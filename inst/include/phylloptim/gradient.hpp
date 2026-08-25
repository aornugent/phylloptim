// -*-c++-*-
#ifndef PHYLLOPTIM_GRADIENT_HPP_
#define PHYLLOPTIM_GRADIENT_HPP_

// Trait gradients, composed here rather than in R (issue #4, PLAN 11d stage 2).
//
// WHAT THIS IS. A transcription of R/gradient.R's `.gradient_ift()` and
// `.gradient_fd()` into C++, plus a loop over observations. It computes exactly
// what `leaf_gradient()` computes and returns the same six things; the whole
// point of it is that a four-parameter gradient crosses the R boundary ONCE
// instead of 112 times.
//
// WHY, measured. After #69 and #70 a four-parameter gradient from R costs
// ~237 us per observation, of which the C++ model work -- two solves -- is
// 6 us, or 1.5%. The other 98.5% is dispatch and the R interpreter.
// `tests/cpp/bench_gradient.cpp` times the same composite here at 1.8 us per
// trait. For a 1,327-observation MCMC that is the difference between hours and
// minutes, and it is why the entry point is vectorised over OBSERVATIONS: one
// crossing per likelihood evaluation, not one per parameter per observation.
//
// ⚠️ WHAT IT DELIBERATELY IS NOT: a likelihood. The likelihood is the caller's
// model -- sigma, robustness, hierarchy -- and baking one in would commit this
// package to it. This returns dY/dtheta for the FOUR model parameters a leaf
// has; the caller's parameterisation Jacobian (leaf-calibration maps 40 fitted
// parameters onto 4 model ones) stays in R, vectorised over observations, where
// it is cheap. That split is the `P_fit > P_model` structure
// `vignette("fitting")` identified as where the exact gradient wins at all.
//
// ⚠️⚠️ THIS IS A SECOND IMPLEMENTATION OF AN ALGORITHM THAT ALREADY EXISTS, and
// that is the risk that dominates everything else here. R/gradient.R stays as
// the reference and the two must agree BIT-FOR-BIT, not closely -- a tolerance
// would let a transcription slip hide inside the solver's ~1e-09 floor. Three
// rules follow, and all three are load-bearing:
//
//   1. R's ARITHMETIC ORDER IS KEPT LITERALLY, including where a division could
//      be folded into a neighbouring one. `-((up - dn) / (2 * h)) / H` is not
//      rewritten as `(dn - up) / (2 * h * H)`.
//   2. NO FUSED MULTIPLY-ADD. See `rounded()` below; `a + b * c` written as one
//      expression compiles to `fmadd` on arm64 and differs from R's two
//      roundings on 28% of random triples, which is a 1e-16 disagreement in a
//      test that asserts equality.
//   3. THE CALL ORDER IS SEQUENCED EXPLICITLY. `f(a) - f(b)` has unspecified
//      operand order in C++ and left-to-right order in R, and these `f`s mutate
//      the leaf. Every difference below names its two halves first.
//
// `tests/testthat/test-gradient-batch.R` is what holds that, over both supply
// paths, both methods and the pinned and shut-down rows.

#include <phylloptim/leaf_model.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/util.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <array>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace phylloptim {
namespace gradient {

#include <phylloptim/inputs.hpp>


// Read straight off the members rather than through `operating_point_values()`,
// which is what R has to use. Bit-identical: that reader copies these same five
// fields into positions 3, 5, 0, 1 and 6 of its twelve, and the three columns it
// computes rather than copies (uptake, lambda, g1_eff) are not among them.
//
// The uptake entries are the leaf's own per-layer consumption at the collar this
// evaluation placed -- `find_psi_stem_from_psi_root` writes it on its way past --
// and not a second computation of it.
inline void outputs(const Leaf& l, OutputValues& y) {
  y[out_assim] = l.assim_colimited_;
  y[out_stom_cond] = l.stom_cond_CO2_;
  y[out_psi_stem] = l.opt_psi_stem_;
  y[out_collar] = l.opt_root_psi_;
  y[out_profit] = l.profit_;
  // A shorter buffer than the request means the drivers' potentials and their
  // depths disagree in length; reading past it would report the layer above.
  if (y.n_uptake() > int(l.soil_consumption_.size())) {
    util::stop("uptake was asked for " + std::to_string(y.n_uptake()) +
               " layers, and this leaf's soil profile has " +
               std::to_string(l.soil_consumption_.size()) + ".");
  }
  for (int i = 0; i < y.n_uptake(); ++i) {
    y[out_uptake_first + i] = l.soil_consumption_[std::size_t(i)];
  }
}

// The outputs with the collar held at `psi` rather than optimised. False where the
// state does not admit that collar -- and REFUSED there rather than evaluated.
//
// ⚠️ ASKING BEFORE EVALUATING IS THE REQUIREMENT, NOT AN ECONOMY. The clamping
// entry point moves the target to the nearest end of the feasible interval and
// answers there, which is silently a one-sided difference over a shorter interval;
// below the wet bound it answers from the model's zero-flux substitution, one dark
// respiration away and of the same type. So a check applied AFTER evaluating has
// to separate two branches by their output alone. `profit_at_fixed_collar`
// establishes feasibility first and does not evaluate when the answer is no.
//
// Where the collar IS admissible this is the same arithmetic at the same collar:
// the clamp is a min/max, so `clamped == psi` is exactly
// `bound_a <= psi <= bound_b`, and both routes then run the same expressions.
// What it removes is the evaluation that used to happen anyway on the way to being
// rejected -- and with it one case the comparison could not see, where a perturbed
// state shuts down at the collar the base point was pinned to and its shut-down
// outputs come back indistinguishable from a held evaluation.
inline bool outputs_at(Leaf& l, double psi, OutputValues& y) {
  if (!l.profit_at_fixed_collar(psi).feasible) {
    return false;
  }
  outputs(l, y);
  return true;
}

// The product, rounded, so that `a + rounded(b * c)` is two IEEE operations and
// never one fused multiply-add.
//
// ⚠️ THIS IS NOT DEFENSIVE, IT IS REQUIRED. `direct + dY_dpsi * dpsi_dtheta` as
// one expression compiles to a single `fmadd` on arm64 -- clang contracts within
// an expression by default, gcc contracts across statements -- and the fused
// result differs from R's two roundings. Measured over 2,000,000 random triples:
// the contracted form disagrees on 565,762 of them, i.e. 28%. Every one of those
// is a last-bit difference in a test whose whole value is that it asserts
// equality. A named intermediate is enough under clang and not under gcc's
// `-ffp-contract=fast`, so the barrier is a volatile store, which the standard
// guarantees rounds. It costs a store and a load per output per parameter,
// against ~2.5 us of solving.
inline double rounded(double x) {
  volatile double v = x;
  return v;
}

// The step: relative to the parameter for values above 1, and plain `step`
// below it. The floor is at 1, which is deliberate rather than an epsilon --
// traits here span `a` = 0.3 to `jmax_25` = 157.44, and a strictly relative step
// would perturb the small ones so little that the difference is dominated by the
// solve's ~1e-09 noise.
//
// ⚠️ That floor is WRONG for a parameter whose natural magnitude is far below 1,
// and `leaf_specific_conductance_max` is: it defaults to 3.14e-05, so flooring
// at 1 would perturb it by 3% and measure a secant across a range over which the
// model is visibly nonlinear. Those two get a plain relative step. For
// `resistance` (~1e3 and up) the two rules coincide; it is listed for the reason
// rather than for the arithmetic.
//
// ⚠️ Root carbon is the third, and there the floor is not merely wrong but
// unsafe: a layer holding less carbon than the floor would be stepped past zero,
// where the architecture model refuses a negative mass and the layer leaves the
// network on one arm only.
inline double step_for(int par, double value, double step, int n_layers) {
  const bool carbon = decode(par, n_layers).block == par_ref::Block::RootCarbon;
  const double floor =
      (par == par_kmax || par == par_resistance || carbon) ? 0.0 : 1.0;
  return std::max(std::abs(value), floor) * step;
}

// --- one observation's drivers ------------------------------------------------
//
// Everything `set_physiology` takes except the two entries a perturbation can
// move, which come out of `theta`. Held C++-side and resolved once, because a
// `RootNetwork` costs 60-100 us to hand across the R boundary -- 60 times a
// trivial `.Call` -- so building one per observation per likelihood evaluation
// would cost more than the gradient it was carrying.
struct Drivers {
  RootNetwork root_network;
  double PPFD = 0.0;
  std::vector<double> psi_soil;
  std::vector<double> soil_depth;
  double atm_vpd = 0.0;
  double ca = 0.0;
  double leaf_temp = 0.0;
  double atm_o2_kpa = 0.0;
  double atm_kpa = 0.0;
};

// How many soil-water rows this observation has. The multi-layer path reads the
// whole profile; the single-potential path reads element 0 and nothing else, so
// it has exactly one soil row however long the caller's vector happens to be.
inline int n_soil_layers(const Drivers& d, bool single) {
  return single ? 1 : static_cast<int>(d.psi_soil.size());
}

// A layer's root carbon, recovered from the network that was built out of it, or
// NA where that layer has none.
//
// The network holds the carbon split three ways rather than the carbon, and the
// vertical third is the half to read it back from -- the same inversion
// `duptake_droot_carbon` performs to answer its own analytic rows. NA rather than
// zero for an unrooted layer, and the distinction is the whole point: the network
// is sized to the deepest rooted layer, so a layer below that one has no carbon
// to move AND no slot to move it in. A zero there would say the outputs are
// insensitive to carbon that could be put there, which is the opposite of true.
inline double root_carbon_of(const Drivers& d, int layer) {
  const std::vector<double>& c = d.root_network.c_r_V;
  return layer >= 0 && layer < int(c.size()) && c[std::size_t(layer)] > 0.0
             ? 3.0 * c[std::size_t(layer)]
             : util::na_value;
}

// The value a parameter currently holds: in `theta` for the sixteen, in the
// observation's drivers for the environment rows. See the environment block above
// for why those two are different places rather than one.
inline double par_value(const double* theta, const Drivers& d, bool single,
                        int par) {
  const par_ref r = decode(par, n_soil_layers(d, single));
  switch (r.block) {
  case par_ref::Block::Parameter:     return theta[r.index];
  case par_ref::Block::Radiation:     return d.PPFD;
  case par_ref::Block::SoilPotential: return d.psi_soil[std::size_t(r.index)];
  case par_ref::Block::RootCarbon:    break;
  }
  return root_carbon_of(d, r.index);
}

// Every requested index names an input this observation has and this package can
// move.
//
// ⚠️ CHECKED PER CALL, because with the environment rows the valid range depends
// on the OBSERVATION rather than on this header: a five-layer row and a one-layer
// row in the same batch do not have the same number of parameters. An index past
// the end used to be impossible (R validates against the fixed sixteen); now it
// is an out-of-bounds read of `psi_soil`, so it is a per-row error.
//
// ⚠️ THE SINGLE-POTENTIAL PATH HAS NO ROOT CARBON, and a silent zero there would
// be the worst available answer. That path builds its own one-element network out
// of `resistance` and never reads the observation's, so a carbon perturbation
// would change nothing the solve consumes and every carbon row would come back
// exactly zero -- indistinguishable from true insensitivity. Refused by name.
// A layer the multi-layer network has no roots in is a per-input refusal instead,
// raised where the rows are assembled.
inline void check_pars(const int* pars, std::size_t npars, int n_layers,
                       bool single, const std::string& caller) {
  for (std::size_t k = 0; k < npars; ++k) {
    if (pars[k] < 0 || pars[k] >= n_pars_total(n_layers)) {
      util::stop(caller + ": parameter index " + std::to_string(pars[k]) +
                 " is out of range; this observation has " +
                 std::to_string(n_layers) + " soil layer(s), so there are " +
                 std::to_string(n_pars_total(n_layers)) + " parameters.");
    }
    if (single && decode(pars[k], n_layers).block == par_ref::Block::RootCarbon) {
      util::stop(caller + ": `" + par_name(pars[k], n_layers) +
                 "` has no row on the single-potential path, which takes one "
                 "series resistance and no root architecture.");
    }
  }
}

enum class Method { Auto, Ift, Fd };

// What the OPERATING POINT is, not whether the call worked -- except for `Error`,
// which is the batch's per-row failure and is the reason this is a status rather
// than an exception (see `batch` below).
enum class Status { Interior, Pinned, NoGradient, Error };

inline std::string status_name(Status s) {
  switch (s) {
    case Status::Interior:   return "interior";
    case Status::Pinned:     return "pinned";
    case Status::NoGradient: return "no-gradient";
    default:                 return "error";
  }
}

struct Settings {
  double step = 1e-6;
  // ⚠️ THE COLLAR'S STEP IS SEPARATE AND MUST STAY SMALL. A pinned point sits one
  // step-in fraction from its bound, so a collar step that crosses a narrow
  // bracket is what detects it; raising this to an input step large enough to
  // move total uptake above the solve's own floor would stop it detecting one.
  double collar = 1e-6;
  double stationarity_tol = 1e-8;
  Method method = Method::Auto;
  bool fast_stem_curve = true;
};

struct Result {
  // The five outputs the gradient is taken at.
  double value[n_outputs];
  // npars * n_outputs, parameter-major: d(output j)/d(pars[k]) at
  // [k * n_outputs + j].
  std::vector<double> grad;
  Status status = Status::Error;
  bool used_ift = false;
  double H = util::na_value;
  double stationarity = util::na_value;
  // dY/dpsi at the base point: the composite's psi-channel, and the third of its
  // three ingredients alongside `H` and the per-parameter M. Reported for the
  // same reason those two are -- it is what the composite stands on -- and NA
  // unless `used_ift`, because the fallback never forms it. `transpose_at` below
  // is its first consumer.
  OutputValues dY_dpsi;
  std::string message;

  void reset(std::size_t npars) {
    for (int j = 0; j < n_outputs; ++j) {
      value[j] = util::na_value;
      dY_dpsi[j] = util::na_value;
    }
    grad.assign(npars * n_outputs, util::na_value);
    status = Status::Error;
    used_ift = false;
    H = util::na_value;
    stationarity = util::na_value;
    message.clear();
  }
};

// --- pushing a parameter vector back onto the leaf ----------------------------
//
// ⚠️ ORDER IS LOAD-BEARING, which is why this is one function rather than three
// lines at each call site. `set_traits()` returns the leaf to its
// just-constructed state, so the drivers have to be re-supplied AFTER it -- and
// `set_physiology()` is what re-derives vcmax_/jmax_/R_d_ behind the temperature
// cache `set_traits()` has just cleared. Setting the traits and then the drivers
// in the other order runs the whole gradient at the first vcmax the object ever
// saw, and reports plausible numbers throughout.
//
// `only` names the single parameter that has moved, or -1 for "all of them".
//
// `PPFD`, `psi_soil` and the network are passed rather than read off `d`, so that
// an environment perturbation is applied by handing over a different environment
// rather than by editing the observation's drivers in place. Every other caller
// takes the two-argument overload below and gets `d`'s own values.
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve, double PPFD,
                  const std::vector<double>& psi_soil,
                  const RootNetwork& root_network) {
  // THE FAST PATH FOR stem_b, which is the whole of PLAN 11f. The stem
  // cumulative-vulnerability integral is homogeneous of degree 1 in stem_b, so
  // the spline for a perturbed stem_b is the existing one with its argument
  // rescaled and the 11.9 us of incomplete gammas a rebuild spends is
  // unnecessary -- 24.5x on that parameter's gradient.
  //
  // ⚠️ Sound only because `only` names a SINGLE parameter, so everything else in
  // `theta` is still what the object was last set to. That is true here and
  // nowhere else, which is why the argument exists rather than the function
  // guessing.
  //
  // `stem_c` is not here because it has no such identity -- it reshapes the curve
  // rather than scaling it, so it rebuilds. The reason once given beside that, that
  // reading the curve from its closed form differentiates a slightly different
  // model and disagrees by 3e-4, is measured and gone: the interpolant carries the
  // closed-form slope now, and the closed-form trait derivative agrees with a
  // difference across the rebuild to 3e-07.
  if (fast_stem_curve && only == par_stem_b) {
    l.perturb_stem_b(theta[par_stem_b]);
    return;
  }
  l.set_traits(theta);
  if (single) {
    // R's `series_resistance()`: a default-constructed network carrying one
    // series resistance in `r_R_V_sum`, which is that field's own meaning with
    // one layer and no vulnerability-weighted term. Built here rather than
    // copied from `d.root_network` so that the two implementations cannot
    // disagree about what the other four fields hold.
    RootNetwork net;
    net.r_R_V_sum.assign(1, theta[par_resistance]);
    l.set_physiology(net, PPFD, psi_soil, d.soil_depth, theta[par_kmax],
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  } else {
    l.set_physiology(root_network, PPFD, psi_soil, d.soil_depth,
                     theta[par_kmax], d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  }
}

// At the observation's own environment, which is what every parameter except the
// environment rows is differentiated in.
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve) {
  apply(l, theta, d, single, only, fast_stem_curve, d.PPFD, d.psi_soil,
        d.root_network);
}

// The network layer `layer`'s carbon would have built, at `value`.
//
// Exact, and it needs neither of the architecture model's two constants nor the
// layer thickness: both resistances the solve reads are proportional to 1/carbon,
// so each one's own constant is recoverable from the network's stored halves and
// re-dividing at the moved carbon reproduces the expression the architecture model
// itself evaluates. That is the same proportionality `duptake_droot_carbon`
// differentiates through, so a differenced carbon row and the analytic bound row
// describe one model rather than two.
//
// The cumulative sum is rebuilt rather than adjusted, in the builder's own
// accumulation order. A zero-carbon layer contributes an exact zero to it either
// way, so a profile with an interior gap is unaffected.
inline void perturb_root_carbon(const RootNetwork& base, int layer, double value,
                                RootNetwork& out) {
  out = base;
  const std::size_t k = std::size_t(layer);
  const double horizontal = base.r_R_H_min[k] * base.c_r_H[k];
  const double vertical = base.r_R_V[k] * base.c_r_V[k];
  out.c_r_V[k] = value / 3.0;
  out.c_r_H[k] = value * 2.0 / 3.0;
  out.r_R_H_min[k] = horizontal / out.c_r_H[k];
  out.r_R_V[k] = vertical / out.c_r_V[k];
  double sum = 0.0;
  for (std::size_t i = 0; i < out.r_R_V.size(); ++i) {
    sum += out.r_R_V[i];
    out.r_R_V_sum[i] = sum;
  }
}

// The buffers a gradient's perturbations reuse, owned by the caller so that a
// copy allocates once per gradient rather than once per perturbation.
struct Scratch {
  std::vector<double> psi_soil;
  RootNetwork root_network;
};

// Put parameter `par` at `value` and everything else at base. One function so
// that both routes below address the sixteen and the environment rows the same
// way, and so that the water and light channels are written where they can be
// read side by side.
//
// ⚠️ THE TWO ENVIRONMENT CHANNELS ARE NOT ALIKE, and the difference is structural
// rather than a matter of which field gets assigned.
//
//   LIGHT is read before the plant's decision. `PPFD` enters assimilation
//   directly through the photosynthesis kernel; at a fixed collar the hydraulics
//   do not move at all, because E and psi_stem are set by the collar and the soil
//   and not by how much carbon the leaf gains. Its row is a plain partial
//   derivative, carrying the intercellular-CO2 root-find's implicit-function
//   term and nothing else.
//
//   WATER is written as a consequence of that decision. `psi_soil` reaches profit
//   only through uptake -- the hydraulic cost depends on stem tension, not on how
//   wet the soil is -- so at a fixed collar potential
//
//       dprofit/dpsi_soil_j = price * dE_up/dpsi_soil_j,
//
//   ONE scalar price times L layer conductances. The whole soil block is rank one
//   across layers, which is a free self-test with real content: over the 198
//   interior points of the golden grid, dprofit/dpsi_soil_j divided by
//   dE_up/dpsi_soil_j varies across layers by at most 2.1e-05. Asserted in
//   `test_leaf.cpp: environment rows: the water channel is rank one across
//   layers`.
//
//   ⚠️ That price is NOT `marginal_cost_water_multilayer()`, though that is the
//   lambda `find_root_collar_psi` equalises and so the obvious thing to reach
//   for. The multilayer lambda is dA/dE along the COLLAR direction, where the
//   stem's own cost is already netted against dA/dE by the first-order condition.
//   Along the soil direction the collar is held FIXED, so water arriving from a
//   wetter layer still has to leave through the stem, and the leaf pays
//   `marginal_cost_water()` for that:
//
//       price = marginal_cost_water_multilayer() - marginal_cost_water().
//
//   Measured against a differenced solve over the same 198 points: the difference
//   agrees to 6.3e-05 relative, while the bare multilayer lambda overstates the
//   price by a factor of 1.20 to 2.37.
//
//   ⚠️ LAYERS BELOW THE DEEPEST ROOTED ONE GET EXACTLY ZERO UPTAKE, so their
//   rows are exactly 0.0 on both routes. That zero is CORRECT and is asserted
//   (`environment rows: an unrooted layer's rows are exactly zero`). It is called
//   out because in this codebase an exact zero is otherwise the signature of a
//   missing accumulator -- and equally so that a genuinely missing row is not
//   read as this.
//
//   ROOT CARBON is the third channel and it reaches the solve the way the soil
//   does -- only through uptake -- but it moves the conductance rather than the
//   driving potential, so it is applied by handing over a network built at the
//   moved carbon. A layer with none has no row: see `root_carbon_of`.
//
// The scratch buffers are the caller's rather than locals so that the copies
// allocate once per gradient instead of once per perturbation.
inline void set_one(Leaf& l, double* th, const double* theta, const Drivers& d,
                    bool single, int par, double value, bool fast_stem_curve,
                    Scratch& scratch) {
  std::copy(theta, theta + n_pars, th);
  const par_ref r = decode(par, n_soil_layers(d, single));
  switch (r.block) {
  case par_ref::Block::Parameter:
    th[par] = value;
    apply(l, th, d, single, par, fast_stem_curve);
    return;
  case par_ref::Block::Radiation:
    apply(l, th, d, single, par, fast_stem_curve, value, d.psi_soil,
          d.root_network);
    return;
  case par_ref::Block::SoilPotential:
    scratch.psi_soil = d.psi_soil;
    scratch.psi_soil[std::size_t(r.index)] = value;
    apply(l, th, d, single, par, fast_stem_curve, d.PPFD, scratch.psi_soil,
          d.root_network);
    return;
  case par_ref::Block::RootCarbon:
    break;
  }
  perturb_root_carbon(d.root_network, r.index, value, scratch.root_network);
  apply(l, th, d, single, par, fast_stem_curve, d.PPFD, d.psi_soil,
        scratch.root_network);
}

// --- the two routes ----------------------------------------------------------

// ⚠️ THE INVARIANT BOTH LOOPS BELOW MAINTAIN, AND IT IS THE FIX FOR #72: every
// parameter's gradient is taken from the BASE point.
//
// The loops restore base once at the END, not between parameters, because every
// parameter's setter normally goes through the full `set_traits()` +
// `set_physiology()` path and restores everything on the way. The `stem_b`
// shortcut is the one that does not: `perturb_stem_b()` rescales the stem spline
// and touches nothing else, which is sound only if the rest of the object is
// already at base. So a `stem_b` that is not the FIRST entry of `pars` was
// differentiated at a point displaced by one step in whichever parameter
// preceded it -- up to 3.4e-5 relative, four orders above the ~1e-9 this is
// supposed to deliver, on both routes.
//
// Cheap in the case that matters: `set_traits()` decides the two spline rebuilds
// by comparing the pairs it is given, so after a parameter that owns no
// vulnerability curve nothing is rebuilt and this costs one trait write plus one
// driver write. After `stem_c`/`root_b`/`root_c` it does rebuild, and those are
// the parameters whose own gradients cost a rebuild per side anyway.
//
// ⚠️ Written as "this parameter takes a shortcut", not as "this parameter is
// stem_b". `root_b` obeys the same homogeneity identity and would get the same
// treatment, at which point a name-based test would silently stop covering it.
inline bool takes_shortcut(int par, const Settings& s) {
  return s.fast_stem_curve && par == par_stem_b;
}

// --- what every route shares --------------------------------------------------
//
// The base point, the collar channel, and one input's two perturbed evaluations.
// Both entry points below and `rows_at` stand on these, so a change to the
// algebra reaches all of them or none.

// The branch a solve took: the kind of operating point, and which of the two
// limits won the dry bound. The two travel together because a difference taken
// across either is a difference of two functions rather than of one.
struct Branch {
  Leaf::OperatingPointKind kind = Leaf::OperatingPointKind::Unsolved;
  Leaf::DryBoundArm arm = Leaf::DryBoundArm::None;
  bool operator==(const Branch& o) const {
    return kind == o.kind && arm == o.arm;
  }
};

inline Branch branch_here(const Leaf& l) {
  return {l.operating_point_kind(), l.dry_bound_arm()};
}

// The solve, the five outputs, and the marginal profit at the collar the solve
// returned.
//
// ⚠️ THE CURVATURE IS NOT HERE, and it used to be. Taking it differenced the
// marginal profit on both sides of p*, which left the collar one step below the
// point and cost every later read a re-placing; `differenced_curvature` below is
// that difference, and only the route that has no closed form for it pays.
struct BasePoint {
  double psi_star = util::na_value;
  double resid = util::na_value;
  // Whether the evaluation at p* admits a derivative at all, from `dprofit`'s own
  // out-parameter. It is the same evaluation, so a later reader needs no second
  // one to find out.
  bool placed = false;
  // ⚠️ Read the moment the solve ends, because the first evaluation at a held
  // collar overwrites it: `evaluate_root_collar_psi` tags the point Prescribed.
  Branch branch;
};

// How far a row that may recover shrinks its step before giving up. Two decades
// and no more: the largest step that is still local is the best-conditioned one,
// and below that floor a step stops moving the solve above its own noise, so the
// difference measures rounding rather than the model. Both routes that shrink read
// it here, so neither can drift to a different floor.
inline constexpr int shrink_decades = 2;

// The step in the collar potential, floored at 1 MPa for `step_for`'s reason.
inline double collar_step(double psi_star, const Settings& s) {
  return std::max(std::abs(psi_star), 1.0) * s.collar;
}

// ⚠️ `n_uptake` sizes the value buffer, and the outputs must be read HERE. The
// uptake block is variable-length, so a default-sized buffer silently records
// none of it, and every read past this point describes whatever the last
// evaluation left.
//
// The solve closes on the condition at the collar it returns, so the coefficients
// this reads describe the point and nothing here evaluates anything.
inline BasePoint base_point(const Leaf& l) {
  BasePoint b;
  b.branch = branch_here(l);
  b.psi_star = l.opt_root_psi_;
  b.resid = l.collar_resid_;
  b.placed = l.collar_resid_placed_;
  return b;
}

// The curvature of profit at the point, by differencing the marginal profit over
// the collar's own step. `Leaf::condition_collar_slope` is the same quantity in
// closed form, and this is what answers where that refuses.
//
// ⚠️ IT LEAVES THE COLLAR ONE STEP BELOW p*, so a caller that reads the point's
// own coefficients afterwards is reading a neighbouring state -- 4.5e-06 relative,
// which looks like nothing. `base_point` closed with this until the closed form
// existed, and every reader after it had to place it again.
inline double differenced_curvature(Leaf& l, double psi_star, double resid,
                                    const Settings& s) {
  const double h_psi = collar_step(psi_star, s);
  // Named halves: `f(a) - f(b)` has unspecified operand order in C++ and
  // left-to-right order in R, and `dprofit_droot_collar_psi` mutates the leaf.
  const double d_hi = l.dprofit_droot_collar_psi(psi_star + h_psi);
  const double d_lo = l.dprofit_droot_collar_psi(psi_star - h_psi);
  // dprofit returns a bare, exact 0.0 SENTINEL rather than a derivative where
  // the collar is shut down or the ci solve is infeasible, and the `feasible`
  // out-parameter that would distinguish it is not carried through the R
  // binding. So test for it EXACTLY, the same argument `outputs_at` makes with
  // `util::identical` for the clamp: an unclamped evaluation reaches 0.0 only at
  // a genuine stationary point, so an exact zero on one arm of the difference is
  // the sentinel. Left in, it does not make H small, it makes H wrong -- |H| out
  // by a median factor of 8.4e04 over the 288-point grid.
  //
  // One bad arm: difference the good arm against `resid` at psi* instead. That
  // is a real one-sided second derivative of profit, it is continuous across the
  // feasibility boundary where the centred one jumps, and it leaves the reported
  // classification unchanged everywhere. Both arms bad: there is no curvature to
  // report and the 0.0 below is the shut-down signature.
  const bool hi_sentinel = util::identical(d_hi, 0.0);
  const bool lo_sentinel = util::identical(d_lo, 0.0);
  return (hi_sentinel && lo_sentinel) ? 0.0
         : hi_sentinel               ? (resid - d_lo) / h_psi
         : lo_sentinel               ? (d_hi - resid) / h_psi
                                     : (d_hi - d_lo) / (2.0 * h_psi);
}

// dY/dpsi at fixed traits, and a SECOND, INDEPENDENT detector of a pinned
// optimum. At a pinned point psi* sits one step-in fraction (1e-06 of the
// bracket width) from its bound, so a step of `step * psi` crosses it whenever
// the bracket is narrower than psi.
//
// ⚠️ IT DETECTS 18 OF THE GRID'S 42 PINNED POINTS, not all of them, and the
// difference is the shrink below. This comment claimed all 42 and that was
// measured before the shrink existed: two decades of it recover 24 of the 42, so
// the number of points where this refuses -- which is what `at` reads as a pin --
// fell by more than half when the recovery landed. It was never a substitute for
// the stationarity test and it is less of one now.
//
// False where the difference cannot be centred on psi*. `Leaf::collar_rows` is
// the same channel read rather than differenced and it answers at all 42; this
// route is kept because `at` is refereed bit for bit against a captured
// reference, and the two must run the same arithmetic there.
inline bool collar_response(Leaf& l, double psi_star, const Settings& s,
                           OutputValues& dY_dpsi) {
  OutputValues hi(dY_dpsi.n_uptake());
  OutputValues lo(dY_dpsi.n_uptake());
  // ⚠️ IT SHRINKS, for the reason `held_row` does, and it did not until a stand
  // showed why. The operating point can sit closer to a bound than one fixed step:
  // below the wet bound uptake is negative and the state has no stem potential, so
  // an arm there does not merely lose accuracy, it does not exist. Over the golden
  // grid no arm ever landed outside one, which is why a fixed step looked settled;
  // on a drought stand an interior node landed within it and every water row at
  // that node was lost -- reported as an interior point refusing, which is a
  // branch that answers everywhere else.
  //
  // Shrinking rather than going one-sided keeps the difference centred on the
  // point, which is what makes the objective's channel exactly zero and the
  // point's exactly one. Same floor as the other two routes that shrink, for the
  // same reason: two decades below this and a difference measures rounding.
  for (int decade = 0; decade <= shrink_decades; ++decade) {
    const double h_psi = collar_step(psi_star, s) * std::pow(0.1, decade);
    // Both, unconditionally, before the test -- R computes `hi` and `lo` on
    // consecutive lines and only then checks either, and each call moves the leaf.
    const bool hi_ok = outputs_at(l, psi_star + h_psi, hi);
    const bool lo_ok = outputs_at(l, psi_star - h_psi, lo);
    if (!hi_ok || !lo_ok) {
      continue;
    }
    for (int j = 0; j < dY_dpsi.size(); ++j) {
      dY_dpsi[j] = (hi[j] - lo[j]) / (2.0 * h_psi);
    }
    return true;
  }
  return false;
}

// One input's two perturbed evaluations at a HELD collar: the requested outputs'
// direct rows, and dR/du -- the collar derivative of marginal profit. Neither
// evaluation re-solves the model. How many outputs is `direct`'s own length.
//
// `at_base` carries the invariant above: it comes in true only where the leaf is
// at base parameters, and goes out false.
// False where no step keeps both arms inside the perturbed feasible interval;
// `direct` and `dresidual` are not written then.
//
// ⚠️ `may_recover` is ONE fact with two consequences, which is why it is a
// permission and not a count. The reference implementation this file is checked
// against stops at the first infeasible arm, so a route that recovered from one
// would answer states the reference does not and the two would stop being one
// function checked against another. A route answering for itself takes both
// recoveries -- the one-sided arm first, then the shrink.
inline bool held_row(Leaf& l, const double* theta, const Drivers& d,
                       bool single, int par, double psi_star, const Settings& s,
                       bool may_recover, bool& at_base,
                       Scratch& scratch, OutputValues& direct,
                       double& dresidual) {
  if (takes_shortcut(par, s) && !at_base) {
    apply(l, theta, d, single, -1, s.fast_stem_curve);
  }
  at_base = false;
  double th[n_pars];
  OutputValues up(direct.n_uptake());
  OutputValues dn(direct.n_uptake());
  double up_resid = 0.0;
  double dn_resid = 0.0;
  bool up_real = false;
  bool dn_real = false;
  const double base = par_value(theta, d, single, par);
  double h = step_for(par, base, s.step, n_soil_layers(d, single));
  // One evaluation of each output at the held collar with nothing moved, for the
  // one-sided arm below. Taken lazily, because the centred arm is what almost
  // every point takes and it does not need this.
  OutputValues here_out(direct.n_uptake());
  double here_resid = 0.0;
  bool here_resid_real = false;
  bool have_here = false;
  auto at_base_here = [&]() -> bool {
    if (have_here) {
      return true;
    }
    apply(l, theta, d, single, -1, s.fast_stem_curve);
    if (!outputs_at(l, psi_star, here_out)) {
      return false;
    }
    here_resid = l.dprofit_droot_collar_psi(psi_star, &here_resid_real);
    have_here = true;
    return true;
  };
  for (int decade = 0;; ++decade) {
    bool up_ok = false;
    bool dn_ok = false;
    for (int side = 0; side < 2; ++side) {
      // Up first, then down: R evaluates `up <- side(1)` before `dn <- side(-1)`
      // and both mutate the leaf.
      set_one(l, th, theta, d, single, par, side == 0 ? base + h : base - h,
              s.fast_stem_curve, scratch);
      OutputValues& dst = side == 0 ? up : dn;
      // Evaluate first, then read dprofit at the same fixed collar -- R's order,
      // and `evaluate_root_collar_psi` is what places the state `dprofit` reads.
      // Both arms are taken before either is tested, so which one moved the
      // interval does not change the order the leaf is moved in.
      const bool ok = outputs_at(l, psi_star, dst);
      (side == 0 ? up_ok : dn_ok) = ok;
      double& resid = side == 0 ? up_resid : dn_resid;
      // ⚠️ THE `feasible` OUT-PARAMETER IS WHAT SEPARATES A DERIVATIVE FROM A
      // SENTINEL. `dprofit` returns a bare, exact 0.0 in the no-flow or
      // infeasible state, and differencing it produces a plausible finite number
      // that is not a mixed partial -- the same mechanism `base_point` guards the
      // curvature against, where leaving it in put |H| out by a median factor of
      // 8.4e04. It reaches here and not the base point because a pinned collar
      // sits one step-in from where uptake is exactly zero.
      //
      // The OUTPUTS' rows survive it: they come from the evaluation, not from
      // dprofit. So a sentinel costs `dresidual` alone, and the consumer that
      // reads a non-finite one already knows to keep the objective's row and drop
      // the rest -- which is the envelope theorem, not a degradation.
      resid = l.dprofit_droot_collar_psi(psi_star, side == 0 ? &up_real
                                                            : &dn_real);
    }
    if (up_ok && dn_ok) {
      // M = d2profit/dpsi dtheta, with psi held FIXED at psi*.
      dresidual = (up_real && dn_real) ? (up_resid - dn_resid) / (2.0 * h)
                                       : util::na_value;
      for (int j = 0; j < direct.size(); ++j) {
        direct[j] = (up[j] - dn[j]) / (2.0 * h);
      }
      return true;
    }
    // ONE ARM INSIDE: take it one-sided rather than shrinking, and take it to
    // SECOND order. An interior point can sit within a step of a bound, and there
    // the held collar is feasible on one side only -- so this is the arm that
    // answers a drying stand, not an edge case. First order would be 150x worse
    // and would carry a sign that follows which bound the point drifted toward
    // rather than the input, so the same row would differ by side.
    //
    // Shrinking is LAST because the reads carry the solve's own floor: dividing
    // them by a step two decades smaller costs more than this truncation.
    if ((up_ok != dn_ok) && may_recover && at_base_here()) {
      const int sign = up_ok ? 1 : -1;
      OutputValues& near = up_ok ? up : dn;
      const double near_resid = up_ok ? up_resid : dn_resid;
      OutputValues far(direct.n_uptake());
      set_one(l, th, theta, d, single, par, base + sign * 2.0 * h,
              s.fast_stem_curve, scratch);
      if (outputs_at(l, psi_star, far)) {
        bool far_real = false;
        const double far_resid =
            l.dprofit_droot_collar_psi(psi_star, &far_real);
        const bool near_real = up_ok ? up_real : dn_real;
        const double scale = sign / (2.0 * h);
        dresidual = (here_resid_real && near_real && far_real)
                        ? scale *
                              (-3.0 * here_resid + 4.0 * near_resid - far_resid)
                        : util::na_value;
        for (int j = 0; j < direct.size(); ++j) {
          direct[j] = scale * (-3.0 * here_out[j] + 4.0 * near[j] - far[j]);
        }
        return true;
      }
    }
    // A step that carries the held collar out of the perturbed feasible interval
    // is the state a pin sits in, where the interval closes around the point. The
    // shrink is `solved_row`'s and stops at the same floor for the same reason.
    if (!may_recover || decade == shrink_decades) {
      return false;
    }
    h *= 0.1;
  }
}

// For the two entry points that answer a whole gradient rather than one column
// of one: a row they cannot take is the call's answer and not an input's. No
// shrinking, and the message is the R implementation's word for word.
inline void held_row_or_stop(Leaf& l, const double* theta, const Drivers& d,
                             bool single, int par, double psi_star,
                             const Settings& s, const std::string& caller,
                             bool& at_base, Scratch& scratch,
                             OutputValues& direct, double& dresidual) {
  if (!held_row(l, theta, d, single, par, psi_star, s, false, at_base, scratch,
                direct, dresidual)) {
    util::stop(caller + ": perturbing `" +
               par_name(par, n_soil_layers(d, single)) +
               "` moved the feasible collar interval past psi*, so the "
               "operating point could not be evaluated there. This point is "
               "on an active-set boundary; lower `stationarity_tol` or "
               "difference the solve directly.");
  }
}

// Where a pinned point's arms are placed: the bound that defines it, and the
// step-in the solve left between that bound and the collar.
struct FollowBound {
  Leaf::WhichBound which;
  double step_in;
};

// One input's central difference of the WHOLE solve. Correct at a pinned optimum
// because it differences the CONSTRAINED answer, which is exactly what the
// composite cannot do.
//
// `stay` is the branch both arms have to land on, or null for a caller that
// differences whatever they land on. Where it is given, the step shrinks by a
// decade at a time looking for two arms that agree with it, and the row is refused
// rather than taken across a change of branch, down to `shrink_decades`.
//
// False only where `stay` was given and no step reached it; `row` is not written
// then.
inline bool solved_row(Leaf& l, const double* theta, const Drivers& d,
                       bool single, int par, const Settings& s,
                       const Branch* stay, const FollowBound* follow,
                       bool& at_base, Scratch& scratch, OutputValues& row) {
  if (takes_shortcut(par, s) && !at_base) {
    apply(l, theta, d, single, -1, s.fast_stem_curve);
  }
  at_base = false;
  double th[n_pars];
  OutputValues up(row.n_uptake());
  OutputValues dn(row.n_uptake());
  const double base = par_value(theta, d, single, par);
  double h = step_for(par, base, s.step, n_soil_layers(d, single));
  auto difference = [&]() -> void {
    for (int j = 0; j < row.size(); ++j) {
      row[j] = (up[j] - dn[j]) / (2.0 * h);
    }
  };
  for (int decade = 0;; ++decade) {
    bool on_branch = true;
    for (int side = 0; side < 2; ++side) {
      set_one(l, th, theta, d, single, par, side == 0 ? base + h : base - h,
              s.fast_stem_curve, scratch);
      l.find_root_collar_psi();
      // Both arms are taken before either is tested, so which one crossed does
      // not change the order the leaf is moved in.
      on_branch = on_branch && (stay == nullptr || branch_here(l) == *stay);
      outputs(l, side == 0 ? up : dn);
    }
    if (on_branch) {
      difference();
      return true;
    }
    // THE ARMS LEFT THE BRANCH. Where a bound is what defines the point, place
    // them on it rather than letting each solve choose: a pinned point IS its
    // bound, one step-in from it, so an arm at the PERTURBED bound plus that
    // same step-in is on the branch by construction. A drying stand refuses
    // `stem_c`, `psi_crit` and its soil rows without this.
    //
    // Tried only after the re-solve, and that ordering is not cosmetic: where
    // both arms stay on the branch the two placements are not the same row, and
    // the re-solved one is what every other entry point reports.
    //
    // ⚠️ The step-in is load-bearing. At the wet bound total uptake is exactly
    // zero, so an arm at the bare bound is degenerate: one measured 268.561 for
    // a derivative of -0.0814.
    if (follow != nullptr) {
      bool placed = true;
      for (int side = 0; side < 2; ++side) {
        set_one(l, th, theta, d, single, par, side == 0 ? base + h : base - h,
                s.fast_stem_curve, scratch);
        // ⚠️ The solve runs for its SUPPLY STATE and not for its answer: setting
        // traits leaves that stale, and the bound then reads as the wettest
        // layer. The collar it lands on is what is discarded here.
        l.find_root_collar_psi();
        const Leaf::BoundRow moved = l.bound_row(follow->which);
        const bool here =
            moved.finite &&
            outputs_at(l, moved.bound + follow->step_in, side == 0 ? up : dn);
        placed = placed && here;
      }
      if (placed) {
        difference();
        return true;
      }
    }
    if (decade == shrink_decades) {
      return false;
    }
    h *= 0.1;
  }
}

// The implicit-function composite. Two perturbed evaluations per parameter,
// neither of which re-solves the model: `dprofit` at the UNPERTURBED psi* gives
// the mixed partial, and the outputs at that same psi* give the direct term.
inline void gradient_ift(Leaf& l, const double* theta, const Drivers& d,
                         bool single, const int* pars, std::size_t npars,
                         double psi_star, double H,
                         const OutputValues& dY_dpsi, const Settings& s,
                         double* out) {
  Scratch scratch;
  OutputValues direct;
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    double M = 0.0;
    held_row_or_stop(l, theta, d, single, pars[k], psi_star, s,
                     "leaf_gradient()", at_base, scratch, direct, M);
    const double dpsi_dtheta = -(M / H);
    // Each output once, by what it is to the point. `output_role` is where those two
    // facts live; what belongs here is why they are BRANCHES rather than a
    // channel multiplied in.
    //
    // ⚠️ THE OBJECTIVE'S ZERO MUST NOT BE MULTIPLIED, and this is a second reason
    // beyond the envelope theorem itself. `dY_dpsi[out_profit]` is a central
    // difference of a FLAT maximum: ~1e-09 of cancellation in profit over a
    // ~1e-06 step in psi. Measured at psi_soil = 2.0 the exact dprofit at psi* is
    // 8.7e-11, while that difference reads -1.8e-04 -- so composing it in would
    // move this column by up to 6.6e-05 relative, five orders above the solve's
    // ~1e-09 floor and entirely plausible-looking.
    //
    // ⚠️ THIS ROUTE RUNS AT AN INTERIOR OPTIMUM ONLY, which is what makes the
    // objective's zero true here. At a pin the same output's channel is the
    // constraint's shadow price, and the row comes from differencing the solve --
    // which `gradient_fd` does for every output with no special case at all.
    // `at()`'s `status` is the one place that choice is made.
    for (int j = 0; j < n_outputs; ++j) {
      double& into = out[k * n_outputs + j];
      switch (output_role(j)) {
      case OutputRole::Point:
        // It IS psi*, held fixed, so its direct term is zero by construction and
        // the composite reduces to dpsi*/dtheta.
        into = dpsi_dtheta;
        break;
      case OutputRole::Objective:
        into = direct[j];
        break;
      default:
        into = direct[j] + rounded(dY_dpsi[j] * dpsi_dtheta);
      }
    }
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);
}

// The fallback: a central difference of the whole solve. Correct at a pinned
// optimum because it differences the CONSTRAINED answer, which is exactly what
// the composite cannot do.
//
// ⚠️ ALL FIVE OUTPUTS COME OUT OF THE ONE DIFFERENCE, `profit` INCLUDED, and the
// absence of a special case here is deliberate. `gradient_ift`'s envelope
// shortcut assumes dprofit/dpsi = 0, which is false at the pinned optimum this
// function exists for; differencing the solve needs no such premise, because the
// solve it differences is the constrained maximum.
inline void gradient_fd(Leaf& l, const double* theta, const Drivers& d,
                        bool single, const int* pars, std::size_t npars,
                        const Settings& s, double* out) {
  Scratch scratch;
  OutputValues row;
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    // No branch to hold: this route differences whatever the two arms land on and
    // reports the base point's own kind beside the answer, so it cannot refuse.
    // The copy is what leaves a refused row NA rather than partly written.
    if (solved_row(l, theta, d, single, pars[k], s, nullptr, nullptr, at_base,
                   scratch, row)) {
      for (int j = 0; j < n_outputs; ++j) {
        out[k * n_outputs + j] = row[j];
      }
    }
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);
}

// --- one observation ---------------------------------------------------------
//
// Fills `out` progressively and THROWS on the three conditions `leaf_gradient()`
// stops on, so that a partially-determined row still carries its diagnostics
// when `batch` catches. See `batch` for why the batch does not propagate.
inline void at(Leaf& l, const double* theta, const Drivers& d, bool single,
               const int* pars, std::size_t npars, const Settings& s,
               Result& out) {
  out.reset(npars);
  check_pars(pars, npars, n_soil_layers(d, single), single, "leaf_gradient()");

  // This route owns its solve, because it is an entry point rather than a read of
  // a leaf a consumer already solved.
  apply(l, theta, d, single, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  const BasePoint b = base_point(l);
  OutputValues base_values(0);
  outputs(l, base_values);
  const double psi_star = b.psi_star;
  const double resid = b.resid;
  // ⚠️ TAKEN HERE AND DIFFERENCED, not read from the leaf's closed form. This
  // route is refereed bit for bit against a captured reference, so the two must
  // run the same arithmetic; and it is taken immediately, because `base_point`
  // leaves the leaf at p* and this is the same three reads in the same order it
  // used to make itself.
  const double H = differenced_curvature(l, psi_star, resid, s);
  for (int j = 0; j < n_outputs; ++j) {
    out.value[j] = base_values[j];
  }

  // Is the composite's premise true HERE? Stationarity is what the whole
  // derivation rests on and it fails at a pinned optimum, where psi* is a bound,
  // dprofit is not zero at the answer, and -M/H is not the bound's derivative.
  // The formula does not fail loudly: at a wet-pinned point the true gradient is
  // ~1e-08 and the bare composite returns O(1).
  //
  // So the premise is TESTED. The test is the implied Newton step
  // |dprofit(psi*) / H|, a distance in MPa that needs no scale of its own: over
  // this package's 288-point grid the worst interior point is 4.8e-11 and the
  // mildest pinned one 1.2e-02, so the 1e-08 default sits in an empty band six
  // orders wide above it. That lower edge used to be recorded here as 6.3e-06,
  // which was the unguarded central difference below inflating |H| at pinned
  // points and so shrinking |resid / H|; the tolerance was right, the reason
  // given for it was not.
  //
  // H == 0 with resid == 0 is the shut-down signature: dprofit returns a
  // sentinel zero there rather than a derivative, so the ratio would be 0/0.
  // H > 0 would not be a maximum. Both mean the composite has nothing to stand
  // on.
  const bool usable = std::isfinite(H) && H < 0.0 && std::isfinite(resid);
  out.H = H;
  out.stationarity = usable ? std::abs(resid / H)
                            : std::numeric_limits<double>::infinity();
  out.status = !usable ? Status::NoGradient
               : (out.stationarity > s.stationarity_tol ? Status::Pinned
                                                        : Status::Interior);

  // `status` describes the POINT and is reported whichever route runs;
  // `use_ift` is the route. They differ only when the caller has forced one.
  bool use_ift = s.method == Method::Auto ? out.status == Status::Interior
                                          : s.method == Method::Ift;
  if (use_ift && !usable) {
    util::stop("leaf_gradient(): method = \"ift\" was asked for at a point with "
               "no usable curvature (H = " + util::to_string(H) + "), so -M/H "
               "has nothing to stand on. This is a shut-down or otherwise "
               "determined operating point; use method = \"auto\".");
  }

  // Stationarity is the composite's premise and a pinned optimum has none: psi*
  // is a bound, dprofit is not zero at the answer, and -M/H is not the bound's
  // derivative. Refused rather than returned, because the wrong answer here is
  // O(1) against a truth of ~1e-08 and so reads as a gradient. The narrow-bracket
  // test below catches most of these as a side effect of the step not centring;
  // this is the predicate itself, so a pinned optimum on a wide bracket is caught
  // too.
  if (use_ift && out.status != Status::Interior) {
    util::stop("leaf_gradient(): method = \"ift\" was asked for at a pinned "
               "operating point (stationarity = " +
               util::to_string(out.stationarity) + " against a tolerance of " +
               util::to_string(s.stationarity_tol) + "), where psi* is a bound "
               "and -M/H is not its derivative. Use method = \"auto\".");
  }

  // Written straight into the result rather than into a local, so that the
  // transpose can read the same numbers instead of measuring them again.
  OutputValues& dY_dpsi = out.dY_dpsi;
  if (use_ift && !collar_response(l, psi_star, s, dY_dpsi)) {
    if (s.method == Method::Ift) {
      util::stop("leaf_gradient(): method = \"ift\" was asked for at a point "
                 "whose feasible collar interval is narrower than one step, "
                 "so dY/dpsi cannot be centred on psi*. Use "
                 "method = \"auto\".");
    }
    use_ift = false;
    out.status = Status::Pinned;
  }

  out.used_ift = use_ift;
  if (use_ift) {
    gradient_ift(l, theta, d, single, pars, npars, psi_star, H, dY_dpsi, s,
                 out.grad.data());
  } else {
    gradient_fd(l, theta, d, single, pars, npars, s, out.grad.data());
  }
}

// --- the batch ---------------------------------------------------------------
//
// `theta` is a COLUMN-MAJOR matrix of `theta_nrow` x n_pars, as R hands one
// over: either one row per observation, or exactly one row shared by all of
// them.
//
// ⚠️ PER-ROW STATUS, NOT AN EXCEPTION, and this is the design decision the batch
// exists to make. A proposal during a fit WILL reach operating points the solve
// cannot handle -- that is what a proposal distribution does -- and that has to
// cost those rows rather than the whole dataset. Throwing would take out a
// likelihood evaluation, and with it the draw, for one observation the sampler
// was entitled to reject on its own. `leaf_predict()` isolates per row for the
// same reason.
//
// A failed row's gradient is ALL NA rather than partially filled. A parameter
// loop that threw halfway has valid entries before the throw, and returning them
// alongside `status == "error"` would be an invitation to use them: a caller who
// checks the status per row and not per cell would be reading a gradient with a
// hole in it. Its `value`, `H` and `stationarity` are kept where they were
// determined, because those describe the point rather than the derivative.
inline std::vector<Result> batch(Leaf& l, const double* theta,
                                std::size_t theta_nrow,
                                const std::vector<Drivers>& drivers,
                                bool single, const int* pars,
                                std::size_t npars, const Settings& s) {
  const std::size_t n = drivers.size();
  std::vector<Result> out(n);
  double th[n_pars];
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t row = theta_nrow == 1 ? 0 : i;
    for (int j = 0; j < n_pars; ++j) {
      th[j] = theta[row + std::size_t(j) * theta_nrow];
    }
    try {
      at(l, th, drivers[i], single, pars, npars, s, out[i]);
    } catch (const std::exception& e) {
      out[i].status = Status::Error;
      out[i].used_ift = false;
      out[i].message = e.what();
      out[i].grad.assign(npars * n_outputs, util::na_value);
      // Put the leaf back at this row's base parameters before the next one.
      // The next row's own `apply(only = -1)` would do it -- `set_traits` forces
      // the vulnerability rebuild that `perturb_stem_b` displaced -- but the
      // LAST row has no next one, and a batch that ended on the fast path would
      // hand back a leaf quietly running on a rescaled stem curve (hazard 8).
      try {
        apply(l, th, drivers[i], single, -1, s.fast_stem_curve);
      } catch (const std::exception&) {  // NOLINT: nothing better to do here
      }
    }
  }
  return out;
}

// --- one observation, in parts -------------------------------------------------
//
// `at` above returns TOTAL rows: it forms the quotient and the composite itself,
// so a consumer recording them has `n_output * n_input` terms to tape. This
// returns the same algebra unassembled -- the held partials, whichever
// condition defines the operating point together with that condition's own
// slope, and each output's sensitivity to the point -- so the point is ONE node
// and the tape holds `n_output + n_input`. At 5 outputs and 26 inputs that is 32
// terms against 130.
//
// ⚠️ THE QUOTIENT IS NOT TAKEN HERE, and neither is the refusal that goes with
// it. `dresidual` and `residual_slope` come back separately because dividing
// them is a property of the implicit function theorem rather than of leaves, and
// because at a pin it is the same division on the bound's condition.
//
// ⚠️ A CONSTRAINED POINT RETURNS PARTS TOO, and it used to return totals. p*
// sits one step-in fraction from the bound, so a step that moves the bound
// carries p* out of the perturbed feasible interval and there is no held
// evaluation to DIFFERENCE -- which is an argument about steps, and a closed form
// takes none. So the same held partial that answers at an interior point answers
// here, the bound's own row is what `dresidual` carries, and the point's movement
// is reported once rather than folded into every output's row.
//
// What still returns a total is what still has no closed form: a shut collar,
// where no condition defines the point at all, and the single path's series
// resistance. `dresidual` is zero for those, so the consumer's one assembly
// returns the total unchanged and nothing subtracts.
//
// THIS IS THE ONE ROUTE THAT CAN BE ASKED FOR PER-LAYER UPTAKE. `output` indexes
// the enumeration above, so an index at or past `out_uptake_first` names a layer's
// consumption, and an index past the last layer this observation has is refused by
// name. The five stay reportable, and the consumer this exists for asks for
// profit and the uptake block.

using OperatingPointKind = Leaf::OperatingPointKind;

// Which bound a pinned point is sitting on. The dry end is a min of two limits
// that are DIFFERENT FUNCTIONS of the inputs, so this is a three-way question
// rather than "pinned or not".
inline bool pinned_bound(OperatingPointKind kind, Leaf::WhichBound& bound) {
  switch (kind) {
  case OperatingPointKind::PinnedWet:
    bound = Leaf::WhichBound::Wet;
    return true;
  case OperatingPointKind::PinnedDryRootCrit:
    bound = Leaf::WhichBound::DryRootCrit;
    return true;
  case OperatingPointKind::PinnedDryRootPsiCrit:
    bound = Leaf::WhichBound::DryRootPsiCrit;
    return true;
  default:
    return false;
  }
}


// --- a route that does not work, recorded so it is not tried again -----------
//
// A stand adjoint needs per-layer uptake rows, and those are TOTAL derivatives:
// uptake is set as a side effect at the operating point, so it consumes the
// argmax rather than being it and the operating point's own movement is part of
// the answer,
//
//     dE_i/du  =  dE_i/du|_p  +  (dE_i/dp) * (dp*/du),
//     dp*/du   = -(d2profit/dp du) / (d2profit/dp2).
//
// The cheap-looking route is to notice that d2profit/dp du is the collar
// derivative of dprofit/du, which `profit_env_derivatives` already computes
// analytically, and so to difference THAT in the collar: two placements, no
// re-driving, every column at once.
//
// ⚠️ IT IS WRONG, AND MEASURABLY SO -- the collar row comes back at -1.40x the
// answer a re-solve gives. The reason is worth keeping. `marginal_price_water()`
// is `marginal_cost_water_multilayer() - marginal_cost_water()`, and the
// multilayer lambda is the one `find_root_collar_psi` EQUALISES -- it is defined
// by the first-order condition. So the analytic soil row is a valid expression
// for dprofit/dpsi_soil AT the operating point and not away from it, and
// differencing it in the collar differentiates something that is not the profit
// row off the optimum.
//
// So d2profit/dp du has to come from the marginal profit itself, and it does:
// `supply_row` above takes it from the condition's own two coefficients, which
// are written in the stem potential and its collar response rather than in total
// uptake and uptake's collar slope. THAT DISTINCTION IS THE WHOLE OF IT -- the
// two pairs span the same directions, so a rank test passes for both, and it was
// a coefficient derived in the second pairing that cost this corpus two reverted
// attempts.

}  // namespace gradient
}  // namespace phylloptim

#endif
