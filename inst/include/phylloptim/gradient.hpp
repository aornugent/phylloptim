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
#include <limits>
#include <string>
#include <vector>

namespace phylloptim {
namespace gradient {

// --- the parameter enumeration, which R indexes into --------------------------
//
// The thirteen traits in `Leaf::set_traits`' argument order, then the two
// quantities a calibration fits that are not traits: the conductance driver and
// the single-potential path's series resistance.
//
// ⚠️ APPENDING IS SAFE AND REORDERING IS NOT. R names these positions in
// `.gradient_par_names` and passes integer indices, so a swap here silently
// differentiates the wrong parameter. `test-gradient-batch.R` reads the names
// back out of C++ and compares them with R's, so the two cannot drift apart
// without a failure.
inline constexpr int n_traits = 13;
inline constexpr int n_pars = 15;

// The two indices the code below has to know by name: one takes the fast
// homogeneity path, and the two non-traits take a relative step.
inline constexpr int par_stem_b = 2;
inline constexpr int par_kmax = 13;
inline constexpr int par_resistance = 14;

inline const std::vector<std::string>& par_names() {
  static const std::vector<std::string> names{
      "vcmax_25",  "stem_c",              "stem_b",
      "psi_crit",  "root_c",              "root_b",
      "root_psi_crit", "beta2",           "jmax_25",
      "a",         "curv_fact_elec_trans", "curv_fact_colim",
      "cost_scale_TF24",
      "leaf_specific_conductance_max",
      "resistance"};
  return names;
}

// --- the environment rows: the two things plants share ------------------------
//
// Everything above differentiates the leaf with respect to what it IS. These rows
// differentiate it with respect to what it EXPERIENCES. Plants in this model are
// coupled through exactly two quantities -- the light they cast on each other and
// the soil water they draw from -- so a gradient that stops at the traits has
// nothing to say about competition: one plant's uptake is the next plant's
// `psi_soil`, one plant's canopy is the next plant's `PPFD`, and neither had a
// row.
//
// THE ENVIRONMENT ROWS ARE NOT IN `theta`, and that is the design decision. Their
// values already have exactly one home -- the observation's `Drivers` -- and
// copying them into `theta` would give one solve two soil states free to
// disagree. So `theta` stays fifteen wide, R's fifteen-column contract and
// `gradient_par_names()` are untouched, and an environment row is addressed by an
// index PAST the end of `theta` whose value is read from the drivers instead
// (`par_value` below). Nothing else about the two routes changes: `psi_soil` and
// `PPFD` are drivers exactly as `leaf_specific_conductance_max` is, so the
// existing perturbation loops differentiate them with no new algebra.
//
// ARITY IS THE NEW THING. Everything in `theta` is a scalar; `psi_soil` is a
// vector of length L, so the row count depends on the layer configuration and is
// not a compile-time constant. It is handled by keeping the FIXED part fixed:
// `PPFD` sits at `par_PPFD` whatever L is, and the L soil rows follow it. An
// index therefore means the same parameter across observations with different
// layer counts, which it would not if the variable-length block came first. The
// soil rows are named `psi_soil_1..psi_soil_L`, one-based to match the layer
// numbering the caller already has, and generated per call rather than cached
// because L is a property of the observation and not of this header.
inline constexpr int par_PPFD = n_pars;
inline constexpr int par_psi_soil_first = n_pars + 1;

// Total rows for a given layer count, and the full row names in index order.
inline constexpr int n_pars_total(int n_layers) {
  return par_psi_soil_first + n_layers;
}

inline std::vector<std::string> par_names(int n_layers) {
  std::vector<std::string> out = par_names();
  out.reserve(std::size_t(n_pars_total(n_layers)));
  out.emplace_back("PPFD");
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("psi_soil_" + std::to_string(i + 1));
  }
  return out;
}

// One name, for a diagnostic message. Does not build the whole vector.
inline std::string par_name(int par, int n_layers) {
  if (par >= 0 && par < n_pars) {
    return par_names()[std::size_t(par)];
  }
  if (par == par_PPFD) {
    return "PPFD";
  }
  const int layer = par - par_psi_soil_first;
  if (layer >= 0 && layer < n_layers) {
    return "psi_soil_" + std::to_string(layer + 1);
  }
  return "parameter " + std::to_string(par);
}

// --- the four differentiated outputs -----------------------------------------
//
// A, gc, psi_stem and collar, in that order, which is R's
// `.gradient_output_names`. `collar` is psi* itself, which is what makes
// dcollar/dtheta equal dpsi*/dtheta and lets the two routes below compute the
// same quantity by different means.
inline constexpr int n_outputs = 4;
inline constexpr int out_collar = 3;

inline const std::vector<std::string>& output_names() {
  static const std::vector<std::string> names{"A", "gc", "psi_stem", "collar"};
  return names;
}

// Read straight off the members rather than through `operating_point_values()`,
// which is what R has to use. Bit-identical: that reader copies these same four
// fields into positions 3, 5, 0 and 1 of its twelve, and the three columns it
// computes rather than copies (uptake, lambda, g1_eff) are not among them.
inline void outputs(const Leaf& l, double* y) {
  y[0] = l.assim_colimited_;
  y[1] = l.stom_cond_CO2_;
  y[2] = l.opt_psi_stem_;
  y[3] = l.opt_root_psi_;
}

// The outputs with the collar held at `psi` rather than optimised. False when
// the clamp moved the target, because then this is not the evaluation that was
// asked for.
//
// ⚠️ EXACT EQUALITY IS THE RIGHT TEST AND THE ONLY ONE THAT WORKS.
// `evaluate_root_collar_psi` CLAMPS its target into the feasible interval, so a
// clamped evaluation is silently a one-sided difference over a shorter interval
// -- the same class of error as differentiating at a pinned optimum, and just as
// plausible-looking. The clamp is a min/max, so an unclamped target comes back
// bit-identical and a tolerance would only blur the detector.
inline bool outputs_at(Leaf& l, double psi, double* y) {
  l.evaluate_root_collar_psi(psi);
  if (!util::identical(l.opt_root_psi_, psi)) {
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
inline double step_for(int par, double value, double step) {
  const double floor =
      (par == par_kmax || par == par_resistance) ? 0.0 : 1.0;
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

// The value a parameter currently holds: in `theta` for the fifteen, in the
// observation's drivers for the environment rows. See the environment block above
// for why those two are different places rather than one.
inline double par_value(const double* theta, const Drivers& d, int par) {
  if (par < n_pars) {
    return theta[par];
  }
  if (par == par_PPFD) {
    return d.PPFD;
  }
  return d.psi_soil[std::size_t(par - par_psi_soil_first)];
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
  double stationarity_tol = 1e-8;
  Method method = Method::Auto;
  bool fast_stem_curve = true;
};

struct Result {
  // The four outputs the gradient is taken at.
  double value[n_outputs];
  // npars * n_outputs, parameter-major: d(output j)/d(pars[k]) at [k * 4 + j].
  std::vector<double> grad;
  Status status = Status::Error;
  bool used_ift = false;
  double H = util::na_value;
  double stationarity = util::na_value;
  std::string message;

  void reset(std::size_t npars) {
    for (int j = 0; j < n_outputs; ++j) {
      value[j] = util::na_value;
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
// `PPFD` and `psi_soil` are passed rather than read off `d`, so that an
// environment perturbation is applied by handing over a different environment
// rather than by editing the observation's drivers in place. Every other caller
// takes the two-argument overload below and gets `d`'s own values.
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve, double PPFD,
                  const std::vector<double>& psi_soil) {
  // THE FAST PATH FOR stem_b, which is the whole of PLAN 11f. The stem
  // cumulative-vulnerability integral is homogeneous of degree 1 in stem_b, so
  // the spline for a perturbed stem_b is the existing one with its argument
  // rescaled and the 11.9 us of incomplete gammas a rebuild spends is
  // unnecessary -- 24.5x on that parameter's gradient.
  //
  // ⚠️ Sound only because `only` names a SINGLE parameter, so everything else in
  // `theta` is still what the object was last set to. That is true here and
  // nowhere else, which is why the argument exists rather than the function
  // guessing. `stem_c` is deliberately not here: it has no such identity, and
  // reading the curve from its closed form instead differentiates a slightly
  // different model and disagrees by 3e-4 (PLAN 11f).
  if (fast_stem_curve && only == par_stem_b) {
    l.perturb_stem_b(theta[par_stem_b]);
    return;
  }
  l.set_traits(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5],
               theta[6], theta[7], theta[8], theta[9], theta[10], theta[11],
               theta[12]);
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
    l.set_physiology(d.root_network, PPFD, psi_soil, d.soil_depth,
                     theta[par_kmax], d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  }
}

// At the observation's own environment, which is what every parameter except the
// environment rows is differentiated in.
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve) {
  apply(l, theta, d, single, only, fast_stem_curve, d.PPFD, d.psi_soil);
}

// Put parameter `par` at `value` and everything else at base. One function so
// that both routes below address the fifteen and the environment rows the same
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
//   ⚠️ LAYERS BELOW THE DEEPEST ROOTED ONE GET EXACTLY ZERO UPTAKE, so their four
//   rows are exactly 0.0 on both routes. That zero is CORRECT and is asserted
//   (`environment rows: an unrooted layer's rows are exactly zero`). It is called
//   out because in this codebase an exact zero is otherwise the signature of a
//   missing accumulator -- and equally so that a genuinely missing row is not
//   read as this.
//
// `psi_scratch` is the caller's buffer rather than a local so that the copy
// allocates once per gradient instead of once per perturbation.
inline void set_one(Leaf& l, double* th, const double* theta, const Drivers& d,
                    bool single, int par, double value, bool fast_stem_curve,
                    std::vector<double>& psi_scratch) {
  std::copy(theta, theta + n_pars, th);
  if (par < n_pars) {
    th[par] = value;
    apply(l, th, d, single, par, fast_stem_curve);
    return;
  }
  if (par == par_PPFD) {
    apply(l, th, d, single, par, fast_stem_curve, value, d.psi_soil);
    return;
  }
  psi_scratch = d.psi_soil;
  psi_scratch[std::size_t(par - par_psi_soil_first)] = value;
  apply(l, th, d, single, par, fast_stem_curve, d.PPFD, psi_scratch);
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

// The implicit-function composite. Two perturbed evaluations per parameter,
// neither of which re-solves the model: `dprofit` at the UNPERTURBED psi* gives
// the mixed partial, and the outputs at that same psi* give the direct term.
inline void gradient_ift(Leaf& l, const double* theta, const Drivers& d,
                         bool single, const int* pars, std::size_t npars,
                         double psi_star, double H, const double* dY_dpsi,
                         const Settings& s, double* out) {
  double th[n_pars];
  double up[1 + n_outputs];
  double dn[1 + n_outputs];
  std::vector<double> psi_scratch;
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    const int p = pars[k];
    if (takes_shortcut(p, s) && !at_base) {
      apply(l, theta, d, single, -1, s.fast_stem_curve);
    }
    at_base = false;
    const double base = par_value(theta, d, p);
    const double h = step_for(p, base, s.step);
    for (int side = 0; side < 2; ++side) {
      // Up first, then down: R evaluates `up <- side(1)` before `dn <- side(-1)`
      // and both mutate the leaf.
      set_one(l, th, theta, d, single, p, side == 0 ? base + h : base - h,
              s.fast_stem_curve, psi_scratch);
      double* dst = side == 0 ? up : dn;
      // Evaluate first, then read dprofit at the same fixed collar -- R's order,
      // and `evaluate_root_collar_psi` is what seats the state `dprofit` reads.
      if (!outputs_at(l, psi_star, dst + 1)) {
        util::stop("leaf_gradient(): perturbing `" +
                   par_name(p, n_soil_layers(d, single)) +
                   "` moved the feasible collar interval past psi*, so the "
                   "operating point could not be evaluated there. This point is "
                   "on an active-set boundary; lower `stationarity_tol` or "
                   "difference the solve directly.");
      }
      dst[0] = l.dprofit_droot_collar_psi(psi_star);
    }
    // M = d2profit/dpsi dtheta, with psi held FIXED at psi*.
    const double dpsi_dtheta = -((up[0] - dn[0]) / (2.0 * h)) / H;
    for (int j = 0; j < n_outputs; ++j) {
      const double direct = (up[1 + j] - dn[1 + j]) / (2.0 * h);
      out[k * n_outputs + j] = direct + rounded(dY_dpsi[j] * dpsi_dtheta);
    }
    // `collar` is not an output of the evaluation -- it IS psi*, held fixed, so
    // its direct term is zero by construction and the composite reduces to
    // dpsi*/dtheta. Set explicitly rather than left as the difference of two
    // identical numbers.
    out[k * n_outputs + out_collar] = dpsi_dtheta;
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);
}

// The fallback: a central difference of the whole solve. Correct at a pinned
// optimum because it differences the CONSTRAINED answer, which is exactly what
// the composite cannot do.
inline void gradient_fd(Leaf& l, const double* theta, const Drivers& d,
                        bool single, const int* pars, std::size_t npars,
                        const Settings& s, double* out) {
  double th[n_pars];
  double up[n_outputs];
  double dn[n_outputs];
  std::vector<double> psi_scratch;
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    const int p = pars[k];
    if (takes_shortcut(p, s) && !at_base) {
      apply(l, theta, d, single, -1, s.fast_stem_curve);
    }
    at_base = false;
    const double base = par_value(theta, d, p);
    const double h = step_for(p, base, s.step);
    for (int side = 0; side < 2; ++side) {
      set_one(l, th, theta, d, single, p, side == 0 ? base + h : base - h,
              s.fast_stem_curve, psi_scratch);
      l.find_root_collar_psi();
      outputs(l, side == 0 ? up : dn);
    }
    for (int j = 0; j < n_outputs; ++j) {
      out[k * n_outputs + j] = (up[j] - dn[j]) / (2.0 * h);
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

  // ⚠️ CHECKED HERE, because with the environment rows the valid range depends on
  // the OBSERVATION rather than on this header: a five-layer row and a one-layer
  // row in the same batch do not have the same number of parameters. An index
  // past the end used to be impossible (R validates against the fixed fifteen);
  // now it is an out-of-bounds read of `psi_soil`, so it is a per-row error.
  const int n_layers = n_soil_layers(d, single);
  for (std::size_t k = 0; k < npars; ++k) {
    if (pars[k] < 0 || pars[k] >= n_pars_total(n_layers)) {
      util::stop("leaf_gradient(): parameter index " + std::to_string(pars[k]) +
                 " is out of range; this observation has " +
                 std::to_string(n_layers) + " soil layer(s), so there are " +
                 std::to_string(n_pars_total(n_layers)) + " parameters.");
    }
  }

  apply(l, theta, d, single, -1, s.fast_stem_curve);
  l.find_root_collar_psi();

  const double psi_star = l.opt_root_psi_;
  outputs(l, out.value);

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
  const double h_psi = std::max(std::abs(psi_star), 1.0) * s.step;
  const double resid = l.dprofit_droot_collar_psi(psi_star);
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
  // feasibility boundary where the centred one jumps, and it leaves `status`
  // unchanged everywhere. Both arms bad: there is no curvature to report and the
  // 0.0 below is the shut-down signature the next comment describes.
  const bool hi_sentinel = util::identical(d_hi, 0.0);
  const bool lo_sentinel = util::identical(d_lo, 0.0);
  const double H =
      (hi_sentinel && lo_sentinel) ? 0.0
      : hi_sentinel                ? (resid - d_lo) / h_psi
      : lo_sentinel                ? (d_hi - resid) / h_psi
                                   : (d_hi - d_lo) / (2.0 * h_psi);
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

  double dY_dpsi[n_outputs];
  if (use_ift) {
    // dY/dpsi at fixed traits, and a SECOND, INDEPENDENT detector of a pinned
    // optimum. At a pinned point psi* sits one step-in fraction (1e-06 of the
    // bracket width) from its bound, so a step of `step * psi` crosses it
    // whenever the bracket is narrower than psi -- which every pinned row in
    // this package's grid is. Measured, that catches all 42 pinned rows and all
    // 48 shut-down ones on its own.
    //
    // It is NOT a substitute for the stationarity test: it fires only when the
    // bracket is narrow, so a pinned optimum on a wide bracket would pass it.
    double hi[n_outputs];
    double lo[n_outputs];
    // Both, unconditionally, before the test -- R computes `hi` and `lo` on
    // consecutive lines and only then checks either, and each call moves the
    // leaf.
    const bool hi_ok = outputs_at(l, psi_star + h_psi, hi);
    const bool lo_ok = outputs_at(l, psi_star - h_psi, lo);
    if (!hi_ok || !lo_ok) {
      if (s.method == Method::Ift) {
        util::stop("leaf_gradient(): method = \"ift\" was asked for at a point "
                   "whose feasible collar interval is narrower than one step, "
                   "so dY/dpsi cannot be centred on psi*. Use "
                   "method = \"auto\".");
      }
      use_ift = false;
      out.status = Status::Pinned;
    } else {
      for (int j = 0; j < n_outputs; ++j) {
        dY_dpsi[j] = (hi[j] - lo[j]) / (2.0 * h_psi);
      }
    }
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

}  // namespace gradient
}  // namespace phylloptim

#endif
