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
// The fourteen traits in `Leaf::set_traits`' argument order, then the two
// quantities a calibration fits that are not traits: the conductance driver and
// the single-potential path's series resistance.
//
// ⚠️ R INDEXES THESE POSITIONS, so a reordering silently differentiates the wrong
// parameter. `test-gradient-batch.R` reads the names back out of C++ and compares
// them with R's, so the two cannot drift apart without a failure.
inline constexpr int n_traits = 14;
inline constexpr int n_pars = 16;

// Every index by name, so nothing below indexes `theta` with a bare integer.
// The first `n_traits` are `set_traits`' arguments in its order, which is also
// `leaf_traits()`'; the two non-traits follow and take a relative step.
inline constexpr int par_vcmax_25 = 0;
inline constexpr int par_stem_c = 1;
inline constexpr int par_stem_b = 2;
inline constexpr int par_psi_crit = 3;
inline constexpr int par_root_c = 4;
inline constexpr int par_root_b = 5;
inline constexpr int par_root_psi_crit = 6;
inline constexpr int par_beta2 = 7;
inline constexpr int par_jmax_25 = 8;
inline constexpr int par_a = 9;
inline constexpr int par_curv_fact_elec_trans = 10;
inline constexpr int par_curv_fact_colim = 11;
inline constexpr int par_cost_scale_TF24 = 12;
inline constexpr int par_R_d_25 = 13;
inline constexpr int par_kmax = 14;
inline constexpr int par_resistance = 15;

inline const std::vector<std::string>& par_names() {
  static const std::vector<std::string> names{
      "vcmax_25",  "stem_c",              "stem_b",
      "psi_crit",  "root_c",              "root_b",
      "root_psi_crit", "beta2",           "jmax_25",
      "a",         "curv_fact_elec_trans", "curv_fact_colim",
      "cost_scale_TF24", "R_d_25",
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
// disagree. So `theta` stays sixteen wide, R's sixteen-column contract and
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

// One root-carbon row per layer, after the soil block: a consumer's input vector
// carries two entries per layer, so the arity is n_pars + 1 + 2L. With two
// variable-length blocks an index past `PPFD` names one input only at a fixed L.
inline constexpr int par_root_carbon_first(int n_layers) {
  return par_psi_soil_first + n_layers;
}

// Total rows for a given layer count, and the full row names in index order.
inline constexpr int n_pars_total(int n_layers) {
  return par_root_carbon_first(n_layers) + n_layers;
}

inline std::vector<std::string> par_names(int n_layers) {
  std::vector<std::string> out = par_names();
  out.reserve(std::size_t(n_pars_total(n_layers)));
  out.emplace_back("PPFD");
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("psi_soil_" + std::to_string(i + 1));
  }
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("root_carbon_" + std::to_string(i + 1));
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
  const int carbon = par - par_root_carbon_first(n_layers);
  if (carbon >= 0 && carbon < n_layers) {
    return "root_carbon_" + std::to_string(carbon + 1);
  }
  return "parameter " + std::to_string(par);
}

// --- the five differentiated outputs -----------------------------------------
//
// A, gc, psi_stem, collar and profit, in that order. THIS LIST IS THE ONE
// DEFINITION OF IT: R does not keep a second copy, it reads this one back through
// `gradient_output_names()` in src/gradient.cpp, so the two routes cannot
// disagree about what they are reporting. `profit` is appended rather than
// inserted, because R hands `pars` over as positions and a caller reads gradient
// columns by position after that.
//
// Two of the five are not the generic composite in `gradient_ift` below, and for
// two DIFFERENT reasons, which is why both are named here:
//
//   * `collar` is psi* itself, which is what makes dcollar/dtheta equal
//     dpsi*/dtheta and lets the two routes compute the same quantity by
//     different means.
//   * `profit` is the OBJECTIVE rather than an output read at the argmax, which
//     is what brings the envelope theorem into play.
inline constexpr int n_outputs = 5;
inline constexpr int out_collar = 3;
inline constexpr int out_profit = 4;

inline const std::vector<std::string>& output_names() {
  static const std::vector<std::string> names{"A", "gc", "psi_stem", "collar",
                                              "profit"};
  return names;
}

// --- the per-layer uptake outputs ---------------------------------------------
//
// The five above are a calibration's, and R reads gradient columns by position,
// so they stay fixed: `output_names()` and `gradient_output_names()` report those
// five whatever the layer configuration is, and `at`, `batch` and `transpose_at`
// report those five. A stand adjoint reads a different list -- of the five it
// wants profit, and it wants the water each soil layer gave up -- so the uptake
// entries are APPENDED, one per layer, and `rows_at` is the one route that can
// request them.
//
// ⚠️ AN OUTPUT THAT IS NOT REPORTED HAS NO ROUTE TO ANY INPUT. At a pin the
// condition's gradient is zero for every input and each reported output carries a
// total instead, so an unreported one comes back exactly zero at every pinned and
// shut point -- and zero and absent are the same number on a consumer's tape.
// That is why uptake is an enumeration entry rather than a second call.
//
// Uptake is ORDINARY: it is set as a side effect at the operating point, so it
// consumes the collar rather than being it, and it composes the way A and gc do.
// Fixed first, variable second, for the reason the soil rows follow `PPFD` on the
// input side: an index means the same output across observations with different
// layer counts.
inline constexpr int out_uptake_first = n_outputs;

inline constexpr int n_outputs_total(int n_layers) {
  return out_uptake_first + n_layers;
}

inline std::vector<std::string> output_names(int n_layers) {
  std::vector<std::string> out = output_names();
  out.reserve(std::size_t(n_outputs_total(n_layers)));
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("uptake_" + std::to_string(i + 1));
  }
  return out;
}

// One name, for a diagnostic message. Does not build the whole vector.
inline std::string output_name(int out, int n_layers) {
  if (out >= 0 && out < n_outputs) {
    return output_names()[std::size_t(out)];
  }
  const int layer = out - out_uptake_first;
  if (layer >= 0 && layer < n_layers) {
    return "uptake_" + std::to_string(layer + 1);
  }
  return "output " + std::to_string(out);
}

// One evaluation's outputs: the five, then one uptake per layer where the caller
// asked for them. A fixed-five caller leaves `uptake` empty, which allocates
// nothing.
struct OutputValues {
  double fixed[n_outputs];
  std::vector<double> uptake;

  OutputValues() = default;
  explicit OutputValues(int n_uptake) : uptake(std::size_t(n_uptake)) {}

  int n_uptake() const { return int(uptake.size()); }
  int size() const { return n_outputs + n_uptake(); }
  double& operator[](int j) {
    return j < n_outputs ? fixed[j] : uptake[std::size_t(j - out_uptake_first)];
  }
  double operator[](int j) const {
    return j < n_outputs ? fixed[j] : uptake[std::size_t(j - out_uptake_first)];
  }
};

// Read straight off the members rather than through `operating_point_values()`,
// which is what R has to use. Bit-identical: that reader copies these same five
// fields into positions 3, 5, 0, 1 and 6 of its twelve, and the three columns it
// computes rather than copies (uptake, lambda, g1_eff) are not among them.
//
// The uptake entries are the leaf's own per-layer consumption at the collar this
// evaluation seated -- `find_psi_stem_from_psi_root` writes it on its way past --
// and not a second computation of it.
inline void outputs(const Leaf& l, OutputValues& y) {
  y[0] = l.assim_colimited_;
  y[1] = l.stom_cond_CO2_;
  y[2] = l.opt_psi_stem_;
  y[3] = l.opt_root_psi_;
  y[4] = l.profit_;
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
inline bool outputs_at(Leaf& l, double psi, OutputValues& y) {
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

// The value a parameter currently holds: in `theta` for the sixteen, in the
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

// Every requested index names an input this observation has and this package can
// move.
//
// ⚠️ CHECKED PER CALL, because with the environment rows the valid range depends
// on the OBSERVATION rather than on this header: a five-layer row and a one-layer
// row in the same batch do not have the same number of parameters. An index past
// the end used to be impossible (R validates against the fixed sixteen); now it
// is an out-of-bounds read of `psi_soil`, so it is a per-row error.
//
// ⚠️ ROOT CARBON HAS NO ROW ON EITHER ROUTE. It is in the enumeration because a
// consumer's input vector carries it, but `set_physiology` is handed a
// `RootNetwork` that is already built, so nothing here can move the carbon it was
// built from. Refused by name rather than left to read past `psi_soil`.
inline void check_pars(const int* pars, std::size_t npars, int n_layers,
                       const std::string& caller) {
  for (std::size_t k = 0; k < npars; ++k) {
    if (pars[k] < 0 || pars[k] >= n_pars_total(n_layers)) {
      util::stop(caller + ": parameter index " + std::to_string(pars[k]) +
                 " is out of range; this observation has " +
                 std::to_string(n_layers) + " soil layer(s), so there are " +
                 std::to_string(n_pars_total(n_layers)) + " parameters.");
    }
    if (pars[k] >= par_root_carbon_first(n_layers)) {
      util::stop(caller + ": `" + par_name(pars[k], n_layers) +
                 "` has no row here. The root network arrives already built, so "
                 "nothing in this package can move the carbon behind it.");
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
  l.set_traits(theta[par_vcmax_25], theta[par_stem_c], theta[par_stem_b],
               theta[par_psi_crit], theta[par_root_c], theta[par_root_b],
               theta[par_root_psi_crit], theta[par_beta2], theta[par_jmax_25],
               theta[par_a], theta[par_curv_fact_elec_trans],
               theta[par_curv_fact_colim], theta[par_cost_scale_TF24],
               theta[par_R_d_25]);
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

// The solve, the five outputs, the marginal profit at the collar the solve
// returned, and the curvature of profit there.
struct BasePoint {
  double psi_star = util::na_value;
  OutputValues value;
  double resid = util::na_value;
  double H = util::na_value;
  // ⚠️ Read the moment the solve ends, because the first evaluation at a held
  // collar overwrites it: `evaluate_root_collar_psi` tags the point Prescribed.
  Branch branch;
};

// The step in the collar potential, floored at 1 MPa for `step_for`'s reason.
inline double collar_step(double psi_star, const Settings& s) {
  return std::max(std::abs(psi_star), 1.0) * s.step;
}

inline BasePoint base_point(Leaf& l, const double* theta, const Drivers& d,
                            bool single, const Settings& s) {
  apply(l, theta, d, single, -1, s.fast_stem_curve);
  l.find_root_collar_psi();

  BasePoint b;
  b.branch = branch_here(l);
  b.psi_star = l.opt_root_psi_;
  outputs(l, b.value);

  const double h_psi = collar_step(b.psi_star, s);
  b.resid = l.dprofit_droot_collar_psi(b.psi_star);
  // Named halves: `f(a) - f(b)` has unspecified operand order in C++ and
  // left-to-right order in R, and `dprofit_droot_collar_psi` mutates the leaf.
  const double d_hi = l.dprofit_droot_collar_psi(b.psi_star + h_psi);
  const double d_lo = l.dprofit_droot_collar_psi(b.psi_star - h_psi);
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
  b.H = (hi_sentinel && lo_sentinel) ? 0.0
        : hi_sentinel                ? (b.resid - d_lo) / h_psi
        : lo_sentinel                ? (d_hi - b.resid) / h_psi
                                     : (d_hi - d_lo) / (2.0 * h_psi);
  return b;
}

// dY/dpsi at fixed traits, and a SECOND, INDEPENDENT detector of a pinned
// optimum. At a pinned point psi* sits one step-in fraction (1e-06 of the
// bracket width) from its bound, so a step of `step * psi` crosses it whenever
// the bracket is narrower than psi -- which every pinned row in this package's
// grid is. Measured, that catches all 42 pinned rows and all 48 shut-down ones
// on its own.
//
// It is NOT a substitute for the stationarity test: it fires only when the
// bracket is narrow, so a pinned optimum on a wide bracket would pass it. False
// where the difference cannot be centred on psi*.
inline bool collar_channel(Leaf& l, double psi_star, const Settings& s,
                           OutputValues& dY_dpsi) {
  const double h_psi = collar_step(psi_star, s);
  OutputValues hi(dY_dpsi.n_uptake());
  OutputValues lo(dY_dpsi.n_uptake());
  // Both, unconditionally, before the test -- R computes `hi` and `lo` on
  // consecutive lines and only then checks either, and each call moves the leaf.
  const bool hi_ok = outputs_at(l, psi_star + h_psi, hi);
  const bool lo_ok = outputs_at(l, psi_star - h_psi, lo);
  if (!hi_ok || !lo_ok) {
    return false;
  }
  for (int j = 0; j < dY_dpsi.size(); ++j) {
    dY_dpsi[j] = (hi[j] - lo[j]) / (2.0 * h_psi);
  }
  return true;
}

// One input's two perturbed evaluations at a HELD collar: the requested outputs'
// direct rows, and dR/du -- the collar derivative of marginal profit. Neither
// evaluation re-solves the model. How many outputs is `direct`'s own length.
//
// `at_base` carries the invariant above: it comes in true only where the leaf is
// at base parameters, and goes out false.
inline void held_row(Leaf& l, const double* theta, const Drivers& d,
                       bool single, int par, double psi_star, const Settings& s,
                       const std::string& caller, bool& at_base,
                       std::vector<double>& psi_scratch, OutputValues& direct,
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
  const double base = par_value(theta, d, par);
  const double h = step_for(par, base, s.step);
  for (int side = 0; side < 2; ++side) {
    // Up first, then down: R evaluates `up <- side(1)` before `dn <- side(-1)`
    // and both mutate the leaf.
    set_one(l, th, theta, d, single, par, side == 0 ? base + h : base - h,
            s.fast_stem_curve, psi_scratch);
    OutputValues& dst = side == 0 ? up : dn;
    // Evaluate first, then read dprofit at the same fixed collar -- R's order,
    // and `evaluate_root_collar_psi` is what seats the state `dprofit` reads.
    if (!outputs_at(l, psi_star, dst)) {
      util::stop(caller + ": perturbing `" +
                 par_name(par, n_soil_layers(d, single)) +
                 "` moved the feasible collar interval past psi*, so the "
                 "operating point could not be evaluated there. This point is "
                 "on an active-set boundary; lower `stationarity_tol` or "
                 "difference the solve directly.");
    }
    double& resid = side == 0 ? up_resid : dn_resid;
    resid = l.dprofit_droot_collar_psi(psi_star);
  }
  // M = d2profit/dpsi dtheta, with psi held FIXED at psi*.
  dresidual = (up_resid - dn_resid) / (2.0 * h);
  for (int j = 0; j < direct.size(); ++j) {
    direct[j] = (up[j] - dn[j]) / (2.0 * h);
  }
}

// One input's central difference of the WHOLE solve. Correct at a pinned optimum
// because it differences the CONSTRAINED answer, which is exactly what the
// composite cannot do.
//
// `stay` is the branch both arms have to land on, or null for a caller that
// differences whatever they land on. Where it is given, the step shrinks by a
// decade at a time looking for two arms that agree with it, and the row is refused
// rather than taken across a change of branch. Two decades and no more: the
// largest step that is still local is the best-conditioned one, and below that
// floor a step stops moving the solve above its own noise, so the difference
// measures rounding rather than the model.
//
// False only where `stay` was given and no step reached it; `row` is not written
// then.
inline bool solved_row(Leaf& l, const double* theta, const Drivers& d,
                       bool single, int par, const Settings& s,
                       const Branch* stay, bool& at_base,
                       std::vector<double>& psi_scratch, OutputValues& row) {
  if (takes_shortcut(par, s) && !at_base) {
    apply(l, theta, d, single, -1, s.fast_stem_curve);
  }
  at_base = false;
  double th[n_pars];
  OutputValues up(row.n_uptake());
  OutputValues dn(row.n_uptake());
  const double base = par_value(theta, d, par);
  double h = step_for(par, base, s.step);
  for (int decade = 0;; ++decade) {
    bool on_branch = true;
    for (int side = 0; side < 2; ++side) {
      set_one(l, th, theta, d, single, par, side == 0 ? base + h : base - h,
              s.fast_stem_curve, psi_scratch);
      l.find_root_collar_psi();
      // Both arms are taken before either is tested, so which one crossed does
      // not change the order the leaf is moved in.
      on_branch = on_branch && (stay == nullptr || branch_here(l) == *stay);
      outputs(l, side == 0 ? up : dn);
    }
    if (on_branch) {
      for (int j = 0; j < row.size(); ++j) {
        row[j] = (up[j] - dn[j]) / (2.0 * h);
      }
      return true;
    }
    if (decade == 2) {
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
  std::vector<double> psi_scratch;
  OutputValues direct;
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    double M = 0.0;
    held_row(l, theta, d, single, pars[k], psi_star, s, "leaf_gradient()",
               at_base, psi_scratch, direct, M);
    const double dpsi_dtheta = -(M / H);
    for (int j = 0; j < n_outputs; ++j) {
      out[k * n_outputs + j] = direct[j] + rounded(dY_dpsi[j] * dpsi_dtheta);
    }
    // TWO OF THE FIVE ARE NOT THAT COMPOSITE, AND THE TWO ARGUMENTS ARE
    // DIFFERENT. Both overwrite the generic line rather than skipping it, which
    // is R's shape too: `.gradient_ift` computes `g` for all five and then
    // assigns these two.
    //
    // `collar` is not an output of the evaluation -- it IS psi*, held fixed, so
    // its DIRECT term is zero by construction and the composite reduces to
    // dpsi*/dtheta. Set explicitly rather than left as the difference of two
    // identical numbers.
    out[k * n_outputs + out_collar] = dpsi_dtheta;
    // `profit` is the objective rather than an output read at the argmax, so for
    // it the OTHER factor of the second term is the one that vanishes:
    //
    //   dprofit*/dtheta = (dprofit/dpsi)(dpsi*/dtheta) + dprofit/dtheta|_psi
    //
    // and dprofit/dpsi = 0 at an interior optimum. That is the envelope theorem,
    // and its premise is exactly the stationarity `at()` has already tested. So
    // the direct term IS the answer: no M, no H, no dpsi*/dtheta.
    //
    // ⚠️ AND THE ARITHMETIC NEEDS THIS LINE EVEN SO, which is a second reason and
    // not a restatement of the first. `dY_dpsi[out_profit]` is a central
    // difference of a FLAT maximum: ~1e-09 of cancellation in profit over a
    // ~1e-06 step in psi. Measured at psi_soil = 2.0 the exact dprofit at psi* is
    // 8.7e-11, while that difference reads -1.8e-04 -- so leaving the generic
    // line to stand would move this column by up to 6.6e-05 relative. That is
    // five orders above the solve's ~1e-09 floor and entirely plausible-looking.
    //
    // ⚠️ THIS ROUTE ONLY. At a pinned optimum dprofit/dpsi is NOT zero, the
    // envelope theorem does not hold, and the profit row has to come from
    // differencing the solve -- which is what `gradient_fd` already does for
    // every output, with no special case at all. `at()`'s `status` chooses
    // between the two functions and is the only place that decision is made, so
    // there is no way for this line to reach the fallback.
    out[k * n_outputs + out_profit] = direct[out_profit];
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
  std::vector<double> psi_scratch;
  OutputValues row;
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    // No branch to hold: this route differences whatever the two arms land on and
    // reports the base point's own kind beside the answer, so it cannot refuse.
    // The copy is what leaves a refused row NA rather than partly written.
    if (solved_row(l, theta, d, single, pars[k], s, nullptr, at_base,
                   psi_scratch, row)) {
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
  check_pars(pars, npars, n_soil_layers(d, single), "leaf_gradient()");

  const BasePoint b = base_point(l, theta, d, single, s);
  const double psi_star = b.psi_star;
  const double resid = b.resid;
  const double H = b.H;
  for (int j = 0; j < n_outputs; ++j) {
    out.value[j] = b.value[j];
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

  // Written straight into the result rather than into a local, so that the
  // transpose can read the same numbers instead of measuring them again.
  OutputValues& dY_dpsi = out.dY_dpsi;
  if (use_ift && !collar_channel(l, psi_star, s, dY_dpsi)) {
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
// ⚠️ AN INPUT THAT MOVES THE BOUND A CONSTRAINED POINT SITS ON GETS ITS ROWS BY
// FOLLOWING IT, and they come back in the same two fields. p* sits one step-in
// fraction from the bound, so a step that moves the bound carries p* out of the
// perturbed feasible interval and there is no held evaluation to take: the whole
// solve is differenced instead, the arms land on the moved point, and what the
// rows carry is the TOTAL. `dresidual` is zero for such an input, so the
// consumer's one assembly returns that total unchanged and nothing subtracts.
//
// THIS IS THE ONE ROUTE THAT CAN BE ASKED FOR PER-LAYER UPTAKE. `output` indexes
// the enumeration above, so an index at or past `out_uptake_first` names a layer's
// consumption, and an index past the last layer this observation has is refused by
// name. The five stay reportable, and the consumer this exists for asks for
// profit and the uptake block.

using OperatingPointKind = Leaf::OperatingPointKind;

// How the operating point reaches an output. At most one Objective and at most
// one Point; a request declaring two of either is refused.
enum class Role { Objective, Point, Ordinary };

// Why an entry is zero, where it is. A bare zero cannot say, and an exact zero is
// the signature of a missing row more often than of true insensitivity.
enum class Zero { none, slack, structural };

struct RowRequest {
  const int* output;  const Role* role;  std::size_t n_output;
  const int* input;   std::size_t n_input;
};

struct Rows {
  // ⚠️ THE BRANCH THE SOLVE TOOK, never a reading of the numbers. `Status` above
  // derives its own classification from the curvature's sign and the residual's
  // size, which cannot separate a stationary point from the hard 0.0 the no-flow
  // state returns: that reads as stationary, and the curvature taken off the same
  // sentinel confirms it.
  OperatingPointKind kind = OperatingPointKind::Unsolved;
  // Set where the kind is one that has no rows to give, and where an input's own
  // row was refused. It names the input in the second case.
  std::string message;

  double point = util::na_value;
  // R_p at an interior point, the bound's own slope at a pin, NA where no
  // condition defines the point.
  double residual_slope = util::na_value;
  // n_input: grad of the condition that defines p*, and zero for an input whose
  // rows followed the point rather than being taken at a held one.
  std::vector<double> dresidual;
  std::vector<double> dy_dp;      // n_output
  // n_output * n_input, output-major: the row at a held point, or the total where
  // that input's rows followed the point. NA where the input's row was refused.
  std::vector<double> held;
  std::vector<Zero> zero;         // n_output * n_input, output-major
  double amplification = util::na_value;
};

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

// The bound's movement in one input, as `bound_row` reports it -- already divided
// by the bound's slope. Zero for an input neither residual reads: no bound reads
// radiation or any photosynthetic trait, and the wet bound is total uptake, so
// the stem does not enter it.
inline double bound_dpoint(const Leaf::BoundRow& b, int par, int n_layers) {
  switch (par) {
  case par_psi_crit:      return b.d_dpsi_crit;
  case par_root_psi_crit: return b.d_droot_psi_crit;
  case par_stem_b:        return b.d_dstem_b;
  case par_root_b:        return b.d_droot_b;
  case par_kmax:          return b.d_dkappa;
  default:                break;
  }
  const int layer = par - par_psi_soil_first;
  if (layer >= 0 && layer < n_layers) {
    return b.d_dpsi_soil[std::size_t(layer)];
  }
  return 0.0;
}

// The inputs a bound moves with and `bound_row` carries no entry for: the two
// curve traits that RESHAPE a vulnerability grid rather than scaling it, and the
// single path's series resistance.
inline bool bound_moves_by_re_solving(int par, bool single) {
  return par == par_stem_c || par == par_root_c ||
         (single && par == par_resistance);
}

// The bound's movement in one of those, by re-solving it at a perturbed state.
// The forward model rebuilds the grid when a curve trait moves, so the motion is
// the model and a row taken on a held grid differentiates a different function.
inline double differenced_bound(Leaf& l, const double* theta, const Drivers& d,
                                bool single, int par, Leaf::WhichBound bound,
                                const Settings& s, bool& at_base,
                                std::vector<double>& psi_scratch) {
  at_base = false;
  double th[n_pars];
  const double base = par_value(theta, d, par);
  const double h = step_for(par, base, s.step);
  set_one(l, th, theta, d, single, par, base + h, s.fast_stem_curve, psi_scratch);
  const double up = l.bound_row(bound).bound;
  set_one(l, th, theta, d, single, par, base - h, s.fast_stem_curve, psi_scratch);
  const double dn = l.bound_row(bound).bound;
  return (up - dn) / (2.0 * h);
}

// An exact zero, classified where it is written. The root's own critical
// potential is read by nothing but the dry bound, so at an interior optimum its
// zero is complementary slackness; the stem's psi_crit is not the same case,
// since the flux integrates the stem curve up to it. Every other exact zero here
// is an input the branch taken reads through nothing.
inline Zero zero_kind(double entry, int par, bool interior) {
  if (entry != 0.0) {
    return Zero::none;
  }
  if (interior && par == par_root_psi_crit) {
    return Zero::slack;
  }
  return Zero::structural;
}

inline Rows rows_at(Leaf& l, const double* theta, const Drivers& d,
                    const RowRequest& r, const Settings& s) {
  // Which supply path is in force is the leaf's own state, so it is read here
  // rather than passed: an argument beside it could disagree with it.
  const bool single = l.supply_kind_ == Leaf::SupplyKind::SinglePotential;
  const int n_layers = n_soil_layers(d, single);
  check_pars(r.input, r.n_input, n_layers, "leaf_rows()");

  int objective = -1;
  int point = -1;
  // One uptake entry per layer, or none: a single evaluation writes the whole
  // consumption profile, so a request naming any layer carries all of them.
  int n_uptake = 0;
  for (std::size_t j = 0; j < r.n_output; ++j) {
    if (r.output[j] < 0 || r.output[j] >= n_outputs_total(n_layers)) {
      // Named where it can be, because an index past the end of the uptake block
      // is a layer count that does not match the drivers rather than a typo.
      const std::string what =
          r.output[j] >= out_uptake_first
              ? "`uptake_" + std::to_string(r.output[j] - out_uptake_first + 1) +
                    "`"
              : "output index " + std::to_string(r.output[j]);
      util::stop("leaf_rows(): " + what +
                 " is not an output here; this observation has " +
                 std::to_string(n_layers) + " soil layer(s), so there are " +
                 std::to_string(n_outputs_total(n_layers)) +
                 " outputs, the last of them `" +
                 output_name(n_outputs_total(n_layers) - 1, n_layers) + "`.");
    }
    if (r.output[j] >= out_uptake_first) {
      n_uptake = n_layers;
    }
    if (r.role[j] == Role::Objective) {
      if (objective >= 0) {
        util::stop("leaf_rows(): two outputs are declared Objective, and the "
                   "operating point maximises one.");
      }
      objective = int(j);
    } else if (r.role[j] == Role::Point) {
      if (point >= 0) {
        util::stop("leaf_rows(): two outputs are declared Point, and the "
                   "operating point is one number.");
      }
      point = int(j);
    }
  }

  Rows out;
  out.dresidual.assign(r.n_input, util::na_value);
  out.dy_dp.assign(r.n_output, util::na_value);
  out.held.assign(r.n_output * r.n_input, util::na_value);
  out.zero.assign(r.n_output * r.n_input, Zero::none);

  const BasePoint b = base_point(l, theta, d, single, s);
  out.kind = b.branch.kind;
  out.point = b.psi_star;

  Leaf::WhichBound bound = Leaf::WhichBound::Wet;
  const bool pinned = pinned_bound(b.branch.kind, bound);
  const bool interior = b.branch.kind == OperatingPointKind::Interior;
  const bool shut = b.branch.kind == OperatingPointKind::HydraulicShutdown;
  // Hydraulic shutdown holds the stem at psi_crit and moves no water, so no
  // condition defines a collar to differentiate and every row below is a total.
  // Every other kind either never solved or could not choose, shade death
  // included: it seats both potentials at the wet bound, where the objective's
  // sensitivity to the collar is the cost's slope and not the `nu` a pin reads.
  if (!interior && !pinned && !shut) {
    out.message = std::string("no rows at an operating point that is ") +
                  Leaf::operating_point_kind_name(b.branch.kind);
    apply(l, theta, d, single, -1, s.fast_stem_curve);
    return out;
  }

  Leaf::BoundRow condition;
  if (pinned) {
    condition = l.bound_row(bound);
    if (!condition.finite) {
      out.message = "the bound this point is pinned to has no derivative here";
      apply(l, theta, d, single, -1, s.fast_stem_curve);
      return out;
    }
    out.residual_slope = condition.residual_slope;
  } else if (interior) {
    out.residual_slope = b.H;
  }

  OutputValues dY_dpsi(n_uptake);
  const bool have_channel = collar_channel(l, b.psi_star, s, dY_dpsi);
  for (std::size_t j = 0; j < r.n_output; ++j) {
    if (int(j) == point) {
      out.dy_dp[j] = 1.0;
    } else if (int(j) == objective) {
      // Zero at an interior optimum by the envelope theorem; at a pin the point
      // is the bound rather than a maximum, and this is the constraint's shadow
      // price.
      out.dy_dp[j] = interior ? 0.0 : pinned ? b.resid : util::na_value;
    } else if (have_channel) {
      out.dy_dp[j] = dY_dpsi[r.output[j]];
    }
  }

  bool at_base = true;
  std::vector<double> psi_scratch;
  OutputValues direct(n_uptake);
  for (std::size_t i = 0; i < r.n_input; ++i) {
    const int p = r.input[i];
    // Does this input move the point? Where a bound is what defines it, the
    // bound's own row says so; a shut collar sits at a critical potential with no
    // condition to read, so nothing there can be shown to leave the point where
    // it is.
    double dpoint = 0.0;
    if (pinned) {
      dpoint = bound_moves_by_re_solving(p, single)
                   ? differenced_bound(l, theta, d, single, p, bound, s, at_base,
                                       psi_scratch)
                   : bound_dpoint(condition, p, n_layers);
    }
    const bool follows = shut || (pinned && dpoint != 0.0);
    if (follows) {
      // ⚠️ THE ROWS THIS WRITES ARE TOTALS, AND `dresidual` IS ZERO BECAUSE OF
      // IT. Both arms re-solve and land on the moved point, so the movement is
      // already inside every row; a live entry here would have the consumer's
      // assembly count it a second time. The held evaluation is not available to
      // take instead -- p* sits one step-in fraction from the bound, so a step
      // that moves the bound carries p* out of the perturbed feasible interval.
      if (!solved_row(l, theta, d, single, p, s, &b.branch, at_base, psi_scratch,
                      direct)) {
        // This input only. The row would be a difference across a change of
        // branch, which is a difference of two functions.
        for (std::size_t j = 0; j < r.n_output; ++j) {
          out.held[j * r.n_input + i] = util::na_value;
        }
        out.message += out.message.empty() ? "" : "; ";
        out.message += "no row for `" + par_name(p, n_layers) +
                       "`: no step within two decades keeps both arms on the " +
                       Leaf::operating_point_kind_name(b.branch.kind) + " branch";
        continue;
      }
      if (pinned) {
        out.dresidual[i] = 0.0;
      }
    } else {
      double dR = 0.0;
      held_row(l, theta, d, single, p, b.psi_star, s, "leaf_rows()", at_base,
                 psi_scratch, direct, dR);
      if (interior) {
        out.dresidual[i] = dR;
      } else {
        // `bound_row` has already divided by the bound's slope; what comes back
        // from here is the condition's own gradient, so that undoing it is the
        // consumer's one division rather than a second convention.
        out.dresidual[i] = -dpoint * out.residual_slope;
      }
    }
    for (std::size_t j = 0; j < r.n_output; ++j) {
      const std::size_t at = j * r.n_input + i;
      // The point's own held row is zero by construction: `held_row` has
      // just asserted both perturbed evaluations sat at exactly p*. Where the
      // arms followed the point instead, its row is the movement they measured.
      out.held[at] = (int(j) == point && !follows) ? 0.0 : direct[r.output[j]];
      out.zero[at] = zero_kind(out.held[at], p, interior);
    }
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);

  if (std::isfinite(out.residual_slope)) {
    double worst = 0.0;
    for (const double g : out.dresidual) {
      worst = std::max(worst, std::abs(g / out.residual_slope));
    }
    out.amplification = worst;
  }
  return out;
}

// --- the transpose ------------------------------------------------------------
//
// WHAT IT IS FOR. Everything above runs FORWARD: one perturbation per parameter,
// so the cost is however many parameters you ask about, which is why
// `leaf_gradient()` warns that asking for all fourteen is the most expensive
// thing you can do. That is the right direction for a leaf, where the parameter
// count is small. It is the wrong direction for a STAND. A forest model wants
// the sensitivity of a few whole-stand summaries -- total leaf area, basal area,
// biomass -- to EVERY trait at once, over a trajectory carrying thousands of
// individuals. Done forwards that is one full model run per trait; done
// backwards it is one run and one sweep whatever the trait count is.
//
// A backward sweep cannot use a forward derivative. When it reaches an
// individual it is holding the sensitivity of the final answer to that
// individual's OUTPUTS, and what it needs is the sensitivity to that
// individual's INPUTS. That is the transpose of what everything above computes,
// and this is it: one output adjoint in, one input adjoint out.
//
// THE CONTRACTION, which is all this is. The forward composite is
//
//   dY/dtheta  =  dY/dtheta|_psi  +  (dY/dpsi) * (-M/H)
//
// so for an output adjoint v the transpose is two scalars and one scaled row:
//
//   s    = v . (dY/dpsi)                a scalar
//   m    = -s / H                       a scalar
//   row  = v . (dY/dtheta|_psi)  +  m * M
//
// THE MATRIX (dY/dpsi)(dpsi*/dtheta) IS NEVER FORMED. It is rank one, because
// the collar potential is a single number, and in this direction it collapses to
// s and m. That is the whole economy and it is why the transpose is cheap: the
// combination drops from npars * n_outputs multiply-adds to npars + n_outputs,
// and nothing npars * n_outputs is ever allocated.
//
// ⚠️ WHAT IS *NOT* CHEAPER, AND SAYING SO HERE SAVES A DISAPPOINTMENT. There is
// no tape and no reverse mode through the model. The per-parameter quantities M
// and dY/dtheta|_psi are the same two perturbed evaluations `gradient_ift` takes,
// so AT THE LEAF the transpose costs what the forward gradient over the same
// parameters costs. Measured interleaved at one interior point: 1.00x over
// fourteen parameters, 1.00x over the nine that rebuild no spline, and a marginal
// cost per parameter within 0.4% of the forward path's. The saving is at the
// STAND, where this is the primitive one adjoint sweep calls per individual and
// the forward alternative is a whole trajectory per trait.
//
// THE IDENTITY THAT CHECKS IT, and it is exact:
//
//   <v, J u>  ==  <J^T v, u>
//
// for any output adjoint v and any input direction u. Both sides multiply the
// same numbers and sum them in a different order, so they agree to
// reassociation and to nothing looser --
// `test_gradient_transpose_matches_the_forward_jacobian` in tests/cpp requires
// that over the whole golden grid and over randomised v and u. Measured worst
// 1.41e-14, five orders below the solve's ~1e-09 floor, at a tolerance of 1e-12.

// The index of the OBJECTIVE among the reported outputs, or -1 where it is not
// reported at all. Profit is reported here, so it is that index.
//
// ⚠️ THE OBJECTIVE'S PSI-CHANNEL MUST BE EXCLUDED FROM `s`. At an interior
// optimum dprofit/dpsi = 0 by the envelope theorem, so profit's psi-channel
// contributes nothing and adding a measured central difference of it would
// double-count -- and that difference is not small noise either, being a
// difference of a FLAT maximum. #8 appends `profit` as a fifth output and gives
// it exactly that special case in the forward composite. The dot-product
// identity is what makes that one line rather than a hazard: it compares the
// transpose against the forward path column by column, so a forward path that
// special-cases profit and a transpose that does not fails immediately -- which
// is what it did when the two were first merged with this constant still unset.
inline constexpr int out_objective = out_profit;

// dY_j/dpsi as the TRANSPOSE must weight it, which is not always the measured
// central difference. Two outputs are exceptions, for two different reasons, and
// both are exceptions in `gradient_ift` too -- a transpose that disagreed with
// the forward path about either would be the transpose of a different operator.
//
//   * `collar` IS psi, so its psi-channel is exactly 1 and its direct term
//     exactly 0; `gradient_ift` writes dpsi*/dtheta into that row rather than
//     composing it. Using the measured difference instead is wrong by about
//     eps * psi / h_psi, i.e. ~1e-10 relative -- small, and still two orders
//     above the residual this transpose is required to leave.
//   * the objective's psi-channel is zero by the envelope theorem, above.
inline double psi_channel(int j, const OutputValues& dY_dpsi) {
  if (j == out_collar) {
    return 1.0;
  }
  if (j == out_objective) {
    return 0.0;
  }
  return dY_dpsi[j];
}

struct TransposeResult {
  // The outputs the transpose is taken at, and the operating-point
  // diagnostics, all of them `at()`'s own -- see `transpose_at`.
  double value[n_outputs];
  // One entry per requested parameter: v . dY/d(pars[k]).
  std::vector<double> adjoint;
  Status status = Status::Error;
  bool used_ift = false;
  double H = util::na_value;
  double stationarity = util::na_value;
  // `s` above: the sensitivity of the weighted output to the collar potential,
  // and the one number the rank-one term collapses to. NA off the composite.
  double psi_adjoint = util::na_value;
  std::string message;

  // Cleared BEFORE anything can throw, for `Result::reset`'s reason: this is a
  // caller's struct and may be a reused one, so a call that stops partway must
  // not leave the previous point's numbers sitting in it looking current.
  void reset(std::size_t npars) {
    for (int j = 0; j < n_outputs; ++j) {
      value[j] = util::na_value;
    }
    adjoint.assign(npars, util::na_value);
    status = Status::Error;
    used_ift = false;
    H = util::na_value;
    stationarity = util::na_value;
    psi_adjoint = util::na_value;
    message.clear();
  }
};

// The transposed gradient at one operating point.
//
// ⚠️ THE OPERATING-POINT CLASSIFICATION IS `at()`'s, NOT A SECOND COPY OF IT, and
// that is why this starts by calling `at()` with no parameters. The composite is
// valid only where stationarity holds; at a pinned optimum the collar follows a
// bound and -M/H is not its derivative. Rather than repeat that test -- and with
// it the curvature, the sentinel handling and the narrow-bracket detector -- this
// runs the forward entry point with `npars == 0`, which does the solve, the
// classification and dY/dpsi and then loops over nothing. So `status`,
// `used_ift`, `H` and `stationarity` are the same numbers the forward path would
// report at this point BY CONSTRUCTION, this route refuses wherever that one
// refuses (`at()` throws and the throw propagates), and a change to the
// classification cannot reach one route without reaching the other.
//
// The fallback is transposed too. Where `at()` declines the composite it
// differences the whole solve, and the transpose of that is the same contraction
// against a plainer Jacobian -- no s, no m, no H. So the identity below holds at
// pinned and shut-down points as well, which is most of what makes it worth
// having.
inline void transpose_at(Leaf& l, const double* theta, const Drivers& d,
                         bool single, const int* pars, std::size_t npars,
                         const double* v, const Settings& s,
                         TransposeResult& out) {
  out.reset(npars);
  check_pars(pars, npars, n_soil_layers(d, single), "leaf_gradient_transpose()");

  Result point;
  at(l, theta, d, single, nullptr, 0, s, point);

  for (int j = 0; j < n_outputs; ++j) {
    out.value[j] = point.value[j];
  }
  out.status = point.status;
  out.used_ift = point.used_ift;
  out.H = point.H;
  out.stationarity = point.stationarity;

  // psi* as `at()` recorded it. `outputs()` copies `opt_root_psi_` into the
  // collar slot, so this is that member and not a re-derivation of it -- the
  // leaf itself is back at base parameters and unsolved by now.
  const double psi_star = out.value[out_collar];

  std::vector<double> psi_scratch;
  bool at_base = true;

  if (out.used_ift) {
    // s, the first of the two scalars. `rounded()` here and below is not for
    // agreement with R -- there is no R implementation of this to agree with --
    // but so that no compiler contracts a multiply-add and the identity's
    // residual is the same number under gcc and clang. It costs a store per
    // term against a model evaluation per parameter.
    double s_psi = 0.0;
    for (int j = 0; j < n_outputs; ++j) {
      s_psi += rounded(v[j] * psi_channel(j, point.dY_dpsi));
    }
    out.psi_adjoint = s_psi;
    const double m = -s_psi / point.H;

    OutputValues direct;
    for (std::size_t k = 0; k < npars; ++k) {
      double M = 0.0;
      held_row(l, theta, d, single, pars[k], psi_star, s,
                 "leaf_gradient_transpose()", at_base, psi_scratch, direct, M);
      // v . dY/dtheta|_psi. `collar` is skipped rather than summed: its direct
      // term is zero by construction, `outputs_at` having just asserted that
      // both sides sit at exactly psi*.
      double row = 0.0;
      for (int j = 0; j < n_outputs; ++j) {
        if (j == out_collar) {
          continue;
        }
        row += rounded(v[j] * direct[j]);
      }
      out.adjoint[k] = row + rounded(m * M);
    }
  } else {
    // The fallback, transposed: difference the solve and contract. No output is
    // exceptional here, `collar` included, because nothing is being composed --
    // which is exactly why `gradient_fd` has no special case either.
    OutputValues solved;
    for (std::size_t k = 0; k < npars; ++k) {
      solved_row(l, theta, d, single, pars[k], s, nullptr, at_base, psi_scratch,
                 solved);
      double row = 0.0;
      for (int j = 0; j < n_outputs; ++j) {
        row += rounded(v[j] * solved[j]);
      }
      out.adjoint[k] = row;
    }
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);
}

// --- the environment rows, contracted -----------------------------------------
//
// `transpose_at` above is a transpose in its STRUCTURE -- one `s`, one `m`, the
// rank-one channel collapsed once rather than per parameter -- but every column
// is still two perturbed evaluations. That is the right trade for the fourteen
// traits, which have no closed form. It is the wrong one for the environment,
// which now has.
//
// So this is the other half: given the adjoint carried on profit, the rows for
// the light and for each soil layer, analytically. A stand sweep calls it once
// per cohort per stage, where differencing would be 2(L+1) solves.
//
// THE TWO CHANNELS ARE NOT ONE LOOP, and the reason is structural rather than
// tidiness. A soil row is a PRICE times a supply derivative: the layer moves
// water, and `marginal_price_water` converts water into carbon. The light row
// has no supply derivative to multiply -- at a fixed collar PPFD moves no water
// at all -- and is a direct kernel partial closed by the ci root-find. Writing
// them as one loop over "environment parameters" would force one of the two into
// the other's shape.
//
// Profit only. The other four outputs' environment rows are still `transpose_at`'s
// to difference; this contracts the one output a census functional reads.
struct ProfitEnvDerivatives {
  // Named for what they are rather than for their place in a matrix: a reader
  // should not have to know which way round "row" and "column" run here.
  double dprofit_dlight = util::na_value;
  // One per layer, in the caller's layer order.
  std::vector<double> dprofit_dpsi_soil;
  // False where the operating point has no analytic rows: a branch kink leaves
  // the supply derivative undefined, and a pinned or shut-down point is not an
  // interior optimum, so the envelope step these rows take does not hold. A
  // caller that ignores this reads sentinels as numbers.
  bool usable = false;
  std::string message;
  // True where the rows above carry the bound's movement as well as the direct
  // partial. At an interior optimum the argmax contributes nothing to profit by
  // the envelope theorem; at a pin the operating point IS the bound, so it moves
  // with the inputs and the profit row picks the term back up. Reported because
  // the two are different derivations, and a caller reading a pinned row as an
  // envelope row would be reading the wrong theory's number.
  bool pinned = false;
  // Which bound, where `pinned`. Reported rather than re-derived so a consumer
  // that follows the bound cannot pick a different arm from the one these rows
  // were built on; meaningless where `pinned` is false.
  Leaf::WhichBound bound = Leaf::WhichBound::Wet;

  void reset(std::size_t n_layers) {
    dprofit_dlight = util::na_value;
    dprofit_dpsi_soil.assign(n_layers, util::na_value);
    usable = false;
    pinned = false;
    message.clear();
  }
};

// The three pieces a soil row is built from, at whatever collar the leaf is
// currently seated at: the radiation row, the price that turns water into carbon,
// and each layer's supply derivative. False where a piece of it is undefined.
//
// The price is handed back rather than already multiplied in because a pinned
// point needs a DIFFERENT price from an interior one, and what differs is a term
// inside the price rather than a factor on the row.
inline bool profit_env_rows_here(Leaf& l, double& light, double& price,
                                 std::vector<double>& duptake_dpsi_soil,
                                 std::string& message) {
  const std::vector<double>& psi_soil = l.supply_psi_soil();
  price = l.marginal_price_water();
  light = l.dprofit_dPPFD();
  if (!std::isfinite(price) || !std::isfinite(light)) {
    message = "the supply derivative or the ci residual is undefined here";
    return false;
  }
  l.dE_from_soil_dpsi_soil(l.opt_root_psi_, psi_soil, duptake_dpsi_soil);
  for (const double v : duptake_dpsi_soil) {
    if (!std::isfinite(v)) {
      message = "a layer sits on a branch kink, so its row does not exist";
      return false;
    }
  }
  return true;
}

// Requires a solved operating point on `l` -- the caller's own solve, not one
// taken here, so this reads state rather than moving it.
inline void profit_env_derivatives(Leaf& l, ProfitEnvDerivatives& out) {
  const std::vector<double>& psi_soil = l.supply_psi_soil();
  out.reset(psi_soil.size());

  using Kind = Leaf::OperatingPointKind;
  const Kind kind = l.operating_point_kind();
  Leaf::WhichBound bound = Leaf::WhichBound::Wet;
  const bool pinned = pinned_bound(kind, bound);

  if (kind != Kind::Interior && !pinned) {
    out.message =
        "the environment rows need an interior optimum or a bound to follow; "
        "this point is " + std::string(Leaf::operating_point_kind_name(kind));
    return;
  }

  double light = 0.0;
  double price = 0.0;
  std::vector<double> duptake;
  if (!profit_env_rows_here(l, light, price, duptake, out.message)) {
    return;
  }
  std::vector<double> soil(duptake.size(), 0.0);

  if (pinned) {
    // The envelope theorem fails at a boundary: the operating point is not
    // stationary there, it IS the bound, so it moves with the inputs and the
    // profit row carries that movement.
    //
    //   dProfit*/du = dProfit/du|_p  +  (dProfit/dp) * dB/du
    //
    // ⚠️ AND THE FIRST TERM IS NOT THE INTERIOR ROW. `marginal_price_water` is
    // λ·kmax·f(p)/S, which is what dProfit/dE_up REDUCES TO once dProfit/dp is
    // zero -- the stationarity condition is inside it. Where dProfit/dp is `nu`
    // instead, the held-collar price is that value plus nu/S, S being the
    // soil-to-collar conductance the price is already built on. Dropping the
    // correction leaves both terms finite and the row wrong by 15% at a dry pin;
    // at a WET pin it cancels the second term exactly and reads 8% out.
    const Leaf::BoundRow b = l.bound_row(bound);
    if (!b.finite) {
      out.message = "the bound this point is pinned to has no derivative here";
      return;
    }
    const double nu = l.dprofit_droot_collar_psi(l.opt_root_psi_);
    if (!std::isfinite(nu)) {
      out.message = "the marginal profit at the bound is not finite";
      return;
    }
    if (b.d_dpsi_soil.size() != soil.size()) {
      out.message = "the bound's row and the soil rows disagree on the layer count";
      return;
    }
    // dmarginal_profit_duptake_slope builds dProfit/dE_up from the cost and
    // assimilation kernels directly, with no stationarity anywhere in it, so it
    // is the held-collar price wherever the point sits. marginal_price_water
    // agrees with it at an interior optimum and only there.
    //
    // The interior branch keeps the price it had: the two are the same number by
    // the first-order condition, they differ in the last bits, and a gradient
    // that already answers must not move.
    price = l.dmarginal_profit_duptake_slope();
    if (!std::isfinite(price)) {
      out.message = "the held-collar price of water is undefined at this pin";
      return;
    }
    for (std::size_t j = 0; j < soil.size(); ++j) {
      soil[j] = duptake[j] * price + nu * b.d_dpsi_soil[j];
    }
    // The light row gains nothing: neither bound reads radiation, so dB/dlight
    // is exactly zero. It needs no price correction either -- it is closed by the
    // ci residual at a held collar rather than by the stationarity condition.
    out.pinned = true;
    out.bound = bound;
  } else {
    for (std::size_t j = 0; j < soil.size(); ++j) {
      soil[j] = duptake[j] * price;
    }
  }

  out.dprofit_dlight = light;
  out.dprofit_dpsi_soil = soil;
  out.usable = true;
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
// derivative of dprofit/du, which `profit_env_derivatives` already computes analytically,
// and so to difference THAT in the collar: two seatings, no re-driving, every
// column at once.
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
// So d2profit/dp du has to come from the marginal profit itself, and that is
// what report 05 s7.3's rank-two factorisation over the state directions is for:
// R = dprofit/dp reads the state only through total uptake and uptake's own
// collar sensitivity, so dR/du = a * dE_up/du + b * d(dE_up/dp)/du with two
// scalars shared across every direction. Both vectors are closed form --
// `MultiLayerRoots::duptake_dpsi_soil` and `d2uptake_dpsi_dpsi_soil` -- and the
// pair is recovered from two perturbed evaluations in directions of different
// families, which is verified out of sample in
// plant-dev/docs/design/verify-factorisation.R.

}  // namespace gradient
}  // namespace phylloptim

#endif
