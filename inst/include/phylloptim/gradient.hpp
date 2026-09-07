// -*-c++-*-
#ifndef PHYLLOPTIM_GRADIENT_HPP_
#define PHYLLOPTIM_GRADIENT_HPP_

// Trait gradients, composed here rather than in R (#4).
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
// has; the caller's parameterisation Jacobian (one study maps 40 fitted
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

// Re-exported here under the names this namespace's callers already use, so the
// move is not a rename. One definition, in the model; two ways to reach it.
using phylloptim::n_traits;
using phylloptim::par_vcmax_25;   using phylloptim::par_stem_c;
using phylloptim::par_stem_P50;   using phylloptim::par_root_c;
using phylloptim::par_root_P50;   using phylloptim::par_TF24_beta2;
using phylloptim::par_jmax_25;    using phylloptim::par_a;
using phylloptim::par_curv_fact_elec_trans;
using phylloptim::par_curv_fact_colim;
using phylloptim::par_TF24_cost_scale; using phylloptim::par_R_d_25;
using phylloptim::par_JS22_gamma; using phylloptim::par_CMax_a;
using phylloptim::par_CMax_b;     using phylloptim::par_kmax;
using phylloptim::par_resistance; using phylloptim::par_CF77_lambda;
using phylloptim::par_TF24_floor_lambda_o;
using phylloptim::par_PPFD;

// The parameter enumeration now lives in leaf_model.hpp, beside the set_traits
// whose argument order it IS. It was here while only this file indexed it; the
// differentiable surface indexes it too, and a list two headers can disagree
// about is the hazard the comment above it warns of.

// How long `theta` is: one entry per `par_names()` below, indexed by the same
// enumeration.
//
// ⚠️ NOT `n_pars`, WHICH IS LONGER. The pack the leaf's kernels read also holds
// slots a calibration never proposes -- the light a cohort stands in is seated
// from the drivers -- so reading `theta` out to `n_pars` runs past its end. The
// proposed slots come FIRST in the enumeration and `par_PPFD` is where they stop,
// which is why this is spelled as that index rather than as a count to keep in
// step by hand.
inline constexpr int n_theta = par_PPFD;

// ⚠️ THESE ARE THE SLOTS A CALIBRATION PROPOSES, which is fewer than the pack
// holds: `par_PPFD` is seated per cohort or per observation from the drivers and
// is never fitted, so it has a slot in the enumeration and no name here. R's
// `.gradient_par_names()` is compared against this list by
// `test-gradient-batch.R`, so the two cannot drift apart without a failure.
inline const std::vector<std::string>& par_names() {
  static const std::vector<std::string> names{
      "vcmax_25",  "stem_c",              "stem_P50",
      "root_c",    "root_P50",            "TF24_beta2",
      "jmax_25",   "a",                   "curv_fact_elec_trans",
      "curv_fact_colim", "TF24_cost_scale", "R_d_25",
      "JS22_gamma", "CMax_a", "CMax_b",
      "leaf_specific_conductance_max",
      "resistance",
      "CF77_lambda_",
      "TF24_floor_lambda_o"};
  return names;
}

// --- the five differentiated outputs ------------------------------------------
//
// A, gc, psi_stem, collar and profit, in that order, which is R's
// `.gradient_output_names`. `collar` is psi* itself, which is what makes
// dcollar/dtheta equal dpsi*/dtheta and lets the two routes below compute the
// same quantity by different means.
//
// ⚠️ APPENDING IS SAFE AND REORDERING IS NOT, exactly as for `par_names` above.
//
// The first four are what a gas-exchange calibration OBSERVES. `profit` is here
// because it is what plant CONSUMES: `leaf.profit_`, not `assim_colimited_`, is
// the carbon that reaches its mass budget, so until #87 the two sets were
// disjoint and no gradient this package produced reached a demographic model.
inline constexpr int n_outputs = 5;
inline constexpr int out_psi_stem = 2;
inline constexpr int out_collar = 3;
inline constexpr int out_profit = 4;

inline const std::vector<std::string>& output_names() {
  static const std::vector<std::string> names{"A", "gc", "psi_stem", "collar",
                                              "profit"};
  return names;
}

// Read straight off the members rather than through `operating_point_values()`,
// which is what R has to use. Bit-identical: that reader copies these same five
// fields into positions 3, 5, 0, 1 and 6 of its twelve, and the three columns it
// computes rather than copies (uptake, lambda, g1_eff) are not among them.
//
// ⚠️ That is why `profit` was cheap to add here and `uptake` would not be. Every
// output in this list has to be a field R COPIES; `uptake` is one R sums over
// the finite soil layers, so adding it means reproducing that summation -- and
// its order -- on this side too.
// --- which model's optimum is being differentiated ---------------------------
//
// `kCollar` is the production path: the TF24 cost maximised over the root-collar
// potential, which is what this file did unconditionally. Anything >= 0 is a
// `Leaf::CostCurve` value and selects a STEM route -- psi_stem with the upstream
// potential pinned at psi_soil.
//
// ⚠️ THE COLLAR BRANCH IS TEXTUALLY WHAT IT ALWAYS WAS. These wrappers exist so
// the six call sites below stop naming a solver, not to change one: the collar
// arms call exactly the same functions in the same order, which is what keeps the
// R-versus-C++ bit-for-bit contract and the golden files intact.
inline constexpr int kCollar = -1;

// Seat the model on the leaf ONCE, at the top of a gradient, rather than naming a
// curve at each of the calls below. `Settings::curve` is still an int here --
// nothing outside this file sees it -- but the Leaf's own API takes names only, so
// this is the single place the two representations meet.
//
// ⚠️ SEAT IT EVEN ON THE PRESCRIBED PATH, which never calls `route_solve`. The
// stem entry points read the seated curve, so without this a prescribed gradient
// would silently differentiate whatever model the leaf was last used for.
inline void route_seat(Leaf& l, int curve) {
  l.set_model(curve < 0 ? Leaf::CostCurve::TF24
                        : static_cast<Leaf::CostCurve>(curve),
              curve < 0);
}

inline void route_solve(Leaf& l, int curve) {
  if (curve < 0) {
    l.find_root_collar_psi();
  } else {
    l.optimise();
  }
}

inline double route_psi(const Leaf& l, int curve) {
  return curve < 0 ? l.opt_root_psi_ : l.opt_psi_stem_;
}

inline void route_evaluate(Leaf& l, int curve, double psi) {
  if (curve < 0) {
    l.evaluate_root_collar_psi(psi);
  } else {
    l.evaluate_psi_stem_at(psi);
  }
}

// ⚠️ THE FLAG IS THE POINT; the unchecked form is this with a null. Both routes
// return a hard 0.0 SENTINEL on their shut-down and reversed-gradient exits, and a
// bare zero is indistinguishable from a stationary point -- see
// Leaf::dprofit_droot_collar_psi's header. There used to be two functions here,
// the shorter one being a call to the longer with `nullptr`; there is one now.
inline double route_dprofit_checked(Leaf& l, int curve, double psi,
                                    bool* feasible = nullptr) {
  if (curve < 0) {
    return l.dprofit_droot_collar_psi(psi, feasible);
  }
  const std::vector<double> r = l.dprofit_dpsi_stem_checked(psi);
  if (feasible != nullptr) {
    *feasible = r[1] != 0.0;
  }
  return r[0];
}

inline double route_dprofit(Leaf& l, int curve, double psi) {
  return route_dprofit_checked(l, curve, psi, nullptr);
}

// ⚠️ WHICH OUTPUT IS THE DECISION VARIABLE MOVES WITH THE ROUTE. On the collar
// route `collar` IS psi*, so `dcollar/dtheta` is `dpsi* / dtheta` and assigned
// rather than differenced. On a stem route `psi_stem` plays that part and the
// collar is never solved for.
inline int route_decision(int curve) {
  return curve < 0 ? out_collar : out_psi_stem;
}

inline void outputs(const Leaf& l, double* y) {
  y[0] = l.assim_colimited_;
  y[1] = l.stom_cond_CO2_;
  y[2] = l.opt_psi_stem_;
  y[3] = l.opt_root_psi_;
  y[4] = l.profit_;
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
inline bool outputs_at(Leaf& l, double psi, double* y, int curve = kCollar) {
  route_evaluate(l, curve, psi);
  if (!util::identical(route_psi(l, curve), psi)) {
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

enum class Method { Auto, Ift, Fd };

// What the OPERATING POINT is, not whether the call worked -- except for `Error`,
// which is the batch's per-row failure and is the reason this is a status rather
// than an exception (see `batch` below).
// `Prescribed` and `Clamped` are the two the caller-supplied-psi path can end in
// (#88): the collar was imposed, and in the second case the feasible interval
// moved it, so the derivative belongs to the BOUND and is not reported.
// ⚠️ APPEND AFTER `Error`, DO NOT INSERT BEFORE IT. The names are
// source-compatible either way, but the VALUES are not: a consumer holding the
// integer would break quietly if `Error` moved from 3 to 5. `Prescribed` and
// `Clamped` are therefore last, out of narrative order, on purpose.
enum class Status { Interior, Pinned, NoGradient, Error, Prescribed, Clamped };

// ⚠️ NO `default:`. An exhaustive switch turns the next member added here into a
// compiler diagnostic; a `default: return "error"` would silently LABEL it
// "error" instead, which is exactly the kind of quiet mislabelling `Prescribed`
// and `Clamped` were added to avoid.
inline std::string status_name(Status s) {
  switch (s) {
    case Status::Interior:   return "interior";
    case Status::Pinned:     return "pinned";
    case Status::NoGradient: return "no-gradient";
    case Status::Prescribed: return "prescribed";
    case Status::Clamped:    return "clamped";
    case Status::Error:      return "error";
  }
  return "error";  // unreachable; silences -Wreturn-type on a bad cast
}

struct Settings {
  double step = 1e-6;
  double stationarity_tol = 1e-8;
  Method method = Method::Auto;
  bool fast_stem_curve = true;
  // Which model's optimum to differentiate. `kCollar` is the production path and
  // the default, so an existing caller is unaffected.
  int curve = kCollar;
  // ⚠️ THE STEP THE FINITE-DIFFERENCE ARM USES, separately from `step`. That arm
  // differences the SOLVE, so it has to clear the solver's own resolution: the
  // collar route root-finds `dprofit == 0` to ~1e-15, while the stem routes scan
  // and refine to ~1e-06 of the bracket. Differencing a stem solve at 1e-06
  // returns quantisation -- measured, the wrong SIGN for JW26. `step` still drives
  // the composite, which perturbs at fixed psi and never re-solves.
  double fd_step = 1e-6;
  // |A|max held fixed across perturbations, for ProfitMax's `Scaled` link. Zero
  // means "not pinned"; the caller seeds it from the base solve.
  // ProfitMax only: re-solve `|A|max` after each perturbation, because `apply()`
  // clears it. Named for what it does rather than for a value to hold, since
  // holding one is exactly what made this a partial derivative.
  bool resolve_profitmax = false;
};

// A collar potential the caller imposes, in place of the one `at` would solve
// for, plus how that collar responds to each parameter (#88).
//
// ⚠️ `dpsi_dtheta` is npars long and in `pars` order, NOT n_theta: it is one value
// per parameter ASKED FOR, because that is the vector a caller integrating its
// own sensitivity state is carrying. Null means zero -- the partial at fixed
// collar -- which is a different statement from the solving path's "derive it".
struct Prescribed {
  double psi = 0.0;
  const double* dpsi_dtheta = nullptr;
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
  // The mixed partials d2profit/dpsi dtheta, npars long, and dY/dpsi at fixed
  // traits. Kept rather than consumed: `-M/H` is what the composite needs, but M
  // and H are also the coefficients of a caller's own sensitivity ODE for psi
  // (traitecoevo/plant#614), and they cannot be recovered once divided.
  std::vector<double> M;
  double dY_dpsi[n_outputs];
  // The collar the outputs were evaluated at -- psi* when solved, and when
  // prescribed the value actually USED, which differs from the requested one
  // exactly when `status` is Clamped.
  double psi = util::na_value;
  std::string message;

  void reset(std::size_t npars) {
    for (int j = 0; j < n_outputs; ++j) {
      value[j] = util::na_value;
      dY_dpsi[j] = util::na_value;
    }
    grad.assign(npars * n_outputs, util::na_value);
    M.assign(npars, util::na_value);
    status = Status::Error;
    used_ift = false;
    H = util::na_value;
    stationarity = util::na_value;
    psi = util::na_value;
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
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve) {
  // THE FAST PATH FOR stem_b. The stem
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
  // different model and disagrees by 3e-4.
  if (fast_stem_curve && only == par_stem_P50) {
    l.perturb_stem_P50(theta[par_stem_P50]);
    return;
  }
  // AND THE WAY BACK OUT OF IT (#74), which is what makes the shortcut worth its
  // 24.5x through a batch instead of 2.4x. `set_traits()` below rebuilds the stem
  // curve whenever `stem_b != stem_b_spline_`, and that third clause is what
  // returns a shortcut-displaced leaf to a rebuilt one -- so a restore that
  // followed a `perturb_stem_P50()` always paid for a rebuild, once per observation,
  // whatever `pars` contained. Undoing the displacement WITH the shortcut leaves
  // the clause false and the splines alone.
  //
  // Bit-identical by construction rather than by measurement. Everything
  // `perturb_stem_P50()` writes -- `stem_P50`, and `stem_b`/`psi_crit` derived
  // from it -- is a pure function of `(stem_P50, stem_c)` computed by the same
  // expression `set_traits()` uses, so restoring the base P50 restores all three
  // to the base bit pattern. The splines here ARE the ones built at
  // `stem_b_spline_`, not a rescaled copy, and at `stem_b == stem_b_spline_` all
  // four `stem_curve_*` accessors take their scale == 1 branch and read them
  // directly.
  //
  // ⚠️ The two guards are both load-bearing, and neither is an optimisation.
  // Displacement can only be created by `perturb_stem_P50()`, which is sound only
  // when everything else is already at base -- so a displaced leaf is one whose
  // stem_c/root_* are the values `stem_b_spline_` was validated against. Without
  // the equality test that argument is gone: in a batch with a theta MATRIX the
  // next row's restore moves stem_b somewhere new, and pushing it through the
  // shortcut would rescale off a spline built for a different curve.
  //
  // The equality is on the DERIVED scale, not on the trait, because
  // `stem_b_spline_` records a `b`. Same expression as the one inside
  // `perturb_stem_P50()`, so a genuine round trip compares exactly equal.
  if (fast_stem_curve && l.stem_b != l.stem_b_spline_ &&
      Leaf::weibull_b_from_P50(theta[par_stem_P50], l.stem_c) ==
          l.stem_b_spline_) {
    l.perturb_stem_P50(theta[par_stem_P50]);
  }
  l.set_traits(theta[par_vcmax_25], theta[par_stem_c], theta[par_stem_P50],
               theta[par_root_c], theta[par_root_P50],
               theta[par_TF24_beta2], theta[par_jmax_25],
               theta[par_a], theta[par_curv_fact_elec_trans],
               theta[par_curv_fact_colim], theta[par_TF24_cost_scale],
               theta[par_R_d_25], theta[par_JS22_gamma],
               theta[par_CMax_a], theta[par_CMax_b]);
  // Not traits, so they are set directly rather than through set_traits -- and
  // nothing is derived from either, so a bare write leaves no stale state behind.
  l.CF77_lambda_ = theta[par_CF77_lambda];
  l.TF24_floor_lambda_o = theta[par_TF24_floor_lambda_o];
  if (single) {
    // R's `series_resistance()`: a default-constructed network carrying one
    // series resistance in `r_R_V_sum`, which is that field's own meaning with
    // one layer and no vulnerability-weighted term. Built here rather than
    // copied from `d.root_network` so that the two implementations cannot
    // disagree about what the other four fields hold.
    RootNetwork net;
    net.r_R_V_sum.assign(1, theta[par_resistance]);
    l.set_physiology(net, d.PPFD, d.psi_soil, d.soil_depth, theta[par_kmax],
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  } else {
    l.set_physiology(d.root_network, d.PPFD, d.psi_soil, d.soil_depth,
                     theta[par_kmax], d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  }
}

// --- the two routes ----------------------------------------------------------

// ⚠️ THE INVARIANT BOTH LOOPS BELOW MAINTAIN, AND IT IS THE FIX FOR #72: every
// parameter's gradient is taken from the BASE point.
//
// The loops restore base once at the END, not between parameters, because every
// parameter's setter normally goes through the full `set_traits()` +
// `set_physiology()` path and restores everything on the way. The `stem_b`
// shortcut is the one that does not: `perturb_stem_P50()` rescales the stem spline
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
  return s.fast_stem_curve && par == par_stem_P50;
}

// The implicit-function composite. Two perturbed evaluations per parameter,
// neither of which re-solves the model: `dprofit` at the UNPERTURBED psi* gives
// the mixed partial, and the outputs at that same psi* give the direct term.
//
// ONE composite for both paths. `dpsi_dtheta` null means "derive it by the
// implicit function theorem", which is right where the collar was solved for;
// non-null means the caller imposed the collar and knows how it moves. Keeping
// them in one function is what makes the two paths' agreement at psi* a real
// assertion rather than two implementations that happen to line up.
inline void gradient_ift(Leaf& l, const double* theta, const Drivers& d,
                         bool single, const int* pars, std::size_t npars,
                         double psi_star, double H, const double* dY_dpsi,
                         const Settings& s, const double* dpsi_dtheta,
                         bool envelope, double* M_out, double* out) {
  double th[n_theta];
  double up[1 + n_outputs];
  double dn[1 + n_outputs];
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    const int p = pars[k];
    if (takes_shortcut(p, s) && !at_base) {
      apply(l, theta, d, single, -1, s.fast_stem_curve);
    }
    at_base = false;
    const double h = step_for(p, theta[p], s.step);
    for (int side = 0; side < 2; ++side) {
      // Up first, then down: R evaluates `up <- side(1)` before `dn <- side(-1)`
      // and both mutate the leaf.
      std::copy(theta, theta + n_theta, th);
      th[p] = side == 0 ? theta[p] + h : theta[p] - h;
      apply(l, th, d, single, p, s.fast_stem_curve);
      double* dst = side == 0 ? up : dn;
      // Evaluate first, then read dprofit at the same fixed collar -- R's order,
      // and `evaluate_root_collar_psi` is what seats the state `dprofit` reads.
      if (s.resolve_profitmax) {
        // `apply()` above cleared ProfitMax's normaliser, so it has to be put back
        // before anything reads the objective. RE-SOLVED, not pinned at the base
        // point: |A|max is found by a root-find now rather than by a scan, so it is
        // differentiable in the traits and letting it follow them makes this the
        // TOTAL derivative -- which is what a fit needs. Must match R/gradient.R's
        // `pinned$reseat`, which calls the same thing.
        l.prepare_profitmax();
      }
      if (!outputs_at(l, psi_star, dst + 1, s.curve)) {
        util::stop_infeasible(
            "gradient_active_set",
            "leaf_gradient(): perturbing `" + par_names()[std::size_t(p)] +
                   "` moved the feasible collar interval past psi*, so the "
                   "operating point could not be evaluated there. This point is "
                   "on an active-set boundary; lower `stationarity_tol` or "
                   "difference the solve directly.");
      }
      // ⚠️ AGAIN BEFORE THE DERIVATIVE READ, because R re-seats inside EVERY
      // accessor rather than once per perturbation, and `prepare_profitmax()`
      // mutates the leaf's rate fields while it maximises. Once here and twice
      // there leaves the two reading `dprofit` at different states: measured, the
      // batch and the one-observation route then agreed only to 3e-05 where they
      // are required to agree bit-for-bit.
      if (s.resolve_profitmax) {
        l.prepare_profitmax();
      }
      dst[0] = route_dprofit(l, s.curve, psi_star);
    }
    // M = d2profit/dpsi dtheta, with psi held FIXED at psi*.
    const double M = (up[0] - dn[0]) / (2.0 * h);
    M_out[k] = M;
    const double d_psi = dpsi_dtheta == nullptr ? -M / H : dpsi_dtheta[k];
    for (int j = 0; j < n_outputs; ++j) {
      const double direct = (up[1 + j] - dn[1 + j]) / (2.0 * h);
      out[k * n_outputs + j] = direct + rounded(dY_dpsi[j] * d_psi);
    }
    // `collar` is not an output of the evaluation -- it IS psi*, held fixed, so
    // its direct term is zero by construction and the composite reduces to
    // dpsi*/dtheta. Set explicitly rather than left as the difference of two
    // identical numbers.
    out[k * n_outputs + route_decision(s.curve)] = d_psi;
    // The envelope theorem, ASSIGNED for the same reason `collar` is: profit's
    // indirect term is identically zero at a stationary point, so stating that
    // beats multiplying a measured near-zero by dpsi/dtheta. It is also immune
    // to a non-finite `d_psi` -- which a CALLER supplies on the prescribed path
    // -- where `0 * x` would be NaN in this one column while the other four
    // carried +-Inf.
    if (envelope) {
      out[k * n_outputs + out_profit] =
          (up[1 + out_profit] - dn[1 + out_profit]) / (2.0 * h);
    }
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);
}

// The fallback: a central difference of the whole solve. Correct at a pinned
// optimum because it differences the CONSTRAINED answer, which is exactly what
// the composite cannot do.
inline void gradient_fd(Leaf& l, const double* theta, const Drivers& d,
                        bool single, const int* pars, std::size_t npars,
                        const Settings& s, double* out) {
  double th[n_theta];
  double up[n_outputs];
  double dn[n_outputs];
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    const int p = pars[k];
    if (takes_shortcut(p, s) && !at_base) {
      apply(l, theta, d, single, -1, s.fast_stem_curve);
    }
    at_base = false;
    const double h = step_for(p, theta[p], s.fd_step);
    for (int side = 0; side < 2; ++side) {
      std::copy(theta, theta + n_theta, th);
      th[p] = side == 0 ? theta[p] + h : theta[p] - h;
      apply(l, th, d, single, p, s.fast_stem_curve);
      route_solve(l, s.curve);
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
               Result& out, const Prescribed* prescribed = nullptr) {
  out.reset(npars);

  apply(l, theta, d, single, -1, s.fast_stem_curve);
  // Seat the model before anything reads it. `apply` runs `set_traits`, which
  // clears the solved state but NOT the seated curve, so this could be hoisted
  // out of the per-observation loop -- it is here because `at()` is a public
  // entry point in its own right and must not depend on `batch` having run.
  route_seat(l, s.curve);

  // Two ways in. The default SOLVES for the collar potential; `prescribed`
  // IMPOSES one, which is what a caller tracking the optimum rather than
  // finding it has (#88).
  //
  // ⚠️ `psi_star` keeps its name on both paths and is no longer always the
  // argmax. It is "the collar the outputs were evaluated at", which is what
  // every use of it below actually means.
  bool clamped = false;
  if (prescribed != nullptr) {
    route_evaluate(l, s.curve, prescribed->psi);
    // Exact equality, for the reason `outputs_at` documents above.
    clamped = !util::identical(route_psi(l, s.curve), prescribed->psi);
  } else {
    route_solve(l, s.curve);
  }
  const double psi_star = route_psi(l, s.curve);
  out.psi = psi_star;
  outputs(l, out.value);

  // Is the composite's premise true HERE? Stationarity is what the whole
  // derivation rests on and it fails at a pinned optimum, where psi* is a bound,
  // dprofit is not zero at the answer, and -M/H is not the bound's derivative.
  // The formula does not fail loudly: at a boundary-soil point the true gradient is
  // ~1e-08 and the bare composite returns O(1).
  //
  // So the premise is TESTED. The test is the implied Newton step
  // |dprofit(psi*) / H|, a distance in MPa that needs no scale of its own: over
  // this package's 288-point grid the worst interior point is 4.8e-11 and the
  // mildest pinned one 6.3e-06, so the 1e-08 default sits in an empty band four
  // orders wide on each side.
  const double h_psi = std::max(std::abs(psi_star), 1.0) * s.step;
  // ⚠️ WITH ITS FEASIBILITY, not bare. `dprofit_droot_collar_psi` returns a hard
  // 0.0 SENTINEL on its shut-down and reversed-gradient exits, and a bare zero is
  // indistinguishable from a stationary point. The solving path got away with the
  // value alone because `H` collapses to zero too and `usable` catches the pair;
  // the prescribed path does not divide by `H`, so it would adopt the sentinel as
  // if it were dprofit/dpsi.
  bool resid_feasible = false;
  const double resid =
      route_dprofit_checked(l, s.curve, psi_star, &resid_feasible);
  // Named halves: `f(a) - f(b)` has unspecified operand order in C++ and
  // left-to-right order in R, and `dprofit_droot_collar_psi` mutates the leaf.
  const double d_hi = route_dprofit(l, s.curve, psi_star + h_psi);
  const double d_lo = route_dprofit(l, s.curve, psi_star - h_psi);
  const double H = (d_hi - d_lo) / (2.0 * h_psi);
  // H == 0 with resid == 0 is the shut-down signature: dprofit returns a
  // sentinel zero there rather than a derivative, so the ratio would be 0/0.
  // H > 0 would not be a maximum. Both mean the composite has nothing to stand
  // on.
  const bool usable = std::isfinite(H) && H < 0.0 && std::isfinite(resid);
  out.H = H;
  out.stationarity = usable ? std::abs(resid / H)
                            : std::numeric_limits<double>::infinity();
  if (prescribed != nullptr) {
    // ⚠️ THE STATIONARITY TEST DOES NOT ROUTE HERE, AND IS STILL WORTH TAKING.
    // What it decides on the solving path -- composite or fallback -- is
    // meaningless at a collar the caller chose: there is no argmax to be pinned
    // against, and differencing the solve would answer about the optimum
    // instead. So `Method` is refused at the R boundary.
    //
    // The NUMBER keeps its meaning and gains a better one: it is how far the
    // point handed over sits from the optimum, in MPa. Reported, and used for
    // exactly one thing -- see the envelope below.
    //
    // ⚠️ NoGradient REACHES THIS PATH TOO, and it is not the solving path's
    // condition. There `usable` also demands `H < 0` -- a MAXIMUM test, which a
    // caller-chosen collar has no business satisfying: a prescribed psi away from
    // the optimum may sit where profit is convex, and that is fine because
    // nothing here divides by `H`. What disables the point is INFEASIBILITY.
    // Almost every such point is already Clamped -- the shut-down state seats a
    // collar of its own choosing -- but a caller can pass exactly that collar
    // back, so "almost" is not a guarantee.
    out.status = clamped ? Status::Clamped
                 : (resid_feasible ? Status::Prescribed : Status::NoGradient);
  } else {
    out.status = !usable ? Status::NoGradient
                 : (out.stationarity > s.stationarity_tol ? Status::Pinned
                                                          : Status::Interior);
  }

  // `status` describes the POINT and is reported whichever route runs;
  // `use_ift` is the route. They differ only when the caller has forced one.
  bool use_ift = prescribed != nullptr
                     ? (!clamped && resid_feasible)
                     : (s.method == Method::Auto
                            ? out.status == Status::Interior
                            : s.method == Method::Ift);
  if (use_ift && prescribed == nullptr && !usable) {
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
    const bool hi_ok = outputs_at(l, psi_star + h_psi, hi, s.curve);
    const bool lo_ok = outputs_at(l, psi_star - h_psi, lo, s.curve);
    if (!hi_ok || !lo_ok) {
      if (prescribed != nullptr) {
        // The same reasoning as the clamp on psi itself, one step out: a psi
        // inside the interval but within h_psi of an end cannot have dY/dpsi
        // centred on it, and a one-sided difference over a shortened interval is
        // exactly what this detector exists to refuse. No fallback to offer --
        // differencing the solve would answer about the optimum -- so the row is
        // reported clamped.
        clamped = true;
        use_ift = false;
        out.status = Status::Clamped;
      } else if (s.method == Method::Ift) {
        util::stop("leaf_gradient(): method = \"ift\" was asked for at a point "
                   "whose feasible collar interval is narrower than one step, "
                   "so dY/dpsi cannot be centred on psi*. Use "
                   "method = \"auto\".");
      } else {
        use_ift = false;
        out.status = Status::Pinned;
      }
    } else {
      for (int j = 0; j < n_outputs; ++j) {
        dY_dpsi[j] = (hi[j] - lo[j]) / (2.0 * h_psi);
      }
      // `dprofit_droot_collar_psi` is EXACT in psi -- forward AD plus the IFT at
      // the ci root-find -- so for profit alone there is something better than a
      // difference of the same quantity, and it is already computed. The other
      // four have no such route and must be differenced. One rule, both paths,
      // which is what keeps a prescribed psi* reproducing the solve exactly.
      //
      // ⚠️ THIS IS THE READER #87 SAID DID NOT EXIST YET. Until the prescribed
      // path landed, the only consumer was `gradient_ift` with `envelope` false --
      // a forced Method::Ift at a pinned point, which throws at all 42 pinned rows
      // of the grid. A prescribed psi away from the optimum is not stationary, so
      // it takes this branch for real, and the exactness now matters.
      dY_dpsi[out_profit] = resid;
    }
  }

  // ⚠️ THE ENVELOPE THEOREM, and the only place this package uses it. At a
  // STATIONARY point dprofit/dpsi is analytically zero, so profit's indirect
  // term vanishes identically and dprofit/dtheta is the direct partial alone.
  // `gradient_ift` is told to ASSIGN that column rather than reach it by
  // multiplying a near-zero dY/dpsi -- the same treatment `collar` gets, for the
  // same reason: an identity is stated, not arrived at.
  //
  // ⚠️ Conditional on stationarity, not on `use_ift`. The identity comes from
  // dprofit/dpsi == 0; at a pinned optimum psi* is a theta-dependent BOUND,
  // dprofit/dpsi is not zero there, and the indirect term survives. Someone
  // forcing Method::Ift there already gets a confidently wrong number and should
  // not get a differently wrong one for this column alone.
  const bool envelope = usable && out.stationarity <= s.stationarity_tol;

  out.used_ift = use_ift;
  if (use_ift) {
    for (int j = 0; j < n_outputs; ++j) {
      out.dY_dpsi[j] = dY_dpsi[j];
    }
    // ⚠️ NULL MEANS TWO DIFFERENT THINGS AND ONLY ONE OF THEM REACHES
    // `gradient_ift`. There, null is "derive dpsi/dtheta by the implicit
    // function theorem", which is right for a collar that was SOLVED for. A
    // prescribed collar with no dpsi_dtheta means the opposite -- it does not
    // move with theta, so the answer is the partial at fixed collar -- and
    // passing the null straight through would silently give the caller the
    // solving path's indirect term on top of it. Zeros are materialised here so
    // that every C++ caller gets that, not just `batch`.
    std::vector<double> zeros;
    const double* d_psi = nullptr;
    if (prescribed != nullptr) {
      d_psi = prescribed->dpsi_dtheta;
      if (d_psi == nullptr) {
        zeros.assign(npars, 0.0);
        d_psi = zeros.data();
      }
    }
    gradient_ift(l, theta, d, single, pars, npars, psi_star, H, dY_dpsi, s,
                 d_psi, envelope, out.M.data(), out.grad.data());
  } else if (prescribed != nullptr) {
    // ⚠️ A CLAMPED PRESCRIBED PSI GETS NO GRADIENT, RATHER THAN THE DIRECT TERM.
    // Not a failure -- the outputs at the clamped collar are perfectly good, and
    // TF24f relies on the clamp to pull an out-of-range tracked state back
    // inside. It is that the derivative is not the one this can compute: the
    // collar used is min(max(psi, a(theta)), b(theta)), so it moves with the
    // BOUND, and dY/dtheta picks up the bound's derivative rather than the
    // caller's dpsi_dtheta. That is the active-set problem arriving through the
    // clamp instead of through the optimiser, and the direct term alone would be
    // plausible and wrong in the documented way.
    //
    // `reset` already left grad, M and dY_dpsi as NA; value, H, stationarity and
    // psi stay, because they describe the point rather than the derivative.
    //
    // The same holds for an INFEASIBLE one, where `dprofit` is a sentinel rather
    // than a derivative -- see the status block above.
  } else {
    gradient_fd(l, theta, d, single, pars, npars, s, out.grad.data());
  }
}

// --- the batch ---------------------------------------------------------------
//
// `theta` is a COLUMN-MAJOR matrix of `theta_nrow` x n_theta, as R hands one
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
//
// `psi` is null to solve, or one collar potential per observation to impose
// (#88); `dpsi_dtheta` is then a COLUMN-MAJOR n x npars matrix, again as R hands
// one over, and null means zero throughout.
inline std::vector<Result> batch(Leaf& l, const double* theta,
                                std::size_t theta_nrow,
                                const std::vector<Drivers>& drivers,
                                bool single, const int* pars,
                                std::size_t npars, const Settings& s,
                                const double* psi = nullptr,
                                const double* dpsi_dtheta = nullptr) {
  const std::size_t n = drivers.size();
  std::vector<Result> out(n);
  double th[n_theta];
  std::vector<double> dpsi(npars);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t row = theta_nrow == 1 ? 0 : i;
    for (int j = 0; j < n_theta; ++j) {
      th[j] = theta[row + std::size_t(j) * theta_nrow];
    }
    Prescribed p;
    if (psi != nullptr) {
      p.psi = psi[i];
      if (dpsi_dtheta != nullptr) {
        for (std::size_t k = 0; k < npars; ++k) {
          dpsi[k] = dpsi_dtheta[i + k * n];
        }
        p.dpsi_dtheta = dpsi.data();
      }
    }
    try {
      at(l, th, drivers[i], single, pars, npars, s, out[i],
         psi == nullptr ? nullptr : &p);
    } catch (const std::exception& e) {
      out[i].status = Status::Error;
      out[i].used_ift = false;
      out[i].message = e.what();
      out[i].grad.assign(npars * n_outputs, util::na_value);
      // Put the leaf back at this row's base parameters before the next one.
      // The next row's own `apply(only = -1)` would do it -- `set_traits` forces
      // the vulnerability rebuild that `perturb_stem_P50` displaced -- but the
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
