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
  // dY/dpsi at the base point: the composite's psi-channel, and the third of its
  // three ingredients alongside `H` and the per-parameter M. Reported for the
  // same reason those two are -- it is what the composite stands on -- and NA
  // unless `used_ift`, because the fallback never forms it. `transpose_at` below
  // is its first consumer.
  double dY_dpsi[n_outputs];
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
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve) {
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
      std::copy(theta, theta + n_pars, th);
      th[p] = side == 0 ? theta[p] + h : theta[p] - h;
      apply(l, th, d, single, p, s.fast_stem_curve);
      double* dst = side == 0 ? up : dn;
      // Evaluate first, then read dprofit at the same fixed collar -- R's order,
      // and `evaluate_root_collar_psi` is what seats the state `dprofit` reads.
      if (!outputs_at(l, psi_star, dst + 1)) {
        util::stop("leaf_gradient(): perturbing `" + par_names()[std::size_t(p)] +
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
  bool at_base = true;
  for (std::size_t k = 0; k < npars; ++k) {
    const int p = pars[k];
    if (takes_shortcut(p, s) && !at_base) {
      apply(l, theta, d, single, -1, s.fast_stem_curve);
    }
    at_base = false;
    const double h = step_for(p, theta[p], s.step);
    for (int side = 0; side < 2; ++side) {
      std::copy(theta, theta + n_pars, th);
      th[p] = side == 0 ? theta[p] + h : theta[p] - h;
      apply(l, th, d, single, p, s.fast_stem_curve);
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
  // mildest pinned one 6.3e-06, so the 1e-08 default sits in an empty band four
  // orders wide on each side.
  const double h_psi = std::max(std::abs(psi_star), 1.0) * s.step;
  const double resid = l.dprofit_droot_collar_psi(psi_star);
  // Named halves: `f(a) - f(b)` has unspecified operand order in C++ and
  // left-to-right order in R, and `dprofit_droot_collar_psi` mutates the leaf.
  const double d_hi = l.dprofit_droot_collar_psi(psi_star + h_psi);
  const double d_lo = l.dprofit_droot_collar_psi(psi_star - h_psi);
  const double H = (d_hi - d_lo) / (2.0 * h_psi);
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
  double* const dY_dpsi = out.dY_dpsi;
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
// reported at all -- which is the case here, `n_outputs` being four.
//
// ⚠️ THE OBJECTIVE'S PSI-CHANNEL MUST BE EXCLUDED FROM `s`. At an interior
// optimum dprofit/dpsi = 0 by the envelope theorem, so profit's psi-channel
// contributes nothing and adding a measured central difference of it would
// double-count -- and that difference is not small noise either, being a
// difference of a FLAT maximum. #8 appends `profit` as a fifth output and gives
// it exactly that special case in the forward composite; when it lands this
// constant becomes `out_profit` and nothing else here changes. The dot-product
// identity is what makes that one line rather than a hazard: it compares the
// transpose against the forward path column by column, so a forward path that
// special-cases profit and a transpose that does not fails immediately.
inline constexpr int out_objective = -1;

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
inline double psi_channel(int j, const double* dY_dpsi) {
  if (j == out_collar) {
    return 1.0;
  }
  if (j == out_objective) {
    return 0.0;
  }
  return dY_dpsi[j];
}

struct TransposeResult {
  // The four outputs the transpose is taken at, and the operating-point
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

  double th[n_pars];
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

    double up[1 + n_outputs];
    double dn[1 + n_outputs];
    for (std::size_t k = 0; k < npars; ++k) {
      // The base-point invariant and the shortcut restore are `gradient_ift`'s,
      // for `gradient_ift`'s reasons -- see `takes_shortcut`.
      const int p = pars[k];
      if (takes_shortcut(p, s) && !at_base) {
        apply(l, theta, d, single, -1, s.fast_stem_curve);
      }
      at_base = false;
      const double h = step_for(p, theta[p], s.step);
      for (int side = 0; side < 2; ++side) {
        std::copy(theta, theta + n_pars, th);
        th[p] = side == 0 ? theta[p] + h : theta[p] - h;
        apply(l, th, d, single, p, s.fast_stem_curve);
        double* dst = side == 0 ? up : dn;
        if (!outputs_at(l, psi_star, dst + 1)) {
          util::stop("leaf_gradient_transpose(): perturbing `" +
                     par_names()[std::size_t(p)] +
                     "` moved the feasible collar interval past psi*, so the "
                     "operating point could not be evaluated there. This point "
                     "is on an active-set boundary; lower `stationarity_tol` or "
                     "difference the solve directly.");
        }
        dst[0] = l.dprofit_droot_collar_psi(psi_star);
      }
      const double M = (up[0] - dn[0]) / (2.0 * h);
      // v . dY/dtheta|_psi. `collar` is skipped rather than summed: its direct
      // term is zero by construction, `outputs_at` having just asserted that
      // both sides sit at exactly psi*.
      double row = 0.0;
      for (int j = 0; j < n_outputs; ++j) {
        if (j == out_collar) {
          continue;
        }
        row += rounded(v[j] * ((up[1 + j] - dn[1 + j]) / (2.0 * h)));
      }
      out.adjoint[k] = row + rounded(m * M);
    }
  } else {
    // The fallback, transposed: difference the solve and contract. No output is
    // exceptional here, `collar` included, because nothing is being composed --
    // which is exactly why `gradient_fd` has no special case either.
    double up[n_outputs];
    double dn[n_outputs];
    for (std::size_t k = 0; k < npars; ++k) {
      const int p = pars[k];
      if (takes_shortcut(p, s) && !at_base) {
        apply(l, theta, d, single, -1, s.fast_stem_curve);
      }
      at_base = false;
      const double h = step_for(p, theta[p], s.step);
      for (int side = 0; side < 2; ++side) {
        std::copy(theta, theta + n_pars, th);
        th[p] = side == 0 ? theta[p] + h : theta[p] - h;
        apply(l, th, d, single, p, s.fast_stem_curve);
        l.find_root_collar_psi();
        outputs(l, side == 0 ? up : dn);
      }
      double row = 0.0;
      for (int j = 0; j < n_outputs; ++j) {
        row += rounded(v[j] * ((up[j] - dn[j]) / (2.0 * h)));
      }
      out.adjoint[k] = row;
    }
  }
  apply(l, theta, d, single, -1, s.fast_stem_curve);
}

}  // namespace gradient
}  // namespace phylloptim

#endif
