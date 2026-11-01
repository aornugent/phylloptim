// -*-c++-*-
#ifndef PHYLLOPTIM_GRADIENT_HPP_
#define PHYLLOPTIM_GRADIENT_HPP_

// Trait gradients composed in C++: R/gradient.R's `.gradient_ift()` and
// `.gradient_fd()` transcribed, plus a loop over observations. From R, dispatch
// rather than the model dominates a gradient's cost, so this crosses the
// boundary once per likelihood evaluation, not once per parameter per
// observation. It returns dY/dtheta for the model parameters; the likelihood and
// the caller's parameterisation Jacobian stay with the caller.
//
// R/gradient.R is the reference, and the two must agree bit for bit: a
// tolerance would let a transcription slip hide inside the solver's ~1e-9
// floor. Three rules make that possible:
//
//   1. R's arithmetic order is kept literally: `-((up - dn) / (2 * h)) / H` is
//      not rewritten as `(dn - up) / (2 * h * H)`.
//   2. No fused multiply-add; see `rounded()`.
//   3. The call order is sequenced explicitly. `f(a) - f(b)` has unspecified
//      operand order in C++ and left-to-right in R, and these `f`s mutate the
//      leaf, so every difference names its two halves first.
//
// `tests/testthat/test-gradient-batch.R` holds this over both supply paths,
// both methods and the pinned and shut-down rows.

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

// The parameter indices, defined in leaf_model.hpp beside the set_traits whose
// argument order they are, reachable under this namespace too.
using phylloptim::n_traits;
using phylloptim::par_vcmax_25;   using phylloptim::par_stem_c;
using phylloptim::par_stem_P50;   using phylloptim::par_root_c;
using phylloptim::par_root_P50;   using phylloptim::par_TF24_beta2;
using phylloptim::par_jmax_25;    using phylloptim::par_phi_J;
using phylloptim::par_curv_fact_elec_trans;
using phylloptim::par_curv_fact_colim;
using phylloptim::par_TF24_cost_scale; using phylloptim::par_R_d_25;
using phylloptim::par_JS22_gamma; using phylloptim::par_CMax_a;
using phylloptim::par_CMax_b;     using phylloptim::par_kmax;
using phylloptim::par_resistance; using phylloptim::par_CF77_lambda;
using phylloptim::par_TF24_floor_lambda_o;
using phylloptim::par_LeastCost_beta;
using phylloptim::par_A;          using phylloptim::par_I;
using phylloptim::par_D;          using phylloptim::par_topt_offset;
using phylloptim::par_dTopt_max;  using phylloptim::par_K_A;
using phylloptim::par_k_i;        using phylloptim::par_k_i_baseline;
using phylloptim::par_k_rec;      using phylloptim::par_k_rec_Ea;
using phylloptim::par_c_acclim_maint;
using phylloptim::par_c_repair_maint; using phylloptim::par_c_protect_maint;
using phylloptim::par_c_repair_flux;
using phylloptim::par_PPFD;

// How long `theta` is: one entry per par_names(), the slots a calibration
// proposes, which come first in the enumeration and stop at par_PPFD.
// ⚠️ Not n_leaf_pars, which also counts par_PPFD (hazard 12).
inline constexpr int n_theta = par_PPFD;

// The proposed slots' names. par_PPFD is seated from the drivers and never
// fitted, so it has none. Append, never reorder (hazard 13);
// test-gradient-batch.R compares R's `.gradient_par_names()` with this list.
inline const std::vector<std::string>& par_names() {
  static const std::vector<std::string> names{
      "vcmax_25",  "stem_c",              "stem_P50",
      "root_c",    "root_P50",            "TF24_beta2",
      "jmax_25",   "phi_J",               "curv_fact_elec_trans",
      "curv_fact_colim", "TF24_cost_scale", "R_d_25",
      "JS22_gamma", "CMax_a", "CMax_b", "g_min",
      "leaf_specific_conductance_max",
      "resistance",
      "CF77_lambda_",
      "TF24_floor_lambda_o",
      "LeastCost_beta",
      "A_", "I_", "D_", "topt_offset_", "dTopt_max_", "K_A_", "k_i_",
      "k_i_baseline_", "k_rec_", "k_rec_Ea_", "c_acclim_maint_",
      "c_repair_maint_", "c_protect_maint_", "c_repair_flux_"};
  return names;
}

// --- the five differentiated outputs ------------------------------------------
//
// A, gc, psi_stem, collar and profit, in R's `.gradient_output_names` order:
// append, never reorder. `collar` is psi* itself, so dcollar/dtheta is
// dpsi*/dtheta. The first four are what a gas-exchange calibration observes;
// `profit` is the carbon plant's mass budget consumes.
inline constexpr int n_outputs = 5;
inline constexpr int out_psi_stem = 2;
inline constexpr int out_collar = 3;
inline constexpr int out_profit = 4;

inline const std::vector<std::string>& output_names() {
  static const std::vector<std::string> names{"A", "gc", "psi_stem", "collar",
                                              "profit"};
  return names;
}

// --- which model's optimum is being differentiated ---------------------------
//
// `kCollar` is the production path, the TF24 cost maximised over the collar
// potential; anything >= 0 is a `Leaf::CostCurve` and selects a stem route
// (psi_stem with the upstream potential pinned at psi_soil). The wrappers below
// only choose a route: changing what a collar arm calls, or its order, breaks
// the bit-for-bit contract with R and the golden files.
inline constexpr int kCollar = -1;

// Seat the model once, at the top of a gradient: the one place the int
// `Settings::curve` meets the Leaf's named API. Seated on the prescribed path
// too, which never calls route_solve, because the stem entry points read the
// seated curve and would otherwise differentiate the leaf's last model.
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

// The same, false where the leaf refuses `psi` as outside its feasible interval.
// Only that refusal is caught: anything else is a failure of the evaluation and
// propagates.
inline bool route_evaluate_feasible(Leaf& l, int curve, double psi) {
  try {
    route_evaluate(l, curve, psi);
  } catch (const util::infeasible_error& e) {
    const std::string token = util::infeasible_token("prescribed_potential");
    if (std::string(e.what()).compare(0, token.size(), token) != 0) {
      throw;
    }
    return false;
  }
  return true;
}

// dprofit with its feasibility: both routes return a hard 0.0 sentinel on their
// shut-down and reversed-gradient exits, which a bare value cannot tell from a
// stationary point (Leaf::dprofit_droot_collar_psi).
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

// The output that is the decision variable, whose derivative is dpsi*/dtheta,
// assigned rather than differenced: `collar` on the collar route, `psi_stem` on
// a stem route, where the collar is never solved for.
inline int route_decision(int curve) {
  return curve < 0 ? out_collar : out_psi_stem;
}

// Read straight off the members, where R reads operating_point_values(), which
// copies these same fields: bit-identical. An output R computes rather than
// copies (uptake, a sum over layers) would need that sum reproduced here.
inline void outputs(const Leaf& l, double* y) {
  y[0] = l.assim_colimited_;
  y[1] = l.stom_cond_CO2_;
  y[2] = l.opt_psi_stem_;
  y[3] = l.opt_root_psi_;
  y[4] = l.profit_;
}

// The outputs with the collar held at `psi`. False where the leaf refuses `psi`
// as infeasible, so a stencil reaching past an end is detected rather than
// differenced one-sidedly over a shorter interval.
inline bool outputs_at(Leaf& l, double psi, double* y, int curve = kCollar) {
  if (!route_evaluate_feasible(l, curve, psi)) {
    return false;
  }
  outputs(l, y);
  return true;
}

// The product, rounded, so that `a + rounded(b * c)` is two IEEE operations and
// never one fused multiply-add. Required: `direct + dY_dpsi * dpsi_dtheta`
// contracts to one `fmadd` on arm64 (clang within an expression, gcc across
// statements), which differs from R's two roundings in the last bit on 28% of
// random triples. A named intermediate does not stop gcc's
// `-ffp-contract=fast`; a volatile store does, because the standard makes it round.
inline double rounded(double x) {
  volatile double v = x;
  return v;
}

// The step: relative to the parameter above 1 and plain `step` below it, so a
// small trait (phi_J = 0.35) is not perturbed into the solve's ~1e-9 noise.
// kmax and resistance take a plain relative step: floored at 1, kmax (3.14e-5)
// would move 3% and measure a secant across visible curvature.
inline double step_for(int par, double value, double step) {
  const double floor =
      (par == par_kmax || par == par_resistance) ? 0.0 : 1.0;
  return std::max(std::abs(value), floor) * step;
}

// Where the stencil is centred: on the value, except that `g_min` is refused
// below zero and defaults to zero, so there the stencil shifts up one step and
// differentiates at g_min = h. The twin of R's `.gradient_centre()`, bit for bit.
inline double centre_for(int par, double value, double h) {
  return (par == par_g_min && value - h < 0.0) ? value + h : value;
}

// --- one observation's drivers ------------------------------------------------
//
// Everything set_physiology takes except kmax and resistance, which come out of
// `theta`. Held C++-side and resolved once, because handing a RootNetwork
// across the R boundary costs more than the gradient it carries.
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

// What the operating point is, not whether the call worked, except `Error`, the
// batch's per-row failure (see `batch`). `Prescribed` and `Refused` end the
// caller-supplied-psi path: the collar was imposed, and for `Refused` it, or
// the stencil about it, fell outside the feasible interval. Append only: a
// consumer holding the integer breaks quietly if a value moves.
enum class Status { Interior, Pinned, NoGradient, Error, Prescribed, Refused };

// No `default:`, so a new Status is a compiler diagnostic rather than silently
// labelled "error".
inline std::string status_name(Status s) {
  switch (s) {
    case Status::Interior:   return "interior";
    case Status::Pinned:     return "pinned";
    case Status::NoGradient: return "no-gradient";
    case Status::Prescribed: return "prescribed";
    case Status::Refused:    return "refused";
    case Status::Error:      return "error";
  }
  return "error";  // unreachable; silences -Wreturn-type on a bad cast
}

struct Settings {
  double step = 1e-6;
  double stationarity_tol = 1e-8;
  Method method = Method::Auto;
  bool fast_stem_curve = true;
  // Which model's optimum to differentiate; `kCollar`, the production path, by
  // default.
  int curve = kCollar;
  // The finite-difference arm's step, separate from `step`: that arm re-solves,
  // so its step must clear the solver's resolution (hazard 15). `step` drives
  // the composite, which perturbs at fixed psi and never re-solves.
  double fd_step = 1e-6;
  // ProfitMax only: re-solve |A|max after each perturbation, because `apply()`
  // clears it; holding it fixed would make the gradient a partial derivative.
  bool resolve_profitmax = false;
};

// A collar potential the caller imposes in place of the solved one, and how it
// responds to each parameter. `dpsi_dtheta` is npars long in `pars` order, not
// n_theta; null means zero (the partial at fixed collar), which is not the
// solving path's "derive it".
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
  // d2profit/dpsi dtheta (npars long) and dY/dpsi at fixed traits: kept, being
  // also the coefficients of a caller's own sensitivity ODE for psi, which
  // cannot be recovered from -M/H.
  std::vector<double> M;
  double dY_dpsi[n_outputs];
  // The collar the outputs were evaluated at -- psi* when solved, the requested
  // one when prescribed (refused or not).
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
// ⚠️ Order is load-bearing (hazard 10): set_traits() returns the leaf to its
// just-constructed state, so the drivers are re-supplied after it, and
// set_physiology() re-derives vcmax_/jmax_/R_d_ behind the temperature cache
// set_traits() cleared. The other order runs the whole gradient at the first
// vcmax the object saw. `only` is the single parameter that moved, or -1 for all.
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single,
                  int only, bool fast_stem_curve) {
  // The stem_P50 fast path: the stem table depends on stem_c alone, so
  // perturb_stem_P50 moves the trait and its derived potentials and leaves the
  // physiology seated. Sound only because `only` names a single parameter, so
  // the rest of `theta` is what the object was last set to; the restore through
  // set_traits recomputes stem_b by the same expression, to the bit.
  if (fast_stem_curve && only == par_stem_P50) {
    l.perturb_stem_P50(theta[par_stem_P50]);
    return;
  }
  l.set_traits(theta[par_vcmax_25], theta[par_stem_c], theta[par_stem_P50],
               theta[par_root_c], theta[par_root_P50],
               theta[par_TF24_beta2], theta[par_jmax_25],
               theta[par_phi_J], theta[par_curv_fact_elec_trans],
               theta[par_curv_fact_colim], theta[par_TF24_cost_scale],
               theta[par_R_d_25], theta[par_JS22_gamma],
               theta[par_CMax_a], theta[par_CMax_b], theta[par_g_min]);
  // Not traits and nothing is derived from them, so written directly. Only a
  // set slot is written, as R does (`!is.na`): on a route that does not own a
  // slot it is NA, and writing it would overwrite the leaf's aligned default.
  if (!std::isnan(theta[par_CF77_lambda])) {
    l.CF77_lambda_ = theta[par_CF77_lambda];
  }
  if (!std::isnan(theta[par_TF24_floor_lambda_o])) {
    l.TF24_floor_lambda_o = theta[par_TF24_floor_lambda_o];
  }
  if (!std::isnan(theta[par_LeastCost_beta])) {
    l.LeastCost_beta = theta[par_LeastCost_beta];
  }
  // The thermal-damage layer's inputs, ONLY with its gate on -- off, the slots
  // hold whatever a caller zero-filled and are never read -- and BEFORE the
  // drivers, because set_physiology is what seats the block that reads them.
  if (l.use_thermal_damage_) {
    l.A_ = theta[par_A];
    l.I_ = theta[par_I];
    l.D_ = theta[par_D];
    l.topt_offset_ = theta[par_topt_offset];
    l.dTopt_max_ = theta[par_dTopt_max];
    l.K_A_ = theta[par_K_A];
    l.k_i_ = theta[par_k_i];
    l.k_i_baseline_ = theta[par_k_i_baseline];
    l.k_rec_ = theta[par_k_rec];
    l.k_rec_Ea_ = theta[par_k_rec_Ea];
    l.c_acclim_maint_ = theta[par_c_acclim_maint];
    l.c_repair_maint_ = theta[par_c_repair_maint];
    l.c_protect_maint_ = theta[par_c_protect_maint];
    l.c_repair_flux_ = theta[par_c_repair_flux];
  }
  if (single) {
    // R's `series_resistance()`: one series resistance in `r_R_V_sum`, built
    // here rather than copied from `d.root_network` so the two implementations
    // cannot disagree about the other four fields.
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

// Both loops below take every parameter's gradient from the base point. A
// setter normally restores everything through set_traits() + set_physiology(),
// so base is restored only at the end; a shortcut (perturb_stem_P50) restores
// nothing and is sound only from base, so one that follows another parameter
// restores base first. Cheap: set_traits() rebuilds a table only for a new `c`.
// Asked as "takes a shortcut", not "is stem_P50", so a second shortcut is covered.
inline bool takes_shortcut(int par, const Settings& s) {
  return s.fast_stem_curve && par == par_stem_P50;
}

// The implicit-function composite: per parameter, two perturbed evaluations at
// the unperturbed psi*, neither re-solving; `dprofit` there gives the mixed
// partial M and the outputs the direct term. `dpsi_dtheta` null means derive it
// as -M/H (a solved collar); non-null is an imposed collar's own. One function
// for both, so their agreement at psi* is a real assertion.
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
    const double centre = centre_for(p, theta[p], h);
    for (int side = 0; side < 2; ++side) {
      // Up first, then down: R evaluates `up <- side(1)` before `dn <- side(-1)`
      // and both mutate the leaf.
      std::copy(theta, theta + n_theta, th);
      th[p] = side == 0 ? centre + h : centre - h;
      apply(l, th, d, single, p, s.fast_stem_curve);
      double* dst = side == 0 ? up : dn;
      // Evaluate first, then read dprofit at the same fixed collar -- R's order,
      // and `evaluate_root_collar_psi` is what seats the state `dprofit` reads.
      if (s.resolve_profitmax) {
        // apply() cleared ProfitMax's normaliser. Re-solved, not pinned, so the
        // derivative is total; must match R/gradient.R's `pinned$reseat`.
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
      // Again before the derivative read, because R re-seats inside every
      // accessor and prepare_profitmax() mutates the leaf's rate fields.
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
    // The decision variable is psi* held fixed, so its direct term is zero and
    // the composite is dpsi*/dtheta: assigned, not differenced.
    out[k * n_outputs + route_decision(s.curve)] = d_psi;
    // The envelope theorem, assigned for the same reason: profit's indirect term
    // is identically zero at a stationary point. Also immune to a non-finite
    // caller-supplied d_psi, where 0 * x would be NaN.
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
    const double centre = centre_for(p, theta[p], h);
    for (int side = 0; side < 2; ++side) {
      std::copy(theta, theta + n_theta, th);
      th[p] = side == 0 ? centre + h : centre - h;
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

// Both routes carry the thermal-damage layer: they difference the solve's own
// dprofit and outputs, needing only its inputs in `theta` (seated by apply()).
// The derivative surface (operating_collar, billed_outputs_at) refuses it,
// because the parameter pack's kernels do not apply the layer.

// --- one observation ---------------------------------------------------------
//
// Fills `out` progressively and throws on the conditions leaf_gradient() stops
// on, so a partly determined row keeps its diagnostics when `batch` catches.
inline void at(Leaf& l, const double* theta, const Drivers& d, bool single,
               const int* pars, std::size_t npars, const Settings& s,
               Result& out, const Prescribed* prescribed = nullptr) {
  out.reset(npars);

  apply(l, theta, d, single, -1, s.fast_stem_curve);
  // Seated here, not only by `batch`, because at() is a public entry point.
  route_seat(l, s.curve);

  // The default solves for the collar; `prescribed` imposes one. Either way
  // `psi_star` is the collar the outputs were evaluated at, not always the argmax.
  if (prescribed != nullptr) {
    // A prescribed psi outside the feasible interval is reported, not thrown: a
    // fit visits such points routinely, and `status` tells a batch which rows to
    // distrust. Everything but `psi` stays NA.
    if (!route_evaluate_feasible(l, s.curve, prescribed->psi)) {
      out.status = Status::Refused;
      out.psi = prescribed->psi;
      return;
    }
  } else {
    route_solve(l, s.curve);
  }
  const double psi_star = route_psi(l, s.curve);
  out.psi = psi_star;
  outputs(l, out.value);

  // The composite rests on stationarity, which fails at a pinned optimum, and
  // fails quietly (O(1) where the truth is ~1e-8). So it is tested, by the
  // implied Newton step |dprofit(psi*) / H| in MPa: over the 288-point grid
  // interior points reach 4.8e-11 and pinned ones start at 6.3e-06, so the 1e-8
  // default sits in an empty band.
  const double h_psi = std::max(std::abs(psi_star), 1.0) * s.step;
  // With its feasibility: the prescribed path does not divide by H, so it would
  // adopt dprofit's 0.0 sentinel as dprofit/dpsi.
  bool resid_feasible = false;
  const double resid =
      route_dprofit_checked(l, s.curve, psi_star, &resid_feasible);
  // Named halves (rule 3 in the header).
  const double d_hi = route_dprofit(l, s.curve, psi_star + h_psi);
  const double d_lo = route_dprofit(l, s.curve, psi_star - h_psi);
  const double H = (d_hi - d_lo) / (2.0 * h_psi);
  // H == 0 with resid == 0 is the shut-down sentinel, and H > 0 is not a
  // maximum: either way the composite has nothing to stand on.
  const bool usable = std::isfinite(H) && H < 0.0 && std::isfinite(resid);
  out.H = H;
  out.stationarity = usable ? std::abs(resid / H)
                            : std::numeric_limits<double>::infinity();
  if (prescribed != nullptr) {
    // At an imposed collar the stationarity test routes nothing (`Method` is
    // refused at the R boundary); it is reported as the point's distance from
    // the optimum in MPa and decides the envelope below. NoGradient here means
    // infeasible, not H >= 0: an imposed collar may sit where profit is convex.
    out.status = resid_feasible ? Status::Prescribed : Status::NoGradient;
  } else {
    out.status = !usable ? Status::NoGradient
                 : (out.stationarity > s.stationarity_tol ? Status::Pinned
                                                          : Status::Interior);
  }

  // `status` describes the POINT and is reported whichever route runs;
  // `use_ift` is the route. They differ only when the caller has forced one.
  bool use_ift = prescribed != nullptr
                     ? resid_feasible
                     : (s.method == Method::Auto
                            ? out.status == Status::Interior
                            : s.method == Method::Ift);
  if (use_ift && prescribed == nullptr && !usable) {
    util::stop("leaf_gradient(): method = \"ift\" was asked for at a point with "
               "no usable curvature (H = " + util::to_string(H) + "), so -M/H "
               "has nothing to stand on. This is a shut-down or otherwise "
               "determined operating point; use method = \"auto\".");
  }
  // Refused by name, as R/gradient.R does and in the same place: the step check
  // below does not catch every pinned row.
  if (use_ift && prescribed == nullptr && out.status == Status::Pinned) {
    util::stop("leaf_gradient(): method = \"ift\" was asked for at a pinned "
               "optimum, where psi* is held by a bound and the "
               "implicit-function composite does not apply. Use "
               "method = \"auto\".");
  }

  double dY_dpsi[n_outputs];
  if (use_ift) {
    // dY/dpsi at fixed traits. It catches some pinned optima (psi* within h_psi
    // of a bound), but nothing relies on that.
    double hi[n_outputs];
    double lo[n_outputs];
    // Both before either is tested, as in R: each call moves the leaf.
    const bool hi_ok = outputs_at(l, psi_star + h_psi, hi, s.curve);
    const bool lo_ok = outputs_at(l, psi_star - h_psi, lo, s.curve);
    if (!hi_ok || !lo_ok) {
      if (prescribed != nullptr) {
        // A psi within h_psi of an end cannot have dY/dpsi centred on it, and
        // differencing the solve would answer about the optimum, so the row is
        // refused, with the point's value but no gradient.
        use_ift = false;
        out.status = Status::Refused;
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
      // dprofit is exact in psi (forward AD plus the IFT at the ci root-find) and
      // already computed, so profit's dY/dpsi is that, not a difference. One rule
      // on both paths keeps a prescribed psi* reproducing the solve exactly; it
      // matters where the envelope does not apply, as off the optimum.
      dY_dpsi[out_profit] = resid;
    }
  }

  // The envelope theorem: at a stationary point dprofit/dpsi is zero, so
  // gradient_ift assigns profit's column the direct partial alone. Conditional
  // on stationarity, not on use_ift: at a pinned optimum psi* is a bound and the
  // indirect term survives.
  const bool envelope = usable && out.stationarity <= s.stationarity_tol;

  out.used_ift = use_ift;
  if (use_ift) {
    for (int j = 0; j < n_outputs; ++j) {
      out.dY_dpsi[j] = dY_dpsi[j];
    }
    // To gradient_ift a null dpsi_dtheta means "derive it", right only for a
    // solved collar. A prescribed collar with none does not move with theta, so
    // zeros are materialised here, for every C++ caller and not just `batch`.
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
    // A prescribed psi within one step of an end, or infeasible, gets no
    // gradient: `reset` left grad, M and dY_dpsi NA, while value, H,
    // stationarity and psi describe the point and stay.
  } else {
    gradient_fd(l, theta, d, single, pars, npars, s, out.grad.data());
  }
}

// --- the batch ---------------------------------------------------------------
//
// `theta` is a column-major `theta_nrow` x n_theta matrix, as R hands one over:
// one row per observation, or one row shared by all.
//
// Per-row status, not an exception: a proposal during a fit will reach points
// the solve cannot handle, and that must cost those rows, not the likelihood
// evaluation. A failed row's gradient is all NA rather than partly filled, so a
// caller checking status per row cannot read a gradient with a hole in it; its
// value, H and stationarity stay where they were determined.
//
// `psi` is null to solve, or one collar per observation to impose;
// `dpsi_dtheta` is then a column-major n x npars matrix, null meaning zero.
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
      // Back to this row's base: the last row has no next apply() to do it, and
      // a batch ending on the fast path would hand back a leaf on a rescaled stem
      // curve (hazard 8).
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
