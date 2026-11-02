// -*-c++-*-
#ifndef PHYLLOPTIM_CLOSED_FORM_HPP_
#define PHYLLOPTIM_CLOSED_FORM_HPP_

// The closed form: an alternative method for the stem-route optimum, seated with
// `set_model(curve, "stem", "closed")` and run by `optimise()`. It exists for
// speed; `bench_solve` prices it against the exact stem solve, fallbacks included.
//
// Every curve here satisfies dA/dE = lambda at the optimum, and given lambda the
// optimum is the Medlyn USO form
//
//     ci/ca = xi/(xi + sqrt(D)),     xi = sqrt(Q/lambda),     Q from `uso_group`
//
// Where lambda is a power law psi^n at the wet end, the circularity inverts
// explicitly, and one Newton step on the supply-minus-demand residual finishes.
// That needs the benefit link's h'(A) (`Leaf::BenefitLink`) to be a solve-time
// constant: true of Identity (TF24, CF77, JS22, CMax, TF24_floor) and of Scaled
// (ProfitMax, whose |A|max is fixed per driver set); false of Log (SOX, JW26,
// LeastCost), where lambda carries A and the problem is a fixed point in A.
// Implemented: TF24 (n = stem_c*TF24_beta2 - 1), CF77 (n = 0) and TF24_floor
// (n = 0, its wet-end lambda being lambda_o). `optimise_into` refuses the rest.
//
// The USO collapse, not n, limits accuracy: its 3*Gstar is the electron-transport
// dA/dci and it assumes ci >> Gstar. Given the exact lambda at the exact optimum,
// ci is off by 1-2% at 25 C with ci/ca > 0.6, 4-10% at 40 C and up to 139% at
// ci/ca 0.16, which is why the guard is on ci/ca. A Newton polish in ci would
// remove this but is a second solver (`.claude/CLAUDE.md`, "One solver").
//
// Five things not to get wrong:
//
//  1. `newton_steps = 1` is deliberate: two steps are worse in the tail.
//  2. The guard (`within_guard`) tests an output, so it runs after solving and a
//     failure falls back to the exact solve; `optimise_into` counts fallbacks.
//  3. Quote the realised speedup 1/[phi + (1-phi)/S], with phi from
//     `Leaf::closed_form_fallback_fraction()`, never the ceiling S.
//  4. ⚠️ Hazard 3: the argmax jumps where a smooth input carries the leaf across
//     the guard, so this method must not sit under a differentiated path, and
//     `gradient::route_seat` seats the exact solve. Over a 201-point +/-30% `kmax`
//     sweep the mean |second difference| of `opt_psi_stem_` is 1.2x the exact
//     solve's where the sweep never falls back, and 100-153x where it crosses.
//  5. psi* is not E/kmax: the supply is the vulnerability integral, so
//     psi* = G^-1(E/kmax + G(psi_upstream)), which is `psi_from_supply`.
//
// Refused rather than approximated (`Leaf::optimise_closed`): more than one layer,
// the collar route, and the energy balance, under which D = `vpd_leaf_` depends on
// the E being solved for. `within_guard` cannot see that wrong input (A overstated
// by up to 178% at 40 C with every row passing), so each solver stops on it too.

#include <phylloptim/constants.hpp>
#include <phylloptim/leaf_model.hpp>

#include <algorithm>
#include <cmath>

namespace phylloptim {
namespace closed_form {

// The group converting a marginal cost of water into the USO slope:
//     xi = sqrt(Q/lambda),   Q = 3*Gstar*kg_to_mol_h2o/(H2O_CO2_stom_diff_ratio_*1e-3)
// with Gstar in Pa and lambda in umol CO2 (kg H2O)^-1. The 3*Gstar is the
// electron-transport-limited dA/dci, so under Rubisco-bound colimitation the USO
// form approximates the optimum even where n = 0 (CF77).
inline double uso_group(const Leaf &l) {
  const double gstar_Pa = l.gamma_ * l.umol_per_mol_to_Pa_;
  return 3.0 * gstar_Pa * kg_to_mol_h2o / (l.H2O_CO2_stom_diff_ratio_ * 1e-3);
}

// Transpiration implied by assimilation on the demand side, kg H2O m-2 s-1; the
// closed form drives its difference from `Leaf::transpiration` to zero. Reads
// `vpd_leaf_`, which is `atm_vpd_` on every path that reaches here.
inline double transpiration_from_assim(const Leaf &l, double assim, double ci) {
  return l.H2O_CO2_stom_diff_ratio_ * 1e-3 * assim * l.vpd_leaf_ /
         ((l.ca_ - ci) * kg_to_mol_h2o);
}

// Stomatal conductance to CO2 implied by a transpiration, mol CO2 m-2 s-1: the
// expression `Leaf::stom_cond_CO2` uses, so the two cannot drift.
inline double stom_cond_from_E(const Leaf &l, double E) {
  return l.atm_kpa_ * E * kg_to_mol_h2o / l.vpd_leaf_ /
         l.H2O_CO2_stom_diff_ratio_;
}

// d(lambda)/d(psi) of `Leaf::lambda_TF24`, analytic:
//     lambda = K * (1-f)^(TF24_beta2-1) * p^(stem_c-1),  p = psi/stem_b,  f = exp(-p^stem_c)
// It is TF24_floor's too, since lambda_o differentiates away. JS22, CMax and
// ProfitMax each need their own before the closed form can serve them.
inline double dlambda_TF24(const Leaf &l, double psi) {
  const double K =
      l.TF24_cost_scale * l.TF24_beta2 * l.stem_c / (l.stem_b * l.leaf_specific_conductance_max_);
  const double p = psi / l.stem_b;
  const double f = std::exp(-std::pow(p, l.stem_c));
  return K * std::pow(p, l.stem_c - 2.0) *
         ((l.TF24_beta2 - 1.0) * std::pow(1.0 - f, l.TF24_beta2 - 2.0) * f * l.stem_c *
              std::pow(p, l.stem_c) +
          std::pow(1.0 - f, l.TF24_beta2 - 1.0) * (l.stem_c - 1.0)) /
         l.stem_b;
}

// dA/dci of the colimitation quadratic, analytic.
inline double dassim_dci(const Leaf &l, double ci, double electron_transport) {
  const double gstar = l.gamma_ * l.umol_per_mol_to_Pa_;
  const double ar = l.vcmax_ * (ci - gstar) / (ci + l.km_);
  const double ae =
      electron_transport / 4.0 * (ci - gstar) / (ci + 2.0 * gstar);
  const double dar =
      l.vcmax_ * (l.km_ + gstar) / ((ci + l.km_) * (ci + l.km_));
  const double dae = electron_transport / 4.0 * 3.0 * gstar /
                     ((ci + 2.0 * gstar) * (ci + 2.0 * gstar));
  const double s = ar + ae, ds = dar + dae;
  const double cv = l.curv_fact_colim;
  const double disc = std::sqrt(s * s - 4.0 * cv * ar * ae);
  const double ddisc =
      (2.0 * s * ds - 4.0 * cv * (dar * ae + ar * dae)) / (2.0 * disc);
  return (ds - ddisc) / (2.0 * cv);
}

// The energy-balance refusal, called by `Leaf::optimise_closed` and by each solver
// so that a direct call cannot skip it. A stop rather than a NaN: a NaN fails
// `within_guard`, which falls back to the exact solve and hides the cause.
inline void require_prescribed_temperature(const Leaf &l, const char *who) {
  if (l.use_energy_balance_) {
    util::stop(std::string(who) +
               " is refused with use_energy_balance_ on, and this is a statement "
               "about the form rather than about its accuracy (#116). Fick's law "
               "divides by vpd_leaf_, which depends on Tleaf, which depends on E "
               "-- and E is what the inversion solves for. So the deficit D is not "
               "known at the moment of inversion and the form is not closed; "
               "forced through, it overstates assimilation by up to 178% and "
               "within_guard does not see it. Use method \"exact\", or prescribe "
               "the leaf temperature by setting use_energy_balance_ to false.");
  }
}

struct Solution {
  double psi_stem;      // MPa, positive magnitude
  double ci;            // Pa
  double assim;         // umol CO2 m-2 s-1
  double transpiration; // kg H2O m-2 s-1
  double stom_cond_CO2; // mol CO2 m-2 s-1
  double g1_eff;        // the USO slope xi, kPa^0.5
};

// The stem potential that carries transpiration E from psi_upstream (note 5):
//
//     E = kmax * [G(psi) - G(psi_upstream)]   =>   psi = G^-1(E/kmax + G(psi_up))
//
// A demand above what the stem carries at psi_crit returns psi_crit, which
// `within_guard` rejects, handing the row to the exact solve.
inline double psi_from_supply(Leaf &l, double E, double psi_upstream) {
  if (!(E > 0.0)) {
    return psi_upstream;
  }
  if (!(E < l.transpiration(l.psi_crit, psi_upstream))) {
    return l.psi_crit;
  }
  return l.transpiration_to_psi_stem(E, psi_upstream);
}

// The outputs implied by a stem potential and the lambda there. The curve enters
// only through `lambda`, so one body serves every curve with a constant h'.
inline Solution evaluate_at_lambda(Leaf &l, double psi, double lambda, double Q,
                                   double sqrt_D) {
  const double xi = std::sqrt(Q / lambda);
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  return Solution{psi, ci, assim, E, stom_cond_from_E(l, E), xi};
}

// `evaluate_at_lambda` at TF24's lambda.
inline Solution evaluate_at(Leaf &l, double psi, double Q, double sqrt_D) {
  return evaluate_at_lambda(l, psi, l.lambda_TF24(psi), Q, sqrt_D);
}

// The wet-end starting potential in units of stem_b: the root of
//     h(p) = (p - p_up) * p^(n/2) - C
// the power-law demand against a supply linearised at the wet end. The p_up term
// is psi_soil on the stem route, and dropping it overstates the supply by the
// whole wet-end integral. At p_up == 0 the root is C^(2/(n+2)) exactly.
inline double wet_end_p(double C, double p_up, double n, double p_crit) {
  const double p0 = std::pow(C, 2.0 / (n + 2.0));
  if (!(p_up > 0.0)) {
    return p0;
  }
  // Newton from above the root: h increases on p > p_up for n >= 0, so the
  // iteration walks down. Four fixed steps suffice for this scalar equation.
  double p = p_up + p0;
  for (int i = 0; i < 4; ++i) {
    const double q = std::pow(p, n / 2.0);
    const double h = (p - p_up) * q - C;
    const double dh = q * (1.0 + 0.5 * n * (p - p_up) / p);
    if (!(dh > 0.0) || !std::isfinite(h)) {
      break;
    }
    const double next = p - h / dh;
    if (!std::isfinite(next)) {
      break;
    }
    // Strictly inside (p_up, p_crit): outside it there is no operating point, and
    // `pow` of a negative base is a NaN that would propagate silently.
    p = std::min(std::max(next, p_up * (1.0 + 1e-12) + 1e-12), p_crit);
  }
  return p;
}

// General TF24_beta2: power-law leading order, then `newton_steps` Newton steps on
// the supply-minus-demand residual. Requires set_physiology. Leave newton_steps
// at 1 (note 1).
inline Solution solve(Leaf &l, double psi_upstream, int newton_steps = 1) {
  require_prescribed_temperature(l, "closed_form::solve");
  const double kmax = l.leaf_specific_conductance_max_;
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double Q = uso_group(l);
  const double K_lambda = l.TF24_cost_scale * l.TF24_beta2 * l.stem_c / (l.stem_b * kmax);
  const double Xi = std::sqrt(Q / K_lambda);
  const double n = l.stem_c * l.TF24_beta2 - 1.0;
  const double electron_transport = l.electron_transport();

  // Leading order, taking kappa from the wet-end limit A(ci -> ca).
  const double assim_wet = l.assim_colimited(l.ca_ * (1.0 - 1e-9));
  const double kappa = l.H2O_CO2_stom_diff_ratio_ * 1e-3 * assim_wet /
                       (l.ca_ * kg_to_mol_h2o);
  const double C = kappa * sqrt_D * Xi / (kmax * l.stem_b);
  const double p_crit = l.psi_crit / l.stem_b;
  double p = wet_end_p(C, psi_upstream / l.stem_b, n, p_crit);

  for (int k = 0; k < newton_steps; ++k) {
    const double psi = p * l.stem_b;
    // The bracket, not just positivity: below psi_upstream no water moves and the
    // residual's supply term is negative, which Newton would chase outward.
    if (!(psi > psi_upstream) || psi >= l.psi_crit) {
      break;
    }
    const double lambda = l.lambda_TF24(psi);
    const double xi = std::sqrt(Q / lambda);
    const double ci = l.ca_ * xi / (xi + sqrt_D);
    const double assim = l.assim_colimited(ci);
    const double dassim = dassim_dci(l, ci, electron_transport);
    const double u = l.ca_ - ci;
    const double E = transpiration_from_assim(l, assim, ci);
    const double dE_dci = l.H2O_CO2_stom_diff_ratio_ * 1e-3 * l.vpd_leaf_ /
                          kg_to_mol_h2o * (dassim * u + assim) / (u * u);
    const double dci_dxi = l.ca_ * sqrt_D / ((xi + sqrt_D) * (xi + sqrt_D));
    const double dxi_dpsi = -0.5 * xi / lambda * dlambda_TF24(l, psi);
    // Residual: hydraulic supply minus stomatal demand, kg H2O m-2 s-1. Measuring
    // the supply from psi_upstream moves its value but not its slope kmax*f(psi).
    const double R = l.transpiration(psi, psi_upstream) - E;
    const double dR =
        kmax * std::exp(-std::pow(p, l.stem_c)) - dE_dci * dci_dxi * dxi_dpsi;
    if (dR == 0.0 || !std::isfinite(R) || !std::isfinite(dR)) {
      break;
    }
    double psi_next = psi - R / dR;
    psi_next = std::max(psi_upstream,
                        std::min(psi_next, l.psi_crit * (1.0 - 1e-9)));
    p = psi_next / l.stem_b;
  }

  return evaluate_at(l, p * l.stem_b, Q, sqrt_D);
}

// TF24_beta2 == 1/stem_c: lambda is independent of psi (n = 0), so xi is constant
// and only the supply inversion remains. An algebraic special case and a test of
// the n = 0 arm, not a preferred value: 1/stem_c is 0.373 at the default stem_c.
inline Solution solve_exact_beta2(Leaf &l, double psi_upstream) {
  require_prescribed_temperature(l, "closed_form::solve_exact_beta2");
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double Q = uso_group(l);
  const double xi = std::sqrt(Q * l.stem_b * l.leaf_specific_conductance_max_ /
                              (l.TF24_cost_scale * l.TF24_beta2 * l.stem_c));
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  // The potential comes from the supply; `opt_psi_stem_` is set from it.
  return Solution{psi_from_supply(l, E, psi_upstream),
                  ci, assim, E, stom_cond_from_E(l, E), xi};
}

inline bool beta2_is_exact(const Leaf &l, double tol = 1e-12) {
  return std::abs(l.TF24_beta2 - 1.0 / l.stem_c) <= tol * std::max(1.0, 1.0 / l.stem_c);
}

// CF77, n = 0: lambda is `cf77_price()`, so xi, ci, A and E follow in one pass and
// the potential from the supply inversion. No iteration is not exactness: the USO
// collapse (`uso_group`) still approximates ci under Rubisco-bound colimitation.
inline Solution solve_CF77(Leaf &l, double psi_upstream) {
  require_prescribed_temperature(l, "closed_form::solve_CF77");
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double xi = std::sqrt(uso_group(l) / l.cf77_price());
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  return Solution{psi_from_supply(l, E, psi_upstream),
                  ci, assim, E, stom_cond_from_E(l, E), xi};
}

// TF24_floor: the leading order is CF77's at price lambda_o, because the TF24 term
// vanishes at the wet end while stem_c*TF24_beta2 > 1 (below that it diverges and
// this start is wrong; see `Leaf::lambda_TF24`). That start is the lowest price,
// so it lands too dry, and the Newton step walks back toward the soil restoring
// the TF24 term. The start is clamped strictly inside psi_crit, where the loop
// would take no step and the row would fall back.
//
// Its tail is worse than TF24's and CF77's, and `within_guard` cannot see it: ci/ca
// comes from the closed form's own lambda. Worst served |psi* - psi*_exact|/psi*_exact
// on a 150-row grid at 25 C: TF24 5.3e-02, CF77 2.4e-01, and TF24_floor 2.9e-01 at
// TF24_cost_scale 0.5 but 5.2e-01 at the default 7.5.
inline Solution solve_TF24_floor(Leaf &l, double psi_upstream,
                              int newton_steps = 1) {
  require_prescribed_temperature(l, "closed_form::solve_TF24_floor");
  const double kmax = l.leaf_specific_conductance_max_;
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double Q = uso_group(l);
  const double electron_transport = l.electron_transport();

  // Leading order: lambda == lambda_o, i.e. `solve_CF77` at this curve's price.
  const double xi0 = std::sqrt(Q / l.TF24_floor_lambda_o);
  const double ci0 = l.ca_ * xi0 / (xi0 + sqrt_D);
  const double assim0 = l.assim_colimited(ci0);
  double psi = std::min(psi_from_supply(l, transpiration_from_assim(l, assim0, ci0),
                                        psi_upstream),
                        l.psi_crit * (1.0 - 1e-9));

  for (int k = 0; k < newton_steps; ++k) {
    // The bracket, not just positivity, for `solve()`'s reason.
    if (!(psi > psi_upstream) || psi >= l.psi_crit) {
      break;
    }
    const double lambda = l.lambda_TF24_floor(psi);
    const double xi = std::sqrt(Q / lambda);
    const double ci = l.ca_ * xi / (xi + sqrt_D);
    const double assim = l.assim_colimited(ci);
    const double dassim = dassim_dci(l, ci, electron_transport);
    const double u = l.ca_ - ci;
    const double E = transpiration_from_assim(l, assim, ci);
    const double dE_dci = l.H2O_CO2_stom_diff_ratio_ * 1e-3 * l.vpd_leaf_ /
                          kg_to_mol_h2o * (dassim * u + assim) / (u * u);
    const double dci_dxi = l.ca_ * sqrt_D / ((xi + sqrt_D) * (xi + sqrt_D));
    // lambda_o differentiates away, so TF24's derivative is this curve's.
    const double dxi_dpsi = -0.5 * xi / lambda * dlambda_TF24(l, psi);
    const double R = l.transpiration(psi, psi_upstream) - E;
    const double dR = kmax * l.proportion_of_conductivity(psi) -
                      dE_dci * dci_dxi * dxi_dpsi;
    if (dR == 0.0 || !std::isfinite(R) || !std::isfinite(dR)) {
      break;
    }
    const double psi_next = psi - R / dR;
    if (!std::isfinite(psi_next)) {
      break;
    }
    psi = std::max(psi_upstream,
                   std::min(psi_next, l.psi_crit * (1.0 - 1e-9)));
  }

  return evaluate_at_lambda(l, psi, l.lambda_TF24_floor(psi), Q, sqrt_D);
}

// The validity guard (note 2): finite, ci/ca > 0.5, and psi_stem strictly inside
// (psi_soil, psi_crit (1 - 1e-9)). Below ci/ca ~ 0.3 the fixed point stops
// existing, and an answer on either bound is a clamp, not an optimum; both go to
// the exact solve. psi_crit bounds the expansion, not the model. It reads no
// driver, so `true` never means the answer is right (the energy balance, above).
inline bool within_guard(const Leaf &l, const Solution &s) {
  return std::isfinite(s.psi_stem) && std::isfinite(s.ci) &&
         std::isfinite(s.assim) && s.ci / l.ca_ > 0.5 &&
         s.psi_stem > l.supply_psi_soil_scalar() &&
         s.psi_stem < l.psi_crit * (1.0 - 1e-9);
}

// ---------------------------------------------------------------------------
// The entry point: solve, guard, and write the operating point -- or fall back.
// ---------------------------------------------------------------------------
//
// ⚠️ Hazard 8: every output the exact solve writes is written here, by the same
// functions. The closed form only locates psi*; `profit_psi_stem_for<K>` and
// `lambda_for<K>` write the rest, as at the end of `optimise_psi_stem_single`. So
// `ci_`, `assim_colimited_` and `transpiration_` are the exact state at the
// approximate psi*, not the closed form's USO estimates.
template <Leaf::CostCurve K>
inline void optimise_into(Leaf &l) {
  using CostCurve = Leaf::CostCurve;
  if constexpr (Curve<K>::closed_form == ClosedForm::Implemented) {
    // Frozen: these three have per-curve solvers below, and no new curve adds one.
    static_assert(K == CostCurve::TF24 || K == CostCurve::CF77 ||
                      K == CostCurve::TF24_floor,
                  "the closed form's per-curve solvers are TF24's, CF77's and "
                  "TF24_floor's alone");
    l.clear_collar_solve_state();
    l.check_cost_parameters<K>();
    ++l.closed_form_calls_;

    // At lambda_o == 0 TF24_floor is TF24, whose wet-end lambda is zero, so xi
    // diverges. Refused here and not in `check_cost_parameters`: zero is valid for
    // the model and invalid only for this method.
    if constexpr (K == CostCurve::TF24_floor) {
      if (!(l.TF24_floor_lambda_o > 0.0)) {
        util::stop("the closed form needs TF24_floor_lambda_o > 0; got " +
                   util::to_string(l.TF24_floor_lambda_o) + ". At lambda_o = 0 this "
                   "curve IS TF24 -- its wet-end marginal cost of water is zero "
                   "rather than a constant, so the USO slope xi = sqrt(Q/lambda) "
                   "diverges and the leading order this method expands about does "
                   "not exist. Seat TF24 instead, whose own power-law start "
                   "(n = stem_c*TF24_beta2 - 1) is implemented, or use method "
                   "\"exact\", which takes lambda_o = 0 without complaint.");
      }
    }

    const double psi_soil = l.supply_psi_soil_scalar();
    l.opt_psi_stem_ = psi_soil;

    // Soil drier than the stem can reach leaves one feasible potential: the exact
    // solve's branch, verbatim.
    if (psi_soil > l.dry_end<K>()) {
      l.profit_ = l.profit_psi_stem_for<K>(psi_soil, psi_soil);
      l.lambda_emergent_ = l.lambda_for<K>(psi_soil, psi_soil);
      return;
    }

    // A soil at or past psi_crit empties the expansion's window while a curve
    // unbounded there still has an interval, so the exact solve answers.
    if (!(psi_soil < l.psi_crit)) {
      ++l.closed_form_fallbacks_;
      l.last_solve_fell_back_ = true;
      l.optimise_psi_stem_single<K>();
      return;
    }

    Solution s{};
    if constexpr (K == CostCurve::CF77) {
      s = solve_CF77(l, psi_soil);
    } else if constexpr (K == CostCurve::TF24_floor) {
      s = solve_TF24_floor(l, psi_soil);
    } else {
      s = beta2_is_exact(l) ? solve_exact_beta2(l, psi_soil)
                            : solve(l, psi_soil);
    }

    if (!within_guard(l, s)) {
      // The whole exact solve, not a refinement of s: the guard fires where the
      // fixed point may not exist.
      ++l.closed_form_fallbacks_;
      l.last_solve_fell_back_ = true;
      l.optimise_psi_stem_single<K>();
      return;
    }

    // Into the exact solve's bracket, so the clamp is a property of the answer
    // and not of one code path.
    const double psi =
        std::max(psi_soil, std::min(s.psi_stem, l.psi_crit));
    l.opt_psi_stem_ = psi;
    l.profit_ = l.profit_psi_stem_for<K>(psi, psi_soil);
    l.lambda_emergent_ = l.lambda_for<K>(psi, psi_soil);
  } else if constexpr (Curve<K>::closed_form == ClosedForm::Unwritten) {
    // Feasible, not implemented: each needs a `dlambda_*` for the Newton step, and
    // ProfitMax the |A|max rescaling of lambda (`profitmax_A_max_` and
    // `profitmax_k_span_` are already computed).
    util::stop("the closed form is not implemented for the " +
               Leaf::curve_name(static_cast<int>(K)) + " cost curve. It EXISTS "
               "for this curve -- its benefit-link derivative does not depend on "
               "the solution, and its wet-end lambda is a power law (n = 1 for "
               "JS22; 0 for CMax unless CMax_b is zero, then 1; stem_c - 1 for "
               "ProfitMax) -- but the analytic dlambda/dpsi the Newton step needs "
               "has not been written, and ProfitMax additionally needs the |A|max "
               "rescaling. Use method \"exact\", or add the derivative to "
               "closed_form.hpp.");
  } else {
    // ClosedForm::Impossible, not missing: the Log link makes lambda carry A, a
    // fixed point. LeastCost's published closed form (Prentice et al. 2014) solves the
    // first-order condition in ci, where A cancels; this routine inverts lambda
    // into a potential, which the fixed point blocks. The exact solve is checked
    // against that form in `test_leastcost_matches_prentice_closed_form`.
    util::stop("the closed form does not exist for the " +
               Leaf::curve_name(static_cast<int>(K)) + " cost curve, and it is "
               "not a matter of unwritten algebra. SOX, JW26 and LeastCost are the "
               "three curves with the LOG benefit link, where h'(A) = 1/A -- so "
               "lambda carries the assimilation the solve is for, and given "
               "lambda -> xi -> ci -> A the dependency closes into a fixed point "
               "rather than an inversion. Use method \"exact\".");
  }
}

} // namespace closed_form

// The closed-form arm of `Leaf::optimise()`, declared in `leaf.hpp`. Defined
// here because it needs `closed_form`, which needs the complete `Leaf`.
// `<phylloptim.hpp>` includes this header and `leaf_model.hpp` does not.
//
// The refusals are configuration checks, so they run at solve time: route, supply,
// energy balance, thermal damage and g_min are all settable after `set_model`.
inline void Leaf::optimise_closed() {
  if (route_is_collar_) {
    util::stop("the closed form is a STEM-route method: it inverts psi_stem "
               "against a supply measured from psi_soil, and says nothing about "
               "the soil-to-collar path. The leaf is seated on the collar route; "
               "call set_model(\"" +
               curve_name(static_cast<int>(cost_curve_)) +
               "\", \"stem\", \"closed\"), or use method \"exact\".");
  }
  if (!supply_is_single_layer()) {
    util::stop("the closed form is single-layer only: the multi-layer case is a "
               "three-level nest and its closed form does not exist. psi_soil has "
               + util::to_string(static_cast<int>(supply_n_layers())) +
               " layers; use one layer, or method \"exact\".");
  }
  // The solvers' own stop, so the two cannot drift.
  closed_form::require_prescribed_temperature(*this, "the closed form");
  if (use_thermal_damage_) {
    util::stop("the closed form is refused with use_thermal_damage_ on: it has "
               "not been checked against the exact solve with the thermal-damage "
               "layer's shifted capacities, functional fraction and costs in place, "
               "and its validity guard was measured without them. Use method "
               "\"exact\", or set use_thermal_damage_ to false.");
  }
  // The solvers invert the USO form for the stomatal flux, which is the whole
  // flux only at g_min = 0; the residual-conductance forms are not derived.
  if (g_min != 0.0) {
    util::stop("the closed form is refused with a residual conductance: it "
               "assumes the stomatal flux is the whole flux, and with g_min = " +
               util::to_string(g_min) +
               " the cost is read at the total. Use method \"exact\", or set "
               "g_min to 0.");
  }
  with_curve(cost_curve_, [&](auto tag) {
    closed_form::optimise_into<tag.value>(*this);
  });
}

} // namespace phylloptim

#endif
