// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_COST_CURVES_HPP_
#define PHYLLOPTIM_LEAF_COST_CURVES_HPP_

// The cost curves: each one's cost, profit and lambda, and the per-curve
// tables the solves dispatch on. A new curve needs an arm in every table
// here (each ends in a static_assert that names it), in curve_name() here,
// in with_curve() in leaf.hpp, and in closed_form::optimise_into.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

// dC/dpsi_stem for one cost curve: the only place a cost curve enters a
// derivative, shared by both routes (the collar route multiplies it by
// dpsi_stem/dpsi and adds cost_deriv_upstream<K>). Every arm is explicit and the
// last asserts, so a new curve is a compile error rather than a plausible number.
template <Leaf::CostCurve K>
inline double Leaf::cost_deriv(double psi_stem, double psi_upstream) {
  using AD = xad::fwd<double>::active_type;
  (void)psi_upstream;
  double C_prime;
  if constexpr (K == CostCurve::TF24) {
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    C_prime = xad::derivative(hydraulic_cost_TF_kernel(ps_ad));
  } else if constexpr (K == CostCurve::CF77) {
    // lambda * E, and dE/dpsi_stem is kmax times the conductivity fraction,
    // which stem_curve_integral_deriv returns.
    C_prime = cf77_price() * leaf_specific_conductance_max_ *
              stem_curve_integral_deriv(psi_stem);
  } else if constexpr (K == CostCurve::JS22) {
    // C = gamma*(psi - psi_up)^2.
    C_prime = 2.0 * JS22_gamma * (psi_stem - psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // TF24's arm plus CF77's, each in its own operand order, so both reductions
    // stay bit-exact (see `hydraulic_cost_TF24_floor`).
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    C_prime = xad::derivative(hydraulic_cost_TF_kernel(ps_ad)) +
              TF24_floor_lambda_o * leaf_specific_conductance_max_ *
                  stem_curve_integral_deriv(psi_stem);
  } else if constexpr (K == CostCurve::SOX) {
    // Under the log link the cost is `-log g`, so dC/dpsi is `-g'/g`; `g` is zero
    // only at psi_crit, which no optimiser returns.
    C_prime = -sox_reduction_deriv(psi_stem) / sox_reduction(psi_stem);
  } else if constexpr (K == CostCurve::JW26) {
    // The same, and here `-g'/g` collapses: with `g = 1 - psi/psi_crit` it is
    // exactly `1/(psi_crit - psi)`.
    C_prime = 1.0 / (psi_crit - psi_stem);
  } else if constexpr (K == CostCurve::LeastCost) {
    // Under the log link the cost is `+log W`, `W = E_mol + beta*Vcmax_mol`
    // (`leastcost_flux`), so dC/dpsi is `W'/W` with `W' = kg_to_mol_h2o * kmax * f`:
    // the opposite sign to SOX's `-g'/g`, as W rises where g falls. Vcmax's slope
    // in Tleaf is taken in `dprofit_energy_balance_term`.
    C_prime = kg_to_mol_h2o * leaf_specific_conductance_max_ *
              stem_curve_integral_deriv(psi_stem) /
              leastcost_flux(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::ProfitMax) {
    // C = [k(psi_s) - k(psi)]/k_span, so dC/dpsi = kmax*|f'|/k_span. The thermal
    // half's slope (energy balance on) is taken in dprofit_energy_balance_term.
    const double f = proportion_of_conductivity(psi_stem);
    const double abs_fprime =
        f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
    C_prime = leaf_specific_conductance_max_ * abs_fprime / profitmax_k_span_;
  } else {
    static_assert(K == CostCurve::CMax, "unhandled CostCurve in cost_deriv");
    // Linear in the absolute potential: the form Anderegg et al. state.
    C_prime = CMax_a * psi_stem + CMax_b;
  }
  return C_prime;
}

// Whether a curve's cost reads the upstream potential as well as psi_stem, which
// matters on the collar route, where upstream is the decision variable. A curve
// that reads it but says false here solves the wrong first-order condition with
// no error; `test_collar_foc_matches_profit` checks every curve. Every arm is
// explicit for the reason `cost_deriv`'s are.
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_reads_upstream() {
  if constexpr (K == CostCurve::CF77 || K == CostCurve::JS22 ||
                K == CostCurve::CMax || K == CostCurve::TF24_floor ||
                K == CostCurve::LeastCost) {
    // LeastCost: `log W` reads upstream through E, the only log-link cost that does.
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::SOX ||
                      K == CostCurve::JW26 || K == CostCurve::ProfitMax,
                  "unhandled CostCurve in cost_reads_upstream");
    return false;
  }
}

// Whether a curve's own form ends at psi_crit. For the others the marginal cost
// diverges as the stem dries, so the dry end is the end of the stem curve
// (vulnerability_psi_end). For these four the bound belongs to the curve:
//
//   SOX, JW26  the gain factor reaches zero at psi_crit and is negative past it
//   CF77       a constant price need not overtake A'(E) before the stem's ceiling
//   LeastCost  sits at psi_crit on most rows, dark and light
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_domain_ends_at_psi_crit() {
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26 ||
                K == CostCurve::CF77 || K == CostCurve::LeastCost) {
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::JS22 ||
                      K == CostCurve::CMax || K == CostCurve::ProfitMax ||
                      K == CostCurve::TF24_floor,
                  "unhandled CostCurve in cost_domain_ends_at_psi_crit");
    return false;
  }
}

// SOX and JW26 maximise A g, g falling to 0 at psi_crit. With A <= 0 at every
// opening the supremum is 0 as g -> 0, driving a dark leaf to psi_crit, so a
// leaf with A <= 0 at the dry end (A rises with opening) is shut at the wet
// bound instead. The argmax jumps at the light compensation point regardless.
template <Leaf::CostCurve K>
constexpr bool Leaf::objective_is_product() {
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26) {
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::CF77 ||
                      K == CostCurve::JS22 || K == CostCurve::CMax ||
                      K == CostCurve::ProfitMax || K == CostCurve::TF24_floor ||
                      K == CostCurve::LeastCost,
                  "unhandled CostCurve in objective_is_product");
    return false;
  }
}

template <Leaf::CostCurve K>
inline bool Leaf::dark_under_product(double psi_dry, double psi_upstream) {
  static_assert(objective_is_product<K>(),
                "dark_under_product is only defined for a product objective");
  set_leaf_states_rates_from_psi_stem(psi_dry, psi_upstream);
  return !(assim_colimited_ > 0.0);
}

template <Leaf::CostCurve K>
inline double Leaf::dry_end() const {
  if constexpr (cost_domain_ends_at_psi_crit<K>()) {
    return psi_crit;
  } else {
    return vulnerability_psi_end(stem_b, stem_c);
  }
}

// dC/dpsi_upstream at fixed psi_stem, for the curves `cost_reads_upstream` names.
// The collar FOC subtracts it with a unit chain factor, since upstream is the
// collar there.
template <Leaf::CostCurve K>
inline double Leaf::cost_deriv_upstream(double psi_stem, double psi_upstream) {
  static_assert(cost_reads_upstream<K>(),
                "cost_deriv_upstream is only defined where the cost reads upstream");
  if constexpr (K == CostCurve::CF77) {
    // lambda * kmax * (G(psi_stem) - G(psi_up)), so -lambda * kmax * G'(psi_up),
    // with G' the table's own slope: the solve ran on the table.
    (void)psi_stem;
    return -(cf77_price() * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream));
  } else if constexpr (K == CostCurve::JS22) {
    return -2.0 * JS22_gamma * (psi_stem - psi_upstream);
  } else if constexpr (K == CostCurve::CMax) {
    // d/du of (psi - u)*(a*(psi + u)/2 + b) is -(a*u + b): the integrand at the
    // lower limit, as for any integral of the price.
    (void)psi_stem;
    return -(CMax_a * psi_upstream + CMax_b);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // TF24's half reads psi_stem alone; the floor's lambda_o*E is CF77's arm.
    (void)psi_stem;
    return -(TF24_floor_lambda_o * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream));
  } else {
    static_assert(K == CostCurve::LeastCost, "unhandled CostCurve");
    // `(dW/du)/W` with `dW/du = -kg_to_mol_h2o * kmax * G'(u)`: CF77's arm with
    // kg-to-mol in place of the price, over `cost_deriv`'s `W`.
    return -(kg_to_mol_h2o * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream)) /
           leastcost_flux(psi_stem, psi_upstream);
  }
}

// Hydraulic cost equations

// --- Marginal cost of water -------------------------------------------------

inline double Leaf::lambda_TF24(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  return TF24_cost_scale * TF24_beta2 * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0) *
         pow(1.0 - f, TF24_beta2 - 1.0) / leaf_specific_conductance_max_;
}

// JS22: dC/dpsi = 2*gamma*(psi - psi_up) over dE/dpsi = kmax*f(psi). Unlike
// TF24's, f does not cancel, so this divides by it: safe inside the feasible
// bracket, where f never reaches zero.
inline double Leaf::lambda_JS22(double psi_stem, double psi_upstream) const {
  const double f = proportion_of_conductivity(psi_stem);
  return 2.0 * JS22_gamma * (psi_stem - psi_upstream) /
         (leaf_specific_conductance_max_ * f);
}

// CMax: the marginal cost does not vanish at the wet end unless CMax_b is zero,
// so it is the one curve that can be boundary-soil with a non-zero cost there.
inline double Leaf::lambda_CMax(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double f = proportion_of_conductivity(psi_stem);
  return (CMax_a * psi_stem + CMax_b) /
         (leaf_specific_conductance_max_ * f);
}

// TF24_floor: TF24's own cost plus a linear price of water, so lambda is
// lambda_TF24 + lambda_o. Called, not repeated: at lambda_o == 0 it is
// lambda_TF24 bit for bit, since that ends in a division no FMA can absorb.
inline double Leaf::lambda_TF24_floor(double psi_stem) const {
  return lambda_TF24(psi_stem) + TF24_floor_lambda_o;
}

// SOX's objective is a product: at its optimum A'g + Ag' = 0, so
//
//     lambda = A'/E' = -A*(g'/g) / (kmax*f)
//
// It carries a factor of A, since a product prices water as a share of
// assimilation (the-models.Rmd), and diverges at psi_crit, where g = 0: never
// pinned there, and not a number to evaluate at the dry end.
inline double Leaf::lambda_SOX(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double g = sox_reduction(psi_stem);
  const double gp = sox_reduction_deriv(psi_stem);
  const double f = proportion_of_conductivity(psi_stem);
  return -assim_colimited_ * (gp / g) / (leaf_specific_conductance_max_ * f);
}

// Jones: g = 1 - psi/psi_crit makes -g'/g = 1/(psi_crit - psi), so lambda =
// A / [(psi_crit - psi) kmax f(psi)], as the-models.Rmd quotes; it carries A and
// diverges at psi_crit, as lambda_SOX's does.
inline double Leaf::lambda_JW26(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double f = proportion_of_conductivity(psi_stem);
  return assim_colimited_ /
         ((psi_crit - psi_stem) * leaf_specific_conductance_max_ * f);
}

// ProfitMax: cost (k_soil - kmax f)/k_span + TC(Tleaf), in carbon units by |A|max:
//
//     lambda = |A|max * [ |f'|/(f*k_span) + (dTC/dT)(dT/dE) ]
//
// The thermal term vanishes with either gate off.
inline double Leaf::lambda_ProfitMax(double psi_stem) {
  const double f = proportion_of_conductivity(psi_stem);
  if (!(f > 0.0) || !(profitmax_k_span_ > 0.0) || !(profitmax_A_max_ > 0.0)) {
    return util::na_value;
  }
  // |f'| = (c/b)(psi/b)^(c-1) f, positive because f is decreasing.
  const double abs_fprime =
      (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0) * f;

  double thermal = 0.0;
  if (use_thermal_cost_ && use_energy_balance_) {
    const double E = operating_point_flux(psi_stem, supply_psi_soil_scalar());
    double dT_dE = 0.0;
    const double Tleaf = leaf_temp_from_E(E, &dT_dE);
    // The same slope dprofit_energy_balance_term subtracts from dJ/dpsi, so the
    // reported lambda and the first-order condition cannot disagree about it.
    thermal = thermal_cost_temp_deriv(Tleaf) * dT_dE;
  }
  return profitmax_A_max_ *
         (abs_fprime / (f * profitmax_k_span_) + thermal);
}

inline double Leaf::marginal_cost_water() const {
  return lambda_TF24(opt_psi_stem_);
}

inline double Leaf::marginal_cost_water_molar() const {
  // umol CO2 (kg H2O)^-1 -> mol CO2 (mol H2O)^-1
  return marginal_cost_water() * umol_to_mol / kg_to_mol_h2o;
}

inline double Leaf::marginal_cost_water_multilayer() {
  // S is the soil->collar conductance, positive by construction (hazard 2).
  const double S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
  if (!std::isfinite(S) || S <= 0.0) {
    return util::na_value;
  }
  // f is the STEM vulnerability curve at the collar suction: kmax*f(psi_r) is
  // dE_stem/dpsi_r, the stem-side conductance where the two paths meet.
  const double f_r = proportion_of_conductivity(opt_root_psi_);
  return lambda_TF24(opt_psi_stem_) *
         (1.0 + leaf_specific_conductance_max_ * f_r / S);
}

// Pure: no write to hydraulic_cost_, so the AD pass cannot scribble model state
// while probing. The caching is the double entry point's job.
template <typename T, typename U>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem, const U& scale,
                                        const U& b, const U& c,
                                        const U& beta2) const {
  return scale *
         pow((1 - proportion_of_conductivity_kernel<T, U>(psi_stem, b, c)), beta2);
}

template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(T psi_stem) const {
  return hydraulic_cost_TF_kernel<T, double>(psi_stem, TF24_cost_scale, stem_b,
                                             stem_c, TF24_beta2);
}

template <typename T>
inline T Leaf::cost_slope_at(const T& sigma, const leaf_pars<T>& pars) const {
  // Named rather than passed inline: a row borrows its input by reference, and a
  // temporary is not a thing to hand one.
  const T b =
      phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], pars[par_stem_c]);
  return kernel_slope_at<T>(
      [this]<class U>(const U& s, const U& scale, const U& bb, const U& c,
                      const U& beta2) -> U {
        return hydraulic_cost_TF_kernel<U, U>(s, scale, bb, c, beta2);
      },
      "hydraulic cost in stem potential", sigma, pars[par_TF24_cost_scale], b,
      pars[par_stem_c], pars[par_TF24_beta2]);
}

// The curve itself is closed form, so the (P50, c) chain comes out of the one
// derivation rather than a supplied row: only the pre-integrated G(psi) is tabulated and
// needs its rows supplied (stem_integral_at).
template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem,
                                        const leaf_pars<T>& pars) const {
  return hydraulic_cost_TF_kernel<T>(
      psi_stem, pars[par_TF24_cost_scale],
      phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], pars[par_stem_c]),
      pars[par_stem_c], pars[par_TF24_beta2]);
}

inline double Leaf::hydraulic_cost_TF(double psi_stem) {

  hydraulic_cost_ = hydraulic_cost_TF_kernel(psi_stem);

return hydraulic_cost_;
}

// Sicangco et al. (2026) Eqns 8-9: the fraction of maximum PSII damage at this
// leaf temperature, in [0,1]; `r` follows gsthermal, not the printed Eqn 8 (see
// use_thermal_cost_). Exactly 0.0 with the gate off, so callers add it freely.
inline double Leaf::thermal_cost_at(double leaf_temp) const {
  if (!use_thermal_cost_) {
    return 0.0;
  }
  const double span = T50_ - Tcrit_;
  if (!(span > 0.0)) {
    util::stop("thermal cost needs T50 > Tcrit; got T50 = " +
               util::format_double(T50_) + ", Tcrit = " +
               util::format_double(Tcrit_));
  }
  return 1.0 / (1.0 + std::exp(-(2.0 / span) * (leaf_temp - T50_)));
}

inline double Leaf::thermal_cost_temp_deriv(double leaf_temp) const {
  if (!use_thermal_cost_) {
    return 0.0;
  }
  const double span = T50_ - Tcrit_;
  const double TC = thermal_cost_at(leaf_temp);
  return (2.0 / span) * TC * (1.0 - TC);
}

// Profit functions

inline double Leaf::profit_psi_stem_TF(double psi_stem, double psi_upstream) {
set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF(psi_stem);

  return benefit_ - cost;
}

// `CF77_lambda_ * E`, written into `hydraulic_cost_` because it IS a carbon flux
// (umol C per kg water times kg m^-2 s^-1). Straight after
// set_leaf_states_rates_from_psi_stem the memo returns that point's E.
inline double Leaf::hydraulic_cost_CF77(double psi_stem,
                                                double psi_upstream) {
  hydraulic_cost_ = cf77_price() * operating_point_flux(psi_stem, psi_upstream);
  return hydraulic_cost_;
}

// The CF77 price with the optional Medlyn soil-moisture factor, rising without
// bound toward the wilting point:
//
//   beta = (theta - theta_w) / (theta_fc - theta_w),   price = lambda / beta
//
// beta is clamped at 1, so CF77_lambda_ is a floor and theta >= theta_fc is plain
// CF77 exactly. The off path returns the member, not a division by 1.0, so it
// stays bit-identical to plain CF77.
inline double Leaf::cf77_price() const {
  if (!CF77_soil_beta_) { return CF77_lambda_; }
  const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
  return CF77_lambda_ / (beta > 1.0 ? 1.0 : beta);
}

inline double Leaf::profit_psi_stem_CF77(double psi_stem,
                                                 double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_CF77(psi_stem, psi_upstream);

  return benefit_ - cost;
}

// gamma * (dpsi)^2, in carbon units like TF24 and CF77, so it goes in
// `hydraulic_cost_`. `dpsi*dpsi`, not `pow`: the optimiser's inner loop (hazard 5).
inline double Leaf::hydraulic_cost_JS22(double psi_stem,
                                        double psi_upstream) {
  const double dpsi = psi_stem - psi_upstream;
  hydraulic_cost_ = JS22_gamma * dpsi * dpsi;
  return hydraulic_cost_;
}

inline double Leaf::profit_psi_stem_JS22(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_JS22(psi_stem, psi_upstream);

  return benefit_ - cost;
}

// The integral of `CMax_a*psi + CMax_b` from psi_upstream, zero at no flow:
//
//     C = (psi - psi_up) * [ CMax_a*(psi + psi_up)/2 + CMax_b ]
//
// Factored: the difference of squares cancels as the potentials converge (off by
// 1.2e-04 at a 1e-12 drop), which a gradient's step reaches near no flow. For the
// same reason JS22, this family's (2 gamma, -2 gamma psi_up) member, keeps its own.
inline double Leaf::hydraulic_cost_CMax(double psi_stem,
                                        double psi_upstream) {
  hydraulic_cost_ = (psi_stem - psi_upstream) *
                    (0.5 * CMax_a * (psi_stem + psi_upstream) + CMax_b);
  return hydraulic_cost_;
}

inline double Leaf::profit_psi_stem_CMax(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_CMax(psi_stem, psi_upstream);

  return benefit_ - cost;
}

// TF24's cost plus CF77's price, both carbon fluxes. Each term is its parent's own
// expression in its operand order, so lambda_o = 0 and TF24_cost_scale = 0 reduce
// bit-exactly, under FMA contraction too (hazards.md, "One solver: the long
// form"). Reassociating breaks that; test_tf24_floor_reduces guards it.
inline double Leaf::hydraulic_cost_TF24_floor(double psi_stem,
                                           double psi_upstream) {
  hydraulic_cost_ = hydraulic_cost_TF_kernel(psi_stem) +
                    TF24_floor_lambda_o * operating_point_flux(psi_stem, psi_upstream);
  return hydraulic_cost_;
}

// See the declaration for why zero on the other curves is a statement about the
// curve and not about the cost.
inline double Leaf::shadow_cost() const {
  return cost_curve_ == CostCurve::TF24_floor
           ? TF24_floor_lambda_o * transpiration_
           : 0.0;
}

inline double Leaf::profit_psi_stem_TF24_floor(double psi_stem,
                                            double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF24_floor(psi_stem, psi_upstream);

  return benefit_ - cost;
}

// --- the product-objective family --------------------------------------------
//
// Eller's reduction factor: the conductivity fraction rescaled to 1 at zero
// tension and exactly 0 at `psi_crit`, through the k_crit_fraction psi_crit and
// Sperry's k_crit also read, so the three cannot disagree.
// Clamped at zero because the no-flow branch evaluates past psi_crit, where g < 0
// and A = -R_d made a shut leaf's product positive. `g < 0`, not `std::max`, so a
// NaN propagates.
inline double Leaf::sox_reduction(double psi_stem) const {
  const double g = (proportion_of_conductivity(psi_stem) - k_crit_fraction) /
                   (1.0 - k_crit_fraction);
  return g < 0.0 ? 0.0 : g;
}

// Its slope, negative: `|f'| = f*(c/b)*(psi/b)^(c-1)`, lambda_TF24's analytic form.
inline double Leaf::sox_reduction_deriv(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  const double abs_fprime =
      f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
  return -abs_fprime / (1.0 - k_crit_fraction);
}

// A product, not in carbon units: the argmax is log A + log g's, but the value
// must not be compared with a TF24 profit.
// ⚠️ hydraulic_cost_ is poisoned (hazard 8): with no subtracted cost, leaving it
// would report whichever curve ran last on a reused Leaf.
inline double Leaf::profit_psi_stem_SOX(double psi_stem,
                                        double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * sox_reduction(psi_stem);
}

// Jones' linear reduction factor on the derived `psi_crit`, clamped at zero, NaN
// propagating, as `sox_reduction` is.
inline double Leaf::jw26_reduction(double psi_stem) const {
  const double g = 1.0 - psi_stem / psi_crit;
  return g < 0.0 ? 0.0 : g;
}

// Constant, the whole difference from SOX, whose slope follows the vulnerability
// curve. Zero once the clamp binds.
inline double Leaf::jw26_reduction_deriv(double psi_stem) const {
  return psi_stem > psi_crit ? 0.0 : -1.0 / psi_crit;
}

inline double Leaf::profit_psi_stem_JW26(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * jw26_reduction(psi_stem);
}

// --- least-cost / cost-per-unit-carbon ---------------------------------------

// `W = E_mol + beta*Vcmax_mol`, mol m^-2 s^-1; see the declaration for why a
// dimensionless beta needs both conversions. It reads `vcmax_`, not `vcmax_25`, so
// with the energy balance on it moves with the solved leaf temperature, as the
// ratio Prentice et al. minimise; cost_deriv leaves that motion out and
// dprofit_energy_balance_term adds it.
inline double Leaf::leastcost_flux(double psi_stem, double psi_upstream) {
  return operating_point_flux(psi_stem, psi_upstream) * kg_to_mol_h2o +
         LeastCost_beta * vcmax_ * umol_to_mol;
}

// A ratio, not a carbon flux; `hydraulic_cost_` is poisoned as for SOX (hazard 8).
inline double Leaf::profit_psi_stem_LeastCost(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ / leastcost_flux(psi_stem, psi_upstream);
}

// `lambda = A*kg_to_mol_h2o/W`, reading the solved `assim_colimited_` as
// `lambda_SOX` does: valid only at the potential last evaluated.
inline double Leaf::lambda_LeastCost(double psi_stem, double psi_upstream) {
  if (!use_energy_balance_) {
    return assim_colimited_ * kg_to_mol_h2o /
           leastcost_flux(psi_stem, psi_upstream);
  }
  // W's capacity term moves with E through Tleaf, adding beta*Vcmax_mol'(T)*dT/dE:
  // the slope dprofit_energy_balance_term takes in the first-order condition.
  double dT_dE = 0.0;
  const double Tleaf =
      leaf_temp_from_E(operating_point_flux(psi_stem, psi_upstream), &dT_dE);
  double dW_dE = kg_to_mol_h2o;
  if (dT_dE != 0.0) {
    dW_dE += LeastCost_beta * umol_to_mol * vcmax_temp_deriv(Tleaf) * dT_dE;
  }
  return assim_colimited_ * dW_dE / leastcost_flux(psi_stem, psi_upstream);
}

// ===========================================================================
// Sperry et al. (2017) ProfitMax, as Sicangco et al. (2026) implement it
// ---------------------------------------------------------------------------
// Sperry maximises `Profit = CG - HC`, both normalised. Times |A|max that is
// `A - lambda*(k(psi_soil)-k(psi))` with the same argmax, but lambda* moves with
// every driver, so a prescribed lambda is not this model:
//
//     CG      = A(psi)/|A|max
//     HC      = [k(psi_soil)-k(psi)] / [k(psi_soil)-kcrit]
//     lambda* = |A|max / [k(psi_soil) - kcrit]
//
// The profit is in Sperry's units, so a figure overlays the paper's. HC runs 0 to
// 1 however wide the vulnerability curve, where the TF24 cost becomes negligible
// as stem_b widens. kcrit = k_crit_fraction * kmax, so psi_crit is P95, as Sperry.
// HC measures loss from psi_soil on one soil potential, and on a layered profile
// from the collar at which uptake is zero (root_zero_E, prepare_collar_solve's
// wet bracket end).
inline double Leaf::profitmax_reference_psi() {
  if (supply_is_single_layer()) {
    return supply_psi_soil_scalar();
  }
  double bound_a = 0.0, bound_b = 0.0;
  if (!prepare_collar_solve<CostCurve::TF24>(bound_a, bound_b)) {
    return util::na_value;
  }
  return bound_a;
}

inline void Leaf::prepare_profitmax_norms() {
  // A non-finite reference is propagated, not refused: this runs on leaves whose
  // drivers are not seated yet (the batch route seeds |A|max first). The
  // optimiser's `!(profitmax_k_span_ > 0)` guard is the one place that decides it.
  const double psi_soil = profitmax_reference_psi();
  profitmax_k_soil_ =
      leaf_specific_conductance_max_ * proportion_of_conductivity(psi_soil);
  profitmax_k_span_ =
      profitmax_k_soil_ -
      leaf_specific_conductance_max_ * proportion_of_conductivity(psi_crit);
}

inline void Leaf::prepare_profitmax_at(double A_max) {
  prepare_profitmax_norms();
  if (!(std::isfinite(A_max) && A_max > 0.0)) {
    util::stop("prepare_profitmax_at needs a finite positive |A|max -- it is the "
               "normaliser held fixed across a perturbation, so an unset one "
               "means the base point was never solved; got " +
               util::to_string(A_max));
  }
  profitmax_A_max_ = A_max;
}

inline void Leaf::prepare_profitmax() {
  const double psi_soil = profitmax_reference_psi();
  prepare_profitmax_norms();

  // |A|max by the one solver, not a scan, whose argmax is piecewise constant in
  // the parameters (hazard 3). One maximisation of A suffices: |A| peaks at a
  // minimum of A only where A < 0, at closure or along a whole stream that cannot
  // fix carbon (below). dA/dpsi is dJ/dpsi plus the cost derivative, for an
  // identity link. Closure is excluded, as gsthermal drops it (`A[E == 0] <- NA`):
  // the search starts one probe inside the wet bound.
  const auto A_of = [&](double psi) {
    set_leaf_states_rates_from_psi_stem(psi, psi_soil);
    return assim_colimited_;
  };
  const auto dA_of = [&](double psi, bool* ok) {
    const double dJ = dprofit_dpsi_stem<CostCurve::TF24>(psi, psi_soil, ok);
    return dJ + cost_deriv<CostCurve::TF24>(psi, psi_soil);
  };
  const double psi_wet = stem_wet_bound(psi_soil);
  const double width = psi_crit - psi_wet;
  const double lo = psi_wet + 1e-9 * width;   // the open end: closure excluded
  double a_hi = -std::numeric_limits<double>::infinity();
  if (width > 0.0) {
    util::maximise_over_closed_interval_foc(A_of, dA_of, lo, psi_crit, 0,
                                            collar_root_tol,
                                            static_cast<size_t>(ci_niter), &a_hi);
  }
  double a_max = (std::isfinite(a_hi) && a_hi > 0.0) ? a_hi : 0.0;
  // The reachable case where max|A| is not max A: A negative along the whole
  // stream, a leaf too hot to cover respiration. Normalising by zero would leave
  // no objective, so maximise -A, and only then.
  if (!(a_max > 0.0) && width > 0.0) {
    double neg = -std::numeric_limits<double>::infinity();
    util::maximise_over_closed_interval_foc(
        [&](double psi) { return -A_of(psi); },
        [&](double psi, bool* ok) { return -dA_of(psi, ok); }, lo, psi_crit, 0,
        collar_root_tol, static_cast<size_t>(ci_niter), &neg);
    if (std::isfinite(neg) && neg > 0.0) a_max = neg;
  }
  profitmax_A_max_ = a_max;
}

// CG - (HC + TC) at one candidate potential, with the normalisers prepare_profitmax
// seeded. Writes carbon_gain_, hydraulic_cost_norm_ and thermal_cost_ so a caller
// plotting the paper's Figure 2 can read the three components off the object.
inline double Leaf::profit_psi_stem_ProfitMax(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  carbon_gain_ = (profitmax_A_max_ > 0.0)
                     ? assim_colimited_ / profitmax_A_max_
                     : 0.0;
  hydraulic_cost_norm_ =
      (profitmax_k_soil_ -
       leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem)) /
      profitmax_k_span_;
  // At the operating point's leaf temperature; with the energy balance off TC is
  // constant along the stream, as in Sicangco's ProfitMaxTC without it.
  thermal_cost_ = thermal_cost_at(use_energy_balance_
                                      ? leaf_temp_from_E(transpiration_)
                                      : leaf_temp_);

  return carbon_gain_ - (hydraulic_cost_norm_ + thermal_cost_);
}

inline std::vector<double> Leaf::profitmax_curve(int n) {
  if (n < 2) {
    util::stop("profitmax_curve needs at least 2 points");
  }
  prepare_profitmax();
  const double psi_soil = supply_psi_soil_scalar();
  const double psi_wet = stem_wet_bound(psi_soil);
  const double step = (psi_crit - psi_wet) / double(n - 1);
  std::vector<double> out(static_cast<std::size_t>(5 * n));
  for (int i = 0; i < n; ++i) {
    const double p = psi_wet + step * double(i);
    const double profit = profit_psi_stem_ProfitMax(p, psi_soil);
    const std::size_t k = static_cast<std::size_t>(i);
    const std::size_t N = static_cast<std::size_t>(n);
    out[k] = p;
    out[N + k] = carbon_gain_;
    out[2 * N + k] = hydraulic_cost_norm_;
    out[3 * N + k] = thermal_cost_;
    out[4 * N + k] = profit;
  }
  return out;
}

inline void Leaf::set_profitmax_degenerate_state(double psi_stem,
                                                 double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  profit_ = 0;
  carbon_gain_ = 0;
  hydraulic_cost_norm_ = 0;
  thermal_cost_ = 0;
  lambda_emergent_ = util::na_value;
}

// --- runtime curve selection --------------------------------------------------
//
// One switch, so the integer-to-curve mapping is written once. `-Werror=switch`
// makes a curve added to the enum and forgotten here a build failure.
inline std::string Leaf::curve_name(int curve) {
  if (curve < 0 || curve >= n_cost_curves) {
    return "unknown";
  }
  switch (static_cast<CostCurve>(curve)) {
    case CostCurve::TF24: return "TF24";
    case CostCurve::CF77: return "CF77";
    case CostCurve::JS22: return "JS22";
    case CostCurve::CMax: return "CMax";
    case CostCurve::SOX:  return "SOX";
    case CostCurve::JW26: return "JW26";
    case CostCurve::ProfitMax: return "ProfitMax";
    case CostCurve::TF24_floor: return "TF24_floor";
    case CostCurve::LeastCost: return "LeastCost";
  }
  return "unknown";   // unreachable past the bounds check above
}

// Which link each curve composes with assimilation. One table, and the last arm
// asserts, so a curve added to the enum without a link is a build failure.
template <Leaf::CostCurve K>
constexpr Leaf::BenefitLink Leaf::benefit_link() {
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26 ||
                K == CostCurve::LeastCost) {
    // SOX and JW26 are products `A*g`, LeastCost a ratio `A/W`: the link is the
    // same `h = log`, and the sign (-log g against +log W) lands in `cost_deriv`.
    return BenefitLink::Log;
  } else if constexpr (K == CostCurve::ProfitMax) {
    return BenefitLink::Scaled;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::CF77 ||
                  K == CostCurve::JS22 || K == CostCurve::CMax ||
                  K == CostCurve::TF24_floor,
                  "unhandled CostCurve in benefit_link");
    return BenefitLink::Identity;
  }
}

// `h'(A)`, all the derivative needs from the link. The Log arm differentiates
// `log A + log g`, which shares `A*g`'s argmax. `log A` needs `A > 0`, and at
// full closure `A = -R_d`: a non-positive `A` returns NaN, which the caller's
// feasibility flag carries outward.
template <Leaf::CostCurve K>
inline double Leaf::benefit_link_deriv(double A) const {
  if constexpr (benefit_link<K>() == BenefitLink::Log) {
    return A > 0.0 ? 1.0 / A : util::na_value;
  } else if constexpr (benefit_link<K>() == BenefitLink::Scaled) {
    return profitmax_A_max_ > 0.0 ? 1.0 / profitmax_A_max_
                                   : util::na_value;
  } else {
    (void)A;
    return 1.0;
  }
}

// --- the three things that vary, each in ONE place ---------------------------
//
// A new single-layer cost curve is a `CostCurve` member plus one arm in each of
// these three. Every arm is explicit and the last asserts, so a curve forgotten
// here is a compile error, not a fall-through to the `else`.

// What must be set before the search. TF24's parameters are traits with
// defaults; SOX and JW26 take none.
template <Leaf::CostCurve K>
inline void Leaf::check_cost_parameters() {
  if constexpr (K == CostCurve::CF77) {
    // An input with an aligned default, so non-finite means a caller cleared it:
    // a NaN objective's argmax is a property of the bracket and looks plausible.
    if (!std::isfinite(CF77_lambda_)) {
      util::stop("the CF77 cost curve needs CF77_lambda_ set: it is the "
                 "PRESCRIBED marginal value of water in umol C (kg H2O)^-1, and "
                 "it was set to " + util::to_string(CF77_lambda_) + ". Its "
                 "default is aligned to TF24 (see align_curves()); "
                 "set_physiology and set_traits never touch it.");
    }
    // Only with the option on, so a caller who never sets it is not refused for
    // soil values the constant never reads. beta <= 0 prices water at infinity
    // or below zero, which pays the leaf to transpire to the wet bound.
    if (CF77_soil_beta_) {
      if (!(theta_fc_ > theta_w_)) {
        util::stop("CF77_soil_beta_ needs theta_fc_ > theta_w_; got theta_fc_ = " +
                   util::to_string(theta_fc_) + ", theta_w_ = " +
                   util::to_string(theta_w_));
      }
      const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
      if (!std::isfinite(beta) || beta <= 0.0) {
        util::stop("CF77_soil_beta_ is on and the soil-moisture factor beta = "
                   "(theta_ - theta_w_)/(theta_fc_ - theta_w_) is " +
                   util::to_string(beta) + ", which prices water at infinity or "
                   "below zero. theta_ must exceed theta_w_; got theta_ = " +
                   util::to_string(theta_) + ", theta_w_ = " +
                   util::to_string(theta_w_));
      }
    }
  } else if constexpr (K == CostCurve::JS22) {
    // A trait with a default, still checked: a NaN maximises a NaN objective.
    if (!std::isfinite(JS22_gamma) || JS22_gamma < 0.0) {
      util::stop("JS22_gamma must be a finite, non-negative hydraulic unit cost "
                 "in umol C m^-2 s^-1 MPa^-2; got " + util::to_string(JS22_gamma));
    }
  } else if constexpr (K == CostCurve::CMax) {
    // Only CMax_a gets the sign check: CMax_b is signed, negative in the paper's
    // own convention (see the member).
    if (!std::isfinite(CMax_a) || CMax_a < 0.0) {
      util::stop("CMax_a must be a finite, non-negative slope for the CMax "
                 "marginal cost in umol C m^-2 s^-1 MPa^-2; got " +
                 util::to_string(CMax_a));
    }
    if (!std::isfinite(CMax_b)) {
      util::stop("CMax_b must be finite; got " + util::to_string(CMax_b));
    }
  } else if constexpr (K == CostCurve::TF24_floor) {
    // No default, unlike CF77_lambda_: at lambda_o = 0 this curve IS TF24, so an
    // unset price would run TF24 under this curve's name. An explicit zero is
    // accepted (the reduction test's case); the hydraulic half is TF24's traits.
    if (!std::isfinite(TF24_floor_lambda_o) || TF24_floor_lambda_o < 0.0) {
      util::stop("the TF24_floor cost curve needs TF24_floor_lambda_o set: it is the "
                 "PRESCRIBED marginal cost of water as transpiration goes to "
                 "zero, in umol C (kg H2O)^-1, non-negative, and NA until you "
                 "assign one -- set_physiology and set_traits never touch it. It "
                 "has no default because 0 would silently make this curve TF24. "
                 "Set `$TF24_floor_lambda_o` on the leaf, or pass "
                 "`TF24_floor_lambda_o =` to leaf_solve() / leaf_batch(); got " +
                 util::to_string(TF24_floor_lambda_o));
    }
  } else if constexpr (K == CostCurve::LeastCost) {
    // An aligned default, refused when cleared for CF77's reason. Strictly
    // positive, unlike TF24_floor_lambda_o: zero leaves `A/E`, which rises without
    // bound as E -> 0 and puts every optimum at the wet bound.
    if (!std::isfinite(LeastCost_beta) || LeastCost_beta <= 0.0) {
      util::stop("the LeastCost cost curve needs LeastCost_beta set: it is "
                 "Prentice et al. (2014)'s beta, the unit cost of maintaining "
                 "carboxylation capacity relative to that of transpiration, "
                 "DIMENSIONLESS (both fluxes are weighed in mol m^-2 s^-1), "
                 "strictly positive; its default is aligned to TF24 (see "
                 "align_curves()), and set_physiology and set_traits never touch "
                 "it. The published value is 146. Set `$LeastCost_beta` on the "
                 "leaf, or pass "
                 "`LeastCost_beta =` to leaf_solve() / leaf_batch(); got " +
                 util::to_string(LeastCost_beta));
    }
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::SOX ||
                  K == CostCurve::JW26 || K == CostCurve::ProfitMax,
                  "unhandled CostCurve");
  }
}

// Which objective. Not all one kind of number: five are carbon, two products
// (`profit_psi_stem_SOX`) and one a ratio (`profit_psi_stem_LeastCost`).
template <Leaf::CostCurve K>
inline double Leaf::profit_psi_stem_for(double psi_stem, double psi_upstream) {
  if constexpr (K == CostCurve::TF24) {
    return profit_psi_stem_TF(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::CF77) {
    return profit_psi_stem_CF77(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::JS22) {
    return profit_psi_stem_JS22(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::CMax) {
    return profit_psi_stem_CMax(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::SOX) {
    return profit_psi_stem_SOX(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::JW26) {
    return profit_psi_stem_JW26(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    return profit_psi_stem_TF24_floor(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::LeastCost) {
    return profit_psi_stem_LeastCost(psi_stem, psi_upstream);
  } else {
    static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
    return profit_psi_stem_ProfitMax(psi_stem, psi_upstream);
  }
}

// Which lambda to report, normalised to two arguments (TF24's reads psi_stem
// alone). CF77's is exact by construction: its cost is lambda*E, so dC/dpsi over
// dE/dpsi is lambda identically.
template <Leaf::CostCurve K>
inline double Leaf::lambda_for(double psi_stem, double psi_upstream) {
  if constexpr (K == CostCurve::TF24) {
    (void)psi_upstream;
    return lambda_TF24(psi_stem);
  } else if constexpr (K == CostCurve::CF77) {
    // Still independent of BOTH potentials -- the soil-beta option keys on theta,
    // a driver, not on where in the bracket the optimiser is.
    (void)psi_stem;
    (void)psi_upstream;
    return cf77_price();
  } else if constexpr (K == CostCurve::JS22) {
    return lambda_JS22(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::CMax) {
    return lambda_CMax(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::SOX) {
    return lambda_SOX(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::JW26) {
    return lambda_JW26(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // The hydraulic half reads psi_stem alone, as TF24's does.
    (void)psi_upstream;
    return lambda_TF24_floor(psi_stem);
  } else if constexpr (K == CostCurve::LeastCost) {
    return lambda_LeastCost(psi_stem, psi_upstream);
  } else {
    static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
    (void)psi_upstream;
    return lambda_ProfitMax(psi_stem);
  }
}

// lambda_emergent on the COLLAR route: dC/dE along the path the solve moves on,
// the collar psi_r free and psi_stem following. Uptake gives dE/dpsi_r = S and the
// stem path dpsi_stem/dpsi_r = (S + kmax f(psi_r)) / (kmax f(psi_stem)), so the
// part of a cost reading psi_stem alone takes the series correction
//
//     lambda = lambda_stem * (1 + kmax f(psi_r) / S),
//
// the part reading only E carries over unchanged, and a cost that also reads the
// collar (JS22's drop, CMax's lower limit) adds its partial over S. Each row must
// equal dA/dE at an interior optimum: test_collar_lambda_emergent_is_dA_dE.
template <Leaf::CostCurve K>
inline double Leaf::lambda_collar_for() {
  if constexpr (K == CostCurve::TF24) {
    return marginal_cost_water_multilayer();
  } else if constexpr (K == CostCurve::CF77) {
    return cf77_price();
  } else if constexpr (K == CostCurve::LeastCost) {
    return lambda_LeastCost(opt_psi_stem_, opt_root_psi_);
  } else {
    const double S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
    if (!std::isfinite(S) || S <= 0.0) {
      return util::na_value;
    }
    const double series =
        1.0 + leaf_specific_conductance_max_ *
                  proportion_of_conductivity(opt_root_psi_) / S;
    if constexpr (K == CostCurve::TF24_floor) {
      return lambda_TF24(opt_psi_stem_) * series + TF24_floor_lambda_o;
    } else if constexpr (K == CostCurve::JS22 || K == CostCurve::CMax) {
      return lambda_for<K>(opt_psi_stem_, opt_root_psi_) * series +
             cost_deriv_upstream<K>(opt_psi_stem_, opt_root_psi_) / S;
    } else if constexpr (K == CostCurve::SOX || K == CostCurve::JW26) {
      return lambda_for<K>(opt_psi_stem_, opt_root_psi_) * series;
    } else {
      static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
      // The conductance term reads psi_stem alone; the thermal term, where the
      // gates put one in, reads the leaf temperature and so E alone.
      const double f = proportion_of_conductivity(opt_psi_stem_);
      if (!(f > 0.0) || !(profitmax_k_span_ > 0.0) || !(profitmax_A_max_ > 0.0)) {
        return util::na_value;
      }
      const double abs_fprime = (stem_c / stem_b) *
                                pow(opt_psi_stem_ / stem_b, stem_c - 1.0) * f;
      double thermal = 0.0;
      if (use_thermal_cost_ && use_energy_balance_) {
        double dT_dE = 0.0;
        const double Tleaf = leaf_temp_from_E(
            operating_point_flux(opt_psi_stem_, opt_root_psi_), &dT_dE);
        thermal = thermal_cost_temp_deriv(Tleaf) * dT_dE;
      }
      return profitmax_A_max_ *
             (abs_fprime / (f * profitmax_k_span_) * series + thermal);
    }
  }
}

} // namespace phylloptim

#endif
