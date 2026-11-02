// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_ENERGY_BALANCE_HPP_
#define PHYLLOPTIM_LEAF_ENERGY_BALANCE_HPP_

// The leaf energy balance, off unless use_energy_balance_: leaf temperature
// from transpiration, the deficit at the leaf, the boundary layer, and the
// temperature terms of the profit's derivative.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

// See the declaration. The step is 1e-3 K, matching dprofit_energy_balance_term's
// dA/dT difference; the full restore makes calling this inside a derivative safe.
template <class F>
inline double Leaf::temp_deriv(double leaf_temp, F read) {
  const double h = 1e-3;
  const double vc0 = vcmax_, jm0 = jmax_, ga0 = gamma_, ko0 = ko_,
               kc0 = kc_, rd0 = R_d_, km0 = km_, J0 = electron_transport_;
  const double kimp0 = k_imp_eff_, krec0 = k_rec_eff_, func0 = functional_;
  const double pt0 = photo_temp_;
  update_temperature_dependent_params(leaf_temp + h);
  const double up = read();
  update_temperature_dependent_params(leaf_temp - h);
  const double dn = read();
  vcmax_ = vc0; jmax_ = jm0; gamma_ = ga0; ko_ = ko0; kc_ = kc0;
  R_d_ = rd0; km_ = km0; electron_transport_ = J0;
  k_imp_eff_ = kimp0; k_rec_eff_ = krec0; functional_ = func0;
  photo_temp_ = pt0;
  return (up - dn) / (2.0 * h);
}

inline double Leaf::vcmax_temp_deriv(double leaf_temp) {
  return temp_deriv(leaf_temp, [&] { return vcmax_; });
}

// The energy-balance correction to dprofit/dpsi; zero with the gate off. Out of
// line on purpose (hazard 5): a branch inside the inlined cost core can change
// FMA contraction and so the gate-off golden file (`nm -C test_golden | grep dprofit`).
//
// With the gate on psi also reaches A through Tleaf, directly and through the
// ci residual, which gains an explicit T. By the implicit-function theorem
//
//   dA/dT (total) = [A_T*gc*inv_atm + A_prime*dgc_dT*(ca-ci)*inv_atm] / g_ci
//
// times tau = dT/dE * dE/dpsi. Checks on it: at dgc_dT = 0 the damping factor
// (gc*inv_atm)/g_ci lies in (0,1); tau < 0, so above the thermal optimum
// (A_T < 0) the optimum moves drier, toward evaporative cooling. dgc_dT is not
// zero: a form putting it under A_T rather than A_prime passes every test at
// dgc_dT = 0 and is wrong. Derived in vignette("the-models"), Appendix B, "The
// energy-balance term".
template <Leaf::CostCurve K>
inline double Leaf::dprofit_energy_balance_term(
    double ci, double gc, double g_ci, double inv_atm, double gc_const,
    double A_prime, double dgc_dpsistem, double dgc_dpsi, double dpsistem_dpsi,
    double dT_dE, double Tleaf, double E_s, double r_star, double drstar_dpsi,
    double dgc_dT_residual, double bl_slope) {
  if (!use_energy_balance_ || dT_dE == 0.0) {
    return 0.0;
  }

  // dA/dTleaf by a central difference over the model's own temperature block,
  // so it cannot drift from the evaluated function. Not a templated block: that
  // changes what shares an inlined body, and moved the gate-off results. ⚠️ h is
  // fixed and absolute: relative to T, the estimator becomes a function of psi
  // and puts a staircase in the argmax (hazard 3).
  const double h = 1e-3;
  const double vcmax0 = vcmax_, jmax0 = jmax_, gamma0 = gamma_, ko0 = ko_,
               kc0 = kc_, Rd0 = R_d_, km0 = km_, J0 = electron_transport_;
  const double kimp0 = k_imp_eff_, krec0 = k_rec_eff_, func0 = functional_;
  const double pt0 = photo_temp_;
  // Tleaf, not leaf_temp_, which the energy-balance path reads as air temperature.
  const double T0 = Tleaf;

  update_temperature_dependent_params(T0 + h);
  const double A_up = assim_colimited(ci);
  const double vcmax_up = vcmax_;
  update_temperature_dependent_params(T0 - h);
  const double A_dn = assim_colimited(ci);
  const double vcmax_dn = vcmax_;

  // Restored by assignment: exact, and cheaper than a third update.
  vcmax_ = vcmax0; jmax_ = jmax0; gamma_ = gamma0; ko_ = ko0; kc_ = kc0;
  R_d_ = Rd0; km_ = km0; electron_transport_ = J0;
  k_imp_eff_ = kimp0; k_rec_eff_ = krec0; functional_ = func0;
  photo_temp_ = pt0;

  const double A_T = (A_up - A_dn) / (2.0 * h);

  // dE/dpsi from the gc partials: dividing out gc_const and the boundary layer's
  // slope leaves the stomatal flux's slope, and the residual's is added back,
  // since Tleaf reads the total flux.
  const double dE_dpsi =
      (dgc_dpsistem * dpsistem_dpsi + dgc_dpsi) / (gc_const * bl_slope) +
      drstar_dpsi;
  const double tau = dT_dE * dE_dpsi;

  // dgc/dTleaf at fixed E. Stomatal part: -gc_s * esat'(Tleaf)/D, negative (a
  // hotter leaf has a larger deficit) and zero on the deficit's floor; the
  // residual's part is appended. Both pass through the boundary layer's split,
  // which reads no leaf temperature.
  const double gc_s = gc_const * E_s;
  const double dgc_dT =
      (((vpd_leaf_ > vpd_leaf_min)
            ? -gc_s * saturation_vapour_pressure_slope(Tleaf) / vpd_leaf_
            : 0.0) +
       dgc_dT_residual) *
      bl_slope;

  // dA/dTleaf total times tau: dA/dpsi along the temperature path, benefit-side
  // only. The costs that read Tleaf (ProfitMax's thermal cost, LeastCost's
  // capacity term) are taken below. The identity arm keeps the additive curves'
  // operation order.
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    (void)vcmax_up;
    (void)vcmax_dn;
    return tau *
           (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) / g_ci;
  } else {
    // Under a link, J = h(A(psi)) - C(psi) with A the total dependence, so this
    // path takes the call site's single h'(A); h'(A) is never differentiated.
    const double dA_dpsi_T =
        tau * (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) /
        g_ci;
    // The block is restored, so this is the A the call site's `base` used.
    const double A = assim_colimited_kernel(ci);
    double term = benefit_link_deriv<K>(A) * dA_dpsi_T;
    if constexpr (K == CostCurve::ProfitMax) {
      // dTC/dpsi = TC'(Tleaf) * tau; cost-side, so no h' factor, and the slope
      // lambda_ProfitMax reports.
      term -= thermal_cost_temp_deriv(Tleaf) * tau;
    } else if constexpr (K == CostCurve::LeastCost) {
      // C = log W with `vcmax_` at Tleaf, so d(log W)/dpsi gains
      // beta*umol_to_mol*dVcmax/dT*tau/W (`cost_deriv` holds only the E half).
      // Cost-side, so no h'. W reads the total flux, stomatal part from gc.
      const double dVcmax_dT = (vcmax_up - vcmax_dn) / (2.0 * h);
      const double W = gc_s / gc_const * kg_to_mol_h2o +
                       LeastCost_beta * vcmax_ * umol_to_mol +
                       r_star * kg_to_mol_h2o;
      term -= LeastCost_beta * umol_to_mol * dVcmax_dT * tau / W;
    } else {
      (void)vcmax_up;
      (void)vcmax_dn;
    }
    return term;
  }
}

// Saturation vapour pressure es(T) in kPa (Tetens), T in deg C.
inline double Leaf::saturation_vapour_pressure(double temp) const {
  return 0.6108 * exp(17.27 * temp / (temp + 237.3));
}

// Slope of the saturation vapour pressure curve Delta(T) in kPa K^-1, T in deg C:
// the exact derivative of the Tetens form above. Not FAO-56's rounded 4098
// (4.2e-5 low), which reaches the energy-balance first-order condition and put
// dprofit/dpsi_stem 2e-6 to 2e-5 off a central difference of profit.
inline double Leaf::saturation_vapour_pressure_slope(double temp) const {
  return (17.27 * 237.3) * saturation_vapour_pressure(temp) /
         ((temp + 237.3) * (temp + 237.3));
}

// The deficit driving diffusion out of the leaf, from the cavity saturated at
// Tleaf to the air (the driver `atm_vpd_` is at Tair):
//
//   vpd_leaf = esat(Tleaf) - e_air = atm_vpd + [esat(Tleaf) - esat(Tair)]
//
// In that association, so off the energy-balance path the bracket is exactly
// +0.0 and the result bit-identical. Held at `vpd_leaf_min`, a kept exception to
// the no-clamp rule: below it the leaf condenses, and a conductance formed from
// the deficit would be infinite or negative.
inline double Leaf::leaf_vpd_at(double leaf_temp) const {
  if (!use_energy_balance_) {
    return atm_vpd_;
  }
  const double d =
      atm_vpd_ + (saturation_vapour_pressure(leaf_temp) -
                  saturation_vapour_pressure(Tair_));
  return std::max(d, vpd_leaf_min);
}

// Also seats the residual conductance at the same temperature: it is read
// beside the deficit wherever a conductance is formed, and the two must
// belong to one leaf temperature.
inline void Leaf::set_leaf_vpd(double leaf_temp) {
  vpd_leaf_ = leaf_vpd_at(leaf_temp);
  g_min_leaf_ = g_min_at(leaf_temp);
}

// Explicit leaf energy balance: Tleaf = Tair + (Rn - lambda*E) / (rho*cp/ra + g_rad),
// with E the hydraulically pinned transpiration (kg H2O m^-2 s^-1). The one place
// Tleaf and dTleaf/dE are computed: a second copy is how a radiative term gets
// dropped from a derivative while the value keeps it.
inline double Leaf::leaf_temp_from_E(double E, double* dT_dE) const {
  // Sensible heat and the linearised longwave as one summed conductance, so the
  // radiative term cannot be dropped by editing one factor.
  const double g_total = vol_heat_cap_air / ra_ + g_rad_;
  const double Tleaf = Tair_ + (Rn_ - latent_heat_vap * E) / g_total;
  if (dT_dE != nullptr) {
    // Linear in E, so the slope is constant and negative.
    *dT_dE = -latent_heat_vap / g_total;
  }
  return Tleaf;
}

// Refuses, never clamps, a leaf temperature range the model cannot evaluate:
// R_d's Q10, intercept - slope*(T + 25)/2, must be positive for its power to be
// real, and the Tetens form has a pole at -237.3 C. Q10 is linear in T, so
// testing both ends tests the range.
inline void Leaf::check_leaf_temperature_range(double coldest,
                                               double hottest) const {
  const auto q10_at = [&](double t) {
    return rd_q10_intercept_ - rd_q10_slope_ * (t + 25.0) / 2.0;
  };
  if (!(coldest > -237.3) || !(q10_at(coldest) > 0.0) ||
      !(q10_at(hottest) > 0.0)) {
    util::stop("set_physiology: the leaf temperature can run from " +
               util::format_double(coldest) + " to " +
               util::format_double(hottest) +
               " C over the fluxes this solve can reach, and the model is "
               "evaluable only above the Tetens pole (-237.3 C) and where R_d's "
               "Q10, rd_q10_intercept_ - rd_q10_slope_*(T + 25)/2, is positive "
               "(Q10 " + util::format_double(q10_at(coldest)) + " at the cold "
               "end, " + util::format_double(q10_at(hottest)) + " at the hot)");
  }
}

inline double Leaf::boundary_layer_cond_H2O() const {
  if (!use_energy_balance_) {
    return std::numeric_limits<double>::infinity();
  }
  // Molar density of air at air temperature, mol m^-3, as plantecophys's CMOLAR.
  const double air_molar_density =
      atm_kpa_ * kPa_to_Pa / (gas_constant * (Tair_ + C_to_K));
  return bl_vapour_heat_ratio * air_molar_density / ra_;
}

inline double Leaf::boundary_layer_co2_excess() const {
  const double r_s = H2O_CO2_stom_diff_ratio_;
  const double r_b = std::cbrt(r_s * r_s);
  return (r_s - r_b) / boundary_layer_cond_H2O();
}

template <class T>
inline std::pair<T, T> Leaf::boundary_layer_split_CO2(const T& gc_hat) const {
  using std::exp;
  using std::log1p;
  using odelia::util::to_passive;
  const double k = boundary_layer_co2_excess();
  const double s = 1.0 / bl_co2_split_floor_width;
  const T d = T(1.0) - T(k) * gc_hat;
  // The soft maximum of d and the floor, d + softplus(d0 - d) at width w, and
  // its slope in d, the logistic of (d - d0)/w. Written so neither overflows:
  // e = exp(-|z|) is at most 1 whichever side of the floor d is.
  const T z = T(s) * (T(bl_co2_split_floor) - d);
  const bool above = to_passive(z) <= 0.0;
  const T e = above ? T(exp(z)) : T(exp(-z));
  const T d_s = above ? T(d + log1p(e) / T(s))
                      : T(T(bl_co2_split_floor) + log1p(e) / T(s));
  const T dds_dd = above ? T(T(1.0) / (T(1.0) + e)) : T(e / (T(1.0) + e));
  const T gc = gc_hat / d_s;
  const T slope = T(1.0) / d_s + gc_hat * T(k) * dds_dd / (d_s * d_s);
  return {gc, slope};
}

template <class T>
inline T Leaf::boundary_layer_split_slope_at(const T& gc) const {
  using odelia::util::to_passive;
  const double k = boundary_layer_co2_excess();
  const T phi = T(1.0) + T(k) * gc;
  // d = 1/phi above the floor; within 30 widths of it the floor has moved the
  // value and this form would not be its slope.
  if (!(1.0 / to_passive(phi) >
        bl_co2_split_floor + 30.0 * bl_co2_split_floor_width)) {
    util::stop("the boundary layer's split of the CO2 path is on its soft floor "
               "at this operating point (g_c = " + util::to_string(to_passive(gc)) +
               "): the transpiration is more than the boundary layer can carry, "
               "so the point is not a physical one to differentiate.");
  }
  return phi * phi;
}

inline double Leaf::stom_cond_stomatal_H2O() const {
  const double g_total = stom_cond_H2O();
  if (!use_energy_balance_) {
    return g_total;
  }
  // A shut leaf has no conductance on either side of the series, and a
  // non-finite one has nothing to subtract from.
  if (!(g_total > 0.0) || !std::isfinite(g_total)) {
    return g_total;
  }
  const double g_b = boundary_layer_cond_H2O();
  if (!(g_total < g_b)) {
    util::stop("no stomatal conductance exists at this operating point: the total "
               "leaf-to-air conductance " + util::to_string(g_total) +
               " mol m^-2 s^-1 is not below the boundary layer's own " +
               util::to_string(g_b) + ", so the transpiration is more than the "
               "boundary layer can carry at this deficit (ra_ = " +
               util::to_string(ra_) + " s m^-1). Check wind_speed_ and d_.");
  }
  return 1.0 / (1.0 / g_total - 1.0 / g_b);
}

// The surface carries the boundary layer's split but not the leaf-temperature
// terms the solve's first-order condition holds (dprofit_energy_balance_term), so
// with the gate on its rows would differentiate a different model.
inline void Leaf::refuse_energy_balance(const char* who) const {
  if (use_energy_balance_) {
    util::stop(std::string(who) +
               ": refused with use_energy_balance_ on. The derivative surface "
               "holds leaf temperature at the value the solve placed and carries "
               "none of its terms, while the solve's first-order condition "
               "includes them, so its rows would describe a different model from "
               "the one the solve ran. Set use_energy_balance_ to false.");
  }
}

} // namespace phylloptim

#endif
