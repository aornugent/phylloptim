// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_PHOTOSYNTHESIS_HPP_
#define PHYLLOPTIM_LEAF_PHOTOSYNTHESIS_HPP_

// Photosynthesis: the temperature responses of Vcmax, Jmax and R_d, the
// Farquhar kernels on any scalar, ci at a conductance (solve_ci), and the
// operating point's outputs at a stem potential in one pass.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

template <typename T>
inline T Leaf::arrh_curve(double Ea, const T& ref_value, double leaf_temp) const {
  return ref_value*std::exp(Ea*((leaf_temp+C_to_K) - (25 + C_to_K))/((25 + C_to_K)*gas_constant*(leaf_temp+C_to_K)));
}

template <typename T>
inline T Leaf::peak_arrh_curve(double Ea, const T& ref_value, double leaf_temp,
                               double H_d, double d_S) const {
  T arrh = arrh_curve<T>(Ea, ref_value, leaf_temp);
  double arg2 = 1 + std::exp((d_S*(25 + C_to_K) - H_d)/(gas_constant*(25 + C_to_K)));
  double arg3 = 1 + std::exp((d_S*(leaf_temp + C_to_K) - H_d)/(gas_constant*(leaf_temp + C_to_K)));

  return arrh * arg2/arg3;
}

// Every scalar update_temperature_dependent_params() reads: its cache key
// (hazard 10). ⚠️ Keep it beside that function and in step with it: an input
// missing here is served stale with no error; `photo_temp_key_size` makes a
// mismatched list a compile error. Compared by ==, so one ULP invalidates.
inline std::array<double, Leaf::photo_temp_key_size> Leaf::photo_temp_key() const {
  return {leaf_temp_, atm_o2_kpa_,
          vcmax_25, vcmax_ha_, vcmax_H_d_, vcmax_d_S_,
          jmax_25, jmax_ha_, jmax_H_d_, jmax_d_S_,
          gamma_25_, gamma_ha_,
          kc_25_, kc_ha_,
          ko_25_, ko_ha_,
          R_d_25, rd_q10_intercept_, rd_q10_slope_,
          // The thermal cost reaches jmax_; a bool widens to 0.0/1.0.
          double(use_thermal_cost_), T50_, Tcrit_,
          // The thermal-damage gate is here and its inputs (A_, I_, D_) are not:
          // set_physiology bypasses this cache with the gate on; drop that bypass
          // and every input the block reads must join this key. Without the gate,
          // switching it off would be served the gate-on block.
          double(use_thermal_damage_)};
}

// The temperature-dependent photosynthetic parameters at `leaf_temp`, in the
// operation order the golden files pin. electron_transport_ is recomputed from
// the new jmax_ and the call's PPFD_ (set_physiology refreshes it on a cache
// hit). R_d comes from R_d_25 and a declining Q10, so it rises with temperature.
inline void Leaf::update_temperature_dependent_params(double leaf_temp) {
  // First: two kernels called below read the block's temperature from here.
  photo_temp_ = leaf_temp;
  if (use_thermal_damage_) {
    // The ATLS layer (see use_thermal_damage_). One shift moves both capacity
    // optima; the impairment has its own. A itself, for delta_impair's reason.
    const double A = A_;
    const double delta_topt = topt_offset_ + dTopt_max_ * A / (A + K_A_);
    const double vcmax_dS =
        d_S_shifted(vcmax_ha_, vcmax_H_d_, vcmax_d_S_, delta_topt);
    const double jmax_dS = d_S_shifted(jmax_ha_, jmax_H_d_, jmax_d_S_, delta_topt);
    vcmax_ = peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_dS);
    jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_dS);
    // Only the functional area F photosynthesises, so both capacities scale by
    // F; reversible deactivation is already in the peaked Arrhenius. F is held
    // above zero by a softplus floor of width damage_surv_floor, exact in F where
    // F > ~3.4e-5: a min/max floor puts a corner in every output.
    {
      const double F = 1.0 - I_ - D_;
      const double f = damage_surv_floor;
      functional_ = F > 0.0 ? F + f * std::log1p(std::exp(-F / f))
                            : f * std::log1p(std::exp(F / f));
    }
    vcmax_ *= functional_;
    jmax_ *= functional_;
    k_imp_eff_ = k_imp_eff(leaf_temp, delta_impair_at(A));
    k_rec_eff_ = k_rec_eff(leaf_temp);
  } else {
    vcmax_ =
        peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_d_S_);
    jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_d_S_);
    k_imp_eff_ = 0.0;
    k_rec_eff_ = 0.0;
    functional_ = 1.0;
  }
  // Sicangco et al. (2026) Eqn 12: PSII damage scales Jmax by (1 - TC).
  // Branched, not multiplied by a gate-off zero, so gate-off arithmetic is untouched.
  if (use_thermal_cost_) {
    jmax_ *= 1.0 - thermal_cost_at(leaf_temp);
  }
  gamma_ = arrh_curve(gamma_ha_, gamma_25_, leaf_temp);
  ko_ = arrh_curve(ko_ha_, ko_25_, leaf_temp);
  kc_ = arrh_curve(kc_ha_, kc_25_, leaf_temp);
  // Exactly R_d_25 at 25 C, so a check at 25 C alone is blind to this response;
  // test_rd_temperature_response covers it. An unset R_d_25 is a caller error,
  // not a cue to derive one from vcmax.
  if (!std::isfinite(R_d_25) || R_d_25 < 0.0) {
    util::stop("R_d_25 must be a finite, non-negative dark respiration at 25 C; "
               "got " + util::format_double(R_d_25));
  }
  const double q10 =
      rd_q10_intercept_ - rd_q10_slope_ * (leaf_temp + 25.0) / 2.0;
  R_d_ = R_d_25 * std::pow(q10, (leaf_temp - 25.0) / 10.0);
  if (use_thermal_damage_) {
    // The layer's maintenance costs, added to R_d rather than folded into it so
    // gate-off arithmetic is untouched: standing costs for acclimation, repair
    // capacity and protection, and an activity cost for repair done. Protection
    // is signed: the default k_i_ equals k_i_baseline_, so a max(0, .) would put
    // the defaults on a corner.
    const double acclim_load = A_;
    const double protection = k_i_baseline_ - k_i_;
    const double repair_flux = k_rec_eff_ * I_;
    R_d_ += c_acclim_maint_ * acclim_load + c_repair_maint_ * k_rec_ +
            c_protect_maint_ * protection + c_repair_flux_ * repair_flux;
  }
  km_ = (kc_*umol_per_mol_to_Pa_)*(1 + (atm_o2_kpa_*kPa_to_Pa)/(ko_*umol_per_mol_to_Pa_));
  electron_transport_ = electron_transport();
}

// The intercellular CO2 at which supply meets demand,
// `A(ci) umol_to_mol = gc (ca - ci) / P`, root-found on [Gamma* P, hi]. The
// target is negative at Gamma* (A = -R_d) and rises with ci. With A(ca) < 0 the
// leaf is a net source and the root lies above ca, below
// 2 |A(ca)| umol_to_mol P / gc above it, where hi is placed. With no conductance
// the residual is A(ci) = 0: its root where one exists (the g_min -> 0 limit of
// the shut leaf), else ci = +inf.
//
// ci_root_tol (1e-13) sets the floor of what every output means; it is not the
// settable `ci_abs_tol`, which reaches only the Medlyn route. Worst golden-grid
// difference from a converged solve, by TOMS748 tolerance (the plateaus are its
// discrete iterations):
//
//   1e-7: 1.82e-07   1e-8: 3.15e-08   1e-10, 1e-11: 5.41e-10   1e-13: 7.48e-13
//
// The Newton solve below stops at rounding; TOMS748 is its fallback. Why 1e-13
// and not 1e-10, and why not to loosen it: .claude/measurements.md, "Magnitudes:
// what a difference means".
inline double Leaf::solve_ci(double gc) {
  // NA in, NA out, without entering the solver.
  if (!std::isfinite(gc)) {
    return util::na_value;
  }
  auto target = [&](double x) mutable -> double {
    const double assim_colimited_x_ = assim_colimited(x);
    return assim_colimited_x_ * umol_to_mol -
      (gc * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
  };
  const double lo = gamma_ * umol_per_mol_to_Pa_;
  const double A_ca = assim_colimited(ca_);
  double hi = ca_;
  // The target at ca is A(ca) exactly (the supply term is gc * 0), handed over
  // rather than re-evaluated: bit-identical, since x - 0.0 is x.
  if (A_ca < 0.0) {
    if (gc > 0.0) {
      hi = ca_ + 2.0 * (-A_ca) * umol_to_mol * (atm_kpa_ * kPa_to_Pa) / gc;
    } else {
      // No conductance at all: A(ci) = 0, or nothing nets zero.
      if (!(assim_saturated() > 0.0)) {
        return std::numeric_limits<double>::infinity();
      }
      // A rises with ci towards a positive limit, so doubling reaches it.
      do {
        hi *= 2.0;
      } while (assim_colimited(hi) <= 0.0);
    }
  }
  try {
    const double f_hi = hi == ca_ ? A_ca * umol_to_mol -
                                        (gc * (ca_ - hi) / (atm_kpa_ * kPa_to_Pa))
                                  : target(hi);
    const double f_lo = target(lo);
    // Safeguarded Newton first. The residual is increasing and concave in ci, so
    // Newton from below the root climbs to it and from above overshoots at most
    // once; a step leaving the bracket is replaced by the Newton step from its
    // lower end, or a bisection where that slope is not finite. Stopping on a
    // step below ci_root_tol (1 + |ci|) places the root to rounding, not in a
    // band, so there is no state-dependent noise for a caller to difference.
    if (f_lo < 0.0 && f_hi >= 0.0) {
      if (f_hi == 0.0) {
        return hi;
      }
      const PhotoCapacity<double> cap = photo_capacity<double>();
      const double supply_per_pa = gc / (atm_kpa_ * kPa_to_Pa);
      // The value is `target`'s and the slope is the kernel's own tangent
      // (assim_slope_at), so the Newton root and TOMS748's are the same
      // function's to the last bit and the slope cannot drift from it.
      auto residual = [&](double x, double& slope) -> double {
        slope = assim_slope_at<double>(x, cap) * umol_to_mol + supply_per_pa;
        return target(x);
      };
      double a = lo, b = hi;
      double f_a = f_lo, slope_a = std::numeric_limits<double>::quiet_NaN();
      double x = lo - f_lo * (hi - lo) / (f_hi - f_lo);
      for (int i = 0; i < static_cast<int>(ci_niter); ++i) {
        double slope;
        const double f = residual(x, slope);
        if (f == 0.0) {
          return x;
        }
        if (f < 0.0) {
          a = x;
          f_a = f;
          slope_a = slope;
        } else {
          b = x;
        }
        const double newton = x - f / slope;
        // Convergence is tested before the safeguard: at the root `newton` lands
        // on the bracket end x just became, which the safeguard would bisect.
        if (std::abs(newton - x) <= ci_root_tol * (1.0 + std::abs(x))) {
          return newton;
        }
        if (newton > a && newton < b) {
          x = newton;
          continue;
        }
        // Overshot from above: the step from a, since hi can sit orders of
        // magnitude above the root (a net source with almost no conductance).
        // Where a's slope is not finite (Gamma*) a wide bracket splits geometrically.
        if (!std::isfinite(slope_a)) {
          residual(a, slope_a);
        }
        const double from_a = a - f_a / slope_a;
        x = (from_a > a && from_a < b)   ? from_a
            : (a > 0.0 && b > 4.0 * a) ? std::sqrt(a * b)
                                       : 0.5 * (a + b);
      }
    }
    return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, ci_root_tol,
                                ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("ci_solve", "solve_ci failed: " + std::string(e.what()) +
               "; min=" + util::to_string(lo) +
               "; max=" + util::to_string(hi) +
               "; gc=" + util::to_string(gc));
  }
}

// The TOTAL conductance to CO2, mol CO2 m^-2 LA s^-1: the stomata's, from the
// flux the residual leaves them, plus the cuticle's. Reads the deficit as seated
// (set_leaf_vpd), so on the energy-balance path seat the leaf temperature first.
inline double Leaf:: stom_cond_CO2(double psi_stem, double psi_upstream) {
  const double E = transpiration(psi_stem, psi_upstream);
  double E_s = E - residual_flux(psi_upstream, E);
  if (E_s < 0.0) {
    clamps.note(CLAMP_STOMATAL_FLUX);
    E_s = 0.0;
  }
  const double gc =
      atm_kpa_ * E_s * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_ +
      cuticular_cond_CO2();
  if (use_energy_balance_) [[unlikely]] {  // see dprofit_at_collar_psi
    return boundary_layer_split_CO2(gc).first;
  }
  return gc;
}

// Electron transport: the non-rectangular hyperbola in phi_J * APAR, APAR =
// alpha_PAR_ * PPFD_ (intercepted photon flux, umol m^-2 s^-1).
// electron_transport_at is the same expression in the same order.
inline double Leaf::electron_transport() {
  const double I = phi_J * alpha_PAR_ * PPFD_;
  double electron_transport_ = (I + jmax_ - sqrt(pow(I + jmax_, 2) -
  4 * curv_fact_elec_trans * I * jmax_)) / (2 * curv_fact_elec_trans);
    return electron_transport_;
}

// The scalar-generic photosynthesis cores, umol m^-2 s^-1; the golden files pin
// the double instantiation's association. U is independent of T on purpose: a
// double times a tangent is one multiply, a tangent times a tangent three, so
// wrapping the capacity in T measured slower per solve.
template <typename T, typename U>
inline T Leaf::assim_rubisco_limited_kernel(const T& ci, const U& vcmax) const {
  return (vcmax * (ci - gamma_ * umol_per_mol_to_Pa_)) / (ci + km_);
}

template <typename T>
inline T Leaf::assim_rubisco_limited_kernel(T ci) const {
  return assim_rubisco_limited_kernel<T, double>(ci, vcmax_);
}

template <typename T, typename U>
inline T Leaf::assim_electron_limited_kernel(const T& ci,
                                             const U& transport) const {
  return transport / 4 *
  ((ci - gamma_ * umol_per_mol_to_Pa_) / (ci + 2 * gamma_ * umol_per_mol_to_Pa_));
}

template <typename T>
inline T Leaf::assim_electron_limited_kernel(T ci) const {
  return assim_electron_limited_kernel<T, double>(ci, electron_transport_);
}

template <typename T, typename U>
inline T Leaf::assim_colimited_kernel(const T& ci,
                                      const PhotoCapacity<U>& cap) const {
  const T rubisco = assim_rubisco_limited_kernel<T, U>(ci, cap.vcmax);
  const T electron = assim_electron_limited_kernel<T, U>(ci, cap.transport);

  return (rubisco + electron -
          sqrt(pow(rubisco + electron, 2) -
               4 * cap.curvature * rubisco * electron)) /
             (2 * cap.curvature) -
         cap.respiration;
}

// The colimited form with both limitations at their c_i -> infinity values:
// Rubisco's at vcmax, the electron-limited one at J/4. The same association
// as assim_colimited_kernel, so the two agree where both are finite.
template <typename T, typename U>
inline T Leaf::assim_saturated_kernel(const PhotoCapacity<U>& cap) const {
  const T rubisco = T(cap.vcmax);
  const T electron = T(cap.transport / 4);
  return (rubisco + electron -
          sqrt(pow(rubisco + electron, 2) -
               4 * cap.curvature * rubisco * electron)) /
             (2 * cap.curvature) -
         cap.respiration;
}

inline double Leaf::assim_saturated() const {
  return assim_saturated_kernel<double, double>(photo_capacity<double>());
}

template <typename T>
inline T Leaf::assim_colimited_kernel(T ci) const {
  return assim_colimited_kernel<T, double>(ci, photo_capacity<double>());
}

inline double Leaf::assim_rubisco_limited(double ci_) {
  return assim_rubisco_limited_kernel(ci_);
}

//calculate the light-limited assimilation rate, returns umol m^-2 s^-1
inline double Leaf::assim_electron_limited(double ci_) {
  return assim_electron_limited_kernel(ci_);
}

// returns co-limited assimilation umol m^-2 s^-1, NET of dark respiration (the
// trailing `- R_d_`), so gross assimilation is this value + R_d_.
inline double Leaf::assim_colimited(double ci_) {
  return assim_colimited_kernel(ci_);
}

// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
inline double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream) {

  double assim_colimited_x_ = assim_colimited(x);

  double stom_cond_CO2_x_ = stom_cond_CO2(psi_stem, psi_upstream);
  return assim_colimited_x_ * umol_to_mol -
         (stom_cond_CO2_x_ * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
}

// The ci at which A(ci) = gc(ca - ci). The target is smooth and strictly
// monotone in ci (Ar and Ae vanish linearly at gamma*), so a superlinear
// bracketing solver suits it; the bracket and the degenerate case are solve_ci's.
inline double Leaf::psi_stem_to_ci(double psi_stem, double psi_upstream) {
  return ci_ = solve_ci(stom_cond_CO2(psi_stem, psi_upstream));
}

// The operating point at psi_stem in one pass, with no branch on the flow: the
// flux (zero for a reversed pair), the leaf temperature, the residual the cuticle
// loses, and ci on the total conductance. A shut leaf is the same expressions
// with the stomatal term zero; with no conductance at all solve_ci takes the
// limit. Water moves whatever photosynthesis does: no exit on assim_max_ < 0.
// ⚠️ Tleaf_ is written on every path (hazard 8): with the energy balance on, a
// leaf transpiring nothing is the hottest. The temperature block and the deficit
// are seated at the stored Tleaf_, so the three temperatures are one number.
inline void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream) {
  transpiration_ = operating_point_flux(psi_stem, psi_upstream);
  Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(transpiration_) : leaf_temp_;
  if (use_energy_balance_) {
    update_temperature_dependent_params(Tleaf_);
    set_leaf_vpd(Tleaf_);
  }
  residual_flux_ = residual_flux(psi_upstream, transpiration_);
  // Below the wet bound the stomatal flux reads as zero, counted, so a collar
  // that inverts the pair by rounding is the wet bound, not a refusal.
  double E_s = transpiration_ - residual_flux_;
  if (E_s < 0.0) {
    clamps.note(CLAMP_STOMATAL_FLUX);
    E_s = 0.0;
  }
  stom_cond_CO2_ = atm_kpa_ * E_s * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_ +
                   cuticular_cond_CO2();
  if (use_energy_balance_) [[unlikely]] {  // see dprofit_at_collar_psi
    stom_cond_CO2_ = boundary_layer_split_CO2(stom_cond_CO2_).first;
  }
  ci_ = solve_ci(stom_cond_CO2_);
  // With no conductance at all A is zero by definition, not the root's ~1e-15
  // residue, which the product and ratio curves' lambda turn into +-inf and NaN.
  // On the tape this A takes its rows through the conductance's g_min row.
  assim_colimited_ = !std::isfinite(ci_)          ? assim_saturated()
                     : stom_cond_CO2_ == 0.0      ? 0.0
                                                  : assim_colimited(ci_);
}

// The kernels' arguments, read off a caller's parameter pack. At double they
// equal the members bit for bit (passive_pars builds the pack; responses are at
// photo_temp_), which test_leaf asserts on both temperature paths.

// These three responses do not apply the thermal-damage layer, so each refuses
// a leaf with it on rather than disagree with the member it twins.
template <typename T>
inline T Leaf::vcmax_at(const leaf_pars<T>& pars) const {
  refuse_thermal_damage("Leaf::vcmax_at");
  return peak_arrh_curve<T>(vcmax_ha_, pars[par_vcmax_25], photo_temp_,
                            vcmax_H_d_, vcmax_d_S_);
}

template <typename T>
inline T Leaf::jmax_at(const leaf_pars<T>& pars) const {
  refuse_thermal_damage("Leaf::jmax_at");
  T j = peak_arrh_curve<T>(jmax_ha_, pars[par_jmax_25], photo_temp_, jmax_H_d_,
                           jmax_d_S_);
  // The thermal cost reads no trait, so it is a factor, not a term with rows.
  if (use_thermal_cost_) {
    // Not `*=`, whose overload is ambiguous at a nested scalar.
    j = j * T(1.0 - thermal_cost_at(photo_temp_));
  }
  return j;
}

template <typename T>
inline T Leaf::respiration_at(const leaf_pars<T>& pars) const {
  refuse_thermal_damage("Leaf::respiration_at");
  const double q10 =
      rd_q10_intercept_ - rd_q10_slope_ * (photo_temp_ + 25.0) / 2.0;
  return pars[par_R_d_25] *
         T(std::pow(q10, (photo_temp_ - 25.0) / 10.0));
}

template <typename T>
inline T Leaf::electron_transport_at(const leaf_pars<T>& pars) const {
  // alpha_PAR_ is a leaf property, not a slot, so it enters as a constant.
  const T I = pars[par_phi_J] * alpha_PAR_ * pars[par_PPFD];
  const T jm = jmax_at<T>(pars);
  const T curv = pars[par_curv_fact_elec_trans];
  return (I + jm -
          sqrt(pow(I + jm, 2) - 4 * curv * I * jm)) /
         (2 * curv);
}

// gamma_ and km_ are temperature responses of fixed physiology -- no trait
// reaches either -- so they stay members and the capacity carries the four that
// do.
template <typename T>
inline Leaf::PhotoCapacity<T> Leaf::photo_capacity_at(
    const leaf_pars<T>& pars) const {
  return {vcmax_at<T>(pars), electron_transport_at<T>(pars),
          pars[par_curv_fact_colim], respiration_at<T>(pars)};
}

template <typename T>
inline Leaf::PhotoCapacity<T> Leaf::photo_capacity() const {
  return {T(vcmax_), T(electron_transport_), T(curv_fact_colim), T(R_d_)};
}

template <class S, class Kernel, std::same_as<S>... Args>
inline S Leaf::kernel_slope_at(Kernel&& kernel, const char* what, const S& x,
                               const Args&... args) const {
  using odelia::util::to_passive;
  // Both levels are tangents at double, so nothing here reaches a tape. The inner
  // direction is x's (the slope); the outer carries x, then each argument (rows).
  using inner = odelia::ode::tangent_scalar<double>;
  constexpr std::size_t n_dir = 1 + sizeof...(Args);
  using outer = typename xad::fwd<inner, n_dir>::active_type;

  // At double, the inner level only: there are no rows, and the outer walk would
  // cost n_dir directions on the path that places every operating point.
  if constexpr (std::is_same_v<S, double>) {
    inner xd(x);
    odelia::ode::set_direction(xd, 1.0);
    return odelia::ode::derivative_along(kernel(xd, inner(args)...));
  } else {
    std::array<outer, n_dir> lift;
    std::size_t d = 0;
    const auto set = [&](const S& v) {
      lift[d] = outer(to_passive(v));
      xad::value(xad::derivative(lift[d])[d]) = 1.0;
      ++d;
    };
    set(x);
    (set(args), ...);
    odelia::ode::set_direction(xad::value(lift[0]), 1.0);

    const outer y = std::apply(kernel, lift);
    const double slope = odelia::ode::derivative_along(xad::value(y));

    // A braced list is sequenced left to right, so the counter follows the walk;
    // a stack array because the rows cross as a span, which borrows its storage.
    std::size_t r = 0;
    const auto row = [&](const S& v) -> odelia::input_and_derivative<S> {
      return {v, odelia::ode::derivative_along(xad::derivative(y)[r++])};
    };
    const odelia::input_and_derivative<S> rows[n_dir]{row(x), row(args)...};
    S out;
    const odelia::supply_report report =
        odelia::supply_derivatives<S>(slope, rows, out);
    if (!report.whole) {
      util::stop(std::string("Leaf::kernel_slope_at: ") + what + ": " +
                 report.why);
    }
    return out;
  }
}

template <typename T>
inline T Leaf::assim_slope_at(const T& ci, const PhotoCapacity<T>& cap) const {
  return kernel_slope_at<T>(
      [this]<class U>(const U& c, const U& vcmax, const U& transport,
                      const U& curvature, const U& respiration) -> U {
        return assim_colimited_kernel<U>(
            c, PhotoCapacity<U>{vcmax, transport, curvature, respiration});
      },
      "assimilation in intercellular CO2", ci, cap.vcmax, cap.transport,
      cap.curvature, cap.respiration);
}

} // namespace phylloptim

#endif
