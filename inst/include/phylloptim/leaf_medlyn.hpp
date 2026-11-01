// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_MEDLYN_HPP_
#define PHYLLOPTIM_LEAF_MEDLYN_HPP_

// MEDLYN STOMATAL-CONDUCTANCE MODEL
//
// Standalone, R-callable coupling of the Medlyn (2011) optimal stomatal model
// to colimited photosynthesis, not on the optimality routes. beta_ is a
// soil-moisture factor in [0,1] from theta_, theta_w_ and theta_fc_.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

inline double Leaf::medlyn_model_gs(double assim_colimited_){

  const double beta_ = medlyn_beta();

  // The USO expression as published: a WATER VAPOUR conductance, which `g1` is
  // fitted against. The zero-deficit guard reads the driven `atm_vpd_`, the value
  // the expression reads, not the `atm_vpd` input field.
  double gs_h2o;
  if (atm_vpd_ == 0) {
     gs_h2o = g0;
  } else {
     gs_h2o = g0 + medlyn_H2O_CO2_ratio*(1 + (g1*beta_)/sqrt(atm_vpd_))*(assim_colimited_/(ca_*(1/umol_per_mol_to_Pa_)));
  }

  // ...converted to the CO2 basis every other conductance here is on, so
  // `A == stom_cond_CO2_*(ca_ - ci_)` holds on this route too; `g0` included,
  // since Medlyn's `g0` is a water-basis conductance as well.
  medlyn_model_gs_ = gs_h2o / medlyn_H2O_CO2_ratio;
  return medlyn_model_gs_;
}

// Supply==demand residual for the Medlyn solver in ci (x, Pa), multiplied
// through by (ca_ - x) so it stays finite at ca_; for x < ca_ its root is the
// raw conductance difference's:
//   gs_medlyn*(ca_-x) - gs_coupled*(ca_-x),  where
//   gs_coupled*(ca_-x) = assim * (atm_kpa_*kPa_to_Pa) / 1e6.
// Both terms are on the CO2 basis; scaling one alone would move the root.
inline double Leaf::medlyn_stom_cond_minus_coupled_stom_cond(double x) {
  const double assim_colimited_x_ = assim_colimited(x);
  medlyn_model_gs_ = medlyn_model_gs(assim_colimited_x_);
  return medlyn_model_gs_ * (ca_ - x)
         - assim_colimited_x_ * (atm_kpa_ * kPa_to_Pa) / 1e6;
}

// The ci at which the Medlyn conductance balances the diffusion-implied one. The
// residual is negative at both ends of [gamma*, ca_] with a hump between, so it
// can hold a spurious sub-compensation root: locate its maximum first, then
// root-find on [argmax, ca_] for the high-ci operating point (g0 == 0 included).
inline void Leaf::solve_medlyn_ci_numerical(){
  auto target = [&](double x) -> double {
    return medlyn_stom_cond_minus_coupled_stom_cond(x);
  };
  const double lo = gamma_ * umol_per_mol_to_Pa_;
  const double hi = ca_;

  double neg_peak = 0.0;
  const double ci_peak =
      util::brent_fmin([&](double x) { return -target(x); }, lo, hi, ci_abs_tol,
                       &neg_peak);
  const double residual_peak = -neg_peak;

  if (residual_peak <= 0.0) {
    // Supply never reaches demand: no feasible Medlyn operating point. Report
    // the closest approach (the residual maximum) rather than failing.
    ci_ = ci_peak;
  } else {
    try {
      ci_ = util::uniroot(target, ci_peak, hi, ci_abs_tol, ci_niter);
    } catch (const std::exception& e) {
      util::stop_infeasible("ci_solve", "solve_medlyn_ci_numerical failed: " + std::string(e.what()) +
                 "; ci_peak=" + util::to_string(ci_peak) +
                 "; max=" + util::to_string(hi));
    }
  }
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  set_medlyn_shutdown_if_closed();
  return;
}

// Medlyn's closed form for chi = ci/ca: the numerical route's limit as g0 -> 0,
// dropping g0, Gamma* and respiration. beta multiplies g1 here as in
// `medlyn_model_gs`, so the two routes are one model.
inline void Leaf::solve_medlyn_ci_analytical(){

  const double beta_ = medlyn_beta();
  const double g1_beta = g1 * beta_;
  ci_ = ca_ * (g1_beta/(g1_beta + sqrt(atm_vpd_)));
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  set_medlyn_shutdown_if_closed();
  return;
}

// The linear soil-moisture factor, clamped to [0, 1]: above field capacity it
// would scale g1 past its fitted value, below the wilting point make it negative.
// CF77's soil price clamps it at 1 only, since it divides by it.
inline double Leaf::medlyn_beta() const {
  const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
  return std::min(std::max(beta, 0.0), 1.0);
}

// A negative conductance (the USO `A/ca` term below the compensation point
// outrunning `g0`) is not an operating point, so it becomes a closed state: no
// water, ci at the compensation point, A = -R_d_. That is not the optimality
// routes' shut state (the evaluator at the wet bound), so in light this route
// steps to -R_d_ as it shuts. Clamped, not refused: a hot or dry leaf reaching
// it is a legitimate state.
inline void Leaf::set_medlyn_shutdown_if_closed(){
  if (medlyn_model_gs_ < 0.0 || stom_cond_CO2_ < 0.0) {
    stom_cond_CO2_ = 0.0;
    medlyn_model_gs_ = 0.0;
    transpiration_ = 0.0;
    ci_ = gamma_ * umol_per_mol_to_Pa_;
    assim_colimited_ = -R_d_;
  }
}

inline double Leaf::g1_eff() const {
  const double chi = ci_ / ca_;
  if (!std::isfinite(chi) || chi >= 1.0) {
    return util::na_value;
  }
  return chi * std::sqrt(vpd_leaf_) / (1.0 - chi);
}

} // namespace phylloptim

#endif
