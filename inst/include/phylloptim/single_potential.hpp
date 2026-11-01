// -*-c++-*-
#ifndef PHYLLOPTIM_SINGLE_POTENTIAL_HPP_
#define PHYLLOPTIM_SINGLE_POTENTIAL_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/util.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace phylloptim {

// The other supply path: one soil water potential and a constant series
// resistance to the root collar, the setting plantecophys and tealeaves users have
// and that Medlyn, least-cost and Cowan-Farquhar are formulated in. It provides
// every operation `Leaf` dispatches on (hazard 6), every potential a positive
// magnitude in MPa. It has no kinks, so `duptake_dpsi` never returns the NaN that
// asks for a finite-difference fallback.
class SinglePotential {
public:
  // Soil water potential as a POSITIVE magnitude, -MPa.
  double psi_soil_ = util::na_value;
  // Series soil-to-collar resistance per unit leaf area (hazard 4),
  // MPa s m2 leaf (mol H2O)^-1. A per-call driver from set_physiology's
  // RootNetwork, as the multi-layer path's resistances are.
  double resistance_ = util::na_value;
  // The smallest series resistance the collar route resolves: the flux is
  // (collar - soil)/r, so a collar placed to collar_root_tol (1e-12) errs in E by
  // about 1e-12/(E*r) relative, worst where E is smallest. Five orders below the
  // nominal 1e3, so no physical root path reaches it.
  static constexpr double resistance_min = 1e-2;
  // Gravitational head to lift water to the collar, MPa. Configuration, zero by
  // default: this path has no depth profile to derive one from (hazard 6).
  double grav_head_ = 0.0;

  // psi_soil_ as a one-element vector, so Leaf's vector-taking soil machinery
  // serves this path unchanged. begin_solve is the only writer.
  std::vector<double> psi_soil_vec_{0.0};

  void clear() {
    psi_soil_ = util::na_value;
    psi_soil_vec_.assign(1, util::na_value);
    // A driver, so cleared too: set_physiology always re-supplies it.
    resistance_ = util::na_value;
  }

  // One layer, always: the consumption buffer gets exactly one entry.
  int n_layers() const { return 1; }

  void set_soil_state(double psi_soil) { psi_soil_ = psi_soil; }

  // The per-call resistance, from the RootNetwork the multi-layer path also takes:
  // `r_R_V_sum[0]` is the series resistance of one layer. A non-zero `r_R_H_min`
  // means a network built for the other path, refused rather than silently dropped.
  void set_supply_resistances(const RootNetwork& network) {
    if (network.r_R_V_sum.size() != 1) {
      util::stop("the single-potential path needs exactly one series resistance "
                 "in r_R_V_sum; got " +
                 std::to_string(network.r_R_V_sum.size()) + " entries");
    }
    for (double h : network.r_R_H_min) {
      if (h != 0.0) {
        util::stop("the single-potential path has no vulnerability-weighted "
                   "horizontal term, so r_R_H_min must be empty or zero; got " +
                   util::to_string(h) + ". Build the leaf with "
                   "set_supply_multilayer() for a root network.");
      }
    }
    const double r = network.r_R_V_sum[0];
    if (!std::isfinite(r) || !(r > 0.0)) {
      util::stop("the single-potential path needs a positive, finite series "
                 "resistance per unit leaf area; got " + util::to_string(r));
    }
    if (r < resistance_min) {
      util::stop("the single-potential path needs a series resistance of at "
                 "least " + util::to_string(resistance_min) +
                 " MPa s m2 (mol H2O)-1, below which the collar route cannot "
                 "resolve the flux; got " + util::to_string(r) +
                 ". The stem routes are the zero-resistance model.");
    }
    resistance_ = r;
  }

  // Per-solve entry point, as MultiLayerRoots::begin_solve: publish the soil
  // state as a vector and return the wettest (here the only) suction.
  double begin_solve() {
    psi_soil_vec_.assign(1, psi_soil_);
    return psi_soil_;
  }

  // E_up at a collar suction, by Ohm's law:
  //
  //   E = (T_collar - T_soil - grav) / r
  //
  // An unset or non-positive resistance is refused rather than returning an
  // infinity.
  void uptake(double T_collar, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake_from(T_collar, psi_soil_, soil_consumption, E_up);
  }

  void uptake_from(double T_collar, double T_soil,
                   std::vector<double>& soil_consumption, double& E_up) const {
    if (!std::isfinite(T_collar)) {
      util::stop_infeasible("uptake", "SinglePotential::uptake invalid input; T_collar=" +
                 util::to_string(T_collar));
    }
    if (!(resistance_ > 0.0)) {
      util::stop("SinglePotential::uptake needs a positive resistance_; got " +
                 util::to_string(resistance_));
    }
    const double E_i = (T_collar - T_soil - grav_head_) / resistance_;
    if (!soil_consumption.empty()) {
      soil_consumption[0] = E_i;  // mol, as MultiLayerRoots leaves it
    }
    E_up = E_i * kg_per_mol_h2o;
    if (!std::isfinite(E_up)) {
      util::stop_infeasible("uptake",
          "SinglePotential::uptake non-finite E_up; T_collar=" +
                 util::to_string(T_collar));
    }
  }

  // Vector-taking forms, so Leaf calls both paths identically. The suction is
  // element 0, not psi_soil_, so Leaf::E_from_Soil_to_Root_Collar, which may be
  // handed any vector, means the same here as for MultiLayerRoots.
  void uptake_at(double T_collar, const std::vector<double>& psi_soil,
                 std::vector<double>& soil_consumption, double& E_up) const {
    if (psi_soil.empty()) {
      util::stop("SinglePotential::uptake_at needs at least one potential");
    }
    uptake_from(T_collar, psi_soil[0], soil_consumption, E_up);
  }

  double duptake_dpsi(double T_collar,
                      const std::vector<double>& psi_soil) const {
    static_cast<void>(T_collar);
    static_cast<void>(psi_soil);
    return duptake_dpsi();
  }

  // Constant, the analytic derivative of uptake: a positive conductance, as
  // MultiLayerRoots::duptake_dpsi is.
  double duptake_dpsi() const { return kg_per_mol_h2o / resistance_; }
};

}  // namespace phylloptim

#endif
