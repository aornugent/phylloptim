// -*-c++-*-
// The root-architecture model the C++ suite drives the leaf with.
//
// Before #33 these two constants were `Leaf` constructor arguments with defaults,
// and the carbon -> resistance map ran inside `set_physiology`. Now the map is a
// helper the CALLER invokes (`phylloptim::root_network_from_carbon`) and the
// constants belong to whoever owns that model -- in plant, `TF24_Strategy`; here,
// this file.
//
// It exists as a shared header rather than a copy per test program for one
// reason: the golden file's bit-exactness depends on these two numbers. Four
// copies would mean a change to one of them fails `test_golden` while the two
// bench programs quietly go on measuring a different root system.
#ifndef PHYLLOPTIM_TESTS_ROOT_NETWORK_HPP_
#define PHYLLOPTIM_TESTS_ROOT_NETWORK_HPP_

#include <phylloptim.hpp>

#include <vector>

namespace fixture {

// The values that used to be `Leaf`'s beta_R_H / beta_R_V defaults, and are
// still plant's TF24 defaults.
inline constexpr double beta_R_H = 3.4e2; // MPa s (mol C) / (mol H2O)
inline constexpr double beta_R_V = 9.4e3; // MPa (mol C) s / (mol H2O) / m^2

// Carbon -> resistances, exactly as set_physiology used to do it internally.
// The suite drives the leaf from a cumulative depth profile, so it differences
// that profile into per-layer widths here (#626) rather than holding widths of
// its own. Bit-exact for the golden grid, whose layers are exactly 1 m.
inline phylloptim::RootNetwork
root_network(const std::vector<double>& root_carbon_per_leaf_area,
             const std::vector<double>& soil_depth) {
  return phylloptim::root_network_from_carbon(
      root_carbon_per_leaf_area, phylloptim::layer_thicknesses(soil_depth),
      beta_R_H, beta_R_V);
}

// The single-potential path's driver: one series resistance, and no
// vulnerability-weighted horizontal term. The SAME RootNetwork the multi-layer
// path takes -- `r_R_V_sum` already means "series resistance to the collar", so
// this is that field with one layer, which is why one set_physiology argument
// serves both paths.
inline phylloptim::RootNetwork series_resistance(double r) {
  phylloptim::RootNetwork out;
  out.r_R_V_sum.assign(1, r);
  return out;
}

}  // namespace fixture

#endif
