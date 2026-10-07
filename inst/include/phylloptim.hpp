// -*-c++-*-
#ifndef PHYLLOPTIM_HPP_
#define PHYLLOPTIM_HPP_

// phylloptim: a header-only leaf gas-exchange and hydraulics model.
//
// One include gets you phylloptim::Leaf, which couples Farquhar-von Caemmerer-Berry
// photosynthesis to an explicit soil-root-stem-leaf hydraulic path and picks the
// operating point by gain-risk profit maximisation. There is nothing to link
// against.
//
// Typical use:
//
//   #include <phylloptim.hpp>
//
//   phylloptim::Leaf l;                 // default traits; vulnerability tables built
//   l.set_physiology(root_network, PPFD, psi_soil, soil_depth,
//                    leaf_specific_conductance_max, atm_vpd, ca,
//                    leaf_temp, atm_o2_kpa, atm_kpa);
//   l.find_root_collar_psi();           // solve; results in l.opt_psi_stem_ etc.
//
// `root_network` is a phylloptim::RootNetwork: the per-layer root hydraulic
// resistances per unit leaf area. If you have root carbon rather than
// resistances, phylloptim::root_network_from_carbon is the architecture model
// that maps one to the other -- but it is a helper you call, not something the
// solve does for you.
//
// To differentiate the operating point, see the derivative surface at the top of
// phylloptim/leaf_model.hpp.

#include <phylloptim/closed_form.hpp>
#include <phylloptim/constants.hpp>
#include <phylloptim/gradient.hpp>
#include <phylloptim/leaf_model.hpp>
#include <phylloptim/optimize.hpp>
#include <phylloptim/quadrature.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/single_potential.hpp>
#include <phylloptim/uniroot.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/vulnerability.hpp>

#endif
