// -*-c++-*-
#ifndef PHYLLOPTIM_HPP_
#define PHYLLOPTIM_HPP_

// phylloptim: a header-only leaf gas-exchange and hydraulics model. One include
// gives phylloptim::Leaf, which couples Farquhar-von Caemmerer-Berry photosynthesis
// to a soil-root-stem-leaf hydraulic path and picks the operating point by
// gain-risk profit maximisation. There is nothing to link against.
//
//   #include <phylloptim.hpp>
//
//   phylloptim::Leaf l;                 // default traits; vulnerability tables built
//   l.set_physiology(root_network, PPFD, psi_soil, soil_depth,
//                    leaf_specific_conductance_max, atm_vpd, ca,
//                    leaf_temp, atm_o2_kpa, atm_kpa);
//   l.find_root_collar_psi();           // solve; results in l.opt_psi_stem_ etc.
//
// `root_network` is a phylloptim::RootNetwork of per-layer root resistances per
// unit leaf area (hazard 4). phylloptim::root_network_from_carbon builds one from
// root carbon; the solve never calls it. To differentiate the operating point,
// see the derivative surface in phylloptim/leaf_derivative_surface.hpp.
//
// closed_form.hpp defines Leaf::optimise_closed, which leaf.hpp declares, so
// it must stay in this list or a call to optimise() fails to link.

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
