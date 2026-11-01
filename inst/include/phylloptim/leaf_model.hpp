// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_MODEL_HPP_
#define PHYLLOPTIM_LEAF_MODEL_HPP_

// `Leaf`, whole: the class (leaf.hpp) and its members' definitions, one header
// per subsystem, in an order their dependencies allow. The header to include
// for the leaf, except that Leaf::optimise_closed is defined in closed_form.hpp,
// which <phylloptim.hpp> adds.

#include <phylloptim/leaf.hpp>
#include <phylloptim/leaf_photosynthesis.hpp>
#include <phylloptim/leaf_hydraulics.hpp>
#include <phylloptim/leaf_energy_balance.hpp>
#include <phylloptim/leaf_damage.hpp>
#include <phylloptim/leaf_cost_curves.hpp>
#include <phylloptim/leaf_collar_solve.hpp>
#include <phylloptim/leaf_stem_routes.hpp>
#include <phylloptim/leaf_derivative_surface.hpp>
#include <phylloptim/leaf_medlyn.hpp>

#endif
