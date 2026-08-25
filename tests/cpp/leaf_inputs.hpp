// -*-c++-*-
// Shared probe fixture: the leaf's inputs at one scalar, with one of them moving.
//
// Named rather than indexed. An index is a second vocabulary for the same fact,
// and two vocabularies can disagree while both compile -- which is the whole
// reason the boundary this fixture drives stopped carrying one.
#ifndef PHYLLOPTIM_TESTS_LEAF_INPUTS_HPP_
#define PHYLLOPTIM_TESTS_LEAF_INPUTS_HPP_

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <type_traits>
#include <vector>

namespace fixture {

// Which input is moving, and where a per-layer one is moving.
enum class Input {
  None,
  vcmax_25, jmax_25, quantum_yield, curv_elec, curv_colim, respiration, ppfd,
  kmax, stem_b, stem_c, beta2, cost_scale, psi_crit,
  root_b, root_c, root_psi_crit,
  psi_soil, root_carbon
};

inline const char* name_of(Input which) {
  switch (which) {
    case Input::None: return "none";
    case Input::vcmax_25: return "vcmax_25";
    case Input::jmax_25: return "jmax_25";
    case Input::quantum_yield: return "a";
    case Input::curv_elec: return "curv_elec";
    case Input::curv_colim: return "curv_colim";
    case Input::respiration: return "R_d_25";
    case Input::ppfd: return "PPFD";
    case Input::kmax: return "kmax";
    case Input::stem_b: return "stem_b";
    case Input::stem_c: return "stem_c";
    case Input::beta2: return "beta2";
    case Input::cost_scale: return "cost_scale";
    case Input::psi_crit: return "psi_crit";
    case Input::root_b: return "root_b";
    case Input::root_c: return "root_c";
    case Input::root_psi_crit: return "root_psi_crit";
    case Input::psi_soil: return "psi_soil";
    case Input::root_carbon: return "root_carbon";
  }
  return "unknown";
}

// The soil the caller drove the leaf with, which the carbon -> resistance map
// needs and the leaf does not hold: it takes resistances and knows nothing about
// carbon.
struct Soil {
  std::vector<double> carbon;
  std::vector<double> depth;
};

// A seed is not a thing on a plain double, and these are built at one to make the
// passive inputs an implicit node takes.
template <class T>
void seed(T& x, double direction) {
  if constexpr (!std::is_same_v<T, double>) {
    phylloptim::seed_direction(x, direction);
  } else {
    (void)x;
    (void)direction;
  }
}

// Everything the leaf answers for, read off a solved leaf, with `which` seeded.
//
// The two _25 traits and dark respiration enter as the temperature-adjusted
// values the kernels take; the Arrhenius is linear in its reference value, so
// the chain from the trait is the factor between them.
//
// Root carbon is seeded through the architecture model, exactly as a caller that
// owns that model does it -- the resistances are what the leaf reads.
template <class T>
phylloptim::LeafInputs<T> leaf_inputs(const phylloptim::Leaf& l, const T& collar,
                                      Input which, int layer, const Soil& soil) {
  phylloptim::LeafInputs<T> in;
  in.profit = phylloptim::ProfitInputs<T>{
      T(l.vcmax_),          T(l.jmax_),
      T(l.a),               T(l.curv_fact_elec_trans),
      T(l.curv_fact_colim), T(l.PPFD_),
      T(l.R_d_),            T(l.leaf_specific_conductance_max_),
      T(l.stem_b),          T(l.stem_c),
      T(l.beta2),           T(l.cost_scale_TF24),
      collar,               T(0.0)};
  in.psi_crit = T(l.psi_crit);
  in.root_psi_crit = T(l.supply_psi_crit());
  in.supply.root_b = T(l.roots_.root_b);
  in.supply.root_c = T(l.roots_.root_c);
  for (double v : l.supply_psi_soil()) in.supply.psi_soil.push_back(T(v));

  switch (which) {
    case Input::vcmax_25: seed(in.profit.vcmax, l.vcmax_ / l.vcmax_25); break;
    case Input::jmax_25:
      seed(in.profit.transport_jmax, l.jmax_ / l.jmax_25); break;
    case Input::quantum_yield: seed(in.profit.quantum_yield, 1.0); break;
    case Input::curv_elec: seed(in.profit.curv_elec, 1.0); break;
    case Input::curv_colim: seed(in.profit.curv_colim, 1.0); break;
    case Input::respiration:
      seed(in.profit.respiration, l.R_d_ / l.R_d_25); break;
    case Input::ppfd: seed(in.profit.ppfd, 1.0); break;
    case Input::kmax: seed(in.profit.kmax, 1.0); break;
    case Input::stem_b: seed(in.profit.stem_b, 1.0); break;
    case Input::stem_c: seed(in.profit.stem_c, 1.0); break;
    case Input::beta2: seed(in.profit.beta2, 1.0); break;
    case Input::cost_scale: seed(in.profit.cost_scale, 1.0); break;
    case Input::psi_crit: seed(in.psi_crit, 1.0); break;
    case Input::root_b: seed(in.supply.root_b, 1.0); break;
    case Input::root_c: seed(in.supply.root_c, 1.0); break;
    case Input::root_psi_crit: seed(in.root_psi_crit, 1.0); break;
    case Input::psi_soil:
      if (layer >= 0 && layer < static_cast<int>(in.supply.psi_soil.size())) {
        seed(in.supply.psi_soil[std::size_t(layer)], 1.0);
      }
      break;
    default: break;
  }

  std::vector<T> carbon;
  for (std::size_t i = 0; i < soil.carbon.size(); ++i) {
    carbon.push_back(T(soil.carbon[i]));
  }
  if (which == Input::root_carbon && layer >= 0 &&
      layer < static_cast<int>(carbon.size())) {
    seed(carbon[std::size_t(layer)], 1.0);
  }
  phylloptim::root_resistances_from_carbon<T>(
      carbon, phylloptim::layer_thickness(soil.depth), beta_R_H, beta_R_V,
      in.supply.r_R_H_min, in.supply.r_R_V_sum);
  return in;
}

}  // namespace fixture

#endif
