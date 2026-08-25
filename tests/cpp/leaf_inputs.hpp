// Shared probe fixture: the leaf's inputs at one scalar, with one of them seeded.
//
// The three probes that score the tape path against the row layer all need the
// same two things -- the inputs read off a solved leaf, and a way to name which
// one is moving -- so they come from here rather than from three copies that can
// drift apart.
#ifndef PHYLLOPTIM_PROBE_LEAF_INPUTS_HPP_
#define PHYLLOPTIM_PROBE_LEAF_INPUTS_HPP_

#include <phylloptim.hpp>

#include <type_traits>
#include <vector>

namespace fixture {

// A seed is not a thing on a plain double, and these are built at one to make the
// passive inputs collar_condition takes.
template <class T>
void seed(T& x, double direction) {
  if constexpr (!std::is_same_v<T, double>) {
    phylloptim::seed_direction(x, direction);
  } else {
    (void)x;
    (void)direction;
  }
}

// Everything the leaf answers for, read off a solved leaf, with `par` seeded.
// The two _25 traits and dark respiration enter as the temperature-adjusted
// values the kernels take; the Arrhenius is linear in its reference value, so
// their chain is the ratio between the two.
template <class T>
phylloptim::LeafInputs<T> leaf_inputs(const phylloptim::Leaf& l, const T& collar,
                                      int par) {
  namespace grad = phylloptim::gradient;
  phylloptim::LeafInputs<T> in;
  in.profit = phylloptim::ProfitInputs<T>{
      T(l.vcmax_),          T(l.jmax_),
      T(l.a),               T(l.curv_fact_elec_trans),
      T(l.curv_fact_colim), T(l.PPFD_),
      T(l.R_d_),            T(l.leaf_specific_conductance_max_),
      T(l.stem_b),          T(l.stem_c),
      T(l.beta2),           T(l.cost_scale_TF24),
      collar,               T(0.0)};
  const phylloptim::SupplyAt<double> held = l.held_supply();
  for (double v : held.psi_soil) in.supply.psi_soil.push_back(T(v));
  for (double v : held.r_R_H_min) in.supply.r_R_H_min.push_back(T(v));
  for (double v : held.r_R_V_sum) in.supply.r_R_V_sum.push_back(T(v));
  in.supply.root_b = T(held.root_b);
  in.supply.root_c = T(held.root_c);
  in.psi_crit = T(l.psi_crit);
  in.root_psi_crit = T(l.supply_psi_crit());

  switch (par) {
    case grad::par_vcmax_25: seed(in.profit.vcmax, l.vcmax_ / l.vcmax_25); break;
    case grad::par_jmax_25:
      seed(in.profit.transport_jmax, l.jmax_ / l.jmax_25); break;
    case grad::par_a: seed(in.profit.quantum_yield, 1.0); break;
    case grad::par_curv_fact_elec_trans: seed(in.profit.curv_elec, 1.0); break;
    case grad::par_curv_fact_colim: seed(in.profit.curv_colim, 1.0); break;
    case grad::par_R_d_25: seed(in.profit.respiration, l.R_d_ / l.R_d_25); break;
    case grad::par_PPFD: seed(in.profit.ppfd, 1.0); break;
    case grad::par_stem_b: seed(in.profit.stem_b, 1.0); break;
    case grad::par_stem_c: seed(in.profit.stem_c, 1.0); break;
    case grad::par_beta2: seed(in.profit.beta2, 1.0); break;
    case grad::par_cost_scale_TF24: seed(in.profit.cost_scale, 1.0); break;
    case grad::par_kmax: seed(in.profit.kmax, 1.0); break;
    case grad::par_root_b: seed(in.supply.root_b, 1.0); break;
    case grad::par_root_c: seed(in.supply.root_c, 1.0); break;
    case grad::par_psi_crit: seed(in.psi_crit, 1.0); break;
    case grad::par_root_psi_crit: seed(in.root_psi_crit, 1.0); break;
    default: break;
  }
  const int layer = par - grad::par_psi_soil_first;
  if (layer >= 0 && layer < static_cast<int>(in.supply.psi_soil.size())) {
    seed(in.supply.psi_soil[std::size_t(layer)], 1.0);
  }
  return in;
}

}  // namespace fixture

#endif
