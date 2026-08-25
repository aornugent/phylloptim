// PROBE (not a test, not shipped): Leaf::profit_at against the row layer it is
// meant to replace.
//
// The row layer reaches profit's response to an input through eight producers,
// each re-deriving its own chain: photo_trait_rows, cost_trait_rows,
// transport_trait_rows, hydraulic_cost_row and the rest. profit_at reaches it
// through the two residuals the model actually solves, seeded in one input, and
// nothing in it knows which input that is.
//
//   T1  kmax * (G(sigma) - G(collar)) - E_up = 0
//   T2  A(ci) * k1 - gc * (ca - ci) * k2   = 0
//       profit = A(ci) - C(sigma)
//
// The collar is held, which at an interior point is the whole answer: profit is
// stationary in it. So this scores interior points, and reports the pinned ones
// separately -- there the collar is a bound that moves with the inputs, and the
// row carries a term this holds at zero.
//
//   make CXX=g++ probe_profit_at && ./probe_profit_at

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace grad = phylloptim::gradient;
using phylloptim::Leaf;
using phylloptim::ProfitInputs;
using phylloptim::tangent;
using phylloptim::seed_direction;
using phylloptim::derivative_along;

namespace {

const double kRd25 = 1.44;
const double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,      157.44,   0.30,
                         0.7,      0.99,     7.5,      kRd25,
                         1.0 * 0.000157 / 5.0, 1e3};

grad::Drivers drivers(double psi_soil, double ppfd, int layers) {
  grad::Drivers d;
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers / 0.05;
  }
  d.root_network = fixture::root_network(root, depth);
  d.PPFD = ppfd;
  d.psi_soil = ps;
  d.soil_depth = depth;
  d.atm_vpd = 2.0;
  d.ca = 40.0;
  d.leaf_temp = 25.0;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

struct Named { int par; const char* name; };

const Named kInputs[] = {
    {grad::par_vcmax_25, "vcmax_25"},
    {grad::par_jmax_25, "jmax_25"},
    {grad::par_a, "a"},
    {grad::par_curv_fact_elec_trans, "curv_elec"},
    {grad::par_curv_fact_colim, "curv_colim"},
    {grad::par_R_d_25, "R_d_25"},
    {grad::par_PPFD, "PPFD"},
    {grad::par_stem_b, "stem_b"},
    {grad::par_stem_c, "stem_c"},
    {grad::par_beta2, "beta2"},
    {grad::par_cost_scale_TF24, "cost_scale"},
    {grad::par_kmax, "kmax"},
};

// The inputs at the leaf's current state, with one of them seeded. The two _25
// traits and dark respiration are carried as the temperature-adjusted value the
// kernels take, so the chain from the trait is the ratio.
ProfitInputs<tangent> inputs_at(const Leaf& l, int par) {
  ProfitInputs<tangent> p{
      tangent(l.vcmax_),        tangent(l.jmax_),
      tangent(l.a),             tangent(l.curv_fact_elec_trans),
      tangent(l.curv_fact_colim), tangent(l.PPFD_),
      tangent(l.R_d_),          tangent(l.leaf_specific_conductance_max_),
      tangent(l.stem_b),        tangent(l.stem_c),
      tangent(l.beta2),         tangent(l.cost_scale_TF24),
      tangent(l.transpiration_)};
  switch (par) {
    case grad::par_vcmax_25: seed_direction(p.vcmax, l.vcmax_ / l.vcmax_25); break;
    case grad::par_jmax_25:
      seed_direction(p.transport_jmax, l.jmax_ / l.jmax_25); break;
    case grad::par_a: seed_direction(p.quantum_yield, 1.0); break;
    case grad::par_curv_fact_elec_trans: seed_direction(p.curv_elec, 1.0); break;
    case grad::par_curv_fact_colim: seed_direction(p.curv_colim, 1.0); break;
    case grad::par_R_d_25: seed_direction(p.respiration, l.R_d_ / l.R_d_25); break;
    case grad::par_PPFD: seed_direction(p.ppfd, 1.0); break;
    case grad::par_stem_b: seed_direction(p.stem_b, 1.0); break;
    case grad::par_stem_c: seed_direction(p.stem_c, 1.0); break;
    case grad::par_beta2: seed_direction(p.beta2, 1.0); break;
    case grad::par_cost_scale_TF24: seed_direction(p.cost_scale, 1.0); break;
    case grad::par_kmax: seed_direction(p.kmax, 1.0); break;
    default: break;
  }
  return p;
}

}  // namespace

int main() {
  std::printf("Leaf::profit_at against the row layer, one input at a time\n\n");
  std::printf("   %-8s %-7s %s  %-10s %2s %-12s %14s %14s %11s\n", "psi_soil",
              "PPFD", "L", "point", "n", "worst input", "rows", "profit_at",
              "worst rel");

  double worst_interior = 0.0, worst_pinned = 0.0;
  int compared = 0, states = 0;

  for (double psi_soil : {0.2, 0.8, 1.5, 2.5, 3.5, 4.5}) {
    for (double ppfd : {30.0, 300.0, 1000.0, 2000.0}) {
      for (int layers : {1, 3}) {
        const grad::Drivers d = drivers(psi_soil, ppfd, layers);
        Leaf l;
        const grad::Settings s;
        grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
        try {
          l.find_root_collar_psi();
        } catch (const std::exception&) {
          continue;
        }
        ++states;
        const bool interior =
            l.operating_point_kind() == Leaf::OperatingPointKind::Interior;

        std::vector<int> want;
        for (const Named& in : kInputs) want.push_back(in.par);
        std::vector<int> out{grad::out_profit};
        const grad::RowRequest req{out.data(), out.size(), want.data(),
                                   want.size()};
        grad::Rows rows;
        try {
          rows = grad::rows_differenced(l, kTheta, d, req, s);
        } catch (const std::exception&) {
          continue;
        }
        // The read moves the leaf, so put it back at the point before asking
        // profit_at about it.
        l.find_root_collar_psi();

        // The scale is profit's own response over this input set: a row that is
        // small beside its neighbours is not a disagreement.
        double scale = 0.0;
        for (std::size_t i = 0; i < want.size(); ++i) {
          const double h = rows.held[i];
          if (std::isfinite(h) && std::abs(h) > scale) scale = std::abs(h);
        }
        if (!(scale > 0.0)) continue;

        double worst = 0.0;
        int at = 0;
        const char* worst_name = "-";
        double worst_rows = 0.0, worst_got = 0.0;
        for (std::size_t i = 0; i < want.size(); ++i) {
          const double hand = rows.held[i];
          if (!std::isfinite(hand)) continue;
          const ProfitInputs<tangent> p = inputs_at(l, kInputs[i].par);
          const double got = derivative_along(l.profit_at<tangent>(
              l.opt_root_psi_, l.opt_psi_stem_, l.ci_, p));
          ++compared;
          ++at;
          const double rel = std::abs(got - hand) / scale;
          if (rel > worst) {
            worst = rel; worst_name = kInputs[i].name;
            worst_rows = hand; worst_got = got;
          }
        }
        std::printf("   %-8.3g %-7.4g %d  %-10s %2d %-12s %14.6e %14.6e %11.3e\n",
                    psi_soil, ppfd, layers,
                    Leaf::operating_point_kind_name(l.operating_point_kind()),
                    at, worst_name, worst_rows, worst_got, worst);
        double& into = interior ? worst_interior : worst_pinned;
        if (worst > into) into = worst;
      }
    }
  }
  std::printf("\n   %d rows over %d states.\n", compared, states);
  std::printf("   interior: worst %.3e\n", worst_interior);
  std::printf("   pinned  : worst %.3e   (the collar moves there; profit_at holds it)\n",
              worst_pinned);
  return 0;
}
