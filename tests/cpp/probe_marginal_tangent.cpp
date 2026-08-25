// PROBE (not a test, not shipped): the collar's own condition, from profit_at.
//
// The marginal profit M(p; theta) = dprofit/dp is what places an interior collar,
// and every water output reads the collar. Today M comes from
// dprofit_at_collar_psi (93 code lines of hand-written chain rule) and its
// derivatives from condition_curvature (68 more, a hand-transcribed second-order
// theorem). This asks whether both are one seeded pass over profit_at.
//
// A tangent above a tangent gives all four at once. Seed the collar in the inner
// layer and one input in the outer:
//
//   value().value()            profit
//   value().derivative()       M          = dprofit/dp
//   derivative().value()       dprofit/dtheta at a held collar
//   derivative().derivative()  dM/dtheta  -- the condition's own gradient
//
// and with the collar seeded in BOTH layers, derivative().derivative() is dM/dp,
// the curvature the theorem divides by.
//
// ⚠️ TWO SCORING REGIMES, AND THE SPLIT IS RULE 6. At the optimum M IS zero -- it
// is the condition the solve drove to nothing -- so scoring M there compares two
// zeros. M is scored ACROSS the feasible interval against dprofit_at_collar_psi.
// dM/dtheta and dM/dp are not zero there, so they are scored AT the point,
// against what the row layer reports.
//
//   make CXX=g++ ODELIA_INC=... probe_marginal_tangent && ./probe_marginal_tangent

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
using phylloptim::tangent2;
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

// The supply at one scalar, owned so the view has something to point at.
template <class T>
struct SupplyStore {
  std::vector<T> psi_soil, r_R_H_min, r_R_V_sum;
  T root_b, root_c;
  phylloptim::SupplyAt<T> view() const {
    return {psi_soil, r_R_H_min, r_R_V_sum, root_b, root_c};
  }
};

template <class T>
SupplyStore<T> supply_of(const Leaf& l, int par, int n_layers) {
  const phylloptim::SupplyAt<double> held = l.held_supply();
  SupplyStore<T> out;
  for (double v : held.psi_soil) out.psi_soil.push_back(T(v));
  for (double v : held.r_R_H_min) out.r_R_H_min.push_back(T(v));
  for (double v : held.r_R_V_sum) out.r_R_V_sum.push_back(T(v));
  out.root_b = T(held.root_b);
  out.root_c = T(held.root_c);
  if (par == grad::par_root_b) seed_direction(out.root_b, 1.0);
  if (par == grad::par_root_c) seed_direction(out.root_c, 1.0);
  const int layer = par - grad::par_psi_soil_first;
  if (layer >= 0 && layer < n_layers &&
      layer < static_cast<int>(out.psi_soil.size())) {
    seed_direction(out.psi_soil[std::size_t(layer)], 1.0);
  }
  return out;
}

const Named kInputs[] = {
    {grad::par_vcmax_25, "vcmax_25"}, {grad::par_jmax_25, "jmax_25"},
    {grad::par_a, "a"},               {grad::par_curv_fact_elec_trans, "curv_elec"},
    {grad::par_curv_fact_colim, "curv_colim"}, {grad::par_R_d_25, "R_d_25"},
    {grad::par_PPFD, "PPFD"},         {grad::par_stem_b, "stem_b"},
    {grad::par_stem_c, "stem_c"},     {grad::par_beta2, "beta2"},
    {grad::par_cost_scale_TF24, "cost_scale"}, {grad::par_kmax, "kmax"},
    {grad::par_root_b, "root_b"},     {grad::par_root_c, "root_c"},
};

// The inputs at the leaf's current state and a given collar, at whatever scalar.
// `seed` names the outer direction; the collar's own is set by the caller, which
// is what separates dM/dtheta from dM/dp.
template <class T>
ProfitInputs<T> inputs_at(const Leaf& l, const T& collar, int seed_par) {
  ProfitInputs<T> p{
      T(l.vcmax_),   T(l.jmax_),   T(l.a),      T(l.curv_fact_elec_trans),
      T(l.curv_fact_colim), T(l.PPFD_), T(l.R_d_),
      T(l.leaf_specific_conductance_max_), T(l.stem_b), T(l.stem_c),
      T(l.beta2),    T(l.cost_scale_TF24),  collar,     T(0.0)};
  switch (seed_par) {
    case grad::par_vcmax_25: seed_direction(p.vcmax, l.vcmax_ / l.vcmax_25); break;
    case grad::par_jmax_25: seed_direction(p.transport_jmax, l.jmax_ / l.jmax_25); break;
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
  std::printf("the collar's condition, from profit_at at a nested tangent\n\n");

  // --- M, away from the optimum -------------------------------------------
  std::printf("  M = dprofit/dp, across the feasible interval\n");
  std::printf("   %-8s %-7s %s %3s %15s %15s %11s\n", "psi_soil", "PPFD", "L",
              "n", "hand", "profit_at", "worst rel");
  double worst_M = 0.0;
  int M_compared = 0;

  for (double psi_soil : {0.2, 0.8, 1.5, 2.5, 3.5, 4.5}) {
    for (double ppfd : {30.0, 300.0, 1000.0, 2000.0}) {
      for (int layers : {1, 3}) {
        const grad::Drivers d = drivers(psi_soil, ppfd, layers);
        Leaf l;
        const grad::Settings s;
        grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
        try { l.find_root_collar_psi(); } catch (const std::exception&) { continue; }
        double a = 0.0, b = 0.0;
        if (!l.prepare_collar_solve(a, b)) continue;

        double worst = 0.0, wh = 0.0, wg = 0.0;
        int at = 0;
        for (double frac : {0.05, 0.2, 0.4, 0.6, 0.8, 0.95}) {
          const double p = a + frac * (b - a);
          bool feasible = false;
          const double hand = l.dprofit_at_collar_psi(p, &feasible);
          if (!feasible || !std::isfinite(hand)) continue;
          // Place the point at p so sigma* and ci* describe it.
          if (!l.profit_at_fixed_collar(p).feasible) continue;
          const double sigma = l.opt_psi_stem_, ci = l.ci_;
          tangent pc = p;
          seed_direction(pc, 1.0);
          const ProfitInputs<tangent> in = inputs_at<tangent>(l, pc, -1);
          const SupplyStore<tangent> sup = supply_of<tangent>(l, -1, layers);
          const double got = derivative_along(
              l.profit_at<tangent>(sigma, ci, in, sup.view()));
          ++M_compared; ++at;
          const double rel = std::abs(hand) > 0.0
                                 ? std::abs(got - hand) / std::abs(hand)
                                 : std::abs(got);
          if (rel > worst) { worst = rel; wh = hand; wg = got; }
        }
        if (at == 0) continue;
        std::printf("   %-8.3g %-7.4g %d %3d %15.7e %15.7e %11.3e\n", psi_soil,
                    ppfd, layers, at, wh, wg, worst);
        if (worst > worst_M) worst_M = worst;
      }
    }
  }
  std::printf("\n   %d collars. worst %.3e\n\n", M_compared, worst_M);

  // --- dM/dp and dM/dtheta, at the point -----------------------------------
  std::printf("  dM/dp and dM/dtheta, at the operating point, against the row layer\n");
  std::printf("   %-8s %-7s %s %-10s %11s %-12s %11s\n", "psi_soil", "PPFD", "L",
              "point", "curvature", "worst input", "worst rel");
  double worst_curv = 0.0, worst_grad = 0.0;
  int grad_compared = 0;
  // An input whose reference row is zero everywhere agrees perfectly and proves
  // nothing, so each one's live count is reported beside the worst.
  std::vector<int> live(sizeof(kInputs) / sizeof(kInputs[0]) + 1, 0);

  for (double psi_soil : {0.2, 0.8, 1.5, 2.5, 3.5, 4.5}) {
    for (double ppfd : {30.0, 300.0, 1000.0, 2000.0}) {
      for (int layers : {1, 3}) {
        const grad::Drivers d = drivers(psi_soil, ppfd, layers);
        Leaf l;
        const grad::Settings s;
        grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
        try { l.find_root_collar_psi(); } catch (const std::exception&) { continue; }
        if (l.operating_point_kind() != Leaf::OperatingPointKind::Interior) continue;

        std::vector<int> want;
        for (const Named& in : kInputs) want.push_back(in.par);
        want.push_back(grad::par_psi_soil_first);
        std::vector<const char*> names;
        for (const Named& in : kInputs) names.push_back(in.name);
        names.push_back("psi_soil_1");
        std::vector<int> out{grad::out_profit};
        const grad::RowRequest req{out.data(), out.size(), want.data(), want.size()};
        grad::Rows rows;
        try { rows = grad::rows_differenced(l, kTheta, d, req, s); }
        catch (const std::exception&) { continue; }
        l.find_root_collar_psi();
        const double p = l.opt_root_psi_;
        if (!l.profit_at_fixed_collar(p).feasible) continue;
        const double sigma = l.opt_psi_stem_, ci = l.ci_;

        // dM/dp: the collar seeded in both layers.
        tangent2 pc2 = p;
        pc2.value().derivative() = 1.0;
        pc2.derivative().value() = 1.0;
        const SupplyStore<tangent2> sup2 = supply_of<tangent2>(l, -1, layers);
        const tangent2 rr = l.profit_at<tangent2>(
            sigma, ci, inputs_at<tangent2>(l, pc2, -1), sup2.view());
        const double curv = rr.derivative().derivative();
        const double curv_rel =
            std::abs(rows.residual_slope) > 0.0
                ? std::abs(curv - rows.residual_slope) / std::abs(rows.residual_slope)
                : std::abs(curv);
        if (curv_rel > worst_curv) worst_curv = curv_rel;

        // dM/dtheta: the collar inner, one input outer.
        double scale = 0.0;
        for (std::size_t i = 0; i < want.size(); ++i) {
          const double v = rows.dresidual[i];
          if (std::isfinite(v) && std::abs(v) > scale) scale = std::abs(v);
        }
        double worst = 0.0; const char* wn = "-";
        if (scale > 0.0) {
          for (std::size_t i = 0; i < want.size(); ++i) {
            const double hand = rows.dresidual[i];
            if (!std::isfinite(hand)) continue;
            tangent2 pin = p;
            pin.value().derivative() = 1.0;   // inner: d/dp
            const SupplyStore<tangent2> sp =
                supply_of<tangent2>(l, want[i], layers);
            const tangent2 r = l.profit_at<tangent2>(
                sigma, ci, inputs_at<tangent2>(l, pin, want[i]), sp.view());
            const double got = r.derivative().derivative();
            ++grad_compared;
            if (hand != 0.0) ++live[i];
            const double rel = std::abs(got - hand) / scale;
            if (rel > worst) { worst = rel; wn = names[i]; }
          }
        }
        if (worst > worst_grad) worst_grad = worst;
        std::printf("   %-8.3g %-7.4g %d %-10s %11.3e %-12s %11.3e\n", psi_soil,
                    ppfd, layers,
                    Leaf::operating_point_kind_name(l.operating_point_kind()),
                    curv_rel, wn, worst);
      }
    }
  }
  std::printf("\n   states in which each condition row is non-zero:\n    ");
  {
    std::vector<const char*> nm;
    for (const Named& in : kInputs) nm.push_back(in.name);
    nm.push_back("psi_soil_1");
    for (std::size_t i = 0; i < nm.size(); ++i) std::printf("%s=%d ", nm[i], live[i]);
  }
  std::printf("\n\n   %d condition rows.\n", grad_compared);
  std::printf("   dM/dp     worst %.3e\n", worst_curv);
  std::printf("   dM/dtheta worst %.3e\n", worst_grad);
  return 0;
}
