// PROBE (not a test, not shipped): the bounds that pin a collar, against
// bound_row.
//
// A pinned point sits where a limit stops it, and its row is that limit's own.
// bound_row is 98 code lines of hand-written quotient across three arms, one of
// which carries a transcription of Euler's identity for the stem curve. This
// asks whether all three are ordinary implicit values over functions the model
// already has:
//
//   Wet             E_up(x) = 0                            uptake vanishes
//   DryRootCrit     kmax*(G(psi_crit) - G(x)) - E_up(x) = 0  T1 at the stem's limit
//   DryRootPsiCrit  x = root_psi_crit                       a registered constant
//
// ⚠️ WHAT IS SCORED IS THE ROW, NOT THE RESIDUAL. Each residual is zero at its
// own bound by construction, so comparing residuals there compares two zeros --
// rule 6. The rows are not zero there, and they are what a consumer reads.
//
//   make CXX=g++ ODELIA_INC=... probe_bound_tangent && ./probe_bound_tangent

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace grad = phylloptim::gradient;
using phylloptim::Leaf;
using phylloptim::LeafInputs;
using phylloptim::tangent;
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
    {grad::par_kmax, "kmax"},       {grad::par_stem_b, "stem_b"},
    {grad::par_stem_c, "stem_c"},   {grad::par_psi_crit, "psi_crit"},
    {grad::par_root_b, "root_b"},   {grad::par_root_c, "root_c"},
    {grad::par_root_psi_crit, "root_psi_crit"},
};

// The entry of bound_row that answers for one input.
double hand_row(const Leaf::BoundRow& r, int par, int n_layers) {
  switch (par) {
    case grad::par_kmax: return r.d_dkappa;
    case grad::par_stem_b: return r.d_dstem_b;
    case grad::par_stem_c: return r.d_dstem_c;
    case grad::par_psi_crit: return r.d_dpsi_crit;
    case grad::par_root_b: return r.d_droot_b;
    case grad::par_root_c: return r.d_droot_c;
    case grad::par_root_psi_crit: return r.d_droot_psi_crit;
    default: break;
  }
  const int layer = par - grad::par_psi_soil_first;
  if (layer >= 0 && layer < n_layers &&
      layer < static_cast<int>(r.d_dpsi_soil.size())) {
    return r.d_dpsi_soil[std::size_t(layer)];
  }
  return 0.0;
}

const char* bound_name(Leaf::WhichBound w) {
  switch (w) {
    case Leaf::WhichBound::Wet: return "wet";
    case Leaf::WhichBound::DryRootCrit: return "dry-root-crit";
    default: return "dry-root-psi-crit";
  }
}

}  // namespace

int main() {
  std::printf("the bounds that pin a collar, against bound_row\n\n");
  std::printf("   %-8s %-7s %s %-16s %8s %11s %-14s %11s\n", "psi_soil", "PPFD",
              "L", "bound", "value", "slope rel", "worst input", "worst rel");

  double worst[3] = {0.0, 0.0, 0.0};
  double worst_value = 0.0, worst_slope = 0.0;
  int compared = 0, states = 0;
  std::vector<int> live(sizeof(kInputs) / sizeof(kInputs[0]) + 1, 0);

  for (double psi_soil : {0.2, 0.8, 1.5, 2.5, 3.5, 4.5}) {
    for (double ppfd : {30.0, 300.0, 1000.0, 2000.0}) {
      for (int layers : {1, 3}) {
        const grad::Drivers d = drivers(psi_soil, ppfd, layers);
        Leaf l;
        const grad::Settings s;
        grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
        try { l.find_root_collar_psi(); } catch (const std::exception&) { continue; }
        ++states;

        const Leaf::WhichBound arms[] = {Leaf::WhichBound::Wet,
                                         Leaf::WhichBound::DryRootCrit,
                                         Leaf::WhichBound::DryRootPsiCrit};
        for (int a = 0; a < 3; ++a) {
          Leaf::BoundRow row;
          try { row = l.bound_row(arms[a]); }
          catch (const std::exception&) { continue; }
          if (!row.finite || !std::isfinite(row.bound)) continue;

          std::vector<int> want;
          for (const Named& in : kInputs) want.push_back(in.par);
          for (int j = 0; j < layers; ++j) {
            want.push_back(grad::par_psi_soil_first + j);
          }
          std::vector<const char*> names;
          for (const Named& in : kInputs) names.push_back(in.name);
          names.push_back("psi_soil");

          double scale = 0.0;
          for (int par : want) {
            const double v = hand_row(row, par, layers);
            if (std::isfinite(v) && std::abs(v) > scale) scale = std::abs(v);
          }
          if (!(scale > 0.0)) continue;

          // The value first: a bound that comes back at a different place is a
          // different bound, and every row below would then be about it.
          const LeafInputs<double> passive =
              fixture::leaf_inputs<double>(l, row.bound, -1);
          const double got_value =
              l.bound_at<double>(arms[a], row.bound, passive);
          const double vrel = std::abs(got_value - row.bound) /
                              (std::abs(row.bound) + 1e-300);
          if (vrel > worst_value) worst_value = vrel;

          double w = 0.0;
          const char* wn = "-";
          for (std::size_t i = 0; i < want.size(); ++i) {
            const double hand = hand_row(row, want[i], layers);
            if (!std::isfinite(hand)) continue;
            if (hand != 0.0 && i < live.size()) ++live[i];
            const LeafInputs<tangent> in =
                fixture::leaf_inputs<tangent>(l, tangent(row.bound), want[i]);
            const double got =
                derivative_along(l.bound_at<tangent>(arms[a], row.bound, in));
            ++compared;
            const double rel = std::abs(got - hand) / scale;
            if (rel > w) { w = rel; wn = names[i < names.size() ? i : names.size() - 1]; }
          }
          if (w > worst[a]) worst[a] = w;
          std::printf("   %-8.3g %-7.4g %d %-16s %8.4f %11.3e %-14s %11.3e\n",
                      psi_soil, ppfd, layers, bound_name(arms[a]), row.bound,
                      vrel, wn, w);
        }
      }
    }
  }
  std::printf("\n   states in which each row is non-zero on some bound:\n    ");
  {
    std::vector<const char*> nm;
    for (const Named& in : kInputs) nm.push_back(in.name);
    nm.push_back("psi_soil");
    for (std::size_t i = 0; i < nm.size(); ++i) std::printf("%s=%d ", nm[i], live[i]);
  }
  std::printf("\n\n   %d rows over %d states.\n", compared, states);
  std::printf("   bound value    worst %.3e\n", worst_value);
  std::printf("   wet            worst %.3e\n", worst[0]);
  std::printf("   dry-root-crit  worst %.3e\n", worst[1]);
  std::printf("   dry-psi-crit   worst %.3e\n", worst[2]);
  return 0;
}
