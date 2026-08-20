// PROBE (not a test, not shipped): the wet bound at a shade death whose profile
// straddles psi_crit. `rows_at` falls back to differencing the whole solve wherever
// this row refuses, and that fallback is the last consumer of a step in the row
// layer -- so whether it refuses, and why, decides the signature.
#include <phylloptim.hpp>
#include "root_network.hpp"
#include <cstdio>
#include <exception>
#include <vector>
namespace grad = phylloptim::gradient;
namespace {
const double kRd25 = 1.44;
double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
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
  d.PPFD = ppfd; d.psi_soil = ps; d.soil_depth = depth; d.atm_vpd = 2.0;
  d.ca = 40.0; d.leaf_temp = 25.0; d.atm_o2_kpa = 21.0; d.atm_kpa = 101.3;
  return d;
}
void one(double ps, double ppfd, int L, double pc) {
  double th[16];
  for (int i = 0; i < 16; ++i) th[i] = kTheta[i];
  th[grad::par_psi_crit] = pc;
  const grad::Drivers d = drivers(ps, ppfd, L);
  phylloptim::Leaf l;
  const grad::Settings s;
  grad::apply(l, th, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  std::printf("psi_soil=%.2f..%.2f psi_crit=%.4f L=%d  %-16s seat=%.10g\n", ps,
              ps + 0.25 * (L - 1), pc, L,
              phylloptim::Leaf::operating_point_kind_name(l.operating_point_kind()),
              l.opt_root_psi_);
  try {
    const phylloptim::Leaf::BoundRow w =
        l.bound_row(phylloptim::Leaf::WhichBound::Wet);
    std::printf("   bound_row(Wet): bound=%.10g finite=%d slope=%.6g "
                "seat-match=%d\n",
                w.bound, int(w.finite), w.residual_slope,
                int(w.bound == l.opt_root_psi_));
  } catch (const std::exception& e) {
    std::printf("   bound_row(Wet) THREW: %s\n", e.what());
  }
}
}  // namespace
int main() {
  // The fifteen: five layers spanning psi_crit, shade-dead.
  one(6.00, 0.0, 5, 6.5);
  one(6.00, 10.0, 5, 6.5);
  one(6.00, 30.0, 5, 6.5);
  // Neighbours, to see where it starts refusing.
  for (double ps = 5.50; ps <= 6.51; ps += 0.25) one(ps, 10.0, 5, 6.5);
  return 0;
}
