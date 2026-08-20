// PROBE (not a test, not shipped): two disagreements, each refereed by a
// SEQUENCE of steps rather than by one.
//
//   A. the closed-form curvature against a differenced marginal profit, over a
//      decade sweep of the collar step. A central difference divides the ci
//      root-find's floor by the step, so a single step reports its own noise; a
//      plateau, if there is one, says which number the difference is converging on.
//   B. the wet bound's own soil row against a differenced find_root_psi, the same
//      way. The wet bound is where total uptake vanishes, so the quotient there
//      divides by a small slope.
//
//   make CXX=g++ probe_plateau && stdbuf -oL ./probe_plateau

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace grad = phylloptim::gradient;
using Kind = grad::OperatingPointKind;

namespace {

const double kRd25 = 1.44;
const double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,      157.44,   0.30,
                         0.7,      0.99,     7.5,      kRd25,
                         1.0 * 0.000157 / 5.0, 1e3};

grad::Drivers drivers(double psi_soil, double ppfd, double vpd, int layers) {
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
  d.atm_vpd = vpd;
  d.ca = 40.0;
  d.leaf_temp = 25.0;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

void curvature_plateau(double psi_soil, double ppfd, double vpd, int layers) {
  const grad::Drivers d = drivers(psi_soil, ppfd, vpd, layers);
  phylloptim::Leaf l;
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  const Kind k = l.operating_point_kind();
  const double p = l.opt_root_psi_;
  const double r0 = l.dprofit_droot_collar_psi(p);
  double closed = 0.0;
  const bool ok = l.condition_collar_slope(closed);
  std::printf("\nA. psi_soil=%.2f ppfd=%.0f vpd=%.1f L=%d  %s  p*=%.10g "
              "resid=%.3g\n",
              psi_soil, ppfd, vpd, layers,
              phylloptim::Leaf::operating_point_kind_name(k), p, r0);
  std::printf("   closed form Pi_pp = %.12g (ok=%d)\n", closed, int(ok));
  std::printf("   %-12s %-18s %-12s\n", "collar step", "differenced", "rel");
  for (int e = 3; e <= 10; ++e) {
    const double h = std::max(std::abs(p), 1.0) * std::pow(10.0, -double(e));
    l.dprofit_droot_collar_psi(p);
    const double resid = l.dprofit_droot_collar_psi(p);
    const double d_hi = l.dprofit_droot_collar_psi(p + h);
    const double d_lo = l.dprofit_droot_collar_psi(p - h);
    const double H = (d_hi - d_lo) / (2.0 * h);
    std::printf("   1e-%-9d %-18.12g %-12.3g%s\n", e, H,
                std::abs(H / closed - 1.0), e == 6 ? "   <- the one in use" : "");
    static_cast<void>(resid);
  }
}

void bound_plateau(double psi_soil, double ppfd, double vpd, int layers) {
  const grad::Drivers d = drivers(psi_soil, ppfd, vpd, layers);
  phylloptim::Leaf l;
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  const Kind k = l.operating_point_kind();
  const phylloptim::Leaf::BoundRow wet =
      l.bound_row(phylloptim::Leaf::WhichBound::Wet);
  std::printf("\nB. psi_soil=%.2f ppfd=%.0f vpd=%.1f L=%d  %s\n", psi_soil, ppfd,
              vpd, layers, phylloptim::Leaf::operating_point_kind_name(k));
  std::printf("   wet bound = %.12g   dR/dx = %.6g   collar_root_tol = %.3g\n",
              wet.bound, wet.residual_slope, phylloptim::Leaf::collar_root_tol);
  std::printf("   closed form d(bound)/d(psi_soil_1) = %.12g\n",
              wet.d_dpsi_soil[0]);
  std::printf("   %-12s %-18s %-12s\n", "step", "differenced", "rel");
  for (int e = 3; e <= 9; ++e) {
    const double h = std::pow(10.0, -double(e));
    double b[2];
    for (int side = 0; side < 2; ++side) {
      grad::Drivers moved = d;
      moved.psi_soil[0] += side == 0 ? h : -h;
      grad::apply(l, kTheta, moved, false, -1, s.fast_stem_curve);
      b[side] = l.find_root_psi(moved.psi_soil[0], moved.psi_soil, 0);
    }
    const double got = (b[0] - b[1]) / (2.0 * h);
    std::printf("   1e-%-9d %-18.12g %-12.3g%s\n", e, got,
                std::abs(got / wet.d_dpsi_soil[0] - 1.0),
                e == 6 ? "   <- the one in use" : "");
  }
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
}

}  // namespace

int main() {
  std::printf("=== A. THE CURVATURE ===\n");
  // The census's worst point, then two others for shape.
  curvature_plateau(0.5, 900.0, 2.0, 3);
  curvature_plateau(2.0, 900.0, 2.0, 5);
  curvature_plateau(0.5, 1500.0, 4.0, 1);

  std::printf("\n=== B. THE WET BOUND'S SOIL ROW ===\n");
  // A wet pin from the golden grid's own soil range, and an interior point for
  // contrast -- the bound exists whether or not the optimum is pinned to it.
  bound_plateau(3.5, 900.0, 2.0, 5);
  bound_plateau(4.0, 100.0, 4.0, 5);
  bound_plateau(2.0, 900.0, 2.0, 5);
  return 0;
}
