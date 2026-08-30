// The stem curve's round trip: G^-1(G(psi)) against psi. At zero flux this IS the
// question transpiration_to_psi_stem asks, and the sign of its error decided a
// feasibility classification.
#include <phylloptim.hpp>
#include "root_network.hpp"
#include <cstdio>
#include <cmath>
int main() {
  phylloptim::Leaf l;
  double worst = 0.0, worst_at = 0.0;
  int wrong_side = 0, n = 0;
  const double lo = 0.05, hi = l.psi_crit * 0.98;
  for (int i = 1; i < 4000; ++i) {
    const double psi = lo + (hi - lo) * i / 4000.0;
    const double w = l.stem_curve_integral(psi, "probe");
    const double back = l.stem_curve_integral_inverse(w, "probe");
    const double rel = std::abs(back - psi) / psi;
    if (rel > worst) { worst = rel; worst_at = psi; }
    if (back > psi) ++wrong_side;   // the side that reads as feasible
    ++n;
  }
  std::printf("round trip over %d potentials: worst relative %.3e at psi = %.4f\n",
              n, worst, worst_at);
  std::printf("  came back ABOVE the potential asked for: %d of %d\n", wrong_side, n);
  return 0;
}
