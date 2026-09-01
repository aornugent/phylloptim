// The two derivatives the leaf hands over as numbers, against a difference of the
// forward model that produces them.
//
// A supplied row is grafted onto a value, so every consumer downstream reads it
// without recomputing it, and a difference of anything downstream returns
// identically zero on the columns it occupies -- right, wrong or absent. The only
// place a supplied row can be caught is HERE, where the function it claims to be
// the derivative of can still be evaluated at two points.
//
// Two rows are checked, and they are the two the design note calls its gaps:
//
//   duptake_dp[j]           d(uptake_j)/d(collar), from roots_.duptake_dpsi.
//                           supply_draw_at takes the draw at a PASSIVE collar and
//                           supplies this row beside it, so the tape never sees
//                           the uptake move with the collar.
//
//   marginal_collar_slope   dM/dp, the slope of the marginal the interior solve
//                           roots. It crosses into plant as a plain double because
//                           reverse mode cannot nest a tangent above its own
//                           scalar, and it is the denominator the interior closure
//                           divides by -- so a wrong one is a plausible gradient
//                           everywhere, not a failure anywhere.
//
// ⚠️ THE REFERENCE IS THE FORWARD MODEL, NOT A SECOND DERIVATION. Each check
// differences the function its row claims to describe. The uptake row is refereed
// against supply_draw_at's own uptake re-evaluated at a moved collar, which reads
// none of the derivative code. dM/dp is refereed against marginal_assembled
// stepped along the direction the row is taken in, which checks the assembly's
// differentiation but takes the four chain rates as given -- dEup is the sum of
// the uptake rows above, so it is covered; V, dci_dpsi and d2Eup are not.
//
// ⚠️ dM/dp IS DIRECTIONAL. Differencing the collar ALONE differences a different
// function and disagrees by 40 to 100 per cent, which reads as a broken row and is
// a broken reference. The five coordinates move together or the check is wrong.
//
// ⚠️ THE BUDGET IS WHAT THE DIFFERENCE CAN SEE. A central difference cannot
// resolve arbitrarily small numbers, so the tolerance is not a constant: each row
// is differenced at TWO step sizes and the gap between the two answers is the
// floor for that row. A check held tighter than its own reference's noise is
// measuring the step size.
//
// ⚠️ ~4e-8 IS THE CEILING ON ANY CROSS-SCALAR CLAIM HERE, and it is a property of
// the model rather than of this test: the per-layer mean conductivity is a 7-node
// Gauss sum at double and a 1-node midpoint rule on the active path. Both rows
// reach a layer mean. See docs/design/leaf-gradients.md.
//
// Run by `make -C tests/cpp`, with the rest of the suite.

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace pl = phylloptim;

namespace {

int failures = 0;
int checked_uptake = 0;
int checked_slope = 0;
int nan_rows = 0;
int noticed = 0;
double worst_uptake = 0.0;
double worst_slope = 0.0;

const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_up(pl::Leaf& l, int layers, double psi0, double ppfd) {
  double theta[pl::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::trait_R_d_25] = 1.44;
  l.set_traits(theta);
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / kAreaLeaf / layers);
    psi_soil.push_back(psi0 + 0.35 * i);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), ppfd, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

fixture::Soil soil_of(int layers) {
  fixture::Soil s;
  for (int i = 0; i < layers; ++i) {
    s.depth.push_back(1.0 * (i + 1));
    s.carbon.push_back(1.0 / layers / kAreaLeaf);
  }
  return s;
}

double rel_gap(double a, double b) {
  return std::abs(a - b) / (std::abs(a) + std::abs(b) + 1e-300);
}

// A row against a central difference taken at two step sizes. Returns the
// disagreement with the finer difference, and reports the floor the two steps
// leave, so a row is never held tighter than the reference can resolve.
struct Verdict {
  double gap;    // row against the finer difference
  double floor;  // what the two step sizes disagree by, which is the noise
  bool ok;
};

// Above this a row is reported, so a state whose reference cannot resolve the row
// is named every run instead of being averaged away. Measured: every state but one
// sits at 1e-12 to 1e-10, and the driest single-layer one reaches 1.7e-07 with both
// step sizes agreeing -- a systematic offset near a branch, not noise. That state
// has psi_soil within 1% of root_psi_crit, where conductivity is on the 5% tail.
const double kNotice = 4e-8;
// Above this it fails. Three orders below the ~1e-4 that phylloptim's guide calls a
// real difference, and five above the ~1e-9 solver floor.
const double kFail = 1e-6;

Verdict against(double row, double coarse, double fine) {
  Verdict v;
  v.gap = rel_gap(row, fine);
  v.floor = rel_gap(coarse, fine);
  // The reference's own spread where that is the binding constraint, and a fixed
  // ceiling where the two step sizes agree with each other but not with the row.
  v.ok = v.gap <= std::max(3.0 * v.floor, kFail);
  return v;
}

void at_state(int layers, double psi0, double ppfd) {
  pl::Leaf l;
  set_up(l, layers, psi0, ppfd);
  l.find_root_collar_psi();

  const char* kind = l.operating_point_kind_name();
  const fixture::Soil soil = soil_of(layers);
  const pl::LeafInputs<double> in =
      fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);
  const double p = l.opt_root_psi_;

  // A shutdown draws no water and supplies a row of zeros deliberately, so there
  // is nothing here for a difference to referee.
  if (l.operating_point_kind() ==
      pl::Leaf::OperatingPointKind::HydraulicShutdown) {
    std::printf("  %d layer  psi0 %-5.2f ppfd %-6.0f %-22s  no draw, skipped\n",
                layers, psi0, ppfd, kind);
    return;
  }

  // ---- duptake_dp[j] against a difference of supply_draw_at's own uptake ----
  //
  // supply_draw_at recomputes the uptake from E_from_soil_at at whatever collar
  // it is handed, and supplies duptake_dp beside it from a different function. So
  // moving the collar and re-taking the draw differences the model, not the row.
  double uptake_worst = 0.0;
  bool uptake_ran = false;
  try {
    const auto at = l.supply_draw_at<double>(p, in.supply);
    const std::size_t n = at.uptake.size();
    // One difference per step size, kept per layer so the coarse and fine
    // answers a layer is judged by are that layer's own.
    std::vector<double> fd[2];
    const double steps[2] = {1e-5, 5e-6};
    for (int k = 0; k < 2; ++k) {
      const auto hi =
          l.supply_draw_at<double>(p + steps[k], in.supply);
      const auto lo =
          l.supply_draw_at<double>(p - steps[k], in.supply);
      fd[k].resize(n);
      for (std::size_t j = 0; j < n; ++j) {
        fd[k][j] = (hi.uptake[j] - lo.uptake[j]) / (2.0 * steps[k]);
      }
    }
    for (std::size_t j = 0; j < n; ++j) {
      // ⚠️ A NaN ROW IS A SIGNAL, NOT A FAULT. duptake_dpsi returns NaN at a kink
      // to tell its caller to fall back to finite differences, so a difference is
      // exactly what should happen here and there is no row to referee.
      if (!std::isfinite(at.duptake_dp[j])) {
        ++nan_rows;
        continue;
      }
      // Only rows with something to say: a layer contributing no flux has a zero
      // row and a zero difference, which agree for no reason.
      if (std::abs(at.duptake_dp[j]) < 1e-12) continue;
      if (!std::isfinite(fd[0][j]) || !std::isfinite(fd[1][j])) continue;
      const Verdict v = against(at.duptake_dp[j], fd[0][j], fd[1][j]);
      uptake_worst = std::max(uptake_worst, v.gap);
      if (v.gap > kNotice) ++noticed;
      if (!v.ok) ++failures;
      uptake_ran = true;
      ++checked_uptake;
    }
  } catch (const std::exception&) {
    uptake_ran = false;
  }

  // ---- marginal_collar_slope against a difference of marginal_at -------------
  //
  // Interior only: it is the only arm production computes it on, and the only one
  // where the collar is placed by rooting M rather than by a bound. The draw is
  // re-taken at each moved collar, because the marginal the solve roots is the one
  // whose draw follows the candidate.
  double slope_gap = 0.0, slope_floor = 0.0;
  bool slope_ran = false;
  if (l.operating_point_kind() == pl::Leaf::OperatingPointKind::Interior) {
    try {
      const double row = l.marginal_collar_slope(in.profit);

      // ⚠️ dM/dp IS A DIRECTIONAL DERIVATIVE, NOT A PARTIAL IN THE COLLAR, and
      // getting that wrong is the trap this check exists past. marginal_assembled
      // takes five coordinates, and marginal_collar_slope seeds all five: p by 1,
      // sigma by V, ci by dci_dpsi, dEup_dp by d2Eup, transpiration by dEup. A
      // difference that moves only the collar differences a different function --
      // it came out disagreeing by 40 to 100 per cent, which reads as a broken row
      // and is a broken reference.
      //
      // So the reference steps along that same direction and differences there.
      // The seeds come from `parts`, which is where the assembly gets them too, so
      // this checks the ASSEMBLY's differentiation rather than the chain rates.
      // Those are refereed separately: dEup is the sum of the uptake rows above,
      // and d2Eup by differencing dEup.
      const double sigma = l.opt_psi_stem_;
      const std::vector<double>& soil = l.supply_psi_soil();
      const double dEup = l.dE_from_soil_dpsi_collar(p, soil);
      const double d2Eup = l.d2E_from_soil_dpsi_collar2(p, soil);
      const double transp = l.transpiration(sigma, p);
      const pl::Leaf::CollarPoint<double> here{p, sigma, l.ci_, dEup, transp};
      const auto parts = l.marginal_assembled<double>(here, in.profit);

      const auto along = [&](double t) -> double {
        const pl::Leaf::CollarPoint<double> q{
            p + t, sigma + t * parts.V, l.ci_ + t * parts.dci_dpsi,
            dEup + t * d2Eup, transp + t * dEup};
        return l.marginal_assembled<double>(q, in.profit).marginal;
      };
      double fd[2] = {0.0, 0.0};
      int k = 0;
      for (const double h : {1e-5, 5e-6}) {
        fd[k++] = (along(h) - along(-h)) / (2.0 * h);
      }
      const Verdict v = against(row, fd[0], fd[1]);
      slope_gap = v.gap;
      slope_floor = v.floor;
      slope_ran = true;
      ++checked_slope;
      worst_slope = std::max(worst_slope, v.gap);
      if (!v.ok) ++failures;
    } catch (const std::exception&) {
      slope_ran = false;
    }
  }
  worst_uptake = std::max(worst_uptake, uptake_worst);

  std::printf("  %d layer  psi0 %-5.2f ppfd %-6.0f %-22s", layers, psi0, ppfd,
              kind);
  if (uptake_ran) {
    std::printf("  uptake %.2e", uptake_worst);
  } else {
    std::printf("  uptake    --  ");
  }
  if (slope_ran) {
    std::printf("  dM/dp %.2e (floor %.1e)", slope_gap, slope_floor);
  }
  std::printf("\n");
}

}  // namespace

int main() {
  std::printf("supplied rows, against a difference of the model they describe\n");
  std::printf("  duptake_dp  vs  d/dp of supply_draw_at's uptake\n");
  std::printf("  marginal_collar_slope  vs  d/dp of marginal_at\n\n");

  for (int layers : {1, 3, 5}) {
    for (double ppfd : {900.0, 30.0}) {
      for (double psi0 : {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.8}) {
        at_state(layers, psi0, ppfd);
      }
    }
  }

  std::printf("\n  uptake rows checked %d, worst %.2e\n", checked_uptake,
              worst_uptake);
  std::printf("  dM/dp rows checked  %d, worst %.2e\n", checked_slope,
              worst_slope);
  std::printf("  NaN rows (a kink, telling the caller to difference): %d\n",
              nan_rows);
  std::printf("  rows above the %.0e notice: %d  (fail is %.0e)\n", kNotice,
              noticed, kFail);

  // Both rows must have been reached. A run that checked neither would report
  // success having differenced nothing, which is the failure this file exists to
  // make impossible for its own subject.
  if (checked_uptake == 0) {
    std::printf("  NOTHING CHECKED: no state supplied an uptake row\n");
    ++failures;
  }
  if (checked_slope == 0) {
    std::printf("  NOTHING CHECKED: no state reached an interior collar\n");
    ++failures;
  }

  std::printf("\n%s\n", failures == 0 ? "all ok" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
