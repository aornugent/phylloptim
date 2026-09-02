// The row the leaf hands over as a number, against a difference of the forward
// model that produces it.
//
// A supplied row is grafted onto a value, so every consumer downstream reads it
// without recomputing it -- and a difference of anything downstream returns
// IDENTICALLY ZERO on the column it occupies, right, wrong or absent. The only
// place it can be caught is here, where the function it claims to be the
// derivative of can still be evaluated at two points.
//
// ONE row is supplied on this surface now:
//
//   duptake_dp[j]   d(uptake_j)/d(collar), from roots_.duptake_dpsi.
//                   supply_draw_at takes the draw at a PASSIVE collar and
//                   supplies this row beside it, so the tape never sees a
//                   layer's uptake move with the collar. outputs_at grafts it
//                   as `uptake[j] + duptake_dp[j] * step`.
//
// ⚠️ THE AGGREGATE IS NOT THE ROWS. duptake_dpsi's total is refereed against a
// difference of uptake_at elsewhere, and it passed at 6e-12 while a per-layer
// term could still be wrong in a way the sum hides -- two layers off by equal and
// opposite amounts cancel exactly. Each layer is differenced separately here for
// that reason.
//
// ⚠️ THE ROWS ARE IN MOL AND THE AGGREGATE IS IN KG. Their sum times
// kg_per_mol_h2o is the total, and asserting that is what caught the dropped
// constant once: every per-layer value stayed exact and only the aggregate moved.
//
// ⚠️ dM/dp IS NO LONGER A SUPPLIED ROW and is not checked here. It is a central
// difference of the model's own marginal, so differencing it against the model
// would compare a difference with itself. The conditioning measurement that made
// it a difference is at marginal_collar_slope.
//
// ⚠️ THE BUDGET IS WHAT THE DIFFERENCE CAN SEE. A central difference cannot
// resolve arbitrarily small numbers, so each row is differenced at TWO steps and
// the gap between those answers is the floor for that row. A check held tighter
// than its own reference's noise is measuring the step size.
//
// Run by `make -C tests/cpp`, with the rest of the suite.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace pl = phylloptim;

namespace {

int failures = 0;
int compared = 0;
double worst_ratio = 0.0;   // |analytic - differenced| / that row's own floor
double worst_sum = 0.0;

void ok(bool cond, const std::string& what) {
  if (!cond) {
    ++failures;
    printf("FAIL %s\n", what.c_str());
  }
}

struct Drivers {
  double PPFD = 1500, area_leaf = 1.0, K_s = 2.0, theta = 0.5, h = 5.0;
  double atm_vpd = 2.0, ca = 40.0, leaf_temp = 25.0, atm_o2_kpa = 21.0,
         atm_kpa = 101.3;
};

pl::Leaf make_leaf(const Drivers& d, const std::vector<double>& psi_soil,
                   const std::vector<double>& soil_depth, std::size_t rooted) {
  pl::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> mass(psi_soil.size(),
                           1.0 / double(psi_soil.size()) / d.area_leaf);
  // Roots that stop short of the profile, which is the case a fixture with
  // carbon in every layer cannot see.
  for (std::size_t i = rooted; i < mass.size(); ++i) mass[i] = 0.0;
  l.set_physiology(fixture::root_network(mass, soil_depth), d.PPFD, psi_soil,
                   soil_depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                   d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// Layer j's uptake at a collar, which reads none of the derivative code.
double uptake_of(const pl::Leaf& l, double collar, std::size_t j) {
  std::vector<double> lay(l.roots_.psi_soil_.size(), 0.0);
  double total = 0.0;
  l.roots_.uptake_at(collar, l.roots_.psi_soil_, lay, total);
  return lay[j];
}

void check(const std::vector<double>& psi_soil,
           const std::vector<double>& soil_depth, double collar,
           std::size_t rooted) {
  Drivers d;
  pl::Leaf l = make_leaf(d, psi_soil, soil_depth, rooted);
  l.find_root_collar_psi();

  std::vector<double> rows;
  const double total_kg = l.roots_.duptake_dpsi(collar, l.roots_.psi_soil_, rows);
  // ⚠️ ONE ENTRY PER SOIL LAYER, matching the uptake it is grafted onto. Sized by
  // the deepest ROOTED layer instead, a shallow-rooted plant hands back a shorter
  // vector than its own uptake -- which plant's gradient ladder refuses a whole
  // sweep on, and which every fixture whose layers are all rooted is blind to.
  ok(rows.size() == l.roots_.psi_soil_.size(),
     "one collar row per soil layer, not per rooted layer");
  ++compared;
  if (!std::isfinite(total_kg)) {
    return;  // a kink: duptake_dpsi refuses by contract, and rows is all NaN
  }

  double sum_mol = 0.0;
  for (std::size_t j = 0; j < rows.size(); ++j) {
    sum_mol += rows[j];
    // Two steps: the answer, and the floor that says how well it is known.
    const double h1 = 1e-5, h2 = 1e-6;
    const double d1 = (uptake_of(l, collar + h1, j) -
                       uptake_of(l, collar - h1, j)) / (2.0 * h1);
    const double d2 = (uptake_of(l, collar + h2, j) -
                       uptake_of(l, collar - h2, j)) / (2.0 * h2);
    const double floor_j = std::max(std::abs(d1 - d2), 1e-14);
    const double gap = std::abs(rows[j] - d1);
    ++compared;
    const double ratio = gap / floor_j;
    worst_ratio = std::max(worst_ratio, ratio);
    // Held against the REFERENCE's own noise rather than a constant: a tolerance
    // tighter than that is measuring the step size.
    ok(ratio < 50.0,
       "d(uptake_" + std::to_string(j) + ")/d(collar) at " +
           std::to_string(collar) + ": supplied " + std::to_string(rows[j]) +
           ", differenced " + std::to_string(d1) + ", floor " +
           std::to_string(floor_j));
  }

  // The rows are the aggregate, in the units the aggregate is reported in.
  const double from_rows = sum_mol * pl::kg_per_mol_h2o;
  const double rel = std::abs(from_rows - total_kg) /
                     std::max(std::abs(total_kg), 1e-300);
  worst_sum = std::max(worst_sum, rel);
  ok(rel < 1e-12,
     "the per-layer rows sum to the aggregate at " + std::to_string(collar));
}

}  // namespace

int main() {
  printf("the supplied per-layer collar rows, against a difference of the uptake\n");
  const std::vector<std::vector<double>> soils = {
      {2.0}, {1.0, 2.0, 3.0}, {0.5, 1.0, 1.5, 2.0, 2.5}};
  const std::vector<std::vector<double>> depths = {
      {1.0}, {0.5, 1.0, 1.5}, {0.4, 0.8, 1.2, 1.6, 2.0}};
  for (std::size_t k = 0; k < soils.size(); ++k) {
    for (double collar : {0.25, 0.75, 1.0, 1.25, 1.5, 2.5, 3.5, 4.5, 5.5}) {
      check(soils[k], depths[k], collar, soils[k].size());
      // and the same profile with roots reaching only the top layers
      if (soils[k].size() > 1) check(soils[k], depths[k], collar, 1);
    }
  }
  printf("    %d layer rows compared | worst gap %.2fx the reference's own floor"
         " | rows-to-aggregate %.3e\n", compared, worst_ratio, worst_sum);
  // A row not reached is a row not checked, and the sweep is what reaches them.
  ok(compared >= 40, "the sweep reached at least 40 layer rows");
  printf("%d checks, %d failures\n", compared + 1, failures);
  return failures == 0 ? 0 : 1;
}
