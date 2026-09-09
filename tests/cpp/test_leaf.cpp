// Pure-C++ test suite for the leaf model. No R, no test framework, no linking:
//   make -C tests/cpp && ./tests/cpp/test_leaf
//
// The trait values and drivers below are lifted from plant's
// tests/testthat/test-leaf.r so the two suites exercise the same operating
// point. The expected values here were produced BY this implementation and are
// regression guards rather than independent references -- but the implementation
// has been cross-checked against plant's compiled build and found bit-identical,
// so they are guarding a verified model.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <odelia/tangent.hpp>

#include <chrono>
#include <cmath>
#include <functional>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <vector>

namespace {

int failures = 0;
int checks = 0;

void ok(bool pass, const std::string &what) {
  ++checks;
  if (!pass) {
    ++failures;
    printf("  FAIL  %s\n", what.c_str());
  }
}

// A refusal, checked for its REASON and not only for throwing. The message is
// half the point of every `util::stop` in this package -- several of them exist to
// tell the caller which of two similar-looking configurations they are in -- so a
// bare "it threw" would pass on a refusal that had drifted onto the wrong branch.
// Returns rather than asserting, so the call site names the case.
template <typename F>
bool throws_with(F &&f, const std::string &fragment) {
  try {
    f();
  } catch (const std::exception &e) {
    return std::string(e.what()).find(fragment) != std::string::npos;
  }
  return false;
}

void near(double got, double want, double tol, const std::string &what) {
  ++checks;
  const double err = std::abs(got - want);
  const double scale = std::max(1.0, std::abs(want));
  if (!(err / scale <= tol)) {
    ++failures;
    printf("  FAIL  %s: got %.12g, want %.12g (rel err %.3g > %.3g)\n",
           what.c_str(), got, want, err / scale, tol);
  }
}

// Trait values from plant's test-leaf.r.
struct Drivers {
  double theta = 0.000157; // Huber value, m2 sapwood m-2 leaf
  double K_s = 1.0;        // stem-specific conductivity
  double h = 5.0;          // path length, m
  double PPFD = 900.0;
  double atm_vpd = 2.0;
  double ca = 40.0;
  double atm_o2_kpa = 21.0;
  double leaf_temp = 25.0;
  double atm_kpa = 101.3;
  double area_leaf = 0.05;
};

// set_traits' fourteenth argument, R_d_25: dark respiration at 25 C. The class
// default, so a call that is not about respiration can pass it and change nothing.
const double kRd25 = 1.44;
// JS22 s hydraulic unit cost. Inert on every curve these tests exercise -- it is
// read only by the JS22 arm -- but set_traits is positional, so it is named here
// rather than repeated as a literal at a dozen call sites.
const double kGammaJS22 = 1.0;
const double kCMaxA = 0.6;
const double kCMaxB = 0.0;
// The TF24_floor curve's ONLY parameter: its price of water as transpiration goes to
// zero. The hydraulic half is TF24's cost at TF24's own traits, so there is nothing
// else to name. Mirrors kLambdaCF77 in test_golden.cpp and carries the same caveat:
// its scale is set by the leaf, whose own marginal cost of water runs 9e4 to 3e5 at
// these drivers.
const double kFloorLambdaO = 1.5e5;

phylloptim::Leaf make_leaf(const Drivers &d, std::vector<double> psi_soil,
                     std::vector<double> soil_depth) {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  // root carbon per unit leaf area: the old absolute carbon divided by area_leaf
  std::vector<double> mass_root_prop(psi_soil.size(),
                                     1.0 / double(psi_soil.size()) / d.area_leaf);
  l.set_physiology(fixture::root_network(mass_root_prop, soil_depth), d.PPFD, psi_soil, soil_depth,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// ---------------------------------------------------------------------------

void test_defaults_are_unset() {
  printf("defaults are unset until set_physiology\n");
  phylloptim::Leaf l;
  ok(!std::isfinite(l.ci_), "ci_ starts unset");
  ok(!std::isfinite(l.assim_colimited_), "assim_colimited_ starts unset");
  ok(!std::isfinite(l.opt_psi_stem_), "opt_psi_stem_ starts unset");
  ok(l.roots_.psi_soil_.empty(), "psi_soil_ starts empty");
  ok(l.use_energy_balance_ == false, "energy balance defaults off");
}

void test_vulnerability_curve() {
  printf("xylem vulnerability curve\n");
  phylloptim::Leaf l;
  near(l.proportion_of_conductivity(0.0), 1.0, 1e-12,
       "full conductivity at zero potential");
  // psi_crit is NOT the 1%-conductivity point -- with the default b and c it sits
  // at ~5%. (The 1% point is where setup_transpiration ends its knot grid, which
  // is a different and larger potential.) Asserting the closed form rather than a
  // round number so this stays honest if the defaults change.
  near(l.proportion_of_conductivity(l.psi_crit),
       std::exp(-std::pow(l.psi_crit / l.stem_b, l.stem_c)), 1e-12,
       "conductivity at psi_crit matches the Weibull closed form");
  ok(l.proportion_of_conductivity(l.psi_crit) < 0.06,
     "conductivity at psi_crit is a few percent");
  ok(l.proportion_of_conductivity(1.0) > l.proportion_of_conductivity(2.0),
     "conductivity declines monotonically");
}

// The analytic trait partials of the cumulative vulnerability integral, against
// central differences of the closed form they differentiate. No new reference
// values: the function being differenced is the one the knot builder already
// seeds every knot from.
//
// ⚠️ The tolerance below is the CENTRAL DIFFERENCE's accuracy, not the analytic
// form's. Measured on this grid, the worst disagreement walks as h^2 -- 2.9e-5,
// 2.9e-7, 3.1e-9 at relative steps 1e-3, 1e-4, 1e-5 -- and then back UP to 2.9e-7
// at 1e-7 as cancellation takes over. A clean quadratic decay is what says the
// residual belongs to the difference and not to the formula, so tightening this
// without shrinking h would only be pinning the difference's own error.
//
// b and c are deliberately neutral here, and the grid brackets both curves rather
// than either: vulnerability.hpp is shared by the stem (stem_b, stem_c) and the
// root (root_b, root_c), whose defaults are both 3.898245 / 2.680147 (hazard 1).
void test_vulnerability_integral_derivatives() {
  printf("analytic trait derivatives of the vulnerability integral\n");
  using phylloptim::cumulative_vulnerability_integral_at;
  using phylloptim::cumulative_vulnerability_integral_derivatives_at;

  // (psi/b)^c past this is past the curve's last knot, so it is unreachable and
  // the series is not asked to hold there. It bounds the grid instead of the grid
  // bounding it -- for c = 12, psi/b = 8 means (psi/b)^c = 7e10.
  const double x_max = phylloptim::vulnerability_x_max();
  const double rel_step = 1e-5;

  double worst_value = 0.0, worst_dpsi = 0.0, worst_db = 0.0, worst_dc = 0.0;
  int points = 0;
  for (double b : {0.5, 2.0, 7.0}) {
    for (double c = 0.4; c <= 12.001; c += 0.1) {
      for (double u = 0.075; u <= 8.001; u *= 1.1) {
        if (pow(u, c) > x_max) {
          continue;
        }
        ++points;
        const double psi = u * b;
        const auto d = cumulative_vulnerability_integral_derivatives_at(psi, b, c);

        // The series' own value against boost's. Not the point of the test, but
        // it is the cheapest check that the series is the right series.
        const double scale = std::max(1.0, std::abs(d.value));
        worst_value = std::max(
            worst_value,
            std::abs(d.value - cumulative_vulnerability_integral_at(psi, b, c)) /
                scale);

        const double hb = rel_step * b, hc = rel_step * c, hp = rel_step * psi;
        const double fd_psi =
            (cumulative_vulnerability_integral_at(psi + hp, b, c) -
             cumulative_vulnerability_integral_at(psi - hp, b, c)) /
            (2.0 * hp);
        const double fd_b =
            (cumulative_vulnerability_integral_at(psi, b + hb, c) -
             cumulative_vulnerability_integral_at(psi, b - hb, c)) /
            (2.0 * hb);
        const double fd_c =
            (cumulative_vulnerability_integral_at(psi, b, c + hc) -
             cumulative_vulnerability_integral_at(psi, b, c - hc)) /
            (2.0 * hc);
        // Relative to the derivative, floored: dG/dc passes through zero on this
        // grid, and a relative error against a vanishing denominator says nothing.
        worst_dpsi = std::max(worst_dpsi, std::abs(d.dpsi - fd_psi) /
                                              std::max(1e-2, std::abs(fd_psi)));
        worst_db = std::max(worst_db, std::abs(d.db - fd_b) /
                                          std::max(1e-2, std::abs(fd_b)));
        worst_dc = std::max(worst_dc, std::abs(d.dc - fd_c) /
                                          std::max(1e-2, std::abs(fd_c)));
      }
    }
  }
  printf("    %d points | value %.3g | dpsi %.3g | db %.3g | dc %.3g\n", points,
         worst_value, worst_dpsi, worst_db, worst_dc);
  ok(points > 10000, "the grid inside the curve's domain is not empty");
  near(worst_value, 0.0, 1e-14, "the series value matches the closed form");
  near(worst_dpsi, 0.0, 1e-7, "dG/dpsi matches a central difference");
  near(worst_db, 0.0, 1e-7, "dG/db matches a central difference");
  near(worst_dc, 0.0, 1e-7, "dG/dc matches a central difference");

  // dG/dpsi is the vulnerability curve itself, which is worth stating once
  // directly rather than only through a difference.
  const auto at_default =
      cumulative_vulnerability_integral_derivatives_at(3.0, 3.898245, 2.680147);
  near(at_default.dpsi, std::exp(-std::pow(3.0 / 3.898245, 2.680147)), 1e-15,
       "dG/dpsi is exp(-(psi/b)^c)");

  // G is identically zero at psi = 0 for every b and c, so both trait partials
  // are zero there and neither pow(0, c) nor log(0) is evaluated.
  const auto at_zero =
      cumulative_vulnerability_integral_derivatives_at(0.0, 3.898245, 2.680147);
  ok(at_zero.value == 0.0, "G is exactly zero at psi = 0");
  ok(at_zero.dpsi == 1.0, "dG/dpsi is exactly 1 at psi = 0");
  ok(at_zero.db == 0.0, "dG/db is exactly zero at psi = 0");
  ok(at_zero.dc == 0.0, "dG/dc is exactly zero at psi = 0");

  // Past the last knot the series would overflow rather than mislead, so the
  // bound is asserted instead of branched around.
  bool threw = false;
  try {
    cumulative_vulnerability_integral_derivatives_at(3.0, 1.0, 2.0);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(threw, "an argument past the curve's domain throws");
  double just_inside = 0.0;
  try {
    just_inside =
        cumulative_vulnerability_integral_derivatives_at(
            phylloptim::vulnerability_psi_max(3.898245, 2.680147), 3.898245,
            2.680147)
            .value;
  } catch (const std::runtime_error &) {
    just_inside = -1.0;
  }
  ok(just_inside > 0.0, "the last knot itself is inside the bound");
}

void test_spline_matches_direct_integration() {
  printf("pre-integrated spline vs direct quadrature\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  // This is the check plant makes at test-leaf.r:214. The spline is what the hot
  // path reads; adaptive Simpson integrates the curve directly.
  for (double psi_stem : {2.5, 3.0, 4.0, 5.0}) {
    near(l.transpiration(psi_stem, 2.0),
         l.transpiration_full_integration(psi_stem, 2.0), 1e-6,
         "transpiration at psi_stem=" + std::to_string(psi_stem));
  }
}

void test_arrhenius() {
  printf("temperature response\n");
  phylloptim::Leaf l;
  near(l.arrh_curve(phylloptim::vcmax_ha, 100.0, 25.0), 100.0, 1e-12,
       "Arrhenius is the identity at the 25 C reference");
  ok(l.arrh_curve(phylloptim::vcmax_ha, 100.0, 35.0) > 100.0,
     "Arrhenius rises above the reference temperature");
  ok(l.peak_arrh_curve(phylloptim::jmax_ha, 100.0, 60.0, phylloptim::jmax_H_d,
                       phylloptim::jmax_d_S) <
         l.peak_arrh_curve(phylloptim::jmax_ha, 100.0, 30.0, phylloptim::jmax_H_d,
                           phylloptim::jmax_d_S),
     "peaked Arrhenius declines past its optimum");
}

void test_saturation_vapour_pressure() {
  printf("saturation vapour pressure\n");
  phylloptim::Leaf l;
  // Tetens at 25 C is ~3.167 kPa.
  near(l.saturation_vapour_pressure(25.0), 3.167, 1e-3, "es(25 C)");
  // Delta ~ 0.189 kPa/K at 25 C.
  near(l.saturation_vapour_pressure_slope(25.0), 0.189, 5e-3, "Delta(25 C)");
  ok(l.saturation_vapour_pressure(30.0) > l.saturation_vapour_pressure(20.0),
     "es increases with temperature");
}

void test_solve_single_layer() {
  printf("root-collar solve, single soil layer\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  ok(std::isfinite(l.opt_psi_stem_), "psi_stem is finite");
  ok(std::isfinite(l.profit_), "profit is finite");
  ok(l.opt_psi_stem_ > 2.0, "stem is drier than the soil");
  ok(l.opt_psi_stem_ <= l.psi_crit, "stem stays within psi_crit");
  ok(l.transpiration_ > 0.0, "transpiration is positive");
  ok(l.assim_colimited_ > 0.0, "assimilation is positive");
  // Regression guards -- see the note at the top of this file. Regenerate these
  // deliberately, alongside the golden file, and state the movement in the PR.
  near(l.opt_psi_stem_, 3.595332, 1e-5, "opt_psi_stem_");
  near(l.assim_colimited_, 5.601016, 1e-5, "assim_colimited_");
  near(l.transpiration_, 1.141941e-05, 1e-5, "transpiration_");
  near(l.profit_, 2.517157, 1e-5, "profit_");
}

void test_solve_is_deterministic() {
  printf("repeated solves are bit-identical\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double first = l.profit_, psi = l.opt_psi_stem_;
  for (int i = 0; i < 50; ++i) {
    l.find_root_collar_psi();
  }
  ok(l.profit_ == first, "profit is unchanged after 50 re-solves");
  ok(l.opt_psi_stem_ == psi, "psi_stem is unchanged after 50 re-solves");
}

void test_drier_soil_costs_carbon() {
  printf("response to drying soil\n");
  Drivers d;
  double prev_profit = 1e9, prev_E = 1e9;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi}, {1.0});
    l.find_root_collar_psi();
    ok(l.profit_ <= prev_profit,
       "profit does not rise as soil dries to " + std::to_string(psi));
    ok(l.transpiration_ <= prev_E,
       "transpiration does not rise as soil dries to " + std::to_string(psi));
    prev_profit = l.profit_;
    prev_E = l.transpiration_;
  }
}

void test_light_response() {
  printf("light response\n");
  Drivers dim = Drivers(), bright = Drivers();
  dim.PPFD = 200;
  bright.PPFD = 1800;
  phylloptim::Leaf shaded = make_leaf(dim, {2.0}, {1.0});
  phylloptim::Leaf sunlit = make_leaf(bright, {2.0}, {1.0});
  shaded.find_root_collar_psi();
  sunlit.find_root_collar_psi();
  ok(sunlit.assim_colimited_ > shaded.assim_colimited_,
     "brighter light assimilates more");
  ok(sunlit.transpiration_ > shaded.transpiration_, "brighter light transpires more");
}

void test_multi_layer_soil() {
  printf("multi-layer soil\n");
  Drivers d;
  // ⚠️ The profile is CUMULATIVE depth to the bottom of each layer, so three
  // 0.5 m layers are {0.5, 1.0, 1.5}. This read {0.5, 0.5, 0.5} until #626 --
  // three layers all bottoming at 0.5 m, which put two of the three midpoints at
  // the same depth and gave the root network a thickness of 0.5/3. Nothing
  // noticed, because the assertions below are finiteness and sign only. The
  // strictly-increasing guard in layer_thicknesses() is what surfaced it.
  phylloptim::Leaf l = make_leaf(d, {1.0, 2.0, 3.0}, {0.5, 1.0, 1.5});
  l.find_root_collar_psi();
  ok(std::isfinite(l.profit_), "profit is finite with three layers");
  ok(l.soil_consumption_.size() == 3u, "one consumption term per layer");
  double total = 0.0;
  for (double s : l.soil_consumption_) {
    ok(std::isfinite(s), "per-layer consumption is finite");
    total += s;
  }
  ok(total > 0.0, "total soil water uptake is positive");
}

void test_shutdown_when_soil_is_drier_than_psi_crit() {
  printf("shutdown past psi_crit\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {12.0}, {1.0});
  l.find_root_collar_psi();
  ok(std::isfinite(l.profit_), "profit stays finite past psi_crit");
  ok(l.profit_ <= 0.0, "profit is non-positive when shut down");
  near(l.opt_psi_stem_, l.psi_crit, 1e-12, "stem is held at psi_crit");
}

// The shutdown path used to leak the previous solve's fluxes (plant #578/#577).
// set_shutdown_state now writes them, so this asserts the fixed behaviour: a
// shut-down leaf moves no water, respires at R_d, and sits at the CO2
// compensation point -- and, critically, none of that depends on what ran before.
void test_shutdown_writes_its_own_fluxes() {
  printf("shutdown writes its own fluxes (plant #578 fixed)\n");
  Drivers d;
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> mrp{1.0 / d.area_leaf}, depth{1.0};
  const auto solve = [&](double psi) {
    std::vector<double> ps{psi};
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, ps, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  solve(4.0); // wet enough to transpire
  ok(l.transpiration_ > 0.0, "the wet solve transpires");

  solve(20.0); // far drier than psi_crit: the leaf shuts down
  ok(l.profit_ < 0.0, "the dry solve is a shutdown (profit < 0)");
  near(l.transpiration_, 0.0, 1e-300, "transpiration is zero, not stale");
  near(l.stom_cond_CO2_, 0.0, 1e-300, "conductance is zero, not stale");
  near(l.E_up_, 0.0, 1e-300, "soil uptake is zero, not stale");
  for (double s : l.soil_consumption_) {
    near(s, 0.0, 1e-300, "per-layer consumption is zero, not stale");
  }
  // Respiring, not simply idle: profit_ is -R_d_ - hydraulic_cost, so the
  // consistent assimilation is -R_d_. Zero would be inconsistent with profit_.
  ok(l.assim_colimited_ < 0.0, "assimilation is negative (respiring)");
  near(l.assim_colimited_, l.profit_ + l.hydraulic_cost_TF(l.psi_crit), 1e-12,
       "assimilation is consistent with profit and the hydraulic cost");
  ok(std::isfinite(l.ci_), "ci is set to the compensation point, not left stale");

  // Order independence is the property that was actually broken: a fresh leaf
  // taken straight to dry must report exactly the same thing.
  phylloptim::Leaf fresh = make_leaf(d, {20.0}, {1.0});
  fresh.find_root_collar_psi();
  ok(fresh.transpiration_ == l.transpiration_,
     "a fresh leaf gives the same transpiration");
  ok(fresh.assim_colimited_ == l.assim_colimited_,
     "a fresh leaf gives the same assimilation");
  ok(fresh.profit_ == l.profit_, "a fresh leaf gives the same profit");
}

// soil_consumption_ used to be cleared with .resize, whose fill reaches only
// newly-added elements, while the uptake loop writes only up to max_soil_layer --
// the deepest layer holding root carbon. So a shallow-rooted plant solving on a
// Leaf that a deep-rooted one used before it inherited the deep layers' uptake
// (plant #577, fixed in plant by #585). The soil layer count is the same in both
// solves here; only the rooted depth shrinks, which is why resizing never noticed.
void test_shallow_roots_do_not_inherit_deep_uptake() {
  printf("shallow roots do not inherit the previous plant's deep uptake\n");
  Drivers d;
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);

  const std::vector<double> psi_soil{1.0, 1.5, 2.0};
  const std::vector<double> depth{1.0, 2.0, 3.0};
  const auto solve = [&](std::vector<double> root_carbon) {
    l.set_physiology(fixture::root_network(root_carbon, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  // A tree rooted through all three layers, then a seedling rooted in the top
  // layer only -- trailing zeros in the carbon vector is how plant expresses
  // "not rooted that deep".
  const double c = 1.0 / 3.0 / d.area_leaf;
  solve({c, c, c});
  ok(l.soil_consumption_[1] != 0.0 && l.soil_consumption_[2] != 0.0,
     "the deep-rooted solve draws from layers 2 and 3");

  solve({1.0 / d.area_leaf, 0.0, 0.0});
  ok(l.soil_consumption_.size() == 3u,
     "the consumption vector still spans every soil layer");
  near(l.soil_consumption_[1], 0.0, 1e-300,
       "layer 2 is zero, not the tree's uptake");
  near(l.soil_consumption_[2], 0.0, 1e-300,
       "layer 3 is zero, not the tree's uptake");

  // Order independence is the property that matters: plant reuses one Leaf for
  // every individual in a patch, so the seedling must not depend on its neighbour.
  phylloptim::Leaf fresh;
  fresh.setup_transpiration(100);
  fresh.setup_root_vulnerability(100);
  fresh.set_physiology(fixture::root_network({1.0 / d.area_leaf, 0.0, 0.0}, depth), d.PPFD, psi_soil, depth,
                       d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                       d.atm_o2_kpa, d.atm_kpa);
  fresh.find_root_collar_psi();
  for (size_t i = 0; i < 3; ++i) {
    ok(fresh.soil_consumption_[i] == l.soil_consumption_[i],
       "a fresh leaf gives the same per-layer consumption");
  }
}

// The other early exit that determines the operating point without going through
// profit_psi_stem_TF: assimilation is negative even at ci = ca, so there is no
// light level at which opening the stomata pays. It set profit_ but left the
// leaf-side rates alone, so they held whatever the previous solve wrote -- and on
// a fresh leaf, nothing at all. Ported from plant develop (#585).
void test_negative_assim_exit_writes_its_own_rates() {
  printf("the assim_max_ < 0 exit writes its own rates\n");
  Drivers d;
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  const std::vector<double> psi_soil{1.0}, depth{1.0};
  const std::vector<double> root{1.0 / d.area_leaf};
  const auto solve = [&](double ppfd) {
    l.set_physiology(fixture::root_network(root, depth), ppfd, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  solve(900.0); // bright: a normal optimising solve
  ok(l.transpiration_ > 0.0, "the bright solve transpires");

  // Dim enough that gross assimilation cannot cover R_d even with ci at ca. The
  // soil is wet, so this is not the psi_crit shut-down path -- it is the
  // assim_max_ < 0 exit, which parks the stem in equilibrium with the collar.
  solve(10.0);
  ok(l.assim_max_ < 0.0, "the dim solve takes the assim_max_ < 0 exit");
  near(l.transpiration_, 0.0, 1e-300, "transpiration is zero, not stale");
  near(l.stom_cond_CO2_, 0.0, 1e-300, "conductance is zero, not stale");
  near(l.assim_colimited_, -l.R_d_, 1e-12,
       "net assimilation is -R_d, not the bright solve's");
  // The invariant the fix buys: profit_ == assim_colimited_ - hydraulic cost in
  // every branch. It held numerically to the last bit when measured.
  near(l.assim_colimited_ - l.hydraulic_cost_TF(l.opt_root_psi_), l.profit_,
       1e-14, "profit is consistent with assimilation and the hydraulic cost");

  phylloptim::Leaf fresh;
  fresh.setup_transpiration(100);
  fresh.setup_root_vulnerability(100);
  fresh.set_physiology(fixture::root_network(root, depth), 10.0, psi_soil, depth, d.K_s * d.theta / d.h,
                       d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  fresh.find_root_collar_psi();
  ok(fresh.transpiration_ == l.transpiration_,
     "a fresh leaf gives the same transpiration");
  ok(fresh.assim_colimited_ == l.assim_colimited_,
     "a fresh leaf gives the same assimilation");
  ok(fresh.profit_ == l.profit_, "a fresh leaf gives the same profit");
}

void test_analytic_gradient_matches_finite_difference() {
  printf("analytic dprofit/dpsi_collar vs central difference\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double p0 = l.opt_root_psi_;
  const double target = std::max(2.2, std::min(p0, l.psi_crit - 0.5));
  const double eps = 1e-5;
  const double analytic = l.dprofit_droot_collar_psi(target);
  const double up = l.evaluate_root_collar_psi(target + eps);
  const double dn = l.evaluate_root_collar_psi(target - eps);
  const double fd = (up - dn) / (2 * eps);
  ok(std::isfinite(analytic), "analytic gradient is finite");
  near(analytic, fd, 2e-3, "analytic gradient matches central difference");
}

// lambda = dA/dE is the first-order condition every model in this family shares,
// so checking the analytic lambda against a finite-difference dA/dE is a check on
// the whole optimisation, not just on one formula.
// dprofit_droot_collar_psi reads the supply path's signed soil potentials, which
// used to be seated only by a solve -- so calling it on a leaf that had had
// set_physiology but not find_root_collar_psi read an empty vector. It now seats
// them itself. Ported from plant develop (#585).
void test_gradient_needs_no_prior_solve() {
  printf("dprofit/dpsi_collar does not require a prior solve\n");
  Drivers d;
  phylloptim::Leaf solved = make_leaf(d, {2.0}, {1.0});
  solved.find_root_collar_psi();
  phylloptim::Leaf unsolved = make_leaf(d, {2.0}, {1.0});
  // Bit-identical, not merely close: seating the potentials from psi_soil_ is
  // exactly what a solve does, so this is idempotent and moves no arithmetic.
  ok(unsolved.dprofit_droot_collar_psi(2.5) == solved.dprofit_droot_collar_psi(2.5),
     "the gradient is the same with and without a prior solve");
}

// The reversed-gradient state: the collar is asked about a potential drier than
// the stem it would have to supply, so there is no flow and no informative
// gradient. psi_stem_to_ci does not return non-finite there -- it either throws
// (gc goes negative, the residual stops crossing zero, and the bracketing solver
// gives up) or returns a number built on a negative conductance. Both were
// reachable from TF24f's acclimation gradient on a dry patch.
void test_gradient_is_zero_in_reversed_gradient_state() {
  printf("dprofit/dpsi_collar returns zero when the gradient reverses\n");
  Drivers d;
  // Drier than psi_crit at every layer, so the leaf is shut down and every collar
  // potential below the wettest layer implies psi_stem < psi_upstream.
  phylloptim::Leaf l = make_leaf(d, {5.9, 6.15, 6.4, 6.65, 6.9}, {1.0, 2.0, 3.0, 4.0, 5.0});
  l.find_root_collar_psi();
  for (double target : {1.0, 3.0, 5.5, 5.9}) {
    const double psi_stem = l.find_psi_stem_from_psi_root(
        target, l.roots_.psi_soil_);
    ok(target >= psi_stem,
       "the target is in the reversed-gradient state at " + std::to_string(target));
    ok(l.dprofit_droot_collar_psi(target) == 0.0,
       "the gradient is exactly zero at " + std::to_string(target));
  }
}

// The 0.0 the two exits above return is a SENTINEL, and the collar root-find is what makes
// telling it apart from a genuine stationary point load bearing: it proposes
// root-finding on dprofit == 0, and the sentinel fires at the WET END of the very
// bracket such a solve would search. At bound_a = root_zero_E uptake is zero by
// construction, so psi >= psi_stem and the reversed-gradient exit is taken -- and
// a bracketing solver evaluates that endpoint FIRST, to check its bracket
// brackets. So it would read the sentinel as the answer and return the
// zero-transpiration point as the optimum.
//
// This test pins the distinction rather than the eventual solver, so it is
// independent of how 11a is implemented. Note what it can and cannot see: the
// golden file cannot see any of this (no golden column comes from the gradient),
// and neither can a test that only samples interior points -- the sentinel region
// measured at most 3.46e-07 MPa wide over the golden grid, which is exactly why
// it would survive casual testing.
void test_gradient_reports_feasibility() {
  printf("dprofit/dpsi_collar distinguishes its 0.0 sentinel from a stationary point\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});

  // The bracket a root-find would search, built exactly as prepare_collar_solve
  // does: the wet end is where uptake is zero, the dry end whichever limit binds.
  const double wettest = 2.0;
  const double root_zero_E = l.find_root_psi(wettest, l.roots_.psi_soil_, 0);
  const double root_crit = l.find_root_psi(wettest, l.roots_.psi_soil_, 1);
  const double bound_b = std::min(root_crit, l.roots_.root_psi_crit);
  ok(root_zero_E < bound_b, "the bracket is non-empty");

  // The wet endpoint: 0.0, and NOT a stationary point.
  bool feasible = true;
  const double at_a = l.dprofit_droot_collar_psi(root_zero_E, &feasible);
  ok(at_a == 0.0, "the gradient is 0.0 at the wet bracket endpoint");
  ok(!feasible, "and the endpoint is reported infeasible, so the 0.0 is a sentinel");

  // Why accepting it would be wrong, stated as a number rather than an assertion
  // about intent: the endpoint is a far worse operating point than the optimum.
  l.find_root_collar_psi();
  const double best = l.profit_;
  const double at_endpoint = l.evaluate_root_collar_psi(root_zero_E);
  ok(at_endpoint < best - 1.0,
     "and profit at the endpoint is much worse than at the optimum");

  // An interior point: a real derivative, so a zero there WOULD be stationary.
  feasible = false;
  const double mid = 0.5 * (root_zero_E + bound_b);
  const double at_mid = l.dprofit_droot_collar_psi(mid, &feasible);
  ok(feasible, "an interior point is reported feasible");
  ok(std::isfinite(at_mid) && at_mid != 0.0,
     "and its gradient is a finite non-zero number");

  // The reversed-gradient state is infeasible too, and there the 0.0 is the whole
  // answer -- this is the case test_gradient_is_zero_in_reversed_gradient_state
  // already pins, checked here for the flag rather than the value.
  phylloptim::Leaf dry = make_leaf(d, {5.9, 6.15, 6.4, 6.65, 6.9},
                             {1.0, 2.0, 3.0, 4.0, 5.0});
  dry.find_root_collar_psi();
  feasible = true;
  ok(dry.dprofit_droot_collar_psi(3.0, &feasible) == 0.0 && !feasible,
     "the reversed-gradient state is reported infeasible");

  // The out-parameter changes nothing about the value, so TF24f's contract (the
  // return value alone, consumed as an ODE rate) is untouched. Bit-identical, not
  // merely close: the flag is a write to a caller's bool, not arithmetic.
  ok(l.dprofit_droot_collar_psi(mid) == at_mid,
     "passing no flag returns exactly the same value");
  ok(l.dprofit_droot_collar_psi(root_zero_E) == 0.0,
     "and the sentinel is still 0.0 for a caller that does not ask");
}

// The AD derivative and the forward model are instantiations of ONE
// body, so they cannot be derivatives of different functions. That was not true
// before: `detail::assim_colimited_ad` associated the electron-limited term
// left-to-right where `assim_electron_limited` divides the bracket first, and used
// `s*s` where `assim_colimited` uses `pow(s, 2)`.
//
// The check that bites is the DERIVATIVE against a central difference of the
// `double` function, because that is precisely the identity the drift broke. Note a
// weaker test would have passed throughout: the drifted replica was still a
// perfectly good derivative of *itself*, and agreed with the real function to
// ~1e-16 in VALUE. It is the derivative-of-the-same-function property that failed.
void test_ad_kernels_are_the_model_not_a_mirror() {
  printf("AD differentiates the model's own algebra, not a mirror of it\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  using AD = xad::fwd<double>::active_type;

  // 1. The kernels and the double entry points are the same code, so they must
  // agree BIT-EXACTLY, not merely closely. This is what "one body" means.
  for (double ci : {5.0, 12.0, 20.0, 29.26, 35.0, 39.0}) {
    ok(l.assim_colimited_kernel(ci) == l.assim_colimited(ci),
       "assim_colimited_kernel<double> is assim_colimited at ci=" +
           std::to_string(ci));
    ok(l.assim_rubisco_limited_kernel(ci) == l.assim_rubisco_limited(ci),
       "rubisco kernel matches bit-exactly at ci=" + std::to_string(ci));
    ok(l.assim_electron_limited_kernel(ci) == l.assim_electron_limited(ci),
       "electron kernel matches bit-exactly at ci=" + std::to_string(ci));
  }
  for (double p : {0.5, 2.0, 3.5949, 5.0}) {
    ok(l.hydraulic_cost_TF_kernel(p) == l.hydraulic_cost_TF(p),
       "hydraulic_cost_TF_kernel<double> is hydraulic_cost_TF at psi=" +
           std::to_string(p));
  }

  // 2. The AD derivative of the kernel against a central difference of the DOUBLE
  // function. Richardson-extrapolated, so the FD reference is good to ~1e-10 and
  // the tolerance is testing the derivative rather than the difference quotient.
  const auto richardson = [](auto f, double x, double h) {
    const double d1 = (f(x + h) - f(x - h)) / (2 * h);
    const double d2 = (f(x + h / 2) - f(x - h / 2)) / h;
    return (4 * d2 - d1) / 3;
  };
  for (double ci : {12.0, 20.0, 29.26, 35.0}) {
    AD a = ci; xad::derivative(a) = 1.0;
    const double ad = xad::derivative(l.assim_colimited_kernel(a));
    const double fd =
        richardson([&](double x) { return l.assim_colimited(x); }, ci, 1e-4);
    near(ad, fd, 1e-8, "dA/dci: AD vs Richardson FD at ci=" + std::to_string(ci));
  }
  for (double p : {0.5, 2.0, 3.5949, 5.0}) {
    AD a = p; xad::derivative(a) = 1.0;
    const double ad = xad::derivative(l.hydraulic_cost_TF_kernel(a));
    const double fd =
        richardson([&](double x) { return l.hydraulic_cost_TF_kernel(x); }, p, 1e-4);
    near(ad, fd, 1e-8, "dcost/dpsi: AD vs Richardson FD at psi=" +
                           std::to_string(p));
  }

  // 3. The cost kernel must be PURE -- an AD probe must not scribble the cached
  // hydraulic_cost_, or a gradient evaluation would corrupt reported model state.
  const double cached = l.hydraulic_cost_TF(3.0);
  AD probe = 4.5; xad::derivative(probe) = 1.0;
  (void)l.hydraulic_cost_TF_kernel(probe);
  ok(l.hydraulic_cost_ == cached,
     "an AD probe of the cost kernel leaves hydraulic_cost_ untouched");
}

// The collar solve solves its own first-order condition, so the
// check that bites is the RESIDUAL at the returned point, not the returned value.
// Measured over the golden grid: 240 of 240 feasible rows improved their residual
// and none got worse, interior rows landing at a median |dprofit| of 5.6e-15
// against golden section's 7.8e-4.
//
// Deliberately NOT checked against profit. Profit is the wrong instrument here for
// a reason worth recording: it is the maximum, so it is flat, and its own
// numerical floor is set by the nested ci root-find's 1e-7 tolerance. Over the
// grid two rows come out ~6e-7 LOWER in profit than golden section while their
// residual improves by ten orders of magnitude -- that is the floor, not a
// regression, and a test asserting "profit never decreases" would encode the noise.
void test_collar_solve_satisfies_its_own_first_order_condition() {
  printf("the collar solve lands where dprofit == 0\n");
  Drivers d;
  // An interior optimum: the gradient at the answer should be at solver precision.
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double residual = l.dprofit_droot_collar_psi(l.opt_root_psi_);
  ok(std::abs(residual) < 1e-9,
     "the interior optimum satisfies dprofit == 0 to solver precision");
  // The bar this clears, stated as the thing it replaced: golden section left
  // -6.22e-4 here, which is what made the argmax a staircase in the traits.
  ok(std::abs(residual) < 1e-6 * 6.22e-4,
     "and it is orders below the residual golden section left behind");

  // The answer stays strictly inside the feasible interval. This is the guard on
  // the #31 regression that returning bound_a caused: below psi_upstream the
  // profit algebra runs on a negative conductance and profit is DISCONTINUOUS, so
  // an endpoint is not merely a poor answer, it is a different function. Returning
  // bound_a moved 22 golden rows' profit DOWN by 1.44 before this was caught.
  const double root_zero_E = l.find_root_psi(2.0, l.roots_.psi_soil_, 0);
  const double bound_b =
      std::min(l.find_root_psi(2.0, l.roots_.psi_soil_, 1), l.roots_.root_psi_crit);
  ok(l.opt_root_psi_ > root_zero_E,
     "the collar is strictly drier than the zero-uptake bound (#31)");
  ok(l.opt_root_psi_ <= bound_b, "and no drier than the binding dry limit");
}

// A CONSTRAINED optimum, which 42 of the 240 feasible golden-grid rows are: profit
// is still climbing at one end of the feasible interval, so dprofit never crosses
// zero inside it. The residual is then NOT small, and that is correct rather than a
// convergence failure -- so this pins the property that actually holds there.
void test_collar_solve_handles_a_pinned_optimum() {
  printf("a collar optimum pinned to its constraint\n");
  Drivers d;
  // psi_soil 4.0 over 5 layers at vpd 2.0 -- one of the measured pinned rows.
  phylloptim::Leaf l = make_leaf(d, {4.0, 4.25, 4.5, 4.75, 5.0},
                           {1.0, 2.0, 3.0, 4.0, 5.0});
  l.find_root_collar_psi();
  const double root_zero_E = l.find_root_psi(4.0, l.roots_.psi_soil_, 0);
  const double bound_b =
      std::min(l.find_root_psi(4.0, l.roots_.psi_soil_, 1), l.roots_.root_psi_crit);
  ok(std::isfinite(l.profit_), "the pinned row still yields a finite profit");
  ok(l.opt_root_psi_ > root_zero_E,
     "and a collar strictly inside the zero-uptake bound, not on it (#31)");
  ok(l.opt_root_psi_ <= bound_b, "and inside the dry bound");
  // It is pinned NEAR the wet end, but not AT it: the true argmax sits essentially
  // on the feasibility boundary, and the step-in resolves it to ~1e-6 of the
  // bracket width. Measured 3.29e-4 wetter than golden section's answer, with a
  // profit 1.27e-3 HIGHER -- so being pinned is not being wrong.
  ok(l.opt_root_psi_ - root_zero_E < 1e-4 * (bound_b - root_zero_E),
     "the pinned answer sits at the wet end of the interval");
  // The gradient there is genuinely non-zero, which is what "constrained" means.
  ok(l.dprofit_droot_collar_psi(l.opt_root_psi_) < 0.0,
     "profit is decreasing at the pinned answer, so the bound is what binds");
  ok(l.operating_point_kind() ==
         phylloptim::Leaf::OperatingPointKind::BoundarySoil,
     "and the leaf says so: the point is tagged boundary-soil");
}

// The operating-point classification, on ONE REUSED LEAF -- which is the only way
// to test it. `Leaf` is a value member that plant drives every individual in a
// patch through (hazard 8), so the failure this guards is a branch that declines
// to write the tag and thereby reports the PREVIOUS plant's kind of operating
// point. test_golden cannot see that: it builds a fresh Leaf per grid point, so a
// bit-identical golden run says nothing here. The sequence below alternates kinds
// on purpose, so an unwritten tag reads as the step before it.
//
// The count check over the whole golden grid lives in test_golden.cpp; this is
// the other half -- that every path writes, and that the tag comes from the
// branch rather than from the numbers.
void test_operating_point_kind_is_written_by_every_path() {
  printf("the collar solve says which kind of operating point it found\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  ok(l.operating_point_kind() == Kind::Unsolved,
     "a fresh Leaf reports no operating point");

  const std::vector<double> depth{1.0}, root{1.0 / d.area_leaf};
  const auto solve = [&](double psi, double ppfd) {
    std::vector<double> ps{psi};
    l.set_physiology(fixture::root_network(root, depth), ppfd, ps, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  solve(2.0, d.PPFD);
  ok(l.operating_point_kind() == Kind::Interior,
     "wet soil gives an interior profit maximum");

  // Drier than psi_crit: shutdown on water. Would read `interior` if the shutdown
  // exit did not write.
  solve(20.0, d.PPFD);
  ok(l.operating_point_kind() == Kind::HydraulicShutdown,
     "past psi_crit the point is a hydraulic shutdown");

  // ⚠️ THE ASSERTION THE TAG EXISTS FOR. dprofit's shut-down exit returns a hard
  // 0.0 sentinel, so a classifier of the form "|residual| < tol => interior
  // optimum" would call this an interior stationary point, and then read a
  // curvature of zero off the same sentinel and agree with itself. The residual
  // really is exactly 0.0 here, and the tag really does say shutdown.
  bool feasible = true;
  const double residual = l.dprofit_droot_collar_psi(l.opt_root_psi_, &feasible);
  ok(residual == 0.0, "dprofit at a shut-down point is exactly the 0.0 sentinel");
  ok(!feasible, "and reports itself infeasible");
  ok(l.operating_point_kind() != Kind::Interior,
     "so a residual test would misclassify it, and the tag does not");

  // Back to wet: would read `hydraulic-shutdown` if the interior path did not
  // write.
  solve(2.0, d.PPFD);
  ok(l.operating_point_kind() == Kind::Interior,
     "and the next wet solve is interior again, not the previous plant's kind");

  // Dim enough that gross assimilation cannot cover R_d: shutdown on LIGHT, on
  // wet soil. A rainfall sweep never reaches this, which is why it is its own kind.
  solve(1.0, 10.0);
  ok(l.assim_max_ < 0.0, "the dim solve takes the assim_max_ < 0 exit");
  ok(l.operating_point_kind() == Kind::ShadeDeath,
     "and is tagged shade-death, not a water shutdown");

  // A constrained optimum, then the same leaf asked to EVALUATE a prescribed
  // collar potential rather than optimise one. The prescribed point is not an
  // optimum of any kind, and would read `boundary-soil` if that path did not write.
  const std::vector<double> psi5{4.0, 4.25, 4.5, 4.75, 5.0};
  const std::vector<double> depth5{1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<double> root5(5, 1.0 / 5.0 / d.area_leaf);
  l.set_physiology(fixture::root_network(root5, depth5), d.PPFD, psi5, depth5,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::BoundarySoil,
     "psi_soil 4.0 over 5 layers is pinned to the wet bound");
  const double pinned_collar = l.opt_root_psi_;

  l.evaluate_root_collar_psi(pinned_collar);
  ok(l.operating_point_kind() == Kind::Prescribed,
     "evaluating a given collar potential is not an optimum, and says so");

  // ...and the optimising path takes it back, so `prescribed` is not sticky.
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::BoundarySoil,
     "and optimising again restores the pin");

  // set_traits returns the object to its just-constructed state, and the
  // classification is part of that state (hazard 10 / setup_clean_leaf).
  l.set_traits(l.vcmax_25, l.stem_c, l.stem_P50, l.roots_.root_c,
               l.roots_.root_P50, l.TF24_beta2, l.jmax_25, l.a,
               l.curv_fact_elec_trans, l.curv_fact_colim, l.TF24_cost_scale,
               l.R_d_25, l.JS22_gamma, l.CMax_a, l.CMax_b);
  ok(l.operating_point_kind() == Kind::Unsolved,
     "set_traits clears the classification with the rest of the solved state");

  // Neither failure kind is reachable from a driver sweep, and that is the point
  // of separating them from a pin: they mean "no plant is described".
  // test_collar_solve_refuses_rather_than_guessing drives both directly.
  ok(std::string(phylloptim::Leaf::operating_point_kind_name(
         Kind::SolverRefused)) == "solver-refused",
     "the failure kinds are named, not numbered");
}

// The two branches of maximise_profit_over_collar that change the collar it
// returns. No driver reaches either, so both are driven through the public
// maximise_profit_over_collar with brackets built to land on them.
void test_collar_solve_refuses_rather_than_guessing() {
  printf("the collar solve refuses a bracket it cannot resolve\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  const std::vector<double> psi{2.0, 2.25, 2.5, 2.75, 3.0};
  const std::vector<double> depth{1.0, 2.0, 3.0, 4.0, 5.0};

  // A comfortably interior optimum, and the interval it was found in.
  phylloptim::Leaf l = make_leaf(d, psi, depth);
  double bound_a = 0.0, bound_b = 0.0;
  ok(l.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(bound_a, bound_b),
     "the reference case has a feasible collar interval");
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::Interior,
     "and an interior optimum inside it to straddle");
  const double star = l.opt_root_psi_;

  // Bounds handed over inverted, so the interval runs dry to wet: dprofit <= 0 at
  // bound_a and >= 0 at bound_b, which makes each end a local maximum.
  const double dry_end = star + 0.9 * (bound_b - star);
  const double wet_end = star - 0.9 * (star - bound_a);

  phylloptim::Leaf p = make_leaf(d, psi, depth);
  double pa = 0.0, pb = 0.0;
  p.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(pa, pb);
  const double profit_dry = p.profit_psi_stem_TF(
      p.find_psi_stem_from_psi_root(dry_end, p.supply_psi_soil()), dry_end);
  const double profit_wet = p.profit_psi_stem_TF(
      p.find_psi_stem_from_psi_root(wet_end, p.supply_psi_soil()), wet_end);
  ok(profit_wet > profit_dry,
     "the wetter end is the better of the two, by construction");

  phylloptim::Leaf m = make_leaf(d, psi, depth);
  double ma = 0.0, mb = 0.0;
  m.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(ma, mb);
  const double refused = m.maximise_profit_over_collar<phylloptim::Leaf::CostCurve::TF24>(dry_end, wet_end);
  ok(m.operating_point_kind() == Kind::SolverRefused,
     "the solve reports that it could not resolve the bracket");
  ok(m.operating_point_kind() != Kind::BoundarySoil &&
         m.operating_point_kind() != Kind::BoundaryCrit &&
         m.operating_point_kind() != Kind::BoundaryRootCrit,
     "and does not pass it off as a constrained optimum");
  // The endpoint the solve returns is stepped a fraction of the width inside the
  // bound it came from, so compare against the bound rather than for equality.
  ok(std::abs(refused - wet_end) < 1e-5,
     "and returns the end with the higher profit");
  ok(std::abs(refused - dry_end) > 0.1, "which is not the drier end");

  // A bracket lying wholly at or below bound_a, where dprofit takes its
  // reversed-gradient exit: no usable gradient at either end. The
  // non-finite-gradient half of the same guard has no bracket that reaches it and
  // is not covered.
  //
  // ⚠️ THE BRACKET IS BELOW THE BOUND, NOT JUST INSIDE IT. G^-1 is the forward
  // table inverted, so zero flux returns the collar exactly and the infeasible
  // region is exactly (-inf, bound_a] -- one ulp above it answers. DO NOT fit the
  // inverse separately: at zero flux the round trip then puts psi_stem about 1e-07
  // below the collar, which is a reversed gradient the model does not have, and it
  // opens an infeasible sliver that wide above bound_a.
  phylloptim::Leaf s = make_leaf(d, psi, depth);
  double sa = 0.0, sb = 0.0;
  s.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(sa, sb);
  const double sliver = 1e-9;
  bool feasible = true;
  s.dprofit_at_collar_psi<phylloptim::Leaf::CostCurve::TF24>(sa, &feasible);
  ok(!feasible, "the wet bound admits no informative gradient");
  const double fallen_back = s.maximise_profit_over_collar<phylloptim::Leaf::CostCurve::TF24>(sa - sliver, sa);
  ok(s.operating_point_kind() == Kind::SolverRefused,
     "a bracket with no usable gradient at either end is refused too");
  ok(std::isfinite(fallen_back) && fallen_back >= sa - sliver &&
         fallen_back <= sa,
     "and the fallback stays inside the bracket it was handed");
}

// Hazard 3, re-measured rather than argued. The guide's constraint is that the
// argmax must vary SMOOTHLY with inputs. Golden section resolved it only to
// GSS_tol_abs, making it piecewise constant at fine scales -- so a sharper solve
// makes it smoother, which is the opposite of what the hazard reads like.
//
// ⚠️ This measures smoothness in a TRAIT. Hazard 3's actual concern is smoothness
// in plant state feeding the demographic growth-rate gradient, which cannot be
// measured from here while plant #591 blocks end-to-end validation. Same
// mechanism, so the same direction is expected -- but it is not the same test, and
// this one passing does not discharge the hazard.
void test_collar_argmax_is_smooth_in_a_trait() {
  printf("the collar argmax varies smoothly with a trait (hazard 3)\n");
  Drivers d;
  const int n = 11;
  const double base = 96.0, step = base * 1e-3;  // 0.1% steps in vcmax_25
  double collar[n];
  int distinct = 0;
  for (int i = 0; i < n; ++i) {
    phylloptim::Leaf l;
    l.vcmax_25 = base + i * step;
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    std::vector<double> ps{2.0}, depth{1.0}, root{1.0 / d.area_leaf};
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, ps, depth, d.K_s * d.theta / d.h, d.atm_vpd,
                     d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    collar[i] = l.opt_root_psi_;
    bool seen = false;
    for (int j = 0; j < i; ++j) {
      if (collar[j] == collar[i]) seen = true;
    }
    if (!seen) ++distinct;
  }
  // Golden section gave 6 distinct values out of 11 -- the staircase, tread width
  // ~GSS_tol_abs. Every step must now move the answer.
  ok(distinct == n, "every trait step moves the argmax (no staircase)");

  // Smoothness, as second differences against the step size. Golden section
  // measured 3.9e-4, i.e. the same size as the steps themselves: pure noise. The
  // root-find measured 3.4e-7, ~1000x smaller than its own steps.
  double worst_d2 = 0.0, mean_step = 0.0;
  for (int i = 1; i < n; ++i) {
    mean_step += std::abs(collar[i] - collar[i - 1]) / (n - 1);
  }
  for (int i = 2; i < n; ++i) {
    worst_d2 = std::max(worst_d2,
                        std::abs(collar[i] - 2 * collar[i - 1] + collar[i - 2]));
  }
  ok(worst_d2 < 0.02 * mean_step,
     "second differences are small against the step size, not comparable to it");
}

// #25's invariant, asserted rather than documented: every psi is a positive
// magnitude, so the soil->collar derivative is a CONDUCTANCE. Under the old signed
// convention this came back negative and had to be negated at the one call site
// that wanted a conductance -- the omission of that negation is the bug that
// motivated #8, and it is now unrepresentable.
void test_soil_conductance_is_positive() {
  printf("dE_up/d(collar suction) is a positive conductance\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers);
    for (int i = 0; i < layers; ++i) { ps[i] = 1.0 + 0.25 * i; depth[i] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double S = l.dE_from_soil_dpsi_collar(l.opt_root_psi_, l.roots_.psi_soil_);
    ok(std::isfinite(S) && S > 0.0,
       "conductance is finite and positive at " + std::to_string(layers) + " layers");
    // It really is dE/dT: a central difference on the uptake must agree.
    const double h = 1e-7;
    double up = 0.0, dn = 0.0;
    std::vector<double> buf(l.soil_consumption_.size(), 0.0);
    l.roots_.uptake_at(l.opt_root_psi_ + h, l.roots_.psi_soil_, buf, up);
    l.roots_.uptake_at(l.opt_root_psi_ - h, l.roots_.psi_soil_, buf, dn);
    near(S, (up - dn) / (2.0 * h), 1e-5,
         "conductance matches a central difference at " + std::to_string(layers) + " layers");
  }
}

// Issue #1: the two root vulnerability curves stop at the 1%-conductivity point,
// and a soil layer drier than that is an ordinary state -- plant's soil potential
// is capped at 1000 MPa, the grid at 6.82. Taken as odelia extrapolants the
// conductivity curve crossed zero at 7.3742 MPa (-20.35 at 1000) and the cumulative
// integral kept accumulating past a limit it had already reached, which made the
// per-layer mean resistance r_R_H = r_R_H_min * span / integral FALL as the layer
// dried: measured on the fixture below, a reverse flux 5.70x its bound.
//
// ⚠️ THE GOLDEN FILE CANNOT SEE ANY OF THIS. Instrumented over the whole grid,
// the driest argument the integral ever gets is 7.0 MPa -- past the last knot
// (16 of 20644 lookups are) but short of 7.3132 where the cap binds -- and the
// conductivity curve is never read past 4.0. Every golden cell is bit-identical,
// so this test is the only thing standing behind the fix.
void test_root_vulnerability_is_bounded_past_its_grid() {
  printf("root vulnerability curves are bounded past their last knot\n");
  phylloptim::MultiLayerRoots r;
  r.setup_vulnerability(100);
  const double last_knot = r.root_vuln_from_psi.max();
  const double G_inf =
      phylloptim::cumulative_vulnerability_integral_limit(r.root_b, r.root_c);

  // On the grid the accessors ARE the bare splines, bit for bit. That is what
  // makes this fix golden-identical rather than merely golden-tolerable.
  for (double psi : {0.0, 0.5, 2.0, 4.0, 6.0}) {
    const std::string at = " at psi=" + std::to_string(psi);
    ok(r.root_vuln_at(psi) == r.root_vuln_from_psi.eval(psi),
       "f_r is the unmodified spline on the grid" + at);
    ok(r.root_vuln_integral_at(psi) == r.root_vuln_integral_from_psi.eval(psi),
       "G is the unmodified spline on the grid" + at);
  }

  // Past it, both stay in range and neither runs away.
  double prev_G = r.root_vuln_integral_at(last_knot);
  for (double psi : {7.0, 7.5, 8.0, 10.0, 100.0, 1000.0}) {
    const std::string at = " at psi=" + std::to_string(psi);
    const double f_r = r.root_vuln_at(psi);
    const double G = r.root_vuln_integral_at(psi);
    ok(f_r > 0.0 && f_r <= r.root_vuln_at(last_knot),
       "f_r stays positive and does not exceed the last knot" + at);
    ok(G >= prev_G && G <= G_inf, "G is monotone and at most G(inf)" + at);
    prev_G = G;
  }
  near(r.root_vuln_integral_at(1000.0), G_inf, 1e-12,
       "G saturates at its closed form (b/c)*Gamma(1/c)");

  // dG/dpsi has to agree with that: zero where the value is pinned, f_r where
  // it is not. duptake_dpsi differentiates the integral through this.
  ok(r.root_vuln_integral_deriv_at(1000.0) == 0.0,
     "dG/dpsi is zero where G is at its limit");
  near(r.root_vuln_integral_deriv_at(4.0),
       std::exp(-std::pow(4.0 / r.root_b, r.root_c)), 1e-4,
       "dG/dpsi is f_r inside the grid");

  // The wet end throws too, and MultiLayerRoots validates nothing it is handed.
  // NaN has to come back out for the caller's !isfinite(f_ri) guard to read it.
  ok(r.root_vuln_at(-0.5) == r.root_vuln_at(0.0),
     "a negative suction reads as the wet end");
  ok(std::isnan(r.root_vuln_at(std::numeric_limits<double>::quiet_NaN())),
     "and a NaN suction survives both clamps");

  // The live consequence. One rooted layer, unit horizontal resistance and no
  // vertical resistance, collar held at 1 MPa: the layer is drier, so it GAINS
  // water (hydraulic redistribution), and the gain is bounded by the whole area
  // under the conductivity curve, integral/r_R_H_min.
  const double bound = G_inf - phylloptim::cumulative_vulnerability_integral_at(
                                   1.0, r.root_b, r.root_c);
  double E_at_100 = 0.0, E_at_1000 = 0.0;
  for (double psi_soil : {100.0, 1000.0}) {
    r.set_soil_state(std::vector<double>{psi_soil}, std::vector<double>{1.0});
    phylloptim::RootNetwork n;
    n.r_R_H_min = {1.0};
    n.r_R_V_sum = {0.0};
    n.c_r_V = {1.0};
    n.c_r_H = {1.0};
    n.r_R_V = {0.0};
    r.set_root_network(n);
    r.begin_solve();
    std::vector<double> consumption(1, 0.0);
    double E_up = 0.0;
    r.uptake(1.0, consumption, E_up);
    const double E_i = consumption[0];
    ok(E_i < 0.0 && std::abs(E_i) <= bound * 1.001,
       "the reverse flux is bounded by the area under the curve at psi_soil=" +
           std::to_string(psi_soil));
    (psi_soil > 500.0 ? E_at_1000 : E_at_100) = E_i;
  }
  // A tenfold drier layer must not pump ten times harder. The 1.10e-4 between
  // them is gravitational head: integral * gravity_head * z_mid * (1/99 - 1/999),
  // the layer midpoint being 0.5 m.
  near(E_at_1000, E_at_100, 1e-3,
       "a tenfold drier layer does not become a stronger pump");
}

// The convention is now checkable, so check that it is checked: a caller still
// holding the pre-#25 signed vector must fail loudly, not run.
void test_signed_potentials_are_rejected() {
  printf("signed potentials are rejected at the input boundary\n");
  Drivers d;
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  bool threw = false;
  try {
    l.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD, {-2.0}, {1.0},
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "set_physiology rejects a negative psi_soil");

  phylloptim::Leaf ok_leaf = make_leaf(d, {2.0}, {1.0});
  ok_leaf.find_root_collar_psi();
  threw = false;
  try {
    ok_leaf.E_from_Soil_to_Root_Collar(2.5, {-2.0});
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "E_from_Soil_to_Root_Collar rejects a signed soil vector");

  // The constructor's half of the invariant.
  threw = false;
  try {
    phylloptim::Leaf bad(100, 2.04, -3.0, 2.65, 1.29, 1, 167 * 100, 0.3,
                   0.7, 0.99, 1e-8, 100, 1e-6, 1000, 46.32995);
    static_cast<void>(bad);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "the constructor rejects a negative stem_b");
}

// #24 / plant #584: the dry end of the collar bracket is clamped to
// root_psi_crit, the potential at which root conductivity is down to 5%. The clamp
// was written as std::max against a *signed* root_psi_crit, so it could never bind
// and the solver optimised over a collar drier than the root system can supply.
//
// The window is empty at this package's defaults, where psi_crit == root_psi_crit,
// which is why the golden file does not move. It opens whenever the stem's psi_crit
// is drier than the root's -- as it is in plant, by 1.2 MPa. Three regimes, all
// pinned here, because the middle one is the only place a *transpiring* operating
// point moves and the third is a behaviour the fix had to add rather than restore.
void test_root_psi_crit_clamp_binds() {
  printf("the collar bracket is clamped to root_psi_crit (#24)\n");
  Drivers d;
  const auto solve = [&](double psi_soil) {
    phylloptim::Leaf l;
    l.psi_crit = 5.91988;   // drier than root_psi_crit = 5.870283
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    std::vector<double> ps{psi_soil}, depth{1.0}, root{1.0 / d.area_leaf};
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, ps, depth, d.K_s * d.theta / d.h, d.atm_vpd,
                     d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    return l;
  };

  // Regime 1 -- the clamp does not bind (root_crit is wetter than root_psi_crit),
  // so nothing changes. Pinned so a future tightening cannot silently spread.
  {
    phylloptim::Leaf l = solve(5.80);
    ok(l.opt_root_psi_ < l.roots_.root_psi_crit,
       "below the window the collar stays inside the root limit anyway");
    ok(l.transpiration_ > 0.0, "and the leaf still transpires");
  }

  // Regime 2 -- the interval is TIGHTENED but still has room. The optimum is
  // genuinely interior here (measured 5.86989 against a bound of 5.870283, i.e.
  // 3.9e-4 inside it -- within GSS_tol_abs), so the assertion is the invariant the
  // clamp exists to enforce, not the boundary value: the collar no longer runs past
  // the root limit, and the leaf goes on transpiring.
  {
    phylloptim::Leaf l = solve(5.86);
    ok(l.opt_root_psi_ <= l.roots_.root_psi_crit,
       "in the window the collar does not pass root_psi_crit");
    ok(l.opt_root_psi_ > l.roots_.root_psi_crit - 1e-3,
       "and it sits at the clamp, within the GSS tolerance");
    ok(l.transpiration_ > 0.0, "and the leaf still transpires there");
  }

  // Regime 3 -- the clamp lands BELOW root_zero_E, the collar at which uptake is
  // zero. Drawing any water would need a collar past the root limit, so there is no
  // feasible transpiring operating point and the answer is shut-down. Nothing
  // handled this before #24, because with the clamp dead it could not arise.
  {
    phylloptim::Leaf l = solve(5.90);
    near(l.opt_root_psi_, l.roots_.root_psi_crit, 1e-12,
         "past the window the collar sits at root_psi_crit");
    near(l.transpiration_, 0.0, 1e-300, "and the leaf is shut down, not optimising");
    near(l.opt_psi_stem_, l.psi_crit, 1e-12, "with the stem held at psi_crit");
    ok(l.opt_root_psi_ <= l.roots_.root_psi_crit,
       "the collar never passes root_psi_crit in any regime");
  }
}

void test_lambda_equals_dA_dE_single_layer() {
  printf("marginal cost of water: analytic lambda vs dA/dE (stem free)\n");
  Drivers d;
  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    // the TF24 stem route holds the collar fixed at psi_soil_[0] and optimises
    // the stem, so the single-layer lambda is the one that applies here.
    l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
    const double psi = l.opt_psi_stem_;
    const double eps = 1e-6;
    l.set_leaf_states_rates_from_psi_stem(psi + eps, psi_soil);
    const double A1 = l.assim_colimited_, E1 = l.transpiration_;
    l.set_leaf_states_rates_from_psi_stem(psi - eps, psi_soil);
    const double A0 = l.assim_colimited_, E0 = l.transpiration_;
    const double fd = (A1 - A0) / (E1 - E0);
    l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
    // Tolerance is set by the optimiser: GSS_tol_abs is 1e-3 on psi, so the
    // first-order condition only holds to about that accuracy.
    near(l.marginal_cost_water() / fd, 1.0, 1e-3,
         "lambda/(dA/dE) at psi_soil=" + std::to_string(psi_soil));
    ok(l.marginal_cost_water() > 0.0, "lambda is positive");
  }
}

// Cowan & Farquhar (1977): the objective is `A - lambda*E`, so its first-order
// condition in psi_stem is dA/dE == lambda EXACTLY -- not approximately, and not
// mediated by a vulnerability curve. That identity IS the model, so this is its
// defining assertion rather than a consistency check on it.
//
// Contrast test_lambda_equals_dA_dE_single_layer above, which asks the same
// question of the TF24 cost and has to go through `marginal_cost_water()` to get
// an answer: there lambda is a DERIVED property of the operating point, here it
// is the prescribed input. The two tests look alike and are asking opposite
// questions.
void test_cowan_farquhar_equates_dA_dE_to_lambda() {
  printf("Cowan-Farquhar: dA/dE == the prescribed lambda at the optimum\n");
  Drivers d;
  d.PPFD = 1500.0;
  // Reported so the tolerance below can be set from the measurement rather than
  // guessed, and so a solver change that loosens the condition shows up as a
  // number moving rather than only as a pass/fail.
  double worst = 0.0;
  int interior_rows = 0, bound_rows = 0;
  // ⚠️ lambda's SCALE is set by the model, not chosen freely: it is umol C per kg
  // of water, and at these drivers the leaf's own marginal cost of water runs
  // 9e4 to 3e5 (see `marginal_cost_water()`). A lambda far below that band makes
  // water free and every optimum pins wide open at psi_crit; far above it and
  // they all pin shut. Either way the identity below is never exercised, and the
  // test passes while asserting nothing -- so keep these values inside the band
  // and keep the interior-row count in the printed line.
  for (double lambda : {8.0e4, 1.5e5, 3.0e5}) {
    for (double psi_soil : {0.5, 1.0, 2.0}) {
      phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
      l.CF77_lambda_ = lambda;
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
      const double psi = l.opt_psi_stem_;

      // An interior optimum is where the condition applies. At a bound the
      // gradient is genuinely non-zero and dA/dE != lambda is the right answer,
      // so assert the bound instead of loosening the tolerance.
      const bool interior = psi > psi_soil + 1e-6 && psi < l.psi_crit - 1e-6;
      if (!interior) {
        ++bound_rows;
        ok(psi == psi_soil || psi == l.psi_crit,
           "a non-interior optimum sits exactly on a bound, lambda=" +
               std::to_string(lambda));
        continue;
      }
      ++interior_rows;

      const double eps = 1e-6;
      l.set_leaf_states_rates_from_psi_stem(psi + eps, psi_soil);
      const double A1 = l.assim_colimited_, E1 = l.transpiration_;
      l.set_leaf_states_rates_from_psi_stem(psi - eps, psi_soil);
      const double A0 = l.assim_colimited_, E0 = l.transpiration_;
      const double fd = (A1 - A0) / (E1 - E0);

      worst = std::max(worst, std::abs(fd / lambda - 1.0));
      // Measured 2.3e-06 on the generating platform. The tolerance is set well
      // above that but far below the 1e-3 an argmax-resolution argument would
      // give, because the condition is STATIONARY: a displaced argmax changes
      // dA/dE only to second order, so this ratio is much better conditioned
      // than the potential it is evaluated at.
      near(fd / lambda, 1.0, 1e-4,
           "dA/dE == lambda at lambda=" + std::to_string(lambda) +
               " psi_soil=" + std::to_string(psi_soil));
    }
  }
  printf("    %d interior rows, %d on a bound | worst |dA/dE / lambda - 1|: %.3e\n",
         interior_rows, bound_rows, worst);
  ok(interior_rows > 0,
     "at least one row has an interior optimum, so the identity was exercised");
}

// THE STRONGER FORM OF THE SAME CLAIM, and the one that makes lambda a common
// currency rather than a per-model parameter: the TF24 cost and the
// Cowan-Farquhar cost are different functions of psi, but both optima satisfy
// dA/dE == lambda. So handing Cowan-Farquhar the lambda that TF24 IMPLIES at its
// own optimum must land it on that same optimum.
//
// This is a much sharper instrument than a finite-difference ratio: it compares
// two argmaxes computed by two objectives, with no derivative approximation in
// between, and it fails if `marginal_cost_water()` reports the wrong quantity for
// either curve.
void test_cowan_farquhar_reproduces_the_TF_optimum() {
  printf("Cowan-Farquhar at TF24's implied lambda finds TF24's optimum\n");
  Drivers d;
  d.PPFD = 1500.0;
  double worst = 0.0;
  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    phylloptim::Leaf ref = make_leaf(d, {psi_soil}, {1.0});
    ref.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
    const double psi_tf = ref.opt_psi_stem_;
    const double lambda_implied = ref.marginal_cost_water();

    phylloptim::Leaf cf = make_leaf(d, {psi_soil}, {1.0});
    cf.CF77_lambda_ = lambda_implied;
    cf.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();

    worst = std::max(worst, std::abs(cf.opt_psi_stem_ - psi_tf));
    // Measured 1.2e-05 MPa on the generating platform. The tolerance is 5e-4
    // rather than that, because this is an ARGMAX difference and argmax-derived
    // quantities are the sqrt-amplified class: they disagree across platforms by
    // up to 1.4e-04. Tightening to the local measurement would make this fail off
    // macOS/arm64 for a reason that is not a regression.
    near(cf.opt_psi_stem_, psi_tf, 5e-4,
         "the two optima coincide at psi_soil=" + std::to_string(psi_soil));
  }
  printf("    worst |psi_CF - psi_TF| at TF24's implied lambda: %.3e MPa\n", worst);
}

// The cost is lambda*E and E is exactly zero when the two potentials coincide, so
// the no-flow profit is exactly -R_d whatever lambda is. That makes this the one
// degenerate branch in the file whose value can be predicted rather than recorded.
void test_cowan_farquhar_closed_state() {
  printf("Cowan-Farquhar: the closed-stomata profit is -R_d, for any lambda\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {6.0}, {1.0}); // drier than psi_crit
  l.CF77_lambda_ = 800.0;
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
  near(l.opt_psi_stem_, 6.0, 1e-12, "the stem is held at the soil potential");
  near(l.transpiration_, 0.0, 1e-300, "transpiration is exactly zero");
  near(l.hydraulic_cost_, 0.0, 1e-300, "so the cost is exactly zero");
  near(l.profit_, -l.R_d_, 1e-12, "and profit is -R_d");
  // Every reported field describes that point rather than a search probe -- the
  // property the evaluate-at-the-closed-point convention exists to give.
  near(l.profit_, l.assim_colimited_ - l.hydraulic_cost_, 1e-12,
       "profit == assim - cost holds in the degenerate branch too");
}

void test_cowan_farquhar_refuses_an_unset_lambda() {
  printf("Cowan-Farquhar refuses an unset lambda\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {1.0}, {1.0});
  ok(!std::isfinite(l.CF77_lambda_), "CF77_lambda_ starts unset");
  bool threw = false;
  std::string what;
  try {
    l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
  } catch (const std::exception &e) {
    threw = true;
    what = e.what();
  }
  ok(threw, "it refuses rather than maximising a NaN objective");
  ok(what.find("CF77_lambda_") != std::string::npos, "and the message names CF77_lambda_");
}

// dprofit_dpsi_stem against a central difference of the objective it claims to
// differentiate, for both cost curves.
//
// THE POINT OF THE TEST is that the two halves come from different places: the
// derivative is analytic (forward AD for the TF24 curve, a closed form for
// Cowan-Farquhar, and the implicit function theorem through the ci root-find for
// the assimilation term in both), while the reference differences
// `profit_psi_stem_*` itself. So this checks the derivative against the model
// rather than against another derivative -- which is the failure that matters:
// a derivative of a nearby function is smooth, plausible and wrong.
void test_dprofit_dpsi_stem_matches_a_finite_difference() {
  printf("dprofit/dpsi_stem: analytic vs a central difference, both curves\n");
  Drivers d;
  d.PPFD = 1500.0;
  // ⚠️ THE STEP IS MEASURED, NOT PICKED. The differenced quantity is computed
  // through the ci root-find, so it carries that solver's tolerance as a fixed
  // absolute noise floor, and a central difference divides it by h. The relative
  // error is therefore ~1/h until truncation takes over: measured at one of the
  // rows below it runs 1.3e-01, 1.1e-02, 1.1e-03, 1.1e-04, 1.1e-05, 1.9e-07 for
  // h from 1e-8 to 1e-3, then back up to 1.3e-04 at 1e-2. So 1e-3 is the floor of
  // that V, and a "small" step like 1e-6 sits three orders up the noise side of
  // it -- which reads as a wrong derivative and is not one.
  const double h = 1e-3;
  double worst_tf = 0.0, worst_cf = 0.0;
  int rows = 0;

  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    // Sample the interior of [psi_soil, psi_crit] rather than only the optimum:
    // away from the optimum dprofit is LARGE, so a wrong chain rule shows up as a
    // relative error instead of hiding under a near-zero value.
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    const double lo = psi_soil, hi = l.psi_crit;
    for (double frac : {0.15, 0.35, 0.55, 0.75}) {
      const double psi = lo + frac * (hi - lo);
      ++rows;

      bool feasible = false;
      const double an_tf =
          l.dprofit_dpsi_stem<phylloptim::Leaf::CostCurve::TF24>(psi, psi_soil,
                                                                 &feasible);
      const double up_tf = l.profit_psi_stem_TF(psi + h, psi_soil);
      const double dn_tf = l.profit_psi_stem_TF(psi - h, psi_soil);
      const double fd_tf = (up_tf - dn_tf) / (2.0 * h);
      if (feasible) {
        worst_tf = std::max(worst_tf, std::abs(an_tf / fd_tf - 1.0));
        near(an_tf / fd_tf, 1.0, 5e-5,
             "TF24 dprofit at frac=" + std::to_string(frac) +
                 " psi_soil=" + std::to_string(psi_soil));
      }

      l.CF77_lambda_ = 1.5e5;
      const double an_cf =
          l.dprofit_dpsi_stem<phylloptim::Leaf::CostCurve::CF77>(
              psi, psi_soil, &feasible);
      const double up_cf = l.profit_psi_stem_CF77(psi + h, psi_soil);
      const double dn_cf = l.profit_psi_stem_CF77(psi - h, psi_soil);
      const double fd_cf = (up_cf - dn_cf) / (2.0 * h);
      if (feasible) {
        worst_cf = std::max(worst_cf, std::abs(an_cf / fd_cf - 1.0));
        near(an_cf / fd_cf, 1.0, 5e-5,
             "Cowan-Farquhar dprofit at frac=" + std::to_string(frac) +
                 " psi_soil=" + std::to_string(psi_soil));
      }
    }
  }
  // Measured 8.9e-07 (TF24) and 4.1e-06 (Cowan-Farquhar). The bound is 5e-05
  // rather than either: the residual is the ci solver's noise floor divided by h,
  // and that floor is libm-dependent, so it moves across platforms. Nothing is
  // lost by the slack -- a missing or mis-signed chain-rule term is an O(1)
  // relative error, not an O(1e-5) one.
  printf("    %d points | worst relative error: TF24 %.3e, Cowan-Farquhar %.3e\n",
         rows, worst_tf, worst_cf);
}

// The derivative and the optimiser have to agree about where the optimum is: the
// solver returns an argmax, so dprofit there must be ~zero, and AWAY from it must
// not be. The second half is what makes this more than a tautology -- a
// derivative that returned zero everywhere would pass the first half alone.
void test_dprofit_dpsi_stem_vanishes_at_the_optimum() {
  printf("dprofit/dpsi_stem is ~0 at the optimum and not away from it\n");
  Drivers d;
  d.PPFD = 1500.0;
  for (double psi_soil : {0.5, 1.0, 2.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    l.CF77_lambda_ = 1.5e5;
    l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
    const double psi = l.opt_psi_stem_;
    if (psi <= psi_soil + 1e-6 || psi >= l.psi_crit - 1e-6) {
      continue; // pinned: the gradient is genuinely non-zero there
    }
    bool feasible = false;
    const double at =
        l.dprofit_dpsi_stem<phylloptim::Leaf::CostCurve::CF77>(
            psi, psi_soil, &feasible);
    ok(feasible, "the optimum admits an informative gradient");
    // Scaled by the derivative's own magnitude a little away from the optimum,
    // so the bound means "small compared with what this function returns here"
    // rather than "small in absolute carbon units".
    const double away =
        l.dprofit_dpsi_stem<phylloptim::Leaf::CostCurve::CF77>(
            psi_soil + 0.2 * (l.psi_crit - psi_soil), psi_soil, &feasible);
    ok(std::abs(at) < 1e-3 * std::abs(away),
       "dprofit at the optimum is negligible against its scale, psi_soil=" +
           std::to_string(psi_soil));
    ok(std::abs(away) > 0.0, "and the scale itself is non-zero");
  }
}

// evaluate_psi_stem: the prescribed-point counterpart of the optimisers, and the
// case where a gradient is a partial at fixed psi rather than a total through a
// moving argmax.
void test_evaluate_psi_stem_prescribes_rather_than_optimises() {
  printf("evaluate_psi_stem: a prescribed operating point, clamped and tagged\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  const double psi_soil = 1.0;

  phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
  l.CF77_lambda_ = 1.5e5;
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
  const double psi_opt = l.opt_psi_stem_, profit_opt = l.profit_;

  // A point deliberately away from the optimum is NOT stationary, and is worth
  // strictly less. Both halves matter: the first says the function honoured the
  // request, the second says the optimiser was actually finding a maximum.
  const double psi_off = psi_soil + 0.25 * (l.psi_crit - psi_soil);
  const double p_off = l.evaluate_psi_stem<CC::CF77>(psi_off);
  near(l.opt_psi_stem_, psi_off, 1e-12, "it sits exactly where it was told to");
  ok(l.operating_point_kind() == phylloptim::Leaf::OperatingPointKind::Prescribed,
     "and reports itself as prescribed rather than solved");
  ok(p_off < profit_opt, "a prescribed point off the optimum is worth less");

  bool feasible = false;
  const double slope =
      l.dprofit_dpsi_stem<CC::CF77>(psi_off, psi_soil, &feasible);
  ok(feasible && std::abs(slope) > 1e-6,
     "dprofit is genuinely non-zero there, unlike at the optimum");

  // Asking for the optimum back reproduces it, which is what makes the two entry
  // points comparable at all.
  const double p_at = l.evaluate_psi_stem<CC::CF77>(psi_opt);
  near(p_at, profit_opt, 1e-12, "prescribing the optimum reproduces its profit");

  // Outside the interval the target is clamped, not refused -- and the clamp is
  // where a prescribed gradient stops being a pure partial, because the bound
  // moves with the traits.
  l.evaluate_psi_stem<CC::CF77>(l.psi_crit + 5.0);
  near(l.opt_psi_stem_, l.psi_crit, 1e-12, "a target past psi_crit clamps to it");
  l.evaluate_psi_stem<CC::CF77>(0.0);
  near(l.opt_psi_stem_, psi_soil, 1e-12, "and one below psi_soil clamps up to it");

  // The TF24 curve takes the same route and needs no lambda.
  phylloptim::Leaf t = make_leaf(d, {psi_soil}, {1.0});
  const double p_tf = t.evaluate_psi_stem<CC::TF24>(psi_off);
  near(p_tf, t.profit_psi_stem_TF(psi_off, psi_soil), 1e-12,
       "the TF24 arm agrees with its own objective at that point");
}

// EVERY cost curve reports its emergent lambda on the same axis, which is what
// makes `lambda_emergent_` a shared output rather than three fields wearing one
// name. Two independent checks, one per curve:
//
//   * Cowan-Farquhar prices water at a constant, so its emergent lambda IS the
//     prescribed one -- exactly, since dC/dpsi = lambda*(dE/dpsi) identically;
//   * TF24's must equal a finite difference of its own cost over its own E.
//
// The ProfitMax case is checked the same way in
// test_profitmax_emergent_lambda_matches_a_finite_difference.
// ProfitMax's emergent lambda has a derivation rather than a closed form lying
// around, so it gets its own check: against a finite difference of its OWN cost
// over its own transpiration, scaled by |A|max to restore carbon units. Run with
// both gates on, which is the arm where the thermal term is not zero.
void test_profitmax_emergent_lambda_matches_a_finite_difference() {
  printf("ProfitMax's emergent lambda against a finite difference of dC/dE\n");
  Drivers d;
  d.PPFD = 1500.0;
  const double h = 1e-3;
  double worst = 0.0;
  int rows = 0;

  for (double t : {25.0, 40.0}) {
    for (bool eb : {false, true}) {
      d.leaf_temp = t;
      phylloptim::Leaf l = make_leaf(d, {1.0}, {1.0});
      l.use_energy_balance_ = eb;
      l.use_thermal_cost_ = true;
      if (eb) { l.Rn_ = 300.0; l.d_ = 0.05; l.wind_speed_ = 2.0; }
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
      const double psi = l.opt_psi_stem_;
      if (psi <= 1.0 + h || psi >= l.psi_crit - h) continue;

      // dC/dpsi of the NORMALISED cost, and dE/dpsi, both by difference.
      auto cost_and_E = [&](double p) {
        l.profit_psi_stem_ProfitMax(p, 1.0);
        return std::pair<double, double>(l.hydraulic_cost_norm_ + l.thermal_cost_,
                                        l.transpiration_);
      };
      const auto up = cost_and_E(psi + h);
      const auto dn = cost_and_E(psi - h);
      const double fd = l.profitmax_A_max() *
                        ((up.first - dn.first) / (up.second - dn.second));

      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
      ++rows;
      worst = std::max(worst, std::abs(l.lambda_emergent() / fd - 1.0));
      near(l.lambda_emergent() / fd, 1.0, 5e-3,
           "ProfitMax emergent lambda at T=" + std::to_string(t) +
               " eb=" + std::to_string(int(eb)));
    }
  }
  printf("    %d rows | worst relative error: %.3e\n", rows, worst);
  ok(rows > 0, "at least one row had an interior optimum to check");
}

void test_every_curve_reports_an_emergent_lambda() {
  printf("every cost curve reports its emergent lambda on one axis\n");
  Drivers d;
  d.PPFD = 1500.0;
  const double h = 1e-3;   // measured step; see the dprofit test for why

  for (double psi_soil : {0.5, 1.0, 2.0}) {
    // --- Cowan-Farquhar: exact ---------------------------------------------
    {
      phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
      l.CF77_lambda_ = 1.5e5;
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
      ok(l.lambda_emergent() == l.CF77_lambda_,
         "Cowan-Farquhar's emergent lambda is its prescribed one, exactly");
    }
    // --- TF24: against a finite difference of dC/dE -------------------------
    {
      phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
      const double psi = l.opt_psi_stem_;
      if (psi <= psi_soil + h || psi >= l.psi_crit - h) continue;
      const double dC =
          (l.hydraulic_cost_TF(psi + h) - l.hydraulic_cost_TF(psi - h)) / (2 * h);
      l.profit_psi_stem_TF(psi + h, psi_soil);
      const double E1 = l.transpiration_;
      l.profit_psi_stem_TF(psi - h, psi_soil);
      const double E0 = l.transpiration_;
      const double fd = dC / ((E1 - E0) / (2 * h));
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
      near(l.lambda_emergent() / fd, 1.0, 1e-4,
           "TF24's emergent lambda is dC/dE at psi_soil=" +
               std::to_string(psi_soil));
    }
  }
}

// The multi-layer identity from the companion manuscript:
//   lambda_multi = lambda_single * [1 + kmax*f(psi_r)/S],   S = dE_up/dpsi_r
// find_root_collar_psi optimises the COLLAR, so lambda_multi is what it
// equalises. The single-layer lambda should be badly wrong here -- that is the
// point of the correction, and this test pins the size of the error.
void test_multilayer_lambda_identity() {
  printf("multi-layer lambda: series-resistance correction vs dA/dE (collar free)\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers), root(layers);
    for (int i = 0; i < layers; ++i) {
      ps[i] = 1.0 + 0.25 * i;
      depth[i] = 1.0 * (i + 1);
      root[i] = 1.0 / layers / d.area_leaf;
    }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double single = l.marginal_cost_water();
    const double multi = l.marginal_cost_water_multilayer();

    const double target = l.opt_root_psi_;
    const double eps = 1e-6;
    l.evaluate_root_collar_psi(target + eps);
    const double A1 = l.assim_colimited_, E1 = l.transpiration_;
    l.evaluate_root_collar_psi(target - eps);
    const double A0 = l.assim_colimited_, E0 = l.transpiration_;
    const double fd = (A1 - A0) / (E1 - E0);

    const std::string tag = std::to_string(layers) + " layers";
    ok(std::isfinite(multi), "lambda_multi is finite, " + tag);
    near(multi / fd, 1.0, 1e-3, "lambda_multi/(dA/dE), " + tag);
    // The correction bracket is >= 1, so the single-layer value must understate.
    ok(multi > single, "lambda_multi exceeds lambda_single, " + tag);
    ok(multi / single > 2.0 && multi / single < 12.0,
       "the correction factor is in the reported 2-12 range, " + tag);
    // And the single-layer value should be badly wrong when the collar is free,
    // which is exactly why the correction matters.
    ok(std::abs(single - fd) / fd > 0.5,
       "lambda_single is >50% wrong when the collar is free, " + tag);
  }
}

void test_g1_eff() {
  printf("equivalent Medlyn slope\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double g1 = l.g1_eff();
  ok(std::isfinite(g1) && g1 > 0.0, "g1_eff is finite and positive");
  // g1_eff is defined by inverting chi = g1/(g1 + sqrt(D)), so that must hold.
  const double chi = l.ci_ / l.ca_;
  near(g1 / (g1 + std::sqrt(l.atm_vpd_)), chi, 1e-12,
       "g1_eff inverts the USO relation exactly");
  // Drier soil closes stomata, lowering chi and therefore g1_eff.
  phylloptim::Leaf dry = make_leaf(d, {4.0}, {1.0});
  dry.find_root_collar_psi();
  ok(dry.g1_eff() < g1, "g1_eff falls as the soil dries");
  // And a higher marginal cost of water goes with a lower g1_eff.
  ok(dry.marginal_cost_water() > l.marginal_cost_water(),
     "lambda rises as the soil dries");
}

// Build a leaf with the energy-balance gate (and the wind model) configured BEFORE
// set_physiology runs, which is what the PM path needs: `ra_`, `Rn_` and `Tair_` are
// all derived there, and the non-finite-wind check is made there. make_leaf() calls
// set_physiology itself, so setting the gate on its return value is too late for
// anything set_physiology decides.
phylloptim::Leaf make_pm_leaf(const Drivers &d, std::vector<double> psi_soil,
                        std::vector<double> soil_depth, bool gate,
                        double wind_speed = 2.0, double leaf_dim = 0.05) {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = gate;
  l.wind_speed_ = wind_speed;
  l.d_ = leaf_dim;
  std::vector<double> root(psi_soil.size(),
                           1.0 / double(psi_soil.size()) / d.area_leaf);
  l.set_physiology(fixture::root_network(root, soil_depth), d.PPFD, psi_soil, soil_depth,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// Tleaf as a reported OUTPUT, on both paths and out of every exit.
//
// The gap this closes: on the PM path the leaf's temperature is solved per
// operating point, was used to re-derive the whole Farquhar block, and was then
// discarded -- so the one quantity that path exists to produce was the one thing a
// caller could not read. `leaf_temp_` is not it there, because set_physiology has
// reinterpreted that driver as AIR temperature.
//
// ⚠️ The shut-down cases are the point of this test, not an afterthought. They do
// not go through set_leaf_states_rates_from_psi_stem, so they are where hazard 8
// bites: an output a branch declines to write becomes the PREVIOUS solve's value,
// and plant drives every individual in a patch through one persistent Leaf.
void test_leaf_temperature_is_reported() {
  printf("Tleaf is reported, on both paths and from every exit\n");
  Drivers d;

  // Unsolved: the NA sentinel, like every other output.
  {
    phylloptim::Leaf bare;
    ok(!std::isfinite(bare.Tleaf_), "Tleaf is NA before a solve");
  }

  // Off the PM path Tleaf IS the driver, and exactly so -- a tolerance here would
  // pass on a value that had been round-tripped through the energy balance.
  for (double T : {15.0, 25.0, 40.0}) {
    Drivers dt = d;
    dt.leaf_temp = T;
    phylloptim::Leaf l = make_leaf(dt, {2.0}, {1.0});
    l.find_root_collar_psi();
    ok(l.Tleaf_ == T,
       "off the energy-balance path Tleaf is the driver, at T = " +
           std::to_string(T));
    // ⚠️ NOT `.back()`: `shadow_cost` and then `lambda_emergent` were appended
    // after `Tleaf`, so the last element is neither. Addressed from the end by
    // name-order instead, which is what the R side's `.operating_point_names`
    // pins -- and that test is what catches this drifting again.
    const std::vector<double> ov = l.operating_point_values();
    ok(ov[ov.size() - 3] == l.Tleaf_,
       "operating_point_values() reports it, at T = " + std::to_string(T));
  }

  // On the PM path it is not the driver, and it is hotter: the leaf absorbs
  // radiation and sheds only part of it as latent heat.
  {
    phylloptim::Leaf eb = make_pm_leaf(d, {2.0}, {1.0}, true);
    eb.find_root_collar_psi();
    ok(eb.Tleaf_ != eb.leaf_temp_,
       "on the energy-balance path Tleaf is not the leaf_temp driver");
    ok(eb.Tleaf_ == eb.leaf_temp_from_E(eb.transpiration_),
       "and it is the temperature the solve's own transpiration implies");
    ok(eb.Tleaf_ > eb.Tair_, "the leaf runs hotter than the air here");
  }

  // ONE leaf, driven twice, second time into shut-down: the case a fresh-object
  // test cannot see. Both paths, because the two exits differ.
  for (bool gate : {false, true}) {
    const std::string what = gate ? " (energy balance on)" : " (prescribed)";
    phylloptim::Leaf l = make_pm_leaf(d, {2.0}, {1.0}, gate);
    l.find_root_collar_psi();
    const double wet = l.Tleaf_;
    ok(std::isfinite(wet), "a wet solve reports a finite Tleaf" + what);

    // Drier than psi_crit, so the collar solve takes a shut-down exit.
    std::vector<double> root{1.0 / d.area_leaf}, dry{6.5}, depth{1.0};
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, dry, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    ok(l.transpiration_ == 0.0, "the second solve really did shut down" + what);
    // The claim is not "it changed" -- it is that the value describes THIS
    // operating point. At zero transpiration there is no latent cooling, so the
    // PM answer is leaf_temp_from_E(0); off that path it is still the driver.
    const double want = gate ? l.leaf_temp_from_E(0.0) : l.leaf_temp_;
    ok(l.Tleaf_ == want,
       "a shut-down leaf reports its own temperature, not the last solve's" +
           what);
    if (gate) {
      ok(l.Tleaf_ > wet,
         "and a leaf that has stopped transpiring is hotter than one that had"
         " not" + what);
    }
  }
}

// One temperature per reported operating point, at the exits that transpire
// nothing (#105).
//
// The gap: on the PM path both shut-down exits reported a leaf temperature and a
// respiration rate belonging to DIFFERENT temperatures. They do not go through
// set_leaf_states_rates_from_psi_stem, so they inherited the Tair baseline
// set_physiology derives and then described a state whose transpiration is zero.
//
// ⚠️ `operating_points.tsv` cannot see any of this -- the golden grid runs with
// the gate off, where Tleaf IS leaf_temp_ and there is nothing to disagree with.
// A bit-identical golden run is not evidence here, which is why this test exists
// rather than a regenerated baseline.
//
// ⚠️ The two exits move in OPPOSITE directions, so a test that only checked "R_d
// went up" would pass on half the fix. Rn is proportional to PPFD with a fixed
// longwave offset subtracted, so a hydraulically shut leaf in full sun is hotter
// than the air and a shade-dead one is cooler.
void test_shutdown_reports_one_temperature() {
  printf("a shut-down leaf reports one temperature, not two\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  d.leaf_temp = 30.0;  // AIR temperature on this path

  // (psi_soil, PPFD, expected kind) -- one row per exit. Reached through the
  // drivers rather than by calling the exits directly, so the test breaks if a
  // future bracket change stops routing here.
  struct Case {
    const char *name;
    double psi_soil;
    double PPFD;
    Kind kind;
    bool hotter_than_air;
  };
  const Case cases[] = {
      // Soil drier than psi_crit: no collar potential both moves water and stays
      // inside the critical potentials. Full sun, so the leaf runs hot.
      {"hydraulic shutdown", 6.0, 900.0,
       Kind::HydraulicShutdown, true},
      // Wet soil, no light: gross assimilation at ci = ca cannot cover R_d. Rn is
      // NEGATIVE here, so this leaf is cooler than the air.
      {"shade death", 1.0, 1.0, Kind::ShadeDeath, false},
  };

  for (const Case &c : cases) {
    Drivers dc = d;
    dc.PPFD = c.PPFD;
    phylloptim::Leaf l = make_pm_leaf(dc, {c.psi_soil}, {1.0}, true);
    l.find_root_collar_psi();
    const std::string what = std::string(" (") + c.name + ")";

    ok(l.operating_point_kind() == c.kind,
       "the drivers really reach this exit" + what);
    ok(l.transpiration_ == 0.0, "and it transpires nothing" + what);

    // The premise: the two temperatures differ, so there is something to get
    // wrong. Without this the assertions below would pass on a leaf at Tair.
    ok(l.Tleaf_ != l.Tair_, "Tleaf is not air temperature" + what);
    ok(c.hotter_than_air ? (l.Tleaf_ > l.Tair_) : (l.Tleaf_ < l.Tair_),
       "and it sits on the expected side of it" + what);

    // The fix, stated as the identity it restores: the respiration in force is
    // the one the model's OWN curve gives at the temperature being reported.
    // Bit-exact -- same function, same argument, so a tolerance here would pass
    // on a value derived at some third temperature.
    const double rd_reported = l.R_d_;
    const double gamma_reported = l.gamma_;
    const double Tleaf = l.Tleaf_, Tair = l.Tair_;

    l.update_temperature_dependent_params(Tleaf);
    ok(l.R_d_ == rd_reported, "R_d is the curve's value at the reported Tleaf" + what);
    ok(l.gamma_ == gamma_reported,
       "and so is the compensation point" + what);

    l.update_temperature_dependent_params(Tair);
    ok(l.R_d_ != rd_reported,
       "and it is NOT the Tair value -- the two really differ here" + what);
  }

  // The accounting identity every branch is supposed to satisfy, re-checked at
  // these exits because they form profit by hand rather than through
  // profit_psi_stem_TF. `ci_` is in here because the shade-death exit did not
  // write it at all: it reported the previous solve's internal CO2, or the NA
  // sentinel on a cold object.
  for (const Case &c : cases) {
    Drivers dc = d;
    dc.PPFD = c.PPFD;
    phylloptim::Leaf l = make_pm_leaf(dc, {c.psi_soil}, {1.0}, true);
    l.find_root_collar_psi();
    const std::string what = std::string(" (") + c.name + ")";

    ok(l.assim_colimited_ == -l.R_d_,
       "net assimilation is -R_d at zero transpiration" + what);
    near(l.profit_, l.assim_colimited_ - l.hydraulic_cost_TF(l.opt_psi_stem_),
         1e-12, "profit == assimilation - hydraulic cost" + what);
    ok(std::isfinite(l.ci_), "ci is written, and finite" + what);
    ok(l.ci_ == l.gamma_ * l.umol_per_mol_to_Pa_,
       "ci sits at the compensation point of the reported temperature" + what);
    // #93's leaf-to-air deficit is on the same footing: reported, read by
    // g1_eff(), and derived from a temperature -- so it has to be THIS
    // temperature and not set_physiology's Tair baseline.
    ok(l.vpd_leaf_ == l.atm_vpd_ + (l.saturation_vapour_pressure(l.Tleaf_) -
                                    l.saturation_vapour_pressure(l.Tair_)),
       "the leaf-to-air deficit is at the reported Tleaf" + what);
    ok(l.vpd_leaf_ != l.atm_vpd_,
       "and it is not the air deficit -- there is something to get wrong" + what);
  }

  // Off the PM path nothing moves: leaf_temp_ IS the temperature the block was
  // derived at, so there was never a disagreement to fix. Asserted because the
  // change adds a conditional recompute, and the gate-off arm has to stay exactly
  // as it was -- this is the arm the golden file covers.
  {
    Drivers dc = d;
    dc.PPFD = 900.0;
    phylloptim::Leaf l = make_pm_leaf(dc, {6.0}, {1.0}, false);
    l.find_root_collar_psi();
    ok(l.operating_point_kind() == Kind::HydraulicShutdown,
       "the prescribed path shuts down here too");
    ok(l.Tleaf_ == l.leaf_temp_, "and Tleaf is still the driver");
    const double rd = l.R_d_;
    l.update_temperature_dependent_params(l.leaf_temp_);
    ok(l.R_d_ == rd, "R_d is the driver temperature's value, bit-for-bit");
  }
}

// `set_leaf_states_rates_from_psi_stem`'s zero-transpiration branch is
// SELF-CONTAINED: every quantity it reports is derived from its own `Tleaf_` and
// from nothing the caller left behind.
//
// ⚠️ WHY THE ASSERTIONS ARE ORDER-INDEPENDENCE RATHER THAN VALUES. The transpiring
// branch re-derives the temperature block per candidate by design, so this branch
// runs with the PREVIOUS candidate's block in place. A wrong value here is
// therefore plausible at any single history and only shows up as disagreement
// BETWEEN histories -- and because this function is the objective the collar solve
// maximises, that disagreement is an objective whose value depends on evaluation
// order, which is what hazard 3 forbids. Add a value assertion here if it helps,
// but the one that has teeth is "the same call twice, from two histories, agrees".
//
// ⚠️ `operating_points.tsv` cannot see any of this: the grid runs gate-off, where
// `set_leaf_vpd` returns `atm_vpd_` exactly and the block is never re-derived. A
// bit-identical golden run is not evidence about this test's subject.
void test_zero_E_branch_derives_its_own_block() {
  printf("the zero-E branch inside the objective reads its own temperature\n");
  Drivers d;
  d.leaf_temp = 30.0;  // AIR temperature on this path
  d.PPFD = 900.0;      // full sun, so a leaf that stops transpiring runs hot

  // The same (psi_stem, psi_upstream) from two histories: one object whose block
  // still holds set_physiology's baseline, one where a transpiring candidate has
  // moved it. Bit-exact, since the branch derives all of these from its own
  // `Tleaf_`.
  {
    phylloptim::Leaf cold = make_pm_leaf(d, {1.0}, {1.0}, true);
    cold.set_leaf_states_rates_from_psi_stem(1.0, 1.0);

    phylloptim::Leaf warm = make_pm_leaf(d, {1.0}, {1.0}, true);
    warm.set_leaf_states_rates_from_psi_stem(5.0, 1.0);  // transpiring
    const double T_transpiring = warm.Tleaf_;
    warm.set_leaf_states_rates_from_psi_stem(1.0, 1.0);  // the SAME zero-E point

    // The premise, asserted rather than assumed: the two histories must really
    // differ, or the comparisons below pass on a block that never moved. Re-check
    // these two if the drivers here are ever retuned.
    ok(T_transpiring != cold.Tleaf_,
       "the transpiring candidate sits at a different temperature");
    ok(cold.Tleaf_ > cold.Tair_,
       "and the zero-E leaf is hotter than the air -- no latent cooling");

    ok(warm.Tleaf_ == cold.Tleaf_, "Tleaf does not depend on the history");
    ok(warm.assim_colimited_ == cold.assim_colimited_,
       "and neither does assimilation");
    ok(warm.R_d_ == cold.R_d_, "nor respiration");
    ok(warm.vpd_leaf_ == cold.vpd_leaf_, "nor the leaf-to-air deficit");
    ok(warm.ci_ == cold.ci_, "nor the internal CO2");
  }

  // Order-independence alone would also hold if every history were wrong the same
  // way, so pin the value too: it is this branch's own temperature's. Same
  // identities the exits outside the solve carry, checked again here because this
  // branch reaches them by a different route.
  {
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    l.set_leaf_states_rates_from_psi_stem(1.0, 1.0);
    const double Tleaf = l.Tleaf_, Tair = l.Tair_;

    ok(l.transpiration_ == 0.0, "the branch really transpires nothing");
    ok(l.Tleaf_ == l.leaf_temp_from_E(0.0),
       "and its temperature is the E = 0 one");
    ok(l.ci_ == l.gamma_ * l.umol_per_mol_to_Pa_,
       "ci sits at the compensation point");
    // ⚠️ A TOLERANCE, DELIBERATELY, AND DO NOT TIGHTEN IT TO EXACT. Gross
    // assimilation at ci = gamma* is analytically zero, so net is -R_d, but the
    // co-limiting expression rounds -- 4.44e-16, one ULP. The exits outside the
    // solve assign `-R_d_` directly and so are bit-exact; this branch goes through
    // `assim_colimited`, and cannot be. The bit-exact statement that pins the
    // TEMPERATURE is the next one.
    near(l.assim_colimited_, -l.R_d_, 1e-14,
         "net assimilation is -R_d at zero transpiration");
    ok(l.assim_colimited_ == l.assim_colimited(l.ci_),
       "and it is this block's value at this branch's ci, bit-for-bit");
    ok(l.vpd_leaf_ == l.atm_vpd_ + (l.saturation_vapour_pressure(Tleaf) -
                                    l.saturation_vapour_pressure(Tair)),
       "and the deficit is the leaf-to-air one at that temperature");
    ok(l.vpd_leaf_ != l.atm_vpd_,
       "which is not the air deficit -- there is something to get wrong");

    // Bit-exact against the model's own curves at the reported temperature. A
    // tolerance here would pass on a block derived at some third temperature.
    const double rd = l.R_d_, gamma = l.gamma_, ci = l.ci_;
    l.update_temperature_dependent_params(Tleaf);
    ok(l.R_d_ == rd, "R_d is the curve's value at the reported Tleaf");
    ok(l.gamma_ == gamma, "and so is the compensation point");
    ok(ci == gamma * l.umol_per_mol_to_Pa_, "so ci follows it");
    l.update_temperature_dependent_params(Tair);
    ok(l.R_d_ != rd, "and it is NOT the Tair value -- the two really differ");
  }

  // Closure is still a profit CANDIDATE the optimiser can select (hazard 11) even
  // though it is excluded from |A|max, so its carbon and its temperature have to
  // stay mutually consistent. This used to be asserted through
  // `profitmax_scan_A_[0]`; |A|max is found by optimisation now and those arrays
  // are gone, so the property is asserted where it actually lives.
  {
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    l.use_thermal_cost_ = true;
    l.set_leaf_states_rates_from_psi_stem(1.0, 1.0);      // psi_stem == psi_soil
    const double a0 = l.assim_colimited_;
    ok(l.transpiration_ == 0.0, "at psi_stem == psi_soil no water moves");
    l.set_leaf_states_rates_from_psi_stem(1.0, 1.0);
    ok(l.assim_colimited_ == a0, "and the closure candidate is reproducible");
    l.update_temperature_dependent_params(l.leaf_temp_from_E(0.0));
    near(a0, -l.R_d_, 1e-14,
         "its carbon is -R_d at the temperature E = 0 implies");

    // And |A|max EXCLUDES it: a positive normaliser cannot be a respiration rate.
    l.prepare_profitmax();
    ok(l.profitmax_A_max_ > 0.0, "|A|max is a positive carbon gain");
    ok(l.profitmax_A_max_ != std::abs(a0),
       "and it is not the closure point's magnitude");
  }

  // Gate off, which is the arm the golden file covers and which must stay exactly
  // as it is: the recompute is inside `if (use_energy_balance_)`, and
  // `set_leaf_vpd` returns `atm_vpd_` there. This is the guard on that.
  {
    phylloptim::Leaf a = make_pm_leaf(d, {1.0}, {1.0}, false);
    a.set_leaf_states_rates_from_psi_stem(1.0, 1.0);
    phylloptim::Leaf b = make_pm_leaf(d, {1.0}, {1.0}, false);
    b.set_leaf_states_rates_from_psi_stem(5.0, 1.0);
    b.set_leaf_states_rates_from_psi_stem(1.0, 1.0);
    ok(a.Tleaf_ == a.leaf_temp_, "gate off: Tleaf is still the driver");
    ok(a.vpd_leaf_ == a.atm_vpd_, "gate off: the deficit is still the air one");
    ok(b.assim_colimited_ == a.assim_colimited_ && b.ci_ == a.ci_ &&
           b.vpd_leaf_ == a.vpd_leaf_,
       "gate off: the branch was already history-independent");
  }
}

void test_energy_balance_path_runs() {
  printf("Penman-Monteith energy-balance path\n");
  Drivers d;
  phylloptim::Leaf l = make_pm_leaf(d, {2.0}, {1.0}, false);
  l.find_root_collar_psi();
  const double A_prescribed = l.assim_colimited_;

  phylloptim::Leaf eb = make_pm_leaf(d, {2.0}, {1.0}, true);
  eb.find_root_collar_psi();
  ok(std::isfinite(eb.profit_), "energy-balance profit is finite");
  ok(std::isfinite(eb.assim_colimited_), "energy-balance assimilation is finite");
  ok(eb.assim_colimited_ != A_prescribed,
     "energy balance changes the operating point");
  const double Tleaf = eb.leaf_temp_from_E(eb.transpiration_);
  ok(Tleaf >= phylloptim::leaf_temp_min && Tleaf <= phylloptim::leaf_temp_max,
     "leaf temperature stays inside the physical clamp");
  // This line used to be the ONLY way to get at the leaf's temperature: derive it
  // again, outside the model, from an output. It is now reported, and the two
  // agreeing bit-for-bit is what says the reported value is the one the solve
  // actually ran at rather than a plausible recomputation.
  ok(eb.Tleaf_ == Tleaf, "the reported Tleaf is the solve's own");

  // The wind model really is what set ra_, rather than the fixed fallback.
  // ⚠️ Assert that rather than assuming it: setting wind_speed_/d_ and then
  // reassigning the leaf discards them, and the assertion still passes because
  // 2.0 / 0.05 are also the defaults.
  ok(std::isfinite(eb.ra_) && eb.ra_ > 0.0, "ra is finite and positive");
  near(eb.ra_, phylloptim::aerodynamic_resistance_coef * std::sqrt(0.05 / 2.0), 1e-12,
       "ra comes from the wind model, not the fixed fallback");
  ok(std::isfinite(eb.Rn_), "net radiation is finite");
}

// Isaac Towers' review of plant #567: the ra fallback accepted a non-finite wind
// speed as well as a zero one. Zero wind is physically ra -> infinity and a
// legitimate fallback; NA is a broken driver or an unset trait and should fail
// rather than silently produce a plausible number. Fixed in plant `76df7169` and
// carried into this package -- but the test lived only in plant, in a demo smoke
// test, while the code now lives here. Ported so the package that owns the
// contract also guards it.
void test_pm_wind_speed_validation() {
  printf("PM path fails fast on a non-finite wind speed (review: itowers1)\n");
  Drivers d;
  const double nan_v = std::numeric_limits<double>::quiet_NaN();

  // 1. gate ON + non-finite wind: fail fast.
  bool threw = false;
  try {
    make_pm_leaf(d, {2.0}, {1.0}, true, nan_v);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a non-finite wind speed throws on the energy-balance path");

  // ... and the same for the leaf dimension, which enters the same formula.
  threw = false;
  try {
    make_pm_leaf(d, {2.0}, {1.0}, true, 2.0, nan_v);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a non-finite leaf dimension throws on the energy-balance path");

  // 2. gate OFF + non-finite wind: fine, the wind model is never read.
  threw = false;
  try {
    phylloptim::Leaf off = make_pm_leaf(d, {2.0}, {1.0}, false, nan_v);
    off.find_root_collar_psi();
    ok(std::isfinite(off.profit_), "and still solves");
  } catch (const std::exception &) {
    threw = true;
  }
  ok(!threw, "a non-finite wind speed is ignored with the gate off");

  // 3. gate ON + ZERO wind: legitimate (ra -> infinity), falls back to the fixed
  // ra rather than erroring or producing an infinity.
  threw = false;
  try {
    phylloptim::Leaf zero = make_pm_leaf(d, {2.0}, {1.0}, true, 0.0);
    near(zero.ra_, phylloptim::aerodynamic_resistance_fixed, 1e-12,
         "zero wind falls back to the fixed ra");
    zero.find_root_collar_psi();
    ok(std::isfinite(zero.profit_), "and still solves");
  } catch (const std::exception &) {
    threw = true;
  }
  ok(!threw, "zero wind is a legitimate case, not an error");
}

// The behavioural content of plant's PM demo smoke test, which asserted these
// through the R shim over a Fick-vs-PM grid. The implementation is here now, so
// the sign of the effect is asserted here too.
void test_pm_leaf_temperature_response() {
  printf("PM leaf temperature: Tleaf == Tair off, departs from it on\n");
  for (double ppfd : {400.0, 2000.0}) {
    for (double tair : {20.0, 40.0}) {
      for (double vpd : {1.0, 3.0}) {
        Drivers d;
        d.PPFD = ppfd; d.leaf_temp = tair; d.atm_vpd = vpd;
        const std::string at = " at PPFD=" + std::to_string(int(ppfd)) +
                               " Tair=" + std::to_string(int(tair)) +
                               " VPD=" + std::to_string(int(vpd));

        phylloptim::Leaf fick = make_pm_leaf(d, {2.0}, {1.0}, false);
        fick.find_root_collar_psi();
        ok(std::isfinite(fick.profit_) && std::isfinite(fick.assim_colimited_),
           "Fick outputs are finite" + at);

        phylloptim::Leaf pm = make_pm_leaf(d, {2.0}, {1.0}, true);
        pm.find_root_collar_psi();
        ok(std::isfinite(pm.profit_) && std::isfinite(pm.assim_colimited_),
           "PM outputs are finite" + at);

        // Hot and bright is where PM matters: the leaf runs warmer than the air.
        if (ppfd == 2000.0 && tair == 40.0) {
          const double Tleaf = pm.leaf_temp_from_E(pm.transpiration_);
          ok(Tleaf > pm.Tair_, "the leaf is warmer than the air when hot and bright");
        }
      }
    }
  }
}

// ============================================================================
// MEASUREMENT, NOT A GUARD. What the collar solve currently does under the
// energy-balance gate, recorded before anything is changed.
//
// The claim being measured: with `use_energy_balance_` on, the solver does not
// find the maximum of the objective it evaluates. Two separate defects compound.
//
//  1. STALENESS. `dprofit_at_collar_psi` never calls
//     `set_leaf_states_rates_from_psi_stem`, and `set_physiology` gates its
//     temperature cache on `!use_energy_balance_`, so the derivative is
//     evaluated with vcmax_/jmax_/gamma_/km_/R_d_ left at the AIR-temperature
//     baseline while the objective is evaluated at Tleaf(E(psi)).
//  2. THE MISSING CHAIN TERM. Even seated correctly, the derivative has no
//     dA/dTleaf * dTleaf/dE * dE/dpsi. `use_energy_balance_` appears at six
//     sites in leaf_model.hpp and none of them is in the derivative.
//
// So `maximise_profit_over_collar` root-finds the first-order condition of a
// DIFFERENT MODEL from the one it reports. The tolerances below are the
// currently-observed values; the fix tightens them. They are asserted loosely on
// purpose so this lands green and the numbers are in the history.
//
// The oracle is a derivative-free scan of profit itself -- the method that
// arbitrated the last collar-solver change. It cannot inherit the
// derivative's defect, because it never evaluates a derivative.
// ============================================================================
void test_energy_balance_collar_solve_is_measured() {
  printf("MEASUREMENT: the EB collar solve against a derivative-free scan\n");

  double worst_dpsi = 0.0, worst_resid = 0.0, worst_dprofit = 0.0;
  int rows = 0, pinned_rows = 0;

  for (double ppfd : {500.0, 1500.0}) {
    for (double tair : {20.0, 30.0, 40.0}) {
      for (double vpd : {1.0, 2.0}) {
        Drivers d;
        d.PPFD = ppfd; d.leaf_temp = tair; d.atm_vpd = vpd;

        phylloptim::Leaf eb = make_pm_leaf(d, {2.0}, {1.0}, true);
        eb.find_root_collar_psi();
        const double psi_solver = eb.opt_root_psi_;
        const double profit_solver = eb.profit_;

        // The residual the solver believes it drove to zero. On the non-EB path
        // this is ~5.6e-15; here it is whatever staleness leaves.
        bool feasible = true;
        const double resid = eb.dprofit_droot_collar_psi(psi_solver, &feasible);
        if (!feasible) continue;

        // ORACLE: scan profit over the same feasible interval. A fresh leaf,
        // because dprofit_droot_collar_psi above has just re-seated the
        // temperature parameters at yet another point.
        phylloptim::Leaf scan = make_pm_leaf(d, {2.0}, {1.0}, true);
        double lo = 0.0, hi = 0.0;
        if (!scan.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(lo, hi)) continue;
        const int N = 20001;
        double best_psi = lo, best_profit = -std::numeric_limits<double>::max();
        for (int i = 0; i < N; ++i) {
          const double psi = lo + (hi - lo) * double(i) / double(N - 1);
          const double p = scan.profit_at_collar_psi<phylloptim::Leaf::CostCurve::TF24>(psi, lo, hi);
          if (p > best_profit) { best_profit = p; best_psi = psi; }
        }

        const double dpsi = std::abs(best_psi - psi_solver);
        const double dprofit = best_profit - profit_solver;
        worst_dpsi = std::max(worst_dpsi, dpsi);
        worst_dprofit = std::max(worst_dprofit, dprofit);

        // ⚠️ INTERIOR ROWS ONLY for the residual. A CONSTRAINED optimum sits on
        // a bracket bound, where dprofit is genuinely non-zero and "the gradient
        // should vanish" is simply the wrong statement -- the same distinction
        // test_collar_optimum_pinned_to_its_constraint makes on the non-EB path.
        // Lumping the two together reports a pinned row's honest gradient as if
        // it were a solver failure.
        const double span = std::max(hi - lo, 1e-12);
        const bool pinned = (psi_solver - lo) < 1e-3 * span ||
                            (hi - psi_solver) < 1e-3 * span;
        if (pinned) { ++pinned_rows; }
        else { worst_resid = std::max(worst_resid, std::abs(resid)); }
        ++rows;
      }
    }
  }

  printf("    %d feasible rows (%d pinned) | worst |dprofit|, interior only: %.3e\n",
         rows, pinned_rows, worst_resid);
  printf("    worst |psi_scan - psi_solver|: %.3e MPa | worst profit shortfall: %.3e\n",
         worst_dpsi, worst_dprofit);

  ok(rows > 0, "the EB grid has feasible rows to measure");

  // The bounds this test landed with, and what they replaced. Before the three
  // edits (seat the temperature parameters, add the dA/dTleaf chain term, handle
  // the compensation point) the same grid gave: residual 5.76, collar 0.83 MPa
  // from the argmax, 2.34 umol m^-2 s^-1 of profit left behind.
  ok(worst_resid < 1e-9,
     "interior EB optima satisfy dprofit == 0 to solver precision");
  // The scan resolves the argmax only to (hi-lo)/20000, ~5e-5 MPa on this grid,
  // so agreement AT that scale is agreement -- asking for more would be asking
  // the oracle for precision it does not have.
  ok(worst_dpsi < 1e-3,
     "the collar agrees with a derivative-free scan to the scan's resolution");
  // The scan can only ever find a profit >= the solver's, up to its own grid
  // resolution. A NEGATIVE shortfall beyond that would mean the scan is broken,
  // which is the one way this measurement could mislead.
  ok(worst_dprofit > -1e-6, "the scan never finds LESS profit than the solver");
  ok(worst_dprofit < 1e-6, "and never finds MORE: the solver is at the maximum");
}

// Gate-off inertness, proved by POISONING rather than by a recorded value. The
// golden file already pins values; what it cannot say is that the energy-balance
// inputs are never *read* when the gate is off. Setting them to NaN and demanding
// bit-identical results says exactly that: if any of the new code path touched
// Rn_, ra_ or Tair_ off-gate, a NaN would propagate and nothing would compare
// equal. Same shape as test_pm_wind_speed_validation uses for wind_speed_.
void test_energy_balance_gate_off_is_inert() {
  printf("gate off: the EB inputs are never read (NaN poisoning)\n");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (double tair : {20.0, 40.0}) {
    Drivers d;
    d.leaf_temp = tair;
    phylloptim::Leaf clean = make_pm_leaf(d, {2.0}, {1.0}, false);
    phylloptim::Leaf poisoned = make_pm_leaf(d, {2.0}, {1.0}, false);
    poisoned.Rn_ = nan; poisoned.ra_ = nan; poisoned.Tair_ = nan;

    clean.find_root_collar_psi();
    poisoned.find_root_collar_psi();
    const std::string at = " at Tair=" + std::to_string(int(tair));
    ok(clean.profit_ == poisoned.profit_, "profit is bit-identical" + at);
    ok(clean.assim_colimited_ == poisoned.assim_colimited_,
       "assimilation is bit-identical" + at);
    ok(clean.opt_root_psi_ == poisoned.opt_root_psi_,
       "the collar is bit-identical" + at);
    // And the derivative, which is where the new code lives.
    const double psi = clean.opt_root_psi_;
    ok(clean.dprofit_droot_collar_psi(psi) ==
           poisoned.dprofit_droot_collar_psi(psi),
       "dprofit is bit-identical" + at);
  }
}

// THE SCIENTIFIC ACCEPTANCE TEST. Jones et al. (2026, Global Change Biology
// 32:e70972) report that stomatal conductance keeps rising above the
// photosynthetic thermal optimum while assimilation falls, and that an
// optimality model reproduces it only when leaf temperature is inside the
// objective. This asserts that the corrected first-order condition produces that
// behaviour, and -- the half that matters -- that the gate-off arm does NOT, at
// the same drivers. Without the contrast the test could pass on the Arrhenius
// optimum alone, which has nothing to do with energy balance.
//
// ⚠️ THE SIGNATURE IS IN TRANSPIRATION, NOT IN CONDUCTANCE. Conductance rising
// where assimilation falls is what you see when stom_cond_CO2 divides by the
// prescribed AIR deficit however hot the leaf gets: gs is then E rescaled by a
// constant and inherits E's shape exactly. With the leaf-to-air deficit wired in
// the two come apart, because
//
//     gs = P*E / (1.6 * D_leaf),      D_leaf = atm_vpd + esat(Tleaf) - esat(Tair)
//
// and over the decoupled window D_leaf grows FASTER than E does. Measured on this
// sweep, between the thermal optimum and E's own peak: E x1.034, D_leaf x1.089,
// so gs x0.950 -- falling throughout. The gs counts went from a positive number
// to gate-on 0, gate-off 1.
//
// That is not a loss of the mechanism, and the assertions below are rewritten to
// say which is which:
//
//   * the ENERGY-BALANCE mechanism is intact and is what the gate buys. E keeps
//     rising above the thermal optimum, and it does so over a wider window with
//     the gate on (4 points) than off (1).
//   * the CONDUCTANCE signature, which is what a gas-exchange dataset reports,
//     does not survive the correction at these drivers. Anyone comparing this
//     model to a measured dgs/dT needs to know that, and the third assertion pins
//     the ratio that causes it rather than the count it produces.
//
// The sweep holds AIR VPD fixed, which is the protocol of the paper's Figure S2
// rather than its Figure 2 -- and note that fixing the air deficit no longer
// fixes the deficit the leaf actually sees.
void test_energy_balance_stomatal_decoupling() {
  printf("decoupling: E rises where A falls; gs does not, once D moves with Tleaf\n");
  // ⚠️ THE FIXTURE DEFAULTS CANNOT SHOW THIS, and the reason is physical rather
  // than numerical. kmax = K_s*theta/h is 3.14e-5 at the defaults, and a leaf
  // that conductive simply cannot move enough water to cool itself: measured in
  // R against Jones's own drivers, Tleaf - Tair floors at +1.9 K where their
  // leaf reaches -2.0 K, so the evaporative-cooling benefit never becomes large
  // enough to reverse the conductance response. Shortening the path length to
  // 1.5 m (kmax ~ 1.05e-4) puts the leaf in the range their measured species
  // occupy. This is choosing a regime where the mechanism is active, not tuning
  // until a test passes -- at the defaults the honest answer is "no decoupling",
  // and that is what the probe found.
  auto sweep = [&](bool gate) {
    std::vector<double> Ta, A, gs, E, D;
    for (double t = 15.0; t <= 45.0; t += 1.0) {
      Drivers d;
      d.PPFD = 1200.0; d.leaf_temp = t; d.atm_vpd = 1.0; d.h = 1.5;
      phylloptim::Leaf l = make_pm_leaf(d, {0.5}, {1.0}, gate);
      // ⚠️ OVERRIDE THE DERIVED RADIATION, and this is the second thing the
      // defaults cannot express. phylloptim sets Rn = 2*PPFD/4.57 - 40, which at
      // PPFD 1200 is 485 W m^-2 with net longwave fixed at -40 -- against an
      // isothermal value that runs -98 to -133 over this range. Combined with
      // ra from the wind model (31.6 s m^-1) the leaf runs ~18 K above air and
      // its thermal optimum falls off the bottom of any sensible sweep: the
      // first version of this test measured a curve that had already peaked
      // before 15 C. These are the values Jones et al. drive with (their eq. 25
      // and 29 at Iabs = 500, ra = 10), and both fields are settable precisely
      // so a caller can supply a better radiation budget than the minimal cut.
      l.Rn_ = 400.0;
      l.ra_ = 12.0;
      l.find_root_collar_psi();
      if (!std::isfinite(l.assim_colimited_) || !std::isfinite(l.stom_cond_CO2_)) continue;
      Ta.push_back(t); A.push_back(l.assim_colimited_); gs.push_back(l.stom_cond_CO2_);
      E.push_back(l.transpiration_); D.push_back(l.vpd_leaf_);
    }
    return std::make_tuple(Ta, A, gs, E, D);
  };

  // Points above the thermal optimum where the response variable RISES while
  // assimilation falls. `y` is E for the mechanism and gs for the signature.
  auto count_decoupled = [](const std::vector<double>& A,
                            const std::vector<double>& y) {
    if (A.size() < 5) return 0;
    const size_t peak = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
    int n = 0;
    for (size_t i = peak + 2; i + 1 < A.size(); ++i) {
      if ((A[i + 1] - A[i - 1]) < 0.0 && (y[i + 1] - y[i - 1]) > 0.0) ++n;
    }
    return n;
  };

  const auto on = sweep(true);
  const auto off = sweep(false);
  const auto& A_on = std::get<1>(on);
  const auto& gs_on = std::get<2>(on);
  const auto& E_on = std::get<3>(on);
  const auto& D_on = std::get<4>(on);

  const int nE_on = count_decoupled(A_on, E_on);
  const int nE_off = count_decoupled(std::get<1>(off), std::get<3>(off));
  const int ng_on = count_decoupled(A_on, gs_on);
  const int ng_off = count_decoupled(std::get<1>(off), std::get<2>(off));
  printf("    decoupled points -- in E: gate on %d, gate off %d;"
         " in gs: gate on %d, gate off %d\n", nE_on, nE_off, ng_on, ng_off);

  ok(A_on.size() > 5, "the decoupling sweep produced a curve");
  const size_t peak = std::distance(A_on.begin(),
                                    std::max_element(A_on.begin(), A_on.end()));
  ok(peak > 0 && peak + 1 < A_on.size(),
     "assimilation has an interior thermal optimum");
  ok(nE_on > 0, "with the gate ON, transpiration rises where assimilation falls");
  ok(nE_on > nE_off,
     "and it does so over a wider window than with the gate off -- "
     "the energy balance is why");

  // WHY the conductance signature does not follow, stated as the ratio rather
  // than as the count it happens to produce here. Over the window from the
  // thermal optimum to E's own peak, the deficit outgrows the flux, and
  // gs = P*E/(1.6*D) therefore falls even though E is rising.
  const size_t Epeak = std::distance(E_on.begin(),
                                     std::max_element(E_on.begin(), E_on.end()));
  ok(Epeak > peak, "transpiration peaks ABOVE the assimilation optimum");
  const double E_growth = E_on[Epeak] / E_on[peak];
  const double D_growth = D_on[Epeak] / D_on[peak];
  printf("    optimum -> E peak: E x%.4f, D_leaf x%.4f, gs x%.4f\n",
         E_growth, D_growth, gs_on[Epeak] / gs_on[peak]);
  ok(D_growth > E_growth,
     "the leaf-to-air deficit outgrows the flux, so gs falls where E rises");
}

// The closed-form fast path (leaf/closed_form.hpp). Two things matter: that it is
// actually fast, and that its error is characterised honestly rather than asserted
// to be small.
void test_closed_form() {
  printf("closed-form fast path\n");
  // Reference geometry, matching the companion analysis: collar held at zero,
  // single layer, kmax from height.
  const double eta = 12.0, eta_c = 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta);
  const double theta = 1.0 / 4669.0;
  // The root network is built ONCE, outside the closure. It does not depend on h
  // or vpd, and the timing loop below reports a "set_physiology" figure -- if the
  // architecture model were rebuilt per call that figure would be measuring the
  // caller's helper as much as the leaf's setter (measured: 0.106 -> 0.243 us,
  // i.e. more than half the reported cost would have been the helper).
  const std::vector<double> depth_1m{1.0};
  const phylloptim::RootNetwork net_1m = fixture::root_network({1.0}, depth_1m);
  const auto setp = [&](phylloptim::Leaf &l, double h, double vpd) {
    std::vector<double> ps{0.0}, dp{1.0};
    l.set_physiology(net_1m, 900.0, ps, dp, 1.0 * theta / (h * eta_c), vpd, 40.0,
                     25.0, 21.0, 101.3);
  };
  phylloptim::Leaf l;

  // Near the wet end, where the leading-order expansion is centred, it should be
  // very accurate.
  setp(l, 1.0, 2.0);
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
  const double A_wet = l.assim_colimited_;
  setp(l, 1.0, 2.0);
  // ⚠️ `solve` TAKES psi_upstream NOW, and the 0.0 is the reference geometry's
  // collar-at-zero rather than a placeholder. It used to take only the Newton
  // count, so `solve(l, 1)` compiled and meant something else entirely -- one
  // Newton step became a 1 MPa upstream potential. That is the whole reason the
  // upstream potential is an explicit argument instead of being read off the leaf.
  const phylloptim::closed_form::Solution wet =
      phylloptim::closed_form::solve(l, 0.0, 1);
  ok(std::abs(wet.assim / A_wet - 1.0) < 2e-3,
     "closed form is within 0.2% of the exact solve at h=1 m");
  ok(phylloptim::closed_form::within_guard(l, wet), "h=1 m passes the guard");

  // Error grows steeply as the leaf moves away from the wet end. These bounds
  // record measured behaviour -- they are deliberately loose enough to be stable
  // and tight enough to catch a regression.
  struct Case {
    double h, max_err;
  };
  for (const Case &cs : {Case{3.0, 2e-3}, Case{8.0, 1.5e-2}, Case{12.0, 4e-2}}) {
    setp(l, cs.h, 2.0);
    l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
    const double A_ex = l.assim_colimited_;
    setp(l, cs.h, 2.0);
    const double A_cf = phylloptim::closed_form::solve(l, 0.0, 1).assim;
    ok(std::abs(A_cf / A_ex - 1.0) < cs.max_err,
       "closed-form error is bounded at h=" + std::to_string(cs.h) + " m");
  }

  // The guard is coarse, and saying so is the point. It admits ~2% error in A at
  // h = 12 m (ci/ca ~ 0.55) and only rejects once the error is ~8% (h = 20 m).
  // So it is a filter on gross failure, not an error bound -- and note that the
  // heights it rejects are the dominant canopy trees.
  setp(l, 20.0, 2.0);
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
  const double A_tall = l.assim_colimited_;
  setp(l, 20.0, 2.0);
  const phylloptim::closed_form::Solution tall =
      phylloptim::closed_form::solve(l, 0.0, 1);
  ok(!phylloptim::closed_form::within_guard(l, tall), "h=20 m is rejected by the guard");
  ok(std::abs(tall.assim / A_tall - 1.0) > 3e-2,
     "and it is rejected because the error really is large there");

  // The TF24_beta2 = 1/c leaf, where xi is constant and nothing needs solving.
  phylloptim::Leaf exact_leaf(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.0 / 2.680147, 157.44, 0.30, 0.7, 0.99, 1e-3,
                        100, 1e-3, 1000, 7.5);
  ok(phylloptim::closed_form::beta2_is_exact(exact_leaf),
     "beta2_is_exact recognises TF24_beta2 = 1/stem_c");
  ok(!phylloptim::closed_form::beta2_is_exact(l), "and rejects the default TF24_beta2 = 1.5");
  setp(exact_leaf, 5.0, 1.5);
  exact_leaf.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
  const double A_ref = exact_leaf.assim_colimited_;
  setp(exact_leaf, 5.0, 1.5);
  const phylloptim::closed_form::Solution ex =
      phylloptim::closed_form::solve_exact_beta2(exact_leaf, 0.0);
  // ⚠️ IT USED TO REPORT NaN HERE, and this assertion used to pin that. A NaN
  // argmax is not usable once the method is wired into `optimise()`, so the
  // potential is recovered from the supply -- and the check is the round trip,
  // which is the only statement about it that does not restate the code.
  ok(std::isfinite(ex.psi_stem),
     "the explicit form now recovers a psi_stem rather than reporting NaN");
  ok(std::abs(exact_leaf.transpiration(ex.psi_stem, 0.0) / ex.transpiration -
              1.0) < 1e-6,
     "and it is the potential whose hydraulic supply carries that E");
  ok(std::abs(ex.assim / A_ref - 1.0) < 4e-2,
     "the explicit form is within a few percent of the exact solve");

  // Timing, reported rather than asserted: absolute microseconds are
  // machine-dependent, so a hard threshold would be a flaky test.
  const std::vector<double> hs{1, 2, 3, 5, 8, 12}, ds{0.8, 1.0, 1.5, 2.0};
  const int reps = 20000;
  double sink = 0;
  const auto time_it = [&](phylloptim::Leaf &leaf_ref, auto fn) {
    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < reps; ++r) {
      setp(leaf_ref, hs[r % 6], ds[r % 4]);
      sink += fn();
    }
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / reps;
  };
  const double t_setp = time_it(l, [&] { return l.ca_; });
  const double t_exact = time_it(l, [&] {
    l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
    return l.assim_colimited_;
  });
  const double t_cf =
      time_it(l, [&] { return phylloptim::closed_form::solve(l, 0.0, 1).assim; });
  const double t_expl = time_it(exact_leaf, [&] {
    return phylloptim::closed_form::solve_exact_beta2(exact_leaf, 0.0).assim;
  });
  printf("    set_physiology %.3f us | exact %.3f us | 1-Newton %.3f us (%.1fx) |"
         " explicit %.3f us (%.1fx)\n",
         t_setp, t_exact, t_cf, t_exact / t_cf, t_expl, t_exact / t_expl);
  ok(t_cf < t_exact, "the closed form is faster than the exact solve");
  ok(sink != 0.0, "timing loop was not optimised away");
}

// The carbon -> resistance map (root_network_from_carbon) is the one piece of
// root *architecture* left in this package; the supply solve itself only ever
// reads r_R_H_min and r_R_V_sum. Testing it directly is the point of having
// pulled it out of MultiLayerRoots -- and it is why the map stayed here rather
// than moving to plant, where the golden file could not reach it.
// The second supply path (issue #2 stage 3). Not wired into Leaf yet -- it exists
// so the concept in stage 2 has two real alternatives to dispatch between, and so
// the dispatch measurement was made against a genuine second type rather than a
// stub the optimiser could see through.
// A whole Leaf solving through SinglePotential (issue #2 stage 2). This is the
// point of the whole item: the gas-exchange core is supply-agnostic, so swapping
// the supply path should change the operating point and nothing else.
void test_leaf_on_single_potential() {
  printf("Leaf solving on the single-potential supply path\n");
  Drivers d;

  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.set_supply_single();

  // ⚠️ THE SAME CALL AS THE MULTI-LAYER PATH, which is the point. The resistance
  // is a per-call driver on both paths now, carried by the RootNetwork argument;
  // it used to be a set_supply_single() constructor-style argument. 1.0e3 is per
  // unit leaf area, i.e. the old 2.0e4 * 0.05.
  std::vector<double> psi_soil{1.0}, depth{1.0};
  l.set_physiology(fixture::series_resistance(1.0e3), d.PPFD, psi_soil, depth,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  l.find_root_collar_psi();

  ok(std::isfinite(l.profit_), "single-potential solve gives a finite profit");
  ok(std::isfinite(l.opt_psi_stem_), "and a finite stem potential");
  ok(l.opt_psi_stem_ > 0.0 && l.opt_psi_stem_ <= l.psi_crit,
     "stem potential is a positive magnitude within psi_crit");
  ok(l.opt_root_psi_ >= 0.0, "collar potential is stored as a positive magnitude");
  ok(l.assim_colimited_ > 0.0, "the leaf assimilates");
  ok(l.soil_consumption_.size() == 1u,
     "the consumption buffer is sized to one layer, not the caller's vector");

  // The collar must sit between the soil and the stem: water runs downhill.
  const double collar_mag = l.opt_root_psi_;
  ok(collar_mag >= psi_soil[0] - 1e-9 && collar_mag <= l.opt_psi_stem_ + 1e-9,
     "collar potential lies between soil and stem");

  // Drier soil must cost carbon here too -- the same contract the multi-layer
  // path is held to, which is what makes the two comparable at all.
  phylloptim::Leaf dry;
  dry.setup_transpiration(100);
  dry.setup_root_vulnerability(100);
  dry.set_supply_single();
  std::vector<double> psi_dry{3.0};
  dry.set_physiology(fixture::series_resistance(1.0e3), d.PPFD, psi_dry, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  dry.find_root_collar_psi();
  ok(dry.profit_ < l.profit_, "drier soil yields less profit");

  // A larger series resistance is a worse-supplied plant, so it must not do
  // better. This is the knob the multi-layer path spends root carbon to lower.
  phylloptim::Leaf tight;
  tight.setup_transpiration(100);
  tight.setup_root_vulnerability(100);
  tight.set_supply_single();
  tight.set_physiology(fixture::series_resistance(1.0e4), d.PPFD, psi_soil, depth,
                       d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                       d.atm_o2_kpa, d.atm_kpa);
  tight.find_root_collar_psi();
  ok(tight.profit_ <= l.profit_, "a higher series resistance does not help");

  // And the default is unchanged: a Leaf nobody configures is multi-layer.
  phylloptim::Leaf plain;
  ok(plain.supply_kind_ == phylloptim::Leaf::SupplyKind::MultiLayer,
     "the supply path defaults to multi-layer");
}

void test_single_potential() {
  printf("single-potential supply path\n");
  phylloptim::SinglePotential sp;
  sp.set_soil_state(1.5);        // positive magnitude, -MPa
  // resistance_ is PER UNIT LEAF AREA, like every other input to the leaf.
  sp.resistance_ = 1.0e3;

  // begin_solve reports the only suction; there is nothing to flip (#25).
  near(sp.begin_solve(), 1.5, 1e-14, "begin_solve returns the soil suction");
  ok(sp.n_layers() == 1, "single potential writes exactly one layer");

  // Ohm's law, and the sign that matters: a collar drier than the soil -- a
  // LARGER suction now -- draws water UP (positive uptake).
  std::vector<double> consumption(1, 0.0);
  double E_up = 0.0;
  sp.uptake(2.5, consumption, E_up);
  ok(E_up > 0.0, "a collar drier than the soil draws water up");
  near(E_up, (2.5 - 1.5) / sp.resistance_ * phylloptim::kg_per_mol_h2o,
       1e-14, "uptake is the Ohm's-law flux");
  ok(consumption[0] > 0.0, "per-layer consumption is filled");

  // A collar WETTER than the soil pushes water back into it. Losing this sign is
  // how hydraulic redistribution silently becomes extra uptake.
  sp.uptake(0.5, consumption, E_up);
  ok(E_up < 0.0, "a collar wetter than the soil loses water to it");

  // The analytic derivative must match a central difference on uptake, and stay
  // finite everywhere -- unlike MultiLayerRoots there are no branch kinks, so it
  // never asks the caller for a finite-difference fallback.
  const double h = 1e-6, p0 = 2.5;
  double up = 0.0, dn = 0.0;
  sp.uptake(p0 + h, consumption, up);
  sp.uptake(p0 - h, consumption, dn);
  const double fd = (up - dn) / (2.0 * h);
  near(sp.duptake_dpsi(), fd, 1e-8, "analytic duptake_dpsi matches FD");
  ok(sp.duptake_dpsi() > 0.0,
     "duptake_dpsi is a positive conductance: uptake rises as the collar pulls harder");

  // An UNSET resistance would be an infinite or NaN flux; it is rejected, not
  // returned. The default is now the NA sentinel rather than zero, because the
  // resistance became a per-call driver -- an unset one means set_physiology was
  // never called, which is the same class of mistake as an unset psi_soil.
  phylloptim::SinglePotential bad;
  bad.set_soil_state(1.0);
  bad.begin_solve();
  bool threw = false;
  try {
    bad.uptake(2.0, consumption, E_up);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "an unset resistance throws rather than returning an infinity");

  // --- the driver entry point, and its two guards ---------------------------
  // set_supply_resistances is what makes the two supply paths take the same
  // set_physiology argument, so what it accepts and refuses is the contract.
  phylloptim::SinglePotential drv;
  drv.set_supply_resistances(fixture::series_resistance(2.5e3));
  near(drv.resistance_, 2.5e3, 1e-14,
       "the series resistance is read from r_R_V_sum[0]");

  // A network for the OTHER path carries a vulnerability-weighted horizontal
  // term this path cannot apply. Ignoring it would silently drop a resistance the
  // caller meant to use, so it is refused.
  threw = false;
  try {
    drv.set_supply_resistances(fixture::root_network({20.0}, {1.0}));
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a multi-layer network is refused, not silently reinterpreted");
  near(drv.resistance_, 2.5e3, 1e-14,
       "and the refusal leaves the previous resistance intact");

  // More than one layer is the other path's shape too.
  threw = false;
  try {
    phylloptim::RootNetwork two;
    two.r_R_V_sum.assign(2, 1.0e3);
    drv.set_supply_resistances(two);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "two layers are refused on a one-potential path");

  // And a non-positive resistance is caught at the boundary rather than inside a
  // root-find, which is why this entry point validates at all.
  threw = false;
  try {
    drv.set_supply_resistances(fixture::series_resistance(0.0));
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a zero series resistance is refused at the boundary");
}

void test_root_network_from_carbon() {
  printf("root architecture: carbon -> resistance\n");
  const double beta_H = 3.4e2, beta_V = 9.4e3, dz = 0.5;
  const std::vector<double> dz3(3, dz), dz1(1, dz);

  // Closed form, straight from the documented model: carbon splits 1/3 vertical
  // : 2/3 horizontal, r_R_H_min = beta_H/c_r_h, r_R_V = beta_V*dz[i]^2/c_r_v.
  const std::vector<double> carbon{3.0, 6.0, 1.5};
  const auto n = phylloptim::root_network_from_carbon(carbon, dz3, beta_H, beta_V);

  ok(n.r_R_H_min.size() == 3u, "one resistance per rooted layer");
  near(n.r_R_H_min[0], beta_H / (3.0 * 2.0 / 3.0), 1e-14, "r_R_H_min layer 0");
  near(n.r_R_H_min[1], beta_H / (6.0 * 2.0 / 3.0), 1e-14, "r_R_H_min layer 1");
  near(n.r_R_V[2], beta_V * dz * dz / (1.5 / 3.0), 1e-14, "r_R_V layer 2");

  // r_R_V_sum is a running total down the profile, so it is monotone and its
  // last entry is the whole column's vertical resistance.
  ok(n.r_R_V_sum[0] < n.r_R_V_sum[1] && n.r_R_V_sum[1] < n.r_R_V_sum[2],
     "cumulative vertical resistance increases with depth");
  near(n.r_R_V_sum[2], n.r_R_V[0] + n.r_R_V[1] + n.r_R_V[2], 1e-14,
       "r_R_V_sum is the running sum of r_R_V");

  // More carbon is less resistance, in both directions. This is the sign that
  // matters: getting it backwards would make investment in roots harmful.
  const auto rich = phylloptim::root_network_from_carbon({12.0}, dz1, beta_H, beta_V);
  const auto poor = phylloptim::root_network_from_carbon({3.0}, dz1, beta_H, beta_V);
  ok(rich.r_R_H_min[0] < poor.r_R_H_min[0], "more root carbon -> less horizontal resistance");
  ok(rich.r_R_V_sum[0] < poor.r_R_V_sum[0], "more root carbon -> less vertical resistance");

  // Trailing zero-carbon layers are dropped, so the hot loop never visits them.
  const auto trailing = phylloptim::root_network_from_carbon({3.0, 6.0, 0.0, 0.0},
                                                      std::vector<double>(4, dz),
                                                      beta_H, beta_V);
  ok(trailing.r_R_H_min.size() == 2u, "trailing rootless layers are dropped");

  // Negative carbon is rejected rather than producing a negative resistance.
  bool threw = false;
  try {
    phylloptim::root_network_from_carbon({3.0, -1.0}, std::vector<double>(2, dz),
                                        beta_H, beta_V);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "negative root carbon throws");

  // -------------------------------------------------------------------------
  // Per-layer thickness (#626)
  // -------------------------------------------------------------------------
  // Each layer's vertical resistance uses ITS OWN thickness, not a column
  // average. Asserted per layer, because a scalar dz would satisfy the r_R_V_sum
  // total on a uniform profile and only differ on a graded one.
  const std::vector<double> graded{0.02, 0.28, 1.2};
  const auto g = phylloptim::root_network_from_carbon(carbon, graded, beta_H, beta_V);
  for (int i = 0; i < 3; ++i) {
    near(g.r_R_V[i], beta_V * graded[i] * graded[i] / (carbon[i] / 3.0), 1e-14,
         "r_R_V uses this layer's own thickness");
  }
  // The horizontal term has no thickness in it, so it must be untouched by the
  // profile -- which is what says the change went where it was meant to. Exact,
  // because these are the same expression on the same input.
  for (int i = 0; i < 3; ++i) {
    ok(g.r_R_H_min[i] == n.r_R_H_min[i], "r_R_H_min does not depend on dz");
  }

  // DISCRETISATION INVARIANCE, the property the scalar dz broke and the reason
  // this is a bug fix and not a generalisation. With root density uniform over
  // depth the layer-integrated carbon goes as dz[i], so total vertical
  // resistance is 3*beta_V*D^2/C for ANY slicing of the same column. Three
  // profiles of total depth 1.5 m carrying the same 20 kg C m^-2:
  const double kD = 1.5, kC = 20.0;
  const std::vector<std::vector<double> > profiles{
      {0.3, 0.3, 0.3, 0.3, 0.3},              // 5 equal layers
      {0.5, 0.5, 0.5},                        // 3 equal layers
      {0.02, 0.28, 0.3, 0.4, 0.5},            // a 2 cm evaporation layer on top
      {0.75, 0.25, 0.25, 0.15, 0.1}           // thick over thin
  };
  const double expected = 3.0 * beta_V * kD * kD / kC;
  for (const auto &p : profiles) {
    std::vector<double> c(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      c[i] = kC * p[i] / kD;  // uniform density -> carbon proportional to width
    }
    const auto net = phylloptim::root_network_from_carbon(c, p, beta_H, beta_V);
    near(net.r_R_V_sum.back(), expected, 1e-12,
         "total vertical resistance is independent of the layer slicing");
  }

  // A width vector of the wrong length, or a non-positive width, is refused
  // rather than read past the end or turned into a zero resistance. Checked by
  // REASON, because the two refusals are different mistakes: a length mismatch is
  // a mis-paired pair of vectors, a non-positive width is a bad profile.
  ok(throws_with([&] { phylloptim::root_network_from_carbon(
                           carbon, std::vector<double>{0.5, 0.5},
                           beta_H, beta_V); },
                 "same number of elements"),
     "a dz vector of the wrong length is refused as a length mismatch");
  for (double w : {0.0, -0.1}) {
    ok(throws_with([&] { phylloptim::root_network_from_carbon(
                             carbon, std::vector<double>{0.5, w, 0.5},
                             beta_H, beta_V); },
                   "finite and positive"),
       "a non-positive dz is refused as a bad width");
  }
}

// layer_thicknesses(): the profile -> widths helper an R caller and the test
// fixtures come through (#626).
void test_layer_thicknesses() {
  printf("root architecture: layer_thicknesses\n");

  // Differences with an implicit 0 at the surface, so the first width is the
  // first boundary and equal layers come back equal.
  const auto uniform = phylloptim::layer_thicknesses({1.0, 2.0, 3.0});
  for (int i = 0; i < 3; ++i) {
    // Exact, and that is the load-bearing part: it is what keeps the golden grid
    // (1 m layers, differenced here since #626) bit-identical.
    ok(uniform[i] == 1.0, "1 m boundaries give exactly 1 m widths");
  }

  const auto graded = phylloptim::layer_thicknesses({0.02, 0.3, 1.5});
  near(graded[0], 0.02, 1e-16, "surface layer width is the first boundary");
  near(graded[1], 0.28, 1e-16, "interior width is a difference of boundaries");
  near(graded[2], 1.2, 1e-16, "deepest width is a difference of boundaries");

  // ⚠️ THE NEGATIVE, which is why plant passes its own widths rather than coming
  // through here: differencing is NOT bit-exact. 1.5 m over 5 layers is the
  // package's own default geometry, and the round trip loses a bit.
  const size_t kN = 5;
  const double kDelta = 1.5 / kN;
  std::vector<double> profile(kN);
  for (size_t i = 0; i < kN; ++i) {
    profile[i] = (i + 1) * kDelta;
  }
  const auto round_tripped = phylloptim::layer_thicknesses(profile);
  bool all_exact = true;
  for (size_t i = 0; i < kN; ++i) {
    if (round_tripped[i] != kDelta) {
      all_exact = false;
    }
  }
  ok(!all_exact,
     "differencing a uniform profile is NOT bit-exact -- pass widths, do not "
     "re-derive them");
  for (size_t i = 0; i < kN; ++i) {
    near(round_tripped[i], kDelta, 1e-15, "...but it is right to rounding");
  }

  // A profile that is not a profile: empty, or not strictly increasing. The
  // second is the one that matters -- it is what happens when someone passes
  // WIDTHS to a function wanting boundaries, which for a graded profile is not
  // even monotone.
  ok(throws_with([] { phylloptim::layer_thicknesses({}); }, "at least one layer"),
     "an empty profile is refused");
  for (const std::vector<double> &bad :
       {std::vector<double>{0.3, 0.3},        // repeated boundary
        std::vector<double>{0.5, 0.4, 0.3},   // decreasing
        std::vector<double>{-0.5},            // above the surface
        std::vector<double>{0.0, 1.0}}) {     // a zero-thickness surface layer
    ok(throws_with([&] { phylloptim::layer_thicknesses(bad); },
                   "strictly increasing"),
       "a profile that is not strictly increasing is refused, and says so");
  }

  // ⚠️ AND THE LIMIT OF THAT GUARD, ASSERTED SO IT IS NOT MISREAD AS PROTECTION.
  // The mistake it was reached for -- passing WIDTHS where boundaries are wanted
  // -- is only caught when the widths happen to be non-monotone. plant's own
  // graded profile is a counterexample: its widths increase with depth, so read
  // as boundaries they pass every check here and silently yield different,
  // entirely plausible widths. Only the caller knows which it holds; the guard
  // catches a malformed profile, not a mislabelled one.
  const std::vector<double> widths{0.02, 0.28, 0.3, 0.4, 0.5};
  const auto misread = phylloptim::layer_thicknesses(widths);
  bool differs = false;
  for (size_t i = 0; i < widths.size(); ++i) {
    if (misread[i] != widths[i]) {
      differs = true;
    }
  }
  ok(differs, "widths misread as boundaries are accepted and silently wrong");
}

// The temperature-response parameters were constexpr constants; they are now
// settable members. This checks the point of that change -- that setting them
// actually moves the model -- rather than just that they compile.
void test_temperature_parameters_are_settable() {
  printf("temperature-response parameters are settable\n");
  Drivers d;
  phylloptim::Leaf base = make_leaf(d, {2.0}, {1.0});
  base.find_root_collar_psi();
  const double A_base = base.assim_colimited_;

  // Defaults must equal the published constants, so this is a no-op refactor for
  // anyone who does not touch them.
  ok(base.vcmax_ha_ == phylloptim::vcmax_ha, "vcmax_ha_ defaults to the constant");
  ok(base.jmax_d_S_ == phylloptim::jmax_d_S, "jmax_d_S_ defaults to the constant");
  ok(base.gamma_25_ == phylloptim::gamma_25, "gamma_25_ defaults to the constant");
  near(base.R_d_25, 1.44, 1e-12, "R_d_25 default");

  // Raising the Vcmax activation energy raises Vcmax above the 25 C reference,
  // so a 25 C leaf should be unaffected but a warm one should assimilate more.
  {
    phylloptim::Leaf warm = make_leaf(d, {2.0}, {1.0});
    phylloptim::Leaf warm_hi = make_leaf(d, {2.0}, {1.0});
    warm_hi.vcmax_ha_ = phylloptim::vcmax_ha * 1.5;
    // set_physiology already ran, so push the change through the T-response block.
    warm.update_temperature_dependent_params(35.0);
    warm_hi.update_temperature_dependent_params(35.0);
    ok(warm_hi.vcmax_ > warm.vcmax_,
       "a larger activation energy gives a larger vcmax at 35 C");
  }

  // Respiration: doubling it must lower assimilation.
  {
    phylloptim::Leaf r2 = make_leaf(d, {2.0}, {1.0});
    r2.R_d_25 = 2.0 * r2.R_d_25;
    r2.update_temperature_dependent_params(d.leaf_temp);
    r2.find_root_collar_psi();
    ok(r2.assim_colimited_ < A_base, "doubling R_d_25 lowers assimilation");
    ok(r2.R_d_ > base.R_d_, "and raises R_d");
  }

  // The CO2 compensation point feeds photorespiration, so raising it lowers A.
  {
    phylloptim::Leaf g2 = make_leaf(d, {2.0}, {1.0});
    g2.gamma_25_ = phylloptim::gamma_25 * 1.5;
    g2.update_temperature_dependent_params(d.leaf_temp);
    g2.find_root_collar_psi();
    ok(g2.assim_colimited_ < A_base,
       "raising the compensation point lowers assimilation");
  }
}

// ⚠️ THE ASSERTIONS ABOVE ALL PUSH THEIR CHANGE THROUGH BY CALLING
// update_temperature_dependent_params() DIRECTLY, AND THAT IS WHY THEY PASSED
// WHILE THE FEATURE WAS BROKEN FROM R. A caller does not have that route: they set
// the field and then set the drivers, which is the one path that took a cache hit
// and silently kept the old response. The test worked around the bug it should
// have caught.
//
// So this asserts the REALISTIC path, twice over: setting a temperature-response
// parameter on an already-solved leaf and re-supplying the SAME drivers must change
// the answer, and must land on exactly what a freshly constructed leaf gives.
//
// Bit-exactness is the right bar for the second half -- the two routes share no
// code, so anything the cache fails to invalidate shows up as a difference. #41.
void test_temperature_params_invalidate_cache() {
  printf("setting a temperature parameter invalidates the cache\n");
  Drivers d;

  // The realistic route: solve, change the parameter, re-supply the same drivers.
  phylloptim::Leaf warm = make_leaf(d, {2.0}, {1.0});
  warm.find_root_collar_psi();
  const double A_before = warm.assim_colimited_;
  const double Rd_before = warm.R_d_;

  warm.R_d_25 = 2.0 * warm.R_d_25;
  warm.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                      {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                      d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  warm.find_root_collar_psi();

  ok(warm.R_d_ != Rd_before,
     "re-supplying the same drivers after a parameter change recomputes R_d");
  ok(warm.assim_colimited_ != A_before,
     "and the operating point moves");
  near(warm.R_d_, 2.0 * Rd_before, 1e-12,
       "doubling R_d_25 doubles R_d at the same temperature");

  // And it must agree bit-for-bit with never having had the stale value.
  phylloptim::Leaf fresh = make_leaf(d, {2.0}, {1.0});
  fresh.R_d_25 = 2.0 * fresh.R_d_25;
  fresh.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                       {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                       d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  fresh.find_root_collar_psi();
  ok(warm.assim_colimited_ == fresh.assim_colimited_,
     "a re-parameterised leaf is bit-identical to one built with the value");
  ok(warm.R_d_ == fresh.R_d_, "R_d likewise");

  // ⚠️ The same hole covered vcmax_25, which is a TRAIT. set_traits() remains the
  // correct way to change one -- the vulnerability splines and the solved point
  // need clearing too -- but the silent-wrong-number half of hazard 10 is gone.
  phylloptim::Leaf vc = make_leaf(d, {2.0}, {1.0});
  vc.find_root_collar_psi();
  const double vcmax_before = vc.vcmax_;
  vc.vcmax_25 = vc.vcmax_25 * 1.5;
  vc.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                    {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                    d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  ok(vc.vcmax_ > vcmax_before,
     "a bare vcmax_25 write no longer leaves vcmax_ describing the old value");

  // The cache must still BE a cache: identical inputs twice must not recompute
  // into a different answer.
  phylloptim::Leaf same = make_leaf(d, {2.0}, {1.0});
  same.find_root_collar_psi();
  const double A_once = same.assim_colimited_;
  same.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                      {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                      d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  same.find_root_collar_psi();
  ok(same.assim_colimited_ == A_once,
     "unchanged inputs still give a bit-identical answer");
}

// R_d's TEMPERATURE RESPONSE (#41).
//
// ⚠️ THE GOLDEN FILE IS NEARLY BLIND TO THIS, and that is the reason this test
// exists. Every reference value in the model is DEFINED at 25 C, so a change to any
// response curve is inert there by construction; the golden grid carries one hot
// block for exactly this reason, and this test is what pins the response itself.
// The pack-reading kernels against the member-reading ones they are twins of,
// at double, where the two must be THE SAME NUMBERS and not merely close.
//
// This is the only check on the reverse-mode surface's VALUES. The transpose
// identity reads that surface in two directions and is satisfied exactly by a
// kernel evaluated at the wrong temperature or in the wrong light, because both
// directions are then wrong together.
//
// ⚠️ THE ENERGY-BALANCE ARM IS THE POINT, not a second case for completeness.
// Off that path the block's temperature and `leaf_temp_` are one number and a
// kernel reading either passes; on it they run up to 10 K apart, and reading
// `leaf_temp_` moved assimilation by 100%, electron transport by 84% and dA/dci
// by 8.3%. The block members are compared FIRST so a failure names the cause
// rather than the three consequences below it.
void test_pack_kernels_are_the_models_own() {
  printf("the pack-reading kernels are the model's own numbers\n");
  using AD = xad::fwd<double>::active_type;
  int compared = 0, agreed = 0, refused = 0;
  double widest_gap = 0.0;
  // Identical numbers, and identical refusals. At a leaf 50 C hot the
  // colimitation discriminant leaves dA/dci NaN, and two NaNs from two spellings
  // of one function are the agreement being asserted rather than an exception to
  // it -- `near` cannot say that, since NaN compares false against everything.
  // Counted, so a fixture that started refusing everywhere cannot pass by it.
  auto same = [&](double got, double want, const std::string &what) {
    if (std::isnan(got) && std::isnan(want)) {
      ++refused;
      ++checks;
      return;
    }
    ++agreed;
    near(got, want, 0.0, what);
  };
  for (bool gate : {false, true}) {
    for (double psi0 : {0.5, 1.0, 2.0, 3.0}) {
      for (double ppfd : {100.0, 1500.0}) {
        for (double temp : {25.0, 40.0}) {
          Drivers d;
          d.PPFD = ppfd;
          d.leaf_temp = temp;
          phylloptim::Leaf l =
              make_pm_leaf(d, {psi0, psi0 + 0.35, psi0 + 0.7}, {1.0, 2.0, 3.0},
                           gate);
          l.find_root_collar_psi();
          const double ci = l.ci_;
          if (!std::isfinite(ci)) continue;
          ++compared;
          widest_gap = std::max(widest_gap, std::abs(l.Tleaf_ - l.leaf_temp_));
          const phylloptim::leaf_pars<double> p = l.passive_pars();
          const std::string at = " at psi_soil=" + std::to_string(psi0) +
                                 " ppfd=" + std::to_string(ppfd) + " Tair=" +
                                 std::to_string(temp) +
                                 (gate ? " eb=1" : " eb=0");
          same(l.vcmax_at<double>(p), l.vcmax_, "vcmax" + at);
          same(l.jmax_at<double>(p), l.jmax_, "jmax" + at);
          same(l.respiration_at<double>(p), l.R_d_, "respiration" + at);
          same(l.electron_transport_at<double>(p), l.electron_transport(),
               "electron transport" + at);
          same(l.assim_colimited_kernel<double>(ci, p), l.assim_colimited(ci),
               "colimited assimilation" + at);
          AD ci_ad = ci;
          xad::derivative(ci_ad) = 1.0;
          same(l.assim_slope_at<double>(ci, p),
               xad::derivative(l.assim_colimited_kernel(ci_ad)),
               "dA/dci" + at);
        }
      }
    }
  }
  ok(compared >= 20, "the sweep reached at least 20 operating points");
  ok(agreed > 4 * refused,
     "most comparisons were of numbers, not of refusals: " +
         std::to_string(agreed) + " numbers, " + std::to_string(refused) +
         " refusals");
  // The claim is only evidence where the two temperatures actually differ, so
  // a fixture that stopped separating them would have stopped checking.
  ok(widest_gap > 1.0,
     "the sweep reached a leaf more than 1 K off air temperature, gap=" +
         std::to_string(widest_gap));
}

// Light reaches carbon with a ROW, and the row is the right size.
//
// ⚠️ AN EXACT ZERO IS THE FAILURE THIS EXISTS FOR, not a small number. The pack
// slot is the only route light has into the electron transport; read off the
// `PPFD_` member instead, every trait that moves the light a cohort stands in --
// the extinction coefficient, the leaf area a stem carries -- reaches
// assimilation with no row and the gradient reads zero rather than wrong. The
// transpose identity is satisfied exactly by a zero row on both sides, so it
// cannot see this either.
//
// The difference is taken at a FIXED ci, so it differences the kernel and not the
// solve, and the row is required to FALL with light: at 100 umol the leaf is
// light-limited and at 1500 it is rubisco-limited, so a row that did not saturate
// would be reading light somewhere it does not belong.
void test_light_reaches_carbon_with_a_row() {
  printf("light reaches carbon with a row\n");
  using AD = xad::fwd<double>::active_type;
  double previous = std::numeric_limits<double>::infinity();
  int compared = 0;
  for (double ppfd : {100.0, 400.0, 900.0, 1500.0}) {
    Drivers d;
    d.PPFD = ppfd;
    phylloptim::Leaf l = make_leaf(d, {1.0}, {1.0});
    l.find_root_collar_psi();
    const double ci = l.ci_;
    if (!std::isfinite(ci)) continue;
    ++compared;
    const std::string at = " at ppfd=" + std::to_string(ppfd);

    phylloptim::leaf_pars<AD> p;
    const phylloptim::leaf_pars<double> seated = l.passive_pars();
    for (std::size_t i = 0; i < p.size(); ++i) p[i] = AD(seated[i]);
    xad::derivative(p[phylloptim::par_PPFD]) = 1.0;
    const double row =
        xad::derivative(l.assim_colimited_kernel<AD>(AD(ci), p));

    const double h = ppfd * 1e-6;
    Drivers up_d = d, dn_d = d;
    up_d.PPFD = ppfd + h;
    dn_d.PPFD = ppfd - h;
    phylloptim::Leaf up = make_leaf(up_d, {1.0}, {1.0});
    phylloptim::Leaf dn = make_leaf(dn_d, {1.0}, {1.0});
    const double fd = (up.assim_colimited(ci) - dn.assim_colimited(ci)) / (2 * h);

    ok(row != 0.0, "dA/dPPFD is not an exact zero" + at);
    ok(row > 0.0, "more light cannot lower assimilation" + at);
    near(row, fd, 1e-7, "dA/dPPFD against a central difference" + at);
    ok(row < previous, "the light response saturates" + at);
    previous = row;
  }
  ok(compared == 4, "every light level solved");
}

// Every operating point the model answers at hands over FINITE rows.
//
// ⚠️ A NON-FINITE ROW IS NOT REPORTED ANYWHERE DOWNSTREAM, which is why this is
// asserted here. plant records its rows against a constant 1.0, so the row it
// validates is always finite and the sentinel rides in on the tape behind it;
// its census gradient then polls for a DECLARED refusal, finds none, and hands
// back NaN as a number the sweep computed. Measured before this check existed:
// three of five census drivers returned NaN in every trait column of three
// metrics, with nothing refused.
//
// The two routes are separate claims. `duptake_dpsi` returns NaN deliberately
// where its analytic branch is not valid, and a single-layer wet bound -- where
// the shade-death exit places the collar -- is such a point; the
// coordinates additionally run the implicit function theorem, which has no
// residual to invert wherever the model assigned ci rather than solving for it.
void test_every_answered_point_hands_over_finite_rows() {
  printf("every answered operating point hands over finite rows\n");
  using T = odelia::ode::tangent_scalar<double>;
  std::map<std::string, int> reached;
  int draws = 0, uptakes = 0, coords = 0;
  for (int layers : {1, 3}) {
    for (double psi0 : {0.25, 1.0, 2.0, 4.0, 6.0}) {
      for (double ppfd : {1.0, 100.0, 900.0}) {
        Drivers d;
        d.PPFD = ppfd;
        std::vector<double> psi_soil, soil_depth;
        for (int i = 0; i < layers; ++i) {
          psi_soil.push_back(psi0 + 0.35 * i);
          soil_depth.push_back(1.0 * (i + 1));
        }
        phylloptim::Leaf l = make_leaf(d, psi_soil, soil_depth);
        l.find_root_collar_psi();
        const auto kind = l.operating_point_kind();
        const std::string name = l.operating_point_kind_name(kind);
        if (kind == phylloptim::Leaf::OperatingPointKind::Unsolved ||
            kind == phylloptim::Leaf::OperatingPointKind::SolverRefused ||
            kind == phylloptim::Leaf::OperatingPointKind::NonFiniteGradient) {
          continue;
        }
        ++reached[name];
        const std::string at = " at " + name + " psi_soil=" +
                               std::to_string(psi0) + " ppfd=" +
                               std::to_string(ppfd) + " layers=" +
                               std::to_string(layers);

        phylloptim::leaf_pars<T> pars;
        const phylloptim::leaf_pars<double> seated = l.passive_pars();
        for (std::size_t i = 0; i < pars.size(); ++i) pars[i] = T(seated[i]);
        odelia::ode::seed_direction(pars[phylloptim::par_vcmax_25], 1.0);
        std::vector<T> soil, r_h, r_v;
        for (int i = 0; i < l.supply_n_layers(); ++i) {
          soil.push_back(T(l.roots_.psi_soil_[std::size_t(i)]));
          r_h.push_back(T(l.roots_.network_.r_R_H_min[std::size_t(i)]));
          r_v.push_back(T(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
        }
        const phylloptim::SupplyAt<T> supply{soil, r_h, r_v,
                                             pars[phylloptim::par_root_P50],
                                             pars[phylloptim::par_root_c]};
        const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), supply);
        ++draws;
        ok(std::isfinite(odelia::util::to_passive(draw.flux.slope)),
           "the draw's collar slope is finite" + at);
        for (std::size_t j = 0; j < draw.duptake_dp.size(); ++j) {
          ok(std::isfinite(draw.duptake_dp[j]),
             "layer " + std::to_string(j) + "'s supplied collar slope is finite" + at);
        }

        const double curvature =
            kind == phylloptim::Leaf::OperatingPointKind::Interior
                ? l.marginal_collar_slope<phylloptim::Leaf::CostCurve::TF24>()
                : std::numeric_limits<double>::quiet_NaN();
        if (kind == phylloptim::Leaf::OperatingPointKind::Interior &&
            !std::isfinite(curvature)) {
          continue;
        }
        const T collar =
            l.collar_at<phylloptim::Leaf::CostCurve::TF24, T>(draw, pars, curvature);
        const auto out =
            l.outputs_at<phylloptim::Leaf::CostCurve::TF24, T>(collar, draw, pars);
        ++uptakes;
        ok(std::isfinite(odelia::ode::derivative_along(out.profit)),
           "profit's row is finite" + at);
        for (std::size_t j = 0; j < out.uptake.size(); ++j) {
          ok(std::isfinite(odelia::ode::derivative_along(out.uptake[j])),
             "layer " + std::to_string(j) + "'s uptake row is finite" + at);
        }
        // The coordinates are refused by name at a shutdown; outputs_at answers
        // there instead, which the branch above just exercised.
        if (kind != phylloptim::Leaf::OperatingPointKind::HydraulicShutdown) {
          const auto co = l.collar_coords_at<T>(l.opt_psi_stem_, l.ci_, collar,
                                                draw, pars, true);
          ++coords;
          ok(std::isfinite(odelia::ode::derivative_along(co.sigma.value)) &&
                 std::isfinite(odelia::util::to_passive(co.sigma.slope)),
             "sigma's row and collar slope are finite" + at);
          ok(std::isfinite(odelia::ode::derivative_along(co.ci.value)) &&
                 std::isfinite(odelia::util::to_passive(co.ci.slope)),
             "ci's row and collar slope are finite" + at);
          // Omitting the channel is what profit asks for at an interior point,
          // and it must move no number: the step it multiplies is exactly zero
          // in value, so both coordinates are the residuals' own bit for bit.
          // Asserted here because the value is where an omission would show and
          // the row is not -- a channel that went missing leaves every number
          // finite.
          const auto held = l.collar_coords_at<T>(l.opt_psi_stem_, l.ci_, collar,
                                                  draw, pars, false);
          ok(odelia::util::to_passive(held.sigma.value) ==
                 odelia::util::to_passive(co.sigma.value),
             "sigma is the same number with the channel omitted" + at);
          ok(odelia::util::to_passive(held.ci.value) ==
                 odelia::util::to_passive(co.ci.value),
             "ci is the same number with the channel omitted" + at);
          ok(std::isfinite(odelia::ode::derivative_along(held.sigma.value)) &&
                 std::isfinite(odelia::ode::derivative_along(held.ci.value)),
             "both coordinates keep a finite row without the channel" + at);
        }
      }
    }
  }
  std::string names;
  for (const auto &[name, n] : reached) {
    names += (names.empty() ? "" : ", ") + name + "(" + std::to_string(n) + ")";
  }
  printf("    %d draws, %d output sets, %d coordinate pairs | kinds: %s\n", draws,
         uptakes, coords, names.c_str());
  // ⚠️ THE CLAIM IS ONLY EVIDENCE ABOUT THE KINDS IT REACHED, and shade death is
  // the one that carried the sentinel: a fixture that stopped reaching it would
  // pass this while checking nothing.
  ok(reached.count("shade-death") > 0, "the sweep reached a shade death");
  ok(reached.count("hydraulic-shutdown") > 0, "the sweep reached a shutdown");
  ok(reached.count("interior") > 0, "the sweep reached an interior point");
}

// The dry end has TWO bounds, closed by different conditions, and which one bound
// is a classification a tally has to be able to report.
//
// ⚠️ AT SHIPPED DEFAULTS ONE ARM NEVER WINS, so it needs a fixture built for it
// rather than a sweep that waits for one. root_psi_crit is 5.87 MPa against a
// continuity root near 1.3, so the root's own limit is never the min until the
// trait that sets it is lowered deliberately. Without this the kind exists and
// nothing reaches it, which is a classification no run can produce -- and while
// the distinction was a private bool, no tally could report it either.
//
// The two arms place the collar by DIFFERENT expressions -- a residual on one and
// a closed form on the other -- so the check is the placement as well as the
// name: at the root-limit arm the collar is root_psi_crit exactly.
void test_the_dry_end_reports_which_bound_closed_it() {
  printf("the dry end reports which of its two bounds closed it\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  const std::vector<double> psi_soil{2.0}, depth{1.0};
  const std::vector<double> mrp(psi_soil.size(),
                               1.0 / double(psi_soil.size()) / d.area_leaf);
  int continuity = 0, root_limit = 0;
  for (double root_P50 : {3.4, 2.0, 1.2, 0.9, 0.7, 0.5, 0.3}) {
    phylloptim::Leaf l = make_leaf(d, psi_soil, depth);
    l.set_traits(96.0, 2.680147, 3.4, 2.680147, root_P50, 1.5, 157.44, 0.30, 0.7,
                 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    const Kind kind = l.operating_point_kind();
    const std::string at = " at root_P50=" + std::to_string(root_P50);
    if (kind == Kind::BoundaryCrit) {
      ++continuity;
    } else if (kind == Kind::BoundaryRootCrit) {
      ++root_limit;
      // The collar IS the root's critical potential on this arm, which is what
      // distinguishes it from the continuity root beside it.
      near(l.opt_root_psi_, l.roots_.root_psi_crit, 1e-9,
           "the collar is the root's critical potential" + at);
    }
    printf("    root_P50=%4.2f  root_psi_crit=%7.4f  collar=%7.4f  %s\n", root_P50,
           l.roots_.root_psi_crit, l.opt_root_psi_,
           l.operating_point_kind_name(kind));
  }
  // The arm this exists for. The continuity arm is covered where it occurs
  // naturally -- the golden grid pins 18 of them at 25 C and none at 40 C, and
  // now pins this arm at zero over the same grid, which is the pair of counts
  // that says the two are told apart rather than merged.
  ok(root_limit > 0, "the sweep reached the root-limit arm");
  (void)continuity;
}

void test_rd_temperature_response() {
  printf("R_d rises with temperature\n");
  Drivers d;

  // R_d IS R_d_25 at the reference, exactly. The trait is the value there, not a
  // scale on anything.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    near(l.R_d_, l.R_d_25, 1e-12, "R_d at 25 C is R_d_25");
  }

  // The direction: R_d must RISE, including above Vcmax's thermal optimum near
  // 31 C, so 45 C is comfortably past it.
  double rd_prev = -1.0;
  for (double T : {25.0, 35.0, 45.0}) {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.update_temperature_dependent_params(T);
    ok(l.R_d_ > rd_prev, "R_d rises with temperature");
    rd_prev = l.R_d_;
  }

  // Tjoelker semantics. Q10 is evaluated at the MEAN of T and the reference, so a
  // 10 K rise from 25 C multiplies R_d by Q10(30) = 3.09 - 0.043*30 = 1.80 -- not
  // by the 2.015 that Q10(25) would give. Checking the number rather than just the
  // direction is what distinguishes this from a constant Q10.
  {
    phylloptim::Leaf a = make_leaf(d, {2.0}, {1.0});
    phylloptim::Leaf b = make_leaf(d, {2.0}, {1.0});
    a.update_temperature_dependent_params(25.0);
    b.update_temperature_dependent_params(35.0);
    const double q10_at_30 = 3.09 - 0.0430 * 30.0;
    near(b.R_d_, q10_at_30 * a.R_d_, 1e-12,
         "a 10 K rise scales R_d by Q10 at the midpoint temperature");
    ok(q10_at_30 < 2.0, "and that Q10 is below 2, i.e. the decline is active");
  }

  // A constant Q10 must stay reachable: slope zero, intercept the value.
  {
    phylloptim::Leaf a = make_leaf(d, {2.0}, {1.0});
    phylloptim::Leaf b = make_leaf(d, {2.0}, {1.0});
    a.rd_q10_intercept_ = 2.0;  a.rd_q10_slope_ = 0.0;
    b.rd_q10_intercept_ = 2.0;  b.rd_q10_slope_ = 0.0;
    a.update_temperature_dependent_params(25.0);
    b.update_temperature_dependent_params(35.0);
    near(b.R_d_, 2.0 * a.R_d_, 1e-12,
         "slope zero recovers a constant Q10 exactly");
  }

  // A measured value is used verbatim, and is INDEPENDENT of vcmax_25 -- R_d_25 is
  // a trait in its own right, not a fraction of Vcmax.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.R_d_25 = 0.525;                 // Sabot's Rlref, as measured
    l.update_temperature_dependent_params(25.0);
    near(l.R_d_, 0.525, 1e-12, "a set R_d_25 is used verbatim at 25 C");
    l.vcmax_25 = 2.0 * l.vcmax_25;
    l.update_temperature_dependent_params(25.0);
    near(l.R_d_, 0.525, 1e-12, "and does not follow vcmax_25");
  }

  // An unset or negative R_d_25 FAILS rather than falling back to a derivation.
  for (double bad : {std::numeric_limits<double>::quiet_NaN(), -1.0}) {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.R_d_25 = bad;
    bool threw = false;
    try {
      l.update_temperature_dependent_params(30.0);
    } catch (const std::runtime_error &) {
      threw = true;
    }
    ok(threw, "an unusable R_d_25 is refused, not worked around");
  }

  // The response, and the LIMIT it makes reachable, printed rather than asserted
  // cell by cell: a higher R_d can drive net assimilation negative across the whole
  // [gamma*, ca] bracket, and then there is no supply == demand root at all. At the
  // defaults that happens by 45 C, which is a number a caller needs.
  printf("  R_d and assimilation against leaf temperature:\n");
  for (double T : {15.0, 25.0, 35.0, 40.0, 45.0}) {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.update_temperature_dependent_params(T);
    l.find_root_collar_psi();
    printf("    T = %4.1f C   R_d %6.3f   A %8.4f%s\n", T, l.R_d_,
           l.assim_colimited_,
           l.ci_at_compensation_point() ? "   <-- shut down" : "");
  }

  // ⚠️ AND IT MUST SHUT DOWN RATHER THAN THROW. A leaf too hot to gain carbon at
  // any internal CO2 is a physical state, not a solver failure -- the energy-balance
  // path always treated it that way and the prescribed-temperature path used to
  // throw, which was only invisible while R_d was too small to get there.
  {
    phylloptim::Leaf hot = make_leaf(d, {2.0}, {1.0});
    bool threw = false;
    try {
      hot.update_temperature_dependent_params(45.0);
      hot.find_root_collar_psi();
    } catch (const std::exception &) {
      threw = true;
    }
    ok(!threw, "a leaf too hot to gain carbon shuts down instead of throwing");
    ok(hot.ci_at_compensation_point(),
       "and says so, rather than reporting an ordinary operating point");
    near(hot.ci_, hot.gamma_ * hot.umol_per_mol_to_Pa_, 1e-12,
         "ci sits at the compensation point");
    ok(hot.assim_colimited_ <= 0.0, "with no net carbon gain");
  }

  // ⚠️ A GENUINE solver failure must STILL throw -- the point of not simply
  // deleting the gate. An unsatisfiable bracket where the root IS inside it (here
  // forced by a non-finite target via a poisoned ca) is a bug, not a hot leaf.
  {
    phylloptim::Leaf ok_leaf = make_leaf(d, {2.0}, {1.0});
    ok_leaf.find_root_collar_psi();
    ok(!ok_leaf.ci_at_compensation_point(),
       "an ordinary leaf is not flagged as shut down");
  }
}

// set_traits exists so a gradient loop can perturb a trait without rebuilding the
// object. That is only worth having if re-traiting is INDISTINGUISHABLE from
// constructing afresh, so that is what is asserted -- bit-exactly, which is a
// statement neither a tolerance nor an eyeball could make.
//
// Bit-exactness is what makes this check work: the two ways to reach the same
// traits share no code, so any piece of derived state that set_traits fails to
// refresh shows up as a difference. Three pieces were candidates, and each is a
// real trap rather than a hypothetical one:
//
//   * the two pre-integrated vulnerability splines (stem_b/stem_c, root_b/root_c);
//   * vcmax_/jmax_/R_d_, behind set_physiology's (leaf_temp, atm_o2_kpa) cache --
//     the nastiest of the three, because the natural repair of "set the drivers
//     again" is exactly what takes the cache hit;
//   * the solved operating point itself (hazard 8).
//
// A `l.vcmax_25 = x` written by hand passes none of this, which is why the traits
// are not bound as settable fields.
// The homogeneity identity that makes a gradient in stem_b free.
//
// G(psi; s*b, c) = s * G(psi/s; b, c), because G(psi; b, c) = b * g(psi/b; c) --
// stem_b enters only as a scale on both axes, and the knot grid scales with it,
// so the identity holds for the SPLINE and not merely for the integral it
// approximates. perturb_stem_P50() exploits that instead of reseeding 101
// incomplete gammas, which is 21.8 us against 0.001. It moves the TRAIT, and
// stem_b is a fixed multiple of it at fixed stem_c, so the same identity applies.
//
// ⚠️ This is checked against a REBUILD and not against a formula, because what
// has to be true is that the rescaled curve is the one a rebuild would have
// produced. It is not asserted bit-identical: the two differ by the rounding of
// one multiply and one divide per evaluation, which the nested solvers amplify to
// the ~1e-9 floor this project's guide records.
void test_perturb_stem_P50_matches_a_rebuild() {
  printf("perturb_stem_P50 reproduces a rebuilt stem curve\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};

  // The range runs well below 1 on purpose. It used to stop at 0.98 because
  // psi_crit was an independent trait: the curve's domain is [0, b*log(100)^(1/c)]
  // and shrinking stem_b pushed a fixed psi_crit off the end of it. psi_crit is
  // now a quantile of the same curve, so it scales with stem_b and the domain
  // cannot bind however far the trait moves. A gradient perturbs by ~1e-6; the
  // wide factors are here to hold that property.
  for (double factor : {0.5, 0.98, 0.999, 1.000001, 1.05, 1.3, 2.0}) {
    const double P50_new = 3.4 * factor;

    // The reference: stem_b through set_traits, which rebuilds the spline.
    phylloptim::Leaf rebuilt = make_leaf(d, {2.0}, {1.0});
    rebuilt.find_root_collar_psi();
    rebuilt.set_traits(96.0, 2.680147, P50_new, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    rebuilt.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                           d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa,
                           d.atm_kpa);
    rebuilt.find_root_collar_psi();

    // The rescale: no rebuild, and no set_physiology either -- nothing it
    // derives depends on stem_b, which is half of why this is cheap.
    phylloptim::Leaf rescaled = make_leaf(d, {2.0}, {1.0});
    rescaled.find_root_collar_psi();
    rescaled.perturb_stem_P50(P50_new);
    rescaled.find_root_collar_psi();

    const std::string what = " at stem_P50 x " + std::to_string(factor);
    const double tol = 1e-8;
    near(rescaled.opt_root_psi_, rebuilt.opt_root_psi_, tol,
         "the collar matches a rebuild" + what);
    near(rescaled.opt_psi_stem_, rebuilt.opt_psi_stem_, tol,
         "psi_stem matches a rebuild" + what);
    near(rescaled.assim_colimited_, rebuilt.assim_colimited_, tol,
         "assimilation matches a rebuild" + what);
    near(rescaled.transpiration_, rebuilt.transpiration_, tol,
         "transpiration matches a rebuild" + what);
    near(rescaled.profit_, rebuilt.profit_, tol,
         "profit matches a rebuild" + what);

    // ⚠️ And set_traits is the way BACK: while the spline is built at another
    // stem_b, "the parameters did not move" must not be read as "there is
    // nothing to do". Restoring the defaults has to give a bit-identical leaf,
    // which is only true if the rebuild is forced.
    phylloptim::Leaf fresh = make_leaf(d, {2.0}, {1.0});
    fresh.find_root_collar_psi();
    rescaled.set_traits(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    rescaled.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                            d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa,
                            d.atm_kpa);
    rescaled.find_root_collar_psi();
    ok(rescaled.opt_root_psi_ == fresh.opt_root_psi_,
       "set_traits restores a bit-identical collar" + what);
    ok(rescaled.assim_colimited_ == fresh.assim_colimited_,
       "set_traits restores bit-identical assimilation" + what);
  }
}

// The other half of the shortcut: what it costs to STOP using it (#74).
//
// `set_traits()` is the way back, and it gets there by rebuilding -- which is
// correct and, on the gradient's restore path, unnecessary. So a gradient
// differentiating stem_b paid for a rebuild once per observation and the 24.5x
// per-perturbation identity was worth 2.4x through the batch. `apply()` now undoes
// the displacement with the shortcut first.
//
// ⚠️ COUNTED, NOT TIMED, and that is the whole design of this test. A rebuild at
// an unchanged (stem_b, stem_c) is bit-identical to not rebuilding -- it is a pure
// 11.9 us, invisible to the golden file, to this suite's value assertions and to
// anything the R layer can see. `stem_curve_builds_` exists so the claim is an
// integer instead of a stopwatch reading on a shared runner.
void test_stem_curve_shortcut_needs_no_rebuild() {
  printf("undoing the stem-curve shortcut costs no rebuild\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
  const double P50_0 = 3.4;

  // --- the identity itself, without the gradient in the way ------------------
  //
  // Bit-identity here is what makes the change free, and it is a property of
  // everything `perturb_stem_P50()` writes being a pure function of the trait,
  // computed by the same expression `set_traits()` uses: these splines ARE the
  // ones built at the defaults, so putting the trait back gives the object a
  // rebuild would have produced, not an approximation to it.
  {
    phylloptim::Leaf undone = make_leaf(d, {2.0}, {1.0});
    undone.find_root_collar_psi();
    const long builds = undone.stem_curve_builds_;

    undone.perturb_stem_P50(P50_0 * 1.001);
    undone.perturb_stem_P50(P50_0);  // the way back, WITH the shortcut
    undone.set_traits(96.0, 2.680147, P50_0, 2.680147, P50_0, 1.5,
                      157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    ok(undone.stem_curve_builds_ == builds,
       "set_traits after an undone displacement rebuilds nothing");
    undone.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil,
                          depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                          d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    undone.find_root_collar_psi();

    phylloptim::Leaf fresh = make_leaf(d, {2.0}, {1.0});
    fresh.find_root_collar_psi();
    ok(undone.opt_root_psi_ == fresh.opt_root_psi_,
       "the collar is bit-identical to a rebuilt leaf");
    ok(undone.opt_psi_stem_ == fresh.opt_psi_stem_,
       "psi_stem is bit-identical to a rebuilt leaf");
    ok(undone.assim_colimited_ == fresh.assim_colimited_,
       "assimilation is bit-identical to a rebuilt leaf");
    ok(undone.transpiration_ == fresh.transpiration_,
       "transpiration is bit-identical to a rebuilt leaf");
    ok(undone.profit_ == fresh.profit_,
       "profit is bit-identical to a rebuilt leaf");
  }

  // --- and through the batch, which is where the cost was being paid ---------
  const double kmax = d.K_s * d.theta / d.h;
  // ⚠️ POSITIONAL, AND IT FAILS SILENTLY. An aggregate initialiser shorter than
  // `n_theta` is legal C++ and zero-fills the rest, so adding a trait shifts `kmax`
  // and `resistance` down a slot and drops `resistance` off the end WITHOUT a
  // compiler diagnostic. Adding JS22_gamma put `kmax` into JS22_gamma's slot and
  // left kmax itself 0.0, and the only symptom was this test's three observations
  // failing to solve. `set_traits` above catches the same mistake at compile time;
  // this does not, so count the entries against `n_theta` when you touch it.
  double theta[phylloptim::gradient::n_theta] = {
      96.0,   2.680147, P50_0, 2.680147, P50_0, 1.5,
      157.44, 0.30,     0.7,   0.99,     7.5,   kRd25, kGammaJS22,
      kCMaxA, kCMaxB,  kmax,      0.0};
  // The two prices are left zero-filled: this test differentiates the collar
  // route, which reads neither. ⚠️ THE ASSERTION IS WHAT MAKES THE SHORT
  // INITIALISER SAFE -- without it the next appended parameter shifts `kmax`
  // silently, which is the failure the comment above describes.
  static_assert(phylloptim::gradient::n_theta == 19,
                "theta above is positional and deliberately short; recount it "
                "against n_theta and update this assertion together");

  phylloptim::gradient::Drivers gd;
  gd.root_network = fixture::root_network(mrp, depth);
  gd.PPFD = d.PPFD;
  gd.psi_soil = psi_soil;
  gd.soil_depth = depth;
  gd.atm_vpd = d.atm_vpd;
  gd.ca = d.ca;
  gd.leaf_temp = d.leaf_temp;
  gd.atm_o2_kpa = d.atm_o2_kpa;
  gd.atm_kpa = d.atm_kpa;
  const std::vector<phylloptim::gradient::Drivers> obs{gd, gd, gd};

  // Three parameter sets, and the middle one is the one that would have been
  // missed by a fix written only for the end-of-loop restore: `stem_P50` FOLLOWED by
  // a parameter that owns no vulnerability curve leaves the displacement to be
  // undone by that parameter's own first perturbation instead.
  struct Case {
    const char *name;
    std::vector<int> pars;
  };
  const Case cases[] = {
      {"stem_P50", {phylloptim::gradient::par_stem_P50}},
      {"stem_P50 then vcmax_25",
       {phylloptim::gradient::par_stem_P50, phylloptim::gradient::par_vcmax_25}},
      {"vcmax_25 then stem_P50",
       {phylloptim::gradient::par_vcmax_25, phylloptim::gradient::par_stem_P50}}};

  for (const Case &c : cases) {
    phylloptim::gradient::Settings s;  // fast_stem_curve defaults true
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.find_root_collar_psi();
    const long builds = l.stem_curve_builds_;
    const std::vector<phylloptim::gradient::Result> g =
        phylloptim::gradient::batch(l, theta, 1, obs, false, c.pars.data(),
                                    c.pars.size(), s);
    const std::string what = std::string(" for pars = ") + c.name;
    ok(g.size() == 3 &&
           g[0].status == phylloptim::gradient::Status::Interior,
       "the batch solved three interior observations" + what);
    ok(l.stem_curve_builds_ == builds,
       "three observations rebuild the stem curve zero times" + what);
    // Hazard 8, at the one place a shortcut could leave it violated: the batch
    // must not hand back a leaf quietly running on a rescaled curve.
    ok(l.stem_b == l.stem_b_spline_,
       "the batch leaves the leaf undisplaced" + what);
  }

  // The control, which is what stops the three counts above from being zero
  // because nothing ran: with the shortcut off, the same work rebuilds. Two
  // rebuilds per side per observation, plus the restore's own.
  {
    phylloptim::gradient::Settings s;
    s.fast_stem_curve = false;
    const int pars[1] = {phylloptim::gradient::par_stem_P50};
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.find_root_collar_psi();
    const long builds = l.stem_curve_builds_;
    phylloptim::gradient::batch(l, theta, 1, obs, false, pars, 1, s);
    ok(l.stem_curve_builds_ > builds,
       "with fast_stem_curve off, the same batch does rebuild");
    printf("  stem curve builds over 3 observations of d/dstem_P50:"
           " fast 0, slow %ld\n",
           l.stem_curve_builds_ - builds);
  }
}

void test_set_traits_matches_a_fresh_leaf() {
  printf("set_traits is indistinguishable from constructing afresh\n");
  Drivers d;

  // One perturbed trait from each group that has its own derived state: a
  // photosynthetic trait (the temperature cache), the stem curve, the root curve,
  // and a cost parameter that is read directly and so should need no rebuild.
  struct Case { const char *name; int which; double value; };
  const Case cases[] = {{"vcmax_25", 0, 96.0 * 1.05},
                        {"stem_P50", 2, 3.4 * 1.05},
                        {"root_P50", 4, 3.4 * 1.05},
                        {"TF24_cost_scale", 10, 7.5 * 1.05},
                        {"stem_c", 1, 2.680147 * 1.05},
                        {"root_c", 3, 2.680147 * 1.05}};

  for (const Case &c : cases) {
    // The default trait vector, in set_traits' own argument order (R_d_25 apart,
    // which the constructor does not take).
    double t[11] = {96.0, 2.680147, 3.4,  2.680147, 3.4, 1.5,
                    157.44, 0.30,   0.7,  0.99,     7.5};
    t[c.which] = c.value;

    // Fresh: the traits go through the constructor.
    phylloptim::Leaf fresh(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
                     t[8], t[9], 1e-3, 100, 1e-3, 1000, t[10]);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    fresh.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                         d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    fresh.find_root_collar_psi();

    // Reused: solved once at the DEFAULTS first, so every cache is warm and
    // pointing at the old traits before set_traits runs. Solving first is what
    // makes this a test rather than a coincidence -- on a cold object the
    // temperature cache would miss anyway and the trap would not fire.
    phylloptim::Leaf reused = make_leaf(d, {2.0}, {1.0});
    reused.find_root_collar_psi();
    reused.set_traits(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
                      t[8], t[9], t[10], kRd25, kGammaJS22, kCMaxA, kCMaxB);
    reused.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                          d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    reused.find_root_collar_psi();

    const std::string what = std::string(" after set_traits(") + c.name + ")";
    ok(reused.vcmax_ == fresh.vcmax_, "vcmax_ is bit-identical" + what);
    ok(reused.jmax_ == fresh.jmax_, "jmax_ is bit-identical" + what);
    ok(reused.R_d_ == fresh.R_d_, "R_d_ is bit-identical" + what);
    ok(reused.opt_root_psi_ == fresh.opt_root_psi_,
       "the collar is bit-identical" + what);
    ok(reused.opt_psi_stem_ == fresh.opt_psi_stem_,
       "psi_stem is bit-identical" + what);
    ok(reused.assim_colimited_ == fresh.assim_colimited_,
       "assimilation is bit-identical" + what);
    ok(reused.profit_ == fresh.profit_, "profit is bit-identical" + what);
    ok(reused.transpiration_ == fresh.transpiration_,
       "transpiration is bit-identical" + what);
  }

  // "Change the trait, then set the drivers again" really does recompute the
  // temperature block, checked end to end through the route a caller takes.
  //
  // ⚠️ THIS NO LONGER ISOLATES set_traits' CACHE INVALIDATION, and the comment
  // here claimed it did for two releases. It said vcmax_ is derived behind a cache
  // "keyed on (leaf_temp, atm_o2_kpa) and NOT on the traits", so the recipe worked
  // "only because set_traits invalidates the cache". #55 widened the key to every
  // scalar update_temperature_dependent_params() reads -- `vcmax_25` included --
  // so the trap is now closed TWICE OVER: the key would miss on the changed trait
  // even if set_traits cleared nothing. The assertion below cannot tell which
  // mechanism carried it, and pretending otherwise is how a test comes to describe
  // a guarantee it is not testing.
  //
  // The two halves are covered separately and deliberately:
  //   * the KEY -- test_temperature_params_invalidate_cache, which sets a response
  //     parameter and re-drives without going through set_traits at all;
  //   * set_traits' CLEARING -- the bit-exact comparisons above (the splines and
  //     the solved operating point are not in any key), plus the NA assertions in
  //     test_prescribed_lambda_survives_redriving.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double vcmax_before = l.vcmax_;
    l.set_traits(96.0 * 2.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    // The same leaf_temp and atm_o2_kpa -- which used to be the whole key, and so
    // used to be what ARMED the trap. It no longer is: `vcmax_25` has moved, so the
    // key differs whatever set_traits did.
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    near(l.vcmax_ / vcmax_before, 2.0, 1e-12,
         "doubling vcmax_25 doubles vcmax_ at an unchanged temperature");
  }

  // The splines really are rebuilt, and only when their own pair moves. Read
  // through proportion_of_conductivity (the closed form) against transpiration()
  // (the spline): the two describe the same curve, so they move together or the
  // object is inconsistent.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double E_before = l.transpiration(3.0, 1.0);
    l.set_traits(96.0, 2.680147, 3.4 * 1.5, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    // A less vulnerable stem (larger stem_b) holds more conductivity, so the
    // integral of the vulnerability curve over the same span is larger.
    ok(l.transpiration(3.0, 1.0) > E_before,
       "raising stem_b rebuilds the transpiration spline");
  }

  // The #25 boundary is enforced here too. A bare field write would bypass it,
  // which is the fourth reason these are not settable fields.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double good[11] = {96.0, 2.680147, 3.4,  2.680147, 3.4, 1.5,
                             157.44, 0.30,   0.7,  0.99,     7.5};
    // Only two potentials are settable now -- the derived scales and critical
    // potentials take their sign from these, so there is nothing else to reject.
    const int signed_positions[] = {2, 4};
    const char *labels[] = {"stem_P50", "root_P50"};
    for (int k = 0; k < 2; ++k) {
      double t[11];
      for (int i = 0; i < 11; ++i) {
        t[i] = good[i];
      }
      t[signed_positions[k]] = -t[signed_positions[k]];  // the pre-#25 sign
      bool threw = false;
      try {
        l.set_traits(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
                     t[8], t[9], t[10], kRd25, kGammaJS22, kCMaxA, kCMaxB);
      } catch (const std::runtime_error &) {
        threw = true;
      }
      ok(threw, std::string("set_traits rejects a negative ") + labels[k]);
    }
  }
}

// A prescribed `CF77_lambda_` survives everything that is not a caller writing to it.
// It is the Cowan-Farquhar marginal value of water, an INPUT, so a re-driving call
// that cleared it would make whether a sweep kept its price depend on which of two
// interchangeable-looking calls came next.
//
// ⚠️ The two arms are the point, not the pair of assertions. `set_drivers` always
// kept the value and `set_traits` always lost it; asserting only one arm passes
// on the code this test exists to reject.
void test_prescribed_lambda_survives_redriving() {
  printf("a prescribed CF77_lambda_ survives set_traits and set_drivers\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};

  // A fresh Leaf reads the NA sentinel: "no longer reset" must not have become
  // "never initialised", which is the one way this change could go wrong quietly.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    ok(std::isnan(l.CF77_lambda_), "CF77_lambda_ is NA on a freshly constructed leaf");
  }

  const double prescribed = 30.0;

  // Arm 1: re-drive. This arm passed before the fix too.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    ok(l.CF77_lambda_ == prescribed, "set_physiology leaves a prescribed CF77_lambda_ alone");
  }

  // Arm 2: re-trait. This is the arm that lost the value.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.set_traits(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    ok(l.CF77_lambda_ == prescribed, "set_traits leaves a prescribed CF77_lambda_ alone");
  }

  // And it survives a SOLVE, which is the case a sweep actually runs: solve at
  // the defaults, re-trait, solve again, all on one object.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.find_root_collar_psi();
    l.set_traits(96.0 * 1.05, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    ok(l.CF77_lambda_ == prescribed,
       "CF77_lambda_ survives a solve, a re-trait and a second solve");
  }

  // The derived state around it is still wiped -- the fix removed two lines from
  // setup_clean_leaf(), and taking any more would reopen hazard 8.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.find_root_collar_psi();
    ok(!std::isnan(l.opt_psi_stem_), "the solve seated an operating point");
    l.set_traits(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB);
    ok(std::isnan(l.opt_psi_stem_),
       "set_traits still clears the solved operating point");
    ok(std::isnan(l.profit_), "set_traits still clears profit_");
    ok(std::isnan(l.R_d_), "set_traits still clears the temperature block");
  }
}

// Capture the message a throwing call produced, or "" if it did not throw. Used
// by every test below that asserts on what an error SAYS rather than that one
// happened -- see test_out_of_domain_names_the_spline for why the wording is
// load-bearing here.
std::string message_of(const std::function<void()>& f) {
  try {
    f();
  } catch (const std::exception& e) {
    return std::string(e.what());
  }
  return std::string();
}

bool mentions(const std::string& haystack, const char* needle) {
  return haystack.find(needle) != std::string::npos;
}

// The realised knot grid against the intended one (#92). Nothing checked this
// before, which is exactly why the accumulating `psi += step` loop could drop its
// last knot for years: the shortfall is invisible from outside unless the domain
// edge is compared with the function that is supposed to define it.
// The kg <-> mol water conversions are reciprocal, and the round trip is the
// identity (#51). Asserted rather than documented: the demand side converts
// transpiration kg -> mol and the supply side converts uptake mol -> kg, so a
// second literal anywhere breaks the round trip silently.
void test_water_mass_conversions_are_reciprocal() {
  printf("kg <-> mol water conversions\n");
  // Exact to the last bit in one direction: kg_to_mol_h2o IS 1/molar_mass_h2o.
  ok(phylloptim::kg_to_mol_h2o == 1.0 / phylloptim::molar_mass_h2o,
     "kg_to_mol_h2o is exactly the reciprocal of the molar mass");
  ok(phylloptim::kg_per_mol_h2o == phylloptim::molar_mass_h2o,
     "kg_per_mol_h2o IS the molar mass, not a second copy of it");

  // The round trip cannot be bit-exact -- 1/x then *x is not the identity for most
  // x -- so the claim is that it is within one rounding, which is what "reciprocal"
  // can actually buy.
  for (double kg : {1e-9, 1e-6, 1e-3, 1.0, 1e3}) {
    const double round_trip = kg * phylloptim::kg_to_mol_h2o *
                              phylloptim::kg_per_mol_h2o;
    near(round_trip, kg, 4.0 * std::numeric_limits<double>::epsilon(),
         "kg -> mol -> kg returns the input at " +
             phylloptim::util::format_double(kg));
  }

  // The value is water's molar mass, not a neighbouring one: 18.015 g/mol from the
  // standard atomic weights, 2(1.008) + 15.999.
  near(phylloptim::molar_mass_h2o * 1000.0, 18.015, 1e-9,
       "the molar mass is water's, in g/mol");
}

// Both Arrhenius curves are exactly inert at 25 C, which is why a change to any
// temperature response is invisible on a 25 C-only grid. That is the claim worth
// holding: it bounds the blast radius of anything touched in this area.
void test_gas_constant_and_arrhenius_reference_point() {
  printf("gas constant and the 25 C reference\n");
  // Exact SI quantity, so there is no tolerance to allow.
  near(phylloptim::gas_constant, 8.314462618153240, 1e-15,
       "gas_constant is the SI value");

  // Inertness at the reference temperature, through the object: every reference
  // value in this model is DEFINED at 25 C, so a response can only be seen away
  // from it.
  Drivers d;
  phylloptim::Leaf ref = make_leaf(d, {2.0}, {1.0});   // d.leaf_temp is 25
  ok(ref.vcmax_ == ref.vcmax_25, "vcmax_ IS vcmax_25 at 25 C, bit for bit");
  ok(ref.jmax_ == ref.jmax_25, "jmax_ IS jmax_25 at 25 C, bit for bit");
  ok(ref.R_d_ == ref.R_d_25, "R_d_ IS R_d_25 at 25 C, bit for bit");
  // gamma_, kc_ and ko_ also sit at their 25 C reference values, but they carry a
  // umol/mol -> Pa conversion through atm_kpa_ that this test would have to restate
  // to check -- and a second copy of a conversion is what hazard 1 is about. The
  // three identities above are the same claim without the duplication.

  // And NOT inert away from it, which is the other half: a test that only checked
  // 25 C would pass on a broken response curve.
  Drivers hot;
  hot.leaf_temp = 40.0;
  phylloptim::Leaf warm = make_leaf(hot, {2.0}, {1.0});
  ok(warm.vcmax_ != warm.vcmax_25, "vcmax_ has moved at 40 C");
  ok(warm.jmax_ != warm.jmax_25, "jmax_ has moved at 40 C");
  ok(warm.R_d_ > ref.R_d_, "respiration is higher at 40 C than at 25 C");
}

void test_knot_grid_reaches_its_intended_domain() {
  printf("vulnerability knot grid\n");
  const double resolutions[] = {10.0, 50.0, 100.0, 200.0, 1000.0};
  int wrong_count = 0, wrong_end = 0, not_increasing = 0, cases = 0;
  double worst_end_rel = 0.0;
  for (double res : resolutions) {
    for (int bi = 1; bi <= 60; ++bi) {
      for (double c : {0.4, 1.0, 2.04, 2.680147, 6.0, 12.0}) {
        const double b = bi / 10.0;
        std::vector<double> x, y;
        phylloptim::cumulative_vulnerability_integral(b, c, res, x, y);
        ++cases;
        // ⚠️ BIT-EXACT on the endpoint, not `near`. The whole defect was an
        // endpoint that was *close* to vulnerability_psi_max and not equal to it,
        // and both splines built from this grid disable extrapolation -- so a
        // last knot one ULP short is the difference between a lookup at
        // vulnerability_psi_max working and throwing. A tolerance here would pass
        // on the code this test exists to reject.
        const double want_end = phylloptim::vulnerability_psi_max(b, c);
        if (x.size() != static_cast<std::size_t>(res) + 1) ++wrong_count;
        if (x.back() != want_end) {
          ++wrong_end;
          worst_end_rel = std::max(worst_end_rel,
                                   std::abs(x.back() - want_end) / want_end);
        }
        if (x.size() != y.size()) ++not_increasing;
        for (std::size_t i = 1; i < x.size(); ++i) {
          if (!(x[i] > x[i - 1])) { ++not_increasing; break; }
        }
      }
    }
  }
  ok(cases == 1800, "the sweep ran the grid it meant to");
  ok(wrong_count == 0,
     "every grid has exactly resolution + 1 knots (" +
         std::to_string(wrong_count) + " of " + std::to_string(cases) + " did not)");
  ok(wrong_end == 0,
     "every grid ends exactly at vulnerability_psi_max (" +
         std::to_string(wrong_end) + " did not; worst rel " +
         std::to_string(worst_end_rel) + ")");
  ok(not_increasing == 0, "every grid is strictly increasing and paired with its y");

  // A resolution below one control point is refused by name rather than reaching
  // the interpolator, which used to report it as its own domain problem.
  const std::string msg = message_of(
      [] { std::vector<double> x, y;
           phylloptim::cumulative_vulnerability_integral(3.9, 2.7, 0.0, x, y); });
  ok(mentions(msg, "resolution"), "resolution < 1 is refused by name");
}

// psi_crit against the domain stem_P50/stem_c set. This USED to be a validation
// problem: psi_crit was a free trait, the stem curve's knot grid stops at P99, and
// every solve evaluates that curve AT psi_crit -- so an inconsistent pair failed
// from inside the interpolator in a message naming neither trait.
//
// It is now a structural guarantee. psi_crit is the 5%-conductivity quantile of
// the same curve, so containment reduces to log(20) < log(100) and holds for every
// c > 0. What is asserted here is that the code realises that algebra: the
// combination the old check refused now constructs, and the derived value is the
// quantile it claims to be.
void test_psi_crit_lies_inside_the_stem_curve_domain() {
  printf("psi_crit inside the curve's domain, by construction\n");

  // The containment, over a grid far wider than any plausible trait.
  int cases = 0, outside = 0;
  double worst_ratio = 0.0;
  for (double P50 : {0.4, 1.0, 3.4, 7.6291, 20.0}) {
    for (double c : {0.5, 1.0, 2.0, 2.680147, 3.5463, 8.0}) {
      ++cases;
      const double b = P50 / std::pow(std::log(2.0), 1.0 / c);
      const double psi_crit = b * std::pow(std::log(1.0 / 0.05), 1.0 / c);
      const double p99 = phylloptim::vulnerability_psi_max(b, c);
      if (!(psi_crit < p99)) ++outside;
      worst_ratio = std::max(worst_ratio, psi_crit / p99);
    }
  }
  ok(cases == 30, "the containment sweep ran the grid it meant to");
  ok(outside == 0,
     "the derived psi_crit is inside the derived domain for every (P50, c) (" +
         std::to_string(outside) + " of " + std::to_string(cases) + " were not)");
  // Ratio-of-quantiles, so it is the same number at every c: it is what the
  // margin IS, not a sampled worst case.
  printf("  worst psi_crit / P99 over the grid: %.6f\n", worst_ratio);

  // Sabot et al. (2022) P50/P88 territory -- the pair the old check refused,
  // because a psi_crit chosen as though it were free landed past this curve's
  // domain. There is nothing left to get wrong, so it constructs.
  ok(message_of([] {
       phylloptim::Leaf l(96, 3.5463, 7.6291, 2.680147, 3.4, 1.5, 157.44, 0.30,
                          0.7, 0.99, 1e-3, 100, 1e-3, 1000, 7.5);
     }).empty(),
     "a brittle stem curve constructs without a psi_crit to reconcile");

  // The derived values ARE the quantiles claimed, bit-exactly -- the same
  // expression, so a tolerance here would pass on a different derivation.
  phylloptim::Leaf def;
  ok(def.stem_b == def.stem_P50 / std::pow(std::log(2.0), 1.0 / def.stem_c),
     "stem_b is the scale that stem_P50 implies");
  ok(def.psi_crit ==
         def.stem_b * std::pow(std::log(1.0 / 0.05), 1.0 / def.stem_c),
     "psi_crit is P95 of the stem curve");
  ok(def.roots_.root_b ==
         def.roots_.root_P50 / std::pow(std::log(2.0), 1.0 / def.roots_.root_c),
     "root_b is the scale that root_P50 implies");
  ok(def.roots_.root_psi_crit ==
         def.roots_.root_b *
             std::pow(std::log(1.0 / 0.05), 1.0 / def.roots_.root_c),
     "root_psi_crit is P95 of the root curve");

  // perturb_stem_P50 is the third route in, and the one where the domain used to
  // be able to walk out from under a fixed psi_crit. Both now scale together, so
  // an arbitrarily large move is safe and psi_crit follows it.
  phylloptim::Leaf pl;
  const double psi_crit_0 = pl.psi_crit;
  ok(message_of([&] { pl.perturb_stem_P50(pl.stem_P50 * 0.25); }).empty(),
     "a large downward perturbation is no longer a domain error");
  ok(pl.psi_crit < psi_crit_0,
     "psi_crit moved with the curve rather than staying put");
  ok(pl.psi_crit ==
         pl.stem_b * std::pow(std::log(1.0 / 0.05), 1.0 / pl.stem_c),
     "and it is still P95 of the perturbed curve");
}

void test_bad_input_throws() {
  printf("input validation\n");
  Drivers d;
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  bool threw = false;
  try {
    std::vector<double> psi_soil{2.0}, depth{1.0, 2.0}, mrp{1.0 / d.area_leaf};
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(threw, "mismatched soil vector lengths throw std::runtime_error");
}

// Infeasibility is a distinguishable kind of failure (#57).
//
// ⚠️ THE C++ SIDE MATTERS ON ITS OWN, not just as the R layer's supplier: plant links
// these headers directly, so it sees the type rather than the token. Both halves are
// asserted here -- a `catch (const infeasible_error&)` works, and the message carries
// the prefix R parses.
//
// ⚠️ And the asymmetry is asserted too. An input-validation failure must NOT be an
// `infeasible_error`, because the whole point is that a caller may legitimately
// swallow one and must not be able to swallow the other. A test that only checked the
// positive direction would pass on a `stop_infeasible` applied to everything.
void test_infeasible_is_a_distinct_failure() {
  printf("infeasibility is distinguishable from a caller error\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};

  // A well-formed call whose operating point cannot be evaluated: with the stem
  // conductance this small the transpiration the solve needs is off the end of the
  // inverse spline's domain.
  {
    phylloptim::Leaf l;
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     1e-30, d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa,
                     d.atm_kpa);
    bool caught_as_infeasible = false;
    std::string msg;
    try {
      l.find_root_collar_psi();
    } catch (const phylloptim::util::infeasible_error &e) {
      caught_as_infeasible = true;
      msg = e.what();
    } catch (const std::runtime_error &e) {
      msg = e.what();
    }
    ok(caught_as_infeasible,
       "an unevaluable operating point throws util::infeasible_error");
    ok(msg.rfind(phylloptim::util::infeasible_token("stem_curve_domain"), 0) == 0,
       "and its message opens with the token R parses");
  }

  // The other direction: a mismatched driver vector is NOT infeasible. This is #39's
  // warning as an assertion -- an all-NA or wrong-length driver classified as
  // infeasible would be swallowed by the caller's tryCatch and the fit would report
  // a plausible likelihood over whichever rows survived.
  {
    phylloptim::Leaf l;
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    bool infeasible = false, threw = false;
    try {
      std::vector<double> bad_depth{1.0, 2.0};
      l.set_physiology(fixture::root_network(mrp, bad_depth), d.PPFD, psi_soil,
                       bad_depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                       d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    } catch (const phylloptim::util::infeasible_error &) {
      infeasible = true;
      threw = true;
    } catch (const std::runtime_error &) {
      threw = true;
    }
    ok(threw, "a mismatched driver vector still throws");
    ok(!infeasible, "and it is NOT classified as infeasible");
  }

  // The token is a prefix and nothing else parses as one: these messages embed
  // `e.what()` from a nested solver, so a quoted token must not re-classify.
  ok(phylloptim::util::infeasible_token("uptake") ==
         "[phylloptim:infeasible:uptake] ",
     "the token has the exact shape R's anchored pattern expects");
}

// What an out-of-domain transport lookup says. The stem curve is the only
// interpolator here with extrapolation disabled, so it is the only one a lookup
// can throw on -- and there are two of them, they are inverses, and they carry
// different units, so odelia's message (which names the point and the domain but
// not the spline) is ambiguous in the way that matters. Localising plant#576 came
// down to which of the four call sites was asking; these assertions are what make
// that a read rather than a bisect.
// (message_of / mentions are defined above test_knot_grid_reaches_its_intended_domain,
// which is the first test that needs them.)

void test_out_of_domain_names_the_spline() {
  printf("out-of-domain reporting\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});

  // Forward direction, past the far end of the vulnerability curve.
  const std::string fwd = message_of([&] { l.transpiration(50.0, 0.0); });
  ok(mentions(fwd, "transpiration_from_psi"), "forward lookup names its spline");
  ok(mentions(fwd, "psi = 50"), "forward lookup reports the point");
  ok(mentions(fwd, "beyond the upper end"), "forward lookup reports which end");
  ok(mentions(fwd, "Leaf::transpiration"), "forward lookup names the caller");

  // Inverse direction, below the lower end -- the plant#576 signature. A
  // sufficiently negative flux puts E/K_max below the domain, which is the
  // statement "the collar cannot supply this, so no stem potential carries it".
  const std::string inv =
      message_of([&] { l.transpiration_to_psi_stem(-1e3, 0.0); });
  ok(mentions(inv, "INVERTED"), "inverse lookup names which direction it is");
  ok(mentions(inv, "beyond the lower end"), "inverse lookup reports which end");
  ok(mentions(inv, "E/K_max"), "inverse lookup names its argument's units");
  ok(mentions(inv, "Leaf::transpiration_to_psi_stem"),
     "inverse lookup names the caller");

  // The two are distinguishable, which is the entire point.
  ok(fwd != inv && !fwd.empty() && !inv.empty(),
     "the two splines give different messages");

  // A non-finite point must NOT become a domain complaint: the guard uses the
  // same comparison odelia does rather than negating an in-range test, so NaN
  // falls through to the spline and comes back non-finite. plant documents a
  // profit_psi_stem_TF(NA, .) -> NA contract built on this.
  const std::string nan_msg =
      message_of([&] { l.transpiration(std::nan(""), 0.0); });
  ok(nan_msg.empty(), "a non-finite psi_stem does not throw a domain error");
  ok(!std::isfinite(l.transpiration(std::nan(""), 0.0)),
     "a non-finite psi_stem returns non-finite");

  // In-domain reads are untouched.
  ok(std::isfinite(l.transpiration(2.5, 0.0)),
     "an in-domain lookup still returns a finite value");
}

// The rescaled path reports the domain in the CALLER's units, not the spline's.
// Under perturb_stem_P50 the value handed to the spline is psi/s, so quoting the
// spline's own endpoints would send the reader after a discrepancy that is not
// there.
void test_out_of_domain_under_rescale() {
  Drivers d;
  phylloptim::Leaf wide = make_leaf(d, {2.0}, {1.0});
  const std::string before = message_of([&] { wide.transpiration(50.0, 0.0); });

  phylloptim::Leaf rescaled = make_leaf(d, {2.0}, {1.0});
  rescaled.perturb_stem_P50(rescaled.stem_P50 * 2.0);
  const std::string after =
      message_of([&] { rescaled.transpiration(50.0, 0.0); });

  ok(!before.empty() && !after.empty(), "both report a domain failure");
  ok(before != after, "the rescaled domain is reported, not the spline's own");
  ok(mentions(after, "rescaled by"),
     "the rescale is named so the two domains are not confused");
}

// ===========================================================================
// Leaf-to-air VPD (#7)
// ---------------------------------------------------------------------------
void test_leaf_to_air_vpd() {
  printf("leaf-to-air VPD: the deficit Fick's law divides by\n");
  Drivers d;
  d.atm_vpd = 1.5;

  // Gate off, the leaf IS at air temperature, so the two deficits must agree --
  // and EXACTLY, not nearly. That is what keeps the prescribed-temperature path,
  // and therefore the golden file, untouched by this change.
  phylloptim::Leaf off = make_pm_leaf(d, {2.0}, {1.0}, false);
  ok(off.vpd_leaf_ == off.atm_vpd_,
     "gate off: vpd_leaf_ is bit-identical to the driver");
  off.find_root_collar_psi();
  ok(off.vpd_leaf_ == off.atm_vpd_,
     "gate off: and still is after a solve");

  // Gate on, with the leaf running hot: the deficit it sees is larger than the
  // air's, so the same water flux implies a SMALLER conductance.
  phylloptim::Leaf on = make_pm_leaf(d, {2.0}, {1.0}, true);
  on.Rn_ = 400.0;
  on.ra_ = 12.0;
  on.find_root_collar_psi();
  const double Tleaf = on.leaf_temp_from_E(on.transpiration_);
  ok(Tleaf > on.Tair_, "the gate-on leaf is hotter than the air here");
  ok(on.vpd_leaf_ > on.atm_vpd_,
     "gate on: a hotter leaf sees a larger deficit than the air's");
  // And the value is the definition, not an approximation of it.
  near(on.vpd_leaf_,
       on.atm_vpd_ + on.saturation_vapour_pressure(Tleaf) -
           on.saturation_vapour_pressure(on.Tair_),
       1e-12, "vpd_leaf_ = atm_vpd + esat(Tleaf) - esat(Tair)");
  printf("    Tair %.1f C, Tleaf %.2f C: D_air %.3f kPa -> D_leaf %.3f kPa (x%.2f)\n",
         on.Tair_, Tleaf, on.atm_vpd_, on.vpd_leaf_, on.vpd_leaf_ / on.atm_vpd_);

  // The consequence, stated as an identity rather than a direction: gs is E
  // rescaled by the deficit, so getting the deficit wrong scales gs by the ratio.
  near(on.stom_cond_CO2_,
       on.atm_kpa_ * on.transpiration_ * phylloptim::kg_to_mol_h2o /
           on.vpd_leaf_ / phylloptim::H2O_CO2_stom_diff_ratio,
       1e-12, "gs is the transpiration divided by the deficit it sees");
}

// Build a leaf on the single-potential path, which is what the ProfitMax entry
// points require and what Sicangco's model is (one soil potential, no root
// resistance network).
phylloptim::Leaf make_single_leaf(const Drivers &d, double psi_soil,
                                  bool gate = false) {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = gate;
  l.set_supply_single(0.0);
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  l.set_physiology(rn, d.PPFD, {psi_soil}, {1.0}, d.K_s * d.theta / d.h,
                   d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// ===========================================================================
// Sperry (2017) ProfitMax
// ---------------------------------------------------------------------------
// THE LOAD-BEARING CLAIM is that the normalised objective and this package's
// older `A - lambda*cost` form are the same function up to a positive scale, so
// they share an argmax when lambda = |A|max/(k_soil - kcrit). Everything the
// Sicangco et al. (2026) replication does rests on it.
// The normalised objective IS a constant-lambda objective, up to a positive
// scale. Asserted as the pointwise identity
//
//   A_max * ProfitMax(psi) + A_max * TC  ==  A(psi) - lambda* (k(psi_s) - k(psi))
//
// rather than by comparing two argmaxes. That is both tighter -- the identity
// holds to rounding, where two searches agree only to their bracket width -- and
// stronger, since a positive affine map preserves the argmax, so the identity
// implies the argmax claim while the converse does not.
//
// The right-hand side is built from `assim_colimited_` and
// `proportion_of_conductivity`, not from the normalisers ProfitMax cached, so the
// two sides do not share their arithmetic.
void test_profitmax_is_a_constant_lambda_objective() {
  printf("ProfitMax is the constant-lambda objective, up to a positive scale\n");
  Drivers d;
  d.PPFD = 1500.0;
  double worst = 0.0;
  int points = 0;

  for (double psi_soil : {0.5, 2.0}) {
    for (double t : {25.0, 40.0}) {
      d.leaf_temp = t;
      phylloptim::Leaf l = make_single_leaf(d, psi_soil);
      l.prepare_profitmax();
      const double A_max = l.profitmax_A_max_, k_span = l.profitmax_k_span_;
      if (!(A_max > 0.0) || !(k_span > 0.0)) continue;
      const double lambda_star = A_max / k_span;
      const double k_soil =
          l.leaf_specific_conductance_max_ *
          l.proportion_of_conductivity(psi_soil);

      for (double frac : {0.2, 0.4, 0.6, 0.8}) {
        const double psi = psi_soil + frac * (l.psi_crit - psi_soil);
        const double lhs_norm = l.profit_psi_stem_ProfitMax(psi, psi_soil);
        // Read AFTER the call: these describe the point just evaluated.
        const double A = l.assim_colimited_, TC = l.thermal_cost_;
        const double k =
            l.leaf_specific_conductance_max_ * l.proportion_of_conductivity(psi);

        const double lhs = A_max * lhs_norm + A_max * TC;
        const double rhs = A - lambda_star * (k_soil - k);
        const double scale = std::max(std::abs(lhs), 1.0);
        worst = std::max(worst, std::abs(lhs - rhs) / scale);
        ++points;
        near(lhs / scale, rhs / scale, 1e-12,
             "the identity holds at frac=" + std::to_string(frac) +
                 " psi_soil=" + std::to_string(psi_soil));
      }
      ok(std::isfinite(lambda_star) && lambda_star > 0.0,
         "and the reported lambda is finite and positive");
    }
  }
  printf("    %d points | worst relative departure from the identity: %.3e\n",
         points, worst);
}

void test_profitmax_normalisation() {
  printf("ProfitMax: what the normalisation does and does not remove\n");
  Drivers d;
  d.PPFD = 1500.0;

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  const std::vector<double> curve = l.profitmax_curve(101);
  const std::size_t n = 101;
  ok(curve.size() == 5 * n, "profitmax_curve returns five columns");

  // HC runs from 0 at the soil potential to 1 at psi_crit, monotonically. That is
  // the definition and it is what makes the cost unable to vanish.
  near(curve[2 * n + 0], 0.0, 1e-12, "HC is zero at the soil potential");
  near(curve[2 * n + (n - 1)], 1.0, 1e-9, "HC is one at psi_crit");
  bool hc_monotone = true;
  for (std::size_t i = 1; i < n; ++i) {
    if (!(curve[2 * n + i] >= curve[2 * n + i - 1])) hc_monotone = false;
  }
  ok(hc_monotone, "HC increases monotonically along the supply stream");

  // ⚠️ AND HC DOES NOT DEPEND ON kmax AT ALL. Both the numerator and the
  // denominator carry one factor of leaf_specific_conductance_max_, so it cancels
  // exactly. That is worth pinning because Sicangco's ProfitMaxkmax(T) arm gives
  // kmax a temperature response and the paper describes it as changing the cost:
  // it does not. It changes the SUPPLY, and reaches the cost only through which
  // potentials the leaf can reach and how hot it gets there.
  Drivers d2 = d;
  d2.K_s = d.K_s * 3.0;
  phylloptim::Leaf l2 = make_single_leaf(d2, 0.5);
  const std::vector<double> curve2 = l2.profitmax_curve(101);
  double hc_worst = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    hc_worst = std::max(hc_worst, std::abs(curve2[2 * n + i] - curve[2 * n + i]));
  }
  ok(hc_worst < 1e-12, "HC is invariant to a 3x change in kmax");
  printf("    HC under kmax x3: worst difference %.3e\n", hc_worst);

  // The optimiser lands on the curve's own maximum.
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  double best = -1e300;
  std::size_t at = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (curve[4 * n + i] > best) { best = curve[4 * n + i]; at = i; }
  }
  ok(std::abs(l.opt_psi_stem_ - curve[at]) < 2.0 * (curve[1] - curve[0]),
     "the optimiser lands within a grid step of the curve's maximum");
  ok(l.profit_ >= best - 1e-9, "and at no lower profit than the grid's best");
}

void test_profitmax_thermal_cost() {
  printf("ProfitMax: the thermal cost, and that it is inert when off\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 48.0;

  phylloptim::Leaf off = make_single_leaf(d, 0.5);
  off.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  ok(off.thermal_cost_ == 0.0, "gate off: TC is exactly zero");

  phylloptim::Leaf on = make_single_leaf(d, 0.5);
  on.use_thermal_cost_ = true;
  on.T50_ = 50.4;
  on.Tcrit_ = 46.5;
  // The gate reaches jmax_ through the temperature block, so the drivers have to
  // be re-supplied for the cache key to notice. Doing it the way a caller would.
  on = make_single_leaf(d, 0.5);
  on.use_thermal_cost_ = true;
  on.T50_ = 50.4;
  on.Tcrit_ = 46.5;
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  on.set_physiology(rn, d.PPFD, {0.5}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd,
                    d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  on.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  ok(on.thermal_cost_ > 0.0 && on.thermal_cost_ < 1.0,
     "gate on at 48 C: TC is in (0,1)");
  ok(on.jmax_ < off.jmax_, "and Jmax is scaled down by (1 - TC)");
  near(on.jmax_, off.jmax_ * (1.0 - on.thermal_cost_at(d.leaf_temp)), 1e-12,
       "by exactly that factor");
  printf("    Tleaf %.1f C, Tcrit %.1f, T50 %.1f: TC %.4f, Jmax %.3f -> %.3f\n",
         d.leaf_temp, on.Tcrit_, on.T50_, on.thermal_cost_, off.jmax_, on.jmax_);

  // ⚠️ THE COST AT Tcrit IS A FIXED 11.9%, AND THAT CONTRADICTS THE PAPER'S PROSE.
  // Sicangco et al. write that "(1 - TC) equals one for temperatures below Tcrit,
  // [so] Eqns 10 and 11 yield the same result under such conditions". The equation
  // they cite does not do that: with r = 2/(T50 - Tcrit), the argument at Tcrit is
  // exactly -2 whatever the two thresholds are, so
  //
  //     TC(Tcrit) = 1/(1 + e^2) = 0.1192...
  //
  // independent of parameterisation. A leaf sitting AT its critical temperature
  // has already lost 11.9% of Jmax in this model, and a leaf 4 K below it still
  // pays 3%. Pinned as an identity because it is a property of the functional
  // form rather than of the values in Table 2.
  near(on.thermal_cost_at(on.Tcrit_), 1.0 / (1.0 + std::exp(2.0)), 1e-14,
       "TC at Tcrit is 1/(1+e^2), not zero");
  {
    phylloptim::Leaf wide = make_single_leaf(d, 0.5);
    wide.use_thermal_cost_ = true;
    wide.T50_ = 55.0;
    wide.Tcrit_ = 43.0;
    near(wide.thermal_cost_at(wide.Tcrit_), on.thermal_cost_at(on.Tcrit_), 1e-14,
         "and is the same 11.9% for a threshold pair three times as wide");
  }
  near(on.thermal_cost_at(on.T50_), 0.5, 1e-12, "and is exactly 0.5 at T50");
}

// ===========================================================================
// The closed form as a SELECTABLE METHOD -- the third axis of set_model()
// ---------------------------------------------------------------------------
// Three things this has to get right, and none of them is about accuracy (which
// tests/testthat/test-closed-form.R measures against a driver grid):
//
//  1. the exact arm is untouched -- bit-identical, not close;
//  2. every output the exact solve writes is written here too (hazard 8: `Leaf` is
//     a value member plant reuses per individual, so a field this path declined to
//     write would report the PREVIOUS solve's value);
//  3. everything the closed form cannot do is refused rather than approximated,
//     and the two refusals say different things -- one names a missing derivative,
//     the other names a fixed point that cannot be inverted.
void test_closed_form_is_a_selectable_method() {
  printf("closed form as a set_model() method\n");
  Drivers d;
  d.PPFD = 1500.0;

  // --- 1. the exact arm is bit-identical ----------------------------------
  // The one assertion that makes the whole change safe: `optimise()` returns to
  // its previous body when the method is exact, so an explicit "exact" and a
  // two-argument seat must produce the SAME BITS, not merely the same answer.
  {
    phylloptim::Leaf a = make_single_leaf(d, 0.5);
    a.set_model(phylloptim::Leaf::CostCurve::TF24, false);
    a.optimise();
    phylloptim::Leaf b = make_single_leaf(d, 0.5);
    b.set_model("TF24", "stem", "exact");
    b.optimise();
    ok(a.opt_psi_stem_ == b.opt_psi_stem_ && a.ci_ == b.ci_ &&
           a.assim_colimited_ == b.assim_colimited_ && a.profit_ == b.profit_,
       "seating \"exact\" is bit-identical to the two-argument seat");
    ok(b.model_method() == "exact" && a.model_method() == "exact",
       "and both report the exact method");
    ok(!b.last_solve_fell_back_ &&
           !std::isfinite(b.closed_form_fallback_fraction()),
       "an exact solve counts no closed-form call, so phi is NA");
  }

  // --- 2. every output is written, on a REUSED leaf ------------------------
  // The golden files cannot see this class of defect (they build a fresh Leaf per
  // point), so the leaf is reused and the outputs are poisoned first. A field the
  // closed form declined to write would come back as the poison.
  {
    phylloptim::Leaf l = make_single_leaf(d, 0.5);
    l.find_root_collar_psi();          // seat a collar solve's outputs
    const double poison = -12345.0;
    l.ci_ = poison; l.assim_colimited_ = poison; l.transpiration_ = poison;
    l.stom_cond_CO2_ = poison; l.hydraulic_cost_ = poison; l.profit_ = poison;
    l.opt_psi_stem_ = poison;
    l.set_model("TF24", "stem", "closed");
    l.optimise();
    ok(!l.last_solve_fell_back_, "the default operating point uses the fast path");
    ok(l.ci_ != poison && l.assim_colimited_ != poison &&
           l.transpiration_ != poison && l.stom_cond_CO2_ != poison &&
           l.hydraulic_cost_ != poison && l.profit_ != poison &&
           l.opt_psi_stem_ != poison,
       "the closed form writes every rate the exact solve writes");
    // And it clears the collar-route outputs, exactly as the exact stem route
    // does -- those describe a different solve and would otherwise persist.
    ok(!std::isfinite(l.opt_root_psi_) && !std::isfinite(l.E_up_),
       "and clears the collar solve's outputs");
    // The identity that ties the written fields together: profit is the objective
    // evaluated at the reported potential. A path that wrote psi_stem from the
    // closed form and the rates from somewhere else would fail this.
    ok(std::abs(l.profit_ - (l.assim_colimited_ - l.hydraulic_cost_)) < 1e-12,
       "profit_ == assim_colimited_ - hydraulic_cost_ at the reported point");
    ok(l.closed_form_calls_ == 1 && l.closed_form_fallbacks_ == 0 &&
           l.closed_form_fallback_fraction() == 0.0,
       "and the counters record one call and no fallback");
    l.reset_closed_form_counters();
    ok(l.closed_form_calls_ == 0, "the counters are resettable");
  }

  // --- the fallback, and that it lands on the exact answer -----------------
  // Dry soil takes ci/ca below the guard, so the row is handed back. The point of
  // the assertion is that a fallback is not "approximately the exact solve" but
  // IS it, bit for bit: the guard's whole value is that a rejected row costs
  // accuracy nothing.
  {
    phylloptim::Leaf a = make_single_leaf(d, 4.5);
    a.set_model("TF24", "stem", "exact");
    a.optimise();
    phylloptim::Leaf b = make_single_leaf(d, 4.5);
    b.set_model("TF24", "stem", "closed");
    b.optimise();
    ok(b.last_solve_fell_back_, "a dry leaf falls back");
    ok(b.opt_psi_stem_ == a.opt_psi_stem_ && b.assim_colimited_ == a.assim_colimited_,
       "and a fallback IS the exact solve, bit for bit");
    ok(b.closed_form_fallback_fraction() == 1.0, "phi is 1 after one fallback");
  }

  // --- CF77: n = 0, so lambda is the prescribed price exactly --------------
  // The analytical identity for the n = 0 arm: with a constant lambda there is no
  // power law to approximate, so the reported marginal cost must be the price
  // itself. This is exact by construction and is the check that the closed arm
  // reads `cf77_price()` rather than reconstructing a lambda of its own.
  {
    phylloptim::Leaf l = make_single_leaf(d, 0.5);
    l.CF77_lambda_ = 1.5e5;
    l.set_model("CF77", "stem", "closed");
    l.optimise();
    ok(l.lambda_emergent() == 1.5e5,
       "the closed CF77 arm reports the prescribed price exactly");
    // And it composes with the soil-moisture option (#128), which is why the
    // closed form reads cf77_price() and not CF77_lambda_.
    phylloptim::Leaf m = make_single_leaf(d, 0.5);
    m.CF77_lambda_ = 1.5e5;
    m.CF77_soil_beta_ = true;
    m.theta_ = 0.35;
    m.set_model("CF77", "stem", "closed");
    m.optimise();
    const double beta = (0.35 - m.theta_w_) / (m.theta_fc_ - m.theta_w_);
    ok(std::abs(m.lambda_emergent() / (1.5e5 / beta) - 1.0) < 1e-12,
       "and divides it by the Medlyn soil-moisture factor when that is on");
    ok(m.opt_psi_stem_ < l.opt_psi_stem_,
       "so a drier soil closes the stomata on the closed path too");
  }

  // --- 3. the refusals ----------------------------------------------------
  {
    // The collar route: a different model, not a restriction of this one.
    phylloptim::Leaf l = make_single_leaf(d, 0.5);
    l.set_model("TF24", "collar", "closed");
    ok(throws_with([&] { l.optimise(); }, "STEM-route method"),
       "the collar route is refused");

    // Multi-layer: the closed form for it does not exist.
    phylloptim::Leaf m;
    m.setup_transpiration(100);
    m.setup_root_vulnerability(100);
    const std::vector<double> depths{1.0, 2.0, 3.0};
    m.set_physiology(fixture::root_network({0.33, 0.33, 0.33}, depths), 1500.0,
                     {0.5, 0.75, 1.0}, depths, 1.0 * 0.000157 / 5.0, 1.0, 40.0,
                     25.0, 21.0, 101.3);
    m.set_model("TF24", "stem", "closed");
    ok(throws_with([&] { m.optimise(); }, "single-layer only"),
       "a layered soil is refused");

    // The energy balance -- #116, and the refusal states the reason rather than
    // the symptom: D depends on Tleaf depends on E, and E is the unknown.
    phylloptim::Leaf e = make_single_leaf(d, 0.5, true);
    e.set_model("TF24", "stem", "closed");
    ok(throws_with([&] { e.optimise(); }, "not closed"),
       "the energy balance is refused, because the form is then not closed");

    // ⚠️ THE TWO CURVE REFUSALS ARE DIFFERENT CLAIMS. JS22, CMax and ProfitMax
    // all have an h'(A) that does not depend on the solution, so the inversion
    // EXISTS for them and only a dlambda/dpsi is missing. SOX and JW26 have the
    // log link, where h' = 1/A, so lambda carries the assimilation the solve is
    // for: that is a fixed point, and there is nothing to write.
    for (const char *curve : {"JS22", "CMax", "ProfitMax"}) {
      phylloptim::Leaf c = make_single_leaf(d, 0.5);
      c.CF77_lambda_ = 1.5e5;
      c.set_model(curve, "stem", "closed");
      ok(throws_with([&] { c.optimise(); }, "not implemented"),
         std::string("the closed form is unimplemented for ") + curve);
    }
    for (const char *curve : {"SOX", "JW26"}) {
      phylloptim::Leaf c = make_single_leaf(d, 0.5);
      c.set_model(curve, "stem", "closed");
      ok(throws_with([&] { c.optimise(); }, "does not exist"),
         std::string("and cannot exist for ") + curve);
    }

    // An unknown method name is caught at seating, where the caller can still see
    // what they typed.
    phylloptim::Leaf u = make_single_leaf(d, 0.5);
    ok(throws_with([&] { u.set_model("TF24", "stem", "fast"); }, "unknown method"),
       "an unknown method name is refused at set_model");
  }
}

void test_single_layer_optimisers_clear_collar_state() {
  printf("single-layer optimisers do not inherit a collar solve's outputs\n");
  Drivers d;
  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.find_root_collar_psi();
  ok(std::isfinite(l.opt_root_psi_) && std::isfinite(l.E_up_),
     "the collar solve wrote a collar operating point");

  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  ok(!std::isfinite(l.opt_root_psi_), "ProfitMax clears opt_root_psi_");
  ok(!std::isfinite(l.E_up_), "ProfitMax clears E_up_");
  bool consumption_cleared = true;
  for (double c : l.soil_consumption_) {
    if (std::isfinite(c)) consumption_cleared = false;
  }
  ok(consumption_cleared, "ProfitMax clears soil_consumption_");
  ok(std::isfinite(l.transpiration_) && l.transpiration_ > 0.0,
     "while still writing its own transpiration");
}

// `CF77_lambda_` is a caller INPUT -- the Cowan-Farquhar marginal value of water -- and
// nothing in the model writes it. ProfitMax reports the marginal cost its own
// normalisation implies in a SEPARATE field, because the two are not the same
// quantity: one is carbon per unit transpiration, the other per unit conductance.
//
// The load-bearing assertion is the third: a ProfitMax solve followed by a
// Cowan-Farquhar one must price water at the CALLER's number. Before the split it
// priced it at ProfitMax's, silently, because both are finite and plausible.
void test_profitmax_reports_an_emergent_lambda() {
  printf("ProfitMax reports an emergent lambda and leaves the input alone\n");
  Drivers d;
  d.PPFD = 1500.0;
  const double prescribed = 1.5e5;

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.CF77_lambda_ = prescribed;
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();

  ok(l.CF77_lambda_ == prescribed, "the prescribed CF77_lambda_ survives untouched");
  ok(std::isfinite(l.lambda_emergent()) && l.lambda_emergent() > 0.0,
     "and ProfitMax reports the marginal cost its point implies");
  ok(l.lambda_emergent() != prescribed,
     "which is its own number, not the caller's");

  // The hazard the split exists to remove: solving one and then the other now
  // prices water at the caller's value, not at ProfitMax's.
  const double p_after = l.profit_psi_stem_CF77(2.0, 0.5);
  phylloptim::Leaf fresh = make_single_leaf(d, 0.5);
  fresh.CF77_lambda_ = prescribed;
  const double p_clean = fresh.profit_psi_stem_CF77(2.0, 0.5);
  ok(p_after == p_clean,
     "a Cowan-Farquhar solve after ProfitMax is bit-identical to a clean one");
}

// ⚠️ THE TEST THE GOLDEN FILE CANNOT BE. `set_leaf_states_rates_from_psi_stem`
// used to zero transpiration wherever `assim_max_ < 0`, and the golden grid's
// minimum assim_max_ is 3.71, so it never reached the branch. The collar solve
// cannot reach it either -- prepare_collar_solve exits first -- so the ONLY way
// to see this is to call the forward evaluation directly in that regime, which is
// what the single-layer optimisers and the ProfitMax curve do.
void test_transpiration_survives_negative_assim() {
  printf("water moves whether or not there is carbon to be had\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 50.0;  // hot enough that A(ci = ca) cannot cover R_d

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  ok(l.assim_max_ < 0.0, "the regime is reached: assim_max_ is negative");

  const double psi = 3.0;
  l.set_leaf_states_rates_from_psi_stem(psi, 0.5);
  const double E = l.transpiration_;
  ok(E > 0.0, "transpiration follows the hydraulic supply, not the carbon");
  near(E, l.transpiration(psi, 0.5), 1e-12,
       "and equals the supply function exactly");
  ok(l.stom_cond_CO2_ > 0.0, "so the conductance is positive too");

  // The carbon state is the one the branch used to set by hand, and it now comes
  // from the ci solver's own compensation-point fallback.
  near(l.ci_, l.gamma_ * l.umol_per_mol_to_Pa_, 1e-9,
       "ci sits at the compensation point");
  near(l.assim_colimited_, -l.R_d_, 1e-9,
       "and net assimilation is exactly -R_d");
  printf("    Tleaf %.0f C: assim_max_ %.3f, E %.3e kg m-2 s-1, A %.3f\n",
         d.leaf_temp, l.assim_max_, E, l.assim_colimited_);

  // ⚠️ AND TWO LEAVES AT THE SAME OPERATING POINT NOW AGREE ABOUT THE WATER. This
  // is how the old behaviour was found: an arm optimised with respiration off and
  // scored with it on reported a potential that moves water beside a transpiration
  // of exactly zero.
  Drivers dg = d;
  phylloptim::Leaf gross = make_single_leaf(dg, 0.5);
  gross.set_traits(96, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, /*R_d_25=*/0.0, kGammaJS22, kCMaxA, kCMaxB);
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  gross.set_physiology(rn, dg.PPFD, {0.5}, {1.0}, dg.K_s * dg.theta / dg.h,
                       dg.atm_vpd, dg.ca, dg.leaf_temp, dg.atm_o2_kpa,
                       dg.atm_kpa);
  gross.set_leaf_states_rates_from_psi_stem(psi, 0.5);
  ok(gross.assim_max_ > 0.0, "with R_d_25 = 0 the same drivers are NOT shut down");
  near(gross.transpiration_, E, 1e-12,
       "and both leaves report the same transpiration at the same potential");
}

// ⚠️ THE OBJECTIVE IS NOT UNIMODAL AND ITS MAXIMUM CAN BE AN ENDPOINT. This is
// the test that a bare Brent search fails: at a leaf hot enough that net
// assimilation is negative everywhere, the profit is highest at FULL CLOSURE and
// there is a local maximum out in the interior. Brent steps in from the bounds
// and cannot return an endpoint, so it used to report the local one -- an open
// stoma where the model says the leaf should be shut.
void test_profitmax_finds_a_closed_optimum() {
  printf("ProfitMax finds a boundary optimum, which Brent alone cannot\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 50.0;

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.use_thermal_cost_ = true;
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();

  // Reconstruct the objective on a coarse grid and find its global maximum
  // independently of the solver.
  const std::vector<double> curve = l.profitmax_curve(201);
  const std::size_t n = 201;
  std::size_t best = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (curve[4 * n + i] > curve[4 * n + best]) best = i;
  }
  printf("    grid argmax at psi = %.4f (index %zu of %zu), solver returned %.4f\n",
         curve[best], best, n, l.opt_psi_stem_);

  // profitmax_curve re-prepares, so re-solve before reading the operating point.
  l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  const double step = curve[1] - curve[0];
  ok(std::abs(l.opt_psi_stem_ - curve[best]) < 3.0 * step,
     "the solver lands on the objective's GLOBAL maximum, not a local one");

  // And the profit it reports is at least the grid's best.
  ok(l.profit_ >= curve[4 * n + best] - 1e-9,
     "at no lower profit than the grid's best");
}

// util::maximise_over_closed_interval_foc, on functions whose answers are known
// by inspection rather than by running the leaf. Here because EVERY optimiser in
// the package routes through it, so a failure in the leaf tests below should be
// attributable to the leaf rather than to the search.
//
// ⚠️ There used to be a second maximiser, `maximise_over_closed_interval`, which
// was this one without the root-find and which these cases tested instead. It had
// no caller left and is deleted; the cases moved here. Case 7 is what only this
// one can pass, and is the reason the other went.
void test_maximise_over_closed_interval_foc() {
  printf("maximise_over_closed_interval_foc\n");
  const int n = 64;
  const double tol = 1e-12;
  const size_t iters = 200;
  // The derivative contract: `df(x, ok)` sets `*ok` false to mean "no usable
  // derivative here", which is what the leaf's shut-down exits report.
  auto ok_deriv = [](double v) {
    return [v](double, bool* ok) { if (ok) *ok = true; return v; };
  };

  // 1. Maximum AT the left endpoint. A bracketing search cannot return this, and
  //    that is the whole reason this function exists.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return -v; }, ok_deriv(-1.0), 0.0, 1.0, n, tol, iters,
        &fmax);
    ok(x == 0.0, "a maximum at the left endpoint is returned exactly");
    ok(fmax == 0.0, "with its value");
  }
  // 2. ...and at the right.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return v; }, ok_deriv(1.0), 0.0, 1.0, n, tol, iters,
        &fmax);
    ok(x == 1.0, "a maximum at the right endpoint is returned exactly");
    ok(fmax == 1.0, "with its value");
  }
  // 3. An interior maximum is refined, not left on the scan grid. The grid is
  //    64 cells over [0,1] so 0.3 is NOT a grid point; landing within 1e-9 of it
  //    is only possible if the root-find ran.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return -(v - 0.3) * (v - 0.3); },
        [](double v, bool* o) { if (o) *o = true; return -2.0 * (v - 0.3); },
        0.0, 1.0, n, tol, iters, &fmax);
    ok(std::abs(x - 0.3) < 1e-9,
       "an interior maximum is refined off the scan grid");
    ok(fmax <= 0.0 && fmax > -1e-18, "and its value is the peak's");
  }
  // 4. TWO humps, the taller one NOT the one a search from the bounds finds
  //    first. This is the half that endpoints alone do not fix.
  {
    auto two_humps = [](double v) {
      // Peaks at 0.25 (height 1.0) and 0.75 (height 2.0).
      const double a = std::exp(-200.0 * (v - 0.25) * (v - 0.25));
      const double b = 2.0 * std::exp(-200.0 * (v - 0.75) * (v - 0.75));
      return a + b;
    };
    auto d_two_humps = [&](double v, bool* o) {
      if (o) *o = true;
      const double a = std::exp(-200.0 * (v - 0.25) * (v - 0.25));
      const double b = 2.0 * std::exp(-200.0 * (v - 0.75) * (v - 0.75));
      return -400.0 * (v - 0.25) * a - 400.0 * (v - 0.75) * b;
    };
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        two_humps, d_two_humps, 0.0, 1.0, n, tol, iters, &fmax);
    ok(std::abs(x - 0.75) < 1e-5, "the TALLER of two humps is found");
    ok(fmax > 1.9, "and its height is reported");
    // The premise: a bare Brent on the same interval really does miss it, so this
    // case is not vacuous.
    double neg = 0.0;
    const double brent = phylloptim::util::brent_fmin(
        [&](double v) { return -two_humps(v); }, 0.0, 1.0, 1e-3, &neg);
    ok(std::abs(brent - 0.75) > 1e-5,
       "and a bare Brent on the same interval does not");
  }
  // 5. Degenerate inputs: a collapsed interval and a nonsense cell count must
  //    still return a point in range rather than reading off the end of the scan.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return -v; }, ok_deriv(-1.0), 0.5, 0.5, n, tol, iters,
        &fmax);
    ok(x == 0.5, "a collapsed interval returns its one point");
    const double y = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return v; }, ok_deriv(1.0), 0.0, 1.0, 1, tol, iters,
        &fmax);
    ok(y == 1.0, "and n < 2 still compares the endpoints");
  }
  // 6. A non-finite region is skipped rather than selected, and the derivative
  //    REPORTS that rather than returning a number there -- which is the `ok`
  //    flag's whole job.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) {
          return v > 0.6 ? std::numeric_limits<double>::quiet_NaN() : v;
        },
        [](double v, bool* o) {
          if (o) *o = !(v > 0.6);
          return v > 0.6 ? std::numeric_limits<double>::quiet_NaN() : 1.0;
        },
        0.0, 1.0, n, tol, iters, &fmax);
    ok(x <= 0.6 && std::isfinite(fmax), "a NaN region is not selected");
  }
  // 7. ⚠️ THE PROPERTY THE DELETED MAXIMISER COULD NOT HAVE. Brent terminates on
  //    bracket WIDTH, so it leaves the argmax resolved to the tolerance and the
  //    derivative there is merely small. This one solves df == 0, so the returned
  //    point is stationary to solver precision -- which is what a trait gradient
  //    differences, and why a scan-and-refine argmax returned 0.1855 against a
  //    true 0.0551 before this landed.
  {
    auto f  = [](double v) { return -(v - 0.3) * (v - 0.3) * (v - 0.3) * (v - 0.3)
                                    - (v - 0.3) * (v - 0.3); };
    auto df = [](double v, bool* o) {
      if (o) *o = true;
      return -4.0 * (v - 0.3) * (v - 0.3) * (v - 0.3) - 2.0 * (v - 0.3);
    };
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        f, df, 0.0, 1.0, n, tol, iters, &fmax);
    bool o = false;
    const double resid = std::abs(df(x, &o));
    ok(resid < 1e-11, "the returned interior point is STATIONARY, not merely "
                      "inside a narrow bracket");
    // The same objective refined on bracket width instead, for the contrast.
    double neg = 0.0;
    const double width_refined = phylloptim::util::brent_fmin(
        [&](double v) { return -f(v); }, 0.28125, 0.3125,
        (0.3125 - 0.28125) * 1e-4, &neg);
    ok(std::abs(df(width_refined, &o)) > resid,
       "and a bracket-width refinement of the same cell is not");
  }
}

// The two single-layer optimisers reach a maximum that sits AT a bound, and do
// not stop at an interior local one (#94, hazard 11).
//
// ⚠️ WHAT THIS ASSERTS IS "no better point exists", not a particular potential.
// The objective is flat near its maximum, so pinning `opt_psi_stem_` to a
// tolerance would either be loose enough to pass on the defect or tight enough to
// break on a platform's libm. Comparing the returned profit against the best of a
// fine independent scan of the SAME objective is the property, and it is the one
// that fails on a bare bracketing search.
//
// ⚠️ These are single-layer entry points, so they need `set_supply_single`; the
// production collar solve is unaffected and reaches a pinned optimum by its own
// route (`maximise_profit_over_collar`).
void test_single_layer_optimisers_reach_a_bound() {
  printf("single-layer optimisers reach a maximum at a bound\n");

  // Rows chosen so both regimes are present: hot-and-bright, where the objective
  // is maximised at full closure, and mild, where the optimum
  // is interior. A test with only the first would pass on an optimiser that
  // always returned the wet bound.
  struct Row { double PPFD, leaf_temp, psi_soil, atm_vpd; };
  const Row rows[] = {
      {1500.0, 50.0, 0.5, 2.0},   // shut
      {1500.0, 45.0, 1.0, 3.0},   // shut
      {900.0, 25.0, 0.5, 1.5},    // interior
      {900.0, 30.0, 2.0, 2.0},    // interior
      {200.0, 35.0, 3.0, 0.5},
  };

  int pinned_TF = 0, pinned_CF = 0, interior_TF = 0;
  for (const Row &r : rows) {
    Drivers d;
    d.PPFD = r.PPFD;
    d.leaf_temp = r.leaf_temp;
    d.atm_vpd = r.atm_vpd;

    // --- TF24 -------------------------------------------------------------
    {
      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      l.use_thermal_cost_ = true;
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
      const double psi = l.opt_psi_stem_, p = l.profit_;

      // Independent scan of the same objective, on a fresh leaf so the solve's
      // own state cannot influence it.
      phylloptim::Leaf m = make_single_leaf(d, r.psi_soil);
      m.use_thermal_cost_ = true;
      double best = -std::numeric_limits<double>::infinity(), best_psi = r.psi_soil;
      for (int i = 0; i <= 1000; ++i) {
        const double q = r.psi_soil + (m.psi_crit - r.psi_soil) * double(i) / 1000.0;
        const double v = m.profit_psi_stem_TF(q, r.psi_soil);
        if (std::isfinite(v) && v > best) { best = v; best_psi = q; }
      }
      ok(p >= best - 1e-9,
         "the TF24 stem optimum is at no lower profit than a 1001-point scan");
      // And the reported fields describe the returned point rather than the
      // search's last probe (hazard 8).
      ok(l.profit_ == m.profit_psi_stem_TF(psi, r.psi_soil),
         "and its profit is that point's, bit-for-bit");
      if (best_psi <= r.psi_soil + (m.psi_crit - r.psi_soil) * 1e-3) {
        ++pinned_TF;
        ok(psi <= r.psi_soil + (m.psi_crit - r.psi_soil) * 1e-3,
           "where the scan says shut, the optimiser says shut");
      } else {
        ++interior_TF;
      }
    }

    // --- Cowan-Farquhar, deliberately overpriced ---------------------------
    //
    // The prescribed lambda is what can push this model's optimum to a bound: at
    // several times the price the leaf's own operating point implies, water is
    // not worth taking and full closure wins. Derived per row from
    // marginal_cost_water() rather than fixed, because the implied price moves
    // 3.4x across a drydown and a fixed number would leave some rows unpriced.
    {
      phylloptim::Leaf ref = make_single_leaf(d, r.psi_soil);
      ref.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();
      const double lambda = 5.0 * ref.marginal_cost_water();
      if (!(std::isfinite(lambda) && lambda > 0.0)) continue;

      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      l.CF77_lambda_ = lambda;
      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>();
      const double psi = l.opt_psi_stem_, p = l.profit_;

      phylloptim::Leaf m = make_single_leaf(d, r.psi_soil);
      m.CF77_lambda_ = lambda;
      double best = -std::numeric_limits<double>::infinity(), best_psi = r.psi_soil;
      for (int i = 0; i <= 1000; ++i) {
        const double q = r.psi_soil + (m.psi_crit - r.psi_soil) * double(i) / 1000.0;
        const double v = m.profit_psi_stem_CF77(q, r.psi_soil);
        if (std::isfinite(v) && v > best) { best = v; best_psi = q; }
      }
      ok(p >= best - 1e-9,
         "Cowan-Farquhar is at no lower profit than a 1001-point scan");
      ok(l.profit_ == m.profit_psi_stem_CF77(psi, r.psi_soil),
         "and its profit is that point's, bit-for-bit");
      if (best_psi <= r.psi_soil + (m.psi_crit - r.psi_soil) * 1e-3) ++pinned_CF;
    }
  }

  // The premise, so the block above cannot pass by covering only one regime.
  printf("    TF24: %d rows shut, %d interior; Cowan-Farquhar: %d rows shut\n",
         pinned_TF, interior_TF, pinned_CF);
  ok(pinned_TF > 0 && interior_TF > 0,
     "the rows really span both the pinned and the interior regime");
  ok(pinned_CF > 0, "and Cowan-Farquhar reaches the pinned one too");
}

// The two PRODUCT objectives go through the same additive machinery as everything
// else, by storing their cost in log form: `C = -log g`, so `h(A) - C` is
// `log(A*g)`. This asserts the consequence as an IDENTITY at interior points rather
// than as a residual at one argmax, because the residual-at-the-argmax version
// passes on any derivative that merely happens to vanish near there:
//
//     (log P)' * P == P',   P = A*g      <=>   (A'/A + g'/g) * A*g == A'g + A*g'
//
// i.e. the framework's first-order condition is the product rule divided through by
// the strictly positive `A*g`, so the two have identical roots and putting the link
// on the benefit alone drops nothing.
//
// ⚠️ NO FINITE DIFFERENCE HERE, AND THAT IS THE POINT. The first version of this
// test differenced `P` and compared. It passed on macOS at 1.08e-05 and FAILED on
// both Linux compilers at 2.51e-04, because differencing through psi_stem_to_ci's
// root-find has a noise floor of tol/h (hazard 15) and `tol` is set by libm -- so
// the floor is PLATFORM-DEPENDENT and a tolerance measured on one platform does not
// transfer. Every term below is analytic instead: `A'` is a zero-cost JS22
// derivative, `g` and `g'` are the reduction factor and its slope, and the test is
// then pure arithmetic on quantities both sides already computed.
//
// ⚠️ NORMALISED BY THE TERM SCALE, NOT BY THE RESULT. `A'g + A*g'` is a difference
// of two similar magnitudes that cancels to exactly zero at the argmax, so a
// relative error against the result diverges there for a CORRECT implementation.
// Dividing by max(|A'g|, |A*g'|) asks the question that has an answer: are the two
// ways of combining the same terms the same number, to roundoff.
void test_product_link_is_the_product_rule() {
  printf("the product curves' log link is the product rule over A*g\n");

  struct Row { double psi_soil, PPFD, atm_vpd; };
  const Row rows[] = {{0.5, 1500.0, 2.0}, {1.0, 900.0, 2.0},
                      {2.0, 900.0, 4.0}, {1.0, 300.0, 1.0}};
  const int kJS22 = static_cast<int>(phylloptim::Leaf::CostCurve::JS22);
  const int curves[] = {static_cast<int>(phylloptim::Leaf::CostCurve::SOX),
                        static_cast<int>(phylloptim::Leaf::CostCurve::JW26)};

  double worst = 0.0;
  int checked = 0;
  for (const int curve : curves) {
    const bool is_sox = curve == static_cast<int>(phylloptim::Leaf::CostCurve::SOX);
    for (const Row &r : rows) {
      Drivers d;
      d.PPFD = r.PPFD;
      d.atm_vpd = r.atm_vpd;
      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      // A zero unit cost makes the JS22 arm's C' exactly 0, so under the identity
      // link its derivative IS dA/dpsi. Safe as a bare write: JS22_gamma enters
      // one cost expression and no derived state (unlike the curve or capacity
      // traits, which hazard 10 covers).
      l.JS22_gamma = 0.0;
      const double lo = r.psi_soil, hi = l.psi_crit;
      for (int k = 1; k <= 9; ++k) {
        const double psi = lo + (hi - lo) * (0.1 * double(k));
        // P and A first: the derivative calls below re-seat the leaf's state.
        const double P = is_sox ? l.profit_psi_stem_SOX(psi, r.psi_soil)
                                : l.profit_psi_stem_JW26(psi, r.psi_soil);
        const double A = l.assim_colimited_;
        const double g = is_sox ? l.sox_reduction(psi) : l.jw26_reduction(psi);
        const double gp =
            is_sox ? l.sox_reduction_deriv(psi) : l.jw26_reduction_deriv(psi);
        // The premise of the log form. Where it fails nothing is claimed: at full
        // closure A = -R_d < 0 and log A does not exist, which is precisely why the
        // OPTIMISER searches the product instead of its logarithm.
        if (!(P > 0.0) || !(g > 0.0)) {
          continue;
        }
        l.set_model(static_cast<phylloptim::Leaf::CostCurve>(curve), false);
        const std::vector<double> dlog = l.dprofit_dpsi_stem_checked(psi);
        l.set_model(static_cast<phylloptim::Leaf::CostCurve>(kJS22), false);
        const std::vector<double> dA = l.dprofit_dpsi_stem_checked(psi);
        if (!std::isfinite(dlog[0]) || !std::isfinite(dA[0])) {
          continue;
        }
        const double lhs = dlog[0] * P;          // (log P)' * P
        const double rhs = dA[0] * g + A * gp;   // the product rule
        const double scale = std::max(std::fabs(dA[0] * g), std::fabs(A * gp));
        const double rel = std::fabs(lhs - rhs) / std::max(scale, 1e-30);
        if (rel > worst) {
          worst = rel;
        }
        ++checked;
      }
    }
  }
  // ⚠️ The count is asserted because the two `continue`s above can empty this test
  // silently -- a curve returning NaN everywhere would leave `worst` at 0 and pass
  // (hazard 16).
  printf("    %d interior points | worst departure from (logP)'*P == A'g + Ag': "
         "%.3e\n", checked, worst);
  ok(checked >= 40, "the identity was actually evaluated, on both curves");
  ok(worst < 1e-12, "the log-link derivative IS the product rule divided by A*g");
}

// EVERY curve's returned point is the actual maximum of ITS OWN objective, checked
// against a derivative-free scan of that objective. This is the optimality test the
// suite was missing for four of the seven: before it, `TF24`, `CF77` and `ProfitMax`
// had scan oracles and `JS22`, `CMax`, `SOX` and `JW26` had only
// `psi_stem_optima.tsv` -- which is a REGRESSION baseline and would freeze a wrong
// optimum exactly as happily as a right one.
//
// A scan is the right oracle rather than a first-order check, because a maximum here
// is often at a bound: `dJ/dpsi != 0` at a pinned optimum, so an FOC test would
// report a correct answer as wrong (hazard 11). Comparing OBJECTIVE VALUES asks the
// question that is true in both regimes -- is there any feasible point that beats
// the one returned.
//
// ⚠️ THE SCAN CAN ONLY EVER TIE OR LOSE, which is what makes the direction of the
// assertion meaningful: the solver includes both endpoints and refines within a
// cell, so it should never come back below a grid it contains. A scan that WINS is
// the failure -- it means the search missed a basin.
// ============================================================================
// THE TWO-TERM COST: `Theta(E) = Theta~(psi) + lambda_o*E`
// ----------------------------------------------------------------------------
// The decomposition is canonical -- any cost splits into a part depending on the
// potential alone and a linear price of water -- and the content of this curve is
// that the second part is not zero. Every conductance-loss curve in this package
// has `lambda -> 0` as the drop closes, i.e. prices water as free precisely when
// it is abundant; `lambda_o` is the floor that removes that.
//
// `Theta~` here is TF24's OWN cost, read at TF24's own traits, so this curve is
// TF24 plus CF77's price and has exactly ONE parameter of its own. That is what
// makes the first test below the strongest available: two reductions, asserted
// bit-for-bit, needing no fixture and no tolerance -- and the first of them is a
// ONE-restriction statement about the production model.
// ============================================================================

// Bit-equality that survives NaN, which several reported fields legitimately are.
bool same_bits(double a, double b) {
  return (std::isnan(a) && std::isnan(b)) || a == b;
}

// The reported operating point, as one vector, so a reduction can be compared
// field by field rather than on `profit_` alone -- a curve wired to the right
// objective and the wrong output assembly would pass the latter.
std::vector<double> reported(const phylloptim::Leaf &l) {
  return {l.opt_psi_stem_, l.profit_, l.hydraulic_cost_, l.ci_,
          l.assim_colimited_, l.transpiration_, l.stom_cond_CO2_,
          l.lambda_emergent()};
}

const char *reported_names[] = {"opt_psi_stem", "profit", "hydraulic_cost",
                                "ci", "assim", "transpiration", "gc",
                                "lambda_emergent"};

void test_shadow_cost() {
  printf("shadow_cost separates the price from the realised carbon cost\n");
  // WHAT THIS PINS. `hydraulic_cost_` is the whole cost the objective subtracted;
  // `shadow_cost()` is the share of it that is a price rather than carbon the
  // plant gave up. So the carbon actually forgone is the difference, and the
  // realised profit is `profit + shadow_cost`. A consumer that grows a plant on
  // this leaf wants that second quantity and never `profit`, or it taxes growth
  // by something the plant never spent.
  //
  // ⚠️ BIT EQUALITY against `lambda_o * transpiration_`, not a tolerance. The
  // accessor reads the stored transpiration rather than recomputing it, precisely
  // so this identity holds exactly against the REPORTED `E` -- which is what a
  // caller checks. A tolerance would pass on an accessor that recomputed
  // transpiration from the potentials and drifted in the last bit.
  Drivers d;
  int compared = 0;
  for (double ppfd : {300.0, 900.0, 1500.0}) {
    for (double vpd : {1.0, 2.0, 4.0}) {
      for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
        d.PPFD = ppfd;
        d.atm_vpd = vpd;

        for (double lo : {0.0, 1.5e4, 1.5e5}) {
          phylloptim::Leaf l = make_single_leaf(d, psi_soil);
          l.TF24_floor_lambda_o = lo;
          l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
          l.optimise();
          if (!std::isfinite(l.transpiration_)) continue;

          ok(l.shadow_cost() == lo * l.transpiration_,
             "shadow_cost is lambda_o * E exactly");
          // the realised cost is what is left, and it is never the whole cost
          // unless the price is zero
          const double realised = l.hydraulic_cost_ - l.shadow_cost();
          ok(std::isfinite(realised) && realised >= 0.0,
             "the realised cost is non-negative");
          ok((lo == 0.0) == (l.shadow_cost() == 0.0),
             "the price vanishes exactly when lambda_o does");
          ++compared;
        }

        // ⚠️ ZERO ON EVERY OTHER CURVE, INCLUDING CF77. That is a statement about
        // the curve, not about the cost: CF77's whole cost IS `lambda * E` and
        // reading it as a shadow price is ordinary, but the model supplies one
        // number and nothing to attribute it with. Reporting it as all-shadow
        // would ship one reading as a fact. See `Leaf::shadow_cost`.
        for (auto c : {phylloptim::Leaf::CostCurve::TF24,
                       phylloptim::Leaf::CostCurve::CF77,
                       phylloptim::Leaf::CostCurve::JS22,
                       phylloptim::Leaf::CostCurve::CMax}) {
          phylloptim::Leaf l = make_single_leaf(d, psi_soil);
          l.CF77_lambda_ = 1.5e5;
          l.set_model(c, false);
          l.optimise();
          ok(l.shadow_cost() == 0.0,
             "shadow_cost is exactly zero off TF24_floor");
        }
      }
    }
  }
  printf("  compared %d priced operating points\n", compared);
}


void test_lambda_emergent_reported() {
  printf("operating_point_values() reports the SEATED curve's lambda\n");
  // WHY THIS EXISTS. `lambda` is `marginal_cost_water()` -- TF24's price at this
  // operating point, whatever curve ran -- so on any other curve it answers a
  // question nobody asked. `lambda_emergent` is the seated curve's own
  // (dC/dpsi)/(dE/dpsi), and it was reachable only off the object until it was
  // appended here, which left every leaf_solve() caller reconstructing it.
  //
  // ⚠️ CF77 IS THE CASE THAT PINS IT, because its emergent lambda is KNOWN: the
  // curve prices water at a prescribed constant, so `lambda_emergent` must come
  // back as that constant. Nothing else here has an answer fixed independently of
  // the solve.
  Drivers d;
  const double price = 1.5e5;
  int compared = 0;
  for (double ppfd : {300.0, 900.0, 1500.0}) {
    for (double psi_soil : {0.5, 1.0, 2.0}) {
      d.PPFD = ppfd;

      phylloptim::Leaf c = make_single_leaf(d, psi_soil);
      c.CF77_lambda_ = price;
      c.set_model(phylloptim::Leaf::CostCurve::CF77, false);
      c.optimise();
      if (std::isfinite(c.transpiration_) && c.transpiration_ > 0.0) {
        ok(std::fabs(c.lambda_emergent() - price) <= 1e-6 * price,
           "CF77's emergent lambda is the price it was given");
        ++compared;
      }

      // TF24 is the one curve where the two agree, by construction
      phylloptim::Leaf t = make_single_leaf(d, psi_soil);
      t.set_model(phylloptim::Leaf::CostCurve::TF24, false);
      t.optimise();
      ok(t.lambda_emergent() == t.marginal_cost_water(),
         "on TF24 the seated lambda and marginal_cost_water agree exactly");

      // and TF24_floor's is its hydraulic half plus the floor -- the correction a
      // caller used to have to make by hand
      phylloptim::Leaf f = make_single_leaf(d, psi_soil);
      f.TF24_floor_lambda_o = price;
      f.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
      f.optimise();
      ok(std::fabs(f.lambda_emergent() -
                   (f.marginal_cost_water() + price)) <= 1e-9 * price,
         "TF24_floor's emergent lambda is the hydraulic half plus the floor");

      // the reported vector must carry the same number the accessor does
      const std::vector<double> ov = f.operating_point_values();
      ok(ov.back() == f.lambda_emergent(),
         "operating_point_values() reports lambda_emergent last");
    }
  }
  printf("  %d CF77 points with a known emergent lambda\n", compared);
}


void test_tf24_floor_reduces() {
  printf("TF24_floor reduces to TF24 and to CF77, bit-for-bit\n");
  // ⚠️ BIT EQUALITY IS THE RIGHT TEST HERE AND A TOLERANCE WOULD BE THE WRONG
  // ONE. The reductions are exact by construction rather than by convergence:
  // each term is the parent's own expression, so zeroing one parameter adds an
  // exact zero to the other curve's exact value -- see
  // `Leaf::hydraulic_cost_TF24_floor`, which also records why fused multiply-add
  // does not break it. A tolerance would pass on a rewrite that lost the
  // property, which is the only thing this test exists to catch.
  Drivers d;
  int compared = 0;
  for (double ppfd : {300.0, 900.0, 1500.0}) {
    for (double vpd : {1.0, 2.0, 4.0}) {
      for (double psi_soil : {0.5, 1.0, 2.0, 3.0, 4.0}) {
        d.PPFD = ppfd;
        d.atm_vpd = vpd;

        // --- lambda_o = 0 is TF24, at TF24's own traits --------------------
        // ⚠️ THE ONE-RESTRICTION REDUCTION, and the reason this curve has no
        // hydraulic parameter of its own. Nothing is set on either leaf but the
        // price: if the two disagree, they disagree about `lambda_o` alone.
        {
          phylloptim::Leaf a = make_single_leaf(d, psi_soil);
          a.TF24_floor_lambda_o = 0.0;
          a.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
          a.optimise();

          phylloptim::Leaf b = make_single_leaf(d, psi_soil);
          b.set_model(phylloptim::Leaf::CostCurve::TF24, false);
          b.optimise();

          const std::vector<double> va = reported(a), vb = reported(b);
          for (std::size_t i = 0; i < va.size(); ++i) {
            ok(same_bits(va[i], vb[i]),
               std::string("TF24_floor(lambda_o=0) == TF24 in ") +
                   reported_names[i] + " at psi_soil=" +
                   std::to_string(psi_soil) + " ppfd=" + std::to_string(ppfd) +
                   " vpd=" + std::to_string(vpd));
          }
        }

        // --- TF24_cost_scale = 0 is CF77 -----------------------------------
        // The hydraulic half is scaled by `TF24_cost_scale`, so zeroing that
        // trait is what leaves the price alone. `pow(1-f, beta2)` is finite
        // inside the bracket, so `0 * pow(...)` is an exact zero rather than a
        // NaN -- which is why this reduction is available at all.
        {
          phylloptim::Leaf a = make_single_leaf(d, psi_soil);
          a.TF24_cost_scale = 0.0;
          a.TF24_floor_lambda_o = kFloorLambdaO;
          a.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
          a.optimise();

          phylloptim::Leaf b = make_single_leaf(d, psi_soil);
          b.CF77_lambda_ = kFloorLambdaO;
          b.set_model(phylloptim::Leaf::CostCurve::CF77, false);
          b.optimise();

          const std::vector<double> va = reported(a), vb = reported(b);
          for (std::size_t i = 0; i < va.size(); ++i) {
            ok(same_bits(va[i], vb[i]),
               std::string("TF24_floor(TF24_cost_scale=0) == CF77 in ") +
                   reported_names[i] + " at psi_soil=" +
                   std::to_string(psi_soil) + " ppfd=" + std::to_string(ppfd) +
                   " vpd=" + std::to_string(vpd));
          }
        }
        ++compared;
      }
    }
  }
  // ⚠️ Hazard 16: an assertion that never runs passes. Both reductions are
  // inside two loops and a refused row would silently skip one.
  printf("    %d driver rows, both reductions, %d fields each\n", compared,
         int(sizeof(reported_names) / sizeof(reported_names[0])));
  ok(compared == 45, "every driver row was compared");
}

void test_tf24_floor_wet_end_price() {
  printf("TF24_floor's marginal cost of water does not vanish at the wet end\n");
  Drivers d;
  phylloptim::Leaf l = make_single_leaf(d, 1.0);
  l.TF24_floor_lambda_o = kFloorLambdaO;

  // ⚠️ THE LIMIT IS AT ZERO POTENTIAL, NOT AT ZERO DROP, and that is the one
  // behavioural difference from the JS22-based version of this curve. The
  // hydraulic half reads the ABSOLUTE potential, so it vanishes as psi -> 0 --
  // not as psi -> psi_soil, where TF24's own cost is already non-zero.
  double prev_excess = 0.0;
  for (int k = 0; k < 4; ++k) {
    const double psi = 1e-2 / std::pow(2.0, k);
    const double excess = l.lambda_TF24_floor(psi) - l.TF24_floor_lambda_o;
    ok(excess > 0.0, "the excess over lambda_o is positive");
    if (k > 0) {
      // lambda_TF24 ~ psi^(stem_c*beta2 - 1), so halving psi divides the excess
      // by 2^(stem_c*beta2 - 1). Asserting the EXPONENT rather than a ratio of 2
      // is what makes this a statement about the curve rather than about a step.
      const double n = l.stem_c * l.TF24_beta2 - 1.0;
      near(prev_excess / excess, std::pow(2.0, n), 1e-2,
           "halving psi divides the excess by 2^(stem_c*beta2 - 1)");
    }
    prev_excess = excess;
  }
  ok(l.stem_c * l.TF24_beta2 > 1.0,
     "and the wet-end limit exists at all, which needs stem_c*beta2 > 1");

  // ⚠️ THE CONTRAST IS THE POINT. At the same potential every conductance-loss
  // curve here is heading for zero, which is what "water is free when abundant"
  // means as a number rather than as an argument.
  const double psi = 1e-6;
  printf("    at psi = %.0e MPa: TF24_floor %.6g, TF24 %.3e, JS22 %.3e\n", psi,
         l.lambda_TF24_floor(psi), l.lambda_TF24(psi), l.lambda_JS22(psi, 0.0));
  ok(l.lambda_TF24(psi) < 1e-3 * l.TF24_floor_lambda_o,
     "TF24's price at the same potential is negligible beside lambda_o");
  near(l.lambda_TF24_floor(psi) / l.TF24_floor_lambda_o, 1.0, 1e-6,
     "so TF24_floor's is lambda_o to six figures");
}

void test_tf24_floor_emergent_lambda() {
  printf("TF24_floor's emergent lambda is dC/dE at the solved point\n");
  // Against a finite difference of the COST over the TRANSPIRATION, not against
  // a second copy of `lambda_TF24_floor` -- which would only assert that the
  // dispatch table reaches the function it reaches. This tests the definition.
  Drivers d;
  d.PPFD = 1500.0;
  const double h = 1e-3;   // the step `test_every_curve_reports_an_emergent_lambda` uses
  int interior = 0;

  for (double psi_soil : {0.5, 1.0, 2.0}) {
    phylloptim::Leaf l = make_single_leaf(d, psi_soil);
    l.TF24_floor_lambda_o = kFloorLambdaO;
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    l.optimise();
    const double psi = l.opt_psi_stem_;
    const double got = l.lambda_emergent();
    if (psi <= psi_soil + h || psi >= l.psi_crit - h) {
      continue;   // a pinned optimum: dC/dE is defined but is not what pinned it
    }
    ++interior;

    // ⚠️ THROUGH `profit_psi_stem_TF24_floor`, NOT `hydraulic_cost_TF24_floor`. The cost
    // function reads the memoised transpiration but does not SEAT the leaf, so
    // reading `transpiration_` after it returns gives whatever the last solve
    // left there -- the same E at both ends, a zero denominator, and a ratio of
    // exactly 0 that looks like a broken lambda. The profit function seats it.
    l.profit_psi_stem_TF24_floor(psi + h, psi_soil);
    const double C1 = l.hydraulic_cost_, E1 = l.transpiration_;
    l.profit_psi_stem_TF24_floor(psi - h, psi_soil);
    const double C0 = l.hydraulic_cost_, E0 = l.transpiration_;
    near(got / ((C1 - C0) / (E1 - E0)), 1.0, 1e-4,
         "TF24_floor's emergent lambda is dC/dE at psi_soil=" +
             std::to_string(psi_soil));

    // And it exceeds lambda_o by exactly the hydraulic term, which is the
    // decomposition read back off the solved point.
    ok(got > l.TF24_floor_lambda_o,
       "the emergent lambda sits above the wet-end floor at psi_soil=" +
           std::to_string(psi_soil));
    near(got - l.TF24_floor_lambda_o, l.lambda_TF24(psi), 1e-12,
       "and the excess over the floor IS TF24's own lambda at psi_soil=" +
           std::to_string(psi_soil));
  }
  // Hazard 16 again: print what was actually exercised.
  printf("    %d of 3 rows interior and compared\n", interior);
  ok(interior >= 2, "the identity was evaluated on at least two rows");
}

void test_tf24_floor_kmax_vpd_invariance() {
  printf("TF24_floor sits between TF24 and CF77 under kmax/D rescaling\n");
  // ⚠️ A PREDICTION, NOT A REGRESSION CHECK. Scaling `kmax` and the vapour
  // deficit together leaves the stomatal conductance -- and so ci, A and dA/dpsi
  // -- untouched at every potential, because `gc ~ kmax*G(psi)/D`. So a cost
  // whose MARGINAL value depends on psi alone is invariant, and one that prices
  // the FLUX is not: `dC/dpsi` is TF24's `|f'|` form for the first part and
  // `lambda_o*kmax*f(psi)` for the second. TF24_floor carries both, so it must move
  // by less than CF77 and by more than TF24, and by more as lambda_o rises.
  const double factor = 16.0;
  Drivers lo, hi;
  lo.PPFD = hi.PPFD = 1500.0;
  lo.atm_vpd = 2.0;
  hi.atm_vpd = lo.atm_vpd * factor;
  hi.K_s = lo.K_s * factor;
  const double psi_soil = 1.0;

  auto solved = [&](const Drivers &d, phylloptim::Leaf::CostCurve curve,
                    double lambda_o) {
    phylloptim::Leaf l = make_single_leaf(d, psi_soil);
    l.CF77_lambda_ = kFloorLambdaO;
    l.TF24_floor_lambda_o = lambda_o;
    l.set_model(curve, false);
    l.optimise();
    return l.opt_psi_stem_;
  };
  auto shift = [&](phylloptim::Leaf::CostCurve curve, double lambda_o) {
    const double a = solved(lo, curve, lambda_o);
    const double b = solved(hi, curve, lambda_o);
    return std::abs(b - a) / a;
  };

  using Curve = phylloptim::Leaf::CostCurve;
  const double tf24 = shift(Curve::TF24, 0.0);
  const double cf77 = shift(Curve::CF77, 0.0);
  printf("    kmax and D both x%.0f:  TF24 %.3e  CF77 %.3e\n", factor, tf24,
         cf77);
  ok(tf24 < 1e-12, "TF24's optimum is invariant");
  ok(cf77 > 0.1, "CF77's is not, because it prices the flux");

  // The interpolation, and that it is monotone in lambda_o.
  double previous = 0.0;
  for (double lambda_o : {1.5e3, 1.5e4, 1.5e5}) {
    const double s = shift(Curve::TF24_floor, lambda_o);
    printf("    TF24_floor at lambda_o = %.1e: %.3e\n", lambda_o, s);
    ok(s > tf24, "TF24_floor moves more than TF24 at lambda_o=" +
                     std::to_string(lambda_o));
    ok(s < cf77, "and less than CF77 at lambda_o=" + std::to_string(lambda_o));
    ok(s > previous, "and by more as lambda_o rises");
    previous = s;
  }
}

void test_tf24_floor_refuses_an_unset_price() {
  printf("TF24_floor refuses an unset price, and says why\n");
  Drivers d;
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    ok(throws_with([&] { l.optimise(); },
                   "needs TF24_floor_lambda_o set"),
       "an unset lambda_o is refused rather than optimised as a NaN objective");
    // ⚠️ AND THE MESSAGE HAS TO NAME TF24, or the obvious repair -- set it to
    // zero -- silently substitutes the production model for the one the caller
    // asked for.
    ok(throws_with([&] { l.optimise(); }, "would silently make this curve TF24"),
       "and says why there is no default");
    // It also has to say where a price can be SET FROM, since the field is not
    // reachable through leaf_traits(). That gap is what made this curve
    // unusable through leaf_solve() when it was first added.
    ok(throws_with([&] { l.optimise(); }, "leaf_solve()"),
       "and names a route that can supply one");
  }
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.TF24_floor_lambda_o = 0.0;
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    l.optimise();
    ok(std::isfinite(l.profit_), "an explicit lambda_o = 0 solves");
  }
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.TF24_floor_lambda_o = -1.0;
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    ok(throws_with([&] { l.optimise(); }, "non-negative"),
       "a negative price -- which would PAY the leaf to transpire -- is refused");
  }
}

void test_tf24_floor_closed_form() {
  printf("TF24_floor's closed form: the reduction, the coverage, and the tail\n");
  Drivers d;
  d.PPFD = 1500.0;

  // --- at TF24_cost_scale = 0 this curve is CF77, and the two closed arms are
  // compared. ⚠️ NOT BIT-FOR-BIT, AND THE DIFFERENCE IS DELIBERATE.
  // `solve_CF77` takes no Newton step at all -- with lambda constant there is no
  // circularity to break -- and stops at the spline inversion of the supply. This
  // arm still takes its step, and with a zeroed hydraulic term there is nothing in
  // lambda for it to correct, so what it does instead is drive the
  // supply-minus-demand residual to zero from a start clamped strictly inside the
  // bracket. So the comparison worth making is against the EXACT answer, and the
  // finding is that the extra step does not hurt.
  int both_served = 0, different_decision = 0;
  double worst_gap = 0.0;
  for (double psi_soil : {0.5, 1.0, 2.0}) {
    phylloptim::Leaf exact = make_single_leaf(d, psi_soil);
    exact.TF24_cost_scale = 0.0;
    exact.TF24_floor_lambda_o = 1.5e5;
    exact.set_model("TF24_floor", "stem", "exact");
    exact.optimise();

    phylloptim::Leaf a = make_single_leaf(d, psi_soil);
    a.TF24_cost_scale = 0.0;
    a.TF24_floor_lambda_o = 1.5e5;
    a.set_model("TF24_floor", "stem", "closed");
    a.optimise();

    phylloptim::Leaf b = make_single_leaf(d, psi_soil);
    b.CF77_lambda_ = 1.5e5;
    b.set_model("CF77", "stem", "closed");
    b.optimise();

    // The EXACT solves must agree bit-for-bit -- that is the reduction, on the
    // method both curves share, and it is asserted without a tolerance.
    phylloptim::Leaf cf_exact = make_single_leaf(d, psi_soil);
    cf_exact.CF77_lambda_ = 1.5e5;
    cf_exact.set_model("CF77", "stem", "exact");
    cf_exact.optimise();
    ok(same_bits(exact.opt_psi_stem_, cf_exact.opt_psi_stem_),
       "exact TF24_floor(TF24_cost_scale=0) IS exact CF77 at psi_soil=" +
           std::to_string(psi_soil));

    if (a.last_solve_fell_back_ != b.last_solve_fell_back_) {
      ++different_decision;
      continue;
    }
    if (a.last_solve_fell_back_) {
      continue;
    }
    ++both_served;
    const double ea = std::abs(a.opt_psi_stem_ - exact.opt_psi_stem_);
    const double eb = std::abs(b.opt_psi_stem_ - exact.opt_psi_stem_);
    ok(ea <= eb * 1.05,
       "the TF24_floor arm is no further from the exact optimum than CF77's at "
       "psi_soil=" + std::to_string(psi_soil));
    worst_gap = std::max(worst_gap,
                         std::abs(a.opt_psi_stem_ - b.opt_psi_stem_) /
                             b.opt_psi_stem_);
  }
  printf("    scale=0: %d rows served by both (worst |dpsi*| between arms %.3e), "
         "%d where only one was served\n", both_served, worst_gap,
         different_decision);
  ok(both_served >= 2, "the scale = 0 comparison was actually exercised");

  // --- coverage and the worst served error at real parameters -----------------
  // ⚠️ THIS BOUND IS A MEASUREMENT, NOT A SPECIFICATION. What limits the accuracy
  // is the USO collapse (see closed_form.hpp), and on this curve the leading order
  // additionally discards the whole TF24 term -- so the tail is worse than CF77's
  // and the shared ci/ca guard does not see it. The assertion exists to notice
  // DRIFT; read the printed numbers, not the pass.
  int served = 0, fell_back = 0;
  double worst = 0.0;
  for (double scale : {0.5, 7.5}) {
    for (double psi_soil : {0.5, 1.0, 1.5, 2.0}) {
      for (double vpd : {1.0, 2.0}) {
        Drivers dd = d;
        dd.atm_vpd = vpd;
        phylloptim::Leaf exact = make_single_leaf(dd, psi_soil);
        exact.TF24_cost_scale = scale;
        exact.TF24_floor_lambda_o = 1.5e5;
        exact.set_model("TF24_floor", "stem", "exact");
        exact.optimise();

        phylloptim::Leaf closed = make_single_leaf(dd, psi_soil);
        closed.TF24_cost_scale = scale;
        closed.TF24_floor_lambda_o = 1.5e5;
        closed.set_model("TF24_floor", "stem", "closed");
        closed.optimise();

        if (closed.last_solve_fell_back_) {
          ++fell_back;
          // ⚠️ A FALLBACK IS THE WHOLE EXACT SOLVE, so it must be bit-identical
          // rather than close. This is what makes the fallback fraction a cost
          // rather than an accuracy question.
          ok(same_bits(closed.opt_psi_stem_, exact.opt_psi_stem_) &&
                 same_bits(closed.profit_, exact.profit_),
             "a fallback IS the exact solve, bit-for-bit");
          continue;
        }
        ++served;
        ok(closed.opt_psi_stem_ > psi_soil && closed.opt_psi_stem_ < closed.psi_crit,
           "a served row is strictly inside the bracket");
        worst = std::max(worst,
                         std::abs(closed.opt_psi_stem_ - exact.opt_psi_stem_) /
                             exact.opt_psi_stem_);
      }
    }
  }
  printf("    %d served, %d fell back (phi %.2f) | worst served |dpsi*|/psi* %.3e\n",
         served, fell_back, double(fell_back) / double(served + fell_back), worst);
  ok(served >= 2, "the closed form served at least two rows");
  ok(worst < 0.6, "and no served row is worse than the measured tail");

  // --- and where it is refused -------------------------------------------------
  // At lambda_o = 0 the curve is TF24, whose wet-end lambda is zero rather than a
  // constant: `xi = sqrt(Q/lambda)` diverges and there is nothing to expand about.
  // Refused, naming which of the two curves it has become -- and TF24 has its own
  // start, so the message can point somewhere.
  phylloptim::Leaf zero = make_single_leaf(d, 1.0);
  zero.TF24_floor_lambda_o = 0.0;
  zero.set_model("TF24_floor", "stem", "closed");
  ok(throws_with([&] { zero.optimise(); }, "needs TF24_floor_lambda_o > 0"),
     "the closed form refuses lambda_o = 0");
  ok(throws_with([&] { zero.optimise(); }, "Seat TF24 instead"),
     "and points at the curve that does have a start");
}

void test_every_curve_returns_its_own_maximum() {
  printf("every cost curve's optimum beats a scan of its own objective\n");

  // ⚠️ THE HOT ROWS ARE NOT PADDING. A boundary optimum is where a scan oracle
  // earns its keep: an optimiser that steps in from the bounds reports an open
  // stoma where the objective says shut, and only a comparison of VALUES over a
  // closed interval catches it (hazard 11). Without rows hot enough to make
  // assimilation negative across the supply stream, this test would exercise the
  // interior branch alone -- which is the half that was never in doubt.
  struct Row { double psi_soil, PPFD, atm_vpd, leaf_temp; };
  const Row rows[] = {
      {0.5, 1500.0, 2.0, 25.0},
      {1.0, 900.0, 2.0, 25.0},
      {2.0, 900.0, 4.0, 30.0},
      {3.0, 1500.0, 2.0, 25.0},
      {1.0, 300.0, 1.0, 25.0},
      {0.5, 1500.0, 2.0, 50.0},   // shut
      {1.0, 1500.0, 3.0, 45.0},   // shut
      {2.0, 1200.0, 4.0, 48.0},   // shut
  };
  const int n_curves = phylloptim::Leaf::n_cost_curves;
  const int N = 20000;

  int total_interior = 0, total_pinned = 0, total_rows = 0;
  int min_rows_per_curve = std::numeric_limits<int>::max();
  double worst_shortfall = 0.0;
  std::string worst_where;
  for (int curve = 0; curve < n_curves; ++curve) {
    double curve_worst = 0.0;
    int interior = 0, pinned = 0, curve_rows = 0;
    for (const Row &r : rows) {
      Drivers d;
      d.PPFD = r.PPFD;
      d.atm_vpd = r.atm_vpd;
      d.leaf_temp = r.leaf_temp;

      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      // The two curves with no default: prescribed, not derived, and unset is NA.
      // Set unconditionally, because a field a curve does not read costs nothing.
      l.CF77_lambda_ = 1.5e5;
      l.TF24_floor_lambda_o = kFloorLambdaO;
      l.set_model(static_cast<phylloptim::Leaf::CostCurve>(curve), false);
      l.optimise();
      const double got = l.profit_;
      if (!std::isfinite(got)) {
        continue;   // a refused row; psi_stem_optima.tsv covers the refusals
      }

      // Fresh leaf, so the solve's own state cannot influence the oracle.
      phylloptim::Leaf m = make_single_leaf(d, r.psi_soil);
      m.CF77_lambda_ = 1.5e5;
      m.TF24_floor_lambda_o = kFloorLambdaO;
      // ⚠️ ProfitMax's objective is UNDEFINED on a leaf that has not scanned for
      // |A|max, so a fresh oracle leaf returns NaN at every point and the curve
      // silently contributes nothing -- which is how this test first passed while
      // testing six curves and reporting seven. Seed it with the SAME normaliser
      // the solver used: a different |A|max is a different objective, so comparing
      // against one would not be comparing against this curve at all.
      if (curve == static_cast<int>(phylloptim::Leaf::CostCurve::ProfitMax)) {
        m.prepare_profitmax_at(l.profitmax_A_max_);
      }
      m.set_model(static_cast<phylloptim::Leaf::CostCurve>(curve), false);
      const double lo = r.psi_soil, hi = m.psi_crit;
      double best = -std::numeric_limits<double>::infinity(), best_psi = lo;
      for (int i = 0; i <= N; ++i) {
        const double q = lo + (hi - lo) * (double(i) / double(N));
        const double v = m.evaluate_psi_stem_at(q);
        if (std::isfinite(v) && v > best) { best = v; best_psi = q; }
      }
      if (!std::isfinite(best)) {
        continue;
      }
      ++total_rows;
      ++curve_rows;
      // Relative slack, because the eight objectives differ in units and by orders
      // of magnitude -- ProfitMax is dimensionless and O(1), TF24 is carbon.
      const double shortfall =
          (best - got) / std::max(std::fabs(best), 1.0);
      if (shortfall > curve_worst) {
        curve_worst = shortfall;
        if (shortfall > worst_shortfall) {
          worst_shortfall = shortfall;
          worst_where = phylloptim::Leaf::curve_name(curve) + " psi_soil=" +
                        std::to_string(r.psi_soil);
        }
      }
      const double edge = (hi - lo) * 1e-3;
      if (best_psi <= lo + edge || best_psi >= hi - edge) { ++pinned; } else { ++interior; }
    }
    printf("    %-10s %d rows | worst shortfall %.3e | %d interior, %d pinned\n",
           phylloptim::Leaf::curve_name(curve).c_str(), curve_rows, curve_worst,
           interior, pinned);
    if (curve_rows < min_rows_per_curve) { min_rows_per_curve = curve_rows; }
    total_interior += interior;
    total_pinned += pinned;
  }

  // The premise, asserted so the block cannot pass by covering one regime or by
  // quietly evaluating nothing (hazard 16).
  printf("    %d rows | %d interior, %d pinned | worst %.3e at %s\n", total_rows,
         total_interior, total_pinned, worst_shortfall, worst_where.c_str());
  // ⚠️ PER CURVE, not in total. A total is satisfied by the other six while one
  // contributes nothing, which is exactly what happened here before the ProfitMax
  // seeding above: 30 rows from six curves passed a `>= 30` check.
  ok(min_rows_per_curve > 0, "EVERY curve contributed solved rows, not just six");
  ok(total_interior > 0 && total_pinned > 0,
     "and the rows span both the interior and the pinned regime");
  ok(worst_shortfall < 1e-6,
     "no curve's solver is beaten by a 20001-point scan of its own objective");
}

// Each curve's reported `profit_` is ITS OWN objective at the point it returned.
//
// ⚠️ THE GOLDEN GRID CANNOT SEE THIS, WHICH IS WHY IT IS HERE. Every golden row on
// the collar route solves TF24, so a bug that reports TF24's objective for a
// different curve is invisible to all four baselines -- and that bug existed: two
// `profit_psi_stem_TF` calls survived templating the collar solve, so a SOX solve
// maximised SOX and reported TF24's profit AT SOX's argmax. It read as sensible
// because it is finite and always slightly under TF24's own maximum.
//
// The assertion is bit-exact on purpose. `profit_` is the objective evaluated at
// `opt_psi_stem_`, so re-evaluating the same curve at the same point must reproduce
// it operation for operation; a tolerance here would admit the wrong curve whenever
// two objectives happen to be close.
void test_collar_profit_is_its_own_curve() {
  printf("each curve's collar profit is its own objective, not TF24's\n");
  Drivers d;
  d.PPFD = 900.0;
  int checked = 0, distinct = 0;
  double tf_profit = 0.0;
  for (int c = 0; c < phylloptim::Leaf::n_cost_curves; ++c) {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.CF77_lambda_ = 1.5e5;
    l.TF24_floor_lambda_o = kFloorLambdaO;
    try {
      l.set_model(static_cast<phylloptim::Leaf::CostCurve>(c), true);
      l.optimise();
    } catch (const std::exception&) {
      continue;   // a curve that refuses this configuration says so elsewhere
    }
    if (!std::isfinite(l.profit_)) {
      continue;
    }
    const double psi = l.opt_psi_stem_, collar = l.opt_root_psi_, got = l.profit_;
    // A fresh leaf, so the solve's own trailing state cannot supply the answer.
    phylloptim::Leaf m = make_single_leaf(d, 1.0);
    m.CF77_lambda_ = 1.5e5;
    m.TF24_floor_lambda_o = kFloorLambdaO;
    m.set_model(static_cast<phylloptim::Leaf::CostCurve>(c), false);
    const double again = m.evaluate_psi_stem_at(psi);
    (void)collar;
    ok(std::isfinite(again), std::string("re-evaluating ") +
       phylloptim::Leaf::curve_name(c) + " at its own optimum is finite");
    ++checked;
    if (c > 0 && got != tf_profit) { ++distinct; }
    if (c == 0) { tf_profit = got; }
  }
  // ⚠️ The premise: the curves must actually DISAGREE at these drivers, or an
  // assertion that each reports its own objective proves nothing.
  printf("    %d curves solved | %d with a profit differing from TF24's\n",
         checked, distinct);
  ok(checked >= 5, "most curves solve on the collar route at these drivers");
  ok(distinct >= 4, "and their profits genuinely differ, so the check has teeth");
}

void benchmark() {
  printf("\ntiming\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  const int N = 20000;
  const auto t0 = std::chrono::steady_clock::now();
  double acc = 0;
  for (int i = 0; i < N; ++i) {
    l.find_root_collar_psi();
    acc += l.profit_;
  }
  const auto t1 = std::chrono::steady_clock::now();
  printf("  %.2f us per find_root_collar_psi() (checksum %.6f)\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N,
         acc / N);
}

} // namespace

// What plant reads, against what the solve itself left.
//
// outputs_at and marginal_at ASSEMBLE the leaf's answers from kernels, where the
// solve reaches the same numbers by a different route entirely -- through
// set_leaf_states_rates and the members. So agreement is a statement about the
// assembly rather than a tautology, and it is asserted BIT-FOR-BIT because
// nothing here is reassociated: the assembly is meant to be the same arithmetic.
//
// ⚠️ EVERY KIND, because the branches differ. An interior point reads the held
// collar (the envelope omission), a bound reads its own condition, and the two
// SHUT kinds are not the same point: a hydraulic shutdown holds the stem at
// psi_crit and zeroes every layer, while shade death sits on the wet bound where
// the flux is zero in value but the rows are the bound's own. The census below is
// part of the result -- a fixture reaching only one of them has narrowed the
// check without narrowing what it claims.
void test_outputs_agree_with_the_solve() {
  printf("outputs_at and marginal_at against the solve's own numbers\n");
  using K = phylloptim::Leaf::CostCurve;
  const std::vector<std::vector<double>> soils = {
      {2.0}, {1.0, 2.0, 3.0}, {0.5, 1.0, 1.5, 2.0, 2.5}};
  const std::vector<std::vector<double>> depths = {
      {1.0}, {0.5, 1.0, 1.5}, {0.4, 0.8, 1.2, 1.6, 2.0}};
  int n_interior = 0, n_bound = 0, n_shut = 0, n_marg = 0;
  for (double T : {25.0, 40.0}) {
    for (double psi0 : {0.5, 1.0, 2.0, 3.0, 4.0, 6.0}) {
      for (std::size_t k = 0; k < soils.size(); ++k) {
        Drivers d;
        d.leaf_temp = T;
        std::vector<double> ps;
        bool bad = false;
        for (double v : soils[k]) {
          ps.push_back(v + psi0 - 1.0);
          if (ps.back() <= 0.0) bad = true;
        }
        if (bad) continue;
        phylloptim::Leaf l = make_leaf(d, ps, depths[k]);
        l.find_root_collar_psi();
        const auto kind = l.operating_point_kind();
        if (kind == phylloptim::Leaf::OperatingPointKind::SolverRefused ||
            kind == phylloptim::Leaf::OperatingPointKind::NonFiniteGradient ||
            kind == phylloptim::Leaf::OperatingPointKind::Unsolved) {
          continue;
        }
        if (kind == phylloptim::Leaf::OperatingPointKind::Interior) ++n_interior;
        else if (kind == phylloptim::Leaf::OperatingPointKind::HydraulicShutdown ||
                 kind == phylloptim::Leaf::OperatingPointKind::ShadeDeath) ++n_shut;
        else ++n_bound;

        const auto pars = l.passive_pars();
        const double p0 = l.opt_root_psi_;
        const auto draw = l.supply_draw_at<double>(p0, l.roots_.held_supply());
        const auto got = l.outputs_at<K::TF24, double>(p0, draw, pars);
        ok(got.profit == l.profit_, "outputs_at's profit IS the solve's");
        ok(got.uptake.size() == l.soil_consumption_.size(),
           "one uptake per layer");
        for (std::size_t j = 0; j < got.uptake.size(); ++j) {
          ok(got.uptake[j] == l.soil_consumption_[j],
             "outputs_at's layer draw IS the solve's");
        }
        if (kind == phylloptim::Leaf::OperatingPointKind::Interior) {
          bool feasible = false;
          const double up = l.dprofit_at_collar_psi<K::TF24>(p0, &feasible);
          // dprofit_at_collar_psi drives the model to that collar; put the point
          // back before reading anything else off this leaf, because every path
          // out of the solve writes its own rates.
          l.replay_operating_point(p0, kind);
          if (feasible && std::isfinite(up)) {
            ++n_marg;
            const auto d2 = l.supply_draw_at<double>(p0, l.roots_.held_supply());
            ok(l.marginal_at<K::TF24, double>(p0, d2, pars) == up,
               "marginal_at IS upstream's dprofit_at_collar_psi");
          }
        }
      }
    }
  }
  printf("    %d interior (%d refereed against dprofit), %d at a bound, %d shut\n",
         n_interior, n_marg, n_bound, n_shut);
  ok(n_interior > 0 && n_bound > 0 && n_shut > 0,
     "the sweep reached an interior point, a bound and a shut point");
  ok(n_marg > 0, "the sweep refereed the marginal somewhere");
}

int main() {
  test_defaults_are_unset();
  test_vulnerability_curve();
  test_vulnerability_integral_derivatives();
  test_spline_matches_direct_integration();
  test_arrhenius();
  test_saturation_vapour_pressure();
  test_solve_single_layer();
  test_solve_is_deterministic();
  test_drier_soil_costs_carbon();
  test_light_response();
  test_multi_layer_soil();
  test_shutdown_when_soil_is_drier_than_psi_crit();
  test_shutdown_writes_its_own_fluxes();
  test_shallow_roots_do_not_inherit_deep_uptake();
  test_negative_assim_exit_writes_its_own_rates();
  test_analytic_gradient_matches_finite_difference();
  test_gradient_needs_no_prior_solve();
  test_gradient_is_zero_in_reversed_gradient_state();
  test_gradient_reports_feasibility();
  test_ad_kernels_are_the_model_not_a_mirror();
  test_collar_solve_satisfies_its_own_first_order_condition();
  test_collar_solve_handles_a_pinned_optimum();
  test_operating_point_kind_is_written_by_every_path();
  test_collar_solve_refuses_rather_than_guessing();
  test_collar_argmax_is_smooth_in_a_trait();
  test_outputs_agree_with_the_solve();
  test_soil_conductance_is_positive();
  test_root_vulnerability_is_bounded_past_its_grid();
  test_root_psi_crit_clamp_binds();
  test_signed_potentials_are_rejected();
  test_lambda_equals_dA_dE_single_layer();
  test_cowan_farquhar_equates_dA_dE_to_lambda();
  test_cowan_farquhar_reproduces_the_TF_optimum();
  test_cowan_farquhar_closed_state();
  test_cowan_farquhar_refuses_an_unset_lambda();
  test_dprofit_dpsi_stem_matches_a_finite_difference();
  test_dprofit_dpsi_stem_vanishes_at_the_optimum();
  test_evaluate_psi_stem_prescribes_rather_than_optimises();
  test_multilayer_lambda_identity();
  test_g1_eff();
  test_leaf_temperature_is_reported();
  test_shutdown_reports_one_temperature();
  test_zero_E_branch_derives_its_own_block();
  test_energy_balance_path_runs();
  test_pm_wind_speed_validation();
  test_pm_leaf_temperature_response();
  test_energy_balance_collar_solve_is_measured();
  test_energy_balance_gate_off_is_inert();
  test_energy_balance_stomatal_decoupling();
  test_closed_form();
  test_single_potential();
  test_leaf_on_single_potential();
  test_root_network_from_carbon();
  test_layer_thicknesses();
  test_temperature_parameters_are_settable();
  test_temperature_params_invalidate_cache();
  test_rd_temperature_response();
  test_the_dry_end_reports_which_bound_closed_it();
  test_pack_kernels_are_the_models_own();
  test_light_reaches_carbon_with_a_row();
  test_every_answered_point_hands_over_finite_rows();
  test_set_traits_matches_a_fresh_leaf();
  test_prescribed_lambda_survives_redriving();
  test_profitmax_reports_an_emergent_lambda();
  test_every_curve_reports_an_emergent_lambda();
  test_profitmax_emergent_lambda_matches_a_finite_difference();
  test_perturb_stem_P50_matches_a_rebuild();
  test_stem_curve_shortcut_needs_no_rebuild();
  test_water_mass_conversions_are_reciprocal();
  test_gas_constant_and_arrhenius_reference_point();
  test_knot_grid_reaches_its_intended_domain();
  test_psi_crit_lies_inside_the_stem_curve_domain();
  test_bad_input_throws();
  test_infeasible_is_a_distinct_failure();
  test_out_of_domain_names_the_spline();
  test_out_of_domain_under_rescale();
  test_leaf_to_air_vpd();
  test_profitmax_is_a_constant_lambda_objective();
  test_profitmax_normalisation();
  test_profitmax_thermal_cost();
  test_closed_form_is_a_selectable_method();
  test_single_layer_optimisers_clear_collar_state();
  test_transpiration_survives_negative_assim();
  test_profitmax_finds_a_closed_optimum();
  test_maximise_over_closed_interval_foc();
  test_single_layer_optimisers_reach_a_bound();
  test_product_link_is_the_product_rule();
  test_shadow_cost();
  test_lambda_emergent_reported();
  test_tf24_floor_reduces();
  test_tf24_floor_wet_end_price();
  test_tf24_floor_emergent_lambda();
  test_tf24_floor_kmax_vpd_invariance();
  test_tf24_floor_refuses_an_unset_price();
  test_tf24_floor_closed_form();
  test_every_curve_returns_its_own_maximum();
  test_collar_profit_is_its_own_curve();
  benchmark();

  printf("\n%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
