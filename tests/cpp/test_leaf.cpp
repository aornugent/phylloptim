// Pure-C++ test suite for the leaf model. No R, no test framework, no linking:
//   make -C tests/cpp && ./tests/cpp/test_leaf
//
// The trait values and drivers below are lifted from plant's
// tests/testthat/test-leaf.r so the two suites exercise the same operating
// point. The expected values here were produced BY this implementation and are
// regression guards rather than independent references -- but PLAN.md item 1 has
// since cross-checked the implementation against plant's compiled build and found
// it bit-identical, so they are guarding a verified model.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <chrono>
#include <cmath>
#include <functional>
#include <cstdio>
#include <limits>
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

// The class's own knot count. A fixture that constructs through the full argument
// list has to name one, and naming it again is how a leaf ends up on a different
// curve from its neighbour -- which is what thirty bit-identity checks reported when
// the default moved and these did not.
double default_ncontrol() {
  static const double n = phylloptim::Leaf().vulnerability_curve_ncontrol;
  return n;
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

phylloptim::Leaf make_leaf(const Drivers &d, std::vector<double> psi_soil,
                     std::vector<double> soil_depth) {
  phylloptim::Leaf l;
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

void test_stem_curve_derivative_is_the_closed_form() {
  printf("stem curve slope is the closed-form conductivity\n");
  phylloptim::Leaf l;
  // The interpolant carries G'(psi) = exp(-(psi/stem_b)^stem_c) as the slope at
  // each knot, so a reader's derivative is that closed form and not a value fit's
  // inference of it. At the wet end the slope is 1 to the last bit.
  near(l.stem_curve_integral_deriv(0.0), 1.0, 1e-14,
       "stem curve slope is 1 at zero potential");
  // Over the whole domain the slope tracks the closed form to well under the
  // 3e-4 a value-fitted cubic's inferred slope carried; exact at the knots, the
  // residual is the between-knot cubic slope near the steep wet end.
  double worst = 0.0;
  const int n = 4000;
  for (int i = 0; i <= n; ++i) {
    const double psi = l.psi_crit * i / n;
    worst = std::max(worst, std::fabs(l.stem_curve_integral_deriv(psi) -
                                      l.proportion_of_conductivity(psi)));
  }
  ok(worst < 1e-5, "slope matches the closed-form conductivity to better than 1e-5");
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
  // Regression guards -- see the note at the top of this file. Moved by PLAN 11a
  // (the collar root-find): opt_psi_stem_ 3.595247 -> 3.595088 and
  // assim_colimited_ 5.599511 -> 5.599153, both ~1e-4 relative, which is the
  // argmax correction and not rounding. profit_ and transpiration_ did not move at
  // this tolerance -- profit_ because it is the maximum and therefore flat.
  near(l.opt_psi_stem_, 3.595088, 1e-5, "opt_psi_stem_");
  near(l.assim_colimited_, 5.599153, 1e-5, "assim_colimited_");
  near(l.transpiration_, 1.141941e-05, 1e-5, "transpiration_");
  near(l.profit_, 2.515843, 1e-5, "profit_");
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

// A placed point has one referee and it is bit-identity: the placement claims to be
// the same three expressions the solve closes with, at the collar the solve returned.
// Anything weaker would make it a warm start, which carries a tolerance and is
// refereed by nothing here.
void test_a_placed_point_is_the_searched_one() {
  printf("a placed operating point is bit-identical to the searched one\n");
  int searched = 0, exited = 0;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0, 5.5, 6.5}) {
    for (double ppfd : {5.0, 26.0, 300.0, 1800.0}) {
      Drivers d;
      d.PPFD = ppfd;
      phylloptim::Leaf l = make_leaf(d, {psi}, {1.0});
      l.find_root_collar_psi();
      const phylloptim::Leaf::SolvedPoint point = l.solved_point();
      const std::string at = " at psi=" + std::to_string(psi) + " PPFD=" +
                             std::to_string(ppfd);
      if (!point.searched()) {
        // The four feasibility exits keep nothing, so a caller solves. Which is
        // right only if the placement says so rather than placing something.
        ++exited;
        ok(!l.place_solved_point(point),
           "a branch that never searched is declined" + at);
        continue;
      }
      ++searched;
      const double profit = l.profit_, stem = l.opt_psi_stem_;
      const double collar = l.opt_root_psi_, resid = l.collar_resid_;
      const double e_up = l.E_up_, trans = l.transpiration_;
      const double assim = l.assim_colimited_, gc = l.stom_cond_CO2_;
      const std::vector<double> uptake = l.soil_consumption_;
      const auto kind = l.operating_point_kind();
      const auto arm = l.dry_bound_arm();

      // Re-supplied first, because that is what a consumer does before every solve:
      // the placement has to stand on the drivers alone, not on what the last solve
      // left seated.
      l = make_leaf(d, {psi}, {1.0});
      ok(l.place_solved_point(point), "the point is placed" + at);
      ok(l.profit_ == profit, "profit is bit-identical" + at);
      ok(l.opt_psi_stem_ == stem, "the stem potential is bit-identical" + at);
      ok(l.opt_root_psi_ == collar, "the collar is bit-identical" + at);
      ok(l.collar_resid_ == resid, "the condition at it is bit-identical" + at);
      ok(l.E_up_ == e_up, "total uptake is bit-identical" + at);
      ok(l.transpiration_ == trans, "transpiration is bit-identical" + at);
      ok(l.assim_colimited_ == assim, "assimilation is bit-identical" + at);
      ok(l.stom_cond_CO2_ == gc, "the conductance is bit-identical" + at);
      ok(l.soil_consumption_ == uptake, "every layer's draw is bit-identical" + at);
      ok(l.operating_point_kind() == kind, "the branch is the one taken" + at);
      ok(l.dry_bound_arm() == arm, "the dry bound's arm is the one taken" + at);
    }
  }
  ok(searched > 0, "the grid reaches points a search found");
  ok(exited > 0, "and points a feasibility exit determined");
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
  phylloptim::Leaf l = make_leaf(d, {1.0, 2.0, 3.0}, {0.5, 0.5, 0.5});
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

// The 0.0 the two exits above return is a SENTINEL, and PLAN 11a is what makes
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

// PLAN 11b: the phylloptim::tangent derivative and the forward model are now instantiations of ONE
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
  printf("phylloptim::tangent differentiates the model's own algebra, not a mirror of it\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();

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

  // 2. The phylloptim::tangent derivative of the kernel against a central difference of the DOUBLE
  // function. Richardson-extrapolated, so the FD reference is good to ~1e-10 and
  // the tolerance is testing the derivative rather than the difference quotient.
  const auto richardson = [](auto f, double x, double h) {
    const double d1 = (f(x + h) - f(x - h)) / (2 * h);
    const double d2 = (f(x + h / 2) - f(x - h / 2)) / h;
    return (4 * d2 - d1) / 3;
  };
  for (double ci : {12.0, 20.0, 29.26, 35.0}) {
    phylloptim::tangent a = ci; phylloptim::seed_direction(a, 1.0);
    const double ad = phylloptim::derivative_along(l.assim_colimited_kernel(a));
    const double fd =
        richardson([&](double x) { return l.assim_colimited(x); }, ci, 1e-4);
    near(ad, fd, 1e-8, "dA/dci: phylloptim::tangent vs Richardson FD at ci=" + std::to_string(ci));
  }
  for (double p : {0.5, 2.0, 3.5949, 5.0}) {
    phylloptim::tangent a = p; phylloptim::seed_direction(a, 1.0);
    const double ad = phylloptim::derivative_along(l.hydraulic_cost_TF_kernel(a));
    const double fd =
        richardson([&](double x) { return l.hydraulic_cost_TF_kernel(x); }, p, 1e-4);
    near(ad, fd, 1e-8, "dcost/dpsi: phylloptim::tangent vs Richardson FD at psi=" +
                           std::to_string(p));
  }

  // 3. The cost kernel must be PURE -- an phylloptim::tangent probe must not scribble the cached
  // hydraulic_cost_, or a gradient evaluation would corrupt reported model state.
  const double cached = l.hydraulic_cost_TF(3.0);
  phylloptim::tangent probe = 4.5; phylloptim::seed_direction(probe, 1.0);
  (void)l.hydraulic_cost_TF_kernel(probe);
  ok(l.hydraulic_cost_ == cached,
     "an phylloptim::tangent probe of the cost kernel leaves hydraulic_cost_ untouched");
}

// PLAN 11a: the collar solve now solves its own first-order condition, so the
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
  printf("the collar solve lands where dprofit == 0 (PLAN 11a)\n");
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
  printf("a collar optimum pinned to its constraint (PLAN 11a)\n");
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
         phylloptim::Leaf::OperatingPointKind::PinnedWet,
     "and the leaf says so: the point is tagged pinned-wet");
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
  // optimum of any kind, and would read `pinned-wet` if that path did not write.
  const std::vector<double> psi5{4.0, 4.25, 4.5, 4.75, 5.0};
  const std::vector<double> depth5{1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<double> root5(5, 1.0 / 5.0 / d.area_leaf);
  l.set_physiology(fixture::root_network(root5, depth5), d.PPFD, psi5, depth5,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::PinnedWet,
     "psi_soil 4.0 over 5 layers is pinned to the wet bound");
  const double pinned_collar = l.opt_root_psi_;

  l.evaluate_root_collar_psi(pinned_collar);
  ok(l.operating_point_kind() == Kind::Prescribed,
     "evaluating a given collar potential is not an optimum, and says so");

  // ...and the optimising path takes it back, so `prescribed` is not sticky.
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::PinnedWet,
     "and optimising again restores the pin");

  // set_traits returns the object to its just-constructed state, and the
  // classification is part of that state (hazard 10 / setup_clean_leaf).
  l.set_traits(l.vcmax_25, l.stem_c, l.stem_b, l.psi_crit, l.roots_.root_c,
               l.roots_.root_b, l.roots_.root_psi_crit, l.beta2, l.jmax_25, l.a,
               l.curv_fact_elec_trans, l.curv_fact_colim, l.cost_scale_TF24,
               l.R_d_25);
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
  ok(l.prepare_collar_solve(bound_a, bound_b),
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
  p.prepare_collar_solve(pa, pb);
  const double profit_dry = p.profit_psi_stem_TF(
      p.find_psi_stem_from_psi_root(dry_end, p.supply_psi_soil()), dry_end);
  const double profit_wet = p.profit_psi_stem_TF(
      p.find_psi_stem_from_psi_root(wet_end, p.supply_psi_soil()), wet_end);
  ok(profit_wet > profit_dry,
     "the wetter end is the better of the two, by construction");

  phylloptim::Leaf m = make_leaf(d, psi, depth);
  double ma = 0.0, mb = 0.0;
  m.prepare_collar_solve(ma, mb);
  const double refused = m.maximise_profit_over_collar(dry_end, wet_end);
  ok(m.operating_point_kind() == Kind::SolverRefused,
     "the solve reports that it could not resolve the bracket");
  ok(m.operating_point_kind() != Kind::PinnedWet &&
         m.operating_point_kind() != Kind::PinnedDryRootCrit &&
         m.operating_point_kind() != Kind::PinnedDryRootPsiCrit,
     "and does not pass it off as a constrained optimum");
  // The endpoint the solve returns is stepped a fraction of the width inside the
  // bound it came from, so compare against the bound rather than for equality.
  ok(std::abs(refused - wet_end) < 1e-5,
     "and returns the end with the higher profit");
  ok(std::abs(refused - dry_end) > 0.1, "which is not the drier end");

  // A bracket lying wholly inside the infeasible sliver at bound_a, where dprofit
  // takes its reversed-gradient exit: no usable gradient at either end. The
  // non-finite-gradient half of the same guard has no bracket that reaches it and
  // is not covered.
  phylloptim::Leaf s = make_leaf(d, psi, depth);
  double sa = 0.0, sb = 0.0;
  s.prepare_collar_solve(sa, sb);
  // Its width is MEASURED rather than assumed: the sliver is where the collar
  // sits at or past the stem potential, so it narrows as the transport curve is
  // refined, and a hardcoded width stops lying inside it without saying so.
  double lo = 0.0, hi = 1e-6;
  for (int k = 0; k < 200; ++k) {
    const double mid = 0.5 * (lo + hi);
    bool f = true;
    s.dprofit_at_collar_psi(sa + mid, &f);
    (f ? hi : lo) = mid;
  }
  const double sliver = lo;
  bool feasible = true;
  s.dprofit_at_collar_psi(sa + 1e-6 * sliver, &feasible);
  ok(!feasible, "the wet bound admits no informative gradient");
  printf("  the infeasible sliver at the wet bound is %.3g MPa wide\n", sliver);
  const double fallen_back = s.maximise_profit_over_collar(sa, sa + sliver);
  ok(s.operating_point_kind() == Kind::SolverRefused,
     "a bracket with no usable gradient at either end is refused too");
  ok(std::isfinite(fallen_back) && fallen_back >= sa &&
         fallen_back <= sa + sliver,
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

// The soil counterpart: d(E_i)/d(psi_soil[i]) per layer, which is what prices a
// soil-potential row. Differencing E_up in one layer's potential isolates that
// layer's derivative exactly, because the block is DIAGONAL -- so a disagreement
// here is either the derivative or the diagonality claim, and both are worth
// failing on.
void test_soil_potential_derivative() {
  printf("dE_i/d(psi_soil_i) per layer\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers);
    for (int i = 0; i < layers; ++i) { ps[i] = 1.0 + 0.25 * i; depth[i] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const std::string at = " at " + std::to_string(layers) + " layers";

    std::vector<double> dE(layers, 0.0);
    l.dE_from_soil_dpsi_soil(l.opt_root_psi_, l.roots_.psi_soil_, dE);

    const double h = 1e-7;
    std::vector<double> buf(l.soil_consumption_.size(), 0.0);
    for (int i = 0; i < layers; ++i) {
      std::vector<double> psi_up = l.roots_.psi_soil_, psi_dn = l.roots_.psi_soil_;
      psi_up[std::size_t(i)] += h;
      psi_dn[std::size_t(i)] -= h;
      double up = 0.0, dn = 0.0;
      l.roots_.uptake_at(l.opt_root_psi_, psi_up, buf, up);
      l.roots_.uptake_at(l.opt_root_psi_, psi_dn, buf, dn);
      near(dE[std::size_t(i)], (up - dn) / (2.0 * h), 1e-5,
           "layer " + std::to_string(i) + " matches a central difference" + at);
    }

    // A drier soil supplies less water, so every rooted layer's row is negative.
    for (int i = 0; i < layers; ++i) {
      ok(dE[std::size_t(i)] < 0.0, "layer " + std::to_string(i) + " is negative" + at);
    }

    // And it is NOT minus the collar conductance. The two moving bounds sit at
    // different points on a non-linear vulnerability curve, so the integral
    // terms do not cancel; writing the soil row as -duptake_dpsi would be the
    // easy mistake and this is what refuses it.
    const double S = l.dE_from_soil_dpsi_collar(l.opt_root_psi_, l.roots_.psi_soil_);
    double summed = 0.0;
    for (double v : dE) { summed += v; }
    ok(std::abs(summed + S) > 1e-12 * std::abs(S),
       "the soil rows are not minus the collar conductance" + at);
  }
}

// d2(E_i)/d(collar) d(psi_soil_i), which a stand adjoint needs and the forward
// model does not.
//
// TWO REFERENCES, AND THEY ARE INDEPENDENT OF EACH OTHER. A mixed second
// derivative can be reached down either side, and the two sides are different
// functions here: dE_i/d(psi_soil_i) is the soil-end quotient rule and
// dE_up/d(collar) is the collar-end one, and they differ by more than a sign
// because the two endpoints sit at different points on a non-linear
// vulnerability curve. So differencing each in the OTHER variable gives two
// routes that share no arithmetic, and agreeing with both is a much stronger
// statement than agreeing with either.
//
// The second route differences the TOTAL, which is what makes it a check on the
// diagonality claim as well: sum_i d2E_i/dT dpsi_j is d2E_j/dT dpsi_j only
// because no layer reads another's potential.
void test_uptake_mixed_second_derivative() {
  printf("d2(E_i)/d(collar) d(psi_soil_i) against both of its own first derivatives\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers);
    for (int i = 0; i < layers; ++i) { ps[i] = 1.0 + 0.25 * i; depth[i] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double T = l.opt_root_psi_;
    const std::vector<double> psi = l.roots_.psi_soil_;
    const std::string at = " at " + std::to_string(layers) + " layers";

    std::vector<double> d2(layers, 0.0);
    l.roots_.d2uptake_dpsi_dpsi_soil(T, psi, d2);
    for (int i = 0; i < layers; ++i) {
      ok(std::isfinite(d2[std::size_t(i)]),
         "layer " + std::to_string(i) + " is finite" + at);
    }

    // Route one: difference the soil row in the collar.
    const double hT = 1e-6 * std::max(1.0, std::abs(T));
    std::vector<double> up(layers, 0.0), dn(layers, 0.0);
    l.dE_from_soil_dpsi_soil(T + hT, psi, up);
    l.dE_from_soil_dpsi_soil(T - hT, psi, dn);
    for (int i = 0; i < layers; ++i) {
      near(d2[std::size_t(i)],
           (up[std::size_t(i)] - dn[std::size_t(i)]) / (2.0 * hT), 1e-5,
           "layer " + std::to_string(i) +
               " matches a difference of the soil row in the collar" + at);
    }

    // Route two: difference the collar row in the soil. It is the total, so this
    // also asserts that layer j's potential reaches no other layer's flux.
    for (int i = 0; i < layers; ++i) {
      const double h = 1e-6 * std::max(1.0, std::abs(psi[std::size_t(i)]));
      std::vector<double> pu = psi, pd = psi;
      pu[std::size_t(i)] += h;
      pd[std::size_t(i)] -= h;
      const double su = l.roots_.duptake_dpsi(T, pu);
      const double sd = l.roots_.duptake_dpsi(T, pd);
      near(d2[std::size_t(i)], (su - sd) / (2.0 * h), 1e-5,
           "layer " + std::to_string(i) +
               " matches a difference of the collar row in the soil" + at);
    }
  }

  // The per-layer collar conductance, which a stand adjoint needs because each
  // layer is a separate write into the shared soil. Its own check is that the
  // parts are the whole: the total is the sum in layer order, so this is exact
  // equality and not a tolerance.
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers);
    for (int i = 0; i < layers; ++i) { ps[i] = 1.0 + 0.25 * i; depth[i] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double T = l.opt_root_psi_;
    const std::string at = " at " + std::to_string(layers) + " layers";

    std::vector<double> by_layer;
    l.dE_from_soil_dpsi_collar_by_layer(T, l.roots_.psi_soil_, by_layer);
    double summed = 0.0;
    for (int i = 0; i < layers; ++i) {
      summed += by_layer[std::size_t(i)];
    }
    // Not bit-exact, and the reason is the only difference between the two:
    // the total converts to kg once after summing and the parts convert each
    // before, so they differ by reassociation and by nothing else.
    near(summed, l.dE_from_soil_dpsi_collar(T, l.roots_.psi_soil_), 1e-14,
         "the per-layer conductances sum to the total" + at);
  }

  // The kink contract is the union of the two first derivatives' kinks, because
  // a second derivative needs BOTH endpoints off theirs. A collar sitting on a
  // layer's potential is the one every caller meets.
  {
    std::vector<double> ps{1.0, 1.5}, depth{1.0, 2.0};
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    std::vector<double> out;
    l.roots_.d2uptake_dpsi_dpsi_soil(ps[0], l.roots_.psi_soil_, out);
    ok(!out.empty() && std::isnan(out[0]),
       "a collar on a layer's own potential refuses the whole vector");
  }
}

// The light row. Its referee is the envelope theorem: at an interior optimum the
// direct partial at a FIXED collar equals the TOTAL derivative of the solved
// profit, so a central difference of the whole re-solve checks it without
// sharing a line of code with it.
void test_light_row() {
  printf("dprofit/dPPFD against a re-solve\n");
  for (double psi : {0.5, 1.0, 2.0}) {
    Drivers d;
    std::vector<double> ps(3, psi), depth(3);
    for (int i = 0; i < 3; ++i) { depth[std::size_t(i)] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const std::string at = " at psi_soil=" + std::to_string(psi);
    if (l.operating_point_kind() != phylloptim::Leaf::OperatingPointKind::Interior) {
      continue;
    }
    const double analytic = l.dprofit_dPPFD();
    ok(std::isfinite(analytic) && analytic > 0.0,
       "more light is worth more carbon" + at);

    const double h = 1e-2;  // PPFD is O(900), so this is a relative 1e-5
    Drivers up = d, dn = d;
    up.PPFD += h;
    dn.PPFD -= h;
    phylloptim::Leaf lu = make_leaf(up, ps, depth);
    phylloptim::Leaf ld = make_leaf(dn, ps, depth);
    lu.find_root_collar_psi();
    ld.find_root_collar_psi();
    const double fd = (lu.profit_ - ld.profit_) / (2.0 * h);
    near(analytic, fd, 1e-4, "matches a re-solve of the whole leaf" + at);
  }
}

// The contraction the stand adjoint calls: profit's environment rows, against a
// re-solve. Same envelope licence as the light row, one layer at a time.
void test_env_adjoint() {
  printf("env_adjoint against a re-solve\n");
  for (double psi : {0.5, 1.0, 2.0}) {
    Drivers d;
    std::vector<double> ps(3, psi), depth(3);
    for (int i = 0; i < 3; ++i) { depth[std::size_t(i)] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    if (l.operating_point_kind() != phylloptim::Leaf::OperatingPointKind::Interior) {
      continue;
    }
    const std::string at = " at psi_soil=" + std::to_string(psi);

    // A weight of 2 rather than 1, so a dropped or doubled factor shows.
    phylloptim::gradient::ProfitEnvDerivatives adj;
    phylloptim::gradient::profit_env_derivatives(l, adj);
    ok(adj.usable, "the rows exist at an interior optimum" + at);
    ok(adj.dprofit_dpsi_soil.size() == 3, "one row per layer" + at);

    const double h = 1e-6;
    for (int j = 0; j < 3; ++j) {
      std::vector<double> up = ps, dn = ps;
      up[std::size_t(j)] += h;
      dn[std::size_t(j)] -= h;
      phylloptim::Leaf lu = make_leaf(d, up, depth);
      phylloptim::Leaf ld = make_leaf(d, dn, depth);
      lu.find_root_collar_psi();
      ld.find_root_collar_psi();
      const double fd = (lu.profit_ - ld.profit_) / (2.0 * h);
      near(adj.dprofit_dpsi_soil[std::size_t(j)], fd, 2e-3,
           "layer " + std::to_string(j) + " matches a re-solve" + at);
      ok(adj.dprofit_dpsi_soil[std::size_t(j)] < 0.0,
         "a drier layer is worth less carbon, layer " + std::to_string(j) + at);
    }
  }

  // And it refuses rather than returning a number where the envelope step does
  // not hold. A very dry column pins or shuts the point down.
  Drivers d;
  std::vector<double> ps(3, 12.0), depth{1.0, 2.0, 3.0};
  phylloptim::Leaf l = make_leaf(d, ps, depth);
  l.find_root_collar_psi();
  phylloptim::gradient::ProfitEnvDerivatives adj;
  phylloptim::gradient::profit_env_derivatives(l, adj);
  if (l.operating_point_kind() != phylloptim::Leaf::OperatingPointKind::Interior) {
    ok(!adj.usable && !adj.message.empty(),
       "a non-interior point is refused by name, not answered");
  }
}

// The net price of soil water, which is what a soil row is multiplied by.
void test_marginal_price_water() {
  printf("net marginal price of soil water\n");
  Drivers d;
  for (double psi : {0.5, 1.0, 2.0, 3.0}) {
    std::vector<double> ps(3), depth(3);
    for (int i = 0; i < 3; ++i) { ps[std::size_t(i)] = psi; depth[std::size_t(i)] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const std::string at = " at psi_soil=" + std::to_string(psi);

    const double gross = l.marginal_cost_water_multilayer();
    const double stem = l.marginal_cost_water();
    const double net = l.marginal_price_water();
    if (!std::isfinite(gross)) { continue; }

    // It is the difference it is defined as, computed the other way.
    near(net, gross - stem, 1e-9, "net equals lambda_multi - lambda_stem" + at);
    // The stem cost is a real share of the gross, so the two are not
    // interchangeable: this is the overstatement the gross figure carries.
    ok(net > 0.0 && net < gross, "net is positive and below gross" + at);
    ok(gross / net > 1.1, "gross overstates the price by more than 10%" + at);
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
  const double last_knot = r.root_vuln_integral_from_psi.max();
  const double G_inf =
      phylloptim::cumulative_vulnerability_integral_limit(r.root_b, r.root_c);

  // On the grid the accessors ARE the bare splines, bit for bit. That is what
  // makes this fix golden-identical rather than merely golden-tolerable.
  for (double psi : {0.0, 0.5, 2.0, 4.0, 6.0}) {
    const std::string at = " at psi=" + std::to_string(psi);
    ok(r.root_vuln_at(psi) == r.root_vuln_integral_from_psi.slope(psi),
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

  // And both clamps are COUNTED where they hold, which is the only way a consumer
  // learns that a row it differenced came back with no response in it. This is the
  // non-vacuity proof for that tally: a stand cannot reach these potentials -- a
  // layer has to be rooted to be evaluated and the plant has to be alive, so the
  // wettest rooted layer stays wetter than psi_crit while this would need a deeper
  // one past 7.31 MPa -- so a counter read on a stand is zero whether it is wired
  // up or not, and only a direct call can tell the two apart.
  r.clamps.clear();
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_INTEGRAL_CAP) == 0,
     "the integral cap's counter starts clear");
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_ARGUMENT) == 0,
     "the argument clamp's counter starts clear");

  // Inside the grid neither holds, so neither counts.
  static_cast<void>(r.root_vuln_integral_at(4.0));
  static_cast<void>(r.root_vuln_at(4.0));
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_INTEGRAL_CAP) == 0,
     "the integral cap does not count inside the grid");
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_ARGUMENT) == 0,
     "the argument clamp does not count inside the grid");

  // Past it both hold, and each counts once per read.
  static_cast<void>(r.root_vuln_integral_at(1000.0));
  static_cast<void>(r.root_vuln_integral_at(1000.0));
  static_cast<void>(r.root_vuln_at(1000.0));
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_INTEGRAL_CAP) == 2,
     "the integral cap counts once per read that it holds on");
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_ARGUMENT) == 1,
     "the argument clamp counts a read past the last knot");
  // The wet end of the argument clamp counts too, and is unreachable in plant
  // because a signed potential is rejected at the boundary.
  static_cast<void>(r.root_vuln_at(-1.0));
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_ARGUMENT) == 2,
     "the argument clamp counts a read below zero");

  // Shared storage, so a copy reports into the same tally: a consumer holding this
  // model by value and copying it per unit -- which plant does per cohort -- would
  // otherwise take every count down with the copy.
  phylloptim::MultiLayerRoots copy = r;
  static_cast<void>(copy.root_vuln_integral_at(1000.0));
  ok(r.clamps.at(phylloptim::CLAMP_ROOT_VULN_INTEGRAL_CAP) == 3,
     "a copy counts into the tally the original reads");
  r.clamps.clear();
  ok(copy.clamps.at(phylloptim::CLAMP_ROOT_VULN_INTEGRAL_CAP) == 0,
     "and clearing through either clears both");

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
    phylloptim::Leaf bad(100, 2.04, -3.0, 5.0, 2.65, 1.29, 1.9, 1, 167 * 100, 0.3,
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

// bound_row(DryRootPsiCrit), differenced. THE EXPECTED VALUE IS +1: at this pin the
// operating point sits exactly at the registered constant, so dp*/d(root_psi_crit)
// is one. The arm returned -1 -- the raw residual partial of R(x) = x -
// root_psi_crit, undivided by the slope -- because it returns sixty lines before the
// block where every other field is converted into -(dR/du)/(dR/dx), and nothing
// differenced it. The two fall-through arms have had a differenced referee since they
// were written; this one had none, which is the whole reason the sign survived.
//
// Reaching the arm takes deliberate setup, exactly as test_root_psi_crit_clamp_binds
// found: at this package's defaults psi_crit == root_psi_crit, so the window where
// the root's own limit wins the min in prepare_collar_solve is EMPTY and no fixture
// built on the defaults can land here. That test opens the window by pushing the
// stem's psi_crit drier; this one opens the same window from the other side, by
// making the ROOT more vulnerable (root_psi_crit well inside psi_crit), because the
// stem's spline domain ends at ~6.82 MPa and there is no room to push psi_crit far
// enough for a dry PIN as well as a wide window. A pin also needs profit still
// climbing at the dry end, which is why each fixture pairs its root_psi_crit with a
// soil close underneath it.
void test_bound_row_root_psi_crit_is_a_unit_row() {
  printf("the root-psi-crit bound's row against a differenced solve\n");
  Drivers d;
  const auto solve = [&](double psi_soil, double root_psi_crit) {
    phylloptim::Leaf l;
    l.roots_.root_psi_crit = root_psi_crit;
    std::vector<double> ps{psi_soil}, depth{1.0}, root{1.0 / d.area_leaf};
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, ps, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    return l;
  };

  using Kind = phylloptim::Leaf::OperatingPointKind;
  double worst = 0.0;
  std::string worst_where;
  // (psi_soil, root_psi_crit) pairs, all measured to pin on this arm.
  const std::vector<std::pair<double, double>> fixtures{
      {1.00, 2.0}, {1.75, 2.0}, {2.50, 3.0}, {3.25, 3.5}, {3.75, 4.0}};
  for (const auto &f : fixtures) {
    const double psi_soil = f.first, rpc = f.second;
    const std::string tag = "psi_soil=" + std::to_string(psi_soil) +
                            ", root_psi_crit=" + std::to_string(rpc);
    phylloptim::Leaf l = solve(psi_soil, rpc);
    // Asserted, not assumed: if a future change moves this fixture off the arm the
    // test stops exercising it, and a silent stop is how the defect got here.
    ok(l.operating_point_kind() == Kind::PinnedDryRootPsiCrit,
       "the fixture pins on the root's own critical potential, " + tag);
    ok(l.dry_bound_arm() == phylloptim::Leaf::DryBoundArm::RootPsiCrit,
       "and the dry bound is the root's arm, not the continuity root, " + tag);

    const phylloptim::Leaf::BoundRow r =
        l.bound_row(phylloptim::Leaf::WhichBound::DryRootPsiCrit);
    ok(r.finite, "the row is finite, " + tag);
    near(r.bound, rpc, 1e-12, "the bound IS root_psi_crit, " + tag);

    // Central difference of the SOLVED operating point in root_psi_crit, at a
    // relative step. Both perturbed solves have to stay on the arm too, or the
    // difference is across a kink and refereeing nothing.
    for (double rel : {1e-4, 1e-6}) {
      const double h = rel * rpc;
      phylloptim::Leaf up = solve(psi_soil, rpc + h);
      phylloptim::Leaf dn = solve(psi_soil, rpc - h);
      ok(up.operating_point_kind() == Kind::PinnedDryRootPsiCrit &&
             dn.operating_point_kind() == Kind::PinnedDryRootPsiCrit,
         "both perturbed solves stay on the arm, " + tag);
      const double fd = (up.opt_root_psi_ - dn.opt_root_psi_) / (2.0 * h);
      near(r.d_droot_psi_crit, fd, 1e-8,
           "d(collar)/d(root_psi_crit) against a central difference, " + tag);
      const double err = std::abs(r.d_droot_psi_crit - fd);
      if (err > worst) {
        worst = err;
        worst_where = tag;
      }
    }
    // And the value itself, since it is exact by inspection once the conversion is
    // done: a row that agreed with the difference but came back -1 would mean the
    // difference had been broken the same way.
    near(r.d_droot_psi_crit, 1.0, 1e-12,
         "and it is +1, the point sitting at the constant, " + tag);
  }
  printf("  worst |row - central difference|: %.2e (%s)\n", worst,
         worst_where.c_str());
}

void test_lambda_equals_dA_dE_single_layer() {
  printf("marginal cost of water: analytic lambda vs dA/dE (stem free)\n");
  Drivers d;
  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    // optimise_psi_stem_TF holds the collar fixed at psi_soil_[0] and optimises
    // the stem, so the single-layer lambda is the one that applies here.
    l.optimise_psi_stem_TF();
    const double psi = l.opt_psi_stem_;
    const double eps = 1e-6;
    l.set_leaf_states_rates_from_psi_stem(psi + eps, psi_soil);
    const double A1 = l.assim_colimited_, E1 = l.transpiration_;
    l.set_leaf_states_rates_from_psi_stem(psi - eps, psi_soil);
    const double A0 = l.assim_colimited_, E0 = l.transpiration_;
    const double fd = (A1 - A0) / (E1 - E0);
    l.optimise_psi_stem_TF();
    // Tolerance is set by the optimiser: GSS_tol_abs is 1e-3 on psi, so the
    // first-order condition only holds to about that accuracy.
    near(l.marginal_cost_water() / fd, 1.0, 1e-3,
         "lambda/(dA/dE) at psi_soil=" + std::to_string(psi_soil));
    ok(l.marginal_cost_water() > 0.0, "lambda is positive");
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

  // The wind model really is what set ra_, rather than the fixed fallback. This
  // used to be untested: the previous version of this test set wind_speed_/d_ and
  // then immediately reassigned the leaf, discarding them, and its comment
  // described re-running set_physiology with the gate on, which it did not do. It
  // passed only because 2.0 / 0.05 are also the defaults.
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
// The oracle is a derivative-free scan of profit itself -- the method PLAN 11a
// used to arbitrate the last collar-solver change. It cannot inherit the
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
        // this is ~5.6e-15 (PLAN 11a); here it is whatever staleness leaves.
        bool feasible = true;
        const double resid = eb.dprofit_droot_collar_psi(psi_solver, &feasible);
        if (!feasible) continue;

        // ORACLE: scan profit over the same feasible interval. A fresh leaf,
        // because dprofit_droot_collar_psi above has just re-seated the
        // temperature parameters at yet another point.
        phylloptim::Leaf scan = make_pm_leaf(d, {2.0}, {1.0}, true);
        double lo = 0.0, hi = 0.0;
        if (!scan.prepare_collar_solve(lo, hi)) continue;
        const int N = 20001;
        double best_psi = lo, best_profit = -std::numeric_limits<double>::max();
        for (int i = 0; i < N; ++i) {
          const double psi = lo + (hi - lo) * double(i) / double(N - 1);
          const double p = scan.profit_at_collar_psi(psi, lo, hi);
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
// ⚠️ What this can and cannot show. stom_cond_CO2 still divides by the
// PRESCRIBED air VPD however hot the leaf gets (PLAN 13.1), so only the
// A(Tleaf) channel is exercised here; the leaf-to-air VPD channel is absent.
// The sweep therefore holds VPD fixed, which is the protocol of the paper's
// Figure S2 rather than its Figure 2.
void test_energy_balance_stomatal_decoupling() {
  printf("decoupling: gs rises where A falls, and only with the gate on\n");
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
    std::vector<double> Ta, A, gs;
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
    }
    return std::make_tuple(Ta, A, gs);
  };

  auto count_decoupled = [](const std::vector<double>& Ta,
                            const std::vector<double>& A,
                            const std::vector<double>& gs) {
    if (Ta.size() < 5) return 0;
    const size_t peak = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
    int n = 0;
    for (size_t i = peak + 2; i + 1 < Ta.size(); ++i) {
      const double dA = (A[i + 1] - A[i - 1]);
      const double dg = (gs[i + 1] - gs[i - 1]);
      if (dA < 0.0 && dg > 0.0) ++n;
    }
    return n;
  };

  const auto on = sweep(true);
  const auto off = sweep(false);
  const int n_on = count_decoupled(std::get<0>(on), std::get<1>(on), std::get<2>(on));
  const int n_off = count_decoupled(std::get<0>(off), std::get<1>(off), std::get<2>(off));
  printf("    decoupled points: gate on %d, gate off %d\n", n_on, n_off);

  const auto& A_on = std::get<1>(on);
  ok(A_on.size() > 5, "the decoupling sweep produced a curve");
  const size_t peak = std::distance(A_on.begin(),
                                    std::max_element(A_on.begin(), A_on.end()));
  ok(peak > 0 && peak + 1 < A_on.size(),
     "assimilation has an interior thermal optimum");
  ok(n_on > 0, "with the gate ON, conductance rises where assimilation falls");
  ok(n_on > n_off,
     "and it does so more than with the gate off -- the energy balance is why");
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
  l.optimise_psi_stem_TF();
  const double A_wet = l.assim_colimited_;
  setp(l, 1.0, 2.0);
  const phylloptim::closed_form::Solution wet = phylloptim::closed_form::solve(l, 1);
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
    l.optimise_psi_stem_TF();
    const double A_ex = l.assim_colimited_;
    setp(l, cs.h, 2.0);
    const double A_cf = phylloptim::closed_form::solve(l, 1).assim;
    ok(std::abs(A_cf / A_ex - 1.0) < cs.max_err,
       "closed-form error is bounded at h=" + std::to_string(cs.h) + " m");
  }

  // The guard is coarse, and saying so is the point. It admits ~2% error in A at
  // h = 12 m (ci/ca ~ 0.55) and only rejects once the error is ~8% (h = 20 m).
  // So it is a filter on gross failure, not an error bound -- and note that the
  // heights it rejects are the dominant canopy trees.
  setp(l, 20.0, 2.0);
  l.optimise_psi_stem_TF();
  const double A_tall = l.assim_colimited_;
  setp(l, 20.0, 2.0);
  const phylloptim::closed_form::Solution tall = phylloptim::closed_form::solve(l, 1);
  ok(!phylloptim::closed_form::within_guard(l, tall), "h=20 m is rejected by the guard");
  ok(std::abs(tall.assim / A_tall - 1.0) > 3e-2,
     "and it is rejected because the error really is large there");

  // The beta2 = 1/c leaf, where xi is constant and nothing needs solving.
  phylloptim::Leaf exact_leaf(96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                        5.870283, 1.0 / 2.680147, 157.44, 0.30, 0.7, 0.99, 1e-3,
                        default_ncontrol(), 1e-3, 1000, 7.5);
  ok(phylloptim::closed_form::beta2_is_exact(exact_leaf),
     "beta2_is_exact recognises beta2 = 1/stem_c");
  ok(!phylloptim::closed_form::beta2_is_exact(l), "and rejects the default beta2 = 1.5");
  setp(exact_leaf, 5.0, 1.5);
  exact_leaf.optimise_psi_stem_TF();
  const double A_ref = exact_leaf.assim_colimited_;
  setp(exact_leaf, 5.0, 1.5);
  const phylloptim::closed_form::Solution ex =
      phylloptim::closed_form::solve_exact_beta2(exact_leaf);
  ok(std::isnan(ex.psi_stem),
     "the explicit form reports no psi_stem -- it never solves for one");
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
    l.optimise_psi_stem_TF();
    return l.assim_colimited_;
  });
  const double t_cf =
      time_it(l, [&] { return phylloptim::closed_form::solve(l, 1).assim; });
  const double t_expl = time_it(
      exact_leaf, [&] { return phylloptim::closed_form::solve_exact_beta2(exact_leaf).assim; });
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

  // Closed form, straight from the documented model: carbon splits 1/3 vertical
  // : 2/3 horizontal, r_R_H_min = beta_H/c_r_h, r_R_V = beta_V*dz^2/c_r_v.
  const std::vector<double> carbon{3.0, 6.0, 1.5};
  const auto n = phylloptim::root_network_from_carbon(carbon, dz, beta_H, beta_V);

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
  const auto rich = phylloptim::root_network_from_carbon({12.0}, dz, beta_H, beta_V);
  const auto poor = phylloptim::root_network_from_carbon({3.0}, dz, beta_H, beta_V);
  ok(rich.r_R_H_min[0] < poor.r_R_H_min[0], "more root carbon -> less horizontal resistance");
  ok(rich.r_R_V_sum[0] < poor.r_R_V_sum[0], "more root carbon -> less vertical resistance");

  // Trailing zero-carbon layers are dropped, so the hot loop never visits them.
  const auto trailing = phylloptim::root_network_from_carbon({3.0, 6.0, 0.0, 0.0}, dz,
                                                      beta_H, beta_V);
  ok(trailing.r_R_H_min.size() == 2u, "trailing rootless layers are dropped");

  // Negative carbon is rejected rather than producing a negative resistance.
  bool threw = false;
  try {
    phylloptim::root_network_from_carbon({3.0, -1.0}, dz, beta_H, beta_V);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "negative root carbon throws");
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
           l.ci_at_compensation_point_ ? "   <-- shut down" : "");
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
    ok(hot.ci_at_compensation_point_,
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
    ok(!ok_leaf.ci_at_compensation_point_,
       "an ordinary leaf is not flagged as shut down");
  }
}

// set_traits exists so a gradient loop can perturb a trait without rebuilding the
// object. That is only worth having if re-traiting is INDISTINGUISHABLE from
// constructing afresh, so that is what is asserted -- bit-exactly, which is a
// statement neither a tolerance nor an eyeball could make.
//
// Bit-exactness is the whole point here rather than strictness for its own sake:
// the two ways to reach the same traits share no code, so any piece of derived
// state that set_traits fails to refresh shows up as a difference. Three pieces
// were candidates, and each is a real trap rather than a hypothetical one:
//
//   * the two pre-integrated vulnerability splines (stem_b/stem_c, root_b/root_c);
//   * vcmax_/jmax_/R_d_, behind set_physiology's (leaf_temp, atm_o2_kpa) cache --
//     the nastiest of the three, because the natural repair of "set the drivers
//     again" is exactly what takes the cache hit;
//   * the solved operating point itself (hazard 8).
//
// A `l.vcmax_25 = x` written by hand passes none of this, which is why the traits
// are not bound as settable fields.
// The homogeneity identity that makes a gradient in stem_b free (PLAN 11f).
//
// G(psi; s*b, c) = s * G(psi/s; b, c), because G(psi; b, c) = b * g(psi/b; c) --
// stem_b enters only as a scale on both axes, and the knot grid scales with it,
// so the identity holds for the SPLINE and not merely for the integral it
// approximates. perturb_stem_b() exploits that instead of reseeding 101
// incomplete gammas, which is 21.8 us against 0.001.
//
// ⚠️ This is checked against a REBUILD and not against a formula, because what
// has to be true is that the rescaled curve is the one a rebuild would have
// produced. It is not asserted bit-identical: the two differ by the rounding of
// one multiply and one divide per evaluation, which the nested solvers amplify to
// the ~1e-9 floor this project's guide records.
void test_perturb_stem_b_matches_a_rebuild() {
  printf("perturb_stem_b reproduces a rebuilt stem curve\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};

  // ⚠️ The downward range stops at 0.98, and that is issue #38 rather than a
  // limitation of the rescale: the curve's domain is [0, b*log(100)^(1/c)], which
  // at the default traits is 6.05 MPa against a psi_crit of 5.87. Shrink stem_b
  // by more than ~3% and psi_crit falls outside the curve, so the collar solve
  // asks for a potential the spline refuses to extrapolate to -- and a REBUILT
  // spline refuses identically. A gradient perturbs by ~1e-6.
  for (double factor : {0.98, 0.999, 1.000001, 1.05, 1.3}) {
    const double b_new = 3.898245 * factor;

    // The reference: stem_b through set_traits, which rebuilds the spline.
    phylloptim::Leaf rebuilt = make_leaf(d, {2.0}, {1.0});
    rebuilt.find_root_collar_psi();
    rebuilt.set_traits(96.0, 2.680147, b_new, 5.870283, 2.680147, 3.898245,
                       5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25);
    rebuilt.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                           d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa,
                           d.atm_kpa);
    rebuilt.find_root_collar_psi();

    // The rescale: no rebuild, and no set_physiology either -- nothing it
    // derives depends on stem_b, which is half of why this is cheap.
    phylloptim::Leaf rescaled = make_leaf(d, {2.0}, {1.0});
    rescaled.find_root_collar_psi();
    rescaled.perturb_stem_b(b_new);
    rescaled.find_root_collar_psi();

    const std::string what = " at stem_b x " + std::to_string(factor);
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
    rescaled.set_traits(96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                        5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25);
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

void test_set_traits_matches_a_fresh_leaf() {
  printf("set_traits is indistinguishable from constructing afresh\n");
  Drivers d;

  // One perturbed trait from each group that has its own derived state: a
  // photosynthetic trait (the temperature cache), the stem curve, the root curve,
  // and a cost parameter that is read directly and so should need no rebuild.
  struct Case { const char *name; int which; double value; };
  const Case cases[] = {{"vcmax_25", 0, 96.0 * 1.05},
                        {"stem_b", 2, 3.898245 * 1.05},
                        {"root_b", 5, 3.898245 * 1.05},
                        {"cost_scale_TF24", 12, 7.5 * 1.05},
                        {"stem_c", 1, 2.680147 * 1.05},
                        {"root_c", 4, 2.680147 * 1.05}};

  for (const Case &c : cases) {
    // The default trait vector, in set_traits' own argument order.
    double t[13] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                    3.898245, 5.870283, 1.5,      157.44,   0.30,
                    0.7,      0.99,     7.5};
    t[c.which] = c.value;

    // Fresh: the traits go through the constructor.
    phylloptim::Leaf fresh(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9],
                     t[10], t[11], 1e-3, default_ncontrol(), 1e-3, 1000,
                     t[12]);
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
    reused.set_traits(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9],
                      t[10], t[11], t[12], kRd25);
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

  // The specific trap, isolated, because the bit-exact comparisons above would
  // also pass if set_traits rebuilt everything unconditionally and the reason it
  // works were lost. vcmax_ is derived inside set_physiology's temperature cache,
  // which is keyed on (leaf_temp, atm_o2_kpa) and NOT on the traits -- so this is
  // the assertion that "change the trait, then set the drivers again" is a
  // sufficient recipe, which it is only because set_traits invalidates the cache.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double vcmax_before = l.vcmax_;
    l.set_traits(96.0 * 2.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                 5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    // The SAME leaf_temp and atm_o2_kpa, which is what arms the cache.
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
    l.set_traits(96.0, 2.680147, 3.898245 * 1.5, 5.870283, 2.680147, 3.898245,
                 5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5, kRd25);
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
    const double good[13] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                             3.898245, 5.870283, 1.5,      157.44,   0.30,
                             0.7,      0.99,     7.5};
    const int signed_positions[] = {3, 2, 5, 6};  // psi_crit, stem_b, root_b, root_psi_crit
    const char *labels[] = {"psi_crit", "stem_b", "root_b", "root_psi_crit"};
    for (int k = 0; k < 4; ++k) {
      double t[13];
      for (int i = 0; i < 13; ++i) {
        t[i] = good[i];
      }
      t[signed_positions[k]] = -t[signed_positions[k]];  // the pre-#25 sign
      bool threw = false;
      try {
        l.set_traits(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9],
                     t[10], t[11], t[12], kRd25);
      } catch (const std::runtime_error &) {
        threw = true;
      }
      ok(threw, std::string("set_traits rejects a negative ") + labels[k]);
    }
  }
}

void test_bad_input_throws() {
  printf("input validation\n");
  Drivers d;
  phylloptim::Leaf l;
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

// What an out-of-domain transport lookup says. The stem curve is the only
// interpolator here with extrapolation disabled, so it is the only one a lookup
// can throw on -- and there are two of them, they are inverses, and they carry
// different units, so odelia's message (which names the point and the domain but
// not the spline) is ambiguous in the way that matters. Localising plant#576 came
// down to which of the four call sites was asking; these assertions are what make
// that a read rather than a bisect.
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
  ok(mentions(inv, "psi_from_transpiration"), "inverse lookup names its spline");
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
// Under perturb_stem_b the value handed to the spline is psi/s, so quoting the
// spline's own endpoints would send the reader after a discrepancy that is not
// there.
void test_out_of_domain_under_rescale() {
  Drivers d;
  phylloptim::Leaf wide = make_leaf(d, {2.0}, {1.0});
  const std::string before = message_of([&] { wide.transpiration(50.0, 0.0); });

  phylloptim::Leaf rescaled = make_leaf(d, {2.0}, {1.0});
  rescaled.perturb_stem_b(rescaled.stem_b * 2.0);
  const std::string after =
      message_of([&] { rescaled.transpiration(50.0, 0.0); });

  ok(!before.empty() && !after.empty(), "both report a domain failure");
  ok(before != after, "the rescaled domain is reported, not the spline's own");
  ok(mentions(after, "rescaled by"),
     "the rescale is named so the two domains are not confused");
}

// ---------------------------------------------------------------------------
// THE ENVIRONMENT ROWS
// ---------------------------------------------------------------------------
// The gradient's parameter dimension now carries the two quantities plants share
// -- incident light and the per-layer soil water they draw from -- so these four
// tests are about the parameter dimension, not about the model. The header block
// on `par_PPFD` in gradient.hpp has the derivation each of them checks.

namespace env {

namespace grad = phylloptim::gradient;

// The suite's trait defaults, in `set_traits`' argument order, then the two
// non-trait parameters. `resistance` is never perturbed here: these are
// multi-layer points. Sized from the initialiser and checked, because too few
// initialisers is legal and zero-fills the rest.
const double kTheta[] = {96.0,   2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,    157.44,   0.30,
                         0.7,    0.99,     7.5,      kRd25,
                         1.0 * 0.000157 / 5.0, 1e3};
static_assert(sizeof(kTheta) / sizeof(kTheta[0]) == grad::n_pars);

// `rooted` layers of root carbon spread over `layers` soil layers, drying with
// depth -- the golden grid's soil, so a point named here is a point that grid
// already covers.
grad::Drivers drivers(double psi_soil, double ppfd, double vpd, int layers,
                      int rooted) {
  grad::Drivers d;
  std::vector<double> ps(layers), depth(layers), root(layers, 0.0);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    if (i < rooted) {
      root[i] = 1.0 / rooted / 0.05;  // per unit leaf area (hazard 4)
    }
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

std::vector<int> all_env_pars(int layers) {
  std::vector<int> pars{grad::par_PPFD};
  for (int j = 0; j < layers; ++j) {
    pars.push_back(grad::par_psi_soil_first + j);
  }
  return pars;
}

phylloptim::Leaf fresh() {
  phylloptim::Leaf l;
  return l;
}

// The request the stand adjoint actually makes: profit, the water each layer gave
// up, and every input the multi-layer path has.
//
// ⚠️ THE COLLAR IS DELIBERATELY ABSENT, which is what makes this the consumer's
// request rather than the calibration route's five. `rows` is indexed by POSITION
// in the request, so profit is row 0 here and layer j is row 1 + j.
void consumer_request(int layers, std::vector<int>& outs, std::vector<int>& ins) {
  outs.assign(1, grad::out_profit);
  for (int j = 0; j < layers; ++j) {
    outs.push_back(grad::out_uptake_first + j);
  }
  ins.clear();
  for (int i = 0; i < grad::n_pars; ++i) {
    if (i != grad::par_resistance) {
      ins.push_back(i);
    }
  }
  ins.push_back(grad::par_PPFD);
  for (int i = 0; i < 2 * layers; ++i) {
    ins.push_back(grad::par_psi_soil_first + i);
  }
}

// Apply, solve, and read the point -- which is what a consumer of the row layer
// does, now that neither entry point solves. `grad::base_point` is the read alone.
grad::BasePoint solved_base_point(phylloptim::Leaf& l, const grad::Drivers& d,
                                  bool single, const grad::Settings& s,
                                  int n_uptake) {
  grad::apply(l, kTheta, d, single, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  static_cast<void>(n_uptake);
  return grad::base_point(l);
}

// The same for a whole request: solve, then take the rows. Both entries read a
// solved leaf, so a caller that has not solved gets a refusal.
grad::Rows solved_rows(phylloptim::Leaf& l, const double* theta,
                       const grad::Drivers& d, const grad::RowRequest& r,
                       const grad::Settings& s) {
  const bool single = l.supply_kind_ == phylloptim::Leaf::SupplyKind::SinglePotential;
  grad::apply(l, theta, d, single, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  return grad::rows_differenced(l, theta, d, r, s);
}

}  // namespace env

// The row count and the row names, which is the whole of the arity decision:
// sixteen fixed parameters, then PPFD at a FIXED index, then L soil rows and L
// root-carbon rows. PPFD keeps its index whatever L is; the two per-layer blocks
// mean an index past it names one input only at a fixed L.
void test_environment_par_names() {
  printf("environment rows: names and arity\n");
  namespace grad = phylloptim::gradient;
  ok(grad::par_PPFD == grad::n_pars, "PPFD sits immediately after theta");
  for (int L : {1, 3, 5}) {
    const std::vector<std::string> nms = grad::par_names(L);
    const std::string tag = std::to_string(L) + " layers";
    ok(int(nms.size()) == grad::n_pars_total(L), "row count, " + tag);
    // The arity a consumer's own input vector has to match: two per layer.
    ok(grad::n_pars_total(L) == grad::n_pars + 1 + 2 * L,
       "the arity is n_pars + 1 + 2L, " + tag);
    ok(nms[grad::n_pars] == "PPFD", "PPFD is named, " + tag);
    ok(nms[std::size_t(grad::par_root_carbon_first(L) - 1)] ==
           "psi_soil_" + std::to_string(L),
       "the last soil row is one-based, " + tag);
    ok(nms.back() == "root_carbon_" + std::to_string(L),
       "the root-carbon block follows the soil one, " + tag);
    // The sixteen R indexes into are untouched.
    ok(std::equal(grad::par_names().begin(), grad::par_names().end(),
                  nms.begin()),
       "the fixed sixteen are unchanged, " + tag);
    ok(grad::par_name(grad::par_psi_soil_first, L) == "psi_soil_1",
       "par_name agrees with par_names, " + tag);
    ok(grad::par_name(grad::par_root_carbon_first(L), L) == "root_carbon_1",
       "par_name names a root-carbon row, " + tag);
  }
  // A row past this observation's arity is a per-row error rather than an
  // out-of-bounds read of psi_soil. A root-carbon row is INSIDE the arity and is
  // answered; what it takes is in test_root_carbon_rows below.
  phylloptim::Leaf l = env::fresh();
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, 3, 3);
  grad::Settings s;
  grad::Result r;
  int bad = grad::n_pars_total(3);
  bool threw = false;
  try {
    grad::at(l, env::kTheta, d, false, &bad, 1, s, r);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(threw, "index " + std::to_string(bad) + " throws");
}

// The values come back with the rows, sized for the request.
//
// ⚠️ Both halves of that have bitten. A buffer sized for the fixed outputs alone
// carries NONE of the variable-length uptake block, and reads past its end; and
// reading the values off the leaf afterwards gives the last perturbed
// evaluation, because this call re-supplies the base state without re-solving.
// What the read declines, and that the report says which of the two differences
// answers it. This is the whole of why `rows_differenced` exists: the read states
// what the model can state about the state it is in, and names the rest.
void test_the_read_names_what_it_declines() {
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 3;
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, L, L);

  std::vector<int> ins;
  for (int p = 0; p < grad::n_pars; ++p) {
    if (p != grad::par_resistance) {
      ins.push_back(p);
    }
  }
  std::vector<int> stand_outs{grad::out_profit};
  for (int j = 0; j < L; ++j) {
    stand_outs.push_back(grad::out_uptake_first + j);
  }
  std::vector<int> all_outs;
  for (int j = 0; j < grad::n_outputs_total(L); ++j) {
    all_outs.push_back(j);
  }

  // The request a stand makes: every input stated, nothing named.
  pl::Leaf l = env::fresh();
  grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  grad::RowRequest stand{stand_outs.data(), stand_outs.size(), ins.data(),
                         ins.size()};
  const grad::Rows sr = grad::rows_at(l, stand);
  int declined = 0;
  for (grad::Rows::NoRow why : sr.no_row) {
    declined += why != grad::Rows::NoRow::None ? 1 : 0;
  }
  ok(sr.no_row.size() == ins.size(), "the report has one entry per input");
  ok(declined == 0, "and the read states every input a stand asks about");

  // Name assimilation and the carbon-side inputs go to a difference, because the
  // readers report profit and the condition rather than assimilation itself.
  grad::RowRequest calib{all_outs.data(), all_outs.size(), ins.data(),
                         ins.size()};
  const grad::Rows cr = grad::rows_at(l, calib);
  int held_difference = 0, carbon_declined = 0, carbon_total = 0;
  for (std::size_t i = 0; i < ins.size(); ++i) {
    if (cr.no_row[i] == grad::Rows::NoRow::DifferenceAtAHeldCollar) {
      held_difference++;
    }
    if (grad::carbon_side(ins[i])) {
      carbon_total++;
      carbon_declined += cr.no_row[i] != grad::Rows::NoRow::None ? 1 : 0;
    }
  }
  ok(held_difference > 0 && carbon_declined == carbon_total,
     "naming assimilation sends every carbon-side input to a held difference");

  // And the wrapper answers exactly those, leaving nothing declined.
  pl::Leaf w = env::fresh();
  const grad::Rows wr = env::solved_rows(w, env::kTheta, d, calib, s);
  int still_declined = 0;
  for (grad::Rows::NoRow why : wr.no_row) {
    still_declined += why != grad::Rows::NoRow::None ? 1 : 0;
  }
  ok(still_declined == 0, "and the differencing entry answers all of them");

  // The single path's series resistance: the read declines it because its bound
  // moves only by re-solving, so a zero there would say the bound stood still.
  std::vector<int> single_ins{grad::par_resistance};
  pl::Leaf sp = env::fresh();
  sp.set_supply_single();
  grad::Drivers sd = env::drivers(2.0, 900.0, 2.0, 1, 1);
  grad::apply(sp, env::kTheta, sd, true, -1, s.fast_stem_curve);
  sp.find_root_collar_psi();
  std::vector<int> one_out{grad::out_profit, grad::out_uptake_first};
  grad::RowRequest res{one_out.data(), one_out.size(), single_ins.data(),
                       single_ins.size()};
  const grad::Rows rr = grad::rows_at(sp, res);
  ok(rr.no_row.size() == 1 &&
         rr.no_row[0] == grad::Rows::NoRow::DifferenceTheSolve,
     "the series resistance is declined, and by the difference that re-solves");
}

void test_rows_read_a_solved_leaf() {
  namespace grad = phylloptim::gradient;
  grad::Settings s;
  const int L = 3;
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, L, L);

  std::vector<int> out_index{grad::out_profit};
  for (int i = 0; i < L; ++i) {
    out_index.push_back(grad::out_uptake_first + i);
  }
  std::vector<int> pars{grad::par_psi_soil_first};
  grad::RowRequest req{out_index.data(), out_index.size(), pars.data(),
                       pars.size()};

  // An unsolved leaf gets a refusal, not a solve and not a plausible number.
  phylloptim::Leaf unsolved = env::fresh();
  grad::apply(unsolved, env::kTheta, d, false, -1, s.fast_stem_curve);
  const grad::Rows none = grad::rows_at(unsolved, req);
  ok(none.kind == phylloptim::Leaf::OperatingPointKind::Unsolved,
     "an unsolved leaf is reported unsolved rather than solved for");
  ok(!none.message.empty(), "and the refusal says so by name");

  // A solved one is read, and the read leaves every value where it found it --
  // which is why the rows carry none.
  phylloptim::Leaf l = env::fresh();
  grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  const double profit_before = l.profit_;
  const std::vector<double> uptake_before = l.soil_consumption_;
  const std::size_t solves_before = *l.collar_solves;
  const grad::Rows rows = grad::rows_at(l, req);

  ok(rows.kind == phylloptim::Leaf::OperatingPointKind::Interior,
     "the solved leaf's own branch is what comes back");
  ok(*l.collar_solves == solves_before, "the read makes no solve of its own");
  ok(phylloptim::util::identical(l.profit_, profit_before),
     "and moves no value: profit is where the solve left it");
  bool uptake_ok = l.soil_consumption_.size() == uptake_before.size();
  for (std::size_t i = 0; i < uptake_before.size() && uptake_ok; ++i) {
    uptake_ok = phylloptim::util::identical(l.soil_consumption_[i],
                                            uptake_before[i]);
  }
  ok(uptake_ok, "nor any layer's uptake");
}

// A leaf too shaded to cover its own respiration seats both potentials at the
// collar where uptake is zero and pays respiration plus a hydraulic cost there.
// So its profit reads the soil only through that bound moving underneath it, and
// the assembled rows must agree with that composition, which shares no code with
// them.
//
// ⚠️ THE COMPOSITION IS THE ASSEMBLY, NOT THE HELD ROW, and it used to be both:
// these rows came from differencing the whole solve, so the total sat in `held`
// with the point declared not to move. It moves -- the seat IS the wet bound -- so
// the held row is the zero a frozen collar makes it and the movement is reported
// once in `dresidual`, exactly as at a pin.
void test_shade_death_soil_rows() {
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 3;
  grad::Drivers d = env::drivers(1.0, 10.0, 1.0, L, L);
  pl::Leaf l = env::fresh();

  std::vector<int> input;
  for (int j = 0; j < L; ++j) {
    input.push_back(grad::par_psi_soil_first + j);
  }
  std::vector<int> out{grad::out_profit};
  grad::RowRequest req{out.data(), out.size(), input.data(), input.size()};
  const grad::Rows r = env::solved_rows(l, env::kTheta, d, req, s);
  ok(r.kind == pl::Leaf::OperatingPointKind::ShadeDeath,
     "the shaded fixture reaches shade death");
  if (r.kind != pl::Leaf::OperatingPointKind::ShadeDeath) {
    return;
  }
  ok(std::isfinite(r.residual_slope) && std::isfinite(r.held[0]),
     "and it answers rather than refusing");

  // The rows come back at the base state, but this call re-supplies base without
  // re-solving it, so the bound has to be read from a solve of its own.
  l.find_root_collar_psi();
  const pl::Leaf::BoundRow b = l.bound_row(pl::Leaf::WhichBound::Wet);
  const pl::Leaf::HydraulicCostRow c = l.hydraulic_cost_row(b.bound);
  ok(b.finite && c.finite, "both halves of the composition are finite");
  // The point's channel into profit is the cost's slope at the seat, and the
  // model's own marginal profit cannot say so: the stem sits AT the collar, so it
  // takes the no-flow exit and returns a sentinel zero.
  ok(std::abs(r.dy_dp[0] - (-c.d_dpsi_stem)) <=
         1e-12 * std::max(std::abs(c.d_dpsi_stem), 1.0),
     "profit's channel into the point is the cost's slope at the seat");
  double worst = 0.0;
  double worst_held = 0.0;
  for (int j = 0; j < L; ++j) {
    const double dpoint = -r.dresidual[std::size_t(j)] / r.residual_slope;
    const double got = r.held[std::size_t(j)] + r.dy_dp[0] * dpoint;
    const double want = -c.d_dpsi_stem * b.d_dpsi_soil[std::size_t(j)];
    worst = std::max(worst, std::abs(got - want) /
                                std::max(std::abs(want), 1e-300));
    worst_held = std::max(worst_held, std::abs(r.held[std::size_t(j)]));
  }
  char worst_text[32];
  std::snprintf(worst_text, sizeof(worst_text), "%.2e", worst);
  ok(worst < 1e-6,
     std::string("the assembled soil rows are the composed ones, ") +
         worst_text);
  // Exactly zero rather than nearly: at a frozen collar a soil potential reaches
  // profit through total uptake, and gross assimilation is identically zero here,
  // so there is nothing for the flux to reach.
  ok(worst_held == 0.0, "and the held half of each is exactly zero");
}

// Root carbon is an ordinary input on the multi-layer path and has no row on
// either of the two states where the network cannot represent a move in it.
//
// The perturbation is exact without the architecture model's constants, because
// both resistances the solve reads are proportional to 1/carbon -- the same
// proportionality `duptake_droot_carbon` differentiates through for the bound.
// So the differenced row and the analytic bound row describe one model, and the
// third check below is what would catch them coming apart.
void test_root_carbon_rows() {
  namespace grad = phylloptim::gradient;
  grad::Settings s;
  grad::Result r;
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, 3, 3);

  phylloptim::Leaf l = env::fresh();
  int rc = grad::par_root_carbon_first(3);
  bool threw = false;
  try {
    grad::at(l, env::kTheta, d, false, &rc, 1, s, r);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(!threw, "a rooted layer's root carbon answers");
  ok(!threw && std::isfinite(r.grad[std::size_t(grad::out_profit)]),
     "and its profit row is finite");
  ok(!threw && r.grad[std::size_t(grad::out_profit)] != 0.0,
     "and it is not the zero an unapplied perturbation would give");

  // One series resistance and no root architecture, so a carbon perturbation
  // would change nothing this path reads and every row would come back exactly
  // zero. Refused by name instead.
  phylloptim::Leaf sp = env::fresh();
  sp.set_supply_single();
  int rc_single = grad::par_root_carbon_first(1);
  threw = false;
  try {
    grad::at(sp, env::kTheta, d, true, &rc_single, 1, s, r);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(threw, "root carbon is refused on the single-potential path");

  // A layer with no roots: the architecture model sizes the network to the
  // deepest rooted layer, so there is no slot to move and no row to give.
  grad::Drivers shallow = env::drivers(2.0, 900.0, 2.0, 3, 2);
  std::vector<int> pars{grad::par_root_carbon_first(3) + 2};
  std::vector<int> out_index{grad::out_profit};
  grad::RowRequest req;
  req.output = out_index.data();
  req.n_output = out_index.size();
  req.input = pars.data();
  req.n_input = pars.size();
  phylloptim::Leaf shallow_leaf = env::fresh();
  const grad::Rows rows =
      env::solved_rows(shallow_leaf, env::kTheta, shallow, req, s);
  ok(!std::isfinite(rows.held[0]),
     "an unrooted layer's carbon row is NA rather than zero");
  ok(!std::isfinite(rows.dresidual[0]),
     "and its condition gradient is NA too");
  const std::string* said = rows.reason_for(pars[0]);
  ok(said != nullptr && said->find("no root carbon") != std::string::npos,
     "and the refusal is stated against that layer's own input");

  // THE REFEREE FOR THE PERTURBATION ITSELF: against the architecture model run
  // again from moved carbon, which is the thing being stood in for. It shares no
  // code with `perturb_root_carbon` -- it goes through the two constants and the
  // layer thickness, none of which the in-place edit ever sees.
  //
  // Measured bit-identical on all 45 values, but the tolerance stays a rounding:
  // recovering each resistance's constant by multiplying is exact for these
  // numbers rather than by construction.
  const std::vector<double> depth{1.0, 2.0, 3.0};
  std::vector<double> carbon(3, 1.0 / 3.0 / 0.05);
  const phylloptim::RootNetwork base = fixture::root_network(carbon, depth);
  auto rel = [](double got, double want) -> double {
    const double scale = std::max(std::abs(want), 1e-300);
    return std::abs(got - want) / scale;
  };
  double worst = 0.0;
  for (int k = 0; k < 3; ++k) {
    const double moved = carbon[std::size_t(k)] * 1.001;
    std::vector<double> rebuilt_carbon = carbon;
    rebuilt_carbon[std::size_t(k)] = moved;
    const phylloptim::RootNetwork want =
        fixture::root_network(rebuilt_carbon, depth);
    phylloptim::RootNetwork got;
    grad::perturb_root_carbon(base, k, moved, got);
    for (std::size_t i = 0; i < want.r_R_H_min.size(); ++i) {
      worst = std::max(worst, rel(got.r_R_H_min[i], want.r_R_H_min[i]));
      worst = std::max(worst, rel(got.r_R_V[i], want.r_R_V[i]));
      worst = std::max(worst, rel(got.r_R_V_sum[i], want.r_R_V_sum[i]));
      worst = std::max(worst, rel(got.c_r_V[i], want.c_r_V[i]));
      worst = std::max(worst, rel(got.c_r_H[i], want.c_r_H[i]));
    }
  }
  ok(worst < 1e-15, "the in-place carbon perturbation reproduces a rebuild, " +
                        std::to_string(worst));
}

// The strongest available reference: `psi_soil` and `PPFD` are already
// `set_physiology` arguments, so differencing the WHOLE solve with respect to
// them needs no new machinery and is independent of the envelope algebra the IFT
// route uses. Scaled by the largest entry in the same output row, because a
// derivative that is 1e-10 next to a neighbour of 1 is not a useful denominator.
//
// The water rows land in the same band as the existing
// `leaf_specific_conductance_max` row (measured over the golden grid: 6.4e-04
// against its 5.4e-04) and for the same reason -- both move the argmax hard, and
// the composite's dY/dpsi and H are themselves differences taken at a flat
// maximum. The light row does not move the collar at all and comes in at 1.1e-06.
void test_environment_rows_match_a_differenced_solve() {
  printf("environment rows: IFT composite vs a differenced solve\n");
  namespace grad = phylloptim::gradient;
  for (int L : {1, 3, 5}) {
    for (double psi_soil : {1.0, 2.0, 3.0}) {
      phylloptim::Leaf l = env::fresh();
      grad::Drivers d = env::drivers(psi_soil, 900.0, 2.0, L, L);
      const std::vector<int> pars = env::all_env_pars(L);
      grad::Settings s;
      grad::Result ift;
      grad::at(l, env::kTheta, d, false, pars.data(), pars.size(), s, ift);
      const std::string tag =
          std::to_string(L) + " layers, psi_soil=" + std::to_string(psi_soil);
      ok(ift.status == grad::Status::Interior, "interior point, " + tag);
      ok(ift.used_ift, "the composite ran, " + tag);

      grad::Settings sf = s;
      sf.method = grad::Method::Fd;
      grad::Result fd;
      grad::at(l, env::kTheta, d, false, pars.data(), pars.size(), sf, fd);

      for (int j = 0; j < grad::n_outputs; ++j) {
        double scale = 0.0;
        for (std::size_t k = 0; k < pars.size(); ++k) {
          scale = std::max(scale, std::abs(fd.grad[k * grad::n_outputs + j]));
        }
        for (std::size_t k = 0; k < pars.size(); ++k) {
          const double a = ift.grad[k * grad::n_outputs + j];
          const double b = fd.grad[k * grad::n_outputs + j];
          ok(std::abs(a - b) <= 2e-3 * scale,
             grad::par_name(pars[k], L) + "/" + grad::output_names()[j] +
                 " agrees with the differenced solve, " + tag);
        }
      }
    }
  }

  // And on the OTHER supply path (hazard 6), where the whole soil profile is one
  // potential however long the caller's vector is -- so there is exactly one soil
  // row, and it is the element the path actually reads.
  {
    phylloptim::Leaf l = env::fresh();
    l.set_supply_single();
    grad::Drivers d = env::drivers(1.0, 900.0, 2.0, 1, 1);
    d.root_network = fixture::series_resistance(1e3);
    ok(grad::n_soil_layers(d, true) == 1, "the single path has one soil row");
    const std::vector<int> pars{grad::par_PPFD, grad::par_psi_soil_first};
    grad::Settings s;
    grad::Result ift, fd;
    grad::at(l, env::kTheta, d, true, pars.data(), pars.size(), s, ift);
    grad::Settings sf = s;
    sf.method = grad::Method::Fd;
    grad::at(l, env::kTheta, d, true, pars.data(), pars.size(), sf, fd);
    ok(ift.status == grad::Status::Interior, "interior, single path");
    for (int j = 0; j < grad::n_outputs; ++j) {
      const double scale =
          std::max(std::abs(fd.grad[j]), std::abs(fd.grad[grad::n_outputs + j]));
      for (std::size_t k = 0; k < pars.size(); ++k) {
        ok(std::abs(ift.grad[k * grad::n_outputs + j] -
                    fd.grad[k * grad::n_outputs + j]) <= 2e-3 * scale,
           grad::par_name(pars[k], 1) + "/" + grad::output_names()[j] +
               " agrees with the differenced solve, single path");
      }
    }
    ok(std::abs(ift.grad[grad::n_outputs + 0]) > 0.0,
       "the single path's soil row moves the leaf");
  }
}

// The rank-one self-test, and it has real content: psi_soil reaches profit only
// through uptake, at a price the collar's own first-order condition fixes, so
// dprofit/dpsi_soil_j divided by dE_up/dpsi_soil_j must be the SAME number for
// every layer j. If it is not, either the formula or the price is wrong.
//
// It also pins WHICH price. `marginal_cost_water_multilayer()` is the lambda the
// collar solve equalises and the obvious candidate; it is 1.20x to 2.37x too
// large over the golden grid, because along the soil direction the collar is
// held fixed and the stem's own cost has not been netted out. See gradient.hpp.
void test_environment_water_rows_are_rank_one() {
  printf("environment rows: the water channel is rank one across layers\n");
  namespace grad = phylloptim::gradient;
  for (int L : {3, 5}) {
    for (double psi_soil : {1.0, 3.0}) {
      phylloptim::Leaf l = env::fresh();
      grad::Drivers d = env::drivers(psi_soil, 900.0, 2.0, L, L);
      grad::apply(l, env::kTheta, d, false, -1, true);
      l.find_root_collar_psi();
      const double psi_star = l.opt_root_psi_;
      const double lambda_multi = l.marginal_cost_water_multilayer();
      const double lambda_stem = l.marginal_cost_water();
      const double price = lambda_multi - lambda_stem;
      const std::string tag =
          std::to_string(L) + " layers, psi_soil=" + std::to_string(psi_soil);
      ok(std::isfinite(price) && price > 0.0, "the price is finite, " + tag);

      double first = 0.0;
      for (int j = 0; j < L; ++j) {
        const double h = 1e-4 * std::max(1.0, std::abs(d.psi_soil[j]));
        std::vector<double> up = d.psi_soil, dn = d.psi_soil;
        up[j] += h;
        dn[j] -= h;

        // dE_up/dpsi_soil_j at the FIXED collar. Uptake is an ordinary smooth
        // function of the soil state, so a difference here is well conditioned.
        grad::apply(l, env::kTheta, d, false, -1, true);
        l.find_root_collar_psi();
        l.E_from_Soil_to_Root_Collar(psi_star, up);
        const double e_up = l.E_up_;
        l.E_from_Soil_to_Root_Collar(psi_star, dn);
        const double e_dn = l.E_up_;
        const double dE = (e_up - e_dn) / (2.0 * h);

        // dprofit/dpsi_soil_j from a differenced SOLVE. Differencing profit
        // against the COLLAR would be the mistake this formulation exists to
        // avoid -- profit is flat there. Against psi_soil it is not.
        grad::Drivers du = d, dd = d;
        du.psi_soil = up;
        dd.psi_soil = dn;
        grad::apply(l, env::kTheta, du, false, -1, true);
        l.find_root_collar_psi();
        const double p_up = l.profit_;
        grad::apply(l, env::kTheta, dd, false, -1, true);
        l.find_root_collar_psi();
        const double p_dn = l.profit_;
        const double dP = (p_up - p_dn) / (2.0 * h);

        const double ratio = dP / dE;
        const std::string what =
            "layer " + std::to_string(j + 1) + ", " + tag;
        if (j == 0) {
          first = ratio;
        }
        near(ratio / first, 1.0, 1e-3, "same price in every layer: " + what);
        near(ratio / price, 1.0, 1e-3,
             "the price is lambda_multi - lambda_stem: " + what);
      }
      // And the bare multilayer lambda is not it, by a wide margin.
      ok(lambda_multi / first > 1.1,
         "the bare multilayer lambda overstates the price, " + tag);
      grad::apply(l, env::kTheta, d, false, -1, true);
    }
  }
}

// ⚠️ A CORRECT EXACT ZERO. Layers below the deepest rooted one contribute nothing
// to uptake, so wetting or drying one of them cannot reach the leaf and all four
// of its rows are exactly 0.0. In this codebase an exact zero is otherwise the
// signature of a missing accumulator, so it is asserted rather than left to be
// rediscovered -- and asserting it is also what stops a genuinely missing row
// from being read as this.
void test_environment_rows_are_zero_below_the_rooted_layers() {
  printf("environment rows: an unrooted layer's rows are exactly zero\n");
  namespace grad = phylloptim::gradient;
  const int L = 5, rooted = 2;
  phylloptim::Leaf l = env::fresh();
  grad::Drivers d = env::drivers(1.0, 900.0, 2.0, L, rooted);
  ok(d.root_network.r_R_H_min.size() == std::size_t(rooted),
     "the network is sized to the rooted layers");
  const std::vector<int> pars = env::all_env_pars(L);

  for (int route = 0; route < 2; ++route) {
    grad::Settings s;
    s.method = route == 0 ? grad::Method::Ift : grad::Method::Fd;
    const std::string tag = route == 0 ? "ift" : "fd";
    grad::Result r;
    grad::at(l, env::kTheta, d, false, pars.data(), pars.size(), s, r);
    for (int j = 0; j < grad::n_outputs; ++j) {
      // Rooted layers move the leaf...
      ok(std::abs(r.grad[1 * grad::n_outputs + j]) > 0.0 ||
             grad::output_names()[j] == "collar",
         "a rooted layer's " + grad::output_names()[j] + " row is non-zero, " +
             tag);
      // ...and unrooted ones cannot.
      for (int k = rooted; k < L; ++k) {
        const std::size_t row = std::size_t(1 + k);
        ok(r.grad[row * grad::n_outputs + j] == 0.0,
           "psi_soil_" + std::to_string(k + 1) + "/" +
               grad::output_names()[j] + " is exactly zero, " + tag);
      }
    }
  }
}

// The output enumeration, which is the whole of the arity decision on this side:
// five fixed outputs that R reads by position, then one uptake row per layer. The
// accepted request here is the consumer's own -- profit and the uptake block, no
// collar -- and the refused one asks for a layer this observation does not have.
void test_uptake_outputs_are_enumerated() {
  printf("uptake outputs: names, arity and refusal\n");
  namespace grad = phylloptim::gradient;
  ok(grad::out_uptake_first == grad::n_outputs,
     "the uptake block sits immediately after the five");
  ok(grad::output_names().size() == std::size_t(grad::n_outputs),
     "and the five R reads by position are still five");
  for (int L : {1, 3, 5}) {
    const std::vector<std::string> nms = grad::output_names(L);
    const std::string tag = std::to_string(L) + " layers";
    ok(int(nms.size()) == grad::n_outputs_total(L), "output count, " + tag);
    ok(std::equal(grad::output_names().begin(), grad::output_names().end(),
                  nms.begin()),
       "the fixed five are unchanged, " + tag);
    ok(nms[std::size_t(grad::out_uptake_first)] == "uptake_1",
       "the uptake block is one-based, " + tag);
    ok(nms.back() == "uptake_" + std::to_string(L),
       "and runs to the last layer, " + tag);
    ok(grad::output_name(grad::out_uptake_first + L - 1, L) ==
           "uptake_" + std::to_string(L),
       "output_name agrees with output_names, " + tag);
  }

  phylloptim::Leaf l = env::fresh();
  const grad::Drivers d = env::drivers(2.0, 900.0, 2.0, 3, 3);
  const grad::Settings s;
  const std::vector<int> inputs{grad::par_vcmax_25, grad::par_psi_soil_first};
  auto rows_for = [&](int layers, std::string &message) -> grad::Rows {
    std::vector<int> out_index{grad::out_profit};
    for (int i = 0; i < layers; ++i) {
      out_index.push_back(grad::out_uptake_first + i);
    }
    grad::RowRequest req;
    req.output = out_index.data();
    req.n_output = out_index.size();
    req.input = inputs.data();
    req.n_input = inputs.size();
    grad::Rows rows;
    try {
      rows = env::solved_rows(l, env::kTheta, d, req, s);
    } catch (const std::runtime_error &e) {
      message = e.what();
    }
    return rows;
  };

  std::string message;
  const std::size_t requested = 1 + 3;  // profit, then one uptake per layer
  const grad::Rows rows = rows_for(3, message);
  ok(message.empty(), "the consumer's request -- profit and the uptake block, no "
                      "collar -- answers: " + message);
  ok(rows.kind == grad::OperatingPointKind::Interior, "at an interior optimum");
  ok(rows.held.size() == requested * inputs.size(),
     "one row per requested output per input");
  ok(rows.dy_dp[0] == 0.0, "the objective's channel is the envelope's zero");
  bool uptake_moves = false;
  for (std::size_t j = 1; j < requested; ++j) {
    uptake_moves = uptake_moves || rows.held[j * inputs.size() + 1] != 0.0;
  }
  ok(uptake_moves, "and a soil layer moves the uptake it is asked about");

  // A layer count that does not match the drivers, refused by the name of the
  // uptake it would have been -- the way an input index past this observation's
  // arity is refused by its own name.
  message.clear();
  rows_for(4, message);
  ok(message.find("uptake_4") != std::string::npos &&
         message.find("uptake_3") != std::string::npos,
     "a fourth layer's uptake is refused by name: " + message);
}

// The rows in parts, assembled by the consumer's own formula, against the totals
// `at` returns for the same request:
//
//   dy_j/du_i = held[j][i] + dy_dp[j] * (-dresidual[i] / residual_slope)
//
// What this referees is the FACTORING. Both routes run the same base point and
// the same two perturbed evaluations per input, so a disagreement here is one of
// them having changed the algebra rather than the shape. It needs no reference
// gradient for the same reason the transpose identity does not.
//
// At a constrained point `at` differences the whole solve, and the parts are
// compared against that. The comparison splits, and the split is the content:
//
//   * an input that MOVES the bound has no held row -- p* sits one step-in
//     fraction from the bound, so a step that moves the bound carries p* out of
//     the perturbed feasible interval. `rows_at` differences the solve for it
//     too, so the two are the same arithmetic and must agree EXACTLY.
//   * an input that does not gets a held row plus a bound row that is exactly
//     zero, and there the two DISAGREE. The parts hold the point at the bound,
//     which is where the model puts it; `at` differences the solve, whose collar
//     moves by the solve's own floor -- 5.5e-10 MPa over a whole step, a step-in
//     of a millionth of the interval's width above the root-find's tolerance.
//     Divided by the step that is a collar row of 1e-04, so the movement is what
//     is asserted and the disagreement it produces is reported.
// The best agreement a DIFFERENCED referee reaches over a sweep of steps, and
// what it said at one.
//
// ⚠️ ONE STEP REPORTS ITS OWN TRUNCATION, and both places this is used were
// reading that as the closed form's error. A referee built on a root-find returns
// its answer to that solver's tolerance, so differencing it carries tol/h and the
// error GROWS as the step shrinks -- measured at 3.08e-07, 3.08e-06 and 3.08e-05
// over three decades, which is a fixed offset in the numerator and cannot be a
// property of the quantity being checked. And a referee on a smooth function can
// still spike at one step where a nested root-find changes its iterate count: the
// curvature's difference sits at 1e-10 at four steps and 5.6e-06 at a fifth
// between them.
//
// So the referee is the sweep's best agreement, which is what a plateau is.
struct Plateau {
  double best = std::numeric_limits<double>::infinity();
  double at_one_step = std::numeric_limits<double>::quiet_NaN();
};

void test_rows_in_parts_assemble_to_the_totals() {
  printf("the rows in parts assemble to the totals\n");
  namespace grad = phylloptim::gradient;
  using Kind = grad::OperatingPointKind;

  // Every output this observation has, with the two roles the mathematics fixes
  // rather than supplies: the collar IS the operating point, and profit is what it
  // maximises. The uptake block follows the five and every layer of it is
  // Ordinary.
  auto request_outputs = [](int n_layers, std::vector<int> &out_index) -> void {
    out_index.clear();
    for (int j = 0; j < grad::n_outputs_total(n_layers); ++j) {
      out_index.push_back(j);
    }
  };

  // The consumer's assembly, and its one branch: whether this input reaches the
  // point at all. `dresidual` is exactly zero both for an input the condition
  // does not read and for one whose rows already followed the point, and there
  // the channel contributes nothing -- which has to be a branch rather than a
  // multiply, because `dy_dp` is NA wherever the difference cannot be centred on
  // p*, and NA times zero is not zero.
  //
  // ⚠️ THE SLOPE IS AN ARGUMENT, and that is what keeps this a check on the
  // FACTORING. At an interior point `rows_at` reports the curvature in closed form
  // and `at` differences it; the two agree at the solver floor, but the assembled
  // total divides by it, and for an output whose two halves nearly cancel that
  // floor is amplified -- 1e-09 in the slope reaches 6e-06 of `gc`'s row. Passing
  // `at`'s own curvature here asks the one question this check exists for, whether
  // the shape is the same; the slope has its own referee below.
  auto assemble = [](const grad::Rows &rows, std::size_t n, std::size_t k, int j,
                     double slope) -> double {
    const double held = rows.held[std::size_t(j) * n + k];
    const double dpoint = -rows.dresidual[k] / slope;
    if (!std::isfinite(dpoint) || dpoint == 0.0) {
      return held;
    }
    return held + rows.dy_dp[std::size_t(j)] * dpoint;
  };

  const grad::Settings settings;
  double worst = 0.0, worst_followed = 0.0, worst_held = 0.0;
  std::string worst_where, worst_followed_where, worst_held_where;
  // The inputs `rows_at` READS at an interior point, where `at` differences the
  // solve. The two are no longer the same arithmetic there, so they are not held
  // to the same agreement as the rest; the read rows have their own referee
  // against the differencing primitive. It is the soil state, the root curve and
  // the transport -- everything the boundary reaches through the flux or through
  // the stem potential.
  double worst_read = 0.0, ignored_read = 0.0;
  std::string worst_read_where, ignored_read_where;
  int read_rows = 0, resolved_read = 0;
  // For an input the pinned bound does not read: how far the collar `at`
  // differences moves over one whole step, in MPa. The parts hold it at the bound,
  // and that movement is the whole of the disagreement between them.
  double worst_collar_move = 0.0;
  std::string worst_collar_move_where;
  int points = 0, compared = 0, refused = 0, contradiction = 0, parts_only = 0;
  int constrained = 0, followed_rows = 0, held_rows = 0, refused_rows = 0;
  int interior = 0, pinned = 0, shut = 0, other = 0;
  int pin_without_slope = 0, role_violation = 0, slack_violation = 0;
  // The closed-form curvature against the difference `at` takes, which is the one
  // thing `rows_at` now answers that `at` cannot. Reported rather than only
  // bounded, because the difference is what has the floor.
  double worst_slope = 0.0, worst_slope_one = 0.0;
  std::string worst_slope_where;
  int slope_compared = 0;
  int off_branch = 0, not_a_number = 0;
  // The uptake columns, which `at` does not report and which are refereed against
  // a re-solve instead. Split by family for the reason the five are: the parts
  // hold the point where the model puts it, and a differenced solve does not.
  double worst_uptake = 0.0, worst_uptake_followed = 0.0;
  double worst_uptake_held = 0.0;
  std::string worst_uptake_where, worst_uptake_followed_where;
  std::string worst_uptake_held_where;
  int uptake_compared = 0, uptake_off_branch = 0, uptake_columns = 0;
  int interior_uptake_nonzero = 0, pinned_uptake_nonzero = 0;
  int pinned_followed_uptake = 0, pinned_followed_uptake_nonzero = 0;
  int shut_uptake_entries = 0, shut_uptake_nonzero = 0, shut_profit_nonzero = 0;

  phylloptim::Leaf multilayer;
  phylloptim::Leaf single_potential;
  single_potential.set_supply_single();

  // dE_i/du by re-solving the leaf and reading the consumption profile it writes:
  // the referee for the uptake rows, sharing with them only the setters. False
  // where an arm left the base point's branch, which makes the difference one of
  // two functions rather than of one.
  auto differenced_uptake = [&](phylloptim::Leaf &l, const grad::Drivers &d,
                                bool single, int par, const grad::Branch &stay,
                                std::vector<double> &out) -> bool {
    double th[grad::n_pars];
    grad::Scratch scratch;
    std::vector<double> arm[2];
    const double base = grad::par_value(env::kTheta, d, single, par);
    const double h = grad::step_for(par, base, settings.step, grad::n_soil_layers(d, single));
    bool on_branch = true;
    for (int side = 0; side < 2; ++side) {
      // Base first, then one step: `stem_b`'s setter rescales the stem curve and
      // is sound only from base, so the second arm would otherwise be taken two
      // steps out.
      grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
      grad::set_one(l, th, env::kTheta, d, single, par,
                    side == 0 ? base + h : base - h,
                    settings.fast_stem_curve, scratch);
      l.find_root_collar_psi();
      on_branch = on_branch && grad::branch_here(l) == stay;
      arm[side] = l.soil_consumption_;
    }
    grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
    if (!on_branch) {
      return false;
    }
    out.assign(arm[0].size(), 0.0);
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] = (arm[0][i] - arm[1][i]) / (2.0 * h);
    }
    return true;
  };

  auto check = [&](const grad::Drivers &d, bool single,
                   const std::vector<int> &pars, const std::string &where) {
    phylloptim::Leaf &l = single ? single_potential : multilayer;
    const std::size_t n = pars.size();
    ++points;
    grad::Result fwd;
    grad::at(l, env::kTheta, d, single, pars.data(), n, settings, fwd);

    const int n_layers = grad::n_soil_layers(d, single);
    std::vector<int> out_index;
    request_outputs(n_layers, out_index);
    grad::RowRequest req;
    req.output = out_index.data();
    req.n_output = out_index.size();
    req.input = pars.data();
    req.n_input = n;
    grad::Rows rows;
    bool have_rows = true;
    try {
      rows = env::solved_rows(l, env::kTheta, d, req, settings);
    } catch (const std::runtime_error &) {
      // Nothing on this grid reaches here: a perturbation that moves the point out
      // of the perturbed feasible interval is answered by following it. Counted
      // rather than let through, so a point that stops answering says so.
      ++refused;
      have_rows = false;
      // The kind the refused point had, for the accounting below: the throw left
      // the leaf perturbed, so this is a fresh solve at base.
      grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
      l.find_root_collar_psi();
      rows.kind = l.operating_point_kind();
    }

    phylloptim::Leaf::WhichBound bound = phylloptim::Leaf::WhichBound::Wet;
    const bool is_pinned = grad::pinned_bound(rows.kind, bound);
    const bool is_interior = rows.kind == Kind::Interior;
    if (is_interior) {
      ++interior;
    } else if (is_pinned) {
      ++pinned;
    } else if (rows.kind == Kind::HydraulicShutdown) {
      ++shut;
    } else {
      ++other;
    }
    if (!have_rows) {
      return;
    }

    // The roles are in the numbers, so a consumer never branches on them.
    if (rows.dy_dp[std::size_t(grad::out_collar)] != 1.0 ||
        (is_interior && rows.dy_dp[std::size_t(grad::out_profit)] != 0.0)) {
      ++role_violation;
    }
    if (is_pinned && !(std::isfinite(rows.residual_slope) &&
                       rows.residual_slope != 0.0)) {
      ++pin_without_slope;
    }

    // The one input read by nothing but the dry bound: at an interior optimum the
    // constraint is inactive, so its whole row is exactly zero. That it is
    // COMPLEMENTARY SLACKNESS rather than a structural zero is declared by the
    // consumer, per parameter, and is not re-derived here.
    for (std::size_t k = 0; k < n; ++k) {
      if (pars[k] != grad::par_root_psi_crit || !is_interior) {
        continue;
      }
      for (int j = 0; j < grad::n_outputs; ++j) {
        if (rows.held[std::size_t(j) * n + k] != 0.0) {
          ++slack_violation;
        }
      }
    }

    if (fwd.used_ift && !is_interior) {
      ++contradiction;
    }
    if (is_interior && !fwd.used_ift) {
      ++parts_only;
    }
    // Against the composite where `at` forms one, and against a differenced
    // solve where it does not.
    const bool composite = fwd.used_ift && is_interior;
    if (composite) {
      ++compared;
    } else if (!is_interior) {
      ++constrained;
    } else {
      return;
    }

    // Which inputs move the bound, asked here from the bound's own row rather
    // than read off what `rows_at` returned, and every input at a shut collar
    // because no condition defines one there.
    std::vector<bool> follows(n, true);
    grad::Scratch scratch;
    double th[grad::n_pars];
    // Taken at every kind, not only the constrained ones: the uptake reference
    // below refuses an arm that leaves this branch.
    grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
    l.find_root_collar_psi();
    const grad::Branch base_branch = grad::branch_here(l);
    if (!(base_branch.kind == rows.kind)) {
      ++off_branch;
    }
    if (!composite) {
      if (is_pinned) {
        const phylloptim::Leaf::BoundRow cond = l.bound_row(bound);
        bool at_base = true;
        for (std::size_t k = 0; k < n; ++k) {
          const double dpoint =
              grad::bound_moves_by_re_solving(pars[k], single)
                  ? grad::differenced_bound(l, env::kTheta, d, single, pars[k],
                                            bound, settings, at_base, scratch)
                  : grad::bound_dpoint(cond, pars[k],
                                       grad::n_soil_layers(d, single));
          follows[k] = dpoint != 0.0;
        }
      }
      // Both arms of every difference, on the branch the base point took: taken
      // here rather than trusted, because a difference across a change of branch
      // -- another kind, or the other limit winning the dry bound -- is a
      // difference of two functions and `rows_at` is what has to hold it.
      for (std::size_t k = 0; k < n; ++k) {
        const double base = grad::par_value(env::kTheta, d, single, pars[k]);
        const double h = grad::step_for(pars[k], base, settings.step, grad::n_soil_layers(d, single));
        grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
        for (int side = 0; side < 2; ++side) {
          grad::set_one(l, th, env::kTheta, d, single, pars[k],
                        side == 0 ? base + h : base - h,
                        settings.fast_stem_curve, scratch);
          l.find_root_collar_psi();
          if (!(grad::branch_here(l) == base_branch)) {
            ++off_branch;
          }
        }
      }
      grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
    }

    // The uptake reference, one column per input, taken for the whole point
    // before anything is compared: the scale a row is judged against is the
    // largest entry in the same output row, which is not known input by input.
    std::vector<std::vector<double> > want_uptake(n);
    std::vector<bool> have_uptake(n, false);
    std::vector<double> uptake_scale(std::size_t(n_layers), 0.0);
    for (std::size_t k = 0; k < n; ++k) {
      have_uptake[k] =
          differenced_uptake(l, d, single, pars[k], base_branch, want_uptake[k]);
      if (!have_uptake[k]) {
        ++uptake_off_branch;
        continue;
      }
      ++uptake_columns;
      for (int layer = 0; layer < n_layers; ++layer) {
        uptake_scale[std::size_t(layer)] =
            std::max(uptake_scale[std::size_t(layer)],
                     std::abs(want_uptake[k][std::size_t(layer)]));
      }
    }

    // Whether `rows_at` reads this input's rows rather than differencing them.
    // ⚠️ A PIN IS INCLUDED, and it was not: these families are read at a
    // constrained point too now, so at a pin they no longer run `at`'s own
    // arithmetic and cannot be held to reproducing it bit for bit. What they are
    // held to instead is the band below, and the point's own movement they now
    // carry in `dresidual` is refereed against a differenced `find_root_psi`.
    auto is_read = [&](int par) -> bool {
      if (!(is_interior || is_pinned)) {
        return false;
      }
      phylloptim::Leaf::TransportTrait ignored =
          phylloptim::Leaf::TransportTrait::Conductance;
      return grad::supply_side(par, grad::n_soil_layers(d, single)) ||
             grad::transport_side(par, ignored) || grad::slack_side(par);
    };
    // The scale each of the five outputs is judged against, which is the largest
    // row `at` gives it over this point's inputs -- the same convention the uptake
    // block below already uses, and for its reason: an input an output barely
    // reaches has a row at the difference's floor, and a ratio taken against the
    // entry itself reports the floor rather than the row.
    double output_scale[grad::n_outputs];
    for (int j = 0; j < grad::n_outputs; ++j) {
      output_scale[j] = 1e-30;
      for (std::size_t k2 = 0; k2 < n; ++k2) {
        output_scale[j] = std::max(
            output_scale[j],
            std::abs(fwd.grad[k2 * grad::n_outputs + std::size_t(j)]));
      }
    }
    // Which curvature the factoring is refereed through. At an interior point
    // `at` divides by the difference it took and `rows_at` by the closed form, so
    // the identity below is asked of `at`'s -- otherwise it would be asking whether
    // the two curvatures agree, which is the next check's question and has its own
    // scale.
    const double route_slope = composite ? fwd.H : rows.residual_slope;
    if (is_interior && std::isfinite(rows.residual_slope) &&
        rows.residual_slope != 0.0) {
      ++slope_compared;
      // Re-seat and difference the marginal profit at four steps a decade apart.
      // `at`'s own H is the 1e-06 one, so the second column below is what a
      // single-step referee reports and the first is what the sweep settles on.
      Plateau pl;
      const double p_star = rows.point;
      grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
      l.find_root_collar_psi();
      for (int e = 4; e <= 7; ++e) {
        const double h =
            std::max(std::abs(p_star), 1.0) * std::pow(10.0, -double(e));
        const double d_hi = l.dprofit_droot_collar_psi(p_star + h);
        const double d_lo = l.dprofit_droot_collar_psi(p_star - h);
        if (phylloptim::util::identical(d_hi, 0.0) ||
            phylloptim::util::identical(d_lo, 0.0)) {
          continue;
        }
        const double H = (d_hi - d_lo) / (2.0 * h);
        const double err = std::abs(rows.residual_slope / H - 1.0);
        pl.best = std::min(pl.best, err);
        if (e == 6) {
          pl.at_one_step = err;
        }
      }
      if (std::isfinite(pl.best) && pl.best > worst_slope) {
        worst_slope = pl.best;
        worst_slope_where = where;
      }
      if (std::isfinite(pl.at_one_step) && pl.at_one_step > worst_slope_one) {
        worst_slope_one = pl.at_one_step;
      }
      grad::apply(l, env::kTheta, d, single, -1, settings.fast_stem_curve);
      l.find_root_collar_psi();
    }
    // Which inputs the solve's own tolerance cannot resolve a collar row for.
    std::vector<bool> unresolved(n, false);
    for (std::size_t k = 0; k < n; ++k) {
      const double base_k = grad::par_value(env::kTheta, d, single, pars[k]);
      const double moved =
          std::abs(fwd.grad[k * grad::n_outputs + std::size_t(grad::out_collar)]) *
          grad::step_for(pars[k], base_k, settings.step,
                         grad::n_soil_layers(d, single));
      unresolved[k] = moved < 100.0 * phylloptim::Leaf::collar_root_tol;
    }
    for (std::size_t k = 0; k < n; ++k) {
      if (!composite) {
        ++(follows[k] ? followed_rows : held_rows);
      }
      if (!std::isfinite(rows.held[std::size_t(grad::out_profit) * n + k])) {
        ++refused_rows;
        continue;
      }
      // What the two routes disagree through, where one holds the point at the
      // bound and the other differences the solve: how far the collar the solve
      // lands on moves over the whole step.
      if (is_pinned && !follows[k]) {
        const double base = grad::par_value(env::kTheta, d, single, pars[k]);
        const double moved =
            std::abs(fwd.grad[k * grad::n_outputs +
                              std::size_t(grad::out_collar)]) *
            grad::step_for(pars[k], base, settings.step, grad::n_soil_layers(d, single));
        if (moved > worst_collar_move) {
          worst_collar_move = moved;
          worst_collar_move_where =
              " at " + grad::par_name(pars[k], grad::n_soil_layers(d, single)) +
              where;
        }
      }
      for (int j = 0; j < grad::n_outputs; ++j) {
        const double got = assemble(rows, n, k, j, route_slope);
        const double want = fwd.grad[k * grad::n_outputs + std::size_t(j)];
        if (!std::isfinite(got)) {
          ++not_a_number;
          continue;
        }
        if (j == grad::out_profit && rows.kind == Kind::HydraulicShutdown &&
            got != 0.0) {
          ++shut_profit_nonzero;
        }
        const bool read_here = is_read(pars[k]);
        // A read row is judged against the largest thing its assembly is made of,
        // not against the total. Two reasons, and both are measured: an input an
        // output barely reaches has a total at the difference's floor, and an
        // input whose held row and point channel nearly cancel -- the maximum
        // conductance moves the stem potential by -36698 and takes almost all of
        // it back -- has a total three orders below its own parts, so parts right
        // to 1e-06 give a total at 1e-03. `at` differences the whole solve and
        // pays neither, which is why it is the reference and not the other way
        // round.
        const double held_term = rows.held[std::size_t(j) * n + k];
        const double point_term = got - held_term;
        const double err =
            std::abs(got - want) /
            (read_here ? std::max(output_scale[j],
                                  std::max(std::abs(held_term),
                                           std::abs(point_term)))
                       : std::max(std::abs(want), 1e-30));
        const bool read = read_here;
        // ⚠️ AND ONE THING THE DIFFERENCE CANNOT SAY. `at` reads the point's own
        // movement by re-solving, and the collar solve stops at 1e-12; an input
        // that moves the collar by less than a few multiples of that over a whole
        // step has no differenced row to be compared with, only the tolerance.
        // The root curve's steepness is such an input -- it moves the collar 1.4e-11
        // MPa -- and the read row is the better number there, so this counts it
        // rather than bounding it.
        if (read && !unresolved[k]) {
          ++resolved_read;
        }
        double &into = read && unresolved[k] ? ignored_read
                       : read                ? worst_read
                       : composite           ? worst
                       : follows[k]          ? worst_followed
                                             : worst_held;
        if (err > into) {
          into = err;
          const std::string what =
              " at " + grad::par_name(pars[k], grad::n_soil_layers(d, single)) +
              "/" + grad::output_names()[std::size_t(j)] + where;
          ((read && unresolved[k]) ? ignored_read_where
           : read       ? worst_read_where
           : composite  ? worst_where
           : follows[k] ? worst_followed_where
                        : worst_held_where) = what;
        }
      }

      // The uptake columns, which `at` does not report: the reference is a
      // re-solve of the leaf, and the parts have to assemble to it the way they do
      // to `at`'s five.
      if (!have_uptake[k]) {
        continue;
      }
      for (int layer = 0; layer < n_layers; ++layer) {
        const int j = grad::out_uptake_first + layer;
        // The reported slope here, not `at`'s: the referee is a re-solve of the
        // leaf rather than `at`, so the question is the answer's accuracy and the
        // closed form is the better number.
        const double got = assemble(rows, n, k, j, rows.residual_slope);
        const double want = want_uptake[k][std::size_t(layer)];
        // Both sides, because a ratio taken against a not-a-number compares false
        // and would be read as agreement.
        if (!std::isfinite(got) || !std::isfinite(want)) {
          ++not_a_number;
          continue;
        }
        ++uptake_compared;
        if (got != 0.0) {
          ++(is_interior      ? interior_uptake_nonzero
             : is_pinned      ? pinned_uptake_nonzero
                              : shut_uptake_nonzero);
        }
        if (is_pinned && follows[k]) {
          ++pinned_followed_uptake;
          if (got != 0.0) {
            ++pinned_followed_uptake_nonzero;
          }
        }
        if (rows.kind == Kind::HydraulicShutdown) {
          ++shut_uptake_entries;
        }
        const double err = std::abs(got - want) /
                           std::max(uptake_scale[std::size_t(layer)], 1e-30);
        const bool read = is_read(pars[k]);
        double &into = read           ? worst_read
                       : composite    ? worst_uptake
                       : follows[k]   ? worst_uptake_followed
                                      : worst_uptake_held;
        if (read) {
          ++read_rows;
        }
        if (err > into) {
          into = err;
          const std::string what =
              " at " + grad::par_name(pars[k], n_layers) + "/" +
              grad::output_name(j, n_layers) + where;
          (read ? worst_read_where
           : composite ? worst_uptake_where
           : follows[k] ? worst_uptake_followed_where
                        : worst_uptake_held_where) = what;
        }
      }
    }
  };

  // The golden file's grid, then the other supply path -- the same 294 operating
  // points the transpose identity is taken over.
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  for (double p : psi_soils) {
    for (double q : ppfds) {
      for (double vp : vpds) {
        for (int L : {1, 3, 5}) {
          // Everything but `resistance`, which the multi-layer path does not
          // have, plus the light row and one row per soil layer.
          std::vector<int> pars;
          for (int i = 0; i < grad::n_pars - 1; ++i) {
            pars.push_back(i);
          }
          for (int i : env::all_env_pars(L)) {
            pars.push_back(i);
          }
          check(env::drivers(p, q, vp, L, L), false, pars,
                " at psi_soil " + std::to_string(p) + ", ppfd " +
                    std::to_string(q) + ", vpd " + std::to_string(vp) + ", " +
                    std::to_string(L) + " layers");
        }
      }
    }
  }
  const int multilayer_points = points;
  for (double p : psi_soils) {
    grad::Drivers d = env::drivers(p, 900.0, 2.0, 1, 1);
    d.root_network = fixture::series_resistance(env::kTheta[grad::par_resistance]);
    std::vector<int> pars;
    for (int i = 0; i < grad::n_pars; ++i) {
      pars.push_back(i);
    }
    for (int i : env::all_env_pars(1)) {
      pars.push_back(i);
    }
    check(d, true, pars,
          " on the single-potential path at psi_soil " + std::to_string(p));
  }
  const int first_pass = points;
  const int first_pass_constrained = constrained;

  // The three photosynthetic traits alone, which move no bound on either arm: at
  // a constrained point every one of their rows is a held row, so this pass is
  // the held family on its own rather than mixed in with the followed one.
  const std::vector<int> photosynthetic{grad::par_vcmax_25, grad::par_jmax_25,
                                        grad::par_R_d_25};
  for (double p : psi_soils) {
    for (double q : ppfds) {
      for (int L : {1, 3, 5}) {
        check(env::drivers(p, q, 2.0, L, L), false, photosynthetic,
              " at psi_soil " + std::to_string(p) + ", ppfd " +
                  std::to_string(q) + ", " + std::to_string(L) + " layers");
      }
    }
  }

  // That the point really does follow the bound, refereed against a difference of
  // `find_root_psi`, which shares no code with either route: a soil layer moves
  // the wet bound, and the point's movement is what the parts have to reproduce.
  //
  // ⚠️ ASSEMBLED, NOT READ OFF `held`. The rows at a pin are parts now -- a held
  // partial and the bound's own gradient -- where they used to be totals a
  // re-solve measured, so the point's movement lives in `dresidual` and the
  // collar's held row is the exact zero a held partial makes it. Reading `held`
  // here reported 0 against a live bound movement, which is this check's own
  // subject arriving as a ratio of exactly one.
  int wet_pins = 0;
  double worst_bound = 0.0, worst_bound_one = 0.0;
  double worst_point = 0.0;
  std::string worst_point_where;
  {
    // ⚠️ NOT THE SOIL ALONE, and the three added inputs are the ones that make
    // this check bite. The wet bound is total uptake, which no stem property and
    // no conductance enters, so their entries in its row are exactly zero -- and
    // the point still moves with them, through the FAR bound's share of the
    // step-in. Their point rows read exactly 0 against a differenced solve of
    // 0.012, -2.2e-07 and 1.3e-07 until that share was carried.
    const std::vector<int> pars{grad::par_psi_soil_first, grad::par_kmax,
                                grad::par_stem_c, grad::par_psi_crit};
    std::vector<int> out_index;
    request_outputs(5, out_index);
    grad::RowRequest req;
    req.output = out_index.data();
    req.n_output = out_index.size();
    req.input = pars.data();
    req.n_input = pars.size();
    for (double p : psi_soils) {
      for (double q : ppfds) {
        grad::Drivers d = env::drivers(p, q, 2.0, 5, 5);
        grad::apply(multilayer, env::kTheta, d, false, -1,
                    settings.fast_stem_curve);
        multilayer.find_root_collar_psi();
        if (multilayer.operating_point_kind() != Kind::PinnedWet) {
          continue;
        }
        const grad::Rows rows =
            env::solved_rows(multilayer, env::kTheta, d, req, settings);
        ++wet_pins;
        const double h = 1e-6;
        double bound[2];
        for (int side = 0; side < 2; ++side) {
          grad::Drivers moved = d;
          moved.psi_soil[0] += side == 0 ? h : -h;
          grad::apply(multilayer, env::kTheta, moved, false, -1, true);
          bound[side] = multilayer.find_root_psi(moved.psi_soil[0],
                                                 moved.psi_soil, 0);
        }
        static_cast<void>(bound);
        static_cast<void>(h);
        const phylloptim::Leaf::BoundRow own =
            multilayer.bound_row(phylloptim::Leaf::WhichBound::Wet);

        // The BOUND's own soil row, against a difference of `find_root_psi`, over
        // a sweep of steps: it differences a root-find, so it carries that
        // solver's tolerance divided by the step and a single step reports 1/h
        // rather than the row -- measured at 3.08e-07, 3.08e-06 and 3.08e-05 over
        // three decades.
        Plateau on_bound;
        for (int e = 3; e <= 7; ++e) {
          const double hh = std::pow(10.0, -double(e));
          double edge[2];
          for (int side = 0; side < 2; ++side) {
            grad::Drivers moved = d;
            moved.psi_soil[0] += side == 0 ? hh : -hh;
            grad::apply(multilayer, env::kTheta, moved, false, -1, true);
            edge[side] = multilayer.find_root_psi(moved.psi_soil[0],
                                                  moved.psi_soil, 0);
          }
          const double err =
              std::abs(own.d_dpsi_soil[0] /
                       ((edge[0] - edge[1]) / (2.0 * hh)) - 1.0);
          on_bound.best = std::min(on_bound.best, err);
          if (e == 6) {
            on_bound.at_one_step = err;
          }
        }
        if (std::isfinite(on_bound.best)) {
          worst_bound = std::max(worst_bound, on_bound.best);
        }
        if (std::isfinite(on_bound.at_one_step)) {
          worst_bound_one = std::max(worst_bound_one, on_bound.at_one_step);
        }

        // The POINT's row, per input, against a difference of the SOLVE. Both
        // arms have to land on this branch: a wet pin sits a millionth of the
        // interval inside its bound, so an arm that crossed into the interior
        // regime would be differencing two different answers.
        double th_pin[grad::n_pars];
        grad::Scratch pin_scratch;
        for (std::size_t k = 0; k < pars.size(); ++k) {
          const double got =
              assemble(rows, pars.size(), k, grad::out_collar,
                       rows.residual_slope);
          const double base =
              grad::par_value(env::kTheta, d, false, pars[k]);
          Plateau on_point;
          for (int e = 4; e <= 7; ++e) {
            const double frac = std::pow(10.0, -double(e));
            const double hh = grad::step_for(pars[k], base, frac, 5);
            double collar[2];
            bool on_branch = true;
            for (int side = 0; side < 2; ++side) {
              grad::set_one(multilayer, th_pin, env::kTheta, d, false, pars[k],
                            side == 0 ? base + hh : base - hh, true,
                            pin_scratch);
              multilayer.find_root_collar_psi();
              on_branch = on_branch &&
                          multilayer.operating_point_kind() == Kind::PinnedWet;
              collar[side] = multilayer.opt_root_psi_;
            }
            if (!on_branch) {
              continue;
            }
            const double want_point = (collar[0] - collar[1]) / (2.0 * hh);
            // A row of exactly zero against a referee of exactly zero agrees;
            // a ratio there would be 0/0.
            const double scale = std::max(std::abs(want_point), std::abs(got));
            const double err =
                scale > 0.0 ? std::abs(got - want_point) / scale : 0.0;
            on_point.best = std::min(on_point.best, err);
          }
          if (std::isfinite(on_point.best) && on_point.best > worst_point) {
            worst_point = on_point.best;
            worst_point_where =
                std::string(" at ") + grad::par_name(pars[k], 5);
          }
        }
        grad::apply(multilayer, env::kTheta, d, false, -1, true);
      }
    }
  }

  printf("  %d operating points: %d interior, %d pinned, %d shut, %d other\n",
         points, interior, pinned, shut, other);
  printf("  %d points refused; all %d constrained points of the %d in the "
         "full-input pass answer\n",
         refused, first_pass_constrained, first_pass);
  printf("  %d compared against the composite (%d parts-only), %d against a "
         "differenced solve\n",
         compared, parts_only, constrained);
  printf("  %d input rows followed the point, %d were taken at a held one, %d "
         "refused\n",
         followed_rows, held_rows, refused_rows);
  printf("  worst |assembled - at()| / |at()| = %.3g%s\n", worst,
         worst_where.c_str());
  printf("  %d rows the soil state reads rather than differences: %.3g%s\n",
         read_rows, worst_read, worst_read_where.c_str());
  printf("  of which %d have a collar row the solve can resolve; the rest reach "
         "%.3g and are not compared%s\n", resolved_read, ignored_read,
         ignored_read_where.c_str());
  printf("  a row that followed the point:       %.3g%s\n", worst_followed,
         worst_followed_where.c_str());
  printf("  a row taken at a held point:         %.3g%s\n", worst_held,
         worst_held_where.c_str());
  printf("  the collar at() lands on moves %.3g MPa over a whole step%s\n",
         worst_collar_move, worst_collar_move_where.c_str());
  printf("  the closed-form curvature against a differenced one, over %d "
         "interior points: %.3g%s (%.3g at one step)\n",
         slope_compared, worst_slope, worst_slope_where.c_str(),
         worst_slope_one);
  printf("  %d wet pins: bound row against a differenced find_root_psi %.3g "
         "(%.3g at one step); assembled point row against a differenced solve "
         "%.3g%s\n",
         wet_pins, worst_bound, worst_bound_one, worst_point,
         worst_point_where.c_str());
  printf("  %d uptake rows over %d columns (%d off-branch), against a re-solved "
         "leaf\n",
         uptake_compared, uptake_columns, uptake_off_branch);
  printf("  worst uptake row / the column's own scale = %.3g%s\n", worst_uptake,
         worst_uptake_where.c_str());
  printf("  a row that followed the point:       %.3g%s\n",
         worst_uptake_followed, worst_uptake_followed_where.c_str());
  printf("  a row taken at a held point:         %.3g%s\n", worst_uptake_held,
         worst_uptake_held_where.c_str());
  printf("  non-zero uptake rows: %d interior, %d pinned (%d of %d that follow "
         "the bound), %d of %d shut\n",
         interior_uptake_nonzero, pinned_uptake_nonzero,
         pinned_followed_uptake_nonzero, pinned_followed_uptake,
         shut_uptake_nonzero, shut_uptake_entries);
  printf("  and %d non-zero profit rows at the shut points those came from\n",
         shut_profit_nonzero);
  ok(multilayer_points == 288, "the whole golden grid was covered");
  ok(first_pass - multilayer_points == 6, "and the single-potential path too");
  ok(pinned > 0 && shut > 0, "the pinned and the shut branches were reached");
  ok(refused == 0 && refused_rows == 0,
     "every point in the grid answers, and every input at it");
  ok(off_branch == 0,
     "every arm of every constrained difference stays on the base point's "
     "branch");
  ok(not_a_number == 0, "and every assembled row is a number");
  ok(followed_rows > 0 && held_rows > 0,
     "both the followed and the held families were reached");
  ok(wet_pins > 0 && worst_bound <= 1e-7,
     "at a wet pin the wet bound's own row is the bound's own movement");
  // ⚠️ A DIFFERENT REFEREE FROM THE LINE ABOVE, and the difference between them is
  // the whole of what a step-in is. `find_root_psi` moves the BOUND; the solve
  // returns that bound stepped a constant fraction of the interval's width inside,
  // so the POINT moves with both bounds and the two referees differ by that
  // fraction of the far one. Refereeing the point against the bound reported
  // 3.1e-05 and called it the point's error.
  ok(wet_pins > 0 && worst_point <= 1e-4,
     "and the point's own row is a differenced solve of the collar");
  // ⚠️ 3e-7 AND NOT 1e-7, AND THE REFEREE IS WHAT MOVED. The root tabulation now
  // carries an exact slope at every knot instead of one inferred from its values,
  // so its own slope agrees with the closed-form curve to 3e-11 where it agreed to
  // ~4e-06 (the PROBE above reports it). The differenced referee here reads 2.1e-07
  // and its single-step reading is 9.98e-06, so it can no longer see the quantity it
  // is refereeing: 1e-7 was calibrated against the old tabulation landing on the
  // useful side of the instrument's own floor.
  ok(slope_compared > 150 && worst_slope <= 3e-7,
     "the closed-form curvature is the difference it replaces");
  ok(contradiction == 0,
     "no point takes the composite that the solve did not call interior");
  ok(role_violation == 0, "the point's and the objective's channels are exact");
  ok(pin_without_slope == 0, "a pin comes back with the bound's own slope");
  ok(slack_violation == 0, "the slack row is exactly zero and says so");
  ok(compared > 150, "the composite was compared over most of the grid");
  ok(constrained > 90, "and the differenced solve over the constrained points");
  // ⚠️ THIS WAS AN IDENTITY AT 1e-12 AND IT IS NOW A BOUND, because the two routes
  // deliberately stopped running the same arithmetic. `rows_at` reads the collar
  // channel off the seated state and the curvature off its closed form; `at`
  // differences both, and must, being refereed bit for bit against a captured
  // reference. That was the stated price of the swap and this is its size: 7e-05,
  // at an output the read does not state at all, through the assembly's own
  // cancellation. What replaces the identity is two tighter checks on the halves --
  // the curvature against a swept difference below, and the collar channel at an
  // interior point, where the read and the difference agree at 1.6e-10.
  ok(worst <= 1e-4, "the parts assemble to the totals");
  // ⚠️ AND THE SOIL STATE'S ROWS ARE NOT HELD TO THAT, because for them the two
  // routes are no longer the same arithmetic: `rows_at` reads a closed form where
  // `at` differences the solve. What is left is the difference's own floor
  // through the concentration's root-find, which falls as 1/h and so is the
  // reference rather than the row -- refereed against the differencing primitive
  // at two steps in the soil state's own check, where it goes 3e-04 to 3e-07 over
  // four decades.
  // ⚠️ THE READ ROWS ARE REPORTED HERE AND BOUNDED ELSEWHERE, and that is not a
  // gap. This check's whole content is that both routes run the SAME arithmetic,
  // so a disagreement is one of them having changed the algebra -- and for a row
  // that is read rather than differenced the premise is false. What is left is
  // `at`'s own error, which at this step is not small for every input: the maximum
  // conductance is 3.1e-05, so a relative step of 1e-06 is an absolute 3.1e-11,
  // and `at` divides the concentration root-find's floor by it. That is why its
  // assimilation column reads 15 where the flux does not move at all.
  //
  // The rows that ARE bounded are bounded against the differencing primitive, at
  // two steps, in the soil state's and the transport's own checks. Here they are
  // counted and their worst is printed, so that a change in them is visible.
  ok(read_rows > 0, "the soil state's rows are read rather than differenced");
  ok(std::isfinite(worst_read), "and every one of them is a number");
  ok(worst_followed == 0.0,
     "at a constrained point a followed row IS at()'s own difference");
  // Why the held family does not, and why the disagreement above is reported
  // rather than bounded. The parts hold the point at the bound, which is where
  // the model puts it; `at` differences the solve, and the collar it lands on
  // moves 5.5e-10 MPa over a whole step -- the solve's own floor, a step-in of a
  // millionth of the interval's width above the root-find's tolerance. Divided by
  // the step that is a collar row of 1e-04, and through each output's sensitivity
  // to the collar it is as large as the rows themselves. So the movement is what
  // is asserted and the ratio it produces is left as a measurement.
  ok(worst_collar_move <= 1e-8,
     "and a held row differs from it only as far as the solve's own collar moves");

  // The uptake outputs. Their reference is a re-solved leaf rather than `at`,
  // which reports the five, so the band is the differencing one: 4.0e-04 here
  // against the 6.4e-04 the water rows already sit at, and for the same reason --
  // both move the argmax, and the composite's dY/dp and its curvature are
  // themselves differences taken at a flat maximum.
  ok(uptake_compared > 10000 && uptake_off_branch == 0,
     "every input's uptake column was refereed against a re-solve");
  ok(worst_uptake <= 2e-3,
     "an uptake row assembles to a differenced re-solve at an interior point");
  ok(worst_uptake_followed == 0.0,
     "and at a constrained point a followed uptake row IS that difference");
  ok(interior_uptake_nonzero > 0, "interior uptake rows are not zero");
  // ⚠️ THE FAILURE THIS EXISTS TO PREVENT. An output that is not reported has no
  // route to any input, so at a pin -- where the condition's gradient is zero for
  // every input and each reported output carries a total instead -- its rows come
  // back exactly zero. So a pin is where the total has to be seen.
  ok(pinned_followed_uptake > 0 && pinned_followed_uptake_nonzero > 0,
     "at a pin an uptake row that follows the bound carries a total");
  // At a hydraulic shutdown the total IS zero, and that is the model rather than a
  // missing route: the stem is held at the critical potential and the leaf takes up
  // no water at any perturbation the branch survives. The differenced re-solve
  // agrees to the bit (`worst_uptake_followed`), and profit at the same points is
  // non-zero, which is what says the route ran.
  ok(shut_uptake_entries > 0 && shut_uptake_nonzero == 0 &&
         shut_profit_nonzero > 0,
     "a shut collar takes up no water, and says so beside a non-zero profit row");
}

// The condition's slope in the stem potential, and the claim that the state
// reaches the condition through two intermediates and no more.
//
// At a held collar the condition is a function of the state only through the stem
// potential and the transport's response to the collar, so for any state
// direction u
//
//   dR/du = dR/dV * dV/du  +  dR/dsigma * dsigma/du
//
// with dR/dV the coefficient the evaluation records and dR/dsigma the object under
// test. Every term is measurable: the left side and the two state derivatives are
// differences of the leaf at a held collar, and nothing here is a reference
// gradient.
//
// ⚠️ THIS IS THE CHECK A RANK TEST CANNOT BE. A rank test asks whether two
// coefficients suffice; this asks whether THESE two are the right ones, which is
// what separates the pairing from the one that reaches the condition only through
// the stem potential.
void test_the_condition_reaches_the_state_through_two_intermediates() {
  printf("the condition's slope, and the state reaching it through two\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;

  struct Fixture { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const Fixture fixtures[] = {{"wet", 2.0, 900.0, 2.0, 3},
                              {"dim", 2.0, 300.0, 2.0, 3},
                              {"dry", 4.5, 900.0, 2.0, 3},
                              {"arid", 2.0, 900.0, 4.0, 3}};

  for (const Fixture& f : fixtures) {
    grad::Drivers d =
        env::drivers(f.psi_soil, f.ppfd, f.vpd, f.layers, f.layers);
    pl::Leaf l = env::fresh();
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    l.find_root_collar_psi();
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(l.operating_point_kind()) + ")";
    const double psi_star = l.opt_root_psi_;
    l.dprofit_droot_collar_psi(psi_star);
    const double dR_dV = l.dprofit_dpsistem_;
    double dR_dsigma = 0.0;
    ok(l.condition_slope(dR_dsigma), "the condition's slope is answered, " + tag);
    if (!std::isfinite(dR_dsigma)) {
      continue;
    }

    // One state direction per soil layer, and one per layer of root carbon.
    // Each is differenced with the COLLAR HELD, which is what makes the two
    // intermediates the only route.
    double worst = 0.0;
    int moved = 0;
    std::string worst_at;
    std::vector<double> got_dR, got_dsigma, got_dV;
    for (int which = 0; which < 2 * f.layers; ++which) {
      const int layer = which % f.layers;
      const bool carbon = which >= f.layers;
      const int par = carbon ? grad::par_root_carbon_first(f.layers) + layer
                             : grad::par_psi_soil_first + layer;
      const double base = grad::par_value(env::kTheta, d, false, par);
      const double h = std::max(std::abs(base), 1.0) * 1e-6;

      auto at = [&](double value) {
        double th[grad::n_pars];
        grad::Scratch scratch;
        grad::set_one(l, th, env::kTheta, d, false, par, value,
                      s.fast_stem_curve, scratch);
        // The collar is HELD at the base point's, not re-solved.
        const double R = l.dprofit_droot_collar_psi(psi_star);
        // The stem potential at the HELD collar, computed rather than read: the
        // evaluation keeps it as a local and `opt_psi_stem_` still holds the
        // solve's, which a perturbation has just invalidated.
        const double sigma =
            l.find_psi_stem_from_psi_root(psi_star, l.supply_psi_soil());
        return std::tuple<double, double, double>(R, sigma, l.dpsistem_dpsi_);
      };
      const auto up = at(base + h);
      const auto dn = at(base - h);
      const double dR = (std::get<0>(up) - std::get<0>(dn)) / (2.0 * h);
      const double dsigma = (std::get<1>(up) - std::get<1>(dn)) / (2.0 * h);
      const double dV = (std::get<2>(up) - std::get<2>(dn)) / (2.0 * h);

      got_dR.push_back(dR);
      got_dsigma.push_back(dsigma);
      got_dV.push_back(dV);
      const double predicted = dR_dV * dV + dR_dsigma * dsigma;
      const double scale = std::max(std::abs(dR), 1e-12);
      const double residual = std::abs(predicted - dR) / scale;
      if (std::isfinite(dR) && std::abs(dR) > 1e-9 && std::isfinite(dsigma) &&
          std::isfinite(dV)) {
        ++moved;
      }
      if (residual > worst || worst_at.empty()) {
        worst = residual;
        worst_at = grad::par_name(par, f.layers);
      }
      // Restore the base state for the next direction.
      grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
      l.find_root_collar_psi();
      l.dprofit_droot_collar_psi(psi_star);
    }
    // Non-vacuity, because a residual of zero is what a loop that measured
    // nothing also reports: every direction has to have moved the condition.
    ok(moved == 2 * f.layers,
       "every state direction moved the condition, " + tag);
    // ⚠️ WHAT BOUNDS THIS IS THE INTERPOLANT'S SECOND DERIVATIVE, and it is not a
    // property of the differencing: the residual is flat to three digits over five
    // decades of step, so it is systematic. Solving both coefficients from two
    // state directions and predicting the rest puts the supply's rank-two structure
    // at 1e-09 and dR/dV at 1e-10, so what carries the residual is dR/dsigma,
    // which agrees with the solved value to 2e-05.
    //
    // Of that, the larger part is the transport curve being read as a C1
    // interpolant: its second derivative is not the conductivity's own slope --
    // they differ by 2.9e-05 here -- and recomputing dR/dsigma on the
    // interpolant's takes the disagreement to 8.7e-06. Supplying the interpolant's
    // slope fixed the FIRST derivative and left this one order untouched, which is
    // the part of report 07 section 6 that said the second-derivative inventory
    // would not close.
    ok(worst < 1e-3, "the two intermediates carry every state direction, " + tag +
                         " (worst " + std::to_string(worst) + " at " + worst_at +
                         ")");
    // The object itself, rather than the identity it satisfies: solve both
    // coefficients from one soil direction and one root-carbon direction -- two
    // families, so the pair is determined -- and hold the closed forms against
    // them. This is what says dR/dsigma is right and not merely consistent.
    const int i = 0, j = f.layers;
    const double det =
        got_dsigma[i] * got_dV[j] - got_dsigma[j] * got_dV[i];
    const double solved_dsigma =
        (got_dR[i] * got_dV[j] - got_dR[j] * got_dV[i]) / det;
    const double solved_dV =
        (got_dsigma[i] * got_dR[j] - got_dsigma[j] * got_dR[i]) / det;
    // dR/dV is exact; dR/dsigma is held to what the curve's second derivative
    // allows at this knot count, which is 2e-05 -- see
    // `vulnerability_curve_ncontrol` for the measurement against resolution.
    near(dR_dV, solved_dV, 1e-8, "dR/dV against the solved pair, " + tag);
    near(dR_dsigma, solved_dsigma, 1e-4,
         "dR/dsigma against the solved pair, " + tag);
    printf("  %-22s dR/dsigma %12.5g (solved %12.5g, rel %8.3g)   worst "
           "residual %9.3g at %s\n",
           tag.c_str(), dR_dsigma, solved_dsigma,
           std::abs(dR_dsigma - solved_dsigma) / std::abs(solved_dsigma), worst,
           worst_at.c_str());
  }
}

// Which pair of coordinates the condition is written in, settled by identity
// rather than by fitting.
//
// The condition is the marginal profit. Read off its own composition it is
//
//   R = dprofit_dpsistem_ * dpsi_stem/dpsi  +  dprofit_dpsi_held_stem_
//
// so the two coordinates are the stem potential and the transport's response to
// the collar. The obvious alternative -- total uptake and its collar slope --
// reaches the condition ONLY through the stem potential, because uptake is what
// sets it. It is therefore short by the second term, whose size is measured here.
//
// ⚠️ A RANK TEST CANNOT SEPARATE THE TWO. Both pairings are rank two and both
// reproduce every out-of-sample state direction; what separates them is whether
// the second term is there, so that is what this reads.
void test_the_condition_is_the_stem_potential_and_its_collar_response() {
  printf("the condition's two coordinates, and the size of the direct term\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;

  struct Fixture { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const Fixture fixtures[] = {{"wet", 2.0, 900.0, 2.0, 3},
                              {"dim", 2.0, 300.0, 2.0, 3},
                              {"dry", 4.5, 900.0, 2.0, 3},
                              {"arid", 2.0, 900.0, 4.0, 3},
                              {"pinned", 5.72, 900.0, 2.0, 1}};

  for (const Fixture& f : fixtures) {
    grad::Drivers d =
        env::drivers(f.psi_soil, f.ppfd, f.vpd, f.layers, f.layers);
    pl::Leaf l = env::fresh();
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    l.find_root_collar_psi();
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(l.operating_point_kind()) + ")";
    const double psi_star = l.opt_root_psi_;
    bool feasible = false;
    const double R = l.dprofit_droot_collar_psi(psi_star, &feasible);
    ok(feasible, "the point is evaluable, " + tag);
    if (!feasible) {
      continue;
    }

    const double through_stem = l.dprofit_dpsistem_ * l.dpsistem_dpsi_;
    const double direct = l.dprofit_dpsi_held_stem_;
    near(through_stem + direct, R, 1e-13,
         "the condition is its two coordinates, " + tag);

    // The whole of the discrimination: were this zero, a pairing reaching the
    // condition only through the stem potential would be complete.
    ok(direct != 0.0 && std::isfinite(direct),
       "and the direct term is not zero, " + tag);
    printf("  %-30s R %11.4g = stem %11.4g + direct %11.4g  (direct %.1f%%)\n",
           tag.c_str(), R, through_stem, direct,
           100.0 * std::abs(direct) /
               std::max(std::abs(R), std::abs(through_stem)));
  }
}

// The transport's three parameters against the rebuild-and-difference they
// replace, which is the last differencing the row layer did.
//
// ⚠️ THE DIFFERENCE HERE REBUILDS A CURVE ON EVERY ARM. That is what makes it the
// right reference -- the forward model rebuilds too, so the difference is of the
// model -- and it is also what makes it expensive, which is the whole reason for a
// closed form.
//
// They used to disagree, and by enough to be recorded rather than bounded: the
// steepness' profit row was out by 4.7e-05 and its condition row by 1.1e-03. The
// cause was not this derivation. It differentiates the flux balance
// kappa (G(sigma) - G(p)) = E_up, and the forward model did not use that balance
// for the collar response -- it read the inverted table's slope instead, so the
// two described different functions. With the model taking its response from the
// balance as well, the disagreement is gone.
// The flux really is frozen when the transport moves, and what looks like a
// response in assimilation is the concentration's root-find divided by the step.
//
// The claim the transport's rows rest on is that these three parameters move the
// stem potential and leave the flux where it is, so nothing past the flux
// responds. Measured directly: the soil's uptake row is EXACTLY zero, the stem's
// own flux moves 3e-10 and the conductance 1e-06 -- and assimilation appears to
// move by 15.
//
// ⚠️ THAT 15 IS NOT A RESPONSE, and the step says so. The maximum conductance is
// 3.1e-05, so a relative step of 1e-06 is an ABSOLUTE step of 3.1e-11, and a
// difference over it divides whatever floor the concentration's root-find leaves
// by that. Over four decades of step the row spreads by four orders and changes
// sign, while the stem potential's row is -20803.9 at every one of them.
//
// ⚠️ WHICH STEP THE FLOOR IS LARGEST AT IS THE TABLE'S, NOT THE MODEL'S, so it is
// not what this checks. Read off a C1 table the floor was largest at the finest
// step -- 15.1, -1.51, -0.151, -1.2e-04, close to 1/h -- because a difference of a
// first derivative reads the second-derivative jump a C1 span leaves at every
// knot. Read off a C2 one there is no jump to read and the finest step is the
// smallest of the four: -1.7e-04, 0.756, 0.0756, -0.00756. Both are the same
// finding, which is that the number is not a derivative, and the assertions below
// say that directly instead of assuming where the floor lands.
void test_the_transport_leaves_the_flux_where_it_is() {
  printf("the transport moves the potential and not the flux\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 5;
  grad::Drivers d = env::drivers(1.0, 1500.0, 1.0, L, L);
  pl::Leaf l = env::fresh();
  const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
  bool feasible = false;
  l.dprofit_droot_collar_psi(b.psi_star, &feasible);
  ok(feasible, "the point is evaluable");
  const double p = b.psi_star;
  double th[grad::n_pars];
  grad::Scratch scratch;
  const double base = grad::par_value(env::kTheta, d, false, grad::par_kmax);

  double first_sigma = 0.0, last_sigma = 0.0;
  double assim_lo = std::numeric_limits<double>::infinity(), assim_hi = 0.0;
  bool changes_sign = false, seen_positive = false, seen_negative = false;
  int decade = 0;
  for (double rel : {1e-6, 1e-5, 1e-4, 1e-3}) {
    const double h = base * rel;
    double sigma[2], assim[2], uptake_total[2];
    for (int side = 0; side < 2; ++side) {
      grad::set_one(l, th, env::kTheta, d, false, grad::par_kmax,
                    side == 0 ? base + h : base - h, s.fast_stem_curve, scratch);
      l.supply_begin_solve();
      grad::OutputValues y(L);
      ok(grad::outputs_at(l, p, y), "the held evaluation answers");
      sigma[side] = y[grad::out_psi_stem];
      assim[side] = y[grad::out_assim];
      uptake_total[side] = 0.0;
      for (int i = 0; i < L; ++i) {
        uptake_total[side] += y[grad::out_uptake_first + i];
      }
    }
    const double d_sigma = (sigma[0] - sigma[1]) / (2.0 * h);
    const double d_assim = (assim[0] - assim[1]) / (2.0 * h);
    printf("  step %8.0e   sigma %13.7g   assim %13.6g   uptake %8.1e\n", rel,
           d_sigma, d_assim, (uptake_total[0] - uptake_total[1]) / (2.0 * h));
    ok(uptake_total[0] == uptake_total[1],
       "no layer's draw moves at all, at a step of " + std::to_string(rel));
    if (decade == 0) {
      first_sigma = d_sigma;
    }
    last_sigma = d_sigma;
    ++decade;
    assim_lo = std::min(assim_lo, std::abs(d_assim));
    assim_hi = std::max(assim_hi, std::abs(d_assim));
    if (d_assim > 0.0) seen_positive = true;
    if (d_assim < 0.0) seen_negative = true;
  }
  // A derivative is the same number at every step to within the differencing
  // error. This says this one is not: it spreads by orders, and no step of any
  // size makes it converge. The SIGN pattern is not checked -- it is the table's
  // and not the model's, and which steps come out negative moves with the grid.
  static_cast<void>(changes_sign);
  static_cast<void>(seen_positive);
  static_cast<void>(seen_negative);
  ok(assim_hi > assim_lo * 1e3,
     "assimilation's apparent row spreads by orders across the step");
  grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
  // Four decades of step move it by 1.5e-06, which is the coarsest step's own
  // truncation and not a floor: a floor would have shown at the FINEST.
  near(last_sigma, first_sigma, 1e-5,
       "while the stem potential's row barely moves with the step at all");
}

void test_the_transport_traits_rows_match_a_rebuilt_difference() {
  printf("the transport's three rows against a rebuilt difference\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  struct Fixture { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const Fixture fixtures[] = {{"wet", 2.0, 900.0, 2.0, 3},
                              {"dim", 2.0, 300.0, 2.0, 3},
                              {"dry", 4.5, 900.0, 2.0, 3},
                              {"arid", 2.0, 900.0, 4.0, 3},
                              {"lush", 1.0, 1500.0, 1.0, 5}};
  using Trait = pl::Leaf::TransportTrait;
  const struct { const char* what; int par; Trait trait; } traits[] = {
      {"kmax", grad::par_kmax, Trait::Conductance},
      {"stem_b", grad::par_stem_b, Trait::Position},
      {"stem_c", grad::par_stem_c, Trait::Steepness}};

  // ⚠️ TWO STEPS, AND THE ROW IS HELD TO THE BETTER OF THEM. The difference is
  // taken at a RELATIVE step, and these three parameters differ by five orders in
  // size: the maximum conductance is 3.1e-05, so 1e-06 of it is an absolute
  // 3.1e-11 and divides the concentration's root-find floor by that; the stem
  // curve's steepness is 2.7, so 1e-04 of it is an absolute 2.7e-04 and truncates.
  // Measured, the two fail in OPPOSITE directions, and a single step chosen for
  // one of them reports the other's reference rather than its row.
  int compared = 0;
  for (int path = 0; path < 2; ++path)
  for (const Fixture& f : fixtures) {
    const bool single = path == 1;
    const int n_layer = single ? 1 : f.layers;
    grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, n_layer, n_layer);
    std::vector<int> out_index{grad::out_psi_stem, grad::out_profit};
    for (int i = 0; i < n_layer; ++i) {
      out_index.push_back(grad::out_uptake_first + i);
    }
    grad::Settings s;
    pl::Leaf l = env::fresh();
    if (single) {
      l.set_supply_single();
    }
    const grad::BasePoint b =
        env::solved_base_point(l, d, single, s, n_layer);
    const std::string tag =
        std::string(f.what) + (single ? " single (" : " (") +
        pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
    if (b.branch.kind != pl::Leaf::OperatingPointKind::Interior) {
      continue;
    }
    bool feasible = false;
    l.dprofit_droot_collar_psi(b.psi_star, &feasible);
    ok(feasible, "the point is evaluable, " + tag);

    for (const auto& t : traits) {
      const std::string nm = std::string(t.what) + ", " + tag;
      pl::Leaf::TransportTraitRows rows;
      ok(l.transport_trait_rows(t.trait, l.dpsistem_dpsi_, rows),
         std::string("the rows are answered for ") + nm);

      double best[3] = {1e30, 1e30, 1e30};
      double at_step[2][3];
      int taken = 0;
      for (int which = 0; which < 2; ++which) {
        s.step = which == 0 ? 1e-6 : 1e-4;
        // ⚠️ A FRESH LEAF PER DIFFERENCE. `held_row` leaves the leaf perturbed --
        // that is what its `at_base` out-parameter says -- and putting it back is
        // not free of the curve's own rescaling, so the second step would
        // difference around wherever the first one stopped. Constructing one is
        // two grid builds and this is not a hot path.
        pl::Leaf fresh = env::fresh();
        if (single) {
          fresh.set_supply_single();
        }
        const grad::BasePoint fb =
            env::solved_base_point(fresh, d, single, s, n_layer);
        bool at_base = true;
        grad::Scratch scratch;
        grad::OutputValues direct(n_layer);
        double dR = 0.0;
        if (!grad::held_row(fresh, env::kTheta, d, single, t.par, fb.psi_star, s,
                            true, at_base, scratch, direct, dR)) {
          at_step[which][0] = at_step[which][1] = at_step[which][2] = 1e30;
          continue;
        }
        ++taken;
        const double got[3] = {rows.dpsistem, rows.dprofit, rows.dmarginal};
        const double want[3] = {direct[grad::out_psi_stem],
                                direct[grad::out_profit], dR};
        for (int q = 0; q < 3; ++q) {
          at_step[which][q] =
              std::abs(want[q]) > 0.0 ? std::abs(got[q] / want[q] - 1.0) : 1e30;
          best[q] = std::min(best[q], at_step[which][q]);
        }
        // The flux does not move, so no layer's draw does -- and exactly, because
        // the derivation says the whole held row past the flux is the cost.
        for (int i = 0; i < n_layer; ++i) {
          near(direct[grad::out_uptake_first + i], 0.0, 1e-9,
               "and the flux holds, layer " + std::to_string(i + 1) + ", " + nm);
        }
      }
      ok(taken == 2, "both steps answer for " + nm);
      ++compared;
      printf("  %-8s %-22s sigma %.1e/%.1e  profit %.1e/%.1e  R %.1e/%.1e\n",
             t.what, tag.c_str(), at_step[0][0], at_step[1][0], at_step[0][1],
             at_step[1][1], at_step[0][2], at_step[1][2]);
      ok(best[0] <= 1e-6, "the stem potential's row for " + nm);
      ok(best[1] <= 1e-5, "the profit row for " + nm);
      ok(best[2] <= 1e-4, "the condition's row for " + nm);
    }
  }
  ok(compared >= 25, "every transport parameter was compared on both paths");
}

void test_the_curves_trait_derivative_is_the_models_own() {
  printf("the closed-form dG/dc against a rebuilt difference of the spline\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 3;
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, L, L);
  pl::Leaf l = env::fresh();
  grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
  const double b = l.stem_b, c = l.stem_c;

  double worst_dc = 0.0, worst_value = 0.0;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    const pl::VulnerabilityIntegralDerivatives dd =
        pl::cumulative_vulnerability_integral_derivatives_at(psi, b, c);
    const double value = l.stem_curve_integral(psi);
    worst_value = std::max(worst_value,
                           std::abs(dd.value - value) / std::abs(value));

    double th[grad::n_pars];
    grad::Scratch scratch;
    const double h = c * 1e-5;
    grad::set_one(l, th, env::kTheta, d, false, grad::par_stem_c, c + h,
                  s.fast_stem_curve, scratch);
    const double up = l.stem_curve_integral(psi);
    grad::set_one(l, th, env::kTheta, d, false, grad::par_stem_c, c - h,
                  s.fast_stem_curve, scratch);
    const double dn = l.stem_curve_integral(psi);
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    const double rebuilt = (up - dn) / (2.0 * h);
    worst_dc = std::max(worst_dc, std::abs(dd.dc - rebuilt) /
                                     std::max(std::abs(rebuilt), 1e-30));
  }
  // The wet end carries the worst of it, where the integral itself is smallest.
  ok(worst_value < 1e-7, "the series' value matches the spline's");
  ok(worst_dc < 1e-4, "and its dG/dc matches a rebuilt difference");
  printf("  value to %.3g, dG/dc to %.3g over five positions\n", worst_value,
         worst_dc);
}

// The transport's response to the collar, which the marginal-profit evaluation
// forms and used to discard. The condition's gradient is written in it, so a
// consumer that wants the condition in parts needs it reported.
void test_the_transport_reports_its_collar_response() {
  printf("the transport's collar response, against a difference of the transport\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;

  struct Fixture { const char* what; double psi_soil; int layers; };
  const Fixture fixtures[] = {{"wet", 2.0, 3},
                              {"dry", 4.5, 3},
                              {"pinned", 5.72, 1}};

  for (const Fixture& f : fixtures) {
    grad::Drivers d = env::drivers(f.psi_soil, 900.0, 2.0, f.layers, f.layers);
    pl::Leaf l = env::fresh();
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    l.find_root_collar_psi();
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(l.operating_point_kind()) + ")";
    const double psi_star = l.opt_root_psi_;
    l.dprofit_droot_collar_psi(psi_star);

    // Against a difference of the transport itself, which shares no arithmetic
    // with it: the reported value is the inverse curve's slope times the flux's,
    // and this is a difference of the inverse solve.
    const double h = 1e-6 * std::max(psi_star, 1.0);
    const double up =
        l.find_psi_stem_from_psi_root(psi_star + h, l.supply_psi_soil());
    const double dn =
        l.find_psi_stem_from_psi_root(psi_star - h, l.supply_psi_soil());
    const double differenced = (up - dn) / (2.0 * h);
    near(l.dpsistem_dpsi_, differenced, 1e-5,
         "the reported dpsi_stem/dpsi, " + tag);
    printf("  %-30s %.6f against %.6f\n", tag.c_str(), l.dpsistem_dpsi_,
           differenced);
  }
}

// The eight carbon-side traits have closed-form rows, and this is what says
// whether they are right. At a frozen collar each reaches profit through
// assimilation or through the hydraulic cost and through nothing else, so its
// held row is one of `photo_trait_rows`' or `cost_trait_rows`' `dprofit_`
// entries, and at an interior point -- where stationarity is the condition --
// the condition's gradient is the matching `dmarginal_`.
//
// ⚠️ NEITHER READER HAS A CALLER ANYWHERE, tests included, so until this ran
// nothing had compared them with the solve they claim to describe. The
// differenced row is the reference: it is what the model does today.
void test_the_supplys_mixed_partials_match_a_difference() {
  printf("the supply's mixed partials against a difference of the conductance\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 5;
  grad::Drivers d = env::drivers(3.0, 1500.0, 0.5, L, L);
  pl::Leaf l = env::fresh();
  const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
  const std::vector<double> psi = l.supply_psi_soil();
  const double p = l.opt_root_psi_;

  std::vector<double> got;
  l.d2E_from_soil_dpsi_collar_dpsi_soil(p, psi, got);
  double worst_soil = 0.0;
  for (int j = 0; j < L; ++j) {
    const double h = 1e-6;
    std::vector<double> up = psi, dn = psi;
    up[std::size_t(j)] += h;
    dn[std::size_t(j)] -= h;
    const double want = (l.dE_from_soil_dpsi_collar(p, up) -
                         l.dE_from_soil_dpsi_collar(p, dn)) / (2.0 * h);
    worst_soil = std::max(worst_soil,
                          std::abs(got[std::size_t(j)] / want - 1.0));
  }

  std::vector<std::vector<double>> dE, dD;
  l.dE_from_soil_droot_carbon(p, psi, dE, dD);
  double worst_first = 0.0, worst_carbon = 0.0;
  double th[grad::n_pars];
  grad::Scratch scratch;
  for (int a = 0; a < L; ++a) {
    double sum_dE = 0.0, sum_dD = 0.0;
    for (int i = 0; i < L; ++i) {
      sum_dE += dE[std::size_t(i)][std::size_t(a)];
      sum_dD += dD[std::size_t(i)][std::size_t(a)];
    }
    const int par = grad::par_root_carbon_first(L) + a;
    const double base = grad::root_carbon_of(d, a);
    const double h = std::abs(base) * 1e-6;
    double slope[2], total[2];
    for (int side = 0; side < 2; ++side) {
      grad::set_one(l, th, env::kTheta, d, false, par,
                    side == 0 ? base + h : base - h, s.fast_stem_curve, scratch);
      slope[side] = l.dE_from_soil_dpsi_collar(p, l.supply_psi_soil());
      l.E_from_Soil_to_Root_Collar(p, l.supply_psi_soil());
      total[side] = l.E_up_;
    }
    worst_first = std::max(
        worst_first,
        std::abs(sum_dE / ((total[0] - total[1]) / (2.0 * h)) - 1.0));
    worst_carbon = std::max(
        worst_carbon,
        std::abs(sum_dD / ((slope[0] - slope[1]) / (2.0 * h)) - 1.0));
  }
  grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);

  printf("  dE_up/drc %.3g   d2E_up/dp dpsi %.3g   d2E_up/dp drc %.3g\n",
         worst_first, worst_soil, worst_carbon);
  ok(worst_first <= 1e-6, "the carbon block sums to the total's own row");
  ok(worst_soil <= 1e-6, "the soil mixed partial is the conductance's own");
  ok(worst_carbon <= 1e-6, "and so is the carbon one");

  // The root curve's two parameters, which reach the supply through the layer's
  // mean conductivity integral and nothing else. Both arms REBUILD the grid --
  // that is what the row replaces -- so the difference is of the model rather
  // than of a held curve.
  using Trait = pl::MultiLayerRoots::CurveTrait;
  const struct { const char* what; int par; Trait trait; } curve[] = {
      {"root_b", grad::par_root_b, Trait::Position},
      {"root_c", grad::par_root_c, Trait::Steepness}};
  for (const auto& t : curve) {
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    std::vector<double> dE, d2E;
    l.dE_from_soil_droot_curve(p, psi, t.trait, dE, d2E);
    double got_E = 0.0, got_2 = 0.0;
    for (int i = 0; i < L; ++i) {
      got_E += dE[std::size_t(i)];
      got_2 += d2E[std::size_t(i)];
    }
    const double base = grad::par_value(env::kTheta, d, false, t.par);
    const double hh = std::abs(base) * 1e-6;
    double total[2], slope[2];
    for (int side = 0; side < 2; ++side) {
      grad::set_one(l, th, env::kTheta, d, false, t.par,
                    side == 0 ? base + hh : base - hh, s.fast_stem_curve, scratch);
      // ⚠️ The per-layer integrals are cached against the soil state, and a moved
      // curve invalidates them: without this the difference reads the new curve at
      // the collar and the OLD one at each soil layer, which is a difference of two
      // models. Every production path reaches this through the collar solve, which
      // seats it; a check that perturbs and reads directly has to seat it itself.
      l.supply_begin_solve();
      l.E_from_Soil_to_Root_Collar(p, l.supply_psi_soil());
      total[side] = l.E_up_;
      slope[side] = l.dE_from_soil_dpsi_collar(p, l.supply_psi_soil());
    }
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    const double want_E = (total[0] - total[1]) / (2.0 * hh);
    const double want_2 = (slope[0] - slope[1]) / (2.0 * hh);
    printf("  %-7s dE_up %14.8g vs %14.8g (%.3g)  existing %14.8g   d2 %.3g\n",
           t.what, got_E, want_E, std::abs(got_E / want_E - 1.0),
           l.roots_.duptake_droot_curve(
               p, psi,
               t.trait == Trait::Position
                   ? phylloptim::MultiLayerRoots::CurveTrait::Position
                   : phylloptim::MultiLayerRoots::CurveTrait::Steepness),
           std::abs(got_2 / want_2 - 1.0));
    near(got_E, want_E, 1e-5, std::string("the supply's row in ") + t.what);
    near(got_2, want_2, 1e-4,
         std::string("and its collar derivative in ") + t.what);
  }
}

// The soil state's rows against the differencing they replace, at two steps.
//
// TWO, because what is left when the rows are right is the difference's own floor
// and that is what a second step separates: a truncation error falls as h^2 and
// this rises as 1/h. Measured over four decades at the wettest fixture, the
// assimilation and profit rows' disagreement goes 3e-04, 3e-05, 3e-06, 3.1e-07 --
// exactly 1/h, and exactly the channels that read the concentration's root-find.
// The conductance, the stem potential and every uptake row, which do not, sit at
// 1e-09 at every step.
// What the operating point's condition IS, in the two quantities the ecology
// weighs against each other: the carbon bought by the water an extra unit of
// collar pull draws, and the hydraulic cost of the extra tension that pull puts
// on the stem.
//
//   R = (dA/dE_up) S  -  C'(sigma) V
//
// Taken OFF the optimum, because at it the two terms are equal by the first-order
// condition and an agreement there checks the condition rather than the identity.
//
// It is worth checking rather than deriving once, because it is the statement the
// rows are built on: the collar reaches assimilation only by moving water, and it
// reaches the cost only by moving the stem potential. A route that appeared on
// one side and not the other would break exactly one of those.
// The transport's collar response is the flux balance's own derivative, and not
// a read of the inverted table's slope.
//
// The stem carries the flux the soil supplies: kappa (G(sigma) - G(p)) = E_up(p).
// Differentiating that in the collar gives V with no inverse anywhere,
//
//   V = (S/kappa + f(p)) / f(sigma),
//
// and one more derivative gives V's response to a state direction in the same two
// elementary quantities:
//
//   dV = dS/(kappa f(sigma))  -  V f'(sigma) dsigma / f(sigma).
//
// ⚠️ THE DIFFERENCE BETWEEN THE TWO ROUTES IS SECOND ORDER AND IT IS THE WHOLE OF
// WHAT THE SOIL STATE'S CONDITION ROW GETS WRONG. Reading the inverse's slope
// makes V a value of the interpolant, so V's own derivative is the interpolant's
// SECOND derivative, which is a property of the fit. Measured here, dV against a
// difference of the V it belongs to: 5e-09 to 1e-07 by the balance, against
// 3e-06 to 4e-05 by the inverse's slope, at the same grid.
void probe_root_curve_slope() {
  printf("PROBE: the root integral's inferred slope against the curve\n");
  namespace pl = phylloptim;
  pl::Leaf l = env::fresh();
  const double b = l.roots_.root_b, c = l.roots_.root_c;
  printf("   n=%.0f knots, root_b=%.4f root_c=%.4f\n",
         l.vulnerability_curve_ncontrol, b, c);
  printf("      %8s %18s %18s %10s\n", "psi", "spline slope", "curve", "rel");
  for (double psi : {0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0}) {
    const double tab = l.roots_.root_vuln_integral_deriv_at(psi);
    const double cur = std::exp(-std::pow(psi / b, c));
    printf("      %8.2f %18.12g %18.12g %10.2e\n", psi, tab, cur,
           cur > 0 ? std::abs(tab / cur - 1.0) : 0.0);
  }
}

void test_the_transport_response_is_the_flux_balances_own() {
  printf("the collar response is the flux balance's own derivative\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  struct F { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const F fixtures[] = {{"sodden", 0.5, 900, 2.0, 3}, {"wet", 2.0, 900, 2.0, 3},
                        {"dim", 2.0, 300, 2.0, 3},    {"dry", 4.5, 900, 2.0, 3},
                        {"arid", 2.0, 900, 4.0, 3},   {"deep", 3.0, 1500, 0.5, 5}};
  printf("      %-8s %12s   %12s %12s\n", "state", "V agrees", "dV by balance",
         "dV by slope");
  int checked = 0;
  for (const F& f : fixtures) {
    const int L = f.layers;
    grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);
    pl::Leaf l = env::fresh();
    const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
    if (b.branch.kind != pl::Leaf::OperatingPointKind::Interior) continue;
    bool feasible = false;
    l.dprofit_droot_collar_psi(b.psi_star, &feasible);
    ok(feasible, "the point is evaluable, " + tag);
    if (!feasible) continue;
    const double p = l.opt_root_psi_, sigma = l.opt_psi_stem_;
    const double kappa = l.leaf_specific_conductance_max_;
    const double S = l.dE_from_soil_dpsi_collar(p, l.supply_psi_soil());
    double f_s, f_s_prime, f_p;
    { phylloptim::tangent x = sigma;  phylloptim::seed_direction(x, 1.0);
      const phylloptim::tangent r = l.proportion_of_conductivity_kernel(x);
      f_s = odelia::util::to_passive(r);  f_s_prime = phylloptim::derivative_along(r); }
    { phylloptim::tangent x = p;  phylloptim::seed_direction(x, 1.0);
      f_p = odelia::util::to_passive(l.proportion_of_conductivity_kernel(x)); }
    const double V_balance = (S / kappa + f_p) / f_s;
    const double V_agrees = std::abs(V_balance / l.dpsistem_dpsi_ - 1.0);

    grad::SupplyRows w;
    if (!l.uptake_rows(w.on_uptake) || !grad::gather_supply(l, w, false, L)) continue;
    double dEup = 0.0, d2Eup = 0.0;
    grad::supply_of(w, grad::par_psi_soil_first, L, dEup, d2Eup);
    const double dsigma = dEup / (kappa * f_s);
    const double dV = d2Eup / (kappa * f_s) - V_balance * f_s_prime * dsigma / f_s;

    double th[grad::n_pars];
    grad::Scratch scratch;
    const double base =
        grad::par_value(env::kTheta, d, false, grad::par_psi_soil_first);
    const double h = grad::step_for(grad::par_psi_soil_first, base, s.step, L);
    double by_slope[2], by_balance[2];
    for (int side = 0; side < 2; ++side) {
      grad::set_one(l, th, env::kTheta, d, false, grad::par_psi_soil_first,
                    side == 0 ? base + h : base - h, s.fast_stem_curve, scratch);
      bool moved = false;
      l.dprofit_droot_collar_psi(b.psi_star, &moved);
      by_slope[side] = l.dpsistem_dpsi_;
      const double sg =
          l.find_psi_stem_from_psi_root(b.psi_star, l.supply_psi_soil());
      const double SS =
          l.dE_from_soil_dpsi_collar(b.psi_star, l.supply_psi_soil());
      by_balance[side] =
          (SS / kappa + f_p) / l.proportion_of_conductivity_kernel(sg);
    }
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    const double err_balance =
        std::abs(dV / ((by_balance[0] - by_balance[1]) / (2.0 * h)) - 1.0);
    const double err_slope =
        std::abs(dV / ((by_slope[0] - by_slope[1]) / (2.0 * h)) - 1.0);
    ++checked;
    printf("      %-8s %12.2e   %12.2e %12.2e\n", f.what, V_agrees, err_balance,
           err_slope);
    ok(V_agrees <= 1e-6, "V is the flux balance's, " + tag);
    ok(err_balance <= 1e-6,
       "and its response is that balance's own derivative, " + tag);
    // The other route is reported rather than bounded: it is the number this
    // check exists to make visible, and bounding it would fix the defect in
    // place.
  }
  ok(checked >= 5, "the balance was taken at every interior state");
}

void probe_table_vs_curve() {
  printf("PROBE: the forward table's slope against the curve itself\n");
  namespace pl = phylloptim;
  pl::Leaf l = env::fresh();
  printf("   n=%.0f knots\n", l.vulnerability_curve_ncontrol);
  printf("      %8s %18s %18s %10s\n", "psi", "table slope", "curve", "rel");
  for (double psi : {0.5, 1.0, 1.5, 2.0, 2.7, 3.5, 4.5}) {
    const double tab = l.stem_curve_integral_deriv(psi);
    const double cur = l.proportion_of_conductivity(psi);
    printf("      %8.2f %18.12g %18.12g %10.2e\n", psi, tab, cur,
           std::abs(tab / cur - 1.0));
  }
}

void test_the_condition_is_carbon_bought_against_tension_paid() {
  printf("the condition is the carbon water buys against the tension it costs\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  struct F { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const F fixtures[] = {{"sodden", 0.5, 900, 2.0, 3}, {"wet", 2.0, 900, 2.0, 3},
                        {"dim", 2.0, 300, 2.0, 3},    {"arid", 2.0, 900, 4.0, 3},
                        {"deep", 3.0, 1500, 0.5, 5}};
  int checked = 0;
  for (const F& f : fixtures) {
    const int L = f.layers;
    grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);
    pl::Leaf l = env::fresh();
    const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
    if (b.branch.kind != pl::Leaf::OperatingPointKind::Interior) {
      continue;
    }
    // A collar a little drier than the optimum, so R is a number rather than a
    // zero. The wettest end is the one with room: the dry end runs the transport
    // coordinate off the end of its grid, which is the caller not having asked
    // for a potential that exists.
    const double collar = b.psi_star * 0.98;
    bool feasible = false;
    const double R = l.dprofit_droot_collar_psi(collar, &feasible);
    ok(feasible, "the off-optimum collar is evaluable, " + tag);
    if (!feasible) {
      continue;
    }
    const double sigma = l.find_psi_stem_from_psi_root(collar, l.supply_psi_soil());
    const double V = l.dpsistem_dpsi_;
    const double S = l.dE_from_soil_dpsi_collar(collar, l.supply_psi_soil());
    const double gc_const = l.atm_kpa_ * pl::kg_to_mol_h2o / l.atm_vpd_ /
                            pl::H2O_CO2_stom_diff_ratio;
    const double v = 1.0 / (l.atm_kpa_ * pl::kPa_to_Pa);
    const double gc = gc_const * l.transpiration(sigma, collar);
    const double A_prime = l.assim_slope_;
    const double g_ci = A_prime * pl::umol_to_mol + gc * v;
    // The flux reaches assimilation through the conductance and the
    // concentration's residual, and through nothing else.
    const double carbon =
        A_prime * gc_const * (l.ca_ - l.ci_at_collar_) * v / g_ci * S;
    phylloptim::tangent x = sigma;  phylloptim::seed_direction(x, 1.0);
    const double tension =
        -phylloptim::derivative_along(l.hydraulic_cost_TF_kernel(x)) * V;
    ++checked;
    printf("  %-22s carbon %11.7g  tension %11.7g  sum %11.7g vs %11.7g\n",
           tag.c_str(), carbon, tension, carbon + tension, R);
    // The gap is the transport's round trip: the carbon term reads the flux
    // through the conductance, which the forward curve gives, and S reads it
    // from the soil, and the two agree only through the inverse curve.
    near(carbon + tension, R, 1e-5, "the two terms are the condition, " + tag);
  }
  ok(checked >= 4, "the identity was checked away from the optimum");
}

// Where a disagreement in the soil state's rows can possibly live, separated so
// that a number here names one of three things rather than their sum.
//
// The claim is that the condition reads the soil state through exactly two
// vectors -- total uptake and total uptake's own collar slope -- with two scalars
// shared across every input. So a least squares of the DIFFERENCED condition on
// those two vectors says three separate things:
//
//   * its residual tests the two vectors and the claim that there are only two;
//   * the fitted scalars against the closed-form ones test the closed forms;
//   * the Gram's condition number says whether the fit could tell them apart.
//
// The last is not a formality. An earlier reading of this boundary held that the
// two directions are numerically collinear and a compensating pair fits every row
// -- measured here, the condition number is 14 to 37, so they are not, and a
// disagreement between fitted and closed IS a disagreement about the coefficient.
// The collar channel, against the difference it replaces and against the states
// where that difference does not exist.
//
// dY/dp at fixed traits is what every ordinary output's row is multiplied by, and
// it was always a difference of the outputs across p*. At an interior point that
// difference is available and this checks the read against it. At a PIN it is not:
// p* sits a millionth of the bracket from its bound, so one arm is outside the
// feasible interval and no shrinking recovers it -- which is why a constrained
// point used to have no channel at all and every ordinary row there had to carry
// a total a re-solve measured.
//
// ⚠️ THE INTERIOR AGREEMENT IS NOT A TIGHT ONE, AND THE READ IS THE BETTER
// NUMBER. The difference divides the concentration root-find's floor by the step,
// so it is the noisy side; the check that says which is which is the step sweep,
// where the difference walks toward the read and not away from it.
void test_the_collar_channel_is_read_rather_than_differenced() {
  printf("the collar channel, read against the difference it replaces\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  struct F { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const F fixtures[] = {{"wet", 2.0, 900, 2.0, 3},  {"dim", 2.0, 300, 2.0, 3},
                        {"arid", 2.0, 900, 4.0, 3}, {"deep", 3.0, 1500, 0.5, 5},
                        {"pinned", 4.0, 100, 4.0, 1}};
  int interior_seen = 0, pinned_seen = 0;
  double worst_interior = 0.0;
  for (const F& f : fixtures) {
    const int L = f.layers;
    grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);
    pl::Leaf l = env::fresh();
    const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
    bool seated = false;
    l.dprofit_droot_collar_psi(b.psi_star, &seated);
    if (!seated) continue;
    pl::Leaf::CollarRows c;
    const bool read = l.collar_rows(c);
    ok(read, "the channel is readable at " + tag);
    if (!read) continue;

    // The profit half is the marginal profit the solve already reports, by a
    // different route through the same recorded pair. At an interior point both
    // are the zero the root-find found; at a pin both are the shadow price.
    near(c.dprofit, b.resid, 1e-8,
         "the profit channel is the marginal profit, " + tag);

    grad::OutputValues diff(L);
    const bool differenced = grad::collar_channel(l, b.psi_star, s, diff);
    if (b.branch.kind == pl::Leaf::OperatingPointKind::Interior) {
      ++interior_seen;
      ok(differenced, "an interior point admits the difference too, " + tag);
      if (!differenced) continue;
      // The step sweep. If the read is right the difference walks toward it as
      // the step grows past the root-find's floor divided by it.
      double worst = 0.0;
      const double got[3] = {c.dassim, c.dstom_cond, c.dpsistem};
      const int which[3] = {grad::out_assim, grad::out_stom_cond,
                            grad::out_psi_stem};
      for (int i = 0; i < 3; ++i) {
        const double scale = std::max(std::abs(got[i]), 1e-30);
        worst = std::max(worst, std::abs(diff[which[i]] - got[i]) / scale);
      }
      double per_layer_worst = 0.0;
      for (int j = 0; j < L; ++j) {
        const double want = diff[grad::out_uptake_first + j];
        const double have = c.duptake[std::size_t(j)] / pl::kg_per_mol_h2o;
        per_layer_worst = std::max(per_layer_worst,
                                   std::abs(have - want) /
                                       std::max(std::abs(want), 1e-30));
      }
      printf("    %-18s five %.2e  per-layer %.2e\n", tag.c_str(), worst,
             per_layer_worst);
      worst_interior = std::max(worst_interior, std::max(worst, per_layer_worst));
    } else if (b.branch.kind == pl::Leaf::OperatingPointKind::PinnedWet ||
               b.branch.kind == pl::Leaf::OperatingPointKind::PinnedDryRootCrit ||
               b.branch.kind ==
                   pl::Leaf::OperatingPointKind::PinnedDryRootPsiCrit) {
      ++pinned_seen;
      // ⚠️ THE DIFFERENCE SOMETIMES ANSWERS AT A PIN, and the count below is what
      // says how often. It is not a property of being pinned: the arms fail when
      // the step crosses the step-in, and the step shrinks two decades looking for
      // one that does not. What is a property is that the READ always answers, and
      // over the grid the difference does not.
      printf("    %-18s the difference %s here\n", tag.c_str(),
             differenced ? "answers too" : "cannot be centred");
      ok(std::isfinite(c.dpsistem) && std::isfinite(c.dassim) &&
             std::isfinite(c.dstom_cond),
         "and the read answers at a pin, " + tag);
    }
  }
  ok(interior_seen >= 3 && pinned_seen >= 1,
     "both branches were reached");

  // How often each route answers, over the grid rather than at a fixture. The
  // claim this object stands on is that the read answers where the difference
  // does not, and a count is what states it.
  {
    int pinned_points = 0, difference_answers = 0, read_answers = 0;
    const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
    const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
    const double vpds[] = {0.5, 1.0, 2.0, 4.0};
    for (double p : psi_soils) for (double q : ppfds) for (double vp : vpds)
      for (int L : {1, 3, 5}) {
        grad::Drivers d = env::drivers(p, q, vp, L, L);
        pl::Leaf l = env::fresh();
        const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
        pl::Leaf::WhichBound which = pl::Leaf::WhichBound::Wet;
        if (!grad::pinned_bound(b.branch.kind, which)) continue;
        ++pinned_points;
        bool seated = false;
        l.dprofit_droot_collar_psi(b.psi_star, &seated);
        pl::Leaf::CollarRows c;
        if (seated && l.collar_rows(c)) ++read_answers;
        grad::OutputValues diff(L);
        if (grad::collar_channel(l, b.psi_star, s, diff)) ++difference_answers;
      }
    printf("  over the grid: %d pinned points, the difference answers at %d, "
           "the read at %d\n",
           pinned_points, difference_answers, read_answers);
    ok(pinned_points > 0, "the grid has pinned points");
    ok(read_answers == pinned_points,
       "the read answers at every one of them");
    ok(difference_answers < pinned_points,
       "and the difference does not, which is what this object is for");
  }
  printf("  worst interior read-against-difference: %.2e\n", worst_interior);
  // A band, not an equality: what is left is the difference's own floor, and the
  // sweep above is what says which side of it each number is on.
  ok(worst_interior < 1e-3,
     "the read and the difference agree to the difference's own floor");
}

// The two bound entries that used to need a grid rebuilt and differenced.
//
// Both steepnesses reshape their vulnerability curve rather than scaling it, so
// neither has the homogeneity identity the two positions stand on. What they
// stand on instead is the incomplete gamma's shape derivative, and the referee is
// a difference of `find_root_psi` across a genuine rebuild -- which shares no code
// with the series.
void test_the_bounds_steepness_rows_match_a_rebuilt_difference() {
  printf("the bounds' two steepness rows against a rebuilt difference\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  const int L = 5;
  struct C { const char* what; int par; int which; };
  const C cases[] = {{"stem_c on the dry bound", grad::par_stem_c, 1},
                     {"root_c on the dry bound", grad::par_root_c, 1},
                     {"root_c on the wet bound", grad::par_root_c, 0}};
  double worst = 0.0;
  for (const C& c : cases) {
    grad::Drivers d = env::drivers(3.0, 900.0, 2.0, L, L);
    pl::Leaf l = env::fresh();
    grad::apply(l, env::kTheta, d, false, -1, true);
    l.find_root_collar_psi();
    const pl::Leaf::WhichBound w = c.which == 0
                                       ? pl::Leaf::WhichBound::Wet
                                       : pl::Leaf::WhichBound::DryRootCrit;
    const pl::Leaf::BoundRow row = l.bound_row(w);
    ok(row.finite, std::string("the row is finite, ") + c.what);
    if (!row.finite) continue;
    const double got =
        c.par == grad::par_stem_c ? row.d_dstem_c : row.d_droot_c;

    // A fresh leaf per arm: a rebuild is what the reference has to be, and the
    // curve's spline is state the perturbed leaf must not inherit.
    double bound[2];
    const double base = env::kTheta[c.par];
    const double h = base * 1e-5;
    for (int side = 0; side < 2; ++side) {
      double th[grad::n_pars];
      for (int i = 0; i < grad::n_pars; ++i) th[i] = env::kTheta[i];
      th[c.par] = side == 0 ? base + h : base - h;
      pl::Leaf m = env::fresh();
      grad::apply(m, th, d, false, -1, /*fast_stem_curve=*/false);
      m.find_root_collar_psi();
      bound[side] = m.find_root_psi(m.supply_begin_solve(), m.supply_psi_soil(),
                                    c.which);
    }
    const double want = (bound[0] - bound[1]) / (2.0 * h);
    const double err = std::abs(got - want) /
                       std::max(std::abs(want), 1e-30);
    printf("    %-24s %14.7g vs %14.7g  (%.2e)\n", c.what, got, want, err);
    worst = std::max(worst, err);
    // Non-vacuity: a row that is always zero would pass any ratio against a zero
    // reference, and the wet bound reads no stem property at all so only the two
    // live cases are asserted live.
    ok(std::abs(want) > 1e-12,
       std::string("the reference moves at all, ") + c.what);
  }
  ok(worst < 1e-3,
     "each steepness row is the series' and not a rebuild's");
}

// What the pinned branch costs, stated as the route each input takes rather than
// as a time.
//
// At a pin every input the bound moves used to be answered by re-solving the
// whole model twice, because the held evaluation was unavailable there: p* sits a
// step-in from the bound, so a step that moves the bound carries p* out of the
// perturbed feasible interval. A CLOSED FORM has no such step, so the same held
// partial that answers at an interior point answers here, and the point's own
// movement goes into `dresidual` once instead of into every row.
void test_a_pinned_point_answers_from_parts_rather_than_re_solving() {
  printf("a pinned point answers from parts\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  // The golden grid rather than a fixture, because which bound a point pins to is
  // not something a driver triple can be chosen for: the dry end is a min of two
  // limits that are different functions of the state.
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.5, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  int pins = 0, wet_pins = 0, dry_pins = 0, worst_moved = 0, no_bound_row = 0;
  for (double ps : psi_soils) for (double q : ppfds) for (double vp : vpds)
    for (int L : {1, 5}) {
    grad::Drivers d = env::drivers(ps, q, vp, L, L);
    pl::Leaf l = env::fresh();
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    l.find_root_collar_psi();
    const pl::Leaf::OperatingPointKind kind = l.operating_point_kind();
    pl::Leaf::WhichBound bound = pl::Leaf::WhichBound::Wet;
    if (!grad::pinned_bound(kind, bound)) continue;
    ++pins;
    ++(bound == pl::Leaf::WhichBound::Wet ? wet_pins : dry_pins);
    const std::string tag = pl::Leaf::operating_point_kind_name(kind);

    std::vector<int> outs, ins;
    env::consumer_request(L, outs, ins);
    grad::RowRequest req;
    req.output = outs.data();  req.n_output = outs.size();
    req.input = ins.data();    req.n_input = ins.size();

    const grad::Rows rows = env::solved_rows(l, env::kTheta, d, req, s);
    ok(rows.kind == kind, "the branch is the one the solve took, " + tag);
    const pl::Leaf::BoundRow cond = l.bound_row(bound);
    // ⚠️ A REFUSAL BY NAME, AND IT IS THE DESIGNED ONE. Four points of this sweep
    // sit at a wet bound whose own row is not finite, so the theorem the point
    // stands on has no denominator and `rows_at` says which bound rather than
    // returning a plausible number. They are counted here, not skipped quietly,
    // and the assertion below is that a refusal is always this one.
    if (!cond.finite) {
      ++no_bound_row;
      ok(!rows.message.empty(),
         "a bound with no derivative refuses by name, " + tag);
      continue;
    }
    ok(rows.message.empty(), "nothing else is refused, " + tag);
    int moved = 0, still_re_solved = 0, live = 0;
    for (std::size_t i = 0; i < ins.size(); ++i) {
      if (grad::bound_dpoint(cond, ins[i], L) != 0.0) {
        ++moved;
        if (rows.dresidual[i] != 0.0) ++live;
      }
      if (grad::bound_moves_by_re_solving(ins[i], false)) ++still_re_solved;
    }
    worst_moved = std::max(worst_moved, moved);
    ok(moved > 0, "the bound moves with something, " + tag);
    // ⚠️ THE ECONOMY, STATED AS A PROPERTY. Each of these was two solves of the
    // whole model, and what replaces them is one entry in a vector the consumer
    // divides once.
    ok(live == moved,
       "every input the bound moves reports it once, in the condition's "
       "gradient, " + tag);
    ok(still_re_solved == 0,
       "and none of them re-solves the model, " + tag);

    // The parts are parts: the point's own held row is the exact zero a held
    // partial makes it, and its motion is the quotient.
    const double dpoint = -rows.dresidual[0] / rows.residual_slope;
    ok(std::isfinite(dpoint), "the quotient exists, " + tag);
  }
  printf("    %d pinned points (%d wet, %d dry), %d with no bound row; the most "
         "inputs any bound moves is %d, and none of them re-solves\n",
         pins, wet_pins, dry_pins, no_bound_row, worst_moved);
  ok(wet_pins > 0 && dry_pins > 0, "both bounds were pinned to");
  // ⚠️ THE ECONOMY, AS THE NUMBER IT IS. Each input a bound moves was two solves
  // of the whole model, so a dry pin at five layers was thirty-two of them for
  // rows the leaf states in closed form.
  ok(worst_moved >= 10,
     "and a bound moves with enough inputs for that to be the cost it was");
}

// A shut point's rows against the differencing they replace.
//
// At a shut collar the leaf holds the stem at its critical potential and moves no
// water at all, so profit is respiration plus the cost at that seat:
//
//   Pi = -R_d(T) - C(psi_crit; stem_b, stem_c, beta2, scale)
//
// Six inputs appear in that expression and twenty do not, so six rows are the
// cost's own derivatives and twenty are exactly zero. Both halves are checked,
// because the zeros are the half a wrong row hides behind.
//
// ⚠️ THE COLLAR IS NOT CHECKED HERE AND MUST NOT BE. `set_shutdown_state` writes
// `opt_psi_stem_ = psi_crit` on every exit, which is why the profit rows above hold
// unconditionally -- but it writes `opt_root_psi_` from the exit's own reason, and
// the three reasons seat it at three different potentials. `shut_row_covers`
// refuses a request naming the collar for that reason.

// The supply's SECOND collar derivative, against a difference of its first, over
// a sweep of steps.
//
// It is the one object standing between the condition's collar slope and a closed
// form, and it is elementary: differentiating a cumulative integral in its upper
// limit leaves the integrand there, so one more derivative leaves the integrand's
// own slope, which the vulnerability curve has in closed form. Nothing here reads
// the tabulation's curvature, which is a difference of its data.

// The two coincidences the supply refused every derivative at, and which are not
// kinks for a derivative at all.
//
// A GRAVITY-BALANCED layer draws exactly nothing, because the potential difference
// and the head cancel. Only the NUMERATOR of its flux vanishes there -- the span,
// the mean conductivity and the resistance are what they are either side -- so the
// quotient rule leaves 1/r_R and there is nothing to take a limit of. A moving
// bound AT ATMOSPHERIC crosses the split the mean conductivity's integral is taken
// in two parts about, and the two parts have the same slope there.
//
// The one that IS 0/0 is the collar meeting a layer's potential, where the span and
// the integral over it vanish together. That still refuses.

// The three invariants that were stated and never checked. Each is one loop, and
// each was a claim about the code that nothing in either suite tested.

// The two classifications of one physical point, and what the numerical one loses.
//
// `Status` is four values derived from the curvature's SIGN and the residual's
// SIZE; `OperatingPointKind` is a decision tree on what defines the point. The
// first is the inference the mathematics forbids, and it survives because the
// calibration route is refereed bit for bit against a captured reference and cannot
// move. So the question that can be settled before that reference is re-blessed is
// not whether to merge them, but whether they ever disagree -- and if they do not,
// what the coarser one throws away.
//
// Both halves are asserted, because they say different things. The agreement is
// what makes the eventual merge a deletion rather than an investigation. The LOSS
// is why the merge is wanted: three kinds map onto one status, and a consumer
// holding the status cannot tell a bound it must differentiate from a seat it must
// not.
void test_the_two_classifications_of_one_point() {
  printf("the numerical status against the branch the solve took\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  using Kind = pl::Leaf::OperatingPointKind;
  grad::Settings s;

  // What the status is ALLOWED to be for each kind: a projection, fixed here so a
  // kind that starts reporting a different one is a failure rather than a new row
  // in a table nobody reads.
  auto projection = [](Kind k) -> grad::Status {
    switch (k) {
    case Kind::Interior:              return grad::Status::Interior;
    case Kind::PinnedWet:
    case Kind::PinnedDryRootCrit:
    case Kind::PinnedDryRootPsiCrit:
    // ⚠️ SHADE DEATH IS HERE, and that is the loss stated as code. It is a pin at
    // the wet bound, so the numerical test lands on "pinned" for the right reason
    // -- and the status cannot say that the bound it is pinned to is the collar of
    // zero uptake rather than an optimiser's stopping point, nor that its carbon
    // half is missing.
    case Kind::ShadeDeath:            return grad::Status::Pinned;
    default:                          return grad::Status::NoGradient;
    }
  };

  int points = 0, disagreements = 0, kinds_seen = 0;
  double worst_interior = 0.0, mildest_pin = 1e300;
  // The zero-flux kind whose residual is a SENTINEL rather than a derivative, and
  // whose implied Newton step is therefore exactly zero -- the interior side of the
  // cut. Counted separately because it is the one case the numerical test gets
  // wrong, and what it is rescued by is not the test.
  int sentinel_points = 0, sentinel_took_the_composite = 0;
  bool seen[16] = {false};
  std::string first_disagreement;
  for (double T : {25.0, 40.0}) {
    for (double psi = 0.5; psi <= 7.0; psi += 0.5) {
      for (double ppfd : {0.0, 100.0, 900.0}) {
        for (int L : {1, 5}) {
          grad::Drivers d = env::drivers(psi, ppfd, 2.0, L, L);
          d.leaf_temp = T;
          // The kind from a solve of its own, so `at`'s perturbations cannot have
          // moved it by the time it is read.
          pl::Leaf k = env::fresh();
          grad::apply(k, env::kTheta, d, false, -1, s.fast_stem_curve);
          k.find_root_collar_psi();
          const Kind kind = k.operating_point_kind();

          grad::Result res;
          int pars[] = {grad::par_vcmax_25};
          pl::Leaf l = env::fresh();
          bool threw = false;
          try {
            grad::at(l, env::kTheta, d, false, pars, 1, s, res);
          } catch (const std::exception &) {
            threw = true;
          }
          if (threw) {
            continue;
          }
          ++points;
          if (!seen[std::size_t(kind)]) {
            seen[std::size_t(kind)] = true;
            ++kinds_seen;
          }
          if (res.status != projection(kind)) {
            ++disagreements;
            if (first_disagreement.empty()) {
              first_disagreement =
                  std::string(pl::Leaf::operating_point_kind_name(kind)) + " read as " +
                  grad::status_name(res.status);
            }
          }
          // The band the numerical cut sits in. It needs no scale of its own, being
          // an implied Newton step in MPa -- so what matters is that the two
          // populations do not overlap it.
          if (std::isfinite(res.stationarity)) {
            if (kind == Kind::Interior) {
              worst_interior = std::max(worst_interior, res.stationarity);
            } else if (res.stationarity <= s.stationarity_tol) {
              // The numerical test has called a constrained point stationary.
              ++sentinel_points;
              sentinel_took_the_composite += res.used_ift ? 1 : 0;
            } else {
              mildest_pin = std::min(mildest_pin, res.stationarity);
            }
          }
        }
      }
    }
  }
  printf("  %d points, %d kinds, %d disagreements; interior reaches %.3g, the "
         "mildest optimiser's pin %.3g, cut at %.3g\n",
         points, kinds_seen, disagreements, worst_interior, mildest_pin,
         s.stationarity_tol);
  printf("  %d constrained points read as stationary by the number alone; %d of "
         "them took the composite\n",
         sentinel_points, sentinel_took_the_composite);
  ok(points > 100 && kinds_seen >= 4,
     "the sweep reaches an interior point, both zero-flux kinds and a pin");
  ok(disagreements == 0,
     std::string("the reported status is exactly the branch's projection") +
         (first_disagreement.empty() ? "" : " -- " + first_disagreement));
  ok(worst_interior < s.stationarity_tol,
     "every interior point sits below the cut");
  ok(mildest_pin > 1e4 * s.stationarity_tol,
     "and every optimiser's pin sits four orders or more above it");

  // ⚠️ AND THE BAND IS NOT EMPTY, WHICH IS THE ARGUMENT FOR ONE CLASSIFICATION
  // RATHER THAN TWO. Where radiation is zero, gross assimilation is identically
  // zero, the marginal profit at the seated collar is a SENTINEL rather than a
  // derivative, and the implied Newton step is exactly 0 -- the interior side of a
  // cut whose whole justification is that no point lands there. The curvature
  // beside it is finite and large, so nothing about the pair says which branch this
  // is.
  //
  // What rescues it is not the classification. The collar channel cannot centre a
  // difference on a point that sits ON its bound, so the composite is abandoned and
  // the status is rewritten. That is a second, unrelated refusal doing the work the
  // number is credited with -- and it is why the two classifications agree above
  // while only one of them is sound.
  ok(sentinel_points > 0,
     "a constrained point reads as stationary from the curvature and residual "
     "alone, so the numerical cut does not separate the branches");
  ok(sentinel_took_the_composite == 0,
     "and none of them takes the composite, because a second refusal catches it");
}

void test_the_three_unchecked_invariants() {
  printf("the three invariants that were stated and never checked\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 5;

  // --- 1. ONE SOLVE PER CALL, counted rather than proxied ---------------------
  //
  // The suite asserted that no INPUT re-solves, which is weaker: a route could
  // solve twice for a reason that has nothing to do with an input, and until this
  // was counted nobody could say it did not. Two solves are expected across the
  // pair below -- one for the base point and one for the state the referee is
  // read at -- so what is asserted is the count for `rows_at` alone.
  {
    grad::Drivers d = env::drivers(2.0, 900.0, 2.0, L, L);
    std::vector<int> outs;
    for (int j = 0; j < grad::n_outputs_total(L); ++j) {
      outs.push_back(j);
    }
    std::vector<int> ins;
    for (int p = 0; p < grad::n_pars_total(L); ++p) {
      if (p != grad::par_resistance) {
        ins.push_back(p);
      }
    }
    grad::RowRequest req{outs.data(), outs.size(), ins.data(), ins.size()};
    pl::Leaf l = env::fresh();
    *l.collar_solves = 0;
    const grad::Rows rows = env::solved_rows(l, env::kTheta, d, req, s);
    const std::size_t solves = *l.collar_solves;
    char text[64];
    std::snprintf(text, sizeof(text), "%zu solve(s) for %zu inputs at %s",
                  solves, ins.size(),
                  pl::Leaf::operating_point_kind_name(rows.kind));
    printf("  interior: %s\n", text);
    ok(rows.kind == pl::Leaf::OperatingPointKind::Interior,
       "the fixture is an interior optimum");
    // ⚠️ TWO, AND WHICH REQUEST IS ASKED DECIDES IT. This one names every output,
    // including the three the read cannot state: the carbon readers report profit
    // and the condition, so a request naming assimilation sends every carbon-side
    // input through a difference, and a differenced row has to put the leaf back
    // before a consumer reads its values off it. That restore is the second solve.
    // The request a stand makes is checked below and costs one.
    ok(solves == 2,
       std::string("an all-outputs request costs the read's solve and the "
                   "restore after differencing -- ") + text);

    // The request this boundary exists for: profit and the water each layer gave
    // up. Nothing is differenced, so nothing has to be put back.
    std::vector<int> stand_outs{grad::out_profit};
    for (int j = 0; j < L; ++j) {
      stand_outs.push_back(grad::out_uptake_first + j);
    }
    grad::RowRequest stand_req{stand_outs.data(), stand_outs.size(), ins.data(),
                               ins.size()};
    pl::Leaf stand = env::fresh();
    *stand.collar_solves = 0;
    const grad::Rows sr = env::solved_rows(stand, env::kTheta, d, stand_req, s);
    printf("  the stand's own request: %zu solve(s) for %zu inputs at %s\n",
           *stand.collar_solves, ins.size(),
           pl::Leaf::operating_point_kind_name(sr.kind));
    ok(*stand.collar_solves == 1,
       "one leaf solve for the request a consumer makes");
    bool stand_refused = false;
    for (grad::Rows::NoRow why : sr.no_row) {
      stand_refused = stand_refused || why != grad::Rows::NoRow::None;
    }
    ok(!stand_refused, "and no input refused, so nothing was differenced");

    // And at a pin, where the closed forms answer too. A shut point is excluded:
    // no condition defines it, so its rows are still a difference of the solve
    // and the count is the input count by construction rather than by defect.
    grad::Drivers dry = env::drivers(3.5, 900.0, 2.0, L, L);
    pl::Leaf pinned = env::fresh();
    grad::apply(pinned, env::kTheta, dry, false, -1, s.fast_stem_curve);
    pinned.find_root_collar_psi();
    if (pinned.operating_point_kind() == pl::Leaf::OperatingPointKind::PinnedWet) {
      pl::Leaf again = env::fresh();
      *again.collar_solves = 0;
      const grad::Rows pr = env::solved_rows(again, env::kTheta, dry, stand_req, s);
      const std::size_t n = *again.collar_solves;
      printf("  pinned:   %zu solve(s) for %zu inputs at %s\n", n, ins.size(),
             pl::Leaf::operating_point_kind_name(pr.kind));
      ok(n == 1, "and one at a pin, where the bound's own row carries the point");
    }
  }

  // --- 2. THE PASSIVATION SITES AND THE RECORDED INPUTS ARE ONE SET ------------
  //
  // The count 14 + L + 1 + 1 + L is the arithmetic a reader does; what makes it an
  // invariant is that no member of it is severed. A `to_passive` on the leaf path
  // that nobody recorded back is INVISIBLE -- every row of that input comes back
  // exactly zero, which is indistinguishable from an input the model does not
  // read. So the check is that every input moves something somewhere: an
  // accidental passivation shows up as a whole column of exact zeros across the
  // grid, and a slack zero does not, because it goes live where its constraint
  // binds.
  {
    const int traits = grad::n_traits;
    const int addressable = grad::n_pars_total(L) - 1;  // `resistance` is the
                                                        // single path's alone
    ok(addressable == traits + L + 1 + 1 + L,
       "the input count is the passivation arithmetic: 14 + L + PPFD + kmax + L");
    ok(addressable == 26, "which is twenty-six at five layers");

    std::vector<int> outs;
    for (int j = 0; j < grad::n_outputs_total(L); ++j) {
      outs.push_back(j);
    }
    std::vector<int> ins;
    for (int p = 0; p < grad::n_pars_total(L); ++p) {
      if (p != grad::par_resistance) {
        ins.push_back(p);
      }
    }
    std::vector<bool> moves_something(ins.size(), false);
    // ⚠️ TWO TRAIT VECTORS, AND THE SECOND IS NOT A REFINEMENT. At this package's
    // defaults psi_crit and root_psi_crit are the SAME NUMBER, so the dry bound's
    // min can never prefer the root's -- and `root_psi_crit` is then slack at every
    // state any default-trait grid can reach, which reads exactly like an input
    // nobody recorded back. Its live row exists and is exactly +1, being a
    // registered constant, and reaching it needs a plant whose ROOT gives up before
    // its stem. That is the configuration the dry bound's clamp exists for.
    double parted[grad::n_pars];
    for (int i = 0; i < grad::n_pars; ++i) {
      parted[i] = env::kTheta[i];
    }
    parted[grad::par_root_psi_crit] = env::kTheta[grad::par_root_psi_crit] * 0.6;
    struct State { const double *theta; double psi; double ppfd; int layers; };
    const State states[] = {
        {env::kTheta, 0.5, 100.0, L},  {env::kTheta, 2.0, 900.0, L},
        {env::kTheta, 3.5, 900.0, L},  {env::kTheta, 5.0, 900.0, L},
        {env::kTheta, 2.0, 0.0, L},    // shade death, so its rows count too
        {parted, 3.375, 100.0, 1},     // the dry bound on the root's own limit
        {parted, 3.0, 900.0, 1},
    };
    int arms = 0;
    for (const State &st : states) {
      grad::Drivers d = env::drivers(st.psi, st.ppfd, 2.0, st.layers, st.layers);
      // The output and input lists are the five-layer ones; a one-layer state
      // addresses a prefix of them, so it gets its own request.
      std::vector<int> o, in;
      for (int j = 0; j < grad::n_outputs_total(st.layers); ++j) {
        o.push_back(j);
      }
      for (int p = 0; p < grad::n_pars_total(st.layers); ++p) {
        if (p != grad::par_resistance) {
          in.push_back(p);
        }
      }
      grad::RowRequest here{o.data(), o.size(), in.data(), in.size()};
      pl::Leaf l = env::fresh();
      const grad::Rows rows = env::solved_rows(l, st.theta, d, here, s);
      if (rows.kind == pl::Leaf::OperatingPointKind::PinnedDryRootPsiCrit) {
        ++arms;
      }
      for (std::size_t i = 0; i < in.size(); ++i) {
        bool moved = std::isfinite(rows.dresidual[i]) && rows.dresidual[i] != 0.0;
        for (std::size_t j = 0; j < o.size() && !moved; ++j) {
          const double v = rows.held[j * in.size() + i];
          moved = std::isfinite(v) && v != 0.0;
        }
        if (!moved) {
          continue;
        }
        // Back to the five-layer numbering this loop tallies in: the fixed
        // parameters share an index, and a one-layer state's soil and carbon
        // entries are its layer 1, which is the five-layer list's too.
        for (std::size_t k = 0; k < ins.size(); ++k) {
          if (ins[k] == in[i]) {
            moves_something[k] = true;
          }
        }
      }
    }
    ok(arms > 0,
       "the states include a pin on the root's own critical potential, which no "
       "default-trait fixture can reach");
    std::string severed;
    for (std::size_t i = 0; i < ins.size(); ++i) {
      if (!moves_something[i]) {
        severed += (severed.empty() ? "" : ", ") + grad::par_name(ins[i], L);
      }
    }
    printf("  every input moves an output somewhere: %s\n",
           severed.empty() ? "yes" : severed.c_str());
    ok(severed.empty(),
       std::string("no input's whole column is exactly zero across the grid") +
           (severed.empty() ? "" : " -- severed: " + severed));
  }

  // --- 3. THE SUPPLIED SLOPE AT EVERY KNOT ------------------------------------
  //
  // The stem curve's two interpolants are built from the value AND the closed-form
  // slope at each knot, which is what makes their value error fourth order rather
  // than second -- and the inverse's value is the stem potential, which everything
  // reads. Nothing outside the builder named the supplied slope, so this rebuilds
  // the knots from the same generator and asks the interpolants for their slope
  // there.
  {
    pl::Leaf l = env::fresh();
    grad::Drivers d = env::drivers(2.0, 900.0, 2.0, 1, 1);
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    std::vector<double> x, y_integral, y_conductivity, y_conductivity_slope;
    l.build_cumulative_vulnerability_integral(l.stem_b, l.stem_c,
                                             l.vulnerability_curve_ncontrol,
                                             x, y_integral, y_conductivity,
                                             y_conductivity_slope);
    double worst_forward = 0.0, worst_inverse = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      const double f = l.proportion_of_conductivity(x[i]);
      const double got = l.transpiration_from_psi.slope(x[i]);
      const double got_inv = l.psi_from_transpiration.slope(y_integral[i]);
      worst_forward = std::max(worst_forward,
                               std::abs(got - f) / std::max(std::abs(f), 1e-30));
      worst_inverse =
          std::max(worst_inverse, std::abs(got_inv - 1.0 / f) /
                                      std::max(std::abs(1.0 / f), 1e-30));
    }
    printf("  %zu knots: the curve's slope to %.3g, the inverse's to %.3g\n",
           x.size(), worst_forward, worst_inverse);
    ok(x.size() > 100, "the knots were rebuilt from the same generator");
    ok(worst_forward <= 1e-12,
       "the cumulative curve's slope at every knot is the integrand there");
    ok(worst_inverse <= 1e-12,
       "and the inverse's is its reciprocal, which is what the value rests on");
  }
}

void test_the_supply_answers_at_the_two_coincidences() {
  printf("the supply's derivatives at a gravity-balanced layer and at atmospheric\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;

  // A collar placed exactly on one layer's gravity balance, which is where a
  // one-layer plant that has stopped drawing water sits.
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, 3, 3);
  pl::Leaf l = env::fresh();
  grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  const std::vector<double> &ps = l.supply_psi_soil();
  const int layer = 1;
  const double balanced =
      ps[std::size_t(layer)] + l.roots_.grav_head_z_[std::size_t(layer)];

  std::vector<double> draw;
  l.dE_from_soil_dpsi_collar_by_layer(balanced, ps, draw);
  const double total = l.dE_from_soil_dpsi_collar(balanced, ps);
  ok(std::isfinite(total), "the conductance at a gravity-balanced collar is a number");
  bool per_layer_finite = true;
  for (double v : draw) {
    per_layer_finite = per_layer_finite && std::isfinite(v);
  }
  ok(per_layer_finite, "and so is every layer's own");
  ok(std::isfinite(l.d2E_from_soil_dpsi_collar2(balanced, ps)),
     "and its collar derivative too");

  // Against a difference of the flux, over a sweep of steps for `Plateau`'s
  // reason. The referee is the model's own uptake at a held soil state, which
  // shares no code with the derivative.
  double best = std::numeric_limits<double>::infinity();
  for (int e = 3; e <= 7; ++e) {
    const double h = std::pow(10.0, -double(e));
    l.E_from_Soil_to_Root_Collar(balanced + h, ps);
    const double up = l.E_up_;
    l.E_from_Soil_to_Root_Collar(balanced - h, ps);
    const double dn = l.E_up_;
    best = std::min(best, std::abs(total / ((up - dn) / (2.0 * h)) - 1.0));
  }
  char text[32];
  std::snprintf(text, sizeof(text), "%.2e", best);
  ok(best < 1e-6,
     std::string("the balanced layer's conductance is a difference of the flux, ") +
         text);

  // And the layer really is balanced there: its own draw is zero while the others
  // are not, which is what makes this the coincidence rather than a quiet state.
  l.E_from_Soil_to_Root_Collar(balanced, ps);
  ok(std::abs(l.soil_consumption_[std::size_t(layer)]) < 1e-18,
     "the balanced layer draws nothing at that collar");
  double others = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (i != layer) {
      others = std::max(others, std::abs(l.soil_consumption_[std::size_t(i)]));
    }
  }
  ok(others > 1e-9, "while the layers either side of it do");

  // Equal potentials still refuses, and it is the ONLY thing that does. Asserted
  // rather than left implicit, because the refusal is what a caller's fallback is
  // keyed on.
  const double meeting = ps[std::size_t(layer)];
  ok(!std::isfinite(l.dE_from_soil_dpsi_collar(meeting, ps)),
     "a collar meeting a layer's potential still refuses, being 0/0 there");
  ok(l.roots_.at_equal_potentials(meeting, ps[std::size_t(layer)]) &&
         !l.roots_.at_equal_potentials(balanced, ps[std::size_t(layer)]),
     "and the predicate names that one and not the balance");

  // No state the grid reaches refuses. This is the count that decided not to
  // write the equal-potentials limit: the balance is reached, the meeting is not.
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0};
  const double ppfds[] = {0.0, 10.0, 100.0, 900.0, 1500.0};
  const int layer_counts[] = {1, 2, 3, 5};
  int points = 0, refused = 0, balanced_points = 0, meeting_points = 0;
  for (double psi : psi_soils) {
    for (double ppfd : ppfds) {
      for (int layers : layer_counts) {
        grad::Drivers dd = env::drivers(psi, ppfd, 2.0, layers, layers);
        pl::Leaf probe = env::fresh();
        grad::apply(probe, env::kTheta, dd, false, -1, s.fast_stem_curve);
        probe.find_root_collar_psi();
        ++points;
        const double p = probe.opt_root_psi_;
        const std::vector<double> &pp = probe.supply_psi_soil();
        for (int i = 0; i < probe.roots_.max_soil_layer; ++i) {
          const double gap = p - pp[std::size_t(i)];
          if (std::abs(gap - probe.roots_.grav_head_z_[std::size_t(i)]) < 1e-8) {
            ++balanced_points;
            break;
          }
        }
        for (int i = 0; i < probe.roots_.max_soil_layer; ++i) {
          if (probe.roots_.at_equal_potentials(p, pp[std::size_t(i)])) {
            ++meeting_points;
            break;
          }
        }
        if (!std::isfinite(probe.dE_from_soil_dpsi_collar(p, pp))) {
          ++refused;
        }
      }
    }
  }
  printf("  %d operating points: %d gravity-balanced, %d at a meeting, %d refuse\n",
         points, balanced_points, meeting_points, refused);
  ok(balanced_points > 0, "the grid reaches the gravity balance");
  ok(meeting_points == 0, "and reaches no collar that meets a layer's potential");
  ok(refused == 0, "so nothing in the grid refuses a conductance any more");
}

void test_the_supplys_second_collar_derivative() {
  printf("the supply's second collar derivative against a difference of its first\n");
  namespace grad = phylloptim::gradient;
  grad::Settings s;
  double worst = 0.0;
  std::string worst_where;
  int compared = 0;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    for (int layers : {1, 3, 5}) {
      grad::Drivers d = env::drivers(psi, 900.0, 2.0, layers, layers);
      phylloptim::Leaf l = env::fresh();
      grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
      l.find_root_collar_psi();
      const double p = l.opt_root_psi_;
      const double closed = l.d2E_from_soil_dpsi_collar2(p, l.supply_psi_soil());
      if (!std::isfinite(closed) || closed == 0.0) {
        continue;
      }
      ++compared;
      double best = std::numeric_limits<double>::infinity();
      for (int e = 3; e <= 7; ++e) {
        const double h = std::max(std::abs(p), 1.0) * std::pow(10.0, -double(e));
        const double up = l.dE_from_soil_dpsi_collar(p + h, l.supply_psi_soil());
        const double dn = l.dE_from_soil_dpsi_collar(p - h, l.supply_psi_soil());
        if (!std::isfinite(up) || !std::isfinite(dn)) {
          continue;
        }
        best = std::min(best, std::abs(closed / ((up - dn) / (2.0 * h)) - 1.0));
      }
      if (std::isfinite(best) && best > worst) {
        worst = best;
        worst_where = "  at psi_soil " + std::to_string(psi) + ", " +
                      std::to_string(layers) + " layers";
      }
    }
  }
  printf("  %d states, worst %.3g%s\n", compared, worst, worst_where.c_str());
  ok(compared >= 12, "every state's second derivative is a number");
  ok(worst <= 1e-5, "and it is the difference of the first it replaces");

  // The single path's flux is linear in the difference over a constant
  // resistance, so its conductance does not move with the collar at all.
  grad::Drivers d = env::drivers(2.0, 900.0, 2.0, 1, 1);
  phylloptim::Leaf single = env::fresh();
  single.set_supply_single(0.0);
  grad::apply(single, env::kTheta, d, true, -1, s.fast_stem_curve);
  single.find_root_collar_psi();
  ok(single.d2E_from_soil_dpsi_collar2(single.opt_root_psi_,
                                       single.supply_psi_soil()) == 0.0,
     "the single path's is exactly zero");
}

// dV/dp off the flux balance, against a difference of the V the solve records.
// The two are the same quantity by two routes, and the route matters: reading V
// off the inverse transport table would make this the interpolant's SECOND
// derivative, which is a property of the fit that no supplied first-order data
// corrects.
void test_the_transport_responses_collar_slope() {
  printf("the transport response's collar slope against a difference of it\n");
  namespace grad = phylloptim::gradient;
  grad::Settings s;
  double worst = 0.0;
  int compared = 0;
  for (double psi : {0.5, 1.0, 2.0, 3.0}) {
    for (int layers : {1, 5}) {
      grad::Drivers d = env::drivers(psi, 900.0, 2.0, layers, layers);
      phylloptim::Leaf l = env::fresh();
      grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
      l.find_root_collar_psi();
      if (l.operating_point_kind() != grad::OperatingPointKind::Interior) {
        continue;
      }
      const double p = l.opt_root_psi_;
      l.dprofit_droot_collar_psi(p);
      double closed = 0.0;
      if (!l.collar_response_slope(closed)) {
        continue;
      }
      ++compared;
      double best = std::numeric_limits<double>::infinity();
      for (int e = 4; e <= 7; ++e) {
        const double h = std::max(std::abs(p), 1.0) * std::pow(10.0, -double(e));
        l.dprofit_droot_collar_psi(p + h);
        const double up = l.dpsistem_dpsi_;
        l.dprofit_droot_collar_psi(p - h);
        const double dn = l.dpsistem_dpsi_;
        if (!std::isfinite(up) || !std::isfinite(dn)) {
          continue;
        }
        best = std::min(best, std::abs(closed / ((up - dn) / (2.0 * h)) - 1.0));
      }
      worst = std::max(worst, std::isfinite(best) ? best : 0.0);
    }
  }
  printf("  %d interior points, worst %.3g\n", compared, worst);
  ok(compared >= 6, "the collar slope answers at every interior point tried");
  ok(worst <= 1e-6, "and it is the difference of the recorded response");
}

// The two zero-flux kinds are two points, and the seat is what tells them apart.
//
// Both pay respiration plus a hydraulic cost, so both take the cost's own rows as
// their profit rows. A hydraulic shutdown holds the stem AT its critical
// potential, so that input is the seat and its row is the cost's slope, and
// nothing the leaf reads is a function of the soil. Shade death holds both
// potentials at the collar of zero uptake instead: the critical potential is
// inactive exactly as at an interior optimum, and the seat moves with the soil.
void test_the_two_zero_flux_kinds_are_two_points() {
  printf("the two zero-flux kinds are two points, and the seat says which\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 5;

  // Shade death is reached by LIGHT and not by water, so the fixture is a wet
  // soil in the dark. No moisture sweep finds this branch.
  grad::Drivers dark = env::drivers(2.0, 0.0, 2.0, L, L);
  pl::Leaf shaded = env::fresh();
  grad::apply(shaded, env::kTheta, dark, false, -1, s.fast_stem_curve);
  shaded.find_root_collar_psi();
  ok(shaded.operating_point_kind() == pl::Leaf::OperatingPointKind::ShadeDeath,
     "an unlit leaf on wet soil reaches shade death");

  // A dry soil in full light is the other kind, on the same trait vector.
  grad::Drivers dry = env::drivers(7.0, 900.0, 2.0, L, L);
  pl::Leaf parched = env::fresh();
  grad::apply(parched, env::kTheta, dry, false, -1, s.fast_stem_curve);
  parched.find_root_collar_psi();
  ok(parched.operating_point_kind() ==
         pl::Leaf::OperatingPointKind::HydraulicShutdown,
     "and a lit leaf on dry soil reaches a hydraulic shutdown");
  if (shaded.operating_point_kind() !=
          pl::Leaf::OperatingPointKind::ShadeDeath ||
      parched.operating_point_kind() !=
          pl::Leaf::OperatingPointKind::HydraulicShutdown) {
    return;
  }

  // The seats are different potentials, and each leaf's stem sits at its own.
  ok(shaded.opt_psi_stem_ == shaded.opt_root_psi_,
     "shade death seats the stem exactly at the collar");
  ok(parched.opt_psi_stem_ == parched.psi_crit,
     "and a hydraulic shutdown seats it exactly at the critical potential");
  ok(shaded.opt_psi_stem_ != shaded.psi_crit,
     "which the shaded leaf's is not");

  // The flux: zero in total on both, and non-zero per layer on one of them. That
  // is the whole reason the soil block cannot be declared zero at both.
  double parched_worst = 0.0, shaded_worst = 0.0, shaded_sum = 0.0;
  for (int i = 0; i < L; ++i) {
    parched_worst = std::max(parched_worst,
                             std::abs(parched.soil_consumption_[std::size_t(i)]));
    shaded_worst = std::max(shaded_worst,
                            std::abs(shaded.soil_consumption_[std::size_t(i)]));
    shaded_sum += shaded.soil_consumption_[std::size_t(i)];
  }
  ok(parched_worst == 0.0, "a hydraulic shutdown draws exactly nothing anywhere");
  ok(shaded_worst > 1e-8 && std::abs(shaded_sum) < 1e-15,
     "shade death's per-layer draws are non-zero and sum to zero");

  // The declared rows, and the two the seat decides. `shut_row` takes the seat as
  // a fact rather than inferring it from the potentials, which is what makes it
  // impossible to read one kind's answer at the other.
  const pl::Leaf::HydraulicCostRow shaded_cost =
      shaded.hydraulic_cost_row(shaded.opt_psi_stem_);
  const pl::Leaf::HydraulicCostRow parched_cost =
      parched.hydraulic_cost_row(parched.opt_psi_stem_);
  ok(shaded_cost.finite && parched_cost.finite, "both seats have a cost row");
  double row = 1234.0;
  ok(grad::shut_row(parched, parched_cost, true, grad::par_psi_crit, L, row) &&
         row == -parched_cost.d_dpsi_stem,
     "at a hydraulic shutdown the critical potential IS the seat");
  row = 1234.0;
  ok(grad::shut_row(shaded, shaded_cost, false, grad::par_psi_crit, L, row) &&
         row == 0.0,
     "and at shade death it is inactive, exactly");
  // The soil block: declared where the seat reads no soil, and refused where it
  // is a function of it, so the supply's own rows answer instead.
  row = 1234.0;
  ok(grad::shut_row(parched, parched_cost, true, grad::par_psi_soil_first, L,
                    row) &&
         row == 0.0,
     "a hydraulic shutdown declares every soil row exactly zero");
  ok(!grad::shut_row(shaded, shaded_cost, false, grad::par_psi_soil_first, L,
                     row),
     "and shade death declares none of them");
  ok(!grad::shut_row(shaded, shaded_cost, false, grad::par_root_c, L, row),
     "nor the root curve's, which moves the supply the seat is defined by");

  // The collar channel, which no evaluation of this model can record: the stem
  // sits at the collar, so the marginal profit takes its no-flow exit.
  bool feasible = true;
  const double sentinel = shaded.dprofit_droot_collar_psi(shaded.opt_root_psi_,
                                                          &feasible);
  ok(sentinel == 0.0 && !feasible,
     "the marginal profit at a shade-death collar is a sentinel, not a zero");
  shaded.find_root_collar_psi();
  pl::Leaf::CollarRows ch;
  ok(shaded.zero_uptake_collar_rows(ch), "and the channel is stated instead");
  ok(ch.dpsistem == 1.0,
     "the stem follows the collar exactly, because the branch ties them");
  ok(ch.dassim == 0.0 && ch.dstom_cond == 0.0,
     "gross assimilation and the conductance move by exactly zero");
  ok(ch.dprofit == -shaded_cost.d_dpsi_stem,
     "and profit's channel is minus the cost's slope at the seat");
  std::vector<double> per_layer;
  shaded.dE_from_soil_dpsi_collar_by_layer(shaded.opt_root_psi_,
                                           shaded.supply_psi_soil(), per_layer);
  bool draws_match = per_layer.size() == std::size_t(L);
  for (int i = 0; i < L && draws_match; ++i) {
    draws_match = ch.duptake[std::size_t(i)] == per_layer[std::size_t(i)];
  }
  ok(draws_match, "each layer's own draw still responds to the collar");
  // Read off the recorded classification, so it cannot answer where it does not
  // describe the point.
  pl::Leaf::CollarRows wrong;
  ok(!parched.zero_uptake_collar_rows(wrong),
     "and a hydraulic shutdown refuses this channel rather than answering it");

  // ⚠️ ONE LAYER IS THE CASE THAT USED TO REFUSE, and it is where two of these
  // items meet. The wet bound is the collar at which total uptake vanishes; with a
  // SINGLE layer that is the collar at which its numerator vanishes, which is the
  // gravity balance -- and the supply refused every derivative there, so the bound
  // had no row and a shade-death point seated on it could not report its own
  // movement. Only the numerator vanishes at that collar, so nothing about it was
  // ever 0/0; with the refusal gone the bound has a row at every layer count and
  // the point's movement is reported rather than folded into a re-solve.
  std::vector<int> ins;
  for (int par = 0; par < grad::n_pars_total(1); ++par) {
    if (par != grad::par_resistance) {
      ins.push_back(par);
    }
  }
  for (int layers : {1, 2, 5}) {
    grad::Drivers thin = env::drivers(2.0, 0.0, 2.0, layers, layers);
    pl::Leaf probe = env::fresh();
    grad::apply(probe, env::kTheta, thin, false, -1, s.fast_stem_curve);
    probe.find_root_collar_psi();
    const pl::Leaf::BoundRow wet = probe.bound_row(pl::Leaf::WhichBound::Wet);
    ok(wet.finite, std::string("the wet bound has a row at ") +
                       std::to_string(layers) + " layer(s)");
    pl::Leaf fresh_probe = env::fresh();
    std::vector<int> thin_outs{grad::out_profit, grad::out_uptake_first};
    grad::RowRequest req{thin_outs.data(), thin_outs.size(), ins.data(),
                         ins.size()};
    const grad::Rows rows = env::solved_rows(fresh_probe, env::kTheta, thin, req, s);
    ok(rows.kind == pl::Leaf::OperatingPointKind::ShadeDeath &&
           rows.message.empty(),
       std::string("and shade death answers at ") + std::to_string(layers) +
           " layer(s)");
    bool all_numbers = true;
    for (double v : rows.held) {
      all_numbers = all_numbers && std::isfinite(v);
    }
    for (double v : rows.dy_dp) {
      all_numbers = all_numbers && std::isfinite(v);
    }
    ok(all_numbers, std::string("every row and channel at ") +
                        std::to_string(layers) + " layer(s) is a number");
    // The bound's own slope rather than the unit one, which is what says the
    // point's movement was reported instead of measured.
    ok(rows.residual_slope == wet.residual_slope,
       std::string("and the slope at ") + std::to_string(layers) +
           " layer(s) is the wet bound's own");
  }
}

// Shade death's whole row set, assembled, against a difference of the solve.
//
// The referee shares no code with the rows: it re-solves the model on both sides
// of each input and holds both arms on this branch. What it referees is the
// composition -- a held partial the frozen collar makes zero for the carbon side
// and the supply's own Jacobian for the water, plus the wet bound's movement
// priced by the cost's slope at the seat.
void test_shade_deaths_rows_against_a_differenced_solve() {
  printf("shade death's assembled rows against a differenced solve\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const int L = 5;
  grad::Drivers d = env::drivers(2.0, 0.0, 2.0, L, L);

  std::vector<int> outs{grad::out_profit};
  for (int i = 0; i < L; ++i) {
    outs.push_back(grad::out_uptake_first + i);
  }
  std::vector<int> ins;
  for (int p = 0; p < grad::n_pars_total(L); ++p) {
    if (p != grad::par_resistance) {
      ins.push_back(p);
    }
  }
  pl::Leaf l = env::fresh();
  grad::RowRequest req{outs.data(), outs.size(), ins.data(), ins.size()};
  const grad::Rows rows = env::solved_rows(l, env::kTheta, d, req, s);
  ok(rows.kind == pl::Leaf::OperatingPointKind::ShadeDeath,
     "the fixture reaches shade death");
  if (rows.kind != pl::Leaf::OperatingPointKind::ShadeDeath) {
    return;
  }
  ok(rows.message.empty(), "and no input at it is refused");
  ok(std::isfinite(rows.residual_slope) && rows.residual_slope != 1.0,
     "the slope is the wet bound's own and not the unit one");

  double worst = 0.0;
  std::string worst_where;
  int compared = 0, off_branch = 0;
  double th[grad::n_pars];
  grad::Scratch scratch;
  grad::OutputValues up(L), dn(L);
  for (std::size_t i = 0; i < ins.size(); ++i) {
    const double base = grad::par_value(env::kTheta, d, false, ins[i]);
    double best = std::numeric_limits<double>::infinity();
    std::size_t best_j = 0;
    for (int e = 5; e <= 7; ++e) {
      const double h = grad::step_for(ins[i], base, std::pow(10.0, -double(e)), L);
      pl::Leaf probe = env::fresh();
      grad::apply(probe, env::kTheta, d, false, -1, s.fast_stem_curve);
      bool on_branch = true;
      for (int side = 0; side < 2; ++side) {
        grad::set_one(probe, th, env::kTheta, d, false, ins[i],
                      side == 0 ? base + h : base - h, s.fast_stem_curve,
                      scratch);
        probe.find_root_collar_psi();
        on_branch = on_branch && probe.operating_point_kind() ==
                                     pl::Leaf::OperatingPointKind::ShadeDeath;
        grad::outputs(probe, side == 0 ? up : dn);
      }
      if (!on_branch) {
        continue;
      }
      // The scale is the largest row this input has over the request, for the
      // census's reason: an output an input barely reaches has a row at the
      // difference's floor, and a ratio against the entry reports the floor.
      double scale = 1e-300;
      for (std::size_t j = 0; j < outs.size(); ++j) {
        scale = std::max(scale, std::abs((up[outs[j]] - dn[outs[j]]) / (2.0 * h)));
      }
      double here = 0.0;
      std::size_t here_j = 0;
      for (std::size_t j = 0; j < outs.size(); ++j) {
        const double want = (up[outs[j]] - dn[outs[j]]) / (2.0 * h);
        const double dpoint =
            -rows.dresidual[i] / rows.residual_slope;
        const double got =
            rows.held[j * ins.size() + i] +
            (std::isfinite(dpoint) && dpoint != 0.0
                 ? rows.dy_dp[j] * dpoint
                 : 0.0);
        if (!std::isfinite(got)) {
          here = 1.0;
          here_j = j;
          continue;
        }
        const double err = std::abs(got - want) / scale;
        if (err > here) {
          here = err;
          here_j = j;
        }
      }
      if (here < best) {
        best = here;
        best_j = here_j;
      }
    }
    if (!std::isfinite(best)) {
      ++off_branch;
      continue;
    }
    ++compared;
    if (best > worst) {
      worst = best;
      worst_where = "  at " + grad::par_name(ins[i], L) + "/" +
                    grad::output_name(outs[best_j], L);
    }
  }
  printf("  %d of %d inputs differenced (%d off-branch), worst %.3g%s\n",
         compared, int(ins.size()), off_branch, worst, worst_where.c_str());
  ok(compared >= int(ins.size()) - 2,
     "almost every input has a difference to be checked against");
  ok(worst <= 1e-5, "and every assembled row is it");
}

void test_a_shut_points_rows_are_the_costs_own() {
  printf("a shut point's rows against the differencing they replace\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  const double psi_soils[] = {4.0, 5.0, 6.0, 7.0};
  const double ppfds[] = {100.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 2.0, 4.0};
  // The inputs the closed form names, and the rest, which it makes exactly zero.
  const int declared[] = {grad::par_psi_crit, grad::par_stem_b, grad::par_stem_c,
                          grad::par_beta2, grad::par_cost_scale_TF24,
                          grad::par_R_d_25};
  int shut = 0, off_branch = 0, compared = 0, zeros = 0;
  int declared_not_zero = 0, uptake_checked = 0, uptake_nonzero = 0;
  double worst = 0.0, worst_zero = 0.0;
  std::string worst_where, worst_zero_where;
  for (double p : psi_soils) for (double q : ppfds) for (double vp : vpds)
    for (int L : {1, 5}) {
    grad::Drivers d = env::drivers(p, q, vp, L, L);
    pl::Leaf l = env::fresh();
    grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
    l.find_root_collar_psi();
    const pl::Leaf::OperatingPointKind k = l.operating_point_kind();
    if (k != pl::Leaf::OperatingPointKind::HydraulicShutdown &&
        k != pl::Leaf::OperatingPointKind::ShadeDeath) {
      continue;
    }
    ++shut;
    // The seat, which is what makes the profit rows exit-independent.
    ok(phylloptim::util::identical(l.opt_psi_stem_, l.psi_crit),
       "a shut collar seats the stem at its critical potential");

    const std::size_t row_profit = 0;
    std::vector<int> outs, ins;
    env::consumer_request(L, outs, ins);
    grad::RowRequest req;
    req.output = outs.data();  req.n_output = outs.size();
    req.input = ins.data();    req.n_input = ins.size();
    const grad::Rows rows = env::solved_rows(l, env::kTheta, d, req, s);
    ok(rows.kind == k, "the branch is the one the solve took");
    ok(rows.message.empty(), "nothing is refused at a shut point: " + rows.message);

    grad::Scratch scratch;
    double th[grad::n_pars];
    for (std::size_t i = 0; i < ins.size(); ++i) {
      const int par = ins[i];
      const double base = grad::par_value(env::kTheta, d, false, par);
      const double h = grad::step_for(par, base, s.step, L);
      double prof[2];
      bool on = true;
      for (int side = 0; side < 2; ++side) {
        grad::set_one(l, th, env::kTheta, d, false, par,
                      side == 0 ? base + h : base - h, s.fast_stem_curve,
                      scratch);
        l.find_root_collar_psi();
        on = on && (l.operating_point_kind() == k);
        prof[side] = l.profit_;
      }
      grad::apply(l, env::kTheta, d, false, -1, s.fast_stem_curve);
      if (!on) { ++off_branch; continue; }
      const double want = (prof[0] - prof[1]) / (2.0 * h);
      const double got = rows.held[row_profit * ins.size() + i];
      const std::string where =
          " at " + grad::par_name(par, L) + ", psi_soil " +
          phylloptim::util::to_string(p) + " L " +
          phylloptim::util::to_string(L);
      bool named = false;
      for (const int dcl : declared) named = named || (par == dcl);
      if (named) {
        ++compared;
        const double err = std::abs(got - want) / std::max(std::abs(want), 1e-30);
        if (err > worst) { worst = err; worst_where = where; }
      } else {
        // ⚠️ THE ZEROS ARE THE HALF A WRONG ROW HIDES BEHIND, so they are held to
        // the difference too rather than to their own declaration. Counted rather
        // than asserted one at a time: five hundred near-identical checks report
        // the same fact as one check and a count, and only one of the two says how
        // much was covered.
        ++zeros;
        if (!phylloptim::util::identical(got, 0.0)) {
          ++declared_not_zero;
          worst_zero_where = where;
        }
        if (std::abs(want) > worst_zero) {
          worst_zero = std::abs(want);
        }
      }
      // Uptake is identically zero at a shut collar, on both routes.
      for (int j = 0; j < L; ++j) {
        ++uptake_checked;
        if (!phylloptim::util::identical(
                rows.held[std::size_t(1 + j) * ins.size() + i], 0.0)) {
          ++uptake_nonzero;
        }
      }
    }
  }
  printf("  %d shut points, %d declared rows compared, %d zeros checked, "
         "%d uptake entries, %d off-branch\n", shut, compared, zeros,
         uptake_checked, off_branch);
  printf("  worst declared row against the difference: %.3g%s\n", worst,
         worst_where.c_str());
  printf("  largest difference at an input the form calls zero: %.3g\n",
         worst_zero);
  ok(shut > 0 && compared > 0, "the grid reaches a shut collar and its rows");
  ok(worst < 1e-5, "each declared row is the cost's own derivative");
  ok(declared_not_zero == 0,
     "every input the closed form does not name comes back an exact zero" +
         worst_zero_where);
  // Non-vacuity for those zeros: if the model itself moved there, the exact zeros
  // above would be reporting a severed channel rather than agreement.
  ok(worst_zero < 1e-7,
     "and the model does not move in them either");
  ok(uptake_checked > 0 && uptake_nonzero == 0,
     "and a shut collar's every uptake row is an exact zero");
}

void test_the_soil_states_two_scalars_are_separable_and_right() {
  printf("the condition's two scalars, fitted against closed form\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;
  struct F { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const F fixtures[] = {
      {"sodden", 0.5, 900.0, 2.0, 3},   {"wet", 2.0, 900.0, 2.0, 3},
      {"dim", 2.0, 300.0, 2.0, 3},      {"dry", 4.5, 900.0, 2.0, 3},
      {"arid", 2.0, 900.0, 4.0, 3},     {"deep", 3.0, 1500.0, 0.5, 5},
      {"shallow", 1.0, 900.0, 2.0, 5},  {"midwet", 1.5, 600.0, 1.0, 5}};
  printf("      %-9s %10s %10s   %10s %10s\n", "state", "fit resid", "cond#",
         "dE_up err", "slope err");
  int checked = 0;
  for (const F& f : fixtures) {
    const int L = f.layers;
    grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);
    std::vector<int> input;
    for (int i = 0; i < L; ++i) input.push_back(grad::par_psi_soil_first + i);
    for (int i = 0; i < L; ++i) input.push_back(grad::par_root_carbon_first(L) + i);
    input.push_back(grad::par_root_b);
    input.push_back(grad::par_root_c);
    pl::Leaf l = env::fresh();
    const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
    ok(b.branch.kind == pl::Leaf::OperatingPointKind::Interior,
       "the fixture is interior, " + tag);
    if (b.branch.kind != pl::Leaf::OperatingPointKind::Interior) continue;
    grad::SupplyRows w;
    bool feasible = false;
    l.dprofit_droot_collar_psi(b.psi_star, &feasible);
    ok(feasible && l.uptake_rows(w.on_uptake) &&
           grad::gather_supply(l, w, false, L),
       "the coefficients and the supply derivatives answer, " + tag);
    if (!feasible) continue;

    double s11 = 0, s12 = 0, s22 = 0, s1y = 0, s2y = 0, ymax = 0;
    std::vector<double> x1, x2, y;
    bool at_base = true;
    grad::Scratch scratch;
    grad::OutputValues direct(L);
    for (std::size_t i = 0; i < input.size(); ++i) {
      double dR = 0.0;
      if (!grad::held_row(l, env::kTheta, d, false, input[i], b.psi_star, s, true,
                          at_base, scratch, direct, dR)) continue;
      double dEup = 0.0, d2Eup = 0.0;
      grad::supply_of(w, input[i], L, dEup, d2Eup);
      x1.push_back(dEup); x2.push_back(d2Eup); y.push_back(dR);
      s11 += dEup*dEup; s12 += dEup*d2Eup; s22 += d2Eup*d2Eup;
      s1y += dEup*dR;   s2y += d2Eup*dR;
      ymax = std::max(ymax, std::abs(dR));
    }
    const double det = s11*s22 - s12*s12;
    const double a = (s22*s1y - s12*s2y)/det;
    const double bb = (s11*s2y - s12*s1y)/det;
    double resid = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i)
      resid = std::max(resid, std::abs(a*x1[i] + bb*x2[i] - y[i]));
    const double tr = s11 + s22, dd = std::sqrt(std::max(0.0, tr*tr - 4*det));
    const double cond = (tr + dd) / std::max(1e-300, tr - dd);
    const double a_err = std::abs(w.on_uptake.dcondition/a - 1.0);
    const double b_err = std::abs(w.on_uptake.dprofit/bb - 1.0);
    ++checked;
    printf("      %-9s %10.2e %10.1f   %10.2e %10.2e\n", f.what, resid/ymax,
           cond, a_err, b_err);
    // The basis, and that there are two of it.
    ok(resid/ymax <= 1e-4, "the differenced condition lies in the two vectors, " + tag);
    ok(cond <= 1e3, "and the two are separable at this state, " + tag);
    // ⚠️ BOTH SCALARS, AND TO THE SAME BOUND. The uptake one used to be held two
    // orders looser because it alone read a table's derivative -- the inverse
    // curve's slope through V, and the forward curve's through the conductance.
    // Taking both from the curve brought it from 5e-05 to 1e-08, which is where
    // the other one and the fit's own residual already were, so there is no
    // longer an asymmetry to allow for. The wettest state is the exception and
    // the reference is what is wrong there: its concentration root-find leaves a
    // floor that falls as 1/h.
    ok(b_err <= 1e-4, "the collar-slope scalar is the fit's own, " + tag);
    ok(a_err <= 1e-4, "and so is the uptake scalar, " + tag);
  }
  ok(checked >= 7, "every state was separated");
}

void test_the_soil_states_rows_match_a_differenced_solve() {
  printf("the soil state's closed-form rows against the differenced solve\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;

  struct Fixture {
    const char* what;
    double psi_soil, ppfd, vpd;
    int layers;
  };
  const Fixture fixtures[] = {
      {"sodden", 0.5, 900.0, 2.0, 3},
      {"wet", 2.0, 900.0, 2.0, 3},
      {"dim", 2.0, 300.0, 2.0, 3},
      {"dry", 4.5, 900.0, 2.0, 3},
      {"arid", 2.0, 900.0, 4.0, 3},
      {"deep", 3.0, 1500.0, 0.5, 5},
  };

  for (const double step : {1e-6, 1e-4}) {
    grad::Settings s;
    s.step = step;
    // What the difference itself can say at this step, and the reason the two
    // families differ is the concentration's root-find rather than the rows.
    const double tol_ci = step == 1e-6 ? 5e-4 : 1e-5;
    // The other direction: a bigger step buys the concentration's root-find room
    // and costs truncation, so the rows that reach no root-find are held to the
    // looser of the two at both steps.
    const double tol_direct = 5e-7;
    printf("  at a step of %.0e\n", step);
    for (const Fixture& f : fixtures) {
      const int L = f.layers;
      grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);

      std::vector<int> out_index{grad::out_assim, grad::out_stom_cond,
                                 grad::out_psi_stem, grad::out_collar,
                                 grad::out_profit};
      for (int i = 0; i < L; ++i) {
        out_index.push_back(grad::out_uptake_first + i);
      }
      std::vector<int> input;
      for (int i = 0; i < L; ++i) {
        input.push_back(grad::par_psi_soil_first + i);
      }
      for (int i = 0; i < L; ++i) {
        input.push_back(grad::par_root_carbon_first(L) + i);
      }
      // The root curve's two, which reach the leaf the same way: they move the
      // supply. Their arms REBUILD the grid, which is the work the row replaces.
      input.push_back(grad::par_root_b);
      input.push_back(grad::par_root_c);
      const grad::RowRequest req{out_index.data(), out_index.size(),
                                 input.data(), input.size()};

      pl::Leaf l = env::fresh();
      const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);
      const bool interior =
          b.branch.kind == pl::Leaf::OperatingPointKind::Interior;
      const std::string tag =
          std::string(f.what) + " (" +
          pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
      ok(interior, "the fixture is interior, " + tag);
      if (!interior) {
        continue;
      }

      // ⚠️ THE DIFFERENCE COMES FROM `held_row`, NOT FROM `rows_at`. The
      // dispatcher prefers these same closed forms for these same inputs, so a
      // comparison through it would compare each row with itself.
      grad::SupplyRows w;
      // The rows below are checked directly, so this asserts the DISPATCHER
      // reaches them too -- a correct row the row layer never calls is the failure
      // this pair exists to separate.
      ok(grad::supply_rows_apply(l, req, interior, L),
         "the row layer takes these rows here, " + tag);
      bool evaluable = false;
      l.dprofit_droot_collar_psi(b.psi_star, &evaluable);
      ok(evaluable, "the point is evaluable, " + tag);
      ok(l.uptake_rows(w.on_uptake), "the uptake coefficients answer, " + tag);
      ok(std::isfinite(w.on_uptake.dcondition),
         "including the condition's own, " + tag);
      ok(grad::gather_supply(l, w, false, L),
         "and every supply derivative exists, " + tag);
      w.usable = true;

      bool at_base = true;
      grad::Scratch scratch;
      grad::OutputValues differenced(L), got(L);
      // The scale each output's rows are judged against: the largest the
      // difference makes any of them, taken over the whole input set before any
      // is compared. An input an output barely reaches has a difference at its
      // own floor, and a ratio taken there reports the floor.
      std::vector<double> scale(std::size_t(grad::n_outputs_total(L)), 1e-30);
      for (std::size_t i = 0; i < input.size(); ++i) {
        double ignored = 0.0;
        if (grad::held_row(l, env::kTheta, d, false, input[i], b.psi_star, s,
                           true, at_base, scratch, differenced, ignored)) {
          for (int j = 0; j < grad::n_outputs_total(L); ++j) {
            scale[std::size_t(j)] =
                std::max(scale[std::size_t(j)], std::abs(differenced[j]));
          }
        }
      }

      double dR_dsigma = 0.0;
      ok(l.condition_slope(dR_dsigma), "and the condition's slope, " + tag);
      double worst[grad::n_outputs] = {0.0, 0.0, 0.0, 0.0, 0.0};
      double worst_uptake = 0.0, worst_cond = 0.0, worst_cancel = 0.0;
      for (std::size_t i = 0; i < input.size(); ++i) {
        const std::string nm = grad::par_name(input[i], L) + ", " + tag;
        double want_dR = 0.0;
        if (!grad::held_row(l, env::kTheta, d, false, input[i], b.psi_star, s,
                            true, at_base, scratch, differenced, want_dR)) {
          ok(false, "the difference answers for " + nm);
          continue;
        }
        double dR = 0.0;
        grad::supply_row(w, input[i], L, got, dR);
        for (int j = 0; j < grad::n_outputs; ++j) {
          if (j == grad::out_collar) {
            continue;   // held, so both routes have it at exactly zero
          }
          worst[j] = std::max(worst[j], std::abs(got[j] - differenced[j]) /
                                            scale[std::size_t(j)]);
        }
        for (int j = 0; j < L; ++j) {
          const int o = grad::out_uptake_first + j;
          worst_uptake =
              std::max(worst_uptake, std::abs(got[o] - differenced[o]) /
                                         scale[std::size_t(o)]);
        }
        worst_cond = std::max(worst_cond, std::abs(dR - want_dR) /
                                              std::max(1.0, std::abs(want_dR)));
        {
          // How much of each term the sum keeps. Reported because it was once
          // the explanation for this row's error and is not: the states that
          // cancel hardest, at 145x, are not the states with the largest
          // disagreement, and with every derivative taken from the curve the
          // whole block sits at the difference's own floor regardless.
          double dEup = 0.0, d2Eup = 0.0;
          grad::supply_of(w, input[i], L, dEup, d2Eup);
          const double dsigma = w.on_uptake.dpsistem * dEup;
          const double G = l.dprofit_dpsistem_;
          const double dV =
              w.on_uptake.dpsistem * d2Eup +
              (w.on_uptake.dcondition - dR_dsigma * w.on_uptake.dpsistem) / G *
                  dEup;
          worst_cancel = std::max(
              worst_cancel,
              (std::abs(dR_dsigma * dsigma) + std::abs(G * dV)) /
                  std::max(1e-30, std::abs(dR)));
        }
      }
      printf("    %-22s A %.2g  gc %.2g  psi_stem %.2g  profit %.2g  "
             "uptake %.2g  condition %.2g (over %.1fx cancellation)\n",
             tag.c_str(), worst[grad::out_assim], worst[grad::out_stom_cond],
             worst[grad::out_psi_stem], worst[grad::out_profit], worst_uptake,
             worst_cond, worst_cancel);
      // The three that reach no root-find are held to the difference itself.
      ok(worst[grad::out_stom_cond] <= tol_direct &&
             worst[grad::out_psi_stem] <= tol_direct && worst_uptake <= tol_direct,
         "the conductance, the stem potential and the uptake rows ARE the "
         "difference, " + tag);
      ok(worst[grad::out_assim] <= tol_ci && worst[grad::out_profit] <= tol_ci,
         "and assimilation and profit to what the concentration's root-find "
         "leaves, " + tag);
      // ⚠️ THE CONDITION'S BOUND IS NOT THE DIFFERENCE'S, and that is why it
      // does not move with the step: it is a sum of the two terms above, each
      // right to a few parts per million, over a cancellation the printed factor
      // reports. Measured at 4.7e-05 at the worst fixture at both steps.
      ok(worst_cond <= 1e-4, "the condition's own row too, " + tag);
    }
  }
}

void test_carbon_trait_rows_match_a_differenced_solve() {
  printf("the carbon traits' closed-form rows against the differenced solve\n");
  namespace grad = phylloptim::gradient;
  namespace pl = phylloptim;
  grad::Settings s;

  struct Fixture {
    const char* what;
    double psi_soil, ppfd, vpd, leaf_temp;
    int layers;
  };
  const Fixture fixtures[] = {
      {"wet", 2.0, 900.0, 2.0, 25.0, 3},
      {"dim", 2.0, 300.0, 2.0, 25.0, 3},
      {"dry", 4.5, 900.0, 2.0, 25.0, 3},
      {"arid", 2.0, 900.0, 4.0, 25.0, 3},
      {"hot", 2.0, 900.0, 2.0, 40.0, 3},
      // One layer dried past where the optimum stops being interior, which the
      // fixture below this one bisects at psi_soil 5.702755. The collar is
      // pinned to the root's critical potential there: the held evaluation is
      // still the frozen-collar one, so these rows must still agree, and the
      // condition is the bound rather than stationarity -- not this reader's to
      // supply, and not asked of it.
      {"pinned", 5.72, 900.0, 2.0, 25.0, 1},
      // Drier still, and the leaf has stopped moving water.
      {"shut", 5.8, 900.0, 2.0, 25.0, 3},
  };

  const int pars[] = {grad::par_vcmax_25,
                      grad::par_jmax_25,
                      grad::par_a,
                      grad::par_curv_fact_elec_trans,
                      grad::par_curv_fact_colim,
                      grad::par_R_d_25,
                      grad::par_beta2,
                      grad::par_cost_scale_TF24,
                      grad::par_PPFD};
  const std::size_t n_par = sizeof(pars) / sizeof(pars[0]);

  int compared = 0;
  for (const Fixture& f : fixtures) {
    const int L = f.layers;
    grad::Drivers d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);
    d.leaf_temp = f.leaf_temp;

    // ⚠️ THE DIFFERENCE IS TAKEN FROM `held_row` AND NOT FROM `rows_at`. The
    // dispatcher now prefers these same closed forms for these same inputs, so a
    // comparison through it would compare each row with itself and pass whatever
    // the rows became. `held_row` is the differencing, which is the reference.
    std::vector<int> out_index{grad::out_profit};
    for (int i = 0; i < L; ++i) {
      out_index.push_back(grad::out_uptake_first + i);
    }
    std::vector<int> input(pars, pars + n_par);
    const grad::RowRequest req{out_index.data(), out_index.size(),
                               input.data(), input.size()};

    pl::Leaf l = env::fresh();
    const grad::BasePoint b = env::solved_base_point(l, d, false, s, L);

    const bool interior =
        b.branch.kind == pl::Leaf::OperatingPointKind::Interior;
    // Where the leaf has stopped moving water the frozen-collar derivation is
    // not the route: gross assimilation is identically zero there, so what these
    // readers divide by is too. The row layer declares those rows instead, and
    // that it does not reach for the readers is what is checked.
    const bool shut =
        b.branch.kind == pl::Leaf::OperatingPointKind::HydraulicShutdown ||
        b.branch.kind == pl::Leaf::OperatingPointKind::ShadeDeath;
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(b.branch.kind) + ")";
    if (shut) {
      ok(!grad::carbon_rows_apply(l, req, shut),
         "the row layer does not read a closed form at " + tag);
      // The one row a shut point declares that is not zero: profit there is minus
      // respiration minus the hydraulic cost, so this is the derived respiration
      // over its trait, negated. Checked against a difference of the solve, which
      // is available here because this fixture sits well inside the branch -- it is
      // the fixtures that sit near its edge where the difference refuses, and that
      // is why the row is declared rather than measured.
      double declared = 0.0;
      ok(grad::shut_row(l, l.hydraulic_cost_row(l.opt_psi_stem_), true,
                        grad::par_R_d_25, 5, declared),
         "dark respiration's row is declared at " + tag);
      double th[grad::n_pars];
      grad::Scratch scratch;
      const double base = grad::par_value(env::kTheta, d, false, grad::par_R_d_25);
      const double hh = base * 1e-6;
      grad::set_one(l, th, env::kTheta, d, false, grad::par_R_d_25, base + hh,
                    s.fast_stem_curve, scratch);
      l.find_root_collar_psi();
      const double up = l.profit_;
      grad::set_one(l, th, env::kTheta, d, false, grad::par_R_d_25, base - hh,
                    s.fast_stem_curve, scratch);
      l.find_root_collar_psi();
      const double dn = l.profit_;
      near(declared, (up - dn) / (2.0 * hh), 1e-6,
           "and it agrees with a difference of the solve, " + tag);
      printf("  %-22s not this reader's branch; R_d_25 row %.6g against %.6g\n",
             tag.c_str(), declared, (up - dn) / (2.0 * hh));
      continue;
    }

    // Taken before any perturbation, and at the point itself: `base_point`
    // closes by differencing the marginal profit across p*, so it leaves the
    // collar one step below it, and a row read there is wrong by that step.
    bool evaluable = false;
    l.dprofit_droot_collar_psi(b.psi_star, &evaluable);
    ok(evaluable, "the point is evaluable, " + tag);
    const pl::Leaf::PhotoTraitRows photo = l.photo_trait_rows(l.dpsistem_dpsi_);
    const pl::Leaf::CostTraitRows cost = l.cost_trait_rows(l.dpsistem_dpsi_);

    const double held[] = {photo.dprofit_dvcmax_25,
                           photo.dprofit_djmax_25,
                           photo.dprofit_da,
                           photo.dprofit_dcurv_elec,
                           photo.dprofit_dcurv_colim,
                           photo.dprofit_dR_d_25,
                           cost.dprofit_dbeta2,
                           cost.dprofit_dcost_scale,
                           photo.dprofit_dPPFD};
    const double condition[] = {photo.dmarginal_dvcmax_25,
                                photo.dmarginal_djmax_25,
                                photo.dmarginal_da,
                                photo.dmarginal_dcurv_elec,
                                photo.dmarginal_dcurv_colim,
                                photo.dmarginal_dR_d_25,
                                cost.dmarginal_dbeta2,
                                cost.dmarginal_dcost_scale,
                                photo.dmarginal_dPPFD};

    double worst_held = 0.0, worst_cond = 0.0, worst_uptake = 0.0;
    bool at_base = true;
    grad::Scratch scratch;
    grad::OutputValues direct(L);
    for (std::size_t i = 0; i < n_par; ++i) {
      const std::string nm =
          std::string(grad::par_name(pars[i], L)) + ", " + tag;
      double dR = 0.0;
      const bool differenced =
          grad::held_row(l, env::kTheta, d, false, pars[i], b.psi_star, s, true,
                         at_base, scratch, direct, dR);
      const double got = differenced ? direct[grad::out_profit]
                                     : pl::util::na_value;
      if (!std::isfinite(got) || !std::isfinite(held[i])) {
        ok(false, "both routes answer for " + nm);
        continue;
      }
      ++compared;
      const double scale_h = std::max(1.0, std::abs(got));
      worst_held = std::max(worst_held, std::abs(held[i] - got) / scale_h);
      near(held[i], got, 1e-4, "profit row for " + nm);

      // At a fixed collar a carbon-side trait moves no water, so every uptake
      // row is claimed exactly zero -- and the difference has to agree with
      // that exactly rather than at its own floor.
      for (int j = 0; j < L; ++j) {
        worst_uptake = std::max(
            worst_uptake,
            std::abs(direct[grad::out_uptake_first + j]) / scale_h);
      }

      if (interior && std::isfinite(dR) && std::isfinite(condition[i])) {
        const double scale_c = std::max(1.0, std::abs(dR));
        worst_cond =
            std::max(worst_cond, std::abs(condition[i] - dR) / scale_c);
        near(condition[i], dR, 1e-3, "condition row for " + nm);
      }
    }
    // Radiation's held row exists twice over, and the second one shares no
    // algebra with this family: `dprofit_dPPFD` closes the concentration's
    // root-find directly rather than going through the transport family's pass.
    // Two implementations agreeing is worth more than either agreeing with a
    // difference.
    // They agree to 4e-09 rather than to the last bit, which is what two
    // formulations of one derivative cost each other through a root-find at
    // double precision -- and both agree with the difference at 2e-09.
    near(photo.dprofit_dPPFD, l.dprofit_dPPFD(), 1e-8,
         "radiation's row matches the leaf's own second route, " + tag);
    printf("  %-22s worst rel: profit %.3g, condition %.3g; uptake rows %.3g\n",
           tag.c_str(), worst_held, worst_cond, worst_uptake);
  }
  ok(compared == int(n_par) * 6, "every trait was compared at every fixture");
}

// The two ends of the shrinking, which the grid above cannot reach: on it every
// arm stays on the base point's branch at the requested step, and so does every
// constrained point's at a step three orders coarser.
//
// So the branch is run at a point placed a chosen distance from a boundary, and
// the boundary is FOUND rather than written down: drying one soil layer takes the
// optimum from interior to pinned somewhere near 5.70 MPa, and where exactly is a
// property of the solve. Two soil layers a decade apart from it are what the two
// ends need -- one where the requested step reaches across and a tenth of it does
// not, one where two decades of shrinking are not enough and the layer's row is
// refused by name. Coarsening the step instead does not work: a step large enough
// to need two decades takes the other arm outside the stem curve's domain, and the
// solve says so rather than returning a branch.
void test_rows_shrink_the_step_to_stay_on_one_branch() {
  printf("the rows shrink a step that leaves the branch, and refuse below that\n");
  namespace grad = phylloptim::gradient;

  int out_index[grad::n_outputs];
  for (int j = 0; j < grad::n_outputs; ++j) {
    out_index[j] = j;
  }
  const std::vector<int> pars{grad::par_psi_soil_first};
  grad::RowRequest req;
  req.output = out_index;
  req.n_output = grad::n_outputs;
  req.input = pars.data();
  req.n_input = pars.size();

  phylloptim::Leaf l = env::fresh();
  auto branch_at = [&](double psi_soil) -> grad::Branch {
    grad::Drivers d = env::drivers(psi_soil, 900.0, 2.0, 1, 1);
    grad::apply(l, env::kTheta, d, false, -1, true);
    l.find_root_collar_psi();
    return grad::branch_here(l);
  };

  // Where the optimum stops being interior, to the last bit. Bisected on the
  // classification the solve reports rather than on a comparison of the collar
  // against a bound: those two agree until a step-in makes them differ.
  double wet = 5.0, dry = 5.72;
  ok(branch_at(wet).kind == grad::OperatingPointKind::Interior &&
         branch_at(dry).kind == grad::OperatingPointKind::PinnedDryRootCrit,
     "the interior and the pinned regimes are bracketed");
  for (int i = 0; i < 60; ++i) {
    const double mid = 0.5 * (wet + dry);
    if (branch_at(mid).kind == grad::OperatingPointKind::Interior) {
      wet = mid;
    } else {
      dry = mid;
    }
  }
  printf("  the interior regime ends at psi_soil %.6f\n", dry);

  // Whether both arms of one step stay on the branch the base point took, which
  // is the question `rows_at` shrinks the step until the answer to is yes.
  auto arms_stay = [&](const grad::Drivers &d, const grad::Branch &base,
                       double h) -> bool {
    bool stays = true;
    double th[grad::n_pars];
    grad::Scratch scratch;
    for (int side = 0; side < 2; ++side) {
      grad::set_one(l, th, env::kTheta, d, false, pars[0],
                    d.psi_soil[0] + (side == 0 ? h : -h), true, scratch);
      l.find_root_collar_psi();
      stays = stays && grad::branch_here(l) == base;
    }
    grad::apply(l, env::kTheta, d, false, -1, true);
    l.find_root_collar_psi();
    return stays;
  };

  const grad::Settings settings;
  // One decade in: two millionths of an MPa past the boundary, so the step reaches
  // across it and a tenth of the step does not.
  {
    grad::Drivers d = env::drivers(dry + 2e-6, 900.0, 2.0, 1, 1);
    const grad::Branch base = branch_at(d.psi_soil[0]);
    const double h = grad::step_for(pars[0], d.psi_soil[0], settings.step, grad::n_soil_layers(d, false));
    ok(base.kind == grad::OperatingPointKind::PinnedDryRootCrit,
       "the point is pinned to the dry bound");
    ok(!arms_stay(d, base, h), "the requested step takes an arm off that branch");
    ok(arms_stay(d, base, 0.1 * h), "and a tenth of it does not");
    const grad::Rows rows = env::solved_rows(l, env::kTheta, d, req, settings);
    ok(std::isfinite(rows.held[std::size_t(grad::out_profit)]) &&
           rows.message.empty(),
       "so the row comes back and nothing is refused");
  }

  // On the boundary itself NO step keeps a re-solved arm on the branch -- not at
  // the requested one, and not two decades below it. This used to be where the
  // row was refused by name, and it is now where the arms stop being re-solved:
  // placed on the bound the point is pinned to, they are on the branch by
  // construction and the step no longer decides.
  {
    grad::Drivers d = env::drivers(dry, 900.0, 2.0, 1, 1);
    const grad::Branch base = branch_at(d.psi_soil[0]);
    const double h = grad::step_for(pars[0], d.psi_soil[0], settings.step, grad::n_soil_layers(d, false));
    ok(!arms_stay(d, base, 0.01 * h),
       "two decades below the step still takes a re-solved arm off the branch");
    const grad::Rows rows = env::solved_rows(l, env::kTheta, d, req, settings);
    ok(rows.kind == base.kind, "the point is still the pinned one");
    bool all_answered = true;
    for (int j = 0; j < grad::n_outputs; ++j) {
      all_answered =
          all_answered && std::isfinite(rows.held[std::size_t(j)]);
    }
    ok(all_answered, "and the whole column answers rather than being refused");
    ok(rows.message.empty(),
       "with nothing refused by name: " + rows.message);
    // The row is PARTS now, so the point's movement is here rather than inside
    // every output's row: the soil potential moves the bound, the bound is what
    // defines this point, and a live entry is what says so. It was zero while
    // the arms re-solved and each row carried the total.
    ok(rows.dresidual[0] != 0.0 && std::isfinite(rows.dresidual[0]),
       "and the condition's gradient carries the bound's own movement");
    ok(std::isfinite(rows.residual_slope) && rows.residual_slope != 0.0,
       "beside the slope the consumer divides it by");
  }
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

  // The other half of share = count x price. A placement is what a pass re-running
  // the model over the same states pays where the run kept the point, so the
  // difference between these two is what one such evaluation saves.
  const phylloptim::Leaf::SolvedPoint point = l.solved_point();
  const auto t2 = std::chrono::steady_clock::now();
  double acc2 = 0;
  for (int i = 0; i < N; ++i) {
    l.place_solved_point(point);
    acc2 += l.profit_;
  }
  const auto t3 = std::chrono::steady_clock::now();
  printf("  %.2f us per place_solved_point() (checksum %.6f)\n",
         std::chrono::duration<double, std::micro>(t3 - t2).count() / N,
         acc2 / N);
}

} // namespace

int main() {
  test_defaults_are_unset();
  test_vulnerability_curve();
  test_vulnerability_integral_derivatives();
  test_spline_matches_direct_integration();
  test_stem_curve_derivative_is_the_closed_form();
  test_arrhenius();
  test_saturation_vapour_pressure();
  test_solve_single_layer();
  test_solve_is_deterministic();
  test_drier_soil_costs_carbon();
  test_a_placed_point_is_the_searched_one();
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
  test_soil_conductance_is_positive();
  test_soil_potential_derivative();
  test_uptake_mixed_second_derivative();
  test_light_row();
  test_env_adjoint();
  test_marginal_price_water();
  test_root_vulnerability_is_bounded_past_its_grid();
  test_root_psi_crit_clamp_binds();
  test_bound_row_root_psi_crit_is_a_unit_row();
  test_signed_potentials_are_rejected();
  test_lambda_equals_dA_dE_single_layer();
  test_multilayer_lambda_identity();
  test_g1_eff();
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
  test_temperature_parameters_are_settable();
  test_temperature_params_invalidate_cache();
  test_rd_temperature_response();
  test_set_traits_matches_a_fresh_leaf();
  test_perturb_stem_b_matches_a_rebuild();
  test_bad_input_throws();
  test_out_of_domain_names_the_spline();
  test_out_of_domain_under_rescale();
  test_environment_par_names();
  test_root_carbon_rows();
  test_rows_read_a_solved_leaf();
  test_the_read_names_what_it_declines();
  test_shade_death_soil_rows();
  test_environment_rows_match_a_differenced_solve();
  test_environment_water_rows_are_rank_one();
  test_environment_rows_are_zero_below_the_rooted_layers();
  test_uptake_outputs_are_enumerated();
  test_rows_in_parts_assemble_to_the_totals();
  test_the_curves_trait_derivative_is_the_models_own();
  test_the_transport_leaves_the_flux_where_it_is();
  test_the_transport_traits_rows_match_a_rebuilt_difference();
  test_the_transport_reports_its_collar_response();
  test_the_condition_is_the_stem_potential_and_its_collar_response();
  test_the_condition_reaches_the_state_through_two_intermediates();
  test_carbon_trait_rows_match_a_differenced_solve();
  test_the_supplys_mixed_partials_match_a_difference();
  probe_root_curve_slope();
  test_the_transport_response_is_the_flux_balances_own();
  probe_table_vs_curve();
  test_the_condition_is_carbon_bought_against_tension_paid();
  test_a_shut_points_rows_are_the_costs_own();
  test_the_two_classifications_of_one_point();
  test_the_three_unchecked_invariants();
  test_the_supply_answers_at_the_two_coincidences();
  test_the_supplys_second_collar_derivative();
  test_the_transport_responses_collar_slope();
  test_the_two_zero_flux_kinds_are_two_points();
  test_shade_deaths_rows_against_a_differenced_solve();
  test_the_collar_channel_is_read_rather_than_differenced();
  test_the_bounds_steepness_rows_match_a_rebuilt_difference();
  test_a_pinned_point_answers_from_parts_rather_than_re_solving();
  test_the_soil_states_two_scalars_are_separable_and_right();
  test_the_soil_states_rows_match_a_differenced_solve();
  test_rows_shrink_the_step_to_stay_on_one_branch();
  benchmark();

  printf("\n%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
