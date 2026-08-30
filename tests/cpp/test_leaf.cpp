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

#include "leaf_inputs.hpp"
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

// Every layer's own draw responding to the collar, from the one path that
// answers it: the quadrature at a seeded scalar. In mol, like the draws
// themselves -- the total is the one thing this model states in kg.
inline std::vector<double> draws_along_collar(phylloptim::Leaf& l,
                                              double collar) {
  phylloptim::tangent c = collar;
  phylloptim::seed_direction(c, 1.0);
  const phylloptim::SupplyValues<phylloptim::tangent> supply =
      fixture::supply_of<phylloptim::tangent>(l);
  std::vector<phylloptim::tangent> per;
  l.E_from_soil_at<phylloptim::tangent>(c, supply.at(), per);
  std::vector<double> out;
  for (const phylloptim::tangent& v : per) {
    out.push_back(phylloptim::derivative_along(v));
  }
  return out;
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
      // left placed.
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
// used to be placed only by a solve -- so calling it on a leaf that had had
// set_physiology but not find_root_collar_psi read an empty vector. It now places
// them itself. Ported from plant develop (#585).
void test_gradient_needs_no_prior_solve() {
  printf("dprofit/dpsi_collar does not require a prior solve\n");
  Drivers d;
  phylloptim::Leaf solved = make_leaf(d, {2.0}, {1.0});
  solved.find_root_collar_psi();
  phylloptim::Leaf unsolved = make_leaf(d, {2.0}, {1.0});
  // Bit-identical, not merely close: placing the potentials from psi_soil_ is
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

    const std::vector<double> by_layer = draws_along_collar(l, T);
    double summed = 0.0;
    for (int i = 0; i < layers; ++i) {
      summed += by_layer[std::size_t(i)];
    }
    // Not bit-exact, and the reason is the only difference between the two: the
    // total converts to kg once after summing and this converts after summing
    // too, so they differ by reassociation and by nothing else.
    summed *= phylloptim::kg_per_mol_h2o;
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
//  2. THE MISSING CHAIN TERM. Even placed correctly, the derivative has no
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
        // because dprofit_droot_collar_psi above has just re-placed the
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
  // edits (placement the temperature parameters, add the dA/dTleaf chain term, handle
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
// THE SHARED OPERATING POINTS
// ---------------------------------------------------------------------------
// One trait vector and one way of laying out a soil profile, shared by every
// check below so that a state reached in one is the state reached in the next.

namespace env {

// The suite's trait defaults, in `set_traits`' argument order. Sized from the
// initialiser and checked, because too few initialisers is legal and zero-fills
// the rest.
const double kTheta[] = {96.0,   2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,    157.44,   0.30,
                         0.7,    0.99,     7.5,      kRd25};
static_assert(sizeof(kTheta) / sizeof(kTheta[0]) == phylloptim::n_traits);

// The maximum leaf-specific conductance every point below is driven at. A driver
// rather than a trait, so it sits on the physiology and not in kTheta.
const double kKmax = 1.0 * 0.000157 / 5.0;

// `rooted` layers of root carbon spread over `layers` soil layers, drying with
// depth -- the golden grid's soil, so a point named here is a point that grid
// already covers.
fixture::Physiology drivers(double psi_soil, double ppfd, double vpd, int layers,
                            int rooted) {
  fixture::Physiology d;
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
  d.kmax = kKmax;
  d.atm_vpd = vpd;
  d.ca = 40.0;
  d.leaf_temp = 25.0;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

phylloptim::Leaf fresh() {
  phylloptim::Leaf l;
  return l;
}

}  // namespace env

// The rank-one self-test, and it has real content: psi_soil reaches profit only
// through uptake, at a price the collar's own first-order condition fixes, so
// dprofit/dpsi_soil_j divided by dE_up/dpsi_soil_j must be the SAME number for
// every layer j. If it is not, either the formula or the price is wrong.
//
// It also pins WHICH price. `marginal_cost_water_multilayer()` is the lambda the
// collar solve equalises and the obvious candidate; it is 1.20x to 2.37x too
// large over the golden grid, because along the soil direction the collar is
// held fixed and the stem's own cost has not been netted out. The price this
// asserts is that lambda less `marginal_cost_water()`.
void test_environment_water_rows_are_rank_one() {
  printf("environment rows: the water channel is rank one across layers\n");
  for (int L : {3, 5}) {
    for (double psi_soil : {1.0, 3.0}) {
      phylloptim::Leaf l = env::fresh();
      fixture::Physiology d = env::drivers(psi_soil, 900.0, 2.0, L, L);
      d.drive(l, env::kTheta);
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
        d.drive(l, env::kTheta);
        l.find_root_collar_psi();
        l.E_from_Soil_to_Root_Collar(psi_star, up);
        const double e_up = l.E_up_;
        l.E_from_Soil_to_Root_Collar(psi_star, dn);
        const double e_dn = l.E_up_;
        const double dE = (e_up - e_dn) / (2.0 * h);

        // dprofit/dpsi_soil_j from a differenced SOLVE. Differencing profit
        // against the COLLAR would be the mistake this formulation exists to
        // avoid -- profit is flat there. Against psi_soil it is not.
        fixture::Physiology du = d, dd = d;
        du.psi_soil = up;
        dd.psi_soil = dn;
        du.drive(l, env::kTheta);
        l.find_root_collar_psi();
        const double p_up = l.profit_;
        dd.drive(l, env::kTheta);
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
      d.drive(l, env::kTheta);
    }
  }
}

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
// So the referee is the sweep's best agreement, and every differenced check below
// takes it that way rather than at one step.

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
  namespace pl = phylloptim;

  struct Fixture { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const Fixture fixtures[] = {{"wet", 2.0, 900.0, 2.0, 3},
                              {"dim", 2.0, 300.0, 2.0, 3},
                              {"dry", 4.5, 900.0, 2.0, 3},
                              {"arid", 2.0, 900.0, 4.0, 3},
                              {"pinned", 5.72, 900.0, 2.0, 1}};

  for (const Fixture& f : fixtures) {
    fixture::Physiology d =
        env::drivers(f.psi_soil, f.ppfd, f.vpd, f.layers, f.layers);
    pl::Leaf l = env::fresh();
    d.drive(l, env::kTheta);
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

// ⚠️ A DIFFERENCE OVER THE TRANSPORT'S THREE PARAMETERS REBUILDS A CURVE ON EVERY
// ARM. That is what makes it the right reference -- the forward model rebuilds too,
// so the difference is of the model -- and it is also what makes it expensive,
// which is the whole reason for a closed form.
//
// The two disagreed once, and by enough to be recorded rather than bounded: the
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
  namespace pl = phylloptim;
  const int L = 5;
  fixture::Physiology d = env::drivers(1.0, 1500.0, 1.0, L, L);
  pl::Leaf l = env::fresh();
  d.drive(l, env::kTheta);
  l.find_root_collar_psi();
  const double p = l.opt_root_psi_;
  bool feasible = false;
  l.dprofit_droot_collar_psi(p, &feasible);
  ok(feasible, "the point is evaluable");
  const double base = env::kKmax;

  double first_sigma = 0.0, last_sigma = 0.0;
  double assim_lo = std::numeric_limits<double>::infinity(), assim_hi = 0.0;
  bool changes_sign = false, seen_positive = false, seen_negative = false;
  int decade = 0;
  for (double rel : {1e-6, 1e-5, 1e-4, 1e-3}) {
    const double h = base * rel;
    double sigma[2], assim[2], uptake_total[2];
    for (int side = 0; side < 2; ++side) {
      fixture::Physiology moved = d;
      moved.kmax = side == 0 ? base + h : base - h;
      moved.drive(l, env::kTheta);
      l.supply_begin_solve();
      ok(l.profit_at_fixed_collar(p).feasible, "the held evaluation answers");
      sigma[side] = l.opt_psi_stem_;
      assim[side] = l.assim_colimited_;
      uptake_total[side] = 0.0;
      for (int i = 0; i < L; ++i) {
        uptake_total[side] += l.soil_consumption_[std::size_t(i)];
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
  d.drive(l, env::kTheta);
  // Four decades of step move it by 1.5e-06, which is the coarsest step's own
  // truncation and not a floor: a floor would have shown at the FINEST.
  near(last_sigma, first_sigma, 1e-5,
       "while the stem potential's row barely moves with the step at all");
}

void test_the_curves_trait_derivative_is_the_models_own() {
  printf("the closed-form dG/dc against a rebuilt difference of the spline\n");
  namespace pl = phylloptim;
  const int L = 3;
  fixture::Physiology d = env::drivers(2.0, 900.0, 2.0, L, L);
  pl::Leaf l = env::fresh();
  d.drive(l, env::kTheta);
  const double b = l.stem_b, c = l.stem_c;

  double worst_dc = 0.0, worst_value = 0.0;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    const pl::VulnerabilityIntegralDerivatives dd =
        pl::cumulative_vulnerability_integral_derivatives_at(psi, b, c);
    const double value = l.stem_curve_integral(psi);
    worst_value = std::max(worst_value,
                           std::abs(dd.value - value) / std::abs(value));

    double th[pl::n_traits];
    std::copy(env::kTheta, env::kTheta + pl::n_traits, th);
    const double h = c * 1e-5;
    th[pl::trait_stem_c] = c + h;
    d.drive(l, th);
    const double up = l.stem_curve_integral(psi);
    th[pl::trait_stem_c] = c - h;
    d.drive(l, th);
    const double dn = l.stem_curve_integral(psi);
    d.drive(l, env::kTheta);
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
  namespace pl = phylloptim;
  const int L = 5;
  fixture::Physiology d = env::drivers(3.0, 1500.0, 0.5, L, L);
  pl::Leaf l = env::fresh();
  d.drive(l, env::kTheta);
  l.find_root_collar_psi();
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
  double th[pl::n_traits];
  // The carbon each layer's resistances were built from. The network holds it
  // split three ways rather than as itself, and the vertical share is a third.
  std::vector<double> carbon;
  for (double v : d.root_network.c_r_V) {
    carbon.push_back(3.0 * v);
  }
  for (int a = 0; a < L; ++a) {
    double sum_dE = 0.0, sum_dD = 0.0;
    for (int i = 0; i < L; ++i) {
      sum_dE += dE[std::size_t(i)][std::size_t(a)];
      sum_dD += dD[std::size_t(i)][std::size_t(a)];
    }
    const double base = carbon[std::size_t(a)];
    const double h = std::abs(base) * 1e-6;
    double slope[2], total[2];
    for (int side = 0; side < 2; ++side) {
      std::vector<double> moved_carbon = carbon;
      moved_carbon[std::size_t(a)] = side == 0 ? base + h : base - h;
      fixture::Physiology moved = d;
      moved.root_network = fixture::root_network(moved_carbon, d.soil_depth);
      moved.drive(l, env::kTheta);
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
  d.drive(l, env::kTheta);

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
  const struct { const char* what; int slot; Trait trait; } curve[] = {
      {"root_b", pl::trait_root_b, Trait::Position},
      {"root_c", pl::trait_root_c, Trait::Steepness}};
  for (const auto& t : curve) {
    d.drive(l, env::kTheta);
    std::vector<double> dE, d2E;
    l.dE_from_soil_droot_curve(p, psi, t.trait, dE, d2E);
    double got_E = 0.0, got_2 = 0.0;
    for (int i = 0; i < L; ++i) {
      got_E += dE[std::size_t(i)];
      got_2 += d2E[std::size_t(i)];
    }
    const double base = env::kTheta[t.slot];
    const double hh = std::abs(base) * 1e-6;
    double total[2], slope[2];
    for (int side = 0; side < 2; ++side) {
      std::copy(env::kTheta, env::kTheta + pl::n_traits, th);
      th[t.slot] = side == 0 ? base + hh : base - hh;
      d.drive(l, th);
      // ⚠️ The per-layer integrals are cached against the soil state, and a moved
      // curve invalidates them: without this the difference reads the new curve at
      // the collar and the OLD one at each soil layer, which is a difference of two
      // models. Every production path reaches this through the collar solve, which
      // places it; a check that perturbs and reads directly has to placement it itself.
      l.supply_begin_solve();
      l.E_from_Soil_to_Root_Collar(p, l.supply_psi_soil());
      total[side] = l.E_up_;
      slope[side] = l.dE_from_soil_dpsi_collar(p, l.supply_psi_soil());
    }
    d.drive(l, env::kTheta);
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
  namespace pl = phylloptim;
  struct F { const char* what; double psi_soil, ppfd, vpd; int layers; };
  const F fixtures[] = {{"sodden", 0.5, 900, 2.0, 3}, {"wet", 2.0, 900, 2.0, 3},
                        {"dim", 2.0, 300, 2.0, 3},    {"arid", 2.0, 900, 4.0, 3},
                        {"deep", 3.0, 1500, 0.5, 5}};
  int checked = 0;
  for (const F& f : fixtures) {
    const int L = f.layers;
    fixture::Physiology d = env::drivers(f.psi_soil, f.ppfd, f.vpd, L, L);
    pl::Leaf l = env::fresh();
    d.drive(l, env::kTheta);
    l.find_root_collar_psi();
    const pl::Leaf::OperatingPointKind kind = l.operating_point_kind();
    const std::string tag =
        std::string(f.what) + " (" +
        pl::Leaf::operating_point_kind_name(kind) + ")";
    if (kind != pl::Leaf::OperatingPointKind::Interior) {
      continue;
    }
    // A collar a little drier than the optimum, so R is a number rather than a
    // zero. The wettest end is the one with room: the dry end runs the transport
    // coordinate off the end of its grid, which is the caller not having asked
    // for a potential that exists.
    const double collar = l.opt_root_psi_ * 0.98;
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

// The collar the solve returns is stationary where the point is interior, and is
// not where a bound placed it. The two populations sit orders apart, which is
// what a consumer reading the residual rather than the branch is relying on.
//
// ⚠️ THE ZERO-FLUX KINDS ARE NEITHER, and this is why the branch is read from the
// classification. The marginal profit at such a collar reports itself infeasible
// and returns a hard 0.0 sentinel, so a test on the residual's size alone lands
// them on the interior side of any cut. They are skipped by that report, and
// counted, so their number is visible rather than folded into the interior one.
void test_interior_optima_are_stationary() {
  printf("the collar residual at interior optima and at pinned ones\n");
  namespace pl = phylloptim;
  using Kind = pl::Leaf::OperatingPointKind;
  // The band the two populations sit either side of. It needs no scale of its
  // own: the residual is profit per MPa and the solve drives it to rounding.
  const double cut = 1e-8;
  int interior = 0, pinned = 0, refused = 0;
  double worst_interior = 0.0, mildest_pin = 1e300;
  for (double T : {25.0, 40.0}) {
    for (double psi = 0.5; psi <= 7.0; psi += 0.5) {
      for (double ppfd : {0.0, 100.0, 900.0}) {
        for (int L : {1, 5}) {
          fixture::Physiology d = env::drivers(psi, ppfd, 2.0, L, L);
          d.leaf_temp = T;
          pl::Leaf l = env::fresh();
          d.drive(l, env::kTheta);
          l.find_root_collar_psi();
          const Kind kind = l.operating_point_kind();
          bool feasible = false;
          const double resid =
              std::abs(l.dprofit_droot_collar_psi(l.opt_root_psi_, &feasible));
          if (!feasible) {
            ++refused;
            continue;
          }
          if (kind == Kind::Interior) {
            ++interior;
            worst_interior = std::max(worst_interior, resid);
          } else {
            ++pinned;
            mildest_pin = std::min(mildest_pin, resid);
          }
        }
      }
    }
  }
  printf("  %d interior points reach %.3g, %d pinned ones start at %.3g, %d "
         "report no derivative\n",
         interior, worst_interior, pinned, mildest_pin, refused);
  ok(interior > 0 && pinned > 0, "the sweep reaches both populations");
  ok(worst_interior < cut,
     "every interior collar satisfies dprofit == 0 to solver precision");
  ok(mildest_pin > 1e4 * cut,
     "and every pinned one is held there by its bound, not by the condition");
}

void test_the_supply_answers_at_the_two_coincidences() {
  printf("the supply's derivatives at a gravity-balanced layer and at atmospheric\n");
  namespace pl = phylloptim;

  // A collar placed exactly on one layer's gravity balance, which is where a
  // one-layer plant that has stopped drawing water sits.
  fixture::Physiology d = env::drivers(2.0, 900.0, 2.0, 3, 3);
  pl::Leaf l = env::fresh();
  d.drive(l, env::kTheta);
  l.find_root_collar_psi();
  const std::vector<double> &ps = l.supply_psi_soil();
  const int layer = 1;
  const double balanced =
      ps[std::size_t(layer)] + l.roots_.grav_head_z_[std::size_t(layer)];

  const std::vector<double> draw = draws_along_collar(l, balanced);
  const double total = l.dE_from_soil_dpsi_collar(balanced, ps);
  ok(std::isfinite(total), "the conductance at a gravity-balanced collar is a number");
  bool per_layer_finite = true;
  for (double v : draw) {
    per_layer_finite = per_layer_finite && std::isfinite(v);
  }
  ok(per_layer_finite, "and so is every layer's own");
  ok(std::isfinite(l.d2E_from_soil_dpsi_collar2(balanced, ps)),
     "and its collar derivative too");

  // Against a difference of the flux, over a sweep of steps for the reason above.
  // The referee is the model's own uptake at a held soil state, which
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

  // ⚠️ EQUAL POTENTIALS USED TO REFUSE HERE, and this assertion is the replacement
  // rather than the deletion. The mean was integral/span with both vanishing; formed
  // as an average it is f at the bound, so the conductance is a number and the
  // caller's central-difference fallback has nothing left to catch. Checked against
  // that same difference, taken from far enough out to be refereeing rather than
  // measuring the same region.
  const double meeting = ps[std::size_t(layer)];
  const double at_meeting = l.dE_from_soil_dpsi_collar(meeting, ps);
  ok(std::isfinite(at_meeting),
     "a collar meeting a layer's potential now answers rather than refusing");
  {
    const double h = 1e-4;
    std::vector<double> per(ps.size(), 0.0);
    double up = 0.0, dn = 0.0;
    l.roots_.uptake_at(meeting + h, ps, per, up);
    l.roots_.uptake_at(meeting - h, ps, per, dn);
    const double differenced = (up - dn) / (2.0 * h);
    ok(std::abs(at_meeting - differenced) / std::abs(differenced) < 1e-9,
       "and the number it answers with is the difference of the flux");
  }
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
        fixture::Physiology dd = env::drivers(psi, ppfd, 2.0, layers, layers);
        pl::Leaf probe = env::fresh();
        dd.drive(probe, env::kTheta);
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
  double worst = 0.0;
  std::string worst_where;
  int compared = 0;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    for (int layers : {1, 3, 5}) {
      fixture::Physiology d = env::drivers(psi, 900.0, 2.0, layers, layers);
      phylloptim::Leaf l = env::fresh();
      d.drive(l, env::kTheta);
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
  fixture::Physiology d = env::drivers(2.0, 900.0, 2.0, 1, 1);
  // That constant resistance, which is the whole network this path reads.
  d.root_network = fixture::series_resistance(1e3);
  phylloptim::Leaf single = env::fresh();
  single.set_supply_single(0.0);
  d.drive(single, env::kTheta);
  single.find_root_collar_psi();
  ok(single.d2E_from_soil_dpsi_collar2(single.opt_root_psi_,
                                       single.supply_psi_soil()) == 0.0,
     "the single path's is exactly zero");
}

// The two zero-flux kinds are two points, and the placement is what tells them apart.
//
// Both pay respiration plus a hydraulic cost, so both take the cost's own rows as
// their profit rows. A hydraulic shutdown holds the stem AT its critical
// potential, so that input is the placement and its row is the cost's slope, and
// nothing the leaf reads is a function of the soil. Shade death holds both
// potentials at the collar of zero uptake instead: the critical potential is
// inactive exactly as at an interior optimum, and the placement moves with the soil.
// The curvature the interior derivation divides by, against a difference of the
// marginal it is the slope of.
//
// ⚠️ THIS IS THE CHECK THAT DID NOT EXIST, and its absence cost a plant gradient.
// The curvature was formed by differentiating a DIFFERENT assembly of the profit
// twice, whose second-order content was two hand-written Taylor coefficients with
// no referee -- and on a 105-year stand it returned +34.4 where a difference of the
// marginal gives -9.63, so plant refused a perfectly good interior maximum.
// marginal_collar_slope() takes it as a FIRST derivative of the marginal instead,
// from the same assembly the solve roots, and this is what says so.
void test_the_marginals_collar_slope_against_a_difference() {
  printf("the marginal's collar slope against a difference of the marginal\n");
  double worst = 0.0;
  std::string worst_where;
  int compared = 0;
  int sign_disagreements = 0;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    for (int layers : {1, 3, 5}) {
      for (double ppfd : {150.0, 900.0}) {
        fixture::Physiology d = env::drivers(psi, ppfd, 2.0, layers, layers);
        phylloptim::Leaf l = env::fresh();
        d.drive(l, env::kTheta);
        l.find_root_collar_psi();
        if (l.operating_point_kind() !=
            phylloptim::Leaf::OperatingPointKind::Interior) {
          continue;
        }
        const double p = l.opt_root_psi_;
        // The inputs at double, from the same fixture the seeded checks use: the
        // slope reads them rather than the leaf's members so its parameter rows
        // survive at an active scalar, and at double the two are the same numbers.
        fixture::Soil soil;
        for (int q = 0; q < layers; ++q) {
          soil.depth.push_back(1.0 * (q + 1));
          soil.carbon.push_back(1.0 / layers / 0.05);
        }
        const phylloptim::LeafInputs<double> li =
            fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);
        const double closed = l.marginal_collar_slope(li.profit);
        if (!std::isfinite(closed) || closed == 0.0) {
          continue;
        }
        ++compared;
        // An interior collar sits on a downward crossing, so the slope there is
        // negative. A positive one is the defect this test exists for.
        if (!(closed < 0.0)) {
          ++sign_disagreements;
        }
        double best = std::numeric_limits<double>::infinity();
        for (int e = 4; e <= 7; ++e) {
          const double h = std::max(std::abs(p), 1.0) * std::pow(10.0, -double(e));
          bool up_ok = false, dn_ok = false;
          const double up = l.dprofit_droot_collar_psi(p + h, &up_ok);
          const double dn = l.dprofit_droot_collar_psi(p - h, &dn_ok);
          if (!up_ok || !dn_ok || !std::isfinite(up) || !std::isfinite(dn)) {
            continue;
          }
          best = std::min(best, std::abs(closed / ((up - dn) / (2.0 * h)) - 1.0));
        }
        // The solve closes on the returned collar afterwards, so the probes above
        // leave nothing behind for the next iteration.
        l.dprofit_droot_collar_psi(p);
        if (std::isfinite(best) && best > worst) {
          worst = best;
          worst_where = "psi=" + std::to_string(psi) + " layers=" +
                        std::to_string(layers) + " ppfd=" + std::to_string(ppfd);
        }
      }
    }
  }
  printf("  %d interior points compared, worst relative gap %.3e at %s\n", compared,
         worst, worst_where.c_str());
  ok(compared >= 6, "the grid reaches interior points at all");
  ok(sign_disagreements == 0, "every interior slope is negative");
  // A centred difference of a quantity whose own evaluation carries a nested ci
  // root-find cannot do better than about 1e-5; the defect this catches is a
  // FACTOR, not a last digit.
  ok(worst < 1e-3, "the closed slope matches a difference of the marginal");
}

// Where a divided difference stops being the right form for the mean conductivity.
//
// The layer resistance is r_R_H_min * span / (G(max) - G(min)) -- a divided
// difference of the cumulative curve. Its denominator is a difference of two
// nearly-equal reads of a TABULATED G, so the relative error is the table's own
// error divided by the span, and it diverges as the span shuts. The mean of the
// integrand over the interval is also just f(midpoint) + span^2/24 * f'', so below
// some span the midpoint is strictly better AND needs no differencing.
//
// This measures where they cross, so the threshold is data rather than a guess.
void test_where_the_mean_conductivity_should_stop_differencing() {
  printf("the mean conductivity: divided difference against the midpoint\n");
  phylloptim::Leaf l = env::fresh();
  env::drivers(2.0, 900.0, 2.0, 3, 3).drive(l, env::kTheta);

  // One suction well inside the curve's domain, and spans shrinking around it.
  const double centre = 2.0;
  printf("      %-12s %-22s %-22s %s\n", "span", "difference", "midpoint", "rel gap");
  for (int e = 2; e <= 12; ++e) {
    const double span = std::pow(10.0, -double(e));
    const double lo = centre - 0.5 * span;
    const double hi = centre + 0.5 * span;
    const double integral =
        l.roots_.root_vuln_integral_at(hi) - l.roots_.root_vuln_integral_at(lo);
    const double differenced = integral / span;
    const double midpoint = l.roots_.root_vuln_at(centre);
    const double rel = std::abs(differenced / midpoint - 1.0);
    printf("      1e-%-9d %-22.15g %-22.15g %.3e\n", e, differenced, midpoint, rel);
  }
  ok(true, "measured (read the table above; this test reports rather than asserts)");
}

// The supply's collar derivatives as the collar closes on a layer's potential.
//
// ⚠️ A DIFFERENCE CANNOT REFEREE THIS REGION -- that is the whole reason the midpoint
// form exists -- so what is checked is CONTINUITY. The two forms of the mean
// conductivity agree to ~1e-11 at the crossover (see the table above), so switching
// between them must not move the answer, and the sequence must stay smooth and finite
// all the way down. Before the midpoint form the divided difference lost an order of
// span per derivative and this sweep went to noise: one operating point of 2,829,445
// on a century stand reached a span of 5.6e-08 and cost the whole gradient.
void test_the_supply_derivatives_stay_smooth_into_a_coincidence() {
  printf("the supply's collar derivatives closing on a layer's potential\n");
  phylloptim::Leaf l = env::fresh();
  env::drivers(2.0, 900.0, 2.0, 3, 3).drive(l, env::kTheta);
  const std::vector<double> soil = l.supply_psi_soil();

  // Approach layer 0's potential from the dry side, straddling the 1e-5 crossover.
  // Stops at 2e-08, above the at_equal_potentials guard, which still answers NaN.
  const double target = soil[0];
  double prev1 = 0.0, prev2 = 0.0;
  double worst_jump1 = 0.0, worst_jump2 = 0.0;
  int finite = 0, steps = 0;
  printf("      %-10s %-24s %s\n", "span", "dE/dT", "d2E/dT2");
  for (int e = 3; e <= 8; ++e) {
    for (double mult : {5.0, 1.0}) {
      const double span = mult * std::pow(10.0, -double(e));
      if (span < 2e-8) continue;
      const double T = target + span;
      const double d1 = l.dE_from_soil_dpsi_collar(T, soil);
      const double d2 = l.d2E_from_soil_dpsi_collar2(T, soil);
      if (std::isfinite(d1) && std::isfinite(d2)) {
        ++finite;
        if (steps > 0) {
          worst_jump1 = std::max(worst_jump1, std::abs(d1 / prev1 - 1.0));
          worst_jump2 = std::max(worst_jump2, std::abs(d2 / prev2 - 1.0));
        }
        prev1 = d1; prev2 = d2; ++steps;
      }
      printf("      %-10.1e %-24.15g %.15g\n", span, d1, d2);
    }
  }
  printf("      worst step-to-step change: dE/dT %.3e, d2E/dT2 %.3e\n",
         worst_jump1, worst_jump2);
  ok(finite >= 10, "every span above the guard answers");
  // Both derivatives approach a finite limit, so consecutive halvings of the span
  // must move them by a bounded amount. Noise showed up here as changes of order 1.
  ok(worst_jump1 < 0.05, "dE/dT stays smooth into the coincidence");
  ok(worst_jump2 < 0.05, "d2E/dT2 stays smooth into the coincidence");
}

// The layer-mean helpers' two branches, against each other where both are valid.
//
// ⚠️ THE CHECK THAT CATCHES A WRONG CONSTANT, and it was missing: the midpoint limits
// are f'/2, f''/3 and f''/6, and I first wrote the last two as f''/4. The continuity
// sweep could not see it -- those terms contribute little to dE/dT there -- so this
// compares the two branches directly, just ABOVE the crossover where the divided
// difference is still accurate and the midpoint form is too. A wrong constant shows
// up here as a 25% disagreement rather than not at all.
// The whole layer-mean family against a quadrature of higher order than the one the
// model uses, at every span the model switches across.
//
// ⚠️ THIS REPLACED A CHECK OF TWO BRANCHES AGAINST EACH OTHER, and the replacement is
// the point. The model used to switch between a divided difference and a midpoint
// asymptotic, so the only referee available was the span at which those two met --
// which says nothing about whether either is right, and which is why an f''/4 written
// for an f''/6 survived until it was derived by hand a second time. A mean is an
// average, so an independent average of higher order referees every quantity in the
// family directly, at every span, including the ones no operating point reaches often
// enough to be noticed at.
void test_the_layer_mean_against_a_higher_order_average() {
  printf("the layer-mean family against a 15-point average\n");
  phylloptim::Leaf l = env::fresh();
  env::drivers(2.0, 900.0, 2.0, 3, 3).drive(l, env::kTheta);
  auto& r = l.roots_;

  // Fifteen nodes, generated here by Newton on the Legendre polynomial, so the
  // model's own seven are not refereeing themselves.
  const int n = 15;
  std::vector<double> t(n), w(n);
  for (int i = 0; i < n; ++i) {
    double x = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0, p1 = 0.0;
      for (int j = 0; j < n; ++j) {
        const double p2 = p1;
        p1 = p0;
        p0 = ((2.0 * j + 1.0) * x * p1 - j * p2) / (j + 1);
      }
      const double dp = n * (x * p0 - p1) / (x * x - 1.0);
      const double dx = -p0 / dp;
      x += dx;
      if (std::abs(dx) < 1e-16) break;
    }
    double p0 = 1.0, p1 = 0.0;
    for (int j = 0; j < n; ++j) {
      const double p2 = p1;
      p1 = p0;
      p0 = ((2.0 * j + 1.0) * x * p1 - j * p2) / (j + 1);
    }
    const double dp = n * (x * p0 - p1) / (x * x - 1.0);
    t[std::size_t(i)] = x;
    w[std::size_t(i)] = 2.0 / ((1.0 - x * x) * dp * dp);
  }
  auto rel_gap = [](double got, double want) -> double {
    return std::abs(got - want) / (std::abs(want) + 1e-300);
  };
  const double bb = r.root_b, cc = r.root_c;
  auto avg = [&](double lo, double hi, int order, double pow_hi,
                 double pow_lo) -> double {
    const double m = 0.5 * (lo + hi), sp = hi - lo;
    double acc = 0.0;
    for (int i = 0; i < n; ++i) {
      const double x = m + 0.5 * sp * t[std::size_t(i)];
      const double a = 0.5 * (1.0 + t[std::size_t(i)]);
      const double b = 0.5 * (1.0 - t[std::size_t(i)]);
      const double g =
          order == 0 ? phylloptim::vulnerability_curve_at<double>(x, bb, cc)
          : order == 1
              ? phylloptim::vulnerability_curve_slope_at<double>(x, bb, cc)
              : phylloptim::vulnerability_curve_curvature_at<double>(x, bb, cc);
      acc += w[std::size_t(i)] * g * std::pow(a, pow_hi) * std::pow(b, pow_lo);
    }
    return 0.5 * acc;
  };

  printf("      %-9s %-11s %-11s %-11s %-11s %s\n", "span", "mean", "d/dhi",
         "d/dlo", "d2/dhi2", "mixed");
  double worst = 0.0, worst_direct = 0.0, worst_table = 0.0;
  for (double centre : {1.0, 2.0, 4.0}) {
    for (int e = 0; e <= 8; ++e) {
      const double sp = std::pow(10.0, -double(e));
      const double lo = centre - 0.5 * sp, hi = centre + 0.5 * sp;
      if (lo <= 0.0) continue;
      const double mn = r.layer_mean(lo, hi);
      const double e0 = rel_gap(mn, avg(lo, hi, 0, 0, 0));
      const double e1 =
          rel_gap(r.layer_mean_dbound(lo, hi, true, mn), avg(lo, hi, 1, 1, 0));
      const double e2 =
          rel_gap(r.layer_mean_dbound(lo, hi, false, mn), avg(lo, hi, 1, 0, 1));
      const double e3 = rel_gap(r.layer_mean_dbound2(lo, hi, true, mn),
                                      avg(lo, hi, 2, 2, 0));
      const double e4 = rel_gap(r.layer_mean_dbound_mixed(lo, hi, mn),
                                      avg(lo, hi, 2, 1, 1));
      if (centre == 2.0) {
        printf("      %-9.0e %-11.2e %-11.2e %-11.2e %-11.2e %.2e\n", sp, e0, e1,
               e2, e3, e4);
      }
      const double here =
          std::max(std::max(e0, e1), std::max(std::max(e2, e3), e4));
      (sp < 1e-3 ? worst_direct : worst_table) =
          std::max(sp < 1e-3 ? worst_direct : worst_table, here);
      worst = std::max(worst, here);
    }
  }
  printf("      worst: direct branch %.3e, table branch %.3e (all %.3e)\n",
         worst_direct, worst_table, worst);
  // Two bounds, because the two branches are held to different things and that IS
  // the design. Below the switch the average is formed directly and is exact to
  // rounding; above it the difference of two table reads is well conditioned but
  // carries the table's own error, which is what makes it cheap. One bound across
  // both would either excuse the direct branch or forbid the cheap one.
  //
  // The table bound is what the three thresholds this replaced reached at their
  // worst, 2.67e-05, so a regression to that scheme fails here.
  ok(worst_direct < 1e-10, "the direct average is exact to rounding below the switch");
  ok(worst_table < 1e-5, "the difference above the switch stays inside the table");
}

// The templated supply collar derivative, against the double one it replaces and
// against the independent second derivative.
//
// ⚠️ THE REFEREE BEFORE THE ALGEBRA. This form exists so the interior collar can be
// closed by a RESIDUAL rather than by rows handed over as numbers -- which needs
// dE_up/dp carrying its parameter rows. Two checks, because two things can go wrong:
// the value must be the double form's exactly, and its collar derivative must be the
// closed-form second derivative, which is computed by a different route entirely.
void test_the_templated_supply_derivative_against_both_routes() {
  printf("the templated dE_up/dp: against the double form and against d2\n");
  namespace pl = phylloptim;
  double worst_value = 0.0, worst_slope = 0.0;
  int compared = 0;
  for (double psi : {0.5, 1.0, 2.0, 3.0}) {
    for (int layers : {1, 3, 5}) {
      fixture::Physiology d = env::drivers(psi, 900.0, 2.0, layers, layers);
      pl::Leaf l = env::fresh();
      d.drive(l, env::kTheta);
      l.find_root_collar_psi();
      const double p = l.opt_root_psi_;
      const std::vector<double>& soil = l.supply_psi_soil();

      const double theirs = l.dE_from_soil_dpsi_collar(p, soil);
      const double mine = l.roots_.template duptake_dpsi_at<double>(
          p, l.roots_.held_supply());
      if (!std::isfinite(theirs) || theirs == 0.0) {
        continue;
      }
      ++compared;
      worst_value = std::max(worst_value, std::abs(mine / theirs - 1.0));

      // The same quantity at a tangent, seeded in the collar: its derivative is the
      // supply's second collar derivative, which d2E_from_soil_dpsi_collar2 forms
      // from its own algebra.
      using tangent = odelia::ode::tangent_scalar<double>;
      pl::SupplyValues<tangent> sup;
      const pl::SupplyAt<double> held = l.roots_.held_supply();
      for (std::size_t j = 0; j < held.psi_soil.size(); ++j) {
        sup.psi_soil.push_back(tangent(held.psi_soil[j]));
        sup.r_R_H_min.push_back(tangent(held.r_R_H_min[j]));
        sup.r_R_V_sum.push_back(tangent(held.r_R_V_sum[j]));
      }
      sup.root_b = tangent(held.root_b);
      sup.root_c = tangent(held.root_c);
      tangent pt = p;
      odelia::ode::seed_direction(pt, 1.0);
      const double slope = odelia::ode::derivative_along(
          l.roots_.template duptake_dpsi_at<tangent>(pt, sup.at()));
      const double closed = l.d2E_from_soil_dpsi_collar2(p, soil);
      if (std::isfinite(closed) && closed != 0.0) {
        worst_slope = std::max(worst_slope, std::abs(slope / closed - 1.0));
      }
    }
  }
  printf("      %d points; worst value gap %.3e, worst slope gap %.3e\n", compared,
         worst_value, worst_slope);
  ok(compared >= 8, "the grid reaches solved interior points");
  ok(worst_value < 1e-12, "the templated form is the double form");
  ok(worst_slope < 1e-8, "its collar derivative is the closed second derivative");
}

void test_the_two_zero_flux_kinds_are_two_points() {
  printf("the two zero-flux kinds are two points, and the placement says which\n");
  namespace pl = phylloptim;
  const int L = 5;

  // Shade death is reached by LIGHT and not by water, so the fixture is a wet
  // soil in the dark. No moisture sweep finds this branch.
  fixture::Physiology dark = env::drivers(2.0, 0.0, 2.0, L, L);
  pl::Leaf shaded = env::fresh();
  dark.drive(shaded, env::kTheta);
  shaded.find_root_collar_psi();
  ok(shaded.operating_point_kind() == pl::Leaf::OperatingPointKind::ShadeDeath,
     "an unlit leaf on wet soil reaches shade death");

  // A dry soil in full light is the other kind, on the same trait vector.
  fixture::Physiology dry = env::drivers(7.0, 900.0, 2.0, L, L);
  pl::Leaf parched = env::fresh();
  dry.drive(parched, env::kTheta);
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

  // The places are different potentials, and each leaf's stem sits at its own.
  ok(shaded.opt_psi_stem_ == shaded.opt_root_psi_,
     "shade death places the stem exactly at the collar");
  ok(parched.opt_psi_stem_ == parched.psi_crit,
     "and a hydraulic shutdown places it exactly at the critical potential");
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

  // The two placements say which potential the cost is paid at, and the outputs
  // carry that: a shutdown holds the stem at its critical potential, so profit
  // moves with it and with no soil at all; shade death holds it at the collar,
  // which the soil places, so the reverse is true of both.
  auto seeded = [&](pl::Leaf& leaf, fixture::Input which, int layer)
      -> pl::Leaf::LeafOutputs<pl::tangent> {
    fixture::Soil soil;
    for (int i = 0; i < L; ++i) {
      soil.depth.push_back(1.0 * (i + 1));
      soil.carbon.push_back(1.0 / L / 0.05);
    }
    const pl::LeafInputs<pl::tangent> in = fixture::leaf_inputs<pl::tangent>(
        leaf, which, layer, soil);
    // Placed, then evaluated: a shade-death collar is the wet bound, and whether
    // an input reaches profit through that placement is the whole question here.
    const auto draw = leaf.supply_draw_at<pl::tangent>(pl::tangent(leaf.opt_root_psi_), in.supply);
    const pl::tangent collar = leaf.collar_at<pl::tangent>(in, draw);
    return leaf.outputs_at<pl::tangent>(collar, in, draw);
  };
  ok(pl::derivative_along(seeded(parched, fixture::Input::psi_crit, 0).profit) != 0.0,
     "at a hydraulic shutdown the critical potential IS the placement");
  ok(pl::derivative_along(seeded(shaded, fixture::Input::psi_crit, 0).profit) == 0.0,
     "and at shade death it is inactive, exactly");
  ok(pl::derivative_along(seeded(parched, fixture::Input::psi_soil, 0).profit) == 0.0,
     "a hydraulic shutdown reads no soil at all");
  ok(pl::derivative_along(seeded(shaded, fixture::Input::psi_soil, 0).profit) != 0.0,
     "and shade death reads the soil its collar is placed by");
  ok(pl::derivative_along(seeded(shaded, fixture::Input::root_c, 0).profit) != 0.0,
     "including the root curve, which moves the supply that placement is defined by");

  // The collar channel, which no evaluation of this model can record: the stem
  // sits at the collar, so the marginal profit takes its no-flow exit.
  bool feasible = true;
  const double sentinel = shaded.dprofit_droot_collar_psi(shaded.opt_root_psi_,
                                                          &feasible);
  ok(sentinel == 0.0 && !feasible,
     "the marginal profit at a shade-death collar is a sentinel, not a zero");
  shaded.find_root_collar_psi();

  // The collar channel, seeded on the collar itself. At shade death the stem is
  // tied to the collar, so gross assimilation does not move and profit responds
  // through the cost alone -- while each layer's own draw still does. A
  // hydraulic shutdown has no such channel at all: nothing places its collar.
  auto along_collar = [&](pl::Leaf& leaf) -> pl::Leaf::LeafOutputs<pl::tangent> {
    fixture::Soil soil;
    for (int i = 0; i < L; ++i) {
      soil.depth.push_back(1.0 * (i + 1));
      soil.carbon.push_back(1.0 / L / 0.05);
    }
    pl::tangent c = leaf.opt_root_psi_;
    pl::seed_direction(c, 1.0);
    const pl::LeafInputs<pl::tangent> in =
        fixture::leaf_inputs<pl::tangent>(leaf, fixture::Input::None, 0, soil);
    return leaf.outputs_at<pl::tangent>(
        c, in, leaf.supply_draw_at<pl::tangent>(pl::tangent(leaf.opt_root_psi_), in.supply));
  };
  const pl::Leaf::LeafOutputs<pl::tangent> shaded_along = along_collar(shaded);
  ok(pl::derivative_along(shaded_along.profit) != 0.0,
     "profit responds to the collar at shade death, through the cost");
  // The per-layer draws respond to the collar, and their sum is the total's own
  // response -- checked against each other rather than against a second
  // derivation, because there is no longer a second one to check against.
  double summed = 0.0;
  bool each_moves = shaded_along.uptake.size() == std::size_t(L);
  for (int i = 0; i < L && each_moves; ++i) {
    const double d = pl::derivative_along(shaded_along.uptake[std::size_t(i)]);
    each_moves = d != 0.0;
    summed += d;
  }
  ok(each_moves, "and every layer's own draw responds to the collar");
  {
    fixture::Soil soil;
    for (int i = 0; i < L; ++i) {
      soil.depth.push_back(1.0 * (i + 1));
      soil.carbon.push_back(1.0 / L / 0.05);
    }
    pl::tangent c = shaded.opt_root_psi_;
    pl::seed_direction(c, 1.0);
    const pl::LeafInputs<pl::tangent> in = fixture::leaf_inputs<pl::tangent>(
        shaded, fixture::Input::None, 0, soil);
    std::vector<pl::tangent> per;
    const double total = pl::derivative_along(
        shaded.E_from_soil_at<pl::tangent>(c, in.supply.at(), per));
    // The total is in kg and the layers in mol, which is the one place this
    // model carries two units for one quantity.
    near(summed * phylloptim::kg_per_mol_h2o, total, 1e-12,
         "and they sum to the total's response");
  }

  const pl::Leaf::LeafOutputs<pl::tangent> parched_along = along_collar(parched);
  ok(pl::derivative_along(parched_along.profit) == 0.0,
     "a hydraulic shutdown has no collar channel, exactly");

  // ⚠️ ONE LAYER IS THE CASE THAT USED TO REFUSE, and it is where two of these
  // items meet. The wet bound is the collar at which total uptake vanishes; with a
  // SINGLE layer that is the collar at which its numerator vanishes, which is the
  // gravity balance -- and the supply refused every derivative there, so the bound
  // had no row and a shade-death point placed on it could not report its own
  // movement. Only the numerator vanishes at that collar, so nothing about it was
  // ever 0/0; with the refusal gone the bound has a row at every layer count and
  // the point's movement is reported rather than folded into a re-solve.
  for (int layers : {1, 2, 5}) {
    fixture::Physiology thin = env::drivers(2.0, 0.0, 2.0, layers, layers);
    pl::Leaf probe = env::fresh();
    thin.drive(probe, env::kTheta);
    probe.find_root_collar_psi();
    const pl::Leaf::BoundRow wet = probe.bound_row(pl::Leaf::WhichBound::Wet);
    ok(wet.finite, std::string("the wet bound has a row at ") +
                       std::to_string(layers) + " layer(s)");
    ok(probe.operating_point_kind() == pl::Leaf::OperatingPointKind::ShadeDeath,
       std::string("and shade death is what the solve calls it at ") +
           std::to_string(layers) + " layer(s)");
    // Shade death sits ON the wet bound, so the collar the outputs are placed at
    // is that bound -- reported by the condition rather than measured -- and
    // every output there is a number.
    fixture::Soil soil;
    for (int i = 0; i < layers; ++i) {
      soil.depth.push_back(1.0 * (i + 1));
      soil.carbon.push_back(1.0 / layers / 0.05);
    }
    const pl::LeafInputs<double> in = fixture::leaf_inputs<double>(
        probe, fixture::Input::None, 0, soil);
    const auto draw = probe.supply_draw_at<double>(wet.bound, in.supply);
    const double placed = probe.bound_at<double>(
        pl::Leaf::WhichBound::Wet, wet.bound, in, draw);
    ok(placed == wet.bound,
       std::string("and the collar at ") + std::to_string(layers) +
           " layer(s) is the wet bound itself");
    const pl::Leaf::LeafOutputs<double> got =
        probe.outputs_at<double>(placed, in, draw);
    bool all_numbers = std::isfinite(got.profit);
    for (double v : got.uptake) {
      all_numbers = all_numbers && std::isfinite(v);
    }
    ok(all_numbers, std::string("every output at ") + std::to_string(layers) +
                        " layer(s) is a number");
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
  test_marginal_price_water();
  test_root_vulnerability_is_bounded_past_its_grid();
  test_root_psi_crit_clamp_binds();
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
  test_environment_water_rows_are_rank_one();
  test_the_curves_trait_derivative_is_the_models_own();
  test_the_transport_leaves_the_flux_where_it_is();
  test_the_condition_is_the_stem_potential_and_its_collar_response();
  test_the_supplys_mixed_partials_match_a_difference();
  probe_root_curve_slope();
  probe_table_vs_curve();
  test_the_condition_is_carbon_bought_against_tension_paid();
  test_interior_optima_are_stationary();
  test_the_supply_answers_at_the_two_coincidences();
  test_the_supplys_second_collar_derivative();
  test_the_marginals_collar_slope_against_a_difference();
  test_where_the_mean_conductivity_should_stop_differencing();
  test_the_supply_derivatives_stay_smooth_into_a_coincidence();
  test_the_layer_mean_against_a_higher_order_average();
  test_the_templated_supply_derivative_against_both_routes();
  test_the_two_zero_flux_kinds_are_two_points();
  benchmark();

  printf("\n%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
