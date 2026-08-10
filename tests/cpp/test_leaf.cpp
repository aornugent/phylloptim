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
#include <cstdio>
#include <limits>
#include <string>
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

// PLAN 11b: the AD derivative and the forward model are now instantiations of ONE
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
// The count check over the whole 288-point grid lives in test_golden.cpp; this is
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
               l.curv_fact_elec_trans, l.curv_fact_colim, l.cost_scale_TF24);
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
         m.operating_point_kind() != Kind::PinnedDry,
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
  const double sliver = 1e-9;
  bool feasible = true;
  s.dprofit_at_collar_psi(sa + 1e-6 * sliver, &feasible);
  ok(!feasible, "the wet bound admits no informative gradient");
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
                        100, 1e-3, 1000, 7.5);
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
  near(base.rd_to_vcmax_ratio_, 0.015, 1e-12, "rd_to_vcmax_ratio_ default");

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

  // Respiration fraction: doubling it must lower assimilation.
  {
    phylloptim::Leaf r2 = make_leaf(d, {2.0}, {1.0});
    r2.rd_to_vcmax_ratio_ = 0.030;
    r2.update_temperature_dependent_params(d.leaf_temp);
    r2.find_root_collar_psi();
    ok(r2.assim_colimited_ < A_base,
       "doubling the respiration fraction lowers assimilation");
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
                       5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5);
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
                        5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5);
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
                     t[10], t[11], 1e-3, 100, 1e-3, 1000, t[12]);
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
                      t[10], t[11], t[12]);
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
                 5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5);
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
                 5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5);
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
                     t[10], t[11], t[12]);
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

int main() {
  test_defaults_are_unset();
  test_vulnerability_curve();
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
  test_soil_conductance_is_positive();
  test_root_psi_crit_clamp_binds();
  test_signed_potentials_are_rejected();
  test_lambda_equals_dA_dE_single_layer();
  test_multilayer_lambda_identity();
  test_g1_eff();
  test_energy_balance_path_runs();
  test_pm_wind_speed_validation();
  test_pm_leaf_temperature_response();
  test_closed_form();
  test_single_potential();
  test_leaf_on_single_potential();
  test_root_network_from_carbon();
  test_temperature_parameters_are_settable();
  test_set_traits_matches_a_fresh_leaf();
  test_perturb_stem_b_matches_a_rebuild();
  test_bad_input_throws();
  benchmark();

  printf("\n%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
