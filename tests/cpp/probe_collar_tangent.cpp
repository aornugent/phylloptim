// PROBE (not a test, not shipped): the marginal profit at the collar, which is
// the condition the interior operating point solves.
//
// `dprofit_at_collar_psi` is ninety-odd lines in which every step but two is a
// hand-written chain rule or an implicit-function quotient -- dgc_dpsistem,
// dgc_dpsi, g_ci, dci_dpsistem, dci_dpsi_expl, dci_dpsi, base. The two that are
// not are A'(ci) and C'(psi_stem), which are already tangents of the model's own
// kernels, and its comment says why: "derivatives of the function actually
// evaluated rather than of a hand-kept mirror of it". This asks what is left if
// that reasoning is applied to the whole composition rather than to two terms in
// it.
//
// The claim under test is that the leaf has ONE evaluation at a held collar, and
// that every derivative it must supply is that evaluation seeded in a different
// direction. probe_carbon_tangent seeds a trait; this seeds the collar. Same two
// residuals, same kernels, same theorem.
//
//   T1   kappa*(G(sigma) - G(p)) - E_up(p) = 0     places the stem potential
//   T2   A(ci)*k1 - gc(sigma,p)*(ca - ci)*k2 = 0   places the intercellular CO2
//        profit = A(ci) - C(sigma)
//
// A table read is lifted with the slope named beside it, which is the decision
// probe_transport_tangent priced: G carries its value because it has no cheap
// closed form, and every derivative of it is taken from the curve f exactly.
//
//   make CXX=g++ probe_collar_tangent && ./probe_collar_tangent

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace grad = phylloptim::gradient;
using phylloptim::Leaf;
using phylloptim::tangent;
using phylloptim::seed_direction;
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

// A tabulated read on the tangent scalar: the value at the query's passive part,
// carrying the query's own derivative through the slope named here. This is
// odelia's `with_query_derivative` with the knots left in double, which is what
// a table nothing seeds needs.
tangent lifted(double value, double slope, const tangent& u) {
  const double up = odelia::util::to_passive(u);
  return tangent(value) + slope * (u - tangent(up));
}

// d(profit)/d(collar) at the solved point, from the two residuals and nothing
// else. `l` must have solved: the values are read off it and only the
// derivatives are formed here.
double dprofit_dcollar_generic(Leaf& l, double p) {
  // The point at this collar, in double. Both root-finds run off the tape --
  // that is the refusal the whole design rests on -- and only their derivatives
  // are formed below.
  const double sigma = l.find_psi_stem_from_psi_root(p, l.supply_psi_soil());
  const double ci = l.psi_stem_to_ci(sigma, p);
  const double kappa = l.leaf_specific_conductance_max_;

  tangent p_ad = p;
  seed_direction(p_ad, 1.0);

  // T1. G is tabulated and its derivative is the curve, exactly.
  const auto& G = l.transpiration_from_psi;
  auto G_at = [&](const tangent& x) -> tangent {
    const double xv = odelia::util::to_passive(x);
    return lifted(G.eval(xv), l.proportion_of_conductivity(xv), x);
  };
  // The soil-to-collar uptake, with the conductance the roots model supplies.
  const tangent E_up =
      lifted(l.transpiration(sigma, p),
             l.dE_from_soil_dpsi_collar(p, l.supply_psi_soil()), p_ad);
  (void)0;
  // sigma solves T1; the theorem's denominator is the residual's own slope in it.
  const double R1_sigma = kappa * l.proportion_of_conductivity(sigma);
  const tangent R1_at_held = kappa * (tangent(G.eval(sigma)) - G_at(p_ad)) - E_up;
  tangent sigma_ad = sigma;
  seed_direction(sigma_ad, -derivative_along(R1_at_held) / R1_sigma);

  // T2. The conductance is the stem integral BETWEEN the two potentials, so the
  // collar reaches it directly as well as through sigma -- carried here by both
  // being live rather than by a second coefficient.
  const double gc_const =
      l.atm_kpa_ * phylloptim::kg_to_mol_h2o / l.atm_vpd_ /
      phylloptim::H2O_CO2_stom_diff_ratio;
  const tangent gc = gc_const * kappa * (G_at(sigma_ad) - G_at(p_ad));
  const double inv_atm = 1.0 / (l.atm_kpa_ * phylloptim::kPa_to_Pa);

  auto R2 = [&](const tangent& c) -> tangent {
    return l.assim_colimited_kernel(c) * phylloptim::umol_to_mol -
           gc * (tangent(l.ca_) - c) * inv_atm;
  };
  // Denominator in double, with everything upstream passive.
  tangent ci_seed = ci;
  seed_direction(ci_seed, 1.0);
  const double R2_ci =
      derivative_along(l.assim_colimited_kernel(ci_seed)) *
          phylloptim::umol_to_mol +
      odelia::util::to_passive(gc) * inv_atm;
  tangent ci_ad = ci;
  seed_direction(ci_ad, -derivative_along(R2(tangent(ci))) / R2_ci);

  const tangent profit =
      l.assim_colimited_kernel(ci_ad) - l.hydraulic_cost_TF_kernel(sigma_ad);
  return derivative_along(profit);
}

double one_state(double psi_soil, double ppfd, int layers, int& compared,
                 bool& interior) {
  const grad::Drivers d = drivers(psi_soil, ppfd, layers);
  Leaf l;
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  try {
    l.find_root_collar_psi();
  } catch (const std::exception&) {
    std::printf("   %-8.3g %-8.4g %d   solve failed\n", psi_soil, ppfd, layers);
    return 0.0;
  }
  interior = l.operating_point_kind() == Leaf::OperatingPointKind::Interior;

  // ⚠️ NOT AT THE OPTIMUM. There this quantity IS zero -- it is the stationarity
  // condition the solve drove to nothing -- so the two forms agree at machine
  // noise whether or not they agree about anything. Taken across the feasible
  // interval instead, which is also where the optimiser calls it: a bracketing
  // solve evaluates the marginal at a sequence of collars and reads its sign.
  double bound_a = 0.0, bound_b = 0.0;
  if (!l.prepare_collar_solve(bound_a, bound_b)) {
    std::printf("   %-8.3g %-8.4g %d   %-10s point already determined\n",
                psi_soil, ppfd, layers,
                Leaf::operating_point_kind_name(l.operating_point_kind()));
    return 0.0;
  }

  double worst = 0.0;
  int at = 0;
  double worst_hand = 0.0, worst_got = 0.0;
  for (double frac : {0.05, 0.2, 0.4, 0.6, 0.8, 0.95}) {
    const double p = bound_a + frac * (bound_b - bound_a);
    bool feasible = false;
    const double hand = l.dprofit_at_collar_psi(p, &feasible);
    if (!feasible || !std::isfinite(hand)) {
      continue;
    }
    const double got = dprofit_dcollar_generic(l, p);
    ++compared;
    ++at;
    // Relative to the marginal itself, which is a real number away from the
    // optimum -- so this is the ordinary measure and needs no special scale.
    const double rel = std::abs(hand) > 0.0
                           ? std::abs(got - hand) / std::abs(hand)
                           : std::abs(got);
    if (rel > worst) {
      worst = rel;
      worst_hand = hand;
      worst_got = got;
    }
  }
  std::printf("   %-8.3g %-8.4g %d   %-10s %2d %15.7e %15.7e %11.3e\n",
              psi_soil, ppfd, layers,
              Leaf::operating_point_kind_name(l.operating_point_kind()), at,
              worst_hand, worst_got, worst);
  return worst;
}

}  // namespace

int main() {
  std::printf("d(profit)/d(collar): dprofit_at_collar_psi against two residuals "
              "and the kernels\n\n");
  std::printf("   %-8s %-8s %s   %-10s %2s %15s %15s %11s\n", "psi_soil",
              "PPFD", "L", "point", "n", "hand", "generic", "worst rel");

  double worst_interior = 0.0, worst_pinned = 0.0;
  int compared = 0, pinned = 0;
  for (double psi_soil : {0.2, 0.8, 1.5, 2.5, 3.5, 4.5}) {
    for (double ppfd : {30.0, 300.0, 1000.0, 2000.0}) {
      for (int layers : {1, 3}) {
        bool interior = true;
        const double w = one_state(psi_soil, ppfd, layers, compared, interior);
        double& into = interior ? worst_interior : worst_pinned;
        if (!interior) ++pinned;
        if (w > into) into = w;
      }
    }
  }
  std::printf("\n   %d collars compared over 48 states.\n", compared);
  std::printf("   interior: worst %.3e\n", worst_interior);
  std::printf("   pinned  : worst %.3e\n", worst_pinned);
  return 0;
}
