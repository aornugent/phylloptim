// PROBE (not a test, not shipped): the transport family, and the question its
// hand-written row rests on.
//
// `transport_trait_rows` never differentiates the inverse table the stem
// potential's VALUE comes from. It applies the theorem to the forward relation
//
//     kappa * (G(sigma) - G(p)) - E_up  =  0
//
// and takes dG/dsigma from the closed-form curve f rather than from the table --
// `transport_slope`'s comment says so outright, and the header above
// `stem_curve_integral_inverse` records that polishing the inverse was tried and
// abandoned BECAUSE every analytic derivative here was derived for the unpolished
// map. Value and derivative are two implementations kept consistent by hand.
//
// Two measurements, and the first decides the second.
//
//   1. The table's own slope against the closed form it is supposed to carry.
//      A generic pass reads whatever the table says; the hand row reads f. If
//      those disagree materially, a generic design must be told which slope a
//      table carries -- which is a one-line decision per table, and odelia's
//      `with_query_derivative(value, dydu, u, up)` already takes `dydu` as a
//      parameter for exactly this reason.
//
//   2. The generic row against the hand row. One residual, one application of
//      the theorem, no knowledge of which input is seeded.
//
//   make CXX=g++ probe_transport_tangent && ./probe_transport_tangent

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

// --- 1. what slope does the table carry ---------------------------------------
//
// G is tabulated because it has no cheap closed form; its derivative f does have
// one, and the quintic is built with f as an exact channel. So the two should
// agree to the interpolation's own order, and the question is what that is.
void slope_agreement(const Leaf& l) {
  const auto& G = l.transpiration_from_psi;
  const double lo = G.min(), hi = G.max();
  double worst = 0.0, worst_at = 0.0, scale = 0.0;
  int n = 0;
  // Off-knot on purpose: at a knot a Hermite table reproduces the slope it was
  // given, so sampling there would measure the build and not the interpolation.
  for (int i = 0; i <= 2000; ++i) {
    const double psi = lo + (hi - lo) * (i + 0.37) / 2001.0;
    if (psi >= hi) break;
    const double from_table = G.slope(psi);
    const double closed = l.proportion_of_conductivity(psi);
    if (std::abs(closed) > scale) scale = std::abs(closed);
    const double d = std::abs(from_table - closed);
    if (d > worst) { worst = d; worst_at = psi; }
    ++n;
  }
  std::printf("1. dG/dpsi: the table's own slope against the closed-form curve\n");
  std::printf("   %d off-knot samples over [%.4g, %.4g]\n", n, lo, hi);
  std::printf("   worst absolute %.3e at psi %.6g; curve peaks at %.4g, so "
              "relative %.3e\n\n",
              worst, worst_at, scale, scale > 0.0 ? worst / scale : 0.0);
}

// --- 2. the transport row, generically ----------------------------------------
//
// The relation the potential satisfies, at one scalar. E_up and the collar are
// held: no stem property reaches the soil-to-collar path, so a transport trait
// moves the potential and nothing upstream of it.
//
// `slope_of_G` is where the decision from measurement 1 lands -- one argument,
// not a per-trait derivation.
template <typename T>
T transport_residual(const Leaf& l, T sigma, T kappa, double p, double E_up) {
  const double sv = odelia::util::to_passive(sigma);
  const T G_sigma = T(l.transpiration_from_psi.eval(sv)) +
                    l.proportion_of_conductivity(sv) * (sigma - sv);
  const double G_p = l.transpiration_from_psi.eval(p);
  return kappa * (G_sigma - T(G_p)) - T(E_up);
}

// d(psi_stem)/d(kappa) at a held collar, from the theorem applied once.
double dpsistem_generic(const Leaf& l) {
  const double sigma = l.opt_psi_stem_;
  const double p = l.opt_root_psi_;
  const double kappa = l.leaf_specific_conductance_max_;
  const double E_up = l.transpiration(sigma, p);

  tangent s_seed = sigma;
  seed_direction(s_seed, 1.0);
  const double R_sigma = derivative_along(
      transport_residual(l, s_seed, tangent(kappa), p, E_up));

  tangent k_seed = kappa;
  seed_direction(k_seed, 1.0);
  const double R_kappa = derivative_along(
      transport_residual(l, tangent(sigma), k_seed, p, E_up));

  return -R_kappa / R_sigma;
}

double one_state(double psi_soil, double ppfd, int layers, int& compared,
                 bool& interior) {
  const grad::Drivers d = drivers(psi_soil, ppfd, layers);
  Leaf l;
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  try {
    l.find_root_collar_psi();
  } catch (const std::exception& e) {
    std::printf("   %-8.3g %-8.4g %d   solve failed\n", psi_soil, ppfd, layers);
    return 0.0;
  }
  interior = l.operating_point_kind() == Leaf::OperatingPointKind::Interior;

  Leaf::TransportTraitRows ref;
  if (!l.transport_trait_rows(Leaf::TransportTrait::Conductance, 0.0, ref)) {
    std::printf("   %-8.3g %-8.4g %d   %-10s row refused\n", psi_soil, ppfd,
                layers, Leaf::operating_point_kind_name(l.operating_point_kind()));
    return 0.0;
  }
  const double got = dpsistem_generic(l);
  ++compared;
  const double rel = ref.dpsistem == 0.0
                         ? std::abs(got)
                         : std::abs(got - ref.dpsistem) / std::abs(ref.dpsistem);
  std::printf("   %-8.3g %-8.4g %d   %-10s %16.9e %16.9e %11.3e\n", psi_soil,
              ppfd, layers,
              Leaf::operating_point_kind_name(l.operating_point_kind()),
              ref.dpsistem, got, rel);
  return rel;
}

}  // namespace

int main() {
  {
    const grad::Drivers d = drivers(1.5, 1000.0, 3);
    Leaf l;
    const grad::Settings s;
    grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
    slope_agreement(l);
  }

  std::printf("2. d(psi_stem)/d(kappa): transport_trait_rows against one "
              "generic pass\n\n");
  std::printf("   %-8s %-8s %s   %-10s %16s %16s %11s\n", "psi_soil", "PPFD",
              "L", "point", "hand", "generic", "rel");

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
  std::printf("\n   %d rows compared.\n", compared);
  std::printf("   interior: worst %.3e\n", worst_interior);
  std::printf("   pinned  : worst %.3e\n", worst_pinned);
  return 0;
}
