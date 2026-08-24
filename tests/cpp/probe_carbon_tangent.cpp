// PROBE (not a test, not shipped): does one generic tangent pass reproduce the
// carbon-side profit rows that `photo_trait_rows` assembles by hand?
//
// `photo_trait_rows` composes the leaf's own scalar-generic kernels and then
// differentiates the composition itself: four second-order passes for A's mixed
// partials, four first-order passes for the transport, an implicit-function
// quotient written out, and one hand-threaded chain rule per trait. Every kernel
// in it is already generic; only the composition is not.
//
// This asks whether the composition can be generic too. It evaluates profit at
// the held collar on the tangent scalar with one input seeded, applies the
// theorem once at the ci root-find, and reads the derivative off the result --
// no per-trait threading and no second-order pass.
//
//   make CXX=g++ probe_carbon_tangent && ./probe_carbon_tangent

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

bool util_finite(double x) { return phylloptim::util::is_finite(x); }

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

// The inputs that reach profit through assimilation at a held collar, at one
// scalar. The two _25 traits are carried as the temperature-adjusted value they
// produce, seeded in the direction their ratio gives.
template <typename T>
struct carbon_inputs {
  T vcmax, transport_jmax, quantum_yield, curv_elec, curv_colim, ppfd,
      respiration;
};

// The condition that places the intercellular concentration: assimilation demand
// against the stomatal supply at a conductance the collar fixes.
template <typename T>
T ci_residual(const Leaf& l, T ci, const carbon_inputs<T>& p, double gc,
              double ca, double atm_kpa) {
  const T J = l.electron_transport_kernel(p.ppfd, p.quantum_yield, p.curv_elec,
                                          p.transport_jmax);
  const T A =
      l.assim_colimited_kernel(ci, p.vcmax, J, p.curv_colim, p.respiration);
  return A * phylloptim::umol_to_mol -
         gc * (T(ca) - ci) / (atm_kpa * phylloptim::kPa_to_Pa);
}

// Profit at the held collar, differentiated in one seeded input. `which` names
// the input to seed; everything else rides passive. The point is that nothing
// below knows which input that is.
double dprofit_generic(const Leaf& l, int which) {
  const double ci = l.ci_;
  const double gc = l.stom_cond_CO2_;

  auto inputs = [&](int at) {
    carbon_inputs<tangent> p{tangent(l.vcmax_),
                             tangent(l.jmax_),
                             tangent(l.a),
                             tangent(l.curv_fact_elec_trans),
                             tangent(l.curv_fact_colim),
                             tangent(l.PPFD_),
                             tangent(l.R_d_)};
    // The chain from a _25 trait to its temperature-adjusted value is their
    // ratio, so seeding the adjusted value in that direction seeds the trait.
    switch (at) {
    case 0: seed_direction(p.vcmax, l.vcmax_ / l.vcmax_25); break;
    case 1: seed_direction(p.transport_jmax, l.jmax_ / l.jmax_25); break;
    case 2: seed_direction(p.quantum_yield, 1.0); break;
    case 3: seed_direction(p.curv_elec, 1.0); break;
    case 4: seed_direction(p.curv_colim, 1.0); break;
    case 5: seed_direction(p.ppfd, 1.0); break;
    // Dark respiration is the fourth trait the kernel now carries. The hand
    // assembly reaches it by a special case -- "it shifts A by exactly -1 per
    // unit and leaves A' untouched. No pass needed" -- which is a derivation
    // rather than a reading, and is what a parameter the kernel does not take
    // costs. Seeded like any other here.
    case 6: seed_direction(p.respiration, l.R_d_ / l.R_d_25); break;
    default: break;
    }
    return p;
  };

  const carbon_inputs<tangent> p = inputs(which);

  // ⚠️ FOLLOW THE BRANCH THE FORWARD MODEL TOOK. Where the leaf is not moving
  // water it does not place ci by the residual at all -- it assigns the
  // compensation point, which no carbon input reads. Applying the theorem to a
  // condition the model did not solve gives a number rather than an error, and
  // the number is a cancellation of two equal terms: measured at 5.9e-07 where
  // the row is zero, against 1e-13 everywhere the condition does place ci.
  double dci = 0.0;
  if (!l.ci_at_compensation_point_) {
    // The theorem's denominator: the residual's own slope in the variable it
    // places, taken with every input passive.
    tangent ci_seed = ci;
    seed_direction(ci_seed, 1.0);
    const carbon_inputs<tangent> passive = inputs(-1);
    const double R_ci = derivative_along(
        ci_residual(l, ci_seed, passive, gc, l.ca_, l.atm_kpa_));
    // The numerator: the same residual at the placed ci, with the input seeded.
    const double R_theta =
        derivative_along(ci_residual(l, tangent(ci), p, gc, l.ca_, l.atm_kpa_));
    dci = -R_theta / R_ci;
  }

  // Lift ci onto the tangent carrying what the theorem gives, and read profit
  // off the same kernels. The hydraulic cost is at a held stem potential and
  // reads no carbon input, so it contributes nothing.
  tangent ci_ad = ci;
  seed_direction(ci_ad, dci);
  const tangent J = l.electron_transport_kernel(p.ppfd, p.quantum_yield,
                                                p.curv_elec, p.transport_jmax);
  const tangent A =
      l.assim_colimited_kernel(ci_ad, p.vcmax, J, p.curv_colim, p.respiration);
  return derivative_along(A);
}

}  // namespace

const char* kName[] = {"vcmax_25", "jmax_25", "a", "curv_fact_elec_trans",
                       "curv_fact_colim", "PPFD", "R_d_25"};
inline constexpr int kInputs = 7;

// One state: the six rows both ways, and the worst relative disagreement. A row
// that is exactly zero both ways is reported and not scored -- there is nothing
// for a relative measure to say about it.
double one_state(double psi_soil, double ppfd, int layers, int& compared,
                 bool& interior) {
  const grad::Drivers d = drivers(psi_soil, ppfd, layers);
  Leaf l;
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  try {
    l.find_root_collar_psi();
  } catch (const std::exception& e) {
    std::printf("  %-8.3g %-8.4g %d   solve failed: %s\n", psi_soil, ppfd,
                layers, e.what());
    return 0.0;
  }
  interior = l.operating_point_kind() == Leaf::OperatingPointKind::Interior;
  // dpsistem_dp reaches only the marginal rows, so the profit rows compared here
  // are the same whatever it is.
  const Leaf::PhotoTraitRows ref = l.photo_trait_rows(0.0);
  const double hand[] = {ref.dprofit_dvcmax_25,   ref.dprofit_djmax_25,
                         ref.dprofit_da,          ref.dprofit_dcurv_elec,
                         ref.dprofit_dcurv_colim, ref.dprofit_dPPFD,
                         ref.dprofit_dR_d_25};

  // ⚠️ SCORED AGAINST THE ROW SET'S OWN SCALE, not each row's own magnitude.
  // At a wet pin the conductance is ~1e-08 and vcmax_25's row is ~1e-10 while
  // the rest of the same set are ~1e-05, so a per-row relative measure reports a
  // 2e-08 disagreement that is 2e-18 in absolute terms. What a consumer sums is
  // the set, so the set's largest entry is what a term in it is small against.
  double scale = 0.0;
  for (int i = 0; i < kInputs; ++i) {
    if (util_finite(hand[i]) && std::abs(hand[i]) > scale) {
      scale = std::abs(hand[i]);
    }
  }
  double worst = 0.0;
  int at = -1;
  for (int i = 0; i < kInputs; ++i) {
    if (!util_finite(hand[i])) {
      continue;
    }
    const double got = dprofit_generic(l, i);
    if (hand[i] == 0.0 && got == 0.0) {
      continue;
    }
    ++compared;
    const double rel = std::abs(got - hand[i]) / (scale > 0.0 ? scale : 1.0);
    if (rel > worst) {
      worst = rel;
      at = i;
    }
  }
  std::printf("  %-8.3g %-8.4g %d   %-10s %12.3e  %s\n", psi_soil, ppfd, layers,
              Leaf::operating_point_kind_name(l.operating_point_kind()), worst,
              at < 0 ? "(nothing to compare)" : kName[at]);
  if (worst > 1e-11) {
    std::printf("      gc %.6e  ci %.6e  compensation %d  psi_stem %.6e "
                " collar %.6e\n",
                l.stom_cond_CO2_, l.ci_, int(l.ci_at_compensation_point_),
                l.opt_psi_stem_, l.opt_root_psi_);
    for (int i = 0; i < kInputs; ++i) {
      std::printf("      %-22s hand %18.10e  tangent %18.10e\n", kName[i],
                  hand[i], dprofit_generic(l, i));
    }
  }
  return worst;
}

int main() {
  std::printf("The carbon-side profit rows: `photo_trait_rows` against one "
              "generic tangent pass.\n\n");
  std::printf("  %-8s %-8s %s   %-10s %12s  %s\n", "psi_soil", "PPFD", "L",
              "point", "worst rel.", "at");

  // Two claims, because the two branches have different floors and one bound
  // covering both would be the looser one wearing the tighter one's name. Where
  // the condition places ci the two forms are a reassociation of one expression.
  // At a wet pin the conductance is ~1e-08 and the theorem's quotient is a
  // difference of two near-equal terms, so the floor there is the cancellation's
  // and not the arithmetic's.
  double worst_interior = 0.0, worst_pinned = 0.0;
  int compared = 0, pinned = 0;
  for (double psi_soil : {0.2, 0.8, 1.5, 2.5, 3.5, 4.5}) {
    for (double ppfd : {30.0, 300.0, 1000.0, 2000.0}) {
      for (int layers : {1, 3}) {
        bool interior = true;
        const double w = one_state(psi_soil, ppfd, layers, compared, interior);
        double& into = interior ? worst_interior : worst_pinned;
        if (!interior) {
          ++pinned;
        }
        if (w > into) {
          into = w;
        }
      }
    }
  }
  std::printf("\n%d rows over %d states.\n", compared, 6 * 4 * 2);
  std::printf("  interior (%d states): worst %.3e\n", 6 * 4 * 2 - pinned,
              worst_interior);
  std::printf("  pinned   (%d states): worst %.3e\n", pinned, worst_pinned);
  return (worst_interior <= 1e-12 && worst_pinned <= 1e-8) ? 0 : 1;
}
