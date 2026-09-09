// The transpose identity, which needs no reference gradient and no differencing.
//
//     <v, J u>  ==  <J^T v, u>
//
// J is the map from the leaf's active inputs -- the parameter pack and the soil
// state -- to the two things plant reads, profit and the per-layer draws, plus
// the two coordinates they are assembled from. `J u` is one forward tangent
// seeded with u; `J^T v` is one reverse sweep seeded with v.
//
// Neither side is a reference for the other. They are the same operator read in
// two directions, and only a consistent J satisfies both -- which is what makes
// this the referee for a construction that SUPPLIES rows rather than recording
// them. A difference of the forward model cannot see a wrong supplied row at
// all, and at the shipped spline resolution a difference cannot referee the
// collar responses either: sigma comes from an inverse spline whose round-trip
// dominates dsigma/dcollar - 1. This can.
//
// It cannot referee the VALUES, though. Both sides read one surface, so a kernel
// evaluated at the wrong temperature or the wrong light satisfies this identity
// exactly -- test_leaf's member-against-pack comparison is what covers that.
//
// ⚠️ THE RESIDUAL IS SCALED BY THE TERMS, NOT BY THE PAIRING. Both sides are a
// sum over inputs of u_k times a row, and the sum can cancel: over this sweep the
// pairing runs from 1.0 to 1.6e+03 times smaller than the sum of its terms'
// magnitudes. Dividing by the pairing therefore measures rounding in a
// cancellation wherever it is small, and it measures it against a threshold set
// for the cases where it is not -- one point read 4.0e-09 that way with every
// individual row agreeing to 3.3e-12. Scaling by the terms gives a wrong row of
// relative size e a residual of order e whatever the pairing does, which is the
// sensitivity this is here for.
//
// ⚠️ THE CONTROL IS THE POINT. The collar is LIVE -- placed by collar_at, so it
// carries the rows of whatever condition pins it -- and the second arm holds it
// passive on the REVERSE side only, which is exactly the envelope mistake. That
// arm MUST fail: a test that passes both ways is measuring nothing, so the
// control is asserted rather than reported.
//
// Where the collar genuinely does not move -- a hydraulic shutdown holds the stem
// at psi_crit and nothing defines the collar at all -- holding it must change
// NOTHING, and that is asserted too. Either expectation broken is a failure.
//
// Run by `make -C tests/cpp`, with the rest of the suite.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <odelia/ode_interface.hpp>
#include <odelia/tangent.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <stdexcept>

namespace pl = phylloptim;
using T = odelia::ode::tangent_scalar<double>;
using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;

namespace {

int failures = 0;
int compared = 0;
double worst = 0.0;
double worst_cancel = 0.0;

void ok(bool cond, const std::string& what) {
  if (!cond) {
    ++failures;
    printf("FAIL %s\n", what.c_str());
  }
}

// A reproducible stream, so a failure can be re-run. Deliberately not <random>:
// the sequence has to be the same on every platform for a reported number to
// mean anything.
struct Stream {
  unsigned long s;
  double next() {
    s = s * 6364136223846793005UL + 1442695040888963407UL;
    return double((s >> 11) & 0xFFFFFFUL) / double(0x1000000UL) * 2.0 - 1.0;
  }
};

const double kTheta = 0.5, kKs = 2.0, kH = 5.0, kAreaLeaf = 1.0;

pl::Leaf set_up(int layers, double psi0, double ppfd, double leaf_temp) {
  pl::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / kAreaLeaf / layers);
    psi_soil.push_back(psi0 + 0.35 * i);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), ppfd, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, leaf_temp, 21.0, 101.3);
  return l;
}

// Every active input, in one order both scalars agree on. Written as a template
// so the forward and reverse sides cannot pair different inputs by position.
template <class S>
struct Inputs {
  pl::leaf_pars<S> pars;
  std::vector<S> psi_soil, r_R_H_min, r_R_V_sum;

  std::vector<S*> flat() {
    std::vector<S*> out;
    for (S& v : pars) out.push_back(&v);
    for (S& v : psi_soil) out.push_back(&v);
    for (S& v : r_R_H_min) out.push_back(&v);
    for (S& v : r_R_V_sum) out.push_back(&v);
    return out;
  }
  pl::SupplyAt<S> supply() const {
    return {psi_soil, r_R_H_min, r_R_V_sum, pars[pl::par_root_P50],
            pars[pl::par_root_c]};
  }
};

template <class S>
Inputs<S> inputs_of(const pl::Leaf& l) {
  Inputs<S> in;
  const pl::leaf_pars<double> p = l.passive_pars();
  for (std::size_t i = 0; i < p.size(); ++i) in.pars[i] = S(p[i]);
  const int n = l.supply_n_layers();
  for (int i = 0; i < n; ++i) {
    in.psi_soil.push_back(S(l.roots_.psi_soil_[std::size_t(i)]));
    in.r_R_H_min.push_back(S(l.roots_.network_.r_R_H_min[std::size_t(i)]));
    in.r_R_V_sum.push_back(S(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
  }
  return in;
}

// <v, J u>: one tangent, every input seeded at once.
double forward_side(const pl::Leaf& l, const std::vector<double>& u,
                    const std::vector<double>& v, double curvature,
                    bool has_coords) {
  Inputs<T> in = inputs_of<T>(l);
  std::vector<T*> f = in.flat();
  for (std::size_t k = 0; k < f.size(); ++k) {
    odelia::ode::seed_direction(*f[k], u[k]);
  }
  const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), in.supply());
  const T collar =
      l.collar_at<pl::Leaf::CostCurve::TF24, T>(draw, in.pars, curvature);
  const auto o = l.outputs_at<pl::Leaf::CostCurve::TF24, T>(collar, draw, in.pars);

  double acc = v[2] * odelia::ode::derivative_along(draw.flux.value) +
               v[3] * odelia::ode::derivative_along(o.profit);
  if (has_coords) {
    const auto co =
        l.collar_coords_at<T>(l.opt_psi_stem_, l.ci_, collar, draw, in.pars,
                              l.photo_capacity_at<T>(in.pars), true);
    acc += v[0] * odelia::ode::derivative_along(co.sigma.value) +
           v[1] * odelia::ode::derivative_along(co.ci.value);
  }
  for (std::size_t i = 0; i < o.uptake.size(); ++i) {
    acc += v[4 + i] * odelia::ode::derivative_along(o.uptake[i]);
  }
  return acc;
}

// <J^T v, u> and the size of the terms it was summed from.
struct Pairing {
  double value = 0.0;
  double terms = 0.0;
};

// One recording, one sweep.
Pairing reverse_side(const pl::Leaf& l, const std::vector<double>& u,
                     const std::vector<double>& v, double curvature,
                     bool hold_collar, bool has_coords) {
  Tape tape;
  Inputs<A> in = inputs_of<A>(l);
  std::vector<A*> f = in.flat();
  for (A* p : f) tape.registerInput(*p);
  tape.newRecording();

  auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), in.supply());
  A collar = l.collar_at<pl::Leaf::CostCurve::TF24, A>(draw, in.pars, curvature);
  // The control: strip the collar's rows on THIS side only.
  if (hold_collar) collar = A(odelia::util::to_passive(collar));
  auto o = l.outputs_at<pl::Leaf::CostCurve::TF24, A>(collar, draw, in.pars);

  if (has_coords) {
    auto co = l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, collar, draw,
                                    in.pars, l.photo_capacity_at<A>(in.pars),
                                    true);
    tape.registerOutput(co.sigma.value);
    xad::derivative(co.sigma.value) = v[0];
    tape.registerOutput(co.ci.value);
    xad::derivative(co.ci.value) = v[1];
  }
  tape.registerOutput(draw.flux.value);
  xad::derivative(draw.flux.value) = v[2];
  tape.registerOutput(o.profit);
  xad::derivative(o.profit) = v[3];
  for (std::size_t i = 0; i < o.uptake.size(); ++i) {
    tape.registerOutput(o.uptake[i]);
    xad::derivative(o.uptake[i]) = v[4 + i];
  }
  tape.computeAdjoints();

  Pairing out;
  for (std::size_t k = 0; k < f.size(); ++k) {
    const double term = u[k] * xad::derivative(*f[k]);
    out.value += term;
    out.terms += std::abs(term);
  }
  return out;
}

// Whether the collar carries rows on this arm, which is what decides what the
// control MUST do. Holding it passive has to change the answer wherever the
// collar is placed by a rule that reads the inputs, and has to change nothing
// where the collar does not move at all.
//
// No `default:`, so -Werror=switch refuses a new kind until someone says which
// of the three it is.
enum class Arm { CarriesRows, NoRows, Refuses };

Arm arm_of(pl::Leaf::OperatingPointKind kind) {
  switch (kind) {
    // Placed by a rule that reads the inputs: the interior condition, or one of
    // the two bounds.
    case pl::Leaf::OperatingPointKind::Interior:
    case pl::Leaf::OperatingPointKind::BoundarySoil:
    case pl::Leaf::OperatingPointKind::BoundaryCrit:
    case pl::Leaf::OperatingPointKind::BoundaryRootCrit:
    case pl::Leaf::OperatingPointKind::ShadeDeath:
      return Arm::CarriesRows;
    // The stem holds at psi_crit and nothing defines the collar at all.
    case pl::Leaf::OperatingPointKind::HydraulicShutdown:
      return Arm::NoRows;
    case pl::Leaf::OperatingPointKind::Unsolved:
    case pl::Leaf::OperatingPointKind::Determined:
    case pl::Leaf::OperatingPointKind::Prescribed:
    case pl::Leaf::OperatingPointKind::SolverRefused:
    case pl::Leaf::OperatingPointKind::NonFiniteGradient:
      return Arm::Refuses;
  }
  return Arm::Refuses;
}

// A control below this is "silent": holding the collar passive changed nothing.
// The identity itself passes below 1e-10, so a live control at this floor still
// sits three orders above what the check can resolve. Scaled by the terms too,
// for the reason the identity is.
const double kControlFloor = 1e-7;
int seen_carry = 0, seen_norows = 0;

void check(int layers, double psi0, double ppfd, double leaf_temp,
           unsigned long seed) {
  pl::Leaf l = set_up(layers, psi0, ppfd, leaf_temp);
  l.find_root_collar_psi();
  Stream rng{seed};
  const std::size_t nin = pl::n_pars + 3 * std::size_t(layers);
  const std::size_t nout = 4 + std::size_t(layers);
  std::vector<double> u(nin), v(nout);
  for (double& x : u) x = rng.next();
  for (double& x : v) x = rng.next();

  const Arm arm = arm_of(l.operating_point_kind());
  if (arm == Arm::Refuses) return;
  const double curvature =
      l.operating_point_kind() == pl::Leaf::OperatingPointKind::Interior
          ? l.marginal_collar_slope<pl::Leaf::CostCurve::TF24>()
          : std::numeric_limits<double>::quiet_NaN();
  if (l.operating_point_kind() == pl::Leaf::OperatingPointKind::Interior &&
      !std::isfinite(curvature)) {
    return;  // no curvature to divide by; collar_at would refuse too
  }

  const bool has_coords =
      l.operating_point_kind() != pl::Leaf::OperatingPointKind::HydraulicShutdown;
  const double fwd = forward_side(l, u, v, curvature, has_coords);
  const Pairing rev = reverse_side(l, u, v, curvature, false, has_coords);

  // The control, on the reverse side only.
  const Pairing held = reverse_side(l, u, v, curvature, true, has_coords);
  const double moved =
      std::abs(rev.value - held.value) / std::max(rev.terms, 1e-300);
  if (arm == Arm::CarriesRows) {
    ++seen_carry;
    ok(moved > kControlFloor,
       "CONTROL: holding the collar must change the answer where it is placed "
       "by a rule -- kind=" +
           std::string(l.operating_point_kind_name(l.operating_point_kind())) +
           " moved=" + std::to_string(moved));
  } else {
    ++seen_norows;
    ok(moved <= kControlFloor,
       "CONTROL: holding the collar must change nothing where it does not move "
       "-- kind=" +
           std::string(l.operating_point_kind_name(l.operating_point_kind())) +
           " moved=" + std::to_string(moved));
  }
  const double rel = (rev.terms > 0.0)
                         ? std::abs(fwd - rev.value) / rev.terms
                         : 0.0;
  const double cancel = rev.terms / std::max(std::abs(rev.value), 1e-300);
  ++compared;
  if (rel > worst) worst = rel;
  if (cancel > worst_cancel) worst_cancel = cancel;
  ok(rel < 1e-10,
     "transpose identity: layers=" + std::to_string(layers) + " psi_soil=" +
         std::to_string(psi0) + " ppfd=" + std::to_string(ppfd) + " T=" +
         std::to_string(leaf_temp) + " rel=" + std::to_string(rel) +
         " cancellation=" + std::to_string(cancel));
}

}  // namespace

int main() {
  printf("the transpose identity, with the collar live\n");
  unsigned long seed = 12345;
  for (int layers : {1, 3, 5}) {
    // Dry enough to reach a hydraulic shutdown, which is the arm the control's
    // other half needs -- the sweep asserts below that it got there.
    for (double psi0 : {0.5, 1.0, 2.0, 3.0, 5.0, 7.0}) {
      for (double ppfd : {100.0, 1500.0}) {
        for (double leaf_temp : {25.0, 40.0}) {
          check(layers, psi0, ppfd, leaf_temp, seed += 7919);
        }
      }
    }
  }
  printf("    %d operating points compared | worst relative %.3e\n", compared,
         worst);
  // Reported because it is what the scaling above is for: a sweep that stopped
  // reaching a cancelling pairing would have stopped exercising the reason.
  printf("    worst cancellation in the pairing: %.3gx\n", worst_cancel);
  // The identity is only evidence about the points it ran on, so a fixture that
  // stops reaching them has narrowed the check without narrowing the claim.
  printf("    arms reached: %d placed by a rule, %d with no collar to move\n",
         seen_carry, seen_norows);
  // The identity is only evidence about the arms it ran on, and the control is
  // only evidence where both expectations were exercised.
  ok(compared >= 20, "the sweep reached at least 20 operating points");
  ok(seen_carry > 0, "the sweep reached an arm whose collar carries rows");
  ok(seen_norows > 0, "the sweep reached an arm whose collar does not move");
  printf("%d checks, %d failures\n", compared + seen_carry + seen_norows + 3, failures);
  return failures == 0 ? 0 : 1;
}
