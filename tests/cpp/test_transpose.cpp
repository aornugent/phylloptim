// The transpose identity, which needs no reference gradient and no differencing.
//
//     <v, J u>  ==  <J^T v, u>
//
// J is the map from the leaf's active inputs -- the parameter pack and the soil
// state -- to the two things plant reads, profit and the per-layer draws, plus
// the two coordinates they are assembled from. `J u` is one forward tangent
// along u; `J^T v` is one reverse sweep from v.
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
#include <map>
#include <cstdio>
#include <string>
#include <vector>
#include <stdexcept>

namespace pl = phylloptim;
using T = odelia::ode::tangent_scalar<double>;
using A = odelia::ode::adjoint_scalar<double>;
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

pl::Leaf set_up(int layers, double psi0, double ppfd, double leaf_temp,
                double g_min = 0.0) {
  pl::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  // The residual conductance is a trait (#182): seated through set_traits so
  // the g_min rows the pack carries are exercised at a non-zero value as well
  // as at the default.
  l.set_traits(l.vcmax_25, l.stem_c, l.stem_P50, l.roots_.root_c,
               l.roots_.root_P50, l.TF24_beta2, l.jmax_25, l.phi_J,
               l.curv_fact_elec_trans, l.curv_fact_colim, l.TF24_cost_scale,
               l.R_d_25, l.JS22_gamma, l.CMax_a, l.CMax_b, g_min);
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

// <v, J u>: one tangent, every input given its direction at once.
double forward_side(pl::Leaf l, const std::vector<double>& u,
                    const std::vector<double>& v, double curvature,
                    bool has_coords, bool handoff) {
  Inputs<T> in = inputs_of<T>(l);
  std::vector<T*> f = in.flat();
  for (std::size_t k = 0; k < f.size(); ++k) {
    odelia::ode::set_direction(*f[k], u[k]);
  }
  const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), in.supply());
  pl::Leaf::HeldCoords<T> held;
  pl::Leaf::HeldCoords<T>* h = handoff ? &held : nullptr;
  const T collar =
      l.collar_at<pl::Leaf::CostCurve::TF24, T>(draw, in.pars, curvature, h);
  const auto o =
      l.outputs_at<pl::Leaf::CostCurve::TF24, T>(collar, draw, in.pars, h);

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
  acc += v[4 + o.uptake.size()] * odelia::ode::derivative_along(o.residual_flux);
  return acc;
}

// <J^T v, u> and the size of the terms it was summed from.
struct Pairing {
  double value = 0.0;
  double terms = 0.0;
};

// One recording, one sweep.
Pairing reverse_side(pl::Leaf l, const std::vector<double>& u,
                     const std::vector<double>& v, double curvature,
                     bool hold_collar, bool has_coords, bool handoff) {
  Tape tape;
  Inputs<A> in = inputs_of<A>(l);
  std::vector<A*> f = in.flat();
  for (A* p : f) tape.registerInput(*p);
  tape.newRecording();

  auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), in.supply());
  pl::Leaf::HeldCoords<A> held;
  pl::Leaf::HeldCoords<A>* h = handoff ? &held : nullptr;
  A collar =
      l.collar_at<pl::Leaf::CostCurve::TF24, A>(draw, in.pars, curvature, h);
  // The control: strip the collar's rows on THIS side only.
  if (hold_collar) collar = A(odelia::util::to_passive(collar));
  auto o = l.outputs_at<pl::Leaf::CostCurve::TF24, A>(collar, draw, in.pars, h);

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
  tape.registerOutput(o.residual_flux);
  xad::derivative(o.residual_flux) = v[4 + o.uptake.size()];
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
// collar is placed by a rule that reads the inputs. Every solved kind is one
// now: a shut leaf sits on its wet bound (#182), which reads the soil state and
// the traits through the residual flux, so there is no arm left whose collar
// does not move; the sweep still asserts it reaches a shutdown.
//
// No `default:`, so -Werror=switch refuses a new kind until someone says which
// of the two it is.
enum class Arm { CarriesRows, Refuses };

Arm arm_of(pl::Leaf::OperatingPointKind kind) {
  switch (kind) {
    // Placed by a rule that reads the inputs: the interior condition, or one of
    // the bounds, the shut leaf's wet bound among them.
    case pl::Leaf::OperatingPointKind::Interior:
    case pl::Leaf::OperatingPointKind::BoundarySoil:
    case pl::Leaf::OperatingPointKind::BoundaryCrit:
    case pl::Leaf::OperatingPointKind::BoundaryRootCrit:
    case pl::Leaf::OperatingPointKind::HydraulicShutdown:
      return Arm::CarriesRows;
    case pl::Leaf::OperatingPointKind::Unsolved:
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
int seen_carry = 0, seen_shut = 0, seen_residual = 0;
// The handoff against solving the coordinates again in outputs_at: the same rows
// summed in another order on the reverse side, so equal to rounding there and
// bit for bit on the forward side, where the arithmetic is the same.
double worst_handoff = 0.0;
int handoff_forward_differs = 0;
// Compared points per kind, by name. The identity is evidence only about the
// kinds it ran on, and collar_at places each by a different condition, so every
// kind that carries rows is asserted reached, not just the arm.
std::map<std::string, int> seen_kind;

void check_solved(pl::Leaf l, const std::string& at, unsigned long seed) {
  const std::size_t layers = std::size_t(l.supply_n_layers());
  Stream rng{seed};
  const std::size_t nin = pl::n_pars + 3 * layers;
  // The residual flux is the last output, and its weight is drawn last of all,
  // for the reason the comment below gives about draw order.
  const std::size_t nout = 5 + layers;
  std::vector<double> u(nin), v(nout);
  // ⚠️ DRAWN IN A FIXED ORDER, NOT SLOT ORDER: the parameters before the
  // thermal-damage layer, then PPFD, the soil state and v, and the layer's
  // slots last. Drawn in slot order, every parameter added to the enumeration
  // reweights every term, so the check moves with the parameter count rather
  // than with the leaf. That is not academic: the kmax row disagrees between the
  // two directions by up to 1e-4 relative, small enough to weigh ~1e-11 in the
  // pairing as drawn here, and the layer's sixteen extra slots once reshuffled
  // it to 7e-11, over the limit on Linux. The layer is off in this sweep, so its
  // rows are exact zeros and their draws change nothing.
  for (int k = 0; k < pl::par_A; ++k) u[k] = rng.next();
  u[pl::par_PPFD] = rng.next();
  for (std::size_t k = pl::n_pars; k < nin; ++k) u[k] = rng.next();
  for (std::size_t k = 0; k + 1 < nout; ++k) v[k] = rng.next();
  for (int k = pl::par_A; k < pl::par_PPFD; ++k) u[k] = rng.next();
  v[nout - 1] = rng.next();

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
  // Through the coordinates collar_at hands to outputs_at, which is how plant
  // calls the surface.
  const double fwd = forward_side(l, u, v, curvature, has_coords, true);
  const Pairing rev = reverse_side(l, u, v, curvature, false, has_coords, true);

  const double fwd_solved = forward_side(l, u, v, curvature, has_coords, false);
  const Pairing rev_solved =
      reverse_side(l, u, v, curvature, false, has_coords, false);
  if (fwd != fwd_solved) ++handoff_forward_differs;
  const double handoff_rel =
      std::abs(rev.value - rev_solved.value) / std::max(rev.terms, 1e-300);
  if (handoff_rel > worst_handoff) worst_handoff = handoff_rel;
  // Two orders below the identity's own threshold, which is what rounding in a
  // reordered sum reaches here (1.1e-13 on this sweep).
  char moved_by[32];
  std::snprintf(moved_by, sizeof moved_by, "%.3e", handoff_rel);
  ok(handoff_rel < 1e-12,
     "handing the coordinates over moves the reverse pairing past rounding: " +
         at + " rel=" + moved_by);

  // The control, on the reverse side only.
  const Pairing held = reverse_side(l, u, v, curvature, true, has_coords, true);
  const double moved =
      std::abs(rev.value - held.value) / std::max(rev.terms, 1e-300);
  ++seen_carry;
  if (l.operating_point_kind() == pl::Leaf::OperatingPointKind::HydraulicShutdown) {
    ++seen_shut;
  }
  if (l.g_min > 0.0) ++seen_residual;
  ok(moved > kControlFloor,
     "CONTROL: holding the collar must change the answer where it is placed "
     "by a rule -- kind=" +
         std::string(l.operating_point_kind_name(l.operating_point_kind())) +
         " g_min=" + std::to_string(l.g_min) +
         " moved=" + std::to_string(moved));
  const double rel = (rev.terms > 0.0)
                         ? std::abs(fwd - rev.value) / rev.terms
                         : 0.0;
  const double cancel = rev.terms / std::max(std::abs(rev.value), 1e-300);
  ++compared;
  if (rel > worst) worst = rel;
  if (cancel > worst_cancel) worst_cancel = cancel;
  ++seen_kind[l.operating_point_kind_name(l.operating_point_kind())];
  ok(rel < 1e-10, "transpose identity: " + at + " rel=" + std::to_string(rel) +
                      " cancellation=" + std::to_string(cancel));
}

void check(int layers, double psi0, double ppfd, double leaf_temp,
           unsigned long seed, double g_min) {
  pl::Leaf l = set_up(layers, psi0, ppfd, leaf_temp, g_min);
  l.find_root_collar_psi();
  check_solved(l,
               "layers=" + std::to_string(layers) + " psi_soil=" +
                   std::to_string(psi0) + " ppfd=" + std::to_string(ppfd) +
                   " T=" + std::to_string(leaf_temp) +
                   " g_min=" + std::to_string(g_min),
               seed);
}

// The kinds the sweep above does not reach, each on test_leaf's drivers and
// asserted to land on its kind before it is compared: a fixture that drifted to
// another kind would narrow the check without failing it. A root_P50 of zero
// keeps the default traits.
void check_kind(const std::vector<double>& psi_soil, double ppfd, double vpd,
                double root_P50, pl::Leaf::OperatingPointKind want,
                unsigned long seed) {
  pl::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> depth, root;
  for (std::size_t i = 0; i < psi_soil.size(); ++i) {
    depth.push_back(1.0 * double(i + 1));
    root.push_back(1.0 / 0.05 / double(psi_soil.size()));
  }
  if (root_P50 > 0.0) {
    l.set_traits(96.0, 2.680147, 3.4287, 2.680147, root_P50, 1.5, 157.44, 0.30,
                 0.7, 0.99, 7.5, 1.44, 1.0, 0.6, 0.0, /*g_min=*/0.0);
  }
  l.set_physiology(fixture::root_network(root, depth), ppfd, psi_soil, depth,
                   0.000157 / 5.0, vpd, 40.0, 25.0, 21.0, 101.3);
  l.find_root_collar_psi();
  const std::string at = std::string(pl::Leaf::operating_point_kind_name(want)) +
                         " fixture psi_soil=" + std::to_string(psi_soil[0]) +
                         " ppfd=" + std::to_string(ppfd) +
                         " vpd=" + std::to_string(vpd) +
                         " root_P50=" + std::to_string(root_P50);
  ok(l.operating_point_kind() == want,
     "the " + at + " is " +
         std::string(pl::Leaf::operating_point_kind_name(want)) + ", got " +
         l.operating_point_kind_name(l.operating_point_kind()));
  check_solved(l, at, seed);
}

}  // namespace

int main() {
  printf("the transpose identity, with the collar live\n");
  unsigned long seed = 12345;
  // At g_min = 0 and at a measured cuticular conductance (0.5 mmol), so the
  // g_min row and the residual's rows in kmax and the stem pair are covered.
  for (double g_min : {0.0, 5e-4}) {
    for (int layers : {1, 3, 5}) {
      // Dry enough to reach a hydraulic shutdown, which the sweep asserts
      // below that it got to.
      for (double psi0 : {0.5, 1.0, 2.0, 3.0, 5.0, 7.0}) {
        for (double ppfd : {100.0, 1500.0}) {
          for (double leaf_temp : {25.0, 40.0}) {
            check(layers, psi0, ppfd, leaf_temp, seed += 7919, g_min);
          }
        }
      }
    }
  }
  using K = pl::Leaf::OperatingPointKind;
  // The two dry bounds collar_at places by different conditions -- the
  // continuity root and the root's own critical potential.
  check_kind({4.0, 4.25, 4.5}, 900.0, 0.5, 0.0, K::BoundaryCrit, seed += 7919);
  check_kind({2.0}, 900.0, 2.0, 1.2, K::BoundaryRootCrit, seed += 7919);
  printf("    %d operating points compared | worst relative %.3e\n", compared,
         worst);
  // Reported because it is what the scaling above is for: a sweep that stopped
  // reaching a cancelling pairing would have stopped exercising the reason.
  printf("    worst cancellation in the pairing: %.3gx\n", worst_cancel);
  printf("    coordinates handed over against solved again: reverse worst %.3e, "
         "forward differs at %d points\n",
         worst_handoff, handoff_forward_differs);
  ok(handoff_forward_differs == 0,
     "handing the coordinates over moves no forward tangent, bit for bit");
  // The identity is only evidence about the points it ran on, so a fixture that
  // stops reaching them has narrowed the check without narrowing the claim.
  printf("    arms reached: %d placed by a rule, %d of them shut, %d with a "
         "residual conductance\n",
         seen_carry, seen_shut, seen_residual);
  printf("    kinds reached:");
  for (const auto& [name, count] : seen_kind) printf(" %s %d,", name.c_str(), count);
  printf("\n");
  // The identity is only evidence about the arms and kinds it ran on.
  ok(compared >= 20, "the sweep reached at least 20 operating points");
  ok(seen_shut > 0, "the sweep reached a shut leaf");
  ok(seen_residual > 0, "the sweep reached a leaf with a residual conductance");
  for (const char* kind : {"interior", "boundary-soil", "boundary-crit",
                           "boundary-root-crit", "hydraulic-shutdown"}) {
    ok(seen_kind[kind] > 0, std::string("the identity was checked at a ") + kind +
                                " point");
  }
  printf("%d checks, %d failures\n", compared + seen_carry + 8, failures);
  return failures == 0 ? 0 : 1;
}
