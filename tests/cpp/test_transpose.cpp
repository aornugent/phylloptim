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
// ⚠️ THE COLLAR IS HELD. Every quantity here is placed AT a collar rather than
// deciding one, so the collar carries no rows and there is no envelope reasoning
// to check yet. The live-collar arm, and the control that holds it passive on
// one side only, arrive with collar_at -- until then this checks the two
// residual lifts, the grafts underneath them, the supply path and the objective,
// which is what exists to be checked.
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

namespace pl = phylloptim;
using T = odelia::ode::tangent_scalar<double>;
using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;

namespace {

int failures = 0;
int compared = 0;
double worst = 0.0;

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
                    const std::vector<double>& v) {
  Inputs<T> in = inputs_of<T>(l);
  std::vector<T*> f = in.flat();
  for (std::size_t k = 0; k < f.size(); ++k) {
    odelia::ode::seed_direction(*f[k], u[k]);
  }
  const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), in.supply());
  const auto co = l.collar_coords_at<T>(l.opt_psi_stem_, l.ci_,
                                        T(l.opt_root_psi_), draw, in.pars);
  const auto o = l.outputs_at<pl::Leaf::CostCurve::TF24, T>(
      T(l.opt_root_psi_), draw, in.pars);

  double acc = v[0] * odelia::ode::derivative_along(co.sigma.value) +
               v[1] * odelia::ode::derivative_along(co.ci.value) +
               v[2] * odelia::ode::derivative_along(draw.flux.value) +
               v[3] * odelia::ode::derivative_along(o.profit);
  for (std::size_t i = 0; i < o.uptake.size(); ++i) {
    acc += v[4 + i] * odelia::ode::derivative_along(o.uptake[i]);
  }
  return acc;
}

// <J^T v, u>: one recording, one sweep.
double reverse_side(const pl::Leaf& l, const std::vector<double>& u,
                    const std::vector<double>& v) {
  Tape tape;
  Inputs<A> in = inputs_of<A>(l);
  std::vector<A*> f = in.flat();
  for (A* p : f) tape.registerInput(*p);
  tape.newRecording();

  auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), in.supply());
  auto co = l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, A(l.opt_root_psi_),
                                  draw, in.pars);
  auto o = l.outputs_at<pl::Leaf::CostCurve::TF24, A>(A(l.opt_root_psi_), draw,
                                                      in.pars);

  tape.registerOutput(co.sigma.value);
  xad::derivative(co.sigma.value) = v[0];
  tape.registerOutput(co.ci.value);
  xad::derivative(co.ci.value) = v[1];
  tape.registerOutput(draw.flux.value);
  xad::derivative(draw.flux.value) = v[2];
  tape.registerOutput(o.profit);
  xad::derivative(o.profit) = v[3];
  for (std::size_t i = 0; i < o.uptake.size(); ++i) {
    tape.registerOutput(o.uptake[i]);
    xad::derivative(o.uptake[i]) = v[4 + i];
  }
  tape.computeAdjoints();

  double acc = 0.0;
  for (std::size_t k = 0; k < f.size(); ++k) {
    acc += u[k] * xad::derivative(*f[k]);
  }
  return acc;
}

void check(int layers, double psi0, double ppfd, double leaf_temp,
           unsigned long seed) {
  pl::Leaf l = set_up(layers, psi0, ppfd, leaf_temp);
  l.find_root_collar_psi();
  if (l.operating_point_kind() == pl::Leaf::OperatingPointKind::SolverRefused ||
      l.operating_point_kind() ==
          pl::Leaf::OperatingPointKind::NonFiniteGradient) {
    return;
  }

  Stream rng{seed};
  const std::size_t nin = pl::n_pars + 3 * std::size_t(layers);
  const std::size_t nout = 4 + std::size_t(layers);
  std::vector<double> u(nin), v(nout);
  for (double& x : u) x = rng.next();
  for (double& x : v) x = rng.next();

  const double fwd = forward_side(l, u, v);
  const double rev = reverse_side(l, u, v);
  const double scale = std::max(std::abs(fwd), std::abs(rev));
  const double rel = (scale > 0.0) ? std::abs(fwd - rev) / scale : 0.0;
  ++compared;
  if (rel > worst) worst = rel;
  ok(rel < 1e-10,
     "transpose identity: layers=" + std::to_string(layers) + " psi_soil=" +
         std::to_string(psi0) + " ppfd=" + std::to_string(ppfd) + " T=" +
         std::to_string(leaf_temp) + " rel=" + std::to_string(rel));
}

}  // namespace

int main() {
  printf("the transpose identity at a held collar\n");
  unsigned long seed = 12345;
  for (int layers : {1, 3, 5}) {
    for (double psi0 : {0.5, 1.0, 2.0, 3.0}) {
      for (double ppfd : {100.0, 1500.0}) {
        for (double leaf_temp : {25.0, 40.0}) {
          check(layers, psi0, ppfd, leaf_temp, seed += 7919);
        }
      }
    }
  }
  printf("    %d operating points compared | worst relative %.3e\n", compared,
         worst);
  // The identity is only evidence about the points it ran on, so a fixture that
  // stops reaching them has narrowed the check without narrowing the claim.
  ok(compared >= 40, "the sweep reached at least 40 operating points");
  printf("%d checks, %d failures\n", compared + 1, failures);
  return failures == 0 ? 0 : 1;
}
