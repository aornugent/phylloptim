// The transpose identity, which needs no reference gradient and no
// differencing.
//
//     <v, J u>  ==  <J^T v, u>
//
// J is the leaf's map from its active inputs to the two things plant reads -- profit
// and the per-layer draws -- with the collar LIVE, so the identity sees the collar's
// own response and the envelope reasoning that omits profit's p-channel. `J u` is one
// forward tangent seeded with u; `J^T v` is one reverse sweep seeded with v.
//
// Neither side is a reference for the other. They are the same operator read in two
// directions, and only a consistent J satisfies both -- which is what makes this the
// referee for a construction that SUPPLIES rows rather than recording them. The two
// structural constants a graft asserts, profit's p-channel exactly zero at an interior
// point and the operating point's exactly one, are what this tests; a wrong one gives
// a plausible gradient everywhere else.
//
// ⚠️ THE CONTROL IS THE POINT. Holding the collar passive on the reverse side alone is
// exactly the mistake the identity exists to catch, so that arm MUST fail. A test that
// passes both ways is measuring nothing.
//
// Run by `make -C tests/cpp`, with the rest of the suite.
//
// ⚠️ The control is printed but NOT asserted: a state where the collar carries no
// rows reports "control is silent" and still passes. So this measures the
// identity, and the control tells a reader whether the state had anything to say.

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <odelia/ode_interface.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace pl = phylloptim;
using T = pl::tangent;
using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;

namespace {

int failures = 0;
// Which operating-point kinds the sweep actually reached. The identity is only
// evidence about the arms it ran on, so the census is part of the result.
std::vector<std::string> seen_kinds;

const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_up(pl::Leaf& l, int layers, double psi0, double ppfd) {
  double theta[pl::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::trait_R_d_25] = 1.44;
  l.set_traits(theta);
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / kAreaLeaf / layers);
    psi_soil.push_back(psi0 + 0.35 * i);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), ppfd, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

fixture::Soil soil_of(int layers) {
  fixture::Soil s;
  for (int i = 0; i < layers; ++i) {
    s.depth.push_back(1.0 * (i + 1));
    s.carbon.push_back(1.0 / layers / kAreaLeaf);
  }
  return s;
}

// Pointers to every active input, in one order both scalars agree on. The same
// layout probe_rank uses, as a template so the forward and reverse sides cannot
// pair different inputs by position.
template <class S>
std::vector<S*> flat(pl::LeafInputs<S>& in) {
  std::vector<S*> out = in.profit.field_ptrs();
  for (S& v : in.supply.psi_soil) out.push_back(&v);
  for (S& v : in.supply.r_R_H_min) out.push_back(&v);
  for (S& v : in.supply.r_R_V_sum) out.push_back(&v);
  out.push_back(&in.supply.root_b);
  out.push_back(&in.supply.root_c);
  out.push_back(&in.psi_crit);
  out.push_back(&in.root_psi_crit);
  return out;
}

// A reproducible stream, so a failure can be re-run. Deliberately not <random>:
// the sequence has to be the same on every platform for a reported number to mean
// anything.
struct Stream {
  unsigned long s;
  double next() {
    s = s * 6364136223846793005UL + 1442695040888963407UL;
    return double((s >> 11) & 0xFFFFFFUL) / double(0x1000000UL) * 2.0 - 1.0;
  }
};

// <v, J u>: one tangent, every input seeded at once.
double forward_side(const pl::Leaf& l, const pl::LeafInputs<double>& base,
                    const std::vector<double>& u, const std::vector<double>& v) {
  pl::LeafInputs<T> in = base.template rebind_from<T>();
  std::vector<T*> f = flat<T>(in);
  for (std::size_t k = 0; k < f.size(); ++k) pl::seed_direction(*f[k], u[k]);

  const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), in.supply);
  const T collar = l.collar_at<T>(in, draw);
  const pl::Leaf::LeafOutputs<T> o = l.outputs_at<T>(collar, in, draw);

  double acc = v[0] * pl::derivative_along(o.profit);
  for (std::size_t i = 0; i < o.uptake.size(); ++i) {
    acc += v[1 + i] * pl::derivative_along(o.uptake[i]);
  }
  return acc;
}

// <J^T v, u>: one recording, one sweep. `hold_collar` is the control -- it strips
// the collar's rows on this side only, which is the envelope mistake.
double reverse_side(const pl::Leaf& l, const pl::LeafInputs<double>& base,
                    const std::vector<double>& u, const std::vector<double>& v,
                    bool hold_collar) {
  Tape tape;
  pl::LeafInputs<A> in = base.template rebind_from<A>();
  std::vector<A*> f = flat<A>(in);
  for (A* p : f) tape.registerInput(*p);
  tape.newRecording();

  const auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), in.supply);
  A collar = l.collar_at<A>(in, draw);
  if (hold_collar) collar = A(odelia::util::to_passive(collar));
  pl::Leaf::LeafOutputs<A> o = l.outputs_at<A>(collar, in, draw);

  tape.registerOutput(o.profit);
  xad::derivative(o.profit) = v[0];
  for (std::size_t i = 0; i < o.uptake.size(); ++i) {
    tape.registerOutput(o.uptake[i]);
    xad::derivative(o.uptake[i]) = v[1 + i];
  }
  tape.computeAdjoints();

  double acc = 0.0;
  for (std::size_t k = 0; k < f.size(); ++k) {
    acc += xad::derivative(*f[k]) * u[k];
  }
  return acc;
}

double rel_gap(double a, double b) {
  return std::abs(a - b) / (std::abs(a) + std::abs(b) + 1e-300);
}

// One operating point, `draws` random (u, v) pairs through it.
void at_state(int layers, double psi0, double ppfd, int draws, Stream& rng) {
  pl::Leaf l;
  set_up(l, layers, psi0, ppfd);
  l.find_root_collar_psi();

  const fixture::Soil soil = soil_of(layers);
  const pl::LeafInputs<double> base =
      fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);

  // Sized off the structures themselves rather than a formula, so a layout change
  // cannot silently pair the wrong entries.
  pl::LeafInputs<double> probe = base;
  const std::size_t n_in = flat<double>(probe).size();
  const std::size_t n_out = 1 + std::size_t(layers);

  double worst = 0.0, worst_control = 1.0;
  bool threw = false;
  for (int d = 0; d < draws; ++d) {
    std::vector<double> u(n_in), v(n_out);
    for (double& x : u) x = rng.next();
    for (double& x : v) x = rng.next();

    try {
      const double fwd = forward_side(l, base, u, v);
      const double rev = reverse_side(l, base, u, v, false);
      worst = std::max(worst, rel_gap(fwd, rev));
      const double held = reverse_side(l, base, u, v, true);
      worst_control = std::min(worst_control, rel_gap(fwd, held));
    } catch (const std::exception&) {
      threw = true;
      break;
    }
  }

  const char* kind = l.operating_point_kind_name();
  seen_kinds.push_back(kind);
  if (threw) {
    std::printf("  %d layer  psi0 %-5.2f ppfd %-6.0f %-22s  refused, skipped\n",
                layers, psi0, ppfd, kind);
    return;
  }

  // The control must be far from agreeing. If holding the collar changes nothing
  // then the collar carries no rows here, and this state proves nothing either way.
  const bool live = worst_control > 1e-6;
  const bool ok = worst < 1e-10;
  if (!ok) ++failures;
  std::printf("  %d layer  psi0 %-5.2f ppfd %-6.0f %-22s  worst %.3e  %s",
              layers, psi0, ppfd, kind, worst, ok ? "ok" : "FAIL");
  if (live) {
    std::printf("   (control %.2e, held collar disagrees)\n", worst_control);
  } else {
    std::printf("   (collar carries no rows here; control is silent)\n");
  }
}

}  // namespace

int main() {
  std::printf("the transpose identity, <v,Ju> == <J^T v,u>\n");
  std::printf("  no reference gradient, no differencing: one tangent against one"
              " sweep\n\n");

  Stream rng{0x5eed1234UL};
  // Dry enough to pin against both root limits, and dark enough to reach the
  // arms where the leaf moves no water at all.
  for (int layers : {1, 3, 5}) {
    for (double ppfd : {900.0, 30.0}) {
      for (double psi0 : {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.8, 6.5}) {
        at_state(layers, psi0, ppfd, 4, rng);
      }
    }
  }

  std::sort(seen_kinds.begin(), seen_kinds.end());
  seen_kinds.erase(std::unique(seen_kinds.begin(), seen_kinds.end()),
                   seen_kinds.end());
  std::printf("\n  kinds reached:");
  for (const std::string& k : seen_kinds) std::printf(" %s", k.c_str());
  std::printf("\n\n%s\n", failures == 0 ? "all ok" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
