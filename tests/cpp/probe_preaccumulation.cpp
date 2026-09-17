// PROBE: what preaccumulating the leaf region costs in reproducibility.
//
// plant records the leaf per cohort per stage per step and sweeps the recording
// once per census metric. Taking the region off that tape -- recording it,
// sweeping it locally, and keeping the rows -- is a large saving and changes one
// thing this package's consumer asserts: the ORDER the rows are summed in.
//
// Measured here, against recording the region plainly on the same tape:
//
//   one region                        rows bit-identical, 0.000e+00
//   two sweeps of one recording       stable, 0.000e+00
//   TWO regions sharing inputs        9.8e-10 to 1.4e-08 relative
//
// The last is not an error term. A missing or double-counted row is O(1)
// relative; this is floating point, and the values agree to eight figures. When
// two regions read the same input, the plain route accumulates their
// contributions through the tape and this route adds two partial rows computed
// earlier -- the same number, summed in a different order.
//
// ⚠️ SO IT CANNOT BE FIXED, ONLY CHOSEN. plant's gradient ladder asserts
// BIT-identity in several places -- two consecutive sweeps, a permuted sweep
// order, a sweep split at an interior step -- and those fail by construction
// against a 1e-10 tolerance, not because a row is wrong.
//
//   make CXX=g++ probe_preaccumulation && ./probe_preaccumulation
#include <phylloptim.hpp>
#include "root_network.hpp"
#include <odelia/implicit_node.hpp>
#include <odelia/ode_interface.hpp>
#include <cmath>
#include <cstdio>
#include <span>
#include <vector>

namespace pl = phylloptim;
using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;
constexpr auto K = pl::Leaf::CostCurve::TF24;

namespace {
pl::Leaf set_up(int layers, double psi0) {
  pl::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  const double th = 0.5, ks = 2.0, h = 5.0, area = 1.0;
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / area / layers);
    psi_soil.push_back(psi0 + 0.35 * i);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   ks * th / h, 2.0, 40.0, 25.0, 21.0, 101.3);
  return l;
}

struct Result { double profit; std::vector<double> rows; bool ok; double resweep = 0.0; };

// `which` selects the output to seed: 0 = profit, 1+ = that uptake layer.
Result run(int layers, double psi0, bool preacc, std::size_t which) {
  pl::Leaf l = set_up(layers, psi0);
  l.find_root_collar_psi();
  if (l.operating_point_kind() != pl::Leaf::OperatingPointKind::Interior)
    return {0, {}, false};
  double curv = 0.0;
  try { curv = l.marginal_collar_slope<K>(); } catch (...) { return {0, {}, false}; }
  if (!(curv < 0.0)) return {0, {}, false};

  Tape tape;
  pl::leaf_pars<A> in;
  const pl::leaf_pars<double> p = l.passive_pars();
  for (std::size_t i = 0; i < p.size(); ++i) in[i] = A(p[i]);
  const int n = l.supply_n_layers();
  std::vector<A> psi_soil, rh, rv;
  for (int i = 0; i < n; ++i) {
    psi_soil.push_back(A(l.roots_.psi_soil_[std::size_t(i)]));
    rh.push_back(A(l.roots_.network_.r_R_H_min[std::size_t(i)]));
    rv.push_back(A(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
  }
  for (A& v : in) tape.registerInput(v);
  for (A& v : psi_soil) tape.registerInput(v);
  for (A& v : rh) tape.registerInput(v);
  for (A& v : rv) tape.registerInput(v);
  tape.newRecording();

  const pl::SupplyAt<A> supply{psi_soil, rh, rv, in[pl::par_root_P50],
                               in[pl::par_root_c]};
  pl::Leaf::LeafOutputs<A> got;
  static std::vector<A*> outs;
  std::vector<double> scratch;
  auto body = [&]() -> std::span<A* const> {
    const auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), supply);
    const A collar = l.collar_at<K, A>(draw, in, curv);
    got = l.outputs_at<K, A>(collar, draw, in);
    outs.clear();
    outs.push_back(&got.profit);
    for (A& e : got.uptake) outs.push_back(&e);
    return std::span<A* const>(outs.data(), outs.size());
  };
  if (preacc) {
    odelia::preaccumulate<A>(body, scratch, in, psi_soil, rh, rv);
  } else {
    body();
  }

  A* target = (which == 0) ? &got.profit : &got.uptake[which - 1];
  tape.registerOutput(*target);
  xad::derivative(*target) = 1.0;
  tape.computeAdjoints();
  Result r;
  r.ok = true;
  // The consumer sweeps ONE recording once per census metric, and the ladder
  // asserts two sweeps of it agree. That is the one thing this probe did not do.
  {
    std::vector<double> first;
    for (A& v : in) first.push_back(xad::derivative(v));
    tape.clearDerivatives();
    xad::derivative(*target) = 1.0;
    tape.computeAdjoints();
    double worst = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i) {
      const double now = xad::derivative(in[i]);
      const double sc = std::max(std::fabs(first[i]), std::fabs(now));
      if (sc > 0) worst = std::max(worst, std::fabs(first[i] - now) / sc);
    }
    r.resweep = worst;
  }
  r.profit = odelia::util::to_passive(got.profit);
  for (A& v : in) r.rows.push_back(xad::derivative(v));
  for (A& v : psi_soil) r.rows.push_back(xad::derivative(v));
  for (A& v : rh) r.rows.push_back(xad::derivative(v));
  for (A& v : rv) r.rows.push_back(xad::derivative(v));
  return r;
}

// Two regions on ONE tape, with a consumer reading both -- which is what plant
// does per cohort per stage per step, and the one thing the single-region checks
// above cannot see.
Result run_two(int layers, double psi0, bool preacc) {
  pl::Leaf l = set_up(layers, psi0);
  l.find_root_collar_psi();
  if (l.operating_point_kind() != pl::Leaf::OperatingPointKind::Interior)
    return {0, {}, false};
  double curv = 0.0;
  try { curv = l.marginal_collar_slope<K>(); } catch (...) { return {0, {}, false}; }
  if (!(curv < 0.0)) return {0, {}, false};

  Tape tape;
  pl::leaf_pars<A> in;
  const pl::leaf_pars<double> p = l.passive_pars();
  for (std::size_t i = 0; i < p.size(); ++i) in[i] = A(p[i]);
  const int n = l.supply_n_layers();
  std::vector<A> psi_soil, rh, rv;
  for (int i = 0; i < n; ++i) {
    psi_soil.push_back(A(l.roots_.psi_soil_[std::size_t(i)]));
    rh.push_back(A(l.roots_.network_.r_R_H_min[std::size_t(i)]));
    rv.push_back(A(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
  }
  for (A& v : in) tape.registerInput(v);
  for (A& v : psi_soil) tape.registerInput(v);
  for (A& v : rh) tape.registerInput(v);
  for (A& v : rv) tape.registerInput(v);
  tape.newRecording();

  const pl::SupplyAt<A> supply{psi_soil, rh, rv, in[pl::par_root_P50],
                               in[pl::par_root_c]};
  std::vector<double> scratch;
  A total = A(0.0);
  for (int rep = 0; rep < 2; ++rep) {
    pl::Leaf::LeafOutputs<A> got;
    static std::vector<A*> outs;
    auto body = [&]() -> std::span<A* const> {
      const auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), supply);
      const A collar = l.collar_at<K, A>(draw, in, curv);
      got = l.outputs_at<K, A>(collar, draw, in);
      outs.clear();
      outs.push_back(&got.profit);
      for (A& e : got.uptake) outs.push_back(&e);
      return std::span<A* const>(outs.data(), outs.size());
    };
    if (preacc) odelia::preaccumulate<A>(body, scratch, in, psi_soil, rh, rv);
    else body();
    // A consumer reading BOTH regions, as plant's rates do.
    total = total + got.profit * double(rep + 1);
    for (A& e : got.uptake) total = total + e * 1e6;
  }
  tape.registerOutput(total);
  xad::derivative(total) = 1.0;
  tape.computeAdjoints();
  Result r; r.ok = true; r.profit = odelia::util::to_passive(total);
  for (A& v : in) r.rows.push_back(xad::derivative(v));
  for (A& v : psi_soil) r.rows.push_back(xad::derivative(v));
  for (A& v : rh) r.rows.push_back(xad::derivative(v));
  for (A& v : rv) r.rows.push_back(xad::derivative(v));
  return r;
}

}  // namespace

int main() {
  std::printf("=== TWO regions on one tape, one consumer reading both ===\n");
  for (int layers : {1, 3}) {
    for (double psi0 = 0.10; psi0 < 1.2; psi0 += 0.25) {
      const Result a = run_two(layers, psi0, false);
      const Result b = run_two(layers, psi0, true);
      if (!a.ok || !b.ok) continue;
      double worst = 0.0; std::size_t at = 0;
      for (std::size_t i = 0; i < a.rows.size(); ++i) {
        const double sc = std::max(std::fabs(a.rows[i]), std::fabs(b.rows[i]));
        const double rel = sc > 0 ? std::fabs(a.rows[i] - b.rows[i]) / sc : 0.0;
        if (rel > worst) { worst = rel; at = i; }
      }
      std::printf("layers=%d psi0=%.2f  value plain=%.8e preacc=%.8e  worst row rel %.3e at %zu"
                  "  plain=%.6e preacc=%.6e  %s\n",
                  layers, psi0, a.profit, b.profit, worst, at, a.rows[at],
                  b.rows[at], worst < 1e-10 ? "ok" : "*** DIFFER ***");
    }
  }
  std::printf("\n=== single region (control) ===\n");
  const char* names[] = {"profit", "uptake[0]", "uptake[1]", "uptake[2]"};
  for (int layers : {1, 3}) {
    for (double psi0 = 0.10; psi0 < 1.2; psi0 += 0.25) {
      for (std::size_t w = 0; w < std::size_t(1 + layers) && w < 4; ++w) {
        const Result a = run(layers, psi0, false, w);
        const Result b = run(layers, psi0, true, w);
        if (!a.ok || !b.ok) continue;
        double worst = 0.0; std::size_t at = 0; std::size_t nz = 0;
        for (std::size_t i = 0; i < a.rows.size(); ++i) {
          const double s = std::max(std::fabs(a.rows[i]), std::fabs(b.rows[i]));
          if (a.rows[i] != 0.0) ++nz;
          const double rel = s > 0 ? std::fabs(a.rows[i] - b.rows[i]) / s : 0.0;
          if (rel > worst) { worst = rel; at = i; }
        }
        std::printf("layers=%d psi0=%.2f seed=%-9s rows=%zu nonzero=%zu "
                    "worst rel %.3e at %zu   plain=%.6e preacc=%.6e  %s\n",
                    layers, psi0, names[w], a.rows.size(), nz, worst, at,
                    a.rows[at], b.rows[at], worst < 1e-10 ? "ok" : "*** DIFFER ***");
        std::printf("    resweep: plain %.3e   preacc %.3e   %s\n", a.resweep,
                    b.resweep,
                    (b.resweep < 1e-12) ? "stable" : "*** SECOND SWEEP DIFFERS ***");
      }
    }
  }
  return 0;
}
