// Timing harness for the root-collar solve.
//
//   make -C tests/cpp bench && ./bench_solve [reps]
//
// Written to settle issue #2's open question -- whether the soil/root supply
// path has to be a *template* policy so its calls inline, or whether a plain
// composed class is fine. The answer is recorded in PLAN.md 7b; the harness is
// kept because the question recurs (issue #3 asks the same thing about lambda,
// and gets a different answer).
//
// Workload is 288 operating points over psi_soil, PPFD, VPD and 1/3/5 soil layers,
// all at a leaf temperature of 25 C, shutdown corner included -- the golden grid's
// state axes, at its reference temperature.
//
// ⚠️ It deliberately does NOT follow the golden grid's temperature axis. A timing
// baseline is only useful against its own history (tools/cost-baseline.tsv,
// tools/bench_history.sh), so doubling the workload would end that history to
// measure nothing new about the solve. If a temperature-dependent cost ever needs
// measuring, add a second workload rather than growing this one.
//
// Reports min-of-N, not mean. This is a deterministic computation on a noisy
// machine, so the minimum is the least-contaminated estimate of its cost; the
// mean mostly measures what else the laptop was doing. Run-to-run reproducibility
// at reps=2000 is about +/-0.01 us, and at reps=40 it is +/-0.5 us -- so use
// enough reps before believing a small difference.
//
// The printed checksum folds in every reported output. It exists so that a
// change which perturbed the arithmetic is caught here rather than being
// mistaken for a speed difference.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifndef PHYLLOPTIM_BENCH_LABEL
#define PHYLLOPTIM_BENCH_LABEL "solve"
#endif

namespace {

// Trait values and fixed drivers from plant's tests/testthat/test-leaf.r, the
// same ones test_golden.cpp uses.
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0;
const double kAreaLeaf = 0.05;
const double kCa = 40.0, kO2 = 21.0, kTleaf = 25.0, kPatm = 101.3;

struct Point {
  double psi_soil, ppfd, vpd;
  int layers;
  std::vector<double> ps, depth, root;
};

std::vector<Point> grid() {
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  const int layer_counts[] = {1, 3, 5};

  std::vector<Point> pts;
  for (double p : psi_soils) {
    for (double q : ppfds) {
      for (double d : vpds) {
        for (int n : layer_counts) {
          Point pt{p, q, d, n, {}, {}, {}};
          pt.ps.resize(n);
          pt.depth.resize(n);
          pt.root.resize(n);
          for (int i = 0; i < n; ++i) {
            pt.ps[i] = p + 0.25 * i;
            pt.depth[i] = 1.0 * (i + 1);
            // root carbon PER UNIT LEAF AREA -- set_physiology no longer takes
            // area_leaf, so the ratio is the input. Matches test_golden.cpp.
            pt.root[i] = 1.0 / n / kAreaLeaf;
          }
          pts.push_back(pt);
        }
      }
    }
  }
  return pts;
}

// One pass over the grid. A separate Leaf per point, as test_golden.cpp does --
// the shutdown-state leak (PLAN.md item 2) makes a reused Leaf order-dependent,
// and this harness should measure what the golden file pins.
//
// Spline setup is hoisted out of the timed region: setup_transpiration and
// setup_root_vulnerability are per-strategy work in plant, not per-solve, so
// timing them would dilute the signal. set_physiology is left in, because plant
// does call it once per solve and it invalidates the caches -- which is what
// keeps each timed solve a genuine cold solve rather than a memo hit.
double pass(std::vector<phylloptim::Leaf> &leaves, const std::vector<Point> &pts) {
  double checksum = 0.0;
  // ONE reused RootNetwork, refilled in place per point. This is deliberately
  // what plant does since #33 -- it holds the network as a TF24_Strategy member,
  // exactly as it used to hold the root-carbon buffer, and refills it through the
  // in-place overload. Building a fresh network per point instead measures
  // +0.074 us/solve of allocation that plant does not pay, which would make an
  // interleaved before/after comparison of this harness meaningless.
  phylloptim::RootNetwork net;
  for (size_t i = 0; i < pts.size(); ++i) {
    const Point &pt = pts[i];
    phylloptim::Leaf &l = leaves[i];
    phylloptim::root_network_from_carbon(pt.root,
                                        phylloptim::layer_thickness(pt.depth),
                                        fixture::beta_R_H, fixture::beta_R_V, net);
    l.set_physiology(net, pt.ppfd, pt.ps, pt.depth, kKs * kTheta / kH,
                     pt.vpd, kCa, kTleaf, kO2, kPatm);
    l.find_root_collar_psi();
    for (double v : {l.opt_psi_stem_, l.opt_root_psi_, l.ci_,
                     l.assim_colimited_, l.transpiration_, l.stom_cond_CO2_,
                     l.profit_, l.E_up_}) {
      if (std::isfinite(v)) {
        checksum += v;
      }
    }
  }
  return checksum;
}

} // namespace

int main(int argc, char **argv) {
  const int reps = argc > 1 ? std::atoi(argv[1]) : 2000;
  const std::vector<Point> pts = grid();

  std::vector<phylloptim::Leaf> leaves(pts.size());
  for (phylloptim::Leaf &l : leaves) {
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
  }

  double checksum = 0.0;
  double best = 1e300;
  for (int r = 0; r < reps; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    checksum = pass(leaves, pts);
    const auto t1 = std::chrono::steady_clock::now();
    best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
  }

  printf("%-14s  %8.2f us/solve   (%zu points, best of %d)   checksum %.17g\n",
         PHYLLOPTIM_BENCH_LABEL, best / pts.size() * 1e6, pts.size(), reps, checksum);
  return 0;
}
