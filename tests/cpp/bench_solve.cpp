// Timing harness for the root-collar solve.
//
//   make -C tests/cpp bench && ./bench_solve [reps]
//
// Written to settle issue #2's open question -- whether the soil/root supply
// path has to be a *template* policy so its calls inline, or whether a plain
// composed class is fine. A plain composed class is; the harness is
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
// the shutdown-state leak makes a reused Leaf order-dependent,
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

// ---------------------------------------------------------------------------
// The single-layer optimisers
// ---------------------------------------------------------------------------
//
// A SECOND workload rather than a bigger one, which is what the warning at the
// top of this file asks for: the collar arm's 288 points and its `us/solve` line
// are a history (tools/cost-baseline.tsv), and growing them would end that
// history to measure something else.
//
// It runs on the 1-LAYER subset of the same grid, because every one of these
// refuses a multi-layer supply. Same drivers, same traits, so the arms are
// directly comparable per call, and against the collar line above.
//
// ⚠️ THESE LINES MUST NOT PRINT `us/solve`. tools/bench_history.sh greps every
// occurrence of that exact string out of this program's whole output and assigns
// the result to ONE TSV field, so a second matching line silently corrupts the
// history file it is building. Hence `us/call` here, and hence this comment
// rather than a tidier-looking unit.

// ⚠️ THE `Closed` ARMS ARE IN THE SAME BINARY AS THEIR OWN BASELINES, WHICH IS
// WHAT MAKES THE COMPARISON LEGITIMATE. Hazard 5: between processes the noise on
// this harness is ~+/-0.1 us against +/-0.01 us within one, so timing the exact
// solve in one build and the closed form in another is not good enough to price a
// method against the solve it replaces. Both run in one process, on one grid, in
// one loop -- so the ratio of two printed lines is a controlled A/B.
enum class Arm { TF, ProfitMax, CF77, TFClosed, CF77Closed };

const char *arm_label(Arm a) {
  switch (a) {
    case Arm::TF:            return "psi_stem:TF";
    case Arm::ProfitMax:     return "psi_stem:ProfitMax";
    case Arm::CF77: return "psi_stem:CF77";
    case Arm::TFClosed:      return "psi_stem:TF-closed";
    case Arm::CF77Closed:    return "psi_stem:CF77-closed";
  }
  return "psi_stem:?";
}

// Cowan-Farquhar consumes a PRESCRIBED lambda and throws without one. Fixed here
// rather than taken from a preceding solve, which would time two solves and call
// it one.
const double kLambdaCF77 = 1.5e5;

double pass_optimiser(Arm arm, std::vector<phylloptim::Leaf> &leaves,
                      const std::vector<Point> &pts) {
  double checksum = 0.0;
  phylloptim::RootNetwork net;
  for (size_t i = 0; i < pts.size(); ++i) {
    const Point &pt = pts[i];
    phylloptim::Leaf &l = leaves[i];
    phylloptim::root_network_from_carbon(pt.root,
                                        phylloptim::layer_thickness(pt.depth),
                                        fixture::beta_R_H, fixture::beta_R_V, net);
    l.set_physiology(net, pt.ppfd, pt.ps, pt.depth, kKs * kTheta / kH,
                     pt.vpd, kCa, kTleaf, kO2, kPatm);
    switch (arm) {
      case Arm::TF:        l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();        break;
      case Arm::ProfitMax: l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>(); break;
      case Arm::CF77:
                           l.CF77_lambda_ = kLambdaCF77;
                           l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>(); break;
      // ⚠️ THROUGH `optimise()`, NOT THE FREE FUNCTION, deliberately: what a caller
      // pays for the closed form INCLUDES the fallback to the exact solve on the
      // rows its guard rejects, and pricing the fast path alone would quote the
      // ceiling rather than the realised cost. The printed fallback fraction is
      // what makes the two readable apart.
      case Arm::TFClosed:
                           l.set_model(phylloptim::Leaf::CostCurve::TF24, false, true);
                           l.optimise(); break;
      case Arm::CF77Closed:
                           l.CF77_lambda_ = kLambdaCF77;
                           l.set_model(phylloptim::Leaf::CostCurve::CF77, false, true);
                           l.optimise(); break;
    }
    for (double v : {l.opt_psi_stem_, l.ci_, l.assim_colimited_,
                     l.transpiration_, l.stom_cond_CO2_, l.profit_}) {
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

  printf("%-22s  %8.2f us/solve   (%zu points, best of %d)   checksum %.17g\n",
         PHYLLOPTIM_BENCH_LABEL, best / pts.size() * 1e6, pts.size(), reps, checksum);

  // --- the single-layer optimisers, on the 1-layer subset -------------------
  std::vector<Point> one;
  for (const Point &pt : pts) {
    if (pt.layers == 1) {
      one.push_back(pt);
    }
  }
  std::vector<phylloptim::Leaf> one_leaves(one.size());
  for (phylloptim::Leaf &l : one_leaves) {
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
  }
  for (Arm arm : {Arm::TF, Arm::ProfitMax, Arm::CF77, Arm::TFClosed,
                  Arm::CF77Closed}) {
    double arm_checksum = 0.0;
    double arm_best = 1e300;
    for (int r = 0; r < reps; ++r) {
      for (phylloptim::Leaf &l : one_leaves) {
        l.reset_closed_form_counters();
      }
      const auto t0 = std::chrono::steady_clock::now();
      arm_checksum = pass_optimiser(arm, one_leaves, one);
      const auto t1 = std::chrono::steady_clock::now();
      arm_best = std::min(arm_best, std::chrono::duration<double>(t1 - t0).count());
    }
    // The fallback fraction phi, summed over the grid's leaves, because the
    // realised speedup is 1/[phi + (1-phi)/ceiling] and the ceiling on its own is
    // not a number anyone should quote (closed_form.hpp note 3).
    int calls = 0, backs = 0;
    for (const phylloptim::Leaf &l : one_leaves) {
      calls += l.closed_form_calls_;
      backs += l.closed_form_fallbacks_;
    }
    char phi[48] = "";
    if (calls > 0) {
      snprintf(phi, sizeof phi, "   phi %.3f (%d/%d)",
               double(backs) / double(calls), backs, calls);
    }
    printf("%-22s  %8.2f us/call    (%zu points, best of %d)   checksum %.17g%s\n",
           arm_label(arm), arm_best / one.size() * 1e6, one.size(), reps,
           arm_checksum, phi);
  }
  return 0;
}
