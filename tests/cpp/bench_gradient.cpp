// Timing harness for a TRAIT GRADIENT, with no R in the way.
//
//   make -C tests/cpp bench_gradient && ./bench_gradient [reps]
//
// Written because a projection and a measurement disagreed. Computing dA/dtheta
// by the implicit function theorem -- one solve plus 2N cheap gradient
// evaluations instead of 2N re-solves -- projects 12x. Measured in the R layer
// the composite came out 6% SLOWER, because an R call costs ~1.8 us against the
// 0.26 us of C++ work it wraps.
//
// That is a statement about the R boundary and NOT about the composite, and the
// distinction matters because **plant links these headers directly**. So the
// arithmetic has to be done again with the boundary removed, which is what this
// measures. The composite does win here, by rather less than 12x, and for four
// of the traits by almost nothing at all -- for a reason neither the projection
// nor the R measurement anticipated.
//
// WHAT THE TWO ARMS ARE. Both perturb one trait by a relative 1e-6 either side
// and both pay set_traits + set_physiology per perturbation, because a trait
// change genuinely requires re-deriving the temperature-dependent block:
//
//   FD   2 x [set_traits, set_physiology, find_root_collar_psi]
//   IFT  2 x [set_traits, set_physiology, dprofit at the UNPERTURBED psi*,
//             evaluate_root_collar_psi at that same psi*]
//
// The second evaluation in the IFT arm is the direct term dA/dtheta|_psi, and it
// is the line item 11d's route-2 table omitted -- it costed a perturbation as one
// `dprofit` and nothing else. It is not optional: the composite is
// dA/dtheta = dA/dtheta|_psi + (dA/dpsi)(dpsi*/dtheta), and dropping the first
// term is not an approximation, it is a different quantity.
//
// This harness does NOT check the gradient. tests/testthat/test-gradient.R does
// that, against arbitrated references and over the whole golden grid, including
// the active-set cases this ignores entirely -- it runs at one interior point
// because it is measuring cost, not correctness. Do not read a speedup here as
// evidence that the composite is usable at a given operating point; that is what
// leaf_gradient()'s stationarity test is for.
//
// Reports min-of-N over interleaved rounds, for the reasons bench_solve.cpp
// gives at greater length: the minimum is the least-contaminated estimate, and
// sequential A-then-B is not good enough at this effect size. The printed
// checksum folds in every value both arms compute, so a change that perturbed
// the arithmetic is caught here rather than read as a speed difference.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

// The default trait vector, in set_traits' argument order.
//
// Eleven of the fifteen traits. `R_d_25`, `JS22_gamma`, `CMax_a` and `CMax_b` are
// left at their defaults because none needs a spline rebuild, so all four land in
// the cheap bucket this table already has examples of -- adding them would lengthen
// the output without changing the finding. `beta_R_H` and `beta_R_V` are not here at all -- they left
// with the root-architecture model in #33, and a gradient in either is now taken
// where the network is built, exactly, because root_network_from_carbon is
// homogeneous of degree 1 in each.
//
// ⚠️ ONE COUNT FOR ALL THREE USERS OF IT, and that is the fix for a real bug
// rather than a tidy-up. `Traits::v`, `kNames` and the timing loop each carried
// the number separately, and the loop's copy was a literal `13` against an
// 11-element array. So it read `kBase.v[11..12]` and `kNames[11..12]` out of
// bounds and handed the latter to `printf("%s")`.
//
// It survived every trait-count change since -- the vector has been 14, then 12,
// now 15 -- because nothing compared the three numbers. And it was invisible in
// the ordinary way: the read is only fatal when the address space happens to put
// something unmapped after the globals, so it ran clean under `lldb` (which
// disables ASLR) and segfaulted when run directly. It also *looked* like a crash
// during initialisation, because a piped stdout is block-buffered and the whole
// report was still sitting in the buffer when the process died.
//
// `-fsanitize=address` names it in one line: "global-buffer-overflow ... 0 bytes
// after global variable 'kBase'". Reach for that before reading a stack trace.
constexpr int kNTraits = 11;

struct Traits {
  double v[kNTraits];
};

const Traits kBase{{96.0, 2.680147, 3.4, 2.680147, 3.4, 1.5,
                    157.44, 0.30, 0.7, 0.99, 7.5}};

const char *kNames[kNTraits] = {"vcmax_25",  "stem_c",     "stem_P50",
                          "root_c",    "root_P50",   "TF24_beta2",
                          "jmax_25",   "a",          "curv_elec",
                          "curv_colim", "cost_scale"};

static_assert(sizeof(kBase.v) / sizeof(kBase.v[0]) == kNTraits,
              "kBase must have exactly kNTraits entries");
static_assert(sizeof(kNames) / sizeof(kNames[0]) == kNTraits,
              "kNames must have exactly kNTraits entries");

// Whether perturbing this trait forces a vulnerability spline to be rebuilt:
// stem_c/stem_P50 own the transpiration pair, root_c/root_P50 own the root curve.
// This turns out to be the single most important fact about the results.
bool rebuilds_a_spline(int i) {
  return i == 1 || i == 2 || i == 3 || i == 4;
}

void apply_traits(phylloptim::Leaf &l, const Traits &t) {
  // R_d_25 last, at the class default, so the benchmark measures the default
  // configuration.
  l.set_traits(t.v[0], t.v[1], t.v[2], t.v[3], t.v[4], t.v[5],
               t.v[6], t.v[7], t.v[8], t.v[9], t.v[10], 1.44,
               /*JS22_gamma=*/1.0, /*CMax_a=*/0.6, /*CMax_b=*/0.0);
}

// One interior operating point. Drivers from plant's tests/testthat/test-leaf.r,
// as everywhere else in this suite.
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_drivers(phylloptim::Leaf &l) {
  std::vector<double> root{1.0 / kAreaLeaf}, psi_soil{2.0}, depth{1.0};
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth, kKs * kTheta / kH, 2.0, 40.0,
                   25.0, 21.0, 101.3);
}

using clock_type = std::chrono::steady_clock;

double us_per(clock_type::time_point a, clock_type::time_point b, long n) {
  return std::chrono::duration<double, std::micro>(b - a).count() / double(n);
}

// The perturbation is varied by an inconsequential amount per iteration so that
// no part of either loop can be hoisted out of the timed region. bench_solve's
// header explains why that matters; the tell, when it goes wrong, is a
// suspiciously round zero.
double perturbed(const Traits &t, int idx, int sign, long r) {
  return t.v[idx] * (1.0 + double(sign) * 1e-6 * (1.0 + double(r) * 1e-9));
}

double fd_arm(phylloptim::Leaf &l, int idx, long reps, double &sink) {
  const auto t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    for (int sign = -1; sign <= 1; sign += 2) {
      Traits t = kBase;
      t.v[idx] = perturbed(kBase, idx, sign, r);
      apply_traits(l, t);
      set_drivers(l);
      l.find_root_collar_psi();
      sink += l.assim_colimited_;
    }
  }
  return us_per(t0, clock_type::now(), reps);
}

double ift_arm(phylloptim::Leaf &l, int idx, double psi_star, long reps,
               double &sink) {
  const auto t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    for (int sign = -1; sign <= 1; sign += 2) {
      Traits t = kBase;
      t.v[idx] = perturbed(kBase, idx, sign, r);
      apply_traits(l, t);
      set_drivers(l);
      // The mixed partial, and then the direct term at the same fixed collar.
      sink += l.dprofit_droot_collar_psi(psi_star);
      sink += l.evaluate_root_collar_psi(psi_star);
      sink += l.assim_colimited_;
    }
  }
  return us_per(t0, clock_type::now(), reps);
}

// The same composite, but moving the stem curve by the homogeneity rescale
// instead of through set_traits -- so no spline is rebuilt and set_physiology is
// not needed either, since nothing it derives depends on stem_b.
double ift_arm_rescaled_stem_curve(phylloptim::Leaf &l, double psi_star,
                                   long reps, double &sink) {
  const auto t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    for (int sign = -1; sign <= 1; sign += 2) {
      l.perturb_stem_P50(perturbed(kBase, 2, sign, r));
      sink += l.dprofit_droot_collar_psi(psi_star);
      sink += l.evaluate_root_collar_psi(psi_star);
      sink += l.assim_colimited_;
    }
  }
  const double us = us_per(t0, clock_type::now(), reps);
  apply_traits(l, kBase);  // set_traits is the way back; it forces the rebuild
  set_drivers(l);
  return us;
}

}  // namespace

int main(int argc, char **argv) {
  const long reps = argc > 1 ? std::atol(argv[1]) : 3000;
  const long per_round = reps / 3 > 0 ? reps / 3 : 1;

  phylloptim::Leaf l;
  apply_traits(l, kBase);
  set_drivers(l);
  l.find_root_collar_psi();
  const double psi_star = l.opt_root_psi_;
  double sink = 0.0;

  // Component costs first, because they are what explain the table below rather
  // than merely accompanying it.
  {
    auto t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      set_drivers(l);
      sink += l.assim_max_;
    }
    printf("  set_physiology              %8.3f us\n", us_per(t0, clock_type::now(), reps));

    t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      l.find_root_collar_psi();
      sink += l.opt_root_psi_;
    }
    printf("  find_root_collar_psi       %8.3f us\n", us_per(t0, clock_type::now(), reps));

    t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      sink += l.dprofit_droot_collar_psi(psi_star + double(r) * 1e-15);
    }
    printf("  dprofit_droot_collar_psi   %8.3f us\n", us_per(t0, clock_type::now(), reps));

    t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      sink += l.evaluate_root_collar_psi(psi_star + double(r) * 1e-15);
    }
    printf("  evaluate_root_collar_psi   %8.3f us\n", us_per(t0, clock_type::now(), reps));

    Traits t = kBase;
    t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      t.v[0] = kBase.v[0] + double(r) * 1e-9;  // vcmax_25: no spline
      apply_traits(l, t);
      sink += l.vcmax_25;
    }
    printf("  set_traits, no rebuild     %8.3f us\n", us_per(t0, clock_type::now(), reps));

    t = kBase;
    t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      t.v[2] = kBase.v[2] + double(r) * 1e-9;  // stem_b: rebuilds the stem curve
      apply_traits(l, t);
      sink += l.stem_b;
    }
    printf("  set_traits, spline rebuild %8.3f us   <-- the finding\n",
           us_per(t0, clock_type::now(), reps));

    // And the way round it, for the stem curve only: it is homogeneous of degree
    // 1 in stem_b, which is a fixed multiple of stem_P50, so the spline for a
    // perturbed trait is this one with its argument rescaled.
    t0 = clock_type::now();
    for (long r = 0; r < reps; ++r) {
      l.perturb_stem_P50(kBase.v[2] + double(r) * 1e-9);
      sink += l.stem_b;
    }
    printf("  perturb_stem_P50, no rebuild %8.3f us  <-- and the way round it\n",
           us_per(t0, clock_type::now(), reps));
    apply_traits(l, kBase);
    set_drivers(l);
    l.find_root_collar_psi();
  }

  printf("\n%-15s %10s %10s %9s  %s\n", "trait", "FD", "IFT", "speedup",
         "rebuilds a spline");
  double fd_all = 0.0, ift_all = 0.0, fd_cheap = 0.0, ift_cheap = 0.0;
  int n_cheap = 0;
  for (int i = 0; i < kNTraits; ++i) {
    double fd = 1e300, ift = 1e300;
    for (int round = 0; round < 3; ++round) {
      // Reset to the base state before each arm so neither inherits the other's
      // caches, then interleave -- see bench_solve.cpp on why sequential is not
      // good enough.
      apply_traits(l, kBase);
      set_drivers(l);
      l.find_root_collar_psi();
      const double a = fd_arm(l, i, per_round, sink);
      apply_traits(l, kBase);
      set_drivers(l);
      l.find_root_collar_psi();
      const double b = ift_arm(l, i, psi_star, per_round, sink);
      if (a < fd) {
        fd = a;
      }
      if (b < ift) {
        ift = b;
      }
    }
    printf("%-15s %10.3f %10.3f %8.2fx  %s\n", kNames[i], fd, ift, fd / ift,
           rebuilds_a_spline(i) ? "yes" : "");
    fd_all += fd;
    ift_all += ift;
    if (!rebuilds_a_spline(i)) {
      fd_cheap += fd;
      ift_cheap += ift;
      ++n_cheap;
    }
  }

  printf("\n  %d of the %d traits  FD %8.1f us   IFT %8.1f us   %.2fx\n",
         kNTraits, phylloptim::gradient::n_traits,
         fd_all, ift_all, fd_all / ift_all);
  printf("  the %2d with no rebuild  FD %8.1f us   IFT %8.1f us   %.2fx\n",
         n_cheap, fd_cheap, ift_cheap, fd_cheap / ift_cheap);
  printf("  the  4 that rebuild     FD %8.1f us   IFT %8.1f us   %.2fx\n",
         fd_all - fd_cheap, ift_all - ift_cheap,
         (fd_all - fd_cheap) / (ift_all - ift_cheap));
  // stem_b again, moved by the rescale rather than by a rebuild. Interleaved with
  // its own rebuild arm rather than compared against the table above, because the
  // two have to be measured in one run to be comparable at this effect size.
  double reb = 0.0, resc = 0.0;
  for (int round = 0; round < 3; ++round) {
    const double a = ift_arm(l, 2, psi_star, per_round, sink);
    const double b = ift_arm_rescaled_stem_curve(l, psi_star, per_round, sink);
    reb = round == 0 ? a : std::min(reb, a);
    resc = round == 0 ? b : std::min(resc, b);
  }
  printf("\n  stem_b, IFT: rebuild %8.3f us   rescaled %8.3f us   %.2fx\n",
         reb, resc, reb / resc);

  // --- and the same identity per OBSERVATION, which is the number to quote ----
  //
  // ⚠️ THE TWO FIGURES ARE NOT THE SAME MEASUREMENT, and #74 exists because one
  // was quoted for the other. Above is per PERTURBATION: two of them, and nothing
  // else. Here is per observation through `gradient::batch`, which is what a
  // calibration actually pays -- the perturbations plus one solve, one restore and
  // the fixed cost of the composite. Say which one you mean.
  {
    const double kmax = kKs * kTheta / kH;
    double theta[phylloptim::gradient::n_pars] = {
        kBase.v[0], kBase.v[1], kBase.v[2],  kBase.v[3], kBase.v[4],
        kBase.v[5], kBase.v[6], kBase.v[7],  kBase.v[8], kBase.v[9],
        // ⚠️ POSITIONAL and silently short-fillable -- see test_leaf.cpp's copy of
        // this array for what adding a trait did to it. Count against `n_pars`.
        kBase.v[10], 1.44,      /*JS22_gamma=*/1.0,
        /*CMax_a=*/0.6, /*CMax_b=*/0.0, kmax, 0.0};

    std::vector<double> root{1.0 / kAreaLeaf}, psi_soil{2.0}, depth{1.0};
    phylloptim::gradient::Drivers gd;
    gd.root_network = fixture::root_network(root, depth);
    gd.PPFD = 900.0;
    gd.psi_soil = psi_soil;
    gd.soil_depth = depth;
    gd.atm_vpd = 2.0;
    gd.ca = 40.0;
    gd.leaf_temp = 25.0;
    gd.atm_o2_kpa = 21.0;
    gd.atm_kpa = 101.3;
    const std::size_t n_obs = 20;
    const std::vector<phylloptim::gradient::Drivers> obs(n_obs, gd);

    struct Arm {
      const char *label;
      int par;
      bool fast;
    };
    const Arm arms[] = {
        {"vcmax_25 (no spline anywhere)", 0, true},
        {"stem_P50, fast_stem_curve on", phylloptim::gradient::par_stem_P50, true},
        {"stem_P50, fast_stem_curve off", phylloptim::gradient::par_stem_P50, false},
        {"stem_c (no identity to use)", phylloptim::gradient::par_stem_c, true}};

    printf("\n%-32s %10s %12s\n", "pars", "us/obs", "rebuilds/obs");
    const long batch_reps = std::max(1L, per_round / 20);
    for (const Arm &arm : arms) {
      phylloptim::gradient::Settings s;
      s.fast_stem_curve = arm.fast;
      const int pars[1] = {arm.par};
      double best = 1e300;
      long builds = 0;
      for (int round = 0; round < 3; ++round) {
        apply_traits(l, kBase);
        set_drivers(l);
        l.find_root_collar_psi();
        const long before = l.stem_curve_builds_;
        const auto t0 = clock_type::now();
        for (long r = 0; r < batch_reps; ++r) {
          // Vary theta per rep for the reason `perturbed()` documents -- but in
          // `TF24_beta2`, NOT in the parameter being differentiated. Moving the base
          // `stem_b` per rep moves it off `stem_b_spline_`, which is a rebuild the
          // model genuinely owes and which showed up as 0.05 rebuilds/obs (one per
          // rep over 20 observations) while reading like a leftover of the bug
          // this arm exists to measure. `TF24_beta2` owns no spline.
          theta[phylloptim::gradient::par_TF24_beta2] =
              kBase.v[7] * (1.0 + double(r) * 1e-12);
          const std::vector<phylloptim::gradient::Result> g =
              phylloptim::gradient::batch(l, theta, 1, obs, false, pars, 1, s);
          sink += g[0].grad[0];
        }
        const double us =
            us_per(t0, clock_type::now(), long(n_obs) * batch_reps);
        best = std::min(best, us);
        builds = l.stem_curve_builds_ - before;
      }
      printf("%-32s %10.2f %12.2f\n", arm.label, best,
             double(builds) / double(long(n_obs) * batch_reps));
    }
    printf(
        "\nstem_b's per-observation cost used to carry one rebuild that no\n"
        "perturbation asked for: the end-of-observation restore went through\n"
        "set_traits while the shortcut still had stem_b displaced, and that is\n"
        "exactly the case set_traits is required to rebuild for. `apply()` now\n"
        "undoes the displacement with the shortcut first, so the restore has\n"
        "nothing to rebuild -- rebuilds/obs above is the assertion. #74.\n");
  }

  printf(
      "\nThe four that rebuild a spline are stem_b, stem_c, root_b and root_c --\n"
      "and they are exactly the traits for which the argmax-mediated term is\n"
      "100%% of the gradient. So the composite's advantage is smallest precisely\n"
      "where the composite is most necessary.\n"
      "\nstem_b escapes it: G is homogeneous of degree 1 in stem_b, so the spline\n"
      "for a perturbed stem_b is this one with its argument rescaled and no\n"
      "rebuild is needed. stem_c has no such identity.\n");
  printf("\nchecksum %.6f\n", sink);
  return 0;
}
