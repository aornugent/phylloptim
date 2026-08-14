// Golden-file regression test over a grid of operating points.
//
//   make -C tests/cpp golden        # regenerate golden/operating_points.tsv
//   make -C tests/cpp && ./test_golden
//
// Why this exists. The refactors in PLAN.md items 7-11 are meant to be
// behaviour-preserving, and the only way to know is to pin the behaviour first.
// PLAN.md item 1 -- cross-checking against plant's compiled build -- is the real
// validation, and it is DONE: the swap was bit-identical, including 78 of 78 SCM
// nodes. This file is the cheaper, always-available version: it freezes what THIS
// implementation produces so that a
// refactor which changes any of it fails loudly.
//
// Comparison is bit-exact by default. Values are written with %.17g, which
// round-trips an IEEE double exactly, so a passing run means the refactor did not
// perturb a single floating-point operation. If a change is *meant* to alter
// results, regenerate deliberately and say so in the commit.
//
// Bit-exactness holds on the platform that generated the file (macOS/arm64) and
// cannot hold on any other -- `--rtol` is for those. See the comment above main().
//
// Note: a fresh Leaf is constructed for every grid point. That is not for tidiness
// -- it WAS required, because the shutdown-state leak (PLAN.md item 2) made a
// reused Leaf order-dependent, which would have made this file ill-defined. That
// leak is fixed, so the construction is now belt-and-braces rather than load
// bearing -- and it is kept, because a fresh Leaf per point is also what makes
// this file blind to stale state by construction (hazard 8).

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

const char *kGoldenPath = "golden/operating_points.tsv";

struct Row {
  // inputs
  double psi_soil, ppfd, vpd, leaf_temp;
  int layers;
  // outputs
  double psi_stem, opt_root_psi, ci, assim, transpiration, gc, profit, e_up, uptake;
  // WHICH KIND of operating point the solve found. Deliberately NOT written to
  // the golden file: the file is compared bit-exactly and a new column would
  // force a regeneration, which is the one thing this file must not need for a
  // pure addition. It is checked by count instead (see check_operating_kinds).
  phylloptim::Leaf::OperatingPointKind kind;
};

// Trait values and fixed drivers from plant's tests/testthat/test-leaf.r.
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0;
const double kAreaLeaf = 0.05;
const double kCa = 40.0, kO2 = 21.0, kPatm = 101.3;
// ⚠️ LEAF TEMPERATURE IS A GRID AXIS, NOT A FIXED DRIVER. It was fixed at 25 C, and
// that made this file blind to every temperature response in the model: the
// reference values of Vcmax, Jmax and R_d are DEFINED at 25 C, so a change to any
// response curve is inert there BY CONSTRUCTION.
//
// Two temperatures are enough to see one: the reference, and 40 C where the response
// bites hardest. The extra pinned rows at the hot end are a feature rather than
// noise.
const double kLeafTemps[] = {25.0, 40.0};

Row solve(double psi_soil, double ppfd, double vpd, int layers,
          double leaf_temp) {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);

  // Spread the soil profile over `layers` equal 1 m layers, drying with depth so
  // that multi-layer runs are not just a repeated single layer, and split root
  // carbon evenly.
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers / kAreaLeaf;
  }

  l.set_physiology(fixture::root_network(root, depth), ppfd, ps, depth, kKs * kTheta / kH, vpd, kCa,
                   leaf_temp, kO2, kPatm);
  l.find_root_collar_psi();

  double uptake = 0.0;
  for (double s : l.soil_consumption_) {
    if (std::isfinite(s)) {
      uptake += s;
    }
  }
  return Row{psi_soil,        ppfd,   vpd,        leaf_temp,  layers,
             l.opt_psi_stem_,
             l.opt_root_psi_, l.ci_, l.assim_colimited_, l.transpiration_,
             l.stom_cond_CO2_,   l.profit_, l.E_up_,        uptake,
             l.operating_point_kind()};
}

std::vector<Row> run_grid() {
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  const int layer_counts[] = {1, 3, 5};

  // ⚠️ TEMPERATURE IS THE OUTERMOST LOOP ON PURPOSE. It makes the file contiguous
  // copies of one 288-point state grid, so a regeneration that only adds a
  // temperature is checkable as an ADDITION: the 25 C block must come out
  // byte-identical. Reordering these loops rewrites every row without changing any
  // value.
  std::vector<Row> rows;
  for (double t : kLeafTemps) {
    for (double p : psi_soils) {
      for (double q : ppfds) {
        for (double d : vpds) {
          for (int n : layer_counts) {
            rows.push_back(solve(p, q, d, n, t));
          }
        }
      }
    }
  }
  return rows;
}

// ---------------------------------------------------------------------------
// The classification of the grid, by count
// ---------------------------------------------------------------------------
//
// `Leaf::operating_point_kind()` says which branch of the collar solve produced
// the point, and the grid's split between those branches is already an
// established number AT 25 C: 240 of those 288 points are feasible, of which 198 sit
// at an interior profit maximum and 42 at a constrained optimum pinned to a bracket
// bound (24 wet, 18 dry); the remaining 48 shut down. Those figures are quoted in
// the developer guide and in maximise_profit_over_collar's own comment. The 40 C
// counts were measured when the axis was added.
//
// So this asserts the classification against numbers that already exist rather
// than inventing reference values for it. It is a cheap test with a specific
// job: the tag must come from the branch that was TAKEN, and a plausible-looking
// wrong implementation -- classifying by |dprofit| against a tolerance -- would
// not reproduce this split, because dprofit's shut-down sentinel is exactly 0.0
// and would move all 48 shutdown points into `interior`.
//
// Three kinds are expected to be EMPTY here: `shade-death` because assim_max_ never
// gets there (it is reached by light, not by drying, and not by heat either at 40 C);
// `solver-refused` and `non-finite-gradient` because both bracket endpoints admit a
// usable gradient on every feasible row.
//
// ⚠️ The hot end does NOT reach the compensation-point exit either, which is worth
// knowing before reading the zeros as coverage: at the defaults that needs about
// 45 C. `test_rd_temperature_response` covers it.
//
// A zero counts for this grid, whose psi_soil values are {0.5, 1, 2, 3, 4, 6},
// and says nothing about the model. What the two failure branches do is checked
// by test_collar_solve_refuses_rather_than_guessing in test_leaf.cpp.
// ⚠️ COUNTED PER TEMPERATURE, not over the whole grid: a single total would let
// points move between branches as the leaf warms and still add up.
struct KindCount {
  phylloptim::Leaf::OperatingPointKind kind;
  int expected;
};

// One row per entry of kLeafTemps, in that order.
struct TempKinds {
  double leaf_temp;
  int interior, pinned_wet, pinned_dry, shutdown;
};
//
// The split MOVES with temperature, in a direction that is physical: a hot leaf
// assimilates less and respires more, so it has less to gain from water and its
// optimum presses against the WET bound instead of the dry one. The 48 shut-down
// rows do not move, because there it is hydraulics rather than heat that forbids
// transpiration.
const TempKinds kExpectedKinds[] = {
    {25.0, 198, 24, 18, 48},
    {40.0, 160, 80, 0, 48},
};

int check_operating_kinds(const std::vector<Row> &rows) {
  using Kind = phylloptim::Leaf::OperatingPointKind;
  int failures = 0;
  int classified = 0;
  for (const TempKinds &e : kExpectedKinds) {
    const KindCount expected[] = {
        {Kind::Interior, e.interior},
        {Kind::PinnedWet, e.pinned_wet},
        {Kind::PinnedDry, e.pinned_dry},
        {Kind::HydraulicShutdown, e.shutdown},
        {Kind::Determined, 0},        {Kind::ShadeDeath, 0},
        {Kind::Prescribed, 0},        {Kind::SolverRefused, 0},
        {Kind::NonFiniteGradient, 0}, {Kind::Unsolved, 0},
    };
    for (const KindCount &k : expected) {
      int got = 0;
      for (const Row &r : rows) {
        if (r.leaf_temp == e.leaf_temp && r.kind == k.kind) {
          ++got;
        }
      }
      classified += got;
      if (got != k.expected) {
        ++failures;
        fprintf(stderr, "FAIL kinds at T=%g: %s got %d, expected %d\n",
                e.leaf_temp,
                phylloptim::Leaf::operating_point_kind_name(k.kind), got,
                k.expected);
      }
    }
  }
  if (classified != static_cast<int>(rows.size())) {
    ++failures;
    fprintf(stderr, "FAIL kinds: %d of %zu points classified\n", classified,
            rows.size());
  }
  if (failures == 0) {
    printf("kinds: %zu points at %zu leaf temperatures\n", rows.size(),
           sizeof kExpectedKinds / sizeof kExpectedKinds[0]);
    for (const TempKinds &e : kExpectedKinds) {
      printf("  T = %4.1f C  %3d interior, %2d pinned (%2d wet, %2d dry), "
             "%2d shutdown\n",
             e.leaf_temp, e.interior, e.pinned_wet + e.pinned_dry, e.pinned_wet,
             e.pinned_dry, e.shutdown);
    }
  }
  return failures == 0 ? 0 : 1;
}

const char *kHeader =
    "psi_soil\tppfd\tvpd\tleaf_temp\tlayers\tpsi_stem\topt_root_psi\tci\tassim\t"
    "transpiration\tgc\tprofit\te_up\tuptake\n";

void write_row(FILE *f, const Row &r) {
  fprintf(f, "%.17g\t%.17g\t%.17g\t%.17g\t%d", r.psi_soil, r.ppfd, r.vpd,
          r.leaf_temp, r.layers);
  for (double v : {r.psi_stem, r.opt_root_psi, r.ci, r.assim, r.transpiration, r.gc,
                   r.profit, r.e_up, r.uptake}) {
    fprintf(f, "\t%.17g", v);
  }
  fprintf(f, "\n");
}

int generate() {
  const std::vector<Row> rows = run_grid();
  FILE *f = fopen(kGoldenPath, "w");
  if (f == nullptr) {
    fprintf(stderr, "cannot write %s (run from tests/cpp/)\n", kGoldenPath);
    return 1;
  }
  fputs(kHeader, f);
  for (const Row &r : rows) {
    write_row(f, r);
  }
  fclose(f);
  printf("wrote %zu operating points to %s\n", rows.size(), kGoldenPath);
  return 0;
}

// Exact equality, with NaN treated as equal to NaN -- some grid points shut down
// and legitimately produce the NA sentinel (see PLAN.md item 2).
bool same(double got, double want) {
  if (std::isnan(got) && std::isnan(want)) {
    return true;
  }
  return got == want;
}

// Relative difference, used both to decide the tolerant comparison and to report
// how far apart the two builds actually are. Returns 0 for the NaN/NaN and
// exactly-equal cases, and infinity when one side is NaN and the other is not.
double rel_diff(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) {
    return 0.0;
  }
  if (std::isnan(a) != std::isnan(b)) {
    return std::numeric_limits<double>::infinity();
  }
  if (a == b) {
    return 0.0;
  }
  const double scale = std::max(std::abs(a), std::abs(b));
  if (scale == 0.0) {
    return 0.0;
  }
  return std::abs(a - b) / scale;
}

// Cross-platform tolerances, per field. The two classes are three orders of
// magnitude apart and the reason is structural, not arbitrary -- see the comment
// above main(). Negative means bit-exact.
struct Tolerance {
  double maximum; // profit: the value AT the optimum
  double argmax;  // everything else: evaluated at the optimum's LOCATION
};
const Tolerance kExact{-1.0, -1.0};
// Measured over the full grid on Linux (gcc worst; clang is smaller): profit
// 1.85e-6, argmax-derived 5.53e-4. That leaves 5.4x and 9.0x headroom.
//
// The profit tolerance is deliberately not loosened further despite the modest
// headroom: 1e-4 is the scale at which a real behavioural change shows (PLAN.md
// item 1), so a profit tolerance approaching it would gate nothing. If a future
// toolchain pushes profit past 1e-5, that is worth looking at rather than
// auto-passing.
const Tolerance kCrossPlatform{1e-5, 5e-3};

// `profit` is the only reported field that is the maximum itself; every other
// field is evaluated at the argmax and so inherits its displacement.
bool is_maximum_field(const char *name) {
  return std::strcmp(name, "profit") == 0;
}

int compare(Tolerance tol) {
  FILE *f = fopen(kGoldenPath, "r");
  if (f == nullptr) {
    fprintf(stderr,
            "FAIL: %s is missing. Generate it with `make golden` and commit it.\n",
            kGoldenPath);
    return 1;
  }
  char line[4096];
  if (fgets(line, sizeof line, f) == nullptr) {
    fprintf(stderr, "FAIL: %s is empty\n", kGoldenPath);
    fclose(f);
    return 1;
  }

  const std::vector<Row> rows = run_grid();
  int failures = 0;
  int inexact = 0;
  double worst_rel = 0.0;
  double worst_max_rel = 0.0, worst_argmax_rel = 0.0;
  double worst_got = 0.0, worst_want = 0.0;
  const char *worst_desc = "";
  Row worst_row{};
  size_t i = 0;
  for (; i < rows.size(); ++i) {
    if (fgets(line, sizeof line, f) == nullptr) {
      fprintf(stderr, "FAIL: golden file has only %zu rows, grid has %zu\n", i,
              rows.size());
      ++failures;
      break;
    }
    Row g{};
    // 5 inputs then 9 outputs
    const int n = sscanf(
        line,
        "%lg\t%lg\t%lg\t%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
        &g.psi_soil, &g.ppfd, &g.vpd, &g.leaf_temp, &g.layers, &g.psi_stem,
        &g.opt_root_psi, &g.ci, &g.assim, &g.transpiration, &g.gc, &g.profit,
        &g.e_up, &g.uptake);
    if (n != 14) {
      fprintf(stderr, "FAIL: row %zu is malformed (%d fields)\n", i, n);
      ++failures;
      continue;
    }
    const Row &r = rows[i];
    struct Field {
      const char *name;
      double got, want;
    };
    const Field fields[] = {{"psi_stem", r.psi_stem, g.psi_stem},
                            {"opt_root_psi", r.opt_root_psi, g.opt_root_psi},
                            {"ci", r.ci, g.ci},
                            {"assim", r.assim, g.assim},
                            {"transpiration", r.transpiration, g.transpiration},
                            {"gc", r.gc, g.gc},
                            {"profit", r.profit, g.profit},
                            {"e_up", r.e_up, g.e_up},
                            {"uptake", r.uptake, g.uptake}};
    for (const Field &fd : fields) {
      const double rd = rel_diff(fd.got, fd.want);
      const bool is_max = is_maximum_field(fd.name);
      const double rtol = is_max ? tol.maximum : tol.argmax;
      // Track the worst deviation in each class over the whole grid regardless
      // of pass/fail, so the summary reports magnitudes rather than leaving the
      // reader to infer them from a list truncated at 20 lines.
      double &worst_in_class = is_max ? worst_max_rel : worst_argmax_rel;
      if (rd > worst_in_class) {
        worst_in_class = rd;
      }
      if (rd > worst_rel) {
        worst_rel = rd;
        worst_desc = fd.name;
        worst_row = r;
        worst_got = fd.got;
        worst_want = fd.want;
      }
      if (!same(fd.got, fd.want)) {
        ++inexact;
      }
      const bool ok = rtol < 0.0 ? same(fd.got, fd.want) : (rd <= rtol);
      if (!ok) {
        if (failures < 20) {
          fprintf(stderr,
                  "FAIL psi_soil=%g ppfd=%g vpd=%g T=%g layers=%d %s: got %.17g, "
                  "want %.17g  (rel %.3g)\n",
                  r.psi_soil, r.ppfd, r.vpd, r.leaf_temp, r.layers, fd.name,
                  fd.got, fd.want, rd);
        }
        ++failures;
      }
    }
  }
  if (fgets(line, sizeof line, f) != nullptr) {
    fprintf(stderr, "FAIL: golden file has more rows than the grid (%zu)\n", i);
    ++failures;
  }
  fclose(f);

  // Always report the worst deviation. On a bit-exact pass it is 0 and says so;
  // on a tolerant pass it is the number that matters, and watching it drift is
  // the point of running the tolerant mode at all.
  char worst[256] = "none (bit-identical)";
  if (worst_rel > 0.0) {
    snprintf(worst, sizeof worst,
             "%.3g  at psi_soil=%g ppfd=%g vpd=%g T=%g layers=%d %s "
             "(got %.17g, want %.17g)",
             worst_rel, worst_row.psi_soil, worst_row.ppfd, worst_row.vpd,
             worst_row.leaf_temp, worst_row.layers, worst_desc, worst_got,
             worst_want);
  }

  const bool exact_mode = tol.argmax < 0.0;

  if (failures == 0) {
    if (exact_mode) {
      printf("golden: %zu operating points, all bit-identical\n", rows.size());
    } else {
      printf("golden: %zu operating points within cross-platform tolerance\n"
             "  %d of %zu values differ. Worst by class:\n"
             "    profit  (the maximum)      %.3g   tolerance %.1g\n"
             "    others  (from the argmax)  %.3g   tolerance %.1g\n",
             rows.size(), inexact, rows.size() * 9, worst_max_rel, tol.maximum,
             worst_argmax_rel, tol.argmax);
    }
    return 0;
  }
  fprintf(stderr,
          "\ngolden: %d mismatches over %zu operating points%s.\n"
          "  worst relative difference %s\n",
          failures, rows.size(),
          exact_mode ? "" : " beyond cross-platform tolerance", worst);
  if (!exact_mode) {
    fprintf(stderr,
            "  by class: profit %.3g (tol %.1g), argmax-derived %.3g (tol %.1g)\n",
            worst_max_rel, tol.maximum, worst_argmax_rel, tol.argmax);
  }
  if (exact_mode) {
    fprintf(stderr,
            "If this change was intended, regenerate with `make golden` and say "
            "so in the commit message.\n"
            "If you are on a platform other than the one that generated the "
            "file, do NOT regenerate -- use --cross-platform.\n");
  }
  return 1;
}

} // namespace

// Usage:
//   test_golden                   bit-exact. The default, and the mode that
//                                 matters on the platform that generated the file.
//   test_golden --cross-platform  per-field tolerances, for any OTHER platform.
//   test_golden --generate        overwrite the golden file (see the warning above).
//
// ---------------------------------------------------------------------------
// Why a cross-platform mode is needed, and why it is per-field
// ---------------------------------------------------------------------------
//
// The golden file is a *drift* detector: it pins what this implementation
// produces so that a refactor which perturbs it fails loudly. That needs
// bit-exactness, and gets it, on the platform that generated the file
// (macOS/arm64).
//
// It cannot be bit-exact anywhere else, and never could have been. libm's
// exp/pow are not bit-reproducible between glibc on x86-64 and Apple's libm on
// arm64, and FMA contraction differs too. Regenerating the file on a second
// platform is NOT the fix -- it just moves the failure to the first one.
//
// The interesting part is the SIZE of the disagreement, because it is not one
// number. Read off this program's own summary line on Linux CI, the nine reported
// fields split into two classes five orders of magnitude apart, and gcc and clang
// agree exactly:
//
//                                   gcc        clang
//     profit                        2.14e-09   2.14e-09
//     psi_stem, collar, ci, assim,
//     transpiration, gc, e_up,
//     uptake                        1.4e-04    1.4e-04
//
// ⚠️ NOTHING ASSERTS THAT BLOCK, so it rots in silence -- it has been wrong three
// times. **If you need a magnitude, read the summary line, not this comment and not
// the FAIL lines**: the failure list is truncated at 20 and biased toward whichever
// rows come first.
//
// That split is structural, not luck. `find_root_collar_psi` maximises profit
// over the collar potential, and the maximum is FLAT: curvature k measured
// directly at the two worst operating points gives k = 1.0 and 0.9 in
// profit ~ p* - k(psi_stem - x*)^2. For a flat maximum, an error dp in the
// profit VALUE displaces the ARGMAX by sqrt(dp/k), which is why a
// well-conditioned profit sits three orders below an ill-conditioned argmax.
//
// The relation was checked pointwise, not as a global identity: at the two worst
// points the residuals imply k = 1.01 and 1.70, against the 1.0 and 0.9 measured
// from the curvature directly. Do NOT expect sqrt(worst profit) to equal the
// worst argmax difference -- the two maxima occur at different operating points.
//
// Two consequences worth keeping in mind beyond this file:
//
//   * profit is the only reported field that is well-conditioned across
//     platforms, because it is the maximum itself rather than its location.
//   * 8 of the 9 fields are pinned to 17 digits but only *determined* to about
//     GSS_tol_abs (1e-3). Bit-exactness on one platform still works as a drift
//     detector -- any code change moves the arbitrary choice within that window,
//     and that shows up -- but it is reproducibility, not determinacy.
//
// So the tolerances below are per class, with the measured worst case and
// headroom noted at their definition. Both are far below the ~1e-4-and-up scale
// at which a real behavioural change shows (PLAN.md item 1), except that for the
// argmax fields the noise floor IS that scale -- which is exactly why the mode
// reports the worst value per class on every run rather than only on failure.
int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--generate") == 0) {
    return generate();
  }
  Tolerance tol = kExact;
  if (argc > 1) {
    if (std::strcmp(argv[1], "--cross-platform") == 0) {
      tol = kCrossPlatform;
    } else {
      fprintf(stderr,
              "unknown argument '%s'\n"
              "usage: test_golden [--cross-platform | --generate]\n",
              argv[1]);
      return 2;
    }
  }
  // Two independent checks over the same grid: the nine reported numbers against
  // the golden file, and the operating-point classification against the counts
  // recorded in the developer guide. Both run even if the first fails -- a
  // classification change and a numeric change have different causes, and seeing
  // only one of them is how a mixed diff gets misread. The grid is solved twice,
  // which costs ~4 ms.
  const int golden_status = compare(tol);
  const int kind_status = check_operating_kinds(run_grid());
  return (golden_status != 0 || kind_status != 0) ? 1 : 0;
}
