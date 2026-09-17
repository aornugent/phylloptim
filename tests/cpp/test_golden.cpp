// Golden-file regression test over a grid of operating points.
//
//   make -C tests/cpp golden        # regenerate golden/operating_points.tsv
//   make -C tests/cpp && ./test_golden
//
// Why this exists. A behaviour-preserving refactor can only be checked against
// behaviour that was pinned first. Cross-checking against plant's compiled build
// is the real validation, and it is DONE: the swap was bit-identical, including
// 78 of 78 SCM nodes. This file is the cheaper, always-available version -- it
// freezes what THIS implementation produces so that a refactor which changes any
// of it fails loudly.
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
// -- it WAS required, because the shutdown-state leak made a
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
  int interior, boundary_soil, boundary_crit, shutdown;
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
        {Kind::BoundarySoil, e.boundary_soil},
        {Kind::BoundaryCrit, e.boundary_crit},
        // The dry end's other bound. At these traits root_psi_crit is 5.87 MPa
        // and the continuity root never reaches it, so the arm is zero over the
        // whole grid -- asserted rather than left out, because a kind absent from
        // this list is a kind the golden file stops covering.
        {Kind::BoundaryRootCrit, 0},
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
      printf("  T = %4.1f C  %3d interior, %2d at a bound (%2d soil, %2d crit), "
             "%2d shutdown\n",
             e.leaf_temp, e.interior, e.boundary_soil + e.boundary_crit,
             e.boundary_soil, e.boundary_crit, e.shutdown);
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

// One row's identity, for a diagnostic.
//
// ⚠️ ROTATING BUFFERS, AND THE COUNT IS LOAD-BEARING. The summary prints four of
// these in ONE printf, and all four arguments are evaluated before it runs -- so
// with fewer buffers than call sites the later calls overwrite the earlier ones and
// the summary reports the same row two or four times. That is exactly what a single
// buffer did. Keep kSlots >= the largest number of calls in any one printf.
const char *where(const Row &r, const char *field) {
  constexpr int kSlots = 4;
  static char buf[kSlots][160];
  static int slot = 0;
  slot = (slot + 1) % kSlots;
  if (field[0] == '\0') {
    // No field in this class differed at all, so there is no row to name. On the
    // platform that generated the file that is every run, and printing a
    // zero-filled Row there reads like a real operating point at psi_soil = 0.
    std::snprintf(buf[slot], sizeof buf[slot], "%s", "no difference");
  } else {
    std::snprintf(buf[slot], sizeof buf[slot],
                  "psi_soil=%g ppfd=%g vpd=%g T=%g layers=%d %s", r.psi_soil,
                  r.ppfd, r.vpd, r.leaf_temp, r.layers, field);
  }
  return buf[slot];
}

// Exact equality, with NaN treated as equal to NaN -- some grid points shut down
// and legitimately produce the NA sentinel.
bool same(double got, double want) {
  if (std::isnan(got) && std::isnan(want)) {
    return true;
  }
  return got == want;
}

// Relative difference, used both to decide the tolerant comparison and to report
// how far apart the two builds actually are. Returns 0 for the NaN/NaN and
// exactly-equal cases, and infinity when one side is NaN and the other is not.
// Below this, a value in psi_stem_optima.tsv is zero in every sense that
// matters and a relative comparison against it is meaningless. Justified by a
// measured empty band -- see the use site in compare_optima.
constexpr double kOptNegligible = 1e-12;

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
// headroom: 1e-4 is the scale at which a real behavioural change shows, so a
// profit tolerance approaching it would gate nothing. If a future
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
  // ⚠️ PER CLASS, AND NAMED. The two class magnitudes are what the developer
  // guide quotes, and they have been wrong four times -- most recently the pair
  // that stood through a solver merge which moved one of them 87x. A magnitude
  // with no row attached can only be watched drift; one that names its row can be
  // reproduced, which is the difference between a number and a diagnosis.
  Row worst_max_row{}, worst_argmax_row{};
  const char *worst_max_field = "";
  const char *worst_argmax_field = "";
  // ⚠️ AND THE ABSOLUTE DIFFERENCE, because the relative one alone has misled this
  // project's own documentation for years. A relative metric divides by the value,
  // so wherever a reported quantity passes through zero it is amplified by the
  // cancellation -- and BOTH class figures turn out to land on such a row. profit's
  // worst sits at the grid's MINIMUM |profit| (1.05e-03 against a median of 2.47,
  // 2348x), and the argmax class's at an assimilation of 4.05e-07 against a median
  // of 1.44. So "profit is well-conditioned and the argmax is sqrt-amplified", which
  // this file's guide has said for a long time, is a statement about the metric at
  // two singular rows rather than about the model at an operating point.
  //
  // Reporting both is the fix. It loosens nothing -- the tolerance is still the
  // relative one -- and it makes a cancellation row distinguishable from real drift,
  // which is what watching this number is supposed to be for.
  double worst_max_abs = 0.0, worst_argmax_abs = 0.0;
  Row worst_max_abs_row{}, worst_argmax_abs_row{};
  const char *worst_max_abs_field = "";
  const char *worst_argmax_abs_field = "";
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
        if (is_max) {
          worst_max_row = r;
          worst_max_field = fd.name;
        } else {
          worst_argmax_row = r;
          worst_argmax_field = fd.name;
        }
      }
      const double ad = (std::isnan(fd.got) || std::isnan(fd.want))
                            ? 0.0
                            : std::fabs(fd.got - fd.want);
      double &worst_abs_in_class = is_max ? worst_max_abs : worst_argmax_abs;
      if (ad > worst_abs_in_class) {
        worst_abs_in_class = ad;
        if (is_max) {
          worst_max_abs_row = r;
          worst_max_abs_field = fd.name;
        } else {
          worst_argmax_abs_row = r;
          worst_argmax_abs_field = fd.name;
        }
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
             "    profit  (the maximum)      rel %.3g  (tol %.1g)  at %s\n"
             "                               abs %.3g              at %s\n"
             "    others  (from the argmax)  rel %.3g  (tol %.1g)  at %s\n"
             "                               abs %.3g              at %s\n"
             "  ⚠️ a worst RELATIVE row with a near-zero value is the metric "
             "cancelling, not drift -- compare the two lines in each class.\n",
             rows.size(), inexact, rows.size() * 9,
             worst_max_rel, tol.maximum, where(worst_max_row, worst_max_field),
             worst_max_abs, where(worst_max_abs_row, worst_max_abs_field),
             worst_argmax_rel, tol.argmax,
             where(worst_argmax_row, worst_argmax_field),
             worst_argmax_abs,
             where(worst_argmax_abs_row, worst_argmax_abs_field));
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
            "  by class: profit %.3g (tol %.1g) at %s, argmax-derived %.3g "
            "(tol %.1g) at %s\n",
            worst_max_rel, tol.maximum, where(worst_max_row, worst_max_field),
            worst_argmax_rel, tol.argmax,
            where(worst_argmax_row, worst_argmax_field));
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

// ===========================================================================
// The psi_stem optimiser fixture (golden/psi_stem_optima.tsv)
// ---------------------------------------------------------------------------
// WHY A SECOND FILE. `operating_points.tsv` above solves with
// `find_root_collar_psi()` and nothing else, so it is blind BY CONSTRUCTION to
// the three `optimise_psi_stem_*` entry points -- and so is the R suite, whose
// 112 test_that blocks call none of them. Every existing C++ test that touches
// them asserts a RELATION (a first-order condition, a monotonicity, a thrown
// message) at spot-checked settings; not one records a number. So today it is
// possible to detect that a refactor broke a relation, and impossible to say
// what it changed.
//
// ⚠️ THIS FILE RECORDS TODAY'S BEHAVIOUR, INCLUDING ITS KNOWN DEFECTS, ON
// PURPOSE. Do not read a row as a statement about what the model SHOULD do.
// Specifically recorded here as-is:
//   * two different degenerate conventions -- `_TF` and `_CF77`
//     evaluate the real objective at the closed point, `_ProfitMax` zeroes seven
//     fields instead;
//   * `hydraulic_cost_` written in carbon units by `_TF` and `_CF77`
//     and left STALE by `_ProfitMax`, which does not write it;
//   * `CF77_lambda_` as a caller input that `_ProfitMax` overwrites.
//
// ⚠️ ADDING A SOLVER IS AN ADDITION IN THE "fresh" PASS AND NOT IN THE "reuse"
// ONE, so check the two separately when regenerating. The reuse pass drives one
// Leaf through a fixed SEQUENCE, so inserting a solver changes what every solver
// after it inherits. Adding Cowan-Farquhar moved 232 reuse rows and zero fresh
// rows, and every one of those 232 moved in the `lambda` column alone: it writes
// `CF77_lambda_`, so the solvers that follow it and do not write one report its value
// instead of the previous solver's. That is the leak this pass exists to catch,
// not a defect in the new model -- but it means "regenerated and the diff is
// pure addition" is the wrong check here. The right one is: the fresh pass is
// additive, and every reuse-pass difference is attributable to a named cause.
//
// TWO PASSES, and the second is the one the file exists for. Pass "fresh"
// constructs a Leaf per row, which makes it reproducible and blind to stale
// state -- the same trade the grid above makes. Pass "reuse" drives ONE Leaf
// through a fixed sequence of rows, poisoning every stale-prone output first, so
// a value that survives is provably not written by the path under test. Four of
// the divergences above are stale-state bugs that pass "fresh" cannot see.
//
// The two axes that make a different MODEL rather than a different driver:
//   solver    which objective, and which variable it optimises. `collar`
//             optimises the COLLAR potential with the root path in series;
//             the other three optimise psi_stem with upstream pinned at
//             psi_soil, ignoring that path -- which is what their own refusal
//             message means by "non-root-based". These are different problems,
//             and the gap does not close as the root resistance goes to zero --
//             the collar loses its freedom there and reports `determined`
//             instead of optimising. Neither arm is a check on the other.
//   topology  single-potential supply, or a multi-layer profile at 1 or 3
//             layers. The three psi_stem solvers refuse 3 layers; the refusal
//             is RECORDED rather than skipped, so a change that lifts it shows
//             up as a row that stopped throwing.

const char *kOptimaPath = "golden/psi_stem_optima.tsv";

// A value no solve should ever produce, written into every stale-prone output
// before each solve in the reuse pass. If it comes back, the path under test did
// not write that field.
const double kPoison = -12345.0;

// The prescribed marginal water cost, for the two curves that consume one.
// Fixed rather than taken from a prior ProfitMax solve: an ordering dependence is
// exactly what this file is for detecting, not something to bake into it.
//
// ⚠️ ITS SCALE IS SET BY THE LEAF, not chosen freely: lambda is umol C per kg of
// water, and at these drivers the leaf's own marginal cost of water runs 9e4 to
// 3e5. A value far outside that band pins the optimum against a bound, where the
// row records the bracket rather than the model.
const double kLambdaCF77 = 1.5e5;

// The TF24_floor curve's price of water as transpiration goes to zero. Same units,
// same scale argument, and a separate constant rather than a reuse of the one
// above: the two are independent inputs and sharing them would make a stale write
// to either invisible in the other's rows. The curve's hydraulic half is TF24's
// cost at TF24's own traits, so there is nothing else to set.
const double kLambdaFloorO = 1.5e5;

// ⚠️ APPEND ONLY. The reuse pass drives ONE Leaf through kSolvers in order, so
// inserting a solver changes what every solver after it inherits -- see the header.
// Appending leaves every existing row byte-identical, which is the check when this
// file is regenerated.
enum class Solver { TF, ProfitMax, CF77, Collar, JS22, CMax, SOX, JW26, TF24_floor };

// ⚠️ EVERY PER-SOLVER ARRAY IS SIZED FROM THIS, never from a literal. Those arrays
// are indexed by the enum VALUE, so a hardcoded bound is an out-of-bounds write the
// moment a solver is appended -- which is exactly the bug bench_gradient carried
// through three trait-count changes before ASan named it. Keep it last-member + 1.
constexpr int kNSolvers = static_cast<int>(Solver::TF24_floor) + 1;
enum class Topology { Single, Multi1, Multi3 };

const char *solver_name(Solver s) {
  switch (s) {
    case Solver::TF:            return "TF";
    case Solver::ProfitMax:     return "ProfitMax";
    case Solver::CF77: return "CF77";
    case Solver::Collar:        return "collar";
    case Solver::JS22:          return "JS22";
    case Solver::CMax:          return "CMax";
    case Solver::SOX:           return "SOX";
    case Solver::JW26:          return "JW26";
    case Solver::TF24_floor:       return "TF24_floor";
  }
  return "unknown";
}

const char *topology_name(Topology t) {
  switch (t) {
    case Topology::Single: return "single";
    case Topology::Multi1: return "multi1";
    case Topology::Multi3: return "multi3";
  }
  return "unknown";
}

int topology_layers(Topology t) {
  switch (t) {
    case Topology::Single: return 1;
    case Topology::Multi1: return 1;
    case Topology::Multi3: return 3;
  }
  return 0;
}

struct OptRow {
  // inputs
  Solver solver;
  Topology topology;
  double psi_soil, ppfd, leaf_temp;
  bool energy_balance, thermal_cost;
  bool reuse_pass;
  // what happened
  const char *status; // "ok" or "threw"
  const char *kind;   // operating_point_kind_name
  // outputs
  double opt_psi_stem, opt_root_psi, profit, ci, assim, transpiration, gc,
      hydraulic_cost, lambda, lambda_emergent, tleaf, carbon_gain,
      hydraulic_cost_norm, thermal_cost_out;
};

// The single-potential path's series resistance. Positive and finite is a
// requirement of that path, so there is no r = 0 arm to compare against.
const double kSeriesResistance = 1.0e3;

void configure_supply(phylloptim::Leaf &l, Topology t, double psi_soil,
                      double ppfd, double leaf_temp) {
  const int layers = topology_layers(t);
  if (t == Topology::Single) {
    l.set_supply_single(0.0);
    phylloptim::RootNetwork rn;
    rn.r_R_V_sum = std::vector<double>{kSeriesResistance};
    rn.r_R_H_min = std::vector<double>{0.0};
    rn.r_R_V = std::vector<double>{kSeriesResistance};
    rn.c_r_V = std::vector<double>{0.0};
    rn.c_r_H = std::vector<double>{0.0};
    l.set_physiology(rn, ppfd, {psi_soil}, {1.0}, kKs * kTheta / kH, 2.0, kCa,
                     leaf_temp, kO2, kPatm);
    return;
  }
  l.set_supply_multilayer();
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers / kAreaLeaf;
  }
  l.set_physiology(fixture::root_network(root, depth), ppfd, ps, depth,
                   kKs * kTheta / kH, 2.0, kCa, leaf_temp, kO2, kPatm);
}

void poison(phylloptim::Leaf &l) {
  l.ci_ = kPoison;
  l.assim_colimited_ = kPoison;
  l.transpiration_ = kPoison;
  l.stom_cond_CO2_ = kPoison;
  l.hydraulic_cost_ = kPoison;
  l.opt_psi_stem_ = kPoison;
  l.profit_ = kPoison;
  l.carbon_gain_ = kPoison;
  l.hydraulic_cost_norm_ = kPoison;
  l.thermal_cost_ = kPoison;
}

void read_outputs(const phylloptim::Leaf &l, OptRow &r) {
  r.kind = phylloptim::Leaf::operating_point_kind_name(l.operating_point_kind());
  r.opt_psi_stem = l.opt_psi_stem_;
  r.opt_root_psi = l.opt_root_psi_;
  r.profit = l.profit_;
  r.ci = l.ci_;
  r.assim = l.assim_colimited_;
  r.transpiration = l.transpiration_;
  r.gc = l.stom_cond_CO2_;
  r.hydraulic_cost = l.hydraulic_cost_;
  r.lambda = l.CF77_lambda_;
  r.lambda_emergent = l.lambda_emergent();
  r.tleaf = l.Tleaf_;
  r.carbon_gain = l.carbon_gain_;
  r.hydraulic_cost_norm = l.hydraulic_cost_norm_;
  r.thermal_cost_out = l.thermal_cost_;
}

void blank_outputs(OptRow &r) {
  const double n = std::numeric_limits<double>::quiet_NaN();
  r.kind = "-";
  r.opt_psi_stem = r.opt_root_psi = r.profit = r.ci = r.assim = n;
  r.transpiration = r.gc = r.hydraulic_cost = r.lambda = n;
  r.lambda_emergent = r.tleaf = n;
  r.carbon_gain = r.hydraulic_cost_norm = r.thermal_cost_out = n;
}

// Run one solver on an already-configured leaf. Separated from configuration so
// the reuse pass can drive one Leaf through many rows.
void dispatch(phylloptim::Leaf &l, Solver s) {
  switch (s) {
    case Solver::TF:        l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24>();        break;
    case Solver::ProfitMax: l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::ProfitMax>(); break;
    case Solver::CF77:
                            l.CF77_lambda_ = kLambdaCF77;
                            l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CF77>(); break;
    case Solver::Collar:    l.find_root_collar_psi();        break;
    // These three need nothing set first: their parameters are TRAITS with
    // defaults, where CF77_lambda_ above is an input with none.
    case Solver::JS22:      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::JS22>();      break;
    case Solver::CMax:      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::CMax>();      break;
    case Solver::SOX:       l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::SOX>();       break;
    case Solver::JW26:      l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::JW26>();      break;
    // ⚠️ ONE PARAMETER, AND IT IS AN INPUT WITH NO DEFAULT, like CF77_lambda_
    // above: this curve refuses to solve without it. The hydraulic half reads
    // TF24's own traits, so these rows differ from the TF row above by the price
    // alone.
    case Solver::TF24_floor:
                            l.TF24_floor_lambda_o = kLambdaFloorO;
                            l.optimise_psi_stem_single<phylloptim::Leaf::CostCurve::TF24_floor>(); break;
  }
}

OptRow solve_one(Solver s, Topology t, double psi_soil, double ppfd,
                 double leaf_temp, bool eb, bool tc) {
  OptRow r{s, t, psi_soil, ppfd, leaf_temp, eb, tc, false, "ok", "-",
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = eb;
  l.use_thermal_cost_ = tc;
  try {
    configure_supply(l, t, psi_soil, ppfd, leaf_temp);
    dispatch(l, s);
  } catch (...) {
    r.status = "threw";
    blank_outputs(r);
    return r;
  }
  read_outputs(l, r);
  return r;
}

const double kOptPsiSoils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
const double kOptTemps[] = {25.0, 40.0, 50.0};
const double kOptPPFDs[] = {0.0, 1500.0};
const Solver kSolvers[] = {Solver::TF, Solver::ProfitMax,
                           Solver::CF77, Solver::Collar,
                           Solver::JS22, Solver::CMax, Solver::SOX,
                           Solver::JW26, Solver::TF24_floor};
static_assert(sizeof(kSolvers) / sizeof(kSolvers[0]) == kNSolvers,
              "kSolvers must list every Solver exactly once");
const Topology kTopologies[] = {Topology::Single, Topology::Multi1,
                                Topology::Multi3};

// Pass "reuse" drives ONE Leaf per (topology, gates) through a reduced driver
// set, in a FIXED order, poisoning the outputs before each solve. The driver set
// is reduced because staleness does not need the full cross -- it needs a wet
// row, a dry row, a shut row and a hot row, which is what these three potentials
// and two temperatures give.
const double kReusePsiSoils[] = {0.5, 3.0, 6.0};
const double kReuseTemps[] = {25.0, 50.0};

std::vector<OptRow> run_optima_grid() {
  std::vector<OptRow> rows;

  // --- pass "fresh": a new Leaf per row -----------------------------------
  for (Topology t : kTopologies)
    for (bool eb : {false, true})
      for (bool tc : {false, true})
        for (double ps : kOptPsiSoils)
          for (double T : kOptTemps)
            for (double ppfd : kOptPPFDs)
              for (Solver s : kSolvers) {
                rows.push_back(solve_one(s, t, ps, ppfd, T, eb, tc));
              }

  // --- pass "reuse": one Leaf, poisoned between solves ---------------------
  for (Topology t : kTopologies)
    for (bool eb : {false, true})
      for (bool tc : {false, true}) {
        phylloptim::Leaf l;
        l.setup_transpiration(100);
        l.setup_root_vulnerability(100);
        l.use_energy_balance_ = eb;
        l.use_thermal_cost_ = tc;
        for (double ps : kReusePsiSoils)
          for (double T : kReuseTemps)
            for (double ppfd : kOptPPFDs)
              for (Solver s : kSolvers) {
                OptRow r{s,  t,  ps, ppfd, T, eb, tc, true, "ok", "-",
                         0,  0,  0,  0,    0, 0,  0,  0,    0,    0, 0, 0, 0,
                         0};
                try {
                  configure_supply(l, t, ps, ppfd, T);
                  poison(l);
                  dispatch(l, s);
                } catch (...) {
                  r.status = "threw";
                  blank_outputs(r);
                  rows.push_back(r);
                  continue;
                }
                read_outputs(l, r);
                rows.push_back(r);
              }
      }
  return rows;
}

const char *kOptimaHeader =
    "pass\tsolver\ttopology\tpsi_soil\tppfd\tleaf_temp\teb\ttc\tstatus\tkind\t"
    "opt_psi_stem\topt_root_psi\tprofit\tci\tassim\ttranspiration\tgc\t"
    "hydraulic_cost\tlambda\tlambda_emergent\ttleaf\tcarbon_gain\t"
    "hydraulic_cost_norm\tthermal_cost\n";

void write_opt_row(FILE *f, const OptRow &r) {
  fprintf(f, "%s\t%s\t%s\t%.17g\t%.17g\t%.17g\t%d\t%d\t%s\t%s",
          r.reuse_pass ? "reuse" : "fresh", solver_name(r.solver),
          topology_name(r.topology), r.psi_soil, r.ppfd, r.leaf_temp,
          r.energy_balance ? 1 : 0, r.thermal_cost ? 1 : 0, r.status, r.kind);
  for (double v : {r.opt_psi_stem, r.opt_root_psi, r.profit, r.ci, r.assim,
                   r.transpiration, r.gc, r.hydraulic_cost, r.lambda,
                   r.lambda_emergent, r.tleaf, r.carbon_gain,
                   r.hydraulic_cost_norm, r.thermal_cost_out}) {
    fprintf(f, "\t%.17g", v);
  }
  fprintf(f, "\n");
}

int generate_optima() {
  const std::vector<OptRow> rows = run_optima_grid();
  FILE *f = fopen(kOptimaPath, "w");
  if (f == nullptr) {
    fprintf(stderr, "cannot write %s (run from tests/cpp/)\n", kOptimaPath);
    return 1;
  }
  fputs(kOptimaHeader, f);
  for (const OptRow &r : rows) {
    write_opt_row(f, r);
  }
  fclose(f);
  printf("wrote %zu optimiser rows to %s\n", rows.size(), kOptimaPath);
  return 0;
}

// The thirteen numeric outputs, in the order write_opt_row emits them.
const char *kOptFieldNames[] = {
    "opt_psi_stem", "opt_root_psi", "profit", "ci", "assim", "transpiration",
    "gc", "hydraulic_cost", "lambda", "lambda_emergent", "tleaf",
    "carbon_gain", "hydraulic_cost_norm", "thermal_cost"};

void opt_row_values(const OptRow &r, double *out) {
  out[0] = r.opt_psi_stem;   out[1] = r.opt_root_psi;
  out[2] = r.profit;         out[3] = r.ci;
  out[4] = r.assim;          out[5] = r.transpiration;
  out[6] = r.gc;             out[7] = r.hydraulic_cost;
  out[8] = r.lambda;         out[9] = r.lambda_emergent;
  out[10] = r.tleaf;         out[11] = r.carbon_gain;
  out[12] = r.hydraulic_cost_norm;
  out[13] = r.thermal_cost_out;
}

// ⚠️ THE CROSS-PLATFORM TOLERANCES HERE ARE INHERITED, NOT MEASURED. `kExact` is
// right on the generating platform and needs no justification. The `--cross-platform`
// arm reuses the collar grid's per-field pair, which was measured for THAT grid's
// nine fields; this file has thirteen, three of which (`lambda`, `carbon_gain`,
// `thermal_cost`) have no counterpart there and so no measured magnitude at all.
// The pair is therefore a starting guess. Two things follow:
//
//   * the two NON-numeric columns -- `status` and `kind` -- are compared exactly
//     in every mode, because a refusal that stopped refusing or a branch that
//     changed is never a rounding difference, and those are the changes this file
//     most needs to catch;
//   * the worst observed relative difference is printed even when the comparison
//     PASSES, so the real magnitude can be read off the first non-macOS CI run
//     rather than bisected out of a red one.
//
// Replace the inherited pair with what that run reports, and say in the commit
// which platform it came from.
int compare_optima(Tolerance tol) {
  FILE *f = fopen(kOptimaPath, "r");
  if (f == nullptr) {
    fprintf(stderr,
            "FAIL: %s is missing. Generate it with `make psi-stem-golden` and "
            "commit it.\n",
            kOptimaPath);
    return 1;
  }
  char line[8192];
  if (fgets(line, sizeof line, f) == nullptr) {
    fprintf(stderr, "FAIL: %s is empty\n", kOptimaPath);
    fclose(f);
    return 1;
  }

  const std::vector<OptRow> rows = run_optima_grid();
  const bool exact_mode = tol.maximum < 0.0;
  int failures = 0;
  size_t negligible = 0;
  double worst_rel = 0.0;
  std::string worst_desc;

  size_t i = 0;
  for (; i < rows.size(); ++i) {
    if (fgets(line, sizeof line, f) == nullptr) {
      fprintf(stderr, "FAIL: %s has only %zu rows, grid has %zu\n", kOptimaPath,
              i, rows.size());
      ++failures;
      break;
    }
    char pass[16], solver[24], topo[16], status[16], kind[32];
    double psi_soil, ppfd, leaf_temp;
    int eb, tc;
    double want[14];
    const int n = sscanf(line,
                         "%15s %23s %15s %lf %lf %lf %d %d %15s %31s"
                         " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                         pass, solver, topo, &psi_soil, &ppfd, &leaf_temp, &eb,
                         &tc, status, kind, &want[0], &want[1], &want[2],
                         &want[3], &want[4], &want[5], &want[6], &want[7],
                         &want[8], &want[9], &want[10], &want[11], &want[12],
                         &want[13]);
    if (n != 24) {
      fprintf(stderr, "FAIL: %s line %zu has %d fields, expected 24\n",
              kOptimaPath, i + 2, n);
      ++failures;
      continue;
    }

    const OptRow &r = rows[i];
    char tag[192];
    snprintf(tag, sizeof tag, "%s/%s/%s psi_soil=%g ppfd=%g T=%g eb=%d tc=%d",
             r.reuse_pass ? "reuse" : "fresh", solver_name(r.solver),
             topology_name(r.topology), r.psi_soil, r.ppfd, r.leaf_temp,
             r.energy_balance ? 1 : 0, r.thermal_cost ? 1 : 0);

    // The three non-numeric columns are compared exactly in every mode: a
    // status or classification change is never a rounding difference.
    if (std::strcmp(status, r.status) != 0) {
      fprintf(stderr, "FAIL %s: status %s, expected %s\n", tag, r.status,
              status);
      ++failures;
    }
    if (std::strcmp(kind, r.kind) != 0) {
      fprintf(stderr, "FAIL %s: kind %s, expected %s\n", tag, r.kind, kind);
      ++failures;
    }

    double got[14];
    opt_row_values(r, got);
    for (int j = 0; j < 14; ++j) {
      if (same(got[j], want[j])) {
        continue;
      }
      // ⚠️ A RELATIVE TOLERANCE IS UNDEFINED WHERE THE TRUTH IS ZERO, and two
      // columns here are analytically zero on some rows: `hydraulic_cost_norm`
      // at a boundary-soil ProfitMax optimum, and `profit` likewise. One platform
      // then computes rounding noise where another computes exactly 0.0, and
      // rel_diff correctly reports 1 for a difference of 1e-16.
      //
      // The threshold sits in an EMPTY BAND, which is what makes it safe rather
      // than convenient. Measured over the whole file: the noise cluster tops
      // out at 1.13e-16 (76 values) and the smallest genuine magnitude in either
      // column is >= 1e-4 (796 values) -- NOTHING lies between 1e-12 and 1e-4.
      // So this cannot mask a real difference; it can only stop a zero being
      // compared as though it had a scale. Re-measure the band before moving it.
      if (!exact_mode && std::fabs(got[j]) < kOptNegligible &&
          std::fabs(want[j]) < kOptNegligible) {
        ++negligible;
        continue;
      }
      const double rel = rel_diff(got[j], want[j]);
      const double limit =
          is_maximum_field(kOptFieldNames[j]) ? tol.maximum : tol.argmax;
      if (!exact_mode && rel <= limit) {
        if (rel > worst_rel) {
          worst_rel = rel;
          worst_desc = std::string(tag) + " " + kOptFieldNames[j];
        }
        continue;
      }
      if (rel > worst_rel) {
        worst_rel = rel;
        worst_desc = std::string(tag) + " " + kOptFieldNames[j];
      }
      if (failures < 20) {
        fprintf(stderr, "FAIL %s: %s got %.17g, expected %.17g (rel %.3g)\n",
                tag, kOptFieldNames[j], got[j], want[j], rel);
      }
      ++failures;
    }
  }
  if (failures == 0 && fgets(line, sizeof line, f) != nullptr) {
    fprintf(stderr, "FAIL: %s has more rows than the grid (%zu)\n", kOptimaPath,
            rows.size());
    ++failures;
  }
  fclose(f);

  if (failures == 0) {
    printf("psi_stem optima: %zu rows, all %s\n", rows.size(),
           exact_mode ? "bit-identical"
                      : "within cross-platform tolerance");
    // ⚠️ PRINTED ON SUCCESS TOO, AND ON PURPOSE. The tolerances this file is
    // compared with off the generating platform are INHERITED from the collar
    // grid, not measured here -- see the note above `compare_optima`. Until
    // someone reads a real number off a non-macOS/arm64 run, this line is the
    // only place that number exists, and a summary that only appears on failure
    // cannot supply it. Reading it off a green run is the cheap way; bisecting a
    // red one is not.
    if (!exact_mode && !worst_desc.empty()) {
      printf("  worst relative difference %.3g at %s (tolerances %.1g maximum,"
             " %.1g argmax -- INHERITED, not measured for this file)\n",
             worst_rel, worst_desc.c_str(), tol.maximum, tol.argmax);
    }
    // Never silent: an exempted cell is a cell nothing checked.
    if (!exact_mode && negligible > 0) {
      printf("  %zu values compared as zero (both sides below %.1g)\n",
             negligible, kOptNegligible);
    }
    return 0;
  }
  fprintf(stderr, "\npsi_stem optima: %d mismatches over %zu rows.\n",
          failures, rows.size());
  if (!worst_desc.empty()) {
    fprintf(stderr, "  worst relative difference %.3g at %s\n", worst_rel,
            worst_desc.c_str());
  }
  if (exact_mode) {
    fprintf(stderr,
            "If this change was intended, regenerate with `make psi-stem-golden` "
            "and say what moved in the commit message.\n");
  }
  return 1;
}

// A standing summary of what the fixture contains, printed on every run. It is
// not an assertion -- the counts are recorded in the file itself -- but a
// refusal count that moves is the fastest signal that a supply restriction was
// lifted or added, and it costs nothing to look at.
void report_optima_shape(const std::vector<OptRow> &rows) {
  int threw = 0;
  int by_solver_threw[kNSolvers] = {0};
  int poisoned = 0;
  for (const OptRow &r : rows) {
    if (std::strcmp(r.status, "threw") == 0) {
      ++threw;
      by_solver_threw[static_cast<int>(r.solver)]++;
      continue;
    }
    double v[14];
    opt_row_values(r, v);
    for (int j = 0; j < 14; ++j) {
      if (v[j] == kPoison) {
        ++poisoned;
        break;
      }
    }
  }
  printf("  %zu rows: %d refused", rows.size(), threw);
  for (int s = 0; s < kNSolvers; ++s) {
    if (by_solver_threw[s] != 0) {
      printf(" (%s %d)", solver_name(static_cast<Solver>(s)),
             by_solver_threw[s]);
    }
  }
  printf(", %d rows carry an unwritten field into their output\n", poisoned);
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
// at which a real behavioural change shows, except that for the
// argmax fields the noise floor IS that scale -- which is exactly why the mode
// reports the worst value per class on every run rather than only on failure.
int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--generate") == 0) {
    return generate();
  }
  if (argc > 1 && std::strcmp(argv[1], "--generate-optima") == 0) {
    return generate_optima();
  }
  Tolerance tol = kExact;
  if (argc > 1) {
    if (std::strcmp(argv[1], "--cross-platform") == 0) {
      tol = kCrossPlatform;
    } else {
      fprintf(stderr,
              "unknown argument '%s'\n"
              "usage: test_golden [--cross-platform | --generate"
              " | --generate-optima]\n",
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
  // The psi_stem optimiser fixture is a third INDEPENDENT check, and it runs
  // even when the two above fail: the collar grid and the optimiser grid have
  // disjoint causes, so seeing only one of them is how a mixed diff gets
  // misread. See the block above kOptimaPath for what this file is and is not.
  const int optima_status = compare_optima(tol);
  report_optima_shape(run_optima_grid());
  return (golden_status != 0 || kind_status != 0 || optima_status != 0) ? 1 : 0;
}
