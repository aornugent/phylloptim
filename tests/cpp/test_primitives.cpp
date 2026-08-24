// Golden-file regression test over the model's PRIMITIVES, in call-tree order.
//
//   make -C tests/cpp                     builds and runs it
//   make -C tests/cpp primitives-golden   regenerate golden/primitives.tsv
//   ./test_primitives --cross-platform    per-tier tolerances, off macOS/arm64
//
// ---------------------------------------------------------------------------
// Why this exists, when golden/operating_points.tsv already does
// ---------------------------------------------------------------------------
//
// That file compares the SOLVE. This one compares the functions the solve calls,
// and the difference is localisation. The nested solvers amplify a perturbation up
// to whichever tolerance is loosest, so by the time a last-bit difference reaches a
// reported output it has been through a collar root-find and a ci iteration and no
// longer points anywhere. The operating-point file can say "3496 cells moved"; it
// cannot say which function moved them.
//
// This is issue #64's option 2, and it is filed there with the two occasions that
// made the case. Both were #92, the vulnerability-domain fix:
//
//   * The DEFECT was invisible to the operating-point file. The knot grid dropped
//     its last knot at the package defaults for years -- 99 knots ending at 6.8229
//     against a psi_max of 6.8918, one full step short -- and the golden file
//     recorded the consequences as the correct answer. Only comparing the realised
//     domain against vulnerability_psi_max finds that, which is a primitive-level
//     check.
//   * The FIX's blast radius could not be attributed from it either. Deciding
//     whether 3496 moved cells were a knot displacement (~1e-15) or the added knot
//     needed a standalone program over the primitives, written in a scratch
//     directory and thrown away. Twice.
//
// ---------------------------------------------------------------------------
// What the reference is, which is the question #64 treats as the blocker
// ---------------------------------------------------------------------------
//
// A checked-in table of primitive values, generated deliberately from this
// implementation and compared bit-exactly thereafter. Its authority comes from
// exactly where operating_points.tsv's does.
//
// It is NOT a comparison against plant, and it cannot be: plant has no independent
// copy of this model any more, it consumes these headers, so "compare against
// plant" would be comparing this file with itself. The two `compare_*.R` harnesses
// that did compare against plant are deleted rather than revived, for that reason.
// What they established -- bit-identical, including 78 of 78 SCM nodes -- is
// on record here and cannot be re-run.
//
// ⚠️ ONE THING THE ORIGINAL HARNESS NEEDED AND THIS DOES NOT: hex floats. It emitted
// every value twice, decimal and `%a`, because its counterpart read the file into R
// and R's decimal parser is not correctly rounded -- `as.numeric`, `scan` and
// `read.delim` all return a double one ULP off for about 18% of inputs, which
// manufactures exactly the disagreement it was built to explain. Both sides here are
// C++ and the reader is `strtod`, which IS correctly rounded, so `%.17g` round-trips
// exactly. Do not "restore" the hex column: it would be a second representation of
// the same number for no reader.
//
// ---------------------------------------------------------------------------
// Tiers, and why they are the whole design
// ---------------------------------------------------------------------------
//
// Rows are emitted in CALL-TREE order: pure arithmetic first, then things that call
// it. So the FIRST tier that disagrees is the answer, and the tiers below it are
// consequences rather than clues. That is the method the dead harness encoded and
// the reason #64 judged it worth reviving rather than deleting.
//
// The tolerant mode reports the worst deviation PER TIER on every run, not one
// number for the file. A cross-platform diff that is 1e-16 in tier 1 and 1e-4 in
// tier 5 is libm plus solver amplification and is expected; the same 1e-4 appearing
// in tier 1 is a real change in the arithmetic. One pooled figure cannot make that
// distinction, which is the thing operating_points.tsv cannot do.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

const char *kGoldenPath = "golden/primitives.tsv";

// The tiers, in call-tree order. Named here rather than at each emit site so the
// per-tier reporting and the file's own column agree by construction.
enum Tier {
  kArithmetic = 1,  // no state, no splines, no iteration: libm and FMA only
  kVulnerability,   // the Weibull curve and its pre-integrated grid
  kAssimilation,    // rational arithmetic, then a quadratic root
  kSpline,          // the first functions that READ an interpolator
  kIterative,       // primitives containing their own root-find
  kTierCount
};

const char *tier_name(int t) {
  switch (t) {
  case kArithmetic:
    return "arithmetic";
  case kVulnerability:
    return "vulnerability";
  case kAssimilation:
    return "assimilation";
  case kSpline:
    return "spline";
  case kIterative:
    return "iterative";
  default:
    return "?";
  }
}

struct Row {
  int tier;
  std::string fn;
  double a1, a2;
  double value;
};

// Identical to tests/cpp/test_golden.cpp's fixture. Spelled out rather than left to
// defaults for the same reason that file does: a default that moved would move this
// baseline without anything naming the change.
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0;
const double kAreaLeaf = 0.05;
const double kCa = 40.0, kO2 = 21.0, kTleaf = 25.0, kPatm = 101.3;

phylloptim::Leaf fixture_leaf() {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  // ONE single-layer configuration at psi_soil = 2 MPa. Single-layer on purpose:
  // the multi-layer supply path is a separate question and would confound this
  // one, which is about the leaf-side functions.
  std::vector<double> ps{2.0}, depth{1.0}, root{1.0 / kAreaLeaf};
  l.set_physiology(fixture::root_network(root, depth), /*PPFD*/ 900.0, ps, depth,
                   kKs * kTheta / kH, /*vpd*/ 2.0, kCa, kTleaf, kO2, kPatm);
  return l;
}

std::vector<Row> run_grid() {
  phylloptim::Leaf l = fixture_leaf();
  std::vector<Row> rows;
  const auto row = [&rows](int tier, const char *fn, double a1, double a2,
                           double v) { rows.push_back({tier, fn, a1, a2, v}); };

  // --- tier 1: pure arithmetic ----------------------------------------------
  //
  // If anything here disagrees, the cause is libm or FMA contraction and the
  // search is over.
  //
  // ⚠️ SWEPT OVER TEMPERATURE, AND NOT EVALUATED AT 25 C. Both curves are
  // normalised to their 25 C reference, so at 25 C they return that reference
  // unchanged and would compare equal no matter how broken the exp() underneath
  // was. This is the same blindness the operating-point grid had before it gained
  // a second temperature.
  const double eas[] = {58550.0, 62500.0, 36380.0};
  const double refs[] = {96.0, 157.44, 40.0};
  for (double t = 5.0; t <= 45.0 + 1e-12; t += 2.5) {
    for (double ea : eas) {
      row(kArithmetic, "arrh_curve", ea, t, l.arrh_curve(ea, refs[0], t));
    }
    for (double ref : refs) {
      row(kArithmetic, "peak_arrh_curve", ref, t,
          l.peak_arrh_curve(43900.0, ref, t, 200000.0, 640.0));
    }
    // R_d's Tjoelker declining-Q10 response (#41), which replaced a peaked
    // Arrhenius and is therefore NOT covered by the two curves above. Written out
    // rather than read off R_d_ so that this row is the response and not the
    // response composed with whatever set_physiology last cached.
    const double q10 = 3.09 - 0.0430 * (t + 25.0) / 2.0;
    row(kArithmetic, "rd_q10", t, 0.0, q10);
    row(kArithmetic, "rd_response", t, 0.0, std::pow(q10, (t - 25.0) / 10.0));
    // The PM energy balance, which is arithmetic in E and has no other primitive
    // of its own. Its slope is the thing the collar FOC differentiates (#84).
    double dT_dE = 0.0;
    const double E = 1.0e-5 * (t - 4.0) / 41.0;
    const double T_leaf = l.leaf_temp_from_E(E, &dT_dE);
    row(kArithmetic, "leaf_temp_from_E", E, 0.0, T_leaf);
    row(kArithmetic, "dleaf_temp_dE", E, 0.0, dT_dE);
  }

  // --- tier 2: the vulnerability curve and its grid --------------------------
  //
  // Calls exp/pow; no spline read and no solve. The (b, c) sweep is #92's
  // instrument: the knot grid's realised endpoint is a function of the pair, and
  // pinning `psi_max` beside the integral is what would have caught a domain that
  // fell one step short of the one it claimed.
  for (double psi = 0.0; psi <= 6.0 + 1e-12; psi += 0.25) {
    row(kVulnerability, "proportion_of_conductivity", psi, 0.0,
        l.proportion_of_conductivity(psi));
  }
  for (double b : {2.5, 3.898245, 5.0}) {
    for (double c : {1.5, 2.680147, 4.0}) {
      row(kVulnerability, "vulnerability_psi_max", b, c,
          phylloptim::vulnerability_psi_max(b, c));
      row(kVulnerability, "cumulative_vulnerability_integral_limit", b, c,
          phylloptim::cumulative_vulnerability_integral_limit(b, c));
      // The integral itself, at fractions of the curve's own domain so the sweep
      // stays inside it for every pair rather than for the default pair only.
      const double psi_max = phylloptim::vulnerability_psi_max(b, c);
      for (int k = 1; k <= 4; ++k) {
        const double psi = psi_max * double(k) / 5.0;
        row(kVulnerability, "cumulative_vulnerability_integral_at", psi,
            b, phylloptim::cumulative_vulnerability_integral_at(psi, b, c));
      }
      // And the realised knot grid: its LAST knot, which is the value #92 found
      // was not the endpoint the model claimed. Pinned as a number rather than as
      // an equality so that a regression says by how much.
      std::vector<double> x, y;
      phylloptim::cumulative_vulnerability_integral(b, c, 100, x, y);
      row(kVulnerability, "knot_grid_last_psi", b, c, x.back());
      row(kVulnerability, "knot_grid_last_value", b, c, y.back());
      row(kVulnerability, "knot_grid_size", b, c, double(x.size()));
    }
  }

  // --- tier 3: assimilation --------------------------------------------------
  for (double ci = 1.0; ci <= 39.0 + 1e-12; ci += 1.0) {
    row(kAssimilation, "assim_rubisco_limited", ci, 0.0,
        l.assim_rubisco_limited(ci));
    row(kAssimilation, "assim_electron_limited", ci, 0.0,
        l.assim_electron_limited(ci));
    row(kAssimilation, "assim_colimited", ci, 0.0, l.assim_colimited(ci));
  }

  // --- tier 4: the first functions that read a spline ------------------------
  //
  // This is where the interesting money is: the spline is built by odelia, and
  // odelia's interpolator is the one dependency shared by header rather than by
  // copy. `stem_curve_integral` is included alongside `transpiration` because it
  // is the raw spline read and `transpiration` is that scaled -- so a tier-4
  // failure that appears in one and not the other separates the interpolator from
  // the conductance scaling.
  for (double psi_stem = 0.25; psi_stem <= 5.5 + 1e-12; psi_stem += 0.25) {
    row(kSpline, "stem_curve_integral", psi_stem, 0.0,
        l.stem_curve_integral(psi_stem));
    row(kSpline, "stem_curve_integral_deriv", psi_stem, 0.0,
        l.stem_curve_integral_deriv(psi_stem));
    row(kSpline, "transpiration", psi_stem, 0.0, l.transpiration(psi_stem, 0.0));
    row(kSpline, "stom_cond_CO2", psi_stem, 0.0, l.stom_cond_CO2(psi_stem, 0.0));
    row(kSpline, "hydraulic_cost_TF", psi_stem, 0.0,
        l.hydraulic_cost_TF(psi_stem));
  }

  // --- tier 5: primitives containing their own iteration --------------------
  //
  // A disagreement that FIRST appears here, with tiers 1-4 clean, means the
  // primitives agree and an iteration count differs -- which is a different bug
  // from an arithmetic one and wants a different fix.
  for (double psi_stem = 0.25; psi_stem <= 5.5 + 1e-12; psi_stem += 0.25) {
    row(kIterative, "psi_stem_to_ci", psi_stem, 0.0,
        l.psi_stem_to_ci(psi_stem, 0.0));
  }
  // The inverse spline's domain ends at the transpiration reached at psi_crit, and
  // asking past it throws rather than extrapolating. Derived from the model rather
  // than hard-coded, so a trait change cannot silently push this out of range.
  const double e_max = l.transpiration(l.supply_psi_crit(), 0.0);
  for (int i = 1; i < 20; ++i) {
    const double e = e_max * i / 20.0;
    row(kIterative, "transpiration_to_psi_stem", e, 0.0,
        l.transpiration_to_psi_stem(e, 0.0));
  }
  return rows;
}

// --- tolerances -------------------------------------------------------------
//
// Per tier, and a negative value means bit-exact.
//
// ⚠️ THESE ARE CEILINGS DERIVED FROM THE MECHANISM, NOT MEASURED WORST CASES, and
// saying so is the point of this paragraph -- the developer guide records a
// cross-platform table that has been wrong three times because a reading was
// written down as a fact. Tier 1 is libm's exp/pow alone, which disagrees between
// glibc/x86-64 and Apple's libm at a few ULP; every tier above inherits that and
// adds its own amplification, and tier 5 additionally carries a root-find whose own
// tolerance is 1e-10. Each ceiling is orders above what that implies, so a failure
// here is a real change rather than a tight bound being grazed.
//
// The REAL figures come from this program's own summary line, which prints on every
// tolerant run whether or not it passes -- never from the FAIL list, which is
// truncated at 20 and biased toward whichever rows come first.
//
// FIRST READING, ubuntu-latest against a macOS/arm64 file, one CI run, and labelled
// as a reading rather than promoted to a fact:
//
//     tier            g++        clang++    ceiling    worst row
//     arithmetic      3.97e-15   4.21e-15   1e-13      peak_arrh_curve(96, 45)
//     vulnerability   7.01e-16   7.01e-16   1e-12      cumulative_..._at(4.1351, 3.9)
//     assimilation    1.75e-15   1.75e-15   1e-11      assim_colimited(7)
//     spline          3.07e-14   3.07e-14   1e-10      stem_curve_integral_deriv(3.75)
//     iterative       8.01e-11   8.01e-11   1e-6       psi_stem_to_ci(2.75)
//
// The two compilers agree on four tiers exactly and differ only in tier 1, at the
// same operating point -- so what is being measured is the platform's libm, not the
// compiler, which is the same conclusion test_golden reached.
//
// ⚠️ READ THE SHAPE, NOT JUST THE SIZES, because the shape is the result. Four tiers
// sit at 1e-14 to 1e-16 -- libm and FMA contraction, nothing else -- and the
// iterative tier sits three orders higher at 8.01e-11, which is `psi_stem_to_ci`'s
// own 1e-10 tolerance showing through. So the 1.4e-04 that operating_points.tsv
// reports for its argmax-derived fields is NOT present in the primitives at all: it
// is amplification by the solve, and this file is what establishes that rather than
// assuming it.
//
// The ceilings are deliberately NOT tightened to these numbers yet. One run on one
// runner image is not a distribution, and a bound set to the first reading fails on
// the second. Tighten when several agree.
struct Tolerance {
  double per_tier[kTierCount];
};

const Tolerance kExact = {{-1.0, -1.0, -1.0, -1.0, -1.0, -1.0}};

const Tolerance kCrossPlatform = {{-1.0,
                                   /*arithmetic*/ 1e-13,
                                   /*vulnerability*/ 1e-12,
                                   /*assimilation*/ 1e-11,
                                   /*spline*/ 1e-10,
                                   /*iterative*/ 1e-6}};

int generate() {
  const std::vector<Row> rows = run_grid();
  FILE *f = std::fopen(kGoldenPath, "w");
  if (f == nullptr) {
    std::fprintf(stderr, "cannot write %s (run from tests/cpp)\n", kGoldenPath);
    return 2;
  }
  std::fprintf(f, "tier\tfunction\targ1\targ2\tvalue\n");
  for (const Row &r : rows) {
    std::fprintf(f, "%s\t%s\t%.17g\t%.17g\t%.17g\n", tier_name(r.tier),
                 r.fn.c_str(), r.a1, r.a2, r.value);
  }
  std::fclose(f);
  std::printf("wrote %s: %zu primitive values\n", kGoldenPath, rows.size());
  return 0;
}

int compare(const Tolerance &tol) {
  const std::vector<Row> rows = run_grid();
  FILE *f = std::fopen(kGoldenPath, "r");
  if (f == nullptr) {
    std::fprintf(stderr,
                 "cannot read %s -- run `make -C tests/cpp primitives-golden` "
                 "to create it\n",
                 kGoldenPath);
    return 2;
  }
  char line[512];
  if (std::fgets(line, sizeof line, f) == nullptr) {  // header
    std::fprintf(stderr, "%s is empty\n", kGoldenPath);
    std::fclose(f);
    return 2;
  }

  int failures = 0, inexact = 0;
  double worst[kTierCount] = {0};
  std::string worst_where[kTierCount];
  std::size_t i = 0;
  for (; i < rows.size(); ++i) {
    if (std::fgets(line, sizeof line, f) == nullptr) {
      std::fprintf(stderr,
                   "FAIL: golden file has %zu rows, the grid has %zu. A tier was "
                   "added or removed -- regenerate deliberately.\n",
                   i, rows.size());
      ++failures;
      break;
    }
    char tier_buf[64], fn_buf[128];
    double a1 = 0.0, a2 = 0.0, want = 0.0;
    if (std::sscanf(line, "%63s\t%127s\t%lf\t%lf\t%lf", tier_buf, fn_buf, &a1,
                    &a2, &want) != 5) {
      std::fprintf(stderr, "FAIL: cannot parse %s row %zu\n", kGoldenPath, i + 1);
      ++failures;
      continue;
    }
    const Row &r = rows[i];
    // The KEY is checked before the value. A row that has drifted out of
    // alignment would otherwise be reported as a numeric failure, which sends the
    // reader after arithmetic that never changed.
    if (r.fn != fn_buf || std::string(tier_name(r.tier)) != tier_buf) {
      std::fprintf(stderr,
                   "FAIL row %zu: golden says %s/%s, grid says %s/%s. The grid "
                   "was reordered -- regenerate deliberately.\n",
                   i + 1, tier_buf, fn_buf, tier_name(r.tier), r.fn.c_str());
      ++failures;
      continue;
    }
    if (r.value == want) {
      continue;
    }
    // Both non-finite and the same kind counts as agreement: NaN != NaN, and a
    // primitive that returns the NA sentinel is reporting a state, not a value.
    if (!std::isfinite(r.value) && !std::isfinite(want) &&
        std::isnan(r.value) == std::isnan(want)) {
      continue;
    }
    ++inexact;
    const double scale = std::max(1.0, std::fabs(want));
    const double rel = std::fabs(r.value - want) / scale;
    if (rel > worst[r.tier]) {
      worst[r.tier] = rel;
      char buf[256];
      std::snprintf(buf, sizeof buf, "%s(%.6g, %.6g): got %.17g want %.17g",
                    r.fn.c_str(), r.a1, r.a2, r.value, want);
      worst_where[r.tier] = buf;
    }
    const double t = tol.per_tier[r.tier];
    if (t < 0.0 || rel > t) {
      if (failures < 20) {
        std::fprintf(stderr, "FAIL [%s] %s(%.6g, %.6g): got %.17g, want %.17g  (rel %.3g)\n",
                     tier_name(r.tier), r.fn.c_str(), r.a1, r.a2, r.value, want,
                     rel);
      }
      ++failures;
    }
  }
  if (std::fgets(line, sizeof line, f) != nullptr) {
    std::fprintf(stderr, "FAIL: golden file has more rows than the grid (%zu)\n", i);
    ++failures;
  }
  std::fclose(f);

  const bool exact_mode = tol.per_tier[kArithmetic] < 0.0;
  if (failures == 0 && exact_mode) {
    std::printf("primitives: %zu values across %d tiers, all bit-identical\n",
                rows.size(), kTierCount - 1);
    return 0;
  }

  // The per-tier table, printed on success AND on failure. It is the instrument
  // this program exists to provide: which tier moved is the answer, and one pooled
  // number cannot give it.
  std::FILE *out = failures == 0 ? stdout : stderr;
  std::fprintf(out,
               "primitives: %zu values, %d differ%s. Worst by tier:\n", rows.size(),
               inexact, failures == 0 ? " (within tolerance)" : "");
  for (int t = kArithmetic; t < kTierCount; ++t) {
    std::fprintf(out, "  %-14s %9.3g   tolerance %-8.1g %s\n", tier_name(t),
                 worst[t], tol.per_tier[t],
                 worst_where[t].empty() ? "" : worst_where[t].c_str());
  }
  if (failures == 0) {
    return 0;
  }
  std::fprintf(stderr,
               "\n%d mismatches. READ THE TIER TABLE ABOVE FIRST: the LOWEST tier "
               "that moved is the cause, and the ones above it are consequences.\n"
               "  arithmetic moving at all means libm, FMA contraction, or a real "
               "change to a temperature curve.\n"
               "  iterative moving alone means the primitives agree and an "
               "iteration count differs.\n",
               failures);
  if (exact_mode) {
    std::fprintf(stderr,
                 "If the change was intended, regenerate with `make -C tests/cpp "
                 "primitives-golden` and say so in the commit.\n"
                 "If you are not on the platform that generated the file "
                 "(macOS/arm64), do NOT regenerate -- use --cross-platform.\n");
  }
  return 1;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--generate") == 0) {
    return generate();
  }
  Tolerance tol = kExact;
  if (argc > 1) {
    if (std::strcmp(argv[1], "--cross-platform") == 0) {
      tol = kCrossPlatform;
    } else {
      std::fprintf(stderr,
                   "unknown argument '%s'\n"
                   "usage: test_primitives [--cross-platform | --generate]\n",
                   argv[1]);
      return 2;
    }
  }
  return compare(tol);
}
