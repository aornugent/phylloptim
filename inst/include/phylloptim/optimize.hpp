// -*-c++-*-
#ifndef PHYLLOPTIM_OPTIMIZE_HPP_
#define PHYLLOPTIM_OPTIMIZE_HPP_

// Robust 1-D function minimiser: Brent's method (golden-section search with
// parabolic interpolation). This is the algorithm behind R's optimize() /
// Forsythe-Malcolm-Moler `fmin`. It converges super-linearly on a smooth
// objective near its optimum but falls back to a golden-section step whenever
// the parabolic fit is unreliable, so it keeps the bracketing robustness that
// the leaf hydraulic solvers depend on (cf. the TOMS748 caveat in uniroot.h:
// here we only ever evaluate STRICTLY INTERIOR points, never the clamped
// bracket endpoints).

#include <phylloptim/uniroot.hpp>
#include <cfloat>  // DBL_EPSILON
#include <cmath>
#include <limits>
#include <vector>

namespace phylloptim {
namespace util {

// Minimise f over [ax, bx] (requires ax <= bx). Returns the argmin; on return
// fmin (if non-null) holds f at the argmin. `tol` is the absolute tolerance on
// the location of the minimum; values below ~sqrt(DBL_EPSILON)*|x| are not
// useful. A direct transcription of R's Brent_fmin (src/appl/fmin.c).
template <typename Function>
double brent_fmin(Function f, double ax, double bx, double tol,
                  double* fmin = nullptr) {
  // c is the squared inverse of the golden ratio (~0.3819660).
  const double c = (3.0 - std::sqrt(5.0)) * 0.5;
  const double eps = std::sqrt(DBL_EPSILON);

  double a = ax, b = bx;
  double v = a + c * (b - a);
  double w = v, x = v;
  double d = 0.0, e = 0.0;
  double fx = f(x);
  double fv = fx, fw = fx;
  const double tol3 = tol / 3.0;

  for (;;) {
    const double xm = (a + b) * 0.5;
    const double tol1 = eps * std::abs(x) + tol3;
    const double tol2 = tol1 * 2.0;

    // Convergence check.
    if (std::abs(x - xm) <= tol2 - (b - a) * 0.5)
      break;

    double p = 0.0, q = 0.0, r = 0.0;
    bool use_golden = true;

    if (std::abs(e) > tol1) {
      // Fit a parabola through (x,fx), (v,fv), (w,fw).
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.0;
      if (q > 0.0)
        p = -p;
      else
        q = -q;
      r = e;
      e = d;
      // Accept the parabolic step only if it is well inside (a,b) and shrinks.
      if (std::abs(p) < std::abs(0.5 * q * r) &&
          p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        const double u = x + d;
        // Keep the new point away from the bracket endpoints.
        if (u - a < tol2 || b - u < tol2)
          d = (x < xm) ? tol1 : -tol1;
        use_golden = false;
      }
    }

    if (use_golden) {
      e = (x < xm) ? (b - x) : (a - x);
      d = c * e;
    }

    // Evaluate f at a point at least tol1 away from x.
    double u;
    if (std::abs(d) >= tol1)
      u = x + d;
    else
      u = (d > 0.0) ? (x + tol1) : (x - tol1);
    const double fu = f(u);

    // Update the bracket and the three best points.
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }

  if (fmin != nullptr)
    *fmin = fx;
  return x;
}

// Golden-section search for the MAXIMUM of a unimodal f over [ax, bx]. Returns
// the argmax (midpoint of the final bracket); terminates when the bracket width
// falls to `tol`. Reuses one interior golden point per iteration, so it costs a
// single new f() evaluation per step after the initial two.
//
// Why this exists alongside brent_fmin: brent_fmin converges faster but its
// parabolic step makes the argmax a *non-smooth* function of the inputs. Where
// the argmax does not feed a gradient (the single-layer leaf optimisers), prefer
// brent_fmin, which was measured ~2.3-2.6x faster there.
//
// ⚠️ **This is NOT the production collar solver.** That is a safeguarded
// root-find on the first-order condition (Leaf::maximise_profit_over_collar);
// this is only its fallback, for the case where neither bracket endpoint has a
// usable gradient. The collar argmax feeds the demographic growth-rate gradient
// and so must vary smoothly with plant state -- and a comparison-based search is
// not the way to get that:
//
//   * A fixed iteration COUNT is not smoothness. Golden section terminates on
//     bracket WIDTH, so it resolves the argmax only to `tol` and the residual
//     offset wanders discontinuously as the comparison sequence flips. Measured:
//     6 distinct answers across 11 trait steps, tread width ~GSS_tol_abs.
//   * That made the argmax piecewise constant at fine scales, so trait
//     derivatives came back exactly zero -- or, for traits in the hydraulic path,
//     smooth, plausible and SIGN-INVERTED.
//   * Solving dprofit == 0 instead resolves the argmax to solver precision and
//     measured ~1000x smoother second differences in a trait, and 24.5% faster
//     (2.65 vs 3.51 us/solve, interleaved at reps=2000).
//
// Keep that constraint in mind before changing the collar solver, and do not
// read this function's existence as evidence that a comparison-based search is
// the safe choice here.
template <typename Function>
double golden_section_max(Function f, double ax, double bx, double tol) {
  const double gr = (std::sqrt(5.0) + 1.0) / 2.0;  // ~1.6180339...
  double a = ax, b = bx;
  double c = b - (b - a) / gr;
  double d = a + (b - a) / gr;
  double fc = f(c);
  double fd = f(d);
  while (std::abs(b - a) > tol) {
    if (fc > fd) {
      b  = d;
      d  = c;
      fd = fc;                 // reuse
      c  = b - (b - a) / gr;
      fc = f(c);               // 1 new eval
    } else {
      a  = c;
      c  = d;
      fc = fd;                 // reuse
      d  = a + (b - a) / gr;
      fd = f(d);               // 1 new eval
    }
  }
  return (a + b) / 2.0;
}

// Maximise f over the CLOSED interval [lo, hi], where the maximum may sit at an
// endpoint and f need not be unimodal, refining by a ROOT-FIND on the first-order
// condition rather than by a bracket-width search. Returns the argmax; `fmax`, if
// non-null, receives f there.
//
// THIS IS THE ONE SOLVER, and the difference between its callers is only which
// derivative is handed in. There was a second, `maximise_over_closed_interval`,
// identical bar the refinement; it went unused when the collar and stem routes
// both moved onto this one and was deleted rather than left as a trap.
//
// ⚠️ `Leaf::maximise_profit_over_collar` IS NOT A THIRD ONE, and the reason it is
// written out rather than calling this is worth knowing before someone folds it
// in. It is this algorithm at `n = 0` PLUS a classification: it reports which of
// five `OperatingPointKind` tags the answer is (interior, pinned wet, pinned dry,
// refused, non-finite gradient), and 42 of 240 feasible golden-grid rows are
// pinned and read that tag. This function has no concept of one -- it returns an
// argmax and nothing about how it got there. Folding them together means either
// giving this an out-parameter every other caller ignores, or losing the
// classification.
//
// ⚠️ WHY A BARE brent_fmin IS NOT THIS FUNCTION. Brent steps in from the bounds,
// so it can return neither endpoint, and it follows one basin, so it cannot see
// past a local maximum. A bracketing optimiser answers "where is the interior
// maximum", which is a different question from "where is the maximum".
//
// ⚠️ WHY THIS EXISTS AND `brent_fmin` REFINEMENT DOES NOT SUFFICE. Brent terminates
// on bracket WIDTH, so it has no stationarity guarantee at all -- nothing in it
// references df. At the package defaults the scan cell is 0.068 MPa and the
// tolerance `(b - a) * 1e-4` is 6.8e-06 MPa, so the argmax is resolved to about
// that. Harmless for the objective, which is flat at its maximum: 6.8e-06 MPa of
// displacement costs ~5e-11 of profit. Fatal for a DERIVATIVE, because a
// derivative divides by a step, and a 1e-06 relative parameter step moves psi* by
// ~3e-06 MPa -- below the resolution. That is where the 0.1855-against-0.0551
// quantisation and a sign-flipped gradient came from.
//
// Three properties preserved from the bracket-width version, each load-bearing:
//   1. both endpoints are evaluated, so a constrained optimum is reachable;
//   2. the scan locates the basin, so a second interior hump cannot hide the
//      global one -- a root-find alone cannot see past a local maximum;
//   3. the grid argmax is ALWAYS kept as a candidate, so this can only improve on
//      the scan and never return something worse than the grid it contains.
//
// `df` is called as `df(x, &ok)`: `ok` false means the derivative is a sentinel
// rather than a value (a shut-down or reversed-gradient exit), and a sentinel must
// never be read as a stationary point. Where the cell ends do not bracket a
// maximum -- df(a) > 0 > df(b) -- the refinement falls back to the width search,
// which is what the scan already guaranteed.
template <typename Function, typename Deriv>
double maximise_over_closed_interval_foc(Function f, Deriv df, double lo, double hi,
                                         int n, double tol, size_t max_iterations,
                                         double* fmax = nullptr) {
  double best_x = lo;
  double best_f = -std::numeric_limits<double>::infinity();
  auto consider = [&](double x, double fx) {
    if (std::isfinite(fx) && fx > best_f) { best_f = fx; best_x = x; }
  };

  if (!(hi > lo)) {
    consider(lo, f(lo));
    if (fmax != nullptr) *fmax = best_f;
    return best_x;
  }

  // ⚠️ `n < 2` MEANS "NO BASIN SCAN", NOT "NO SEARCH". The cell to refine is then
  // the whole interval, which makes this exactly the endpoints-plus-root-find
  // method -- what the collar route has always done, and what a UNIMODAL objective
  // needs and nothing more. The scan is an argument because multi-modality is a
  // property of the configuration, measured: over a 1728-row sweep the only stem
  // objectives with two prominent interior basins are TF24 (21 rows) and JS22, and
  // EVERY one of those rows has the energy balance on. With it off the objective is
  // unimodal, so scanning would cost ~n extra evaluations to confirm what geometry
  // already guarantees -- on the path plant calls millions of times.
  double a = lo, b = hi;
  if (n >= 2) {
    int arg = 0;
    std::vector<double> xs(static_cast<std::size_t>(n) + 1);
    for (int i = 0; i <= n; ++i) {
      xs[static_cast<std::size_t>(i)] = lo + (hi - lo) * double(i) / double(n);
      const double fx = f(xs[static_cast<std::size_t>(i)]);
      if (std::isfinite(fx) && fx > best_f) arg = i;
      consider(xs[static_cast<std::size_t>(i)], fx);
    }
    if (arg == 0 || arg == n) {
      if (fmax != nullptr) *fmax = best_f;
      return best_x;   // a constrained optimum, returned unrefined on purpose
    }
    a = xs[static_cast<std::size_t>(arg - 1)];
    b = xs[static_cast<std::size_t>(arg + 1)];
  } else {
    consider(lo, f(lo));
    consider(hi, f(hi));
  }
  // ⚠️ PROBE INWARD FOR A USABLE DERIVATIVE; DO NOT EVALUATE AT THE ENDS. At the wet
  // end transpiration is ~0, so gc ~ 0 and the implicit-function quotient behind
  // dJ/dpsi is 0/0; at the dry end the supply integral is at its domain edge. Both
  // hand back a SENTINEL rather than a value, and a sentinel read as a derivative
  // makes the bracket test fail, which silently drops the refinement back to a
  // width search -- measured, that leaves |dJ/dpsi| at 4.4e-05 where the root-find
  // reaches the solver floor. The collar route has always stepped in for this
  // reason; it is not incidental and this is the same loop.
  const double width = b - a;
  bool ok_a = false, ok_b = false;
  double da = 0.0, db = 0.0, xa = a, xb = b;
  for (double frac = 0.0; frac < 0.5; frac = (frac == 0.0 ? 1e-6 : frac * 10.0)) {
    xa = a + frac * width;
    da = df(xa, &ok_a);
    if (ok_a && std::isfinite(da)) break;
  }
  for (double frac = 0.0; frac < 0.5; frac = (frac == 0.0 ? 1e-6 : frac * 10.0)) {
    xb = b - frac * width;
    db = df(xb, &ok_b);
    if (ok_b && std::isfinite(db)) break;
  }
  if (ok_a && ok_b && std::isfinite(da) && std::isfinite(db) &&
      da > 0.0 && db < 0.0) {
    const double x = util::uniroot_smooth(
        [&](double v) { bool ok = false; return df(v, &ok); }, xa, xb, da, db, tol,
        max_iterations);
    consider(x, f(x));
  } else {
    double neg = 0.0;
    const double x = brent_fmin([&](double v) { return -f(v); }, a, b,
                                (b - a) * 1e-4, &neg);
    consider(x, -neg);
  }
  if (fmax != nullptr) *fmax = best_f;
  return best_x;
}


}
}

#endif
