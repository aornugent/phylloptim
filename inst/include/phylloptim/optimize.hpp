// -*-c++-*-
#ifndef PHYLLOPTIM_OPTIMIZE_HPP_
#define PHYLLOPTIM_OPTIMIZE_HPP_

// One-dimensional optimisers. `maximise_over_closed_interval_foc` is the one
// solver every route uses. `brent_fmin` and `golden_section_max` terminate on
// bracket width and evaluate only strictly interior points.

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

// Golden-section search for the maximum of a unimodal f over [ax, bx]: returns
// the midpoint of the final bracket once its width falls to `tol`, at one new f()
// per step after the first two.
//
// It exists only as the fallback of `Leaf::maximise_profit_over_collar` where a
// bracket end has no usable gradient, so that case still gets an answer. ⚠️ It
// is not a safe default: terminating on width resolves the argmax only to `tol`,
// so the argmax is piecewise constant in the inputs and a trait derivative through
// it comes back zero or sign-inverted (hazard 3). The collar solve root-finds
// dprofit == 0 instead.
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

// Maximise f over the closed interval [lo, hi], where the maximum may sit at an
// endpoint and f need not be unimodal, refining by a root-find on df == 0. Returns
// the argmax; `fmax`, if non-null, receives f there. This is the one solver, and
// `Leaf::maximise_profit_over_collar` is it at n = 0 plus an `OperatingPointKind`
// tag, which this cannot report.
//
// The root-find is what makes the argmax differentiable: a width-terminated search
// resolves it more coarsely than a 1e-6 parameter step moves it (hazard 15), and
// a bare `brent_fmin` can return neither endpoint nor see past one basin (hazard
// 11). Three properties are load-bearing:
//   1. both endpoints are evaluated, so a constrained optimum is reachable;
//   2. the scan locates the basin, so a second interior hump cannot hide the
//      global one;
//   3. the grid argmax is always a candidate, so the result is never worse than
//      the scan.
//
// `df(x, &ok)` with `ok` false is a sentinel (a shut-down or reversed-gradient
// exit), never a stationary point. Where the cell ends do not bracket a maximum,
// df(a) > 0 > df(b), the refinement falls back to `brent_fmin`.
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

  // `n < 2` means no basin scan, not no search: the cell is the whole interval,
  // endpoints plus a root-find, which is all a unimodal objective needs. The stem
  // routes pass `Leaf::basin_scan_cells()`, 0 unless the energy balance is on.
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
  // Probe inward for a usable derivative, never at the ends: at the wet end gc ~ 0
  // makes dJ/dpsi 0/0, and at the dry end the supply integral is at its domain
  // edge. A sentinel there fails the bracket test and silently drops the
  // refinement to `brent_fmin`, leaving the argmax width-resolved.
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
