// -*-c++-*-
#ifndef PHYLLOPTIM_QUADRATURE_HPP_
#define PHYLLOPTIM_QUADRATURE_HPP_

// Adaptive Simpson quadrature, used only by Leaf::transpiration_full_integration:
// the independent check that the pre-integrated table `transpiration()` reads is
// faithful (`test_spline_matches_direct_integration`). Off the hot path, and kept
// header-only so the package needs no compiled quadrature library. It agrees with
// other quadratures to tolerance, not bit for bit.

#include <cmath>
#include <cstddef>

namespace phylloptim {
namespace quadrature {

namespace internals {

// One adaptive Simpson step over [a, b]. `fa`/`fm`/`fb` are f at the endpoints
// and midpoint (passed in so they are never recomputed), `whole` the Simpson
// estimate over the interval from those three. Recurses until the Richardson
// error estimate on the two halves falls under `tol`, or the depth budget runs
// out.
template <typename Function>
double simpson_step(Function f, double a, double b, double fa, double fm,
                    double fb, double whole, double tol, int depth) {
  const double m = 0.5 * (a + b);
  const double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
  const double flm = f(lm), frm = f(rm);
  const double h = (b - a) / 12.0;
  const double left = h * (fa + 4.0 * flm + fm);
  const double right = h * (fm + 4.0 * frm + fb);
  const double both = left + right;
  // Richardson: the Simpson error scales as h^4, so (both - whole)/15 estimates
  // the error remaining in `both`.
  if (depth <= 0 || std::abs(both - whole) <= 15.0 * tol) {
    return both + (both - whole) / 15.0;
  }
  return simpson_step(f, a, m, fa, flm, fm, left, 0.5 * tol, depth - 1) +
         simpson_step(f, m, b, fm, frm, fb, right, 0.5 * tol, depth - 1);
}

} // namespace internals

// Integrate f over [a, b] to absolute tolerance `tol`. Reversed limits are
// handled by the usual sign flip; a degenerate interval integrates to zero.
template <typename Function>
double adaptive_simpson(Function f, double a, double b, double tol = 1e-8,
                        int max_depth = 30) {
  if (a == b) {
    return 0.0;
  }
  if (a > b) {
    return -adaptive_simpson(f, b, a, tol, max_depth);
  }
  const double m = 0.5 * (a + b);
  const double fa = f(a), fm = f(m), fb = f(b);
  const double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
  return internals::simpson_step(f, a, b, fa, fm, fb, whole, tol, max_depth);
}

} // namespace quadrature
} // namespace phylloptim

#endif
