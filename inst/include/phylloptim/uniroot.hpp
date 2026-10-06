// -*-c++-*-
#ifndef PHYLLOPTIM_UNIROOT_HPP_
#define PHYLLOPTIM_UNIROOT_HPP_

// Really simple wrapper around Boost's 1d root finding with bisection
// method.

#include <boost/math/tools/roots.hpp>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <phylloptim/util.hpp>

namespace phylloptim {
namespace util {

namespace internals {
// True once the bracket's ends agree to 4 DBL_EPSILON relative, or are adjacent
// doubles where it closes on zero, so a root moves with its inputs to roundoff.
inline bool within_roundoff(double a, double b) {
  return std::abs(a - b) <= 4 * DBL_EPSILON * std::max(std::abs(a), std::abs(b)) ||
         std::nextafter(a, b) == b;
}
}

// Wrapper around boost's root finder as a black-box function.
template <typename Function>
double uniroot(Function f, double min, double max, size_t max_iterations) {
  using boost::math::tools::bisect;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root =
      bisect(f, min, max, internals::within_roundoff, it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

// Faster root finder for SMOOTH, monotonic functions (TOMS748 / Brent-like).
//
// Drop-in replacement for uniroot() with the same [min, max] bracketing
// contract (f(min) and f(max) must have opposite signs). For smooth functions
// TOMS748 converges super-linearly (~5-8 evals) versus bisection's ~one bit per
// iteration, so it is attractive for deeply nested, expensive-per-eval solvers.
//
// One gotcha vs bisect: this validates its bracket and THROWS on non-finite or
// same-sign endpoints where boost::bisect returned NaN silently (guard upstream
// if a finite-but-degenerate bracket can occur; see psi_stem_to_ci's NA guard).
template <typename Function>
double uniroot_smooth(Function f, double min, double max,
                      size_t max_iterations) {
  using boost::math::tools::toms748_solve;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root =
      toms748_solve(f, min, max, internals::within_roundoff, it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

// Same solver, for a caller that has ALREADY evaluated both endpoints.
//
// This exists for the collar solve (Leaf::maximise_profit_over_collar), which has
// to evaluate f at both ends anyway -- it must know the signs there to tell an
// interior optimum from one pinned to a constraint, and it must step the wet end
// past dprofit's infeasibility sentinel. Handing those two values back saves 2 of
// ~12 evaluations per solve, and a collar gradient evaluation is not cheap: each
// one carries a nested ci root-find. Worth the overload on a path plant runs
// millions of times.
template <typename Function>
double uniroot_smooth(Function f, double min, double max, double f_min,
                      double f_max, size_t max_iterations) {
  using boost::math::tools::toms748_solve;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = toms748_solve(
      f, min, max, f_min, f_max, internals::within_roundoff, it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

}
}

#endif
