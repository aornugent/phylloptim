// -*-c++-*-
#ifndef PHYLLOPTIM_UNIROOT_HPP_
#define PHYLLOPTIM_UNIROOT_HPP_

// Boost's one-dimensional root-finders: `uniroot` (bisection) and `uniroot_smooth`
// (TOMS748). Each takes a bracket with opposite signs at its ends, and running out
// of max_iterations is infeasible ("root_find_iterations"), not a caller error.

#include <boost/math/tools/roots.hpp>
#include <phylloptim/util.hpp>

namespace phylloptim {
namespace util {

namespace internals {
struct uniroot_tol {
  uniroot_tol(double atol_, double rtol_) : atol(atol_), rtol(rtol_) {}
  bool operator()(double a, double b) {
    return std::abs(a - b) < atol + rtol * std::min(std::abs(a), std::abs(b));
  }
  double atol;
  double rtol;
};
}

// Bisection, for a target with no smoothness to exploit. Returns the midpoint of
// the final bracket.
template <typename Function>
double uniroot(Function f, double min, double max, double tol,
               size_t max_iterations) {
  using boost::math::tools::bisect;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = bisect(f, min, max,
                                          internals::uniroot_tol(tol, tol),
                                          it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

// TOMS748, for a target smooth and monotone across the whole bracket: same
// contract as uniroot(), in a handful of evaluations where bisection gains one bit
// per step. Judge smoothness target by target, on the bracket actually handed in:
// the leaf's ci, soil-to-collar and collar residuals qualify. Unlike boost::bisect it
// throws on a non-finite or same-sign bracket, so guard upstream wherever a
// degenerate bracket can occur.
template <typename Function>
double uniroot_smooth(Function f, double min, double max, double tol,
                      size_t max_iterations) {
  using boost::math::tools::toms748_solve;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = toms748_solve(
      f, min, max, internals::uniroot_tol(tol, tol), it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

// Same solver, given f at both ends, returning the final bracket rather than its
// midpoint for a caller that needs the root from one side: f(first) and f(second)
// differ in sign (or one is zero).
template <typename Function>
std::pair<double, double> uniroot_smooth_bracket(Function f, double min,
                                                 double max, double f_min,
                                                 double f_max, double tol,
                                                 size_t max_iterations) {
  using boost::math::tools::toms748_solve;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = toms748_solve(
      f, min, max, f_min, f_max, internals::uniroot_tol(tol, tol), it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return root;
}

// Same solver, for a caller that has already evaluated both ends, as the collar
// solve must to classify its optimum: saves two evaluations, each carrying a
// nested ci root-find.
template <typename Function>
double uniroot_smooth(Function f, double min, double max, double f_min,
                      double f_max, double tol, size_t max_iterations) {
  using boost::math::tools::toms748_solve;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = toms748_solve(
      f, min, max, f_min, f_max, internals::uniroot_tol(tol, tol), it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop_infeasible("root_find_iterations",
                          "exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

}
}

#endif
