// -*-c++-*-
#ifndef PHYLLOPTIM_UTIL_HPP_
#define PHYLLOPTIM_UTIL_HPP_

// Minimal, R-free utilities for the leaf model (hazard 9): `stop()` throws a
// std::runtime_error, which Rcpp turns into an R error at the package boundary,
// and the missing-value sentinel is a quiet NaN rather than NA_REAL.

#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace phylloptim {
namespace util {

// Throws; never returns, so callers whose catch blocks end in util::stop() raise
// no -Wreturn-type warning.
[[noreturn]] inline void stop(const std::string &msg) {
  throw std::runtime_error(msg);
}

// --- infeasibility, as distinct from a caller error --------------------------
//
// A fit's proposals reach operating points the solve cannot handle, and those must
// cost their rows, not the likelihood. Rcpp carries only `what()` into R, where a
// derived C++ type is invisible, so R/conditions.R re-signals on the token:
//
//     [phylloptim:infeasible:collar_bracket] find_root_psi(...) failed: ...
//
// Set here because the R boundary is RcppR6-generated and this header stays R-free
// (hazard 9). Marking a caller error infeasible is the one mistake that matters:
// the fit's tryCatch would swallow it and report a likelihood over the surviving
// rows. So the default is `stop()`, and a site is infeasible only when a
// well-formed call with in-range data reaches it; R/conditions.R lists the codes.
struct infeasible_error : std::runtime_error {
  explicit infeasible_error(const std::string &what_arg)
      : std::runtime_error(what_arg) {}
};

inline std::string infeasible_token(const std::string &code) {
  return "[phylloptim:infeasible:" + code + "] ";
}

[[noreturn]] inline void stop_infeasible(const std::string &code,
                                         const std::string &msg) {
  throw infeasible_error(infeasible_token(code) + msg);
}

template <typename T> std::string to_string(T x) { return std::to_string(x); }

// A double in an error message, to six significant figures as odelia's and
// plant's util::format_double do; `to_string` gives six decimals, so 1e-22 reads
// "0.000000". For reading, not for reconstructing a double.
inline std::string format_double(double x) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.6g", x);
  return std::string(buf);
}

inline bool is_finite(double x) { return std::isfinite(x); }

// Use this to be explicit when a deliberate floating-point equality test is
// being made (e.g. checking a value is exactly the one we set).
inline bool identical(double a, double b) { return a == b; }

template <typename T> T clamp(T x, T min_val, T max_val) {
  return std::max(std::min(x, max_val), min_val);
}

// log(1 + exp(x)) and its slope, the logistic, in the forms that do not
// overflow. The smooth replacement for a max(0, x), at a width the caller
// scales x by; every soft floor and soft cap in the leaf is one of these.
inline double softplus(double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}
inline double logistic(double x) {
  return x > 0.0 ? 1.0 / (1.0 + std::exp(-x))
                 : std::exp(x) / (1.0 + std::exp(x));
}

// Missing-value sentinel for uninitialised leaf state, in place of R's NA_REAL:
// nothing reads NA_REAL's payload, only !isfinite(), so a quiet NaN is exact.
inline constexpr double na_value = std::numeric_limits<double>::quiet_NaN();

// Integer counterpart, for max_soil_layer before any layer is seen. Matches
// R's NA_INTEGER so the value round-trips unchanged through plant's bindings.
inline constexpr int na_value_int = -2147483648;

} // namespace util
} // namespace phylloptim

#endif
