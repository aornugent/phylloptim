// -*-c++-*-
#ifndef PHYLLOPTIM_UTIL_HPP_
#define PHYLLOPTIM_UTIL_HPP_

// Minimal, R-free utilities for the leaf model.
//
// This replaces plant/util.h, which reached into Rcpp for `stop()` and into
// R_ext/Arith.h for the NA_REAL sentinel. Neither is needed: `stop()` throws a
// std::runtime_error, which Rcpp converts into an ordinary R error at the
// package boundary, and the sentinel is a quiet NaN. Keeping this header free
// of R is what lets the leaf model compile and run as plain C++ (see
// tests/cpp/), which in turn is what makes it testable and profilable without
// an R session.

#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace phylloptim {
namespace util {

// Throws; never returns. The attribute lets the compiler see that callers whose
// catch blocks end in util::stop() do not fall through, which avoids spurious
// -Wreturn-type warnings.
[[noreturn]] inline void stop(const std::string &msg) {
  throw std::runtime_error(msg);
}

// --- infeasibility, as distinct from a caller error (#57) --------------------
//
// A parameter proposal during a fit WILL reach operating points the solve cannot
// handle. That is not an edge case, it is what an optimiser does, and it has to
// cost those rows rather than the whole likelihood evaluation. So a caller needs to
// tell "this operating point cannot be solved" apart from "you handed me a
// mismatched driver vector", and until now both arrived as a plain R error and the
// only way to separate them was to match on message text.
//
// ⚠️ THE TOKEN IS THE LOAD-BEARING HALF, NOT THE TYPE. Rcpp converts any
// std::exception into an ordinary R error carrying only `what()`, so a derived C++
// type is invisible from R -- and R is where the caller is. The token is a stable,
// machine-readable prefix inside the message:
//
//     [phylloptim:infeasible:collar_bracket] find_root_psi(...) failed: ...
//
// which R's thin wrapper reads to re-signal with a condition class. The type is here
// for a C++ consumer, who can catch it directly, and because it costs nothing.
//
// ⚠️ WHY THIS IS NOT AT THE R BOUNDARY, WHERE IT WOULD BE TIDIER. `util.hpp` must
// not reach for Rcpp -- that is hazard 9, and it is enforced by CI building this on
// runners with no R -- and the boundary itself is RcppR6-generated and must not be
// hand-edited. So the classification happens where the knowledge is (here) and the
// translation happens in hand-written R.
//
// ⚠️ AND WHY THE MESSAGE TEXT IS NOT ALREADY ENOUGH. It reads like it would be. It
// is not: the out-of-domain wording has been rewritten twice recently (#65, #79) and
// #92 added another, so a caller matching on prose is silently broken by an
// improvement to an error message. The token is ours and is asserted by tests.
//
// ⚠️ CLASSIFYING A CALLER ERROR AS INFEASIBLE IS THE ONE MISTAKE THAT MATTERS.
// #39 rejected silent NA rows because an all-NA driver column once hid behind "every
// prediction is NA". An input-validation failure marked infeasible would be
// swallowed by the very tryCatch this exists to enable, and the fit would report a
// plausible likelihood over the rows that survived. So the default is `stop()`, and a
// site becomes infeasible only when it is reachable at a WELL-FORMED call with
// in-range data because the parameters or state make the operating point unsolvable.
// The list of codes that have passed that test is in R/conditions.R.
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

// A double in an error message. to_string above is std::to_string, which for a
// double is fixed-point with six DECIMALS -- it renders a transpiration of 1e-22
// as "0.000000" and a hydraulic potential of 6.8918 with trailing zeros. Six
// significant figures instead, matching odelia's and plant's util::format_double
// so the family renders numbers the same way. For reading, not for reconstructing
// a double.
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

// Missing-value sentinel for uninitialised leaf state. Stands in for R's
// NA_REAL, which is itself a NaN with a payload; nothing in the leaf model
// inspects the payload, it only ever tests !isfinite(), so a quiet NaN is an
// exact substitute.
inline constexpr double na_value = std::numeric_limits<double>::quiet_NaN();

// Integer counterpart, for max_soil_layer before any layer is seen. Matches
// R's NA_INTEGER so the value round-trips unchanged through plant's bindings.
inline constexpr int na_value_int = -2147483648;

} // namespace util
} // namespace phylloptim

#endif
