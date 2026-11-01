// -*-c++-*-
#ifndef PHYLLOPTIM_DAMAGE_HPP_
#define PHYLLOPTIM_DAMAGE_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/util.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace phylloptim {

// The thermal-damage dynamics, written once on the scalar type, so the solve,
// damage_rates(), damage_trajectory() and its gradient share one set of expressions.
namespace damage {

// A forward-mode dual number along one input; comparisons read the value.
struct Dual {
  double v = 0.0;
  double d = 0.0;
  Dual() = default;
  Dual(double value) : v(value), d(0.0) {}  // NOLINT: implicit by design
  Dual(double value, double slope) : v(value), d(slope) {}
};
inline double val(double x) { return x; }
inline double val(const Dual& x) { return x.v; }
inline Dual operator+(const Dual& a, const Dual& b) { return {a.v + b.v, a.d + b.d}; }
inline Dual operator-(const Dual& a, const Dual& b) { return {a.v - b.v, a.d - b.d}; }
inline Dual operator-(const Dual& a) { return {-a.v, -a.d}; }
inline Dual operator*(const Dual& a, const Dual& b) {
  return {a.v * b.v, a.d * b.v + a.v * b.d};
}
inline Dual operator/(const Dual& a, const Dual& b) {
  return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)};
}
inline Dual operator+(const Dual& a, double b) { return {a.v + b, a.d}; }
inline Dual operator+(double a, const Dual& b) { return {a + b.v, b.d}; }
inline Dual operator-(const Dual& a, double b) { return {a.v - b, a.d}; }
inline Dual operator-(double a, const Dual& b) { return {a - b.v, -b.d}; }
inline Dual operator*(const Dual& a, double b) { return {a.v * b, a.d * b}; }
inline Dual operator*(double a, const Dual& b) { return {a * b.v, a * b.d}; }
inline Dual operator/(const Dual& a, double b) { return {a.v / b, a.d / b}; }
inline Dual operator/(double a, const Dual& b) {
  return {a / b.v, -a * b.d / (b.v * b.v)};
}
inline Dual exp(const Dual& a) {
  const double e = std::exp(a.v);
  return {e, e * a.d};
}
inline Dual log1p(const Dual& a) { return {std::log1p(a.v), a.d / (1.0 + a.v)}; }
inline double exp(double a) { return std::exp(a); }
inline double log1p(double a) { return std::log1p(a); }

// R's leaf_damage_gradient() names these by position through
// `Leaf::damage_par_names()`; the three initial states follow them.
enum : int {
  dp_k_i, dp_k_i_Ea, dp_k_i_Tiso, dp_k_rec, dp_k_rec_Ea, dp_k_mat, dp_k_mat_Ea,
  dp_acclim_alpha, dp_acclim_beta, dp_acclim_T, dp_acclim_s, dp_dTimpair_max,
  dp_K_A, dp_topt_offset,
  n_dyn_pars
};
template <class T> using dyn_pars = std::array<T, n_dyn_pars>;

// Reference temperature of repair and maturation, deg C; impairment's is k_i_Tiso_.
inline constexpr double damage_Tref = 25.0;

// The smooth max(0, x) at sharpness s, ln(1 + e^(s x)) / s, unclamped: e^(s x)
// underflows harmlessly below zero, and above 30 a non-overflowing form takes over.
template <class T> inline T softplus(const T& x, const T& s) {
  const T z = s * x;
  return (val(z) > 30.0 ? z + log1p(exp(-z)) : log1p(exp(z))) / s;
}

// Rate k at Tref_K (kelvin) on a centred Arrhenius, Ea in J mol^-1; exactly k at
// Ea = 0. Unclamped: a clamp binds in ordinary use (high Ea below 25 C) and zeroes
// the derivative there. Only an overflow is refused.
template <class T, class U>
inline T centred_arrhenius(const T& k, const T& Ea, const U& Tref_K,
                           double leaf_temp) {
  const T arg =
      (Ea / gas_constant) * (1.0 / Tref_K - 1.0 / (leaf_temp + C_to_K));
  const T rate = k * exp(arg);
  if (std::isinf(val(rate))) {
    util::stop("thermal damage: a centred Arrhenius rate overflowed at leaf "
               "temperature " + util::format_double(leaf_temp) + " C (Ea " +
               util::format_double(val(Ea)) + " J/mol, reference " +
               util::format_double(val(Tref_K)) + " K)");
  }
  return rate;
}

// The acclimation shift of the impairment's isokinetic temperature, smooth through
// A = 0; its pole, A = -K_A, is refused where a state is seated.
template <class T> inline T delta_impair(const T& A, const dyn_pars<T>& p) {
  return p[dp_topt_offset] + p[dp_dTimpair_max] * A / (A + p[dp_K_A]);
}

// Impairment, functional -> impaired: Bison & Michaletz (2026)'s PSII inactivation
// about their isokinetic temperature, where k_i_ is the rate and k_i_Ea_ alone
// sets the steepness. That temperature, shifted, must lie in [250, 350] K.
template <class T>
inline T k_imp(double leaf_temp, const T& delta, const dyn_pars<T>& p) {
  const T Tiso_K = p[dp_k_i_Tiso] + C_to_K + delta;
  if (val(Tiso_K) < 250.0 || val(Tiso_K) > 350.0) {
    util::stop("thermal damage: the impairment's isokinetic temperature, "
               "k_i_Tiso_ shifted by the acclimation, must lie in [250, 350] K; "
               "got " + util::format_double(val(Tiso_K)) + " K");
  }
  return centred_arrhenius(p[dp_k_i], p[dp_k_i_Ea], Tiso_K, leaf_temp);
}

// Repair, impaired -> functional.
template <class T> inline T k_rec_eff(double leaf_temp, const dyn_pars<T>& p) {
  return centred_arrhenius(p[dp_k_rec], p[dp_k_rec_Ea], damage_Tref + C_to_K,
                           leaf_temp);
}

// Maturation, impaired -> damaged.
template <class T> inline T k_mat_eff(double leaf_temp, const dyn_pars<T>& p) {
  return centred_arrhenius(p[dp_k_mat], p[dp_k_mat_Ea], damage_Tref + C_to_K,
                           leaf_temp);
}

// {dA/dt, dI/dt, dD/dt}, per day.
template <class T>
inline std::array<T, 3> rates(double leaf_temp, const T& A, const T& I,
                              const T& D, const dyn_pars<T>& p) {
  const T F = 1.0 - I - D;
  const T k_i = k_imp(leaf_temp, delta_impair(A, p), p);
  const T k_rec = k_rec_eff(leaf_temp, p);
  const T k_mat = k_mat_eff(leaf_temp, p);
  const T dA = p[dp_acclim_alpha] *
                   softplus(T(leaf_temp) - p[dp_acclim_T], p[dp_acclim_s]) -
               p[dp_acclim_beta] * A;
  const T dI = k_i * F - (k_rec + k_mat) * I;
  const T dD = k_mat * I;
  return {dA, dI, dD};
}

// RK4 over a piecewise-linear leaf temperature, calling `put(i, y)` at every
// time. The sub-step count is chosen from values, so a Dual run takes the double
// run's steps and differentiates the discrete map it evaluates.
template <class T, class Put>
inline void trajectory(const std::vector<double>& times,
                       const std::vector<double>& leaf_temp,
                       std::array<T, 3> y, const dyn_pars<T>& p, Put put) {
  const std::size_t n = times.size();
  put(0, y);
  for (std::size_t i = 1; i < n; ++i) {
    const double t0 = times[i - 1], dt_all = times[i] - t0;
    const double T0 = leaf_temp[i - 1], T1 = leaf_temp[i];
    // The stiffest rate on the interval bounds the sub-step: every rate is monotone
    // in temperature, and impairment falls as its isokinetic shift rises.
    const T delta_min =
        p[dp_topt_offset] +
        (val(p[dp_dTimpair_max]) < 0.0 ? p[dp_dTimpair_max] : T(0.0));
    double fastest = val(p[dp_acclim_beta]);
    for (double Tx : {T0, T1}) {
      fastest = std::max(fastest, val(k_imp(Tx, delta_min, p)) +
                                      val(k_rec_eff(Tx, p)) +
                                      val(k_mat_eff(Tx, p)) +
                                      val(p[dp_acclim_beta]));
    }
    // Fastest rate times step below 0.1, and at most 0.25 K of temperature per step:
    // the impairment rate is steep enough to be under-resolved across a larger one.
    const int m = std::max({1, static_cast<int>(std::ceil(dt_all * fastest / 0.1)),
                            static_cast<int>(std::ceil(std::abs(T1 - T0) / 0.25))});
    const double h = dt_all / m;
    auto temp_at = [&](double t) { return T0 + (T1 - T0) * (t - t0) / dt_all; };
    for (int k = 0; k < m; ++k) {
      const double t = t0 + k * h;
      const std::array<T, 3> k1 = rates(temp_at(t), y[0], y[1], y[2], p);
      const std::array<T, 3> k2 =
          rates(temp_at(t + h / 2), y[0] + h / 2 * k1[0], y[1] + h / 2 * k1[1],
                y[2] + h / 2 * k1[2], p);
      const std::array<T, 3> k3 =
          rates(temp_at(t + h / 2), y[0] + h / 2 * k2[0], y[1] + h / 2 * k2[1],
                y[2] + h / 2 * k2[2], p);
      const std::array<T, 3> k4 = rates(temp_at(t + h), y[0] + h * k3[0],
                                        y[1] + h * k3[1], y[2] + h * k3[2], p);
      for (int j = 0; j < 3; ++j) {
        y[j] = y[j] + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
      }
    }
    put(i, y);
  }
}

}  // namespace damage

}  // namespace phylloptim

#endif
