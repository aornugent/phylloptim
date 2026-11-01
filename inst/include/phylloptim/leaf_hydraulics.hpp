// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_HYDRAULICS_HPP_
#define PHYLLOPTIM_LEAF_HYDRAULICS_HPP_

// Water through the stem: the residual (cuticular) conductance and its soft
// cap, the stem's vulnerability curve as the four operations performed on
// it, and transpiration, at double and on any scalar.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

// --- the residual (cuticular) conductance ------------------------------------

// The temperature factor on g_min: see the member. Exactly 1 at 25 C, because
// both terms of the exponent are differences of identical expressions there.
inline double Leaf::g_min_factor_compute(double leaf_temp) const {
  const double w = g_min_T_width_;
  if (!(w > 0.0)) {
    util::stop("g_min_T_width_ must be a positive width in K; got " +
               util::to_string(w));
  }
  const double arg =
      g_min_k_low_ * (leaf_temp - 25.0) +
      w * g_min_k_high_ *
          (util::softplus((leaf_temp - g_min_T_trans_) / w) -
           util::softplus((25.0 - g_min_T_trans_) / w));
  const double factor = std::exp(arg);
  g_min_factor_memo_ = {leaf_temp, g_min_k_low_, g_min_k_high_,
                        g_min_T_trans_, g_min_T_width_, factor};
  return factor;
}

inline double Leaf::g_min_factor_temp_deriv(double leaf_temp) const {
  return g_min_factor(leaf_temp) *
         (g_min_k_low_ +
          g_min_k_high_ *
              util::logistic((leaf_temp - g_min_T_trans_) / g_min_T_width_));
}

// kmax Q(psi_up) = kmax (G(inf) - G(psi_up)): the most the stem can ever
// carry from psi_up, which it approaches as its own potential runs away.
inline double Leaf::supply_ceiling(double psi_upstream) const {
  return leaf_specific_conductance_max_ *
         stem_curve_complement(psi_upstream, "Leaf::supply_ceiling");
}

// Past q = residual_cap_q_limit the cap is its limit phi E_inf to within
// q^-p / p, below double precision, and the power would overflow; the limit is
// taken there, on every path, so none of them can form inf * 0.
inline constexpr double residual_cap_q_limit = 1e4;

// Below q^p = 2^-54, 1 + q^p rounds to 1 exactly, so the cap returns r bit for
// bit with slope 1: nearly every solve, and skipping the arithmetic there is
// exact, not an approximation.
inline constexpr double residual_cap_inert_qp = 0x1p-54;

// (1 + q^p)^(-1/p) at p = 8 by three squarings and three square roots, each
// correctly rounded, in place of two pow() calls, which profiled as the largest
// cost a residual conductance added. One expression for the double path and the
// tape, so the two cannot drift.
template <class T>
inline T residual_cap_factor(const T& q, T* one_plus_qp = nullptr) {
  using std::sqrt;
  const T q2 = q * q;
  const T q4 = q2 * q2;
  const T s = T(1.0) + q4 * q4;
  if (one_plus_qp != nullptr) {
    *one_plus_qp = s;
  }
  return T(1.0) / sqrt(sqrt(sqrt(s)));
}

static_assert(Leaf::residual_cap_power == 8.0,
              "residual_cap_factor takes the eighth root by three square roots; "
              "another power needs another expression");

// r* = r (1 + (r/E_inf)^p)^(-1/p), the smooth minimum of the demand and the
// ceiling (see residual_cap_power). Exactly zero at r = 0, and the ceiling
// itself where the ceiling has closed.
inline double Leaf::residual_cap(double r, double E_inf) {
  if (r == 0.0) {
    return 0.0;
  }
  if (!(E_inf > 0.0)) {
    return 0.0;
  }
  const double E_cap = residual_cap_fraction * E_inf;
  const double q = r / E_cap;
  if (q > residual_cap_q_limit) {
    return E_cap;
  }
  const double q2 = q * q, q4 = q2 * q2;
  if (q4 * q4 < residual_cap_inert_qp) {
    return r;
  }
  return r * residual_cap_factor(q);
}

// With q = r / (phi E_inf): dr*/dr = (1 + q^p)^(-1/p - 1), and
// dr*/dE_inf = phi q^(p+1) (1 + q^p)^(-1/p - 1).
inline double Leaf::residual_cap_ddemand(double r, double E_inf) {
  if (r == 0.0 || !(E_inf > 0.0)) {
    return r == 0.0 ? 1.0 : 0.0;
  }
  const double q = r / (residual_cap_fraction * E_inf);
  if (q > residual_cap_q_limit) {
    return 0.0;
  }
  const double q2 = q * q, q4 = q2 * q2;
  if (q4 * q4 < residual_cap_inert_qp) {
    return 1.0;
  }
  double s = 0.0;
  const double f = residual_cap_factor(q, &s);
  return f / s;
}

inline double Leaf::residual_cap_dceiling(double r, double E_inf) {
  if (r == 0.0 || !(E_inf > 0.0)) {
    return 0.0;
  }
  const double q = r / (residual_cap_fraction * E_inf);
  if (q > residual_cap_q_limit) {
    return residual_cap_fraction;
  }
  const double q2 = q * q, q4 = q2 * q2, q8 = q4 * q4;
  if (q8 < residual_cap_inert_qp) {
    return residual_cap_fraction * q8 * q;
  }
  double s = 0.0;
  const double f = residual_cap_factor(q, &s);
  return residual_cap_fraction * q8 * q * f / s;
}

// g_min(T) D(T) / P, in kg: the cuticle's loss at a leaf temperature.
inline double Leaf::residual_demand_at(double leaf_temp) const {
  return g_min_at(leaf_temp) * leaf_vpd_at(leaf_temp) / atm_kpa_ /
         kg_to_mol_h2o;
}

inline double Leaf::residual_demand_dT(double leaf_temp) const {
  const double D = leaf_vpd_at(leaf_temp);
  // The deficit's slope is the saturation curve's, and zero on its floor --
  // the same rule dprofit_energy_balance_term applies to the stomatal part.
  const double dD_dT = (use_energy_balance_ && D > vpd_leaf_min)
                           ? saturation_vapour_pressure_slope(leaf_temp)
                           : 0.0;
  return (g_min * g_min_factor_temp_deriv(leaf_temp) * D +
          g_min_at(leaf_temp) * dD_dT) /
         atm_kpa_ / kg_to_mol_h2o;
}

inline double Leaf::residual_flux_active(double psi_upstream, double E_total) const {
  if (g_min == 0.0) {
    return 0.0;
  }
  const double T = use_energy_balance_ ? leaf_temp_from_E(E_total) : leaf_temp_;
  const double r = residual_demand_at(T);
  return residual_root_cap(residual_cap(r, supply_ceiling(psi_upstream)), r);
}

inline double Leaf::residual_flux_dpsi_upstream_active(double psi_upstream,
                                                       double leaf_temp) const {
  if (g_min == 0.0) {
    return 0.0;
  }
  // dE_inf/dpsi_up = -kmax f(psi_up), the table's slope. The roots' ceiling
  // does not move with the collar, so the second cap contributes only its
  // slope in what the first passes it.
  const double r = residual_demand_at(leaf_temp);
  const double E_inf = supply_ceiling(psi_upstream);
  return residual_root_cap_ddemand(residual_cap(r, E_inf), r) *
         residual_cap_dceiling(r, E_inf) *
         (-leaf_specific_conductance_max_ *
          stem_curve_integral_deriv(psi_upstream));
}

inline double Leaf::residual_cond_temp_deriv_active(double psi_upstream,
                                                    double leaf_temp,
                                                    double gc_const) const {
  if (g_min == 0.0) {
    return 0.0;
  }
  const double r = residual_demand_at(leaf_temp);
  const double E_inf = supply_ceiling(psi_upstream);
  const double r1 = residual_cap(r, E_inf);
  // dr*/dr: through the stem's cap into the roots', plus the roots' ceiling's
  // own dependence on r.
  const double drstar_dr = residual_root_cap_ddemand(r1, r) *
                               residual_cap_ddemand(r, E_inf) +
                           residual_root_cap_dceiling_dr(r1, r);
  return -gc_const * drstar_dr * residual_demand_dT(leaf_temp) +
         g_min_co2_ratio_ * g_min * g_min_factor_temp_deriv(leaf_temp);
}

template <class S>
inline S Leaf::residual_flux_at(const S& collar, const leaf_pars<S>& pars,
                                const SupplyDraw<S>& draw) const {
  // The demand: g_min at the seated leaf temperature, at the seated deficit.
  // Both are the operating point's and passive; the trait carries the row.
  const double per_g = g_min_factor(Tleaf_) * vpd_leaf_ / atm_kpa_ / kg_to_mol_h2o;
  const S r = pars[par_g_min] * S(per_g);
  const S E_inf = pars[par_kmax] * stem_complement_at<S>(collar, pars);
  const S E_cap = S(residual_cap_fraction) * E_inf;
  // The same limits as the double path's residual_cap, so the two agree where
  // the ceiling has closed: the cap is the ceiling's fraction there, and a
  // ceiling of zero (the stem curve at its limit) carries no flux.
  if (!(odelia::util::to_passive(E_inf) > 0.0)) {
    return S(0.0) * r;
  }
  const S q = r / E_cap;
  // The same three regimes as residual_cap, chosen at double: the limit, the
  // inert cap (whose rows in the ceiling are below double precision), and the
  // soft minimum itself.
  auto inert = [](double qd) { const double q2 = qd * qd, q4 = q2 * q2; return q4 * q4 < residual_cap_inert_qp; };
  const S r1 = odelia::util::to_passive(q) > residual_cap_q_limit
                   ? E_cap
                   : inert(odelia::util::to_passive(q))
                         ? r
                         : r * residual_cap_factor<S>(q);
  // The roots' ceiling, the draw's, with its rows; the same limits again.
  if (!draw.has_uptake_limit) {
    return r1;
  }
  if (!(odelia::util::to_passive(draw.uptake_limit) > 0.0)) {
    return S(0.0) * r1;
  }
  // The ceiling shrinks as U (U / (U + r))^2, as residual_root_ceiling_at says.
  const S& U = draw.uptake_limit;
  const S w = U / (U + r);
  const S U_cap = S(residual_cap_fraction) * (U * w * w);
  const S q_root = r1 / U_cap;
  if (odelia::util::to_passive(q_root) > residual_cap_q_limit) {
    return U_cap;
  }
  if (inert(odelia::util::to_passive(q_root))) {
    return r1;
  }
  return r1 * residual_cap_factor<S>(q_root);
}

// transpiration supply functions

// returns proportion of conductance taken from hydraulic vulnerability curve (unitless)
template <typename T, typename U>
inline T Leaf::proportion_of_conductivity_kernel(const T& psi, const U& b,
                                                 const U& c) const {
  return exp(-pow((psi / b), c));
}

template <typename T>
inline T Leaf::proportion_of_conductivity_kernel(T psi) const {
  return proportion_of_conductivity_kernel<T, double>(psi, stem_b, stem_c);
}

inline double Leaf::proportion_of_conductivity(double psi) const {
  return proportion_of_conductivity_kernel(psi);
}

// Seat the stem's vulnerability table at stem_c.
inline void Leaf::setup_transpiration(double resolution) {
  // Recorded because `set_traits` rebuilds both curves at this member: otherwise
  // a re-traited leaf silently differs from a fresh one at those traits (hazard
  // 10). Two curves built at different resolutions record the last.
  vulnerability_curve_ncontrol = resolution;
  // The root's table is the same table wherever the two c agree (the defaults),
  // so it is shared rather than built twice.
  stem_curve_ = vulnerability_table(stem_c, resolution, roots_.vulnerability_);
  ++stem_curve_builds_;
}

// --- the stem curve, as the four operations performed on it ------------------

// Cold and out of line: reached only on a refusal, so its string building stays
// off the hot path.
[[noreturn, gnu::cold, gnu::noinline]] inline void stem_curve_out_of_domain(
    double lo, double hi, double u, const char* what, const char* arg_name,
    const char* caller) {
  const bool below = u < lo;
  util::stop_infeasible("stem_curve_domain",
             std::string("Leaf hydraulics: ") + what +
             " evaluated outside its domain: " + arg_name + " = " +
             util::format_double(u) + " lies " +
             util::format_double(below ? lo - u : u - hi) + " beyond the " +
             (below ? "lower" : "upper") + " end of [" +
             util::format_double(lo) + ", " + util::format_double(hi) + "]" +
             (caller == nullptr ? std::string()
                                : std::string("; asked by ") + caller) +
             ".");
}

inline double Leaf::stem_curve_limit() const {
  return stem_b * stem_curve_->G_inf_over_b();
}

// A NaN fails the comparison and falls through to the read, which returns NaN:
// plant's profit_psi_stem_TF(NA, .) -> NA contract rests on this.
inline double Leaf::stem_curve_integral(double psi, const char* caller) const {
  if (psi < 0.0) {
    stem_curve_out_of_domain(0.0, std::numeric_limits<double>::infinity(), psi,
                             "the cumulative xylem conductivity integral G "
                             "(psi in +MPa)",
                             "psi", caller);
  }
  return vulnerability_integral(*stem_curve_, psi, stem_b, stem_inv_b_.of(stem_b));
}

inline double Leaf::stem_curve_integral_deriv(double psi) const {
  return vulnerability_integral_deriv(*stem_curve_, psi, stem_inv_b_.of(stem_b));
}

inline double Leaf::stem_curve_integral_inverse(double w, const char* caller,
                                                bool) const {
  // A potential carries w exactly when 0 <= w <= G(inf), and G(inf) maps to the
  // end of the curve. The slack is rounding: a flux kept below the ceiling (the
  // residual under its soft cap, a continuity root) can sum a few ulps past it.
  const double G_inf = stem_curve_limit();
  if (w < 0.0 || w > G_inf * (1.0 + 16.0 * std::numeric_limits<double>::epsilon())) {
    stem_curve_out_of_domain(0.0, G_inf, w,
                             "the cumulative xylem conductivity integral G, "
                             "INVERTED (G^-1, argument in E/K_max), whose upper "
                             "end is its limit G(inf)",
                             "E/K_max", caller);
  }
  return vulnerability_integral_inverse(*stem_curve_, std::min(w, G_inf), stem_b);
}

inline double Leaf::stem_curve_integral_inverse_deriv(double w) const {
  // dpsi/dG = 1/G'(psi) at the inverse: the forward slope's reciprocal, so there
  // is one definition of this number.
  return 1.0 / stem_curve_integral_deriv(stem_curve_integral_inverse(w));
}

inline double Leaf::stem_curve_complement(double psi, const char* caller) const {
  if (psi < 0.0) {
    stem_curve_out_of_domain(0.0, std::numeric_limits<double>::infinity(), psi,
                             "what the cumulative xylem conductivity integral "
                             "has left, G(inf) - G (psi in +MPa)",
                             "psi", caller);
  }
  return vulnerability_complement(*stem_curve_, psi, stem_b,
                                  stem_inv_b_.of(stem_b));
}

inline double Leaf::stem_curve_complement_inverse(double Q) const {
  return vulnerability_complement_inverse(*stem_curve_, Q, stem_b);
}

inline bool Leaf::stem_curve_in_tail(double psi) const {
  return vulnerability_in_tail(*stem_curve_, psi, stem_inv_b_.of(stem_b));
}

// Before the tail in G; from an upstream potential in the tail in Q = G(inf) - G,
// where what is left after E is Q(upstream) - E/kmax (hazard 1).
inline bool Leaf::stem_carries(double E, double psi_upstream) const {
  constexpr double slack = 1.0 + 16.0 * std::numeric_limits<double>::epsilon();
  const double w = E / leaf_specific_conductance_max_;
  if (!stem_curve_in_tail(psi_upstream)) {
    return w + stem_curve_integral(psi_upstream) <= stem_curve_limit() * slack;
  }
  return w - stem_curve_complement(psi_upstream) <=
         stem_curve_limit() * (slack - 1.0);
}

inline double Leaf::psi_stem_carrying(double E, double psi_upstream,
                                      const char* caller) const {
  const double w = E / leaf_specific_conductance_max_;
  if (!stem_curve_in_tail(psi_upstream)) {
    return stem_curve_integral_inverse(w + stem_curve_integral(psi_upstream, caller),
                                       caller);
  }
  const double Q_up = stem_curve_complement(psi_upstream, caller);
  const double left = Q_up - w;
  if (std::isnan(left)) {
    return left;  // a NaN flux places no stem, and is not a flux out of range
  }
  // The slack is a few ulps of G(inf), not of Q: a root-find in G's arithmetic
  // can overshoot a remaining capacity of ~1e-15 G(inf), and that flux is the
  // whole capacity. More than G(inf) left is a stem wetter than saturation.
  if (!(left >= -16.0 * std::numeric_limits<double>::epsilon() * stem_curve_limit()) ||
      !(left <= stem_curve_limit())) {
    stem_curve_out_of_domain(0.0, Q_up, w,
                             "the stem's remaining capacity from psi_upstream, "
                             "Q(psi_upstream) = G(inf) - G(psi_upstream), "
                             "which bounds the flux it carries",
                             "E/K_max", caller);
  }
  return stem_curve_complement_inverse(left > 0.0 ? left : 0.0);
}

inline void Leaf::perturb_stem_P50(double stem_P50_new) {
  check_psi_magnitudes(stem_P50_new, stem_c, roots_.root_P50, roots_.root_c);
  stem_P50 = stem_P50_new;
  stem_b = weibull_b_from_P50(stem_P50_new, stem_c);
  psi_crit = weibull_psi_at_fraction(stem_b, stem_c, k_crit_fraction);
  // The transpiration memo is keyed on the potentials, not the curve, so without
  // this a gradient's perturbation returns a plausible value from the old stem_b.
  transpiration_cached_ = false;
}

// Direct quadrature of the vulnerability curve: an independent check on the table
// transpiration() reads, not used in production. Convergent with plant's QAG, not
// bit-identical; see quadrature.hpp.
inline double Leaf::transpiration_full_integration(double psi_stem, double psi_upstream) {
  const auto f = [&](double psi) -> double {
    return proportion_of_conductivity(psi);
  };
  return leaf_specific_conductance_max_ *
         quadrature::adaptive_simpson(f, psi_upstream, psi_stem,
                                      integration_tol_);
}

// transpiration()'s refusal, cold and out of line: its string building inlined
// into that hot path cost the solve measurably.
[[noreturn, gnu::cold, gnu::noinline]] inline void
stop_reversed_flow(double psi_stem, double psi_upstream) {
  util::stop("Leaf::transpiration: psi_stem = " + util::format_double(psi_stem) +
             " MPa is wetter than psi_upstream = " +
             util::format_double(psi_upstream) +
             " MPa, which would move water from the leaf to the soil; "
             "psi_stem must be at least psi_upstream");
}

// Supply-side transpiration through the stem, kg H2O s^-1 m^-2 LA, from two
// positive magnitudes in MPa.
inline double Leaf::transpiration(double psi_stem, double psi_upstream) {

  // One-entry memo: a profit evaluation asks for the same pair several times.
  // Anything that changes the curve must clear it.
  if (transpiration_cached_ &&
      psi_stem == transpiration_cache_psi_stem_ &&
      psi_upstream == transpiration_cache_psi_upstream_) {
    return transpiration_cache_value_;
  }

  // A difference of G, or from an upstream potential in the tail of Q (hazard 1).
  const double E =
      stem_curve_in_tail(psi_upstream)
          ? leaf_specific_conductance_max_ *
                (stem_curve_complement(psi_upstream,
                                       "Leaf::transpiration, at psi_upstream") -
                 stem_curve_complement(psi_stem, "Leaf::transpiration, at psi_stem"))
          : leaf_specific_conductance_max_ *
                (stem_curve_integral(psi_stem, "Leaf::transpiration, at psi_stem") -
                 stem_curve_integral(psi_upstream,
                                     "Leaf::transpiration, at psi_upstream"));

  // A reversed pair is a caller error (solves read it as no flow, through
  // operating_point_flux) that would yield a negative conductance; equality is
  // shut-down and legal. Tested after the lookups, so a negative potential gets
  // the table's message, and before the memo, so the memo never holds one.
  if (psi_stem < psi_upstream) [[unlikely]] {
    stop_reversed_flow(psi_stem, psi_upstream);
  }

  transpiration_cache_psi_stem_ = psi_stem;
  transpiration_cache_psi_upstream_ = psi_upstream;
  transpiration_cache_value_ = E;
  transpiration_cached_ = true;
  return E;
}

inline double Leaf::operating_point_flux(double psi_stem, double psi_upstream) {
  return psi_upstream >= psi_stem ? 0.0 : transpiration(psi_stem, psi_upstream);
}

// The inverse of transpiration(): the psi_stem that carries a known flux from
// psi_upstream, both positive magnitudes.
inline double Leaf::transpiration_to_psi_stem(double transpiration_, double psi_upstream) {
  // Past the end of the curve G^-1 cannot place the stem: a zero flux has it AT
  // the upstream potential, and any positive flux there is infeasible.
  if (transpiration_ == 0.0 &&
      psi_upstream >= vulnerability_psi_end(stem_b, stem_c)) {
    return psi_upstream;
  }

  // A refusal below the table's lower end means a NEGATIVE demanded flux: no
  // stem potential carries it, so the caller should not ask, and widening the
  // domain is not the fix. The message names E/K_max so that reads.
  return psi_stem_carrying(transpiration_, psi_upstream,
                           "Leaf::transpiration_to_psi_stem, inverting for psi_stem");
}

// ⚠️ Hazard 1: the rows are in (P50, c), not (b, c); stem_b is derived, so a
// row in b reaches nothing. With L = ln 2 and b = P50 * L^(-1/c),
//
//     db/dP50 = b / P50            db/dc = b * ln(L) / c^2
//
// The value is the table's, which the solve ran on; the closed form gives only
// the rows.
template <class S>
inline S Leaf::stem_integral_at(const S& psi, const leaf_pars<S>& pars) const {
  // One site serves both curves (closed_form_rows.hpp).
  const double at = odelia::util::to_passive(psi);
  return integral_with_closed_form_rows<S>(
      stem_curve_integral(at, "Leaf::stem_integral_at"),
      stem_curve_integral_deriv(at), psi, pars[par_stem_P50],
      pars[par_stem_c]);
}

template <class S>
inline S Leaf::stem_complement_at(const S& psi, const leaf_pars<S>& pars) const {
  const double at = odelia::util::to_passive(psi);
  return complement_with_closed_form_rows<S>(
      stem_curve_complement(at, "Leaf::stem_complement_at"),
      -stem_curve_integral_deriv(at), psi, pars[par_stem_P50],
      pars[par_stem_c]);
}

template <typename T>
inline T Leaf::transpiration_at(const T& sigma, const T& collar,
                                const leaf_pars<T>& pars) const {
  // The difference transpiration() forms, on the same test, so the value is
  // its value bit for bit and the rows are those of what it differences.
  if (stem_curve_in_tail(odelia::util::to_passive(collar))) {
    return pars[par_kmax] * (stem_complement_at<T>(collar, pars) -
                             stem_complement_at<T>(sigma, pars));
  }
  return pars[par_kmax] * (stem_integral_at<T>(sigma, pars) -
                           stem_integral_at<T>(collar, pars));
}

} // namespace phylloptim

#endif
