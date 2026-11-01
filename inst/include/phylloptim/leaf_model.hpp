// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_MODEL_HPP_
#define PHYLLOPTIM_LEAF_MODEL_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/uniroot.hpp>
#include <phylloptim/optimize.hpp>
#include <phylloptim/quadrature.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/single_potential.hpp>
#include <phylloptim/vulnerability.hpp>
#include <phylloptim/closed_form_rows.hpp>
#include <phylloptim/clamp_sites.hpp>

#include <odelia/interpolator.hpp>
#include <odelia/with_slope.hpp>
#include <odelia/tangent.hpp>
#include <odelia/implicit_node.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <XAD/XAD.hpp>

// The leaf.
//
// `Leaf` holds a leaf's traits, its environment and its water-supply path, and
// find_root_collar_psi() finds its operating point: the root-collar water
// potential that maximises carbon profit (assimilation through a benefit link,
// minus a hydraulic cost) subject to continuity from soil to stem. The answer is
// left in members (opt_root_psi_, opt_psi_stem_, ci_, profit_, E_up_, ...). The
// cost curves share that solve and differ only in their cost (see CostCurve).
//
// THE DERIVATIVE SURFACE. After the solve has run at double, a consumer takes
// the operating point's derivative, at whatever scalar it holds, from three
// member templates:
//
//   SupplyDraw<S>  draw = leaf.supply_draw_at<S>(S(leaf.opt_root_psi_), supply);
//   HeldCoords<S>  held;
//   S              psi  = leaf.operating_collar<K, S>(draw, pars, curvature, &held);
//   BilledOutputs<S> out  = leaf.billed_outputs_at<K, S>(psi, draw, pars, &held);
//
// `out.profit` and each `out.uptake[j]` then carry their derivatives in `pars`
// and `supply`, from the implicit function theorem at the converged point;
// nothing the root-find did is recorded. S is double, a tangent
// (odelia::ode::tangent_scalar) or an adjoint (odelia::ode::adjoint_scalar).
// `held` is optional: it spares billed_outputs_at solving the stem potential and
// ci again, and changes no value or tangent (adjoint rows only to rounding).
// Preconditions, checked where marked:
//
//   - the draw is taken at the solved collar (checked: operating_collar and
//     billed_outputs_at refuse another);
//   - `held`, if passed, comes from operating_collar at this point (checked:
//     billed_outputs_at refuses one placed at another collar or kind);
//   - nothing moves the leaf between the solve and billed_outputs_at;
//   - `pars` and `supply` hold, in value, the seated traits (passive_pars()) and
//     the seated soil state and root network;
//   - K is the curve the solve used, TF24 or TF24_floor (checked against the
//     seated curve; see check_surface_curve for the one solve it cannot see);
//   - the supply path is multi-layer (checked: supply_draw_at refuses another);
//   - `curvature` is profit_curvature_at_collar<K>() at an Interior point, and is
//     not read at any other kind or at double; operating_collar refuses one that
//     is not negative or lies inside Leaf::interior_curvature_floor;
//   - the thermal-damage layer is off (checked): the kernels read do not apply it.
//
// operating_point_kind() says which condition placed the point. Branch on it,
// never on a returned value. What operating_collar does at each:
//
//   Interior            stationarity of profit in the collar
//   BoundarySoil        uptake equal to the residual flux at the wet bound
//                       (polished root)
//   BoundaryDry         the stem reaching its dry end: psi_crit for the curves
//                       bounded there, the end of the stem curve for the rest
//   HydraulicShutdown   the same wet bound, reached with no interval to
//                       optimise over
//   Unsolved, Prescribed, SolverRefused, NonFiniteGradient: refused
//
// replay_operating_point(collar, kind) puts a recorded point back with the rows
// its solve produced (hazard 17). profit_row_at_current_point in
// tests/cpp/test_leaf.cpp is a complete, compiled use of the three calls; the
// vignette "Differentiating the operating point" explains the design.
//
// "Hazard N" below is the numbered hazard in the repository's .claude/hazards.md,
// which is not shipped with the package.

namespace phylloptim {

// --- the parameter enumeration, which R indexes into -------------------------
//
// The sixteen traits in `Leaf::set_traits`' argument order, then the non-traits:
// `kmax`, `resistance`, three model-owned slots, the thermal-damage inputs that
// reach the solve, and `PPFD`. Everything before `PPFD` is what a calibration
// fits, `gradient::n_theta`; `n_leaf_pars` also counts `PPFD`, so size `theta` by
// `n_theta` (hazard 12).
//
// ⚠️ R indexes these positions, so a reordering silently differentiates the wrong
// parameter (hazard 13); `test-gradient-batch.R` compares the names with R's.
inline constexpr int n_traits = 16;
inline constexpr int n_leaf_pars = 36;

// The first `n_traits` are `set_traits`' arguments in its order, which is also
// `leaf_traits()`'; the non-traits after them take a relative step.
inline constexpr int par_vcmax_25 = 0;
inline constexpr int par_stem_c = 1;
inline constexpr int par_stem_P50 = 2;
inline constexpr int par_root_c = 3;
inline constexpr int par_root_P50 = 4;
inline constexpr int par_TF24_beta2 = 5;
inline constexpr int par_jmax_25 = 6;
inline constexpr int par_phi_J = 7;
inline constexpr int par_curv_fact_elec_trans = 8;
inline constexpr int par_curv_fact_colim = 9;
inline constexpr int par_TF24_cost_scale = 10;
inline constexpr int par_R_d_25 = 11;
inline constexpr int par_JS22_gamma = 12;
inline constexpr int par_CMax_a = 13;
inline constexpr int par_CMax_b = 14;
// Residual (cuticular) conductance to water vapour at 25 C, mol H2O m^-2 s^-1.
inline constexpr int par_g_min = 15;
// The non-traits. A new fitted slot goes before `par_PPFD`: `gradient::n_theta`
// is `par_PPFD`, so a slot appended last is outside every caller's theta.
inline constexpr int par_kmax = 16;
inline constexpr int par_resistance = 17;
// The model-owned slots, each one curve's parameter: R's `.gradient_model_pars()`
// offers each for its curve alone and refuses it by name elsewhere (hazard 13).
// CF77's marginal value of water.
inline constexpr int par_CF77_lambda = 18;
// TF24_floor's price of water at zero transpiration.
inline constexpr int par_TF24_floor_lambda_o = 19;
// Least-cost's unit-cost ratio beta (Prentice et al. 2014): dimensionless, not a
// price. `gradient::step_for` takes a relative step only while |beta| >= 1; below
// that its absolute step is a large perturbation (2x at beta = 0.5).
inline constexpr int par_LeastCost_beta = 20;
// The thermal-damage inputs that reach the solve, refused unless
// use_thermal_damage_ is on. Not the dynamics parameters: they move only
// damage_rates(), so a slot would read an exact zero; their gradient is
// damage_trajectory_gradient().
inline constexpr int par_A = 21;
inline constexpr int par_I = 22;
inline constexpr int par_D = 23;
inline constexpr int par_topt_offset = 24;
inline constexpr int par_dTopt_max = 25;
inline constexpr int par_K_A = 26;
inline constexpr int par_k_i = 27;
inline constexpr int par_k_i_baseline = 28;
inline constexpr int par_k_rec = 29;
inline constexpr int par_k_rec_Ea = 30;
inline constexpr int par_c_acclim_maint = 31;
inline constexpr int par_c_repair_maint = 32;
inline constexpr int par_c_protect_maint = 33;
inline constexpr int par_c_repair_flux = 34;
// The intercepted PPFD per unit leaf area (see PPFD_): a driver, never fitted, and
// light's only route to carbon, so a kernel reading the `PPFD_` member instead
// gives every trait that moves light a zero row.
inline constexpr int par_PPFD = 35;
static_assert(n_leaf_pars == par_PPFD + 1,
              "the pack ends at par_PPFD; move n_leaf_pars with it");


// The differentiable parameters, as one array rather than a struct, which would
// spell set_traits' order a fourth time beside leaf_traits() and par_names(). It
// carries the traits, not what derives from them: stem_b and psi_crit are derived
// where used, so a direction in stem_P50 reaches them by the chain rule.
template <class S> using leaf_pars = std::array<S, n_leaf_pars>;

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

class Leaf {
public:
  Leaf();
  
  Leaf(double vcmax_25, 
       double stem_c, 
       double stem_P50, 
       double root_c,
       double root_P50,
       double TF24_beta2, 
       double jmax_25, 
       double phi_J, 
       double curv_fact_elec_trans, 
       double curv_fact_colim,
       double GSS_tol_abs,
       double vulnerability_curve_ncontrol,
       double ci_abs_tol,
       double ci_niter,
      double TF24_cost_scale);

  // The stem's cumulative vulnerability integral G, a table in x at stem_c alone
  // (VulnerabilityTable), shared between copies; read only through stem_curve_*.
  std::shared_ptr<const VulnerabilityTable> stem_curve_;
  InverseScale stem_inv_b_;

  // Instrumentation: how many times the stem curve has been built, never reset, so
  // a caller takes differences. A rebuild is otherwise invisible: it is bit-identical.
  long stem_curve_builds_ = 0;

  // The soil -> root-collar supply. Held by value because Leaf must stay copyable,
  // and public because plant's RcppR6 bindings reach its fields (hazard 7).
  MultiLayerRoots roots_;

  // The alternative supply path. supply_kind_ selects between the two members by
  // a branch, not a std::variant, which measured slower (hazard 6).
  SinglePotential single_;
  enum class SupplyKind { MultiLayer, SinglePotential };
  // MultiLayer by default, which is what plant runs.
  SupplyKind supply_kind_ = SupplyKind::MultiLayer;

  // --- choosing the supply path ----------------------------------------------
  //
  // Two entry points rather than a settable tag, which would leave the other path's
  // state configured, ignored and then stale. Each leaves tag and supply agreeing
  // and clears the solved state, so call set_physiology() again.
  void set_supply_multilayer() {
    supply_kind_ = SupplyKind::MultiLayer;
    single_.clear();
    setup_clean_leaf();
  }

  // No resistance: that is a per-call driver on both paths, from set_physiology's
  // RootNetwork. `gravity_head`, the head to lift water to the collar in MPa, is
  // configuration, because this path has no depth profile to derive it from
  // (hazard 6); one layer of thickness d under the multi-layer rule is
  // `<gravity head> * d / 2`.
  void set_supply_single(double gravity_head = 0.0) {
    if (!std::isfinite(gravity_head) || gravity_head < 0.0) {
      util::stop("set_supply_single needs a finite, non-negative gravity_head in "
                 "MPa; got " + util::to_string(gravity_head));
    }
    supply_kind_ = SupplyKind::SinglePotential;
    setup_clean_leaf();
    // After setup_clean_leaf(): clear() spares grav_head_, which is configuration.
    single_.grav_head_ = gravity_head;
  }

  // Which path is in force, as a string: the enum has no R representation.
  std::string supply_kind_name() const {
    return supply_kind_ == SupplyKind::MultiLayer ? "multilayer" : "single";
  }

  // --- supply dispatch -------------------------------------------------------
  // Where the two paths differ; the rest of the solve reads the suction vector both
  // provide. `default:` means the single path, so a third path needs an explicit arm
  // in every switch (hazard 6).
  double supply_begin_solve() {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.begin_solve();
      default:                     return single_.begin_solve();
    }
  }
  // The current soil state, as positive suction magnitudes.
  const std::vector<double>& supply_psi_soil() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.psi_soil_;
      default:                     return single_.psi_soil_vec_;
    }
  }
  // The wet bound's floor (MultiLayerRoots::wet_floor_layer): the lowest balance
  // psi_i + z_i among conducting layers, else the wettest layer's. `conducts` and
  // `layer` say which, for the derivative surface.
  double supply_wet_floor(bool* conducts = nullptr, int* layer = nullptr) const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: {
        bool c = false;
        const int i = roots_.wet_floor_layer(roots_.psi_soil_, c);
        if (conducts != nullptr) *conducts = c;
        if (layer != nullptr) *layer = i;
        const std::size_t k = static_cast<std::size_t>(i);
        return roots_.psi_soil_[k] + roots_.grav_head_z_[k];
      }
      default:
        if (conducts != nullptr) *conducts = true;
        if (layer != nullptr) *layer = 0;
        return single_.psi_soil_ + single_.grav_head_;
    }
  }
  // Its partner (MultiLayerRoots::wet_ceiling): the highest such balance, where
  // uptake is at least zero. Equal to the floor on one layer or a single potential.
  double supply_wet_ceiling() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.wet_ceiling(roots_.psi_soil_);
      default:                     return single_.psi_soil_ + single_.grav_head_;
    }
  }
  // The one soil potential the single-layer solvers read. It checks nothing: on a
  // multi-layer profile it returns layer zero's, so test supply_is_single_layer().
  double supply_psi_soil_scalar() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.psi_soil_[0];
      default:                     return single_.psi_soil_;
    }
  }
  bool supply_is_single_layer() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.psi_soil_.size() == 1;
      default:                     return true;
    }
  }
  int supply_n_layers() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        return static_cast<int>(roots_.soil_number_of_depths_);
      default:
        return single_.n_layers();
    }
  }


  double vcmax_25;
  // Named stem_*: there are two Weibull curves, stem and root (hazard 1).
  double stem_c;
  // P50 is the scale parameter and the trait: the curve is 2^-(psi/P50)^c.
  // ⚠️ Do not add stem_b back as settable (hazard 1): two scale parameters fix the
  // shape by c = ln(ln 2)/ln(P50/b), which diverges as P50 -> b.
  double stem_P50;
  // Derived from (stem_P50, stem_c) and read-only from R.
  double stem_b;
  double psi_crit;
  double TF24_beta2;
  double jmax_25;
  double phi_J;
  double curv_fact_elec_trans; // unitless - obtained from Smith and Keenan (2020)
  double curv_fact_colim;
  // Tolerance of maximise_profit_over_collar's golden-section fallback, taken when a
  // bracket end has no usable gradient; not of the collar root-find or the stem routes.
  double GSS_tol_abs;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
  double TF24_cost_scale;

  // Tolerance on the collar for the profit-maximising root-find, fixed by the solve:
  // ~8 orders below the 1e-4 at which this package calls a difference real.
  static constexpr double collar_root_tol = 1e-12;
  // Tolerance on ci for every supply == demand root-find; solve_ci says why it
  // is 1e-13.
  static constexpr double ci_root_tol = 1e-13;
  // Tolerance to which find_root_psi places the collar interval's two bounds,
  // and that tolerance at x, absolute plus relative.
  static constexpr double bracket_tol = 1e-4;
  static double bracket_tol_at(double x) {
    return 1.5 * bracket_tol * (1.0 + std::abs(x));
  }
  // A stationary point this close to the dry end, in MPa of collar, is reported as
  // the pin at the end: the stem's inverse resolves its rows to ~1e-6 here and
  // nothing within ~1e-10. ⚠️ Measured on the collar, not on the stem's remaining
  // capacity, which is rounding noise there and flips under a tangent (hazard 19).
  static constexpr double dry_pin_resolution = 1e-8;
  // How far inside the dry end, in MPa, the marginal follows the log of the distance
  // to it; an interior point this close takes its curvature in that logarithm, which
  // is why the solve and a replay solve the end to collar precision.
  static constexpr double dry_end_neighbourhood = 1e-2;


  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
  std::vector<double> soil_consumption_;
  double profit_;
  // ⚠️ Kept only because plant's RcppR6 bindings name it (hazard 7); nothing here
  // reads or writes it, and the initialiser makes a read through them a defined NA.
  double psi_stem = util::na_value;
  // The Cowan-Farquhar marginal value of water, carbon per unit transpiration. An
  // input: no model code assigns it and setup_clean_leaf() does not clear it.
  //
  // ALIGNED DEFAULTS. Each level parameter (`CF77_lambda_`, `JS22_gamma`, `CMax_a`
  // with `CMax_b` 0, `LeastCost_beta`) defaults to where its stem optimum has
  // TF24's conductance at one reference: single-potential supply, psi_soil 1 MPa,
  // PPFD 900, VPD 2 kPa, ca 40 Pa, leaf at 25 C, energy balance off, other traits at
  // default. R's `align_curves()` solves them and `test-align-curves.R` fails when a
  // change moves an optimum there: re-align every copy rather than loosen it.
  double CF77_lambda_ = 107750.8846;            // umol C (kg H2O)^-1
  // Optional Medlyn (2011) soil-moisture factor, making the CF77 price
  // `CF77_lambda_ / beta(theta)`. Off must stay bit-identical (the golden files pin
  // CF77 rows), so cf77_price() then returns CF77_lambda_ rather than dividing by 1.
  bool CF77_soil_beta_ = false;
  // TF24_floor's price of water as transpiration goes to zero, and its only
  // parameter: an input like `CF77_lambda_`, with no default. Unset is refused,
  // because at zero the curve is TF24 under another name. Its scale is the leaf's
  // lambda, 9e4-3e5 at the defaults; far outside it the optimum pins (hazard 16).
  double TF24_floor_lambda_o = util::na_value;     // umol C (kg H2O)^-1
  // Least-cost's unit-cost ratio beta = c_V / c_E (Prentice et al. 2014), the curve's
  // only parameter; an input like `CF77_lambda_`. Dimensionless, so leastcost_flux()
  // puts E and vcmax_ in the same molar units, and it weighs vcmax_, so this is
  // Prentice's beta/eta* at the leaf's temperature. The default is aligned, not the
  // published 146, which at the default leaf pins the optimum at psi_crit.
  double LeastCost_beta = 51.59667042;             // dimensionless, c_V / c_E
  double hydraulic_cost_;
  
  double electron_transport_;
  double gamma_;
  double ko_;
  double kc_;
  double km_;
  double R_d_;
  double leaf_specific_conductance_max_;
  double k_s_;
  double vcmax_;
  double jmax_;
  double lma_; //kg m^-2
  
  double leaf_temp_;

  // The leaf's own temperature at the operating point, deg C: an output. Solved from
  // transpiration on the energy-balance path, where `leaf_temp_` is air temperature;
  // off it, a copy of `leaf_temp_` rather than NA.
  double Tleaf_;

  // Energy-balance state, re-seated by every set_physiology: a caller supplying its
  // own radiation budget sets these after each call.
  double Tair_;  // air temperature, deg C (reinterprets the leaf_temp driver)
  // Isothermal net radiation, W m^-2: what a leaf at air temperature would absorb
  // net. Its departure from Tair is g_rad_'s (see the longwave block in constants.hpp).
  double Rn_;
  // Radiative conductance, W m^-2 K^-1: 4*eps_leaf*sigma*Tair_K^3, the linearised
  // longwave response to leaf temperature, beside vol_heat_cap_air/ra_.
  double g_rad_;
  // Boundary-layer resistance to sensible heat, s m^-1, convective only. A combined
  // convective-plus-radiative resistance here needs g_rad_ = 0, or radiation is
  // counted twice and Tleaf - Tair comes out 5-7% small.
  double ra_;
  // Gate for the energy balance. Off by default, with leaf_temp prescribed.
  bool use_energy_balance_ = false;
  // --- the residual (cuticular) conductance ----------------------------------
  //
  // The cuticle passes water and CO2 after the stomata shut. The leaf optimises
  // its stomatal flux and the xylem carries the total:
  //
  //     E   = E_s + r,            r  = g_min(T_leaf) D_leaf / P_atm
  //     g_c = g_cs + rho g_min,   g_cs = (E - r) P_atm / (1.67 D_leaf)
  //
  // with g_c the total conductance to CO2 and rho = g_min_co2_ratio_. Every cost
  // that reads the flux reads the total. With g_c > 0 the supply == demand residual
  // has a root for every leaf, so a shut leaf is the ordinary evaluator at E_s = 0
  // (hazard 18). A trait, default 0, where the term is an exact zero appended to
  // every expression (hazard 18). Measured cuticular conductance, not a fitted g0,
  // which is mostly stomata not yet closed (Duursma et al. 2018).
  double g_min = 0.0;  // mol H2O m^-2 s^-1, at 25 C
  // The cuticle's CO2:H2O conductance ratio, 1/20 to 1/40 measured (Boyer 2015);
  // 0.05 is the upper end.
  double g_min_co2_ratio_ = 0.05;
  // The temperature response of g_min, a softplus blend of two exponentials,
  //
  //     g_min(T) = g_min * exp( k_lo (T - 25) + w k_hi [ sp((T - T_t)/w) -
  //                                                     sp((25 - T_t)/w) ] )
  //
  // with sp(x) = log(1 + exp(x)): slope k_lo below T_t, k_lo + k_hi above, exactly 1
  // at 25 C, and no kink. The defaults are placeholders, inert at 25 C and at
  // g_min = 0; a run that leans on the hot branch needs measured values.
  double g_min_k_low_ = 0.03;     // K^-1, the slope below the transition
  double g_min_k_high_ = 0.15;    // K^-1, the EXTRA slope above it
  double g_min_T_trans_ = 40.0;   // deg C, the transition temperature
  double g_min_T_width_ = 2.5;    // K, the width over which the slopes blend
  // The soft cap on the realised residual flux, r* = r (1 + (r/E_cap)^p)^(-1/p),
  // with p = residual_cap_power and E_cap = phi kmax (G(inf) - G(psi_up)), phi =
  // residual_cap_fraction. Smooth, so psi_soil has no kink; phi < 1 keeps the shut
  // stem's potential finite (hazard 18).
  static constexpr double residual_cap_power = 8.0;
  static constexpr double residual_cap_fraction = 0.5;
  // Outputs, written by every path out of a solve (hazard 8). `transpiration_` is
  // the total; the stomatal flux is `transpiration_ - residual_flux_`.
  double g_min_leaf_ = 0.0;        // mol H2O m^-2 s^-1
  double residual_flux_ = 0.0;     // kg H2O m^-2 s^-1
  // The collar route's second cap: residual_cap_fraction of what the roots can
  // deliver (MultiLayerRoots::uptake_limit), so the wet bound always exists and
  // moves continuously with the soil. Infinite, and exactly inert, off the collar
  // route and at g_min = 0; seated by prepare_collar_solve, cleared by stem_wet_bound.
  double residual_root_ceiling_ = std::numeric_limits<double>::infinity();
  // g_min_factor's one-entry memo: its five inputs, then its value.
  mutable std::array<double, 6> g_min_factor_memo_{
      std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0, 0.0, 0.0};

  // The residual conductance's temperature factor, exactly 1 at 25 C, and its slope.
  // Memoised on all five inputs, since a solve with the energy balance off asks at
  // one temperature throughout.
  double g_min_factor(double leaf_temp) const {
    if (leaf_temp == g_min_factor_memo_[0] && g_min_k_low_ == g_min_factor_memo_[1] &&
        g_min_k_high_ == g_min_factor_memo_[2] && g_min_T_trans_ == g_min_factor_memo_[3] &&
        g_min_T_width_ == g_min_factor_memo_[4]) {
      return g_min_factor_memo_[5];
    }
    return g_min_factor_compute(leaf_temp);
  }
  double g_min_factor_compute(double leaf_temp) const;
  double g_min_factor_temp_deriv(double leaf_temp) const;
  // Zero at g_min = 0 without forming the factor.
  double g_min_at(double leaf_temp) const {
    return g_min == 0.0 ? 0.0 : g_min * g_min_factor(leaf_temp);
  }
  // The leaf-to-air deficit at a leaf temperature, as set_leaf_vpd seats it.
  double leaf_vpd_at(double leaf_temp) const;
  // The most the stem can ever carry from an upstream potential:
  // kmax (G(inf) - G(psi_up)).
  double supply_ceiling(double psi_upstream) const;
  // The soft cap above, as a function of the demand and the ceiling, and its
  // slope in the demand and in the ceiling.
  static double residual_cap(double r, double E_inf);
  static double residual_cap_ddemand(double r, double E_inf);
  static double residual_cap_dceiling(double r, double E_inf);
  // The cap at the roots' ceiling, applied to the stem-capped flux r1; the identity
  // while that ceiling is infinite. The ceiling is U (U / (U + r))^2, not U, which
  // keeps the shut leaf C1 through the roots' cap (hazard 18).
  static double residual_root_ceiling_at(double U, double r) {
    const double w = U / (U + r);
    return U * w * w;
  }
  double residual_root_cap(double r1, double r) const {
    return std::isinf(residual_root_ceiling_)
               ? r1
               : residual_cap(r1, residual_root_ceiling_at(residual_root_ceiling_, r));
  }
  // Its slope in r1, and its slope in r through the ceiling.
  double residual_root_cap_ddemand(double r1, double r) const {
    return std::isinf(residual_root_ceiling_)
               ? 1.0
               : residual_cap_ddemand(
                     r1, residual_root_ceiling_at(residual_root_ceiling_, r));
  }
  double residual_root_cap_dceiling_dr(double r1, double r) const {
    if (std::isinf(residual_root_ceiling_)) {
      return 0.0;
    }
    const double U = residual_root_ceiling_;
    const double w = U / (U + r);
    const double dC_dr = -2.0 * w * w * w;
    return residual_cap_dceiling(r1, residual_root_ceiling_at(U, r)) * dC_dr;
  }
  // The realised residual flux at an upstream potential and a total flux,
  // kg H2O m^-2 s^-1, exactly zero at g_min = 0. Pure, so a bracket can ask freely.
  // The g_min = 0 test is in the class, here and below, so it inlines (hazard 18).
  double residual_flux(double psi_upstream, double E_total) const {
    return g_min == 0.0 ? 0.0 : residual_flux_active(psi_upstream, E_total);
  }
  double residual_flux_active(double psi_upstream, double E_total) const;
  // The demand before the cap at a leaf temperature, and its slope in temperature.
  double residual_demand_at(double leaf_temp) const;
  double residual_demand_dT(double leaf_temp) const;
  // The realised flux's slope in the upstream potential, and the temperature slope of
  // the residual's share of the CO2 conductance. Both exact zeros at g_min = 0.
  double residual_flux_dpsi_upstream(double psi_upstream, double leaf_temp) const {
    return g_min == 0.0 ? 0.0
                        : residual_flux_dpsi_upstream_active(psi_upstream, leaf_temp);
  }
  double residual_flux_dpsi_upstream_active(double psi_upstream,
                                            double leaf_temp) const;
  double residual_cond_temp_deriv(double psi_upstream, double leaf_temp,
                                  double gc_const) const {
    return g_min == 0.0 ? 0.0
                        : residual_cond_temp_deriv_active(psi_upstream, leaf_temp,
                                                          gc_const);
  }
  double residual_cond_temp_deriv_active(double psi_upstream, double leaf_temp,
                                         double gc_const) const;
  // The same flux on the tape, at a collar, with rows in g_min, kmax and the stem
  // curve's pair; leaf temperature and deficit are the seated point's, passive.
  template <class S> struct SupplyDraw;
  template <class S>
  S residual_flux_at(const S& collar, const leaf_pars<S>& pars,
                     const SupplyDraw<S>& draw) const;
  // The cuticle's conductance to CO2 at the seated leaf temperature.
  double cuticular_cond_CO2() const { return g_min_co2_ratio_ * g_min_leaf_; }
  // Net assimilation as c_i -> infinity, the colimited capacity less respiration:
  // the carbon of a leaf with no conductance that cannot net zero at any c_i.
  double assim_saturated() const;
  // c_i from a total conductance to CO2: the supply == demand root, bracketed
  // above c_a where A(c_a) < 0. See the definition for the g_c = 0 limit.
  double solve_ci(double gc);
  // Boundary-layer inputs for ra = C_ra*sqrt(d/U0), read only on the energy-balance
  // path; where either is unusable, set_physiology falls back to the fixed ra.
  double d_ = 0.05;          // characteristic leaf dimension, m
  double wind_speed_ = 2.0;  // above-canopy wind speed U0, m s^-1

  // Absorptances, fractions of intercepted radiation, applied here and nowhere
  // upstream; alpha_NIR_ reaches only the energy balance. A gas-exchange fit
  // identifies only phi_J * alpha_PAR_: fix alpha_PAR_ and phi_J is per absorbed photon.
  double alpha_PAR_ = 0.85;  // PAR absorptance, unitless
  double alpha_NIR_ = 0.20;  // near-infrared absorptance, unitless
  // Intercepted near-infrared per unit leaf area, W m^-2, for the energy balance.
  // No default: it depends on the light source (about PPFD_ / umol_par_per_joule in
  // daylight, zero under an LED), and set_physiology refuses it non-finite there.
  double NIR_int_ = util::na_value;

  // --- Sicangco et al. (2026) thermal cost, default OFF ----------------------
  // An instantaneous cost of PSII damage, sigmoid in leaf temperature,
  //
  //   TC(Tleaf) = 1 / (1 + exp(-r*(Tleaf - T50))),   r = 2/(T50 - Tcrit)
  //
  // with Jmax scaled by (1 - TC), as in the paper's ProfitMaxTC. r is gsthermal's,
  // which made the figures; the paper's Eqn 8 prints (T50 - Tcrit)/2.
  bool use_thermal_cost_ = false;
  double T50_ = 50.4;    // deg C; Sicangco Table 2, heatwave treatment
  double Tcrit_ = 46.5;  // deg C; Sicangco Table 2, heatwave treatment

  // --- ATLS thermal damage and acclimation, default OFF ----------------------
  // A lasting damage model on three fractions of leaf area:
  //
  //   functional F  <=>  impaired I  -->  damaged D,      F = 1 - I - D
  //
  //   dA/dt = acclim_alpha_ * softplus(T - acclim_T_; acclim_s_) - acclim_beta_ * A
  //   dI/dt = k_imp(T) * F - (k_rec(T) + k_mat(T)) * I
  //   dD/dt = k_mat(T) * I
  //   k(T)  = k * exp[(Ea/R)(1/T_ref - 1/T)]   for each of the three rates
  //
  // with T the leaf temperature and rates per day. damage_rates() evaluates this and
  // damage_trajectory() integrates it; a caller integrating the states itself seats
  // A_, I_ and D_ before each set_physiology. Behind the gate,
  // update_temperature_dependent_params shifts the Vcmax and Jmax optima by
  // topt_offset_ + dTopt_max_ * A/(A + K_A_), scales both capacities by F, and adds
  // four maintenance costs to R_d_, which F does not scale. With the gate off
  // nothing here is read and every number is bit-identical. The states enter the
  // benefit only, and the rates no objective. The derivative surface refuses the
  // layer; the dynamics' sensitivities are damage_trajectory_gradient().
  bool use_thermal_damage_ = false;
  // Constitutive thermostability shift, deg C, of both the photosynthetic optimum
  // and the impairment's isokinetic temperature.
  double topt_offset_ = 0.0;
  // Impairment, k_imp = k_i_ * exp[(k_i_Ea_/R)(1/T_iso - 1/T)] about the shifted
  // isokinetic temperature T_iso = k_i_Tiso_ + topt_offset_ + dTimpair_max_ A/(A + K_A_):
  // k_i_ in day^-1, k_i_Ea_ in J mol^-1. Defaults are Bison & Michaletz (2026)'s
  // 177-species fit: exp(b) = 3.2e-3 s^-1 (276 day^-1), T_iso 50.5 C, mean Ea
  // 316 kJ mol^-1. Their A is per minute in the paper as first published,
  // mislabelled per second (corrected in Zenodo 10.5281/zenodo.23752859 v3).
  double k_i_ = 276.0;
  double k_i_Ea_ = 316.0e3;
  double k_i_Tiso_ = 50.5;
  // The baseline k_i_ is costed against, day^-1: protection bought is
  // k_i_baseline_ - k_i_, signed.
  double k_i_baseline_ = 276.0;
  // Repair and maturation: rates at 25 C, day^-1, each with its own Ea, J mol^-1.
  double k_rec_ = 0.25;
  double k_rec_Ea_ = 0.0;
  double k_mat_ = 0.02;
  double k_mat_Ea_ = 0.0;
  // Acclimation kinetics, above: gain day^-1 deg C^-1, decay day^-1, threshold
  // deg C, sharpness deg C^-1.
  double acclim_alpha_ = 0.02;
  double acclim_beta_ = 0.01;
  double acclim_T_ = 30.0;
  double acclim_s_ = 1.0;
  // What acclimation buys, deg C: the most each optimum can rise, half-saturating
  // in A at K_A_.
  double dTopt_max_ = 6.0;
  double dTimpair_max_ = 6.0;
  double K_A_ = 1.0;
  // State, seated by the caller before each set_physiology; zero is unacclimated
  // and intact. set_physiology refuses I_ < 0, D_ < 0 and A_ <= -K_A_.
  double A_ = 0.0;
  double I_ = 0.0;
  double D_ = 0.0;
  // Maintenance costs added to R_d_, umol C m^-2 s^-1 per unit of A_, of k_rec_, of
  // protection bought, and of the repair flux k_rec_eff * I_. Zero by default.
  double c_acclim_maint_ = 0.0;
  double c_repair_maint_ = 0.0;
  double c_protect_maint_ = 0.0;
  double c_repair_flux_ = 0.0;
  // Width of F's soft floor; positive, so the capacities stay positive as I_ + D_ -> 1.
  static constexpr double damage_surv_floor = 1e-6;
  // Outputs, at the temperature the block was last seated at.
  double k_imp_eff_ = 0.0;   // impairment rate there, day^-1; zero with the gate off
  double k_rec_eff_ = 0.0;   // repair rate there, day^-1; zero with the gate off
  double functional_ = 1.0;  // the F that scaled the capacities; one with the gate off

  // Cells a stem route scans for the basin, where basin_scan_cells() asks for a
  // scan. 64 is where the answer stops moving: against a 2001-point reference over
  // a 1728-row sweep, 64 matches every row and 32 leaves 6 rows short by 3.9e-04.
  // GSS_tol_abs does not sharpen a stem argmax: the refine is a root-find.
  int boundary_scan_n_ = 64;

  // Zero unless the energy balance is on: with it off every stem objective in that
  // sweep is unimodal. ⚠️ Re-measure before widening this: a missed basin is a
  // large error (hazard 11).
  int basin_scan_cells() const {
    return use_energy_balance_ ? boundary_scan_n_ : 0;
  }
  // Intercepted PPFD per unit leaf area, umol m^-2 s^-1: incident, not absorbed.
  double PPFD_;
  double atm_vpd_;
  // The deficit the diffusion equations use, kPa: `atm_vpd_` off the energy
  // balance, the leaf-to-air deficit at the operating leaf temperature on it.
  double vpd_leaf_;
  double atm_o2_kpa_;

  // --- Temperature-response parameters ---------------------------------------
  // plantecophys's EaV / EdVC / delsC / EaJ / EdVJ / delsJ, settable because thermal
  // acclimation modifies them; defaults in constants.hpp.
  double vcmax_ha_ = phylloptim::vcmax_ha;    // activation energy, J mol^-1
  double vcmax_H_d_ = phylloptim::vcmax_H_d;  // deactivation energy, J mol^-1
  double vcmax_d_S_ = phylloptim::vcmax_d_S;  // entropy term, J mol^-1 K^-1
  double jmax_ha_ = phylloptim::jmax_ha;
  double jmax_H_d_ = phylloptim::jmax_H_d;
  double jmax_d_S_ = phylloptim::jmax_d_S;
  // Bernacchi kinetics, value at 25 C and activation energy.
  double gamma_25_ = phylloptim::gamma_25;    // CO2 compensation point, umol mol^-1
  double gamma_ha_ = phylloptim::gamma_ha;
  double kc_25_ = phylloptim::kc_25;          // Rubisco Km for CO2, umol mol^-1
  double kc_ha_ = phylloptim::kc_ha;
  double ko_25_ = phylloptim::ko_25;          // Rubisco Km for O2, umol mol^-1
  double ko_ha_ = phylloptim::ko_ha;
  // Dark respiration at 25 C, umol m^-2 s^-1, a trait in its own right. Initialised
  // here because plant's RcppR6 bindings pin the constructor by arity.
  double R_d_25 = 1.44;
  // Joshi & Stocker (2022)'s hydraulic unit cost, `gamma` in
  // `F = A - alpha*Jmax - gamma*(dpsi)^2`, initialised here for R_d_25's reason.
  // Their capacity term is not implemented, so this is their cost, not their model.
  // The default is aligned (see `CF77_lambda_`).
  double JS22_gamma = 0.4717887793;    // umol C m^-2 s^-1 MPa^-2
  // Wolf, Anderegg & Pacala (2016) carbon maximisation as Sabot's TractLSM
  // implements it (`SPAC/fregulate.py`), a marginal cost linear in leaf potential,
  //
  //     dC/dpsi = CMax_a*psi + CMax_b
  //
  // CMax_b is signed, negative in the paper's convention, and TractLSM stores it
  // sign-inverted: pin the sign by reproducing a TractLSM number. CMax_a is aligned
  // with CMax_b = 0 (see `CF77_lambda_`).
  double CMax_a = 0.6385596982;        // umol C m^-2 s^-1 MPa^-2
  double CMax_b = 0.0;                 // umol C m^-2 s^-1 MPa^-1
  // Dark respiration's temperature response, Tjoelker et al. (2001):
  //
  //     R_d(T) = R_d_25 * Q10(T)^((T - 25) / 10)
  //     Q10(T) = rd_q10_intercept_ - rd_q10_slope_ * (T + 25) / 2
  //
  // A declining Q10: a constant 2 leaves the solve no operating point by 45 C.
  double rd_q10_intercept_ = 3.09;
  double rd_q10_slope_ = 0.0430;
  double atm_kpa_;
  // Mixing ratio (umol mol^-1) to partial pressure (Pa), 1e-6 * P: derived from
  // atm_kpa_, so every CO2 quantity sees the pressure the stomatal side does.
  double umol_per_mol_to_Pa_;
  double ca_;
  double opt_root_psi_;
  double assim_max_;

  // --- Sperry (2017) ProfitMax outputs, all unitless -------------------------
  // Written only by the ProfitMax path; the normalisers are seeded by
  // prepare_profitmax() and constant over one solve.
  double profitmax_A_max_;   // |A|max over the supply stream, umol m^-2 s^-1
  double profitmax_k_soil_;  // k(psi_soil), kg m^-2 s^-1 MPa^-1
  double profitmax_k_span_;  // k(psi_soil) - k_crit, the HC denominator
  // The normaliser ratio profitmax_A_max_ / profitmax_k_span_ is per unit
  // conductance lost, not a lambda. The lambda every curve reports on one axis,
  // with the `CF77_lambda_` input, is
  //
  //     lambda_emergent = (dC/dpsi) / (dE/dpsi)      umol C (kg H2O)^-1
  //
  // written by whichever optimiser ran; NA before one has.
  double lambda_emergent_ = util::na_value;
  double carbon_gain_;         // CG = A/|A|max over the supply stream
  double hydraulic_cost_norm_; // HC = [k(psi_soil)-k(psi)]/[k(psi_soil)-k_crit]
  double thermal_cost_;        // TC, zero unless use_thermal_cost_

  
  double opt_psi_stem_;
  // Kept for plant's bindings (hazard 7); the operating point's ci is ci_.
  double opt_ci_;
  // ⚠️ Kept only for plant's bindings, like psi_stem (hazard 7).
  double count = util::na_value;
  double E_up_;

  // --- Medlyn stomatal-conductance model -------------------------------------
  // A standalone, R-callable alternative to the profit optimisation
  // (solve_medlyn_ci_*), not used by the TF24 path.
  double g0 = 0.022;        // residual stomatal conductance (umol m^-2 s^-1)
  double g1 = 2.57;         // sensitivity to vpd (kPa^0.5)

  // The H2O:CO2 stomatal diffusion ratio, a convention (Medlyn et al. 2011 use 1.6).
  // It reaches the solve, so re-solve after changing it.
  double H2O_CO2_stom_diff_ratio_ = phylloptim::H2O_CO2_stom_diff_ratio;

  // The Medlyn route's ratio, 1.6 because published g1 is fitted with it. The USO
  // expression, the residual and the CO2 conversion must all use it, or
  // A = gs (ca - ci) fails by 4.2%.
  static constexpr double medlyn_H2O_CO2_ratio = 1.6;

  // CO2 basis: medlyn_model_gs divides the USO water-vapour conductance by
  // medlyn_H2O_CO2_ratio before storing.
  double medlyn_model_gs_;  // mol CO2 m^-2 s^-1
  double theta_w_;          // current soil water content at wilting point (m^3 m^-3)
  double theta_fc_;         // current soil water content at field capacity (m^3 m^-3)
  double theta_;            // current soil water content (m^3 m^-3)

  // One-entry memo for transpiration(), keyed on (psi_stem, psi_upstream): the rest
  // is fixed within a solve. Invalidated in set_physiology.
  bool   transpiration_cached_ = false;
  double transpiration_cache_psi_stem_ = 0.0;
  double transpiration_cache_psi_upstream_ = 0.0;
  double transpiration_cache_value_ = 0.0;

  // Cache of the temperature- and O2-dependent block seated in set_physiology,
  // reused while photo_temp_key() is unchanged. electron_transport_ reads PPFD_
  // and is not cached.
  bool   photo_temp_cached_ = false;
  // The leaf temperature the block holds. Kernels reading the pack must read this,
  // not `leaf_temp_` (air temperature on the energy-balance path); test_leaf's
  // member-against-pack comparison catches the difference.
  double photo_temp_;
  // The key is every input of update_temperature_dependent_params(), traits
  // included (hazard 10). ⚠️ A member added to that block and left out of
  // photo_temp_key() is served stale, silently.
  static constexpr int photo_temp_key_size = 23;
  std::array<double, photo_temp_key_size> photo_temp_cache_key_{};
  std::array<double, photo_temp_key_size> photo_temp_key() const;
  // ⚠️ Kept only for plant's bindings, like psi_stem (hazard 7).
  std::vector<double> f_r;

  double atm_vpd = 2.0; //kPa
  double ca = 40.0; // Pa
  double atm_kpa = 101.3; //kPa
  //partial pressure o2 (kPa)
  double atm_o2_kpa = 21;
  //leaf temperature (deg C)
  double leaf_temp = 25;
  // Default soil-moisture content for the Medlyn beta_ factor.
  double theta_w = 0.2;  //m^3 m^-3
  double theta_fc = 0.5; //m^3 m^-3
  double theta = 0.3;    //m^3 m^-3


  
  // `root_network`: per-layer root hydraulic resistances; deriving them from root
  // carbon is the caller's business (root_network_from_carbon). ⚠️ Per unit leaf
  // area (hazard 4): resistances from absolute root carbon give a silently wrong
  // E_up. On the energy-balance path, seat NIR_int_ before this call.
  void set_physiology(const RootNetwork& root_network, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double leaf_temp, double atm_o2_kpa, double atm_kpa);

  // Replace the traits, in `leaf_traits()`' order, leaving the numerical controls
  // alone: the one way to change a trait (hazard 10). Rebuilds a table only when
  // its c moved and ends with setup_clean_leaf(), so call set_physiology() again.
  void set_traits(double vcmax_25, double stem_c, double stem_P50,
                  double root_c, double root_P50,
                  double TF24_beta2, double jmax_25, double phi_J,
                  double curv_fact_elec_trans, double curv_fact_colim,
                  double TF24_cost_scale, double R_d_25, double JS22_gamma,
                  double CMax_a, double CMax_b, double g_min);

  // The critical point is a quantile of the stem curve: the conductivity fraction
  // left there, 0.05 giving P95, also Sperry's k_crit = f * kmax. Moving it is not
  // local: it sets SOX's and JW26's costs, ProfitMax's normaliser, the P50 chain
  // rule (hazard 14) and the dry end of the curves cost_domain_ends_at_psi_crit names.
  static constexpr double k_crit_fraction = 0.05;
  // Knot intervals in each vulnerability table, over x = (psi/b)^c; a consumer seats
  // the same number. 800 puts the table within 1.1e-14 of G(inf) for c from 1 to 12,
  // two orders below the solve's floor (psi_stem_to_ci's 1e-13).
  static constexpr double ncontrol_default = 800.0;

  // The two derivations, one place each; `f` is the remaining-conductivity fraction.
  static double weibull_b_from_P50(double P50, double c) {
    return phylloptim::weibull_b_from_P50(P50, c);
  }
  static double weibull_psi_at_fraction(double b, double c, double f) {
    return b * std::pow(std::log(1.0 / f), 1.0 / c);
  }

  // Seats both curves' derived scale and critical potential. Constructors only:
  // set_traits computes these before assigning, to decide which tables to rebuild.
  void derive_vulnerability_curves_from_P50() {
    stem_b = weibull_b_from_P50(stem_P50, stem_c);
    psi_crit = weibull_psi_at_fraction(stem_b, stem_c, k_crit_fraction);
    roots_.root_b = weibull_b_from_P50(roots_.root_P50, roots_.root_c);
    roots_.root_psi_crit =
        weibull_psi_at_fraction(roots_.root_b, roots_.root_c, k_crit_fraction);
  }

  // The positive-magnitude checks (hazard 2), for the constructor and set_traits.
  static void check_psi_magnitudes(double stem_P50, double stem_c,
                                   double root_P50, double root_c);

  void setup_transpiration(double resolution);

  // The four operations on the stem's G; every read goes through these. G and G'
  // come from one polynomial, G^-1 is Newton on that pair, (G^-1)' is 1/G' there.
  // Past the curve G is G(inf); only a negative potential or a flux past G(inf) is
  // refused, naming `caller`.
  double stem_curve_integral(double psi, const char* caller = nullptr) const;
  double stem_curve_integral_deriv(double psi) const;
  // `converge` has no effect: the inverse always converges to a few ulps.
  double stem_curve_integral_inverse(double w, const char* caller = nullptr,
                                     bool converge = false) const;
  double stem_curve_integral_inverse_deriv(double w) const;
  // G(inf) = stem_b Gamma(1/c) / c, the most the stem carries per unit conductance.
  double stem_curve_limit() const;
  // Q = G(inf) - G(psi), to its own relative precision, and its inverse; zero past
  // the end of the curve.
  double stem_curve_complement(double psi, const char* caller = nullptr) const;
  double stem_curve_complement_inverse(double Q) const;
  // Whether psi is in the stem table's tail (VulnerabilityTable::tail_share), where a
  // flux is a difference of Q: a difference of G loses the digits (hazard 1).
  bool stem_curve_in_tail(double psi) const;
  // Whether the stem can carry E from psi_upstream, to a few ulps, and the stem
  // potential that carries it.
  bool stem_carries(double E, double psi_upstream) const;
  double psi_stem_carrying(double E, double psi_upstream,
                           const char* caller) const;

  // Move stem_P50, stem_b and psi_crit without touching the stem table. Derivative
  // work only: the solved point is not cleared, so write the outputs before reading
  // them (hazard 8); set_traits() is the way back.
  void perturb_stem_P50(double stem_P50_new);
  // Kept on Leaf for the construction sequence and plant's bindings.
  void setup_root_vulnerability(double resolution) {
    // Recorded because set_traits rebuilds both curves at this resolution.
    vulnerability_curve_ncontrol = resolution;
    roots_.setup_vulnerability(resolution, stem_curve_);
  }
  void setup_clean_leaf();

  // Tolerance of transpiration_full_integration's direct quadrature; off the
  // production path.
  double integration_tol_ = 1e-8;

  // Kept for plant's R bindings. Sets transpiration_full_integration's tolerance;
  // `integration_rule` is ignored.
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-8) {
    static_cast<void>(integration_rule);
    integration_tol_ = integration_tol;
  }

  // Medlyn stomatal conductance; R-callable, standalone.
  double medlyn_model_gs(double assim_colimited_);
  double medlyn_stom_cond_minus_coupled_stom_cond(double x);
  void solve_medlyn_ci_numerical();
  void solve_medlyn_ci_analytical();
  // Collapses a negative-conductance result to the class's closed state. Called by
  // both solve routes, so neither can report a negative gs. See the definition.
  void set_medlyn_shutdown_if_closed();
  // The soil-moisture factor both Medlyn routes multiply g1 by, clamped to
  // [0, 1]. See the definition.
  double medlyn_beta() const;

  // Refuses a psi_soil argument that is not positive magnitudes in MPa (hazard 2).
  // This object's own vector, already validated, is recognised by address and
  // skips the scan, which would otherwise sit on the hot path (hazard 5).
  void require_suction_vector(const std::vector<double>& psi_soil,
                              const char* who) const {
    if (&psi_soil == &supply_psi_soil()) {
      return;
    }
    for (size_t i = 0; i < psi_soil.size(); ++i) {
      if (!(psi_soil[i] >= 0.0)) {
        util::stop(std::string(who) +
                   ": psi_soil must be positive magnitudes in MPa, not signed "
                   "potentials (#25); got psi_soil[" + std::to_string(i) + "]=" +
                   util::to_string(psi_soil[i]));
      }
    }
  }

  // Uptake at a collar suction against any vector of layer suctions. The buffers
  // stay on Leaf because plant writes into them by name.
  void E_from_Soil_to_Root_Collar(double T_collar, const std::vector<double>& psi_soil) {
    require_suction_vector(psi_soil, "E_from_Soil_to_Root_Collar");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        roots_.uptake_at(T_collar, psi_soil, soil_consumption_, E_up_);
        break;
      default:
        single_.uptake_at(T_collar, psi_soil, soil_consumption_, E_up_);
        break;
    }
  }

  // What the leaf hands over, at any scalar: the objective plant bills, the water
  // drawn from each layer, and the cuticle's share, each with rows matching its
  // own value. The stomatal flux is the total less `residual_flux`.
  template <class T>
  struct BilledOutputs {
    T profit{};
    std::vector<T> uptake;
    T residual_flux{};
  };

  // The soil's draw at one collar, with that collar, so check_draw can refuse a
  // flux and a slope from different ones. `E_up` and its slope are an active pair;
  // the per-layer slopes are supplied at double (a live collar would record the
  // supply twice), so the per-layer draws have implicit_value's shape.
  template <class S>
  struct SupplyDraw {
    double at = 0.0;
    // Zero rather than default-initialised, so no path leaves it indeterminate.
    odelia::with_slope<S> E_up{};
    std::vector<S> uptake;
    std::vector<double> duptake_dp;
    // What the roots can deliver, capping the residual flux on the collar route;
    // carried only with a residual conductance.
    bool has_uptake_limit = false;
    S uptake_limit{};
    // Where no layer conducts, the wet bound is the wettest layer's gravity balance
    // and takes its rows: the residual is zero at every collar there.
    bool no_conduction = false;
    S gravity_balance{};

    // Everything that carries a row. ⚠️ A member not visited here gives exact-zero
    // columns, not an error (hazard 17); the static_assert below counts members.
    template <class F>
    void for_each_active(F&& f) const {
      f(E_up.value);
      f(E_up.slope);
      for (const S& u : uptake) {
        f(u);
      }
      if (has_uptake_limit) {
        f(uptake_limit);
      }
      if (no_conduction) {
        f(gravity_balance);
      }
    }
  };
  // Each optional row is a flag and a value; at double the flag pads to a double.
  static_assert(sizeof(SupplyDraw<double>) ==
                    sizeof(double) + sizeof(odelia::with_slope<double>) +
                        sizeof(std::vector<double>) * 2 +
                        2 * (alignof(double) + sizeof(double)),
                "SupplyDraw gained or lost a member: update for_each_active, "
                "then this count");
  static_assert(sizeof(odelia::with_slope<double>) == 2 * sizeof(double),
                "odelia::with_slope is no longer a value and a slope: re-check "
                "that SupplyDraw::for_each_active visits all of it");

  // The soil's draw at `collar`, which is passive whatever the caller passes: where
  // the point moves is the supplied row's business, and a live collar would record
  // the search that placed it.
  template <class S>
  SupplyDraw<S> supply_draw_at(const S& collar,
                               const SupplyInputs<S>& supply) const {
    // The walk is MultiLayerRoots', which a single-potential solve never read.
    if (supply_kind_ != SupplyKind::MultiLayer) {
      util::stop("Leaf::supply_draw_at: the derivative surface reads the "
                 "multi-layer supply, and this leaf is on the " +
                 supply_kind_name() + " path");
    }
    const std::size_t n = static_cast<std::size_t>(supply_n_layers());
    SupplyDraw<S> d;
    d.at = odelia::util::to_passive(collar);

    // A shut leaf sits on its wet bound and is drawn like any other point.

    // One walk for all three; a live collar would record the supply a second time.
    const S held = S(d.at);
    d.uptake.assign(n, S(0.0));
    CollarConductance<S> conductance;
    roots_.template uptake_at<S>(held, supply, d.uptake, d.E_up.value,
                                 &conductance);
    d.E_up.slope = conductance.total;
    d.duptake_dp = std::move(conductance.per_layer);
    if (!std::isinf(residual_root_ceiling_)) {
      d.uptake_limit = roots_.template uptake_limit<S>(supply);
      d.has_uptake_limit = true;
    }
    bool conducts = true;
    const int floor_layer = roots_.wet_floor_layer(roots_.psi_soil_, conducts);
    if (!conducts) {
      const std::size_t k = static_cast<std::size_t>(floor_layer);
      d.no_conduction = true;
      d.gravity_balance = supply.psi_soil[k] + S(roots_.grav_head_z_[k]);
    }

    // A kink is answered here, once: `duptake_dpsi` returns NaN where its analytic
    // slope is not valid (hazard 6), as at a wet bound, and then every slope is
    // differenced from its own function. A NaN left in one layer reaches plant as a
    // NaN in every trait column sharing the adjoint accumulator.
    const bool at_a_kink =
        !std::isfinite(odelia::util::to_passive(d.E_up.slope)) ||
        std::any_of(d.duptake_dp.begin(), d.duptake_dp.end(),
                    [](double x) { return !std::isfinite(x); });
    if (at_a_kink) {
      const double step = 1e-6;
      std::vector<double> soil_at;
      soil_at.reserve(supply.psi_soil.size());
      for (const S& v : supply.psi_soil) {
        soil_at.push_back(odelia::util::to_passive(v));
      }
      std::vector<double> up_layer(d.duptake_dp.size(), 0.0);
      std::vector<double> down_layer(d.duptake_dp.size(), 0.0);
      double up = 0.0, down = 0.0;
      roots_.uptake_at(d.at + step, soil_at, up_layer, up);
      roots_.uptake_at(d.at - step, soil_at, down_layer, down);
      d.E_up.slope = S((up - down) / (2.0 * step));
      for (std::size_t j = 0; j < d.duptake_dp.size(); ++j) {
        d.duptake_dp[j] = (up_layer[j] - down_layer[j]) / (2.0 * step);
      }
    }
    return d;
  }

  // A draw from the wrong collar gives a wrong number, not a missing row, so this
  // refuses it, exactly: the draw stores the collar it was taken at.
  template <class S>
  void check_draw(double expect, const SupplyDraw<S>& draw) const {
    if (draw.at != expect) {
      util::stop("Leaf: the supply draw was taken at a collar of " +
                 util::to_string(draw.at) + " and is being read at " +
                 util::to_string(expect));
    }
  }

  // The surface's coordinates: each a value and its response in the collar, which
  // travel together because every coordinate must move with the collar.

  // What profit needs: the stem potential the collar determines, and its ci.
  template <class T>
  struct CollarCoords {
    odelia::with_slope<T> sigma;
    odelia::with_slope<T> ci;
  };

  // The cost curves. A curve enters the derivative only through dC/dpsi_stem.
  enum class CostCurve { TF24, CF77, JS22, CMax, SOX, JW26, ProfitMax, TF24_floor,
                         LeastCost };

  // Every curve maximises `h(A(psi)) - C(psi)`, a benefit link h composed with
  // assimilation minus a cost, so one derivative serves all and only h' varies:
  //
  //     d/dpsi [ h(A) - C ] = h'(A) * dA/dpsi - dC/dpsi
  //
  //     Identity  h = A          h' = 1         TF24, CF77, JS22, CMax, TF24_floor
  //     Log       h = log A      h' = 1/A       SOX, JW26, LeastCost
  //     Scaled    h = A/|A|max   h' = 1/|A|max  ProfitMax
  //
  // Log serves the products (log A + log g has A g's argmax), whose cost is stored
  // as -log g, and LeastCost's ratio. The energy-balance path through Tleaf(E) is
  // part of dA/dpsi and takes the same h'; cost-side Tleaf terms take none (see
  // `dprofit_energy_balance_term`). `Scaled` holds |A|max fixed, so ProfitMax's trait
  // gradients are partials: |A|max comes from a scan.
  enum class BenefitLink { Identity, Log, Scaled };

  // Names the enum's last member: move it when a curve is appended, or the new
  // curve falls outside every bounds check.
  static constexpr int n_cost_curves = static_cast<int>(CostCurve::LeastCost) + 1;

  // The model is configuration, because its parameters already are: set_model()
  // seats it and optimise() takes nothing. TF24 on the collar route by default.
  CostCurve cost_curve_ = CostCurve::TF24;
  bool route_is_collar_ = true;
  // The method: exact (the root-find of the first-order condition) or the closed
  // form in closed_form.hpp. A method, not a model: the two must agree where the
  // closed form is valid, which `test-closed-form.R` measures.
  bool method_is_closed_ = false;

  // Whether the last optimise() fell back from the closed form, and the counts that
  // measure the fallback fraction phi: the guard tests an output, so phi cannot be
  // predicted. Written by every optimise(), exact ones included (hazard 8).
  bool last_solve_fell_back_ = false;
  int closed_form_calls_ = 0;
  int closed_form_fallbacks_ = 0;

  void find_root_collar_psi();
  // The same solve for any cost curve; find_root_collar_psi() is the TF24 one.
  template <CostCurve K> void find_root_collar_psi_for();
  // The entry point: curve, route and method by name, never a curve index. `route`
  // is "collar" (soil to collar, any supply) or "stem" (upstream pinned at psi_soil,
  // one soil potential); `method` is "exact" (the default) or "closed".
  void set_model(const std::string& curve, const std::string& route,
                 const std::string& method = "exact");
  // The same, pre-parsed. The method defaults to exact, so seating resets it, and
  // `gradient::route_seat` always differentiates the exact solve.
  void set_model(CostCurve curve, bool collar_route, bool closed_method = false);
  // The curve, route and method currently seated, by name.
  std::string model_curve() const { return curve_name(static_cast<int>(cost_curve_)); }
  std::string model_route() const { return route_is_collar_ ? "collar" : "stem"; }
  std::string model_method() const { return method_is_closed_ ? "closed" : "exact"; }
  // Solve for the operating point under the seated model. THE entry point.
  void optimise();
  // The closed-form arm of optimise(), with its refusals and its fallback. Defined
  // in closed_form.hpp, which <phylloptim.hpp> includes after this header; including
  // only leaf_model.hpp and calling it is a link error.
  void optimise_closed();
  // Fraction of closed-form optimise() calls that fell back since the counters were
  // reset; NaN before any, so "never asked" is not "never fell back".
  double closed_form_fallback_fraction() const {
    return closed_form_calls_ > 0
               ? static_cast<double>(closed_form_fallbacks_) /
                     static_cast<double>(closed_form_calls_)
               : util::na_value;
  }
  void reset_closed_form_counters() {
    closed_form_calls_ = 0;
    closed_form_fallbacks_ = 0;
  }
  // Shared setup for the root-collar solve: builds the soil-side caches and settles
  // every feasibility exit itself, returning false; otherwise true, with the feasible
  // collar interval [bound_a, bound_b]. `first_bracket` lets a curve not bounded at
  // psi_crit try the bracket that ends at psi_crit first; maximise_profit_over_collar
  // extends it. A prescribed collar needs the whole domain, so does not ask.
  template <CostCurve K>
  bool prepare_collar_solve(double& bound_a, double& bound_b,
                            bool first_bracket = false);
  // Move the bracket's dry end from psi_crit to the end of the stem curve.
  double extend_collar_bracket(double bound_a);
  // At a dry pin of a curve psi_crit does not bound, seat the stem at the end of
  // its curve, which carries the uptake, rather than where the inverse lands
  // (anywhere G rounds to G(inf)). Returns whether it did.
  template <CostCurve K> bool seat_stem_at_dry_pin();
  // The leaf at a prescribed collar potential (MPa, positive magnitude) instead
  // of the optimum: refuses a target outside the feasible interval (code
  // `prescribed_potential`), writes find_root_collar_psi's outputs, returns profit_.
  double evaluate_root_collar_psi(double target_opt_root_psi);
  template <CostCurve K>
  double evaluate_root_collar_psi_for(double target_opt_root_psi);
  // evaluate_root_collar_psi once prepare_collar_solve has run this step, so
  // evaluations share one prepare; refuses a target outside [bound_a, bound_b].
  template <CostCurve K>
  double profit_at_collar_psi(double target_opt_root_psi,
                              double bound_a, double bound_b);
  // Exact d(profit)/d(collar) at a given collar (MPa, positive magnitude), by
  // forward-mode AD and the implicit-function theorem; seats its own caches.
  // The 0.0 on the shut and reversed-gradient exits is a sentinel, not a
  // stationary point, and sits at the wet bracket end a root-find tries first;
  // `feasible`, when non-null, is false there. Not NaN: plant integrates the
  // return as a rate where 0 means "do not acclimate", and NaN means "use finite
  // differences" on duptake_dpsi (hazard 6).
  double dprofit_droot_collar_psi(double opt_root_psi, bool* feasible = nullptr);
  template <CostCurve K>
  double dprofit_droot_collar_psi_for(double opt_root_psi,
                                      bool* feasible = nullptr);
  // The same, returning `{dprofit, feasible}`, because R cannot reach a `bool*`:
  // RcppR6 drops the argument, and an R caller would lose the sentinel's flag.
  std::vector<double> dprofit_droot_collar_psi_checked(double opt_root_psi);
  // dprofit_droot_collar_psi with the supply path's per-solve caches already
  // seated, so the collar solve seats them once for all its evaluations.
  template <CostCurve K>
  double dprofit_at_collar_psi(double opt_root_psi, bool* feasible = nullptr);


  // dprofit/dpsi_stem with the upstream potential held, for the stem routes:
  // dprofit_at_collar_psi with dpsi_stem/dx = 1. `feasible` false means the 0.0
  // returned is a sentinel, not a stationary point.
  template <CostCurve K>
  double dprofit_dpsi_stem(double psi_stem, double psi_upstream,
                           bool* feasible = nullptr);

  static std::string curve_name(int curve);

  // ⚠️ Hazard 17: the derivative surface read on another curve is a wrong,
  // finite number. Checked against the SEATED curve; find_root_collar_psi()
  // solves TF24 whatever is seated, so a leaf solved through it on another curve
  // passes this check wrongly.
  template <CostCurve K>
  void check_surface_curve(const char* who) const {
    if (K != cost_curve_) {
      util::stop(std::string(who) + ": asked for the " +
                 curve_name(static_cast<int>(K)) + " curve, and the leaf is "
                 "seated on " + curve_name(static_cast<int>(cost_curve_)) +
                 "; read the surface on the curve the solve used");
    }
  }

  // The only runtime curve dispatch: calls `f(tag)` with `tag` an
  // `std::integral_constant<CostCurve, K>`, so the body is compiled per curve.
  // A new curve needs an arm here and in curve_name(); neither switch has a
  // `default:`, so `-Werror=switch` makes a forgotten curve a build failure.
  template <typename F>
  static decltype(auto) with_curve(CostCurve curve, F&& f) {
    switch (curve) {
      case CostCurve::TF24:
        return f(std::integral_constant<CostCurve, CostCurve::TF24>{});
      case CostCurve::CF77:
        return f(std::integral_constant<CostCurve, CostCurve::CF77>{});
      case CostCurve::JS22:
        return f(std::integral_constant<CostCurve, CostCurve::JS22>{});
      case CostCurve::CMax:
        return f(std::integral_constant<CostCurve, CostCurve::CMax>{});
      case CostCurve::SOX:
        return f(std::integral_constant<CostCurve, CostCurve::SOX>{});
      case CostCurve::JW26:
        return f(std::integral_constant<CostCurve, CostCurve::JW26>{});
      case CostCurve::ProfitMax:
        return f(std::integral_constant<CostCurve, CostCurve::ProfitMax>{});
      case CostCurve::TF24_floor:
        return f(std::integral_constant<CostCurve, CostCurve::TF24_floor>{});
      case CostCurve::LeastCost:
        return f(std::integral_constant<CostCurve, CostCurve::LeastCost>{});
    }
    // Unreachable, and not a `default:` arm, which would satisfy -Werror=switch.
    util::stop("unknown cost curve index " +
               util::to_string(static_cast<int>(curve)));
    return f(std::integral_constant<CostCurve, CostCurve::TF24>{});
  }
  // The leaf at a prescribed psi_stem under the curve set_model() seated; returns
  // profit_. No curve argument: the curve's constants live on the object.
  double evaluate_psi_stem_at(double target_psi_stem);
  // Refuse a stem-route method on a leaf seated on the collar route. See the
  // definition for why the route and not only the curve has to be checked.
  void require_stem_route(const char* who) const;
  // Returns {dprofit, feasible}. No `psi_upstream` argument: on a stem route it
  // is always psi_soil, read through the accessor optimise_psi_stem_single uses
  // (R's `psi_soil_` is the multi-layer vector, empty on the single path). The
  // curve is the seated one, as for evaluate_psi_stem_at.
  std::vector<double> dprofit_dpsi_stem_checked(double psi_stem);

  template <CostCurve K> static constexpr BenefitLink benefit_link();
  template <CostCurve K> double benefit_link_deriv(double A) const;
  template <CostCurve K> void optimise_psi_stem_single();
  // ProfitMax's shut state: with |A|max <= 0 its objective A/|A|max is undefined,
  // so it is zeroed, and the leaf's own state is the evaluator's at the wet end,
  // written rather than left stale (hazard 8).
  void set_profitmax_degenerate_state(double psi_stem, double psi_upstream);
  template <CostCurve K> void check_cost_parameters();
  template <CostCurve K> double cost_deriv(double psi_stem,
                                           double psi_upstream);
  template <CostCurve K> static constexpr bool cost_reads_upstream();
  template <CostCurve K> static constexpr bool cost_domain_ends_at_psi_crit();
  // Whether the objective is A times a gain factor, which has no maximum
  // where A <= 0 throughout, so the leaf is shut there (dark_under_product).
  template <CostCurve K> static constexpr bool objective_is_product();
  // For such a curve, whether A <= 0 with the stem at its dry end, the most
  // open the leaf can be; this leaves the leaf evaluated there.
  template <CostCurve K> bool dark_under_product(double psi_dry,
                                                 double psi_upstream);
  // The stem potential past which a solve does not look: psi_crit for a curve
  // whose own domain ends there, and the end of the stem curve for every other.
  template <CostCurve K> double dry_end() const;
  template <CostCurve K> double cost_deriv_upstream(double psi_stem,
                                                    double psi_upstream);
  template <CostCurve K> double profit_psi_stem_for(double psi_stem,
                                                    double psi_upstream);
  template <CostCurve K> double lambda_for(double psi_stem,
                                           double psi_upstream);
  // The collar route's lambda_emergent, for the seated curve. See the definition.
  template <CostCurve K> double lambda_collar_for();

  template <CostCurve K>
  double evaluate_psi_stem(double target_psi_stem);
  // d/dTleaf of what `read` takes from the temperature block, by a central
  // difference of the block. It restores every member the block writes; a
  // member added to update_temperature_dependent_params must be added here and
  // to dprofit_energy_balance_term, or the leaf stays 1e-3 degrees off.
  template <class F> double temp_deriv(double leaf_temp, F read);
  // dVcmax/dTleaf, the same way and with the same restore, for lambda_LeastCost
  // with the energy balance on (LeastCost's capacity cost reads vcmax_).
  double vcmax_temp_deriv(double leaf_temp);
  // The energy-balance contribution to dJ/dpsi through Tleaf(E(psi)), zero with
  // the gate off; out of line so it cannot change FMA contraction on the gate-off
  // path. `E_s`, `r_star`, `drstar_dpsi` and `dgc_dT_residual` are the residual
  // conductance's terms, exact zeros at g_min = 0; derivation at the definition.
  template <CostCurve K>
  double dprofit_energy_balance_term(double ci, double gc, double g_ci,
                                     double inv_atm, double gc_const,
                                     double A_prime, double dgc_dpsistem,
                                     double dgc_dpsi, double dpsistem_dpsi,
                                     double dT_dE, double Tleaf, double E_s,
                                     double r_star, double drstar_dpsi,
                                     double dgc_dT_residual, double bl_slope);
  // The profit-maximising collar in [bound_a, bound_b] by a safeguarded root-find
  // on dprofit == 0, or a bound when the optimum is pinned there; golden section
  // if an endpoint has no usable gradient.
  template <CostCurve K>
  double maximise_profit_over_collar(double bound_a, double bound_b);
  // d(E_up_)/d(collar), a conductance, positive. Forwards to the supply path's
  // duptake_dpsi, which returns NaN at a kink (hazard 6).
  double dE_from_soil_dpsi_collar(double T_collar, const std::vector<double>& psi_soil) {
    require_suction_vector(psi_soil, "dE_from_soil_dpsi_collar");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        return roots_.duptake_dpsi(T_collar, psi_soil);
      default:
        return single_.duptake_dpsi(T_collar, psi_soil);
    }
  }
  // The shut operating point the collar solve's early exits place: the collar
  // where uptake meets the residual flux, the stem where it carries it, and every
  // output from the ordinary evaluator there (hazard 18).
  template <CostCurve K>
  void set_shut_state();
  // The wet bound of every collar solve: the collar where uptake equals the
  // realised residual flux, bracket widened past psi_crit for a dry soil.
  double find_wet_bound(double wettest_soil_layer, double E_up_at_wettest,
                        double E_up_at_crit);
  // The same bound for the shut leaf, bracketed from its floor (hazard 18).
  double find_shut_wet_bound();
  // The wet end of every stem route: the stem potential carrying the residual
  // flux from a pinned upstream potential (that potential at g_min = 0).
  double stem_wet_bound(double psi_upstream);
  // The root on [wettest, end], with the stem held at `end` for the continuity
  // root (find_root_crit = 1) and the zero-uptake root otherwise.
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit,
                       double f_wettest, double f_crit, double end);
  // The same root with the stem at psi_crit, for a caller outside the solve.
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit) {
    const auto f = [&](double x) {
      return find_root_crit == 1 ? E_column(x, psi_soil, psi_crit)
                                 : E_column_zero(x, psi_soil);
    };
    const double f_wettest = f(wettest_soil_layer);
    return find_root_psi(wettest_soil_layer, psi_soil, find_root_crit, f_wettest,
                         f(psi_crit), psi_crit);
  }
  // A bracket bound find_root_psi placed at its 1e-4 tolerance, re-solved to
  // collar_root_tol; called only where the bound is the answer.
  double polish_bracket_bound(double x0, const std::vector<double>& psi_soil,
                              int find_root_crit);
  double find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil);
  double E_column(double x, const std::vector<double>& psi_soil, double psi_leaf);
  double E_column_zero(double x, const std::vector<double>& psi_soil);
  
  // The two temperature responses, exactly linear in ref_value. Only ref_value
  // is templated: the other arguments do not move, and lifting one to T makes
  // every operation carry a derivative vector of zeros. The operation order is
  // load-bearing: reordering moves the double results by an ULP.
  template <typename T>
  T arrh_curve(double Ea, const T& ref_value, double leaf_temp) const;
  template <typename T>
  T peak_arrh_curve(double Ea, const T& ref_value, double leaf_temp, double H_d,
                    double d_S) const;

  // --- Penman-Monteith leaf energy balance (minimal core) ---------------------
  // Recompute the temperature-dependent photosynthetic parameters (vcmax_,
  // jmax_, gamma_, ko_, kc_, R_d_, km_, electron_transport_, and the damage
  // layer's k_imp_eff_, k_rec_eff_, functional_) at a leaf temperature. With the
  // energy balance on it runs per operating point, bypassing the photo_temp cache.
  void update_temperature_dependent_params(double leaf_temp);
  // Saturation vapour pressure es(T) (kPa) and its slope (kPa K^-1), Tetens.
  double saturation_vapour_pressure(double temp) const;
  double saturation_vapour_pressure_slope(double temp) const;

  // Seat `vpd_leaf_` at a leaf temperature: the one place the leaf-to-air deficit
  // is defined. Not in update_temperature_dependent_params, whose differences at
  // T +/- h must not move it; the VPD route into the derivative is dgc_dT's.
  void set_leaf_vpd(double leaf_temp);
  // Explicit energy balance, Tleaf = Tair + (Rn - lambda*E) / (rho*cp/ra + g_rad),
  // E in kg H2O m^-2 s^-1: affine in E, so `dT_dE` (when non-null) is a constant.
  // Not clamped: check_leaf_temperature_range refuses drivers that would need it.
  double leaf_temp_from_E(double E, double* dT_dE = nullptr) const;
  // Refuses a leaf whose temperature can run, over the fluxes the solve can
  // reach, from `coldest` to `hottest` deg C, if either end is past what the
  // model can evaluate: R_d's Q10 not positive, or the Tetens pole.
  void check_leaf_temperature_range(double coldest, double hottest) const;

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;

  // Supply-side transpiration through the stem (kg H2O s^-1 m^-2 LA). Refuses
  // psi_stem < psi_upstream, a flow from the leaf to the soil.
  double transpiration(double psi_stem, double psi_upstream);
  // The flux an operating point uses: zero where psi_upstream >= psi_stem, else
  // transpiration(). A cost or lambda that re-reads the flux must use this, or it
  // refuses a point the evaluator accepted (a prescribed collar at the wet bound).
  double operating_point_flux(double psi_stem, double psi_upstream);
  // Supply-side transpiration, integrated internally (kg H2O s^-1 m^-2 LA).
  double transpiration_full_integration(double psi_stem, double psi_upstream);                    
  // stomatal conductance to CO2 (mol CO2 m^-2 s^-1)
  double stom_cond_CO2(double psi_stem, double psi_upstream);
  // the stem potential carrying a transpiration rate (MPa, positive magnitude)
  double transpiration_to_psi_stem(double transpiration_, double psi_upstream);
  
  // assimilation functions

  double assim_rubisco_limited(double ci_);
  double electron_transport();
  double assim_electron_limited(double ci_);
  double assim_colimited(double ci_);

  // Photosynthetic capacity at this leaf's temperature and light: the two
  // limiting capacities, the colimiting curvature and dark respiration. None
  // moves with ci, so dA/dci walks four directions rather than one per pack slot.
  template <class T>
  struct PhotoCapacity {
    T vcmax;
    T transport;
    T curvature;
    T respiration;

    // ⚠️ Hazard 17: every member that carries a row is visited here, or its
    // columns read exactly zero. The member count is asserted below.
    template <class F>
    void for_each_active(F&& f) const {
      f(vcmax);
      f(transport);
      f(curvature);
      f(respiration);
    }
  };
  static_assert(sizeof(PhotoCapacity<double>) == 4 * sizeof(double),
                "PhotoCapacity gained or lost a member: update for_each_active, "
                "then this count");

  // The coordinates operating_collar placed at an interior collar, with their
  // rows, for billed_outputs_at to read instead of solving again. They live on
  // the caller's tape: a reset or new recording in between leaves them dangling.
  // `at` is the collar they were placed at; NaN at double and any kind but interior.
  template <class S>
  struct HeldCoords {
    double at = std::numeric_limits<double>::quiet_NaN();
    PhotoCapacity<S> cap{};
    S sigma{};
    S gc{};
    S ci{};

    bool placed() const { return !std::isnan(at); }

    // ⚠️ Hazard 17: the marginal's residual reads every one of these, and one
    // left out is a column of exact zeros in the collar's row.
    template <class F>
    void for_each_active(F&& f) const {
      cap.for_each_active(f);
      f(sigma);
      f(gc);
      f(ci);
    }
  };

  // --- the scalar-generic cores -------------------------------------------
  // The single definition of the photosynthesis and cost algebra, instantiated
  // at double for the model and at AD types for its derivative, so the two
  // cannot disagree. A kernel's arguments are its read set (fixed physiology
  // stays a member), and it must not write a member, or an AD probe moves the model.
  template <typename T, typename U> T assim_rubisco_limited_kernel(const T& ci, const U& vcmax) const;
  template <typename T, typename U> T assim_electron_limited_kernel(const T& ci, const U& transport) const;
  template <typename T, typename U> T assim_colimited_kernel(const T& ci,
                                                 const PhotoCapacity<U>& cap) const;
  // The curve parameters take their own scalar U, so pow(tangent, double) pays
  // no log for an exponent that does not move.
  template <typename T, typename U> T proportion_of_conductivity_kernel(
      const T& psi, const U& b, const U& c) const;
  template <typename T, typename U> T hydraulic_cost_TF_kernel(
      const T& psi_stem, const U& scale, const U& b, const U& c,
      const U& beta2) const;

  // The same kernels reading the members, for the forward model.
  template <typename T> T assim_rubisco_limited_kernel(T ci) const;
  template <typename T> T assim_electron_limited_kernel(T ci) const;
  template <typename T> T assim_colimited_kernel(T ci) const;
  template <typename T> T proportion_of_conductivity_kernel(T psi) const;
  template <typename T> T hydraulic_cost_TF_kernel(T psi_stem) const;

  // The kernels' arguments read off a caller's parameter pack. From the pack,
  // never the members: at double they hold the same numbers, so a member read is
  // right in value and a silent zero in every trait's row. The temperature
  // responses are re-run, linear in the 25 C reference, which carries the row.
  template <typename T> T vcmax_at(const leaf_pars<T>& pars) const;
  template <typename T> T jmax_at(const leaf_pars<T>& pars) const;
  template <typename T> T respiration_at(const leaf_pars<T>& pars) const;
  template <typename T> T electron_transport_at(const leaf_pars<T>& pars) const;

  // The one place that says which pack slots each kernel reads.
  template <typename T> PhotoCapacity<T> photo_capacity_at(const leaf_pars<T>& pars) const;
  template <typename T> PhotoCapacity<T> photo_capacity() const;
  // Net assimilation as ci -> infinity: the colimited capacity less respiration,
  // what a leaf with no conductance nets where A reaches zero at no ci (the
  // g_min = 0 limit of a shut leaf in the dark). At any scalar.
  template <typename T, typename U>
  T assim_saturated_kernel(const PhotoCapacity<U>& cap) const;
  template <typename T>
  T hydraulic_cost_TF_kernel(const T& psi_stem, const leaf_pars<T>& pars) const;

  // The slope of a kernel in its first argument, carrying rows: a tangent of a
  // tangent at double through the kernel alone, supplied as rows, since nesting
  // a tangent above the caller's adjoint is forbidden. One direction per argument,
  // including `x`, whose point moves with the traits. `kernel` must declare its
  // return type (`-> U`): a deduced one returns an XAD expression over dead temporaries.
  template <class S, class Kernel, std::same_as<S>... Args>
  S kernel_slope_at(Kernel&& kernel, const char* what, const S& x,
                    const Args&... args) const;

  // dA/dci and dC/dsigma, each a tangent through its own kernel, never a
  // hand-written twin that could disagree with it by a finite, plausible amount.
  template <typename T>
  T assim_slope_at(const T& ci, const PhotoCapacity<T>& cap) const;
  template <typename T>
  T cost_slope_at(const T& sigma, const leaf_pars<T>& pars) const;

  // The STEM's flux at a held stem potential and collar, kmax * (G(sigma) -
  // G(collar)). Not the soil draw: the collar solve meets T1 only to its
  // tolerance, and the model places ci from this flux, so a ci derived from the
  // draw differentiates a residual the model did not solve (finite, and wrong).
  template <typename T>
  T transpiration_at(const T& sigma, const T& collar,
                     const leaf_pars<T>& pars) const;

  // G(psi) at any scalar: the table's value at the passive point, carrying the
  // query slope and the two TRAIT rows the closed form gives exactly.
  template <class S>
  S stem_integral_at(const S& psi, const leaf_pars<S>& pars) const;
  // Q(psi) = G(inf) - G(psi) the same way: the table's value, to its relative
  // precision however small, with the closed form's rows for what is left.
  template <class S>
  S stem_complement_at(const S& psi, const leaf_pars<S>& pars) const;

  // The two conditions that place the operating point at a held collar:
  //
  //   T1  kmax * (G(sigma) - G(collar)) - E_up = 0   places the stem potential
  //   T2  A(ci) * k1 - gc * (ca - ci) * k2   = 0     places the intercellular CO2
  //
  // gc is derived from the draw. `collar_moves` is false at an interior point
  // (the envelope theorem), and the slopes are then absent rather than zero.
  template <class S>
  CollarCoords<S> collar_coords_at(double sigma_star, double ci_star,
                                   const S& collar, const SupplyDraw<S>& draw,
                                   const leaf_pars<S>& pars,
                                   const PhotoCapacity<S>& cap,
                                   bool collar_moves) const;

  // collar_coords_at in two halves, the coordinates at the held collar and then
  // the collar's channel, so a caller holding the first need not solve again.
  template <class S>
  void held_coords_at(double sigma_star, double ci_star, const S& collar,
                      const SupplyDraw<S>& draw, const leaf_pars<S>& pars,
                      const PhotoCapacity<S>& cap, S& sigma_h, S& gc,
                      S& ci_h) const;
  template <class S>
  CollarCoords<S> collar_channel_at(double sigma_star, const S& collar,
                                    const SupplyDraw<S>& draw,
                                    const leaf_pars<S>& pars,
                                    const PhotoCapacity<S>& cap,
                                    const S& sigma_h, const S& gc,
                                    const S& ci_h, bool collar_moves) const;
  // Whether sigma_star is the stem held at a dry pin at the end of its curve
  // (seat_stem_at_dry_pin), where sigma is a function of the traits, not the collar.
  bool stem_held_at_dry_pin(double sigma_star) const {
    return operating_point_kind_ == OperatingPointKind::BoundaryDry &&
           sigma_star == vulnerability_psi_end(stem_b, stem_c);
  }
  // Whether sigma_star is a shut leaf's stem seated at its collar
  // (set_shut_state). There sigma is the collar, rows and all.
  bool stem_held_at_collar(double sigma_star, double collar) const {
    return operating_point_kind_ == OperatingPointKind::HydraulicShutdown &&
           sigma_star == collar;
  }

  // The stem's critical potential at any scalar, derived from the trait pair as
  // at double, so the chain differentiates it.
  template <class S>
  S psi_crit_at(const leaf_pars<S>& pars) const;
  // dry_end<K>() at any scalar: psi_crit_at for a curve bounded there, and the
  // end of the stem curve, b x_end^(1/c), for every other.
  template <CostCurve K, class S>
  S dry_end_at(const leaf_pars<S>& pars) const;
  template <class S>
  S stem_curve_end_at(const leaf_pars<S>& pars) const;

  // The objective at a placed operating point, and everything plant reads.
  // Implemented for TF24 and TF24_floor; a static_assert refuses the rest.
  template <CostCurve K, class S>
  S profit_at(const S& collar, const SupplyDraw<S>& draw,
              const leaf_pars<S>& pars, bool collar_moves) const;
  // The objective at coordinates already placed.
  template <CostCurve K, class S>
  S profit_from_coords(const S& collar, const leaf_pars<S>& pars,
                       const PhotoCapacity<S>& cap, const S& sigma,
                       const S& ci) const;

  // `placed`, when operating_collar placed it at this point, is read instead of
  // solving sigma and ci again. Refused at any kind but interior or at any other
  // collar, either of which means the leaf has moved (hazard 17).
  template <CostCurve K, class S>
  BilledOutputs<S> billed_outputs_at(const S& collar, const SupplyDraw<S>& draw,
                                     const leaf_pars<S>& pars,
                                     const HeldCoords<S>* placed = nullptr) const;

  // M = dprofit/dcollar, a scalar; the collar responses are the coordinates'
  // own slope halves.
  template <CostCurve K, class S>
  S dprofit_dcollar_at(const S& collar, const SupplyDraw<S>& draw,
                       const leaf_pars<S>& pars) const;
  // The same marginal on coordinates operating_collar placed, adding the collar's channel.
  template <CostCurve K, class S>
  S dprofit_dcollar_at(const S& collar, const SupplyDraw<S>& draw,
                       const leaf_pars<S>& pars, const HeldCoords<S>& held) const;
  template <CostCurve K, class S>
  S dprofit_dcollar_from_coords(const S& collar, const leaf_pars<S>& pars,
                                const PhotoCapacity<S>& cap,
                                const CollarCoords<S>& at) const;

  // dM/dp at the solved point, the curvature the interior derivation divides by,
  // as a first difference of M: a tangent is wrong by a factor of eight, because
  // dci/dcollar cancels below the table's fit error. Not const: it moves the
  // model to p0 +/- h and restores the operating point (hazard 8).
  template <CostCurve K>
  double profit_curvature_at_collar();

  // The collar this solve left, at the caller's scalar; the kind chooses which
  // condition closes the system. Not const: the interior arm takes the curvature
  // by moving the model and back, unless passed as `interior_curvature`. `placed`
  // receives the interior arm's coordinates for billed_outputs_at.
  template <CostCurve K, class S>
  S operating_collar(const SupplyDraw<S>& draw, const leaf_pars<S>& pars,
                     double interior_curvature =
                         std::numeric_limits<double>::quiet_NaN(),
                     HeldCoords<S>* placed = nullptr);

  // The smallest |dM/dp| operating_collar's interior arm divides by; nearer
  // zero, every collar row is amplification. Equal to plant's default
  // gradient_curvature_floor, so the two guards refuse the same points.
  static constexpr double interior_curvature_floor = 1e-3;

  // The conditions that pin a collar, one per bound, each an ordinary implicit
  // value. The wet bound: the collar where uptake meets the residual flux (zero
  // uptake at g_min = 0); it reads no cost, so it takes no curve.
  template <class S>
  S wet_bound_at(double bound_x, const SupplyDraw<S>& draw,
                 const leaf_pars<S>& pars) const;
  // The dry end: the continuity root with the stem held at dry_end_at<K>(),
  // which is where the curve enters.
  template <CostCurve K, class S>
  S dry_end_collar_at(double bound_x, const SupplyDraw<S>& draw,
                      const leaf_pars<S>& pars) const;


  // The pack as this leaf stands, so the double path reaches the kernels the way
  // the active path does, with one calling convention.
  leaf_pars<double> passive_pars() const;
  double assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream);
  double psi_stem_to_ci(double psi_stem, double psi_upstream);
  void set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream);


// --- Marginal cost of water ------------------------------------------------
// lambda = dA/dE: every stomatal optimality model satisfies dA/dE = lambda at
// its optimum and differs only in lambda(state). Accessors, not stored state
// (hazard 5). CF77_lambda_ is an input; these are emergent outputs.

  // (dC/dpsi) / (dE/dpsi) for TF24 at a stem potential (MPa), umol CO2 (kg H2O)^-1,
  // with the conductivity cancelled analytically so it stays finite at psi_crit.
  // With TF24_beta2 < 1 it diverges as psi -> 0, a property of the cost.
  double lambda_TF24(double psi_stem) const;
  // Takes both potentials: JS22's cost is a function of the drop, so its
  // marginal cost depends on the upstream potential.
  double lambda_JS22(double psi_stem, double psi_upstream) const;
  // Also both potentials, but only because `E` is measured from the upstream one --
  // the COST here is a function of the absolute potential, unlike JS22's.
  double lambda_CMax(double psi_stem, double psi_upstream) const;
  // lambda_TF24(psi) + lambda_o: the one lambda here that does not vanish at the
  // wet end. One argument, because the hydraulic half reads the absolute
  // potential; lambda_for<TF24_floor> discards the upstream one.
  double lambda_TF24_floor(double psi_stem) const;

  // --- the product-objective family --------------------------------
  // SOX (Eller 2018, 2020) maximises `A * g(psi)` rather than `A - C(psi)`, with
  // `g` the conductivity fraction rescaled onto [0, 1] (TractLSM's `kcost`). No
  // new parameter: `g` reads proportion_of_conductivity and k_crit_fraction.
  double sox_reduction(double psi_stem) const;
  double sox_reduction_deriv(double psi_stem) const;
  double profit_psi_stem_SOX(double psi_stem, double psi_upstream);
  double lambda_SOX(double psi_stem, double psi_upstream) const;

  // --- least-cost / cost-per-unit-carbon --------------------------------------
  // Prentice et al. (2014) and Dong et al. (2017) minimise the cost per unit
  // carbon, `(c_E*E + c_V*Vcmax)/A`; its log reciprocal is `log A - log(c_E*E +
  // c_V*Vcmax)`, the log benefit link of SOX and JW26 with a penalty on capacity
  // as well as flux (theory manuscript S1.4). The only penalty here carrying
  // something other than water.

  // The weighted resource use `W = E_mol + beta*Vcmax_mol`, mol m^-2 s^-1: the
  // one place both unit conversions live, read by the cost, its derivative, the
  // objective and the price so they agree on the units beta is dimensionless in.
  double leastcost_flux(double psi_stem, double psi_upstream);

  // A ratio, `A/W` (umol C per mol of weighted resource), not a carbon flux:
  // monotone in the published cost per unit carbon, so the argmax is the
  // published model's, but the value is not comparable with a TF24 profit.
  double profit_psi_stem_LeastCost(double psi_stem, double psi_upstream);

  // The marginal cost of water, umol C (kg H2O)^-1, `A * kg_to_mol_h2o / W` (eq.
  // S1.4), plus `beta*umol_to_mol*dVcmax/dT*dT/dE / W` with the energy balance
  // on. Unlike SOX's and JW26's, it does not vanish as E -> 0.
  double lambda_LeastCost(double psi_stem, double psi_upstream);

  // Jones et al. (2026): the same product objective with a linear reduction,
  // `1 - psi/psi_crit` (TractLSM's `phiLWP`). No parameter: psi_crit is derived
  // from stem_P50 and stem_c, so this and SOX interpolate between the same two
  // anchors (1 at zero tension, 0 at psi_crit).
  double jw26_reduction(double psi_stem) const;
  double jw26_reduction_deriv(double psi_stem) const;
  double profit_psi_stem_JW26(double psi_stem, double psi_upstream);
  double lambda_JW26(double psi_stem, double psi_upstream) const;

  // ProfitMax's marginal cost of water, umol C (kg H2O)^-1: the normalised dC/dE
  // times |A|max. Requires prepare_profitmax() to have seeded the normalisers.
  double lambda_ProfitMax(double psi_stem);

  // lambda at the operating point as the stem sees it: dA/dE only with the
  // collar fixed (the optimise_psi_stem_* solvers).
  double marginal_cost_water() const;

  // The same in mol CO2 (mol H2O)^-1, the optimality literature's units.
  double marginal_cost_water_molar() const;

  // lambda with the root network's series resistance, which find_root_collar_psi
  // equalises: lambda_TF24 * [1 + kmax*f(psi_r)/S], S = dE_up/dpsi_r. NA where S
  // is NaN (a kink) or zero. Requires a collar solve to have seated the caches.
  double marginal_cost_water_multilayer();

  // Equivalent Medlyn USO slope, kPa^0.5: chi * sqrt(D) / (1 - chi) from the
  // solved chi = ci/ca, comparable with fitted g1 (plantecophys::fitBB).
  double g1_eff() const;

  // The part of the cost that is a price, not a realised carbon loss, umol C
  // m^-2 s^-1: `TF24_floor_lambda_o * E` on that curve, bit-exact against the
  // reported E, and zero on every other. Zero says the curve does not separate
  // the two, not that its cost is all realised carbon (CF77 could be read either way).
  double shadow_cost() const;

  // Conductance to water vapour, mol H2O m^-2 s^-1, from the solved CO2
  // conductance; an accessor (hazard 5). With the energy balance on it is the
  // total leaf-to-air conductance, recovered past the boundary layer's split as
  // g_c / (1 + k g_c): compare a measured gs with stom_cond_stomatal_H2O().
  double stom_cond_H2O() const {
    if (!use_energy_balance_) {
      return stom_cond_CO2_ * H2O_CO2_stom_diff_ratio_;
    }
    return stom_cond_CO2_ /
           (1.0 + boundary_layer_co2_excess() * stom_cond_CO2_) *
           H2O_CO2_stom_diff_ratio_;
  }
  // Boundary-layer conductance to water vapour, mol H2O m^-2 s^-1: the molar
  // form of ra_ at air temperature, scaled by bl_vapour_heat_ratio. +Inf off
  // the energy-balance path, where ra_ is a fallback nothing validates.
  double boundary_layer_cond_H2O() const;
  // The boundary layer's share of the CO2 path, diffusing at r_b = r_s^(2/3)
  // (Sabot et al. 2022 eq 10): g_c = g_c,hat / (1 - k g_c,hat), k = (r_s - r_b)/g_bw,
  // under bl_co2_split_floor; the split returns g_c and dg_c/dg_c,hat. Called
  // only with the energy balance on, so the gate-off path stays bit-identical.
  double boundary_layer_co2_excess() const;
  template <class T> std::pair<T, T> boundary_layer_split_CO2(const T& gc_hat) const;
  // The split's slope at a solved conductance, (1 + k g_c)^2, for the derivative
  // surface; exact above the floor, where every solved point is, refused below.
  template <class T> T boundary_layer_split_slope_at(const T& gc) const;
  // Stomatal conductance to water vapour, 1/g_s = 1/g_total - 1/g_b, mol H2O
  // m^-2 s^-1. Refuses where g_total >= g_b; the refusal is here, not in the
  // solve, so it cannot abort a plant run.
  double stom_cond_stomatal_H2O() const;

  // Read-only: it is derived from the operating point.
  double lambda_emergent() const { return lambda_emergent_; }
  // ProfitMax's normalisers, so a caller can form its equivalence ratio without
  // this class pretending that ratio is a lambda.
  double profitmax_A_max() const { return profitmax_A_max_; }
  double profitmax_k_span() const { return profitmax_k_span_; }

  // Every reported output in one call, for the R boundary. R's
  // .operating_point_names names these by position (test-surface.R), so append,
  // never insert. Index 9 sums finite layers; the rest is as-is, sentinels included.
  std::vector<double> operating_point_values() const {
    double uptake = 0.0;
    for (double c : soil_consumption_) {
      if (std::isfinite(c)) {
        uptake += c;
      }
    }
    return {opt_psi_stem_,   // MPa, positive magnitude
            opt_root_psi_,   // MPa, positive magnitude
            ci_,             // Pa
            assim_colimited_,// umol C m^-2 s^-1
            transpiration_,  // kg H2O m^-2 s^-1
            stom_cond_CO2_,  // mol CO2 m^-2 s^-1
            profit_,         // umol C m^-2 s^-1
            hydraulic_cost_,
            E_up_,
            uptake,
            marginal_cost_water(),
            g1_eff(),
            Tleaf_,          // deg C
            shadow_cost(),
            // The seated curve's marginal cost of water; `lambda` above is TF24's.
            lambda_emergent(),
            // The realised residual flux, kg H2O m^-2 s^-1; `E` above is the total.
            residual_flux_};
  }

// leaf economics functions
  double hydraulic_cost_TF(double psi_stem);

  // Cowan & Farquhar (1977): `CF77_lambda_ * E`, umol C m^-2 s^-1, water priced
  // at a constant marginal value; maximising `A - lambda*E` gives dA/dE == lambda.
  double hydraulic_cost_CF77(double psi_stem, double psi_upstream);

  // The CF77 price, umol C (kg H2O)^-1, and the only place the soil-beta option
  // is applied: the cost, its derivative and lambda_for all read it, so they agree.
  double cf77_price() const;

  // Joshi & Stocker (2022): `gamma*(dpsi)^2`, quadratic in the soil-to-leaf drop.
  // See `JS22_gamma` for what is missing from it.
  double hydraulic_cost_JS22(double psi_stem, double psi_upstream);

  double profit_psi_stem_TF(double psi_stem, double psi_upstream);
  double profit_psi_stem_CF77(double psi_stem, double psi_upstream);
  double profit_psi_stem_JS22(double psi_stem, double psi_upstream);

  // Wolf/Anderegg CMax: the antiderivative of `CMax_a*psi + CMax_b` taken from the
  // upstream potential, so the cost of moving no water is zero.
  double hydraulic_cost_CMax(double psi_stem, double psi_upstream);
  double profit_psi_stem_CMax(double psi_stem, double psi_upstream);

  // The TF24_floor cost, TF24's `gamma*(1-f)^beta2` plus `lambda_o*E`, in carbon
  // units and written into `hydraulic_cost_` like both.
  double hydraulic_cost_TF24_floor(double psi_stem, double psi_upstream);
  double profit_psi_stem_TF24_floor(double psi_stem, double psi_upstream);

  // The instantaneous thermal cost at a leaf temperature, in [0,1]. Zero when the
  // gate is off, so callers need not branch.
  double thermal_cost_at(double leaf_temp) const;
  // Its slope in leaf temperature, r*TC*(1 - TC) with r = 2/(T50 - Tcrit): the
  // sigmoid's own derivative. Zero when the gate is off, like thermal_cost_at.
  double thermal_cost_temp_deriv(double leaf_temp) const;

  // --- the ATLS thermal-damage layer's pieces (see use_thermal_damage_) --------
  // No argument is clamped: an exponential that would overflow is refused, and
  // a shifted optimum or isokinetic temperature outside [250, 350] K likewise.
  //
  // The smooth max(0, x), log1p(exp(s*x))/s. Static: a caller's acclimation
  // forcing uses it too.
  static double softplus_(double x, double s);
  // The entropy term that moves a peaked Arrhenius optimum,
  // Topt = H_d / (d_S - R*ln(Ea/(H_d - Ea))), by delta_topt degrees.
  double d_S_shifted(double Ea, double H_d, double d_S_base,
                     double delta_topt) const;
  // The impairment rate, functional -> impaired, day^-1: k_i_ at the isokinetic
  // temperature k_i_Tiso_ moved by delta_impair degrees. Reads nothing of the
  // Vcmax or Jmax curves.
  double k_imp_eff(double leaf_temp, double delta_impair) const;
  // The repair rate, impaired -> functional, day^-1.
  double k_rec_eff(double leaf_temp) const;
  // The acclimation shift of the impairment's isokinetic temperature at a state A.
  double delta_impair_at(double A) const;

  // The layer's parameters, packed for the templates in `damage`, and the
  // refusal of a negative k_i_Ea_ every entry point applies.
  damage::dyn_pars<double> damage_pars() const;
  void check_damage_pars() const;
  void check_damage_series(const std::vector<double>& times,
                           const std::vector<double>& leaf_temp) const;
  void check_damage_initial_states(double A0, double I0, double D0) const;
  // The names of damage_trajectory_gradient()'s inputs, in index order: the
  // fourteen dynamics parameters, then the three initial states.
  std::vector<std::string> damage_par_names() const;
  // d(state)/d(input) along the trajectory, forward-mode through the same RK4:
  // time fastest, then state (A, I, D), then input (indexing damage_par_names()).
  std::vector<double> damage_trajectory_gradient(
      const std::vector<double>& times, const std::vector<double>& leaf_temp,
      double A0, double I0, double D0, const std::vector<int>& which) const;

  // {dA/dt, dI/dt, dD/dt} per day at a leaf temperature and state, from the
  // layer's parameters alone (no set_physiology, members or gate). dD/dt is
  // maturation only; a caller whose leaf area turns over subtracts the dilution.
  std::vector<double> damage_rates(double leaf_temp, double A, double I,
                                   double D) const;
  // damage_rates integrated by RK4 over a temperature series linear between the
  // times (days, increasing) from (A0, I0, D0), with rate*step < 0.1 and at most
  // 0.25 K per step. Returns 3 * times.size() states, column-major (A, I, D).
  std::vector<double> damage_trajectory(const std::vector<double>& times,
                                        const std::vector<double>& leaf_temp,
                                        double A0, double I0, double D0) const;
  // Throws, naming `who`, when use_thermal_damage_ is on: the paths that read
  // capacities off the parameter pack (vcmax_at and the rest) do not apply it.
  void refuse_thermal_damage(const char* who) const;

  // Sperry (2017) ProfitMax, both terms normalised as the paper defines them
  // (notes above optimise_psi_stem_single): seeds |A|max and the conductance span.
  void prepare_profitmax();

  // The conductance normalisers alone, `k(psi_soil)` and the span to `k_crit`:
  // analytic in the traits, so refreshed whenever the traits move.
  double profitmax_reference_psi();
  void prepare_profitmax_norms();

  // The same with |A|max supplied, not scanned: the scan is piecewise constant
  // in the traits, so a gradient holds |A|max at the base point and refreshes
  // only the analytic normalisers after set_traits clears them.
  void prepare_profitmax_at(double A_max);
  double profit_psi_stem_ProfitMax(double psi_stem, double psi_upstream);

  // The cost/gain/profit curve over [psi_soil, psi_crit] in one R crossing: n
  // rows of (psi, CG, HC, TC, profit), column-major (Sicangco et al. Figs 2, 3, S4).
  std::vector<double> profitmax_curve(int n);

  // Clear the collar-solve outputs a single-layer optimiser does not write, which
  // would otherwise report a point assembled from two solves (hazard 8).
  void clear_collar_solve_state();

  // --- WHICH KIND of operating point the collar solve found -------------------
  // The kind changes what a derivative there means. Recorded by the branch
  // TAKEN, never inferred from a residual: dprofit_at_collar_psi's 0.0 sentinel
  // is indistinguishable from a stationary point.
  enum class OperatingPointKind {
    // No solve since construction, set_traits, or an off-path psi_stem
    // optimiser; the operating point reads NA sentinels.
    Unsolved,
    // Interior profit maximum: dprofit == 0 strictly inside the feasible interval.
    Interior,
    // Constrained optimum at the wet (psi_soil) end, just inside the collar of
    // zero uptake, with a non-zero gradient. Named for the bound that binds, not
    // the conditions: a stressed leaf declines to open here.
    BoundarySoil,
    // Constrained optimum at the dry end, the continuity root with the stem at
    // dry_end<K>(), including a stationary point within dry_pin_resolution of it
    // (hazard 19).
    BoundaryDry,
    // Shutdown on water: no collar both carries more than the residual flux and
    // keeps the stem short of its dry end, so the leaf sits shut at the wet bound
    // (hazard 18). Light and heat never shut a leaf.
    HydraulicShutdown,
    // The collar was imposed (evaluate_root_collar_psi / profit_at_collar_psi),
    // not optimised; dprofit there is the acclimation rate, not a residual.
    Prescribed,
    // The solver could not move: an endpoint had no feasible gradient, or the
    // stationary point is a minimum. Not a pin, which would be plausible and wrong.
    SolverRefused,
    // A partial derivative came back non-finite where feasibility said it could not.
    NonFiniteGradient,
  };

  // This leaf's clamp sites plus the supply model's, one list, summed on read so
  // that rebuilding the root network cannot detach the tally.
  std::vector<std::size_t> clamp_counts() const {
    std::vector<std::size_t> out = *clamps.counts;
    const std::vector<std::size_t>& supply = *roots_.clamps.counts;
    for (std::size_t i = 0; i < out.size() && i < supply.size(); ++i) {
      out[i] += supply[i];
    }
    return out;
  }
  std::size_t clamp_count(int site) const {
    return clamps.at(site) + roots_.clamps.at(site);
  }
  void clear_clamp_counts() const {
    clamps.clear();
    roots_.clamps.clear();
  }

  // Counted from the last enumerator, so an appended kind is inside it; size a
  // per-kind tally from here, never from a repeated literal.
  static constexpr std::size_t operating_point_kind_count =
      static_cast<std::size_t>(OperatingPointKind::NonFiniteGradient) + 1;

  // Read-only: the kind is an output of the solve; replay_operating_point is the
  // one way to put a recorded kind back.
  OperatingPointKind operating_point_kind() const {
    return operating_point_kind_;
  }
  // Put back a recorded operating point, collar and kind in one call:
  // evaluate_root_collar_psi alone tags it `prescribed`, and operating_collar
  // switches on the kind. A BoundarySoil point re-takes the polish of its
  // zero-flux bound, so the rows are the fresh solve's bit for bit (hazard 17).
  double replay_operating_point(double collar, OperatingPointKind kind) {
    return with_curve(cost_curve_, [&](auto tag) {
      return replay_operating_point_for<tag.value>(collar, kind);
    });
  }
  template <CostCurve K>
  double replay_operating_point_for(double collar, OperatingPointKind kind);
  static const char* operating_point_kind_name(OperatingPointKind kind);

private:
  // Where this model replaces a value by a bound: a row severed by a clamp and an
  // honest zero are the same number, so the only way to see one is to count it.
  clamp_counter clamps;

  // Written by every path out of the collar solve and reset to Unsolved at the
  // top of prepare_collar_solve, so a path that forgets reports "unsolved"
  // rather than the previous solve's kind (hazard 8).
  OperatingPointKind operating_point_kind_ = OperatingPointKind::Unsolved;
  // dry_end<K>() of the seated curve for the solve in progress, seated at the
  // top of prepare_collar_solve for the bracket helpers not templated on K.
  double dry_end_ = std::numeric_limits<double>::quiet_NaN();
  // Whether dry_end_ is psi_crit standing in as a first bracket for a curve it
  // does not bound (prepare_collar_solve's first_bracket).
  bool dry_end_provisional_ = false;
  // The wettest rooted layer and the uptake there, from the top of
  // prepare_collar_solve, for extend_collar_bracket.
  double bracket_wettest_ = std::numeric_limits<double>::quiet_NaN();
  double bracket_E_up_at_wettest_ = std::numeric_limits<double>::quiet_NaN();
  // The zero-flux collar a wet pin was solved against, to collar precision, or
  // NaN. wet_bound_at corrects the residual from the pinned point (a step
  // inside) back to it, so the row is the root's derivative, not the step's.
  double wet_bound_root_ = std::numeric_limits<double>::quiet_NaN();
  // The dry end polished to collar precision where the solve read it, or NaN.
  // profit_curvature_at_collar measures from it, as the marginal steepens there.
  double dry_bound_root_ = std::numeric_limits<double>::quiet_NaN();
};

// Human-readable tag. The switch has no default, so a missing name is a -Wswitch
// warning -- an error under the test suite's -Werror=switch, "unknown" elsewhere.
inline const char* Leaf::operating_point_kind_name(OperatingPointKind kind) {
  switch (kind) {
    case OperatingPointKind::Unsolved:          return "unsolved";
    case OperatingPointKind::Interior:          return "interior";
    case OperatingPointKind::BoundarySoil:      return "boundary-soil";
    case OperatingPointKind::BoundaryDry:       return "boundary-dry";
    case OperatingPointKind::HydraulicShutdown: return "hydraulic-shutdown";
    case OperatingPointKind::Prescribed:        return "prescribed";
    case OperatingPointKind::SolverRefused:     return "solver-refused";
    case OperatingPointKind::NonFiniteGradient: return "non-finite-gradient";
  }
  return "unknown";
}


inline Leaf::Leaf()
    :
    vcmax_25(96), // umol m^-2 s^-1 
    stem_c(2.680147), //unitless
    stem_P50(3.4), // MPa, positive magnitude -- THE trait; stem_b/psi_crit derive
    TF24_beta2(1.5), //exponent for effect of hydraulic risk (unitless)
    jmax_25(157.44), // maximum electron transport rate umol m^-2 s^-1
    phi_J(0.30 / 0.85), // quantum yield of electron transport per ABSORBED photon (mol mol^-1); 0.30 per incident photon at alpha_PAR_ = 0.85
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosynthesis equation
    GSS_tol_abs(1e-3),
    vulnerability_curve_ncontrol(ncontrol_default),
    ci_abs_tol(1e-3),
    ci_niter(1000),
    TF24_cost_scale(7.5) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      // The root traits and the beta_R_* constants keep their defaults in
      // MultiLayerRoots, which owns them; a second copy here is hazard 1.
      derive_vulnerability_curves_from_P50();
      // The member, never a literal: set_traits rebuilds at the member's resolution,
      // and curves built at another break the rescale-against-rebuild identity.
      setup_transpiration(vulnerability_curve_ncontrol);
      setup_root_vulnerability(vulnerability_curve_ncontrol);
      setup_clean_leaf();
}

inline Leaf::Leaf(double vcmax_25, double stem_c, double stem_P50,
           double root_c,
           double root_P50,
           double TF24_beta2, double jmax_25,
           double phi_J, double curv_fact_elec_trans, double curv_fact_colim, 
           double GSS_tol_abs,
           double vulnerability_curve_ncontrol,
           double ci_abs_tol,
           double ci_niter,
           double TF24_cost_scale)
    : vcmax_25(vcmax_25), // umol m^-2 s^-1 
    stem_c(stem_c), //unitless
    stem_P50(stem_P50), // MPa, positive magnitude
    TF24_beta2(TF24_beta2), //exponent for effect of hydraulic risk (unitless)
    jmax_25(jmax_25), // maximum electron transport rate umol m^-2 s^-1
    phi_J(phi_J), // quantum yield of electron transport per ABSORBED photon (mol mol^-1)
    curv_fact_elec_trans(curv_fact_elec_trans), //curvature factor for the light response curve (unitless)
    curv_fact_colim(curv_fact_colim), //curvature factor for the colimited photosynthesis equation
    GSS_tol_abs(GSS_tol_abs),
    vulnerability_curve_ncontrol(vulnerability_curve_ncontrol),
    ci_abs_tol(ci_abs_tol),
    ci_niter(ci_niter),
    TF24_cost_scale(TF24_cost_scale) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      // Every psi is a positive magnitude (hazard 2); set_traits shares this check.
      check_psi_magnitudes(stem_P50, stem_c, root_P50, root_c);

      // The root traits belong to the supply path, handed over before its curve is
      // built; root_b and root_psi_crit derive from them there.
      roots_.root_c = root_c;      //unitless
      roots_.root_P50 = root_P50;  // MPa, positive magnitude
      derive_vulnerability_curves_from_P50();

      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_root_vulnerability(vulnerability_curve_ncontrol);
      setup_clean_leaf();
}

inline void Leaf::check_psi_magnitudes(double stem_P50, double stem_c,
                                      double root_P50, double root_c) {
  // Every psi is a positive magnitude in MPa (hazard 2).
  if (!(stem_P50 > 0.0)) {
    util::stop("stem_P50 must be a positive magnitude in MPa; got " +
               util::to_string(stem_P50));
  }
  if (!(root_P50 > 0.0)) {
    util::stop("root_P50 must be a positive magnitude in MPa; got " +
               util::to_string(root_P50));
  }
  // A shape below about 1 (no safe plateau) usually flags a measurement
  // artefact, but it is representable, so only the sign is refused.
  if (!(stem_c > 0.0)) {
    util::stop("stem_c must be positive; got " + util::to_string(stem_c));
  }
  if (!(root_c > 0.0)) {
    util::stop("root_c must be positive; got " + util::to_string(root_c));
  }
}


// The one route to change traits (hazard 10).
inline void Leaf::set_traits(double vcmax_25_, double stem_c_, double stem_P50_,
                             double root_c_, double root_P50_, double TF24_beta2_,
                             double jmax_25_, double phi_J_,
                             double curv_fact_elec_trans_,
                             double curv_fact_colim_,
                             double TF24_cost_scale_, double R_d_25_,
                             double JS22_gamma_, double CMax_a_,
                             double CMax_b_, double g_min_) {
  check_psi_magnitudes(stem_P50_, stem_c_, root_P50_, root_c_);
  // A negative residual conductance would move water into the leaf through
  // its cuticle; refused rather than clamped, like every other sign here.
  if (!(g_min_ >= 0.0)) {
    util::stop("g_min must be a finite, non-negative residual conductance in "
               "mol H2O m^-2 s^-1; got " + util::to_string(g_min_));
  }

  // Both scale parameters and both critical potentials derive from the traits.
  const double stem_b_ = weibull_b_from_P50(stem_P50_, stem_c_);
  const double root_b_ = weibull_b_from_P50(root_P50_, root_c_);
  const double psi_crit_ =
      weibull_psi_at_fraction(stem_b_, stem_c_, k_crit_fraction);
  const double root_psi_crit_ =
      weibull_psi_at_fraction(root_b_, root_c_, k_crit_fraction);

  // Which tables to rebuild, decided before the assignment. A table depends on
  // its curve's c alone, so exact equality is right: an unchanged c gives a
  // bit-identical table, and a gradient's 1e-08 perturbation compares unequal.
  const bool stem_curve_moved = (stem_c_ != stem_c);
  const bool root_curve_moved = (root_c_ != roots_.root_c);

  vcmax_25 = vcmax_25_;
  stem_c = stem_c_;
  stem_P50 = stem_P50_;
  stem_b = stem_b_;
  psi_crit = psi_crit_;
  TF24_beta2 = TF24_beta2_;
  jmax_25 = jmax_25_;
  phi_J = phi_J_;
  curv_fact_elec_trans = curv_fact_elec_trans_;
  curv_fact_colim = curv_fact_colim_;
  TF24_cost_scale = TF24_cost_scale_;
  R_d_25 = R_d_25_;
  JS22_gamma = JS22_gamma_;
  CMax_a = CMax_a_;
  CMax_b = CMax_b_;
  g_min = g_min_;

  roots_.root_c = root_c_;
  roots_.root_P50 = root_P50_;
  roots_.root_b = root_b_;
  roots_.root_psi_crit = root_psi_crit_;

  if (stem_curve_moved) {
    setup_transpiration(vulnerability_curve_ncontrol);
  }
  if (root_curve_moved) {
    setup_root_vulnerability(vulnerability_curve_ncontrol);
  }

  // Last: back to the just-constructed state, outputs NA and both caches reset
  // (hazards 8 and 10). The caller inputs CF77_lambda_, TF24_floor_lambda_o and
  // LeastCost_beta stay standing.
  setup_clean_leaf();
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

inline leaf_pars<double> Leaf::passive_pars() const {
  leaf_pars<double> p{};
  p[par_vcmax_25] = vcmax_25;
  p[par_stem_c] = stem_c;
  p[par_stem_P50] = stem_P50;
  p[par_root_c] = roots_.root_c;
  p[par_root_P50] = roots_.root_P50;
  p[par_TF24_beta2] = TF24_beta2;
  p[par_jmax_25] = jmax_25;
  p[par_phi_J] = phi_J;
  p[par_curv_fact_elec_trans] = curv_fact_elec_trans;
  p[par_curv_fact_colim] = curv_fact_colim;
  p[par_TF24_cost_scale] = TF24_cost_scale;
  p[par_R_d_25] = R_d_25;
  p[par_JS22_gamma] = JS22_gamma;
  p[par_CMax_a] = CMax_a;
  p[par_CMax_b] = CMax_b;
  p[par_g_min] = g_min;
  p[par_kmax] = leaf_specific_conductance_max_;
  // The single path's series resistance; NA on the layered route (per layer).
  p[par_resistance] = roots_.network_.r_R_V_sum.empty() ? util::na_value
                                    : roots_.network_.r_R_V_sum.front();
  p[par_CF77_lambda] = CF77_lambda_;
  p[par_TF24_floor_lambda_o] = TF24_floor_lambda_o;
  p[par_LeastCost_beta] = LeastCost_beta;
  p[par_A] = A_;
  p[par_I] = I_;
  p[par_D] = D_;
  p[par_topt_offset] = topt_offset_;
  p[par_dTopt_max] = dTopt_max_;
  p[par_K_A] = K_A_;
  p[par_k_i] = k_i_;
  p[par_k_i_baseline] = k_i_baseline_;
  p[par_k_rec] = k_rec_;
  p[par_k_rec_Ea] = k_rec_Ea_;
  p[par_c_acclim_maint] = c_acclim_maint_;
  p[par_c_repair_maint] = c_repair_maint_;
  p[par_c_protect_maint] = c_protect_maint_;
  p[par_c_repair_flux] = c_repair_flux_;
  p[par_PPFD] = PPFD_;
  return p;
}

// Reset every derived state and output to the NA sentinel (hazard 8).
inline void Leaf::setup_clean_leaf() {
  ci_ = util::na_value; // Pa
  stom_cond_CO2_= util::na_value; //mol Co2 m^-2 s^-1 
  assim_colimited_= util::na_value; // umol C m^-2 s^-1 
  transpiration_= util::na_value; // kg m^-2 s^-1
  residual_flux_ = util::na_value; // kg m^-2 s^-1
  g_min_leaf_ = util::na_value;    // mol m^-2 s^-1, seated with the deficit
  residual_root_ceiling_ = std::numeric_limits<double>::infinity();
  profit_= util::na_value; // umol C m^-2 s^-1
  // Not CF77_lambda_, TF24_floor_lambda_o or LeastCost_beta: caller inputs, which
  // a reset here would make survive set_drivers() and vanish on set_traits().
  carbon_gain_= util::na_value;
  hydraulic_cost_norm_= util::na_value;
  thermal_cost_= util::na_value;
  profitmax_A_max_= util::na_value;
  profitmax_k_soil_= util::na_value;
  profitmax_k_span_= util::na_value;
  hydraulic_cost_= util::na_value; // umol C m^-2 s^-1
  electron_transport_= util::na_value; //electron transport rate umol m^-2 s^-1
  gamma_= util::na_value;
  ko_= util::na_value;
  kc_= util::na_value;
  km_= util::na_value;
  R_d_= util::na_value;
  leaf_specific_conductance_max_= util::na_value; //kg m^-2 s^-1 MPa^-1 
  vcmax_= util::na_value; // umol m^-2 s^-1
  jmax_= util::na_value; // umol m^-2 s^-1
  opt_root_psi_ = util::na_value; // MPa, positive magnitude
  // No solve has run, so there is no kind of point to report.
  operating_point_kind_ = OperatingPointKind::Unsolved;
  leaf_temp_= util::na_value; // deg C
  photo_temp_= util::na_value; // deg C -- no block seated yet
  Tleaf_= util::na_value; // deg C -- an output, so NA until a solve writes it
  Tair_= util::na_value; // deg C
  Rn_= util::na_value; // W m^-2 (isothermal)
  g_rad_= util::na_value; // W m^-2 K^-1
  ra_= util::na_value; // s m^-1
  PPFD_= util::na_value; //umol m^-2 s^-1
  atm_vpd_= util::na_value; //kPa
  vpd_leaf_= util::na_value; //kPa 
  atm_o2_kpa_= util::na_value; // kPa
  atm_kpa_= util::na_value; // kPa
  umol_per_mol_to_Pa_ = util::na_value;
  ca_= util::na_value; //Pa
  opt_psi_stem_= util::na_value; // MPa, positive magnitude
  opt_ci_= util::na_value; //Pa
  E_up_ = util::na_value;
  medlyn_model_gs_ = util::na_value; // mol CO2 m^-2 s^-1 (Medlyn model)
  theta_w_ = util::na_value;
  theta_fc_ = util::na_value;
  theta_ = util::na_value;
  roots_.clear(); // soil state, geometry and the root resistance network
  // Both supply paths, or a Leaf switched between them gets stale soil state back
  // (hazard 8). single_.clear() keeps grav_head_ (configuration, so
  // set_supply_single is not order-dependent) and clears resistance_ (a driver).
  single_.clear();
  soil_consumption_.clear(); // soil consumption mol  m^-2 s^-1;

  transpiration_cached_ = false; // invalidate transpiration() memo
  photo_temp_cached_ = false;    // members above set to NA; force recompute
}

// Seat the per-timestep physiology that stays constant through the collar or
// stem optimisation: temperature-dependent parameters, assim_max_ (ci = ca),
// the soil state and the supply resistances.
inline void Leaf::set_physiology(const RootNetwork& root_network, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double leaf_temp, double atm_o2_kpa, double atm_kpa) {
  if (psi_soil.size() != soil_depth.size()) {
    util::stop("soil_depth and psi_soil must have the same number of elements");
  }
  if (!std::isfinite(PPFD) || !std::isfinite(leaf_specific_conductance_max) ||
      !std::isfinite(atm_vpd) || !std::isfinite(ca) ||
      !std::isfinite(leaf_temp) ||
      !std::isfinite(atm_o2_kpa) || !std::isfinite(atm_kpa)) {
    util::stop("set_physiology received non-finite scalar input");
  }
  // A real atmosphere's deficit: positive and short of saturation at the air
  // temperature (off the energy balance, leaf_temp), so e_s(Tair) - D > 0.
  // Refused rather than floored, which would bend the leaf temperature.
  if (!(atm_vpd > 0.0) || !(atm_vpd < saturation_vapour_pressure(leaf_temp))) {
    util::stop("set_physiology: atm_vpd must lie in (0, e_s(T_air)), a relative "
               "humidity between 0 and 1; got atm_vpd = " +
               util::format_double(atm_vpd) + " kPa against e_s = " +
               util::format_double(saturation_vapour_pressure(leaf_temp)) +
               " kPa at " + util::format_double(leaf_temp) + " C");
  }
  // The thermal-damage states: a negative fraction of leaf area or the
  // acclimation shift's pole is refused as infeasible, so an integrating
  // consumer can reject the trial step; a NaN passes to poison the outputs.
  if (use_thermal_damage_ && (I_ < 0.0 || D_ < 0.0 || A_ <= -K_A_)) {
    util::stop_infeasible(
        "damage_state",
        "set_physiology: the thermal-damage states need I_ >= 0, D_ >= 0 and "
        "A_ > -K_A_; got A_ = " + util::format_double(A_) + ", I_ = " +
            util::format_double(I_) + ", D_ = " + util::format_double(D_) +
            ", K_A_ = " + util::format_double(K_A_));
  }
  for (size_t i = 0; i < psi_soil.size(); ++i) {
    if (!std::isfinite(psi_soil[i])) {
      util::stop("set_physiology received non-finite psi_soil at layer=" + std::to_string(i) +
                 "; psi_soil=" + util::to_string(psi_soil[i]));
    }
    // Hazard 2: a signed potential fails here rather than running with the soil
    // and the collar on opposite sides of zero.
    if (psi_soil[i] < 0.0) {
      util::stop("set_physiology: psi_soil must be positive magnitudes in MPa, "
                 "not signed potentials (#25); got psi_soil[" +
                 std::to_string(i) + "]=" + util::to_string(psi_soil[i]));
    }
    if (!std::isfinite(soil_depth[i])) {
      util::stop("set_physiology received non-finite soil_depth at layer=" + std::to_string(i) +
                 "; soil_depth=" + util::to_string(soil_depth[i]));
    }
  }
  // The root network is validated in set_root_network below, against the soil
  // profile: it is sized to the deepest rooted layer, which may be shallower.
   atm_vpd_ = atm_vpd;
   leaf_temp_ = leaf_temp;
   atm_kpa_ = atm_kpa;
   umol_per_mol_to_Pa_ = atm_kpa_ * kPa_to_Pa * umol_to_mol;
   atm_o2_kpa_ = atm_o2_kpa;
   PPFD_ = PPFD;
   switch (supply_kind_) {
     case SupplyKind::MultiLayer:
       roots_.set_soil_state(psi_soil, soil_depth);
       break;
     default:
       // One potential: the depth and root-mass profiles are not read, nor any
       // element of psi_soil past the first.
       single_.set_soil_state(psi_soil[0]);
       break;
   }

   leaf_specific_conductance_max_ = leaf_specific_conductance_max;
   // conductance changed -> invalidate the transpiration() memo
   transpiration_cached_ = false;
   ca_ = ca;
   // Energy-balance inputs: the leaf_temp driver is the air temperature, and Rn_
   // is the absorbed shortwave of both bands less the sky's longwave.
   if (use_energy_balance_ && !std::isfinite(NIR_int_)) {
     util::stop("set_physiology: non-finite NIR_int_ on the energy-balance path. "
                "Set the intercepted near-infrared (W m^-2) before set_physiology: "
                "about PPFD / umol_par_per_joule in daylight, 0 under an LED; "
                "NIR_int_=" + util::to_string(NIR_int_));
   }
   Tair_ = leaf_temp_;
   {
     // Rn_ and g_rad_ depend only on the drivers, so they are seated once here;
     // Tleaf enters the balance linearly through g_rad_.
     const double Tair_K = Tair_ + C_to_K;
     // Actual vapour pressure in Pa, as the emissivity coefficient requires;
     // positive because a deficit at or past saturation is refused above.
     const double ea_Pa = 1000.0 * (saturation_vapour_pressure(Tair_) - atm_vpd_);
     // Emissivity held below 1 by a soft minimum: Brutsaert's fit passes unity
     // once ea > 22.25 Pa/K * Tair_K, and a hard min() puts a corner in every
     // output's slope in Tair and humidity.
     const double ema_fit =
         atmos_emissivity_coef *
         std::pow(ea_Pa / Tair_K, atmos_emissivity_exponent);
     const double w = atmos_emissivity_softmin_width;
     const double ema = ema_fit - w * util::softplus((ema_fit - 1.0) / w);
     const double sw_abs =
         alpha_PAR_ * PPFD_ / umol_par_per_joule + alpha_NIR_ * NIR_int_;
     Rn_ = sw_abs - (1.0 - ema) * stefan_boltzmann * Tair_K * Tair_K * Tair_K *
                        Tair_K;
     g_rad_ = 4.0 * leaf_emissivity * stefan_boltzmann * Tair_K * Tair_K * Tair_K;
   }
   // Aerodynamic resistance ra = C_ra*sqrt(d/U). On the energy-balance path a
   // non-finite wind_speed_ or d_ is a broken driver and is refused; zero wind
   // or d (ra -> infinity) falls back to the fixed ra.
   if (use_energy_balance_ &&
       (!std::isfinite(wind_speed_) || !std::isfinite(d_))) {
     util::stop("set_physiology: non-finite wind_speed_/d_ on the energy-balance "
                "path; wind_speed_=" + util::to_string(wind_speed_) +
                "; d_=" + util::to_string(d_));
   }
   ra_ = (std::isfinite(d_) && std::isfinite(wind_speed_) &&
          d_ > 0.0 && wind_speed_ > 0.0)
             ? aerodynamic_resistance_coef * std::sqrt(d_ / wind_speed_)
             : aerodynamic_resistance_fixed;
   // The leaf temperatures the solve can reach: with the balance on, from the
   // shut leaf (E = 0, hottest) to k_max G(inf) (coolest); off it, the driver.
   if (use_energy_balance_) {
     const double E_most = stem_curve_
                               ? leaf_specific_conductance_max_ * stem_curve_limit()
                               : 0.0;
     check_leaf_temperature_range(leaf_temp_from_E(E_most), leaf_temp_from_E(0.0));
   } else {
     check_leaf_temperature_range(leaf_temp_, leaf_temp_);
   }

   // The deficit at the baseline temperature, so assim_max_ and a bare
   // stom_cond_CO2() are defined before a solve re-seats it.
   set_leaf_vpd(leaf_temp_);

   // Temperature/O2 block, cached on photo_temp_key() (hazard 10). Bypassed with
   // the energy balance on (re-seated per candidate at Tleaf) and with the damage
   // gate on (it reads A_, I_ and D_, which the key does not hold).
   const std::array<double, photo_temp_key_size> photo_key = photo_temp_key();
   if (!use_energy_balance_ &&
       !use_thermal_damage_ &&
       photo_temp_cached_ &&
       photo_key == photo_temp_cache_key_) {
     // Cache hit: only electron_transport_ reads the per-call PPFD_.
     electron_transport_ = electron_transport();
   } else {
     update_temperature_dependent_params(leaf_temp_);
     photo_temp_cache_key_ = photo_key;
     photo_temp_cached_ = true;
   }

  // The supply resistances, a validated copy. Both paths read the one argument:
  // the single path takes `r_R_V_sum[0]` as its series resistance. After
  // set_soil_state, because the multi-layer length check needs the soil profile.
  switch (supply_kind_) {
    case SupplyKind::MultiLayer:
      roots_.set_root_network(root_network);
      break;
    default:
      single_.set_supply_resistances(root_network);
      break;
  }

  // Root water uptake per layer; plant writes the crown-integrated value back by
  // name. `.assign`, not `.resize`, so a shallower-rooted plant does not inherit
  // the previous one's deep-layer uptake (hazard 8).
  soil_consumption_.assign(supply_n_layers(), 0.0);

  // Soil moisture for the Medlyn beta_ stress factor, so the standalone Medlyn
  // methods are defined; the root-water path does not read these.
  theta_w_ = theta_w;
  theta_fc_ = theta_fc;
  theta_ = theta;

  // Find maximum assimilation assuming ci = ca
  assim_max_ = assim_colimited(ca_);
}

// ===========================================================================
// ONE REPRESENTATION FOR WATER POTENTIAL (hazard 2)
// ---------------------------------------------------------------------------
// Every psi here, members and arguments alike, is a positive magnitude in MPa,
// asserted at the entry points. Where one potential opposes another the sign
// is in the equation, as in the soil -> collar flux
// E_i = (T_collar - T_soil_i - gravity_head * z_i) / r_R, which is < 0 where a
// layer gains water. Do not call it "tension": that holds only while there is
// no osmotic term and gravity is carried separately in `gravity_head`.
// ===========================================================================
//
// The continuity residual at collar potential `x`: soil uptake less the stem's
// flux to `psi_leaf`, zero where the continuum balances.
inline double Leaf::E_column(double x, const std::vector<double>& psi_soil, double psi_leaf) {

  E_from_Soil_to_Root_Collar(x, psi_soil);
  double E_root_to_leaf = transpiration(psi_leaf, x);
  return E_up_ - E_root_to_leaf;
}

// Soil uptake at collar potential `x`; zero at the zero-uptake collar.
inline double Leaf::E_column_zero(double x, const std::vector<double>& psi_soil) {

  E_from_Soil_to_Root_Collar(x, psi_soil);

  return E_up_;
}

// The collar root on [wettest_soil_layer, end] of E_column with the stem held
// at `end` (find_root_crit == 1) or of E_column_zero (0), to bracket_tol. Both
// are smooth and monotone there, and the caller's early exits guarantee the
// sign change, so a superlinear bracketing solver (TOMS748) is safe.
inline double Leaf::find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit,
                                  double f_wettest, double f_crit, double end) {
  if (find_root_crit == 1) {
    auto target = [&](double x) -> double {
      return E_column(x, psi_soil, end);
    };
    try {
      // The wet end of the final bracket, not its midpoint: the target rises
      // with the collar, and past the root the soil can deliver more than
      // G(inf), which no stem potential carries and the inverse refuses.
      return util::uniroot_smooth_bracket(target, wettest_soil_layer, end,
                                          f_wettest, f_crit, bracket_tol,
                                          ci_niter)
          .first;
    } catch (const std::exception& e) {
      util::stop_infeasible("collar_bracket",
                 "find_root_psi(find_root_crit=1) failed: " + std::string(e.what()) +
                 "; min=" + util::to_string(wettest_soil_layer) +
                 "; max=" + util::to_string(end));
    }
  }

  auto target = [&](double x) -> double {
    return E_column_zero(x, psi_soil);
  };
  try {
    return util::uniroot_smooth(target, wettest_soil_layer, end, f_wettest,
                                f_crit, bracket_tol, ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
               "find_root_psi(find_root_crit=0) failed: " + std::string(e.what()) +
               "; min=" + util::to_string(wettest_soil_layer) +
               "; max=" + util::to_string(end));
  }

}

// Re-solves find_root_psi's 1e-4 MPa bound `x0` to 1e-15, for a bound that is
// the answer. Where a loose root-find stops moves with the parameters through
// its iterates, while the derivative rows (wet_bound_at) differentiate the root
// itself, so an unpolished pinned optimum's root-trait rows were 100x off a
// difference of the solve. Only the answering bound is polished: polishing
// every bound slows every solve.
inline double Leaf::polish_bracket_bound(double x0,
                                         const std::vector<double>& psi_soil,
                                         int find_root_crit) {
  // The wet bound's target is uptake less the residual flux (the uptake itself
  // at g_min = 0).
  auto target = [&](double x) -> double {
    if (find_root_crit == 1) {
      return E_column(x, psi_soil, dry_end_);
    }
    const double up = E_column_zero(x, psi_soil);
    return up - residual_flux(x, up);
  };
  // Tighter than collar_root_tol (~6e-12 at a 5 MPa collar): near the dry end
  // the steep stem curve turns that into ~2e-11 of psi_stem, a 2e-6 row error
  // to a finite difference at a 1e-6 step.
  constexpr double polish_tol = 1e-15;
  const double f0 = target(x0);
  if (f0 == 0.0 || !std::isfinite(f0)) {
    return x0;
  }
  // Both targets rise with the collar, so f0 > 0 puts the root below x0: that
  // side is probed first, starting at the bracket's own tolerance at x0. The
  // continuity root's probes stop at the dry end, past which the collar is
  // drier than the stem and transpiration() refuses the reversed gradient.
  const double hi = find_root_crit == 1 ? dry_end_
                                        : std::numeric_limits<double>::infinity();
  const double sign = f0 > 0.0 ? -1.0 : 1.0;
  for (double d = bracket_tol_at(x0); d <= 0.1; d *= 4.0) {
    for (double side : {sign, -sign}) {
      const double x = std::min(x0 + side * d, hi);
      if (x == x0) {
        continue;
      }
      const double fx = target(x);
      if (std::isfinite(fx) && fx * f0 <= 0.0) {
        const std::pair<double, double> r =
            side < 0.0
                ? util::uniroot_smooth_bracket(target, x, x0, fx, f0, polish_tol,
                                               static_cast<size_t>(ci_niter))
                : util::uniroot_smooth_bracket(target, x0, x, f0, fx, polish_tol,
                                               static_cast<size_t>(ci_niter));
        // The continuity root from its feasible side, for the reason
        // find_root_psi takes it; the wet bound from the middle of its bracket.
        return find_root_crit == 1 ? r.first : 0.5 * (r.first + r.second);
      }
    }
  }
  // No sign change near x0: keep the bracket's answer rather than guess.
  return x0;
}

// The stem potential that carries the soil's uptake at collar `psi_root`.
inline double Leaf::find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil){
  E_from_Soil_to_Root_Collar(psi_root, psi_soil);

  double psi_stem = transpiration_to_psi_stem(E_up_, psi_root);
  return psi_stem;
}

// ---------------------------------------------------------------------------
// MASTER SOLVER: optimal root-collar (and stem) water potential
// ---------------------------------------------------------------------------
// find_root_collar_psi solves soil -> root -> stem -> leaf once per call and
// stores the operating point in opt_psi_stem_, opt_root_psi_ and profit_. The
// inner level is CONTINUITY (find_root_psi on E_column: soil uptake equals the
// stem's flux); the outer is OPTIMISATION, a root-find of dprofit/dcollar = 0
// over [wet bound, collar_at_dry_end] (maximise_profit_over_collar). Where no
// collar in that interval moves water the leaf is shut (set_shut_state).
//
// find_wet_bound: the collar at which soil uptake equals the residual flux the
// leaf realises there, the wet bound and the shut leaf's collar. The target
// rises with the collar, and can change sign past the dry end, so the bracket
// is widened until it does, which it must: the residual flux is capped below
// the roots' ceiling. At g_min = 0 it is the zero-uptake collar, bit-identical
// (hazard 18). `E_up_at_crit` is the uptake at the dry end, NaN if not in hand.
inline double Leaf::find_wet_bound(double wettest_soil_layer,
                                   double E_up_at_wettest,
                                   double E_up_at_crit) {
  const std::vector<double>& soil = supply_psi_soil();
  auto target = [&](double x) -> double {
    E_from_Soil_to_Root_Collar(x, soil);
    return E_up_ - residual_flux(x, E_up_);
  };
  const double f_lo0 =
      E_up_at_wettest - residual_flux(wettest_soil_layer, E_up_at_wettest);
  double lo = wettest_soil_layer, f_lo = f_lo0;
  // The first bracket is [wettest, psi_crit] whatever this solve's dry end, so
  // the wet bound is the same root from the same bracket on every route to it
  // and a replay re-takes it bit for bit.
  double hi = std::max(psi_crit, wettest_soil_layer +
                                     bracket_tol_at(wettest_soil_layer));
  double f_hi = (hi == psi_crit && dry_end_ == psi_crit && std::isfinite(E_up_at_crit))
                    ? E_up_at_crit - residual_flux(psi_crit, E_up_at_crit)
                    : target(hi);
  for (int i = 0; f_hi < 0.0; ++i) {
    if (i >= 64) {
      util::stop_infeasible("collar_bracket",
                            "find_wet_bound: the soil delivers less than the "
                            "residual flux at every collar potential tried; "
                            "wettest=" + util::to_string(wettest_soil_layer) +
                            "; last=" + util::to_string(hi));
    }
    lo = hi;
    f_lo = f_hi;
    hi = hi + (hi - wettest_soil_layer) + bracket_tol_at(hi);
    f_hi = target(hi);
  }
  try {
    return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, bracket_tol,
                                ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
                          "find_wet_bound failed: " + std::string(e.what()) +
                          "; min=" + util::to_string(lo) +
                          "; max=" + util::to_string(hi));
  }
}

// The shut leaf's wet bound, bracketed from below by supply_wet_floor(), the
// lowest gravity balance among conducting layers: the target is <= 0 there and
// rises above, so one root-find to collar precision places it. Where no layer
// conducts the floor is returned, the bound's limit from the wet side; the soil
// potential instead would leave a step of the gravity head in psi_stem. The
// open leaf keeps find_wet_bound, whose bracket keeps g_min = 0 bit-identical.
inline double Leaf::find_shut_wet_bound() {
  const std::vector<double>& soil = supply_psi_soil();
  auto target = [&](double x) -> double {
    E_from_Soil_to_Root_Collar(x, soil);
    return E_up_ - residual_flux(x, E_up_);
  };
  const double floor = supply_wet_floor();
  double lo = floor, f_lo = target(lo);
  if (!(f_lo < 0.0)) {
    return lo;  // nothing conducts, or the floor is the bound itself
  }
  double hi = std::max(supply_wet_ceiling(), lo + bracket_tol_at(lo));
  double f_hi = target(hi);
  for (int i = 0; f_hi < 0.0; ++i) {
    if (i >= 64) {
      util::stop_infeasible("collar_bracket",
                            "find_shut_wet_bound: the soil delivers less than "
                            "the residual flux at every collar tried; floor=" +
                            util::to_string(floor) +
                            "; last=" + util::to_string(hi));
    }
    lo = hi;
    f_lo = f_hi;
    hi = hi + (hi - floor) + bracket_tol_at(hi);
    f_hi = target(hi);
  }
  if (f_hi == 0.0) {
    return hi;
  }
  // polish_bracket_bound's 1e-15, for its reason: the bound is the answer.
  constexpr double shut_tol = 1e-15;
  try {
    return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, shut_tol,
                                static_cast<size_t>(ci_niter));
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
                          "find_shut_wet_bound failed: " + std::string(e.what()) +
                          "; min=" + util::to_string(lo) +
                          "; max=" + util::to_string(hi));
  }
}

// The stem potential carrying exactly the residual flux from `psi_upstream`:
// the wet end of a stem route. It exists because the residual flux is strictly
// below the stem's ceiling, and may lie past psi_crit. Solved to 1e-15 because
// a leaf that declines to open sits exactly here.
inline double Leaf::stem_wet_bound(double psi_upstream) {
  // A stem route has no roots upstream of psi_soil, so it clears the roots'
  // ceiling a previous collar solve seated.
  residual_root_ceiling_ = std::numeric_limits<double>::infinity();
  if (g_min == 0.0) {
    return psi_upstream;
  }
  auto target = [&](double psi) -> double {
    const double E = transpiration(psi, psi_upstream);
    return E - residual_flux(psi_upstream, E);
  };
  double lo = psi_upstream, f_lo = -residual_flux(psi_upstream, 0.0);
  double hi = std::max(psi_crit, psi_upstream + bracket_tol_at(psi_upstream));
  double f_hi = target(hi);
  for (int i = 0; f_hi < 0.0; ++i) {
    if (i >= 64) {
      util::stop_infeasible("stem_wet_bound",
                            "the stem cannot carry the residual flux from "
                            "psi_upstream=" + util::to_string(psi_upstream));
    }
    lo = hi;
    f_lo = f_hi;
    hi = hi + (hi - psi_upstream) + bracket_tol_at(hi);
    f_hi = target(hi);
  }
  constexpr double wet_tol = 1e-15;
  return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, wet_tol,
                              static_cast<size_t>(ci_niter));
}

// The shut operating point every early exit of the collar solve places
// (hazard 18): the collar where uptake meets the residual flux, the stem where
// it carries it, and every output from the ordinary evaluator there, so profit
// is continuous across the exit. Tagged here, not at the call sites, so no
// exit can shut the leaf unclassified.
template <Leaf::CostCurve K>
inline void Leaf::set_shut_state() {
  const std::vector<double>& soil = supply_psi_soil();
  const double collar = find_shut_wet_bound();
  wet_bound_root_ = collar;
  opt_root_psi_ = collar;
  E_from_Soil_to_Root_Collar(collar, soil);
  // The stem carries the uptake. Where that is zero (g_min = 0) or the bound's
  // root-find leaves it at or below zero (a ~1e-50 residual flux at the end of
  // the stem curve), the stem sits at the collar, the flux's limit from above:
  // inverting a negative flux would put the stem wetter than the collar.
  opt_psi_stem_ = (residual_flux(collar, E_up_) == 0.0 || !(E_up_ > 0.0))
                      ? collar
                      : find_psi_stem_from_psi_root(collar, soil);
  profit_ = profit_psi_stem_for<K>(opt_psi_stem_, collar);
  // NA where the supply's conductance is unavailable (a one-layer wet bound is
  // its gravity-balance kink).
  lambda_emergent_ = lambda_collar_for<K>();
  operating_point_kind_ = OperatingPointKind::HydraulicShutdown;
  if (!std::isfinite(profit_)) {
    util::stop_infeasible("collar_solve",
                          "profit is not finite at the shut state; collar=" +
                              util::to_string(collar) +
                              "; psi_stem=" + util::to_string(opt_psi_stem_));
  }
}

// Setup and feasibility for every collar solve. Returns false when the leaf is
// shut, with its operating point already placed; otherwise true with
// [bound_a, bound_b] the feasible collar interval. `first_bracket` lets a curve
// not bounded at psi_crit try psi_crit as a provisional dry end.
template <Leaf::CostCurve K>
inline bool Leaf::prepare_collar_solve(double& bound_a, double& bound_b,
                                       bool first_bracket){

  // ⚠️ Cleared first, so a path that declines to classify reports Unsolved
  // rather than the previous plant's kind (hazard 8).
  operating_point_kind_ = OperatingPointKind::Unsolved;
  wet_bound_root_ = std::numeric_limits<double>::quiet_NaN();
  dry_bound_root_ = std::numeric_limits<double>::quiet_NaN();

  // Builds the supply's per-solve caches; returns the wettest rooted layer, the
  // smallest suction and the lower bracket end.
  const double wettest_soil_layer = supply_begin_solve();
  // Seated before anything reads the residual flux. Only the layered supply
  // bounds its uptake; a series resistance delivers without limit.
  residual_root_ceiling_ =
      (g_min > 0.0 && supply_kind_ == SupplyKind::MultiLayer)
          ? roots_.uptake_limit<double>(roots_.held_supply())
          : std::numeric_limits<double>::infinity();
  const double E_up_at_wettest =
      E_column_zero(wettest_soil_layer, supply_psi_soil());

  bracket_wettest_ = wettest_soil_layer;
  bracket_E_up_at_wettest_ = E_up_at_wettest;

  // One bracket with the stem held at `end` (seated in dry_end_ first, which the
  // helpers read): false if no collar short of it both supplies the cuticle and
  // leaves the stem flow, else [wet bound, continuity root]. The two root-finds
  // share endpoint values, handed over rather than re-taken: at the end the
  // column term is exactly zero, and at the wettest layer both read one uptake.
  auto bracket = [&](double end) -> bool {
    dry_end_ = end;
    if (wettest_soil_layer >= end) {
      return false;
    }
    // The most the soil delivers with the stem short of `end`.
    const double E_at_end = E_column(end, supply_psi_soil(), end);
    if (E_at_end < residual_flux(end, E_at_end)) {
      return false;
    }
    const double collar_at_dry_end = find_root_psi(
        wettest_soil_layer, supply_psi_soil(), 1,
        E_up_at_wettest - transpiration(end, wettest_soil_layer), E_at_end, end);
    if (collar_at_dry_end >= end) {
      return false;
    }
    bound_a = find_wet_bound(wettest_soil_layer, E_up_at_wettest, E_at_end);
    // The root's own critical potential bounds nothing: a clamp to it puts a
    // kink in the stem and root traits wherever the two tie.
    bound_b = collar_at_dry_end;
    // Ends of an interval narrower than bracket_tol cannot be trusted, so they
    // are solved to collar precision before anything is decided on them.
    if (bound_b - bound_a < bracket_tol_at(bound_b)) {
      bound_a = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
      bound_b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    }
    return bound_b - bound_a > collar_root_tol;
  };

  // For a curve not bounded at psi_crit, psi_crit is only a first bracket: it
  // holds the optimum wherever profit already falls there, and
  // maximise_profit_over_collar extends it otherwise. Where it holds no interval
  // the end of the stem curve is tried before the leaf is called shut.
  if constexpr (!cost_domain_ends_at_psi_crit<K>()) {
    if (first_bracket && bracket(psi_crit)) {
      dry_end_provisional_ = true;
      return true;
    }
  }
  dry_end_provisional_ = false;
  if (!bracket(dry_end<K>())) {
    set_shut_state<K>();
    return false;
  }
  return true;
}

template <Leaf::CostCurve K>
inline bool Leaf::seat_stem_at_dry_pin() {
  if constexpr (cost_domain_ends_at_psi_crit<K>()) {
    return false;
  } else {
    if (operating_point_kind_ != OperatingPointKind::BoundaryDry) {
      return false;
    }
    opt_psi_stem_ = vulnerability_psi_end(stem_b, stem_c);
    return true;
  }
}

inline double Leaf::extend_collar_bracket(double bound_a) {
  const double end = vulnerability_psi_end(stem_b, stem_c);
  dry_end_ = end;
  dry_end_provisional_ = false;
  const double E_at_end = E_column(end, supply_psi_soil(), end);
  const double collar_at_dry_end = find_root_psi(
      bracket_wettest_, supply_psi_soil(), 1,
      bracket_E_up_at_wettest_ - transpiration(end, bracket_wettest_), E_at_end,
      end);
  // A narrow interval's dry end is polished, as in prepare_collar_solve, or a
  // root within bracket_tol_at of the wet bound can land on it and report shut
  // a leaf whose profit still rises. The wet end is not: moving it would strand
  // the caller's `lo`, already stepped in from it.
  const double bound_b = std::max(collar_at_dry_end, bound_a);
  if (bound_b - bound_a < bracket_tol_at(bound_b)) {
    return std::max(polish_bracket_bound(bound_b, supply_psi_soil(), 1),
                    bound_a);
  }
  return bound_b;
}

// The profit-maximising collar on [bound_a, bound_b], by a bracketing
// root-find of dprofit/dpsi == 0; prepare_collar_solve must have seated the
// soil caches. The first-order condition, not a search of profit: a search's
// argmax wanders within its tolerance as parameters move, an offset that
// swamps, and can reverse, the hydraulic traits' true response (hazard 3).
// It is util::maximise_over_closed_interval_foc at n = 0 plus the
// OperatingPointKind every exit sets, which pinned golden rows read: the one
// solver plus a classification (.claude/CLAUDE.md "One solver").
template <Leaf::CostCurve K>
inline double Leaf::maximise_profit_over_collar(double bound_a, double bound_b) {
  double width = bound_b - bound_a;
  const auto dprofit = [&](double psi, bool* feasible = nullptr) {
    return dprofit_at_collar_psi<K>(psi, feasible);
  };

  // Places lo at the wet bound: the answer where profit falls from the wet end,
  // and for a dark leaf under a product objective (below).
  double lo = bound_a;
  double hi = bound_b;
  const auto place_at_wet_bound = [&]() {
    const double a = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
    wet_bound_root_ = a;
    if (g_min > 0.0) {
      // With a residual conductance the bound is a regular point (no stomatal
      // flux, the cuticle conducts), and a leaf that declines to open sits on it.
      lo = a;
    } else if (a < hi) {
      // With nothing conducting at the bound a dark leaf has no finite ci there
      // (ci = ca + R_d/g_c), so the state is taken a step inside, always from
      // the polished root and on the bound's own scale: a fraction of the width
      // moves with the dry bound, a dependence the rows do not have. Not below
      // 1e-6: the rows' error grows with the step, but at 1e-7 the stem
      // inversion reaches its rounding floor. Only the no-flow test is repeated:
      // f_lo's sign cannot flip across the ~1e-4 the bound moved.
      const double scale = std::max(a, 1e-3);
      const double limit = 0.5 * (bound_b - a);
      for (double frac = 1e-6; frac * scale < limit; frac *= 10.0) {
        const double x = a + frac * scale;
        const double s = find_psi_stem_from_psi_root(x, supply_psi_soil());
        if (std::isfinite(s) && x < s) {
          lo = x;
          break;
        }
      }
    }
  };

  // A product objective with A <= 0 at the dry end has no maximum, and the leaf
  // is shut (objective_is_product). Decided before either end's gradient is
  // read: through the log link at a negative A they say nothing.
  if constexpr (objective_is_product<K>()) {
    if (dark_under_product<K>(
            find_psi_stem_from_psi_root(bound_b, supply_psi_soil()), bound_b)) {
      // Where the open path would start: the bound at g_min > 0, a step inside
      // the no-flow sentinel at g_min = 0, should placement find no step.
      lo = g_min > 0.0 ? bound_a : bound_a + 1e-6 * width;
      place_at_wet_bound();
      operating_point_kind_ = OperatingPointKind::BoundarySoil;
      return lo;
    }
  }

  // The WET endpoint. With g_min > 0 it is an ordinary point, tried as-is. At
  // g_min = 0 it is the zero-uptake collar, where dprofit's 0.0 is a sentinel a
  // bracketing solver would report as the optimum, so the read steps inside
  // until the gradient is real (the infeasible region is under 1e-6 of the
  // bracket, so the first step lands). Starting at the bound there instead
  // moves every interior row by the solver floor.
  bool ok_lo = false;
  double f_lo = 0.0;
  for (double frac = (g_min > 0.0 ? 0.0 : 1e-6); frac < 0.5;
       frac = (frac == 0.0 ? 1e-6 : frac * 10.0)) {
    lo = bound_a + frac * width;
    f_lo = dprofit(lo, &ok_lo);
    if (ok_lo && std::isfinite(f_lo)) {
      break;
    }
  }

  // The DRY endpoint. For a curve bounded at psi_crit it is tried as-is. For any
  // other it is the continuity root with the stem at the end of its curve,
  // placed only to bracket_tol_at, and a read there can fall past the root,
  // where transpiration() refuses; so the first read is a step inside.
  // Profit can still rise at that end, which is the model: where the stem has
  // little left to carry and A climbs steeply in E, the leaf transpires all the
  // stem can carry, a pin at the dry end. The probes only say where to look;
  // the pin is decided at the end itself, solved to collar precision (hazard 19).
  // On a first bracket (stem at psi_crit) the end is an ordinary point: profit
  // falling there keeps the bracket, anything else extends it.
  constexpr bool bounded = cost_domain_ends_at_psi_crit<K>();
  bool ok_hi = false;
  bool rising = false;
  double f_hi = 0.0;
  bool settled = false;
  if (dry_end_provisional_) {
    f_hi = dprofit(hi, &ok_hi);
    settled = ok_hi && std::isfinite(f_hi) && f_hi < 0.0;
    if (!settled) {
      bound_b = extend_collar_bracket(bound_a);
      width = bound_b - bound_a;
    }
  }
  // The bracket the interior root-find runs on: the step inside the wet bound,
  // and `hi` wherever the probes leave it.
  double root_lo = lo;
  double f_root_lo = f_lo;
  // Whether the dry end has decided the point: an interior maximum between
  // the first probe and the end, or the pin. Neither needs the probes below.
  bool decided = false;
  if (!settled) {
    hi = bounded ? bound_b : bound_b - 1e-6 * width;
    f_hi = dprofit(hi, &ok_hi);
  }
  // ⚠️ Where the first probe reads a rising profit the end decides (hazard 19):
  // the marginal falls steeply over the last micro-MPa, so a probe's sign is
  // not the end's. Solved to collar precision, falling at the end is an
  // interior maximum between probe and end, rising is the pin; deciding on the
  // probes made the solve jump. Read before the deeper probes, each a ci solve.
  // Only where profit rises from the wet end too: falling from it is two
  // maxima, which the refusal below compares.
  if (!bounded && !settled && f_lo > 0.0 && ok_hi && std::isfinite(f_hi) &&
      f_hi >= 0.0) {
    dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    bool ok_end = false;
    const double f_end = dprofit(dry_bound_root_, &ok_end);
    if (ok_end && std::isfinite(f_end) && f_end < 0.0 && dry_bound_root_ > hi) {
      root_lo = hi;
      f_root_lo = f_hi;
      hi = dry_bound_root_;
      f_hi = f_end;
      decided = true;
    } else if (ok_end && std::isfinite(f_end) && f_end >= 0.0) {
      rising = true;
      decided = true;
    }
    bound_b = dry_bound_root_;
  }
  if (!settled && !decided &&
      !(ok_hi && std::isfinite(f_hi) && (bounded || f_hi < 0.0))) {
    rising = ok_hi && std::isfinite(f_hi);
    for (double frac = bounded ? 1e-6 : 1e-5; frac < 0.5; frac *= 10.0) {
      hi = bound_b - frac * width;
      f_hi = dprofit(hi, &ok_hi);
      if (ok_hi && std::isfinite(f_hi) && (bounded || f_hi < 0.0)) {
        rising = false;
        break;
      }
      rising = rising || (ok_hi && std::isfinite(f_hi));
    }
  }
  if (!bounded && rising) {
    // Rising at the end, or at every probe where the end was not read: the pin,
    // on the end solved to collar precision.
    if (!std::isfinite(dry_bound_root_)) {
      dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      bound_b = dry_bound_root_;
    }
    hi = dry_bound_root_;
    f_hi = 1.0;
    ok_hi = true;
  }

  // No usable gradient at one end, which no golden-grid row reaches: a signal,
  // answered by golden section rather than a failure. `!ok` is the supply
  // reporting no informative gradient (SolverRefused); `ok && !isfinite` is a
  // partial feasibility said would be finite (NonFiniteGradient). A NaN would
  // otherwise fail both pin tests below and reach the root-find unbracketed.
  if (!ok_lo || !ok_hi || !std::isfinite(f_lo) || !std::isfinite(f_hi)) {
    operating_point_kind_ =
        ((ok_lo && !std::isfinite(f_lo)) || (ok_hi && !std::isfinite(f_hi)))
            ? OperatingPointKind::NonFiniteGradient
            : OperatingPointKind::SolverRefused;
    return util::golden_section_max(
        [&](double bound) {
          const double psi_stem =
              find_psi_stem_from_psi_root(bound, supply_psi_soil());
          return profit_psi_stem_for<K>(psi_stem, bound);
        },
        bound_a, bound_b, GSS_tol_abs);
  }

  // A CONSTRAINED optimum: profit still climbs at one end, so the answer is that
  // end, which is common on dry soil. Return lo/hi, never the raw bound_a/bound_b:
  // at g_min = 0 the wet bound is the sentinel stepped over above. Where a row
  // crosses between interior and pinned the argmax has a real kink in traits.
  // A pin's bound is the answer, so only a side that would pin is polished
  // (polish_bracket_bound says why); an interior solve pays nothing.
  // dprofit <= 0 at the wet end and >= 0 at the dry end is not a third pin: the
  // interior stationary point is a minimum, so it is SolverRefused and the
  // larger of the two ends' profits is returned.
  if (f_lo <= 0.0) {
    place_at_wet_bound();
  }
  if (f_hi >= 0.0 && hi == bound_b && !std::isfinite(dry_bound_root_)) {
    const double b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    // f_hi's sign cannot flip across the ~1e-4 the bound moved.
    if (b != bound_b && b > lo) {
      hi = b;
    }
    dry_bound_root_ = hi;
  }

  if (f_lo <= 0.0 && f_hi >= 0.0) {
    operating_point_kind_ = OperatingPointKind::SolverRefused;
    const double p_lo = profit_psi_stem_for<K>(
        find_psi_stem_from_psi_root(lo, supply_psi_soil()), lo);
    const double p_hi = profit_psi_stem_for<K>(
        find_psi_stem_from_psi_root(hi, supply_psi_soil()), hi);
    if (!std::isfinite(p_hi)) {
      return lo;
    }
    if (!std::isfinite(p_lo)) {
      return hi;
    }
    return p_hi > p_lo ? hi : lo;
  }
  if (f_lo <= 0.0) {
    operating_point_kind_ = OperatingPointKind::BoundarySoil;
    return lo;
  }
  if (f_hi >= 0.0) {
    operating_point_kind_ = OperatingPointKind::BoundaryDry;
    return hi;
  }

  // Interior stationary point: f_lo > 0 > f_hi, a valid bracket.
  operating_point_kind_ = OperatingPointKind::Interior;
  double root;
  if (!bounded && decided && std::isfinite(dry_bound_root_)) {
    // Near the end the marginal is close to linear in the log of the distance
    // to it (dry_end_neighbourhood), so the root is found in that coordinate: a
    // few superlinear steps where the potential takes a dozen ci solves. A root
    // within dry_pin_resolution of the end is the pin.
    const double end = dry_bound_root_;
    bool ok_in = false;
    const double f_in = end - root_lo > dry_pin_resolution
                            ? dprofit(end - dry_pin_resolution, &ok_in)
                            : 0.0;
    if (!(ok_in && std::isfinite(f_in) && f_in < 0.0)) {
      operating_point_kind_ = OperatingPointKind::BoundaryDry;
      return end;
    }
    const auto in_log_distance = [&](double s) {
      return dprofit(end - std::exp(s));
    };
    // 1e-10 in the logarithm is 1e-12 MPa at the neighbourhood's outer edge.
    const double s = util::uniroot_smooth(
        in_log_distance, std::log(dry_pin_resolution), std::log(end - root_lo),
        f_in, f_root_lo, 1e-10, static_cast<size_t>(ci_niter));
    root = end - std::exp(s);
  } else {
    root = util::uniroot_smooth(dprofit, root_lo, hi, f_root_lo, f_hi,
                                collar_root_tol, static_cast<size_t>(ci_niter));
  }
  if constexpr (!bounded) {
    // A root within dry_end_neighbourhood of the dry end reads that end for its
    // curvature, so the end is polished; within dry_pin_resolution it is the
    // pin. Only where the bracket's end is the domain's: a first bracket's end
    // (stem at psi_crit) is an ordinary point, and a replay, bracketing the
    // whole domain, would take a different stencil.
    if (!dry_end_provisional_ &&
        bound_b - root < dry_end_neighbourhood + bracket_tol_at(bound_b)) {
      if (!std::isfinite(dry_bound_root_)) {
        dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      }
      if (dry_bound_root_ - root < dry_pin_resolution) {
        operating_point_kind_ = OperatingPointKind::BoundaryDry;
        return dry_bound_root_;
      }
    }
  }
  return root;
}

template <Leaf::CostCurve K>
inline void Leaf::find_root_collar_psi_for(){
  // The stem route's parameter check: without it CF77 with an unset lambda
  // maximises a NaN objective and returns a bracket end as an operating point.
  check_cost_parameters<K>();
  // ProfitMax seeds its normalisers before searching, on either topology: the
  // reference potential is psi_soil where there is one and the zero-flow collar
  // where there is not (see profitmax_reference_psi).
  if constexpr (K == CostCurve::ProfitMax) {
    prepare_profitmax();
  }
    double bound_a, bound_b;
    if (!prepare_collar_solve<K>(bound_a, bound_b, true)) {
      return;
    }

    const double opt_root_psi = maximise_profit_over_collar<K>(bound_a, bound_b);

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());
    seat_stem_at_dry_pin<K>();

    opt_root_psi_ = opt_root_psi;
    profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi);
    // On the single-layer solvers' axis: the collar is the free variable, so the
    // marginal cost carries the soil-to-collar series-resistance correction.
    lambda_emergent_ = lambda_collar_for<K>();

    if(!std::isfinite(profit_)){
        util::stop_infeasible("collar_solve", "non-finite profit; opt_psi_stem_=" + util::to_string(opt_psi_stem_) +
             "; opt_root_psi_=" + util::to_string(opt_root_psi_) +
             "; bound_a=" + util::to_string(bound_a) +
             "; bound_b=" + util::to_string(bound_b) +
             "; E_up_=" + util::to_string(E_up_) +
             "; assim_colimited_=" + util::to_string(assim_colimited_) +
             "; hydraulic_cost_=" + util::to_string(hydraulic_cost_));
    }
}

// The collar solve on TF24 whatever curve is seated; optimise() solves the
// seated model.
inline void Leaf::find_root_collar_psi() {
  find_root_collar_psi_for<CostCurve::TF24>();
}


// Names in, model seated: this parses, and the typed overload alone assigns, so
// R and C++ cannot hold two ideas of what seating means.
inline void Leaf::set_model(const std::string& curve,
                           const std::string& route,
                           const std::string& method) {
  int k = -1;
  for (int i = 0; i < n_cost_curves; ++i) {
    if (curve_name(i) == curve) { k = i; break; }
  }
  if (k < 0) {
    std::string names;
    for (int i = 0; i < n_cost_curves; ++i) {
      names += (i == 0 ? "" : ", ") + curve_name(i);
    }
    util::stop("unknown cost curve \"" + curve + "\". Available: " + names + ".");
  }
  if (route != "collar" && route != "stem") {
    util::stop("unknown route \"" + route + "\": use \"collar\" for the full "
               "soil-to-collar path or \"stem\" for psi_stem with the upstream "
               "potential pinned at psi_soil.");
  }
  // Only the name is checked here: what a method can do depends on the curve,
  // route and energy balance, all settable after this call, so the configuration
  // is checked at solve time.
  if (method != "exact" && method != "closed") {
    util::stop("unknown method \"" + method + "\": use \"exact\" for the "
               "first-order-condition root-find over the closed interval, or "
               "\"closed\" for the closed form. Available: exact, closed.");
  }
  set_model(static_cast<CostCurve>(k), route == "collar", method == "closed");
}

// The same, already parsed, for a C++ caller holding a `CostCurve`
// (`gradient::route_seat` seats once per observation).
inline void Leaf::set_model(CostCurve curve, bool collar_route,
                           bool closed_method) {
  cost_curve_ = curve;
  route_is_collar_ = collar_route;
  method_is_closed_ = closed_method;
  // The curve's constants, and whether the method can serve this configuration,
  // are checked at solve time: a caller may seat the model, then set `CF77_lambda_`.
}


// The seated model, solved: the closed form, or the exact solve on the seated
// curve (with_curve) and route.
inline void Leaf::optimise() {
  last_solve_fell_back_ = false;
  if (method_is_closed_) {
    optimise_closed();
    return;
  }
  const bool collar = route_is_collar_;
  with_curve(cost_curve_, [&](auto tag) {
    constexpr CostCurve K = tag.value;
    if (collar) {
      find_root_collar_psi_for<K>();
    } else {
      optimise_psi_stem_single<K>();
    }
  });
}


// The operating point at a prescribed collar, on the seated curve. A collar
// outside the feasible interval is refused as infeasible (profit_at_collar_psi),
// so a caller tracking the collar as a state can reject the step.
inline double Leaf::evaluate_root_collar_psi(double target_opt_root_psi) {
  return with_curve(cost_curve_, [&](auto tag) {
    return evaluate_root_collar_psi_for<tag.value>(target_opt_root_psi);
  });
}


inline double Leaf::dprofit_droot_collar_psi(double opt_root_psi, bool* feasible) {
  return with_curve(cost_curve_, [&](auto tag) {
    return dprofit_droot_collar_psi_for<tag.value>(opt_root_psi, feasible);
  });
}


template <Leaf::CostCurve K>
inline double Leaf::evaluate_root_collar_psi_for(double target_opt_root_psi){
    double bound_a, bound_b;
    if (!prepare_collar_solve<K>(bound_a, bound_b)) {
      // Shut: profit_ is already set.
      return profit_;
    }
    // The bracket's dry end is placed only to its tolerance, and an interior
    // maximum can sit between it and the continuity root: a target past the
    // bracket is read against the root, solved to collar precision.
    if (target_opt_root_psi > bound_b) {
      const double b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      if (b > bound_b) {
        dry_bound_root_ = b;
        bound_b = b;
      }
    }

    return profit_at_collar_psi<K>(target_opt_root_psi, bound_a, bound_b);
}

template <Leaf::CostCurve K>
inline double Leaf::replay_operating_point_for(double collar,
                                               OperatingPointKind kind) {
    double bound_a, bound_b;
    double out;
    if (!prepare_collar_solve<K>(bound_a, bound_b)) {
      out = profit_;
    } else {
      if (kind == OperatingPointKind::BoundarySoil) {
        wet_bound_root_ = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
      }
      // A dry pin sits on the polished end, and an interior point near it reads
      // that end for its curvature: seat it as the solve did, under the solve's
      // own test, bracket tolerance included, so the two take the same stencil.
      if (kind == OperatingPointKind::BoundaryDry ||
          (kind == OperatingPointKind::Interior &&
           bound_b - collar < dry_end_neighbourhood + bracket_tol_at(bound_b))) {
        dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      }
      // A pinned solve's collar is a polished bound, up to ~1e-4 outside the
      // bracket returned here; clamped back it would replay a different point,
      // so the interval is widened to hold it.
      out = profit_at_collar_psi<K>(collar, std::min(bound_a, collar),
                                    std::max(bound_b, collar));
    }
    operating_point_kind_ = kind;
    if (seat_stem_at_dry_pin<K>()) {
      profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi_);
      lambda_emergent_ = lambda_collar_for<K>();
      out = profit_;
    }
    return out;
}

// The operating point at a prescribed collar, after prepare_collar_solve: a
// caller evaluating several collars in one step reuses the soil-side caches.
// A collar outside [bound_a, bound_b] is refused, not projected: projected, the
// outputs follow the bound rather than the collar asked for. Infeasible-class,
// so a consumer integrating the collar can reject the step, not end the run.
template <Leaf::CostCurve K>
inline double Leaf::profit_at_collar_psi(double target_opt_root_psi,
                                  double bound_a, double bound_b){
    if (!(target_opt_root_psi >= bound_a && target_opt_root_psi <= bound_b)) {
      util::stop_infeasible(
          "prescribed_potential",
          "the prescribed collar potential " +
              util::format_double(target_opt_root_psi) +
              " MPa is outside the feasible collar interval [" +
              util::format_double(bound_a) + ", " +
              util::format_double(bound_b) + "] MPa");
    }
    const double opt_root_psi = target_opt_root_psi;

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());
    opt_root_psi_ = opt_root_psi;
    profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi);
    // dC/dE is a property of the point, defined off the optimum, where it is not dA/dE.
    lambda_emergent_ = lambda_collar_for<K>();
    // Here a non-zero dprofit is the answer (TF24f's acclimation rate), not a
    // residual to check against zero.
    operating_point_kind_ = OperatingPointKind::Prescribed;

    if(!std::isfinite(profit_)){
        util::stop_infeasible("collar_solve",
             "non-finite profit in evaluate_root_collar_psi; "
             "target=" + util::to_string(target_opt_root_psi) +
             "; opt_root_psi=" + util::to_string(opt_root_psi) +
             "; bound_a=" + util::to_string(bound_a) +
             "; bound_b=" + util::to_string(bound_b));
    }
    return profit_;
}

// Exact d(profit)/d(collar): A'(ci) dci/dpsi - C'(psi_stem) dpsi_stem/dpsi,
// with dci/dpsi from the implicit-function theorem on the ci residual and A', C'
// by forward AD of the model's own kernels; vignette("the-models"), "The
// benefit, and the root-find inside it". This wrapper seats the supply caches on
// psi_soil_; the collar solve calls dprofit_at_collar_psi, seated once.
template <Leaf::CostCurve K>
inline double Leaf::dprofit_droot_collar_psi_for(double opt_root_psi, bool* feasible) {
  supply_begin_solve();
  return dprofit_at_collar_psi<K>(opt_root_psi, feasible);
}

inline std::vector<double> Leaf::dprofit_droot_collar_psi_checked(
    double opt_root_psi) {
  bool feasible = false;
  const double d = dprofit_droot_collar_psi(opt_root_psi, &feasible);
  return {d, feasible ? 1.0 : 0.0};
}

template <Leaf::CostCurve K>
inline double Leaf::dprofit_at_collar_psi(double opt_root_psi, bool* feasible) {
  using AD = xad::fwd<double>::active_type;
  const double psi = opt_root_psi;
  // Infeasible until the exits below are passed: their 0.0 must be
  // distinguishable from a stationary point.
  if (feasible != nullptr) {
    *feasible = false;
  }

  // Past the continuity root with the stem at the end of its curve the soil
  // delivers more than the stem can carry: no operating point, reported like the
  // other no-gradient exits rather than thrown (a difference across a dry pin reaches it).
  E_from_Soil_to_Root_Collar(psi, supply_psi_soil());
  if (!stem_carries(E_up_, psi)) {
    return 0.0;
  }
  const double psi_stem = transpiration_to_psi_stem(E_up_, psi);
  // psi >= psi_stem is the no-flow / reversed case, caught before the ci solve:
  // there psi_stem_to_ci throws (transpiration() refuses the pair) rather than
  // returning non-finite, and the throw would end a consumer's run.
  if (!std::isfinite(psi_stem) || psi >= psi_stem) {
    return 0.0;  // shut-down / infeasible: no informative gradient
  }

  // With the energy balance on, the temperature parameters are seated at this
  // candidate's Tleaf(E), as the objective is: otherwise this is the first-order
  // condition of the air-temperature model. Same order as
  // set_leaf_states_rates_from_psi_stem: transpiration, temperature, ci solve.
  const double E = transpiration(psi_stem, psi);
  double dT_dE = 0.0;
  double Tleaf_here = leaf_temp_;
  if (use_energy_balance_) {
    Tleaf_here = leaf_temp_from_E(E, &dT_dE);
    update_temperature_dependent_params(Tleaf_here);
    set_leaf_vpd(Tleaf_here);
  }

  // Below the wet bound the stomatal flux is zero, the corner the solve keeps
  // (E_s >= 0); with nothing conducting (the g_min = 0 wet bound) there is no
  // derivative to take.
  const double r_star = residual_flux(psi, E);
  const double E_s = E - r_star < 0.0 ? 0.0 : E - r_star;
  const double cuticle = cuticular_cond_CO2();
  if (E_s == 0.0 && cuticle == 0.0) {
    return 0.0;
  }

  const double ci = psi_stem_to_ci(psi_stem, psi);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  // dpsi_stem/dpsi, computed once here and shared by both branches below.
  const double dEup_dpsi = dE_from_soil_dpsi_collar(psi, supply_psi_soil());
  double dpsistem_dpsi;
  if (std::isfinite(dEup_dpsi)) {
    E_from_Soil_to_Root_Collar(psi, supply_psi_soil());  // refresh E_up_ at psi
    const double E_psi_stem =
        E_up_ / leaf_specific_conductance_max_ +
        stem_curve_integral(psi, "Leaf::dprofit_at_collar_psi, forming "
                                 "dpsi_stem/dpsi at the operating point");
    const double dEpsistem_dpsi =
        dEup_dpsi / leaf_specific_conductance_max_ + stem_curve_integral_deriv(psi);
    // dpsi_stem/dE_psi_stem = 1/G'(G^-1(E_psi_stem)), and G^-1(E_psi_stem) is
    // the psi_stem inverted above from the same E_up_ and psi: read the slope there.
    (void)E_psi_stem;
    dpsistem_dpsi = (1.0 / stem_curve_integral_deriv(psi_stem)) * dEpsistem_dpsi;
  } else {
    // NaN from the supply means "use finite differences" (hazard 6).
    const double h = 1e-6;
    dpsistem_dpsi =
        (find_psi_stem_from_psi_root(psi + h, supply_psi_soil()) -
         find_psi_stem_from_psi_root(psi - h, supply_psi_soil())) / (2.0 * h);
  }


  // Past both exits: whatever is returned below is a real derivative, so a zero
  // from here IS a stationary point.
  if (feasible != nullptr) {
    *feasible = true;
  }

  // A'(ci) by forward AD of the kernel assim_colimited() instantiates, so it is
  // the derivative of the function evaluated, not of a hand-kept mirror.
  AD ci_ad = ci;            xad::derivative(ci_ad) = 1.0;
  const double A_prime = xad::derivative(assim_colimited_kernel(ci_ad));
  const double C_prime = cost_deriv<K>(psi_stem, psi);

  // gc and its partials, on the stem table's own slope. The residual flux reads
  // psi_stem only through Tleaf (the energy-balance term) and the collar through
  // its ceiling (below); both are exact zeros at g_min = 0.
  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  double gc = gc_const * E_s + cuticle;
  double dgc_dpsistem =
      gc_const * leaf_specific_conductance_max_ * stem_curve_integral_deriv(psi_stem);
  double dgc_dpsi =
      gc_const * leaf_specific_conductance_max_ * (-stem_curve_integral_deriv(psi));
  const double drstar_dpsi = residual_flux_dpsi_upstream(psi, Tleaf_here);
  dgc_dpsi = dgc_dpsi - gc_const * drstar_dpsi;
  // The boundary layer's share of the CO2 path: the value and every partial
  // through the split's one slope. [[unlikely]] is measured: without it this
  // branch slows the gate-off collar solve, and moving it out of line does not help.
  double bl_slope = 1.0;
  if (use_energy_balance_) [[unlikely]] {
    const std::pair<double, double> split = boundary_layer_split_CO2(gc);
    gc = split.first;
    bl_slope = split.second;
    dgc_dpsistem *= bl_slope;
    dgc_dpsi *= bl_slope;
  }

  // IFT on g(ci; psi_stem, psi): dci/dp = -(dg/dp)/(dg/dci).
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;      // dg/dci
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;
  const double dci_dpsi_expl = -(-dgc_dpsi * (ca_ - ci) * inv_atm) / g_ci;


  const double dci_dpsi = dci_dpsistem * dpsistem_dpsi + dci_dpsi_expl;
  // The identity arm keeps every difference objective's operation order, which
  // the golden files pin; only a non-identity link takes the general path.
  double base;
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    base = A_prime * dci_dpsi - C_prime * dpsistem_dpsi;
  } else {
    const double A = assim_colimited_kernel(ci);
    base = benefit_link_deriv<K>(A) * (A_prime * dci_dpsi) - C_prime * dpsistem_dpsi;
    // A log link has no derivative where the leaf is a net source (A <= 0), an
    // ordinary point once the cuticle conducts CO2: reported as absent, never NaN.
    if (!std::isfinite(base)) {
      if (feasible != nullptr) {
        *feasible = false;
      }
      return 0.0;
    }
  }
  // The cost's direct dependence on the collar, appended so a curve that does not
  // read upstream keeps the arm above operation for operation.
  if constexpr (cost_reads_upstream<K>()) {
    base = base - cost_deriv_upstream<K>(psi_stem, psi);
  }
  // Gated at the CALL SITE, not just inside the callee: the block is out of line
  // (deliberately, so adding it cannot change FMA contraction in this inlined
  // body), and an out-of-line call costs even when it returns 0.0 immediately.
  if (!use_energy_balance_) {
    return base;
  }
  return base + dprofit_energy_balance_term<K>(
                    ci, gc, g_ci, inv_atm, gc_const, A_prime, dgc_dpsistem,
                    dgc_dpsi, dpsistem_dpsi, dT_dE, Tleaf_here, E_s, r_star,
                    drstar_dpsi,
                    residual_cond_temp_deriv(psi, Tleaf_here, gc_const), bl_slope);
}


// dC/dpsi_stem for one cost curve: the only place a cost curve enters a
// derivative, shared by both routes (the collar route multiplies it by
// dpsi_stem/dpsi and adds cost_deriv_upstream<K>). Every arm is explicit and the
// last asserts, so a new curve is a compile error rather than a plausible number.
template <Leaf::CostCurve K>
inline double Leaf::cost_deriv(double psi_stem, double psi_upstream) {
  using AD = xad::fwd<double>::active_type;
  (void)psi_upstream;
  double C_prime;
  if constexpr (K == CostCurve::TF24) {
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    C_prime = xad::derivative(hydraulic_cost_TF_kernel(ps_ad));
  } else if constexpr (K == CostCurve::CF77) {
    // lambda * E, and dE/dpsi_stem is kmax times the conductivity fraction,
    // which stem_curve_integral_deriv returns.
    C_prime = cf77_price() * leaf_specific_conductance_max_ *
              stem_curve_integral_deriv(psi_stem);
  } else if constexpr (K == CostCurve::JS22) {
    // C = gamma*(psi - psi_up)^2.
    C_prime = 2.0 * JS22_gamma * (psi_stem - psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // TF24's arm plus CF77's, each in its own operand order, so both reductions
    // stay bit-exact (see `hydraulic_cost_TF24_floor`).
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    C_prime = xad::derivative(hydraulic_cost_TF_kernel(ps_ad)) +
              TF24_floor_lambda_o * leaf_specific_conductance_max_ *
                  stem_curve_integral_deriv(psi_stem);
  } else if constexpr (K == CostCurve::SOX) {
    // Under the log link the cost is `-log g`, so dC/dpsi is `-g'/g`; `g` is zero
    // only at psi_crit, which no optimiser returns.
    C_prime = -sox_reduction_deriv(psi_stem) / sox_reduction(psi_stem);
  } else if constexpr (K == CostCurve::JW26) {
    // The same, and here `-g'/g` collapses: with `g = 1 - psi/psi_crit` it is
    // exactly `1/(psi_crit - psi)`.
    C_prime = 1.0 / (psi_crit - psi_stem);
  } else if constexpr (K == CostCurve::LeastCost) {
    // Under the log link the cost is `+log W`, `W = E_mol + beta*Vcmax_mol`
    // (`leastcost_flux`), so dC/dpsi is `W'/W` with `W' = kg_to_mol_h2o * kmax * f`:
    // the opposite sign to SOX's `-g'/g`, as W rises where g falls. Vcmax's slope
    // in Tleaf is taken in `dprofit_energy_balance_term`.
    C_prime = kg_to_mol_h2o * leaf_specific_conductance_max_ *
              stem_curve_integral_deriv(psi_stem) /
              leastcost_flux(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::ProfitMax) {
    // C = [k(psi_s) - k(psi)]/k_span, so dC/dpsi = kmax*|f'|/k_span. The thermal
    // half's slope (energy balance on) is taken in dprofit_energy_balance_term.
    const double f = proportion_of_conductivity(psi_stem);
    const double abs_fprime =
        f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
    C_prime = leaf_specific_conductance_max_ * abs_fprime / profitmax_k_span_;
  } else {
    static_assert(K == CostCurve::CMax, "unhandled CostCurve in cost_deriv");
    // Linear in the absolute potential: the form Anderegg et al. state.
    C_prime = CMax_a * psi_stem + CMax_b;
  }
  return C_prime;
}


// Whether a curve's cost reads the upstream potential as well as psi_stem, which
// matters on the collar route, where upstream is the decision variable. A curve
// that reads it but says false here solves the wrong first-order condition with
// no error; `test_collar_foc_matches_profit` checks every curve. Every arm is
// explicit for the reason `cost_deriv`'s are.
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_reads_upstream() {
  if constexpr (K == CostCurve::CF77 || K == CostCurve::JS22 ||
                K == CostCurve::CMax || K == CostCurve::TF24_floor ||
                K == CostCurve::LeastCost) {
    // LeastCost: `log W` reads upstream through E, the only log-link cost that does.
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::SOX ||
                      K == CostCurve::JW26 || K == CostCurve::ProfitMax,
                  "unhandled CostCurve in cost_reads_upstream");
    return false;
  }
}


// Whether a curve's own form ends at psi_crit. For the others the marginal cost
// diverges as the stem dries, so the dry end is the end of the stem curve
// (vulnerability_psi_end). For these four the bound belongs to the curve:
//
//   SOX, JW26  the gain factor reaches zero at psi_crit and is negative past it
//   CF77       a constant price need not overtake A'(E) before the stem's ceiling
//   LeastCost  sits at psi_crit on most rows, dark and light
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_domain_ends_at_psi_crit() {
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26 ||
                K == CostCurve::CF77 || K == CostCurve::LeastCost) {
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::JS22 ||
                      K == CostCurve::CMax || K == CostCurve::ProfitMax ||
                      K == CostCurve::TF24_floor,
                  "unhandled CostCurve in cost_domain_ends_at_psi_crit");
    return false;
  }
}

// SOX and JW26 maximise A g, g falling to 0 at psi_crit. With A <= 0 at every
// opening the supremum is 0 as g -> 0, driving a dark leaf to psi_crit, so a
// leaf with A <= 0 at the dry end (A rises with opening) is shut at the wet
// bound instead. The argmax jumps at the light compensation point regardless.
template <Leaf::CostCurve K>
constexpr bool Leaf::objective_is_product() {
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26) {
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::CF77 ||
                      K == CostCurve::JS22 || K == CostCurve::CMax ||
                      K == CostCurve::ProfitMax || K == CostCurve::TF24_floor ||
                      K == CostCurve::LeastCost,
                  "unhandled CostCurve in objective_is_product");
    return false;
  }
}

template <Leaf::CostCurve K>
inline bool Leaf::dark_under_product(double psi_dry, double psi_upstream) {
  static_assert(objective_is_product<K>(),
                "dark_under_product is only defined for a product objective");
  set_leaf_states_rates_from_psi_stem(psi_dry, psi_upstream);
  return !(assim_colimited_ > 0.0);
}

template <Leaf::CostCurve K>
inline double Leaf::dry_end() const {
  if constexpr (cost_domain_ends_at_psi_crit<K>()) {
    return psi_crit;
  } else {
    return vulnerability_psi_end(stem_b, stem_c);
  }
}

// dC/dpsi_upstream at fixed psi_stem, for the curves `cost_reads_upstream` names.
// The collar FOC subtracts it with a unit chain factor, since upstream is the
// collar there.
template <Leaf::CostCurve K>
inline double Leaf::cost_deriv_upstream(double psi_stem, double psi_upstream) {
  static_assert(cost_reads_upstream<K>(),
                "cost_deriv_upstream is only defined where the cost reads upstream");
  if constexpr (K == CostCurve::CF77) {
    // lambda * kmax * (G(psi_stem) - G(psi_up)), so -lambda * kmax * G'(psi_up),
    // with G' the table's own slope: the solve ran on the table.
    (void)psi_stem;
    return -(cf77_price() * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream));
  } else if constexpr (K == CostCurve::JS22) {
    return -2.0 * JS22_gamma * (psi_stem - psi_upstream);
  } else if constexpr (K == CostCurve::CMax) {
    // d/du of (psi - u)*(a*(psi + u)/2 + b) is -(a*u + b): the integrand at the
    // lower limit, as for any integral of the price.
    (void)psi_stem;
    return -(CMax_a * psi_upstream + CMax_b);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // TF24's half reads psi_stem alone; the floor's lambda_o*E is CF77's arm.
    (void)psi_stem;
    return -(TF24_floor_lambda_o * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream));
  } else {
    static_assert(K == CostCurve::LeastCost, "unhandled CostCurve");
    // `(dW/du)/W` with `dW/du = -kg_to_mol_h2o * kmax * G'(u)`: CF77's arm with
    // kg-to-mol in place of the price, over `cost_deriv`'s `W`.
    return -(kg_to_mol_h2o * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream)) /
           leastcost_flux(psi_stem, psi_upstream);
  }
}

// dprofit/dpsi_stem = h'(A) A'(ci) dci/dpsi_stem - C'(psi_stem), upstream held
// fixed, for every curve through its `BenefitLink` h. A deliberate copy of
// dprofit_at_collar_psi's algebra: equal in value at dpsistem_dpsi = 1 and
// dgc_dpsi = 0, but not in operation order, which the golden files pin.
template <Leaf::CostCurve K>
inline double Leaf::dprofit_dpsi_stem(double psi_stem, double psi_upstream,
                                      bool* feasible) {
  using AD = xad::fwd<double>::active_type;
  if (feasible != nullptr) {
    *feasible = false;
  }
  // No flow (upstream >= stem), so no informative gradient: the collar
  // version's sentinel and contract.
  if (!std::isfinite(psi_stem) || psi_upstream >= psi_stem) {
    return 0.0;
  }

  const double E = transpiration(psi_stem, psi_upstream);
  double dT_dE = 0.0;
  double Tleaf_here = leaf_temp_;
  if (use_energy_balance_) {
    Tleaf_here = leaf_temp_from_E(E, &dT_dE);
    update_temperature_dependent_params(Tleaf_here);
    set_leaf_vpd(Tleaf_here);
  }

  // As on the collar route: the stomatal flux the residual leaves, read as
  // zero below the wet bound, and no derivative where nothing conducts.
  const double r_star = residual_flux(psi_upstream, E);
  const double E_s = E - r_star < 0.0 ? 0.0 : E - r_star;
  const double cuticle = cuticular_cond_CO2();
  if (E_s == 0.0 && cuticle == 0.0) {
    return 0.0;
  }

  const double ci = psi_stem_to_ci(psi_stem, psi_upstream);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  const double C_prime = cost_deriv<K>(psi_stem, psi_upstream);

  if (feasible != nullptr) {
    *feasible = true;
  }

  AD ci_ad = ci;
  xad::derivative(ci_ad) = 1.0;
  const double A_prime = xad::derivative(assim_colimited_kernel(ci_ad));

  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  double gc = gc_const * E_s + cuticle;
  double dgc_dpsistem =
      gc_const * leaf_specific_conductance_max_ *
      stem_curve_integral_deriv(psi_stem);
  // The boundary layer's share of the CO2 path, as on the collar route.
  double bl_slope = 1.0;
  if (use_energy_balance_) [[unlikely]] {
    const std::pair<double, double> split = boundary_layer_split_CO2(gc);
    gc = split.first;
    bl_slope = split.second;
    dgc_dpsistem *= bl_slope;
  }
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;

  // The identity arm is written out: a literal 1.0 factor is exact but can change
  // FMA contraction, and the golden files compare the last bit.
  double base;
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    base = A_prime * dci_dpsistem - C_prime;
  } else {
    const double A = assim_colimited_kernel(ci);
    base = benefit_link_deriv<K>(A) * (A_prime * dci_dpsistem) - C_prime;
    // As on the collar route: no derivative at a net source under a log link,
    // reported rather than returned as NaN.
    if (!std::isfinite(base)) {
      if (feasible != nullptr) {
        *feasible = false;
      }
      return 0.0;
    }
  }
  if (!use_energy_balance_) {
    return base;
  }
  // dgc_dpsi = 0 and dpsistem_dpsi = 1, so the term's dE/dpsi is kmax f(psi_stem).
  // The term applies the curve's link itself, templated on K: an unscaled term
  // added here would be plausible and wrong.
  return base + dprofit_energy_balance_term<K>(
                    ci, gc, g_ci, inv_atm, gc_const, A_prime, dgc_dpsistem, 0.0,
                    1.0, dT_dE, Tleaf_here, E_s, r_star, 0.0,
                    residual_cond_temp_deriv(psi_upstream, Tleaf_here, gc_const),
                    bl_slope);
}

// The operating point at a prescribed psi_stem. See the declaration for why a
// gradient here differs from one at an optimum.
template <Leaf::CostCurve K>
inline double Leaf::evaluate_psi_stem(double target_psi_stem) {
  clear_collar_solve_state();

  if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }
  // The optimiser's checks and ProfitMax's normalisers, so a prescribed point
  // and a solved one agree on what a curve requires and whether its objective exists.
  check_cost_parameters<K>();
  if constexpr (K == CostCurve::ProfitMax) {
    prepare_profitmax();
  }
  const double psi_soil = supply_psi_soil_scalar();
  // The wet end: where the stem carries exactly the residual flux, psi_soil
  // itself at g_min = 0.
  const double psi_wet = stem_wet_bound(psi_soil);

  // A soil too dry for the stem to carry more than the cuticle's loss short of
  // its dry end leaves one feasible potential, and the target is irrelevant to it.
  const double psi_dry = dry_end<K>();
  if (psi_wet > psi_dry) {
    opt_psi_stem_ = psi_wet;
    profit_ = profit_psi_stem_for<K>(psi_wet, psi_soil);
    lambda_emergent_ = lambda_for<K>(psi_wet, psi_soil);
    operating_point_kind_ = OperatingPointKind::Prescribed;
    return profit_;
  }

  // Refused, not projected onto [psi_wet, psi_dry], for profit_at_collar_psi's reason.
  if (!(target_psi_stem >= psi_wet && target_psi_stem <= psi_dry)) {
    util::stop_infeasible(
        "prescribed_potential",
        "the prescribed stem potential " + util::format_double(target_psi_stem) +
            " MPa is outside the feasible interval [" +
            util::format_double(psi_wet) + ", " + util::format_double(psi_dry) +
            "] MPa");
  }
  const double psi = target_psi_stem;

  opt_psi_stem_ = psi;
  profit_ = profit_psi_stem_for<K>(psi, psi_soil);
  lambda_emergent_ = lambda_for<K>(psi, psi_soil);
  operating_point_kind_ = OperatingPointKind::Prescribed;
  return profit_;
}

// See the declaration. The step is 1e-3 K, matching dprofit_energy_balance_term's
// dA/dT difference; the full restore makes calling this inside a derivative safe.
template <class F>
inline double Leaf::temp_deriv(double leaf_temp, F read) {
  const double h = 1e-3;
  const double vc0 = vcmax_, jm0 = jmax_, ga0 = gamma_, ko0 = ko_,
               kc0 = kc_, rd0 = R_d_, km0 = km_, J0 = electron_transport_;
  const double kimp0 = k_imp_eff_, krec0 = k_rec_eff_, func0 = functional_;
  const double pt0 = photo_temp_;
  update_temperature_dependent_params(leaf_temp + h);
  const double up = read();
  update_temperature_dependent_params(leaf_temp - h);
  const double dn = read();
  vcmax_ = vc0; jmax_ = jm0; gamma_ = ga0; ko_ = ko0; kc_ = kc0;
  R_d_ = rd0; km_ = km0; electron_transport_ = J0;
  k_imp_eff_ = kimp0; k_rec_eff_ = krec0; functional_ = func0;
  photo_temp_ = pt0;
  return (up - dn) / (2.0 * h);
}

inline double Leaf::vcmax_temp_deriv(double leaf_temp) {
  return temp_deriv(leaf_temp, [&] { return vcmax_; });
}

// The energy-balance correction to dprofit/dpsi; zero with the gate off. Out of
// line on purpose (hazard 5): a branch inside the inlined cost core can change
// FMA contraction and so the gate-off golden file (`nm -C test_golden | grep dprofit`).
//
// With the gate on psi also reaches A through Tleaf, directly and through the
// ci residual, which gains an explicit T. By the implicit-function theorem
//
//   dA/dT (total) = [A_T*gc*inv_atm + A_prime*dgc_dT*(ca-ci)*inv_atm] / g_ci
//
// times tau = dT/dE * dE/dpsi. Checks on it: at dgc_dT = 0 the damping factor
// (gc*inv_atm)/g_ci lies in (0,1); tau < 0, so above the thermal optimum
// (A_T < 0) the optimum moves drier, toward evaporative cooling. dgc_dT is not
// zero: a form putting it under A_T rather than A_prime passes every test at
// dgc_dT = 0 and is wrong. Derived in vignette("the-models"), Appendix B, "The
// energy-balance term".
template <Leaf::CostCurve K>
inline double Leaf::dprofit_energy_balance_term(
    double ci, double gc, double g_ci, double inv_atm, double gc_const,
    double A_prime, double dgc_dpsistem, double dgc_dpsi, double dpsistem_dpsi,
    double dT_dE, double Tleaf, double E_s, double r_star, double drstar_dpsi,
    double dgc_dT_residual, double bl_slope) {
  if (!use_energy_balance_ || dT_dE == 0.0) {
    return 0.0;
  }

  // dA/dTleaf by a central difference over the model's own temperature block,
  // so it cannot drift from the evaluated function. Not a templated block: that
  // changes what shares an inlined body, and moved the gate-off results. ⚠️ h is
  // fixed and absolute: relative to T, the estimator becomes a function of psi
  // and puts a staircase in the argmax (hazard 3).
  const double h = 1e-3;
  const double vcmax0 = vcmax_, jmax0 = jmax_, gamma0 = gamma_, ko0 = ko_,
               kc0 = kc_, Rd0 = R_d_, km0 = km_, J0 = electron_transport_;
  const double kimp0 = k_imp_eff_, krec0 = k_rec_eff_, func0 = functional_;
  const double pt0 = photo_temp_;
  // Tleaf, not leaf_temp_, which the energy-balance path reads as air temperature.
  const double T0 = Tleaf;

  update_temperature_dependent_params(T0 + h);
  const double A_up = assim_colimited(ci);
  const double vcmax_up = vcmax_;
  update_temperature_dependent_params(T0 - h);
  const double A_dn = assim_colimited(ci);
  const double vcmax_dn = vcmax_;

  // Restored by assignment: exact, and cheaper than a third update.
  vcmax_ = vcmax0; jmax_ = jmax0; gamma_ = gamma0; ko_ = ko0; kc_ = kc0;
  R_d_ = Rd0; km_ = km0; electron_transport_ = J0;
  k_imp_eff_ = kimp0; k_rec_eff_ = krec0; functional_ = func0;
  photo_temp_ = pt0;

  const double A_T = (A_up - A_dn) / (2.0 * h);

  // dE/dpsi from the gc partials: dividing out gc_const and the boundary layer's
  // slope leaves the stomatal flux's slope, and the residual's is added back,
  // since Tleaf reads the total flux.
  const double dE_dpsi =
      (dgc_dpsistem * dpsistem_dpsi + dgc_dpsi) / (gc_const * bl_slope) +
      drstar_dpsi;
  const double tau = dT_dE * dE_dpsi;

  // dgc/dTleaf at fixed E. Stomatal part: -gc_s * esat'(Tleaf)/D, negative (a
  // hotter leaf has a larger deficit) and zero on the deficit's floor; the
  // residual's part is appended. Both pass through the boundary layer's split,
  // which reads no leaf temperature.
  const double gc_s = gc_const * E_s;
  const double dgc_dT =
      (((vpd_leaf_ > vpd_leaf_min)
            ? -gc_s * saturation_vapour_pressure_slope(Tleaf) / vpd_leaf_
            : 0.0) +
       dgc_dT_residual) *
      bl_slope;

  // dA/dTleaf total times tau: dA/dpsi along the temperature path, benefit-side
  // only. The costs that read Tleaf (ProfitMax's thermal cost, LeastCost's
  // capacity term) are taken below. The identity arm keeps the additive curves'
  // operation order.
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    (void)vcmax_up;
    (void)vcmax_dn;
    return tau *
           (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) / g_ci;
  } else {
    // Under a link, J = h(A(psi)) - C(psi) with A the total dependence, so this
    // path takes the call site's single h'(A); h'(A) is never differentiated.
    const double dA_dpsi_T =
        tau * (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) /
        g_ci;
    // The block is restored, so this is the A the call site's `base` used.
    const double A = assim_colimited_kernel(ci);
    double term = benefit_link_deriv<K>(A) * dA_dpsi_T;
    if constexpr (K == CostCurve::ProfitMax) {
      // dTC/dpsi = TC'(Tleaf) * tau; cost-side, so no h' factor, and the slope
      // lambda_ProfitMax reports.
      term -= thermal_cost_temp_deriv(Tleaf) * tau;
    } else if constexpr (K == CostCurve::LeastCost) {
      // C = log W with `vcmax_` at Tleaf, so d(log W)/dpsi gains
      // beta*umol_to_mol*dVcmax/dT*tau/W (`cost_deriv` holds only the E half).
      // Cost-side, so no h'. W reads the total flux, stomatal part from gc.
      const double dVcmax_dT = (vcmax_up - vcmax_dn) / (2.0 * h);
      const double W = gc_s / gc_const * kg_to_mol_h2o +
                       LeastCost_beta * vcmax_ * umol_to_mol +
                       r_star * kg_to_mol_h2o;
      term -= LeastCost_beta * umol_to_mol * dVcmax_dT * tau / W;
    } else {
      (void)vcmax_up;
      (void)vcmax_dn;
    }
    return term;
  }
}

template <typename T>
inline T Leaf::arrh_curve(double Ea, const T& ref_value, double leaf_temp) const {
  return ref_value*std::exp(Ea*((leaf_temp+C_to_K) - (25 + C_to_K))/((25 + C_to_K)*gas_constant*(leaf_temp+C_to_K)));
}

template <typename T>
inline T Leaf::peak_arrh_curve(double Ea, const T& ref_value, double leaf_temp,
                               double H_d, double d_S) const {
  T arrh = arrh_curve<T>(Ea, ref_value, leaf_temp);
  double arg2 = 1 + std::exp((d_S*(25 + C_to_K) - H_d)/(gas_constant*(25 + C_to_K)));
  double arg3 = 1 + std::exp((d_S*(leaf_temp + C_to_K) - H_d)/(gas_constant*(leaf_temp + C_to_K)));

  return arrh * arg2/arg3;
}


// Every scalar update_temperature_dependent_params() reads: its cache key
// (hazard 10). ⚠️ Keep it beside that function and in step with it: an input
// missing here is served stale with no error; `photo_temp_key_size` makes a
// mismatched list a compile error. Compared by ==, so one ULP invalidates.
inline std::array<double, Leaf::photo_temp_key_size> Leaf::photo_temp_key() const {
  return {leaf_temp_, atm_o2_kpa_,
          vcmax_25, vcmax_ha_, vcmax_H_d_, vcmax_d_S_,
          jmax_25, jmax_ha_, jmax_H_d_, jmax_d_S_,
          gamma_25_, gamma_ha_,
          kc_25_, kc_ha_,
          ko_25_, ko_ha_,
          R_d_25, rd_q10_intercept_, rd_q10_slope_,
          // The thermal cost reaches jmax_; a bool widens to 0.0/1.0.
          double(use_thermal_cost_), T50_, Tcrit_,
          // The thermal-damage gate is here and its inputs (A_, I_, D_) are not:
          // set_physiology bypasses this cache with the gate on; drop that bypass
          // and every input the block reads must join this key. Without the gate,
          // switching it off would be served the gate-on block.
          double(use_thermal_damage_)};
}

// The temperature-dependent photosynthetic parameters at `leaf_temp`, in the
// operation order the golden files pin. electron_transport_ is recomputed from
// the new jmax_ and the call's PPFD_ (set_physiology refreshes it on a cache
// hit). R_d comes from R_d_25 and a declining Q10, so it rises with temperature.
inline void Leaf::update_temperature_dependent_params(double leaf_temp) {
  // First: two kernels called below read the block's temperature from here.
  photo_temp_ = leaf_temp;
  if (use_thermal_damage_) {
    // The ATLS layer (see use_thermal_damage_). One shift moves both capacity
    // optima; the impairment has its own. A itself, for delta_impair's reason.
    const double A = A_;
    const double delta_topt = topt_offset_ + dTopt_max_ * A / (A + K_A_);
    const double vcmax_dS =
        d_S_shifted(vcmax_ha_, vcmax_H_d_, vcmax_d_S_, delta_topt);
    const double jmax_dS = d_S_shifted(jmax_ha_, jmax_H_d_, jmax_d_S_, delta_topt);
    vcmax_ = peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_dS);
    jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_dS);
    // Only the functional area F photosynthesises, so both capacities scale by
    // F; reversible deactivation is already in the peaked Arrhenius. F is held
    // above zero by a softplus floor of width damage_surv_floor, exact in F where
    // F > ~3.4e-5: a min/max floor puts a corner in every output.
    {
      const double F = 1.0 - I_ - D_;
      const double f = damage_surv_floor;
      functional_ = F > 0.0 ? F + f * std::log1p(std::exp(-F / f))
                            : f * std::log1p(std::exp(F / f));
    }
    vcmax_ *= functional_;
    jmax_ *= functional_;
    k_imp_eff_ = k_imp_eff(leaf_temp, delta_impair_at(A));
    k_rec_eff_ = k_rec_eff(leaf_temp);
  } else {
    vcmax_ =
        peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_d_S_);
    jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_d_S_);
    k_imp_eff_ = 0.0;
    k_rec_eff_ = 0.0;
    functional_ = 1.0;
  }
  // Sicangco et al. (2026) Eqn 12: PSII damage scales Jmax by (1 - TC).
  // Branched, not multiplied by a gate-off zero, so gate-off arithmetic is untouched.
  if (use_thermal_cost_) {
    jmax_ *= 1.0 - thermal_cost_at(leaf_temp);
  }
  gamma_ = arrh_curve(gamma_ha_, gamma_25_, leaf_temp);
  ko_ = arrh_curve(ko_ha_, ko_25_, leaf_temp);
  kc_ = arrh_curve(kc_ha_, kc_25_, leaf_temp);
  // Exactly R_d_25 at 25 C, so a check at 25 C alone is blind to this response;
  // test_rd_temperature_response covers it. An unset R_d_25 is a caller error,
  // not a cue to derive one from vcmax.
  if (!std::isfinite(R_d_25) || R_d_25 < 0.0) {
    util::stop("R_d_25 must be a finite, non-negative dark respiration at 25 C; "
               "got " + util::format_double(R_d_25));
  }
  const double q10 =
      rd_q10_intercept_ - rd_q10_slope_ * (leaf_temp + 25.0) / 2.0;
  R_d_ = R_d_25 * std::pow(q10, (leaf_temp - 25.0) / 10.0);
  if (use_thermal_damage_) {
    // The layer's maintenance costs, added to R_d rather than folded into it so
    // gate-off arithmetic is untouched: standing costs for acclimation, repair
    // capacity and protection, and an activity cost for repair done. Protection
    // is signed: the default k_i_ equals k_i_baseline_, so a max(0, .) would put
    // the defaults on a corner.
    const double acclim_load = A_;
    const double protection = k_i_baseline_ - k_i_;
    const double repair_flux = k_rec_eff_ * I_;
    R_d_ += c_acclim_maint_ * acclim_load + c_repair_maint_ * k_rec_ +
            c_protect_maint_ * protection + c_repair_flux_ * repair_flux;
  }
  km_ = (kc_*umol_per_mol_to_Pa_)*(1 + (atm_o2_kpa_*kPa_to_Pa)/(ko_*umol_per_mol_to_Pa_));
  electron_transport_ = electron_transport();
}

// --- the ATLS thermal-damage layer ------------------------------------------


inline double Leaf::softplus_(double x, double s) {
  return damage::softplus(x, s);
}

// Inverts the peaked Arrhenius optimum Topt = H_d/(d_S - corr), corr =
// R*ln(Ea/(H_d - Ea)), for the d_S that puts it delta_topt higher. The optimum
// exists only for 0 < Ea < H_d, so outside it is refused, not a NaN capacity.
inline double Leaf::d_S_shifted(double Ea, double H_d, double d_S_base,
                                double delta_topt) const {
  if (!(Ea > 0.0 && H_d > Ea)) {
    util::stop("the thermal-damage layer moves the peaked Arrhenius optimum, "
               "which exists only for 0 < Ea < H_d; got Ea = " +
               util::format_double(Ea) + ", H_d = " + util::format_double(H_d));
  }
  const double corr = gas_constant * std::log(Ea / (H_d - Ea));
  const double Topt0_K = H_d / (d_S_base - corr);
  // Refused outside [250, 350] K rather than clamped there, as the impairment's
  // isokinetic temperature is (damage::k_imp).
  const double Topt_K = Topt0_K + delta_topt;
  if (Topt_K < 250.0 || Topt_K > 350.0) {
    util::stop("the thermal-damage layer moves the peaked Arrhenius optimum to " +
               util::format_double(Topt_K) + " K, outside [250, 350] K; "
               "topt_offset_ + dTopt_max_ * A/(A + K_A_) = " +
               util::format_double(delta_topt) + " K");
  }
  return H_d / Topt_K + corr;
}

// The layer's helpers on this leaf: its parameters, packed for the templates in
// `damage` above.
inline damage::dyn_pars<double> Leaf::damage_pars() const {
  damage::dyn_pars<double> p{};
  p[damage::dp_k_i] = k_i_;
  p[damage::dp_k_i_Ea] = k_i_Ea_;
  p[damage::dp_k_i_Tiso] = k_i_Tiso_;
  p[damage::dp_k_rec] = k_rec_;
  p[damage::dp_k_rec_Ea] = k_rec_Ea_;
  p[damage::dp_k_mat] = k_mat_;
  p[damage::dp_k_mat_Ea] = k_mat_Ea_;
  p[damage::dp_acclim_alpha] = acclim_alpha_;
  p[damage::dp_acclim_beta] = acclim_beta_;
  p[damage::dp_acclim_T] = acclim_T_;
  p[damage::dp_acclim_s] = acclim_s_;
  p[damage::dp_dTimpair_max] = dTimpair_max_;
  p[damage::dp_K_A] = K_A_;
  p[damage::dp_topt_offset] = topt_offset_;
  return p;
}

// A negative k_i_Ea_ is refused: impairment would slow as the leaf heats. Zero
// is allowed. A NaN passes, poisoning the outputs where a test can see it.
inline void Leaf::check_damage_pars() const {
  if (k_i_Ea_ < 0.0) {
    util::stop("the thermal-damage layer needs k_i_Ea_ >= 0, so that "
               "impairment does not slow with temperature; got k_i_Ea_ = " +
               util::format_double(k_i_Ea_));
  }
}

// Its own parameters only: reading jmax_H_d_ or jmax_d_S_ here would make
// recalibrating photosynthesis silently move the impairment rate.
inline double Leaf::k_imp_eff(double leaf_temp, double delta_impair) const {
  check_damage_pars();
  return damage::k_imp(leaf_temp, delta_impair, damage_pars());
}

inline double Leaf::delta_impair_at(double A) const {
  return damage::delta_impair(A, damage_pars());
}

inline double Leaf::k_rec_eff(double leaf_temp) const {
  return damage::k_rec_eff(leaf_temp, damage_pars());
}

inline std::vector<double> Leaf::damage_rates(double leaf_temp, double A,
                                              double I, double D) const {
  check_damage_pars();
  if (A <= -K_A_) {
    util::stop("damage_rates: the acclimation shift A/(A + K_A_) has its pole "
               "at A = -K_A_; got A = " + util::format_double(A) +
               ", K_A_ = " + util::format_double(K_A_));
  }
  const std::array<double, 3> r = damage::rates(leaf_temp, A, I, D, damage_pars());
  return {r[0], r[1], r[2]};
}

inline void Leaf::check_damage_initial_states(double A0, double I0,
                                              double D0) const {
  if (I0 < 0.0 || D0 < 0.0 || A0 <= -K_A_) {
    util::stop("damage_trajectory: the initial states need I0 >= 0, D0 >= 0 "
               "and A0 > -K_A_; got A0 = " + util::format_double(A0) +
               ", I0 = " + util::format_double(I0) + ", D0 = " +
               util::format_double(D0));
  }
}

inline void Leaf::check_damage_series(const std::vector<double>& times,
                                      const std::vector<double>& leaf_temp) const {
  const std::size_t n = times.size();
  if (n == 0 || leaf_temp.size() != n) {
    util::stop("damage_trajectory: times and leaf_temp must be non-empty and "
               "the same length; got " + std::to_string(n) + " and " +
               std::to_string(leaf_temp.size()));
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(times[i]) || !std::isfinite(leaf_temp[i]) ||
        (i > 0 && !(times[i] > times[i - 1]))) {
      util::stop("damage_trajectory: times must be finite and strictly "
                 "increasing, and leaf_temp finite; failed at index " +
                 std::to_string(i));
    }
  }
  check_damage_pars();
}

inline std::vector<double> Leaf::damage_trajectory(
    const std::vector<double>& times, const std::vector<double>& leaf_temp,
    double A0, double I0, double D0) const {
  check_damage_series(times, leaf_temp);
  check_damage_initial_states(A0, I0, D0);
  const std::size_t n = times.size();
  std::vector<double> out(3 * n);
  damage::trajectory<double>(times, leaf_temp, {A0, I0, D0}, damage_pars(),
                             [&](std::size_t i, const std::array<double, 3>& y) {
                               out[i] = y[0];
                               out[n + i] = y[1];
                               out[2 * n + i] = y[2];
                             });
  return out;
}

inline std::vector<std::string> Leaf::damage_par_names() const {
  return {"k_i_", "k_i_Ea_", "k_i_Tiso_", "k_rec_", "k_rec_Ea_", "k_mat_",
          "k_mat_Ea_", "acclim_alpha_", "acclim_beta_", "acclim_T_", "acclim_s_",
          "dTimpair_max_", "K_A_", "topt_offset_", "A0", "I0", "D0"};
}

// One forward-mode pass per requested input (a parameter or initial state)
// through the double trajectory's RK4: exact for the discrete map, at the
// double run's sub-step count.
inline std::vector<double> Leaf::damage_trajectory_gradient(
    const std::vector<double>& times, const std::vector<double>& leaf_temp,
    double A0, double I0, double D0, const std::vector<int>& which) const {
  check_damage_series(times, leaf_temp);
  check_damage_initial_states(A0, I0, D0);
  const std::size_t n = times.size(), np = which.size();
  const int n_inputs = damage::n_dyn_pars + 3;
  std::vector<double> out(3 * n * np);
  const damage::dyn_pars<double> base = damage_pars();
  for (std::size_t k = 0; k < np; ++k) {
    const int w = which[k];
    if (w < 0 || w >= n_inputs) {
      util::stop("damage_trajectory_gradient: input index " + std::to_string(w) +
                 " is out of range; see damage_par_names()");
    }
    damage::dyn_pars<damage::Dual> p;
    for (int j = 0; j < damage::n_dyn_pars; ++j) {
      p[std::size_t(j)] = damage::Dual(base[std::size_t(j)], j == w ? 1.0 : 0.0);
    }
    const int s0 = w - damage::n_dyn_pars;
    std::array<damage::Dual, 3> y0{damage::Dual(A0, s0 == 0 ? 1.0 : 0.0),
                                   damage::Dual(I0, s0 == 1 ? 1.0 : 0.0),
                                   damage::Dual(D0, s0 == 2 ? 1.0 : 0.0)};
    damage::trajectory<damage::Dual>(
        times, leaf_temp, y0, p,
        [&](std::size_t i, const std::array<damage::Dual, 3>& y) {
          for (int st = 0; st < 3; ++st) {
            out[i + n * (std::size_t(st) + 3 * k)] = y[std::size_t(st)].d;
          }
        });
  }
  return out;
}

inline void Leaf::refuse_thermal_damage(const char* who) const {
  if (use_thermal_damage_) {
    util::stop(std::string(who) +
               ": refused with use_thermal_damage_ on. The derivative surface "
               "takes Vcmax, Jmax and R_d from the traits through a temperature "
               "response that does not apply the thermal-damage layer -- its "
               "shifted optimum, its functional fraction and its maintenance costs -- "
               "so its rows would describe a different model from the one the "
               "solve ran. Set use_thermal_damage_ to false.");
  }
}

// Saturation vapour pressure es(T) in kPa (Tetens), T in deg C.
inline double Leaf::saturation_vapour_pressure(double temp) const {
  return 0.6108 * exp(17.27 * temp / (temp + 237.3));
}

// Slope of the saturation vapour pressure curve Delta(T) in kPa K^-1, T in deg C:
// the exact derivative of the Tetens form above. Not FAO-56's rounded 4098
// (4.2e-5 low), which reaches the energy-balance first-order condition and put
// dprofit/dpsi_stem 2e-6 to 2e-5 off a central difference of profit.
inline double Leaf::saturation_vapour_pressure_slope(double temp) const {
  return (17.27 * 237.3) * saturation_vapour_pressure(temp) /
         ((temp + 237.3) * (temp + 237.3));
}

// The deficit driving diffusion out of the leaf, from the cavity saturated at
// Tleaf to the air (the driver `atm_vpd_` is at Tair):
//
//   vpd_leaf = esat(Tleaf) - e_air = atm_vpd + [esat(Tleaf) - esat(Tair)]
//
// In that association, so off the energy-balance path the bracket is exactly
// +0.0 and the result bit-identical. Held at `vpd_leaf_min`, a kept exception to
// the no-clamp rule: below it the leaf condenses, and a conductance formed from
// the deficit would be infinite or negative.
inline double Leaf::leaf_vpd_at(double leaf_temp) const {
  if (!use_energy_balance_) {
    return atm_vpd_;
  }
  const double d =
      atm_vpd_ + (saturation_vapour_pressure(leaf_temp) -
                  saturation_vapour_pressure(Tair_));
  return std::max(d, vpd_leaf_min);
}

// Also seats the residual conductance at the same temperature: it is read
// beside the deficit wherever a conductance is formed, and the two must
// belong to one leaf temperature.
inline void Leaf::set_leaf_vpd(double leaf_temp) {
  vpd_leaf_ = leaf_vpd_at(leaf_temp);
  g_min_leaf_ = g_min_at(leaf_temp);
}

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

// The intercellular CO2 at which supply meets demand,
// `A(ci) umol_to_mol = gc (ca - ci) / P`, root-found on [Gamma* P, hi]. The
// target is negative at Gamma* (A = -R_d) and rises with ci. With A(ca) < 0 the
// leaf is a net source and the root lies above ca, below
// 2 |A(ca)| umol_to_mol P / gc above it, where hi is placed. With no conductance
// the residual is A(ci) = 0: its root where one exists (the g_min -> 0 limit of
// the shut leaf), else ci = +inf.
//
// ci_root_tol (1e-13) sets the floor of what every output means; it is not the
// settable `ci_abs_tol`, which reaches only the Medlyn route. Worst golden-grid
// difference from a converged solve, by TOMS748 tolerance (the plateaus are its
// discrete iterations):
//
//   1e-7: 1.82e-07   1e-8: 3.15e-08   1e-10, 1e-11: 5.41e-10   1e-13: 7.48e-13
//
// The Newton solve below stops at rounding; TOMS748 is its fallback. Why 1e-13
// and not 1e-10, and why not to loosen it: .claude/measurements.md, "Magnitudes:
// what a difference means".
inline double Leaf::solve_ci(double gc) {
  // NA in, NA out, without entering the solver.
  if (!std::isfinite(gc)) {
    return util::na_value;
  }
  auto target = [&](double x) mutable -> double {
    const double assim_colimited_x_ = assim_colimited(x);
    return assim_colimited_x_ * umol_to_mol -
      (gc * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
  };
  const double lo = gamma_ * umol_per_mol_to_Pa_;
  const double A_ca = assim_colimited(ca_);
  double hi = ca_;
  // The target at ca is A(ca) exactly (the supply term is gc * 0), handed over
  // rather than re-evaluated: bit-identical, since x - 0.0 is x.
  if (A_ca < 0.0) {
    if (gc > 0.0) {
      hi = ca_ + 2.0 * (-A_ca) * umol_to_mol * (atm_kpa_ * kPa_to_Pa) / gc;
    } else {
      // No conductance at all: A(ci) = 0, or nothing nets zero.
      if (!(assim_saturated() > 0.0)) {
        return std::numeric_limits<double>::infinity();
      }
      // A rises with ci towards a positive limit, so doubling reaches it.
      do {
        hi *= 2.0;
      } while (assim_colimited(hi) <= 0.0);
    }
  }
  try {
    const double f_hi = hi == ca_ ? A_ca * umol_to_mol -
                                        (gc * (ca_ - hi) / (atm_kpa_ * kPa_to_Pa))
                                  : target(hi);
    const double f_lo = target(lo);
    // Safeguarded Newton first. The residual is increasing and concave in ci, so
    // Newton from below the root climbs to it and from above overshoots at most
    // once; a step leaving the bracket is replaced by the Newton step from its
    // lower end, or a bisection where that slope is not finite. Stopping on a
    // step below ci_root_tol (1 + |ci|) places the root to rounding, not in a
    // band, so there is no state-dependent noise for a caller to difference.
    if (f_lo < 0.0 && f_hi >= 0.0) {
      if (f_hi == 0.0) {
        return hi;
      }
      const PhotoCapacity<double> cap = photo_capacity<double>();
      const double supply_per_pa = gc / (atm_kpa_ * kPa_to_Pa);
      // The value is `target`'s and the slope is the kernel's own tangent
      // (assim_slope_at), so the Newton root and TOMS748's are the same
      // function's to the last bit and the slope cannot drift from it.
      auto residual = [&](double x, double& slope) -> double {
        slope = assim_slope_at<double>(x, cap) * umol_to_mol + supply_per_pa;
        return target(x);
      };
      double a = lo, b = hi;
      double f_a = f_lo, slope_a = std::numeric_limits<double>::quiet_NaN();
      double x = lo - f_lo * (hi - lo) / (f_hi - f_lo);
      for (int i = 0; i < static_cast<int>(ci_niter); ++i) {
        double slope;
        const double f = residual(x, slope);
        if (f == 0.0) {
          return x;
        }
        if (f < 0.0) {
          a = x;
          f_a = f;
          slope_a = slope;
        } else {
          b = x;
        }
        const double newton = x - f / slope;
        // Convergence is tested before the safeguard: at the root `newton` lands
        // on the bracket end x just became, which the safeguard would bisect.
        if (std::abs(newton - x) <= ci_root_tol * (1.0 + std::abs(x))) {
          return newton;
        }
        if (newton > a && newton < b) {
          x = newton;
          continue;
        }
        // Overshot from above: the step from a, since hi can sit orders of
        // magnitude above the root (a net source with almost no conductance).
        // Where a's slope is not finite (Gamma*) a wide bracket splits geometrically.
        if (!std::isfinite(slope_a)) {
          residual(a, slope_a);
        }
        const double from_a = a - f_a / slope_a;
        x = (from_a > a && from_a < b)   ? from_a
            : (a > 0.0 && b > 4.0 * a) ? std::sqrt(a * b)
                                       : 0.5 * (a + b);
      }
    }
    return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, ci_root_tol,
                                ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("ci_solve", "solve_ci failed: " + std::string(e.what()) +
               "; min=" + util::to_string(lo) +
               "; max=" + util::to_string(hi) +
               "; gc=" + util::to_string(gc));
  }
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

// Explicit leaf energy balance: Tleaf = Tair + (Rn - lambda*E) / (rho*cp/ra + g_rad),
// with E the hydraulically pinned transpiration (kg H2O m^-2 s^-1). The one place
// Tleaf and dTleaf/dE are computed: a second copy is how a radiative term gets
// dropped from a derivative while the value keeps it.
inline double Leaf::leaf_temp_from_E(double E, double* dT_dE) const {
  // Sensible heat and the linearised longwave as one summed conductance, so the
  // radiative term cannot be dropped by editing one factor.
  const double g_total = vol_heat_cap_air / ra_ + g_rad_;
  const double Tleaf = Tair_ + (Rn_ - latent_heat_vap * E) / g_total;
  if (dT_dE != nullptr) {
    // Linear in E, so the slope is constant and negative.
    *dT_dE = -latent_heat_vap / g_total;
  }
  return Tleaf;
}

// Refuses, never clamps, a leaf temperature range the model cannot evaluate:
// R_d's Q10, intercept - slope*(T + 25)/2, must be positive for its power to be
// real, and the Tetens form has a pole at -237.3 C. Q10 is linear in T, so
// testing both ends tests the range.
inline void Leaf::check_leaf_temperature_range(double coldest,
                                               double hottest) const {
  const auto q10_at = [&](double t) {
    return rd_q10_intercept_ - rd_q10_slope_ * (t + 25.0) / 2.0;
  };
  if (!(coldest > -237.3) || !(q10_at(coldest) > 0.0) ||
      !(q10_at(hottest) > 0.0)) {
    util::stop("set_physiology: the leaf temperature can run from " +
               util::format_double(coldest) + " to " +
               util::format_double(hottest) +
               " C over the fluxes this solve can reach, and the model is "
               "evaluable only above the Tetens pole (-237.3 C) and where R_d's "
               "Q10, rd_q10_intercept_ - rd_q10_slope_*(T + 25)/2, is positive "
               "(Q10 " + util::format_double(q10_at(coldest)) + " at the cold "
               "end, " + util::format_double(q10_at(hottest)) + " at the hot)");
  }
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

// The TOTAL conductance to CO2, mol CO2 m^-2 LA s^-1: the stomata's, from the
// flux the residual leaves them, plus the cuticle's. Reads the deficit as seated
// (set_leaf_vpd), so on the energy-balance path seat the leaf temperature first.
inline double Leaf:: stom_cond_CO2(double psi_stem, double psi_upstream) {
  const double E = transpiration(psi_stem, psi_upstream);
  double E_s = E - residual_flux(psi_upstream, E);
  if (E_s < 0.0) {
    clamps.note(CLAMP_STOMATAL_FLUX);
    E_s = 0.0;
  }
  const double gc =
      atm_kpa_ * E_s * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_ +
      cuticular_cond_CO2();
  if (use_energy_balance_) [[unlikely]] {  // see dprofit_at_collar_psi
    return boundary_layer_split_CO2(gc).first;
  }
  return gc;
}


// Electron transport: the non-rectangular hyperbola in phi_J * APAR, APAR =
// alpha_PAR_ * PPFD_ (intercepted photon flux, umol m^-2 s^-1).
// electron_transport_at is the same expression in the same order.
inline double Leaf::electron_transport() {
  const double I = phi_J * alpha_PAR_ * PPFD_;
  double electron_transport_ = (I + jmax_ - sqrt(pow(I + jmax_, 2) -
  4 * curv_fact_elec_trans * I * jmax_)) / (2 * curv_fact_elec_trans);
    return electron_transport_;
}

// The scalar-generic photosynthesis cores, umol m^-2 s^-1; the golden files pin
// the double instantiation's association. U is independent of T on purpose: a
// double times a tangent is one multiply, a tangent times a tangent three, so
// wrapping the capacity in T measured slower per solve.
template <typename T, typename U>
inline T Leaf::assim_rubisco_limited_kernel(const T& ci, const U& vcmax) const {
  return (vcmax * (ci - gamma_ * umol_per_mol_to_Pa_)) / (ci + km_);
}

template <typename T>
inline T Leaf::assim_rubisco_limited_kernel(T ci) const {
  return assim_rubisco_limited_kernel<T, double>(ci, vcmax_);
}

template <typename T, typename U>
inline T Leaf::assim_electron_limited_kernel(const T& ci,
                                             const U& transport) const {
  return transport / 4 *
  ((ci - gamma_ * umol_per_mol_to_Pa_) / (ci + 2 * gamma_ * umol_per_mol_to_Pa_));
}

template <typename T>
inline T Leaf::assim_electron_limited_kernel(T ci) const {
  return assim_electron_limited_kernel<T, double>(ci, electron_transport_);
}

template <typename T, typename U>
inline T Leaf::assim_colimited_kernel(const T& ci,
                                      const PhotoCapacity<U>& cap) const {
  const T rubisco = assim_rubisco_limited_kernel<T, U>(ci, cap.vcmax);
  const T electron = assim_electron_limited_kernel<T, U>(ci, cap.transport);

  return (rubisco + electron -
          sqrt(pow(rubisco + electron, 2) -
               4 * cap.curvature * rubisco * electron)) /
             (2 * cap.curvature) -
         cap.respiration;
}

// The colimited form with both limitations at their c_i -> infinity values:
// Rubisco's at vcmax, the electron-limited one at J/4. The same association
// as assim_colimited_kernel, so the two agree where both are finite.
template <typename T, typename U>
inline T Leaf::assim_saturated_kernel(const PhotoCapacity<U>& cap) const {
  const T rubisco = T(cap.vcmax);
  const T electron = T(cap.transport / 4);
  return (rubisco + electron -
          sqrt(pow(rubisco + electron, 2) -
               4 * cap.curvature * rubisco * electron)) /
             (2 * cap.curvature) -
         cap.respiration;
}

inline double Leaf::assim_saturated() const {
  return assim_saturated_kernel<double, double>(photo_capacity<double>());
}

template <typename T>
inline T Leaf::assim_colimited_kernel(T ci) const {
  return assim_colimited_kernel<T, double>(ci, photo_capacity<double>());
}

inline double Leaf::assim_rubisco_limited(double ci_) {
  return assim_rubisco_limited_kernel(ci_);
}

//calculate the light-limited assimilation rate, returns umol m^-2 s^-1
inline double Leaf::assim_electron_limited(double ci_) {
  return assim_electron_limited_kernel(ci_);
}

// returns co-limited assimilation umol m^-2 s^-1, NET of dark respiration (the
// trailing `- R_d_`), so gross assimilation is this value + R_d_.
inline double Leaf::assim_colimited(double ci_) {
  return assim_colimited_kernel(ci_);
}


// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
inline double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream) {

  double assim_colimited_x_ = assim_colimited(x);

  double stom_cond_CO2_x_ = stom_cond_CO2(psi_stem, psi_upstream);
  return assim_colimited_x_ * umol_to_mol -
         (stom_cond_CO2_x_ * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
}

// The ci at which A(ci) = gc(ca - ci). The target is smooth and strictly
// monotone in ci (Ar and Ae vanish linearly at gamma*), so a superlinear
// bracketing solver suits it; the bracket and the degenerate case are solve_ci's.
inline double Leaf::psi_stem_to_ci(double psi_stem, double psi_upstream) {
  return ci_ = solve_ci(stom_cond_CO2(psi_stem, psi_upstream));
}

// The operating point at psi_stem in one pass, with no branch on the flow: the
// flux (zero for a reversed pair), the leaf temperature, the residual the cuticle
// loses, and ci on the total conductance. A shut leaf is the same expressions
// with the stomatal term zero; with no conductance at all solve_ci takes the
// limit. Water moves whatever photosynthesis does: no exit on assim_max_ < 0.
// ⚠️ Tleaf_ is written on every path (hazard 8): with the energy balance on, a
// leaf transpiring nothing is the hottest. The temperature block and the deficit
// are seated at the stored Tleaf_, so the three temperatures are one number.
inline void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream) {
  transpiration_ = operating_point_flux(psi_stem, psi_upstream);
  Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(transpiration_) : leaf_temp_;
  if (use_energy_balance_) {
    update_temperature_dependent_params(Tleaf_);
    set_leaf_vpd(Tleaf_);
  }
  residual_flux_ = residual_flux(psi_upstream, transpiration_);
  // Below the wet bound the stomatal flux reads as zero, counted, so a collar
  // that inverts the pair by rounding is the wet bound, not a refusal.
  double E_s = transpiration_ - residual_flux_;
  if (E_s < 0.0) {
    clamps.note(CLAMP_STOMATAL_FLUX);
    E_s = 0.0;
  }
  stom_cond_CO2_ = atm_kpa_ * E_s * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_ +
                   cuticular_cond_CO2();
  if (use_energy_balance_) [[unlikely]] {  // see dprofit_at_collar_psi
    stom_cond_CO2_ = boundary_layer_split_CO2(stom_cond_CO2_).first;
  }
  ci_ = solve_ci(stom_cond_CO2_);
  // With no conductance at all A is zero by definition, not the root's ~1e-15
  // residue, which the product and ratio curves' lambda turn into +-inf and NaN.
  // On the tape this A takes its rows through the conductance's g_min row.
  assim_colimited_ = !std::isfinite(ci_)          ? assim_saturated()
                     : stom_cond_CO2_ == 0.0      ? 0.0
                                                  : assim_colimited(ci_);
}


// Hydraulic cost equations

// --- Marginal cost of water -------------------------------------------------

inline double Leaf::lambda_TF24(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  return TF24_cost_scale * TF24_beta2 * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0) *
         pow(1.0 - f, TF24_beta2 - 1.0) / leaf_specific_conductance_max_;
}


// JS22: dC/dpsi = 2*gamma*(psi - psi_up) over dE/dpsi = kmax*f(psi). Unlike
// TF24's, f does not cancel, so this divides by it: safe inside the feasible
// bracket, where f never reaches zero.
inline double Leaf::lambda_JS22(double psi_stem, double psi_upstream) const {
  const double f = proportion_of_conductivity(psi_stem);
  return 2.0 * JS22_gamma * (psi_stem - psi_upstream) /
         (leaf_specific_conductance_max_ * f);
}


// CMax: the marginal cost does not vanish at the wet end unless CMax_b is zero,
// so it is the one curve that can be boundary-soil with a non-zero cost there.
inline double Leaf::lambda_CMax(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double f = proportion_of_conductivity(psi_stem);
  return (CMax_a * psi_stem + CMax_b) /
         (leaf_specific_conductance_max_ * f);
}


// TF24_floor: TF24's own cost plus a linear price of water, so lambda is
// lambda_TF24 + lambda_o. Called, not repeated: at lambda_o == 0 it is
// lambda_TF24 bit for bit, since that ends in a division no FMA can absorb.
inline double Leaf::lambda_TF24_floor(double psi_stem) const {
  return lambda_TF24(psi_stem) + TF24_floor_lambda_o;
}


// SOX's objective is a product: at its optimum A'g + Ag' = 0, so
//
//     lambda = A'/E' = -A*(g'/g) / (kmax*f)
//
// It carries a factor of A, since a product prices water as a share of
// assimilation (the-models.Rmd), and diverges at psi_crit, where g = 0: never
// pinned there, and not a number to evaluate at the dry end.
inline double Leaf::lambda_SOX(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double g = sox_reduction(psi_stem);
  const double gp = sox_reduction_deriv(psi_stem);
  const double f = proportion_of_conductivity(psi_stem);
  return -assim_colimited_ * (gp / g) / (leaf_specific_conductance_max_ * f);
}


// Jones: g = 1 - psi/psi_crit makes -g'/g = 1/(psi_crit - psi), so lambda =
// A / [(psi_crit - psi) kmax f(psi)], as the-models.Rmd quotes; it carries A and
// diverges at psi_crit, as lambda_SOX's does.
inline double Leaf::lambda_JW26(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double f = proportion_of_conductivity(psi_stem);
  return assim_colimited_ /
         ((psi_crit - psi_stem) * leaf_specific_conductance_max_ * f);
}


// ProfitMax: cost (k_soil - kmax f)/k_span + TC(Tleaf), in carbon units by |A|max:
//
//     lambda = |A|max * [ |f'|/(f*k_span) + (dTC/dT)(dT/dE) ]
//
// The thermal term vanishes with either gate off.
inline double Leaf::lambda_ProfitMax(double psi_stem) {
  const double f = proportion_of_conductivity(psi_stem);
  if (!(f > 0.0) || !(profitmax_k_span_ > 0.0) || !(profitmax_A_max_ > 0.0)) {
    return util::na_value;
  }
  // |f'| = (c/b)(psi/b)^(c-1) f, positive because f is decreasing.
  const double abs_fprime =
      (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0) * f;

  double thermal = 0.0;
  if (use_thermal_cost_ && use_energy_balance_) {
    const double E = operating_point_flux(psi_stem, supply_psi_soil_scalar());
    double dT_dE = 0.0;
    const double Tleaf = leaf_temp_from_E(E, &dT_dE);
    // The same slope dprofit_energy_balance_term subtracts from dJ/dpsi, so the
    // reported lambda and the first-order condition cannot disagree about it.
    thermal = thermal_cost_temp_deriv(Tleaf) * dT_dE;
  }
  return profitmax_A_max_ *
         (abs_fprime / (f * profitmax_k_span_) + thermal);
}

inline double Leaf::marginal_cost_water() const {
  return lambda_TF24(opt_psi_stem_);
}

inline double Leaf::marginal_cost_water_molar() const {
  // umol CO2 (kg H2O)^-1 -> mol CO2 (mol H2O)^-1
  return marginal_cost_water() * umol_to_mol / kg_to_mol_h2o;
}

inline double Leaf::marginal_cost_water_multilayer() {
  // S is the soil->collar conductance, positive by construction (hazard 2).
  const double S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
  if (!std::isfinite(S) || S <= 0.0) {
    return util::na_value;
  }
  // f is the STEM vulnerability curve at the collar suction: kmax*f(psi_r) is
  // dE_stem/dpsi_r, the stem-side conductance where the two paths meet.
  const double f_r = proportion_of_conductivity(opt_root_psi_);
  return lambda_TF24(opt_psi_stem_) *
         (1.0 + leaf_specific_conductance_max_ * f_r / S);
}

inline double Leaf::g1_eff() const {
  const double chi = ci_ / ca_;
  if (!std::isfinite(chi) || chi >= 1.0) {
    return util::na_value;
  }
  return chi * std::sqrt(vpd_leaf_) / (1.0 - chi);
}

// Pure: no write to hydraulic_cost_, so the AD pass cannot scribble model state
// while probing. The caching is the double entry point's job.
template <typename T, typename U>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem, const U& scale,
                                        const U& b, const U& c,
                                        const U& beta2) const {
  return scale *
         pow((1 - proportion_of_conductivity_kernel<T, U>(psi_stem, b, c)), beta2);
}

template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(T psi_stem) const {
  return hydraulic_cost_TF_kernel<T, double>(psi_stem, TF24_cost_scale, stem_b,
                                             stem_c, TF24_beta2);
}

// The kernels' arguments, read off a caller's parameter pack. At double they
// equal the members bit for bit (passive_pars builds the pack; responses are at
// photo_temp_), which test_leaf asserts on both temperature paths.

// These three responses do not apply the thermal-damage layer, so each refuses
// a leaf with it on rather than disagree with the member it twins.
template <typename T>
inline T Leaf::vcmax_at(const leaf_pars<T>& pars) const {
  refuse_thermal_damage("Leaf::vcmax_at");
  return peak_arrh_curve<T>(vcmax_ha_, pars[par_vcmax_25], photo_temp_,
                            vcmax_H_d_, vcmax_d_S_);
}

template <typename T>
inline T Leaf::jmax_at(const leaf_pars<T>& pars) const {
  refuse_thermal_damage("Leaf::jmax_at");
  T j = peak_arrh_curve<T>(jmax_ha_, pars[par_jmax_25], photo_temp_, jmax_H_d_,
                           jmax_d_S_);
  // The thermal cost reads no trait, so it is a factor, not a term with rows.
  if (use_thermal_cost_) {
    // Not `*=`, whose overload is ambiguous at a nested scalar.
    j = j * T(1.0 - thermal_cost_at(photo_temp_));
  }
  return j;
}

template <typename T>
inline T Leaf::respiration_at(const leaf_pars<T>& pars) const {
  refuse_thermal_damage("Leaf::respiration_at");
  const double q10 =
      rd_q10_intercept_ - rd_q10_slope_ * (photo_temp_ + 25.0) / 2.0;
  return pars[par_R_d_25] *
         T(std::pow(q10, (photo_temp_ - 25.0) / 10.0));
}

template <typename T>
inline T Leaf::electron_transport_at(const leaf_pars<T>& pars) const {
  // alpha_PAR_ is a leaf property, not a slot, so it enters as a constant.
  const T I = pars[par_phi_J] * alpha_PAR_ * pars[par_PPFD];
  const T jm = jmax_at<T>(pars);
  const T curv = pars[par_curv_fact_elec_trans];
  return (I + jm -
          sqrt(pow(I + jm, 2) - 4 * curv * I * jm)) /
         (2 * curv);
}

// gamma_ and km_ are temperature responses of fixed physiology -- no trait
// reaches either -- so they stay members and the capacity carries the four that
// do.
template <typename T>
inline Leaf::PhotoCapacity<T> Leaf::photo_capacity_at(
    const leaf_pars<T>& pars) const {
  return {vcmax_at<T>(pars), electron_transport_at<T>(pars),
          pars[par_curv_fact_colim], respiration_at<T>(pars)};
}

template <typename T>
inline Leaf::PhotoCapacity<T> Leaf::photo_capacity() const {
  return {T(vcmax_), T(electron_transport_), T(curv_fact_colim), T(R_d_)};
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

template <class S, class Kernel, std::same_as<S>... Args>
inline S Leaf::kernel_slope_at(Kernel&& kernel, const char* what, const S& x,
                               const Args&... args) const {
  using odelia::util::to_passive;
  // Both levels are tangents at double, so nothing here reaches a tape. The inner
  // direction is x's (the slope); the outer carries x, then each argument (rows).
  using inner = odelia::ode::tangent_scalar<double>;
  constexpr std::size_t n_dir = 1 + sizeof...(Args);
  using outer = typename xad::fwd<inner, n_dir>::active_type;

  // At double, the inner level only: there are no rows, and the outer walk would
  // cost n_dir directions on the path that places every operating point.
  if constexpr (std::is_same_v<S, double>) {
    inner xd(x);
    odelia::ode::set_direction(xd, 1.0);
    return odelia::ode::derivative_along(kernel(xd, inner(args)...));
  } else {
    std::array<outer, n_dir> lift;
    std::size_t d = 0;
    const auto set = [&](const S& v) {
      lift[d] = outer(to_passive(v));
      xad::value(xad::derivative(lift[d])[d]) = 1.0;
      ++d;
    };
    set(x);
    (set(args), ...);
    odelia::ode::set_direction(xad::value(lift[0]), 1.0);

    const outer y = std::apply(kernel, lift);
    const double slope = odelia::ode::derivative_along(xad::value(y));

    // A braced list is sequenced left to right, so the counter follows the walk;
    // a stack array because the rows cross as a span, which borrows its storage.
    std::size_t r = 0;
    const auto row = [&](const S& v) -> odelia::input_and_derivative<S> {
      return {v, odelia::ode::derivative_along(xad::derivative(y)[r++])};
    };
    const odelia::input_and_derivative<S> rows[n_dir]{row(x), row(args)...};
    S out;
    const odelia::supply_report report =
        odelia::supply_derivatives<S>(slope, rows, out);
    if (!report.whole) {
      util::stop(std::string("Leaf::kernel_slope_at: ") + what + ": " +
                 report.why);
    }
    return out;
  }
}

template <typename T>
inline T Leaf::assim_slope_at(const T& ci, const PhotoCapacity<T>& cap) const {
  return kernel_slope_at<T>(
      [this]<class U>(const U& c, const U& vcmax, const U& transport,
                      const U& curvature, const U& respiration) -> U {
        return assim_colimited_kernel<U>(
            c, PhotoCapacity<U>{vcmax, transport, curvature, respiration});
      },
      "assimilation in intercellular CO2", ci, cap.vcmax, cap.transport,
      cap.curvature, cap.respiration);
}

template <typename T>
inline T Leaf::cost_slope_at(const T& sigma, const leaf_pars<T>& pars) const {
  // Named rather than passed inline: a row borrows its input by reference, and a
  // temporary is not a thing to hand one.
  const T b =
      phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], pars[par_stem_c]);
  return kernel_slope_at<T>(
      [this]<class U>(const U& s, const U& scale, const U& bb, const U& c,
                      const U& beta2) -> U {
        return hydraulic_cost_TF_kernel<U, U>(s, scale, bb, c, beta2);
      },
      "hydraulic cost in stem potential", sigma, pars[par_TF24_cost_scale], b,
      pars[par_stem_c], pars[par_TF24_beta2]);
}

// The curve itself is closed form, so the (P50, c) chain comes out of the one
// derivation rather than a supplied row: only the pre-integrated G(psi) is tabulated and
// needs its rows supplied (stem_integral_at).
template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem,
                                        const leaf_pars<T>& pars) const {
  return hydraulic_cost_TF_kernel<T>(
      psi_stem, pars[par_TF24_cost_scale],
      phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], pars[par_stem_c]),
      pars[par_stem_c], pars[par_TF24_beta2]);
}

template <class S>
inline Leaf::CollarCoords<S> Leaf::collar_coords_at(
    double sigma_star, double ci_star, const S& collar,
    const SupplyDraw<S>& draw, const leaf_pars<S>& pars,
    const PhotoCapacity<S>& cap, bool collar_moves) const {
  S sigma_h, gc, ci_h;
  held_coords_at<S>(sigma_star, ci_star, collar, draw, pars, cap, sigma_h, gc,
                    ci_h);
  return collar_channel_at<S>(sigma_star, collar, draw, pars, cap, sigma_h, gc,
                              ci_h, collar_moves);
}

template <class S>
inline void Leaf::held_coords_at(double sigma_star, double ci_star,
                                 const S& collar, const SupplyDraw<S>& draw,
                                 const leaf_pars<S>& pars,
                                 const PhotoCapacity<S>& cap, S& sigma_h,
                                 S& gc, S& ci_h) const {
  using odelia::util::to_passive;
  check_draw(to_passive(collar), draw);

  // The residuals see a HELD collar and the collar's channel is added back as one
  // supplied slope (closed_form_rows.hpp): sigma is placed by an inverse table
  // whose slope is not 1/G'(sigma), so a live collar would carry two values of
  // dsigma/dcollar. The residuals carry the trait rows, the slopes the collar's.
  const S held(to_passive(collar));
  const double gc_per_flux =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const Leaf& leaf = *this;

  // dT1/dsigma from the table, the derivative the model forms, with kmax off the
  // pack the residual multiplies by: the slope of the residual as written.
  const double dT1_dsigma =
      to_passive(pars[par_kmax]) * stem_curve_integral_deriv(sigma_star);
  // ⚠️ Every active input the residual reads must be in implicit_value's list:
  // the residual is rewound, so one left out is a row that reads as an exact
  // zero (hazard 18). `held` is passive and carries none.
  // At a dry pin the stem is held at the curve's end, a function of the traits:
  // T1's slope there is kmax e^-40, and dividing by it amplifies rounding. A shut
  // leaf's stem seated at its collar IS the collar (V = 1 in the channel):
  // inverting T1 there divides by a conductivity of ~1e-20.
  if (stem_held_at_dry_pin(sigma_star)) {
    sigma_h = stem_curve_end_at<S>(pars);
  } else if (stem_held_at_collar(sigma_star, to_passive(collar))) {
    sigma_h = held;
  } else {
    sigma_h = odelia::implicit_value<S>(
        sigma_star, dT1_dsigma,
        [&](const S& sg) -> S {
          return leaf.template transpiration_at<S>(sg, held, pars) -
                 draw.E_up.value;
        },
        pars, draw.E_up.value);
  }

  // gc is built from the stem's flux at the placed sigma, which the ci solve
  // reads, not from the draw (they agree only to the table's round-trip). At
  // g_min = 0 the residual and the cuticle are exact zeros carrying the g_min row.
  const S stem_flux = transpiration_at<S>(sigma_h, held, pars);
  const S r_star = residual_flux_at<S>(held, pars, draw);
  const S cuticle = pars[par_g_min] * S(g_min_co2_ratio_ * g_min_factor(Tleaf_));
  // Where nothing conducts (a shut leaf at g_min = 0, on its wet bound) the
  // stomatal term carries no row: the collar's channel is skipped there, so
  // nothing would cancel the held collar's partial in stem_flux.
  gc = stom_cond_CO2_ == 0.0
           ? cuticle
           : S(gc_per_flux) * (stem_flux - r_star) + cuticle;
  // Through the boundary layer's split as the evaluator formed it, on the tape.
  if (use_energy_balance_) {
    gc = boundary_layer_split_CO2<S>(gc).first;
  }

  // Follow the forward model's branch: with no conductance and no root of A(ci)
  // it reported ci = +inf and the saturated assimilation, nothing to place. A'
  // is a tangent through the kernel the residual calls, not a hand-written twin.
  const bool no_coordinate = !std::isfinite(ci_star);
  double dT2_dci = 0.0;
  if (!no_coordinate) {
    const PhotoCapacity<double> at_cap{
        to_passive(cap.vcmax), to_passive(cap.transport),
        to_passive(cap.curvature), to_passive(cap.respiration)};
    dT2_dci = assim_slope_at<double>(ci_star, at_cap) * umol_to_mol +
              to_passive(gc) * inv_atm;
  }
  // The capacity and gc are the active things the residual reads, so they are
  // the list; ca_ and inv_atm are double.
  ci_h =
      no_coordinate
          ? S(ci_star)
          : odelia::implicit_value<S>(
                ci_star, dT2_dci,
                [&](const S& c) -> S {
                  const S A = leaf.template assim_colimited_kernel<S>(c, cap);
                  return A * umol_to_mol - gc * (S(leaf.ca_) - c) * S(inv_atm);
                },
                cap, gc);
}

template <class S>
inline Leaf::CollarCoords<S> Leaf::collar_channel_at(
    double sigma_star, const S& collar, const SupplyDraw<S>& draw,
    const leaf_pars<S>& pars, const PhotoCapacity<S>& cap, const S& sigma_h,
    const S& gc, const S& ci_h, bool collar_moves) const {
  using odelia::util::to_passive;
  const S held(to_passive(collar));
  const double gc_per_flux =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  // A ci with no root was handed over as +inf, as collar_coords_at's no-root
  // branch places it; the channel has nothing to carry there.
  const bool no_coordinate = !std::isfinite(to_passive(ci_h));

  // ⚠️ A held collar's slopes come back NaN on purpose: spelled 0.0, a missing
  // channel leaves every number finite (hazard 17), and supply_derivatives
  // refuses NaN by name. Only dprofit_dcollar_at reads one, and it asks for it.
  if (!collar_moves) {
    const S absent(std::numeric_limits<double>::quiet_NaN());
    return CollarCoords<S>{odelia::with_slope<S>{sigma_h, absent},
                           odelia::with_slope<S>{ci_h, absent}};
  }

  // The collar's own channel. The conductivity's value is the table's, since the
  // solve drove dprofit_at_collar_psi's table assembly to zero (a closed-form
  // value makes M non-zero there); its rows are the curve's (closed_form_rows.hpp).
  const S& kmax = pars[par_kmax];
  const S f_p = curve_with_closed_form_rows<S>(stem_curve_integral_deriv(to_passive(collar)),
                               held, pars[par_stem_P50], pars[par_stem_c]);
  const S f_sigma = curve_with_closed_form_rows<S>(stem_curve_integral_deriv(sigma_star),
                                   sigma_h, pars[par_stem_P50],
                                   pars[par_stem_c]);

  // dsigma/dcollar: the inverse's own slope times the transport's, as
  // dprofit_at_collar_psi forms it, read where psi_stem_carrying placed the stem
  // for the soil draw. Formed only where read, since the inverse can refuse.
  const auto dsigma_du = [&]() -> S {
    const double sigma_at = psi_stem_carrying(
        to_passive(draw.E_up.value), to_passive(collar),
        "Leaf::collar_channel_at, placing the stem for the soil draw");
    const S recip = S(1.0) / f_sigma;
    return S(1.0 / stem_curve_integral_deriv(sigma_at)) +
           (recip - S(to_passive(recip)));
  };
  // psi_crit does not move with the collar, nor does a stem held at a dry pin:
  // its channel is an exact zero, which formed would be e^40 times terms of
  // order e^-40, platform-dependent rounding. A stem seated at a shut leaf's
  // collar moves with it one for one.
  const S V = stem_held_at_dry_pin(sigma_star)
                  ? S(0.0)
                  : stem_held_at_collar(sigma_star, to_passive(collar))
                        ? S(1.0)
                        : dsigma_du() * (draw.E_up.slope / kmax + f_p);

  const S gc_const = S(gc_per_flux);
  const S inv = S(inv_atm);
  // The conductance partials are kmax * f at each end. The residual flux moves
  // with the collar through its ceiling, appended as dprofit_at_collar_psi
  // appends it; an exact zero at g_min = 0.
  S dgc_dsigma = gc_const * kmax * f_sigma;
  S dgc_dp = -gc_const * kmax * f_p;
  dgc_dp = dgc_dp - gc_const * S(residual_flux_dpsi_upstream(to_passive(collar),
                                                             Tleaf_));

  // The same condition governs the slope: with no conductance at all ci is
  // placed by A(ci) = 0, which the collar does not enter.
  S dci_dp = S(0.0);
  if (!no_coordinate && stom_cond_CO2_ != 0.0) {
    // The partials above are before the boundary layer's split; gc is after it.
    if (use_energy_balance_) {
      const S bl_slope = boundary_layer_split_slope_at<S>(gc);
      dgc_dsigma = dgc_dsigma * bl_slope;
      dgc_dp = dgc_dp * bl_slope;
    }
    const S A_prime = assim_slope_at<S>(ci_h, cap);
    const S g_ci = A_prime * umol_to_mol + gc * inv;
    const S ca_minus_ci = S(ca_) - ci_h;
    const S dci_dsigma = (dgc_dsigma * ca_minus_ci * inv) / g_ci;
    const S dci_dp_expl = (dgc_dp * ca_minus_ci * inv) / g_ci;
    dci_dp = dci_dsigma * V + dci_dp_expl;
  }

  // The step is exactly zero in VALUE, so both coordinates are the residuals'
  // bit for bit and only the collar's channel is added.
  const S step = collar - held;
  return CollarCoords<S>{odelia::with_slope<S>{sigma_h + V * step, V},
                         odelia::with_slope<S>{ci_h + dci_dp * step, dci_dp}};
}

template <class S>
inline S Leaf::psi_crit_at(const leaf_pars<S>& pars) const {
  const S& c = pars[par_stem_c];
  const S b = phylloptim::weibull_b_from_P50<S>(pars[par_stem_P50], c);
  return b * pow(S(std::log(1.0 / k_crit_fraction)), 1.0 / c);
}

template <class S>
inline S Leaf::stem_curve_end_at(const leaf_pars<S>& pars) const {
  const S& c = pars[par_stem_c];
  const S b = phylloptim::weibull_b_from_P50<S>(pars[par_stem_P50], c);
  return b * pow(S(vulnerability_x_end), 1.0 / c);
}

template <Leaf::CostCurve K, class S>
inline S Leaf::dry_end_at(const leaf_pars<S>& pars) const {
  if constexpr (cost_domain_ends_at_psi_crit<K>()) {
    return psi_crit_at<S>(pars);
  } else {
    return stem_curve_end_at<S>(pars);
  }
}

template <Leaf::CostCurve K, class S>
inline S Leaf::profit_at(const S& collar, const SupplyDraw<S>& draw,
                         const leaf_pars<S>& pars, bool collar_moves) const {
  static_assert(K == CostCurve::TF24 || K == CostCurve::TF24_floor,
                "the reverse-mode surface covers TF24 and TF24_floor; another "
                "curve is a row in profit_at and dprofit_dcollar_at");
  // One capacity for the point and for the assimilation read off it, so the ci
  // the residual placed and the A evaluated there are the same function.
  const PhotoCapacity<S> cap = photo_capacity_at<S>(pars);
  const CollarCoords<S> at = collar_coords_at<S>(opt_psi_stem_, ci_, collar,
                                                 draw, pars, cap, collar_moves);
  return profit_from_coords<K, S>(collar, pars, cap, at.sigma.value,
                                  at.ci.value);
}

template <Leaf::CostCurve K, class S>
inline S Leaf::profit_from_coords(const S& collar, const leaf_pars<S>& pars,
                                  const PhotoCapacity<S>& cap, const S& sigma,
                                  const S& ci) const {
  static_assert(K == CostCurve::TF24 || K == CostCurve::TF24_floor,
                "the reverse-mode surface covers TF24 and TF24_floor; another "
                "curve is a row in profit_at and dprofit_dcollar_at");
  // With no conductance and no root, ci = +inf and the saturated assimilation,
  // as the forward reported. With no conductance and a root, A is zero exactly,
  // carrying the kernel's rows, which are what move it off zero.
  S A;
  if (!std::isfinite(ci_)) {
    A = assim_saturated_kernel<S, S>(cap);
  } else {
    A = assim_colimited_kernel<S>(ci, cap);
    if (stom_cond_CO2_ == 0.0) {
      A = A - S(odelia::util::to_passive(A));
    }
  }
  S cost = hydraulic_cost_TF_kernel<S>(sigma, pars);
  if constexpr (K == CostCurve::TF24_floor) {
    // TF24's own cost plus lambda_o*E, as the parent's expression plus a term,
    // so the reduction at lambda_o = 0 is exact, fused multiply-add included.
    cost = cost + pars[par_TF24_floor_lambda_o] *
                      transpiration_at<S>(sigma, collar, pars);
  }
  return A - cost;
}

template <Leaf::CostCurve K, class S>
inline Leaf::BilledOutputs<S> Leaf::billed_outputs_at(const S& collar,
                                                      const SupplyDraw<S>& draw,
                                                      const leaf_pars<S>& pars,
                                                      const HeldCoords<S>* placed) const {
  using odelia::util::to_passive;
  refuse_thermal_damage("Leaf::billed_outputs_at");
  check_surface_curve<K>("Leaf::billed_outputs_at");
  check_draw(to_passive(collar), draw);
  // Placed coordinates hold only at an interior point and at the collar
  // operating_collar placed them at; either check failing means the leaf moved
  // between the calls, and they would answer for another point.
  const bool reuse = placed != nullptr && placed->placed();
  if (reuse) {
    if (operating_point_kind_ != OperatingPointKind::Interior) {
      util::stop(std::string("Leaf::billed_outputs_at: the held coordinates were "
                             "placed at an interior point and are being read "
                             "at one that is ") +
                 operating_point_kind_name(operating_point_kind_));
    }
    if (placed->at != to_passive(collar)) {
      util::stop("Leaf::billed_outputs_at: the held coordinates were placed at a "
                 "collar of " + util::to_string(placed->at) +
                 " and are being read at " + util::to_string(to_passive(collar)));
    }
  }
  BilledOutputs<S> out;
  const std::size_t n = static_cast<std::size_t>(supply_n_layers());

  // The draws carry supplied rows, not a second supply: the draw holds each
  // layer's uptake at the passive collar with the state's rows, and the collar's
  // movement is one supplied slope per layer. The step is zero in value, so each
  // uptake is the draw's bit for bit and the supply reaches the tape once.
  const S step = collar - S(draw.at);
  odelia::util::check_length(draw.duptake_dp.size(), draw.uptake.size());
  // Built, not copied and overwritten: a copy records a statement per layer.
  out.uptake.reserve(draw.uptake.size());
  for (std::size_t j = 0; j < draw.uptake.size(); ++j) {
    out.uptake.push_back(draw.uptake[j] + S(draw.duptake_dp[j]) * step);
  }
  (void)n;

  // Profit reads the HELD collar at an interior point: the envelope theorem as
  // an omission. The marginal is zero there by the solve, so carrying the collar
  // would add M dp/dtheta, of the size of the solve's tolerance. A shut leaf is
  // pinned: the wet bound places its collar and profit carries its movement.
  const bool interior = operating_point_kind_ == OperatingPointKind::Interior;
  const S collar_held(to_passive(collar));
  if (reuse) {
    // The same coordinates profit_at would solve for at the held collar: the
    // residuals see a held collar there too, so these are its numbers and rows.
    out.profit = profit_from_coords<K, S>(collar_held, pars, placed->cap,
                                          placed->sigma, placed->ci);
  } else {
    out.profit = profit_at<K, S>(interior ? collar_held : collar, draw, pars,
                                 !interior);
  }
  // The residual flux at the LIVE collar: the envelope theorem is profit's,
  // and the cuticle's loss moves with the collar through the stem's ceiling.
  // An exact zero carrying the g_min row at g_min = 0.
  out.residual_flux = residual_flux_at<S>(collar, pars, draw);
  return out;
}

template <Leaf::CostCurve K, class S>
inline S Leaf::dprofit_dcollar_at(const S& collar, const SupplyDraw<S>& draw,
                                  const leaf_pars<S>& pars) const {
  static_assert(K == CostCurve::TF24 || K == CostCurve::TF24_floor,
                "the reverse-mode surface covers TF24 and TF24_floor; another "
                "curve is a row in profit_at and dprofit_dcollar_at");
  // Every factor is a value, so this is arithmetic at S and dM/dp a first
  // derivative of it: a table's value with the true function's derivative
  // disagree by the interpolation error, enough to flip a curvature's sign.
  const PhotoCapacity<S> cap = photo_capacity_at<S>(pars);
  const CollarCoords<S> at =
      collar_coords_at<S>(opt_psi_stem_, ci_, collar, draw, pars, cap, true);
  return dprofit_dcollar_from_coords<K, S>(collar, pars, cap, at);
}

template <Leaf::CostCurve K, class S>
inline S Leaf::dprofit_dcollar_at(const S& collar, const SupplyDraw<S>& draw,
                                  const leaf_pars<S>& pars,
                                  const HeldCoords<S>& held) const {
  using odelia::util::to_passive;
  check_draw(to_passive(collar), draw);
  if (held.at != to_passive(collar)) {
    util::stop("Leaf::dprofit_dcollar_at: the held coordinates were placed at a "
               "collar of " + util::to_string(held.at) +
               " and are being read at " + util::to_string(to_passive(collar)));
  }
  const CollarCoords<S> at =
      collar_channel_at<S>(opt_psi_stem_, collar, draw, pars, held.cap,
                           held.sigma, held.gc, held.ci, true);
  return dprofit_dcollar_from_coords<K, S>(collar, pars, held.cap, at);
}

template <Leaf::CostCurve K, class S>
inline S Leaf::dprofit_dcollar_from_coords(const S& collar, const leaf_pars<S>& pars,
                                           const PhotoCapacity<S>& cap,
                                           const CollarCoords<S>& at) const {
  static_assert(K == CostCurve::TF24 || K == CostCurve::TF24_floor,
                "the reverse-mode surface covers TF24 and TF24_floor; another "
                "curve is a row in profit_at and dprofit_dcollar_at");
  const S A_prime = assim_slope_at<S>(at.ci.value, cap);
  // dC/dpsi_stem as a tangent through the cost kernel, not a slope beside it.
  S C_prime = cost_slope_at<S>(at.sigma.value, pars);
  if constexpr (K == CostCurve::TF24_floor) {
    // d/dsigma of lambda_o kmax (G(sigma) - G(collar)) is lambda_o kmax f, with
    // the table's conductivity because the solve ran on the table.
    C_prime = C_prime +
              pars[par_TF24_floor_lambda_o] * pars[par_kmax] *
                  curve_with_closed_form_rows<S>(stem_curve_integral_deriv(opt_psi_stem_),
                                 at.sigma.value, pars[par_stem_P50],
                                 pars[par_stem_c]);
  }

  S M = A_prime * at.ci.slope - C_prime * at.sigma.slope;
  if constexpr (cost_reads_upstream<K>()) {
    using odelia::util::to_passive;
    // -dC/dcollar at fixed sigma: the floor's lambda_o*E reads the collar as its
    // upstream end, as cost_deriv_upstream bills it in the solve. At the LIVE
    // collar, so the term's own dM/dcollar reaches the implicit row.
    M = M + pars[par_TF24_floor_lambda_o] * pars[par_kmax] *
                curve_with_closed_form_rows<S>(stem_curve_integral_deriv(to_passive(collar)),
                                               collar, pars[par_stem_P50], pars[par_stem_c]);
  }
  return M;
}


template <class S>
inline S Leaf::wet_bound_at(double bound_x, const SupplyDraw<S>& draw,
                            const leaf_pars<S>& pars) const {
  using odelia::util::to_passive;
  check_draw(bound_x, draw);
  const Leaf& leaf = *this;

  // The bound's slope is the supply's own conductance, carried as one number:
  // where a kink makes the analytic branch invalid, supply_draw_at has already
  // differenced it.
  const double dflux_dx = to_passive(draw.E_up.slope);

  // Nothing conducts: the bound is the gravity balance the draw carries,
  // with its soil row (find_shut_wet_bound placed the solve there).
  if (draw.no_conduction) {
    return draw.gravity_balance;
  }
  // The residual is the draw less the residual flux; the draw was taken at the
  // passive bound, so the supply is not recorded twice. At g_min = 0 the draw
  // sits a step inside the root (hazard 18), which puts slope * step in every
  // row: taking it back off is the Taylor correction. At g_min > 0 it is zero.
  const double step = std::isfinite(wet_bound_root_)
                          ? bound_x - wet_bound_root_
                          : 0.0;
  const double dres_dx =
      dflux_dx - residual_flux_dpsi_upstream(bound_x, Tleaf_);
  const S held(bound_x);
  // ⚠️ Every active input the residual reads goes in the list, the roots'
  // ceiling among them when in force, or the reverse sweep drops its rows while
  // a forward tangent keeps them (hazard 18; test_transpose).
  if (step != 0.0) {
    const S at_root = draw.E_up.value - step * draw.E_up.slope;
    auto residual = [&](const S&) -> S {
      return at_root - leaf.template residual_flux_at<S>(held, pars, draw);
    };
    return draw.has_uptake_limit
               ? odelia::implicit_value<S>(bound_x, dres_dx, residual, at_root,
                                           draw.uptake_limit, pars)
               : odelia::implicit_value<S>(bound_x, dres_dx, residual, at_root,
                                           pars);
  }
  auto residual = [&](const S&) -> S {
    return draw.E_up.value -
           leaf.template residual_flux_at<S>(held, pars, draw);
  };
  const S wet = draw.has_uptake_limit
                    ? odelia::implicit_value<S>(bound_x, dres_dx, residual,
                                                draw.E_up.value,
                                                draw.uptake_limit, pars)
                    : odelia::implicit_value<S>(bound_x, dres_dx, residual,
                                                draw.E_up.value, pars);
  return wet;
}

template <Leaf::CostCurve K, class S>
inline S Leaf::dry_end_collar_at(double bound_x, const SupplyDraw<S>& draw,
                                 const leaf_pars<S>& pars) const {
  using odelia::util::to_passive;
  check_draw(bound_x, draw);
  const Leaf& leaf = *this;
  // The supply's own conductance at the bound, as wet_bound_at reads it.
  const double dflux_dx = to_passive(draw.E_up.slope);

  // T1 with the stem held at dry_end_at: the collar at which the soil delivers
  // what the column can still carry, so only the collar's integral and the flux
  // move. At the end of the stem curve the stem's rows are zero to e^-40.
  const double dT_dx =
      -to_passive(pars[par_kmax]) * stem_curve_integral_deriv(bound_x) -
      dflux_dx;
  const S dry = odelia::implicit_value<S>(
      bound_x, dT_dx,
      [&](const S& x) -> S {
        return leaf.template transpiration_at<S>(dry_end_at<K, S>(pars), x,
                                                 pars) -
               draw.E_up.value;
      },
      pars, draw.E_up.value);
  return dry;
}

// An interior point is the only kind whose condition is a second derivative, so
// it is the only one that needs a curvature handed to it; every other kind's
// condition is first order and composes here.
template <Leaf::CostCurve K, class S>
inline S Leaf::operating_collar(const SupplyDraw<S>& draw, const leaf_pars<S>& pars,
                                double interior_curvature, HeldCoords<S>* placed) {
  // Refused at every kind, not only where the kernels are read: a surface that
  // answers at some kinds and not others is one a caller cannot plan around.
  refuse_thermal_damage("Leaf::operating_collar");
  check_surface_curve<K>("Leaf::operating_collar");
  check_draw(opt_root_psi_, draw);
  // Unplaced until the interior arm places it, so coordinates left from an
  // earlier point can never be read at this one.
  if (placed != nullptr) {
    placed->at = std::numeric_limits<double>::quiet_NaN();
  }
  S collar = S(opt_root_psi_);
  switch (operating_point_kind_) {
    case OperatingPointKind::Interior: {
      // Closed by the residual like every other kind: dM/dp is what the theorem
      // divides by, and every parameter row is taped from M itself.
      // ⚠️ Not at double (hazard 5), where this IS opt_root_psi_ but the arguments
      // are evaluated first: the curvature costs two model evaluations a solve.
      if constexpr (!std::is_same_v<S, double>) {
        const double curvature =
            std::isnan(interior_curvature)
                ? this->template profit_curvature_at_collar<K>()
                : interior_curvature;
        // implicit_value refuses only a zero or non-finite denominator, so these
        // are this arm's: -1e-300 would return rows of ~1e296 unflagged. A
        // non-negative curvature indicts itself; a tiny negative one is flat.
        if (!(curvature < 0.0)) {
          util::stop("Leaf::operating_collar: the profit curvature at this interior "
                     "point is " + util::format_double(curvature) +
                     ", which is not negative, so the implicit function "
                     "theorem has nothing to divide by");
        }
        if (curvature > -interior_curvature_floor) {
          util::stop("Leaf::operating_collar: the profit curvature at this interior "
                     "point is " + util::format_double(curvature) +
                     ", inside the floor of " +
                     util::format_double(interior_curvature_floor) +
                     ", so the collar rows would divide by too little to answer");
        }
        // The coordinates are placed before the residual, which is rewound, so
        // they stay on the caller's tape for billed_outputs_at. The residual
        // reads them, so they join the list: one left out is a column of zeros.
        if (placed != nullptr) {
          const S at = S(opt_root_psi_);
          placed->cap = photo_capacity_at<S>(pars);
          held_coords_at<S>(opt_psi_stem_, ci_, at, draw, pars, placed->cap,
                            placed->sigma, placed->gc, placed->ci);
          placed->at = opt_root_psi_;
          const HeldCoords<S>& h = *placed;
          collar = odelia::implicit_value<S>(
              opt_root_psi_, curvature,
              [&](const S& y) -> S {
                return dprofit_dcollar_at<K, S>(y, draw, pars, h);
              },
              draw, pars, h);
        } else {
          // The whole draw, because dprofit_dcollar_at reads more of it than
          // E_up: SupplyDraw names its own active members for this.
          collar = odelia::implicit_value<S>(
              opt_root_psi_, curvature,
              [&](const S& y) -> S {
                return dprofit_dcollar_at<K, S>(y, draw, pars);
              },
              draw, pars);
        }
      }
      break;
    }
    // A shut leaf sits ON the wet bound: the collar at which uptake meets the
    // residual flux, so the bound is the point and it moves with the soil.
    case OperatingPointKind::BoundarySoil:
    case OperatingPointKind::HydraulicShutdown:
      collar = wet_bound_at<S>(opt_root_psi_, draw, pars);
      break;
    case OperatingPointKind::BoundaryDry:
      // The continuity root: T1 with the stem held at its dry end.
      collar = dry_end_collar_at<K, S>(opt_root_psi_, draw, pars);
      break;
    default:
      util::stop(std::string("Leaf::operating_collar: no operating point to place at "
                             "one that is ") +
                 operating_point_kind_name(operating_point_kind_));
  }
  return collar;
}

template <Leaf::CostCurve K>
inline double Leaf::profit_curvature_at_collar() {
  const double p0 = opt_root_psi_;
  // An absolute step: the collar runs 0.5 to 6 MPa, and steps 1e-04 to 1e-06
  // agree to seven figures (1e-02 is truncation, 1e-07 round-off). Within
  // dry_end_neighbourhood of the dry end a step across the distance reads the
  // wrong decade's slope; the marginal is near linear in the log of the
  // distance, so it is differenced there, extrapolated against the half step,
  // with both points inside the end (hazard 19).
  const OperatingPointKind kind0 = operating_point_kind_;
  const double wet0 = wet_bound_root_;
  const double to_dry = std::isfinite(dry_bound_root_)
                            ? dry_bound_root_ - p0
                            : std::numeric_limits<double>::infinity();
  if (to_dry > 0.0 && to_dry < dry_end_neighbourhood) {
    auto slope_in_log_distance = [&](double hs, bool* ok) -> double {
      bool in_ok = false, out_ok = false;
      const double in = dprofit_at_collar_psi<K>(dry_bound_root_ - to_dry * std::exp(-hs), &in_ok);
      const double out = dprofit_at_collar_psi<K>(dry_bound_root_ - to_dry * std::exp(hs), &out_ok);
      *ok = in_ok && out_ok && std::isfinite(in) && std::isfinite(out);
      return (in - out) / (2.0 * hs);
    };
    bool ok1 = false, ok2 = false;
    const double coarse = slope_in_log_distance(0.2, &ok1);
    const double fine = slope_in_log_distance(0.1, &ok2);
    // Back to the point the solve placed, with the bounds it polished: the
    // restore prepares the solve again, which clears them.
    const double dry0 = dry_bound_root_;
    evaluate_root_collar_psi_for<K>(p0);
    operating_point_kind_ = kind0;
    dry_bound_root_ = dry0;
    wet_bound_root_ = wet0;
    if (ok1 && ok2) {
      // The quotient above is -dM/ds, its inner point lying at the smaller
      // distance, and dM/dpsi = dM/ds * ds/dpsi with s = ln(end - psi), so
      // ds/dpsi = -1/(end - psi): the two signs cancel.
      return ((4.0 * fine - coarse) / 3.0) / to_dry;
    }
  }
  const double h = 1e-5;
  bool up_ok = false, down_ok = false, here_ok = false;
  const double up = dprofit_at_collar_psi<K>(p0 + h, &up_ok);
  const double down = dprofit_at_collar_psi<K>(p0 - h, &down_ok);
  // One side is enough: an interior collar within h of a reachable bound has a
  // neighbour where dprofit refuses, and a NaN makes plant refuse the whole
  // sweep. A one-sided difference is the same derivative, one order lower.
  const double here =
      (up_ok && down_ok) ? 0.0 : dprofit_at_collar_psi<K>(p0, &here_ok);
  // Back to the point the solve placed. evaluate_root_collar_psi restores every
  // number but tags the point `prescribed`, so the kind is put back by hand:
  // operating_collar switches on it and would throw.
  const double dry0 = dry_bound_root_;
  evaluate_root_collar_psi_for<K>(p0);
  operating_point_kind_ = kind0;
  dry_bound_root_ = dry0;
  wet_bound_root_ = wet0;
  const bool have_up = up_ok && std::isfinite(up);
  const bool have_down = down_ok && std::isfinite(down);
  const bool have_here = here_ok && std::isfinite(here);
  if (have_up && have_down) {
    return (up - down) / (2.0 * h);
  }
  if (have_here && have_up) {
    return (up - here) / h;
  }
  if (have_here && have_down) {
    return (here - down) / h;
  }
  // The same convention duptake_dpsi uses at a kink: a caller that cannot get a
  // curvature here has to know that, not be handed a zero it will divide by.
  return std::numeric_limits<double>::quiet_NaN();
}


inline double Leaf::hydraulic_cost_TF(double psi_stem) {

  hydraulic_cost_ = hydraulic_cost_TF_kernel(psi_stem);

return hydraulic_cost_;
}

// Sicangco et al. (2026) Eqns 8-9: the fraction of maximum PSII damage at this
// leaf temperature, in [0,1]; `r` follows gsthermal, not the printed Eqn 8 (see
// use_thermal_cost_). Exactly 0.0 with the gate off, so callers add it freely.
inline double Leaf::thermal_cost_at(double leaf_temp) const {
  if (!use_thermal_cost_) {
    return 0.0;
  }
  const double span = T50_ - Tcrit_;
  if (!(span > 0.0)) {
    util::stop("thermal cost needs T50 > Tcrit; got T50 = " +
               util::format_double(T50_) + ", Tcrit = " +
               util::format_double(Tcrit_));
  }
  return 1.0 / (1.0 + std::exp(-(2.0 / span) * (leaf_temp - T50_)));
}

inline double Leaf::thermal_cost_temp_deriv(double leaf_temp) const {
  if (!use_thermal_cost_) {
    return 0.0;
  }
  const double span = T50_ - Tcrit_;
  const double TC = thermal_cost_at(leaf_temp);
  return (2.0 / span) * TC * (1.0 - TC);
}

// Profit functions

inline double Leaf::profit_psi_stem_TF(double psi_stem, double psi_upstream) {
set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF(psi_stem);

  return benefit_ - cost;
}


// `CF77_lambda_ * E`, written into `hydraulic_cost_` because it IS a carbon flux
// (umol C per kg water times kg m^-2 s^-1). Straight after
// set_leaf_states_rates_from_psi_stem the memo returns that point's E.
inline double Leaf::hydraulic_cost_CF77(double psi_stem,
                                                double psi_upstream) {
  hydraulic_cost_ = cf77_price() * operating_point_flux(psi_stem, psi_upstream);
  return hydraulic_cost_;
}


// The CF77 price with the optional Medlyn soil-moisture factor, rising without
// bound toward the wilting point:
//
//   beta = (theta - theta_w) / (theta_fc - theta_w),   price = lambda / beta
//
// beta is clamped at 1, so CF77_lambda_ is a floor and theta >= theta_fc is plain
// CF77 exactly. The off path returns the member, not a division by 1.0, so it
// stays bit-identical to plain CF77.
inline double Leaf::cf77_price() const {
  if (!CF77_soil_beta_) { return CF77_lambda_; }
  const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
  return CF77_lambda_ / (beta > 1.0 ? 1.0 : beta);
}


inline double Leaf::profit_psi_stem_CF77(double psi_stem,
                                                 double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_CF77(psi_stem, psi_upstream);

  return benefit_ - cost;
}


// gamma * (dpsi)^2, in carbon units like TF24 and CF77, so it goes in
// `hydraulic_cost_`. `dpsi*dpsi`, not `pow`: the optimiser's inner loop (hazard 5).
inline double Leaf::hydraulic_cost_JS22(double psi_stem,
                                        double psi_upstream) {
  const double dpsi = psi_stem - psi_upstream;
  hydraulic_cost_ = JS22_gamma * dpsi * dpsi;
  return hydraulic_cost_;
}


inline double Leaf::profit_psi_stem_JS22(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_JS22(psi_stem, psi_upstream);

  return benefit_ - cost;
}


// The integral of `CMax_a*psi + CMax_b` from psi_upstream, zero at no flow:
//
//     C = (psi - psi_up) * [ CMax_a*(psi + psi_up)/2 + CMax_b ]
//
// Factored: the difference of squares cancels as the potentials converge (off by
// 1.2e-04 at a 1e-12 drop), which a gradient's step reaches near no flow. For the
// same reason JS22, this family's (2 gamma, -2 gamma psi_up) member, keeps its own.
inline double Leaf::hydraulic_cost_CMax(double psi_stem,
                                        double psi_upstream) {
  hydraulic_cost_ = (psi_stem - psi_upstream) *
                    (0.5 * CMax_a * (psi_stem + psi_upstream) + CMax_b);
  return hydraulic_cost_;
}


inline double Leaf::profit_psi_stem_CMax(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_CMax(psi_stem, psi_upstream);

  return benefit_ - cost;
}


// TF24's cost plus CF77's price, both carbon fluxes. Each term is its parent's own
// expression in its operand order, so lambda_o = 0 and TF24_cost_scale = 0 reduce
// bit-exactly, under FMA contraction too (hazards.md, "One solver: the long
// form"). Reassociating breaks that; test_tf24_floor_reduces guards it.
inline double Leaf::hydraulic_cost_TF24_floor(double psi_stem,
                                           double psi_upstream) {
  hydraulic_cost_ = hydraulic_cost_TF_kernel(psi_stem) +
                    TF24_floor_lambda_o * operating_point_flux(psi_stem, psi_upstream);
  return hydraulic_cost_;
}


inline double Leaf::boundary_layer_cond_H2O() const {
  if (!use_energy_balance_) {
    return std::numeric_limits<double>::infinity();
  }
  // Molar density of air at air temperature, mol m^-3, as plantecophys's CMOLAR.
  const double air_molar_density =
      atm_kpa_ * kPa_to_Pa / (gas_constant * (Tair_ + C_to_K));
  return bl_vapour_heat_ratio * air_molar_density / ra_;
}

inline double Leaf::boundary_layer_co2_excess() const {
  const double r_s = H2O_CO2_stom_diff_ratio_;
  const double r_b = std::cbrt(r_s * r_s);
  return (r_s - r_b) / boundary_layer_cond_H2O();
}

template <class T>
inline std::pair<T, T> Leaf::boundary_layer_split_CO2(const T& gc_hat) const {
  using std::exp;
  using std::log1p;
  using odelia::util::to_passive;
  const double k = boundary_layer_co2_excess();
  const double s = 1.0 / bl_co2_split_floor_width;
  const T d = T(1.0) - T(k) * gc_hat;
  // The soft maximum of d and the floor, d + softplus(d0 - d) at width w, and
  // its slope in d, the logistic of (d - d0)/w. Written so neither overflows:
  // e = exp(-|z|) is at most 1 whichever side of the floor d is.
  const T z = T(s) * (T(bl_co2_split_floor) - d);
  const bool above = to_passive(z) <= 0.0;
  const T e = above ? T(exp(z)) : T(exp(-z));
  const T d_s = above ? T(d + log1p(e) / T(s))
                      : T(T(bl_co2_split_floor) + log1p(e) / T(s));
  const T dds_dd = above ? T(T(1.0) / (T(1.0) + e)) : T(e / (T(1.0) + e));
  const T gc = gc_hat / d_s;
  const T slope = T(1.0) / d_s + gc_hat * T(k) * dds_dd / (d_s * d_s);
  return {gc, slope};
}

template <class T>
inline T Leaf::boundary_layer_split_slope_at(const T& gc) const {
  using odelia::util::to_passive;
  const double k = boundary_layer_co2_excess();
  const T phi = T(1.0) + T(k) * gc;
  // d = 1/phi above the floor; within 30 widths of it the floor has moved the
  // value and this form would not be its slope.
  if (!(1.0 / to_passive(phi) >
        bl_co2_split_floor + 30.0 * bl_co2_split_floor_width)) {
    util::stop("the boundary layer's split of the CO2 path is on its soft floor "
               "at this operating point (g_c = " + util::to_string(to_passive(gc)) +
               "): the transpiration is more than the boundary layer can carry, "
               "so the point is not a physical one to differentiate.");
  }
  return phi * phi;
}

inline double Leaf::stom_cond_stomatal_H2O() const {
  const double g_total = stom_cond_H2O();
  if (!use_energy_balance_) {
    return g_total;
  }
  // A shut leaf has no conductance on either side of the series, and a
  // non-finite one has nothing to subtract from.
  if (!(g_total > 0.0) || !std::isfinite(g_total)) {
    return g_total;
  }
  const double g_b = boundary_layer_cond_H2O();
  if (!(g_total < g_b)) {
    util::stop("no stomatal conductance exists at this operating point: the total "
               "leaf-to-air conductance " + util::to_string(g_total) +
               " mol m^-2 s^-1 is not below the boundary layer's own " +
               util::to_string(g_b) + ", so the transpiration is more than the "
               "boundary layer can carry at this deficit (ra_ = " +
               util::to_string(ra_) + " s m^-1). Check wind_speed_ and d_.");
  }
  return 1.0 / (1.0 / g_total - 1.0 / g_b);
}

// See the declaration for why zero on the other curves is a statement about the
// curve and not about the cost.
inline double Leaf::shadow_cost() const {
  return cost_curve_ == CostCurve::TF24_floor
           ? TF24_floor_lambda_o * transpiration_
           : 0.0;
}


inline double Leaf::profit_psi_stem_TF24_floor(double psi_stem,
                                            double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF24_floor(psi_stem, psi_upstream);

  return benefit_ - cost;
}


// --- the product-objective family --------------------------------------------
//
// Eller's reduction factor: the conductivity fraction rescaled to 1 at zero
// tension and exactly 0 at `psi_crit`, through the k_crit_fraction psi_crit and
// Sperry's k_crit also read, so the three cannot disagree.
// Clamped at zero because the no-flow branch evaluates past psi_crit, where g < 0
// and A = -R_d made a shut leaf's product positive. `g < 0`, not `std::max`, so a
// NaN propagates.
inline double Leaf::sox_reduction(double psi_stem) const {
  const double g = (proportion_of_conductivity(psi_stem) - k_crit_fraction) /
                   (1.0 - k_crit_fraction);
  return g < 0.0 ? 0.0 : g;
}


// Its slope, negative: `|f'| = f*(c/b)*(psi/b)^(c-1)`, lambda_TF24's analytic form.
inline double Leaf::sox_reduction_deriv(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  const double abs_fprime =
      f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
  return -abs_fprime / (1.0 - k_crit_fraction);
}


// A product, not in carbon units: the argmax is log A + log g's, but the value
// must not be compared with a TF24 profit.
// ⚠️ hydraulic_cost_ is poisoned (hazard 8): with no subtracted cost, leaving it
// would report whichever curve ran last on a reused Leaf.
inline double Leaf::profit_psi_stem_SOX(double psi_stem,
                                        double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * sox_reduction(psi_stem);
}


// Jones' linear reduction factor on the derived `psi_crit`, clamped at zero, NaN
// propagating, as `sox_reduction` is.
inline double Leaf::jw26_reduction(double psi_stem) const {
  const double g = 1.0 - psi_stem / psi_crit;
  return g < 0.0 ? 0.0 : g;
}


// Constant, the whole difference from SOX, whose slope follows the vulnerability
// curve. Zero once the clamp binds.
inline double Leaf::jw26_reduction_deriv(double psi_stem) const {
  return psi_stem > psi_crit ? 0.0 : -1.0 / psi_crit;
}


inline double Leaf::profit_psi_stem_JW26(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * jw26_reduction(psi_stem);
}


// --- least-cost / cost-per-unit-carbon ---------------------------------------

// `W = E_mol + beta*Vcmax_mol`, mol m^-2 s^-1; see the declaration for why a
// dimensionless beta needs both conversions. It reads `vcmax_`, not `vcmax_25`, so
// with the energy balance on it moves with the solved leaf temperature, as the
// ratio Prentice et al. minimise; cost_deriv leaves that motion out and
// dprofit_energy_balance_term adds it.
inline double Leaf::leastcost_flux(double psi_stem, double psi_upstream) {
  return operating_point_flux(psi_stem, psi_upstream) * kg_to_mol_h2o +
         LeastCost_beta * vcmax_ * umol_to_mol;
}


// A ratio, not a carbon flux; `hydraulic_cost_` is poisoned as for SOX (hazard 8).
inline double Leaf::profit_psi_stem_LeastCost(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ / leastcost_flux(psi_stem, psi_upstream);
}


// `lambda = A*kg_to_mol_h2o/W`, reading the solved `assim_colimited_` as
// `lambda_SOX` does: valid only at the potential last evaluated.
inline double Leaf::lambda_LeastCost(double psi_stem, double psi_upstream) {
  if (!use_energy_balance_) {
    return assim_colimited_ * kg_to_mol_h2o /
           leastcost_flux(psi_stem, psi_upstream);
  }
  // W's capacity term moves with E through Tleaf, adding beta*Vcmax_mol'(T)*dT/dE:
  // the slope dprofit_energy_balance_term takes in the first-order condition.
  double dT_dE = 0.0;
  const double Tleaf =
      leaf_temp_from_E(operating_point_flux(psi_stem, psi_upstream), &dT_dE);
  double dW_dE = kg_to_mol_h2o;
  if (dT_dE != 0.0) {
    dW_dE += LeastCost_beta * umol_to_mol * vcmax_temp_deriv(Tleaf) * dT_dE;
  }
  return assim_colimited_ * dW_dE / leastcost_flux(psi_stem, psi_upstream);
}


//optimisation functions

// Everything a single-layer optimiser leaves untouched, cleared rather than
// inherited (hazard 8).
inline void Leaf::clear_collar_solve_state() {
  operating_point_kind_ = OperatingPointKind::Unsolved;
  opt_root_psi_ = util::na_value;
  E_up_ = util::na_value;
  soil_consumption_.assign(soil_consumption_.size(), util::na_value);
}

// A single-layer optimiser maximises over the CLOSED interval from the wet bound
// to dry_end<K>(), endpoints included, through
// util::maximise_over_closed_interval_foc (hazard 11); maximise_profit_over_collar
// tests the gradient's sign at each end instead.

// ===========================================================================
// Sperry et al. (2017) ProfitMax, as Sicangco et al. (2026) implement it
// ---------------------------------------------------------------------------
// Sperry maximises `Profit = CG - HC`, both normalised. Times |A|max that is
// `A - lambda*(k(psi_soil)-k(psi))` with the same argmax, but lambda* moves with
// every driver, so a prescribed lambda is not this model:
//
//     CG      = A(psi)/|A|max
//     HC      = [k(psi_soil)-k(psi)] / [k(psi_soil)-kcrit]
//     lambda* = |A|max / [k(psi_soil) - kcrit]
//
// The profit is in Sperry's units, so a figure overlays the paper's. HC runs 0 to
// 1 however wide the vulnerability curve, where the TF24 cost becomes negligible
// as stem_b widens. kcrit = k_crit_fraction * kmax, so psi_crit is P95, as Sperry.
// HC measures loss from psi_soil on one soil potential, and on a layered profile
// from the collar at which uptake is zero (root_zero_E, prepare_collar_solve's
// wet bracket end).
inline double Leaf::profitmax_reference_psi() {
  if (supply_is_single_layer()) {
    return supply_psi_soil_scalar();
  }
  double bound_a = 0.0, bound_b = 0.0;
  if (!prepare_collar_solve<CostCurve::TF24>(bound_a, bound_b)) {
    return util::na_value;
  }
  return bound_a;
}


inline void Leaf::prepare_profitmax_norms() {
  // A non-finite reference is propagated, not refused: this runs on leaves whose
  // drivers are not seated yet (the batch route seeds |A|max first). The
  // optimiser's `!(profitmax_k_span_ > 0)` guard is the one place that decides it.
  const double psi_soil = profitmax_reference_psi();
  profitmax_k_soil_ =
      leaf_specific_conductance_max_ * proportion_of_conductivity(psi_soil);
  profitmax_k_span_ =
      profitmax_k_soil_ -
      leaf_specific_conductance_max_ * proportion_of_conductivity(psi_crit);
}


inline void Leaf::prepare_profitmax_at(double A_max) {
  prepare_profitmax_norms();
  if (!(std::isfinite(A_max) && A_max > 0.0)) {
    util::stop("prepare_profitmax_at needs a finite positive |A|max -- it is the "
               "normaliser held fixed across a perturbation, so an unset one "
               "means the base point was never solved; got " +
               util::to_string(A_max));
  }
  profitmax_A_max_ = A_max;
}


inline void Leaf::prepare_profitmax() {
  const double psi_soil = profitmax_reference_psi();
  prepare_profitmax_norms();

  // |A|max by the one solver, not a scan, whose argmax is piecewise constant in
  // the parameters (hazard 3). One maximisation of A suffices: |A| peaks at a
  // minimum of A only where A < 0, at closure or along a whole stream that cannot
  // fix carbon (below). dA/dpsi is dJ/dpsi plus the cost derivative, for an
  // identity link. Closure is excluded, as gsthermal drops it (`A[E == 0] <- NA`):
  // the search starts one probe inside the wet bound.
  const auto A_of = [&](double psi) {
    set_leaf_states_rates_from_psi_stem(psi, psi_soil);
    return assim_colimited_;
  };
  const auto dA_of = [&](double psi, bool* ok) {
    const double dJ = dprofit_dpsi_stem<CostCurve::TF24>(psi, psi_soil, ok);
    return dJ + cost_deriv<CostCurve::TF24>(psi, psi_soil);
  };
  const double psi_wet = stem_wet_bound(psi_soil);
  const double width = psi_crit - psi_wet;
  const double lo = psi_wet + 1e-9 * width;   // the open end: closure excluded
  double a_hi = -std::numeric_limits<double>::infinity();
  if (width > 0.0) {
    util::maximise_over_closed_interval_foc(A_of, dA_of, lo, psi_crit, 0,
                                            collar_root_tol,
                                            static_cast<size_t>(ci_niter), &a_hi);
  }
  double a_max = (std::isfinite(a_hi) && a_hi > 0.0) ? a_hi : 0.0;
  // The reachable case where max|A| is not max A: A negative along the whole
  // stream, a leaf too hot to cover respiration. Normalising by zero would leave
  // no objective, so maximise -A, and only then.
  if (!(a_max > 0.0) && width > 0.0) {
    double neg = -std::numeric_limits<double>::infinity();
    util::maximise_over_closed_interval_foc(
        [&](double psi) { return -A_of(psi); },
        [&](double psi, bool* ok) { return -dA_of(psi, ok); }, lo, psi_crit, 0,
        collar_root_tol, static_cast<size_t>(ci_niter), &neg);
    if (std::isfinite(neg) && neg > 0.0) a_max = neg;
  }
  profitmax_A_max_ = a_max;
}

// CG - (HC + TC) at one candidate potential, with the normalisers prepare_profitmax
// seeded. Writes carbon_gain_, hydraulic_cost_norm_ and thermal_cost_ so a caller
// plotting the paper's Figure 2 can read the three components off the object.
inline double Leaf::profit_psi_stem_ProfitMax(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  carbon_gain_ = (profitmax_A_max_ > 0.0)
                     ? assim_colimited_ / profitmax_A_max_
                     : 0.0;
  hydraulic_cost_norm_ =
      (profitmax_k_soil_ -
       leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem)) /
      profitmax_k_span_;
  // At the operating point's leaf temperature; with the energy balance off TC is
  // constant along the stream, as in Sicangco's ProfitMaxTC without it.
  thermal_cost_ = thermal_cost_at(use_energy_balance_
                                      ? leaf_temp_from_E(transpiration_)
                                      : leaf_temp_);

  return carbon_gain_ - (hydraulic_cost_norm_ + thermal_cost_);
}

inline std::vector<double> Leaf::profitmax_curve(int n) {
  if (n < 2) {
    util::stop("profitmax_curve needs at least 2 points");
  }
  prepare_profitmax();
  const double psi_soil = supply_psi_soil_scalar();
  const double psi_wet = stem_wet_bound(psi_soil);
  const double step = (psi_crit - psi_wet) / double(n - 1);
  std::vector<double> out(static_cast<std::size_t>(5 * n));
  for (int i = 0; i < n; ++i) {
    const double p = psi_wet + step * double(i);
    const double profit = profit_psi_stem_ProfitMax(p, psi_soil);
    const std::size_t k = static_cast<std::size_t>(i);
    const std::size_t N = static_cast<std::size_t>(n);
    out[k] = p;
    out[N + k] = carbon_gain_;
    out[2 * N + k] = hydraulic_cost_norm_;
    out[3 * N + k] = thermal_cost_;
    out[4 * N + k] = profit;
  }
  return out;
}

inline void Leaf::set_profitmax_degenerate_state(double psi_stem,
                                                 double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  profit_ = 0;
  carbon_gain_ = 0;
  hydraulic_cost_norm_ = 0;
  thermal_cost_ = 0;
  lambda_emergent_ = util::na_value;
}


// --- runtime curve selection --------------------------------------------------
//
// One switch, so the integer-to-curve mapping is written once. `-Werror=switch`
// makes a curve added to the enum and forgotten here a build failure.
inline std::string Leaf::curve_name(int curve) {
  if (curve < 0 || curve >= n_cost_curves) {
    return "unknown";
  }
  switch (static_cast<CostCurve>(curve)) {
    case CostCurve::TF24: return "TF24";
    case CostCurve::CF77: return "CF77";
    case CostCurve::JS22: return "JS22";
    case CostCurve::CMax: return "CMax";
    case CostCurve::SOX:  return "SOX";
    case CostCurve::JW26: return "JW26";
    case CostCurve::ProfitMax: return "ProfitMax";
    case CostCurve::TF24_floor: return "TF24_floor";
    case CostCurve::LeastCost: return "LeastCost";
  }
  return "unknown";   // unreachable past the bounds check above
}

// The route is checked, not just the curve: these pin the upstream potential at
// psi_soil, so on a collar-route leaf they would answer for a different model
// than `model_route()` reports (hazard 11).
inline void Leaf::require_stem_route(const char* who) const {
  if (route_is_collar_) {
    util::stop(std::string(who) + " is a STEM-route method: it pins the "
               "upstream potential at psi_soil and ignores the soil-to-collar "
               "path. The leaf is seated on the collar route; call "
               "set_model(\"" + curve_name(static_cast<int>(cost_curve_)) +
               "\", \"stem\") first, or use the collar equivalent.");
  }
}

inline double Leaf::evaluate_psi_stem_at(double target_psi_stem) {
  require_stem_route("evaluate_psi_stem_at");
  return with_curve(cost_curve_, [&](auto tag) {
    return evaluate_psi_stem<tag.value>(target_psi_stem);
  });
}

// Returns (value, feasible), as the collar version does: `dprofit`'s 0.0 on its
// shut-down and reversed exits is indistinguishable from a stationary point.
inline std::vector<double> Leaf::dprofit_dpsi_stem_checked(double psi_stem) {
  require_stem_route("dprofit_dpsi_stem_checked");
  if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }
  const double psi_upstream = supply_psi_soil_scalar();
  bool feasible = false;
  const double v = with_curve(cost_curve_, [&](auto tag) {
    return dprofit_dpsi_stem<tag.value>(psi_stem, psi_upstream, &feasible);
  });
  return std::vector<double>{v, feasible ? 1.0 : 0.0};
}


// Which link each curve composes with assimilation. One table, and the last arm
// asserts, so a curve added to the enum without a link is a build failure.
template <Leaf::CostCurve K>
constexpr Leaf::BenefitLink Leaf::benefit_link() {
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26 ||
                K == CostCurve::LeastCost) {
    // SOX and JW26 are products `A*g`, LeastCost a ratio `A/W`: the link is the
    // same `h = log`, and the sign (-log g against +log W) lands in `cost_deriv`.
    return BenefitLink::Log;
  } else if constexpr (K == CostCurve::ProfitMax) {
    return BenefitLink::Scaled;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::CF77 ||
                  K == CostCurve::JS22 || K == CostCurve::CMax ||
                  K == CostCurve::TF24_floor,
                  "unhandled CostCurve in benefit_link");
    return BenefitLink::Identity;
  }
}

// `h'(A)`, all the derivative needs from the link. The Log arm differentiates
// `log A + log g`, which shares `A*g`'s argmax. `log A` needs `A > 0`, and at
// full closure `A = -R_d`: a non-positive `A` returns NaN, which the caller's
// feasibility flag carries outward.
template <Leaf::CostCurve K>
inline double Leaf::benefit_link_deriv(double A) const {
  if constexpr (benefit_link<K>() == BenefitLink::Log) {
    return A > 0.0 ? 1.0 / A : util::na_value;
  } else if constexpr (benefit_link<K>() == BenefitLink::Scaled) {
    return profitmax_A_max_ > 0.0 ? 1.0 / profitmax_A_max_
                                   : util::na_value;
  } else {
    (void)A;
    return 1.0;
  }
}


// --- the three things that vary, each in ONE place ---------------------------
//
// A new single-layer cost curve is a `CostCurve` member plus one arm in each of
// these three. Every arm is explicit and the last asserts, so a curve forgotten
// here is a compile error, not a fall-through to the `else`.

// What must be set before the search. TF24's parameters are traits with
// defaults; SOX and JW26 take none.
template <Leaf::CostCurve K>
inline void Leaf::check_cost_parameters() {
  if constexpr (K == CostCurve::CF77) {
    // An input with an aligned default, so non-finite means a caller cleared it:
    // a NaN objective's argmax is a property of the bracket and looks plausible.
    if (!std::isfinite(CF77_lambda_)) {
      util::stop("the CF77 cost curve needs CF77_lambda_ set: it is the "
                 "PRESCRIBED marginal value of water in umol C (kg H2O)^-1, and "
                 "it was set to " + util::to_string(CF77_lambda_) + ". Its "
                 "default is aligned to TF24 (see align_curves()); "
                 "set_physiology and set_traits never touch it.");
    }
    // Only with the option on, so a caller who never sets it is not refused for
    // soil values the constant never reads. beta <= 0 prices water at infinity
    // or below zero, which pays the leaf to transpire to the wet bound.
    if (CF77_soil_beta_) {
      if (!(theta_fc_ > theta_w_)) {
        util::stop("CF77_soil_beta_ needs theta_fc_ > theta_w_; got theta_fc_ = " +
                   util::to_string(theta_fc_) + ", theta_w_ = " +
                   util::to_string(theta_w_));
      }
      const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
      if (!std::isfinite(beta) || beta <= 0.0) {
        util::stop("CF77_soil_beta_ is on and the soil-moisture factor beta = "
                   "(theta_ - theta_w_)/(theta_fc_ - theta_w_) is " +
                   util::to_string(beta) + ", which prices water at infinity or "
                   "below zero. theta_ must exceed theta_w_; got theta_ = " +
                   util::to_string(theta_) + ", theta_w_ = " +
                   util::to_string(theta_w_));
      }
    }
  } else if constexpr (K == CostCurve::JS22) {
    // A trait with a default, still checked: a NaN maximises a NaN objective.
    if (!std::isfinite(JS22_gamma) || JS22_gamma < 0.0) {
      util::stop("JS22_gamma must be a finite, non-negative hydraulic unit cost "
                 "in umol C m^-2 s^-1 MPa^-2; got " + util::to_string(JS22_gamma));
    }
  } else if constexpr (K == CostCurve::CMax) {
    // Only CMax_a gets the sign check: CMax_b is signed, negative in the paper's
    // own convention (see the member).
    if (!std::isfinite(CMax_a) || CMax_a < 0.0) {
      util::stop("CMax_a must be a finite, non-negative slope for the CMax "
                 "marginal cost in umol C m^-2 s^-1 MPa^-2; got " +
                 util::to_string(CMax_a));
    }
    if (!std::isfinite(CMax_b)) {
      util::stop("CMax_b must be finite; got " + util::to_string(CMax_b));
    }
  } else if constexpr (K == CostCurve::TF24_floor) {
    // No default, unlike CF77_lambda_: at lambda_o = 0 this curve IS TF24, so an
    // unset price would run TF24 under this curve's name. An explicit zero is
    // accepted (the reduction test's case); the hydraulic half is TF24's traits.
    if (!std::isfinite(TF24_floor_lambda_o) || TF24_floor_lambda_o < 0.0) {
      util::stop("the TF24_floor cost curve needs TF24_floor_lambda_o set: it is the "
                 "PRESCRIBED marginal cost of water as transpiration goes to "
                 "zero, in umol C (kg H2O)^-1, non-negative, and NA until you "
                 "assign one -- set_physiology and set_traits never touch it. It "
                 "has no default because 0 would silently make this curve TF24. "
                 "Set `$TF24_floor_lambda_o` on the leaf, or pass "
                 "`TF24_floor_lambda_o =` to leaf_solve() / leaf_batch(); got " +
                 util::to_string(TF24_floor_lambda_o));
    }
  } else if constexpr (K == CostCurve::LeastCost) {
    // An aligned default, refused when cleared for CF77's reason. Strictly
    // positive, unlike TF24_floor_lambda_o: zero leaves `A/E`, which rises without
    // bound as E -> 0 and puts every optimum at the wet bound.
    if (!std::isfinite(LeastCost_beta) || LeastCost_beta <= 0.0) {
      util::stop("the LeastCost cost curve needs LeastCost_beta set: it is "
                 "Prentice et al. (2014)'s beta, the unit cost of maintaining "
                 "carboxylation capacity relative to that of transpiration, "
                 "DIMENSIONLESS (both fluxes are weighed in mol m^-2 s^-1), "
                 "strictly positive; its default is aligned to TF24 (see "
                 "align_curves()), and set_physiology and set_traits never touch "
                 "it. The published value is 146. Set `$LeastCost_beta` on the "
                 "leaf, or pass "
                 "`LeastCost_beta =` to leaf_solve() / leaf_batch(); got " +
                 util::to_string(LeastCost_beta));
    }
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::SOX ||
                  K == CostCurve::JW26 || K == CostCurve::ProfitMax,
                  "unhandled CostCurve");
  }
}


// Which objective. Not all one kind of number: five are carbon, two products
// (`profit_psi_stem_SOX`) and one a ratio (`profit_psi_stem_LeastCost`).
template <Leaf::CostCurve K>
inline double Leaf::profit_psi_stem_for(double psi_stem, double psi_upstream) {
  if constexpr (K == CostCurve::TF24) {
    return profit_psi_stem_TF(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::CF77) {
    return profit_psi_stem_CF77(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::JS22) {
    return profit_psi_stem_JS22(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::CMax) {
    return profit_psi_stem_CMax(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::SOX) {
    return profit_psi_stem_SOX(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::JW26) {
    return profit_psi_stem_JW26(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    return profit_psi_stem_TF24_floor(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::LeastCost) {
    return profit_psi_stem_LeastCost(psi_stem, psi_upstream);
  } else {
    static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
    return profit_psi_stem_ProfitMax(psi_stem, psi_upstream);
  }
}


// Which lambda to report, normalised to two arguments (TF24's reads psi_stem
// alone). CF77's is exact by construction: its cost is lambda*E, so dC/dpsi over
// dE/dpsi is lambda identically.
template <Leaf::CostCurve K>
inline double Leaf::lambda_for(double psi_stem, double psi_upstream) {
  if constexpr (K == CostCurve::TF24) {
    (void)psi_upstream;
    return lambda_TF24(psi_stem);
  } else if constexpr (K == CostCurve::CF77) {
    // Still independent of BOTH potentials -- the soil-beta option keys on theta,
    // a driver, not on where in the bracket the optimiser is.
    (void)psi_stem;
    (void)psi_upstream;
    return cf77_price();
  } else if constexpr (K == CostCurve::JS22) {
    return lambda_JS22(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::CMax) {
    return lambda_CMax(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::SOX) {
    return lambda_SOX(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::JW26) {
    return lambda_JW26(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // The hydraulic half reads psi_stem alone, as TF24's does.
    (void)psi_upstream;
    return lambda_TF24_floor(psi_stem);
  } else if constexpr (K == CostCurve::LeastCost) {
    return lambda_LeastCost(psi_stem, psi_upstream);
  } else {
    static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
    (void)psi_upstream;
    return lambda_ProfitMax(psi_stem);
  }
}


// lambda_emergent on the COLLAR route: dC/dE along the path the solve moves on,
// the collar psi_r free and psi_stem following. Uptake gives dE/dpsi_r = S and the
// stem path dpsi_stem/dpsi_r = (S + kmax f(psi_r)) / (kmax f(psi_stem)), so the
// part of a cost reading psi_stem alone takes the series correction
//
//     lambda = lambda_stem * (1 + kmax f(psi_r) / S),
//
// the part reading only E carries over unchanged, and a cost that also reads the
// collar (JS22's drop, CMax's lower limit) adds its partial over S. Each row must
// equal dA/dE at an interior optimum: test_collar_lambda_emergent_is_dA_dE.
template <Leaf::CostCurve K>
inline double Leaf::lambda_collar_for() {
  if constexpr (K == CostCurve::TF24) {
    return marginal_cost_water_multilayer();
  } else if constexpr (K == CostCurve::CF77) {
    return cf77_price();
  } else if constexpr (K == CostCurve::LeastCost) {
    return lambda_LeastCost(opt_psi_stem_, opt_root_psi_);
  } else {
    const double S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
    if (!std::isfinite(S) || S <= 0.0) {
      return util::na_value;
    }
    const double series =
        1.0 + leaf_specific_conductance_max_ *
                  proportion_of_conductivity(opt_root_psi_) / S;
    if constexpr (K == CostCurve::TF24_floor) {
      return lambda_TF24(opt_psi_stem_) * series + TF24_floor_lambda_o;
    } else if constexpr (K == CostCurve::JS22 || K == CostCurve::CMax) {
      return lambda_for<K>(opt_psi_stem_, opt_root_psi_) * series +
             cost_deriv_upstream<K>(opt_psi_stem_, opt_root_psi_) / S;
    } else if constexpr (K == CostCurve::SOX || K == CostCurve::JW26) {
      return lambda_for<K>(opt_psi_stem_, opt_root_psi_) * series;
    } else {
      static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
      // The conductance term reads psi_stem alone; the thermal term, where the
      // gates put one in, reads the leaf temperature and so E alone.
      const double f = proportion_of_conductivity(opt_psi_stem_);
      if (!(f > 0.0) || !(profitmax_k_span_ > 0.0) || !(profitmax_A_max_ > 0.0)) {
        return util::na_value;
      }
      const double abs_fprime = (stem_c / stem_b) *
                                pow(opt_psi_stem_ / stem_b, stem_c - 1.0) * f;
      double thermal = 0.0;
      if (use_thermal_cost_ && use_energy_balance_) {
        double dT_dE = 0.0;
        const double Tleaf = leaf_temp_from_E(
            operating_point_flux(opt_psi_stem_, opt_root_psi_), &dT_dE);
        thermal = thermal_cost_temp_deriv(Tleaf) * dT_dE;
      }
      return profitmax_A_max_ *
             (abs_fprime / (f * profitmax_k_span_) * series + thermal);
    }
  }
}


// The body, written once, over the CLOSED interval (hazard 11).
//
// ============================================================================
// WHAT EACH CURVE DOES DIFFERENTLY IN HERE
// ----------------------------------------------------------------------------
// TF24 -- the reference; wet-pinning is its dominant pinned class.
//
// CF77 -- the objective is the only difference. At a soil drier than psi_crit
//   the objective is EVALUATED at the no-flow point, not zeroed: profit is -R_d.
//
// JS22 -- no second hump: its marginal cost rises from zero, so the first-order
//   condition crosses once (hazard 11's double crossing needs a marginal cost
//   that falls at the dry end), but the endpoints are still candidates. The wet
//   end is interior whenever there is a bracket (dC/dpsi -> 0 there, dA/dpsi > 0);
//   with none, psi* is psi_soil exactly through the no-flow branch, not a pin.
//
// CMax -- can be wet-pinned: its marginal cost tends to CMax_a*psi_soil + CMax_b,
//   not zero. A pinned answer is not psi_soil exactly but the bracket's step-in,
//   which does not move with the parameters: detect closure against the bound
//   with the step-in tolerance, or by the gradient's sign there.
//
// SOX -- a product, never dry-pinned since g(psi_crit) = 0; the wet bound is
//   reachable when A is negative throughout (hazard 11). No parameter to
//   validate. With C(psi) = -log g and h = log, C reads psi only, so the
//   first-order conditions have identical roots:
//
//       h(A) - C          = log A + log g = log(A*g)
//       d/dpsi[h(A) - C]  = A'/A + g'/g   = (A'*g + A*g') / (A*g) = P' / (A*g)
//   The optimiser still maximises the product: log A does not exist at full
//   closure, the endpoint hazard 11 requires.
//
// JW26 -- SOX's objective with a linear reduction factor, on this package's
//   hydraulics (psi_crit derived, supply the curve's integral): not a
//   reproduction of the paper.
//
// ProfitMax -- needs setup before the search, and its degenerate exit ZEROES
//   where every other curve EVALUATES (the branch below).
//
// TF24_floor -- TF24 at TF24's own traits plus lambda_o*E, so its marginal cost
//   does not vanish at the wet end: it can be wet-pinned as CMax can, with the
//   same step-in. It inherits TF24's second hump, so basin_scan_cells() matters
//   as for TF24. Under kmax/D rescaling its optimum sits between TF24 and CF77
//   (test_tf24_floor_kmax_vpd_invariance).
//
// `profit_` carries three kinds of number and nothing enforces it: carbon for the
// difference objectives, a normalised profit for ProfitMax, carbon times a factor
// for the products. Compare argmaxes across curves, never `profit_`.
// ============================================================================
template <Leaf::CostCurve K>
inline void Leaf::optimise_psi_stem_single() {
  clear_collar_solve_state();

  if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }
  check_cost_parameters<K>();

  const double psi_soil = supply_psi_soil_scalar();
  // The wet end: where the stem carries exactly the residual flux, psi_soil
  // itself at g_min = 0.
  const double psi_wet = stem_wet_bound(psi_soil);
  opt_psi_stem_ = psi_wet;

  // The dry end: psi_crit for a curve whose own domain ends there, the end of
  // the stem curve for every other.
  const double psi_dry = dry_end<K>();

  // ProfitMax divides by |A|max, so it is seeded here, and with unusable
  // normalisers there is no objective at the wet end to evaluate. A dark leaf is
  // tested first, as its |A|max is non-positive anyway. Its normalisers are
  // Sperry's, over the stream to psi_crit, so a wet end past psi_crit has none.
  if constexpr (K == CostCurve::ProfitMax) {
    if ((PPFD_ < 1.5e-8) | (psi_wet > psi_crit)) {
      set_profitmax_degenerate_state(psi_wet, psi_soil);
      return;
    }
    prepare_profitmax();
    if (!(profitmax_k_span_ > 0.0) || !(profitmax_A_max_ > 0.0)) {
      set_profitmax_degenerate_state(psi_wet, psi_soil);
      return;
    }
  } else if (psi_wet > psi_dry || !(stem_curve_complement(psi_soil) > 0.0)) {
    // Drier soil than the stem can carry more than the cuticle's loss from
    // short of its dry end, or than it can carry anything from: one feasible
    // potential, and this is it.
    profit_ = profit_psi_stem_for<K>(psi_wet, psi_soil);
    lambda_emergent_ = lambda_for<K>(psi_wet, psi_soil);
    return;
  }
  // A product objective with A <= 0 at the dry end has no maximum, and the leaf
  // is shut at the wet end (objective_is_product).
  if constexpr (objective_is_product<K>()) {
    if (dark_under_product<K>(psi_dry, psi_soil)) {
      profit_ = profit_psi_stem_for<K>(psi_wet, psi_soil);
      lambda_emergent_ = lambda_for<K>(psi_wet, psi_soil);
      return;
    }
  }

  double profit_opt = 0.0;
  const auto objective = [&](double psi_stem) {
    return profit_psi_stem_for<K>(psi_stem, psi_soil);
  };
  // One solver, with no caller-side branch on whether a derivative exists: the
  // closure reports that through `ok` or a non-finite value and the solver falls
  // back. A missing derivative must be reported, never thrown, or the whole
  // solve is lost. basin_scan_cells() scans only where two basins can occur.
  const auto dobjective = [&](double psi_stem, bool* ok) -> double {
    return dprofit_dpsi_stem<K>(psi_stem, psi_soil, ok);
  };
  // For a curve psi_crit does not bound, psi_crit is still tried first, as on
  // the collar route: where profit already falls there the optimum lies inside
  // the narrower bracket. Not under a basin scan, whose cells the interval sets.
  double hi = psi_dry;
  if constexpr (!cost_domain_ends_at_psi_crit<K>()) {
    if (basin_scan_cells() == 0 && psi_wet < psi_crit) {
      bool ok = false;
      const double d = dobjective(psi_crit, &ok);
      if (ok && std::isfinite(d) && d < 0.0) {
        hi = psi_crit;
      }
    }
  }
  // A scan keeps its cell width, not its count: past psi_crit the interval runs
  // to the end of the stem curve, and wider cells step over a dark leaf's narrow
  // basin near its wet end. So the count grows with the interval, to at most 4x.
  int cells = basin_scan_cells();
  if (cells > 0 && hi > psi_crit && psi_crit > psi_wet) {
    cells = std::min(4 * cells,
                     static_cast<int>(std::ceil(cells * (hi - psi_wet) /
                                                (psi_crit - psi_wet))));
  }
  opt_psi_stem_ = util::maximise_over_closed_interval_foc(
      objective, dobjective, psi_wet, hi, cells, collar_root_tol,
      static_cast<size_t>(ci_niter), &profit_opt);
  // An argmax at psi_dry for a curve that is not bounded at psi_crit is the end
  // of the stem curve: the leaf transpires all the stem can carry, as on the
  // collar route (see maximise_profit_over_collar).
  profit_ = profit_psi_stem_for<K>(opt_psi_stem_, psi_soil);
  lambda_emergent_ = lambda_for<K>(opt_psi_stem_, psi_soil);
  (void)profit_opt;
}


// ===========================================================================
// MEDLYN STOMATAL-CONDUCTANCE MODEL
// ---------------------------------------------------------------------------
// Standalone, R-callable coupling of the Medlyn (2011) optimal stomatal model
// to colimited photosynthesis, not on the optimality routes. beta_ is a
// soil-moisture factor in [0,1] from theta_, theta_w_ and theta_fc_.
// ===========================================================================
inline double Leaf::medlyn_model_gs(double assim_colimited_){

  const double beta_ = medlyn_beta();

  // The USO expression as published: a WATER VAPOUR conductance, which `g1` is
  // fitted against. The zero-deficit guard reads the driven `atm_vpd_`, the value
  // the expression reads, not the `atm_vpd` input field.
  double gs_h2o;
  if (atm_vpd_ == 0) {
     gs_h2o = g0;
  } else {
     gs_h2o = g0 + medlyn_H2O_CO2_ratio*(1 + (g1*beta_)/sqrt(atm_vpd_))*(assim_colimited_/(ca_*(1/umol_per_mol_to_Pa_)));
  }

  // ...converted to the CO2 basis every other conductance here is on, so
  // `A == stom_cond_CO2_*(ca_ - ci_)` holds on this route too; `g0` included,
  // since Medlyn's `g0` is a water-basis conductance as well.
  medlyn_model_gs_ = gs_h2o / medlyn_H2O_CO2_ratio;
  return medlyn_model_gs_;
}

// Supply==demand residual for the Medlyn solver in ci (x, Pa), multiplied
// through by (ca_ - x) so it stays finite at ca_; for x < ca_ its root is the
// raw conductance difference's:
//   gs_medlyn*(ca_-x) - gs_coupled*(ca_-x),  where
//   gs_coupled*(ca_-x) = assim * (atm_kpa_*kPa_to_Pa) / 1e6.
// Both terms are on the CO2 basis; scaling one alone would move the root.
inline double Leaf::medlyn_stom_cond_minus_coupled_stom_cond(double x) {
  const double assim_colimited_x_ = assim_colimited(x);
  medlyn_model_gs_ = medlyn_model_gs(assim_colimited_x_);
  return medlyn_model_gs_ * (ca_ - x)
         - assim_colimited_x_ * (atm_kpa_ * kPa_to_Pa) / 1e6;
}

// The ci at which the Medlyn conductance balances the diffusion-implied one. The
// residual is negative at both ends of [gamma*, ca_] with a hump between, so it
// can hold a spurious sub-compensation root: locate its maximum first, then
// root-find on [argmax, ca_] for the high-ci operating point (g0 == 0 included).
inline void Leaf::solve_medlyn_ci_numerical(){
  auto target = [&](double x) -> double {
    return medlyn_stom_cond_minus_coupled_stom_cond(x);
  };
  const double lo = gamma_ * umol_per_mol_to_Pa_;
  const double hi = ca_;

  double neg_peak = 0.0;
  const double ci_peak =
      util::brent_fmin([&](double x) { return -target(x); }, lo, hi, ci_abs_tol,
                       &neg_peak);
  const double residual_peak = -neg_peak;

  if (residual_peak <= 0.0) {
    // Supply never reaches demand: no feasible Medlyn operating point. Report
    // the closest approach (the residual maximum) rather than failing.
    ci_ = ci_peak;
  } else {
    try {
      ci_ = util::uniroot(target, ci_peak, hi, ci_abs_tol, ci_niter);
    } catch (const std::exception& e) {
      util::stop_infeasible("ci_solve", "solve_medlyn_ci_numerical failed: " + std::string(e.what()) +
                 "; ci_peak=" + util::to_string(ci_peak) +
                 "; max=" + util::to_string(hi));
    }
  }
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  set_medlyn_shutdown_if_closed();
  return;
}

// Medlyn's closed form for chi = ci/ca: the numerical route's limit as g0 -> 0,
// dropping g0, Gamma* and respiration. beta multiplies g1 here as in
// `medlyn_model_gs`, so the two routes are one model.
inline void Leaf::solve_medlyn_ci_analytical(){

  const double beta_ = medlyn_beta();
  const double g1_beta = g1 * beta_;
  ci_ = ca_ * (g1_beta/(g1_beta + sqrt(atm_vpd_)));
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  set_medlyn_shutdown_if_closed();
  return;
}

// The linear soil-moisture factor, clamped to [0, 1]: above field capacity it
// would scale g1 past its fitted value, below the wilting point make it negative.
// CF77's soil price clamps it at 1 only, since it divides by it.
inline double Leaf::medlyn_beta() const {
  const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
  return std::min(std::max(beta, 0.0), 1.0);
}

// A negative conductance (the USO `A/ca` term below the compensation point
// outrunning `g0`) is not an operating point, so it becomes a closed state: no
// water, ci at the compensation point, A = -R_d_. That is not the optimality
// routes' shut state (the evaluator at the wet bound), so in light this route
// steps to -R_d_ as it shuts. Clamped, not refused: a hot or dry leaf reaching
// it is a legitimate state.
inline void Leaf::set_medlyn_shutdown_if_closed(){
  if (medlyn_model_gs_ < 0.0 || stom_cond_CO2_ < 0.0) {
    stom_cond_CO2_ = 0.0;
    medlyn_model_gs_ = 0.0;
    transpiration_ = 0.0;
    ci_ = gamma_ * umol_per_mol_to_Pa_;
    assim_colimited_ = -R_d_;
  }
}

} // namespace phylloptim

#endif

