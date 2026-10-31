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
// left in members (opt_root_psi_, opt_psi_stem_, ci_, profit_, E_up_, ...).
// The eight cost curves share that solve and differ only in their cost (see
// CostCurve).
//
// THE DERIVATIVE SURFACE. A consumer differentiating a model that contains this
// leaf takes the operating point's derivative from three member templates, at
// whatever scalar it holds, after the solve has run at double:
//
//   SupplyDraw<S>  draw = leaf.supply_draw_at<S>(S(leaf.opt_root_psi_), supply);
//   HeldCoords<S>  held;
//   S              psi  = leaf.operating_collar<K, S>(draw, pars, curvature, &held);
//   BilledOutputs<S> out  = leaf.billed_outputs_at<K, S>(psi, draw, pars, &held);
//
// `out.profit` and each layer's `out.uptake[j]` then carry their derivatives in
// `pars` and `supply`, supplied from the implicit function theorem at the
// converged point; nothing the root-find did is recorded. S is double, a tangent
// (odelia::ode::tangent_scalar) or an adjoint (odelia::ode::adjoint_scalar).
// `held` is optional: it carries the stem potential and ci operating_collar solved at
// an interior point to billed_outputs_at, which otherwise solves them again. Without
// it every value and tangent is the same bit for bit, and every adjoint row the
// same to rounding, because it is summed in another order. Preconditions,
// checked where marked:
//
//   - the draw is taken at the solved collar (checked: operating_collar and billed_outputs_at
//     refuse another);
//   - `held`, if passed, comes from operating_collar at this point (checked: billed_outputs_at
//     refuses one placed at another collar or kind);
//   - nothing moves the leaf between the solve and billed_outputs_at;
//   - `pars` and `supply` hold, in value, the seated traits (passive_pars()) and
//     the seated soil state and root network;
//   - K is the curve the solve used, TF24 or TF24_floor (checked against the
//     seated curve; see check_surface_curve for the one solve it cannot see);
//   - the supply path is multi-layer (checked: supply_draw_at refuses another);
//   - `curvature` is profit_curvature_at_collar<K>() at an Interior point, and is
//     not read at any other kind or at double; operating_collar refuses one that is
//     not negative or lies inside Leaf::interior_curvature_floor;
//   - the thermal-damage layer is off: operating_collar and billed_outputs_at refuse
//     use_thermal_damage_, because the kernels they read do not apply it.
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
// A recorded point is put back with replay_operating_point(collar, kind), which
// hands over the same rows the solve did. profit_row_at_current_point in
// tests/cpp/test_leaf.cpp is a complete, compiled use of the three calls; the
// vignette "Differentiating the operating point" explains the design.
//
// "Hazard N" in comments below refers to the numbered hazards in the
// repository's .claude/CLAUDE.md, which is not shipped with the package.

namespace phylloptim {

// --- the parameter enumeration, which R indexes into --------------------------
//
// The sixteen traits in `Leaf::set_traits`' argument order, then the slots that
// are not traits: `kmax` (the conductance driver), `resistance` (the
// single-potential path's series resistance), the three model-owned slots
// (`CF77_lambda`, `TF24_floor_lambda_o`, `LeastCost_beta`), the fourteen inputs of
// the thermal-damage layer that reach the solve, and `PPFD`. Everything before
// `PPFD` is what a calibration fits, `gradient::n_theta` (35); `PPFD` is seated
// per observation from the drivers and is never fitted, and `n_leaf_pars` (36),
// the length of the pack the kernels read, counts it. Size `theta` by
// `n_theta`; read it out to `n_leaf_pars` and you run past its end.
//
// ⚠️ R INDEXES THESE POSITIONS, so a reordering silently differentiates the wrong
// parameter. `test-gradient-batch.R` reads the names back out of C++ and compares
// them with R's, so the two cannot drift apart without a failure.
inline constexpr int n_traits = 16;
inline constexpr int n_leaf_pars = 36;

// Every index by name, so nothing below indexes `theta` with a bare integer.
// The first `n_traits` are `set_traits`' arguments in its order, which is also
// `leaf_traits()`'; the non-traits follow and take a relative step.
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
// The residual (cuticular) conductance to water vapour at 25 C, mol H2O m^-2
// s^-1 (#182). The sixteenth trait and the last one: it is appended to the
// trait block so every slot after it moves by one, which is the renumbering
// the comment below describes and every consumer survives by reading the pack
// by name.
inline constexpr int par_g_min = 15;
// ⚠️ THESE MOVE WHENEVER A TRAIT IS ADDED, and bumping them is the whole cost.
// They are the non-traits and they sit AFTER the contiguous trait block, which
// is a readability convention rather than a constraint: R's
// `.gradient_theta_matrix()` addresses EVERY column by name, including these.
// What IS load-bearing is the ORDER ITSELF -- R passes integer positions into this
// enumeration, so reordering silently differentiates the wrong parameter, and
// `test-gradient-batch.R` compares this enumeration against R's copy. A new
// FITTED slot goes before `par_PPFD`, not after it: `gradient::n_theta` is
// `par_PPFD`, so a slot appended last is outside every caller's theta.
inline constexpr int par_kmax = 16;
inline constexpr int par_resistance = 17;
// Cowan-Farquhar's prescribed marginal value of water. A pure APPEND after the two
// existing non-traits, which is only safe because R addresses theta's non-trait
// columns by NAME rather than by position -- a positional rule ("everything but
// the last two") reads this as `resistance`.
//
// ⚠️ AVAILABLE FOR ONE MODEL. It is CF77's only parameter and every other curve's
// lambda is EMERGENT, derived from that curve's own parameters rather than set. So
// `.gradient_available_pars()` offers it only for CF77 and refuses it elsewhere,
// naming the model -- the same treatment `resistance` gets on the wrong supply path.
inline constexpr int par_CF77_lambda = 18;
// TF24_floor's price of water at zero transpiration, on exactly the same footing:
// an append after the non-traits, available for ONE model, and refused elsewhere
// by `.gradient_available_pars()` with the model named.
//
// ⚠️ IT IS THE SECOND MODEL-SPECIFIC SLOT, so "the CF77 one" has stopped being a
// safe way to talk about this class. R's `.gradient_model_pars()` is the single
// table that says which model owns which slot; there is no second copy here.
inline constexpr int par_TF24_floor_lambda_o = 19;
// Least-cost's unit-cost ratio (Prentice et al. 2014's beta), the THIRD
// model-owned slot. Same treatment as the two prices above:
// `.gradient_model_pars()` offers it for `LeastCost` alone and refuses it
// elsewhere, naming the model. Inserted here rather than appended after the
// thermal-damage block so the three slots one curve each owns stay contiguous,
// which costs the slots below a renumbering -- safe, because every
// consumer reads them by name (plant's `tf24_strategy.h` included).
//
// ⚠️ IT IS NOT A PRICE, unlike the other two model-owned slots, and the message
// R raises for it says so. `CF77_lambda_` and `TF24_floor_lambda_o` are marginal
// values of water in umol C (kg H2O)^-1; this is a DIMENSIONLESS ratio of two unit
// costs, and least-cost's price of water is emergent from it rather than equal to
// it. What the three share is the reason they live here: a parameter of one curve
// that every other curve has no slot for.
//
// ⚠️ `gradient::step_for`'s RULE IS RIGHT FOR IT BY LUCK RATHER THAN BY DESIGN.
// Its published value is 146, so the `max(|value|, 1)` floor never binds and it
// takes a relative step. A caller sweeping beta below 1 would get the absolute
// step instead, which at beta = 0.5 is a 2x perturbation -- the same trap
// `leaf_specific_conductance_max` is exempted from, and not exempted here because
// no such use exists yet.
inline constexpr int par_LeastCost_beta = 20;
// The thermal-damage layer's inputs that REACH THE SOLVE (use_thermal_damage_):
// its three states, the acclimated shift of the Vcmax/Jmax optima, and every
// quantity a maintenance cost multiplies. Model-owned in the same sense as the
// three slots above -- available only on a leaf with the gate on, and refused by
// name elsewhere -- and seated BEFORE the drivers, because set_physiology is what
// reads them.
//
// ⚠️ THE LAYER'S DYNAMICS PARAMETERS ARE NOT HERE, AND MUST NOT BE ADDED. The
// impairment rate's Arrhenius (k_i_Ea_, k_i_Tiso_, dTimpair_max_), maturation
// (k_mat_, k_mat_Ea_) and the acclimation kinetics move only damage_rates(),
// never a solved output, so a slot here would report an exact zero that reads as
// an answer. Their gradient is the trajectory's (damage_trajectory_gradient).
// Repair is here because the repair-flux cost reads k_rec_eff(T), so both of its
// parameters reach R_d_.
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
// The light the leaf is standing in: the INTERCEPTED PPFD per unit leaf area,
// the photon flux incident on the leaf before its absorptance. A DRIVER rather
// than a trait, on the same footing as `kmax`: the consumer seats it per cohort
// and it carries that consumer's rows.
//
// ⚠️ IT IS THE ONLY ROUTE LIGHT HAS TO CARBON. Read off the `PPFD_` member the
// forward pass seats, the electron transport is a constant in it, so every trait
// that moves the light a cohort stands in -- the extinction coefficient, the leaf
// area a stem carries -- reaches assimilation with no row at all and reads as an
// exact zero, which is the signature of a missing accumulator.
inline constexpr int par_PPFD = 35;
static_assert(n_leaf_pars == par_PPFD + 1,
              "the pack ends at par_PPFD; move n_leaf_pars with it");


// The differentiable parameters, as one array indexed by the enumeration above.
//
// THIS IS THE INPUT PACK, and it is deliberately not a struct of named fields.
// A struct would be a fourth spelling of an order that set_traits, leaf_traits()
// and par_names() already spell three times, and every field would have to be
// listed again in a walk that sets directions -- which is how a parameter gets added to
// the model and silently omitted from its gradient. Indexed by `par_*` the
// access reads the same (`pars[par_kmax]`, not `pars[15]`) and there is one list.
//
// ⚠️ IT CARRIES THE TRAITS, NOT THE QUANTITIES DERIVED FROM THEM. stem_b and
// psi_crit are functions of (stem_P50, stem_c) and are derived where they are
// used, so a direction in stem_P50 reaches them by the chain rule rather than by a
// caller remembering to move three numbers together. Passing the derived pair
// instead is what lets a curve be stated twice and the two disagree.
template <class S> using leaf_pars = std::array<S, n_leaf_pars>;

// THE LAYER'S DYNAMICS, written ONCE as templates on the scalar type, so the
// rates the solve's block reports, damage_rates(), damage_trajectory() and the
// trajectory's exact sensitivities (damage_trajectory_gradient(), on a forward
// dual number) are one set of expressions. At double each is the same expression,
// in the same order, as the member it replaced, so nothing at double moves.
namespace damage {

// A forward-mode dual number: a value and its derivative along one input.
// Only what the dynamics need. Comparisons read the value.
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

// The dynamics parameters, by name. ⚠️ R's leaf_damage_gradient() names these by
// POSITION through `Leaf::damage_par_names()`, and `test-thermal-damage.R` reads
// the names back, so a reordering fails there rather than differentiating the
// wrong parameter. The three initial states follow them in a gradient request.
enum : int {
  dp_k_i, dp_k_i_Ea, dp_k_i_Tiso, dp_k_rec, dp_k_rec_Ea, dp_k_mat, dp_k_mat_Ea,
  dp_acclim_alpha, dp_acclim_beta, dp_acclim_T, dp_acclim_s, dp_dTimpair_max,
  dp_K_A, dp_topt_offset,
  n_dyn_pars
};
template <class T> using dyn_pars = std::array<T, n_dyn_pars>;

// The reference temperature of the repair and maturation rates, deg C: the
// temperature at which k_rec_ and k_mat_ are stated. Impairment has its own,
// the isokinetic temperature k_i_Tiso_.
inline constexpr double damage_Tref = 25.0;

// The smooth max(0, x) at sharpness s, ln(1 + e^(s x)) / s, with no clamp on
// its argument: below zero e^(s x) underflows to a finite zero, and above 30 the
// form that cannot overflow takes over. Up to 30 it is the plain form, which
// keeps every value the argument clamp to [-30, 30] used to produce, and past
// either end the function rather than a constant: the clamp held the forcing
// fixed below acclim_T - 30/s, which a leaf below 0 C reaches at the defaults.
template <class T> inline T softplus(const T& x, const T& s) {
  const T z = s * x;
  return (val(z) > 30.0 ? z + log1p(exp(-z)) : log1p(exp(z))) / s;
}

// Every rate in the layer: k at the reference temperature Tref_K (kelvin), on a
// centred Arrhenius with activation energy Ea (J mol^-1). At Ea = 0 the argument
// is an exact zero and the rate exactly k.
//
// ⚠️ NO CLAMP ON THE ARGUMENT. It used to be held to [-30, 30], and the lower
// end binds in ordinary use -- a high-Ea species (619 kJ mol^-1) below about
// 25 C, any species near 0 C -- where it held the rate at k e^-30 and its
// derivative at zero. e^arg underflows harmlessly, so the rate is the function
// itself; only an overflow, which needs an absurd temperature or Ea, is refused.
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

// The acclimation shift of the impairment's isokinetic temperature, read off A
// itself. A >= 0 is invariant under the dynamics (the forcing is a softplus and
// the decay is linear), so this used to read max(A, 0) and never bound; without
// it the shift is smooth through zero, which an integrator's trial stage may
// cross by a hair. Its one singularity, A = -K_A, is refused where a state is
// seated.
template <class T> inline T delta_impair(const T& A, const dyn_pars<T>& p) {
  return p[dp_topt_offset] + p[dp_dTimpair_max] * A / (A + p[dp_K_A]);
}

// Impairment, functional -> impaired: Bison & Michaletz (2026)'s first-order
// PSII inactivation, written about the isokinetic temperature their
// Meyer-Neldel relation makes common to every species, so k_i_ is the rate there
// and k_i_Ea_ alone sets the steepness. Unbounded in temperature by design.
//
// The isokinetic temperature, shifted, must lie in [250, 350] K: it used to be
// clamped there, and never was on any forcing measured, so outside it is
// refused rather than held.
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

// RK4 over a piecewise-linear leaf temperature; calls `put(i, y)` at every time.
// The sub-step count is chosen from VALUES, so on a dual run it is the same
// count as on the double run and the sensitivity is that of the discrete map
// the double run evaluates.
template <class T, class Put>
inline void trajectory(const std::vector<double>& times,
                       const std::vector<double>& leaf_temp,
                       std::array<T, 3> y, const dyn_pars<T>& p, Put put) {
  const std::size_t n = times.size();
  put(0, y);
  for (std::size_t i = 1; i < n; ++i) {
    const double t0 = times[i - 1], dt_all = times[i] - t0;
    const double T0 = leaf_temp[i - 1], T1 = leaf_temp[i];
    // The stiffest rate on the interval bounds the sub-step. Both ends, because
    // every rate is monotone in temperature; impairment at the smallest shift
    // its isokinetic temperature can take, because it falls as that rises.
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
    // Two bounds on the sub-step: the fastest rate times the step below 0.1, and
    // at most 0.25 K of temperature change per step, because the impairment rate
    // is steep in temperature and a forcing swept across it in one step is
    // under-resolved however slow the rates are.
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

  // The stem's cumulative vulnerability integral G, as a table of H(x) at
  // stem_c -- see VulnerabilityTable for why the table is in x and why it is C2.
  // It depends on stem_c alone, so moving stem_P50 (and with it stem_b) rebuilds
  // nothing; and it is immutable and shared, so copying a Leaf copies a pointer.
  // Read only through the four stem_curve_* accessors below.
  std::shared_ptr<const VulnerabilityTable> stem_curve_;
  InverseScale stem_inv_b_;

  // INSTRUMENTATION, not state: how many times the stem curve has been built.
  // Nothing in the model reads it, and `set_traits()` deliberately does NOT reset
  // it -- it counts events over an object's whole life, so a caller takes
  // differences.
  //
  // It exists because a rebuild is invisible to every other instrument: the table
  // is a pure function of (stem_c, knot count), so a rebuild at an unchanged pair
  // is bit-identical to not rebuilding, and the only thing it costs is time. This
  // is the countable form (#74), which is what the cost rules in the developer
  // guide ask for wherever one exists.
  long stem_curve_builds_ = 0;

  // The soil -> root-collar water supply (issue #2). Everything the leaf needs
  // from the soil enters through this object: `uptake` and its derivative. It is
  // held by value -- Leaf must stay copyable, because plant's
  // make_strategy_ptr(TF24_Strategy) takes the strategy by value and TF24_Strategy
  // holds a Leaf member.
  //
  // Public because plant's RcppR6 bindings reach the moved fields by name; they
  // now spell them `roots_.psi_soil_` etc.
  MultiLayerRoots roots_;

  // The alternative supply path (issue #2 stage 2/3). Both alternatives are held
  // as members and selected by `supply_kind_` -- measured free, where
  // std::variant costs +1.0%, because a predictable branch in front of an
  // already-out-of-line call disappears into it. ⚠️ Re-measure before "tidying"
  // this into a variant, and beware that a compile-time-known tag folds the
  // branch away and reports a false zero. SinglePotential is four doubles plus a one-element vector, so carrying it
  // unused costs tens of bytes per Leaf.
  SinglePotential single_;
  enum class SupplyKind { MultiLayer, SinglePotential };
  // Default MultiLayer: every existing caller, plant included, keeps today's
  // behaviour bit-for-bit without knowing this exists.
  SupplyKind supply_kind_ = SupplyKind::MultiLayer;

  // --- choosing the supply path (issue #32) ----------------------------------
  //
  // TWO ENTRY POINTS, NOT A SETTABLE TAG, and the difference matters. Assigning
  // supply_kind_ on its own leaves the other path's state configured and
  // silently ignored; assign it back and that state is now stale rather than
  // absent. This is the footgun to design around before
  // exposing any of it to R, where a settable field is the obvious thing to
  // reach for. Each of these leaves the object in a state where the tag and the
  // supply agree, and there is no intermediate state in which they do not.
  //
  // Both CLEAR the solved state, so set_physiology() must be called again
  // afterwards. That is not an inconvenience being papered over -- the two paths
  // read different inputs, so any state carried across would be answering a
  // question about the other model.
  void set_supply_multilayer() {
    supply_kind_ = SupplyKind::MultiLayer;
    single_.clear();
    setup_clean_leaf();
  }

  // ⚠️ THIS TAKES NO RESISTANCE, and that is deliberate. The soil-to-collar
  // resistance is a per-call DRIVER on
  // both supply paths now: it arrives through `set_physiology`, out of the same
  // `RootNetwork` the multi-layer path is given (see
  // SinglePotential::set_supply_resistances). Before, the multi-layer path took its
  // resistances per call and this one took its resistance at construction, so the
  // same quantity arrived at two different times depending on which path was in
  // force -- and `resistance` was the only differentiable parameter whose setter
  // reset the whole object, because it had to come back through here.
  //
  // `gravity_head` is the head to lift water to the collar in MPa, and IS still
  // configuration. That is the one asymmetry left, and it is not laziness: the
  // multi-layer path derives a per-layer head from the depth profile it is handed
  // (gravity_head * z_soil_mid), and this path has no depth profile to derive one
  // from. A bare leaf also wants zero rather than a geometric default, which the
  // multi-layer rule cannot express. A caller who does want the multi-layer rule
  // for one layer of thickness d passes `gravity_head = <gravity head> * d / 2`.
  void set_supply_single(double gravity_head = 0.0) {
    if (!std::isfinite(gravity_head) || gravity_head < 0.0) {
      util::stop("set_supply_single needs a finite, non-negative gravity_head in "
                 "MPa; got " + util::to_string(gravity_head));
    }
    supply_kind_ = SupplyKind::SinglePotential;
    setup_clean_leaf();
    // grav_head_ is configuration and clear() spares it, so this must come after
    // setup_clean_leaf. resistance_ used to need the same treatment and no longer
    // does -- it is a driver, clear() resets it, and set_physiology re-supplies it.
    single_.grav_head_ = gravity_head;
  }

  // Which path is in force, as a string, because the enum has no R
  // representation and a bare integer would be a worse one.
  std::string supply_kind_name() const {
    return supply_kind_ == SupplyKind::MultiLayer ? "multilayer" : "single";
  }

  // --- supply dispatch -------------------------------------------------------
  // The four points where the two paths differ. Everything else in the solve is
  // supply-agnostic and goes through the vector of signed potentials below,
  // which both paths provide -- that is what keeps this stage off the three
  // R-facing signatures that thread it (find_root_psi, find_psi_stem_from_psi_root,
  // E_from_Soil_to_Root_Collar).
  //
  // ⚠️ Those three take the soil state as an argument, and #25 changed what the
  // argument MEANS (positive suctions, not signed potentials) without changing
  // any signature. An R caller passing the old `-psi_soil` would get a silently
  // wrong answer, so each of them validates the vector is non-negative and stops
  // if not -- see require_suction_vector.
  double supply_begin_solve() {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.begin_solve();
      default:                     return single_.begin_solve();
    }
  }
  // The current soil state, as positive suction magnitudes. Threaded through
  // E_column, find_root_psi and find_psi_stem_from_psi_root.
  const std::vector<double>& supply_psi_soil() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.psi_soil_;
      default:                     return single_.psi_soil_vec_;
    }
  }
  // The wet bound's floor (MultiLayerRoots::wet_floor_layer): the lowest
  // gravity balance psi_i + z_i among the layers that still conduct, or the
  // wettest layer's where none does. The single potential's balance is its
  // soil plus its head, and it always conducts. `conducts` and `layer` report
  // which case and which layer, for the derivative surface.
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
  // Its partner (MultiLayerRoots::wet_ceiling): the highest balance among
  // conducting layers, where uptake is at least zero. Equal to the floor on
  // one layer and on the single potential.
  double supply_wet_ceiling() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.wet_ceiling(roots_.psi_soil_);
      default:                     return single_.psi_soil_ + single_.grav_head_;
    }
  }
  // The single soil potential the psi_soil_[0]-style solvers (optimise_psi_stem_*)
  // work against, positive magnitude. Those solvers already require exactly one
  // layer, so this is the same value either way -- it just stops them reaching
  // into MultiLayerRoots for it.
  // ⚠️ THIS CHECKS NOTHING, AND A COMMENT ELSEWHERE USED TO CLAIM IT DID. On a
  // multi-layer profile it returns LAYER ZERO's potential rather than refusing, so a
  // caller that needs one soil potential must test `supply_is_single_layer()`
  // itself. Relying on this to enforce the topology silently runs a single-potential
  // model on layer 0 of a layered soil -- which is exactly what happened when a
  // refusal was moved out of `prepare_profitmax`.
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

  // psi_from_E

  double vcmax_25;
  // NAMED stem_* deliberately. There are TWO Weibull vulnerability curves in this
  // model -- stem (stem_b, stem_c), which drives hydraulic_cost_TF, and root
  // (root_b, root_c), which drives uptake -- and these used to be the unmarked
  // default `b` and `c`. That is not a style question: the companion analysis used
  // the root parameters for the stem cost and carried lambda ~ psi^3.02 into a
  // manuscript draft where it should have been psi^0.64. Never leave an unmarked
  // default for a parameter that exists in two versions.
  double stem_c;
  // THE SCALE PARAMETER IS P50, the potential at 50% loss of conductivity, and it
  // is the trait. Written this way the curve is `2^-(psi/P50)^c`, the same
  // function as `exp(-(psi/b)^c)` but with scale and shape cleanly separated --
  // P50 is measured and reported, `b` (the 63% point) is not.
  //
  // ⚠️ Do not add `stem_b` back as a settable parameter alongside it. Both are
  // scale parameters in the same units, so specifying both means specifying two
  // points on one curve and recovering the shape from their ratio,
  // `c = ln(ln 2)/ln(P50/b)`, which diverges as P50 -> b. That is also why
  // vulnerability curves are fitted on a Px and the slope at Px rather than (b, c).
  double stem_P50;
  // DERIVED, and read-only from R. Kept as members because every internal use is
  // in terms of them: the curve reads, the brackets, the cost kernels.
  double stem_b;
  double psi_crit;
  double TF24_beta2;
  double jmax_25;
  double phi_J;
  double curv_fact_elec_trans; // unitless - obtained from Smith and Keenan (2020)
  double curv_fact_colim;
  // Despite its name, no longer the tolerance of the collar solve, which is a
  // root-find on dprofit = 0. It is the tolerance of the golden-section FALLBACK
  // maximise_profit_over_collar drops to when either bracket endpoint has no
  // usable gradient. It does NOT reach the stem route, which refines by a
  // root-find at collar_root_tol, and it no longer sets how well the reported
  // operating point is determined -- collar_root_tol does.
  double GSS_tol_abs;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
  double TF24_cost_scale;

  // Tolerance on the collar potential for the profit-maximising root-find.
  // Hard-coded rather than a control field, for the same reason
  // find_root_psi's and psi_stem_to_ci's tolerances are: it is a property of the
  // solve, not a knob. 1e-12 sits ~8 orders below the 1e-4 at which this package
  // calls a difference real, and the cost of going there from GSS_tol_abs's 1e-3
  // is a handful of evaluations, because TOMS748 is superlinear.
  static constexpr double collar_root_tol = 1e-12;
  // Tolerance on ci for every supply == demand root-find; psi_stem_to_ci says
  // why it is 1e-13.
  static constexpr double ci_root_tol = 1e-13;
  // Tolerance to which find_root_psi places the collar interval's two bounds,
  // and that tolerance at x, absolute plus relative.
  static constexpr double bracket_tol = 1e-4;
  static double bracket_tol_at(double x) {
    return 1.5 * bracket_tol * (1.0 + std::abs(x));
  }


  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
  std::vector<double> soil_consumption_;
  double profit_;
  // ⚠️ KEPT ONLY BECAUSE plant's RcppR6 bindings name it (hazard 7), and nothing
  // in this package writes or reads it. Initialised so a read through those
  // bindings is a defined NA rather than whatever the memory held; this package
  // no longer binds it.
  double psi_stem = util::na_value;
  // ⚠️ THE ONLY TWO INPUTS IN THIS BLOCK, and the reason they carry their default
  // here rather than in setup_clean_leaf(). Everything around them is derived
  // state or a solved output, which setup_clean_leaf() exists to wipe; these are
  // the Cowan-Farquhar marginal value of water, supplied by the CALLER and read
  // by profit_psi_stem_CF77. Wiping an input on the caller's behalf
  // makes whether a prescribed value survives depend on which of two
  // interchangeable-looking re-driving calls comes next, so neither clears it.
  //
  // An in-class initialiser, not a line in setup_clean_leaf(), so a fresh Leaf
  // reads the default: the member is otherwise uninitialised, and "never reset"
  // must not become "never initialised".
  // ⚠️ IT IS AN INPUT AND NOTHING ELSE. No model code assigns it. The units are
  // carbon per unit TRANSPIRATION, which is the Cowan-Farquhar marginal value of
  // water and the convention the lambda literature uses -- NOT the same quantity
  // as ProfitMax's normaliser ratio, which is per unit CONDUCTANCE.
  //
  // ALIGNED DEFAULTS (#156). Every cost curve with a level parameter defaults to
  // the value at which its stem optimum has TF24's conductance at one reference
  // environment: single-potential supply, psi_soil 1 MPa, PPFD 900, VPD 2 kPa,
  // ca 40 Pa, leaf at 25 C with the energy balance off, and every other trait at
  // its default. The curves then agree at the reference and differ only in how
  // their optimum moves away from it, which is what a comparison across curves
  // is after. TF24 is the anchor and does not move. Four parameters carry it:
  // `CF77_lambda_`, `JS22_gamma`, `CMax_a` (with `CMax_b` 0) and `LeastCost_beta`.
  //
  // ⚠️ THE FOUR NUMBERS ARE SOLVED, NOT CHOSEN, and R's `align_curves()` is the
  // solver. `test-align-curves.R` re-solves them and fails when a model change
  // moves an optimum at the reference; re-align then, rather than loosening it.
  // `leaf_traits()` and `.level_defaults` carry the same numbers on the R side.
  //
  // What is left unaligned, and why. SOX, JW26 and ProfitMax have no level
  // parameter: the shared supply sets their optimum, and aligning them through
  // `kmax` would break the one supply every curve shares. TF24_floor's hydraulic
  // half IS TF24's cost at TF24's traits, so the only lambda_o that reaches
  // TF24's conductance is zero, where the curve is TF24; it keeps no default.
  //
  // CF77's aligned price is exactly TF24's emergent lambda at the reference: at
  // that price the two first-order conditions coincide, so the optima do.
  double CF77_lambda_ = 107750.8846;            // umol C (kg H2O)^-1
  // Optional Medlyn (2011) soil-moisture factor on the CF77 price: with this on,
  // the marginal value of water becomes `CF77_lambda_ / beta(theta)` rather than
  // the bare constant, where beta is the same expression medlyn_model_gs uses.
  //
  // WHY THIS EXISTS. Cowan-Farquhar prices water at a rate that is invariant to
  // soil state BY CONSTRUCTION -- cost_deriv<CF77> discards psi_upstream, and
  // lambda_for<CF77> returns the bare constant -- so soil water reaches the
  // optimum only through the feasible bracket and the cost's VALUE, never through
  // its price. That makes CF77 the one curve here with no soil-moisture shutdown,
  // and this field is what supplies one without a second CostCurve enumerator.
  //
  // ⚠️ DEFAULT off, and off must stay BIT-IDENTICAL to the bare constant: the
  // golden files pin CF77 rows, so cf77_price() returns CF77_lambda_ itself on
  // this path rather than dividing it by a computed 1.0.
  bool CF77_soil_beta_ = false;
  // THE PRICE OF WATER AS TRANSPIRATION GOES TO ZERO, umol C (kg H2O)^-1 -- the
  // `lambda_o` of the TF24_floor cost curve, and an INPUT on the same footing as
  // `CF77_lambda_` above: same units, same meaning, never assigned by model code,
  // and not cleared by setup_clean_leaf() -- but, unlike it, with no default.
  //
  // WHY IT HAS NO DEFAULT, and this is the sharpest case in the class. Zero is a
  // defensible number and it is the WRONG default, because at `lambda_o = 0` this
  // curve IS TF24 -- the hydraulic half is TF24's cost, read at TF24's own traits
  // -- so a caller who seated "TF24_floor" and never set a price would silently be
  // running the production model under a new name. Every conductance-loss curve in
  // this package already has `lambda_o = 0` by construction; this one exists to
  // make the price explicit, so leaving it unset is refused rather than defaulted.
  //
  // ⚠️ IT IS THE CURVE'S ONLY PARAMETER. There is deliberately no `TF24_floor_a`: the
  // hydraulic half reads `TF24_cost_scale` and `TF24_beta2`, so TF24 is TF24_floor at
  // `lambda_o = 0` AT IDENTICAL PARAMETER VALUES and the comparison between them is
  // a one-restriction test rather than a fit of two differently-parameterised
  // models. An earlier version gave it its own scale; that made "is TF24 missing a
  // price of water" a three-restriction question for no gain.
  //
  // ⚠️ ITS SCALE IS SET BY THE LEAF, not chosen freely. lambda here runs 9e4 to 3e5
  // umol C (kg H2O)^-1 at this package's defaults, so a lambda_o far outside that
  // band pins the optimum against a bound and the answer describes the bracket
  // rather than the model. See `kLambdaCF77` in tests/cpp/test_golden.cpp, which
  // makes the same point about the same units.
  double TF24_floor_lambda_o = util::na_value;     // umol C (kg H2O)^-1
  // THE UNIT-COST RATIO OF LEAST-COST THEORY (@Prentice-2014's beta), and the
  // LeastCost curve's ONLY parameter. An INPUT on the same footing as
  // `CF77_lambda_`: never assigned by model code, not cleared by
  // setup_clean_leaf(), and defaulted to its ALIGNED value (see `CF77_lambda_`).
  //
  // Least-cost minimises `C_P = (c_E*E + c_V*Vcmax)/A`, so it has two unit costs
  // rather than one. ⚠️ ONE FIELD IS ENOUGH, AND THIS IS NOT A SIMPLIFICATION.
  // Scaling both costs by k multiplies C_P by k, which moves neither the argmax nor
  // the marginal cost of water: dividing through by `c_E` leaves
  // `C_P/c_E = (E + beta*Vcmax)/A` with `beta = c_V/c_E`, and lambda (below) has
  // `c_E` cancelling top and bottom. The ONLY quantity the missing factor reaches
  // is the reported `profit_`, which on this curve -- as on SOX and JW26 -- is a
  // ratio rather than a carbon flux and is not comparable with a subtracted curve's
  // profit anyway. So the curve is parameterised by the published quantity: beta.
  //
  // ⚠️ DIMENSIONLESS, WHICH FIXES THE UNITS OF BOTH TERMS IT WEIGHS. beta is a ratio
  // of two unit costs, so `E` and `Vcmax` must be in the SAME molar units for it to
  // be the published number. This package carries `E` in kg H2O m^-2 s^-1 and
  // `vcmax_` in umol CO2 m^-2 s^-1, so `leastcost_flux()` converts both to
  // mol m^-2 s^-1 before weighing them. Feeding this field a beta derived against
  // kg-and-umol would be off by a factor of 5.6e7.
  //
  // ⚠️ IT WEIGHS `vcmax_`, THE TEMPERATURE-ADJUSTED CAPACITY, not `vcmax_25`. That
  // is what appears in the assimilation the ratio divides by, and it is what makes
  // the first-order condition here the one @Prentice-2014 write down. It also means
  // the two costs move apart as the leaf warms, which is the acclimation channel
  // S1.4 of the theory paper flags: their beta carries a viscosity term (eta*) that
  // this package does not model, so beta here is beta/eta* at the leaf's own
  // temperature. At 25 C eta* is 1 by definition and the two agree.
  //
  // @Prentice-2014 report beta = 146 and @Dong-2017 fit 146 +- a factor of about
  // two; below ~10 the water term dominates and the leaf closes, above ~1e4 the
  // capacity term does and it opens to the wet bound.
  //
  // ⚠️ THE DEFAULT IS NOT THE PUBLISHED 146. At this package's default leaf, 146
  // pins the stem optimum at psi_crit from the reference psi_soil of 1 MPa
  // drier (the leaf opens fully at the reference), so it cannot be compared with
  // the other curves there. The aligned value, 146/2.83, sits below Dong's
  // factor-of-two band about 146.
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

  // THE LEAF'S OWN TEMPERATURE AT THE OPERATING POINT, deg C -- an OUTPUT, not a
  // driver, and the distinction is the whole reason this member exists.
  //
  // ⚠️ `leaf_temp_` is not it on the energy-balance path. There `set_physiology`
  // reinterprets the `leaf_temp` driver as AIR temperature (`Tair_ = leaf_temp_`)
  // and the leaf's temperature is solved per operating point from its own
  // transpiration. That value was computed, used to re-derive the whole Farquhar
  // temperature block, and then thrown away -- so the one quantity the PM path
  // exists to produce was the one thing a caller could not read. It reaches 62 C
  // at plant's drivers, which is not a detail a consumer should have to infer
  // from an assimilation it cannot explain.
  //
  // Off the PM path this is a copy of `leaf_temp_`, deliberately rather than NA:
  // an output column that is sometimes the driver and sometimes NA cannot be
  // plotted against anything.
  //
  // ⚠️ IT IS THE LEAF'S TEMPERATURE AT THE OPERATING POINT, which is NOT always
  // the temperature the reported fluxes were derived at -- and the difference is
  // #105, not a looseness in this definition. At the two shut-down exits the
  // Farquhar block still holds its Tair baseline, so the reported `R_d_` there is
  // 36% low against this temperature. Reporting Tair instead would make the pair
  // agree by reporting a temperature the leaf is not at, which is what kept that
  // inconsistency invisible.
  double Tleaf_;

  // Penman-Monteith leaf energy balance state (#523), only meaningful on the
  // use_energy_balance_ path. Set once per set_physiology; Tleaf itself is a
  // per-operating-point quantity computed in set_leaf_states_rates_from_psi_stem.
  //
  // ⚠️ Tair_, Rn_, g_rad_ and ra_ are all RE-SEATED by every set_physiology, so a
  // caller supplying its own radiation budget overrides them after each call.
  double Tair_;  // air temperature, deg C (reinterprets the leaf_temp driver)
  // ISOTHERMAL net radiation, W m^-2: absorbed shortwave minus the net longwave a
  // leaf AT AIR TEMPERATURE would exchange. The leaf's own departure from Tair is
  // carried by g_rad_ below rather than by this term, which is what keeps the
  // balance explicit -- see the longwave block in constants.hpp.
  double Rn_;
  // Radiative conductance, W m^-2 K^-1: 4*eps_leaf*sigma*Tair_K^3, the linearised
  // sensitivity of the leaf's longwave emission to its own temperature. Adds to
  // the sensible-heat conductance vol_heat_cap_air/ra_ in the balance, so a
  // radiatively coupled leaf sits CLOSER to air temperature than ra_ alone
  // implies: 13-18% of the total at the default leaf size and wind.
  double g_rad_;
  // AERODYNAMIC (boundary-layer) resistance to sensible heat, s m^-1. Convective
  // ONLY: the radiative path is g_rad_, added beside it in leaf_temp_from_E.
  // ⚠️ A caller injecting a COMBINED convective-plus-radiative resistance here
  // (Jones et al.'s r_HR, or 1/(Gbh + Grad)) must also set g_rad_ = 0, or the
  // radiative conductance is counted twice and every Tleaf - Tair and dTleaf/dE
  // comes out 5-7% small with nothing to flag it. g_rad_ = 0 recovers the
  // single-resistance balance Tair + (Rn_ - lambda*E)*ra_/vol_heat_cap_air.
  double ra_;
  // Gate for the PM leaf energy balance. Default OFF: today's path runs
  // (prescribed leaf_temp, single-shot cached Arrhenius). R-settable (#523 full
  // cut) so TF24 (via pars.use_energy_balance) and the leaf-level demo can
  // turn PM on; default preserves backward compatibility.
  bool use_energy_balance_ = false;
  // --- the residual (cuticular) conductance (#182) ----------------------------
  //
  // A leaf keeps losing water through its cuticle after its stomata have shut,
  // and the cuticle lets CO2 through as well. The residual flux is not selected:
  // the leaf optimises its STOMATAL flux and the xylem carries the total,
  //
  //     E   = E_s + r,            r  = g_min(T_leaf) D_leaf / P_atm
  //     g_c = g_cs + rho g_min,   g_cs = (E - r) P_atm / (1.67 D_leaf)
  //
  // with g_c the total conductance to CO2 the supply == demand residual reads,
  // and rho the cuticle's CO2:H2O permeability ratio (g_min_co2_ratio_). The
  // stem potential (or the collar) fixes E through the supply relation as it
  // always has, so only the conductance changes: g_cs is formed from E - r
  // instead of E, and every cost that reads the flux (CF77, TF24_floor,
  // LeastCost, ProfitMax's normalisers) reads the TOTAL. The first-order
  // condition keeps its form, which is the manuscript's eq-g0shift.
  //
  // ⚠️ THE CO2 PATHWAY IS WHAT REMOVES THE NO-FLOW FAMILY (#183). With g_c >=
  // rho g_min > 0 the residual A(c_i) = g_c (c_a - c_i) has a root for every
  // leaf: in the dark c_i sits above c_a and respiratory CO2 leaves through the
  // cuticle, and a leaf too hot to fix carbon does the same. So there is no
  // `min(0, A(c_a))`, no c_i assigned at c_a, no shade-death kind and no
  // hand-written shut state: a shut leaf is the ordinary evaluator at E_s = 0.
  // At g_min = 0 exactly the conductance of a shut leaf is zero and the
  // residual degenerates to A(c_i) = 0; solve_ci takes that limit.
  //
  // A TRAIT, default 0, so that every existing caller and all four golden
  // files are untouched: at g_min = 0 the term is an exact zero appended to
  // every expression it enters (hazard 13 names the cost of a sixteenth slot).
  // An in-class initialiser rather than a constructor argument, for R_d_25's
  // reason: plant's generated glue pins the constructor by arity.
  //
  // ⚠️ MEASURED CUTICULAR CONDUCTANCE, NOT A FITTED g0. A g0 fitted to gas
  // exchange is an order of magnitude larger and is mostly stomata that have
  // not finished closing (Duursma et al. 2018); the Medlyn route's `g0` below is
  // that quantity and is left alone. De Kauwe et al. (2022) put measured
  // cuticular conductance at 0.25-0.8 mmol m^-2 s^-1.
  double g_min = 0.0;  // mol H2O m^-2 s^-1, at 25 C
  // The cuticle's CO2:H2O conductance ratio. Boyer (2015, J. Exp. Bot.,
  // doi:10.1093/jxb/erv065) measured the cuticle's conductance to CO2 at 1/20
  // to 1/40 of its conductance to water vapour in grape and sunflower, in line
  // with Boyer et al. (1997); 0.05 is the upper end. Stomata are at 1/1.67.
  // Held constant: the temperature response below scales both pathways.
  double g_min_co2_ratio_ = 0.05;
  // The temperature response of g_min, ONE smooth expression and never a
  // breakpoint (#182): a low-temperature exponential, with a second exponential
  // blended in above a transition temperature by a softplus of stated width,
  //
  //     g_min(T) = g_min * exp( k_lo (T - 25) + w k_hi [ sp((T - T_t)/w) -
  //                                                     sp((25 - T_t)/w) ] )
  //
  // with sp(x) = log(1 + exp(x)). Below T_t the slope is k_lo, above it k_lo +
  // k_hi, and the two branches meet over about w degrees. The normalisation
  // makes the factor exactly 1 at 25 C, so the reference value is the trait.
  //
  // WHY THIS SHAPE. Cuticular permeance rises with temperature, steeply above a
  // phase transition of the cuticular waxes near 35-45 C (Riederer & Schreiber
  // 2001; Schuster et al. 2016, both reviewed in Duursma et al. 2018), and that
  // is what carries a hot leaf to hydraulic failure. The measured response is
  // usually drawn as two slopes meeting at a point, which is a kink in exactly
  // the regime where this package is removing kinks; the softplus is the C-inf
  // form of the same two slopes.
  //
  // ⚠️ THE DEFAULTS ARE PLACEHOLDERS, not fitted values: a Q10 near 1.35 below
  // the transition and a tenfold rise per 15 K above it, at a transition of
  // 40 C. They are inert at 25 C and at g_min = 0, which is where every golden
  // number and every manuscript claim is taken; a run that leans on the hot
  // branch needs measured values (registry: placeholder).
  double g_min_k_low_ = 0.03;     // K^-1, the slope below the transition
  double g_min_k_high_ = 0.15;    // K^-1, the EXTRA slope above it
  double g_min_T_trans_ = 40.0;   // deg C, the transition temperature
  double g_min_T_width_ = 2.5;    // K, the width over which the slopes blend
  // The soft cap that limits the realised residual flux to what the xylem can
  // carry: r* = r (1 + (r/E_cap)^p)^(-1/p), with E_cap = phi E_inf and E_inf =
  // kmax (G(inf) - G(psi_up)) the most the stem can ever deliver from its
  // upstream potential. Inert while r is well inside the cap ((r/E_cap)^p / p:
  // at half the cap 5e-4 of r), and a smooth minimum where the soil cannot
  // supply the cuticle -- the one regime a steady-state leaf has no answer for
  // (#183 leaves it to a later step; this form is a placeholder for it). Never
  // a min(): that would put a kink into psi_soil at the crossing.
  //
  // ⚠️ phi < 1 IS WHAT KEEPS THE SHUT STEM FINITE. At the cap the stem carries
  // phi of the flux it could carry at infinite tension, so its potential is
  // the one at which a fraction 1 - phi of the remaining conductance survives;
  // at phi = 1 that potential is infinite and the inverse table has no answer
  // (measured: a 7 MPa soil with 0.5 mmol of cuticle stopped the solve there).
  static constexpr double residual_cap_power = 8.0;
  static constexpr double residual_cap_fraction = 0.5;
  // OUTPUTS, written by every path out of a solve (hazard 8): the residual
  // conductance at the operating point's leaf temperature, and the residual
  // flux the leaf realised there, kg H2O m^-2 s^-1. `transpiration_` is the
  // TOTAL; the stomatal flux is `transpiration_ - residual_flux_`.
  double g_min_leaf_ = 0.0;        // mol H2O m^-2 s^-1
  double residual_flux_ = 0.0;     // kg H2O m^-2 s^-1
  // The second cap, on the collar route only: the realised residual flux is
  // also held below residual_cap_fraction of the most the ROOTS can deliver from
  // this soil (MultiLayerRoots::uptake_limit, kg H2O m^-2 s^-1), so the wet
  // bound -- the collar at which uptake meets the residual flux -- always
  // exists inside the uptake's range and moves continuously with the soil
  // (#184). Without it, as a dry layer approached the root integral's cap the
  // uptake's limit fell to zero under a residual flux that only the stem
  // bounded, and the collar ran off to where the stem's ceiling underflows: a
  // jump of ~10 MPa in psi_stem across 0.005 MPa of soil. Infinite (inert,
  // bit-exactly) off the collar route and at g_min = 0; seated per solve by
  // prepare_collar_solve and cleared by stem_wet_bound, because the stem routes
  // have no roots upstream of psi_soil.
  double residual_root_ceiling_ = std::numeric_limits<double>::infinity();
  // g_min_factor's one-entry memo: its five inputs, then its value. A cache
  // of a pure function, so it is mutable and never part of the leaf's state.
  mutable std::array<double, 6> g_min_factor_memo_{
      std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0, 0.0, 0.0};

  // The temperature factor of the residual conductance, 1 at 25 C exactly, and
  // its slope in temperature.
  // Served from the memo when nothing it reads has moved: with the energy
  // balance off every evaluation in a solve asks at the driver temperature,
  // and the exp and two softpluses were a third of what a residual conductance
  // cost the solve. Keyed on all five inputs, so a parameter changed between
  // solves cannot be answered with the old curve; a NaN key never matches.
  // The check is here, in the class, so it inlines (1.3% at g_min > 0).
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
  // At g_min = 0 the factor is never formed: set_leaf_vpd seats this at every
  // leaf temperature a solve establishes, and the exp it costs is waste there.
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
  // The cap at the roots' ceiling, applied to the stem-capped flux r1; the
  // identity, exactly, while that ceiling is infinite. The ceiling is not the
  // roots' limit U itself but C = U (U / (U + r))^2, with r the cuticle's
  // demand. C ~ U wherever the roots can deliver far more than the cuticle
  // asks, which is every ordinary soil; as U -> 0, C ~ U^3 / r^2, so the
  // realised flux falls much faster than the supply. Two measured reasons:
  //
  //   - With C = U the bound solved uptake = phi U as the last layer stopped
  //     conducting, and converged to psi + z / (1 - phi): a step of the gravity
  //     head (4.9e-3 MPa) against a dead soil, whose bound is the balance psi + z.
  //     Any C shrinking faster than U closes it.
  //   - U itself falls to zero LINEARLY at the root integral's cap (the root
  //     table holds G(inf) past it, #85), and with C = U^2 / (U + r) the bound's
  //     offset from the balance followed U, leaving a slope change there: the
  //     stem's slope in psi_soil went -2.3 on one side and 1 on the other. With
  //     the cube the offset follows U^2, and both one-sided slopes converge to
  //     the same value (1.0000 for the stem, -0.0866 for profit).
  //
  // So the shut leaf is continuous with a continuous slope through the cap.
  // Its second derivative is not, and that, like the cap itself, goes with a
  // closed-form root tail (#85, #184 section 2).
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
  // The residual flux the leaf realises at an upstream potential and a TOTAL
  // flux, kg H2O m^-2 s^-1. Pure: the leaf temperature is taken from the total
  // flux on the energy-balance path and from the driver off it, so a bracket
  // root-find can ask about a collar without seating anything. Exactly zero at
  // g_min = 0.
  //
  // ⚠️ THE g_min = 0 TEST IS IN THE CLASS ON PURPOSE, as for the two slopes
  // below. The full functions are too large to inline, and these are called on
  // every evaluation of the solve: a call returning zero cost ~2.5% of an open
  // leaf's solve at g_min = 0 (profiled), which the default must not pay.
  double residual_flux(double psi_upstream, double E_total) const {
    return g_min == 0.0 ? 0.0 : residual_flux_active(psi_upstream, E_total);
  }
  double residual_flux_active(double psi_upstream, double E_total) const;
  // The demand before the cap, kg H2O m^-2 s^-1 at a leaf temperature, and its
  // slope in that temperature (through g_min(T) and the deficit).
  double residual_demand_at(double leaf_temp) const;
  double residual_demand_dT(double leaf_temp) const;
  // The realised flux's slope in the upstream potential, through the ceiling,
  // at a leaf temperature; and the slope in leaf temperature of the TOTAL
  // conductance to CO2 that the residual contributes, -gc_const dr*/dT + rho
  // dg_min/dT, for the energy-balance term. Both exact zeros at g_min = 0.
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
  // The same quantities on the tape: G(inf) with its trait rows, and the
  // realised residual flux at a collar, carrying rows in g_min, kmax and the
  // stem curve's pair. The leaf temperature and deficit are the seated
  // operating point's and are passive, as the deficit is in collar_coords_at.
  template <class S> S stem_integral_limit_at(const leaf_pars<S>& pars) const;
  // Declared ahead of its definition with the derivative surface below; the
  // draw carries the roots' ceiling on the residual flux.
  template <class S> struct SupplyDraw;
  template <class S>
  S residual_flux_at(const S& collar, const leaf_pars<S>& pars,
                     const SupplyDraw<S>& draw) const;
  // The cuticle's conductance to CO2 at the seated leaf temperature.
  double cuticular_cond_CO2() const { return g_min_co2_ratio_ * g_min_leaf_; }
  // Net assimilation as c_i -> infinity: the colimited capacity itself, less
  // respiration. Where a leaf with no conductance cannot net zero at any c_i
  // this is its carbon, and c_i is reported as +inf.
  double assim_saturated() const;
  // c_i from a total conductance to CO2: the supply == demand root, bracketed
  // above c_a where A(c_a) < 0. See the definition for the g_c = 0 limit.
  double solve_ci(double gc);
  // Boundary-layer inputs for ra = C_ra*sqrt(d/U0) (doc 4.1). d is a per-strategy
  // trait (set from pars.d in prepare_strategy); wind_speed_ is the per-timestep
  // above-canopy driver (set from the environment before set_physiology). Both
  // R-settable so a bare Leaf can exercise the wind model; if either is unusable
  // set_physiology falls back to the fixed ra. Only read on the PM path.
  double d_ = 0.05;          // characteristic leaf dimension, m
  double wind_speed_ = 2.0;  // above-canopy wind speed U0, m s^-1

  // Leaf absorptances, the fractions of INTERCEPTED radiation the leaf absorbs.
  // Light reaches the leaf as intercepted quantities (PPFD_ and NIR_int_); these
  // are the leaf's own optical properties, so the conversion to absorbed happens
  // here and nowhere upstream. alpha_PAR_ reaches both photosynthesis (APAR =
  // alpha_PAR_ * PPFD_, which phi_J multiplies) and the energy balance;
  // alpha_NIR_ only the energy balance. Leaf properties, not fitted traits.
  //
  // ⚠️ IN A GAS-EXCHANGE FIT ONLY phi_J * alpha_PAR_ IS IDENTIFIED: the two enter
  // J as a product. Fix alpha_PAR_ (at this default or a measured absorptance)
  // and the fitted phi_J is per absorbed photon.
  double alpha_PAR_ = 0.85;  // PAR absorptance, unitless
  double alpha_NIR_ = 0.20;  // near-infrared absorptance, unitless
  // Intercepted near-infrared (700-3000 nm) per unit leaf area, W m^-2: the
  // non-PAR part of the shortwave, which heats the leaf but drives no
  // photosynthesis. A per-call driver like wind_speed_, set before
  // set_physiology, and read on the energy-balance path ONLY.
  //
  // ⚠️ NO DEFAULT, BECAUSE NONE IS RIGHT. It depends on the light source, which
  // the leaf cannot see: about PPFD_ / umol_par_per_joule in daylight, zero under
  // an LED. A default of either would be silently wrong for the other, so
  // set_physiology refuses a non-finite value on the energy-balance path.
  double NIR_int_ = util::na_value;

  // --- Sicangco et al. (2026) thermal cost, default OFF -----------------------
  // An instantaneous cost of PSII damage, sigmoid in LEAF temperature:
  //
  //   TC(Tleaf) = 1 / (1 + exp(-r*(Tleaf - T50))),   r = 2/(T50 - Tcrit)
  //
  // and, in the same paper's ProfitMaxTC, Jmax is scaled by (1 - TC) so electron
  // transport is destroyed at the same rate the cost is incurred. Both halves are
  // behind this one gate.
  //
  // ⚠️ IT IS INSTANTANEOUS, and that is a modelling position rather than an
  // omission. plant's ATLS layer models thermal damage as a LASTING ratchet --
  // damage accumulated at one time step is still there at the next. This is the
  // other thing: a cost paid at the temperature the leaf is at right now, with no
  // memory. They are not interchangeable and neither is a port of the other.
  //
  // ⚠️ AND THE PAPER AND ITS CODE DISAGREE ON r. The published Eqn 8 area reads
  // "r = (T50 - Tcrit)/2"; gsthermal's F0_func and TJmax_updated both compute
  // `r = 2/(T50 - Tcrit)`. The code is what produced the figures, so the code is
  // what is implemented here.
  bool use_thermal_cost_ = false;
  double T50_ = 50.4;    // deg C; Sicangco Table 2, heatwave treatment
  double Tcrit_ = 46.5;  // deg C; Sicangco Table 2, heatwave treatment

  // --- ATLS thermal damage and acclimation, default OFF (traitecoevo/plant#566) --
  // A LASTING damage model on three states, fractions of leaf area named for what
  // they do rather than for a mechanism:
  //
  //   functional F  <=>  impaired I  -->  damaged D,      F = 1 - I - D
  //
  //   dA/dt = acclim_alpha_ * softplus(T - acclim_T_; acclim_s_) - acclim_beta_ * A
  //   dI/dt = k_imp(T) * F - (k_rec(T) + k_mat(T)) * I
  //   dD/dt = k_mat(T) * I
  //   k(T)  = k * exp[(Ea/R)(1/T_ref - 1/T)]   for each of the three rates
  //
  // with T the LEAF temperature throughout and every rate per DAY. Every rate has
  // the one form, a centred Arrhenius, and nothing switches on or off: repair and
  // maturation are stated at 25 C, and impairment at its isokinetic temperature
  // k_i_Tiso_ moved by the acclimation shift. Impaired area returns to functional
  // only by repair; damaged area does not return. damage_rates() evaluates the
  // right-hand side and damage_trajectory() integrates it over a
  // leaf-temperature series, so a leaf experiment can be
  // fitted with this class alone; a caller that integrates the states itself (plant's
  // TF24t, which adds dilution by leaf turnover and growth to dD/dt) calls
  // damage_rates() and seats A_, I_ and D_ before each set_physiology.
  //
  // Behind the gate, update_temperature_dependent_params also:
  //
  //   * shifts the entropy term of the Vcmax and Jmax peaked Arrhenius so their
  //     optimum moves by topt_offset_ + dTopt_max_ * A/(A + K_A_) degrees;
  //   * scales both capacities by the functional fraction F = 1 - I_ - D_, held
  //     above zero by a soft floor of width damage_surv_floor: impaired and
  //     damaged area does not photosynthesise;
  //   * reports k_imp_eff_, k_rec_eff_ and functional_ at that temperature;
  //   * adds four maintenance costs to R_d_, after the Q10 response.
  //
  // ⚠️ IMPAIRMENT IS NOT THE PEAKED ARRHENIUS'S DEACTIVATION. The Vcmax and Jmax
  // curves are the INSTANTANEOUS, fully reversible temperature response, and
  // nothing about them persists past the call; impairment is rate-and-state and
  // lasts until repaired. They share no parameter, so recalibrating the
  // photosynthetic response cannot move the impairment rate. A hot leaf loses
  // capacity both ways, by different mechanisms on different timescales, and that
  // is intended rather than a double count.
  //
  // With the gate off none of these inputs is read by the solve, k_imp_eff_ and
  // k_rec_eff_ are written as zero and functional_ as one, and every other number
  // is bit-identical to a leaf without the layer.
  //
  // ⚠️ DAMAGE ENTERS THE BENEFIT, NOT THE HYDRAULIC COST: both capacities carry
  // F and R_d_ carries the costs, so the optimiser sees the states only through
  // assimilation. The rates enter no objective: a leaf cannot trade carbon now
  // against impairment later.
  //
  // ⚠️ R_d_ IS NOT SCALED BY F. Respiration is R_d_25 * Q10(T) whatever the
  // state, plus the costs below: impaired and damaged area costs carbon and
  // gains none.
  //
  // ⚠️ ONLY THE TRAIT GRADIENT CARRIES THIS LAYER. leaf_gradient() and the
  // batch difference the solve's own dprofit and outputs, which read the members
  // this block writes, and take A_, I_, D_ and the parameters that reach the
  // solve as inputs (gradient.hpp). The derivative surface (operating_collar,
  // billed_outputs_at) reads capacities off the parameter pack, whose kernels do not
  // apply the layer, and the closed form has not been checked with it, so those
  // refuse a leaf with the gate on. The dynamics' own sensitivities are
  // damage_trajectory_gradient().
  bool use_thermal_damage_ = false;
  // Constitutive thermostability shift, deg C. Moves BOTH the photosynthetic
  // optimum (through the entropy term) and the impairment's isokinetic
  // temperature by the same amount; the acclimated shifts below are separate.
  double topt_offset_ = 0.0;
  // Impairment, functional -> impaired, Bison & Michaletz (2026)'s first-order
  // PSII inactivation k = A exp(-Ea/RT) under their Meyer-Neldel relation
  // ln A = m Ea + b, which pins every species' rate to one value exp(b) at the
  // isokinetic temperature T_iso = 1/(m R). Written about that point,
  //
  //     k_imp = k_i_ * exp[(k_i_Ea_/R)(1/T_iso - 1/T)],
  //     T_iso = k_i_Tiso_ + topt_offset_ + dTimpair_max_ * A/(A + K_A_)  (deg C),
  //
  // k_i_ (day^-1) is that common rate and k_i_Ea_ (J mol^-1) the species'
  // steepness, tenfold per ln(10)*R*T^2/k_i_Ea_ degrees (6.1 K at the default near
  // 45 C). A different Ea pivots the curve about T_iso rather than moving it.
  //
  // The defaults are their 177-species fit: exp(b) = 3.2e-3 s^-1, T_iso = 50.5 C
  // (Extended Data Table 2) and the mean Ea, 316 kJ mol^-1 (Extended Data Table 1;
  // species span 116-619). ⚠️ THEIR A IS PER MINUTE IN THE PAPER AS FIRST
  // PUBLISHED, mislabelled per second (author correction, Zenodo
  // 10.5281/zenodo.23752859 v3 corrects it): exp(b) here is the per-second value
  // times 86400. Their Ea is kJ mol^-1, and these are J mol^-1.
  double k_i_ = 276.0;
  double k_i_Ea_ = 316.0e3;
  double k_i_Tiso_ = 50.5;
  // The no-extra-protection baseline k_i_ is costed against, day^-1: protection
  // bought is k_i_baseline_ - k_i_, signed, so a leaf less protected than the
  // baseline saves what a more protected one pays.
  double k_i_baseline_ = 276.0;
  // Repair, impaired -> functional, and maturation, impaired -> damaged: each a
  // rate at 25 C (day^-1) on a centred Arrhenius with its own activation energy
  // (J mol^-1). Neither shuts off when hot: impairment outruns the default repair
  // above about 32 C and keeps rising tenfold every ~6 K, so a hot leaf is
  // impaired whatever repair does.
  double k_rec_ = 0.25;
  double k_rec_Ea_ = 0.0;
  double k_mat_ = 0.02;
  double k_mat_Ea_ = 0.0;
  // Acclimation kinetics: dA/dt = acclim_alpha_ * softplus(T - acclim_T_;
  // acclim_s_) - acclim_beta_ * A, with T the LEAF temperature. acclim_alpha_
  // (day^-1 deg C^-1) is the gain, acclim_beta_ (day^-1) the decay, so 1/beta is
  // the timescale and (alpha/beta)*softplus the equilibrium; acclim_T_ (deg C) is
  // the forcing threshold and acclim_s_ (deg C^-1) the softplus sharpness.
  double acclim_alpha_ = 0.02;
  double acclim_beta_ = 0.01;
  double acclim_T_ = 30.0;
  double acclim_s_ = 1.0;
  // What acclimation buys: the most the photosynthetic optimum can rise, and the
  // most the impairment's isokinetic temperature can rise, deg C, both
  // saturating in A with half-saturation K_A_. Neither is a hard cap on anything but the shift.
  double dTopt_max_ = 6.0;
  double dTimpair_max_ = 6.0;
  double K_A_ = 1.0;
  // STATE, set by the caller before each set_physiology: the acclimation state,
  // and the impaired and damaged fractions of leaf area. Zero is an unacclimated,
  // intact leaf. set_physiology refuses I_ < 0, D_ < 0 and A_ <= -K_A_ (the
  // acclimation shift's pole); a hair below zero in A_ is read as itself.
  double A_ = 0.0;
  double I_ = 0.0;
  double D_ = 0.0;
  // Maintenance costs added to R_d_, umol C m^-2 s^-1 per unit of what each
  // multiplies: the acclimation state A_; the repair capacity k_rec_ (day^-1);
  // the protection bought, k_i_baseline_ - k_i_ (day^-1); and the realised
  // repair flux k_rec_eff * I_ (day^-1). Zero by default, so the layer costs
  // nothing until a caller prices it.
  double c_acclim_maint_ = 0.0;
  double c_repair_maint_ = 0.0;
  double c_protect_maint_ = 0.0;
  double c_repair_flux_ = 0.0;
  // The width of the functional fraction's soft floor. Strictly positive so the
  // capacities stay positive and the solve stays well-conditioned as
  // I_ + D_ -> 1.
  static constexpr double damage_surv_floor = 1e-6;
  // OUTPUTS, at the temperature the block was last seated at -- the operating
  // point's leaf temperature after a solve on the energy-balance path, the
  // leaf_temp driver off it.
  double k_imp_eff_ = 0.0;   // impairment rate there, day^-1; zero with the gate off
  double k_rec_eff_ = 0.0;   // repair rate there, day^-1; zero with the gate off
  double functional_ = 1.0;  // the F that scaled the capacities; one with the gate off

  // Cells a stem route scans to locate the basin before refining, WHERE ONE IS
  // ASKED FOR -- see `basin_scan_cells()` below, which returns 0 with the energy
  // balance off. This is only the size of the scan, not whether there is one.
  //
  // ⚠️ THE REFINEMENT TOLERANCE IS NOT `GSS_tol_abs` AND IS NOT SETTABLE. The
  // winning cell is refined by a ROOT-FIND on dJ/dpsi == 0 at `collar_root_tol`;
  // the `(cell width) * 1e-4` figure belongs to
  // maximise_over_closed_interval_foc's golden-section FALLBACK, which fires only
  // where neither cell end has a usable derivative. `GSS_tol_abs` reaches neither,
  // so tightening it to sharpen a stem argmax is inert -- a downstream project has
  // already tried exactly that.
  //
  // Costs n+1 objective evaluations per solve where it runs, so it is not free,
  // and 64 is where the answer stops moving: measured over a 1728-row driver sweep
  // against a 2001-point reference, 64 matches it on every row of both objectives
  // while 32 leaves 6 rows short by 3.9e-04. A member rather than a constructor
  // argument because that arity is pinned by plant's generated glue and by the CI
  // consumer program.
  int boundary_scan_n_ = 64;

  // Cells to scan before refining. ZERO unless the energy balance is on, and that
  // is measured rather than assumed: over a 1728-row sweep (8 air temperatures x 4
  // gate combinations x 6 soil potentials x 3 deficits x 3 light levels) the only
  // stem objectives carrying TWO prominent interior basins are TF24 (21 rows) and
  // JS22 (66), and EVERY one of those rows has `use_energy_balance_` true. With the
  // gate off the objective is unimodal, so endpoints plus a root-find on the
  // first-order condition is not an approximation -- it is the whole answer, and it
  // costs `n` fewer objective evaluations on the path plant calls millions of times.
  //
  // ⚠️ Re-measure this before widening the gate predicate. The valleys involved are
  // 7.0e-01 (TF24) and 3.5e+00 (JS22) deep, so a missed basin is a large error, not
  // a rounding one -- this is a cheap scan to skip and an expensive one to skip
  // wrongly.
  int basin_scan_cells() const {
    return use_energy_balance_ ? boundary_scan_n_ : 0;
  }
  // INTERCEPTED photosynthetic photon flux density per unit leaf area, umol m^-2
  // s^-1: incident on the leaf, NOT absorbed. Absorption is the leaf's own
  // business (alpha_PAR_), so the caller hands over the light that arrives and
  // never the light that is kept. Seated by set_physiology's PPFD argument.
  double PPFD_;
  double atm_vpd_;
  // The vapour pressure deficit the DIFFUSION equations use, kPa. Off the
  // energy-balance path it is exactly `atm_vpd_`; on it, the leaf-to-air deficit
  // at the operating-point leaf temperature (#7). See set_leaf_vpd.
  double vpd_leaf_;
  double atm_o2_kpa_;

  // --- Temperature-response parameters ---------------------------------------
  // Settable, NOT constexpr, and that is a deliberate correction. These are
  // exactly plantecophys's EaV / EdVC / delsC / EaJ / EdVJ / delsJ, which are
  // *user parameters* there -- and whose defaults that package revised at v1.4
  // after a literature review. They are also precisely what thermal acclimation
  // modifies (Kattge & Knorr), so TF24t needs them mutable. `constexpr` asserted
  // that they were not modelling choices, which was wrong.
  //
  // Defaults are the constants in leaf/constants.hpp, so there is still one source
  // of truth for the published values.
  double vcmax_ha_ = phylloptim::vcmax_ha;    // activation energy, J mol^-1
  double vcmax_H_d_ = phylloptim::vcmax_H_d;  // deactivation energy, J mol^-1
  double vcmax_d_S_ = phylloptim::vcmax_d_S;  // entropy term, J mol^-1 K^-1
  double jmax_ha_ = phylloptim::jmax_ha;
  double jmax_H_d_ = phylloptim::jmax_H_d;
  double jmax_d_S_ = phylloptim::jmax_d_S;
  // Bernacchi kinetics: reference value at 25 C and activation energy. Weaker case
  // for being settable than the six above -- enzyme kinetics vary less among
  // species -- but they are literature values that get revised, and bigleaf
  // exposes them.
  double gamma_25_ = phylloptim::gamma_25;    // CO2 compensation point, umol mol^-1
  double gamma_ha_ = phylloptim::gamma_ha;
  double kc_25_ = phylloptim::kc_25;          // Rubisco Km for CO2, umol mol^-1
  double kc_ha_ = phylloptim::kc_ha;
  double ko_25_ = phylloptim::ko_25;          // Rubisco Km for O2, umol mol^-1
  double ko_ha_ = phylloptim::ko_ha;
  // Dark respiration at the 25 C reference, umol m^-2 s^-1. Many datasets report
  // it directly (Sabot et al. give `Rlref` per site), and across their 16 species
  // the implied fraction of vcmax_25 spans 0.0046 to 0.0302, so it is a trait in
  // its own right rather than a fixed multiple of anything.
  //
  // The default is 0.015 * 96, the value the old vcmax-derived form gave at the
  // default vcmax_25. Initialised here rather than in the constructors' init lists
  // because plant's RcppR6 bindings pin the 17-argument constructor by arity, so
  // this cannot become an 18th argument without breaking plant's generated glue.
  double R_d_25 = 1.44;
  // Joshi & Stocker (2022)'s hydraulic unit cost, `gamma` in
  // `F = A - alpha*Jmax - gamma*(dpsi)^2`. A TRAIT, initialised here rather than
  // in the constructors' init lists for exactly the reason R_d_25 is: RcppR6 pins
  // the constructor by arity.
  //
  // ⚠️ THE CAPACITY TERM IS NOT IMPLEMENTED, so this is Joshi's HYDRAULIC cost
  // and not their model -- theirs optimises Jmax jointly, and dropping that term
  // is a real omission rather than a simplification. Do not describe a run of this
  // curve as "Joshi et al. (2022)".
  //
  // ⚠️ THE DEFAULT IS NOT A LITERATURE VALUE. It is ALIGNED to TF24 at the
  // reference environment (see `CF77_lambda_`), and no single value matches more
  // than one soil potential: the gamma reproducing TF24's cost at its own
  // optimum runs 0.287 (psi_soil 0.5 MPa) to 1.761 (3.0 MPa), a 6.1x range,
  // because TF24's cost tracks the ABSOLUTE potential through (1 - f) while this
  // one tracks the DROP -- and across a drydown those move in opposite directions.
  // That spread is the models' structural difference, not a fitting problem.
  double JS22_gamma = 0.4717887793;    // umol C m^-2 s^-1 MPa^-2
  // Wolf, Anderegg & Pacala (2016) carbon maximisation, in the form Anderegg et al.
  // (2018) gave it and Sabot's `TractLSM` implements (`SPAC/fregulate.py`,
  // `dcost_dpsi = Alpha*|P| + Beta`): a marginal cost LINEAR in the absolute leaf
  // water potential,
  //
  //     dC/dpsi = CMax_a*psi + CMax_b
  //
  // Wolf et al. argued the cost should be concave-up without specifying it; the
  // linear marginal form is Anderegg's. Read the form off TractLSM rather than the
  // 2016 paper.
  //
  // ⚠️ CMax_b IS SIGNED, and negative in the paper's own convention. TractLSM's
  // docstring says it stores Beta with the sign inverted "so as to get a positive
  // parameter value" and that a reported value should carry the minus sign, so the
  // sign here has to be pinned by reproducing a TractLSM number and not by reading
  // the formula. It is therefore NOT range-checked the way JS22_gamma is.
  //
  // ⚠️ NEITHER DEFAULT IS A LITERATURE VALUE. CMax_a is ALIGNED to TF24 at the
  // reference environment with CMax_b = 0 (see `CF77_lambda_`). Measured: the
  // CMax_a reproducing TF24's cost at TF24's own optimum runs 0.406 (psi_soil
  // 0.5 MPa) to 0.807 (3.0 MPa) -- a 2.0x range, against JS22's 6.1x, which is the
  // absolute-versus-drop distinction showing up as a number.
  double CMax_a = 0.6385596982;        // umol C m^-2 s^-1 MPa^-2
  double CMax_b = 0.0;                 // umol C m^-2 s^-1 MPa^-1
  // Dark respiration's temperature response, Tjoelker et al. (2001):
  //
  //     R_d(T) = R_d_25 * Q10(T)^((T - 25) / 10)
  //     Q10(T) = rd_q10_intercept_ - rd_q10_slope_ * (T + 25) / 2
  //
  // The Q10 is evaluated at the MEAN of the measurement and reference
  // temperatures, which is the form the land-surface literature implements, and it
  // DECLINES with temperature: a constant Q10 of 2 is too aggressive at the top of
  // the range, putting R_d 2.8x its 25 C value over 15 K and leaving the solve no
  // operating point at all by 45 C.
  //
  // Set the slope to zero and the intercept IS a constant Q10, for anyone who
  // wants the conventional form.
  double rd_q10_intercept_ = 3.09;
  double rd_q10_slope_ = 0.0430;
  double atm_kpa_;
  // Conversion from a mixing ratio (umol mol^-1) to a partial pressure (Pa).
  // DERIVED from atm_kpa_, not a constant: it is 1e-6 * P, so the old
  // `umol_per_mol_to_Pa_ = 0.1013` was atmospheric pressure of 101.3 kPa in
  // disguise. With atm_kpa_ a settable input, hard-coding it made the model
  // internally inconsistent away from sea level -- the stomatal side responded to
  // atm_kpa while Gamma*, Kc, Ko, Km and the ci root-find bounds all kept assuming
  // 101.3 kPa. This is the trap plantecophys warns about in bold for its own
  // `Patm`. Bit-identical to the old literal at 101.3 kPa (checked).
  double umol_per_mol_to_Pa_;
  double ca_;
  double opt_root_psi_;
  double assim_max_;

  // --- Sperry (2017) ProfitMax outputs, all unitless -------------------------
  // Written ONLY by the ProfitMax path (optimise_psi_stem_single<ProfitMax> and
  // profit_psi_stem_ProfitMax). Separate members rather than reusing
  // `hydraulic_cost_`, which is in carbon units on the TF24 path and in
  // conductance units on the Sperry-cost one: three meanings behind one name is
  // hazard 8 waiting to happen.
  // Normalisers, seeded by prepare_profitmax() and read by the two functions
  // above. Cached rather than recomputed per candidate because they are constant
  // over one solve and each costs an exp+pow.
  double profitmax_A_max_;   // |A|max over the supply stream, umol m^-2 s^-1
  double profitmax_k_soil_;  // k(psi_soil), kg m^-2 s^-1 MPa^-1
  double profitmax_k_span_;  // k(psi_soil) - k_crit, the HC denominator
  // ⚠️ profitmax_A_max_ / profitmax_k_span_ is NOT a lambda, whatever it gets
  // called: its units are carbon per unit CONDUCTANCE lost, umol C MPa (kg H2O)^-1,
  // because it multiplies `k(psi_soil) - k(psi)`. It is ProfitMax's normaliser and
  // nothing more. The comparable quantity is `lambda_emergent_` below.
  // The scan prepare_profitmax() already runs, kept so the objective can be
  // rebuilt on it without evaluating the model a second time.
  // THE MARGINAL COST OF WATER THE OPERATING POINT IMPLIES, and the one output
  // every cost curve can report on the same axis:
  //
  //     lambda_emergent = (dC/dpsi) / (dE/dpsi)      umol C (kg H2O)^-1
  //
  // Same units as the `CF77_lambda_` INPUT, and directly comparable with it -- that is
  // the point. Cowan-Farquhar prices water at a constant, so its emergent lambda
  // IS that constant; every other curve has one that moves with the drivers, and
  // this is where you read it. Written by whichever optimiser ran; NA before one
  // has.
  double lambda_emergent_ = util::na_value;
  double carbon_gain_;         // CG = A/|A|max over the supply stream
  double hydraulic_cost_norm_; // HC = [k(psi_soil)-k(psi)]/[k(psi_soil)-k_crit]
  double thermal_cost_;        // TC, zero unless use_thermal_cost_

  
  double opt_psi_stem_;
  // Kept for plant's bindings (hazard 7) and set to NA by setup_clean_leaf;
  // nothing writes it, and the operating point's ci is ci_.
  double opt_ci_;
  // ⚠️ KEPT ONLY BECAUSE plant's RcppR6 bindings name it (hazard 7), and nothing
  // in this package writes or reads it. Initialised so a read through those
  // bindings is a defined NA rather than whatever the memory held; this package
  // no longer binds it.
  double count = util::na_value;
  double E_up_;

  // --- Medlyn stomatal-conductance model (from develop #450) ------------------
  // Standalone, R-callable alternative to the root-collar profit optimisation
  // (solve_medlyn_ci_*); NOT used by the TF24 compute path, which optimises
  // psi_stem directly. g0/g1 default to the published values and are exposed as
  // settable fields. beta_ (soil-moisture stress) uses theta_/theta_w_/theta_fc_,
  // set from the default soil-moisture members in set_physiology.
  double g0 = 0.022;        // residual stomatal conductance (umol m^-2 s^-1)
  double g1 = 2.57;         // sensitivity to vpd (kPa^0.5)

  // The H2O:CO2 stomatal diffusion ratio. Settable because it encodes a CONVENTION,
  // not a property of this leaf: 1.67 here, 1.6 in Medlyn et al. (2011) and the `g1`
  // literature built on it, from the binary diffusivities (#50).
  //
  // ⚠️ IT REACHES THE SOLVE, not just a reported output -- the conductance conversion
  // and `dprofit`'s `gc_const` -- so a solved leaf must be re-solved after changing
  // it. Nothing caches a conductance, so there is no cache to invalidate.
  double H2O_CO2_stom_diff_ratio_ = phylloptim::H2O_CO2_stom_diff_ratio;

  // The diffusion ratio the MEDLYN route uses, deliberately 1.6 rather than
  // `H2O_CO2_stom_diff_ratio_`. It appears in three places that must agree -- the
  // USO expression, the solver residual, and the conversion of the result to a CO2
  // basis -- so it is named once here rather than written as a literal three times.
  //
  // ⚠️ WHY NOT 1.67. Published `g1` values are fitted with 1.6, so the USO
  // expression has to use it or the parameter means something else. Given that, the
  // conversion back has to use the SAME constant: the residual below equates
  // `gs/1.6` to the CO2 flux, so dividing the result by 1.67 instead reports a
  // conductance inconsistent with the `ci` the solve returned -- measured as a
  // systematic -4.2% violation of `A = gs*(ca-ci)`, where 1.6 satisfies it. This is
  // a property of the published parameterisation, not of the leaf, which is why it
  // is a constant here while `H2O_CO2_stom_diff_ratio_` is a settable field.
  static constexpr double medlyn_H2O_CO2_ratio = 1.6;

  // ⚠️ mol CO2 m^-2 s^-1, and it is the CO2 basis that is load-bearing: the USO
  // expression is a WATER VAPOUR conductance, so `medlyn_model_gs` divides by
  // `medlyn_H2O_CO2_ratio` before storing. Omitting that division is #129, which
  // reported this field 1.6x high on both Medlyn entry points.
  double medlyn_model_gs_;  // mol CO2 m^-2 s^-1
  double theta_w_;          // current soil water content at wilting point (m^3 m^-3)
  double theta_fc_;         // current soil water content at field capacity (m^3 m^-3)
  double theta_;            // current soil water content (m^3 m^-3)

  // 1-entry memo for transpiration(). Within a single root-collar/profit solve
  // leaf_specific_conductance_max_ and the stem curve are fixed, so
  // supply-side transpiration depends only on (psi_stem, psi_upstream). That
  // pair is queried repeatedly with identical values per profit evaluation
  // (psi_stem_to_ci -> stom_cond_CO2 -> transpiration, then transpiration again
  // for transpiration_). Caching the last result avoids the redundant table
  // lookups; it is invalidated in set_physiology when conductance/soil change.
  bool   transpiration_cached_ = false;
  double transpiration_cache_psi_stem_ = 0.0;
  double transpiration_cache_psi_upstream_ = 0.0;
  double transpiration_cache_value_ = 0.0;

  // Cache for the temperature/O2-dependent photosynthesis parameters set in
  // set_physiology (vcmax_, jmax_, gamma_, ko_, kc_, R_d_, km_). They are pure
  // functions of (leaf_temp_, atm_o2_kpa_) and constants, so when the key is
  // unchanged the Arrhenius transcendentals are skipped and the members reused
  // (bit-identical: same inputs -> same outputs). In the current driver
  // leaf_temp_/atm_o2_kpa_ are constant across the run, so this fires once.
  // NOTE: electron_transport_ is deliberately NOT cached here -- it also depends
  // on the per-call PPFD_ and is recomputed every call.
  bool   photo_temp_cached_ = false;
  // The leaf temperature the block above currently holds.
  //
  // ⚠️ NOT `leaf_temp_`, AND THE PACK-READING KERNELS MUST READ THIS ONE. On the
  // energy-balance path `leaf_temp_` is AIR temperature and the block is re-seated
  // per candidate psi at the leaf's own -- measured 10 K apart at Tair 40 -- so a
  // pack kernel reading `leaf_temp_` there evaluated the Farquhar block at air
  // temperature while the model evaluated it at the leaf's: 100% of assimilation,
  // 84% of electron transport, 8.3% of dA/dci. The transpose identity cannot see
  // it, because it reads one surface in two directions and both directions are
  // wrong together; test_leaf's member-against-pack comparison is what does.
  double photo_temp_;
  // ⚠️ THE KEY MUST COVER EVERY SCALAR THE BLOCK READS, NOT JUST THE DRIVERS.
  // It used to be (leaf_temp_, atm_o2_kpa_) alone, and the argument for that was
  // "same inputs -> bit-identical outputs, so reusing is exact". That argument was
  // true while the temperature-response parameters were unreachable C++ members
  // and became FALSE the moment they were bound to R: `l$rd_q10_slope_ <- 0`
  // followed by `set_drivers()` at the same temperature took a cache HIT and
  // silently kept the old response, with A unchanged to every digit.
  //
  // So the key is now every input of update_temperature_dependent_params(). Two
  // consequences worth knowing:
  //
  //   * it also covers `vcmax_25` and `jmax_25`, which closes the third and least
  //     visible half of hazard 10 -- a bare `l$vcmax_25 <- x` write no longer
  //     leaves `vcmax_`/`jmax_`/`R_d_` describing the old value. `set_traits()` is
  //     still the right way to change a trait (the vulnerability tables and the
  //     solved point need clearing too), but the silent-wrong-number failure mode
  //     is gone.
  //   * the cost is 17 double comparisons per set_physiology() call, i.e. per
  //     driver set, NOT per inner solve iteration. Measured: within run-to-run
  //     noise on bench_solve.
  static constexpr int photo_temp_key_size = 23;
  std::array<double, photo_temp_key_size> photo_temp_cache_key_{};
  std::array<double, photo_temp_key_size> photo_temp_key() const;
  // ⚠️ KEPT ONLY BECAUSE plant's RcppR6 bindings name it (hazard 7), and nothing
  // in this package writes or reads it. Initialised so a read through those
  // bindings is a defined NA rather than whatever the memory held; this package
  // no longer binds it.
  std::vector<double> f_r;
  // TODO: move into environment?

  // TODO: atm_vpd - now set in set_physiology although ideally should be moved to enviroment
  double atm_vpd = 2.0; //kPa
  double ca = 40.0; // Pa
  double atm_kpa = 101.3; //kPa
  //partial pressure o2 (kPa)
  double atm_o2_kpa = 21;
  //leaf temperature (deg C)
  double leaf_temp = 25;
  // default soil-moisture content used by the Medlyn beta_ stress factor
  // (matches the fixed values develop's TF24 caller passed to set_physiology)
  double theta_w = 0.2;  //m^3 m^-3
  double theta_fc = 0.5; //m^3 m^-3
  double theta = 0.3;    //m^3 m^-3
  // density of water


  
  // set-up functions
  //
  // `root_network`: the per-layer root hydraulic RESISTANCES, PER UNIT LEAF AREA.
  // Not the root carbon they may have been derived from (#33). The solve reads
  // exactly two of the five fields -- r_R_H_min and r_R_V_sum -- and nothing in it
  // touches root carbon, the 1/3 : 2/3 vertical/horizontal split, the layer
  // thickness or either beta_R_* constant. Those are a root-ARCHITECTURE model,
  // and which one is in force is no more this package's business than which
  // conductance-versus-height model produced `leaf_specific_conductance_max`,
  // which plant has always passed in already reduced. `root_network_from_carbon`
  // is the architecture model that used to run in here; it is still in this
  // package, still tested, and still what the golden grid calls -- but the CALL is
  // the caller's now.
  //
  // ⚠️ PER UNIT LEAF AREA, which is where the intensiveness lives (hazard 4). The
  // leaf once took absolute root carbon plus area_leaf and divided, but the two
  // only ever appeared as that ratio -- r_R = beta/c_r is exactly linear in root
  // carbon, so E_i depends on (root carbon / area_leaf) and nothing else. That
  // reduction now happens on the far side of the boundary, so it is the
  // resistances that arrive already scaled, and a caller passing resistances built
  // from absolute carbon gets a silently wrong E_up rather than an error: no
  // check here can see the difference, because both are just positive numbers.
  //
  // `PPFD` is the INTERCEPTED photon flux per unit leaf area (see PPFD_). On the
  // energy-balance path the leaf also needs the intercepted near-infrared,
  // NIR_int_, set before this call.
  void set_physiology(const RootNetwork& root_network, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double leaf_temp, double atm_o2_kpa, double atm_kpa);

  // Replace the fourteen traits on an existing Leaf, leaving the four numerical
  // controls alone. Same arguments, same order, as the constructor's trait subset,
  // plus R_d_25 which the constructor does not take.
  //
  // `beta_R_H` and `beta_R_V` are NOT traits here since #33: they left with the root
  // architecture model, so there is no route to d(output)/d(beta_R_*) through this
  // object. A caller who needs one differences the NETWORK, which is now an input:
  // root_network_from_carbon is
  // homogeneous of degree 1 in each constant (r_R_H_min proportional to beta_R_H,
  // r_R_V to beta_R_V), so the perturbed network is a scaling of the base one and
  // costs no rebuild -- but the two solves either side of it are still two solves.
  // ⚠️ Do not read that as "the gradient is free"; only the perturbation is.
  //
  // ONE ENTRY POINT, NOT FIFTEEN SETTABLE FIELDS, and this is the set_supply_*
  // judgement again (issue #32) rather than a stylistic preference. The traits are
  // public plain doubles, so assigning one *compiles* -- and three separate pieces
  // of derived state go stale when it does:
  //
  //   * the two vulnerability CURVES. Each c owns a pre-integrated table and each
  //     P50 sets the scale (b, and the memoised 1/b) it is read at, so a bare
  //     `l.stem_c = x` leaves the hot path reading the previous curve's integral
  //     while proportion_of_conductivity reports the new one -- two answers to
  //     the same question.
  //   * the solved OPERATING POINT, which is hazard 8: an output that no code path
  //     rewrites goes on reading as though it belonged to the new traits.
  //   * vcmax_, jmax_ and R_d_, which are derived from vcmax_25/jmax_25 inside
  //     set_physiology's TEMPERATURE CACHE. ⚠️ THIS ONE IS CLOSED (#55) and the
  //     entry stays only so the repair is not undone: the key is now every scalar
  //     update_temperature_dependent_params() reads, vcmax_25 and jmax_25 among
  //     them, so "change the trait, then set the drivers again" DOES recompute
  //     them. It used to be keyed on (leaf_temp_, atm_o2_kpa_) alone and took a
  //     cache hit, which ran a whole trait sweep at the first vcmax the object
  //     ever saw and reported plausible numbers throughout. The guarantee lives
  //     in photo_temp_key(), not in the field: a member added to the temperature
  //     block and left out of that key reopens it silently.
  //
  // A fourth, which is not derived state but is the same argument: a bare write
  // bypasses the #25 positive-magnitude checks below entirely.
  //
  // So: the invariant checks run, the tables
  // are rebuilt only when the c that owns them actually moved, and
  // setup_clean_leaf() puts the object back in its just-constructed state -- which
  // resets both caches and requires set_physiology() before the next solve, exactly
  // as a fresh Leaf would. That last part is not conservatism: the derived
  // photosynthetic parameters really are unknown until the drivers are re-supplied.
  void set_traits(double vcmax_25, double stem_c, double stem_P50,
                  double root_c, double root_P50,
                  double TF24_beta2, double jmax_25, double phi_J,
                  double curv_fact_elec_trans, double curv_fact_colim,
                  double TF24_cost_scale, double R_d_25, double JS22_gamma,
                  double CMax_a, double CMax_b, double g_min);

  // THE CRITICAL POINT IS A QUANTILE OF THE CURVE, not a parameter: the fraction
  // of maximum conductivity remaining there. 0.05 gives P95. Sperry's reference
  // conductance `k_crit = f * kmax` reads the same constant, so the two cannot
  // disagree -- held as a free trait plus a hard-coded 0.05 they needed a
  // consistency check to police them.
  //
  // A convention rather than a knob: 0.05 everywhere it appears across this
  // family, and `sicangco-2026` parameterises it only to mirror the paper's own
  // `calc_Pcrit(b, c, ratiocrit = 0.05)` signature.
  //
  // ⚠️ MOVING IT IS NOT A LOCAL CHANGE, and the tempting argument for moving it is
  // that the bound it sets rarely binds -- measured, 18 of 240 feasible reference
  // rows at 25 C and NONE at 40 C, since heat drives leaves to the other bound.
  // That is true and it is only an argument about the BRACKET. The same 0.05 also
  // sets, at every operating point rather than at the rare ones:
  //
  //   * SOX's cost, `g = (f - 0.05)/(1 - 0.05)`
  //   * JW26's cost, `g = 1 - psi/psi_crit`
  //   * ProfitMax's normaliser, `k_span = k(psi_soil) - 0.05*kmax`, so its lambda
  //   * the P50 chain rule, through `dpsi_crit/dP50`
  //
  //   * the dry end of the four curves whose own domain ends at psi_crit
  //     (cost_domain_ends_at_psi_crit)
  //
  // So three of the nine cost curves change shape, not just a bound. For the
  // others psi_crit is a labelled point on a curve that continues past it (#183):
  // it bounds nothing, because their marginal cost of water diverges as the stem
  // dries, and that, not a bound, is what keeps their optimum inside.
  static constexpr double k_crit_fraction = 0.05;
  // Knot intervals in each vulnerability table (VulnerabilityTable), over
  // x = (psi/b)^c in [0, vulnerability_x_end]. Named because a CONSUMER has to
  // seat the same number -- plant's Control carries it -- and a literal repeated
  // there would drift from this one silently.
  //
  // 800 puts the table within 1.1e-14 of G(inf) of the closed form for every c
  // from 1 to 12 (1600: 1.3e-15, at boost's own accuracy). That is two orders
  // below the solve's own floor, set by psi_stem_to_ci's 1e-13, and the table's
  // derivatives are not read at all -- G' and G'' are closed forms. At 800 a
  // Leaf's two tables hold the same number of doubles the cubic splines did at
  // 400, and they are shared between copies.
  static constexpr double ncontrol_default = 800.0;

  // The least share of G(inf) the stem must still be able to carry from the
  // wettest layer for the collar solve to have an interval. G is read to
  // ~2e-16 of G(inf), so a stem whose remaining capacity G(inf) - G(psi) is
  // 1e-11 of it carries a flux resolved to ~2e-5 of itself; below ~1e-13 every
  // gradient is rounding, and the solve was measured refusing there (soil
  // 13.4-14 MPa at the defaults, where x = (psi/b)^c > 27). The leaf shuts
  // instead, at a soil of ~12.4 MPa at the defaults (x ~ 22), giving up a flux
  // below 1e-11 of the stem's capacity: continuous to that order.
  //
  // ⚠️ A RESOLUTION FLOOR, NOT A MODELLED BOUND, and it binds only where psi_crit
  // does not (the curves cost_domain_ends_at_psi_crit names leave ~1e-2 there).
  // Carrying the curve's tail as G(inf) - G would resolve it and remove this.
  static constexpr double stem_capacity_floor = 1e-11;

  // The two derivations, one place each, so a caller cannot get them
  // inconsistent. `f` is the remaining-conductivity fraction.
  // Forwards: the derivation lives in vulnerability.hpp, beside the curve it
  // scales and the bundles that differentiate it.
  static double weibull_b_from_P50(double P50, double c) {
    return phylloptim::weibull_b_from_P50(P50, c);
  }
  static double weibull_psi_at_fraction(double b, double c, double f) {
    return b * std::pow(std::log(1.0 / f), 1.0 / c);
  }

  // Seats both curves' derived scale and critical potential from the traits
  // already in place. For the CONSTRUCTORS only: set_traits has to compute the
  // same four values BEFORE it assigns anything, because it decides which
  // tables to rebuild by comparing them against the ones standing.
  void derive_vulnerability_curves_from_P50() {
    stem_b = weibull_b_from_P50(stem_P50, stem_c);
    psi_crit = weibull_psi_at_fraction(stem_b, stem_c, k_crit_fraction);
    roots_.root_b = weibull_b_from_P50(roots_.root_P50, roots_.root_c);
    roots_.root_psi_crit =
        weibull_psi_at_fraction(roots_.root_b, roots_.root_c, k_crit_fraction);
  }

  // The #25 boundary: the potentials that must be positive magnitudes. One copy,
  // called from both the constructor and set_traits.
  static void check_psi_magnitudes(double stem_P50, double stem_c,
                                   double root_P50, double root_c);

  void setup_transpiration(double resolution);

  // The stem cumulative-vulnerability integral G and its inverse, as the FOUR
  // operations the model actually performs on them. Every read goes through
  // these. G and G' are the table's value and slope, read from one polynomial so
  // that G' is the derivative of the G the solve reads; G^-1 is Newton on that
  // pair; and (G^-1)' is 1/G' at the inverse. None of them has a
  // domain end: past the end of the curve G is G(inf), and only a negative
  // potential, or a flux at or past G(inf), is refused.
  //
  // The two non-derivative reads take an optional `caller` label, which is
  // appended to an out-of-domain failure, so a refusal names the read that made
  // it (localising plant#576 came down to which).
  double stem_curve_integral(double psi, const char* caller = nullptr) const;
  double stem_curve_integral_deriv(double psi) const;
  // `converge` is kept for its callers; the inverse always converges to a few
  // ulps now.
  double stem_curve_integral_inverse(double w, const char* caller = nullptr,
                                     bool converge = false) const;
  double stem_curve_integral_inverse_deriv(double w) const;
  // G(inf) = stem_b Gamma(1/c) / c, the most the stem can carry per unit
  // conductance.
  double stem_curve_limit() const;

  // Move stem_P50, and with it stem_b and psi_crit, without touching the stem
  // table, which depends on stem_c alone.
  //
  // ⚠️ FOR DERIVATIVE WORK ONLY. It deliberately does NOT clear the solved
  // operating point -- re-seating the physiology is exactly the cost being
  // avoided -- so whatever reads the outputs afterwards must write them first.
  // `set_traits()` is the way back.
  void perturb_stem_P50(double stem_P50_new);
  // Forwards to roots_.setup_vulnerability. Kept on Leaf because it is part of
  // the published construction sequence (see the umbrella header) and plant's
  // bindings name it.
  void setup_root_vulnerability(double resolution) {
    // Recorded for the reason setup_transpiration records it: `set_traits`
    // rebuilds both curves at this member.
    vulnerability_curve_ncontrol = resolution;
    roots_.setup_vulnerability(resolution, stem_curve_);
  }
  void setup_clean_leaf();

  // Absolute tolerance for the direct quadrature in
  // transpiration_full_integration (leaf/quadrature.hpp). Not used anywhere on
  // the production path, which reads the pre-integrated table instead.
  double integration_tol_ = 1e-8;

  // Retained for API compatibility with plant's R bindings, which expose this
  // method. In plant it configured a compiled adaptive Gauss-Kronrod integrator
  // (plant::quadrature::QAG); that dependency is gone, so `integration_tol` now
  // sets the tolerance of the header-only adaptive Simpson quadrature and
  // `integration_rule` is accepted but ignored -- Simpson has no rule order to
  // choose. Only affects transpiration_full_integration.
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-8) {
    static_cast<void>(integration_rule);
    integration_tol_ = integration_tol;
  }

  // Medlyn stomatal-conductance model (from develop #450); R-callable, standalone.
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
  // std::vector<double> root_collar_psi(std::vector<double> soil_moist_);

  // Every psi crossing this boundary is a POSITIVE MAGNITUDE in MPa (#25). The
  // entry points that take the soil state as an *argument* are reachable from R,
  // where the pre-#25 signed vector would compile, run, and be silently wrong --
  // so they check.
  //
  // The check is skipped when the caller handed back this object's own vector,
  // which set_physiology already validated. That is not a micro-optimisation:
  // E_from_Soil_to_Root_Collar runs ~10^3 times per collar solve (hazard 5), and
  // an O(layers) scan there would be paid on every one of them. The address
  // compare costs nothing and is the same idiom MultiLayerRoots::uptake_at uses
  // to select its cached path.
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

  // Uptake at a collar suction against an arbitrary vector of layer suctions.
  // Thin forwarder to roots_.uptake_at; the E_up_ / soil_consumption_ buffers
  // stay on Leaf and are handed over by reference, because plant writes back
  // into them by name after crown integration.
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

  // What the leaf hands over: the objective plant bills, the water it drew
  // from each layer, and the share of that water the cuticle lost. At any
  // scalar, because all three re-enter plant's chain.
  //
  // The three fluxes are named so a consumer never takes one's value with
  // another's rows (plant#655): the soil draw is the draw's E_up, the stem
  // carries the same total on the collar route (the water budget closes,
  // E_up = E), the residual flux is `residual_flux` here, and the stomatal flux
  // is the total less it. Each carries rows that match its own value.
  template <class T>
  struct BilledOutputs {
    T profit{};
    std::vector<T> uptake;
    T residual_flux{};
  };

  // The soil's draw at ONE collar, with the collar it was taken at, so a flux
  // from one and a slope from another cannot be paired. check_draw enforces that
  // at every read rather than trusting the caller.
  //
  // ⚠️ TWO DIFFERENT PAIRINGS, DELIBERATELY. `E_up` (the total uptake) and its
  // collar slope are both active, so they are a pair. The per-layer draws are
  // active but their collar slopes are SUPPLIED at double -- taking the supply
  // again at a live collar would record the whole of it a second time -- so
  // those are implicit_value's shape, not the pair's, and flattening both into
  // one would lose that.
  template <class S>
  struct SupplyDraw {
    double at = 0.0;
    // Zero, not default-initialised: at a hydraulic shutdown supply_draw_at
    // returns before writing it, and a default-initialised double is
    // indeterminate.
    odelia::with_slope<S> E_up{};
    std::vector<S> uptake;
    std::vector<double> duptake_dp;
    // The most the roots can deliver from this soil, which caps the residual
    // flux on the collar route (#182, #184). Carried only where that cap is in
    // force -- a residual conductance on the layered supply -- so a draw at
    // g_min = 0 records nothing more than it did.
    bool has_uptake_limit = false;
    S uptake_limit{};
    // Where no layer conducts (every rooted potential past the root integral's
    // cap) the wet bound is the wettest layer's gravity balance, psi_i + z_i,
    // and its rows are that layer's potential's: the residual is zero at every
    // collar there, so the implicit function theorem has no slope to divide by.
    // Carried only in that case.
    bool no_conduction = false;
    S gravity_balance{};

    // ⚠️ WHAT CARRIES A ROW, NAMED FOR visit_active. A shape that dispatch does
    // not open is skipped in silence, and here that arrives as an exact zero in
    // a gradient column rather than as an error -- so this exists to make the
    // theorem's inputs reachable, not for convenience. `at` and `duptake_dp` are
    // double and carry no row. Const because the harvest reads slots.
    // A member added above and not visited here would be such a column, so the
    // layout is asserted below this struct: adding a member fails to compile
    // until both are updated.
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
  // The two optional rows (#182), the roots' ceiling and the dead soil's
  // gravity balance, are each a flag and a value; at double the flag pads to
  // the value's alignment, so each pair counts as two doubles.
  static_assert(sizeof(SupplyDraw<double>) ==
                    sizeof(double) + sizeof(odelia::with_slope<double>) +
                        sizeof(std::vector<double>) * 2 +
                        2 * (alignof(double) + sizeof(double)),
                "SupplyDraw gained or lost a member: update for_each_active, "
                "then this count");
  static_assert(sizeof(odelia::with_slope<double>) == 2 * sizeof(double),
                "odelia::with_slope is no longer a value and a slope: re-check "
                "that SupplyDraw::for_each_active visits all of it");

  // How a draw is meant to be made: it records the collar it was taken at, so its
  // values cannot come from different ones, and operating_collar refuses a draw whose
  // collar is not the solve's.
  //
  // ⚠️ THE COLLAR IS PASSIVE HERE, whatever the caller hands in. This is the
  // supply at a POINT; where the point itself moves is the supplied row's business, and
  // taking it live would record the search that placed it.
  template <class S>
  SupplyDraw<S> supply_draw_at(const S& collar,
                               const SupplyInputs<S>& supply) const {
    // The walk below is MultiLayerRoots', whatever path the leaf solved on, so
    // a single-potential leaf would be handed a draw from a supply its solve
    // never read.
    if (supply_kind_ != SupplyKind::MultiLayer) {
      util::stop("Leaf::supply_draw_at: the derivative surface reads the "
                 "multi-layer supply, and this leaf is on the " +
                 supply_kind_name() + " path");
    }
    const std::size_t n = static_cast<std::size_t>(supply_n_layers());
    SupplyDraw<S> d;
    d.at = odelia::util::to_passive(collar);

    // A shut leaf sits on its wet bound -- the collar at which uptake meets the
    // residual flux -- so it is drawn like any other point: the uptake is the
    // residual flux (zero at g_min = 0) and the rows are the bound's own.

    // ONE WALK for all three. The per-layer collar slopes are SUPPLIED at
    // double: the draws already carry their own rows, and what is missing is
    // only the channel through the collar moving -- taking the supply again at
    // a live collar would record the whole of it a second time for the same
    // numbers.
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

    // ⚠️ A KINK IS ANSWERED HERE, ONCE, so that no consumer can read the sentinel.
    // `duptake_dpsi` returns NaN to mean "the analytic branch is not valid
    // across this collar, difference it instead". A wet bound is one -- for
    // a single rooted layer it is exactly the gravity-balance point, and that is
    // where a shut leaf's collar sits at g_min = 0.
    //
    // Honoured at the total and NOT at the per-layer slopes, the sentinel became
    // `NaN * step` in every layer's uptake, NaN in the value and in the row. That
    // reached plant, where one NaN row takes every trait column that shares the
    // adjoint accumulator -- and nothing named it, because the row plant hands to
    // `supply_derivatives` is the constant 1.0 and the sentinel is already
    // on the tape behind it.
    //
    // The rows go with the sentinel, which is right rather than a loss: the
    // analytic form is not valid here, so neither are its rows. Each slope is
    // differenced from the function its own values come from, so the total keeps
    // its kg and the layers keep their mol without a conversion to remember.
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

  // ⚠️ A DRAW FROM THE WRONG COLLAR IS A WRONG NUMBER, not a missing row: profit
  // reads the flux at whatever collar it is evaluated at, so a stale draw answers
  // with the uptake from somewhere else and every value stays finite. Exact,
  // because the draw stores the point it was taken at rather than deriving it.
  template <class S>
  void check_draw(double expect, const SupplyDraw<S>& draw) const {
    if (draw.at != expect) {
      util::stop("Leaf: the supply draw was taken at a collar of " +
                 util::to_string(draw.at) + " and is being read at " +
                 util::to_string(expect));
    }
  }

  // The differentiable surface's coordinates.
  //
  // A coordinate is a value AND its response in the collar potential, so the two
  // travel as one. Written as separate scalars they were seeded from separate
  // places, and a first attempt that moved only the collar disagreed with a
  // difference by 40 to 100 per cent -- dM/dp is directional, and five
  // coordinates move together or the assembly differentiates nothing.

  // What profit needs: the stem potential the collar determines, and the
  // intercellular CO2 that goes with it.
  template <class T>
  struct CollarCoords {
    odelia::with_slope<T> sigma;
    odelia::with_slope<T> ci;
  };

  // Which cost curve a psi_stem derivative differentiates. The cost enters the
  // chain through exactly ONE quantity -- dC/dpsi_stem -- so this selects that
  // and nothing else.
  enum class CostCurve { TF24, CF77, JS22, CMax, SOX, JW26, ProfitMax, TF24_floor,
                         LeastCost };

  // ⚠️ THE ONE ABSTRACTION THAT MAKES EVERY MODEL THE SAME MODEL.
  //
  // All eight maximise `h(A(psi)) - C(psi)`: a BENEFIT LINK `h` composed with
  // assimilation, minus a cost curve. That is not a convenience -- it is why the
  // derivative below is one expression rather than three, and `h'` is the only
  // thing that varies:
  //
  //     d/dpsi [ h(A) - C ] = h'(A) * dA/dpsi - dC/dpsi
  //     ----------------------------------------------------------------------
  //     Identity  h(A) = A            h' = 1          TF24, CF77, JS22, CMax,
  //                                                   TF24_floor
  //     Log       h(A) = log A        h' = 1/A        SOX, JW26  (products: A*g
  //                                                   and log A + log g share
  //                                                   an argmax), LeastCost
  //                                                   (a RATIO: minimising
  //                                                   (c_E E + c_V V)/A is
  //                                                   maximising log A -
  //                                                   log(c_E E + c_V V))
  //     Scaled    h(A) = A/|A|max     h' = 1/|A|max   ProfitMax
  //
  // With the energy balance on, dA/dpsi also runs through Tleaf(E(psi)). That
  // path is part of dA/dpsi and takes the same single h': there is no h'' in a
  // first derivative, because h'(A) is never itself differentiated. Two costs
  // read Tleaf, ProfitMax's thermal term and LeastCost's capacity term through
  // `vcmax_`. Both are cost-side and take none. All three slopes live in
  // `dprofit_energy_balance_term`.
  //
  // ⚠️ THE LINK IS NOT ONLY ON THE BENEFIT, AND READING IT THAT WAY IS THE NATURAL
  // MISTAKE. For the two product curves the COST is stored in log form as well --
  // `C = -log g` -- so `h(A) - C` is `log(A*g)`, the transform of the whole
  // objective rather than of half of it. That asymmetry in the NAMING is real
  // though: `C` is a function of psi alone, so its transform is just another cost
  // curve and needs no declaration, while `A` comes out of the `ci` root-find and
  // cannot be reparameterised that way. The per-curve notes above
  // `optimise_psi_stem_single` carry the two lines of algebra showing the
  // resulting first-order condition is the product rule divided by `A*g`, and so
  // has identical roots.
  //
  // ⚠️ `Scaled` treats `|A|max` as CONSTANT, which makes ProfitMax's trait
  // gradients PARTIALS at fixed normaliser rather than total derivatives. That is
  // deliberate and it is what the code computes: `|A|max` comes from a scan over
  // the supply stream, so its argmax is piecewise constant in the traits and its
  // derivative is a sequence of zeros and jumps. A total derivative would need
  // that chain rule; a partial at fixed |A|max is well defined and is the honest
  // thing to report until someone needs the other.
  enum class BenefitLink { Identity, Log, Scaled };

  // ⚠️ ONE BOUND, derived from the enum's last member rather than named. The
  // bounds checks below were written as `> CostCurve::JW26` and went stale the
  // moment ProfitMax was appended -- the curve existed, had a link and a
  // derivative, and was simply invisible to R because `curve_name()` called it
  // unknown. Keep this the last member.
  static constexpr int n_cost_curves = static_cast<int>(CostCurve::LeastCost) + 1;

  // ⚠️ THE MODEL IS CONFIGURATION, NOT A CALL ARGUMENT, and the reason is that its
  // PARAMETERS already are. `CF77_lambda_` and `TF24_floor_lambda_o` are fields;
  // `JS22_gamma`, `CMax_a` and `CMax_b` are traits. Naming the curve
  // at the call site while the constants it
  // needs live on the object is incoherent -- you cannot configure half a model.
  // So the curve is seated here, beside them, and `optimise()` takes nothing.
  //
  // Defaults to TF24 and the collar route, which is what plant runs, so an existing
  // caller that never sets a model gets the production path unchanged.
  CostCurve cost_curve_ = CostCurve::TF24;
  bool route_is_collar_ = true;
  // THE THIRD AXIS: which METHOD reaches the operating point the seated curve and
  // route define. `false` is the exact solve -- the root-find of the first-order
  // condition over the closed interval -- and `true` is the closed form in
  // `closed_form.hpp`, which inverts the Medlyn USO form given a wet-end power law
  // for lambda(psi).
  //
  // ⚠️ IT IS A METHOD, NOT A MODEL, and that is the whole reason it is a third axis
  // rather than an eighth curve: the objective is identical and only the route to
  // its argmax differs. So the two must agree wherever the closed form is valid,
  // which is what `tests/testthat/test-closed-form.R` measures, and a disagreement
  // is an error in the approximation rather than a different answer to a different
  // question.
  //
  // Defaults to the exact solve. Every existing caller therefore gets a
  // bit-identical answer, and that is asserted rather than asserted-in-a-comment:
  // `optimise()` returns to its previous body untouched when this is false.
  bool method_is_closed_ = false;

  // Whether the LAST `optimise()` fell back from the closed form to the exact
  // solve, and the running counts that make the fallback fraction phi measurable.
  //
  // ⚠️ WHY THIS IS INSTRUMENTED RATHER THAN INFERRED. The closed form's validity
  // guard tests an OUTPUT (ci/ca > 0.5), so it can only be applied after solving
  // and the fallback cannot be predicted from the drivers. The realised speedup is
  // therefore not the ceiling but `1 / [phi + (1-phi)/speedup]` -- 3.6x at
  // phi = 0.2 against a 10.8x ceiling -- and quoting the ceiling for a
  // water-limited scenario is the mistake these three fields exist to prevent.
  // `closed_form.hpp`'s open item 3 was that phi had never been measured on a real
  // scenario; it is measurable now.
  //
  // ⚠️ `last_solve_fell_back_` IS WRITTEN BY EVERY `optimise()`, INCLUDING THE
  // EXACT ONE (hazard 8: an output a branch declines to write becomes the previous
  // solve's value, and plant reuses one Leaf per individual). It is a bool, so
  // writing it on the exact path costs nothing numeric and every golden baseline
  // is untouched.
  bool last_solve_fell_back_ = false;
  int closed_form_calls_ = 0;
  int closed_form_fallbacks_ = 0;

  void find_root_collar_psi();
  // The same solve for ANY cost curve. `find_root_collar_psi()` is the TF24
  // instantiation, kept under that name because plant's generated glue and the R
  // bindings call it.
  template <CostCurve K> void find_root_collar_psi_for();
  // ⚠️ THE ONE ENTRY POINT. Takes the curve BY NAME and the route by name, so a
  // caller never handles a curve index -- an integer that silently means a
  // different model when the enumeration grows is the kind of API this package has
  // been removing. Everything else here is either an alias over this or the
  // compile-time form the hot path uses.
  //
  // `route` is "collar" (the production formulation: the full soil-to-collar path,
  // any supply) or "stem" (psi_stem with the upstream potential pinned at psi_soil,
  // the form the literature is written in, single soil potential only).
  //
  // `method` is "exact" (the root-find of the first-order condition, the default and
  // the only thing any existing caller gets) or "closed" (the closed form in
  // `closed_form.hpp`). The third argument is DEFAULTED rather than required so that
  // every two-argument call site -- and the R binding, which carries the same
  // default -- is untouched.
  void set_model(const std::string& curve, const std::string& route,
                 const std::string& method = "exact");
  // The same, pre-parsed, for a C++ caller. See the definition.
  //
  // ⚠️ THE THIRD ARGUMENT DEFAULTS TO THE EXACT SOLVE, WHICH MEANS SEATING RESETS
  // IT. `gradient::route_seat` seats once per observation through this overload, so
  // a gradient always differentiates the exact solve even on a leaf the caller had
  // put on the closed form. That is deliberate -- the implicit-function-theorem
  // composite is derived from the exact first-order condition and says nothing
  // about an approximation to its root -- but it is silent, so it is stated here.
  void set_model(CostCurve curve, bool collar_route, bool closed_method = false);
  // The curve, route and method currently seated, by name.
  std::string model_curve() const { return curve_name(static_cast<int>(cost_curve_)); }
  std::string model_route() const { return route_is_collar_ ? "collar" : "stem"; }
  std::string model_method() const { return method_is_closed_ ? "closed" : "exact"; }
  // Solve for the operating point under the seated model. THE entry point.
  void optimise();
  // The closed-form arm of `optimise()`, with every refusal the closed form needs.
  //
  // ⚠️ DEFINED IN `closed_form.hpp`, NOT HERE, and the include order is what makes
  // that legal: `closed_form.hpp` includes this header (it needs the complete
  // `Leaf`), and `<phylloptim.hpp>` includes `closed_form.hpp`, so in any
  // translation unit that takes the umbrella the definition follows this
  // declaration. It is a plain `inline` member rather than a template, so there is
  // no instantiation subtlety: a translation unit that includes only
  // `leaf_model.hpp` and calls this gets a link error, which is a legible failure
  // rather than a silent one.
  //
  // Not called directly. `optimise()` reaches it when the seated method is
  // "closed"; the fallback to the exact solve is inside it, because the validity
  // guard tests an output.
  void optimise_closed();
  // Fraction of closed-form `optimise()` calls that fell back to the exact solve
  // since the counters were last reset -- the phi in `1/[phi + (1-phi)/speedup]`.
  // NaN before any closed-form call, so a caller cannot mistake "never asked" for
  // "never fell back".
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
  // Shared setup for the root-collar solve: builds the soil-side caches, handles
  // every feasibility early-exit (shutdown / closed interval) by
  // setting the final operating point itself, and otherwise returns the feasible
  // collar-potential interval [bound_a, bound_b] (positive magnitudes). Returns
  // false when the operating point is already fully determined (caller is done),
  // true when there is a real interval to choose a collar potential within.
  //
  // `first_bracket` lets a curve not bounded at psi_crit try [wet bound, the
  // continuity root with the stem at psi_crit] first, which holds the optimum
  // on most rows at the old bracket's cost; maximise_profit_over_collar extends
  // it to the end of the stem curve when profit is still rising there. Only the
  // optimising solve asks for it: a prescribed collar needs the whole domain.
  template <CostCurve K>
  bool prepare_collar_solve(double& bound_a, double& bound_b,
                            bool first_bracket = false);
  // Move the bracket's dry end from psi_crit to the end of the stem curve.
  double extend_collar_bracket(double bound_a);
  // At a dry pin of a curve psi_crit does not bound, the stem is held at the end
  // of its curve, which carries the uptake (G(inf) - G(collar), to the
  // continuity root's precision): seat it there, rather than where the inverse
  // lands, which is anywhere G rounds to G(inf). Returns whether it did.
  template <CostCurve K> bool seat_stem_at_dry_pin();
  // Evaluate the leaf at a *given* root-collar potential (positive magnitude)
  // rather than optimising it: reuses prepare_collar_solve, refuses a target
  // outside the feasible interval (infeasibility code `prescribed_potential`),
  // and evaluates there (no golden-section search). Leaves
  // exactly the same outputs as find_root_collar_psi and returns profit_. Used by
  // TF24f's gradient-ascent acclimation (#525).
  double evaluate_root_collar_psi(double target_opt_root_psi);
  template <CostCurve K>
  double evaluate_root_collar_psi_for(double target_opt_root_psi);
  // Evaluate profit at a given root-collar potential (positive magnitude)
  // *assuming prepare_collar_solve has already run this step* (soil-side caches
  // built, feasible interval [bound_a, bound_b] known). Refuses a target outside
  // the interval and otherwise sets the operating point (opt_psi_stem_, opt_root_psi_,
  // profit_), returning profit_. This is the post-prepare body of
  // evaluate_root_collar_psi, factored out so the centred finite-difference leaf
  // solve can share one prepare_collar_solve across its three profit evals (#530).
  template <CostCurve K>
  double profit_at_collar_psi(double target_opt_root_psi,
                              double bound_a, double bound_b);
  // Exact d(profit)/d(opt_root_psi) at a given root-collar potential (positive
  // magnitude), for TF24f's acclimation tracking (#525/#527). Combines
  // forward-mode AD for the analytic photosynthesis/cost algebra, the
  // implicit-function theorem at the psi_stem_to_ci root-find, and the
  // vulnerability tables' own slopes for the smooth transport. Replaces
  // the noisy finite-difference gradient. Seats the soil-side caches itself, so a
  // solve need not have run first.
  //
  // ⚠️ **The 0.0 returned on the shut-down / reversed-gradient exits is a
  // SENTINEL, not a stationary point**, and the distinction only became load
  // bearing because the collar solve root-finds on `dprofit == 0`. It matters
  // because the sentinel fires at `prepare_collar_solve`'s WET bracket endpoint
  // -- at `root_zero_E` uptake is zero by construction, so `psi >= psi_stem` --
  // which is the first point a bracketing solver evaluates when it checks that
  // its bracket brackets. Measured at the default operating point: profit there
  // is -1.897 against 2.516 at the true optimum, so a solver that reads the
  // sentinel as a root returns the zero-transpiration point as the answer. The
  // region is narrow (at most 3.46e-07 MPa into the bracket over the golden grid,
  // median 1.22e-08), which is exactly why it would survive casual testing.
  //
  // `feasible`, when non-null, reports whether the point admits an informative
  // gradient at all, so a caller searching for a zero can tell the two kinds of
  // 0.0 apart. The out-parameter rather than a NaN return is deliberate three
  // ways: TF24f consumes the return value directly as an ODE rate
  // (`dpsi/dt = k_acclim * dprofit`, plant/src/tf24f_strategy.cpp), so a NaN
  // would propagate into plant's state vector where 0.0 correctly means "do not
  // acclimate this step"; NaN already carries a *different* contract on the
  // neighbouring derivative (hazard 6: `duptake_dpsi` returning NaN means "fall
  // back to finite differences"), and overloading it would make the two
  // indistinguishable; and defaulting to nullptr leaves every existing call site
  // and the generated R binding untouched.
  double dprofit_droot_collar_psi(double opt_root_psi, bool* feasible = nullptr);
  template <CostCurve K>
  double dprofit_droot_collar_psi_for(double opt_root_psi,
                                      bool* feasible = nullptr);
  // The same, with the flag as a return rather than an out-parameter, because
  // R cannot reach a `bool*` -- RcppR6 has no form for one, so the generated
  // binding drops it and every R-side composite silently gets the "a composite
  // that ignores it inherits the bug" case this header warns about above.
  // `{dprofit, feasible}`, in that order.
  std::vector<double> dprofit_droot_collar_psi_checked(double opt_root_psi);
  // Post-prepare body of dprofit_droot_collar_psi, with the same `feasible`
  // contract. Assumes the supply path's per-solve caches are already seated, so
  // the collar solve can share ONE supply_begin_solve across all ~10 of its
  // gradient evaluations instead of re-seating per call -- the same saving #530
  // made for the finite-difference path, and it matters here because
  // begin_solve() is a vulnerability read per soil layer.
  template <CostCurve K>
  double dprofit_at_collar_psi(double opt_root_psi, bool* feasible = nullptr);


  // dprofit/dpsi_stem for the solvers that optimise psi_stem directly with the
  // upstream potential held fixed.
  //
  // The same chain as dprofit_at_collar_psi, and strictly simpler: there is no
  // collar-to-stem map, so dpsi_stem/dx is 1 and every term in the upstream
  // potential drops out. A template rather than a runtime switch because it is
  // called inside a root-find, and because a compile-time choice cannot disagree
  // with the objective the caller is maximising.
  //
  // ⚠️ `feasible` carries the same contract as the collar version: false means
  // "no informative gradient here", and the 0.0 returned alongside it is a
  // SENTINEL rather than a stationary point. A caller that ignores it and hands
  // the zero to a root-find gets the shut-down state reported as the optimum.
  template <CostCurve K>
  double dprofit_dpsi_stem(double psi_stem, double psi_upstream,
                           bool* feasible = nullptr);

  static std::string curve_name(int curve);

  // ⚠️ THE DERIVATIVE SURFACE READ ON ANOTHER CURVE IS A WRONG NUMBER, like a
  // draw from another collar (check_draw): operating_collar and billed_outputs_at evaluate
  // whatever K names at the solve's collar, and every value stays finite.
  // Checked against the SEATED curve, which is what optimise() and
  // replay_operating_point() solve. find_root_collar_psi() is the exception --
  // it solves TF24 whatever is seated -- so a leaf seated on another curve and
  // solved through it passes this check for the wrong reason.
  template <CostCurve K>
  void check_surface_curve(const char* who) const {
    if (K != cost_curve_) {
      util::stop(std::string(who) + ": asked for the " +
                 curve_name(static_cast<int>(K)) + " curve, and the leaf is "
                 "seated on " + curve_name(static_cast<int>(cost_curve_)) +
                 "; read the surface on the curve the solve used");
    }
  }

  // ⚠️ THE ONLY RUNTIME CURVE DISPATCH IN THE CLASS, and the reason it is a
  // template rather than four switches. Calls `f(tag)` where `tag` is an
  // `std::integral_constant<CostCurve, K>`, so the body it hands back is
  // COMPILE-TIME specialised on the curve -- the same instantiation the hot path
  // gets -- while the choice is made at run time.
  //
  // There were four of these switches (solve-collar, solve-stem, evaluate-stem,
  // dprofit-stem), all with the same seven arms, so adding a curve meant
  // remembering six places. Now it is two: an arm here and an arm in
  // `curve_name`. Both are switches over the enum with no `default`, so
  // `-Werror=switch` makes a forgotten curve a build failure -- which is the
  // whole reason neither is an `if`/`else` chain.
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
    // Unreachable for any enumerator, and not a `default:` arm -- a `default`
    // would satisfy `-Werror=switch` and so remove the check this exists for.
    util::stop("unknown cost curve index " +
               util::to_string(static_cast<int>(curve)));
    return f(std::integral_constant<CostCurve, CostCurve::TF24>{});
  }
  // Evaluate the leaf at a prescribed psi_stem under the SEATED cost curve (see
  // set_model). Returns profit_.
  //
  // ⚠️ NO CURVE ARGUMENT, deliberately, and it used to take an integer INDEX.
  // The curve's own constants live on the object, so naming it at the call site
  // configured half a model in each place -- and an index silently means a
  // different model the moment the enumeration grows. Seat the model once with
  // set_model(), then call this.
  double evaluate_psi_stem_at(double target_psi_stem);
  // Refuse a stem-route method on a leaf seated on the collar route. See the
  // definition for why the route and not only the curve has to be checked.
  void require_stem_route(const char* who) const;
  // ⚠️ NO `psi_upstream` ARGUMENT, deliberately. On a stem route the upstream
  // potential is ALWAYS psi_soil -- that is what "non-root-based" means -- so it is
  // read here through the same accessor `optimise_psi_stem_single` uses. Taking it
  // from the caller invited a mismatch, and the first caller got it wrong: R read
  // `psi_soil_`, which is the MULTI-LAYER roots' vector and is empty on the single
  // path, so the subscript was out of bounds rather than merely inconsistent.
  //
  // The CURVE is read from the seated model too, for the reason
  // `evaluate_psi_stem_at` gives.
  std::vector<double> dprofit_dpsi_stem_checked(double psi_stem);

  template <CostCurve K> static constexpr BenefitLink benefit_link();
  template <CostCurve K> double benefit_link_deriv(double A) const;
  template <CostCurve K> void optimise_psi_stem_single();
  // ProfitMax's shut state. Its objective is A/|A|max minus a normalised cost,
  // so with |A|max <= 0 there is nothing to evaluate -- unlike every other curve,
  // whose objective is defined everywhere. The objective is zeroed; the leaf's
  // own state (its fluxes, conductance, ci and carbon) is the evaluator's at
  // the wet end, so it is written rather than left stale (hazard 8).
  void set_profitmax_degenerate_state(double psi_stem, double psi_upstream);
  template <CostCurve K> void check_cost_parameters();
  template <CostCurve K> double cost_deriv(double psi_stem,
                                           double psi_upstream);
  template <CostCurve K> static constexpr bool cost_reads_upstream();
  template <CostCurve K> static constexpr bool cost_domain_ends_at_psi_crit();
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
  // difference of the block. ⚠️ It restores every member the block writes; a
  // member added to update_temperature_dependent_params must be added here and
  // to dprofit_energy_balance_term, or the leaf stays 1e-3 degrees off.
  template <class F> double temp_deriv(double leaf_temp, F read);
  // dVcmax/dTleaf, the same way and with the same save-and-restore. LeastCost's
  // capacity cost reads `vcmax_`, so `lambda_LeastCost` needs this slope with the
  // energy balance on. dprofit_energy_balance_term reads the same two
  // evaluations inline, beside A_T.
  double vcmax_temp_deriv(double leaf_temp);
  // The energy-balance contribution to dJ/dpsi, zero when the gate is off: the
  // benefit's motion through Tleaf(E(psi)) under the curve's link, plus the
  // cost's for the two curves whose cost reads Tleaf (ProfitMax's thermal cost,
  // LeastCost's Vcmax). Kept out of line so that adding it cannot
  // change FMA contraction in the inlined gate-off path; the derivation and the
  // two sign checks are at the definition.
  //
  // `E_s` and `r_star` are the stomatal and the residual flux at the point,
  // `drstar_dpsi` the residual's slope in the collar (zero on a stem route),
  // and `dgc_dT_residual` the residual's contribution to dgc/dTleaf; see the
  // definition for where each enters. All four are exact zeros at g_min = 0.
  template <CostCurve K>
  double dprofit_energy_balance_term(double ci, double gc, double g_ci,
                                     double inv_atm, double gc_const,
                                     double A_prime, double dgc_dpsistem,
                                     double dgc_dpsi, double dpsistem_dpsi,
                                     double dT_dE, double Tleaf, double E_s,
                                     double r_star, double drstar_dpsi,
                                     double dgc_dT_residual);
  // The profit-maximising collar potential within [bound_a, bound_b], by a
  // safeguarded root-find on dprofit == 0. Returns a bound when the
  // optimum is pinned to it, which is the case on 42 of the 240 feasible
  // golden-grid rows and so is a branch that has to be written rather than a
  // corner. Falls back to golden section if either endpoint has no usable
  // gradient -- see the definition for why that fallback should never fire.
  template <CostCurve K>
  double maximise_profit_over_collar(double bound_a, double bound_b);
  // Analytic d(E_up_)/d(collar suction) -- a CONDUCTANCE, positive by
  // construction now that both sides are magnitudes (#25). Thin forwarder to
  // roots_.duptake_dpsi; see there for the derivation and for the NaN-at-a-kink
  // contract. Used only on the TF24f acclimation gradient path, not the base
  // TF24 value path.
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
  // at which uptake meets the residual flux, the stem where it carries that
  // flux, and every output from the ordinary evaluator there. See the
  // definition.
  template <CostCurve K>
  void set_shut_state();
  // The collar at which uptake equals the realised residual flux, on a bracket
  // that is widened past the stem's critical potential when the soil is too
  // dry to deliver it inside; the wet bound of every collar solve.
  double find_wet_bound(double wettest_soil_layer, double E_up_at_wettest,
                        double E_up_at_crit);
  // The same bound for the shut leaf, bracketed from its floor: see the
  // definition for why this is both cheaper and continuous through the root
  // integral's cap.
  double find_shut_wet_bound();
  // The stem potential at which the stem carries exactly the residual flux
  // from a pinned upstream potential; the wet end of every stem route. The
  // upstream potential itself at g_min = 0.
  double stem_wet_bound(double psi_upstream);
  // The root on [wettest, end], with the stem held at `end` for the continuity
  // root (find_root_crit = 1) and the zero-uptake root otherwise.
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit,
                       double f_wettest, double f_crit, double end);
  // The same root, taking its own endpoint values, with the stem at psi_crit:
  // for a caller outside the solve.
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit) {
    const auto f = [&](double x) {
      return find_root_crit == 1 ? E_column(x, psi_soil, psi_crit)
                                 : E_column_zero(x, psi_soil);
    };
    const double f_wettest = f(wettest_soil_layer);
    return find_root_psi(wettest_soil_layer, psi_soil, find_root_crit, f_wettest,
                         f(psi_crit), psi_crit);
  }
  // A bracket bound find_root_psi placed at its 1e-4 tolerance, solved again to
  // collar_root_tol from where it stopped. Called only where the bound IS the
  // answer; see its definition for why the answer needs it.
  double polish_bracket_bound(double x0, const std::vector<double>& psi_soil,
                              int find_root_crit);
  double find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil);
  double E_column(double x, const std::vector<double>& psi_soil, double psi_leaf);
  double E_column_zero(double x, const std::vector<double>& psi_soil);
  
  // The two temperature responses, at any scalar. Both are exactly LINEAR in
  // ref_value and every other argument is fixed physiology, so on the active path
  // the response is a passive factor and the reference value is what carries the
  // row. Templated rather than factored into that product because the operation
  // order is what keeps the double instantiation bit-identical to the member it
  // replaced.
  // ⚠️ ONLY `ref_value` IS TEMPLATED, AND THE REST ARE double BECAUSE THEY DO
  // NOT MOVE. An activation energy, a leaf temperature and the two deactivation
  // constants are the block's own numbers, not the caller's traits -- lifting
  // one into T gives it a derivative vector of zeros that every operation then
  // carries. Measured at a nested tangent, the two curves below were 25.7% of
  // the leaf's instructions with the constants lifted and are a single active
  // multiply without.
  template <typename T>
  T arrh_curve(double Ea, const T& ref_value, double leaf_temp) const;
  template <typename T>
  T peak_arrh_curve(double Ea, const T& ref_value, double leaf_temp, double H_d,
                    double d_S) const;

  // --- Penman-Monteith leaf energy balance (minimal core; #523) ---------------
  // Recompute the temperature-dependent photosynthetic parameters (vcmax_,
  // jmax_, gamma_, ko_, kc_, R_d_, km_, electron_transport_, and the
  // thermal-damage layer's k_imp_eff_, k_rec_eff_ and functional_) at a given leaf
  // temperature. Extracted verbatim from the inline block in set_physiology so
  // the non-PM path is bit-identical; on the PM path it is called per
  // operating point with the energy-balance Tleaf (defeating the photo_temp
  // cache, as intended).
  void update_temperature_dependent_params(double leaf_temp);
  // Saturation vapour pressure es(T) (kPa) and its slope Delta(T) (kPa K^-1),
  // Tetens formula. Not wired into the minimal-cut solve (prescribed atm_vpd is
  // kept); provided and unit-tested for the leaf-to-air VPD in the full cut.
  double saturation_vapour_pressure(double temp) const;
  double saturation_vapour_pressure_slope(double temp) const;

  // Seat `vpd_leaf_` at a leaf temperature. THE ONE PLACE the leaf-to-air deficit
  // is defined, called wherever a leaf temperature is established: set_physiology
  // (at Tair, so off-path callers get atm_vpd_ back exactly),
  // set_leaf_states_rates_from_psi_stem, and dprofit_at_collar_psi.
  //
  // ⚠️ NOT folded into update_temperature_dependent_params, even though that is
  // the other function taking a leaf temperature. The two finite-difference
  // blocks call that one at T +/- h and restore the photosynthesis members by
  // assignment; `vpd_leaf_` must NOT move with them, because dA/dTleaf there is
  // taken at fixed ci and assim_colimited() reads no VPD at all. The VPD route
  // into the derivative is carried separately, by dgc_dT.
  void set_leaf_vpd(double leaf_temp);
  // Explicit leaf energy balance:
  // Tleaf = Tair + (Rn - lambda*E) / (rho*cp/ra + g_rad), with Rn the ISOTHERMAL
  // net radiation and g_rad the linearised longwave (constants.hpp).
  // E is the hydraulically-pinned transpiration (kg H2O m^-2 s^-1); no PM
  // inversion and no A->E feedback, so this is a single algebraic forward pass.
  // `dT_dE`, when non-null, receives dTleaf/dE at the same point: one constant,
  // since the balance is affine in E. Nothing clamps the result; set_physiology
  // refuses drivers under which the leaf could leave the model's range
  // (check_leaf_temperature_range), so every E the solve can reach is inside it.
  double leaf_temp_from_E(double E, double* dT_dE = nullptr) const;
  // Refuses a leaf whose temperature can run, over the fluxes the solve can
  // reach, from `coldest` to `hottest` deg C, if either end is past what the
  // model can evaluate: R_d's Q10 not positive, or the Tetens pole.
  void check_leaf_temperature_range(double coldest, double hottest) const;

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;

  // supply-side transpiration for a given water potential gradient between leaves and soil, 
  // references setup_transpiraiton for values (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  // Stops on psi_stem < psi_upstream, a flow from the leaf to the soil.
  double transpiration(double psi_stem, double psi_upstream);
  // The flux an OPERATING-POINT evaluation uses: zero where psi_upstream >=
  // psi_stem, the same no-flow rule set_leaf_states_rates_from_psi_stem applies,
  // and transpiration() otherwise. A cost or lambda that re-reads the flux after
  // that call must use this, or it refuses a point the evaluator accepted: a
  // prescribed collar at the wet bound puts the stem up to ~1e-9
  // MPa wetter.
  double operating_point_flux(double psi_stem, double psi_upstream);
  // supply-side transpiration for a given water potential gradient between leaves and soil, integrated internally (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration_full_integration(double psi_stem, double psi_upstream);                    
  // stomatal conductance rate of c02 (return: mol CO2 m^-2 s^-1)
  double stom_cond_CO2(double psi_stem, double psi_upstream); // define as a constant
  // converts transpiration in kg h20 s^-1 m^-2 LA to psi_stem (return: -MPa)
  double transpiration_to_psi_stem(double transpiration_, double psi_upstream);
  
  // assimilation functions

  //
  double assim_rubisco_limited(double ci_);
  double electron_transport();
  double assim_electron_limited(double ci_);
  double assim_colimited(double ci_);

  // Photosynthesis at this leaf's temperature and light, as the assimilation
  // kernel reads it: the two limiting capacities, the curvature that colimits
  // them, and the dark respiration. Nothing in it moves with intercellular CO2,
  // which is what makes dA/dci a slope in FOUR directions rather than in one per
  // pack slot -- and the chain from these back to the traits is the tape's,
  // because they are recorded before the slope is taken.
  template <class T>
  struct PhotoCapacity {
    T vcmax;
    T transport;
    T curvature;
    T respiration;

    // ⚠️ WHAT CARRIES A ROW, NAMED FOR visit_active -- see SupplyDraw. A shape
    // dispatch does not open is skipped in silence, and the ci residual reads
    // every one of these. The member count is asserted below this struct.
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

  // The coordinates at an interior collar, held, as operating_collar placed them: the
  // capacity, the stem potential sigma, the stem's flux and ci, each carrying
  // its rows in the traits and the soil. billed_outputs_at reads them rather than
  // solving sigma and ci again. The caller owns them, as it owns the draw.
  //
  // ⚠️ THEY LIVE ON THE CALLER'S TAPE, SO THEIR LIFETIME IS THE DRAW'S. They are
  // recorded before operating_collar's residual is rewound, which is why they survive
  // it; a tape reset or a new recording between operating_collar and billed_outputs_at
  // leaves them pointing at slots that no longer hold them.
  //
  // `at` is the collar they were placed at, NaN until operating_collar places them.
  // operating_collar places them only at an interior point and only at an active
  // scalar, so at double and at every other kind billed_outputs_at solves its own.
  template <class S>
  struct HeldCoords {
    double at = std::numeric_limits<double>::quiet_NaN();
    PhotoCapacity<S> cap{};
    S sigma{};
    S gc{};
    S ci{};

    bool placed() const { return !std::isnan(at); }

    // ⚠️ WHAT CARRIES A ROW, NAMED FOR visit_active -- see SupplyDraw. The
    // marginal's residual reads every one of these, and one left out is a
    // column of exact zeros in the collar's row.
    template <class F>
    void for_each_active(F&& f) const {
      cap.for_each_active(f);
      f(sigma);
      f(gc);
      f(ci);
    }
  };

  // --- the scalar-generic cores -------------------------------------------
  //
  // These are the SINGLE definition of the photosynthesis and cost algebra. The
  // five `double` entry points above and beside them are `T = double`
  // instantiations, and `dprofit_droot_collar_psi` differentiates the same code
  // with `T = xad::fwd<double>::active_type`. So the forward model and its
  // derivative cannot disagree: there is one body, not two.
  //
  // ⚠️ **They replace a hand-maintained `namespace detail`, which HAD ALREADY
  // DRIFTED** -- the deleted `assim_colimited_ad` associated the electron-limited
  // term left-to-right where `assim_electron_limited` divides the bracket first,
  // and used `s*s` where `assim_colimited` uses `pow(s, 2)`. Both are pure
  // reassociation, so the two functions were mathematically identical and
  // numerically not: the AD derivative was the derivative of a *slightly
  // different function* than the model evaluated. Harmless while the gradient only
  // set TF24f's acclimation rate; load-bearing because the collar solve
  // root-find on it. Measured cost of the fix: 4.98e-07 on the golden grid.
  //
  // ⚠️ A KERNEL TAKES THE SCALARS THAT MOVE, AS ARGUMENTS. That is what makes the
  // read set checkable: the parameter list IS the declaration of what the kernel
  // differentiates against, so a row cannot go missing because someone maintained
  // a list wrongly, and a slope through the kernel walks as many directions as it
  // has arguments rather than as many as the pack is wide. Fixed physiology that
  // no trait reaches -- `gamma_`, `km_`, the atmosphere -- stays a member, because
  // nothing differentiates it.
  //
  // Keep them PURE -- no writes to members. `hydraulic_cost_TF` caches into
  // `hydraulic_cost_`; its kernel must not, or the AD pass would write model state
  // while probing.
  template <typename T, typename U> T assim_rubisco_limited_kernel(const T& ci, const U& vcmax) const;
  template <typename T, typename U> T assim_electron_limited_kernel(const T& ci, const U& transport) const;
  template <typename T, typename U> T assim_colimited_kernel(const T& ci,
                                                 const PhotoCapacity<U>& cap) const;
  // The curve parameters take their own scalar U: the forward model and its
  // tangents pass the members as plain doubles, and pow(tangent, double) needs
  // no log where pow(tangent, tangent) differentiates the exponent too.
  template <typename T, typename U> T proportion_of_conductivity_kernel(
      const T& psi, const U& b, const U& c) const;
  template <typename T, typename U> T hydraulic_cost_TF_kernel(
      const T& psi_stem, const U& scale, const U& b, const U& c,
      const U& beta2) const;

  // The same kernels reading the members, for the forward model, which has one
  // set of these and no reason to thread them through every call.
  template <typename T> T assim_rubisco_limited_kernel(T ci) const;
  template <typename T> T assim_electron_limited_kernel(T ci) const;
  template <typename T> T assim_colimited_kernel(T ci) const;
  template <typename T> T proportion_of_conductivity_kernel(T psi) const;
  template <typename T> T hydraulic_cost_TF_kernel(T psi_stem) const;

  // The kernels' arguments, read off a caller's parameter pack, for the path
  // that needs rows in the traits.
  //
  // ⚠️ FROM THE PACK, NOT THE MEMBERS. At double the two are the same numbers, so
  // a member read answers correctly and differentiates to a silent ZERO in every
  // trait it touches -- which is a plausible gradient, in every column at once.
  // The atmospheric constants and the fixed physiology stay members, because
  // nothing differentiates them.
  //
  // The temperature responses are re-run rather than read off vcmax_ / jmax_ /
  // R_d_, and that is what carries a trait's row from its reference value at 25 C
  // to the value at this leaf's temperature. Each is exactly linear in that
  // reference value, so the response contributes a factor and the trait carries
  // the row.
  template <typename T> T vcmax_at(const leaf_pars<T>& pars) const;
  template <typename T> T jmax_at(const leaf_pars<T>& pars) const;
  template <typename T> T respiration_at(const leaf_pars<T>& pars) const;
  template <typename T> T electron_transport_at(const leaf_pars<T>& pars) const;

  // The one place that says which pack slots each kernel reads.
  template <typename T> PhotoCapacity<T> photo_capacity_at(const leaf_pars<T>& pars) const;
  template <typename T> PhotoCapacity<T> photo_capacity() const;
  // Net assimilation as c_i -> infinity on the capacity `cap`: the colimited
  // capacity less respiration, which is what a leaf with no conductance nets
  // where it cannot reach zero at any c_i (the g_min = 0 limit of a shut leaf
  // in the dark). At any scalar, so the tape differentiates it too.
  template <typename T, typename U>
  T assim_saturated_kernel(const PhotoCapacity<U>& cap) const;
  template <typename T>
  T hydraulic_cost_TF_kernel(const T& psi_stem, const leaf_pars<T>& pars) const;

  // The slope of a kernel in its first argument, carrying rows.
  //
  // A slope like dA/dci is SECOND order in the traits, and taking it at the
  // caller's scalar means nesting a tangent above an adjoint -- which
  // `odelia/tangent.hpp` forbids by name, and which measures 521 tape statements
  // against the 17 of the kernel it differentiates. Taken here through the
  // kernel alone, over a tangent of a tangent at double, and handed over as
  // supplied rows: nothing is recorded, and the slope still cannot disagree with
  // the function because it IS a tangent through it.
  //
  // ⚠️ ONE DIRECTION PER ARGUMENT, AND THE ARGUMENTS ARE THE WHOLE READ SET. The
  // walk costs the kernel's arithmetic times the direction count, so handing it
  // a twenty-slot pack to read four numbers out of prices seventeen columns of
  // structural zero: measured, the pack walk was 34.7% of the boundary's
  // instructions and `xad::FReal`'s constructor alone was 18.9% of the run. Every
  // argument's row is supplied whether or not it came back non-zero, because a
  // slot the kernel cannot read is a slot it cannot be handed.
  //
  // ⚠️ THE ROW AGAINST `x` IS AS REAL AS THE OTHERS. The slope is read at a point
  // that MOVES with the traits, so the kernel's second derivative in its own
  // argument is a first-order channel exactly like every argument's. Left out, it
  // is an exact zero in a column rather than an error.
  //
  // `kernel` must DECLARE its return type (`-> U`): a deduced one hands back an
  // XAD expression template over temporaries that are dead on return.
  template <class S, class Kernel, std::same_as<S>... Args>
  S kernel_slope_at(Kernel&& kernel, const char* what, const S& x,
                    const Args&... args) const;

  // dA/dci and dC/dsigma, each a tangent through its own kernel rather than a
  // slope written out beside it: a hand-written twin is a second definition of
  // one function, free to disagree with it by an amount that stays finite and
  // plausible in every column.
  template <typename T>
  T assim_slope_at(const T& ci, const PhotoCapacity<T>& cap) const;
  template <typename T>
  T cost_slope_at(const T& sigma, const leaf_pars<T>& pars) const;

  // The STEM's flux at a held stem potential and collar: kmax * (G(sigma) -
  // G(collar)), which is upstream's `transpiration()` given its parameters.
  //
  // ⚠️ THIS IS NOT THE SOIL DRAW, and the difference is the whole reason ci has
  // to read this one. T1 says the stem carries what the soil delivers, but the
  // collar solve satisfies it only to its own tolerance -- measured 6.7e-04
  // apart at a converged interior point -- and the model places ci from
  // stom_cond_CO2, which is built from THIS flux. Differentiating a residual
  // the model did not solve returns a number rather than an error, and the
  // number was 3 to 10 per cent wrong in dci/dcollar with every value finite.
  template <typename T>
  T transpiration_at(const T& sigma, const T& collar,
                     const leaf_pars<T>& pars) const;

  // G(psi) at any scalar: the table's value at the passive point, carrying the
  // query slope and the two TRAIT rows the closed form gives exactly.
  template <class S>
  S stem_integral_at(const S& psi, const leaf_pars<S>& pars) const;

  // The two conditions that place the operating point at a held collar, and both
  // halves of each coordinate.
  //
  //   T1  kmax * (G(sigma) - G(collar)) - E_up = 0   places the stem potential
  //   T2  A(ci) * k1 - gc * (ca - ci) * k2   = 0     places the intercellular CO2
  //
  // gc is proportional to the flux, so it is derived from the draw here rather
  // than passed: two spellings of one fact is a place they can disagree. The DRAW
  // rather than a bare flux, so a flux and a slope from different collars cannot
  // be paired -- check_draw refuses that before anything reads it.
  //
  // `collar_moves` says whether profit carries the collar's channel here. At an
  // interior point it does not -- billed_outputs_at evaluates profit at a HELD collar,
  // which is the envelope omission -- and the slopes are then absent rather than
  // zero. dprofit_dcollar_at, the only reader of a slope, passes true.
  template <class S>
  CollarCoords<S> collar_coords_at(double sigma_star, double ci_star,
                                   const S& collar, const SupplyDraw<S>& draw,
                                   const leaf_pars<S>& pars,
                                   const PhotoCapacity<S>& cap,
                                   bool collar_moves) const;

  // collar_coords_at in its two halves, so a caller holding the first can add
  // the second without solving sigma and ci again: the coordinates at the HELD
  // collar, carrying the trait and soil rows, and then the collar's channel on
  // top of them. collar_coords_at is exactly the one followed by the other.
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
  // (seat_stem_at_dry_pin). There sigma is that end as a function of the traits,
  // and it does not move with the collar: both halves above read this.
  bool stem_held_at_dry_pin(double sigma_star) const {
    return operating_point_kind_ == OperatingPointKind::BoundaryDry &&
           sigma_star == vulnerability_psi_end(stem_b, stem_c);
  }

  // The stem's critical potential at any scalar. DERIVED from the trait pair, as
  // it is at double, which is what makes the old `no_gradient` question moot: it
  // is an ordinary expression the chain differentiates, not a parameter that has
  // to be declared not to move.
  template <class S>
  S psi_crit_at(const leaf_pars<S>& pars) const;
  // dry_end<K>() at any scalar: psi_crit_at for a curve bounded there, and the
  // end of the stem curve, b x_end^(1/c), for every other.
  template <CostCurve K, class S>
  S dry_end_at(const leaf_pars<S>& pars) const;
  template <class S>
  S stem_curve_end_at(const leaf_pars<S>& pars) const;

  // The objective at a placed operating point, and everything plant reads.
  //
  // ⚠️ TEMPLATED ON THE CURVE, IMPLEMENTED FOR TF24 AND TF24_floor. A
  // static_assert refuses the other seven, so another curve is a branch here rather
  // than a redesign -- and until something differentiates one, a branch written
  // for it would have no caller and no reference to check it against.
  template <CostCurve K, class S>
  S profit_at(const S& collar, const SupplyDraw<S>& draw,
              const leaf_pars<S>& pars, bool collar_moves) const;
  // The objective at coordinates already placed, which is all profit_at does
  // once it has them.
  template <CostCurve K, class S>
  S profit_from_coords(const S& collar, const leaf_pars<S>& pars,
                       const PhotoCapacity<S>& cap, const S& sigma,
                       const S& ci) const;

  // `placed`, when operating_collar placed it at this point, is read instead of
  // solving sigma and ci again; otherwise, and when it is null, this solves its
  // own. ⚠️ Placed coordinates are refused at any kind but interior and at any
  // collar but the one they were placed at, because either means the leaf has
  // moved.
  template <CostCurve K, class S>
  BilledOutputs<S> billed_outputs_at(const S& collar, const SupplyDraw<S>& draw,
                                     const leaf_pars<S>& pars,
                                     const HeldCoords<S>* placed = nullptr) const;

  // M = dprofit/dcollar, as a SCALAR. What used to travel beside it -- the two
  // collar responses -- are the coordinates' own slope halves now, so there is
  // nothing left to return with it.
  template <CostCurve K, class S>
  S dprofit_dcollar_at(const S& collar, const SupplyDraw<S>& draw,
                       const leaf_pars<S>& pars) const;
  // The same marginal on coordinates operating_collar has already placed, adding only
  // the collar's channel.
  template <CostCurve K, class S>
  S dprofit_dcollar_at(const S& collar, const SupplyDraw<S>& draw,
                       const leaf_pars<S>& pars, const HeldCoords<S>& held) const;
  template <CostCurve K, class S>
  S dprofit_dcollar_from_coords(const S& collar, const leaf_pars<S>& pars,
                                const PhotoCapacity<S>& cap,
                                const CollarCoords<S>& at) const;

  // dM/dp at the operating point the solve placed -- the curvature the interior
  // derivation divides by. A FIRST difference of the marginal, which is what this
  // has always meant: the failure it replaces was a SECOND difference of the
  // OBJECTIVE, whose second-order content was two hand-written Taylor
  // coefficients and which put a POSITIVE curvature at an interior maximum on a
  // century stand.
  //
  // ⚠️ NOT A TANGENT THROUGH dprofit_dcollar_at, AND THE REASON IS CONDITIONING RATHER
  // THAN EFFORT. dci/dcollar is a cancellation of two terms of 8.1e+04 that
  // leaves 0.64, and its collar derivative is dominated by dci/dsigma * dV/dp
  // where dV/dp is ITSELF a cancellation -- 0.50 against 0.50, leaving -1.1e-04.
  // So an assembly has to hold four significant digits through that, and then
  // multiply by 8.1e+04. It cannot: the conductivity's VALUE is the table's and
  // its response is the curve's, and those disagree by the fit error, which is
  // larger than the residual being cancelled to. A tangent returns an eighth of
  // the answer, with the right sign.
  //
  // The difference does not have that problem, because M itself is O(1) away
  // from the optimum and is not formed by cancelling anything. Measured over
  // three decades of step it is stable to SEVEN significant figures.
  //
  // ⚠️ NOT const: it drives the model to two neighbouring collars and restores
  // the operating point afterwards. Leaving it at p0 +/- h would break the rule
  // that every path out of the solve writes its own rates, and this is a path out.
  template <CostCurve K>
  double profit_curvature_at_collar();

  // The collar this solve left, at whatever scalar the caller wants, from
  // whatever pins it. The kind chooses which condition closes the system and
  // nothing else.
  //
  // A caller that has already taken profit_curvature_at_collar() -- to refuse on a
  // curvature it cannot divide by -- passes it as `interior_curvature` rather
  // than letting this take it again: it is the same call with the same argument,
  // and it costs two model evaluations.
  // NOT const. The interior arm takes the curvature by differencing, which
  // drives the model off the operating point and puts it back -- so this mutates
  // the leaf even though it leaves it where it found it. Declared const it
  // const_cast itself, which hid that from every caller.
  //
  // `placed`, when given, receives the coordinates the interior arm solved,
  // for billed_outputs_at to read instead of solving them again. It is left unplaced
  // at double and at every other kind, which have none to hand over.
  template <CostCurve K, class S>
  S operating_collar(const SupplyDraw<S>& draw, const leaf_pars<S>& pars,
                     double interior_curvature =
                         std::numeric_limits<double>::quiet_NaN(),
                     HeldCoords<S>* placed = nullptr);

  // The smallest |dM/dp| operating_collar's interior arm divides by. Every collar row
  // there is a quotient over the curvature, so one near zero returns
  // amplification rather than a derivative -- finite, and huge. The value is
  // plant's default gradient_curvature_floor, so plant's own guard and this one
  // refuse the same points at its defaults; the smallest |curvature| on
  // test_transpose's grid is 0.14.
  static constexpr double interior_curvature_floor = 1e-3;

  // The conditions that pin a collar, one function per bound. Neither is a
  // second derivative, so each is an ordinary implicit value: the residual's own
  // slope at the bound, and the residual at S for the rest.
  //
  // The wet bound: the collar at which uptake meets the residual flux (at
  // g_min = 0, where uptake is zero). It reads no cost, so it takes no curve.
  template <class S>
  S wet_bound_at(double bound_x, const SupplyDraw<S>& draw,
                 const leaf_pars<S>& pars) const;
  // The dry end: the continuity root with the stem held at dry_end_at<K>(),
  // which is where the curve enters.
  template <CostCurve K, class S>
  S dry_end_collar_at(double bound_x, const SupplyDraw<S>& draw,
                      const leaf_pars<S>& pars) const;


  // The pack as this leaf currently stands. The double path passes this where the
  // active path passes a copy carrying directions, so both reach the kernels the same way and
  // there is one call convention rather than a member-reading twin.
  leaf_pars<double> passive_pars() const;
  double assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream);
  double psi_stem_to_ci(double psi_stem, double psi_upstream);
  void set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream);


// --- Marginal cost of water ------------------------------------------------
// lambda = dA/dE, the marginal carbon gain per unit water lost. This is the
// quantity that unifies the stomatal optimality models: they all maximise a
// profit, so they all satisfy dA/dE = lambda at the optimum and differ only in
// lambda(state). Reporting it makes this model directly comparable to the
// others, and to fitted Medlyn g1 values via g1_eff() below.
//
// Deliberately accessors rather than stored state: they are wanted once per
// solve, not once per candidate potential, and computing them inside
// set_leaf_states_rates_from_psi_stem would put a pow() on a path that runs
// ~10^3 times per solve.
//
// NOTE the distinction to be careful of: the member `CF77_lambda_` is an *input*, the
// Cowan-Farquhar marginal value of water, never set by set_physiology. The lambda
// here is an *emergent output* -- the marginal cost the operating point implies.

  // lambda at an arbitrary stem water potential (positive magnitude, MPa), in
  // umol CO2 (kg H2O)^-1. Analytic, from the TF24 cost function:
  //
  //   C(psi)   = TF24_cost_scale * (1 - f)^TF24_beta2,   f(psi) = exp(-(psi/stem_b)^stem_c)
  //   dC/dpsi  = TF24_cost_scale * TF24_beta2 * (1-f)^(TF24_beta2-1) *
  //              (stem_c/stem_b)(psi/stem_b)^(stem_c-1) * f
  //   E(psi)     = kmax * integral of f,  so  dE/dpsi = kmax * f
  //   lambda     = (dC/dpsi) / (dE/dpsi)
  //
  // The f cancels exactly; it is cancelled here rather than divided out, so the
  // expression stays finite as psi -> psi_crit where f -> 0.
  //
  // Caveat: with TF24_beta2 < 1 the (1-f)^(TF24_beta2-1) factor diverges as psi -> 0
  // (f -> 1), which is a property of the cost function, not of this code.
  double lambda_TF24(double psi_stem) const;
  // ⚠️ TAKES BOTH POTENTIALS, where lambda_TF24 takes one. That is not a
  // convenience: JS22's cost is a function of the DROP, so its marginal cost
  // genuinely depends on the upstream potential, and no `f` cancels out of the
  // ratio the way it does for TF24.
  double lambda_JS22(double psi_stem, double psi_upstream) const;
  // Also both potentials, but only because `E` is measured from the upstream one --
  // the COST here is a function of the absolute potential, unlike JS22's.
  double lambda_CMax(double psi_stem, double psi_upstream) const;
  // ⚠️ THE ONE lambda IN THIS CLASS THAT DOES NOT VANISH AT THE WET END, and that
  // is the whole content of the TF24_floor curve. It is `lambda_TF24` plus a constant:
  //
  //     lambda = lambda_TF24(psi) + lambda_o
  //
  // so as the drop closes the first term goes to zero and the price does not.
  // Every conductance-loss curve here has `lambda -> 0` there by construction,
  // which prices water as free precisely when it is abundant; this is the repair.
  //
  // ⚠️ ONE ARGUMENT, like `lambda_TF24` and unlike `lambda_JS22`, because the
  // hydraulic half reads the ABSOLUTE potential. `lambda_for<TF24_floor>` discards
  // the upstream potential accordingly.
  double lambda_TF24_floor(double psi_stem) const;

  // --- the product-objective family --------------------------------
  //
  // Eller (2018, 2020)'s SOX maximises `A * g(psi)` rather than `A - C(psi)`, with
  // `g` the conductivity fraction rescaled onto [0, 1]. Sabot's `TractLSM` calls it
  // `kcost` and Jones et al. (2026) is the same objective with a different `g`.
  //
  // ⚠️ IT TAKES NO NEW PARAMETER. `g` is built from `proportion_of_conductivity`
  // and `k_crit_fraction`, both of which the object already has, which is why this
  // curve is cheaper than every subtracted one.
  double sox_reduction(double psi_stem) const;
  double sox_reduction_deriv(double psi_stem) const;
  double profit_psi_stem_SOX(double psi_stem, double psi_upstream);
  double lambda_SOX(double psi_stem, double psi_upstream) const;

  // Jones et al. (2026): the same product objective with a LINEAR reduction factor,
  // `1 - psi/psi_crit`, where SOX's follows the vulnerability curve. Sabot's
  // `TractLSM` has this one too, as `phiLWP`, crediting Dewar.
  //
  // ⚠️ ALSO TAKES NO PARAMETER, because `psi_crit` is DERIVED here -- the 95%
  // quantile of the stem curve, from `stem_P50` and `stem_c`. The paper supplies its
  // own `Pcrit`; deriving it instead is what makes this and SOX two interpolations
  // between the SAME two anchors (1 at zero tension, 0 at `psi_crit`) and so
  // comparable by construction rather than by matching a parameter.
  // --- least-cost / cost-per-unit-carbon --------------------------------------
  //
  // Prentice et al. (2014) and Dong et al. (2017) minimise the cost of acquiring a
  // unit of carbon, `C_P = (c_E*E + c_V*Vcmax)/A`, rather than maximising a
  // difference. Minimising a ratio is maximising its reciprocal, and `log` turns
  // that into `log A - log(c_E*E + c_V*Vcmax)` -- the same benefit link SOX and
  // JW26 use, with the penalty now carrying CAPACITY as well as flux. That is the
  // reduction written out in S1.4 of the theory manuscript; the family form is
  // `h = log`, `Theta = log(c_E*E + c_V*Vcmax)`.
  //
  // ⚠️ THE ONLY CURVE HERE WHOSE PENALTY CARRIES SOMETHING OTHER THAN WATER. Every
  // other cost in this class is a function of psi (or of E) alone; this one weighs
  // photosynthetic capacity against transpiration, which is why its first-order
  // condition reproduces Prentice et al.'s `chi = xi/(xi + sqrt(D))` and why a
  // light sweep on it measures acclimation rather than the price of water.

  // The weighted resource use `W = E_mol + beta*Vcmax_mol`, mol m^-2 s^-1, and the
  // ONE place the two unit conversions live -- E from kg to mol, Vcmax from umol to
  // mol. The cost is `log W`, its derivative `W'/W`, the objective `A/W` and the
  // price `A*kg_to_mol_h2o/W`, so all four read this and cannot disagree about the
  // units beta is dimensionless in.
  double leastcost_flux(double psi_stem, double psi_upstream);

  // ⚠️ A RATIO, so what comes back is NOT in carbon units -- the same caveat
  // `profit_psi_stem_SOX` carries, and for a related reason. It is `A/W`, umol C
  // per mol of weighted resource, and it is `c_E/C_P`: the reciprocal of the
  // published cost per unit carbon, times a constant this curve does not carry
  // (see `LeastCost_beta`). Monotone in `-C_P`, so the argmax is the published
  // model's, and the VALUE must not be compared with a TF24 profit.
  double profit_psi_stem_LeastCost(double psi_stem, double psi_upstream);

  // The marginal cost of water this curve implies, umol C (kg H2O)^-1, on the same
  // axis as every other `lambda_*` here:
  //
  //     lambda = h'(A)^-1 * dTheta/dE = A * kg_to_mol_h2o / W
  //
  // at a fixed leaf temperature. With the energy balance on, Vcmax moves with E
  // through Tleaf and dTheta/dE gains `beta*umol_to_mol*dVcmax/dT*dT/dE / W`.
  //
  // which is eq. S1.4 of the theory manuscript, `lambda = c_E*A/(c_E*E +
  // c_V*Vcmax)`, with `c_E` cancelled and the kg-to-mol factor supplied.
  //
  // ⚠️ IT DOES NOT VANISH AT THE WET END, and it is the only LOG-link curve here of
  // which that is true. As E -> 0 this tends to `A*kg_to_mol_h2o/(beta*Vcmax_mol)`,
  // finite and set by the capacity term; SOX's and JW26's both carry `A` to zero
  // through their own `g`. It does carry `A`, so it still falls as assimilation
  // does -- the wet end is finite in psi, not in light.
  double lambda_LeastCost(double psi_stem, double psi_upstream);

  double jw26_reduction(double psi_stem) const;
  double jw26_reduction_deriv(double psi_stem) const;
  double profit_psi_stem_JW26(double psi_stem, double psi_upstream);
  double lambda_JW26(double psi_stem, double psi_upstream) const;

  // The marginal cost of water the ProfitMax cost implies at `psi_stem`, in the
  // same umol C (kg H2O)^-1 as lambda_TF24 -- so the two are comparable.
  //
  // The normalised objective is dimensionless, so its own dC/dE is not in carbon
  // units; multiplying by |A|max restores them, which is the same scaling that
  // makes the objective a constant-lambda one. Requires prepare_profitmax() to
  // have seeded the normalisers.
  double lambda_ProfitMax(double psi_stem);

  // lambda at the current operating point. This is the single-layer value: it is
  // the marginal cost seen by the stem, and equals dA/dE only when the collar
  // potential is fixed (i.e. for the optimise_psi_stem_* solvers).
  double marginal_cost_water() const;

  // The same in molar units, mol CO2 (mol H2O)^-1, which is the convention the
  // optimality literature states lambda in.
  double marginal_cost_water_molar() const;

  // lambda including the root-network series resistance, which is what
  // find_root_collar_psi actually equalises because it optimises over the collar
  // potential with the stem following from continuity:
  //
  //   lambda_multi = lambda_TF24 * [1 + kmax*f(psi_r)/S],    S = dE_up/dpsi_r
  //
  // S comes from dE_from_soil_dpsi_collar. The bracket is >= 1, so the
  // single-layer lambda always UNDERSTATES the true marginal cost. Returns the
  // NA sentinel where S is unavailable (the branch kinks that
  // dE_from_soil_dpsi_collar reports as NaN) or zero.
  //
  // Requires that a collar solve has run, so that the supply path's per-solve
  // caches are current.
  double marginal_cost_water_multilayer();

  // Equivalent Medlyn USO slope implied by the operating point, in kPa^0.5.
  // Defined operationally from the solved chi = ci/ca, by inverting the USO
  // relation chi = g1/(g1 + sqrt(D)):
  //
  //   g1_eff = chi * sqrt(D) / (1 - chi)
  //
  // Operational rather than predicted-from-lambda on purpose: this form is exact
  // by construction and free of unit conventions, and it is directly comparable
  // to g1 values fitted by plantecophys::fitBB or tabulated by Lin et al. (2015).
  // Reconciling it against the theoretical sqrt(3*Gstar*P/(1.6*lambda)) is the
  // companion manuscript's job, not this header's.
  double g1_eff() const;

  // THE PART OF THE COST THAT IS A PRICE RATHER THAN A REALISED CARBON LOSS, in
  // umol C m^-2 s^-1. `TF24_floor_lambda_o * E` on that curve and exactly zero on
  // every other one.
  //
  // ⚠️ ZERO IS NOT A CLAIM THAT THE CURVE'S COST IS ALL REALISED CARBON. It says
  // the curve does not separate the two. CF77 is the case that matters: its whole
  // cost is `lambda * E`, and reading that as a shadow price is a perfectly
  // ordinary thing to do -- but the model supplies ONE number and nothing to
  // attribute it with, so it behaves identically whether that number is read as
  // carbon lost to water already taken or as the value of water withheld for
  // later. Reporting it as all-shadow would ship one reading as a fact; reporting
  // it as zero-shadow ships the other. Zero plus this comment is the honest
  // choice, because the field is defined only where the CURVE defines it.
  //
  // That is also the argument for the two-term form: it is the first curve here in
  // which the question "is lambda a cost to deduct, or a price that only shapes
  // behaviour?" can be posed at all.
  //
  // Reads `transpiration_` rather than recomputing it, so
  // `shadow_cost == TF24_floor_lambda_o * E` holds bit-exactly against the
  // reported `E` -- which is the identity a caller checks.
  //
  // An accessor, not stored state (hazard 5), and it deliberately does NOT touch
  // `hydraulic_cost_`: that member is the cost the OBJECTIVE subtracted, its
  // expression carries the two bit-exact reductions, and splitting it would
  // reassociate the arithmetic those depend on.
  double shadow_cost() const;

  // Stomatal conductance to WATER VAPOUR, mol H2O m^-2 s^-1 (#56). The model solves
  // for conductance to CO2; every data source and the g1 literature record it to
  // water, so the conversion belongs here rather than in each caller.
  //
  // An ACCESSOR, not a member, per hazard 5: nothing in the solve reads it, so it
  // does not earn storage.
  double stom_cond_H2O() const {
    return stom_cond_CO2_ * H2O_CO2_stom_diff_ratio_;
  }
  // ⚠️ ON THE ENERGY-BALANCE PATH stom_cond_H2O() IS THE TOTAL LEAF-TO-AIR
  // CONDUCTANCE, not a stomatal one. It is E*P/D_leaf, and E crosses the leaf
  // boundary layer as well as the stoma, so comparing it with a measured gs
  // compares different quantities. The two accessors below separate them (#122).
  //
  // Boundary-layer conductance to WATER VAPOUR, mol H2O m^-2 s^-1: the molar form
  // of the heat resistance ra_ the energy balance reads, at air temperature,
  // scaled by bl_vapour_heat_ratio. +Inf off the energy-balance path, where the
  // leaf is at air temperature and coupled to it perfectly, and where ra_ is a
  // fallback nothing validates.
  double boundary_layer_cond_H2O() const;
  // Stomatal conductance to water vapour with the boundary layer taken out in
  // series, 1/g_s = 1/g_total - 1/g_b, mol H2O m^-2 s^-1. Exactly
  // stom_cond_H2O() off the energy-balance path.
  //
  // ⚠️ REFUSES where g_total >= g_b. Then E is more than the boundary layer can
  // carry at this deficit whatever the stoma does, and no finite g_s exists.
  // Hydraulics set E here, so nothing in the solve prevents it; the refusal sits
  // in this accessor rather than in the solve so that it cannot abort a plant
  // run, and test_boundary_layer_bounds_the_total_conductance measures how far
  // a realistic driver grid stays from it.
  double stom_cond_stomatal_H2O() const;

  // Read-only on purpose: a settable one would be a way to disagree with the
  // operating point it is derived from.
  double lambda_emergent() const { return lambda_emergent_; }
  // ProfitMax's normalisers, so a caller can form its equivalence ratio without
  // this class pretending that ratio is a lambda.
  double profitmax_A_max() const { return profitmax_A_max_; }
  double profitmax_k_span() const { return profitmax_k_span_; }

  // Every reported output of a solved point, in one call.
  //
  // ⚠️ THIS EXISTS FOR THE R BOUNDARY AND FOR NOTHING ELSE. A C++ caller should
  // read the members; they are public and free. From R each read is a separate
  // call through an RcppR6 active binding at ~1.1 us, so the thirteen of them cost
  // ~15 us against a ~3 us solve -- five times the model, to report it (#39).
  // One call is ~1.5 us. That is the whole justification, and it is why this is
  // a flat vector rather than a struct: it crosses the boundary as a numeric
  // vector with no glue.
  //
  // The ORDER is the interface. R's .operating_point_names names these positions
  // and test-surface.R checks the two agree by reading all thirteen fields
  // individually and comparing -- so a field inserted here without a
  // corresponding R change fails rather than silently shifting a column.
  //
  // ⚠️ APPEND, NEVER INSERT. `Tleaf` was added at the END for that reason, even
  // though it reads more naturally beside the other state variables: R names
  // these by POSITION, and anything already reading position 6 as `profit`
  // includes the golden comparison and every consumer's saved output.
  //
  // Index 9 is the total soil uptake, summed over FINITE layers only: before a
  // solve the layers hold the NA sentinel, and R's own reader did the same
  // filtering. Everything else is read as-is, sentinels included, because
  // "unsolved" is a state the caller is entitled to see.
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
            Tleaf_,          // deg C -- APPENDED, see below
            shadow_cost(),   // APPENDED likewise, for the same reason
            // THE SEATED CURVE'S marginal cost of water, as against `lambda`
            // above, which is TF24's whatever ran. Appended rather than placed
            // beside its namesake, because these names are positions.
            lambda_emergent(),
            // The residual (cuticular) flux the leaf realised, kg H2O m^-2
            // s^-1 (#182). APPENDED: `E` above stays the TOTAL, and the
            // stomatal flux is their difference.
            residual_flux_};
  }

// leaf economics functions
  double hydraulic_cost_TF(double psi_stem);

  // Cowan & Farquhar (1977): the cost of water is the water itself, priced at a
  // constant marginal value. `CF77_lambda_ * E`, in umol C m^-2 s^-1 -- a carbon flux,
  // because `CF77_lambda_` is carbon per unit transpiration and E is a mass flux.
  //
  // This is the model `CF77_lambda_`'s units belong to, and the reason its first-order
  // condition is the one the whole optimality literature is written in:
  // maximising `A - lambda*E` over psi_stem gives dA/dE == lambda at the optimum.
  // `marginal_cost_water()` reports the same quantity for the OTHER cost curves,
  // which is what puts them all on one axis.
  double hydraulic_cost_CF77(double psi_stem, double psi_upstream);

  // THE CF77 PRICE, and the only place the soil-beta option is applied. All three
  // CF77 arms read it -- the cost, its derivative, and lambda_for -- so the value,
  // the marginal value, and the reported lambda cannot disagree about whether the
  // option is on. Same units as `CF77_lambda_`: umol C (kg H2O)^-1.
  double cf77_price() const;

  // Joshi & Stocker (2022)'s hydraulic term, `gamma*(dpsi)^2`, quadratic in the
  // soil-to-leaf DROP where the other two curves are functions of the absolute
  // potential. See `JS22_gamma` for what is missing from it.
  double hydraulic_cost_JS22(double psi_stem, double psi_upstream);

  double profit_psi_stem_TF(double psi_stem, double psi_upstream);
  double profit_psi_stem_CF77(double psi_stem, double psi_upstream);
  double profit_psi_stem_JS22(double psi_stem, double psi_upstream);

  // Wolf/Anderegg CMax: the antiderivative of `CMax_a*psi + CMax_b` taken from the
  // upstream potential, so the cost of moving no water is zero.
  double hydraulic_cost_CMax(double psi_stem, double psi_upstream);
  double profit_psi_stem_CMax(double psi_stem, double psi_upstream);

  // The TF24_floor cost: TF24's `gamma*(1-f)^beta2` plus CF77's `lambda_o*E`, in
  // carbon units like both of them and so written into `hydraulic_cost_` on the
  // same footing.
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
  // temperature k_i_Tiso_ moved by delta_impair degrees (refused outside
  // [250, 350] K).
  // Its own parameters: it reads nothing of the Vcmax or Jmax curves.
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
  // d(state)/d(input) along the trajectory, exact by forward-mode differentiation
  // of the same RK4: 3 * times.size() * which.size() numbers, time fastest, then
  // state (A, I, D), then input. `which` indexes damage_par_names().
  std::vector<double> damage_trajectory_gradient(
      const std::vector<double>& times, const std::vector<double>& leaf_temp,
      double A0, double I0, double D0, const std::vector<int>& which) const;

  // THE LAYER'S DYNAMICS, per DAY, at a leaf temperature and a state:
  // {dA/dt, dI/dt, dD/dt}. Reads the layer's parameters and nothing of the
  // members A_, I_, D_ or the drivers, and does not need set_physiology, so a
  // caller can evaluate it at a measured leaf temperature. dD/dt is maturation
  // only: a caller whose leaf area turns over (plant) subtracts that dilution
  // itself. Does not read the gate either: asking for the rates is the request.
  std::vector<double> damage_rates(double leaf_temp, double A, double I,
                                   double D) const;
  // Integrates damage_rates over a leaf-temperature series, linear between the
  // given times (days, strictly increasing), from (A0, I0, D0) at times[0]. Fixed-
  // step RK4, with each interval cut so that the fastest rate times the step stays
  // below 0.1 and the temperature moves at most 0.25 K per step. Returns the
  // states at every time, column-major (A, I, D): 3 * times.size() numbers.
  std::vector<double> damage_trajectory(const std::vector<double>& times,
                                        const std::vector<double>& leaf_temp,
                                        double A0, double I0, double D0) const;
  // Throws, naming `who`, when use_thermal_damage_ is on. Called by every path
  // that reads the capacities off the parameter pack (vcmax_at and the rest),
  // which do not apply the layer.
  void refuse_thermal_damage(const char* who) const;

  // Sperry (2017) ProfitMax, with BOTH terms normalised as the paper defines them
  // -- see the ProfitMax notes above optimise_psi_stem_single for what that buys
  // over the lambda form.
  // Seeds |A|max and the conductance span; the two below read what it seeds.
  void prepare_profitmax();

  // The conductance normalisers alone: `k(psi_soil)` and the span down to
  // `k_crit`. Both are analytic in the traits, so they are differentiable and must
  // be refreshed whenever the traits move.
  double profitmax_reference_psi();
  void prepare_profitmax_norms();

  // ⚠️ THE SAME PREPARATION WITH `|A|max` SUPPLIED RATHER THAN SCANNED, and the
  // reason it exists is a gradient.
  //
  // `prepare_profitmax()` finds `|A|max` by scanning the supply stream, so its
  // argmax is piecewise constant in the traits: the derivative is a sequence of
  // zeros and jumps, and a total derivative through it is not a quantity worth
  // computing. What IS well defined is the PARTIAL at fixed normaliser, which is
  // what this package reports and what the code already computes elsewhere.
  //
  // A gradient's perturbation loop calls `set_traits` + `set_physiology`, which
  // clears the whole profitmax block -- measured, `|A|max` goes from 16.757 to
  // NaN. So a perturbation must refresh the analytic normalisers and HOLD `|A|max`
  // at the base point. Passing it in makes that pinning visible in the signature
  // instead of a field write somebody has to notice.
  void prepare_profitmax_at(double A_max);
  double profit_psi_stem_ProfitMax(double psi_stem, double psi_upstream);

  // The whole cost/gain/profit curve over [psi_soil, psi_crit] in ONE crossing of
  // the R boundary: n rows of (psi, CG, HC, TC, profit), flattened column-major.
  // This is Sicangco et al.'s Figures 2, 3 and S4, and building it row by row from
  // R would pay ~1.8 us of call overhead against a ~3 us model evaluation.
  std::vector<double> profitmax_curve(int n);

  // Clear the outputs a single-layer optimiser does NOT write. Hazard 8: these
  // three describe a ROOT-COLLAR solve, and optimise_psi_stem_* never runs one,
  // so leaving the last find_root_collar_psi()'s values standing reports an
  // operating point whose parts came from two different solves. Measured before
  // this existed: E = 9.216e-5 from the Sperry solve sitting beside E_up =
  // 2.626e-5 from a collar solve several calls earlier.
  void clear_collar_solve_state();

  // --- WHICH KIND of operating point the collar solve found -------------------
  //
  // The operating point is not one kind of thing, and which kind it is changes
  // what a derivative taken from it means. Along one drydown the leaf passes
  // through these in order: an interior profit maximum while the soil is wet,
  // then a constrained optimum pinned against a hydraulic limit, then shutdown.
  // At an interior optimum a trait moves a FREE optimum; at a pin it moves the
  // LIMIT itself; at shutdown the water response is zero while the carbon
  // response is not. A consumer of `dprofit_droot_collar_psi` needs to know
  // which of those it is holding, and nothing in the returned numbers says.
  //
  // ⚠️ **This is recorded by the branch that was TAKEN, never inferred after the
  // fact from a residual**, and that is the whole point of it.
  // `dprofit_at_collar_psi` returns a hard 0.0 SENTINEL in the no-flow /
  // infeasible state (see its `feasible` out-parameter). So a classifier of the
  // form "|dprofit| within tolerance => interior optimum" records such a state as
  // a stationary point, and then reads a curvature of zero off the same
  // sentinel -- twice-confirmed and wrong. No residual test can separate the two,
  // because the two return the same number.
  enum class OperatingPointKind {
    // No solve has run on this object since it was constructed, re-traited, or
    // driven through one of the off-path psi_stem optimisers. Reading an
    // operating point in this state reads NA sentinels.
    Unsolved,
    // Interior profit maximum: dprofit == 0 was solved for, strictly inside the
    // feasible collar interval. 232 of the 336 golden grid points at 25 C -- 172 of
    // them at 40 C, where the optimum presses against the WET bound instead.
    Interior,
    // Constrained optimum at the psi_soil end of the feasible interval, just
    // inside root_zero_E (the collar at which uptake is exactly zero). The
    // objective is still climbing toward that bound, so the bound is the answer:
    // the leaf declines to open. 45 golden points at 25 C. The gradient here is
    // genuinely non-zero.
    //
    // ⚠️ THESE TWO ARE NAMED FOR THE BOUND THAT BINDS, NOT FOR THE CONDITIONS,
    // and an earlier `PinnedWet`/`PinnedDry` pair got read backwards for exactly
    // that reason. `psi_stem == psi_soil` is the least-tension end of the
    // bracket, so it reads as the comfortable case -- but it is what a leaf does
    // under STRESS, when no interior point earns its water. The golden grid says
    // so: 45 -> 116 rows as the leaf goes 25 -> 40 C, while the dry end goes
    // 11 -> 0. Warming drives this one and removes the other. Do not reintroduce
    // a name whose reading inverts with temperature.
    BoundarySoil,
    // Constrained optimum at the dry end: the continuity root with the stem
    // held at dry_end<K>(). For a curve whose own domain ends at psi_crit
    // (cost_domain_ends_at_psi_crit) that is psi_crit, the most tension its
    // form allows. For every other curve psi_crit bounds nothing (#183), and
    // the dry end is the end of the stem curve: profit still rising there means
    // the leaf transpires all the stem can carry, to within e^-40 -- a dry soil
    // where A climbs steeply at the little flux left. Named for the end, not
    // for psi_crit, because for most curves psi_crit is not where it is.
    BoundaryDry,
    // Shutdown on water: no collar potential both carries more than the
    // residual flux and keeps the stem short of its dry end, so there is no
    // interval to optimise over. The leaf sits at the wet bound -- the collar at
    // which uptake equals the residual flux, the stem where it carries it --
    // with its stomata shut and every output from the ordinary evaluator there.
    // For a curve bounded at psi_crit this is a soil drier than psi_crit can
    // supply; for every other curve the dry end is the end of the stem curve,
    // where it carries G(inf) to rounding, so a leaf reaches this only when the
    // soil is drier than that and its flux is zero to rounding anyway -- the
    // shutdown is asymptotic. Light and heat do not shut a leaf: with the
    // cuticle's CO2 pathway the residual has a root for every leaf (#183), and a
    // leaf in the dark is an ordinary optimum, usually at the wet bound.
    HydraulicShutdown,
    // The collar potential was IMPOSED by the caller (evaluate_root_collar_psi /
    // profit_at_collar_psi) inside the feasible interval, not optimised.
    // dprofit is generally non-zero and is the whole point on that path -- it is
    // the acclimation rate, not a residual.
    Prescribed,
    // The solver could not move, and NO PLANT IS DESCRIBED. Either an endpoint
    // admitted no feasible gradient (so the golden-section fallback ran), or the
    // gradient was non-positive at the wet end AND non-negative at the dry end,
    // which makes the interior stationary point a MINIMUM and leaves the maximum
    // at one of the two ends with nothing to say which. This is deliberately NOT
    // reported as a pin: attributing it to whichever bound is nearer would return
    // a plausible, finite, wrong answer.
    SolverRefused,
    // A partial derivative came back non-finite where feasibility said it should
    // not have. Distinct from SolverRefused because it is a different fault, and
    // there is nowhere downstream to put the test.
    NonFiniteGradient,
  };

  // This leaf's own sites plus the supply model's, which are ONE list. Summed on
  // read rather than shared on construction, so rebuilding the root network
  // cannot silently detach the tally.
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

  // How many kinds there are, counted from the LAST enumerator so a kind added
  // to the list is inside this by construction. A consumer tallying by kind
  // sizes its array from here; a number repeated on that side would silently
  // drop the new one.
  static constexpr std::size_t operating_point_kind_count =
      static_cast<std::size_t>(OperatingPointKind::NonFiniteGradient) + 1;

  // Read-only on purpose: the kind is an output of the solve, and a settable one
  // would be a way to disagree with it. replay_operating_point is the one way to
  // put a recorded kind back.
  OperatingPointKind operating_point_kind() const {
    return operating_point_kind_;
  }
  // Put back an operating point a previous pass found: the collar it placed, and
  // the ARM it placed it on.
  //
  // ⚠️ BOTH, IN ONE CALL, BECAUSE evaluate_root_collar_psi RESTORES EVERY NUMBER
  // AND THEN TAGS THE POINT `prescribed`. That is correct from its side -- a
  // caller handed it a target -- and it is wrong for a replay, where the recorded
  // arm is what a replay needs: operating_collar switches on the kind, so a replayed
  // interior point tagged `prescribed` has no condition to place it at. Splitting
  // this into two calls is the mistake, and it costs a throw in a caller that
  // never asked for a collar.
  //
  // This is the replay half of a recorded decision: the collar MOVES and is
  // re-derived from, the arm SELECTS and is put back. Differentiating the choice
  // would manufacture a jump the model does not have.
  //
  // A BoundarySoil point also needs the polished root of its zero-flux bound,
  // which wet_bound_at's step-back correction reads (wet_bound_root_). The solve
  // finds it inside the optimiser, which a replay does not run, so the replay
  // re-takes the same polish from the same bracket: the rows it hands over are
  // then the fresh solve's, bit for bit.
  double replay_operating_point(double collar, OperatingPointKind kind) {
    return with_curve(cost_curve_, [&](auto tag) {
      return replay_operating_point_for<tag.value>(collar, kind);
    });
  }
  template <CostCurve K>
  double replay_operating_point_for(double collar, OperatingPointKind kind);
  static const char* operating_point_kind_name(OperatingPointKind kind);

private:
  // Where this model replaces a value by a bound. A row severed that way and a
  // row that is honestly zero are the same number, and a consumer DIFFERENCING
  // through a clamp gets a finite value with no response in it and nothing to
  // notice by -- so the derivative kill is invisible unless it is counted.
  clamp_counter clamps;

  // Written by every path out of the collar solve, and reset to Unsolved at the
  // top of prepare_collar_solve. A `Leaf` is often reused for one solve after
  // another, so a branch that declined to write this would leave the PREVIOUS
  // solve's kind: a plausible answer about a different point. Resetting to
  // Unsolved means a path that forgets reports "unclassified" instead.
  OperatingPointKind operating_point_kind_ = OperatingPointKind::Unsolved;
  // The stem's dry end for the solve in progress, dry_end<K>() of the seated
  // curve: psi_crit for a curve bounded there, the end of the stem curve for
  // every other. Seated at the top of prepare_collar_solve, so the bracket
  // helpers that are not templated on the curve read the same number.
  double dry_end_ = std::numeric_limits<double>::quiet_NaN();
  // Whether dry_end_ is psi_crit standing in as a FIRST bracket for a curve it
  // does not bound (prepare_collar_solve's first_bracket), so the dry end of the
  // bracket is not yet the dry end of the domain.
  bool dry_end_provisional_ = false;
  // The wettest rooted layer and the uptake there, from the top of
  // prepare_collar_solve, for extend_collar_bracket.
  double bracket_wettest_ = std::numeric_limits<double>::quiet_NaN();
  double bracket_E_up_at_wettest_ = std::numeric_limits<double>::quiet_NaN();
  // The zero-flux collar a wet pin was solved against, to collar precision, or NaN
  // where no solve placed one. A pinned solve returns a point a step INSIDE it,
  // and wet_bound_at reads the flux at that point; knowing the root lets it correct
  // the residual back to the root, so the row is the root's derivative rather
  // than one proportional to the step.
  double wet_bound_root_ = std::numeric_limits<double>::quiet_NaN();
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


// `namespace detail` used to live here, holding templated REPLICAS of the profit
// algebra whose comment claimed they "mirror Leaf::assim_colimited and
// Leaf::hydraulic_cost_TF exactly". One of them did not -- see the `_kernel`
// declarations in the class. They are deleted: the real functions are now
// templated on their scalar type, so AD differentiates the model itself and the
// mirror cannot drift because there is no mirror.
inline Leaf::Leaf()
    :
    vcmax_25(96), // umol m^-2 s^-1 
    stem_c(2.680147), //unitless
    stem_P50(3.4), // MPa, positive magnitude -- THE trait; stem_b/psi_crit derive
    TF24_beta2(1.5), //exponent for effect of hydraulic risk (unitless)
    jmax_25(157.44), // maximum electron transport rate umol m^-2 s^-1
    phi_J(0.30 / 0.85), // quantum yield of electron transport per ABSORBED photon (mol mol^-1); 0.30 per incident photon at alpha_PAR_ = 0.85
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosythnthesis equatiom
    GSS_tol_abs(1e-3),
    vulnerability_curve_ncontrol(ncontrol_default),
    ci_abs_tol(1e-3),
    ci_niter(1000),
    TF24_cost_scale(7.5) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      // The root traits (root_c/root_P50) and the two beta_R_* resistance
      // constants keep their defaults in MultiLayerRoots, which owns them --
      // including the root's own derived root_b/root_psi_crit. Deliberately not
      // restated here: a second copy of the root Weibull pair is the exact shape
      // of hazard 1 in the developer guide.
      derive_vulnerability_curves_from_P50();
      // ⚠️ THE MEMBER, NEVER A LITERAL. Both of these were literal `100` while
      // ncontrol_default was 100, so a default-constructed leaf built its curves
      // at one resolution and let `set_traits` rebuild them at another the moment
      // the default moved. Raising it to 400 put the two 4x apart and the
      // rescale-against-rebuild identity, asserted to 1e-08, came apart at
      // 1.4e-07 on 52 checks -- at 200, 400 and 800 knots alike, which is the
      // signature of a resolution mismatch rather than of interpolation error.
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
    curv_fact_colim(curv_fact_colim), //curvature factor for the colimited photosythnthesis equation
    GSS_tol_abs(GSS_tol_abs),
    vulnerability_curve_ncontrol(vulnerability_curve_ncontrol),
    ci_abs_tol(ci_abs_tol),
    ci_niter(ci_niter),
    TF24_cost_scale(TF24_cost_scale) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      // The single convention, asserted at the one place it enters (#25). Before
      // there was no global statement about psi's sign to assert -- psi_soil_ was
      // >= 0, psi_soil_inverted_ was <= 0, psi_crit was >= 0 and the collar was
      // <= 0 -- which is precisely why the convention had to live in comments.
      // Now it is checkable, so it is checked. set_traits shares this check.
      check_psi_magnitudes(stem_P50, stem_c, root_P50, root_c);

      // The root traits belong to the supply path, so hand them over before its
      // vulnerability curve is built. Its root_b and root_psi_crit are DERIVED
      // from the pair below, in the one place that derivation lives.
      roots_.root_c = root_c;      //unitless
      roots_.root_P50 = root_P50;  // MPa, positive magnitude
      derive_vulnerability_curves_from_P50();

      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_root_vulnerability(vulnerability_curve_ncontrol);
      setup_clean_leaf();
}

inline void Leaf::check_psi_magnitudes(double stem_P50, double stem_c,
                                      double root_P50, double root_c) {
  // Every psi here is a POSITIVE magnitude in MPa. There is one representation, so
  // this can be asserted globally rather than described in comments.
  if (!(stem_P50 > 0.0)) {
    util::stop("stem_P50 must be a positive magnitude in MPa; got " +
               util::to_string(stem_P50));
  }
  if (!(root_P50 > 0.0)) {
    util::stop("root_P50 must be a positive magnitude in MPa; got " +
               util::to_string(root_P50));
  }
  // The shape parameters set how steeply conductivity is lost. Below about 1 there
  // is no safe plateau at all, which usually flags a measurement artefact rather
  // than a species trait -- but it is representable, so only the sign is refused.
  if (!(stem_c > 0.0)) {
    util::stop("stem_c must be positive; got " + util::to_string(stem_c));
  }
  if (!(root_c > 0.0)) {
    util::stop("root_c must be positive; got " + util::to_string(root_c));
  }
}


// See the header for why this exists rather than fourteen settable fields.
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

  // Both scale parameters and both critical potentials fall out of the traits, so
  // there is nothing left for a consistency check to police.
  const double stem_b_ = weibull_b_from_P50(stem_P50_, stem_c_);
  const double root_b_ = weibull_b_from_P50(root_P50_, root_c_);
  const double psi_crit_ =
      weibull_psi_at_fraction(stem_b_, stem_c_, k_crit_fraction);
  const double root_psi_crit_ =
      weibull_psi_at_fraction(root_b_, root_c_, k_crit_fraction);

  // Which tables have to be rebuilt, decided BEFORE the assignment. A table
  // depends on its curve's c alone (VulnerabilityTable), so a P50 moves nothing.
  // Exact equality is the right test and not a sloppy one: the table is a pure
  // function of c, so an unchanged c gives a bit-identical table and rebuilding
  // it is pure cost. A trait perturbed by a relative 1e-08 -- what a gradient
  // loop does -- compares unequal and rebuilds, which is the case that must not
  // be missed.
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

  // Last, and it is the whole safety argument: back to the just-constructed
  // state. Clears the solved operating point to the NA sentinels (hazard 8: an
  // output left unwritten becomes the previous solve's value), and resets both
  // the transpiration memo and the photosynthesis temperature cache -- the one
  // that would otherwise hand back the old vcmax_.
  //
  // "The just-constructed state" is exact for the derived state and the outputs,
  // which is all of it bar `CF77_lambda_`, `TF24_floor_lambda_o` and
  // `LeastCost_beta` -- caller inputs, left standing on purpose (#96, see their
  // declarations).
  setup_clean_leaf();
}



// ⚠️ THE ROWS ARE IN (P50, c), NOT (b, c). stem_b is derived, so a row in it
// describes a parameter this model does not have -- and set_traits re-derives b
// from P50 on every call, so a direction in b reaches nothing at all. The chain is
// exact: with L = ln 2 and b = P50 * L^(-1/c),
//
//     db/dP50 = b / P50            db/dc = b * ln(L) / c^2
//
// so dG/dP50 is dG/db times the first, and dG/dc is the partial at fixed b plus
// dG/db times the second.
//
// The VALUE is the table's, not the closed form's: the solve ran on the table, so
// a value from the curve would place the operating point somewhere else. Every
// bracket below is exactly zero at the recording point, so the number is
// untouched and only the tape sees the rows.
template <class S>
inline S Leaf::stem_integral_at(const S& psi, const leaf_pars<S>& pars) const {
  // The value is the stem table's; the rows are the closed form's. One site
  // serves both curves -- see closed_form_rows.hpp for why that is not a convenience.
  const double at = odelia::util::to_passive(psi);
  return integral_with_closed_form_rows<S>(
      stem_curve_integral(at, "Leaf::stem_integral_at"),
      stem_curve_integral_deriv(at), psi, pars[par_stem_P50],
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
  // The single-potential path's series resistance. Empty on the layered route,
  // where the resistances are per layer and this slot has no scalar to hold.
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

// set various states and physiology parameters obtained from TF24 to NA to clean leaf object
inline void Leaf::setup_clean_leaf() {
  ci_ = util::na_value; // Pa
  stom_cond_CO2_= util::na_value; //mol Co2 m^-2 s^-1 
  assim_colimited_= util::na_value; // umol C m^-2 s^-1 
  transpiration_= util::na_value; // kg m^-2 s^-1
  residual_flux_ = util::na_value; // kg m^-2 s^-1
  g_min_leaf_ = util::na_value;    // mol m^-2 s^-1, seated with the deficit
  residual_root_ceiling_ = std::numeric_limits<double>::infinity();
  profit_= util::na_value; // umol C m^-2 s^-1
  // CF77_lambda_, TF24_floor_lambda_o and LeastCost_beta are deliberately NOT
  // here: they are the caller's inputs, not derived state, and carry their NA
  // defaults at the declarations instead (#96). Adding any back makes a
  // prescribed price survive set_drivers() and vanish on set_traits().
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
  vcmax_= util::na_value; //kg m^-3
  jmax_= util::na_value; //kg m^-3
  opt_root_psi_ = util::na_value; // MPa, positive magnitude
  // The NA sentinel's counterpart for the operating-point classification: no
  // solve has run, so there is no kind of point to report.
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
  opt_psi_stem_= util::na_value; //-MPa 
  opt_ci_= util::na_value; //Pa
  E_up_ = util::na_value;
  medlyn_model_gs_ = util::na_value; // mol CO2 m^-2 s^-1 (Medlyn model, develop #450)
  theta_w_ = util::na_value;
  theta_fc_ = util::na_value;
  theta_ = util::na_value;
  roots_.clear(); // soil state, geometry and the root resistance network
  // BOTH supply paths, not just the active one. Hazard 8 is that an output a
  // code path declines to write becomes the previous solve's value, and a Leaf
  // that has been switched between paths (set_supply_single / _multilayer) is
  // exactly the case where the inactive one's stale soil state could come back.
  // clear() leaves grav_head_ alone -- it is the single path's one remaining piece
  // of configuration, and wiping it here would make set_supply_single
  // order-dependent. resistance_ IS cleared, because it is a driver now.
  single_.clear();
  soil_consumption_.clear(); // soil consumption mol  m^-2 s^-1;

  transpiration_cached_ = false; // invalidate transpiration() memo
  photo_temp_cached_ = false;    // members above set to NA; force recompute
}

// Set the per-individual, per-timestep physiology that stays constant during
// the subsequent root-collar/stem optimisation.
//
// Two groups of quantities are computed here:
//   1. Temperature-dependent photosynthetic parameters (vcmax_, jmax_,
//      electron_transport_, gamma_, ko_, kc_, km_, R_d_) via Arrhenius/peaked
//      Arrhenius functions, plus assim_max_ (assimilation at ci = ca).
//   2. The soil state and the root hydraulic-resistance network, both of which
//      now belong to roots_ (MultiLayerRoots::set_soil_state and
//      ::set_root_network). The two calls are left at their original positions
//      in this function, straddling the temperature block, rather than merged:
//      the blocks are numerically independent, but keeping the order means the
//      only thing this refactor has to argue about is where the code lives.
//
// NOTE: the temperature-dependent block (group 1) is a function of one
// temperature, recorded in photo_temp_. Off the energy-balance path that is
// leaf_temp_ and constant across the run in the current driver setup, which is
// what the cache below exploits; on it the block is re-seated per candidate psi
// at the leaf's own temperature and the cache is bypassed.
//
//sets various parameters which are constant for a given node at a given time
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
  // The deficit of a real atmosphere: positive, and short of saturation at the
  // air temperature, so the actual vapour pressure e_s(Tair) - D is positive.
  // Refused rather than floored: the longwave term used to floor e_a at
  // 1 Pa, which bent the leaf temperature by 1.5 C across RH = 0 and was reached
  // only by deficits past saturation. Off the energy balance leaf_temp IS the
  // air temperature the deficit is stated at.
  if (!(atm_vpd > 0.0) || !(atm_vpd < saturation_vapour_pressure(leaf_temp))) {
    util::stop("set_physiology: atm_vpd must lie in (0, e_s(T_air)), a relative "
               "humidity between 0 and 1; got atm_vpd = " +
               util::format_double(atm_vpd) + " kPa against e_s = " +
               util::format_double(saturation_vapour_pressure(leaf_temp)) +
               " kPa at " + util::format_double(leaf_temp) + " C");
  }
  // The thermal-damage states the caller seated. The model reads A_ and I_ as
  // themselves, so a hair outside their range is smooth; what is refused is a
  // negative fraction of leaf area and the acclimation shift's pole. Infeasible
  // rather than a plain error: a consumer integrating these states reaches them
  // from a trial stage, and can reject the step. A NaN passes, to poison the
  // outputs where a test can see it.
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
    // The input boundary for the one representation (#25). A caller that still
    // has the pre-#25 signed vector fails here rather than silently running a
    // model with the soil and the collar on opposite sides of zero.
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
  // The root network is NOT validated here. It is checked in set_root_network
  // below, which is where the length agreement with the soil profile can actually
  // be tested -- the network is sized to the deepest ROOTED layer, so it is
  // shorter than psi_soil whenever the plant is shallower than the soil, and a
  // check up here would have to duplicate that rule.
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
       // One potential; the depth profile and the root-mass profile are not
       // this path's business, and the caller's vectors are simply not read
       // beyond element 0.
       single_.set_soil_state(psi_soil[0]);
       break;
   }

   leaf_specific_conductance_max_ = leaf_specific_conductance_max;
   // conductance changed -> invalidate the transpiration() memo
   transpiration_cached_ = false;
   ca_ = ca;
   // Penman-Monteith leaf energy-balance inputs (#523). Reinterpret the incoming
   // leaf_temp driver as air temperature. The isothermal net radiation is the
   // ABSORBED shortwave -- each intercepted band times its own absorptance, the
   // PAR band converted from photons to W m^-2 at the daylight ratio -- less the
   // sky's longwave (#97). Only read on the use_energy_balance_ path.
   if (use_energy_balance_ && !std::isfinite(NIR_int_)) {
     util::stop("set_physiology: non-finite NIR_int_ on the energy-balance path. "
                "Set the intercepted near-infrared (W m^-2) before set_physiology: "
                "about PPFD / umol_par_per_joule in daylight, 0 under an LED; "
                "NIR_int_=" + util::to_string(NIR_int_));
   }
   Tair_ = leaf_temp_;
   {
     // Isothermal net radiation and the radiative conductance (#97). Both depend
     // only on the DRIVERS, so both are seated once here rather than per candidate
     // potential -- Tleaf enters the balance linearly through g_rad_ and not
     // through either of these.
     const double Tair_K = Tair_ + C_to_K;
     // Actual vapour pressure, IN PASCALS as the emissivity coefficient requires
     // (constants.hpp). atm_vpd_ and esat are both kPa, hence the 1000. Positive:
     // a deficit at or past saturation is refused above.
     const double ea_Pa = 1000.0 * (saturation_vapour_pressure(Tair_) - atm_vpd_);
     // ⚠️ HELD BELOW 1, A PHYSICAL BOUND AND NOT PADDING. Brutsaert's fit
     // extrapolates past unity once ea > 22.25 Pa/K * Tair_K (about 7 kPa): above
     // Tair 38.9 C in saturated air, 43.9 C at a 2 kPa deficit, 48.0 C at 4 kPa.
     // Unbounded, the isothermal longwave turns POSITIVE there (+100 W m^-2 by
     // Tair 60 at 2 kPa): a sky brighter than a blackbody heating the leaf, the
     // right order of magnitude with the wrong sign. The bound is a SOFT minimum:
     // a min() at 1 put a corner in every output's slope in Tair and
     // humidity on a heatwave sweep. Past the fit's unity the term tends to 0,
     // the limit for a sky radiating as a blackbody at Tair, and Rn_ responds to
     // Tair and humidity only weakly there. It is computed once per call, so the
     // soft form costs nothing per candidate.
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
   // Aerodynamic resistance from leaf boundary-layer theory: ra = C_ra*sqrt(d/U)
   // (doc 4.1), with the per-strategy leaf dimension d_ and the above-canopy wind
   // wind_speed_. On the PM path a non-finite wind_speed_/d_ is a broken driver /
   // unset trait, not a modelling choice, so fail fast rather than silently using
   // the fixed fallback (review: itowers1). Zero wind or zero d (physically
   // ra -> infinity) is a legitimate case and falls back to the fixed ra.
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
   // The leaf temperatures the solve can reach. With the balance on, Tleaf is
   // affine in E, so the range runs from the shut leaf (E = 0, the hottest) to
   // the most the stem can carry, k_max G(inf) (the coolest; every supply path
   // carries less, because G >= 0 upstream). Off it the driver is the leaf's.
   if (use_energy_balance_) {
     const double E_most = stem_curve_
                               ? leaf_specific_conductance_max_ * stem_curve_limit()
                               : 0.0;
     check_leaf_temperature_range(leaf_temp_from_E(E_most), leaf_temp_from_E(0.0));
   } else {
     check_leaf_temperature_range(leaf_temp_, leaf_temp_);
   }

   // The diffusion deficit at the baseline temperature. On the PM path the solve
   // re-seats this per candidate psi; here it makes assim_max_ and any bare
   // stom_cond_CO2() call well defined before one has run.
   set_leaf_vpd(leaf_temp_);

   // Temperature/O2-dependent block. Off the PM path this is recomputed only
   // when (leaf_temp_, atm_o2_kpa_) changes from the previous call (see
   // photo_temp_cache_ in the header); same inputs -> bit-identical outputs, so
   // reusing is exact. On the PM path the cache is bypassed: leaf_temp_ here is
   // only Tair, and the operating-point Tleaf (hence these params) is set per
   // candidate psi in set_leaf_states_rates_from_psi_stem, so we always recompute
   // the Tair baseline (reported as assim_max_) and let the solve
   // override it.
   //
   // ⚠️ AND IT IS BYPASSED WITH THE THERMAL-DAMAGE GATE ON. The block then reads
   // A_, I_ and D_, which a caller moves between calls and which the key
   // does not hold (see photo_temp_key).
   const std::array<double, photo_temp_key_size> photo_key = photo_temp_key();
   if (!use_energy_balance_ &&
       !use_thermal_damage_ &&
       photo_temp_cached_ &&
       photo_key == photo_temp_cache_key_) {
     // Cache hit (non-PM): every input unchanged; only electron_transport_
     // depends on the per-call PPFD_, so refresh just that (as before).
     electron_transport_ = electron_transport();
   } else {
     update_temperature_dependent_params(leaf_temp_);
     photo_temp_cache_key_ = photo_key;
     photo_temp_cached_ = true;
   }

  // The supply resistances, handed over as given. #33 moved the carbon ->
  // resistance step out to the caller, so this is a validated copy rather than a
  // model evaluation. BOTH paths are served from the one argument: the
  // single-potential path reads `r_R_V_sum[0]` as its series resistance, which is
  // that field's own meaning with one layer and no horizontal term. That is what
  // makes the calling convention independent of which path is in force.
  //
  // It stays HERE, after set_soil_state, because the multi-layer length check
  // needs the soil profile seated first.
  switch (supply_kind_) {
    case SupplyKind::MultiLayer:
      roots_.set_root_network(root_network);
      break;
    default:
      single_.set_supply_resistances(root_network);
      break;
  }

  // Set up vector of root water uptake from layer. Stays on Leaf: plant writes
  // the crown-integrated value back into leaf.soil_consumption_ by name.
  //
  // .assign, not .resize: the uptake loop writes only up to the deepest *rooted*
  // layer, and resize's fill reaches only newly added elements, so a solve with
  // fewer rooted layers than the last one on this Leaf would leave the tail
  // holding the previous plant's values -- which plant then bills to the patch
  // water balance.
  soil_consumption_.assign(supply_n_layers(), 0.0);

  // Soil-moisture state for the Medlyn beta_ stress factor (develop #450). The
  // root-water compute path does not use these; they make the standalone,
  // R-callable Medlyn methods well-defined with the default soil-moisture values.
  theta_w_ = theta_w;
  theta_fc_ = theta_fc;
  theta_ = theta;

  // Find maximum assimilation assuming ci = ca
  assim_max_ = assim_colimited(ca_);
}

// ===========================================================================
// ONE REPRESENTATION FOR WATER POTENTIAL  [#25, superseding #7's two-convention map]
// ---------------------------------------------------------------------------
// EVERY psi in this package is a POSITIVE MAGNITUDE in MPa. There is no second
// representation, no flip, and no bridge point. `psi_soil_`, `psi_crit`,
// `root_psi_crit`, `stem_b`, `root_b`, `opt_psi_stem_`, `opt_root_psi_`, the
// find_root_psi / E_column root variable, find_psi_stem_from_psi_root's
// psi_root, transpiration's psi_upstream, and the vulnerability tables' arguments are the
// same kind of number. `set_physiology` asserts psi_soil >= 0 and the
// constructor asserts psi_crit > 0 and stem_b > 0, so the convention is
// checkable rather than merely documented -- which is the point: before #25
// there was no global statement to assert, because the answer depended on which
// variable you asked about.
//
// Where an equation needs one potential to oppose another, THE MINUS SIGN IS IN
// THE EQUATION, not in the storage. The soil -> collar flux is
//
//     E_i = (T_collar - T_soil_i - gravity_head * z_i) / r_R
//
// which reads directly: to draw water you must pull harder than the soil holds
// it, plus enough to lift it. E_i < 0 still means the layer is gaining water.
//
// What this buys, beyond tidiness: dE_up/dT_collar is now +1/r, a conductance,
// positive by construction. Under the signed convention the same derivative came
// back negative where a conductance was wanted, and that produced a negative
// lambda in marginal_cost_water_multilayer -- a bug that had to be patched with a
// negation. That whole class of error cannot occur here.
//
// ⚠️ DO NOT call this convention "tension" in the code. It is exact only because
// there is no osmotic term anywhere and the gravitational component is carried
// separately as `gravity_head`, so psi here is a pure pressure potential. Add
// solutes and "tension" becomes wrong while "magnitude of psi" stays right.
//
// TWO THINGS #7 GOT RIGHT AND #25 KEEPS:
//
//   * transpiration() and its inverse transpiration_to_psi_stem() used to read
//     eval(psi_upstream) and eval(-psi_upstream) respectively. #7 recorded that
//     as "NOT a bug -- they are called with psi_upstream of OPPOSITE sign". True,
//     and it stopped being necessary: they now take the same convention, and the
//     inverse's negation is gone.
//   * opt_root_psi_ (exported as the opt_root_psi aux) and opt_psi_stem_ agree in
//     ALL branches of find_root_collar_psi. #7 established that for opt_psi_stem_
//     and for the collar's sign-consistency across exits; #25 makes them the same
//     kind of number as each other, and as plant's state variable of the same
//     name, whose two compensating negations are deleted.
//
// The member is called `opt_root_psi_`, not `root_collar_psi_`. Renamed
// deliberately: keeping the old name with a flipped sign is the one genuinely
// dangerous outcome here, because an old analysis would silently read the wrong
// sign. A rename gives a binding error instead.
// ===========================================================================
//
// This function is used to find root collar suction which equilibrates the
// soil-root-stem water continuum. `x` is the candidate collar suction; it and
// psi_leaf are the same kind of number now, so the mid-solve scratch write into
// the collar member that #7 flagged (a magnitude parked in a member documented as
// signed) is simply gone -- it never needed to be a flip.
inline double Leaf::E_column(double x, const std::vector<double>& psi_soil, double psi_leaf) {

  E_from_Soil_to_Root_Collar(x, psi_soil);
  double E_root_to_leaf = transpiration(psi_leaf, x);
  return E_up_ - E_root_to_leaf;
}

// This function is used to find the root collar suction where water from soil is zero
inline double Leaf::E_column_zero(double x, const std::vector<double>& psi_soil) {

  E_from_Soil_to_Root_Collar(x, psi_soil);

  return E_up_;
}

// find root psi based on required condition, i.e. equilibrated continuum, zero water from soil
//
// #486: both targets (E_column / E_column_zero, the soil->collar continuity
// residual over the collar suction x in [wettest_soil_layer, dry_end_]) are
// smooth and strictly monotone in their *normal operating regime* -- a clean
// single sign-change with derivatives continuous across every x == psi_soil[i]
// layer crossing (the per-layer branch switches in E_from_Soil_to_Root_Collar
// are bit-level kinks, relative slope jump ~1e-6, not real corners). So a
// superlinear bracketing solver is safe and faster here: TOMS748 reaches the
// same root in ~6-8 E_from_Soil evals vs bisection's ~15-16 at the same 1e-4
// tol (see the test-leaf.r "find_root_psi soil->collar continuity solve"
// contract block). This directly attacks the dominant E_from_Soil per-layer
// arithmetic hot-spot, evaluated ~270x per collar solve through this finder.
//
// SCOPE/CAVEAT: the brackets here are guaranteed valid (opposite-sign, finite
// endpoints) by find_root_collar_psi's preceding early-exits -- the crit=1
// dry_end_ endpoint is exactly the E_column(dry_end_) < 0 shutdown test, and
// crit=0 is only reached for soil wetter than dry_end_. The genuinely
// non-smooth failure mode in the earlier blanket-swap rejection (the root
// vulnerability spline extrapolating negative beyond its domain, before #183)
// lived in E_from_Soil_to_Root_Collar itself and bit both solvers identically;
// it did not arise on the brackets this finder is actually handed. Like
// psi_stem_to_ci (Phase 6) this is a same-tolerance method swap, NOT a tolerance
// loosening: same root, fewer evals.
inline double Leaf::find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit,
                                  double f_wettest, double f_crit, double end) {
  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  if (find_root_crit == 1) {
    auto target = [&](double x) -> double {
      return E_column(x, psi_soil, end);
    };
    try {
      // Ends swapped, not just renamed: a bracketing solver needs opposite signs
      // at the two endpoints, but the LOWER bound must come first, and in
      // magnitudes the wettest layer is the smallest suction (#25).
      //
      // ⚠️ THE WET END OF THE FINAL BRACKET, NOT ITS MIDPOINT. The target rises
      // with the collar, so its wet end is a collar at which the soil delivers
      // no more than the stem carries with its potential at the dry end. The
      // midpoint can sit up to bracket_tol past the root, where it delivers
      // more; with the dry end at the end of the stem curve, that is more than
      // G(inf), which no stem potential carries, and the inverse refuses it.
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

// ⚠️ WHY A PINNED OPTIMUM NEEDS ITS BOUND SOLVED TIGHTLY. find_root_psi stops at
// 1e-4 MPa, which is plenty for a BRACKET -- the interior root-find refines to
// collar_root_tol inside it whatever the ends are. But at a pinned optimum the
// bound is the answer, and where a 1e-4 root-find stops moves with the parameters
// through its own iterates, not the way the root does. The derivative rows
// (operating_collar -> wet_bound_at) differentiate the root itself, so the two disagreed:
// at psi_soil = 4 the zero-flux bound sat 1.3e-4 MPa from the true root, its
// d/droot_c was 1.2e-4 against the root's 1.1e-6, and the profit rows for the
// root traits came out 100x off a finite difference of the solve.
//
// Solving every bound this tightly costs +2% on every solve; polishing only the
// bound that answers costs nothing on an interior one. The root is within ~1e-4
// of x0, so the first bracket tried almost always holds it.
inline double Leaf::polish_bracket_bound(double x0,
                                         const std::vector<double>& psi_soil,
                                         int find_root_crit) {
  // The wet bound's target is uptake less the residual flux, which at g_min = 0
  // is the uptake exactly.
  auto target = [&](double x) -> double {
    if (find_root_crit == 1) {
      return E_column(x, psi_soil, dry_end_);
    }
    const double up = E_column_zero(x, psi_soil);
    return up - residual_flux(x, up);
  };
  // Tighter than collar_root_tol, which at a collar of ~5 MPa is ~6e-12: near
  // the dry end the stem curve is steep enough to turn that into ~2e-11 of psi_stem,
  // which a finite difference at a 1e-6 step reads as a 2e-6 error in a row.
  constexpr double polish_tol = 1e-15;
  const double f0 = target(x0);
  if (f0 == 0.0 || !std::isfinite(f0)) {
    return x0;
  }
  // Both targets rise with the collar, so f0 > 0 puts the root below x0; that side
  // is tried first and the other only if it holds no sign change. The first width
  // is the bracket's own tolerance at x0, which holds the
  // root it was placed around, so one probe usually suffices.
  //
  // ⚠️ THE CONTINUITY ROOT'S PROBES STOP AT THE DRY END, as find_root_psi's
  // bracket does: past it the collar is drier than the stem held there, and
  // transpiration() refuses a reversed gradient.
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

// When root pressure is known, find E from soil, then use E from soil to find psi stem
inline double Leaf::find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil){
  E_from_Soil_to_Root_Collar(psi_root, psi_soil);

  double psi_stem = transpiration_to_psi_stem(E_up_, psi_root);
  return psi_stem;
}

// ---------------------------------------------------------------------------
// MASTER SOLVER: optimal root-collar (and stem) water potential
// ---------------------------------------------------------------------------
// find_root_collar_psi is the entry point, called once per solve. It solves the whole
// soil -> root -> stem -> leaf hydraulic continuum and stores the optimal
// operating point in opt_psi_stem_, opt_root_psi_ and profit_.
//
// The solve has two nested levels:
//
//   1. CONTINUITY (find_root_psi / E_column): for any candidate root-collar
//      potential, water supplied from the soil (E_from_Soil_to_Root_Collar)
//      must equal water transpired through the stem (transpiration()). This is
//      a 1-D root-find on the collar potential.
//
//   2. OPTIMISATION: among feasible collar potentials, choose the one that
//      maximises carbon profit = assimilation - hydraulic cost. This is a
//      root-find (TOMS748) on dprofit/dcollar = 0, with golden-section search
//      only as a fallback when an endpoint has no usable gradient. The collar
//      potential is bracketed between the wet bound (the collar at which the
//      soil delivers the residual flux; at g_min = 0, where uptake is zero) and
//      `collar_at_dry_end` (the collar at which the stem must sit at its dry
//      end to carry what the soil delivers).
//
// Several early exits skip the optimisation when no meaningful choice is
// possible. In each case the plant is effectively
// shut down (moving no water, paying the hydraulic cost where flow stopped):
//   * wettest soil layer is already drier than psi_crit -> no transpiration;
//   * even at psi_crit the soil cannot supply the demanded flux (E_column<0);
//   * the continuity root would require the collar drier than psi_crit;
//   * maximum possible assimilation (at ci = ca) is negative.
//
// psi_soil is a positive magnitude here as everywhere else, so nothing is
// flipped.
//
// The collar at which the soil's uptake equals the residual flux the leaf
// realises there: the wet bound of the collar solve, and the shut leaf's own
// collar. At g_min = 0 it is the zero-uptake collar, found by the call
// find_root_psi(.., 0, ..) always made, with the same bracket and the same
// endpoint values, so that path is bit-identical.
//
// With g_min > 0 the target is uptake minus the residual flux, which rises
// with the collar: uptake does, and the residual flux falls with it (the
// ceiling tightens, and on the energy-balance path a larger uptake cools the
// leaf). It is negative at the wettest layer, where uptake is at most the
// gravity term, and positive once the soil delivers the cuticle's loss. That
// can be past the dry end, so the bracket is widened out from there until the
// target changes sign, which it must: the residual flux is capped below the
// roots' ceiling.
//
// `E_up_at_crit` is the uptake at the dry end where the caller has it (NaN
// otherwise), handed over rather than re-taken for the reason the two bracket
// root-finds share their endpoint values.
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
  // The first bracket is [wettest, psi_crit], whenever it is a bracket at all,
  // whatever the dry end of this solve: the wet bound is then the same root
  // from the same bracket on every route to it, and a replay re-takes it bit for
  // bit. The uptake there is handed over where the caller has it.
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

// The shut leaf's wet bound, bracketed from below by the supply's floor
// (supply_wet_floor): at the lowest gravity balance among conducting layers
// every layer draws nothing or less, so uptake less the residual flux is at
// most zero there, and above it the target rises. Two things follow:
//
//   - at g_min = 0 on one layer the floor IS the bound (the draw is exactly
//     zero at the balance), and on several the bound lies between the lowest
//     and highest balance, so ONE root-find to collar precision on that bracket
//     replaces a walk out from the wettest soil potential, a coarse solve and
//     a polish -- most of what made a shut solve cost 3.6x an old held one.
//     With a residual conductance the bound lies higher, and the walk starts
//     from the highest balance;
//   - where no layer conducts the target is zero everywhere and the floor is
//     returned as it stands: the wettest layer's balance, the limit the bound
//     reaches from the wet side. Returning the soil potential instead, as the
//     walk did, left a step of the layer's gravity head in psi_stem (#184).
//
// The open leaf keeps find_wet_bound and its bracket, which is what keeps its
// rows bit-identical at g_min = 0.
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
  // To collar precision directly, polish_bracket_bound's 1e-15, for the reason
  // it gives: the bound IS the shut leaf's answer.
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

// The stem potential from which the stem carries exactly the residual flux,
// with the upstream potential pinned: the wet end of a stem route. The flux
// through the stem rises with the stem potential towards its ceiling, and the
// realised residual flux is strictly below that ceiling, so the root exists;
// it is beyond psi_crit when the soil is too dry to supply the cuticle inside
// it. Solved tightly, because a leaf that declines to open sits exactly here.
inline double Leaf::stem_wet_bound(double psi_upstream) {
  // Every stem route starts here, and a stem route has no roots upstream of
  // psi_soil, so the roots' ceiling a previous collar solve seated does not
  // apply to it.
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

// The shut operating point the collar solve's early exits place (#182).
//
// With its stomata shut the leaf still loses the residual flux, so water
// still moves: the collar sits where the soil's uptake meets that flux, the
// stem where it carries it, and every output -- the flux, the conductance the
// cuticle leaves, the intercellular CO2 and the carbon -- comes from the
// ordinary evaluator at that point. Nothing is written by hand, so the shut
// state and the open state next to it are the same function, which is what
// makes profit continuous across the exit (#175, #183). At g_min = 0 the
// collar is the zero-uptake collar and the stem sits at it, the E -> 0 limit
// of the open leaf; the old convention of holding the stem at a critical
// potential went with the no-flow family.
//
// Tagged here rather than at the call sites, so a fifth reason to shut down
// cannot arrive without a classification. All of them are one statement: no
// collar potential both carries more than the residual flux and stays inside
// the critical potentials, so there is nothing to optimise over.
template <Leaf::CostCurve K>
inline void Leaf::set_shut_state() {
  const std::vector<double>& soil = supply_psi_soil();
  // The bound IS the answer here, so find_shut_wet_bound solves it to collar
  // precision itself, for the reason polish_bracket_bound gives.
  const double collar = find_shut_wet_bound();
  wet_bound_root_ = collar;
  opt_root_psi_ = collar;
  E_from_Soil_to_Root_Collar(collar, soil);
  // The stem carries the uptake, as it does at every collar. At g_min = 0 that
  // uptake is zero to the polish tolerance, and E = 0 puts the stem AT the
  // collar exactly rather than a round-trip of the inverse table away from it.
  opt_psi_stem_ = residual_flux(collar, E_up_) == 0.0
                      ? collar
                      : find_psi_stem_from_psi_root(collar, soil);
  profit_ = profit_psi_stem_for<K>(opt_psi_stem_, collar);
  // The seated curve's marginal cost of water at this point, NA where the
  // supply's conductance is unavailable (a one-layer wet bound is its
  // gravity-balance kink).
  lambda_emergent_ = lambda_collar_for<K>();
  operating_point_kind_ = OperatingPointKind::HydraulicShutdown;
  if (!std::isfinite(profit_)) {
    util::stop_infeasible("collar_solve",
                          "profit is not finite at the shut state; collar=" +
                              util::to_string(collar) +
                              "; psi_stem=" + util::to_string(opt_psi_stem_));
  }
}

// Shared setup + feasibility handling for the root-collar solve. Extracted
// verbatim from find_root_collar_psi (no reordering of floating-point ops, so
// the TF24 optimisation path stays bit-identical) and reused by
// evaluate_root_collar_psi. Returns false when the final operating point is
// already determined here (shutdown / closed interval) and the
// caller should stop; returns true with [bound_a, bound_b] set to the feasible
// collar-potential interval (positive magnitudes) otherwise.
template <Leaf::CostCurve K>
inline bool Leaf::prepare_collar_solve(double& bound_a, double& bound_b,
                                       bool first_bracket){

  // Clear the classification FIRST, so that a path which declines to write it
  // reports "unclassified" rather than the previous plant's kind of operating
  // point (hazard 8). Every exit below, and both callers, write it again.
  operating_point_kind_ = OperatingPointKind::Unsolved;
  wet_bound_root_ = std::numeric_limits<double>::quiet_NaN();

  // Hand the supply path the start of a solve: it builds its per-solve caches and
  // reports back the wettest rooted layer -- the SMALLEST suction, and the lower
  // bracket endpoint everything below needs.
  const double wettest_soil_layer = supply_begin_solve();
  // The roots' ceiling on the residual flux, seated before anything reads the
  // residual flux. Only the layered supply bounds its uptake; a constant series
  // resistance delivers without limit as the collar dries.
  residual_root_ceiling_ =
      (g_min > 0.0 && supply_kind_ == SupplyKind::MultiLayer)
          ? roots_.uptake_limit<double>(roots_.held_supply())
          : std::numeric_limits<double>::infinity();
  const double E_up_at_wettest =
      E_column_zero(wettest_soil_layer, supply_psi_soil());

  bracket_wettest_ = wettest_soil_layer;
  bracket_E_up_at_wettest_ = E_up_at_wettest;

  // One bracket, with the stem held at `end`: whether any collar short of it
  // can both supply the cuticle and leave the stem something to carry, and if
  // so the interval [wet bound, continuity root]. The bracket helpers read the
  // end through dry_end_, so it is seated first.
  //
  // ⚠️ THE TWO BRACKET ROOT-FINDS SHARE THEIR ENDPOINT VALUES, and they are
  // handed over rather than re-taken. At the end the column term is exactly
  // zero, so E_column there IS the uptake E_column_zero reads; at the wettest
  // layer both read one uptake.
  auto bracket = [&](double end) -> bool {
    dry_end_ = end;
    if (wettest_soil_layer >= end) {
      return false;
    }
    // Nor where the stem can carry less, from the wettest layer, than the curve
    // resolves (stem_capacity_floor): no gradient there is more than rounding.
    const double G_inf = stem_curve_limit();
    if (!(G_inf - stem_curve_integral(wettest_soil_layer) >
          stem_capacity_floor * G_inf)) {
      return false;
    }
    // Uptake at a collar held at the end is the most the soil delivers with the
    // stem short of it; below the residual flux the cuticle cannot be supplied
    // (at g_min = 0, a negative uptake).
    const double E_at_end = E_column(end, supply_psi_soil(), end);
    if (E_at_end < residual_flux(end, E_at_end)) {
      return false;
    }
    const double collar_at_dry_end = find_root_psi(
        wettest_soil_layer, supply_psi_soil(), 1,
        E_up_at_wettest - transpiration(end, wettest_soil_layer), E_at_end, end);
    // The continuity root at or past the end: no collar is left.
    if (collar_at_dry_end >= end) {
      return false;
    }
    // The wet bound: the collar at which the soil delivers the residual flux,
    // and at g_min = 0 the collar at which uptake is zero.
    bound_a = find_wet_bound(wettest_soil_layer, E_up_at_wettest, E_at_end);
    // The dry end of the interval: the collar at which the stem must sit at
    // `end` to carry what the soil delivers. The root's own critical potential
    // bounds nothing (#183); the clamp to it, and the tie with the stem's that
    // made every default-trait shutdown a kink in the stem and root traits, are
    // gone.
    bound_b = collar_at_dry_end;
    // An interval narrower than bracket_tol has ends that cannot be trusted, so
    // they are solved to collar precision before anything is decided on it.
    if (bound_b - bound_a < bracket_tol_at(bound_b)) {
      bound_a = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
      bound_b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    }
    return bound_b - bound_a > collar_root_tol;
  };

  // For a curve not bounded at psi_crit, psi_crit is only a first bracket: it
  // holds the optimum wherever profit already falls there, at the old cost, and
  // maximise_profit_over_collar extends it otherwise. Where it holds no interval
  // the end of the stem curve is tried before the leaf is called shut.
  if constexpr (!cost_domain_ends_at_psi_crit<K>()) {
    if (first_bracket && bracket(psi_crit)) {
      dry_end_provisional_ = true;
      return true;
    }
  }
  dry_end_provisional_ = false;
  // Every exit below places the leaf at its wet bound with the stomata shut,
  // through one function: no collar short of the dry end both supplies the
  // cuticle and moves water.
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
  // ⚠️ THE DRY END IS POLISHED HERE TOO, as prepare_collar_solve polishes a
  // narrow first bracket. find_root_psi places the root only to within
  // bracket_tol_at of it, so where the whole interval is narrower than that (a
  // small E times the series resistance) the unpolished root can fall on the
  // wet bound, and a leaf whose profit still rises is reported shut. The wet
  // end is not re-polished: it was polished already wherever the first bracket
  // was narrow, a wide first bracket only widens when its dry end moves drier,
  // and moving the wet end here would strand the caller's `lo`, stepped in
  // from it before this call.
  const double bound_b = std::max(collar_at_dry_end, bound_a);
  if (bound_b - bound_a < bracket_tol_at(bound_b)) {
    return std::max(polish_bracket_bound(bound_b, supply_psi_soil(), 1),
                    bound_a);
  }
  return bound_b;
}

// The profit-maximising collar potential, by solving the first-order condition
// dprofit/dpsi == 0 instead of searching profit itself. Assumes
// prepare_collar_solve has run, so the soil-side caches are seated and every
// gradient evaluation below can go straight to dprofit_at_collar_psi.
//
// WHY this replaced golden section, in one line each:
//
//  * Golden section resolves the argmax only to GSS_tol_abs (1e-3), leaving the
//    returned collar ~1e-4 MPa from the true stationary point with an offset that
//    wanders as the comparison sequence flips. That offset IS the staircase that
//    made trait derivatives unusable.
//  * The damage was not merely a staircase. For traits in the hydraulic path the
//    wandering offset DOMINATED the real trait response, so the argmax came back
//    smooth, plausible and SIGN-INVERTED (root_b: -2.6e-03 where the truth is
//    +2.6e-04). Arbitrated against a 20001-point derivative-free scan of profit.
//  * It cannot be fixed downstream. No finite-difference step size and no amount
//    of differentiating through the iterations recovers it, because the error is
//    in the solved argmax rather than in how the argmax is differentiated.
//  * Hazard 3 -- "the argmax must vary smoothly with inputs" -- is IMPROVED, not
//    threatened: measured ~800x better second differences in the traits. The
//    hazard's own concern is smoothness in plant state, which cannot be measured
//    here while plant #591 blocks end-to-end validation; same mechanism, so the
//    same direction is expected, but it is not the same measurement.
//
// ⚠️ AND WHY IT IS NOT `util::maximise_over_closed_interval_foc` AT `n = 0`,
// which is the algorithm it runs. What it has and that function does not is the
// CLASSIFICATION: every exit below sets an `OperatingPointKind`, and 42 of the 240
// feasible golden-grid rows are pinned and read it. Folding the two together means
// either handing the shared solver an out-parameter no other caller wants, or
// dropping the tag. Do not read this as a second solver -- it is the one solver
// plus a report of which branch answered.
//
// Safeguarded, not Newton: TOMS748 keeps a bracket throughout, so a non-monotone
// or kinked case degrades instead of diverging. Monotonicity was measured on all
// 240 feasible golden-grid rows, but that is a grid, not a theorem. This also
// makes the collar the third leaf solver on this method -- see uniroot.hpp's
// history note, which is about exactly this judgement call.
template <Leaf::CostCurve K>
inline double Leaf::maximise_profit_over_collar(double bound_a, double bound_b) {
  double width = bound_b - bound_a;
  const auto dprofit = [&](double psi, bool* feasible = nullptr) {
    return dprofit_at_collar_psi<K>(psi, feasible);
  };

  // The WET endpoint. With a residual conductance it is an ordinary point: the
  // stomatal flux is zero there but the cuticle still conducts, so dprofit is a
  // real derivative and the bound is tried as-is first. At g_min = 0 it is
  // infeasible by construction -- bound_a is the zero-uptake collar, where
  // psi_stem == bound_a and dprofit takes its reversed-gradient exit. Its 0.0
  // is a sentinel, and a bracketing solver handed it would report the
  // zero-transpiration point as the optimum: profit -1.897 there against 2.516
  // at the true optimum, measured at the default drivers.
  //
  // So step inside until the gradient is real. Measured over the golden grid the
  // infeasible region is at most 3.46e-07 MPa wide (median 1.22e-08, never more
  // than 6.2e-07 of the bracket), so the first step lands and the loop is a guard.
  //
  // ⚠️ THE BOUND IS TRIED ONLY WHEN THE CUTICLE CONDUCTS. At g_min = 0 the
  // bracket's root sits within its 1e-4 tolerance of the true zero-uptake
  // collar, on either side, so the bound sometimes admits a tiny real flux and
  // a real derivative with it; starting the interior root-find there rather
  // than a step inside moves every interior row by the solver floor. Keeping
  // the old first step at g_min = 0 is what keeps the golden files bit-identical.
  bool ok_lo = false;
  double lo = bound_a;
  double f_lo = 0.0;
  for (double frac = (g_min > 0.0 ? 0.0 : 1e-6); frac < 0.5;
       frac = (frac == 0.0 ? 1e-6 : frac * 10.0)) {
    lo = bound_a + frac * width;
    f_lo = dprofit(lo, &ok_lo);
    if (ok_lo && std::isfinite(f_lo)) {
      break;
    }
  }

  // The DRY endpoint. For a curve bounded at psi_crit it is feasible everywhere
  // measured, so it is tried as-is first and only stepped inside if that fails.
  //
  // For every other curve the dry end is the end of the stem curve, where G is
  // flat to rounding: the stem potential there is wherever G first rounds to
  // G(inf), and every collar derivative divides by a G' of ~1e-16, so the
  // gradient read at the bound is noise (measured: +0.99 at the bound, -35 one
  // step of 1e-6 of the bracket inside it, -14.6 at 1e-2). So the first try is
  // that step inside, and a gradient that is not yet negative steps further in
  // rather than being read as a pin: the curve's marginal cost of water diverges,
  // so a falling profit is what the dry end has.
  //
  // ⚠️ AND PROFIT CAN STILL BE RISING THERE, which is the model and not noise.
  // These curves' marginal cost of water diverges only as x^((c-1)/c) along the
  // stem curve, so where the stem has little left to carry (a dry soil) and A
  // still climbs steeply in E, the optimum lies past the end of the curve: the
  // leaf transpires everything the stem can carry, to within e^-40 of it. Every
  // output is defined there; only the stem potential is not (any potential past
  // the end carries the same flux). That is a pin at the dry end, reported as
  // one: the continuity root with the stem at the end of its curve.
  //
  // On a first bracket (prepare_collar_solve's first_bracket) the dry end is
  // the continuity root with the stem at psi_crit, an ordinary point: read as
  // it stands, a falling profit there means the optimum lies inside and the
  // bracket stays, and anything else extends it to the end of the stem curve.
  constexpr bool bounded = cost_domain_ends_at_psi_crit<K>();
  bool ok_hi = false;
  bool rising = false;
  double hi = bound_b;
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
  if (!settled) {
    hi = bounded ? bound_b : bound_b - 1e-6 * width;
    f_hi = dprofit(hi, &ok_hi);
  }
  if (!settled &&
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
    // Rising wherever the gradient could be read: the dry end is the answer.
    hi = bound_b;
    f_hi = 1.0;
    ok_hi = true;
  }

  // No usable gradient at one end. No driver sweep reaches this -- both ends are
  // feasible on all 240 feasible golden-grid rows -- so treat it as a signal
  // rather than a routine path if it ever fires. Falling back to the search this
  // replaced means the change cannot make a previously-working case fail, which
  // is worth six lines on a solve plant runs millions of times.
  //
  // Two distinguishable faults reach here, and the tag separates them because
  // they are different things to have gone wrong. `!ok` is the supply/feasibility
  // report: no point near that end admits an informative gradient at all.
  // `ok && !isfinite` is a partial that came back non-finite where feasibility
  // said it would not -- checked here because there is nowhere downstream to put
  // that test. A NaN endpoint fails both pin comparisons below, since NaN
  // compares false against everything, so without this it reaches the root-find
  // on a bracket it cannot have whenever the other endpoint also fails its pin
  // test.
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

  // A CONSTRAINED optimum: profit is still climbing at one end of the feasible
  // interval, so dprofit never crosses zero inside it and the answer is that end.
  // Not a corner case -- 42 of the 240 feasible golden-grid rows land here (24
  // pinned wet, 18 pinned dry), all at psi_soil of 3 to 4, where profit is
  // negative and the best the leaf can do is barely transpire.
  //
  // ⚠️ **Return lo/hi, NOT bound_a/bound_b.** Returning the raw bound was tried
  // and is wrong, because of #31: below psi_upstream the profit algebra is built
  // on a negative conductance, so profit is DISCONTINUOUS across the feasibility
  // boundary rather than merely steep. Measured at psi_soil=4, vpd=2, 5 layers:
  // profit is -4.696 just inside the boundary and -6.136 at bound_a itself, a jump
  // of exactly the 1.44 that showed up on 22 golden rows as a profit DECREASE. The
  // #31 region is the reason the wet endpoint has to be stepped over rather than
  // evaluated, and the same reason it must not be returned.
  //
  // What lo buys, stated honestly: on a boundary-soil row the true argmax sits
  // essentially ON the feasibility boundary -- a fine scan puts it at 6.4e-07 and
  // 3.2e-09 of the bracket width from bound_a on the two worst rows -- so lo
  // resolves it to the step-in scale (~1e-6 of the width), not to
  // collar_root_tol. That is still ~100x tighter than the GSS_tol_abs it replaces,
  // and it beat golden section on profit at both of those rows. Bisecting for the
  // exact feasibility boundary would cost ~40 gradient evaluations on the hot path
  // to gain ~1e-07 of profit on a flat maximum, which is not a trade worth making.
  //
  // Worth knowing downstream: where a row crosses between interior and pinned the
  // argmax has a genuine KINK. That belongs to the constrained problem rather than
  // to this solver -- golden section has it too -- but it is a real
  // non-differentiable set in trait space, and it sits where a calibration is most
  // likely to wander.
  //
  // ⚠️ The two pin tests are not exhaustive, and the leftover case is a FAILURE
  // rather than a third pin. dprofit <= 0 at the wet end AND >= 0 at the dry end
  // means profit falls away from both ends into the interval: the interior
  // stationary point is a MINIMUM and the maximum sits at one of the two bounds.
  // Tagged SolverRefused, because reporting it as a pin -- to whichever bound the
  // ordering of the tests happens to reach first -- would hand back a plausible,
  // finite, wrong answer with a gradient that is genuinely non-zero in the wrong
  // direction.
  //
  // The gradients cannot say which of the two bounds it is; the profits can,
  // since f_lo <= 0 makes lo a local maximum and f_hi >= 0 makes hi one.
  // A pin's bound IS the answer, so it is solved to collar_root_tol before the
  // tests below read it (polish_bracket_bound says why the rows need this). Only
  // a side that would pin is touched, so an interior solve pays nothing; and a
  // polish whose re-evaluation fails keeps the bracket's answer.
  if (f_lo <= 0.0) {
    const double a = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
    wet_bound_root_ = a;
    if (g_min > 0.0) {
      // With a residual conductance the bound is a regular point -- the
      // stomatal flux is zero and the cuticle still conducts -- so the leaf
      // that declines to open sits exactly on it. That is what removes the
      // boundary-soil artefact of #183: a reported state set by the step
      // below rather than by the model.
      lo = a;
    } else if (a < hi) {
      // With nothing conducting at the bound itself, a shut leaf in the dark
      // has no finite ci there (ci = ca + R_d/g_c), so the reported state is
      // taken a step inside it. Stepped in on the bound's own scale, not the
      // width's: a step that is a fraction of the width moves with the DRY
      // bound, a parameter dependence the rows do not have -- and for root_c the
      // same size as the answer.
      //
      // ⚠️ ALWAYS FROM THE POLISHED ROOT, even where the polish left the bound
      // where the bracket had it. Keeping the bracket's own lower end there
      // instead put two different steps on either side of a driver value, and
      // the reported A and E jumped by 3e-6 of their range where a sweep
      // crossed between them.
      //
      // ⚠️ Do not start smaller than 1e-6. The rows read the flux at the stepped
      // point rather than at the root, an error proportional to the step (~1e-6
      // of profit per unit of a root trait, visible where the true derivative
      // is zero: a one-layer wet bound is the gravity balance). But at 1e-7 and
      // 1e-8 the forward model's stem-potential inversion reaches its rounding
      // floor, and a finite-difference referee found 56 rows wrong against 3.
      //
      // Only the no-flow test is repeated here, not the whole gradient: the sign
      // f_lo already has cannot flip across the ~1e-4 the bound moved at a pin,
      // and re-taking dprofit on every pinned row cost 1.5% of every solve.
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
  }
  if (f_hi >= 0.0 && hi == bound_b) {
    const double b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    // The dry end is feasible wherever it has been measured, and the sign f_hi
    // has cannot flip across the ~1e-4 the bound moved.
    if (b != bound_b && b > lo) {
      hi = b;
    }
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

  // Interior stationary point. f_lo > 0 > f_hi, so the bracket is valid and
  // uniroot_smooth's throw-on-bad-bracket cannot fire; passing the two endpoint
  // values it already has saves it re-evaluating them.
  operating_point_kind_ = OperatingPointKind::Interior;
  return util::uniroot_smooth(dprofit, lo, hi, f_lo, f_hi, collar_root_tol,
                              static_cast<size_t>(ci_niter));
}

template <Leaf::CostCurve K>
inline void Leaf::find_root_collar_psi_for(){
  // ⚠️ THE SAME PARAMETER CHECK THE STEM ROUTE MAKES. Without it, CF77 with an
  // unset lambda maximises a NaN objective here and returns a bracket property
  // that looks entirely like an operating point -- the exact failure
  // check_cost_parameters exists to refuse.
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
    // collar_at_dry_end / root_zero_E were consumed inside prepare_collar_solve; recover
    // them for the diagnostic message only if the profit check below fails.

    // Maximise carbon profit over the feasible collar-potential interval by
    // solving its first-order condition, dprofit/dpsi == 0, rather than by
    // searching the objective.
    const double opt_root_psi = maximise_profit_over_collar<K>(bound_a, bound_b);

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());
    seat_stem_at_dry_pin<K>();

    opt_root_psi_ = opt_root_psi;
    // ⚠️ `<K>`, NOT TF24. This read TF24's objective for every curve after the
    // collar route was templated -- so a SOX solve maximised SOX and then reported
    // TF24's profit at SOX's argmax, which is finite, plausible, and always a
    // little below TF24's own maximum. Two of these survived the templating pass.
    profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi);
    // Reported on the same axis as the single-layer solvers: the collar is the
    // free variable here, so the marginal cost carries the series-resistance
    // correction for the soil-to-collar path. For the SEATED curve, `<K>`.
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

// Evaluate the operating point at a given collar potential rather than
// optimising it (see header). A target outside the feasible interval is refused
// as infeasible (see profit_at_collar_psi), so a caller tracking the collar as a
// state rejects the step rather than reading an operating point at the bound.
inline void Leaf::find_root_collar_psi() {
  find_root_collar_psi_for<CostCurve::TF24>();
}


// Names in, model seated. The parsing half of set_model; the typed overload below
// is what actually assigns, so R and C++ cannot end up with two ideas of what
// seating means.
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
  // ⚠️ THE METHOD IS PARSED HERE AND VALIDATED NOWHERE ELSE, deliberately: what a
  // given method can and cannot do depends on the curve, the route and the energy
  // balance, and all three are settable after this call. So the name is checked now
  // and the CONFIGURATION is checked at solve time, exactly as the curve's own
  // constants are (see the typed overload).
  if (method != "exact" && method != "closed") {
    util::stop("unknown method \"" + method + "\": use \"exact\" for the "
               "first-order-condition root-find over the closed interval, or "
               "\"closed\" for the closed form. Available: exact, closed.");
  }
  set_model(static_cast<CostCurve>(k), route == "collar", method == "closed");
}

// The same, already parsed. For a C++ caller that has a `CostCurve` in hand --
// notably `gradient::route_seat`, which seats once per observation and would
// otherwise pay seven `std::string` constructions per gradient to look up a name
// it started from.
inline void Leaf::set_model(CostCurve curve, bool collar_route,
                           bool closed_method) {
  cost_curve_ = curve;
  route_is_collar_ = collar_route;
  method_is_closed_ = closed_method;
  // ⚠️ THE CURVE'S OWN CONSTANTS ARE NOT CHECKED HERE, deliberately: a caller may
  // reasonably seat the model and then set `CF77_lambda_`. They are checked at
  // solve time, where the ordering cannot be got wrong. The same goes for whether
  // the seated METHOD can serve the seated curve, route and energy-balance gate --
  // all of which are settable after this returns.
}


// The seated model, solved. Both routes in one expression, because the curve is
// the only thing `with_curve` has to choose and the route is a plain bool.
//
// ⚠️ THE CLOSED-FORM ARM RETURNS BEFORE THE ORIGINAL BODY, WHICH IS WHY THE EXACT
// PATH IS BIT-IDENTICAL BY CONSTRUCTION RATHER THAN BY TESTING. Everything below
// the early return is the function as it stood; the one thing the exact path
// gained is the bool write, which no golden baseline reads.
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


// The seated curve, like every other entry point that names one. Pinned to TF24
// these answered on a curve the leaf was not on: TF24_floor carries the price, so
// a floor leaf priced its collar with the price absent and reported a state that
// did not respond to it.
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
      // Operating point fully determined by feasibility handling (shutdown /
      // closed interval); profit_ is already set.
      return profit_;
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
      // The recorded collar is where this solve put the point, and a pinned
      // solve puts it at a POLISHED bound, which can lie up to ~1e-4 outside the
      // bracket prepare_collar_solve returns. Clamping it back into the bracket
      // would replay a different point, so the interval is widened to hold it.
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

// Post-prepare body of evaluate_root_collar_psi (see header). Kept as a separate
// entry point so callers that evaluate several collar potentials within one step
// (a centred finite difference) can run prepare_collar_solve once and reuse the
// soil-side caches across every profit eval.
//
// ⚠️ A COLLAR OUTSIDE [bound_a, bound_b] IS REFUSED, NOT PROJECTED. Projecting
// it onto the bound made the outputs move with the bound rather than with the
// collar asked for: a kink in every output of a tracked collar, and a derivative
// that answered a different question. The refusal is infeasible-class, because
// a tracked state reaches it at a well-formed call, so a consumer integrating
// the collar can reject the step rather than end the run.
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
    // Well defined off the optimum too: dC/dE is a property of the point, not of
    // its being stationary. At a prescribed point it simply is not equal to dA/dE.
    lambda_emergent_ = lambda_collar_for<K>();
    // Not an optimum of any kind: the collar potential came from the caller. On
    // this path a non-zero dprofit is the answer (TF24f's acclimation rate), not
    // a residual to be checked against zero.
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

// Exact d(profit)/d(opt_root_psi). profit(psi) = assim_colimited(ci) -
// hydraulic_cost_TF(psi_stem), with psi_stem = find_psi_stem_from_psi_root(psi)
// (smooth transport) and ci = psi_stem_to_ci(psi_stem, psi) (root-find).
// Chain rule:
//   dprofit/dpsi = A'(ci) dci/dpsi - C'(psi_stem) dpsi_stem/dpsi
// where dci/dpsi = (dci/dpsi_stem) dpsi_stem/dpsi + (dci/dpsi)|_explicit, and the
// dci/d* terms come from the implicit-function theorem on the residual
//   g(ci; psi_stem, psi) = A(ci) umol_to_mol - gc(psi_stem,psi) (ca-ci)/(atm kPa)
// with gc = const * transpiration(psi_stem,psi). A'/C' are obtained by forward
// AD; the gc partials use the stem table's own slope (stem_curve_integral_deriv); dpsi_stem/dpsi by a tight central difference on the smooth transport.
template <Leaf::CostCurve K>
inline double Leaf::dprofit_droot_collar_psi_for(double opt_root_psi, bool* feasible) {
  // Every transport evaluation below reads the supply path's per-solve caches, so
  // seat them on the current psi_soil_ here rather than depending on whatever the
  // caller's last solve left cached. Keeping this in the wrapper is what lets the
  // collar solve call the body directly and seat them once for the whole solve.
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
  // gstar_Pa used to be precomputed here and threaded into the AD replicas. The
  // kernels read gamma_ and umol_per_mol_to_Pa_ themselves, which is one fewer
  // place for the two sides to disagree.
  // Infeasible until the two exits below have been passed; see the header for why
  // the 0.0 they return must be distinguishable from a genuine stationary point.
  if (feasible != nullptr) {
    *feasible = false;
  }

  // Operating point in double. Past the continuity root with the stem at the
  // end of its curve the soil delivers more than the stem can carry at any
  // potential: no operating point, and so no gradient -- reachable by a caller
  // differencing across a dry pin, as the gradient's own stationarity check
  // does. Reported like the other no-gradient exits rather than thrown.
  E_from_Soil_to_Root_Collar(psi, supply_psi_soil());
  if (!(E_up_ / leaf_specific_conductance_max_ + stem_curve_integral(psi) <=
        stem_curve_limit() * (1.0 + 16.0 * std::numeric_limits<double>::epsilon()))) {
    return 0.0;
  }
  const double psi_stem = transpiration_to_psi_stem(E_up_, psi);
  // Shut down before the ci solve, not after. psi and psi_stem are both positive
  // magnitudes here, so psi >= psi_stem is the no-flow / reversed-gradient case --
  // the same condition set_leaf_states_rates_from_psi_stem treats as zero
  // transpiration. It has to be caught *here* because psi_stem_to_ci does not
  // return non-finite on a reversed pair, it throws: transpiration() refuses one
  // (#31). The isfinite check below was written to cover
  // shut-down but cannot see a thrown exception, so a dry patch killed a whole
  // plant run: reproduced on TF24f at 5 layers, theta = 0.005-0.03 with 1 m/yr
  // rainfall, at psi_stem = 1.23 against psi_upstream = 5.92 MPa.
  if (!std::isfinite(psi_stem) || psi >= psi_stem) {
    return 0.0;  // shut-down / infeasible: no informative gradient
  }

  // ⚠️ E1: SEAT THE TEMPERATURE PARAMETERS AT THIS CANDIDATE. Without this the
  // whole derivative below is evaluated at the AIR-temperature baseline, because
  // nothing on this path calls set_leaf_states_rates_from_psi_stem and
  // set_physiology gates its temperature cache on !use_energy_balance_. The
  // objective meanwhile is evaluated at Tleaf(E(psi)). So the solver was
  // root-finding the first-order condition of a different model from the one it
  // reported: measured before this landed, |dprofit| at the returned collar
  // reached 5.76 against the non-EB path's 5.6e-15, and the collar sat up to
  // 0.83 MPa from the true argmax.
  //
  // Same call order as set_leaf_states_rates_from_psi_stem uses, deliberately:
  // transpiration first, then the temperature update, then the ci solve, so the
  // ci root-find sees the same parameters the objective's does.
  const double E = transpiration(psi_stem, psi);
  double dT_dE = 0.0;
  double Tleaf_here = leaf_temp_;
  if (use_energy_balance_) {
    Tleaf_here = leaf_temp_from_E(E, &dT_dE);
    update_temperature_dependent_params(Tleaf_here);
    set_leaf_vpd(Tleaf_here);
  }

  // The stomatal flux the residual leaves, and the conductance the cuticle
  // adds (#182). Below the wet bound the stomatal flux reads as zero, the
  // corner the solve keeps (E_s >= 0); with nothing conducting at all -- the
  // g_min = 0 wet bound -- there is no derivative to take.
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

  // dpsi_stem/dpsi. Computed HERE, before the assigned-ci branch, so both
  // branches share ONE evaluation.
  //
  // ⚠️ Do not argue from inlining here. Factored into a separate member function
  // and inlined measure indistinguishably, with both binaries interleaved. Inline
  // is kept because sharing one evaluation between the two branches is less work
  // regardless, not because the compiler stops inlining a factored version.
  //
  // ⚠️ AND THE GATE-OFF COST OF THIS WHOLE CHANGE IS 3.1%, NOT THE ~1% FIRST
  // REPORTED. The first figure was taken at load average 8.9 with another
  // job's eight workers on the machine, where the arms straddled the noise.
  // Re-measured at load 5.9, six interleaved rounds, each arm stable to
  // 0.03 us: 3.132 us before, 3.230 us after. It is NOT the factoring and it is
  // NOT the out-of-line EB term (gated at the call site, so gate-off never
  // calls it); dprofit_at_collar_psi is out of line in both builds. Unexplained,
  // and worth explaining before this reaches plant, which runs this millions of
  // times. Object layout is EXCLUDED: sizeof(Leaf) is 2008 both before and
  // after, so the extra bool packed into existing padding. What remains
  // untested is whether the compensation branch's code enlarges
  // dprofit_at_collar_psi enough to change how it is scheduled -- it is out of
  // line in both builds, so this is about the body, not about the call.
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
    // psi_stem, already inverted at the top of this function from the same
    // E_up_ and psi. Read the slope there rather than inverting a second time.
    (void)E_psi_stem;
    dpsistem_dpsi = (1.0 / stem_curve_integral_deriv(psi_stem)) * dEpsistem_dpsi;
  } else {
    // Near a branch kink the analytic conductance returns NaN; fall back to a
    // central difference on the transport, as this path has always done.
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

  // A'(ci) and C'(psi_stem) by forward-mode AD of THE MODEL'S OWN algebra -- the
  // same kernels assim_colimited() and hydraulic_cost_TF() are instantiations of,
  // so these are derivatives of the function actually evaluated rather than of a
  // hand-kept mirror of it.
  AD ci_ad = ci;            xad::derivative(ci_ad) = 1.0;
  const double A_prime = xad::derivative(assim_colimited_kernel(ci_ad));
  const double C_prime = cost_deriv<K>(psi_stem, psi);

  // Stomatal-conductance supply coefficient gc and its partials. gc =
  // gc_const * (E - r*) + rho g_min, with E = transpiration(psi_stem, psi) =
  // conductance_max * (transp_from_psi(psi_stem) - transp_from_psi(psi)), so
  // the partials use the stem table's own slope. The residual flux reads
  // psi_stem only through the leaf temperature (the energy-balance term's) and
  // the collar through its ceiling, appended below; both are exact zeros at
  // g_min = 0.
  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double gc = gc_const * E_s + cuticle;
  const double dgc_dpsistem =
      gc_const * leaf_specific_conductance_max_ * stem_curve_integral_deriv(psi_stem);
  double dgc_dpsi =
      gc_const * leaf_specific_conductance_max_ * (-stem_curve_integral_deriv(psi));
  const double drstar_dpsi = residual_flux_dpsi_upstream(psi, Tleaf_here);
  dgc_dpsi = dgc_dpsi - gc_const * drstar_dpsi;

  // IFT on g(ci; psi_stem, psi): dci/dp = -(dg/dp)/(dg/dci).
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;      // dg/dci
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;
  const double dci_dpsi_expl = -(-dgc_dpsi * (ca_ - ci) * inv_atm) / g_ci;

  // dpsi_stem/dpsi: psi_stem = P(E_psi_stem) with
  //   E_psi_stem = E_up_(psi)/k_max + S(psi),
  // S = transpiration_from_psi and P its inverse, taken on that same table, and
  // E_up_(psi) the soil->collar uptake at collar suction psi. The collar variable
  // IS psi now, so there is no dr/dpsi = -1 factor to carry and the two terms add:
  //   dE_psi_stem/dpsi = E_up_'(psi)/k_max + S'(psi)
  //   dpsi_stem/dpsi   = P'(E_psi_stem) * dE_psi_stem/dpsi.
  // E_up_'(psi) is the analytic conductance (dE_from_soil_dpsi_collar), positive;
  // near a branch kink it returns NaN and we fall back to a central difference on
  // the transport.

  const double dci_dpsi = dci_dpsistem * dpsistem_dpsi + dci_dpsi_expl;
  // ⚠️ THE IDENTITY ARM IS TEXTUALLY UNCHANGED, so every difference objective's
  // arithmetic order survives operation for operation and the golden files stay
  // bit-identical. Only a non-identity link takes the general path.
  double base;
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    base = A_prime * dci_dpsi - C_prime * dpsistem_dpsi;
  } else {
    const double A = assim_colimited_kernel(ci);
    base = benefit_link_deriv<K>(A) * (A_prime * dci_dpsi) - C_prime * dpsistem_dpsi;
    // A log link has no derivative where the leaf is a net source (A <= 0,
    // its ci above ambient), which with the cuticle's CO2 pathway is an
    // ordinary point of the solve (#182). The objective is still defined
    // there; the derivative is reported as absent, never returned as NaN.
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
                    residual_cond_temp_deriv(psi, Tleaf_here, gc_const));
}


// dprofit/dpsi_stem, with the upstream potential held fixed.
//
// The same chain as above, with two simplifications that both follow from the
// decision variable being psi_stem itself rather than a collar potential mapped
// onto it: dpsi_stem/dx is 1, and every term carrying the upstream potential's
// motion is zero. What is left is
//
//     dprofit/dpsi_stem = A'(ci) dci/dpsi_stem - C'(psi_stem)
//
// and the cost curve enters through C' alone.
//
// ⚠️ THIS IS A DELIBERATE SECOND COPY OF `dprofit_at_collar_psi`'s ALGEBRA, and
// merging the two is a change to make knowingly or not at all. Setting
// `dpsistem_dpsi = 1` and `dgc_dpsi = 0` in the collar body does give this one
// mathematically -- but not TEXTUALLY, and the Identity arm's exact operation
// sequence is what keeps all four golden baselines bit-identical (multiplying by a
// literal 1.0 is exact and can still change FMA contraction; see the arm itself).
// The shared parts that carry no such constraint ARE factored out:
// `cost_deriv<K>`, `benefit_link_deriv<K>`, `temp_deriv` and
// `dprofit_energy_balance_term`. What is left duplicated is the arithmetic the
// golden files pin.
// dC/dpsi_stem for one cost curve: THE ONLY PLACE THE COST CURVE ENTERS a
// derivative, and shared by both routes. `dprofit_dpsi_stem` uses it directly;
// `dprofit_at_collar_psi` multiplies it by dpsi_stem/dpsi_collar, and adds
// `cost_deriv_upstream<K>` for a cost that also reads the collar. Factored out so
// there is one dispatch table rather than one per solver -- adding a curve must not
// mean remembering a second place.
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
    // lambda * E. E is kmax times the integral of the conductivity fraction over
    // [psi_upstream, psi_stem], so its derivative in psi_stem is kmax times that
    // fraction, which is what stem_curve_integral_deriv returns. No AD needed.
    C_prime = cf77_price() * leaf_specific_conductance_max_ *
              stem_curve_integral_deriv(psi_stem);
  } else if constexpr (K == CostCurve::JS22) {
    // C = gamma*(psi - psi_up)^2, so dC/dpsi = 2*gamma*(psi - psi_up). Analytic,
    // no AD and no table -- the simplest derivative of the four.
    C_prime = 2.0 * JS22_gamma * (psi_stem - psi_upstream);
  } else if constexpr (K == CostCurve::TF24_floor) {
    // C = TF24's cost + lambda_o*E, so dC/dpsi is the TF24 arm's expression plus
    // the CF77 arm's, each written in that arm's own operand order -- the AD pass
    // over the same kernel, and the same analytic `kmax*f`. See
    // `hydraulic_cost_TF24_floor` for why that makes both reductions bit-exact, FMA
    // included; the argument applies here term for term.
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    C_prime = xad::derivative(hydraulic_cost_TF_kernel(ps_ad)) +
              TF24_floor_lambda_o * leaf_specific_conductance_max_ *
                  stem_curve_integral_deriv(psi_stem);
  } else if constexpr (K == CostCurve::SOX) {
    // Under the log link the cost is `-log g`, so dC/dpsi is `-g'/g`. Both are
    // analytic; `g` is bounded below by zero only AT psi_crit, where the objective
    // is zero anyway and no optimiser returns it.
    C_prime = -sox_reduction_deriv(psi_stem) / sox_reduction(psi_stem);
  } else if constexpr (K == CostCurve::JW26) {
    // The same, and here `-g'/g` collapses: with `g = 1 - psi/psi_crit` it is
    // exactly `1/(psi_crit - psi)`.
    C_prime = 1.0 / (psi_crit - psi_stem);
  } else if constexpr (K == CostCurve::LeastCost) {
    // Under the log link the cost is `+log W`, where `W = E_mol + beta*Vcmax_mol`
    // is the weighted resource use of `leastcost_flux`. So dC/dpsi is `W'/W`, and
    // `W' = kg_to_mol_h2o * kmax * f(psi)` -- the capacity term is constant in
    // psi_stem at a fixed leaf temperature, so it survives only in the
    // denominator. With the energy balance on, Vcmax moves with Tleaf, and that
    // slope is taken in `dprofit_energy_balance_term` beside the benefit's own
    // temperature path, not here.
    //
    // ⚠️ THE SIGN IS THE OPPOSITE OF SOX's AND JW26's, and that is the one place a
    // ratio differs from a product. Theirs is `-g'/g` with `g` falling; this is
    // `+W'/W` with `W` rising. Both are costs that grow with tension; they get
    // there from different halves of the objective.
    C_prime = kg_to_mol_h2o * leaf_specific_conductance_max_ *
              stem_curve_integral_deriv(psi_stem) /
              leastcost_flux(psi_stem, psi_upstream);
  } else if constexpr (K == CostCurve::ProfitMax) {
    // C = [k(psi_s) - k(psi)]/k_span, so dC/dpsi = kmax*|f'|/k_span: the
    // hydraulic half alone. The thermal half TC(Tleaf) is constant in psi with
    // the energy balance off. With it on, Tleaf moves with E(psi) and the slope
    // needs dT/dE*dE/dpsi, which the benefit's own temperature path needs too,
    // so it is taken in dprofit_energy_balance_term beside that path, not here.
    const double f = proportion_of_conductivity(psi_stem);
    const double abs_fprime =
        f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
    C_prime = leaf_specific_conductance_max_ * abs_fprime / profitmax_k_span_;
  } else {
    // ⚠️ EVERY ARM IS EXPLICIT AND THE LAST ONE ASSERTS. This used to be a two-way
    // `if constexpr (TF24) {...} else {CF77}`, which would have routed a third cost
    // curve into the CF77 branch silently -- and `-Werror=switch` cannot see it,
    // because this is not a switch. The static_assert is what makes a FIFTH curve a
    // compile error instead of a plausible number.
    static_assert(K == CostCurve::CMax, "unhandled CostCurve in cost_deriv");
    // Linear in the ABSOLUTE potential by definition -- this IS the form Anderegg
    // et al. state, so no antiderivative and no cancellation is involved here.
    C_prime = CMax_a * psi_stem + CMax_b;
  }
  return C_prime;
}


// Whether a curve's cost reads the UPSTREAM potential as well as psi_stem. On the
// stem routes upstream is psi_soil and fixed, so this never matters there; on the
// collar route upstream IS the decision variable, and a cost that reads it moves
// with the collar directly as well as through psi_stem.
//
// ⚠️ A CURVE THAT SAYS FALSE HERE AND READS UPSTREAM ANYWAY SOLVES THE WRONG
// FIRST-ORDER CONDITION ON THE COLLAR ROUTE, and nothing else notices: the root-find
// still drives its own dprofit to ~1e-15, and every reported number is plausible.
// The symptom is `dprofit_droot_collar_psi` disagreeing with a central difference
// of `profit_` from `evaluate_root_collar_psi` -- which is what
// `test_collar_foc_matches_profit` asserts, for every curve. Every arm is explicit
// for the reason `cost_deriv`'s are.
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_reads_upstream() {
  if constexpr (K == CostCurve::CF77 || K == CostCurve::JS22 ||
                K == CostCurve::CMax || K == CostCurve::TF24_floor ||
                K == CostCurve::LeastCost) {
    // LeastCost: `log W` with `W = E_mol + beta*Vcmax_mol`, and E is measured
    // from the upstream potential -- the only log-link cost here that is.
    return true;
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::SOX ||
                      K == CostCurve::JW26 || K == CostCurve::ProfitMax,
                  "unhandled CostCurve in cost_reads_upstream");
    return false;
  }
}


// Whether a curve's own published form ends at psi_crit (#183).
//
// For most curves psi_crit is a labelled point on a vulnerability curve that
// continues past it. Their marginal cost of water diverges as the stem dries,
// so the optimum stays inside without a bound, and the dry end of the solve is
// the end of the stem curve itself (vulnerability_psi_end), where the stem
// carries all it ever can.
//
// Four curves are different, and for each the bound belongs to the curve, not
// to the hydraulics. Measured on psi_stem_optima.tsv, each sits AT psi_crit on
// many rows:
//
//   SOX, JW26  their published gain factor reaches zero at psi_crit. With A < 0
//              (every dark row) A*g is largest where g = 0, and past psi_crit
//              g < 0 and profit grows without bound, so psi_crit is where these
//              forms are defined to stop.
//   CF77       a constant price need not overtake A'(E) before E reaches the
//              stem's ceiling, so the supremum can lie at psi_stem -> infinity.
//   LeastCost  sits at psi_crit on most rows, dark and light.
//
// What these should do instead is open (a follow-up to #183); until then they
// keep the bound, a dry pin at psi_crit and the shut exit they had.
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
    // `d(log W)/du = (dW/du)/W`, and `W` reads upstream only through E, so
    // `dW/du = -kg_to_mol_h2o * kmax * G'(u)` -- CF77's arm with the kg-to-mol
    // factor in place of the price, over the same `W` `cost_deriv` divides by.
    // The capacity term is constant in both potentials and survives only in
    // that denominator.
    return -(kg_to_mol_h2o * leaf_specific_conductance_max_ *
             stem_curve_integral_deriv(psi_upstream)) /
           leastcost_flux(psi_stem, psi_upstream);
  }
}

template <Leaf::CostCurve K>
inline double Leaf::dprofit_dpsi_stem(double psi_stem, double psi_upstream,
                                      bool* feasible) {
  // ⚠️ EVERY CURVE, through the benefit link. This returns
  // `h'(A)*dA/dpsi - dC/dpsi` for whichever `h` the curve composes with
  // assimilation, so a product objective is not a special case -- it is the `Log`
  // link, and ProfitMax is the `Scaled` one. See `BenefitLink`.
  using AD = xad::fwd<double>::active_type;
  if (feasible != nullptr) {
    *feasible = false;
  }
  // No flow, so no informative gradient -- the same sentinel and the same
  // contract as the collar version. These are magnitudes, so upstream >= stem is
  // the reversed case.
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
  const double gc = gc_const * E_s + cuticle;
  const double dgc_dpsistem =
      gc_const * leaf_specific_conductance_max_ *
      stem_curve_integral_deriv(psi_stem);
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;

  // ⚠️ THE IDENTITY ARM IS WRITTEN OUT SEPARATELY AND THAT IS DELIBERATE. It is
  // textually what this line always was, so the four additive curves stay
  // bit-identical -- multiplying by a literal 1.0 would be exact but could still
  // change FMA contraction, and the golden files compare at the last bit.
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
  // dgc_dpsi = 0 and dpsistem_dpsi = 1: the term derives dE/dpsi from those two,
  // and with the upstream potential fixed that reduces to kmax * f(psi_stem).
  // The term applies the curve's link itself, so a non-identity curve's
  // temperature path takes the same h'(A) as its explicit path above, and
  // ProfitMax's thermal cost takes its own slope. An unscaled term added here
  // would be plausible and wrong, which is why the term is templated on K
  // rather than scaled by its callers.
  return base + dprofit_energy_balance_term<K>(
                    ci, gc, g_ci, inv_atm, gc_const, A_prime, dgc_dpsistem, 0.0,
                    1.0, dT_dE, Tleaf_here, E_s, r_star, 0.0,
                    residual_cond_temp_deriv(psi_upstream, Tleaf_here, gc_const));
}

// Evaluate at a prescribed psi_stem rather than optimising one. See the
// declaration for why a gradient here is a different derivative from a gradient
// at an optimum.
template <Leaf::CostCurve K>
inline double Leaf::evaluate_psi_stem(double target_psi_stem) {
  clear_collar_solve_state();

  if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }
  // Same three chains the optimiser uses, so a prescribed point and a solved one
  // cannot disagree about what a curve requires or which objective it is.
  check_cost_parameters<K>();
  // ⚠️ INCLUDING ProfitMax's NORMALISERS. The optimisers seed these and this did
  // not, so evaluating ProfitMax at a prescribed potential returned NaN while
  // solving for one worked -- a prescribed point disagreeing with a solved point
  // about whether the objective exists, which is exactly what the comment above
  // promises cannot happen.
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
    operating_point_kind_ = OperatingPointKind::Prescribed;
    return profit_;
  }

  // Refused, not projected onto [psi_wet, psi_dry], for the reason
  // profit_at_collar_psi gives: projected, the outputs moved with the bound
  // rather than with the potential asked for.
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
  operating_point_kind_ = OperatingPointKind::Prescribed;
  return profit_;
}

// The energy-balance correction to dprofit/dpsi, and zero when the gate is off.
//
// ⚠️ OUT OF LINE ON PURPOSE. Everything around it -- profit_psi_stem_TF,
// hydraulic_cost_TF, assim_colimited -- is fully inlined with no out-of-line
// symbol at all (hazard 5), and adding a branch inside such a body can change
// which expressions land in one inlined block and therefore change FMA
// contraction. That would move the GATE-OFF path, i.e. the golden file, for no
// reason anyone could read from the diff. Keeping the block behind its own
// symbol is what makes "bit-identical with the gate off" structural rather than
// lucky. Check with `nm -C test_golden | grep dprofit` before and after.
//
// THE DERIVATION. When the gate is on, psi reaches profit by two further routes
// beyond the two already accounted for:
//
//  * DIRECT:   psi -> psi_stem -> E -> Tleaf -> theta(Tleaf) -> A
//  * INDIRECT: psi -> psi_stem -> E -> Tleaf -> the ci residual -> ci -> A
//
// The second is the subtle one. psi_stem_to_ci root-finds
// g(ci; psi_stem, psi, T) = A(ci,T)*umol_to_mol - gc*(ca-ci)*inv_atm = 0, and the
// demand side reads the temperature-dependent members, so g gains an EXPLICIT T
// argument. Differentiating g = 0 totally in psi adds g_T * dT/dpsi to the
// existing terms, with g_T = A_T * umol_to_mol, so `dci/dpsi` gains
// `-(A_T * umol_to_mol * tau) / g_ci`.
//
// Adding that to the direct term `A_T * tau` and collecting gives
// `Delta = A_T * tau * (1 - A_prime*umol_to_mol/g_ci)`, which since
// `g_ci = A_prime*umol_to_mol + gc*inv_atm` is the same as
// `Delta = A_T * tau * (gc*inv_atm) / g_ci`.
//
// Two consequences worth keeping, because each is a free check on the algebra:
//
//  * the damping factor (gc*inv_atm)/g_ci lies strictly in (0,1). The direct
//    thermal effect on A is partly cancelled because ci re-equilibrates. If an
//    implementation ever produces |Delta| > |A_T*tau|, it is wrong.
//  * dE/dpsi > 0 and dT/dE < 0, so tau < 0. Above the thermal optimum A_T < 0
//    and therefore Delta > 0: the corrected condition pushes the optimum DRIER,
//    i.e. toward more transpiration, because the extra flux now buys evaporative
//    cooling back toward the optimum. That is the decoupling mechanism, and its
//    sign is the sharpest available test that this term is right.
//
// ⚠️ `dgc_dT` IS NON-ZERO, and placeholder algebra written while it was zero does
// NOT generalise -- switching the term on is a re-derivation, not a one-line
// change. A damping factor of the form `(gc*inv_atm - dgc_dT*(ca-ci)*inv_atm)/g_ci`
// puts the new term under `A_T` where the derivation puts it under `A_prime`. The
// two agree only at dgc_dT = 0, so such a form passes every test until the term
// becomes non-zero. Derived
// again, in full, below.
//
// g(ci; psi_stem, psi, T) = A(ci,T)*umol_to_mol - gc(T)*(ca-ci)*inv_atm = 0, so
//
//   g_T    = A_T*umol_to_mol - dgc_dT*(ca-ci)*inv_atm
//   dci/dT = -g_T/g_ci
//   dA/dT (total) = A_T + A_prime*dci/dT
//                 = [A_T*gc*inv_atm + A_prime*dgc_dT*(ca-ci)*inv_atm] / g_ci
//
// using g_ci = A_prime*umol_to_mol + gc*inv_atm to collect. At dgc_dT = 0 this is
// `A_T * (gc*inv_atm)/g_ci`, i.e. exactly what the old expression returned, so
// the prescribed-VPD behaviour is unchanged. The two checks in the paragraph
// above still hold for the first term; the second is new and has its own sign
// argument at the assignment.
// See the declaration. The step is 1e-3 degrees, matching dprofit_energy_balance_
// term's dA/dT difference, and the full restore is what makes calling this inside
// a derivative safe.
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

template <Leaf::CostCurve K>
inline double Leaf::dprofit_energy_balance_term(
    double ci, double gc, double g_ci, double inv_atm, double gc_const,
    double A_prime, double dgc_dpsistem, double dgc_dpsi, double dpsistem_dpsi,
    double dT_dE, double Tleaf, double E_s, double r_star, double drstar_dpsi,
    double dgc_dT_residual) {
  if (!use_energy_balance_ || dT_dE == 0.0) {
    return 0.0;
  }

  // dA/dTleaf by a CENTRAL DIFFERENCE over the temperature block only.
  //
  // Chosen over templating the Arrhenius block on the scalar type, and the
  // reason is the golden file rather than accuracy. Re-expressing
  // update_temperature_dependent_params as the T=double instantiation of a
  // templated kernel changes which expressions share an inlined body, and
  // That exact move was measured changing results even where the algebra
  // was identical. That risks the gate-off path to buy precision nobody can
  // observe: at h = 1e-3 K the truncation error is ~2e-9 relative against a
  // model whose own floor is set by psi_stem_to_ci's tolerance.
  //
  // ⚠️ h is FIXED and ABSOLUTE, not relative to T. A relative step would make
  // the estimator itself a function of psi, which is exactly the kind of thing
  // that puts a staircase back into the argmax (hazard 3).
  //
  // Differencing the model's own update + kernel means this cannot drift from
  // the function actually evaluated, and it picks up all eight temperature
  // parameters -- including R_d_ and electron_transport_, which a hand-derived
  // expression is precisely the sort of thing to forget.
  const double h = 1e-3;
  const double vcmax0 = vcmax_, jmax0 = jmax_, gamma0 = gamma_, ko0 = ko_,
               kc0 = kc_, Rd0 = R_d_, km0 = km_, J0 = electron_transport_;
  const double kimp0 = k_imp_eff_, krec0 = k_rec_eff_, func0 = functional_;
  const double pt0 = photo_temp_;
  // ⚠️ Tleaf, NOT leaf_temp_. On the energy-balance path set_physiology
  // reinterprets leaf_temp_ as AIR temperature (Tair_ = leaf_temp_), so
  // differencing around it evaluates dA/dT at the wrong point entirely. The
  // first version of this did exactly that; the resulting A_T was wrong enough
  // to flip the sign of dprofit at a bracket endpoint, and the collar solve
  // aborted with "parameters a and b do not bracket the root".
  const double T0 = Tleaf;

  update_temperature_dependent_params(T0 + h);
  const double A_up = assim_colimited(ci);
  const double vcmax_up = vcmax_;
  update_temperature_dependent_params(T0 - h);
  const double A_dn = assim_colimited(ci);
  const double vcmax_dn = vcmax_;

  // Restore by assignment rather than by a third update call: exact, and cheaper
  // than re-running the Arrhenius block.
  vcmax_ = vcmax0; jmax_ = jmax0; gamma_ = gamma0; ko_ = ko0; kc_ = kc0;
  R_d_ = Rd0; km_ = km0; electron_transport_ = J0;
  k_imp_eff_ = kimp0; k_rec_eff_ = krec0; functional_ = func0;
  photo_temp_ = pt0;

  const double A_T = (A_up - A_dn) / (2.0 * h);

  // dE/dpsi from the gc partials already in hand -- the stomatal part of gc is
  // gc_const * (E - r*), so dividing by gc_const recovers dE/dpsi less the
  // residual's own slope in the collar, which is put back. Tleaf reads the
  // TOTAL flux, so this is the slope the balance needs.
  const double dE_dpsi =
      (dgc_dpsistem * dpsistem_dpsi + dgc_dpsi) / gc_const + drstar_dpsi;
  const double tau = dT_dE * dE_dpsi;

  // dgc/dTleaf at fixed E. The stomatal part is K*(E - r*)/D(T) with E
  // hydraulically pinned, so its slope is -gc_s * D'(T)/D with D'(T) =
  // esat'(Tleaf) since e_air is fixed. NEGATIVE: a hotter leaf has a larger
  // deficit, so the same water flux implies a SMALLER conductance. Zero where
  // the deficit is on its floor, because the reported value no longer responds
  // to the input there. The residual's part -- a hotter cuticle loses more, which takes
  // conductance from the stomata and adds it to the cuticle -- is appended.
  const double gc_s = gc_const * E_s;
  const double dgc_dT =
      ((vpd_leaf_ > vpd_leaf_min)
           ? -gc_s * saturation_vapour_pressure_slope(Tleaf) / vpd_leaf_
           : 0.0) +
      dgc_dT_residual;

  // What the bracket holds is dA/dTleaf TOTAL, the implicit-function theorem
  // applied to the supply==demand residual with Tleaf moving: the direct A_T and
  // the deficit's dgc/dT combine into exactly this quotient, and tau turns it
  // into dA/dpsi along the temperature path. It is purely benefit-side. Two
  // costs read Tleaf, and both are taken below: ProfitMax's thermal cost, and
  // LeastCost's capacity term, through `vcmax_`. No identity-link cost does.
  //
  // ⚠️ THE IDENTITY ARM IS TEXTUALLY WHAT THIS RETURN ALWAYS WAS, so the five
  // additive curves keep their arithmetic operation for operation.
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    (void)vcmax_up;
    (void)vcmax_dn;
    return tau *
           (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) / g_ci;
  } else {
    // Under a link, J = h(A(psi)) - C(psi) with A(psi) the TOTAL dependence,
    // temperature path included, so this path takes the same single h'(A) the
    // call site put on the explicit dci/dpsi path. There is no h'' in a first
    // derivative: h'(A) is never itself differentiated (#124).
    const double dA_dpsi_T =
        tau * (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) /
        g_ci;
    // A at the operating point. The block was restored by assignment above, so
    // this is the value the call site's own `base` used.
    const double A = assim_colimited_kernel(ci);
    double term = benefit_link_deriv<K>(A) * dA_dpsi_T;
    if constexpr (K == CostCurve::ProfitMax) {
      // The thermal cost TC(Tleaf) reads the leaf temperature, so with the gate
      // on it moves with psi through E: dTC/dpsi = TC'(Tleaf) * dT/dE * dE/dpsi.
      // COST-SIDE, so it takes no h' factor -- the link sits on the benefit
      // alone -- and its slope is the one lambda_ProfitMax reports.
      term -= thermal_cost_temp_deriv(Tleaf) * tau;
    } else if constexpr (K == CostCurve::LeastCost) {
      // C = log W with W = E_mol + beta*Vcmax_mol (`leastcost_flux`), and
      // `vcmax_` is the value at Tleaf, so with the gate on the capacity term
      // moves with psi: d(log W)/dpsi gains beta*umol_to_mol*dVcmax/dT*tau/W.
      // `cost_deriv` holds only the E half, at fixed temperature. Cost-side, so
      // no h'. dVcmax/dT comes from the same two temperature-block evaluations
      // as A_T. W is rebuilt from gc, since gc = gc_const*E; the block is
      // restored, so `vcmax_` is the operating point's.
      const double dVcmax_dT = (vcmax_up - vcmax_dn) / (2.0 * h);
      // W reads the TOTAL flux: the stomatal part recovered from gc, plus the
      // residual, appended.
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

// Recompute the temperature-dependent photosynthetic parameters at a given leaf
// temperature. The arithmetic (and order) is exactly the inline block that used
// to live in set_physiology, so the non-PM path is bit-identical; extracting it
// lets the PM path recompute per operating-point Tleaf. electron_transport_ also
// depends on the per-call PPFD_ and is (re)computed here from the just-updated
// jmax_ -- on the non-PM cache-hit path set_physiology refreshes it separately.
//
// R_d comes from R_d_25 and the declining Q10, so it RISES with temperature.
// It used to be a fraction of vcmax_(T) and therefore fell above the thermal
// optimum, which was the wrong direction (#41).

// Every scalar update_temperature_dependent_params() below reads, in one place so
// the cache key and the computation cannot drift apart.
//
// ⚠️ KEEP THIS IMMEDIATELY ABOVE THAT FUNCTION AND IN STEP WITH IT. Adding an input
// there without adding it here reintroduces exactly the silent staleness this key
// exists to remove -- the fit still converges and the numbers stay plausible.
// `photo_temp_key_size` is the guard: the compiler rejects a mismatched list.
//
// Bit equality (`==`) is the right comparison. The question is "did any input
// change", not "did it change materially": an input that moved by one ULP produces
// a different response and must invalidate.
inline std::array<double, Leaf::photo_temp_key_size> Leaf::photo_temp_key() const {
  return {leaf_temp_, atm_o2_kpa_,
          vcmax_25, vcmax_ha_, vcmax_H_d_, vcmax_d_S_,
          jmax_25, jmax_ha_, jmax_H_d_, jmax_d_S_,
          gamma_25_, gamma_ha_,
          kc_25_, kc_ha_,
          ko_25_, ko_ha_,
          R_d_25, rd_q10_intercept_, rd_q10_slope_,
          // The thermal cost reaches jmax_ above, so it belongs in the key. A
          // bool widens to 0.0/1.0, which compares by == exactly like the rest.
          double(use_thermal_cost_), T50_, Tcrit_,
          // ⚠️ THE THERMAL-DAMAGE GATE IS HERE AND ITS INPUTS ARE NOT, because
          // set_physiology never serves this cache while the gate is on: A_, I_
          // and D_ are state a caller moves every step. The gate itself must be
          // here, or switching it off at unchanged drivers would be served the
          // gate-on block. Drop that bypass and every thermal input has to join
          // this list.
          double(use_thermal_damage_)};
}

inline void Leaf::update_temperature_dependent_params(double leaf_temp) {
  // FIRST: the pack-reading kernels take the block's temperature from here, and
  // two of them are called below.
  photo_temp_ = leaf_temp;
  if (use_thermal_damage_) {
    // The ATLS layer (see use_thermal_damage_). One shift moves both capacity
    // optima; the impairment's isokinetic temperature has its own, sharing only
    // the constitutive offset and the acclimation state.
    // A itself, for the reason delta_impair gives.
    const double A = A_;
    const double delta_topt = topt_offset_ + dTopt_max_ * A / (A + K_A_);
    const double vcmax_dS =
        d_S_shifted(vcmax_ha_, vcmax_H_d_, vcmax_d_S_, delta_topt);
    const double jmax_dS = d_S_shifted(jmax_ha_, jmax_H_d_, jmax_d_S_, delta_topt);
    vcmax_ = peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_dS);
    jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_dS);
    // Only the functional area photosynthesises, so BOTH capacities are scaled
    // by F. The instantaneous, reversible deactivation is inside the peaked
    // Arrhenius already; impairment and damage are lasting losses on top of it,
    // so this is not a second count of it.
    //
    // F is held above zero by a SOFT floor of width damage_surv_floor, f:
    // F + f ln(1 + e^(-F/f)), the softplus written so that F itself is kept
    // exactly wherever the correction is below half an ulp of F (F > ~3.4e-5). A
    // min/max floor put a corner in every output at F = f; sustained lethal
    // heat can drive F that low, and nothing below it is refused. I_ and D_
    // are refused below zero, so F <= 1 needs no cap.
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
  // Sicangco et al. (2026) Eqn 12: irreversible PSII damage takes electron
  // transport down with it, so the peaked Arrhenius Jmax is scaled by (1 - TC).
  // Branched rather than multiplied by a gate-off zero, so the gate-off
  // arithmetic is untouched.
  if (use_thermal_cost_) {
    jmax_ *= 1.0 - thermal_cost_at(leaf_temp);
  }
  gamma_ = arrh_curve(gamma_ha_, gamma_25_, leaf_temp);
  ko_ = arrh_curve(ko_ha_, ko_25_, leaf_temp);
  kc_ = arrh_curve(kc_ha_, kc_25_, leaf_temp);
  // Respiration RISES with temperature, on the declining Q10 above. ⚠️ It is
  // exactly R_d_25 at 25 C and larger above, so any check taken at 25 C alone is
  // blind to this whole response -- test_rd_temperature_response is what covers it.
  //
  // Fail rather than fall back: R_d_25 is a trait, and an unset one is a caller
  // error, not a cue to derive something from vcmax.
  if (!std::isfinite(R_d_25) || R_d_25 < 0.0) {
    util::stop("R_d_25 must be a finite, non-negative dark respiration at 25 C; "
               "got " + util::format_double(R_d_25));
  }
  const double q10 =
      rd_q10_intercept_ - rd_q10_slope_ * (leaf_temp + 25.0) / 2.0;
  R_d_ = R_d_25 * std::pow(q10, (leaf_temp - 25.0) / 10.0);
  if (use_thermal_damage_) {
    // The layer's maintenance costs, added to the respiration above rather than
    // folded into it, so the gate-off arithmetic is untouched. Standing costs for
    // the acclimation held and for the repair capacity and protection kept, paid
    // whether or not they are used; an activity cost for the repair done.
    //
    // Each is read as itself. A_ and I_ are non-negative states (see
    // delta_impair). Protection is SIGNED: the default k_i_ equals
    // k_i_baseline_, so a max(0, .) put the defaults exactly on a corner, and a
    // leaf less protected than the baseline saves the cost one more protected
    // pays.
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
// R*ln(Ea/(H_d - Ea)), for the d_S that puts it delta_topt higher.
//
// ⚠️ THE OPTIMUM EXISTS ONLY FOR 0 < Ea < H_d, and both are settable members, so
// that is refused here rather than left to return a NaN capacity.
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

// ⚠️ A NEGATIVE k_i_EA_ IS REFUSED: impairment would then slow as the leaf
// heats. Zero is a temperature-independent rate, and allowed. A NaN passes
// through rather than throwing, so an unset input poisons the outputs where a
// test can see it.
inline void Leaf::check_damage_pars() const {
  if (k_i_Ea_ < 0.0) {
    util::stop("the thermal-damage layer needs k_i_Ea_ >= 0, so that "
               "impairment does not slow with temperature; got k_i_Ea_ = " +
               util::format_double(k_i_Ea_));
  }
}

// ⚠️ ITS OWN PARAMETERS, AND NOTHING OF THE JMAX CURVE. Do not read jmax_H_d_ or
// jmax_d_S_ here: that would make recalibrating the instantaneous photosynthetic
// response silently move the impairment rate.
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

// One forward-mode pass per requested input, each setting a dual number's
// direction in that parameter (or initial state) and running the same RK4 the
// double trajectory runs. Exact for the discrete map, at the sub-step count the
// double run uses.
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
// the derivative of the Tetens form above, so its numerator is the product of
// that form's two coefficients.
//
// ⚠️ KEEP IT THE EXACT DERIVATIVE, NOT FAO-56's QUOTED 4098. That constant rounds
// 17.27 * 237.3 = 4098.171, 4.2e-5 low, and the slope reaches the energy-balance
// first-order condition through dgc/dT = -gc * Delta(Tleaf) / D_leaf. Rounded, the
// analytic dprofit/dpsi_stem missed a central difference of the profit by 2e-6 to
// 2e-5 relative at every step size; exact, it agrees to the difference's own floor.
inline double Leaf::saturation_vapour_pressure_slope(double temp) const {
  return (17.27 * 237.3) * saturation_vapour_pressure(temp) /
         ((temp + 237.3) * (temp + 237.3));
}

// The deficit that drives diffusion out of the leaf (#7).
//
// Fick's law needs the deficit between the SUB-STOMATAL cavity, saturated at the
// LEAF temperature, and the air. The driver `atm_vpd_` is the deficit at AIR
// temperature, so the two agree only when the leaf is at air temperature:
//
//   e_air     = esat(Tair) - atm_vpd
//   vpd_leaf  = esat(Tleaf) - e_air = atm_vpd + [esat(Tleaf) - esat(Tair)]
//
// ⚠️ WRITTEN AS `atm_vpd_ + (esat(T) - esat(Tair_))` FOR AN ARITHMETIC REASON.
// Off the energy-balance path the leaf IS at air temperature, so the bracket is
// a value minus itself: exactly +0.0, and `x + 0.0 == x` for every finite x. The
// prescribed-temperature path is therefore bit-identical, which is what keeps
// the golden file untouched. Computing it as `esat(T) - esat(Tair) + atm_vpd` in
// a different association would not guarantee that.
//
// ⚠️ AND IT CAN GO NON-POSITIVE. A leaf transpiring hard enough to sit well below
// air temperature has esat(Tleaf) < e_air, i.e. condensation rather than
// evaporation, which Fick's law in this direction does not describe -- and a
// deficit passing through zero is an infinite conductance. Clamped to
// `vpd_leaf_min` rather than allowed to invert the flux: the model would
// otherwise report a NEGATIVE stomatal conductance for a positive transpiration.
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

// --- the residual (cuticular) conductance (#182) ------------------------------

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

// kmax (G(inf) - G(psi_up)): the most the stem can ever carry from psi_up,
// which it approaches as its own potential runs away.
inline double Leaf::supply_ceiling(double psi_upstream) const {
  return leaf_specific_conductance_max_ *
         (stem_curve_limit() -
          stem_curve_integral(psi_upstream, "Leaf::supply_ceiling"));
}

// r* = r (1 + (r/E_inf)^p)^(-1/p), the smooth minimum of the demand and the
// ceiling (see residual_cap_power). Exactly zero at r = 0, and the ceiling
// itself where the ceiling has closed.
// Past q = residual_cap_q_limit the cap is its limit phi E_inf to within
// q^-p / p, below double precision, and the power would overflow; the limit is
// taken there, on every path, so none of them can form inf * 0.
inline constexpr double residual_cap_q_limit = 1e4;
// Below q^p = 2^-54 the sum 1 + q^p rounds to 1 exactly, so the cap returns r
// bit for bit and its slope in the demand is 1. That is nearly every solve --
// a cuticle loses ~1e-3 of what the stem can carry -- and skipping the two
// pow() calls there is most of what a residual conductance costs the solve.
// Exact, not an approximation: the arithmetic below would give the same bits.
inline constexpr double residual_cap_inert_qp = 0x1p-54;

// (1 + q^p)^(-1/p) at p = 8, by squarings and square roots: p is a power of
// two, so q^8 is three squarings and the eighth root three square roots, each
// correctly rounded. Two pow() calls here were the largest single cost a
// residual conductance added to the solve (a seventh of it, profiled). One
// expression for the double path and the tape, so the two cannot drift.
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

// The intercellular CO2 at which diffusive supply meets demand,
// `A(ci) umol_to_mol = gc (ca - ci) / P`, root-found on [Gamma* cP, hi]. The target is negative at Gamma* (A = -R_d
// there) and rises with ci. With A(ca) >= 0 it is non-negative at ca, the
// bracket every solve used to take. With A(ca) < 0 the root is ABOVE ca -- the
// leaf is a net source and CO2 leaves it -- and hi is placed where the supply
// term exceeds the demand's shortfall, 2 |A(ca)| umol_to_mol P / gc above ca,
// since A only rises past ca. With no conductance at all (g_min = 0 and the
// stomata shut) the residual degenerates to A(ci) = 0: its root where one
// exists, which is the limit g_min -> 0 of the shut leaf, and otherwise the
// leaf nets its saturated assimilation at ci = +inf.
//
// ⚠️ **The 1e-13 is load-bearing and was chosen by measurement, not taste. It
// sets the floor of what every reported output of this model MEANS.** This is
// the model's dominant amplifier: a last-bit change anywhere upstream shifts
// which point TOMS748 lands on inside its tolerance band, and that surfaces in
// ci, assim, gc and profit. Measured by perturbing the model identically at each
// tolerance and reading the golden grid:
//
//   ci tol   worst diff vs a converged (1e-15) solve   us/solve (indicative)
//   1e-7                       1.82e-07                      2.655
//   1e-8                       3.15e-08                      2.640
//   1e-10                      5.41e-10                      2.695
//   1e-11                      5.41e-10                      2.768
//   1e-13                      7.48e-13                      2.885   <- here
//   1e-15                      0                             2.945
//
// The plateaus (1e-8 = 1e-9, 1e-10 = 1e-11) are TOMS748's discrete iterations,
// so tightening *between* them buys nothing and still costs. The timing column
// is only indicative: those binaries were built in one loop in a scratch tree,
// and code-layout luck moves this benchmark by ~2%, the size of the effect. If
// you re-derive the curve, build every point from one tree and interleave.
//
// ⚠️ **What chose 1e-13 over 1e-10 is a consumer that differences the solve
// (#172).** Where TOMS748 stops inside its band depends on the state, so at
// 1e-10 the reported outputs are noisy along any smooth path through the
// drivers -- up to 3.7e-11 relative over a grid of plant heights and light
// levels, 1e-15 at others. plant takes finite differences of the leaf at steps
// of 1e-7, where that noise is a 1e-4 error. At 1e-13 the worst over the same
// grid is 2.3e-14, and 1e-12 still leaves 1.2e-12. Do not loosen it without
// re-reading the guide's rounding-versus-bug magnitudes, which this figure sets.
//
// ⚠️ Not the same knob as `ci_abs_tol` (the settable control, default 1e-3).
// `ci_abs_tol` is read in exactly one place,
// solve_medlyn_ci_numerical, i.e. the empirical Medlyn-conductance route. It
// reaches NEITHER family of optimality solver -- both go through this function
// and so through the 1e-13 above. A caller tightening it gets no extra precision
// anywhere in the optimality model.
inline double Leaf::solve_ci(double gc) {
  // Propagate non-finite inputs as NA rather than entering the solver: the
  // NA-in -> NA-out contract (see test-leaf.r "Basic functions").
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
  // The target at ca is A(ca) exactly (the supply term is gc * 0), and it is
  // already in hand: handed to the root-finder rather than evaluated again,
  // which saves one assimilation kernel per ci solve, the cost the bracket
  // choice above added. Bit-identical: x - 0.0 is x.
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
    return util::uniroot_smooth(target, lo, hi, target(lo), f_hi, ci_root_tol,
                                ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("ci_solve", "solve_ci failed: " + std::string(e.what()) +
               "; min=" + util::to_string(lo) +
               "; max=" + util::to_string(hi) +
               "; gc=" + util::to_string(gc));
  }
}

template <class S>
inline S Leaf::stem_integral_limit_at(const leaf_pars<S>& pars) const {
  using odelia::util::to_passive;
  const double P50_0 = to_passive(pars[par_stem_P50]);
  const double c_0 = to_passive(pars[par_stem_c]);
  const double b = weibull_b_from_P50(P50_0, c_0);
  const VulnerabilityIntegralDerivatives d =
      cumulative_vulnerability_integral_limit_derivatives(b, c_0);
  const TraitRows rows = rows_in_P50(d.db, d.dc, b, c_0);
  return S(d.value) + S(rows.dP50) * (pars[par_stem_P50] - S(P50_0)) +
         S(rows.dc) * (pars[par_stem_c] - S(c_0));
}

template <class S>
inline S Leaf::residual_flux_at(const S& collar, const leaf_pars<S>& pars,
                                const SupplyDraw<S>& draw) const {
  // The demand: g_min at the seated leaf temperature, at the seated deficit.
  // Both are the operating point's and passive; the trait carries the row.
  const double per_g = g_min_factor(Tleaf_) * vpd_leaf_ / atm_kpa_ / kg_to_mol_h2o;
  const S r = pars[par_g_min] * S(per_g);
  const S E_inf = pars[par_kmax] * (stem_integral_limit_at<S>(pars) -
                                    stem_integral_at<S>(collar, pars));
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

// Explicit leaf energy balance (#523, #97):
// Tleaf = Tair + (Rn - lambda*E) / (rho*cp/ra + g_rad).
// E is the hydraulically-pinned transpiration (kg H2O m^-2 s^-1), so lambda*E is
// the latent heat flux (W m^-2) and (Rn - lambda*E) is what sensible heat and the
// leaf's own extra longwave emission carry away between them.
//
// ⚠️ THE ONE PLACE Tleaf AND dTleaf/dE ARE COMPUTED. Every energy-balance route
// (the collar and stem first-order conditions, the thermal cost, the reported
// Tleaf_) calls this; a second copy of the expression is how a radiative term
// gets dropped from a derivative while the value keeps it.
inline double Leaf::leaf_temp_from_E(double E, double* dT_dE) const {
  // Sensible heat plus the linearised longwave, in one conductance. Written as a
  // division by the SUM rather than a multiplication by ra_/vol_heat_cap_air so
  // the radiative term cannot be dropped by editing one factor.
  const double g_total = vol_heat_cap_air / ra_ + g_rad_;
  const double Tleaf = Tair_ + (Rn_ - latent_heat_vap * E) / g_total;
  if (dT_dE != nullptr) {
    // The balance is linear in E, so the slope is a constant and negative: more
    // transpiration, more latent heat, cooler leaf.
    *dT_dE = -latent_heat_vap / g_total;
  }
  return Tleaf;
}

// ⚠️ NO CLAMP, AND THESE ARE NOT ONE. Tleaf used to be clamped to [-40, 70] C,
// and a heatwave leaf shut in still air (Tair 45, PPFD 2000, wind 0.3 m/s, d 0.2
// m) solved AT 70.00: a clamp artefact reported as the leaf's temperature.
// Nothing in the leaf needs that range. The two limits here are where it stops
// being evaluable: R_d's Q10, intercept - slope*(T + 25)/2, must be positive for
// its power to be real (about 119 C at the defaults), and the Tetens form has a
// pole at -237.3 C. Q10 is linear in T, so testing both ends tests the range.
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
  // ⚠️ RECORDED, BECAUSE `set_traits` REBUILDS AT THE MEMBER. Without this the
  // object holds a curve at one resolution and re-traits it at another, so a
  // re-traited leaf is not the leaf a fresh one at those traits would be -- which
  // is the guarantee set_traits rests on, and it fails silently: every number
  // stays plausible and only the curve under the solve changes. A caller that
  // builds the two curves at DIFFERENT resolutions gets the last one recorded,
  // which is the same one member both rebuilds have always read.
  vulnerability_curve_ncontrol = resolution;
  // The root's table is the same table wherever the two c agree (the defaults),
  // so it is shared rather than built twice.
  stem_curve_ = vulnerability_table(stem_c, resolution, roots_.vulnerability_);
  ++stem_curve_builds_;
}

// --- the stem curve, as the four operations performed on it ------------------

// Off the hot path: only reached when a read is about to be refused, so the
// string building costs nothing in production.
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

// ⚠️ The comparisons are written so that a NaN fails every one and falls
// through to the read, which returns NaN: plant documents a
// profit_psi_stem_TF(NA, .) -> NA contract built on exactly this.
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
  // G rises to G(inf), which it reaches to rounding at the end of the curve,
  // so a potential gives w exactly when 0 <= w <= G(inf), and G(inf) itself is
  // the end of the curve. Past G(inf) no stem potential carries the flow.
  //
  // ⚠️ "PAST" MEANS PAST ROUNDING. A flux the model keeps below the stem's
  // ceiling -- the residual flux under its soft cap, a continuity root taken
  // from its feasible side -- arrives here as w = E/kmax + G(upstream), and the
  // sum can round a few ulps past G(inf) when the ceiling is nearly all used.
  // That is G(inf), and the end of the curve is its potential.
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
  // ⚠️ ONE SPELLING. dpsi/dG = 1/G'(psi) is an identity, so the reciprocal of
  // the forward slope at the psi this inverts to IS the answer -- and taking it
  // any other way is a second definition of one number.
  return 1.0 / stem_curve_integral_deriv(stem_curve_integral_inverse(w));
}

inline void Leaf::perturb_stem_P50(double stem_P50_new) {
  check_psi_magnitudes(stem_P50_new, stem_c, roots_.root_P50, roots_.root_c);
  stem_P50 = stem_P50_new;
  stem_b = weibull_b_from_P50(stem_P50_new, stem_c);
  psi_crit = weibull_psi_at_fraction(stem_b, stem_c, k_crit_fraction);
  // ⚠️ The transpiration memo is keyed on (psi_stem, psi_upstream) and NOT on the
  // curve, so a perturbation that leaves the potentials alone -- which is every
  // perturbation a gradient makes -- would hit an entry computed from the old
  // stem_b and return it. That is the bug this line prevents, and it is invisible:
  // the returned value is a plausible transpiration.
  transpiration_cached_ = false;
}

// Direct integration of the xylem vulnerability curve, as an independent check
// on the pre-integrated table that transpiration() uses on the hot path. Not
// used in production. Uses leaf/quadrature.hpp rather than plant's compiled QAG,
// so values are convergent-but-not-bit-identical to plant's; see that header.
inline double Leaf::transpiration_full_integration(double psi_stem, double psi_upstream) {
  const auto f = [&](double psi) -> double {
    return proportion_of_conductivity(psi);
  };
  return leaf_specific_conductance_max_ *
         quadrature::adaptive_simpson(f, psi_upstream, psi_stem,
                                      integration_tol_);
}

// transpiration()'s refusal, kept out of line: the string building inlined into
// that hot path cost the solve measurably, where a call to a cold function does
// not.
[[noreturn, gnu::cold, gnu::noinline]] inline void
stop_reversed_flow(double psi_stem, double psi_upstream) {
  util::stop("Leaf::transpiration: psi_stem = " + util::format_double(psi_stem) +
             " MPa is wetter than psi_upstream = " +
             util::format_double(psi_upstream) +
             " MPa, which would move water from the leaf to the soil; "
             "psi_stem must be at least psi_upstream");
}

//calculates supply-side transpiration from psi_stem and opt_root_psi_, returns kg h20 s^-1 m^-2 LA
// SIGN: psi_stem and psi_upstream are POSITIVE magnitudes here (passed straight
// to the table). Contrast transpiration_to_psi_stem below. See the sign-
// conventions block above.
inline double Leaf::transpiration(double psi_stem, double psi_upstream) {

  // 1-entry memo: identical (psi_stem, psi_upstream) is requested several times
  // per profit evaluation; return the cached value (bit-identical) to skip the
  // redundant table lookups. Cache invalidated in set_physiology.
  if (transpiration_cached_ &&
      psi_stem == transpiration_cache_psi_stem_ &&
      psi_upstream == transpiration_cache_psi_upstream_) {
    return transpiration_cache_value_;
  }

  // integration of proportion_of_conductivity over [opt_root_psi_, psi_stem]
  const double E = leaf_specific_conductance_max_ *
    (stem_curve_integral(psi_stem, "Leaf::transpiration, at psi_stem") -
     stem_curve_integral(psi_upstream, "Leaf::transpiration, at psi_upstream"));
  // return (transpiration_full_integration(psi_stem));

  // A stem wetter than its upstream end would move water from the leaf to the
  // soil. The integral above is then negative, and stom_cond_CO2 and
  // psi_stem_to_ci would hand back a negative conductance and a ci solved on it,
  // as finite numbers. No solve asks for this flux (#31): the operating-point
  // evaluator reads psi_upstream >= psi_stem as no flow without calling here, and
  // a cost that re-reads the flux goes through operating_point_flux. So a reversed
  // pair here is a caller error. Equality is the shut-down state and stays legal;
  // NaN falls through to a non-finite E. Tested after the lookups so a negative
  // potential still gets the table's domain message, and before the memo is
  // written so the memo never holds a reversed pair.
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

// converts a known transpiration to its corresponding psi_stem; both potentials
// are positive magnitudes, as in transpiration() above, of which this is the
// inverse.
inline double Leaf::transpiration_to_psi_stem(double transpiration_, double psi_upstream) {
  // integration of proportion_of_conductivity over [opt_root_psi_, psi_stem].
  // psi_upstream is a positive magnitude, same as in transpiration() -- #7 noted
  // that this function read eval(-psi_upstream) while its forward direction read
  // eval(psi_upstream), and that this was consistent because the two were called
  // with opposite-sign arguments. Under one representation they are not, and the
  // negation is deleted (#25).

  // Past the end of the curve G is G(inf) to rounding, so every potential there
  // maps to one value and G^-1 cannot say which. A zero flux has the stem AT the
  // upstream potential, which is the answer the inverse would give if it could;
  // any positive flux there needs more than G(inf) and is infeasible.
  if (transpiration_ == 0.0 &&
      psi_upstream >= vulnerability_psi_end(stem_b, stem_c)) {
    return psi_upstream;
  }

  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +
    stem_curve_integral(psi_upstream, "Leaf::transpiration_to_psi_stem, at psi_upstream");

  // The lookup that killed plant#576, and the reason the inverse names its units:
  // E_psi_stem below the table's lower end means the demanded flux is NEGATIVE,
  // i.e. the collar cannot supply it and the stem potential that would carry it is
  // wetter than saturation. There is no such potential, so widening the domain is
  // not the fix -- the caller should not have asked. Naming the argument E/K_max
  // rather than a bare `u` is what makes that readable from the message.
  return stem_curve_integral_inverse(
      E_psi_stem, "Leaf::transpiration_to_psi_stem, inverting for psi_stem");
}

// The TOTAL conductance to CO2 at a candidate pair, mol CO2 m^-2 LA s^-1: the
// stomata's, from the flux the residual leaves them, plus the cuticle's. Reads
// the deficit and the residual conductance as seated (set_leaf_vpd), so a
// caller on the energy-balance path seats the leaf temperature first, as the
// evaluator and the first-order conditions do. Below the wet bound the
// stomatal flux reads as zero: the corner E_s >= 0, counted.
inline double Leaf:: stom_cond_CO2(double psi_stem, double psi_upstream) {
  const double E = transpiration(psi_stem, psi_upstream);
  double E_s = E - residual_flux(psi_upstream, E);
  if (E_s < 0.0) {
    clamps.note(CLAMP_STOMATAL_FLUX);
    E_s = 0.0;
  }
  return atm_kpa_ * E_s * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_ +
         cuticular_cond_CO2();
}


// biochemical photosynthesis model equations
// Electron transport: the non-rectangular hyperbola in the electron flux the
// ABSORBED light can drive, phi_J * APAR with APAR = alpha_PAR_ * PPFD_ (PPFD_ is
// the intercepted photon flux, umol m^-2 s^-1). electron_transport_at is the same
// expression on the parameter pack, in the same order.
inline double Leaf::electron_transport() {
  const double I = phi_J * alpha_PAR_ * PPFD_;
  double electron_transport_ = (I + jmax_ - sqrt(pow(I + jmax_, 2) -
  4 * curv_fact_elec_trans * I * jmax_)) / (2 * curv_fact_elec_trans);
    return electron_transport_;
}

//calculate the rubisco-limited assimilation rate, returns umol m^-2 s^-1
// The scalar-generic cores. Arithmetic and call structure are IDENTICAL to the
// `double` bodies these replaced -- same association, same `pow(s, 2)`, same
// two-call structure in the colimited form -- so the `double` instantiation is not
// meant to move results, and the golden file is the check on that. What moves is
// the DERIVATIVE, which now comes from this code instead of from the drifted
// replica.
// The capacity's scalar U is independent of the kernel's T on purpose: on the
// value and tangent paths the capacity is four cached doubles, and a double
// times a tangent costs one multiply where a tangent times a tangent costs
// three. Wrapping the constants in T doubled the kernel's cost per solve
// (measured 5% to 10% of a collar solve).
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

// converts psi stem to ci, used to find ci which makes A(ci) = gc(ca - ci)
//
// #486: this target (assim_colimited demand minus the linear gc supply) is
// smooth and strictly monotone in ci -- no singularity at the bracket ends
// (Ar, Ae vanish linearly at gamma* so sqrt(disc) is linear, not singular) and
// its only sharp feature is the colimitation elbow far below the operating
// root. So it is a well-behaved case for a superlinear bracketing solver:
// TOMS748 reaches the same root in ~9 evals vs bisection's ~29. The hydraulic
// find_root_psi path keeps its own solver (see the warning on
// util::uniroot_smooth). The bracket and the degenerate case are solve_ci's.
inline double Leaf::psi_stem_to_ci(double psi_stem, double psi_upstream) {
  return ci_ = solve_ci(stom_cond_CO2(psi_stem, psi_upstream));
}

// given psi_stem, find assimilation, transpiration and stomatal conductance to CO2
//
// ONE PASS, NO BRANCH ON THE FLOW (#182). Transpiration is the hydraulic
// supply at the candidate potential, zero where the pair is reversed (#31);
// the leaf temperature follows it; the residual flux the cuticle loses at
// that temperature is taken off it; what is left drives the stomata, the
// cuticle adds its own conductance to CO2, and ci is the root of supply ==
// demand on the total. A shut leaf (E_s = 0) is the same expressions with the
// stomatal term zero, so its outputs cannot disagree with the open leaf's next
// to it; with nothing conducting at all, which is only g_min = 0, solve_ci
// takes the limit.
//
// ⚠️ AN `assim_max_ < 0` EARLY EXIT USED TO SIT HERE, AND IT ZEROED THE WATER.
// Transpiration on this path does not ask whether there is carbon to be had --
// water moves down a potential gradient whatever photosynthesis is doing.
// Zeroing it made two leaves at the SAME operating point disagree about whether
// water was moving purely because one of them had respiration switched on.
//
// ⚠️ `Tleaf_` IS WRITTEN HERE ON EVERY PATH, which is hazard 8 applied to the
// one output whose value is a driver on one path and a solved quantity on the
// other: on the PM path a leaf that transpires nothing is the HOTTEST one, so
// leaving it to a transpiring candidate would report the previous candidate's
// temperature for exactly the operating points where the answer is most
// extreme (#106). The temperature block and the deficit are seated at the
// STORED `Tleaf_`, so the temperature the parameters were derived at, the one
// the deficit was computed at, and the one reported cannot be three different
// numbers -- and it costs no extra `leaf_temp_from_E` on the hot path.
inline void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream) {
  transpiration_ = operating_point_flux(psi_stem, psi_upstream);
  Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(transpiration_) : leaf_temp_;
  if (use_energy_balance_) {
    update_temperature_dependent_params(Tleaf_);
    set_leaf_vpd(Tleaf_);
  }
  residual_flux_ = residual_flux(psi_upstream, transpiration_);
  // The corner the solve keeps (#183): below the wet bound the stomatal flux
  // reads as zero, counted, so a prescribed collar that inverts the pair by
  // rounding is the wet bound rather than a refusal.
  double E_s = transpiration_ - residual_flux_;
  if (E_s < 0.0) {
    clamps.note(CLAMP_STOMATAL_FLUX);
    E_s = 0.0;
  }
  stom_cond_CO2_ = atm_kpa_ * E_s * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_ +
                   cuticular_cond_CO2();
  ci_ = solve_ci(stom_cond_CO2_);
  // With no conductance at all (g_min = 0 and the stomata shut) solve_ci's
  // residual IS A(ci), so A at its root is zero by definition, and the root's
  // own rounding residue (~1e-15, either sign) is not a value: the product and
  // ratio curves' emergent lambda scales with or divides by A, and turned that
  // residue into platform-dependent +-1e-9, +-inf and NaN. On the tape the
  // theorem gives this A its rows through the conductance's g_min row.
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


// The same quantity for JS22's quadratic-in-the-drop cost:
//
//   dC/dpsi = 2*gamma*(psi - psi_up),   dE/dpsi = kmax*f(psi)
//
// ⚠️ `f` does NOT cancel here, where it does for TF24 (whose |f'| carries a factor
// of `f`). So this divides by the conductivity fraction. That is safe inside the
// feasible bracket -- `f` runs 1 down to `k_crit_fraction` = 0.05 at psi_crit, and
// never reaches zero -- but it is the reason this cannot be written in TF24's form.
inline double Leaf::lambda_JS22(double psi_stem, double psi_upstream) const {
  const double f = proportion_of_conductivity(psi_stem);
  return 2.0 * JS22_gamma * (psi_stem - psi_upstream) /
         (leaf_specific_conductance_max_ * f);
}


// And for CMax, whose marginal cost does not vanish at the wet end unless CMax_b is
// zero -- which is what makes it the one curve here that CAN be boundary-soil with a
// non-zero cost at the bound. `psi_upstream` enters only through dE/dpsi.
inline double Leaf::lambda_CMax(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double f = proportion_of_conductivity(psi_stem);
  return (CMax_a * psi_stem + CMax_b) /
         (leaf_specific_conductance_max_ * f);
}


// And for TF24_floor, which is the canonical decomposition of any cost into a part
// depending on the potential alone and a linear price of water:
//
//     Theta(E) = Theta~(psi) + lambda_o*E   =>   lambda = Theta~'(psi)/K(psi) + lambda_o
//
// with `K = kmax*f` the stem conductance. `Theta~` here is TF24's OWN cost, read at
// TF24's own traits, so the first term is `lambda_TF24` exactly and this is that
// function plus a constant. The `f` that does not cancel for JS22 cancels here,
// because it cancels for TF24 -- whose `|f'|` carries a factor of `f`.
//
// ⚠️ CALLED, NOT REPEATED, and the bit-exactness argument still holds. At
// `TF24_floor_lambda_o == 0` this must reproduce `lambda_TF24` bit for bit; it does,
// because the addend is an exact zero and `lambda_TF24`'s body ends in a DIVISION,
// which no fused multiply-add can absorb. (The JS22-based version of this curve
// had to repeat its parent's body instead, because that one ended in a multiply.)
inline double Leaf::lambda_TF24_floor(double psi_stem) const {
  return lambda_TF24(psi_stem) + TF24_floor_lambda_o;
}


// And for SOX, where the objective is a PRODUCT. At its optimum `A'g + Ag' = 0`, so
// `A' = -A(g'/g)` and
//
//     lambda = A'/E' = -A*(g'/g) / (kmax*f)
//
// ⚠️ IT CARRIES A FACTOR OF `A`, and that is what a product objective means rather
// than an artefact of this `g`: a difference objective prices water in absolute
// carbon, a product prices it as a share of current assimilation. `the-models.Rmd`
// derives it in one line and gives the scale-invariance test that pins it.
//
// ⚠️ It DIVERGES at `psi_crit`, where `g` is exactly zero. That is the mechanism by
// which a product objective can never be pinned at psi_crit -- water becomes infinitely
// expensive before the bound is reached -- and it means this is not a number to
// evaluate at the dry end of the bracket.
inline double Leaf::lambda_SOX(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double g = sox_reduction(psi_stem);
  const double gp = sox_reduction_deriv(psi_stem);
  const double f = proportion_of_conductivity(psi_stem);
  return -assim_colimited_ * (gp / g) / (leaf_specific_conductance_max_ * f);
}


// The same expression for Jones, where `-g'/g` collapses: with `g = 1 - psi/psi_crit`
// it is `1/(psi_crit - psi)`, so
//
//     lambda = A / [(psi_crit - psi) * kmax * f(psi)]
//
// which is the form `the-models.Rmd` quotes. It carries a factor of `A` and diverges
// at `psi_crit`, both for the product-objective reasons `lambda_SOX` records.
inline double Leaf::lambda_JW26(double psi_stem, double psi_upstream) const {
  (void)psi_upstream;
  const double f = proportion_of_conductivity(psi_stem);
  return assim_colimited_ /
         ((psi_crit - psi_stem) * leaf_specific_conductance_max_ * f);
}


// The same quantity for the normalised ProfitMax cost. Writing `f` for the
// conductivity fraction and `f'` for its slope, the cost is
// `HC + TC = (k_soil - kmax*f)/k_span + TC(Tleaf)`, so -- dividing through and
// restoring carbon units with |A|max:
//
//     dC/dpsi = kmax*|f'|/k_span + (dTC/dT)(dT/dE)(dE/dpsi),  dE/dpsi = kmax*f
//     lambda  = |A|max * [ |f'|/(f*k_span) + (dTC/dT)(dT/dE) ]
//
// The thermal term vanishes with either gate off -- with no energy balance
// dT/dE is zero, and with no thermal cost dTC/dT is.
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
  // S is the soil->collar conductance, and dE_from_soil_dpsi_collar now returns
  // exactly that: with both sides magnitudes it differentiates uptake with respect
  // to how hard the collar pulls, so it is positive by construction (#25). The
  // negation this line used to carry -- and the comment explaining why a
  // conductance came back negative -- are gone. That negation was the patch for
  // the bug that motivated #8: getting it backwards made lambda_multi negative,
  // which is how it was caught. It is now unrepresentable rather than guarded.
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

// The kernels' arguments, read off a caller's parameter pack.
//
// At double the pack holds exactly what the members hold (passive_pars builds it
// from them) and the temperature responses are taken at the temperature the
// block itself holds (photo_temp_), so filling a kernel from the pack and
// filling it from the members give the same numbers bit for bit -- which
// test_leaf asserts on both temperature paths rather than leaving to the reader.

// ⚠️ THE THREE RESPONSES BELOW DO NOT APPLY THE THERMAL-DAMAGE LAYER, so each
// refuses a leaf with it on rather than disagree with the member it twins.
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
  // The thermal cost reads leaf temperature and two thresholds, none of which is
  // a trait, so it is a factor here rather than a term carrying rows.
  if (use_thermal_cost_) {
    // Written out rather than `*=`, whose overload is ambiguous at a nested
    // scalar.
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
  return pars[par_kmax] * (stem_integral_at<T>(sigma, pars) -
                           stem_integral_at<T>(collar, pars));
}

template <class S, class Kernel, std::same_as<S>... Args>
inline S Leaf::kernel_slope_at(Kernel&& kernel, const char* what, const S& x,
                               const Args&... args) const {
  using odelia::util::to_passive;
  // Both levels are tangents and both sit at double, so none of this reaches a
  // tape whatever S is. The inner direction is x's, which makes the value a
  // slope; the outer walk carries x first and then each argument, which is where
  // the rows come from.
  using inner = odelia::ode::tangent_scalar<double>;
  constexpr std::size_t n_dir = 1 + sizeof...(Args);
  using outer = typename xad::fwd<inner, n_dir>::active_type;

  // ⚠️ THE FORWARD SOLVE TAKES THE INNER LEVEL ONLY. At double there are no rows
  // to carry, and the outer walk would cost n_dir directions for numbers nothing
  // reads -- on the path that places every operating point.
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

    // Every argument gets a row, in the order the walk set them. A braced
    // list is sequenced left to right, so the counter follows the walk; a stack
    // array rather than a vector because the rows cross as a span, which borrows
    // storage instead of owning it.
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

  // ⚠️ THE RESIDUALS SEE A HELD COLLAR, AND THE COLLAR'S CHANNEL IS ADDED BACK AS
  // ONE SUPPLIED SLOPE. That is the split closed_form_rows.hpp states, and here it
  // is forced rather than preferred: sigma is placed by an INVERSE TABLE, so its
  // response to the collar is that table's own slope -- while the implicit
  // function theorem on T1 divides by G'(sigma). The two agree only if the two
  // stem reads are exact mutual inverses, and they are not. Letting the
  // residual see a live collar therefore puts TWO values of dsigma/dcollar into
  // one object: measured, the lift said 0.6913 for dci/dcollar where the slope
  // beside it said 0.6382, and dM/dcollar came out an eighth of its size.
  //
  // Held, the residuals carry the TRAIT rows (which the inverse table has none
  // of) and the slopes below carry the collar's. One spelling of each.
  const S held(to_passive(collar));
  const double gc_per_flux =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const Leaf& leaf = *this;

  // dT1/dsigma: the stem curve at the operating point, scaled by kmax, and taken
  // from THE TABLE because that is the derivative the model itself forms. kmax
  // is read off the PACK the residual below multiplies by, not off the member
  // holding the same number, so the supplied slope is the slope of the residual
  // as written.
  const double dT1_dsigma =
      to_passive(pars[par_kmax]) * stem_curve_integral_deriv(sigma_star);
  // The inputs are handed over rather than left for the outer sweep to find:
  // the theorem's row is written here and the residual is rewound, so a walk of
  // it never reaches the caller's tape. ⚠️ EVERY ACTIVE THING THE RESIDUAL READS
  // MUST BE IN THIS LIST -- one left out is a row that never arrives, which
  // reads as an exact zero in that column rather than as an error. `held` is
  // to_passive and carries none.
  //
  // At a dry pin at the end of the stem curve the stem is HELD there
  // (seat_stem_at_dry_pin), so sigma is that end as a function of the traits,
  // not the root of T1: T1's slope there is kmax e^-40, and dividing by it hands
  // over amplified rounding (measured: forward and reverse disagreed by 0.4%).
  if (stem_held_at_dry_pin(sigma_star)) {
    sigma_h = stem_curve_end_at<S>(pars);
  } else {
    sigma_h = odelia::implicit_value<S>(
        sigma_star, dT1_dsigma,
        [&](const S& sg) -> S {
          return pars[par_kmax] *
                     (leaf.template stem_integral_at<S>(sg, pars) -
                      leaf.template stem_integral_at<S>(held, pars)) -
                 draw.E_up.value;
        },
        pars, draw.E_up.value);
  }

  // The stem's flux at the placed sigma. gc is built from THIS, not from the
  // draw, because that is what the model's ci solve reads -- the two agree only
  // to the stem table's round-trip. The residual flux comes off it and the
  // cuticle's conductance goes on, as the evaluator forms them (#182); at
  // g_min = 0 both are exact zeros carrying the g_min row.
  const S stem_flux = transpiration_at<S>(sigma_h, held, pars);
  const S r_star = residual_flux_at<S>(held, pars, draw);
  const S cuticle = pars[par_g_min] * S(g_min_co2_ratio_ * g_min_factor(Tleaf_));
  // ⚠️ WHERE NOTHING CONDUCTS, THE STOMATAL TERM CARRIES NO ROW. That is a shut
  // leaf at g_min = 0: it sits on its wet bound, where the stomatal flux is
  // zero identically -- the collar moves to keep it so -- and the collar's
  // channel below is skipped there, so nothing would cancel the held collar's
  // partial that stem_flux carries. Left in, it reached ci and A: the root_P50
  // row of a three-layer shut leaf's profit read 9.5e-2 where a difference of
  // the solve gives 1.3e-2. The cuticle keeps its g_min row, which is the
  // derivative of a leaf that starts to conduct.
  gc = stom_cond_CO2_ == 0.0
           ? cuticle
           : S(gc_per_flux) * (stem_flux - r_star) + cuticle;

  // ⚠️ FOLLOW THE BRANCH THE FORWARD MODEL TOOK. With nothing conducting at all
  // (g_min = 0 and the stomata shut) the residual is A(ci) = 0, which the
  // theorem still inverts where it has a root; where it has none the forward
  // reported ci = +inf and the saturated assimilation, and there is no
  // coordinate to place.
  //
  // A' comes from a TANGENT THROUGH THE SAME KERNEL the residual calls, which is
  // upstream's own way of taking it -- a hand-written slope twin would be a
  // second definition of one function, free to disagree with it while staying
  // finite.
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

  // Everything below is the collar's channel, and where the collar is held it is
  // multiplied by a step that is zero in value AND carries no derivative -- so
  // all of it records and none of it reaches a row. probe_tape_regions measures
  // 592 statements at a live collar against 591 at a held one: the channel is
  // worth one of them and the rest were recorded either way.
  //
  // ⚠️ THE SLOPES COME BACK NOT-A-NUMBER, WHICH IS DELIBERATE. They are absent
  // here, and an absent derivative spelled 0.0 is the one failure this codebase
  // cannot see -- a channel that has gone missing with every number still
  // finite. A caller that reads one gets NaN instead, which supply_derivatives
  // refuses by name. Only dprofit_dcollar_at reads a slope, and it asks for the channel.
  if (!collar_moves) {
    const S absent(std::numeric_limits<double>::quiet_NaN());
    return CollarCoords<S>{odelia::with_slope<S>{sigma_h, absent},
                           odelia::with_slope<S>{ci_h, absent}};
  }

  // The collar's own channel.
  //
  // ⚠️ THE CONDUCTIVITY HERE IS THE TABLE'S, WITH THE CURVE'S ROWS. Upstream's
  // dprofit_at_collar_psi forms every one of these from stem_curve_integral_deriv,
  // and the solve drove THAT assembly to zero -- so a closed-form value here
  // makes M non-zero at the collar the solve placed, and the envelope omission
  // stops being true. closed_form_rows.hpp is the same rule: the value is the table's
  // because the solve ran on the table, the rows are the curve's because a table
  // carries none.
  const S& kmax = pars[par_kmax];
  const S f_p = curve_with_closed_form_rows<S>(stem_curve_integral_deriv(to_passive(collar)),
                               held, pars[par_stem_P50], pars[par_stem_c]);
  const S f_sigma = curve_with_closed_form_rows<S>(stem_curve_integral_deriv(sigma_star),
                                   sigma_h, pars[par_stem_P50],
                                   pars[par_stem_c]);

  // dsigma/dcollar: THE INVERSE'S OWN SLOPE times the transport's, which
  // is how upstream's dprofit_at_collar_psi forms it and what actually placed
  // sigma. Read at E_up/kmax + G(collar) -- THE SOIL DRAW, the point the inverse
  // table was queried at, not the stem flux.
  const double u_at = to_passive(draw.E_up.value) / to_passive(kmax) +
                      to_passive(stem_integral_at<S>(held, pars));
  const S recip = S(1.0) / f_sigma;
  const S dsigma_du =
      S(stem_curve_integral_inverse_deriv(u_at)) + (recip - S(to_passive(recip)));
  // psi_crit does not move with the collar.
  //
  // ⚠️ NOR DOES A STEM HELD AT A DRY PIN. sigma_h is already the end of the curve
  // as a function of the traits, and the collar is placed to keep the stem
  // there, so its channel is an exact zero. Formed anyway it is e^40 (the
  // inverse's slope at the end) times terms of order e^-40 -- the cost's slope
  // and f_sigma -- which carry no information past rounding: measured, it moved
  // the profit rows by up to 8e-5 on macOS/arm64 and by up to 60% on Linux/x86,
  // where a difference of the solve agrees with the held stem's rows to 1e-9.
  const S V = stem_held_at_dry_pin(sigma_star)
                  ? S(0.0)
                  : dsigma_du * (draw.E_up.slope / kmax + f_p);

  const S gc_const = S(gc_per_flux);
  const S inv = S(inv_atm);
  // The conductance partials are the transport slope at each end, which is
  // kmax * f -- already formed above, so read rather than recomputed. The
  // residual flux moves with the collar through its ceiling, appended as the
  // solve's dprofit_at_collar_psi appends it; an exact zero at g_min = 0.
  const S dgc_dsigma = gc_const * kmax * f_sigma;
  S dgc_dp = -gc_const * kmax * f_p;
  dgc_dp = dgc_dp - gc_const * S(residual_flux_dpsi_upstream(to_passive(collar),
                                                             Tleaf_));

  // ⚠️ THE SAME CONDITION AGAIN, AND IT GOVERNS THE SLOPE AS WELL AS THE VALUE.
  // With no coordinate there is no residual to invert, and a leaf with no
  // conductance at all places ci by A(ci) = 0, which the collar does not enter.
  S dci_dp = S(0.0);
  if (!no_coordinate && stom_cond_CO2_ != 0.0) {
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
  // A leaf with no conductance and no ci at which it nets zero reported
  // ci = +inf and its saturated assimilation; the same here. One with no
  // conductance and a root reported A = 0 exactly (the residual there IS A),
  // so its value here is that zero, carrying the kernel's rows: the theorem
  // places ci where A stays zero, so the rows are what move it off zero.
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
    // TF24's own cost plus a part LINEAR IN THE FLUX: Theta(E) = Theta~(psi) +
    // lambda_o*E. Written as the parent's expression plus a term, which is what
    // makes the reduction at lambda_o = 0 exact rather than approximate -- an
    // exact zero added to TF24's exact value, fused multiply-add included.
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
  // Placed coordinates exist only at an interior point and only at the collar
  // operating_collar placed them at, so either check failing means the leaf moved
  // between the two calls -- and the coordinates would answer for another point
  // with every number finite.
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

  // The draws first. Profit reads the same flux wherever its collar is the live
  // one, so taking them here is what stops the supply reaching the tape twice.
  //
  // ⚠️ THE DRAWS CARRY SUPPLIED ROWS, NOT A SECOND SUPPLY. The draw already holds every
  // layer's uptake at the passive collar carrying the state's own rows; what
  // is missing is the channel through the collar MOVING, and that is one
  // supplied slope per layer. Taking the supply again at the live collar
  // records the whole of it a second time for the same numbers -- the step is
  // exactly zero in value, so every uptake here is the draw's, bit for bit.
  const S step = collar - S(draw.at);
  odelia::util::check_length(draw.duptake_dp.size(), draw.uptake.size());
  // Built rather than copied and then overwritten: copying the draw's vector
  // records one statement per layer for numbers the next line replaces.
  out.uptake.reserve(draw.uptake.size());
  for (std::size_t j = 0; j < draw.uptake.size(); ++j) {
    out.uptake.push_back(draw.uptake[j] + S(draw.duptake_dp[j]) * step);
  }
  (void)n;

  // ⚠️ PROFIT READS THE HELD COLLAR AT AN INTERIOR POINT, and that is the
  // envelope theorem written as an omission rather than as a term that has to
  // come out to zero. The marginal is zero there by the condition the solve
  // drove to nothing, so carrying the collar's movement through profit would add
  // M * dp/dtheta -- a term whose size is the solve's tolerance rather than the
  // model's. It is also why a collar that cannot be recorded costs the water
  // rows and not the objective. A shut leaf is a pinned one: its collar is
  // placed by the wet bound and profit carries its movement.
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
  // ⚠️ EVERY FACTOR IS A VALUE, so this is arithmetic at S and dM/dp is a FIRST
  // derivative of it. The lift this replaces took its value and first derivative
  // from a TABULATION and its second from the TRUE function, and the two disagree
  // by the interpolation error -- which is harmless at first order and is what
  // flipped the sign of a curvature where the collar came within 5.6e-08 of a
  // soil layer's potential.
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
  // dC/dpsi_stem, the same way: a tangent through the cost kernel rather than a
  // slope written out beside it.
  S C_prime = cost_slope_at<S>(at.sigma.value, pars);
  if constexpr (K == CostCurve::TF24_floor) {
    // d/dsigma of lambda_o * kmax * (G(sigma) - G(collar)) is lambda_o*kmax*f,
    // and the conductivity is the TABLE's for the reason every other one here
    // is: the solve ran on the table.
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

  // The bound moves water, so its slope is the supply's own conductance
  // there. The draw carries it as ONE number: where a kink makes the
  // analytic branch invalid -- a single-layer wet bound is exactly the
  // gravity-balance point -- supply_draw_at has already differenced it, so there
  // is no sentinel to test for here and no second definition of the fallback.
  const double dflux_dx = to_passive(draw.E_up.slope);

  // Nothing conducts: the bound is the gravity balance the draw carries,
  // with its soil row (find_shut_wet_bound placed the solve there).
  if (draw.no_conduction) {
    return draw.gravity_balance;
  }
  // The residual is the draw less the residual flux: uptake meets the
  // cuticle's loss at this collar, and vanishes there at g_min = 0.
  // implicit_value evaluates at the passive bound, and the draw was taken
  // there, so re-recording the supply would put the same expression on the
  // tape twice. The residual flux carries its rows (g_min, kmax, the stem
  // pair through the ceiling) and its own slope in the collar.
  //
  // ⚠️ AND AT g_min = 0 THE DRAW WAS TAKEN A STEP INSIDE THE ROOT, NOT AT IT.
  // A pinned solve returns root + step there, so the flux is flux(root) +
  // slope * step, and its parameter dependence carries slope's times the
  // step: an error in the row proportional to the step, visible wherever the
  // root's own derivative is zero (a one-layer wet bound is the gravity
  // balance, independent of the root curve). Taking the step's worth of slope
  // back off is the Taylor correction to the root, and the slope is already
  // on the tape. With a residual conductance the bound is a regular point
  // and the solve returns it exactly, so the step is zero.
  const double step = std::isfinite(wet_bound_root_)
                          ? bound_x - wet_bound_root_
                          : 0.0;
  const double dres_dx =
      dflux_dx - residual_flux_dpsi_upstream(bound_x, Tleaf_);
  const S held(bound_x);
  // ⚠️ EVERY ACTIVE INPUT THE RESIDUAL READS GOES IN THE LIST: the rewinding
  // form differentiates only those, and the roots' ceiling on the residual
  // flux is one when it is in force. Left out, the reverse sweep drops its
  // rows while a forward tangent keeps them; test_transpose sees exactly that.
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

  // The dry end is T1 with the stem held at its dry end (dry_end_at): the collar
  // at which the soil delivers exactly what the column can still carry. The stem
  // is held there, so only the second integral and the flux move with the collar.
  // Where the dry end is the end of the stem curve, the stem integral there is
  // G(inf) and its rows are G(inf)'s: the leaf transpires all the stem can carry.
  const double dT_dx =
      -to_passive(pars[par_kmax]) * stem_curve_integral_deriv(bound_x) -
      dflux_dx;
  const S dry = odelia::implicit_value<S>(
      bound_x, dT_dx,
      [&](const S& x) -> S {
        return pars[par_kmax] *
                   (leaf.template stem_integral_at<S>(dry_end_at<K, S>(pars),
                                                      pars) -
                    leaf.template stem_integral_at<S>(x, pars)) -
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
  // At every kind, not only where the kernels below are read: a bound is placed
  // by the supply alone, but billed_outputs_at refuses beside it, and a surface that
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
      // Closed by the RESIDUAL, like every other kind here. dM/dp is what the
      // theorem divides by; every parameter row is taped from M itself rather
      // than handed over as numbers from a second-order pass over a different
      // assembly.
      //
      // ⚠️ NOT ON THE DOUBLE PATH. At double every term of implicit_value
      // vanishes and this IS opt_root_psi_ -- but C++ evaluates the arguments
      // first, and the curvature costs two model evaluations. Charging the
      // forward model for a derivative it discards is what made a century stand
      // 191 s against 32 s.
      if constexpr (!std::is_same_v<S, double>) {
        const double curvature =
            std::isnan(interior_curvature)
                ? this->template profit_curvature_at_collar<K>()
                : interior_curvature;
        // ⚠️ odelia's implicit_value refuses only a zero or non-finite
        // denominator, so these two are this arm's to refuse: a curvature of
        // -1e-300 otherwise returns rows of ~1e296 with nothing flagged. Two
        // messages because they are two findings: an interior collar sits on a
        // downward crossing of the marginal, so a curvature that is not negative
        // indicts the curvature, while a tiny negative one is a flat optimum.
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
        // ⚠️ THE COORDINATES ARE PLACED BEFORE THE RESIDUAL, NOT INSIDE IT. The
        // residual is recorded, swept and rewound, so anything it solves is gone
        // when this returns; placed out here they stay on the caller's tape for
        // billed_outputs_at to read, and the residual adds only the collar's channel.
        // They are then inputs the residual reads, so they join the list: one
        // left out is a column of exact zeros in the collar's row.
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
      // The continuity root: T1 with the stem held at its critical potential.
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
  // Absolute, not relative: the collar runs 0.5 to 6 MPa, so one step serves the
  // whole range, and the error-versus-step curve has a floor three decades wide
  // around it. Read off dprofit_at_collar_psi's own values: 1e-04, 1e-05 and
  // 1e-06 agree to seven figures, 1e-02 is truncation and 1e-07 is round-off.
  const double h = 1e-5;
  bool up_ok = false, down_ok = false, here_ok = false;
  const double up = dprofit_at_collar_psi<K>(p0 + h, &up_ok);
  const double down = dprofit_at_collar_psi<K>(p0 - h, &down_ok);
  // ⚠️ ONE SIDE IS ENOUGH, AND AN INTERIOR COLLAR CAN BE A STEP FROM A BOUND.
  // The feasible interval is closed and its ends are reachable: a collar sitting
  // within h of one has a neighbour outside it, where dprofit refuses. Returning
  // NaN there hands the caller nothing to divide by, and plant's interior
  // derivation refuses the whole SWEEP on it -- so a point perfectly well placed
  // in the interior costs a gradient because of where its neighbour fell. A
  // one-sided difference at the same step is a derivative of the same function,
  // to one order lower, which is the right answer rather than no answer.
  const double here =
      (up_ok && down_ok) ? 0.0 : dprofit_at_collar_psi<K>(p0, &here_ok);
  // Back to the point the solve placed, because this drove the model away from
  // it.
  //
  // ⚠️ AND THE KIND HAS TO BE PUT BACK BY HAND. evaluate_root_collar_psi restores
  // every NUMBER -- collar, stem potential, ci, profit, every layer's draw, bit
  // for bit -- and then tags the point `prescribed`, because from its side that
  // is what it is. operating_collar switches on the kind, so leaving it there turns a
  // solved interior point into one with no condition to place it at: the failure
  // is a throw, not a wrong number, and it surfaces in a caller that never asked
  // for a curvature.
  const OperatingPointKind kind0 = operating_point_kind_;
  evaluate_root_collar_psi_for<K>(p0);
  operating_point_kind_ = kind0;
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

// Sicangco et al. (2026) Eqns 8-9: the fraction of maximum PSII damage sustained
// at this leaf temperature, in [0,1]. `r` follows gsthermal, not the printed Eqn
// 8 -- see the note on use_thermal_cost_.
//
// Exactly 0.0 with the gate off, so every caller can add it unconditionally.
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


// `CF77_lambda_ * E`. Written into `hydraulic_cost_` like the TF24 cost and unlike the
// Sperry one, because this product IS a carbon flux: `CF77_lambda_` is umol C per kg of
// water and E is kg m^-2 s^-1.
//
// `operating_point_flux()` reads the memoised `transpiration()`, so calling it here
// straight after set_leaf_states_rates_from_psi_stem costs nothing and returns
// exactly the E of the operating point that call established -- including the
// no-flow case, which both read as zero by the same rule.
inline double Leaf::hydraulic_cost_CF77(double psi_stem,
                                                double psi_upstream) {
  hydraulic_cost_ = cf77_price() * operating_point_flux(psi_stem, psi_upstream);
  return hydraulic_cost_;
}


// The CF77 price, with the optional Medlyn soil-moisture factor.
//
//   beta = (theta - theta_w) / (theta_fc - theta_w),   price = lambda / beta
//
// so the price is `CF77_lambda_` in soil at or above field capacity and rises
// without bound as theta approaches the wilting point. That direction is the whole
// point: it gives Cowan-Farquhar the soil-moisture shutdown it structurally lacks.
//
// ⚠️ BETA IS CLAMPED AT 1, so theta above field capacity does not price water
// BELOW `CF77_lambda_`. Two reasons. It matches the convention the beta-factor
// literature uses, and it makes `CF77_lambda_` a floor rather than a value the
// wet end can undercut -- which is what the parameter is for. The consequence
// worth knowing when reading output: for any theta >= theta_fc this curve is
// plain CF77, exactly, not approximately.
//
// The off path returns the member untouched rather than dividing by a computed
// 1.0, so every existing CF77 number is bit-identical. Do not "simplify" that.
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


// gamma * (dpsi)^2, in carbon units like the TF24 and CF77 costs -- so it goes in
// `hydraulic_cost_` on the same footing as those two and not in ProfitMax's
// separate normalised members.
//
// `dpsi*dpsi` rather than `pow(dpsi, 2.0)`: this runs inside the optimiser's inner
// loop, and hazard 5 is about exactly that.
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


// The integral of `CMax_a*psi + CMax_b` from psi_upstream to psi_stem, so no flow
// costs nothing:
//
//     C = (psi - psi_up) * [ CMax_a*(psi + psi_up)/2 + CMax_b ]
//
// ⚠️ FACTORED, NOT `a*(psi^2 - psi_up^2)/2 + b*(psi - psi_up)`, and that is a
// measured requirement rather than tidiness. The difference-of-squares form
// subtracts two nearly equal numbers as the potentials converge: against the
// factored value it is wrong by 1.6e-08 at a 1e-08 drop and 1.2e-04 at 1e-12,
// the latter being what this guide calls a real difference rather than rounding.
// A gradient perturbs by a relative 1e-06, so an operating point near no flow
// reaches that regime.
//
// ⚠️ THIS IS ALSO WHY IT DOES NOT SHARE A KERNEL WITH JS22, tempting as that is.
// JS22 is this family's `(a, b) = (2*gamma, -2*gamma*psi_up)` member, and
// putting it through the form above computes `gamma*(psi + psi_up) - 2*gamma*psi_up`
// -- reintroducing exactly the cancellation that `gamma*dpsi*dpsi` does not have.
// The two share a family, not an implementation.
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


// TF24's cost plus CF77's price. Both are carbon fluxes, so the sum goes into
// `hydraulic_cost_` on the same footing as each of them alone.
//
// ⚠️ THE TWO REDUCTIONS ARE BIT-EXACT BY CONSTRUCTION, WHICH IS WHY EACH TERM IS
// THE PARENT'S OWN EXPRESSION. `hydraulic_cost_TF_kernel(psi_stem)` is the
// function `hydraulic_cost_TF` assigns from, and `TF24_floor_lambda_o *
// transpiration(...)` is `hydraulic_cost_CF77`'s product in that operand order.
// So `lambda_o = 0` adds an exact zero to TF24's exact value, and
// `TF24_cost_scale = 0` -- which zeroes the kernel, since `pow(1-f, beta2)` is
// finite inside the bracket -- adds an exact zero to CF77's.
//
// ⚠️ AND IT SURVIVES FUSED MULTIPLY-ADD, which is the part that is not obvious and
// the reason no `rounded()` barrier is needed here (contrast gradient.hpp, where
// one is REQUIRED). Contraction has two candidates -- fuse the kernel's own
// `scale * pow(...)` with the final add, or fuse `lambda_o * E` with it -- and the
// zeroed parameter makes the fused addend an exact zero either way:
// `fma(x, y, 0)` is `round(x*y)`, which is what the unfused form computes, and
// `fma(0, y, q)` is `q`. So both reductions hold under `-ffp-contract=fast` as
// well as under `off`. Reassociating would break that, so `test_tf24_floor_reduces`
// is what to re-run after touching this line.
//
// `operating_point_flux()` reads the memoised `transpiration()`, so calling it
// here straight after set_leaf_states_rates_from_psi_stem costs nothing -- see
// hydraulic_cost_CF77, which relies on the same thing.
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
// Eller's reduction factor: the conductivity fraction rescaled so it is 1 at zero
// tension and exactly 0 at `psi_crit`, where `f` reaches `k_crit_fraction` by
// construction. That last equality is why this needs no parameter of its own -- the
// 0.05 is the same one `psi_crit` is derived from and the same one Sperry's
// `k_crit` reads, so the three cannot disagree.
// ⚠️ CLAMPED AT ZERO, and that is not defensive padding -- without it this curve
// reports a shut-down leaf as PROFITABLE. Past `psi_crit` the conductivity fraction
// falls below `k_crit_fraction`, so `g` goes negative; the objective is `A*g` and a
// shut-down leaf has `A = -R_d < 0`, so the product of two negatives came back
// POSITIVE and grew as the soil dried. Measured before the clamp, at PPFD 1500:
// psi_soil 6.0 gave +0.0125 and 7.0 gave +0.0633, against TF24's -8.48 and -8.85.
// Sabot's `phiLWP` clamps for the same reason; her `kcost` does not, because
// `TractLSM` never evaluates outside `[Ps, Pcrit]` and the no-flow branch here does.
//
// Written as `g < 0` rather than `std::max(0.0, g)` so a NaN propagates instead of
// being silently turned into a zero: `NaN < 0.0` is false, so the NaN is returned.
inline double Leaf::sox_reduction(double psi_stem) const {
  const double g = (proportion_of_conductivity(psi_stem) - k_crit_fraction) /
                   (1.0 - k_crit_fraction);
  return g < 0.0 ? 0.0 : g;
}


// Its slope. `|f'| = f*(c/b)*(psi/b)^(c-1)` is the same analytic form `lambda_TF24`
// uses, so this needs neither AD nor a table slope; `f` decreases, so this is
// negative.
inline double Leaf::sox_reduction_deriv(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  const double abs_fprime =
      f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
  return -abs_fprime / (1.0 - k_crit_fraction);
}


// ⚠️ A PRODUCT, so what comes back is NOT in carbon units, unlike every other
// `profit_psi_stem_*` here. The argmax is unaffected -- `A*g` and `log A + log g`
// share it because `log` is monotone, and the paper's own first-order condition is
// the log-differentiated form -- but the VALUE is a different kind of number and
// must not be compared with a TF24 profit.
//
// ⚠️ `hydraulic_cost_` IS EXPLICITLY POISONED, which is hazard 8 rather than
// tidiness. There is no subtracted cost on this curve, so leaving the field alone
// would report whichever curve ran last -- and on one reused `Leaf` that is exactly
// the stale-state leak the `psi_stem_optima` golden file's second pass exists to
// catch.
inline double Leaf::profit_psi_stem_SOX(double psi_stem,
                                        double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * sox_reduction(psi_stem);
}


// Jones' linear reduction factor, on the DERIVED `psi_crit`. Clamped at zero for
// exactly the reason `sox_reduction` is, and the NaN-propagating form for the same
// reason; here the clamp binds for every `psi > psi_crit` rather than only past the
// point where `f` drops under `k_crit_fraction`.
inline double Leaf::jw26_reduction(double psi_stem) const {
  const double g = 1.0 - psi_stem / psi_crit;
  return g < 0.0 ? 0.0 : g;
}


// Constant, which is the whole difference from SOX: a linear `g` prices each further
// MPa of tension identically, where SOX's follows the vulnerability curve and so
// charges little until the curve turns over. Zero once the clamp binds.
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

// `W = E_mol + beta*Vcmax_mol`, mol m^-2 s^-1. Both conversions in one place; see
// the declaration for why beta's dimensionlessness requires them.
//
// ⚠️ IT READS `vcmax_`, NOT `vcmax_25`, so with the energy balance on it moves with
// the solved leaf temperature -- which is what makes `A/W` the ratio Prentice et
// al. minimise rather than a ratio to a reference-temperature capacity. `vcmax_` is
// re-derived per candidate potential by `set_leaf_states_rates_from_psi_stem`, and
// this function is called from `cost_deriv` (which does not call it) as well as
// from the profit (which does). At a fixed leaf temperature the two agree because
// `vcmax_` does not move. With the energy balance ON, both derivatives read it at
// the leaf temperature they solved for, and `dprofit_energy_balance_term` adds
// the capacity term's own motion with that temperature, which `cost_deriv` leaves
// out.
inline double Leaf::leastcost_flux(double psi_stem, double psi_upstream) {
  return operating_point_flux(psi_stem, psi_upstream) * kg_to_mol_h2o +
         LeastCost_beta * vcmax_ * umol_to_mol;
}


// ⚠️ A RATIO, not a carbon flux, and `hydraulic_cost_` IS EXPLICITLY POISONED for
// the reason `profit_psi_stem_SOX` gives: there is no subtracted cost on this
// curve, so leaving the field alone reports whichever curve ran last on a reused
// `Leaf`.
inline double Leaf::profit_psi_stem_LeastCost(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ / leastcost_flux(psi_stem, psi_upstream);
}


// `lambda = A*kg_to_mol_h2o/W`. Reads the solved `assim_colimited_` rather than
// re-deriving it, exactly as `lambda_SOX` does, so it means what it says only at
// the potential the leaf was last evaluated at -- which is how the optimisers call
// it, immediately after writing `profit_`.
inline double Leaf::lambda_LeastCost(double psi_stem, double psi_upstream) {
  if (!use_energy_balance_) {
    return assim_colimited_ * kg_to_mol_h2o /
           leastcost_flux(psi_stem, psi_upstream);
  }
  // With the energy balance on, W's capacity term moves with E through Tleaf,
  // so dTheta/dE gains beta*Vcmax_mol'(T)*dT/dE -- the slope the first-order
  // condition takes in dprofit_energy_balance_term, as lambda_ProfitMax carries
  // its thermal term.
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
// inherited. See the declaration for what it cost when this did not exist.
inline void Leaf::clear_collar_solve_state() {
  operating_point_kind_ = OperatingPointKind::Unsolved;
  opt_root_psi_ = util::na_value;
  E_up_ = util::na_value;
  soil_consumption_.assign(soil_consumption_.size(), util::na_value);
}

// A single-layer optimiser maximises over the CLOSED interval
// [psi_soil, psi_crit]. The endpoints have to be candidates: the objective is
// maximised at full closure whenever water is priced above what the carbon is
// worth, and a bracketing search steps in from the bounds and so can never
// return one. Some cost curves also carry a second interior hump, which the same
// search cannot see past. See util::maximise_over_closed_interval_foc, which is
// the one solver both routes reach.
//
// The root-based solve reaches the same conclusion by a different route:
// maximise_profit_over_collar tests the gradient's sign at each end and reports a
// pinned optimum explicitly.

// ===========================================================================
// Sperry et al. (2017) ProfitMax, as Sicangco et al. (2026) implement it
// ---------------------------------------------------------------------------
// Sperry maximises `Profit = CG - HC` with both terms normalised. Multiplying
// that objective by |A|max turns it into `A - lambda*(k(psi_soil)-k(psi))`, the
// same function up to a positive scale factor and so the same argmax, with lambda
// taking the value below. That is a property worth knowing rather than a second
// entry point: the lambda is not a constant, so prescribing one is not this
// model.
//
//     CG      = A(psi)/|A|max
//     HC      = [k(psi_soil)-k(psi)] / [k(psi_soil)-kcrit]
//     lambda* = |A|max / [k(psi_soil) - kcrit]
//
// Checked numerically on a 4001-point grid at gross assimilation, at net, and at
// a leaf temperature of 48 C where CG is negative throughout: same grid point
// every time. So the lambda form is not wrong -- it is unusable, because lambda*
// is not a constant. |A|max moves with every driver and k(psi_soil) with the
// soil, so a caller has to rescan the supply stream and recompute lambda at every
// observation, from R, at ~1.8 us of call overhead per crossing against a ~3 us
// model. This does it in C++ and reports the profit in Sperry's own units, so a
// figure from this package can be laid over one from the paper.
//
// THE NORMALISATION IS NOT COSMETIC, which is the other reason to have it. HC
// runs 0 to 1 across the operating range HOWEVER WIDE the vulnerability curve is,
// because the denominator rescales with the curve -- so a ProfitMax cost cannot
// become numerically negligible the way the TF24 cost does as stem_b widens
// (measured at 1.3e-5 at stem_b = 200). Both terms are scale-free, which also
// means the plant's willingness to spend water does not depend on how much carbon
// is at stake in absolute terms. That is exactly what TF24_cost_scale decides on
// the TF24 path, and it is the parameter a calibration cannot identify without
// leaf water potential.
//
// ⚠️ kcrit IS 5% OF kmax, and that is the same convention this package already
// has rather than a new one. Sperry, Sabot et al. (2020) and Sicangco all set
// kcrit = 0.05*kmax, i.e. psi_crit is P95; at this package's defaults
// f(psi_crit) = 0.0500 exactly, so nothing has to be reconciled.
// The potential Sperry's normalisation measures conductance loss FROM. On one soil
// potential that is psi_soil. On a layered profile there is no scalar psi_soil, and
// the right generalisation is the collar potential at which uptake is exactly zero
// (`root_zero_E`, which is prepare_collar_solve's wet bracket end): the potential at
// which the soil delivers no water is what the plant experiences as "the soil".
// Single-layer callers get psi_soil unchanged, so nothing there moves.
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
  // ⚠️ A NON-FINITE REFERENCE IS PROPAGATED, NOT REFUSED, and that is deliberate.
  // This is called on leaves whose drivers are not seated yet -- the batch route
  // seeds its |A|max capture before applying an observation -- and a hard error
  // there breaks a path that has always worked. NaN norms flow into the optimiser's
  // `!(profitmax_k_span_ > 0)` guard, which reports the degenerate shut state. That
  // guard is the one place this is decided; do not add a second.
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

  // |A|max BY OPTIMISATION, NOT BY A SCAN. Sperry's carbon gain is normalised by
  // the largest |A| on the supply stream, and this used to find it with a
  // 500-point sweep. Two reasons that was the wrong instrument:
  //
  //   * the argmax of a scan is piecewise constant in the parameters, so
  //     d|A|max/dtheta is a staircase -- which is the whole reason ProfitMax's
  //     gradient had a partial/total split and needed a coarse finite-difference
  //     step where every other route uses 1e-06;
  //   * it is 500 evaluations to locate something that measurement puts at the DRY
  //     BOUND on 1318 of 1320 driver rows, interior on 2, with the interior case
  //     beating the bound by at most 1.1e-03.
  //
  // So: endpoints plus a root-find on dA/dpsi == 0, the same method the objective
  // itself uses. ONE maximisation of A suffices, and that is a statement about the
  // model rather than a shortcut: on this interval the only way |A| could peak at a
  // MINIMUM of A is where A is negative, which happens at full closure (excluded
  // below) and along the whole stream where the leaf cannot fix carbon even at
  // ambient CO2. Elsewhere A is the carbon the leaf gains, so |A|max is max A.
  //
  // dA/dpsi needs no new algebra: for an identity link dJ/dpsi = A'*dci/dpsi - C',
  // so adding the cost derivative back recovers the benefit term exactly.
  //
  // ⚠️ THE CLOSURE POINT IS EXCLUDED, as it was from the scan. At the wet end
  // the stomata are shut: a closed stoma, not a point on the stream, and
  // gsthermal drops it too (`A[E == 0] <- NA`). The scan expressed that as "skip
  // grid point 0"; with no grid, the search starts one probe inside the bound.
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
  // ⚠️ THE ONE CASE WHERE max|A| IS NOT max A, and it is reachable. If A is
  // NEGATIVE across the whole stream -- a leaf hot enough that gross assimilation
  // cannot cover respiration anywhere, measured at 48 C -- then the largest
  // magnitude is the most negative value, and normalising by zero would leave
  // ProfitMax with no objective at all. Solved the same way, by maximising -A, and
  // only in that case: one extra root-find on the rows that need it and none on the
  // rows that do not.
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
  // At the operating-point leaf temperature, which on the energy-balance path is
  // the one set_leaf_states_rates_from_psi_stem just solved for. Off it the leaf
  // is at the prescribed temperature and TC is constant across the stream --
  // which is Sicangco's ProfitMaxTC run with the energy balance disabled, and is
  // why their thermal cost does nothing without it.
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

// ⚠️ THE ROUTE IS CHECKED, NOT JUST THE CURVE. Both of these optimise or
// evaluate psi_stem with the upstream potential pinned at psi_soil, so calling
// one on a leaf seated on the collar route would answer a question about a
// different model while `model_route()` said otherwise -- the same "half a model
// configured in each place" incoherence `set_model` exists to remove.
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

// ⚠️ RETURNS (value, feasible) LIKE THE COLLAR VERSION, and for the same reason:
// `dprofit` hands back a hard 0.0 SENTINEL on its shut-down and reversed-gradient
// exits, and a bare zero is indistinguishable from a stationary point. A composite
// that reads the value without the flag inherits that bug.
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
    // ⚠️ LeastCost REACHES THIS ARM BY A DIFFERENT ROUTE FROM THE OTHER TWO, and
    // the difference is worth keeping in view. SOX and JW26 are products `A*g`,
    // whose log is `log A + log g`; LeastCost is a RATIO `A/W`, whose log is
    // `log A - log W`. A product's log-cost is `-log g` and a ratio's is `+log W`,
    // so the sign lands in `cost_deriv` rather than here -- the link itself is the
    // same `h = log`, and `h' = 1/A` is all this table reports.
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

// `h'(A)`, the only thing the derivative below needs from the link.
//
// ⚠️ The Log arm is why a product objective's gradient exists at all: `A*g` and
// `log A + log g` share an argmax, and differentiating the second is the same
// additive shape as every other curve with `1/A` on the benefit term.
//
// ⚠️ It is also where a product curve's ONE genuine hazard lives. `log A` needs
// `A > 0`, and at full closure `A = -R_d < 0`. The optimiser sidesteps this by
// maximising the product directly; a DERIVATIVE cannot, so a non-positive `A`
// returns NaN rather than a number, and the caller's feasibility flag is what
// carries that outward.
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
// these three, and nothing else. Every arm is explicit and the last asserts, so a
// curve added to the enum and forgotten here is a compile error rather than a
// silent fall-through to whichever arm happened to be the `else`.

// What must be set before the search. TF24, SOX and JW26 have nothing: TF24's
// parameters are traits with defaults, and the two product curves take no
// parameter at all -- their `g` is built from the vulnerability curve and
// `k_crit_fraction`.
template <Leaf::CostCurve K>
inline void Leaf::check_cost_parameters() {
  if constexpr (K == CostCurve::CF77) {
    // CF77_lambda_ is an input with an aligned default, so a non-finite value
    // means a caller cleared it. That maximises a NaN objective, whose argmax is
    // a property of the bracket and looks entirely plausible. Refuse instead.
    if (!std::isfinite(CF77_lambda_)) {
      util::stop("the CF77 cost curve needs CF77_lambda_ set: it is the "
                 "PRESCRIBED marginal value of water in umol C (kg H2O)^-1, and "
                 "it was set to " + util::to_string(CF77_lambda_) + ". Its "
                 "default is aligned to TF24 (see align_curves()); "
                 "set_physiology and set_traits never touch it.");
    }
    // ⚠️ Only when the option is ON, so a caller who never touches it cannot be
    // refused for soil-moisture values the bare constant never reads. beta <= 0
    // means theta at or below the wilting point: the price is infinite or
    // negative, and a negative price PAYS the leaf to transpire, which the
    // optimiser will happily take to the wet bound. Refuse rather than return a
    // bracket property that looks like an operating point.
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
    // A trait with a default, so an unset one is a programming error rather than
    // a caller omission -- but still checked, because a NaN here maximises a NaN
    // objective and returns a bracket property that looks like an operating point.
    if (!std::isfinite(JS22_gamma) || JS22_gamma < 0.0) {
      util::stop("JS22_gamma must be a finite, non-negative hydraulic unit cost "
                 "in umol C m^-2 s^-1 MPa^-2; got " + util::to_string(JS22_gamma));
    }
  } else if constexpr (K == CostCurve::CMax) {
    // ⚠️ CMax_a only gets the sign check. CMax_b is SIGNED -- negative in the
    // paper's own convention, see the member -- so rejecting a negative would
    // reject the literature.
    if (!std::isfinite(CMax_a) || CMax_a < 0.0) {
      util::stop("CMax_a must be a finite, non-negative slope for the CMax "
                 "marginal cost in umol C m^-2 s^-1 MPa^-2; got " +
                 util::to_string(CMax_a));
    }
    if (!std::isfinite(CMax_b)) {
      util::stop("CMax_b must be finite; got " + util::to_string(CMax_b));
    }
  } else if constexpr (K == CostCurve::TF24_floor) {
    // ⚠️ ONE PARAMETER, AND IT IS AN INPUT WITH NO DEFAULT -- unlike
    // CF77_lambda_, which has an aligned one -- because zero is a perfectly good
    // value here, and it is the value at which this curve stops being itself. The hydraulic
    // half is TF24's cost at TF24's own traits, so at lambda_o = 0 TF24_floor IS
    // TF24, and a caller who seated this curve and never named a price would be
    // running the production model under the new model's name. Refuse the
    // omission; an explicit zero is still accepted and is the reduction test's own
    // case. The hydraulic half needs no check, for TF24's reason: its parameters
    // are traits with defaults.
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
    // An input with an aligned default, refused when cleared for CF77's reason
    // rather than TF24_floor's: there is no other curve this one collapses into,
    // so a non-finite beta is a NaN objective, whose argmax is a property of the
    // bracket and looks entirely plausible.
    //
    // ⚠️ AND STRICTLY POSITIVE, WHICH IS NOT THE SAME CHECK `TF24_floor_lambda_o`
    // GETS. Zero is accepted there and is the reduction test's own case; here it
    // deletes the capacity term, leaving `A/E` -- water-use efficiency, which rises
    // without bound as transpiration goes to zero, so the optimiser returns the wet
    // bound for every driver set. That is a degenerate objective rather than a
    // limiting case of this one.
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


// Which objective. ⚠️ These are NOT all the same kind of number: five are carbon
// (TF24, CF77, JS22, CMax, TF24_floor), two are products (see
// `profit_psi_stem_SOX`) and one is a ratio (see `profit_psi_stem_LeastCost`).
// That is already true of `profit_` across the optimisers and is not introduced
// here.
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


// Which lambda to report. ⚠️ TF24's reads psi_stem ALONE where the rest need both
// potentials -- JS22 because its cost is a function of the drop, the others
// because E is measured from the upstream one. Normalised to two arguments here so
// the body below never has to know which. CF77's is exact by construction rather
// than derived: its cost is lambda*E, so dC/dpsi over dE/dpsi is lambda
// identically.
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
    // The hydraulic half reads the absolute potential, so like TF24's this needs
    // only one of the two potentials.
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
// where the collar psi_r is the free variable and psi_stem follows it. Uptake
// gives dE/dpsi_r = S, and the stem path E = kmax [G(psi_stem) - G(psi_r)] gives
// dpsi_stem/dpsi_r = (S + kmax f(psi_r)) / (kmax f(psi_stem)). So, for the part
// of a cost that reads psi_stem alone, with its stem-route lambda C_s/(kmax f),
//
//     lambda = lambda_stem * (1 + kmax f(psi_r) / S),
//
// the series correction. The part of a cost that reads the potentials only
// through E (CF77's lambda*E, LeastCost's log W, the floor's lambda_o*E,
// ProfitMax's thermal term) moves with E directly, so its lambda carries over
// unchanged. A cost that reads the collar explicitly as well (JS22's drop,
// CMax's lower limit) adds its partial dC/dpsi_r over S.
//
// ⚠️ THIS USED TO BE marginal_cost_water_multilayer() FOR EVERY CURVE, which is
// the TF24 row of this table: every collar solve reported TF24's lambda whatever
// it had maximised. CF77 at a price of 1e5 reported 1.05e6 where its emergent
// lambda is, by construction, its price. TF24's row is still that function, so
// its number is bit-identical. The check that every row is right is that at an
// interior optimum lambda_emergent equals dA/dE along the collar path, for every
// curve (test_collar_lambda_emergent_is_dA_dE).
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


// And the body, written once. Maximises over the CLOSED interval for the reason
// hazard 11 gives: the objective is highest at full closure whenever water is
// priced above what the carbon is worth, and a bracketing search that steps in
// from the bounds can never return one.
//
// ============================================================================
// WHAT EACH CURVE DOES DIFFERENTLY IN HERE, and it is less than you would guess
// ----------------------------------------------------------------------------
// This used to be seven one-line functions -- `optimise_psi_stem_TF()` and six
// siblings -- each `{ optimise_psi_stem_single<K>(); }` and each carrying one of
// the notes below. The functions are gone; the notes are the part worth keeping,
// so they are collected here rather than attached to an alias apiece.
//
// TF24 -- the reference. Wet-pinning is its dominant pinned class: 24 rows at
//   25 C and 80 at 40 C over the golden grid.
//
// CF77 -- same shape, objective the only difference. Shares TF24's degenerate
//   convention: at a soil potential drier than psi_crit the objective is
//   EVALUATED at the no-flow point rather than zeroed, so every reported field
//   describes that point. Profit there is `-R_d`, since E is exactly zero and so
//   is the cost.
//
// JS22 -- ⚠️ NO SECOND HUMP, and the scan is kept anyway. Its marginal cost rises
//   monotonically from zero, so against a saturating marginal benefit the
//   first-order condition crosses exactly once; hazard 11's double crossing needs
//   a marginal cost that FALLS at the dry end, which |f'| does and a quadratic
//   does not. The endpoints still have to be candidates.
//   ⚠️ THE WET END IS INTERIOR WHENEVER THERE IS A BRACKET AT ALL. dC/dpsi -> 0 as
//   psi -> psi_soil while dA/dpsi > 0 there, so dprofit/dpsi is strictly positive
//   at the wet bound. Measured at gamma = 20, PPFD 200, VPD 4 -- a regime built to
//   make closure attractive -- psi* still clears psi_soil by 1.8e-02 at psi_soil
//   5.5 and 1.1e-02 at 5.84.
//   ⚠️ Stated that way rather than as "can never be boundary-soil", which the same
//   measurement falsifies: at psi_soil 6.0, above the default psi_crit of 5.870,
//   psi* comes back EXACTLY psi_soil. That is the no-flow branch below -- there is
//   no feasible interval to be interior in -- not a pinned optimum. A fit on this
//   curve should therefore reach the exact-gradient route more often than TF24.
//
// CMax -- ⚠️ UNLIKE JS22, CAN BE WET-PINNED, and `CMax_b` is what decides it.
//   JS22's marginal cost vanishes as the drop closes; this one approaches
//   `CMax_a*psi_soil + CMax_b`, which is non-zero. So the closed interval is doing
//   real work here rather than guarding a case that cannot arise.
//   ⚠️ AND A WET-PINNED ANSWER IS NOT `psi_soil` EXACTLY. Measured at psi_soil 3,
//   PPFD 200, VPD 4: `psi* - psi_soil` is 5.930e-06 at CMax_b = 6 and **the same
//   5.930e-06** at 8, 12 and 20. A value that does not move with the parameter is
//   the tell -- it is the bracket's step-in fraction (~1e-06 of a 2.870 MPa width),
//   so the answer is determined by the bound rather than by the objective. A
//   caller testing `psi* == psi_soil` to detect closure will conclude this curve
//   never pins. Compare against the bound with the step-in tolerance, or read the
//   gradient's sign at the bound.
//
// SOX -- the first PRODUCT objective, and correct through the same maximisation
//   because `A*g` and `log A + log g` share an argmax.
//   ⚠️ NEVER DRY-PINNED, by construction rather than by luck: `g(psi_crit)` is
//   exactly zero, so the objective is zero at the dry bound while any interior
//   point with positive `A` beats it. The wet bound IS reachable -- at full
//   closure `A = -R_d` and the product is `-R_d * g(psi_soil)`, which can be the
//   maximum when assimilation is negative throughout (hazard 11).
//   ⚠️ NO PARAMETER TO VALIDATE, which is the appeal: `g` is built from
//   `proportion_of_conductivity` and `k_crit_fraction`. Nothing here can be unset.
//   ⚠️ THE BENEFIT LINK IS WHAT MAKES THIS A ROW OF THE COST TABLE RATHER THAN AN
//   EXCEPTION TO IT, exactly rather than approximately. Write `C(psi) = -log g`
//   and take `h = log`; then
//
//       h(A) - C          = log A + log g = log(A*g)
//       d/dpsi[h(A) - C]  = A'/A + g'/g   = (A'*g + A*g') / (A*g) = P' / (A*g)
//
//   so the two first-order conditions have IDENTICAL roots and nothing is dropped
//   by the link sitting on the benefit alone. It can sit there because `C` is a
//   function of psi ONLY, so `-log g` is itself just another cost curve;
//   transforming the cost buys no generality. `A` is not free that way -- it comes
//   out of the `ci` root-find -- so its transform has to be declared, and `h'` is
//   the whole of what the derivative needs from it. Verified numerically over six
//   driver rows on both product curves: the log-form derivative at a 400k-point
//   scan argmax of the PRODUCT is 1.6e-07 to 4.0e-06, which is that scan's own
//   resolution times the curvature, and the solver lands within 3.2e-06 MPa of it.
//   ⚠️ THE OPTIMISER STILL MAXIMISES THE PRODUCT, deliberately. `log` is monotone
//   so the argmax is identical, and the product avoids the one place the log form
//   breaks: at full closure `A = -R_d < 0` and `log A` does not exist, which is
//   exactly the endpoint hazard 11 requires to be a candidate. The link is the
//   right way to STATE this cost and the wrong way to SEARCH for it.
//
// JW26 -- the same product objective with a LINEAR reduction factor. Everything
//   under SOX applies unchanged.
//   ⚠️ NOT THEIR FULL MODEL. Their psi_crit is a free parameter and their supply
//   is linear in the potential; here psi_crit is DERIVED from the stem curve and
//   the supply is that curve integral. So this is their objective on our
//   hydraulics -- which is what makes it comparable with SOX, the two `g` being
//   interpolations between the same two anchors -- and is NOT a reproduction of
//   the paper.
//
// ProfitMax -- the one curve needing setup before the search, and the one whose
//   degenerate exit ZEROES where every other curve EVALUATES. See the branch below.
//
// TF24_floor -- TF24's cost plus CF77's price, and the only DIFFERENCE objective
//   here whose marginal cost of water does NOT vanish at the wet end (LeastCost's
//   does not either, through its capacity term; see `lambda_LeastCost`). That is
//   the whole reason it exists: every conductance-loss curve in this package has
//   `lambda -> 0` as the drop closes, which prices water as free precisely when
//   it is abundant.
//   ⚠️ IT IS TF24 AT `lambda_o = 0`, AT TF24's OWN TRAITS. There is no second
//   hydraulic parameter: the cost reads `TF24_cost_scale` and `TF24_beta2`, so the
//   two curves differ by exactly one restriction and `lambda_o` is the whole of
//   what is being tested.
//   ⚠️ SO IT CAN BE WET-PINNED WHERE TF24 CANNOT, and the mechanism is CMax's.
//   dC/dpsi at the wet bound is `lambda_o*kmax*f(psi_soil)`, strictly positive for
//   any lambda_o > 0, where TF24's `|f'|` vanishes there. The closed interval is
//   doing real work here; read the note under CMax on comparing a pinned answer
//   against the bound, since the same step-in fraction applies.
//   ⚠️ IT INHERITS TF24's SECOND HUMP, which JS22 does not have. The marginal cost
//   is TF24's -- which FALLS at the dry end, the condition hazard 11 names -- plus
//   a constant, so `basin_scan_cells()` matters here for the same rows it matters
//   for TF24. Re-measure with the energy balance on before assuming otherwise.
//   ⚠️ ITS TWO TERMS ARE NOT INTERCHANGEABLE UNDER kmax/D RESCALING, and that is a
//   prediction rather than a fitting nuisance. `Theta~'/K` is invariant when kmax
//   and the vapour deficit move together; `lambda_o` is not. So this curve's
//   optimum sits BETWEEN TF24 (invariant) and CF77, by an amount that grows with
//   lambda_o -- which is what `test_tf24_floor_kmax_vpd_invariance` measures.
//
// ⚠️ `profit_` CARRIES THREE KINDS OF NUMBER because of these, and nothing
// enforces it: carbon for the five difference objectives, a dimensionless
// normalised profit for ProfitMax, and carbon times a dimensionless factor for the
// two products. A caller comparing `profit_` ACROSS curves is exposed; comparing
// each curve's argmax is fine. The guide's rule for `hydraulic_cost_` ("a third
// meaning behind the same name is how a reader quotes the wrong number") applies
// here and would be worth a separate field per objective kind.
// ============================================================================
template <Leaf::CostCurve K>
inline void Leaf::optimise_psi_stem_single() {
  clear_collar_solve_state();

  if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }
  check_cost_parameters<K>();

  const double psi_soil = supply_psi_soil_scalar();
  // The wet end of the interval: where the stem carries exactly the residual
  // flux, which is psi_soil itself at g_min = 0 (#182).
  const double psi_wet = stem_wet_bound(psi_soil);
  opt_psi_stem_ = psi_wet;

  // The dry end of the interval: psi_crit for a curve whose own domain ends
  // there, the end of the stem curve for every other (#183).
  const double psi_dry = dry_end<K>();

  // ⚠️ ONE CURVE NEEDS SETUP BEFORE THE SEARCH, AND ITS DEGENERATE EXIT DIFFERS.
  // ProfitMax divides by |A|max, so it is seeded here -- and where the normalisers
  // come out unusable there is no normalised objective to evaluate at the wet
  // end, which is why that arm ZEROES where every other curve EVALUATES. `PPFD_`
  // joins the test because a dark leaf has A = -R_d everywhere, so |A|max is
  // non-positive and the second bail-out would catch it anyway; testing it first
  // avoids seeding in order to discover that. Its normalisers are Sperry's own,
  // taken over the stream to psi_crit (k_crit is part of the published form),
  // so a wet end past psi_crit leaves k_span <= 0 and no objective, though the
  // search itself is no longer bounded there.
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
  } else if (psi_wet > psi_dry ||
             !(stem_curve_limit() - stem_curve_integral(psi_soil) >
               stem_capacity_floor * stem_curve_limit())) {
    // Drier soil than the stem can carry more than the cuticle's loss from
    // short of its dry end, or than the curve resolves a flux from
    // (stem_capacity_floor): one feasible potential, and this is it.
    profit_ = profit_psi_stem_for<K>(psi_wet, psi_soil);
    lambda_emergent_ = lambda_for<K>(psi_wet, psi_soil);
    return;
  }

  double profit_opt = 0.0;
  const auto objective = [&](double psi_stem) {
    return profit_psi_stem_for<K>(psi_stem, psi_soil);
  };
  // ONE solver, and NO caller-side branch on whether a derivative exists: the
  // derivative closure reports that itself, through `ok` at a shut-down or
  // compensation-point exit and as a non-finite value where a log link meets a
  // non-positive A, and the solver's own fallback handles it. Every curve has
  // dJ/dpsi on both gate states, so nothing here decides by curve. ⚠️ A
  // derivative that does not exist must be REPORTED that way, never thrown: the
  // objective is well defined wherever the derivative is not, and a throw from
  // inside the optimiser loses the whole solve. `basin_scan_cells()` returns a
  // scan whenever the gate is on, which is where the objective can carry two
  // basins.
  const auto dobjective = [&](double psi_stem, bool* ok) -> double {
    return dprofit_dpsi_stem<K>(psi_stem, psi_soil, ok);
  };
  // For a curve psi_crit does not bound, psi_crit is still the first bracket
  // tried, as on the collar route: where profit already falls there the
  // optimum lies inside, and the old interval finds it at the old cost. Not
  // under a basin scan, whose cells are set by the interval.
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
  // ⚠️ A SCAN KEEPS ITS CELL WIDTH, NOT ITS CELL COUNT. Its cells span the
  // interval, and past psi_crit the interval runs on to the end of the stem
  // curve, 2.6x as wide at the defaults; at the old count the cells widened with
  // it and stepped over the narrow basin a dark leaf has near its wet end
  // (measured: JS22 and TF dropped to the wet end at lower profit on every dark
  // energy-balance row). So the count grows with the interval against
  // [psi_wet, psi_crit], to at most 4x.
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
// MEDLYN STOMATAL-CONDUCTANCE MODEL (from develop #450)
// ---------------------------------------------------------------------------
// Standalone, R-callable coupling of the Medlyn (2011) optimal stomatal model
// to colimited photosynthesis. NOT invoked by the TF24 compute path, which
// optimises psi_stem directly via find_root_collar_psi; provided so the model
// remains available for R-level experimentation/comparison. beta_ is a soil-
// moisture stress factor in [0,1] from theta_/theta_w_/theta_fc_ (set in
// set_physiology from the default soil-moisture values).
// ===========================================================================
inline double Leaf::medlyn_model_gs(double assim_colimited_){

  const double beta_ = medlyn_beta();

  // The USO expression as published, which is a WATER VAPOUR conductance -- that is
  // what its diffusion ratio is doing there, and what `g1` is fitted against.
  //
  // ⚠️ The zero-deficit guard reads `atm_vpd_`, the DRIVEN value, not the `atm_vpd`
  // input field. It used to read the input while the expression below read the
  // working copy, so a leaf driven at zero deficit took the else-branch and
  // returned Inf from `sqrt(0)` -- the input field is still at its 2.0 default and
  // so never looked like zero. Same shape as the theta_/theta hazard.
  double gs_h2o;
  if (atm_vpd_ == 0) {
     gs_h2o = g0;
  } else {
     gs_h2o = g0 + medlyn_H2O_CO2_ratio*(1 + (g1*beta_)/sqrt(atm_vpd_))*(assim_colimited_/(ca_*(1/umol_per_mol_to_Pa_)));
  }

  // ...converted to the CO2 basis every other conductance in this class is on, so
  // that `A == stom_cond_CO2_*(ca_ - ci_)` holds for this route as it does for the
  // optimality routes. Storing the water-basis value here was #129. Dividing the
  // whole expression, `g0` included, is deliberate: Medlyn's `g0` is a water-basis
  // residual conductance too.
  medlyn_model_gs_ = gs_h2o / medlyn_H2O_CO2_ratio;
  return medlyn_model_gs_;
}

// Supply==demand residual for the Medlyn solver, as a function of ci (x, Pa).
// It is the difference between the Medlyn optimal stomatal conductance and the
// diffusion-implied conductance, *multiplied through by (ca_ - x)* so it stays
// finite across the whole [gamma*, ca_] bracket -- the raw gs difference has a
// 1/(ca_-x) singularity at the upper end. The root (zero crossing) is identical
// to that of the raw difference for x < ca_:
//   gs_medlyn*(ca_-x) - gs_coupled*(ca_-x),  where
//   gs_coupled*(ca_-x) = assim * (atm_kpa_*kPa_to_Pa) / 1e6.
//
// ⚠️ BOTH TERMS ARE ON THE CO2 BASIS. `medlyn_model_gs` returns a CO2 conductance
// now, so the CO2 flux on the right carries no diffusion ratio. This residual used
// to compare a water-basis conductance against a ratio-scaled flux, which is the
// same equation scaled by that ratio -- so THE ROOT IS UNCHANGED and `ci_`,
// `assim_colimited_` and everything downstream keep their values. Only the stored
// conductance moves (#129). Scaling both sides is what makes that true; changing
// one of them alone would move the operating point.
inline double Leaf::medlyn_stom_cond_minus_coupled_stom_cond(double x) {
  const double assim_colimited_x_ = assim_colimited(x);
  medlyn_model_gs_ = medlyn_model_gs(assim_colimited_x_);
  return medlyn_model_gs_ * (ca_ - x)
         - assim_colimited_x_ * (atm_kpa_ * kPa_to_Pa) / 1e6;
}

// Solve for the leaf-internal CO2 (ci) at which the Medlyn optimal stomatal
// conductance balances the diffusion-implied conductance, with Brent's method
// (util::uniroot). This replaces an earlier golden-section search on
// 1/|gs difference|, which could lock onto a spurious interior maximum
// (observed at high VPD).
//
// The residual is not monotone on [gamma*, ca_]: it has a hump and is negative
// at BOTH ends (near gamma* assimilation is below the compensation point;
// near ca_ the diffusion conductance diverges), so the interval can contain a
// spurious sub-compensation root as well as the meaningful Medlyn root. We
// therefore first locate the residual's maximum (Brent minimiser on -residual),
// then root-find on [argmax, ca_], which isolates the physically meaningful
// high-ci operating point in every case (including g0 == 0).
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

// Medlyn's CLOSED FORM for chi = ci/ca, which is the same model as the numerical
// route above rather than a different one: the coupled solve converges on this
// expression exactly as g0 -> 0 (measured, chi agreeing to every printed digit at
// g0 == 0, against a gap of 0.036 at the default g0 = 0.022). The gap IS the g0
// term, which the closed form drops along with Gamma* and respiration. Use this
// when a non-iterative chi is wanted and that approximation is acceptable.
//
// ⚠️ beta MULTIPLIES g1 HERE TOO. It did not, so this route was invariant to soil
// moisture while the coupled route moved chi from 0.681 to 0.446 over
// beta in [0.2, 1] -- the two entry points were different models rather than two
// forms of one. The substitution is the same one `medlyn_model_gs` makes, so at
// beta == 1 this reproduces the previous value exactly and only a soil-moisture
// sweep moves.
inline void Leaf::solve_medlyn_ci_analytical(){

  const double beta_ = medlyn_beta();
  const double g1_beta = g1 * beta_;
  ci_ = ca_ * (g1_beta/(g1_beta + sqrt(atm_vpd_)));
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  set_medlyn_shutdown_if_closed();
  return;
}

// The soil-moisture factor, the linear beta of the beta-factor literature,
// CLAMPED TO [0, 1]. Unclamped, a soil wetter than field capacity scaled g1 ABOVE
// its fitted value (theta 0.6 with field capacity 0.5 and wilting point 0.2 gave
// beta 1.33), and one drier than the wilting point gave a negative beta, so a
// negative g1 term the shutdown exit then had to catch. Between the two the
// factor is unchanged, so every value inside [theta_w_, theta_fc_] is
// bit-identical. CF77's soil price clamps the same factor at 1 for the same
// reason, and divides by it, so it has no use for the lower clamp.
inline double Leaf::medlyn_beta() const {
  const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
  return std::min(std::max(beta, 0.0), 1.0);
}

// A negative conductance is not an operating point. It arises where the USO
// expression's `A/ca` term goes negative faster than `g0` holds the stomata open --
// a leaf below its compensation point with a small or zero `g0` -- and the
// numerical route's "closest approach" exit would otherwise return it silently. At
// g0 == 0 and beta = 3.3e-4 that was gs = -3.95e-06, A = -9.7e-04.
//
// The answer is a closed state: no water moves, ci sits at the compensation point,
// and the leaf still respires, so assimilation is -R_d_ rather than zero.
// ⚠️ THIS IS NOT THE OPTIMALITY ROUTES' SHUT STATE. Those evaluate the leaf at
// the wet bound through the ordinary evaluator (#182), netting zero in light
// at g_min = 0, so in light this route still steps from a near-zero A to -R_d_
// as it shuts. Clamped rather than
// refused because a hot or dry leaf reaching it is a legitimate model state, unlike
// the parameter errors `check_cost_parameters` rejects.
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

