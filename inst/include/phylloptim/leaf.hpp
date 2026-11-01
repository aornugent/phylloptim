// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_HPP_
#define PHYLLOPTIM_LEAF_HPP_

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
#include <phylloptim/damage.hpp>

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
// The class and its construction are here; the other members are defined
// in the leaf_*.hpp headers by subsystem, and optimise_closed in
// closed_form.hpp. Include <phylloptim/leaf_model.hpp>, which includes the
// leaf_*.hpp headers: a member used and not defined in a translation unit is
// a link error.

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

// "Hazard N" in these headers is the numbered hazard in the repository's
// .claude/hazards.md, which is not shipped with the package.

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
  // in closed_form.hpp, which <phylloptim.hpp> includes and leaf_model.hpp does
  // not: including only leaf_model.hpp and calling it is a link error.
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
  // E_from_Soil_to_Root_Collar and dE_from_soil_dpsi_collar from one walk where
  // the path has one, each bit for bit what the two calls give; returns the slope.
  double E_and_slope_from_soil(double T_collar, const std::vector<double>& psi_soil) {
    require_suction_vector(psi_soil, "E_and_slope_from_soil");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        return roots_.uptake_and_slope_at(T_collar, psi_soil, soil_consumption_,
                                          E_up_);
      default:
        single_.uptake_at(T_collar, psi_soil, soil_consumption_, E_up_);
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
  // (notes above optimise_psi_stem_single in leaf_stem_routes.hpp): seeds |A|max
  // and the conductance span.
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

} // namespace phylloptim

#endif
