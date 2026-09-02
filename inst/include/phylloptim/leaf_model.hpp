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
#include <phylloptim/graft.hpp>

#include <odelia/interpolator.hpp>
#include <odelia/with_slope.hpp>
#include <odelia/tangent.hpp>
#include <odelia/implicit_node.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <XAD/XAD.hpp>

namespace phylloptim {

// --- the parameter enumeration, which R indexes into --------------------------
//
// The fifteen traits in `Leaf::set_traits`' argument order, then the two
// quantities a calibration fits that are not traits: the conductance driver and
// the single-potential path's series resistance.
//
// ⚠️ R INDEXES THESE POSITIONS, so a reordering silently differentiates the wrong
// parameter. `test-gradient-batch.R` reads the names back out of C++ and compares
// them with R's, so the two cannot drift apart without a failure.
inline constexpr int n_traits = 15;
inline constexpr int n_pars = 19;

// Every index by name, so nothing below indexes `theta` with a bare integer.
// The first `n_traits` are `set_traits`' arguments in its order, which is also
// `leaf_traits()`'; the two non-traits follow and take a relative step.
inline constexpr int par_vcmax_25 = 0;
inline constexpr int par_stem_c = 1;
inline constexpr int par_stem_P50 = 2;
inline constexpr int par_root_c = 3;
inline constexpr int par_root_P50 = 4;
inline constexpr int par_TF24_beta2 = 5;
inline constexpr int par_jmax_25 = 6;
inline constexpr int par_a = 7;
inline constexpr int par_curv_fact_elec_trans = 8;
inline constexpr int par_curv_fact_colim = 9;
inline constexpr int par_TF24_cost_scale = 10;
inline constexpr int par_R_d_25 = 11;
inline constexpr int par_JS22_gamma = 12;
inline constexpr int par_CMax_a = 13;
inline constexpr int par_CMax_b = 14;
// ⚠️ THESE MOVE WHENEVER A TRAIT IS ADDED, and bumping them is the whole cost.
// They are the non-traits and they sit AFTER the contiguous trait block, which is
// a readability convention rather than a constraint now: R's
// `.gradient_theta_matrix()` addresses EVERY column by name, including these.
// It used to take the traits as "everything but the last two", which is what made
// the ordering load-bearing; that is fixed. What is still load-bearing is the
// ORDER ITSELF -- R passes integer positions into this enumeration, so appending
// is safe and reordering silently differentiates the wrong parameter, and
// `test-gradient-batch.R` compares this enumeration against R's copy.
inline constexpr int par_kmax = 15;
inline constexpr int par_resistance = 16;
// Cowan-Farquhar's prescribed marginal value of water. A pure APPEND after the two
// existing non-traits, which is only safe because R addresses theta's non-trait
// columns by NAME rather than by position -- a positional rule ("everything but
// the last two") reads this as `resistance`.
//
// ⚠️ AVAILABLE FOR ONE MODEL. It is CF77's only parameter and every other curve's
// lambda is EMERGENT, derived from that curve's own parameters rather than set. So
// `.gradient_available_pars()` offers it only for CF77 and refuses it elsewhere,
// naming the model -- the same treatment `resistance` gets on the wrong supply path.
inline constexpr int par_CF77_lambda = 17;
// TF24_floor's price of water at zero transpiration, on exactly the same footing:
// an append after the non-traits, available for ONE model, and refused elsewhere
// by `.gradient_available_pars()` with the model named.
//
// ⚠️ IT IS THE SECOND MODEL-SPECIFIC SLOT, so "the CF77 one" has stopped being a
// safe way to talk about this class. R's `.gradient_model_pars()` is the single
// table that says which model owns which slot; there is no second copy here.
inline constexpr int par_TF24_floor_lambda_o = 18;


// The differentiable parameters, as one array indexed by the enumeration above.
//
// THIS IS THE INPUT PACK, and it is deliberately not a struct of named fields.
// A struct would be a fourth spelling of an order that set_traits, leaf_traits()
// and par_names() already spell three times, and every field would have to be
// listed again in a walk for seeding -- which is how a parameter gets added to
// the model and silently omitted from its gradient. Indexed by `par_*` the
// access reads the same (`pars[par_kmax]`, not `pars[15]`) and there is one list.
//
// ⚠️ IT CARRIES THE TRAITS, NOT THE QUANTITIES DERIVED FROM THEM. stem_b and
// psi_crit are functions of (stem_P50, stem_c) and are derived where they are
// used, so a seed in stem_P50 reaches them by the chain rule rather than by a
// caller remembering to move three numbers together. Passing the derived pair
// instead is what lets a curve be stated twice and the two disagree.
template <class S> using leaf_pars = std::array<S, n_pars>;

class Leaf {
public:
  //anonymous Leaf function as in canopy.h
  Leaf();
  
  Leaf(double vcmax_25, 
       double stem_c, 
       double stem_P50, 
       double root_c,
       double root_P50,
       double TF24_beta2, 
       double jmax_25, 
       double a, 
       double curv_fact_elec_trans, 
       double curv_fact_colim,
       double GSS_tol_abs,
       double vulnerability_curve_ncontrol,
       double ci_abs_tol,
       double ci_niter,
      double TF24_cost_scale);

  odelia::interpolator::Interpolator transpiration_from_psi;
  odelia::interpolator::Interpolator psi_from_transpiration;

  // The `stem_b` the two splines above were built at, which is normally just
  // `stem_b` -- and is not, while a gradient is perturbing it.
  //
  // ⚠️ THIS IS A HOMOGENEITY RESULT, AND IT IS WHAT MAKES A GRADIENT IN stem_b
  // FREE. The cumulative integral obeys
  //
  //     G(psi; s*b, c) = s * G(psi/s; b, c)
  //
  // because G(psi; b, c) = b * g(psi/b; c) for g(u; c) = integral of
  // exp(-u^c) -- b enters only as a scale on both axes. The knot grid scales
  // with it too (psi_max = b*log(100)^(1/c)), so the identity holds for the
  // SPLINE and not merely for the integral it approximates: measured, the
  // rescaled spline reproduces a rebuilt one to 0-3e-16.
  //
  // So moving stem_b needs no rebuild at all -- just the existing spline
  // evaluated at a rescaled argument. That matters because a rebuild is 11.9 us
  // of incomplete gammas plus 3.1 us per interpolator, which is the entire cost
  // of a gradient in stem_b.
  //
  // ⚠️ There is NO SUCH IDENTITY FOR stem_c, and the obvious substitute -- read
  // G from its closed form instead of the spline -- was built, measured and
  // rejected: it differentiates a slightly different model and disagrees with the
  // spline's own derivative by 3e-4. stem_c rebuilds.
  //
  // Two rules keep this from becoming a stale-state bug of the kind hazard 8
  // records: ONLY the four stem_curve_* accessors may read the splines, so one
  // scale factor is sufficient; and `setup_transpiration()` resets it, with
  // `set_traits()` forcing a rebuild while it is displaced, so no route out of
  // the perturbed state leaves it set.
  double stem_b_spline_ = 0.0;

  // INSTRUMENTATION, not state: how many times the stem curve has been built.
  // Nothing in the model reads it, and `set_traits()` deliberately does NOT reset
  // it -- it counts events over an object's whole life, so a caller takes
  // differences.
  //
  // It exists because a rebuild is invisible to every other instrument. The spline
  // is a pure function of (stem_b, stem_c), so a rebuild at an unchanged pair is
  // bit-identical to not rebuilding: the only thing it costs is 11.9 us, and a
  // timing assertion for one rebuild per observation would drift with the machine
  // and still pass on a fast day. This is the countable form (#74), which is what
  // the cost rules in the developer guide ask for wherever one exists.
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
  // Driest collar potential the supply path can be asked about, positive
  // magnitude. For roots it is the root vulnerability limit; a constant-
  // conductance path has no such limit, so the stem's psi_crit binds instead.
  double supply_psi_crit() const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.root_psi_crit;
      default:                     return psi_crit;
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
  // in terms of them: the splines, the brackets, the cost kernels.
  double stem_b;
  double psi_crit;
  double TF24_beta2;
  double jmax_25;
  double a;
  double curv_fact_elec_trans; // unitless - obtained from Smith and Keenan (2020)
  double curv_fact_colim;
  // Still a settable control, and it still has two jobs after the root-find replaced
  // the collar golden-section search: prepare_collar_solve's "this interval is too
  // narrow to solve over" threshold, and the tolerance of the golden-section
  // FALLBACK maximise_profit_over_collar drops to when neither bracket endpoint
  // has a usable gradient. It does NOT reach the stem route, which refines by a
  // root-find at collar_root_tol, and it no longer sets how well the reported
  // operating point is determined -- collar_root_tol does.
  double GSS_tol_abs;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
  double TF24_cost_scale;

  // Tolerance on the collar potential for the profit-maximising root-find.
  // Hard-coded rather than a control field, for the same reason
  // find_root_psi's 1e-4 and psi_stem_to_ci's 1e-7 are: it is a property of the
  // solve, not a knob. 1e-12 sits ~8 orders below the 1e-4 at which this package
  // calls a difference real, and the cost of going there from GSS_tol_abs's 1e-3
  // is a handful of evaluations, because TOMS748 is superlinear.
  static constexpr double collar_root_tol = 1e-12;


  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
  std::vector<double> soil_consumption_;
  double profit_;
  double psi_stem;
  // ⚠️ THE ONLY TWO INPUTS IN THIS BLOCK, and the reason they carry their default
  // here rather than in setup_clean_leaf(). Everything around them is derived
  // state or a solved output, which setup_clean_leaf() exists to wipe; these are
  // the Cowan-Farquhar marginal value of water, supplied by the CALLER and read
  // by profit_psi_stem_CF77. Wiping an input on the caller's behalf
  // makes whether a prescribed value survives depend on which of two
  // interchangeable-looking re-driving calls comes next, so neither clears it.
  //
  // An in-class initialiser, not a line in setup_clean_leaf(), so a fresh Leaf
  // still reads the NA sentinel: the member is otherwise uninitialised, and
  // "never reset" must not become "never initialised".
  // ⚠️ IT IS AN INPUT AND NOTHING ELSE. No model code assigns it. The units are
  // carbon per unit TRANSPIRATION, which is the Cowan-Farquhar marginal value of
  // water and the convention the lambda literature uses -- NOT the same quantity
  // as ProfitMax's normaliser ratio, which is per unit CONDUCTANCE.
  double CF77_lambda_ = util::na_value;         // umol C (kg H2O)^-1
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
  // `lambda_o` of the TF24_floor cost curve, and an INPUT on exactly the same footing
  // as `CF77_lambda_` above: same units, same meaning, no default, never assigned
  // by model code, and not cleared by setup_clean_leaf().
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
  double Tair_;  // air temperature, deg C (reinterprets the leaf_temp driver)
  double Rn_;    // net radiation at the leaf, W m^-2
  double ra_;    // aerodynamic (boundary-layer) resistance, s m^-1
  // Gate for the PM leaf energy balance. Default OFF: today's path runs
  // (prescribed leaf_temp, single-shot cached Arrhenius). R-settable (#523 full
  // cut) so TF24 (via pars.use_energy_balance) and the leaf-level demo can
  // turn PM on; default preserves backward compatibility.
  bool use_energy_balance_ = false;
  // Set by psi_stem_to_ci when the energy-balance fallback fires: the leaf is so
  // hot that assimilation is negative across the whole [gamma*, ca] bracket, so
  // there is no supply==demand root and ci is placed AT the compensation point
  // instead. The residual g is then NOT zero, which voids the implicit-function
  // theorem dprofit_at_collar_psi is built on -- see the branch there.
  bool ci_at_compensation_point_ = false;
  // Boundary-layer inputs for ra = C_ra*sqrt(d/U0) (doc 4.1). d is a per-strategy
  // trait (set from pars.d in prepare_strategy); wind_speed_ is the per-timestep
  // above-canopy driver (set from the environment before set_physiology). Both
  // R-settable so a bare Leaf can exercise the wind model; if either is unusable
  // set_physiology falls back to the fixed ra. Only read on the PM path.
  double d_ = 0.05;          // characteristic leaf dimension, m
  double wind_speed_ = 2.0;  // above-canopy wind speed U0, m s^-1

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
  // ⚠️ THE DEFAULT IS NOT A LITERATURE VALUE. It is order-of-magnitude matched to
  // TF24 at this package's defaults, and there is no single value that matches
  // more than one soil potential: the gamma reproducing TF24's cost at its own
  // optimum runs 0.287 (psi_soil 0.5 MPa) to 1.761 (3.0 MPa), a 6.1x range,
  // because TF24's cost tracks the ABSOLUTE potential through (1 - f) while this
  // one tracks the DROP -- and across a drydown those move in opposite directions.
  // That spread is the models' structural difference, not a fitting problem.
  double JS22_gamma = 1.0;             // umol C m^-2 s^-1 MPa^-2
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
  // ⚠️ NEITHER DEFAULT IS A LITERATURE VALUE. Measured: with CMax_b = 0 the CMax_a
  // reproducing TF24's cost at TF24's own optimum runs 0.406 (psi_soil 0.5 MPa) to
  // 0.807 (3.0 MPa) -- a 2.0x range, against JS22's 6.1x, which is the absolute-
  // versus-drop distinction showing up as a number. 0.6 is mid-range of that band.
  double CMax_a = 0.6;                 // umol C m^-2 s^-1 MPa^-2
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
  double opt_ci_;
  double count;
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
  double medlyn_model_gs_;  // mol CO2 m^-2 s^-1
  double theta_w_;          // current soil water content at wilting point (m^3 m^-3)
  double theta_fc_;         // current soil water content at field capacity (m^3 m^-3)
  double theta_;            // current soil water content (m^3 m^-3)

  // 1-entry memo for transpiration(). Within a single root-collar/profit solve
  // leaf_specific_conductance_max_ and the transpiration spline are fixed, so
  // supply-side transpiration depends only on (psi_stem, psi_upstream). That
  // pair is queried repeatedly with identical values per profit evaluation
  // (psi_stem_to_ci -> stom_cond_CO2 -> transpiration, then transpiration again
  // for transpiration_). Caching the last result avoids the redundant spline
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
  //     still the right way to change a trait (the vulnerability splines and the
  //     solved point need clearing too), but the silent-wrong-number failure mode
  //     is gone.
  //   * the cost is 17 double comparisons per set_physiology() call, i.e. per
  //     driver set, NOT per inner solve iteration. Measured: within run-to-run
  //     noise on bench_solve.
  static constexpr int photo_temp_key_size = 22;
  std::array<double, photo_temp_key_size> photo_temp_cache_key_{};
  std::array<double, photo_temp_key_size> photo_temp_key() const;
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
  //   * the two vulnerability SPLINES. stem_b/stem_c build transpiration_from_psi
  //     and psi_from_transpiration; root_b/root_c build the root curve. Both are
  //     pre-integrated at construction, so a bare `l.stem_b = x` leaves the hot
  //     path reading the previous curve's integral while proportion_of_conductivity
  //     reports the new one -- two answers to the same question.
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
  // So: the invariant checks run, the splines
  // are rebuilt only when the curve that owns them actually moved, and
  // setup_clean_leaf() puts the object back in its just-constructed state -- which
  // resets both caches and requires set_physiology() before the next solve, exactly
  // as a fresh Leaf would. That last part is not conservatism: the derived
  // photosynthetic parameters really are unknown until the drivers are re-supplied.
  void set_traits(double vcmax_25, double stem_c, double stem_P50,
                  double root_c, double root_P50,
                  double TF24_beta2, double jmax_25, double a,
                  double curv_fact_elec_trans, double curv_fact_colim,
                  double TF24_cost_scale, double R_d_25, double JS22_gamma,
                  double CMax_a, double CMax_b);

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
  // So three of the eight cost curves change shape, not just a bound that seldom
  // binds. And there is a ceiling: the vulnerability splines are built out to P99
  // (`vulnerability_conductivity_floor`), and the collar path recovers psi_stem by
  // INVERTING the stem integral, whose range ends there -- measured, asking 0.03%
  // past it throws `stem_curve_domain`. A threshold at P99 leaves that inversion no
  // headroom, and past P99 the spline domain has to move too, which redistributes
  // every knot and therefore every number the model reports.
  //
  // Smoothness is the wrong reason to reach for this in any case: moving a bound
  // moves the kink rather than removing it. What removes it is a cost that DIVERGES
  // at the bound, which is exactly why no product objective can be boundary-crit.
  static constexpr double k_crit_fraction = 0.05;

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
  // splines to rebuild by comparing them against the ones standing.
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

  // The #38 boundary, and the STRONGER of the two: `psi_crit` must lie inside the
  // stem vulnerability spline's domain, which `stem_b`/`stem_c` set and `psi_crit`
  // does not enter. The knot grid stops at P99 =
  // vulnerability_psi_max(stem_b, stem_c) and setup_transpiration disables
  // extrapolation, while every solve evaluates the stem curve AT `psi_crit` -- the
  // dry bracket bound, the shutdown exit's cost, and E_column's feasibility test
  // all do. So `psi_crit > P99` is not a configuration that sometimes works: it
  // throws, out of the interpolator, in a message naming neither trait.
  //
  // ⚠️ A caller is much more likely to violate this than the positivity checks,
  // because `psi_crit` LOOKS independent of the curve and is not. What the
  // defaults actually say is that it is P95:
  //
  //   stem_b = 3.898245, stem_c = 2.680147  ->  P95 = 5.870283 = psi_crit
  //
  // to six decimal places, against P99 = 6.891842. That relationship is the thing
  // a caller fitting measured vulnerability curves needs and cannot find anywhere
  // else, so the message quotes the P95 that WOULD work rather than only the bound
  // that failed.
  //
  // Separate from check_psi_magnitudes rather than folded into it, because that
  // one is `static` and reached by name from plant's generated glue and from the
  // CI consumer program: extending its signature is an API break where adding a
  // function is not.
  //
  // ⚠️ NOT applied to the root curve, and that is a finding rather than an
  // omission. #38 assumed `root_psi_crit` carried the same latent constraint;
  // since #77 bounded the root curves past their last knot it does not throw --
  // `root_vuln_at` CLAMPS its argument to the last knot instead. So a
  // `root_psi_crit` past the root P99 silently reports the floor conductivity
  // rather than failing, which is a different defect and is #85's question, not
  // this check's.


  void setup_transpiration(double resolution);

  // The stem cumulative-vulnerability integral G and its inverse, as the FOUR
  // operations the model actually performs on them. Every read of
  // transpiration_from_psi / psi_from_transpiration goes through these, which is
  // what lets `stem_curve_closed_form_` be a single flag rather than a condition
  // repeated at eight call sites.
  //
  // At stem_b == stem_b_spline_ they are the splines, and bit-identically so --
  // the scale is exactly 1.0, and dividing by it is the identity. Otherwise they
  // apply the homogeneity identity documented at stem_b_spline_, with
  // s = stem_b / stem_b_spline_:
  //
  //   G(psi)     = s * G_spline(psi/s)
  //   G'(psi)    =     G'_spline(psi/s)
  //   G^-1(w)    = s * G^-1_spline(w/s)
  //   (G^-1)'(w) =     (G^-1)'_spline(w/s)
  //
  // The two derivative lines have no leading `s` on purpose: differentiating
  // s*G(psi/s) with respect to psi cancels it.
  //
  // The two non-derivative reads take an optional `caller` label, which is
  // appended to an out-of-domain failure. See eval_stem_curve for why the two
  // splines cannot identify themselves and why the caller has to.
  double stem_curve_integral(double psi, const char* caller = nullptr) const;
  double stem_curve_integral_deriv(double psi) const;
  double stem_curve_integral_inverse(double w, const char* caller = nullptr) const;
  double stem_curve_integral_inverse_deriv(double w) const;

  // Domain-guarded read behind the two accessors above. The stem curve is the
  // only interpolator in this file built with extrapolation DISABLED (the root
  // vulnerability pair clamps instead -- see setup_root_vulnerability), so it is
  // the only one a lookup can throw on. odelia's message names the point and the
  // domain but cannot name WHICH spline, because it does not know: there are two
  // here, they are inverses of each other, and they carry different units, so
  // "u = 7.5 beyond the upper end" is ambiguous in exactly the way that matters.
  // Nor can it name the caller -- the same spline is read from four places, and
  // localising plant#576 came down to which.
  static double eval_stem_curve(const odelia::interpolator::Interpolator& spline,
                                double u, double scale, const char* spline_name,
                                const char* arg_name, const char* caller);

  // Move stem_b WITHOUT rebuilding the stem vulnerability spline, by the
  // homogeneity identity above. stem_c is not accepted: it has no such identity.
  //
  // ⚠️ FOR DERIVATIVE WORK ONLY. It leaves the splines describing a different
  // stem_b, which is sound only because every read goes through the four
  // accessors, and it deliberately does NOT clear the solved operating point --
  // re-seating the physiology is exactly the cost being avoided -- so whatever
  // reads the outputs afterwards must write them first. `set_traits()` is the
  // way back, and forces a rebuild.
  void perturb_stem_P50(double stem_P50_new);
  // Forwards to roots_.setup_vulnerability. Kept on Leaf because it is part of
  // the published construction sequence (see the umbrella header) and plant's
  // bindings name it.
  void setup_root_vulnerability(double resolution) {
    roots_.setup_vulnerability(resolution);
  }
  // Forwards to phylloptim::cumulative_vulnerability_integral, which now lives in
  // vulnerability.hpp because it is shared by the stem and the root curves and
  // so belongs to neither. The parameters are deliberately neutral names, not
  // stem_*: this is called with (stem_b, stem_c) from setup_transpiration and
  // with (root_b, root_c) from setup_root_vulnerability.
  void build_cumulative_vulnerability_integral(double weibull_b, double weibull_c,
                                               double resolution,
                                               std::vector<double>& x,
                                               std::vector<double>& y_integral) {
    cumulative_vulnerability_integral(weibull_b, weibull_c, resolution, x,
                                      y_integral);
  }
  void setup_clean_leaf();

  // Absolute tolerance for the direct quadrature in
  // transpiration_full_integration (leaf/quadrature.hpp). Not used anywhere on
  // the production path, which reads the pre-integrated spline instead.
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
  // Which cost curve a psi_stem derivative differentiates. The cost enters the
  // chain through exactly ONE quantity -- dC/dpsi_stem -- so this selects that
  // and nothing else.


  // A value and its response in the collar potential, travelling as one pair.
  template <class T> using pair = odelia::with_slope<T>;

  // What the leaf hands over: the objective plant bills, and the water it drew
  // from each layer. At any scalar, because both re-enter plant's chain.
  template <class T>
  struct LeafOutputs {
    T profit{};
    std::vector<T> uptake;
  };

  // The same two, at double, plus whether the collar could be evaluated at all.
  // DERIVES rather than repeating them.
  //
  // The flag is not an error channel and should not become one: an infeasible
  // collar is an EXPECTED outcome on a sweep -- 1536 of the golden file's 5184
  // rows refuse -- so unwinding for it would be the wrong shape. This is the same
  // convention the `_checked` readers use, under a name that says which collar.
  struct FixedCollarEval : LeafOutputs<double> {
    bool feasible = false;
  };

  // The soil's draw at ONE collar, with the collar it was taken at, so a flux
  // from one and a slope from another cannot be paired. check_draw enforces that
  // at every read rather than trusting the caller.
  //
  // ⚠️ TWO DIFFERENT PAIRINGS, DELIBERATELY. `flux` and its collar slope are both
  // active, so they are a pair. The per-layer draws are active but their collar
  // slopes are SUPPLIED at double -- taking the supply again at a live collar
  // would record the whole of it a second time -- so those are implicit_value's
  // shape, not the pair's, and flattening both into one would lose that.
  template <class S>
  struct SupplyDraw {
    double at = 0.0;
    pair<S> flux;
    std::vector<S> uptake;
    std::vector<double> duptake_dp;
  };

  // The only way to make one, so the draw records the collar it was taken at and
  // its values cannot come from different ones.
  //
  // ⚠️ THE COLLAR IS PASSIVE HERE, whatever the caller hands in. This is the
  // supply at a POINT; where the point itself moves is the graft's business, and
  // taking it live would record the search that placed it.
  template <class S>
  SupplyDraw<S> supply_draw_at(const S& collar,
                               const SupplyAt<S>& supply) const {
    const std::size_t n = static_cast<std::size_t>(supply_n_layers());
    SupplyDraw<S> d;
    d.at = odelia::util::to_passive(collar);

    // ⚠️ NO DRAW AT A SHUTDOWN, and not because it would be zero. The collar is
    // held at the stem's critical potential there, which is past where the uptake
    // model answers at all -- the integral is exactly zero and E_up is not a
    // number, so asking costs a stop rather than a wrong row. Every consumer
    // reads zero at that kind anyway. ShadeDeath is NOT this case: it sits on the
    // wet bound, where E_up is zero in value and its rows are the bound's own.
    if (operating_point_kind_ == OperatingPointKind::HydraulicShutdown) {
      d.uptake.assign(n, S(0.0));
      d.duptake_dp.assign(n, 0.0);
      return d;
    }

    const S held = S(d.at);
    d.uptake.assign(n, S(0.0));
    roots_.template uptake_at<S>(held, supply, d.uptake, d.flux.value);
    d.flux.slope = roots_.template duptake_dpsi<S>(held, supply);

    // The per-layer collar slopes, SUPPLIED at double: the draws already carry
    // their own rows, and what is missing is only the channel through the collar
    // moving. Taking the supply again at a live collar would record the whole of
    // it a second time for the same numbers.
    std::vector<double> soil_at;
    soil_at.reserve(supply.psi_soil.size());
    for (const S& v : supply.psi_soil) {
      soil_at.push_back(odelia::util::to_passive(v));
    }
    roots_.duptake_dpsi(d.at, soil_at, d.duptake_dp);
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

  // ---- the differentiable surface's coordinates ------------------------------
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
    pair<T> sigma;
    pair<T> ci;
  };

  // What the marginal needs: those two, plus the collar itself and the stem's
  // flux. DERIVES rather than repeating the pair above -- written out separately
  // they were two structs differing by two members, and a sigma from one with a
  // ci from another was a thing that compiled.
  //
  // ⚠️ FOUR MEMBERS, NOT FIVE. What used to be a separate `dEup_dp` field IS
  // `transpiration.slope`: the stem's flux and the soil-to-collar conductance are
  // one quantity at two orders, which a flat struct of bare scalars could only say
  // in a comment. p.slope is 1 by construction.
  template <class T>
  struct CollarPoint : CollarCoords<T> {
    pair<T> p;
    pair<T> transpiration;
  };

  enum class CostCurve {
 TF24, CF77, JS22, CMax, SOX, JW26, ProfitMax, TF24_floor };

  // ⚠️ THE ONE ABSTRACTION THAT MAKES EVERY MODEL THE SAME MODEL.
  //
  // All seven maximise `h(A(psi)) - C(psi)`: a BENEFIT LINK `h` composed with
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
  //                                                   an argmax)
  //     Scaled    h(A) = A/|A|max     h' = 1/|A|max   ProfitMax
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
  static constexpr int n_cost_curves = static_cast<int>(CostCurve::TF24_floor) + 1;

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
  // every feasibility early-exit (shutdown / assim<0 / collapsed interval) by
  // setting the final operating point itself, and otherwise returns the feasible
  // collar-potential interval [bound_a, bound_b] (positive magnitudes). Returns
  // false when the operating point is already fully determined (caller is done),
  // true when there is a real interval to choose a collar potential within.
  template <CostCurve K>
  bool prepare_collar_solve(double& bound_a, double& bound_b);
  // Evaluate the leaf at a *given* root-collar potential (positive magnitude)
  // rather than optimising it: reuses prepare_collar_solve, clamps the target to
  // the feasible interval, and evaluates there (no golden-section search). Leaves
  // exactly the same outputs as find_root_collar_psi and returns profit_. Used by
  // TF24f's gradient-ascent acclimation (#525).
  double evaluate_root_collar_psi(double target_opt_root_psi);
  template <CostCurve K>
  double evaluate_root_collar_psi_for(double target_opt_root_psi);
  // Evaluate profit at a given root-collar potential (positive magnitude)
  // *assuming prepare_collar_solve has already run this step* (soil-side caches
  // built, feasible interval [bound_a, bound_b] known). Clamps the target into
  // the interval and sets the operating point (opt_psi_stem_, opt_root_psi_,
  // profit_), returning profit_. This is the post-prepare body of
  // evaluate_root_collar_psi, factored out so the centred finite-difference leaf
  // solve can share one prepare_collar_solve across its three profit evals (#530).
  template <CostCurve K>
  double profit_at_collar_psi(double target_opt_root_psi,
                              double bound_a, double bound_b);
  // Exact d(profit)/d(opt_root_psi) at a given root-collar potential (positive
  // magnitude), for TF24f's acclimation tracking (#525/#527). Combines
  // forward-mode AD for the analytic photosynthesis/cost algebra, the
  // implicit-function theorem at the psi_stem_to_ci root-find, and analytic
  // spline derivatives (Interpolator::deriv) for the smooth transport. Replaces
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
  // begin_solve() is a spline evaluation per soil layer.
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
  // whose objective is defined everywhere. Hence zeroed rather than evaluated.
  void set_profitmax_degenerate_state();
  template <CostCurve K> void check_cost_parameters();
  template <CostCurve K> double cost_deriv(double psi_stem,
                                           double psi_upstream);
  template <CostCurve K> double profit_psi_stem_for(double psi_stem,
                                                    double psi_upstream);
  template <CostCurve K> double lambda_for(double psi_stem,
                                           double psi_upstream);

  template <CostCurve K>
  double evaluate_psi_stem(double target_psi_stem);
  // dR_d/dT at a leaf temperature, by differencing the model's OWN temperature
  // block so the two cannot drift apart. Used on the compensation-point branch of
  // both dprofit routes, which is the whole reason it is a function: the two
  // carried an identical copy, and the copy is nine lines of save-and-restore
  // around a call with side effects, which is not a thing to keep two of.
  //
  // ⚠️ NOT const, and the caller sees no change. `update_temperature_dependent_
  // params` writes eight members; this saves them, perturbs, reads R_d_, and puts
  // all eight back. Adding a member to that block without adding it here leaves
  // the leaf perturbed by 1e-3 degrees after a gradient evaluation, silently --
  // the same failure mode `photo_temp_key()` guards for the cache.
  double respiration_temp_deriv(double leaf_temp);
  // The energy-balance correction to the above, zero when the gate is off. Kept
  // out of line so that adding it cannot change FMA contraction in the inlined
  // gate-off path; the derivation and the two sign checks are at the definition.
  double dprofit_energy_balance_term(double ci, double gc, double g_ci,
                                     double inv_atm, double gc_const,
                                     double A_prime, double dgc_dpsistem,
                                     double dgc_dpsi, double dpsistem_dpsi,
                                     double dT_dE, double Tleaf);
  // The profit-maximising collar potential within [bound_a, bound_b], by a
  // safeguarded root-find on dprofit == 0. Returns a bound when the
  // optimum is pinned to it, which is the case on 42 of the 240 feasible
  // golden-grid rows and so is a branch that has to be written rather than a
  // corner. Falls back to golden section if either endpoint has no usable
  // gradient -- see the definition for why that fallback should never fire.
  template <CostCurve K>
  double maximise_profit_over_collar(double bound_a, double bound_b);
  // The objective's value at a no-flow point, for the shut-down exits.
  template <CostCurve K> double no_flow_profit(double psi_stem);
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
  // Shut-down operating point used by the find_root_collar_psi early-exits: stem
  // held at psi_crit (no transpiration), paying only respiration + hydraulic
  // cost. Only opt_root_psi_ differs between the cases, so it is the argument
  // (a positive magnitude, like every other psi here).
  template <CostCurve K> void set_shutdown_state(double root_collar);
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit);
  double find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil);
  double E_column(double x, const std::vector<double>& psi_soil, double psi_leaf);
  double E_column_zero(double x, const std::vector<double>& psi_soil);
  
  // The two temperature responses, at any scalar. Both are exactly LINEAR in
  // ref_value and every other argument is fixed physiology, so on the active path
  // the response is a passive factor and the reference value is what carries the
  // row. Templated rather than factored into that product because the operation
  // order is what keeps the double instantiation bit-identical to the member it
  // replaced.
  template <typename T> T arrh_curve(T Ea, T ref_value, T leaf_temp) const;
  template <typename T>
  T peak_arrh_curve(T Ea, T ref_value, T leaf_temp, T H_d, T d_S) const;

  // --- Penman-Monteith leaf energy balance (minimal core; #523) ---------------
  // Recompute the temperature-dependent photosynthetic parameters (vcmax_,
  // jmax_, gamma_, ko_, kc_, R_d_, km_, electron_transport_) at a given leaf
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
  // Explicit leaf energy balance: Tleaf = Tair + (Rn - lambda*E) * ra / (rho*cp).
  // E is the hydraulically-pinned transpiration (kg H2O m^-2 s^-1); no PM
  // inversion and no A->E feedback, so this is a single algebraic forward pass.
  // `dT_dE`, when non-null, receives dTleaf/dE at the same point. It is an
  // out-parameter rather than a second function so that ONE clamp test decides
  // both the value and the slope: a separate slope function would return the
  // interior slope while the value was clamped, which is a silently wrong
  // derivative rather than an imprecise one, and the collar solve would then
  // converge to a point where dprofit is genuinely non-zero.
  double leaf_temp_from_E(double E, double* dT_dE = nullptr) const;

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;
  // Scalar-generic core of the above. See the `_kernel` block below the assim
  // declarations for why these exist.
  template <typename T> T proportion_of_conductivity_kernel(T psi) const;

  // supply-side transpiration for a given water potential gradient between leaves and soil, 
  // references setup_transpiraiton for values (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration(double psi_stem, double psi_upstream);
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
  // Why members templated on the scalar type, rather than free functions taking
  // their parameters explicitly (which is what #4 comment 1 proposed): a free
  // function needs seven or eight arguments threaded from the members at each call
  // site, and a transposed pair there is exactly the class of silent error the
  // replicas were. Reading the members directly makes that unrepresentable. It is
  // also the shape item 11's `Leaf<T>` wants, so it is a step rather than a
  // detour.
  //
  // Keep them PURE -- no writes to members. `hydraulic_cost_TF` caches into
  // `hydraulic_cost_`; its kernel must not, or the AD pass would write model state
  // while probing.
  template <typename T> T assim_rubisco_limited_kernel(T ci) const;
  template <typename T> T assim_electron_limited_kernel(T ci) const;
  template <typename T> T assim_colimited_kernel(T ci) const;
  template <typename T> T hydraulic_cost_TF_kernel(T psi_stem) const;

  // The same kernels against a caller's parameter pack, for the path that needs
  // rows in the traits.
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
  template <typename T>
  T assim_colimited_kernel(const T& ci, const leaf_pars<T>& pars) const;
  // dA/dci, by a TANGENT THROUGH THE KERNEL ABOVE rather than a slope written
  // out beside it. Upstream takes A' this way at double for exactly the reason it
  // matters here: a hand-written twin is a second definition of one function,
  // free to disagree with it by an amount that stays finite and plausible in
  // every column. At a nested tangent this is the same statement one order up.
  template <typename T>
  T assim_slope_at(const T& ci, const leaf_pars<T>& pars) const;

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
  template <typename T>
  T proportion_of_conductivity_kernel(const T& psi,
                                      const leaf_pars<T>& pars) const;
  template <typename T>
  T hydraulic_cost_TF_kernel(const T& psi_stem, const leaf_pars<T>& pars) const;

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
  template <class S>
  CollarCoords<S> collar_coords_at(double sigma_star, double ci_star,
                                   const S& collar, const SupplyDraw<S>& draw,
                                   const leaf_pars<S>& pars) const;

  // The stem's critical potential at any scalar. DERIVED from the trait pair, as
  // it is at double, which is what makes the old `no_gradient` question moot: it
  // is an ordinary expression the chain differentiates, not a parameter that has
  // to be declared not to move.
  template <class S>
  S psi_crit_at(const leaf_pars<S>& pars) const;

  // The objective at a placed operating point, and everything plant reads.
  //
  // ⚠️ TEMPLATED ON THE CURVE, IMPLEMENTED FOR TF24. The `if constexpr` chain
  // ends in a static_assert naming the rest, so a second curve is a row here
  // rather than a redesign -- and until something differentiates one, an arm
  // written for it would be a path with no caller and no reference to check it
  // against. plant's TF24_Strategy is the only consumer of reverse mode.
  template <CostCurve K, class S>
  S profit_at(const S& collar, const SupplyDraw<S>& draw,
              const leaf_pars<S>& pars) const;

  template <CostCurve K, class S>
  LeafOutputs<S> outputs_at(const S& collar, const SupplyDraw<S>& draw,
                            const leaf_pars<S>& pars) const;

  // M = dprofit/dcollar, as a SCALAR. What used to travel beside it -- the two
  // collar responses -- are the coordinates' own slope halves now, so there is
  // nothing left to return with it.
  template <CostCurve K, class S>
  S marginal_at(const S& collar, const SupplyDraw<S>& draw,
                const leaf_pars<S>& pars) const;

  // The pack as this leaf currently stands. The double path passes this where the
  // active path passes a seeded copy, so both reach the kernels the same way and
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
            lambda_emergent()};
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
    // feasible collar interval. 198 of the 288 golden grid points at 25 C -- 160 of
    // them at 40 C, where the optimum presses against the WET bound instead.
    Interior,
    // Constrained optimum at the psi_soil end of the feasible interval, just
    // inside root_zero_E (the collar at which uptake is exactly zero). The
    // objective is still climbing toward that bound, so the bound is the answer:
    // the leaf declines to open. 24 golden points at 25 C. The gradient here is
    // genuinely non-zero.
    //
    // ⚠️ THESE TWO ARE NAMED FOR THE BOUND THAT BINDS, NOT FOR THE CONDITIONS,
    // and an earlier `PinnedWet`/`PinnedDry` pair got read backwards for exactly
    // that reason. `psi_stem == psi_soil` is the least-tension end of the
    // bracket, so it reads as the comfortable case -- but it is what a leaf does
    // under STRESS, when no interior point earns its water. The golden grid says
    // so: 24 -> 80 rows as the leaf goes 25 -> 40 C, while the psi_crit end goes
    // 18 -> 0. Warming drives this one and removes the other. Do not reintroduce
    // a name whose reading inverts with temperature.
    BoundarySoil,
    // Constrained optimum at the psi_crit end, min(root_crit,
    // supply_psi_crit()) -- the most tension the stem is allowed. Worth pulling
    // that hard only when water is cheap relative to carbon, so this is the
    // UNstressed pin. 18 golden points at 25 C, none at 40 C. Gradient genuinely
    // non-zero.
    BoundaryCrit,
    // The feasible interval collapsed to a point (width <= GSS_tol_abs), so
    // feasibility DETERMINED the collar potential and nothing was optimised.
    // There is no free variable left to differentiate.
    Determined,
    // Shutdown on water: no collar potential both moves water and stays inside
    // the stem's and the root's critical potentials. The stem holds at psi_crit,
    // transpiration is zero, and the leaf pays respiration plus the hydraulic
    // cost there. 48 of the 288 golden grid points at every temperature, since it is
    // hydraulics rather than heat that forbids transpiration. The water response is zero;
    // the carbon response is not.
    HydraulicShutdown,
    // Shutdown on light: assim_max_ < 0, so gross assimilation at ci = ca cannot
    // cover dark respiration. Governed by PPFD and temperature, not by water, so
    // no rainfall sweep reaches it -- and the golden grid does not either (its
    // minimum assim_max_ is 3.71).
    ShadeDeath,
    // The collar potential was IMPOSED by the caller (evaluate_root_collar_psi /
    // profit_at_collar_psi), clamped into the feasible interval, not optimised.
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

  // Read-only on purpose: the tag is an output of the solve, and a settable one
  // would be a way to disagree with it. Same argument as `supply_kind_`'s two
  // entry points above, one step further.
  OperatingPointKind operating_point_kind() const {
    return operating_point_kind_;
  }
  static const char* operating_point_kind_name(OperatingPointKind kind);

private:
  // Written by every path out of the collar solve, and reset to Unsolved at the
  // top of prepare_collar_solve. Hazard 8 in the developer guide is why: `Leaf`
  // is a value member that plant reuses for every individual in a patch, so a
  // branch that declines to write this would leave the PREVIOUS plant's
  // classification -- a plausible answer about a different plant, which is the
  // worst failure shape available here. Defaulting the reset to Unsolved means a
  // path that forgets reports "unclassified" instead.
  OperatingPointKind operating_point_kind_ = OperatingPointKind::Unsolved;
};

// Human-readable tag. The switch has no default, so a missing name is a -Wswitch
// warning -- an error under the test suite's -Werror=switch, "unknown" elsewhere.
inline const char* Leaf::operating_point_kind_name(OperatingPointKind kind) {
  switch (kind) {
    case OperatingPointKind::Unsolved:          return "unsolved";
    case OperatingPointKind::Interior:          return "interior";
    case OperatingPointKind::BoundarySoil:      return "boundary-soil";
    case OperatingPointKind::BoundaryCrit:      return "boundary-crit";
    case OperatingPointKind::Determined:        return "determined";
    case OperatingPointKind::HydraulicShutdown: return "hydraulic-shutdown";
    case OperatingPointKind::ShadeDeath:        return "shade-death";
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
    a(0.30), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosythnthesis equatiom
    GSS_tol_abs(1e-3),
    vulnerability_curve_ncontrol(100),
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
      setup_transpiration(100); // arg: num control points for integration
      setup_root_vulnerability(100);
      setup_clean_leaf();
}

inline Leaf::Leaf(double vcmax_25, double stem_c, double stem_P50,
           double root_c,
           double root_P50,
           double TF24_beta2, double jmax_25,
           double a, double curv_fact_elec_trans, double curv_fact_colim, 
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
    a(a), //quantum yield of photosynthetic electron transport (mol mol^-1)
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
                             double jmax_25_, double a_,
                             double curv_fact_elec_trans_,
                             double curv_fact_colim_,
                             double TF24_cost_scale_, double R_d_25_,
                             double JS22_gamma_, double CMax_a_,
                             double CMax_b_) {
  check_psi_magnitudes(stem_P50_, stem_c_, root_P50_, root_c_);

  // Both scale parameters and both critical potentials fall out of the traits, so
  // there is nothing left for a consistency check to police.
  const double stem_b_ = weibull_b_from_P50(stem_P50_, stem_c_);
  const double root_b_ = weibull_b_from_P50(root_P50_, root_c_);
  const double psi_crit_ =
      weibull_psi_at_fraction(stem_b_, stem_c_, k_crit_fraction);
  const double root_psi_crit_ =
      weibull_psi_at_fraction(root_b_, root_c_, k_crit_fraction);

  // Which splines have to be rebuilt, decided BEFORE the assignment. Exact
  // equality is the right test and not a sloppy one: the spline is a pure
  // function of the pair, so an unchanged pair gives a bit-identical spline and
  // rebuilding it is pure cost. A trait perturbed by a relative 1e-08 -- what a
  // gradient loop does -- compares unequal and rebuilds, which is the case that
  // must not be missed.
  // ⚠️ The third clause is what makes set_traits the way back from
  // perturb_stem_P50. While the splines are built at a different stem_b, "the
  // parameters did not move" is not a reason to keep them -- it is exactly the
  // case where the equality test would conclude there is nothing to do and leave
  // the object rescaling forever.
  const bool stem_curve_moved = (stem_b_ != stem_b) || (stem_c_ != stem_c) ||
                                (stem_b != stem_b_spline_);
  const bool root_curve_moved =
      (root_b_ != roots_.root_b) || (root_c_ != roots_.root_c);

  vcmax_25 = vcmax_25_;
  stem_c = stem_c_;
  stem_P50 = stem_P50_;
  stem_b = stem_b_;
  psi_crit = psi_crit_;
  TF24_beta2 = TF24_beta2_;
  jmax_25 = jmax_25_;
  a = a_;
  curv_fact_elec_trans = curv_fact_elec_trans_;
  curv_fact_colim = curv_fact_colim_;
  TF24_cost_scale = TF24_cost_scale_;
  R_d_25 = R_d_25_;
  JS22_gamma = JS22_gamma_;
  CMax_a = CMax_a_;
  CMax_b = CMax_b_;

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
  // which is all of it bar `CF77_lambda_` and `TF24_floor_lambda_o` -- caller inputs,
  // left standing on purpose (#96, see their declarations).
  setup_clean_leaf();
}



// ⚠️ THE ROWS ARE IN (P50, c), NOT (b, c). stem_b is derived, so a row in it
// describes a parameter this model does not have -- and set_traits re-derives b
// from P50 on every call, so a seed in b reaches nothing at all. The chain is
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
  // The value is the stem table's; the rows are the closed form's. One graft
  // serves both curves -- see graft.hpp for why that is not a convenience.
  return graft_integral<S>(
      stem_curve_integral(odelia::util::to_passive(psi),
                          "Leaf::stem_integral_at"),
      psi, pars[par_stem_P50], pars[par_stem_c]);
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
  p[par_a] = a;
  p[par_curv_fact_elec_trans] = curv_fact_elec_trans;
  p[par_curv_fact_colim] = curv_fact_colim;
  p[par_TF24_cost_scale] = TF24_cost_scale;
  p[par_R_d_25] = R_d_25;
  p[par_JS22_gamma] = JS22_gamma;
  p[par_CMax_a] = CMax_a;
  p[par_CMax_b] = CMax_b;
  p[par_kmax] = leaf_specific_conductance_max_;
  // The single-potential path's series resistance. Empty on the layered route,
  // where the resistances are per layer and this slot has no scalar to hold.
  p[par_resistance] = roots_.network_.r_R_V_sum.empty() ? util::na_value
                                    : roots_.network_.r_R_V_sum.front();
  p[par_CF77_lambda] = CF77_lambda_;
  p[par_TF24_floor_lambda_o] = TF24_floor_lambda_o;
  return p;
}

// set various states and physiology parameters obtained from TF24 to NA to clean leaf object
inline void Leaf::setup_clean_leaf() {
  ci_at_compensation_point_ = false;  // hazard 8: every exit writes its own state
  ci_ = util::na_value; // Pa
  stom_cond_CO2_= util::na_value; //mol Co2 m^-2 s^-1 
  assim_colimited_= util::na_value; // umol C m^-2 s^-1 
  transpiration_= util::na_value; // kg m^-2 s^-1 
  profit_= util::na_value; // umol C m^-2 s^-1
  // CF77_lambda_ and TF24_floor_lambda_o are deliberately NOT here: they are the
  // caller's inputs, not derived state, and carry their NA defaults at the
  // declarations instead (#96). Adding either back makes a prescribed price
  // survive set_drivers() and vanish on set_traits().
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
  Tleaf_= util::na_value; // deg C -- an output, so NA until a solve writes it
  Tair_= util::na_value; // deg C
  Rn_= util::na_value; // W m^-2
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
// NOTE: the temperature-dependent block (group 1) depends only on leaf_temp_
// (constant across the run in the current driver setup) yet is recomputed on
// every call; see the optimisation notes / caching opportunity.
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
   // leaf_temp driver as air temperature; derive net radiation from the absorbed
   // PAR (PPFD_, umol m^-2 s^-1): shortwave ~= 2*PAR converted to W m^-2, plus a
   // fixed clear-sky longwave cooling offset (doc 3.3). ra fixed for the minimal
   // cut (doc fallback). Only read on the use_energy_balance_ path.
   Tair_ = leaf_temp_;
   Rn_ = sw_abs_per_par * PPFD_ / umol_par_per_joule + longwave_net_offset;
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
   // the Tair baseline (used for assim_max_ / feasibility) and let the solve
   // override it.
   const std::array<double, photo_temp_key_size> photo_key = photo_temp_key();
   if (!use_energy_balance_ &&
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
// psi_root, transpiration's psi_upstream, and all four spline domains are the
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
// residual over the collar suction x in [wettest_soil_layer, psi_crit]) are
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
// psi_crit endpoint is exactly the E_column(psi_crit) < 0 shutdown test, and
// crit=0 is only reached for soil wetter than psi_crit. The genuinely
// non-smooth failure mode in the earlier blanket-swap rejection (the root
// vulnerability spline extrapolating negative beyond its ~root_psi_crit domain)
// lives in E_from_Soil_to_Root_Collar itself and bites both solvers identically;
// it does not arise on the brackets this finder is actually handed. Like
// psi_stem_to_ci (Phase 6) this is a same-tolerance method swap, NOT a tolerance
// loosening: same root, fewer evals.
inline double Leaf::find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit) {
  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  if (find_root_crit == 1) {
    auto target = [&](double x) -> double {
      return E_column(x, psi_soil, psi_crit);
    };
    try {
      // Ends swapped, not just renamed: a bracketing solver needs opposite signs
      // at the two endpoints, but the LOWER bound must come first, and in
      // magnitudes the wettest layer is the smallest suction (#25).
      return util::uniroot_smooth(target, wettest_soil_layer, psi_crit, 1e-4, ci_niter);
    } catch (const std::exception& e) {
      util::stop_infeasible("collar_bracket",
                 "find_root_psi(find_root_crit=1) failed: " + std::string(e.what()) +
                 "; min=" + util::to_string(wettest_soil_layer) +
                 "; max=" + util::to_string(psi_crit));
    }
  }

  auto target = [&](double x) -> double {
    return E_column_zero(x, psi_soil);
  };
  try {
    return util::uniroot_smooth(target, wettest_soil_layer, psi_crit, 1e-4, ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
               "find_root_psi(find_root_crit=0) failed: " + std::string(e.what()) +
               "; min=" + util::to_string(wettest_soil_layer) +
               "; max=" + util::to_string(psi_crit));
  }

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
// This is the entry point called once per individual per environment update
// (from TF24_Strategy::net_mass_production_dt). It solves the whole
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
//   2. OPTIMISATION (Golden-Section Search): among feasible collar potentials,
//      choose the one that maximises carbon profit = assimilation - hydraulic
//      cost (profit_psi_stem_TF). The collar potential is bracketed between
//      `root_zero_E` (collar where soil uptake is zero, the wettest feasible
//      point) and `root_crit` (collar at which the stem reaches psi_crit, the
//      driest feasible point), clamped to root_psi_crit.
//
// Several early-exit short-circuits avoid the (expensive) GSS loop when no
// meaningful optimisation is possible. In each case the plant is effectively
// shut down (operating at psi_crit, paying only respiration + hydraulic cost):
//   * wettest soil layer is already drier than psi_crit -> no transpiration;
//   * even at psi_crit the soil cannot supply the demanded flux (E_column<0);
//   * the continuity root would require the collar drier than psi_crit;
//   * maximum possible assimilation (at ci = ca) is negative.
//
// Implementation note: psi_soil is a positive magnitude here as everywhere else
// (#25), so nothing is flipped. The GSS reuses one
// profit evaluation per iteration (golden ratio) to halve function calls, and
// a collapsed-interval branch handles the degenerate single-feasible-point case.
// Shut-down operating point shared by find_root_collar_psi's early-exits: the
// stem is held at psi_crit (transpiration not possible), so the plant pays only
// respiration (R_d_) plus the hydraulic cost at psi_crit. Only the recorded
// root-collar potential differs between the calling cases.
template <Leaf::CostCurve K>
inline void Leaf::set_shutdown_state(double root_collar) {
  opt_root_psi_ = root_collar;
  opt_psi_stem_ = psi_crit;
  // FIRST, because every number below is read off the temperature block: on the
  // PM path this exit's own leaf temperature is the E = 0 one, and the block
  // arrives holding the Tair baseline set_physiology derived. Forming profit
  // before this line reported a respiration rate belonging to a different
  // temperature than the Tleaf beside it -- 35.8% low at Tair 30, in the
  // direction that makes shutting down look cheaper than it is (#105).
  //
  // Zero transpiration is where the gap is LARGEST, not smallest: no latent
  // cooling, so this is the hottest the leaf gets.
  Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(0.0) : leaf_temp_;
  // The leaf-to-air deficit belongs to the same temperature, for the same reason
  // and by the same argument #93 makes at the zero-transpiration branch of
  // set_leaf_states_rates_from_psi_stem: `vpd_leaf_` is reported, and `g1_eff()`
  // reads it. Left alone it held the Tair deficit set_physiology seated. A no-op
  // off the PM path -- `set_leaf_vpd` returns `atm_vpd_` there exactly.
  set_leaf_vpd(Tleaf_);
  if (use_energy_balance_) {
    // Leaves the block at Tleaf rather than Tair. Safe because set_physiology
    // bypasses the photo_temp cache entirely when the gate is on, so the next
    // solve re-derives the baseline instead of taking a hit on stale members.
    update_temperature_dependent_params(Tleaf_);
  }
  profit_ = no_flow_profit<K>(psi_crit);
  // Tagged here rather than at the four call sites, so a fifth reason to shut
  // down cannot arrive without a classification. All four are the same
  // ecological statement -- the soil is too dry for any collar potential that
  // both moves water and stays inside the critical potentials.
  operating_point_kind_ = OperatingPointKind::HydraulicShutdown;

  // Write the flux outputs too. Before this they were left holding whatever the
  // PREVIOUS solve on this object put there, and set_physiology did not reset them
  // either (setup_clean_leaf runs from the constructors only) -- so a reused Leaf
  // reported the previous plant's water and carbon use after a shutdown. plant
  // holds one persistent Leaf per TF24_Strategy and drives every node, height and
  // timestep through it, and soil_consumption_ feeds the patch water balance, so
  // this was a live water-balance error on the dry margin. plant #578, #577.
  //
  // The values are the ones consistent with the profit already set above: the leaf
  // is holding at psi_crit, so it moves no water at all, but it IS respiring --
  // profit_ is -R_d_ - hydraulic_cost, so assimilation is -R_d_, not zero. ci sits
  // at the CO2 compensation point, matching the two zero-transpiration branches in
  // set_leaf_states_rates_from_psi_stem.
  transpiration_ = 0.0;
  stom_cond_CO2_ = 0.0;
  assim_colimited_ = -R_d_;
  ci_ = gamma_ * umol_per_mol_to_Pa_;
  E_up_ = 0.0;
  // `Tleaf_` is written at the top of this function, not here, because the
  // temperature block is derived from it and `R_d_`/`gamma_` above are read out
  // of that block. The two used to be set at opposite ends and disagree.
  std::fill(soil_consumption_.begin(), soil_consumption_.end(), 0.0);
  // ⚠️ AND THE EMERGENT LAMBDA, for exactly the reason the block above exists. A
  // leaf that moves no water has no marginal cost of moving it, so the honest
  // value is NA -- and without this line it kept whatever the PREVIOUS solve on
  // this object left, which on a reused Leaf is a finite, plausible number
  // belonging to a different plant. `leaf_solve()` reuses one Leaf across rows,
  // so a shut-down row reported the previous row's price: 3.49e+05 against the NA
  // the same drivers give on a fresh leaf. Found when `lambda_emergent` was added
  // to the reported vector; before that the field was reachable only off an object
  // the caller had driven themselves, so nobody met it.
  //
  // `shadow_cost()` needs no equivalent -- it is computed from `transpiration_`,
  // which is set to zero just above, so it follows to zero on its own.
  lambda_emergent_ = util::na_value;
  // Invalidate the transpiration memo: it is keyed on (psi_stem, psi_upstream) and
  // we have just written transpiration_ without going through transpiration().
  transpiration_cached_ = false;
}

// Shared setup + feasibility handling for the root-collar solve. Extracted
// verbatim from find_root_collar_psi (no reordering of floating-point ops, so
// the TF24 optimisation path stays bit-identical) and reused by
// evaluate_root_collar_psi. Returns false when the final operating point is
// already determined here (shutdown / assim<0 / collapsed interval) and the
// caller should stop; returns true with [bound_a, bound_b] set to the feasible
// collar-potential interval (positive magnitudes) otherwise.
template <Leaf::CostCurve K>
inline bool Leaf::prepare_collar_solve(double& bound_a, double& bound_b){

  // Clear the classification FIRST, so that a path which declines to write it
  // reports "unclassified" rather than the previous plant's kind of operating
  // point (hazard 8). Every exit below, and both callers, write it again.
  operating_point_kind_ = OperatingPointKind::Unsolved;

  // Hand the supply path the start of a solve: it builds its per-solve caches and
  // reports back the wettest rooted layer -- the SMALLEST suction, and the lower
  // bracket endpoint everything below needs.
  const double wettest_soil_layer = supply_begin_solve();

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to
  // shut down

  if (wettest_soil_layer >= psi_crit){
    set_shutdown_state<K>(psi_crit);
    return false;
  }

if(E_column(psi_crit, supply_psi_soil(), psi_crit) < 0){
      set_shutdown_state<K>(supply_psi_crit());
      return false;
}

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to
  // shut down
double root_crit = find_root_psi(wettest_soil_layer, supply_psi_soil(), 1);

// If root crit would have to be larger than psi crit, also avoid loop as above

    if (root_crit >= psi_crit){
    set_shutdown_state<K>(root_crit);
    return false;
  }

// Find root collar where transpiration from soil is 0
double root_zero_E = find_root_psi(wettest_soil_layer, supply_psi_soil(), 0);

// If assimilation would be less than 0 even at Ca, also end loop
if(assim_max_ < 0){
    // At zero transpiration the stem equilibrates with the collar (no flux, no
    // gradient), so the operating point is root_zero_E for both -- and now they
    // are literally the same number, where #7 had to pair a magnitude with its
    // negation to say the same thing.
    opt_psi_stem_ = root_zero_E;
    opt_root_psi_ = root_zero_E;
    E_from_Soil_to_Root_Collar(opt_root_psi_, supply_psi_soil());

    // As at the shut-down exits, and for the same reason: the reported state is
    // the E = 0 one, so the temperature block has to be at the E = 0 temperature
    // before `R_d_` and `gamma_` are read out of it (#105).
    //
    // ⚠️ THE TEST ABOVE STAYS AT THE Tair BASELINE, and that is a decision rather
    // than an oversight. `assim_max_` asks "can gross assimilation at ci = ca
    // cover respiration", i.e. whether any OPEN operating point is worth taking --
    // and an open leaf transpires, so it is cooler than this exit's leaf. Re-taking
    // the test at the zero-transpiration temperature would judge the coolest option
    // at the hottest temperature and shut down leaves that would have paid. The
    // asymmetry is therefore: the branch is chosen on the transpiring baseline, the
    // state it reports is self-consistent at the temperature it describes.
    Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(0.0) : leaf_temp_;
    set_leaf_vpd(Tleaf_);
    if (use_energy_balance_) {
      update_temperature_dependent_params(Tleaf_);
    }

    profit_ = no_flow_profit<K>(opt_root_psi_);
    // As on the shut-down exits: transpiration is zero here, so gross
    // assimilation is zero and the reported net rate is -R_d_. Set it
    // explicitly -- this branch does not go through profit_psi_stem_TF, so
    // assim_colimited_ would otherwise keep whatever the last probe wrote,
    // and it is reported. Keeps profit_ == assim_colimited_ -
    // hydraulic_cost_TF() in every branch.
    assim_colimited_ = -R_d_;
    // Shutdown on LIGHT, not on water: gross assimilation at ci = ca cannot
    // cover dark respiration. Tagged separately from the hydraulic exits
    // because the two respond to different drivers -- a rainfall sequence
    // traverses the hydraulic ones and never reaches this.
    operating_point_kind_ = OperatingPointKind::ShadeDeath;
    // E_up_ and soil_consumption_ are already correct: the
    // E_from_Soil_to_Root_Collar call above evaluates them at root_zero_E, the
    // collar potential at which uptake is zero. The leaf-side pair is set
    // nowhere on this path, though, so zero it here rather than leave the
    // previous solve's values -- see set_shutdown_state for why that matters.
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;
    // And `ci_`, which this exit never wrote at all -- so it reported the last
    // solve's internal CO2 beside a state that moves no gas, or the NA sentinel
    // on a cold object. The compensation point, matching `set_shutdown_state` and
    // both zero-transpiration branches of set_leaf_states_rates_from_psi_stem.
    ci_ = gamma_ * umol_per_mol_to_Pa_;
    // And the emergent lambda, for the same reason and by the same argument as in
    // set_shutdown_state: no water moves here, so there is no marginal cost of
    // moving it, and leaving it unwritten hands a reused Leaf the previous row's
    // price. This exit is shade death rather than hydraulic shutdown, but the
    // stale-state hazard is identical.
    lambda_emergent_ = util::na_value;

        if(std::isnan(profit_)){
          util::stop_infeasible("collar_solve", "profit is not finite at the shade-death "
                                                "exit");
    }

    return false;
}
// opt_psi_stem_ = psi_soil_;


  // optimise for stem water potential
    bound_a = root_zero_E;
    // The dry end of the feasible interval is whichever limit binds FIRST: the
    // continuity root, or the potential at which root conductivity is down to 5%.
    // Both are positive magnitudes, so that is a min (#24, plant #584).
    //
    // This line used to read std::max(-root_crit, -root_psi_crit) -- a magnitude
    // compared against a signed potential, so the second term was always negative
    // and the clamp could never bind. The master solver's comment has claimed the
    // bracket is "clamped to root_psi_crit" throughout; it now is. In magnitudes
    // the correct form is the obvious one and the trap is gone, which is the
    // clearest argument for #25 there is: the bug was a property of having two
    // representations, not of this line.
    bound_b = std::min(root_crit, supply_psi_crit());

    // ⚠️ The clamp can INVERT the interval, and nothing handled that before,
    // because with the clamp dead it could not happen. bound_a is root_zero_E, the
    // collar suction at which uptake is exactly zero: to draw any water at all the
    // collar must pull harder than that. So when root_psi_crit lands *below*
    // root_zero_E there is no operating point that both moves water and stays
    // inside the root vulnerability limit, and the answer is shut-down rather than
    // an optimisation over an empty interval.
    //
    // Measured at psi_crit = 5.91988 (plant's TF24 value, drier than
    // root_psi_crit = 5.870283), single layer:
    //
    //   psi_soil   root_zero_E   root_crit   bound_b     inverted
    //     5.850      5.854900     5.863944   5.863944    no
    //     5.880      5.884900     5.889736   5.870283    YES
    //     5.910      5.914900     5.915587   5.870283    YES
    //
    // Without this exit, golden_section_max is handed bound_a > bound_b and returns
    // a point between them -- past the limit the clamp exists to enforce, which is
    // the very failure #24 is about, reintroduced by its own fix.
    if (bound_b < bound_a) {
      set_shutdown_state<K>(supply_psi_crit());
      return false;
    }

    // If no interval exists (single feasible root-collar value), use that
    // point directly as the alternative solution instead of running GSS.
    if (std::abs(bound_b - bound_a) <= GSS_tol_abs) {
      const double opt_root_psi = 0.5 * (bound_a + bound_b);
      const double psi_stem_single = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());

      if (!std::isfinite(psi_stem_single)) {
        util::stop_infeasible("collar_solve",
                   "non-finite psi_stem_single in collapsed-root interval; "
                   "opt_root_psi=" + util::to_string(opt_root_psi) +
                   "; bound_a=" + util::to_string(bound_a) +
                   "; bound_b=" + util::to_string(bound_b) +
                   "; root_crit=" + util::to_string(root_crit) +
                   "; root_zero_E=" + util::to_string(root_zero_E) +
                   "; E_up_=" + util::to_string(E_up_));
      }

      // ⚠️ THE COLLAR IS DETERMINED; psi_stem IS NOT. This branch used to DERIVE
      // psi_stem from continuity at the collapsed collar and stop, which threw the
      // remaining degree of freedom away: the leaf can still choose how far to let
      // psi_stem fall, it is only the COLLAR that feasibility has pinned. Measured
      // at psi_soil 1.0, PPFD 900, single potential: deriving gave psi_stem 1.6465
      // and profit 8.255 where optimising gives 3.0935 and 13.389 -- a 38% loss,
      // and the reason the collar route's r -> 0 limit did not converge onto the
      // stem route's answer.
      //
      // So optimise psi_stem with the upstream potential pinned at the determined
      // collar. That is exactly the problem `optimise_psi_stem_single` solves, and
      // the same shared solver does it -- `profit_psi_stem_for` and
      // `dprofit_dpsi_stem` both already take the upstream potential as an
      // argument, so nothing new is derived here.
      //
      // ⚠️ COSTS NOTHING ON THE NORMAL PATH. This runs only where the interval has
      // collapsed, which is where the solve previously gave up; a bracket wider
      // than GSS_tol_abs never reaches it.
      double collapsed_profit = 0.0;
      opt_psi_stem_ = util::maximise_over_closed_interval_foc(
          [&](double psi_stem) {
            return profit_psi_stem_for<K>(psi_stem, opt_root_psi);
          },
          [&](double psi_stem, bool* ok) -> double {
            if constexpr (benefit_link<K>() != BenefitLink::Identity) {
              if (use_energy_balance_) {
                if (ok != nullptr) *ok = false;
                return 0.0;
              }
            }
            return dprofit_dpsi_stem<K>(psi_stem, opt_root_psi, ok);
          },
          opt_root_psi, psi_crit, basin_scan_cells(), collar_root_tol,
          static_cast<size_t>(ci_niter), &collapsed_profit);
      profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi);
      (void)psi_stem_single;
      opt_root_psi_ = opt_root_psi;
      // Feasibility DETERMINED this point; no maximisation happened, and there is
      // no free variable left for a derivative to move.
      operating_point_kind_ = OperatingPointKind::Determined;

      if (!std::isfinite(profit_)) {
        util::stop_infeasible("collar_solve",
                   "non-finite profit in collapsed-root interval; "
                   "opt_psi_stem_=" + util::to_string(opt_psi_stem_) +
                   "; opt_root_psi_=" + util::to_string(opt_root_psi_) +
                   "; bound_a=" + util::to_string(bound_a) +
                   "; bound_b=" + util::to_string(bound_b) +
                   "; root_crit=" + util::to_string(root_crit) +
                   "; root_zero_E=" + util::to_string(root_zero_E) +
                   "; E_up_=" + util::to_string(E_up_) +
                   "; assim_colimited_=" + util::to_string(assim_colimited_) +
                   "; hydraulic_cost_=" + util::to_string(hydraulic_cost_));
      }
      return false;
    }

    return true;
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
  const double width = bound_b - bound_a;
  const auto dprofit = [&](double psi, bool* feasible = nullptr) {
    return dprofit_at_collar_psi<K>(psi, feasible);
  };

  // The WET endpoint is infeasible by construction: bound_a is root_zero_E, where
  // uptake is exactly zero, so psi_stem == bound_a and dprofit takes its
  // reversed-gradient exit. Its 0.0 is a sentinel, and a bracketing solver handed
  // it would report the zero-transpiration point as the optimum -- profit -1.897
  // there against 2.516 at the true optimum, measured at the default drivers.
  //
  // So step inside until the gradient is real. Measured over the golden grid the
  // infeasible region is at most 3.46e-07 MPa wide (median 1.22e-08, never more
  // than 6.2e-07 of the bracket), so the first try lands and the loop is a guard.
  bool ok_lo = false;
  double lo = bound_a;
  double f_lo = 0.0;
  for (double frac = 1e-6; frac < 0.5; frac *= 10.0) {
    lo = bound_a + frac * width;
    f_lo = dprofit(lo, &ok_lo);
    if (ok_lo && std::isfinite(f_lo)) {
      break;
    }
  }

  // The DRY endpoint is feasible everywhere measured, so it is tried as-is first
  // and only stepped inside if that fails.
  bool ok_hi = false;
  double hi = bound_b;
  double f_hi = dprofit(hi, &ok_hi);
  if (!(ok_hi && std::isfinite(f_hi))) {
    for (double frac = 1e-6; frac < 0.5; frac *= 10.0) {
      hi = bound_b - frac * width;
      f_hi = dprofit(hi, &ok_hi);
      if (ok_hi && std::isfinite(f_hi)) {
        break;
      }
    }
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
    operating_point_kind_ = OperatingPointKind::BoundaryCrit;
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
    if (!prepare_collar_solve<K>(bound_a, bound_b)) {
      return;
    }
    // root_crit / root_zero_E were consumed inside prepare_collar_solve; recover
    // them for the diagnostic message only if the profit check below fails.

    // Maximise carbon profit over the feasible collar-potential interval by
    // solving its first-order condition, dprofit/dpsi == 0, rather than by
    // searching the objective.
    const double opt_root_psi = maximise_profit_over_collar<K>(bound_a, bound_b);

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());

    opt_root_psi_ = opt_root_psi;
    // ⚠️ `<K>`, NOT TF24. This read TF24's objective for every curve after the
    // collar route was templated -- so a SOX solve maximised SOX and then reported
    // TF24's profit at SOX's argmax, which is finite, plausible, and always a
    // little below TF24's own maximum. Two of these survived the templating pass.
    profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi);
    // Reported on the same axis as the single-layer solvers, using the MULTI-LAYER
    // lambda: the collar is the free variable here, so the marginal cost carries
    // the series-resistance correction for the soil-to-collar path that the
    // single-layer form has no term for.
    lambda_emergent_ = marginal_cost_water_multilayer();

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
// optimising it (see header). Clamps the target into the feasible interval so a
// tracked state that has drifted outside it still yields a finite operating
// point; the gradient computed by the caller then pulls it back inside.
// The objective's value where no water moves. The TF24 arm is spelled out because
// these exits also hand-write their own rate fields (hazard 8) and the golden files
// pin the result; every other curve goes through the shared evaluator, which is what
// gets a product objective's `A*g` right rather than an `A - C` that does not apply.
template <Leaf::CostCurve K>
inline double Leaf::no_flow_profit(double psi_stem) {
  if constexpr (K == CostCurve::TF24) {
    return -R_d_ - hydraulic_cost_TF(psi_stem);
  } else {
    return profit_psi_stem_for<K>(psi_stem, psi_stem);
  }
}


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


inline double Leaf::evaluate_root_collar_psi(double target_opt_root_psi) {
  return evaluate_root_collar_psi_for<CostCurve::TF24>(target_opt_root_psi);
}


inline double Leaf::dprofit_droot_collar_psi(double opt_root_psi, bool* feasible) {
  return dprofit_droot_collar_psi_for<CostCurve::TF24>(opt_root_psi, feasible);
}


template <Leaf::CostCurve K>
inline double Leaf::evaluate_root_collar_psi_for(double target_opt_root_psi){
    double bound_a, bound_b;
    if (!prepare_collar_solve<K>(bound_a, bound_b)) {
      // Operating point fully determined by feasibility handling (shutdown /
      // assim<0 / collapsed interval); profit_ is already set.
      return profit_;
    }

    return profit_at_collar_psi<K>(target_opt_root_psi, bound_a, bound_b);
}

// Post-prepare body of evaluate_root_collar_psi (see header). Kept as a separate
// entry point so callers that evaluate several collar potentials within one step
// (the centred finite difference, #530) can run prepare_collar_solve once and
// reuse the soil-side caches across every profit eval. The clamp into
// [bound_a, bound_b] is identical to evaluate_root_collar_psi's, so near a
// boundary a perturbed potential collapses onto the boundary -- which is exactly
// how the FD path degrades gracefully to a one-sided difference.
template <Leaf::CostCurve K>
inline double Leaf::profit_at_collar_psi(double target_opt_root_psi,
                                  double bound_a, double bound_b){
    const double opt_root_psi =
        std::min(std::max(target_opt_root_psi, bound_a), bound_b);

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());
    opt_root_psi_ = opt_root_psi;
    profit_ = profit_psi_stem_for<K>(opt_psi_stem_, opt_root_psi);
    // Well defined off the optimum too: dC/dE is a property of the point, not of
    // its being stationary. At a prescribed point it simply is not equal to dA/dE.
    lambda_emergent_ = marginal_cost_water_multilayer();
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
// (smooth spline transport) and ci = psi_stem_to_ci(psi_stem, psi) (root-find).
// Chain rule:
//   dprofit/dpsi = A'(ci) dci/dpsi - C'(psi_stem) dpsi_stem/dpsi
// where dci/dpsi = (dci/dpsi_stem) dpsi_stem/dpsi + (dci/dpsi)|_explicit, and the
// dci/d* terms come from the implicit-function theorem on the residual
//   g(ci; psi_stem, psi) = A(ci) umol_to_mol - gc(psi_stem,psi) (ca-ci)/(atm kPa)
// with gc = const * transpiration(psi_stem,psi). A'/C' are obtained by forward
// AD; the gc partials use the analytic spline derivative (transpiration_from_psi
// .deriv); dpsi_stem/dpsi by a tight central difference on the smooth transport.
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

  // Operating point in double.
  const double psi_stem = find_psi_stem_from_psi_root(psi, supply_psi_soil());
  // Shut down before the ci solve, not after. psi and psi_stem are both positive
  // magnitudes here, so psi >= psi_stem is the no-flow / reversed-gradient case --
  // the same condition set_leaf_states_rates_from_psi_stem treats as zero
  // transpiration. It has to be caught *here* because psi_stem_to_ci does not
  // return non-finite in that state, it throws: gc = const * transpiration goes
  // negative, which flips the sign of the supply term so the residual no longer
  // crosses zero over (gamma*, ca] and the bracketing root-find reports that its
  // endpoints do not bracket a root. The isfinite check below was written to cover
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
  double dT_dE = 0.0;
  double Tleaf_here = leaf_temp_;
  if (use_energy_balance_) {
    Tleaf_here = leaf_temp_from_E(transpiration(psi_stem, psi), &dT_dE);
    update_temperature_dependent_params(Tleaf_here);
    set_leaf_vpd(Tleaf_here);
  }

  const double ci = psi_stem_to_ci(psi_stem, psi);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  // dpsi_stem/dpsi. Computed HERE, before the compensation-point branch, so both
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
    dpsistem_dpsi = stem_curve_integral_inverse_deriv(E_psi_stem) * dEpsistem_dpsi;
  } else {
    // Near a branch kink the analytic conductance returns NaN; fall back to a
    // central difference on the transport, as this path has always done.
    const double h = 1e-6;
    dpsistem_dpsi =
        (find_psi_stem_from_psi_root(psi + h, supply_psi_soil()) -
         find_psi_stem_from_psi_root(psi - h, supply_psi_soil())) / (2.0 * h);
  }


  // ⚠️ E3: THE COMPENSATION-POINT BRANCH, where the implicit function theorem
  // below is void. When the leaf is hot enough that assimilation is negative
  // across the whole [gamma*, ca] bracket, psi_stem_to_ci cannot find a
  // supply==demand root and places ci AT gamma* instead. The residual g is then
  // not zero, so differentiating "g = 0" is meaningless -- and worse, it is
  // silently meaningless: at the wet bracket endpoint transpiration is ~0 so
  // gc ~ 0, and at the compensation point A'(ci) ~ 0, which makes
  // g_ci = A'*umol_to_mol + gc*inv_atm ~ 0 and the IFT quotient a 0/0. That NaN
  // then propagated into maximise_profit_over_collar as f_lo = f_hi = NaN, where
  // both sign tests are false for NaN, and TOMS748 aborted the whole run with
  // "parameters a and b do not bracket the root". Found exactly that way.
  //
  // On this branch gross assimilation is identically zero, so A = -R_d(T) and
  // the only surviving temperature dependence is respiration's:
  // `dprofit/dpsi = -R_d'(T) * tau - C'(psi_stem) * dpsi_stem/dpsi`.
  //
  // R_d' is obtained the same way A_T is, by differencing the model's own
  // temperature block, so the two cannot drift apart.
  if (ci_at_compensation_point_) {
    const double C_prime0 = cost_deriv<K>(psi_stem, psi);
    double dprofit = -C_prime0 * dpsistem_dpsi;
    if (use_energy_balance_ && dT_dE != 0.0) {
      const double Rd_T = respiration_temp_deriv(Tleaf_here);
      // dE/dpsi from the same spline derivatives the main branch uses;
      // recomputed here because the main branch's locals are below this early
      // return. It used to be written as the two gc partials divided back by
      // their shared constant, which cancelled a coefficient that is no longer a
      // single number now that the deficit moves with Tleaf. Stated directly.
      const double dE_dpsi =
          leaf_specific_conductance_max_ *
          (stem_curve_integral_deriv(psi_stem) * dpsistem_dpsi -
           stem_curve_integral_deriv(psi));
      dprofit += -Rd_T * dT_dE * dE_dpsi;
    }
    return std::isfinite(dprofit) ? dprofit : 0.0;
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
  // gc_const * transpiration(psi_stem, psi); transpiration is conductance_max *
  // (transp_from_psi(psi_stem) - transp_from_psi(psi)), so the partials use the
  // analytic spline derivative.
  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double gc = gc_const * transpiration(psi_stem, psi);
  const double dgc_dpsistem =
      gc_const * leaf_specific_conductance_max_ * stem_curve_integral_deriv(psi_stem);
  const double dgc_dpsi =
      gc_const * leaf_specific_conductance_max_ * (-stem_curve_integral_deriv(psi));

  // IFT on g(ci; psi_stem, psi): dci/dp = -(dg/dp)/(dg/dci).
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;      // dg/dci
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;
  const double dci_dpsi_expl = -(-dgc_dpsi * (ca_ - ci) * inv_atm) / g_ci;

  // dpsi_stem/dpsi: psi_stem = P(E_psi_stem) with
  //   E_psi_stem = E_up_(psi)/k_max + S(psi),
  // S = transpiration_from_psi, P = psi_from_transpiration (both C2 splines), and
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
  }
  // Gated at the CALL SITE, not just inside the callee: the block is out of line
  // (deliberately, so adding it cannot change FMA contraction in this inlined
  // body), and an out-of-line call costs even when it returns 0.0 immediately.
  if (!use_energy_balance_) {
    return base;
  }
  return base + dprofit_energy_balance_term(ci, gc, g_ci, inv_atm, gc_const,
                                            A_prime, dgc_dpsistem, dgc_dpsi,
                                            dpsistem_dpsi, dT_dE, Tleaf_here);
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
// `cost_deriv<K>`, `benefit_link_deriv<K>`, `respiration_temp_deriv` and
// `dprofit_energy_balance_term`. What is left duplicated is the arithmetic the
// golden files pin.
// dC/dpsi_stem for one cost curve: THE ONLY PLACE THE COST CURVE ENTERS a
// derivative, and shared by both routes. `dprofit_dpsi_stem` uses it directly;
// `dprofit_at_collar_psi` multiplies it by dpsi_stem/dpsi_collar. Factored out so
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
    // no AD and no spline -- the simplest derivative of the four.
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
  } else if constexpr (K == CostCurve::ProfitMax) {
    // C = [k(psi_s) - k(psi)]/k_span, so dC/dpsi = kmax*|f'|/k_span. The thermal
    // term is constant in psi with the energy balance off, and the guard below
    // refuses the case where it is not.
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

  double dT_dE = 0.0;
  double Tleaf_here = leaf_temp_;
  if (use_energy_balance_) {
    Tleaf_here = leaf_temp_from_E(transpiration(psi_stem, psi_upstream), &dT_dE);
    update_temperature_dependent_params(Tleaf_here);
    set_leaf_vpd(Tleaf_here);
  }

  const double ci = psi_stem_to_ci(psi_stem, psi_upstream);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  const double C_prime = cost_deriv<K>(psi_stem, psi_upstream);

  // ci pinned at the compensation point: the supply==demand residual is not zero
  // there, so the implicit function theorem below does not hold and the
  // assimilation term is dropped rather than approximated. `feasible` stays
  // false, which is what tells a caller the difference.
  if (ci_at_compensation_point_) {
    double dprofit = -C_prime;
    if (use_energy_balance_ && dT_dE != 0.0) {
      const double Rd_T = respiration_temp_deriv(Tleaf_here);
      // dpsi_stem/dpsi is 1 here and the upstream potential is fixed, so the
      // collar version's two-term dE/dpsi collapses to this one.
      const double dE_dpsi =
          leaf_specific_conductance_max_ * stem_curve_integral_deriv(psi_stem);
      dprofit += -Rd_T * dT_dE * dE_dpsi;
    }
    return std::isfinite(dprofit) ? dprofit : 0.0;
  }

  if (feasible != nullptr) {
    *feasible = true;
  }

  AD ci_ad = ci;
  xad::derivative(ci_ad) = 1.0;
  const double A_prime = xad::derivative(assim_colimited_kernel(ci_ad));

  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double gc = gc_const * transpiration(psi_stem, psi_upstream);
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
  }
  if (!use_energy_balance_) {
    return base;
  }
  if constexpr (benefit_link<K>() != BenefitLink::Identity) {
    // ⚠️ REFUSED RATHER THAN APPROXIMATED. The energy-balance term below is a
    // further contribution to dA/dpsi through Tleaf(E), so under a non-identity
    // link it needs the same `h'` factor -- and it also carries cost-side pieces
    // (ProfitMax's thermal term moves with Tleaf) that are not separated in it.
    // Returning it unscaled would be plausible and wrong.
    util::stop("dprofit/dpsi_stem with a non-identity benefit link and the "
               "energy balance ON is not implemented: the temperature term needs "
               "the same h'(A) factor and, for ProfitMax, its thermal cost term "
               "separated out. Solve with the optimiser, or run with the energy "
               "balance off.");
  }
  // dgc_dpsi = 0 and dpsistem_dpsi = 1: the term derives dE/dpsi from those two,
  // and with the upstream potential fixed that reduces to kmax * f(psi_stem).
  return base + dprofit_energy_balance_term(ci, gc, g_ci, inv_atm, gc_const,
                                            A_prime, dgc_dpsistem, 0.0, 1.0,
                                            dT_dE, Tleaf_here);
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

  // Drier soil than the stem can reach is the no-flow case, and the target is
  // irrelevant to it: there is one feasible potential and this is it.
  if (psi_soil > psi_crit) {
    opt_psi_stem_ = psi_soil;
    profit_ = profit_psi_stem_for<K>(psi_soil, psi_soil);
    operating_point_kind_ = OperatingPointKind::Prescribed;
    return profit_;
  }

  // Clamp, rather than refuse: a caller carrying psi_stem as a tracked state can
  // arrive here just outside the interval, and a finite operating point plus a
  // gradient that points back inside is more useful than an exception. The tag
  // does not distinguish a clamped target from an interior one -- compare the
  // returned opt_psi_stem_ with what you asked for if you need to know.
  const double psi = std::min(std::max(target_psi_stem, psi_soil), psi_crit);

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
// term's dA/dT difference, and the eight-member restore is what makes calling this
// inside a derivative safe.
inline double Leaf::respiration_temp_deriv(double leaf_temp) {
  const double h = 1e-3;
  const double vc0 = vcmax_, jm0 = jmax_, ga0 = gamma_, ko0 = ko_,
               kc0 = kc_, rd0 = R_d_, km0 = km_, J0 = electron_transport_;
  update_temperature_dependent_params(leaf_temp + h);
  const double Rd_up = R_d_;
  update_temperature_dependent_params(leaf_temp - h);
  const double Rd_dn = R_d_;
  vcmax_ = vc0; jmax_ = jm0; gamma_ = ga0; ko_ = ko0; kc_ = kc0;
  R_d_ = rd0; km_ = km0; electron_transport_ = J0;
  return (Rd_up - Rd_dn) / (2.0 * h);
}

inline double Leaf::dprofit_energy_balance_term(
    double ci, double gc, double g_ci, double inv_atm, double gc_const,
    double A_prime, double dgc_dpsistem, double dgc_dpsi, double dpsistem_dpsi,
    double dT_dE, double Tleaf) {
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
  // model whose own floor (psi_stem_to_ci at 1e-10) is ~1e-9.
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
  // ⚠️ Tleaf, NOT leaf_temp_. On the energy-balance path set_physiology
  // reinterprets leaf_temp_ as AIR temperature (Tair_ = leaf_temp_), so
  // differencing around it evaluates dA/dT at the wrong point entirely. The
  // first version of this did exactly that; the resulting A_T was wrong enough
  // to flip the sign of dprofit at a bracket endpoint, and the collar solve
  // aborted with "parameters a and b do not bracket the root".
  const double T0 = Tleaf;

  update_temperature_dependent_params(T0 + h);
  const double A_up = assim_colimited(ci);
  update_temperature_dependent_params(T0 - h);
  const double A_dn = assim_colimited(ci);

  // Restore by assignment rather than by a third update call: exact, and cheaper
  // than re-running the Arrhenius block.
  vcmax_ = vcmax0; jmax_ = jmax0; gamma_ = gamma0; ko_ = ko0; kc_ = kc0;
  R_d_ = Rd0; km_ = km0; electron_transport_ = J0;

  const double A_T = (A_up - A_dn) / (2.0 * h);

  // dE/dpsi from the gc partials already in hand -- gc = gc_const * E, so
  // dividing by gc_const recovers it with no new spline evaluation.
  const double dE_dpsi = (dgc_dpsistem * dpsistem_dpsi + dgc_dpsi) / gc_const;
  const double tau = dT_dE * dE_dpsi;

  // dgc/dTleaf through the deficit alone. gc = K*E/D(T) with E hydraulically
  // pinned, so dgc/dT = -gc * D'(T)/D, and D'(T) = esat'(Tleaf) since e_air is
  // fixed. NEGATIVE: a hotter leaf has a larger deficit, so the same water flux
  // implies a SMALLER conductance. Zero where the deficit is on its floor, for
  // the same reason leaf_temp_from_E returns dT_dE = 0 on its clamp -- the
  // reported value no longer responds to the input.
  const double dgc_dT =
      (vpd_leaf_ > vpd_leaf_min)
          ? -gc * saturation_vapour_pressure_slope(Tleaf) / vpd_leaf_
          : 0.0;

  return tau *
         (A_T * gc * inv_atm + A_prime * dgc_dT * (ca_ - ci) * inv_atm) / g_ci;
}

template <typename T>
inline T Leaf::arrh_curve(T Ea, T ref_value, T leaf_temp) const {
  return ref_value*exp(Ea*((leaf_temp+C_to_K) - (25 + C_to_K))/((25 + C_to_K)*gas_constant*(leaf_temp+C_to_K)));
}

template <typename T>
inline T Leaf::peak_arrh_curve(T Ea, T ref_value, T leaf_temp, T H_d, T d_S) const {
  T arrh = arrh_curve<T>(Ea, ref_value, leaf_temp);
  T arg2 = 1 + exp((d_S*(25 + C_to_K) - H_d)/(gas_constant*(25 + C_to_K)));
  T arg3 = 1 + exp((d_S*(leaf_temp + C_to_K) - H_d)/(gas_constant*(leaf_temp + C_to_K)));

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
          double(use_thermal_cost_), T50_, Tcrit_};
}

inline void Leaf::update_temperature_dependent_params(double leaf_temp) {
  vcmax_ =
      peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_d_S_);
  jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_d_S_);
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
  km_ = (kc_*umol_per_mol_to_Pa_)*(1 + (atm_o2_kpa_*kPa_to_Pa)/(ko_*umol_per_mol_to_Pa_));
  electron_transport_ = electron_transport();
}

// Saturation vapour pressure es(T) in kPa (Tetens), T in deg C.
inline double Leaf::saturation_vapour_pressure(double temp) const {
  return 0.6108 * exp(17.27 * temp / (temp + 237.3));
}

// Slope of the saturation vapour pressure curve Delta(T) in kPa K^-1, T in deg C.
inline double Leaf::saturation_vapour_pressure_slope(double temp) const {
  return 4098.0 * saturation_vapour_pressure(temp) / ((temp + 237.3) * (temp + 237.3));
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
inline void Leaf::set_leaf_vpd(double leaf_temp) {
  if (!use_energy_balance_) {
    vpd_leaf_ = atm_vpd_;
    return;
  }
  const double d =
      atm_vpd_ + (saturation_vapour_pressure(leaf_temp) -
                  saturation_vapour_pressure(Tair_));
  vpd_leaf_ = std::max(d, vpd_leaf_min);
}

// Explicit leaf energy balance (#523): Tleaf = Tair + (Rn - lambda*E)*ra/(rho*cp).
// E is the hydraulically-pinned transpiration (kg H2O m^-2 s^-1), so lambda*E is
// the latent heat flux (W m^-2) and (Rn - lambda*E) the sensible heat flux H.
inline double Leaf::leaf_temp_from_E(double E, double* dT_dE) const {
  const double Tleaf = Tair_ + (Rn_ - latent_heat_vap * E) * ra_ / vol_heat_cap_air;
  // Clamp to a physical range so an extreme (non-equilibrium) E cannot drive the
  // Arrhenius block non-finite; see leaf_temp_min/max in the header.
  const double clamped = std::min(std::max(Tleaf, leaf_temp_min), leaf_temp_max);
  if (dT_dE != nullptr) {
    // Inside the clamp the balance is linear in E, so the slope is a constant
    // and negative: more transpiration, more latent heat, cooler leaf. ON the
    // clamp it is zero, because the returned temperature no longer responds to
    // E at all. The two branches share the one comparison above deliberately --
    // see the note on the declaration.
    *dT_dE = (clamped == Tleaf)
                 ? -latent_heat_vap * ra_ / vol_heat_cap_air
                 : 0.0;
  }
  return clamped;
}


// transpiration supply functions

// returns proportion of conductance taken from hydraulic vulnerability curve (unitless)
template <typename T>
inline T Leaf::proportion_of_conductivity_kernel(T psi) const {
  return exp(-pow((psi / stem_b), stem_c));
}

inline double Leaf::proportion_of_conductivity(double psi) const {
  return proportion_of_conductivity_kernel(psi);
}

// set spline for proportion of conductivity
inline void Leaf::setup_transpiration(double resolution) {
  std::vector<double> x_psi_, y_cumulative_transpiration_;
  build_cumulative_vulnerability_integral(stem_b, stem_c, resolution, x_psi_,
                                          y_cumulative_transpiration_);

  // THE CLOSED FORM SUPPLIES BOTH HALVES. G' is the surviving conductivity,
  // which this class already computes -- so the interpolant is handed a value and
  // a slope at every knot rather than choosing slopes from the values.
  //
  // ⚠️ THIS IS NOT AN OPTIMISATION, and letting the interpolant choose is a
  // measured regression. A shape-preserving rule (Fritsch-Carlson) limits slopes
  // to stop overshoot, which is the wrong objective when a gradient READS the
  // slope: at psi = 2 it lands 2,360x further from exp(-(psi/b)^c) than the fit it
  // replaces. Measured on the suites -- the spline tier moves 8.18e-05 against a
  // 1e-10 tolerance, test_leaf takes 7 failures, and the R suite goes 1484/0 to
  // 1458/36. With the slopes supplied it is 1484/0 and the spline tier is 2.57e-08.
  //
  // The inverse carries the reciprocal, d(psi)/dG = 1/G', finite everywhere the
  // forward slope is.
  std::vector<double> g_(x_psi_.size()), g_inv_(x_psi_.size());
  for (std::size_t i = 0; i < x_psi_.size(); ++i) {
    g_[i] = proportion_of_conductivity_kernel(x_psi_[i]);
    g_inv_[i] = 1.0 / g_[i];
  }
  transpiration_from_psi.init(x_psi_, y_cumulative_transpiration_, g_);
  transpiration_from_psi.set_extrapolate(false);

  psi_from_transpiration.init(y_cumulative_transpiration_, x_psi_, g_inv_);
  psi_from_transpiration.set_extrapolate(false);

  // The splines now describe the current stem_b, so the rescaling is over.
  // Recording it HERE rather than at each caller is what makes it impossible to
  // rebuild and forget.
  stem_b_spline_ = stem_b;
  ++stem_curve_builds_;
}

// --- the stem curve, as the four operations performed on it ------------------
//
// See the declarations for the identity these implement. Each returns the spline
// unchanged when nothing is rescaled, which is the production path.

// Named, scale-aware domain check in front of the spline read. Checked here
// rather than by catching odelia's throw and rethrowing, so the hot path carries
// no exception machinery: one comparison on a path that already does a spline
// solve.
//
// The domain is reported in the CALLER's units, not the spline's. Under a stem_b
// rescale (#46) the value handed to the spline is u/s, so the range the caller
// can actually reach is [s*min, s*max]; quoting the spline's own endpoints would
// send the reader after a discrepancy that is not there. At s == 1 -- the
// production path -- the two coincide.
//
// ⚠️ The comparison is `v < lo || v > hi` and NOT the negation of an in-range
// test, matching odelia's Interpolator::eval for the same reason: every
// comparison against NaN is false, so a non-finite u must keep falling through to
// the spline and coming back non-finite. Callers rely on it (plant documents a
// profit_psi_stem_TF(NA, .) -> NA contract built on exactly this), and negating
// an in-range test here would read as a tightening while being a behaviour
// change.
// Off the hot path: only reached when the lookup is about to fail, so the string
// building costs nothing in production. Kept out of line from eval_stem_curve so
// that function stays small enough to inline.
[[noreturn]] inline void stem_curve_out_of_domain(
    const odelia::interpolator::Interpolator& spline, double u, double v,
    double scale, const char* spline_name, const char* arg_name,
    const char* caller);

inline double Leaf::eval_stem_curve(const odelia::interpolator::Interpolator& spline,
                                    double u, double scale,
                                    const char* spline_name,
                                    const char* arg_name, const char* caller) {
  // scale == 1.0 is the production path (no stem_b rescale), and it must not pay
  // for the rescaled one: dividing and re-multiplying by 1.0 is exact but not
  // free, and collapsing the two cases into one measured +5% on
  // find_root_collar_psi. The check itself is a duplicate of the one inside
  // odelia's eval and costs nothing measurable -- it exists only to raise a
  // message that names this spline and its caller.
  if (scale == 1.0) {
    if (u < spline.min() || u > spline.max()) {
      stem_curve_out_of_domain(spline, u, u, scale, spline_name, arg_name, caller);
    }
    return spline.eval(u);
  }
  const double v = u / scale;
  if (v < spline.min() || v > spline.max()) {
    stem_curve_out_of_domain(spline, u, v, scale, spline_name, arg_name, caller);
  }
  return scale * spline.eval(v);
}

inline void stem_curve_out_of_domain(
    const odelia::interpolator::Interpolator& spline, double u, double v,
    double scale, const char* spline_name, const char* arg_name,
    const char* caller) {
    const bool below = v < spline.min();
    const double lo = scale * spline.min();
    const double hi = scale * spline.max();
    util::stop_infeasible("stem_curve_domain",
               std::string("Leaf hydraulics: ") + spline_name +
               " evaluated outside its domain: " + arg_name + " = " +
               util::format_double(u) + " lies " +
               util::format_double(below ? lo - u : u - hi) + " beyond the " +
               (below ? "lower" : "upper") + " end of [" +
               util::format_double(lo) + ", " + util::format_double(hi) + "]" +
               (scale == 1.0 ? std::string()
                             : "; the spline is rescaled by stem_b / "
                               "stem_b_spline = " +
                                   util::format_double(scale) +
                                   ", so it was read at " +
                                   util::format_double(v)) +
               (caller == nullptr ? std::string()
                                  : std::string("; asked by ") + caller) +
               ".");
}

inline double Leaf::stem_curve_integral(double psi, const char* caller) const {
  // x/x is exactly 1.0 for any finite non-zero x, so the equal case would fall
  // out of the division anyway; the branch is kept because stem_b_spline_ is 0.0
  // before the first setup_transpiration and 0/0 is not 1.
  const double s = (stem_b == stem_b_spline_) ? 1.0 : stem_b / stem_b_spline_;
  return eval_stem_curve(transpiration_from_psi, psi, s,
                         "transpiration_from_psi (the cumulative xylem "
                         "conductivity integral G, psi in +MPa)",
                         "psi", caller);
}

inline double Leaf::stem_curve_integral_deriv(double psi) const {
  if (stem_b == stem_b_spline_) {
    return transpiration_from_psi.deriv(psi);
  }
  return transpiration_from_psi.deriv(psi / (stem_b / stem_b_spline_));
}

inline double Leaf::stem_curve_integral_inverse(double w, const char* caller) const {
  const double s = (stem_b == stem_b_spline_) ? 1.0 : stem_b / stem_b_spline_;
  return eval_stem_curve(psi_from_transpiration, w, s,
                         "psi_from_transpiration (the INVERSE cumulative xylem "
                         "conductivity integral G^-1, argument in E/K_max)",
                         "E/K_max", caller);
}

inline double Leaf::stem_curve_integral_inverse_deriv(double w) const {
  if (stem_b == stem_b_spline_) {
    return psi_from_transpiration.deriv(w);
  }
  return psi_from_transpiration.deriv(w / (stem_b / stem_b_spline_));
}

inline void Leaf::perturb_stem_P50(double stem_P50_new) {
  check_psi_magnitudes(stem_P50_new, stem_c, roots_.root_P50, roots_.root_c);
  // psi_crit is a quantile of the same curve, so it is a fixed multiple of
  // stem_b at fixed stem_c: the ratio psi_crit/stem_b depends only on stem_c.
  // Both therefore scale by the same factor here, the rescaled spline evaluated
  // at the rescaled bound returns the same conductivity fraction, and the
  // domain cannot newly bind -- so the rescaling shortcut below stays exact and
  // no curve is rebuilt.
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
// on the pre-integrated spline that transpiration() uses on the hot path. Not
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

//calculates supply-side transpiration from psi_stem and opt_root_psi_, returns kg h20 s^-1 m^-2 LA
// SIGN: psi_stem and psi_upstream are POSITIVE magnitudes here (passed straight
// to the spline). Contrast transpiration_to_psi_stem below. See the sign-
// conventions block above.
inline double Leaf::transpiration(double psi_stem, double psi_upstream) {

  // 1-entry memo: identical (psi_stem, psi_upstream) is requested several times
  // per profit evaluation; return the cached value (bit-identical) to skip the
  // redundant spline lookups. Cache invalidated in set_physiology.
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

  transpiration_cache_psi_stem_ = psi_stem;
  transpiration_cache_psi_upstream_ = psi_upstream;
  transpiration_cache_value_ = E;
  transpiration_cached_ = true;
  return E;
}

// converts a known transpiration to its corresponding psi_stem, returns -MPa
// SIGN: unlike transpiration() above, psi_upstream here is a SIGNED (negative)
// potential, so it is flipped with a leading `-` before the spline lookup. The
// two functions are inverses called with opposite-sign psi_upstream.
inline double Leaf::transpiration_to_psi_stem(double transpiration_, double psi_upstream) {
  // integration of proportion_of_conductivity over [opt_root_psi_, psi_stem].
  // psi_upstream is a positive magnitude, same as in transpiration() -- #7 noted
  // that this function read eval(-psi_upstream) while its forward direction read
  // eval(psi_upstream), and that this was consistent because the two were called
  // with opposite-sign arguments. Under one representation they are not, and the
  // negation is deleted (#25).

  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +
    stem_curve_integral(psi_upstream, "Leaf::transpiration_to_psi_stem, at psi_upstream");

  // The lookup that killed plant#576, and the reason the inverse names its units:
  // E_psi_stem below the spline's lower end means the demanded flux is NEGATIVE,
  // i.e. the collar cannot supply it and the stem potential that would carry it is
  // wetter than saturation. There is no such potential, so widening the domain is
  // not the fix -- the caller should not have asked. Naming the argument E/K_max
  // rather than a bare `u` is what makes that readable from the message.
  return stem_curve_integral_inverse(
      E_psi_stem, "Leaf::transpiration_to_psi_stem, inverting for psi_stem");
  }

// returns stomatal conductance to CO2, mol C m^-2 LA s^-1
inline double Leaf:: stom_cond_CO2(double psi_stem, double psi_upstream) {
  double transpiration_ = transpiration(psi_stem, psi_upstream);
  return atm_kpa_ * transpiration_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
}


// biochemical photosynthesis model equations
//ensure that units of PPFD_ actually correspond to something real.
// electron trnansport rate based on light availability and vcmax assuming co-limitation hypothesis
inline double Leaf::electron_transport() {



  double electron_transport_ = (a * PPFD_ + jmax_ - sqrt(pow(a * PPFD_ + jmax_, 2) - 
  4 * curv_fact_elec_trans * a * PPFD_ * jmax_)) / (2 * curv_fact_elec_trans); // check brackets are correct

  // double electron_transport_ = (4*a*PPFD_)/sqrt(pow(4*a*PPFD_/jmax_,2)+ 1);
    return electron_transport_;           
}

//calculate the rubisco-limited assimilation rate, returns umol m^-2 s^-1
// The scalar-generic cores. Arithmetic and call structure are IDENTICAL to the
// `double` bodies these replaced -- same association, same `pow(s, 2)`, same
// two-call structure in the colimited form -- so the `double` instantiation is not
// meant to move results, and the golden file is the check on that. What moves is
// the DERIVATIVE, which now comes from this code instead of from the drifted
// replica.
template <typename T>
inline T Leaf::assim_rubisco_limited_kernel(T ci) const {
  return (vcmax_ * (ci - gamma_ * umol_per_mol_to_Pa_)) / (ci + km_);
}

template <typename T>
inline T Leaf::assim_electron_limited_kernel(T ci) const {
  return electron_transport_ / 4 *
  ((ci - gamma_ * umol_per_mol_to_Pa_) / (ci + 2 * gamma_ * umol_per_mol_to_Pa_));
}

template <typename T>
inline T Leaf::assim_colimited_kernel(T ci) const {
  T assim_rubisco_limited_ = assim_rubisco_limited_kernel(ci);
  T assim_electron_limited_ = assim_electron_limited_kernel(ci);

  return (assim_rubisco_limited_ + assim_electron_limited_ - sqrt(pow(assim_rubisco_limited_ + assim_electron_limited_, 2) - 4 * curv_fact_colim * assim_rubisco_limited_ * assim_electron_limited_)) /
             (2 * curv_fact_colim)- R_d_;
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
inline double Leaf::psi_stem_to_ci(double psi_stem, double psi_upstream) {
  const double stom_cond_CO2_fixed = stom_cond_CO2(psi_stem, psi_upstream);

  // Propagate non-finite inputs as NA rather than entering the solver. A
  // non-finite psi_stem (e.g. NA from profit_psi_stem_TF(NA, .)) makes gc and
  // hence the whole target non-finite. The previous bisection returned NaN
  // silently in this case; the bracketing TOMS748 solver below instead throws
  // ("a and b do not bracket the root"), so guard explicitly to preserve the
  // NA-in -> NA-out contract (see test-leaf.r "Basic functions").
  if (!std::isfinite(stom_cond_CO2_fixed)) {
    return ci_ = util::na_value;
  }

  auto target = [&](double x) mutable -> double {
    const double assim_colimited_x_ = assim_colimited(x);
    return assim_colimited_x_ * umol_to_mol -
      (stom_cond_CO2_fixed * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
  };

  // #486: this target (assim_colimited demand minus the linear gc supply) is
  // smooth and strictly monotone over (gamma*, ca] -- no singularity at the
  // bracket ends (Ar,Ae vanish linearly at gamma* so sqrt(disc) is linear, not
  // singular) and its only sharp feature is the colimitation elbow far below the
  // operating root. So it is a well-behaved case for a superlinear bracketing
  // solver: TOMS748 reaches the same root in ~9 evals vs bisection's ~29. This is
  // deliberately scoped to psi_stem_to_ci ONLY; the hydraulic find_root_psi path
  // keeps bisection (its target is not smooth -- see the warning on
  // util::uniroot_smooth).
  //
  // ⚠️ **The 1e-10 is load-bearing and was chosen by measurement, not taste. It
  // sets the floor of what every reported output of this model MEANS.** It used to
  // be 1e-7, which was invisible while golden section's GSS_tol_abs (1e-3)
  // dominated; with that gone, this became the model's dominant
  // amplifier -- a last-bit change anywhere upstream shifts which point TOMS748
  // lands on inside its tolerance band, and that surfaces in ci, assim, gc and
  // profit. Measured by perturbing the model identically at each tolerance and
  // reading the golden grid:
  //
  //   ci tol   worst diff vs a converged (1e-15) solve   us/solve (indicative)
  //   1e-7                       1.82e-07                      2.655
  //   1e-8                       3.15e-08                      2.640
  //   1e-10                      5.41e-10                      2.695   <- here
  //   1e-11                      5.41e-10                      2.768
  //   1e-13                      7.48e-13                      2.885
  //   1e-15                      0                             2.945
  //
  // 1e-10 is the knee: it lands **335x closer to a converged solve** than 1e-7
  // did, where going on to 1e-13 buys ~700x more for roughly three times the extra
  // time. The plateaus (1e-8 = 1e-9, 1e-10 = 1e-11) are TOMS748's discrete
  // iterations, so tightening *between* them buys nothing and still costs.
  //
  // ⚠️ **The precision column is solid; the timing column is only indicative and
  // reads ~2x too cheap.** Those seven binaries were built in one loop in a
  // scratch tree, and code-layout luck moves this benchmark by ~2%. The controlled
  // figure for the step actually taken -- both binaries built from the same tree,
  // interleaved x5 -- is **+3.4%** (2.65 -> 2.75 us/solve), against the +1.5% the
  // table implies. Still a good trade next to the root-find's 24.5%, and still
  // 21.7% faster than the 3.51 us this package ran at before 11a. If you re-derive
  // the curve, build every point from one tree.
  //
  // Do not tighten further without a use that needs it, and do not loosen it back
  // without re-reading the guide's rounding-versus-bug magnitudes, which this
  // figure sets.
  //
  // ⚠️ Not the same knob as `ci_abs_tol` (the settable control, default 1e-3).
  // `ci_abs_tol` is read in exactly one place,
  // solve_medlyn_ci_numerical, i.e. the empirical Medlyn-conductance route. It
  // reaches NEITHER family of optimality solver -- both go through this function
  // and so through the 1e-10 above. A caller tightening it gets no extra precision
  // anywhere in the optimality model.
  ci_at_compensation_point_ = false;
  try {
    return ci_ = util::uniroot_smooth(target, gamma_ * umol_per_mol_to_Pa_, ca_, 1e-10, ci_niter);
  } catch (const std::exception& e) {
    // Assimilation can be negative across the WHOLE [gamma*, ca] bracket -- the
    // leaf is too hot to gain carbon at any internal CO2 -- and then there is no
    // supply==demand root at all. That is a physically-meaningful shut-down, not a
    // solver failure: operate at the compensation point (ci = gamma*, gross A = 0,
    // net A = -R_d) and let the profit optimiser move away from it.
    //
    // ⚠️ THIS USED TO BE GATED ON use_energy_balance_, on the grounds that the
    // prescribed-temperature path "never reaches here". That was true and #41 made
    // it false: raising R_d to a Q10 response is enough to drive assimilation
    // negative throughout by ~45 C at the defaults, and the prescribed path then
    // threw where the PM path shut down -- the same physics, two different
    // outcomes, decided by how leaf_temp happened to be obtained.
    //
    // ⚠️ THE GATE IS NOT SIMPLY REMOVED, because that would mask real solver
    // failures as shut-downs: this `catch` sees ANY exception from the root-find,
    // and fail-fast has value. The target is strictly monotone over (gamma*, ca]
    // (#486), so "no root exists" is distinguishable from "the solver broke" by
    // the sign at the two ends -- equal signs means the root is genuinely outside
    // the bracket. Only that case shuts down; anything else still stops.
    const double lo = gamma_ * umol_per_mol_to_Pa_;
    const double t_lo = target(lo);
    const double t_hi = target(ca_);
    const bool no_root_in_bracket =
        std::isfinite(t_lo) && std::isfinite(t_hi) &&
        ((t_lo > 0.0 && t_hi > 0.0) || (t_lo < 0.0 && t_hi < 0.0));
    if (use_energy_balance_ || no_root_in_bracket) {
      ci_at_compensation_point_ = true;
      return ci_ = lo;
    }
    util::stop_infeasible("ci_solve", "psi_stem_to_ci failed: " + std::string(e.what()) +
               "; min=" + util::to_string(gamma_ * umol_per_mol_to_Pa_) +
               "; max=" + util::to_string(ca_) +
               "; psi_stem=" + util::to_string(psi_stem) +
               "; psi_upstream=" + util::to_string(psi_upstream));
  }
}

// given psi_stem, find assimilation, transpiration and stomal conductance to c02
inline void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream) {

  // ⚠️ AN `assim_max_ < 0` EARLY EXIT USED TO SIT HERE, AND IT ZEROED THE WATER.
  // It read: if gross assimilation at ci = ca cannot cover respiration, put ci at
  // the compensation point and set transpiration and conductance to zero. The
  // first half is right and is now the ci solver's own job; the second half is
  // wrong, and wrong in a way that matters.
  //
  // Transpiration on this path is the HYDRAULIC SUPPLY at the candidate potential.
  // It does not ask whether there is carbon to be had -- water moves down a
  // potential gradient whatever photosynthesis is doing. Zeroing it made two leaves
  // at the SAME operating point disagree about whether water was moving purely
  // because one of them had respiration switched on, which is how it was found.
  //
  // Deleting it changes nothing about the shut-down STATE it was reaching for.
  // Where `assim_max_ < 0` the ci root-find has no supply==demand root anywhere in
  // [gamma*, ca] -- assimilation is monotone in ci and negative at both ends -- so
  // psi_stem_to_ci takes its compensation-point fallback and returns ci = gamma*
  // with gross assimilation zero and net = -R_d, which is exactly what this branch
  // set by hand. That fallback did not exist when the branch was written; it was
  // added with the Q10 respiration in #41, which is also what made this regime
  // reachable at ordinary light.
  //
  // ⚠️ AND IT IS UNREACHABLE FROM THE COLLAR SOLVE, which is why plant is
  // unaffected: prepare_collar_solve has its OWN `assim_max_ < 0` exit and returns
  // false before any candidate potential is evaluated. The golden grid's minimum
  // assim_max_ is 3.71, so the file never reaches this branch either and is
  // bit-identical across the change -- which is a statement about the grid, not
  // evidence that the change is inert. test_transpiration_survives_negative_assim
  // is the test that can see it.
  //
  // ⚠️ `Tleaf_` IS WRITTEN ON BOTH PATHS BELOW, which is hazard 8 applied to the
  // one output whose value is a driver on one path and a solved quantity on the
  // other. The remaining shut-down branch transpires nothing, and on the PM path
  // a leaf that transpires nothing is the HOTTEST one -- so leaving `Tleaf_` to
  // the transpiring branch would report the previous candidate's temperature for
  // exactly the operating points where the answer is most extreme. (#106 wrote
  // this on all THREE branches it found; deleting the `assim_max_ < 0` exit
  // leaves two.)
  if (psi_upstream >= psi_stem){
    // ⚠️ THE TEMPERATURE BLOCK MUST BE SEATED BEFORE ANYTHING BELOW READS IT, so
    // these three lines stay above the assignments and in this order. `ci_` reads
    // `gamma_`, `assim_colimited_` at the bottom of this function reads the whole
    // block, and `set_leaf_vpd` writes a field that is both reported and read by
    // `g1_eff()`. This branch transpires nothing, so its temperature is the E = 0
    // one -- on the PM path the HOTTEST the leaf gets, since nothing is cooling it.
    //
    // ⚠️ THE BLOCK ARRIVES HOLDING THE PREVIOUS CANDIDATE'S TEMPERATURE, not the
    // baseline: the transpiring branch below re-derives it per candidate by design.
    // So a read here without a preceding write is order-dependent, not merely
    // offset, and hazard 3 (the argmax must vary smoothly with inputs) is what
    // that costs. Anything added to this branch that reads the block must sit
    // below these lines.
    Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(0.0) : leaf_temp_;
    if (use_energy_balance_) {
      update_temperature_dependent_params(Tleaf_);
      set_leaf_vpd(Tleaf_);
    }
    ci_ = gamma_*umol_per_mol_to_Pa_;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    } else{
      {
      // Transpiration is the hydraulic supply, independent of ci; compute it
      // first so the PM path can derive the operating-point leaf temperature.
      // Off the PM path this is a memoised no-op reorder (psi_stem_to_ci ->
      // stom_cond_CO2 requests the same (psi_stem, psi_upstream), returning the
      // bit-identical cached value), so the non-PM result is unchanged.
      transpiration_ = transpiration(psi_stem, psi_upstream);
      // Tleaf = f(E) is explicit (no PM inversion, no A->E feedback), so this is
      // a single forward pass: recompute the Farquhar temperature params at this
      // candidate's Tleaf before solving for ci. Defeats the photo_temp cache by
      // design -- Tleaf varies per operating point.
      //
      // Stored rather than passed straight through, so the temperature the
      // parameters were derived at and the temperature reported to the caller
      // cannot be two different numbers. No extra `leaf_temp_from_E` call: this
      // is the hot path, ~10^3 candidates per solve.
      Tleaf_ = use_energy_balance_ ? leaf_temp_from_E(transpiration_) : leaf_temp_;
      if (use_energy_balance_) {
        update_temperature_dependent_params(Tleaf_);
        // ...and the deficit Fick's law divides by, which moves with it. Both
        // read the STORED `Tleaf_` rather than re-deriving it, so the temperature
        // the Farquhar parameters were computed at, the one the deficit was
        // computed at, and the one reported to the caller cannot be three
        // different numbers -- and it costs no extra `leaf_temp_from_E`.
        set_leaf_vpd(Tleaf_);
      }
      ci_ = psi_stem_to_ci(psi_stem, psi_upstream);
      stom_cond_CO2_ = atm_kpa_ * transpiration_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
      }
    }
  assim_colimited_ = assim_colimited(ci_);
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
// which a product objective can never be boundary-crit -- water becomes infinitely
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
    const double E = transpiration(psi_stem, supply_psi_soil_scalar());
    double dT_dE = 0.0;
    const double Tleaf = leaf_temp_from_E(E, &dT_dE);
    const double span = T50_ - Tcrit_;
    const double TC = thermal_cost_at(Tleaf);
    // Sigmoid derivative: dTC/dT = r*TC*(1-TC), r = 2/(T50 - Tcrit).
    thermal = (2.0 / span) * TC * (1.0 - TC) * dT_dE;
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
template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(T psi_stem) const {
  return TF24_cost_scale * pow((1 - proportion_of_conductivity_kernel(psi_stem)), TF24_beta2);
}

// ---- the same kernels, against a caller's parameter pack -------------------
//
// Each is its member-reading twin above with every differentiated quantity taken
// from `pars` instead. At double the pack holds exactly what the members hold
// (passive_pars builds it from them), and the operation order is unchanged, so
// these instantiate to the same numbers bit for bit.

template <typename T>
inline T Leaf::vcmax_at(const leaf_pars<T>& pars) const {
  return peak_arrh_curve<T>(T(vcmax_ha_), pars[par_vcmax_25],
                            T(leaf_temp_), T(vcmax_H_d_), T(vcmax_d_S_));
}

template <typename T>
inline T Leaf::jmax_at(const leaf_pars<T>& pars) const {
  T j = peak_arrh_curve<T>(T(jmax_ha_), pars[par_jmax_25],
                           T(leaf_temp_), T(jmax_H_d_), T(jmax_d_S_));
  // The thermal cost reads leaf temperature and two thresholds, none of which is
  // a trait, so it is a factor here rather than a term carrying rows.
  if (use_thermal_cost_) {
    // Written out rather than `*=`, whose overload is ambiguous at a nested
    // scalar.
    j = j * T(1.0 - thermal_cost_at(leaf_temp_));
  }
  return j;
}

template <typename T>
inline T Leaf::respiration_at(const leaf_pars<T>& pars) const {
  const double q10 =
      rd_q10_intercept_ - rd_q10_slope_ * (leaf_temp_ + 25.0) / 2.0;
  return pars[par_R_d_25] *
         T(std::pow(q10, (leaf_temp_ - 25.0) / 10.0));
}

template <typename T>
inline T Leaf::electron_transport_at(const leaf_pars<T>& pars) const {
  const T a_at = pars[par_a];
  const T jm = jmax_at<T>(pars);
  const T curv = pars[par_curv_fact_elec_trans];
  return (a_at * PPFD_ + jm -
          sqrt(pow(a_at * PPFD_ + jm, 2) - 4 * curv * a_at * PPFD_ * jm)) /
         (2 * curv);
}

template <typename T>
inline T Leaf::assim_colimited_kernel(const T& ci,
                                      const leaf_pars<T>& pars) const {
  // gamma_ and km_ are temperature responses of fixed physiology -- no trait
  // reaches either -- so they stay members.
  const T vc = vcmax_at<T>(pars);
  const T J = electron_transport_at<T>(pars);
  const T curv = pars[par_curv_fact_colim];
  const T rubisco = (vc * (ci - gamma_ * umol_per_mol_to_Pa_)) / (ci + km_);
  const T electron =
      J / 4 *
      ((ci - gamma_ * umol_per_mol_to_Pa_) /
       (ci + 2 * gamma_ * umol_per_mol_to_Pa_));
  return (rubisco + electron -
          sqrt(pow(rubisco + electron, 2) - 4 * curv * rubisco * electron)) /
             (2 * curv) -
         respiration_at<T>(pars);
}

template <typename T>
inline T Leaf::transpiration_at(const T& sigma, const T& collar,
                                const leaf_pars<T>& pars) const {
  return pars[par_kmax] * (stem_integral_at<T>(sigma, pars) -
                           stem_integral_at<T>(collar, pars));
}

template <typename T>
inline T Leaf::assim_slope_at(const T& ci, const leaf_pars<T>& pars) const {
  // ⚠️ A NESTED ACTIVE IS BUILT BY ASSIGNING INTO ITS VALUE, not by converting.
  // The outer scalar's constructor takes the INNERMOST value type, so
  // `nested(inner)` either refuses to compile or strips the inner rows -- and a
  // stripped row here is a slope that is right at double and zero in every trait.
  using nested = odelia::ode::tangent_scalar<T>;
  auto lift = [](const T& x) -> nested {
    nested out;
    xad::value(out) = x;
    return out;
  };
  leaf_pars<nested> np;
  for (std::size_t i = 0; i < np.size(); ++i) np[i] = lift(pars[i]);
  nested c = lift(ci);
  odelia::ode::seed_direction(c, 1.0);
  return odelia::ode::derivative_along(
      assim_colimited_kernel<nested>(c, np));
}

// The curve itself is closed form, so the (P50, c) chain comes out of the one
// derivation rather than a graft: only the pre-integrated G(psi) is tabulated and
// needs its rows supplied (stem_integral_at).
template <typename T>
inline T Leaf::proportion_of_conductivity_kernel(const T& psi,
                                                 const leaf_pars<T>& pars) const {
  const T c = pars[par_stem_c];
  const T b = phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], c);
  return exp(-pow((psi / b), c));
}

template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem,
                                        const leaf_pars<T>& pars) const {
  return pars[par_TF24_cost_scale] *
         pow((1 - proportion_of_conductivity_kernel<T>(psi_stem, pars)),
             pars[par_TF24_beta2]);
}

template <class S>
inline Leaf::CollarCoords<S> Leaf::collar_coords_at(
    double sigma_star, double ci_star, const S& collar,
    const SupplyDraw<S>& draw, const leaf_pars<S>& pars) const {
  using odelia::util::to_passive;
  check_draw(to_passive(collar), draw);

  // ⚠️ THE LEAF-TO-AIR DEFICIT AND THE SETTABLE RATIO, which is what
  // stom_cond_CO2 divides by. With the energy balance off vpd_leaf_ IS atm_vpd_,
  // so a check with the gate off cannot tell them apart -- and with it on they
  // differ by up to 1.52x.
  const double gc_per_flux =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const Leaf& leaf = *this;
  const leaf_pars<double> at_pars = passive_pars();

  // The flux at the LIVE collar. The draw holds it at the passive one carrying
  // its own rows; what is missing is the channel through the collar moving, and
  // that is the one recorded slope. At the operating point the step is exactly
  // zero, so this is the draw's flux bit for bit.
  const S step = collar - S(draw.at);
  const S flux = draw.flux.value + draw.flux.slope * step;

  // dT1/dsigma is one term of the residual: the stem curve at the operating
  // point, scaled by kmax. Taken at the passive point, which is where the theorem
  // divides.
  const double dT1_dsigma =
      at_pars[par_kmax] *
      proportion_of_conductivity_kernel<double>(sigma_star, at_pars);
  const S sigma = odelia::implicit_value<S>(
      sigma_star, dT1_dsigma, [&](const S& sg) -> S {
        return pars[par_kmax] * (leaf.template stem_integral_at<S>(sg, pars) -
                                 leaf.template stem_integral_at<S>(collar, pars)) -
               flux;
      });

  // ⚠️ FOLLOW THE BRANCH THE FORWARD MODEL TOOK. Where the leaf is not moving
  // water it does not place ci by the residual at all -- it assigns the CO2
  // compensation point, which is a function of temperature and which no carbon
  // input reaches. Applying the theorem to a condition the model did not solve
  // returns a number rather than an error, and the number is a cancellation of
  // two vanishing terms.
  //
  // dT2/dci is the assimilation's own slope plus the supply term's, at the
  // passive point. A' comes from a TANGENT THROUGH THE SAME KERNEL the residual
  // calls, which is upstream's own way of taking it -- a hand-written slope twin
  // would be a second definition of one function, free to disagree with it in a
  // way that stays finite and plausible.
  // The stem's flux at the placed sigma. gc is built from THIS, not from the
  // draw, because that is what the model's ci solve reads.
  const S stem_flux = transpiration_at<S>(sigma, collar, pars);

  double dT2_dci = 0.0;
  if (!ci_at_compensation_point_) {
    dT2_dci = assim_slope_at<double>(ci_star, at_pars) * umol_to_mol +
              gc_per_flux * to_passive(stem_flux) * inv_atm;
  }
  const S ci =
      ci_at_compensation_point_
          ? S(ci_star)
          : odelia::implicit_value<S>(
                ci_star, dT2_dci, [&](const S& c) -> S {
                  const S A = leaf.template assim_colimited_kernel<S>(c, pars);
                  return A * umol_to_mol -
                         S(gc_per_flux) * stem_flux * (S(leaf.ca_) - c) *
                             S(inv_atm);
                });

  // ---- the other half of each pair: the response in the collar --------------
  //
  // Both come out of the same two conditions, differentiated in the collar with
  // the parameters held still. They are returned WITH their values because a
  // consumer seeded from separate places can pair a value from one point with a
  // slope from another, and a first attempt that moved only the collar disagreed
  // with a difference by 40 to 100 per cent.
  // ⚠️ THE CONDUCTIVITY HERE IS THE TABLE'S, WITH THE CURVE'S ROWS. Upstream's
  // dprofit_at_collar_psi forms every one of these from stem_curve_integral_deriv,
  // and the solve drove THAT assembly to zero -- so a closed-form value here
  // makes M non-zero at the collar the solve placed, and the envelope omission
  // stops being true. graft.hpp is the same rule and the same one place: the
  // value is the table's because the solve ran on the table, the rows are the
  // curve's because a derivative off a table is a derivative of the fit.
  const S kmax = pars[par_kmax];
  const S f_p = graft_curve<S>(stem_curve_integral_deriv(to_passive(collar)),
                               collar, pars[par_stem_P50], pars[par_stem_c]);
  const S f_sigma = graft_curve<S>(stem_curve_integral_deriv(sigma_star), sigma,
                                   pars[par_stem_P50], pars[par_stem_c]);

  // dsigma/dcollar: THE INVERSE SPLINE'S OWN SLOPE times the transport's, which
  // is how upstream's dprofit_at_collar_psi forms it and what actually placed
  // sigma -- transpiration_to_psi_stem inverts through psi_from_transpiration,
  // a separate table.
  //
  // ⚠️ 1/f(sigma) IS THAT SLOPE ONLY IF THE TWO STEM SPLINES ARE EXACT MUTUAL
  // INVERSES, AND THEY ARE NOT. Writing it that way puts sigma's response ~8 per
  // cent out at the shipped resolution, which the marginal then amplifies by its
  // own cancellation: M read 6.4e-02 at a collar where the solve had driven
  // upstream's dprofit to 1e-12, so the envelope omission -- dropping
  // M * dp/dtheta at an interior point BECAUSE M is zero there -- stopped being
  // true of the number this function returns.
  //
  // So the VALUE is the table's, which is what makes M vanish at the collar the
  // solve placed, and the ROWS are the closed form's, which is graft.hpp's rule:
  // a derivative read off a table is a derivative of the fit.
  // ⚠️ THE SOIL DRAW, NOT THE STEM FLUX. transpiration_to_psi_stem inverts at
  // E_up/kmax + G(collar), so this is the point the inverse table was read at.
  // The two fluxes differ by that table's own round-trip, which is the whole
  // reason this argument has to be the one sigma was actually placed from.
  const double u_at = to_passive(draw.flux.value) / to_passive(kmax) +
                      to_passive(stem_integral_at<S>(collar, pars));
  const S recip = S(1.0) / f_sigma;
  const S dsigma_du =
      S(stem_curve_integral_inverse_deriv(u_at)) + (recip - S(to_passive(recip)));
  const S V = dsigma_du * (draw.flux.slope / kmax + f_p);

  const S gc_const = S(gc_per_flux);
  const S inv = S(inv_atm);
  const S gc = gc_const * stem_flux;
  // The conductance partials are the transport slope at each end, which is the
  // closed form kmax * f -- already formed above, so read rather than recomputed.
  const S dgc_dsigma = gc_const * kmax * f_sigma;
  const S dgc_dp = -gc_const * kmax * f_p;

  const S A_prime = assim_slope_at<S>(ci, pars);
  const S g_ci = A_prime * umol_to_mol + gc * inv;
  const S ca_minus_ci = S(ca_) - ci;
  const S dci_dsigma = (dgc_dsigma * ca_minus_ci * inv) / g_ci;
  const S dci_dp_expl = (dgc_dp * ca_minus_ci * inv) / g_ci;

  return CollarCoords<S>{pair<S>{sigma, V},
                         pair<S>{ci, dci_dsigma * V + dci_dp_expl}};
}

template <class S>
inline S Leaf::psi_crit_at(const leaf_pars<S>& pars) const {
  const S c = pars[par_stem_c];
  const S b = phylloptim::weibull_b_from_P50<S>(pars[par_stem_P50], c);
  return b * pow(S(std::log(1.0 / k_crit_fraction)), 1.0 / c);
}

template <Leaf::CostCurve K, class S>
inline S Leaf::profit_at(const S& collar, const SupplyDraw<S>& draw,
                         const leaf_pars<S>& pars) const {
  static_assert(K == CostCurve::TF24,
                "the reverse-mode surface covers TF24; a second curve is a row "
                "in profit_at, outputs_at and marginal_at");
  const CollarCoords<S> at =
      collar_coords_at<S>(opt_psi_stem_, ci_, collar, draw, pars);
  const S A = assim_colimited_kernel<S>(at.ci.value, pars);
  const S cost = hydraulic_cost_TF_kernel<S>(at.sigma.value, pars);
  return A - cost;
}

template <Leaf::CostCurve K, class S>
inline Leaf::LeafOutputs<S> Leaf::outputs_at(const S& collar,
                                             const SupplyDraw<S>& draw,
                                             const leaf_pars<S>& pars) const {
  using odelia::util::to_passive;
  check_draw(to_passive(collar), draw);
  LeafOutputs<S> out;
  const std::size_t n = static_cast<std::size_t>(supply_n_layers());

  // The draws first. Profit reads the same flux wherever its collar is the live
  // one, so taking them here is what stops the supply reaching the tape twice.
  if (operating_point_kind_ == OperatingPointKind::HydraulicShutdown) {
    // Every flux is zero and stays zero: nothing the soil holds reaches a leaf
    // that has stopped drawing on it.
    out.uptake.assign(n, S(0.0));
  } else {
    // ⚠️ THE DRAWS ARE GRAFTED, NOT TAKEN AGAIN. The draw already holds every
    // layer's uptake at the passive collar carrying the state's own rows; what
    // is missing is the channel through the collar MOVING, and that is one
    // supplied slope per layer. Taking the supply again at the live collar
    // records the whole of it a second time for the same numbers -- the step is
    // exactly zero in value, so every uptake here is the draw's, bit for bit.
    const S step = collar - S(draw.at);
    out.uptake = draw.uptake;
    odelia::util::check_length(draw.duptake_dp.size(), out.uptake.size());
    for (std::size_t j = 0; j < out.uptake.size(); ++j) {
      out.uptake[j] = out.uptake[j] + S(draw.duptake_dp[j]) * step;
    }
  }

  // ⚠️ PROFIT READS THE HELD COLLAR AT AN INTERIOR POINT, and that is the
  // envelope theorem written as an omission rather than as a term that has to
  // come out to zero. The marginal is zero there by the condition the solve
  // drove to nothing, so carrying the collar's movement through profit would add
  // M * dp/dtheta -- a term whose size is the solve's tolerance rather than the
  // model's. It is also why a collar that cannot be recorded costs the water
  // rows and not the objective.
  const bool interior = operating_point_kind_ == OperatingPointKind::Interior;
  const S collar_held(to_passive(collar));

  if (operating_point_kind_ == OperatingPointKind::HydraulicShutdown ||
      operating_point_kind_ == OperatingPointKind::ShadeDeath) {
    // No water moves, so no carbon is fixed: assimilation is minus the dark
    // respiration, and the cost is the one at the potential being held.
    const S held = operating_point_kind_ == OperatingPointKind::HydraulicShutdown
                       ? psi_crit_at<S>(pars)
                       : collar;
    out.profit = -respiration_at<S>(pars) - hydraulic_cost_TF_kernel<S>(held, pars);
  } else {
    out.profit = profit_at<K, S>(interior ? collar_held : collar, draw, pars);
  }
  return out;
}

template <Leaf::CostCurve K, class S>
inline S Leaf::marginal_at(const S& collar, const SupplyDraw<S>& draw,
                           const leaf_pars<S>& pars) const {
  static_assert(K == CostCurve::TF24,
                "the reverse-mode surface covers TF24; a second curve is a row "
                "in profit_at, outputs_at and marginal_at");
  // ⚠️ EVERY FACTOR IS A VALUE, so this is arithmetic at S and dM/dp is a FIRST
  // derivative of it. The lift this replaces took its value and first derivative
  // from a TABULATION and its second from the TRUE function, and the two disagree
  // by the interpolation error -- which is harmless at first order and is what
  // flipped the sign of a curvature where the collar came within 5.6e-08 of a
  // soil layer's potential.
  const CollarCoords<S> at =
      collar_coords_at<S>(opt_psi_stem_, ci_, collar, draw, pars);
  const S A_prime = assim_slope_at<S>(at.ci.value, pars);
  // dC/dpsi_stem from upstream's own arm, at a scalar: a tangent through the
  // cost kernel rather than a slope written out beside it.
  using nested = odelia::ode::tangent_scalar<S>;
  auto lift = [](const S& x) -> nested {
    nested out;
    xad::value(out) = x;
    return out;
  };
  leaf_pars<nested> np;
  for (std::size_t i = 0; i < np.size(); ++i) np[i] = lift(pars[i]);
  nested sg = lift(at.sigma.value);
  odelia::ode::seed_direction(sg, 1.0);
  const S C_prime =
      odelia::ode::derivative_along(hydraulic_cost_TF_kernel<nested>(sg, np));

  return A_prime * at.ci.slope - C_prime * at.sigma.slope;
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
// `transpiration()` is memoised on (psi_stem, psi_upstream), so calling it here
// straight after set_leaf_states_rates_from_psi_stem costs nothing and returns
// exactly the E of the operating point that call established -- including the
// no-flow case, where the two potentials are equal and the integral difference is
// exactly zero.
inline double Leaf::hydraulic_cost_CF77(double psi_stem,
                                                double psi_upstream) {
  hydraulic_cost_ = cf77_price() * transpiration(psi_stem, psi_upstream);
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
// `transpiration()` is memoised on (psi_stem, psi_upstream), so calling it here
// straight after set_leaf_states_rates_from_psi_stem costs nothing -- see
// hydraulic_cost_CF77, which relies on the same thing.
inline double Leaf::hydraulic_cost_TF24_floor(double psi_stem,
                                           double psi_upstream) {
  hydraulic_cost_ = hydraulic_cost_TF_kernel(psi_stem) +
                    TF24_floor_lambda_o * transpiration(psi_stem, psi_upstream);
  return hydraulic_cost_;
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
// uses, so this needs neither AD nor a spline derivative; `f` decreases, so this is
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
  // MINIMUM of A is where A is negative, and A is negative only at full closure,
  // where E = 0 and A = -R_d -- the one point excluded below. Away from it A is the
  // carbon the leaf actually gains, so |A|max is max A.
  //
  // dA/dpsi needs no new algebra: for an identity link dJ/dpsi = A'*dci/dpsi - C',
  // so adding the cost derivative back recovers the benefit term exactly.
  //
  // ⚠️ THE CLOSURE POINT IS EXCLUDED, as it was from the scan. At exactly psi_soil
  // E is zero and A is -R_d: a closed stoma rather than a point on the stream, and
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
  const double width = psi_crit - psi_soil;
  const double lo = psi_soil + 1e-9 * width;   // the open end: closure excluded
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
  const double step = (psi_crit - psi_soil) / double(n - 1);
  std::vector<double> out(static_cast<std::size_t>(5 * n));
  for (int i = 0; i < n; ++i) {
    const double p = psi_soil + step * double(i);
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

inline void Leaf::set_profitmax_degenerate_state() {
  profit_ = 0;
  transpiration_ = 0;
  stom_cond_CO2_ = 0;
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
  if constexpr (K == CostCurve::SOX || K == CostCurve::JW26) {
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
    // ⚠️ CF77_lambda_ IS AN INPUT AND HAS NO DEFAULT. It is the marginal value of
    // water -- the one parameter this model is defined by -- and no model code
    // supplies one, so a caller who has not set it is maximising a NaN objective.
    // The potential that comes back from that is a property of the bracket rather
    // than of the leaf, and it looks entirely plausible. Refuse instead.
    if (!std::isfinite(CF77_lambda_)) {
      util::stop("the CF77 cost curve needs CF77_lambda_ set: it is the "
                 "PRESCRIBED marginal value of water in umol C (kg H2O)^-1, NA "
                 "until you assign one, and never set by set_physiology or "
                 "set_traits.");
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
    // ⚠️ ONE PARAMETER, AND IT IS AN INPUT WITH NO DEFAULT -- exactly like
    // CF77_lambda_ and for a sharper reason: zero is a perfectly good value here,
    // and it is the value at which this curve stops being itself. The hydraulic
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
  } else {
    static_assert(K == CostCurve::TF24 || K == CostCurve::SOX ||
                  K == CostCurve::JW26 || K == CostCurve::ProfitMax,
                  "unhandled CostCurve");
  }
}


// Which objective. ⚠️ These are NOT all the same kind of number: five are carbon
// (TF24, CF77, JS22, CMax, TF24_floor) and two are products (see
// `profit_psi_stem_SOX`). That is
// already true of `profit_` across the optimisers and is not introduced here.
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
  } else {
    static_assert(K == CostCurve::ProfitMax, "unhandled CostCurve");
    (void)psi_upstream;
    return lambda_ProfitMax(psi_stem);
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
// TF24_floor -- TF24's cost plus CF77's price, and the only curve here whose marginal
//   cost of water does NOT vanish at the wet end. That is the whole reason it
//   exists: every conductance-loss curve in this package has `lambda -> 0` as the
//   drop closes, which prices water as free precisely when it is abundant.
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
  opt_psi_stem_ = psi_soil;

  // ⚠️ ONE CURVE NEEDS SETUP BEFORE THE SEARCH, AND ITS DEGENERATE EXIT DIFFERS.
  // ProfitMax divides by |A|max, so it is seeded here -- and where the normalisers
  // come out unusable there is no normalised objective to evaluate at the no-flow
  // point, which is why that arm ZEROES where every other curve EVALUATES. `PPFD_`
  // joins the test because a dark leaf has A = -R_d everywhere, so |A|max is
  // non-positive and the second bail-out would catch it anyway; testing it first
  // avoids seeding in order to discover that.
  if constexpr (K == CostCurve::ProfitMax) {
    if ((PPFD_ < 1.5e-8) | (psi_soil > psi_crit)) {
      set_profitmax_degenerate_state();
      return;
    }
    prepare_profitmax();
    if (!(profitmax_k_span_ > 0.0) || !(profitmax_A_max_ > 0.0)) {
      set_profitmax_degenerate_state();
      return;
    }
  } else if (psi_soil > psi_crit) {
    // Drier soil than the stem can reach: one feasible potential, and this is it.
    profit_ = profit_psi_stem_for<K>(psi_soil, psi_soil);
    lambda_emergent_ = lambda_for<K>(psi_soil, psi_soil);
    return;
  }

  double profit_opt = 0.0;
  const auto objective = [&](double psi_stem) {
    return profit_psi_stem_for<K>(psi_stem, psi_soil);
  };
  // ⚠️ WHETHER A DERIVATIVE EXISTS IS DECIDED HERE, NOT BY CATCHING A THROW. Under
  // a non-identity benefit link with the energy balance on, dJ/dpsi is REFUSED
  // rather than approximated (the temperature term needs the same h'(A) factor).
  // The OBJECTIVE is still perfectly well defined there, so that configuration
  // keeps the bracket-width refinement instead of losing the solve to an
  // exception -- letting the refusal fire from inside the optimiser turned 3236 of
  // 4608 golden rows into throws.
  // ONE solver, and NO caller-side branch on whether a derivative exists: the
  // derivative closure reports that itself through `ok`, and the solver's own
  // fallback handles it. A non-identity link with the energy balance on refuses
  // dJ/dpsi (the temperature term needs the same h'(A) factor), and that refusal is
  // reported rather than thrown, because the OBJECTIVE is perfectly well defined
  // there. `basin_scan_cells()` already returns a scan whenever the gate is on,
  // which is exactly when the fallback needs one.
  opt_psi_stem_ = util::maximise_over_closed_interval_foc(
      objective,
      [&](double psi_stem, bool* ok) -> double {
        if constexpr (benefit_link<K>() != BenefitLink::Identity) {
          if (use_energy_balance_) {
            if (ok != nullptr) *ok = false;
            return 0.0;
          }
        }
        return dprofit_dpsi_stem<K>(psi_stem, psi_soil, ok);
      },
      psi_soil, psi_crit, basin_scan_cells(), collar_root_tol,
      static_cast<size_t>(ci_niter), &profit_opt);
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

  double beta_ = (theta_ - theta_w_)/(theta_fc_ - theta_w_);

  if(atm_vpd == 0){
     medlyn_model_gs_ = g0;
  } else{
     medlyn_model_gs_ = g0 + 1.6*(1 + (g1*beta_)/sqrt(atm_vpd_))*(assim_colimited_/(ca_*(1/umol_per_mol_to_Pa_)));
  }
  return medlyn_model_gs_;
}

// Supply==demand residual for the Medlyn solver, as a function of ci (x, Pa).
// It is the difference between the Medlyn optimal stomatal conductance and the
// diffusion-implied conductance, *multiplied through by (ca_ - x)* so it stays
// finite across the whole [gamma*, ca_] bracket -- the raw gs difference has a
// 1/(ca_-x) singularity at the upper end. The root (zero crossing) is identical
// to that of the raw difference for x < ca_:
//   gs_medlyn*(ca_-x) - gs_coupled*(ca_-x),  where
//   gs_coupled*(ca_-x) = assim * (atm_kpa_*kPa_to_Pa) * 1.6 / 1e6.
inline double Leaf::medlyn_stom_cond_minus_coupled_stom_cond(double x) {
  const double assim_colimited_x_ = assim_colimited(x);
  medlyn_model_gs_ = medlyn_model_gs(assim_colimited_x_);
  return medlyn_model_gs_ * (ca_ - x)
         - assim_colimited_x_ * (atm_kpa_ * kPa_to_Pa) * 1.6 / 1e6;
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
  return;
}

inline void Leaf::solve_medlyn_ci_analytical(){

  ci_ = ca_ * (g1/(g1 + sqrt(atm_vpd_)));
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  return;
}

} // namespace phylloptim

#endif

