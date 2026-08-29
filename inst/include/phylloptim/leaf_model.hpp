// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_MODEL_HPP_
#define PHYLLOPTIM_LEAF_MODEL_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/inputs.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/uniroot.hpp>
#include <phylloptim/optimize.hpp>
#include <phylloptim/quadrature.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/single_potential.hpp>
#include <phylloptim/vulnerability.hpp>

#include <odelia/implicit_node.hpp>
#include <odelia/interpolator.hpp>
#include <odelia/tangent.hpp>

#include <algorithm>
#include <array>
#include <cmath>
// cstdio/cstdlib for the environment-gated scan of a non-monotone marginal in
// maximise_profit_over_collar, and nothing else in this header.
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace phylloptim {

// The scalars a derivative here runs on, and the two things done to one: seed a
// direction, read the derivative it produced. Every derivative in this package
// is a forward tangent -- there is no tape and the build sets none of the flags
// one would need -- and a curvature is a tangent of a tangent.
//
// Taken from odelia, whose forward-mode header names no tape, so this package
// stays visibly without one while the AD library is named in one place for the
// whole family.
using tangent = odelia::ode::tangent_scalar<double>;
using odelia::ode::derivative_along;
using odelia::ode::seed_direction;

// The parameters profit answers for, at one scalar. WHERE the leaf is operating
// is not among them: the collar, the stem potential and the intercellular CO2
// are what the solve found, and these are what it found them for. They arrive as
// arguments, so nothing can hold a stale one and nothing has to copy this to
// override one.
//
// The two _25 traits enter as the temperature-adjusted values the kernels take,
// so a caller seeds them in the direction their ratio gives and the chain from
// the trait is one factor rather than a pass of its own.
template <typename T>
struct ProfitInputs {
  using value_type = T;
  T vcmax, transport_jmax, quantum_yield, curv_elec, curv_colim, ppfd,
      respiration;
  T kmax, stem_b, stem_c, beta2, cost_scale;

  // Every field, in one order, so a gradient over them cannot be matched to them
  // by hand. The const form carries the list and the mutable one takes the
  // constness back off, so a field is named in one place.
  std::vector<const T*> field_ptrs() const {
    return {&vcmax,  &transport_jmax, &quantum_yield, &curv_elec,
            &curv_colim, &ppfd,       &respiration,   &kmax,
            &stem_b, &stem_c,         &beta2,         &cost_scale};
  }
  std::vector<T*> field_ptrs() {
    std::vector<T*> out;
    for (const T* q : std::as_const(*this).field_ptrs()) {
      out.push_back(const_cast<T*>(q));
    }
    return out;
  }

  // The same inputs with the derivatives taken off, which is what an implicit
  // node asks for to take a slope with the parameters held still.
  template <class U>
  ProfitInputs<U> rebind_from() const {
    using odelia::util::to_passive;
    return {U(to_passive(vcmax)),         U(to_passive(transport_jmax)),
            U(to_passive(quantum_yield)), U(to_passive(curv_elec)),
            U(to_passive(curv_colim)),    U(to_passive(ppfd)),
            U(to_passive(respiration)),   U(to_passive(kmax)),
            U(to_passive(stem_b)),        U(to_passive(stem_c)),
            U(to_passive(beta2)),         U(to_passive(cost_scale))};
  }
};

// `field_ptrs` and `rebind_from` above both list the fields by hand, and a field
// added without a line in either gets no gradient row -- a zero that reads as an
// answer rather than as an omission. The compiler already counts them; this is
// the two agreeing.
static_assert(sizeof(ProfitInputs<double>) == 12 * sizeof(double),
              "ProfitInputs gained or lost a field, and the two lists inside it "
              "are written by hand.");

// Everything the leaf answers for, at one scalar. Two halves because the supply
// is per layer and the rest is not, and one struct because a caller holding two
// would have to keep them in step.
//
// `psi_crit` and `root_psi_crit` sit here rather than in either half: profit
// does not read them at all, and the bounds are the only thing that does.
template <typename T>
struct LeafInputs {
  using value_type = T;
  ProfitInputs<T> profit;
  SupplyValues<T> supply;
  T psi_crit{}, root_psi_crit{};

  template <class U>
  LeafInputs<U> rebind_from() const {
    using odelia::util::to_passive;
    LeafInputs<U> out;
    out.profit = profit.template rebind_from<U>();
    out.supply = supply.template rebind_from<U>();
    out.psi_crit = U(to_passive(psi_crit));
    out.root_psi_crit = U(to_passive(root_psi_crit));
    return out;
  }

  // Every entry a caller can seed, in one order, so a gradient over them cannot
  // be matched to them by hand.
  std::vector<const T*> field_ptrs() const {
    std::vector<const T*> out = profit.field_ptrs();
    for (const T* q : supply.field_ptrs()) out.push_back(q);
    out.push_back(&psi_crit);
    out.push_back(&root_psi_crit);
    return out;
  }
  std::vector<T*> field_ptrs() {
    std::vector<T*> out;
    for (const T* q : std::as_const(*this).field_ptrs()) {
      out.push_back(const_cast<T*>(q));
    }
    return out;
  }
};

class Leaf {
public:
  //anonymous Leaf function as in canopy.h
  Leaf();
  
  Leaf(double vcmax_25, 
       double stem_c, 
       double stem_b, 
       double psi_crit,
       double root_c,
       double root_b,
       double root_psi_crit,
       double beta2, 
       double jmax_25, 
       double a, 
       double curv_fact_elec_trans, 
       double curv_fact_colim,
       double GSS_tol_abs,
       double vulnerability_curve_ncontrol,
       double ci_abs_tol,
       double ci_niter,
      double cost_scale_TF24);

  odelia::interpolator::hermite_interpolator<double, 5> transpiration_from_psi;
  odelia::interpolator::hermite_interpolator<double, 5> psi_from_transpiration;

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
  // of a gradient in stem_b (PLAN 11f).
  //
  // ⚠️ There is NO SUCH IDENTITY FOR stem_c, so a move in stem_c rebuilds. The
  // interpolant carries the closed-form derivative G'(psi) = exp(-(psi/stem_b)^stem_c)
  // as the slope at each knot, so the derivative a reader gets is that closed form,
  // not a value fit's inference of it -- the 3e-4 the two once disagreed by was the
  // fit's error, and it is gone.
  //
  // Two rules keep this from becoming a stale-state bug of the kind hazard 8
  // records: ONLY the four stem_curve_* accessors may read the splines, so one
  // scale factor is sufficient; and `setup_transpiration()` resets it, with
  // `set_traits()` forcing a rebuild while it is displaced, so no route out of
  // the perturbed state leaves it set.
  double stem_b_spline_ = 0.0;

  // The soil -> root-collar water supply (issue #2). Everything the leaf needs
  // from the soil enters through this object: `uptake` and its derivative. It is
  // held by value -- Leaf must stay copyable, because plant's
  // make_strategy_ptr(TF24_Strategy) takes the strategy by value and TF24_Strategy
  // holds a Leaf member.
  //
  // Public because plant's RcppR6 bindings reach the moved fields by name; they
  // now spell them `roots_.psi_soil_` etc. (see PLAN 7b-iii stage 4).
  MultiLayerRoots roots_;

  // The alternative supply path (issue #2 stage 2/3). Both alternatives are held
  // as members and selected by `supply_kind_` -- measured free, where
  // std::variant costs +1.0%; see PLAN 7b-iii stage 2 for the numbers and for why
  // a predictable branch in front of an already-out-of-line call disappears into
  // it. SinglePotential is four doubles plus a one-element vector, so carrying it
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
  // absent. PLAN 7b-iii flagged this as the footgun to design around before
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

  // ⚠️ THIS TAKES NO RESISTANCE, and that is the point of the change that
  // introduced this comment. The soil-to-collar resistance is a per-call DRIVER on
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
  // A layer's root carbon, off the leaf's own network. The single path has none:
  // it carries a series resistance and no root architecture.
  double supply_root_carbon(int layer) const {
    switch (supply_kind_) {
      case SupplyKind::MultiLayer: return roots_.root_carbon(layer);
      default:                     return util::na_value;
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
  double stem_b;
  double psi_crit;  // derived from b and c
  double beta2;
  double jmax_25;
  double a;
  double curv_fact_elec_trans; // unitless - obtained from Smith and Keenan (2020)
  double curv_fact_colim;
  // The single-layer optimisers (optimise_psi_stem_TF / _Sperry) only. They are
  // off the production path and keep brent_fmin because their argmax feeds no
  // gradient. This does not set how well the reported operating point is
  // determined -- collar_root_tol does.
  double GSS_tol_abs;
  // Below this width prepare_collar_solve stops optimising and substitutes the
  // interval's midpoint. That is a decision about which kind of operating point
  // the leaf reports, not a tolerance on a search, so it is a separate number
  // from GSS_tol_abs however similar the two look: a caller loosening a search
  // must not thereby widen the set of points it declines to optimise.
  double collar_interval_min_width;
  // Knots in each vulnerability curve, and the ONE place the count is chosen: both
  // constructors and set_traits' rebuild read it, so a leaf cannot be built on one
  // curve and rebuilt onto another. The default constructor used to hardcode 100
  // beside this member, which is exactly how that happened.
  //
  // ⚠️ IT IS NO LONGER SET BY THE GRID'S OWN ERROR, and what moved is the read
  // rather than the count. Every channel of both tables is a closed form -- G, its
  // slope f, and f's own slope -- so both are read as quintics: h^6 in the value
  // and h^5 in the slope. What the rows follow is the SLOPE, which measured h^3
  // against an h^4 value read, so the slope column is the one that decides.
  // Measured against the closed form, worst over every span:
  //
  //   knots   cubic value   quint value   cubic slope   quint slope   quint KB
  //      25     1.899e-06     3.532e-09     5.252e-05     1.099e-07        0.6
  //      50     1.201e-07     5.555e-11     6.636e-06     3.455e-09        1.2
  //     100     7.506e-09     8.717e-13     8.296e-07     1.084e-10        2.4
  //     200     4.695e-10     1.421e-14     1.038e-07     3.417e-12        4.7
  //     400     2.934e-11     2.220e-15     1.297e-08     4.827e-13        9.4
  //     800     1.834e-12     2.220e-15     1.622e-09     1.251e-12       18.8
  //    1600     1.155e-13     2.220e-15     2.033e-10     2.144e-12       37.5
  //
  // 1600 was the count a CUBIC read needed to bring the rows to the differenced
  // reference's own floor of 2e-09 to 4e-09, and it reached it on the slope
  // column's 2.033e-10. A quintic passes that at 100 knots and sits 60x below it at
  // 200, so by the old rule this would now be 200.
  //
  // ⚠️ IT IS NOT, AND WHAT SETS IT IS THE ROWS RATHER THAN EITHER READ. Reducing it
  // moves the knot SPACING, and two row checks follow the spacing and not G's own
  // error: a wet pin's assembled point row and shade death's differenced rows both
  // pass at 1600, one fails at 800, and two fail at 400 and below -- monotonically,
  // and while the slope column above is sixty times inside its budget the whole
  // way. So the rows are limited by the grid's geometry rather than by how well a
  // span reproduces G, and refining the read does not buy a coarser grid.
  //
  // ⚠️ AND A SECOND THING MOVES WITH THE SPACING, which is why the counts above are
  // read from the row checks and not from the whole suite. At zero flux
  // transpiration_to_psi_stem round-trips G through its own inverse, so the wet end
  // of the collar bracket reads as feasible or not according to the SIGN of the two
  // tables' disagreement -- and the suite's failure count then runs 0, 1, 9, 8, 8
  // over the counts above rather than monotonically. stem_curve_integral_inverse
  // records what fixing that takes and why it is not this accessor's to fix.
  //
  // ⚠️ AND MORE KNOTS ARE NOT SAFER EITHER. The quintic's slope error bottoms out at
  // 400 and then gets WORSE -- 4.8e-13, 1.3e-12, 2.1e-12 -- because the divided
  // differences that placement the top three coefficients lose to roundoff as the span
  // shrinks. A count chosen by refining until it stops moving would overshoot.
  //
  // The solve does not pay for any of it: a read is O(1) on a uniform grid.
  //
  // ⚠️ AND IT IS ONE NAME ACROSS THE PACKAGE BOUNDARY, not one per package. This
  // was three defaults -- Leaf's, phylloptim's R control, and plant's Control --
  // and plant's was the one every stand ran on, so two rounds of refining this
  // member reached nothing plant does. plant's Control now reads this constant;
  // phylloptim's own R defaults restate it, because an R default cannot read it.
  static constexpr double ncontrol_default = 1600.0;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
  double cost_scale_TF24;

  // Tolerance on the collar potential for the profit-maximising root-find (PLAN
  // 11a). Hard-coded rather than a control field, for the same reason
  // find_root_psi's 1e-4 and psi_stem_to_ci's 1e-7 are: it is a property of the
  // solve, not a knob. 1e-12 sits ~8 orders below the 1e-4 at which this package
  // calls a difference real, and the cost of going there from GSS_tol_abs's 1e-3
  // is a handful of evaluations, because TOMS748 is superlinear.
  static constexpr double collar_root_tol = 1e-12;

  // How far dry of the returned root the marginal is probed, to find out which way
  // it crosses there. Deliberately NOT collar_root_tol: the root is located to
  // 1e-12, so a probe that close reads the solver's own noise rather than the
  // function's slope. 1e-6 relative is the scale gradient.hpp's differenced
  // curvature already works at, and it is a probe rather than a difference -- only
  // the SIGN is read, so nothing here needs the step to be optimal.
  static constexpr double collar_probe_frac = 1e-6;


  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
  std::vector<double> soil_consumption_;
  double profit_;
  double psi_stem;
  double lambda_;
  double lambda_analytical_;
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
  double PPFD_;
  double atm_vpd_;
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

  
  double opt_psi_stem_;
  double opt_ci_;
  double count;
  double E_up_;

  // dpsi_stem/dpsi at the collar the marginal profit was last evaluated at. The
  // transport's whole response to the collar, formed there and previously
  // discarded; the condition's gradient is written in it.
  double dpsistem_dpsi_ = util::na_value;

  // The marginal profit's two coefficients at that same collar: its response to
  // the stem potential, and its response to the collar at a HELD stem potential.
  // The condition is the first times dpsistem_dpsi_ plus the second, so a
  // consumer differentiating the condition differentiates these.
  double dprofit_dpsistem_ = util::na_value;
  double dprofit_dpsi_held_stem_ = util::na_value;

  // The rest of the first-order picture there, kept because `condition_slope`
  // would otherwise re-solve the concentration to recover it and because these
  // are the quantities its second derivatives multiply.
  double ci_at_collar_ = util::na_value;
  double assim_slope_ = util::na_value;         // dA/dci
  double dci_dpsistem_ = util::na_value;        // dci/dpsi_stem
  double dci_dpsi_held_stem_ = util::na_value;  // dci/dpsi at a held stem

  // The condition's own value at the collar the solve returned, and whether that
  // evaluation admitted a derivative at all.
  //
  // ⚠️ THE SOLVE CLOSES ON THIS, AND THE BLOCK ABOVE IS WHY. A search's last probe
  // is its answer only at an interior point; at a pin it evaluates both ends and
  // returns one, so a reader taking the block afterwards reads the OTHER end --
  // dpsi_stem/dpsi came back 7.67 where the point's value is 1.46. Evaluating the
  // condition at the collar being returned costs one marginal profit and leaves
  // every coefficient describing the point, so a reader needs no evaluation of its
  // own and the leaf it reads can be const.
  //
  // The two are one channel and one report of whether it exists: the value is the
  // model's own sentinel zero where no condition defines the point, and the flag is
  // what separates that from stationarity.
  double collar_resid_ = util::na_value;
  bool collar_resid_placed_ = false;

  // --- Medlyn stomatal-conductance model (from develop #450) ------------------
  // Standalone, R-callable alternative to the root-collar profit optimisation
  // (solve_medlyn_ci_*); NOT used by the TF24 compute path, which optimises
  // psi_stem directly. g0/g1 default to the published values and are exposed as
  // settable fields. beta_ (soil-moisture stress) uses theta_/theta_w_/theta_fc_,
  // set from the default soil-moisture members in set_physiology.
  double g0 = 0.022;        // residual stomatal conductance (umol m^-2 s^-1)
  double g1 = 2.57;         // sensitivity to vpd (kPa^0.5)
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
  // A memo, so it is mutable: the value it returns is a function of its two
  // arguments alone and the cache is bit-identical to recomputing. Without this
  // every reader that prices a flux is non-const for a stored copy of a number
  // it just derived.
  mutable bool   transpiration_cached_ = false;
  mutable double transpiration_cache_psi_stem_ = 0.0;
  mutable double transpiration_cache_psi_upstream_ = 0.0;
  mutable double transpiration_cache_value_ = 0.0;

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
  static constexpr int photo_temp_key_size = 19;
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
  //     set_physiology's TEMPERATURE CACHE. That cache is keyed on (leaf_temp_,
  //     atm_o2_kpa_) and on nothing else, so calling set_physiology again after a
  //     bare trait write at the SAME temperature takes the cache hit and never
  //     recomputes them. The obvious repair -- "change the trait, then set the
  //     drivers again" -- therefore does not work, which is what makes a settable
  //     field actively dangerous here rather than merely untidy.
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
  // The traits, in the order gradient::par_table numbers them. One argument, so
  // that order is not also spelled by fourteen parameter positions here and
  // fourteen argument positions at every caller.
  void set_traits(const double* theta);
  void set_traits(const std::vector<double>& theta) {
    odelia::util::check_length(theta.size(), std::size_t(gradient::n_traits));
    set_traits(theta.data());
  }

  // The fourteen as separate values, for a caller that holds them that way: the R
  // binding and the tests. Each is PLACED BY THE NAME OF ITS INDEX, so this is the
  // one spelling of the order that a reader can check against par_table without
  // counting positions -- and the only one left.
  void set_traits(double vcmax_25, double stem_c, double stem_b, double psi_crit,
                  double root_c, double root_b, double root_psi_crit,
                  double beta2, double jmax_25, double a,
                  double curv_fact_elec_trans, double curv_fact_colim,
                  double cost_scale_TF24, double R_d_25) {
    double theta[gradient::n_traits];
    theta[gradient::par_vcmax_25] = vcmax_25;
    theta[gradient::par_stem_c] = stem_c;
    theta[gradient::par_stem_b] = stem_b;
    theta[gradient::par_psi_crit] = psi_crit;
    theta[gradient::par_root_c] = root_c;
    theta[gradient::par_root_b] = root_b;
    theta[gradient::par_root_psi_crit] = root_psi_crit;
    theta[gradient::par_beta2] = beta2;
    theta[gradient::par_jmax_25] = jmax_25;
    theta[gradient::par_a] = a;
    theta[gradient::par_curv_fact_elec_trans] = curv_fact_elec_trans;
    theta[gradient::par_curv_fact_colim] = curv_fact_colim;
    theta[gradient::par_cost_scale_TF24] = cost_scale_TF24;
    theta[gradient::par_R_d_25] = R_d_25;
    set_traits(theta);
  }

  // The #25 boundary: the four potentials that must be positive magnitudes. One
  // copy, called from both the constructor and set_traits -- the alternative is
  // two copies that agree until one of them is edited.
  static void check_psi_magnitudes(double psi_crit, double stem_b, double root_b,
                                   double root_psi_crit);

  void setup_transpiration(double resolution);

  // A built vulnerability curve, kept against the pair that determines it.
  //
  // The curve is a pure function of (b, c, resolution) and of nothing else --
  // read cumulative_vulnerability_integral, which takes exactly those three --
  // so the key is everything the value depends on rather than a subset of it,
  // which is the one condition a cache here has to meet.
  //
  // It exists because a trait gradient rebuilds the SAME few pairs over and
  // over: the curve belongs to the traits, and the traits do not depend on which
  // plant is being differentiated, so a caller driving one trait across a stand
  // was paying for a hundred rebuilds of five distinct curves. A perturbation
  // loop touches five pairs per curve -- the base and two sides of each of its
  // two traits -- so a handful of entries holds all of them.
  struct StemCurveCache {
    double b = 0.0, c = 0.0, resolution = 0.0;
    odelia::interpolator::hermite_interpolator<double, 5> from_psi, to_psi;
  };
  // Bounded: past the pairs a perturbation loop visits an entry is never read
  // again, so growing the store without limit would be a leak rather than a hit.
  // The root curve keeps its own beside its splines, in roots_.
  static constexpr std::size_t curve_cache_size = 32;
  // One store per thread, shared by every Leaf on it, because the key determines
  // the value completely: the splines are a pure function of (b, c, resolution)
  // and hold no state of the plant being solved. Per-object stores start empty,
  // so constructing a leaf to overwrite it -- which rebinding to another scalar
  // does, per cohort per stage -- rebuilt both curves every time.
  static std::shared_ptr<std::vector<StemCurveCache>> stem_curve_store() {
    static thread_local std::shared_ptr<std::vector<StemCurveCache>> store =
        std::make_shared<std::vector<StemCurveCache>>();
    return store;
  }
  std::shared_ptr<std::vector<StemCurveCache>> stem_curve_cache_ =
      stem_curve_store();

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
  // d(transpiration)/d(psi) at one end of the column. The flux is the curve
  // integrated BETWEEN the two potentials, so this is the maximum conductance
  // times the curve there -- and the curve, not the table's slope: the
  // tabulation carries G's value because G has no cheap closed form, where f
  // does, so every derivative of G is taken from f exactly.
  double transport_slope(double psi) const {
    return leaf_specific_conductance_max_ * proportion_of_conductivity(psi);
  }
  double stem_curve_integral_inverse(double w, const char* caller = nullptr) const;

  // dG/d(stem_b) at fixed psi, from the same homogeneity and with NO rebuild.
  // G is homogeneous of degree one in (psi, stem_b), so Euler's theorem gives
  //
  //   psi * G'(psi) + stem_b * dG/dstem_b  ==  G(psi)
  //
  // identically, and rearranging it is the whole derivation. Cheap enough to be
  // free beside the two reads it is built from.
  //
  // stem_c has NO counterpart, for the reason perturb_stem_b refuses it: it
  // reshapes the curve rather than scaling it, so its row needs the rebuild.
  double stem_curve_integral_dstem_b(double psi,
                                     const char* caller = nullptr) const;

  // The ROOT curve's three, forwarded so the homogeneity they rest on can be
  // refereed from outside C++ exactly as the stem's is. Hazard 1: these are a
  // different curve in different parameters, and reading one for the other is the
  // mistake that has already reached a manuscript.
  double root_curve_integral(double psi) const {
    return roots_.root_vuln_integral_at(psi);
  }
  double root_curve_integral_deriv(double psi) const {
    return roots_.root_vuln_integral_deriv_at(psi);
  }
  double root_curve_integral_droot_b(double psi) const {
    return roots_.root_vuln_integral_droot_b(psi);
  }

  // Domain-guarded read behind the two accessors above. The stem curve is the
  // only interpolator in this file built with extrapolation DISABLED (the root
  // vulnerability pair clamps instead -- see setup_root_vulnerability), so it is
  // the only one a lookup can throw on. odelia's message names the point and the
  // domain but cannot name WHICH spline, because it does not know: there are two
  // here, they are inverses of each other, and they carry different units, so
  // "u = 7.5 beyond the upper end" is ambiguous in exactly the way that matters.
  // Nor can it name the caller -- the same spline is read from four places, and
  // localising plant#576 came down to which.
  static double eval_stem_curve(
      const odelia::interpolator::hermite_interpolator<double, 5>& spline, double u,
      double scale, const char* spline_name, const char* arg_name,
      const char* caller);

  // Move stem_b WITHOUT rebuilding the stem vulnerability spline, by the
  // homogeneity identity above. stem_c is not accepted: it has no such identity.
  //
  // ⚠️ FOR DERIVATIVE WORK ONLY. It leaves the splines describing a different
  // stem_b, which is sound only because every read goes through the four
  // accessors, and it deliberately does NOT clear the solved operating point --
  // re-placing the physiology is exactly the cost being avoided -- so whatever
  // reads the outputs afterwards must write them first. `set_traits()` is the
  // way back, and forces a rebuild.
  void perturb_stem_b(double stem_b_new);
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
                                               std::vector<double>& y_integral,
                                               std::vector<double>& y_conductivity,
                                               std::vector<double>& y_conductivity_slope) {
    cumulative_vulnerability_integral(weibull_b, weibull_c, resolution, x,
                                      y_integral, y_conductivity,
                                      y_conductivity_slope);
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
  // into them by name after crown integration (PLAN 7b-ii trap 1).
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
  void find_root_collar_psi();
  // Shared setup for the root-collar solve: builds the soil-side caches, handles
  // every feasibility early-exit (shutdown / assim<0 / collapsed interval) by
  // setting the final operating point itself, and otherwise returns the feasible
  // collar-potential interval [bound_a, bound_b] (positive magnitudes). Returns
  // false when the operating point is already fully determined (caller is done),
  // true when there is a real interval to choose a collar potential within.
  bool prepare_collar_solve(double& bound_a, double& bound_b);
  // Evaluate the leaf at a *given* root-collar potential (positive magnitude)
  // rather than optimising it: reuses prepare_collar_solve, clamps the target to
  // the feasible interval, and evaluates there (no golden-section search). Leaves
  // exactly the same outputs as find_root_collar_psi and returns profit_. Used by
  // TF24f's gradient-ascent acclimation (#525).
  double evaluate_root_collar_psi(double target_opt_root_psi);

  // The other consumer of a given collar potential, and it wants the opposite
  // treatment. evaluate_root_collar_psi CLAMPS, which is right for a variant
  // whose collar is a tracked state being nudged: a target outside the feasible
  // interval means the state has drifted out and the projection back in is the
  // model. It is wrong for a partial derivative taken at a FROZEN collar, where
  // the clamp silently replaces the collar the caller asked about with a
  // different one, and near a bound it does so on essentially every arm.
  //
  // So this one does not clamp and does not project. A collar outside the
  // perturbed interval comes back `feasible = false` and IS NOT EVALUATED --
  // below the wet bound the profit algebra runs on a negative conductance and
  // returns a plausible number, so evaluating anyway would answer the question
  // with a number from a different model. The caller refuses the row instead.
  //
  // Not const: it re-runs prepare_collar_solve, which is what re-establishes
  // feasibility at the perturbed state, and that writes the soil-side caches and
  // the classification.
  struct FixedCollarEval {
    double profit;
    std::vector<double> uptake;   // per layer, at the collar as given
    bool feasible;
  };
  FixedCollarEval profit_at_fixed_collar(double collar);

  // A feasibility bound's own derivative. Both bounds are roots of residuals the
  // leaf already evaluates, so the implicit function theorem gives the row
  // without differentiating the search that found it:
  //
  //   dB/du = -(dR/du) / (dR/dx)
  //
  // WHICH bound decides the residual, and they are different functions:
  //   Wet             R0(x) = E_up(x, psi)                    -- no stem terms
  //   DryRootCrit     R(x)  = E_up(x, psi) - kappa*[G(psi_crit) - G(x)]
  //   DryRootPsiCrit  the bound IS a registered constant       -- the row is +1
  //
  // Both steepnesses are here, and for a while neither was: they reshape their
  // vulnerability curve rather than scaling it, so they have no homogeneity
  // identity and their rows were taken by rebuilding the grid and differencing
  // it. What replaced that is the incomplete gamma's shape series, which comes
  // out of the same loop as the integral's value.
  enum class WhichBound { Wet, DryRootCrit, DryRootPsiCrit };
  struct BoundRow {
    std::vector<double> d_dpsi_soil;      // per layer
    std::vector<double> d_droot_carbon;   // per layer
    double d_dkappa = 0.0;
    double d_dpsi_crit = 0.0;             // the STEM's
    double d_droot_psi_crit = 0.0;        // the ROOT's, and only the dry arm has it
    double d_dstem_b = 0.0;
    // The root curve's two parameters. They enter BOTH bounds by the same route
    // -- the layer mean-conductivity integral inside total uptake -- because the
    // stem half of the dry residual does not read either.
    double d_droot_b = 0.0;
    double d_droot_c = 0.0;
    // The stem curve's steepness, and only the dry arm has it: the wet bound is
    // total uptake, which no stem property enters.
    double d_dstem_c = 0.0;
    // dR/dx: the theorem's denominator, and the guard. For the dry arm it is a
    // sum of two strictly positive terms so it cannot change sign, but it can
    // approach zero in deep drought -- so a consumer guards on the amplification
    // it produces rather than on its sign, and gets it back here to do that
    // without recomputing.
    double residual_slope = 0.0;
    double bound = 0.0;                   // where the row was taken
    bool finite = false;
  };



  // The same, flattened for the R boundary: [feasible, profit, uptake_1 ...].
  // Flat for operating_point_values' reason, and with the same obligation --
  // the position of each value here IS the interface.
  std::vector<double> profit_at_fixed_collar_values(double collar) {
    const FixedCollarEval e = profit_at_fixed_collar(collar);
    std::vector<double> out;
    out.reserve(2 + e.uptake.size());
    out.push_back(e.feasible ? 1.0 : 0.0);
    out.push_back(e.profit);
    out.insert(out.end(), e.uptake.begin(), e.uptake.end());
    return out;
  }
  // Evaluate profit at a given root-collar potential (positive magnitude)
  // *assuming prepare_collar_solve has already run this step* (soil-side caches
  // built, feasible interval [bound_a, bound_b] known). Clamps the target into
  // the interval and sets the operating point (opt_psi_stem_, opt_root_psi_,
  // profit_), returning profit_. This is the post-prepare body of
  // evaluate_root_collar_psi, factored out so the centred finite-difference leaf
  // solve can share one prepare_collar_solve across its three profit evals (#530).
  double profit_at_collar_psi(double target_opt_root_psi,
                              double bound_a, double bound_b);
  // Exact d(profit)/d(opt_root_psi) at a given root-collar potential (positive
  // magnitude), for TF24f's acclimation tracking (#525/#527). Combines
  // forward-mode tangent for the analytic photosynthesis/cost algebra, the
  // implicit-function theorem at the psi_stem_to_ci root-find, and analytic
  // spline derivatives (the stem curve's slope) for the smooth transport. Replaces
  // the noisy finite-difference gradient. Places the soil-side caches itself, so a
  // solve need not have run first.
  //
  // ⚠️ **The 0.0 returned on the shut-down / reversed-gradient exits is a
  // SENTINEL, not a stationary point**, and the distinction only became load
  // bearing when PLAN 11a proposed root-finding on `dprofit == 0`. It matters
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
  // Post-prepare body of dprofit_droot_collar_psi, with the same `feasible`
  // contract. Assumes the supply path's per-solve caches are already placed, so
  // the collar solve can share ONE supply_begin_solve across all ~10 of its
  // gradient evaluations instead of re-placing per call -- the same saving #530
  // made for the finite-difference path, and it matters here because
  // begin_solve() is a spline evaluation per soil layer.
  double dprofit_at_collar_psi(double opt_root_psi, bool* feasible = nullptr);

  // ---- the marginal profit, assembled once and read at any scalar ------------
  //
  // The operating point's coordinates, each carrying its OWN derivative in the
  // collar. A structure rather than five parameters because the five must be
  // seeded consistently or the assembly below is a derivative of nothing, and a
  // struct is where that requirement can be stated.
  template <class T>
  struct CollarPoint {
    T p;              // the collar itself; dp/dp = 1
    T sigma;          // the stem potential;  dsigma/dp = V
    T ci;             // intercellular CO2;   dci/dp
    T dEup_dp;        // the soil->collar conductance; its own slope is d2E/dp2
    T transpiration;  // the stem's flux; its collar slope IS dEup_dp, by V's
                      // definition -- the identity dprofit_at_collar_psi states
  };

  // M, the marginal profit this class roots, from those coordinates.
  //
  // ⚠️ ONE ASSEMBLY, TWO SCALARS. `dprofit_at_collar_psi` is this at double and
  // `marginal_collar_slope` is this at a tangent, so the function the solve roots
  // and the slope the implicit function theorem divides by cannot drift apart.
  // They did: plant's curvature came from differentiating a DIFFERENT assembly
  // (profit_at) twice, whose second-order content was two hand-written Taylor
  // coefficients, and on a century stand that put a POSITIVE curvature at an
  // interior maximum -- +34.4 where a difference of this function gives -9.63.
  //
  // A' and C' come from the slope primitives written beside their own kernels, so
  // every factor here is a value and the whole assembly is arithmetic at T. That is
  // what makes dM/dp a FIRST derivative of this function: there is no order above T
  // anywhere in the marginal.

  // What the assembly forms on the way. V and dci_dpsi are the collar responses
  // the coordinates above have to be SEEDED with, so they are returned rather than
  // recomputed by a caller -- the seeds and the assembly then come from one place
  // and cannot disagree.
  template <class T>
  struct MarginalParts {
    T marginal, V, dci_dpsi;
  };

  template <class T>
  MarginalParts<T> marginal_assembled(const CollarPoint<T>& at,
                                      const ProfitInputs<T>& in) const {
    // ⚠️ FROM THE INPUTS, NOT THE MEMBERS. At double the two are the same numbers,
    // but a member read carries no row -- so an assembly that reached for members
    // would answer correctly and differentiate to a silent zero in every parameter
    // it touched. The atmospheric constants below stay members because nothing
    // differentiates them.
    const T kmax = in.kmax;
    const T f_p = proportion_of_conductivity_kernel<T>(at.p, in.stem_b, in.stem_c);
    const T f_sigma =
        proportion_of_conductivity_kernel<T>(at.sigma, in.stem_b, in.stem_c);

    // dpsi_stem/dpsi, from the same closed form the double path uses. At a tangent
    // its own slope falls out of this arithmetic, so no quotient rule is written.
    const T V = (at.dEup_dp / kmax + f_p) / f_sigma;

    // A' and C' read off the primitives that carry them, at T. Every factor of the
    // marginal is now a value, so nothing here is differentiated to be assembled.
    const T J0 = electron_transport_kernel<T>(in.ppfd, in.quantum_yield,
                                              in.curv_elec, in.transport_jmax);
    const T A_prime =
        assim_colimited_slope_kernel<T>(at.ci, in.vcmax, J0, in.curv_colim);
    const T C_prime = hydraulic_cost_TF_slope_kernel<T>(
        at.sigma, in.stem_b, in.stem_c, in.beta2, in.cost_scale);

    const T gc_const =
        T(atm_kpa_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio);
    const T inv_atm = T(1.0 / (atm_kpa_ * kPa_to_Pa));
    const T gc = gc_const * at.transpiration;
    // The conductance partials are the transport slope at each end, which is the
    // closed form kmax * f -- already formed above, so read rather than recomputed.
    const T dgc_dpsistem = gc_const * kmax * f_sigma;
    const T dgc_dpsi = -gc_const * kmax * f_p;

    const T g_ci = A_prime * umol_to_mol + gc * inv_atm;
    const T ca_minus_ci = T(ca_) - at.ci;
    const T dci_dpsistem = (dgc_dpsistem * ca_minus_ci * inv_atm) / g_ci;
    const T dci_dpsi_expl = (dgc_dpsi * ca_minus_ci * inv_atm) / g_ci;
    const T dci_dpsi = dci_dpsistem * V + dci_dpsi_expl;

    return MarginalParts<T>{A_prime * dci_dpsi - C_prime * V, V, dci_dpsi};
  }

  // dM/dp at the operating point the solve placed -- the curvature the interior
  // derivation divides by, taken as a FIRST derivative of the marginal rather than
  // a second of the objective.
  //
  // ⚠️ ONE DEFINITION PER FACTOR AND BOTH ORDERS FROM IT. Every leaf quantity is
  // seeded with its own closed-form slope, so each pair is a Taylor series of a
  // single function. That is the whole finding: the lift this replaces took its
  // value and first derivative from a TABULATION and its second from the TRUE
  // function, and the two disagree by the interpolation error -- which is harmless
  // at first order and is what flipped the sign of a curvature where the collar
  // came within 5.6e-08 of a soil layer's potential.
  //
  // The flux's second order is d2E_from_soil_dpsi_collar2: a closed form that
  // existed, was tested against a difference, and had no production caller.
  double marginal_collar_slope(const ProfitInputs<double>& in) const;


  // What the soil state reaches, and it reaches all of it through total uptake:
  // at a frozen collar the stem potential is the transport read of
  // E_up/kappa + G(p), and the concentration, the conductance, the cost and
  // profit all sit downstream of that. So a soil potential's or a layer carbon's
  // row is one of these times its own supply derivative, and the two supply
  // derivatives it needs are closed form.
  //
  // The condition takes two of them, because it reads the state through the stem
  // potential AND through the stem potential's collar response:
  //
  //   dR/du = dcondition * dE_up/du  +  dprofit * d2E_up/dp du
  //
  // ⚠️ WRITTEN IN (psi_stem, dpsi_stem/dpsi), NOT IN (E_up, dE_up/dpsi). The two
  // pairs span the same directions, so no rank test separates them and a
  // coefficient derived in the second reproduces every supply direction while
  // being short by the collar's direct route into the stomatal conductance.
  // Written this way that route is inside condition_slope, where it is one term
  // of a number that can be checked against a difference of the condition.
  //
  // dcondition is the only second-order field and is left missing where
  // condition_slope refuses; the rest describe the branch wherever the last
  // marginal-profit evaluation recorded one.
  struct UptakeRows {
    double dassim = util::na_value;
    double dstom_cond = util::na_value;
    double dpsistem = util::na_value;
    // dProfit/dE_up at the held collar. It is ALSO dR/d(dE_up/dpsi): uptake's
    // collar slope reaches the condition only through dpsi_stem/dpsi, and the
    // condition's response to that is the stem's own marginal profit, so both
    // are that times dpsistem.
    double dprofit = util::na_value;
    double dcondition = util::na_value;
  };

  // Every output's response to the collar potential, at the traits and the
  // drivers the solve was given. Differencing the outputs across p* measures
  // this; here it is read off what the marginal-profit evaluation recorded.
  //
  // ⚠️ A DIFFERENCE CANNOT BE CENTRED AT A PINNED POINT, which is where the
  // constrained rows need this most: p* sits a millionth of the bracket from its
  // bound, so one arm is outside the feasible interval and no shrinking brings it
  // back. Read this way the response exists wherever the point does.
  //
  // The collar's own response is not a field. It is 1 by identity, and a stored
  // copy of it is a number that could disagree with what it is.
  //
  // `duptake` is in kg, as every supply derivative in this class is, where the
  // reported per-layer output is the soil's own consumption in mol.
  struct CollarRows {
    double dassim = util::na_value;
    double dstom_cond = util::na_value;
    double dpsistem = util::na_value;
    double dprofit = util::na_value;
    std::vector<double> duptake;
  };


  // The three parameters of the TRANSPORT, whose rows were the last taken by
  // rebuilding a curve and differencing it.
  //
  // At a frozen collar the flux through the stem IS the flux the soil supplies --
  // kappa (G(sigma) - G(p)) = E_up(p) -- and the soil reads no stem property, so
  // all three of these move the stem potential and leave the flux where it is.
  // Everything downstream of the flux therefore holds: measured over a step, the
  // conductance, the concentration and assimilation hold to 5e-10 and the uptake
  // to exactly zero, while sigma and profit move. So the whole held row is the
  // hydraulic cost, which is the one thing past the flux that reads sigma.
  //
  // The condition's row is the cost too, for the reason at the definition.
  //
  // False with the energy-balance gate on, for the reason the other readers refuse
  // there.
  enum class TransportTrait { Conductance, Position, Steepness };
  struct TransportTraitRows {
    double dpsistem = util::na_value;
    double dprofit = util::na_value;
    double dmarginal = util::na_value;
  };


  // The three photosynthesis traits, whose rows come off two second-order passes
  // rather than six re-solves.
  //
  // At a frozen collar these move nothing but assimilation: the stem potential,
  // the stomatal conductance and the hydraulic cost are all fixed, and the only
  // thing that responds is the intercellular CO2 the residual places. So with
  // g_ci = A'*umol_to_mol + gc*inv_atm and K = (ca - ci)*inv_atm/g_ci,
  //
  //   dci/dtheta  = -(dA/dtheta) * umol_to_mol / g_ci
  //   dprofit     =  (dA/dtheta) * gc * inv_atm / g_ci
  //   dmarginal   =  (dA'/dtheta) * D * K + A' * D * dK/dtheta
  //
  // with D the collar's route into gc, which carries no trait. dA'/dtheta needs
  // the mixed second partial and A'', and both come from the assimilation kernel
  // seeded twice -- the model's own algebra, not a hand-kept second copy.
  //
  // Six traits, four passes, because the family shares its intermediates.
  // `a`, curv_fact_elec_trans and jmax_25 reach assimilation only through the
  // electron transport, so ONE pass in that direction serves all three and they
  // enter as dJ/dtheta. curv_fact_colim reaches it only through the
  // colimitation and vcmax_25 only through the rubisco-limited rate, one pass
  // each. R_d_25 needs no pass at all: dark respiration is subtracted from the
  // colimitation, so dA/dR_d_ is exactly -1 and d2A/dR_d_ dci exactly 0.
  //
  // The three that carry a temperature-derived scalar -- vcmax_, jmax_, R_d_ --
  // reach it from their _25 trait through a factor that does not depend on the
  // trait, since both Arrhenius forms are the reference value times a function
  // of temperature. So the chain is the derived value over the trait, exactly.
  //
  // Their frozen-collar uptake rows are exactly zero, for the cost traits'
  // reason: at a fixed collar a carbon-side trait moves no water.
  //
  // On the compensation-point branch gross assimilation is identically zero, so
  // none of them reaches profit at all and every row here is zero -- except
  // R_d_25, which is what net assimilation is reduced BY, and so still carries
  // its own row there.
  struct PhotoTraitRows {
    double dprofit_da, dprofit_dcurv_elec, dprofit_dcurv_colim;
    double dprofit_dvcmax_25, dprofit_djmax_25, dprofit_dR_d_25;
    double dmarginal_da, dmarginal_dcurv_elec, dmarginal_dcurv_colim;
    double dmarginal_dvcmax_25, dmarginal_djmax_25, dmarginal_dR_d_25;
    // Radiation is not a trait, and it belongs here anyway: it reaches
    // assimilation through the electron transport and through nothing else, which
    // is the family `a` and the transport curvature are in, so its rows come off
    // the pass already taken in that direction for the cost of one more seed.
    double dprofit_dPPFD, dmarginal_dPPFD;
  };
  // The profit-maximising collar potential within [bound_a, bound_b], by a
  // safeguarded root-find on dprofit == 0 (PLAN 11a). Returns a bound when the
  // optimum is pinned to it, which is the case on 42 of the 240 feasible
  // golden-grid rows and so is a branch that has to be written rather than a
  // corner. Falls back to golden section if either endpoint has no usable
  // gradient -- see the definition for why that fallback should never fire.
  double maximise_profit_over_collar(double bound_a, double bound_b);
  // Analytic d(E_up_)/d(collar suction) -- a CONDUCTANCE, positive by
  // construction now that both sides are magnitudes (#25). Thin forwarder to
  // roots_.duptake_dpsi; see there for the derivation and for the NaN-at-a-kink
  // contract. Used only on the TF24f acclimation gradient path, not the base
  // TF24 value path.
  double dE_from_soil_dpsi_collar(double T_collar,
                                  const std::vector<double>& psi_soil) const {
    require_suction_vector(psi_soil, "dE_from_soil_dpsi_collar");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        return roots_.duptake_dpsi(T_collar, psi_soil);
      default:
        return single_.duptake_dpsi(T_collar, psi_soil);
    }
  }
  // The collar derivative of that conductance. It is what turns the condition's
  // slope in the collar into a statement rather than a difference: the transport
  // response V is (E_up'(p)/kappa + f(p))/f(sigma), so differentiating it in the
  // collar reads this and the vulnerability curve's own slope and nothing else.
  // Same NaN-at-a-kink contract as the conductance.
  double d2E_from_soil_dpsi_collar2(double T_collar,
                                    const std::vector<double>& psi_soil) const {
    require_suction_vector(psi_soil, "d2E_from_soil_dpsi_collar2");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        return roots_.d2uptake_dpsi2(T_collar, psi_soil);
      default:
        return single_.d2uptake_dpsi2(T_collar, psi_soil);
    }
  }
  // Per-layer d(E_i)/d(psi_soil[i]), the soil counterpart of the conductance
  // above. Diagonal, so one entry per layer is the whole of it. Same NaN-at-a-
  // kink contract; see MultiLayerRoots::duptake_dpsi_soil.
  void dE_from_soil_dpsi_soil(double T_collar,
                              const std::vector<double>& psi_soil,
                              std::vector<double>& out) {
    require_suction_vector(psi_soil, "dE_from_soil_dpsi_soil");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        roots_.duptake_dpsi_soil(T_collar, psi_soil, out);
        break;
      default:
        single_.duptake_dpsi_soil(T_collar, psi_soil, out);
        break;
    }
  }
  // d2(E_i)/d(collar suction) d(psi_soil[i]), the collar derivative of the
  // conductance above. Diagonal for its reason. The single path's flux is linear
  // in the difference over a constant resistance, so there it is exactly zero.
  void d2E_from_soil_dpsi_collar_dpsi_soil(double T_collar,
                                           const std::vector<double>& psi_soil,
                                           std::vector<double>& out) {
    require_suction_vector(psi_soil, "d2E_from_soil_dpsi_collar_dpsi_soil");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        roots_.d2uptake_dpsi_dpsi_soil(T_collar, psi_soil, out);
        break;
      default:
        out.assign(psi_soil.size(), 0.0);
        break;
    }
  }
  // d(E_i)/d(root carbon in layer a) and its collar derivative, both lower
  // triangular. Root carbon is the multi-layer architecture's input; the single
  // path has no carbon profile to move and both blocks come back empty.
  void dE_from_soil_droot_carbon(double T_collar,
                                 const std::vector<double>& psi_soil,
                                 std::vector<std::vector<double>>& dE_drc,
                                 std::vector<std::vector<double>>& dD_drc) {
    require_suction_vector(psi_soil, "dE_from_soil_droot_carbon");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        roots_.duptake_droot_carbon(T_collar, psi_soil, dE_drc, dD_drc);
        break;
      default:
        dE_drc.clear();
        dD_drc.clear();
        break;
    }
  }
  using SupplyCurveTrait = MultiLayerRoots::CurveTrait;
  // d(E_i)/d(a root curve parameter) and its collar derivative, per layer and in
  // kg. Both are closed form: the curve reaches the supply through one integral,
  // and that integral's trait derivatives are Euler's identity and the incomplete
  // gamma's shape series. The single path has no root curve to move.
  void dE_from_soil_droot_curve(double T_collar,
                                const std::vector<double>& psi_soil,
                                SupplyCurveTrait trait,
                                std::vector<double>& dE,
                                std::vector<double>& d2E) {
    require_suction_vector(psi_soil, "dE_from_soil_droot_curve");
    switch (supply_kind_) {
      case SupplyKind::MultiLayer:
        roots_.duptake_droot_curve_by_layer(T_collar, psi_soil, trait, dE);
        roots_.d2uptake_dpsi_droot_curve(T_collar, psi_soil, trait, d2E);
        break;
      default:
        dE.assign(psi_soil.size(), 0.0);
        d2E.assign(psi_soil.size(), 0.0);
        break;
    }
  }
  // Shut-down operating point used by the find_root_collar_psi early-exits: stem
  // held at psi_crit (no transpiration), paying only respiration + hydraulic
  // cost. Only opt_root_psi_ differs between the cases, so it is the argument
  // (a positive magnitude, like every other psi here).
  void set_shutdown_state(double root_collar);
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit);
  double find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil);
  double E_column(double x, const std::vector<double>& psi_soil, double psi_leaf);
  double E_column_zero(double x, const std::vector<double>& psi_soil);
  
  double arrh_curve(double Ea, double ref_value, double leaf_temp) const;
  double peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const;

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

  // How often each clamp in this model held. Sites are in clamp_sites.hpp, and
  // every one of them pairs its clamp with a matching derivative -- so this is a
  // measure of DISTANCE from a defect rather than of one.
  //
  // It exists because two of the sites bind within a few hundredths of a soil
  // moisture unit of the wilt point, and because a consumer DIFFERENCING through a
  // clamp gets a finite number with no thermal or hydraulic response in it and
  // nothing else to notice by: the derivative kill that keeps the analytic rows
  // honest is invisible to a difference.
  clamp_counter clamps;

  // How many times the operating point has been SOLVED on this object.
  //
  // A cost rather than a defect, which is why it is counted here and not beside
  // the clamps. A row layer's whole claim is that it makes ONE solve however many
  // inputs are asked for, and nothing outside could see that -- so it was asserted
  // by the proxy that no input re-solves, which is weaker: a route could solve
  // twice for reasons that have nothing to do with an input.
  //
  // Behind a pointer for clamp_counter's reason and no other: a consumer holds this
  // model by value and copies it per unit, so a plain member would take the count
  // down with the copy.
  std::shared_ptr<std::size_t> collar_solves =
      std::make_shared<std::size_t>(0);


  // This leaf's own sites plus the supply model's, which are one list. Summed on
  // read rather than shared on construction, so rebuilding the root network cannot
  // silently detach the tally.
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

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;
  // Scalar-generic core of the above. See the `_kernel` block below the assim
  // declarations for why these exist.
  template <typename T> T proportion_of_conductivity_kernel(T psi) const;
  // The same curve with its two parameters as arguments rather than read off
  // the object, so forward mode can be seeded on THEM. The one-argument form
  // above reaches only d/dpsi, which is why a trait row through this curve had
  // no analytic route.
  template <typename T>
  T proportion_of_conductivity_kernel(T psi, T b, T c) const;
  // d/dpsi of the curve above, in one definition beside it, so a caller wanting the
  // slope differentiates nothing.
  template <typename T>
  T proportion_of_conductivity_slope_kernel(T psi, T b, T c) const;

  // supply-side transpiration for a given water potential gradient between leaves and soil,
  // references setup_transpiraiton for values (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration(double psi_stem, double psi_upstream) const;
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
  // The same expression with the light as an argument, so it can be
  // differentiated in PPFD by the same forward-mode route the assimilation and
  // cost kernels already use. Cold: reached once per set_physiology and once per
  // light row, never from inside the collar solve.
  template <typename T> T electron_transport_kernel(T ppfd) const;
  // The same expression again with the three traits that reach it on the
  // scalar, which is the only route into them: each appears here and nowhere
  // else. jmax_ arrives from jmax_25 through the temperature block, which is
  // exactly linear in it.
  template <typename T> T electron_transport_kernel(T ppfd, T quantum_yield,
                                                    T curvature, T jmax) const;
  double assim_electron_limited(double ci_);
  double assim_colimited(double ci_);

  // --- the scalar-generic cores -------------------------------------------
  //
  // These are the SINGLE definition of the photosynthesis and cost algebra. The
  // five `double` entry points above and beside them are `T = double`
  // instantiations, and `dprofit_droot_collar_psi` differentiates the same code
  // with `T = tangent`. So the forward model and its
  // derivative cannot disagree: there is one body, not two.
  //
  // ⚠️ **They replace a hand-maintained `namespace detail`, which HAD ALREADY
  // DRIFTED** -- the deleted `assim_colimited_ad` associated the electron-limited
  // term left-to-right where `assim_electron_limited` divides the bracket first,
  // and used `s*s` where `assim_colimited` uses `pow(s, 2)`. Both are pure
  // reassociation, so the two functions were mathematically identical and
  // numerically not: the tangent derivative was the derivative of a *slightly
  // different function* than the model evaluated. Harmless while the gradient only
  // set TF24f's acclimation rate; load-bearing once PLAN 11a made the collar solve
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
  // `hydraulic_cost_`; its kernel must not, or the tangent pass would write model state
  // while probing.
  //
  // Three of them take a trait on the scalar as well as the variable, and the
  // two-argument forms forward to those with the member. `a` and
  // `curv_fact_elec_trans` reach assimilation through the electron transport and
  // through nothing else, and `curv_fact_colim` through the colimitation and
  // through nothing else -- each appears in exactly one expression in this file
  // -- so carrying them here is what lets photo_trait_rows differentiate the
  // family instead of moving a member and re-solving.
  template <typename T> T assim_rubisco_limited_kernel(T ci) const;
  template <typename T> T assim_rubisco_limited_kernel(T ci, T vcmax) const;
  template <typename T> T assim_electron_limited_kernel(T ci) const;
  template <typename T> T assim_electron_limited_kernel(T ci, T transport) const;
  template <typename T> T assim_colimited_kernel(T ci) const;
  template <typename T> T assim_colimited_kernel(T ci, T vcmax, T transport,
                                                 T curvature,
                                                 T respiration) const;
  template <typename T> T colimit_kernel(T assim_rubisco_limited_,
                                         T assim_electron_limited_) const;
  template <typename T> T colimit_kernel(T assim_rubisco_limited_,
                                         T assim_electron_limited_,
                                         T curvature, T respiration) const;

  // Each kernel's slope in ci, beside the kernel it belongs to and decomposed the
  // same way, so every pair is checkable against the two lines above it.
  // Respiration is absent from the colimited slope because it is subtracted, and a
  // constant in ci contributes nothing.
  template <typename T> T assim_rubisco_limited_slope_kernel(T ci, T vcmax) const;
  template <typename T> T assim_electron_limited_slope_kernel(T ci,
                                                              T transport) const;
  template <typename T> T assim_colimited_slope_kernel(T ci, T vcmax, T transport,
                                                       T curvature) const;
  // The colimitation's slope along a direction the two limited rates move in, so
  // the chain rule is written once here rather than at each caller.
  template <typename T> T colimit_slope_kernel(T assim_rubisco_limited_,
                                               T assim_electron_limited_,
                                               T d_rubisco, T d_electron,
                                               T curvature) const;
  // G(psi) at one scalar: the value the table holds, carrying the three slopes
  // the closed form gives exactly.
  //
  // ⚠️ THE VALUE IS THE TABLE'S AND EVERY SLOPE IS THE CURVE'S. The solve ran on
  // the table, so a value from the closed form would move the operating point;
  // the table's own slope is a fit of the curve and differs from it by 1.15e-11,
  // so a derivative taken from the table is a derivative of the fit.
  template <typename T>
  T stem_integral_at(const T& psi, const ProfitInputs<T>& in) const;

  // The water the soil delivers at a collar, at one scalar, against the soil
  // state the solve cached. Both supply paths answer the same question, so a
  // caller differentiating the collar asks one thing rather than branching.
  //
  // This is the ONLY route the collar reaches the carbon side by: the stomatal
  // conductance is proportional to the flux and the stem potential is what
  // carries it.
  template <typename T>
  T E_from_soil_at(const T& collar, const SupplyAt<T>& supply,
                   std::vector<T>& per_layer) const {
    per_layer.assign(static_cast<std::size_t>(supply_n_layers()), T(0.0));
    T out = T(0.0);
    if (supply_kind_ == SupplyKind::MultiLayer) {
      roots_.uptake(collar, supply, per_layer, out);
      return out;
    }
    if constexpr (std::is_same_v<T, double>) {
      single_.uptake_at(collar, supply.psi_soil, per_layer, out);
      return out;
    } else {
      // Refused rather than answered with zeros. The single-potential path's
      // resistance is a parameter of its own and no differentiated caller is on
      // that path; a silent zero row here would read as an insensitivity.
      util::stop("E_from_soil_at: the single-potential supply path has no "
                 "differentiated form; the multi-layer one is what a gradient "
                 "runs on");
    }
  }

  // This leaf's own supply state, as the view the quadrature takes. No copies.
  SupplyAt<double> held_supply() const { return roots_.held_supply(); }

  // Everything a consumer asks the leaf for, at one scalar.
  //
  // ⚠️ `point` USED TO BE HERE, and it went because the distinction it carried was
  // one nothing acted on. It reported a collar the theorem could not place, and it
  // was filled by ONE of collar_at's five arms -- the interior one -- while the
  // other four threw. Its single reader turned it into `refuse(...)`, which is
  // exactly what the catch around the throws does. Two mechanisms, one outcome, and
  // a field on every LeafOutputs to carry the difference between them.
  template <class S>
  struct LeafOutputs {
    S profit{};
    std::vector<S> uptake;
  };

  // What the soil delivers at the operating point the solve placed, recorded ONCE.
  //
  // The state reaches everything above this through the two waists and the draws
  // and through nothing else, so a consumer holding one of these needs no supply at
  // all. Taken at the PASSIVE collar, which is where every residual below is
  // evaluated anyway: implicit_value asks its residual for the value at the point
  // the solve found and the rows in the inputs, so the collar it reads carries no
  // direction.
  //
  // ⚠️ THE COLLAR IS PART OF THE VALUE. A flux and a collar that disagree are a
  // caller's mistake that no arithmetic here can see, so the two travel together
  // and `at` is checked against the operating point before anything reads it.
  template <class S>
  struct SupplyDraw {
    double at = 0.0;        // the collar this was taken at
    S flux;                 // E_up there
    S slope;                // dE_up/dp there
    std::vector<S> uptake;  // the per-layer draws
  };

  // The only way to make one, so the draw records the collar it was taken at and
  // the three values cannot come from different ones. The collar is passive here
  // whatever the caller hands in: this is the supply at a POINT, and where the
  // point itself moves is the graft's business rather than the draw's.
  template <class S>
  SupplyDraw<S> supply_draw_at(const S& collar,
                               const SupplyValues<S>& supply) const {
    SupplyDraw<S> d;
    d.at = odelia::util::to_passive(collar);
    // ⚠️ NO DRAW AT A SHUTDOWN, and not because it would be zero. The collar is
    // held at the stem's critical potential there, which is past where the uptake
    // model answers at all -- the integral is exactly zero and E_up is not a
    // number, so asking costs a stop rather than a wrong row. Every consumer reads
    // zero at that kind anyway. ShadeDeath is NOT this case: it sits on the wet
    // bound, where E_up is zero in value and its rows are the bound's own theorem.
    if (operating_point_kind_ == OperatingPointKind::HydraulicShutdown) {
      d.uptake.assign(static_cast<std::size_t>(supply_n_layers()), S(0.0));
      return d;
    }
    const S held = S(d.at);
    d.flux = E_from_soil_at<S>(held, supply.at(), d.uptake);
    d.slope = roots_.template duptake_dpsi_at<S>(held, supply.at());
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

  // The collar this solve left, at whatever scalar the caller wants, from
  // whatever pins it. The kind chooses which condition closes the system and
  // nothing else.
  template <typename S>
  S collar_at(const LeafInputs<S>& in, const SupplyDraw<S>& draw) const;

  // The outputs at a collar, whatever placed it. Underneath the kind switch this
  // is one composition: the two residuals, the kernels and the quadrature.
  template <typename S>
  LeafOutputs<S> outputs_at(const S& collar, const LeafInputs<S>& in,
                            const SupplyDraw<S>& draw) const;

  // Profit at a collar the solve already placed, differentiable in every input
  // that reaches it AND in the collar itself. `sigma_star` and `ci_star` are the
  // two the solve found there, so nothing here searches: they enter the two
  // residuals that define them and come back carrying the implicit function
  // theorem's quotient.
  //
  // The collar is an argument at the working scalar, which is what makes the
  // marginal profit a seeded pass over this rather than a second derivation:
  // seed the collar and read the direction. At an interior point profit's own
  // row needs none of it -- the collar is where the marginal is zero, which is
  // the envelope theorem -- but the marginal itself is what places the collar,
  // and every water output reads that.
  // The operating point's other coordinates at scalar S, from the anchors the solve
  // left: the flux the soil delivers at this collar, and the stem potential and the
  // intercellular CO2 the two residuals define there. profit_at composes the profit
  // from these and marginal_at composes the marginal from them, so ONE construction
  // serves both and they cannot disagree about where the point is.
  template <class S>
  struct CollarCoords {
    S sigma, ci;
  };

  // The flux is an ARGUMENT rather than a local, which is the waist: the soil
  // reaches these two residuals through it and through nothing else, so a caller
  // holding one recording of it can serve every consumer here.
  template <typename S>
  CollarCoords<S> collar_coords_at(double sigma_star, double ci_star,
                                   const S& collar, const S& flux,
                                   const ProfitInputs<S>& in) const;

  template <typename S>
  S profit_at(double sigma_star, double ci_star, const S& collar, const S& flux,
              const ProfitInputs<S>& in) const;

  // M at scalar S: the marginal profit as a function of the ACTIVE INPUTS, with the
  // collar wherever the caller puts it. This is the residual the interior collar is
  // closed by -- `implicit_value_reported(p*, dM/dp, marginal_at)` -- so that its
  // parameter rows are TAPED from the condition rather than handed over as numbers
  // from a second-order pass nothing referees.
  template <typename S>
  S marginal_at(const S& collar, const SupplyDraw<S>& draw,
                const ProfitInputs<S>& in) const;

  // The collar where a pinned point sits, at one scalar, from the condition that
  // pins it. Three arms because the dry end is a `min` of two limits that are
  // different functions of the inputs -- and none of them needs a second
  // derivative, so unlike the interior condition this composes on the caller's
  // tape with nothing handed over.
  //
  //   Wet           E_up(x) = 0                       uptake vanishes
  //   DryRootCrit   kmax*(G(psi_crit) - G(x)) - E_up(x) = 0
  //                                                   T1 with the stem at its limit
  //   DryRootPsiCrit  x = root_psi_crit                a registered constant
  template <typename S>
  S bound_at(WhichBound which, double bound_x, const LeafInputs<S>& in,
             const SupplyDraw<S>& draw) const;

  template <typename T> T hydraulic_cost_TF_kernel(T psi_stem) const;
  template <typename T>
  T hydraulic_cost_TF_kernel(T psi_stem, T b, T c, T beta, T scale) const;
  // dC/dpsi_stem, one line over the curve's own slope.
  template <typename T>
  T hydraulic_cost_TF_slope_kernel(T psi_stem, T b, T c, T beta, T scale) const;

  // dC/d(psi_stem, stem_b, stem_c, beta2, cost_scale) at a stem potential.
  //
  // A zero-flux operating point pays only respiration and this cost -- profit is
  // -R_d - C(psi_crit), reading neither soil nor light -- so these ARE its trait
  // rows, and every environment row there is exactly zero.
  struct HydraulicCostRow {
    double d_dpsi_stem = 0.0;
    double d_dstem_b = 0.0;
    double d_dstem_c = 0.0;
    double d_dbeta2 = 0.0;
    double d_dcost_scale = 0.0;
    bool finite = false;
  };
  double assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream);
  // The energy-balance correction to the above, zero when the gate is off. Kept
  // out of line so that adding it cannot change FMA contraction in the inlined
  // gate-off path; the derivation and the two sign checks are at the definition.
  double dprofit_energy_balance_term(double ci, double gc, double g_ci,
                                     double inv_atm, double gc_const,
                                     double dgc_dpsistem, double dgc_dpsi,
                                     double dpsistem_dpsi, double dT_dE,
                                     double Tleaf);
  BoundRow bound_row(WhichBound which);
  double psi_stem_to_ci(double psi_stem, double psi_upstream);
  // `ci_known` is the intercellular CO2 an earlier solve found at this same
  // state. Every no-flux branch fixes ci itself, so it is read at the one place a
  // root-find would otherwise run, and NA there means solve for it.
  void set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream,
                                           double ci_known = util::na_value);


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
// NOTE the name clash to be careful of: the member `lambda_` is an *input* to
// profit_psi_stem_Sperry (Sperry's prescribed marginal water cost, never set by
// set_physiology). The lambda here is an *emergent output*. See PLAN.md 10a.

  // lambda at an arbitrary stem water potential (positive magnitude, MPa), in
  // umol CO2 (kg H2O)^-1. Analytic, from the TF24 cost function:
  //
  //   C(psi)   = cost_scale_TF24 * (1 - f)^beta2,   f(psi) = exp(-(psi/stem_b)^stem_c)
  //   dC/dpsi  = cost_scale_TF24 * beta2 * (1-f)^(beta2-1) *
  //              (stem_c/stem_b)(psi/stem_b)^(stem_c-1) * f
  //   E(psi)     = kmax * integral of f,  so  dE/dpsi = kmax * f
  //   lambda     = (dC/dpsi) / (dE/dpsi)
  //
  // The f cancels exactly; it is cancelled here rather than divided out, so the
  // expression stays finite as psi -> psi_crit where f -> 0.
  //
  // Caveat: with beta2 < 1 the (1-f)^(beta2-1) factor diverges as psi -> 0
  // (f -> 1), which is a property of the cost function, not of this code.
  double lambda_TF24(double psi_stem) const;

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

  // The NET marginal value of water arriving from the soil: lambda_multi minus
  // lambda_stem.
  //
  // This, and not lambda_multi, is what prices a change in soil potential. Water
  // a wetter layer supplies still has to be pushed through the stem, and
  // lambda_stem is by its own definition "the marginal cost seen by the stem",
  // so that cost nets off; lambda_multi is the GROSS value and already includes
  // the soil->collar transport the wetter soil is relieving you of. Taking the
  // gross figure overstates the soil rows by 1.201-2.369x over the golden grid.
  //
  // Computed as the product rather than as the subtraction it is defined by.
  // The two are algebraically identical -- lambda_multi = lambda_stem*(1 + k/S)
  // so the difference is lambda_stem*k/S exactly -- and the product form does
  // not cancel two nearby numbers to get a small one. Same reason the soil rows
  // exist at all rather than being differenced out of profit.
  //
  // Returns the NA sentinel wherever marginal_cost_water_multilayer does, and
  // for the same reason: S unavailable at a branch kink, or non-positive.
  double marginal_price_water();


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

  // Every reported output of a solved point, in one call.
  //
  // ⚠️ THIS EXISTS FOR THE R BOUNDARY AND FOR NOTHING ELSE. A C++ caller should
  // read the members; they are public and free. From R each read is a separate
  // call through an RcppR6 active binding at ~1.1 us, so the twelve of them cost
  // ~15 us against a ~3 us solve -- five times the model, to report it (#39).
  // One call is ~1.5 us. That is the whole justification, and it is why this is
  // a flat vector rather than a struct: it crosses the boundary as a numeric
  // vector with no glue.
  //
  // The ORDER is the interface. R's .operating_point_names names these positions
  // and test-surface.R checks the two agree by reading all twelve fields
  // individually and comparing -- so a field inserted here without a
  // corresponding R change fails rather than silently shifting a column.
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
            g1_eff()};
  }

// leaf economics functions
  double hydraulic_cost_Sperry(double psi_stem, double psi_upstream);
  double hydraulic_cost_TF(double psi_stem);

  double profit_psi_stem_Sperry(double psi_stem, double psi_upstream);
  double profit_psi_stem_TF(double psi_stem, double psi_upstream,
                            double ci_known = util::na_value);

// optimiser functions
  void optimise_psi_stem_Sperry();
  void optimise_psi_stem_TF();

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
    // Constrained optimum pinned at the WET end of the feasible interval, just
    // inside root_zero_E (the collar at which uptake is exactly zero). profit is
    // still climbing toward the wet bound, so the bound is the answer. 24 golden
    // points. The gradient here is genuinely non-zero.
    PinnedWet,
    // Constrained optimum pinned at the DRY end, min(root_crit,
    // supply_psi_crit()). 18 golden points. Gradient genuinely non-zero.
    // The dry bound is a min of two limits that are DIFFERENT FUNCTIONS of the
    // inputs, so a consumer forming the bound's row needs to know which won. The
    // stem's continuity root is a search result and its row is dense; the root's
    // own critical potential is a registered constant, and its row is exactly
    // minus the unit vector in its own direction and zero everywhere else.
    PinnedDryRootCrit,
    PinnedDryRootPsiCrit,
    // The feasible interval collapsed to a point (width below
    // collar_interval_min_width), so feasibility DETERMINED the collar
    // potential and nothing was optimised.
    // There is no free variable left to differentiate.
    Determined,
    // Shutdown on water: no collar potential both moves water and stays inside
    // the stem's and the root's critical potentials. The stem holds at psi_crit,
    // transpiration is zero, and the leaf pays respiration plus the hydraulic
    // cost there. 48 of the 288 golden grid points at every temperature, since it is
    // hydraulics rather than heat that forbids transpiration. The water response is zero;
    // the carbon response is not.
    HydraulicShutdown,
    // The feasible interval INVERTED: the dry bound landed wetter than the
    // collar at which uptake is exactly zero, so no operating point both moves
    // water and stays inside the critical potentials. Reported apart from
    // HydraulicShutdown because it is reached by a different branch, and whether
    // it is a dry soil or a parameterisation the model cannot represent is then
    // a question a counter can answer rather than one this tag has to assume.
    InfeasibleBracket,
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
  // Which limit won the dry bound. Recorded where the min is taken, because a
  // min is the one operation that destroys the information a consumer needs
  // afterwards: the two arms are different functions of the inputs.
  enum class DryBoundArm { None, RootCrit, RootPsiCrit };
  DryBoundArm dry_bound_arm() const { return dry_bound_arm_; }

  // How many kinds there are, derived from the last enumerator rather than
  // written out, so a kind added before it is counted without a second edit.
  // NonFiniteGradient must stay last.
  static constexpr std::size_t operating_point_kind_count =
      static_cast<std::size_t>(OperatingPointKind::NonFiniteGradient) + 1;

  OperatingPointKind operating_point_kind() const {
    return operating_point_kind_;
  }
  static const char* operating_point_kind_name(OperatingPointKind kind);
  // The branch this solve took, by name. A fixture that asserts which regime it
  // is in has to read the classification rather than compare the collar against
  // a bound: the two agree until a step-in or a tolerance makes them differ.
  const char* operating_point_kind_name() const {
    return operating_point_kind_name(operating_point_kind_);
  }

  // True where the solve ended with the stomata shut. No flux, no gross
  // assimilation, and a profit that is respiration plus the hydraulic cost
  // alone -- so nothing about such a point is an argmax, and a consumer that
  // takes an envelope step at one is answering about a maximisation that did
  // not happen.
  //
  // ⚠️ THE TWO ARE NOT THE SAME POINT. HydraulicShutdown holds the stem at
  // psi_crit and moves no water at all, so every environment row there is
  // exactly zero. ShadeDeath places BOTH potentials at the collar of zero
  // uptake -- which is the wet bound -- so its profit reads the soil through
  // that bound, and the per-layer consumptions are not zero either: they sum
  // to zero.
  bool zero_flux_operating_point() const {
    return operating_point_kind_ == OperatingPointKind::HydraulicShutdown ||
           operating_point_kind_ == OperatingPointKind::ShadeDeath;
  }

  // The whole point a search for the operating point found, in the form the
  // placement below takes it back in: the collar it returned, the two quantities
  // the collar determines, the branch it was found on, and which limit won the dry
  // bound.
  //
  // ⚠️ THE STEM POTENTIAL AND THE CONCENTRATION ARE HERE BECAUSE RE-DERIVING THEM
  // IS THE COST. The collar alone leaves a placement inverting the stem curve and
  // running psi_stem_to_ci -- a bracketing root-find to 1e-10 -- to reach numbers
  // the search already had. Two doubles per entry against ~9 evaluations of the
  // colimited assimilation.
  //
  // Only Leaf writes one. A caller can carry one and hand it back and cannot make
  // one up, which is what keeps the branch an output of the solve rather than
  // something a caller can disagree with -- the same argument that makes the tag
  // read-only, with the difference that this can be handed back.
  class SolvedPoint {
  public:
    SolvedPoint() = default;
    // False where no search happened: the four branches that exit on feasibility
    // cost nothing to reach again, so nothing is kept for them.
    bool searched() const { return kind != OperatingPointKind::Unsolved; }

  private:
    friend class Leaf;
    double collar = util::na_value;
    double sigma = util::na_value;  // the stem potential the collar determines
    double ci = util::na_value;     // the intercellular CO2 that goes with it
    OperatingPointKind kind = OperatingPointKind::Unsolved;
    DryBoundArm arm = DryBoundArm::None;
  };

  // The point this solve found, or one holding nothing where it exited before
  // searching for one.
  SolvedPoint solved_point() const {
    SolvedPoint out;
    switch (operating_point_kind_) {
      case OperatingPointKind::Interior:
      case OperatingPointKind::PinnedWet:
      case OperatingPointKind::PinnedDryRootCrit:
      case OperatingPointKind::PinnedDryRootPsiCrit:
        out.collar = opt_root_psi_;
        out.sigma = opt_psi_stem_;
        out.ci = ci_;
        out.kind = operating_point_kind_;
        out.arm = dry_bound_arm_;
        break;
      default:
        break;
    }
    return out;
  }

  // Place a point an earlier solve found instead of searching for it. Returns false
  // for a point holding nothing, which is the caller's cue to solve.
  bool place_solved_point(const SolvedPoint& point);

private:
  // Written by every path out of the collar solve, and reset to Unsolved at the
  // top of prepare_collar_solve. Hazard 8 in the developer guide is why: `Leaf`
  // is a value member that plant reuses for every individual in a patch, so a
  // branch that declines to write this would leave the PREVIOUS plant's
  // classification -- a plausible answer about a different plant, which is the
  // worst failure shape available here. Defaulting the reset to Unsolved means a
  // path that forgets reports "unclassified" instead.
  OperatingPointKind operating_point_kind_ = OperatingPointKind::Unsolved;
  DryBoundArm dry_bound_arm_ = DryBoundArm::None;
};

// Human-readable tag. The switch has no default, so a missing name is a -Wswitch
// warning -- an error under the test suite's -Werror=switch, "unknown" elsewhere.
inline const char* Leaf::operating_point_kind_name(OperatingPointKind kind) {
  switch (kind) {
    case OperatingPointKind::Unsolved:          return "unsolved";
    case OperatingPointKind::Interior:          return "interior";
    case OperatingPointKind::PinnedWet:         return "pinned-wet";
    case OperatingPointKind::PinnedDryRootCrit: return "pinned-dry-root-crit";
    case OperatingPointKind::PinnedDryRootPsiCrit:
      return "pinned-dry-root-psi-crit";
    case OperatingPointKind::Determined:        return "determined";
    case OperatingPointKind::HydraulicShutdown: return "hydraulic-shutdown";
    case OperatingPointKind::InfeasibleBracket: return "infeasible-bracket";
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
// templated on their scalar type, so tangent differentiates the model itself and the
// mirror cannot drift because there is no mirror.
inline Leaf::Leaf()
    :
    vcmax_25(96), // umol m^-2 s^-1 
    stem_c(2.680147), //unitless
    stem_b(3.898245), //-MPa
    psi_crit(5.870283), //-MPa
    beta2(1.5), //exponent for effect of hydraulic risk (unitless)
    jmax_25(157.44), // maximum electron transport rate umol m^-2 s^-1
    a(0.30), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosythnthesis equatiom
    GSS_tol_abs(1e-3),
    collar_interval_min_width(1e-3),
    vulnerability_curve_ncontrol(ncontrol_default),
    ci_abs_tol(1e-3),
    ci_niter(1000),
    cost_scale_TF24(7.5) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      // The root traits (root_c/root_b/root_psi_crit) and the two beta_R_*
      // resistance constants keep their defaults in MultiLayerRoots, which owns
      // them. Deliberately not restated here: a second copy of the root Weibull
      // pair is the exact shape of hazard 1 in the developer guide.
      setup_transpiration(vulnerability_curve_ncontrol);
      setup_root_vulnerability(vulnerability_curve_ncontrol);
      setup_clean_leaf();
}

inline Leaf::Leaf(double vcmax_25, double stem_c, double stem_b,
           double psi_crit, // derived from b and c,
           double root_c,
           double root_b,
           double root_psi_crit,
           double beta2, double jmax_25,
           double a, double curv_fact_elec_trans, double curv_fact_colim, 
           double GSS_tol_abs,
           double vulnerability_curve_ncontrol,
           double ci_abs_tol,
           double ci_niter,
           double cost_scale_TF24)
    : vcmax_25(vcmax_25), // umol m^-2 s^-1 
    stem_c(stem_c), //unitless
    stem_b(stem_b), //-MPa
    psi_crit(psi_crit), //-MPa
    beta2(beta2), //exponent for effect of hydraulic risk (unitless)
    jmax_25(jmax_25), // maximum electron transport rate umol m^-2 s^-1
    a(a), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(curv_fact_elec_trans), //curvature factor for the light response curve (unitless)
    curv_fact_colim(curv_fact_colim), //curvature factor for the colimited photosythnthesis equation
    GSS_tol_abs(GSS_tol_abs),
    collar_interval_min_width(GSS_tol_abs),
    vulnerability_curve_ncontrol(vulnerability_curve_ncontrol),
    ci_abs_tol(ci_abs_tol),
    ci_niter(ci_niter),
    cost_scale_TF24(cost_scale_TF24) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      // The single convention, asserted at the one place it enters (#25). Before
      // there was no global statement about psi's sign to assert -- psi_soil_ was
      // >= 0, psi_soil_inverted_ was <= 0, psi_crit was >= 0 and the collar was
      // <= 0 -- which is precisely why the convention had to live in comments.
      // Now it is checkable, so it is checked. set_traits shares this check.
      check_psi_magnitudes(psi_crit, stem_b, root_b, root_psi_crit);

      // The root traits belong to the supply path, so hand them over before its
      // vulnerability curve is built.
      roots_.root_c = root_c;          //unitless
      roots_.root_b = root_b;          // MPa, positive magnitude
      roots_.root_psi_crit = root_psi_crit; // MPa, positive magnitude

      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_root_vulnerability(vulnerability_curve_ncontrol);
      setup_clean_leaf();
}

inline void Leaf::check_psi_magnitudes(double psi_crit, double stem_b,
                                      double root_b, double root_psi_crit) {
  if (!(psi_crit > 0.0)) {
    util::stop("psi_crit must be a positive magnitude in MPa (#25); got " +
               util::to_string(psi_crit));
  }
  if (!(stem_b > 0.0)) {
    util::stop("stem_b must be a positive magnitude in MPa (#25); got " +
               util::to_string(stem_b));
  }
  if (!(root_b > 0.0)) {
    util::stop("root_b must be a positive magnitude in MPa (#25); got " +
               util::to_string(root_b));
  }
  if (!(root_psi_crit > 0.0)) {
    util::stop("root_psi_crit must be a positive magnitude in MPa (#25); got " +
               util::to_string(root_psi_crit));
  }
}

// See the header for why this exists rather than fourteen settable fields.
inline void Leaf::set_traits(const double* theta) {
  // Each trait taken by the name of its index, so a member and the slot it is
  // read from are written down beside each other.
  const double vcmax_25_ = theta[gradient::par_vcmax_25];
  const double stem_c_ = theta[gradient::par_stem_c];
  const double stem_b_ = theta[gradient::par_stem_b];
  const double psi_crit_ = theta[gradient::par_psi_crit];
  const double root_c_ = theta[gradient::par_root_c];
  const double root_b_ = theta[gradient::par_root_b];
  const double root_psi_crit_ = theta[gradient::par_root_psi_crit];
  const double beta2_ = theta[gradient::par_beta2];
  const double jmax_25_ = theta[gradient::par_jmax_25];
  const double a_ = theta[gradient::par_a];
  const double curv_fact_elec_trans_ = theta[gradient::par_curv_fact_elec_trans];
  const double curv_fact_colim_ = theta[gradient::par_curv_fact_colim];
  const double cost_scale_TF24_ = theta[gradient::par_cost_scale_TF24];
  const double R_d_25_ = theta[gradient::par_R_d_25];

  check_psi_magnitudes(psi_crit_, stem_b_, root_b_, root_psi_crit_);

  // Which splines have to be rebuilt, decided BEFORE the assignment. Exact
  // equality is the right test and not a sloppy one: the spline is a pure
  // function of the pair, so an unchanged pair gives a bit-identical spline and
  // rebuilding it is pure cost. A trait perturbed by a relative 1e-08 -- what a
  // gradient loop does -- compares unequal and rebuilds, which is the case that
  // must not be missed.
  // ⚠️ The third clause is what makes set_traits the way back from
  // perturb_stem_b. While the splines are built at a different stem_b, "the
  // parameters did not move" is not a reason to keep them -- it is exactly the
  // case where the equality test would conclude there is nothing to do and leave
  // the object rescaling forever.
  const bool stem_curve_moved = (stem_b_ != stem_b) || (stem_c_ != stem_c) ||
                                (stem_b != stem_b_spline_);
  const bool root_curve_moved =
      (root_b_ != roots_.root_b) || (root_c_ != roots_.root_c);

  vcmax_25 = vcmax_25_;
  stem_c = stem_c_;
  stem_b = stem_b_;
  psi_crit = psi_crit_;
  beta2 = beta2_;
  jmax_25 = jmax_25_;
  a = a_;
  curv_fact_elec_trans = curv_fact_elec_trans_;
  curv_fact_colim = curv_fact_colim_;
  cost_scale_TF24 = cost_scale_TF24_;
  R_d_25 = R_d_25_;

  roots_.root_c = root_c_;
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
  setup_clean_leaf();
}

// set various states and physiology parameters obtained from TF24 to NA to clean leaf object
inline void Leaf::setup_clean_leaf() {
  ci_at_compensation_point_ = false;  // hazard 8: every exit writes its own state
  ci_ = util::na_value; // Pa
  stom_cond_CO2_= util::na_value; //mol Co2 m^-2 s^-1 
  assim_colimited_= util::na_value; // umol C m^-2 s^-1 
  transpiration_= util::na_value; // kg m^-2 s^-1 
  profit_= util::na_value; // umol C m^-2 s^-1 
  lambda_= util::na_value; // umol C m^-2 s^-1 kg^-1 m^2 s^1
  lambda_analytical_= util::na_value; // umol C m^-2 s^-1 kg^-1 m^2 s^1
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
  Tair_= util::na_value; // deg C
  Rn_= util::na_value; // W m^-2
  ra_= util::na_value; // s m^-1
  PPFD_= util::na_value; //umol m^-2 s^-1
  atm_vpd_= util::na_value; //kPa 
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
  // needs the soil profile placed first.
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
      util::stop("find_root_psi(find_root_crit=1) failed: " + std::string(e.what()) +
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
    util::stop("find_root_psi(find_root_crit=0) failed: " + std::string(e.what()) +
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
inline void Leaf::set_shutdown_state(double root_collar) {
  opt_root_psi_ = root_collar;
  opt_psi_stem_ = psi_crit;
  profit_ = -R_d_ - hydraulic_cost_TF(psi_crit);
  // Tagged here rather than at the call sites, so a further reason to shut down
  // cannot arrive without a classification. The sites that take this default are
  // one ecological statement -- the soil is too dry for any collar potential
  // that both moves water and stays inside the critical potentials. The
  // INVERTED-interval site is not that statement and overwrites this tag
  // immediately; it is the only one that does, and it says so there.
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
  std::fill(soil_consumption_.begin(), soil_consumption_.end(), 0.0);
  // No condition defines this collar, so the pair is declared rather than
  // evaluated -- which is also what the evaluation returns here, since the collar
  // is at least as dry as the stem and the marginal profit takes its no-flow exit.
  // Declaring it costs nothing and evaluating it would move the fluxes just
  // written above.
  collar_resid_ = 0.0;
  collar_resid_placed_ = false;
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
inline bool Leaf::prepare_collar_solve(double& bound_a, double& bound_b){

  // Clear the classification FIRST, so that a path which declines to write it
  // reports "unclassified" rather than the previous plant's kind of operating
  // point (hazard 8). Every exit below, and both callers, write it again.
  //
  // The condition at the returned collar goes with it, and for the same reason: a
  // reader after an exit that declined to write it would otherwise get a plausible
  // number from the previous plant.
  collar_resid_ = util::na_value;
  collar_resid_placed_ = false;
  operating_point_kind_ = OperatingPointKind::Unsolved;
  // Same reason, and it needs saying because the arm is written LATER than the
  // classification is: an exit taken before the dry bound is formed would
  // otherwise leave the previous plant's arm beside this plant's kind.
  dry_bound_arm_ = DryBoundArm::None;

  // Hand the supply path the start of a solve: it builds its per-solve caches and
  // reports back the wettest rooted layer -- the SMALLEST suction, and the lower
  // bracket endpoint everything below needs.
  const double wettest_soil_layer = supply_begin_solve();

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to
  // shut down

  if (wettest_soil_layer >= psi_crit){
    set_shutdown_state(psi_crit);
    return false;
  }

if(E_column(psi_crit, supply_psi_soil(), psi_crit) < 0){
      set_shutdown_state(supply_psi_crit());
      return false;
}

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to
  // shut down
double root_crit = find_root_psi(wettest_soil_layer, supply_psi_soil(), 1);

// If root crit would have to be larger than psi crit, also avoid loop as above

    if (root_crit >= psi_crit){
    set_shutdown_state(root_crit);
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

    profit_ = - R_d_ - hydraulic_cost_TF(opt_root_psi_);
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
    // As on the shut-down exits: the stem is tied to the collar here, so the
    // marginal profit takes its no-flow exit and returns the sentinel. The cost's
    // own slope at the placement is what a consumer needs instead, and it is a read.
    collar_resid_ = 0.0;
    collar_resid_placed_ = false;

        if(std::isnan(profit_)){
          util::stop("Error: profit nan");
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
    // std::min returns its SECOND argument only when that is strictly smaller,
    // so this ternary reproduces the tie-break rather than guessing at it.
    dry_bound_arm_ = supply_psi_crit() < root_crit ? DryBoundArm::RootPsiCrit
                                                   : DryBoundArm::RootCrit;
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
      // The outputs are a shutdown's -- no flux, stem at the critical potential
      // -- so they are written by the same call. The KIND is not: this is the
      // interval inverting, which is a different branch from a soil too dry to
      // pay for water, and the two are worth telling apart before either is
      // answered for.
      set_shutdown_state(supply_psi_crit());
      operating_point_kind_ = OperatingPointKind::InfeasibleBracket;
      return false;
    }

    // If no interval exists (single feasible root-collar value), use that
    // point directly as the alternative solution instead of optimising.
    if (std::abs(bound_b - bound_a) <= collar_interval_min_width) {
      const double opt_root_psi = 0.5 * (bound_a + bound_b);
      const double psi_stem_single = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());

      if (!std::isfinite(psi_stem_single)) {
        util::stop("Error: non-finite psi_stem_single in collapsed-root interval; "
                   "opt_root_psi=" + util::to_string(opt_root_psi) +
                   "; bound_a=" + util::to_string(bound_a) +
                   "; bound_b=" + util::to_string(bound_b) +
                   "; root_crit=" + util::to_string(root_crit) +
                   "; root_zero_E=" + util::to_string(root_zero_E) +
                   "; E_up_=" + util::to_string(E_up_));
      }

      opt_psi_stem_ = psi_stem_single;
      profit_ = profit_psi_stem_TF(opt_psi_stem_, opt_root_psi);
      opt_root_psi_ = opt_root_psi;
      // Feasibility DETERMINED this point; no maximisation happened, and there is
      // no free variable left for a derivative to move -- so there is no condition
      // to evaluate at it either.
      operating_point_kind_ = OperatingPointKind::Determined;
      collar_resid_ = 0.0;
      collar_resid_placed_ = false;

      if (!std::isfinite(profit_)) {
        util::stop("Error: non-finite profit in collapsed-root interval; "
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
// dprofit/dpsi == 0 instead of searching profit itself (PLAN 11a). Assumes
// prepare_collar_solve has run, so the soil-side caches are placed and every
// gradient evaluation below can go straight to dprofit_at_collar_psi.
//
// WHY this replaced golden section, in one line each -- PLAN 11a has the numbers:
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
// Safeguarded, not Newton: TOMS748 keeps a bracket throughout, so a non-monotone
// or kinked case degrades instead of diverging. Monotonicity was measured on all
// 240 feasible golden-grid rows, but that is a grid, not a theorem. This also
// makes the collar the third leaf solver on this method -- see uniroot.hpp's
// history note, which is about exactly this judgement call.
inline double Leaf::maximise_profit_over_collar(double bound_a, double bound_b) {
  const double width = bound_b - bound_a;
  const auto dprofit = [&](double psi, bool* feasible = nullptr) {
    return dprofit_at_collar_psi(psi, feasible);
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
          return profit_psi_stem_TF(psi_stem, bound);
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
  // What lo buys, stated honestly: on a wet-pinned row the true argmax sits
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
    const double p_lo = profit_psi_stem_TF(
        find_psi_stem_from_psi_root(lo, supply_psi_soil()), lo);
    const double p_hi = profit_psi_stem_TF(
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
    operating_point_kind_ = OperatingPointKind::PinnedWet;
    return lo;
  }
  if (f_hi >= 0.0) {
    operating_point_kind_ = dry_bound_arm_ == DryBoundArm::RootPsiCrit
        ? OperatingPointKind::PinnedDryRootPsiCrit
        : OperatingPointKind::PinnedDryRootCrit;
    return hi;
  }

  // Interior stationary point. f_lo > 0 > f_hi, so the bracket is valid and
  // uniroot_smooth's throw-on-bad-bracket cannot fire; passing the two endpoint
  // values it already has saves it re-evaluating them.
  //
  // ⚠️ A NOTE THAT WAS WRONG, KEPT BECAUSE THE MEASUREMENT THAT KILLED IT MATTERS.
  // This said the solve returns a MINIMUM of profit on a century stand, because
  // plant's reverse-mode guard reported a POSITIVE profit curvature at an interior
  // point and an interior point is only reached on a bracket the marginal crosses
  // downward. The reasoning was sound and the premise was false.
  //
  // MEASURED, at the offending point (collar 1.7984860121573381) on that stand:
  //
  //     dprofit at x-h = +1.71e-05      dprofit at x+h = -1.75e-05
  //     centred difference             = -9.633
  //     plant's nested-AD curvature    = +34.414
  //
  // The marginal crosses DOWNWARD. The collar is a genuine maximum, the curvature
  // there is negative, and this solve is right. The probe below -- one evaluation
  // dry of every root it returns -- ran over 2,829,445 interior solves on that
  // stand and found the marginal rising exactly zero times.
  //
  // So the defect was in the CURVATURE, not here: it came from differentiating
  // profit_at twice through a divided difference that had already spent its
  // accuracy on the span. Both are gone -- the mean conductivity no longer
  // differences at a short span, and the curvature is marginal_collar_slope(), a
  // FIRST derivative of the marginal this solve roots.
  operating_point_kind_ = OperatingPointKind::Interior;
  const double root = util::uniroot_smooth(dprofit, lo, hi, f_lo, f_hi,
                                           collar_root_tol,
                                           static_cast<size_t>(ci_niter));

  // ⚠️ THE CROSSING DIRECTION IS NOT MEASURABLE BY A PROBE HERE, and a probe was
  // tried. An interior collar is only reached on a bracket the marginal crosses
  // downward, so the slope at a converged root must be non-positive -- but near a
  // collar/soil coincidence the marginal has a sharp feature, and MEASURED on the
  // century stand: a difference at a 1e-6 relative step gives -9.63, at 1e-8 gives
  // +166, at 1e-10 gives -10588. It does not converge. Above the feature the step
  // straddles it; inside it, the marginal's own ~1e-6 evaluation noise divided by
  // the step swamps the signal. So no step size works, and a probe at 1e-6 reported
  // the crossing as downward at all 2,829,445 interior solves on that stand while
  // two independent analytic routes agreed the slope there is +34.41.
  //
  // Any detection here has to be ANALYTIC for that reason. marginal_collar_slope()
  // is dM/dp in closed form and is what plant's guard now tests -- the same number
  // the interior closure divides by, rather than a second one that could disagree
  // with it.
  return root;
}

// The same operating point, placed rather than searched for. Every line below is
// find_root_collar_psi's own closing, in its own order, so at the collar that solve
// returned this places what it placed -- bit for bit, which is the only referee a
// placement can have.
//
// What it does not do is look for the feasible interval: prepare_collar_solve places
// the soil-side caches on its way to the bounds, so this places them itself and stops
// there. That is the whole saving -- two bound root-finds and the search between
// them, against one evaluation of the condition.
inline bool Leaf::place_solved_point(const SolvedPoint& point) {
    if (!point.searched()) {
      return false;
    }
    // Counted like any other: a placement is a solve of the operating point, so the
    // tally stays one per state whichever way the state got there.
    ++(*collar_solves);
    supply_begin_solve();
    operating_point_kind_ = point.kind;
    dry_bound_arm_ = point.arm;

    collar_resid_ = dprofit_at_collar_psi(point.collar, &collar_resid_placed_);

    opt_psi_stem_ = point.sigma;
    opt_root_psi_ = point.collar;
    profit_ = profit_psi_stem_TF(point.sigma, point.collar, point.ci);

    if(!std::isfinite(profit_)){
        util::stop("Error: non-finite profit at a placed operating point; "
             "opt_psi_stem_=" + util::to_string(opt_psi_stem_) +
             "; opt_root_psi_=" + util::to_string(opt_root_psi_) +
             "; E_up_=" + util::to_string(E_up_) +
             "; assim_colimited_=" + util::to_string(assim_colimited_) +
             "; hydraulic_cost_=" + util::to_string(hydraulic_cost_));
    }
    return true;
}

inline void Leaf::find_root_collar_psi(){
    // Counted before the feasibility exits, not after: a shut-down or determined
    // point is a solve of the operating point whether or not it optimises anything.
    ++(*collar_solves);
    double bound_a, bound_b;
    if (!prepare_collar_solve(bound_a, bound_b)) {
      return;
    }
    // root_crit / root_zero_E were consumed inside prepare_collar_solve; recover
    // them for the diagnostic message only if the profit check below fails.

    // Maximise carbon profit over the feasible collar-potential interval by
    // solving its first-order condition, dprofit/dpsi == 0, rather than by
    // searching the objective. PLAN 11a has the reasoning and the measurements.
    const double opt_root_psi = maximise_profit_over_collar(bound_a, bound_b);

    // ⚠️ CLOSE ON THE CONDITION AT THE COLLAR BEING RETURNED, and do it BEFORE the
    // outputs are placed below. The coefficients a consumer reads off this object
    // describe the last marginal profit evaluated, and the search's last probe is
    // the answer only at an interior point: at a pin it evaluates both ends of the
    // interval and returns one, so a reader afterwards described the other end --
    // dpsi_stem/dpsi came back 7.67 where the point's value is 1.46, finite and
    // plausible, at every pinned point.
    //
    // Ordered before the placement rather than after it because this evaluation
    // moves the fluxes: the placement is what the outputs are, so it goes last and
    // nothing has to be saved and restored. The soil caches are already placed by
    // prepare_collar_solve, so the body is called rather than the wrapper.
    collar_resid_ = dprofit_at_collar_psi(opt_root_psi, &collar_resid_placed_);

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());

    opt_root_psi_ = opt_root_psi;
    profit_ = profit_psi_stem_TF(opt_psi_stem_, opt_root_psi);


    if(!std::isfinite(profit_)){
        util::stop("Error: non-finite profit; opt_psi_stem_=" + util::to_string(opt_psi_stem_) +
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
inline double Leaf::evaluate_root_collar_psi(double target_opt_root_psi){
    double bound_a, bound_b;
    if (!prepare_collar_solve(bound_a, bound_b)) {
      // Operating point fully determined by feasibility handling (shutdown /
      // assim<0 / collapsed interval); profit_ is already set.
      return profit_;
    }

    return profit_at_collar_psi(target_opt_root_psi, bound_a, bound_b);
}


inline Leaf::FixedCollarEval Leaf::profit_at_fixed_collar(double collar) {
  FixedCollarEval out;
  out.profit = std::numeric_limits<double>::quiet_NaN();
  out.feasible = false;

  double bound_a, bound_b;
  if (!prepare_collar_solve(bound_a, bound_b)) {
    // Feasibility fixed the operating point on its own, so there is no interval
    // for a collar to be held inside and no frozen-collar partial to take.
    return out;
  }
  if (collar < bound_a || collar > bound_b) {
    return out;
  }

  // find_psi_stem_from_psi_root places the uptake at this collar on its way past,
  // and the profit is taken at the stem potential it returns, so both outputs
  // describe the collar that was asked for.
  opt_root_psi_ = collar;
  opt_psi_stem_ = find_psi_stem_from_psi_root(collar, supply_psi_soil());
  profit_ = profit_psi_stem_TF(opt_psi_stem_, collar);
  E_from_Soil_to_Root_Collar(collar, supply_psi_soil());
  operating_point_kind_ = OperatingPointKind::Prescribed;

  out.profit = profit_;
  out.uptake = soil_consumption_;
  out.feasible = std::isfinite(profit_);
  return out;
}

inline double Leaf::profit_at_collar_psi(double target_opt_root_psi,
                                  double bound_a, double bound_b){
    if (target_opt_root_psi < bound_a || target_opt_root_psi > bound_b) {
      // Counted, because a caller differencing this function across the bound gets
      // the feasibility jump rather than a derivative. The unclamped form for that
      // purpose is profit_at_fixed_collar, which refuses instead.
      clamps.note(CLAMP_COLLAR_POTENTIAL);
    }
    const double opt_root_psi =
        std::min(std::max(target_opt_root_psi, bound_a), bound_b);

    opt_psi_stem_ = find_psi_stem_from_psi_root(opt_root_psi, supply_psi_soil());
    opt_root_psi_ = opt_root_psi;
    profit_ = profit_psi_stem_TF(opt_psi_stem_, opt_root_psi);
    // Not an optimum of any kind: the collar potential came from the caller. On
    // this path a non-zero dprofit is the answer (TF24f's acclimation rate), not
    // a residual to be checked against zero.
    operating_point_kind_ = OperatingPointKind::Prescribed;

    if(!std::isfinite(profit_)){
        util::stop("Error: non-finite profit in evaluate_root_collar_psi; "
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
// tangent; the gc partials use the analytic spline derivative (transpiration_from_psi
// .deriv); dpsi_stem/dpsi by a tight central difference on the smooth transport.
inline double Leaf::dprofit_droot_collar_psi(double opt_root_psi, bool* feasible) {
  // Every transport evaluation below reads the supply path's per-solve caches, so
  // placement them on the current psi_soil_ here rather than depending on whatever the
  // caller's last solve left cached. Keeping this in the wrapper is what lets the
  // collar solve call the body directly and placement them once for the whole solve.
  supply_begin_solve();
  return dprofit_at_collar_psi(opt_root_psi, feasible);
}

// dM/dp at the placed operating point. See the declaration for why this is a FIRST
// derivative of the marginal rather than a second of the objective.
inline double Leaf::marginal_collar_slope(const ProfitInputs<double>& in) const {
  using odelia::ode::derivative_along;
  using odelia::ode::seed_direction;
  using T = odelia::ode::tangent_scalar<double>;

  const double p = opt_root_psi_;
  const double sigma = opt_psi_stem_;
  const std::vector<double>& soil = supply_psi_soil();

  // The two leaf quantities that are not closed-form arithmetic, each with its OWN
  // slope: the soil-to-collar conductance and its tested second derivative, and
  // the stem's flux, whose collar slope IS that conductance by V's definition.
  const double dEup = dE_from_soil_dpsi_collar(p, soil);
  const double d2Eup = d2E_from_soil_dpsi_collar2(p, soil);
  const double transp = transpiration(sigma, p);

  // One pass at double for the collar responses the tangent seeds are, then one at
  // a tangent to read their slope. Same assembly both times.
  const CollarPoint<double> here{p, sigma, ci_, dEup, transp};
  const MarginalParts<double> parts = marginal_assembled<double>(here, in);

  CollarPoint<T> at{T(p), T(sigma), T(ci_), T(dEup), T(transp)};
  seed_direction(at.p, 1.0);
  seed_direction(at.sigma, parts.V);
  seed_direction(at.ci, parts.dci_dpsi);
  seed_direction(at.dEup_dp, d2Eup);
  seed_direction(at.transpiration, dEup);

  return derivative_along(
      marginal_assembled<T>(at, in.template rebind_from<T>()).marginal);
}

inline double Leaf::dprofit_at_collar_psi(double opt_root_psi, bool* feasible) {
  const double psi = opt_root_psi;
  // gstar_Pa used to be precomputed here and threaded into the tangent replicas. The
  // kernels read gamma_ and umol_per_mol_to_Pa_ themselves, which is one fewer
  // place for the two sides to disagree.
  // Infeasible until the two exits below have been passed; see the header for why
  // the 0.0 they return must be distinguishable from a genuine stationary point.
  if (feasible != nullptr) {
    *feasible = false;
  }
  // ⚠️ The recorded picture below describes THIS evaluation, and every exit
  // before it returns without writing it. Cleared here so a reader after one of
  // them refuses on a missing number instead of answering from the last
  // evaluation that got through -- which is a plausible number at a nearby
  // collar, and so invisible.
  dprofit_dpsistem_ = util::na_value;
  dprofit_dpsi_held_stem_ = util::na_value;
  ci_at_collar_ = util::na_value;
  assim_slope_ = util::na_value;
  dci_dpsistem_ = util::na_value;
  dci_dpsi_held_stem_ = util::na_value;

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

  // ⚠️ E1: PLACE THE TEMPERATURE PARAMETERS AT THIS CANDIDATE. Without this the
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
  }

  const double ci = psi_stem_to_ci(psi_stem, psi);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  // dpsi_stem/dpsi. Computed HERE, before the compensation-point branch, so both
  // branches share ONE evaluation.
  //
  // ⚠️ It was briefly a separate member function, and this comment claimed that
  // cost 1.4% because the compiler stopped inlining it. THAT WAS WRONG, and
  // measuring it properly is what showed so: with both binaries interleaved,
  // the factored and inlined versions are indistinguishable. Kept inline anyway
  // because sharing one evaluation between the two branches is less work
  // regardless -- but not for the reason first given.
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
    // The stem carries the flux the soil supplies -- kappa (G(sigma) - G(p)) =
    // E_up(p) -- so the collar response follows from differentiating THAT, with
    // no inverse curve in it and nothing read off a table's slope:
    //
    //   V = (dE_up/dp / kappa + f(p)) / f(sigma)
    //
    // ⚠️ IT WAS P'(x) TIMES THE BRACKET, AND THE DIFFERENCE IS SECOND ORDER.
    // Reading the inverse's slope makes V a value of an interpolant, so a
    // consumer differentiating V differentiates that interpolant's SECOND
    // derivative -- which is a property of the fit and not of the curve, and no
    // supplied first-order data corrects it. Against a difference of the V it
    // belongs to: 5e-09 to 1e-07 this way, 3e-06 to 4e-05 the other. The two V
    // agree to 1e-08, so this is a change of route rather than of answer.
    const double f_stem = proportion_of_conductivity(psi_stem);
    dpsistem_dpsi =
        (dEup_dpsi / leaf_specific_conductance_max_ +
         proportion_of_conductivity(psi)) / f_stem;
  } else {
    // Near a branch kink the analytic conductance returns NaN; fall back to a
    // central difference on the transport, as this path has always done.
    const double h = 1e-6;
    dpsistem_dpsi =
        (find_psi_stem_from_psi_root(psi + h, supply_psi_soil()) -
         find_psi_stem_from_psi_root(psi - h, supply_psi_soil())) / (2.0 * h);
  }
  dpsistem_dpsi_ = dpsistem_dpsi;


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
    tangent ps_ad0 = psi_stem;  seed_direction(ps_ad0, 1.0);
    const double C_prime0 = derivative_along(hydraulic_cost_TF_kernel(ps_ad0));
    double dprofit = -C_prime0 * dpsistem_dpsi;
    if (use_energy_balance_ && dT_dE != 0.0) {
      const double h = 1e-3;
      const double vc0 = vcmax_, jm0 = jmax_, ga0 = gamma_, ko0 = ko_,
                   kc0 = kc_, rd0 = R_d_, km0 = km_, J0 = electron_transport_;
      update_temperature_dependent_params(Tleaf_here + h);
      const double Rd_up = R_d_;
      update_temperature_dependent_params(Tleaf_here - h);
      const double Rd_dn = R_d_;
      vcmax_ = vc0; jmax_ = jm0; gamma_ = ga0; ko_ = ko0; kc_ = kc0;
      R_d_ = rd0; km_ = km0; electron_transport_ = J0;
      const double Rd_T = (Rd_up - Rd_dn) / (2.0 * h);
      // gc = gc_const * E, so dE/dpsi comes from the same spline derivatives the
      // main branch uses; recomputed here because the main branch's locals are
      // below this early return.
      const double gc_c = atm_kpa_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
      const double dgc_ps = gc_c * transport_slope(psi_stem);
      const double dgc_p = -gc_c * transport_slope(psi);
      const double dE_dpsi = (dgc_ps * dpsistem_dpsi + dgc_p) / gc_c;
      dprofit += -Rd_T * dT_dE * dE_dpsi;
    }
    return std::isfinite(dprofit) ? dprofit : 0.0;
  }
  // Past both exits: whatever is returned below is a real derivative, so a zero
  // from here IS a stationary point.
  if (feasible != nullptr) {
    *feasible = true;
  }

  // A'(ci) and C'(psi_stem) by forward-mode tangent of THE MODEL'S OWN algebra -- the
  // same kernels assim_colimited() and hydraulic_cost_TF() are instantiations of,
  // so these are derivatives of the function actually evaluated rather than of a
  // hand-kept mirror of it.
  tangent ci_ad = ci;            seed_direction(ci_ad, 1.0);
  const double A_prime = derivative_along(assim_colimited_kernel(ci_ad));
  tangent ps_ad = psi_stem;      seed_direction(ps_ad, 1.0);
  const double C_prime = derivative_along(hydraulic_cost_TF_kernel(ps_ad));

  // Stomatal-conductance supply coefficient gc and its partials. gc =
  // gc_const * transpiration(psi_stem, psi); transpiration is conductance_max *
  // (transp_from_psi(psi_stem) - transp_from_psi(psi)), so the partials use the
  // analytic spline derivative.
  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
  const double gc = gc_const * transpiration(psi_stem, psi);
  const double dgc_dpsistem = gc_const * transport_slope(psi_stem);
  const double dgc_dpsi = -gc_const * transport_slope(psi);

  // IFT on g(ci; psi_stem, psi): dci/dp = -(dg/dp)/(dg/dci).
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;      // dg/dci
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;
  const double dci_dpsi_expl = -(-dgc_dpsi * (ca_ - ci) * inv_atm) / g_ci;

  // dpsi_stem/dpsi: psi_stem = P(E_psi_stem) with
  //   E_psi_stem = E_up_(psi)/k_max + S(psi),
  // S = transpiration_from_psi, P = psi_from_transpiration (both C1 interpolants), and
  // E_up_(psi) the soil->collar uptake at collar suction psi. The collar variable
  // IS psi now, so there is no dr/dpsi = -1 factor to carry and the two terms add:
  //   dE_psi_stem/dpsi = E_up_'(psi)/k_max + S'(psi)
  //   dpsi_stem/dpsi   = P'(E_psi_stem) * dE_psi_stem/dpsi.
  // E_up_'(psi) is the analytic conductance (dE_from_soil_dpsi_collar), positive;
  // near a branch kink it returns NaN and we fall back to a central difference on
  // the transport.

  const double dci_dpsi = dci_dpsistem * dpsistem_dpsi + dci_dpsi_expl;
  const double base = A_prime * dci_dpsi - C_prime * dpsistem_dpsi;
  // The same expression read as a pair of coefficients, which is what a consumer
  // wanting this condition's gradient needs:
  //
  //   base = (A' dci/dpsi_stem - C') V  +  A' (dci/dpsi at held stem)
  //        =            dprofit_dpsistem_ * V  +  dprofit_dpsi_held_stem_
  //
  // The second is not a remainder. The collar moves the stomatal conductance
  // directly as well as through the stem potential -- transpiration is the stem
  // integral BETWEEN them -- so profit responds to the collar at a held stem
  // potential, and a coordinate pair that reaches the condition only through the
  // stem potential is short by exactly this.
  dprofit_dpsistem_ = A_prime * dci_dpsistem - C_prime;
  dprofit_dpsi_held_stem_ = A_prime * dci_dpsi_expl;
  ci_at_collar_ = ci;
  assim_slope_ = A_prime;
  dci_dpsistem_ = dci_dpsistem;
  dci_dpsi_held_stem_ = dci_dpsi_expl;
  // Gated at the CALL SITE, not just inside the callee: the block is out of line
  // (deliberately, so adding it cannot change FMA contraction in this inlined
  // body), and an out-of-line call costs even when it returns 0.0 immediately.
  if (!use_energy_balance_) {
    return base;
  }
  return base + dprofit_energy_balance_term(ci, gc, g_ci, inv_atm, gc_const,
                                            dgc_dpsistem, dgc_dpsi,
                                            dpsistem_dpsi, dT_dE, Tleaf_here);
}













inline double Leaf::arrh_curve(double Ea, double ref_value, double leaf_temp) const {


  return ref_value*exp(Ea*((leaf_temp+C_to_K) - (25 + C_to_K))/((25 + C_to_K)*gas_constant*(leaf_temp+C_to_K)));
}

inline double Leaf::peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const {
  double arrh = arrh_curve(Ea, ref_value, leaf_temp);
  double arg2 = 1 + exp((d_S*(25 + C_to_K) - H_d)/(gas_constant*(25 + C_to_K)));
  double arg3 = 1 + exp((d_S*(leaf_temp + C_to_K) - H_d)/(gas_constant*(leaf_temp + C_to_K)));

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
          R_d_25, rd_q10_intercept_, rd_q10_slope_};
}

inline void Leaf::update_temperature_dependent_params(double leaf_temp) {
  vcmax_ =
      peak_arrh_curve(vcmax_ha_, vcmax_25, leaf_temp, vcmax_H_d_, vcmax_d_S_);
  jmax_ = peak_arrh_curve(jmax_ha_, jmax_25, leaf_temp, jmax_H_d_, jmax_d_S_);
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

// Explicit leaf energy balance (#523): Tleaf = Tair + (Rn - lambda*E)*ra/(rho*cp).
// E is the hydraulically-pinned transpiration (kg H2O m^-2 s^-1), so lambda*E is
// the latent heat flux (W m^-2) and (Rn - lambda*E) the sensible heat flux H.
inline double Leaf::leaf_temp_from_E(double E, double* dT_dE) const {
  const double Tleaf = Tair_ + (Rn_ - latent_heat_vap * E) * ra_ / vol_heat_cap_air;
  // Clamp to a physical range so an extreme (non-equilibrium) E cannot drive the
  // Arrhenius block non-finite; see leaf_temp_min/max in the header.
  const double clamped = std::min(std::max(Tleaf, leaf_temp_min), leaf_temp_max);
  if (clamped != Tleaf) {
    // Counted, because the kill below keeps the ANALYTIC rows consistent and can
    // do nothing for a consumer differencing this function: two arms both on the
    // clamp return the same temperature, so a difference across them reports no
    // thermal response at all and reports it as a finite number. The count is what
    // lets that consumer notice.
    clamps.note(CLAMP_LEAF_TEMPERATURE);
  }
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

template <typename T>
inline T Leaf::proportion_of_conductivity_kernel(T psi, T b, T c) const {
  return exp(-pow((psi / b), c));
}

template <typename T>
inline T Leaf::proportion_of_conductivity_slope_kernel(T psi, T b, T c) const {
  const T x = pow((psi / b), c);
  return -exp(-x) * c * x / psi;
}

inline double Leaf::proportion_of_conductivity(double psi) const {
  return proportion_of_conductivity_kernel(psi);
}

// set spline for proportion of conductivity
inline void Leaf::setup_transpiration(double resolution) {
  for (const StemCurveCache& hit : *stem_curve_cache_) {
    if (hit.b == stem_b && hit.c == stem_c && hit.resolution == resolution) {
      transpiration_from_psi = hit.from_psi;
      psi_from_transpiration = hit.to_psi;
      stem_b_spline_ = stem_b;
      return;
    }
  }
  // All three channels of G come back from the one loop that forms it: the value,
  // the conductivity G'(psi) = exp(-(psi/stem_b)^stem_c) which is the slope at each
  // knot, and the conductivity's own slope which is the curvature. The slopes used
  // to be recomputed here through proportion_of_conductivity -- the same expression
  // written twice, and now three times would have been needed.
  std::vector<double> x_psi_, y_cumulative_transpiration_, conductivity,
                      conductivity_slope;
  build_cumulative_vulnerability_integral(stem_b, stem_c, resolution, x_psi_,
                                          y_cumulative_transpiration_,
                                          conductivity, conductivity_slope);

  // The inverse carries the reciprocal, d(psi)/dG = 1 / G' (finite everywhere,
  // since G' = 1 at the wet end and stays positive to the dry end), and its own
  // curvature is d2(psi)/dG2 = -G'' / G'^3, which differentiating the reciprocal
  // gives and which needs nothing beyond the two vectors above.
  std::vector<double> inverse_slope(x_psi_.size()),
                      inverse_curvature(x_psi_.size());
  for (size_t i = 0; i < x_psi_.size(); ++i) {
    const double f = conductivity[i];
    inverse_slope[i] = 1.0 / f;
    inverse_curvature[i] = -conductivity_slope[i] / (f * f * f);
  }
  transpiration_from_psi.init(x_psi_, y_cumulative_transpiration_, conductivity,
                              conductivity_slope);
  psi_from_transpiration.init(y_cumulative_transpiration_, x_psi_, inverse_slope,
                              inverse_curvature);

  if (stem_curve_cache_->size() >= curve_cache_size) {
    stem_curve_cache_->erase(stem_curve_cache_->begin());
  }
  stem_curve_cache_->push_back(StemCurveCache{stem_b, stem_c, resolution,
                                              transpiration_from_psi,
                                              psi_from_transpiration});

  // The splines now describe the current stem_b, so the rescaling is over.
  // Recording it HERE rather than at each caller is what makes it impossible to
  // rebuild and forget.
  stem_b_spline_ = stem_b;
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
    const odelia::interpolator::hermite_interpolator<double, 5>& spline, double u,
    double v,
    double scale, const char* spline_name, const char* arg_name,
    const char* caller);

inline double Leaf::eval_stem_curve(
    const odelia::interpolator::hermite_interpolator<double, 5>& spline, double u,
    double scale,
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
    const odelia::interpolator::hermite_interpolator<double, 5>& spline, double u,
    double v, double scale, const char* spline_name, const char* arg_name,
    const char* caller) {
    const bool below = v < spline.min();
    const double lo = scale * spline.min();
    const double hi = scale * spline.max();
    util::stop(std::string("Leaf hydraulics: ") + spline_name +
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
    return transpiration_from_psi.slope(psi);
  }
  return transpiration_from_psi.slope(psi / (stem_b / stem_b_spline_));
}

inline double Leaf::stem_curve_integral_dstem_b(double psi,
                                                const char* caller) const {
  return (stem_curve_integral(psi, caller) -
          psi * stem_curve_integral_deriv(psi)) / stem_b;
}

inline double Leaf::stem_curve_integral_inverse(double w, const char* caller) const {
  const double s = (stem_b == stem_b_spline_) ? 1.0 : stem_b / stem_b_spline_;
  return eval_stem_curve(psi_from_transpiration, w, s,
                         "psi_from_transpiration (the INVERSE cumulative xylem "
                         "conductivity integral G^-1, argument in E/K_max)",
                         "E/K_max", caller);
}

// ⚠️ POLISHING THIS INVERSE IS NOT A LOCAL CHANGE, and it was tried. One Newton
// step against the forward table makes the round trip exact -- 6.27e-16 over 3999
// potentials against 2.10e-12 raw -- and that matters because at zero flux
// transpiration_to_psi_stem asks for G^-1(G(psi_upstream)), whose answer IS
// psi_upstream: round-tripped through two tables it comes back within their
// disagreement instead, and the SIGN of that disagreement decides whether the wet
// end of the collar bracket reads as feasible. Made exact, the suite's failure
// count over the knot counts 1600, 800, 400, 200, 100 goes from 0, 1, 9, 8, 8 to
// 0, 1, 2, 2, 2 -- from a sentinel appearing and vanishing to a number drifting.
//
// What it cost is what stops it: every analytic derivative of this relation --
// dpsistem_dpsi_, and the collar response slope the profit's curvature is built
// from -- was derived for the unpolished map. Polishing the value alone leaves the
// value and its derivatives describing different functions, and the symptom is not
// a small error: 47 of 47 nodes on a dry stand refuse with "the profit's curvature
// has no closed form at this point" at operating points the solve called interior.
// A polished inverse needs the whole chain re-derived from the polished map, which
// is the refusal channel's own work rather than this accessor's.

inline void Leaf::perturb_stem_b(double stem_b_new) {
  check_psi_magnitudes(psi_crit, stem_b_new, roots_.root_b, roots_.root_psi_crit);
  stem_b = stem_b_new;
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
inline double Leaf::transpiration(double psi_stem, double psi_upstream) const {

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
  return atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
}


// biochemical photosynthesis model equations
//ensure that units of PPFD_ actually correspond to something real.
// electron trnansport rate based on light availability and vcmax assuming co-limitation hypothesis
template <typename T>
inline T Leaf::electron_transport_kernel(T ppfd, T quantum_yield,
                                         T curvature, T jmax) const {
  return (quantum_yield * ppfd + jmax - sqrt(pow(quantum_yield * ppfd + jmax, 2) -
  4 * curvature * quantum_yield * ppfd * jmax)) / (2 * curvature);
}

template <typename T>
inline T Leaf::electron_transport_kernel(T ppfd) const {
  return electron_transport_kernel(ppfd, T(a), T(curv_fact_elec_trans),
                                   T(jmax_));
}

inline double Leaf::electron_transport() {



  double electron_transport_ = electron_transport_kernel(PPFD_);

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
inline T Leaf::assim_rubisco_limited_kernel(T ci, T vcmax) const {
  return (vcmax * (ci - gamma_ * umol_per_mol_to_Pa_)) / (ci + km_);
}

template <typename T>
inline T Leaf::assim_rubisco_limited_kernel(T ci) const {
  return assim_rubisco_limited_kernel(ci, T(vcmax_));
}

template <typename T>
inline T Leaf::assim_rubisco_limited_slope_kernel(T ci, T vcmax) const {
  const T g = gamma_ * umol_per_mol_to_Pa_;
  return vcmax * (km_ + g) / ((ci + km_) * (ci + km_));
}

template <typename T>
inline T Leaf::assim_electron_limited_kernel(T ci, T transport) const {
  return transport / 4 *
  ((ci - gamma_ * umol_per_mol_to_Pa_) / (ci + 2 * gamma_ * umol_per_mol_to_Pa_));
}

template <typename T>
inline T Leaf::assim_electron_limited_kernel(T ci) const {
  return assim_electron_limited_kernel(ci, T(electron_transport_));
}

template <typename T>
inline T Leaf::assim_electron_limited_slope_kernel(T ci, T transport) const {
  const T g = gamma_ * umol_per_mol_to_Pa_;
  return transport / 4 * (3 * g) / ((ci + 2 * g) * (ci + 2 * g));
}

// The colimitation itself, lifted out so it can be differentiated in the
// ELECTRON-LIMITED rate as well as in ci -- which is what the light row needs,
// because PPFD reaches assimilation only through that term. The expression is
// unchanged character for character, so the association is too and the golden
// file is the check on that.
template <typename T>
inline T Leaf::colimit_kernel(T assim_rubisco_limited_,
                              T assim_electron_limited_,
                              T curvature, T respiration) const {
  return (assim_rubisco_limited_ + assim_electron_limited_ - sqrt(pow(assim_rubisco_limited_ + assim_electron_limited_, 2) - 4 * curvature * assim_rubisco_limited_ * assim_electron_limited_)) /
             (2 * curvature)- respiration;
}

template <typename T>
inline T Leaf::colimit_kernel(T assim_rubisco_limited_,
                              T assim_electron_limited_) const {
  return colimit_kernel(assim_rubisco_limited_, assim_electron_limited_,
                        T(curv_fact_colim), T(R_d_));
}

template <typename T>
inline T Leaf::colimit_slope_kernel(T assim_rubisco_limited_,
                                    T assim_electron_limited_, T d_rubisco,
                                    T d_electron, T curvature) const {
  const T sum = assim_rubisco_limited_ + assim_electron_limited_;
  const T d_sum = d_rubisco + d_electron;
  const T root = sqrt(pow(sum, 2) - 4 * curvature * assim_rubisco_limited_ *
                                        assim_electron_limited_);
  const T d_root =
      (sum * d_sum - 2 * curvature * (d_rubisco * assim_electron_limited_ +
                                      assim_rubisco_limited_ * d_electron)) /
      root;
  return (d_sum - d_root) / (2 * curvature);
}

template <typename T>
inline T Leaf::assim_colimited_kernel(T ci, T vcmax, T transport,
                                      T curvature, T respiration) const {
  T assim_rubisco_limited_ = assim_rubisco_limited_kernel(ci, vcmax);
  T assim_electron_limited_ = assim_electron_limited_kernel(ci, transport);

  return colimit_kernel(assim_rubisco_limited_, assim_electron_limited_,
                        curvature, respiration);
}

template <typename T>
inline T Leaf::assim_colimited_kernel(T ci) const {
  return assim_colimited_kernel(ci, T(vcmax_), T(electron_transport_),
                                T(curv_fact_colim), T(R_d_));
}

template <typename T>
inline T Leaf::assim_colimited_slope_kernel(T ci, T vcmax, T transport,
                                            T curvature) const {
  return colimit_slope_kernel(
      assim_rubisco_limited_kernel(ci, vcmax),
      assim_electron_limited_kernel(ci, transport),
      assim_rubisco_limited_slope_kernel(ci, vcmax),
      assim_electron_limited_slope_kernel(ci, transport), curvature);
}

inline double Leaf::assim_rubisco_limited(double ci_) {
  return assim_rubisco_limited_kernel(ci_);
}

//calculate the light-limited assimilation rate, returns umol m^-2 s^-1
inline double Leaf::assim_electron_limited(double ci_) {
  return assim_electron_limited_kernel(ci_);
}

// returns co-limited assimilation umol m^-2 s^-1, NET of dark respiration (the
// trailing `- R_d_`), so gross assimilation is this value + R_d_. The comment
// that used to sit on the return statement said the opposite.
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
  // dominated; once PLAN 11a removed that, this became the model's dominant
  // amplifier -- a last-bit change anywhere upstream shifts which point TOMS748
  // lands on inside its tolerance band, and that surfaces in ci, assim, gc and
  // profit. Measured by perturbing the model identically at each tolerance and
  // reading the golden grid (PLAN 11b):
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
  // table implies. Still a good trade next to the 24.5% PLAN 11a bought, and still
  // 21.7% faster than the 3.51 us this package ran at before 11a. If you re-derive
  // the curve, build every point from one tree.
  //
  // Do not tighten further without a use that needs it, and do not loosen it back
  // without re-reading the guide's rounding-versus-bug magnitudes, which this
  // figure sets.
  //
  // Not the same knob as `ci_abs_tol` (the settable control, default 1e-3), which
  // reaches only the off-path optimise_psi_stem_* solvers. A caller tightening
  // that one gets no extra precision here.
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
    util::stop("psi_stem_to_ci failed: " + std::string(e.what()) +
               "; min=" + util::to_string(gamma_ * umol_per_mol_to_Pa_) +
               "; max=" + util::to_string(ca_) +
               "; psi_stem=" + util::to_string(psi_stem) +
               "; psi_upstream=" + util::to_string(psi_upstream));
  }
}

// given psi_stem, find assimilation, transpiration and stomal conductance to c02
inline void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem,
                                                     double psi_upstream,
                                                     double ci_known) {

  if (psi_upstream >= psi_stem){
    ci_ = gamma_*umol_per_mol_to_Pa_;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    } else{
      if(assim_max_ < 0){
        ci_ = gamma_*umol_per_mol_to_Pa_;
        transpiration_ = 0;
        stom_cond_CO2_ = 0;
        } else{
      // Transpiration is the hydraulic supply, independent of ci; compute it
      // first so the PM path can derive the operating-point leaf temperature.
      // Off the PM path this is a memoised no-op reorder (psi_stem_to_ci ->
      // stom_cond_CO2 requests the same (psi_stem, psi_upstream), returning the
      // bit-identical cached value), so the non-PM result is unchanged.
      transpiration_ = transpiration(psi_stem, psi_upstream);
      if (use_energy_balance_) {
        // Tleaf = f(E) is explicit (no PM inversion, no A->E feedback), so this
        // is a single forward pass: recompute the Farquhar temperature params at
        // this candidate's Tleaf before solving for ci. Defeats the photo_temp
        // cache by design -- Tleaf varies per operating point.
        update_temperature_dependent_params(leaf_temp_from_E(transpiration_));
      }
      ci_ = util::is_finite(ci_known) ? ci_known
                                     : psi_stem_to_ci(psi_stem, psi_upstream);
      stom_cond_CO2_ = atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
      }
    }
  assim_colimited_ = assim_colimited(ci_);
}


// Hydraulic cost equations

// Sperry et al. 2017; Sabot et al. 2020 implementation

inline double Leaf::hydraulic_cost_Sperry(double psi_stem, double psi_upstream) {
  // Cost is definitionally zero when the potentials are equal. Returning it
  // explicitly avoids a tiny non-zero residual from FMA contraction of the
  // k_l_soil_ - k_l_stem_ subtraction (arch-dependent; see arm64 build, #468).
  if (psi_stem == psi_upstream) {
    hydraulic_cost_ = 0.0;
    return hydraulic_cost_;
  }
  double k_l_soil_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_upstream);
  double k_l_stem_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem);
  
  hydraulic_cost_ = k_l_soil_ - k_l_stem_;
  
  return hydraulic_cost_;
}

// --- Marginal cost of water -------------------------------------------------

inline double Leaf::lambda_TF24(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  return cost_scale_TF24 * beta2 * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0) *
         pow(1.0 - f, beta2 - 1.0) / leaf_specific_conductance_max_;
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


inline double Leaf::marginal_price_water() {
  const double S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
  if (!std::isfinite(S) || S <= 0.0) {
    return util::na_value;
  }
  const double f_r = proportion_of_conductivity(opt_root_psi_);
  // lambda_multi - lambda_stem, without forming either.
  return lambda_TF24(opt_psi_stem_) * leaf_specific_conductance_max_ * f_r / S;
}

inline double Leaf::g1_eff() const {
  const double chi = ci_ / ca_;
  if (!std::isfinite(chi) || chi >= 1.0) {
    return util::na_value;
  }
  return chi * std::sqrt(atm_vpd_) / (1.0 - chi);
}

// Pure: no write to hydraulic_cost_, so the tangent pass cannot scribble model state
// while probing. The caching is the double entry point's job.
template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(T psi_stem) const {
  return cost_scale_TF24 * pow((1 - proportion_of_conductivity_kernel(psi_stem)), beta2);
}

template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(T psi_stem, T b, T c, T beta,
                                        T scale) const {
  return scale *
         pow((1 - proportion_of_conductivity_kernel(psi_stem, b, c)), beta);
}

template <typename T>
inline T Leaf::hydraulic_cost_TF_slope_kernel(T psi_stem, T b, T c, T beta,
                                              T scale) const {
  const T f = proportion_of_conductivity_kernel(psi_stem, b, c);
  const T df = proportion_of_conductivity_slope_kernel(psi_stem, b, c);
  return -scale * beta * pow((1 - f), beta - 1) * df;
}


inline double Leaf::hydraulic_cost_TF(double psi_stem) {

  hydraulic_cost_ = hydraulic_cost_TF_kernel(psi_stem);

return hydraulic_cost_;
}

// Profit functions

inline double Leaf::profit_psi_stem_Sperry(double psi_stem, double psi_upstream) {

set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_Sperry(psi_stem, psi_upstream);

  return benefit_ - lambda_ * cost;
}


inline double Leaf::profit_psi_stem_TF(double psi_stem, double psi_upstream,
                                       double ci_known) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream, ci_known);

double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF(psi_stem);

  return benefit_ - cost;
}


//optimisation functions


// need docs on Golden Section Search.
inline void Leaf::optimise_psi_stem_Sperry() {

    // These two optimise psi_stem directly and produce no root-collar operating
    // point, so they have no kind to report -- clear it rather than leave the
    // last collar solve's classification standing over their outputs.
    operating_point_kind_ = OperatingPointKind::Unsolved;

    if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }

  opt_psi_stem_ = supply_psi_soil_scalar();


  if ((PPFD_ < 1.5e-8 )| (supply_psi_soil_scalar() > psi_crit)){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    return;
  }

  // Maximise carbon profit over [psi_soil, psi_crit]. Brent's method (golden-
  // section + parabolic interpolation) converges super-linearly on this smooth
  // objective; we minimise -profit and recover the maximum from neg_profit_opt.
    double neg_profit_opt = 0.0;
    opt_psi_stem_ = util::brent_fmin(
        [&](double psi_stem) { return -profit_psi_stem_Sperry(psi_stem, supply_psi_soil_scalar()); },
        supply_psi_soil_scalar(), psi_crit, GSS_tol_abs, &neg_profit_opt);
    profit_ = -neg_profit_opt;

  }
  

inline void Leaf::optimise_psi_stem_TF() {

  // See optimise_psi_stem_Sperry: no root-collar operating point, so no kind.
  operating_point_kind_ = OperatingPointKind::Unsolved;

  if (!supply_is_single_layer()) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }

  opt_psi_stem_ = supply_psi_soil_scalar();

  if (supply_psi_soil_scalar() > psi_crit){
    profit_ = profit_psi_stem_TF(supply_psi_soil_scalar(), supply_psi_soil_scalar());
    return;
  }

  // Maximise carbon profit over [psi_soil, psi_crit] via Brent's method
  // (minimise -profit), matching find_root_collar_psi's multi-layer solver.
    double neg_profit_opt = 0.0;
    opt_psi_stem_ = util::brent_fmin(
        [&](double psi_stem) { return -profit_psi_stem_TF(psi_stem, supply_psi_soil_scalar()); },
        supply_psi_soil_scalar(), psi_crit, GSS_tol_abs, &neg_profit_opt);
    profit_ = -neg_profit_opt;

    return;
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
      util::stop("solve_medlyn_ci_numerical failed: " + std::string(e.what()) +
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

// G(psi) at one scalar. Every slope comes from the same series evaluation the
// value's own table was built out of, so the three are one read rather than
// three.
template <typename T>
inline T Leaf::stem_integral_at(const T& psi, const ProfitInputs<T>& in) const {
  using odelia::util::to_passive;
  const double at = to_passive(psi);
  const double b = to_passive(in.stem_b);
  const double c = to_passive(in.stem_c);
  const VulnerabilityIntegralDerivatives d =
      cumulative_vulnerability_integral_derivatives_at(at, b, c);
  const T step = psi - T(at);
  // The table's value with the query slope and the two trait rows the closed form
  // gives exactly. `step` and both trait offsets are exactly zero in VALUE, so a
  // term multiplying two of them contributes exactly nothing to a first
  // derivative -- and nothing reads a second, so the curve at the passive point is
  // the whole first-order coefficient.
  return T(stem_curve_integral(at, "Leaf::stem_integral_at")) +
         T(proportion_of_conductivity_kernel<double>(at, b, c)) * step +
         T(d.db) * (in.stem_b - T(b)) + T(d.dc) * (in.stem_c - T(c));
}

// The two residuals that place the operating point at a held collar, and the
// profit that follows from them. Written as free templates taking their
// parameters explicitly, because the theorem's denominator is taken with those
// parameters held still and a closure over active values cannot be asked for
// that.
//
//   T1  kmax * (G(sigma) - G(collar)) - E_up = 0     places the stem potential
//   T2  A(ci) * k1 - gc * (ca - ci) * k2   = 0       places the intercellular CO2
//
// gc is proportional to the flux, so it is derived from the uptake here rather
// than passed: two spellings of one fact is a place they can disagree.
template <typename S>
inline Leaf::CollarCoords<S> Leaf::collar_coords_at(
    double sigma_star, double ci_star, const S& collar, const S& flux,
    const ProfitInputs<S>& in) const {
  using odelia::util::to_passive;
  const double gc_per_flux =
      atm_kpa_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const Leaf& leaf = *this;

  // dT1/dsigma is one term of the residual below: the stem curve at the operating
  // point, scaled by kmax. Taken at the passive point, which is where the theorem
  // divides.
  const double dT1_dsigma =
      to_passive(in.kmax) * proportion_of_conductivity_kernel<double>(
                                sigma_star, to_passive(in.stem_b),
                                to_passive(in.stem_c));
  const S sigma = odelia::implicit_value<S>(
      sigma_star, dT1_dsigma, [&](const S& s) -> S {
        return in.kmax * (leaf.template stem_integral_at<S>(s, in) -
                          leaf.template stem_integral_at<S>(collar, in)) -
               flux;
      });

  // ⚠️ FOLLOW THE BRANCH THE FORWARD MODEL TOOK. Where the leaf is not moving
  // water it does not place ci by the residual at all -- it assigns the CO2
  // compensation point, which is a function of temperature and which no carbon
  // input reaches. Applying the theorem to a condition the model did not solve
  // returns a number rather than an error, and the number is a cancellation of
  // two vanishing terms.
  // dT2/dci: the assimilation's own slope plus the supply term's, at the passive
  // point. The marginal profit forms the same quantity and calls it g_ci -- and
  // takes A' the same way, by a tangent through the one kernel rather than
  // through the whole residual.
  double dT2_dci = 0.0;
  if (!ci_at_compensation_point_) {
    const double J0 = electron_transport_kernel<double>(
        to_passive(in.ppfd), to_passive(in.quantum_yield),
        to_passive(in.curv_elec), to_passive(in.transport_jmax));
    tangent c_at = ci_star;
    odelia::ode::seed_direction(c_at, 1.0);
    const double A_prime = odelia::ode::derivative_along(
        assim_colimited_kernel<tangent>(c_at, tangent(to_passive(in.vcmax)),
                                        tangent(J0),
                                        tangent(to_passive(in.curv_colim)),
                                        tangent(to_passive(in.respiration))));
    dT2_dci = A_prime * umol_to_mol +
              gc_per_flux * to_passive(flux) * inv_atm;
  }
  const S ci =
      ci_at_compensation_point_
          ? S(ci_star)
          : odelia::implicit_value<S>(
      ci_star, dT2_dci, [&](const S& c) -> S {
        const S J = leaf.template electron_transport_kernel<S>(
            in.ppfd, in.quantum_yield, in.curv_elec, in.transport_jmax);
        const S A = leaf.template assim_colimited_kernel<S>(
            c, in.vcmax, J, in.curv_colim, in.respiration);
        return A * umol_to_mol -
               S(gc_per_flux) * flux * (S(leaf.ca_) - c) * S(inv_atm);
      });


  return CollarCoords<S>{sigma, ci};
}

// Profit at that point. The coordinates above are the whole of what it shares with
// the marginal, so they are built once and both read them.
template <typename S>
inline S Leaf::profit_at(double sigma_star, double ci_star, const S& collar,
                         const S& flux, const ProfitInputs<S>& in) const {
  const CollarCoords<S> at =
      collar_coords_at<S>(sigma_star, ci_star, collar, flux, in);
  const S& sigma = at.sigma;
  const S& ci = at.ci;
  const S J = electron_transport_kernel<S>(in.ppfd, in.quantum_yield,
                                           in.curv_elec, in.transport_jmax);
  const S A = assim_colimited_kernel<S>(ci, in.vcmax, J, in.curv_colim,
                                        in.respiration);
  const S cost = hydraulic_cost_TF_kernel<S>(sigma, in.stem_b, in.stem_c,
                                             in.beta2, in.cost_scale);
  return A - cost;
}

// M at scalar S, from the same coordinates profit_at uses plus the supply's collar
// derivative -- which is templated for exactly this, and which the transport response
// V inside the assembly needs. `transpiration` is the flux: the T1 residual says the
// stem carries what the soil delivers, so there is no second quantity to form.
template <typename S>
inline S Leaf::marginal_at(const S& collar, const SupplyDraw<S>& draw,
                           const ProfitInputs<S>& in) const {
  const CollarCoords<S> c =
      collar_coords_at<S>(opt_psi_stem_, ci_, collar, draw.flux, in);
  const CollarPoint<S> at{collar, c.sigma, c.ci, draw.slope, draw.flux};
  return marginal_assembled<S>(at, in).marginal;
}


// Every input paired with its own entry in a gradient over the same struct. One
// order, from one table, so the two cannot be matched up by hand.
template <class S>
inline std::vector<odelia::input_and_derivative<S>> against(
    const LeafInputs<S>& in, const LeafInputs<double>& gradient) {
  std::vector<const S*> active = in.field_ptrs();
  std::vector<const double*> rows = gradient.field_ptrs();
  odelia::util::check_length(active.size(), rows.size());
  std::vector<odelia::input_and_derivative<S>> out;
  out.reserve(active.size());
  for (std::size_t k = 0; k < active.size(); ++k) {
    out.push_back({*active[k], *rows[k]});
  }
  return out;
}

// The three conditions that pin a collar. None is a second derivative, so each
// is an ordinary implicit value: the residual at a tangent for its own slope,
// and at S for the rest.
template <typename S>
inline S Leaf::bound_at(WhichBound which, double bound_x,
                        const LeafInputs<S>& in,
                        const SupplyDraw<S>& draw) const {
  check_draw(bound_x, draw);
  const Leaf& leaf = *this;
  if (which == WhichBound::DryRootPsiCrit) {
    // The bound IS the trait, so it moves with that and with nothing else. No
    // theorem to apply and none applied.
    return in.root_psi_crit;
  }
  // Both arms move water, so both slopes are the supply's own conductance at the
  // bound, which the model states analytically rather than differencing.
  using odelia::util::to_passive;
  std::vector<double> soil_at;
  soil_at.reserve(in.supply.psi_soil.size());
  for (const S& v : in.supply.psi_soil) {
    soil_at.push_back(to_passive(v));
  }
  const double dflux_dx = dE_from_soil_dpsi_collar(bound_x, soil_at);

  if (which == WhichBound::Wet) {
    // The residual IS the draw: implicit_value evaluates at the passive bound, and
    // the draw was taken there, so re-recording the supply would put the same
    // expression on the tape a second time.
    return odelia::implicit_value<S>(
        bound_x, dflux_dx, [&](const S&) -> S { return draw.flux; });
  }
  // The dry end is T1 with the stem held at its critical potential: the collar
  // at which the soil delivers exactly what the column can still carry.
  // The stem is held at its critical potential here, so only the second integral
  // and the flux move with the collar.
  const double dT_dx =
      -to_passive(in.profit.kmax) *
          proportion_of_conductivity_kernel<double>(
              bound_x, to_passive(in.profit.stem_b),
              to_passive(in.profit.stem_c)) -
      dflux_dx;
  return odelia::implicit_value<S>(
      bound_x, dT_dx, [&](const S& x) -> S {
        return in.profit.kmax *
                   (leaf.template stem_integral_at<S>(in.psi_crit, in.profit) -
                    leaf.template stem_integral_at<S>(x, in.profit)) -
               draw.flux;
      });
}

// An interior point is the only kind whose condition is a second derivative, so
// it is the only one whose gradient had to be taken in forward mode and handed
// over; every other kind's condition is first order and composes here.
template <typename S>
inline S Leaf::collar_at(const LeafInputs<S>& in,
                         const SupplyDraw<S>& draw) const {
  check_draw(opt_root_psi_, draw);
  S collar = S(opt_root_psi_);
  switch (operating_point_kind_) {
    case OperatingPointKind::Interior: {
      // Closed by the RESIDUAL, like every other kind here. dM/dp is what the
      // theorem divides by; every parameter row is taped from M itself rather than
      // handed over as numbers from a second-order pass over a different assembly.
      // Reported rather than stopped, because at a fold the point is still the point
      // and an output the envelope theorem spares does not read the collar at all.
      // ⚠️ NOT ON THE DOUBLE PATH. implicit_value_reported short-circuits at double
      // -- there is no tape to record against -- but C++ evaluates its arguments
      // first, and dM/dp costs a closed-form assembly of the whole marginal. The
      // forward model was paying for it per interior operating point and discarding
      // it: measured at 191 s against 32 s on a century stand before this line.
      // Recording-only, so both arms place the same collar.
      if constexpr (std::is_same_v<S, double>) {
        collar = S(opt_root_psi_);
      } else {
        collar = odelia::implicit_value<S>(
            opt_root_psi_,
            marginal_collar_slope(in.profit.template rebind_from<double>()),
            [&](const S& y) -> S {
              return marginal_at<S>(y, draw, in.profit);
            });
      }
      break;
    }
    case OperatingPointKind::PinnedWet:
    // Shade death sits ON the wet bound: both potentials at the collar where
    // uptake vanishes, so the bound is the point and it moves with the soil.
    case OperatingPointKind::ShadeDeath:
      collar = bound_at<S>(WhichBound::Wet, opt_root_psi_, in, draw);
      break;
    case OperatingPointKind::PinnedDryRootCrit:
      collar = bound_at<S>(WhichBound::DryRootCrit, opt_root_psi_, in, draw);
      break;
    case OperatingPointKind::PinnedDryRootPsiCrit:
      collar = bound_at<S>(WhichBound::DryRootPsiCrit, opt_root_psi_, in, draw);
      break;
    case OperatingPointKind::HydraulicShutdown:
      // The stem holds at its critical potential and nothing defines the collar
      // at all, so it does not move.
      break;
    default:
      util::stop(std::string("Leaf::collar_at: no operating point to place at "
                             "one that is ") +
                 operating_point_kind_name(operating_point_kind_));
  }
  return collar;
}

template <typename S>
inline Leaf::LeafOutputs<S> Leaf::outputs_at(const S& collar,
                                             const LeafInputs<S>& in,
                                             const SupplyDraw<S>& draw) const {
  check_draw(odelia::util::to_passive(collar), draw);
  LeafOutputs<S> out;

  // The draws first. Profit reads the same flux wherever its collar is the live
  // one, so taking them here is what stops the supply reaching the tape twice.
  S live_flux = S(0.0);
  if (operating_point_kind_ == OperatingPointKind::HydraulicShutdown) {
    // Every flux is zero and stays zero: nothing the soil holds reaches a leaf
    // that has stopped drawing on it.
    out.uptake.assign(static_cast<std::size_t>(supply_n_layers()), S(0.0));
  } else {
    live_flux = E_from_soil_at<S>(collar, in.supply.at(), out.uptake);
  }

  // ⚠️ PROFIT READS THE HELD COLLAR AT AN INTERIOR POINT, and that is the
  // envelope theorem written as an omission rather than as a term that has to
  // come out to zero. The marginal is zero there by the condition the solve
  // drove to nothing, so carrying the collar's movement through profit would add
  // M * dp/dtheta -- a term whose size is the solve's tolerance rather than the
  // model's. It is also why a collar that cannot be recorded costs the water
  // rows and not the objective.
  //
  // The flux follows the collar, which is the whole reason the draw is taken at
  // the passive point: at an interior point it IS profit's flux.
  const bool interior = operating_point_kind_ == OperatingPointKind::Interior;
  const S collar_held(odelia::util::to_passive(collar));
  const S& profit_collar = interior ? collar_held : collar;
  const S& profit_flux = interior ? draw.flux : live_flux;

  if (zero_flux_operating_point()) {
    // No water moves, so no carbon is fixed: assimilation is minus the dark
    // respiration, and the cost is the one at the potential being held.
    const S& held = operating_point_kind_ == OperatingPointKind::HydraulicShutdown
                        ? in.psi_crit
                        : collar;
    out.profit = -in.profit.respiration -
                 hydraulic_cost_TF_kernel<S>(held, in.profit.stem_b,
                                             in.profit.stem_c, in.profit.beta2,
                                             in.profit.cost_scale);
  } else {
    out.profit = profit_at<S>(opt_psi_stem_, ci_, profit_collar, profit_flux,
                              in.profit);
  }
  return out;
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
// `dgc_dT` is a named zero rather than an omission: stom_cond_CO2 currently
// divides by the prescribed AIR vpd, so gc does not depend on Tleaf. Wiring
// leaf-to-air VPD (PLAN 13.1) makes it non-zero and this becomes a one-line
// change instead of a re-derivation.
inline double Leaf::dprofit_energy_balance_term(
    double ci, double gc, double g_ci, double inv_atm, double gc_const,
    double dgc_dpsistem, double dgc_dpsi, double dpsistem_dpsi, double dT_dE,
    double Tleaf) {
  if (!use_energy_balance_ || dT_dE == 0.0) {
    return 0.0;
  }

  // dA/dTleaf by a CENTRAL DIFFERENCE over the temperature block only.
  //
  // Chosen over templating the Arrhenius block on the scalar type, and the
  // reason is the golden file rather than accuracy. Re-expressing
  // update_temperature_dependent_params as the T=double instantiation of a
  // templated kernel changes which expressions share an inlined body, and
  // PLAN 11b measured that exact move changing results even where the algebra
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

  const double dgc_dT = 0.0;  // see the note above; non-zero once PLAN 13.1 lands
  const double damping = (gc * inv_atm - dgc_dT * (ca_ - ci) * inv_atm) / g_ci;

  return A_T * tau * damping;
}

// Post-prepare body of evaluate_root_collar_psi (see header). Kept as a separate
// entry point so callers that evaluate several collar potentials within one step
// (the centred finite difference, #530) can run prepare_collar_solve once and
// reuse the soil-side caches across every profit eval. The clamp into
// [bound_a, bound_b] is identical to evaluate_root_collar_psi's, so near a
// boundary a perturbed potential collapses onto the boundary -- which is exactly
// how the FD path degrades gracefully to a one-sided difference.
inline Leaf::BoundRow Leaf::bound_row(WhichBound which) {
  const std::vector<double>& psi_soil = supply_psi_soil();
  const std::size_t n = psi_soil.size();
  BoundRow row;
  row.d_dpsi_soil.assign(n, 0.0);
  row.d_droot_carbon.assign(n, 0.0);

  if (which == WhichBound::DryRootPsiCrit) {
    // The bound is a registered constant, so it moves with nothing except
    // itself. Exact, and the cheapest row here.
    //
    // ⚠️ THIS ARM RETURNS BEFORE THE SHARED CONVERSION at the end of the
    // function, so the value stored here must ALREADY be the quotient
    // -(dR/du)/(dR/dx) that every other field only becomes down there. The
    // residual is R(x) = x - root_psi_crit, so the raw partial is -1 and the
    // slope is 1, and the row is -(-1)/1 = +1. It read -1 -- the raw partial,
    // never divided or negated -- for as long as nothing differenced it, which
    // is the sign of this row backwards for every consumer. A future arm added
    // above the conversion block has to do its own division here too, or return
    // through it.
    row.bound = supply_psi_crit();
    row.d_droot_psi_crit = 1.0;
    row.residual_slope = 1.0;
    row.finite = true;
    return row;
  }

  // find_root_psi probes collars, and every probe writes E_up_ and
  // soil_consumption_ on its way past. Those describe the OPERATING POINT for
  // whoever solved it, so this restores them: a row is a read, and a read that
  // moves the outputs is the cross-plant channel hazard 8 is about.
  const std::vector<double> saved_consumption = soil_consumption_;
  const double saved_E_up = E_up_;
  const double wettest = supply_begin_solve();
  const double x =
      find_root_psi(wettest, psi_soil, which == WhichBound::Wet ? 0 : 1);
  row.bound = x;
  auto restore = [&]() -> void {
    soil_consumption_ = saved_consumption;
    E_up_ = saved_E_up;
  };
  if (!std::isfinite(x)) {
    restore();
    return row;
  }

  // Common to both: the residual is total uptake, so its state partials are
  // uptake's. dE_from_soil_dpsi_soil is diagonal, so entry j IS dE_up/dpsi_j.
  const double dEup_dx = dE_from_soil_dpsi_collar(x, psi_soil);
  // The root curve reaches total uptake through the layer mean-conductivity
  // integral and nothing else, so each of these is the whole of that parameter's
  // residual partial for either bound. The single path has no root curve.
  const bool multi = supply_kind_ == SupplyKind::MultiLayer;
  const double dEup_droot_b =
      multi ? roots_.duptake_droot_curve(x, psi_soil,
                                         MultiLayerRoots::CurveTrait::Position)
            : 0.0;
  const double dEup_droot_c =
      multi ? roots_.duptake_droot_curve(x, psi_soil,
                                         MultiLayerRoots::CurveTrait::Steepness)
            : 0.0;
  std::vector<double> dEup_dpsi;
  dE_from_soil_dpsi_soil(x, psi_soil, dEup_dpsi);
  // Root carbon is the multi-layer architecture's input; the single-potential
  // path is given a series resistance and has no carbon profile to move.
  std::vector<std::vector<double>> dE_drc, dD_drc;
  const bool has_root_carbon = supply_kind_ == SupplyKind::MultiLayer;
  if (has_root_carbon) {
    roots_.duptake_droot_carbon(x, psi_soil, dE_drc, dD_drc);
  }

  double slope = dEup_dx;
  if (which == WhichBound::DryRootCrit) {
    // The stem's half of the residual. G'(x) is positive, so this can only make
    // the denominator larger -- which is why the sign is safe and the magnitude
    // is not.
    slope += transport_slope(x);
    row.d_dkappa =
        -(stem_curve_integral(psi_crit, "Leaf::bound_row") -
          stem_curve_integral(x, "Leaf::bound_row"));
    row.d_dpsi_crit =
        -transport_slope(psi_crit);
    row.d_dstem_b =
        -leaf_specific_conductance_max_ *
        (stem_curve_integral_dstem_b(psi_crit, "Leaf::bound_row") -
         stem_curve_integral_dstem_b(x, "Leaf::bound_row"));
    // The steepness reaches the same difference of cumulative integrals, from
    // the series rather than from Euler's identity. Taken at the traits the
    // spline was built at, so it is the derivative of the value this residual
    // reads and not of a nearby curve.
    const VulnerabilityIntegralDerivatives G_crit =
        cumulative_vulnerability_integral_derivatives_at(psi_crit, stem_b,
                                                         stem_c);
    const VulnerabilityIntegralDerivatives G_x =
        cumulative_vulnerability_integral_derivatives_at(x, stem_b, stem_c);
    row.d_dstem_c =
        -leaf_specific_conductance_max_ * (G_crit.dc - G_x.dc);
  }
  row.residual_slope = slope;
  if (!std::isfinite(slope) || slope == 0.0) {
    restore();
    return row;
  }

  // Every entry is the same quotient: minus the residual's partial over the
  // residual's slope. The fields above still hold the raw partials at this
  // point, so the conversion happens once, here.
  row.d_droot_b = dEup_droot_b;
  row.d_droot_c = dEup_droot_c;
  bool ok = std::isfinite(row.d_dkappa) && std::isfinite(row.d_dpsi_crit) &&
            std::isfinite(row.d_dstem_b) && std::isfinite(row.d_dstem_c) &&
            std::isfinite(row.d_droot_b) && std::isfinite(row.d_droot_c);
  for (std::size_t j = 0; j < n && ok; ++j) {
    double dEup_drc = 0.0;
    if (has_root_carbon) {
      for (std::size_t i = 0; i < n; ++i) {
        dEup_drc += dE_drc[i][j];
      }
    }
    if (!std::isfinite(dEup_dpsi[j]) || !std::isfinite(dEup_drc)) {
      ok = false;
      break;
    }
    row.d_dpsi_soil[j] = -dEup_dpsi[j] / slope;
    row.d_droot_carbon[j] = -dEup_drc / slope;
  }
  row.d_dkappa = -row.d_dkappa / slope;
  row.d_dpsi_crit = -row.d_dpsi_crit / slope;
  row.d_dstem_b = -row.d_dstem_b / slope;
  row.d_dstem_c = -row.d_dstem_c / slope;
  row.d_droot_b = -row.d_droot_b / slope;
  row.d_droot_c = -row.d_droot_c / slope;
  row.finite = ok;
  restore();
  return row;
}

} // namespace phylloptim

#endif
