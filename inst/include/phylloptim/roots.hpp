// -*-c++-*-
#ifndef PHYLLOPTIM_ROOTS_HPP_
#define PHYLLOPTIM_ROOTS_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/vulnerability.hpp>

#include <odelia/interpolator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace phylloptim {

// The per-layer root hydraulic resistances the supply solve actually consumes.
// Two of the five fields are load-bearing; see set_root_network on why the other
// three are here at all.
struct RootNetwork {
  // Minimum (fully-hydrated) horizontal, intra-layer soil->root resistance.
  // Divided by the vulnerability-weighted mean conductivity at the operating
  // potential to get the actual horizontal resistance.
  std::vector<double> r_R_H_min;
  // Cumulative vertical, inter-layer resistance from the surface down to layer i.
  std::vector<double> r_R_V_sum;
  // Diagnostics only -- read by nothing in the model.
  std::vector<double> c_r_V, c_r_H, r_R_V;
};

// The root-architecture model: how carbon invested in roots becomes hydraulic
// resistance.
//
// ⚠️ THIS IS A HELPER, NOT PART OF THE SUPPLY PATH. Nothing in this package
// calls it. `Leaf::set_physiology` takes the resistances themselves (#33), so
// which root-architecture model is in force is the caller's business, exactly as
// which conductance-versus-height model produced `leaf_specific_conductance_max`
// already was. It stays here, public and tested, because the arithmetic is worth
// sharing and because the golden grid calls it -- see the note on the in-place
// overload below.
//
// Each layer's root carbon is split 1/3 vertical : 2/3 horizontal, and
//   network_.r_R_H_min[i] = beta_R_H / c_r_h        (min horizontal resistance, i.e. the
//                                           reciprocal of max conductance)
//   r_R_V[i]     = beta_R_V * dz^2 / c_r_v (vertical; dz^2 because vertical
//                                           conductivity scales with root
//                                           cross-sectional area)
//   network_.r_R_V_sum[i] = cumulative vertical resistance from the surface to layer i.
//
// The returned vectors are sized to the deepest layer with non-zero root carbon,
// so the hot loop only iterates over layers that actually contain roots.
//
// A free function rather than a member, because it is a *model of the plant* and
// not of water transport: keeping it out of MultiLayerRoots is what lets an
// alternative supply path exist without inventing a root system.
//
// NOTE on zero-carbon layers: a layer with no roots currently gets
// r_R_H_min = 0, i.e. *zero* horizontal resistance, which is infinite
// soil-to-root conductance in a layer with no roots -- backwards. It appears
// unreachable from plant today (max_soil_layer truncates at the last non-zero
// layer, and plant's Q() root distribution does not produce an exact interior
// zero), so this preserves the behaviour rather than changing it silently.
//
// Fills `out` in place so its buffers are reused: a per-solve caller runs this
// once per solve, and building five fresh vectors each time measured
// +0.074 us/call (0.061 -> 0.135), about +2% of a whole solve. That is why this
// overload exists and why it is the one plant uses -- it holds a RootNetwork as a
// strategy member and refills it, the same way it already held the root-carbon
// buffer. The value-returning overload below is for tests and one-off callers,
// where that does not matter.
inline void root_network_from_carbon(
    const std::vector<double>& root_carbon_per_layer, double dz,
    double beta_R_H, double beta_R_V, RootNetwork& out) {
  const size_t n_layers = root_carbon_per_layer.size();

  // deepest layer with non-zero root carbon
  int max_soil_layer = 0;
  for (size_t i = 0; i < n_layers; ++i) {
    if (root_carbon_per_layer[i] != 0) {
      max_soil_layer = i + 1;
    }
  }
  out.c_r_V.assign(max_soil_layer, 0.0);
  out.c_r_H.assign(max_soil_layer, 0.0);
  out.r_R_H_min.resize(max_soil_layer);
  out.r_R_V.resize(max_soil_layer);
  out.r_R_V_sum.resize(max_soil_layer);

  const double dz_sq = dz * dz;
  double vertical_resistance_sum = 0.0;
  for (int i = 0; i < max_soil_layer; ++i) {
    if(root_carbon_per_layer[i] < 0){
            util::stop("Root mass lower than 0");
    }
    const double root_mass = root_carbon_per_layer[i];
    if (root_mass == 0.0) {
      out.r_R_H_min[i] = 0.0;
      out.r_R_V[i] = 0.0;
      out.r_R_V_sum[i] = vertical_resistance_sum;
      continue;
    }

    const double c_r_v = root_mass / 3.0;
    const double c_r_h = root_mass * 2.0 / 3.0;
    out.c_r_V[i] = c_r_v;
    out.c_r_H[i] = c_r_h;

    // Set horizantal minimum resistance per soil layer (i.e. reciprocal of maximum conductance).
    out.r_R_H_min[i] = beta_R_H / c_r_h;
    // The vertical conductivity is likely linearly proportional to the root area projected onto the horizontal plane, hence dz^2.
    out.r_R_V[i] = beta_R_V * dz_sq / c_r_v;
    vertical_resistance_sum += out.r_R_V[i];
    out.r_R_V_sum[i] = vertical_resistance_sum;
  }
}

// Same, returning a fresh network. Convenience for tests and standalone callers.
inline RootNetwork root_network_from_carbon(
    const std::vector<double>& root_carbon_per_layer, double dz,
    double beta_R_H, double beta_R_V) {
  RootNetwork out;
  root_network_from_carbon(root_carbon_per_layer, dz, beta_R_H, beta_R_V, out);
  return out;
}

// Layer thickness implied by a cumulative soil-depth profile.
//
// Exported as its own function because since #33 there are TWO callers that must
// agree on it: MultiLayerRoots::set_soil_state, and whoever builds the root
// network before handing it over. root_network_from_carbon scales the vertical
// resistance by dz^2, so two definitions drifting apart would put a silent
// squared factor on every vertical resistance -- and nothing in either package
// would notice, because both halves would still be internally consistent.
inline double layer_thickness(const std::vector<double>& soil_depth) {
  if (soil_depth.empty()) {
    // Guarded because this is a public entry point a caller reaches for directly,
    // where the same expression inside set_soil_state was only ever reachable
    // after set_physiology had validated the profile against psi_soil.
    util::stop("layer_thickness: soil_depth must have at least one layer");
  }
  return soil_depth.back() / soil_depth.size();
}

// ---------------------------------------------------------------------------
// SOIL -> ROOT-COLLAR WATER SUPPLY
// ---------------------------------------------------------------------------
// The water supply side of the model, lifted out of Leaf (issue #2). Everything
// here answers one question: given a collar potential, how much water can the
// root system deliver, and how fast does that change? The gas-exchange core does
// not need to know that soil has layers -- it only ever consumes
//
//   E_up  = uptake(P_collar)        and       dE_up/dP = duptake_dpsi(P_collar)
//
// so this is where a single-potential alternative will plug in (stage 3).
//
// Scientific model (after Potkay et al. 2021; prototyped in plant's
// vignettes/models/root_water_uptake.Rmd as E_from_Soil_to_Root_Collar):
//
// The root system is represented as a set of parallel soil layers, each
// connected to a single root collar (the point where roots join the stem).
// Within each layer i, water flows from soil to collar driven by the water
// suction gradient (T_collar - T_soil[i]), corrected for the gravitational
// head needed to lift water to the layer midpoint (gravity_head * z_soil_mid).
//
// The hydraulic resistance of each layer is the sum of two terms:
//   * r_R_H : horizontal (intra-layer, soil->root) resistance. Set during
//             set_root_network as network_.r_R_H_min[i] / f_r, where r_R_H_min scales
//             with the carbon invested in horizontal roots and f_r is the
//             fractional loss of conductivity from the root vulnerability
//             curve at the operating potential.
//   * r_R_V : vertical (inter-layer, along the root axis to the collar)
//             resistance, accumulated from the surface down to layer i
//             (r_R_V_sum). It scales with dz^2 / carbon-in-vertical-roots.
//
// Because the root vulnerability curve f_r is non-linear in psi, the
// horizontal resistance is evaluated using the *average* fractional
// conductivity over the potential interval spanned between the soil and the
// collar (T_src_min..T_src_max). This mean is obtained as
// (1/(b-a)) * integral_a^b f_r dpsi from a pre-integrated curve
// (root_vuln_integral_from_psi) with two spline evals, the same technique used
// for stem transpiration in Leaf::setup_transpiration.
//
// SIGN CONVENTION: every water potential here is a POSITIVE MAGNITUDE in MPa,
// as everywhere else in this package (#25). psi_soil_ arrives that way and stays
// that way -- there is no second representation and no flip. The vulnerability
// splines are indexed by magnitude, so they are read directly. Where an equation
// needs one potential to oppose another the minus sign is written in the
// equation: the collar draws water when its magnitude exceeds the soil's, so the
// flux numerator is (T_collar - T_soil - gravity_head).
//
// TWO OUTPUTS ARE NOT OWNED HERE. E_up and the per-layer soil_consumption
// buffer are passed in by reference rather than stored, because plant reaches
// into `leaf.E_up_` / `leaf.soil_consumption_[a]` by name and *writes back*
// into them after crown integration (tf24_strategy.cpp:501-508). They are
// plant's buffers, not this object's state. Note the deliberate unit split:
// E_up is kg H2O m^-2 s^-1, soil_consumption[i] is mol, converted downstream.
class MultiLayerRoots {
public:
  // --- root vulnerability trait pair (hazard 1: NOT the stem's b/c) ---------
  double root_c = 2.680147;       // unitless
  double root_b = 3.898245;       // -MPa
  double root_psi_crit = 5.870283; // -MPa

  // NOTE: beta_R_H and beta_R_V used to live here. They are parameters of the
  // root-architecture model, not of water transport, and since #33 this class
  // takes resistances rather than the carbon they were applied to -- so they are
  // arguments to root_network_from_carbon and members of whoever owns that
  // model. In plant that is TF24_Strategy.

  // pre-computed root vulnerability curve f_r(m) = exp(-(m/root_b)^root_c)
  odelia::interpolator::Interpolator root_vuln_from_psi;
  // cumulative integral of it, G(m) = int_0^m f_r(s) ds, indexed by magnitude
  // m = -psi. Lets uptake() obtain the mean conductivity over a potential
  // interval from 2 evals instead of (n+1).
  odelia::interpolator::Interpolator root_vuln_integral_from_psi;
  // What the two splines do past their last knot, cached by setup_vulnerability.
  // Both are read through the accessors below, never as bare .eval() calls --
  // see the accessors for why. util::na_value until a curve is built.
  double root_vuln_last_knot_ = util::na_value;
  double root_vuln_integral_limit_ = util::na_value;

  // --- soil geometry -------------------------------------------------------
  // The four scalars carry the same unset sentinels clear() assigns, so a bare
  // MultiLayerRoots is never indeterminate. Leaf reaches them only after
  // set_soil_state / set_root_network, but the class is public now.
  double soil_number_of_depths_ = util::na_value_int;
  int max_soil_layer = util::na_value_int;  // deepest layer with non-zero root mass
  std::vector<double> soil_depth_;
  std::vector<double> z_soil_mid_;
  // Per-layer gravitational head gravity_head * z_soil_mid_[i], precomputed once
  // per set_soil_state (z_soil_mid_ is fixed across a collar solve). Used three
  // times per layer in uptake()'s hot loop; caching it removes a redundant
  // multiply per layer per (re)evaluation.
  std::vector<double> grav_head_z_;
  bool use_precomputed_z_soil_mid_ = false;
  double dz_ = util::na_value;

  // --- soil state ----------------------------------------------------------
  std::vector<double> psi_soil_;           // positive magnitudes, as supplied
  // Per-layer cache of root_vuln_integral_from_psi.eval(psi_soil_[i]).
  // psi_soil_ is fixed for the whole collar solve, so the soil-side
  // endpoint of the cumulative-integral lookup in uptake() is constant across
  // every (re)evaluation of the nested root-finders. Precomputing it once per
  // solve (alongside the T_collar-side eval, hoisted out of the layer loop)
  // collapses ~2 spline evals per layer to ~1 per call. Rebuilt in begin_solve.
  std::vector<double> root_vuln_integral_soil_;

  // --- root resistance network --------------------------------------------
  // Held as one object rather than five loose vectors so the carbon -> resistance
  // map can fill it in place and reuse its buffers. The solve reads exactly two
  // of its fields.
  RootNetwork network_;

  // -------------------------------------------------------------------------

  // Reset every state member to the unset sentinel. Mirrors Leaf::setup_clean_leaf.
  void clear() {
    psi_soil_.clear();
    soil_depth_.clear();
    z_soil_mid_.clear();
    grav_head_z_.clear();
    use_precomputed_z_soil_mid_ = false;
    network_.c_r_V.clear();
    network_.c_r_H.clear();
    network_.r_R_H_min.clear();
    network_.r_R_V.clear();
    network_.r_R_V_sum.clear();
    soil_number_of_depths_ = util::na_value_int;
    max_soil_layer = util::na_value_int;
  }

  // Pre-compute the root vulnerability curve and its cumulative integral over
  // [0, psi_max_root], the range where conductivity drops to 1%. Avoids repeated
  // exp(pow(...)) inside uptake().
  //
  // A LAYER DRIER THAN THE GRID IS AN ORDINARY STATE, NOT AN ERROR: the grid
  // stops where conductivity has fallen to 1% (~6.8 MPa at the root defaults),
  // while the soil potential a caller may pass is bounded only by its own ceiling
  // -- 1000 MPa in plant. So each curve needs a stated answer out there, and the
  // two answers differ because the two limits do. Both are applied by the
  // accessors below, which are the only way the rest of this class reads a curve.
  // Issue #1 measured what taking odelia's extrapolant unexamined cost instead: a
  // conductivity that turns negative past 7.3742 MPa, and an integral 4.35x its
  // own limit at 1000 MPa.
  void setup_vulnerability(double resolution) {
    std::vector<double> x_psi_root, y_integral;
    cumulative_vulnerability_integral(root_b, root_c, resolution, x_psi_root,
                                      y_integral);

    // f_r conductivity knots on the same grid. f_r(0) = exp(-pow(0,root_c)) = 1.
    std::vector<double> y_f_r(x_psi_root.size());
    for (size_t i = 0; i < x_psi_root.size(); ++i) {
      y_f_r[i] = exp(-pow(x_psi_root[i] / root_b, root_c));
    }
    // Conductivity: NO extrapolation, matching the stem pair in
    // Leaf::setup_transpiration. Its limit is zero, which is unusable as the
    // divisor it becomes, so root_vuln_at clamps the argument to the last knot
    // instead and this setting is what stops any other reading being possible.
    root_vuln_from_psi.init(x_psi_root, y_f_r);
    root_vuln_from_psi.set_extrapolate(false);

    // Integral: extrapolation stays ON, under a ceiling. Its limit is finite and
    // non-zero, and the end-knot polynomial tracks the true G closely over the
    // half-MPa it takes to reach that limit (3.46212 against 3.46176 at 7 MPa), so
    // capping the VALUE is both smooth and tighter than clamping the argument
    // would be. root_vuln_integral_at applies the cap.
    //
    // Do not set this false "to match" the conductivity spline: odelia's deriv()
    // has no extrapolation check where eval() does, so eval would throw while
    // deriv went on extrapolating.
    root_vuln_integral_from_psi.init(x_psi_root, y_integral);
    root_vuln_integral_from_psi.set_extrapolate(true);

    // The last knot, NOT vulnerability_psi_max: the knot loop advances by
    // accumulation and stops one step short of psi_max (6.8229 against 6.8918 at
    // the root defaults), so psi_max is itself outside the domain and clamping to
    // it would throw.
    root_vuln_last_knot_ = root_vuln_from_psi.max();
    root_vuln_integral_limit_ =
        cumulative_vulnerability_integral_limit(root_b, root_c);
  }

  // f_r at a suction, clamped into the knot domain -- set_extrapolate(false)
  // throws at both ends, and the last knot's ~1% is the driest conductivity this
  // curve describes.
  //
  // Operand order is load-bearing: written this way both clamps return psi when
  // psi is NaN, where the reversed forms return the bound. The uptake call site's
  // !isfinite(f_ri) guard is what reads that NaN.
  double root_vuln_at(double psi) const {
    return root_vuln_from_psi.eval(
        std::max(std::min(psi, root_vuln_last_knot_), 0.0));
  }

  // G at a suction, with the closed-form limit G(inf) = (b/c)*Gamma(1/c) as a
  // ceiling. The spline sits below that ceiling at every knot (99.83% of it at
  // the last), so this is continuous, monotone, and identical to a bare eval
  // everywhere on the grid: it binds only past 7.3132 MPa, where the end-knot
  // polynomial would otherwise carry on accumulating for ever.
  //
  // WHY IT MATTERS: the general branch of uptake_impl forms the layer's mean
  // resistance as r_R_H = r_R_H_min * span / integral. Unbounded, the integral
  // makes r_R_H FALL as the layer dries, so the drier a near-embolised layer, the
  // harder the plant supposedly pumps water into it -- and because whole-plant
  // shutdown keys off the wettest layer, nothing stops it. Capped, the reverse
  // flux tends to integral/r_R_H_min <= G(inf)/r_R_H_min as the span grows, which
  // is the whole area under the conductivity curve and is the right limit.
  double root_vuln_integral_at(double psi) const {
    return std::min(root_vuln_integral_from_psi.eval(psi),
                    root_vuln_integral_limit_);
  }

  // dG/dpsi, consistent with root_vuln_integral_at: zero wherever that returns
  // the cap, because there the value no longer depends on psi. Costs a second
  // spline eval, which is affordable here and would not be on the uptake hot path
  // -- duptake_dpsi runs once per gradient, not ~10^3 times per solve.
  double root_vuln_integral_deriv_at(double psi) const {
    if (root_vuln_integral_from_psi.eval(psi) >= root_vuln_integral_limit_) {
      return 0.0;
    }
    return root_vuln_integral_from_psi.deriv(psi);
  }

  // Per-timestep soil state: the layer potentials, the layer depths, and the
  // gravitational head that follows from them.
  void set_soil_state(const std::vector<double>& psi_soil,
                      const std::vector<double>& soil_depth) {
    psi_soil_ = psi_soil;
    soil_depth_ = soil_depth;
    soil_number_of_depths_ = soil_depth_.size();

    if (!(use_precomputed_z_soil_mid_ &&
          z_soil_mid_.size() == static_cast<size_t>(soil_number_of_depths_))) {
      // Fallback for paths that do not provide environment-precomputed midpoints.
      z_soil_mid_.resize(soil_number_of_depths_);
      for (size_t i = 0; i < soil_number_of_depths_; ++i) {
        if (i == 0) {
          z_soil_mid_[i] = (soil_depth_[i] / 2.0);
        } else {
          z_soil_mid_[i] = ((soil_depth_[i - 1] + soil_depth_[i]) / 2.0);
        }
      }
    }

    use_precomputed_z_soil_mid_ = false;

    grav_head_z_.resize(soil_number_of_depths_);
    for (size_t i = 0; i < soil_number_of_depths_; ++i) {
      grav_head_z_[i] = gravity_head * z_soil_mid_[i];
    }

    // Layer thickness is soil geometry, not root architecture, so it is set here
    // rather than alongside the resistance network that consumes it.
    //
    // ⚠️ Since #33 nothing in this package READS dz_: the only thing that did was
    // the carbon -> resistance map, which is now the caller's. It is kept because
    // it is a property of the soil profile this object is given, and because the
    // caller needs the same number -- see layer_thickness, which is the shared
    // definition. It is a removal candidate with the diagnostics (item 6).
    dz_ = layer_thickness(soil_depth_);
  }

  // Per-timestep root resistance network. Takes the resistances themselves, not
  // the root carbon they are derived from.
  //
  // WHY THIS TAKES RESISTANCES. The solve reads exactly two of these vectors --
  // r_R_H_min and r_R_V_sum -- plus grav_head_z_ and max_soil_layer. Nothing in
  // uptake() or duptake_dpsi() touches root carbon, the 1/3 : 2/3 split, dz, or
  // either beta_R_* constant; those are inputs to a *root architecture* model
  // that happens to run just before. Splitting them out is the same move
  // leaf_specific_conductance_max already makes: plant computes
  // kmax = K_s*theta/(h*eta_c) and hands over a scalar, so which
  // conductance-versus-height model is in force is not this package's business.
  // The carbon -> resistance map lives in root_network_from_carbon below, so an
  // alternative supply path can supply resistances any way it likes.
  //
  // c_r_V_, c_r_H_ and r_R_V ride along as diagnostics: nothing in the model
  // reads them, but plant exposes them through RcppR6, so they are carried
  // rather than dropped. They are removal candidates with item 6.
  //
  // ⚠️ TAKES const& AND COPY-ASSIGNS, WHERE IT USED TO TAKE BY VALUE AND MOVE.
  // The move was right when the caller built a throwaway network per call; it is
  // WRONG now that the caller holds one as a member and refills it, because
  // moving would empty the caller's buffers and force root_network_from_carbon to
  // reallocate all five vectors on the next call -- reintroducing exactly the
  // +0.074 us the in-place overload exists to avoid. Copy-assigning into
  // already-sized vectors allocates nothing on either side once both are warm.
  void set_root_network(const RootNetwork& network) {
    if (network.r_R_V_sum.size() != network.r_R_H_min.size()) {
      util::stop("set_root_network: r_R_H_min and r_R_V_sum must have the same "
                 "length; got " + std::to_string(network.r_R_H_min.size()) +
                 " and " + std::to_string(network.r_R_V_sum.size()));
    }
    // The rooted layers index psi_soil_ and grav_head_z_ directly in uptake(),
    // so a network deeper than the soil profile is an out-of-bounds read rather
    // than a wrong number. Before #33 the length agreement came for free, because
    // set_physiology validated root carbon against soil_depth; now the network
    // arrives from outside the package and the check has to be here.
    if (soil_number_of_depths_ > 0 &&
        network.r_R_H_min.size() >
            static_cast<size_t>(soil_number_of_depths_)) {
      util::stop("set_root_network: network has " +
                 std::to_string(network.r_R_H_min.size()) +
                 " rooted layers but the soil profile has only " +
                 std::to_string(static_cast<int>(soil_number_of_depths_)));
    }
    for (size_t i = 0; i < network.r_R_H_min.size(); ++i) {
      // Resistances, so non-negative. Zero is permitted, and means infinite
      // conductance: root_network_from_carbon produces it for a zero-carbon
      // layer, which is backwards but is the behaviour this preserves (see its
      // note on zero-carbon layers).
      if (!std::isfinite(network.r_R_H_min[i]) || network.r_R_H_min[i] < 0.0 ||
          !std::isfinite(network.r_R_V_sum[i]) || network.r_R_V_sum[i] < 0.0) {
        util::stop("set_root_network: root resistances must be finite and "
                   "non-negative; layer=" + std::to_string(i) +
                   "; r_R_H_min=" + util::to_string(network.r_R_H_min[i]) +
                   "; r_R_V_sum=" + util::to_string(network.r_R_V_sum[i]));
      }
    }
    network_ = network;
    max_soil_layer = static_cast<int>(network_.r_R_H_min.size());
  }

  // Per-solve entry point. Builds the soil-side cumulative-integral cache and
  // returns the wettest rooted layer -- which in magnitudes is the layer of
  // SMALLEST suction, a minimum where the signed convention took a maximum. That
  // is the bracket endpoint the collar solve needs. The cache is valid until the
  // next call.
  //
  // This is a single pass on purpose: the cache is the measured hot-path
  // optimisation described on root_vuln_integral_soil_, and the wettest layer
  // falls out of the same loop.
  double begin_solve() {
    root_vuln_integral_soil_.resize(max_soil_layer);
    double wettest_soil_layer = std::numeric_limits<double>::infinity();
    for (int i = 0; i < max_soil_layer; ++i) {
      root_vuln_integral_soil_[i] = root_vuln_integral_at(psi_soil_[i]);
      wettest_soil_layer = std::min(wettest_soil_layer, psi_soil_[i]);
    }
    return wettest_soil_layer;
  }

  // Uptake at a collar suction, against the soil state begin_solve() cached.
  // This is the hot path: ~10^3 calls per collar solve.
  void uptake(double T_collar, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake_impl(T_collar, psi_soil_,
                root_vuln_integral_soil_.size() ==
                    static_cast<size_t>(max_soil_layer),
                soil_consumption, E_up);
  }

  // Uptake against an arbitrary vector of layer suctions, for callers that
  // want to probe the supply function away from the current soil state (the
  // R-facing Leaf::E_from_Soil_to_Root_Collar).
  //
  // The `&psi_soil == &psi_soil_` test is what used to select the
  // cached path for every caller, including the hot one. It is kept here only so
  // this entry point cannot change behaviour for a caller that happens to hand
  // back psi_soil_ itself; the hot path above no longer depends on
  // address identity to be fast. PLAN 7b-ii trap 3.
  void uptake_at(double T_collar, const std::vector<double>& psi_soil,
                 std::vector<double>& soil_consumption, double& E_up) const {
    uptake_impl(T_collar, psi_soil,
                (&psi_soil == &psi_soil_) &&
                    root_vuln_integral_soil_.size() ==
                        static_cast<size_t>(max_soil_layer),
                soil_consumption, E_up);
  }

  // Analytic d(E_up)/d(T_collar): the collar-suction derivative of the uptake,
  // mirroring the general branch of uptake_impl. It is a CONDUCTANCE and is
  // positive by construction -- pulling harder at the collar draws more water --
  // which is the whole reason for working in magnitudes (#25). Per layer, with
  // span = |T_collar - T_soil[i]| and integral = \int f_r over
  // [T_src_min, T_src_max] (root_vuln_integral_from_psi, whose integrand is
  // root_vuln_from_psi):
  //   E_i        = (T_collar - T_soil[i] - grav) / r_R,
  //   r_R        = r_R_H_min[i] * span / integral + r_R_V_sum[i],
  //   dspan/dT   = sign_var   (+1 if T_collar is the upper bound, else -1),
  //   dinteg/dT  = sign_var * f_r(T_collar)  for T_collar>0  (else sign_var, f_r==1),
  // and dE_i/dT follows by the quotient rule.
  //
  // CONTRACT: returns NaN when any layer sits on a branch kink (T_collar ==
  // T_soil[i], the gravity-balance point, or T_collar == 0). That is deliberate,
  // not a failure -- the analytic general-branch derivative is not valid across
  // those, and the caller falls back to a central difference. An implementation
  // that threw, or returned 0, would silently degrade TF24f's acclimation
  // gradient. Any alternative supply path must keep this contract.
  double duptake_dpsi(double T_collar,
                      const std::vector<double>& psi_soil) const {
    const double kink_tol = 1e-8;
    double dEup_dT_mol = 0.0;

    for (int i = 0; i < max_soil_layer; i++) {
      if (std::abs(T_collar - psi_soil[i]) < kink_tol ||
          std::abs((T_collar - psi_soil[i]) - grav_head_z_[i]) < kink_tol ||
          std::abs(T_collar) < kink_tol) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;
      const double sign_var = (T_collar > psi_soil[i]) ? 1.0 : -1.0;  // = dspan/dT_collar

      // integral, replicated bit-for-bit from uptake_impl.
      const double T_pos_lo = std::max(T_src_min, 0.0);
      const double T_neg_hi = std::min(T_src_max, 0.0);
      double integral = 0.0;
      if (T_pos_lo < T_src_max) {
        integral += root_vuln_integral_at(T_src_max) -
                    root_vuln_integral_at(T_pos_lo);
      }
      if (T_src_min < T_neg_hi) {
        integral += (T_neg_hi - T_src_min);
      }

      // d(integral)/d(T_collar): for T_collar>0 the moving bound is in the
      // vulnerable region. The integrand is the derivative of the *same*
      // cumulative curve that produced `integral`
      // (root_vuln_integral_deriv_at), NOT the separate root_vuln_from_psi
      // spline: the two agree on the knot domain but are bounded differently past
      // it -- the conductivity lookup clamps its argument to the last knot, the
      // integral is capped at G(inf) -- so beyond the domain only the integral's
      // own derivative stays consistent with the value used here (issue #1; the
      // reasoning is #527's, the "both clamp-to-last-value" it used to cite was
      // never true of either). For T_collar<0 (an above-atmospheric collar) the
      // moving bound is in the f_r==1 part, contributed linearly, so the slope
      // is 1.
      const double fr_at =
          (T_collar > 0.0) ? root_vuln_integral_deriv_at(T_collar) : 1.0;
      const double dinteg_dT = sign_var * fr_at;

      const double r_R_H = network_.r_R_H_min[i] * span / integral;
      const double r_R = r_R_H + network_.r_R_V_sum[i];
      const double dr_R_H_dT =
          network_.r_R_H_min[i] * (sign_var * integral - span * dinteg_dT) / (integral * integral);
      const double dr_R_dT = dr_R_H_dT;

      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      const double dnum_dT = 1.0;
      // E_i = num / r_R  ->  quotient rule.
      dEup_dT_mol += (dnum_dT * r_R - num * dr_R_dT) / (r_R * r_R);
    }

    return dEup_dT_mol * kg_per_mol_h2o;  // match E_up's kg units
  }

private:
  // Total water drawn from all layers to the collar. Writes E_up (kg H2O m^-2
  // leaf s^-1) and soil_consumption[i] (mol H2O m^-2 leaf s^-1, note the unit
  // split); a negative E_i in a layer means that layer is *gaining* water
  // (hydraulic redistribution).
  //
  // Implementation decisions:
  //   * f_r and its running integral are read from pre-computed splines
  //     (root_vuln_from_psi, root_vuln_integral_from_psi) instead of repeatedly
  //     evaluating exp(-(psi/b)^c); see setup_vulnerability.
  //   * Two special cases are handled exactly to avoid division/round-off
  //     issues: (a) collar potential equals layer potential, and (b) the
  //     gradient exactly balances gravity (E_i = 0).
  //   * The isfinite() guards are present because this is called from within
  //     nested root-finders where bad brackets can produce NaNs; they fail fast
  //     with diagnostic context rather than propagating NaN.
  void uptake_impl(double T_collar, const std::vector<double>& psi_soil,
                   bool use_integral_cache,
                   std::vector<double>& soil_consumption, double& E_up) const {

    if (!std::isfinite(T_collar)) {
      util::stop("E_from_Soil_to_Root_Collar invalid input; T_collar=" + util::to_string(T_collar));
    }

    E_up = 0;

    // Cumulative-integral spline caching (bit-identical fast path). The only two
    // arguments ever passed to root_vuln_integral_from_psi in the loop below are
    // T_src_max and T_pos_lo, each of which resolves to exactly one of
    // {psi_soil[i], T_collar, 0}. T_collar is constant across all layers (compute
    // once), and psi_soil[i] is constant across the whole solve (precomputed in
    // begin_solve).
    const double G_at_T_collar =
        use_integral_cache ? root_vuln_integral_at(T_collar) : 0.0;

    // GUARD POLICY (the per-layer isfinite/stop guards here were added while
    // debugging the #485 drought-NaN, now fixed at source by the soil residual-
    // moisture floor). Most were defensive and redundant, so they have been
    // removed from this hot loop; the remaining two are load-bearing:
    //   * the equal-potentials f_ri <= 0 check below: it USED to be the only
    //     thing standing between a deep-drought layer and the negative
    //     conductivity root_vuln_from_psi's extrapolant produced past its
    //     domain -- negative-but-FINITE r_R, so a wrong-sign E_i the post-loop
    //     isfinite(E_up) net would not catch. That case is now prevented at
    //     source: root_vuln_at clamps its argument to the last knot, so f_ri is
    //     bounded below by the last knot's ~1% (issue #1). The check stays as a
    //     cheap assertion on that, not as the enforcement.
    //   * the post-loop isfinite(E_up) check: any non-finite produced anywhere
    //     in the loop propagates into the sum and is caught there once per call.
    // Everything else is provably safe to drop on the valid path: psi_soil is
    // validated in Leaf::set_physiology; T_src_min<=T_src_max by construction;
    // the general-branch integral comes from a monotone-increasing spline so it
    // is strictly > 0 (span>0), giving r_R>0 and finite E_i; and any stray
    // NaN/Inf still reaches the post-loop net.
    for(int i = 0; i < max_soil_layer; i++){

    // The wetter end of the interval spanned between this layer and the collar --
    // the SMALLER suction, where the signed convention took a minimum.
    double T_src_min = std::min(psi_soil[i], T_collar);

    // The drier end: the LARGER suction.
    double T_src_max = std::max(psi_soil[i], T_collar);

     // If root collar soil water potential equals the soil water potential in a given layer
    if(std::abs(T_collar - psi_soil[i]) < 1e-8){

      // Fraction of conductance in roots in a given layer at the driest suction
      // (which here equals the root collar's).
      // root_vuln_from_psi is a pre-built spline of exp(-(psi/b_root)^c_root),
      // read through root_vuln_at so a layer drier than the grid gets the last
      // knot's conductivity rather than an extrapolated (eventually negative) one.
      double f_ri = root_vuln_at(T_src_max);
      if (!std::isfinite(f_ri) || f_ri <= 0.0) {
        util::stop("E_from_Soil_to_Root_Collar invalid f_ri; layer=" + std::to_string(i) +
                   "; f_ri=" + util::to_string(f_ri) +
                   "; T_src_max=" + util::to_string(T_src_max) +
                   "; T_collar=" + util::to_string(T_collar));
      }

      // Fraction of conductance in roots in a given layer at the driest suction
      double r_R_H = network_.r_R_H_min[i] / f_ri; // [MPa * s * (mol H2O)^-1]

      // Total root resistance (horizantal plus vertical)
      double r_R = r_R_H + network_.r_R_V_sum[i];

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      double E_i = -grav_head_z_[i] / r_R ;

      soil_consumption[i] = E_i;
      E_up += E_i;

    }
    else if(std::abs((T_collar - psi_soil[i]) - grav_head_z_[i]) < 1e-8){
      // If pressure difference perfectly balances gravity transpiration is equal to zero
      double E_i = 0.0; // [mol H2O / m^2 / s]

      soil_consumption[i] = E_i;

      E_up += E_i;

    } else{

      // Mean fractional root conductivity over the suction interval
      // [T_src_min, T_src_max], i.e. (1/(b-a)) * integral_a^b f_r dT.
      // Computed from the pre-integrated curve G(m) = integral_0^m f_r(s) ds
      // (root_vuln_integral_from_psi, indexed by the suction magnitude) with 2
      // evals instead of the old (n+1)-point sample mean. The interval is split
      // at T = 0: for T < 0 (an above-atmospheric potential) vulnerability is 1.
      double T_pos_lo = std::max(T_src_min, 0.0); // wet end of the T>=0 part
      double T_neg_hi = std::min(T_src_max, 0.0); // dry end of the T<0 part

      // Memoised cumulative-integral lookup. Returns the exact same double the
      // spline would (same input -> same output); the comparisons select the
      // precomputed value because T_src_max / T_pos_lo are bit-for-bit equal to
      // one of the cached arguments in the common (T>=0) case.
      auto G_integral = [&](double arg) -> double {
        if (use_integral_cache) {
          if (arg == T_collar) return G_at_T_collar;
          if (arg == psi_soil[i]) return root_vuln_integral_soil_[i];
        }
        return root_vuln_integral_at(arg);
      };

      double integral = 0.0;
      if (T_pos_lo < T_src_max) {
        // T>=0 part: suction runs from T_pos_lo up to T_src_max
        integral += G_integral(T_src_max) - G_integral(T_pos_lo);
      }
      if (T_src_min < T_neg_hi) {
        // T<0 part: f_r == 1 over its length
        integral += (T_neg_hi - T_src_min);
      }

    // span = T_src_max - T_src_min > 0 here (the equal-potentials case is
    // handled in the branch above). integral comes from the monotone-increasing
    // cumulative-vulnerability curve so it is > 0 over a span>0 interval;
    // forming r_R_H as r_R_H_min * span / integral is one division (vs the old
    // f_r_average = integral/span then r_R_H_min/f_r_average two), and needs no
    // per-layer finiteness guard (any stray NaN/Inf propagates to the post-loop
    // isfinite(E_up) net).
    //
    // The exception, since the cap: with BOTH bounds past 7.3132 MPa the integral
    // is exactly 0 and r_R_H +Inf. E_i is then -0 and duptake_dpsi NaN, which its
    // caller reads as "use central differences".
    const double span = T_src_max - T_src_min;

    // Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity) by the fractional loss of conductivity
    double r_R_H = network_.r_R_H_min[i] * span / integral; // [MPa * s * (mol H2O)^-1]

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    double r_R = r_R_H + network_.r_R_V_sum[i]; // [MPa * s * (mol H2O)^-1]

    // Transpiration is equal to the potential gradient between the root collar
    // and the soil, accounting for gravitational potential. In magnitudes the
    // collar has to pull HARDER than the soil holds, plus enough to lift the
    // water -- hence the subtraction order. E_i < 0 still means the layer gains.
    double E_i = (T_collar - psi_soil[i] - grav_head_z_[i]) / r_R; // [mol H2O / m^2 / s]

    soil_consumption[i] = E_i;
    E_up += E_i;

    }
  }
  // Convert the summed uptake to kg H2O m^-2 s^-1, consistent with the rest of
  // the leaf model and environment. NOTE (review #10): only the aggregate E_up
  // is converted to kg here; the per-layer soil_consumption[i] above is left in
  // mol H2O m^-2 s^-1 and converted downstream in TF24_Strategy::compute_rates.
  // The two siblings therefore carry different units by design.
  E_up = E_up * kg_per_mol_h2o;
  if (!std::isfinite(E_up)) {
    util::stop("E_from_Soil_to_Root_Collar non-finite E_up_; T_collar=" + util::to_string(T_collar) +
               "; max_soil_layer=" + std::to_string(max_soil_layer));
  }
  }
};

}  // namespace phylloptim

#endif
