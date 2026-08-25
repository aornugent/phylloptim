// -*-c++-*-
#ifndef PHYLLOPTIM_ROOTS_HPP_
#define PHYLLOPTIM_ROOTS_HPP_

#include <phylloptim/clamp_sites.hpp>
#include <phylloptim/constants.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/vulnerability.hpp>

#include <odelia/interpolator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace phylloptim {

// The per-layer root hydraulic resistances the supply solve actually consumes.
// Two of the five fields are load-bearing; see set_root_network on why the other
// three are here at all.
// The supply state one uptake reads, at one scalar. References rather than
// values because the double path runs it ~10^3 times per collar solve and owns
// the arrays already; the active path builds its own and hands them the same way.
//
// The layer geometry is not here. Depths belong to the architecture model rather
// than to anything plant differentiates, so they stay where they are.
template <typename T>
struct SupplyAt {
  const std::vector<T>& psi_soil;
  const std::vector<T>& r_R_H_min;   // per layer, fully hydrated
  const std::vector<T>& r_R_V_sum;   // per layer, cumulative from the surface
  const T& root_b;
  const T& root_c;
};

// The same, owned. The active path builds one of these and hands the view; a
// gradient in the supply has this shape too, which is what lets a caller read it
// by name rather than by position.
template <typename T>
struct SupplyValues {
  using value_type = T;
  std::vector<T> psi_soil, r_R_H_min, r_R_V_sum;
  T root_b{}, root_c{};

  SupplyAt<T> at() const {
    return {psi_soil, r_R_H_min, r_R_V_sum, root_b, root_c};
  }

  template <class U>
  SupplyValues<U> rebind_from() const {
    using odelia::util::to_passive;
    SupplyValues<U> out;
    for (const T& v : psi_soil) out.psi_soil.push_back(U(to_passive(v)));
    for (const T& v : r_R_H_min) out.r_R_H_min.push_back(U(to_passive(v)));
    for (const T& v : r_R_V_sum) out.r_R_V_sum.push_back(U(to_passive(v)));
    out.root_b = U(to_passive(root_b));
    out.root_c = U(to_passive(root_c));
    return out;
  }

  // The entries a caller can seed, in one order, so a gradient over them cannot
  // be matched to them by hand. The layer arrays come first, layer-major.
  std::vector<T*> field_ptrs() {
    std::vector<T*> out;
    out.reserve(psi_soil.size() + r_R_H_min.size() + r_R_V_sum.size() + 2);
    for (T& v : psi_soil) out.push_back(&v);
    for (T& v : r_R_H_min) out.push_back(&v);
    for (T& v : r_R_V_sum) out.push_back(&v);
    out.push_back(&root_b);
    out.push_back(&root_c);
    return out;
  }
};

struct RootNetwork {
  // Minimum (fully-hydrated) horizontal, intra-layer soil->root resistance.
  // Divided by the vulnerability-weighted mean conductivity at the operating
  // potential to get the actual horizontal resistance.
  std::vector<double> r_R_H_min;
  // Cumulative vertical, inter-layer resistance from the surface down to layer i.
  std::vector<double> r_R_V_sum;
  // c_r_H is a diagnostic. c_r_V and r_R_V are read by
  // duptake_droot_carbon: the first recovers the layer's carbon, the second is
  // the per-layer vertical resistance the cumulative sum is built from.
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

  // G(m) = int_0^m f_r(s) ds, indexed by magnitude m = -psi, carrying a value, a
  // slope and a curvature at every knot -- all three the closed form's. Lets
  // uptake() obtain the mean conductivity over a potential interval from 2 evals
  // instead of (n+1).
  //
  // ONE table, not two. f_r is dG/dpsi, so the curve is this table's slope rather
  // than a second tabulation of it: two tables agreed at their knots and differed
  // between them, so uptake read G from one polynomial and f_r from another.
  odelia::interpolator::hermite_interpolator<double, 5> root_vuln_integral_from_psi;
  // What the two splines do past their last knot, cached by setup_vulnerability.
  // Both are read through the accessors below, never as bare .eval() calls --
  // see the accessors for why. util::na_value until a curve is built.
  double root_vuln_last_knot_ = util::na_value;
  double root_vuln_integral_limit_ = util::na_value;

  // How often each clamp held. Shared storage, so a consumer copying this object
  // per unit still reads what the copy counted.
  clamp_counter clamps;

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
  // The built curve against the pair that determines it; see Leaf's cache for
  // why the key is (b, c, resolution) and why a handful of entries is enough.
  struct CurveCache {
    double b = 0.0, c = 0.0, resolution = 0.0;
    odelia::interpolator::hermite_interpolator<double, 5> integral;
    double last_knot = 0.0, integral_limit = 0.0;
  };
  static constexpr std::size_t curve_cache_size = 32;
  // Held behind a pointer so a copied MultiLayerRoots SHARES the store rather
  // than duplicating it. A caller that rebuilds this object per cohort per stage
  // -- which a reverse sweep does -- otherwise starts empty every time and the
  // cache never reads back a single entry it wrote.
  std::shared_ptr<std::vector<CurveCache>> curve_cache_ =
      std::make_shared<std::vector<CurveCache>>();

  void setup_vulnerability(double resolution) {
    for (const CurveCache& hit : *curve_cache_) {
      if (hit.b == root_b && hit.c == root_c && hit.resolution == resolution) {
        root_vuln_integral_from_psi = hit.integral;
        root_vuln_last_knot_ = hit.last_knot;
        root_vuln_integral_limit_ = hit.integral_limit;
        return;
      }
    }
    // Three knot vectors from one pass: dG/dpsi IS exp(-(psi/root_b)^root_c) and
    // its own slope needs nothing the same series has not already formed, so every
    // channel of this table is the closed form's and none is inferred.
    // f_r(0) = exp(-pow(0,root_c)) = 1.
    std::vector<double> x_psi_root, y_integral, y_f_r, y_f_r_slope;
    cumulative_vulnerability_integral(root_b, root_c, resolution, x_psi_root,
                                      y_integral, y_f_r, y_f_r_slope);

    root_vuln_integral_from_psi.init(x_psi_root, y_integral, y_f_r, y_f_r_slope);

    // Neither carries an extrapolation setting, because neither ever reads past
    // its knots: root_vuln_at clamps the argument, and root_vuln_integral_at caps
    // the value at the closed-form limit. Past the last knot the interpolant
    // extends the end slope, which for the integral is f_r there -- under a
    // hundredth -- so the cap binds smoothly instead of a quadratic running away.

    // The last knot, NOT vulnerability_psi_max: the knot loop advances by
    // accumulation and stops one step short of psi_max (6.8229 against 6.8918 at
    // the root defaults), so psi_max is itself outside the domain and clamping to
    // it would throw.
    root_vuln_last_knot_ = root_vuln_integral_from_psi.max();
    root_vuln_integral_limit_ =
        cumulative_vulnerability_integral_limit(root_b, root_c);

    if (curve_cache_->size() >= curve_cache_size) {
      curve_cache_->erase(curve_cache_->begin());
    }
    curve_cache_->push_back(CurveCache{root_b, root_c, resolution,
                                       root_vuln_integral_from_psi,
                                       root_vuln_last_knot_,
                                       root_vuln_integral_limit_});
  }

  // f_r at a suction, clamped into the knot domain: the last knot holds about 1%
  // of full conductivity, which is the driest this curve describes, and the clamp
  // is what makes the read in-domain rather than a setting on the interpolant.
  //
  // Operand order is load-bearing: written this way both clamps return psi when
  // psi is NaN, where the reversed forms return the bound. The uptake call site's
  // !isfinite(f_ri) guard is what reads that NaN.
  double root_vuln_at(double psi) const {
    if (psi > root_vuln_last_knot_ || psi < 0.0) {
      clamps.note(CLAMP_ROOT_VULN_ARGUMENT);
    }
    // The integral's own slope, so f_r and G are derivatives of one polynomial. A
    // second table of f_r agreed with this at the knots and nowhere between them.
    return root_vuln_integral_from_psi.slope(
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
    const double raw = root_vuln_integral_from_psi.eval(psi);
    if (raw > root_vuln_integral_limit_) {
      clamps.note(CLAMP_ROOT_VULN_INTEGRAL_CAP);
    }
    return std::min(raw, root_vuln_integral_limit_);
  }

  // dG/dpsi, consistent with root_vuln_integral_at: zero wherever that returns
  // the cap, because there the value no longer depends on psi.
  double root_vuln_integral_deriv_at(double psi) const {
    if (root_vuln_integral_from_psi.eval(psi) >= root_vuln_integral_limit_) {
      return 0.0;
    }
    // The interpolant's own slope, which is the exact derivative of the value the
    // line above just read. The closed form f_r is the derivative of G and NOT of
    // this tabulation of G, so returning it would give a value and a derivative
    // taken from two different functions -- they agree at the knots, where the
    // supplied slope IS f_r, and differ between them by the interpolation error.
    // Measured: a differenced curvature moves from 1e-07 to 1.9e-07 that way.
    return root_vuln_integral_from_psi.slope(psi);
  }

  // d(dG/dpsi)/dpsi -- the integrand's own slope, and NOT a second derivative of
  // the tabulation. Differentiating a cumulative integral in its upper limit
  // leaves the integrand there, so one more derivative leaves the integrand's
  // slope, and exp(-(psi/root_b)^root_c) has that in closed form. A C1
  // interpolant's second derivative is a difference of its data; this is data.
  //
  // Zero wherever the integral is at its cap, and zero at or below the surface,
  // for the two reasons the first derivative is: past the cap the value no longer
  // moves with psi, and below it the integrand is the constant 1.
  double root_vuln_integrand_deriv_at(double psi) const {
    if (!(psi > 0.0)) {
      return 0.0;
    }
    if (root_vuln_integral_from_psi.eval(psi) >= root_vuln_integral_limit_) {
      return 0.0;
    }
    return vulnerability_curve_slope(psi, root_b, root_c);
  }

  // The suction gap below which a layer's mean conductivity is 0/0, and the ONE
  // coincidence a derivative kernel here cannot answer at.
  //
  // The mean conductivity is a span over an integral over that span, so as the
  // collar approaches a layer's potential both vanish together. The value's limit
  // is 1/f_r, which uptake() already computes; the collar derivative's is
  // -f_r'/(2 f_r^2), from the reciprocal of G's divided difference. That limit is
  // NOT written, and the reason is a count rather than a difficulty: no state
  // reaches it -- zero of 540 operating points over one to five layers,
  // radiations from dark to full sun, and potentials from 0.5 to 7 MPa.
  //
  // ⚠️ TWO FURTHER COINCIDENCES WERE REFUSED HERE AND NEITHER IS A KINK FOR A
  // DERIVATIVE, and one of them cost a consumer its whole water channel. A
  // GRAVITY-BALANCED layer has only its NUMERATOR vanish -- the span, the integral
  // and the resistance are all untouched -- so dE_i/dT is 1/r_R there, agreeing
  // with a difference to 2e-12. A moving bound AT ATMOSPHERIC crosses the split the
  // integral is taken in two parts about, and both parts have slope 1 there: below
  // the surface the integrand is 1 and contributes linearly, and f_r(0) = 1 above
  // it. The first is reached by one-layer shade death, where the collar of zero
  // uptake IS the gravity balance; 30 of those 540 points came back not-a-number
  // for it.
  bool at_equal_potentials(double T_collar, double psi_soil_i) const {
    return std::abs(T_collar - psi_soil_i) < 1e-8;
  }

  // dG/d(root_b) at a fixed suction, from the homogeneity the curve already has
  // and with NO rebuild. G integrates exp(-(sigma/root_b)^root_c), which is
  // homogeneous of degree one in (psi, root_b), so Euler's theorem gives
  //
  //   psi * dG/dpsi + root_b * dG/droot_b  ==  G(psi)
  //
  // and rearranging is the whole derivation. Same identity as the stem curve's,
  // on the same builder.
  //
  // It stays right under the cap rather than in spite of it: past the ceiling G
  // is the complete-gamma limit (root_b/root_c)*Gamma(1/root_c), LINEAR in
  // root_b, and dG/dpsi is zero there -- so the identity returns G/root_b, which
  // is that limit's own derivative.
  //
  // root_c has no counterpart: it reshapes the curve rather than scaling it, so
  // its row needs the grid rebuilt, exactly as stem_c's does.
  double root_vuln_integral_droot_b(double psi) const {
    return (root_vuln_integral_at(psi) -
            psi * root_vuln_integral_deriv_at(psi)) / root_b;
  }

  // Which parameter of the root curve a row is taken in. The two reach the supply
  // by the same single route -- the layer's mean conductivity integral, and
  // nothing else -- so one loop serves both and all that differs is which
  // derivative of the curve it accumulates.
  enum class CurveTrait { Position, Steepness };

  // dG/dtheta at a fixed suction. Position is Euler's identity on the TABULATED
  // pair, which stays right past the cap because the limit is homogeneous of
  // degree one in root_b too. Steepness has no such identity -- it reshapes the
  // curve rather than scaling it -- and comes off the incomplete gamma's own
  // shape series, which is the same loop the tabulation is built from.
  //
  // ⚠️ CAPPED, LIKE EVERY OTHER READER HERE, and this was the one that was not.
  // The series holds only on the grid, so past the cap it refused -- by throwing,
  // where the value beside it returns its limit. A deep dry layer under a wetter
  // one is the state that reaches it, which is the arrangement a drying profile
  // produces and whole-plant shutdown does not catch, since that keys off the
  // WETTEST layer. It cost a shaded five-layer plant its whole water channel: the
  // wet bound's row threw, and the row layer fell back to differencing the solve.
  // Past the cap the integral IS the limit, so the limit's own derivative is the
  // answer and no series is involved.
  double root_vuln_integral_dtrait(double psi, CurveTrait trait) const {
    if (trait == CurveTrait::Position) {
      return root_vuln_integral_droot_b(psi);
    }
    if (root_vuln_integral_from_psi.eval(psi) >= root_vuln_integral_limit_) {
      return cumulative_vulnerability_integral_limit_dc(root_b, root_c);
    }
    // ⚠️ A BAND BETWEEN THE TWO DOMAINS, AND IT REFUSES RATHER THAN THROWING.
    // The grid stops one step short of psi_max and the spline extrapolates past it
    // with the slope at its last knot, so between psi_max and the potential where
    // that straight line finally exceeds the limit the value is neither on the grid
    // nor at the cap: the series has no domain there and the model's own integral is
    // a linear extrapolation whose steepness derivative is not one of these
    // expressions. Non-finite is the honest answer -- the caller differences a
    // genuine rebuild instead, which moves the grid as the model does.
    //
    // It threw here until now, and a throw is not a refusal: a deep dry layer under
    // a wetter one put a live plant in this band, and the exception took the whole
    // metric's gradient rather than one row's.
    if (!(std::pow(psi / root_b, root_c) <= vulnerability_x_max())) {
      return util::na_value;
    }
    return cumulative_vulnerability_integral_derivatives_at(psi, root_b, root_c)
        .dc;
  }

  // d(dG/dpsi)/dtheta -- the curve's own trait derivative, elementary in all
  // three of its arguments. Zero wherever the integral is at its cap, for the
  // reason its psi-derivative is zero there: past the cap the value no longer
  // moves, so nothing that moves it can either.
  double root_vuln_integrand_dtrait(double psi, CurveTrait trait) const {
    if (!(psi > 0.0)) {
      return 0.0;   // f_r == 1 there whatever the parameters are
    }
    if (root_vuln_integral_from_psi.eval(psi) >= root_vuln_integral_limit_) {
      return 0.0;
    }
    const double x = std::pow(psi / root_b, root_c);
    const double f = std::exp(-x);
    return trait == CurveTrait::Position
               ? f * x * root_c / root_b
               : -f * x * std::log(psi / root_b);
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
  // A layer's root carbon, recovered from the network that was built out of it.
  // The network holds the carbon split three ways rather than the carbon, and the
  // vertical third is the half to read it back from -- the same inversion
  // `duptake_droot_carbon` performs. NA rather than zero for an unrooted layer:
  // the network is sized to the deepest rooted layer, so a layer below it has no
  // carbon to move AND no slot to move it in, and a zero there would say the
  // outputs are insensitive to carbon that could be put there.
  double root_carbon(int layer) const {
    const std::vector<double>& c = network_.c_r_V;
    return layer >= 0 && layer < int(c.size()) && c[std::size_t(layer)] > 0.0
               ? 3.0 * c[std::size_t(layer)]
               : util::na_value;
  }

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
  // At the soil state and the network this object holds, all in double.
  void uptake(double T_collar, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake(T_collar, held_supply(), soil_consumption, E_up);
  }

  // At a supply the caller owns, which is how the active path reaches it: the
  // resistances come from root carbon, and that arithmetic is the caller's.
  template <typename T>
  void uptake(const T& T_collar, const SupplyAt<T>& at_scalar,
              std::vector<T>& soil_consumption, T& E_up) const {
    uptake_impl<T>(T_collar, at_scalar,
                   root_vuln_integral_soil_.size() ==
                       static_cast<size_t>(max_soil_layer),
                   soil_consumption, E_up);
  }

  // This object's own state, as the view above. No copies: the arrays are here.
  SupplyAt<double> held_supply() const {
    return {psi_soil_, network_.r_R_H_min, network_.r_R_V_sum, root_b, root_c};
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
    const SupplyAt<double> at_scalar{psi_soil, network_.r_R_H_min,
                                     network_.r_R_V_sum, root_b, root_c};
    uptake_impl<double>(T_collar, at_scalar,
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
  // [T_src_min, T_src_max] (root_vuln_integral_from_psi, whose integrand is that
  // table's own slope):
  //   E_i        = (T_collar - T_soil[i] - grav) / r_R,
  //   r_R        = r_R_H_min[i] * span / integral + r_R_V_sum[i],
  //   dspan/dT   = sign_var   (+1 if T_collar is the upper bound, else -1),
  //   dinteg/dT  = sign_var * f_r(T_collar)  for T_collar>0  (else sign_var, f_r==1),
  // and dE_i/dT follows by the quotient rule.
  //
  // CONTRACT: returns NaN when any layer's potential equals T_collar, and only
  // then -- see at_equal_potentials for why that one is 0/0 and why the
  // gravity-balance and at-atmospheric coincidences it used to refuse at as well
  // are not. That refusal is deliberate, not a failure: the caller falls back to a
  // central difference, and an implementation that threw, or returned 0, would
  // silently degrade an acclimation gradient. Any alternative supply path must
  // keep it.
  double duptake_dpsi(double T_collar,
                      const std::vector<double>& psi_soil) const {
    std::vector<double> per_layer;
    const double total_mol = duptake_dpsi_impl(T_collar, psi_soil, per_layer);
    return total_mol * kg_per_mol_h2o;  // match E_up's kg units
  }

  // The same, per layer rather than summed, in kg to match E_up.
  //
  // A stand adjoint needs this and the total will not do: uptake reaches the
  // shared soil one layer at a time, so the operating point's movement has to be
  // priced into each layer's flux separately. Diagonality is not the point here
  // -- this is one column, d(E_i)/d(collar), and every layer has one.
  //
  // Summing this in layer order and multiplying once is what duptake_dpsi does,
  // so the two cannot disagree.
  void duptake_dpsi_by_layer(double T_collar,
                             const std::vector<double>& psi_soil,
                             std::vector<double>& out) const {
    duptake_dpsi_impl(T_collar, psi_soil, out);
    for (double& v : out) {
      v *= kg_per_mol_h2o;
    }
  }

  // d2(E_up)/d(T_collar)2, in kg to match the conductance above. Mirrors
  // duptake_dpsi_impl term for term with one more derivative of each moving part:
  // the span is linear in the collar so its second derivative is zero, and the
  // integral's is the integrand's own slope -- not the tabulation's curvature,
  // which is a difference of the data rather than data.
  //
  // Same refusal as duptake_dpsi, and for the same one reason: NaN where a layer's
  // potential equals the collar, because the mean conductivity is 0/0 there.
  double d2uptake_dpsi2(double T_collar,
                        const std::vector<double>& psi_soil) const {
    double d2_mol = 0.0;
    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;
      const double sign_var = (T_collar > psi_soil[i]) ? 1.0 : -1.0;

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

      // Below the surface the moving bound is in the f_r == 1 part, contributed
      // linearly, so the slope is 1 and its own slope is 0.
      const double fr_at =
          (T_collar > 0.0) ? root_vuln_integral_deriv_at(T_collar) : 1.0;
      const double dfr_at =
          (T_collar > 0.0) ? root_vuln_integrand_deriv_at(T_collar) : 0.0;
      const double dinteg_dT = sign_var * fr_at;
      const double d2integ_dT = sign_var * dfr_at;

      // The mean resistance is r_R_H_min * span / integral, and both derivatives
      // of that quotient are taken here rather than once each: the span's second
      // derivative is zero, which is what leaves only two terms.
      const double r_R_H = network_.r_R_H_min[i] * span / integral;
      const double r_R = r_R_H + network_.r_R_V_sum[i];
      const double quotient_dT =
          (sign_var * integral - span * dinteg_dT) / (integral * integral);
      const double quotient_d2T = -span * d2integ_dT / (integral * integral) -
                                  2.0 * quotient_dT * dinteg_dT / integral;
      const double dr_R_dT = network_.r_R_H_min[i] * quotient_dT;
      const double d2r_R_dT = network_.r_R_H_min[i] * quotient_d2T;

      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      // E_i = num / r_R with num linear in the collar, so the same two terms
      // again -- one for the resistance's curvature, one for the product of the
      // two first derivatives.
      const double dE_i = (r_R - num * dr_R_dT) / (r_R * r_R);
      d2_mol += -num * d2r_R_dT / (r_R * r_R) - 2.0 * dE_i * dr_R_dT / r_R;
    }
    return d2_mol * kg_per_mol_h2o;
  }

private:
  // Per layer in mol, returning the sum in mol: the one loop both public forms
  // read, so a change reaches them together and the total is the sum of the
  // parts by construction rather than by agreement.
  double duptake_dpsi_impl(double T_collar,
                           const std::vector<double>& psi_soil,
                           std::vector<double>& per_layer) const {
    double dEup_dT_mol = 0.0;
    per_layer.assign(psi_soil.size(), 0.0);

    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        per_layer.assign(psi_soil.size(),
                         std::numeric_limits<double>::quiet_NaN());
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
      // cumulative curve that produced `integral` (root_vuln_integral_deriv_at),
      // and the two accessors differ only in how they are bounded past the knots:
      // the conductivity read clamps its argument to the last knot, the integral is
      // capped at G(inf) -- so beyond the domain only the integral's
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
      const double dE_i = (dnum_dT * r_R - num * dr_R_dT) / (r_R * r_R);
      per_layer[std::size_t(i)] = dE_i;
      dEup_dT_mol += dE_i;
    }

    return dEup_dT_mol;
  }

  // d(E_i)/d(root_b) at a fixed collar, mirroring duptake_dpsi_impl term for
  // term. root_b enters through ONE quantity -- the layer's mean conductivity
  // integral -- so span, the numerator and both resistances are constants here
  // and the quotient rule has a single moving part.
  //
  // Same kink contract as duptake_dpsi_impl, and for the same reason: the
  // general branch is not valid across them, and a caller mixing an analytic row
  // with a missing one is worse off than one told the whole block is undefined.
  double duptake_droot_curve_impl(double T_collar,
                                  const std::vector<double>& psi_soil,
                                  CurveTrait trait,
                                  std::vector<double>& per_layer) const {
    double dEup_db_mol = 0.0;
    per_layer.assign(psi_soil.size(), 0.0);

    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        per_layer.assign(psi_soil.size(),
                         std::numeric_limits<double>::quiet_NaN());
        return std::numeric_limits<double>::quiet_NaN();
      }

      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;

      // integral, replicated bit-for-bit from uptake_impl.
      const double T_pos_lo = std::max(T_src_min, 0.0);
      const double T_neg_hi = std::min(T_src_max, 0.0);
      double integral = 0.0;
      double dinteg_db = 0.0;
      if (T_pos_lo < T_src_max) {
        integral += root_vuln_integral_at(T_src_max) -
                    root_vuln_integral_at(T_pos_lo);
        dinteg_db += root_vuln_integral_dtrait(T_src_max, trait) -
                     root_vuln_integral_dtrait(T_pos_lo, trait);
      }
      if (T_src_min < T_neg_hi) {
        // The above-atmospheric part contributes its width, with f_r == 1
        // throughout. No curve, so neither parameter of it.
        integral += (T_neg_hi - T_src_min);
      }

      const double r_R_H = network_.r_R_H_min[i] * span / integral;
      const double r_R = r_R_H + network_.r_R_V_sum[i];
      // Only the integral moves, and it is in the denominator.
      const double dr_R_db =
          -network_.r_R_H_min[i] * span * dinteg_db / (integral * integral);

      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      // E_i = num / r_R with num constant in root_b.
      const double dE_i = -num * dr_R_db / (r_R * r_R);
      per_layer[std::size_t(i)] = dE_i;
      dEup_db_mol += dE_i;
    }

    return dEup_db_mol;
  }

public:

  // d(E_up)/d(a root curve parameter) at a fixed collar, summed and in kg to
  // match E_up. Closed form in both: the position by Euler's identity, because
  // the curve is scaled rather than reshaped by it, and the steepness by the
  // incomplete gamma's shape series.
  double duptake_droot_curve(double T_collar,
                             const std::vector<double>& psi_soil,
                             CurveTrait trait) const {
    std::vector<double> per_layer;
    return duptake_droot_curve_impl(T_collar, psi_soil, trait, per_layer) *
           kg_per_mol_h2o;
  }

  // The same in either parameter, per layer and in kg. The steepness used to have
  // no row here at all and its perturbation rebuilt the grid; it reaches the
  // supply exactly where the position does.
  void duptake_droot_curve_by_layer(double T_collar,
                                    const std::vector<double>& psi_soil,
                                    CurveTrait trait,
                                    std::vector<double>& out) const {
    duptake_droot_curve_impl(T_collar, psi_soil, trait, out);
    for (double& v : out) {
      v *= kg_per_mol_h2o;
    }
  }

  // d(E_i)/d(psi_soil[i]) for every rooted layer, in kg to match E_up.
  //
  // DIAGONAL, and that is a property of the model rather than an approximation:
  // layer i's flux reads its own potential and the collar, and no other layer's.
  // So this vector is the whole soil Jacobian of the supply, not its diagonal
  // part.
  //
  // Same quotient rule as duptake_dpsi with the moving bound at the soil end
  // instead of the collar. Two factors flip and nothing else does: the span
  // shrinks as the soil approaches the collar (dspan/dpsi = -dspan/dT_collar),
  // and the vulnerability integral's moving endpoint is psi_i rather than
  // T_collar, so its derivative is read AT psi_i and carries the same minus.
  // Writing it as -duptake_dpsi would be wrong: the two endpoints sit at
  // different points on a non-linear curve, so the integral terms do not cancel.
  //
  // Layers past max_soil_layer carry no roots and are written zero rather than
  // left alone -- the caller's buffer is reused across solves (hazard 8).
  //
  // NaN contract is duptake_dpsi's: where a bound meets a layer the whole vector is NaN,
  // because a caller that used some layers and not others would be mixing an
  // analytic row with a missing one.
  void duptake_dpsi_soil(double T_collar, const std::vector<double>& psi_soil,
                         std::vector<double>& out) const {
    out.assign(psi_soil.size(), 0.0);

    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        out.assign(psi_soil.size(), std::numeric_limits<double>::quiet_NaN());
        return;
      }

      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;
      // dspan/dT_collar was +1 when the collar pulls; the soil end moves the
      // other way.
      const double sign_var = (T_collar > psi_soil[i]) ? 1.0 : -1.0;
      const double dspan_dpsi = -sign_var;

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

      // Read at psi_i, not at the collar: same curve, other endpoint.
      const double fr_at =
          (psi_soil[i] > 0.0) ? root_vuln_integral_deriv_at(psi_soil[i]) : 1.0;
      const double dinteg_dpsi = dspan_dpsi * fr_at;

      const double r_R_H = network_.r_R_H_min[i] * span / integral;
      const double r_R = r_R_H + network_.r_R_V_sum[i];
      const double dr_R_dpsi = network_.r_R_H_min[i] *
                               (dspan_dpsi * integral - span * dinteg_dpsi) /
                               (integral * integral);

      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      const double dnum_dpsi = -1.0;
      out[std::size_t(i)] =
          ((dnum_dpsi * r_R - num * dr_R_dpsi) / (r_R * r_R)) * kg_per_mol_h2o;
    }
  }

  // d(E_i)/d(root carbon in layer a), and the same for d(E_i)/d(T_collar), both
  // in kg to match the accessors above. Two blocks, layer by carbon-layer.
  //
  // The carbon reaches a flux only through the two resistances the architecture
  // model builds from it, and both are proportional to 1/rc: the horizontal one
  // is layer i's own, the vertical one is summed over every layer at or above i.
  // So layer a's carbon reaches every layer i >= a and the blocks are LOWER
  // TRIANGULAR, not diagonal. An implementation that assumes diagonality gets
  // the shallow layers right and loses the deep ones, which on a drying profile
  // is the half that carries the flux.
  //
  // Why this cannot be recovered from duptake_dpsi: there the cumulative
  // vertical resistance drops out entirely, because it has no collar dependence,
  // so that column carries the horizontal term alone. The carbon direction
  // reaches both.
  //
  // NaN contract is duptake_dpsi's -- where the collar meets a layer both blocks are
  // NaN, because a caller using some layers and not others would be mixing an
  // analytic row with a missing one. A layer with no carbon has no resistance
  // and contributes nothing rather than dividing by its carbon.
  void duptake_droot_carbon(double T_collar,
                            const std::vector<double>& psi_soil,
                            std::vector<std::vector<double>>& dE_drc,
                            std::vector<std::vector<double>>& dD_drc) const {
    const std::size_t n = psi_soil.size();
    dE_drc.assign(n, std::vector<double>(n, 0.0));
    dD_drc.assign(n, std::vector<double>(n, 0.0));
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        dE_drc.assign(n, std::vector<double>(n, nan));
        dD_drc.assign(n, std::vector<double>(n, nan));
        return;
      }

      // Everything down to g is duptake_dpsi_impl's loop, unchanged: the same
      // span, the same integral, the same moving-bound derivative.
      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;
      const double sign_var = (T_collar > psi_soil[i]) ? 1.0 : -1.0;

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
      const double fr_at =
          (T_collar > 0.0) ? root_vuln_integral_deriv_at(T_collar) : 1.0;
      const double dinteg_dT = sign_var * fr_at;

      // A and B carry the whole of the carbon dependence; f and g carry the
      // whole of the collar dependence. That separation is what makes the rest
      // a quotient rule rather than a new model.
      const double A = network_.r_R_H_min[i];
      const double f = span / integral;
      const double B = network_.r_R_V_sum[i];
      const double r_R = A * f + B;
      const double g = (sign_var * integral - span * dinteg_dT) /
                       (integral * integral);
      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      const double E_i = num / r_R;
      const double N = r_R - num * A * g;   // the numerator of dE_i/dT_collar
      const double r3 = r_R * r_R * r_R;
      const double dD_dA = ((f - num * g) * r_R - 2.0 * N * f) / r3;
      const double dD_dB = (-r_R + 2.0 * num * A * g) / r3;

      for (int a = 0; a < max_soil_layer; a++) {
        // The network stores the split carbon rather than the carbon; either
        // half recovers it, and c_r_V is the one the vertical term is built on.
        const double rc_a = 3.0 * network_.c_r_V[std::size_t(a)];
        if (!(rc_a > 0.0)) {
          continue;   // no roots in that layer, so no route from its carbon
        }
        const double dA = (i == a) ? -A / rc_a : 0.0;
        const double dB = (a <= i) ? -network_.r_R_V[std::size_t(a)] / rc_a
                                   : 0.0;
        const double dr = f * dA + dB;
        dE_drc[std::size_t(i)][std::size_t(a)] =
            -(E_i / r_R) * dr * kg_per_mol_h2o;
        dD_drc[std::size_t(i)][std::size_t(a)] =
            (dD_dA * dA + dD_dB * dB) * kg_per_mol_h2o;
      }
    }
  }

  // d2(E_i)/d(T_collar) d(theta) for either parameter of the root curve, per
  // layer and in kg.
  //
  // The curve reaches a flux through the layer's mean conductivity integral and
  // nothing else, so this is one more application of the quotient rule to
  // duptake_droot_curve_impl -- with the integral's moving bound now carrying the
  // curve's own trait derivative as well as its slope.
  //
  // NaN contract is duptake_dpsi's, over the union of both first derivatives'
  // kinks: a second derivative needs both of them off theirs.
  void d2uptake_dpsi_droot_curve(double T_collar,
                                 const std::vector<double>& psi_soil,
                                 CurveTrait trait,
                                 std::vector<double>& out) const {
    out.assign(psi_soil.size(), 0.0);

    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        out.assign(psi_soil.size(), std::numeric_limits<double>::quiet_NaN());
        return;
      }

      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;
      const double sign_var = (T_collar > psi_soil[i]) ? 1.0 : -1.0;

      const double T_pos_lo = std::max(T_src_min, 0.0);
      const double T_neg_hi = std::min(T_src_max, 0.0);
      double integral = 0.0;
      double integral_t = 0.0;   // dI/dtheta
      if (T_pos_lo < T_src_max) {
        integral += root_vuln_integral_at(T_src_max) -
                    root_vuln_integral_at(T_pos_lo);
        integral_t += root_vuln_integral_dtrait(T_src_max, trait) -
                      root_vuln_integral_dtrait(T_pos_lo, trait);
      }
      if (T_src_min < T_neg_hi) {
        integral += (T_neg_hi - T_src_min);
      }

      // The integrand at the moving bound, and its own trait derivative there.
      const double fr =
          (T_collar > 0.0) ? root_vuln_integral_deriv_at(T_collar) : 1.0;
      const double fr_t = (T_collar > 0.0)
                              ? root_vuln_integrand_dtrait(T_collar, trait)
                              : 0.0;
      const double dI_dT = sign_var * fr;
      const double dI_t_dT = sign_var * fr_t;

      const double H = network_.r_R_H_min[i];
      const double r_R = H * span / integral + network_.r_R_V_sum[i];
      const double dr_dT =
          H * (sign_var * integral - span * dI_dT) / (integral * integral);
      // dr/dtheta, and its own collar derivative.
      const double dr_dt = -H * span * integral_t / (integral * integral);
      const double d2r =
          -H * (sign_var * integral_t + span * dI_t_dT) / (integral * integral) +
          2.0 * H * span * integral_t * dI_dT / (integral * integral * integral);

      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      // E_i = num / r_R, so dE/dtheta = -num dr/dtheta / r_R^2, and one more
      // collar derivative of THAT is the three terms below.
      out[std::size_t(i)] =
          (-(dr_dt + num * d2r) / (r_R * r_R) +
           2.0 * num * dr_dt * dr_dT / (r_R * r_R * r_R)) * kg_per_mol_h2o;
    }
  }

  // d2(E_i)/d(T_collar) d(psi_soil[i]), diagonal for duptake_dpsi_soil's reason.
  //
  // A stand adjoint needs this and the forward model does not, so it is worth
  // saying what it is for: the marginal profit reads the state only through
  // total uptake and through uptake's own collar sensitivity, so a row of the
  // mixed second derivative of profit is a pair of scalars times this vector and
  // d(E_i)/d(psi_soil[i]). Differencing it instead costs 2(L+1) leaf
  // re-evaluations per cohort per stage, which is about eight collar solves.
  //
  // The two endpoints of the vulnerability integral are independent, so
  // d2(integral)/dT dpsi is zero and so is d2(span)/dT dpsi; that is what keeps
  // this to one more application of the quotient rule rather than a new object.
  //
  // NaN contract is duptake_dpsi's, over the union of both first derivatives'
  // kinks: a second derivative needs both endpoints off theirs.
  void d2uptake_dpsi_dpsi_soil(double T_collar,
                               const std::vector<double>& psi_soil,
                               std::vector<double>& out) const {
    out.assign(psi_soil.size(), 0.0);

    for (int i = 0; i < max_soil_layer; i++) {
      if (at_equal_potentials(T_collar, psi_soil[i])) {
        out.assign(psi_soil.size(), std::numeric_limits<double>::quiet_NaN());
        return;
      }

      const double T_src_min = std::min(psi_soil[i], T_collar);
      const double T_src_max = std::max(psi_soil[i], T_collar);
      const double span = T_src_max - T_src_min;
      const double sign_var = (T_collar > psi_soil[i]) ? 1.0 : -1.0;

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

      // The integrand at each moving endpoint, each read at its own end.
      const double fr_T =
          (T_collar > 0.0) ? root_vuln_integral_deriv_at(T_collar) : 1.0;
      const double fr_psi =
          (psi_soil[i] > 0.0) ? root_vuln_integral_deriv_at(psi_soil[i]) : 1.0;

      const double H = network_.r_R_H_min[i];
      const double dinteg_dT = sign_var * fr_T;
      const double dinteg_dpsi = -sign_var * fr_psi;

      const double r_R = H * span / integral + network_.r_R_V_sum[i];
      const double dr_dT =
          H * (sign_var * integral - span * dinteg_dT) / (integral * integral);
      // A = dspan/dpsi * integral - span * dinteg/dpsi, so dr/dpsi = H A / I^2.
      const double A = -sign_var * integral - span * dinteg_dpsi;
      const double dr_dpsi = H * A / (integral * integral);
      // dA/dT, with both second derivatives of the endpoints vanishing.
      const double dA_dT = fr_psi - fr_T;
      const double d2r_dT_dpsi =
          H * (dA_dT * integral - 2.0 * A * dinteg_dT) /
          (integral * integral * integral);

      const double num = T_collar - psi_soil[i] - grav_head_z_[i];
      // E = num / r_R; dE/dpsi = N / r^2 with N = -r - num dr/dpsi.
      const double N = -r_R - num * dr_dpsi;
      const double dN_dT = -dr_dT - dr_dpsi - num * d2r_dT_dpsi;
      out[std::size_t(i)] = ((dN_dT * r_R - 2.0 * N * dr_dT) /
                             (r_R * r_R * r_R)) * kg_per_mol_h2o;
    }
  }

private:
  // Total water drawn from all layers to the collar. Writes E_up (kg H2O m^-2
  // leaf s^-1) and soil_consumption[i] (mol H2O m^-2 leaf s^-1, note the unit
  // split); a negative E_i in a layer means that layer is *gaining* water
  // (hydraulic redistribution).
  //
  // Implementation decisions:
  //   * f_r and its running integral are read from one pre-computed table
  //     (root_vuln_integral_from_psi, f_r being its slope) instead of repeatedly
  //     evaluating exp(-(psi/b)^c); see setup_vulnerability.
  //   * Two special cases are handled exactly to avoid division/round-off
  //     issues: (a) collar potential equals layer potential, and (b) the
  //     gradient exactly balances gravity (E_i = 0).
  //   * The isfinite() guards are present because this is called from within
  //     nested root-finders where bad brackets can produce NaNs; they fail fast
  //     with diagnostic context rather than propagating NaN.
  // ⚠️ ONE QUADRATURE AT TWO SCALARS, and that is the point of the template. The
  // forward solve runs it at double ~10^3 times per collar solve; the derivative
  // path runs it once at an active scalar. A second implementation for the second
  // scalar is a place the two can disagree, which is the one thing this boundary
  // cannot afford.
  //
  // At double it is the arithmetic it always was, to the bit: the min/max are
  // written as the selects they compile to, and the cumulative-integral cache is
  // taken only there -- at an active scalar the argument carries a derivative and
  // an equality test on it is a test on the value alone.
  //
  // The soil state and the network are held. The collar is what moves, which is
  // what the marginal profit reads.
  template <typename T>
  void uptake_impl(const T& T_collar, const SupplyAt<T>& at_scalar,
                   bool use_integral_cache,
                   std::vector<T>& soil_consumption, T& E_up) const {
    using odelia::util::to_passive;
    const std::vector<T>& psi_soil = at_scalar.psi_soil;
    const double collar_at = to_passive(T_collar);
    const double root_b0 = to_passive(at_scalar.root_b);
    const double root_c0 = to_passive(at_scalar.root_c);

    if (!std::isfinite(collar_at)) {
      util::stop("E_from_Soil_to_Root_Collar invalid input; T_collar=" + util::to_string(collar_at));
    }

    E_up = T(0.0);

    // Cumulative-integral spline caching (bit-identical fast path). The only two
    // arguments ever passed to root_vuln_integral_from_psi in the loop below are
    // T_src_max and T_pos_lo, each of which resolves to exactly one of
    // {psi_soil[i], T_collar, 0}. T_collar is constant across all layers (compute
    // once), and psi_soil[i] is constant across the whole solve (precomputed in
    // begin_solve).
    const double G_at_T_collar =
        use_integral_cache ? root_vuln_integral_at(collar_at) : 0.0;

    // GUARD POLICY (the per-layer isfinite/stop guards here were added while
    // debugging the #485 drought-NaN, now fixed at source by the soil residual-
    // moisture floor). Most were defensive and redundant, so they have been
    // removed from this hot loop; the remaining two are load-bearing:
    //   * the equal-potentials f_ri <= 0 check below: it USED to be the only
    //     thing standing between a deep-drought layer and the negative
    //     conductivity the curve's extrapolant produced past its
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
    // the SMALLER suction, where the signed convention took a minimum. Written as
    // the select std::min compiles to, so the collar can carry a derivative
    // through whichever end it is.
    const T& psi_i = psi_soil[std::size_t(i)];
    const double soil_at = to_passive(psi_i);
    const T T_src_min = (collar_at < soil_at) ? T_collar : psi_i;

    // The drier end: the LARGER suction.
    const T T_src_max = (soil_at < collar_at) ? T_collar : psi_i;

     // If root collar soil water potential equals the soil water potential in a given layer
    if(std::abs(collar_at - soil_at) < 1e-8){

      // Fraction of conductance in roots in a given layer at the driest suction
      // (which here equals the root collar's).
      // f_r is the cumulative table's own slope, read through root_vuln_at so a
      // layer drier than the grid gets the last knot's conductivity rather than an
      // extrapolated (eventually negative) one.
      const double f_ri = root_vuln_at(to_passive(T_src_max));
      if (!std::isfinite(f_ri) || f_ri <= 0.0) {
        util::stop("E_from_Soil_to_Root_Collar invalid f_ri; layer=" + std::to_string(i) +
                   "; f_ri=" + util::to_string(f_ri) +
                   "; T_src_max=" + util::to_string(to_passive(T_src_max)) +
                   "; T_collar=" + util::to_string(collar_at));
      }

      // Fraction of conductance in roots in a given layer at the driest suction
      const T r_R_H = at_scalar.r_R_H_min[std::size_t(i)] / f_ri; // [MPa * s * (mol H2O)^-1]

      // Total root resistance (horizantal plus vertical)
      const T r_R = r_R_H + at_scalar.r_R_V_sum[std::size_t(i)];

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains
      // water). The collar does not appear, so this branch carries the
      // resistances' derivative and none of the collar's -- which is what the
      // model says about it.
      const T E_i = T(-grav_head_z_[i]) / r_R;

      soil_consumption[i] = E_i;
      E_up += E_i;

    }
    else if(std::abs((collar_at - soil_at) - grav_head_z_[i]) < 1e-8){
      // If pressure difference perfectly balances gravity transpiration is equal to zero
      const T E_i = T(0.0); // [mol H2O / m^2 / s]

      soil_consumption[i] = E_i;

      E_up += E_i;

    } else{

      // Mean fractional root conductivity over the suction interval
      // [T_src_min, T_src_max], i.e. (1/(b-a)) * integral_a^b f_r dT.
      // Computed from the pre-integrated curve G(m) = integral_0^m f_r(s) ds
      // (root_vuln_integral_from_psi, indexed by the suction magnitude) with 2
      // evals instead of the old (n+1)-point sample mean. The interval is split
      // at T = 0: for T < 0 (an above-atmospheric potential) vulnerability is 1.
      const T T_pos_lo = (to_passive(T_src_min) < 0.0) ? T(0.0) : T_src_min;
      const T T_neg_hi = (0.0 < to_passive(T_src_max)) ? T(0.0) : T_src_max;

      // Memoised cumulative-integral lookup. Returns the exact same double the
      // spline would (same input -> same output); the comparisons select the
      // precomputed value because T_src_max / T_pos_lo are bit-for-bit equal to
      // one of the cached arguments in the common (T>=0) case.
      // The cumulative curve, carrying its own integrand as the slope: the value
      // is the table's, because the solve ran on the table, and the derivative is
      // the curve's, which is what root_vuln_integral_deriv_at is.
      // The cumulative root curve, carrying every slope the same curve gives:
      // its own integrand in the query, the integrand's slope for the second
      // order the collar's condition takes, and the two trait derivatives with
      // their cross terms in the query. Five reads of one curve, none of them a
      // hand-kept mirror of it.
      auto G_integral = [&](const T& arg) -> T {
        const double q = to_passive(arg);
        if constexpr (std::is_same_v<T, double>) {
          if (use_integral_cache) {
            if (q == collar_at) return G_at_T_collar;
            if (q == soil_at) return root_vuln_integral_soil_[i];
          }
          return root_vuln_integral_at(q);
        } else {
          const T step = arg - T(q);
          const T db = at_scalar.root_b - T(root_b0);
          const T dc = at_scalar.root_c - T(root_c0);
          const T slope =
              T(root_vuln_integral_deriv_at(q)) +
              T(root_vuln_integrand_dtrait(q, CurveTrait::Position)) * db +
              T(root_vuln_integrand_dtrait(q, CurveTrait::Steepness)) * dc;
          return T(root_vuln_integral_at(q)) + slope * step +
                 T(0.5 * root_vuln_integrand_deriv_at(q)) * step * step +
                 T(root_vuln_integral_dtrait(q, CurveTrait::Position)) * db +
                 T(root_vuln_integral_dtrait(q, CurveTrait::Steepness)) * dc;
        }
      };

      T integral = T(0.0);
      if (to_passive(T_pos_lo) < to_passive(T_src_max)) {
        // T>=0 part: suction runs from T_pos_lo up to T_src_max
        integral += G_integral(T_src_max) - G_integral(T_pos_lo);
      }
      if (to_passive(T_src_min) < to_passive(T_neg_hi)) {
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
    const T span = T_src_max - T_src_min;

    // Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity) by the fractional loss of conductivity
    const T r_R_H =
        at_scalar.r_R_H_min[std::size_t(i)] * span / integral; // [MPa * s * (mol H2O)^-1]

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    const T r_R = r_R_H + at_scalar.r_R_V_sum[std::size_t(i)]; // [MPa * s * (mol H2O)^-1]

    // Transpiration is equal to the potential gradient between the root collar
    // and the soil, accounting for gravitational potential. In magnitudes the
    // collar has to pull HARDER than the soil holds, plus enough to lift the
    // water -- hence the subtraction order. E_i < 0 still means the layer gains.
    const T E_i = (T_collar - psi_i - T(grav_head_z_[i])) / r_R; // [mol H2O / m^2 / s]

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
  if (!std::isfinite(to_passive(E_up))) {
    util::stop("E_from_Soil_to_Root_Collar non-finite E_up_; T_collar=" + util::to_string(collar_at) +
               "; max_soil_layer=" + std::to_string(max_soil_layer));
  }
  }
};

}  // namespace phylloptim

#endif
