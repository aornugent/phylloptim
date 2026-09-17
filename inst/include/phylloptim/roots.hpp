// -*-c++-*-
#ifndef PHYLLOPTIM_ROOTS_HPP_
#define PHYLLOPTIM_ROOTS_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/vulnerability.hpp>
#include <phylloptim/closed_form_rows.hpp>
#include <phylloptim/clamp_sites.hpp>

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
//   r_R_V[i]     = beta_R_V * dz[i]^2 / c_r_v (vertical; dz[i]^2 because
//                                           vertical conductivity scales with
//                                           root cross-sectional area, and the
//                                           segment spanning layer i is dz[i]
//                                           long)
//   network_.r_R_V_sum[i] = cumulative vertical resistance from the surface to layer i.
//
// ⚠️ `dz` IS PER LAYER, AND THAT IS LOAD-BEARING RATHER THAN GENERALITY FOR ITS
// OWN SAKE. It was one scalar until #626, which was correct only for a profile of
// equal layers and silently wrong for any other. The property that fails is
// DISCRETISATION INVARIANCE: with a root density uniform over depth, the
// layer-integrated carbon in layer i is c_r_V[i] = (C/3)*dz[i]/D, so the
// per-layer form gives r_R_V[i] = 3*beta_R_V*D*dz[i]/C and the total
// sum(r_R_V) = 3*beta_R_V*D^2/C -- independent of how the column is sliced, as it
// must be, since slicing is a numerical choice and not a property of the plant. A
// single scalar dz = D/n breaks that: total resistance picks up a factor
// (D/n)^2 * sum(1/dz[i]) / D^2, which is 1 only for equal layers. Measured with a
// 2 cm surface layer over 1.5 m in five layers, it inflates total vertical root
// resistance 3.68x and throttles uptake, with nothing anywhere reporting it.
//
// So the reduction test is a total, not a per-layer value, and it is the one
// worth keeping: sum(r_R_V) must not depend on the profile.
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
    const std::vector<double>& root_carbon_per_layer,
    const std::vector<double>& dz,
    double beta_R_H, double beta_R_V, RootNetwork& out) {
  const size_t n_layers = root_carbon_per_layer.size();

  if (dz.size() != n_layers) {
    util::stop("root_network_from_carbon: dz and root_carbon_per_layer must "
               "have the same number of elements");
  }
  for (size_t i = 0; i < n_layers; ++i) {
    // Checked for every layer rather than only the rooted ones, because a
    // caller's profile is wrong as a whole if any width is: the layers below the
    // rooting depth are the ones a shallow plant does not reach *today* and will
    // reach as it grows.
    if (!util::is_finite(dz[i]) || dz[i] <= 0.0) {
      util::stop("root_network_from_carbon: dz must be finite and positive in "
                 "every layer; layer " + util::to_string(i + 1) +
                 " is " + util::format_double(dz[i]));
    }
  }

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
    // The vertical conductivity is likely linearly proportional to the root area
    // projected onto the horizontal plane, hence dz^2 -- and it is THIS layer's
    // thickness, because the vertical root segment being resisted is the one
    // spanning this layer. See the discretisation-invariance note above.
    //
    // ⚠️ SQUARE FIRST, INTO ITS OWN VARIABLE. `beta_R_V * dz[i] * dz[i] / c_r_v`
    // associates left to right, i.e. (beta_R_V*dz)*dz, which is NOT the
    // beta_R_V*(dz*dz) this has always computed and would move the golden file by
    // a last bit for no reason.
    const double dz_sq = dz[i] * dz[i];
    out.r_R_V[i] = beta_R_V * dz_sq / c_r_v;
    vertical_resistance_sum += out.r_R_V[i];
    out.r_R_V_sum[i] = vertical_resistance_sum;
  }
}

// Same, returning a fresh network. Convenience for tests and standalone callers.
// The same map at any scalar, writing the two series the supply view reads.
//
// An OVERLOAD rather than a second name: the arithmetic below is upstream's term
// for term, including squaring dz into its own variable first, so the double
// instantiation is bit-identical to the RootNetwork one. It exists because the
// carbon -> resistance chain belongs to whoever owns the architecture model --
// plant -- and has to land on the same tape as everything else rather than
// crossing as a hand-written d(uptake)/d(carbon).
//
// It writes only r_R_H_min and r_R_V_sum. c_r_V, c_r_H and the per-layer r_R_V
// are intermediates the RootNetwork keeps for the R boundary and nothing on a
// differentiated path reads.
template <class T>
inline void root_network_from_carbon(
    const std::vector<T>& root_carbon_per_layer, const std::vector<double>& dz,
    double beta_R_H, double beta_R_V, std::vector<T>& r_R_H_min,
    std::vector<T>& r_R_V_sum) {
  using odelia::util::to_passive;
  const std::size_t n_layers = root_carbon_per_layer.size();
  if (dz.size() != n_layers) {
    util::stop("root_network_from_carbon: dz and root_carbon_per_layer must "
               "have the same number of elements");
  }
  for (std::size_t i = 0; i < n_layers; ++i) {
    if (!util::is_finite(dz[i]) || dz[i] <= 0.0) {
      util::stop("root_network_from_carbon: dz must be finite and positive in "
                 "every layer; layer " + util::to_string(i + 1) +
                 " is " + util::format_double(dz[i]));
    }
  }

  // The deepest layer carrying carbon. ⚠️ A SELECTOR, so it reads passive: which
  // layers are rooted is piecewise constant in the carbon, and differentiating
  // the comparison would manufacture a jump where a layer starts being reached.
  int max_soil_layer = 0;
  for (std::size_t i = 0; i < n_layers; ++i) {
    if (to_passive(root_carbon_per_layer[i]) != 0.0) {
      max_soil_layer = static_cast<int>(i) + 1;
    }
  }
  r_R_H_min.assign(static_cast<std::size_t>(max_soil_layer), T(0.0));
  r_R_V_sum.assign(static_cast<std::size_t>(max_soil_layer), T(0.0));

  T vertical_resistance_sum = T(0.0);
  for (int i = 0; i < max_soil_layer; ++i) {
    const T& root_mass = root_carbon_per_layer[std::size_t(i)];
    const double at = to_passive(root_mass);
    if (at < 0.0) {
      util::stop("Root mass lower than 0");
    }
    // A layer with no carbon has no resistance to give and none to accumulate;
    // the running sum passes through it unchanged.
    if (at == 0.0) {
      r_R_V_sum[std::size_t(i)] = vertical_resistance_sum;
      continue;
    }
    const T c_r_v = root_mass / 3.0;
    const T c_r_h = root_mass * 2.0 / 3.0;
    r_R_H_min[std::size_t(i)] = T(beta_R_H) / c_r_h;
    // ⚠️ SQUARE FIRST, INTO ITS OWN VARIABLE -- see the double form above.
    const double dz_sq = dz[std::size_t(i)] * dz[std::size_t(i)];
    vertical_resistance_sum += T(beta_R_V) * T(dz_sq) / c_r_v;
    r_R_V_sum[std::size_t(i)] = vertical_resistance_sum;
  }
}

inline RootNetwork root_network_from_carbon(
    const std::vector<double>& root_carbon_per_layer,
    const std::vector<double>& dz,
    double beta_R_H, double beta_R_V) {
  RootNetwork out;
  root_network_from_carbon(root_carbon_per_layer, dz, beta_R_H, beta_R_V, out);
  return out;
}

// Per-layer thicknesses implied by a cumulative soil-depth profile: the
// differences between consecutive boundaries, with an implicit 0 at the surface.
//
// FOR A CALLER WHO HAS A PROFILE AND NOT THE WIDTHS. It is a convenience, and
// deliberately no longer billed as "the shared definition of dz" the way its
// scalar predecessor was: since #626 there is no derivation to share, because the
// widths ARE the geometry a caller states and everything else follows from them.
// A caller who holds its own widths -- plant's TF24_Environment does -- must pass
// those, NOT re-derive them through here.
//
// ⚠️ THE REASON THAT IS NOT PEDANTRY: differencing is not bit-exact. Over 180
// (total depth, n) pairs, `diff` of a profile built as `(i+1)*depth/n` differs
// from `depth/n` in the last bit in 122 of them -- the package defaults among
// them. So re-deriving widths a caller already holds silently perturbs every
// vertical resistance, and through the solve that reaches the reported outputs.
// Exact whenever the boundaries are exactly representable, which is why the
// golden grid (1 m layers) is unaffected.
inline std::vector<double>
layer_thicknesses(const std::vector<double>& soil_depth) {
  if (soil_depth.empty()) {
    // Guarded because this is a public entry point a caller reaches for directly,
    // where the same expression inside set_soil_state was only ever reachable
    // after set_physiology had validated the profile against psi_soil.
    util::stop("layer_thicknesses: soil_depth must have at least one layer");
  }
  std::vector<double> out(soil_depth.size());
  double previous = 0.0;
  for (size_t i = 0; i < soil_depth.size(); ++i) {
    // Strictly increasing, so every width comes out positive. Checked here rather
    // than left to root_network_from_carbon's positivity guard because the
    // diagnosis differs: there the widths are wrong, here the PROFILE is not a
    // profile.
    //
    // ⚠️ IT DOES NOT CATCH WIDTHS PASSED WHERE BOUNDARIES WERE WANTED, and do not
    // read it as though it did. A width vector that happens to increase with
    // depth -- which a graded profile thinning towards the surface does -- passes
    // every check here and yields different, entirely plausible widths. Only the
    // caller knows which of the two it is holding. Asserted in
    // test_layer_thicknesses so the limit is on the record.
    if (!util::is_finite(soil_depth[i]) || !(soil_depth[i] > previous)) {
      util::stop("layer_thicknesses: soil_depth must be finite and strictly "
                 "increasing (cumulative depth to the bottom of each layer); "
                 "layer " + util::to_string(i + 1) + " is " +
                 util::format_double(soil_depth[i]) + " after " +
                 util::format_double(previous));
    }
    out[i] = soil_depth[i] - previous;
    previous = soil_depth[i];
  }
  return out;
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
// into them after crown integration, in TF24_Strategy::net_mass_production_dt.
// They are plant's buffers, not this object's state. Note the deliberate unit split:
// E_up is kg H2O m^-2 s^-1, soil_consumption[i] is mol, converted downstream.

// The active inputs the uptake path reads, as a view.
//
// A VIEW rather than an owned bundle, because there are two owners: the double
// path assembles one from this object's own members (held_supply below), and the
// active path from vectors its caller owns. One uptake_impl then serves both,
// where an owned type would force the double path to copy its own state.
//
// ⚠️ THE CURVE ARRIVES AS (root_P50, root_c), NOT root_b. b is derived and
// set_traits re-derives it, so a row in it reaches nothing -- closed_form_rows.hpp carries
// the chain.
template <class T>
struct SupplyAt {
  const std::vector<T>& psi_soil;
  const std::vector<T>& r_R_H_min;
  const std::vector<T>& r_R_V_sum;
  const T& root_P50;
  const T& root_c;
};

// How the draw responds when the collar pulls harder: dE_up/dT_collar in total,
// and each layer's own term.
//
// ⚠️ THE TOTAL IS AT THE CALLER'S SCALAR AND THE LAYERS ARE AT DOUBLE, which is
// the split SupplyDraw keeps: the layer terms cross as SUPPLIED rows, and
// recording them would take the whole supply a second time for the same numbers.
// Units follow the draw's -- kg for the total, mol per layer.
template <class T>
struct CollarConductance {
  T total{};
  std::vector<double> per_layer;
};

class MultiLayerRoots {
public:
  // How often each clamp held. Shared storage, so a consumer copying this object
  // still reports -- plant copies its whole strategy per cohort and discards it,
  // and a plain member would take every count down with the copy.
  clamp_counter clamps;

  // --- root vulnerability trait pair (hazard 1: NOT the stem's b/c) ---------
  double root_c = 2.680147;       // unitless
  // THE TRAIT is root_P50, the potential at 50% loss of conductivity. root_b and
  // root_psi_crit are DERIVED from it and root_c by Leaf, and are kept here
  // because the supply path's own curve and bracket are written in terms of them.
  double root_P50 = 3.4;          // MPa, positive magnitude
  double root_b = 3.8982451221145307;   // MPa, = P50 / (ln 2)^(1/c)
  double root_psi_crit = 5.8702827267723245; // MPa, = P95 of that curve

  // NOTE: beta_R_H and beta_R_V used to live here. They are parameters of the
  // root-architecture model, not of water transport, and since #33 this class
  // takes resistances rather than the carbon they were applied to -- so they are
  // arguments to root_network_from_carbon and members of whoever owns that
  // model. In plant that is TF24_Strategy.

  // pre-computed root vulnerability curve f_r(m) = exp(-(m/root_b)^root_c)
  odelia::interpolator::hermite_interpolator<double> root_vuln_from_psi;
  // cumulative integral of it, G(m) = int_0^m f_r(s) ds, indexed by magnitude
  // m = -psi. Lets uptake() obtain the mean conductivity over a potential
  // interval from 2 evals instead of (n+1).
  odelia::interpolator::hermite_interpolator<double> root_vuln_integral_from_psi;
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
  // NOTE: `dz_` used to live here, a scalar layer thickness re-derived from the
  // profile on every set_soil_state. Nothing in this package had read it since
  // #33, and #626 removed the last reason to carry it: layer thickness is stated
  // by the caller now, per layer, and is an argument to root_network_from_carbon
  // rather than something this object derives a second opinion about.

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
    // df_r/dpsi in closed form, so the interpolant is handed both halves rather
    // than choosing slopes from the values -- see Leaf::setup_transpiration for
    // what a shape-preserving rule costs a reader of the slope. f_r'(0) = 0 for
    // root_c > 1, which the expression gives exactly at psi = 0.
    std::vector<double> y_f_r_slope(x_psi_root.size());
    for (size_t i = 0; i < x_psi_root.size(); ++i) {
      const double u = x_psi_root[i] / root_b;
      y_f_r_slope[i] = -y_f_r[i] * (root_c / root_b) * std::pow(u, root_c - 1.0);
    }
    root_vuln_from_psi.init(x_psi_root, y_f_r, y_f_r_slope);

    // ⚠️ NEITHER SPLINE REFUSES AN OUT-OF-DOMAIN READ. A hermite interpolator
    // extends linearly from the end knot in eval, value_at and slope alike, and
    // there is no switch to turn that off -- DO NOT reach for one, and do not
    // read a bound into the init call above. What bounds each curve is its
    // accessor below: root_vuln_at clamps the argument, root_vuln_integral_at
    // caps the value, and those are the only way the rest of this class reads a
    // curve.
    //
    // Integral: the extension is wanted, under a ceiling. Its limit is finite and
    // non-zero, and the end-knot polynomial tracks the true G closely over the
    // half-MPa it takes to reach that limit (3.46212 against 3.46176 at 7 MPa), so
    // capping the VALUE is both smooth and tighter than clamping the argument
    // would be. root_vuln_integral_at applies the cap.
    // dG_root/dpsi IS f_r by the fundamental theorem, and y_f_r is already built
    // above -- so this pair costs nothing beyond passing it.
    root_vuln_integral_from_psi.init(x_psi_root, y_integral, y_f_r);

    // Read off the spline rather than recomputed. cumulative_vulnerability_integral
    // sets its final knot to vulnerability_psi_max EXACTLY -- the `i == n` branch
    // is there so the two cannot disagree by a rounding -- so this is that number
    // (6.8918 MPa at the root defaults). Taking it from the curve keeps it one
    // number rather than two spellings of one.
    root_vuln_last_knot_ = root_vuln_from_psi.max();
    root_vuln_integral_limit_ =
        cumulative_vulnerability_integral_limit(root_b, root_c);
  }

  // f_r at a suction, clamped into the knot domain. This clamp is the whole of
  // the bound: past the last knot the spline extends linearly and runs negative,
  // and the last knot's ~1% is the driest conductivity this curve describes.
  //
  // Operand order is load-bearing: written this way both clamps return psi when
  // psi is NaN, where the reversed forms return the bound. The uptake call site's
  // !isfinite(f_ri) guard is what reads that NaN.
  double root_vuln_at(double psi) const {
    const double held = std::max(std::min(psi, root_vuln_last_knot_), 0.0);
    if (held != psi) clamps.note(CLAMP_ROOT_VULN_ARGUMENT);
    return root_vuln_from_psi.eval(held);
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
      return root_vuln_integral_limit_;
    }
    return raw;
  }

  // dG/dpsi, consistent with root_vuln_integral_at: zero wherever that returns
  // the cap, because there the value no longer depends on psi. Costs a second
  // spline eval, which is affordable here and would not be on the uptake hot path
  // -- duptake_dpsi runs once per gradient, not ~10^3 times per solve.
  double root_vuln_integral_deriv_at(double psi) const {
    if (root_vuln_integral_from_psi.eval(psi) >= root_vuln_integral_limit_) {
      return 0.0;
    }
    return root_vuln_integral_from_psi.slope(psi);
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
      //
      // Correct at ANY layer spacing, and left alone by #626 for that reason:
      // `soil_depth` holds the cumulative depth to the BOTTOM of each layer, so
      // averaging consecutive boundaries (with an implicit 0 at the surface) is
      // the midpoint of layer i whether or not the layers are equal.
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
  //
  // ⚠️ THE SAME WALK THE DERIVATIVE TAKES, at T = double and with no conductance
  // asked for. A second body for this path is a second answer to the same
  // question, free to drift from the one the rows are formed against -- which is
  // the whole reason the templated walk exists.
  void uptake(double T_collar, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake_at<double>(T_collar, supply_over(psi_soil_), soil_consumption, E_up);
  }

  // Uptake against an arbitrary vector of layer suctions, for callers that
  // want to probe the supply function away from the current soil state (the
  // R-facing Leaf::E_from_Soil_to_Root_Collar).
  //
  // The `&psi_soil == &psi_soil_` test is what used to select the
  // cached path for every caller, including the hot one. It is kept here only so
  // this entry point cannot change behaviour for a caller that happens to hand
  // back psi_soil_ itself; the hot path above no longer depends on
  // address identity to be fast.

  // The same, at any scalar. The cache is refused on the active path: it returns
  // a table read with no rows, which is a value that has quietly stopped
  // responding to anything -- the one failure this whole surface exists to avoid.
  //
  // ⚠️ THE COLLAR RESPONSE COMES OUT OF THE SAME WALK OR NOT AT ALL. It shares
  // this layer's span, integral and resistance, and a second walk forming them
  // again is a second spelling that has to agree bit for bit -- which is what
  // the deleted `duptake_dpsi_impl` said of itself, and what cost 73 tape
  // statements against this walk's 48 over five layers for one number. Pass
  // nullptr on the path that reads no slope; the forward solve calls this ~10^3
  // times per collar and is that path.
  template <class T>
  void uptake_at(const T& T_collar, const SupplyAt<T>& at,
                 std::vector<T>& soil_consumption, T& E_up,
                 CollarConductance<T>* conductance = nullptr) const {
    const bool cache = std::is_same_v<T, double> &&
                       (&at.psi_soil == reinterpret_cast<const std::vector<T>*>(&psi_soil_)) &&
                       root_vuln_integral_soil_.size() ==
                           static_cast<size_t>(max_soil_layer);
    uptake_impl<T>(T_collar, at, cache, soil_consumption, E_up, conductance);
  }

  void uptake_at(double T_collar, const std::vector<double>& psi_soil,
                 std::vector<double>& soil_consumption, double& E_up) const {
    uptake_at<double>(T_collar, supply_over(psi_soil), soil_consumption, E_up);
  }

  // d(E_up)/d(T_collar), for a caller that wants the conductance and not the
  // draw. It is a CONDUCTANCE and is positive by construction -- pulling harder
  // at the collar draws more water -- which is the whole reason for working in
  // magnitudes. The derivation is at the walk above; the per-layer terms stay in
  // MOL, matching the per-layer draws rather than the aggregate, and the kg
  // conversion is applied once, to the total.
  //
  // CONTRACT: returns NaN when any layer sits on a branch kink (T_collar ==
  // T_soil[i], the gravity-balance point, or T_collar == 0). That is deliberate,
  // not a failure -- the analytic general-branch derivative is not valid across
  // those, and the caller falls back to a central difference. An implementation
  // that threw, or returned 0, would silently degrade TF24f's acclimation
  // gradient. Any alternative supply path must keep this contract.
  //
  // It answers only where the draw does, because the walk that forms the slope
  // is the walk that forms the draw: a collar the uptake refuses as infeasible
  // refuses from here too, rather than handing back a slope for a point with no
  // flux.
  double duptake_dpsi(double T_collar,
                      const std::vector<double>& psi_soil) const {
    std::vector<double> per_layer;
    return duptake_dpsi(T_collar, psi_soil, per_layer);
  }

  double duptake_dpsi(double T_collar, const std::vector<double>& psi_soil,
                      std::vector<double>& per_layer) const {
    CollarConductance<double> conductance;
    std::vector<double> draw(psi_soil.size(), 0.0);
    double E_up = 0.0;
    uptake_impl<double>(T_collar, supply_over(psi_soil), false, draw, E_up,
                        &conductance);
    per_layer = std::move(conductance.per_layer);
    return conductance.total;
  }

private:
  // This object's members as a view over a caller's soil vector, so the double
  // entries reach the same body without copying their own state.
  SupplyAt<double> supply_over(const std::vector<double>& psi_soil) const {
    return {psi_soil, network_.r_R_H_min, network_.r_R_V_sum, root_P50, root_c};
  }

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

public:
  // This object's own state as a view, for the double path.
  SupplyAt<double> held_supply() const {
    return {psi_soil_, network_.r_R_H_min, network_.r_R_V_sum, root_P50, root_c};
  }

  // The uptake, at any scalar. Upstream's function line for line, with one
  // discipline added: ⚠️ EVERY COMPARISON READS PASSIVE. std::min on two active
  // values branches on a taped comparison and manufactures a discontinuity the
  // model does not have -- so the branch is chosen at double and the ACTIVE value
  // is selected, which is what std::min compiles to anyway.
  //
  // At double every supplied row below collapses to its table read and the cache path is
  // taken, so this instantiation is the original arithmetic.
  template <class T>
  void uptake_impl(const T& T_collar, const SupplyAt<T>& at,
                   bool use_integral_cache, std::vector<T>& soil_consumption,
                   T& E_up, CollarConductance<T>* conductance) const {
    using odelia::util::to_passive;
    const std::vector<T>& psi_soil = at.psi_soil;
    const double collar_at = to_passive(T_collar);
    // The band in which two potentials count as equal. It is the uptake's
    // branch width AND the slope's kink width because they are the same
    // condition: where the span vanishes the general branch divides by nothing,
    // and the derivative of the branch that replaces it is not the derivative
    // of this one.
    const double kink_tol = 1e-8;
    if (!std::isfinite(collar_at)) {
      util::stop_infeasible("uptake",
          "E_from_Soil_to_Root_Collar invalid input; T_collar=" +
          util::to_string(collar_at));
    }
    E_up = T(0.0);
    // ⚠️ ONE ENTRY PER SOIL LAYER, NOT PER ROOTED LAYER, which is the convention
    // the uptake already keeps: the loop below writes only as far as
    // max_soil_layer -- the deepest layer carrying roots -- and the layers under
    // it stay zero, because a layer no root reaches draws nothing and its collar
    // slope is nothing. Sized by the rooted count instead, a shallow-rooted plant
    // hands back a shorter vector than its own uptake, and the supplied row in
    // outputs_at pairs a row with the wrong layer or refuses on the length.
    // Measured: plant's gradient ladder refused a whole sweep on
    // "expected 5, received 2".
    //
    // ⚠️ max_soil_layer BELONGS TO THE ROOT NETWORK, NOT TO THE CALLER'S SOIL
    // VECTOR, and the loop below indexes all three by it. A caller probing the
    // supply at a shorter psi_soil than the network was built over -- which
    // duptake_dpsi and the R-facing E_from_Soil_to_Root_Collar both accept --
    // reads past the end of psi_soil and WRITES past the end of
    // soil_consumption. Checked here, once, against the bound the loop actually
    // uses, rather than at each entry point against a length each of them
    // derives differently.
    const std::size_t layers = static_cast<std::size_t>(max_soil_layer);
    if (psi_soil.size() < layers) {
      util::stop("uptake: the root network reaches " + util::to_string(max_soil_layer) +
                 " layers and psi_soil carries " + util::to_string(static_cast<int>(psi_soil.size())) +
                 "; a supply probe must cover every rooted layer");
    }
    // Grown rather than assigned: a caller handing a longer vector keeps its
    // length, which is the per-soil-layer convention above.
    if (soil_consumption.size() < layers) {
      soil_consumption.resize(layers, T(0.0));
    }
    if (conductance != nullptr) {
      conductance->total = T(0.0);
      conductance->per_layer.assign(psi_soil.size(), 0.0);
    }
    bool at_a_kink = false;

    auto curve = [&](const T& x) -> T {
      const double q = to_passive(x);
      if constexpr (std::is_same_v<T, double>) { return root_vuln_at(q); }
      else { return closed_form_curve<T>(root_vuln_at(q), x, at.root_P50, at.root_c); }
    };

    const double G_at_T_collar =
        use_integral_cache ? root_vuln_integral_at(collar_at) : 0.0;

    for (int i = 0; i < max_soil_layer; i++) {
      const T& psi_i = psi_soil[std::size_t(i)];
      const double soil_at = to_passive(psi_i);
      // The three collars the analytic slope is not valid across. Tested here
      // rather than in a second walk, because the layer that fails is the layer
      // whose span, integral and resistance are being formed anyway.
      if (conductance != nullptr &&
          (std::abs(collar_at - soil_at) < kink_tol ||
           std::abs((collar_at - soil_at) - grav_head_z_[i]) < kink_tol ||
           std::abs(collar_at) < kink_tol)) {
        at_a_kink = true;
      }
      const bool want_slope = conductance != nullptr && !at_a_kink;
      // ⚠️ BOUND, NOT COPIED. Both arms are lvalues that outlive the layer, so
      // this selects one; taken by value it would COPY an active, and a copy of
      // one is a recorded statement that the sweep then walks once per census
      // metric for nothing. Four of these a layer was 20 statements a placement
      // at five layers.
      const T& T_src_min = (collar_at < soil_at) ? T_collar : psi_i;
      const T& T_src_max = (soil_at < collar_at) ? T_collar : psi_i;

      if (std::abs(collar_at - soil_at) < kink_tol) {
        const T f_ri = curve(T_src_max);
        if (!std::isfinite(to_passive(f_ri)) || to_passive(f_ri) <= 0.0) {
          util::stop_infeasible("uptake",
              "E_from_Soil_to_Root_Collar invalid f_ri; layer=" +
              std::to_string(i) + "; f_ri=" +
              util::to_string(to_passive(f_ri)));
        }
        const T r_R = at.r_R_H_min[std::size_t(i)] / f_ri +
                      at.r_R_V_sum[std::size_t(i)];
        const T E_i = T(-grav_head_z_[i]) / r_R;
        soil_consumption[std::size_t(i)] = E_i;
        E_up += E_i;
      } else if (std::is_same_v<T, double> &&
                 std::abs((collar_at - soil_at) - grav_head_z_[i]) < kink_tol) {
        soil_consumption[std::size_t(i)] = T(0.0);
      } else {
        // Named so both arms are lvalues and the selection binds rather than
        // copies; T(0.0) records nothing either way.
        const T zero(0.0);
        const T& T_pos_lo = (to_passive(T_src_min) < 0.0) ? zero : T_src_min;
        const T& T_neg_hi = (0.0 < to_passive(T_src_max)) ? zero : T_src_max;

        auto G_integral = [&](const T& arg) -> T {
          const double q = to_passive(arg);
          if constexpr (std::is_same_v<T, double>) {
            if (use_integral_cache) {
              if (q == collar_at) return G_at_T_collar;
              if (q == soil_at) return root_vuln_integral_soil_[std::size_t(i)];
            }
            return root_vuln_integral_at(q);
          } else {
            return closed_form_integral<T>(root_vuln_integral_at(q),
                                     root_vuln_integral_deriv_at(q), arg,
                                     at.root_P50, at.root_c);
          }
        };

        T integral = T(0.0);
        if (to_passive(T_pos_lo) < to_passive(T_src_max)) {
          integral += G_integral(T_src_max) - G_integral(T_pos_lo);
        }
        if (to_passive(T_src_min) < to_passive(T_neg_hi)) {
          integral += (T_neg_hi - T_src_min);
        }

        const T span = T_src_max - T_src_min;
        const T r_R = at.r_R_H_min[std::size_t(i)] * span / integral +
                      at.r_R_V_sum[std::size_t(i)];
        const T num = T_collar - psi_i - T(grav_head_z_[i]);
        const T E_i = num / r_R;
        soil_consumption[std::size_t(i)] = E_i;
        E_up += E_i;

        if (want_slope) {
          // dE_i/dT_collar by the quotient rule, on the span, integral and
          // resistance THIS layer just formed. With
          //   dspan/dT  = sign_var  (+1 where the collar is the upper bound),
          //   dinteg/dT = sign_var * f_r(T_collar)  for T_collar > 0, else
          //               sign_var, since f_r is 1 there,
          //   dnum/dT   = 1,
          // and only r_R_H depends on the collar.
          const double sign_var = (collar_at > soil_at) ? 1.0 : -1.0;
          // ⚠️ THE MOVING BOUND'S INTEGRAND IS THE INTEGRAL'S OWN DERIVATIVE,
          // not the separate conductivity spline. The two agree on the knot
          // domain and are bounded differently past it -- the lookup clamps to
          // the last knot, the integral is capped at G(inf) -- so only this one
          // stays consistent with the `integral` above. closed_form_curve
          // supplies the ROWS; the VALUE stays the derivative table's, which is
          // what keeps that distinction.
          T fr_at = T(1.0);
          if (collar_at > 0.0) {
            const double table = root_vuln_integral_deriv_at(collar_at);
            if constexpr (std::is_same_v<T, double>) {
              fr_at = table;
            } else {
              fr_at = closed_form_curve<T>(table, T_collar, at.root_P50, at.root_c);
            }
          }
          const T dinteg_dT = T(sign_var) * fr_at;
          const T dr_R_dT = at.r_R_H_min[std::size_t(i)] *
                            (T(sign_var) * integral - span * dinteg_dT) /
                            (integral * integral);
          const T term = (r_R - num * dr_R_dT) / (r_R * r_R);
          conductance->per_layer[std::size_t(i)] = to_passive(term);
          conductance->total += term;
        }
      }
    }

    // ⚠️ THE TWO OUTPUTS CARRY DIFFERENT UNITS, and this line is the whole of the
    // difference. E_up leaves in kg H2O m^-2 s^-1, matching the rest of the leaf
    // and the environment; the per-layer draws stay in mol and are converted
    // downstream in TF24_Strategy::compute_rates.
    //
    // Dropping it is invisible layer by layer -- every soil_consumption[i] is
    // still exact -- and shows up only in the aggregate, as a clean factor of
    // 1/0.018015. Which is how it was caught here: a ratio that round is a
    // missing constant, not drift.
    E_up = E_up * T(kg_per_mol_h2o);
    if (conductance != nullptr) {
      if (at_a_kink) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        conductance->per_layer.assign(psi_soil.size(), nan);
        conductance->total = T(nan);
      } else {
        conductance->total = conductance->total * T(kg_per_mol_h2o);
      }
    }
  }
};

}  // namespace phylloptim

#endif
