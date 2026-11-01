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

// The per-layer root hydraulic resistances the supply solve consumes.
struct RootNetwork {
  // Fully hydrated horizontal (soil -> root) resistance, divided by the layer's
  // mean conductivity fraction over its potential span to give the actual one.
  std::vector<double> r_R_H_min;
  // Cumulative vertical, inter-layer resistance from the surface down to layer i.
  std::vector<double> r_R_V_sum;
  // Diagnostics only -- read by nothing in the model.
  std::vector<double> c_r_V, c_r_H, r_R_V;
};

// The root-architecture model: root carbon per layer -> resistances. A helper
// for the caller; nothing in this package calls it, since set_physiology takes
// the resistances. Each layer's carbon is split 1/3 vertical : 2/3 horizontal:
//
//   r_R_H_min[i] = beta_R_H / c_r_h
//   r_R_V[i]     = beta_R_V * dz[i]^2 / c_r_v
//   r_R_V_sum[i] = r_R_V[0] + ... + r_R_V[i]
//
// `dz` is per layer, which is what makes the total vertical resistance
// independent of how the column is sliced (test_root_network_from_carbon); a
// single dz breaks that for unequal layers. Outputs are sized to the deepest
// rooted layer. A zero-carbon layer above it gets r_R_H_min = 0, infinite
// conductance where there are no roots: backwards, unreachable from plant, and
// kept rather than changed silently. Fills `out` in place so a per-solve
// caller reuses its buffers: five fresh vectors per call measured slower.
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
    // Every layer, not only the rooted ones: the layers below the rooting depth
    // are the ones a growing plant will reach.
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

    out.r_R_H_min[i] = beta_R_H / c_r_h;
    // dz^2 of this layer: vertical conductivity scales with root cross-section,
    // and the segment resisted is the one spanning this layer. Squared into its
    // own variable, because `beta_R_V * dz * dz` associates as (beta_R_V*dz)*dz
    // and moves the golden file by a last bit.
    const double dz_sq = dz[i] * dz[i];
    out.r_R_V[i] = beta_R_V * dz_sq / c_r_v;
    vertical_resistance_sum += out.r_R_V[i];
    out.r_R_V_sum[i] = vertical_resistance_sum;
  }
}

// The same map at any scalar, so a caller can put the carbon -> resistance chain
// on its own tape. Writes only the two series the supply reads. Its arithmetic is
// the double form's term for term, so at double it is bit-identical
// (test_root_network_template_is_the_double_overload).
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

  // The deepest layer carrying carbon, read passive: which layers are rooted is
  // piecewise constant in the carbon, and a taped comparison would make a jump.
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
    // Squared first, as in the double form.
    const double dz_sq = dz[std::size_t(i)] * dz[std::size_t(i)];
    vertical_resistance_sum += T(beta_R_V) * T(dz_sq) / c_r_v;
    r_R_V_sum[std::size_t(i)] = vertical_resistance_sum;
  }
}

// Returns a fresh network: for tests and one-off callers.
inline RootNetwork root_network_from_carbon(
    const std::vector<double>& root_carbon_per_layer,
    const std::vector<double>& dz,
    double beta_R_H, double beta_R_V) {
  RootNetwork out;
  root_network_from_carbon(root_carbon_per_layer, dz, beta_R_H, beta_R_V, out);
  return out;
}

// Per-layer thicknesses from a cumulative soil-depth profile, with an implicit 0
// at the surface: for a caller holding a profile and not the widths. A caller
// holding its own widths (plant's TF24_Environment) must pass those: differencing
// is not bit-exact (`diff` of `(i+1)*depth/n` misses `depth/n` by a last bit at
// the package defaults), and that perturbs every vertical resistance.
inline std::vector<double>
layer_thicknesses(const std::vector<double>& soil_depth) {
  if (soil_depth.empty()) {
    util::stop("layer_thicknesses: soil_depth must have at least one layer");
  }
  std::vector<double> out(soil_depth.size());
  double previous = 0.0;
  for (size_t i = 0; i < soil_depth.size(); ++i) {
    // Strictly increasing, so every width is positive. Widths passed as
    // boundaries pass this whenever they increase with depth, and only the
    // caller can tell (test_layer_thicknesses asserts that limit).
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
// Given a collar potential, how much water the root system delivers and how
// fast that changes. The gas-exchange core consumes only
//
//   E_up = uptake(P_collar)   and   dE_up/dP = duptake_dpsi(P_collar)
//
// The model (after Potkay et al. 2021) is parallel soil layers, each joined to
// one root collar. Layer i draws (T_collar - T_soil[i] - gravity_head *
// z_soil_mid[i]) / (r_R_H + r_R_V), where r_R_V is r_R_V_sum[i] and
// r_R_H = r_R_H_min[i] / mean f_r, the mean of the root conductivity fraction
// over the span between soil and collar, (G(hi) - G(lo)) / (hi - lo).
//
// Every water potential is a positive magnitude in MPa (hazard 2): the minus
// signs live in the equations. E_up (kg H2O m^-2 s^-1) and soil_consumption
// (mol per layer) are passed in by reference, because they are plant's buffers:
// it writes back into them by name after crown integration.

// The active inputs the uptake path reads, as a view: the double path builds one
// over this object's members, the active path over its caller's vectors, and one
// uptake_impl serves both without copying.
// ⚠️ The curve arrives as (root_P50, root_c), not root_b, which is derived
// (hazard 1): a row in b reaches nothing. closed_form_rows.hpp carries the chain.
template <class T>
struct SupplyInputs {
  const std::vector<T>& psi_soil;
  const std::vector<T>& r_R_H_min;
  const std::vector<T>& r_R_V_sum;
  const T& root_P50;
  const T& root_c;
};

// dE_up/dT_collar, in total (kg, at the caller's scalar) and per layer (mol, at
// double). The layer terms cross as supplied rows, as SupplyDraw's do, so
// recording them would take the whole supply a second time.
template <class T>
struct CollarConductance {
  T total{};
  std::vector<double> per_layer;
  // False for a caller that reads only the total; per_layer is then left alone.
  bool keep_per_layer = true;
};

class MultiLayerRoots {
public:
  // How often each clamp held. Shared storage, so the counts survive plant
  // copying its strategy per cohort.
  clamp_counter clamps;

  // --- root vulnerability trait pair (hazard 1: NOT the stem's b/c) ---------
  double root_c = 2.680147;       // unitless
  // The trait is root_P50, the potential at 50% loss of conductivity; root_b and
  // root_psi_crit are derived from it and root_c by Leaf, and kept here because
  // the supply path's curve and bracket are written in them.
  double root_P50 = 3.4;          // MPa, positive magnitude
  double root_b = 3.8982451221145307;   // MPa, = P50 / (ln 2)^(1/c)
  double root_psi_crit = 5.8702827267723245; // MPa, = P95 of that curve

  // The root's G as a table at root_c, shared with the stem's where the two c
  // agree. f_r has no table: it is exp(-(psi/root_b)^root_c).
  std::shared_ptr<const VulnerabilityTable> vulnerability_;
  InverseScale inv_b_;

  // --- soil geometry -------------------------------------------------------
  // Initialised to the sentinels clear() assigns, so a bare MultiLayerRoots is
  // never indeterminate.
  double soil_number_of_depths_ = util::na_value_int;
  int max_soil_layer = util::na_value_int;  // deepest layer with non-zero root mass
  std::vector<double> soil_depth_;
  std::vector<double> z_soil_mid_;
  // gravity_head * z_soil_mid_[i], set by set_soil_state for uptake's hot loop.
  std::vector<double> grav_head_z_;
  bool use_precomputed_z_soil_mid_ = false;
  // Scratch for duptake_dpsi(), reused across calls: mutable and unsynchronised,
  // so one solve at a time per object.
  mutable std::vector<double> draw_scratch_;
  mutable CollarConductance<double> conductance_scratch_;

  // --- soil state ----------------------------------------------------------
  std::vector<double> psi_soil_;           // positive magnitudes, as supplied
  // root_vuln_integral_at(psi_soil_[i]), constant across a collar solve and built
  // by begin_solve, so a cached uptake() reads the table once per call rather
  // than twice per layer.
  std::vector<double> root_vuln_integral_soil_;

  // --- root resistance network --------------------------------------------
  // One object, so root_network_from_carbon can refill it in place.
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

  // Seat the root's table at root_c, or share `reuse` (the stem's) when it is
  // already that table. A layer drier than the curve is an ordinary state: both
  // reads below answer at every potential (past the end, G(inf) and exp(-x)).
  void setup_vulnerability(
      double resolution,
      const std::shared_ptr<const VulnerabilityTable>& reuse = nullptr) {
    vulnerability_ = vulnerability_table(root_c, resolution, reuse);
  }

  // f_r at a suction: the surviving root conductivity, exp(-(psi/root_b)^c). It
  // underflows to exactly zero far past the curve (above ~46 MPa at the
  // defaults), and the one reader that divides by it, uptake_impl's
  // equal-potentials branch, takes a layer with none left as drawing nothing.
  double root_vuln_at(double psi) const {
    if (psi <= 0.0) {
      return 1.0;  // the curve starts at full conductivity; NaN falls through
    }
    return exp(-pow(psi / root_b, root_c));
  }

  // G(inf) = root_b Gamma(1/c) / c: the whole area under the conductivity curve.
  double root_vuln_integral_limit() const {
    return root_b * vulnerability_->G_inf_over_b();
  }

  // G at a suction. It tends to G(inf) as the layer dries, so the flux through a
  // layer is bounded by G(inf) / r_R_H_min however wide its span.
  double root_vuln_integral_at(double psi) const {
    return vulnerability_integral(*vulnerability_, psi, root_b, inv_b_.of(root_b));
  }

  // dG/dpsi = f_r, in closed form.
  double root_vuln_integral_deriv_at(double psi) const {
    return vulnerability_integral_deriv(*vulnerability_, psi, inv_b_.of(root_b));
  }

  // Q = G(inf) - G(psi), the area left under the conductivity curve past psi,
  // to its own relative precision however small (the table's tail). Zero past
  // the end of the curve, where a layer conducts nothing.
  double root_vuln_complement_at(double psi) const {
    return vulnerability_complement(*vulnerability_, psi, root_b, inv_b_.of(root_b));
  }

  // Whether psi is in the table's tail, where a layer's integral between two
  // potentials is formed as a difference of Q rather than of G (uptake_impl):
  // as a difference of G it keeps none of the digits G(inf) swamps.
  bool root_in_tail(double psi) const {
    return vulnerability_in_tail(*vulnerability_, psi, inv_b_.of(root_b));
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
      // Where no midpoints were precomputed. `soil_depth` is the cumulative depth
      // to the bottom of each layer, so the mean of consecutive boundaries is
      // the midpoint at any spacing.
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

  // Per-timestep root resistance network: the resistances, not the root carbon
  // they come from, which is the caller's root-architecture model. c_r_V, c_r_H
  // and r_R_V ride along as diagnostics plant exposes through RcppR6.
  // Takes const& and copy-assigns: the caller refills one network per call, and a
  // move would empty its buffers and make it reallocate all five next time.
  void set_root_network(const RootNetwork& network) {
    if (network.r_R_V_sum.size() != network.r_R_H_min.size()) {
      util::stop("set_root_network: r_R_H_min and r_R_V_sum must have the same "
                 "length; got " + std::to_string(network.r_R_H_min.size()) +
                 " and " + std::to_string(network.r_R_V_sum.size()));
    }
    // uptake() indexes psi_soil_ and grav_head_z_ by rooted layer, so a network
    // deeper than the soil profile would be an out-of-bounds read.
    if (soil_number_of_depths_ > 0 &&
        network.r_R_H_min.size() >
            static_cast<size_t>(soil_number_of_depths_)) {
      util::stop("set_root_network: network has " +
                 std::to_string(network.r_R_H_min.size()) +
                 " rooted layers but the soil profile has only " +
                 std::to_string(static_cast<int>(soil_number_of_depths_)));
    }
    for (size_t i = 0; i < network.r_R_H_min.size(); ++i) {
      // Non-negative. Zero means infinite conductance, which
      // root_network_from_carbon gives a zero-carbon layer (see its comment).
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

  // Per-solve entry point: builds root_vuln_integral_soil_ (valid until the next
  // call) and returns the wettest rooted layer's suction, the smallest magnitude,
  // which is the collar solve's bracket end.
  double begin_solve() {
    root_vuln_integral_soil_.resize(max_soil_layer);
    double wettest_soil_layer = std::numeric_limits<double>::infinity();
    for (int i = 0; i < max_soil_layer; ++i) {
      root_vuln_integral_soil_[i] = root_vuln_integral_at(psi_soil_[i]);
      wettest_soil_layer = std::min(wettest_soil_layer, psi_soil_[i]);
    }
    return wettest_soil_layer;
  }

  // Uptake at a collar suction, against the soil state begin_solve() cached: the
  // hot path, called at every collar evaluation. It is the walk the derivative
  // takes, at double, so the draw and the rows formed against it cannot drift.
  void uptake(double T_collar, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake_at<double>(T_collar, supply_over(psi_soil_), soil_consumption, E_up);
  }

  // Uptake against any layer suctions, at any scalar. The integral cache is
  // taken only at double and only when `at` views psi_soil_ itself: on the
  // active path a cached read would carry no rows. The collar response comes out
  // of this same walk, never a second one: pass `conductance` to read it, and
  // nullptr on the forward solve.
  template <class T>
  void uptake_at(const T& T_collar, const SupplyInputs<T>& at,
                 std::vector<T>& soil_consumption, T& E_up,
                 CollarConductance<T>* conductance = nullptr) const {
    const bool cache = std::is_same_v<T, double> &&
                       (&at.psi_soil == reinterpret_cast<const std::vector<T>*>(&psi_soil_)) &&
                       root_vuln_integral_soil_.size() ==
                           static_cast<size_t>(max_soil_layer);
    uptake_impl<T>(T_collar, at, cache, soil_consumption, E_up, conductance);
  }

  // For a caller probing the supply away from the current soil state (the
  // R-facing Leaf::E_from_Soil_to_Root_Collar).
  void uptake_at(double T_collar, const std::vector<double>& psi_soil,
                 std::vector<double>& soil_consumption, double& E_up) const {
    uptake_at<double>(T_collar, supply_over(psi_soil), soil_consumption, E_up);
  }

  // d(E_up)/d(T_collar): a conductance, positive by construction. The total is in
  // kg, the per-layer terms in mol, as the draws are.
  //
  // ⚠️ NaN when any layer sits on a branch kink (T_collar == T_soil[i], the
  // gravity balance, or T_collar == 0), and the caller then takes a central
  // difference (hazard 6); an alternative supply path must keep that rather than
  // throw or return 0. Refuses where the uptake refuses. On the hot path, so its
  // vectors are scratch this object keeps.
  double duptake_dpsi(double T_collar,
                      const std::vector<double>& psi_soil) const {
    double E_up = 0.0;
    conductance_scratch_.keep_per_layer = false;
    uptake_impl<double>(T_collar, supply_over(psi_soil), false, draw_scratch_,
                        E_up, &conductance_scratch_);
    return conductance_scratch_.total;
  }

  double duptake_dpsi(double T_collar, const std::vector<double>& psi_soil,
                      std::vector<double>& per_layer) const {
    double E_up = 0.0;
    conductance_scratch_.keep_per_layer = true;
    uptake_impl<double>(T_collar, supply_over(psi_soil), false, draw_scratch_,
                        E_up, &conductance_scratch_);
    per_layer = conductance_scratch_.per_layer;
    return conductance_scratch_.total;
  }

private:
  // This object's members as a view over a caller's soil vector, so the double
  // entries reach the same body without copying their own state.
  SupplyInputs<double> supply_over(const std::vector<double>& psi_soil) const {
    return {psi_soil, network_.r_R_H_min, network_.r_R_V_sum, root_P50, root_c};
  }

public:
  // This object's own state as a view, for the double path.
  SupplyInputs<double> held_supply() const {
    return {psi_soil_, network_.r_R_H_min, network_.r_R_V_sum, root_P50, root_c};
  }

  // The most the soil can deliver, the uptake's limit as the collar dries without
  // bound: the sum over layers of Q(psi_i) / r_R_H_min[i], Q being the area left
  // under the conductivity curve, with rows from Q's closed form so they survive
  // as Q goes to zero. A layer past the end of the curve adds nothing; a zero
  // r_R_H_min makes the limit infinite. kg H2O m^-2 leaf s^-1, as E_up.
  template <class T>
  T uptake_limit(const SupplyInputs<T>& at) const {
    using odelia::util::to_passive;
    T total = T(0.0);
    for (int i = 0; i < max_soil_layer; ++i) {
      const T& psi_i = at.psi_soil[i];
      const double q = to_passive(psi_i);
      T I_inf;
      if constexpr (std::is_same_v<T, double>) {
        I_inf = root_vuln_complement_at(q);
      } else {
        I_inf = complement_with_closed_form_rows<T>(root_vuln_complement_at(q),
                                                    -root_vuln_integral_deriv_at(q),
                                                    psi_i, at.root_P50, at.root_c);
      }
      if (!(to_passive(I_inf) > 0.0)) {
        continue;
      }
      if (!(to_passive(at.r_R_H_min[i]) > 0.0)) {
        return T(std::numeric_limits<double>::infinity());
      }
      total += I_inf / at.r_R_H_min[i];
    }
    // In kg, as E_up leaves the walk; the per-layer draws are in mol.
    return total * T(kg_per_mol_h2o);
  }

  // The wet bound's floor: the lowest gravity balance psi_i + z_i among layers
  // that still conduct, where every layer draws nothing or less. Where none
  // conducts, the wettest layer's balance, which keeps a shut leaf continuous
  // past the end of the curve. Returns that layer; `conducts` says which case.
  int wet_floor_layer(const std::vector<double>& psi_soil, bool& conducts) const {
    int lowest = -1, wettest = -1;
    for (int i = 0; i < max_soil_layer; ++i) {
      const double s = psi_soil[static_cast<std::size_t>(i)];
      if (wettest < 0 || s < psi_soil[static_cast<std::size_t>(wettest)]) {
        wettest = i;
      }
      if (!(root_vuln_complement_at(s) > 0.0)) {
        continue;
      }
      const double b = s + grav_head_z_[static_cast<std::size_t>(i)];
      if (lowest < 0 ||
          b < psi_soil[static_cast<std::size_t>(lowest)] +
                  grav_head_z_[static_cast<std::size_t>(lowest)]) {
        lowest = i;
      }
    }
    conducts = lowest >= 0;
    return conducts ? lowest : wettest;
  }

  // The highest gravity balance among layers that still conduct, where uptake is
  // at least zero; at g_min = 0 the wet bound lies between this and the floor.
  // The floor's value where none conducts.
  double wet_ceiling(const std::vector<double>& psi_soil) const {
    double highest = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < max_soil_layer; ++i) {
      const double s = psi_soil[static_cast<std::size_t>(i)];
      if (!(root_vuln_complement_at(s) > 0.0)) {
        continue;
      }
      highest = std::max(highest, s + grav_head_z_[static_cast<std::size_t>(i)]);
    }
    if (std::isinf(highest)) {
      bool conducts = false;
      const std::size_t k = static_cast<std::size_t>(wet_floor_layer(psi_soil, conducts));
      return psi_soil[k] + grav_head_z_[k];
    }
    return highest;
  }

  // The uptake walk, at any scalar. Writes E_up (kg H2O m^-2 leaf s^-1) and
  // soil_consumption[i] (mol H2O m^-2 leaf s^-1); a negative entry is a layer
  // gaining water (hydraulic redistribution). At double every supplied row
  // collapses to its table read. Every comparison reads passive: a branch on a
  // taped comparison makes a discontinuity the model does not have, so the
  // branch is chosen at double and the active value selected.
  template <class T>
  void uptake_impl(const T& T_collar, const SupplyInputs<T>& at,
                   bool use_integral_cache, std::vector<T>& soil_consumption,
                   T& E_up, CollarConductance<T>* conductance) const {
    using odelia::util::to_passive;
    const std::vector<T>& psi_soil = at.psi_soil;
    const double collar_psi = to_passive(T_collar);
    // Where two potentials count as equal: the uptake's branch width and the
    // slope's kink width, one condition, since where the span vanishes the
    // general branch divides by nothing and its replacement has another slope.
    const double kink_tol = 1e-8;
    if (!std::isfinite(collar_psi)) {
      util::stop_infeasible("uptake",
          "E_from_Soil_to_Root_Collar invalid input; T_collar=" +
          util::to_string(collar_psi));
    }
    E_up = T(0.0);
    // One entry per soil layer, not per rooted layer: layers below max_soil_layer
    // stay zero. Sized by the rooted count, a shallow-rooted plant's supplied rows
    // in billed_outputs_at pair with the wrong layer or refuse on the length.
    // max_soil_layer belongs to the root network, not to the caller's psi_soil, so
    // a shorter probe would be read past its end: checked here, once, against the
    // bound the loop uses.
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
      if (conductance->keep_per_layer) {
        conductance->per_layer.assign(psi_soil.size(), 0.0);
      }
    }
    bool at_a_kink = false;

    auto curve = [&](const T& x) -> T {
      const double q = to_passive(x);
      if constexpr (std::is_same_v<T, double>) { return root_vuln_at(q); }
      else { return curve_with_closed_form_rows<T>(root_vuln_at(q), x, at.root_P50, at.root_c); }
    };

    const double G_at_T_collar =
        use_integral_cache ? root_vuln_integral_at(collar_psi) : 0.0;

    // Raw pointers taken once: the walk writes soil_consumption, which the
    // compiler cannot prove aliases none of the vectors it reads, so references
    // would reload every data pointer on every layer.
    const T* const psi_p = psi_soil.data();
    const T* const rH_p = at.r_R_H_min.data();
    const T* const rV_p = at.r_R_V_sum.data();
    const double* const zg_p = grav_head_z_.data();
    const double* const Gs_p = root_vuln_integral_soil_.data();
    T* const out_p = soil_consumption.data();
    // At double, plain copies, for the same aliasing reason; at an active scalar
    // a copy is a recorded statement, so references.
    using hold_t = std::conditional_t<std::is_same_v<T, double>, const T, const T&>;
    hold_t Tc = T_collar;
    for (int i = 0; i < max_soil_layer; i++) {
      hold_t psi_i = psi_p[i];
      const double soil_at = to_passive(psi_i);
      // The three collars the analytic slope is not valid across (duptake_dpsi's
      // NaN contract).
      if (conductance != nullptr &&
          (std::abs(collar_psi - soil_at) < kink_tol ||
           std::abs((collar_psi - soil_at) - zg_p[i]) < kink_tol ||
           std::abs(collar_psi) < kink_tol)) {
        at_a_kink = true;
      }
      const bool want_slope = conductance != nullptr && !at_a_kink;
      // Bound, not copied: both arms outlive the layer, and a copy of an active is
      // a recorded statement the sweep walks once per census metric.
      hold_t T_src_min = (collar_psi < soil_at) ? Tc : psi_i;
      hold_t T_src_max = (soil_at < collar_psi) ? Tc : psi_i;

      if (std::abs(collar_psi - soil_at) < kink_tol) {
        const T f_ri = curve(T_src_max);
        if (!std::isfinite(to_passive(f_ri)) || to_passive(f_ri) < 0.0) {
          util::stop_infeasible("uptake",
              "E_from_Soil_to_Root_Collar invalid f_ri; layer=" +
              std::to_string(i) + "; f_ri=" +
              util::to_string(to_passive(f_ri)));
        }
        // -z / (r_H / f + r_V), multiplied through by f so a fully embolised
        // layer (f_r underflowed to zero) draws exactly nothing rather than
        // dividing by zero. With no horizontal resistance it is r_V alone.
        T E_i;
        if (to_passive(rH_p[i]) > 0.0) {
          E_i = T(-zg_p[i]) * f_ri / (rH_p[i] + rV_p[i] * f_ri);
        } else {
          E_i = T(-zg_p[i]) / rV_p[i];
        }
        out_p[i] = E_i;
        E_up += E_i;
      } else if (std::is_same_v<T, double> &&
                 std::abs((collar_psi - soil_at) - zg_p[i]) < kink_tol) {
        out_p[i] = T(0.0);
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
              if (q == collar_psi) return G_at_T_collar;
              if (q == soil_at) return Gs_p[i];
            }
            return root_vuln_integral_at(q);
          } else {
            return integral_with_closed_form_rows<T>(root_vuln_integral_at(q),
                                     root_vuln_integral_deriv_at(q), arg,
                                     at.root_P50, at.root_c);
          }
        };

        // From a lower end in the table's tail the span's integral is a
        // difference of what is left, Q(lo) - Q(hi), which keeps its digits
        // where G(hi) - G(lo) keeps none.
        auto Q_integral = [&](const T& arg) -> T {
          const double q = to_passive(arg);
          if constexpr (std::is_same_v<T, double>) {
            return root_vuln_complement_at(q);
          } else {
            return complement_with_closed_form_rows<T>(root_vuln_complement_at(q),
                                     -root_vuln_integral_deriv_at(q), arg,
                                     at.root_P50, at.root_c);
          }
        };

        T integral = T(0.0);
        if (to_passive(T_pos_lo) < to_passive(T_src_max)) {
          if (root_in_tail(to_passive(T_pos_lo))) {
            integral += Q_integral(T_pos_lo) - Q_integral(T_src_max);
          } else {
            integral += G_integral(T_src_max) - G_integral(T_pos_lo);
          }
        }
        if (to_passive(T_src_min) < to_passive(T_neg_hi)) {
          integral += (T_neg_hi - T_src_min);
        }

        // A layer whose whole span lies past the end of the curve conducts
        // nothing, and is skipped with no flux and no rows: on the tape the
        // arithmetic would be inf * 0. A shut leaf's collar can sit there.
        if (!(to_passive(integral) > 0.0)) {
          out_p[i] = T(0.0);
          continue;
        }
        const T span = T_src_max - T_src_min;
        const T r_R = rH_p[i] * span / integral +
                      rV_p[i];
        const T num = Tc - psi_i - T(zg_p[i]);
        const T E_i = num / r_R;
        out_p[i] = E_i;
        E_up += E_i;

        if (want_slope) {
          // dE_i/dT_collar by the quotient rule on this layer's span, integral
          // and resistance; only r_R_H depends on the collar. With
          //   dspan/dT  = sign_var  (+1 where the collar is the upper bound),
          //   dinteg/dT = sign_var * f_r(T_collar), or sign_var for T_collar <= 0,
          //   dnum/dT   = 1.
          const double sign_var = (collar_psi > soil_at) ? 1.0 : -1.0;
          // The moving bound's integrand is the integral's own derivative,
          // f_r(T_collar), in closed form; curve_with_closed_form_rows supplies its rows.
          T fr_at = T(1.0);
          if (collar_psi > 0.0) {
            const double table = root_vuln_integral_deriv_at(collar_psi);
            if constexpr (std::is_same_v<T, double>) {
              fr_at = table;
            } else {
              fr_at = curve_with_closed_form_rows<T>(table, T_collar, at.root_P50, at.root_c);
            }
          }
          const T dinteg_dT = T(sign_var) * fr_at;
          const T dr_R_dT = rH_p[i] *
                            (T(sign_var) * integral - span * dinteg_dT) /
                            (integral * integral);
          const T term = (r_R - num * dr_R_dT) / (r_R * r_R);
          if (conductance->keep_per_layer) {
            conductance->per_layer[std::size_t(i)] = to_passive(term);
          }
          conductance->total += term;
        }
      }
    }

    // The two outputs' units differ by this line alone: E_up leaves in
    // kg H2O m^-2 s^-1, the per-layer draws stay in mol (plant converts them in
    // TF24_Strategy::compute_rates). Dropped, every layer still looks exact and
    // only the total is off, by a factor of exactly 1/0.018015.
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
