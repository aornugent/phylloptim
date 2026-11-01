// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_STEM_ROUTES_HPP_
#define PHYLLOPTIM_LEAF_STEM_ROUTES_HPP_

// The stem routes: the operating point maximised over the stem potential
// with the upstream potential held at psi_soil, over the closed interval
// (hazard 11). A stem route is a different model from the collar solve.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>
#include <phylloptim/leaf_cost_curves.hpp>

namespace phylloptim {

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

} // namespace phylloptim

#endif
