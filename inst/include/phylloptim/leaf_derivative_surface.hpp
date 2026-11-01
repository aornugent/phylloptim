// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_DERIVATIVE_SURFACE_HPP_
#define PHYLLOPTIM_LEAF_DERIVATIVE_SURFACE_HPP_

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
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

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

} // namespace phylloptim

#endif
