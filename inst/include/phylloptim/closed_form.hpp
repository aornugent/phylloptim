// -*-c++-*-
#ifndef PHYLLOPTIM_CLOSED_FORM_HPP_
#define PHYLLOPTIM_CLOSED_FORM_HPP_

// The closed form: an alternative METHOD for reaching the leaf optimum, selected
// with `set_model(curve, "stem", "closed")` and then `optimise()`.
//
// WHY. The exact solve is the model's whole cost. It root-finds the first-order
// condition over a closed interval, and plant calls it millions of times.
//
// WHAT IT ACTUALLY BUYS, re-measured in `bench_solve` on the 96-point single-layer
// grid with both arms in ONE process (hazard 5). ⚠️ These supersede the 10.8x and
// 47x this header used to quote, which predated the solver merge -- back then the
// baseline was a golden-section search over the objective, where it now root-finds
// the first-order condition, and the old figures also priced the inversion alone
// without the output-writing evaluation a real caller needs:
//
//     exact TF24 stem solve        2.65 us/call   1x
//     closed TF24                  1.16 us/call   2.29x realised, phi = 0.271
//     exact CF77 stem solve        2.74 us/call   1x
//     closed CF77                  1.37 us/call   1.99x realised, phi = 0.385
//
// Inverting `1/[phi + (1-phi)/S]` puts the CEILING at S = 4.4x (TF24) and 5.2x
// (CF77) -- what the method would be worth at phi = 0. Quote the realised column.
//
// HOW. Every model in this family satisfies dA/dE = lambda at the optimum, and
// given lambda the solution collapses to the Medlyn USO form
//
//     ci/ca = xi/(xi + sqrt(D)),     xi = sqrt(Q/lambda)
//
// with Q below. So the whole problem is knowing lambda, which is a function of the
// psi the solve is looking for -- and where lambda is a POWER LAW in psi to leading
// order at the wet end, that circularity inverts explicitly. One Newton step on the
// unapproximated supply-minus-demand residual finishes the job.
//
// ============================================================================
// WHEN A CLOSED FORM EXISTS, WHICH IS A SHARPER CRITERION THAN THE BENEFIT LINK
// ----------------------------------------------------------------------------
// The objective is `h(A(psi)) - C(psi)` for every curve (see `Leaf::BenefitLink`),
// and the first-order condition is `h'(A)*dA/dpsi = dC/dpsi`. The closed form
// exists iff **h'(A) is independent of the solution** -- because that is what lets
// lambda be evaluated before A is known:
//
//   Identity  h' = 1                    TF24, CF77, JS22, CMax    -- qualifies
//   Scaled    h' = 1/|A|max             ProfitMax                 -- qualifies
//   Log       h' = 1/A                  SOX, JW26                 -- DOES NOT
//
// `Scaled` qualifies although h' is not 1, because `|A|max` is computed once per
// driver set by its own maximisation over the supply stream: at solve time it is a
// CONSTANT, and the condition becomes `dA/dpsi = |A|max * C'`, the same structure
// with a rescaled lambda.
//
// `Log` does not, and the failure is structural rather than algebraic: h' = 1/A, so
// lambda carries A, xi carries lambda, ci carries xi and A = A(ci) closes the loop.
// That is a genuine fixed point in A. It can be iterated but not inverted, and the
// two curves with that link are the two the refusal below names.
//
// Given the criterion, what remains is the wet-end exponent
// `n = dln(lambda)/dln(psi)`, because that is what has to be a power law for the
// inversion to be explicit:
//
//   curve      wet-end lambda                                    n
//   -------    ----------------------------------------------    -----------------
//   CF77       cf77_price(), constant                            0    (exact)
//   CMax       (a*psi+b)/(kmax*f) -> b/kmax, or ~psi if b == 0    0 if b != 0 else 1
//   JS22       2*gamma*dpsi/(kmax*f) -> 2*gamma*psi/kmax          1
//   TF24       ~psi^(stem_c*beta2 - 1)                            stem_c*beta2 - 1
//   ProfitMax  |A|max/(k(psi_soil)-kcrit) rescaling on |f'|/f     stem_c - 1
//   SOX, JW26  carries A                                          no power law
//
// IMPLEMENTED HERE: **TF24 and CF77 only.** The other four are refused, and the two
// refusals say different things on purpose -- see `optimise_into` below.
// ============================================================================
//
// ⚠️ HOW ACCURATE IT CAN POSSIBLY BE, WHICH IS NOT SET BY THE POWER LAW.
// ----------------------------------------------------------------------------
// The Newton step converges: 2, 3 and 5 steps give bit-identical answers over a
// 384-row driver grid, and they are not the exact solve's answers. So the residual's
// fixed point is not the optimum, and what separates them is the USO COLLAPSE
// itself, not `n`. Measured directly -- taking the EXACT lambda at the EXACT optimum
// and asking what `ca*xi/(xi + sqrt(D))` makes of it -- the relative error in ci is
//
//     0.8 % to 2 %      at 25 C with ci/ca > 0.6
//     4 % to 10 %       at 40 C, across the range
//     up to 139 %       at 25 C, psi_soil 4, ci/ca 0.16
//
// The reason is in `uso_group`: the 3*Gstar is the electron-transport-limited
// dA/dci, and the collapse additionally drops `(ca - ci)` beside
// `(ci-Gstar)(ci+2Gstar)/(3Gstar)` and assumes ci >> Gstar. Under Rubisco-bound
// colimitation, at high temperature, or at low ci, none of those is free. **So the
// closed form is a few-percent method at 25 C in wet soil and no better than its
// collapse anywhere** -- which is what `within_guard` is really guarding, and why
// the guard is calibrated on ci/ca rather than on anything about `n`.
//
// A polish that would remove this is available and NOT taken: given psi, the true
// first-order condition `lambda = A'(ca-ci)^2 / [kappa*D*(A'(ca-ci) + A)]` is a
// scalar equation in ci with `dassim_dci` already analytic here, so two Newton
// steps on it would make the fixed point the exact optimum. That is a second solver
// for the operating point, which is the one thing the developer guide's "One
// solver. Do not add a second." section exists to refuse -- and it would not fix
// note 4, which is the actual blocker.
//
// FIVE THINGS NOT TO GET WRONG, all learned the expensive way:
//
//  1. `newton_steps = 1` is deliberate. Two steps are *worse* in the tail.
//  2. The validity guard tests an OUTPUT (ci/ca > 0.5), so it can only be applied
//     after the fact, with a fallback to the exact solve. You cannot branch on it
//     up front. `optimise_into` does exactly that, and COUNTS the fallbacks.
//  3. Quote the realised speedup, not the ceiling. With a fallback fraction phi it
//     is 1/[phi + (1-phi)/speedup]. `Leaf::closed_form_fallback_fraction()` is what
//     measures phi; before it existed, phi had never been measured on a real
//     water-limited scenario and the ceiling was the only number anyone had.
//  4. ⚠️ **SMOOTHNESS IS THE BLOCKER, AND IT IS THE FALLBACK RATHER THAN THE
//     ARITHMETIC.** The argmax must vary smoothly with inputs because the
//     demographic growth-rate gradient depends on it (hazard 3), so this was
//     re-measured: a 201-point sweep of `kmax` over +/-30% -- which is how plant's
//     height reaches the leaf -- reporting the mean |second difference| of
//     `opt_psi_stem_`. Three regimes, and only the middle one is new information:
//
//         sweep never falls back      exact 6.1e-06   closed 7.2e-06    1.2x
//         sweep always falls back     identical by construction         1.0x
//         sweep CROSSES the guard     exact 1.0e-05   closed 1.0e-03    100x
//                                     exact 2.2e-05   closed 3.4e-03    153x
//
//     The closed form's own arithmetic is barely rougher than the solve -- 1.2x,
//     which is the register the reference reported (0.0015 against 0.0011). What
//     costs two orders of magnitude is the guard boundary: it tests an OUTPUT, so
//     the set of drivers on which the method switches is a hypersurface in driver
//     space, and the argmax JUMPS as a smooth input carries the leaf across it. A
//     mixture of two methods is discontinuous even where both are individually
//     smooth. That is structural, not a tolerance to tune, and it means **this
//     method must not sit under a differentiated path** -- which is why
//     `gradient::route_seat` seats the exact solve and says so.
//  5. ⚠️ Do not read `psi* = E/kmax`. The supply here is the vulnerability curve's
//     INTEGRAL, so the inversion is `psi* = G^-1(E/kmax + G(psi_upstream))` --
//     `psi_from_supply` below. `E/kmax` is the special case of a leaf with no
//     vulnerability curve at all, and it is not this leaf. The claim that it was
//     ("exact, and reduces to psi* = E/kmax") is withdrawn.
//
// LIMITATIONS, all refused rather than approximated (`Leaf::optimise_closed`):
//
//  - **Single-layer, stem route only.** The multi-layer closed form does not exist;
//    the collar route is a three-level nest with a soil-to-collar path this
//    inversion says nothing about.
//  - **The energy balance is refused, and that is issue #116's resolution.** Every
//    `D` below is `vpd_leaf_`, the leaf-to-air deficit Fick's law actually divides
//    by since #93 -- which is a no-op on this path, because `set_leaf_vpd` returns
//    `atm_vpd_` EXACTLY when `use_energy_balance_` is false. It is not a no-op with
//    the gate on, and there the form is not closed at all: `vpd_leaf_` depends on
//    Tleaf, Tleaf depends on E, and E is what these expressions solve for, so `D`
//    is not known at the moment of inversion. Option 3 of the three in #116, taken
//    because it is the honest one and because the measurement is in that issue: the
//    error with the gate bypassed is reported in the PR.

#include <phylloptim/constants.hpp>
#include <phylloptim/leaf_model.hpp>

#include <algorithm>
#include <cmath>

namespace phylloptim {
namespace closed_form {

// The group that converts a marginal cost of water into the USO slope xi:
//     xi = sqrt(Q/lambda),   Q = 3*Gstar*kg_to_mol_h2o/(1.67e-3)
// with Gstar the CO2 compensation point in Pa and lambda in umol CO2 (kg H2O)^-1.
//
// The 1.67e-3 is the H2O:CO2 diffusion ratio carrying the same 1e-3 scale factor
// that appears in the E-from-A relation below; it is kept in this one place rather
// than spread through the expressions. The ratio is `l.H2O_CO2_stom_diff_ratio_`, a
// SETTABLE field since #50 -- 1.67 by default where the g1 literature uses 1.6, so a
// caller comparing against fitted g1 values may want to set it.
//
// ⚠️ THE 3*Gstar IS WHERE THE EXACTNESS CLAIM LIVES, AND IT IS NARROWER THAN IT
// LOOKS. It comes from differentiating the ELECTRON-TRANSPORT-limited assimilation
// `A = J/4 * (ci - Gstar)/(ci + 2*Gstar)`, whose dA/dci is `J/4 * 3*Gstar/(ci +
// 2*Gstar)^2` -- that 3 is that curve's, not a universal constant. Under
// colimitation with the Rubisco term binding, the USO form is an approximation
// rather than the optimum, which is why even the n = 0 case (CF77) is checked
// against the exact solve rather than asserted equal to it.
inline double uso_group(const Leaf &l) {
  const double gstar_Pa = l.gamma_ * l.umol_per_mol_to_Pa_;
  return 3.0 * gstar_Pa * kg_to_mol_h2o / (l.H2O_CO2_stom_diff_ratio_ * 1e-3);
}

// Transpiration implied by assimilation on the DEMAND side, kg H2O m-2 s-1.
// The supply-side counterpart is Leaf::transpiration; the closed form works by
// driving their difference to zero.
//
// `vpd_leaf_` and not `atm_vpd_` -- see the energy-balance note in the header. The
// two are the same number on every path that reaches here.
inline double transpiration_from_assim(const Leaf &l, double assim, double ci) {
  return l.H2O_CO2_stom_diff_ratio_ * 1e-3 * assim * l.vpd_leaf_ /
         ((l.ca_ - ci) * kg_to_mol_h2o);
}

// Stomatal conductance to CO2 implied by a transpiration, mol CO2 m-2 s-1. The
// same expression `Leaf::stom_cond_CO2` uses, so the closed form's own diagnostic
// and the leaf's reported value cannot drift.
inline double stom_cond_from_E(const Leaf &l, double E) {
  return l.atm_kpa_ * E * kg_to_mol_h2o / l.vpd_leaf_ /
         l.H2O_CO2_stom_diff_ratio_;
}

// d(lambda)/d(psi), analytic. Differentiates Leaf::lambda_TF24, whose form is
//     lambda = K * (1-f)^(TF24_beta2-1) * p^(stem_c-1),  p = psi/stem_b,  f = exp(-p^stem_c)
//
// ⚠️ THE ONE PIECE JS22, CMax AND ProfitMax ARE MISSING. Each needs a `dlambda_*`
// of its own before `solve()` can serve it -- see the taxonomy in the header for
// the n each of them carries.
inline double dlambda_TF24(const Leaf &l, double psi) {
  const double K =
      l.TF24_cost_scale * l.TF24_beta2 * l.stem_c / (l.stem_b * l.leaf_specific_conductance_max_);
  const double p = psi / l.stem_b;
  const double f = std::exp(-std::pow(p, l.stem_c));
  return K * std::pow(p, l.stem_c - 2.0) *
         ((l.TF24_beta2 - 1.0) * std::pow(1.0 - f, l.TF24_beta2 - 2.0) * f * l.stem_c *
              std::pow(p, l.stem_c) +
          std::pow(1.0 - f, l.TF24_beta2 - 1.0) * (l.stem_c - 1.0)) /
         l.stem_b;
}

// dA/dci for the colimitation quadratic, analytic. Replaces the central
// difference an earlier prototype used, saving two assimilation evaluations per
// Newton step.
inline double dassim_dci(const Leaf &l, double ci, double electron_transport) {
  const double gstar = l.gamma_ * l.umol_per_mol_to_Pa_;
  const double ar = l.vcmax_ * (ci - gstar) / (ci + l.km_);
  const double ae =
      electron_transport / 4.0 * (ci - gstar) / (ci + 2.0 * gstar);
  const double dar =
      l.vcmax_ * (l.km_ + gstar) / ((ci + l.km_) * (ci + l.km_));
  const double dae = electron_transport / 4.0 * 3.0 * gstar /
                     ((ci + 2.0 * gstar) * (ci + 2.0 * gstar));
  const double s = ar + ae, ds = dar + dae;
  const double cv = l.curv_fact_colim;
  const double disc = std::sqrt(s * s - 4.0 * cv * ar * ae);
  const double ddisc =
      (2.0 * s * ds - 4.0 * cv * (dar * ae + ar * dae)) / (2.0 * disc);
  return (ds - ddisc) / (2.0 * cv);
}

struct Solution {
  double psi_stem;      // MPa, positive magnitude
  double ci;            // Pa
  double assim;         // umol CO2 m-2 s-1
  double transpiration; // kg H2O m-2 s-1
  double stom_cond_CO2; // mol CO2 m-2 s-1
  double g1_eff;        // the USO slope xi, kPa^0.5
};

// Invert the hydraulic supply for the stem potential that carries a transpiration:
//
//     E = kmax * [G(psi) - G(psi_upstream)]   =>   psi = G^-1(E/kmax + G(psi_up))
//
// ⚠️ NOT `E/kmax`. G is the cumulative vulnerability integral, so the two agree only
// where f == 1 -- i.e. nowhere on a real curve. See note 5 in the header.
//
// The domain check is done against the transpiration AT psi_crit rather than left to
// the spline, because `stem_curve_integral_inverse` throws past its domain and the
// only thing a caller here can do with that is report psi_crit. A demand above what
// the stem can carry at psi_crit IS the dry bound, so it is returned rather than
// raised: the guard downstream then sees the resulting ci/ca and falls back.
inline double psi_from_supply(Leaf &l, double E, double psi_upstream) {
  if (!(E > 0.0)) {
    return psi_upstream;
  }
  if (!(E < l.transpiration(l.psi_crit, psi_upstream))) {
    return l.psi_crit;
  }
  return l.transpiration_to_psi_stem(E, psi_upstream);
}

// Assemble the outputs implied by a stem potential, for TF24.
inline Solution evaluate_at(Leaf &l, double psi, double Q, double sqrt_D) {
  const double lambda = l.lambda_TF24(psi);
  const double xi = std::sqrt(Q / lambda);
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  return Solution{psi, ci, assim, E, stom_cond_from_E(l, E), xi};
}

// The wet-end starting potential, in units of stem_b. Pure arithmetic: no spline
// read, no assimilation evaluation, so this is free relative to the Newton step
// that follows it.
//
// Equating the power-law demand `C * p^(-n/2)` against a supply linearised at the
// wet end, `p - p_up`, gives
//
//     h(p) = (p - p_up) * p^(n/2) - C = 0
//
// ⚠️ THE `p_up` TERM IS NOT COSMETIC, AND THE VERSION OF THIS FILE THAT PREDATED
// BEING WIRED IN DID NOT HAVE IT. It solved the psi_upstream == 0 problem --
// `l.transpiration(psi, 0.0)`, matching the reference analysis, where the collar was
// held at zero. On this package's stem route the upstream potential is psi_soil,
// which the golden grid runs from 0.5 to 6.0 MPa, and ignoring it is not a small
// error: the supply available at a given psi is overstated by the whole wet-end
// integral. At p_up == 0 the Newton iteration below is exact in one step and
// reproduces the old explicit `C^(2/(n+2))`.
inline double wet_end_p(double C, double p_up, double n, double p_crit) {
  const double p0 = std::pow(C, 2.0 / (n + 2.0));
  if (!(p_up > 0.0)) {
    return p0;
  }
  // Newton on h, from above the root: h is increasing on p > p_up for n >= 0, so a
  // start at `p_up + p0` is on the far side of it and the iteration walks down.
  // Four steps, fixed: this is an algebraic scalar equation, not the model, and its
  // cost is four `pow` calls against the spline reads and root-find that follow.
  double p = p_up + p0;
  for (int i = 0; i < 4; ++i) {
    const double q = std::pow(p, n / 2.0);
    const double h = (p - p_up) * q - C;
    const double dh = q * (1.0 + 0.5 * n * (p - p_up) / p);
    if (!(dh > 0.0) || !std::isfinite(h)) {
      break;
    }
    const double next = p - h / dh;
    if (!std::isfinite(next)) {
      break;
    }
    // Strictly inside (p_up, p_crit): outside it the model has no operating point
    // and the `pow` of a negative base is a NaN that would propagate silently.
    p = std::min(std::max(next, p_up * (1.0 + 1e-12) + 1e-12), p_crit);
  }
  return p;
}

// General TF24_beta2. Explicit power-law leading order, then `newton_steps` Newton
// steps on the full supply-minus-demand residual. Requires set_physiology to have
// run. See note 1 in the header: leave newton_steps at 1.
inline Solution solve(Leaf &l, double psi_upstream, int newton_steps = 1) {
  const double kmax = l.leaf_specific_conductance_max_;
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double Q = uso_group(l);
  const double K_lambda = l.TF24_cost_scale * l.TF24_beta2 * l.stem_c / (l.stem_b * kmax);
  const double Xi = std::sqrt(Q / K_lambda);
  const double n = l.stem_c * l.TF24_beta2 - 1.0;
  const double electron_transport = l.electron_transport();

  // Leading order, taking kappa from the wet-end limit A(ci -> ca).
  const double assim_wet = l.assim_colimited(l.ca_ * (1.0 - 1e-9));
  const double kappa = l.H2O_CO2_stom_diff_ratio_ * 1e-3 * assim_wet /
                       (l.ca_ * kg_to_mol_h2o);
  const double C = kappa * sqrt_D * Xi / (kmax * l.stem_b);
  const double p_crit = l.psi_crit / l.stem_b;
  double p = wet_end_p(C, psi_upstream / l.stem_b, n, p_crit);

  for (int k = 0; k < newton_steps; ++k) {
    const double psi = p * l.stem_b;
    // The bracket, not just positivity: below psi_upstream no water moves and the
    // residual's supply term is negative, which Newton would chase outward.
    if (!(psi > psi_upstream) || psi >= l.psi_crit) {
      break;
    }
    const double lambda = l.lambda_TF24(psi);
    const double xi = std::sqrt(Q / lambda);
    const double ci = l.ca_ * xi / (xi + sqrt_D);
    const double assim = l.assim_colimited(ci);
    const double dassim = dassim_dci(l, ci, electron_transport);
    const double u = l.ca_ - ci;
    const double E = transpiration_from_assim(l, assim, ci);
    const double dE_dci = l.H2O_CO2_stom_diff_ratio_ * 1e-3 * l.vpd_leaf_ /
                          kg_to_mol_h2o * (dassim * u + assim) / (u * u);
    const double dci_dxi = l.ca_ * sqrt_D / ((xi + sqrt_D) * (xi + sqrt_D));
    const double dxi_dpsi = -0.5 * xi / lambda * dlambda_TF24(l, psi);
    // Residual: hydraulic supply minus stomatal demand, both kg H2O m-2 s-1. The
    // supply is measured from psi_upstream, so its DERIVATIVE is unchanged --
    // kmax*f(psi) either way -- while its VALUE is not.
    const double R = l.transpiration(psi, psi_upstream) - E;
    const double dR =
        kmax * std::exp(-std::pow(p, l.stem_c)) - dE_dci * dci_dxi * dxi_dpsi;
    if (dR == 0.0 || !std::isfinite(R) || !std::isfinite(dR)) {
      break;
    }
    double psi_next = psi - R / dR;
    psi_next = std::max(psi_upstream,
                        std::min(psi_next, l.psi_crit * (1.0 - 1e-9)));
    p = psi_next / l.stem_b;
  }

  return evaluate_at(l, p * l.stem_b, Q, sqrt_D);
}

// TF24_beta2 == 1/stem_c. The psi dependence cancels out of lambda entirely (n = 0),
// so xi is constant and there is nothing to solve for it -- no power law, no Newton
// step. The stem potential still has to come out of the supply inversion, which is
// one spline solve.
//
// ⚠️ A VALID SPECIAL CASE AT AN UNMOTIVATED PARAMETER VALUE, AND THIS FILE USED TO
// CLAIM MORE. It said 1/stem_c was 0.917, "within a hair" of the beta2 <= 1 argued
// for on independent grounds -- so the algebraic convenience looked biologically
// preferred. That coincidence was an artefact of a vulnerability-curve error: the
// stem curve had been built with *E. saligna*'s measured P12 (1.85 MPa) sitting in
// the P50 slot, giving stem_c = 1.09 where the correct pair (P50 = 3.40, P88 = 5.16)
// gives 2.680147 -- which is this package's default. So **1/stem_c = 0.373**, and
// reaching beta2 ~ 1 through beta2 = 1/stem_c would require stem_c ~ 1, i.e. a
// non-threshold "exponential" vulnerability curve that the eucalypt literature
// falsifies (De Kauwe et al. 2022 Table 1: stem_c 2.86-8.30 over 15 species, none
// below 2.8). The case is kept because it is a real algebraic simplification and a
// clean test of the n = 0 arm; the claim that it is the preferred parameterisation
// is withdrawn.
inline Solution solve_exact_beta2(Leaf &l, double psi_upstream) {
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double Q = uso_group(l);
  const double xi = std::sqrt(Q * l.stem_b * l.leaf_specific_conductance_max_ /
                              (l.TF24_cost_scale * l.TF24_beta2 * l.stem_c));
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  // ⚠️ THIS USED TO RETURN NaN HERE, which was tolerable while nothing consumed the
  // Solution and is not now: `opt_psi_stem_` is the operating point. The potential
  // is recovered from the supply, which is the only place it can come from.
  return Solution{psi_from_supply(l, E, psi_upstream),
                  ci, assim, E, stom_cond_from_E(l, E), xi};
}

inline bool beta2_is_exact(const Leaf &l, double tol = 1e-12) {
  return std::abs(l.TF24_beta2 - 1.0 / l.stem_c) <= tol * std::max(1.0, 1.0 / l.stem_c);
}

// CF77, the n = 0 case: lambda is `cf77_price()`, a constant, so there is no
// circularity to break and no iteration at all. xi, ci, A and E follow in one pass
// and the potential comes out of the supply inversion.
//
// ⚠️ "NO ITERATION" IS NOT THE SAME AS "EXACT", and conflating them is what note 5
// in the header withdraws. Two separate approximations remain: the USO form's 3*Gstar
// is the electron-transport-limited dA/dci (see `uso_group`), so under Rubisco-bound
// colimitation this ci is not the optimal one; and `psi_from_supply` is a spline
// inversion, exact to the spline. The n = 0 exactness is a statement about lambda,
// which is the only thing the power-law machinery approximates.
inline Solution solve_CF77(Leaf &l, double psi_upstream) {
  const double sqrt_D = std::sqrt(l.vpd_leaf_);
  const double xi = std::sqrt(uso_group(l) / l.cf77_price());
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  return Solution{psi_from_supply(l, E, psi_upstream),
                  ci, assim, E, stom_cond_from_E(l, E), xi};
}

// The validity guard. The closed form degrades where the leaf is far from the
// wet-end limit its leading order is expanded about, and ci/ca is the diagnostic:
// the reference reports good agreement while ci/ca > 0.5 and does not claim it
// below. Tests an OUTPUT, so it can only be applied after solving -- see note 2.
//
// ⚠️ BELOW ci/ca ~ 0.3 THE FIXED POINT CEASES TO EXIST, which is why this is a guard
// rather than a tolerance. The exact optimum turns back toward shutdown while the
// closed form runs on to psi_crit; that is the wrong topology, not a loose answer,
// and no amount of Newton refinement recovers it.
// ⚠️ AND IT ALSO REQUIRES THE ANSWER TO BE STRICTLY INSIDE THE BRACKET, which the
// ci/ca test does not imply and which the version of this file that predated being
// wired in did not check. Measured over a 384-row driver grid: on the dim, hot rows
// -- psi_soil 1, PPFD 100, VPD 2, 40 C, where net assimilation at the exact optimum
// is 0.08 umol -- the inversion returns a potential at or below psi_soil while
// reporting ci/ca = 0.975, and clamping that to psi_soil lands on the NO-FLOW branch
// (A = -R_d = -3.17 against the exact 0.08). That was the whole of the grid's worst
// error. A bound is where the two methods are least comparable and where the exact
// solve is cheapest to be right about -- it tests the gradient's sign there
// explicitly -- so a closed form that arrives at one hands the row back.
inline bool within_guard(const Leaf &l, const Solution &s) {
  return std::isfinite(s.psi_stem) && std::isfinite(s.ci) &&
         std::isfinite(s.assim) && s.ci / l.ca_ > 0.5 &&
         s.psi_stem > l.supply_psi_soil_scalar() && s.psi_stem < l.psi_crit;
}

// ---------------------------------------------------------------------------
// The entry point: solve, guard, and write the operating point -- or fall back.
// ---------------------------------------------------------------------------
//
// ⚠️ EVERY OUTPUT THE EXACT SOLVE WRITES IS WRITTEN HERE, AND BY THE SAME
// FUNCTIONS. Hazard 8: `Leaf` is a value member that plant reuses for every
// individual in a patch, so an output this path declined to write would silently
// become the PREVIOUS solve's value -- a plausible number about a different leaf.
// The way that is guaranteed rather than audited is that the closed form is used
// ONLY to locate psi*, and everything else comes from `profit_psi_stem_for<K>` and
// `lambda_for<K>`, which is exactly what `optimise_psi_stem_single` calls at the end
// of its own search. So the two paths share the entire output-assembly code and
// cannot disagree about which fields exist.
//
// The consequence worth knowing when reading output: `ci_`, `assim_colimited_` and
// `transpiration_` are the EXACT values at the approximate psi*, not the closed
// form's own USO estimates of them. That is the right semantics -- the method
// approximates an argmax, and the state at a given potential is not in question --
// and it is what makes the comparison against the exact solve a comparison of
// argmaxes.
template <Leaf::CostCurve K>
inline void optimise_into(Leaf &l) {
  using CostCurve = Leaf::CostCurve;
  if constexpr (K == CostCurve::TF24 || K == CostCurve::CF77) {
    l.clear_collar_solve_state();
    l.check_cost_parameters<K>();
    ++l.closed_form_calls_;

    const double psi_soil = l.supply_psi_soil_scalar();
    l.opt_psi_stem_ = psi_soil;

    // The degenerate branch, and it is the exact solve's branch verbatim: soil
    // drier than the stem can reach leaves exactly one feasible potential, so
    // there is nothing for either method to choose between. Kept identical rather
    // than shared because the shared version would be a fourth caller of a
    // two-line body.
    if (psi_soil > l.psi_crit) {
      l.profit_ = l.profit_psi_stem_for<K>(psi_soil, psi_soil);
      l.lambda_emergent_ = l.lambda_for<K>(psi_soil, psi_soil);
      return;
    }

    Solution s{};
    if constexpr (K == CostCurve::CF77) {
      s = solve_CF77(l, psi_soil);
    } else {
      s = beta2_is_exact(l) ? solve_exact_beta2(l, psi_soil)
                            : solve(l, psi_soil);
    }

    if (!within_guard(l, s)) {
      // ⚠️ THE FALLBACK IS THE WHOLE EXACT SOLVE, not a refinement of s. The guard
      // fires where the fixed point may not exist, so there is nothing to refine.
      ++l.closed_form_fallbacks_;
      l.last_solve_fell_back_ = true;
      l.optimise_psi_stem_single<K>();
      return;
    }

    // Into the bracket the exact solve maximises over. The Newton step clamps
    // already; this is what makes the clamping a property of the answer rather
    // than of one code path.
    const double psi =
        std::max(psi_soil, std::min(s.psi_stem, l.psi_crit));
    l.opt_psi_stem_ = psi;
    l.profit_ = l.profit_psi_stem_for<K>(psi, psi_soil);
    l.lambda_emergent_ = l.lambda_for<K>(psi, psi_soil);
  } else if constexpr (K == CostCurve::JS22 || K == CostCurve::CMax ||
                       K == CostCurve::ProfitMax) {
    // ⚠️ FEASIBLE, NOT IMPLEMENTED -- and the distinction from the arm below is the
    // point. All three have an h'(A) that is independent of the solution (JS22 and
    // CMax by the identity link; ProfitMax because |A|max is a per-driver-set
    // constant), so the inversion exists. What is missing is one `dlambda_*` each
    // for the Newton step, plus, for ProfitMax, picking up the |A|max rescaling of
    // lambda -- its `profitmax_A_max_` and `profitmax_k_span_` are already computed
    // and exposed, so the constants are in hand.
    util::stop("the closed form is not implemented for the " +
               Leaf::curve_name(static_cast<int>(K)) + " cost curve. It EXISTS "
               "for this curve -- its benefit-link derivative does not depend on "
               "the solution, and its wet-end lambda is a power law (n = 1 for "
               "JS22; 0 for CMax unless CMax_b is zero, then 1; stem_c - 1 for "
               "ProfitMax) -- but the analytic dlambda/dpsi the Newton step needs "
               "has not been written, and ProfitMax additionally needs the |A|max "
               "rescaling. Use method \"exact\", or add the derivative to "
               "closed_form.hpp.");
  } else {
    static_assert(K == CostCurve::SOX || K == CostCurve::JW26,
                  "unhandled CostCurve in closed_form::optimise_into");
    // ⚠️ NOT A MISSING PIECE. This is the one place the closed form cannot exist,
    // and saying so precisely matters more than the refusal: the LOG benefit link
    // has h'(A) = 1/A, so lambda carries A, xi carries lambda, ci carries xi and
    // A = A(ci) closes the loop. It is a fixed point in the unknown, which can be
    // iterated but not inverted. Every other link's h' is a constant at solve time.
    util::stop("the closed form does not exist for the " +
               Leaf::curve_name(static_cast<int>(K)) + " cost curve, and it is "
               "not a matter of unwritten algebra. SOX and JW26 are the two "
               "curves with the LOG benefit link, where h'(A) = 1/A -- so lambda "
               "carries the assimilation the solve is for, and given lambda -> xi "
               "-> ci -> A the dependency closes into a fixed point rather than "
               "an inversion. Use method \"exact\".");
  }
}

} // namespace closed_form

// The closed-form arm of `Leaf::optimise()`. Declared in `leaf_model.hpp`, defined
// here because it needs `closed_form` above and `closed_form` needs the complete
// `Leaf` -- see the declaration for why that is legal and what breaks if the
// include order changes.
//
// ⚠️ THE THREE REFUSALS ARE CONFIGURATION CHECKS, SO THEY BELONG AT SOLVE TIME AND
// NOT IN `set_model`. Every one of route, supply topology and energy balance is
// settable after the model is seated, so a check at seating time would either fire
// on a leaf that was about to become valid or miss a leaf that had just stopped
// being.
inline void Leaf::optimise_closed() {
  if (route_is_collar_) {
    util::stop("the closed form is a STEM-route method: it inverts psi_stem "
               "against a supply measured from psi_soil, and says nothing about "
               "the soil-to-collar path. The leaf is seated on the collar route; "
               "call set_model(\"" +
               curve_name(static_cast<int>(cost_curve_)) +
               "\", \"stem\", \"closed\"), or use method \"exact\".");
  }
  if (!supply_is_single_layer()) {
    util::stop("the closed form is single-layer only: the multi-layer case is a "
               "three-level nest and its closed form does not exist. psi_soil has "
               + util::to_string(static_cast<int>(supply_n_layers())) +
               " layers; use one layer, or method \"exact\".");
  }
  if (use_energy_balance_) {
    util::stop("the closed form is refused with use_energy_balance_ on, and this "
               "is a statement about the form rather than about its accuracy "
               "(#116). Fick's law divides by vpd_leaf_, which depends on Tleaf, "
               "which depends on E -- and E is what the inversion solves for. So "
               "the deficit D is not known at the moment of inversion and the "
               "form is not closed. Use method \"exact\", or prescribe the leaf "
               "temperature by setting use_energy_balance_ to false.");
  }
  with_curve(cost_curve_, [&](auto tag) {
    closed_form::optimise_into<tag.value>(*this);
  });
}

} // namespace phylloptim

#endif
