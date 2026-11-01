# The hand-written R surface above the generated `Leaf` class in R/RcppR6.R, and
# the only route from here into C++. It splits traits from numerical settings and
# assigns the traits the 15-argument constructor does not take.

# Eucalyptus saligna; the defaults for every function below. They must equal
# phylloptim::Leaf's default constructor: test-golden.R solves through them and
# compares bit-exactly with a golden file from a default-constructed C++ Leaf.
# The order is the gradient's trait vector: append, never reorder (hazard 13).
# JS22_gamma and CMax_a are solved by align_curves(); test-align-curves.R checks.
# root_c / root_P50 default in MultiLayerRoots, not Leaf (hazard 1); they are
# restated because the R constructor passes them.
.leaf_trait_defaults <- list(
  vcmax_25 = 96,
  stem_c = 2.680147,
  stem_P50 = 3.4,
  root_c = 2.680147,
  root_P50 = 3.4,
  TF24_beta2 = 1.5,
  jmax_25 = 157.44,
  phi_J = 0.30 / 0.85,
  curv_fact_elec_trans = 0.7,
  curv_fact_colim = 0.99,
  TF24_cost_scale = 7.5,
  R_d_25 = 1.44,
  JS22_gamma = 0.4717887793,
  CMax_a = 0.6385596982,
  CMax_b = 0.0,
  g_min = 0.0
)

##' Physiological traits for a leaf
##'
##' The traits, and only the traits. Numerical settings live in
##' [leaf_control()]; the C++ constructor mixes the two and this is the split.
##'
##' Defaults are *Eucalyptus saligna*, and are the same values
##' `phylloptim::Leaf`'s default constructor uses.
##'
##' @section Two vulnerability curves, not one:
##' `stem_P50`/`stem_c` describe the STEM Weibull curve, which drives the
##' hydraulic cost; `root_P50`/`root_c` describe the ROOT curve, which drives
##' uptake. They are separate parameters with separate meanings and they are
##' named accordingly, because they were once an unmarked `b`/`c` pair alongside
##' `root_b`/`root_c` and an analysis used the root parameters for the stem cost.
##'
##' Each curve has exactly TWO settable numbers. `stem_b`, `psi_crit`, `root_b`
##' and `root_psi_crit` are quantiles of the curve those two define -- readable on
##' the resulting `Leaf` and not settable anywhere. See [weibull_p50_c()].
##'
##' @section psi_crit is not a free trait:
##' `psi_crit` is the stem curve's **P95** and is derived from `stem_P50`/`stem_c`,
##' not set: a free `psi_crit` would describe a curve it was not derived from. It
##' is a labelled point on a curve that continues past it, and the dry end of the
##' solve only for SOX, JW26, CF77 and LeastCost, whose own forms end there; for
##' the other curves the stem may run past it (#183).
##'
##' At the defaults:
##'
##' ```
##' stem_b = 3.898245, stem_c = 2.680147
##' 3.898245 * log(1/0.05)^(1/2.680147) = 5.870283 = psi_crit
##' ```
##'
##' to six decimal places. It moves with the pair, so there is nothing left to set
##' inconsistently. `vignette("fitting")` derives the
##' pair from a published P50/P88 pair, which is the right way round.
##'
##' @param vcmax_25 maximum carboxylation rate at 25 C (umol m^-2 s^-1)
##' @param stem_c shape parameter of the stem vulnerability curve (unitless)
##' @param stem_P50 stem potential at 50% loss of conductivity (MPa, positive
##'   magnitude). `stem_b` and `psi_crit` are derived from this pair and are
##'   readable but not settable on the resulting `Leaf`.
##' @param root_c shape parameter of the root vulnerability curve (unitless)
##' @param root_P50 root potential at 50% loss of conductivity (MPa, positive
##'   magnitude), with `root_b` and `root_psi_crit` derived from it.
##' @param TF24_beta2 exponent for the effect of hydraulic risk (unitless)
##' @param jmax_25 maximum electron transport rate at 25 C (umol m^-2 s^-1)
##' @param phi_J quantum yield of electron transport per ABSORBED photon
##'   (mol mol^-1). The leaf absorbs `alpha_PAR_` (0.85 by default) of the
##'   intercepted `PPFD`, so the default 0.30 / 0.85 is the old incident-basis
##'   yield of 0.30 restated. A gas-exchange fit identifies only
##'   `phi_J * alpha_PAR_`; fix `alpha_PAR_` to make the fitted value per absorbed
##'   photon. Was `a` (incident basis) before 0.13.0: `phi_J = a / 0.85`.
##' @param curv_fact_elec_trans curvature of the light response curve (unitless)
##' @param curv_fact_colim curvature of the colimited photosynthesis equation
##' @param TF24_cost_scale cost parameter for the TF24 profit model
##'   (umol m^-2 s^-1)
##' @param R_d_25 dark respiration at 25 C (umol m^-2 s^-1). It is the value at
##'   25 C only: respiration rises from there on Tjoelker's declining-Q10 curve.
##' @param JS22_gamma Joshi & Stocker (2022)'s hydraulic unit cost
##'   (umol C m^-2 s^-1 MPa^-2), read only by the `JS22` cost curve and returning
##'   an exactly zero gradient on every other one.
##'
##'   ⚠️ **The default is not a literature value.** It is aligned to TF24 by
##'   [align_curves()]: the value at which the JS22 stem optimum has TF24's
##'   conductance at that function's reference environment. No single value matches more
##'   than one soil potential: the `JS22_gamma` reproducing TF24's cost at its own
##'   optimum runs 0.287 at `psi_soil` 0.5 MPa to 1.761 at 3.0 MPa. TF24's cost
##'   tracks the absolute potential and this one tracks the drop, and across a
##'   drydown those move in opposite directions.
##' @param CMax_a,CMax_b slope and intercept of the `CMax` marginal cost,
##'   `dC/dpsi = CMax_a * psi + CMax_b`, following Wolf et al. (2016) as
##'   Anderegg et al. (2018) parameterised it. Read only by the `CMax` curve.
##'   `CMax_a`'s default is aligned to TF24 by [align_curves()], with `CMax_b` at
##'   zero; neither default is a literature value.
##'
##'   ⚠️ **`CMax_b` is signed, and negative in the source convention.** Sabot's
##'   `TractLSM` stores it with the sign inverted "so as to get a positive
##'   parameter value", so a value taken from that code or from a calibration
##'   against it must have its sign checked rather than assumed.
##'
##'   ⚠️ **Neither default is a literature value.** With `CMax_b = 0` the `CMax_a`
##'   reproducing TF24's cost at its own optimum runs 0.406 at `psi_soil` 0.5 MPa
##'   to 0.807 at 3.0 MPa — a 2.0x range against `JS22_gamma`'s 6.1x, which is the
##'   absolute-versus-drop distinction as a number.
##' @param g_min the residual (cuticular) conductance to water vapour at 25 C
##'   (mol H2O m^-2 s^-1). A leaf keeps losing water through its cuticle after
##'   its stomata have shut, and the cuticle lets CO2 through as well (at the
##'   ratio `$g_min_co2_ratio_`, 0.05 by default). The leaf optimises its
##'   STOMATAL flux; the xylem carries the total, `E = E_s + E_r` with
##'   `E_r = g_min(T_leaf) D_leaf / P_atm`, and every cost that reads the flux
##'   charges the total. Default 0, at which nothing moves.
##'
##'   ⚠️ **A measured cuticular conductance, not a fitted `g0`.** A `g0` fitted
##'   to gas exchange is an order of magnitude larger and is mostly stomata that
##'   have not finished closing; measured cuticular conductances run
##'   0.25-0.8 mmol m^-2 s^-1 (2.5e-4 to 8e-4 here). The Medlyn route's `$g0` is
##'   the fitted quantity and is untouched by this trait. `g_min` rises with
##'   leaf temperature on a smooth curve whose three parameters are fields on
##'   the `Leaf` (`$g_min_k_low_`, `$g_min_k_high_`, `$g_min_T_trans_`,
##'   `$g_min_T_width_`); their defaults are placeholders, inert at 25 C.
##'
##' @section Where the two root-resistance constants went:
##' `beta_R_H` and `beta_R_V` were traits here until #33. They parameterise the
##' root-ARCHITECTURE model, which the leaf no longer runs -- it takes the
##' resistances. They are arguments to [root_network_from_carbon()] now.
##'
##' The cost is real and worth stating: [leaf_gradient()] can no longer
##' differentiate with respect to either. Perturbing them is cheap, because
##' [root_network_from_carbon()] is homogeneous of degree 1 in each, but the
##' derivative of a solved output still needs two solves either side.
##'
##' @return A `leaf_traits` object; a named list.
##' @seealso [leaf_control()], [leaf_model()], [leaf_solve()]
##' @examples
##' leaf_traits()
##' # A more brittle stem. Only the curve's own two parameters are set: the scale
##' # `stem_b` and the critical potential `psi_crit` are quantiles of that curve
##' # and are derived from them, so they cannot be set inconsistently.
##' leaf_traits(vcmax_25 = 120, stem_P50 = 2.2)
##' @export
leaf_traits <- function(vcmax_25 = 96,
                        stem_c = 2.680147,
                        stem_P50 = 3.4,
                        root_c = 2.680147,
                        root_P50 = 3.4,
                        TF24_beta2 = 1.5,
                        jmax_25 = 157.44,
                        phi_J = 0.30 / 0.85,
                        curv_fact_elec_trans = 0.7,
                        curv_fact_colim = 0.99,
                        TF24_cost_scale = 7.5,
                        R_d_25 = 1.44,
                        JS22_gamma = 0.4717887793,
                        CMax_a = 0.6385596982,
                        CMax_b = 0.0,
                        g_min = 0.0) {
  out <- list(vcmax_25 = vcmax_25, stem_c = stem_c, stem_P50 = stem_P50,
              root_c = root_c, root_P50 = root_P50, TF24_beta2 = TF24_beta2,
              jmax_25 = jmax_25, phi_J = phi_J,
              curv_fact_elec_trans = curv_fact_elec_trans,
              curv_fact_colim = curv_fact_colim,
              TF24_cost_scale = TF24_cost_scale,
              R_d_25 = R_d_25, JS22_gamma = JS22_gamma,
              CMax_a = CMax_a, CMax_b = CMax_b, g_min = g_min)
  .check_scalars(out, "leaf_traits")
  if (R_d_25 < 0) {
    stop("leaf_traits(): R_d_25 must be non-negative", call. = FALSE)
  }
  if (g_min < 0) {
    stop("leaf_traits(): g_min must be non-negative", call. = FALSE)
  }
  structure(out, class = c("leaf_traits", "list"))
}

##' Numerical settings for a leaf solve
##'
##' Tolerances and iteration counts. These are not traits: they control how
##' precisely the solve is carried out, not what is being solved, and separating
##' them from [leaf_traits()] is the point of having two functions.
##'
##' @section What is deliberately not here:
##' A range for the leaf temperature. Nothing clamps it: with the energy balance
##' on it is affine in the transpiration, and [set_drivers()] refuses drivers
##' under which a leaf the solve can reach would be past what the model can
##' evaluate (R_d's Q10 reaching zero, about 119 C at the defaults), rather than
##' holding the temperature at a bound.
##'
##' @param GSS_tol_abs a width on the collar potential, despite its name, which
##'   dates from a golden-section search the solve no longer uses. ⚠️ **This does not set how well the operating point is
##'   determined.** Every solver in the package now root-finds its own first-order
##'   condition instead of searching the objective, so changing this leaves the
##'   answer bit-identical on the production path. What it still does, and it is
##'   less than the name suggests: it is read in exactly two places, both on the
##'   COLLAR route -- the width below which the feasible collar interval is treated
##'   as a single point, and the tolerance of the golden-section FALLBACK used when
##'   either bracket endpoint has no usable gradient (no driver sweep has reached
##'   that fallback).
##'
##'   ⚠️ **It does not reach the stem route at all.** That refines by a root-find
##'   on `dJ/dpsi == 0` at a hardcoded tolerance; the `(cell width) * 1e-4` figure
##'   is its own golden-section fallback, not its normal path. Tightening this to
##'   sharpen a stem-route argmax does nothing.
##' @param vulnerability_curve_ncontrol number of knot intervals in each Weibull
##'   vulnerability curve's table (a C2 quintic of its cumulative integral).
##'   Higher is more accurate and slower to construct; it does not affect solve
##'   speed. The default, 800, is read from the C++ `Leaf::ncontrol_default`, so
##'   the R and C++ defaults cannot differ.
##' @param ci_abs_tol absolute tolerance for root-solving intercellular CO2
##' @param ci_niter maximum iterations for root-solving intercellular CO2
##' @param integration_rule accepted and ignored. It selected a Gauss-Kronrod
##'   rule order in plant's compiled quadrature; the header-only adaptive
##'   Simpson that replaced it has no rule order. Retained because plant's R
##'   bindings expose it.
##' @param integration_tol absolute tolerance for the direct quadrature in
##'   `transpiration_full_integration()`. Not on the production path, which
##'   reads the vulnerability table instead.
##'
##' @return A `leaf_control` object; a named list.
##' @seealso [leaf_traits()], [leaf_model()]
##' @examples
##' leaf_control()
##' leaf_control(GSS_tol_abs = 1e-5)
##' @export
# `vulnerability_curve_ncontrol` is read from C++ `Leaf::ncontrol_default`: a
# literal would be a second copy whose drift the golden file cannot see (a table
# resolution only refines the answer). test-surface.R asserts it; plant's
# test-control.R compares plant's Control with `leaf_control()`.
leaf_control <- function(GSS_tol_abs = 1e-3,
                         vulnerability_curve_ncontrol = .leaf_ncontrol_default(),
                         ci_abs_tol = 1e-3,
                         ci_niter = 1000,
                         integration_rule = 21,
                         integration_tol = 1e-8) {
  out <- list(GSS_tol_abs = GSS_tol_abs,
              vulnerability_curve_ncontrol = vulnerability_curve_ncontrol,
              ci_abs_tol = ci_abs_tol, ci_niter = ci_niter,
              integration_rule = integration_rule,
              integration_tol = integration_tol)
  .check_scalars(out, "leaf_control")
  structure(out, class = c("leaf_control", "list"))
}

##' Choose how water reaches the root collar
##'
##' Two supply paths, and the choice belongs with the model rather than with the
##' drivers, so it is made once in [leaf_model()].
##'
##' `leaf_supply_multilayer()` is the default and is what plant uses: a soil
##' profile of one or more layers, each with its own water potential, feeding a
##' root network whose resistances are derived from a root carbon profile.
##'
##' `leaf_supply_singlelayer()` collapses the whole soil-to-collar path to **one
##' resistance**. This is what a leaf physiologist arriving from `plantecophys`
##' or `tealeaves` actually has: a soil water potential and no root-mass profile.
##' It is also what makes comparison against other optimality models meaningful,
##' since Medlyn, Prentice least-cost and Cowan-Farquhar are all formulated
##' against a single soil potential — comparing them against a multi-layer root
##' network compares two things at once.
##'
##' @section Both paths are configured the same way, and take the same drivers:
##' Neither function takes a resistance. The soil-to-collar resistance is a
##' per-call **driver** on both paths, supplied to [set_drivers()] as
##' `root_network` — [root_network_from_carbon()] on the multi-layer path,
##' [series_resistance()] on this one. It used to be an argument here, so the same
##' quantity arrived at construction on one path and per call on the other, and
##' `resistance` was the only differentiable parameter whose setter reset the whole
##' object.
##'
##' ⚠️ **`gravity_head` is the one asymmetry left, and it is deliberate.** The
##' multi-layer path derives a per-layer head from the depth profile it is handed
##' (`gravity_head * z_soil_mid`); this path has no depth profile to derive one
##' from, and a bare leaf wants **zero** rather than a geometric default — which is
##' precisely the caller `leaf_supply_singlelayer()` exists for. Making it a driver here
##' would mean either inventing a depth for a leaf that has none, or adding a
##' second supply-shaped argument to [set_drivers()] that only one path reads. If
##' you do want the multi-layer rule for a single layer of thickness `d`, pass
##' `gravity_head = 0.00981 * d / 2`.
##'
##' @section Why this is not a settable field:
##' The obvious R interface would be `leaf$supply_kind <- "single"`. It is not
##' offered, because flipping the tag on its own leaves the other path's state
##' configured and silently ignored, and flipping back makes it stale rather than
##' absent. Both entry points here reconfigure the object completely, so the tag
##' and the supply state can never disagree. The cost is that switching clears
##' the solved state and the drivers must be set again — which is honest, since
##' the two paths read different inputs.
##'
##' @param gravity_head head required to lift water to the collar, MPa. Zero by
##'   default, which is right for a bare leaf that is not thinking about rooting
##'   depth. See the section above for why this, alone, is not a driver.
##'
##' @section The two ends of this path are in different unit bases:
##' \strong{Read this before parameterising a whole soil-to-leaf path}, which is
##' exactly what this supply path is for (#56).

##'
##' - `leaf_specific_conductance_max` is **kg** m^-2 s^-1 MPa^-1
##' - [series_resistance()]'s `resistance` is MPa s **mol**^-1 m^2
##'
##' So the two quantities the package presents as the two ends of one series are on
##' opposite sides of a kg-per-mol factor, and a caller carrying the whole path has to
##' supply it. The calibration study recorded dropping that factor of 0.018 as its
##' original error, **worth three orders of magnitude** — which is why it ended up
##' named in their code rather than inlined.
##'
##' The factor is `molar_mass_h2o` = 0.018015 kg/mol. ⚠️ Nothing here can check you
##' applied it: both quantities are just positive numbers, so a path built on the wrong
##' basis returns a plausible operating point, off by ~55x one way or ~0.018x the other.
##'
##' @return A `leaf_supply` object, for the `supply` argument of [leaf_model()]
##'   and [leaf_solve()].
##' @seealso [series_resistance()] for the resistance itself, which is a driver.
##' @examples
##' # a bare leaf: one soil potential, one resistance, no root profile
##' leaf_solve(psi_soil = 1.5, PPFD = 900, supply = leaf_supply_singlelayer(),
##'            root_network = series_resistance(1e3))
##' @export
leaf_supply_singlelayer <- function(gravity_head = 0) {
  out <- list(kind = "single", gravity_head = gravity_head)
  .check_scalars(out["gravity_head"], "leaf_supply_singlelayer")
  if (gravity_head < 0) {
    stop("leaf_supply_singlelayer(): `gravity_head` must be non-negative (MPa)",
         call. = FALSE)
  }
  structure(out, class = c("leaf_supply", "list"))
}

##' A single soil-to-collar series resistance, as a supply driver
##'
##' The single-potential path's counterpart to [root_network_from_carbon()]: it
##' packages one series resistance into the [RootNetwork()] that [set_drivers()]
##' takes, so the two supply paths are driven through the same argument.
##'
##' It goes in `r_R_V_sum` rather than `r_R_H_min`, and that is exact rather than a
##' convention: `r_R_V_sum` already means "series resistance from the surface down
##' to this layer, with no vulnerability weighting", and this is that quantity with
##' one layer. `r_R_H_min` is the vulnerability-weighted horizontal term, which
##' this path does not have — supplying a non-zero one is an error rather than
##' ignored, since it would otherwise silently drop a resistance you meant to use.
##'
##' @param resistance the whole soil-to-collar path as one series resistance,
##'   **per unit leaf area**, MPa s (mol H2O)^-1 m^2 leaf. Must be positive: a
##'   zero resistance is an infinite flux. Per unit leaf area because the leaf is
##'   purely intensive — nothing here may scale with plant size. [set_drivers()]
##'   refuses one below 1e-2, where the collar route can no longer resolve the
##'   flux; the stem routes are the zero-resistance model.
##'
##' @return A [RootNetwork()] with one entry in `r_R_V_sum`.
##' @seealso [leaf_supply_singlelayer()], [root_network_from_carbon()], [set_drivers()]
##' @examples
##' series_resistance(1500)
##' leaf_solve(psi_soil = 1.5, supply = leaf_supply_singlelayer(),
##'            root_network = series_resistance(1500))
##' @export
series_resistance <- function(resistance) {
  .check_scalars(list(resistance = resistance), "series_resistance")
  if (resistance <= 0) {
    stop("series_resistance(): `resistance` must be positive; a zero ",
         "resistance is an infinite flux", call. = FALSE)
  }
  # Not `RootNetwork(r_R_V_sum = resistance)`, whose `Rcpp::wrap` is dear and
  # would be paid per likelihood evaluation by a calibration (test-cost.R).
  # Nor a hand-written `structure()`, which could drift from the C++ struct.
  # Copy-on-write leaves the cached prototype unmutated; test-surface.R asserts it.
  out <- .root_network_prototype()
  out$r_R_V_sum <- resistance
  out
}

##' @rdname leaf_supply_singlelayer
##' @export
leaf_supply_multilayer <- function() {
  structure(list(kind = "multilayer"), class = c("leaf_supply", "list"))
}

##' Build a leaf
##'
##' The recommended way to construct a leaf. [Leaf()] is the raw C++ constructor,
##' with all fifteen arguments positional and no defaults; this splits them into
##' traits and numerical settings, defaults both, assigns the five traits the
##' constructor does not take, and initialises the integrator.
##'
##' The result is a stateful R6 object: set drivers with [set_drivers()], choose a
##' model with `$set_model()` if you want one other than the default, solve with
##' `$optimise()`, then read the operating point off the object or with
##' [operating_point()]. For a single solve you probably want [leaf_solve()]
##' instead, which does all of that in one call.
##'
##' @param traits a [leaf_traits()] object
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_singlelayer()]
##'
##' @return A `Leaf` R6 object.
##' @seealso [leaf_traits()], [leaf_control()], [leaf_supply_singlelayer()],
##'   [set_drivers()], [leaf_solve()]
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$optimise()
##' operating_point(l)
##' @export
leaf_model <- function(traits = leaf_traits(), control = leaf_control(),
                       supply = leaf_supply_multilayer()) {
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  if (!inherits(control, "leaf_control")) {
    stop("`control` must come from leaf_control()", call. = FALSE)
  }
  if (!inherits(supply, "leaf_supply")) {
    stop("`supply` must come from leaf_supply_multilayer() or ",
         "leaf_supply_singlelayer()", call. = FALSE)
  }

  l <- Leaf(
    vcmax_25 = traits$vcmax_25,
    stem_c = traits$stem_c,
    stem_P50 = traits$stem_P50,
    root_c = traits$root_c,
    root_P50 = traits$root_P50,
    TF24_beta2 = traits$TF24_beta2,
    jmax_25 = traits$jmax_25,
    phi_J = traits$phi_J,
    curv_fact_elec_trans = traits$curv_fact_elec_trans,
    curv_fact_colim = traits$curv_fact_colim,
    GSS_tol_abs = control$GSS_tol_abs,
    vulnerability_curve_ncontrol = control$vulnerability_curve_ncontrol,
    ci_abs_tol = control$ci_abs_tol,
    ci_niter = control$ci_niter,
    TF24_cost_scale = traits$TF24_cost_scale
  )
  # The traits the constructor does not take, since plant's bindings pin its
  # arity. Each must be assigned here, or its `leaf_traits()` argument is
  # accepted and silently ignored.
  l$R_d_25 <- traits$R_d_25
  l$JS22_gamma <- traits$JS22_gamma
  l$CMax_a <- traits$CMax_a
  l$CMax_b <- traits$CMax_b
  l$g_min <- traits$g_min
  l$initialize_integrator(control$integration_rule, control$integration_tol)
  # After the integrator, so this does not rely on set_supply_single (which
  # clears the solved state) leaving the integrator tolerance alone.
  if (identical(supply$kind, "single")) {
    l$set_supply_single(supply$gravity_head)
  }
  l
}

# --- default RootNetworks, memoised ------------------------------------------
#
# Materialising a RootNetwork in R is dear, nearly all of it `Rcpp::wrap`; built
# per call it would dominate `set_drivers()` and `leaf_solve()`. test-cost.R
# counts the crossings.
#
# Each memo holds one entry: a sweep holds `soil_depth` fixed, so it hits from the
# second row on and cannot grow. Sharing one object is safe: `set_physiology`
# copies it into C++, and copy-on-write protects it from callers.

.network_memo <- new.env(parent = emptyenv())

# The multi-layer default: 20 kg C m^-2 leaf split evenly over the layers. Keying
# on `soil_depth` alone holds only while the carbon is `rep(20 / n, n)` with
# `n == length(soil_depth)`.
.default_root_network <- function(n, soil_depth) {
  if (!is.null(.network_memo$depth) &&
      identical(.network_memo$depth, soil_depth)) {
    return(.network_memo$net)
  }
  net <- root_network_from_carbon(
    root_carbon_per_leaf_area = rep(20 / n, n),
    soil_depth = soil_depth
  )
  .network_memo$depth <- soil_depth
  .network_memo$net <- net
  net
}

# One default-constructed RootNetwork per session, as a prototype to copy; see
# series_resistance().
.root_network_prototype <- function() {
  if (is.null(.network_memo$proto)) {
    .network_memo$proto <- RootNetwork()
  }
  .network_memo$proto
}

# The single-potential default: a nominal 1e3 MPa s (mol H2O)^-1 m^2 leaf, the
# value the vignettes use. A stand-in, not a recommendation.
.default_series_resistance <- function() {
  if (is.null(.network_memo$series)) {
    .network_memo$series <- series_resistance(1e3)
  }
  .network_memo$series
}

##' Set the drivers for a leaf
##'
##' Named and defaulted, over the C++ `set_physiology()`'s ten positional
##' arguments. Modifies `x` in place and returns it invisibly, so it can be
##' piped.
##'
##' @section Water potentials are positive magnitudes:
##' `psi_soil` is a positive number of MPa. There is one representation for water
##' potential throughout this package and it is asserted at the input boundary
##' rather than documented: a negative entry is an error. That check is the only
##' thing standing between a script written against the old signed convention and
##' a plausible wrong number, so this wrapper surfaces it rather than smoothing
##' it over -- it does not take `abs()`.
##'
##' @section Multiple soil layers:
##' Pass vectors for `psi_soil`, and optionally for `soil_depth`. Both must end up
##' the same length; by default `soil_depth` becomes 1 m layers.
##'
##' @section The supply resistances, on either path:
##' The leaf takes the soil-to-collar hydraulic RESISTANCES, not the root carbon
##' they may have been derived from (#33). The solve reads two vectors of them,
##' `r_R_H_min` and `r_R_V_sum`, and nothing in it knows about root carbon, layer
##' thickness, the 1/3 : 2/3 vertical/horizontal split or either `beta_R_*`
##' constant. This function used to ask for a root carbon profile, which a leaf
##' physiologist does not have, and invented a default of 20 kg C m^-2 leaf.
##'
##' **`root_network` is the same argument on both supply paths**, which is the
##' point. On the multi-layer path it carries one `r_R_H_min` and one `r_R_V_sum`
##' per rooted layer; on the single-potential path it carries one series resistance
##' in `r_R_V_sum` and nothing else, which is that field's own meaning with one
##' layer and no vulnerability-weighted term. [root_network_from_carbon()] builds
##' the first, [series_resistance()] the second. Before this change the
##' single-potential resistance was an argument to [leaf_supply_singlelayer()] instead —
##' so the same quantity arrived at a different *time* depending on which path was
##' in force, and it was the only fitted parameter whose setter reset the object.
##'
##' Both paths default it, because neither should force a caller to own a supply
##' model they do not have, and both defaults are written out in this function's
##' body rather than buried in a signature: a nominal 20 kg C m^-2 leaf through
##' [root_network_from_carbon()] on the multi-layer path, and a nominal 1e3
##' MPa s (mol H2O)^-1 m^2 leaf series resistance on the single-potential one.
##' Both are stand-ins, not recommendations.
##'
##' `soil_depth` is the one argument that is genuinely multi-layer-only, and
##' passing it on the single-potential path is an error rather than ignored: there
##' is no depth profile there for anything to read. That is also why
##' [leaf_supply_singlelayer()] keeps `gravity_head` — see its documentation for the one
##' remaining asymmetry between the two paths, and why it is left in place.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param psi_soil soil water potential, MPa, **positive magnitude**. Length
##'   sets the number of soil layers.
##' @param PPFD INTERCEPTED photosynthetic photon flux density per unit leaf area,
##'   umol m^-2 s^-1: the light incident on the leaf, before its absorptance
##'   (`alpha_PAR_`). On the energy-balance path the leaf also needs `$NIR_int_`,
##'   the intercepted near-infrared in W m^-2, set before this call.
##' @param soil_depth cumulative depth to the bottom of each layer, m. Defaults
##'   to 1 m layers.
##' @param root_network the soil-to-collar hydraulic resistances, a
##'   [RootNetwork()]. **The same argument on both supply paths**: per rooted layer
##'   from [root_network_from_carbon()] on the multi-layer path, or one series
##'   resistance from [series_resistance()] on the single-potential path. Defaults
##'   to a nominal value on either — 20 kg C m^-2 leaf split evenly, or 1e3 —
##'   written out in the body where you can see and replace it. Stand-ins, not
##'   recommendations.
##' @param leaf_specific_conductance_max maximum leaf-specific hydraulic
##'   conductance, kg m^-2 s^-1 MPa^-1
##' @param atm_vpd atmospheric vapour pressure deficit, kPa. Must lie in
##'   `(0, e_s(leaf_temp))`, a relative humidity between 0 and 1; a deficit at or
##'   past saturation is refused.
##' @param ca atmospheric CO2 partial pressure, Pa
##' @param leaf_temp leaf temperature, deg C
##' @param atm_o2_kpa atmospheric O2 partial pressure, kPa
##' @param atm_kpa atmospheric pressure, kPa. Note this is not decorative: the
##'   ppm-to-Pa conversion is derived from it, and getting it wrong moves
##'   assimilation by percent, not by rounding.
##'
##' @return `x`, invisibly.
##' @seealso [leaf_model()], [leaf_solve()]
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900, atm_vpd = 1.5)
##' # three drying layers
##' set_drivers(l, psi_soil = c(1.0, 1.5, 2.0), PPFD = 900)
##' @export
set_drivers <- function(x,
                        psi_soil,
                        PPFD = 900,
                        soil_depth = NULL,
                        root_network = NULL,
                        leaf_specific_conductance_max = 3.14e-5,
                        atm_vpd = 2.0,
                        ca = 40.0,
                        leaf_temp = 25.0,
                        atm_o2_kpa = 21.0,
                        atm_kpa = 101.3) {
  a <- .resolve_drivers(x, psi_soil, PPFD, soil_depth, root_network,
                        leaf_specific_conductance_max, atm_vpd, ca, leaf_temp,
                        atm_o2_kpa, atm_kpa)
  x$set_physiology(a$root_network, a$PPFD, a$psi_soil, a$soil_depth,
                   a$leaf_specific_conductance_max, a$atm_vpd, a$ca,
                   a$leaf_temp, a$atm_o2_kpa, a$atm_kpa)
  invisible(x)
}

# Validate and default the drivers, without applying them. `leaf_gradient()`
# resolves once and re-applies per perturbation, and must call this rather than
# copy it: the defaults here (1 m layers, the nominal networks, the single-path
# placeholder depth) are load-bearing.
.resolve_drivers <- function(x, psi_soil, PPFD, soil_depth, root_network,
                             leaf_specific_conductance_max, atm_vpd, ca,
                             leaf_temp, atm_o2_kpa, atm_kpa) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  psi_soil <- as.numeric(psi_soil)
  n <- length(psi_soil)
  if (n < 1L) {
    stop("`psi_soil` must have at least one layer", call. = FALSE)
  }

  # C++ enforces this too (hazard 2); this check exists to name the fix.
  if (any(psi_soil < 0)) {
    stop("`psi_soil` must be positive magnitudes in MPa, not signed ",
         "potentials. If this came from an older script, drop the minus sign ",
         "rather than negating the result.", call. = FALSE)
  }

  single <- identical(x$supply_kind, "single")

  if (single) {
    # This path reads psi_soil[1] and no depth profile, so extra layers and a
    # `soil_depth` are refused rather than ignored.
    if (n != 1L) {
      stop("`psi_soil` must be a single value on the single-potential supply ",
           "path; got ", n, " layers. Build the leaf with ",
           "leaf_supply_multilayer() for a layered profile.", call. = FALSE)
    }
    if (!is.null(soil_depth)) {
      stop("`soil_depth` is not used on the single-potential supply path -- ",
           "there is no depth profile to read. Drop it, or build the leaf with ",
           "leaf_supply_multilayer().", call. = FALSE)
    }
    # `root_network` is read here too: one series resistance in r_R_V_sum.
    if (is.null(root_network)) {
      root_network <- .default_series_resistance()
    } else if (!inherits(root_network, "RootNetwork")) {
      stop("`root_network` must be a RootNetwork, from series_resistance() or ",
           "RootNetwork(); got ", class(root_network)[[1]], call. = FALSE)
    }
    # set_physiology requires psi_soil and soil_depth to agree in length; this
    # path reads only psi_soil[1], so the depth is a placeholder rather than input.
    soil_depth <- 1.0
  } else {
    soil_depth <- if (is.null(soil_depth)) {
      seq_len(n) * 1.0
    } else {
      as.numeric(soil_depth)
    }
    if (length(soil_depth) != n) {
      stop("`psi_soil` (", n, ") and `soil_depth` (", length(soil_depth),
           ") must have one entry per soil layer", call. = FALSE)
    }

    if (is.null(root_network)) {
      root_network <- .default_root_network(n, soil_depth)
    } else if (!inherits(root_network, "RootNetwork")) {
      stop("`root_network` must be a RootNetwork, from RootNetwork() or ",
           "root_network_from_carbon(); got ", class(root_network)[[1]],
           call. = FALSE)
    }
  }

  # In `set_physiology()`'s argument order: callers pass it positionally.
  list(root_network = root_network,
       PPFD = PPFD,
       psi_soil = psi_soil,
       soil_depth = soil_depth,
       leaf_specific_conductance_max = leaf_specific_conductance_max,
       atm_vpd = atm_vpd,
       ca = ca,
       leaf_temp = leaf_temp,
       atm_o2_kpa = atm_o2_kpa,
       atm_kpa = atm_kpa)
}

##' The solved operating point, as one row
##'
##' Reads the outputs off a leaf that has been solved. Call after `$optimise()`;
##' before that the values are the missing-value sentinels the object was
##' constructed with.
##'
##' Costs about 4 µs, so it is usable in a loop. It was 180 µs until it stopped
##' reading the fifteen outputs through fifteen separate calls into C++ and
##' stopped building its one row with `data.frame()` (#39) -- 45× more than the
##' ~3 µs solve it was reporting on.
##'
##' @section What the cost columns mean:
##'
##' `hydraulic_cost` is the whole cost the objective subtracted. `shadow_cost` is
##' the share of it that is a *price* rather than carbon the plant gave up: the
##' value of water in its best alternative use, which for a leaf is assimilation
##' later. So
##'
##' \preformatted{
##'   realised cost   = hydraulic_cost - shadow_cost
##'   carbon profit   = profit + shadow_cost
##' }
##'
##' and a consumer that grows a plant on this leaf wants the second of those,
##' never `profit`. Deducting a shadow price from a carbon budget taxes growth by
##' something the plant never spent.
##'
##' `shadow_cost` is `TF24_floor_lambda_o * E` on that curve and zero on every
##' other. Zero means *this curve does not separate the two*, not *this curve's
##' cost is all realised carbon*. `CF77` is the case worth stating: its whole cost
##' is `lambda * E`, and reading that as a shadow price is an ordinary thing to do,
##' but the model supplies one number and nothing to attribute it with and behaves
##' identically whichever way it is read. Reporting it as all-shadow would present
##' one reading as a fact. The field is defined only where the curve defines it.
##'
##' @param x a solved `Leaf`
##' @return A one-row data.frame.
##' @seealso [leaf_solve()], which does the whole thing in one call.
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$optimise()
##' operating_point(l)
##' @export
operating_point <- function(x) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  # Built directly: data.frame() costs far more than the solve this reports.
  # test-surface.R asserts the result is `identical()` to data.frame()'s, which
  # needs row.names to be the integer 1L, not 1.0.
  v <- as.list(x$operating_point_values())
  names(v) <- .operating_point_names
  structure(v, class = "data.frame", row.names = 1L)
}

# The names of the values `Leaf::operating_point_values()` returns, in its order,
# shared by `operating_point()` and `leaf_solve()`. The order is an interface the
# types do not enforce: append, never insert, or saved outputs shift. test-surface.R
# checks each against its own binding; test-golden.R checks leaf_solve() bit-exactly.
.operating_point_names <- c(
  "psi_stem",       # MPa, positive magnitude
  "collar",         # MPa, positive magnitude
  "ci",             # Pa
  "A",              # umol C m^-2 s^-1
  "E",              # kg H2O m^-2 s^-1
  "gc",             # mol CO2 m^-2 s^-1
  "profit",         # umol C m^-2 s^-1
  "hydraulic_cost",
  "E_up",
  "uptake",
  # TF24's marginal cost of water at this operating point, whatever curve is
  # seated; the seated curve's own is `lambda_emergent`.
  "lambda",         # dA/dE under the TF24 cost
  "g1_eff",         # the Medlyn g1 this leaf implies
  # deg C. Not the `leaf_temp` driver: on the energy-balance path that is air
  # temperature and this the leaf's own; off it, this is a copy of the driver.
  "Tleaf",
  # umol C m^-2 s^-1: the share of `hydraulic_cost` that is a price, not carbon
  # given up; `TF24_floor_lambda_o * E` on that curve and zero on every other,
  # where zero means the curve does not separate the two. See operating_point().
  "shadow_cost",
  # umol C (kg H2O)^-1: the seated curve's marginal cost of water,
  # `(dC/dpsi)/(dE/dpsi)`, on the same axis for every curve.
  "lambda_emergent",
  # kg H2O m^-2 s^-1: the residual (cuticular) flux. `E` is the total flux; the
  # stomatal flux is `E - residual_flux`.
  "residual_flux"
)

##' Solve a leaf, in one call
##'
##' Drivers in, operating point out. Vectorised over the scalar drivers, so a
##' response curve is one call rather than a loop.
##'
##' This is the entry point for someone who would otherwise reach for
##' `plantecophys::Photosyn()`. The difference in what comes back is the point of
##' the package: `gc` here is not from a fitted conductance model, it is what
##' falls out of maximising profit over an explicit hydraulic path, and `g1_eff`
##' reports the Medlyn `g1` that would have been needed to reproduce it.
##'
##' @section `leaf_temp` is a driver, `Tleaf` is an output:
##' Two temperature columns come back and they are not the same thing.
##' `leaf_temp` is the driver you supplied. `Tleaf` is the temperature the
##' reported fluxes were computed at.
##'
##' By default they are equal, because the default path takes leaf temperature as
##' prescribed. With the Penman-Monteith energy balance on
##' (`$use_energy_balance_ <- TRUE`) they are not: `leaf_temp` is then **air**
##' temperature, and `Tleaf` is solved per operating point from the leaf's own
##' transpiration — so it rises as the stomata close, and it is the column to plot
##' assimilation against.
##'
##' @section Performance, and what the lever actually is:
##' About **20 µs per row**, of which the model is **2.8 µs**. The rest is the R
##' boundary: each call into C++ costs ~1.1 µs and a solved row needs several.
##'
##' Two consequences, both the opposite of what you might expect:
##'
##' * **This is the fast path, not the convenient-but-slow one.** It is within 6%
##'   of building a `Leaf` yourself and looping `set_drivers()` +
##'   `$find_root_collar_psi()` + [operating_point()]. Reaching into the object
##'   saves about 1 µs a row and is worth doing for access to intermediate state,
##'   not for speed. (This was not true before: building the result row by row and
##'   `rbind`ing cost 344 µs a row, 26× the stateful loop.)
##' * **The lever is fewer R calls per row, not a different R function.** Pass all
##'   your drivers to one vectorised call rather than looping in R. If you need
##'   more than that -- a fit's inner loop, say -- what has to go is the boundary
##'   itself, which means C++.
##'
##' `reuse = TRUE` is the default because constructing a `Leaf` from R costs
##' about **180x a trivial `.Call`** (~230 µs where this was measured), or 45
##' solves, and only a small part of that is the vulnerability tables; the rest is
##' R-side object construction. See [set_traits()] for varying traits without
##' reconstructing, and note it costs many solves rather than almost nothing when
##' the trait you change is `stem_c` or `root_c`, which own a table.
##'
##' @section Recycling:
##' `PPFD`, `atm_vpd`, `ca`, `leaf_temp`, `atm_o2_kpa`, `atm_kpa` and
##' `leaf_specific_conductance_max` are recycled to a common length, the usual R
##' way. `psi_soil` is recycled too **when it is a plain numeric vector** --
##' each element is then a separate single-layer solve. For multi-layer solves
##' pass a list of numeric vectors, one per row; `soil_depth` follows the same
##' rule. `root_network` is either one [RootNetwork()], used for every row, or a
##' list of them, one per row.
##'
##' @inheritParams set_drivers
##' @param traits a [leaf_traits()] object
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_singlelayer()]. On the single-potential path
##'   `soil_depth` and `root_network` must be omitted, and each `psi_soil` is one
##'   value rather than a profile.
##' @param model which optimality model to solve, as in [leaf_gradient()]:
##'   `"collar"` (the default, the production path -- the TF24 cost maximised over
##'   the root-collar potential) or any name from [cost_curve_names()], which
##'   solves that curve's **stem** optimum instead: `psi_stem` with the upstream
##'   potential pinned at `psi_soil`, which needs `leaf_supply_singlelayer()`.
##'
##'   ⚠️ **`"collar"` and `"TF24"` are different models, not two ways of asking
##'   one question**, even though both use the TF24 cost. They optimise different
##'   variables over different supply topologies and disagree on 20 of 30 driver
##'   rows.
##'
##'   ⚠️ **On a stem route the `collar` column comes back non-finite**, because
##'   the collar is not solved for at all. An absent number is the truthful answer
##'   rather than a gap -- reading one would be reading whatever the last solve
##'   left behind.
##'
##'   Two curves take a PRESCRIBED price of water rather than deriving one:
##'   `CF77` reads `CF77_lambda` and `TF24_floor` reads `TF24_floor_lambda_o`,
##'   both below. `CF77`'s has a default, aligned to TF24 by [align_curves()];
##'   `TF24_floor`'s has none, and that curve refuses to solve without one.
##' @param LeastCost_beta the `LeastCost` curve's unit-cost ratio (Prentice et
##'   al. 2014's beta): the cost of maintaining carboxylation capacity relative to
##'   the cost of transpiration, dimensionless. Only `model = "LeastCost"` reads
##'   it. `NA_real_` (the default) keeps the leaf's own default, aligned to TF24
##'   by [align_curves()]. Zero is refused, since at zero the capacity term
##'   vanishes and what is left is water-use efficiency, which rises without bound
##'   as transpiration falls. The published value is 146, which at this
##'   package's default leaf pins the stem optimum at `psi_crit` from the
##'   reference `psi_soil` of 1 MPa drier.
##' @param CF77_lambda,TF24_floor_lambda_o the prescribed marginal value of water
##'   for the two curves that take one, in umol C (kg H2O)^-1. `NA_real_` (the
##'   default) leaves the field as it is: `CF77`'s aligned default (see
##'   [align_curves()]), or unset for `TF24_floor`, which refuses to solve.
##'
##'   ⚠️ **Their scale is set by the leaf, not chosen freely.** At this package's
##'   defaults the leaf's own marginal cost of water runs 9e4 to 3e5, and a price
##'   far outside that band pins the optimum against a bound — where the answer
##'   describes the bracket rather than the model. [operating_point()]'s
##'   `lambda_emergent` is how you find a value that means something.
##' @param reuse solve every row with one `Leaf` object (`TRUE`) or construct a
##'   fresh one per row (`FALSE`, the default). Reuse is faster because the
##'   vulnerability tables are built once, and it is safe -- every exit from the
##'   solve writes its own outputs. `FALSE` exists so that can be checked rather
##'   than assumed, and the package's tests do check it.
##'
##' @return A data.frame with one row per driver combination: the drivers, then
##'   the columns of [operating_point()].
##' @seealso [leaf_model()] for the stateful interface.
##' @examples
##' # one point
##' leaf_solve(psi_soil = 2.0, PPFD = 900)
##'
##' # a light response curve
##' leaf_solve(psi_soil = 1.0, PPFD = seq(100, 1800, length.out = 6))
##'
##' # a drought response, and the effective g1 it implies
##' d <- leaf_solve(psi_soil = seq(0.5, 5, length.out = 6), PPFD = 900)
##' d[, c("psi_soil", "A", "gc", "g1_eff")]
##'
##' # two three-layer profiles
##' leaf_solve(psi_soil = list(c(1, 1.5, 2), c(3, 3.5, 4)), PPFD = 900)
##'
##' # a different optimality model, on the single-potential path it needs
##' leaf_solve(psi_soil = 1.5, PPFD = 900, model = "SOX",
##'            supply = leaf_supply_singlelayer())
##' @export
leaf_solve <- function(psi_soil,
                       PPFD = 900,
                       soil_depth = NULL,
                       root_network = NULL,
                       leaf_specific_conductance_max = 3.14e-5,
                       atm_vpd = 2.0,
                       ca = 40.0,
                       leaf_temp = 25.0,
                       atm_o2_kpa = 21.0,
                       atm_kpa = 101.3,
                       traits = leaf_traits(),
                       control = leaf_control(),
                       supply = leaf_supply_multilayer(),
                       reuse = TRUE,
                       model = "collar",
                       CF77_lambda = NA_real_,
                       TF24_floor_lambda_o = NA_real_,
                       LeastCost_beta = NA_real_) {
  .check_model(model, supply)
  # One price for every row: the price is seated once per Leaf.
  long <- c(CF77_lambda = length(CF77_lambda),
            TF24_floor_lambda_o = length(TF24_floor_lambda_o),
            LeastCost_beta = length(LeastCost_beta)) != 1L
  if (any(long)) {
    stop("`", names(long)[long][[1]], "` must be a single number (or NA): ",
         "leaf_solve() seats one price for every row. For a price per ",
         "observation, use leaf_batch().", call. = FALSE)
  }
  prices <- .check_prices(model, CF77_lambda, TF24_floor_lambda_o,
                          LeastCost_beta)
  layered <- .as_layer_list(psi_soil, "psi_soil")
  scalars <- list(PPFD = PPFD, atm_vpd = atm_vpd, ca = ca,
                  leaf_temp = leaf_temp, atm_o2_kpa = atm_o2_kpa,
                  atm_kpa = atm_kpa,
                  leaf_specific_conductance_max = leaf_specific_conductance_max)

  n <- max(length(layered), vapply(scalars, length, integer(1)))
  if (n == 0L) {
    stop("nothing to solve: every driver has length zero", call. = FALSE)
  }
  layered <- .recycle_to(layered, n, "psi_soil")
  scalars <- lapply(seq_along(scalars), function(i) {
    .recycle_to(scalars[[i]], n, names(scalars)[[i]])
  })
  names(scalars) <- c("PPFD", "atm_vpd", "ca", "leaf_temp", "atm_o2_kpa",
                      "atm_kpa", "leaf_specific_conductance_max")

  depths <- if (is.null(soil_depth)) NULL else
    .recycle_to(.as_layer_list(soil_depth, "soil_depth"), n, "soil_depth")
  # One network for every row, or one per row. inherits() first: a RootNetwork
  # is a list of five, and would otherwise read as five rows' networks.
  networks <- if (is.null(root_network)) {
    NULL
  } else if (inherits(root_network, "RootNetwork")) {
    rep(list(root_network), n)
  } else {
    .recycle_to(root_network, n, "root_network")
  }

  shared <- if (reuse) .seat_model(leaf_model(traits, control, supply), model,
                                   prices)
            else NULL

  # Nothing allocates per row, which is most of this function's speed: drivers
  # are already length-n vectors, outputs fill one matrix, and the one
  # data.frame() call, at the end, is a per-call cost.
  outputs <- matrix(NA_real_, nrow = n, ncol = length(.operating_point_names),
                    dimnames = list(NULL, .operating_point_names))

  # One condition handler per call, outside the loop, so classifying a failure
  # costs nothing per row (test-cost.R). The first infeasible row still aborts
  # the call; see with_phylloptim_conditions().
  with_phylloptim_conditions(
  for (i in seq_len(n)) {
    l <- if (reuse) shared else .seat_model(leaf_model(traits, control, supply),
                                            model, prices)
    set_drivers(l,
                psi_soil = layered[[i]],
                PPFD = scalars$PPFD[[i]],
                soil_depth = if (is.null(depths)) NULL else depths[[i]],
                root_network =
                  if (is.null(networks)) NULL else networks[[i]],
                leaf_specific_conductance_max =
                  scalars$leaf_specific_conductance_max[[i]],
                atm_vpd = scalars$atm_vpd[[i]],
                ca = scalars$ca[[i]],
                leaf_temp = scalars$leaf_temp[[i]],
                atm_o2_kpa = scalars$atm_o2_kpa[[i]],
                atm_kpa = scalars$atm_kpa[[i]])
    l$optimise()
    outputs[i, ] <- l$operating_point_values()
  })

  # psi_soil reports the wettest layer, with the layer count alongside, so the
  # column stays numeric and the frame rectangular.
  data.frame(psi_soil = vapply(layered, min, numeric(1)),
             layers = lengths(layered),
             PPFD = scalars$PPFD,
             atm_vpd = scalars$atm_vpd,
             ca = scalars$ca,
             leaf_temp = scalars$leaf_temp,
             atm_kpa = scalars$atm_kpa,
             outputs)
}

# --- internals ---------------------------------------------------------------

# Validate a model name against the supply path, before a Leaf is built. The
# names come from `.gradient_models()`, which `leaf_gradient()` also uses, so a
# curve added in C++ is solvable here and the two entry points cannot disagree.
.check_model <- function(model, supply) {
  if (!(is.character(model) && length(model) == 1L)) {
    stop("`model` must be a single model name; see cost_curve_names()",
         call. = FALSE)
  }
  if (!(model %in% .gradient_models())) {
    stop("unknown model \"", model, "\". Available: ",
         paste(.gradient_models(), collapse = ", "), ".", call. = FALSE)
  }
  if (!identical(model, "collar") && !identical(supply$kind, "single")) {
    stop("the stem routes optimise `psi_stem` with the upstream potential ",
         "pinned at `psi_soil`, so they need the single-potential supply path; ",
         "got the multi-layer one. Use model = \"collar\" for a root network.",
         call. = FALSE)
  }
  invisible(model)
}

# The model and its prescribed prices are configuration that survives
# set_drivers() and set_traits(), so this runs once per Leaf, not per row. It is
# `leaf_solve()`'s only route to the price fields, which `leaf_traits()` cannot carry.
.seat_model <- function(x, model, prices = NULL) {
  x$set_model(if (identical(model, "collar")) "TF24" else model,
              if (identical(model, "collar")) "collar" else "stem")
  for (nm in names(prices)) {
    if (!is.na(prices[[nm]])) {
      x[[nm]] <- prices[[nm]]
    }
  }
  x
}

# Which field belongs to which curve: the solve's twin of `.gradient_owned_pars`,
# so a new priced curve needs a row in both. A price the seated model does not
# read is refused, not ignored; an omitted one is reported by the curve itself.
.solve_price_fields <- c(CF77 = "CF77_lambda_",
                         TF24_floor = "TF24_floor_lambda_o",
                         LeastCost = "LeastCost_beta")

# The defaults of the owned fields that have one: solved by align_curves(), and
# the R copy of the in-class initialisers in leaf.hpp. test-align-curves.R
# checks both against a fresh Leaf and a re-run of the alignment.
# `TF24_floor_lambda_o` has none.
.level_defaults <- c(CF77_lambda_ = 107750.8846,
                     LeastCost_beta = 51.59667042)

.owns_field <- function(model, field) {
  model %in% names(.solve_price_fields) &&
    identical(.solve_price_fields[[model]], field)
}

# The value an owned field takes in a gradient's theta: off the leaf when there
# is one, its default when there is not, and NA on every route that does not own
# it. The NA is what keeps the R reference from writing the field once per
# perturbation on routes that never read it -- a crossing `test-cost.R` counts.
.owned_value <- function(model, field, x = NULL) {
  if (!.owns_field(model, field)) {
    return(NA_real_)
  }
  if (!is.null(x)) {
    return(x[[field]])
  }
  if (field %in% names(.level_defaults)) .level_defaults[[field]] else NA_real_
}

# `LeastCost_beta` is a dimensionless ratio of unit costs, not a price of water,
# so the refusal names it differently. What the three fields share is that one
# curve reads each and no other does.
.check_prices <- function(model, CF77_lambda, TF24_floor_lambda_o,
                          LeastCost_beta = NA_real_) {
  # A list, because leaf_gradient_batch() passes a column of `theta`, one value
  # per row; a field is given when any row sets it.
  given <- list(CF77_lambda_ = CF77_lambda,
                TF24_floor_lambda_o = TF24_floor_lambda_o,
                LeastCost_beta = LeastCost_beta)
  # Not `[[model]]`, which errors on an unpriced model; "" matches no field.
  mine <- if (model %in% names(.solve_price_fields)) {
    .solve_price_fields[[model]]
  } else {
    ""
  }
  set <- vapply(given, function(x) any(!is.na(x)), logical(1))
  wrong <- names(given)[set & names(given) != mine]
  if (length(wrong)) {
    owner <- names(.solve_price_fields)[match(wrong[[1]], .solve_price_fields)]
    what <- if (wrong[[1]] == "LeastCost_beta") {
      paste("the LeastCost curve's unit-cost ratio, the cost of maintaining",
            "carboxylation capacity relative to that of transpiration")
    } else {
      paste("the", owner, "curve's prescribed price of water")
    }
    stop("`", sub("_$", "", wrong[[1]]), "` is ", what,
         ", and model = \"", model,
         "\" does not read it: on every other curve the price is emergent, ",
         "derived from that curve's own parameters rather than set. Pass ",
         "model = \"", owner, "\", or drop the argument.", call. = FALSE)
  }
  unlist(given)
}


.check_scalars <- function(x, what) {
  bad <- names(x)[!vapply(x, function(v) {
    is.numeric(v) && length(v) == 1L && is.finite(v)
  }, logical(1))]
  if (length(bad)) {
    stop(what, "(): each argument must be a single finite number; ",
         "the following ", if (length(bad) == 1L) "is" else "are", " not: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Normalise a driver that may describe soil layers into a list of numeric
# vectors, one per row. A plain numeric vector means "one single-layer solve per
# element"; a list means "one multi-layer solve per element".
.as_layer_list <- function(x, what) {
  if (is.list(x)) {
    out <- lapply(x, as.numeric)
  } else {
    out <- as.list(as.numeric(x))
  }
  if (any(vapply(out, length, integer(1)) == 0L)) {
    stop("`", what, "` has an entry with no layers", call. = FALSE)
  }
  out
}

.recycle_to <- function(x, n, what) {
  len <- length(x)
  if (len == n) {
    return(x)
  }
  if (len == 1L) {
    return(rep(x, n))
  }
  stop("`", what, "` has length ", len,
       ", which cannot be recycled to ", n,
       ". Driver lengths must be 1 or the longest driver's length.",
       call. = FALSE)
}
