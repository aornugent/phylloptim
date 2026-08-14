# The friendly R surface (issue #5, stage 2).
#
# Everything here is hand-written and sits ABOVE the generated glue in
# R/RcppR6.R, which is what makes it rewritable without touching the build, the
# YAML or CI. Nothing in this file reaches into C++ except through the generated
# `Leaf` R6 class.
#
# Two things it is for:
#
#   1. Separating traits from tolerances. The C++ constructor takes 17 positional
#      arguments and four of them -- GSS_tol_abs, ci_abs_tol, ci_niter,
#      vulnerability_curve_ncontrol -- are numerical settings sitting among the
#      physiology. A trait-calibration loop should not have to know which of 19
#      arguments are not traits. `leaf_traits()` and `leaf_control()` split them.
#   2. A one-call entry point. `leaf_solve()` takes drivers and returns a
#      data.frame, which is the shape someone reaching for plantecophys expects.

# The C++ default constructor's values, in one place. These are Eucalyptus
# saligna, and they are the defaults for every function below.
#
# ⚠️ These must equal phylloptim::Leaf's default constructor, and it is not enough to
# believe so: tests/testthat/test-golden.R solves through THESE defaults and
# compares bit-exactly against a golden file generated from a default-constructed
# C++ Leaf. If the two ever drift apart, that file fails.
#
# root_c / root_b / root_psi_crit default in MultiLayerRoots rather than in Leaf,
# because a second copy of the root Weibull pair inside Leaf is exactly hazard 1
# in the developer guide -- it has already cost this project a wrong exponent in a
# manuscript draft. They are restated here because the R constructor has to pass
# them, not because Leaf owns them.
#
# beta_R_H / beta_R_V are NOT here since #33: they parameterise the
# root-architecture model, which the leaf no longer runs. Their defaults live with
# that model, on root_network_from_carbon().
.leaf_trait_defaults <- list(
  vcmax_25 = 96,
  stem_c = 2.680147,
  stem_b = 3.898245,
  psi_crit = 5.870283,
  root_c = 2.680147,
  root_b = 3.898245,
  root_psi_crit = 5.870283,
  beta2 = 1.5,
  jmax_25 = 157.44,
  a = 0.30,
  curv_fact_elec_trans = 0.7,
  curv_fact_colim = 0.99,
  cost_scale_TF24 = 7.5,
  R_d_25 = 1.44
)

.leaf_control_defaults <- list(
  GSS_tol_abs = 1e-3,
  vulnerability_curve_ncontrol = 100,
  ci_abs_tol = 1e-3,
  ci_niter = 1000,
  integration_rule = 21,
  integration_tol = 1e-8
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
##' `stem_b`/`stem_c` describe the STEM Weibull curve, which drives the
##' hydraulic cost; `root_b`/`root_c` describe the ROOT curve, which drives
##' uptake. They are separate parameters with separate meanings and they are
##' named accordingly, because they were once an unmarked `b`/`c` pair alongside
##' `root_b`/`root_c` and an analysis used the root parameters for the stem cost.
##'
##' @param vcmax_25 maximum carboxylation rate at 25 C (umol m^-2 s^-1)
##' @param stem_c shape parameter of the stem vulnerability curve (unitless)
##' @param stem_b sensitivity parameter of the stem vulnerability curve (MPa)
##' @param psi_crit critical stem water potential (MPa, positive magnitude)
##' @param root_c shape parameter of the root vulnerability curve (unitless)
##' @param root_b sensitivity parameter of the root vulnerability curve (MPa)
##' @param root_psi_crit critical root water potential (MPa, positive magnitude)
##' @param beta2 exponent for the effect of hydraulic risk (unitless)
##' @param jmax_25 maximum electron transport rate at 25 C (umol m^-2 s^-1)
##' @param a quantum yield of photosynthetic electron transport (mol mol^-1)
##' @param curv_fact_elec_trans curvature of the light response curve (unitless)
##' @param curv_fact_colim curvature of the colimited photosynthesis equation
##' @param cost_scale_TF24 cost parameter for the TF24 profit model
##'   (umol m^-2 s^-1)
##' @param R_d_25 dark respiration at 25 C (umol m^-2 s^-1). It is the value at
##'   25 C only: respiration rises from there on Tjoelker's declining-Q10 curve.
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
##' leaf_traits(vcmax_25 = 120, stem_b = 2.5)
##' @export
leaf_traits <- function(vcmax_25 = 96,
                        stem_c = 2.680147,
                        stem_b = 3.898245,
                        psi_crit = 5.870283,
                        root_c = 2.680147,
                        root_b = 3.898245,
                        root_psi_crit = 5.870283,
                        beta2 = 1.5,
                        jmax_25 = 157.44,
                        a = 0.30,
                        curv_fact_elec_trans = 0.7,
                        curv_fact_colim = 0.99,
                        cost_scale_TF24 = 7.5,
                        R_d_25 = 1.44) {
  out <- list(vcmax_25 = vcmax_25, stem_c = stem_c, stem_b = stem_b,
              psi_crit = psi_crit, root_c = root_c, root_b = root_b,
              root_psi_crit = root_psi_crit, beta2 = beta2,
              jmax_25 = jmax_25, a = a,
              curv_fact_elec_trans = curv_fact_elec_trans,
              curv_fact_colim = curv_fact_colim,
              cost_scale_TF24 = cost_scale_TF24,
              R_d_25 = R_d_25)
  .check_scalars(out, "leaf_traits")
  if (R_d_25 < 0) {
    stop("leaf_traits(): R_d_25 must be non-negative", call. = FALSE)
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
##' The internal leaf-temperature clamp (`-40 C` to `70 C`) is a compile-time
##' constant and stays one. It is a guard that keeps an extreme non-equilibrium
##' transpiration on the energy-balance path from driving the Arrhenius block
##' non-finite -- not a tolerance anyone tunes. Making it settable would let a
##' caller produce non-finite photosynthetic parameters and get NaNs back with no
##' indication of why.
##'
##' @param GSS_tol_abs absolute tolerance for the golden-section search over
##'   stem water potential. ⚠️ **This no longer sets how well the operating point
##'   is determined**, and this text used to say that it did. The collar solve now
##'   solves its own first-order condition to about `1e-12` instead of searching
##'   profit to `GSS_tol_abs`, so changing this leaves the answer bit-identical on
##'   the production path. What it still does: it is the width below which the
##'   feasible collar interval is treated as a single point, and it is the
##'   tolerance of the two single-layer `optimise_psi_stem_*` optimisers, which
##'   are off that path.
##' @param vulnerability_curve_ncontrol number of control points used to
##'   pre-integrate the two Weibull vulnerability curves into splines. Higher is
##'   more accurate and slower to construct; it does not affect solve speed.
##' @param ci_abs_tol absolute tolerance for root-solving intercellular CO2
##' @param ci_niter maximum iterations for root-solving intercellular CO2
##' @param integration_rule accepted and ignored. It selected a Gauss-Kronrod
##'   rule order in plant's compiled quadrature; the header-only adaptive
##'   Simpson that replaced it has no rule order. Retained because plant's R
##'   bindings expose it.
##' @param integration_tol absolute tolerance for the direct quadrature in
##'   `transpiration_full_integration()`. Not on the production path, which
##'   reads the pre-integrated spline instead.
##'
##' @return A `leaf_control` object; a named list.
##' @seealso [leaf_traits()], [leaf_model()]
##' @examples
##' leaf_control()
##' leaf_control(GSS_tol_abs = 1e-5)
##' @export
leaf_control <- function(GSS_tol_abs = 1e-3,
                         vulnerability_curve_ncontrol = 100,
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
##' `leaf_supply_single()` collapses the whole soil-to-collar path to **one
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
##' precisely the caller `leaf_supply_single()` exists for. Making it a driver here
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
##' @return A `leaf_supply` object, for the `supply` argument of [leaf_model()]
##'   and [leaf_solve()].
##' @seealso [series_resistance()] for the resistance itself, which is a driver.
##' @examples
##' # a bare leaf: one soil potential, one resistance, no root profile
##' leaf_solve(psi_soil = 1.5, PPFD = 900, supply = leaf_supply_single(),
##'            root_network = series_resistance(1e3))
##' @export
leaf_supply_single <- function(gravity_head = 0) {
  out <- list(kind = "single", gravity_head = gravity_head)
  .check_scalars(out["gravity_head"], "leaf_supply_single")
  if (gravity_head < 0) {
    stop("leaf_supply_single(): `gravity_head` must be non-negative (MPa)",
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
##'   purely intensive — nothing here may scale with plant size.
##'
##' @return A [RootNetwork()] with one entry in `r_R_V_sum`.
##' @seealso [leaf_supply_single()], [root_network_from_carbon()], [set_drivers()]
##' @examples
##' series_resistance(1500)
##' leaf_solve(psi_soil = 1.5, supply = leaf_supply_single(),
##'            root_network = series_resistance(1500))
##' @export
series_resistance <- function(resistance) {
  .check_scalars(list(resistance = resistance), "series_resistance")
  if (resistance <= 0) {
    stop("series_resistance(): `resistance` must be positive; a zero ",
         "resistance is an infinite flux", call. = FALSE)
  }
  # ⚠️ NOT `RootNetwork(r_R_V_sum = resistance)`, WHICH COSTS ~58 us. That call
  # reaches C++ for a default-constructed network and pays `Rcpp::wrap` building the
  # five-element named list -- measured at 58 us against 1.1 us for a trivial
  # `.Call`. This function is reached per likelihood evaluation by a calibration
  # whose fitted resistance changes every proposal, so it cannot be hoisted out of
  # the caller's loop the way a fixed network can.
  #
  # Instead: take ONE default-constructed network per session and overwrite the one
  # field. R's copy-on-write means the assignment copies rather than mutating the
  # cached prototype, so callers cannot corrupt each other -- asserted in
  # test-surface.R. The field list still comes from the real constructor, so it
  # cannot drift from the C++ struct the way a hand-written `structure()` would.
  out <- .root_network_prototype()
  out$r_R_V_sum <- resistance
  out
}

##' @rdname leaf_supply_single
##' @export
leaf_supply_multilayer <- function() {
  structure(list(kind = "multilayer"), class = c("leaf_supply", "list"))
}

##' Build a leaf
##'
##' The recommended way to construct a leaf. [Leaf()] is the raw C++ constructor,
##' with all seventeen arguments positional and no defaults; this splits them into
##' traits and numerical settings, defaults both, and initialises the integrator.
##'
##' The result is a stateful R6 object: set drivers with [set_drivers()], solve
##' with `$find_root_collar_psi()`, then read the operating point off the object
##' or with [operating_point()]. For a single solve you probably want
##' [leaf_solve()] instead, which does all of that in one call.
##'
##' @param traits a [leaf_traits()] object
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_single()]
##'
##' @return A `Leaf` R6 object.
##' @seealso [leaf_traits()], [leaf_control()], [leaf_supply_single()],
##'   [set_drivers()], [leaf_solve()]
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$find_root_collar_psi()
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
         "leaf_supply_single()", call. = FALSE)
  }

  # Positional, because that is what the generated constructor takes. The whole
  # value of this function is that the ordering is written down once, here,
  # instead of at every call site.
  l <- Leaf(
    vcmax_25 = traits$vcmax_25,
    stem_c = traits$stem_c,
    stem_b = traits$stem_b,
    psi_crit = traits$psi_crit,
    root_c = traits$root_c,
    root_b = traits$root_b,
    root_psi_crit = traits$root_psi_crit,
    beta2 = traits$beta2,
    jmax_25 = traits$jmax_25,
    a = traits$a,
    curv_fact_elec_trans = traits$curv_fact_elec_trans,
    curv_fact_colim = traits$curv_fact_colim,
    GSS_tol_abs = control$GSS_tol_abs,
    vulnerability_curve_ncontrol = control$vulnerability_curve_ncontrol,
    ci_abs_tol = control$ci_abs_tol,
    ci_niter = control$ci_niter,
    cost_scale_TF24 = traits$cost_scale_TF24
  )
  # ⚠️ AFTER construction, because plant's RcppR6 bindings pin the generated
  # constructor by arity so R_d_25 cannot be an argument to it. Without this line
  # `leaf_traits(R_d_25 = )` would be accepted and silently ignored.
  l$R_d_25 <- traits$R_d_25
  l$initialize_integrator(control$integration_rule, control$integration_tol)
  # After the integrator, because set_supply_single clears the solved state --
  # not the integrator tolerance, but relying on that ordering would be a
  # dependency on an implementation detail of setup_clean_leaf.
  if (identical(supply$kind, "single")) {
    l$set_supply_single(supply$gravity_head)
  }
  l
}

# ---------------------------------------------------------------------------
# The R-boundary tax on a RootNetwork, and why these two are memoised
# ---------------------------------------------------------------------------
#
# ⚠️ MATERIALISING A RootNetwork AS AN R OBJECT COSTS ~58 us, AND ALMOST NONE OF
# IT IS OUR CODE. Measured against the installed package: `RootNetwork__ctor()` --
# the bare `.Call`, before any R wrapper -- is 58 us, where a trivial `.Call` on a
# Leaf field is 1.1 us and `$set_physiology()` handed a ready-made network is
# 4.3 us. It is `Rcpp::wrap` building the five-element named list, which RcppR6
# generates as five successive `ret["name"] = ...` assignments on an initially
# empty `Rcpp::List` -- each one grows the list and its names attribute.
#
# Building the default network per call therefore made `set_drivers()` **9.4 ->
# 73 us** on the multi-layer path and **5.8 -> 66 us** on the single-potential one,
# and `leaf_solve()` **23.7 -> 108 us/row** -- a 4.5x regression on a figure this
# package's own documentation quotes. Caught only by measuring; every test passed.
#
# This is the R boundary dominating a loop again, which is the standing lesson
# here: count R calls, not arithmetic.
#
# Both memos are deliberately SIZE ONE rather than keyed caches. A `leaf_solve()`
# sweep holds `soil_depth` fixed across rows, so a one-entry memo hits on every row
# after the first, and it cannot grow without bound the way a keyed cache would if
# someone swept the soil profile. The comparison is `identical()` on the depth
# vector -- about 1 us -- not a pasted string key, because formatting a numeric
# vector to build one costs more than it saves.
#
# Sharing one R object across calls is safe: `set_physiology` copies it into C++
# through `Rcpp::as`, nothing mutates it, and R's copy-on-write means a caller who
# got hold of it could not mutate ours either.

.network_memo <- new.env(parent = emptyenv())

# The multi-layer default: 20 kg C m^-2 leaf split evenly over the layers, through
# the architecture model the leaf used to run internally (#33) -- so the numbers
# are the pre-#33 ones and the provenance is on the page rather than buried in a
# signature. `soil_depth` determines the whole thing, because the default carbon is
# `rep(20 / n, n)` with `n == length(soil_depth)`, so it is the only key needed.
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

# The single-potential path's default: one nominal series resistance, in the same
# memo and for the same reason. 1e3 MPa s (mol H2O)^-1 m^2 leaf is the value the
# package's own vignettes and the companion calibration study use, and like the
# multi-layer default it is a stand-in rather than a recommendation -- but it means
# `leaf_solve(psi_soil = 2, supply = leaf_supply_single())` means something, exactly
# as it does on the other path. Neither path forces a caller to own a supply model.
# One default-constructed RootNetwork per session, as a prototype to copy. See
# series_resistance() for why, and .network_memo above for the measurement.
.root_network_prototype <- function() {
  if (is.null(.network_memo$proto)) {
    .network_memo$proto <- RootNetwork()
  }
  .network_memo$proto
}

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
##' single-potential resistance was an argument to [leaf_supply_single()] instead —
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
##' [leaf_supply_single()] keeps `gravity_head` — see its documentation for the one
##' remaining asymmetry between the two paths, and why it is left in place.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param psi_soil soil water potential, MPa, **positive magnitude**. Length
##'   sets the number of soil layers.
##' @param PPFD photosynthetic photon flux density, umol m^-2 s^-1
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
##' @param atm_vpd atmospheric vapour pressure deficit, kPa
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

# Validate and default the drivers, without applying them.
#
# WHY THIS IS SPLIT OUT. `leaf_gradient()` re-drives the leaf once per
# perturbation -- eleven times for a four-parameter gradient -- and all of this
# validation and defaulting produces the same answer every time bar the one
# parameter being moved. Resolving once and applying many times is worth ~12% of a
# gradient (see PLAN), but only if there is ONE definition of the rules: the
# defaults here are load-bearing (1 m layers, the nominal networks, the
# single-path placeholder depth) and a second copy in the gradient code would be
# free to drift from this one silently. So the gradient calls this, not a
# reimplementation of it.
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

  # The C++ boundary rejects this too, and it is the enforcement. This check is
  # here only to name the fix, because "psi_soil must be positive magnitudes"
  # arriving from inside a solve is a less useful thing to read than this is.
  if (any(psi_soil < 0)) {
    stop("`psi_soil` must be positive magnitudes in MPa, not signed ",
         "potentials. If this came from an older script, drop the minus sign ",
         "rather than negating the result.", call. = FALSE)
  }

  single <- identical(x$supply_kind, "single")

  if (single) {
    # The single-potential path reads psi_soil[1] and nothing else: no depth
    # profile, no root carbon. Rather than quietly ignoring extra layers -- which
    # would let someone pass a profile and believe it was used -- say so.
    if (n != 1L) {
      stop("`psi_soil` must be a single value on the single-potential supply ",
           "path; got ", n, " layers. Build the leaf with ",
           "leaf_supply_multilayer() for a layered profile.", call. = FALSE)
    }
    # `soil_depth` really is multi-layer-only: nothing on this path reads a depth
    # profile, and silently accepting one someone took the trouble to pass is how
    # a plausible wrong number happens.
    if (!is.null(soil_depth)) {
      stop("`soil_depth` is not used on the single-potential supply path -- ",
           "there is no depth profile to read. Drop it, or build the leaf with ",
           "leaf_supply_multilayer().", call. = FALSE)
    }
    # `root_network` IS used here, and is the same argument the multi-layer path
    # takes: one series resistance in r_R_V_sum. That is the whole point of the
    # consistency change -- the resistance is a driver on both paths now, where it
    # used to be a leaf_supply_single() argument on this one.
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

  # In `set_physiology`'s own argument order, so a caller can splat it positionally
  # without naming ten arguments at a call site that runs in a loop.
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
##' Reads the outputs off a leaf that has been solved. Call after
##' `$find_root_collar_psi()`; before that the values are the missing-value
##' sentinels the object was constructed with.
##'
##' Costs about 4 µs, so it is usable in a loop. It was 180 µs until it stopped
##' reading the twelve outputs through twelve separate calls into C++ and stopped
##' building its one row with `data.frame()` (#39) -- 45× more than the ~3 µs
##' solve it was reporting on.
##'
##' @param x a solved `Leaf`
##' @return A one-row data.frame.
##' @seealso [leaf_solve()], which does the whole thing in one call.
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$find_root_collar_psi()
##' operating_point(l)
##' @export
operating_point <- function(x) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  # One row, from the same C++ reader leaf_solve()'s loop uses, rather than a
  # data.frame() call listing the columns a second time: the names are written
  # down once, in .operating_point_names.
  #
  # ⚠️ Built directly rather than through data.frame(), which costs 158 us
  # against 2 us for this -- on a function called once per solved point, to
  # report a 3 us solve. What data.frame() spends it on is checking and recycling
  # twelve arguments that are already twelve length-1 doubles by construction.
  # The result is `identical()` to what data.frame() returned, which
  # test-surface.R asserts rather than assumes; note row.names has to be the
  # integer 1L and not 1.0, or it would not be.
  v <- as.list(x$operating_point_values())
  names(v) <- .operating_point_names
  structure(v, class = "data.frame", row.names = 1L)
}

# What an operating point IS: the names of the twelve values
# `Leaf::operating_point_values()` returns, in its order. `operating_point()`
# wraps them in a one-row data.frame and `leaf_solve()` fills a matrix row with
# them, so the two cannot disagree about which outputs there are.
#
# ⚠️ THE ORDER IS AN INTERFACE, and it is one nothing in the types enforces --
# the C++ side returns a flat vector, because a flat vector is what crosses the
# boundary for free. test-surface.R asserts these names line up with the
# twelve individual bindings by reading each one and comparing, so a field
# inserted on either side without the other fails there instead of silently
# shifting a column. test-golden.R then compares leaf_solve()'s output
# bit-exactly against a file generated in C++.
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
  "lambda",         # dA/dE
  "g1_eff"          # the Medlyn g1 this leaf implies
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
##' ~204 µs, some 70 solves, and only ~32 µs of that is the two vulnerability
##' splines; the rest is R-side object construction. See [set_traits()] for
##' varying traits without reconstructing, and note it costs 21.8 µs rather than
##' 0.02 µs when the trait you change owns a vulnerability spline.
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
##'   (the default) or [leaf_supply_single()]. On the single-potential path
##'   `soil_depth` and `root_network` must be omitted, and each `psi_soil` is one
##'   value rather than a profile.
##' @param reuse solve every row with one `Leaf` object (`TRUE`) or construct a
##'   fresh one per row (`FALSE`, the default). Reuse is faster because the two
##'   vulnerability splines are built once, and it is safe -- every exit from the
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
                       reuse = TRUE) {
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
  # One network for every row, or one per row. The inherits() test comes FIRST
  # because a RootNetwork *is* a list, so length() on it would be 5 and it would
  # be read as five rows' worth of networks.
  networks <- if (is.null(root_network)) {
    NULL
  } else if (inherits(root_network, "RootNetwork")) {
    rep(list(root_network), n)
  } else {
    .recycle_to(root_network, n, "root_network")
  }

  shared <- if (reuse) leaf_model(traits, control, supply) else NULL

  # ⚠️ COLUMNWISE, AND THAT IS MOST OF THE PERFORMANCE STORY OF THIS FUNCTION.
  # It used to build a one-row data.frame of drivers per row, cbind an
  # `operating_point()` data.frame to it, and rbind the lot at the end. Measured
  # over 32 rows on the single-potential path, that cost **344 us per row**
  # against ~5 us of solving, and 310 us of it was the two data.frame() calls,
  # the cbind and the rbind -- none of which is the model (#39). Filling
  # preallocated storage and assembling once takes it to 34 us; the one-call
  # reader below takes it to **20.5 us**, which is 1.2x the cost of driving the
  # object by hand rather than the 22x it was. Same result to the bit.
  #
  # So: nothing that allocates per row. The driver columns are already length-n
  # vectors from the recycling above, the outputs go into one matrix, and there
  # is exactly one data.frame() call, at the end. Its ~105 us is a per-CALL cost
  # rather than a per-row one, which is why it is still data.frame(): at n = 1 it
  # is invisible beside the 204 us of constructing the Leaf.
  outputs <- matrix(NA_real_, nrow = n, ncol = length(.operating_point_names),
                    dimnames = list(NULL, .operating_point_names))

  for (i in seq_len(n)) {
    l <- if (reuse) shared else leaf_model(traits, control, supply)
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
    l$find_root_collar_psi()
    outputs[i, ] <- l$operating_point_values()
  }

  # The drivers, then the outputs. psi_soil is reported as the wettest layer
  # when there are several, with the layer count alongside, so the column stays
  # numeric and the frame stays rectangular.
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
