# Trait gradients at a solved operating point: set_traits() and leaf_gradient().
#
# The reference implementation: inst/include/phylloptim/gradient.hpp reproduces
# this file's composite bit-for-bit, and test-gradient-batch.R asserts it. So the
# arithmetic order here and the order of calls into the leaf are load-bearing: do
# not reassociate or reorder either without changing gradient.hpp to match.
#
# When the exact gradient beats differencing is ?leaf_gradient's "What it
# costs"; the measured per-call costs are README's "Performance from R".

##' Replace the traits on an existing leaf
##'
##' Changes the traits of a `Leaf` in place, leaving [leaf_control()]'s numerical
##' settings alone. This is the fast path for anything that varies traits --
##' calibration, sensitivity analysis, a trait response curve -- because
##' constructing a `Leaf` from R costs about 45 times as much as solving one.
##'
##' @section You must set the drivers again afterwards:
##' This returns the leaf to its just-constructed state, so [set_drivers()] has to
##' be called before the next solve. That is not caution: `vcmax_`, `jmax_` and
##' `R_d_` are derived from the traits inside `set_physiology()`, so they are
##' genuinely unknown until the drivers are re-supplied.
##'
##' The reason it is a function rather than sixteen assignable fields is that
##' `leaf$vcmax_25 <- x` could not be made correct. Three separate pieces of
##' derived state go stale on a bare trait write: the two pre-integrated
##' vulnerability tables, the solved operating point, and -- least visibly --
##' `vcmax_`/`jmax_`/`R_d_`, which `set_physiology()` computes behind a cache. The
##' cache is keyed on the traits it reads, so setting the drivers again refreshes
##' it, but nothing short of `set_traits()` refreshes the tables or the solved point.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param traits a [leaf_traits()] object
##'
##' @return `x`, invisibly.
##' @seealso [leaf_model()], [set_drivers()], [leaf_gradient()]
##' @examples
##' l <- leaf_model()
##' set_traits(l, leaf_traits(vcmax_25 = 120))
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$optimise()
##' operating_point(l)
##' @export
set_traits <- function(x, traits) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  # Positional, in the C++ setter's argument order; a trait added there is an
  # argument added here.
  x$set_traits(traits$vcmax_25, traits$stem_c, traits$stem_P50,
               traits$root_c, traits$root_P50, traits$TF24_beta2,
               traits$jmax_25, traits$phi_J, traits$curv_fact_elec_trans,
               traits$curv_fact_colim, traits$TF24_cost_scale, traits$R_d_25,
               traits$JS22_gamma, traits$CMax_a, traits$CMax_b, traits$g_min)
  invisible(x)
}

##' Trait gradients of a solved operating point
##'
##' The derivatives of the solved outputs with respect to the traits:
##' `dA/dtheta`, `dgc/dtheta`, `dpsi_stem/dtheta`, `dcollar/dtheta` and
##' `dprofit/dtheta`, at one operating point.
##'
##' @section Choosing which model to differentiate:
##' `model` selects it, and every model in the package has a route. `"collar"`,
##' the default, is the TF24 cost maximised over the root-collar potential -- the
##' production path. Naming any of the seven cost curves
##' (`cost_curve_names()`) differentiates the corresponding **stem** optimum
##' instead: `psi_stem` with the upstream potential pinned at `psi_soil`, which
##' needs the single-potential supply path.
##'
##' The product objectives, `SOX` and `JW26`, need no special case. They maximise
##' `A * g(psi)`, which is the `Log` benefit link applied to the same
##' difference form everything else uses: `d/dpsi[h(A) - C] =
##' h'(A) * dA/dpsi - dC/dpsi` with `h' = 1/A`. One derivative serves all seven.
##'
##' ⚠️ **`"collar"` and `"TF24"` are different models, not two ways of asking one
##' question**, even though both use the TF24 cost. They optimise different
##' variables over different supply topologies, disagree on 20 of 30 driver rows,
##' and the gap does not close as the root resistance goes to zero.
##'
##' ⚠️ **`model` SEATS the leaf, including one you passed as `x`.** The route is
##' applied with `$set_model()` and is left there when this returns, alongside the
##' base-point traits and solve. A loop that alternates models on one `x` is fine;
##' a loop that takes a gradient and then reads `$optimise()` off the same object
##' expecting the default is not.
##'
##' ⚠️ **On a stem route the decision variable is `psi_stem`, so `collar` comes
##' back non-finite.** The collar is not solved for at all there, and an absent
##' number is the truthful answer rather than a gap -- reading one would be reading
##' state from whichever solve ran last.
##'
##' **Every route reaches the composite under `method = "auto"`, and every route
##' uses the same finite-difference step (1e-06).** Both were once untrue.
##' `ProfitMax` normalises by `|A|max`, and while that came from a 500-point scan
##' its argmax was piecewise constant in the parameters -- so the composite had to
##' pin it and returned a PARTIAL (0.0766 against a true 0.0579), and differencing
##' needed a step large enough to cross a scan cell. `|A|max` is found by a
##' root-find now, so it is re-solved per perturbation, the composite returns the
##' TOTAL, and no model carries that split.
##'
##' ⚠️ **A zero column for another curve's parameter means two different things,
##' and they are not distinguishable from the number.** On a `collar` or `TF24`
##' route, `JS22_gamma`, `CMax_a`, `CMax_b` and the rest come back exactly zero,
##' and read one way that is correct: those parameters do not enter the TF24
##' cost, so their derivative is genuinely zero. Read the other way it is a trap:
##' if you are *fitting* one of those curves, the zero is telling you this
##' differentiated a different model, not that the parameter does nothing.
##' Measured on the single-potential path at `psi_soil` 1.5, the TF24 collar
##' optimum sits at 2.09 MPa and the JS22 optimum at 3.15.
##'
##' An optimiser handed such a column concludes the parameter is inert and leaves
##' it at its starting value, which looks like convergence. **Name the curve you
##' are fitting in `model`**; `.gradient_available_pars()` reports which
##' parameters a given route can move.
##'
##' @section The maths, and why it is not just a finite difference:
##' The reported outputs are evaluated at the profit-maximising collar potential
##' `psi*`, so a trait moves them two ways: directly, and by moving `psi*`. At an
##' interior optimum `dprofit/dpsi = 0`, and differentiating that condition gives
##' `psi*`'s response without re-solving the model:
##'
##' \deqn{d\psi^*/d\theta = -M / H, \quad
##'       M = \partial^2 profit/\partial\psi\partial\theta, \quad
##'       H = \partial^2 profit/\partial\psi^2}
##'
##' \deqn{dY/d\theta = \partial Y/\partial\theta|_\psi +
##'                    (\partial Y/\partial\psi)(d\psi^*/d\theta)}
##'
##' Both second derivatives are differences of `dprofit_droot_collar_psi()`, which
##' is exact in `psi` and smooth in the traits -- `psi` is one of its arguments, so
##' no argmax is involved and `M` is a first difference of an exact quantity. It is
##' stable to seven significant figures across five decades of step size.
##'
##' The second term is not a correction. For `TF24_cost_scale`, `TF24_beta2` and
##' `stem_c` it is 100% of the answer, and for `vcmax_25` 52%.
##'
##' @section profit, which is the one output the envelope theorem reaches:
##' The first four outputs are what a gas-exchange calibration observes. `profit`
##' is here because it is what a **demographic** consumer bills: plant's carbon
##' budget reads `leaf.profit_`, not `assim_colimited_`, so before it was added
##' the two sets were disjoint and no gradient from this package reached a
##' demographic model at all.
##'
##' It is also the cheapest of the five, because it is the objective. At an
##' interior optimum `dprofit/dpsi = 0`, so the indirect term above vanishes
##' identically and
##'
##' \deqn{dprofit/d\theta = \partial profit/\partial\theta|_{\psi^*}}
##'
##' -- the direct partial alone, with no `dY/dpsi` and no `-M/H`. That is the
##' envelope theorem, and it is the only place this package uses it.
##'
##' ⚠️ **The term this drops is noise, not an `h^2` truncation, and that is why
##' dropping it matters.** `profit` is the maximum, so it is flat, and a central
##' difference of it at `psi*` divides the solve's ~1e-09 floor by a ~1e-06 step.
##' Over the golden grid's 136 interior rows:
##'
##' | | median | max |
##' |---|---|---|
##' | `\|dprofit/dpsi\|`, exact (forward AD) | 4.8e-15 | 5.4e-10 |
##' | `\|dprofit/dpsi\|`, central difference | 7.8e-10 | 2.1e-04 |
##' | relative move in `dprofit/dtheta` if kept | 2.7e-10 | 8.0e-05 |
##'
##' Eleven orders between the two instruments at the median, and the worst row
##' sits in the band this package calls a real difference. So `dprofit/dtheta` is
##' set from the direct term rather than computed through the composite, exactly
##' as `collar` is set from `dpsi*/dtheta` rather than differenced.
##'
##' ⚠️ **The envelope does not survive the active set.** At a pinned optimum
##' `psi*` is a trait-dependent bound and `dprofit/dpsi` is not zero there, so the
##' indirect term is real. This applies the identity only where `status` is
##' `"interior"`; everywhere else `profit` goes through the same finite-difference
##' fallback as the other four.
##'
##' @section A collar potential you supply, instead of the one this solves for:
##' Pass `psi` and the outputs are evaluated at that collar rather than at the
##' argmax. This is for a caller whose model **tracks** the optimum instead of
##' finding it — plant's TF24f carries the collar potential as an ODE state,
##' `dpsi/dt = k * dprofit/dpsi`, so at finite gain its `dprofit/dpsi` is
##' deliberately non-zero and *is* the acclimation rate.
##'
##' There the derivation above simplifies rather than breaks. `psi` is exogenous,
##' so the indirect term is whatever the caller says it is:
##'
##' \deqn{dY/d\theta = \partial Y/\partial\theta|_\psi +
##'                    (\partial Y/\partial\psi)(d\psi/d\theta)}
##'
##' with `dpsi_dtheta` supplied — defaulting to zero, the partial at fixed
##' collar. Nothing is derived from `-M/H`, so nothing needs stationarity, and
##' `method` is refused because the two routes it chooses between are both about
##' a solved optimum.
##'
##' `M`, `H` and `dY_dpsi` come back in the result, because a caller with a
##' *dynamic* `psi` cannot supply `dpsi_dtheta` as a constant: for the
##' gradient-ascent law above it obeys `ds/dt = k(M + H s)`, and those are its
##' coefficients. Its fixed point is `-M/H`, which is what the solving path
##' returns — the two agree in the large-gain limit.
##'
##' ⚠️ **`stationarity` is still reported, and is used for exactly one thing.**
##' It no longer routes anything; it now measures how far the collar you supplied
##' sits from the optimum. The one decision it makes is `profit`'s: at a
##' stationary point the envelope theorem applies and the analytic zero is used,
##' and away from one the exact `dprofit/dpsi` is used instead of differencing
##' it. So `psi = <the solved psi*>` with `dpsi_dtheta = -M/H` reproduces the
##' solving path **bit-for-bit**, which the tests assert.
##'
##' ⚠️ **An infeasible `psi` gets no gradient either**, reported as
##' `"no-gradient"`. `dprofit/dpsi` is a sentinel `0.0` on the shut-down and
##' reversed-gradient exits rather than a derivative, and this path never divides
##' by `H`, so it would otherwise adopt that zero as the real thing.
##'
##' ⚠️ **A `psi` outside the feasible interval is refused, and gets no
##' gradient.** The leaf does not project it onto the interval's end, which
##' would make the outputs move with the *bound* rather than with `dpsi_dtheta`.
##' `status` reports `"refused"`, `psi` in the result is the one you gave, and
##' everything else is `NA`. This is reported rather than thrown because a fit
##' reaches these points routinely: a trait perturbation moves the interval past
##' a fixed `psi`. A `psi` inside the interval but within one step of an end is
##' `"refused"` too, with its value, because `dY/dpsi` cannot be centred there.
##'
##' @section The active set, which is the reason this function is careful:
##' Stationarity is the premise of the whole derivation, and it fails when the
##' optimum is pinned to an end of the feasible collar interval. There `psi*` is a
##' bound, `dprofit` is not zero at the answer, and `-M/H` is not the bound's
##' derivative. The formula does not fail loudly: measured at a boundary-soil point
##' the true gradient is about `1e-08` and the bare composite returns `O(1)` --
##' plausible-looking and wrong by seven orders of magnitude. About one operating
##' point in six is pinned across this package's test grid, all at soil suctions of
##' 3-4 MPa, which is the dry end a drought calibration will visit.
##'
##' So the premise is **tested rather than assumed**. The test is the implied
##' Newton step `|dprofit(psi*) / H|`, a distance in MPa: at a stationary point it
##' is zero to solver precision, and at a pinned one it is not. Over the package's
##' 288-point grid the two populations are five orders of magnitude apart -- worst
##' interior `5e-11`, smallest pinned `6e-06` -- so `stationarity_tol` sits in an
##' empty band rather than on a judgement call.
##'
##' When the premise does not hold, this falls back to a central finite difference
##' of the whole solve, which respects the constraint by construction because it
##' differences the constrained answer. That is the one thing the crude method does
##' better, and `method` in the result says which was used.
##'
##' It is better in a second way, which is sharper than accuracy. The dry end of
##' the feasible collar interval -- `psi_crit` for the curves it bounds, the end of
##' the stem curve for the rest -- does not appear in the profit function at all,
##' so the composite necessarily returns **exactly zero** for it. At an interior
##' optimum that is correct. At a dry pin it *is* the binding constraint, and a
##' zero is arguably a worse answer than a wildly wrong one: it tells an optimiser
##' the parameter does nothing, and nothing about it looks suspicious. The
##' difference of the solve gets it for free.
##'
##' There is a second guard behind that one, found by measurement rather than
##' designed: a pinned optimum sits so close to its bound that the difference in
##' `psi` cannot be centred there without stepping out of the interval, and that is
##' detected too.
##' Over this package's grid it catches every pinned point on its own, so forcing
##' `method = "ift"` at one fails loudly instead of returning the wrong number.
##' Do not read that as making the stationarity test redundant -- it fires only
##' when the feasible interval is narrow, which is true of these points but not of
##' pinned optima in general.
##'
##' @section Two of these are not traits:
##' `pars` also accepts `leaf_specific_conductance_max` and, on the
##' single-potential path, `resistance`. They are here because a calibration fits
##' them: of a companion calibration study's four free parameters two are traits and two are
##' these, so restricting `pars` to [leaf_traits()] left half of that fit without
##' an exact gradient.
##'
##' Nothing in the derivation cares that the parameter is a trait — the implicit
##' function theorem is applied to `dprofit/dpsi = 0`, and any parameter the profit
##' function depends on goes through it identically. What differs is only which
##' setter applies the perturbation. `leaf_specific_conductance_max` is also how
##' plant's height reaches the leaf, so the same route gives a gradient with
##' respect to plant state and not only with respect to traits.
##'
##' ⚠️ **Both move the feasible collar interval, not just the profit function** —
##' `resistance` directly and `leaf_specific_conductance_max` through the supply
##' curve — so a perturbation can push a bound past `psi*` more readily than a
##' trait perturbation can. That is the case the refusal detector below exists for,
##' and it fails loudly rather than quietly.
##'
##' @section What it costs, and when that beats differencing:
##' Both routes are linear in a parameter count, and it is a **different** count
##' for each. With `N` observations:
##'
##' \deqn{T_{FD} = 2\,P_{fit}\,T_{objective}}
##' \deqn{T_{exact} = N\,(t_{construct} + t_{solve}) + 2\,P_{model}\,N\,t_{eval}}
##'
##' An optimiser differences the `P_fit` parameters it is moving. This function is
##' asked for the `P_model` parameters the leaf has — the length of `pars`. Those
##' are equal only when the fit varies traits directly; a parameterisation in
##' between (pooling, a hierarchy, any shared or derived parameter) makes
##' `P_fit > P_model`, and then this route's cost does not follow the optimiser's
##' dimension.
##'
##' ⚠️ **So "is the exact gradient faster" has no answer independent of your
##' parameterisation, and both regimes are measured.** `vignette("fitting")` fits
##' 72 observations with `P_fit == P_model` and this loses at every count up to 13.
##' A companion calibration study has `N = 1327`, `P_fit = 40`,
##' `P_model = 4`, wins 3.4x, and reaches the same optimum — with its
##' 57-parameter variant costing the same as its 40-parameter one, because
##' `P_model` is 4 in both.
##'
##' ⚠️ **Always pass `pars`.** It is `P_model`, so the default — every trait plus
##' `leaf_specific_conductance_max` on the multi-layer path — is the most expensive
##' thing you can ask for, and a fit that reads four of them pays for the rest.
##'
##' The intercept is a fresh `Leaf` per call. Pass `x` to reuse one across calls
##' (#52); without it, construction is 29% of the exact gradient in the study
##' above.
##'
##' @section Precision:
##' Do not ask for more than about `1e-09` from any of this. That is the
##' achievable precision of the solved outputs themselves, set by the tolerance of
##' the intercellular-CO2 root-find.
##'
##' @section Reusing a leaf across gradients:
##' Constructing a `Leaf` is the largest single cost on this surface, and this
##' function does it once per call: it is **most of a one-parameter gradient** on a
##' fresh leaf, paid once per observation by a fit that differentiates per
##' observation. The measured ratios are in the README's "Performance from R", the
##' one place the package quotes them. Pass `x` to reuse one:
##'
##' ```r
##' l <- leaf_model(traits, control, supply)
##' for (i in seq_len(n)) leaf_gradient(psi_soil = ps[i], x = l, traits = traits, pars = FIT)
##' ```
##'
##' `x` is a **vessel, not the point**: `traits` still says where the gradient is
##' taken and is applied to `x`, rather than read from it (a `Leaf` does not expose
##' its traits). That is why `traits` is required with `x` — defaulting it would
##' differentiate at the package defaults on somebody else's leaf and return a
##' plausible answer for a point nobody asked about.
##'
##' `control` and `supply` are fixed when `x` is built, so passing them alongside it
##' is an error rather than silently resolved.
##'
##' **`x` is left solved at the base point**, not at the last perturbation — the
##' state this function seats it in, restored on the way out including on an error.
##' The restore costs one re-trait and one solve, a small fraction of the
##' construction it saves.
##'
##' @inheritParams set_drivers
##' @param x an existing `Leaf` from [leaf_model()] to reuse instead of building
##'   one, for a loop that takes many gradients. `traits` is required with it, and
##'   `control`/`supply` must be omitted. See the section above.
##' @param traits a [leaf_traits()] object: the point in trait space to
##'   differentiate at
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_singlelayer()]
##' @param model which optimality model to differentiate: `"collar"` (the
##'   default, the production path) or any name from [cost_curve_names()], which
##'   differentiates that curve's stem optimum instead. See *Choosing which model
##'   to differentiate*.
##' @param pars what to differentiate with respect to. Any of the sixteen
##'   [leaf_traits()] names, plus `"leaf_specific_conductance_max"`, on the
##'   single-potential path only `"resistance"`, on `model = "CF77"` only
##'   `"CF77_lambda_"`, on `model = "TF24_floor"` only `"TF24_floor_lambda_o"`,
##'   and on `model = "LeastCost"` only `"LeastCost_beta"`.
##'   On a reused leaf `x` with `use_thermal_damage_` on, also the layer's inputs
##'   that reach the solve: its states `"A_"`, `"I_"`, `"D_"`, the optimum shift
##'   `"topt_offset_"`, `"dTopt_max_"`, `"K_A_"`, and every quantity a maintenance
##'   cost multiplies (`"k_i_"`, `"k_i_baseline_"`, `"k_rec_"`, `"k_rec_Ea_"`,
##'   `"c_acclim_maint_"`, `"c_repair_maint_"`,
##'   `"c_protect_maint_"`, `"c_repair_flux_"`). The layer's dynamics parameters
##'   move no solved output and are refused; their gradient is
##'   [leaf_damage_gradient()]'s. Defaults to all the route allows.
##'
##'   ⚠️ At 25 C the peaked Arrhenius is normalised to its reference, so
##'   `topt_offset_` and `dTopt_max_` have an exactly zero effect there -- a true
##'   zero, not a missing row.
##'   The protection cost is signed, `k_i_baseline_ - k_i_`, so its columns are
##'   the same linear slope on both sides of the default, where the two are
##'   equal.
##'
##'   `dY/dvcmax_25` is a PARTIAL at fixed respiration: `R_d_25` is its own trait,
##'   so a fit that wants respiration to follow Vcmax moves both and adds the two
##'   columns.
##'
##'   `beta_R_H` and `beta_R_V` were here until #33 and are not any more: they
##'   parameterise the root-architecture model, which the leaf no longer runs, so
##'   this function has no way to reach them. Getting a gradient in either means
##'   solving at two networks yourself. The perturbation itself is cheap —
##'   [root_network_from_carbon()] is homogeneous of degree 1 in each constant, so
##'   scaling `beta_R_H` scales `r_R_H_min` by the same factor and needs no
##'   rebuild — but the two solves are still two solves.
##' @param psi a collar water potential to evaluate at, in MPa as a positive
##'   magnitude, instead of solving for the profit-maximising one. `NULL` (the
##'   default) solves. See the section above; `method` cannot be given with it.
##' @param dpsi_dtheta how the prescribed `psi` itself responds to each
##'   parameter, as one value per entry of `pars` (or one value recycled, or a
##'   vector named by `pars`). Defaults to zero, which is the partial derivative
##'   at fixed collar. Must be finite — an infinite value would return an
##'   all-infinite gradient row with a finite, plausible-looking `profit` in it,
##'   because the envelope theorem assigns that column rather than multiplying
##'   through. Only meaningful with `psi`.
##' @param step relative step for the trait difference. The default `1e-06` is
##'   near the middle of the five decades over which the mixed partial was
##'   measured stable; it is also used, relative to the collar potential, for the
##'   two differences in `psi`.
##' @param stationarity_tol how far from stationary, as an implied Newton step in
##'   MPa, the solved point may be before the optimum is called pinned and the
##'   finite-difference fallback is used. The default `1e-08` sits in the empty
##'   band between the two measured populations; there is no operating point in
##'   this package's grid within four orders of magnitude of it.
##' @param fast_stem_curve perturb `stem_P50` through `perturb_stem_P50()`, which
##'   moves `stem_b` and leaves the seated physiology alone (`TRUE`, the
##'   default), rather than through `set_traits()`. Neither route rebuilds the
##'   stem's vulnerability table, which depends on `stem_c` alone; `FALSE` exists
##'   so the equivalence is checkable rather than assumed, and the package's tests
##'   require the two routes to agree.
##' @param method which route to take: `"auto"` (the default) uses the
##'   implicit-function composite where its premise holds and the finite
##'   difference where it does not. `"ift"` and `"fd"` force one route, for
##'   checking one against the other. Forcing `"ift"` at a pinned optimum returns
##'   a confidently wrong number -- that is the whole reason `"auto"` exists --
##'   and `status` still reports the truth about the point.
##'
##' @return A list with
##'   \describe{
##'     \item{`gradient`}{a matrix, one row per trait in `pars`, with columns
##'       `A` (umol C m^-2 s^-1 per trait unit), `gc`, `psi_stem`, `collar` and
##'       `profit` (also umol C m^-2 s^-1 per trait unit)}
##'     \item{`value`}{the solved outputs the gradient is taken at}
##'     \item{`method`}{`"ift"` if the implicit-function composite was used,
##'       `"fd"` if the fallback was, `"prescribed"` if `psi` was given}
##'     \item{`status`}{`"interior"`, `"pinned"` or `"no-gradient"` on the
##'       solving path; `"prescribed"`, `"refused"` or `"no-gradient"` with
##'       `psi`}
##'     \item{`H`}{the curvature of profit in the collar potential}
##'     \item{`stationarity`}{the implied Newton step, in MPa: what `method` was
##'       decided on, or with `psi` how far the collar you gave is from the
##'       optimum}
##'     \item{`M`}{the mixed partial `d2profit/dpsi dtheta`, one per entry of
##'       `pars`, or `NA` where the fallback ran}
##'     \item{`dY_dpsi`}{`dY/dpsi` at fixed traits, one per output. `profit`'s
##'       entry is the analytic zero at a stationary point and the exact
##'       `dprofit/dpsi` elsewhere -- see the envelope section}
##'     \item{`psi`}{the collar the outputs were evaluated at: `psi*` on the
##'       solving path, and the `psi` argument on the prescribed path, refused
##'       or not}
##'   }
##'
##' @seealso [leaf_solve()] for the operating point itself, [set_traits()].
##' @examples
##' g <- leaf_gradient(psi_soil = 2.0, PPFD = 900)
##' g$method
##' g$gradient[c("vcmax_25", "stem_P50"), ]
##'
##' # In dry soil at high VPD the optimum pins against the edge of the feasible
##' # range, the composite's premise fails, and the fallback is used instead.
##' dry <- leaf_gradient(psi_soil = 4.5, PPFD = 900, atm_vpd = 3.0)
##' dry$method       # "fd"
##' dry$status       # "pinned"
##' dry$stationarity # how far from stationary the optimum is, in MPa
##' @export
leaf_gradient <- function(psi_soil,
                          PPFD = 900,
                          soil_depth = NULL,
                          root_network = NULL,
                          leaf_specific_conductance_max = 3.14e-5,
                          atm_vpd = 2.0,
                          ca = 40.0,
                          leaf_temp = 25.0,
                          atm_o2_kpa = 21.0,
                          atm_kpa = 101.3,
                          x = NULL,
                          traits = leaf_traits(),
                          control = leaf_control(),
                          supply = leaf_supply_multilayer(),
                          pars = NULL,
                          psi = NULL,
                          dpsi_dtheta = NULL,
                          step = 1e-6,
                          stationarity_tol = 1e-8,
                          method = c("auto", "ift", "fd"),
                          fast_stem_curve = TRUE,
                          model = "collar") {
  method <- match.arg(method)
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  if (!(is.numeric(step) && length(step) == 1L && step > 0)) {
    stop("`step` must be a single positive number", call. = FALSE)
  }
  # A caller who names `step` means it; otherwise the finite-difference arm picks
  # one matched to the solver whose argmax it is differencing.
  if (!missing(step)) {
    attr(step, "user_set") <- TRUE
  }
  # Use the return: it coerces an integer `psi` to the double C++ receives, so the
  # R reference and the C++ batch see one value.
  psi <- .gradient_check_psi(psi, dpsi_dtheta, method)

  # --- reusing a Leaf --------------------------------------------------------
  # `x` is a vessel: `traits` says where the gradient is taken and is applied to
  # it, so `traits` is required. `control` and `supply` are fixed in `x` and are
  # refused; `supply` is rebuilt from it because only its `kind` is read.
  if (!is.null(x)) {
    if (!inherits(x, "Leaf")) {
      stop("`x` must be a Leaf, from leaf_model(); got ", class(x)[[1]],
           call. = FALSE)
    }
    if (missing(traits)) {
      stop("`traits` must be given with `x`: it is the point the gradient is ",
           "taken at, and it is applied to `x` rather than read from it. Pass ",
           "the same leaf_traits() the leaf was built with.", call. = FALSE)
    }
    if (!missing(control)) {
      stop("`control` comes from `x` and cannot be given alongside it -- the ",
           "tolerances and the integrator were fixed when `x` was built.",
           call. = FALSE)
    }
    if (!missing(supply)) {
      stop("`supply` comes from `x` and cannot be given alongside it. `x` is on ",
           "the ", x$supply_kind, " path.", call. = FALSE)
    }
    supply <- if (identical(x$supply_kind, "single")) {
      leaf_supply_singlelayer()
    } else {
      leaf_supply_multilayer()
    }
  }
  # The thermal-damage layer's inputs are parameters ONLY on a leaf with its gate
  # on, which only a reused leaf can have: one built here has it off.
  damage <- !is.null(x) && isTRUE(x$use_thermal_damage_)

  # Resolved here rather than in set_drivers() because `resistance` can be
  # differentiated, so theta must hold the value the solve uses. The same
  # function set_drivers() calls, so there is one default.
  if (identical(supply$kind, "single") && is.null(root_network)) {
    root_network <- .default_series_resistance()
  }

  drivers <- list(psi_soil = psi_soil, PPFD = PPFD, soil_depth = soil_depth,
                  root_network = root_network,
                  leaf_specific_conductance_max = leaf_specific_conductance_max,
                  atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp,
                  atm_o2_kpa = atm_o2_kpa, atm_kpa = atm_kpa)

  # Model-owned fields are read off the leaf when there is one, so theta holds
  # the value the solve uses; NA on a route that does not own the field, so the
  # setter writes nothing (see .owned_value()).
  theta <- .gradient_theta(traits, leaf_specific_conductance_max, supply,
                           root_network,
                           CF77_lambda = .owned_value(model, "CF77_lambda_", x),
                           TF24_floor_lambda_o =
                             .owned_value(model, "TF24_floor_lambda_o", x),
                           LeastCost_beta =
                             .owned_value(model, "LeastCost_beta", x),
                           damage = if (damage) .gradient_damage_values(x))
  if (is.null(pars)) {
    # What this route can differentiate, not every slot in theta: theta holds
    # other routes' owned slots, which the check below refuses.
    pars <- .gradient_available_pars(identical(supply$kind, "single"), model,
                                     damage)
  }
  # Shared with leaf_gradient_batch(), so the two entry points cannot disagree
  # about which parameters exist or explain a rejection differently.
  .gradient_check_pars(pars, identical(supply$kind, "single"), model, damage)
  # Checked here, unconditionally: the composite reads it only at a usable point,
  # so checking it there would accept a bad value at some psi and not others.
  dpsi_dtheta <- .gradient_dpsi_dtheta(dpsi_dtheta, pars, !is.null(psi))

  # One leaf for the whole gradient (and with `x`, across gradients), re-traited
  # rather than rebuilt: construction is the largest cost on this surface.
  l <- if (is.null(x)) leaf_model(traits, control, supply) else x
  reset <- .gradient_setter(l, traits, drivers, supply, fast_stem_curve, damage)

  if (is.null(x)) {
    do.call(set_drivers, c(list(l), drivers))
  } else {
    # Seat a reused leaf with reset(), not set_drivers() alone: it carries its
    # last user's traits, and set_drivers() never touches traits.
    reset(theta)
    # Leave it solved at the base point, not at the last perturbation (hazard 8),
    # on the model the route seated; `on.exit` because paths below stop().
    on.exit({
      try({
        reset(theta)
        l$optimise()
      }, silent = TRUE)
    }, add = TRUE)
  }
  # --- seat the operating point ---------------------------------------------
  # Solved for by default, imposed by `psi`. Either way `psi_star` is the collar
  # the outputs are evaluated at; it is the argmax only on the solving path.
  route <- .gradient_route(l, model, supply)
  prescribed <- !is.null(psi)
  if (prescribed) {
    # Refused by the leaf and reported, not thrown: a fit visits these points,
    # and `status` is how a batch flags them.
    if (!.gradient_evaluate_feasible(route, psi)) {
      na_out <- stats::setNames(rep(NA_real_, length(.gradient_output_names())),
                                .gradient_output_names())
      return(list(
        gradient = matrix(NA_real_, length(pars), length(na_out),
                          dimnames = list(pars, names(na_out))),
        value = na_out, method = "prescribed", status = "refused",
        H = NA_real_, stationarity = NA_real_, model = route$model,
        M = stats::setNames(rep(NA_real_, length(pars)), pars),
        dY_dpsi = na_out, psi = psi))
    }
    psi_star <- route$psi_star()
  } else {
    route$solve()
    psi_star <- route$psi_star()
  }
  value <- .gradient_outputs(l)

  # Curvature and residual at the seated collar. `route$dprofit` takes psi as an
  # argument, so neither needs the point to be an optimum.
  h_psi <- max(abs(psi_star), 1) * step
  # Read with its feasibility: dprofit is a 0.0 sentinel on the shut-down and
  # reversed-gradient exits, indistinguishable by value from stationarity. The
  # prescribed path never divides by `H`, so it relies on the flag.
  checked <- route$checked(psi_star)
  resid <- checked[[1L]]
  feasible <- checked[[2L]] == 1
  H <- (route$dprofit(psi_star + h_psi) -
        route$dprofit(psi_star - h_psi)) / (2 * h_psi)
  # H == 0 (the shut-down sentinel) makes the ratio 0/0, and H > 0 is not a
  # maximum: either way the composite has nothing to stand on.
  usable <- is.finite(H) && H < 0 && is.finite(resid)
  stationarity <- if (usable) abs(resid / H) else Inf

  if (prescribed) {
    # Stationarity does not route here: there is no argmax to be pinned against,
    # so `method` is refused upstream. It is still reported, as the distance from
    # the optimum in MPa, and decides `envelope` below. Only infeasibility
    # disables the point: a prescribed psi may sit where profit is convex, which
    # is fine because nothing here divides by `H`.
    status <- if (!feasible) "no-gradient" else "prescribed"
  } else {
    status <- if (!usable) "no-gradient" else
      if (stationarity > stationarity_tol) "pinned" else "interior"
  }

  # `status` describes the point and is reported whichever route runs; `use_ift`
  # is the route. They differ only when the caller has forced one.
  auto_ift <- identical(status, "interior")
  use_ift <- if (prescribed) feasible else switch(method,
                    auto = auto_ift,
                    ift = TRUE,
                    fd = FALSE)
  if (use_ift && !prescribed && !usable) {
    stop("leaf_gradient(): method = \"ift\" was asked for at a point with no ",
         "usable curvature (H = ", format(H), "), so -M/H has nothing to stand ",
         "on. This is a shut-down or otherwise determined operating point; use ",
         "method = \"auto\".", call. = FALSE)
  }
  # Refused on the stationarity test: the step check below does not catch every
  # pinned optimum, and at one the composite is wrong by orders of magnitude.
  if (use_ift && !prescribed && identical(status, "pinned")) {
    stop("leaf_gradient(): method = \"ift\" was asked for at a pinned optimum, ",
         "where psi* is held by a bound and the implicit-function composite ",
         "does not apply. Use method = \"auto\".", call. = FALSE)
  }

  if (use_ift) {
    # dY/dpsi at fixed traits, centred. A stencil reaching past an end of the
    # feasible interval is detected rather than taken one-sided; that catches
    # some pinned optima, but the refusal above is what guards them.
    hi <- .gradient_outputs_at(l, psi_star + h_psi, route)
    lo <- .gradient_outputs_at(l, psi_star - h_psi, route)
    if (is.null(hi) || is.null(lo)) {
      if (prescribed) {
        # Within `h_psi` of an end: dY/dpsi cannot be centred, and differencing
        # the solve would answer about the optimum, so the row is refused with
        # its value and no gradient.
        status <- "refused"
        use_ift <- FALSE
      } else if (identical(method, "ift")) {
        stop("leaf_gradient(): method = \"ift\" was asked for at a point whose ",
             "feasible collar interval is narrower than one step, so dY/dpsi ",
             "cannot be centred on psi*. Use method = \"auto\".", call. = FALSE)
      } else {
        use_ift <- FALSE
        status <- "pinned"
      }
    } else {
      dY_dpsi <- (hi - lo) / (2 * h_psi)
      # profit's dY/dpsi is the exact `dprofit` already in hand, not a difference;
      # one rule on both paths keeps `psi = psi*` reproducing the solve exactly.
      dY_dpsi[["profit"]] <- resid
    }
  }

  # The envelope theorem: at a stationary point profit's indirect term is zero,
  # so the composite assigns the direct partial (see ?leaf_gradient). Keyed on
  # `stationarity`, not on `status`, which differs between the two paths at
  # psi = psi* and would break their bit-for-bit match, nor on `use_ift`: at a
  # pinned optimum the indirect term is real.
  envelope <- isTRUE(usable && stationarity <= stationarity_tol)

  if (prescribed && !use_ift) {
    grad <- matrix(NA_real_, length(pars), length(.gradient_output_names()),
                   dimnames = list(pars, .gradient_output_names()))
    M <- stats::setNames(rep(NA_real_, length(pars)), pars)
    dY_dpsi <- stats::setNames(rep(NA_real_, length(.gradient_output_names())),
                               .gradient_output_names())
  } else if (use_ift) {
    # One composite for both paths; only dpsi/dtheta differs, derived by the
    # implicit function theorem when solved and supplied by the caller when
    # imposed.
    fit <- .gradient_ift(l, reset, theta, pars, psi_star, H, dY_dpsi, step,
                         fast_stem_curve, dpsi_dtheta = dpsi_dtheta,
                         envelope = envelope, route = route)
    grad <- fit$gradient
    M <- fit$M
  } else {
    grad <- .gradient_fd(l, reset, theta, pars, step, fast_stem_curve, route)
    M <- stats::setNames(rep(NA_real_, length(pars)), pars)
    dY_dpsi <- stats::setNames(rep(NA_real_, length(.gradient_output_names())),
                               .gradient_output_names())
  }

  list(gradient = grad,
       value = value,
       method = if (prescribed) "prescribed" else if (use_ift) "ift" else "fd",
       status = status,
       H = H,
       stationarity = stationarity,
       model = route$model,
       M = M,
       dY_dpsi = dY_dpsi,
       psi = psi_star)
}

# --- the route: which model's optimum is being differentiated -----------------
#
# Chosen once here rather than branched on at each call site. The collar route
# and a stem route are different models, not one restricted to the other
# (hazard 11). `decision` names the output that IS the decision variable, whose
# gradient is assigned dpsi*/dtheta rather than differenced. On a stem route
# `collar` is not solved for, and its NA gradient is the truthful answer.
.GRADIENT_COLLAR_ROUTE <- "collar"

# Every model this can differentiate: the collar route plus every cost curve's
# stem route. The one list `.gradient_model_pars()` validates against and
# `.gradient_route()` dispatches on. Verify a route against a finite difference
# of its solve by sweeping the step, not at one step (hazard 15).
.gradient_models <- function() {
  c(.GRADIENT_COLLAR_ROUTE, cost_curve_names())
}

.gradient_route <- function(l, model, supply) {
  if (identical(model, .GRADIENT_COLLAR_ROUTE)) {
    # Seated, as the stem routes are below: the closures read the seated curve,
    # and a reused leaf may come seated on another.
    l$set_model("TF24", "collar")
    return(list(
      model = model,
      decision = "collar",
      # A root-find on `dprofit == 0`, stationary to ~1e-15, so a 1e-06 trait step
      # is resolvable: psi* moves far more than the solve's own precision.
      fd_step = 1e-6,
      solve = function() l$find_root_collar_psi(),
      psi_star = function() l$opt_root_psi_,
      checked = function(psi) l$dprofit_droot_collar_psi_checked(psi),
      dprofit = function(psi) l$dprofit_droot_collar_psi(psi),
      evaluate = function(psi) l$evaluate_root_collar_psi(psi)
    ))
  }
  if (!(model %in% cost_curve_names())) {
    stop("unknown model \"", model, "\". Available: ",
         paste(.gradient_models(), collapse = ", "), ".", call. = FALSE)
  }
  if (!identical(supply$kind, "single")) {
    stop("the stem routes optimise `psi_stem` with the upstream potential ",
         "pinned at `psi_soil`, so they need the single-potential supply path; ",
         "got the multi-layer one. Use model = \"collar\" for a root network.",
         call. = FALSE)
  }

  # Seated on construction, not in `solve()`: the closures below read whichever
  # model the leaf is set to, and the prescribed path never calls `solve()`.
  # set_traits() does not clear the seated model.
  l$set_model(model, "stem")

  # ProfitMax's `|A|max` is cleared by every perturbation, so it is re-solved
  # before each read rather than pinned at the base point: that is what makes the
  # composite a total derivative. A no-op on every other curve.
  reseat <- if (identical(model, "ProfitMax")) {
    function() l$prepare_profitmax()
  } else {
    function() invisible(NULL)
  }

  list(
    model = model,
    decision = "psi_stem",
    # The same step as the collar route: the stem routes root-find dJ/dpsi == 0,
    # so psi* is stationary to ~1e-15. The floor moves with the solver; re-sweep
    # before changing this or the solver (hazard 15).
    fd_step = 1e-6,
    solve = function() l$optimise(),
    psi_star = function() l$opt_psi_stem_,
    checked = function(psi) {
      reseat()
      l$dprofit_dpsi_stem_checked(psi)
    },
    dprofit = function(psi) {
      reseat()
      l$dprofit_dpsi_stem_checked(psi)[[1L]]
    },
    evaluate = function(psi) {
      reseat()
      l$evaluate_psi_stem_at(psi)
    }
  )
}

# --- internals ---------------------------------------------------------------

# The differentiated outputs, read from C++ rather than restated here, so adding
# one is one edit in gradient.hpp. `collar` is psi* itself, so dcollar/dtheta is
# dpsi*/dtheta. Cached at first call: calling into the shared library while the
# package is being sourced is not safe to rely on.
.gradient_output_names <- local({
  nms <- NULL
  function() {
    if (is.null(nms)) {
      nms <<- gradient_output_names()
    }
    nms
  }
})

# One C++ reader for every output, not a field read each: every `l$field` is an
# R6 active binding wrapping a `.Call`, and this runs once per perturbation. The
# reader is the accessor the bindings wrap, so the values are the same.
#
# Cached at first call: R collates `R/` alphabetically, so
# `.operating_point_names` (leaf-model.R) does not exist when this is sourced.
.gradient_outputs_idx <- local({
  idx <- NULL
  function() {
    if (is.null(idx)) {
      idx <<- match(.gradient_output_names(), .operating_point_names)
    }
    idx
  }
})

.gradient_outputs <- function(l) {
  v <- l$operating_point_values()[.gradient_outputs_idx()]
  names(v) <- .gradient_output_names()
  v
}

# Outputs with the collar held at `psi` rather than optimised. NULL when the
# leaf refuses `psi` as outside the feasible interval, so that a stencil
# reaching past an end is detected rather than differenced one-sidedly.
.gradient_outputs_at <- function(l, psi, route) {
  if (!.gradient_evaluate_feasible(route, psi)) {
    return(NULL)
  }
  .gradient_outputs(l)
}

# `route$evaluate(psi)`, FALSE where the leaf refuses `psi` as outside its
# feasible interval. Only that refusal is caught, by its infeasibility code;
# anything else is a failure of the evaluation and propagates.
.gradient_evaluate_feasible <- function(route, psi) {
  tryCatch({
    route$evaluate(psi)
    TRUE
  }, error = function(e) {
    if (identical(.infeasible_code(conditionMessage(e)), "prescribed_potential")) {
      return(FALSE)
    }
    stop(e)
  })
}

# The differentiable parameters, in the order C++ indexes them, derived from
# leaf_traits() so they cannot drift from the trait vector. C++ keeps its own
# copy (`gradient::par_names`) and test-gradient-batch.R compares the two: R
# passes positions, so appending is safe and reordering differentiates the wrong
# parameter (hazard 13). Computed at first call, for the collation reason above.
.gradient_par_names <- local({
  nms <- NULL
  function() {
    if (is.null(nms)) {
      nms <<- c(names(.leaf_trait_defaults), .gradient_non_traits)
    }
    nms
  }
})

# The differentiable parameters that are not traits, in enumeration order after
# the traits: the one list of the slots `leaf_traits()` does not supply.
.gradient_non_traits <- c("leaf_specific_conductance_max", "resistance",
                          "CF77_lambda_", "TF24_floor_lambda_o",
                          "LeastCost_beta",
                          # the same as .gradient_damage_pars, spelled out
                          # because that is defined below
                          "A_", "I_", "D_", "topt_offset_", "dTopt_max_", "K_A_",
                          "k_i_", "k_i_baseline_", "k_rec_", "k_rec_Ea_",
                          "c_acclim_maint_", "c_repair_maint_",
                          "c_protect_maint_", "c_repair_flux_")

# The thermal-damage layer's inputs that reach the solve, in enumeration order:
# Leaf fields, read and written by name, available only with
# `use_thermal_damage_` on. Its dynamics parameters move no solved output and are
# absent (their gradient is leaf_damage_gradient()'s), except repair's two, which
# the repair-flux cost reads.
.gradient_damage_pars <- c("A_", "I_", "D_", "topt_offset_", "dTopt_max_", "K_A_",
                           "k_i_", "k_i_baseline_", "k_rec_", "k_rec_Ea_",
                           "c_acclim_maint_", "c_repair_maint_",
                           "c_protect_maint_", "c_repair_flux_")

# Their current values, off a leaf, in that order.
.gradient_damage_values <- function(l) {
  vapply(.gradient_damage_pars, function(nm) l[[nm]], numeric(1))
}

# What `pars` may name, shared by leaf_gradient() and leaf_gradient_batch(). A
# parameter is unavailable on one of two axes, the supply path (`resistance` is
# single-potential only) or the model (an owned parameter belongs to one curve),
# and the refusal names which. Refusing is not a zero: a zero means "in the
# objective but inactive here", which is information.
#
# The one table mapping a model-owned parameter to its owner (hazard 13); its
# twin on the solve side is `.solve_price_fields`. Not all are prices:
# `LeastCost_beta` is a dimensionless cost ratio. A new owned parameter needs a
# row here, a row in `.solve_price_fields` and a write in `.gradient_setter`.
.gradient_owned_pars <- c(CF77 = "CF77_lambda_",
                          TF24_floor = "TF24_floor_lambda_o",
                          LeastCost = "LeastCost_beta")

.gradient_model_pars <- function(model) {
  if (!(model %in% .gradient_models())) {
    return(NULL)
  }
  unname(.gradient_owned_pars[names(.gradient_owned_pars) != model])
}

.gradient_available_pars <- function(single, model = "collar", damage = FALSE) {
  nms <- .gradient_par_names()
  if (!single) nms <- setdiff(nms, "resistance")
  if (!damage) nms <- setdiff(nms, .gradient_damage_pars)
  excluded <- .gradient_model_pars(model)
  if (is.null(excluded)) {
    stop("no gradient route for model \"", model, "\". Available: ",
         paste(.gradient_models(), collapse = ", "), ".", call. = FALSE)
  }
  setdiff(nms, excluded)
}

.gradient_check_pars <- function(pars, single, model = "collar", damage = FALSE) {
  unknown <- setdiff(pars, .gradient_available_pars(single, model, damage))
  if (length(unknown)) {
    off <- intersect(unknown, .gradient_damage_pars)
    if (length(off) && !damage) {
      stop("`pars` names `", off[[1L]], "`, an input of the thermal-damage layer, ",
           "which is off on this leaf. Set use_thermal_damage_ to TRUE on the leaf ",
           "you pass (`x`, or the batch's leaf) to differentiate it.",
           call. = FALSE)
    }
    why <- if (!single && "resistance" %in% unknown) {
      paste(" `resistance` is a parameter of the single-potential path only;",
            "on the multi-layer path the soil-to-collar resistances are derived",
            "from root carbon.")
    } else {
      ""
    }
    owned <- intersect(unknown, .gradient_owned_pars)
    if (length(owned)) {
      par <- owned[[1L]]
      owner <- names(.gradient_owned_pars)[match(par, .gradient_owned_pars)]
      stop("`pars` names `", par, "`, which is not a parameter of the ", model,
           " model: it is the ", owner, " curve's PRESCRIBED price of water, ",
           "and on every other cost curve the price is emergent -- derived from ",
           "that curve's own parameters rather than set. Differentiate those ",
           "instead, or pass model = \"", owner, "\".", call. = FALSE)
    }
    stop("`pars` names things this cannot differentiate: ",
         paste(unknown, collapse = ", "),
         ". Available: the traits from leaf_traits(), plus ",
         "`leaf_specific_conductance_max`",
         if (single) " and `resistance`" else "",
         ".", why, call. = FALSE)
  }
  invisible(pars)
}

# Everything this can differentiate, with its current value, as one named vector:
# the traits, then the non-traits. `leaf_specific_conductance_max` and
# `resistance` are drivers (set_drivers()); `resistance` is read from
# `r_R_V_sum[1]` and exists only on the single-potential path, where it is one
# scalar. The derivation does not care whether a parameter is a trait; only the
# setter differs (.gradient_setter).
.gradient_theta <- function(traits, kmax, supply, root_network,
                            CF77_lambda = NA_real_,
                            TF24_floor_lambda_o = NA_real_,
                            LeastCost_beta = NA_real_,
                            damage = NULL) {
  theta <- c(unlist(traits), leaf_specific_conductance_max = kmax)
  if (identical(supply$kind, "single")) {
    theta <- c(theta, resistance = root_network$r_R_V_sum[[1]])
  }
  # NA unless the caller is on the route that owns it. C++ reads every slot, so
  # each holds the value the solve uses; off its route the value is never read.
  theta <- c(theta, CF77_lambda_ = CF77_lambda,
             TF24_floor_lambda_o = TF24_floor_lambda_o,
             LeastCost_beta = LeastCost_beta)
  # The thermal-damage inputs, read off the leaf when its gate is on and NA
  # otherwise -- in which case neither side reads them.
  if (is.null(damage)) {
    damage <- stats::setNames(rep(NA_real_, length(.gradient_damage_pars)),
                              .gradient_damage_pars)
  }
  c(theta, damage)
}


# The arguments `.resolve_drivers()` takes, in its order, so the drivers list can
# be handed over without naming them.
.resolve_names <- c("psi_soil", "PPFD", "soil_depth", "root_network",
                    "leaf_specific_conductance_max", "atm_vpd", "ca",
                    "leaf_temp", "atm_o2_kpa", "atm_kpa")

# A closure that pushes a parameter vector onto the leaf. Order is load-bearing
# (hazard 10): set_traits() returns the leaf to its just-constructed state, so
# the drivers go on after it, and that write re-derives vcmax_/jmax_/R_d_.
.gradient_setter <- function(l, traits, drivers, supply,
                             fast_stem_curve = TRUE, damage = FALSE) {
  trait_names <- names(traits)
  single <- identical(supply$kind, "single")
  # The drivers are resolved once, by the one definition of their rules; the
  # closure applies them positionally, varying only what a perturbation moves.
  # A second copy of the rules here would drift from set_drivers().
  base <- do.call(.resolve_drivers, c(list(l), drivers[.resolve_names]))

  # Method closures pulled out once: `$` on an R6 object is a lookup per use, and
  # these run per perturbation.
  apply_traits <- l$set_traits
  apply_drivers <- l$set_physiology

  function(theta, only = NULL) {
    # Fast path for stem_P50: moves stem_b and leaves the seated physiology
    # alone. Sound only when `only` names this one parameter and the rest of the
    # leaf is at base, which .gradient_reseat_base() ensures.
    if (fast_stem_curve && identical(only, "stem_P50")) {
      l$perturb_stem_P50(theta[["stem_P50"]])
      return(invisible(l))
    }
    # Straight onto the object, bypassing set_traits()'s checks: sound only
    # because theta came from a caller's leaf_traits() and the C++ setter asserts
    # the positive-magnitude invariants. Do not copy this where values are not
    # known-good.
    #
    # A do.call over the derived trait names, and it must stay one (hazard 13):
    # spelled-out subscripts break at run time when a trait is added, surfacing
    # as an R-versus-C++ disagreement. test-cost.R counts the crossing.
    tv <- theta[trait_names]
    do.call(apply_traits, unname(as.list(tv)))
    # The thermal-damage inputs, ONLY with its gate on, and BEFORE the drivers:
    # set_physiology is what seats the block that reads them. A field write per
    # input costs a `.Call` each, which is why the gate-off path skips the loop.
    if (damage) {
      for (nm in .gradient_damage_pars) {
        l[[nm]] <- theta[[nm]]
      }
    }
    # `resistance` is a driver and goes in with the others, so set_traits() is
    # the only call here that resets the object.
    net <- if (single) {
      series_resistance(theta[["resistance"]])
    } else {
      base$root_network
    }
    apply_drivers(net, base$PPFD, base$psi_soil, base$soil_depth,
                  theta[["leaf_specific_conductance_max"]], base$atm_vpd,
                  base$ca, base$leaf_temp, base$atm_o2_kpa, base$atm_kpa)
    # Model-owned fields, each written only when set: off its route the slot is
    # NA, which would poison the leaf, and each write is a `.Call` test-cost.R
    # counts. Every owned slot needs its write here, or its column is an exact
    # zero in R while C++'s `gradient::apply()` reports a number.
    lam <- theta[["CF77_lambda_"]]
    if (!is.na(lam)) {
      l$CF77_lambda_ <- lam
    }
    lam_o <- theta[["TF24_floor_lambda_o"]]
    if (!is.na(lam_o)) {
      l$TF24_floor_lambda_o <- lam_o
    }
    beta <- theta[["LeastCost_beta"]]
    if (!is.na(beta)) {
      l$LeastCost_beta <- beta
    }
    invisible(l)
  }
}

# Parameters whose step is RELATIVE rather than floored at 1. See .gradient_step.
.gradient_relative_pars <- c("leaf_specific_conductance_max", "resistance")

# The step: relative above 1 and plain `step` below, so a small trait is still
# perturbed above the solve's ~1e-09 noise and a trait at zero is perturbed at
# all. Parameters far below 1 by nature (`leaf_specific_conductance_max`,
# 3.14e-05) take a plain relative step, or the floor would make it a 3% secant.
.gradient_step <- function(par, value, step) {
  floor <- if (par %in% .gradient_relative_pars) 0 else 1
  max(abs(value), floor) * step
}

# Where the stencil is centred. `g_min` is refused below zero and defaults to
# zero, so where its stencil would cross zero it shifts up one step and the
# difference is the derivative at `g_min = h`. Mirrored bit-for-bit in
# `gradient::centre_for()`.
.gradient_centre <- function(par, value, h) {
  if (identical(par, "g_min") && value - h < 0) value + h else value
}

# Parameters whose setter takes a SHORTCUT that leaves the rest of the object
# alone, and which therefore require the object to be at the base point before
# they run. See `.gradient_reseat_base()`.
.gradient_shortcut_pars <- function(fast_stem_curve) {
  if (fast_stem_curve) "stem_P50" else character(0)
}

# Every parameter's gradient is taken from the base point. Every setter but a
# shortcut (stem_P50's fast path) restores base on its way, so base is restored
# only before a shortcut that follows another parameter. Cheap: set_traits()
# rebuilds a table only when a `c` changed.
.gradient_reseat_base <- function(reset, theta, pars, fast_stem_curve) {
  shortcut <- pars %in% .gradient_shortcut_pars(fast_stem_curve)
  at_base <- TRUE
  function(k) {
    if (shortcut[[k]] && !at_base) {
      reset(theta)
    }
    at_base <<- FALSE
    invisible(NULL)
  }
}

# The implicit-function composite. Two perturbed evaluations per parameter,
# neither of which re-solves the model: `dprofit` at the UNPERTURBED psi* gives
# the mixed partial, and the outputs at that same psi* give the direct term.
.gradient_ift <- function(l, reset, theta, pars, psi_star, H, dY_dpsi, step,
                          fast_stem_curve = TRUE, dpsi_dtheta = NULL,
                          envelope = FALSE, route) {
  seat <- .gradient_reseat_base(reset, theta, pars, fast_stem_curve)
  # The mixed partials are returned too: a caller integrating its own sensitivity
  # of psi needs `M` and `H`, which cannot be recovered once divided.
  M <- numeric(length(pars))
  out <- t(vapply(seq_along(pars), function(k) {
    p <- pars[[k]]
    seat(k)
    h <- .gradient_step(p, theta[[p]], step)
    centre <- .gradient_centre(p, theta[[p]], h)
    side <- function(sign) {
      th <- theta
      th[[p]] <- centre + sign * h
      reset(th, only = p)
      y <- .gradient_outputs_at(l, psi_star, route)
      if (is.null(y)) {
        stop("leaf_gradient(): perturbing `", p, "` moved the feasible collar ",
             "interval past psi*, so the operating point could not be ",
             "evaluated there. This point is on an active-set boundary; ",
             "lower `stationarity_tol` or difference the solve directly.",
             call. = FALSE)
      }
      c(dprofit = route$dprofit(psi_star), y)
    }
    up <- side(1)
    dn <- side(-1)
    # M = d2profit/dpsi dtheta, with psi held FIXED at psi*.
    M[[k]] <<- (up[["dprofit"]] - dn[["dprofit"]]) / (2 * h)
    # Solved: psi* moves with theta as the implicit function theorem says.
    # Imposed: it moves only as the caller says, by `dpsi_dtheta`.
    d_psi <- if (is.null(dpsi_dtheta)) -M[[k]] / H else dpsi_dtheta[[k]]
    direct <- (up[-1L] - dn[-1L]) / (2 * h)
    # The decision variable IS psi*, held fixed, so its gradient is assigned
    # dpsi*/dtheta rather than read off a difference of two identical numbers.
    g <- direct + dY_dpsi * d_psi
    g[[route$decision]] <- d_psi
    # The envelope term, assigned for the same reason: profit's indirect term is
    # identically zero at a stationary point.
    if (envelope) {
      g[["profit"]] <- direct[["profit"]]
    }
    g[.gradient_output_names()]
  }, numeric(length(.gradient_output_names()))))
  reset(theta)
  rownames(out) <- pars
  list(gradient = out, M = stats::setNames(M, pars))
}

# `psi` and `dpsi_dtheta` only mean anything together, and neither means anything
# alongside `method`.
.gradient_check_psi <- function(psi, dpsi_dtheta, method) {
  if (is.null(psi)) {
    if (!is.null(dpsi_dtheta)) {
      stop("`dpsi_dtheta` is the trait response of a collar potential YOU ",
           "imposed, so it needs `psi`. Without one the collar is solved for ",
           "and its response is derived by the implicit function theorem.",
           call. = FALSE)
    }
    return(NULL)
  }
  if (!(is.numeric(psi) && length(psi) == 1L && is.finite(psi) && psi > 0)) {
    stop("`psi` must be a single finite positive number: it is a collar water ",
         "potential in MPa, as a positive magnitude (see the sign convention ",
         "in ?leaf_model).", call. = FALSE)
  }
  if (!identical(method, "auto")) {
    stop("`method` cannot be given with `psi`. The two routes it chooses ",
         "between are about a SOLVED operating point -- \"fd\" differences the ",
         "solve, which would answer about the optimum rather than about the ",
         "collar you gave -- so at a prescribed psi neither applies and the ",
         "reported method is \"prescribed\".", call. = FALSE)
  }
  # Coerced, and the caller must use the return, so the R reference and the C++
  # batch see the same double; `.gradient_check_psi_batch()` does the same.
  as.numeric(psi)
}

# dpsi/dtheta as one value per parameter, in `pars` order. NULL on the solving
# path means "derive it"; NULL on the prescribed path means "the collar does not
# move with theta", which is zero and not the same statement.
.gradient_dpsi_dtheta <- function(dpsi_dtheta, pars, prescribed) {
  if (!prescribed) {
    return(NULL)
  }
  if (is.null(dpsi_dtheta)) {
    return(stats::setNames(numeric(length(pars)), pars))
  }
  # Finite, not just non-NA: an infinite value makes four columns +-Inf while
  # `profit`, assigned by the envelope theorem, stays finite and plausible.
  if (!is.numeric(dpsi_dtheta) || !all(is.finite(dpsi_dtheta))) {
    stop("`dpsi_dtheta` must be numeric and finite -- no NA, NaN or Inf. It ",
         "is dpsi/dtheta for a collar YOU imposed; a non-finite value returns ",
         "an all-infinite gradient row with a plausible `profit` in it.",
         call. = FALSE)
  }
  # Names are matched before recycling, or a named value would be applied to
  # whichever parameter recycling put it on.
  if (!is.null(names(dpsi_dtheta))) {
    if (!setequal(names(dpsi_dtheta), pars)) {
      stop("`dpsi_dtheta` is named, so its names must be exactly `pars`. ",
           "Missing: ", paste(setdiff(pars, names(dpsi_dtheta)),
                              collapse = ", "),
           "; unexpected: ", paste(setdiff(names(dpsi_dtheta), pars),
                                   collapse = ", "), call. = FALSE)
    }
    return(dpsi_dtheta[pars])
  }
  if (length(dpsi_dtheta) == 1L) {
    return(stats::setNames(rep(dpsi_dtheta, length(pars)), pars))
  }
  if (length(dpsi_dtheta) != length(pars)) {
    stop("`dpsi_dtheta` must be length 1 or one value per parameter (",
         length(pars), " for this `pars`); got ", length(dpsi_dtheta),
         call. = FALSE)
  }
  stats::setNames(dpsi_dtheta, pars)
}

# The fallback: a central difference of the whole solve. Correct at a pinned
# optimum because it differences the constrained answer, which is exactly what the
# composite cannot do.
.gradient_fd <- function(l, reset, theta, pars, step, fast_stem_curve = TRUE,
                         route) {
  # The route's step unless the caller set one: this arm differences the solve,
  # so its step must clear the solver's resolution (`fd_step`).
  if (is.null(attr(step, "user_set"))) {
    step <- route$fd_step
  }
  seat <- .gradient_reseat_base(reset, theta, pars, fast_stem_curve)
  out <- t(vapply(seq_along(pars), function(k) {
    p <- pars[[k]]
    seat(k)
    h <- .gradient_step(p, theta[[p]], step)
    centre <- .gradient_centre(p, theta[[p]], h)
    side <- function(sign) {
      th <- theta
      th[[p]] <- centre + sign * h
      reset(th, only = p)
      route$solve()
      .gradient_outputs(l)
    }
    ((side(1) - side(-1)) / (2 * h))[.gradient_output_names()]
  }, numeric(length(.gradient_output_names()))))
  rownames(out) <- pars
  reset(theta)
  out
}
