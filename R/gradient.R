# Trait gradients (#4).
#
# dA/dtheta, dgs/dtheta and dpsi_leaf/dtheta at a solved operating point, for the
# traits. The customer is a companion calibration study, whose hierarchical fit spent 97.6%
# of 10,045 model evaluations on finite-difference gradients.
#
# WHAT THIS FILE IS FOR, in one paragraph, because the arithmetic that motivated
# it did not survive being measured. The plan was that the implicit function
# theorem would replace 2N re-solves per gradient with 2N evaluations of the exact
# `dprofit`, at 7% of the cost of a solve, for a 10x win in the R layer. The
# theorem is right and it is implemented below. The 10x is not there, and the
# reason is worth stating once: in R the per-call overhead (~1.8 us) is seven
# times the C++ work it wraps (0.26 us), so replacing a 3.6 us solve with a
# 0.26 us gradient evaluation saves almost nothing measurable from R. What DOES
# dominate is that perturbing a trait used to mean constructing a new Leaf, which
# costs ~180x a trivial `.Call`, or 45 solves. `set_traits()` removed that, and it is where the
# speedup actually came from.
#
# ⚠️ THAT PARAGRAPH IS TRUE AND IT WAS OVER-GENERALISED, INCLUDING HERE. "The
# value of this file is exactness, not speed" held for the case measured -- a fit
# whose parameters ARE traits, one per fitted parameter. It is false for a fit
# where they are not. Both arms are linear in a parameter count, and it is a
# DIFFERENT count for each:
#
#   T_fd    = 2 * P_fit   * T_objective
#   T_exact = N * (t_construct + t_solve)  +  2 * P_model * N * t_eval
#
# `optim` differences the P_fit parameters it is moving; `leaf_gradient()` is asked
# for the P_model parameters the leaf actually has. Those are equal only when the
# fit varies traits directly. Any parameterisation in between -- pooling, a
# hierarchy, a shared or derived parameter -- makes P_fit > P_model, and the
# composite's cost does not follow the optimiser's dimension.
#
# Measured both regimes: vignette("fitting")'s 72-observation trait fit has
# P_fit == P_model and the composite loses at every count up to 13. The companion
# calibration study has N = 1327, P_fit = 40, P_model = 4 and it wins 3.4x,
# converging to the same optimum -- and its 57-parameter variant costs the SAME as
# its 40-parameter one, because P_model is 4 in both. One set of per-observation
# coefficients fitted on the first design predicts the second to ~5%.
#
# So: the value of the code below is EXACTNESS, the ACTIVE-SET CLASSIFICATION, and
# speed WHEN P_fit EXCEEDS P_model. The classification is still the part that took
# the work. Do not quote a speed verdict from this file without saying which of the
# two counts it was denominated in -- that conflation has now cost this family a
# retraction in each direction.

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
##' The reason it is a function rather than fifteen assignable fields is that
##' `leaf$vcmax_25 <- x` could not be made correct. Three separate pieces of
##' derived state go stale on a bare trait write: the two pre-integrated
##' vulnerability splines, the solved operating point, and -- least visibly --
##' `vcmax_`/`jmax_`/`R_d_`, which `set_physiology()` computes behind a cache keyed
##' on leaf temperature and O2 alone. That last one means the obvious repair,
##' "change the trait and then set the drivers again", silently does not work.
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
  # Positional, in the C++ argument order, for the same reason leaf_model() is:
  # the ordering is written down once rather than at every call site.
  x$set_traits(traits$vcmax_25, traits$stem_c, traits$stem_P50,
               traits$root_c, traits$root_P50, traits$TF24_beta2,
               traits$jmax_25, traits$a, traits$curv_fact_elec_trans,
               traits$curv_fact_colim, traits$TF24_cost_scale, traits$R_d_25,
               traits$JS22_gamma, traits$CMax_a, traits$CMax_b)
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
##' The second term is not a correction. For `TF24_cost_scale`, `TF24_beta2`, `stem_b`
##' and `stem_c` it is 100% of the answer, and for `vcmax_25` 52%.
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
##' ⚠️ **A clamped `psi` gets no gradient, and that is the interesting case.**
##' The collar actually used is `psi` clamped into the feasible interval, so it
##' moves with the *bound* rather than with `dpsi_dtheta` — the active-set
##' problem, arriving through the clamp instead of through the optimiser. The
##' direct term alone would be plausible and wrong. `status` reports `"clamped"`
##' and the gradient is all `NA`; `psi` in the result is the collar that was
##' used. This is reported rather than thrown because a tracking model reaches
##' these points routinely — the clamp is how TF24f pulls an out-of-range state
##' back inside — and it fires for a `psi` within one step of an end too, since
##' `dY/dpsi` cannot be centred there.
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
##' It is better in a second way, which is sharper than accuracy. `psi_crit` does
##' not appear in the profit function at all -- it only sets the dry end of the
##' feasible collar interval -- so the composite necessarily returns **exactly
##' zero** for it. At an interior optimum that is correct. At a boundary-crit optimum
##' `psi_crit` *is* the binding constraint, the true `dA/dθ` is about `1.26`, and a
##' zero is arguably a worse answer than a wildly wrong one: it tells an optimiser
##' the parameter does nothing, and nothing about it looks suspicious. The
##' difference of the solve gets it for free.
##'
##' There is a second guard behind that one, found by measurement rather than
##' designed: a pinned optimum sits so close to its bound that the difference in
##' `psi` cannot be centred there without being clamped, and that is detected too.
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
##' trait perturbation can. That is the case the clamp detector below exists for,
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
##' ⚠️ **Always pass `pars`.** It is `P_model`, so the default — all fifteen on
##' the multi-layer path — is the most expensive thing you can ask for, and a fit
##' that reads four of them pays for twelve it discards.
##'
##' The intercept is a fresh `Leaf` per call, ~45 solves' worth. Pass `x` to reuse one
##' across calls (#52); without it, construction is 29% of the exact gradient in
##' the study above.
##'
##' @section Precision:
##' Do not ask for more than about `1e-09` from any of this. That is the
##' achievable precision of the solved outputs themselves, set by the tolerance of
##' the intercellular-CO2 root-find.
##'
##' @section Reusing a leaf across gradients:
##' Constructing a `Leaf` costs about **180x a trivial `.Call`** — 45 solves — and
##' this function does it once per call, which is **61% of a one-parameter
##' gradient**: the largest single term on this surface by some way, and paid once
##' per observation by a fit that differentiates per observation.
##'
##' ⚠️ **That share has GROWN, and the reason is worth knowing.** It was quoted at
##' 40% when the solvers searched their objectives; they root-find their
##' first-order conditions now, so the composite got cheaper while construction did
##' not. Re-measure it as a RATIO in one process — absolute microseconds here moved
##' by 50% between machines, and four places in this package quoted the same
##' absolute with four different solve counts beside it. Pass `x` to reuse one:
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
##' The restore costs one re-trait and one solve, about 4 solves against the 45
##' saved.
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
##' @param pars what to differentiate with respect to. Any of the fifteen
##'   [leaf_traits()] names, plus `"leaf_specific_conductance_max"`, on the
##'   single-potential path only `"resistance"`, on `model = "CF77"` only
##'   `"CF77_lambda_"`, and on `model = "TF24_floor"` only `"TF24_floor_lambda_o"`.
##'   Defaults to all the route allows.
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
##' @param fast_stem_curve perturb `stem_b` without rebuilding the stem
##'   vulnerability spline (`TRUE`, the default). The curve is homogeneous of
##'   degree 1 in `stem_b`, so the spline for a perturbed `stem_b` is the existing
##'   one with its argument rescaled — which makes the rebuild, otherwise most of
##'   the cost of that gradient, unnecessary. `FALSE` rebuilds, and exists so the
##'   equivalence is checkable rather than assumed; the package's tests require
##'   the two routes to agree. `stem_c` has no such identity and always rebuilds.
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
##'       solving path; `"prescribed"`, `"clamped"` or `"no-gradient"` with
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
##'       solving path, and on the prescribed path the value actually used,
##'       which differs from the `psi` argument exactly when `status` is
##'       `"clamped"`}
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
  # ⚠️ THE RETURN IS THE POINT, NOT JUST THE CHECK. `psi` is compared against
  # `opt_root_psi_` to detect the clamp, and a caller writing `psi = 3L` -- or the
  # entirely ordinary `for (p in 2:5)`, since `2:5` is integer -- would otherwise
  # be reported CLAMPED at a collar it was given exactly. The checker coerces;
  # discarding what it returns was the bug.
  psi <- .gradient_check_psi(psi, dpsi_dtheta, method)

  # --- reusing a Leaf (#52) -------------------------------------------------
  # `x` is a VESSEL, not the point being differentiated: `traits` still says where
  # the gradient is taken, and the setter applies it. That is why `traits` becomes
  # REQUIRED here rather than optional -- with `x` supplied and `traits` left at its
  # default, this would differentiate at the package defaults on somebody else's
  # leaf and return a plausible answer for a point nobody asked about.
  #
  # `control` and `supply` are baked into `x`, so accepting them alongside it would
  # invite a silent disagreement (a `supply` naming one path, an `x` built on the
  # other). Rejected rather than resolved. `supply` is still needed internally --
  # only its `kind` is ever read -- so it is derived from the object.
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

  # Resolve the single path's default resistance HERE rather than letting
  # set_drivers() do it, because `resistance` is one of the parameters being
  # differentiated: theta has to hold the value the solve will actually use, and a
  # NULL would leave it unknown. `.default_series_resistance()` is the same
  # function set_drivers() would have called, so there is no second default to
  # drift.
  if (identical(supply$kind, "single") && is.null(root_network)) {
    root_network <- .default_series_resistance()
  }

  drivers <- list(psi_soil = psi_soil, PPFD = PPFD, soil_depth = soil_depth,
                  root_network = root_network,
                  leaf_specific_conductance_max = leaf_specific_conductance_max,
                  atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp,
                  atm_o2_kpa = atm_o2_kpa, atm_kpa = atm_kpa)

  # The differentiable parameters: the fifteen traits, plus the four that are not
  # traits and that a calibration nonetheless fits (#44). See .gradient_theta.
  # ⚠️ READ OFF THE LEAF, not defaulted, for the same reason the resistance is:
  # theta has to hold the value the solve will actually use, and on the two priced
  # routes that is whatever the caller assigned to `x$CF77_lambda_` /
  # `x$TF24_floor_lambda_o` before calling.
  theta <- .gradient_theta(traits, leaf_specific_conductance_max, supply,
                           root_network,
                           CF77_lambda = if (is.null(x)) NA_real_ else
                             x$CF77_lambda_,
                           TF24_floor_lambda_o = if (is.null(x)) NA_real_ else
                             x$TF24_floor_lambda_o)
  if (is.null(pars)) {
    # ⚠️ WHAT IS AVAILABLE, not every slot in `theta`. The enumeration now carries
    # a parameter that belongs to one model, so "all of them" and "all of the ones
    # this route can differentiate" have stopped being the same set -- and taking
    # the former made the default refuse itself.
    pars <- .gradient_available_pars(identical(supply$kind, "single"), model)
  }
  # Shared with leaf_gradient_batch(), so the two entry points cannot disagree
  # about which parameters exist or explain a rejection differently.
  .gradient_check_pars(pars, identical(supply$kind, "single"), model)
  # ⚠️ HERE, NOT WHERE IT IS USED. `dpsi_dtheta` is only READ inside the
  # composite, which a clamped or shut-down operating point never reaches -- so
  # validating it there made a wrong-length vector an error at some psi and
  # silently accepted at others. Argument checking cannot be conditional on what
  # the model does with the arguments.
  dpsi_dtheta <- .gradient_dpsi_dtheta(dpsi_dtheta, pars, !is.null(psi))

  # ONE leaf for the whole gradient, re-traited rather than reconstructed. This is
  # the measurement that decided it: a fresh Leaf costs ~45 solves against ~3.4 to
  # re-trait and re-drive this one, so reconstructing per perturbation would swamp
  # any difference between the two gradient routes below.
  #
  # And with `x`, one leaf across gradients too -- which is #52, and worth ~40% of
  # a call in a per-observation fit. Construction is the single largest term on
  # this surface.
  l <- if (is.null(x)) leaf_model(traits, control, supply) else x
  reset <- .gradient_setter(l, traits, drivers, supply, fast_stem_curve)

  if (is.null(x)) {
    do.call(set_drivers, c(list(l), drivers))
  } else {
    # ⚠️ A REUSED LEAF MUST BE SEATED WITH reset(), NOT set_drivers() ALONE. It
    # arrives carrying whatever traits its last user left on it, and everything
    # below -- `value`, `psi_star`, `H`, and the whole composite -- is supposed to
    # describe `traits`. set_drivers() re-derives the drivers and never touches
    # traits, so a leaf built at other traits would be differentiated at the wrong
    # point: plausibly, and with no symptom. reset() applies both.
    reset(theta)
    # Hand the caller's object back as this function found it -- solved at the base
    # point -- rather than at whichever perturbation happened to run last, which is
    # hazard 8 territory. `on.exit` because several paths below stop().
    on.exit({
      try({
        reset(theta)
        l$find_root_collar_psi()
      }, silent = TRUE)
    }, add = TRUE)
  }
  # --- seat the operating point ---------------------------------------------
  # Two ways in, and which one ran is what everything below branches on. The
  # default SOLVES for the collar potential; `psi` IMPOSES one, which is what a
  # caller tracking the optimum rather than finding it has (#88).
  #
  # ⚠️ `psi_star` keeps its name on both paths and it is no longer always the
  # argmax. It is "the collar the outputs were evaluated at", which is what every
  # use of it below actually means.
  route <- .gradient_route(l, model, supply)
  prescribed <- !is.null(psi)
  if (prescribed) {
    route$evaluate(psi)
    # Exact equality, for the reason `.gradient_outputs_at` documents: the clamp
    # is a min/max, so an unclamped target returns bit-identically.
    clamped <- !identical(route$psi_star(), psi)
    psi_star <- route$psi_star()
  } else {
    route$solve()
    clamped <- FALSE
    psi_star <- route$psi_star()
  }
  value <- .gradient_outputs(l)

  # The curvature and the residual, at whichever collar was seated. Both are
  # differences of `dprofit_droot_collar_psi`, which takes psi as an ARGUMENT --
  # so neither needs the point to be an optimum, and both mean the same thing on
  # the two paths.
  h_psi <- max(abs(psi_star), 1) * step
  # ⚠️ WITH ITS FEASIBILITY, not bare. `dprofit_droot_collar_psi` returns a hard
  # 0.0 SENTINEL on its shut-down and reversed-gradient exits, and a bare zero is
  # indistinguishable from a stationary point -- the header says so, and says a
  # composite that ignores the flag inherits the bug. The solving path got away
  # with reading the value alone because `H` collapses to zero too and `usable`
  # catches the pair; the prescribed path does NOT divide by `H`, so it would
  # have adopted the sentinel as if it were dprofit/dpsi.
  checked <- route$checked(psi_star)
  resid <- checked[[1L]]
  feasible <- checked[[2L]] == 1
  H <- (route$dprofit(psi_star + h_psi) -
        route$dprofit(psi_star - h_psi)) / (2 * h_psi)
  # H == 0 with resid == 0 is the shut-down signature: dprofit returns a sentinel
  # zero there rather than a derivative, so the ratio below would be 0/0. H > 0
  # would not be a maximum. Both mean the composite has nothing to stand on.
  usable <- is.finite(H) && H < 0 && is.finite(resid)
  stationarity <- if (usable) abs(resid / H) else Inf

  if (prescribed) {
    # ⚠️ THE STATIONARITY TEST DOES NOT ROUTE HERE, AND IS STILL WORTH TAKING.
    # What it decides on the solving path -- composite or fallback -- is
    # meaningless at a collar the caller chose: there is no argmax to be pinned
    # against, and differencing the solve would answer a question about the
    # optimum instead of about this point. So `method` is rejected upstream.
    #
    # The NUMBER keeps its meaning, though, and gains a better one: it is how far
    # the point you handed over sits from the optimum, in MPa. It is reported,
    # and it is used for exactly one thing -- see `envelope` below.
    #
    # ⚠️ `no-gradient` REACHES THIS PATH TOO, and it is not the same condition as
    # the solving path's. There, `usable` also demands `H < 0` -- a MAXIMUM test,
    # which is exactly what a caller-chosen collar has no business satisfying: a
    # prescribed psi away from the optimum may sit where profit is convex, and
    # that is fine, because nothing here divides by `H`. What genuinely disables
    # the point is INFEASIBILITY: the shut-down and reversed-gradient exits, where
    # `dprofit` is a sentinel rather than a derivative. Almost every such point is
    # already caught as `clamped` -- the shut-down state seats a collar of its own
    # choosing -- but "almost" is not a guarantee, since a caller can pass exactly
    # that collar back.
    status <- if (clamped) "clamped" else if (!feasible) "no-gradient" else
      "prescribed"
  } else {
    status <- if (!usable) "no-gradient" else
      if (stationarity > stationarity_tol) "pinned" else "interior"
  }

  # `status` describes the POINT and is reported whichever route runs; `use_ift`
  # is the route. They differ only when the caller has forced one.
  # ⚠️ ProfitMax's TWO METHODS ANSWER DIFFERENT QUESTIONS, and `auto` picks the one
  # a caller almost always wants. Its `Scaled` link divides by `|A|max`, which the
  # composite holds FIXED across a perturbation -- so `"ift"` returns the PARTIAL at
  # fixed normaliser. Differencing the solve lets the scan re-run, so `"fd"` returns
  # the TOTAL. Measured at psi_soil 1.5, PPFD 1500, for `vcmax_25`: partial 0.0766,
  # total 0.0579, and `d|A|max/dvcmax_25` = 0.0768 -- the omitted chain term is 32%
  # of the partial, not a rounding difference.
  #
  # A FIT needs the total: the model re-scans `|A|max` at every parameter vector, so
  # an optimiser fed the partial stops where the partial vanishes rather than where
  # the objective does. Hence `auto` differences the solve here, and `"ift"` remains
  # available for anyone who wants the partial deliberately.
  #
  # No other curve has this split -- for TF24 the two agree to 0.0%.
  auto_ift <- identical(status, "interior")
  use_ift <- if (prescribed) !clamped && feasible else switch(method,
                    auto = auto_ift,
                    ift = TRUE,
                    fd = FALSE)
  if (use_ift && !prescribed && !usable) {
    stop("leaf_gradient(): method = \"ift\" was asked for at a point with no ",
         "usable curvature (H = ", format(H), "), so -M/H has nothing to stand ",
         "on. This is a shut-down or otherwise determined operating point; use ",
         "method = \"auto\".", call. = FALSE)
  }
  # ⚠️ A CLAMPED PRESCRIBED PSI GETS NO GRADIENT, RATHER THAN THE DIRECT TERM.
  # It is not a failure -- the outputs at the clamped collar are perfectly good,
  # and TF24f relies on the clamp to pull an out-of-range tracked state back
  # inside. It is that the derivative is not the one this can compute: the collar
  # actually used is `min(max(psi, a(theta)), b(theta))`, so it moves with the
  # BOUND, and dY/dtheta picks up the bound's derivative rather than the caller's
  # `dpsi_dtheta`. That is the active-set problem arriving through the clamp
  # instead of through the optimiser, and the direct term alone would be
  # plausible and wrong in exactly the documented way.
  #
  # Reported rather than thrown, because a fit will visit these points routinely
  # and `status` is how a batch tells its caller which rows to distrust. The
  # all-NA result is assembled below, with the second way of reaching it.

  if (use_ift) {
    # dY/dpsi at fixed traits. evaluate_root_collar_psi CLAMPS its target into the
    # feasible interval, so check it landed where it was asked to: a clamped
    # evaluation would silently make this a one-sided difference over a shorter
    # interval, which is the same class of error as the pinned case.
    #
    # This turns out to be a SECOND, INDEPENDENT detector of a pinned optimum
    # rather than the unreachable guard it was written as, and the measurement is
    # worth recording. At a pinned point psi* sits one step-in fraction (1e-06 of
    # the bracket width) from the bound, so a step of `step` * psi in psi crosses
    # it whenever the bracket is narrower than psi -- which every pinned row in
    # the package's grid is, being at the dry end where the feasible interval has
    # nearly closed. Measured: forcing method = "ift" fails here on all 42 pinned
    # rows and on all 48 shut-down ones, so the composite's silently-wrong answer
    # is not reachable through this function at all.
    #
    # It is NOT a substitute for the stationarity test, and reading it as one
    # would be the mistake: it fires only when the bracket is narrow, so a pinned
    # optimum on a wide bracket would pass it. The stationarity test is the one
    # that is scale-free and the one the default relies on.
    hi <- .gradient_outputs_at(l, psi_star + h_psi, route)
    lo <- .gradient_outputs_at(l, psi_star - h_psi, route)
    if (is.null(hi) || is.null(lo)) {
      if (prescribed) {
        # Same reasoning as the clamp on `psi` itself, one step out: a psi that
        # is inside the interval but within `h_psi` of an end cannot have dY/dpsi
        # centred on it, and a one-sided difference over a shortened interval is
        # the failure this detector exists for. There is no fallback to offer --
        # differencing the solve would answer about the optimum -- so the row is
        # reported as clamped.
        clamped <- TRUE
        status <- "clamped"
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
      # `dprofit_droot_collar_psi` is EXACT in psi -- forward AD plus the IFT at
      # the ci root-find -- so for profit alone the package has something better
      # than a difference of the same quantity, and it is already computed. The
      # other four have no such route and must be differenced. One rule, both
      # paths, which is what keeps `psi = psi*` reproducing the solve exactly.
      #
      # ⚠️ THIS IS THE READER #87 SAID DID NOT EXIST YET. Until the prescribed path
      # landed, the only consumer was `.gradient_ift(envelope = FALSE)` -- a forced
      # method = "ift" at a pinned point, which throws at all 42 pinned rows of the
      # grid. A prescribed psi away from the optimum is not stationary, so it takes
      # this branch for real, and the exactness now matters.
      dY_dpsi[["profit"]] <- resid
    }
  }

  # ⚠️ THE ENVELOPE THEOREM, and the only place this package uses it. At a
  # STATIONARY point dprofit/dpsi is analytically zero, so profit's indirect term
  # vanishes identically and dprofit/dtheta is the direct partial alone. The
  # composite is told to ASSIGN that column rather than reach it by multiplying a
  # near-zero dY/dpsi -- the same treatment `collar` gets, and for the same
  # reason: an identity is stated, not arrived at.
  #
  # ⚠️ The test is `stationarity`, NOT `status`, so the two paths agree at a psi
  # that happens to BE psi*: `status` carries "prescribed" on one and "interior"
  # on the other, and routing on it would break the bit-for-bit equivalence while
  # looking equivalent. It is also not `use_ift`: at a pinned optimum psi* is a
  # theta-dependent BOUND, dprofit/dpsi is not zero, and the indirect term
  # survives -- someone forcing method = "ift" there already gets a confidently
  # wrong number and should not get a differently wrong one for this column.
  #
  # The measured size of what this removes is in ?leaf_gradient.
  envelope <- isTRUE(usable && stationarity <= stationarity_tol)

  if (prescribed && !use_ift) {
    grad <- matrix(NA_real_, length(pars), length(.gradient_output_names()),
                   dimnames = list(pars, .gradient_output_names()))
    M <- stats::setNames(rep(NA_real_, length(pars)), pars)
    dY_dpsi <- stats::setNames(rep(NA_real_, length(.gradient_output_names())),
                               .gradient_output_names())
  } else if (use_ift) {
    # ONE composite for both paths, which is what makes the equivalence above a
    # real assertion rather than two implementations that happen to agree. The
    # only difference is where dpsi/dtheta comes from: derived by the implicit
    # function theorem when the collar was solved for, and supplied by the caller
    # when it was imposed -- because then it is not this function's to know.
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
# `leaf_gradient()` used to hardwire the collar solve, so it always differentiated
# the TF24 cost maximised over the root-collar potential whatever the caller had
# set up. That is now a choice, and it is made ONCE here rather than branched on at
# each of the six call sites.
#
# ⚠️ THE TWO ROUTES ARE DIFFERENT MODELS, not one restricted to the other. The
# collar route optimises the COLLAR potential with the root path in series; a stem
# route optimises `psi_stem` with the upstream potential pinned at `psi_soil` and
# ignores the soil-to-collar path entirely. Measured over 30 driver rows they
# disagree on 20, and the gap does NOT close as the root resistance goes to zero.
# So `model = "TF24"` and `model = "collar"` are not two ways of asking one
# question even though both use the TF24 cost.
#
# ⚠️ WHICH OUTPUT IS THE DECISION VARIABLE MOVES WITH THE ROUTE, and that is what
# the `decision` field is for. On the collar route `collar` IS psi*, so
# `dcollar/dtheta` is `dpsi*/dtheta` and is assigned rather than differenced. On a
# stem route `psi_stem` plays that part and `collar` is not solved for at all --
# `clear_collar_solve_state()` leaves it NA, and an NA gradient for it is the
# truthful answer rather than a gap.
.GRADIENT_COLLAR_ROUTE <- "collar"

# Every model this can differentiate: the collar route plus the seven cost
# curves' stem routes. ONE definition -- `.gradient_model_pars()` validates
# against it and `.gradient_route()` dispatches on it, where there used to be
# three lists that could disagree.
#
# ⚠️ EVERY ONE IS VERIFIED AGAINST A FINITE DIFFERENCE OF ITS OWN SOLVE, and the
# verification is a step SWEEP rather than one step. Skipping the sweep cost a
# false bug report: at a relative step of 1e-06 the FD of a SCAN-resolved argmax
# is pure quantisation -- psi* moves ~3e-06 MPa while the scan resolves it to
# ~4e-06 over a 4.4 MPa bracket. Measured for SOX's `vcmax_25`, the FD read
# 0.1855, 0.0552, 0.0551 at steps of 1e-06, 1e-03, 1e-02 while the composite
# returned a flat 0.0551 -- and for JW26 the 1e-06 step got the SIGN wrong. The
# composite was right in every case and the verification was what was broken. The
# stem routes root-find their first-order condition now, so one step serves all.
.gradient_models <- function() {
  c(.GRADIENT_COLLAR_ROUTE, cost_curve_names())
}

.gradient_route <- function(l, model, supply) {
  if (identical(model, .GRADIENT_COLLAR_ROUTE)) {
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

  # ⚠️ SEATED HERE, NOT IN `solve()`, and that is a correctness fix rather than
  # tidying. The four closures below read whichever model the leaf is set to --
  # they take no curve argument, because an integer index silently means a
  # different model when the enumeration grows. The PRESCRIBED path never calls
  # `solve()`, so seating there would have left `evaluate()` and `dprofit()`
  # answering about whatever model the leaf was last used for. Seating on
  # construction means every closure has one.
  #
  # It survives a perturbation: `set_traits()` clears the splines, the operating
  # point and the temperature cache, but not the seated curve.
  l$set_model(model, "stem")

  # ⚠️ ProfitMax RE-SOLVES ITS NORMALISER RATHER THAN PINNING IT, and that is what
  # makes the composite a TOTAL derivative. `|A|max` is state a perturbation
  # destroys -- `set_traits` + `set_physiology` clears it, measured 16.757 -> NaN --
  # so something has to put it back before every read.
  #
  # It used to be PINNED at the base point, for a reason that has since been fixed:
  # `|A|max` was the argmax of a 500-point scan, so it was piecewise constant in the
  # traits and a total derivative through it was zeros and jumps. It is found by a
  # root-find now, and cheaply (it sits at the dry bound on 1318 of 1320 driver
  # rows), so re-solving it lets it follow the traits -- which is the quantity a fit
  # needs, and removes the partial-versus-total split this route used to carry.
  #
  # A no-op on every other curve, which is why it is one function and not a pair
  # with a `capture` half that did nothing on either branch.
  reseat <- if (identical(model, "ProfitMax")) {
    function() l$prepare_profitmax()
  } else {
    function() invisible(NULL)
  }

  list(
    model = model,
    decision = "psi_stem",
    # THE SAME STEP AS THE COLLAR ROUTE, because there is now one solver. The stem
    # entry points root-find dJ/dpsi == 0 inside a scanned basin, so psi* is
    # stationary to ~1e-15 and a 1e-06 relative perturbation moves it far more than
    # the solve can resolve.
    #
    # ⚠️ THIS USED TO BE 1e-03 AND THE REASON IS WORTH KNOWING BEFORE ANYONE
    # TIGHTENS IT FURTHER. While these entry points refined on bracket WIDTH they
    # resolved psi* only to ~4e-06 MPa, so a 1e-06 step differenced quantisation:
    # measured for SOX's `vcmax_25`, 0.1855 against a converged 0.0551, and the
    # wrong SIGN for JW26. At the same step now: 2.5e-05 and 1.5e-05 relative.
    # The floor moves with the solver, so re-sweep if the solver changes again.
    #
    # ⚠️ ONE STEP FOR EVERY ROUTE, INCLUDING ProfitMax, AND THAT IS NEW. While
    # `|A|max` was the argmax of a 500-point scan it was piecewise constant in the
    # traits, so a step smaller than a scan cell differenced a staircase: 0.014 at
    # 1e-06 against 0.058 at 1e-02. `|A|max` is found by a root-find now, so the
    # difference is flat from 1e-02 to 1e-06 (0.057884 to 0.057885) and the
    # exception is gone. Sweep again if that normaliser ever goes back to a grid.
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

# The differentiated outputs, READ from C++ rather than restated here.
#
# ⚠️ This is deliberately NOT the pattern `.gradient_par_names` uses. That one is a
# second copy kept honest by a test comparing it with `gradient_par_names()`, and
# it has to be, because R builds `theta` from `leaf_traits()` and so needs the
# order before any C++ call. The OUTPUT list has no such constraint, so it is a
# single definition and there is no hazard comment to write: adding an output is
# one edit in `gradient.hpp`, not two that can disagree.
#
# `collar` is psi* itself, which makes dcollar/dtheta equal to dpsi*/dtheta -- so
# the two routes compute the same quantity by different means and a test can
# compare them. `profit` is what plant CONSUMES (`leaf.profit_`, not
# `assim_colimited_`) and is the one output the envelope theorem reaches; see
# ?leaf_gradient.
#
# ⚠️ Cached at FIRST CALL, not at build time, for the reason
# `.gradient_outputs_idx()` below records -- and here there is a second reason:
# calling into the shared library while the package is still being sourced is not
# something to rely on.
.gradient_output_names <- local({
  nms <- NULL
  function() {
    if (is.null(nms)) {
      nms <<- gradient_output_names()
    }
    nms
  }
})

# ⚠️ ONE call, not four field reads. Every `l$field` is an R6 ACTIVE BINDING -- a
# closure call wrapping a `.Call` -- and this function runs once per perturbation, so
# eleven times per four-parameter gradient. Four reads cost 4.65 us against 0.93 us
# for the one C++ reader that returns all fifteen outputs, and `$` was 12.7% of a
# gradient's self time before this. Same numbers: the reader is the same accessor the
# bindings wrap, and test-gradient.R requires bit-identical gradients.
#
# `.operating_point_names` is the order that reader emits; the four wanted are
# selected by position, matched once and cached.
#
# ⚠️ Cached at FIRST CALL, not at build time. R collates `R/` alphabetically, so this
# file is sourced before `leaf-model.R` and `.operating_point_names` does not exist
# yet -- a build-time `match()` here fails the package load with a message that names
# neither file. Deferring costs one `is.null` per call.
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
# clamp moved the target, because then this is not the evaluation that was asked
# for. Exact equality is the right test: the clamp is a min/max, so an unclamped
# target comes back bit-identical.
.gradient_outputs_at <- function(l, psi, route) {
  route$evaluate(psi)
  if (route$psi_star() != psi) {
    return(NULL)
  }
  .gradient_outputs(l)
}

# The differentiable parameters, in the order C++ indexes them (#4 stage 2).
#
# Derived from leaf_traits() rather than written out, so it cannot drift from the
# trait vector. The C++ side has its own copy, in `phylloptim::gradient::par_names`,
# because it has no way to read this one; `test-gradient-batch.R` compares them
# rather than trusting them, since R passes integer POSITIONS into that
# enumeration -- appending to it is safe and reordering it would silently
# differentiate the wrong parameter.
#
# ⚠️ Computed at FIRST CALL, not at build time, for the reason
# `.gradient_outputs_idx` records: R collates `R/` alphabetically, so this file is
# sourced before `leaf-model.R` and `.leaf_trait_defaults` does not exist yet.
.gradient_par_names <- local({
  nms <- NULL
  function() {
    if (is.null(nms)) {
      nms <<- c(names(.leaf_trait_defaults), .gradient_non_traits)
    }
    nms
  }
})

# The differentiable parameters that are NOT traits, in enumeration order after the
# traits. Named once because three places need to know which slots `leaf_traits()`
# does not supply: the enumeration itself, `.gradient_theta_matrix()`'s trait
# subset, and the batch's column fills. It used to be spelled as "everything but
# the last two", which was correct only while there were exactly two.
.gradient_non_traits <- c("leaf_specific_conductance_max", "resistance",
                          "CF77_lambda_", "TF24_floor_lambda_o")

# What `pars` may name, and the message when it names something else. One
# definition, used by `leaf_gradient()` and by `leaf_gradient_batch()`: the
# multi-layer case has to be named specially, because `resistance` is a real
# parameter on the OTHER supply path and "not differentiable" would be misleading
# rather than merely unhelpful.
# ⚠️ TWO AXES, not one. A parameter can be unavailable because of the SUPPLY PATH
# (`resistance` exists only on the single-potential one) or because of the ACTIVE
# MODEL (`CF77_lambda_` is Cowan-Farquhar's PRESCRIBED price of water and
# `TF24_floor_lambda_o` is the two-term curve's; on every other curve the price is
# emergent rather than set). Both refuse, and both name which axis ruled the
# parameter out, because "not differentiable" alone sends a caller looking in the
# wrong place.
#
# Refusing is deliberate and is NOT the same as returning a zero. A zero means "in
# this objective but inactive at this operating point", which is the `psi_crit`
# case and is information. Structurally absent from the objective is a different
# statement and gets an error.
#
# ⚠️ DERIVED FROM `.gradient_models()`, not written out, and that is the point:
# this used to be an eight-row table that said `"CF77_lambda_"` seven times and
# was a THIRD place a model name had to appear. A curve added to the C++
# enumeration now reaches here on its own.
#
# ⚠️ ONE TABLE, `.gradient_owned_pars`, MAPS PRICE TO OWNER, and this is the only
# copy of it. There are two such parameters now, so "the CF77 slot" has stopped
# being a safe way to describe the class: what the function returns is everything
# owned by a model OTHER than the active one.
.gradient_owned_pars <- c(CF77 = "CF77_lambda_", TF24_floor = "TF24_floor_lambda_o")

.gradient_model_pars <- function(model) {
  if (!(model %in% .gradient_models())) {
    return(NULL)
  }
  unname(.gradient_owned_pars[names(.gradient_owned_pars) != model])
}

.gradient_available_pars <- function(single, model = "collar") {
  nms <- .gradient_par_names()
  if (!single) nms <- setdiff(nms, "resistance")
  excluded <- .gradient_model_pars(model)
  if (is.null(excluded)) {
    stop("no gradient route for model \"", model, "\". Available: ",
         paste(.gradient_models(), collapse = ", "), ".", call. = FALSE)
  }
  setdiff(nms, excluded)
}

.gradient_check_pars <- function(pars, single, model = "collar") {
  unknown <- setdiff(pars, .gradient_available_pars(single, model))
  if (length(unknown)) {
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

# Everything this can differentiate, and its current value, as one named vector.
#
# The fifteen traits, plus the four quantities a calibration fits that are NOT
# traits (#44) and that `pars` therefore used to reject:
#
#   * `leaf_specific_conductance_max` -- a DRIVER, set through set_drivers(). It
#     is also how plant's height reaches the leaf, so the same code path answers
#     plant #537's question about differentiating w.r.t. state.
#   * `resistance` -- the single-potential path's whole soil-to-collar
#     resistance. A DRIVER, like leaf_specific_conductance_max above: it reaches
#     the leaf as the `root_network` argument to set_drivers(), and is read back out
#     of `r_R_V_sum[1]`. It used to be a property of the supply object set through
#     $set_supply_single(), which made it the one differentiable parameter whose
#     setter called setup_clean_leaf() and reset the whole object. Only present on
#     that path: on the multi-layer path the resistances are per-layer vectors
#     rather than one scalar, so there is no such parameter -- which is why it
#     appears here conditionally rather than being rejected later with a worse
#     message.
#
# Nothing in the derivation cares that theta is a trait -- `dpsi*/dtheta = -M/H`
# and `dY/dtheta = dY/dtheta|_psi + (dY/dpsi)(dpsi*/dtheta)` hold for any
# parameter profit depends on. What differs per parameter is only which setter
# applies it, which is .gradient_setter's job.
.gradient_theta <- function(traits, kmax, supply, root_network,
                            CF77_lambda = NA_real_,
                            TF24_floor_lambda_o = NA_real_) {
  theta <- c(unlist(traits), leaf_specific_conductance_max = kmax)
  if (identical(supply$kind, "single")) {
    theta <- c(theta, resistance = root_network$r_R_V_sum[[1]])
  }
  # ⚠️ NA unless the caller is on the route that owns it. C++ reads both slots
  # unconditionally, so each must hold a number the solve will actually use -- and
  # on every other route that number is never read, exactly as `resistance` is not
  # read on the multi-layer path.
  c(theta, CF77_lambda_ = CF77_lambda, TF24_floor_lambda_o = TF24_floor_lambda_o)
}


# Push a parameter vector back onto the leaf, in the one order that is correct.
#
# ⚠️ ORDER IS LOAD-BEARING and the reason this is a function rather than three
# lines at each call site. `set_traits()` returns the leaf to its just-constructed
# state, so the drivers have to be re-supplied AFTER it -- and `set_drivers()` is
# what re-derives vcmax_/jmax_/R_d_ behind the temperature cache `set_traits()`
# has just cleared.
#
# It used to be TWO such setters: `$set_supply_single()` reset the object as well,
# because the single path's resistance was a property of the supply object rather
# than a driver. Now it is a driver on both paths and only `set_traits()` resets.
# The arguments `.resolve_drivers()` takes, in its order, so the drivers list can be
# handed over without naming them at the call site.
.resolve_names <- c("psi_soil", "PPFD", "soil_depth", "root_network",
                    "leaf_specific_conductance_max", "atm_vpd", "ca",
                    "leaf_temp", "atm_o2_kpa", "atm_kpa")

.gradient_setter <- function(l, traits, drivers, supply,
                             fast_stem_curve = TRUE) {
  trait_names <- names(traits)
  single <- identical(supply$kind, "single")
  # ⚠️ RESOLVE THE DRIVERS ONCE. This closure runs once per perturbation -- eleven
  # times for a four-parameter gradient -- and `set_drivers()` re-validates and
  # re-defaults the same values every time. `.resolve_drivers()` is the definition
  # of those rules and is called here exactly once; the loop then applies the result
  # positionally, varying only the two entries a perturbation can move. That is
  # ~12% of a gradient, and it is why the rules are not reimplemented here: a second
  # copy of the 1 m-layer default or the single-path placeholder depth would be free
  # to drift from set_drivers() with nothing to notice.
  base <- do.call(.resolve_drivers, c(list(l), drivers[.resolve_names]))

  # Method closures pulled out of the object once. `$` on an R6 object is a lookup
  # per use and was 12.7% of a gradient's self time; these three run per
  # perturbation.
  apply_traits <- l$set_traits
  apply_drivers <- l$set_physiology

  function(theta, only = NULL) {
    # THE FAST PATH FOR stem_b.
    #
    # stem_b owns the pre-integrated stem vulnerability spline, so moving it
    # through set_traits() rebuilds it -- 11.9 us of incomplete gammas plus two
    # interpolator inits, which is essentially the entire cost of its gradient.
    # But the curve is homogeneous of degree 1 in stem_b, so the spline for
    # another stem_b IS this one with its argument rescaled, and no rebuild is
    # needed at all.
    #
    # ⚠️ stem_c is NOT here, and that is a measured decision rather than an
    # omission: it has no such identity, and reading the curve from its closed
    # form instead -- the obvious alternative -- differentiates a slightly
    # different model and disagrees by 3e-4.
    #
    # ⚠️ Sound only because `only` names a SINGLE parameter, so everything else in
    # `theta` is still what the object was last set to. That is true here and
    # nowhere else, which is why the argument exists rather than the function
    # guessing: .gradient_ift and .gradient_fd perturb exactly one parameter at a
    # time, and the next parameter's first call goes through the full path below,
    # rebuilding the spline on its way.
    if (fast_stem_curve && identical(only, "stem_P50")) {
      l$perturb_stem_P50(theta[["stem_P50"]])
      return(invisible(l))
    }
    # Positional, straight onto the object, rather than through `set_traits()` and
    # `set_drivers()`. Those two rebuild a `leaf_traits` object with
    # `structure(as.list(...))`, re-extract fifteen fields by name, and re-run the
    # driver validation -- 3.0% + 3.3% + 4.9% of a gradient's self time between them,
    # for work whose answer cannot change across perturbations.
    #
    # ⚠️ THIS BYPASSES set_traits()'s INVARIANT CHECKS, and that is sound here for
    # one reason only: `theta` came from `.gradient_theta()`, whose values came from a
    # `leaf_traits()` object the caller already handed in, and the C++ setter asserts
    # the #25 positive-magnitude invariants itself. Do not copy this pattern anywhere
    # the values are not already known-good.
    # ⚠️ POSITIONAL, AND THE COUNT IS NO LONGER WRITTEN DOWN HERE. `trait_names`
    # comes from the C++ enumeration with the non-traits removed, so the argument
    # list is however many traits `set_traits()` takes, in its order, derived
    # rather than restated.
    #
    # It used to be fifteen subscripts spelled out, and the comment here said so:
    # "adding a trait breaks here and nowhere else -- at run time, with `argument
    # <name> is missing` raised inside the generated binding, which names neither
    # this line nor the count." That is how #41 broke, and adding `TF24_floor_a` broke
    # it again in exactly the predicted way -- 118 gradient-batch rows reporting
    # `error` where they had reported `interior`, because the R reference threw and
    # the batch did not. A `do.call` over the derived vector costs no `.Call` (the
    # boundary crossing is `apply_traits` either way, which `test-cost.R` counts)
    # and removes the class rather than the instance.
    tv <- theta[trait_names]
    do.call(apply_traits, unname(as.list(tv)))
    # `resistance` is a driver, so it goes in with the others rather than through
    # $set_supply_single() -- which keeps this to ONE object-resetting call.
    net <- if (single) {
      series_resistance(theta[["resistance"]])
    } else {
      base$root_network
    }
    apply_drivers(net, base$PPFD, base$psi_soil, base$soil_depth,
                  theta[["leaf_specific_conductance_max"]], base$atm_vpd,
                  base$ca, base$leaf_temp, base$atm_o2_kpa, base$atm_kpa)
    # ⚠️ ONLY WHEN IT IS IN PLAY, and the guard is about cost as much as
    # correctness. This is a field write, so it costs a `.Call` per perturbation --
    # and `test-cost.R` counts those. On every route but the one that owns a price
    # the slot is NA, so writing it unconditionally would both poison the leaf and
    # add a crossing to every gradient that has nothing to do with this parameter.
    # The `is.na` is pure R and free.
    lam <- theta[["CF77_lambda_"]]
    if (!is.na(lam)) {
      l$CF77_lambda_ <- lam
    }
    lam_o <- theta[["TF24_floor_lambda_o"]]
    if (!is.na(lam_o)) {
      l$TF24_floor_lambda_o <- lam_o
    }
    invisible(l)
  }
}

# Parameters whose step is RELATIVE rather than floored at 1. See .gradient_step.
.gradient_relative_pars <- c("leaf_specific_conductance_max", "resistance")

# The step: relative to the parameter for values above 1, and plain `step` below
# it. Not a relative step with an epsilon floor -- the floor is at 1, which is
# deliberate. Traits here span `a` = 0.3 to `jmax_25` = 157.44, and a strictly
# relative step would perturb the small ones so little that the difference is
# dominated by the solve's ~1e-09 noise. It also means a trait sitting at zero is
# still perturbed rather than not differentiated at all.
#
# ⚠️ That floor is WRONG for a parameter whose natural magnitude is far below 1,
# and `leaf_specific_conductance_max` is: it defaults to 3.14e-05, so flooring at
# 1 would perturb it by **3%** and measure a secant across a range over which the
# model is visibly nonlinear, not a derivative. Those parameters get a plain
# relative step. For `resistance` (~1e3 and up) the two rules coincide; it is
# listed for the reason rather than for the arithmetic.
.gradient_step <- function(par, value, step) {
  floor <- if (par %in% .gradient_relative_pars) 0 else 1
  max(abs(value), floor) * step
}

# Parameters whose setter takes a SHORTCUT that leaves the rest of the object
# alone, and which therefore require the object to be at the base point before
# they run. See `.gradient_reseat_base()`.
.gradient_shortcut_pars <- function(fast_stem_curve) {
  if (fast_stem_curve) "stem_P50" else character(0)
}

# ⚠️ THE FIX FOR #72, AND THE INVARIANT IT RESTORES: every parameter's gradient
# is taken from the BASE point.
#
# `perturb_stem_P50()` rescales the stem vulnerability spline and touches nothing
# else, which is sound only if everything else on the object is already at base.
# The loops below restore base once at the END, not between parameters, because
# every other parameter's setter goes through the full `set_traits()` +
# `set_physiology()` path and restores it on the way. `stem_P50` on the fast path
# is the one that does not -- so a `stem_b` that is not the first entry of `pars`
# was differentiated at a point displaced by one step in whichever parameter
# preceded it. Measured up to **3.4e-5 relative**, four orders above the ~1e-9
# this function documents as achievable, and on BOTH routes.
#
# Restoring base before the shortcut is the whole fix. It is cheap in the case
# that matters, because `set_traits()` decides the two spline rebuilds by
# comparing the pairs it was given: after a parameter that owns no vulnerability
# curve nothing is rebuilt, so this costs one trait write and one driver write.
# After `stem_c`/`root_b`/`root_c` it does rebuild -- but those are the
# parameters whose own gradients cost a rebuild per side anyway.
#
# ⚠️ Not `stem_b`-specific by name at the call site, deliberately. `root_b` obeys
# the same homogeneity identity and would get the same fast path, at which point
# the condition has to be "this parameter takes a shortcut" rather than a name.
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
  # The mixed partials, kept rather than consumed. `-M/H` is what this function
  # needs, but `M` and `H` are also what a caller integrating its own sensitivity
  # of psi needs (traitecoevo/plant#614), and they are not recoverable from the
  # gradient once divided.
  M <- numeric(length(pars))
  out <- t(vapply(seq_along(pars), function(k) {
    p <- pars[[k]]
    seat(k)
    h <- .gradient_step(p, theta[[p]], step)
    side <- function(sign) {
      th <- theta
      th[[p]] <- th[[p]] + sign * h
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
    # Where the collar was SOLVED for, psi* moves with theta and the implicit
    # function theorem says how. Where it was IMPOSED, it does not move unless
    # the caller says it does -- so the caller supplies dpsi/dtheta, defaulting
    # to zero, and this function has no business deriving one.
    d_psi <- if (is.null(dpsi_dtheta)) -M[[k]] / H else dpsi_dtheta[[k]]
    direct <- (up[-1L] - dn[-1L]) / (2 * h)
    # `collar` is not an output of the evaluation -- it IS psi*, held fixed, so
    # its direct term is zero by construction and the composite reduces to
    # dpsi*/dtheta. Setting it explicitly says so, rather than relying on the
    # difference of two identical numbers.
    g <- direct + dY_dpsi * d_psi
    g[[route$decision]] <- d_psi
    # The envelope theorem, ASSIGNED for the same reason `collar` is: profit's
    # indirect term is identically zero at a stationary point, so stating that is
    # better than multiplying a measured near-zero by dpsi/dtheta. It is also
    # immune to a non-finite `d_psi` -- which a CALLER supplies on the prescribed
    # path -- where `0 * x` would be NaN in this one column while the other four
    # carried +-Inf.
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
  # ⚠️ COERCED, AND THE CALLER MUST USE THE RETURN. `psi` is compared against
  # `opt_root_psi_` with `identical()` to detect the clamp, and `identical(3, 3L)`
  # is FALSE -- so an integer `psi` was reported clamped at a collar it had been
  # given exactly, with the "where it was pulled to" diagnostic showing no
  # movement and nothing to notice. `.gradient_check_psi_batch()` already returned
  # `as.numeric(psi)`, and C++'s `util::identical` is `a == b`, so R alone had it.
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
  # ⚠️ FINITE, AND `anyNA` WAS NOT ENOUGH -- `is.finite` covers NA and NaN,
  # so this is one check rather than two, and it adds the case that matters. An
  # infinite `dpsi_dtheta` does not fail loudly. The composite is
  # `direct + dY_dpsi * dpsi_dtheta`, so four columns come back +-Inf -- and
  # `profit`, which the envelope theorem ASSIGNS from the direct term at a
  # stationary psi, comes back FINITE AND PLAUSIBLE beside them. Measured at
  # `psi = psi*` with `dpsi_dtheta = Inf`: A, gc, psi_stem and collar all Inf,
  # profit 0.0105. Reading `profit` alone is plant's own case (#87), so the one
  # column that survives is the one most likely to be believed. The sibling
  # argument is already checked this way -- `.gradient_check_psi()` demands
  # `is.finite(psi)`.
  if (!is.numeric(dpsi_dtheta) || !all(is.finite(dpsi_dtheta))) {
    stop("`dpsi_dtheta` must be numeric and finite -- no NA, NaN or Inf. It ",
         "is dpsi/dtheta for a collar YOU imposed; a non-finite value returns ",
         "an all-infinite gradient row with a plausible `profit` in it.",
         call. = FALSE)
  }
  # ⚠️ NAMES ARE CHECKED BEFORE RECYCLING, and the order matters. Recycling first
  # overwrote whatever the caller wrote with `pars`, so `dpsi_dtheta =
  # c(stem_b = 1)` against `pars = "vcmax_25"` was silently applied to vcmax_25 --
  # a named argument quietly meaning a different parameter, which is the exact
  # failure naming it was supposed to prevent.
  if (!is.null(names(dpsi_dtheta))) {
    # Matched, not assumed aligned: reordering it silently would be the same
    # class of bug as reordering `pars`.
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
  # ⚠️ THE ROUTE'S STEP, not the caller's, unless the caller asked for one. This arm
  # differences the SOLVE, so its step has to clear the solver's own resolution --
  # see `fd_step`. The composite's step is a different quantity and stays as given.
  if (is.null(attr(step, "user_set"))) {
    step <- route$fd_step
  }
  seat <- .gradient_reseat_base(reset, theta, pars, fast_stem_curve)
  out <- t(vapply(seq_along(pars), function(k) {
    p <- pars[[k]]
    seat(k)
    h <- .gradient_step(p, theta[[p]], step)
    side <- function(sign) {
      th <- theta
      th[[p]] <- th[[p]] + sign * h
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
