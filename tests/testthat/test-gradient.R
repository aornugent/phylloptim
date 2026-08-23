# Trait gradients: set_traits() and leaf_gradient() (issue #4, PLAN 11d stage 1).
#
# Two things are being pinned here, and they are not equally hard.
#
#   1. That the implicit-function composite is right where its premise holds.
#      This is the easy half: the reference values below were arbitrated against
#      a least-squares slope of the solved output over a +-2% trait span at
#      n = 41, and the composite agrees to about four digits.
#   2. That it is NOT used where its premise fails. This is the half that took
#      the work. At a wet-pinned optimum the bare composite returns O(1) where
#      the truth is ~1e-08 -- plausible-looking and wrong by seven orders of
#      magnitude -- so the test that matters is the one asserting the guard
#      fires. `method = "ift"` exists so that failure can be provoked on demand
#      rather than described in a comment.

# `grid_drivers()` and `grid_gradient()` moved to helper-gradient.R when
# test-gradient-batch.R needed the same operating points: each test file gets its
# own environment, so a second copy would have been free to drift and the two
# files would have pinned different points while appearing to pin the same ones.

test_that("set_traits() on a used leaf equals a leaf built with those traits", {
  # The R-side statement of the C++ suite's test_set_traits_matches_a_fresh_leaf,
  # made at the boundary a calibration loop actually crosses. Bit-exact, because
  # the two routes share no code: any derived state set_traits fails to refresh
  # -- the two vulnerability splines, or vcmax_/jmax_/R_d_ behind
  # set_physiology's temperature cache -- shows up as a difference.
  traits <- leaf_traits(vcmax_25 = 110, stem_b = 4.2, root_b = 3.5,
                        cost_scale_TF24 = 8.0)
  d <- grid_drivers(2.0)

  fresh <- leaf_model(traits)
  do.call(set_drivers, c(list(fresh), d))
  fresh$find_root_collar_psi()

  # Solved at the DEFAULTS first, so every cache is warm and pointing at the old
  # traits. Solving first is what makes this a test rather than a coincidence.
  reused <- leaf_model()
  do.call(set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()
  set_traits(reused, traits)
  do.call(set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()

  expect_identical(operating_point(reused), operating_point(fresh))
  expect_identical(reused$vcmax_, fresh$vcmax_)
  expect_identical(reused$R_d_, fresh$R_d_)
})

test_that("set_traits() requires the drivers to be set again", {
  # Not a nicety: the derived photosynthetic parameters really are unknown until
  # set_physiology runs, so the operating point must not be readable as if it
  # were still valid. setup_clean_leaf puts the NA sentinels back.
  l <- leaf_model()
  do.call(set_drivers, c(list(l), grid_drivers(2.0)))
  l$find_root_collar_psi()
  expect_true(is.finite(l$assim_colimited_))

  set_traits(l, leaf_traits(vcmax_25 = 110))
  expect_false(is.finite(l$assim_colimited_))
})

test_that("set_traits() enforces the one representation for psi (#25)", {
  # A bare field write would bypass this, which is one of the four reasons the
  # traits are not bound as settable fields.
  l <- leaf_model()
  for (p in c("psi_crit", "stem_b", "root_b", "root_psi_crit")) {
    tr <- leaf_traits()
    tr[[p]] <- -tr[[p]]
    expect_error(set_traits(l, structure(tr, class = class(leaf_traits()))),
                 "positive magnitude", label = p)
  }
})

test_that("the composite reproduces the arbitrated reference gradients", {
  # psi_soil = 2, PPFD = 900, VPD = 2, one layer, default traits. H = -8.9561.
  # Ratios against a least-squares slope over +-2% at n = 41 were 0.9979-1.0000
  # when these were established; the tolerance below is that agreement, not the
  # composite's own precision, which is finer.
  ref <- rbind(
    vcmax_25        = c(collar =  1.7459e-03, A =  1.7209e-02, profit =  8.2143e-03),
    jmax_25         = c(collar =  1.0862e-04, A =  9.1320e-04, profit =  3.5411e-04),
    cost_scale_TF24 = c(collar = -7.6704e-02, A = -3.9520e-01, profit = -4.1114e-01),
    beta2           = c(collar = -4.2610e-02, A = -2.1954e-01, profit =  1.8272e+00),
    stem_b          = c(collar =  5.5247e-01, A =  2.8465e+00, profit =  2.8396e+00),
    stem_c          = c(collar = -1.6390e-01, A = -8.4448e-01, profit =  4.4233e-01),
    root_b          = c(collar =  2.6656e-04, A =  2.3895e-02, profit =  1.3852e-02))
  # beta_R_H used to be an eighth row here (collar 8.1973e-06, A -2.8802e-04).
  # #33 removed it from the trait vector along with the root-architecture model,
  # so there is no longer a `pars` name for it. The recorded values are kept in
  # this comment rather than deleted, because they are what a hand-rolled
  # two-network difference should reproduce if anyone needs to check that route.

  g <- grid_gradient(2.0, pars = rownames(ref))
  expect_identical(g$status, "interior")
  expect_identical(g$method, "ift")
  expect_equal(g$H, -8.9561, tolerance = 1e-4)

  expect_equal(g$gradient[rownames(ref), "collar"], ref[, "collar"],
               tolerance = 5e-3)
  expect_equal(g$gradient[rownames(ref), "A"], ref[, "A"], tolerance = 5e-3)
  # The profit column, arbitrated the same way. Note the SIGNS differ from A's for
  # beta2 and stem_c: raising either lowers assimilation, and lowers the hydraulic
  # cost by more, so A falls while profit rises. That is the whole reason profit is
  # worth reporting -- it is not a rescaling of dA/dtheta.
  #
  # ⚠️ jmax_25's ratio to the arbitrated slope is 0.989, the worst in this table,
  # and it is the FLOOR rather than the composite: its profit gradient is 3.5e-04
  # against a difference quotient whose noise is the solve's ~1e-09 in profit
  # divided by the step, ~3e-06 here. The other six sit within 6e-04. This is the
  # Precision section's "do not ask for more than about 1e-09", read in a column
  # whose values happen to be small.
  expect_equal(g$gradient[rownames(ref), "profit"], ref[, "profit"],
               tolerance = 5e-3)
})

test_that("the profit gradient is the envelope-theorem partial", {
  # The claim this column rests on. profit is the OBJECTIVE, so
  #   dprofit*/dtheta = (dprofit/dpsi)(dpsi*/dtheta) + dprofit/dtheta|_psi
  # and at an interior optimum the first term is zero. Two consequences are
  # checked, and the second is the one worth having:
  #
  #   1. the composite's profit row equals a difference of the WHOLE SOLVE. A
  #      supplied partial cannot be checked by differencing the step that consumes
  #      it, so the reference has to be the constrained maximum itself.
  #   2. the same number comes out of an evaluation with the collar HELD FIXED at
  #      psi* -- no optimisation, no H, no -M/H. That is the envelope theorem
  #      stated as a computation rather than as algebra.
  #
  # `pars` here excludes root_c and jmax_25 on purpose: both have small profit
  # gradients (8.9e-03 and 3.5e-04) against the ~1e-09/step floor of any
  # difference quotient in profit, so the disagreement they show is the floor and
  # not the identity. Measured at psi_soil = 2.0, step 1e-06: root_c 2.2e-02 and
  # jmax_25 9.4e-03, against 1.2e-02 and 3.0e-03 for the SAME two parameters in
  # the long-standing `A` column at the same point. The floor is the existing one.
  pars <- c("vcmax_25", "stem_b", "stem_c", "cost_scale_TF24", "beta2", "root_b")
  d <- grid_drivers(2.0)

  # 2: the partial at a fixed collar potential, built from scratch rather than
  # through leaf_gradient()'s setter, so this shares no code with the thing it is
  # checking beyond the model itself.
  base <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
  do.call(set_drivers, c(list(base), d))
  base$find_root_collar_psi()
  psi_star <- base$opt_root_psi_
  partial_at_fixed_collar <- function(p, h) {
    at <- function(sgn) {
      tt <- leaf_traits()
      tt[[p]] <- tt[[p]] + sgn * h
      l <- leaf_model(tt, leaf_control(), leaf_supply_multilayer())
      do.call(set_drivers, c(list(l), d))
      l$evaluate_root_collar_psi(psi_star)
      # If the clamp moved the target this is not the evaluation asked for.
      expect_identical(l$opt_root_psi_, psi_star)
      operating_point(l)$profit
    }
    (at(1) - at(-1)) / (2 * h)
  }

  ift <- do.call(leaf_gradient, c(d, list(pars = pars, method = "ift")))
  fd <- do.call(leaf_gradient, c(d, list(pars = pars, method = "fd")))
  expect_identical(ift$status, "interior")

  for (p in pars) {
    h <- max(abs(leaf_traits()[[p]]), 1) * 1e-6
    # Per element, and relative: the whole-column tolerance a matrix comparison
    # gives would let one bad row hide behind stem_b, which is 200x the others.
    expect_equal(ift$gradient[p, "profit"], fd$gradient[p, "profit"],
                 tolerance = 5e-3, label = paste("solve difference,", p))
    expect_equal(ift$gradient[p, "profit"], partial_at_fixed_collar(p, h),
                 tolerance = 5e-3, label = paste("fixed-collar partial,", p))
  }
})

test_that("the profit gradient needs no argmax sensitivity, across the wet grid", {
  # The same identity over the soil-moisture range the interior classification
  # covers, and the number to quote for it. Worst per-element relative
  # disagreement between the composite's profit row and a difference of the whole
  # solve, over eleven parameters: 7.8e-07 at psi_soil 0.5, 2.9e-07 at 1.0 and
  # 1.8e-07 at 3.0. psi_soil = 2.0 is the exception at 2.2e-02, and the test above
  # records why -- two small gradients against the difference quotient's floor,
  # exactly as the `A` column behaves there.
  pars <- c("vcmax_25", "stem_b", "cost_scale_TF24", "beta2", "root_b",
            "leaf_specific_conductance_max")
  for (psi_soil in c(0.5, 1.0, 3.0)) {
    ift <- grid_gradient(psi_soil, pars = pars, method = "ift")
    fd <- grid_gradient(psi_soil, pars = pars, method = "fd")
    expect_identical(ift$status, "interior")
    for (p in pars) {
      expect_equal(ift$gradient[p, "profit"], fd$gradient[p, "profit"],
                   tolerance = 1e-5,
                   label = paste("psi_soil", psi_soil, p))
    }
  }
})

test_that("the composite and the finite difference agree at interior points", {
  # The cross-check that needs no external reference: two routes to the same
  # derivative, sharing only the model. They are checked across soil moisture
  # because the indirect term's share of the answer varies with it.
  pars <- c("vcmax_25", "jmax_25", "stem_b", "cost_scale_TF24", "root_b")
  for (psi_soil in c(0.5, 1.0, 2.0, 3.0)) {
    ift <- grid_gradient(psi_soil, pars = pars, method = "ift")
    fd <- grid_gradient(psi_soil, pars = pars, method = "fd")
    expect_identical(ift$status, "interior", label = paste("psi_soil", psi_soil))
    # 1e-3 relative, which is the finite difference's precision rather than the
    # composite's: differencing the solved argmax carries the solve's own ~1e-09
    # floor divided by the step. Do not read agreement here as either route
    # being good to more than three digits.
    expect_equal(ift$gradient, fd$gradient, tolerance = 1e-3,
                 label = paste("psi_soil", psi_soil))
  }
})

test_that("dcollar/dtheta is dpsi*/dtheta, computed two ways", {
  # The `collar` column is psi* itself, so the composite reports -M/H there while
  # the fallback differences the solved argmax. They are different computations of
  # the same quantity, which is why the column is worth reporting at all.
  pars <- c("vcmax_25", "stem_b", "beta2")
  ift <- grid_gradient(2.0, pars = pars, method = "ift")
  fd <- grid_gradient(2.0, pars = pars, method = "fd")
  expect_equal(ift$gradient[, "collar"], fd$gradient[, "collar"],
               tolerance = 1e-3)
})

test_that("a pinned optimum takes the fallback, and the composite would be wrong", {
  # The hazard, pinned to specific operating points. These three are rows of the
  # golden grid; the wet-pinned pair are where the bare composite is worst.
  #
  # ⚠️ The truth at a wet-pinned point is ~1e-07 or smaller -- the leaf is stuck
  # at the wettest feasible collar and barely responds -- so the ASSERTION IS
  # ABSOLUTE, not relative. A ratio test here would be meaningless in both
  # directions: the fallback's relative error is large and does not matter, and
  # the composite's O(1) answer is catastrophic precisely because the truth is
  # not O(1).
  wet <- list(list(psi_soil = 4, vpd = 2.0, layers = 5L),
              list(psi_soil = 3, vpd = 4.0, layers = 3L))
  for (w in wet) {
    lab <- sprintf("psi_soil=%g vpd=%g layers=%d", w$psi_soil, w$vpd, w$layers)
    auto <- grid_gradient(w$psi_soil, vpd = w$vpd, layers = w$layers,
                          pars = c("stem_b", "vcmax_25"))
    expect_identical(auto$status, "pinned", label = lab)
    expect_identical(auto$method, "fd", label = lab)
    # The guard fired with room to spare: measured stationarity at these points
    # is ~6e-06 to ~9e-06 against a 1e-08 threshold.
    expect_gt(auto$stationarity, 1e-7)
    expect_lt(abs(auto$gradient["stem_b", "A"]), 1e-3)

    # Forcing the composite here does not produce the wrong number -- it fails, on
    # the stationarity the composite rests on. Asserted by the refusal rather than
    # by which guard raised it: psi* also sits 1e-06 of a bracket width from its
    # bound here, so the centred step cannot fit either, and pinning that message
    # would pin the narrower of two reasons.
    expect_error(grid_gradient(w$psi_soil, vpd = w$vpd, layers = w$layers,
                               pars = "stem_b", method = "ift"),
                 "pinned operating point", label = lab)
  }

  # Pinned DRY is the milder case and is worth separating: the composite is only
  # a factor of 1.6 out for vcmax_25 there rather than seven orders, so a test
  # that only looked at the wet rows would let a threshold drift past it.
  dry <- grid_gradient(4, vpd = 0.5, layers = 3L,
                       pars = c("stem_b", "vcmax_25", "R_d_25"))
  expect_identical(dry$status, "pinned")
  expect_equal(dry$gradient["vcmax_25", "A"], 0.017548, tolerance = 1e-3)
  expect_equal(dry$gradient["stem_b", "A"], 4.6991, tolerance = 1e-3)

  # `dA/dvcmax_25` is a PARTIAL at fixed respiration, since `R_d_25` is its own
  # trait. A fit that wants respiration to follow Vcmax adds the two columns, so
  # that combination is asserted rather than left implicit.
  expect_equal(dry$gradient["vcmax_25", "A"] +
                 0.015 * dry$gradient["R_d_25", "A"],
               0.014357, tolerance = 1e-3)
})

test_that("at a pinned optimum the BOUND's own trait carries the gradient", {
  # The sharpest statement of why the fallback is not merely a fallback, and it
  # is a stronger claim than "the composite is inaccurate there".
  #
  # `psi_crit` does not appear in the profit function at all -- it only sets the
  # dry end of the feasible collar interval. So at an INTERIOR optimum it has no
  # effect and its gradient is exactly zero, which is right. At a DRY-PINNED
  # optimum it IS the binding constraint, so moving it moves the answer directly.
  #
  # ⚠️ The composite cannot see that, and not by a small margin: psi* enters
  # `dprofit` nowhere, so M = 0, dpsi*/dtheta = 0, the direct term is 0, and the
  # composite returns EXACTLY ZERO where the truth is 1.26. A silent zero is
  # arguably a worse failure than the wet-pinned 3.5e+07 -- it tells an optimiser
  # the parameter does nothing, and nothing about it looks wrong.
  interior <- grid_gradient(2.0, pars = c("psi_crit", "root_psi_crit"))
  expect_identical(interior$status, "interior")
  expect_equal(interior$gradient["psi_crit", "A"], 0)
  expect_equal(interior$gradient["root_psi_crit", "A"], 0)
  # Including `profit`, and there the zero is worth saying out loud rather than
  # merely tolerating: psi_crit is absent from the profit function, so at an
  # interior optimum the envelope-theorem partial is exactly zero because the
  # constraint is not binding. It is a result, not a row that failed to be filled.
  expect_equal(interior$gradient["psi_crit", "profit"], 0)
  expect_equal(interior$gradient["root_psi_crit", "profit"], 0)

  pinned <- grid_gradient(4, vpd = 0.5, layers = 3L, pars = "psi_crit")
  expect_identical(pinned$status, "pinned")
  expect_identical(pinned$method, "fd")
  # Arbitrated against a least-squares slope of the solved A over +-2% of
  # psi_crit at n = 41, which gives 1.2605 (R^2 = 0.998; the fit is imperfect
  # because the response has a kink in it, which is the whole point).
  expect_equal(pinned$gradient["psi_crit", "A"], 1.2605, tolerance = 5e-3)
  expect_gt(pinned$gradient["psi_crit", "collar"], 0)
  # And where the constraint IS binding the fallback recovers it in profit as well.
  # Same arbitration, over the solved profit: 0.5099, R^2 = 0.996 -- the same kink.
  expect_equal(pinned$gradient["psi_crit", "profit"], 0.5099, tolerance = 1e-2)
})

test_that("a shut-down operating point reports no gradient and still differences", {
  # psi_soil = 6 is drier than psi_crit, so there is no optimisation to
  # differentiate: dprofit returns a sentinel zero rather than a derivative, and
  # H comes back zero with it. The composite has nothing to stand on and says so
  # -- but the gradient itself is not all zero, because a shut-down leaf still
  # respires, so the fallback still has work to do.
  #
  # WHICH parameter carries it is the sharpest statement in the suite: `A = -R_d`
  # exactly here, so `dA/dR_d_25` is -1 and `dA/dvcmax_25` is EXACTLY zero --
  # vcmax_25 does not reach A at all at a shut-down point, so both perturbed solves
  # return the same bits. The -1 is a central difference and lands within ~6e-11.
  g <- grid_gradient(6.0, pars = c("vcmax_25", "stem_b", "R_d_25"))
  expect_identical(g$status, "no-gradient")
  expect_identical(g$method, "fd")
  expect_equal(g$value[["A"]], -leaf_traits()$R_d_25)
  expect_equal(g$gradient["R_d_25", "A"], -1, tolerance = 1e-8)
  expect_identical(g$gradient["vcmax_25", "A"], 0)
  expect_equal(g$gradient["stem_b", "A"], 0)
  expect_equal(g$gradient["cost_scale_TF24", "A"], 0)

  for (p in names(analytic)) {
    # 1e-03 relative, and the loose tolerance is the DIFFERENCE QUOTIENT's floor
    # rather than the closed form's: the solved profit carries ~1e-09 and the
    # fallback divides a difference of two of them by a ~1e-06 step, which bounds
    # the error at ~1e-04 relative here. Worst OBSERVED is 5.1e-06 (stem_b) with
    # the other five at 1e-10 to 5e-10 -- but that is one draw from a distribution
    # the floor bounds at 1e-04, so the tolerance is set from the bound. Pinning
    # 5.1e-06 would be pinning the luck.
    expect_equal(g$gradient[[p, "profit"]], analytic[[p]], tolerance = 1e-3,
                 label = paste("shut-down dprofit/d", p, sep = ""))
    # Each is a real response, not a small number that happens to pass above.
    expect_gt(abs(analytic[[p]]), 1e-2)
  }
  expect_equal(g$gradient["root_b", "profit"], 0)
  # vcmax_25 is the one parameter reaching both columns, and the cost term does not
  # depend on it, so the two differences are differences of the same function offset
  # by a constant. NOT bit-identical, and the reason is worth one line: profit is
  # -8.39 where A is -1.44, so differencing profit cancels against an operand six
  # times larger and loses ~6e-10 relative. Measured 9e-11.
  expect_equal(g$gradient[["vcmax_25", "profit"]],
               g$gradient[["vcmax_25", "A"]], tolerance = 1e-9)

  # Forcing the composite here is an error rather than a wrong number: unlike a
  # pinned point, there is no curvature to divide by at all.
  expect_error(grid_gradient(6.0, pars = "vcmax_25", method = "ift"),
               "usable curvature")
})

test_that("the active-set classification splits the golden grid as measured", {
  # The whole 288-point grid, which is what makes `stationarity_tol` a measured
  # threshold rather than a guess. Two claims:
  #
  #   * the counts: 198 interior, 42 pinned, 48 with no gradient at all;
  #   * the GAP. The worst interior point is ~5e-11 and the mildest pinned one
  #     ~1.2e-02, so the default 1e-08 sits in an empty band nine orders wide.
  #     That is the number worth defending, because a threshold in a populated
  #     region would be a judgement call instead.
  #
  #     ⚠️ The pinned edge read ~6e-06 until the curvature was guarded against
  #     dprofit's shut-down sentinel. On 20 of these 42 pinned rows one arm of
  #     the central difference was the sentinel, which inflated |H| and so shrank
  #     |resid/H|. The lower bound below is set an order under the guarded figure
  #     so that reinstating the sentinel fails here rather than passing quietly.
  psi_soils <- c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0)
  vpds <- c(0.5, 1.0, 2.0, 4.0)
  layer_counts <- c(1L, 3L, 5L)
  ppfds <- c(100, 500, 900, 1500)

  rows <- expand.grid(layers = layer_counts, vpd = vpds, ppfd = ppfds,
                      psi_soil = psi_soils)
  res <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    # One trait only: the classification is a property of the operating point,
    # not of which trait is being differentiated.
    g <- grid_gradient(r$psi_soil, ppfd = r$ppfd, vpd = r$vpd,
                       layers = r$layers, pars = "vcmax_25")
    # Does forcing the composite here produce a number or refuse? Recorded per
    # row so the claim below is about the whole grid rather than a sample.
    forced <- tryCatch({
      grid_gradient(r$psi_soil, ppfd = r$ppfd, vpd = r$vpd, layers = r$layers,
                    pars = "vcmax_25", method = "ift")
      TRUE
    }, error = function(e) FALSE)
    data.frame(psi_soil = r$psi_soil, status = g$status,
               stationarity = g$stationarity, forced_ift_ran = forced)
  }))

  expect_identical(nrow(res), 288L)
  expect_identical(sum(res$status == "interior"), 198L)
  expect_identical(sum(res$status == "pinned"), 42L)
  expect_identical(sum(res$status == "no-gradient"), 48L)

  # Pinned points are all at the dry end, which is where a drought calibration
  # wanders and is the reason this is not a curiosity.
  expect_setequal(unique(res$psi_soil[res$status == "pinned"]), c(3, 4))

  worst_interior <- max(res$stationarity[res$status == "interior"])
  mildest_pinned <- min(res$stationarity[res$status == "pinned"])
  expect_lt(worst_interior, 1e-9)
  expect_gt(mildest_pinned, 1e-3)

  # The second guard, over the whole grid: the composite runs on every interior
  # row and refuses on every other one. So a caller who forces `ift` cannot
  # obtain the silently-wrong answer -- they get an error. This is a measured
  # property of where these points sit, not a designed one, and it is asserted
  # here so that it is noticed if it stops being true.
  expect_true(all(res$forced_ift_ran[res$status == "interior"]))
  expect_false(any(res$forced_ift_ran[res$status != "interior"]))
})

test_that("leaf_gradient() works on the single-potential supply path", {
  g <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                     supply = leaf_supply_single(),
                     root_network = series_resistance(1e3),
                     pars = c("vcmax_25", "stem_b"))
  expect_identical(g$status, "interior")
  fd <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                      supply = leaf_supply_single(),
                      root_network = series_resistance(1e3),
                      pars = c("vcmax_25", "stem_b"), method = "fd")
  expect_equal(g$gradient, fd$gradient, tolerance = 1e-3)
})

test_that("leaf_gradient() leaves the traits it was given", {
  # The perturbation loop mutates one Leaf, so it has to put it back. The Leaf is
  # internal, so what this really checks is that two calls agree -- a loop that
  # leaked its last perturbation would drift.
  first <- grid_gradient(2.0, pars = "vcmax_25")
  second <- grid_gradient(2.0, pars = "vcmax_25")
  expect_identical(first$gradient, second$gradient)
  expect_identical(first$value, second$value)
})

test_that("GSS_tol_abs does not reach the production path", {
  # Not a gradient test, but it belongs with the ones above: it is the claim that
  # makes them possible. Since PLAN 11a the collar solve solves dprofit == 0 at
  # 1e-12 instead of searching profit to GSS_tol_abs, so this control has no
  # effect on the answer -- and an argmax determined only to 1e-3 could not be
  # differentiated at all.
  #
  # Asserted because it was documented WRONG in two places until this branch: the
  # vignette demonstrated 1e-3 against 1e-5 and said they "differ in the fourth
  # decimal", and leaf_control()'s roxygen said this sets how well the operating
  # point is determined. Both were true before 11a. Bit-exact, because "no
  # effect" is exactly what is being claimed.
  loose <- leaf_solve(psi_soil = 2.0, PPFD = 900)
  tight <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                      control = leaf_control(GSS_tol_abs = 1e-5))
  expect_identical(loose, tight)
})

test_that("the documented examples classify the way their comments claim", {
  # Written after the first version of ?leaf_gradient offered
  # `psi_soil = c(4, 4.25, 4.5), atm_vpd = 2` as "a pinned optimum" and it came
  # back interior. An example whose comment contradicts its output is worse than
  # no example, and nothing else here would have caught it: examples are run by
  # R CMD check but their *classifications* are not asserted anywhere. So assert
  # them, using exactly the calls in the roxygen block.
  expect_identical(leaf_gradient(psi_soil = 2.0, PPFD = 900)$method, "ift")

  dry <- leaf_gradient(psi_soil = 4.5, PPFD = 900, atm_vpd = 3.0)
  expect_identical(dry$method, "fd")
  expect_identical(dry$status, "pinned")
})

test_that("the two non-trait parameters agree with a resolved reference", {
  # #44. Half of leaf-calibration's free parameters are not traits: `K_total` and
  # `f_plant` reach the leaf as leaf_specific_conductance_max (a driver) and the
  # single-potential resistance. Nothing in the derivation cares -- the implicit
  # function theorem is applied to dprofit/dpsi = 0, and any parameter profit
  # depends on goes through it identically -- so the thing to check is that the
  # plumbing applies the perturbation through the right setter.
  #
  # Arbitrated the way 11c arbitrated the traits: against a central difference of
  # the WHOLE SOLVE, read where the step sweep plateaus rather than at one step.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_single(),
            root_network = series_resistance(1e4))

  solve_at <- function(par, value) {
    a <- d
    # `resistance` is a driver now, so perturbing it means rebuilding the
    # root_network driver rather than rebuilding the supply object -- which is the
    # whole of the consistency change, seen from the caller's side.
    if (identical(par, "resistance")) {
      a$root_network <- series_resistance(value)
    } else {
      a[[par]] <- value
    }
    x <- do.call(leaf_solve, a)
    c(A = x$A, gc = x$gc, psi_stem = x$psi_stem, collar = x$collar,
      profit = x$profit)
  }

  g <- do.call(leaf_gradient,
               c(d, list(pars = c("leaf_specific_conductance_max",
                                  "resistance"), method = "ift")))
  expect_identical(g$status, "interior")

  for (p in c("leaf_specific_conductance_max", "resistance")) {
    v <- if (identical(p, "resistance")) 1e4 else 3.14e-5
    h <- v * 1e-5
    ref <- (solve_at(p, v + h) - solve_at(p, v - h)) / (2 * h)
    # Six digits, which is well inside the reference's own precision. Measured
    # agreement is 8 digits for the driver and 9 for the resistance.
    expect_equal(g$gradient[p, ], ref, tolerance = 1e-6)
    # And not accidentally zero, which is how this would fail if the setter did
    # nothing at all -- the failure mode psi_crit genuinely has (see above).
    expect_gt(abs(g$gradient[p, "A"]), 1e-8)
  }
})

test_that("the non-trait gradients are stable across decades of step", {
  # The mixed partial is a difference of an exact quantity, so it should be
  # step-insensitive over a wide band; a parameter whose step rule was wrong
  # would show up here as drift rather than as a wrong answer.
  #
  # ⚠️ This is the test that would have caught the step-rule trap:
  # leaf_specific_conductance_max defaults to 3.14e-05, so the traits' rule of
  # flooring the step at 1 would perturb it by 3% and measure a secant.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_single(),
            root_network = series_resistance(1e4))
  for (p in c("leaf_specific_conductance_max", "resistance")) {
    vals <- vapply(10^-(4:7), function(st) {
      do.call(leaf_gradient, c(d, list(pars = p, step = st,
                                       method = "ift")))$gradient[p, "A"]
    }, numeric(1))
    expect_lt(diff(range(vals)) / abs(mean(vals)), 1e-6)
  }
})

test_that("the active-set guard covers the non-trait parameters too", {
  # Both move the FEASIBLE COLLAR INTERVAL and not only the profit function --
  # resistance directly, the conductance through the supply curve -- so a
  # perturbation can push a bound past psi* more readily than a trait
  # perturbation can. The guard must not be trait-specific.
  dry <- list(psi_soil = 5.5, PPFD = 900, atm_vpd = 3.0,
              supply = leaf_supply_single(),
              root_network = series_resistance(1e4))
  auto <- do.call(leaf_gradient,
                  c(dry, list(pars = c("leaf_specific_conductance_max",
                                       "resistance"))))
  expect_true(auto$status %in% c("pinned", "no-gradient"))
  expect_identical(auto$method, "fd")
  expect_true(all(is.finite(auto$gradient)))

  # ...and forcing the composite there fails loudly rather than returning a
  # confident number, which is the whole reason `auto` exists.
  expect_error(do.call(leaf_gradient,
                       c(dry, list(pars = "resistance", method = "ift"))))
})

test_that("perturbing stem_b by rescaling equals perturbing it by a rebuild", {
  # PLAN 11f. The stem cumulative-vulnerability integral is homogeneous of degree
  # 1 in stem_b -- G(psi; s*b, c) = s*G(psi/s; b, c) -- so the spline for a
  # perturbed stem_b is the existing one with its argument rescaled, and the
  # 11.9 us of incomplete gammas a rebuild spends is unnecessary. Worth 24x on
  # that parameter's gradient in C++.
  #
  # ⚠️ What has to be true is not the formula but that the two routes give the
  # SAME GRADIENT, so that is what is asserted, at four operating points and on
  # both routes through the composite. `fast_stem_curve = FALSE` exists for this.
  for (psi_soil in c(1.0, 1.5, 2.0, 3.0)) {
    for (m in c("ift", "fd")) {
      d <- list(psi_soil = psi_soil, PPFD = 900, atm_vpd = 2.0)
      fast <- do.call(leaf_gradient, c(d, list(pars = "stem_b", method = m,
                                               fast_stem_curve = TRUE)))
      slow <- do.call(leaf_gradient, c(d, list(pars = "stem_b", method = m,
                                               fast_stem_curve = FALSE)))
      # 1e-3 rather than tighter, and the reason is worth knowing: a rebuild
      # reseeds 101 incomplete gammas per side, so the two sides of the central
      # difference carry UNCORRELATED rounding, which dividing by h ~ 4e-6
      # amplifies. Measured, the rebuild route jitters by up to 9e-05 between
      # neighbouring steps where the rescale route does not -- so the loose
      # tolerance here is the rebuild's noise, not the rescale's error.
      expect_equal(fast$gradient, slow$gradient, tolerance = 1e-3,
                   info = paste(psi_soil, m))
      expect_identical(fast$status, slow$status)
    }
  }
})

test_that("a gradient leaves the leaf's stem curve as it found it", {
  # ⚠️ perturb_stem_b() leaves the splines built at a different stem_b on
  # purpose, so the thing that could go wrong is a gradient call that does not
  # put it back -- after which every later solve is quietly on a rescaled curve.
  # set_traits() forces the rebuild for exactly this reason, and BIT-identity is
  # the right test: a restored leaf shares no code path with a fresh one.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0)
  before <- do.call(leaf_solve, d)

  invisible(do.call(leaf_gradient, c(d, list(pars = "stem_b"))))
  invisible(do.call(leaf_gradient, c(d, list(pars = c("stem_b", "vcmax_25")))))
  # ...and one that ends on the fast path, which is the ordering that would leave
  # it displaced if the restore at the end of the loop were missing.
  invisible(do.call(leaf_gradient, c(d, list(pars = c("vcmax_25", "stem_b")))))

  expect_identical(do.call(leaf_solve, d), before)
})

test_that("leaf_gradient() rejects bad arguments", {
  expect_error(grid_gradient(2.0, pars = "not_a_trait"),
               "cannot differentiate")
  # `resistance` is a real parameter -- on the other supply path -- so the
  # message says which path it belongs to rather than calling it unknown.
  expect_error(grid_gradient(2.0, pars = "resistance"),
               "single-potential path only")
  expect_error(grid_gradient(2.0, step = -1), "positive number")
  expect_error(grid_gradient(2.0, traits = list(vcmax_25 = 96)),
               "leaf_traits")
  expect_error(grid_gradient(2.0, method = "magic"))
})

test_that("the gradient is reported for every output the fit needs", {
  # leaf-calibration fits three responses -- A, gs and psi_leaf -- so all three
  # are differentiated, not just A. `collar` comes along because it is psi*, and
  # `profit` because it is what the leaf maximises and what feeds net production.
  # The order is appended, not inserted, so a caller selecting by position does not
  # move.
  g <- grid_gradient(2.0, pars = "vcmax_25")
  expect_identical(colnames(g$gradient),
                   c("A", "gc", "psi_stem", "collar", "profit"))
  expect_identical(names(g$value),
                   c("A", "gc", "psi_stem", "collar", "profit"))
  expect_true(all(is.finite(g$gradient)))
  # Raising vcmax_25 raises assimilation and opens the stomata, and the leaf pays
  # for it with a more negative water potential (a larger positive magnitude). The
  # carbon it gains exceeds the water it spends, or the optimum would not have
  # moved that way, so profit rises too -- by less than A.
  expect_gt(g$gradient["vcmax_25", "A"], 0)
  expect_gt(g$gradient["vcmax_25", "gc"], 0)
  expect_gt(g$gradient["vcmax_25", "psi_stem"], 0)
  expect_gt(g$gradient["vcmax_25", "profit"], 0)
  expect_lt(g$gradient["vcmax_25", "profit"], g$gradient["vcmax_25", "A"])
  # `value` carries the solved profit alongside the gradient, so a fit reading
  # profit as a response needs no second call.
  expect_identical(g$value[["profit"]],
                   do.call(leaf_solve, grid_drivers(2.0))$profit)
})

# ---------------------------------------------------------------------------
# Reusing a Leaf across gradients (#52)
# ---------------------------------------------------------------------------

test_that("leaf_gradient(x =) matches building a leaf per call", {
  # The whole claim: reuse is a pure cost saving. Bit-identical, not merely close --
  # the setter applies the same traits and drivers either way, so anything else means
  # the reused leaf carried state into the answer.
  # ⚠️ stem_b is raised, not lowered. #38: the vulnerability spline is built over
  # [0, b*log(100)^(1/c)] and never consults psi_crit, so shrinking stem_b below
  # ~3.4 at the default psi_crit = 5.87 puts the solve outside its own domain and it
  # dies naming neither. stem_b = 3.1 does exactly that; 4.2 gives a domain of 7.4.
  tr <- leaf_traits(vcmax_25 = 105, stem_b = 4.2)
  args <- list(psi_soil = 2.0, PPFD = 900, atm_vpd = 1.5, traits = tr,
               pars = c("vcmax_25", "stem_b", "cost_scale_TF24"))
  fresh <- do.call(leaf_gradient, args)

  l <- leaf_model(traits = tr)
  reused <- do.call(leaf_gradient, c(args, list(x = l)))
  expect_identical(reused$gradient, fresh$gradient)
  expect_identical(reused$value, fresh$value)
  expect_identical(reused$status, fresh$status)
  expect_identical(reused$method, fresh$method)
  expect_identical(reused$H, fresh$H)

  # And again on the SAME object, which is the point of it: a second call must not
  # inherit the first call's last perturbation.
  again <- do.call(leaf_gradient, c(args, list(x = l)))
  expect_identical(again$gradient, fresh$gradient)
})

test_that("a reused leaf is differentiated at `traits`, not at its own state", {
  # `x` is a vessel. A leaf built at OTHER traits must give the same answer as a
  # fresh one, because the setter applies `traits` before anything is read. If it
  # did not, `value` and `psi_star` -- and so the whole composite -- would describe
  # the wrong point, plausibly and with no symptom. This is the failure mode that
  # made placing the leaf with reset() rather than set_drivers() load-bearing.
  tr <- leaf_traits(vcmax_25 = 105, stem_b = 4.2)
  args <- list(psi_soil = 2.0, PPFD = 900, traits = tr, pars = "vcmax_25")
  fresh <- do.call(leaf_gradient, args)

  # Deliberately different, and deliberately still inside the #38 domain -- the
  # vessel's traits are overwritten before any solve, but leaf_model() does build
  # its splines, and an invalid vessel would confuse the diagnosis if this failed.
  wrong <- leaf_model(traits = leaf_traits(vcmax_25 = 60, stem_b = 4.6))
  reused <- do.call(leaf_gradient, c(args, list(x = wrong)))
  expect_identical(reused$gradient, fresh$gradient)
  expect_identical(reused$value, fresh$value)
})

test_that("leaf_gradient(x =) leaves the caller's leaf solved at the base point", {
  # hazard 8: an output no path rewrites reads as the previous caller's. The setter
  # leaves the leaf at the last perturbation, so this function restores it.
  tr <- leaf_traits()
  l <- leaf_model(traits = tr)
  set_drivers(l, psi_soil = 2.0, PPFD = 900)
  l$find_root_collar_psi()
  before <- operating_point(l)

  leaf_gradient(psi_soil = 2.0, PPFD = 900, x = l, traits = tr, pars = "vcmax_25")
  expect_identical(operating_point(l), before)

  # Restored on an ERROR path too -- there are several stop()s downstream of the
  # first perturbation, and a half-perturbed leaf handed back would be worse than
  # the error.
  expect_error(
    leaf_gradient(psi_soil = 2.0, PPFD = 900, x = l, traits = tr,
                  pars = "not_a_parameter"),
    "cannot differentiate")
  expect_identical(operating_point(l), before)
})

test_that("leaf_gradient() refuses the combinations that would disagree", {
  tr <- leaf_traits()
  l <- leaf_model(traits = tr)
  expect_error(leaf_gradient(psi_soil = 2, x = l), "`traits` must be given with")
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr,
                             control = leaf_control()), "`control` comes from")
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr,
                             supply = leaf_supply_multilayer()),
               "`supply` comes from")
  expect_error(leaf_gradient(psi_soil = 2, x = leaf_traits(), traits = tr),
               "must be a Leaf")

  # The supply path is read off the object, so `resistance` is differentiable when
  # `x` is on the single-potential path and not when it is not -- without the caller
  # naming the path twice.
  s <- leaf_model(supply = leaf_supply_single())
  g <- leaf_gradient(psi_soil = 1.5, x = s, traits = tr,
                     root_network = series_resistance(1e3), pars = "resistance")
  expect_true(is.finite(g$gradient["resistance", "A"]))
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr, pars = "resistance"),
               "single-potential path only")
})

test_that("the setter's positional trait call cannot drift in arity", {
  # `.gradient_setter()` applies traits POSITIONALLY, straight onto the object, to
  # skip rebuilding a leaf_traits per perturbation. That is a hard-coded FOURTEEN. If
  # a trait is added to leaf_traits() and to the C++ setter, nothing about that call
  # fails to compile -- it would silently pass the wrong value for every argument
  # after the new one. So the arity is asserted here rather than trusted.
  #
  expect_length(leaf_traits(), 14L)
  expect_length(formals(leaf_model()$set_traits), 14L)
  expect_identical(names(leaf_traits()), names(formals(leaf_model()$set_traits)))
})
