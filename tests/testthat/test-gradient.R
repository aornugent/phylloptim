# Trait gradients: set_traits() and leaf_gradient() (#4).
#
# Two things are being pinned here, and they are not equally hard.
#
#   1. That the implicit-function composite is right where its premise holds.
#      This is the easy half: the reference values below were arbitrated against
#      a least-squares slope of the solved output over a +-2% trait span at
#      n = 41, and the composite agrees to about four digits.
#   2. That it is NOT used where its premise fails. This is the half that took
#      the work. At a boundary-soil optimum the bare composite returns O(1) where
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
  traits <- leaf_traits(vcmax_25 = 110, stem_P50 = 4.2, root_P50 = 3.5,
                        TF24_cost_scale = 8.0)
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
  for (p in c("stem_P50", "root_P50")) {
    tr <- leaf_traits()
    tr[[p]] <- -tr[[p]]
    expect_error(set_traits(l, structure(tr, class = class(leaf_traits()))),
                 "positive magnitude", label = p)
  }
})

test_that("the composite reproduces the arbitrated reference gradients", {
  # psi_soil = 2, PPFD = 900, VPD = 2, one layer, default traits.
  # Ratios against a least-squares slope over +-2% at n = 41 were 0.9979-1.0000
  # when these were established; the tolerance below is that agreement, not the
  # composite's own precision, which is finer.
  # ⚠️ THE CURVE ROWS ARE NOT RE-ARBITRATED NUMBERS -- they are the arbitrated
  # b-gradients carried through the reparameterisation by its own chain rule, so
  # this table asserts the chain rule as well as the gradient. With
  # b = P50 / (ln 2)^(1/c),
  #
  #     d/dP50 |_c = (ln 2)^(-1/c) . d/db
  #     d/dc |_P50 = d/dc |_b  +  (db/dc) . d/db,   db/dc = b.ln(ln 2)/c^2
  #
  # and the derived psi_crit contributes NOTHING at this row because it is an
  # INTERIOR optimum, where dY/dpsi_crit is exactly zero. That is why the identity
  # is clean here and is not at a pinned optimum -- see the pinned test below,
  # which asserts the other half.
  tr <- leaf_traits()
  stem_b <- tr$stem_P50 / log(2)^(1 / tr$stem_c)
  d_dP50 <- log(2)^(-1 / tr$stem_c)
  d_dc   <- stem_b * log(log(2)) / tr$stem_c^2
  b_stem <- c(collar = 5.5247e-01, A = 2.8465e+00)   # arbitrated d/dstem_b
  b_root <- c(collar = 2.6656e-04, A = 2.3895e-02)   # arbitrated d/droot_b
  ref <- rbind(
    vcmax_25        = c(collar =  1.7459e-03, A =  1.7209e-02),
    jmax_25         = c(collar =  1.0862e-04, A =  9.1320e-04),
    TF24_cost_scale = c(collar = -7.6704e-02, A = -3.9520e-01),
    TF24_beta2           = c(collar = -4.2610e-02, A = -2.1954e-01),
    stem_P50        = d_dP50 * b_stem,
    stem_c          = c(collar = -1.6390e-01, A = -8.4448e-01) + d_dc * b_stem,
    root_P50        = (log(2)^(-1 / tr$root_c)) * b_root)
  # beta_R_H used to be an eighth row here (collar 8.1973e-06, A -2.8802e-04).
  # #33 removed it from the trait vector along with the root-architecture model,
  # so there is no longer a `pars` name for it. The recorded values are kept in
  # this comment rather than deleted, because they are what a hand-rolled
  # two-network difference should reproduce if anyone needs to check that route.

  g <- grid_gradient(2.0, pars = rownames(ref))
  expect_identical(g$status, "interior")
  expect_identical(g$method, "ift")
  # ⚠️ THE MODEL'S OWN NUMBER, NOT AN ARBITRATED ONE, and it had drifted from the
  # model by more than this tolerance: the literal read -8.9578 and the comment
  # above it -8.9561, where the model gives -8.9554 at 100 knots. It is pinned at
  # the CONVERGED value, which holds at every resolution measured -- -8.955414,
  # -8.955057, -8.955110, -8.955133, -8.955141 at 100, 200, 400, 800 and 1600
  # knots, a spread of 4e-05 relative and so well inside this tolerance. The
  # seven rows below are the arbitrated ones and did not move.
  expect_equal(g$H, -8.9551, tolerance = 1e-4)

  expect_equal(g$gradient[rownames(ref), "collar"], ref[, "collar"],
               tolerance = 5e-3)
  expect_equal(g$gradient[rownames(ref), "A"], ref[, "A"], tolerance = 5e-3)
})

test_that("the composite and the finite difference agree at interior points", {
  # The cross-check that needs no external reference: two routes to the same
  # derivative, sharing only the model. They are checked across soil moisture
  # because the indirect term's share of the answer varies with it.
  pars <- c("vcmax_25", "jmax_25", "stem_P50", "TF24_cost_scale", "root_P50")
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
  pars <- c("vcmax_25", "stem_P50", "TF24_beta2")
  ift <- grid_gradient(2.0, pars = pars, method = "ift")
  fd <- grid_gradient(2.0, pars = pars, method = "fd")
  expect_equal(ift$gradient[, "collar"], fd$gradient[, "collar"],
               tolerance = 1e-3)
})

test_that("a pinned optimum takes the fallback, and the composite would be wrong", {
  # The hazard, pinned to specific operating points. These three are rows of the
  # golden grid; the boundary-soil pair are where the bare composite is worst.
  #
  # ⚠️ The truth at a boundary-soil point is ~1e-07 or smaller -- the leaf is stuck
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
                          pars = c("stem_P50", "vcmax_25"))
    expect_identical(auto$status, "pinned", label = lab)
    expect_identical(auto$method, "fd", label = lab)
    # The guard fired with room to spare: measured stationarity at these points
    # is ~6e-06 to ~9e-06 against a 1e-08 threshold.
    expect_gt(auto$stationarity, 1e-7)
    expect_lt(abs(auto$gradient["stem_P50", "A"]), 1e-3)

    # Forcing the composite here does not produce the wrong number -- it fails.
    # That was not the design and is worth stating as a measurement: psi* sits
    # 1e-06 of a bracket width from its bound at a pinned point, so the step in
    # psi cannot be centred without clamping, and that is checked. The grid test
    # below shows it holds at all 42 pinned rows, which means the composite's
    # O(1) answer is not reachable through this function.
    expect_error(grid_gradient(w$psi_soil, vpd = w$vpd, layers = w$layers,
                               pars = "stem_P50", method = "ift"),
                 "narrower than one step", label = lab)
  }

  # Pinned DRY is the milder case and is worth separating: the composite is only
  # a factor of 1.6 out for vcmax_25 there rather than seven orders, so a test
  # that only looked at the wet rows would let a threshold drift past it.
  dry <- grid_gradient(4, vpd = 0.5, layers = 3L,
                       pars = c("stem_P50", "vcmax_25", "R_d_25"))
  expect_identical(dry$status, "pinned")
  expect_equal(dry$gradient["vcmax_25", "A"], 0.017548, tolerance = 1e-3)
  # ⚠️ NOT the interior identity. At a boundary-crit optimum the derived psi_crit is
  # the binding constraint, so dA/dstem_P50 picks up a second term through it --
  # and both halves are independently arbitrated numbers from before the
  # reparameterisation (4.6991 for d/dstem_b, 1.2605 for d/dpsi_crit):
  #
  #     dA/dP50 = (db/dP50).dA/db + (dpsi_crit/dP50).dA/dpsi_crit
  #
  # Dropping the second term is the plausible-looking mistake this catches: it
  # would give 5.39 rather than 7.56, a 29% error with nothing about it that looks
  # wrong.
  tr <- leaf_traits()
  d_dP50 <- log(2)^(-1 / tr$stem_c)
  bound  <- (log(1 / 0.05) / log(2))^(1 / tr$stem_c)
  expect_equal(dry$gradient["stem_P50", "A"],
               d_dP50 * 4.6991 + bound * 1.2605, tolerance = 2e-3)

  # `dA/dvcmax_25` is a PARTIAL at fixed respiration, since `R_d_25` is its own
  # trait. A fit that wants respiration to follow Vcmax adds the two columns, so
  # that combination is asserted rather than left implicit.
  expect_equal(dry$gradient["vcmax_25", "A"] +
                 0.015 * dry$gradient["R_d_25", "A"],
               0.014357, tolerance = 1e-3)
})

test_that("at a pinned optimum the moving BOUND carries part of the gradient", {
  # The sharpest statement of why the fallback is not merely a fallback, and it
  # is a stronger claim than "the composite is inaccurate there".
  #
  # `psi_crit` sets the dry end of the feasible collar interval and appears in the
  # profit function nowhere else. It is no longer a trait -- it is the
  # 5%-conductivity quantile of the stem curve -- so moving `stem_P50` moves it,
  # and the effect that used to be a column of its own is now a TERM inside the
  # stem gradient. Where the optimum is interior that term is exactly zero; where
  # it is pinned at the dry bound the term is what the answer is made of.
  #
  # ⚠️ The composite cannot see it, and not by a small margin: psi* enters
  # `dprofit` nowhere, so M = 0, dpsi*/dtheta = 0 and the direct term through the
  # bound is 0. A composite forced through here understates the pinned gradient by
  # the whole bound term, and nothing about the number looks wrong.
  tr <- leaf_traits()
  d_dP50 <- log(2)^(-1 / tr$stem_c)
  bound  <- (log(1 / 0.05) / log(2))^(1 / tr$stem_c)

  # Interior: the bound term drops out, so the gradient is the b-gradient scaled.
  interior <- grid_gradient(2.0, pars = "stem_P50")
  expect_identical(interior$status, "interior")
  expect_equal(interior$gradient["stem_P50", "A"], d_dP50 * 2.8465,
               tolerance = 5e-3)

  # Pinned: the bound term is present and is 40% of the total.
  pinned <- grid_gradient(4, vpd = 0.5, layers = 3L, pars = "stem_P50")
  expect_identical(pinned$status, "pinned")
  expect_identical(pinned$method, "fd")
  # Both halves were arbitrated separately against least-squares slopes over +-2%
  # at n = 41, before the reparameterisation: 4.6991 for dA/dstem_b and 1.2605 for
  # dA/dpsi_crit (R^2 = 0.998; the fit is imperfect because the response has a
  # kink in it, which is the whole point).
  expect_equal(pinned$gradient["stem_P50", "A"],
               d_dP50 * 4.6991 + bound * 1.2605, tolerance = 2e-3)
  expect_gt(pinned$gradient["stem_P50", "collar"], 0)
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
  g <- grid_gradient(6.0, pars = c("vcmax_25", "stem_P50", "R_d_25"))
  expect_identical(g$status, "no-gradient")
  expect_identical(g$method, "fd")
  expect_equal(g$value[["A"]], -leaf_traits()$R_d_25)
  expect_equal(g$gradient["R_d_25", "A"], -1, tolerance = 1e-8)
  expect_identical(g$gradient["vcmax_25", "A"], 0)

  # ⚠️ `profit` NEEDS ITS OWN CLOSED FORM HERE, AND HAS ONE. This is the single
  # regime where `profit_` is written by a branch that leaves the other outputs
  # alone (hazard 8), so a recorded hex with nothing saying what it ought to be
  # would pin a number rather than a fact. Two identities close it:
  #
  #   E = 0, so A = -R_d exactly and R_d(25 C) = R_d_25 -- and the hydraulic cost
  #   does not depend on R_d_25 at all. So dprofit/dR_d_25 = dA/dR_d_25 = -1.
  l <- leaf_model()
  do.call(set_drivers, c(list(l), grid_drivers(6.0)))
  l$find_root_collar_psi()
  expect_identical(l$transpiration_, 0)
  expect_identical(l$assim_colimited_, -l$R_d_)
  expect_identical(g$value[["profit"]],
                   l$assim_colimited_ - l$hydraulic_cost_)
  expect_equal(g$gradient["R_d_25", "profit"], -1, tolerance = 1e-8)
  expect_equal(g$gradient["R_d_25", "profit"], g$gradient["R_d_25", "A"],
               tolerance = 1e-8)

  #   The shut-down collar is PINNED AT psi_crit, and psi_crit is the
  #   5%-conductivity quantile of the stem curve -- so dcollar/dstem_P50 is
  #   exactly (log2(1/0.05))^(1/stem_c), a constant at fixed stem_c.
  bound <- (log(1 / 0.05) / log(2))^(1 / leaf_traits()$stem_c)
  expect_equal(g$gradient["stem_P50", "collar"], bound, tolerance = 1e-6)
  expect_equal(g$gradient["stem_P50", "psi_stem"], bound, tolerance = 1e-6)
  expect_equal(g$gradient["stem_P50", "A"], 0)

  #   And PROFIT does not move either, which is a stronger claim than it looks and
  #   is the reparameterisation's doing. The collar moves with the trait, but it
  #   moves TO the same place on the curve: psi_crit is defined as the potential at
  #   5% remaining conductivity, so f(psi_crit) is exactly 0.05 whatever stem_P50
  #   is, and the hydraulic cost of a shut-down leaf is therefore scale-free.
  #   Asserted as bit-identical across three curves rather than as a gradient,
  #   because that is what it actually is.
  costs <- vapply(c(3.0, 3.4, 4.0), function(P50) {
    ll <- leaf_model(leaf_traits(stem_P50 = P50))
    do.call(set_drivers, c(list(ll), grid_drivers(6.0)))
    ll$find_root_collar_psi()
    ll$hydraulic_cost_
  }, numeric(1))
  expect_identical(costs[2:3], costs[c(1L, 1L)])
  expect_equal(g$gradient["stem_P50", "profit"], 0, tolerance = 1e-6)

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
  #     ~6e-06, so the default 1e-08 sits in an empty band four orders wide on
  #     each side. That is the number worth defending, because a threshold in a
  #     populated region would be a judgement call instead.
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
  expect_gt(mildest_pinned, 1e-7)

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
                     supply = leaf_supply_singlelayer(),
                     root_network = series_resistance(1e3),
                     pars = c("vcmax_25", "stem_P50"))
  expect_identical(g$status, "interior")
  fd <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                      supply = leaf_supply_singlelayer(),
                      root_network = series_resistance(1e3),
                      pars = c("vcmax_25", "stem_P50"), method = "fd")
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
  # makes them possible. The collar solve solves dprofit == 0 at
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
  # #44. Half of a companion study's free parameters are not traits: `K_total` and
  # `f_plant` reach the leaf as leaf_specific_conductance_max (a driver) and the
  # single-potential resistance. Nothing in the derivation cares -- the implicit
  # function theorem is applied to dprofit/dpsi = 0, and any parameter profit
  # depends on goes through it identically -- so the thing to check is that the
  # plumbing applies the perturbation through the right setter.
  #
  # Arbitrated the way 11c arbitrated the traits: against a central difference of
  # the WHOLE SOLVE, read where the step sweep plateaus rather than at one step.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_singlelayer(),
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
            supply = leaf_supply_singlelayer(),
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
              supply = leaf_supply_singlelayer(),
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
  # The stem cumulative-vulnerability integral is homogeneous of degree
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
      fast <- do.call(leaf_gradient, c(d, list(pars = "stem_P50", method = m,
                                               fast_stem_curve = TRUE)))
      slow <- do.call(leaf_gradient, c(d, list(pars = "stem_P50", method = m,
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
  # ⚠️ perturb_stem_P50() leaves the splines built at a different stem_b on
  # purpose, so the thing that could go wrong is a gradient call that does not
  # put it back -- after which every later solve is quietly on a rescaled curve.
  # set_traits() forces the rebuild for exactly this reason, and BIT-identity is
  # the right test: a restored leaf shares no code path with a fresh one.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0)
  before <- do.call(leaf_solve, d)

  invisible(do.call(leaf_gradient, c(d, list(pars = "stem_P50"))))
  invisible(do.call(leaf_gradient, c(d, list(pars = c("stem_P50", "vcmax_25")))))
  # ...and one that ends on the fast path, which is the ordering that would leave
  # it displaced if the restore at the end of the loop were missing.
  invisible(do.call(leaf_gradient, c(d, list(pars = c("vcmax_25", "stem_P50")))))

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
  # a companion study fits three responses -- A, gs and psi_leaf -- so all three
  # are differentiated, not just A. `collar` comes along because it is psi*, and
  # `profit` because it is what a DEMOGRAPHIC consumer bills: plant's carbon is
  # `leaf.profit_`, not `assim_colimited_`, so without it the four columns above
  # -- the calibration set -- were disjoint from plant's (#87).
  g <- grid_gradient(2.0, pars = "vcmax_25")
  expect_identical(colnames(g$gradient),
                   c("A", "gc", "psi_stem", "collar", "profit"))
  expect_identical(names(g$value),
                   c("A", "gc", "psi_stem", "collar", "profit"))
  expect_true(all(is.finite(g$gradient)))
  # Raising vcmax_25 raises assimilation and opens the stomata, and the leaf pays
  # for it with a more negative water potential (a larger positive magnitude).
  expect_gt(g$gradient["vcmax_25", "A"], 0)
  expect_gt(g$gradient["vcmax_25", "gc"], 0)
  expect_gt(g$gradient["vcmax_25", "psi_stem"], 0)
  # More photosynthetic capacity is worth having: profit is A minus the hydraulic
  # cost, and the extra carbon exceeds the extra cost. It is BELOW dA/dvcmax_25,
  # which is the statement that the cost is not zero -- the cheapest available
  # check that this column is the objective and not a copy of A.
  expect_gt(g$gradient["vcmax_25", "profit"], 0)
  expect_lt(g$gradient["vcmax_25", "profit"], g$gradient["vcmax_25", "A"])
})

test_that("profit's gradient is the direct term alone at an interior optimum", {
  # The envelope theorem, which is the ONE place this package uses it. At an
  # interior optimum dprofit/dpsi == 0, so the indirect term
  # (dprofit/dpsi)(dpsi*/dtheta) vanishes IDENTICALLY, and dprofit/dtheta is the
  # direct partial at fixed psi. `leaf_gradient()` encodes that by ASSIGNING the
  # profit column from the direct term -- it does not zero dY_dpsi["profit"],
  # which carries the exact dprofit/dpsi for the pinned route -- so the check is
  # that the reported column equals a direct central difference with the collar
  # HELD at psi*.
  #
  # ⚠️ THE OPERATING POINT IS CHOSEN, NOT ARBITRARY, and choosing it is what makes
  # this a test rather than a formality. The assignment only matters where the term
  # it removes is big enough to see, and that term is NOISE rather than an h^2
  # truncation: `profit` is the maximum, so it is flat, and a central difference
  # of it divides the solve's ~1e-9 floor by a ~1e-6 step. So the point has to be
  # one where that noise is LARGE; at the suite's usual `grid_drivers(2.0)` the same
  # term is ~1e-10 and this test would pass either way, which is the version of it
  # written first.
  #
  # ⚠️ THE POINT DOES NOT SURVIVE A RESULTS CHANGE, so RE-DERIVE IT rather than
  # adjusting the bounds whenever this fails on the `fd` assertion: sweep the golden
  # grid's interior rows and take the largest |fd|.
  #
  #   for each (psi_soil, ppfd, vpd, layers) in the golden grid:
  #     g <- leaf_gradient(..., pars = "vcmax_25"); keep if g$status == "interior"
  #     solve, then fd <- central difference of profit_ in the collar at psi*
  #   pick argmax |fd|
  #
  # A failure saying `abs(fd) <= 1e-5` means the point went quiet, NOT that the
  # feature broke. A failure on `abs(exact)` is the real signal.
  #
  # ⚠️ Those figures are macOS/arm64's, and the second half of this test says why
  # that matters. Read it before adding an assertion on a magnitude here.
  d <- grid_drivers(1.0, ppfd = 1500, vpd = 1.0, layers = 1L)
  g <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
  expect_identical(g$status, "interior")

  l <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
  do.call(set_drivers, c(list(l), d))
  l$find_root_collar_psi()
  psi_star <- l$opt_root_psi_

  # The two instruments, both read on a leaf still at base traits.
  hp <- max(abs(psi_star), 1) * 1e-6
  exact <- l$dprofit_droot_collar_psi(psi_star)
  l$evaluate_root_collar_psi(psi_star + hp)
  hi <- l$profit_
  l$evaluate_root_collar_psi(psi_star - hp)
  lo <- l$profit_
  fd <- (hi - lo) / (2 * hp)

  # The answer itself, against a direct difference at fixed collar. THIS IS THE
  # PART THAT HOLDS EVERYWHERE -- it is the claim the feature makes, and it does
  # not depend on how big the term that was dropped happens to be.
  v <- leaf_traits()$vcmax_25
  h <- max(abs(v), 1) * 1e-6
  at <- function(x) {
    set_traits(l, leaf_traits(vcmax_25 = x))
    do.call(set_drivers, c(list(l), d))
    l$evaluate_root_collar_psi(psi_star)
    l$profit_
  }
  direct <- (at(v + h) - at(v - h)) / (2 * h)
  expect_equal(g$gradient["vcmax_25", "profit"], direct, tolerance = 1e-7)

  # ⚠️ THE MAGNITUDES BELOW ARE PLATFORM-SPECIFIC, AND THE FIRST VERSION OF THIS
  # TEST ASSERTED THEM EVERYWHERE. It passed on macOS/arm64 and failed on Linux
  # CI, where the same operating point gives `exact` = 6.6e-11 and `fd` = 8.0e-10
  # rather than 2.4e-15 and 2.1e-04.
  #
  # That is not a different answer, it is a different NOISE FLOOR: which side of
  # the collar solver's tolerance the root-find lands on is set by libm's exp/pow,
  # and those are not reproducible between Apple's arm64 libm and glibc on
  # x86-64. The package already has a name for that -- it is the same reason the
  # golden files compare bit-exactly on one platform and with a tolerance
  # elsewhere -- so the same predicate gates it here.
  #
  # ⚠️ So this test has TEETH only on the platform the measurement was made on.
  # Off it, the check above still holds and the one below is skipped; a Linux-only
  # run would not catch the assignment being removed. Said plainly rather than left
  # for someone to discover from a green CI.
  skip_if_not(golden_bit_exact_platform(),
              "the noise floor these two numbers measure is macOS/arm64's")
  # The exact instrument says stationary to solver precision; the differenced one
  # says nothing of the kind. Asserting BOTH is the point -- it is the difference
  # between "the dropped term is small" and "the dropped term is unmeasurable by
  # the route that would have supplied it".
  #
  # ⚠️ THE BOUND ON `exact` IS THE SOLVER FLOOR (~1e-9), not a landing point. Pinning
  # it tighter pins which side of its own tolerance the collar root-find lands on,
  # which moves for reasons unrelated to this feature.
  #
  # What carries the test is the GAP, and it is asserted as a RATIO rather than as
  # a magnitude on `fd`.
  #
  # ⚠️ `fd`'s absolute size is not a property of the feature -- it is where the
  # collar root-find happened to land inside its own tolerance, divided by the
  # step. The (P50, c) reparameterisation moved every parameter by ~5e-08 and
  # `fd` at this point went from 2.7e-04 to 3.6e-10, four orders, with nothing
  # about the model changed. A bound on `fd` therefore pins a coincidence; the
  # ratio is the claim. Measured across eight operating points it runs 1e+04 to
  # 1e+09, so 1e+03 is inside every one of them.
  expect_lt(abs(exact), 1e-9)
  expect_gt(abs(fd) / abs(exact), 1e3)
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
  tr <- leaf_traits(vcmax_25 = 105, stem_P50 = 4.2)
  args <- list(psi_soil = 2.0, PPFD = 900, atm_vpd = 1.5, traits = tr,
               pars = c("vcmax_25", "stem_P50", "TF24_cost_scale"))
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
  # made seating the leaf with reset() rather than set_drivers() load-bearing.
  tr <- leaf_traits(vcmax_25 = 105, stem_P50 = 4.2)
  args <- list(psi_soil = 2.0, PPFD = 900, traits = tr, pars = "vcmax_25")
  fresh <- do.call(leaf_gradient, args)

  # Deliberately different, and deliberately still inside the #38 domain -- the
  # vessel's traits are overwritten before any solve, but leaf_model() does build
  # its splines, and an invalid vessel would confuse the diagnosis if this failed.
  wrong <- leaf_model(traits = leaf_traits(vcmax_25 = 60, stem_P50 = 4.6))
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
  s <- leaf_model(supply = leaf_supply_singlelayer())
  g <- leaf_gradient(psi_soil = 1.5, x = s, traits = tr,
                     root_network = series_resistance(1e3), pars = "resistance")
  expect_true(is.finite(g$gradient["resistance", "A"]))
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr, pars = "resistance"),
               "single-potential path only")
})

test_that("the setter's positional trait call cannot drift in arity", {
  # `.gradient_setter()` applies traits POSITIONALLY, straight onto the object, to
  # skip rebuilding a leaf_traits per perturbation. It used to name a fixed number
  # of subscripts, and a trait added to leaf_traits() and to the C++ setter would
  # silently pass the wrong value for every argument after the new one -- so this
  # test asserted the arity. It is a `do.call` over the derived trait vector now,
  # which cannot drift; what is still worth asserting is that the three lists agree
  # with each other, since the setter takes its ORDER from the C++ enumeration.
  expect_length(leaf_traits(), 15L)
  expect_length(formals(leaf_model()$set_traits), 15L)
  expect_identical(names(leaf_traits()), names(formals(leaf_model()$set_traits)))
  expect_identical(names(leaf_traits()),
                   gradient_par_names()[seq_along(leaf_traits())])
})

# --- a collar potential the caller supplies (#88) -----------------------------

test_that("a prescribed psi at psi* reproduces the solving path bit-for-bit", {
  # ⚠️ THE LOAD-BEARING TEST OF THE WHOLE FEATURE, and the reason `psi` reuses
  # `.gradient_ift()` rather than getting a composite of its own. Give it back
  # the collar it solved for, and the response of that collar that it derived,
  # and the two paths are the same arithmetic -- so this is `identical()` and not
  # a tolerance. Anything that made the prescribed path a second implementation
  # would show up here as a last-bit difference rather than as a design note.
  d <- grid_drivers(2.0)
  pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale")
  a <- do.call(leaf_gradient, c(d, list(pars = pars)))
  expect_identical(a$status, "interior")

  b <- do.call(leaf_gradient,
               c(d, list(pars = pars, psi = a$psi, dpsi_dtheta = -a$M / a$H)))
  expect_identical(b$status, "prescribed")
  expect_identical(b$method, "prescribed")
  expect_identical(b$gradient, a$gradient)
  expect_identical(b$M, a$M)
  expect_identical(b$dY_dpsi, a$dY_dpsi)
  expect_identical(b$psi, a$psi)

  # `collar` IS psi, so its column must be exactly what was supplied -- the one
  # column whose value the caller controls outright.
  expect_identical(unname(b$gradient[, "collar"]), unname(-a$M / a$H))

  # And the equivalence is not an accident of the default: dpsi_dtheta = 0 gives
  # the direct terms, which differ from the solving path by the indirect term.
  z <- do.call(leaf_gradient, c(d, list(pars = pars, psi = a$psi)))
  expect_identical(unname(z$gradient[, "collar"]), rep(0, length(pars)))
  expect_true(all(abs(z$gradient[, "A"] - a$gradient[, "A"]) > 1e-12))
})

test_that("a prescribed psi off the optimum matches a difference at fixed collar", {
  # The claim the feature makes: at a collar the caller chose, dY/dtheta is the
  # direct partial. Checked against a central difference computed the long way,
  # at a psi deliberately away from psi* so the indirect term the solving path
  # would add is large and its absence is visible.
  d <- grid_drivers(2.0)
  a <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
  psi_off <- a$psi + 0.2
  g <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25", psi = psi_off)))
  expect_identical(g$status, "prescribed")
  # Far from stationary, which is the regime this exists for.
  expect_gt(g$stationarity, 1e-3)

  l <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
  v <- leaf_traits()$vcmax_25
  h <- max(abs(v), 1) * 1e-6
  at <- function(x, field) {
    set_traits(l, leaf_traits(vcmax_25 = x))
    do.call(set_drivers, c(list(l), d))
    l$evaluate_root_collar_psi(psi_off)
    l[[field]]
  }
  fields <- c(A = "assim_colimited_", profit = "profit_")
  for (column in names(fields)) {
    ref <- (at(v + h, fields[[column]]) - at(v - h, fields[[column]])) / (2 * h)
    expect_equal(g$gradient["vcmax_25", column], ref, tolerance = 1e-7,
                 label = column)
  }

  # profit's dY/dpsi is the EXACT dprofit/dpsi here, not a difference of it --
  # the same rule as the envelope, reaching the opposite answer because the point
  # is not stationary. Compare against the primitive directly.
  set_traits(l, leaf_traits())
  do.call(set_drivers, c(list(l), d))
  expect_identical(g$dY_dpsi[["profit"]], l$dprofit_droot_collar_psi(psi_off))
  expect_gt(abs(g$dY_dpsi[["profit"]]), 1)
})

test_that("a clamped psi reports itself and returns no gradient", {
  # ⚠️ NOT AN ERROR, AND NOT THE DIRECT TERM EITHER. The collar used is `psi`
  # clamped into the feasible interval, so it moves with the BOUND rather than
  # with the caller's dpsi_dtheta -- the active-set problem arriving through the
  # clamp. The direct term alone would be plausible and wrong, so it is withheld;
  # a tracking model reaches these points routinely, so it is reported and not
  # thrown.
  d <- grid_drivers(2.0)
  pars <- c("vcmax_25", "stem_P50")

  g <- do.call(leaf_gradient, c(d, list(pars = pars, psi = 99)))
  expect_identical(g$status, "clamped")
  expect_identical(g$method, "prescribed")
  expect_true(all(is.na(g$gradient)))
  expect_true(all(is.na(g$M)))
  expect_true(all(is.na(g$dY_dpsi)))
  # What the point IS still comes back: `psi` is the collar actually used, and
  # `value` describes it. That is what lets a caller see where it was pulled to.
  expect_lt(g$psi, 99)
  expect_true(all(is.finite(g$value)))

  # And it fires one step inside the end too, because dY/dpsi cannot be centred
  # there -- a one-sided difference over a shortened interval is the same failure
  # the solving path's clamp detector refuses.
  a <- do.call(leaf_gradient, c(d, list(pars = pars)))
  expect_identical(
    do.call(leaf_gradient, c(d, list(pars = pars, psi = g$psi)))$status,
    "clamped")
  # ...while a psi comfortably inside it does not.
  expect_identical(
    do.call(leaf_gradient, c(d, list(pars = pars, psi = a$psi)))$status,
    "prescribed")
})

test_that("an integer psi is not mistaken for a clamped one", {
  # ⚠️ THE CLAMP TEST IS `identical()` AGAINST WHAT THE CALLER PASSED, and
  # `identical(3, 3L)` is FALSE. So an integer `psi` -- `psi = 3L`, or the
  # entirely ordinary `for (p in 2:5)`, since `2:5` yields integers -- was
  # reported CLAMPED at a collar it had been given exactly, with an all-NA
  # gradient and a `psi` in the result showing no movement, so there was nothing
  # to notice.
  #
  # Neither sibling had it: `.gradient_check_psi_batch()` already coerced, and
  # C++'s `util::identical` is `a == b`. R alone, and only because the checker's
  # return was discarded.
  d <- grid_drivers(2.0)
  a <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
  whole <- 3                        # inside the feasible interval at this point
  dbl <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25", psi = whole)))
  int <- do.call(leaf_gradient,
                 c(d, list(pars = "vcmax_25", psi = as.integer(whole))))
  expect_identical(dbl$status, "prescribed")
  expect_identical(int$status, "prescribed")
  expect_identical(int$gradient, dbl$gradient)
  expect_identical(int$psi, dbl$psi)

  # And the batch entry point agrees on the same input, which it did not before.
  b <- leaf_batch(psi_soil = 2.0, PPFD = 900)
  expect_identical(
    leaf_gradient_batch(b, pars = "vcmax_25",
                        psi = as.integer(whole))$status, "prescribed")
})

test_that("an infeasible prescribed psi is no-gradient, not a sentinel zero", {
  # ⚠️ `dprofit_droot_collar_psi` RETURNS A HARD 0.0 ON ITS SHUT-DOWN EXIT, and a
  # bare zero is indistinguishable from a stationary point. The solving path got
  # away with reading the value alone because `H` collapses to zero with it and
  # `usable` catches the pair; the prescribed path never divides by `H`, so it
  # would have adopted the sentinel AS dprofit/dpsi -- and with a non-zero
  # `dpsi_dtheta`, which is the case this feature exists for, silently lost
  # profit's indirect term at exactly the dry points a tracking model lives in.
  #
  # The fix reads the feasibility flag the header has always had and R could not
  # reach. Most such points are caught as `clamped` first -- the shut-down state
  # seats a collar of its own choosing -- so the test hands back exactly that
  # collar, which is the case `clamped` cannot see.
  d <- grid_drivers(6.0)                       # drier than psi_crit: shut down
  solved <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
  expect_identical(solved$status, "no-gradient")

  g <- do.call(leaf_gradient,
               c(d, list(pars = "vcmax_25", psi = solved$psi,
                         dpsi_dtheta = 1)))
  # Unclamped -- it is the collar the shut-down state itself seated -- and still
  # refused, which is the whole point of the test.
  expect_identical(g$psi, solved$psi)
  expect_identical(g$status, "no-gradient")
  expect_true(all(is.na(g$gradient)))
  expect_true(all(is.na(g$dY_dpsi)))

  b <- leaf_batch(psi_soil = 6.0, PPFD = 900)
  expect_identical(leaf_gradient_batch(b, pars = "vcmax_25",
                                       psi = solved$psi)$status, "no-gradient")
})

test_that("psi and dpsi_dtheta are validated together", {
  d <- grid_drivers(2.0)
  gr <- function(...) do.call(leaf_gradient, c(d, list(pars = "vcmax_25", ...)))
  expect_error(gr(psi = -1), "positive")
  expect_error(gr(psi = 0), "positive")
  expect_error(gr(psi = c(1, 2)), "single")
  expect_error(gr(psi = NA_real_), "finite")
  # `method` is about a solved optimum, so it cannot come with a collar the
  # caller chose -- "fd" especially, which would difference the solve and answer
  # a question about the optimum instead.
  expect_error(gr(psi = 2, method = "fd"), "cannot be given with `psi`")
  expect_error(gr(psi = 2, method = "ift"), "cannot be given with `psi`")
  # dpsi_dtheta without psi is meaningless: there the collar's response is
  # derived, not supplied.
  expect_error(gr(dpsi_dtheta = 1), "needs `psi`")
  expect_error(gr(psi = 2, dpsi_dtheta = c(1, 2)), "length 1 or one value")
  expect_error(gr(psi = 2, dpsi_dtheta = c(nope = 1)), "names must be exactly")

  # ⚠️ AN INFINITE `dpsi_dtheta` IS REFUSED, AND THE ENVELOPE IS WHY. The
  # composite is `direct + dY_dpsi * dpsi_dtheta`, so four columns come back
  # +-Inf -- but `profit` is ASSIGNED from the direct term at a stationary psi
  # rather than multiplied through, so it came back finite and plausible beside
  # them. Measured at `psi = psi*` with `dpsi_dtheta = Inf`: A, gc, psi_stem and
  # collar all Inf, profit 0.0105. Reading `profit` alone is plant's own case
  # (#87), so the one surviving column was the one most likely to be believed.
  #
  # `anyNA` already caught NA and NaN; `Inf` it did not, and `is.finite` is the
  # single check that covers all three.
  expect_error(gr(psi = 2, dpsi_dtheta = Inf), "finite")
  expect_error(gr(psi = 2, dpsi_dtheta = -Inf), "finite")
  expect_error(gr(psi = 2, dpsi_dtheta = NaN), "finite")
  # At psi* specifically, which is where the envelope fires and where the
  # asymmetry was, rather than only at an arbitrary prescribed collar.
  solved <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
  expect_error(gr(psi = solved$psi, dpsi_dtheta = Inf), "finite")

  # The batch entry point validates the same way, on BOTH shapes it accepts --
  # the matrix form never reaches `.gradient_dpsi_dtheta()`, so one check there
  # would have left it open.
  b <- leaf_batch(psi_soil = c(2.0, 2.5), PPFD = 900)
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25", psi = c(3, 3),
                                   dpsi_dtheta = Inf), "finite")
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25", psi = c(3, 3),
                                   dpsi_dtheta = matrix(c(1, Inf), nrow = 2)),
               "finite")

  # Named is matched by name, not by position -- the same discipline `pars` has,
  # for the same reason.
  pars <- c("vcmax_25", "stem_P50")
  by_name <- do.call(leaf_gradient,
                     c(d, list(pars = pars, psi = 3.0,
                               dpsi_dtheta = c(stem_P50 = 2, vcmax_25 = 1))))
  by_pos <- do.call(leaf_gradient,
                    c(d, list(pars = pars, psi = 3.0, dpsi_dtheta = c(1, 2))))
  expect_identical(by_name$gradient, by_pos$gradient)
})

# --- per-model parameter availability (plan item 4) --------------------------

test_that("CF77_lambda_ is in the enumeration and refused by other models", {
  # It is Cowan & Farquhar's only parameter. Every other curve's lambda is
  # EMERGENT -- derived from that curve's own parameters rather than set -- so the
  # slot exists once in the enumeration and is available for one model.
  expect_true("CF77_lambda_" %in% gradient_par_names())
  # ⚠️ NO LONGER LAST. `TF24_floor_lambda_o` is the second model-owned price and was
  # appended after it, which is exactly the "appending is safe" the enumeration is
  # designed for -- so what is asserted is the pair and its order, not a position.
  expect_identical(tail(gradient_par_names(), 2L),
                   c("CF77_lambda_", "TF24_floor_lambda_o"))

  # ⚠️ REFUSED, not returned as a zero. A zero means "in this objective but
  # inactive at this operating point", which is the psi_crit case and is
  # information. Structurally absent from the objective is a different statement,
  # and an optimiser handed a zero would read the parameter as inert and converge.
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "CF77_lambda_"),
               "not a parameter of the collar model")
  # And the message names the axis that ruled it out, so a caller is not left
  # looking at the supply path.
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "CF77_lambda_"),
               "emergent")

  # The default must not offer what it would then refuse.
  expect_false("CF77_lambda_" %in%
                 phylloptim:::.gradient_available_pars(TRUE, "collar"))
  expect_true("CF77_lambda_" %in%
                phylloptim:::.gradient_available_pars(TRUE, "CF77"))

  # The two axes are independent: the supply path still rules out `resistance`.
  multi <- phylloptim:::.gradient_available_pars(FALSE, "CF77")
  expect_false("resistance" %in% multi)
  expect_true("CF77_lambda_" %in% multi)

  # An unknown model is refused rather than silently offering everything.
  expect_error(phylloptim:::.gradient_available_pars(TRUE, "Medlyn"),
               "no gradient route")
})

test_that("TF24_floor_lambda_o is the second model-owned price, and refused elsewhere", {
  # ⚠️ THE CLASS, NOT THE INSTANCE. There are two prescribed prices in the
  # enumeration now -- Cowan-Farquhar's and the two-term curve's -- so the thing
  # under test is that `.gradient_owned_pars` maps each to its owner rather than
  # that one name is special. On every other curve the price is EMERGENT, derived
  # from that curve's own parameters rather than set.
  expect_true("TF24_floor_lambda_o" %in% gradient_par_names())

  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "TF24_floor_lambda_o"),
               "not a parameter of the collar model")
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "TF24_floor_lambda_o"),
               "emergent")
  # And the refusal names the OWNER, so a caller is told what to pass rather than
  # only what is wrong.
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "TF24_floor_lambda_o"),
               "TF24_floor")

  # The two prices exclude each other, which is what makes this a table rather
  # than a special case.
  cf <- phylloptim:::.gradient_available_pars(TRUE, "CF77")
  tt <- phylloptim:::.gradient_available_pars(TRUE, "TF24_floor")
  expect_true("CF77_lambda_" %in% cf)
  expect_false("TF24_floor_lambda_o" %in% cf)
  expect_true("TF24_floor_lambda_o" %in% tt)
  expect_false("CF77_lambda_" %in% tt)

  # ⚠️ AND THE HYDRAULIC HALF IS NOT IN THIS CLASS AT ALL. `TF24_floor` reads
  # `TF24_cost_scale`, which is a TRAIT and available on every model -- returning
  # an exactly zero gradient off the curves that do not read it, which is the
  # `psi_crit` case and is information rather than an error. The distinction
  # between "zero" and "refused" is the whole design here, so both are checked.
  expect_true("TF24_cost_scale" %in%
                phylloptim:::.gradient_available_pars(TRUE, "TF24_floor"))
  g <- leaf_gradient(psi_soil = 2.0, PPFD = 900, pars = "JS22_gamma")
  expect_identical(unname(g$gradient["JS22_gamma", "A"]), 0)
})

test_that("the curve registry is read from C++, not restated in R", {
  # R selects a curve by POSITION in this vector, so it is compared rather than
  # trusted -- the same discipline gradient_par_names() gets.
  nms <- cost_curve_names()
  expect_identical(nms, c("TF24", "CF77", "JS22", "CMax", "SOX", "JW26",
                          "ProfitMax", "TF24_floor"))

  # ⚠️ THERE IS NO `cost_curve_has_derivative()` ANY MORE, and its absence is the
  # statement. Every curve is `h(A) - C(psi)` and every derivative is
  # `h'(A)*dA/dpsi - dC/dpsi`, so the predicate returned TRUE for all seven and
  # nothing read it -- a constant dressed as a query. Product objectives are the
  # `log` link and ProfitMax the `A/|A|max` one; neither is a special case.
  expect_false(exists("cost_curve_has_derivative",
                      where = asNamespace("phylloptim"), inherits = FALSE))

  # ⚠️ AND THERE IS NO SECOND LIST OF MODEL NAMES EITHER. `.gradient_link_verified`
  # / `.GRADIENT_VERIFIED_LINKS` recorded which routes had been checked against a
  # finite difference of their own solve; once every route was, the predicate was
  # a constant and the `stop()` behind it unreachable -- and that dead branch still
  # named `$optimise_psi_stem_<curve>()`, which no longer exists. Deleted. What
  # replaces it is derived, so a curve added in C++ arrives here on its own.
  expect_identical(phylloptim:::.gradient_models(), c("collar", nms))
  expect_false(exists(".gradient_link_verified",
                      where = asNamespace("phylloptim"), inherits = FALSE))

  # Exactly one curve owns `CF77_lambda_`; every other model refuses it, and
  # nothing else is refused on the single-potential path.
  for (m in phylloptim:::.gradient_models()) {
    got <- phylloptim:::.gradient_available_pars(TRUE, m)
    expect_identical("CF77_lambda_" %in% got, identical(m, "CF77"),
                     info = m)
  }

  # The C++ derivative itself answers for a product curve rather than refusing --
  # that is the point of the link -- and it is reached by NAME, through the model
  # the leaf is seated on, with no curve index anywhere.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(l, psi_soil = 1.5, PPFD = 900)
  l$set_model("SOX", "stem"); l$optimise()
  expect_identical(l$model_curve(), "SOX")
  expect_identical(l$model_route(), "stem")
  d <- l$dprofit_dpsi_stem_checked(l$opt_psi_stem_)
  expect_length(d, 2L)
  expect_true(is.finite(d[[1]]))

  # ⚠️ AND A STEM METHOD REFUSES ON THE COLLAR ROUTE rather than answering about a
  # model the leaf is not seated on. This is what the integer argument used to
  # make unrepresentable-by-accident: an index said which curve, and nothing said
  # which route.
  l$set_model("SOX", "collar")
  expect_error(l$dprofit_dpsi_stem_checked(l$opt_psi_stem_), "STEM-route")
  expect_error(l$evaluate_psi_stem_at(2.0), "STEM-route")
})

# --- differentiating a stem model's optimum (plan item 5) --------------------

test_that("every additive curve's IFT gradient matches differencing its solve", {
  # ⚠️ The IFT is the EXACT one and the finite difference is the noisy one -- a
  # difference of a solve divides that solve's ~1e-09 floor by the step -- so this
  # bounds the FD's error rather than the composite's. What it establishes is that
  # the composite is differentiating the same model the optimiser solved.
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  mk <- function() {
    l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())
    l$CF77_lambda_ <- 1.5e5   # in band: marginal_cost_water runs 9e4-3e5
    l$TF24_floor_lambda_o <- 1.5e5   # the same band, and the same reason
    l
  }
  # ⚠️ `TF24_floor` APPEARS TWICE ON PURPOSE. It is the one curve with a parameter in
  # each half of the decomposition -- `TF24_cost_scale` in `Theta~(psi)` and
  # `TF24_floor_lambda_o` in the price -- and the two reach the first-order
  # condition by different routes: the hydraulic term is a function of the
  # potential alone, while `lambda_o*kmax*f(psi)` carries the supply. A composite
  # that had them the right way round on one and not the other would pass on
  # either alone. Note the first is a parameter the curve SHARES with TF24, which
  # is exactly the point of the shared parameterisation.
  cases <- list(c("TF24", "TF24_cost_scale"), c("CF77", "CF77_lambda_"),
                c("JS22", "JS22_gamma"), c("CMax", "CMax_a"),
                c("TF24_floor", "TF24_cost_scale"),
                c("TF24_floor", "TF24_floor_lambda_o"))
  n_checked <- 0L
  for (case in cases) {
    m <- case[[1]]
    p <- case[[2]]
    g <- function(meth) {
      leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net, x = mk(),
                    traits = tr, pars = p, model = m, method = meth)
    }
    ift <- g("ift")
    fd <- g("fd")
    expect_identical(ift$model, m)
    # ⚠️ THIS ASSERTED `"pinned"` AND A RESIDUAL ABOVE 1e-09, AND IT FAILED ON
    # PURPOSE. It was written to fail the day a stem route gained a
    # first-order-condition polish, and that is what happened: the stem entry
    # points now root-find dJ/dpsi == 0 through the same solver the collar route
    # uses, so an interior point is CLASSIFIED interior instead of being read as a
    # bound by its own solver residual.
    expect_identical(ift$status, "interior")
    expect_lt(ift$stationarity, 1e-9)
    # And `method = "auto"` therefore reaches the composite on a stem route rather
    # than falling back to differencing the solve. Still forced here so the arm
    # under test cannot change silently.
    expect_identical(ift$method, "ift")
    for (out in c("A", "gc", "psi_stem")) {
      expect_equal(ift$gradient[p, out], fd$gradient[p, out],
                   tolerance = 5e-3, label = paste(m, out))
    }
    n_checked <- n_checked + 1L
  }
  expect_identical(n_checked, 6L)
})

test_that("the decision variable follows the route", {
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())

  collar <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, pars = "TF24_cost_scale")
  stem <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net, x = l,
                        traits = tr, pars = "TF24_cost_scale", model = "TF24")

  # On the collar route `collar` IS psi*, so its derivative is assigned.
  expect_identical(collar$model, "collar")
  expect_true(is.finite(collar$gradient[1, "collar"]))
  # On a stem route the collar is never solved for, so NOT finite is the truthful
  # answer -- reading a number there would be reading stale or absent state.
  expect_false(is.finite(stem$gradient[1, "collar"]))
  expect_true(is.finite(stem$gradient[1, "psi_stem"]))

  # ⚠️ And they are DIFFERENT MODELS, not one restricted to the other: same
  # drivers, same cost curve, different optimum.
  expect_false(isTRUE(all.equal(collar$psi, stem$psi)))
})

test_that("a product objective's gradient survives an h-sweep", {
  # ⚠️ THE SWEEP IS THE TEST, not a single step. A finite difference of a
  # SCAN-resolved argmax is quantisation noise at a small step: psi* moves ~3e-06
  # MPa under a 1e-06 relative perturbation while the scan resolves it to ~4e-06
  # over a 4.4 MPa bracket. At that step the FD reported 0.1855 against a true
  # 0.0551 for SOX, and got the SIGN wrong for JW26 -- which is what a false bug
  # report is made of. The composite is flat across six decades; the FD is not.
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  mk <- function(v, p50) {
    l <- leaf_model(leaf_traits(vcmax_25 = v, stem_P50 = p50), leaf_control(),
                    leaf_supply_singlelayer())
    set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
    l
  }
  n_checked <- 0L
  for (m in c("SOX", "JW26")) {
    for (par in c("vcmax_25", "stem_P50")) {
      base <- if (par == "vcmax_25") 96 else 3.4
      solve_at <- function(x) {
        l <- if (par == "vcmax_25") mk(x, 3.4) else mk(96, x)
        l$set_model(if (m == "TF") "TF24" else m, "stem")
        l$optimise()
        l$assim_colimited_
      }
      h <- base * 1e-2               # the floor of the V, not the toe
      fd <- (solve_at(base + h) - solve_at(base - h)) / (2 * h)
      g <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net,
                         x = mk(96, 3.4), traits = tr, pars = par, model = m,
                         method = "ift")
      expect_equal(g$gradient[par, "A"], fd, tolerance = 1e-2,
                   label = paste(m, par))
      n_checked <- n_checked + 1L
    }
  }
  expect_identical(n_checked, 4L)
})

test_that("the remaining route restrictions are the topological ones", {
  tr <- leaf_traits()
  # A stem route needs the single-potential path, and says so rather than
  # returning a number from the wrong topology.
  expect_error(leaf_gradient(psi_soil = 1.5, PPFD = 1500,
                             pars = "TF24_cost_scale", model = "TF24"),
               "single-potential supply path")
  expect_error(leaf_gradient(psi_soil = 1.5, PPFD = 1500,
                             pars = "TF24_cost_scale", model = "Medlyn"),
               "no gradient route")
})

test_that("ProfitMax's normaliser follows the traits, so the composite is a total", {
  net <- series_resistance(1e4)
  mkv <- function(v) {
    l <- leaf_model(leaf_traits(vcmax_25 = v), leaf_control(),
                    leaf_supply_singlelayer())
    set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
    l
  }
  g <- function(meth) {
    leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net, x = mkv(96),
                  traits = leaf_traits(), pars = "vcmax_25", model = "ProfitMax",
                  method = meth)
  }
  # ⚠️ IT USED TO BE NaN, THEN A PARTIAL, AND IS NOW A TOTAL -- three states, and
  # the middle one is worth knowing about. `|A|max` is state a perturbation
  # destroys (`set_traits` + `set_physiology` clears it: 16.757 -> NaN), so
  # something must restore it. Pinning it at the base point fixed the NaN but froze
  # the normaliser, which made the composite a PARTIAL: 0.0766 against a true
  # 0.0579, the omitted term 32% of the answer. That was the right call only while
  # `|A|max` was the argmax of a 500-point scan and so piecewise constant in the
  # traits. It is a root-find now, so it is RE-SOLVED per perturbation and follows
  # the traits.
  expect_true(is.finite(g("ift")$gradient[1, "A"]))

  # The composite and a difference of the solve now agree, which is the statement
  # that the split is gone. Loose because the FD is the noisy side.
  expect_equal(g("ift")$gradient[1, "A"], g("fd")$gradient[1, "A"],
               tolerance = 1e-3)

  # And a step-swept difference of the solve, which is what a fit sees.
  h <- 96 * 1e-2
  solve_A <- function(v) {
    l <- mkv(v)
    l$set_model("ProfitMax", "stem"); l$optimise()
    l$assim_colimited_
  }
  swept <- (solve_A(96 + h) - solve_A(96 - h)) / (2 * h)
  expect_equal(g("ift")$gradient[1, "A"], swept, tolerance = 1e-3)

  # ⚠️ `auto` NO LONGER SPECIAL-CASES THIS CURVE. It used to be forced onto the
  # finite difference because the composite returned the partial; now every route
  # reaches the composite at an interior optimum.
  expect_identical(g("auto")$method, "ift")
  expect_identical(g("auto")$status, "interior")
})
