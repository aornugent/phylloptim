# The closed form as a third axis on set_model(): `method = "closed"`.
#
# WHAT THE METHOD IS. Every model here satisfies dA/dE = lambda at the optimum,
# and given lambda the solution collapses to the Medlyn USO relation ci/ca =
# xi/(xi + sqrt(D)) with xi = sqrt(Q/lambda). Where lambda is a power law in psi
# at the wet end that inverts explicitly, so a leading-order psi plus one Newton
# step on the supply-minus-demand residual replaces a root-find over the whole
# bracket. It is a different METHOD for the same MODEL, so the two must agree
# wherever the approximation holds -- which is what this file measures.
#
# WHAT IT IS NOT. It is not exact and it is not universal. Three separate
# approximations sit inside it (the power law for lambda, the one Newton step,
# and the USO collapse itself, whose 3*Gstar is the electron-transport-limited
# dA/dci), and the last of those is the binding one: see `closed_form.hpp`. So
# the tests below assert analytical IDENTITIES where identities hold -- the n =
# 0 arm's lambda, the supply inversion's round trip, the fallback's bit-equality
# -- and report the accuracy as a table binned by ci/ca rather than pinning a
# number.
#
# ⚠️ NO SNAPSHOTS. There are none anywhere in this package and there should be
# none: a snapshot of a driver grid records what the code does, which is exactly
# what a test of an approximation must not do.

cf_leaf <- function(psi_soil = 0.5, PPFD = 1500, atm_vpd = 1.0,
                    leaf_temp = 25.0, method = "closed", curve = "TF24",
                    lambda = 1.5e5) {
  l <- leaf_model(supply = leaf_supply_singlelayer())
  l$CF77_lambda_ <- lambda
  set_drivers(l, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd,
              leaf_temp = leaf_temp,
              root_network = series_resistance(1500))
  l$set_model(curve, "stem", method)
  l
}

test_that("the method defaults to exact and a two-argument seat works", {
  # The compatibility assertion, and the reason the third argument carries an
  # R-side default rather than being required: every existing call site, every
  # vignette and README block, and `.seat_model()` itself pass two arguments.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  expect_identical(l$model_method(), "exact")

  l$set_model("TF24", "stem")
  expect_identical(l$model_method(), "exact")
  expect_identical(l$model_route(), "stem")
  expect_identical(l$model_curve(), "TF24")

  # And it is settable back and forth, because it is configuration like the
  # other two axes rather than a property of a solve.
  l$set_model("TF24", "stem", "closed")
  expect_identical(l$model_method(), "closed")
  l$set_model("TF24", "stem")
  expect_identical(l$model_method(), "exact")

  expect_error(l$set_model("TF24", "stem", "fast"), "unknown method")
})

test_that("seating exact is bit-identical to not seating a method at all", {
  # THE assertion that makes this change safe. `optimise()` returns to its
  # previous body when the method is exact, so this is bit-identity by
  # construction -- and asserting it here is what would catch someone "tidying"
  # the early return into a branch that reorders the arithmetic.
  a <- cf_leaf(method = "exact")
  a$optimise()
  b <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(b, psi_soil = 0.5, PPFD = 1500, atm_vpd = 1.0,
              root_network = series_resistance(1500))
  b$set_model("TF24", "stem")
  b$optimise()

  expect_identical(a$opt_psi_stem_, b$opt_psi_stem_)
  expect_identical(a$ci_, b$ci_)
  expect_identical(a$assim_colimited_, b$assim_colimited_)
  expect_identical(a$transpiration_, b$transpiration_)
  expect_identical(a$stom_cond_CO2_, b$stom_cond_CO2_)
  expect_identical(a$profit_, b$profit_)
  expect_identical(a$lambda_emergent, b$lambda_emergent)
})

test_that("the closed form writes every output the exact solve writes", {
  # Hazard 8: `Leaf` is a value member that plant reuses for every individual in
  # a patch, so an output a branch declines to write silently becomes the
  # PREVIOUS solve's value. The golden files cannot see this class of defect
  # because they build a fresh leaf per point, so the leaf is REUSED here and
  # the outputs are poisoned first.
  l <- cf_leaf()
  l$set_model("TF24", "collar")
  l$optimise()                       # seat a collar solve's outputs
  poison <- -12345
  for (f in c("ci_", "assim_colimited_", "transpiration_", "stom_cond_CO2_",
              "hydraulic_cost_", "profit_", "opt_psi_stem_")) {
    l[[f]] <- poison
  }
  l$set_model("TF24", "stem", "closed")
  l$optimise()

  expect_false(l$last_solve_fell_back_)
  for (f in c("ci_", "assim_colimited_", "transpiration_", "stom_cond_CO2_",
              "hydraulic_cost_", "profit_", "opt_psi_stem_")) {
    expect_false(identical(l[[f]], poison), info = f)
    expect_true(is.finite(l[[f]]), info = f)
  }
  # The collar-route outputs describe a different solve and are cleared, exactly
  # as the exact stem route clears them.
  expect_false(is.finite(l$opt_root_psi_))
  expect_false(is.finite(l$E_up_))

  # The identity that ties the written fields to each other rather than to a
  # recorded number: profit IS the objective at the reported potential. A path
  # that took psi_stem from the closed form and the rates from anywhere else
  # would fail this.
  expect_equal(l$profit_, l$assim_colimited_ - l$hydraulic_cost_,
               tolerance = 1e-12)
})

test_that("a fallback IS the exact solve, bit for bit", {
  # The guard tests an OUTPUT (ci/ca > 0.5, and strictly inside the bracket), so
  # it can only be applied after solving. What makes that acceptable is that a
  # rejected row costs accuracy nothing at all -- not "is close to the exact
  # answer", but is it.
  dry <- 4.5
  a <- cf_leaf(psi_soil = dry, method = "exact")
  a$optimise()
  b <- cf_leaf(psi_soil = dry, method = "closed")
  b$optimise()

  expect_true(b$last_solve_fell_back_)
  expect_identical(b$opt_psi_stem_, a$opt_psi_stem_)
  expect_identical(b$assim_colimited_, a$assim_colimited_)
  expect_identical(b$profit_, a$profit_)
  expect_identical(b$closed_form_fallback_fraction(), 1)
})

test_that("the fallback fraction is measurable, which it was not before", {
  # `closed_form.hpp`'s open item 3: the realised speedup is
  # 1/[phi + (1-phi)/S], and phi had never been measured on a real water-limited
  # scenario because nothing counted it. One leaf, one sweep from wet to dry.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  # NA, not 0: "never asked" and "never fell back" are different answers.
  expect_true(is.na(l$closed_form_fallback_fraction()))

  l$set_model("TF24", "stem", "closed")
  soils <- seq(0.25, 5.0, length.out = 20)
  for (p in soils) {
    set_drivers(l, psi_soil = p, PPFD = 1500, atm_vpd = 1.0,
                root_network = series_resistance(1500))
    l$optimise()
  }
  expect_identical(l$closed_form_calls_, 20L)
  phi <- l$closed_form_fallback_fraction()
  # Both bounds matter. A phi of 0 would mean the guard never fires on a sweep
  # that reaches psi_crit, i.e. the guard is not doing its job; a phi of 1 would
  # mean the fast path never runs and the method is pure overhead.
  expect_gt(phi, 0)
  expect_lt(phi, 1)
  expect_identical(phi, l$closed_form_fallbacks_ / l$closed_form_calls_)

  l$reset_closed_form_counters()
  expect_identical(l$closed_form_calls_, 0L)
  expect_true(is.na(l$closed_form_fallback_fraction()))
})

test_that("CF77 is the n = 0 arm, so its lambda is the price exactly", {
  # The one analytical identity available at n = 0: with lambda constant there
  # is no power law to approximate, so the reported marginal cost must be the
  # prescribed price itself, to the last bit. This is also the check that the
  # closed arm reads `cf77_price()` rather than reconstructing a lambda -- which
  # is what lets it compose with the soil-moisture option (#128).
  lambda <- 1.5e5
  l <- cf_leaf(curve = "CF77", lambda = lambda)
  l$optimise()
  expect_false(l$last_solve_fell_back_)
  expect_identical(l$lambda_emergent, lambda)

  # With the Medlyn soil-moisture factor on, the price is lambda/beta and the
  # closed form must see the divided price, not the bare field.
  theta <- 0.35
  m <- leaf_model(supply = leaf_supply_singlelayer())
  m$CF77_lambda_ <- lambda
  m$CF77_soil_beta_ <- TRUE
  m$theta <- theta
  set_drivers(m, psi_soil = 0.5, PPFD = 1500, atm_vpd = 1.0,
              root_network = series_resistance(1500))
  m$set_model("CF77", "stem", "closed")
  m$optimise()
  beta <- (theta - m$theta_w_) / (m$theta_fc_ - m$theta_w_)
  expect_equal(m$lambda_emergent, lambda / beta, tolerance = 1e-12)
  # And the behaviour that follows: a dearer price closes the stomata.
  expect_lt(m$opt_psi_stem_, l$opt_psi_stem_)
})

test_that("the reported potential is the one whose supply carries its E", {
  # The supply inversion, checked as a round trip rather than restated. ⚠️ It is
  # `psi = G^-1(E/kmax + G(psi_up))` and NOT `E/kmax`: the supply here is the
  # vulnerability curve's INTEGRAL, so the two agree only on a leaf with no
  # vulnerability curve. The header used to claim the latter and that claim is
  # withdrawn -- this is the assertion that keeps it withdrawn.
  l <- cf_leaf(curve = "CF77")
  l$optimise()
  psi <- l$opt_psi_stem_
  psi_soil <- l$single_psi_soil_

  expect_gt(psi, psi_soil)
  expect_lt(psi, l$psi_crit)
  expect_equal(l$transpiration(psi, psi_soil), l$transpiration_,
               tolerance = 1e-10)
  # And the naive form really is different, by a lot -- so this is a test and
  # not a tautology.
  naive <- l$transpiration_ / l$leaf_specific_conductance_max_
  expect_gt(abs(naive - psi) / psi, 0.1)
})

test_that("dA/dE equals lambda at the closed form's own operating point", {
  # The first-order condition the whole family shares, differenced against the
  # model's own A and E at the potential the closed form returned. It does NOT
  # hold exactly -- the closed form approximates the argmax, so its point is not
  # stationary -- and the tolerance is what states how far off it is. That is
  # the honest form of this assertion: an exact-equality version would be a
  # claim about the method that the method does not make.
  l <- cf_leaf(curve = "CF77")
  l$optimise()
  psi <- l$opt_psi_stem_
  expect_false(l$last_solve_fell_back_)

  h <- 1e-5
  at <- function(p) {
    l$evaluate_psi_stem_at(p)
    c(A = l$assim_colimited_, E = l$transpiration_)
  }
  hi <- at(psi + h)
  lo <- at(psi - h)
  fd <- (hi[["A"]] - lo[["A"]]) / (hi[["E"]] - lo[["E"]])
  # Within 10 %: measured 3 % here. The exact solve satisfies the same identity
  # to 1e-4 (test-cf77-soil-beta.R), and the gap between those two figures IS
  # the closed form's error.
  expect_equal(fd / 1.5e5, 1, tolerance = 0.1)
})

test_that("accuracy against the exact solve degrades as ci/ca falls", {
  # THE ACCURACY TEST, and it is a statement about the SHAPE of the error rather
  # than about any one number. Measured over a driver grid at 25 C, binned by
  # the exact solve's ci/ca:
  #
  #   ci/ca      mean |dA|/A   worst
  #   > 0.9      0.20 %        0.34 %
  #   0.8-0.9    0.26 %        0.68 %
  #   0.7-0.8    0.45 %        1.4 %
  #   0.6-0.7    0.81 %        3.6 %
  #   0.5-0.6    2.2 %         6.5 %
  #   < 0.4      falls back, so 0 by construction
  #
  # ⚠️ The bounds asserted below are an ORDER OF MAGNITUDE looser than those
  # figures, deliberately. What must not silently change is that the error is
  # small at the wet end, grows monotonically as ci/ca falls, and is exactly
  # zero once the guard fires. Pinning 0.20 % would make this a snapshot of the
  # approximation instead of a test of it.
  grid <- expand.grid(psi_soil = c(0.25, 0.5, 1, 2, 3, 4, 5),
                      PPFD = c(100, 500, 1500),
                      atm_vpd = c(0.5, 1, 2, 4))

  rows <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    ex <- cf_leaf(g$psi_soil, g$PPFD, g$atm_vpd, method = "exact")
    ex$optimise()
    cl <- cf_leaf(g$psi_soil, g$PPFD, g$atm_vpd, method = "closed")
    cl$optimise()
    data.frame(frac = ex$ci_ / ex$ca_,
               fell = cl$last_solve_fell_back_,
               # Relative in A with a floor of 1 umol. ⚠️ A bare relative error
               # divides by A, and a dim hot leaf has A ~ 1e-5, so the metric is
               # then amplified by the cancellation and says nothing about the
               # method -- the same trap the golden file's cross-platform
               # figures fell into.
               err = abs(cl$assim_colimited_ - ex$assim_colimited_) /
                 max(1, abs(ex$assim_colimited_)),
               dpsi = abs(cl$opt_psi_stem_ - ex$opt_psi_stem_))
  })
  res <- do.call(rbind, rows)
  res <- res[is.finite(res$frac) & is.finite(res$err), ]
  expect_gt(nrow(res), 50)

  # 1. Where the guard fires, the answer is the exact one, so the error is zero
  #    and not merely small.
  expect_true(all(res$err[res$fell] == 0))
  expect_true(all(res$dpsi[res$fell] == 0))
  expect_gt(sum(res$fell), 0)

  # 2. On the fast path, the wet end is accurate.
  wet <- res[!res$fell & res$frac > 0.8, ]
  expect_gt(nrow(wet), 5)
  expect_lt(max(wet$err), 0.02)

  # 3. And it degrades monotonically in the bin means, which is the property the
  # ci/ca guard is built on. Monotone in the MEANS rather than row by row: the
  # error is a function of more than ci/ca, so individual rows cross over.
  fast <- res[!res$fell, ]
  bins <- cut(fast$frac, breaks = c(0.4, 0.6, 0.7, 0.8, 1.01))
  means <- tapply(fast$err, bins, mean)
  means <- means[!is.na(means)]
  expect_gt(length(means), 2)
  # `means` runs dry -> wet, so it must be DECREASING.
  expect_true(all(diff(means) <= 0))

  # 4. The guard admits nothing catastrophic. It is a filter on gross failure
  # and not an error bound -- saying which is the point of the bound's size.
  expect_lt(max(fast$err), 0.15)
})

test_that("everything the closed form cannot do is refused, with the reason", {
  # ⚠️ THE TWO CURVE REFUSALS ARE DIFFERENT CLAIMS, and that distinction is
  # going into the paper as the existence criterion, so the messages are
  # asserted rather than just the throwing. The closed form exists iff h'(A)
  # does not depend on the solution: it is 1 for the identity link and 1/|A|max
  # for ProfitMax, whose normaliser is a per-driver-set constant -- so those
  # four are FEASIBLE and only a dlambda/dpsi is unwritten. SOX and JW26 have
  # the log link, h' = 1/A, so lambda carries the assimilation the solve is for;
  # that is a fixed point and there is nothing to write.
  for (curve in c("JS22", "CMax", "ProfitMax")) {
    l <- cf_leaf(curve = curve)
    expect_error(l$optimise(), "not implemented", info = curve)
    expect_error(l$optimise(), "EXISTS for this curve", info = curve)
  }
  for (curve in c("SOX", "JW26")) {
    l <- cf_leaf(curve = curve)
    expect_error(l$optimise(), "does not exist", info = curve)
    expect_error(l$optimise(), "LOG benefit link", info = curve)
  }

  # The collar route is a different model rather than a restriction of this one
  # -- it optimises the soil-to-collar path, which this inversion says nothing
  # about.
  collar <- cf_leaf()
  collar$set_model("TF24", "collar", "closed")
  expect_error(collar$optimise(), "STEM-route method")

  # Multi-layer: the closed form for it does not exist, so it is refused rather
  # than run on layer zero.
  ml <- leaf_model()
  set_drivers(ml, psi_soil = c(0.5, 0.75, 1.0), soil_depth = c(1, 2, 3),
              root_network = root_network_from_carbon(rep(1 / 3, 3),
                                                   c(1, 2, 3)))
  ml$set_model("TF24", "stem", "closed")
  expect_error(ml$optimise(), "single-layer only")

  # #116: with the energy balance on, the deficit Fick's law divides by depends
  # on Tleaf, which depends on E, which is what the inversion solves for. So D
  # is not known at the moment of inversion and the form is not closed -- which
  # is a statement about the FORM, not about its accuracy.
  eb <- cf_leaf()
  eb$use_energy_balance_ <- TRUE
  set_drivers(eb, psi_soil = 0.5, PPFD = 1500, atm_vpd = 1.0,
              root_network = series_resistance(1500))
  eb$set_model("TF24", "stem", "closed")
  expect_error(eb$optimise(), "not closed")
  # And the exact solve on the same leaf is untouched, so the refusal is the
  # method's and not the leaf's.
  eb$set_model("TF24", "stem", "exact")
  expect_silent(eb$optimise())
})

test_that("the gradient differentiates the exact solve whatever is seated", {
  # ⚠️ Seating a model resets the method to exact, so `gradient::route_seat` --
  # which seats once per observation -- always differentiates the exact solve.
  # That is deliberate: the implicit-function-theorem composite is derived from
  # the exact first-order condition, and the closed form's argmax is
  # DISCONTINUOUS across its own guard boundary (measured: the mean |second
  # difference| of the argmax over a kmax sweep rises 100x where the sweep
  # crosses it). A gradient through that is not a gradient. This test pins the
  # silent part.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  l$set_model("TF24", "stem", "closed")
  expect_identical(l$model_method(), "closed")

  g <- leaf_gradient(psi_soil = 0.5, PPFD = 1500, x = l,
                     traits = leaf_traits(), model = "TF24",
                     root_network = series_resistance(1500),
                     pars = "vcmax_25")
  # `collar` is NaN on a stem route by construction, so the column that says the
  # gradient ran is `A`.
  expect_true(is.finite(g$gradient[1, "A"]))
  expect_identical(g$status, "interior")
  # And the seat is back to exact, which is the silent part being pinned.
  expect_identical(l$model_method(), "exact")
})
