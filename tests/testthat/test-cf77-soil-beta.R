# The optional Medlyn (2011) soil-moisture factor on the Cowan-Farquhar price.
#
# WHY THE OPTION EXISTS. CF77 is the one curve here whose PRICE of water is
# invariant to soil state by construction: cost_deriv<CF77> discards
# psi_upstream and lambda_for<CF77> returns the bare constant, so soil water
# reaches the optimum only through the feasible bracket and the cost's VALUE.
# That leaves it without a soil-moisture shutdown, which is the one thing a
# Medlyn-style model is usually reached for. `CF77_soil_beta_` supplies one by
# dividing the price by beta = (theta - theta_w)/(theta_fc - theta_w), and does
# it as a FIELD rather than an eighth CostCurve so the enum, the gradient
# parameter indices, R's mirror enumeration and the 4608-row golden file are all
# untouched.
#
# The tests below are in the order the change has to survive them: off is
# bit-identical, on reduces to off in the wet limit, the first-order condition
# still holds with the modified price, and the behaviour it was added for
# actually happens.

cf77_leaf <- function(theta = NULL, beta_on = FALSE, lambda = 1.5e5,
                      psi_soil = 0.5) {
  l <- leaf_model(supply = leaf_supply_singlelayer())
  # In band: marginal_cost_water runs 9e4-3e5 at package defaults, and a lambda
  # far outside it pins the optimum shut or wide open regardless of beta.
  l$CF77_lambda_ <- lambda
  l$CF77_soil_beta_ <- beta_on
  if (!is.null(theta)) l$theta <- theta
  set_drivers(l, psi_soil = psi_soil, root_network = series_resistance(1500))
  l$set_model("CF77", "stem")
  l
}

test_that("the option is off by default and off is bit-identical", {
  # The golden files pin CF77 rows, so this is the test that says the change is
  # additive. cf77_price() returns the member untouched on the off path rather
  # than dividing it by a computed 1.0 -- if someone "simplifies" that to a
  # single expression, this is what should fail.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  expect_false(l$CF77_soil_beta_)

  off <- cf77_leaf(beta_on = FALSE)
  off$optimise()

  # Same leaf, same drivers, option explicitly off and theta set to a value
  # the factor would have reacted to. Bit-identical, not merely close.
  also_off <- cf77_leaf(theta = 0.25, beta_on = FALSE)
  also_off$optimise()

  expect_identical(also_off$opt_psi_stem_, off$opt_psi_stem_)
  expect_identical(also_off$stom_cond_CO2_, off$stom_cond_CO2_)
  expect_identical(also_off$assim_colimited_, off$assim_colimited_)
  expect_identical(also_off$lambda_emergent, off$lambda_emergent)
})

test_that("at or above field capacity it reduces to plain CF77, exactly", {
  off <- cf77_leaf(beta_on = FALSE)
  off$optimise()

  # beta == 1 exactly at theta == theta_fc.
  at_fc <- cf77_leaf(theta = 0.5, beta_on = TRUE)
  at_fc$optimise()
  expect_identical(at_fc$opt_psi_stem_, off$opt_psi_stem_)
  expect_identical(at_fc$lambda_emergent, off$lambda_emergent)

  # And above it, because beta is CLAMPED at 1. This is the behaviour that makes
  # CF77_lambda_ a floor the wet end cannot undercut, and it is a choice rather
  # than an accident: without the clamp, saturated soil would price water
  # BELOW the prescribed marginal value.
  above_fc <- cf77_leaf(theta = 0.8, beta_on = TRUE)
  above_fc$optimise()
  expect_identical(above_fc$opt_psi_stem_, off$opt_psi_stem_)
  expect_identical(above_fc$lambda_emergent, off$lambda_emergent)
})

test_that("dA/dE equals lambda/beta at the optimum", {
  # The analytical identity, and the reason this is a cost curve rather than an
  # empirical gs formula: CF77's first-order condition is dA/dE == lambda, and
  # with the factor on the prescribed value becomes lambda/beta. Differenced
  # against the model's own A and E, so nothing here trusts lambda_emergent.
  theta <- 0.35
  beta <- (theta - 0.2) / (0.5 - 0.2)
  lambda <- 1.5e5

  l <- cf77_leaf(theta = theta, beta_on = TRUE, lambda = lambda)
  l$optimise()
  psi <- l$opt_psi_stem_

  # Interior optimum, or the identity is a statement about a bound instead.
  expect_gt(psi, l$single_psi_soil_ + 1e-6)
  expect_lt(psi, l$psi_crit - 1e-6)

  h <- 1e-5
  at <- function(p) {
    l$evaluate_psi_stem_at(p)
    c(A = l$assim_colimited_, E = l$transpiration_)
  }
  hi <- at(psi + h)
  lo <- at(psi - h)
  fd_dassim_dtransp <- (hi[["A"]] - lo[["A"]]) / (hi[["E"]] - lo[["E"]])

  expect_equal(fd_dassim_dtransp, lambda / beta, tolerance = 1e-4)
  # And the reported lambda agrees with the prescribed price, which is what puts
  # this curve on the same axis as the others.
  expect_equal(l$lambda_emergent, lambda / beta, tolerance = 1e-10)
})

test_that("drying the soil closes stomata, which plain CF77 will not do", {
  thetas <- c(0.50, 0.45, 0.40, 0.35, 0.30)

  solve_at <- function(theta, beta_on) {
    l <- cf77_leaf(theta = theta, beta_on = beta_on)
    l$optimise()
    c(gs = l$stom_cond_CO2_, psi = l$opt_psi_stem_)
  }

  on <- vapply(thetas, solve_at, numeric(2), beta_on = TRUE)
  gs_on <- on["gs", ]

  # Monotone, and not by a rounding margin: the wettest point carries several
  # times the conductance of the driest.
  expect_true(all(diff(gs_on) < 0))
  expect_gt(gs_on[1] / gs_on[length(gs_on)], 2)

  # Drier still and the leaf shuts: beta -> 0 sends the price to infinity, which
  # is the shutdown CF77 structurally lacks.
  shut <- solve_at(0.22, beta_on = TRUE)
  expect_lt(shut[["gs"]], 1e-3)

  # THE CONTROL, and the point of the whole exercise: with the option off, the
  # same sweep does nothing at all, because soil moisture is not an argument to
  # the bare constant.
  off <- vapply(thetas, solve_at, numeric(2), beta_on = FALSE)
  expect_identical(length(unique(off["gs", ])), 1L)
})

test_that("an unusable soil-moisture factor is refused, not optimised", {
  # beta <= 0 is theta at or below the wilting point. The price is then infinite
  # or NEGATIVE, and a negative price pays the leaf to transpire -- which the
  # optimiser takes to the wet bound, returning a bracket property that looks
  # like an operating point. Same reasoning as the NA CF77_lambda_ refusal.
  at_wilting <- cf77_leaf(theta = 0.2, beta_on = TRUE)
  expect_error(at_wilting$optimise(), "beta")

  below_wilting <- cf77_leaf(theta = 0.15, beta_on = TRUE)
  expect_error(below_wilting$optimise(), "must exceed theta_w_")

  # A degenerate range is refused separately, because there the division
  # itself fails rather than the sign of its result.
  degenerate <- cf77_leaf(theta = 0.3, beta_on = TRUE)
  degenerate$theta_fc <- 0.2
  set_drivers(degenerate, psi_soil = 0.5,
              root_network = series_resistance(1500))
  degenerate$set_model("CF77", "stem")
  expect_error(degenerate$optimise(), "theta_fc_ > theta_w_")

  # ⚠️ And none of this fires when the option is off, so a caller who never
  # touches it cannot be refused for values the bare constant never reads.
  ok <- cf77_leaf(theta = 0.15, beta_on = FALSE)
  expect_silent(ok$optimise())
})

test_that("the soil-moisture inputs survive re-driving", {
  # The #96-shaped hazard, and the reason `theta` is exposed and `theta_` is
  # not: set_physiology overwrites all three working copies from the inputs on
  # every driver set, so a bare `l$theta_ <- x` is silently undone by the next
  # set_drivers() and a sweep would run every point at the default.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  l$theta <- 0.25
  l$theta_fc <- 0.45
  l$theta_w <- 0.15

  set_drivers(l, psi_soil = 0.5, root_network = series_resistance(1500))
  expect_identical(l$theta, 0.25)
  expect_identical(l$theta_, 0.25)      # propagated to the working copy
  expect_identical(l$theta_fc_, 0.45)
  expect_identical(l$theta_w_, 0.15)

  # Still there after a re-trait and a second drive, which is the sequence a
  # sweep actually runs.
  set_traits(l, leaf_traits(vcmax_25 = 120))
  set_drivers(l, psi_soil = 0.5, root_network = series_resistance(1500))
  expect_identical(l$theta_, 0.25)

  # The working copy is derived state and IS wiped, so the input is what carries
  # the value across. Asserting this keeps the two roles distinct.
  expect_identical(l$theta, 0.25)
})

test_that("the option touches CF77 and nothing else", {
  # It reads theta_, which every curve now has, so the guard against it leaking
  # is that no other cost curve consults cf77_price(). TF24 at the same drivers
  # must not notice.
  solve_tf24 <- function(beta_on, theta) {
    l <- leaf_model(supply = leaf_supply_singlelayer())
    l$CF77_soil_beta_ <- beta_on
    l$theta <- theta
    set_drivers(l, psi_soil = 0.5, root_network = series_resistance(1500))
    l$set_model("TF24", "stem")
    l$optimise()
    l$opt_psi_stem_
  }
  expect_identical(solve_tf24(TRUE, 0.25), solve_tf24(FALSE, 0.5))
})
