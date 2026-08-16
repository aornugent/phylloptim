# The hydraulic cost's own derivatives, which ARE the trait row of a zero-flux
# operating point.
#
# Under shutdown or at the compensation point the leaf pays only respiration and
# this cost -- profit is -R_d - C(psi_crit), reading neither soil nor light -- so
# every environment row there is exactly zero and these five numbers are the
# whole of the trait row.
#
# The kernel used to take only psi as its scalar argument, reading stem_b and
# stem_c off the object, so forward mode could reach d/dpsi and nothing else. The
# overload takes them as arguments, which is what makes the rest available at
# all.

cost_traits <- function() leaf_traits()

# C(psi) = cost_scale * (1 - proportion_of_conductivity(psi))^beta2, rebuilt from
# the traits each time. Shares no code with the forward-mode row: this evaluates
# the curve, that differentiates it.
cost_value <- function(tr, psi) {
  l <- leaf_model(traits = tr)
  tr$cost_scale_TF24 * (1 - l$proportion_of_conductivity(psi))^tr$beta2
}

cost_row <- function(tr, psi) {
  leaf_model(traits = tr)$hydraulic_cost_row_values(psi)
}

# The row's layout, and the trait each entry belongs to. Position IS the
# interface, so naming it here is what a reordering would break.
cost_row_traits <- c("stem_b", "stem_c", "beta2", "cost_scale_TF24")

test_that("the cost row agrees with a rebuilt difference in every direction", {
  tr <- cost_traits()
  worst <- 0
  for (psi in c(1.0, 3.0, 5.0, 6.5)) {
    r <- cost_row(tr, psi)
    expect_equal(r[[1]], 1)

    h <- 1e-6
    fd_psi <- (cost_value(tr, psi + h) - cost_value(tr, psi - h)) / (2 * h)
    worst <- max(worst, abs(r[[2]] / fd_psi - 1))

    for (k in seq_along(cost_row_traits)) {
      nm <- cost_row_traits[[k]]
      hh <- max(abs(tr[[nm]]), 1) * 1e-6
      up <- tr; up[[nm]] <- up[[nm]] + hh
      dn <- tr; dn[[nm]] <- dn[[nm]] - hh
      fd <- (cost_value(up, psi) - cost_value(dn, psi)) / (2 * hh)
      worst <- max(worst, abs(r[[2 + k]] / fd - 1))
    }
  }
  message(sprintf("  cost row vs rebuilt difference: worst %.2e", worst))
  expect_lt(worst, 1e-7)
})

test_that("the cost row is not vacuously zero, and stem_c changes sign", {
  tr <- cost_traits()
  # Every entry live somewhere: a row of zeros would pass a ratio check against
  # a difference that is also zero.
  mid <- cost_row(tr, 3.0)
  for (k in 2:6) expect_true(abs(mid[[k]]) > 1e-6)

  # The curve is sigmoid, so its steepness pulls the cost in opposite directions
  # either side of the inflexion. That is a property of the model, and a row that
  # got the sign from the value rather than the derivative would not have it.
  below <- cost_row(tr, 3.0)[[4]]
  above <- cost_row(tr, 5.0)[[4]]
  message(sprintf("  d(cost)/d(stem_c): %.4f at psi 3, %.4f at psi 5", below, above))
  expect_lt(below, 0)
  expect_gt(above, 0)
})

test_that("the cost scale enters as a pure prefactor", {
  # C is linear in cost_scale_TF24, so its derivative is C over that scale
  # exactly -- an identity the row must satisfy without being told.
  tr <- cost_traits()
  for (psi in c(1.0, 3.0, 5.0)) {
    r <- cost_row(tr, psi)
    expect_equal(r[[6]], cost_value(tr, psi) / tr$cost_scale_TF24,
                 tolerance = 1e-12)
  }
})
