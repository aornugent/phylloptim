# A feasibility bound's own derivative, by the implicit function theorem on the
# residual that defines it rather than by differentiating the search that found
# it.
#
# The referee is a central difference of find_root_psi itself, which shares no
# code with the row: one solves the residual, the other differentiates it.
#
# ⚠️ THE SOIL PROFILE MUST NOT BE UNIFORM, and that is a property of the wet
# bound rather than a convenience. The wet bound is the collar at which the
# per-layer fluxes SUM TO ZERO, and on a uniform profile that lands essentially
# on top of every layer's own potential -- which is exactly where each layer's
# span |x - psi_i| has its kink. A central difference there straddles the kink
# while the analytic row takes the side the model is on, so the two disagree by
# percent and neither is wrong. Measured: 5.5e-02 on a uniform profile against
# 3.6e-07 on a mild gradient and 4.4e-05 on a strong one.
#
# A drying stand leaves uniformity within weeks -- the profile drains -- so the
# gradient cases below are the states this row is actually taken at.

bound_row_leaf <- function(psi_soil) {
  l <- leaf_model(traits = leaf_traits())
  set_drivers(l, psi_soil = psi_soil)
  l
}

# find_root_psi's own answer, re-solved from scratch at a perturbed profile.
bound_position <- function(psi_soil, which) {
  bound_row_leaf(psi_soil)$find_root_psi(
    min(psi_soil), psi_soil, if (which == 0L) 0L else 1L)
}

bound_row_of <- function(psi_soil, which) {
  bound_row_leaf(psi_soil)$bound_row_values(which)
}

# Layout is documented at bound_row_values; entry 7 + j is layer j's soil row.
bound_row_worst <- function(psi_soil, which, h = 1e-5) {
  r <- bound_row_of(psi_soil, which)
  testthat::expect_equal(r[[1]], 1)          # the row is finite
  worst <- 0
  for (j in seq_along(psi_soil)) {
    up <- psi_soil; up[[j]] <- up[[j]] + h
    dn <- psi_soil; dn[[j]] <- dn[[j]] - h
    fd <- (bound_position(up, which) - bound_position(dn, which)) / (2 * h)
    if (abs(fd) > 1e-12) worst <- max(worst, abs(r[[7 + j]] / fd - 1))
  }
  worst
}

mild_gradient <- c(2.6, 2.8, 3.0, 3.2, 3.4)
strong_gradient <- c(1.0, 2.0, 3.0, 4.0, 5.0)

test_that("the wet bound's row agrees with a differenced solve", {
  # Residual is total uptake, so the row has NO stem terms at all.
  for (profile in list(mild_gradient, strong_gradient)) {
    worst <- bound_row_worst(profile, 0L)
    message(sprintf("  wet bound, psi %s: worst %.2e",
                    paste(range(profile), collapse = "-"), worst))
    expect_lt(worst, 1e-4)
  }
  r <- bound_row_of(mild_gradient, 0L)
  expect_equal(r[[4]], 0)   # d/d(kappa)
  expect_equal(r[[5]], 0)   # d/d(psi_crit)
  expect_equal(r[[7]], 0)   # d/d(stem_b)
  expect_gt(r[[3]], 0)      # the residual slope is a conductance
})

test_that("the dry bound's row carries the stem terms and agrees too", {
  # Residual is E_up - kappa*[G(psi_crit) - G(x)], so the stem enters three ways
  # and dR/dstem_b comes from the homogeneity identity rather than a rebuild.
  for (profile in list(mild_gradient, strong_gradient)) {
    worst <- bound_row_worst(profile, 1L)
    message(sprintf("  dry bound, psi %s: worst %.2e",
                    paste(range(profile), collapse = "-"), worst))
    expect_lt(worst, 1e-4)
  }
  r <- bound_row_of(mild_gradient, 1L)
  # Non-vacuity: these are the entries the wet row does not have, so a dry row
  # that left them at zero would pass every check above.
  expect_true(r[[4]] != 0)
  expect_true(r[[5]] != 0)
  expect_true(r[[7]] != 0)
  # dR/dx is dE_up/dx + kappa*G'(x), a sum of two strictly positive terms, so it
  # is larger than the wet bound's slope rather than merely different.
  expect_gt(r[[3]], bound_row_of(mild_gradient, 0L)[[3]])
})

test_that("the root's own critical potential is an exact unit row", {
  # This arm's bound is a registered constant, so it moves with nothing but
  # itself. No solve, no differencing -- the row is exact by inspection, and it
  # is the cheapest of the three.
  r <- bound_row_of(mild_gradient, 2L)
  expect_equal(r[[1]], 1)
  expect_equal(r[[6]], -1)     # d/d(root_psi_crit)
  expect_equal(r[[3]], 1)      # residual slope
  expect_equal(r[[4]], 0)
  expect_equal(r[[5]], 0)
  expect_equal(r[[7]], 0)
  expect_true(all(r[-(1:7)] == 0))
})

test_that("a uniform profile cannot referee the wet bound", {
  # The fixture requirement stated as a measurement rather than a comment. If
  # this ever stops holding, the kink has moved and the gradient cases above
  # need re-reading -- not this one relaxing.
  uniform <- rep(3.0, 5)
  worst <- bound_row_worst(uniform, 0L)
  message(sprintf("  wet bound on a uniform profile: worst %.2e", worst))
  expect_gt(worst, 1e-3)
  # And the dry bound is unaffected, because its own bound sits well away from
  # the soil potentials -- which is what identifies the cause as the kink.
  expect_lt(bound_row_worst(uniform, 1L), 1e-4)
})

test_that("the bound row's parameter half agrees with a rebuilt difference", {
  # The soil rows above are only part of the row. These three are what a pinned
  # point's kappa, psi_crit and stem_b columns are built from, and consuming them
  # unrefereed would be building on numbers nothing had checked.
  base <- leaf_traits()
  kmax0 <- 3.14e-05
  seat <- function(tr, kmax = kmax0) {
    l <- leaf_model(traits = tr)
    set_drivers(l, psi_soil = mild_gradient, leaf_specific_conductance_max = kmax)
    l
  }
  bound_of <- function(tr, which, kmax = kmax0) {
    seat(tr, kmax)$find_root_psi(min(mild_gradient), mild_gradient,
                                 if (which == 0L) 0L else 1L)
  }

  # Entry 4 is kappa, 5 is the stem's psi_crit, 7 is stem_b.
  dry <- seat(base)$bound_row_values(1L)
  for (e in list(list(5, "psi_crit"), list(7, "stem_b"))) {
    h <- max(abs(base[[e[[2]]]]), 1) * 1e-6
    up <- base; up[[e[[2]]]] <- up[[e[[2]]]] + h
    dn <- base; dn[[e[[2]]]] <- dn[[e[[2]]]] - h
    fd <- (bound_of(up, 1L) - bound_of(dn, 1L)) / (2 * h)
    expect_equal(dry[[e[[1]]]] / fd, 1, tolerance = 1e-4)
  }
  hk <- kmax0 * 1e-4
  fd_k <- (bound_of(base, 1L, kmax0 + hk) - bound_of(base, 1L, kmax0 - hk)) / (2 * hk)
  expect_equal(dry[[4]] / fd_k, 1, tolerance = 1e-4)

  # And the wet bound has none of them, which is the derivation rather than a
  # coincidence: its residual is total uptake and the stem is not in it.
  wet <- seat(base)$bound_row_values(0L)
  expect_equal(wet[[4]], 0)
  expect_equal(wet[[5]], 0)
  expect_equal(wet[[7]], 0)
})
