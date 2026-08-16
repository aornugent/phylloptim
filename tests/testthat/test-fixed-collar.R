# Two consumers want a collar potential treated in opposite ways, and one entry
# point used to serve both.
#
# evaluate_root_collar_psi CLAMPS the target into the feasible interval. That is
# right for a tracked collar being nudged -- a target outside the interval means
# the state has drifted out, and the projection back in is the model. It is wrong
# for a partial derivative taken at a FROZEN collar, where the clamp replaces the
# collar the caller asked about with a different one and says nothing.
#
# profit_at_fixed_collar does not clamp and does not project: a collar outside
# the interval comes back infeasible and IS NOT EVALUATED, because below the wet
# bound the profit algebra runs on a negative conductance and returns a plausible
# number.

fixed_collar_seated <- function(psi, traits = leaf_traits()) {
  l <- leaf_model(traits = traits)
  set_drivers(l, psi_soil = rep(psi, 5))
  l$find_root_collar_psi()
  l
}

# One finite-difference arm, exactly as the stand's gradient takes it: solve,
# hold the collar, move a trait, re-evaluate at the held collar.
fixed_collar_arm <- function(psi, name, side, rel = 1e-3) {
  base <- leaf_traits()
  seated <- fixed_collar_seated(psi, base)
  collar <- seated$opt_root_psi_

  tr <- base
  tr[[name]] <- base[[name]] + side * max(abs(base[[name]]), 1) * rel
  m <- leaf_model(traits = tr)
  set_drivers(m, psi_soil = rep(psi, 5))
  v <- m$profit_at_fixed_collar_values(collar)
  list(collar = collar,
       bound_a = seated$find_root_psi(psi, rep(psi, 5), 0L),
       moved_bound_a = m$find_root_psi(psi, rep(psi, 5), 0L),
       feasible = v[[1]] == 1,
       profit = v[[2]],
       clamped = m$evaluate_root_collar_psi(collar),
       clamped_at = m$opt_root_psi_)
}

test_that("a held collar well inside the interval is evaluated where it is", {
  # Wet: the collar sits far from either bound, so a trait step cannot move a
  # bound across it and the two routes must agree exactly.
  for (side in c(1, -1)) {
    a <- fixed_collar_arm(2.0, "root_b", side)
    expect_gt(a$collar - a$bound_a, 0.1)
    expect_true(a$feasible)
    expect_equal(a$profit, a$clamped)
    expect_equal(a$clamped_at, a$collar)
  }
})

test_that("a trait step across the wet bound is refused, not projected", {
  # Dry: the collar is pinned a millionth of the interval's width off the wet
  # bound, while a 1e-3 trait step moves that bound by two orders more. So one
  # arm crosses, and the clamp on that arm silently answers about the bound
  # rather than about the collar it was given.
  crossed <- FALSE
  for (psi in c(4.2, 5.0)) {
    margin <- local({
      a <- fixed_collar_arm(psi, "root_b", 1)
      a$collar - a$bound_a
    })
    # The fixture must be AT a pin or it is testing nothing.
    expect_lt(margin, 1e-5)

    for (side in c(1, -1)) {
      a <- fixed_collar_arm(psi, "root_b", side)
      moved <- abs(a$moved_bound_a - a$bound_a)
      expect_gt(moved, margin)          # the bound moves further than the margin
      if (!a$feasible) {
        crossed <- TRUE
        # The clamp's tell: it reseats the collar somewhere else and returns a
        # number anyway.
        expect_gt(abs(a$clamped_at - a$collar), 0)
        expect_true(is.finite(a$clamped))
      }
    }
  }
  # Non-vacuity: at least one arm really did cross, or the test above asserted
  # nothing about the branch it exists for.
  expect_true(crossed)
})

test_that("the clamped route's difference is the feasibility jump", {
  # What the row looked like before it was refused. The two arms sit on opposite
  # sides of a discontinuity in profit, so their difference divided by the step
  # is enormous and perfectly finite -- this design's worst failure shape.
  base <- leaf_traits()
  h <- max(abs(base$root_b), 1) * 1e-3
  up <- fixed_collar_arm(5.0, "root_b", 1)
  dn <- fixed_collar_arm(5.0, "root_b", -1)
  jump <- (up$clamped - dn$clamped) / (2 * h)
  message(sprintf("  clamped dprofit/droot_b at a pin: %.6g (against %.6g well inside)",
                  jump,
                  local({
                    u <- fixed_collar_arm(2.0, "root_b", 1)
                    d <- fixed_collar_arm(2.0, "root_b", -1)
                    (u$clamped - d$clamped) / (2 * h)
                  })))
  expect_gt(abs(jump), 100)
  expect_true(is.finite(jump))
})
