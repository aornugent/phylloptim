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

# The arms a row is differenced along when one side is not available.
#
# A collar held at the operating point can fall outside the interval the
# PERTURBED state has -- the bound moves out from under it -- and the stand meets
# this at points classified interior, sitting a hundred-thousandth of an MPa
# inside the dry bound. The step cannot always be reduced enough, so the row is
# taken one-sided on the side that is still inside.
#
# It is taken SECOND-ORDER, at the cost of one extra evaluation, and these
# numbers are why.

fixed_collar_interior <- function() {
  l <- leaf_model(traits = leaf_traits())
  set_drivers(l, psi_soil = c(2.6, 2.8, 3.0, 3.2, 3.4))
  l$find_root_collar_psi()
  l
}

# profit at a HELD collar, with one trait moved. The quantity every driven row is
# a difference of.
fixed_collar_profit <- function(collar, name, delta) {
  tr <- leaf_traits(); tr[[name]] <- tr[[name]] + delta
  m <- leaf_model(traits = tr)
  set_drivers(m, psi_soil = c(2.6, 2.8, 3.0, 3.2, 3.4))
  v <- m$profit_at_fixed_collar_values(collar)
  if (v[[1]] != 1) NA_real_ else v[[2]]
}

test_that("a one-sided second-order arm reproduces the centred one, from either side", {
  # The fixture is a state where BOTH sides are available, because that is the
  # only place the two can be compared at all -- where the one-sided form is
  # actually needed there is nothing to compare it against.
  l <- fixed_collar_interior()
  expect_identical(l$operating_point_kind_name(), "interior")
  collar <- l$opt_root_psi_
  fit_step <- 1e-3

  worst_2nd <- 0
  worst_1st <- 0
  for (nm in c("stem_c", "stem_b", "root_c", "root_b")) {
    h <- max(abs(leaf_traits()[[nm]]), 1) * fit_step
    f0 <- fixed_collar_profit(collar, nm, 0)
    at <- function(d) fixed_collar_profit(collar, nm, d)
    centred <- (at(h) - at(-h)) / (2 * h)
    forward <- (-3 * f0 + 4 * at(h) - at(2 * h)) / (2 * h)
    backward <- (3 * f0 - 4 * at(-h) + at(-2 * h)) / (2 * h)

    worst_2nd <- max(worst_2nd, abs(forward / centred - 1), abs(backward / centred - 1))
    worst_1st <- max(worst_1st, abs((at(h) - f0) / h / centred - 1),
                     abs((f0 - at(-h)) / h / centred - 1))

    # The property that matters most, and the one a first-order form does not
    # have: the two sides agree with EACH OTHER. Where this arm is really used,
    # which side is available is decided by which bound the point drifted toward
    # -- so a formula whose answer depends on the side would make the row a
    # function of the geometry rather than of the trait.
    expect_equal(forward / backward, 1, tolerance = 1e-5)
  }
  message(sprintf("  one-sided vs centred: second order %.2e, first order %.2e",
                  worst_2nd, worst_1st))
  expect_lt(worst_2nd, 1e-4)

  # Non-vacuity, and the reason for the extra evaluation. A first-order arm is
  # two orders worse and its error has a SIGN that follows the side, which is
  # what identifies it as truncation rather than noise.
  expect_gt(worst_1st, 20 * worst_2nd)
})

test_that("the side a one-sided arm must take is the side away from the bound", {
  # Which side is available is not a free choice, and the fixture states the rule
  # the stand relies on: a step that moves the bound AWAY from the held collar
  # stays feasible however far it goes, so the second point at 2h is available
  # whenever the first at h is.
  psi <- rep(5.0, 5)
  l <- leaf_model(traits = leaf_traits())
  set_drivers(l, psi_soil = psi)
  l$find_root_collar_psi()
  collar <- l$opt_root_psi_
  base <- leaf_traits()
  h <- max(abs(base$root_b), 1) * 1e-3

  ok <- function(d) {
    tr <- base; tr$root_b <- tr$root_b + d
    m <- leaf_model(traits = tr); set_drivers(m, psi_soil = psi)
    m$profit_at_fixed_collar_values(collar)[[1]] == 1
  }
  # One side crosses at this pin -- that is what test-fixed-collar's earlier
  # block measured -- and the other carries both of its points.
  expect_false(ok(h) && ok(-h))
  side <- if (ok(h)) 1 else -1
  expect_true(ok(side * h))
  expect_true(ok(side * 2 * h))
})
