# The profit's rows in the environment -- radiation and each soil layer -- at
# the operating point the leaf is seated at.
#
# The referee is a central difference of the profit at a RE-SOLVED operating
# point. That is the relation the rows are the derivative of, and it shares no
# code with the derivation: one solves for the collar and reads profit, the other
# never solves at all.
#
# ⚠️ THE PRICE CARRIES THE STATIONARITY CONDITION, which is why these rows are
# not one derivation. `marginal_price_water` is lambda*kmax*f(p)/S, and that IS
# dProfit/dE_up only once dProfit/dp is zero -- the interior condition is inside
# the expression rather than beside it. At a pin dProfit/dp is `nu`, and the
# frozen-collar price is a DIFFERENT object -- dmarginal_profit_duptake_slope,
# which is built from the cost and assimilation kernels with no stationarity in
# it anywhere. Measured, using the interior price at a pin leaves the row 8% out
# at a wet pin and 15% to 47% out at a dry one, with every intermediate finite.
#
# The interior branch keeps the interior price. The two agree there by the
# first-order condition and differ in the last bits, and a gradient that already
# answers must not move.

env_row_seat <- function(psi_soil, ...) {
  l <- leaf_model(traits = leaf_traits())
  set_drivers(l, psi_soil = psi_soil, ...)
  l$find_root_collar_psi()
  l
}

env_row_profile <- c(2.6, 2.8, 3.0, 3.2, 3.4)

# The three regimes, each named by the branch its solve takes rather than by how
# close the collar sits to a bound.
env_row_states <- list(
  interior = list(psi = env_row_profile, args = list(), kind = "interior"),
  wet      = list(psi = env_row_profile * 1.4, args = list(), kind = "pinned-wet"),
  dry      = list(psi = env_row_profile,
                  args = list(leaf_specific_conductance_max = 3e-4),
                  kind = "pinned-dry-root-crit"),
  dry_hard = list(psi = env_row_profile,
                  args = list(leaf_specific_conductance_max = 1e-3),
                  kind = "pinned-dry-root-crit")
)

env_row_leaf <- function(st, psi = st$psi, extra = list()) {
  do.call(env_row_seat, c(list(psi_soil = psi), st$args, extra))
}

# dProfit*/dpsi_j by re-solving either side. Returns NA where an arm left the
# branch, because a difference across a change of branch is a difference of two
# different functions.
env_row_total <- function(st, j, h) {
  up <- st$psi; up[[j]] <- up[[j]] + h
  dn <- st$psi; dn[[j]] <- dn[[j]] - h
  a <- env_row_leaf(st, up)
  b <- env_row_leaf(st, dn)
  if (a$operating_point_kind_name() != st$kind ||
      b$operating_point_kind_name() != st$kind) {
    return(NA_real_)
  }
  (a$profit_ - b$profit_) / (2 * h)
}

test_that("the environment rows agree with a re-solved difference in every regime", {
  # Tolerances are per regime because the regimes do not have one floor. Interior
  # and the wet pin sit on the model's own consistency floor; the dry pin sits an
  # order above it, and §"the dry pin's residual" below says what was ruled out.
  tol <- c(interior = 1e-4, wet = 1e-6, dry = 1e-3, dry_hard = 1e-3)
  for (nm in names(env_row_states)) {
    st <- env_row_states[[nm]]
    l <- env_row_leaf(st)
    expect_identical(l$operating_point_kind_name(), st$kind)
    r <- profit_env_row_values(l)
    expect_true(r$usable)
    expect_identical(r$pinned, st$kind != "interior")

    # Four steps over three orders, so the plateau can be SEEN rather than
    # assumed -- a single reading cannot tell a converged difference from a
    # step-size artefact. 1e-6 is the loose end and is reported, not asserted:
    # the re-solve's own floor is what widens it, and a check held to the
    # coarsest reading it can get would be a check on the referee.
    worst <- setNames(numeric(4), c("1e-6", "1e-5", "1e-4", "1e-3"))
    for (j in seq_along(st$psi)) {
      for (k in seq_along(worst)) {
        h <- c(1e-6, 1e-5, 1e-4, 1e-3)[[k]]
        fd <- env_row_total(st, j, h)
        expect_false(is.na(fd))
        worst[[k]] <- max(worst[[k]], abs(r$dprofit_dpsi_soil[[j]] / fd - 1))
      }
    }
    message(sprintf("  %-9s (%s): worst by step %s",
                    nm, st$kind, paste(sprintf("%.1e", worst), collapse = " ")))
    expect_lt(max(worst[-1]), tol[[nm]])
  }
})

test_that("a wet pin cannot referee the price correction, and a dry pin can", {
  # ⚠️ THE BLIND SPOT, asserted rather than described. At the wet bound the
  # residual is total uptake alone, so dB/dpsi_j is exactly -(dE_up/dpsi_j)/S --
  # and the two prices differ by exactly nu/S. So the bound term cancels the whole
  # difference between them, the interior formula lands on the right answer, and
  # this fixture passes whichever price is used.
  #
  # The dry bound's slope carries the stem's half as well, so nothing cancels
  # there. That is the only fixture the derivation is evidence from.
  st <- env_row_states$wet
  l <- env_row_leaf(st)
  r <- profit_env_row_values(l)
  nu <- l$dprofit_droot_collar_psi(l$opt_root_psi_)
  wet <- l$bound_row_values(0L)
  expect_true(abs(nu) > 1e-3)          # the correction is not zero for want of nu

  # The operational form of the blind spot: the wet fixture agrees with its
  # reference far better than the dry one does, whatever the derivation, so it is
  # not the fixture that discriminates between them. Two orders, measured -- if
  # this ever stops holding, the cancellation has changed and the dry fixture is
  # no longer carrying the evidence on its own.
  resid <- function(nm) {
    s <- env_row_states[[nm]]
    rr <- profit_env_row_values(env_row_leaf(s))
    max(vapply(seq_along(s$psi),
               function(j) abs(rr$dprofit_dpsi_soil[[j]] /
                                 env_row_total(s, j, 1e-4) - 1), numeric(1)))
  }
  wet_r <- resid("wet"); dry_r <- resid("dry")
  message(sprintf("  wet residual %.1e against dry %.1e", wet_r, dry_r))
  expect_lt(wet_r * 100, dry_r)

  # And the half that WAS shipped on its own -- interior price, plus nu*dB -- is
  # wrong here by a margin nothing about the wet bound makes small.
  half <- vapply(seq_along(st$psi),
                 function(j) r$dprofit_dpsi_soil[[j]] + nu * wet[[8 + j]],
                 numeric(1))
  off <- max(abs(half / unlist(r$dprofit_dpsi_soil) - 1))
  message(sprintf("  wet pin: the bound term alone would move the row by %.1f%%",
                  100 * off))
  expect_gt(off, 0.05)
})

test_that("the dry pin is where both halves of the correction are load-bearing", {
  # Non-vacuity for the correction itself. Reconstruct the two wrong rows the
  # correction sits between and require each to miss the re-solved total by a
  # margin far outside the tolerance the correct row is held to.
  for (nm in c("dry", "dry_hard")) {
    st <- env_row_states[[nm]]
    l <- env_row_leaf(st)
    r <- profit_env_row_values(l)
    nu <- l$dprofit_droot_collar_psi(l$opt_root_psi_)
    dry <- l$bound_row_values(1L)
    expect_equal(dry[[1]], 1)

    j <- 1L
    fd <- env_row_total(st, j, 1e-4)
    row <- r$dprofit_dpsi_soil[[j]]
    # (a) the bound term dropped -- the frozen-collar partial alone.
    without_bound <- row - nu * dry[[8 + j]]
    # (b) the price correction dropped -- what was shipped before this check.
    #     nu/S is recovered from the two rows rather than restated: the bound
    #     row's own slope is S + kmax*G'(x) on this arm.
    S <- dry[[3]] - l$leaf_specific_conductance_max_ *
      l$stem_curve_integral_deriv(dry[[2]])
    duptake <- -dry[[8 + j]] * dry[[3]]
    without_price <- row - duptake * (nu / S)

    message(sprintf("  %-8s: correct %.5f  no bound term %.5f (%.0f%% out)  no price term %.5f (%.0f%% out)",
                    nm, row / fd, without_bound / fd,
                    100 * abs(without_bound / fd - 1), without_price / fd,
                    100 * abs(without_price / fd - 1)))
    expect_gt(abs(without_bound / fd - 1), 0.02)
    expect_gt(abs(without_price / fd - 1), 0.05)
  }
})

test_that("the radiation row is the same derivation at a pin as at an interior point", {
  # Neither bound reads radiation, so dB/dlight is exactly zero and the row picks
  # up no bound term. That is asserted as an EXACT equality of the bound itself
  # rather than as a small derivative, because "exactly zero" is the claim.
  for (nm in names(env_row_states)) {
    st <- env_row_states[[nm]]
    which <- if (st$kind == "pinned-dry-root-crit") 1L else 0L
    l <- env_row_leaf(st)
    b0 <- l$find_root_psi(min(st$psi), st$psi, which)
    bright <- env_row_leaf(st, extra = list(PPFD = 1200))
    expect_identical(b0, bright$find_root_psi(min(st$psi), st$psi, which))
  }

  # And the row itself, where it is big enough to referee. At the wet pin the
  # radiation row is ~5e-13 -- the leaf is barely transpiring there -- so a ratio
  # against a difference is dominated by the difference's own floor, and this
  # says so instead of asserting a tolerance it cannot support.
  for (nm in c("interior", "dry")) {
    st <- env_row_states[[nm]]
    r <- profit_env_row_values(env_row_leaf(st))
    # Around set_drivers' own default, or the row and the difference describe
    # two different radiation levels and the ratio is a reading of that gap.
    h <- 1
    a <- env_row_leaf(st, extra = list(PPFD = 900 + h))
    b <- env_row_leaf(st, extra = list(PPFD = 900 - h))
    expect_identical(a$operating_point_kind_name(), st$kind)
    fd <- (a$profit_ - b$profit_) / (2 * h)
    message(sprintf("  radiation row, %-8s: %.6g against %.6g", nm,
                    r$dprofit_dlight, fd))
    expect_equal(r$dprofit_dlight / fd, 1, tolerance = 1e-3)
  }
  expect_lt(abs(profit_env_row_values(env_row_leaf(env_row_states$wet))$dprofit_dlight),
            1e-9)
})

test_that("a point on no branch these rows serve is refused by name", {
  # The rows exist for an interior optimum and for a bound to follow. Anything
  # else has to come back unusable WITH a message, because a caller reading
  # sentinels as numbers is the failure this flag exists to prevent.
  l <- env_row_seat(rep(3.0, 5), PPFD = 0)
  r <- profit_env_row_values(l)
  expect_false(r$usable)
  expect_false(r$pinned)
  expect_gt(nchar(r$message), 0)
  expect_true(all(is.na(unlist(r$dprofit_dpsi_soil))))
})
