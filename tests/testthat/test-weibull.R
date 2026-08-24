# The (P50, c) converter.
#
# The load-bearing tests here are the ROUND TRIPS. A converter's failure mode is
# not that it errors, it is that it returns a plausible pair that describes a
# different curve -- so every path is checked by converting away and back, against
# the model's own derived quantities rather than against a second copy of the
# formula.

test_that("psi_at_plc agrees with the model's own derived quantiles", {
  # The model derives stem_b and psi_crit from (P50, c) internally. If this
  # function disagrees with those, one of the two is wrong -- and the model's is
  # the one the solve uses, so it is the reference.
  tr <- leaf_traits()
  l <- leaf_model(tr)

  # psi_crit IS P95 by construction: the critical fraction is 0.05 remaining.
  expect_equal(psi_at_plc(tr$stem_P50, tr$stem_c, 0.95), l$psi_crit)
  # And P50 is its own 50% point, exactly.
  expect_equal(psi_at_plc(tr$stem_P50, tr$stem_c, 0.5), tr$stem_P50)
  # The root curve too, which has its own derived critical potential.
  expect_equal(psi_at_plc(tr$root_P50, tr$root_c, 0.95), l$root_psi_crit)

  # Monotone in plc, which is what makes the inversions below single-valued.
  p <- vapply(c(0.1, 0.3, 0.5, 0.88, 0.95, 0.99),
              function(q) psi_at_plc(3.4, 2.68, q), numeric(1))
  expect_true(all(diff(p) > 0))
})

test_that("weibull_s50 matches a finite difference of the curve", {
  # S50 is quoted in % PLC per MPa, so the comparison needs the same units: the
  # derivative of (1 - f) at P50, times 100.
  P50 <- 3.4
  cc <- 2.680147
  f <- function(psi) 2^-((psi / P50)^cc)
  h <- 1e-6
  fd <- 100 * -(f(P50 + h) - f(P50 - h)) / (2 * h)
  expect_equal(weibull_s50(P50, cc), fd, tolerance = 1e-6)
})

test_that("every pair of inputs recovers the curve it came from", {
  P50 <- 3.4
  cc <- 2.680147
  P88 <- psi_at_plc(P50, cc, 0.88)
  P95 <- psi_at_plc(P50, cc, 0.95)
  P99 <- psi_at_plc(P50, cc, 0.99)
  S50 <- weibull_s50(P50, cc)

  same <- function(got) {
    expect_equal(got$P50, P50, tolerance = 1e-8)
    expect_equal(got$c, cc, tolerance = 1e-8)
  }

  same(weibull_p50_c(P50 = P50, c = cc))          # identity
  same(weibull_p50_c(P50 = P50, P88 = P88))       # two quantiles, one is P50
  same(weibull_p50_c(P88 = P88, P95 = P95))       # two quantiles, neither is P50
  same(weibull_p50_c(P95 = P95, P99 = P99))
  same(weibull_p50_c(P50 = P50, S50 = S50))       # the collapsing S50 case
  same(weibull_p50_c(S50 = S50, c = cc))
  same(weibull_p50_c(P95 = P95, c = cc))          # a quantile plus a shape
  same(weibull_p50_c(px = P95, plc = 0.95, c = cc))
  # ⚠️ The one transcendental case, solved by root-find rather than in closed
  # form. It is the reason the root-find exists, so it is asserted rather than
  # assumed to be reachable.
  same(weibull_p50_c(P95 = P95, S50 = S50))
  same(weibull_p50_c(px = P88, plc = 0.88, S50 = S50))
})

test_that("S50 with a non-P50 quantile has two solutions, and says so", {
  # ⚠️ NOT a defensive message: S50(c) = K*c*L^(1/c)/psi_f has an analytic minimum
  # at c = ln(L), so it is U-shaped and any S50 above that minimum is matched by
  # two curves. This is what the round-trip test found by failing with "no root in
  # (0, 20]" -- both bracket ends were positive because the root sat in a dip.
  P50 <- 3.4
  cc <- 2.680147
  P95 <- psi_at_plc(P50, cc, 0.95)
  S50 <- weibull_s50(P50, cc)

  expect_message(got <- weibull_p50_c(P95 = P95, S50 = S50), "two solutions")
  # The upper branch is returned, and it is the one that came in.
  expect_equal(got$c, cc, tolerance = 1e-8)

  # The other root is real, not an artefact of the search: it reproduces the same
  # S50 at the same P95 to full precision.
  L <- log(1 / 0.05) / log(2)
  other_c <- 0.8846
  other_P50 <- P95 * L^(-1 / other_c)
  expect_equal(weibull_s50(other_P50, other_c), S50, tolerance = 1e-3)
  # And it lands below 1, which is why the upper branch is the right default:
  # measured angiosperm stems sit at 1.8-2.6.
  expect_lt(other_c, 1)

  # When the quantile IS P50 the ambiguity disappears -- L = 1 collapses the
  # substitution -- so there is no message and no root-find.
  expect_silent(weibull_p50_c(P50 = P50, S50 = S50))

  # A slope flatter than any curve can manage at that quantile is refused, and the
  # message names the achievable minimum rather than just failing.
  expect_error(weibull_p50_c(P95 = P95, S50 = 1), "smallest slope")
})

test_that("the arbitrary quantile is the same thing as the named ones", {
  P50 <- 2.9
  cc <- 1.8
  for (q in c(0.5, 0.88, 0.95, 0.99)) {
    psi <- psi_at_plc(P50, cc, q)
    a <- weibull_p50_c(px = psi, plc = q, c = cc)
    expect_equal(a$P50, P50, tolerance = 1e-8)
  }
})

test_that("it refuses what does not determine a curve, and says why", {
  expect_error(weibull_p50_c(P50 = 3.4), "exactly two")
  expect_error(weibull_p50_c(), "exactly two")
  expect_error(weibull_p50_c(P50 = 3.4, P88 = 5, c = 2), "exactly two")
  # A quantile is one quantity, so its two arguments must arrive together.
  expect_error(weibull_p50_c(px = 5.0, c = 2), "given together")
  expect_error(weibull_p50_c(plc = 0.9, c = 2), "given together")
  # Two potentials at the same loss carry one piece of information.
  expect_error(weibull_p50_c(P95 = 5.87, px = 5.87, plc = 0.95),
               "one piece of information")
})

test_that("the sign convention is enforced rather than silently repaired", {
  # ⚠️ Published P50 is NEGATIVE. Negating it here would let a caller working in
  # the other convention pass everything signed and get plausible answers, so it
  # is refused instead.
  expect_error(weibull_p50_c(P50 = -3.4, c = 2), "positive magnitude")
  expect_error(psi_at_plc(-3.4, 2, 0.5), "positive magnitude")
  expect_error(weibull_s50(-3.4, 2), "positive magnitude")
  expect_error(psi_at_plc(3.4, -2, 0.5), "positive")
  # plc is a proportion, not a percentage.
  expect_error(psi_at_plc(3.4, 2, 95), "\\(0, 1\\)")
})

test_that("a shape below 1 warns without refusing", {
  # Possible, but usually a native-embolism artefact rather than a trait, so the
  # caller is told and handed the answer.
  expect_warning(got <- weibull_p50_c(P50 = 3.4, c = 0.6), "below 1")
  expect_equal(got$c, 0.6)
})

test_that("the result splices into leaf_traits() and reaches the model", {
  # The point of the converter: a published pair goes in, a solved leaf comes out.
  pars <- weibull_p50_c(P50 = 2.9, P88 = 4.6)
  tr <- leaf_traits(stem_P50 = pars$P50, stem_c = pars$c)
  l <- leaf_model(tr)
  expect_equal(l$stem_P50, pars$P50)
  # And the derived quantities follow from it, which is the property that made
  # psi_crit stop being a trait.
  expect_equal(l$psi_crit, psi_at_plc(pars$P50, pars$c, 0.95))
  expect_equal(l$stem_b, pars$P50 / log(2)^(1 / pars$c))
})

test_that("b <-> P50 round-trips against the model's own derived stem_b", {
  tr <- leaf_traits()
  l <- leaf_model(tr)
  # The model derives stem_b from stem_P50; this inverts it. Asserted against the
  # object rather than against the formula, so the two cannot drift.
  expect_equal(weibull_p50_from_b(l$stem_b, tr$stem_c), tr$stem_P50)
  expect_equal(weibull_p50_from_b(l$root_b, tr$root_c), tr$root_P50)

  # ⚠️ P50 < b always, so a swapped argument is plausible rather than an error.
  # Pinned as a number so the size of that mistake is on record.
  expect_lt(weibull_p50_from_b(3.898245, 2.680147), 3.898245)
  expect_equal(weibull_p50_from_b(3.898245, 2.680147) / 3.898245,
               log(2)^(1 / 2.680147))
  expect_error(weibull_p50_from_b(-3.9, 2.68), "positive magnitude")
})
