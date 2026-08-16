# The ROOT vulnerability curve's cumulative integral and its row in root_b.
#
# ⚠️ HAZARD 1. This is not the stem curve. It is built from root_b and root_c and
# it drives uptake; the stem's is built from stem_b and stem_c and drives the
# hydraulic cost. Reading one for the other is the mistake that carried a wrong
# exponent into a manuscript draft.
#
# root_b scales the curve rather than reshaping it, so G is homogeneous of degree
# one in (psi, root_b) and Euler gives dG/droot_b with no rebuild -- the same
# identity the stem curve uses, on the same builder. root_c has no counterpart.
#
# The row is reached by the bound rows already (test-bound-row.R differences
# find_root_psi in root_b), but only THROUGH the uptake integral and the collar
# solve. This checks the quantity itself, which is the only way to say whether a
# disagreement there belongs to the curve or to what consumes it.

root_curve_leaf <- function(tr = leaf_traits()) leaf_model(traits = tr)

# Past this the integral is the complete-gamma limit and dG/dpsi is zero.
root_curve_capped_psi <- 9.0

test_that("the homogeneity identity holds across the curve, cap included", {
  tr <- leaf_traits()
  l <- root_curve_leaf(tr)
  worst <- 0
  for (psi in c(0.5, 1, 2, 3, 5, 6.5, 7, 8, root_curve_capped_psi)) {
    G <- l$root_curve_integral(psi)
    resid <- psi * l$root_curve_integral_deriv(psi) +
      tr$root_b * l$root_curve_integral_droot_b(psi) - G
    worst <- max(worst, abs(resid) / G)
  }
  message(sprintf("  Euler identity: worst relative %.2e", worst))
  expect_lt(worst, 1e-15)
})

test_that("the identity is blind to the derivative it is assembled from", {
  # ⚠️ WHAT THIS INVARIANT CANNOT SEE, stated as arithmetic rather than as a
  # caveat. dG/droot_b is DEFINED as (G - psi*G')/root_b, so
  #
  #   psi*G' + root_b*((G - psi*G')/root_b)  ==  G
  #
  # for ANY G' whatsoever. Substituting a corrupted derivative below leaves the
  # residual at exactly zero, so a check built on the identity alone would pass
  # while the row it certifies was wrong by whatever factor the corruption was.
  # The rebuilt difference in the next block is what is not blind to it.
  tr <- leaf_traits()
  l <- root_curve_leaf(tr)
  psi <- 3.0
  G <- l$root_curve_integral(psi)
  for (corrupt in c(1.0, 0.5, 2.0, -1.0)) {
    Gp <- corrupt * l$root_curve_integral_deriv(psi)
    Gb <- (G - psi * Gp) / tr$root_b          # the row, rebuilt from that G'
    expect_equal(psi * Gp + tr$root_b * Gb - G, 0)
  }
})

test_that("the row agrees with a rebuilt difference at every step tried", {
  # The referee that is not blind: rebuild the curve at a perturbed root_b and
  # difference the value. Shares no code with the row -- one evaluates the
  # integral, the other differentiates it.
  tr <- leaf_traits()
  l <- root_curve_leaf(tr)
  worst <- 0
  flat <- TRUE
  for (psi in c(1, 3, 5, root_curve_capped_psi)) {
    row <- l$root_curve_integral_droot_b(psi)
    ratios <- numeric(0)
    # FIVE steps over four orders. The composed row in test-bound-row.R needs
    # 1e-5 or coarser; this one does not, and the contrast is the finding: the
    # knot-grid artefact that puts 1e-6 off the plateau there is a property of
    # the collar solve the row passes through, not of the curve.
    for (rel in c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3)) {
      h <- tr$root_b * rel
      up <- tr; up$root_b <- up$root_b + h
      dn <- tr; dn$root_b <- dn$root_b - h
      fd <- (root_curve_leaf(up)$root_curve_integral(psi) -
             root_curve_leaf(dn)$root_curve_integral(psi)) / (2 * h)
      ratios <- c(ratios, row / fd)
      worst <- max(worst, abs(row / fd - 1))
    }
    if (diff(range(ratios)) > 1e-5) flat <- FALSE
  }
  message(sprintf("  row vs rebuilt difference: worst %.2e over 4 potentials x 5 steps",
                  worst))
  expect_lt(worst, 1e-4)
  # The plateau covers every step tried, which is what says the step was not
  # chosen to make this pass.
  expect_true(flat)
})

test_that("the capped branch is reached, and is the limit's own derivative", {
  # Past the ceiling G is (root_b/root_c)*Gamma(1/root_c) -- LINEAR in root_b and
  # flat in psi -- so the identity returns G/root_b, which is that limit's exact
  # derivative. A branch, and one no bound row reaches: the collar potentials the
  # bounds sit at are well inside the curve, so without this the capped arm is
  # code nothing runs.
  tr <- leaf_traits()
  l <- root_curve_leaf(tr)

  # The fixture must be shown to be past the cap rather than assumed to be.
  expect_equal(l$root_curve_integral_deriv(root_curve_capped_psi), 0)
  expect_gt(l$root_curve_integral_deriv(5.0), 0)
  expect_equal(l$root_curve_integral(root_curve_capped_psi),
               l$root_curve_integral(root_curve_capped_psi + 5))

  G <- l$root_curve_integral(root_curve_capped_psi)
  expect_equal(l$root_curve_integral_droot_b(root_curve_capped_psi), G / tr$root_b)

  # And the same at a different root_b -- but the potential has to scale with it.
  # ⚠️ THE CAP MOVES. It binds where the curve has run out, and root_b is what
  # sets that width, so a fixture that holds psi fixed while scaling root_b walks
  # back onto the rising part and reads 1.298 for a ratio that is 1.3 exactly.
  # That is the homogeneity again rather than a second fact about the cap.
  s <- 1.3
  tr2 <- tr; tr2$root_b <- tr$root_b * s
  l2 <- root_curve_leaf(tr2)
  psi2 <- root_curve_capped_psi * s
  expect_equal(l2$root_curve_integral_deriv(psi2), 0)
  expect_equal(l2$root_curve_integral(psi2) / tr2$root_b,
               l2$root_curve_integral_droot_b(psi2))
  expect_equal(l2$root_curve_integral(psi2) / G, s, tolerance = 1e-9)
})

test_that("the root curve has no unit-scale blind fixture, unlike the stem's", {
  # The stem curve is a reference spline RESCALED by stem_b/stem_b_spline_, so on
  # a freshly built leaf that scale is exactly 1.0 and a misread of it is
  # invisible -- which is why test-stem-curve.R has to perturb stem_b first.
  #
  # The root curve is rebuilt and cached against the (root_b, root_c) pair that
  # determines it, so there is no scale to misread and no such fixture
  # requirement. Asserted rather than assumed, because it is the reason this file
  # does not carry the stem's warning.
  tr <- leaf_traits()
  a <- root_curve_leaf(tr)$root_curve_integral(3.0)
  tr2 <- tr; tr2$root_b <- tr$root_b * 1.25
  b <- root_curve_leaf(tr2)$root_curve_integral(3.0)
  expect_gt(abs(b / a - 1), 0.01)
  # A curve that is genuinely rebuilt reproduces a scaled one exactly: G is
  # homogeneous, so G(psi; s*b) == s*G(psi/s; b).
  s <- 1.25
  expect_equal(root_curve_leaf(tr2)$root_curve_integral(3.0 * s),
               s * root_curve_leaf(tr)$root_curve_integral(3.0),
               tolerance = 1e-9)
})
