# The stem curve's stem_b row, and the homogeneity the four accessors rest on.
#
# G(psi; b) = s * Ghat(psi/s) with s = b / b_spline, so G is homogeneous of
# degree one in (psi, b). That buys two things: a moved stem_b needs no spline
# rebuild, and dG/db follows from Euler's theorem rather than from a difference.
#
# ⚠️ EVERY CHECK HERE IS RUN AT A NON-UNIT SCALE AS WELL AS AT ONE, AND THAT IS
# THE POINT. On a freshly built leaf stem_b == stem_b_spline_, so s is exactly
# 1.0 and dividing by it is the identity -- which means a misread scale is
# invisible on any fixture that does not perturb stem_b first. A suite that only
# ever built fresh leaves would pass with the scaling wrong.

stem_curve_traits <- function() leaf_traits()

# The same leaf reached two ways: rebuild the spline at the new stem_b, or scale
# the old one. They must agree, and that agreement IS the homogeneity claim.
rebuilt_leaf <- function(b_new) {
  tr <- stem_curve_traits()
  tr$stem_b <- b_new
  leaf_model(traits = tr)
}

scaled_leaf <- function(b_new) {
  l <- leaf_model(traits = stem_curve_traits())
  if (b_new != stem_curve_traits()$stem_b) l$perturb_stem_b(b_new)
  l
}

# The spline is laid out to the potential at which the vulnerability function
# reaches a fixed small fraction, so its domain moves with stem_b.
stem_curve_domain <- function(b) b * log(100)^(1 / stem_curve_traits()$stem_c)

sample_psi <- function(b, n = 15) {
  seq(0.1 * stem_curve_domain(b), 0.9 * stem_curve_domain(b), length.out = n)
}

test_that("a scaled stem curve reproduces a rebuilt one", {
  b0 <- stem_curve_traits()$stem_b
  for (f in c(0.75, 1.0, 1.4)) {
    b <- f * b0
    scaled <- scaled_leaf(b)
    rebuilt <- rebuilt_leaf(b)
    for (psi in sample_psi(b)) {
      # ⚠️ THE BAND IS THE GRID'S, AND IT MOVES WITH THE KNOT COUNT. The two
      # routes read splines laid out on different abscissae -- the scaled one at
      # psi/s, the rebuilt one at psi -- so what is left is each interpolant's
      # own resolution of the same curve, and a finer grid puts more of the
      # sample points near a span boundary rather than fewer. Measured over the
      # 45 sample points: 1 past 1e-12 at 400 knots, 3 at 800, 6 at 1600, with
      # the worst going 1.0e-12, 2.9e-12, 3.3e-12. The default is 1600, so the
      # bound is set above what that grid supports and stays four orders below
      # the interpolant's own value error.
      expect_equal(scaled$stem_curve_integral(psi),
                   rebuilt$stem_curve_integral(psi), tolerance = 1e-11)
      expect_equal(scaled$stem_curve_integral_deriv(psi),
                   rebuilt$stem_curve_integral_deriv(psi), tolerance = 1e-11)
    }
  }
})

test_that("Euler's identity holds on the stem curve", {
  # psi*G'(psi) + b*dG/db == G(psi), by homogeneity of degree one.
  #
  # What this catches is an expression-level slip in dG/db -- a dropped 1/b moves
  # it by ORDERS. What it cannot catch is an error in G' itself, because dG/db is
  # built from G' and the error cancels; the next test is the one for that, and
  # the two are not interchangeable.
  b0 <- stem_curve_traits()$stem_b
  worst <- 0
  for (f in c(0.75, 1.0, 1.4)) {
    b <- f * b0
    l <- scaled_leaf(b)
    for (psi in sample_psi(b)) {
      G <- l$stem_curve_integral(psi)
      lhs <- psi * l$stem_curve_integral_deriv(psi) +
        b * l$stem_curve_integral_dstem_b(psi)
      worst <- max(worst, abs(lhs - G) / abs(G))
    }
  }
  message(sprintf("  Euler identity: worst relative %.2e", worst))
  expect_lt(worst, 1e-13)
})

test_that("the Euler row agrees with a rebuilt difference", {
  # The referee that shares no code with the row: rebuild the spline either side
  # of stem_b and difference G there. This is what sees a misread scale, and it
  # sees it only at a scale that is not one.
  b0 <- stem_curve_traits()$stem_b
  worst <- 0
  for (f in c(0.75, 1.0, 1.4)) {
    b <- f * b0
    l <- scaled_leaf(b)
    h <- 1e-5 * b0
    for (psi in sample_psi(b, n = 7)) {
      fd <- (rebuilt_leaf(b + h)$stem_curve_integral(psi) -
               rebuilt_leaf(b - h)$stem_curve_integral(psi)) / (2 * h)
      worst <- max(worst, abs(l$stem_curve_integral_dstem_b(psi) / fd - 1))
    }
  }
  message(sprintf("  Euler vs rebuilt difference: worst relative %.2e", worst))
  # The floor here is the central difference's, not the row's.
  expect_lt(worst, 1e-6)
})
