# The aligned defaults (#156). Each curve's level parameter defaults to the value
# at which its stem optimum has TF24's conductance at one reference environment;
# these tests re-solve that alignment and hold the shipped defaults to it.

ref_gc <- function(model, ...) {
  leaf_solve(psi_soil = 1, model = model, supply = leaf_supply_singlelayer(),
             ...)$gc
}

test_that("align_curves() reproduces the shipped defaults", {
  # ⚠️ IF THIS FAILS AFTER A MODEL CHANGE, RE-ALIGN; DO NOT LOOSEN IT. A change
  # that moves TF24 or one of these curves at the reference moves the aligned
  # values, and the defaults in leaf.hpp, leaf_traits(),
  # .leaf_trait_defaults and .level_defaults must follow, ten significant
  # figures each.
  a <- align_curves()
  expect_named(a, c("CF77_lambda", "JS22_gamma", "CMax_a", "LeastCost_beta"))
  expect_equal(a[["CF77_lambda"]], .level_defaults[["CF77_lambda_"]],
               tolerance = 1e-8)
  expect_equal(a[["LeastCost_beta"]], .level_defaults[["LeastCost_beta"]],
               tolerance = 1e-8)
  expect_equal(a[["JS22_gamma"]], leaf_traits()$JS22_gamma, tolerance = 1e-8)
  expect_equal(a[["CMax_a"]], leaf_traits()$CMax_a, tolerance = 1e-8)
  expect_equal(attr(a, "target"), ref_gc("TF24"), tolerance = 0)
})

test_that("at their defaults the aligned curves have TF24's conductance", {
  tf24 <- ref_gc("TF24")
  for (m in c("CF77", "JS22", "CMax", "LeastCost")) {
    expect_equal(ref_gc(m), tf24, tolerance = 1e-8, info = m)
  }
  # And they differ away from the reference, which is the point of aligning:
  # a curve that matched everywhere would be the same curve.
  dry <- vapply(c("TF24", "JS22", "CMax"), function(m) {
    leaf_solve(psi_soil = 2.5, model = m, supply = leaf_supply_singlelayer())$gc
  }, numeric(1))
  expect_gt(diff(range(dry)) / mean(dry), 0.01)
})

test_that("CF77 aligned to TF24 is TF24's emergent lambda", {
  # At TF24's own marginal cost of water the two first-order conditions
  # coincide, so this is an identity rather than a fit.
  tf24 <- leaf_solve(psi_soil = 1, model = "TF24",
                     supply = leaf_supply_singlelayer())
  expect_equal(align_curves("CF77")[["CF77_lambda"]], tf24$lambda_emergent,
               tolerance = 1e-8)
})

test_that("the R defaults are the leaf's own", {
  # .level_defaults is the R copy of two in-class initialisers in
  # leaf.hpp; a fresh leaf reads the C++ ones.
  l <- leaf_model()
  expect_identical(l$CF77_lambda_, .level_defaults[["CF77_lambda_"]])
  expect_identical(l$LeastCost_beta, .level_defaults[["LeastCost_beta"]])
  expect_true(is.na(l$TF24_floor_lambda_o))
  # And leaf_traits()'s two are .leaf_trait_defaults' two.
  expect_identical(leaf_traits()$JS22_gamma, .leaf_trait_defaults$JS22_gamma)
  expect_identical(leaf_traits()$CMax_a, .leaf_trait_defaults$CMax_a)
})

test_that("a batch gradient keeps the batch leaf's aligned defaults", {
  # The C++ batch writes the owned slots from theta, which is NA on every
  # route that does not own them; writing that NA cleared the default, and a
  # later CF77 or LeastCost solve on the same leaf was refused.
  b <- leaf_batch(psi_soil = c(1, 2), PPFD = 900,
                  supply = leaf_supply_singlelayer())
  leaf_gradient_batch(b, model = "TF24", pars = "vcmax_25")
  expect_identical(b$leaf$CF77_lambda_, .level_defaults[["CF77_lambda_"]])
  expect_identical(b$leaf$LeastCost_beta, .level_defaults[["LeastCost_beta"]])
})

test_that("an NA owned slot in a supplied theta takes the default", {
  b <- leaf_batch(psi_soil = c(1, 2), PPFD = 900,
                  supply = leaf_supply_singlelayer())
  built <- leaf_gradient_batch(b, model = "LeastCost", pars = "vcmax_25")
  th <- .gradient_theta_matrix(b, b$traits, model = "LeastCost")
  th[, match("LeastCost_beta", .gradient_par_names())] <- NA_real_
  given <- leaf_gradient_batch(b, model = "LeastCost", pars = "vcmax_25",
                               theta = th)
  expect_true(all(given$status != "error"))
  expect_identical(given$gradient, built$gradient)
})

test_that("align_curves() follows its reference and its target", {
  w <- align_curves(c("JS22", "CMax"), psi_soil = 0.5)
  tf24 <- leaf_solve(psi_soil = 0.5, model = "TF24",
                     supply = leaf_supply_singlelayer())$gc
  expect_equal(attr(w, "target"), tf24, tolerance = 0)
  got <- leaf_solve(psi_soil = 0.5, model = "CMax",
                    supply = leaf_supply_singlelayer(),
                    traits = leaf_traits(CMax_a = w[["CMax_a"]]))$gc
  expect_equal(got, tf24, tolerance = 1e-8)

  t <- align_curves("LeastCost", target = 0.06)
  expect_equal(ref_gc("LeastCost", LeastCost_beta = t[["LeastCost_beta"]]),
               0.06, tolerance = 1e-8)
})

test_that("align_curves() refuses what it cannot align, and says why", {
  for (m in c("SOX", "JW26", "ProfitMax")) {
    expect_error(align_curves(m), "no level parameter|have a level parameter")
  }
  expect_error(align_curves("TF24_floor"), "lambda_o = 0")
  expect_error(align_curves("CF77", target = 10), "cannot reach")
})

test_that("an unset owned price takes the default, and only on its own route", {
  # leaf_gradient() with no leaf seats the default on the owning route ...
  g <- leaf_gradient(psi_soil = 1, PPFD = 900, model = "CF77",
                     supply = leaf_supply_singlelayer(),
                     pars = "CF77_lambda_")
  expect_true(all(is.finite(g$gradient[, c("A", "gc", "psi_stem")])))
  expect_lt(g$gradient[, "gc"], 0)  # a dearer price of water closes the leaf
  # ... and NA everywhere else, so the R reference writes nothing there.
  expect_true(is.na(.owned_value("TF24", "CF77_lambda_")))
  expect_true(is.na(.owned_value("collar", "LeastCost_beta")))
  expect_true(is.na(.owned_value("TF24_floor", "TF24_floor_lambda_o")))
})
