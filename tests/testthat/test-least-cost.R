# The R surface of the least-cost curve (Prentice et al. 2014).
#
# WHAT IS AND IS NOT TESTED HERE. The physics is checked in C++, in
# `tests/cpp/test_leaf.cpp`, where the solved chi is compared against Prentice et
# al.'s first-order condition to 1e-11 and the published closed form's own error
# is predicted. This file covers what that cannot reach: that the curve is
# REACHABLE from R, that `LeastCost_beta` survives the trip through
# `leaf_solve()`, and that the refusals arrive as R conditions rather than as a
# crash. One closed-form check is repeated here anyway, because a curve that is
# wired to the wrong field would still solve and would still return a number.

test_that("LeastCost is in the registry and solves through leaf_solve()", {
  expect_true("LeastCost" %in% cost_curve_names())

  sp <- leaf_supply_singlelayer()
  r <- leaf_solve(psi_soil = 0.1, PPFD = 3000, atm_vpd = 2, ca = 40,
                  leaf_temp = 25, supply = sp, model = "LeastCost",
                  LeastCost_beta = 146)
  expect_true(is.finite(r$A))
  expect_true(is.finite(r$ci))
  expect_true(is.finite(r$lambda))
  # A ratio objective, so `profit` is not a carbon flux -- see
  # `profit_psi_stem_LeastCost`. Finiteness is all that is asserted of it.
  expect_true(is.finite(r$profit))
  # chi between the compensation point and one, which is the only statement
  # about the answer that needs no closed form.
  expect_gt(r$ci / 40, 0.1)
  expect_lt(r$ci / 40, 1.0)
})

test_that("the solved chi is Prentice et al.'s first-order condition", {
  # ⚠️ AGAINST THE FIRST-ORDER CONDITION, NOT THE PUBLISHED `chi = xi/(xi +
  # sqrt(D))`. Their reported form drops `Gamma*` from the numerator of the
  # condition, which understates chi by `(Gamma*/ca)/(1 + r)` -- some 4% at
  # D = 1 kPa and 8% at D = 4 kPa. The C++ test measures that gap; here the
  # exact condition is used, so the tolerance can be tight enough to catch a
  # mis-wiring.
  #
  # ⚠️ AND `xi` CARRIES 1.67, NOT THE PUBLISHED 1.6. The constant in `xi` is the
  # H2O:CO2 stomatal diffusion ratio, and this package's is 1.67 -- a convention,
  # not a property of the leaf. Using 1.6 against a solver that diffuses at 1.67
  # leaves a 0.6% error that looks like a defect and is not.
  sp <- leaf_supply_singlelayer()
  beta <- 146
  ratio <- 1.67
  ca <- 40

  for (vpd in c(1, 2, 3)) {
    r <- leaf_solve(psi_soil = 0.1, PPFD = 3000, atm_vpd = vpd, ca = ca,
                    leaf_temp = 25, supply = sp, model = "LeastCost",
                    # `kmax` raised so the optimum stays interior: least-cost has
                    # no vulnerability curve, so a solve pinned on `psi_crit` is
                    # answering a different question.
                    leaf_specific_conductance_max = 6.3e-3,
                    LeastCost_beta = beta)
    # Km and Gamma* are read off a leaf carrying the SAME drivers rather than
    # restated here, so a change to the Arrhenius curves moves the prediction
    # with the model instead of breaking the test.
    # The multi-layer supply here, not `sp`: Km and Gamma* are functions of
    # temperature, O2 and pressure alone, so the hydraulic path is irrelevant to
    # them, and the single-potential path refuses a carbon-built root network.
    l <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
    l$set_physiology(
      root_network = root_network_from_carbon(1.0, soil_depth = 1.0),
      PPFD = 3000, psi_soil = 0.1, soil_depth = 1.0,
      leaf_specific_conductance_max = 6.3e-3, atm_vpd = vpd, ca = ca,
      leaf_temp = 25, atm_o2_kpa = 21, atm_kpa = 101.3)
    km <- l$km_
    # `gamma_` is a mole fraction in umol mol^-1 and `km_` is already a partial
    # pressure in Pa, so the two have to be put on one axis before they are
    # added. C++ multiplies by `umol_per_mol_to_Pa_`, which is not exposed to R;
    # it is atmospheric pressure in Pa over 1e6, so `atm_kpa/1000` is the same
    # number and moves with the driver rather than being pinned at 0.1013.
    gamma_star <- l$gamma_ * 101.3 / 1000
    xi <- sqrt(beta * (km + gamma_star) / ratio)
    rr <- xi / sqrt(vpd * 1000)
    foc <- (gamma_star / ca + rr) / (1 + rr)
    # Loose against the C++ tolerance of 1e-11, because R is solving co-limited
    # and net of respiration where the C++ test drives the leaf onto the
    # Rubisco-limited branch. 5% is the size of that co-limitation gap.
    expect_equal(r$ci / ca, foc, tolerance = 0.05)
  }
})

test_that("LeastCost refuses a cleared, zero or negative beta", {
  sp <- leaf_supply_singlelayer()
  solve_at <- function(...) {
    leaf_solve(psi_soil = 0.1, PPFD = 900, atm_vpd = 2, supply = sp,
               model = "LeastCost", ...)
  }
  # Unset takes the aligned default (#156) ...
  expect_identical(solve_at(), solve_at(LeastCost_beta =
                                          .level_defaults[["LeastCost_beta"]]))
  # ... and a beta cleared to NA on the leaf is refused rather than solved as a
  # NaN objective. The message has to carry the units: a beta derived against
  # this package's kg H2O and umol CO2 rather than mol and mol is out by 5.6e7
  # and would solve to something entirely plausible.
  cleared <- leaf_model(supply = sp)
  cleared$set_model("LeastCost", "stem")
  cleared$LeastCost_beta <- NA_real_
  set_drivers(cleared, psi_soil = 0.1, PPFD = 900, atm_vpd = 2)
  expect_error(cleared$optimise(), "needs LeastCost_beta set")
  expect_error(cleared$optimise(), "DIMENSIONLESS")
  # ⚠️ ZERO IS REFUSED, WHICH `TF24_floor` DOES NOT DO. Zero is meaningful there
  # -- the curve is `TF24`. Here it deletes the capacity term, leaving `A/E`,
  # which is water-use efficiency and rises without bound as transpiration
  # falls, so the optimiser returns the wet bound for every driver set.
  expect_error(solve_at(LeastCost_beta = 0), "strictly positive")
  expect_error(solve_at(LeastCost_beta = -1), "strictly positive")
})

test_that("LeastCost_beta is refused by every other curve, and named when it is", {
  sp <- leaf_supply_singlelayer()
  expect_error(
    leaf_solve(psi_soil = 0.1, PPFD = 900, atm_vpd = 2, supply = sp,
               model = "TF24", LeastCost_beta = 146),
    "does not read it")
  # ⚠️ AND IT IS NOT DESCRIBED AS A PRICE. The other two model-owned slots are
  # marginal values of water; this one is a ratio of unit costs, and least-cost's
  # price of water is emergent from it rather than equal to it.
  expect_error(
    leaf_solve(psi_soil = 0.1, PPFD = 900, atm_vpd = 2, supply = sp,
               model = "TF24", LeastCost_beta = 146),
    "unit-cost ratio")
})

test_that("LeastCost_beta is in the gradient enumeration, owned by one curve", {
  expect_true("LeastCost_beta" %in% gradient_par_names())
  # Owned, so every other model refuses it -- the same discipline `CF77_lambda_`
  # and `TF24_floor_lambda_o` get.
  expect_error(
    leaf_gradient(psi_soil = 1.0, PPFD = 900, pars = "LeastCost_beta",
                  model = "TF24"),
    "LeastCost")
})

test_that("LeastCost_beta is differentiated, and the R reference agrees with its FD", {
  # ⚠️ THE SLOT HAS TO REACH THE LEAF ON EVERY PERTURBATION. C++'s
  # `gradient::apply()` writes `LeastCost_beta` from theta; R's `.gradient_setter`
  # must do the same, or the R reference re-solves at the seated beta on every
  # step and reports an EXACT ZERO -- finite, plausible, and wrong. The first
  # version of this curve did exactly that. The sign is the check that the write
  # lands where the objective reads it: a dearer capacity (larger beta) buys more
  # water per unit carbon, so the leaf opens and A rises.
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  mk <- function() {
    l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())
    l$LeastCost_beta <- 146
    l
  }
  g <- function(meth) {
    leaf_gradient(psi_soil = 0.5, PPFD = 1500, atm_vpd = 2,
                  leaf_specific_conductance_max = 6.3e-3,
                  root_network = net, x = mk(), traits = tr,
                  pars = "LeastCost_beta", model = "LeastCost", method = meth)
  }
  ift <- g("ift")
  fd <- g("fd")
  expect_identical(ift$model, "LeastCost")
  expect_identical(ift$status, "interior")
  expect_gt(ift$gradient["LeastCost_beta", "A"], 0)
  for (out in c("A", "gc", "psi_stem")) {
    expect_equal(ift$gradient["LeastCost_beta", out],
                 fd$gradient["LeastCost_beta", out],
                 tolerance = 5e-3, label = out)
  }
})

test_that("the closed form is refused for LeastCost, and says why", {
  # ⚠️ NOT A MISSING PIECE. Prentice et al. DO publish a closed form, but they
  # obtain it by solving the condition in `ci`, where the `A` in lambda cancels.
  # The closed form in this package inverts a lambda into a Medlyn xi and then
  # into a potential, and the log benefit link makes that a fixed point rather
  # than an inversion -- the same wall SOX and JW26 hit.
  sp <- leaf_supply_singlelayer()
  # ⚠️ THE REFUSAL ARRIVES AT THE SOLVE, NOT AT `set_model()`. Selecting the
  # method is legal; it is evaluating the objective through it that has nowhere
  # to go. Expecting it one call earlier is how this test first passed while
  # asserting nothing.
  l <- leaf_model(leaf_traits(), leaf_control(), sp)
  l$LeastCost_beta <- 146
  l$set_model("LeastCost", "stem", "closed")
  expect_error(l$optimise(), "closed form does not exist")
})
