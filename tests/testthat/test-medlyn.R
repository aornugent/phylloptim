# The Medlyn (2011) USO route, which until #129 had no tests at all -- which is
# how a conductance reported 1.6x high survived in shipped code, and why the
# regression tests here are the substantive part of that fix rather than the
# two-line change.
#
# The instrument throughout is the CO2 diffusion identity
#
#     A == stom_cond_CO2_ * (ca - ci) / (atm_kpa*1000) * 1e6
#
# checked against the model's OWN outputs. That needs no reference
# implementation and no cross-model calibration, which is what makes it able to
# catch a units error: the pre-#129 code satisfied every internal consistency
# check it had while violating this by exactly the diffusion ratio.

MEDLYN_RATIO <- 1.6   # the constant the USO parameterisation is fitted with

medlyn_leaf <- function(theta = 0.5, vpd = 2, ppfd = 900, ca = 40, g0 = 0.022,
                        g1 = 2.57, ci_tol = 1e-9) {
  # ci_abs_tol matters here: it is the ONLY place in the package that reads it,
  # and it sets how exactly the identity below can be expected to hold. At the
  # 1e-3 default the residual leaves ~3e-4 of relative slack.
  l <- leaf_model(supply = leaf_supply_singlelayer(),
                  control = leaf_control(ci_abs_tol = ci_tol))
  l$g0 <- g0; l$g1 <- g1
  set_drivers(l, psi_soil = 0.5, atm_vpd = vpd, PPFD = ppfd, ca = ca,
              root_network = series_resistance(1500))
  # ⚠️ AFTER set_drivers, not before. These are working copies that
  # set_physiology overwrites from its own soil-moisture members on every driver
  # set, so assigning them first is silently undone -- and a sweep written that
  # way runs every point at the default theta. Re-drive this leaf and the values
  # go back; that is the reason the helper sets them last and does not re-drive.
  l$theta_ <- theta; l$theta_fc_ <- 0.5; l$theta_w_ <- 0.2
  l
}

diffusion_residual <- function(l) {
  implied <- l$stom_cond_CO2_ * (l$ca_ - l$ci_) / (l$atm_kpa_ * 1000) * 1e6
  implied / l$assim_colimited_
}

test_that("the coupled route reports a CO2-basis conductance (#129)", {
  # THE REGRESSION TEST. Before the fix this ratio was 1.5984-1.6004 -- the USO
  # expression is a water vapour conductance and was stored without conversion.
  # A tolerance loose enough to pass 1.6 would defeat the purpose, so it is tight
  # and ci_abs_tol is tightened to make that meaningful.
  for (vpd in c(1, 2, 3)) {
    for (theta in c(0.5, 0.35)) {
      num <- medlyn_leaf(theta = theta, vpd = vpd)
      num$solve_medlyn_ci_numerical()
      expect_equal(diffusion_residual(num), 1, tolerance = 1e-6,
                   info = sprintf("numerical, vpd=%g theta=%g", vpd, theta))
    }
  }
})

test_that("1.6 is the only divisor that satisfies the identity", {
  # Pins WHICH constant, because 1.67 is the plausible wrong answer: it is what
  # every other conductance in the package uses (H2O_CO2_stom_diff_ratio_). The
  # USO form is fitted with 1.6 and the solver residual converts with 1.6, so the
  # stored value must use 1.6 or it is inconsistent with its own ci. Asserting
  # that the alternative FAILS is what stops a later "consistency" tidy-up.
  l <- medlyn_leaf()
  l$solve_medlyn_ci_numerical()
  expect_equal(diffusion_residual(l), 1, tolerance = 1e-6)

  as_if_167 <- diffusion_residual(l) * (MEDLYN_RATIO / l$H2O_CO2_stom_diff_ratio_)
  expect_lt(as_if_167, 0.96)   # a systematic -4.2%, not a rounding difference
})

test_that("the closed form's gs and ci are consistent only as g0 -> 0", {
  # ⚠️ KNOWN APPROXIMATION, asserted so it cannot be mistaken for the bug above.
  # The analytical route takes ci from the closed form and gs from the USO
  # expression, and those two agree only where the closed form is exact. So it
  # does NOT satisfy the diffusion identity at the default g0. Measured, it runs
  # 1.084 (vpd 1, theta 0.50) to 1.268 (vpd 3, theta 0.35) -- worst where the
  # deficit is highest and the soil driest, which is the signature of a dropped
  # g0 term rather than of a units error. Reporting a Fick's-law gs instead would
  # make it exact by construction but would no longer be the USO conductance;
  # that is a modelling call, and this records today's behaviour, not a wish.
  resid <- function(g0, vpd, theta = 0.5) {
    l <- medlyn_leaf(g0 = g0, vpd = vpd, theta = theta)
    l$solve_medlyn_ci_analytical()
    diffusion_residual(l)
  }
  expect_gt(resid(0.022, 3, 0.35), 1.25)   # the approximation, at the default g0
  expect_gt(resid(0.022, 1, 0.50), 1.08)
  expect_lt(abs(resid(0, 2) - 1), 1e-6)    # exact in the limit it is derived in
  expect_lt(abs(resid(0, 3) - 1), 1e-6)
})

test_that("the closed form is the coupled solve in the g0 -> 0 limit", {
  # These two entry points are two forms of ONE model, and this is the assertion
  # that says so. The gap at the default g0 is the g0 term the closed form drops;
  # send g0 to zero and they must agree. Asserted as a LIMIT, not a tolerance at
  # the default, because at g0 = 0.022 they genuinely differ by 0.036 in chi.
  chi <- function(g0, route) {
    l <- medlyn_leaf(g0 = g0)
    if (route == "num") l$solve_medlyn_ci_numerical() else l$solve_medlyn_ci_analytical()
    l$ci_ / l$ca_
  }
  gap <- function(g0) abs(chi(g0, "num") - chi(g0, "ana"))

  expect_gt(gap(0.022), 1e-3)          # the default: genuinely different
  expect_lt(gap(0.002), gap(0.022))    # monotone in g0
  expect_lt(gap(0), 1e-6)              # and the limit closes
})

test_that("beta reaches the analytical route", {
  # It did not, so this route was invariant to soil moisture while the coupled
  # route shut the leaf down -- the two entry points disagreed about whether the
  # model has a soil-moisture response at all.
  chi_at <- function(theta) {
    l <- medlyn_leaf(theta = theta)
    l$solve_medlyn_ci_analytical()
    l$ci_ / l$ca_
  }
  chis <- vapply(c(0.50, 0.40, 0.30, 0.26), chi_at, numeric(1))
  expect_true(all(diff(chis) < 0))            # drying lowers chi
  expect_gt(chis[1] / chis[length(chis)], 2)  # and not by a rounding margin

  # At beta == 1 it must reproduce the bare-g1 closed form exactly, so the
  # substitution is continuous with the published expression rather than a
  # different model.
  l <- medlyn_leaf(theta = 0.5)
  l$solve_medlyn_ci_analytical()
  expect_equal(l$ci_ / l$ca_, l$g1 / (l$g1 + sqrt(l$atm_vpd_)), tolerance = 1e-12)
})

test_that("a zero deficit is refused at the driver, not divided by", {
  # A deficit must lie in (0, e_s(T)): zero is refused where it is seated, so
  # sqrt(D) in the USO expression never sees it.
  expect_error(medlyn_leaf(vpd = 0), "atm_vpd must lie in \\(0, e_s")
})

test_that("neither route reports a negative conductance", {
  # With g0 = 0 and the soil near the wilting point the USO expression goes
  # negative, and the numerical route's "supply never reaches demand, report the
  # closest approach" exit returned it silently: gs = -3.95e-06, A = -9.7e-04.
  # A negative conductance is not an operating point; the closed state is.
  for (route in c("num", "ana")) {
    for (theta in c(0.30, 0.22, 0.2001)) {
      l <- medlyn_leaf(theta = theta, g0 = 0)
      if (route == "num") l$solve_medlyn_ci_numerical() else l$solve_medlyn_ci_analytical()
      expect_gte(l$stom_cond_CO2_, 0)
      expect_gte(l$medlyn_model_gs_, 0)
      # And where it did shut, it shut into the class's closed state rather than
      # leaving a half-written one: respiring, no water moving.
      if (l$stom_cond_CO2_ == 0) {
        expect_identical(l$transpiration_, 0)
        expect_lt(l$assim_colimited_, 0)
      }
    }
  }
})

test_that("the soil-moisture factor is clamped to [0, 1]", {
  # Unclamped, a soil wetter than field capacity scaled g1 above its fitted value
  # and one drier than the wilting point gave a negative g1 term. Clamped, wetter
  # than field capacity is exactly field capacity, and drier than the wilting
  # point is exactly the wilting point, on both routes.
  for (route in c("num", "ana")) {
    solve <- function(theta) {
      l <- medlyn_leaf(theta = theta)
      if (route == "num") l$solve_medlyn_ci_numerical() else l$solve_medlyn_ci_analytical()
      c(ci = l$ci_, gs = l$stom_cond_CO2_)
    }
    expect_identical(solve(0.6), solve(0.5))     # theta_fc_ is 0.5 in medlyn_leaf
    expect_identical(solve(0.1), solve(0.2))     # theta_w_ is 0.2
    expect_false(identical(solve(0.35), solve(0.5)))  # and it still responds between
  }
})
