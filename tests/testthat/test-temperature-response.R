# The temperature-response parameters, now reachable from R.
#
# WHY THIS FILE EXISTS. These fifteen were public C++ members with no entry in
# inst/RcppR6_classes.yml, so from R the SHAPE of the temperature response was
# fixed and unreachable. The eight derived quantities (vcmax_, jmax_, gamma_,
# ko_, kc_, km_, R_d_, electron_transport_) were bound and settable, which made
# the gap easy to miss: an R caller could write vcmax_ and see it take effect for
# exactly one solve, then be silently overwritten by the next set_physiology.
#
# The point of binding them is calibration -- a temperature-response study cannot
# use a model whose thermal optimum it cannot move -- so what these tests assert
# is not "the field exists" but "changing it changes the answer, through the
# documented workflow".

test_that("the temperature-response parameters are readable at their defaults", {
  l <- leaf_model(traits = leaf_traits(), control = leaf_control(),
                  supply = leaf_supply_single())

  expect_equal(l$vcmax_ha_, 60000)
  expect_equal(l$vcmax_H_d_, 200000)
  expect_equal(l$vcmax_d_S_, 650)
  expect_equal(l$jmax_ha_, 30000)
  expect_equal(l$jmax_H_d_, 200000)
  expect_equal(l$jmax_d_S_, 650)
  expect_equal(l$R_d_25, 1.44)
  expect_equal(l$rd_q10_intercept_, 3.09)
  expect_equal(l$rd_q10_slope_, 0.0430)

  for (f in c("gamma_25_", "gamma_ha_", "kc_25_", "kc_ha_", "ko_25_", "ko_ha_")) {
    expect_true(is.finite(l[[f]]), info = f)
  }
})

# The number that motivated binding these at all. Nothing else in the package
# reports a thermal optimum, so a user calibrating against a temperature response
# had no way to discover that the defaults sit 4-7 C below where a fitted study
# puts them -- Jones et al. (2026) fit 35 C for both.
test_that("the defaults imply the documented thermal optima", {
  l <- leaf_model(traits = leaf_traits(), control = leaf_control(),
                  supply = leaf_supply_single())
  topt <- function(ha, H_d, d_S) H_d / (d_S - 8.314 * log(ha / (H_d - ha))) - 273.15

  expect_equal(topt(l$vcmax_ha_, l$vcmax_H_d_, l$vcmax_d_S_), 31.2, tolerance = 0.05)
  expect_equal(topt(l$jmax_ha_, l$jmax_H_d_, l$jmax_d_S_), 27.9, tolerance = 0.05)
})

# ⚠️ THE TRAP THIS TEST EXISTS FOR (hazard 10, in its temperature form).
# set_physiology's temperature block is a cache keyed on (leaf_temp_,
# atm_o2_kpa_) and on NOTHING ELSE. So the obvious workflow -- change the
# parameter, call set_drivers() again -- takes a cache HIT and runs the whole
# study at the response the object was born with. Every number it produces is
# plausible. set_traits() ends in setup_clean_leaf(), which clears the cache;
# that is why it is the documented route and why this test drives both.
test_that("changing a temperature parameter changes the solve, via set_traits", {
  tr <- leaf_traits()
  build <- function() leaf_model(traits = tr, control = leaf_control(),
                                 supply = leaf_supply_single())
  drive <- function(x) {
    set_drivers(x, psi_soil = 0.5, PPFD = 1200,
                root_network = series_resistance(50),
                leaf_specific_conductance_max = 1e-4, atm_vpd = 1.0, ca = 40,
                leaf_temp = 35, atm_o2_kpa = 20.9, atm_kpa = 101.325)
    x$find_root_collar_psi()
    x$assim_colimited_
  }

  base <- drive(build())

  # d_S = 641 moves the Vcmax optimum from 31.2 C to ~35.5 C, so a solve AT 35 C
  # must gain: the leaf is now sitting near its optimum instead of past it.
  shifted <- build()
  shifted$vcmax_d_S_ <- 641
  set_traits(shifted, tr)                     # clears the temperature cache
  moved <- drive(shifted)

  expect_true(is.finite(base) && is.finite(moved))
  expect_false(isTRUE(all.equal(base, moved)))
  expect_gt(moved, base)
})

test_that("R_d_25 reaches the solve, and rises with temperature", {
  tr <- leaf_traits()
  drive <- function(rd_25, temp = 25) {
    l <- leaf_model(traits = leaf_traits(R_d_25 = rd_25),
                    control = leaf_control(), supply = leaf_supply_single())
    set_drivers(l, psi_soil = 0.5, PPFD = 1200,
                root_network = series_resistance(50),
                leaf_specific_conductance_max = 1e-4, atm_vpd = 1.0, ca = 40,
                leaf_temp = temp, atm_o2_kpa = 20.9, atm_kpa = 101.325)
    l$find_root_collar_psi()
    c(A = l$assim_colimited_, R_d = l$R_d_)
  }

  # More respiration, less net assimilation. Sabot's data imply 0.44 to 2.90
  # umol m^-2 s^-1 across their 16 species against the 1.44 default, so this is the
  # range a calibration actually needs, not an extreme.
  expect_lt(drive(2.90)[["A"]], drive(0.44)[["A"]])
  # And the trait is the value AT 25 C: it is used verbatim there and larger above.
  expect_equal(drive(0.525)[["R_d"]], 0.525)
  expect_gt(drive(0.525, temp = 40)[["R_d"]], 0.525)
})
