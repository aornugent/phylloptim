# set_traits(): replacing the traits on a leaf that has already solved.
#
# The hard part is not the assignment, it is the derived state. Three pieces go
# stale on a bare trait write -- the two pre-integrated vulnerability splines, the
# solved operating point, and vcmax_/jmax_/R_d_ behind set_physiology's
# temperature cache -- so the first test solves at the DEFAULTS before re-traiting.
# Solving first is what makes it a test rather than a coincidence.

# The golden grid's drivers, so a row here is the same operating point
# tests/cpp/test_golden.cpp and test-golden.R use. Taken from test_golden.cpp.
grid_drivers <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L) {
  theta <- 0.000157
  area_leaf <- 0.05
  list(psi_soil = psi_soil + 0.25 * (seq_len(layers) - 1),
       PPFD = ppfd,
       soil_depth = 1.0 * seq_len(layers),
       root_network = root_network_from_carbon(
         rep(1 / layers / area_leaf, layers),
         soil_depth = 1.0 * seq_len(layers)),
       leaf_specific_conductance_max = 1.0 * theta / 5.0,
       atm_vpd = vpd, ca = 40, leaf_temp = 25, atm_o2_kpa = 21,
       atm_kpa = 101.3)
}

test_that("set_traits() on a used leaf equals a leaf built with those traits", {
  # The R-side statement of the C++ suite's test_set_traits_matches_a_fresh_leaf,
  # made at the boundary a calibration loop actually crosses. Bit-exact, because
  # the two routes share no code: any derived state set_traits fails to refresh
  # -- the two vulnerability splines, or vcmax_/jmax_/R_d_ behind
  # set_physiology's temperature cache -- shows up as a difference.
  traits <- leaf_traits(vcmax_25 = 110, stem_b = 4.2, root_b = 3.5,
                        cost_scale_TF24 = 8.0)
  d <- grid_drivers(2.0)

  fresh <- leaf_model(traits)
  do.call(set_drivers, c(list(fresh), d))
  fresh$find_root_collar_psi()

  # Solved at the DEFAULTS first, so every cache is warm and pointing at the old
  # traits. Solving first is what makes this a test rather than a coincidence.
  reused <- leaf_model()
  do.call(set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()
  set_traits(reused, traits)
  do.call(set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()

  expect_identical(operating_point(reused), operating_point(fresh))
  expect_identical(reused$vcmax_, fresh$vcmax_)
  expect_identical(reused$R_d_, fresh$R_d_)
})

test_that("set_traits() requires the drivers to be set again", {
  # Not a nicety: the derived photosynthetic parameters really are unknown until
  # set_physiology runs, so the operating point must not be readable as if it
  # were still valid. setup_clean_leaf puts the NA sentinels back.
  l <- leaf_model()
  do.call(set_drivers, c(list(l), grid_drivers(2.0)))
  l$find_root_collar_psi()
  expect_true(is.finite(l$assim_colimited_))

  set_traits(l, leaf_traits(vcmax_25 = 110))
  expect_false(is.finite(l$assim_colimited_))
})

test_that("set_traits() enforces the one representation for psi (#25)", {
  # A bare field write would bypass this, which is one of the four reasons the
  # traits are not bound as settable fields.
  l <- leaf_model()
  for (p in c("psi_crit", "stem_b", "root_b", "root_psi_crit")) {
    tr <- leaf_traits()
    tr[[p]] <- -tr[[p]]
    expect_error(set_traits(l, structure(tr, class = class(leaf_traits()))),
                 "positive magnitude", label = p)
  }
})

