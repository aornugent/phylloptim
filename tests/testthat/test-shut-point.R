# A shut operating point: stomata closed, no flux, no gross assimilation, and
#
#   profit = -R_d - C(psi_stem)
#
# ⚠️ THERE ARE TWO OF THEM AND THEY ARE NOT THE SAME POINT, which is the whole
# reason a consumer cannot treat "zero flux" as "zero rows".
#
#   hydraulic shutdown  the stem is held at psi_crit and every layer is zeroed,
#                       so the environment rows really are exactly zero
#   shade death         BOTH potentials are placed at the collar of zero uptake,
#                       which IS the wet bound -- so profit reads the soil
#                       through it, and the per-layer consumptions are not zero
#                       either: they sum to zero
#
# The rows for the second are therefore -C'(B) * dB/du, and both halves of that
# product are refereed elsewhere in this suite -- the cost row against a rebuilt
# curve, the bound row against a differenced solve. What is checked here is the
# composition, against a difference of the profit at a RE-SOLVED point.

shut_profile <- c(2.6, 2.8, 3.0, 3.2, 3.4)

shut_leaf <- function(psi_soil = shut_profile, PPFD = 10) {
  l <- leaf_model(traits = leaf_traits())
  set_drivers(l, psi_soil = psi_soil, PPFD = PPFD)
  l$find_root_collar_psi()
  l
}

test_that("a shaded leaf places both potentials at the collar of zero uptake", {
  l <- shut_leaf()
  expect_identical(l$operating_point_kind_name(), "shade-death")
  expect_true(l$zero_flux_operating_point())

  # The identity the rows are built on, asserted rather than assumed: the placed
  # potential IS the wet bound, to the bit.
  wet <- l$find_root_psi(min(shut_profile), shut_profile, 0L)
  expect_identical(l$opt_psi_stem_, wet)
  expect_identical(l$opt_root_psi_, wet)

  # And the profit is respiration plus the cost at that potential, with no
  # assimilation in it at all.
  cost <- leaf_traits()$cost_scale_TF24 *
    (1 - l$proportion_of_conductivity(wet))^leaf_traits()$beta2
  expect_equal(l$profit_, -l$R_d_ - cost)
  expect_equal(l$assim_colimited_, -l$R_d_)
  expect_equal(l$transpiration_, 0)
})

test_that("a shut point's profit does not read light, exactly", {
  # Radiation decides WHICH branch the leaf is on and then appears nowhere in the
  # profit. So the light row is zero structurally rather than numerically, and
  # what says so is that the value is bit-identical across a change of light --
  # not that a difference came back small.
  a <- shut_leaf(PPFD = 20)
  b <- shut_leaf(PPFD = 5)
  expect_identical(a$operating_point_kind_name(), "shade-death")
  expect_identical(b$operating_point_kind_name(), "shade-death")
  expect_identical(a$profit_, b$profit_)
  expect_identical(a$opt_psi_stem_, b$opt_psi_stem_)
  expect_identical(a$soil_consumption_, b$soil_consumption_)
})


test_that("a hydraulically shut leaf is the other point, and reads no soil", {
  # The contrast, and the reason the stand differences these rows rather than
  # declaring them: this branch holds the stem at psi_crit and zeroes every
  # layer, so the same "zero flux" classification carries the opposite rows.
  dry <- rep(leaf_traits()$psi_crit + 0.5, 5)
  l <- shut_leaf(dry, PPFD = 1000)
  expect_identical(l$operating_point_kind_name(), "hydraulic-shutdown")
  expect_true(l$zero_flux_operating_point())
  expect_equal(l$opt_psi_stem_, leaf_traits()$psi_crit)
  expect_true(all(l$soil_consumption_ == 0))
  expect_equal(l$E_up_, 0)

  # Its profit is the cost at psi_crit, which no soil potential moves -- so the
  # soil rows here are exactly zero where the shaded leaf's are not.
  before <- l$profit_
  moved <- dry; moved[[1]] <- moved[[1]] + 0.25
  after <- shut_leaf(moved, PPFD = 1000)
  expect_identical(after$operating_point_kind_name(), "hydraulic-shutdown")
  expect_identical(after$profit_, before)
})
