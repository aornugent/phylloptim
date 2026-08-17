# A shut operating point: stomata closed, no flux, no gross assimilation, and
#
#   profit = -R_d - C(psi_stem)
#
# ⚠️ THERE ARE TWO OF THEM AND THEY ARE NOT THE SAME POINT, which is the whole
# reason a consumer cannot treat "zero flux" as "zero rows".
#
#   hydraulic shutdown  the stem is held at psi_crit and every layer is zeroed,
#                       so the environment rows really are exactly zero
#   shade death         BOTH potentials are seated at the collar of zero uptake,
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

test_that("a shaded leaf seats both potentials at the collar of zero uptake", {
  l <- shut_leaf()
  expect_identical(l$operating_point_kind_name(), "shade-death")
  expect_true(l$zero_flux_operating_point())

  # The identity the rows are built on, asserted rather than assumed: the seated
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

test_that("a shaded leaf's soil rows are the cost row times the bound's", {
  # The composition, against a difference of the profit at a re-solved point.
  # Neither factor is checked here -- both have their own file -- so a
  # disagreement is this product's rather than theirs.
  l <- shut_leaf()
  B <- l$opt_psi_stem_
  cost_row <- l$hydraulic_cost_row_values(B)
  bound <- l$bound_row_values(0L)
  expect_equal(cost_row[[1]], 1)
  expect_equal(bound[[1]], 1)

  worst <- 0
  for (j in seq_along(shut_profile)) {
    predicted <- -cost_row[[2]] * bound[[8 + j]]
    for (h in c(1e-5, 1e-4, 1e-3)) {
      up <- shut_profile; up[[j]] <- up[[j]] + h
      dn <- shut_profile; dn[[j]] <- dn[[j]] - h
      a <- shut_leaf(up); b <- shut_leaf(dn)
      # A difference across a change of branch is a difference of two different
      # functions, so the arms are required to stay on this one.
      expect_identical(a$operating_point_kind_name(), "shade-death")
      expect_identical(b$operating_point_kind_name(), "shade-death")
      worst <- max(worst, abs(predicted / ((a$profit_ - b$profit_) / (2 * h)) - 1))
    }
  }
  message(sprintf("  shaded soil rows vs a re-solved difference: worst %.2e", worst))
  expect_lt(worst, 1e-5)

  # Non-vacuity: the rows are live, and they are live because the BOUND moves --
  # a shut leaf that read no soil at all would pass every check above.
  expect_true(all(abs(bound[9:13]) > 1e-6))
  expect_gt(abs(cost_row[[2]]), 1e-6)
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
