# The golden grid's drivers, shared by the two gradient test files.
#
# A row built with these is the same operating point tests/cpp/test_golden.cpp
# and test-golden.R use, which is what lets a gradient test say "the pinned rows"
# and mean the same rows as the rest of the suite. Taken from test_golden.cpp.
#
# In a helper rather than in test-gradient.R, where it used to live, because
# test-gradient-batch.R needs the same points to compare against and each test
# file gets its own environment -- a second copy of these drivers would be free
# to drift, and then the two files would be pinning different operating points
# while appearing to pin the same ones.
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

grid_gradient <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L, ...) {
  do.call(leaf_gradient, c(grid_drivers(psi_soil, ppfd, vpd, layers), list(...)))
}

# The same drivers as ONE observation of a leaf_batch().
#
# ⚠️ `psi_soil` and `soil_depth` are wrapped in list(), which is the multi-layer
# form: a plain numeric vector is N single-layer observations, exactly as in
# leaf_solve(). Getting that backwards is the easiest mistake to make with
# leaf_batch(), and it fails loudly rather than silently -- the network and the
# profile then disagree in length.
batch_drivers <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L) {
  d <- grid_drivers(psi_soil, ppfd, vpd, layers)
  d$psi_soil <- list(d$psi_soil)
  d$soil_depth <- list(d$soil_depth)
  d
}
