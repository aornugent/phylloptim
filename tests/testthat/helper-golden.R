# Helpers for the R-side golden check. See test-golden.R for why it exists.

# This used to restate the C++ default constructor's trait values, because stage
# 1 bound only the 19-argument constructor and the R side had to name them all.
# It no longer does, and that is the improvement: `leaf_model()` carries the
# package's own defaults now, so solving through it here is what CHECKS those
# defaults still equal `phylloptim::Leaf`'s. If the two ever drift apart the bit-exact
# comparisons in test-golden.R fail, rather than quietly measuring a different
# model than the one the golden file was generated from.
default_leaf <- function() leaf_model()

# --- how exactly to compare -------------------------------------------------
#
# The golden file is bit-exact ONLY on the platform that generated it,
# macOS/arm64, and this is not a shortcut being taken -- libm's exp/pow are not
# bit-reproducible between glibc on x86-64 and Apple's libm on arm64, and FMA
# contraction differs too, so cross-platform bit-equality was never achievable.
# tests/cpp/test_golden.cpp has had a `--cross-platform` mode for this since it
# was written; the R side needs the same policy, and got it the hard way -- the
# first version of this file asserted bit-exactness everywhere, passed on macOS,
# and failed on Linux CI at the 16th significant figure.
#
# Tolerances are the C++ suite's, measured over the full grid rather than
# guessed, and they are per-CLASS because the fields split into two groups three
# orders of magnitude apart:
#
#   profit          the maximum itself, and so well conditioned
#   everything else evaluated at the ARGMAX of a flat maximum, so an error dp in
#                   the profit value displaces its location by sqrt(dp/k)
#
# If you need a magnitude, read those numbers, not a failure list.
golden_bit_exact_platform <- function() {
  si <- Sys.info()
  identical(unname(si[["sysname"]]), "Darwin") &&
    grepl("arm64|aarch64", unname(si[["machine"]]))
}

golden_tolerance <- function(field) {
  if (identical(field, "profit")) 1.0e-5 else 1.0e-3
}

# ⚠️ THE RECORDED GRADIENTS NEED THEIR OWN, LOOSER TOLERANCE, and borrowing the one
# above was wrong -- it just happened not to fail until #41 added a column.
#
# A gradient here is a FINITE DIFFERENCE of the solve, so cross-platform it carries
# the solve's own ~1e-9 floor DIVIDED BY THE STEP. The step is RELATIVE, so the
# amplification is set by the differentiated parameter's own magnitude, and the
# smallest parameter in the file decides the tolerance for all of it.
#
# Measured on Linux CI: `dA/dR_d_25` at interior-1layer disagrees by 1.3e-03
# relative, against 1.3e-04 for every column that was here before. `R_d_25` is 1.44
# where `vcmax_25` is 96 -- a 67x smaller absolute step at the same relative one --
# and the arithmetic closes: 1.3e-03 x 0.2487 x 2h = 9.3e-10, i.e. the floor itself.
#
# Confirmed independently by the step sweep: at `step = 1e-7` the `R_d_25` column has
# already moved 1.3% while `vcmax_25` has moved 0.2%. That column goes
# noise-dominated first, and it is the one that sets this number.
#
# 5e-03 leaves ~4x headroom over the one observation there is. It cannot hide a real
# change: #41's own reallocation moved these cells by 16% to 250%, two orders above
# this. ⚠️ It IS one CI observation, so read the worst-difference line the test
# prints on every run rather than assuming the headroom is still there.
gradient_golden_tolerance <- function() 5.0e-3

# expect_identical where it can hold, expect_equal with the measured per-field
# tolerance where it cannot. Wrapped in one place so no individual test has to
# decide, and so the reason lives next to the policy rather than at every call.
expect_golden <- function(actual, expected_hex, field, label,
                          tolerance = golden_tolerance(field)) {
  expected <- as.numeric(expected_hex)
  if (golden_bit_exact_platform()) {
    expect_identical(actual, expected, label = paste(label, field))
  } else {
    expect_equal(actual, expected, tolerance = tolerance,
                 label = paste(label, field))
  }
}

# The fixed drivers from tests/cpp/test_golden.cpp, which took them in turn from
# plant's tests/testthat/test-leaf.r.
golden_theta <- 0.000157
golden_ks <- 1.0
golden_h <- 5.0
golden_area_leaf <- 0.05
golden_ca <- 40.0
golden_o2 <- 21.0
golden_tleaf <- 25.0
golden_patm <- 101.3
golden_kmax <- golden_ks * golden_theta / golden_h

# The single-layer root network the C++ grid's one-layer rows use: all the root
# carbon in one 1 m layer, through the architecture model (#33). Named rather than
# rebuilt at each call site so the three tests that use it cannot drift apart.
golden_network_1m <- root_network_from_carbon(1.0 / golden_area_leaf,
                                             soil_depth = 1.0)

# One grid point, set up exactly as test_golden.cpp's solve() does: the soil
# profile spread over `layers` equal 1 m layers drying with depth, and the root
# carbon split evenly. Returns the same nine outputs the golden file records.
#
# Deliberately calls $set_physiology() directly rather than going through
# set_drivers(). This is the tie-back to the C++ baseline, so it wants as little
# of our own R between it and the binding as possible; test-surface.R covers
# set_drivers() on its own and checks the two agree on the same point.
#
# Since #33 the carbon -> resistance step is a separate call. That is the point of
# the change rather than a cost of it: this file now exercises exactly the two
# steps the C++ grid does, root_network_from_carbon() then $set_physiology(), so
# the R route to the golden numbers goes through the same public helper a user
# would call.
golden_solve <- function(psi_soil, ppfd, vpd, layers) {
  l <- default_leaf()

  i <- seq_len(layers) - 1L
  ps <- psi_soil + 0.25 * i
  depth <- 1.0 * (i + 1)
  root <- rep(1.0 / layers / golden_area_leaf, layers)

  l$set_physiology(
    root_network = root_network_from_carbon(root, soil_depth = depth),
    PPFD = ppfd,
    psi_soil = ps,
    soil_depth = depth,
    leaf_specific_conductance_max = golden_kmax,
    atm_vpd = vpd,
    ca = golden_ca,
    leaf_temp = golden_tleaf,
    atm_o2_kpa = golden_o2,
    atm_kpa = golden_patm
  )
  l$find_root_collar_psi()

  consumption <- l$soil_consumption_
  list(
    psi_stem = l$opt_psi_stem_,
    opt_root_psi = l$opt_root_psi_,
    ci = l$ci_,
    assim = l$assim_colimited_,
    transpiration = l$transpiration_,
    gc = l$stom_cond_CO2_,
    profit = l$profit_,
    e_up = l$E_up_,
    uptake = sum(consumption[is.finite(consumption)])
  )
}
