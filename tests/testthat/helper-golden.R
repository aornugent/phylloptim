# Helpers for the R-side golden check. See test-golden.R for why it exists.

# This used to restate the C++ default constructor's trait values, because stage
# 1 bound only the 19-argument constructor and the R side had to name them all.
# It no longer does, and that is the improvement: `leaf_model()` carries the
# package's own defaults now, so solving through it here is what CHECKS those
# defaults still equal `phylloptim::Leaf`'s, rather than quietly measuring a
# different model than the one the golden file was generated from.
#
# ⚠️ THAT REACHES THE TRAITS, AND leaf_control()'s NUMERICAL SETTINGS ONLY BARELY.
# A trait enters the answer, so a drifted one moves these rows past tolerance; a
# spline resolution only refines it. Measured: vulnerability_curve_ncontrol at 40
# rather than 400 moves every field here by 3.0e-06, and at 20 by 5.0e-06 --
# past `profit`'s tolerance below, but inside the argmax class's. What holds THAT
# default to the C++ one is that leaf_control() reads it from
# `Leaf::ncontrol_default`, which test-surface.R asserts.
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
# Tolerances are the C++ suite's `kCrossPlatform` (tests/cpp/test_golden.cpp),
# measured over its full grid rather than guessed, and per CLASS because the fields
# split into two groups orders of magnitude apart:
#
#   profit          the maximum itself, and so well conditioned
#   everything else evaluated at the ARGMAX, which inherits its displacement
#
# These six points are a subset of that grid and sit well inside it: on Linux CI
# their worst is 7.6e-14 (profit) and 9.6e-13 (the rest), which test-golden.R
# prints on every run. If you need a magnitude, read those numbers, not a
# failure list.
golden_bit_exact_platform <- function() {
  si <- Sys.info()
  identical(unname(si[["sysname"]]), "Darwin") &&
    grepl("arm64|aarch64", unname(si[["machine"]]))
}

golden_tolerance <- function(field) {
  if (identical(field, "profit")) 1.0e-9 else 1.0e-5
}

# The water fluxes' absolute floor, mirroring kFluxNegligible in
# tests/cpp/test_golden.cpp: a shut leaf's uptake is the polished bound's
# residue, ~1e-20 kg m^-2 s^-1 of either sign on either platform (#182), and the
# smallest flux an open leaf reports is 7.6e-12. Below the floor on both sides
# the two are the same zero. Tolerant platforms only.
golden_flux_negligible <- 1.0e-15
golden_water_flux_fields <- c("transpiration", "gc", "e_up", "uptake")

# THE RECORDED GRADIENTS HAVE THEIR OWN TOLERANCE, because their floor is not the
# solved outputs'. Even the IFT route takes its mixed partial and curvature as
# central differences of `dprofit` (`.gradient_ift`), so a gradient carries the
# solve's floor DIVIDED BY THE STEP. The step is RELATIVE, so the amplification is
# set by the differentiated parameter's own magnitude: `R_d_25` (1.44) takes a 67x
# smaller absolute step than `vcmax_25` (96), and the smallest parameter in the
# file decides the tolerance for all of it. `profit` is the smallest-magnitude
# column carrying that floor, so its relative disagreement is the largest.
#
# Measured on Linux CI, the whole-file worst: 1.19e-07. 1e-5 leaves ~80x of
# headroom and still catches what this file exists for: tightening
# `psi_stem_to_ci` from 1e-10 to 1e-13 moved these cells by up to 2.3e-3 (#172),
# and #41's reallocation by 16% to 250%.
#
# ⚠️ THE FLOOR IS THE `ci` ROOT-FIND'S TOLERANCE, NOT libm. At 1e-10 this file
# disagreed across platforms by 2.34e-3 -- the same 2.3e-3 the macOS values were
# themselves wrong by, measured against a Richardson-extrapolated difference of
# the solve. Loosen `psi_stem_to_ci` and this tolerance must be re-measured, not
# widened until it passes. Read the worst-difference line the test prints on every
# run before adding a column.
gradient_golden_tolerance <- function() 1.0e-5

# expect_identical where it can hold, expect_equal with the measured per-field
# tolerance where it cannot. Wrapped in one place so no individual test has to
# decide, and so the reason lives next to the policy rather than at every call.
expect_golden <- function(actual, expected_hex, field, label,
                          tolerance = golden_tolerance(field)) {
  expected <- as.numeric(expected_hex)
  if (golden_bit_exact_platform()) {
    expect_identical(actual, expected, label = paste(label, field))
  } else if (field %in% golden_water_flux_fields &&
             abs(actual) < golden_flux_negligible &&
             abs(expected) < golden_flux_negligible) {
    succeed()
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
# ⚠️ KEEP `leaf_temp` A PARAMETER, and keep pinning at least one point away from its
# 25 C default. Every reference parameter in this model is DEFINED at 25 C, so every
# temperature response is inert there by construction -- a 25 C-only set of pinned
# values cannot see a response curve change at all.
golden_solve <- function(psi_soil, ppfd, vpd, layers, leaf_temp = golden_tleaf) {
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
    leaf_temp = leaf_temp,
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
