# The R layer's regression baseline: four points from the C++ golden grid,
# reproduced through the R API and compared BIT-EXACTLY.
#
# WHY THIS FILE EXISTS. tests/cpp/ is the regression baseline for the model, and
# it is completely blind to the R layer -- it never loads it. A mistranslation in
# inst/RcppR6_classes.yml (two arguments transposed, a field bound to the wrong
# member, a vector silently truncated) produces a green C++ suite and plausible
# numbers in R. Nothing else in this package would notice. So the R side needs
# its own tie-back to the recorded values, and comparing to anything less than
# the full 17 digits would let exactly the errors this is for slip through.
#
# WHY THE EXPECTED VALUES ARE HEX FLOATS, and why you must not "tidy" them into
# decimals. R's decimal string-to-double conversion is not correctly rounded:
# as.numeric, scan and read.delim all share it, and on a random sample about 18%
# of full-precision inputs come back one ULP away from the correctly rounded
# value. Pasting the golden file's %.17g strings in here would therefore make
# roughly a fifth of these comparisons fail against a model that is exactly
# right. That is not a hypothetical -- it is issue #13, where a supposed 1-ULP
# disagreement between this package and plant turned out to be entirely the
# measuring instrument: all 2352 finite values were bit-identical, and 345 of
# the 585 apparent mismatches were the parser. R reads C99 hex exactly (4000 of
# 4000 round-trip, against 3265 of 4000 for %.17g), so hex is what goes here.
#
# To regenerate after a DELIBERATE change to the golden file:
#
#   cc -O2 -o /tmp/tsv_to_hex tests/validate/tsv_to_hex.c
#   awk -F'\t' 'NR==1 || ($4==25 && $1==2 && $2==900 && $3==2 && $5==3)' \
#     tests/cpp/golden/operating_points.tsv | /tmp/tsv_to_hex
#
# and paste the row. The C tool parses with the C library's strtod, which is
# correctly rounded, and prints %a.
#
# ⚠️ The `$4==25` is not decoration. The golden file's columns are psi_soil, ppfd,
# vpd, LEAF_TEMP, layers -- two temperatures, so a tab-anchored `grep` on position
# matches leaf_temp where it means layers and pulls the wrong row. Select on the
# named column, not on position in a remembered layout.
#
# Last regenerated for the (P50, c) reparameterisation, which moved 30 of these 36
# values by ~1e-08 relative. The six that did not are the shut-down row's four
# exact zeros plus its ci and assim: that row never reaches the collar solve, so
# nothing it reports depends on the vulnerability curve -- except psi_stem, which
# IS psi_crit and moves with it. That split is why the row is here.
#
# Before that, #92 (the indexed knot grid) moved 27, the psi_stem_to_ci tolerance
# 18, and the collar root-find 27.

# Four points, chosen to exercise different parts of the model rather than to
# sample the grid evenly:
#
#   1 layer          the single-layer path, which is also the one the
#                    optimise_psi_stem_* solvers require
#   3 and 5 layers   the multi-layer supply path, where the R side has to get a
#                    whole vector across in the right order -- a transposed or
#                    reversed soil profile is invisible at one layer
#   the last row     a SHUT-DOWN point (psi_soil 6 MPa, low light, high VPD):
#                    transpiration and gc are exactly zero, assimilation is
#                    negative and the collar sits at psi_stem. This is the
#                    branch #15 fixed, and it is the one where an output a code
#                    path declines to write becomes the previous leaf's value.
#                    Pinning it here means the R layer is checked on the exit
#                    that has actually gone wrong before, not just the happy one.
golden_rows <- list(
  list(
    inputs = list(psi_soil = 0.5, ppfd = 1500, vpd = 0.5, layers = 1L),
    expected = list(
      psi_stem      = "0x1.709c089a30e11p+1",
      opt_root_psi  = "0x1.066f36de62f92p+1",
      ci            = "0x1.9f4d94f1ca6b8p+4",
      assim         = "0x1.215cb35114283p+4",
      transpiration = "0x1.450056e104a38p-16",
      gc            = "0x1.0b2aebfd57a1cp-3",
      profit        = "0x1.07970cae93315p+4",
      e_up          = "0x1.4500579e60a2bp-16",
      uptake        = "0x1.19e2624e09a31p-10"
    )
  ),
  list(
    inputs = list(psi_soil = 2.0, ppfd = 900, vpd = 2.0, layers = 3L),
    expected = list(
      psi_stem      = "0x1.b1bde93bb54b8p+1",
      opt_root_psi  = "0x1.8217d564af11ep+1",
      ci            = "0x1.0de739697d312p+3",
      assim         = "0x1.b2de947994e03p+1",
      transpiration = "0x1.b29b590bf94b5p-18",
      gc            = "0x1.6544db916f446p-7",
      profit        = "0x1.8a7e18bae99dcp-1",
      e_up          = "0x1.b29b5c3c10411p-18",
      uptake        = "0x1.78f2e3db6c197p-12"
    )
  ),
  list(
    inputs = list(psi_soil = 4.0, ppfd = 500, vpd = 1.0, layers = 5L),
    expected = list(
      psi_stem      = "0x1.77b2b48dfae5ap+2",
      opt_root_psi  = "0x1.37e82d95ec98cp+2",
      ci            = "0x1.0bd5e89326a56p+3",
      assim         = "0x1.9da2a0b6c4983p+1",
      transpiration = "0x1.9c8a7f2e28eaap-19",
      gc            = "0x1.532136ccec936p-7",
      profit        = "-0x1.db45b8023bdd1p+1",
      e_up          = "0x1.9c8a850103713p-19",
      uptake        = "0x1.65cf771fdb8d4p-13"
    )
  ),
  list(
    inputs = list(psi_soil = 6.0, ppfd = 100, vpd = 4.0, layers = 5L),
    expected = list(
      psi_stem      = "0x1.77b2b65270b19p+2",
      opt_root_psi  = "0x1.77b2b65270b19p+2",
      ci            = "0x1.1528240b78034p+2",
      assim         = "-0x1.70a3d70a3d70ap+0",
      transpiration = "0x0p+0",
      gc            = "0x0p+0",
      profit        = "-0x1.0c4e91b92c8b1p+3",
      e_up          = "0x0p+0",
      uptake        = "0x0p+0"
    )
  )
)

test_that("R's hex parser is exact, which is what the expected values rely on", {
  # If this ever fails the rest of the file is measuring R, not the model -- so
  # check the instrument before using it. 0.1 has no exact binary form, which is
  # the point: the hex literal names the double, the decimal only approximates a
  # rule for finding it.
  expect_identical(as.numeric("0x1.999999999999ap-4"), 0.1)
  expect_identical(as.numeric("0x1.91eb851eb851fp+1"), 3.14)
})

test_that("leaf_model()'s defaults are the C++ default constructor's", {
  # tests/cpp/test_golden.cpp solves with a default-constructed phylloptim::Leaf, and
  # every comparison below goes through leaf_model(). So the golden rows are
  # already the strong form of this check. This is the cheap, legible form that
  # says which thing broke when they fail: the stem vulnerability curve is built
  # from stem_b and stem_c at construction, so proportion_of_conductivity is a
  # fingerprint of that pair.
  l <- leaf_model()
  # exp(-(2/3.898245)^2.680147), the Weibull survival at psi = 2 MPa.
  expect_equal(l$proportion_of_conductivity(2.0),
               exp(-(2.0 / 3.898245)^2.680147))
  # A trait the constructor validates rather than merely stores.
  expect_error(
    leaf_model(leaf_traits(stem_P50 = -3.4)),
    "stem_P50 must be a positive magnitude"
  )
})

for (row in golden_rows) {
  local({
    inputs <- row$inputs
    expected <- row$expected
    label <- sprintf("psi_soil=%g ppfd=%g vpd=%g layers=%d",
                     inputs$psi_soil, inputs$ppfd, inputs$vpd, inputs$layers)

    test_that(paste("the R API reproduces the golden point", label), {
      got <- do.call(golden_solve, inputs)
      for (field in names(expected)) {
        # Bit-exact on macOS/arm64, per-field measured tolerance elsewhere --
        # the same policy tests/cpp/test_golden.cpp applies, and for the same
        # reason. See expect_golden() in helper-golden.R.
        expect_golden(got[[field]], expected[[field]], field, label)
      }
    })
  })
}

# ⚠️ TWO POINTS AT 40 C, and keep at least one. Every value above is a 25 C value,
# where every reference parameter in this model is DEFINED -- so every temperature
# response is inert there by construction and a 25 C-only set of pinned values cannot
# see a response curve change at all. The C++ grid carries a second temperature for
# the same reason.
#
# One 1-layer and one 3-layer, both interior at 40 C. Regenerate with the same recipe
# as above, with `$4==40`.
golden_rows_40 <- list(
  list(
    inputs = list(psi_soil = 2.0, ppfd = 900, vpd = 2.0, layers = 1L,
                  leaf_temp = 40.0),
    expected = list(
      psi_stem      = "0x1.3c43e3dfe1b6p+1",
      opt_root_psi  = "0x1.27d39655ffb4ep+1",
      ci            = "0x1.a8bae6cee41dep+4",
      assim         = "0x1.b63ab134ea75cp-1",
      transpiration = "0x1.00e1092662d9fp-18",
      gc            = "0x1.a655d98f08b69p-8",
      profit        = "-0x1.c63097866a68p-4",
      e_up          = "0x1.00e10a54bc9aap-18",
      uptake        = "0x1.bd996962bee33p-13"
    )
  ),
  list(
    inputs = list(psi_soil = 0.5, ppfd = 1500, vpd = 0.5, layers = 3L,
                  leaf_temp = 40.0),
    expected = list(
      psi_stem      = "0x1.ca0d1f625981ap+0",
      opt_root_psi  = "0x1.8cdaba05d728ap+0",
      ci            = "0x1.1428b2c0cf1f6p+5",
      assim         = "0x1.3ba43a8d6e84fp+1",
      transpiration = "0x1.c63fc674096e2p-18",
      gc            = "0x1.756a7b3df666dp-5",
      profit        = "0x1.1563763793734p+1",
      e_up          = "0x1.c63fc67bd843ep-18",
      uptake        = "0x1.89fc33020b6fbp-12"
    )
  )
)

for (row in golden_rows_40) {
  local({
    inputs <- row$inputs
    expected <- row$expected
    label <- sprintf("psi_soil=%g ppfd=%g vpd=%g layers=%d T=%g",
                     inputs$psi_soil, inputs$ppfd, inputs$vpd, inputs$layers,
                     inputs$leaf_temp)

    test_that(paste("the R API reproduces the golden point", label), {
      got <- do.call(golden_solve, inputs)
      for (field in names(expected)) {
        expect_golden(got[[field]], expected[[field]], field, label)
      }
    })
  })
}

# The 40 C rows must not equal the 25 C ones -- otherwise the two blocks above are
# measuring one temperature twice and the coverage is still missing. Asserted rather
# than assumed, because a `leaf_temp` argument that silently failed to reach
# set_physiology would give exactly that.
test_that("the 40 C block is a different operating point from the 25 C block", {
  hot  <- golden_solve(psi_soil = 2.0, ppfd = 900, vpd = 2.0, layers = 1L,
                       leaf_temp = 40.0)
  cool <- golden_solve(psi_soil = 2.0, ppfd = 900, vpd = 2.0, layers = 1L)
  expect_false(isTRUE(all.equal(hot$assim, cool$assim)))
  expect_lt(hot$assim, cool$assim)     # past the thermal optimum at 40 C
  expect_gt(hot$psi_stem, 0)
})

test_that("a shut-down point writes every output, not just the ones it changed", {
  # Hazard 8: Leaf is a value member that plant reuses for every individual in a
  # patch, so an output a branch declines to write silently becomes the previous
  # plant's value. The golden comparison above cannot see that -- it builds a
  # fresh Leaf per point, by design -- so the check has to reuse one, and the R
  # layer is the natural place because reusing an object is what an R user does.
  l <- default_leaf()
  wet <- list(psi_soil = 0.5, ppfd = 1500, vpd = 0.5)
  dry <- list(psi_soil = 6.0, ppfd = 100, vpd = 4.0)

  run <- function(l, drivers) {
    l$set_physiology(
      root_network = golden_network_1m,
      PPFD = drivers$ppfd, psi_soil = drivers$psi_soil, soil_depth = 1.0,
      leaf_specific_conductance_max = golden_ks * golden_theta / golden_h,
      atm_vpd = drivers$vpd, ca = golden_ca, leaf_temp = golden_tleaf,
      atm_o2_kpa = golden_o2, atm_kpa = golden_patm
    )
    l$find_root_collar_psi()
    c(transpiration = l$transpiration_, gc = l$stom_cond_CO2_,
      assim = l$assim_colimited_)
  }

  run(l, wet)
  after_wet_then_dry <- run(l, dry)

  fresh <- default_leaf()
  dry_from_fresh <- run(fresh, dry)

  # The point of the test: the same drivers must give the same answer whether or
  # not this object has solved something wetter first.
  expect_identical(after_wet_then_dry, dry_from_fresh)
})

test_that("psi_soil is rejected as a signed potential, through the R layer too", {
  # #25 made every water potential a positive magnitude and asserted it at the
  # input boundary. That assertion is the only thing standing between a pre-#25
  # script -- which passed -psi_soil -- and a plausible wrong number, so it has
  # to survive the trip through Rcpp as an error rather than as a warning or a
  # silent abs().
  l <- default_leaf()
  expect_error(
    l$set_physiology(
      root_network = golden_network_1m,
      PPFD = 900, psi_soil = -2.0, soil_depth = 1.0,
      leaf_specific_conductance_max = golden_ks * golden_theta / golden_h,
      atm_vpd = 2.0, ca = golden_ca, leaf_temp = golden_tleaf,
      atm_o2_kpa = golden_o2, atm_kpa = golden_patm
    ),
    "positive magnitudes in MPa"
  )
})

test_that("lambda and g1_eff are reachable and read-only", {
  # Issue #5's explicit ask, and the first item on the companion project's list.
  # They are accessors rather than stored state (hazard 5), so the thing worth
  # pinning is that they compute from the solved state and cannot be assigned.
  l <- default_leaf()
  l$set_physiology(
    root_network = golden_network_1m,
    PPFD = 900, psi_soil = 2.0, soil_depth = 1.0,
    leaf_specific_conductance_max = golden_ks * golden_theta / golden_h,
    atm_vpd = 2.0, ca = golden_ca, leaf_temp = golden_tleaf,
    atm_o2_kpa = golden_o2, atm_kpa = golden_patm
  )
  l$find_root_collar_psi()

  expect_true(is.finite(l$lambda))
  expect_true(is.finite(l$lambda_molar))
  expect_true(is.finite(l$g1_eff))
  expect_error(l$lambda <- 1, "read-only")
  expect_error(l$g1_eff <- 1, "read-only")

  # lambda is dA/dE, so it must agree with a finite difference of the profit
  # optimum's own assimilation against its transpiration. Checking the identity
  # rather than a recorded number means this survives a deliberate change to the
  # golden file, and it is what caught the sign error #25 removed: before the
  # single representation, dE_from_soil_dpsi_collar differentiated with respect
  # to the SIGNED potential and lambda came out negative where a conductance was
  # wanted.
  expect_gt(l$lambda, 0)
})
