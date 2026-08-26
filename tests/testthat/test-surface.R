# The stage 2 R surface: leaf_traits / leaf_control / leaf_model / set_drivers /
# operating_point / leaf_solve.
#
# The load-bearing tests here are the ones that check the convenient path and
# the raw path agree. Everything in R/leaf-model.R is a convenience over the
# generated glue, and the failure mode of a convenience layer is not that it
# errors -- it is that it quietly reorders an argument and returns a plausible
# number. test-golden.R pins the raw path against the C++ baseline; these tests
# pin the friendly path against the raw one, so the chain reaches the golden
# file from both.

test_that("leaf_traits() and leaf_control() partition the C++ constructor", {
  # Between them they must cover the 19 constructor arguments exactly once --
  # no argument silently dropped (which would leave it at whatever the
  # positional call happens to supply) and none duplicated.
  ctor_args <- setdiff(names(formals(Leaf)), "")
  covered <- c(names(leaf_traits()), names(leaf_control()))

  expect_setequal(setdiff(ctor_args, covered), character(0))
  # `R_d_25`, `JS22_gamma` and the two `CMax` parameters are traits the constructor
  # does not take: plant's own RcppR6 bindings pin this constructor by arity, so
  # `leaf_model()` assigns the fields afterwards. The test below checks those
  # assignments really happen.
  expect_setequal(setdiff(covered, ctor_args),
                  c("R_d_25", "JS22_gamma", "CMax_a", "CMax_b",
                    "integration_rule", "integration_tol"))
  expect_length(intersect(names(leaf_traits()), names(leaf_control())), 0)

  # And the split is the one the issue asked for: tolerances on the control
  # side, physiology on the trait side.
  expect_true(all(c("GSS_tol_abs", "ci_abs_tol", "ci_niter",
                    "vulnerability_curve_ncontrol") %in% names(leaf_control())))
  expect_false(any(grepl("tol|niter", names(leaf_traits()))))
})

test_that("every trait can be read back from the object (#95)", {
  # The traits were `set_traits()` arguments and nothing else, so from R they could
  # be written and not read. A caller who needed one to compute a derived quantity
  # -- Sperry's cost normalises by `k_crit = kmax * proportion_of_conductivity(psi_crit)`
  # -- had to carry it themselves, and in one downstream study that
  # was a hard-coded 5.870283 in a probe script.
  #
  # This is a COVERAGE test on purpose: it asserts the readable set contains the
  # whole trait list rather than naming four fields, so a trait added to
  # leaf_traits() without a binding fails here instead of being noticed by the next
  # caller who needs it.
  traits <- leaf_traits(vcmax_25 = 111, stem_P50 = 3.1, TF24_beta2 = 1.7,
                        root_P50 = 3.3, a = 0.28)
  l <- leaf_model(traits)
  for (nm in names(traits)) {
    expect_true(nm %in% names(l), label = paste("Leaf binds", nm))
    expect_identical(l[[nm]], traits[[nm]], label = paste("Leaf$", nm, sep = ""))
  }

  # ⚠️ Read-only, and that is hazard 10 rather than tidiness: changing a trait means
  # rebuilding up to two vulnerability splines and clearing the solved operating
  # point, so a bare write would leave the object describing two different curves.
  # `R_d_25` and `JS22_gamma` are the documented exceptions -- both are settable
  # because plant's bindings pin the generated constructor by arity, so
  # `leaf_model()` applies them afterwards.
  #
  # ⚠️ Settable is only SAFE for these two because nothing is derived from either.
  # `JS22_gamma` is read at call time by `hydraulic_cost_JS22` -- no spline, no
  # temperature cache, no precomputation -- so a bare write cannot leave the object
  # describing two different models the way a bare `stem_P50` write would. Do not
  # read this exemption as permission for the next trait.
  for (nm in setdiff(names(traits), c("R_d_25", "JS22_gamma", "CMax_a", "CMax_b"))) {
    expect_error(l[[nm]] <- 1, "read-only", label = paste(nm, "rejects a write"))
  }

  # And the read tracks set_traits(), which is what makes it a read-back rather
  # than a second copy of the constructor arguments.
  set_traits(l, leaf_traits(stem_P50 = 3.9))
  expect_identical(l$stem_P50, 3.9)

  # The DERIVED pair is bound too, and moves with the trait rather than being
  # settable beside it -- which is the whole reason they are no longer traits.
  expect_identical(l$stem_b, 3.9 / log(2)^(1 / l$stem_c))
  expect_identical(l$psi_crit, l$stem_b * log(1 / 0.05)^(1 / l$stem_c))
})

test_that("conductance is reported to water as well as to CO2 (#56)", {
  l <- leaf_model()
  set_drivers(l, psi_soil = 2.0, PPFD = 900)
  l$find_root_collar_psi()

  # Exactly the ratio, bit for bit: one multiply on the solved value, not a second
  # derivation that could drift from it.
  expect_identical(l$gs_H2O, l$stom_cond_CO2_ * l$H2O_CO2_stom_diff_ratio_)
  expect_gt(l$gs_H2O, l$stom_cond_CO2_)   # water diffuses faster than CO2

  # Read-only: it is an accessor over solved state, so a write would be the stale-
  # state trap hazard 8 describes.
  expect_error(l$gs_H2O <- 1, "read-only")
})

test_that("the H2O:CO2 diffusion ratio is settable, and 1.67 changes nothing (#50)", {
  expect_identical(leaf_model()$H2O_CO2_stom_diff_ratio_, 1.67)

  solve_at <- function(ratio) {
    l <- leaf_model()
    l$H2O_CO2_stom_diff_ratio_ <- ratio
    set_drivers(l, psi_soil = 2.0, PPFD = 900)
    l$find_root_collar_psi()
    l
  }

  # The default must stay bit-identical to not touching the field at all.
  base <- leaf_model()
  set_drivers(base, psi_soil = 2.0, PPFD = 900)
  base$find_root_collar_psi()
  expect_identical(solve_at(1.67)$profit_, base$profit_)
  expect_identical(solve_at(1.67)$stom_cond_CO2_, base$stom_cond_CO2_)

  # ⚠️ IT REACHES THE SOLVE, and `g1_eff` does not contain the ratio at all -- it is
  # chi*sqrt(D)/(1-chi) -- so the effect arrives through ci moving. That is why the
  # offset below is not predictable from the ratio and has to be pinned.
  at_16 <- solve_at(1.60)
  expect_false(isTRUE(all.equal(at_16$profit_, base$profit_)))
  expect_gt(at_16$g1_eff, base$g1_eff)
  expect_equal(abs(at_16$g1_eff - base$g1_eff) / base$g1_eff, 0.0367,
               tolerance = 0.02)
})

test_that("leaf_model() and the raw Leaf() constructor agree", {
  # The reason leaf_model() exists is that mapping 11 traits and 4 tolerances
  # onto 15 positional slots is exactly the kind of thing that goes wrong once
  # and is never noticed. So check it against a hand-written positional call
  # with the same values, on a full solve rather than on the arguments.
  raw <- Leaf(96, 2.680147, 3.4, 2.680147, 3.4,
              1.5, 157.44, 0.30, 0.7, 0.99, 1e-3, 100, 1e-3, 1000, 7.5)
  raw$initialize_integrator(21, 1e-8)
  friendly <- leaf_model()

  net <- root_network_from_carbon(20, soil_depth = 1.0)
  for (l in list(raw, friendly)) {
    l$set_physiology(net, 900, 2.0, 1.0, golden_kmax, 2.0, 40.0, 25, 21, 101.3)
    l$find_root_collar_psi()
  }
  expect_identical(operating_point(raw), operating_point(friendly))
})

test_that("a non-default trait reaches the model through leaf_model()", {
  # The previous test would pass even if leaf_model() ignored `traits` entirely
  # and always used the defaults, because the defaults are what it compares. So
  # move one and check it lands in the right slot -- vcmax_25, which raises
  # assimilation, against stem_P50, which moves the vulnerability curve.
  base <- leaf_solve(psi_soil = 2.0, PPFD = 900)
  hi_vcmax <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                         traits = leaf_traits(vcmax_25 = 150))
  expect_gt(hi_vcmax$A, base$A)

  # One trait moves the whole curve: `stem_b` and `psi_crit` are quantiles of it
  # and follow, so a brittle stem is stated once rather than as a pair that has to
  # be kept consistent by hand.
  brittle <- leaf_model(leaf_traits(stem_P50 = 2.0))
  expect_lt(brittle$proportion_of_conductivity(2.0),
            leaf_model()$proportion_of_conductivity(2.0))

  # ⚠️ `R_d_25` NEEDS ITS OWN CASE, because it is the one trait the constructor does
  # not take: if `leaf_model()`'s assignment went away, every test above would still
  # pass while `leaf_traits(R_d_25 = )` was accepted and silently ignored.
  default_rd <- leaf_traits()$R_d_25
  expect_identical(
    leaf_solve(psi_soil = 2.0, PPFD = 900,
               traits = leaf_traits(R_d_25 = default_rd))$A,
    base$A)
  expect_lt(leaf_solve(psi_soil = 2.0, PPFD = 900,
                       traits = leaf_traits(R_d_25 = 2 * default_rd))$A,
            base$A)
  expect_gt(leaf_solve(psi_soil = 2.0, PPFD = 900,
                       traits = leaf_traits(R_d_25 = 0))$A,
            base$A)
  expect_error(leaf_traits(R_d_25 = -1), "non-negative")
  expect_error(leaf_traits(R_d_25 = NA_real_), "finite")
})

test_that("a control setting reaches the model and is not treated as a trait", {
  # vulnerability_curve_ncontrol changes the spline resolution, so it moves the
  # answer slightly without being a trait. That it moves the answer at all is
  # what shows it reached the constructor rather than being dropped.
  coarse <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                       control = leaf_control(vulnerability_curve_ncontrol = 20))
  fine <- leaf_solve(psi_soil = 2.0, PPFD = 900)
  expect_false(identical(coarse$A, fine$A))
  # Still the same model, though -- a resolution change is not a trait change.
  expect_equal(coarse$A, fine$A, tolerance = 1e-2)
})

test_that("set_drivers() agrees with a positional set_physiology()", {
  raw <- leaf_model()
  raw$set_physiology(
    root_network_from_carbon(rep(20 / 3, 3), soil_depth = c(1, 2, 3)),
    900, c(1, 1.25, 1.5), c(1, 2, 3),
    golden_kmax, 2.0, 40.0, 25, 21, 101.3)
  raw$find_root_collar_psi()

  friendly <- leaf_model()
  set_drivers(friendly, psi_soil = c(1, 1.25, 1.5), PPFD = 900,
              leaf_specific_conductance_max = golden_kmax)
  friendly$find_root_collar_psi()

  # The default set_drivers() supplies -- 1 m layers, and 20 kg C m^-2 leaf split
  # evenly through root_network_from_carbon() -- must be exactly what was passed by
  # hand above. This is also the assertion that pins the #33 default to the numbers
  # it had before, since the hand-written call is the pre-#33 one.
  expect_identical(operating_point(raw), operating_point(friendly))
})

test_that("set_drivers() keeps the signed-psi_soil rejection", {
  # #25's input-boundary assertion is the only thing between a pre-#25 script and
  # a plausible wrong number. A friendly wrapper is precisely where someone would
  # be tempted to "helpfully" abs() it.
  l <- leaf_model()
  expect_error(set_drivers(l, psi_soil = -2.0), "positive magnitudes in MPa")
  expect_error(set_drivers(l, psi_soil = c(1, -2)), "positive magnitudes in MPa")
  expect_error(leaf_solve(psi_soil = -2.0), "positive magnitudes in MPa")
})

test_that("set_drivers() rejects mismatched layer counts", {
  l <- leaf_model()
  expect_error(set_drivers(l, psi_soil = c(1, 2), soil_depth = 1),
               "one entry per soil layer")
  # A network deeper than the soil profile is an out-of-bounds read in uptake(),
  # not a wrong number, so the C++ boundary rejects it (#33).
  expect_error(
    set_drivers(l, psi_soil = c(1, 2),
                root_network = root_network_from_carbon(c(1, 2, 3),
                                                        soil_depth = 1:3)),
    "rooted layers but the soil profile has only")
  expect_error(set_drivers(l, psi_soil = 2, root_network = list(a = 1)),
               "must be a RootNetwork")
  expect_error(set_drivers(l, psi_soil = numeric(0)), "at least one layer")
})

test_that("leaf_traits() and leaf_control() reject non-scalars", {
  expect_error(leaf_traits(vcmax_25 = c(90, 100)), "single finite number")
  expect_error(leaf_traits(vcmax_25 = NA_real_), "single finite number")
  expect_error(leaf_control(ci_niter = NULL), "single finite number")
  expect_error(leaf_model(traits = list()), "must come from leaf_traits")
  expect_error(leaf_model(control = list()), "must come from leaf_control")
})

test_that("leaf_solve() reproduces the stateful path exactly", {
  stateful <- leaf_model()
  set_drivers(stateful, psi_soil = 2.0, PPFD = 900,
              leaf_specific_conductance_max = golden_kmax)
  stateful$find_root_collar_psi()

  one_call <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                         leaf_specific_conductance_max = golden_kmax)
  expect_identical(one_call[, names(operating_point(stateful))],
                   operating_point(stateful))
})

# Read the fourteen outputs the slow way -- one active binding at a time, which is
# what operating_point() did before #39 -- so the one-call C++ reader can be
# checked against it.
#
# ⚠️ `uptake` IS NOT A BINDING, it is a sum, so this helper has to derive it -- and
# `sum()` IS THE WRONG WAY TO DERIVE IT for a comparison that demands bit-equality
# with C++. R's `sum()` accumulates in `LDOUBLE`, which is 80-bit on x86-64 Linux and
# 64-bit on arm64 macOS; `Leaf::operating_point_values()` accumulates in `double`. So
# `sum()` agrees with it bit-for-bit on arm64 and can differ by an ULP on x86-64,
# for a sum of more than one term.
#
# That is exactly how it failed: this test passed on macOS and on Linux for years,
# then #92 moved the values and the three-layer case came apart on ubuntu only,
# element [10] of twelve, 4.05328878168362639e-06 against ...724e-06. Nothing about
# the code under test had changed platform behaviour -- the assertion had always been
# platform-dependent and had happened to hold.
#
# `Reduce("+", ...)` accumulates left to right in plain double addition, which is
# what the C++ loop does, so the comparison is bit-exact on every platform rather
# than relaxed on some. Preferred over widening the tolerance because a shifted
# column -- the thing this test exists to catch -- is worth catching at the last bit.
outputs_one_at_a_time <- function(l) {
  consumption <- l$soil_consumption_
  finite <- consumption[is.finite(consumption)]
  c(psi_stem = l$opt_psi_stem_,
    collar = l$opt_root_psi_,
    ci = l$ci_,
    A = l$assim_colimited_,
    E = l$transpiration_,
    gc = l$stom_cond_CO2_,
    profit = l$profit_,
    hydraulic_cost = l$hydraulic_cost_,
    E_up = l$E_up_,
    uptake = Reduce(`+`, finite, 0),
    lambda = l$lambda,
    g1_eff = l$g1_eff,
    Tleaf = l$Tleaf_,
    shadow_cost = l$shadow_cost)
}

test_that("operating_point_values() returns the fourteen fields, in that order", {
  # ⚠️ THE ORDER IS AN INTERFACE AND NOTHING IN THE TYPES ENFORCES IT. The C++
  # method returns a flat vector because that is what crosses the R boundary for
  # free (#39: twelve active bindings cost ~15 us against a ~3 us solve, one call
  # ~1.5 us), and R names its positions in .operating_point_names. A field
  # inserted on one side and not the other would shift a column and keep
  # returning plausible numbers, so it is checked rather than commented.
  l <- leaf_model()

  # Unsolved FIRST, because that is the case a test written after a solve would
  # miss: every output is the NA sentinel and soil_consumption_ is empty, so the
  # uptake sum has nothing to sum. R summed the finite layers only and the C++
  # does the same; sum(numeric(0)) is 0 on both sides.
  expect_identical(l$operating_point_values(),
                   unname(outputs_one_at_a_time(l)))
  expect_identical(names(outputs_one_at_a_time(l)),
                   phylloptim:::.operating_point_names)

  # Then solved, and then in shut-down, where the layer uptakes are written but
  # the flux outputs are zero -- the branch where "every path writes its own
  # outputs" (hazard 8) is load-bearing.
  for (drivers in list(list(psi_soil = 2.0, PPFD = 900),
                       list(psi_soil = 6.0, PPFD = 100, atm_vpd = 4.0))) {
    do.call(set_drivers, c(list(l), drivers))
    l$find_root_collar_psi()
    expect_identical(l$operating_point_values(),
                     unname(outputs_one_at_a_time(l)))
  }
  # ...and the second of those really was the shut-down branch.
  expect_identical(l$transpiration_, 0)

  # Multi-layer, so the uptake sum has more than one term to get wrong.
  ml <- leaf_model()
  set_drivers(ml, psi_soil = c(1.0, 2.0, 3.0), PPFD = 900,
              soil_depth = c(0.3, 0.6, 1.0),
              root_network = root_network_from_carbon(
                c(0.1, 0.05, 0.02), soil_depth = c(0.3, 0.6, 1.0)))
  ml$find_root_collar_psi()
  expect_identical(ml$operating_point_values(),
                   unname(outputs_one_at_a_time(ml)))
  expect_length(ml$soil_consumption_, 3L)
})

test_that("Tleaf is reported, and is not the leaf_temp driver on the PM path", {
  # The gap: on the energy-balance path leaf temperature is SOLVED per operating
  # point, and `leaf_temp` has been reinterpreted as air temperature — so before
  # this column existed the only quantity that path produces was the one thing an
  # R caller could not read. Off that path the two agree, which is why `Tleaf` is
  # a copy of the driver rather than NA.
  d <- leaf_solve(psi_soil = c(1.0, 2.0, 3.0), PPFD = 900, leaf_temp = 30)
  expect_true("Tleaf" %in% names(d))
  expect_identical(d$Tleaf, rep(30, 3L))
  # Last column, not inserted: these names are positions, and a saved output
  # would shift under an insertion.
  expect_identical(names(d)[[length(names(d)) - 1L]], "Tleaf")

  # With the energy balance on, driven by hand because the gate is a field.
  l <- leaf_model()
  l$use_energy_balance_ <- TRUE
  set_drivers(l, psi_soil = 2.0, PPFD = 900, leaf_temp = 30)
  l$find_root_collar_psi()
  op <- operating_point(l)
  expect_false(isTRUE(all.equal(op$Tleaf, 30)))
  expect_gt(op$Tleaf, 30)  # absorbing radiation, shedding only part of it
  expect_identical(op$Tleaf, l$Tleaf_)
})

test_that("operating_point() is the data.frame it replaced", {
  # #39 replaced a twelve-argument data.frame() call -- 158 us, on a function
  # called once per 3 us solve -- with a direct list-to-data.frame construction
  # at 2 us. That is only safe if the result is indistinguishable, so this
  # rebuilds the old call verbatim and demands identical(). row.names is the part
  # that would break silently: it has to be the integer 1L, not 1.0.
  l <- leaf_model()
  set_drivers(l, psi_soil = 2.0, PPFD = 900)
  l$find_root_collar_psi()

  consumption <- l$soil_consumption_
  as_written_before <- data.frame(
    psi_stem = l$opt_psi_stem_,
    collar = l$opt_root_psi_,
    ci = l$ci_,
    A = l$assim_colimited_,
    E = l$transpiration_,
    gc = l$stom_cond_CO2_,
    profit = l$profit_,
    hydraulic_cost = l$hydraulic_cost_,
    E_up = l$E_up_,
    uptake = sum(consumption[is.finite(consumption)]),
    lambda = l$lambda,
    g1_eff = l$g1_eff,
    Tleaf = l$Tleaf_,
    shadow_cost = l$shadow_cost
  )
  expect_identical(operating_point(l), as_written_before)

  # And on an unsolved leaf, where every value is a sentinel.
  expect_identical(operating_point(leaf_model()),
                   operating_point(leaf_model()))
  expect_true(all(is.na(unlist(operating_point(leaf_model())[
    c("psi_stem", "A", "gc")]))))
})

test_that("leaf_solve()'s rows are assembled in the right order", {
  # The columnwise rewrite (#39) fills a preallocated matrix by index instead of
  # rbinding one-row frames, which trades an obvious cost for a
  # non-obvious failure mode: a misaligned row would pair driver i with the
  # answer to driver j and look entirely reasonable. So every row of a vectorised
  # call is checked against its own single-row call.
  #
  # The drivers deliberately vary in a way that makes each row distinguishable
  # and non-monotonic, and they cross into shut-down, which is where a
  # last-written-value bug would show up as a duplicated row.
  drivers <- list(psi_soil = c(1.0, 6.0, 0.5, 3.0, 2.0),
                  PPFD = c(900, 100, 1500, 600, 1200),
                  atm_vpd = c(2.0, 4.0, 0.5, 3.0, 1.0))
  many <- do.call(leaf_solve, drivers)

  for (i in seq_along(drivers$psi_soil)) {
    one <- leaf_solve(psi_soil = drivers$psi_soil[[i]],
                      PPFD = drivers$PPFD[[i]],
                      atm_vpd = drivers$atm_vpd[[i]])
    row <- many[i, ]
    # The extracted row carries its position as its row name; the single-row
    # call is row 1. That is the only licensed difference, so it is removed here
    # rather than by comparing loosely.
    attr(row, "row.names") <- 1L
    expect_identical(row, one, info = paste("row", i))
  }

  # Not vacuous: the rows have to actually differ, and one of them has to be the
  # shut-down branch.
  expect_identical(anyDuplicated(many$A), 0L)
  expect_true(any(many$E == 0))
})

test_that("leaf_solve() vectorises and recycles", {
  d <- leaf_solve(psi_soil = c(0.5, 1, 2), PPFD = 900)
  expect_s3_class(d, "data.frame")
  expect_identical(nrow(d), 3L)
  expect_identical(d$psi_soil, c(0.5, 1, 2))
  expect_identical(unique(d$PPFD), 900)

  # Two vectors of the same length pair up elementwise; a scalar recycles.
  pairs <- leaf_solve(psi_soil = c(1, 2), PPFD = c(400, 1200))
  expect_identical(nrow(pairs), 2L)
  expect_identical(pairs$PPFD, c(400, 1200))

  expect_error(leaf_solve(psi_soil = c(1, 2, 3), PPFD = c(400, 1200)),
               "cannot be recycled")
})

test_that("leaf_solve() takes a list for multi-layer profiles", {
  d <- leaf_solve(psi_soil = list(c(1, 1.5, 2), c(3, 3.5, 4)), PPFD = 900)
  expect_identical(nrow(d), 2L)
  expect_identical(d$layers, c(3L, 3L))
  # psi_soil is reported as the wettest layer, so the column stays numeric.
  expect_identical(d$psi_soil, c(1, 3))
  # Drier profile, less assimilation.
  expect_lt(d$A[[2]], d$A[[1]])
})

test_that("reuse = TRUE and reuse = FALSE give identical answers", {
  # This is hazard 8 checked at the R level, and it is the reason leaf_solve()
  # is allowed to reuse one object for speed. Leaf is a value type that plant
  # reuses for every individual in a patch, and three code paths have historically
  # left an output holding the PREVIOUS solve's value. The drivers below cross
  # from a wet, productive point to a shut-down one and back, which is the
  # transition that exposed those bugs.
  drivers <- list(psi_soil = c(0.5, 6.0, 1.0, 6.0, 0.5),
                  PPFD = c(1500, 100, 900, 1500, 100),
                  atm_vpd = c(0.5, 4.0, 2.0, 4.0, 0.5))
  shared <- do.call(leaf_solve, c(drivers, list(reuse = TRUE)))
  fresh <- do.call(leaf_solve, c(drivers, list(reuse = FALSE)))
  expect_identical(shared, fresh)

  # And the sweep really does reach the shut-down branch, so the check above is
  # not vacuous -- reading a bit-identical result as evidence without checking
  # the branch was reached is a mistake this project has made before.
  expect_true(any(shared$E == 0))
  expect_true(any(shared$E > 0))
})

test_that("operating_point() reports lambda and g1_eff from the solved state", {
  d <- leaf_solve(psi_soil = 1.0, PPFD = 900)
  expect_true(is.finite(d$lambda))
  expect_true(is.finite(d$g1_eff))
  expect_gt(d$lambda, 0)

  # g1_eff is computed from the solve rather than supplied, so it must vary with
  # conditions. A flat curve here would mean it had become disconnected from the
  # solved state -- reported, but not actually reporting anything.
  sweep <- leaf_solve(psi_soil = seq(0.5, 4, length.out = 8), PPFD = 900)
  live <- sweep[sweep$gc > 0 & is.finite(sweep$g1_eff), ]
  expect_gt(nrow(live), 2)
  expect_gt(diff(range(live$g1_eff)) / mean(live$g1_eff), 0.1)
})

test_that("the supply path can be chosen, and reports which is in force", {
  expect_identical(leaf_model()$supply_kind, "multilayer")

  single <- leaf_model(supply = leaf_supply_singlelayer(gravity_head = 0.05))
  expect_identical(single$supply_kind, "single")
  expect_identical(single$single_gravity_head_, 0.05)
  # The resistance is a DRIVER now, so it is unset until set_drivers() runs -- the
  # same way psi_soil is, and the change this test exists to pin.
  expect_true(is.na(single$single_resistance_))
  set_drivers(single, psi_soil = 1.5, root_network = series_resistance(1e3))
  expect_identical(single$single_resistance_, 1e3)
})

test_that("there is no state in which the tag and the supply disagree", {
  # The footgun to design around, and the reason this is two entry points
  # rather than a settable field: assigning the tag alone would leave the other
  # path's state configured and silently ignored. So the tag must not be
  # assignable at all, and the resistance must not be settable behind the tag's
  # back either.
  l <- leaf_model()
  expect_error(l$supply_kind <- "single", "read-only")
  expect_error(l$single_resistance_ <- 1e3, "read-only")

  # Switching after the drivers are set must not leave the previous path's
  # solved state lying around to be read as if it belonged to the new one.
  set_drivers(l, psi_soil = 2.0, PPFD = 900)
  l$find_root_collar_psi()
  expect_true(is.finite(l$profit_))

  l$set_supply_single(0)
  expect_identical(l$supply_kind, "single")
  expect_false(is.finite(l$profit_))
  expect_length(l$psi_soil_, 0L)
})

test_that("the single-potential path solves, and responds to its resistance", {
  solve_at <- function(r) {
    leaf_solve(psi_soil = 1.5, PPFD = 900,
               supply = leaf_supply_singlelayer(),
               root_network = series_resistance(r))
  }
  easy <- solve_at(1e3)
  hard <- solve_at(1e4)

  expect_true(all(is.finite(c(easy$A, easy$gc, easy$psi_stem))))
  expect_gt(easy$A, 0)
  # A tenfold harder path to the collar must cost the leaf something. This is
  # the check that the resistance actually reaches the model rather than being
  # stored and ignored -- the failure mode a read-only field invites.
  expect_lt(hard$A, easy$A)
  expect_lt(hard$gc, easy$gc)

  # Water runs downhill: soil <= collar <= stem, all positive magnitudes.
  expect_gte(easy$collar, 1.5 - 1e-9)
  expect_lte(easy$collar, easy$psi_stem + 1e-9)
})

test_that("the single path refuses inputs it would otherwise ignore", {
  # Silently ignoring a soil profile someone took the trouble to pass is the
  # kind of thing that produces a plausible wrong number, so it errors.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  expect_error(set_drivers(l, psi_soil = c(1, 2)), "single value")
  # `soil_depth` is the one argument that stays multi-layer-only: this path has no
  # depth profile for anything to read.
  expect_error(set_drivers(l, psi_soil = 1, soil_depth = 1),
               "no depth profile to read")

  # But `root_network` is now USED here rather than refused, which is the
  # consistency change: it is the same argument on both paths. What is refused is a
  # network built for the OTHER path -- a vulnerability-weighted horizontal term
  # this path cannot apply, which would otherwise be silently dropped.
  expect_error(
    set_drivers(l, psi_soil = 1,
                root_network = root_network_from_carbon(20, soil_depth = 1)),
    "r_R_H_min must be empty or zero")
  expect_error(
    set_drivers(l, psi_soil = 1,
                root_network = RootNetwork(r_R_V_sum = c(1e3, 2e3))),
    "exactly one series resistance")
  expect_error(set_drivers(l, psi_soil = 1, root_network = list(a = 1)),
               "must be a RootNetwork")

  expect_error(series_resistance(0), "must be positive")
  expect_error(series_resistance(-1), "must be positive")
  expect_error(leaf_supply_singlelayer(gravity_head = -1), "non-negative")
  expect_error(leaf_model(supply = list(kind = "single")),
               "must come from leaf_supply")
})

test_that("gravity_head costs the leaf water, on the single path", {
  flat <- leaf_solve(psi_soil = 1.5, PPFD = 900,
                     supply = leaf_supply_singlelayer(),
                     root_network = series_resistance(1e3))
  uphill <- leaf_solve(psi_soil = 1.5, PPFD = 900,
                       supply = leaf_supply_singlelayer(gravity_head = 0.5),
                       root_network = series_resistance(1e3))
  expect_lt(uphill$A, flat$A)
})

test_that("atm_kpa is not decorative", {
  # 10c: the ppm-to-Pa conversion is derived from atm_kpa. This was hard-coded at
  # 101.3 until #15, and plant's own driver default is 100.5 -- which turned out
  # to move TF24 offspring production by +2.4%, against +0.10% for the entire
  # rest of the swap. So the R surface must actually pass it through, and a test
  # that only ever uses 101.3 would not notice if it stopped.
  a <- leaf_solve(psi_soil = 2.0, PPFD = 900, atm_kpa = 101.3)
  b <- leaf_solve(psi_soil = 2.0, PPFD = 900, atm_kpa = 100.5)
  expect_false(identical(a$A, b$A))
})

# ---------------------------------------------------------------------------
# The root network as an input (#33)
# ---------------------------------------------------------------------------

test_that("root_network_from_carbon() reproduces the leaf's own layer thickness", {
  # The one thing that could go silently wrong when the carbon -> resistance step
  # crossed the package boundary: dz is derived from the soil-depth profile, and
  # the vertical resistance scales with dz^2, so a caller deriving dz differently
  # would be wrong by a squared factor with nothing to catch it. This asserts the
  # R helper agrees with what the leaf computes from the same profile.
  depth <- c(0.4, 0.9, 1.5)
  l <- leaf_model()
  set_drivers(l, psi_soil = c(1, 2, 3), soil_depth = depth)
  expect_equal(l$dz_, depth[[3]] / 3)

  n <- root_network_from_carbon(rep(2, 3), soil_depth = depth)
  # r_R_V[i] = beta_R_V * dz^2 / (carbon/3)
  expect_equal(n$r_R_V, rep(9.4e3 * (depth[[3]] / 3)^2 / (2 / 3), 3))
})

test_that("root_network_from_carbon() is homogeneous of degree 1 in each beta", {
  # The claim made in leaf_traits()' and leaf_gradient()'s documentation, where it
  # is what a caller who wants a beta_R_* gradient is told to rely on. Asserted so
  # the claim cannot rot: scaling one constant scales one vector and leaves the
  # other alone.
  carbon <- c(3, 6, 1)
  depth <- 1:3
  base <- root_network_from_carbon(carbon, soil_depth = depth)
  h2 <- root_network_from_carbon(carbon, soil_depth = depth, beta_R_H = 2 * 3.4e2)
  v2 <- root_network_from_carbon(carbon, soil_depth = depth, beta_R_V = 2 * 9.4e3)

  expect_equal(h2$r_R_H_min, 2 * base$r_R_H_min)
  expect_equal(h2$r_R_V_sum, base$r_R_V_sum)
  expect_equal(v2$r_R_V_sum, 2 * base$r_R_V_sum)
  expect_equal(v2$r_R_H_min, base$r_R_H_min)
})

test_that("a hand-built RootNetwork needs only the two load-bearing vectors", {
  # The R-user argument for #33: someone with a measured or fitted series
  # resistance has no carbon profile, no layer thickness and no opinion about
  # beta_R_V. They should be able to state the resistances and solve.
  l <- leaf_model()
  set_drivers(l, psi_soil = 1.5,
              root_network = RootNetwork(r_R_H_min = 25.5, r_R_V_sum = 1410))
  l$find_root_collar_psi()
  expect_true(is.finite(l$profit_))

  # And it must agree with the carbon route that produces those same numbers --
  # 20 kg C m^-2 leaf in one 1 m layer is r_R_H_min = 25.5, r_R_V_sum = 1410.
  from_carbon <- root_network_from_carbon(20, soil_depth = 1)
  expect_equal(from_carbon$r_R_H_min, 25.5)
  expect_equal(from_carbon$r_R_V_sum, 1410)

  ref <- leaf_model()
  set_drivers(ref, psi_soil = 1.5, root_network = from_carbon)
  ref$find_root_collar_psi()
  expect_identical(operating_point(l), operating_point(ref))
})

test_that("the root network is rejected when it cannot be read safely", {
  # max_soil_layer indexes psi_soil and the per-layer gravity head directly, so a
  # network longer than the soil profile is an out-of-bounds read rather than a
  # wrong number. Before #33 the agreement came for free from a length check on
  # root carbon; now it has to be asserted where the network arrives.
  l <- leaf_model()
  expect_error(
    set_drivers(l, psi_soil = 1,
                root_network = RootNetwork(r_R_H_min = c(1, 2),
                                           r_R_V_sum = c(1, 2))),
    "rooted layers but the soil profile has only")
  expect_error(
    set_drivers(l, psi_soil = 1,
                root_network = RootNetwork(r_R_H_min = 1,
                                           r_R_V_sum = c(1, 2))),
    "must have the same length")
  expect_error(
    set_drivers(l, psi_soil = 1,
                root_network = RootNetwork(r_R_H_min = -1, r_R_V_sum = 1)),
    "finite and non-negative")
})

test_that("the default-root-network memo cannot go stale", {
  # `set_drivers()` memoises the default network because materialising one costs
  # ~58 us at the R boundary (see the note above set_drivers). It is a SIZE-ONE
  # memo keyed on `soil_depth` by identical(), so the failure mode to rule out is a
  # second call with a different profile silently reusing the first profile's
  # resistances -- a wrong number with nothing to signal it, since both are
  # plausible.
  op <- function(depth) {
    l <- leaf_model()
    set_drivers(l, psi_soil = rep(1.5, length(depth)), soil_depth = depth)
    l$find_root_collar_psi()
    list(point = operating_point(l), r_R_V_sum = l$r_R_V_sum)
  }

  shallow <- op(c(0.2, 0.4))     # dz = 0.2
  deep    <- op(c(1.0, 2.0))     # dz = 1.0
  # r_R_V scales with dz^2, so a 5x thicker layer is 25x the vertical resistance.
  expect_equal(deep$r_R_V_sum, 25 * shallow$r_R_V_sum)
  expect_false(isTRUE(all.equal(deep$point$A, shallow$point$A)))

  # Going back must recompute, not return the `deep` entry the memo now holds.
  expect_identical(op(c(0.2, 0.4))$point, shallow$point)

  # And a hit must be indistinguishable from a miss.
  expect_identical(op(c(0.2, 0.4))$point, shallow$point)

  # The same for the single-potential path's cached empty network: reused across
  # calls, and reuse must not carry state.
  s1 <- leaf_model(supply = leaf_supply_singlelayer(1e3))
  set_drivers(s1, psi_soil = 1.5); s1$find_root_collar_psi()
  first <- operating_point(s1)
  s2 <- leaf_model(supply = leaf_supply_singlelayer(1e3))
  set_drivers(s2, psi_soil = 1.5); s2$find_root_collar_psi()
  expect_identical(operating_point(s2), first)
})

test_that("series_resistance() copies its prototype rather than mutating it", {
  # It avoids `RootNetwork()`'s ~58 us C++ round trip by overwriting one field of a
  # session-cached prototype. That is only safe because R copies on write, so the
  # failure mode to rule out is one call corrupting the prototype and every later
  # call inheriting it -- which would be a wrong resistance, silently, in whichever
  # caller happened to run second.
  a <- series_resistance(1500)
  b <- series_resistance(2500)
  expect_identical(a$r_R_V_sum, 1500)
  expect_identical(b$r_R_V_sum, 2500)
  # Third call, after two different values have been through: still correct.
  expect_identical(series_resistance(1500)$r_R_V_sum, 1500)

  # Mutating a returned network must not reach the prototype either.
  a$r_R_V_sum <- 99
  expect_identical(series_resistance(1500)$r_R_V_sum, 1500)

  # And the shape is the real struct's, not a hand-written list that could drift.
  expect_identical(sort(names(a)), sort(names(RootNetwork())))
  expect_s3_class(a, "RootNetwork")
  expect_length(a$r_R_H_min, 0L)

  # It must still drive a solve identically to the constructor route it replaced.
  by_ctor <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(by_ctor, psi_soil = 1.5,
              root_network = RootNetwork(r_R_V_sum = 1500))
  by_ctor$find_root_collar_psi()
  by_helper <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(by_helper, psi_soil = 1.5, root_network = series_resistance(1500))
  by_helper$find_root_collar_psi()
  expect_identical(operating_point(by_helper), operating_point(by_ctor))
})

test_that("a prescribed CF77_lambda_ survives both re-driving calls (#96)", {
  # The issue's own reproduction, from R, which is where it was found: two calls
  # that look interchangeable in a sweep disagreed about whether the caller's
  # `CF77_lambda_` survived. `CF77_lambda_` is Sperry's prescribed marginal water cost --
  # an input -- so the answer is that neither call touches it.
  #
  # Both arms matter. `set_drivers` always kept it; asserting only that one
  # passes on the code this test exists to reject.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  expect_true(is.na(l$CF77_lambda_))          # not "never initialised"

  l$CF77_lambda_ <- 30
  set_drivers(l, psi_soil = 1.5, root_network = series_resistance(1500))
  expect_identical(l$CF77_lambda_, 30)

  set_traits(l, leaf_traits(vcmax_25 = 120))
  expect_identical(l$CF77_lambda_, 30)

  # And it is still there after a solve on the re-traited object, which is the
  # sequence a sweep actually runs.
  set_drivers(l, psi_soil = 1.5, root_network = series_resistance(1500))
  l$find_root_collar_psi()
  expect_identical(l$CF77_lambda_, 30)

  # The derived state around it is still cleared, so the fix did not widen into
  # hazard 8: set_traits must leave nothing describing the old traits.
  set_traits(l, leaf_traits(vcmax_25 = 96))
  expect_true(is.na(l$profit_))
  expect_true(is.na(l$opt_psi_stem_))
  expect_identical(l$CF77_lambda_, 30)
})

test_that("TF24_floor_lambda_o is a caller input on CF77_lambda_'s footing (#96)", {
  # The same three-part contract, because it is the same kind of field: an input
  # the caller supplies, not derived state, so re-driving and re-traiting must both
  # leave it standing while clearing everything around it.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  expect_true(is.na(l$TF24_floor_lambda_o))

  l$TF24_floor_lambda_o <- 1.5e5
  set_drivers(l, psi_soil = 1.5, root_network = series_resistance(1500))
  expect_identical(l$TF24_floor_lambda_o, 1.5e5)

  set_traits(l, leaf_traits(vcmax_25 = 120))
  expect_identical(l$TF24_floor_lambda_o, 1.5e5)
})

test_that("TF24_floor is TF24 plus a price, and reduces to each of them", {
  # ⚠️ THE SAME TWO REDUCTIONS THE C++ SUITE ASSERTS BIT-FOR-BIT, restated at the R
  # surface so that the BINDINGS are covered too: a field bound to the wrong member
  # would give a plausible answer the C++ test cannot see.
  net <- series_resistance(1e4)
  solve_at <- function(model, ..., traits = leaf_traits()) {
    l <- leaf_model(traits, leaf_control(), leaf_supply_singlelayer())
    prices <- list(...)
    for (nm in names(prices)) l[[nm]] <- prices[[nm]]
    set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
    l$set_model(model, "stem")
    l$optimise()
    operating_point(l)
  }

  # ⚠️ THE ONE-RESTRICTION REDUCTION. Nothing is set on either leaf but the price,
  # because the hydraulic half reads TF24's own traits -- which is the whole reason
  # this curve has no scale of its own.
  tf24 <- solve_at("TF24")
  floor_at_zero <- solve_at("TF24_floor", TF24_floor_lambda_o = 0)
  expect_identical(floor_at_zero, tf24)

  cf77 <- solve_at("CF77", CF77_lambda_ = 1.5e5)
  no_hydraulic <- solve_at("TF24_floor", TF24_floor_lambda_o = 1.5e5,
                           traits = leaf_traits(TF24_cost_scale = 0))
  # ⚠️ EVERY COLUMN BUT `lambda`, AND THE EXCEPTION IS THE POINT. `lambda` here is
  # `marginal_cost_water()`, i.e. TF24's own price at this operating point -- and
  # the two leaves reach the same operating point with DIFFERENT
  # `TF24_cost_scale`, so it is 0 on the left and 1.65e+05 on the right. That is
  # the diagnostic disagreeing, not the model: every reported state and flux is
  # bit-identical. Pinned rather than papered over, because this column is the one
  # a reader of this curve is most likely to quote by mistake.
  #
  # ⚠️ AND `shadow_cost`, WHICH IS A SECOND EXCEPTION OF A DIFFERENT KIND. These two
  # leaves are the same model numerically -- same state, same fluxes, bit for bit --
  # and this is the one column that separates them, because only one of them
  # DECLARES its price as a price. TF24_floor reports `lambda_o * E`; CF77 reports
  # zero, not because its cost is realised carbon but because its single `lambda`
  # merges the two readings and supplies nothing to attribute the number with. So
  # the pair below is the sharpest available statement of what the column means: it
  # is about what a curve can attribute, not about what a leaf does.
  cols <- setdiff(names(cf77), c("lambda", "shadow_cost"))
  expect_identical(no_hydraulic[cols], cf77[cols])
  expect_identical(no_hydraulic$lambda, 0)
  expect_gt(cf77$lambda, 0)
  expect_gt(no_hydraulic$shadow_cost, 0)
  expect_identical(cf77$shadow_cost, 0)

  # And with both halves live it is neither of them, which is what makes the two
  # equalities above tests rather than tautologies.
  both <- solve_at("TF24_floor", TF24_floor_lambda_o = 1.5e5)
  expect_false(isTRUE(all.equal(both$psi_stem, tf24$psi_stem)))
  expect_false(isTRUE(all.equal(both$psi_stem, cf77$psi_stem)))

  # ⚠️ THE FLOOR IS ON THE EMERGENT PRICE, and `$lambda_emergent` is where to read
  # it -- NOT `operating_point()$lambda`, which is `marginal_cost_water()` and so
  # reports the TF24 cost's price whatever curve is seated. On THIS curve the two
  # are especially easy to confuse, since the hydraulic half really is TF24's.
  emergent <- function(model, ..., traits = leaf_traits()) {
    l <- leaf_model(traits, leaf_control(), leaf_supply_singlelayer())
    prices <- list(...)
    for (nm in names(prices)) l[[nm]] <- prices[[nm]]
    set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
    l$set_model(model, "stem")
    l$optimise()
    l$lambda_emergent
  }
  expect_gt(emergent("TF24_floor", TF24_floor_lambda_o = 1.5e5), 1.5e5)

  # ⚠️ AND THE FLOOR IS A WET-END STATEMENT, NOT A STATEMENT ABOUT THE OPTIMUM. A
  # curve with no price floor can still be expensive where it settles; what
  # distinguishes the two is the limit as the POTENTIAL closes on zero -- zero for
  # every conductance-loss curve here, `lambda_o` for this one. Read through the
  # bindings, so this covers them too.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  l$TF24_floor_lambda_o <- 1.5e5
  set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
  expect_equal(l$lambda_TF24_floor(1e-6), 1.5e5, tolerance = 1e-6)
  expect_lt(l$lambda_TF24(1e-6), 1)
})

test_that("TF24_floor refuses an unset price at the R surface too", {
  l <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(l, psi_soil = 1.5, PPFD = 1500,
              root_network = series_resistance(1e4))
  l$set_model("TF24_floor", "stem")
  expect_error(l$optimise(), "needs TF24_floor_lambda_o set")
  # The message has to say why there is no default, or the obvious repair -- set it
  # to zero -- silently substitutes TF24 for the model the caller asked for.
  expect_error(l$optimise(), "TF24")
})

test_that("the prescribed prices are reachable from leaf_solve() (#132)", {
  # ⚠️ THE GAP THIS CLOSES WAS REAL AND WAS NOT HYPOTHETICAL. `$CF77_lambda_` and
  # `$TF24_floor_lambda_o` are FIELDS, not traits, so `leaf_traits()` cannot carry
  # them -- and `leaf_solve()` builds its own `Leaf` internally, so before these
  # arguments existed the two priced curves could only ever raise their own
  # refusal from the one-call surface. It survived for CF77 because its docs said
  # "build the leaf yourself"; it stopped being survivable for a curve whose ONLY
  # parameter is a price.
  one <- leaf_solve(psi_soil = 1.5, PPFD = 1500,
                    supply = leaf_supply_singlelayer(),
                    root_network = series_resistance(1e4),
                    model = "TF24_floor", TF24_floor_lambda_o = 1.5e5)
  expect_true(is.finite(one$psi_stem))

  # The same answer the stateful route gives, to the bit.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  l$TF24_floor_lambda_o <- 1.5e5
  set_drivers(l, psi_soil = 1.5, PPFD = 1500,
              root_network = series_resistance(1e4))
  l$set_model("TF24_floor", "stem")
  l$optimise()
  expect_identical(one$psi_stem, operating_point(l)$psi_stem)

  # CF77 too, which is the pre-existing half of the same gap.
  expect_true(is.finite(leaf_solve(psi_soil = 1.5, PPFD = 1500,
                                   supply = leaf_supply_singlelayer(),
                                   root_network = series_resistance(1e4),
                                   model = "CF77",
                                   CF77_lambda = 1.5e5)$psi_stem))

  # Without one, the curve's own refusal reaches the caller rather than a NaN.
  expect_error(leaf_solve(psi_soil = 1.5, PPFD = 1500,
                          supply = leaf_supply_singlelayer(),
                          root_network = series_resistance(1e4),
                          model = "TF24_floor"),
               "needs TF24_floor_lambda_o set")

  # ⚠️ A PRICE THE SEATED MODEL DOES NOT READ IS REFUSED, NOT IGNORED. Silently
  # ignoring it is how someone spends an afternoon wondering why their lambda had
  # no effect, so the refusal names the curve that does read it.
  expect_error(leaf_solve(psi_soil = 1.5, PPFD = 1500,
                          supply = leaf_supply_singlelayer(),
                          root_network = series_resistance(1e4),
                          model = "TF24_floor", CF77_lambda = 1.5e5),
               "CF77 curve's prescribed price")
  expect_error(leaf_solve(psi_soil = 2.0, PPFD = 900,
                          TF24_floor_lambda_o = 1.5e5),
               "model = \"collar\" does not read it")

  # And it survives the `reuse = FALSE` path, which builds a Leaf per row and so
  # applies the price through a different call.
  fresh <- leaf_solve(psi_soil = c(1.0, 1.5), PPFD = 1500,
                      supply = leaf_supply_singlelayer(),
                      root_network = series_resistance(1e4),
                      model = "TF24_floor", TF24_floor_lambda_o = 1.5e5,
                      reuse = FALSE)
  shared <- leaf_solve(psi_soil = c(1.0, 1.5), PPFD = 1500,
                       supply = leaf_supply_singlelayer(),
                       root_network = series_resistance(1e4),
                       model = "TF24_floor", TF24_floor_lambda_o = 1.5e5,
                       reuse = TRUE)
  expect_identical(fresh, shared)
})
