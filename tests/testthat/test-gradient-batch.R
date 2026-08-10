# The batched trait gradient (issue #4, PLAN 11d stage 2).
#
# WHAT IS BEING PINNED HERE, in order of how much it matters.
#
#   1. That the C++ composite computes the SAME THING as `leaf_gradient()`,
#      bit-for-bit, over both supply paths, both methods and the pinned and
#      shut-down rows. This is the whole risk of stage 2: it is a second
#      implementation of an algorithm that already exists, and the two share no
#      code, so a transcription slip has nothing to hide behind. Bit-for-bit
#      rather than close, because a tolerance would let one hide inside the
#      solver's own ~1e-09 floor -- and because the composite's arithmetic order,
#      its avoidance of fused multiply-add and its explicit call sequencing were
#      all written for this test and would otherwise not be checked at all.
#
#   2. That the values are what they were. ⚠️ (1) CANNOT SEE A CHANGE APPLIED
#      CONSISTENTLY TO BOTH. Move the step rule, the stationarity threshold or a
#      solver tolerance and the two implementations move together, the equality
#      test passes, and nothing says the gradients changed. That is the trap #70
#      hit inside R alone. So `gradient_golden.tsv` records them, and its five
#      rows agree with the independently arbitrated references in
#      test-gradient.R -- which is the check that the transcription is RIGHT and
#      not merely self-consistent. ⚠️ It is bit-exact only on macOS/arm64, and
#      compared with the golden file's own per-class tolerance elsewhere, for the
#      reason helper-golden.R gives at length.
#
#   3. That a row the solve cannot handle costs THAT ROW. A proposal during a fit
#      will reach such points; taking out the dataset for one of them would take
#      out the draw.
#
#   4. That the parameter enumeration cannot drift. R passes integer POSITIONS
#      into it, so a reordering silently differentiates the wrong parameter.

# `grid_drivers()` and `batch_drivers()` are in helper-gradient.R, shared with
# test-gradient.R so both files pin the same operating points.

# One row's gradient through the batch, in `leaf_gradient()`'s shape, so the two
# can be compared field by field without every test restating the indexing.
one_row <- function(g) {
  list(gradient = g$gradient[1, , , drop = FALSE],
       value = g$value[1, ],
       status = g$status[[1]],
       method = g$method[[1]],
       H = g$H[[1]],
       stationarity = g$stationarity[[1]])
}

expect_batch_matches_r <- function(args, pars, method, label) {
  batch_args <- args
  ref <- tryCatch(
    do.call(leaf_gradient, c(args, list(pars = pars, method = method))),
    error = function(e) conditionMessage(e))
  if (!is.null(batch_args$psi_soil) && !is.list(batch_args$psi_soil)) {
    batch_args$psi_soil <- list(batch_args$psi_soil)
    if (!is.null(batch_args$soil_depth)) {
      batch_args$soil_depth <- list(batch_args$soil_depth)
    }
  }
  b <- do.call(leaf_batch, batch_args)
  got <- one_row(leaf_gradient_batch(b, pars = pars, method = method))

  if (is.character(ref)) {
    # R stops; the batch has to report the same refusal per row rather than
    # propagate it, and it must not quietly return a number instead.
    expect_identical(got$status, "error", label = label)
    return(invisible(NULL))
  }
  expect_identical(as.vector(got$gradient), as.vector(ref$gradient),
                   label = paste(label, "gradient"))
  expect_identical(unname(got$value), unname(ref$value),
                   label = paste(label, "value"))
  expect_identical(got$status, ref$status, label = paste(label, "status"))
  expect_identical(got$method, ref$method, label = paste(label, "method"))
  expect_identical(got$H, ref$H, label = paste(label, "H"))
  expect_identical(got$stationarity, ref$stationarity,
                   label = paste(label, "stationarity"))
  invisible(ref$status)
}

test_that("the parameter enumeration is the same order in R and in C++", {
  # ⚠️ THE ORDER IS AN INTERFACE AND NOTHING IN THE TYPES ENFORCES IT. R sends
  # integer positions into C++'s `phylloptim::gradient::par_names`, so appending
  # is safe and reordering would differentiate the wrong parameter and report
  # plausible numbers for it. Compared rather than trusted, in both directions.
  r_side <- c(names(leaf_traits()), "leaf_specific_conductance_max",
              "resistance")
  expect_identical(gradient_par_names(), r_side)
  # And the count, which is what a positional trait call would silently break:
  # thirteen traits then the two that are not traits.
  expect_length(gradient_par_names(), 15L)
  expect_identical(gradient_par_names()[1:13], names(leaf_traits()))
})

test_that("the batch reproduces leaf_gradient() bit-for-bit across the grid", {
  # THE ACCEPTANCE TEST. Eighteen operating points spanning wet to shut-down and
  # one to five soil layers, times all three methods -- so it covers the interior
  # rows the composite runs on, the pinned rows it must decline, the shut-down
  # rows where there is no curvature at all, and the forced-route cases where R
  # raises an error and the batch has to report one per row instead.
  #
  # `psi_crit` is in `pars` deliberately: it appears nowhere in the profit
  # function, so the composite must return exactly zero for it and the fallback
  # must not. A transcription that dropped the explicit `collar` assignment, or
  # differenced two identical numbers instead, shows up on that column.
  pars <- c("vcmax_25", "stem_b", "cost_scale_TF24", "psi_crit")
  rows <- expand.grid(layers = c(1L, 3L, 5L), psi_soil = c(0.5, 2.0, 3.0, 4.0,
                                                          6.0),
                      KEEP.OUT.ATTRS = FALSE)
  seen <- character(0)
  for (method in c("auto", "ift", "fd")) {
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      lab <- sprintf("psi_soil=%g layers=%d method=%s", r$psi_soil, r$layers,
                     method)
      # `grid_drivers()`, not `batch_drivers()`: the helper below hands the same
      # arguments to BOTH entry points and puts the profile in a list for the
      # batch itself, since `leaf_gradient()` takes a bare vector.
      st <- expect_batch_matches_r(grid_drivers(r$psi_soil, layers = r$layers),
                                   pars, method, lab)
      if (method == "auto") {
        seen <- c(seen, st)
      }
    }
  }
  # ⚠️ The grid has to REACH all three classes, or the loop above is only testing
  # the happy path while looking thorough -- which is the way a comparison test
  # fails silently. Asserted by WHERE the classes fall rather than by bare counts,
  # so the claim survives adding a row to the grid: the shut-down rows are exactly
  # the ones drier than psi_crit, and the pinned ones are at the dry end a drought
  # calibration visits.
  expect_length(seen, 15L)
  expect_identical(seen[rows$psi_soil == 6], rep("no-gradient", 3L))
  expect_true(all(rows$psi_soil[seen == "pinned"] %in% c(3, 4)))
  expect_gt(sum(seen == "pinned"), 0L)
  expect_identical(sum(seen == "interior"), 11L)
})

test_that("the batch matches leaf_gradient() on the single-potential path", {
  # The other supply path, and with it the two parameters that are not traits.
  # Both of them move the FEASIBLE COLLAR INTERVAL and not only the profit
  # function -- `resistance` directly, the conductance through the supply curve --
  # so a perturbation can push a bound past psi* more readily than a trait
  # perturbation can, and the clamp detector has to be reached through this route
  # too.
  pars <- c("vcmax_25", "stem_b", "leaf_specific_conductance_max", "resistance")
  for (psi_soil in c(0.5, 1.5, 3.5, 5.5, 6.5)) {
    for (resistance in c(1e3, 1e4)) {
      for (method in c("auto", "ift", "fd")) {
        args <- list(psi_soil = psi_soil, PPFD = 900, atm_vpd = 3.0,
                     supply = leaf_supply_single(),
                     root_network = series_resistance(resistance))
        expect_batch_matches_r(
          args, pars, method,
          sprintf("psi_soil=%g r=%g method=%s", psi_soil, resistance, method))
      }
    }
  }
})

test_that("the batch matches leaf_gradient() with the stem curve rebuilt", {
  # `fast_stem_curve = FALSE` rebuilds the vulnerability spline instead of
  # rescaling it (PLAN 11f), so it is a different code path through `apply()` and
  # not merely a slower one. It has to be transcribed too.
  for (psi_soil in c(1.5, 2.0, 4.0)) {
    for (method in c("auto", "fd")) {
      args <- batch_drivers(psi_soil)
      ref <- do.call(leaf_gradient, c(grid_drivers(psi_soil),
                                      list(pars = "stem_b", method = method,
                                           fast_stem_curve = FALSE)))
      b <- do.call(leaf_batch, args)
      got <- leaf_gradient_batch(b, pars = "stem_b", method = method,
                                 fast_stem_curve = FALSE)
      expect_identical(as.vector(got$gradient[1, , ]),
                       as.vector(ref$gradient),
                       label = sprintf("psi_soil=%g method=%s", psi_soil,
                                       method))
    }
  }
})

test_that("the recorded gradients have not moved", {
  # The guard the bit-for-bit comparisons above cannot be: they would pass a
  # change applied to BOTH implementations. See this file's header, and
  # tools/gradient_golden.R for how to regenerate -- deliberately.
  #
  # ⚠️ Hex floats, read exactly. Do not convert them to decimals: R's decimal
  # parser is not correctly rounded and about 18% of full-precision values come
  # back one ULP out (#13).
  gold <- utils::read.delim(test_path("gradient_golden.tsv"),
                           stringsAsFactors = FALSE)
  cases <- list(
    "interior-1layer" = list(args = batch_drivers(2.0),
                             pars = c("vcmax_25", "stem_b", "psi_crit")),
    "interior-5layer" = list(args = batch_drivers(0.5, vpd = 0.5, layers = 5L),
                             pars = c("vcmax_25", "stem_b", "psi_crit")),
    "pinned-dry-3layer" = list(args = batch_drivers(4.0, vpd = 0.5,
                                                    layers = 3L),
                               pars = c("vcmax_25", "stem_b", "psi_crit")),
    "shutdown-1layer" = list(args = batch_drivers(6.0),
                             pars = c("vcmax_25", "stem_b", "psi_crit")),
    "single-potential" = list(
      args = list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
                  supply = leaf_supply_single(),
                  root_network = series_resistance(1e4)),
      pars = c("vcmax_25", "leaf_specific_conductance_max", "resistance")))

  expect_setequal(unique(gold$case), names(cases))
  # ⚠️ BIT-EXACT ONLY ON THE PLATFORM THAT GENERATED IT, macOS/arm64, exactly as
  # tests/cpp/golden/ is -- and the first version of this test asserted equality
  # everywhere, passed here and failed on Linux CI, which is the same mistake
  # test-golden.R's first version made. libm's exp/pow are not bit-reproducible
  # between glibc on x86-64 and Apple's libm on arm64.
  #
  # THE TOLERANCE IS NOT A GUESS AND NOT NEW. These are derivatives of outputs
  # evaluated at the ARGMAX of a flat maximum, so they belong to the golden file's
  # own sqrt-amplified class -- and `golden_tolerance()` is that policy, reused
  # here rather than restated. Measured worst relative disagreement on Linux is
  # 1.3e-4, against the golden file's 1.4e-4 for the same class: a central
  # difference cancels the systematic part of a libm difference but not all of it.
  #
  # ⚠️ Read the SUMMARY LINE below for a magnitude, never the FAIL lines. The
  # figures in this package's guide were wrong twice because a truncated failure
  # list was read as if it were the distribution.
  #
  # ⚠️ WHICH COLUMNS ARE PINNED IS READ OUT OF THE FILE, not written here, and that
  # matters now that there are five outputs and this file records four. It was
  # generated before `profit` existed and CANNOT be regenerated here: hex floats
  # off Linux would be wrong on the platform the other columns came from. So
  # `profit` is unpinned until someone regenerates on macOS/arm64, at which point
  # `tools/gradient_golden.R` emits it and this loop picks it up with no edit.
  # Meanwhile the profit column is held by the two things this file cannot be: the
  # bit-for-bit route comparison above, and test-gradient.R's arbitration against
  # a difference of the whole solve.
  outs <- setdiff(names(gold), c("case", "status", "method", "par"))
  expect_true(all(outs %in% gradient_output_names()))
  worst <- 0
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    rows <- gold[gold$case == nm, ]
    b <- do.call(leaf_batch, cs$args)
    g <- leaf_gradient_batch(b, pars = cs$pars)
    expect_identical(g$status[[1]], rows$status[[1]], label = nm)
    expect_identical(g$method[[1]], rows$method[[1]], label = nm)
    for (i in seq_len(nrow(rows))) {
      for (out in outs) {
        got <- g$gradient[1, rows$par[[i]], out]
        expect_golden(got, rows[[out]][[i]], out,
                      paste(nm, rows$par[[i]]))
        ref <- as.numeric(rows[[out]][[i]])
        if (ref != 0) {
          worst <- max(worst, abs(got - ref) / abs(ref))
        }
      }
    }
  }
  cat(sprintf("\ngradient_golden.tsv: worst relative difference %.3g (%s)\n",
              worst,
              if (golden_bit_exact_platform()) "bit-exact platform" else
                "compared with tolerance"))
  # And the two claims the choice of rows exists to make, stated rather than left
  # implicit in the hex: `psi_crit` is exactly zero at an interior optimum
  # because it is not in the profit function, and it carries the gradient at a
  # dry-pinned one because there it IS the binding constraint.
  interior <- gold[gold$case == "interior-1layer" & gold$par == "psi_crit", ]
  expect_identical(as.numeric(interior$A), 0)
  pinned <- gold[gold$case == "pinned-dry-3layer" & gold$par == "psi_crit", ]
  expect_gt(as.numeric(pinned$A), 1)
})
test_that("the profit column agrees across routes and with a differenced solve", {
  # ⚠️ THE ONE COLUMN WHOSE REPORTED VALUE IS NOT A DIFFERENCE OF WHAT IT REPORTS,
  # and that is why it needs a THIRD reference rather than the two this file uses
  # everywhere else. Both composites take `profit` from the envelope theorem: the
  # partial at a FIXED collar, with (dprofit/dpsi)(dpsi*/dtheta) dropped because
  # dprofit/dpsi is zero at an interior optimum. Neither one differences the
  # constrained maximum it claims to be the derivative of, so two independent
  # implementations agreeing establishes the transcription and says nothing
  # whatever about the identity. `leaf_solve()` is the arbiter: it shares no code
  # with either composite beyond the model itself.
  pars <- c("vcmax_25", "stem_b", "stem_c", "cost_scale_TF24", "beta2", "root_b")

  # A central difference of the WHOLE SOLVE, at the step the gradient uses.
  solve_profit_gradient <- function(bd, p) {
    v <- leaf_traits()[[p]]
    h <- max(abs(v), 1) * 1e-6
    at <- function(x) {
      tt <- leaf_traits()
      tt[[p]] <- x
      do.call(leaf_solve, c(bd, list(traits = tt)))$profit
    }
    (at(v + h) - at(v - h)) / (2 * h)
  }

  # ⚠️ psi_soil = 2.0 IS DELIBERATELY NOT IN THIS LIST, and the reason is the
  # REFERENCE's floor rather than the composite's. There vcmax_25's profit gradient
  # is 8.2e-03 against a step of 1.5e-04, so the solve's ~1e-09 in profit bounds
  # any difference quotient of it at ~4e-04 relative -- and 3.3e-04 is what it
  # shows. The three points below have larger profit gradients and land at 6.5e-07,
  # 2.9e-07 and 8.2e-08.
  for (psi_soil in c(0.5, 1.0, 3.0)) {
    bd <- batch_drivers(psi_soil)
    cpp <- leaf_gradient_batch(do.call(leaf_batch, bd), pars = pars)
    r <- do.call(leaf_gradient, c(grid_drivers(psi_soil), list(pars = pars)))
    expect_identical(cpp$status[[1]], "interior")
    expect_identical(cpp$method[[1]], "ift")
    # The two composites on this column, bit-for-bit, as everywhere else here.
    expect_identical(cpp$gradient[1, , "profit"], r$gradient[, "profit"],
                     label = paste("routes at", psi_soil))
    for (p in pars) {
      expect_equal(cpp$gradient[1, p, "profit"], solve_profit_gradient(bd, p),
                   tolerance = 1e-5,
                   label = paste("solve difference at", psi_soil, p))
    }
  }

  # ⚠️ AND THE ENVELOPE SHORTCUT IS CONFINED TO THE COMPOSITE, which is the other
  # half of the claim and the half that could fail silently. `dprofit/dpsi` is zero
  # at an INTERIOR optimum and nowhere else, so anywhere else the fixed-collar
  # partial is not the total derivative and the profit row has to come from
  # differencing the solve. `gradient_fd()` and `.gradient_fd()` have no profit
  # special case at all; this is what says so. At all three points where the
  # fallback runs, its profit row is BIT-IDENTICAL to a central difference of
  # `leaf_solve()` -- no tolerance, because there is nothing to tolerate: it is the
  # same difference of the same solve.
  #
  # The dry-pinned row is the one a leak could not survive. `psi_crit` does not
  # appear in the profit function, so the fixed-collar partial for it is EXACTLY
  # ZERO -- while it is the binding constraint there and the true derivative is
  # 0.507. A profit special case that reached this path would report the zero, and
  # a zero is the kind of wrong answer nothing about looks suspicious.
  #
  # ⚠️ `fast_stem_curve = FALSE` so that `stem_b` is IN the comparison rather than
  # excluded from it. The shortcut rescales the stem spline instead of rebuilding
  # it, so it differentiates a marginally different model -- 4.7e-05 relative on
  # this parameter, which is PLAN 11f's own trade and has nothing to do with
  # profit.
  fd_pars <- c(pars, "psi_crit", "root_c", "jmax_25", "a")
  points <- list(pinned = batch_drivers(4.0, vpd = 0.5, layers = 3L),
                 shutdown = batch_drivers(6.0),
                 interior = batch_drivers(2.0))
  status <- c(pinned = "pinned", shutdown = "no-gradient",
              interior = "interior")
  for (nm in names(points)) {
    bd <- points[[nm]]
    g <- leaf_gradient_batch(do.call(leaf_batch, bd), pars = fd_pars,
                             method = "fd", fast_stem_curve = FALSE)
    expect_identical(g$status[[1]], status[[nm]], label = nm)
    expect_identical(g$method[[1]], "fd", label = nm)
    for (p in fd_pars) {
      expect_identical(g$gradient[1, p, "profit"], solve_profit_gradient(bd, p),
                       label = paste("fd", nm, p))
    }
  }

  # And under the default `method`, which is the one a caller gets: the dry-pinned
  # point routes itself to the fallback and reports the binding constraint's
  # derivative rather than the envelope theorem's zero.
  auto <- leaf_gradient_batch(do.call(leaf_batch, points$pinned),
                              pars = "psi_crit")
  expect_identical(auto$status[[1]], "pinned")
  expect_identical(auto$method[[1]], "fd")
  expect_equal(auto$gradient[1, "psi_crit", "profit"], 0.5074, tolerance = 1e-3)
})

test_that("an unsolvable row costs that row and not the batch", {
  # ⚠️ THE REASON THE BATCH REPORTS A STATUS RATHER THAN THROWING. A proposal
  # during a fit will reach operating points the solve cannot handle -- that is
  # what a proposal distribution does -- and taking out the whole likelihood
  # evaluation for one observation would take out the draw.
  pars <- c("vcmax_25", "stem_b")
  ok <- c(1.5, 2.0, 2.5)
  b <- leaf_batch(psi_soil = c(ok[[1]], 6.5, ok[[2]], 7.5, ok[[3]]), PPFD = 900)
  g <- leaf_gradient_batch(b, pars = pars, method = "ift")

  expect_identical(g$status, c("interior", "error", "interior", "error",
                               "interior"))
  # BIT-IDENTICAL to being differentiated alone, which is the real claim: a
  # failed row must not leave the leaf carrying state into its neighbours. The
  # displaced stem-vulnerability spline `perturb_stem_b()` leaves behind is
  # exactly the kind of thing that would show up here and nowhere else.
  solo <- vapply(ok, function(p) {
    as.vector(leaf_gradient_batch(leaf_batch(psi_soil = p, PPFD = 900),
                                  pars = pars, method = "ift")$gradient[1, , ])
  }, numeric(length(pars) * length(gradient_output_names())))
  for (i in seq_along(ok)) {
    expect_identical(as.vector(g$gradient[2 * i - 1, , ]), solo[, i],
                     label = paste("row", 2 * i - 1))
  }

  # A failed row's gradient is ALL NA rather than partially filled. The parameter
  # loop fails on the first parameter here, but even where it does not, returning
  # the entries computed before the failure alongside `status == "error"` would
  # invite a caller who checks the status per row to read a gradient with a hole
  # in it.
  expect_true(all(is.na(g$gradient[c(2, 4), , ])))
  expect_true(all(nzchar(g$message[c(2, 4)])))
  expect_identical(g$message[c(1, 3, 5)], rep("", 3L))
  # The point's own diagnostics survive, because they describe the operating
  # point rather than the derivative and were determined before the failure.
  expect_true(all(is.finite(g$stationarity[c(2, 4)]) |
                    is.infinite(g$stationarity[c(2, 4)])))
})

test_that("a batch that failed on the fast stem_b path leaves a usable leaf", {
  # ⚠️ `perturb_stem_b()` leaves the splines built at a different stem_b ON
  # PURPOSE, and `set_traits()` is the way back. A row that threw between the two
  # skipped the restore, so the batch does it in the handler -- and the LAST row
  # is the case that needs it, because there is no next row whose own
  # `set_traits()` would have covered for it (hazard 8).
  b <- leaf_batch(psi_soil = c(2.0, 7.5), PPFD = 900)
  invisible(leaf_gradient_batch(b, pars = "stem_b", method = "ift"))
  after <- leaf_gradient_batch(b, pars = "stem_b")
  fresh <- leaf_gradient_batch(leaf_batch(psi_soil = c(2.0, 7.5), PPFD = 900),
                               pars = "stem_b")
  expect_identical(after$gradient, fresh$gradient)
  expect_identical(after$value, fresh$value)

  # And the leaf itself, which the caller can still reach, solves as a fresh one
  # does -- the strongest form of the claim, since a rescaled curve would be
  # invisible in the gradient columns above if both routes carried it.
  set_drivers(b$leaf, psi_soil = 1.5, PPFD = 900)
  b$leaf$find_root_collar_psi()
  expect_identical(operating_point(b$leaf),
                   leaf_solve(psi_soil = 1.5, PPFD = 900)[
                     , names(operating_point(b$leaf))])
})

test_that("per-observation theta differentiates each row at its own parameters", {
  # THE HIERARCHICAL CASE, which is what the batch is for: `leaf-calibration` maps
  # 40 fitted parameters onto 4 model ones, so the model parameters differ per
  # observation and a shared `theta` would be the wrong shape entirely.
  pars <- c("vcmax_25", "stem_b")
  psv <- c(1.0, 1.5, 2.0)
  # ⚠️ stem_b raised, not lowered. #38: the vulnerability spline is built over
  # [0, b*log(100)^(1/c)] and never consults psi_crit, so a stem_b below ~3.4 at
  # the default psi_crit puts the solve outside its own domain and dies naming
  # neither.
  traits <- lapply(c(96, 110, 120), function(v) {
    leaf_traits(vcmax_25 = v, stem_b = 4.2)
  })

  b <- leaf_batch(psi_soil = psv, PPFD = 900)
  theta <- t(vapply(traits, function(tr) {
    unname(c(unlist(tr)[gradient_par_names()[1:13]], 3.14e-5, NA_real_))
  }, numeric(15)))
  g <- leaf_gradient_batch(b, theta = theta, pars = pars)

  for (i in seq_along(psv)) {
    ref <- leaf_gradient(psi_soil = psv[[i]], PPFD = 900, traits = traits[[i]],
                         pars = pars)
    expect_identical(as.vector(g$gradient[i, , ]), as.vector(ref$gradient),
                     label = paste("row", i))
    expect_identical(unname(g$value[i, ]), unname(ref$value),
                     label = paste("row", i))
  }

  # A one-row theta is shared by every observation, which is the other half of the
  # recycling rule and the shape the default takes.
  shared <- leaf_gradient_batch(b, theta = theta[2, , drop = FALSE], pars = pars)
  all_at_2 <- leaf_gradient_batch(b, traits = traits[[2]], pars = pars)
  expect_identical(shared$gradient, all_at_2$gradient)
})

test_that("the default theta describes the point the batch was built at", {
  # `leaf_specific_conductance_max` is both a driver and a fitted parameter, so
  # `theta` has to hold the value the solve will actually use -- the same thing
  # `.gradient_theta()` does for one observation. A default that quietly used the
  # package default instead would differentiate at a point nobody asked about and
  # report plausible numbers.
  kmax <- 5e-5
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900,
                  leaf_specific_conductance_max = kmax)
  g <- leaf_gradient_batch(b, pars = c("vcmax_25",
                                       "leaf_specific_conductance_max"))
  ref <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                       leaf_specific_conductance_max = kmax,
                       pars = c("vcmax_25", "leaf_specific_conductance_max"))
  expect_identical(as.vector(g$gradient[1, , ]), as.vector(ref$gradient))

  # Non-shared conductances take the n-row path, which is a different branch in
  # `.gradient_theta_matrix()` and has to produce the same answers.
  kv <- c(3.14e-5, 5e-5)
  b2 <- leaf_batch(psi_soil = c(1.5, 1.5), PPFD = 900,
                   leaf_specific_conductance_max = kv)
  expect_identical(b2$theta_nrow, 2L)
  g2 <- leaf_gradient_batch(b2, pars = "leaf_specific_conductance_max")
  for (i in 1:2) {
    ref2 <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                          leaf_specific_conductance_max = kv[[i]],
                          pars = "leaf_specific_conductance_max")
    expect_identical(as.vector(g2$gradient[i, , ]), as.vector(ref2$gradient),
                     label = paste("row", i))
  }
})

test_that("the result is shaped and named for a caller applying a Jacobian", {
  # The layering: C++ returns dY/dtheta for the model parameters and R applies the
  # parameterisation chain rule, vectorised over observations. That only works if
  # the array's dimensions are labelled, so it is asserted rather than assumed.
  pars <- c("vcmax_25", "stem_b", "cost_scale_TF24")
  b <- leaf_batch(psi_soil = c(1.0, 1.5, 2.0, 2.5), PPFD = 900)
  g <- leaf_gradient_batch(b, pars = pars)
  outs <- c("A", "gc", "psi_stem", "collar", "profit")
  # ⚠️ WRITTEN OUT HERE ON PURPOSE, and it is the one place in the package that
  # does. Everything else derives the list from `gradient_output_names()`, which is
  # what makes the two routes unable to disagree about it -- so nothing else would
  # notice a column silently appearing, disappearing or moving. A caller reads
  # these by position after applying a Jacobian, so the ORDER is as much of an
  # interface as the names: `profit` is appended, not inserted.
  expect_identical(gradient_output_names(), outs)
  expect_identical(dim(g$gradient), c(4L, 3L, 5L))
  expect_identical(dimnames(g$gradient)[[2]], pars)
  expect_identical(dimnames(g$gradient)[[3]], outs)
  expect_identical(dim(g$value), c(4L, 5L))
  expect_identical(colnames(g$value), outs)
  for (f in c("status", "method", "message")) {
    expect_length(g[[f]], 4L)
  }
  for (f in c("H", "stationarity")) {
    expect_length(g[[f]], 4L)
  }
  # `pars` selects, and its ORDER is the array's order -- a fit reads columns by
  # position after this, so a permuted request has to come back permuted rather
  # than in the enumeration's order.
  rev_g <- leaf_gradient_batch(b, pars = rev(pars))
  expect_identical(dimnames(rev_g$gradient)[[2]], rev(pars))
  expect_identical(g$gradient[, "vcmax_25", ], rev_g$gradient[, "vcmax_25", ])
})

test_that("`pars` order does not change any gradient (#72)", {
  # ⚠️ THE REGRESSION GUARD FOR #72, AND IT IS THE PRIMARY ONE -- not the golden
  # file. The defect was worth up to 3.4e-5 relative, and the golden file is
  # compared with a 1e-3 tolerance off macOS/arm64, so on Linux it could not
  # catch a recurrence. This test compares two runs IN THE SAME PROCESS, so it is
  # exact everywhere and needs no tolerance at all.
  #
  # What went wrong, because the shape of it is what this asserts: the `stem_b`
  # shortcut calls `perturb_stem_b()`, which rescales the spline and touches
  # nothing else, while the loops only restore base AFTER the whole parameter
  # loop. So a `stem_b` that was not first was differentiated at a point
  # displaced by one step in whichever parameter preceded it.
  #
  # ⚠️ Asserted as the GENERAL invariant -- order changes nothing, for any
  # parameter -- rather than as "stem_b is now reseated". `root_b` obeys the same
  # homogeneity identity and would get the same shortcut, and a stem_b-shaped
  # test would silently stop covering the case it was written for.
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900)
  pars <- c("vcmax_25", "stem_c", "a", "stem_b", "cost_scale_TF24", "root_b")

  for (method in c("auto", "fd")) {
    for (fast in c(TRUE, FALSE)) {
      lab <- paste(method, "fast =", fast)
      ref <- leaf_gradient_batch(b, pars = pars, method = method,
                                 fast_stem_curve = fast)$gradient[1, , ]
      # Reversed, and rotated, so the claim is not "one other ordering agrees".
      for (perm in list(rev(pars), pars[c(4, 1, 6, 2, 5, 3)],
                        c("stem_b", setdiff(pars, "stem_b")))) {
        got <- leaf_gradient_batch(b, pars = perm, method = method,
                                   fast_stem_curve = fast)$gradient[1, , ]
        expect_identical(got[pars, ], ref[pars, ],
                         label = paste(lab, paste(perm, collapse = ",")))
      }
    }
  }

  # ⚠️ And the predecessor that made the defect WORST, which is not the one you
  # would pick by eye: `a` = 0.30 and `curv_fact_colim` = 0.99 take the absolute
  # 1e-6 step floor rather than a relative step, and `stem_c`'s perturbation
  # REBUILDS the very curve the `stem_b` shortcut then rescales. Each of the
  # thirteen is checked as a predecessor rather than a sample of them.
  solo <- leaf_gradient_batch(b, pars = "stem_b")$gradient[1, "stem_b", ]
  for (p in setdiff(names(leaf_traits()), "stem_b")) {
    got <- leaf_gradient_batch(b, pars = c(p, "stem_b"))$gradient[1, "stem_b", ]
    expect_identical(got, solo, label = paste("after", p))
  }

  # The R reference has the same fix, so the two implementations still agree --
  # which is the point of fixing them together in one change rather than letting
  # them drift.
  r_ref <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                         pars = c("cost_scale_TF24", "stem_b"))
  batch <- leaf_gradient_batch(b, pars = c("cost_scale_TF24", "stem_b"))
  expect_identical(as.vector(batch$gradient[1, , ]), as.vector(r_ref$gradient))
})

test_that("leaf_batch() recycles its drivers the way leaf_solve() does", {
  # Same helpers, so the rules cannot drift -- but that is only true while both
  # actually call them, which is what this checks.
  b <- leaf_batch(psi_soil = c(1.0, 1.5, 2.0), PPFD = 900)
  expect_identical(b$n, 3L)
  g <- leaf_gradient_batch(b, pars = "vcmax_25")
  solved <- leaf_solve(psi_soil = c(1.0, 1.5, 2.0), PPFD = 900)
  # EVERY column, not two of them, and against a reader that names its fields
  # independently. `gradient::outputs()` copies five `Leaf` members into an
  # unnamed array by position while `operating_point_values()` copies the same
  # members into a twelve-long one; `gradient_output_names()` is the only thing
  # tying those two orders together, and it is a list of names with nothing in the
  # types to check it against the fields. This is that check.
  for (o in gradient_output_names()) {
    expect_identical(g$value[, o], solved[[o]], label = paste("value", o))
  }

  # A plain numeric psi_soil is N single-layer observations; a list is one
  # multi-layer observation each. Getting that backwards is the easiest mistake to
  # make with this function, and it fails loudly rather than silently, because the
  # network and the profile then disagree in length.
  expect_identical(leaf_batch(psi_soil = list(c(1, 1.5, 2)), PPFD = 900)$n, 1L)
  expect_length(leaf_gradient_batch(
    leaf_batch(psi_soil = list(c(1, 1.5, 2)), PPFD = 900),
    pars = "vcmax_25")$status, 1L)
})

test_that("leaf_gradient_batch() rejects what it cannot do", {
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900)
  s <- leaf_batch(psi_soil = 1.5, supply = leaf_supply_single())

  expect_error(leaf_gradient_batch(b, pars = "not_a_trait"),
               "cannot differentiate")
  # The same message `leaf_gradient()` gives, from the same function: `resistance`
  # is a real parameter on the other supply path, so saying "unknown" would be
  # misleading rather than merely unhelpful.
  expect_error(leaf_gradient_batch(b, pars = "resistance"),
               "single-potential path only")
  expect_true(is.finite(leaf_gradient_batch(s, pars = "resistance")$gradient[
    1, 1, "A"]))

  expect_error(leaf_gradient_batch(b, step = -1), "positive number")
  expect_error(leaf_gradient_batch(b, method = "magic"))
  expect_error(leaf_gradient_batch(leaf_traits()), "must come from leaf_batch")
  expect_error(leaf_batch(psi_soil = 1.5, traits = list(vcmax_25 = 96)),
               "leaf_traits")

  # `traits` and `theta` both say where the gradient is taken, so passing both is
  # refused rather than resolved in favour of one of them.
  th <- matrix(1, nrow = 1, ncol = 15)
  expect_error(leaf_gradient_batch(b, traits = leaf_traits(), theta = th),
               "pass one")
  expect_error(leaf_gradient_batch(b, theta = matrix(1, 1, 14)), "15 columns")
  expect_error(leaf_gradient_batch(b, theta = matrix(1, 3, 15)),
               "1 row or one per observation")
  expect_error(leaf_gradient_batch(b, theta = as.data.frame(th)),
               "numeric matrix")
  # ⚠️ Column names are checked against the enumeration ORDER, because C++
  # indexes by position: a reordered matrix would differentiate the wrong
  # parameters and return plausible numbers rather than fail.
  bad <- th
  colnames(bad) <- rev(gradient_par_names())
  expect_error(leaf_gradient_batch(b, theta = bad), "in that order")
})

test_that("a batch that has lost its C++ pointer says so rather than crashing", {
  # ⚠️ The prepared drivers and the leaf are both external pointers, so a batch
  # that has been through `serialize()` holds NULL ones -- and that is true within
  # the SAME session too, which is worth knowing: `readRDS()` does not restore a
  # pointer, it restores a null. Serialising a batch is a reasonable thing for
  # someone to try, and "external pointer is not valid" -- what the first read of
  # either would otherwise report, from whichever one Rcpp unwrapped first -- is
  # true and gives no hint that the fix is to rebuild.
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900)
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(b, f)
  expect_error(leaf_gradient_batch(readRDS(f), pars = "vcmax_25"),
               "does not survive")
  # The live one still works afterwards: the check is on the handle, not a flag
  # set by serialising, so `saveRDS()` does not poison the object it copied.
  expect_true(is.finite(leaf_gradient_batch(b, pars = "vcmax_25")$gradient[
    1, 1, "A"]))
  # Something that was never a handle at all, which the same check has to catch:
  # dereferencing a double as a pointer is not an error the session survives.
  expect_error(phylloptim:::gradient_batch_check(1.0), "not a prepared batch")
})
