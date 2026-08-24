# Guards against the R boundary quietly getting more expensive.
#
# WHY THIS FILE EXISTS. Nothing else in this package looks at the clock, on the
# stated grounds that absolute microseconds are machine-dependent. That is true, and
# it is also exactly how a 4.5x regression reached a green CI: `set_drivers()` began
# constructing a `RootNetwork` per call, which took `leaf_solve()` from 23.7 to
# 108 us/row, and the whole suite, the golden tie-back and six CI jobs passed.
#
# So the guards here are of two kinds, in order of trustworthiness:
#
#   1. COUNTING, which is deterministic. The regression was a `.Call` appearing in a
#      per-row path, and a count of `.Call`s has no variance at all. This is the real
#      guard; prefer adding to it.
#   2. A RATIO, for what counting cannot foresee. Never an absolute time: a bare
#      `.Call` on one machine moved 0.69 -> 1.10 us between two runs an hour apart on
#      the same build. Cost divided by the cost of a trivial `.Call` in the SAME
#      process is stable across machines, so that is what is bounded -- generously,
#      so it catches a doubling rather than policing a few percent.
#
# `tools/bench_user_cost.R` is the measurement these thresholds came from, and
# `tools/bench_history.sh` runs it against past commits.

# Count how many times a generated `.Call` wrapper is invoked while evaluating
# `expr`, by swapping the namespace binding for a counting one. Deterministic, and
# it observes the thing that actually went wrong rather than a symptom of it.
count_calls <- function(name, expr) {
  ns <- asNamespace("phylloptim")
  orig <- get(name, envir = ns)
  n <- 0L
  unlockBinding(name, ns)
  on.exit({
    assign(name, orig, envir = ns)
    lockBinding(name, ns)
  }, add = TRUE)
  assign(name, function(...) {
    n <<- n + 1L
    orig(...)
  }, envir = ns)
  force(expr)
  n
}

test_that("set_drivers() does not build a supply network per call", {
  # THE REGRESSION THIS FILE WAS WRITTEN FOR. Building the default network per call
  # is correct and 4.5x slower; the memo makes it once per distinct soil profile. So
  # the invariant is a COUNT, not a time: it must not scale with the number of rows.
  #
  # ⚠️ THE TWO PATHS REACH C++ THROUGH DIFFERENT SYMBOLS, and getting that wrong
  # makes this test pass while measuring nothing. The multi-layer default goes through
  # `root_network_from_carbon()`, a plain Rcpp export; the single-potential default
  # goes through `series_resistance()` and so through the RcppR6-generated
  # `RootNetwork__ctor`. Counting the constructor on the multi-layer path returns 0
  # whatever happens, because the constructor is not on that path at all.
  l <- leaf_model()
  set_drivers(l, psi_soil = 2.0)                      # warm the memo
  n <- count_calls("root_network_from_carbon",
                   for (i in 1:20) set_drivers(l, psi_soil = 2.0))
  expect_identical(n, 0L)

  # A genuinely new soil profile MUST build one -- otherwise the memo would be
  # returning a stale network for a different profile, which is the wrong-number
  # failure the memo test covers. One, not twenty.
  n <- count_calls("root_network_from_carbon",
                   for (i in 1:20) set_drivers(l, psi_soil = c(1, 2),
                                               soil_depth = c(0.5, 1.0)))
  expect_identical(n, 1L)

  # The single-potential path, through its own mechanism.
  s <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(s, psi_soil = 1.5)
  n <- count_calls("RootNetwork__ctor",
                   for (i in 1:20) set_drivers(s, psi_soil = 1.5))
  expect_identical(n, 0L)
})

test_that("series_resistance() does not reach C++ after the first call", {
  # A fitted resistance changes every proposal, so this cannot be memoised on its
  # value -- it copies a session-cached prototype instead (#67). The guard is that
  # varying the value does not reconstruct.
  series_resistance(1e3)                              # warm the prototype
  n <- count_calls("RootNetwork__ctor",
                   for (r in seq(1e3, 2e3, length.out = 20)) series_resistance(r))
  expect_identical(n, 0L)
})

test_that("leaf_solve() crosses the boundary a bounded number of times per row", {
  # Not a time: a count of solves. `reuse = TRUE` must build ONE Leaf for the whole
  # sweep, which is the other 155-us-per-call trap on this surface (issue #52 is the
  # same cost in leaf_gradient()).
  n1 <- count_calls("Leaf__ctor", leaf_solve(psi_soil = rep(1.5, 4)))
  n2 <- count_calls("Leaf__ctor", leaf_solve(psi_soil = rep(1.5, 32)))
  expect_identical(n1, n2)      # constant in the number of rows
  expect_lte(n2, 2L)

  # reuse = FALSE is documented to build one per row; if that stops being true the
  # test above is no longer testing anything.
  expect_gt(count_calls("Leaf__ctor",
                        leaf_solve(psi_soil = rep(1.5, 8), reuse = FALSE)), 4L)
})

test_that("leaf_gradient(x =) constructs no Leaf at all", {
  # #52. Construction is ~146 us, about 40% of a one-parameter gradient, and `x`
  # exists to remove it from a per-observation loop. The guard is a COUNT: reuse must
  # construct zero, however many gradients are taken, or the argument is decorative.
  tr <- leaf_traits()
  l <- leaf_model(traits = tr)
  leaf_gradient(psi_soil = 2.0, x = l, traits = tr, pars = "vcmax_25")   # warm

  n <- count_calls("Leaf__ctor", {
    for (i in 1:5) {
      leaf_gradient(psi_soil = 2.0, x = l, traits = tr, pars = "vcmax_25")
    }
  })
  expect_identical(n, 0L)

  # And without `x` it constructs exactly one per call -- if that stopped being true
  # the test above would be measuring nothing.
  n <- count_calls("Leaf__ctor", {
    for (i in 1:5) {
      leaf_gradient(psi_soil = 2.0, traits = tr, pars = "vcmax_25")
    }
  })
  expect_identical(n, 5L)
})

test_that("leaf_gradient_batch() crosses the boundary once, whatever N is", {
  # ⚠️ THIS IS THE WHOLE CLAIM OF #4 STAGE 2, AND IT IS COUNTABLE. A four-parameter
  # gradient through `leaf_gradient()` crosses the boundary 112 times per
  # observation and spends 1.5% of its time in the model. The batch composes the
  # composite in C++ instead, so the crossing count must be CONSTANT IN N -- one
  # call into `gradient_batch_run`, per call and not per row. A timing assertion
  # would not say this: it would drift with the machine, and it would still pass if
  # the count went back to being per-row on a fast day.
  b1 <- leaf_batch(psi_soil = rep(1.5, 4), PPFD = 900)
  b2 <- leaf_batch(psi_soil = rep(1.5, 64), PPFD = 900)
  pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale", "TF24_beta2")
  leaf_gradient_batch(b1, pars = pars)                       # warm anything lazy

  n1 <- count_calls("gradient_batch_run", leaf_gradient_batch(b1, pars = pars))
  n2 <- count_calls("gradient_batch_run", leaf_gradient_batch(b2, pars = pars))
  expect_identical(n1, 1L)
  expect_identical(n2, 1L)                      # 16x the rows, the same one call

  # ⚠️ AND NO MODEL WORK LEAKED BACK TO R, which is the stronger statement and the
  # one that would catch a future "small" convenience added to the R wrapper. Every
  # per-perturbation primitive `leaf_gradient()` reaches through -- the setter, the
  # solve, the two gradient evaluations -- must be called ZERO times by the batch,
  # for any N. A single one of these appearing in the R path is the regression this
  # file exists for, at 11 crossings per parameter per observation.
  for (sym in c("Leaf__set_traits", "Leaf__set_physiology",
                "Leaf__find_root_collar_psi", "Leaf__evaluate_root_collar_psi",
                "Leaf__dprofit_droot_collar_psi", "Leaf__operating_point_values",
                "Leaf__perturb_stem_P50", "Leaf__ctor", "RootNetwork__ctor",
                "root_network_from_carbon")) {
    expect_identical(count_calls(sym, leaf_gradient_batch(b2, pars = pars)), 0L,
                     label = sym)
  }

  # `leaf_batch()` is where the per-observation R work is allowed to live, because
  # it runs once per fit rather than once per likelihood evaluation -- but it must
  # still build exactly ONE Leaf, not one per observation (#52's trap on a new
  # surface), and reach the drivers' C++ side exactly once.
  expect_identical(count_calls("Leaf__ctor",
                               leaf_batch(psi_soil = rep(1.5, 64))), 1L)
  expect_identical(count_calls("gradient_batch_prepare",
                               leaf_batch(psi_soil = rep(1.5, 64))), 1L)
})

test_that("a driven row costs a bounded multiple of a trivial .Call", {
  # The backstop for what counting cannot foresee. Bound is ~2x the measured ratio,
  # so it catches a doubling and ignores machine noise. Measured on macOS/arm64,
  # 2026-08-05, at phylloptim 0.2.0: 14 x .Call for a hand-driven row and 20 x for a
  # leaf_solve() row. See tools/bench_user_cost.R.
  skip_on_cran()
  skip_on_ci()   # shared runners are too noisy for even a 2x bound to be honest

  timeit <- function(f, n) {
    f()
    min(replicate(5, {
      t0 <- Sys.time()
      for (i in seq_len(n)) f()
      1e6 * as.numeric(difftime(Sys.time(), t0, units = "secs")) / n
    }))
  }
  l <- leaf_model(supply = leaf_supply_singlelayer())
  net <- series_resistance(1e3)
  set_drivers(l, psi_soil = 1.5, root_network = net)
  l$find_root_collar_psi()

  ref <- timeit(function() l$ci_, 20000)
  row <- timeit(function() {
    set_drivers(l, psi_soil = 1.5, PPFD = 900, atm_vpd = 1.5, root_network = net)
    l$find_root_collar_psi()
    c(l$assim_colimited_, l$stom_cond_CO2_, l$opt_psi_stem_)
  }, 2000)

  expect_lt(row / ref, 30)      # measured 14
  expect_lt(timeit(function() leaf_solve(psi_soil = rep(1.5, 16),
                                         supply = leaf_supply_singlelayer(),
                                         root_network = net), 30) / 16 / ref, 45)  # measured 20
})
