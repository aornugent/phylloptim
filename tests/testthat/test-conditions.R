# The infeasibility classification (#57).
#
# ⚠️ THE POINT OF TESTING THIS BY TOKEN RATHER THAN BY MESSAGE. The issue exists
# because callers had to match on prose, and the prose moves: #65 and #79 both
# rewrote the out-of-domain text and #92 added another variant. So these assert the
# `code`, never the sentence, and the drift test below is what stops the two ends of
# the mechanism parting company.

test_that("the codes in the headers are exactly the codes documented in R", {
  # ⚠️ THIS IS THE TEST THAT MATTERS MOST, and it is not about any one code. The
  # mechanism has two ends -- a token written in C++, a list documented in R -- and
  # nothing in either language connects them. A `stop_infeasible("something_new")`
  # added to a header would otherwise become an undocumented part of the public API,
  # discovered by a caller from a failure, which is the same class of rot as #64.
  #
  # It reads the INSTALLED headers, so it is also checking that the header carrying a
  # classification actually shipped.
  inc <- system.file("include", "phylloptim", package = "phylloptim")
  skip_if(!nzchar(inc) || !dir.exists(inc), "installed headers not found")
  src <- unlist(lapply(list.files(inc, pattern = "\\.hpp$", full.names = TRUE),
                       readLines, warn = FALSE))
  # ⚠️ COLLAPSED TO ONE STRING FIRST, and this is not tidiness. A line-by-line scan
  # misses `stop_infeasible(` whose code argument wrapped onto the next line, which
  # is what several of them do -- and it misses it SILENTLY, reporting the code as
  # absent from the headers. The first version of this test did exactly that.
  flat <- paste(src, collapse = " ")
  hits <- regmatches(flat, gregexpr('stop_infeasible\\(\\s*"[a-z_]+"', flat))[[1L]]
  codes <- sort(unique(gsub('.*"([a-z_]+)".*', "\\1", hits)))
  expect_gt(length(codes), 0L)
  expect_identical(codes, sort(names(leaf_infeasible_codes())))
})

test_that("every documented code has a description, and none is a placeholder", {
  codes <- leaf_infeasible_codes()
  expect_true(all(nzchar(names(codes))))
  expect_true(all(nchar(codes) > 30L))
  expect_false(anyDuplicated(names(codes)) > 0L)
})

test_that("an untagged error stays an ordinary error", {
  # The safe default, and the half of the design that protects #39's warning: a
  # mismatched driver vector must NOT be catchable as infeasibility, or the tryCatch
  # this mechanism exists to enable would swallow it and the fit would report a
  # plausible likelihood over the rows that survived.
  l <- leaf_model()
  err <- tryCatch(
    with_phylloptim_conditions(
      set_drivers(l, psi_soil = c(1, 2), soil_depth = 1.0, PPFD = 900)),
    condition = function(e) e)
  expect_s3_class(err, "error")
  expect_false(inherits(err, "phylloptim_infeasible"))
  expect_false(inherits(err, "phylloptim_error"))

  # And a plain R error passing through is untouched.
  expect_identical(
    tryCatch(with_phylloptim_conditions(stop("not ours")),
             error = conditionMessage),
    "not ours")
})

test_that("an infeasible exit is catchable by class and carries its code", {
  # A REACHABLE one, end to end, through the ordinary public entry point. The driver
  # is extreme but well formed and in range: with the stem conductance this small the
  # transpiration the solve needs falls outside the inverse spline's domain, so the
  # operating point genuinely cannot be evaluated. That is infeasibility rather than
  # a caller error, which is exactly the distinction being tested.
  caught <- tryCatch(
    leaf_solve(psi_soil = 2.0, PPFD = 900,
               leaf_specific_conductance_max = 1e-30),
    condition = function(e) e)

  expect_s3_class(caught, "phylloptim_infeasible")
  expect_s3_class(caught, "phylloptim_error")
  expect_s3_class(caught, "error")
  expect_identical(caught$code, "stem_curve_domain")
  expect_true(caught$code %in% names(leaf_infeasible_codes()))

  # ⚠️ leaf_solve() classifies WITHOUT the caller wrapping anything -- it routes
  # through the wrapper itself. That is the property a fit depends on.
  expect_identical(
    tryCatch(leaf_solve(psi_soil = 2.0, PPFD = 900,
                        leaf_specific_conductance_max = 1e-30),
             phylloptim_infeasible = function(e) e$code),
    "stem_curve_domain")

  # The token stays in the message on purpose: a log line should still say which
  # code it was, and the parsed `code` is what programs read.
  expect_match(conditionMessage(caught), "phylloptim:infeasible:stem_curve_domain",
               fixed = TRUE)
})

test_that("the same failure driven by hand needs the wrapper, and gets it", {
  # `$find_root_collar_psi()` is generated glue, so it cannot classify on its own --
  # this is the case with_phylloptim_conditions() exists for, and the reason it is
  # exported rather than internal.
  l <- leaf_model()
  set_drivers(l, psi_soil = 2.0, PPFD = 900,
              leaf_specific_conductance_max = 1e-30)

  bare <- tryCatch(l$find_root_collar_psi(), condition = function(e) e)
  expect_s3_class(bare, "error")
  expect_false(inherits(bare, "phylloptim_infeasible"))

  wrapped <- tryCatch(with_phylloptim_conditions(l$find_root_collar_psi()),
                      condition = function(e) e)
  expect_s3_class(wrapped, "phylloptim_infeasible")
  expect_identical(wrapped$code, "stem_curve_domain")
})

test_that("catching phylloptim_error catches every classified failure", {
  # The middle class exists so a caller does not have to enumerate codes, and so a
  # second category later does not rewrite every handler. Asserted on a synthetic
  # condition, because it is a property of the class vector rather than of any
  # particular exit -- and building it from the package's own constructor is what
  # keeps this from testing a hand-written class vector that has drifted.
  cnd <- phylloptim:::.infeasible_condition(
    simpleError("[phylloptim:infeasible:collar_solve] synthetic"), "collar_solve")
  expect_s3_class(cnd, "phylloptim_infeasible")
  expect_s3_class(cnd, "phylloptim_error")
  expect_s3_class(cnd, "error")
  expect_identical(tryCatch(stop(cnd), phylloptim_error = function(e) e$code),
                   "collar_solve")
})

test_that("the token parser is anchored and does not match prose", {
  # A message that merely MENTIONS the token shape mid-sentence is not a
  # classification. The pattern is anchored for that reason: a diagnostic that
  # quotes another error -- which several of these do, they embed `e.what()` -- must
  # not be re-classified by the quote.
  expect_identical(
    phylloptim:::.infeasible_code("[phylloptim:infeasible:uptake] no flux"),
    "uptake")
  expect_true(is.na(
    phylloptim:::.infeasible_code("failed: [phylloptim:infeasible:uptake] no flux")))
  expect_true(is.na(phylloptim:::.infeasible_code("ordinary failure")))
  expect_true(is.na(phylloptim:::.infeasible_code("")))
})

test_that("a solved leaf is unaffected by the wrapper", {
  # The wrapper is one tryCatch per CALL. It must not change any number, and it must
  # not change what a successful call returns.
  a <- leaf_solve(psi_soil = c(1.0, 2.0, 3.0), PPFD = 900)
  b <- with_phylloptim_conditions(
    leaf_solve(psi_soil = c(1.0, 2.0, 3.0), PPFD = 900))
  expect_identical(a, b)
})
