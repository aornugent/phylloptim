# The README's R examples, actually run.
#
# ⚠️ WHY THIS FILE EXISTS. Vignettes are knitted, so a stale example in one fails
# a build. `README.md` is plain markdown -- nothing renders it, nothing parses it,
# and its code is prose as far as every tool in this repo is concerned. So when the
# `stem_b` trait became `stem_P50`, one README example was updated and two were
# not, and both would have errored for the next reader. Neither the C++ suite nor
# `R CMD check` nor CI could see it.
#
# The alternative was a knitted `README.Rmd`, which is the more usual answer and
# was not taken: the README has C++, CMake and Python blocks that no R engine can
# run, so knitting it would mean marking most of the file `eval = FALSE` and the
# guard would cover less than this does.
#
# ⚠️ THIS IS A SYNTAX-AND-API GUARD, NOT A CORRECTNESS ONE. It asserts that every
# ```r block parses and evaluates without error, which is what rots: a renamed
# argument, a removed method, a parameter no longer in `gradient_par_names()`. It
# says nothing about whether the NUMBERS quoted in the surrounding prose are still
# right -- those have to be re-read after a solver change, and this cannot help.
test_that("every r block in the README parses and runs", {
  skip_on_cran()
  # Present in the source tree and in the tarball; absent from the unpacked
  # check directory's test working directory, where this simply skips.
  readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md not reachable from here")

  lines <- readLines(readme, warn = FALSE)
  # Fenced blocks tagged exactly `r`, so the cpp/cmake/python/sh ones are left
  # alone. Openers and closers alternate, which is what makes the pairing safe.
  fences <- grep("^```", lines)
  opens <- fences[seq(1, length(fences), by = 2)]
  closes <- fences[seq(2, length(fences), by = 2)]
  expect_length(closes, length(opens))
  keep <- lines[opens] == "```r"
  expect_gt(sum(keep), 5L)   # the guard is worthless if it matched nothing

  # The names the examples read without defining: the README is written for
  # someone who has their own observations, so it says `obs$psi_soil` without
  # building an `obs`. Supplying them here is not weakening the test -- what it
  # checks is the package's surface, not the prose's self-containedness.
  env <- new.env(parent = globalenv())
  env$obs <- data.frame(psi_soil = c(1, 2, 3), PPFD = c(500, 900, 1200))
  env$traits <- leaf_traits()
  env$control <- leaf_control()

  for (k in which(keep)) {
    code <- lines[(opens[[k]] + 1L):(closes[[k]] - 1L)]
    label <- paste0("README.md:", opens[[k]] + 1L)
    exprs <- tryCatch(parse(text = code), error = function(e) e)
    expect_false(inherits(exprs, "error"),
                 info = paste(label, "does not parse:",
                              conditionMessage(exprs)))
    if (inherits(exprs, "error")) next
    err <- tryCatch({ for (e in exprs) eval(e, envir = env); NULL },
                    error = function(e) e)
    expect_null(err, info = if (is.null(err)) label else
                  paste(label, "failed:", conditionMessage(err)))
  }
})
