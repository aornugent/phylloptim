#!/usr/bin/env Rscript
# Regenerate tests/testthat/gradient_golden.tsv, the recorded trait gradients.
#
#   Rscript tools/gradient_golden.R > tests/testthat/gradient_golden.tsv
#
# WHY THIS FILE EXISTS, and it is not "because the other golden file does".
# `test-gradient-batch.R` compares the C++ composite against `leaf_gradient()`
# bit-for-bit, and that is the strong test of the transcription -- the two share
# no code. What it cannot see is a change applied CONSISTENTLY TO BOTH: move the
# step rule, the active-set threshold or a solver tolerance and the two
# implementations move together, the equality test passes, and nothing anywhere
# says the gradients changed. That is the trap #70 hit inside R alone, and it
# would be far easier to hit across a language boundary. So the values are also
# pinned to a file.
#
# ⚠️ HEX FLOATS, and do not "tidy" them into decimals. R's decimal
# string-to-double conversion is not correctly rounded: on a random sample about
# 18% of full-precision inputs come back one ULP away. R reads C99 hex exactly.
# That is issue #13, where a supposed 1-ULP disagreement with plant turned out to
# be entirely the measuring instrument.
#
# ⚠️ Only regenerate DELIBERATELY, and say in the commit message what moved and by
# how much. Running this after an accidental change rubber-stamps it. The guide to
# magnitudes is in .claude/CLAUDE.md: ~1e-16 is reassociation, ~1e-9 is the
# solver floor, ~1e-4 is a real difference.
#
# ⚠️ RUN THIS ON macOS/arm64 AND NOWHERE ELSE. Like tests/cpp/golden/, the file is
# bit-exact only on the platform that generated it -- libm's exp/pow are not
# bit-reproducible between glibc on x86-64 and Apple's libm on arm64 -- and
# test-gradient-batch.R compares with `golden_tolerance()` elsewhere. Regenerating
# on Linux would just move the failure to the platform the file came from.
#
# These are derivatives of outputs evaluated at the ARGMAX of a flat maximum, so
# they inherit that file's sqrt-amplified class -- and, being finite differences,
# one amplification more. Measured worst cross-platform disagreement 1.3e-3,
# against 1.4e-4 for the solved outputs themselves. A central difference cancels
# the systematic part of a libm difference, not all of it, and what is left is
# divided by the step: the smallest-magnitude parameter in the file therefore sets
# the tolerance for all of it, which since #41 is `R_d_25` at 1.44 rather than
# `vcmax_25` at 96. `gradient_golden_tolerance()` in tests/testthat/helper-golden.R
# carries the arithmetic.
suppressMessages(library(phylloptim))

# ⚠️ AND SAY WHICH BUILD THAT WAS. `library(phylloptim)` resolves against the
# library path, NOT against the tree you are standing in, so running this from a
# working tree with uninstalled changes regenerates the file from the SITE build --
# reproducing the values you were trying to replace, exactly, with no warning. That
# is not hypothetical: it happened while regenerating for #92 and reported "0 cells
# moved" against a suite that was failing on 63 of them. It is the same trap
# `tools/bench_history.sh` guards with --expect-version, and the same one the guide
# records for plant under "Validating against plant".
#
# So the path goes to stderr on every run, where it does not corrupt the TSV on
# stdout. Check it before believing the diff:
#
#   R CMD INSTALL --library=/tmp/lib . && \
#     R_LIBS=/tmp/lib Rscript tools/gradient_golden.R > tests/testthat/gradient_golden.tsv
# The check compares the INSTALLED HEADERS against this tree's, which is the actual
# question. `R CMD INSTALL` copies `inst/include/` into the package directory, so
# `system.file("include", ...)` is a verbatim copy of the headers the loaded build
# was compiled from -- comparing paths cannot work (an installed package is never
# in the working tree) and comparing versions cannot either, since the version is
# rarely bumped for a numerics change.
message("phylloptim from: ", find.package("phylloptim"),
        " (version ", as.character(utils::packageVersion("phylloptim")), ")")
local({
  tree <- "inst/include"
  if (!dir.exists(tree)) return(invisible(NULL))  # not run from the package root
  installed <- system.file("include", package = "phylloptim")
  digest_dir <- function(d) {
    f <- sort(list.files(d, recursive = TRUE, pattern = "[.]hpp?$"))
    vapply(f, function(p) paste(readLines(file.path(d, p), warn = FALSE),
                                collapse = "\n"), character(1))
  }
  if (!nzchar(installed) || !identical(digest_dir(tree), digest_dir(installed))) {
    stop("the loaded phylloptim was built from DIFFERENT headers than ", tree,
         ".\nRegenerating now would reproduce the installed build's values, not ",
         "this tree's.\nInstall first:\n",
         "  R CMD INSTALL --library=/tmp/lib . && \\\n",
         "    R_LIBS=/tmp/lib Rscript tools/gradient_golden.R > ",
         "tests/testthat/gradient_golden.tsv",
         call. = FALSE)
  }
})

# The golden grid's drivers, so every row here is an operating point the rest of
# the suite already knows -- tests/cpp/test_golden.cpp, test-golden.R and
# test-gradient.R all use these.
grid_drivers <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L) {
  theta <- 0.000157
  area_leaf <- 0.05
  list(psi_soil = list(psi_soil + 0.25 * (seq_len(layers) - 1)),
       PPFD = ppfd,
       soil_depth = list(1.0 * seq_len(layers)),
       root_network = root_network_from_carbon(
         rep(1 / layers / area_leaf, layers),
         soil_depth = 1.0 * seq_len(layers)),
       leaf_specific_conductance_max = 1.0 * theta / 5.0,
       atm_vpd = vpd, ca = 40, leaf_temp = 25, atm_o2_kpa = 21,
       atm_kpa = 101.3)
}

# Five rows, chosen for the branch each one takes rather than to sample evenly.
# `stem_P50` and `stem_c` are in `pars` throughout because they are the two
# parameters of one curve and their derivatives are taken by different routes:
# `stem_P50` moves the curve by the homogeneity rescale, `stem_c` has no such
# identity and rebuilds. Both now carry a term through the derived `stem_b` and
# `psi_crit`, which is the part a hand-written chain rule would be most likely to
# drop -- so pinning them here is what would catch a partial masquerading as a
# total derivative.
#
# `R_d_25` is in `pars` throughout because it is the smallest-magnitude parameter
# here, so it takes the smallest absolute step and sets this file's cross-platform
# tolerance -- see `gradient_golden_tolerance()`.
cases <- list(
  list(label = "interior-1layer",
       args = grid_drivers(2.0),
       pars = c("vcmax_25", "stem_P50", "stem_c", "R_d_25")),
  list(label = "interior-5layer",
       args = grid_drivers(0.5, vpd = 0.5, layers = 5L),
       pars = c("vcmax_25", "stem_P50", "stem_c", "R_d_25")),
  list(label = "pinned-dry-3layer",
       args = grid_drivers(4.0, vpd = 0.5, layers = 3L),
       pars = c("vcmax_25", "stem_P50", "stem_c", "R_d_25")),
  list(label = "shutdown-1layer",
       args = grid_drivers(6.0),
       pars = c("vcmax_25", "stem_P50", "stem_c", "R_d_25")),
  list(label = "single-potential",
       args = list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
                   supply = leaf_supply_singlelayer(),
                   root_network = series_resistance(1e4)),
       pars = c("vcmax_25", "leaf_specific_conductance_max", "resistance",
                "R_d_25"))
)

out <- do.call(rbind, lapply(cases, function(cs) {
  b <- do.call(leaf_batch, cs$args)
  g <- leaf_gradient_batch(b, pars = cs$pars)
  do.call(rbind, lapply(cs$pars, function(p) {
    data.frame(case = cs$label, status = g$status[[1]], method = g$method[[1]],
               par = p,
               A = sprintf("%a", g$gradient[1, p, "A"]),
               gc = sprintf("%a", g$gradient[1, p, "gc"]),
               psi_stem = sprintf("%a", g$gradient[1, p, "psi_stem"]),
               collar = sprintf("%a", g$gradient[1, p, "collar"]),
               profit = sprintf("%a", g$gradient[1, p, "profit"]),
               stringsAsFactors = FALSE)
  }))
}))

write.table(out, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
