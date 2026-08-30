#!/usr/bin/env Rscript
# What a user pays, on the R side, for the three things a calibration actually does:
# one solve, N solves, and a gradient. The C++ side of the same three lives in
# tests/cpp/bench_solve.cpp; tools/bench_history.sh
# runs both and puts them in one table.
#
# WHY THIS EXISTS. The package asserts no absolute timings anywhere else, on the
# grounds that they are machine-dependent -- and that is exactly how a 4.5x
# regression in `set_drivers()` reached a green CI. Every test passed, `R CMD check`
# passed, and nothing looked at the clock.
#
# ⚠️ THE RECORDED QUANTITY IS A RATIO, NOT MICROSECONDS. The same build measured an
# hour apart on this machine moved a bare `.Call` from 0.69 to 1.10 us. What is
# stable across machines and runs is a cost DIVIDED BY the cost of a trivial `.Call`
# in the same process, so that is what is compared. Absolute figures are printed for
# orientation; never pin one.
#
# VERSION-AWARE on purpose -- it is meant to run against older commits in a
# worktree, and the driver API changed twice (#63 root carbon -> resistances, #66 the
# single path's resistance -> a driver). It detects the surface rather than assuming.
suppressMessages(library(phylloptim))

args <- commandArgs(TRUE)
tsv  <- "--tsv" %in% args
reps <- if ("--reps" %in% args) as.integer(args[[which(args == "--reps") + 1L]]) else 7L

# ⚠️ REFUSE TO MEASURE THE WRONG PACKAGE. `--expect-version` is not paranoia: this
# script is run against old commits with `R_LIBS` pointed at a throwaway library, and
# if that install did not produce a package called `phylloptim` then `library()` falls
# back to the SITE build and every number below describes the current tree while
# claiming to describe history. That happened: commits before #47 declare
# `Package: leaf`, installed cleanly as `leaf`, and reported the current API and
# version with no indication anything was wrong. A silent fallback is the failure
# mode this whole file exists to catch, so it is an error here.
expect <- if ("--expect-version" %in% args) {
  args[[which(args == "--expect-version") + 1L]]
} else {
  NA_character_
}
got <- as.character(utils::packageVersion("phylloptim"))
if (!is.na(expect) && !identical(expect, got)) {
  stop(sprintf(paste("loaded phylloptim %s but expected %s -- the intended library",
                     "was not used (a pre-#47 commit installs as `leaf`, not",
                     "`phylloptim`). Refusing to report numbers for the wrong build."),
               got, expect), call. = FALSE)
}
cat(sprintf("# phylloptim %s from %s\n", got, find.package("phylloptim")), file = stderr())

timeit <- function(f, n, reps) {
  f()
  min(replicate(reps, {
    t0 <- Sys.time()
    for (i in seq_len(n)) f()
    1e6 * as.numeric(difftime(Sys.time(), t0, units = "secs")) / n
  }))
}

# --- which API surface is this? -------------------------------------------
sd_args <- names(formals(set_drivers))
has_network <- "root_network" %in% sd_args
single_takes_resistance <- "resistance" %in% names(formals(leaf_supply_single))
api <- paste0(if (has_network) "resistance-driver" else "root-carbon",
              if (single_takes_resistance) "/supply-cfg" else "/supply-driver")

supply <- if (single_takes_resistance) {
  leaf_supply_single(resistance = 1e3)
} else {
  leaf_supply_single()
}
net <- if (has_network && !single_takes_resistance) series_resistance(1e3) else NULL

# Extra driver arguments this version wants, built ONCE outside any loop -- which is
# what a caller who has read the docs does, and what makes the comparison fair.
extra <- if (is.null(net)) list() else list(root_network = net)

l <- leaf_model(supply = supply)
drive <- function() {
  do.call(set_drivers, c(list(l, psi_soil = 1.5, PPFD = 900, atm_vpd = 1.5, ca = 40,
                              leaf_temp = 25, atm_o2_kpa = 21, atm_kpa = 101.3,
                              leaf_specific_conductance_max = 3.14e-5), extra))
}
read3 <- function() c(l$assim_colimited_, l$stom_cond_CO2_, l$opt_psi_stem_)
drive(); l$find_root_collar_psi()

# --- the boundary unit ----------------------------------------------------
ref <- timeit(function() l$ci_, 20000, reps)

# --- 1. ONE solve, driven by hand (the calibration inner loop) -------------
one <- timeit(function() { drive(); l$find_root_collar_psi(); read3() }, 3000, reps)
dr  <- timeit(drive, 5000, reps)
sv  <- timeit(function() l$find_root_collar_psi(), 3000, reps)
rd  <- timeit(read3, 20000, reps)

# --- 2. N solves through the vectorised wrapper (amortised per row) --------
N <- 32L
solve_args <- c(list(psi_soil = rep(1.5, N), PPFD = 900, atm_vpd = 1.5,
                     supply = supply), extra)
many <- timeit(function() do.call(leaf_solve, solve_args), 40, reps) / N

# --- 3. constructing a Leaf, which is 3's dominant term --------------------
ctor <- timeit(function() leaf_model(supply = supply), 300, reps)

out <- data.frame(
  version = as.character(utils::packageVersion("phylloptim")),
  api = api,
  ref_us = round(ref, 3),
  one_us = round(one, 2),      one_ratio   = round(one / ref, 1),
  drive_us = round(dr, 2),     drive_ratio = round(dr / ref, 1),
  solve_us = round(sv, 2),
  read3_us = round(rd, 2),
  many_us = round(many, 2),    many_ratio  = round(many / ref, 1),
  ctor_us = round(ctor, 1)
)
if (tsv) {
  write.table(out, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  cat("\nphylloptim ", out$version, "   api: ", out$api, "\n", sep = "")
  cat(sprintf("  trivial .Call (the unit)          %8.3f us\n", out$ref_us))
  cat("\n  R SIDE                              us      x.Call\n")
  cat(sprintf("    1 solve, hand-driven          %8.2f  %8.1f   <-- calibration inner loop\n", out$one_us, out$one_ratio))
  cat(sprintf("      of which set_drivers()      %8.2f  %8.1f\n", out$drive_us, out$drive_ratio))
  cat(sprintf("      of which the solve          %8.2f\n", out$solve_us))
  cat(sprintf("      of which 3 field reads      %8.2f\n", out$read3_us))
  cat(sprintf("    %d solves, leaf_solve(), /row  %8.2f  %8.1f\n", N, out$many_us, out$many_ratio))
  cat(sprintf("    leaf_model() construction     %8.1f\n", out$ctor_us))
}
