#!/usr/bin/env Rscript
# What a user pays, on the R side, for the three things a calibration actually does:
# one solve, N solves, and a gradient. The C++ side of the same three lives in
# tests/cpp/bench_solve.cpp and tests/cpp/bench_gradient.cpp; tools/bench_history.sh
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
single_takes_resistance <- "resistance" %in% names(formals(leaf_supply_singlelayer))
api <- paste0(if (has_network) "resistance-driver" else "root-carbon",
              if (single_takes_resistance) "/supply-cfg" else "/supply-driver")

supply <- if (single_takes_resistance) {
  leaf_supply_singlelayer(resistance = 1e3)
} else {
  leaf_supply_singlelayer()
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

# --- 4. N gradients x P pars: THE CALIBRATION SHAPE -----------------------
# ⚠️ THIS IS THE NUMBER A USER FEELS, and it is neither of the single-call figures
# above. A calibration takes a gradient PER OBSERVATION over SEVERAL fitted
# parameters, so both counts matter and they pull opposite ways: construction is a
# fixed cost per call, while the differentiated work is linear in P. Reuse therefore
# saves a large fraction at P = 1 and a much smaller one at P = 4, and quoting the
# P = 1 figure overstates what a fit gets.
NOBS <- 24L
PSET <- list(p1 = "vcmax_25",
             p4 = c("vcmax_25", "jmax_25", "TF24_cost_scale", "TF24_beta2"))
gN <- list(p1 = c(fresh = NA_real_, reuse = NA_real_),
           p4 = c(fresh = NA_real_, reuse = NA_real_))
if (exists("leaf_gradient")) {
  psv <- seq(0.5, 3.0, length.out = NOBS)
  # Built here rather than reused from section 5, which runs after this one.
  # `supply` is dropped for the reuse arm: it is fixed when the leaf is built and
  # passing both is refused rather than resolved.
  loop_args <- c(list(PPFD = 900, atm_vpd = 1.5, supply = supply), extra)
  loop_reuse <- c(list(PPFD = 900, atm_vpd = 1.5), extra)
  can_reuse <- "x" %in% names(formals(leaf_gradient))
  ll <- if (can_reuse) leaf_model(traits = leaf_traits(), supply = supply) else NULL
  for (nm in names(PSET)) {
    fit <- PSET[[nm]]
    gN[[nm]]["fresh"] <- timeit(function() {
      for (p in psv) {
        do.call(leaf_gradient, c(loop_args, list(psi_soil = p, pars = fit)))
      }
    }, 3, min(reps, 3)) / NOBS
    if (can_reuse) {
      gN[[nm]]["reuse"] <- timeit(function() {
        for (p in psv) {
          do.call(leaf_gradient, c(loop_reuse, list(psi_soil = p, pars = fit,
                                                    x = ll,
                                                    traits = leaf_traits())))
        }
      }, 3, min(reps, 3)) / NOBS
    }
  }
}

# --- 4b. THE SAME SHAPE, COMPOSED IN C++ (#4 stage 2) ----------------------
# ⚠️ THE ONLY FAIR COMPARISON IS AGAINST SECTION 4, not against the single-call
# figures: `leaf_gradient_batch()` amortises its own per-call cost over the
# observations, so quoting it at N = 1 understates it and quoting the R loop at
# N = 1 overstates that. Both arms here run the same NOBS observations over the
# same DIFFERENTIATED parameters (`length(pars)`), and the batch is prepared
# OUTSIDE the timed region --
# which is what a fit does, since the drivers do not change across draws.
gbatch <- c(p1 = NA_real_, p4 = NA_real_)
gprep <- NA_real_
if (exists("leaf_gradient_batch")) {
  batch_args <- c(list(psi_soil = psv, PPFD = 900, atm_vpd = 1.5,
                       supply = supply), extra)
  bb <- do.call(leaf_batch, batch_args)
  for (nm in names(PSET)) {
    fit <- PSET[[nm]]
    gbatch[[nm]] <- timeit(function() leaf_gradient_batch(bb, pars = fit),
                           30, reps) / NOBS
  }
  # Reported separately rather than folded in, because it is paid ONCE per fit
  # against 30,000-odd draws -- and because it is the term that would grow if
  # someone moved driver resolution back inside the hot path.
  gprep <- timeit(function() do.call(leaf_batch, batch_args), 5, reps)
}

# --- 5. a gradient, per call and per parameter -----------------------------
# NA on commits before #42, which is where leaf_gradient() was added. The whole
# point of running this against history is that the surface changes, so absence is
# recorded rather than fatal.
g1 <- g3 <- greuse <- NA_real_
if (exists("leaf_gradient")) {
  gp <- c("vcmax_25", "stem_P50", "TF24_cost_scale")
  grad_args <- c(list(psi_soil = 1.5, PPFD = 900, atm_vpd = 1.5, supply = supply,
                      pars = gp), extra)
  g3 <- timeit(function() do.call(leaf_gradient, grad_args), 30, reps)
  g1 <- timeit(function() do.call(leaf_gradient,
                c(grad_args[setdiff(names(grad_args), "pars")],
                  list(pars = gp[1]))), 30, reps)
  # The reuse path (#52), where this version has one. NA before it, so the column
  # shows when it arrived rather than pretending it was always there.
  if ("x" %in% names(formals(leaf_gradient))) {
    gl <- leaf_model(traits = leaf_traits(), supply = supply)
    # `supply` must come OUT: it is fixed when `gl` is built, and passing both is
    # rejected rather than resolved. `traits` must go IN, since a Leaf does not
    # expose its own.
    reuse_args <- c(grad_args[setdiff(names(grad_args), c("pars", "supply"))],
                    list(pars = gp[1], x = gl, traits = leaf_traits()))
    greuse <- timeit(function() do.call(leaf_gradient, reuse_args), 30, reps)
  }
}

out <- data.frame(
  version = as.character(utils::packageVersion("phylloptim")),
  api = api,
  ref_us = round(ref, 3),
  one_us = round(one, 2),      one_ratio   = round(one / ref, 1),
  drive_us = round(dr, 2),     drive_ratio = round(dr / ref, 1),
  solve_us = round(sv, 2),
  read3_us = round(rd, 2),
  many_us = round(many, 2),    many_ratio  = round(many / ref, 1),
  ctor_us = round(ctor, 1),
  grad1_us = round(g1, 1),     grad1_ratio = round(g1 / ref, 1),
  grad3_us = round(g3, 1),     grad3_ratio = round(g3 / ref, 1),
  grad1_reuse_us = round(greuse, 1),
  grad1_reuse_ratio = round(greuse / ref, 1),
  gradN_p1_fresh_us = round(gN$p1[["fresh"]], 1),
  gradN_p1_reuse_us = round(gN$p1[["reuse"]], 1),
  gradN_p4_fresh_us = round(gN$p4[["fresh"]], 1),
  gradN_p4_reuse_us = round(gN$p4[["reuse"]], 1),
  gradN_p1_batch_us = round(gbatch[["p1"]], 2),
  gradN_p4_batch_us = round(gbatch[["p4"]], 2),
  gradN_p4_batch_ratio = round(gbatch[["p4"]] / ref, 1),
  batch_prepare_us = round(gprep, 0),
  grad_per_par_us = round((g3 - g1) / 2, 1)
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
  cat(sprintf("    gradient, 1 par               %8.1f  %8.1f\n", out$grad1_us, out$grad1_ratio))
  cat(sprintf("    gradient, 3 pars              %8.1f  %8.1f\n", out$grad3_us, out$grad3_ratio))
  if (!is.na(out$grad1_reuse_us))
    cat(sprintf("    gradient, 1 par, x = a Leaf   %8.1f  %8.1f   <-- #52\n",
                out$grad1_reuse_us, out$grad1_reuse_ratio))
  cat(sprintf("\n  THE CALIBRATION SHAPE: %d gradients, us PER OBSERVATION\n", NOBS))
  cat("  (the count below is P_model = length(pars), NOT what an optimiser moves)\n")
  cat("                                  fresh     reused    saved     in C++   vs reused\n")
  for (nm in c("p1", "p4")) {
    f <- out[[paste0("gradN_", nm, "_fresh_us")]]
    r <- out[[paste0("gradN_", nm, "_reuse_us")]]
    bt <- out[[paste0("gradN_", nm, "_batch_us")]]
    cat(sprintf("    length(pars) = %-3d           %8.1f %10s %8s %10s %10s\n",
                length(PSET[[nm]]), f,
                if (is.na(r)) "-" else sprintf("%.1f", r),
                if (is.na(r)) "-" else sprintf("-%.0f%%", 100 * (f - r) / f),
                if (is.na(bt)) "-" else sprintf("%.2f", bt),
                if (is.na(bt) || is.na(r)) "-" else sprintf("%.1fx", r / bt)))
  }
  if (!is.na(out$batch_prepare_us)) {
    cat(sprintf("      leaf_batch() prepare        %8.0f us once, for all %d\n",
                out$batch_prepare_us, NOBS))
  }
  cat(sprintf("      marginal cost per par       %8.1f\n", out$grad_per_par_us))
  if (!is.na(out$grad1_us))
    cat(sprintf("      construction is %.0f%% of a 1-par gradient  <-- issue #52\n",
                100 * out$ctor_us / out$grad1_us))
}
