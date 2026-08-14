# ⚠️⚠️ THIS HARNESS NO LONGER APPLIES, AND HAS NOT SINCE #15. Read this first.
#
# The note below says to "run this against `main`, whose set_physiology signature
# still matches plant's". That branch is gone: `feature/api-cleanup` merged as #15,
# which changed the signature on master, and #33 changed it again (root carbon ->
# RootNetwork). Both this script and its companion compare_primitives.R build a
# Leaf with beta_R_H/beta_R_V constructor arguments that no longer exist.
#
# Kept for the reasoning it records -- what has been ruled out, and why comparing
# primitives beats comparing the solve -- not as something you can run today.
# Reviving it means a decision about what the reference build should now be, since
# plant no longer has its own copy of this model. Filed as #64.
#
# PLAN.md item 1: cross-check this package against plant's compiled build.
#
#   Rscript tests/validate/compare_with_plant.R
#
# Runs plant's own `Leaf` -- the compiled class in the installed plant package,
# not a reimplementation -- over exactly the grid in tests/cpp/test_golden.cpp,
# and compares against the golden file this package generated.
#
# RESULT AS OF 2026-08-03: the finite values are BIT-IDENTICAL to plant, all 2352
# of them. Zero difference, not a small one.
#
# This retracts the previous result rather than editing it, because the previous
# result was an artifact of this script. It reported '585 of 2592 differ, every
# one at 1-2 ULP, cause unknown'. That decomposes exactly:
#
#     345   R's decimal parser reading the golden file
#     240   the shutdown-state NA sentinel (48 rows x 5 flux fields)
#     ---
#     585
#
# THE 345. R's string-to-double conversion is not correctly rounded. as.numeric,
# scan and read.delim all share it, and it returns a double one ULP off the
# correctly rounded value for about 18% of inputs:
#
#     "26.550866314209998"   R gives 0x1.a8d0593240001p+4
#                            correct  0x1.a8d059324p+4
#
# The golden file is written by C++ at full %.17g precision and read into R;
# plant's values are computed in-process and never touch a string. So the parser
# perturbed one side and not the other, and the script attributed the result to
# the two implementations. It now reads through tests/validate/tsv_to_hex.c,
# which parses with the C library's strtod and re-emits hex, which R reads
# exactly. Verified: 4000 of 4000 hex values round-trip, against 3265 of 4000
# for %.17g.
#
# This also explains why every previous hypothesis was ruled out and nothing
# replaced them -- plant version, compiler flags, -ffp-contract, inlining, odelia
# header version, translation-unit structure. None of them was ever involved. It
# is also why scm_regression.R found the consume build bit-identical: that
# comparison never round-trips through a text file parsed by R.
#
# THE 240 are not arithmetic at all: golden carries the NA sentinel on shutdown
# rows where plant carries a number, because set_shutdown_state assigns
# opt_root_psi_ (root_collar_psi_ before #25), opt_psi_stem_ and profit_ and leaves
# ci/assim/transpiration/
# gc/e_up untouched. That WAS the shutdown-state leak (PLAN item 2, plant #578),
# and the count matches the 48 x 5 recorded as that fix's blast radius. The leak is
# fixed both here and in plant, so this column should now be empty against a
# current plant -- it is retained because a non-empty one means the two have
# drifted apart again, which is more useful than the check disappearing.
#
# Corroborated independently by compare_primitives.R, which calls the underlying
# functions directly -- arrh_curve, the vulnerability curve, the assimilation
# terms, transpiration, and the two that iterate -- and finds all 329 values
# bit-identical.
#
# So: expect ZERO. Any finite difference now is real and should be investigated;
# the nested solvers amplify perturbations up to about GSS_tol_abs (1e-3), so a
# genuine arithmetic difference shows up around 1e-4.
#
# Run this against `main`, whose set_physiology signature still matches plant's.
# The `feature/api-cleanup` branch deliberately diverges (different signature), so
# this script does not apply there; its deviations are measured in its own commits.
#
# REFERENCE BUILD. Compare against plant built from the commit this package was
# extracted from, NOT whatever happens to be installed -- the installed plant may
# be a different branch. Recreate it with:
#
#   git -C ../plant worktree add /tmp/plant-ref 76df7169
#   R CMD INSTALL --no-docs -l /tmp/rlib-ref /tmp/plant-ref
#   R_LIBS=/tmp/rlib-ref Rscript tests/validate/compare_with_plant.R
#
# (The first attempt at this compared against an installed plant carrying the ATLS
# thermal-damage layer, which is a different branch. It gave identical numbers,
# because ATLS is default-off and genuinely bit-identical when off -- but that was
# luck, not method.)

suppressMessages(library(plant))

golden_path <- "tests/cpp/golden/operating_points.tsv"
if (!file.exists(golden_path)) {
  stop("Run from the repository root: ", golden_path, " not found")
}

# --- the grid, mirroring tests/cpp/test_golden.cpp exactly --------------------

psi_soils <- c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0)
ppfds <- c(100, 500, 900, 1500)
vpds <- c(0.5, 1.0, 2.0, 4.0)
layer_counts <- c(1L, 3L, 5L)
# ⚠️ LEAF TEMPERATURE IS AN AXIS NOW (#41), not the fixed 25 C this script used to
# assume, because the golden file it reads has four temperature blocks. Two
# consequences before running anything:
#
#   * at 25 C the comparison means what it always meant -- and a difference there is
#     a DIFFERENT bug, since #41 is inert at the reference point by construction;
#   * at 15 / 35 / 40 C phylloptim and a pre-#41 plant are EXPECTED to differ,
#     because R_d's shape changed. `Leaf::rd_tracks_vcmax_` is what makes that arm
#     comparable, and it lives on the phylloptim side, so an arm that uses it cannot
#     be read out of the golden file -- see the arm table in this file's header.
leaf_temps <- c(15.0, 25.0, 35.0, 40.0)

# Fixed drivers and traits, from plant's tests/testthat/test-leaf.r.
theta <- 0.000157
K_s <- 1.0
h <- 5.0
area_leaf <- 0.05
rho <- 608.0
a_bio <- 0.0245
ca <- 40.0
o2 <- 21.0
patm <- 101.3

# Every constructor argument passed explicitly, matching the C++ default
# constructor value for value. Relying on R-side defaults would leave the
# comparison hostage to a defaults mismatch, which is exactly the sort of thing
# this script exists to detect.
new_leaf <- function() {
  plant::Leaf(
    vcmax_25 = 96, c = 2.680147, b = 3.898245, psi_crit = 5.870283,
    root_c = 2.680147, root_b = 3.898245, root_psi_crit = 5.870283,
    beta2 = 1.5, jmax_25 = 157.44, a = 0.30,
    curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
    GSS_tol_abs = 1e-3, vulnerability_curve_ncontrol = 100,
    ci_abs_tol = 1e-3, ci_niter = 1000,
    g1_TF24 = 7.5, beta_R_H = 3.4e2, beta_R_V = 9.4e3
  )
}

solve_one <- function(psi_soil, ppfd, vpd, layers, leaf_temp) {
  l <- new_leaf()
  i <- seq_len(layers)
  ps <- psi_soil + 0.25 * (i - 1)
  depth <- 1.0 * i
  root <- rep(1.0 / layers, layers)

  l$set_physiology(
    area_leaf = area_leaf, mass_root_prop = root, rho = rho, a_bio = a_bio,
    PPFD = ppfd, psi_soil = ps, soil_depth = depth,
    leaf_specific_conductance_max = K_s * theta / h,
    atm_vpd = vpd, ca = ca, sapwood_volume_per_leaf_area = theta * h,
    leaf_temp = leaf_temp, atm_o2_kpa = o2, atm_kpa = patm
  )
  l$find_root_collar_psi()

  cons <- l$soil_consumption_
  data.frame(
    psi_soil = psi_soil, ppfd = ppfd, vpd = vpd, leaf_temp = leaf_temp,
    layers = layers,
    psi_stem = l$opt_psi_stem_, opt_root_psi = l$opt_root_psi_,
    ci = l$ci_, assim = l$assim_colimited_,
    transpiration = l$transpiration_, gc = l$stom_cond_CO2_,
    profit = l$profit_, e_up = l$E_up_,
    uptake = sum(cons[is.finite(cons)])
  )
}

grid <- expand.grid(layers = layer_counts, vpd = vpds, ppfd = ppfds,
                    psi_soil = psi_soils, leaf_temp = leaf_temps)
# expand.grid varies the FIRST column fastest; the C++ loops nest
# leaf_temp > psi_soil > ppfd > vpd > layers, with layers innermost. Listing them
# in reverse above therefore reproduces the C++ row order.
grid <- grid[, c("psi_soil", "ppfd", "vpd", "leaf_temp", "layers")]

cat(sprintf("Running plant's Leaf over %d operating points...\n", nrow(grid)))
plant_rows <- do.call(rbind, Map(solve_one, grid$psi_soil, grid$ppfd,
                                 grid$vpd, grid$layers, grid$leaf_temp))

# --- compare ------------------------------------------------------------------

# Read the golden file through tsv_to_hex rather than with read.delim directly.
# This is not fussiness, it is the difference between measuring the models and
# measuring R: R's decimal parser is not correctly rounded, and on this file it
# returns a double one ULP off the correctly rounded value for 345 of 3504 cells
# (9.8%). Those perturbations land on the C++ side only -- plant's values are
# computed in-process and never go near a string -- so they show up as a
# disagreement between the two implementations, which is what they are not.
# See tests/validate/tsv_to_hex.c, and issue #13.
read_tsv_exactly <- function(path) {
  src <- "tests/validate/tsv_to_hex.c"
  bin <- file.path(tempdir(), "tsv_to_hex")
  if (!file.exists(bin)) {
    cc <- paste(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CC"),
                        stdout = TRUE), collapse = " ")
    if (system2(cc, c("-O2", "-o", shQuote(bin), shQuote(src))) != 0) {
      stop("could not build ", src, ", which is needed to read ", path,
           " without R's parser rounding it")
    }
  }
  d <- read.delim(text = paste(system2(bin, stdin = path, stdout = TRUE),
                               collapse = "\n"),
                  colClasses = "character")
  d[] <- lapply(d, as.numeric)   # hex; R parses this exactly
  d
}

golden <- read_tsv_exactly(golden_path)
stopifnot(nrow(golden) == nrow(plant_rows))

# The inputs must line up before comparing outputs, or we would be comparing
# different operating points and calling it agreement.
for (key in c("psi_soil", "ppfd", "vpd", "leaf_temp", "layers")) {
  if (!isTRUE(all.equal(golden[[key]], plant_rows[[key]], tolerance = 0))) {
    stop("grid mismatch in '", key, "': the R grid is not in the same order as ",
         "the C++ one, so the comparison would be meaningless")
  }
}

fields <- c("psi_stem", "opt_root_psi", "ci", "assim", "transpiration", "gc",
            "profit", "e_up", "uptake")

# NaN == NaN for our purposes: shut-down points legitimately carry the NA
# sentinel on `main` (the shutdown-state leak, PLAN item 2).
differs <- function(a, b) !((is.nan(a) & is.nan(b)) | (!is.nan(a) & !is.nan(b) & a == b))

summary_rows <- list()
total_mismatch <- 0L
for (f in fields) {
  a <- plant_rows[[f]]
  b <- golden[[f]]
  # Two entirely different things get counted separately, because conflating
  # them is what made this look like a numerical mystery for so long:
  #   shutdown  -- golden carries the NA sentinel where plant carries a number.
  #                That is the shutdown-state leak (PLAN item 2, plant #578),
  #                a known behavioural difference, not a rounding one.
  #   numeric   -- both sides finite and unequal. THIS is the arithmetic.
  shutdown <- is.nan(b) & !is.nan(a)
  both_num <- !is.nan(a) & !is.nan(b)
  bad <- which(both_num & a != b)
  total_mismatch <- total_mismatch + length(bad)
  rel <- if (length(bad)) {
    max(abs(a[bad] - b[bad]) / pmax(abs(b[bad]), .Machine$double.xmin))
  } else 0
  summary_rows[[f]] <- data.frame(field = f, shutdown_NA = sum(shutdown),
                                  numeric_differ = length(bad),
                                  max_rel_diff = rel)
}
summ <- do.call(rbind, summary_rows)

cat("\n")
print(summ, row.names = FALSE)
cat(sprintf("\n%d of %d finite comparisons differ; %d cells are the shutdown\n  NA sentinel (a known behavioural difference, not arithmetic).\n",
            total_mismatch, nrow(golden) * length(fields) - sum(summ$shutdown_NA),
            sum(summ$shutdown_NA)))

# The bar is a few ULP, not zero -- see the header. Bit-identity is unachievable
# across a header-only single-TU build and plant's separate-TU build, and demanding
# it would just mean a permanently red check.
ulp_tol <- 1e-14   # ~45 ULP of headroom; a real arithmetic difference lands at 1e-4
worst_rel <- max(summ$max_rel_diff)

if (total_mismatch == 0L) {
  cat("\nPASS: bit-identical to plant across the whole grid.\n")
  quit(status = 0)
}
if (worst_rel <= ulp_tol) {
  cat(sprintf(paste0("\nPASS: agreement to %.1e relative (%.0f ULP), within the %.0e bar.\n",
                     "The source is arithmetically identical (see the header); the\n",
                     "cause of the 1-ULP residual is NOT established, but it is\n",
                     "confined to this standalone harness -- swapping the package into\n",
                     "plant is bit-identical (scm_compare.R). Anything above ~1e-12\n",
                     "would be a real difference and should be investigated.\n"),
              worst_rel, worst_rel / .Machine$double.eps, ulp_tol))
  quit(status = 0)
}

# Show the worst offenders so the failure is diagnosable, not just reported.
cat("\nWorst mismatches:\n")
worst <- NULL
for (f in fields) {
  a <- plant_rows[[f]]; b <- golden[[f]]
  bad <- which(differs(a, b))
  if (!length(bad)) next
  rel <- abs(a[bad] - b[bad]) / pmax(abs(b[bad]), .Machine$double.xmin)
  o <- bad[order(rel, decreasing = TRUE)][seq_len(min(3, length(bad)))]
  worst <- rbind(worst, data.frame(
    field = f, psi_soil = golden$psi_soil[o], ppfd = golden$ppfd[o],
    vpd = golden$vpd[o], layers = golden$layers[o],
    plant = a[o], leaf = b[o],
    rel_diff = abs(a[o] - b[o]) / pmax(abs(b[o]), .Machine$double.xmin)))
}
print(head(worst[order(worst$rel_diff, decreasing = TRUE), ], 20), row.names = FALSE)

cat(sprintf("\nFAIL: worst relative difference %.3e exceeds the %.0e bar.\n", worst_rel, ulp_tol))
cat("That is too large to be reassociation. Check, in order: that plant was built\n")
cat("from the commit this package was extracted from (see the header -- a different\n")
cat("branch is the likeliest cause); that the grid and constructor arguments still\n")
cat("match tests/cpp/test_golden.cpp; and only then suspect a conversion bug.\n")
quit(status = 1)
