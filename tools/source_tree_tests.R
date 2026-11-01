# The testthat files that read the source tree, run against it.
#
# Under R CMD check the tests run from an unpacked copy that has no README.md,
# no R/ sources and no Doxygen XML, so these files skip there and guard nothing.
# R-CMD-check.yml runs this from the checkout, after `doxygen` and an install,
# and here a skip is a failure: this is the one place they can run.
#
# Usage, from the package root: Rscript tools/source_tree_tests.R
files <- c("test-readme.R", "test-registry.R")
Sys.setenv(NOT_CRAN = "true")

suppressPackageStartupMessages({
  library(testthat)
  library(phylloptim)
})
if (!file.exists(file.path("docs", "xml", "index.xml"))) {
  stop("no Doxygen XML in docs/xml; run `doxygen` at the package root first",
       call. = FALSE)
}
results <- do.call(rbind, lapply(files, function(f) {
  as.data.frame(test_file(file.path("tests", "testthat", f),
                          package = "phylloptim", load_package = "installed",
                          reporter = "summary", stop_on_failure = FALSE))
}))
bad <- results[results$failed > 0 | results$error | results$skipped, ]
if (nrow(bad)) {
  stop("failed, errored or skipped: ",
       paste0(bad$file, ": ", bad$test, collapse = "; "), call. = FALSE)
}
cat(sprintf("%d tests from %d files ran against the source tree, none skipped\n",
            nrow(results), length(files)))
