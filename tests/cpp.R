# Run the C++ suite under `R CMD check` (#12).
#
# This package ships headers and nothing else -- no R code, no compiled code --
# so without this file `R CMD check` validates essentially nothing. The point is
# not to duplicate the GitHub Actions workflow: it is that a `LinkingTo: phylloptim`
# consumer's own `R CMD check`, on their machine and their toolchain, tells them
# when a header stops compiling. That is the failure mode that matters for a
# header-only package, and it is invisible to every other check R runs here.
#
# Two deliberate choices:
#
#   * Compile with R's *configured* compiler (`R CMD config CXX20`) rather than
#     whatever `c++` is on PATH, because that is the compiler a consumer builds
#     with.
#   * Compile against the *installed* headers rather than ../../inst/include.
#     Under `R CMD check` those are different directories, and the installed one
#     is what a consumer gets -- so this also catches an install that failed to
#     ship a header the others include.
#
# WHY NOT JUST CALL `make -C tests/cpp`, which is what a developer runs? Because
# `tests/cpp/Makefile` uses GNU extensions (`?=`, `$(shell)`, `$(wildcard)`,
# `ifeq`), and `R CMD check` scans every Makefile in the tarball and warns about
# them. The sanctioned way to silence that is `SystemRequirements: GNU make` --
# but that would be a false statement about this package. Installing it requires
# no make whatsoever; only its own test harness does, and a `LinkingTo: phylloptim`
# consumer would inherit a declared dependency they do not have. So the Makefile
# stays for developers and is left out of the tarball (see .Rbuildignore), and
# the two translation units are compiled directly here. It is two compiler calls.

cpp_dir <- if (dir.exists("cpp")) "cpp" else file.path("tests", "cpp")

if (!dir.exists(cpp_dir)) {
  stop("cannot find the C++ suite; looked for 'cpp' and 'tests/cpp' from ",
       getwd())
}

# --- Is there a toolchain to test with? -------------------------------------
#
# A binary-only R installation with no Rtools / no command-line tools has no
# compiler, and there is nothing this file can meaningfully do there. Skip
# loudly rather than fail: a missing toolchain is not a broken header.

bail <- function(...) {
  message("SKIP: ", ..., "\n  The C++ suite was not run.")
  quit(save = "no", status = 0)
}

rcmd <- function(what) {
  paste(system2(file.path(R.home("bin"), "R"), c("CMD", "config", what),
                stdout = TRUE, stderr = FALSE), collapse = " ")
}

cxx <- rcmd("CXX20")
cxxstd <- rcmd("CXX20STD")
if (!nzchar(trimws(cxx))) {
  bail("R has no C++20 compiler configured (`R CMD config CXX20` is empty).")
}

# --- Where the headers live --------------------------------------------------
#
# BH and odelia are LinkingTo, so they are installed whenever this is checked.

deps <- vapply(c("phylloptim", "BH", "odelia"),
               function(p) system.file("include", package = p),
               character(1))
missing <- names(deps)[!nzchar(deps)]
if (length(missing)) {
  bail("headers not found for: ", paste(missing, collapse = ", "), ".")
}

# --- Bit-exact where that is achievable, tolerant everywhere else -----------
#
# The golden file was generated on macOS/arm64 and is bit-exact only there;
# libm's exp/pow are not bit-reproducible across platforms and never were. This
# mirrors the same decision in .github/workflows/cpp-tests.yml -- see the long
# comment above main() in tests/cpp/test_golden.cpp for why it is per-field and
# why it is not a loosening.

sysinfo <- Sys.info()
on_generating_platform <-
  identical(unname(sysinfo[["sysname"]]), "Darwin") &&
  identical(unname(sysinfo[["machine"]]), "arm64")
golden_args <- if (on_generating_platform) character() else "--cross-platform"

# --- Build and run -----------------------------------------------------------
#
# Flags mirror tests/cpp/Makefile. They are deliberately not shared with it --
# the point of compiling here is to use R's toolchain rather than the Makefile's
# -- so keep the two in step if you change either.

message("Building the C++ suite with:")
message("  CXX     ", cxx, " ", cxxstd)
message("  leaf    ", deps[["phylloptim"]])
message("  BH      ", deps[["BH"]])
message("  odelia  ", deps[["odelia"]])
message("  golden  ", if (length(golden_args)) golden_args else "bit-exact")

owd <- setwd(cpp_dir)   # both golden readers take their file relatively
on.exit(setwd(owd), add = TRUE)

includes <- c("-I", shQuote(deps[["phylloptim"]]),
              "-isystem", shQuote(deps[["odelia"]]),
              "-isystem", shQuote(deps[["BH"]]))

# ⚠️ Every program under tests/cpp/ that has a golden file belongs in this list, and
# in tests/cpp/Makefile, CMakeLists.txt and cpp-tests.yml. #64 is what a program in
# none of them turns into.
for (prog in c("test_leaf", "test_golden", "test_primitives")) {
  src <- paste0(prog, ".cpp")
  message("\n== ", src)
  status <- system2(cxx, c(cxxstd, "-O2", "-Wall", "-Wextra",
                           "-Wno-unused-parameter", includes,
                           shQuote(src), "-o", shQuote(prog)))
  if (status != 0) {
    stop(src, " did not compile (exit ", status, "). A header in the installed\n",
         "  package is broken for anything that LinkingTo's it.")
  }
  status <- system2(file.path(".", prog), golden_args)
  if (status != 0) {
    stop(prog, " failed (exit ", status, ").\n",
         "  If this was test_primitives, read its per-TIER table first: the lowest\n",
         "  tier that moved is the cause and the ones above it are consequences.\n",
         "  If this was the golden comparison, read the worst-relative-difference\n",
         "  line it printed first: ~1e-15 is cross-platform libm and not a\n",
         "  regression, ~1e-4 is a real change in behaviour.")
  }
}
