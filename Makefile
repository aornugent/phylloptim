PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

# Developer entry points. Same shape as plant's and odelia's Makefiles, plus the
# targets this package needs that they do not: it has both an R layer and a C++
# core that must keep building without R (hazard 9 in .claude/CLAUDE.md).
#
# The split worth knowing before you pick a target:
#
#   make test        the R layer. Needs R.
#   make test-cpp    the model, as plain C++. Needs no R.
#   make test-cmake  the model again, through the route a C++ or Python
#                    consumer takes, including install and find_package.
#
# `make check-all` runs all three, and is what to run before pushing.

all: RcppR6 compile

# --- generated code ----------------------------------------------------------
#
# inst/RcppR6_classes.yml is the source of truth for the R bindings. Edit that,
# then run this -- never edit R/RcppR6.R, src/RcppR6.cpp or
# inst/include/phylloptim/RcppR6_*.hpp, which all carry a "do not edit by hand" banner.
# The generated files are committed, so a consumer never needs RcppR6; CI
# regenerates and diffs to catch a commit that forgot to.
RcppR6:
	$(RSCRIPT) -e "library(methods); RcppR6::RcppR6()"

# src/RcppExports.cpp and R/RcppExports.R, from anything marked Rcpp::export.
attributes:
	$(RSCRIPT) -e "Rcpp::compileAttributes()"

compile:
	$(RSCRIPT) -e 'pkgbuild::compile_dll(compile_attributes = FALSE, debug = FALSE)'

# Also regenerates the Rcpp exports, unlike `compile`.
full_compile:
	$(RSCRIPT) -e 'pkgbuild::compile_dll(debug = FALSE)'

roxygen:
	@mkdir -p man
	$(RSCRIPT) -e "library(methods); devtools::document()"

# Everything generated, from scratch, in dependency order.
rebuild: clean RcppR6 full_compile roxygen

# --- tests -------------------------------------------------------------------

test: all
	$(RSCRIPT) -e 'library(methods); devtools::test()'

# The model as plain C++: no R, no Rcpp, no test framework. This is the
# regression baseline, and it is blind to the R layer -- tests/testthat/ ties the
# R side back to the same golden points.
#
# `bench` is built too, because it is NOT part of `make -C tests/cpp` and CI
# builds it: a rename that misses bench_solve.cpp gives a green local run and
# three red CI jobs.
test-cpp:
	$(MAKE) -C tests/cpp
	$(MAKE) -C tests/cpp bench_solve

# The same headers through CMake, which is what a C++ or Python consumer uses.
# Covers three things the Makefile above cannot: that the install rules ship the
# right headers and exclude the Rcpp-dependent ones, that find_package produces a
# usable target, and that none of it needs R.
#
# odelia is not findable by CMake, so point at a checkout. A sibling one is the
# default, which is the layout of the plant-family tree.
ODELIA_INC ?= $(CURDIR)/../odelia/inst/include

test-cmake:
	cmake -B build-cmake -DPHYLLOPTIM_ODELIA_INCLUDE_DIR="$(ODELIA_INC)"
	cmake --build build-cmake
	ctest --test-dir build-cmake --output-on-failure

check-all: test test-cpp test-cmake

# Time a collar solve. Not part of any test target: it measures the machine as
# much as the code.
#
# ⚠️ Interleave, never A/B sequentially. Between-process noise is ~±0.1 us, and a
# sequential pair once reported a refactor as 4% faster when interleaving said it
# was 1.7% slower. Build both binaries, keep both, alternate.
bench:
	$(MAKE) -C tests/cpp bench

# ⚠️ Regenerates the recorded behaviour of the model. Run this ONLY when a change
# to the numbers is intended, and say in the commit message what moved and by how
# much -- running it after an accidental change rubber-stamps the change. Guide to
# magnitudes: ~1e-16 is reassociation, ~1e-4 is a real difference.
golden:
	$(MAKE) -C tests/cpp golden

# --- package -----------------------------------------------------------------

install:
	R CMD INSTALL .

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

vignettes:
	$(RSCRIPT) -e "devtools::build_vignettes()"

# The C++ API. Doxygen rather than roxygen because the headers are the
# substantive documentation; roxygen handles the R side, via `make roxygen`.
docs:
	doxygen

clean:
	rm -f src/*.o src/*.so src/*.o.tmp
	rm -rf build-cmake
	$(MAKE) -C tests/cpp clean

.PHONY: all RcppR6 attributes compile full_compile roxygen rebuild \
        test test-cpp test-cmake check-all bench golden \
        install build check vignettes docs clean
