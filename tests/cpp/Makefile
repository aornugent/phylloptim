# Build and run the leaf test suite as plain C++ -- no R, no linking, and since
# odelia #44 nothing standing in for R either.
#
#   make          build and run
#   make build    build only
#   make clean
#
# Header locations are discovered from R if it is on the PATH (the normal case
# once BH and odelia are installed), and can be overridden:
#
#   make BH_INC=/opt/homebrew/include ODELIA_INC=../../../odelia/inst/include
#
# Falling back to a sibling odelia checkout means this works in the plant-family
# tree before anything is installed.

CXX      ?= c++
CXXSTD   ?= -std=c++20
OPT      ?= -O2
# -Werror=switch only: operating_point_kind_name switches over OperatingPointKind
# with no default, so a kind added without a name stops this build.
WARN     ?= -Wall -Wextra -Wno-unused-parameter -Werror=switch

PHYLLOPTIM_INC   ?= ../../inst/include
BH_INC     ?= $(shell Rscript -e 'cat(system.file("include", package="BH"))' 2>/dev/null)
ODELIA_INC ?= $(shell Rscript -e 'cat(system.file("include", package="odelia"))' 2>/dev/null)

ifeq ($(strip $(BH_INC)),)
  BH_INC := /opt/homebrew/include
endif
ifeq ($(strip $(ODELIA_INC)),)
  ODELIA_INC := ../../../odelia/inst/include
endif

# Dependency include dirs, as -isystem flags.
#
# /usr/include is filtered out deliberately. Passing `-isystem /usr/include`
# moves it earlier in the system-header search order, and libstdc++'s <cmath>
# then cannot resolve its `#include_next <math.h>` -- the build dies with
#
#   /usr/include/c++/13/cmath:47:15: fatal error: math.h: No such file or directory
#
# which reads like a broken toolchain and is not. It bites on Linux, where
# apt's libboost-dev headers land in /usr/include and CI passes BH_INC=/usr/include.
# The compiler already searches /usr/include, so dropping the flag costs nothing.
# ⚠️ PHYLLOPTIM NOW HAS A REVERSE-MODE PATH, and nothing will tell you if these
# disagree. XAD's active tape is a __thread variable whose STORAGE CLASS is
# not part of the mangled name, so a mismatch against odelia's and plant's
# Makevars links cleanly and is undefined behaviour. These two must stay equal to
# theirs. Forward mode is tapeless and needed neither, which is why
# src/Makevars was exempt -- it is not exempt once the R layer reaches the tape.
XADFLAGS = -DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE
# The one translation unit odelia compiles the tape into. The installed package
# ships headers only, so this comes from a checkout -- and ODELIA_INC resolves to
# the INSTALLED headers whenever odelia is installed, which is the normal case
# here, so deriving this from it alone finds no Tape.cpp. Only test_transpose
# needs the tape, so the rest of the suite builds and runs and the transpose
# identity is the one target that silently stops being checked.
ODELIA_SRC ?= $(ODELIA_INC)/../../src
ifeq ($(wildcard $(ODELIA_SRC)/Tape.cpp),)
  ODELIA_SRC := ../../../odelia/src
endif
TAPE_SRC = $(ODELIA_SRC)/Tape.cpp

DEP_INC  = $(ODELIA_INC) $(BH_INC)
SYS_INC  = $(foreach d,$(filter-out /usr/include /usr/include/,$(strip $(DEP_INC))),-isystem $(d))
INCLUDES = -I$(PHYLLOPTIM_INC) $(SYS_INC)

.PHONY: all build clean golden psi-stem-golden primitives-golden bench
# The suite's own headers are prerequisites too: root_network.hpp holds the two
# ex-`Leaf`-default constants the golden file's bit-exactness depends on, so a
# stale rebuild after editing it would be a silently wrong baseline.
HEADERS = $(wildcard $(PHYLLOPTIM_INC)/phylloptim/*.hpp) $(PHYLLOPTIM_INC)/phylloptim.hpp \
          $(wildcard *.hpp)

# Extra arguments for the golden comparison. Empty means bit-exact, which is what
# you want locally and on the platform that generated the file. CI passes
# GOLDEN_ARGS=--cross-platform on platforms other than that one -- see the comment
# above main() in test_golden.cpp for why that is correct and not a loosening.
GOLDEN_ARGS ?=

# ⚠️ EVERY TEST PROGRAM MUST BE NAMED HERE. This Makefile lists its sources
# explicitly, and that is exactly how tests/validate/ rotted through two interface
# changes without anyone noticing (#64): a file no build touches compiles in
# nobody's head. Adding a .cpp under tests/cpp/ without adding it to `build`, to
# `all`, to CMakeLists.txt and to cpp-tests.yml leaves it in the same state, with a
# newer date on it.
all: build
	./test_leaf
	./test_supplied_rows
	./test_transpose
	./test_golden $(GOLDEN_ARGS)
	./test_primitives $(GOLDEN_ARGS)

build: test_leaf test_supplied_rows test_transpose test_golden test_primitives

# Not built or run by `make all` -- it takes a few seconds and measures the
# machine as much as the code. Run it deliberately when a change might have
# moved the cost of a solve.
bench: bench_solve bench_gradient
	./bench_solve
	./bench_gradient

test_leaf: test_leaf.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

# The one row this surface SUPPLIES rather than records, against a difference of
# the function it claims to describe. A difference of anything downstream is
# identically zero on that column, so this is the only place it can be caught.
test_supplied_rows: test_supplied_rows.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

# The transpose identity, which needs no reference gradient: one tangent against
# one sweep. It is the only thing here that can see a wrong SUPPLIED row, and the
# only one that can referee the collar responses at the shipped spline
# resolution. The one program here that links the tape.
test_transpose: test_transpose.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

# Not in `all`: it answers a design question rather than guarding a property, and
# it needs the tape the same way test_transpose does.
probe_preaccumulation: probe_preaccumulation.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

test_golden: test_golden.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

test_primitives: test_primitives.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

bench_solve: bench_solve.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

bench_gradient: bench_gradient.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

# Regenerate the golden file. Only ever run this deliberately: it overwrites the
# recorded behaviour, so a refactor that changed results would be rubber-stamped.
golden: test_golden
	./test_golden --generate

# The psi_stem optimiser baseline, on the same terms as `golden` above: only ever
# deliberately. Separate from `golden` so the two files cannot be regenerated by
# accident together -- they have disjoint causes, and a diff in one is evidence
# about a different part of the model than a diff in the other.
psi-stem-golden: test_golden
	./test_golden --generate-optima

# The primitive-level baseline, regenerated on the same terms: only ever
# deliberately, and only on macOS/arm64, which is the platform it is bit-exact on.
primitives-golden: test_primitives
	./test_primitives --generate

clean:
	rm -f test_leaf test_supplied_rows test_transpose test_golden test_primitives bench_solve bench_gradient
