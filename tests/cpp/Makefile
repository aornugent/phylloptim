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
DEP_INC  = $(ODELIA_INC) $(BH_INC)
SYS_INC  = $(foreach d,$(filter-out /usr/include /usr/include/,$(strip $(DEP_INC))),-isystem $(d))
INCLUDES = -I$(PHYLLOPTIM_INC) $(SYS_INC)

# The collar's condition records a reverse pass, so the suite needs XAD's tape:
# the two storage-class defines, which must match odelia's and plant's Makevars,
# and the one translation unit odelia compiles the tape into. The installed
# package ships headers only, so this comes from a checkout.
XADFLAGS = -DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE
ODELIA_SRC ?= $(ODELIA_INC)/../../src
TAPE_SRC = $(ODELIA_SRC)/Tape.cpp

.PHONY: all build clean golden bench
# The suite's own headers are prerequisites too: root_network.hpp holds the two
# ex-`Leaf`-default constants the golden file's bit-exactness depends on, so a
# stale rebuild after editing it would be a silently wrong baseline.
HEADERS = $(wildcard $(PHYLLOPTIM_INC)/phylloptim/*.hpp) $(PHYLLOPTIM_INC)/phylloptim.hpp \
          $(wildcard *.hpp)

# The golden file is bit-exact only on the platform that generated it, so the
# comparison that means anything depends on where it runs: bit-exact on
# macOS/arm64, per-field tolerances everywhere else. See the comment above main()
# in test_golden.cpp for why the tolerant mode is correct and not a loosening.
#
# Defaulted by platform rather than left to the caller. Empty everywhere, the
# form this had, makes `make` fail on every other platform for a reason that is
# never a regression -- so a reader learns to ignore it, and the one time it IS a
# regression nothing can say so. A stale golden then went unnoticed for ten days.
ifeq ($(shell uname -s)/$(shell uname -m),Darwin/arm64)
  GOLDEN_ARGS ?=
else
  GOLDEN_ARGS ?= --cross-platform
endif

all: build
	./test_leaf
	./test_leaf_gradient
	./test_transpose
	./test_supplied_rows
	./test_golden $(GOLDEN_ARGS)

build: test_leaf test_leaf_gradient test_transpose test_supplied_rows \
       test_golden

# Not built or run by `make all` -- it takes a few seconds and measures the
# machine as much as the code. Run it deliberately when a change might have
# moved the cost of a solve.
bench: bench_solve
	./bench_solve

test_leaf: test_leaf.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

# The gradient against a difference of the solve. Slower than test_leaf -- it
# re-solves the whole leaf four times per input per state -- and it is the one
# thing that referees the derivative against something that is not itself.
test_leaf_gradient: test_leaf_gradient.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

# The transpose identity, which needs no reference gradient: one tangent against
# one sweep. It is the only thing here that can see a wrong SUPPLIED row, since a
# difference of the recorded step cannot. 0.01 s.
test_transpose: test_transpose.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

# The two derivatives the leaf hands over as plain numbers, against a difference of
# the forward model each describes. A difference of anything downstream is
# identically zero on those columns, so this is the only place they can be caught.
test_supplied_rows: test_supplied_rows.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

test_golden: test_golden.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

bench_solve: bench_solve.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

# One-off measurement programs. Not built or run by anything; `make probe_x`
# builds probe_x.cpp with the suite's own flags so a probe measures the same code
# the tests do.
probe_%: probe_%.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(XADFLAGS) $(INCLUDES) $< $(TAPE_SRC) -o $@

# Regenerate the golden file. Only ever run this deliberately: it overwrites the
# recorded behaviour, so a refactor that changed results would be rubber-stamped.
golden: test_golden
	./test_golden --generate

clean:
	rm -f test_leaf test_leaf_gradient test_transpose test_supplied_rows \
	      test_golden bench_solve
