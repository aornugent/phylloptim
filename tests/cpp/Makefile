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

.PHONY: all build clean golden bench
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

all: build
	./test_leaf
	./test_golden $(GOLDEN_ARGS)

build: test_leaf test_golden

# Not built or run by `make all` -- it takes a few seconds and measures the
# machine as much as the code. Run it deliberately when a change might have
# moved the cost of a solve.
bench: bench_solve bench_gradient
	./bench_solve
	./bench_gradient

test_leaf: test_leaf.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

test_golden: test_golden.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

bench_solve: bench_solve.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

bench_gradient: bench_gradient.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

# Regenerate the golden file. Only ever run this deliberately: it overwrites the
# recorded behaviour, so a refactor that changed results would be rubber-stamped.
golden: test_golden
	./test_golden --generate

clean:
	rm -f test_leaf test_golden bench_solve bench_gradient
