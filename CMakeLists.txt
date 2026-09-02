# Consume the leaf model from C++ -- or from anything that can call C++ -- with
# no R involved at any point.
#
# This exists so that "the model does not need R" is a thing you can run rather
# than a claim in DESCRIPTION. The package gained an R layer in issue #5, and the
# guarantee that survived it (PLAN item 6a) is that inst/include/ is plain C++:
# no R headers, no Rcpp, nothing that requires an R installation to build or to
# run. An R package's build system cannot demonstrate that, because R is present
# the whole time. This can, and CI runs it on a runner with no R.
#
#   cmake -B build -DPHYLLOPTIM_ODELIA_INCLUDE_DIR=/path/to/odelia/inst/include
#   cmake --build build
#   ctest --test-dir build
#
# In another project:
#
#   add_subdirectory(phylloptim)            # or find_package(phylloptim)
#   target_link_libraries(mine PRIVATE phylloptim::phylloptim)
#
# The target is INTERFACE -- the model is header-only and there is nothing to
# link against. `src/` is the R binding layer and is deliberately not built here;
# it is the one part of this repository that does need R.
#
# See README.md, "Using phylloptim without R", for the Python route, which is this
# target plus pybind11.

cmake_minimum_required(VERSION 3.16)
project(phylloptim VERSION 0.1.0 LANGUAGES CXX
        DESCRIPTION "Leaf gas exchange with explicit plant hydraulics")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# PROJECT_IS_TOP_LEVEL is CMake 3.21; this file claims 3.16, so derive it.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(PHYLLOPTIM_IS_TOP_LEVEL ON)
else()
  set(PHYLLOPTIM_IS_TOP_LEVEL OFF)
endif()

# --- default build type, which the golden file cares about --------------------
#
# CMake's default is no build type at all, which means no -O flag, which means
# -O0. Measured here on macOS/arm64: the golden file is bit-identical at -O1,
# -O2 and -O3, and NOT at -O0 -- worst relative difference 3.47e-15, about 13
# ULP, squarely in this project's "~1e-16 is reassociation, ~1e-4 is a real
# difference" band. The mechanism is the documented one: -O0 does not contract
# a*b+c into an FMA and every optimised level does, so the unoptimised build is
# the odd one out rather than the reference.
#
# This matters beyond CMake, and the developer guide did not say it: the golden
# file's bit-exactness is conditional on the OPTIMISATION LEVEL as well as on the
# platform. tests/cpp/Makefile defaults to -O2 and so never met this. A debug
# build that fails test_golden by ~1e-15 has found nothing.
if(PHYLLOPTIM_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- dependencies ------------------------------------------------------------
#
# Two, both header-only:
#
#   Boost   the TOMS748 root finder. Whatever your platform packages.
#   odelia  the spline interpolator and the vendored XAD forward-mode AD.
#           https://github.com/traitecoevo/odelia -- an R package by
#           distribution, but the headers used here are plain C++, so a git
#           checkout is enough and nothing needs installing or building.
#
# odelia is not findable by CMake, so its include directory is asked for. A
# sibling checkout is the default because that is the layout of the plant-family
# tree this package lives in.

# Config mode only. The FindBoost module is removed in CMake 4, and asking for
# it emits a policy warning on 3.30+; config mode has shipped with Boost since
# 1.70. Where there is no BoostConfig.cmake -- some distribution packages still
# install only the headers -- fall back to an include directory, which is all
# this uses: one root finder, header-only, nothing to link.
find_package(Boost QUIET CONFIG)

set(PHYLLOPTIM_BOOST_INCLUDE_DIR "" CACHE PATH
    "Path to Boost's headers; only needed if find_package(Boost CONFIG) fails")

if(NOT Boost_FOUND)
  if(NOT EXISTS "${PHYLLOPTIM_BOOST_INCLUDE_DIR}/boost/math/tools/toms748_solve.hpp")
    message(FATAL_ERROR
      "Boost not found by find_package(Boost CONFIG), and "
      "PHYLLOPTIM_BOOST_INCLUDE_DIR (${PHYLLOPTIM_BOOST_INCLUDE_DIR}) does not contain "
      "boost/math/tools/toms748_solve.hpp.\n"
      "Install Boost's headers and either let CMake find them or pass "
      "-DPHYLLOPTIM_BOOST_INCLUDE_DIR=/path/to/boost/headers")
  endif()
endif()

set(PHYLLOPTIM_ODELIA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../odelia/inst/include"
    CACHE PATH "Path to odelia's inst/include")

if(NOT EXISTS "${PHYLLOPTIM_ODELIA_INCLUDE_DIR}/odelia/interpolator.hpp")
  message(FATAL_ERROR
    "odelia headers not found under PHYLLOPTIM_ODELIA_INCLUDE_DIR "
    "(${PHYLLOPTIM_ODELIA_INCLUDE_DIR}).\n"
    "Clone https://github.com/traitecoevo/odelia and point at its "
    "inst/include, e.g.\n"
    "  cmake -B build -DPHYLLOPTIM_ODELIA_INCLUDE_DIR=/path/to/odelia/inst/include")
endif()

# --- the target ---------------------------------------------------------------

add_library(phylloptim INTERFACE)
add_library(phylloptim::phylloptim ALIAS phylloptim)

# C++20 is XAD's floor, not a preference.
target_compile_features(phylloptim INTERFACE cxx_std_20)

target_include_directories(phylloptim INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inst/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# SYSTEM so that a consumer building with -Werror is not held responsible for
# warnings in Boost, XAD or the interpolator.
target_include_directories(phylloptim SYSTEM INTERFACE
  $<BUILD_INTERFACE:${PHYLLOPTIM_ODELIA_INCLUDE_DIR}>)

if(Boost_FOUND)
  target_link_libraries(phylloptim INTERFACE Boost::headers)
elseif(NOT PHYLLOPTIM_BOOST_INCLUDE_DIR MATCHES "^/usr/include/?$")
  target_include_directories(phylloptim SYSTEM INTERFACE
    $<BUILD_INTERFACE:${PHYLLOPTIM_BOOST_INCLUDE_DIR}>)
endif()
# /usr/include is deliberately not added as a SYSTEM include, and dropping it
# costs nothing because the compiler searches it anyway. Passing
# `-isystem /usr/include` moves it earlier in the system-header search order, and
# libstdc++'s <cmath> then cannot resolve its `#include_next <math.h>`; the build
# dies with "math.h: No such file or directory", which reads like a broken
# toolchain and is not. It bites on Linux, where apt's libboost-dev headers land
# in /usr/include. tests/cpp/Makefile filters it for the same reason.

# --- the C++ suite ------------------------------------------------------------
#
# The same programs tests/cpp/Makefile builds, and the same golden files.
# Deliberately not a reimplementation of that Makefile: it stays the developer
# entry point, and this is here so the no-R claim is checkable by a machine that
# has never had R on it. Both golden readers take their file from the working
# directory, hence WORKING_DIRECTORY below.
#
# ⚠️ A test program added under tests/cpp/ has to be added HERE as well as to that
# Makefile and to cpp-tests.yml. #64 is what happens otherwise: a harness no build
# touches rots through two interface changes in silence.
#
# The golden file is bit-exact only on the platform that generated it
# (macOS/arm64), so pass --cross-platform anywhere else. See the comment above
# main() in tests/cpp/test_golden.cpp for why that is not a loosening.

option(PHYLLOPTIM_BUILD_TESTS "Build the C++ test suite" ${PHYLLOPTIM_IS_TOP_LEVEL})

if(PHYLLOPTIM_BUILD_TESTS)
  enable_testing()

  foreach(prog test_leaf test_golden test_primitives bench_solve)
    add_executable(${prog} tests/cpp/${prog}.cpp)
    target_link_libraries(${prog} PRIVATE phylloptim::phylloptim)
    target_compile_options(${prog} PRIVATE -Wall -Wextra -Wno-unused-parameter)
  endforeach()

  add_test(NAME leaf_unit COMMAND test_leaf
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp)

  if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    set(PHYLLOPTIM_GOLDEN_ARGS "")
  else()
    set(PHYLLOPTIM_GOLDEN_ARGS "--cross-platform")
  endif()

  add_test(NAME leaf_golden COMMAND test_golden ${PHYLLOPTIM_GOLDEN_ARGS}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp)

  # The primitive-level baseline (#64). Same per-platform policy, and the same
  # variable, because it is bit-exact on exactly the same platform and for exactly
  # the same reason -- libm's exp/pow.
  add_test(NAME leaf_primitives COMMAND test_primitives ${PHYLLOPTIM_GOLDEN_ARGS}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp)

  # bench_solve is built but not registered as a test: it is a timing harness,
  # so it would be a flaky gate on a shared runner. Building it keeps it from
  # rotting, which is the same reason the GitHub workflow compiles it.
endif()

# --- install ------------------------------------------------------------------
#
# Headers only. RcppR6_pre/post/support.hpp are generated into inst/include/phylloptim/
# by RcppR6 and DO include Rcpp -- they are excluded here, along with leaf.h, the
# R-facing umbrella that includes them. Installing them would hand a non-R
# consumer headers they cannot compile, and would blur the one-way dependency
# that PLAN item 6a rests on. <phylloptim.hpp> is the entry point; <phylloptim.h> is not.

install(DIRECTORY inst/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
          PATTERN "*.hpp"
          PATTERN "RcppR6_*.hpp" EXCLUDE)

install(TARGETS phylloptim EXPORT phylloptimTargets)
install(EXPORT phylloptimTargets
        FILE phylloptimTargets.cmake
        NAMESPACE phylloptim::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phylloptim)

# The exported target carries phylloptim's own headers but NOT odelia's: odelia has no
# CMake package to depend on, so the path is a BUILD_INTERFACE only and an
# installed consumer has to supply it. phylloptimConfig.cmake asks for it by the same
# variable name and fails with the same message rather than letting the error
# surface as a missing <odelia/interpolator.hpp> from inside our headers.
configure_package_config_file(
  cmake/phylloptimConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/phylloptimConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phylloptim)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/phylloptimConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion
  ARCH_INDEPENDENT)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/phylloptimConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/phylloptimConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phylloptim)
