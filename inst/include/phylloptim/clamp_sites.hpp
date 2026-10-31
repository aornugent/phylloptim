// -*-c++-*-
#ifndef PHYLLOPTIM_CLAMP_SITES_HPP_
#define PHYLLOPTIM_CLAMP_SITES_HPP_

#include <cstddef>
#include <memory>
#include <vector>

namespace phylloptim {

// Where this model replaces a value by a bound, counted. One site is left: the
// stomatal corner E_s >= 0, the one inequality the solve keeps. Every other
// bound a value used to be bent to is now a refusal at the entry point or a
// smooth form, so a solved path notes nothing else.
//
// Counted rather than refused because the corner is part of the model: the
// value function is C1 there and the argmax has one-sided derivatives, as at
// any active bound. The count says how often a run sat on it.
//
// The list is this package's own. A consumer that keeps its own tally folds these
// in rather than sharing an enum, because a header-only model cannot depend on the
// consumer that includes it.
enum clamp_site {
  // The stomatal flux read as zero where the hydraulic flux falls short of the
  // residual flux. Measured: the operating-point evaluator reaches it only by
  // rounding at the wet bound; the derivative evaluations taken just wet of
  // that bound, to read the corner's sign, reach it by up to 0.7% of r*.
  CLAMP_STOMATAL_FLUX = 0,
  CLAMP_SITE_COUNT
};

inline const char* clamp_site_name(int site) {
  switch (site) {
  case CLAMP_STOMATAL_FLUX:          return "stomatal_flux";
  }
  return "unknown";
}

// The tally, shared so that a copy still reports. A consumer holds this model by
// value and copies it per unit -- plant copies its whole strategy per cohort and
// discards it -- so a plain member would take every count down with the copy. The
// counts are behind a pointer for that reason and for no other.
struct clamp_counter {
  std::shared_ptr<std::vector<std::size_t>> counts =
    std::make_shared<std::vector<std::size_t>>(CLAMP_SITE_COUNT, 0);

  void note(int site) const { ++(*counts)[site]; }
  std::size_t at(int site) const { return (*counts)[site]; }
  void clear() const { counts->assign(CLAMP_SITE_COUNT, 0); }
};

}  // namespace phylloptim

#endif
