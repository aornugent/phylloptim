// -*-c++-*-
#ifndef PHYLLOPTIM_CLAMP_SITES_HPP_
#define PHYLLOPTIM_CLAMP_SITES_HPP_

#include <cstddef>
#include <memory>
#include <vector>

namespace phylloptim {

// Where this model replaces a value by a bound. Each one holds the quantity it
// masks constant, so whatever that quantity fed stops responding to it -- and a
// row severed that way and a row that is honestly zero are the same number.
//
// Every site here already pairs its clamp with a matching derivative, so what is
// counted is not a defect: it is the DISTANCE from one. Two of these sit within a
// few hundredths of a soil moisture unit of the wilt point, and until they were
// counted nobody could say whether a run reached them.
//
// The list is this package's own. A consumer that keeps its own tally folds these
// in rather than sharing an enum, because a header-only model cannot depend on the
// consumer that includes it.
enum clamp_site {
  // The cumulative-vulnerability integral, capped at its closed-form limit
  // G(inf) = (b/c)*Gamma(1/c). Binds past 7.3132 MPa at the root defaults, where
  // the end-knot polynomial would otherwise accumulate for ever.
  CLAMP_ROOT_VULN_INTEGRAL_CAP = 0,
  // The conductivity lookup's argument, clamped into the knot domain. The upper
  // end is the curve's last knot (6.8229 MPa at the root defaults), the lower is
  // zero.
  CLAMP_ROOT_VULN_ARGUMENT,
  // The energy-balance leaf temperature, held inside a physical range so an
  // extreme non-equilibrium transpiration cannot drive the Arrhenius block
  // non-finite. Reached only with use_energy_balance on.
  CLAMP_LEAF_TEMPERATURE,
  // A collar potential projected into the feasible interval. For a partial taken
  // at a held collar the unclamped alternative is profit_at_fixed_collar, which
  // refuses instead; this one is for a tracked collar being nudged.
  CLAMP_COLLAR_POTENTIAL,
  CLAMP_SITE_COUNT
};

inline const char* clamp_site_name(int site) {
  switch (site) {
  case CLAMP_ROOT_VULN_INTEGRAL_CAP: return "root_vuln_integral_cap";
  case CLAMP_ROOT_VULN_ARGUMENT:     return "root_vuln_argument";
  case CLAMP_LEAF_TEMPERATURE:       return "leaf_temperature";
  case CLAMP_COLLAR_POTENTIAL:       return "collar_potential";
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
