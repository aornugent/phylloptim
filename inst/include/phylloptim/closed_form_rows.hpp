// -*-c++-*-
#ifndef PHYLLOPTIM_CLOSED_FORM_ROWS_HPP_
#define PHYLLOPTIM_CLOSED_FORM_ROWS_HPP_

#include <phylloptim/vulnerability.hpp>
#include <odelia/ode_util.hpp>

namespace phylloptim {

// Putting a closed-form row onto a table read.
//
// THE VALUE IS THE TABLE'S, in both functions and at every call site, because the
// solve ran on the table: a value from elsewhere would place the operating point
// somewhere else, if only by the table's ~1e-14. THE TRAIT ROWS ARE THE CURVE'S,
// in both, because a table carries none.
//
// THE QUERY ROW IS THE ONE THE REST OF THE ASSEMBLY ALREADY FORMS: dG/dpsi is
// formed by stem_curve_integral_deriv and root_vuln_integral_deriv_at, so
// integral_with_closed_form_rows takes it as an argument rather than forming a
// second value of it; df/dpsi is formed nowhere else, so
// curve_with_closed_form_rows uses the curve's own. Both are closed forms now, so the two could not disagree in any
// case, but the argument keeps there being one value of each.
//
// Every bracket below is exactly zero in VALUE at the recording point, so the
// number is untouched and only the tape sees the rows.
//
// ⚠️ ONE PLACE, FOR TWO CURVES. The stem and the root do this identically, and
// the rows go through rows_in_P50 -- so written twice, one copy gets the chain in
// c and the other keeps a partial at fixed b. That is a finite, plausible, wrong
// row that nothing reports, which is exactly the class of defect this file exists
// to make unwritable.

// A table read with the rows `derivatives_at` gives at its point, in (P50, c).
template <class S, class Derivatives>
inline S table_read_with_closed_form_rows(double table_value,
                                          double table_slope, const S& psi,
                                          const S& P50, const S& c,
                                          Derivatives&& derivatives_at) {
  using odelia::util::to_passive;
  const double at = to_passive(psi);
  const double P50_0 = to_passive(P50);
  const double c_0 = to_passive(c);
  const double b = weibull_b_from_P50(P50_0, c_0);
  const VulnerabilityIntegralDerivatives d = derivatives_at(at, b, c_0);
  const TraitRows rows = rows_in_P50(d.db, d.dc, b, c_0);
  return S(table_value) + S(table_slope) * (psi - S(at)) +
         S(rows.dP50) * (P50 - S(P50_0)) + S(rows.dc) * (c - S(c_0));
}

// G(psi): the cumulative vulnerability integral. Past the end of the curve the
// rows are G(inf)'s, which cumulative_vulnerability_integral_derivatives_at
// returns there.
template <class S>
inline S integral_with_closed_form_rows(double table_value, double table_slope,
                                        const S& psi, const S& P50,
                                        const S& c) {
  return table_read_with_closed_form_rows<S>(
      table_value, table_slope, psi, P50, c,
      cumulative_vulnerability_integral_derivatives_at);
}

// Q(psi) = G(inf) - G(psi): what the curve has left past psi, with rows that
// keep their relative precision as Q goes to zero
// (cumulative_vulnerability_complement_derivatives_at). A difference of two
// G rows there would carry G(inf)'s rounding, which is all of Q's size.
template <class S>
inline S complement_with_closed_form_rows(double table_value,
                                          double table_slope, const S& psi,
                                          const S& P50, const S& c) {
  return table_read_with_closed_form_rows<S>(
      table_value, table_slope, psi, P50, c,
      cumulative_vulnerability_complement_derivatives_at);
}

// f(psi): the surviving conductivity, the curve itself.
template <class S>
inline S curve_with_closed_form_rows(double value, const S& psi, const S& P50,
                                     const S& c) {
  using odelia::util::to_passive;
  const double at = to_passive(psi);
  const double P50_0 = to_passive(P50);
  const double c_0 = to_passive(c);
  const double b = weibull_b_from_P50(P50_0, c_0);
  const VulnerabilityDerivatives d = vulnerability_derivatives_at(at, b, c_0);
  const TraitRows rows = rows_in_P50(d.db, d.dc, b, c_0);
  return S(value) + S(d.dpsi) * (psi - S(at)) +
         S(rows.dP50) * (P50 - S(P50_0)) + S(rows.dc) * (c - S(c_0));
}

}

#endif
