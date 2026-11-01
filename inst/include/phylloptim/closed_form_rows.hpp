// -*-c++-*-
#ifndef PHYLLOPTIM_CLOSED_FORM_ROWS_HPP_
#define PHYLLOPTIM_CLOSED_FORM_ROWS_HPP_

#include <phylloptim/vulnerability.hpp>
#include <odelia/ode_util.hpp>

namespace phylloptim {

// Putting a closed-form row onto a table read, for the stem and root curves alike.
//
// The value is the table's at every call site, because the solve ran on the table;
// the trait rows are the curve's, in (P50, c) through rows_in_P50, because a table
// carries none. The query row is supplied where the assembly already forms it
// (dG/dpsi, from stem_curve_integral_deriv and root_vuln_integral_deriv_at) and is
// the curve's own where nothing else does (df/dpsi), so each has one value. Every
// bracket below is exactly zero at the recording point, so only the tape sees the
// rows.
//
// One place for both curves: written twice, one copy would get the chain in c and
// the other keep a partial at fixed b, a plausible wrong row nothing reports.

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

// Q(psi) = G(inf) - G(psi), with rows that keep their relative precision as Q goes
// to zero; a difference of two G rows would carry G(inf)'s rounding, all of Q.
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
