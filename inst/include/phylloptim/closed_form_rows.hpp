// -*-c++-*-
#ifndef PHYLLOPTIM_CLOSED_FORM_ROWS_HPP_
#define PHYLLOPTIM_CLOSED_FORM_ROWS_HPP_

#include <phylloptim/vulnerability.hpp>
#include <odelia/ode_util.hpp>

namespace phylloptim {

// Putting a closed-form row onto a table read.
//
// THE VALUE IS THE TABLE'S, in both functions and at every call site, because the
// solve ran on the table: a value from the curve would place the operating point
// somewhere else. THE TRAIT ROWS ARE THE CURVE'S, in both, because a table carries
// none.
//
// THE QUERY ROW IS WHICHEVER THE REST OF THE ASSEMBLY ALREADY FORMS, and it is the
// only one that differs between the two below. Where the model forms that
// derivative from a table, a closed-form row is a SECOND value of it and the two
// disagree about where the operating point is stationary; where the model forms it
// nowhere, a table row is a row of the fit rather than of the curve. Which case a
// function is in is a grep and not a judgement:
//
//   G  dG/dpsi IS formed -- stem_curve_integral_deriv and
//      root_vuln_integral_deriv_at, both read by the physics -- so
//      closed_form_integral takes the table's slope as an argument.
//   f  df/dpsi is formed nowhere. root_vuln_from_psi is built WITH slopes and
//      .slope() is never called on it, and odelia's interpolator reads no second
//      derivative at all, so closed_form_curve uses the curve's own.
//
// Every bracket below is exactly zero in VALUE at the recording point, so the
// number is untouched and only the tape sees the rows.
//
// ⚠️ ONE PLACE, FOR TWO CURVES. The stem and the root do this identically, and
// the rows go through rows_in_P50 -- so written twice, one copy gets the chain in
// c and the other keeps a partial at fixed b. That is a finite, plausible, wrong
// row that nothing reports, which is exactly the class of defect this file exists
// to make unwritable.

// G(psi): the cumulative vulnerability integral.
//
// ⚠️ THE QUERY SLOPE IS THE CALLER'S, AND IT MUST BE THE TABLE'S -- this is the
// case above where the model forms dG/dpsi itself. Measured with the closed form
// here instead: the residual lift's dci/dcollar read 0.6913 against 0.6382 from
// the explicit slope beside it, and dM/dcollar came out at an eighth of its true
// size. Both are finite and plausible and they differ by the fit error.
template <class S>
inline S closed_form_integral(double table_value, double table_slope, const S& psi,
                        const S& P50, const S& c) {
  using odelia::util::to_passive;
  const double at = to_passive(psi);
  const double P50_0 = to_passive(P50);
  const double c_0 = to_passive(c);
  const double b = weibull_b_from_P50(P50_0, c_0);
  const VulnerabilityIntegralDerivatives d =
      cumulative_vulnerability_integral_derivatives_at(at, b, c_0);
  const TraitRows rows = rows_in_P50(d.db, d.dc, b, c_0);
  return S(table_value) + S(table_slope) * (psi - S(at)) +
         S(rows.dP50) * (P50 - S(P50_0)) + S(rows.dc) * (c - S(c_0));
}

// f(psi): the surviving conductivity, the curve itself. No table-slope argument,
// because nothing else in the assembly forms df/dpsi to disagree with.
template <class S>
inline S closed_form_curve(double table_value, const S& psi, const S& P50,
                     const S& c) {
  using odelia::util::to_passive;
  const double at = to_passive(psi);
  const double P50_0 = to_passive(P50);
  const double c_0 = to_passive(c);
  const double b = weibull_b_from_P50(P50_0, c_0);
  const VulnerabilityDerivatives d = vulnerability_derivatives_at(at, b, c_0);
  const TraitRows rows = rows_in_P50(d.db, d.dc, b, c_0);
  return S(table_value) + S(d.dpsi) * (psi - S(at)) +
         S(rows.dP50) * (P50 - S(P50_0)) + S(rows.dc) * (c - S(c_0));
}

}

#endif
