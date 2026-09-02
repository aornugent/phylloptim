// -*-c++-*-
#ifndef PHYLLOPTIM_GRAFT_HPP_
#define PHYLLOPTIM_GRAFT_HPP_

#include <phylloptim/vulnerability.hpp>
#include <odelia/ode_util.hpp>

namespace phylloptim {

// Putting a closed-form row onto a table read.
//
// Both vulnerability curves are evaluated from a spline and differentiated from
// the closed form, and the split is deliberate: THE VALUE IS THE TABLE'S because
// the solve ran on the table, so a value from the curve would place the operating
// point somewhere else -- while the table's own slope is a fit of the curve and
// differs from it, so a derivative read off the table is a derivative of the fit.
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
// ⚠️ THE QUERY SLOPE IS THE CALLER'S, AND IT MUST BE THE TABLE'S. dG/dpsi is a
// quantity the model itself computes -- from the table, by
// stem_curve_integral_deriv and root_vuln_integral_deriv_at -- so passing the
// closed form here puts a SECOND value of one derivative into the assembly.
// Both are finite and plausible, they differ by the fit error, and the two then
// disagree about where the operating point is stationary: measured, the residual
// lift's dci/dcollar read 0.6913 against 0.6382 from the explicit slope beside
// it, and dM/dcollar came out at an eighth of its true size.
//
// The TRAIT rows are still the curve's, because the table has none at all. That
// is the whole split: the table knows the query direction, only the closed form
// knows the others.
template <class S>
inline S graft_integral(double table_value, double table_slope, const S& psi,
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

// f(psi): the surviving conductivity, the curve itself.
template <class S>
inline S graft_curve(double table_value, const S& psi, const S& P50,
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
