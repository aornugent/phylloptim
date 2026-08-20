// Which slope rule a tabulated vulnerability curve wants, and what the knot
// count is buying. Four fits of the same knots against the closed form:
//
//   natural   the C2 cubic through values alone -- what the root pair reads today
//   exact     C1 Hermite, values and slopes both closed form -- the stem pair
//   bessel    C1 Hermite, slopes from the parabola through three values
//   monotone  C1 Hermite, slopes limited so no span overshoots its data
//
// Reported per knot count: the worst value and slope error, split between the
// interior and the two spans at each end, where a natural spline's zero-curvature
// condition is wrong; the round trip through the inverse, which is what the knot
// count was chosen on; and the build.

#include <phylloptim/vulnerability.hpp>
#include <odelia/interpolator.hpp>
#include <odelia/interpolator.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using phylloptim::cumulative_vulnerability_integral_derivatives_at;
using phylloptim::vulnerability_psi_max;

namespace {

// The knot grid the model builds, accumulation included.
void knots(double b, double c, double resolution, std::vector<double>& x,
           std::vector<double>& g, std::vector<double>& f) {
  x = {0.0}; g = {0.0}; f = {1.0};
  const double psi_max = vulnerability_psi_max(b, c);
  const double step = psi_max / resolution;
  for (double psi = step; psi <= psi_max; psi += step) {
    const auto d = cumulative_vulnerability_integral_derivatives_at(psi, b, c);
    x.push_back(psi); g.push_back(d.value); f.push_back(d.dpsi);
  }
}

std::vector<double> secants(const std::vector<double>& x,
                            const std::vector<double>& y) {
  std::vector<double> d(x.size() - 1);
  for (size_t k = 0; k + 1 < x.size(); ++k) {
    d[k] = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);
  }
  return d;
}

// The parabola through three values, one-sided at the ends.
std::vector<double> bessel_slopes(const std::vector<double>& x,
                                  const std::vector<double>& y) {
  const std::vector<double> d = secants(x, y);
  const size_t n = x.size();
  std::vector<double> m(n);
  m[0] = d[0]; m[n - 1] = d[n - 2];
  for (size_t k = 1; k + 1 < n; ++k) {
    const double h0 = x[k] - x[k - 1], h1 = x[k + 1] - x[k];
    m[k] = (h1 * d[k - 1] + h0 * d[k]) / (h0 + h1);
  }
  return m;
}

// Bessel, then limited into the region where no span leaves its data's range.
std::vector<double> monotone_slopes(const std::vector<double>& x,
                                    const std::vector<double>& y) {
  const std::vector<double> d = secants(x, y);
  std::vector<double> m = bessel_slopes(x, y);
  for (size_t k = 0; k + 1 < x.size(); ++k) {
    if (d[k] == 0.0) { m[k] = 0.0; m[k + 1] = 0.0; continue; }
    const double a = m[k] / d[k], bb = m[k + 1] / d[k];
    if (a < 0.0) m[k] = 0.0;
    if (bb < 0.0) m[k + 1] = 0.0;
    const double s = a * a + bb * bb;
    if (s > 9.0) {
      const double tau = 3.0 / std::sqrt(s);
      m[k] = tau * a * d[k];
      m[k + 1] = tau * bb * d[k];
    }
  }
  return m;
}

struct Err { double value = 0.0, slope = 0.0; };

void worse(Err& e, double v, double s) {
  if (v > e.value) e.value = v;
  if (s > e.slope) e.slope = s;
}

}  // namespace

int main() {
  const double b = 3.898245, c = 2.680147;  // the root defaults
  const int per_span = 8;

  std::printf("Weibull vulnerability curve, b = %.6f  c = %.6f\n", b, c);
  std::printf("value errors relative to G(psi_max); slope errors absolute "
              "(f_r is in [0.01, 1])\n\n");
  std::printf("%6s  %-9s  %10s %10s  %10s %10s  %11s %8s\n",
              "knots", "rule", "val int", "val end", "slp int", "slp end",
              "round trip", "build us");

  for (double res : {50.0, 100.0, 200.0, 400.0, 800.0, 1600.0}) {
    std::vector<double> x, g, f;
    knots(b, c, res, x, g, f);
    const size_t n = x.size();
    const double range = g.back();

    // The four fits, and the build each costs.
    odelia::spline::Spline natural;
    auto t0 = std::chrono::steady_clock::now();
    natural.set_points(x, g);
    auto t1 = std::chrono::steady_clock::now();
    const double us_natural =
        std::chrono::duration<double, std::micro>(t1 - t0).count();

    std::vector<double> m_bessel, m_mono;
    t0 = std::chrono::steady_clock::now();
    m_bessel = bessel_slopes(x, g);
    t1 = std::chrono::steady_clock::now();
    const double us_bessel =
        std::chrono::duration<double, std::micro>(t1 - t0).count();
    m_mono = monotone_slopes(x, g);

    odelia::interpolator::hermite_interpolator<double> exact, bessel, mono;
    t0 = std::chrono::steady_clock::now();
    exact.init(x, g, f);
    t1 = std::chrono::steady_clock::now();
    const double us_exact =
        std::chrono::duration<double, std::micro>(t1 - t0).count();
    bessel.init(x, g, m_bessel);
    mono.init(x, g, m_mono);

    // The inverse each rule implies, on the same relation read the other way.
    std::vector<double> inv_slope_exact(n), inv_slope_bessel(n);
    for (size_t k = 0; k < n; ++k) {
      inv_slope_exact[k] = 1.0 / f[k];
      inv_slope_bessel[k] = 1.0 / m_bessel[k];
    }
    odelia::interpolator::hermite_interpolator<double> inv_exact, inv_bessel;
    inv_exact.init(g, x, inv_slope_exact);
    inv_bessel.init(g, x, inv_slope_bessel);
    odelia::spline::Spline inv_natural;
    inv_natural.set_points(g, x);

    Err e_nat, e_exa, e_bes, e_mon;
    Err d_nat, d_exa, d_bes, d_mon;  // the end spans
    double rt_nat = 0.0, rt_exa = 0.0, rt_bes = 0.0;

    for (size_t k = 0; k + 1 < n; ++k) {
      const bool end = (k < 2) || (k + 3 > n);
      const double h = x[k + 1] - x[k];
      for (int j = 1; j < per_span; ++j) {
        const double u = x[k] + h * j / per_span;
        const auto ref = cumulative_vulnerability_integral_derivatives_at(u, b, c);
        Err& ev_n = end ? d_nat : e_nat;
        Err& ev_e = end ? d_exa : e_exa;
        Err& ev_b = end ? d_bes : e_bes;
        Err& ev_m = end ? d_mon : e_mon;
        worse(ev_n, std::fabs(natural(u) - ref.value) / range,
              std::fabs(natural.deriv(u) - ref.dpsi));
        worse(ev_e, std::fabs(exact.eval(u) - ref.value) / range,
              std::fabs(exact.slope(u) - ref.dpsi));
        worse(ev_b, std::fabs(bessel.eval(u) - ref.value) / range,
              std::fabs(bessel.slope(u) - ref.dpsi));
        worse(ev_m, std::fabs(mono.eval(u) - ref.value) / range,
              std::fabs(mono.slope(u) - ref.dpsi));

        // G forward then back: the two tables' agreement that they are one
        // relation, which is what the knot count was set on.
        const double back_n = inv_natural(natural(u));
        const double back_e = inv_exact.eval(exact.eval(u));
        const double back_b = inv_bessel.eval(bessel.eval(u));
        rt_nat = std::max(rt_nat, std::fabs(back_n - u));
        rt_exa = std::max(rt_exa, std::fabs(back_e - u));
        rt_bes = std::max(rt_bes, std::fabs(back_b - u));
      }
    }

    std::printf("%6.0f  %-9s  %10.2e %10.2e  %10.2e %10.2e  %11.2e %8.1f\n",
                res, "natural", e_nat.value, d_nat.value, e_nat.slope,
                d_nat.slope, rt_nat, us_natural);
    std::printf("%6s  %-9s  %10.2e %10.2e  %10.2e %10.2e  %11.2e %8.1f\n",
                "", "exact", e_exa.value, d_exa.value, e_exa.slope, d_exa.slope,
                rt_exa, us_exact);
    std::printf("%6s  %-9s  %10.2e %10.2e  %10.2e %10.2e  %11.2e %8.1f\n",
                "", "bessel", e_bes.value, d_bes.value, e_bes.slope, d_bes.slope,
                rt_bes, us_bessel);
    std::printf("%6s  %-9s  %10.2e %10.2e  %10.2e %10.2e  %11s %8s\n\n",
                "", "monotone", e_mon.value, d_mon.value, e_mon.slope,
                d_mon.slope, "-", "-");
  }
  return 0;
}
