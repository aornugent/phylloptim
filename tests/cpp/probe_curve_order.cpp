// G is analytic and its first TWO derivatives are closed forms: dG/dpsi is f_r
// and d2G/dpsi2 is f_r's own slope. So the question is what order of read the
// table should carry, not how many knots a C1 read needs.
#include <phylloptim/vulnerability.hpp>
#include <cstdio>
#include <cmath>
#include <vector>
using namespace phylloptim;

// value and slope at both ends -> the unique cubic; plus curvature -> the unique
// quintic. Written out rather than fitted, so the comparison is of orders only.
static double quintic(double t, double h, double y0, double y1, double m0,
                      double m1, double c0, double c1) {
  const double u = t / h;
  const double a0 = y0, a1 = m0 * h, a2 = 0.5 * c0 * h * h;
  const double d = y1 - a0 - a1 - a2, e = m1 * h - a1 - 2 * a2,
               f = c1 * h * h - 2 * a2;
  const double a3 = 10 * d - 4 * e + 0.5 * f;
  const double a4 = -15 * d + 7 * e - f;
  const double a5 = 6 * d - 3 * e + 0.5 * f;
  return a0 + u * (a1 + u * (a2 + u * (a3 + u * (a4 + u * a5))));
}
static double cubic(double t, double h, double y0, double y1, double m0, double m1) {
  const double d = y1 - y0;
  const double c2 = (3 * d - h * (2 * m0 + m1)) / (h * h);
  const double c3 = (-2 * d + h * (m0 + m1)) / (h * h * h);
  return y0 + t * (m0 + t * (c2 + t * c3));
}
// The rows scale as h^3 in the knot spacing where the read's VALUE is h^4, so what
// they follow is the slope. Both slopes here are derivatives of the polynomial the
// value above comes from.
static double cubic_slope(double t, double h, double y0, double y1, double m0, double m1) {
  const double d = y1 - y0;
  const double c2 = (3 * d - h * (2 * m0 + m1)) / (h * h);
  const double c3 = (-2 * d + h * (m0 + m1)) / (h * h * h);
  return m0 + t * (2 * c2 + t * 3 * c3);
}
static double quintic_slope(double t, double h, double y0, double y1, double m0,
                            double m1, double c0, double c1) {
  const double u = t / h;
  const double a1 = m0 * h, a2 = 0.5 * c0 * h * h;
  const double d = y1 - y0 - a1 - a2, e = m1 * h - a1 - 2 * a2,
               f = c1 * h * h - 2 * a2;
  const double a3 = 10 * d - 4 * e + 0.5 * f;
  const double a4 = -15 * d + 7 * e - f;
  const double a5 = 6 * d - 3 * e + 0.5 * f;
  return (a1 + u * (2 * a2 + u * (3 * a3 + u * (4 * a4 + u * 5 * a5)))) / h;
}

int main() {
  const double b = 2.0, c = 5.0;
  std::printf("%7s %12s %12s %12s %12s %9s %9s\n", "knots", "cubic value",
              "quint value", "cubic slope", "quint slope", "cubic KB", "quint KB");
  for (double res : {25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1600.0}) {
    std::vector<double> x, G, fr;
    cumulative_vulnerability_integral(b, c, res, x, G, fr);
    const std::size_t n = x.size();
    std::vector<double> frs(n);
    for (std::size_t k = 0; k < n; ++k) frs[k] = vulnerability_curve_slope(x[k], b, c);
    double e3 = 0, e5 = 0, s3 = 0, s5 = 0;
    for (std::size_t k = 0; k + 1 < n; ++k) {
      const double h = x[k + 1] - x[k];
      for (int i = 1; i < 8; ++i) {
        const double t = h * i / 8.0;
        const double z = x[k] + t;
        const double exact = cumulative_vulnerability_integral_at(z, b, c);
        const double exact_s = vulnerability_curve(z, b, c);
        e3 = std::max(e3, std::abs(cubic(t, h, G[k], G[k+1], fr[k], fr[k+1]) - exact));
        e5 = std::max(e5, std::abs(quintic(t, h, G[k], G[k+1], fr[k], fr[k+1],
                                           frs[k], frs[k+1]) - exact));
        s3 = std::max(s3, std::abs(cubic_slope(t, h, G[k], G[k+1], fr[k], fr[k+1]) - exact_s));
        s5 = std::max(s5, std::abs(quintic_slope(t, h, G[k], G[k+1], fr[k], fr[k+1],
                                                 frs[k], frs[k+1]) - exact_s));
      }
    }
    std::printf("%7.0f %12.3e %12.3e %12.3e %12.3e %9.1f %9.1f\n", res, e3, e5, s3, s5,
                n * 2 * 8 / 1024.0, n * 3 * 8 / 1024.0);
  }
  return 0;
}
