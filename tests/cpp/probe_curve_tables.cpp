// Two questions about the vulnerability tabulations, measured on the root and
// stem defaults.
//
// 1. The root holds two interpolants on one grid: G, whose supplied slopes ARE
//    f_r, and f_r itself. So f_r is available twice -- as the second table's
//    value and as the first table's slope. How far apart are they?
// 2. The stem holds G and its inverse. dG/dpsi is f_r in closed form, so a
//    coarse inverse table plus Newton steps converges to the same relation. What
//    knot count does that need against the 1600 the raw table needs?
#include <phylloptim/vulnerability.hpp>
#include <odelia/interpolator.hpp>
#include <cstdio>
#include <cmath>
#include <vector>

using phylloptim::cumulative_vulnerability_integral;
using phylloptim::cumulative_vulnerability_integral_at;
using phylloptim::vulnerability_curve;
using phylloptim::vulnerability_curve_slope;
using odelia::interpolator::hermite_interpolator;

int main() {
  const double root_b = 3.898245, root_c = 2.680147;
  const double stem_b = 2.0,      stem_c = 5.0;

  std::printf("== 1. f_r read as its own table against G's slope (root b=%.4f c=%.4f)\n",
              root_b, root_c);
  std::printf("%7s %13s %13s %13s\n", "knots", "own table", "G's slope", "the two apart");
  for (double res : {50.0, 100.0, 200.0, 400.0, 800.0, 1600.0}) {
    std::vector<double> x, G, fr;
    cumulative_vulnerability_integral(root_b, root_c, res, x, G, fr);
    std::vector<double> frs(x.size());
    for (std::size_t k = 0; k < x.size(); ++k)
      frs[k] = vulnerability_curve_slope(x[k], root_b, root_c);
    hermite_interpolator<double> own, integral;
    own.init(x, fr, frs);
    integral.init(x, G, fr);
    double e_own = 0, e_slope = 0, e_pair = 0;
    const int m = 20000;
    for (int i = 1; i < m; ++i) {
      const double psi = x.back() * i / m;
      const double exact = vulnerability_curve(psi, root_b, root_c);
      const double a = own.eval(psi), b = integral.slope(psi);
      e_own   = std::max(e_own,   std::abs(a - exact));
      e_slope = std::max(e_slope, std::abs(b - exact));
      e_pair  = std::max(e_pair,  std::abs(a - b));
    }
    std::printf("%7.0f %13.3e %13.3e %13.3e\n", res, e_own, e_slope, e_pair);
  }

  std::printf("\n== 2. the stem's inverse: a raw table against a coarse one polished\n");
  std::printf("   (Newton on G(psi) - target, with dG/dpsi = f_r exact)\n");
  std::printf("%7s %13s %13s %13s %9s\n",
              "knots", "raw round trip", "+1 Newton", "+2 Newton", "table KB");
  for (double res : {17.0, 25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1600.0}) {
    std::vector<double> x, G, fr;
    cumulative_vulnerability_integral(stem_b, stem_c, res, x, G, fr);
    // the inverse tabulated the way the stem holds it: psi against G, with slope
    // dpsi/dG = 1 / f_r
    std::vector<double> inv_slope(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) inv_slope[k] = 1.0 / fr[k];
    hermite_interpolator<double> inverse;
    inverse.init(G, x, inv_slope);
    double e0 = 0, e1 = 0, e2 = 0;
    const int m = 5000;
    for (int i = 1; i < m; ++i) {
      const double target = G.back() * i / m;
      double psi = inverse.eval(target);
      e0 = std::max(e0, std::abs(cumulative_vulnerability_integral_at(psi, stem_b, stem_c) - target));
      for (int step = 0; step < 2; ++step) {
        const double r = cumulative_vulnerability_integral_at(psi, stem_b, stem_c) - target;
        psi -= r / vulnerability_curve(psi, stem_b, stem_c);
        const double e = std::abs(
            cumulative_vulnerability_integral_at(psi, stem_b, stem_c) - target);
        if (step == 0) e1 = std::max(e1, e); else e2 = std::max(e2, e);
      }
    }
    std::printf("%7.0f %13.3e %13.3e %13.3e %9.1f\n", res, e0, e1, e2,
                x.size() * 9 * 2 * 8 / 1024.0);
  }
  return 0;
}
