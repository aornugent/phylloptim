// PROBE: is there a form of the layer mean that beats the three thresholds?
//
// The mean of the root vulnerability curve over a layer's suction interval is
// (G(hi) - G(lo))/span -- a difference of two reads of a TABULATED G, so its
// relative error is the table's divided by the span, and every derivative of it
// divides by the span again. That is why the model carries THREE measured
// thresholds and not one, and below each of them switches to a midpoint
// asymptotic. The bad band is where neither form is good: at the crossover both
// are at their worst by construction.
//
// Four arms, one reference:
//
//   REF  Gauss-Legendre at high order on the CLOSED-FORM integrand. No table and
//        no subtraction, so it is trustworthy exactly where the others are not.
//   A    what the model does now: the thresholds and the asymptotic
//   B    the divided difference alone, no threshold -- what A exists to avoid
//   C    Gauss-Legendre at LOW order on the closed form, with every bound
//        derivative from differentiating the finite sum
//
// C is the interesting one. mean(lo,hi) = 0.5 * sum w_i f(x_i) with
// x_i = m + (s/2) t_i, so d(mean)/dhi = 0.5 * sum w_i f'(x_i) (1+t_i)/2 and the
// second derivatives are the same sum one order up. NOTHING DIVIDES BY THE SPAN,
// at any order, so there is no cancellation to threshold around.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include \
//        ODELIA_SRC=../../../odelia/src probe_layer_mean && ./probe_layer_mean

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace pl = phylloptim;
using clock_type = std::chrono::steady_clock;

namespace {

const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_up(pl::Leaf& l) {
  double theta[pl::gradient::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::gradient::par_R_d_25] = 1.44;
  l.set_traits(theta);
  std::vector<double> root{1.0 / kAreaLeaf}, psi_soil{2.0}, depth{1.0};
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

// Gauss-Legendre nodes and weights on [-1,1], by Newton on the Legendre
// polynomial. Generated rather than tabulated so the order is a parameter and a
// transcription cannot be wrong.
struct Rule {
  std::vector<double> t, w;
};

Rule gauss_legendre(int n) {
  Rule r;
  r.t.assign(std::size_t(n), 0.0);
  r.w.assign(std::size_t(n), 0.0);
  for (int i = 0; i < n; ++i) {
    double x = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0, p1 = 0.0;
      for (int j = 0; j < n; ++j) {
        const double p2 = p1;
        p1 = p0;
        p0 = ((2.0 * j + 1.0) * x * p1 - j * p2) / (j + 1);
      }
      const double dp = n * (x * p0 - p1) / (x * x - 1.0);
      const double dx = -p0 / dp;
      x += dx;
      if (std::abs(dx) < 1e-16) break;
    }
    double p0 = 1.0, p1 = 0.0;
    for (int j = 0; j < n; ++j) {
      const double p2 = p1;
      p1 = p0;
      p0 = ((2.0 * j + 1.0) * x * p1 - j * p2) / (j + 1);
    }
    const double dp = n * (x * p0 - p1) / (x * x - 1.0);
    r.t[std::size_t(i)] = x;
    r.w[std::size_t(i)] = 2.0 / ((1.0 - x * x) * dp * dp);
  }
  return r;
}

double b_root = 0.0, c_root = 0.0;

double f_at(double psi) {
  return pl::vulnerability_curve_at<double>(psi, b_root, c_root);
}
double df_at(double psi) {
  return pl::vulnerability_curve_slope_at<double>(psi, b_root, c_root);
}
// f'' from the slope's own closed form by one tangent -- the same route
// curve_slope2_at takes, and the only second order anywhere in this probe.
double d2f_at(double psi) {
  pl::tangent p = psi;
  pl::seed_direction(p, 1.0);
  return pl::derivative_along(pl::vulnerability_curve_slope_at<pl::tangent>(
      p, pl::tangent(b_root), pl::tangent(c_root)));
}

// The quadrature form and its three bound derivatives. Every one is a weighted
// sum over the same nodes; none divides by the span.
struct Quad {
  double mean, dhi, dhi2, dmixed;
};

Quad quad_at(const Rule& r, double lo, double hi) {
  const double m = 0.5 * (lo + hi), s = hi - lo;
  Quad q{0.0, 0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < r.t.size(); ++i) {
    const double x = m + 0.5 * s * r.t[i];
    const double a_hi = 0.5 * (1.0 + r.t[i]);   // dx_i/dhi
    const double a_lo = 0.5 * (1.0 - r.t[i]);   // dx_i/dlo
    q.mean += r.w[i] * f_at(x);
    q.dhi += r.w[i] * df_at(x) * a_hi;
    q.dhi2 += r.w[i] * d2f_at(x) * a_hi * a_hi;
    q.dmixed += r.w[i] * d2f_at(x) * a_hi * a_lo;
  }
  q.mean *= 0.5;
  q.dhi *= 0.5;
  q.dhi2 *= 0.5;
  q.dmixed *= 0.5;
  return q;
}

double rel(double got, double want) {
  return std::abs(got - want) / (std::abs(want) + 1e-300);
}

const char* mark(double e) {
  return e < 1e-12 ? "" : (e < 1e-6 ? "  ." : "  <-");
}

}  // namespace

int main() {
  pl::Leaf l;
  set_up(l);
  l.find_root_collar_psi();
  b_root = l.roots_.root_b;
  c_root = l.roots_.root_c;

  const Rule ref = gauss_legendre(24);
  const Rule c3 = gauss_legendre(3);
  const Rule c5 = gauss_legendre(5);
  const Rule c7 = gauss_legendre(7);

  std::printf("the layer mean: the thresholds against forms that need none\n");
  std::printf("  root curve b=%.6g c=%.6g; reference is Gauss-Legendre 24 on the"
              " closed form\n", b_root, c_root);
  std::printf("  thresholds now: mean 1e-5, d2 5e-4, mixed 5e-3\n");

  const double centres[3] = {1.0, 2.0, 4.0};
  for (double centre : centres) {
    std::printf("\n=== centre psi = %.1f ===\n", centre);
    std::printf("  %-9s | %-22s | %-22s | %-22s\n", "span",
                "MEAN  A(model)  C(GL5)", "d/dbound  A       C(GL5)",
                "d2/dbound2  A     C(GL5)");
    for (int e = 0; e >= -12; --e) {
      const double s = std::pow(10.0, double(e));
      const double lo = centre - 0.5 * s, hi = centre + 0.5 * s;
      if (lo <= 0.0) continue;

      const Quad R = quad_at(ref, lo, hi);
      const Quad C = quad_at(c5, lo, hi);

      const double a_mean = l.roots_.layer_mean(lo, hi);
      const double a_dhi = l.roots_.layer_mean_dbound(lo, hi, true, a_mean);
      const double a_d2 = l.roots_.layer_mean_dbound2(lo, hi, true, a_mean);

      const double em_a = rel(a_mean, R.mean), em_c = rel(C.mean, R.mean);
      const double e1_a = rel(a_dhi, R.dhi), e1_c = rel(C.dhi, R.dhi);
      const double e2_a = rel(a_d2, R.dhi2), e2_c = rel(C.dhi2, R.dhi2);

      std::printf("  %-9.0e | %8.2e %8.2e%-3s | %8.2e %8.2e%-3s | %8.2e %8.2e%-3s\n",
                  s, em_a, em_c, mark(em_a), e1_a, e1_c, mark(e1_a), e2_a, e2_c,
                  mark(e2_a));
    }
  }

  // How few nodes C needs across the whole range, which is what decides its cost.
  std::printf("\n=== worst relative error over spans 1e0..1e-12, all three centres ===\n");
  struct Arm { const char* name; const Rule* r; };
  const Arm arms[3] = {{"C, Gauss 3", &c3}, {"C, Gauss 5", &c5}, {"C, Gauss 7", &c7}};
  for (const Arm& a : arms) {
    double wm = 0.0, w1 = 0.0, w2 = 0.0, wx = 0.0;
    for (double centre : centres) {
      for (int e = 0; e >= -12; --e) {
        const double s = std::pow(10.0, double(e));
        const double lo = centre - 0.5 * s, hi = centre + 0.5 * s;
        if (lo <= 0.0) continue;
        const Quad R = quad_at(ref, lo, hi), Q = quad_at(*a.r, lo, hi);
        wm = std::max(wm, rel(Q.mean, R.mean));
        w1 = std::max(w1, rel(Q.dhi, R.dhi));
        w2 = std::max(w2, rel(Q.dhi2, R.dhi2));
        wx = std::max(wx, rel(Q.dmixed, R.dmixed));
      }
    }
    std::printf("  %-12s  mean %.2e   d/db %.2e   d2/db2 %.2e   mixed %.2e\n",
                a.name, wm, w1, w2, wx);
  }

  // The model's own worst, over the same grid, for the same four quantities.
  {
    double wm = 0.0, w1 = 0.0, w2 = 0.0, wx = 0.0;
    for (double centre : centres) {
      for (int e = 0; e >= -12; --e) {
        const double s = std::pow(10.0, double(e));
        const double lo = centre - 0.5 * s, hi = centre + 0.5 * s;
        if (lo <= 0.0) continue;
        const Quad R = quad_at(ref, lo, hi);
        const double m = l.roots_.layer_mean(lo, hi);
        wm = std::max(wm, rel(m, R.mean));
        w1 = std::max(w1, rel(l.roots_.layer_mean_dbound(lo, hi, true, m), R.dhi));
        w2 = std::max(w2, rel(l.roots_.layer_mean_dbound2(lo, hi, true, m), R.dhi2));
        wx = std::max(wx, rel(l.roots_.layer_mean_dbound_mixed(lo, hi, m), R.dmixed));
      }
    }
    std::printf("  %-12s  mean %.2e   d/db %.2e   d2/db2 %.2e   mixed %.2e\n",
                "A, the model", wm, w1, w2, wx);
  }

  // And the price, LIKE FOR LIKE. The table read is one quantity, so a quadrature
  // timed while it also forms three derivatives -- one of them through a tangent --
  // is not the comparison anyone would make.
  {
    const long reps = 300000;
    const double lo = 1.999, hi = 2.001;
    double sink = 0.0;

    auto time_it = [&](auto&& fn) -> double {
      const clock_type::time_point t0 = clock_type::now();
      for (long r = 0; r < reps; ++r) sink += fn(lo + 1e-12 * double(r), hi);
      return std::chrono::duration<double, std::micro>(clock_type::now() - t0)
                 .count() /
             double(reps);
    };

    // mean only, both ways
    auto mean_quad = [](const Rule& r, double a, double b) -> double {
      const double m = 0.5 * (a + b), s = b - a;
      double acc = 0.0;
      for (std::size_t i = 0; i < r.t.size(); ++i) {
        acc += r.w[i] * f_at(m + 0.5 * s * r.t[i]);
      }
      return 0.5 * acc;
    };
    const double us_a = time_it([&](double a, double b) {
      return l.roots_.layer_mean(a, b);
    });
    const double us_c5 =
        time_it([&](double a, double b) { return mean_quad(c5, a, b); });
    const double us_c7 =
        time_it([&](double a, double b) { return mean_quad(c7, a, b); });

    // the first bound derivative, both ways, each given what it needs
    const double us_a1 = time_it([&](double a, double b) {
      return l.roots_.layer_mean_dbound(a, b, true, l.roots_.layer_mean(a, b));
    });
    auto dbound_quad = [](const Rule& r, double a, double b) -> double {
      const double m = 0.5 * (a + b), s = b - a;
      double acc = 0.0;
      for (std::size_t i = 0; i < r.t.size(); ++i) {
        acc += r.w[i] * df_at(m + 0.5 * s * r.t[i]) * 0.5 * (1.0 + r.t[i]);
      }
      return 0.5 * acc;
    };
    const double us_c71 =
        time_it([&](double a, double b) { return dbound_quad(c7, a, b); });

    std::printf("\n  cost, like for like (us per call, %ld reps)\n", reps);
    std::printf("    mean        A(table) %.4f   C(GL5) %.4f (%.1fx)   "
                "C(GL7) %.4f (%.1fx)\n",
                us_a, us_c5, us_c5 / us_a, us_c7, us_c7 / us_a);
    std::printf("    d/dbound    A        %.4f   C(GL7) %.4f (%.1fx)\n", us_a1,
                us_c71, us_c71 / us_a1);
    std::printf("    (checksum %.6g)\n", sink);
  }

  return 0;
}
