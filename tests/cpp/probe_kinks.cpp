// PROBE (not a test, not shipped): the three coincidences the supply's derivative
// kernels refuse at, how often each is reached, and what the limit is at each.
//
// Each kernel opens with the same test and returns a whole vector of NaN if any
// layer satisfies it:
//
//   |T_collar - psi_soil[i]| < 1e-8      the collar meets a layer's potential
//   |(T_collar - psi_soil[i]) - grav| < 1e-8   the layer's flux is gravity-balanced
//   |T_collar| < 1e-8                    the collar sits at atmospheric
//
// The claim to test is that none of the three needs to be a refusal:
//
//   * the gravity balance leaves span, the integral and the resistance untouched
//     and only makes the NUMERATOR zero, so the general expression is already
//     valid there -- nothing to derive;
//   * the collar at atmospheric has the same one-sided slope from both sides,
//     because the integrand is 1 below the surface and f_r(0) = 1 above it;
//   * equal potentials is the only one that needs a limit, and it is the
//     reciprocal of G's DIVIDED DIFFERENCE over the interval, whose expansion is
//     elementary.
//
// Writing D(s,p) for the mean of the integrand over [p, s] -- which is what the
// value branch already calls the layer's mean conductivity -- the whole family is
// Q = 1/D with
//
//   D = int_0^1 f(p + t(s-p)) dt,  D -> f,  D_s -> f'/2,  D_p -> f'/2,
//   D_ss -> f''/3,  D_sp -> f''/6,  D_pp -> f''/3.
//
//   make CXX=g++ probe_kinks && stdbuf -oL ./probe_kinks

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

const double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,      157.44,   0.30,
                         0.7,      0.99,     7.5,      1.44};

fixture::Physiology drivers(double psi_soil, double ppfd, double vpd,
                            int layers) {
  fixture::Physiology d;
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers / 0.05;
  }
  d.root_network = fixture::root_network(root, depth);
  d.PPFD = ppfd;
  d.psi_soil = ps;
  d.soil_depth = depth;
  d.kmax = 1.0 * 0.000157 / 5.0;
  d.atm_vpd = vpd;
  d.ca = 40.0;
  d.leaf_temp = 25.0;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

// --- A. how often each coincidence is reached -----------------------------------

void incidence() {
  std::printf("\n=== A. WHICH COINCIDENCE IS REACHED, AND WHERE ===\n");
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0};
  const double ppfds[] = {0.0, 10.0, 100.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 2.0, 4.0};
  const int layer_counts[] = {1, 2, 3, 5};
  long points = 0, equal = 0, gravity = 0, atmospheric = 0, refused = 0;
  std::string first_gravity, first_equal;
  std::vector<long> by_kind_refused(20, 0), by_kind_points(20, 0);
  for (double psi : psi_soils) {
    for (double ppfd : ppfds) {
      for (double vpd : vpds) {
        for (int layers : layer_counts) {
          const fixture::Physiology d = drivers(psi, ppfd, vpd, layers);
          phylloptim::Leaf l;
          d.drive(l, kTheta);
          l.find_root_collar_psi();
          ++points;
          const int kind = int(l.operating_point_kind());
          ++by_kind_points[std::size_t(kind)];
          const double p = l.opt_root_psi_;
          const std::vector<double>& ps = l.supply_psi_soil();
          bool any_equal = false, any_gravity = false, any_atm = false;
          for (int i = 0; i < l.roots_.max_soil_layer; ++i) {
            const double gap = p - ps[std::size_t(i)];
            any_equal = any_equal || std::abs(gap) < 1e-8;
            any_gravity =
                any_gravity ||
                std::abs(gap - l.roots_.grav_head_z_[std::size_t(i)]) < 1e-8;
            any_atm = any_atm || std::abs(p) < 1e-8;
          }
          equal += any_equal;
          gravity += any_gravity;
          atmospheric += any_atm;
          const bool nan_here =
              !std::isfinite(l.dE_from_soil_dpsi_collar(p, ps));
          refused += nan_here;
          if (nan_here) {
            ++by_kind_refused[std::size_t(kind)];
          }
          if (any_gravity && first_gravity.empty()) {
            first_gravity = "psi_soil " + std::to_string(psi) + ", ppfd " +
                            std::to_string(int(ppfd)) + ", " +
                            std::to_string(layers) + " layers, " +
                            phylloptim::Leaf::operating_point_kind_name(
                                l.operating_point_kind());
          }
          if (any_equal && first_equal.empty()) {
            first_equal = "psi_soil " + std::to_string(psi) + ", ppfd " +
                          std::to_string(int(ppfd)) + ", " +
                          std::to_string(layers) + " layers";
          }
        }
      }
    }
  }
  std::printf("  %ld operating points\n", points);
  std::printf("  equal potentials      %ld   %s\n", equal,
              first_equal.empty() ? "(never)" : first_equal.c_str());
  std::printf("  gravity balance       %ld   %s\n", gravity,
              first_gravity.empty() ? "(never)" : first_gravity.c_str());
  std::printf("  collar at atmospheric %ld\n", atmospheric);
  std::printf("  the conductance refuses at %ld of them\n", refused);
  for (std::size_t k = 0; k < by_kind_points.size(); ++k) {
    if (by_kind_refused[k] > 0) {
      std::printf("    %ld of %ld at %s\n", by_kind_refused[k],
                  by_kind_points[k],
                  phylloptim::Leaf::operating_point_kind_name(
                      phylloptim::Leaf::OperatingPointKind(k)));
    }
  }
}

// --- B. the limits, against the general branch approaching them ------------------

// The general branch's own quantities for one layer, at an arbitrary collar.
struct Layer {
  double Q = 0.0;    // span / integral -- the reciprocal of G's divided difference
  double r_R = 0.0;
  double E = 0.0;
};

Layer general(const phylloptim::MultiLayerRoots& r, double s, double p, int i) {
  Layer out;
  const double lo = std::min(p, s), hi = std::max(p, s);
  const double span = hi - lo;
  const double pos_lo = std::max(lo, 0.0), neg_hi = std::min(hi, 0.0);
  double integral = 0.0;
  if (pos_lo < hi) {
    integral += r.root_vuln_integral_at(hi) - r.root_vuln_integral_at(pos_lo);
  }
  if (lo < neg_hi) {
    integral += (neg_hi - lo);
  }
  out.Q = span / integral;
  out.r_R = r.network_.r_R_H_min[std::size_t(i)] * out.Q +
            r.network_.r_R_V_sum[std::size_t(i)];
  out.E = (s - p - r.grav_head_z_[std::size_t(i)]) / out.r_R;
  return out;
}

void limits() {
  std::printf("\n=== B. THE EQUAL-POTENTIALS LIMIT ===\n");
  const fixture::Physiology d = drivers(2.0, 900.0, 2.0, 3);
  phylloptim::Leaf l;
  d.drive(l, kTheta);
  l.find_root_collar_psi();
  const phylloptim::MultiLayerRoots& r = l.roots_;
  const int i = 1;
  const double p = l.supply_psi_soil()[std::size_t(i)];
  // f and its two derivatives at the meeting point, all off the same curve the
  // general branch's own derivative reads.
  const double f = r.root_vuln_integral_deriv_at(p);
  const double fp = r.root_vuln_integrand_deriv_at(p);
  // f'' = f c x / psi^2 * (c x - c + 1) with x = (psi/b)^c, from differentiating
  // f' = -f c x / psi once more.
  const double x = std::pow(p / r.root_b, r.root_c);
  const double fpp = std::exp(-x) * r.root_c * x / (p * p) *
                     (r.root_c * x - r.root_c + 1.0);
  std::printf("  layer %d meets the collar at psi = %.10g\n", i, p);
  std::printf("  f = %.12g   f' = %.12g   f'' = %.12g\n", f, fp, fpp);
  // f'' against a difference of f'.
  {
    double best = 1e300;
    for (int e = 3; e <= 7; ++e) {
      const double h = std::pow(10.0, -double(e));
      const double got = (r.root_vuln_integrand_deriv_at(p + h) -
                          r.root_vuln_integrand_deriv_at(p - h)) / (2.0 * h);
      best = std::min(best, std::abs(fpp / got - 1.0));
    }
    std::printf("  f'' against a difference of f', best over a sweep: %.3g\n",
                best);
  }

  const double k = r.network_.r_R_H_min[std::size_t(i)];
  const double v = r.network_.r_R_V_sum[std::size_t(i)];
  const double grav = r.grav_head_z_[std::size_t(i)];

  // The claimed limits.
  const double Q0 = 1.0 / f;
  const double Qs = -fp / (2.0 * f * f);
  const double Qss = fp * fp / (2.0 * f * f * f) - fpp / (3.0 * f * f);
  const double r_R0 = k * Q0 + v;
  const double num0 = -grav;
  const double dr = k * Qs;
  const double d2r = k * Qss;
  const double E0 = num0 / r_R0;
  const double Es = (1.0 - num0 * dr / r_R0) / r_R0;
  const double Ess = -num0 * d2r / (r_R0 * r_R0) - 2.0 * Es * dr / r_R0;

  std::printf("\n  %-8s %-16s %-16s %-16s %-16s\n", "d", "Q(general)",
              "dQ/ds(diff)", "E(general)", "dE/ds(diff)");
  for (int e = 3; e <= 7; ++e) {
    const double h = std::pow(10.0, -double(e));
    const Layer at = general(r, p + h, p, i);
    const Layer up = general(r, p + 2.0 * h, p, i);
    const Layer dn = general(r, p, p, i);  // exactly the limit point
    static_cast<void>(dn);
    const Layer hi2 = general(r, p + h + h, p, i);
    const double dQ = (hi2.Q - at.Q) / h;      // one-sided, from h to 2h
    const double dE = (hi2.E - at.E) / h;
    static_cast<void>(up);
    std::printf("  1e-%-5d %-16.12g %-16.12g %-16.12g %-16.12g\n", e, at.Q, dQ,
                at.E, dE);
  }
  std::printf("  limit    %-16.12g %-16.12g %-16.12g %-16.12g\n", Q0, Qs, E0, Es);
  std::printf("  and d2E/ds2 at the limit = %.12g\n", Ess);

  // The two splines the value branch and the derivative branch read for the SAME
  // integrand, which is a latent disagreement rather than one of mine.
  std::printf("\n  f from the integral's own slope   %.17g\n", f);
  std::printf("  f from the conductivity spline    %.17g\n", r.root_vuln_at(p));
  std::printf("  they differ by                    %.3g relative\n",
              std::abs(r.root_vuln_at(p) / f - 1.0));
}

// --- C. the gravity balance and the collar at atmospheric ------------------------

void the_two_that_need_nothing() {
  std::printf("\n=== C. THE TWO THAT NEED NO LIMIT ===\n");
  const fixture::Physiology d = drivers(2.0, 900.0, 2.0, 3);
  phylloptim::Leaf l;
  d.drive(l, kTheta);
  l.find_root_collar_psi();
  const phylloptim::MultiLayerRoots& r = l.roots_;
  const int i = 1;
  const double p = l.supply_psi_soil()[std::size_t(i)];
  const double grav = r.grav_head_z_[std::size_t(i)];

  // At the gravity balance the numerator vanishes and nothing else does, so the
  // general expression is continuous and so is its derivative.
  const double s_bal = p + grav;
  std::printf("  gravity balance for layer %d at collar %.12g (grav %.6g)\n", i,
              s_bal, grav);
  const Layer at = general(r, s_bal, p, i);
  std::printf("    span %.6g  Q %.6g  r_R %.6g  E %.3g (exactly zero: %d)\n",
              std::abs(s_bal - p), at.Q, at.r_R, at.E, at.E == 0.0);
  double best = 1e300;
  const double dE_claim = 1.0 / at.r_R -
      0.0;  // num == 0, so the quotient rule leaves 1/r_R
  for (int e = 3; e <= 7; ++e) {
    const double h = std::pow(10.0, -double(e));
    const double got = (general(r, s_bal + h, p, i).E -
                        general(r, s_bal - h, p, i).E) / (2.0 * h);
    best = std::min(best, std::abs(dE_claim / got - 1.0));
  }
  std::printf("    dE/ds there is 1/r_R = %.12g; against a difference, best "
              "%.3g\n",
              dE_claim, best);

  // At the collar at atmospheric the integrand is 1 from both sides.
  std::printf("  collar at atmospheric: f_r(0) from the integral's slope = "
              "%.17g\n",
              r.root_vuln_integral_deriv_at(0.0));
  std::printf("    and the below-surface contribution is linear, slope 1 -- so "
              "the two sides agree\n");
}

}  // namespace

int main() {
  incidence();
  limits();
  the_two_that_need_nothing();
  return 0;
}
