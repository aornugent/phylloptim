// PROBE: is the supply's collar derivative still 0/0 where a layer's potential
// equals the collar?
//
// duptake_dpsi refuses there -- returns NaN, and its caller falls back to a central
// difference -- because the layer's mean conductivity was integral/span with both
// vanishing. Under an average that is no longer true: at span zero every node of the
// rule collapses onto the bound, so the mean is f(lo), its bound derivative is
// f'(lo)/2 and its second is f''(lo)/3, each an exact limit reached by arithmetic
// rather than by a quotient.
//
// This asks the model directly. It reports what duptake_dpsi returns AT the
// coincidence and, whatever that is, what the layer-mean family returns there --
// against a central difference of the uptake taken from far enough away that the
// difference is not itself in the singular region.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include \
//        ODELIA_SRC=../../../odelia/src probe_coincidence && ./probe_coincidence

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace pl = phylloptim;

namespace {

const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

// One layer, so the coincidence is the whole of what the collar sees.
void set_up(pl::Leaf& l, const std::vector<double>& psi_soil) {
  double theta[pl::gradient::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::gradient::par_R_d_25] = 1.44;
  l.set_traits(theta);
  const std::size_t n = psi_soil.size();
  std::vector<double> root(n, 1.0 / kAreaLeaf / double(n)), depth;
  for (std::size_t i = 0; i < n; ++i) depth.push_back(1.0 * double(i + 1));
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

double uptake_total(const pl::Leaf& l, const std::vector<double>& psi_soil,
                    double collar) {
  std::vector<double> per(psi_soil.size(), 0.0);
  double e_up = 0.0;
  l.roots_.uptake_at(collar, psi_soil, per, e_up);
  return e_up;
}

void at(const char* what, std::vector<double> psi_soil, double collar) {
  pl::Leaf l;
  set_up(l, psi_soil);
  l.roots_.begin_solve();

  const double analytic = l.roots_.duptake_dpsi(collar, psi_soil);
  // A wide central difference, deliberately: a step small enough to sit inside the
  // coincidence would be measuring the same region rather than refereeing it.
  const double h = 1e-4;
  const double diff = (uptake_total(l, psi_soil, collar + h) -
                       uptake_total(l, psi_soil, collar - h)) /
                      (2.0 * h);

  // The layer-mean family at the coincidence itself, which is what the refusal was
  // ever about.
  const double m = l.roots_.layer_mean(collar, collar);
  const double d1 = l.roots_.layer_mean_dbound(collar, collar, true, m);
  const double d2 = l.roots_.layer_mean_dbound2(collar, collar, true, m);

  std::printf("  %-22s collar %.6g\n", what, collar);
  std::printf("    duptake_dpsi        %-18.10g  central difference %.10g\n",
              analytic, diff);
  if (std::isfinite(analytic)) {
    std::printf("    rel gap             %.3e\n",
                std::abs(analytic - diff) / (std::abs(diff) + 1e-300));
  } else {
    std::printf("    rel gap             (refused)\n");
  }
  std::printf("    layer mean at span 0  f %.10g   f'/2 %.10g   f''/3 %.10g  %s\n",
              m, d1, d2,
              (std::isfinite(m) && std::isfinite(d1) && std::isfinite(d2))
                  ? "all finite"
                  : "NOT FINITE");
  // The exact limits the rule must reproduce, from the curve's own closed forms.
  const double b = l.roots_.root_b, c = l.roots_.root_c;
  std::printf("    against the closed forms  f %.10g   f'/2 %.10g   f''/3 %.10g\n",
              pl::vulnerability_curve_at<double>(collar, b, c),
              0.5 * pl::vulnerability_curve_slope_at<double>(collar, b, c),
              pl::vulnerability_curve_curvature_at<double>(collar, b, c) / 3.0);
}

}  // namespace

int main() {
  std::printf("the supply at a coincidence: is the mean still 0/0 there?\n\n");
  at("interior, one layer", {2.0}, 2.0);
  at("interior, three layers", {1.5, 2.0, 2.5}, 2.0);
  at("gravity-balanced", {2.0, 3.0}, 2.0);
  at("near the surface", {0.05}, 0.05);
  return 0;
}
