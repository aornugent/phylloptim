// PROBE: the two slopes the model computes instead of having.
//
// Every derivative mechanism in this tree traces to one missing primitive. The
// nested tangent above the adjoint exists because dA/dci is COMPUTED at run time
// rather than REPRESENTED, so a caller that wants A's parameter rows has to
// differentiate the kernel one order past what it holds. The same for dC/dsigma.
//
// Both are elementary. This writes them as sibling closed forms -- the shape
// vulnerability_curve_slope_at already has beside vulnerability_curve_at -- checks
// them against the tangent they replace, and counts what each costs the tape.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include \
//        ODELIA_SRC=../../../odelia/src probe_primitive && ./probe_primitive

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <odelia/ode_interface.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace pl = phylloptim;
using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;
using Tan = pl::tangent;

namespace {

int failures = 0;

void ok(const char* what, double got, double want, double tol) {
  const double rel = std::abs(got - want) / (std::abs(want) + 1e-300);
  std::printf("  %-44s %.15g\n  %-44s %.15g   rel %.3e  %s\n", what, got,
              "   against the tangent", want, rel, rel < tol ? "ok" : "FAIL");
  if (!(rel < tol)) ++failures;
}

// dA/dci for the colimited assimilation. A = (s - q)/(2t) - Rd with
// s = Ac + Aj and q = sqrt(s^2 - 4 t Ac Aj), so the slope is elementary in the two
// limited rates and their own slopes. Respiration drops out: it is subtracted.
template <class T>
T assim_colimited_slope(const T& ci, const T& vcmax, const T& transport,
                        const T& curvature, const T& gamma_pa, const T& km) {
  const T ac = vcmax * (ci - gamma_pa) / (ci + km);
  const T aj = transport / 4.0 * ((ci - gamma_pa) / (ci + 2.0 * gamma_pa));
  // d/dci of (ci - g)/(ci + k) is (k + g)/(ci + k)^2.
  const T dac = vcmax * (km + gamma_pa) / ((ci + km) * (ci + km));
  const T daj = transport / 4.0 * (3.0 * gamma_pa) /
                ((ci + 2.0 * gamma_pa) * (ci + 2.0 * gamma_pa));
  const T s = ac + aj;
  const T ds = dac + daj;
  const T q = sqrt(s * s - 4.0 * curvature * ac * aj);
  const T dq = (s * ds - 2.0 * curvature * (dac * aj + ac * daj)) / q;
  return (ds - dq) / (2.0 * curvature);
}

// dC/dsigma for the hydraulic cost. C = scale * (1 - f)^beta with f the stem's
// vulnerability curve, so this is one line over the slope that curve already has.
template <class T>
T hydraulic_cost_slope(const T& sigma, const T& b, const T& c, const T& beta,
                       const T& scale) {
  const T f = pl::vulnerability_curve_at<T>(sigma, b, c);
  const T df = pl::vulnerability_curve_slope_at<T>(sigma, b, c);
  return -scale * beta * pow(1.0 - f, beta - 1.0) * df;
}

const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_up(pl::Leaf& l) {
  double theta[pl::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::trait_R_d_25] = 1.44;
  l.set_traits(theta);
  std::vector<double> root{1.0 / kAreaLeaf}, psi_soil{2.0}, depth{1.0};
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

}  // namespace

int main() {
  pl::Leaf l;
  set_up(l);
  l.find_root_collar_psi();

  const double ci = l.ci_, sigma = l.opt_psi_stem_;
  const double gpa = l.gamma_ * l.umol_per_mol_to_Pa_;
  const double km = l.km_;

  std::printf("the two slopes the model computes instead of having\n");
  std::printf("  at ci = %.10g, sigma = %.10g\n\n", ci, sigma);

  // dA/dci: the closed form against a tangent through the kernel it belongs to.
  {
    Tan c_ad = Tan(ci);
    pl::seed_direction(c_ad, 1.0);
    const Tan a = l.assim_colimited_kernel<Tan>(
        c_ad, Tan(l.vcmax_), Tan(l.electron_transport_),
        Tan(l.curv_fact_colim), Tan(l.R_d_));
    const double want = pl::derivative_along(a);
    const double got = assim_colimited_slope<double>(
        ci, l.vcmax_, l.electron_transport_, l.curv_fact_colim, gpa, km);
    ok("dA/dci, closed form", got, want, 1e-13);
  }

  // dC/dsigma, the same way.
  {
    Tan s_ad = Tan(sigma);
    pl::seed_direction(s_ad, 1.0);
    const Tan cst = l.hydraulic_cost_TF_kernel<Tan>(
        s_ad, Tan(l.stem_b), Tan(l.stem_c), Tan(l.beta2),
        Tan(l.cost_scale_TF24));
    const double want = pl::derivative_along(cst);
    const double got = hydraulic_cost_slope<double>(
        sigma, l.stem_b, l.stem_c, l.beta2, l.cost_scale_TF24);
    ok("dC/dsigma, closed form", got, want, 1e-13);
  }

  // And the price. The three kernels at a nested tangent above the adjoint cost 566
  // statements; the same information as primitives costs what a value costs.
  {
    Tape tape;
    A ci_a = A(ci), sg_a = A(sigma), vc = A(l.vcmax_), tr = A(l.electron_transport_),
      cv = A(l.curv_fact_colim), rd = A(l.R_d_), bb = A(l.stem_b), cc = A(l.stem_c),
      b2 = A(l.beta2), sc = A(l.cost_scale_TF24), gp = A(gpa), kmA = A(km);
    for (A* p : {&ci_a, &sg_a, &vc, &tr, &cv, &rd, &bb, &cc, &b2, &sc}) {
      tape.registerInput(*p);
    }
    tape.newRecording();

    const std::size_t s0 = tape.getNumStatements();
    const A a_val = l.assim_colimited_kernel<A>(ci_a, vc, tr, cv, rd);
    const A c_val = l.hydraulic_cost_TF_kernel<A>(sg_a, bb, cc, b2, sc);
    const std::size_t s1 = tape.getNumStatements();
    const A a_slope = assim_colimited_slope<A>(ci_a, vc, tr, cv, gp, kmA);
    const A c_slope = hydraulic_cost_slope<A>(sg_a, bb, cc, b2, sc);
    const std::size_t s2 = tape.getNumStatements();

    std::printf("\n  the two kernels' VALUES at the adjoint      %4zu statements\n",
                s1 - s0);
    std::printf("  their SLOPES as primitives, same scalar    %4zu statements\n",
                s2 - s1);
    std::printf("  values and slopes together                 %4zu statements\n",
                s2 - s0);
    std::printf("  the same information by a nested tangent    566 statements"
                "   <- measured by probe_leaf_tape\n");
    std::printf("  so a primitive costs                        %.1fx less\n",
                566.0 / double(s2 - s0));
    (void)a_val; (void)c_val; (void)a_slope; (void)c_slope;
  }

  std::printf("\n%s\n", failures == 0 ? "all ok" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
