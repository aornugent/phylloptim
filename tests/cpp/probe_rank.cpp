// PROBE: is the leaf's supply-side dependence really rank one at a held collar and
// rank two at the condition? The design reports assert it and measure it on an
// earlier tree; this asks THIS model, and it is the claim everything else rests on.
//
// The test is direct rather than statistical. At a HELD collar, if profit reads the
// supply only through total uptake E, then for every input u_k that reaches profit
// through the supply and nowhere else,
//
//     dprofit/du_k  =  c * dE/du_k        one constant c, every k
//
// and if the marginal profit R reads it only through (E, S = dE/dp), then
//
//     dR/du_k       =  a * dE/du_k + b * dS/du_k
//
// with one pair (a, b). Both are falsifiable by least squares over the inputs: an
// exact rank means a residual at rounding.
//
// ⚠️ AND THE CONTROL IS THE POINT. An input that reaches profit DIRECTLY -- the
// photosynthetic traits, the light, the stem parameters -- must FAIL both fits. A
// test that passes for everything is measuring nothing.
//
// Everything runs at FReal<double>: one forward tangent, no tape, no adjoint.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include \
//        ODELIA_SRC=../../../odelia/src probe_rank && ./probe_rank

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace pl = phylloptim;
using T = pl::tangent;  // FReal<double>

namespace {

const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_up(pl::Leaf& l, int layers, double psi0) {
  double theta[pl::gradient::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::gradient::par_R_d_25] = 1.44;
  l.set_traits(theta);
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / kAreaLeaf / layers);
    // Distinct potentials per layer, so the layers are not interchangeable and a
    // rank result cannot come from an accidental symmetry.
    psi_soil.push_back(psi0 + 0.35 * i);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

fixture::Soil soil_of(int layers) {
  fixture::Soil s;
  for (int i = 0; i < layers; ++i) {
    s.depth.push_back(1.0 * (i + 1));
    s.carbon.push_back(1.0 / layers / kAreaLeaf);
  }
  return s;
}

// One input's four directional derivatives, all at a HELD collar.
struct Row {
  std::string name;
  bool through_supply;  // the claim applies to these
  double dprofit = 0.0, dE = 0.0, dS = 0.0, dR = 0.0;
};

// Seeds one field of a LeafInputs<T> by index into a flat list, so the caller does
// not have to know the layout.
std::vector<T*> flat(pl::LeafInputs<T>& in) {
  std::vector<T*> out = in.profit.field_ptrs();
  for (T& v : in.supply.psi_soil) out.push_back(&v);
  for (T& v : in.supply.r_R_H_min) out.push_back(&v);
  for (T& v : in.supply.r_R_V_sum) out.push_back(&v);
  out.push_back(&in.supply.root_b);
  out.push_back(&in.supply.root_c);
  out.push_back(&in.psi_crit);
  out.push_back(&in.root_psi_crit);
  return out;
}

std::vector<std::string> names(int L) {
  std::vector<std::string> n = {"vcmax",  "jmax",   "quantum", "curv_elec",
                                "curv_colim", "ppfd", "resp",  "kmax",
                                "stem_b", "stem_c", "beta2",   "cost_scale"};
  for (int i = 0; i < L; ++i) n.push_back("psi_soil[" + std::to_string(i) + "]");
  for (int i = 0; i < L; ++i) n.push_back("r_R_H_min[" + std::to_string(i) + "]");
  for (int i = 0; i < L; ++i) n.push_back("r_R_V_sum[" + std::to_string(i) + "]");
  n.push_back("root_b");
  n.push_back("root_c");
  n.push_back("psi_crit");
  n.push_back("root_psi_crit");
  return n;
}

// Which inputs the claim is about: everything that reaches profit ONLY by moving the
// water the roots supply. The photosynthetic and stem parameters reach it directly
// and are the control.
bool reaches_only_through_supply(const std::string& n) {
  return n.rfind("psi_soil", 0) == 0 || n.rfind("r_R_H_min", 0) == 0 ||
         n.rfind("r_R_V_sum", 0) == 0 || n == "root_b" || n == "root_c";
}

std::vector<Row> rows_at(pl::Leaf& l, const fixture::Soil& soil, int L) {
  const pl::LeafInputs<double> base =
      fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);
  const std::vector<std::string> nm = names(L);
  std::vector<Row> out;
  const std::size_t n = nm.size();
  for (std::size_t k = 0; k < n; ++k) {
    pl::LeafInputs<T> in = base.template rebind_from<T>();
    std::vector<T*> f = flat(in);
    if (f.size() != n) {
      std::printf("  LAYOUT MISMATCH: %zu pointers, %zu names\n", f.size(), n);
      return {};
    }
    pl::seed_direction(*f[k], 1.0);

    // The collar is HELD: a passive double, carrying no direction. These are the
    // partials at fixed p, which is what the rank claim is about.
    const T collar = T(l.opt_root_psi_);

    const pl::Leaf::LeafOutputs<T> o = l.outputs_at<T>(collar, in);
    T e_total = T(0.0);
    for (const T& u : o.uptake) e_total += u;
    const T s = l.roots_.template duptake_dpsi_at<T>(collar, in.supply.at());
    const T r = l.marginal_at<T>(collar, in);

    Row row;
    row.name = nm[k];
    row.through_supply = reaches_only_through_supply(nm[k]);
    row.dprofit = pl::derivative_along(o.profit);
    row.dE = pl::derivative_along(e_total);
    row.dS = pl::derivative_along(s);
    row.dR = pl::derivative_along(r);
    out.push_back(row);
  }
  return out;
}

// Least squares of y on the columns given, and the worst RELATIVE residual over the
// rows used. A claim of exact rank predicts rounding.
struct Fit {
  double c0 = 0.0, c1 = 0.0, worst_rel = 0.0;
  std::string worst_name;
};

Fit fit1(const std::vector<Row>& rows, bool want_supply,
         double Row::*y, double Row::*x) {
  double xx = 0.0, xy = 0.0;
  for (const Row& r : rows) {
    if (r.through_supply != want_supply) continue;
    xx += r.*x * r.*x;
    xy += r.*x * r.*y;
  }
  Fit f;
  f.c0 = (xx > 0.0) ? xy / xx : 0.0;
  for (const Row& r : rows) {
    if (r.through_supply != want_supply) continue;
    const double res = r.*y - f.c0 * r.*x;
    const double rel = std::abs(res) / (std::abs(r.*y) + 1e-300);
    if (rel > f.worst_rel) { f.worst_rel = rel; f.worst_name = r.name; }
  }
  return f;
}

Fit fit2(const std::vector<Row>& rows, bool want_supply,
         double Row::*y, double Row::*x0, double Row::*x1) {
  double a00 = 0, a01 = 0, a11 = 0, b0 = 0, b1 = 0;
  for (const Row& r : rows) {
    if (r.through_supply != want_supply) continue;
    a00 += r.*x0 * r.*x0; a01 += r.*x0 * r.*x1; a11 += r.*x1 * r.*x1;
    b0 += r.*x0 * r.*y;   b1 += r.*x1 * r.*y;
  }
  Fit f;
  const double det = a00 * a11 - a01 * a01;
  if (std::abs(det) > 0.0) {
    f.c0 = (b0 * a11 - b1 * a01) / det;
    f.c1 = (a00 * b1 - a01 * b0) / det;
  }
  for (const Row& r : rows) {
    if (r.through_supply != want_supply) continue;
    const double res = r.*y - f.c0 * r.*x0 - f.c1 * r.*x1;
    const double rel = std::abs(res) / (std::abs(r.*y) + 1e-300);
    if (rel > f.worst_rel) { f.worst_rel = rel; f.worst_name = r.name; }
  }
  return f;
}

int failures = 0;

void ok(bool cond, const char* what, double got, const char* against) {
  std::printf("  %-52s %10.3e   %s   %s\n", what, got, against,
              cond ? "ok" : "FAIL");
  if (!cond) ++failures;
}

void one_state(int layers, double psi0, const char* label) {
  pl::Leaf l;
  set_up(l, layers, psi0);
  l.find_root_collar_psi();
  std::printf("\n=== %s: %d layers, psi_soil from %.2f, point %s ===\n", label,
              layers, psi0,
              pl::Leaf::operating_point_kind_name(l.operating_point_kind()));
  const std::vector<Row> rows = rows_at(l, soil_of(layers), layers);
  if (rows.empty()) { ++failures; return; }

  const Fit p1 = fit1(rows, true, &Row::dprofit, &Row::dE);
  ok(p1.worst_rel < 1e-10,
     "profit through E alone, supply inputs", p1.worst_rel, "< 1e-10");

  const Fit r2 = fit2(rows, true, &Row::dR, &Row::dE, &Row::dS);
  ok(r2.worst_rel < 1e-8,
     "the condition through (E, S), supply inputs", r2.worst_rel, "< 1e-8 ");

  // THE CONTROL. The same fits over inputs that reach profit directly must fail,
  // or the test above is measuring nothing.
  const Fit p1c = fit1(rows, false, &Row::dprofit, &Row::dE);
  ok(p1c.worst_rel > 1e-3,
     "CONTROL: direct inputs must NOT fit E alone", p1c.worst_rel, "> 1e-3 ");

  const Fit r2c = fit2(rows, false, &Row::dR, &Row::dE, &Row::dS);
  ok(r2c.worst_rel > 1e-3,
     "CONTROL: direct inputs must NOT fit (E, S)", r2c.worst_rel, "> 1e-3 ");

  std::printf("  fitted   dprofit/dE = %.12g\n", p1.c0);
  std::printf("           dR/dE      = %.12g    dR/dS = %.12g\n", r2.c0, r2.c1);
  // ⚠️ CONDITIONING, and the right measure is not a ratio of extremes. A rank-two
  // fit is vacuous when the two columns are collinear, so what matters is how much
  // of the S column is INDEPENDENT of the E column: 1 - corr^2 over the supply
  // inputs. At zero the second coefficient is unidentifiable and the residual below
  // means nothing.
  double a00 = 0, a01 = 0, a11 = 0; int used = 0;
  for (const Row& r : rows) {
    if (!r.through_supply) continue;
    ++used;
    a00 += r.dE * r.dE; a01 += r.dE * r.dS; a11 += r.dS * r.dS;
  }
  const double corr2 = (a00 > 0 && a11 > 0) ? (a01 * a01) / (a00 * a11) : 1.0;
  std::printf("  conditioning: %d supply inputs, S independent of E by %.4f\n",
              used, 1.0 - corr2);
  ok(used >= 5 && (1.0 - corr2) > 1e-3,
     "the two columns are not collinear", 1.0 - corr2, "> 1e-3 ");
  std::printf("  worst supply-side residual at %s / %s\n",
              p1.worst_name.empty() ? "-" : p1.worst_name.c_str(),
              r2.worst_name.empty() ? "-" : r2.worst_name.c_str());
}

}  // namespace

int main() {
  std::printf("is the leaf's supply-side dependence rank one, and its condition rank two?\n");
  // Across the operating-point kinds the solve can reach, because the claim has to
  // hold wherever plant asks, not only at an interior maximum.
  one_state(5, 2.00, "wet");
  one_state(5, 3.40, "drier");
  one_state(5, 4.60, "drier still");
  one_state(5, 5.40, "near the root limit");
  one_state(3, 2.00, "three layers");
  one_state(1, 2.00, "one layer");
  one_state(1, 5.60, "one layer, dry");
  std::printf("\n%s\n", failures == 0 ? "all ok" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
