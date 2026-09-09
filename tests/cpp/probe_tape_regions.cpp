// PROBE: where the leaf boundary's tape statements go, region by region, and
// which of them reach a row.
//
// plant's TF24_Strategy::record_leaf_outputs makes three calls -- the draw, the
// collar, the outputs -- and the century stand's cost is dominated by what they
// record per placement. This measures each region directly off the tape's own
// counters, so a claim about where the statements are is a number rather than a
// profile attribution.
//
// It also splits the one question the region totals cannot answer: profit_at is
// handed a HELD collar at an interior point (the envelope omission), so anything
// collar_coords_at computes for the collar's channel is multiplied by a step
// that is zero in value AND carries no derivative. Recorded against swept is the
// difference between work that reaches a row and work that does not.
//
//   make CXX=g++ probe_tape_regions && ./probe_tape_regions
#include <phylloptim.hpp>
#include "root_network.hpp"
#include <odelia/implicit_node.hpp>
#include <odelia/ode_interface.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

namespace pl = phylloptim;
using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;
constexpr auto K = pl::Leaf::CostCurve::TF24;

namespace {

pl::Leaf set_up(int layers, double psi0) {
  pl::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  const double th = 0.5, ks = 2.0, h = 5.0, area = 1.0;
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / area / layers);
    psi_soil.push_back(psi0 + 0.35 * i);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   ks * th / h, 2.0, 40.0, 25.0, 21.0, 101.3);
  return l;
}

// Everything the tape holds as an input, so a row can be read for each.
struct Inputs {
  pl::leaf_pars<A> in;
  std::vector<A> psi_soil, rh, rv;
  void seed(Tape& tape, const pl::Leaf& l) {
    const pl::leaf_pars<double> p = l.passive_pars();
    for (std::size_t i = 0; i < p.size(); ++i) in[i] = A(p[i]);
    const int n = l.supply_n_layers();
    for (int i = 0; i < n; ++i) {
      psi_soil.push_back(A(l.roots_.psi_soil_[std::size_t(i)]));
      rh.push_back(A(l.roots_.network_.r_R_H_min[std::size_t(i)]));
      rv.push_back(A(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
    }
    for (A& v : in) tape.registerInput(v);
    for (A& v : psi_soil) tape.registerInput(v);
    for (A& v : rh) tape.registerInput(v);
    for (A& v : rv) tape.registerInput(v);
    tape.newRecording();
  }
  pl::SupplyAt<A> supply() const {
    return pl::SupplyAt<A>{psi_soil, rh, rv, in[pl::par_root_P50],
                           in[pl::par_root_c]};
  }
  std::vector<double> rows() const {
    std::vector<double> r;
    for (const A& v : in) r.push_back(xad::derivative(v));
    for (const A& v : psi_soil) r.push_back(xad::derivative(v));
    for (const A& v : rh) r.push_back(xad::derivative(v));
    for (const A& v : rv) r.push_back(xad::derivative(v));
    return r;
  }
};

// A solved interior leaf and the curvature its collar residual divides by;
// every other kind is skipped, because the century fixture reaches only this one.
bool interior(pl::Leaf& l, double& curv) {
  l.find_root_collar_psi();
  if (l.operating_point_kind() != pl::Leaf::OperatingPointKind::Interior)
    return false;
  try { curv = l.marginal_collar_slope<K>(); } catch (...) { return false; }
  return curv < 0.0;
}

struct Regions {
  bool ok = false;
  long draw = 0, collar = 0, outputs = 0;
  long coords_held = 0, coords_live = 0, profit_held = 0, marginal = 0;
  double profit = 0.0;
  std::vector<double> rows;
};

Regions measure(int layers, double psi0) {
  Regions r;
  pl::Leaf l = set_up(layers, psi0);
  double curv = 0.0;
  if (!interior(l, curv)) return r;

  Tape tape;
  Inputs io;
  io.seed(tape, l);
  const pl::SupplyAt<A> supply = io.supply();

  // The three calls plant makes, in its order, each charged what it added.
  const long s0 = long(tape.getNumStatements());
  const auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), supply);
  const long s1 = long(tape.getNumStatements());
  const A collar = l.collar_at<K, A>(draw, io.in, curv);
  const long s2 = long(tape.getNumStatements());
  const pl::Leaf::LeafOutputs<A> got = l.outputs_at<K, A>(collar, draw, io.in);
  const long s3 = long(tape.getNumStatements());
  r.draw = s1 - s0;
  r.collar = s2 - s1;
  r.outputs = s3 - s2;

  // The same coordinates twice: once at the collar profit is actually handed at
  // an interior point, once at the live one a bound hands it. The difference is
  // the collar's channel.
  const A held(odelia::util::to_passive(collar));
  const long c0 = long(tape.getNumStatements());
  volatile double sink = 0.0;
  {
    const auto at = l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, held, draw, io.in, false);
    sink += odelia::util::to_passive(at.sigma.value);
  }
  const long c1 = long(tape.getNumStatements());
  {
    const auto at = l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, collar, draw, io.in, true);
    sink += odelia::util::to_passive(at.sigma.value);
  }
  const long c2 = long(tape.getNumStatements());
  {
    const A p = l.profit_at<K, A>(held, draw, io.in, false);
    sink += odelia::util::to_passive(p);
  }
  const long c3 = long(tape.getNumStatements());
  {
    const A m = l.marginal_at<K, A>(collar, draw, io.in);
    sink += odelia::util::to_passive(m);
  }
  const long c4 = long(tape.getNumStatements());
  r.coords_held = c1 - c0;
  r.coords_live = c2 - c1;
  r.profit_held = c3 - c2;
  r.marginal = c4 - c3;
  (void)sink;

  // The rows plant reads: the objective's, swept once.
  A* target = const_cast<A*>(&got.profit);
  tape.registerOutput(*target);
  xad::derivative(*target) = 1.0;
  tape.computeAdjoints();
  r.profit = odelia::util::to_passive(got.profit);
  r.rows = io.rows();
  r.ok = true;
  return r;
}


// Each kernel priced on its own, at the scalar the boundary calls it at.
void kernels() {
  pl::Leaf l = set_up(1, 0.5);
  double curv = 0.0;
  if (!interior(l, curv)) { printf("\n(kernel prices: fixture not interior)\n"); return; }
  Tape tape;
  Inputs io;
  io.seed(tape, l);
  const pl::SupplyAt<A> supply = io.supply();
  const A collar(l.opt_root_psi_), sigma(l.opt_psi_stem_), ci(l.ci_);
  volatile double sink = 0.0;
  long s0 = 0;
  auto mark = [&]() { s0 = long(tape.getNumStatements()); };
  auto say = [&](const char* what) {
    printf("  %-44s %6ld\n", what, long(tape.getNumStatements()) - s0);
  };
  auto keep = [&](const A& v) { sink += odelia::util::to_passive(v); };

  printf("\nleaf kernels, tape statements per call\n");
  mark(); keep(l.vcmax_at<A>(io.in));                        say("vcmax_at(pars)");
  mark(); keep(l.electron_transport_at<A>(io.in));           say("electron_transport_at(pars)");
  mark(); keep(l.assim_colimited_kernel<A>(ci, io.in));      say("assim_colimited_kernel(ci, pars)");
  mark(); keep(l.stem_integral_at<A>(sigma, io.in));         say("stem_integral_at(psi, pars)");
  mark(); keep(l.transpiration_at<A>(sigma, collar, io.in)); say("transpiration_at(sigma, collar, pars)");
  mark(); keep(l.hydraulic_cost_TF_kernel<A>(sigma, io.in)); say("hydraulic_cost_TF_kernel(sigma, pars)");
  mark(); keep(pl::closed_form_curve<A>(l.stem_curve_integral_deriv(l.opt_psi_stem_),
                                        sigma, io.in[pl::par_stem_P50],
                                        io.in[pl::par_stem_c]));
                                                             say("closed_form_curve(table, psi, P50, c)");
  mark(); keep(l.assim_slope_at<A>(ci, io.in));              say("assim_slope_at(ci, pars)   <- NESTED");
  (void)sink;
}

}  // namespace

int main() {
  printf("leaf boundary, tape statements by region (surviving, after any rewind)\n\n");
  printf("%-8s %-7s %9s %9s %9s | %11s %11s %11s %9s\n", "layers", "psi0",
         "draw", "collar", "outputs", "coords@held", "coords@live", "profit@held",
         "marginal");
  const int layer_set[] = {1, 3, 5};
  const double psi_set[] = {0.5, 1.0, 2.0};
  for (int layers : layer_set) {
    for (double psi0 : psi_set) {
      const Regions r = measure(layers, psi0);
      if (!r.ok) { printf("%-8d %-7.2f  (not interior)\n", layers, psi0); continue; }
      printf("%-8d %-7.2f %9ld %9ld %9ld | %11ld %11ld %11ld %9ld\n", layers,
             psi0, r.draw, r.collar, r.outputs, r.coords_held, r.coords_live,
             r.profit_held, r.marginal);
    }
  }

  // What each kernel costs, which is where a region total has to be explained.
  // ⚠️ A TANGENT NESTED ABOVE THE ACTIVE SCALAR IS WHAT odelia/tangent.hpp
  // FORBIDS, and the last two rows are it: every operand and cached value of a
  // nested BinaryExpr is copied, and a copy of an active scalar is a recorded
  // statement.
  kernels();

  // The rows, at full precision, so an edit that claims to remove dead work can
  // be held to bit-identity rather than to a tolerance.
  printf("\nprofit rows, 3 layers at psi0 = 0.5\n");
  const Regions r = measure(3, 0.5);
  if (!r.ok || r.rows.empty()) { printf("NO ROWS -- fixture not interior\n"); return 1; }
  {
    printf("rows %zu\n", r.rows.size());
    printf("profit %.17g\n", r.profit);
    for (std::size_t i = 0; i < r.rows.size(); ++i)
      printf("row %2zu %.17g\n", i, r.rows[i]);
  }
  return 0;
}
