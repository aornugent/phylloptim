// PROBE: what the leaf boundary costs, region by region and kernel by kernel.
//
// ⚠️ READ THE OPERATIONS, NOT THE STATEMENTS. Cost is operations times the price
// of one; a statement is the tape's bookkeeping over a run of them, and counting
// statements alone prices the bookkeeping and nothing else. It is also blind
// twice over: to arithmetic that never reaches a tape at all, and to how WIDE a
// scalar's derivative is, which is the entire cost of a nested tangent. Both
// blind spots have already cost this project a wrong conclusion -- a curve
// dismissed at four statements turned out to be a quarter of the instructions.
//
// For what the machine actually pays, run the whole probe under callgrind, whose
// count is exact where seconds on this fixture have swung 28% between sittings
// and 81% for one binary:
//
//   valgrind --tool=callgrind --callgrind-out-file=cg.out ./probe_tape_regions
//   callgrind_annotate cg.out
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
#include <cstdlib>
#include <string>
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
  long draw_op = 0, collar_op = 0, outputs_op = 0;
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
  const long s0 = long(tape.getNumStatements()), o0 = long(tape.getNumOperations());
  const auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), supply);
  const long s1 = long(tape.getNumStatements()), o1 = long(tape.getNumOperations());
  const A collar = l.collar_at<K, A>(draw, io.in, curv);
  const long s2 = long(tape.getNumStatements()), o2 = long(tape.getNumOperations());
  const pl::Leaf::LeafOutputs<A> got = l.outputs_at<K, A>(collar, draw, io.in);
  const long s3 = long(tape.getNumStatements()), o3 = long(tape.getNumOperations());
  r.draw = s1 - s0;      r.draw_op = o1 - o0;
  r.collar = s2 - s1;    r.collar_op = o2 - o1;
  r.outputs = s3 - s2;   r.outputs_op = o3 - o2;
  // ⚠️ collar_at REWINDS, and resetTo TRUNCATES BOTH COUNTERS, so its figures are
  // what SURVIVED and not what it pushed. marginal_at below is the same region
  // measured before the rewind, which is why the two disagree by 30-fold.

  // The same coordinates twice: once at the collar profit is actually handed at
  // an interior point, once at the live one a bound hands it. The difference is
  // the collar's channel.
  const A held(odelia::util::to_passive(collar));
  const long c0 = long(tape.getNumStatements());
  volatile double sink = 0.0;
  {
    const auto at = l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, held, draw,
                                          io.in, l.photo_capacity_at<A>(io.in),
                                          false);
    sink += odelia::util::to_passive(at.sigma.value);
  }
  const long c1 = long(tape.getNumStatements());
  {
    const auto at = l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, collar, draw,
                                          io.in, l.photo_capacity_at<A>(io.in),
                                          true);
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

  // ⚠️ EVERY OUTPUT, NOT ONLY PROFIT. At an interior point profit reads the HELD
  // collar, so the collar's own rows reach the UPTAKES and never the objective --
  // and a check that sweeps profit alone is green whatever happens to the collar.
  // One sweep per output, derivatives cleared between, rows concatenated.
  std::vector<A*> outs{const_cast<A*>(&got.profit)};
  for (const A& u : got.uptake) outs.push_back(const_cast<A*>(&u));
  for (A* o : outs) tape.registerOutput(*o);
  for (A* o : outs) {
    tape.clearDerivatives();
    xad::derivative(*o) = 1.0;
    tape.computeAdjoints();
    for (double v : io.rows()) r.rows.push_back(v);
  }
  r.profit = odelia::util::to_passive(got.profit);
  r.ok = true;
  return r;
}


// The three calls, N times on one tape, and nothing else -- so a run under
// callgrind prices a placement in INSTRUCTIONS, which is the only one of the
// three quantities that sees a tangent's arithmetic. plant records about a
// hundred placements between one newRecording and the next, so the tape is
// shared across the loop as it is there.
void placements(int n) {
  pl::Leaf l = set_up(1, 0.5);
  double curv = 0.0;
  if (!interior(l, curv)) { printf("(fixture not interior)\n"); return; }
  Tape tape;
  Inputs io;
  io.seed(tape, l);
  const pl::SupplyAt<A> supply = io.supply();
  volatile double sink = 0.0;
  for (int i = 0; i < n; ++i) {
    const auto draw = l.supply_draw_at<A>(A(l.opt_root_psi_), supply);
    const A collar = l.collar_at<K, A>(draw, io.in, curv);
    const pl::Leaf::LeafOutputs<A> got = l.outputs_at<K, A>(collar, draw, io.in);
    sink += odelia::util::to_passive(got.profit);
  }
  printf("%d placements, %u statements, %u operations, sink %.17g\n", n,
         tape.getNumStatements(), tape.getNumOperations(), double(sink));
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
  long s0 = 0, o0 = 0;
  auto mark = [&]() {
    s0 = long(tape.getNumStatements());
    o0 = long(tape.getNumOperations());
  };
  auto say = [&](const char* what) {
    printf("  %-42s %5ld stmt %6ld op\n", what,
           long(tape.getNumStatements()) - s0, long(tape.getNumOperations()) - o0);
  };
  auto keep = [&](const A& v) { sink += odelia::util::to_passive(v); };

  printf("\nleaf kernels, per call\n");
  mark(); keep(l.vcmax_at<A>(io.in));                        say("vcmax_at(pars)");
  mark(); keep(l.electron_transport_at<A>(io.in));           say("electron_transport_at(pars)");
  const pl::Leaf::PhotoCapacity<A> cap = l.photo_capacity_at<A>(io.in);
  mark(); keep(l.assim_colimited_kernel<A>(ci, cap));        say("assim_colimited_kernel(ci, capacity)");
  mark(); keep(l.stem_integral_at<A>(sigma, io.in));         say("stem_integral_at(psi, pars)");
  mark(); keep(l.transpiration_at<A>(sigma, collar, io.in)); say("transpiration_at(sigma, collar, pars)");
  mark(); keep(l.hydraulic_cost_TF_kernel<A>(sigma, io.in)); say("hydraulic_cost_TF_kernel(sigma, pars)");
  mark(); keep(pl::closed_form_curve<A>(l.stem_curve_integral_deriv(l.opt_psi_stem_),
                                        sigma, io.in[pl::par_stem_P50],
                                        io.in[pl::par_stem_c]));
                                                             say("closed_form_curve(table, psi, P50, c)");
  mark(); keep(l.assim_slope_at<A>(ci, cap));                say("assim_slope_at(ci, capacity)  <- SUPPLIED");
  mark(); keep(l.cost_slope_at<A>(sigma, io.in));            say("cost_slope_at(sigma, pars)    <- SUPPLIED");
  (void)sink;
}

}  // namespace

int main(int argc, char** argv) {
  // `placements N` runs the boundary and nothing else, for callgrind.
  if (argc > 1 && std::string(argv[1]) == "placements") {
    placements(argc > 2 ? std::atoi(argv[2]) : 60);
    return 0;
  }
  printf("leaf boundary per placement, statements/operations by region\n"
         "(collar_at REWINDS, so its pair is what SURVIVED; marginal_at is the\n"
         " same region measured before the rewind)\n\n");
  printf("%-6s %-5s %14s %14s %14s | %9s %9s %9s %9s\n", "layers", "psi0",
         "draw stmt/op", "collar stmt/op", "outputs stmt/op", "cd@held",
         "cd@live", "profit", "marginal");
  const int layer_set[] = {1, 3, 5};
  const double psi_set[] = {0.5, 1.0, 2.0};
  for (int layers : layer_set) {
    for (double psi0 : psi_set) {
      const Regions r = measure(layers, psi0);
      if (!r.ok) { printf("%-6d %-5.2f  (not interior)\n", layers, psi0); continue; }
      printf("%-6d %-5.2f %7ld/%-6ld %7ld/%-6ld %7ld/%-6ld | %9ld %9ld %9ld %9ld\n",
             layers, psi0, r.draw, r.draw_op, r.collar, r.collar_op, r.outputs,
             r.outputs_op, r.coords_held, r.coords_live, r.profit_held,
             r.marginal);
    }
  }

  // What each kernel costs, which is where a region total has to be explained.
  // ⚠️ THE TWO SLOPES ARE THE CHEAPEST ROWS HERE AND THE MOST EXPENSIVE CALLS.
  // Each records one statement whatever its row count and runs a tangent of a
  // tangent that no counter here can see, so read them under callgrind --
  // `placements` is the mode for that -- and not off this table.
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
