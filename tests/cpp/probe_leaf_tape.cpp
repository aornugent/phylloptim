// PROBE: what one leaf placement costs the tape it is recorded on, counted
// rather than sampled. A profile says which primitives are hot; this says how
// many statements and operations a placement writes, which is the number a
// design decision turns on -- plant records one of these per cohort per stage
// per step and sweeps the result once per census metric.
//
// It also prices the alternative in the same units: a dense block, where the
// leaf hands over each output's value with its derivative rows as plain
// doubles, costs one statement per non-zero row and nothing else.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include \
//        ODELIA_SRC=../../../odelia/src probe_leaf_tape && ./probe_leaf_tape

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <odelia/ode_interface.hpp>
#include <odelia/implicit_node.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace pl = phylloptim;

using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;
using clock_type = std::chrono::steady_clock;

namespace {

// The default trait vector and drivers, from bench_solve/bench_gradient, so this
// measures the same interior operating point the rest of the suite does.
const double kBase[13] = {96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                          5.870283, 1.5, 157.44, 0.30, 0.7, 0.99, 7.5};
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

void set_up(pl::Leaf& l, int layers) {
  double theta[pl::gradient::n_traits];
  std::copy(kBase, kBase + 13, theta);
  theta[pl::gradient::par_R_d_25] = 1.44;
  l.set_traits(theta);
  std::vector<double> root, psi_soil, depth;
  for (int i = 0; i < layers; ++i) {
    root.push_back(1.0 / kAreaLeaf / layers);
    psi_soil.push_back(2.0);
    depth.push_back(1.0 * (i + 1));
  }
  l.set_physiology(fixture::root_network(root, depth), 900.0, psi_soil, depth,
                   kKs * kTheta / kH, 2.0, 40.0, 25.0, 21.0, 101.3);
}

// Every active value the leaf reads, in one list, so the recording has real
// inputs -- an active assigned from a double takes no slot and records nothing,
// so an unregistered run would measure an empty tape.
std::vector<A*> inputs_of(pl::LeafInputs<A>& in) {
  std::vector<A*> out = in.profit.field_ptrs();
  for (A& v : in.supply.psi_soil) out.push_back(&v);
  for (A& v : in.supply.r_R_H_min) out.push_back(&v);
  for (A& v : in.supply.r_R_V_sum) out.push_back(&v);
  out.push_back(&in.supply.root_b);
  out.push_back(&in.supply.root_c);
  out.push_back(&in.psi_crit);
  out.push_back(&in.root_psi_crit);
  return out;
}

fixture::Soil soil_of(int layers) {
  fixture::Soil s;
  for (int i = 0; i < layers; ++i) {
    s.depth.push_back(1.0 * (i + 1));
    s.carbon.push_back(1.0 / layers / kAreaLeaf);
  }
  return s;
}

double us_per(clock_type::time_point a, clock_type::time_point b, long n) {
  return std::chrono::duration<double, std::micro>(b - a).count() / double(n);
}

void one_width(int layers) {
  pl::Leaf l;
  set_up(l, layers);
  l.find_root_collar_psi();
  const fixture::Soil soil = soil_of(layers);

  std::printf("\n=== %d soil layers ===\n", layers);
  std::printf("  operating point            %s\n",
              pl::Leaf::operating_point_kind_name(l.operating_point_kind()));
  std::fflush(stdout);

  // Constructing a Tape activates it; the codebase spells the other case
  // `tape(false)`. Activating again raises TapeAlreadyActive.
  Tape tape;
  // Built at double and rebound, which is the route plant takes in reverse: the
  // fixture's seeding helper is a tangent's and refuses an adjoint scalar.
  const pl::LeafInputs<double> in_d =
      fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);
  pl::LeafInputs<A> in = in_d.template rebind_from<A>();
  std::vector<A*> ptrs = inputs_of(in);
  for (A* p : ptrs) {
    tape.registerInput(*p);
  }
  tape.newRecording();

  const std::size_t s0 = tape.getNumStatements(), o0 = tape.getNumOperations();
  const A collar = l.collar_at<A>(in);
  const std::size_t s1 = tape.getNumStatements(), o1 = tape.getNumOperations();
  pl::Leaf::LeafOutputs<A> got = l.outputs_at<A>(collar, in);
  const std::size_t s2 = tape.getNumStatements(), o2 = tape.getNumOperations();

  const std::size_t n_in = ptrs.size();
  const std::size_t n_out = 1 + got.uptake.size();

  // Where collar_at's statements go: marginal_at is the residual implicit_value
  // evaluates, and it is three calls. marginal_assembled reads A' and C' off the
  // slope primitives, so it is arithmetic at the adjoint and nothing above it.
  const A held = A(xad::value(collar));
  const std::size_t c0 = tape.getNumStatements();
  const pl::Leaf::CollarCoords<A> co =
      l.collar_coords_at<A>(l.opt_psi_stem_, l.ci_, held, in);
  const std::size_t c1 = tape.getNumStatements();
  const A dE = l.roots_.duptake_dpsi_at<A>(held, in.supply.at());
  const std::size_t c2 = tape.getNumStatements();
  const pl::Leaf::CollarPoint<A> at{held, co.sigma, co.ci, dE, co.flux};
  const A mm = l.marginal_assembled<A>(at, in.profit).marginal;
  const std::size_t c3 = tape.getNumStatements();
  (void)mm;

  // What a NESTING would cost, kept as the counterfactual the slope primitives
  // replaced. At TT = FReal<A>, a tangent above the adjoint, both halves of every
  // operation record and FReal assigns each separately, so the expression template
  // cannot fuse. The same kernels at the working scalar A are the same arithmetic.
  using TT = typename xad::fwd<A>::active_type;
  const auto lift = [](const A& v) { TT o{}; xad::value(o) = v; return o; };
  // WHAT A SUPPLIED-ROW DESIGN WOULD STILL RECORD. The soil state reaches the leaf
  // only through total uptake and its collar slope, so plant would tape the supply
  // -- closed-form Ohm's law over a tabulated integral, and the only part whose
  // inputs are numerous -- and graft the gas-exchange model's rows onto it.
  const std::size_t w0 = tape.getNumStatements();
  std::vector<A> per_layer;
  const A e_up = l.E_from_soil_at<A>(held, in.supply.at(), per_layer);
  const A s_up = l.roots_.template duptake_dpsi_at<A>(held, in.supply.at());
  const std::size_t w1 = tape.getNumStatements();
  std::printf("  THE SUPPLY ALONE (E and S)   %8zu statements   <- all a graft would tape\n",
              w1 - w0);
  (void)e_up; (void)s_up;

  const std::size_t k0 = tape.getNumStatements();
  const A jA = l.electron_transport_kernel<A>(in.profit.ppfd, in.profit.quantum_yield,
                                              in.profit.curv_elec, in.profit.transport_jmax);
  const A aA = l.assim_colimited_kernel<A>(co.ci, in.profit.vcmax, jA,
                                           in.profit.curv_colim, in.profit.respiration);
  const A cA = l.hydraulic_cost_TF_kernel<A>(co.sigma, in.profit.stem_b, in.profit.stem_c,
                                             in.profit.beta2, in.profit.cost_scale);
  const std::size_t k1 = tape.getNumStatements();
  const std::size_t q0 = tape.getNumStatements();
  const TT jT = l.electron_transport_kernel<TT>(lift(in.profit.ppfd), lift(in.profit.quantum_yield),
                                                lift(in.profit.curv_elec), lift(in.profit.transport_jmax));
  const std::size_t q1 = tape.getNumStatements();
  TT ciT = lift(co.ci); xad::derivative(ciT) = A(1.0);
  const TT aT = l.assim_colimited_kernel<TT>(ciT, lift(in.profit.vcmax), jT,
                                             lift(in.profit.curv_colim), lift(in.profit.respiration));
  const std::size_t q2 = tape.getNumStatements();
  TT sgT = lift(co.sigma); xad::derivative(sgT) = A(1.0);
  const TT cT = l.hydraulic_cost_TF_kernel<TT>(sgT, lift(in.profit.stem_b), lift(in.profit.stem_c),
                                               lift(in.profit.beta2), lift(in.profit.cost_scale));
  const std::size_t k2 = tape.getNumStatements();
  std::printf("    at FReal<A>, by kernel:  electron %zu  colimited %zu  cost %zu\n",
              q1 - q0, q2 - q1, k2 - q2);
  (void)aA; (void)cA; (void)aT; (void)cT;

  std::printf("  active inputs              %zu\n", n_in);
  std::printf("  outputs plant reads        %zu  (profit + %zu layer draws)\n",
              n_out, got.uptake.size());
  std::printf("  collar_at                  %8zu statements %9zu operations\n",
              s1 - s0, o1 - o0);
  std::printf("  outputs_at                 %8zu statements %9zu operations\n",
              s2 - s1, o2 - o1);
  std::printf("  ONE PLACEMENT              %8zu statements %9zu operations\n",
              s2 - s0, o2 - o0);
  std::printf("    inside marginal_at:\n");
  std::printf("      collar_coords_at         %8zu statements\n", c1 - c0);
  std::printf("      duptake_dpsi_at          %8zu statements\n", c2 - c1);
  std::printf("      marginal_assembled       %8zu statements  <- from the slope primitives\n",
              c3 - c2);
  std::printf("    the three kernels at A     %8zu statements\n", k1 - k0);
  std::printf("    the same at FReal<A>       %8zu statements  <- %.1fx\n",
              k2 - k1, double(k2 - k1) / double(k1 - k0));
  std::printf("  tape memory                %8zu bytes\n", tape.getMemory());

  // The alternative, in the same units. record_with_derivatives writes one
  // statement per non-zero row and skips the zeros, so a dense block is at most
  // this and typically less.
  std::printf("  a dense block would be     %8zu statements  (%zu outputs x %zu rows)\n",
              n_out * n_in, n_out, n_in);
  const double ratio = double(s2 - s0) / double(n_out * n_in);
  std::printf("  ratio, recorded : supplied %8.1fx\n", ratio);

  // Sweeping first: the recording it walks is the one just counted, and the
  // timing below clears it. Seeded on profit -- one output, one walk.
  tape.registerOutput(got.profit);
  const long reps = 400;
  clock_type::time_point t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    tape.clearDerivatives();
    xad::derivative(got.profit) = 1.0;
    tape.computeAdjoints();
  }
  const double us_sweep = us_per(t0, clock_type::now(), reps);

  // ⚠️ THE MARGINAL COST OF A PLACEMENT, which is what plant pays: it clears once
  // per STEP recording and records many placements into it, so charging a clear
  // and a re-registration to each one measures a tape cycle nothing performs.
  tape.clearAll();
  for (A* p : ptrs) { *p = A(xad::value(*p)); }
  for (A* p : ptrs) { tape.registerInput(*p); }
  tape.newRecording();
  t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    const A c = l.collar_at<A>(in);
    pl::Leaf::LeafOutputs<A> o = l.outputs_at<A>(c, in);
    (void)o;
  }
  const double us_record = us_per(t0, clock_type::now(), reps);

  // And the cycle a PRIVATE tape per placement would add on top: release,
  // register, clear, new recording. This is the overhead move 2 has to beat.
  t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    tape.clearAll();
    for (A* p : ptrs) { *p = A(xad::value(*p)); }
    for (A* p : ptrs) { tape.registerInput(*p); }
    tape.newRecording();
  }
  const double us_cycle = us_per(t0, clock_type::now(), reps);

  std::printf("  record one placement       %8.2f us   (marginal, no tape cycle)\n",
              us_record);
  std::printf("  sweep it once              %8.2f us\n", us_sweep);
  std::printf("  an empty tape cycle        %8.2f us   (clear+release+register)\n",
              us_cycle);
  std::printf("  record : sweep             %8.2fx\n", us_record / us_sweep);
}


// The OTHER nesting. FReal<AReal> is a tangent WRAPPING an adjoint, and it costs
// 18x because FReal assigns its two halves separately. AReal<FReal> is an adjoint
// whose SCALAR is a tangent -- the expression template sees one scalar type and
// fuses normally, and one sweep carries both d/dtheta and d2/dp dtheta.
//
// This is what collar_condition used and what increment 2 deleted.
void nested_cost(int layers) {
  using D = typename xad::fwd_adj<double>::active_type;   // AReal<FReal<double>>
  using DTape = typename xad::fwd_adj<double>::tape_type; // Tape<FReal<double>>

  pl::Leaf l;
  set_up(l, layers);
  l.find_root_collar_psi();
  const fixture::Soil soil = soil_of(layers);
  const pl::LeafInputs<double> in_d =
      fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);

  // The referee, taken FIRST because two tapes of one scalar type cannot both be
  // active: d2(profit)/dp d(vcmax) against a central difference of
  // d(profit)/d(vcmax) in the collar, each arm its own recording.
  auto row_at = [&](double p_at) -> double {
    DTape t2;
    pl::LeafInputs<D> q = in_d.template rebind_from<D>();
    std::vector<D*> qp = q.profit.field_ptrs();
    for (D& v : q.supply.psi_soil) qp.push_back(&v);
    for (D& v : q.supply.r_R_H_min) qp.push_back(&v);
    for (D& v : q.supply.r_R_V_sum) qp.push_back(&v);
    qp.push_back(&q.supply.root_b); qp.push_back(&q.supply.root_c);
    qp.push_back(&q.psi_crit); qp.push_back(&q.root_psi_crit);
    for (D* z : qp) { t2.registerInput(*z); }
    t2.newRecording();
    pl::Leaf::LeafOutputs<D> o = l.outputs_at<D>(D(p_at), q);
    t2.registerOutput(o.profit);
    t2.clearDerivatives();
    xad::derivative(o.profit) = 1.0;
    t2.computeAdjoints();
    return xad::value(xad::derivative(*qp[0]));
  };
  const double h = 1e-5 * std::max(1.0, std::abs(l.opt_root_psi_));
  const double fd_hi = row_at(l.opt_root_psi_ + h);
  const double fd_lo = row_at(l.opt_root_psi_ - h);
  const double fd = (fd_hi - fd_lo) / (2 * h);

  // Does the inner direction propagate AT ALL? Evaluate away from the optimum,
  // where dprofit/dp is not zero, and compare against a double difference.
  {
    const double pa = l.opt_root_psi_ * 0.97;
    DTape t3;
    pl::LeafInputs<D> q = in_d.template rebind_from<D>();
    D c3 = D(pa);
    xad::derivative(xad::value(c3)) = 1.0;
    t3.newRecording();
    pl::Leaf::LeafOutputs<D> o3 = l.outputs_at<D>(c3, q);
    const double tang = xad::derivative(xad::value(o3.profit));
    const double hh = 1e-6 * pa;
    const pl::LeafInputs<double> qd = in_d;
    const double up = l.outputs_at<double>(pa + hh, qd).profit;
    const double dn = l.outputs_at<double>(pa - hh, qd).profit;
    std::printf("  off-optimum: tangent       %.12g\n", tang);
    std::printf("  off-optimum: differenced   %.12g\n", (up - dn) / (2 * hh));
  }

  DTape tape;
  pl::LeafInputs<D> in = in_d.template rebind_from<D>();
  std::vector<D*> ptrs = in.profit.field_ptrs();
  for (D& v : in.supply.psi_soil) ptrs.push_back(&v);
  for (D& v : in.supply.r_R_H_min) ptrs.push_back(&v);
  for (D& v : in.supply.r_R_V_sum) ptrs.push_back(&v);
  ptrs.push_back(&in.supply.root_b);
  ptrs.push_back(&in.supply.root_c);
  ptrs.push_back(&in.psi_crit);
  ptrs.push_back(&in.root_psi_crit);
  for (D* q : ptrs) { tape.registerInput(*q); }
  tape.newRecording();

  // The collar is the DOUBLE the solve returned, carrying an inner direction of
  // one. Nothing here is solved on the tape.
  D collar = D(l.opt_root_psi_);
  xad::derivative(xad::value(collar)) = 1.0;

  const std::size_t s0 = tape.getNumStatements();
  pl::Leaf::LeafOutputs<D> got = l.outputs_at<D>(collar, in);
  const std::size_t s1 = tape.getNumStatements();

  std::printf("\n--- %d layers, outputs_at at AReal<FReal<double>> ---\n", layers);
  std::printf("  statements                 %8zu\n", s1 - s0);
  std::printf("  operations                 %8zu\n", tape.getNumOperations());

  // One sweep per output, read BEFORE any re-recording. The VALUE of an input's
  // adjoint is d(output)/d(theta); its DERIVATIVE is d2(output)/dp d(theta) --
  // both from the one recording, which is the whole point.
  tape.registerOutput(got.profit);
  tape.clearDerivatives();
  xad::derivative(got.profit) = 1.0;
  tape.computeAdjoints();
  const double d1 = xad::value(xad::derivative(*ptrs[0]));
  const double d2 = xad::derivative(xad::derivative(*ptrs[0]));
  // Stationarity, free from the same pass: the inner tangent of profit is
  // dprofit/dp, which the solve put at zero.
  const double resid = xad::derivative(xad::value(got.profit));



  const long reps = 400;
  clock_type::time_point t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    tape.clearDerivatives();
    xad::derivative(got.profit) = 1.0;
    tape.computeAdjoints();
  }
  const double us_sweep = us_per(t0, clock_type::now(), reps);

  t0 = clock_type::now();
  for (long r = 0; r < reps; ++r) {
    pl::Leaf::LeafOutputs<D> o = l.outputs_at<D>(collar, in);
    (void)o;
  }
  const double us_record = us_per(t0, clock_type::now(), reps);

  std::printf("  dprofit/dp at the optimum  %.3e   (stationarity, free here)\n", resid);
  std::printf("  d(profit)/d(vcmax)         %.12g\n", d1);
  std::printf("  d2(profit)/dp d(vcmax)     %.12g\n", d2);
  std::printf("  the same, differenced      %.12g   rel gap %.3e\n",
              fd, std::abs(d2 - fd) / (std::abs(fd) + 1e-300));
}

}  // namespace

int main() {
  std::printf("what one leaf placement costs the tape\n");
  for (int layers : {1, 3, 5}) {
    one_width(layers);
  }
  nested_cost(5);
  std::printf(
      "\nplant records one placement per cohort per stage per step and sweeps\n"
      "the result once per census metric.\n");
  return 0;
}
