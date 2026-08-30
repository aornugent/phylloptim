// PROBE: can a second tape run to completion INSIDE an outer recording, if the
// outer one is paused rather than nested?
//
// The implicit function theorem does not want a nested tape. It wants a SEQUENCED
// one: pause the stand's recording, run the leaf's own tape to completion, sweep it
// for the block, put the stand's recording back, and splice. `TapeAlreadyActive` is
// about two tapes being active at once, which this never does.
//
// This asks whether XAD permits it, and then prices the three ways of getting a
// leaf-sized block onto the stand's tape against each other in the one unit that
// decides it: STATEMENT WALKS per placement, counting a record and a sweep alike.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include \
//        ODELIA_SRC=../../../odelia/src probe_two_tapes && ./probe_two_tapes

#include <phylloptim.hpp>

#include <odelia/ode_interface.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using A = odelia::ode::active_scalar<double>;
using Tape = odelia::ode::adjoint_tape<double>;

namespace {

int failures = 0;

void ok(const char* what, bool cond) {
  std::printf("  %-58s %s\n", what, cond ? "ok" : "FAIL");
  if (!cond) ++failures;
}

}  // namespace

int main() {
  std::printf("a second tape inside an outer recording, by pausing rather than"
              " nesting\n\n");

  // The outer recording, as the stand's is: inputs registered, work recorded.
  Tape outer;
  A x = 3.0, z = 5.0;
  outer.registerInput(x);
  outer.registerInput(z);
  outer.newRecording();
  A y = x * z;                       // dy/dx = 5, dy/dz = 3

  // ⚠️ PAUSE, DO NOT NEST. deactivate() clears the active tape; setActive only
  // throws where a DIFFERENT tape is already active, so with none active the inner
  // one starts cleanly.
  const std::size_t outer_before = outer.getNumStatements();
  outer.deactivate();
  ok("the outer tape deactivates while holding a recording",
     !outer.isActive() && outer.getNumStatements() == outer_before);

  double d_inner_da = 0.0, d_inner_db = 0.0, inner_value = 0.0;
  {
    // The leaf's own tape: its inputs are the outer values STRIPPED to doubles, so
    // nothing here can touch a slot the outer tape owns.
    Tape inner;
    A a = odelia::util::to_passive(x);
    A b = odelia::util::to_passive(z);
    inner.registerInput(a);
    inner.registerInput(b);
    inner.newRecording();
    A w = a * a * b + sin(b);        // dw/da = 2ab, dw/db = a^2 + cos(b)
    inner.registerOutput(w);
    inner_value = odelia::util::to_passive(w);
    xad::derivative(w) = 1.0;
    inner.computeAdjoints();
    d_inner_da = xad::derivative(a);
    d_inner_db = xad::derivative(b);
    inner.deactivate();
  }

  const double want_da = 2.0 * 3.0 * 5.0;                 // 30
  const double want_db = 3.0 * 3.0 + std::cos(5.0);
  ok("the inner tape sweeps to the right rows",
     std::abs(d_inner_da - want_da) < 1e-12 &&
         std::abs(d_inner_db - want_db) < 1e-12);

  // Resume. The outer recording must still be there, and still be the same one.
  outer.activate();
  ok("the outer tape resumes with its recording intact",
     outer.isActive() && outer.getNumStatements() == outer_before);

  // Splice the block: the inner value carrying the inner rows, against the OUTER
  // actives. This is what a supplied-row boundary does.
  A spliced;
  const odelia::record_report rep = odelia::record_with_derivatives<A>(
      inner_value, {{x, d_inner_da}, {z, d_inner_db}}, spliced);
  ok("the block splices onto the outer tape", rep.whole);

  // One outer sweep over the whole thing.
  A out = y + spliced;
  outer.registerOutput(out);
  xad::derivative(out) = 1.0;
  outer.computeAdjoints();

  // d(out)/dx = z + dw/da = 5 + 30 ; d(out)/dz = x + dw/db = 3 + 9 + cos(5)
  const double gx = xad::derivative(x), gz = xad::derivative(z);
  std::printf("\n  d/dx = %.12g   (want %.12g)\n", gx, 5.0 + want_da);
  std::printf("  d/dz = %.12g   (want %.12g)\n", gz, 3.0 + want_db);
  ok("the outer sweep carries both the recorded and the spliced halves",
     std::abs(gx - (5.0 + want_da)) < 1e-12 &&
         std::abs(gz - (3.0 + want_db)) < 1e-12);

  // ⚠️ AND THE COST MODEL, which is what actually decides the design. A statement
  // walked is a statement walked, whether it is being recorded or swept. Per
  // placement, with T the leaf's statements, k the stand's seeds and m the leaf's
  // outputs:
  //
  //   recorded inline        T + k*T
  //   preaccumulated         T + m*T  + m + k*m       (record inner, sweep it m times)
  //   analytic block             0    + m + k*m       (no recording at all)
  const double T = 360, k = 3, m = 6;
  const double inline_walks = T + k * T;
  const double preacc_walks = T + m * T + m + k * m;
  const double analytic_walks = m + k * m;
  std::printf("\n  statement walks per placement, T=%.0f leaf statements,"
              " k=%.0f stand seeds, m=%.0f outputs\n", T, k, m);
  std::printf("    recorded inline      %8.0f\n", inline_walks);
  std::printf("    preaccumulated       %8.0f   %.2fx\n", preacc_walks,
              preacc_walks / inline_walks);
  std::printf("    analytic block       %8.0f   %.3fx\n", analytic_walks,
              analytic_walks / inline_walks);
  std::printf("\n  Preaccumulation loses because the stand sweeps k=%.0f times and\n"
              "  the block costs m=%.0f inner sweeps to build: it spends %.0f to save\n"
              "  %.0f. It would win only where m < k.\n", k, m, m, k);

  std::printf("\n%s\n", failures == 0 ? "all ok" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
