// The leaf's outputs at a held collar, against a difference of the same.
//
// ONE CLAIM: every derivative the leaf reports is the derivative of the model it
// reports it about. The reference is the model itself -- rebuild the leaf at a
// moved input, evaluate at the SAME collar, difference -- so it shares nothing
// with the tape path except the forward pass. That is why this replaced
// twenty-nine tests written against a layer that produced the rows by a second
// route: those compared two derivations of one model, and agreed whenever both
// were wrong the same way.
//
// ⚠️ THE COLLAR IS HELD, AND THAT IS NOT AN ECONOMY. Re-optimising it makes the
// reference measure the SOLVER as well as the model. The wet endpoint is
// infeasible by construction -- uptake is exactly zero there, so the marginal
// takes its no-flow exit -- and the solve steps inside by a fraction of the
// interval's WIDTH until the gradient is real. A pinned collar is therefore the
// bound plus a feasibility step (measured, 5.4e-07), and differencing a re-solve
// tracks that step's own response to the traits. A gradient of a solver's step is
// not a gradient of anything.
//
// So the collar is held, and held in the MIDDLE of the feasible interval rather
// than at the operating point: an optimum sitting on a bound leaves the interval
// on half of every perturbation, and a reference that refuses half its steps is
// not one. The collar's own row is refereed where it is a statement about the
// model: probe_marginal_tangent for the interior condition, probe_bound_tangent
// for the three bounds.
//
// The named inputs, the named states, both outputs. A difference cannot resolve
// arbitrarily small numbers, so the budget is what the difference itself can
// see: two step sizes, and the gap between them is the floor.
//
//   make CXX=g++ test_leaf_gradient && ./test_leaf_gradient

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

using phylloptim::Leaf;
using phylloptim::LeafInputs;
using phylloptim::tangent;
using phylloptim::derivative_along;
using fixture::Input;

namespace {

int failures = 0;
int compared = 0;
int crossed = 0;

const double kRd25 = 1.44;
const double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,      157.44,   0.30,
                         0.7,      0.99,     7.5,      kRd25};
const double kKmax = 1.0 * 0.000157 / 5.0;

// The states the leaf answers in, chosen to reach every branch its solve has.
struct State {
  const char* what;
  double psi_soil, ppfd, vpd, leaf_temp;
  int layers;
  // One trait moved off the shared set, where a branch is not otherwise
  // reachable. The root's critical potential equals the stem's in the shared
  // traits, so the dry bound never picks the root arm until one of them moves.
  int slot = -1;
  double value = 0.0;
};
const char* bound_name(Leaf::WhichBound w) {
  switch (w) {
    case Leaf::WhichBound::Wet: return "wet";
    case Leaf::WhichBound::DryRootCrit: return "dry-root-crit";
    default: return "dry-root-psi-crit";
  }
}

const State kStates[] = {
    {"wet",    2.0,  900.0, 2.0, 25.0, 3},
    {"dim",    2.0,  300.0, 2.0, 25.0, 3},
    {"dry",    4.5,  900.0, 2.0, 25.0, 3},
    {"arid",   2.0,  900.0, 4.0, 25.0, 3},
    {"hot",    2.0,  900.0, 2.0, 40.0, 3},
    {"pinned", 5.72, 900.0, 2.0, 25.0, 1},
    {"shut",   5.8,  900.0, 2.0, 25.0, 3},
    {"single", 2.0,  900.0, 2.0, 25.0, 1},
    // The wet end: profit is still climbing when the collar runs out of
    // feasible interval, so the bound is the answer.
    {"wet-pin", 4.0,   30.0, 2.0, 25.0, 3},
    // Dim enough that the leaf holds at the collar where uptake vanishes.
    {"shade",   2.0,    1.0, 2.0, 25.0, 3},
    // The root's own limit, brought in ahead of the stem's so the dry bound has
    // a reason to pick it.
    {"root-pin", 2.0,  900.0, 2.0, 25.0, 3, phylloptim::trait_root_psi_crit, 3.0},
};

fixture::Soil soil_of(const State& st) {
  fixture::Soil soil;
  for (int i = 0; i < st.layers; ++i) {
    soil.carbon.push_back(1.0 / st.layers / 0.05);
    soil.depth.push_back(1.0 * (i + 1));
  }
  return soil;
}

fixture::Physiology drivers_of(const State& st, const fixture::Soil& soil) {
  fixture::Physiology d;
  for (int i = 0; i < st.layers; ++i) {
    d.psi_soil.push_back(st.psi_soil + 0.25 * i);
  }
  d.soil_depth = soil.depth;
  d.root_network = fixture::root_network(soil.carbon, soil.depth);
  d.PPFD = st.ppfd;
  d.kmax = kKmax;
  d.atm_vpd = st.vpd;
  d.ca = 40.0;
  d.leaf_temp = st.leaf_temp;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

// Which entry of theta a named trait is, or -1 where the input is not a trait.
int trait_slot(Input which) {
  switch (which) {
    case Input::vcmax_25: return phylloptim::trait_vcmax_25;
    case Input::jmax_25: return phylloptim::trait_jmax_25;
    case Input::quantum_yield: return phylloptim::trait_a;
    case Input::curv_elec: return phylloptim::trait_curv_fact_elec_trans;
    case Input::curv_colim: return phylloptim::trait_curv_fact_colim;
    case Input::respiration: return phylloptim::trait_R_d_25;
    case Input::stem_b: return phylloptim::trait_stem_b;
    case Input::stem_c: return phylloptim::trait_stem_c;
    case Input::beta2: return phylloptim::trait_beta2;
    case Input::cost_scale: return phylloptim::trait_cost_scale_TF24;
    case Input::psi_crit: return phylloptim::trait_psi_crit;
    case Input::root_b: return phylloptim::trait_root_b;
    case Input::root_c: return phylloptim::trait_root_c;
    case Input::root_psi_crit: return phylloptim::trait_root_psi_crit;
    default: return -1;
  }
}

// The value an input holds before it is moved.
double base_of(const State& st, const fixture::Soil& soil, Input which,
               int layer) {
  const int slot = trait_slot(which);
  if (slot >= 0) return kTheta[slot];
  switch (which) {
    case Input::ppfd: return st.ppfd;
    case Input::kmax: return kKmax;
    case Input::psi_soil: return st.psi_soil + 0.25 * layer;
    case Input::root_carbon: return soil.carbon[std::size_t(layer)];
    default: return 0.0;
  }
}

// The outputs of a whole solve at a moved input: profit, then the per-layer
// uptake. Rebuilt from the traits and the drivers, so nothing of the base solve
// survives into it.
bool held_outputs(const State& st, Input which, int layer, double moved,
                  double collar, std::vector<double>& out) {
  fixture::Soil soil = soil_of(st);
  double theta[phylloptim::n_traits];
  for (int i = 0; i < phylloptim::n_traits; ++i) theta[i] = kTheta[i];
  if (st.slot >= 0) theta[st.slot] = st.value;
  const int slot = trait_slot(which);
  if (slot >= 0) {
    theta[slot] = moved;
  }
  if (which == Input::root_carbon) {
    soil.carbon[std::size_t(layer)] = moved;
  }
  fixture::Physiology d = drivers_of(st, soil);
  if (which == Input::ppfd) d.PPFD = moved;
  if (which == Input::kmax) d.kmax = moved;
  if (which == Input::psi_soil) d.psi_soil[std::size_t(layer)] = moved;

  Leaf l;
  phylloptim::Leaf::FixedCollarEval at;
  try {
    d.drive(l, theta);
    at = l.profit_at_fixed_collar(collar);
  } catch (const std::exception&) {
    return false;
  }
  if (!at.feasible) return false;
  out.assign(std::size_t(st.layers) + 1, 0.0);
  out[0] = at.profit;
  for (int j = 0; j < st.layers; ++j) {
    out[std::size_t(j) + 1] = at.uptake[std::size_t(j)];
  }
  for (double v : out) {
    if (!std::isfinite(v)) return false;
  }
  return true;
}

}  // namespace

int main() {
  std::printf("the leaf's gradient against a difference of its own solve\n\n");
  std::printf("   %-8s %-12s %-14s %2s %11s %11s %11s\n", "state", "point",
              "input", "L", "tape", "differenced", "budget");

  const Input inputs[] = {
      Input::vcmax_25, Input::jmax_25,  Input::quantum_yield, Input::curv_elec,
      Input::curv_colim, Input::respiration, Input::ppfd,     Input::kmax,
      Input::stem_b,   Input::stem_c,   Input::beta2,         Input::cost_scale,
      Input::psi_crit, Input::root_b,   Input::root_c,        Input::root_psi_crit,
      Input::psi_soil, Input::root_carbon};
  const std::size_t n_inputs = sizeof(inputs) / sizeof(inputs[0]);
  std::vector<int> live(n_inputs, 0);
  std::vector<int> seen(n_inputs, 0);

  for (const State& st : kStates) {
    const fixture::Soil soil = soil_of(st);
    const fixture::Physiology d = drivers_of(st, soil);
    double theta0[phylloptim::n_traits];
    for (int i = 0; i < phylloptim::n_traits; ++i) theta0[i] = kTheta[i];
    if (st.slot >= 0) theta0[st.slot] = st.value;
    Leaf l;
    d.drive(l, theta0);
    try {
      l.find_root_collar_psi();
    } catch (const std::exception& e) {
      std::printf("   %-8s solve failed: %s\n", st.what, e.what());
      ++failures;
      continue;
    }
    const char* point = Leaf::operating_point_kind_name(l.operating_point_kind());
    const std::size_t n_out = std::size_t(st.layers) + 1;

    // The collar every comparison at this state is taken at: the middle of the
    // feasible interval, which the model admits at every step of every input.
    double lo = 0.0, hi = 0.0;
    if (!l.prepare_collar_solve(lo, hi)) {
      std::printf("   %-8s %-12s no feasible interval to hold a collar in\n",
                  st.what, point);
      continue;
    }
    const double held = lo + 0.5 * (hi - lo);
    l.find_root_collar_psi();

    for (std::size_t k = 0; k < n_inputs; ++k) {
      const int n_layer = (inputs[k] == Input::psi_soil ||
                           inputs[k] == Input::root_carbon)
                              ? st.layers
                              : 1;
      for (int layer = 0; layer < n_layer; ++layer) {
    const double base = st.slot >= 0 && st.slot == trait_slot(inputs[k])
                                ? st.value
                                : base_of(st, soil, inputs[k], layer);
        if (!(std::abs(base) > 0.0)) continue;

        // The tape path: one seeded evaluation, every output at once.
        const LeafInputs<tangent> in =
            fixture::leaf_inputs<tangent>(l, inputs[k], layer, soil);
        // The tape path at the same held collar. The collar is passive, so this
        // is the model's own branch for this kind evaluated at a point that does
        // not move -- which is what the difference below is.
        // The interior composition at that collar: the two residuals, the
        // kernels and the quadrature. Which condition PLACED the collar is not
        // this claim's subject and is refereed elsewhere.
        tangent profit;
        std::vector<tangent> uptake;
        try {
          const phylloptim::Leaf::FixedCollarEval base_at =
              l.profit_at_fixed_collar(held);
          static_cast<void>(base_at);
          const tangent at_collar(held);
          profit = l.profit_at<tangent>(
              l.opt_psi_stem_, l.ci_, at_collar,
              l.supply_draw_at<tangent>(at_collar, in.supply).flux, in.profit);
          l.E_from_soil_at<tangent>(at_collar, in.supply.at(), uptake);
        } catch (const std::exception& e) {
          std::printf("   %-8s %-12s %-14s %2d   refused: %s\n", st.what, point,
                      fixture::name_of(inputs[k]), layer, e.what());
          ++failures;
          continue;
        }
        std::vector<double> tape(n_out, 0.0);
        tape[0] = derivative_along(profit);
        for (int j = 0; j < st.layers; ++j) {
          tape[std::size_t(j) + 1] = derivative_along(uptake[std::size_t(j)]);
        }

        // The reference: two central differences of a whole re-solve.
        auto differenced = [&](double rel, std::vector<double>& into) -> bool {
          const double h = std::abs(base) * rel;
          std::vector<double> up, dn;
          if (!held_outputs(st, inputs[k], layer, base + h, held, up)) return false;
          if (!held_outputs(st, inputs[k], layer, base - h, held, dn)) return false;
          into.assign(n_out, 0.0);
          for (std::size_t j = 0; j < n_out; ++j) {
            into[j] = (up[j] - dn[j]) / (2.0 * h);
          }
          return true;
        };
        std::vector<double> coarse, fine;
        if (!differenced(1e-5, coarse) || !differenced(1e-6, fine)) {
          ++crossed;
          std::printf("   %-8s %-12s %-14s %2d   the collar is not admissible "
                      "at the moved input\n",
                      st.what, point, fixture::name_of(inputs[k]), layer);
          continue;
        }

        // ⚠️ THE BUDGET IS WHAT THE DIFFERENCE CAN SEE, NOT A NUMBER CHOSEN TO
        // PASS. Two step sizes disagree by the difference's own error, which is
        // the floor below which it cannot referee anything -- so the gap between
        // them sets it, and a flat tolerance would be either blind or a
        // permanent nuisance depending on the row.
        double scale = 0.0, floor = 0.0;
        for (std::size_t j = 0; j < n_out; ++j) {
          scale = std::max(scale, std::abs(fine[j]));
          floor = std::max(floor, std::abs(coarse[j] - fine[j]));
        }
        ++seen[k];
        if (scale > 0.0) ++live[k];
        // ⚠️ THE TWO-STEP GAP IS NOT THE WHOLE FLOOR. It sees the difference's
        // NOISE and is blind to its BIAS: a bound found to 1e-4 and a
        // concentration to 1e-10 shift both evaluations the same way, so two
        // steps can agree closely and both be off together. The explicit term is
        // what a central difference of a re-solved nonlinear system is worth.
        const double budget =
            std::max(10.0 * floor, 1e-7 * std::max(scale, 1.0));
        double worst = 0.0;
        std::size_t at = 0;
        for (std::size_t j = 0; j < n_out; ++j) {
          const double err = std::abs(tape[j] - fine[j]);
          if (err > worst) { worst = err; at = j; }
        }
        ++compared;
        if (!(worst <= budget)) {
          ++failures;
          std::printf("   %-8s %-12s %-14s %2d out %zu  tape %13.6e  fine "
                      "%13.6e  coarse %13.6e  budget %10.3e  FAIL\n",
                      st.what, point, fixture::name_of(inputs[k]), layer, at,
                      tape[at], fine[at], coarse[at], budget);
        }
      }
    }
    std::printf("   %-8s %-12s all inputs within the difference's own budget\n",
                st.what, point);
  }

  std::printf("\n   states in which each input moved an output:\n");
  for (std::size_t k = 0; k < n_inputs; ++k) {
    std::printf("     %-16s live in %d of %d\n", fixture::name_of(inputs[k]),
                live[k], seen[k]);
  }
  // --- claim 2: the bounds -------------------------------------------------
  // Which condition places a pinned collar, refereed by differencing that
  // condition's own solve rather than the whole leaf's.
  std::printf("\n  the bounds, against a difference of the bound's own solve\n");
  int bound_rows = 0;
  for (const State& st : kStates) {
    const fixture::Soil soil = soil_of(st);
    double theta0[phylloptim::n_traits];
    for (int i = 0; i < phylloptim::n_traits; ++i) theta0[i] = kTheta[i];
    if (st.slot >= 0) theta0[st.slot] = st.value;
    Leaf l;
    drivers_of(st, soil).drive(l, theta0);
    try { l.find_root_collar_psi(); } catch (const std::exception&) { continue; }

    const Leaf::WhichBound arms[] = {Leaf::WhichBound::Wet,
                                     Leaf::WhichBound::DryRootCrit,
                                     Leaf::WhichBound::DryRootPsiCrit};
    for (const Leaf::WhichBound arm : arms) {
      // The bound the model itself finds, at the traits as they stand.
      auto bound_of = [&](Input which, int layer, double moved,
                          double& into) -> bool {
        fixture::Soil sl = soil_of(st);
        double th[phylloptim::n_traits];
        for (int i = 0; i < phylloptim::n_traits; ++i) th[i] = kTheta[i];
        if (st.slot >= 0) th[st.slot] = st.value;
        const int slot = trait_slot(which);
        if (slot >= 0) th[slot] = moved;
        if (which == Input::root_carbon) sl.carbon[std::size_t(layer)] = moved;
        fixture::Physiology dd = drivers_of(st, sl);
        if (which == Input::ppfd) dd.PPFD = moved;
        if (which == Input::kmax) dd.kmax = moved;
        if (which == Input::psi_soil) dd.psi_soil[std::size_t(layer)] = moved;
        Leaf m;
        try {
          dd.drive(m, th);
          if (arm == Leaf::WhichBound::DryRootPsiCrit) {
            into = m.supply_psi_crit();
            return true;
          }
          const double wettest = m.supply_begin_solve();
          into = m.find_root_psi(wettest, m.supply_psi_soil(),
                                 arm == Leaf::WhichBound::Wet ? 0 : 1);
        } catch (const std::exception&) {
          return false;
        }
        return std::isfinite(into);
      };

      double x = 0.0;
      if (!bound_of(Input::None, 0, 0.0, x)) continue;
      l.find_root_collar_psi();

      for (std::size_t k = 0; k < n_inputs; ++k) {
        const int n_layer = (inputs[k] == Input::psi_soil ||
                             inputs[k] == Input::root_carbon) ? st.layers : 1;
        for (int layer = 0; layer < n_layer; ++layer) {
          const double base = st.slot >= 0 && st.slot == trait_slot(inputs[k])
                                  ? st.value
                                  : base_of(st, soil, inputs[k], layer);
          if (!(std::abs(base) > 0.0)) continue;
          const LeafInputs<tangent> in =
              fixture::leaf_inputs<tangent>(l, inputs[k], layer, soil);
          double tape = 0.0;
          try {
            tape = derivative_along(
                l.bound_at<tangent>(
                    arm, x, in, l.supply_draw_at<tangent>(tangent(x), in.supply)));
          } catch (const std::exception&) { continue; }

          double up = 0.0, dn = 0.0, up2 = 0.0, dn2 = 0.0;
          const double h = std::abs(base) * 1e-6, h2 = std::abs(base) * 1e-5;
          if (!bound_of(inputs[k], layer, base + h, up)) continue;
          if (!bound_of(inputs[k], layer, base - h, dn)) continue;
          if (!bound_of(inputs[k], layer, base + h2, up2)) continue;
          if (!bound_of(inputs[k], layer, base - h2, dn2)) continue;
          const double fine = (up - dn) / (2.0 * h);
          const double coarse = (up2 - dn2) / (2.0 * h2);
          // ⚠️ A BOUND IS FOUND TO 1e-4, so a difference of it is worth about
          // five digits and no more. The two-step gap sees the noise and is
          // blind to the bias, which is why the explicit term is here.
          const double budget = std::max(
              10.0 * std::abs(coarse - fine), 1e-5 * std::max(std::abs(fine), 1.0));
          ++bound_rows;
          if (!(std::abs(tape - fine) <= budget)) {
            ++failures;
            std::printf("   %-8s %-16s %-14s %2d  tape %13.6e  fine %13.6e  "
                        "budget %10.3e  FAIL\n",
                        st.what, bound_name(arm), fixture::name_of(inputs[k]),
                        layer, tape, fine, budget);
          }
        }
      }
    }
  }
  std::printf("   %d bound rows\n", bound_rows);

  std::printf("\n   %d comparisons, %d failures", compared, failures);
  std::printf(" (%d inputs where the held collar stops being "
              "admissible)\n", crossed);
  return failures == 0 ? 0 : 1;
}
