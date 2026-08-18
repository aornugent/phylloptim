// PROBE (not a test, not shipped): where does `gradient::rows_at` return
// non-finite entries, and per operating-point kind, is each one a genuine "no row
// exists" or an unfilled slot that an identity would fill exactly?
//
// Asks, at every point of the golden grid plus some deliberately off-grid states,
// for EVERY legal input and EVERY output (profit as Objective, collar as Point,
// the rest Ordinary), and tabulates the non-finiteness per kind.
//
//   make CXX=g++ probe_rows_na && stdbuf -oL ./probe_rows_na

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace grad = phylloptim::gradient;
using Kind = grad::OperatingPointKind;

namespace {

// test_leaf.cpp's env::kTheta, which is plant's TF24 default trait vector.
const double kRd25 = 1.44;  // test_leaf.cpp:66, and Leaf's own default
const double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,      157.44,   0.30,
                         0.7,      0.99,     7.5,      kRd25,
                         1.0 * 0.000157 / 5.0, 1e3};
static_assert(sizeof(kTheta) / sizeof(kTheta[0]) == grad::n_pars);

const double kAreaLeaf = 0.05;

// The golden grid's soil: `layers` equal 1 m layers drying with depth, carbon
// split evenly over all of them (so every layer is rooted and every root-carbon
// row is answerable).
grad::Drivers drivers(double psi_soil, double ppfd, double vpd, int layers,
                      double leaf_temp) {
  grad::Drivers d;
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers / kAreaLeaf;
  }
  d.root_network = fixture::root_network(root, depth);
  d.PPFD = ppfd;
  d.psi_soil = ps;
  d.soil_depth = depth;
  d.atm_vpd = vpd;
  d.ca = 40.0;
  d.leaf_temp = leaf_temp;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

phylloptim::Leaf fresh() {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  return l;
}

bool nf(double x) { return !std::isfinite(x); }

struct Pair {
  long seen = 0, nonfinite = 0, zero = 0, nonzero_finite = 0;
};

struct KindStats {
  long n_points = 0;
  long n_threw = 0;               // whole-request rows_at threw
  long n_slope_finite = 0;
  long n_point_finite = 0;
  std::map<std::string, long> messages;

  std::map<std::string, Pair> dydp;   // by output name
  std::map<std::string, Pair> dres;   // by input name
  long held_total = 0, held_nonfinite = 0;
  std::map<std::pair<std::string, std::string>, long> held_na;  // (out,in)

  // Point-level cross-tab: is the point ACTIVE (some finite non-zero dresidual)
  // where a dy_dp is missing?
  long n_any_dydp_na = 0;
  long n_ordinary_dydp_na = 0;   // == collar_channel returned false
  long n_objective_dydp_na = 0;
  long n_active = 0;             // some dresidual finite and != 0
  long n_dres_all_na = 0;
  long n_dres_all_zero_or_na = 0;
  long missing_active = 0, missing_passive = 0, ok_active = 0, ok_passive = 0;

  // At a pin: the mix of inputs whose rows FOLLOWED the point against those
  // taken at a held one.
  long n_pin_mix = 0, n_pin_follow_only = 0, n_pin_held_only = 0;
  long pin_follow_inputs = 0, pin_held_inputs = 0;

  // Does the IDENTITY substitution change any assembled total? Branch assembly
  // is the one test_rows_in_parts_assemble_to_the_totals uses; identity assembly
  // replaces NA slope by 1, NA dresidual by 0, NA dy_dp by 0 and never branches.
  long assembled = 0, identity_differs = 0, identity_new_nan = 0;
  // Is the POINT's own row non-zero where the parts call it passive? i.e. does
  // the collar really move with the inputs, with the movement inside `held`?
  long collar_held_nonzero = 0, collar_held_seen = 0;
  long dres_negative_zero = 0;
};

std::map<Kind, KindStats> stats;

// --- the mix at a pin ---------------------------------------------------------
//
// `rows_at` decides `follows` from whether the bound moves in that input, and
// then writes dresidual = 0 EITHER WAY (0 for a follower, -dpoint*slope = -0 for
// an input the bound does not read). So the mix is not readable off `Rows`; it is
// recomputed here from the same two functions rows_at uses.
void pin_mix(Kind kind, const grad::Drivers& d, const std::vector<int>& inputs,
             int n_layers, KindStats& st) {
  phylloptim::Leaf::WhichBound bound = phylloptim::Leaf::WhichBound::Wet;
  if (!grad::pinned_bound(kind, bound)) {
    return;
  }
  phylloptim::Leaf l = fresh();
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  l.find_root_collar_psi();
  const phylloptim::Leaf::BoundRow row = l.bound_row(bound);
  if (!row.finite) {
    return;
  }
  long follow = 0, held = 0;
  for (int p : inputs) {
    double dpoint = 0.0;
    if (grad::bound_moves_by_re_solving(p, false)) {
      bool at_base = false;
      grad::Scratch scratch;
      dpoint = grad::differenced_bound(l, kTheta, d, false, p, bound, s, at_base,
                                       scratch);
      // Put the leaf back where it was: differenced_bound leaves it perturbed.
      grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
      l.find_root_collar_psi();
    } else {
      dpoint = grad::bound_dpoint(row, p, n_layers);
    }
    if (dpoint != 0.0) {
      ++follow;
    } else {
      ++held;
    }
  }
  st.pin_follow_inputs += follow;
  st.pin_held_inputs += held;
  if (follow > 0 && held > 0) {
    ++st.n_pin_mix;
  } else if (follow > 0) {
    ++st.n_pin_follow_only;
  } else {
    ++st.n_pin_held_only;
  }
}

void tabulate(const grad::Rows& rows, const std::vector<int>& inputs,
              const std::vector<int>& outputs,
              const std::vector<grad::Role>& roles, int n_layers,
              const grad::Drivers& d) {
  KindStats& st = stats[rows.kind];
  ++st.n_points;
  st.n_slope_finite += !nf(rows.residual_slope);
  st.n_point_finite += !nf(rows.point);
  st.messages[rows.message.empty() ? "(none)" : rows.message] += 1;

  const std::size_t n = inputs.size();
  bool any_dydp_na = false, ordinary_na = false, objective_na = false;
  for (std::size_t j = 0; j < outputs.size(); ++j) {
    Pair& p = st.dydp[grad::output_name(outputs[j], n_layers)];
    ++p.seen;
    if (nf(rows.dy_dp[j])) {
      ++p.nonfinite;
      any_dydp_na = true;
      if (roles[j] == grad::Role::Ordinary) ordinary_na = true;
      if (roles[j] == grad::Role::Objective) objective_na = true;
    }
  }
  st.n_any_dydp_na += any_dydp_na;
  st.n_ordinary_dydp_na += ordinary_na;
  st.n_objective_dydp_na += objective_na;

  bool active = false;
  long n_dres_na = 0;
  for (std::size_t i = 0; i < n; ++i) {
    Pair& p = st.dres[grad::par_name(inputs[i], n_layers)];
    ++p.seen;
    const double v = rows.dresidual[i];
    if (nf(v)) {
      ++p.nonfinite;
      ++n_dres_na;
    } else if (v == 0.0) {
      ++p.zero;
    } else {
      ++p.nonzero_finite;
      active = true;
    }
  }
  st.n_active += active;
  st.n_dres_all_na += (n_dres_na == long(n));
  st.n_dres_all_zero_or_na += !active;
  if (ordinary_na || objective_na) {
    active ? ++st.missing_active : ++st.missing_passive;
  } else {
    active ? ++st.ok_active : ++st.ok_passive;
  }

  for (std::size_t i = 0; i < n; ++i) {
    const double v = rows.dresidual[i];
    if (v == 0.0 && std::signbit(v)) ++st.dres_negative_zero;
  }

  // The two assemblies, entry by entry.
  const double slope_id =
      std::isfinite(rows.residual_slope) ? rows.residual_slope : 1.0;
  for (std::size_t j = 0; j < outputs.size(); ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      ++st.held_total;
      const double held = rows.held[j * n + i];
      if (nf(held)) {
        ++st.held_nonfinite;
        st.held_na[{grad::output_name(outputs[j], n_layers),
                    grad::par_name(inputs[i], n_layers)}] += 1;
      }
      if (outputs[j] == grad::out_collar) {
        ++st.collar_held_seen;
        if (std::isfinite(held) && held != 0.0) ++st.collar_held_nonzero;
      }
      // Branch assembly (the existing test's).
      const double dpoint = -rows.dresidual[i] / rows.residual_slope;
      const double branch = (!std::isfinite(dpoint) || dpoint == 0.0)
                                ? held
                                : held + rows.dy_dp[j] * dpoint;
      // Identity assembly: no branch on anything.
      const double dres_id = std::isfinite(rows.dresidual[i]) ? rows.dresidual[i] : 0.0;
      const double dydp_id = std::isfinite(rows.dy_dp[j]) ? rows.dy_dp[j] : 0.0;
      const double id = held + dydp_id * (-dres_id / slope_id);
      ++st.assembled;
      const bool same = (std::isnan(branch) && std::isnan(id)) || branch == id;
      if (!same) ++st.identity_differs;
      if (std::isfinite(held) && !std::isfinite(id)) ++st.identity_new_nan;
    }
  }
  pin_mix(rows.kind, d, inputs, n_layers, st);

  // One worked example per kind: the whole dy_dp vector, and the collar's own
  // held row -- which is where the point's movement lives wherever the parts
  // report the point as passive.
  static std::set<Kind> dumped;
  if (dumped.insert(rows.kind).second) {
    printf("  EXAMPLE %s: slope=%g point=%g\n",
           phylloptim::Leaf::operating_point_kind_name(rows.kind),
           rows.residual_slope, rows.point);
    printf("    dy_dp:");
    for (std::size_t j = 0; j < outputs.size(); ++j) {
      printf(" %s=%g", grad::output_name(outputs[j], n_layers).c_str(),
             rows.dy_dp[j]);
    }
    printf("\n    held[collar][i] (non-zero only):");
    for (std::size_t i = 0; i < n; ++i) {
      const double h = rows.held[std::size_t(grad::out_collar) * n + i];
      if (!(h == 0.0)) {
        printf(" %s=%g", grad::par_name(inputs[i], n_layers).c_str(), h);
      }
    }
    printf("\n    dresidual (non-zero only):");
    for (std::size_t i = 0; i < n; ++i) {
      if (!(rows.dresidual[i] == 0.0)) {
        printf(" %s=%g", grad::par_name(inputs[i], n_layers).c_str(),
               rows.dresidual[i]);
      }
    }
    printf("\n");
  }
}

std::map<std::string, long> throw_messages;

// Cross-check against the golden file's published classification: the kinds the
// 288-point state grid falls into, per leaf temperature.
bool in_golden = false;
std::map<std::pair<double, Kind>, long> golden_kinds;

void one_point(double psi_soil, double ppfd, double vpd, int layers,
               double leaf_temp) {
  const grad::Drivers d = drivers(psi_soil, ppfd, vpd, layers, leaf_temp);
  const grad::Settings s;
  const int n_layers = layers;

  std::vector<int> inputs;
  for (int p = 0; p < grad::n_pars_total(n_layers); ++p) {
    inputs.push_back(p);  // nothing is refused by name on the multi-layer path
  }
  std::vector<int> outputs;
  std::vector<grad::Role> roles;
  for (int j = 0; j < grad::n_outputs_total(n_layers); ++j) {
    outputs.push_back(j);
    roles.push_back(j == grad::out_collar    ? grad::Role::Point
                    : j == grad::out_profit  ? grad::Role::Objective
                                             : grad::Role::Ordinary);
  }
  grad::RowRequest req;
  req.output = outputs.data();
  req.n_output = outputs.size();
  req.input = inputs.data();
  req.n_input = inputs.size();

  phylloptim::Leaf l = fresh();
  try {
    const grad::Rows rows = grad::rows_at(l, kTheta, d, req, s);
    if (in_golden) golden_kinds[{leaf_temp, rows.kind}] += 1;
    if (!rows.message.empty() &&
        rows.message.find("shade-death") == std::string::npos) {
      printf("  MESSAGE at psi_soil=%g ppfd=%g vpd=%g L=%d T=%g (%s): %s\n",
             psi_soil, ppfd, vpd, layers, leaf_temp,
             phylloptim::Leaf::operating_point_kind_name(rows.kind),
             rows.message.c_str());
    }
    tabulate(rows, inputs, outputs, roles, n_layers, d);
  } catch (const std::exception& e) {
    // Which kind was it? Solve it plainly, and record the refusal.
    phylloptim::Leaf l2 = fresh();
    grad::apply(l2, kTheta, d, false, -1, s.fast_stem_curve);
    l2.find_root_collar_psi();
    ++stats[l2.operating_point_kind()].n_threw;
    ++stats[l2.operating_point_kind()].n_points;
    throw_messages[e.what()] += 1;
    printf("  THREW at psi_soil=%g ppfd=%g vpd=%g L=%d T=%g (%s): %s\n",
           psi_soil, ppfd, vpd, layers, leaf_temp,
           phylloptim::Leaf::operating_point_kind_name(
               l2.operating_point_kind()),
           e.what());
  }
}

// A profile whose deepest layers hold NO root carbon: the one per-input refusal
// `rows_at` raises by name. Reported on its own, not folded into the table.
void unrooted(double psi_soil, double ppfd, double vpd, int layers, int rooted,
              double leaf_temp) {
  grad::Drivers d;
  std::vector<double> ps(layers), depth(layers), root(layers, 0.0);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    if (i < rooted) root[i] = 1.0 / rooted / kAreaLeaf;
  }
  d.root_network = fixture::root_network(root, depth);
  d.PPFD = ppfd;
  d.psi_soil = ps;
  d.soil_depth = depth;
  d.atm_vpd = vpd;
  d.ca = 40.0;
  d.leaf_temp = leaf_temp;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;

  std::vector<int> inputs;
  for (int p = 0; p < grad::n_pars_total(layers); ++p) inputs.push_back(p);
  std::vector<int> outputs;
  std::vector<grad::Role> roles;
  for (int j = 0; j < grad::n_outputs_total(layers); ++j) {
    outputs.push_back(j);
    roles.push_back(j == grad::out_collar   ? grad::Role::Point
                    : j == grad::out_profit ? grad::Role::Objective
                                            : grad::Role::Ordinary);
  }
  grad::RowRequest req{outputs.data(), outputs.size(), inputs.data(),
                       inputs.size()};
  phylloptim::Leaf l = fresh();
  const grad::Settings s;
  const grad::Rows rows = grad::rows_at(l, kTheta, d, req, s);
  std::map<std::string, long> na_by_input;
  for (std::size_t j = 0; j < outputs.size(); ++j) {
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      if (nf(rows.held[j * inputs.size() + i])) {
        na_by_input[grad::par_name(inputs[i], layers)] += 1;
      }
    }
  }
  printf("  psi_soil=%g ppfd=%g L=%d rooted=%d T=%g -> %s; held NA at:",
         psi_soil, ppfd, layers, rooted, leaf_temp,
         phylloptim::Leaf::operating_point_kind_name(rows.kind));
  for (const auto& kv : na_by_input) {
    printf(" %s(%ld)", kv.first.c_str(), kv.second);
  }
  printf("\n    dresidual NA at:");
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    if (nf(rows.dresidual[i])) {
      printf(" %s", grad::par_name(inputs[i], layers).c_str());
    }
  }
  printf("\n    message: %s\n",
         rows.message.empty() ? "(none)" : rows.message.c_str());
}

void report_kind(Kind k, const KindStats& st) {
  printf("\n=== %s: %ld point(s)%s\n",
         phylloptim::Leaf::operating_point_kind_name(k), st.n_points,
         st.n_threw ? " (some threw)" : "");
  printf("  residual_slope finite %ld/%ld   point finite %ld/%ld\n",
         st.n_slope_finite, st.n_points, st.n_point_finite, st.n_points);
  printf("  held entries non-finite %ld/%ld\n", st.held_nonfinite,
         st.held_total);
  if (!st.held_na.empty()) {
    printf("    non-finite (output, input) pairs:\n");
    // Collapse per-layer names so the table stays short.
    std::map<std::string, long> by_input;
    std::set<std::string> outs;
    for (const auto& kv : st.held_na) {
      by_input[kv.second > 0 ? kv.first.second : kv.first.second] += kv.second;
      outs.insert(kv.first.first);
    }
    printf("      outputs affected:");
    for (const std::string& o : outs) printf(" %s", o.c_str());
    printf("\n");
    for (const auto& kv : by_input) {
      printf("      input %-22s %ld entr(ies) over all outputs\n",
             kv.first.c_str(), kv.second);
    }
  }
  printf("  dy_dp non-finite, by output:\n");
  for (const auto& kv : st.dydp) {
    if (kv.second.nonfinite) {
      printf("      %-12s %ld/%ld non-finite\n", kv.first.c_str(),
             kv.second.nonfinite, kv.second.seen);
    }
  }
  bool any = false;
  for (const auto& kv : st.dydp) any = any || kv.second.nonfinite;
  if (!any) printf("      (none)\n");
  printf("  dresidual, by input: non-finite / exactly zero / finite non-zero\n");
  {
    long na = 0, z = 0, nz = 0, seen = 0;
    std::map<std::string, long> na_by_input, nz_by_input;
    for (const auto& kv : st.dres) {
      na += kv.second.nonfinite;
      z += kv.second.zero;
      nz += kv.second.nonzero_finite;
      seen += kv.second.seen;
      if (kv.second.nonfinite) na_by_input[kv.first] = kv.second.nonfinite;
      if (kv.second.nonzero_finite) nz_by_input[kv.first] = kv.second.nonzero_finite;
    }
    printf("      totals: %ld non-finite, %ld zero, %ld finite non-zero (of %ld)\n",
           na, z, nz, seen);
    if (!na_by_input.empty()) {
      printf("      non-finite at:");
      for (const auto& kv : na_by_input)
        printf(" %s(%ld)", kv.first.c_str(), kv.second);
      printf("\n");
    }
    if (!nz_by_input.empty() && nz_by_input.size() <= 30) {
      printf("      finite non-zero at:");
      for (const auto& kv : nz_by_input)
        printf(" %s(%ld)", kv.first.c_str(), kv.second);
      printf("\n");
    }
  }
  printf("  point-level: any dy_dp NA %ld, ordinary NA (collar_channel false) "
         "%ld, objective NA %ld\n",
         st.n_any_dydp_na, st.n_ordinary_dydp_na, st.n_objective_dydp_na);
  printf("             active (some dresidual finite non-zero) %ld, all "
         "dresidual NA %ld, passive %ld\n",
         st.n_active, st.n_dres_all_na, st.n_dres_all_zero_or_na);
  printf("             CROSS-TAB  dy_dp missing & ACTIVE %ld | missing & "
         "PASSIVE %ld | complete & active %ld | complete & passive %ld\n",
         st.missing_active, st.missing_passive, st.ok_active, st.ok_passive);
  if (st.n_pin_mix || st.n_pin_follow_only || st.n_pin_held_only) {
    printf("  pin mix: both kinds of input at one point %ld, followers only "
           "%ld, held only %ld  (inputs: %ld following, %ld held)\n",
           st.n_pin_mix, st.n_pin_follow_only, st.n_pin_held_only,
           st.pin_follow_inputs, st.pin_held_inputs);
  }
  printf("  identity vs branch assembly: %ld entries, %ld differ, %ld turned a "
         "finite held into NaN\n",
         st.assembled, st.identity_differs, st.identity_new_nan);
  printf("  collar's own held row non-zero %ld/%ld (the point moves, and the "
         "movement is inside `held`)\n",
         st.collar_held_nonzero, st.collar_held_seen);
  printf("  dresidual entries that are NEGATIVE zero: %ld\n",
         st.dres_negative_zero);
  printf("  messages:\n");
  for (const auto& kv : st.messages) {
    printf("      [%ld x] %s\n", kv.second, kv.first.c_str());
  }
}

}  // namespace

int main() {
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  const int layer_counts[] = {1, 3, 5};
  const double temps[] = {25.0, 40.0};

  printf("golden grid\n");
  in_golden = true;
  for (double t : temps)
    for (double p : psi_soils)
      for (double q : ppfds)
        for (double v : vpds)
          for (int n : layer_counts) one_point(p, q, v, n, t);

  in_golden = false;
  printf("\ngolden-grid classification (cross-check against test_golden):\n");
  for (const auto& kv : golden_kinds) {
    printf("  T=%4.1f  %-24s %ld\n", kv.first.first,
           phylloptim::Leaf::operating_point_kind_name(kv.first.second),
           kv.second);
  }

  // Off-grid states, hunting the kinds the grid does not reach.
  printf("off-grid sweeps\n");
  // Shade death: dim enough that gross assimilation cannot cover R_d.
  for (double q : {2.0, 10.0, 30.0, 50.0})
    for (double p : {0.5, 1.0, 2.0})
      for (int n : layer_counts)
        for (double t : temps) one_point(p, q, 1.0, n, t);
  // Deep drought and hot, for InfeasibleBracket / Determined / SolverRefused.
  for (double p : {7.0, 8.0, 10.0, 12.0, 20.0, 40.0})
    for (double v : {0.5, 4.0, 10.0})
      for (int n : layer_counts)
        for (double t : temps) one_point(p, 900.0, v, n, t);
  // Compensation point by heat, at the temperature the guide names.
  for (double t : {45.0, 50.0})
    for (double p : {0.5, 2.0})
      for (int n : layer_counts) one_point(p, 900.0, 1.0, n, t);

  printf("\nunrooted layers (the per-input refusal by name)\n");
  unrooted(2.0, 900.0, 1.0, 5, 2, 25.0);   // interior
  unrooted(4.0, 900.0, 4.0, 5, 2, 40.0);   // pinned, probably
  unrooted(20.0, 900.0, 1.0, 5, 2, 25.0);  // shut down

  printf("\n################ PER-KIND TABLE ################\n");
  long total = 0;
  for (const auto& kv : stats) total += kv.second.n_points;
  printf("%ld points over %zu kind(s) reached\n", total, stats.size());
  for (const auto& kv : stats) report_kind(kv.first, kv.second);

  printf("\nkinds NOT reached:");
  for (std::size_t i = 0;
       i < phylloptim::Leaf::operating_point_kind_count; ++i) {
    const Kind k = static_cast<Kind>(i);
    if (stats.find(k) == stats.end()) {
      printf(" %s", phylloptim::Leaf::operating_point_kind_name(k));
    }
  }
  printf("\n");
  if (!throw_messages.empty()) {
    printf("\nexceptions thrown by rows_at:\n");
    for (const auto& kv : throw_messages) {
      printf("  [%ld x] %s\n", kv.second, kv.first.c_str());
    }
  }
  return 0;
}
