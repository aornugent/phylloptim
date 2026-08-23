// -*-c++-*-
#ifndef PHYLLOPTIM_INPUTS_HPP_
#define PHYLLOPTIM_INPUTS_HPP_

#include <phylloptim/util.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

// What the leaf's inputs and outputs are called and how they are numbered. Below
// both the model and the row layer, because both need it: the row layer indexes
// by these, and the model has to place a value it is handed by index. Held here
// so neither has to spell the order the other uses.

namespace phylloptim {
namespace gradient {

// --- the parameter enumeration, which R indexes into --------------------------
//
// The fourteen traits, in the order the model places them, then the two
// quantities a calibration fits that are not traits: the conductance driver and
// the single-potential path's series resistance.
//
// ⚠️ R INDEXES THESE POSITIONS, so a reordering silently differentiates the wrong
// parameter. `test-gradient-batch.R` reads the names back out of C++ and compares
// them with R's, so the two cannot drift apart without a failure.
// Which physical route an input reaches the leaf by, which is what decides the
// shape of its row. Carried on the entry rather than tested for by four separate
// lists of names, so the four are a partition by construction and an input
// belonging to none says so.
enum class InputRole {
  Carbon,     // reaches profit through assimilation or the hydraulic cost
  Transport,  // moves the stem potential at a frozen flux
  Supply,      // moves the supply, so it reaches the leaf through total uptake
  Slack,      // a limit the point may or may not be sitting on
  None        // no row here: the other supply path's input
};

// THE list, in order. The first `n_traits` are the traits `set_traits` places and
// `leaf_traits()` names; the two non-traits follow and take a relative step.
struct par_entry {
  std::string_view name;
  InputRole role;
};

inline constexpr std::array<par_entry, 16> par_table{{
    {"vcmax_25", InputRole::Carbon},
    {"stem_c", InputRole::Transport},
    {"stem_b", InputRole::Transport},
    {"psi_crit", InputRole::Slack},
    {"root_c", InputRole::Supply},
    {"root_b", InputRole::Supply},
    {"root_psi_crit", InputRole::Slack},
    {"beta2", InputRole::Carbon},
    {"jmax_25", InputRole::Carbon},
    {"a", InputRole::Carbon},
    {"curv_fact_elec_trans", InputRole::Carbon},
    {"curv_fact_colim", InputRole::Carbon},
    {"cost_scale_TF24", InputRole::Carbon},
    {"R_d_25", InputRole::Carbon},
    {"leaf_specific_conductance_max", InputRole::Transport},
    // The single-potential path's series resistance. No closed form here says how
    // it moves the bound, so the read declines it and a difference answers
    // instead -- which is why it is the one input no role claims.
    {"resistance", InputRole::None}}};

inline constexpr int n_pars = static_cast<int>(par_table.size());

// A parameter's index, found in the list that names it. Every constant below is
// this, so an index and the name it stands for cannot be written down in two
// places -- and a name the list does not hold is not a constant expression, so a
// misspelling is a compile error at the constant rather than a -1 nobody checks.
inline constexpr int par_of(std::string_view name) {
  for (std::size_t i = 0; i < par_table.size(); ++i) {
    if (par_table[i].name == name) {
      return static_cast<int>(i);
    }
  }
  util::stop("par_of: `" + std::string(name) + "` is not a parameter");
  return -1;
}

// Every index by name, so nothing below indexes `theta` with a bare integer.
inline constexpr int par_vcmax_25 = par_of("vcmax_25");
inline constexpr int par_stem_c = par_of("stem_c");
inline constexpr int par_stem_b = par_of("stem_b");
inline constexpr int par_psi_crit = par_of("psi_crit");
inline constexpr int par_root_c = par_of("root_c");
inline constexpr int par_root_b = par_of("root_b");
inline constexpr int par_root_psi_crit = par_of("root_psi_crit");
inline constexpr int par_beta2 = par_of("beta2");
inline constexpr int par_jmax_25 = par_of("jmax_25");
inline constexpr int par_a = par_of("a");
inline constexpr int par_curv_fact_elec_trans = par_of("curv_fact_elec_trans");
inline constexpr int par_curv_fact_colim = par_of("curv_fact_colim");
inline constexpr int par_cost_scale_TF24 = par_of("cost_scale_TF24");
inline constexpr int par_R_d_25 = par_of("R_d_25");
inline constexpr int par_kmax = par_of("leaf_specific_conductance_max");
inline constexpr int par_resistance = par_of("resistance");

// The traits are everything before the first non-trait, read off the list
// rather than counted: a trait added ahead of the conductance driver would
// leave a hand-written count naming one of the two that are not traits.
inline constexpr int n_traits = par_kmax;

inline const std::vector<std::string>& par_names() {
  static const std::vector<std::string> names = [] {
    std::vector<std::string> ret;
    ret.reserve(par_table.size());
    for (const par_entry& e : par_table) {
      ret.emplace_back(e.name);
    }
    return ret;
  }();
  return names;
}

// --- the environment rows: the two things plants share ------------------------
//
// Everything above differentiates the leaf with respect to what it IS. These rows
// differentiate it with respect to what it EXPERIENCES. Plants in this model are
// coupled through exactly two quantities -- the light they cast on each other and
// the soil water they draw from -- so a gradient that stops at the traits has
// nothing to say about competition: one plant's uptake is the next plant's
// `psi_soil`, one plant's canopy is the next plant's `PPFD`, and neither had a
// row.
//
// THE ENVIRONMENT ROWS ARE NOT IN `theta`, and that is the design decision. Their
// values already have exactly one home -- the observation's `Drivers` -- and
// copying them into `theta` would give one solve two soil states free to
// disagree. So `theta` stays sixteen wide, R's sixteen-column contract and
// `gradient_par_names()` are untouched, and an environment row is addressed by an
// index PAST the end of `theta` whose value is read from the drivers instead
// (`par_value` below). Nothing else about the two routes changes: `psi_soil` and
// `PPFD` are drivers exactly as `leaf_specific_conductance_max` is, so the
// existing perturbation loops differentiate them with no new algebra.
//
// ARITY IS THE NEW THING. Everything in `theta` is a scalar; `psi_soil` is a
// vector of length L, so the row count depends on the layer configuration and is
// not a compile-time constant. It is handled by keeping the FIXED part fixed:
// `PPFD` sits at `par_PPFD` whatever L is, and the L soil rows follow it. An
// index therefore means the same parameter across observations with different
// layer counts, which it would not if the variable-length block came first. The
// soil rows are named `psi_soil_1..psi_soil_L`, one-based to match the layer
// numbering the caller already has, and generated per call rather than cached
// because L is a property of the observation and not of this header.
inline constexpr int par_PPFD = n_pars;
inline constexpr int par_psi_soil_first = n_pars + 1;

// One root-carbon row per layer, after the soil block: a consumer's input vector
// carries two entries per layer, so the arity is n_pars + 1 + 2L. With two
// variable-length blocks an index past `PPFD` names one input only at a fixed L.
inline constexpr int par_root_carbon_first(int n_layers) {
  return par_psi_soil_first + n_layers;
}

// Total rows for a given layer count, and the full row names in index order.
inline constexpr int n_pars_total(int n_layers) {
  return par_root_carbon_first(n_layers) + n_layers;
}

// Which block an input index falls in, and where in that block. The layout above
// is arithmetic, and every reader of it used to redo the arithmetic: seven
// functions each subtracted par_psi_soil_first and compared against n_layers, so
// moving a block meant finding all seven.
struct par_ref {
  enum class Block { Parameter, Radiation, SoilPotential, RootCarbon };
  Block block;
  int index;
};

// What part one of the sixteen parameters plays. The blocks past them are the
// supply's, which input_role() below adds -- so a caller asking about a role only
// a parameter can have needs no layer count.
inline constexpr InputRole parameter_role(int par) {
  return par >= 0 && par < n_pars ? par_table[std::size_t(par)].role
                                  : InputRole::None;
}

// ⚠️ AN INPUT NO ROLE CLAIMS HAS NO ROW, AND ONE IS EXPECTED. `rows_at`
// defaults every row to NA and the recording refuses a non-finite derivative by name,
// so an unclaimed input fails safe -- but it fails safe by accident unless the
// set of them is known. `resistance` is the one, and this is what says so.
static_assert(
    [] {
      for (int p = 0; p < n_pars; ++p) {
        if (parameter_role(p) == InputRole::None && p != par_resistance) {
          return false;
        }
      }
      return true;
    }(),
    "a parameter has no role, so it has no row, and only `resistance` is "
    "meant to be in that position");

inline constexpr par_ref decode(int par, int n_layers) {
  if (par < n_pars) {
    return {par_ref::Block::Parameter, par};
  }
  if (par == par_PPFD) {
    return {par_ref::Block::Radiation, 0};
  }
  const int layer = par - par_psi_soil_first;
  if (layer < n_layers) {
    return {par_ref::Block::SoilPotential, layer};
  }
  return {par_ref::Block::RootCarbon, layer - n_layers};
}

// What part an input plays in the solve, for ANY input index. The blocks past the
// parameters belong to the supply by construction: a soil potential and a layer's
// carbon both reach the leaf through total uptake and nothing else.
inline constexpr InputRole input_role(int par, int n_layers) {
  if (par < 0 || par >= n_pars_total(n_layers)) {
    return InputRole::None;
  }
  const par_ref r = decode(par, n_layers);
  switch (r.block) {
  case par_ref::Block::Parameter: return parameter_role(r.index);
  case par_ref::Block::Radiation: return InputRole::Carbon;
  default:                       break;
  }
  return InputRole::Supply;
}


inline std::vector<std::string> par_names(int n_layers) {
  std::vector<std::string> out = par_names();
  out.reserve(std::size_t(n_pars_total(n_layers)));
  out.emplace_back("PPFD");
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("psi_soil_" + std::to_string(i + 1));
  }
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("root_carbon_" + std::to_string(i + 1));
  }
  // The arity and the list are two readings of one layout, and this is where
  // they meet. A block added to one and not the other stops here.
  if (out.size() != std::size_t(n_pars_total(n_layers))) {
    util::stop("par_names: " + std::to_string(out.size()) + " names against " +
               std::to_string(n_pars_total(n_layers)) + " parameters");
  }
  return out;
}

// One name, for a diagnostic message. Read out of the whole list rather than
// re-deriving the block layout, which is the one place that layout is written:
// two spellings of it agree until a block moves. A diagnostic can afford the
// allocation.
inline std::string par_name(int par, int n_layers) {
  const std::vector<std::string> names = par_names(n_layers);
  if (par < 0 || par >= static_cast<int>(names.size())) {
    return "parameter " + std::to_string(par);
  }
  return names[std::size_t(par)];
}

// --- the five differentiated outputs -----------------------------------------
//
// A, gc, psi_stem, collar and profit, in that order. THIS LIST IS THE ONE
// DEFINITION OF IT: R does not keep a second copy, it reads this one back through
// `gradient_output_names()` in src/gradient.cpp, so the two routes cannot
// disagree about what they are reporting. `profit` is appended rather than
// inserted, because R hands `pars` over as positions and a caller reads gradient
// columns by position after that.
//
// Two of the five are not the generic composite in `gradient_ift` below, and for
// two DIFFERENT reasons, which is why both are named here:
//
//   * `collar` is psi* itself, which is what makes dcollar/dtheta equal
//     dpsi*/dtheta and lets the two routes compute the same quantity by
//     different means.
//   * `profit` is the OBJECTIVE rather than an output read at the argmax, which
//     is what brings the envelope theorem into play.
inline constexpr std::array<std::string_view, 5> output_table{
    "A", "gc", "psi_stem", "collar", "profit"};

inline constexpr int n_outputs = static_cast<int>(output_table.size());

// An output's index, found in the list that names it, for the reason `par_of`
// gives: an index and the name it stands for are one entry, and a name the list
// does not hold does not compile.
inline constexpr int output_of(std::string_view name) {
  for (std::size_t i = 0; i < output_table.size(); ++i) {
    if (output_table[i] == name) {
      return static_cast<int>(i);
    }
  }
  util::stop("output_of: `" + std::string(name) + "` is not an output");
  return -1;
}

inline constexpr int out_assim = output_of("A");
inline constexpr int out_stom_cond = output_of("gc");
inline constexpr int out_psi_stem = output_of("psi_stem");
inline constexpr int out_collar = output_of("collar");
inline constexpr int out_profit = output_of("profit");

inline const std::vector<std::string>& output_names() {
  static const std::vector<std::string> names(output_table.begin(),
                                              output_table.end());
  return names;
}

// --- the per-layer uptake outputs ---------------------------------------------
//
// The five above are a calibration's, and R reads gradient columns by position,
// so they stay fixed: `output_names()` and `gradient_output_names()` report those
// five whatever the layer configuration is, and `at`, `batch` and `transpose_at`
// report those five. A stand adjoint reads a different list -- of the five it
// wants profit, and it wants the water each soil layer gave up -- so the uptake
// entries are APPENDED, one per layer, and `rows_at` is the one route that can
// request them.
//
// ⚠️ AN OUTPUT THAT IS NOT REPORTED HAS NO ROUTE TO ANY INPUT. At a pin the
// condition's gradient is zero for every input and each reported output carries a
// total instead, so an unreported one comes back exactly zero at every pinned and
// shut point -- and zero and absent are the same number on a consumer's tape.
// That is why uptake is an enumeration entry rather than a second call.
//
// Uptake is ORDINARY: it is set as a side effect at the operating point, so it
// consumes the collar rather than being it, and it composes the way A and gc do.
// Fixed first, variable second, for the reason the soil rows follow `PPFD` on the
// input side: an index means the same output across observations with different
// layer counts.
inline constexpr int out_uptake_first = n_outputs;

inline constexpr int n_outputs_total(int n_layers) {
  return out_uptake_first + n_layers;
}

inline std::vector<std::string> output_names(int n_layers) {
  std::vector<std::string> out = output_names();
  out.reserve(std::size_t(n_outputs_total(n_layers)));
  for (int i = 0; i < n_layers; ++i) {
    out.push_back("uptake_" + std::to_string(i + 1));
  }
  // The arity and the list are two readings of one layout, and this is where
  // they meet.
  if (out.size() != std::size_t(n_outputs_total(n_layers))) {
    util::stop("output_names: " + std::to_string(out.size()) +
               " names against " + std::to_string(n_outputs_total(n_layers)) +
               " outputs");
  }
  return out;
}

// One name, for a diagnostic message. Read out of the whole list rather than
// re-deriving the block layout, for the reason `par_name` gives.
inline std::string output_name(int out, int n_layers) {
  const std::vector<std::string> names = output_names(n_layers);
  if (out < 0 || out >= static_cast<int>(names.size())) {
    return "output " + std::to_string(out);
  }
  return names[std::size_t(out)];
}

// How the operating point reaches an output. A PROPERTY OF THE OUTPUT, which is
// why `output_role` below is a function of the index and not a field of the request:
// which output IS the objective is fixed by the enumeration, and a caller free to
// say otherwise is a caller free to disagree with it.
enum class OutputRole { Objective, Point, Ordinary };

// The collar IS the operating point; profit is what it maximises; everything else
// reads the point without being it. Read here and nowhere else -- `out_collar` and
// `out_profit` are ordinary indices at every other site.
inline OutputRole output_role(int output) {
  return output == out_collar   ? OutputRole::Point
       : output == out_profit   ? OutputRole::Objective
                                : OutputRole::Ordinary;
}

// dy/dp, from what the output is and what defines the point. That pairing is the
// whole sensitivity theory and this is the one place either axis is read.
//
// ⚠️ THE OBJECTIVE'S ZERO IS NOT A PROPERTY OF BEING THE OBJECTIVE. It is the
// interior stationarity condition, so it holds where that condition does; at a
// bound the same output's dy/dp is the constraint's shadow price. An interface
// that let a consumer infer the zero from the output's identity is correct at an
// interior optimum and silently wrong at exactly the states a pin exists for.
//
// The point's own two are genuine identities: one, because it IS the
// point, and a held partial cannot move what it holds. Neither reads the kind.
inline double dy_dp_for(OutputRole role, bool pinned, double marginal_at_bound,
                            double measured) {
  switch (role) {
  case OutputRole::Point:     return 1.0;
  case OutputRole::Objective: return pinned ? marginal_at_bound : 0.0;
  default:              return measured;
  }
}

// One evaluation's outputs: the five, then one uptake per layer where the caller
// asked for them. A fixed-five caller leaves `uptake` empty, which allocates
// nothing.
struct OutputValues {
  double fixed[n_outputs];
  std::vector<double> uptake;

  OutputValues() = default;
  explicit OutputValues(int n_uptake) : uptake(std::size_t(n_uptake)) {}

  int n_uptake() const { return int(uptake.size()); }
  int size() const { return n_outputs + n_uptake(); }
  double& operator[](int j) {
    return j < n_outputs ? fixed[j] : uptake[std::size_t(j - out_uptake_first)];
  }
  double operator[](int j) const {
    return j < n_outputs ? fixed[j] : uptake[std::size_t(j - out_uptake_first)];
  }
};
}
}

#endif
