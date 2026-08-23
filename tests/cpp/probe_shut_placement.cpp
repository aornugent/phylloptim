// PROBE (not a test, not shipped): what a shut collar's own row would have to be.
//
// Three exits write a shut operating point and they placement the collar at three
// different potentials. The question is whether each placement is a NEW object or one
// the bound machinery already computes, because that decides whether the missing
// row is a derivation or a recorded enumerator.
//
// The answer is a count: over a wide sweep every hydraulic shutdown places the
// collar at one of two REGISTERED INPUTS, the searched dry bound is the placement
// nowhere, and profit is bit-for-bit the two-term expression at the STEM
// potential -- so nothing a shut leaf reports is a function of its collar.
//
//   make CXX=g++ probe_shut_placement && stdbuf -oL ./probe_shut_placement

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace grad = phylloptim::gradient;
using Kind = grad::OperatingPointKind;
using Which = phylloptim::Leaf::WhichBound;

namespace {

const double kRd25 = 1.44;
double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                   3.898245, 5.870283, 1.5,      157.44,   0.30,
                   0.7,      0.99,     7.5,      kRd25,
                   1.0 * 0.000157 / 5.0, 1e3};

grad::Drivers drivers(double psi_soil, double ppfd, int layers) {
  grad::Drivers d;
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers / 0.05;
  }
  d.root_network = fixture::root_network(root, depth);
  d.PPFD = ppfd;
  d.psi_soil = ps;
  d.soil_depth = depth;
  d.atm_vpd = 2.0;
  d.ca = 40.0;
  d.leaf_temp = 25.0;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

}  // namespace

namespace {
struct Tally {
  long shut = 0, shade = 0, infeasible = 0;
  long placed_at_psi_crit_bound = 0, placed_at_root_psi_crit = 0, placed_at_neither = 0;
  long placed_at_dry_searched = 0, placed_at_wet_bound = 0;
  long profit_two_term = 0, all_flux_zero = 0;
  double worst_neither = 0.0;
  long wet_refused = 0;
  std::string wet_why;
};

void tally_one(Tally& t, double psi_soil, double ppfd, int layers,
               double psi_crit_theta, double vpd, double temp) {
  double th[16];
  for (int i = 0; i < 16; ++i) th[i] = kTheta[i];
  th[grad::par_psi_crit] = psi_crit_theta;
  grad::Drivers d = drivers(psi_soil, ppfd, layers);
  d.atm_vpd = vpd;
  d.leaf_temp = temp;
  phylloptim::Leaf l;
  const grad::Settings s;
  grad::apply(l, th, d, false, -1, s.fast_stem_curve);
  try { l.find_root_collar_psi(); } catch (const std::exception&) { return; }
  const Kind k = l.operating_point_kind();
  if (k == Kind::HydraulicShutdown) t.shut++;
  else if (k == Kind::ShadeDeath) t.shade++;
  else if (k == Kind::InfeasibleBracket) t.infeasible++;
  else return;

  const double placement = l.opt_root_psi_;
  bool named = false;
  if (placement == l.psi_crit) { t.placed_at_psi_crit_bound++; named = true; }
  if (placement == l.supply_psi_crit()) { t.placed_at_root_psi_crit++; named = true; }
  if (!named) {
    // Is it one of the two bounds the machinery already computes?
    try {
      if (placement == l.bound_row(Which::DryRootCrit).bound) {
        t.placed_at_dry_searched++; named = true;
      }
    } catch (const std::exception&) {}
    bool wet_refused = false;
    try {
      if (placement == l.bound_row(Which::Wet).bound) {
        t.placed_at_wet_bound++; named = true;
      }
    } catch (const std::exception&) { wet_refused = true; }
    if (!named) t.wet_refused += wet_refused ? 1 : 0;
  }
  if (!named) {
    t.placed_at_neither++;
    t.worst_neither = placement;
    // What IS it, then? root_crit is find_root_psi(wettest, psi_soil, 1) -- the
    // same call bound_row(DryRootCrit) makes, so ask that directly.
    double rc = std::nan("");
    try { rc = l.find_root_psi(l.supply_begin_solve(), l.supply_psi_soil(), 1); }
    catch (const std::exception&) {}
    std::printf("  unnamed placement: %s psi_soil=%.3f ppfd=%.0f L=%d psi_crit=%.6f "
                "root_psi_crit=%.6f placement=%.17g root_crit=%.17g match=%d\n",
                phylloptim::Leaf::operating_point_kind_name(k), psi_soil, ppfd,
                layers, l.psi_crit, l.supply_psi_crit(), placement, rc,
                int(placement == rc));
  }
  const double predicted = -l.R_d_ - l.hydraulic_cost_TF(l.opt_psi_stem_);
  if (l.profit_ == predicted) t.profit_two_term++;
  bool zero = (l.E_up_ == 0.0);
  if (k == Kind::ShadeDeath) zero = true;  // per-layer draws are non-zero there
  else for (double u : l.soil_consumption_) if (u != 0.0) zero = false;
  if (zero) t.all_flux_zero++;
}
}  // namespace

int main() {
  Tally t;
  const double pcs[] = {5.870283, 5.0, 4.0, 6.5, 3.0};
  const double vpds[] = {1.0, 2.0, 4.0};
  const double temps[] = {25.0, 40.0};
  for (double pc : pcs)
    for (double vpd : vpds)
      for (double temp : temps)
        for (int L = 1; L <= 5; ++L)
          for (double ppfd : {0.0, 10.0, 30.0, 100.0, 900.0, 1500.0})
            for (double ps = 0.25; ps <= 8.01; ps += 0.25)
              tally_one(t, ps, ppfd, L, pc, vpd, temp);

  std::printf("=== shut and zero-flux points over a wide sweep ===\n");
  std::printf("hydraulic-shutdown %ld   shade-death %ld   infeasible-bracket %ld\n",
              t.shut, t.shade, t.infeasible);
  const long n = t.shut + t.shade + t.infeasible;
  std::printf("\nwhere the placement already is, over %ld points:\n", n);
  std::printf("  == psi_crit (a registered input)          %ld\n", t.placed_at_psi_crit_bound);
  std::printf("  == root_psi_crit (a registered input)     %ld\n", t.placed_at_root_psi_crit);
  std::printf("  == bound_row(DryRootCrit).bound           %ld\n", t.placed_at_dry_searched);
  std::printf("  == bound_row(Wet).bound                   %ld\n", t.placed_at_wet_bound);
  std::printf("  none of the four  ---------------------->  %ld", t.placed_at_neither);
  if (t.placed_at_neither)
    std::printf("   (of which the wet bound REFUSED at %ld; one placement: %.17g)",
                t.wet_refused, t.worst_neither);
  if (!t.wet_why.empty())
    std::printf("\n  the wet bound's refusal says: %s\n", t.wet_why.c_str());
  std::printf("\n\nprofit == -R_d - C(psi_stem) bit-for-bit: %ld of %ld\n",
              t.profit_two_term, n);
  std::printf("every flux written zero (shutdown only):  %ld of %ld\n",
              t.all_flux_zero, n);
  return 0;
}
