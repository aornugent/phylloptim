// PROBE (not a test, not shipped): every field `rows_at` returns, over a grid and
// two request shapes, in a form a diff can referee.
//
// It is the instrument for a refactor of the row layer: capture before, capture
// after, diff. The two shapes are the two consumers -- a stand asks for profit and
// the uptake block, a calibration asks for all five fixed outputs.
//
//   make CXX=g++ probe_rows_ref && ./probe_rows_ref > before.tsv

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cstdio>
#include <exception>
#include <vector>

namespace grad = phylloptim::gradient;

namespace {

const double kRd25 = 1.44;
const double kTheta[] = {96.0,     2.680147, 3.898245, 5.870283, 2.680147,
                         3.898245, 5.870283, 1.5,      157.44,   0.30,
                         0.7,      0.99,     7.5,      kRd25,
                         1.0 * 0.000157 / 5.0, 1e3};

grad::Drivers drivers(double psi_soil, double ppfd, double vpd, double temp,
                      int layers) {
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
  d.atm_vpd = vpd;
  d.ca = 40.0;
  d.leaf_temp = temp;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

void dump(const char* tag, double psi_soil, double ppfd, double vpd, double temp,
          int layers, bool stand_shape) {
  const grad::Drivers d = drivers(psi_soil, ppfd, vpd, temp, layers);
  phylloptim::Leaf l;
  const grad::Settings s;
  grad::apply(l, kTheta, d, false, -1, s.fast_stem_curve);
  // The consumer solves; the row layer reads. Before the row layer stopped
  // solving, this call was inside it.
  try { l.find_root_collar_psi(); } catch (const std::exception&) {}

  std::vector<int> input;
  for (int p = 0; p < grad::n_pars; ++p) {
    if (p == grad::par_resistance) continue;  // single-potential path only
    input.push_back(p);
  }
  input.push_back(grad::par_PPFD);
  for (int j = 0; j < layers; ++j) input.push_back(grad::par_psi_soil_first + j);
  for (int j = 0; j < layers; ++j)
    input.push_back(grad::par_root_carbon_first(layers) + j);

  std::vector<int> output;
  if (stand_shape) {
    output.push_back(grad::out_profit);
    for (int j = 0; j < layers; ++j) output.push_back(grad::out_uptake_first + j);
  } else {
    for (int j = 0; j < grad::n_outputs; ++j) output.push_back(j);
    for (int j = 0; j < layers; ++j) output.push_back(grad::out_uptake_first + j);
  }

  const grad::RowRequest r{output.data(), output.size(), input.data(),
                           input.size()};
  grad::Rows rows;
  std::string err;
  try {
    rows = grad::rows_differenced(l, kTheta, d, r, s);
  } catch (const std::exception& e) {
    err = e.what();
  }
  const char* head = tag;
  std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\tkind\t%s\n", head, psi_soil, ppfd,
              vpd, temp, layers,
              err.empty() ? phylloptim::Leaf::operating_point_kind_name(rows.kind)
                          : "THREW");
  if (!err.empty()) {
    std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\terr\t%s\n", head, psi_soil,
                ppfd, vpd, temp, layers, err.c_str());
    return;
  }
  std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\tsolves\t%zu\n", head, psi_soil,
              ppfd, vpd, temp, layers, *l.collar_solves);
  std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\tpoint\t%.17g\t%.17g\n", head,
              psi_soil, ppfd, vpd, temp, layers, rows.point,
              rows.residual_slope);
  for (std::size_t i = 0; i < rows.dresidual.size(); ++i)
    std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\tdres\t%zu\t%.17g\n", head,
                psi_soil, ppfd, vpd, temp, layers, i, rows.dresidual[i]);
  for (std::size_t j = 0; j < rows.dy_dp.size(); ++j)
    std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\tdydp\t%zu\t%.17g\n", head,
                psi_soil, ppfd, vpd, temp, layers, j, rows.dy_dp[j]);
  for (std::size_t k = 0; k < rows.held.size(); ++k)
    std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\theld\t%zu\t%.17g\n", head,
                psi_soil, ppfd, vpd, temp, layers, k, rows.held[k]);
  if (!rows.message.empty())
    std::printf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%d\tmsg\t%s\n", head, psi_soil,
                ppfd, vpd, temp, layers, rows.message.c_str());
}

}  // namespace

int main() {
  for (double temp : {25.0, 40.0})
    for (double vpd : {0.5, 2.0, 4.0})
      for (int L = 1; L <= 5; ++L)
        for (double ppfd : {0.0, 20.0, 100.0, 900.0, 1500.0})
          for (double ps = 0.5; ps <= 6.01; ps += 0.5) {
            dump("stand", ps, ppfd, vpd, temp, L, true);
            dump("calib", ps, ppfd, vpd, temp, L, false);
          }
  return 0;
}
