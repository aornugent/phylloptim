# leaf

A C++ leaf gas-exchange model, callable from R, in which **stomatal behaviour
emerges from hydraulics** instead of an empirical conductance function.

Farquhar-von Caemmerer-Berry photosynthesis is coupled to an explicit
soil → root → stem → leaf water transport path — Weibull vulnerability curves
for both xylem and roots, multi-layer soil, per-layer root resistance,
gravitational head — and the operating point is chosen by Sperry-style gain-risk
**profit maximisation** over the root-collar water potential. Forward-mode
automatic differentiation (XAD) supplies exact derivatives of profit with respect
to that potential, for models that need to track acclimation.

A full solve costs about **3 µs**, which is what makes it usable inside a
demographic model that calls it millions of times. That is the C++ figure; from R
a solved row costs ~20 µs, nearly all of it the R boundary rather than the model —
see [Performance from R](#performance-from-r) before optimising anything.

This code was developed as the TF24 strategy inside
[traitecoevo/plant](https://github.com/traitecoevo/plant) and is extracted here
so it can be tested, profiled, extended and embedded on its own.

## Why a separate package
**A home for several stomatal models, not just ours.** The package solves **seven**
optimality models through one optimiser, at identical drivers, behind one gradient
entry point: our hydraulic gain-risk formulation (`TF24`), Sperry-style profit
maximisation (`ProfitMax`), Cowan-Farquhar (`CF77`), Joshi & Stocker's quadratic
cost (`JS22`), Wolf-Anderegg-Pacala carbon maximisation (`CMax`), and the two
product objectives `SOX` and `JW26`. The Medlyn et al. (2011) empirical model is
here too, bypassing the hydraulic solve. `vignette("the-models")` sets them side
by side and derives each one.

Making that comparison apples-to-apples needed one level more generality than
swapping cost functions. Every one of these models maximises

    h(A(psi)) - C(psi)

so a model is fully specified by a **cost curve** `C` and a **benefit link** `h`,
and one derivative serves all seven:

    d/dpsi [h(A) - C] = h'(A) * dA/dpsi - dC/dpsi

The marginal cost of water λ — the quantity that distinguishes these models in the
literature — is a *consequence* of that pair rather than a primitive:
λ = (dC/dψ)/(dE/dψ). That is what lets the product objectives (`A·g`) sit
alongside the difference objectives at all: a log link turns one into the other
without moving the argmax. This is the central result of a companion manuscript,
and it is something none of the existing R packages can support, because each
commits to a single hydraulically explicit scheme or to none.

**One solver, so the comparison is not just fair but cheap.** Every model reaches its
operating point through the same call: evaluate both interval endpoints, optionally
scan for the basin, then root-find `dJ/dψ = 0` inside the winning cell. Nothing here
is per-model except a row in two dispatch tables. That is what keeps a seven-model
package as fast as a one-model one — the stem entry points run at **2.9 µs** and
ProfitMax at **5.6 µs**, against 9.3 and 58.2 before the solvers were merged, because
a basin scan is now used only where multi-modality is measured rather than everywhere
by default. The residual `|dJ/dψ|` at the returned optimum has a median of **1.2e-15**.

**Either supply path, so the comparison is fair.** The gas-exchange core is
soil-agnostic: the multi-layer soil and root system enter the solve only as a
supply function `E_up = f(P_collar)`. So the multi-layer root system can be
swapped for a single soil water potential (`leaf_supply_singlelayer()`), which both
lowers the barrier for a bare-leaf user — no root-mass profile to construct — and
is what makes the comparison fair, because the alternative formulations worth
comparing against are all written for one ψ_soil.

**Fast and differentiable, so it is built for calibration.** A full hydraulic
solve costs ~3 µs, and derivatives are analytic rather than finite differences:
forward-mode AD (XAD) for the collar potential, and `leaf_gradient()` for the
traits, by differentiating the optimality condition. That combination is what
calibration wants — gradient-based optimisers and Hamiltonian samplers need many
evaluations *and* clean gradients, and finite-differencing a nested root-find is
exactly the case where numerical gradients are noisiest.

One honest caveat: there is still no calibration vignette in the package, though
the fit that drove the gradient work exists outside it in
a companion calibration study.

## Status

**v0.6.0 — the model is mature and in production use inside plant; the
*packaging* is what is new.**

- **Cross-checked against plant's compiled build**, and the swap was bit-identical
  at the point it was made: plant's full suite 0 fail / 0 error on both builds,
  and the SCM regression identical across 78/78 nodes. The 1-ULP disagreement
  that held this up turned out to be R's decimal parser rather than either model.
- **The shutdown defect is fixed.** On the hydraulic-shutdown path, transpiration,
  assimilation and uptake were left holding the previous solve's values (plant
  #578). That, and three further stale-state exits ported from plant #585, are all
  fixed here.
- **Results now differ from plant's own leaf, deliberately** — see
  [NEWS.md](NEWS.md). The most consequential single change is deriving the
  ppm→Pa conversion from the actual atmospheric pressure instead of a hard-coded
  101.3 kPa, which moves TF24 offspring production by 2.4% in plant, because
  plant's driver default is 100.5.

## Two ways in, and the model does not need R

The model is a set of self-contained C++ headers under `inst/include`. They use
no R and no Rcpp, depend only on Boost and the header-only parts of odelia, and
compile and run with no R installed — so the same model is available to a C++
program, to a Python extension, and to R, with none of those paying for the
others. The R layer (`src/`, `R/`) sits on top of those headers and is never
included by them; the dependency runs one way only. See
`.github/workflows/cpp-tests.yml` for what enforces it —
it builds the whole C++ suite on a runner with no R on it.

## Use from C++

There is nothing to link against — one include is the whole library.

```cpp
#include <phylloptim.hpp>

phylloptim::Leaf l;                   // default Eucalyptus saligna traits
l.setup_transpiration(100);     // build the xylem vulnerability splines
l.setup_root_vulnerability(100);

std::vector<double> psi_soil{2.0};             // positive suction, MPa
std::vector<double> soil_depth{1.0};           // m

// The per-layer root hydraulic RESISTANCES, per unit leaf area: the leaf is purely
// intensive, and it takes the resistances rather than the root carbon they came
// from. If you have carbon, this is the root-architecture model that maps one to
// the other -- a helper you call, not something the solve does for you.
const phylloptim::RootNetwork roots = phylloptim::root_network_from_carbon(
    /*kg C per m2 LEAF*/ {20.0}, phylloptim::layer_thickness(soil_depth),
    /*beta_R_H*/ 3.4e2, /*beta_R_V*/ 9.4e3);

l.set_physiology(roots, /*PPFD*/ 900,
                 psi_soil, soil_depth,
                 /*leaf_specific_conductance_max*/ 3.14e-5,
                 /*atm_vpd*/ 2.0, /*ca*/ 40.0,
                 /*leaf_temp*/ 25.0, /*atm_o2_kpa*/ 21.0, /*atm_kpa*/ 101.3);

l.find_root_collar_psi();       // solve

l.opt_psi_stem_;                // leaf water potential at the optimum, MPa
l.opt_root_psi_;                // root-collar potential, MPa (positive magnitude)
l.assim_colimited_;             // A, umol m-2 s-1
l.transpiration_;               // E, kg H2O m-2 s-1
l.stom_cond_CO2_;               // gc, mol CO2 m-2 s-1
l.profit_;                      // A - hydraulic cost
l.soil_consumption_;            // per-layer water uptake
```

Compile with C++20 and three include paths — this package, `odelia`, and Boost:

```sh
c++ -std=c++20 -O2 \
  -I /path/to/leaf/inst/include \
  -isystem /path/to/odelia/inst/include \
  -isystem /path/to/boost \
  my_program.cpp -o my_program
```

Or use the CMake package, which handles those three paths for you. `odelia` is
distributed as an R package but the headers used here are plain C++, so a git
checkout of it is enough — nothing needs installing or building:

```sh
git clone https://github.com/traitecoevo/odelia
cmake -B build -DPHYLLOPTIM_ODELIA_INCLUDE_DIR=$PWD/odelia/inst/include
cmake --build build
ctest --test-dir build          # runs the C++ suite, including the golden file
cmake --install build --prefix /usr/local
```

```cmake
find_package(phylloptim REQUIRED)
target_link_libraries(my_program PRIVATE phylloptim::phylloptim)
```

`phylloptim::phylloptim` is an INTERFACE target — headers, an include path and `cxx_std_20`,
with nothing to link. `add_subdirectory(leaf)` works the same way if you would
rather vendor it.

⚠️ **Build with optimisation on if you intend to compare against the golden
file.** Measured on macOS/arm64: `-O1`, `-O2` and `-O3` all reproduce
`tests/cpp/golden/operating_points.tsv` bit-for-bit and `-O0` does not, missing
by about 13 ULP because it declines to contract `a*b + c` into an FMA. A debug
build that fails `test_golden` by ~1e-15 has found nothing. The CMake build
therefore defaults to `Release` rather than to CMake's flagless default.

## Use from Python

The same headers, through [pybind11](https://pybind11.readthedocs.io/) — no R
anywhere in the picture. A minimal extension module:

```cpp
// pyleaf.cpp
#include <phylloptim.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(pyleaf, m) {
  py::class_<phylloptim::Leaf>(m, "Leaf")
      .def(py::init<>())
      .def("set_physiology", &phylloptim::Leaf::set_physiology)
      .def("find_root_collar_psi", &phylloptim::Leaf::find_root_collar_psi)
      .def_readonly("profit", &phylloptim::Leaf::profit_)
      .def_readonly("opt_psi_stem", &phylloptim::Leaf::opt_psi_stem_)
      .def_property_readonly("g1_eff", &phylloptim::Leaf::g1_eff);
}
```

```cmake
find_package(phylloptim REQUIRED)
find_package(pybind11 REQUIRED)
pybind11_add_module(pyleaf pyleaf.cpp)
target_link_libraries(pyleaf PRIVATE phylloptim::phylloptim)
```

```python
>>> import pyleaf
>>> l = pyleaf.Leaf()
>>> roots = pyleaf.root_network_from_carbon([20.0], 1.0, 340.0, 9400.0)
>>> l.set_physiology(roots, 900, [2.0], [1.0], 3.14e-5, 2.0, 40.0, 25.0, 21.0, 101.3)
>>> l.find_root_collar_psi()
>>> l.profit
2.5158434915102319
```

`set_physiology` takes a `RootNetwork`, so the binding above needs it exposed too:

```cpp
  py::class_<phylloptim::RootNetwork>(m, "RootNetwork")
      .def(py::init<>())
      .def_readwrite("r_R_H_min", &phylloptim::RootNetwork::r_R_H_min)
      .def_readwrite("r_R_V_sum", &phylloptim::RootNetwork::r_R_V_sum);
  m.def("root_network_from_carbon",
        py::overload_cast<const std::vector<double>&, double, double, double>(
            &phylloptim::root_network_from_carbon));
```

That value is the golden file's `profit` at this operating point, to the last
bit — which is the point of the example. `util::stop` throws
`std::runtime_error`, so pybind11 turns the model's input validation into a
`RuntimeError` with no extra work.

`std::vector<double>` needs `pybind11/stl.h`, as above; swap it for
`pybind11/numpy.h` and an `Eigen`-style binding if you want the soil profile to
arrive as an array without a copy.

## Use from R

Drivers in, operating point out. `leaf_solve()` is vectorised, so a response
curve is one call:

```r
library(phylloptim)

leaf_solve(psi_soil = 2.0, PPFD = 900)
#>   psi_soil layers PPFD atm_vpd ca leaf_temp atm_kpa psi_stem  collar    ci
#> 1        2      1  900       2 40        25   101.3 3.595247 2.92039 10.49
#>          A         E         gc   profit ... lambda    g1_eff
#> 1 5.599511 1.142e-05 0.01921993 2.515843 ... 159884.6 0.5025448

# a drought response
leaf_solve(psi_soil = seq(0.5, 5, length.out = 20), PPFD = 900)
```

`gc` is not from a fitted conductance model — it is what falls out of maximising
profit over the hydraulic path. `lambda` is the marginal cost of water, dA/dE, at
the operating point, and `g1_eff` re-expresses the solved conductance as a Medlyn
`g1`, which is a convenient common scale for comparison.

Traits and numerical settings are separate, so a calibration loop varying traits
never has to know which of the C++ constructor's fifteen arguments are
tolerances:

```r
leaf_solve(psi_soil = 3.0, PPFD = 900,
           traits  = leaf_traits(vcmax_25 = 120, stem_P50 = 2.5),
           control = leaf_control(GSS_tol_abs = 1e-5))
```

For the stateful interface — which is what plant uses, and what you want if you
care about intermediate state:

```r
l <- leaf_model()                          # or leaf_model(traits, control)
set_drivers(l, psi_soil = 2.0, PPFD = 900)
l$optimise()

operating_point(l)   # the same one-row data.frame
l$profit_            # or reach into the object directly
l$lambda             # marginal cost of water, dA/dE
```

**Which model is configuration, not a call argument.** `$set_model()` seats a cost
curve, a route and a method; `$optimise()` takes nothing, because the curves'
constants are already on the object. The defaults are `"TF24"`, `"collar"` and
`"exact"`, the production path, so an existing caller that never sets a model is
unaffected — `$find_root_collar_psi()` is that solve under its own name, which is
what plant's sources spell.

```r
l <- leaf_model(supply = leaf_supply_singlelayer())
set_drivers(l, psi_soil = 1.5, PPFD = 900)
l$set_model("SOX", "stem")     # any of cost_curve_names(); "collar" or "stem"
l$optimise()
l$model_curve(); l$model_route(); l$model_method()
```

The third axis chooses *how* the operating point is reached rather than *what* is
optimised. `"exact"` root-finds the first-order condition; `"closed"` inverts the
Medlyn USO relation given the marginal cost of water, which is explicit where that
cost is a wet-end power law in `psi`. It exists for `TF24` and `CF77` on the stem
route only, is refused outright with the energy balance on, and falls back to the
exact solve where its validity guard fails — about 2x realised on a mixed driver
grid, at a few percent error in wet soil and none at all on a fallback row.

```r
c <- leaf_model(supply = leaf_supply_singlelayer())
c$CF77_lambda_ <- 1.5e5
set_drivers(c, psi_soil = 0.5, PPFD = 1500, atm_vpd = 1)
c$set_model("CF77", "stem", "closed")
c$optimise()
c$closed_form_fallback_fraction()   # phi, which sets the realised speedup
```

`Leaf()` is also exported: it is the raw C++ constructor, fifteen positional
arguments and no defaults — and four traits it does not take at all.
`leaf_model()` is that with the arguments named, defaulted, split into traits
versus tolerances, and the missing four assigned, and is what you should use.

**All water potentials are positive magnitudes in MPa.** One representation
throughout, and it is asserted rather than documented — a negative `psi_soil` is
an error, not a sign convention the model quietly accepts.

A bare leaf needs no root carbon profile at all — collapse the whole
soil-to-collar path to one resistance:

```r
leaf_solve(psi_soil = 1.5, PPFD = 900,
           supply = leaf_supply_singlelayer(),
           root_network = series_resistance(1e3))
```

The path is chosen when the leaf is built and cannot be flipped afterwards: a
settable tag would leave the other path's state configured and silently ignored.

On the **multi-layer** path the leaf takes the per-layer resistances, so a caller
with measured or fitted ones can state them directly:

```r
l <- leaf_model()
set_drivers(l, psi_soil = 1.5,
            root_network = RootNetwork(r_R_H_min = 25.5, r_R_V_sum = 1410))
```

⚠️ The default `root_network` is a nominal 20 kg C m^-2 leaf put through
`root_network_from_carbon()` — a stand-in rather than a recommendation. It is
written out in `set_drivers()`' body so it can be seen and replaced.

### Trait gradients

`leaf_gradient()` gives the derivatives of the solved outputs with respect to the
traits, which is what a gradient-based optimiser or a Hamiltonian sampler wants:

```r
g <- leaf_gradient(psi_soil = 2.0, PPFD = 900,
                   pars = c("vcmax_25", "stem_P50", "TF24_cost_scale"))
g$gradient   # rows: parameters.  columns: A, gc, psi_stem, collar, profit
g$method     # "ift" or "fd" -- see below

# any of the seven models, same call -- and `leaf_solve()` takes `model` too
leaf_gradient(psi_soil = 2.0, PPFD = 900, model = "JS22",
              supply = leaf_supply_singlelayer(), pars = c("vcmax_25", "JS22_gamma"))
```

The first four columns are what a gas-exchange calibration observes. `profit` is
the objective, and it is there for a demographic consumer: `plant` bills carbon
from the leaf's profit rather than its assimilation, so without that column no
gradient from this package reached a demographic model at all.

`pars` is not restricted to traits: `leaf_specific_conductance_max` and, on the
single-potential path, `resistance` are differentiable too, because a calibration
fits them and nothing in the derivation cares whether a parameter is a trait.

```r
leaf_gradient(psi_soil = 1.5, PPFD = 900,
              supply = leaf_supply_singlelayer(),
              root_network = series_resistance(1e4),
              pars = c("leaf_specific_conductance_max", "resistance"))
```

These are not finite differences of the solve. The outputs are evaluated at the
profit-maximising collar potential, so a trait moves them both directly and by
moving that optimum — and for `TF24_cost_scale`, `TF24_beta2`, `stem_P50` and
`stem_c` the second route is **100%** of the answer. Differentiating the optimality
condition rather than the solved output gets both terms exactly.

`profit` is the exception, and it is the cheapest column for the reason it is the
exception: it *is* the objective, so at an interior optimum the second route
contributes nothing and its gradient is the direct partial alone. That is the
envelope theorem, and the only place this package uses it.

If your model **tracks** the optimum instead of finding it, pass the collar
potential you are operating at and the derivation simplifies rather than breaks:

```r
leaf_gradient(psi_soil = 2.0, PPFD = 900, pars = c("vcmax_25", "stem_P50"),
              psi = 2.7)                # evaluate here, do not solve
```

`psi` is exogenous, so the answer is the partial at fixed collar — plus whatever
`dpsi_dtheta` you supply, if the collar you imposed itself moves with the traits.
`M`, `H` and `dY_dpsi` come back so a caller integrating its own sensitivity of
`psi` has the coefficients. This is what plant's TF24f needs, and giving back the
collar the solver found, with the response it derived, reproduces the solving
path exactly.

That derivation assumes the optimum is interior, and at the dry end it often is
not: with the optimum pinned to the edge of the feasible range the formula returns
a confidently wrong number, off by up to seven orders of magnitude. So the
assumption is **tested** at every point and the function falls back to
differencing the solve where it fails. `g$method` reports which route ran and
`g$status` reports why.

Whether this is *faster* than letting your optimiser difference the objective
depends on your parameterisation, and the two counts that decide it are easy to
conflate. Differencing costs `2 ×` the number of parameters **the optimiser is
moving**; this costs one pass plus a term in the number of parameters **the leaf
has** — `length(pars)`. They are equal only if you fit traits directly. Pooling, a
hierarchy, or any derived parameter makes the first much larger than the second,
which is where this route wins; `vignette("fitting")` measures both regimes and
`?leaf_gradient` has the cost model. ⚠️ **Always pass `pars`** — the default is all
sixteen on the multi-layer path, which is the most expensive request there is.

**For a fit, use `leaf_gradient_batch()`.** It is the same gradient, composed in C++
and vectorised over observations, so a likelihood evaluation crosses the R boundary
once instead of 112 times per observation — **363 → 10.6 µs per observation** at four
differentiated parameters (`length(pars)`), 22×.

```r
b <- leaf_batch(psi_soil = obs$psi_soil, PPFD = obs$PPFD)   # once per fit
g <- leaf_gradient_batch(b, traits, pars = c("vcmax_25", "stem_P50"))
g$gradient   # [observation, parameter, output]
g$status     # per observation: "interior", "pinned", "no-gradient" or "error"
```

The likelihood and your parameterisation Jacobian stay in R, vectorised over
observations: the likelihood is your model, and the chain rule belongs where the win
is — this returns `dY/dθ` for the four parameters the leaf has, and you map your own
onto them. ⚠️ Nothing about the gradient got faster; the boundary was removed from
under it, so the figure needs a batch to be realised.

To vary traits yourself, `set_traits()` replaces them on an existing leaf — much
cheaper than rebuilding one, and the only correct way to do it, since a trait
change invalidates derived state that is not obvious from the outside:

```r
l <- leaf_model()
set_traits(l, leaf_traits(vcmax_25 = 120))
set_drivers(l, psi_soil = 2.0, PPFD = 900)   # required: the drivers must be re-set
l$optimise()
```

See `vignette("phylloptim")` for the whole tour.

### Performance from R

The **~3 µs** quoted at the top of this file is the C++ solve. From R the same
solved row costs about **20 µs**, and the difference is not the model — it is that
each call across the R boundary costs ~1.1 µs, and a solved row needs a handful of
them. Measured on 32 rows, one driver combination per row, default multi-layer
supply:

| | µs per row |
|---|---:|
| `leaf_solve()`, vectorised | **21.5** |
| `leaf_model()` once, then `set_drivers()` + `$optimise()` + `operating_point()` per row | 20.3 |
| the same, reading one field instead of `operating_point()` | 17.1 |
| the C++ solve inside all three | **2.8** |

Three things follow, and the first two are corrections to advice this file used to
imply:

- **Use `leaf_solve()`.** It is within 6% of driving the object by hand. It was
  **26× slower** until it stopped building a one-row `data.frame` per row and
  rbinding them (#39) — 344 → 21.5 µs — so if you are reading advice anywhere that
  says to avoid it for inner loops, that advice has expired.
- **Reaching into the object is not the lever either.** The stateful interface is
  for when you want intermediate state, not for speed; it saves ~1 µs a row.
- **The lever is making fewer R calls per row.** ~18 of the 20 µs is R call
  overhead and R-side assembly. A loop that solves the same leaf at many drivers
  should pass them all to one vectorised `leaf_solve()` call rather than looping in
  R, and if you need a fit's inner loop faster than this, the thing to remove is
  the boundary — which means C++, not better R.

Two costs worth knowing because they surprise people:

- **Constructing a `Leaf` from R costs about 180× a trivial `.Call`** — 45 solves — and only ~32 µs of
  that is the two vulnerability splines; the rest is R-side object construction
  over ~60 active bindings. So construct once and reuse. `leaf_solve(reuse = TRUE)`
  is the default for this reason, and [`set_traits()`](#trait-gradients) exists so
  that a trait sweep need not reconstruct either.
- **`set_traits()` is ~0.02 µs unless you change `stem_P50`, `stem_c`, `root_P50`
  or `root_c`, and 21.8 µs if you do**, because those four own the pre-integrated
  vulnerability splines and it rebuilds one. That is 8× a solve, in C++, where
  batching cannot help — worth knowing before writing a sweep over a vulnerability
  curve. Most of it is the incomplete gamma function seeding 101 knots, not the
  spline machinery. `leaf_gradient()` sidesteps it for `stem_P50`, because the
  curve is homogeneous of degree 1 in the scale that pair implies: see
  `fast_stem_curve` in `?leaf_gradient`.

### As a dependency of another R package

Name it in `LinkingTo` to compile against the headers, the way `BH` is used:

```
LinkingTo: BH, odelia (>= 0.2.0), phylloptim (>= 0.1.0)
```

`LinkingTo` is **not** transitive in R, so you must name `BH` and `odelia`
yourself even though it is `leaf` that includes them — including the odelia
version, for the same reason. A `LinkingTo` consumer gets `<phylloptim.hpp>`, which is
R-free; `<phylloptim.h>` is the R binding layer's own umbrella and is not for you.

## Dependencies

Deliberately few. The two the *model* needs are header-only:

| | why | how |
|---|---|---|
| **odelia** (>= 0.2.0) | cubic-spline interpolator for the pre-integrated vulnerability curves, and the vendored **XAD** automatic-differentiation library | `LinkingTo` |
| **BH** (Boost) | TOMS748 root finder, incomplete gamma for the closed-form vulnerability integral | `LinkingTo` |

**Rcpp** and **R6** are needed by the R layer only. They are not in the model's
include graph and a C++ or Python consumer never sees them.

Nothing else, and **neither model dependency needs R**. The leaf model itself does not
touch **Rcpp** or the R C API: `leaf/util.hpp` replaced plant's `util::stop`
with a plain `std::runtime_error` and `NA_REAL` with a quiet NaN. odelia's
solver core was the last R touchpoint in the include graph, via `ode_util.hpp`;
that was removed upstream in traitecoevo/odelia#44, so the test suite now builds
against the real headers with nothing standing in for R at all. **odelia 0.2.0 is
the first release with that fix**, hence the version requirement: an older odelia
would otherwise fail deep in the build with `RcppCommon.h: No such file or
directory`, which does not point at the cause.

Both dependencies are already required by plant, so plant pays nothing new for
depending on this package.

## Tests

Two suites, and the C++ one is the regression baseline.

```sh
make -C tests/cpp        # plain C++: no R, no test framework
```

It discovers BH and odelia through `Rscript` if R is installed, and otherwise
falls back to a sibling `odelia/` checkout and Homebrew Boost. Override with
`make BH_INC=... ODELIA_INC=...`. `ctest --test-dir build` runs the same two
programs through CMake.

At its centre is `tests/cpp/golden/operating_points.tsv`: 288 operating points
recorded at full precision and compared **bit-exactly**, which is what makes a
large refactor of this code checkable rather than hopeful. It is bit-exact on the
platform that generated it (macOS/arm64) and compared with per-field tolerances
elsewhere, because libm's `exp`/`pow` are not bit-reproducible across platforms.

```sh
R CMD check .            # the C++ suite, plus the R layer's own tests
```

`R CMD check` runs the C++ suite compiled with R's own configured compiler
against the installed headers — so a package that `LinkingTo`s this one finds out
from its own check when a header stops compiling. It also runs
`tests/testthat/`, which ties the R layer back to the same golden points. That
tie-back matters more than it looks: the C++ suite never loads the R layer, so a
mistranslation in the bindings would otherwise produce a green suite and
plausible R numbers. Those expected values are written as **C99 hex floats**, on
purpose — R's decimal parser is not correctly rounded and returns a value one ULP
off for roughly 18% of full-precision inputs, so decimals there would fail
against a model that is exactly right.

## API documentation

```sh
doxygen        # docs/html/index.html
```

Doxygen for the C++ API and roxygen for the R one. The headers' comments are the
substantive documentation here, and `tools/doxygen_filter.awk` presents them to
Doxygen without modifying a single source file.

## How this compares to other leaf models

See [COMPARISON.md](COMPARISON.md) for a feature-by-feature comparison against
`plantecophys`, `bigleaf` and `tealeaves`. The short version: those packages are
stronger on empirical stomatal models, leaf energy balance and fitting to
measured data; this one is the only one with an explicit hydraulic architecture
and a profit-maximisation solve, and the only one written to be embedded in a
larger model.

## Licence

AGPL (>= 3), inherited from plant.
