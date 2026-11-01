# Developer guide — phylloptim

A header-only C++ leaf gas-exchange model in which stomatal behaviour emerges from hydraulics, with an R package on top. Extracted from the TF24 strategy in [traitecoevo/plant](https://github.com/traitecoevo/plant).

The guide is three files. **This one holds the standing rules** and is read at the start of every session. **[hazards.md](hazards.md)** holds the nineteen numbered hazards as narratives with their measurements, plus the long form of "One solver". **[measurements.md](measurements.md)** holds the golden files' cross-platform record, the floating-point findings that recur, and the benchmark method. A rule appears here once, as a line with its reason; its evidence lives there.

Read alongside:

- **[README.md](../README.md)** — what it is and how to use it; its "Performance from R" is the one place the R-side and set_traits costs are quoted
- **[NEWS.md](../NEWS.md)** — what changed in each release, and what it broke
- **[`vignettes/the-models.Rmd`](../vignettes/the-models.Rmd)** — the nine optimality models, the cost/link framework that unifies them, the solver's per-call costs, and two appendices: **B** derives the gradients, **C** is the numerical practice (finite-difference steps, golden files, cross-platform magnitudes)
- **[`vignettes/fitting.Rmd`](../vignettes/fitting.Rmd)** — ⚠️ read *"Fitting a different collection of parameters"* before answering anything about `pars`: what a different collection costs, which parameters return exactly zero and why, which are unreachable, and what a new study has to write itself
- **[COMPARISON.md](../COMPARISON.md)** — how this differs from `plantecophys`, `tealeaves`, `bigleaf`
- **[PLAN.md](../PLAN.md)** — the work queue after 0.27.0, in order, with what "done" is for each item; read it before picking up anything from the issues
- **[issues](https://github.com/traitecoevo/phylloptim/issues)** — the work queue, and the *why* behind each open item

Family context lives in [`plant-meta`](https://github.com/traitecoevo/plant-meta).

## Start here

| To | Read |
|---|---|
| find which header holds what | the map below, in dependency order |
| change the solve | "One solver" below, then hazards [3](hazards.md#3-argmax-smoothness-is-a-hard-constraint-not-a-preference), [8](hazards.md#8-every-path-out-of-the-solve-must-write-its-own-rates), [11](hazards.md#11-a-bracketing-optimiser-answers-where-is-the-interior-maximum-which-is-not-where-is-the-maximum) and [18](hazards.md#18-a-shut-leaf-is-the-ordinary-evaluator-at-its-wet-bound-and-nothing-writes-its-outputs-by-hand-182) |
| add a cost curve | "One solver" below, then [its long form](hazards.md#one-solver-the-long-form) |
| change the derivative surface | the top of `leaf_derivative_surface.hpp`, `vignettes/derivative-surface.Rmd`, hazard [17](hazards.md#17-the-derivative-surface-fails-silently-in-three-ways-and-none-is-a-compile-error) |
| change a trait or the parameter set | hazards [10](hazards.md#10-changing-a-trait-is-not-one-assignment-and-the-reason-is-invisible), [12](hazards.md#12-n_theta-is-a-compile-time-constant-and-a-short-aggregate-initialiser-is-legal-c) and [13](hazards.md#13-the-trait-vector-is-bound-positionally-in-four-places), and grep `.github/` |
| move a golden number | "Golden files" below, then [measurements.md](measurements.md#the-golden-files) |
| time anything | "Cost" below, then [the benchmark method](measurements.md#benchmark-method) |
| build and test | "Build and test" |
| add or change an equation, or bind a new field | "The model registry" below |

**Map of the headers**, `inst/include/phylloptim/`, in dependency order:

| header | holds |
|---|---|
| `util.hpp`, `constants.hpp` | R-free `stop()` and sentinels; physical constants |
| `clamp_sites.hpp` | the `clamp_site` enum: the one bound a solved path keeps, the stomatal corner `E_s >= 0`, counted |
| `quadrature.hpp`, `uniroot.hpp`, `optimize.hpp` | adaptive Simpson; Boost root-finders; `maximise_over_closed_interval_foc`, the one solver |
| `vulnerability.hpp` | `VulnerabilityTable` (the C2 table both curves read), the closed-form integral and its inverse, `vulnerability_derivatives_at` |
| `closed_form_rows.hpp` | a table's value with the closed form's trait rows, for both curves |
| `roots.hpp` | `RootNetwork`, `root_network_from_carbon`, `SupplyInputs`, `CollarConductance`, `MultiLayerRoots` (the one uptake walk), `layer_thicknesses()` (soil depths to layer widths). Nothing here calls `root_network_from_carbon`: since #33 the caller runs the root model and hands `set_physiology` the resistances |
| `single_potential.hpp` | `SinglePotential`, the other supply path |
| `damage.hpp` | `namespace damage`: the thermal-damage dynamics, written once on the scalar type; no `Leaf` |
| `leaf.hpp` | the `par_*` enumeration and `leaf_pars`; the class `Leaf`, with `CostCurve`, `OperatingPointKind`, `with_curve()`, `supply_draw_at` and every `switch (supply_kind_)`; its constructors, `set_traits`, `set_physiology`, `set_model` and `optimise()` |
| `leaf_photosynthesis.hpp` | the temperature responses, the Farquhar kernels, `solve_ci`, and `set_leaf_states_rates_from_psi_stem` (the evaluator at a stem potential) |
| `leaf_hydraulics.hpp` | the residual conductance and its soft cap, the stem curve's four operations, transpiration |
| `leaf_energy_balance.hpp` | leaf temperature, the deficit at the leaf, the boundary layer, `dprofit_energy_balance_term` |
| `leaf_damage.hpp` | `Leaf`'s side of the thermal-damage layer |
| `leaf_cost_curves.hpp` | each curve's cost, profit and λ, ProfitMax's normalisers, `curve_name()`, and every per-curve table |
| `leaf_collar_solve.hpp` | the collar solve: continuity, the wet bound, `maximise_profit_over_collar`, `dprofit_at_collar_psi`, `replay_operating_point_for` |
| `leaf_stem_routes.hpp` | `dprofit_dpsi_stem`, `evaluate_psi_stem`, `optimise_psi_stem_single` |
| `leaf_derivative_surface.hpp` | the derivative surface's contract, then `operating_collar`, `billed_outputs_at` and what only they use |
| `leaf_medlyn.hpp` | the Medlyn model, which is not a `CostCurve`, and `g1_eff` |
| `leaf_model.hpp` | `Leaf` whole: the ten headers above, in order. plant and `gradient.hpp` include it |
| `closed_form.hpp` | `Leaf::optimise_closed`, out of line because it needs the complete `Leaf`; `<phylloptim.hpp>` includes it and `leaf_model.hpp` does not |
| `gradient.hpp` | `n_theta`, the IFT composite and `batch`, the bit-for-bit twin of `R/gradient.R` |

**A member is defined in its subsystem's `leaf_*.hpp`, and a consumer includes `leaf_model.hpp`, never a part.** A part alone parses, but a member it calls and does not define is a link error. Order within `leaf_model.hpp` does not matter, because the templates are instantiated at the end of the translation unit; the two solves include `leaf_cost_curves.hpp` only so that each parses alone. Moving a definition between parts changes no code: Apple clang's objects stay byte-identical, and gcc's differ only in how it numbers its constant labels (#216).

`inst/include/phylloptim.hpp` is the R-free umbrella: one include is the whole model. `phylloptim.h` and the generated `RcppR6_*.hpp` are the R layer's, include Rcpp, and only `src/RcppR6.cpp` includes them; the `.h`/`.hpp` split is forced by RcppR6, which hardwires both names.

## Layout

- `inst/RcppR6_classes.yml` — the R bindings' source of truth. Edit it, then `Rscript -e 'RcppR6::RcppR6()'`; never edit the generated `R/RcppR6.R`, `src/RcppR6.cpp` or `RcppR6_*.hpp`.
- `src/`, `R/` — the R layer. Depends downward on `inst/include/`; nothing there depends back on it. `R/leaf-model.R` is the friendly surface (`leaf_traits`, `leaf_control`, `leaf_model`, `set_drivers`, `leaf_solve`). `R/gradient.R` is `set_traits()` and `leaf_gradient()` (their costs are in README's "Performance from R") and ⚠️ **the reference implementation**: `gradient.hpp` must reproduce it bit-for-bit, so its arithmetic order is load-bearing — no reassociation, no fused multiply-add (see `rounded()`), and every difference names its halves first because C++ does not sequence `f(a) - f(b)` and R does. `R/gradient-batch.R` is `leaf_batch()` and `leaf_gradient_batch()`, the same gradient over N observations in one crossing.
- `tests/cpp/` — the plain-C++ suite, no R, no framework. Five programs: `test_leaf` (unit), `test_supplied_rows` (each layer's collar slope against a difference), `test_transpose` (forward tangent against the reverse sweep; the one that links odelia's `Tape.cpp`), `test_golden` and `test_primitives` (the baselines). `bench_solve` and `bench_gradient` are the timing harnesses; `probe_tape_regions.cpp` prices a design and guards nothing. `root_network.hpp` is the suite's root-architecture fixture, in one place because the golden files' bit-exactness depends on it.
- `tests/cpp/golden/` — three bit-exact baselines. `operating_points.tsv` (672 solved points) says *whether* anything moved and is blind to every stem route. `primitives.tsv` (571 values in five call-tree tiers) says *which* function moved: the lowest tier that changed is the cause, so read `test_primitives.cpp`'s tier table first on a failure. `psi_stem_optima.tsv` (5760 rows) covers the nine stem routes and records known defects on purpose: read its header before reading a row as intended behaviour.
- `tests/testthat/` — the R layer's tie-back to the golden points. `gradient_golden.tsv` holds twenty recorded trait gradients in hex floats, the guard the C++-versus-R equality test cannot be (a change applied to both passes that one); regenerate with `tools/gradient_golden.R` on macOS/arm64, only deliberately.
- `tests/validate/` — the SCM regression pair (needs R and a plant build), and `tsv_to_hex.c`, which `test-golden.R`'s regeneration recipe pipes through because R's decimal parser is not correctly rounded. ⚠️ Never write an R-side expected value by pasting a decimal from a golden file.
- `tools/` — `bench_history.sh`, `bench_user_cost.R`, `cost-baseline.tsv` (the recorded cost table), `gradient_golden.R`, `doxygen_filter.awk`.
- `CMakeLists.txt` — the no-R build for C++ and Python consumers, which makes "does not need R" runnable.
- `inst/model-registry.yml` — the model registry: the leaf's equations, parameters and code symbols, which documents are generated from. Its header holds the schema and the rules; see "The model registry" below. `R/registry-locations.R` reads it and resolves its symbols in Doxygen's XML, and is sourced standalone by `docs.yml`, so it may reference nothing else in the package. `R/registry.R` is the registry as data frames, with the parameters' values read off a live `Leaf`. `tools/registry_symbols.R` is what `docs.yml` runs after the render.

## Build and test

```sh
make -C tests/cpp                    # builds and runs the five programs
make -C tests/cpp bench              # bench_solve and bench_gradient (not in `make all`)
make -C tests/cpp smoothness         # #183's acceptance sweep alone; run it on every PR that touches the solve
make -C tests/cpp golden             # regenerate operating_points.tsv -- deliberately only
make -C tests/cpp psi-stem-golden    # regenerate psi_stem_optima.tsv -- same
make -C tests/cpp primitives-golden  # regenerate primitives.tsv -- same
make -C tests/cpp GOLDEN_ARGS=--cross-platform   # tolerance mode, for platforms other than macOS/arm64
R CMD check .                        # the C++ suite via tests/cpp.R, plus tests/testthat/
doxygen                              # the C++ API, to docs/html/index.html
```

Dependency headers are found via `Rscript` if R is installed, else a sibling `odelia/` checkout and Homebrew Boost; override with `make BH_INC=... ODELIA_INC=...`. CI runs gcc and clang on Linux and clang on macOS.

- **Before pushing, run `make -C tests/cpp && make -C tests/cpp bench`.** CI builds the bench programs and `make all` does not, so a rename that misses a bench source passes locally and fails three CI jobs.
- **CI compiles a third C++ program that is not in the source tree.** `cpp-tests.yml`'s "Consume the installed package" step writes `consumer/main.cpp` in a heredoc and builds it against the installed package through `find_package`. It calls `set_physiology`, so a signature change breaks it; and it fills a `theta[gradient::n_theta]` by hand, so a change to the parameter set breaks it *silently* (a short initialiser built, ran and printed `dA/dvcmax = nan`). **Grep `.github/` whenever the parameter set or a public signature changes.** Reproduce it locally by extracting the two heredocs from the workflow rather than retyping them:

  ```sh
  cmake -B /tmp/bc -DPHYLLOPTIM_ODELIA_INCLUDE_DIR=<odelia>/inst/include
  cmake --build /tmp/bc -j3 && cmake --install /tmp/bc --prefix /tmp/cmp
  # then write consumer/{CMakeLists.txt,main.cpp} from the workflow and:
  cmake -B b -DCMAKE_PREFIX_PATH=/tmp/cmp \
    -DPHYLLOPTIM_ODELIA_INCLUDE_DIR=<odelia>/inst/include \
    -DPHYLLOPTIM_BOOST_INCLUDE_DIR=<BH>/include
  cmake --build b && ./b/consumer
  ```

  The two `-D` include paths are needed again at the consumer step: `phylloptimConfig.cmake` re-requires them and fails with "phylloptim needs odelia's headers", which reads like a broken install and is not.
- **The no-R build is a claim CI keeps honest.** `src/` and `R/` include downward into `inst/include/`, never the reverse, and `cpp-tests.yml` enforces it by building on runners with no R installed (`cmake -B build -DPHYLLOPTIM_ODELIA_INCLUDE_DIR=../odelia/inst/include && cmake --build build && ctest --test-dir build`). It also checks that the install rules ship `*.hpp`, exclude `phylloptim.h` and `RcppR6_*.hpp`, and give a usable `phylloptim::phylloptim`. **Never add an R step to that workflow**; R work goes in `R-CMD-check.yml`, because a job with R present cannot tell the difference.
- **A debug build failing `test_golden` by ~1e-15 has found nothing.** Bit-exactness holds at `-O1` to `-O3`, not at `-O0`, which does not contract FMAs; `CMakeLists.txt` forces `Release` for this reason ([measured](measurements.md#optimisation-level)).
- **`rm -f src/*.o src/*.so` before reinstalling after editing `inst/include/`.** R does not track header dependencies, so `R CMD INSTALL` reuses a stale `src/RcppR6.o` and the R layer goes on running the old model with plausible numbers. Then check one value against the C++ suite, which does track headers.
- **`tests/cpp.R` is not a duplicate of the CI workflow.** It compiles the C++ suite with R's configured compiler (`R CMD config CXX20`) against the *installed* headers, which is what a `LinkingTo: phylloptim` consumer gets. It compiles the sources directly rather than calling `make`, and `tests/cpp/Makefile` is `.Rbuildignore`d, because `R CMD check` warns about the GNU extensions the harness uses; do not add `SystemRequirements: GNU make`, which would be false for every consumer.
- **`'LinkingTo' field is unused: package has no 'src' directory` from `R CMD check` is a failure**: the R layer has stopped being built.
- **`man/` is generated from generated code**: edit `inst/RcppR6_classes.yml` or the roxygen blocks in `R/`, run `Rscript -e 'RcppR6::RcppR6()'`, then `roxygen2::roxygenise()`. Every roxygen block is markdown, so `DESCRIPTION` must keep `Roxygen: list(markdown = TRUE)`; if a regenerated `man/` loses `\code{}` everywhere, that line has gone missing.
- **RcppR6 is deliberately not a declared dependency.** The generated files are committed and RcppR6 is not on CRAN, so it is installed by exactly one CI job (`R-CMD-check.yml`), which regenerates and diffs because a committed generated file goes stale silently.
- **Doxygen reads the headers through `tools/doxygen_filter.awk`**, which turns plain `//` comments into Doxygen ones at render time without touching the sources; CI asserts every non-comment line survives byte-for-byte. A new header needs nothing; a `///` comment is left alone, for a literal Doxygen command. ⚠️ **One indented display per comment block**: Doxygen 1.9, which CI installs, drops the second `\verbatim` in a block and reports an error 100+ lines away in unrelated code, while a local Doxygen 1.17 renders it in silence. `docs.yml` asserts the shape and names the block.
- **`docs.yml` builds two sites and publishes them as one**: the pkgdown site (`_pkgdown.yml`, rendered to `docs/site`) at the root and the Doxygen render (`docs/html`) under `cpp/`, because a repository has one Pages site. Publishing is off until the repository variable `PUBLISH_DOCS` is `true`; until then each run uploads the `pkgdown-site` and `api-docs` artifacts. A new exported topic must be listed in `_pkgdown.yml`'s reference index or the pkgdown build fails; `pkgdown::check_pkgdown()` checks that locally in seconds.
- **`README.md` is guarded by `tests/testthat/test-readme.R`**, which parses and evaluates every ```` ```r ```` block. It guards syntax and API only; the numbers in the surrounding prose rot without it noticing.

## One solver. Do not add a second.

**Every optimality model reaches its operating point through `util::maximise_over_closed_interval_foc(f, df, lo, hi, n, tol, iters, &fmax)`**: endpoints, so a constrained optimum is reachable; an `n`-cell scan, so a second basin cannot hide the global one; a root-find of `df == 0` inside the winning cell; and the grid argmax always kept, so refinement can only improve on it. Four duplicates of it have been written and deleted, and deleting each made the code smaller and faster; two were also wrong in a way the one solver was not.

Everything is an argument to that call:

| you might think you need | you need |
|---|---|
| a different cost curve | a `CostCurve` enumerator and its table rows (below) |
| a different objective *shape* (a product, a normalised ratio) | a `BenefitLink` — the cost goes in log form and `h'` carries the rest |
| a different decision variable or supply topology | the bracket and the chain factor `dpsi_stem/dpsi` |
| a scan, or no scan | `n` |

- **A new cost curve touches no solver.** It is one enumerator appended last (move `n_cost_curves`, which names the last one), one arm each in `with_curve()` and `curve_name()`, and one arm in every per-curve `if constexpr` table: those are the ones ending in `static_assert(..., "unhandled CostCurve")`, so `grep -n 'unhandled CostCurve' leaf_cost_curves.hpp closed_form.hpp` lists them and a missed one is a compile error. Today they are `benefit_link`, `check_cost_parameters`, `cost_deriv`, `cost_reads_upstream`, `cost_domain_ends_at_psi_crit`, `objective_is_product`, `profit_psi_stem_for`, `lambda_for` and `lambda_collar_for`, plus `cost_deriv_upstream` for a curve that reads the upstream potential, and `closed_form::optimise_into`, which refuses a curve with no closed form.
- ⚠️ **A cost that reads the upstream potential needs its `cost_reads_upstream<K>()` and `cost_deriv_upstream<K>()` rows**, or the collar solve converges onto a point that is not the maximum, with no error (#141).
- **A cost or λ that reads the flux reads it through `operating_point_flux()`, never `transpiration()`.** `transpiration()` refuses a reversed pair (#31), while the evaluator reads that pair as no flow, and a prescribed collar at the multi-layer wet bound inverts to a stem ~1e-14 MPa wetter; `test_prescribed_wet_bound_reads_no_flow` lists the curves that read E, and a new one belongs in it.
- **Runtime curve dispatch is ONE switch, `Leaf::with_curve()`**, plus `curve_name()`. Neither has a `default:`, deliberately: a `default` satisfies `-Werror=switch` and removes the check they exist for.
- **Test a new curve by a bit-for-bit REDUCTION to an old one** (`TF24_floor` is `TF24` at `lambda_o = 0` and `CF77` at `TF24_cost_scale = 0`) before reaching for a new baseline, and prefer sharing the parent's parameters when the new curve is a parent plus a term.
- **A curve parameter that is a field rather than a trait needs wiring** in `.solve_price_fields`, `.gradient_owned_pars` and `.gradient_setter`, or `leaf_solve()` cannot reach it.
- **`n` is set from measurement: 0 unless the energy balance is on** (`basin_scan_cells()`). Only energy-balance rows have two basins, and a scan by default is a 3× cost on the hot path. A scan is also never a way to find a maximum you will differentiate: its argmax is piecewise constant.
- **Two things look like a second solver and are not**: `Leaf::maximise_profit_over_collar` (the one solver at `n = 0` plus the `OperatingPointKind` classification) and the shared algebra of `dprofit_dpsi_stem` / `dprofit_at_collar_psi` (whose exact operation order keeps the golden files bit-identical). The shape to refuse is a function duplicating the solver's *structure* with a piece missing.

**If you are about to write a new `optimise_*` or `maximise_*`, first say which of the four arguments above you need**; if none, put the measurement that shows your case is different in the PR. The measurements and the deleted duplicates are in [the long form](hazards.md#one-solver-the-long-form).

## Hazards

Each has cost someone real numbers; the narrative and its measurements are in [hazards.md](hazards.md). Code comments cite these numbers, which never change.

1. **There are two Weibull vulnerability curves**, stem (`stem_P50`, `stem_c`) and root (`root_P50`, `root_c`), each with exactly two settable numbers; `stem_b`, `psi_crit` and `k_crit` are derived. Never leave an unmarked parameter that exists in two versions, and never make `stem_b` settable.
2. **Every ψ is a positive magnitude in MPa** (#25), asserted at every entry point. Signs live in equations, not in variables, and the convention is not "tension".
3. **The argmax must vary smoothly with inputs**; any solver change must re-measure it, and the closed form is never under a derivative because its guard switches on an output.
4. **The leaf is purely intensive**: every input is per unit leaf area, including the `RootNetwork` resistances, which nothing here can check.
5. **The collar solve is the hot path**, so measure before arguing about dispatch or inlining: the supply path is out of line, the cost core is fully inlined.
6. **There are two supply paths, selected by `supply_kind_`**; a third needs an arm in every `switch (supply_kind_)` (each has a `default:` that silently means "single"), and `duptake_dpsi` returning NaN means "use finite differences".
7. **Moving a public member breaks plant's generated bindings**; fix it with a dotted `name_cpp` and land both changes together.
8. **Every path out of the solve writes every output**, because plant reuses one `Leaf` per patch; a bit-identical golden run is no evidence about stale state.
9. **Nothing reachable from `<phylloptim.hpp>` includes Rcpp**, enforced only by `cpp-tests.yml` running without R.
10. **Change traits with `set_traits()`**, never by assignment; three pieces of derived state depend on them, and the temperature cache's guarantee lives in `photo_temp_key()`.
11. **A bracketing optimiser finds the interior maximum, not the maximum**; the stem routes maximise over a closed interval, and a stem route is a different model from the collar solve.
12. **`theta` initialisers are sized by `n_theta`**, and a short one compiles; grep `n_theta ==` and update all three asserting places together.
13. **The trait vector is bound positionally in four places**; append, never reorder, and keep `.gradient_owned_pars` the one owner table.
14. **The P50 chain rule holds only at interior optima**; key any assertion of it on `status`.
15. **A finite difference through a solve has a noise floor that moves with the solver**; sweep the step and read the floor.
16. **Print how many rows reached an assertion**; a test whose rows all pinned passed without evaluating anything.
17. **The derivative surface fails silently**: unvisited `for_each_active` members, unchecked seating assumptions, and a replay that must go through `replay_operating_point`.
18. **A shut leaf is the ordinary evaluator at its wet bound**, where uptake meets the cuticle's residual flux (#182): append residual terms, never interleave them, so `g_min = 0` stays bit-identical; never replace a soft cap with a `min()`; and list every active input a wet-bound `implicit_value` reads.
19. **A bound is decided on the bound itself, solved to the precision it is reported at**; decided a step inside it, the solve jumped where the decision flipped and a difference of a whole run integrated the jumps into an offset no step size removed.

## Golden files

- **Only regenerate deliberately**, and say in the PR what moved and by how much. Regenerating after an accidental change rubber-stamps it.
- **Bit-exact on macOS/arm64 only**; elsewhere CI compares with `--cross-platform`. **Never regenerate to make another platform pass** — it moves the failure to the platform the file came from.
- **Magnitudes: ~1e-16 is reassociation, ~1e-12 is the solver floor, ~1e-4 is a real difference.** The floor is set by `psi_stem_to_ci`'s 1e-13 tolerance; do not loosen it without reading [why it is there](measurements.md#magnitudes-what-a-difference-means). The settable `ci_abs_tol` reaches only the Medlyn route.
- **Read the summary line from master's CI**, not the FAIL lines and not a table in a document; a worst *relative* row is usually a near-zero value, so compare it with the worst *absolute* row.
- **Check a collar-solve change by the residual `|dprofit|`, not by profit**, which is flat at its maximum.
- **Before demanding bit-equality across R and C++, match accumulator width and order**: R's `sum()` is `LDOUBLE` ([measured](measurements.md#rs-sum-is-not-cs-)). **Reversing a bracket costs the last two bits** in TOMS748 ([measured](measurements.md#toms748-is-not-sign-symmetric)).

## Cost: the one thing the golden files cannot see

`tests/cpp/bench_solve.cpp` and `bench_gradient.cpp` cover the C++ side; `tools/bench_user_cost.R` covers the R side; `tools/bench_history.sh` runs both against a list of commits; `tools/cost-baseline.tsv` is the recorded table and `tests/testthat/test-cost.R` the guard. The README's "Performance from R" quotes the R-side figures and the solver table in `vignette("the-models")` the per-model ones; everything else links to those two. How to measure is in [the benchmark method](measurements.md#benchmark-method).

1. **Never pin an absolute microsecond figure, in a test or in prose.** Record cost ÷ a trivial `.Call` in the same process; that ratio travels between machines and runs. Quote each figure in one place and link to it from the others, or the copies drift apart.
2. **Guard by COUNTING, not by timing, wherever you can.** A count of `.Call`s has no variance. Timing assertions are a backstop, bounded at ~2× the measurement and `skip_on_ci()`. Count the symbol the path actually calls: the two supply paths reach C++ through different ones.
3. **For the batch gradient a count is always available**: one `gradient_batch_run` per call for any N, and zero per-perturbation primitives.
4. **Measure the entry point the USER calls.** A regression in `leaf_gradient()` once went unseen while the two entry points measured were flat.

## The model registry

`inst/model-registry.yml` is the one place the leaf's equations are written down for documents: each with its LaTeX, the C++ symbol that evaluates it, the parameters it reads and the bound fields it writes; and each parameter with its symbol, units, class, meaning, source and placeholder flag (#166). It covers every bound field: a parameter, an equation's `writes`, or an entry under `diagnostics` (a bound field that is neither; none is needed today). A document's equation list or parameter table is generated from it rather than typed beside it. `model_registry()` reads it, `registry_values()` adds the values, `registry_resolve()` adds file, line and API link.

Three rules keep it from rotting, and each has a mechanism:

- **It stores no line number and no default.** Both go stale silently: the methods document this grew out of cited line numbers that moved twice in one day, and a default's comment that outlived the parameter's rename. Lines come from Doxygen's XML (`GENERATE_XML` in the Doxyfile); values come off a live `Leaf`, `set_drivers()`'s formals or a constexpr's initializer in that XML. The one number the file holds is a `literal:`, a constant written inside an equation (Tetens' 17.27), and docs.yml checks that it appears in the body of the function its `code` names.
- **CI validates it.** `test-registry.R` checks every symbol and field named exists and asserts that every bound field is accounted for, so binding a field without registering it fails by name. `docs.yml` resolves every symbol in the XML, checks that what an equation `reads` is what its function references (`REFERENCES_RELATION`), and checks every literal, through `tools/registry_symbols.R`.
- **One registry per package.** plant's will cover its strategies and point at these ids rather than copy them.

**When you change an equation's code, re-read its entry.** Nothing compares the LaTeX with the C++: CI sees a renamed symbol, a changed `reads` and a changed literal, but not a changed formula. The registry is a record, not a check; the numbers are checked by the golden files and the test suites, and an equation that has no test there gets one there, not in the registry.

That is a decision, not an omission. Measured by planting sixteen one-line bugs in the C++, a base-R replica of every equation caught all sixteen and the existing suites thirteen; the three it alone caught now have their own tests, except a wrong slope in the TF24 closed form's fixed Newton steps (`wet_end_p`), which moves psi* by 2e-5 in a method accurate to a few percent. A second implementation of the whole leaf, touching about 237 bindings, costs more to keep in step than that. Do not reintroduce a per-equation replica check without a measurement that says otherwise.

Doxygen's line numbers are in the stream it parsed, which is the awk filter's output, not the source: the filter inserts a `\file` block and `\verbatim` pairs and drops held blank `//` lines. The lines that are not whole-line comments are byte-identical in both streams, docs.yml asserts that, and every line Doxygen reports is one of them, so `registry_line_map()` maps the k-th code line of one to the k-th of the other and refuses a header whose code lines differ. Do not read a line out of the XML without that map.

YAML 1.1 reads `1e-3` as a string: write a float with a point and a signed exponent (`1.0e-3`), which the reader checks.

## Invariants and what not to touch

- **The one-way include rule** (hazard 9; older comments call it "PLAN item 6a"): `inst/include/` is plain C++ and never includes from `src/`; `src/` and `R/` include downward only. A separate R-binding package was considered and rejected: two repositories for one model, installed twice.
- **The generated RcppR6 files** (old comments: "PLAN item 6c"): never hand-edit them, and keep RcppR6 undeclared.
- **`gradient.hpp`'s arithmetic order**, which must match `R/gradient.R` bit-for-bit.
- **`psi_stem_to_ci`'s tolerance** and **`basin_scan_cells()`'s gate**, both set by measurement.
- **No clamp, `min`, `max` or branch on a physical quantity inside the solved range** (#183, decided). A driver or state outside the model's range is refused at its entry point, by name and limit; a reachable bound gets a smooth form with a registered width. The kept exceptions are the stomatal corner `E_s >= 0`, `vpd_leaf_min`, SOX's and JW26's reduction factors held at zero past psi_crit (without it a dark leaf's negative A times a negative factor reports a positive profit), SOX's and JW26's shut for a leaf with A <= 0 at the dry end (`objective_is_product`: a product has no maximum there, and its argmax jumps at the light compensation point whatever is done), and the published forms (Medlyn's beta, CF77's soil price), each with its reason in the registry. `make -C tests/cpp smoothness` is the check.
- **The `supply_kind_` enum**, which measured free where `std::variant` did not (hazard 6).
- **Public `Leaf` members**, which plant binds by name (hazard 7).
- **`tsv_to_hex.c`**, which is live: `test-golden.R`'s regeneration recipe needs it.
- **The aligned level defaults** (`CF77_lambda_`, `JS22_gamma`, `CMax_a`, `LeastCost_beta`; #156) are solved by `align_curves()`, not chosen, and live in four places: the in-class initialisers in `leaf.hpp`, `leaf_traits()`, `.leaf_trait_defaults` (`JS22_gamma`, `CMax_a`) and `.level_defaults` (`CF77_lambda_`, `LeastCost_beta`). `test-align-curves.R` re-solves them, so a change that moves TF24 or one of these curves at the reference fails it; re-align every copy rather than loosening the test. The golden files pin their own values and do not follow the defaults.

Comments in `inst/include/phylloptim.h`, `inst/RcppR6_classes.yml`, `Doxyfile` and old `NEWS.md` entries that cite "PLAN item N" refer to a `PLAN.md` retired in #126; its last version is `git show 83fbb25^:PLAN.md`. Item 6a is the one-way include rule above, 6c the RcppR6 layering, and 15 the Doxygen site (issue #12).

## Validating against plant

**There is no harness that compares against plant, and there cannot be one**: plant consumes these headers, so the comparison would be this package against itself. What you want instead:

- `tests/cpp/test_primitives.cpp` for *where* a difference comes from;
- `tests/cpp/golden/operating_points.tsv` for *whether* anything moved;
- `tests/validate/scm_regression.R` + `scm_compare.R` for the SCM half, which a bit-identical leaf solve does not settle.

If you do run anything against a plant build: build plant from the commit you mean rather than whatever is installed; print `find.package("plant")` in every script and check it, because a cleaned `/tmp` once made a run fall back to the site library with plausible numbers; and test an installed plant with `test_check("plant")`, not `test_dir()`, which gives 93 spurious "could not find function" errors. Pass `TESTTHAT_PARALLEL=false` when the package is in a non-default library, since the parallel workers do not inherit `R_LIBS` (traitecoevo/plant#586).

## Related work

- **plant** — consumes these headers directly (`plant::Leaf` is an alias of `phylloptim::Leaf`); its TF24 strategies are the only consumer of the derivative surface, and replay every recorded point during a reverse sweep.
- **odelia** — supplies the supplied-derivative nodes (`implicit_node.hpp`) and the vendored XAD. The solve's own derivatives are forward mode and need no tape; the rows the leaf hands a consumer are recorded on that consumer's reverse tape, and `test_transpose` links odelia's `Tape.cpp` for that.
- **The companion manuscript** — `Falster-stomatal_analytical_analysis` in atelier, *"The marginal cost of water as a common currency for stomatal optimality models"*. **It is this package's first customer, not a downstream user**: its blockers are ours, and its central result — that six optimality models differ only in λ(state) — is what the cost/link framework is built on.

## Branches

**Single-branch now.** `feature/api-cleanup` merged as #15 on 2026-08-03, and with it the two-branch rule this section used to state ("results changes live on the branch, `master` stays a drop-in"). Both halves of that rule are spent: `master` is no longer a source-level drop-in for plant's leaf, and there is no longer a branch to quarantine results changes onto.

- **`master`** — everything. **Behaviourally** still validated against plant's own leaf as of the pre-#15 state: full suite 0 fail / 0 error, SCM regression bit-identical across 78/78 nodes. (The pass *count* is environment-dependent — 2364 when first recorded, 2431 on a rerun with `NOT_CRAN` unset, on both arms of a control. Compare against a control run, never against a remembered number.) ⚠️ **#15 deliberately changed results**, so that validation now needs redoing against plant — see hazard 7.

What replaces the old rule: **a change that moves results lands on `master` through a PR that states its measured blast radius against the golden file.** #15 is the worked example — 240 cells from the shutdown fix, 355 at ≤1.2e-13 from reassociation, split by cause because the two classes are four orders apart. A PR that regenerates the golden file without saying what moved and why is the thing to refuse.

The commit message carries **one line** of that: the magnitude and where the breakdown lives ("moves 240 golden cells; split by cause in #15"). The cell counts, the per-cause split and the tolerance bands go in the first PR comment — see "Writing commits here" below.

## Writing in this guide

These files are read at the start of every agent session, as current instructions, by a reader who was not in your session and cannot tell which parts are out of date. When you change them:

- Put a standing rule here, as one line with its reason, and its evidence in `hazards.md` or `measurements.md`. Keep this file near 250 lines; when it grows, move narrative out rather than compressing rules.
- Keep hazard numbers stable: code comments cite them. Append a new hazard as the next number, in both this file's list and `hazards.md`.
- Write what is true now, in the present tense. If something here turns out to be wrong, replace it; the story of how you found out belongs in the PR.
- Give each rule its reason, in plain case. A reason persuades where capitals and warning signs do not.
- Point to where a value lives (`DESCRIPTION`, a header, a Makefile target, `tools/cost-baseline.tsv`) rather than copying the value, and leave out anything a fresh clone cannot resolve, such as local paths or notes from one session.
- Write in complete sentences that say what happens and what it causes, using names that already exist in the code.

The same reader meets the code's comments and NEWS, so two more rules hold there:

- **A comment in code states what must hold and what breaks if it does not**: a contract, a unit or a hazard, in one to three lines. A derivation that needs more goes in a vignette or `hazards.md`, linked by name. No history (what a line used to be, when a slot was inserted, which probe sat here), no retracted claims and no issue numbers: the PR holds them. Keep ⚠️ for a mistake the reader is about to make that the hazard list names, and cite the hazard by number. No timings: where a choice is made for cost, say so in words and name the bench that shows it; the figures live at the entry points, in the solver table of `vignette("the-models")` and README's "Performance from R", because a figure per element churns with every change.
- **A NEWS entry is one to five one-line bullets per release**: what changed for a user, whether results move and roughly by how much, any API break with what to do instead, and the PR. A speed change is quoted only for the main entry points (the collar solve, plant's stands); tolerances, row counts, golden-file moves and other timings belong in the PR.

`tools/check_lengths.sh`, run by `cpp-tests.yml`, holds both: a NEWS entry is at most 12 non-blank lines, and in a header under `inst/include/phylloptim/` or an R file under `R/` of 200 non-blank lines or more, `//` or `#` comment lines are at most 40% of them.

## Writing commits here

**State the measurement, not the intention** — that part of the old rule here stands. What has changed is *where* the rest of it goes.

This repo squash-merges, so a PR's title and body are copied verbatim into permanent history. Commit bodies here reached a median of 450 words and a maximum of 1822, and two commits exist for nothing but correcting an earlier commit message. That is the habit to drop: a reader three years from now is bisecting, and needs what the change does, not the order in which it was understood.

So: keep the commit to why the change was needed and what it does to observable behaviour, in about 10 lines. Everything else — the measurements in full, the explanations that turned out wrong, the alternatives rejected, what was tried first — goes in the **first comment on the PR**, which is just as permanent and is linked from the squashed subject. Never write a commit whose purpose is to correct an earlier commit message; correct the record on the issue or PR instead.

The family guideline, with the split set out as a table, is [`commit-messages.md`](https://github.com/traitecoevo/plant-meta/blob/main/governance/commit-messages.md).
