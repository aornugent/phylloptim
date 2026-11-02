# Measurements — phylloptim

The numerical record that the standing rules in [CLAUDE.md](CLAUDE.md) rest on: what the golden files can and cannot say, how far apart the platforms are, two findings about floating point that recur, and how to time anything. The hazards and their measurements are in [hazards.md](hazards.md).

## The golden files

`tests/cpp/golden/operating_points.tsv` records 672 operating points — one 336-point state grid at 25 and 40 °C — and is compared **bit-exactly** (`%.17g` round-trips a double). It is what makes a large refactor of this code checkable rather than hopeful: it proved three "surely dead" `set_physiology` arguments really were dead, confined the shutdown fix to exactly 48 rows × 5 fields, and showed the `area_leaf` change was 2 ULP.

⚠️ **A grid at one leaf temperature is blind to every temperature response in the model**, because every reference value is *defined* at 25 °C and a change to any response curve is inert there **by construction**. Hence the second temperature, and hence the classification being printed per temperature: points move between branches as the leaf warms (boundary-dry 9 → 0, boundary-soil 42 → 116). Temperature is the OUTERMOST loop, which is what lets a regeneration that adds one be checked as an addition — the 25 °C block must come out byte-identical.

The three files answer different questions, and the layout in CLAUDE.md says which: `operating_points.tsv` says *whether* anything moved and is blind to every stem route; `primitives.tsv` says *which* function moved (the lowest tier that changed is the cause); `psi_stem_optima.tsv` covers every curve on both routes, solved and at a prescribed point, records today's known defects on purpose (read its header before reading a row as intended behaviour), and reuses one `Leaf` with poisoned outputs so stale state is visible.

### Magnitudes: what a difference means

**~1e-16 is reassociation, ~1e-12 is the solver floor, ~1e-4 is a real difference.** Anything above the floor deserves investigation.

A last-bit perturbation does not stay small: the nested solvers amplify it up to whichever tolerance is loosest, because a last-bit change shifts *which point* inside its tolerance band a bracketing solver lands on. Measured by applying **one identical perturbation** (unifying the AD replica with the function it mirrors) at each stage of the solver's history:

| amplifier in play | that perturbation shows up as |
|---|---|
| golden section, `GSS_tol_abs` 1e-3 | **5.53e-04** — the old collar solve |
| `psi_stem_to_ci` at 1e-7 | **4.98e-07** — golden section gone |
| `psi_stem_to_ci` at 1e-10 | **~1e-09** |
| `psi_stem_to_ci` at **1e-13 — the current tolerance** | **7.16e-13** |

- **`psi_stem_to_ci`'s tolerance sets what every reported output MEANS**, so the call site (`solve_ci` in `leaf_photosynthesis.hpp`) points here, and the table above is the precision behind it. 1e-13 was chosen over the cheaper 1e-10 for consumers that DIFFERENCE the solve (#172): where TOMS748 stops inside its band depends on the state, so at 1e-10 outputs were noisy at up to 3.7e-11 relative along a smooth path through the drivers, and the IFT gradient's own central differences of `dprofit` (step ~1e-6) inherited it, up to 2.3e-3 on `gradient_golden.tsv`. It costs ~5% per solve and per gradient (+3.4% in the controlled A/B for the step actually taken). Do not loosen it without re-reading these magnitudes: the plateaus in that curve are TOMS748's discrete iterations, so some changes cost without buying anything, and 1e-12 still leaves 1.2e-12 at some states.
- **It is not the settable `ci_abs_tol`** (default 1e-3). `ci_abs_tol` is read in exactly one place, `solve_medlyn_ci_numerical` — the empirical Medlyn-conductance route. **Both** families of optimality solver reach `ci` through `psi_stem_to_ci`, so both get its fixed tolerance and neither is reachable by the settable control. Tightening it buys no precision anywhere in the optimality model.

### Across platforms

**Bit-exact only on the platform that generated the files — macOS/arm64.** libm's `exp`/`pow` are not bit-reproducible between glibc on x86-64 and Apple's libm on arm64, and FMA contraction differs too, so cross-platform bit-equality was never achievable. CI compares bit-exactly on macOS and with `--cross-platform` elsewhere.

**The size of the cross-platform disagreement**, read off CI's summary line, with gcc and clang identical in every column (4116 of 5184 values differ at 0.19.0), against `operating_points.tsv`'s tolerances, now 1e-9 for `profit` and 5e-6 for the eight argmax-evaluated fields (#125; each about 50x the worst below, and set in `test_golden.cpp`):

| class | worst RELATIVE | worst ABSOLUTE |
|---|---|---|
| `profit` | 2.13e-11 | 9.05e-13 µmol C m⁻² s⁻¹ |
| the other eight | 8.6e-08 (`assim`) | 3.78e-11 Pa (`ci`) |

⚠️ **What sets these figures is the loosest tolerance in the solve, not libm.** At `psi_stem_to_ci`'s former 1e-10, both classes landed on ~1e-9 absolute, the solver floor of the day, and tightening it to 1e-13 (#172) took them down 80× (`ci`) and 1100× (`profit`) with no change to the platforms. So when a cross-platform figure moves, look at the solver tolerances before reaching for a libm explanation.

**Read the baseline off master's CI, not off this table** — `gh run list --branch master --workflow cpp-tests.yml`, then grep the g++ job's log for "values differ". It has to be the summary line from the right commit. `test_golden` prints only the first 20 failing lines, a sample biased toward whichever operating points are listed first, so never read the FAIL lines as the distribution; the summary line carries the full-grid maxima.

**`psi_stem_optima.tsv` has its own tolerances (#125)**, measured for it rather than inherited, in two groups. The original ten solvers: 1e-10 for `profit` (worst 3.87e-12) and 1e-7 for the argmax class (worst 1.24e-9, at `fresh/collar/multi3 psi_soil=6 ppfd=1500 T=25 eb=1 tc=1 assim`). The solvers appended in 0.27.3, every other curve's collar route and LeastCost at beta 2: 1e-7 for `profit` (worst 2.42e-9, at `fresh/collar_ProfitMax/multi3 psi_soil=4 ppfd=1500 T=40 eb=0 tc=0`, a normalised profit of -6e-7 whose two terms, ~6e-6 each, nearly cancel) and 1e-5 for the argmax class (worst 1.97e-7, at `fresh/collar_JS22/multi3 psi_soil=6 ppfd=0 T=50 eb=1 tc=0 hydraulic_cost`, a leaf all but shut whose flux is ~1e-10). The prescribed rows appended in 0.27.4 (`at_<curve>` and `collar_at_<curve>`) share that pair and set neither worst. The original rows keep their pair so that a Linux-only change to the routes they cover still shows. Until 0.19.0 that class's worst was 6.63e-06, a shut leaf's near-zero `assim` cancelling: a root-find's residue around zero. Since #182 that `A` is zero by definition and reported as an exact zero, so the class sits at the solver floor like the rest. Its worst *absolute* figure is no guide: it includes `lambda_emergent`, whose values reach 2.4e22. It prints the worst observed difference per class even when it passes.

`tests/testthat/gradient_golden.tsv` is further apart across platforms than the solved outputs (1.19e-7 measured on Linux, against a tolerance of 1e-5), because a finite difference divides the solver floor by the step, `R_d_25`'s step is the smallest here, and `profit` is the smallest magnitude carrying that same floor.

gcc and clang report *identical* figures. They used to differ by 2–3×, and what differed was which way a golden-section comparison fell; with that search gone, what is left is the platform's arithmetic, not the compiler. **If a compiler-dependent column reappears, something has reintroduced a discrete decision into the solve.**

Two things follow that matter beyond these files:

- **`profit` is the best-conditioned reported field** — 2.13e-11 against the argmax class's 8.6e-08, relative. A portable check of the solve should prefer it over `opt_psi_stem_`; re-read the summary line before leaning on the margin.
- **The eight argmax-evaluated fields are determined, not just pinned.** The collar solve solves `dprofit == 0` rather than searching profit, so on the 198 interior rows the argmax is determined to solver precision — `|dprofit|` at the returned collar has a median of **5.6e-15**, against golden section's 7.8e-4. The remaining 42 rows have a **constrained** optimum pinned to a bracket bound, where the gradient is genuinely non-zero and the answer is determined to the step-in scale (~1e-6 of the bracket width) instead. So bit-exactness is a drift detector *and* the numbers are determined.

### Optimisation level

**Bit-exactness depends on the optimisation level too.** Measured on macOS/arm64: bit-identical at `-O1`, `-O2` and `-O3`; **not** at `-O0`, which misses by 3.47e-15 (~13 ULP) because it does not contract `a*b + c` into an FMA. That is in the reassociation band, so **a debug build failing `test_golden` by ~1e-15 has found nothing.** The Makefile defaults to `-O2`; CMake's default is flagless, i.e. `-O0`, which is why `CMakeLists.txt` forces `Release`.

### A worst RELATIVE row is usually a cancellation artefact of the metric

A relative difference divides by the value, so wherever a reported quantity nearly cancels — a dim, hot leaf whose assimilation barely clears respiration — the metric is amplified by the cancellation and says nothing about the model. The worst relative rows have repeatedly been such rows (an `assim` of 4e-7 against a grid median of 1.44), and the worst absolute row is usually a different one.

So the summary reports the worst **absolute** difference per class alongside the relative one, each naming its row, and a cancellation row is told from real drift by comparing the two lines. **The tolerance stays relative**: a near-zero row is exactly where a genuinely broken value would also show up, so the check stays where it is and only the reporting distinguishes the cases. The sibling file `psi_stem_optima.tsv` handles the same problem the other way, with an absolute floor (`kOptNegligible`, 1e-12) below which it declines to compare — the two files answer this differently on purpose. The one exception in `operating_points.tsv` is the water fluxes (`kFluxNegligible`, 1e-15 kg m⁻² s⁻¹, mirrored in `helper-golden.R`): a shut leaf's uptake is the polished wet bound's residue, ~1e-20 of either sign by platform, and no relative tolerance can pass on that. The floor sits in a measured empty band, four decades below the smallest flux an open leaf reports (7.6e-12), so it cannot hide a broken open leaf.

### Profit is the wrong instrument for checking a collar-solve change

It is the maximum, so it is flat, and its own numerical floor is set by the nested `ci` root-find's tolerance. When the root-find replaced golden section, two of 288 rows came out ~6e-7 *lower* in profit while their residual improved by ten orders of magnitude. **Check the residual `|dprofit|` at the returned point instead** — it improved on 240 of 240 feasible rows with none worse, which is a statement profit could not have made.

## R's `sum()` is not C++'s `+=`

**R's `sum()` accumulates in `LDOUBLE`: 80-bit on x86-64, 64-bit on arm64.** C++ accumulating in `double` therefore agrees with `sum()` bit-for-bit on arm64 and can differ by an ULP on x86-64 — for any sum of more than one term. A test comparing them is platform-dependent.

Found by #92, on `test-surface.R`'s check that `Leaf::operating_point_values()` agrees with reading the twelve outputs one binding at a time. The twelfth, `uptake`, is a sum over soil layers, so the R side derived it with `sum()`. The assertion had passed on macOS and Linux for years by luck: #92 moved the values and the three-layer case came apart on ubuntu only, at element [10] of twelve, while macOS stayed green. The code under test had not changed platform behaviour; the assertion had always been platform-dependent.

`Reduce("+", finite, 0)` accumulates left to right in plain double addition from the same zero the C++ loop starts at, so the comparison is bit-exact everywhere. Prefer that to widening a tolerance: what the test exists to catch is a **shifted column**, which is worth catching at the last bit.

**Before demanding bit-equality across the R/C++ boundary, check that both sides use the same accumulator width and the same association order.** `sum`, `mean`, `cumsum` and `prod` all use `LDOUBLE`; `+` on doubles does not. This is the mirror image of R's decimal parser, which is not correctly rounded (hence `tsv_to_hex.c` and hex-float expected values): both times the model was right and the instrument was not.

## TOMS748 is not sign-symmetric

**Reversing a bracket's orientation costs the last two bits, and no amount of care in the model code buys them back.** Budget for it when planning a change that flips a sign, because it comes up whenever something is "supposed to be bit-identical".

#25 predicted the golden file would be bit-identical bar an exactly-negated collar column. It came out **276/288 exactly negated, with 12 rows off by 1–3 ULP** and 5–11 rows per other column moving at ≤5e-16. Localised: the flux rewrite *is* exactly antisymmetric — `root_zero_E`, which depends only on `E_up`, is exactly negated — while `root_crit`, from the same solver on a **reversed bracket**, is not. Reproduced independently on a synthetic family of smooth monotone targets: 44 of 192 bracket midpoints are not the exact negation.

## Benchmark method

**C++: `make -C tests/cpp bench`** builds `bench_solve` (the collar solve over the whole golden grid, and every single-layer arm on its 1-layer subset, all in ONE process so their ratios are a controlled A/B) and `bench_gradient` (a trait gradient: the IFT composite against differencing the solve). **R: `tools/bench_user_cost.R`** measures one solve, N solves, construction and a gradient, each as a ratio to a trivial `.Call` in the same process. **`tools/bench_history.sh <sha>...`** runs both against a list of commits and prints one row per commit in the format of `tools/cost-baseline.tsv`, which is the recorded table. Append a row for every release and for any change that moves a cost, and quote figures from the newest row; the cost gate is that row plus the interleaved bench a PR touching `inst/include/` quotes in its first comment (the PR template asks for it), not a timing check on CI, whose shared runners are too noisy to bound.

- **Use `reps=2000` (the `bench` default) before believing a small difference.** `bench` reports min-of-N over the golden grid; reproducibility is ±0.01 µs at 2000 reps and ±0.5 µs at `reps=40`, which is wide enough to invent or hide a few-percent effect.
- **Interleave the two binaries; do not run one build and then the other.** The ±0.01 µs is *within* one process. Between processes it is ~±0.1 µs: one sequential pair said a refactor was 4% *faster*, and interleaving the two binaries three times each said it was 1.7% slower, which is what held up.

  ```sh
  cp bench_solve bench_before      # ...switch branches, rebuild...
  for i in 1 2 3; do ./bench_solve | tail -1; ./bench_before | tail -1; done
  ```

- **Build both arms from one tree.** Code layout moves this benchmark ~2%, the same order as many effects worth measuring. A seven-point tolerance sweep built in one loop in a scratch tree read ~2× too cheap — +1.5% for a step where the controlled A/B said +3.4% — and a sweep is not a set of controlled A/Bs.
- **A loaded machine moves the absolute figures and the ratios.** A bare `.Call` moved 0.69 → 1.10 µs between two runs an hour apart on the same build, so record ratios to that `.Call`; and run nothing else heavy while measuring, since a busy machine stretches different arms by different amounts.
- **`bench_history.sh` refuses commits before #47.** They declare `Package: leaf`, install cleanly as `leaf`, and then `library(phylloptim)` silently falls back to the SITE build — so the row that comes out describes the current tree while looking like history. `bench_user_cost.R` takes `--expect-version` and hard-fails on a mismatch; any harness that installs an old commit needs the same check.
- **The new single-layer arms in `bench_solve` print `us/call`, not `us/solve`**, because `bench_history.sh` greps every occurrence of the latter into one TSV field.

### Why the cost rules are what they are

The rules themselves are in CLAUDE.md. The cases behind them:

- **Counting:** the regression that motivated the cost guard was a `.Call` appearing in a per-row path, a 4.5× regression in `set_drivers()` that reached a green CI. A count of `.Call`s has no variance; `test-cost.R` swaps a namespace binding for a counting one. ⚠️ **The two supply paths reach C++ through different symbols** — the multi-layer default through `root_network_from_carbon()`, the single one through `RootNetwork__ctor`. Counting the wrong one returns 0 forever and the test passes while measuring nothing.
- **The batch gradient:** `leaf_gradient_batch()`'s whole claim is that it crosses the boundary ONCE per call rather than 112 times per observation, and `test-cost.R` asserts one `gradient_batch_run` for any N and **zero** calls to any of the nine per-perturbation primitives the R route reaches through.
- **The user's entry point:** #66 nearly doubled `leaf_gradient()` (286× → 539× a `.Call`) while `set_drivers()` and `leaf_solve()`, the two entry points measured, were flat.
- **Documentation rots the same way:** five places in this package once quoted "204 µs" beside four different solve counts, because nobody re-derived the divisor. Hence one home per figure.
