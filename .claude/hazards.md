# Hazards — phylloptim

The eighteen hazards in this code, each of which has cost someone real numbers, with the measurements behind them. [CLAUDE.md](CLAUDE.md) lists them as one-liners and holds the standing rules; [measurements.md](measurements.md) holds the golden-file record and the benchmark method.

⚠️ **The numbers are stable.** Comments in `inst/include/`, `R/` and `tests/` cite "hazard N", so a hazard keeps its number for life. A resolved hazard keeps its slot and says what replaced it; a new one is appended as 18.

The last section, "One solver: the long form", is the evidence behind the rule of that name in CLAUDE.md.

## 1. There are TWO Weibull vulnerability curves

Stem (`stem_P50`, `stem_c`) drives `hydraulic_cost_TF`; root (`root_P50`, `root_c`) drives uptake. They used to be the unmarked `b`/`c` plus `root_b`/`root_c`, and the companion analysis used the root parameters for the stem cost and carried λ ∝ ψ^3.02 into a manuscript draft where it should have been ψ^0.64. Never leave an unmarked default for a parameter that exists in two versions.

⚠️ **Each curve has exactly TWO settable numbers, and everything else about it is derived.** `stem_b = P50/(ln 2)^(1/c)` and `psi_crit = b·ln(1/0.05)^(1/c)`, i.e. P95 of the same curve — and Sperry's `k_crit = 0.05·kmax` reads that same 0.05 from `Leaf::k_crit_fraction`. All four are readable and none is settable.

This is worth defending, because the failure it removes was expensive and quiet. A free `psi_crit` described a curve it was not derived from: the spline was pre-integrated over `[0, P99]`, `psi_crit` never entered that bound, and every solve evaluated the curve *at* `psi_crit` — so anyone fitting a measured vulnerability curve picked a plausible number and got a domain error naming only the interpolator. Deriving it removed the need for a consistency check and a domain check. Since #183 `psi_crit` bounds only the four curves whose own forms end there (`cost_domain_ends_at_psi_crit`); for the rest it is a labelled point on a curve that continues past it.

Do not reintroduce `stem_b` as settable. Solving `c = ln(ln 2)/ln(P50/b)` diverges as `P50 → b`, so a (b, c) interface has a singularity that a (P50, c) one does not.

⚠️ **Both curves are read from one kind of table, and nothing has a domain end but `G(inf)` (#183, #85).** `VulnerabilityTable` holds `g(u) = G(psi)/b` on `u = psi/b` as a quintic Hermite with exact `g`, `g'` and `g''` at every knot, so `G` is C2 across knots, and `G'` is read from the same polynomial: the derivative of the value the solve reads. It depends on `c` alone, so a P50 rebuilds nothing, and it is shared through a pointer, between copies and between the stem and the root where their `c` agree. Past the end of the curve (`x = (psi/b)^c = 40`, 15.44 MPa at the defaults) `G` is `G(inf)` to rounding. The root's conductivity is its closed form everywhere, with no clamp, and a layer it reaches zero in draws nothing. What stays `stem_curve_domain`: a negative potential and a flow past `G(inf)` beyond rounding, which no stem potential carries. The table's accuracy, measured against a 50-digit reference, is in `VulnerabilityTable`'s comment; the boost incomplete gamma it was seeded from is ~50 table reads per call, which is why it is not on the hot path (plant's TF24 stands ran 3.3–3.9× slower reading it).

⚠️ **Near `G(inf)` the curve resolves nothing.** `G` is read to ~2e-16 of `G(inf)`, so where the stem's remaining capacity `G(inf) - G` falls below ~1e-13 of it every gradient is rounding, and the collar solve refused there. `Leaf::stem_capacity_floor` (1e-11) shuts the leaf before that, at ~12.4 MPa at the defaults. A dry pin at the end of the curve holds the stem at `psi_end` (`seat_stem_at_dry_pin`), not where the inverse lands, for the same reason. Carrying the tail as `G(inf) - G` would resolve it and remove both.

⚠️ **The knot grid is indexed, and its last knot is asserted bit-exactly.** A grid built by accumulating `x += step` produced `resolution` *or* `resolution - 1` knots depending on `b` and `c` (#92). The table's last knot is `u_end` exactly and holds `G(inf)/b` exactly, so the end of the curve and the value past it are one number.

## 2. Every ψ is a positive magnitude in MPa (#25)

This slot used to be "signed versus magnitude water potentials". There is now one representation, so the hazard is removed rather than managed. Two conventions held together by comments cost a negative λ (`dE_from_soil_dpsi_collar` differentiated w.r.t. the *signed* potential and returned a negative where a conductance was wanted) and a dead clamp (plant #584, filed here as #24). Tagging the two conventions as types (#23, `Psi` / `AbsPsi`) was tried and closed unmerged: it described the two-convention model instead of removing one, and it could not cover the derivative, which is neither type.

`psi_soil_`, `psi_crit`, `root_psi_crit`, `stem_b`, `root_b`, `opt_psi_stem_`, `opt_root_psi_`, the collar variable inside the solve, the root-finding brackets and the vulnerability tables' arguments are the same kind of number. Where an equation needs one potential to oppose another the minus sign is **in the equation**: `E_i = (T_collar - T_soil_i - gravity_head*z_i)/r_R`.

Three things follow:

- **The invariant is asserted, not documented.** `set_physiology` rejects a negative `psi_soil`, the constructor rejects a non-positive `psi_crit` / `stem_b` / `root_b` / `root_psi_crit`, and the three entry points that take the soil state *as an argument* reject a signed vector.
- **`dE_from_soil_dpsi_collar` is a positive conductance by construction**, so the negative-λ class of bug is unrepresentable.
- **The member is `opt_root_psi_`, not `root_collar_psi_`.** Renamed on purpose: a flipped sign under the old name would let an old analysis silently read the wrong value, where a rename gives a binding error. plant's `opt_root_psi` aux and state agree in sign.

⚠️ **Do not call the convention "tension" in code.** It is exact only because there is no osmotic term and gravity is carried separately as `gravity_head`, so ψ here is a pure pressure potential. Add solutes and "tension" is wrong while "magnitude of ψ" stays right.

**The collar solve's dry end is the continuity root with the stem at `dry_end<K>()`**, taken from its feasible side, and the root's critical potential bounds nothing (#183). That removed the #24 / plant #584 clamp to `root_psi_crit` and, with it, the tie at the default traits (`psi_crit == root_psi_crit`) that put every shutdown on a kink in the stem and root traits.

## 3. Argmax smoothness is a hard constraint, not a preference

The demographic growth-rate gradient depends on the argmax varying *smoothly* with inputs, and any change to a solver must re-measure it. A safeguarded root-find on `dprofit == 0` meets it better than the golden-section search it replaced:

| | golden section | root-find |
|---|---|---|
| distinct argmax values over 11 trait steps | 6 / 11 | **11 / 11** |
| argmax second differences in a trait | 3.9e-04 (= the step size: noise) | **3.4e-07** |
| `\|dprofit\|` at the returned collar, interior rows | median 7.8e-04 | **median 5.6e-15** |
| cost per solve, interleaved at reps=2000 | 1 | **0.76** |

**A fixed iteration count is not smoothness.** Golden section terminates on bracket *width*, so it resolved the argmax only to `GSS_tol_abs` and the leftover offset wandered discontinuously as the comparison sequence flipped. That made the argmax piecewise constant at fine scales, so trait derivatives came back exactly zero — or, for traits in the hydraulic path, **smooth, plausible and sign-inverted** (`root_b`: −2.6e-03 against a true +2.6e-04).

⚠️ **The measurement above is in a TRAIT, and the hazard's concern is plant state.** Smoothness in height and light, feeding the demographic gradient, cannot be measured from here while plant #591 blocks end-to-end validation. Same mechanism, so the same direction is expected — but it is not the same measurement, and the trait result does not discharge the hazard. **Re-measure on the plant side when #591 clears.** `test_leaf.cpp` carries the trait measurement as a test.

`golden_section_max` still exists as the collar solver's fallback and in the single-layer optimisers; `optimize.hpp` says why. Do not read its existence as evidence that a comparison-based search is the safe default here.

⚠️ **The closed form fails this hazard, which is why it is not on a differentiated path.** Its own arithmetic is fine; what breaks is that its validity guard tests an **output**, so the set of drivers on which the method switches is a surface in driver space and the argmax JUMPS across it. **A mixture of two methods is discontinuous even where both are individually smooth**, and no tolerance fixes that. `set_model` resets the method to exact for this reason, so `gradient::route_seat` cannot pick it up. Measurements are in `closed_form.hpp`; the short version is that the mean |second difference| of `opt_psi_stem_` across a parameter sweep — which is the quantity a trait gradient differentiates, and which for a smooth argmax should shrink with the step — rises a hundredfold where the sweep crosses the guard.

## 4. The leaf is purely intensive

Every input is per unit leaf area or a dimensionless/intensive driver. Nothing scales with plant size — whole-plant allometry (`kmax(h)`, root carbon totals) is computed on the plant side and passed in already reduced. Keep it that way; it is a one-sentence contract.

⚠️ **The ROOT RESISTANCES arrive already reduced too (#33), and nothing here can check it.** `set_physiology` takes a `RootNetwork` rather than a root carbon profile, and a network built from *absolute* carbon instead of carbon per unit leaf area is just five vectors of positive numbers — the leaf cannot tell, and `E_up` comes back wrong by the leaf area with no error anywhere. This is the one place the boundary traded a guarantee for a convention.

## 5. Hot-path discipline — but measure it, because the intuition has been wrong

`find_root_collar_psi` runs ~10³ inner evaluations per solve, and plant calls it millions of times. λ and `g1_eff` are *accessors*, not stored state, for this reason, and don't add a `pow()` to `set_leaf_states_rates_from_psi_stem`.

What is *not* true is the blanket "prefer compile-time dispatch, a virtual call in the solve would cost". Measured (`make -C tests/cpp bench`), the code splits into two halves that answer oppositely:

- **The soil/supply path is already out-of-line.** `E_from_Soil_to_Root_Collar` is inlined into none of its 18 call sites at `-O2`, nor at `-inline-threshold=2000`; neither was the per-layer spline `eval` that the vulnerability table replaced in #183. So making it virtual costs **+1.1%**, `std::function` +0.6%, `std::variant` +2.6% — noise. Forcing more inlining is *worse*, by 4%.
- **The cost core is fully inlined.** `profit_psi_stem_TF`, `hydraulic_cost_TF`, `assim_colimited` and `transpiration_to_psi_stem` have no out-of-line symbol at all, so a change there is a different question and needs its own measurement.

Before arguing from inlining, check: `nm -C test_golden | grep <fn>`. No symbol means inlined; a symbol plus `bl` call sites in `objdump -d` means it is not. How to time an A/B is in [measurements.md](measurements.md#benchmark-method).

## 6. There are TWO supply paths, chosen by `Leaf::supply_kind_`

`MultiLayerRoots` (default) and `SinglePotential`. Both are held as members and selected by an enum, which measured **free** where `std::variant` cost +1.0% — a predictable branch in front of an already-out-of-line call disappears into it. Don't "tidy" this into a variant without re-measuring, and note the trap: a compile-time-known tag folds the branch away and reports a false zero.

**The contract is whatever `Leaf` dispatches on, and `grep -n 'switch (supply_kind_)' leaf_model.hpp` lists it.** Today that is nine operations a path must provide: `clear`, `set_soil_state`, `set_supply_resistances` (`set_root_network` on the multi-layer path), `begin_solve`, `uptake_at`, `duptake_dpsi`, the layer count, the soil-potential vector, and the collar-potential bound (`root_psi_crit` for roots, the stem's `psi_crit` for a path with no root limit). The single-soil-potential and single-layer queries are derived from those. ⚠️ **Every one of those switches has a `default:` arm that means "the single path"**, so a third enumerator compiles and silently runs the single path wherever an arm is missing. Add an explicit arm to every switch.

One rule rides with the contract: **`duptake_dpsi` returning NaN means "fall back to finite differences", so returning 0 or throwing silently degrades TF24f's acclimation gradient.** A path with no branch kinks simply never returns it.

`set_physiology`'s `RootNetwork` argument is part of the contract. Both paths take their soil-to-collar resistances from that one argument, so the calling convention does not depend on which path is in force. `MultiLayerRoots` reads one `r_R_H_min` and one `r_R_V_sum` per rooted layer; `SinglePotential` reads `r_R_V_sum[0]` as its series resistance and **refuses** a non-zero `r_R_H_min` rather than ignoring it. A third path must accept the argument and validate what it can use.

⚠️ **`gravity_head` is the one thing that is still configuration on the single path**, via `set_supply_single(gravity_head)`, and that is deliberate rather than unfinished: the multi-layer path derives a per-layer head from the depth profile it is given, and a path with no depth profile has nothing to derive one from — while a bare leaf wants zero rather than a geometric default. Do not "finish the job" by inventing a depth for it.

## 7. Moving a public member is a plant API break, because RcppR6 binds fields by name

plant's `inst/RcppR6_classes.yml` lists most of `Leaf`'s state as `access: field`, and the generator emits `obj_->psi_soil_` — a getter *and* a setter — straight into `src/RcppR6.cpp`. Relocating a member into a sub-object therefore breaks plant's generated glue, not just its own sources.

It is cheap to fix: for `access: field` the template pastes `name_cpp` verbatim after `->`, so a dotted path works and the R-side name (taken from the YAML key, not `name_cpp`) does not move:

```yaml
psi_soil_:  {type: "std::vector<double>", access: field, name_cpp: "roots_.psi_soil_"}
r_R_H_min:  {type: "std::vector<double>", access: field, name_cpp: "roots_.network_.r_R_H_min"}
```

The path can be nested arbitrarily, as the second line shows — so a member can be moved as deep as the design wants without the R API noticing.

Check `git grep 'access: field' inst/RcppR6_classes.yml` in plant before moving anything public, and land the two changes together — plant tracks this package's **master** via `Remotes:`, so a merge here is what breaks it.

## 8. Every path out of the solve must write its own rates

`Leaf` is a value member that plant reuses for every individual in a patch, so any output a branch declines to write silently becomes the *previous plant's* value. Three paths had this defect, and each fix is a pattern to keep:

- **Clear with `.assign`, never `.resize`.** `set_physiology` cleared `soil_consumption_` with `.resize`, whose fill reaches only newly-added elements, while the uptake loop writes only up to `max_soil_layer` — the deepest *rooted* layer, not the deepest soil layer. A shallow-rooted plant therefore inherited the deep layers' uptake, and plant bills that to the patch water balance (plant #577, fixed in plant by #585).
- **Every exit writes every output.** The `assim_max_ < 0` exit set `profit_` but left `transpiration_`, `stom_cond_CO2_` and `assim_colimited_` alone, which made `profit_ == assim_colimited_ - hydraulic_cost_TF()` false in that one branch (7.79 vs -1.47 when measured).
- **A method that reads seated state must check that it is seated.** `dprofit_droot_collar_psi` read soil potentials a solve had to have seated, so calling it without one **segfaulted**; and it missed the reversed-gradient case (`psi >= psi_stem`), where `psi_stem_to_ci` either throws or returns a gradient built on a negative conductance.

⚠️ **The golden file cannot see any of this, by design.** `test_golden` constructs a fresh `Leaf` per grid point (see its header), so `soil_consumption_` is always empty when `set_physiology` runs; the grid's minimum `assim_max_` is 3.71 so it never reaches that exit; and no golden column comes from the gradient. All three fixes are golden-bit-identical. **A bit-identical golden run is not evidence about a stale-state bug** — reach for a test that reuses one `Leaf`, and check whether the grid reaches the branch at all before reading its silence as a result.

## 9. No Rcpp in the leaf, and none in the include graph either

`util.hpp` throws `std::runtime_error` instead of `Rcpp::stop`, and uses a quiet NaN instead of `NA_REAL`. odelia's `ode_util.hpp` was the last R touchpoint; it was fixed upstream (traitecoevo/odelia#44), so the suite compiles against the real headers with no R anywhere (issue #11). Don't reintroduce Rcpp — and note the tests will not tell you if you do it indirectly by leaning on a header that pulls R in, because there is no shim to notice.

⚠️ **The invariant is directional, not "no Rcpp under `inst/include/`".** There is Rcpp inside `inst/include/` — in `phylloptim.h` and the three generated `RcppR6_*.hpp` files. The rule is: nothing reachable from `<phylloptim.hpp>` may include any of those four, and only `src/RcppR6.cpp` includes `<phylloptim.h>`.

What enforces it is that `cpp-tests.yml` builds the suite on runners with **no R installed**, so the moment a model header reaches for Rcpp, three jobs go red with a missing `RcppCommon.h`. That is also why the R job lives in a separate workflow: adding an `R CMD check` step to `cpp-tests.yml` would install R on those runners and quietly destroy the only test that can see this. The CMake install rules encode the same split from the other end — they ship `*.hpp` and exclude `phylloptim.h` and `RcppR6_*.hpp`, so a C++ or Python consumer is never handed a header they cannot compile.

## 10. Changing a trait is not one assignment, and the reason is invisible

The traits are public plain doubles, so `l.vcmax_25 = x` compiles — and leaves derived state describing the old value. Use `set_traits()`, which exists for exactly this and is the only route bound to R. It ends with `setup_clean_leaf()`, so `set_physiology` must genuinely be called again afterwards.

Three pieces of derived state depend on the traits:

- the **vulnerability tables**, one per `c` (`stem_c` builds the stem's, `root_c` the root's; a P50 builds nothing);
- the **solved operating point** (hazard 8);
- `vcmax_`, `jmax_` and `R_d_`, derived inside `set_physiology`'s **temperature cache**.

⚠️ **The cache key is every scalar `update_temperature_dependent_params()` reads**, `vcmax_25` and `jmax_25` among them, so re-driving after a bare trait write recomputes (#55). A key of `(leaf_temp_, atm_o2_kpa_)` alone once let a trait sweep written as "change the trait, call `set_physiology` again" take a cache *hit* and run the whole sweep at the first vcmax it ever saw, with every number plausible. **The guarantee lives in `photo_temp_key()`, not in the field** — a member added to the temperature block and left out of that key reopens this silently, and `photo_temp_key_size` is the only thing that will complain. `test_temperature_params_invalidate_cache` asserts it through the route a caller actually takes.

⚠️ **One deliberate exception: the thermal-damage layer's inputs are NOT in the key.** `set_physiology` never serves the cache while `use_thermal_damage_` is on, because `A_`, `I_` and `D_` are state a caller moves every step; only the gate itself is keyed, so switching it off at unchanged drivers recomputes. Remove that bypass and every input the block reads (the solve-reaching fourteen in the gradient enumeration, plus `k_i_Ea_`, `k_i_Tiso_` and `dTimpair_max_`) must join the key. `test_thermal_damage_bypasses_the_cache` (C++) and `test-thermal-damage.R` fail if either half goes.

The cache key does not clear the tables or the solved operating point, so `set_traits()` stays the way to change a trait. `test_set_traits_matches_a_fresh_leaf` asserts a re-traited leaf is **bit-identical** to a freshly constructed one, which catches all three at once — the two routes share no code, so anything not refreshed shows up as a difference. It solves at the defaults *first*, on purpose: on a cold object the cache would miss anyway and the trap would not fire.

**Measure before optimising a trait loop, in the layer you actually care about.** Two costs dominate, and which one bites depends on where you are:

- **From R**, constructing a `Leaf` is the largest single cost: most of a one-parameter gradient on a fresh leaf, and mostly R-side object construction over the active bindings rather than the tables. That, not the choice of gradient formula, is what dominates a finite-difference gradient in R. `set_traits` exists to avoid it: a re-trait plus a re-drive costs a few solves.
- **In C++**, `set_traits` is almost free unless it rebuilds a vulnerability table, and then costs many whole solves. Only `stem_c` and `root_c` rebuild: a table depends on `c` alone, so a P50 moves nothing in it, and `perturb_stem_P50` moves `stem_b` without even clearing the operating point.

The measured ratios for both are in the README's "Performance from R", the one place they are quoted; `tools/bench_user_cost.R` and `make -C tests/cpp bench_gradient` reproduce them. ⚠️ **Do not carry the R conclusion into C++**: with the boundary removed, the IFT composite beats differencing the solve several-fold on the traits that touch no table (`bench_gradient`'s summary lines), which an R measurement does not show.

## 11. A bracketing optimiser answers "where is the interior maximum", which is not "where is the maximum"

The single-layer profit is neither unimodal nor always interior: at a leaf hot enough that net assimilation is negative across the whole supply stream, the ProfitMax profit is highest at **full closure** and carries a second maximum out in the interior. At Tair 50 °C with the thermal cost on the profit runs −1.5314 at `psi_soil`, −1.5510 at 1.19 and −1.5459 at 1.88, so a search that steps in from the bounds reports an open stoma where the objective says shut.

So all nine stem routes maximise over a CLOSED interval through the one solver: endpoints included, an `n`-cell scan to pick the basin **where `basin_scan_cells()` asks for one** (`boundary_scan_n_` = 64 cells, and only with the energy balance on), and a refine that root-finds the first-order condition inside the winning cell. Against a 20001-point reference over 30 single-potential rows they are exact. The scan costs `n + 1` objective evaluations per solve where it runs, which is why it is gated; `bench_solve` prices every arm, and the current per-call figures are the solver table in `vignette("the-models")`.

⚠️ **The refine tolerance must scale with the CELL, or a finer grid is WORSE.** Brent terminates on bracket width, so a fixed tolerance comparable to a cell leaves the answer at the grid point — and adding grid points then makes it worse, because cells narrow while the tolerance does not. A shortfall count that *rises* as the grid gets finer is the signature.

⚠️ **Endpoints alone are a half-fix, and the row count hides it.** Adding the two bounds removes every boundary row while still leaving a worst-case shortfall of 3.05e-01 on a genuine interior second hump — worse than the 4.71e-01 it fixes. **Read the worst shortfall, not the number of rows short.** This is the standing objection to any derivative-only route: a root-find between the two ends is exactly the half-fix. What keeps both properties is multi-start — several sub-brackets, each yielding an exact stationary point, compared on profit — because a grid argmax is piecewise constant in the traits and so cannot feed a gradient (hazard 3).

The collar solve needs none of this: `maximise_profit_over_collar` tests the gradient's sign at each end and reports a pinned optimum explicitly, which is why 42 of 240 feasible golden rows are pinned and correct.

⚠️ **A STEM ROUTE IS NOT the collar solve with a restriction — it is a different model, and neither is a check on the other.** These optimise ψ_stem with upstream pinned at ψ_soil, ignoring the soil→collar path entirely, which is what their refusal message means by "non-**root-based**". They disagree with the collar solve on 20 of 30 driver rows, and the gap does **not** close as the root resistance goes to zero (1.67 MPa at r = 1e-2): as that resistance vanishes the collar loses its freedom, `[root_zero_E, root_crit]` collapses, and the collar solve loses precision and finally reports hydraulic shutdown (#176).

## 12. `n_theta` is a compile-time constant, and a short aggregate initialiser is legal C++

Several places fill `theta` with a literal list sized by `gradient::n_theta`, which is `par_PPFD`: everything before `PPFD` in the enumeration. `phylloptim::n_pars` is one more, because it counts `par_PPFD`, and must not size `theta`. Read both values off the enumeration in `leaf_model.hpp` rather than from here. An initialiser shorter than the array zero-fills the rest, so adding a parameter shifts `kmax` and `resistance` down a slot and drops `resistance` off the end **with no diagnostic at all**.

This has bitten twice. `test_leaf.cpp`'s batch `theta` put `kmax` into a new trait's slot and left `kmax` itself `0.0`; the only symptom was three observations failing to solve. `test_golden.cpp`'s `by_solver_threw[4]` was indexed by a `Solver` enum that grew to eight, writing past the end. Both now derive their size from one named constant with a `static_assert`. **Count the entries against `n_theta` whenever you touch either, and prefer a named count to a literal** — `bench_gradient.cpp` carried a literal `13` over an 11-element array through three trait-count changes, and it only ever crashed when the address layout happened to be unlucky.

⚠️ **Two of those initialisers are DELIBERATELY SHORT — they leave the model-owned prices zero-filled because the route they exercise reads neither — and that is exactly the case the zero-fill hides.** Both carry `static_assert(gradient::n_theta == N)` beside them, so the next appended parameter is a compile error rather than a shifted `kmax`. Three places assert that constant: `tests/cpp/test_leaf.cpp`, `tests/cpp/bench_gradient.cpp`, and the consumer program inside `.github/workflows/cpp-tests.yml`. Grep for `n_theta ==` and update all three together.

## 13. The trait vector is bound POSITIONALLY in four places

C++ `gradient::apply()`, R's `.gradient_setter`, the batch route's `theta` matrix, and R's derived copy of the enumeration. Two of those fail loudly on a length change; the batch's non-trait columns did not, when they were addressed as `length(par_names)` and `length(par_names) - 1L`, correct only while `kmax` and `resistance` were the final two in that order. They are addressed by NAME now, and `.gradient_non_traits` names the non-trait set once. Appending a parameter is safe; reordering the enumeration differentiates the wrong thing and returns plausible numbers.

⚠️ **`.gradient_setter` is a `do.call` over the derived trait vector, and must stay one.** When it spelled out fifteen subscripts into `l$set_traits`, adding a sixteenth trait broke it: 118 gradient-batch rows came back `status = "error"` where they had been `"interior"`, because the R reference threw and the batch did not, so the failure surfaced as an R-versus-C++ DISAGREEMENT rather than as a missing argument. A comment beside the old call predicted exactly this. **A comment predicting a failure is not a guard against it.** The sixteenth slot is now `g_min`, the residual conductance (#182), which went in through the `do.call` with no positional edit on the R side.

⚠️ **The MODEL-OWNED parameters are a second class: `CF77_lambda_`, `TF24_floor_lambda_o` and `LeastCost_beta`.** Each is a slot in the enumeration available for exactly one curve and refused, by name, on every other. `R/gradient.R`'s `.gradient_owned_pars` is the single table mapping parameter to owner — do not write a second one, and do not reason about "the CF77 slot" as though it were unique. The three sit together at `par_CF77_lambda`..`par_LeastCost_beta`, BEFORE the thermal-damage block: a new owned slot goes after `par_LeastCost_beta` and renumbers that block, which is safe because every consumer (plant included) reads the pack by name.

## 14. A derived `psi_crit` makes the P50 chain rule CONDITIONAL

Since the dry end is a quantile of the curve rather than a trait — `psi_crit` for the curves it bounds, the end of the curve (`b·40^(1/c)`) for the rest (#183) — moving `stem_P50` moves it, and `d/dP50` picks up a term through it. That term is **exactly zero at an interior optimum and order one at a dry pin**, so the clean identity `d/dP50 = (ln 2)^{-1/c} · d/db` holds only where the optimum is interior. Assert it keyed on `status`, never unconditionally: an unconditional assertion fails on pinned rows for a *correct* implementation.

## 15. A finite difference through a solve has a noise floor, and the floor moves with the solver

The error-versus-step curve is a V and only its floor is a derivative; a small step is the wrong instinct. The V's position depends on which solver produced the argmax: the collar solve and the stem routes both root-find the first-order condition and locate psi* to ~1e-15, so ONE step (1e-06) serves every route. When the stem routes refined on bracket WIDTH they resolved psi* only to ~1e-06 of the bracket, and a relative parameter step of 1e-06 differenced quantisation — measured 0.1855 against a true 0.0551, and the wrong SIGN on another curve. Re-sweep whenever a solver changes.

**Sweep the step and read the floor.** Two bugs were reported against a correct composite on the strength of one 1e-06 step. Each gradient route carries its own `fd_step` for this reason. `vignette("the-models")` Appendix C has the numbers.

## 16. Check that an assertion was actually exercised

A first-order-condition test whose λ was 100× off-scale put all nine rows on a bound, so the identity was never evaluated and the test passed. Print the count of rows that reached the assertion. λ's scale is set by the leaf — `marginal_cost_water()` runs 9e4–3e5 at the defaults.

## 17. The derivative surface fails silently in three ways, and none is a compile error

(a) A scalar member of `PhotoCapacity` or `SupplyDraw` that their `for_each_active` does not visit contributes no rows: its columns read exactly zero and every number stays finite. (b) The surface assumes `pars` and `supply` hold the SEATED traits and soil state in value, that `K` is the curve the solve used, and that the leaf has not moved since the solve; none of that is checked. (c) A recorded point must come back through `replay_operating_point`, which re-takes the polished wet bound; re-evaluating at the collar any other way drops the correction a `boundary-soil` row needs, and 258 of 643 non-negligible rows then moved by more than 1e-4 (worst 25x, sign flips on root traits). `test_a_replay_hands_over_the_rows_its_solve_did` pins it.

## 18. A shut leaf is the ordinary evaluator at its wet bound, and nothing writes its outputs by hand (#182)

The wet bound is where the supply delivers exactly the cuticle's residual flux: the collar where uptake equals `r*` on the collar route, `stem_wet_bound` on the stem routes. At `g_min = 0` it is the zero-uptake collar, found by the same call with the same bracket, which is what keeps the four golden files' open rows bit-identical. Three things follow that a later edit could break without a compile error:

- **Append the residual term, never interleave it.** Every expression the residual conductance enters (the conductance, the flux, the energy-balance term, the costs' flux) carries it as an added term that is an exact zero at `g_min = 0`. Re-associating any of those sums moves every golden row by an ULP, which hides a real change among 5000 noise ones.
- **The wet loop starts at `frac = 0` only when `g_min > 0`.** With a residual conductance the bound is a regular point (the cuticle still passes CO2), so a leaf that declines to open sits exactly on it. At `g_min = 0` there is no conductance there, and evaluating dprofit at `frac = 0` moved interior rows by 1e-14; the solve steps 1e-6 of the bracket inside, as before.
- **The realised residual flux is capped twice, by soft minima**: at `residual_cap_fraction` of what the stem can carry from the collar, and on the collar route at the same fraction of what the roots can deliver from this soil (`MultiLayerRoots::uptake_limit`, seated per solve in `residual_root_ceiling_`). Without the first the shut stem's potential goes to infinity and `G^-1` refuses; without the second, as a layer approaches the root integral's cap the wet bound has nowhere continuous to go and jumps ~10 MPa. The roots' ceiling is a member, so every stem route clears it (in `stem_wet_bound`) and the draw carries it onto the tape; an `implicit_value` whose residual reads it must list `draw.uptake_limit` among its inputs or the reverse sweep drops its rows, which `test_transpose` sees. The caps are a placeholder for the regime #183 leaves to a later step; do not replace either with a `min()`, which kinks every output in `psi_soil`.
- **The roots' ceiling is `U (U / (U + r))^2`, not `U`, and the shape is load-bearing.** Past the root integral's cap (#85) no layer conducts and the shut leaf sits at the wettest layer's gravity balance. With `U` the bound converged to `psi + z / (1 - phi)` instead (a 4.9e-3 MPa step); with `U^2 / (U + r)` it was continuous but its slope in `psi_soil` jumped (-2.3 against 1), because `U` itself vanishes linearly at the cap. The cube makes value and slope continuous; the second derivative still is not, and that goes with a closed-form root tail (#85, #184 section 2). `test_the_capped_regime_is_smooth_in_the_soil` asserts both at the cap.
- **The shut leaf's bound is bracketed by gravity balances** (`find_shut_wet_bound`): every layer draws nothing or less at the lowest balance among conducting layers, and at least nothing at the highest. One root-find to collar precision on that bracket replaced a walk, a coarse solve and a polish, and is what puts a dead soil's bound at the balance. The open leaf keeps `find_wet_bound`, which is what keeps its rows bit-identical at `g_min = 0`.
- **The `g_min == 0` test sits in the class, in front of each out-of-line residual function.** A call that only returns zero cost 2.5% of an open leaf's solve at the default. ⚠️ Measure before merging two residual calls into one pass: doing that in `dprofit` measured 3% SLOWER at `g_min > 0`, interleaved.

`test_smoothness_sweep` in `test_leaf.cpp` is the guard for the epic's smoothness claim: it second-differences profit, `A`, `E`, `psi_stem`, `lambda` and two derivative rows across soil, light, deficit and air temperature, with the energy balance off and on, one and three layers and the damage layer. A slope jump it finds is refined at a hundredth of the step before it fails, which separates a kink from curvature, and an interval whose change dwarfs its neighbours' is bisected, which separates a jump from a steep stretch. ⚠️ The refinement alone cannot see a jump: it is centred on a sweep point and a step between two points straddles it, which is how the 10 MPa jump above first passed. Kind changes are corners and are excluded by design. `test_a_shut_leaf_is_a_value_not_an_exception` holds #184's rule that a dry, dark or hot leaf comes back as a value, from wet soil to plant's 1000 MPa.

## One solver: the long form

The rule is in CLAUDE.md. This is what it rests on.

**The worked example is `TF24_floor`**, and it is the shape a new curve should have. What it touched: one enumerator, appended last (`n_cost_curves` is defined from the last enumerator by name, so a new curve appends and moves that definition, as `LeastCost` since has), one arm each in `with_curve()` and `curve_name()`, one arm in each exhaustive `if constexpr` table, a cost/profit/lambda trio beside its parents' — and nothing in any solver. Both routes, both supply topologies, the gradient and the closed form came for free. Its ONE parameter is a caller INPUT with no default (`TF24_floor_lambda_o`); the hydraulic half reads `TF24_cost_scale` and `TF24_beta2`. It is the one level parameter `align_curves()` cannot align, since the only value reaching TF24's conductance is zero.

⚠️ **Prefer sharing the parent's parameters to giving the new curve its own, when the new curve is a parent plus a term.** A curve that exists to be COMPARED with its parent should be the same curve at one parameter's zero, so "is TF24 missing a price of water?" is a one-restriction question instead of a three-restriction one. It also kept the trait vector from growing; hazard 13 has what a new trait cost.

⚠️ **The cheapest test of a new curve is a REDUCTION to an old one.** `TF24_floor` is `TF24` at `lambda_o = 0` and `CF77` at `TF24_cost_scale = 0`, and both are asserted BIT-FOR-BIT rather than to a tolerance — which works only because each term is the parent's own expression, so zeroing a parameter adds an exact zero to the other curve's exact value. It survives fused multiply-add too: with the addend an exact zero, `fma(x, y, 0)` is `round(x*y)` and `fma(0, y, q)` is `q`, so both reductions hold under `-ffp-contract=fast`. That test needs no fixture, no golden file and no tolerance to argue about, and it catches a wiring error anywhere between the enum and the reported outputs. **Reach for it before reaching for a new baseline.**

⚠️ **A new curve whose parameter is a field rather than a trait is unreachable from `leaf_solve()` until someone wires it.** `leaf_traits()` cannot carry a field and `leaf_solve()` builds its own `Leaf`. `leaf_solve()` and `leaf_batch()` take the owned parameters as arguments and `.seat_model()` applies them; `.solve_price_fields` is the one table saying which parameter belongs to which curve, and `.gradient_owned_pars` is its twin on the gradient side. `LeastCost_beta` is the third row in each (not a price: a dimensionless unit-cost ratio, which is why those tables are about *ownership* rather than prices). A fourth owned parameter needs a row in both, plus the write in `.gradient_setter` that pushes a perturbed theta back onto the leaf — the one place `LeastCost_beta` was first missed, which no finiteness check saw and only the R-versus-C++ batch comparison did. If it has a default (#156), it also needs a row in `.level_defaults`, the R copy of its in-class initialiser, which `.owned_value()` and the batch's theta matrix substitute for an NA on its own route.

⚠️ **`n` is set from measurement, and scanning by default is the expensive mistake.** Over a 1728-row sweep (8 air temperatures × 4 gate combinations × 6 soil potentials × 3 deficits × 3 light levels), the only objectives carrying two prominent interior basins are `TF24` (21 rows) and `JS22` (66), and **every one of those rows has the energy balance on**; the collar objective has none in 432 rows across the same range. So `basin_scan_cells()` returns 0 with the gate off. Re-measure before widening that — the valleys are 7.0e-01 and 3.5e+00 deep, so a missed basin is a large error, but a scan added "to be safe" is a 3× cost on the hot path for nothing. 64 cells is where the answer stops moving: against a 2001-point reference, 64 matches it on every row while 32 leaves 6 rows short by 3.9e-04 (the comment on `boundary_scan_n_`).

⚠️ **A scan is not a way to find a maximum you will later differentiate.** A scan argmax is piecewise constant in the parameters, so any derivative through it is a staircase. The stem routes refined on bracket width and left `|dJ/dpsi|` at a median of 4.4e-05, which made a finite-difference gradient return 0.1855 against a true 0.0551 and the wrong SIGN on another curve; and `|A|max` came from a 500-point scan, which is where ProfitMax's partial-versus-total gradient split came from. Both are root-finds now: median residual **1.2e-15**, one finite-difference step (1e-06) serves every route, and removing the scan made ProfitMax **10× faster**.

**The duplicates that were deleted.** `optimise_psi_stem_Sperry`, ProfitMax's grid walk, the collar solve's separate maximiser, and `util::maximise_over_closed_interval` (the call without the root-find, which ended up called only by its own test). Each was claimed to be a genuinely different case, and deleting each made the code both smaller and faster.

**Two things in the tree LOOK like a second solver and are not.** Both say so at the definition:

- **`Leaf::maximise_profit_over_collar`** is the one solver at `n = 0` PLUS a classification: every exit sets an `OperatingPointKind`, and 42 of 240 feasible golden rows are pinned and read it. The shared maximiser returns an argmax and nothing about how it got there, so folding them in means either an out-parameter no other caller wants, or losing the tag.
- **`dprofit_dpsi_stem` and `dprofit_at_collar_psi`** share ~90 lines of algebra. The stem form IS the collar form at `dpsistem_dpsi = 1`, `dgc_dpsi = 0` — mathematically, not TEXTUALLY, and the Identity arm's exact operation sequence is what keeps all four golden baselines bit-identical. Everything carrying no such constraint is already shared: `cost_deriv<K>`, `benefit_link_deriv<K>`, `temp_deriv`, `dprofit_energy_balance_term`.

So the shape to look for is not a function duplicating the solver's *algebra* — that can be load-bearing — but one duplicating its *structure* with a piece missing.
