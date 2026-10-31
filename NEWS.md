# phylloptim 0.24.0

## The boundary layer's share of the CO2 path (#122)

With the energy balance on, the leaf has a boundary layer, and the CO2 conductance treated the whole leaf-to-air path as stomatal (H2O:CO2 ratio `r_s`, 1.67). The boundary layer's share diffuses at `r_b = r_s^(2/3)` (1.407 at the default; the laminar-layer scaling of Sabot et al. 2022, eq 10), so the conductance is now the series `1/g_c = r_s/g_w - (r_s - r_b)/g_bw`, with `g_w` the total conductance to water (`gs_H2O`) and `g_bw` the boundary layer's (`gb_H2O`). The transpiration is set by the hydraulics as before; the boundary layer changes how much CO2 it buys. On the default traits over psi_soil 0.5-4, VPD 0.5-4 kPa, PPFD 500 and 1500, Tair 25 and 35 C, `A` rises by a median 0.3% at 5 m/s of wind and 0.8% at 0.5 m/s, at most 1.4%, and `gc` by at most 4.1%; no solve's profit falls. Off the energy balance nothing changes, bit for bit.

`$stom_cond_CO2_` is the split conductance; `gs_H2O` is still the total to water, recovered from it, so `gs_stomatal_H2O` means what it did. Where a candidate's total conductance passes what the boundary layer can carry (a deficit on `vpd_leaf_min`), the series has a pole; a soft floor (`bl_co2_split_floor`, width `bl_co2_split_floor_width`) keeps it finite and is exact to the last bit at every physical point. Goldens: `operating_points.tsv` and `primitives.tsv` are bit-identical; `psi_stem_optima.tsv` moves on 1,002 of its 2,880 energy-balance rows and none of its others, by at most 2.6% in profit at 25 C, with no status or kind changed. The gate-off solve's cost is unchanged.

# phylloptim 0.23.1

## Leftovers of the 0.23.0 rename (#149)

No number moves. `MultiLayerRoots::uptake_impl`'s local for a trial collar had been swept into the name `operating_collar`, which means the collar the solve placed; it is `collar_psi`. `vignette("derivative-surface")` now gives its sweep's total as the sum of its rows (1,920) and counts eight kinds, and comments that still called `SupplyDraw::E_up` "flux" use its name.

# phylloptim 0.23.0

## The derivative surface's names say what they return (#149)

C++ API only; no number moves, and every golden file is bit-identical. The R layer and every R-visible string (statuses, kind names) are unchanged.

- `n_pars` is now `n_leaf_pars`, so code that sized `theta` by the old name fails to compile rather than overrunning it. A `static_assert` ties it to `par_PPFD`.
- `Leaf::collar_at` is `operating_collar`: it takes no point, it returns the collar the solve placed. `Leaf::outputs_at` and `Leaf::LeafOutputs` are `billed_outputs_at` and `BilledOutputs`, apart from the unrelated `gradient::outputs_at`.
- `Leaf::marginal_at` is `dprofit_dcollar_at` (and `marginal_from_coords` is `dprofit_dcollar_from_coords`), since "marginal" means lambda everywhere else; `marginal_collar_slope` is `profit_curvature_at_collar`, which is what it returns.
- `Leaf::bound_at(bool wet, ...)` is two functions: `wet_bound_at<S>` (uptake meets the residual flux; it reads no cost, so it takes no curve) and `dry_end_collar_at<K, S>`.
- `SupplyDraw::flux` is `SupplyDraw::E_up`, the name the uptake has everywhere else. `SupplyAt` is `SupplyInputs`.
- `closed_form_integral` and `closed_form_curve` are `integral_with_closed_form_rows` and `curve_with_closed_form_rows`: the value is the table's and only the trait rows are the closed form's.
- `vignette("derivative-surface")` tabulates which R status each C++ kind is reported under.

Migration for plant: `tf24_strategy.h` renames `SupplyAt`, `LeafOutputs`, `marginal_collar_slope`, `collar_at`, `outputs_at` and `draw.flux`.

# phylloptim 0.22.2

## A series resistance below 1e-2 is refused (#176)
The collar route reads the single-potential flux as (collar - soil) / r, so as r falls the whole range of flux fits in a collar interval about E*r wide and a collar placed to `collar_root_tol` no longer resolves it. Measured against the stem route, the r -> 0 model, beyond the difference linear in r, the error was 2e-8 at r = 1e-2, 3.6e-5 at 1e-3 and 2e-3 at 1e-4, and below 1e-8 leaves reported shutdown. `set_physiology` (and so `set_drivers()`) now refuses a single-potential series resistance below `SinglePotential::resistance_min`, 1e-2 MPa s m2 (mol H2O)-1, five orders below the nominal 1e3. Nothing at or above it changes; the golden files are bit-identical.

# phylloptim 0.22.1

## A narrow collar bracket keeps its width when extended (#193)

On a curve not bounded at psi_crit, a collar interval narrower than `find_root_psi`'s bracket tolerance (`bracket_tol_at`, 1.5e-4 (1 + |psi|) MPa, so about 1e-3 at these potentials; reached where E times the series resistance is small) was extended to the end of the stem curve onto its own wet bound, and the leaf was reported shut while profit still rose across the interval: at psi_soil 4.8 and a series resistance of 1, profit -5.63 shut against -5.05 open. `extend_collar_bracket` now polishes a narrow extension as `prepare_collar_solve` polishes its first bracket. Measured on the single-potential path, rows open at r <= 10 over psi_soil 4.5-5.8 (profit up by as much as 1.14) and nothing moves at r >= 100; on the multi-layer path nothing moves at the default root carbon or 10x it, and six rows move by 1e-14 at 100x. The golden files are bit-identical.

# phylloptim 0.22.0

## Refuse rather than clamp; no corner left but the stomatal one (#183)

Every site that bent a value was counted on the goldens, on #183's acceptance sweep and on 60-day damage trajectories, and each is now a refusal at the entry point, a smooth form, or a kept corner with its reason recorded.

- **Drivers.** `set_physiology` (and so `set_drivers()`) refuses `atm_vpd` outside `(0, e_s(leaf_temp))`, a relative humidity outside (0, 1), on every path. The longwave term used to floor `e_a` at 1 Pa there; `vapour_pressure_min` is gone.
- **Leaf temperature.** The `[-40, 70]` C clamp is gone, with `leaf_temp_min`/`leaf_temp_max`: a heatwave leaf shut in still air solved *at* 70.00 C, a clamp artefact. `set_physiology` refuses drivers under which a reachable leaf (from `E = 0` to `k_max G(inf)`) would pass the Tetens pole or a non-positive Q10 of `R_d` (about 119 C at the defaults).
- **Sky emissivity** is held below 1 by a soft minimum of width `atmos_emissivity_softmin_width` (0.002) rather than a `min()`, which put a corner in every output's slope in Tair from about 39-48 C. Bit-identical wherever the Brutsaert fit is below 0.937.
- **Prescribed potentials are refused, not projected.** `evaluate_root_collar_psi()` (and `profit_at_collar_psi`) and `evaluate_psi_stem_at()` refuse a potential outside the feasible interval with the new infeasibility code `prescribed_potential`. `leaf_gradient()`'s and `leaf_gradient_batch()`'s status `"clamped"` is now `"refused"`: with `psi` outside the interval every output but `psi` is `NA`; within one step of an end the value is reported and the gradient is `NA`, as before. C++ `gradient::Status::Clamped` is `Status::Refused`.
- **Thermal damage.** The Arrhenius and softplus argument clamps to [-30, 30] are gone (bit-identical inside; the lower ends bound for a 619 kJ/mol species below about 25 C and for any leaf below 0 C); an overflowing rate is refused. A and I are read as themselves, `A <= -K_A`, `I < 0` and `D < 0` are refused at `set_physiology` with the infeasibility code `damage_state`, and `damage_trajectory()` refuses such initial states. The functional fraction has a soft floor of width `damage_surv_floor`, exact for `F > 3.4e-5`. Protection is signed, `k_i_baseline_ - k_i_`, where the defaults sat on a `max(0, .)` corner. Shifted optima and isokinetic temperatures outside [250, 350] K are refused.
- **Kept.** The stomatal corner `E_s >= 0` (the one `clamp_site` left, `CLAMP_STOMATAL_FLUX`), `vpd_leaf_min` (binds at candidate points only, and at a solved point only on SOX, JW26 and LeastCost dark rows pinned at psi_crit), SOX's and JW26's reduction factors held at zero past psi_crit (without it a dark leaf's negative A times a negative factor is a positive profit), Medlyn's beta in [0, 1] and CF77's soil price (published forms, now said so in the registry). `CLAMP_ROOT_VULN_*`, `CLAMP_LEAF_TEMPERATURE` and `CLAMP_COLLAR_POTENTIAL` are gone.

The acceptance sweep (`make -C tests/cpp smoothness`) now runs VPD to 6 kPa at 40 C, g_min zero and positive and the damage layer with the energy balance on, and prints how many triples reached the assertion: 73,070 with the shut-leaf fix below, worst refined slope jump 0.0087 of scale against 0.05. Goldens: `operating_points.tsv`'s 84 RH < 0 cells (VPD 4 kPa at 25 C) move to 3 kPa; `psi_stem_optima.tsv`'s energy-balance rows move by at most 3.6e-9 relative at 40 C and 7.4e-14 at 50 C (the emissivity); `primitives.tsv` is bit-identical. `bench_solve` is unchanged.

Migration for plant: catch `prescribed_potential` and `damage_state` as step rejections (TF24f's collar no longer comes back projected; TF24t's states are refused below zero); expect `"refused"` where `"clamped"` was; any driver set with a deficit at or past saturation is now an error.

## A shut leaf at g_min = 0 is continuous in the drivers (#183)

At the default `g_min = 0` a shut leaf has no conductance at its wet bound, so its reported state is taken 1e-6 of the bound's scale inside it. That step was taken from the polished bound where the polish moved it and from the bracket's coarser lower end where it did not, and a sweep crossing between the two saw A and E jump by 3e-6 of their range (at psi_soil 5.4407 with the energy balance on: A 2.6e-5 -> 2.7e-7 umol m-2 s-1). The step is now always taken from the polished bound. The acceptance sweep runs g_min = 0 with the energy balance on again (73,070 triples asserted). Goldens: 25 of 672 operating points and 274 of 5,760 collar rows of `psi_stem_optima.tsv` move, all shut leaves at g_min = 0, profit by at most 4.1e-6 relative (5e-5 absolute on the stem file's near-zero profits). With `g_min > 0` nothing moves.

Decided with it (#183 step 4): the stomatal corner E_s >= 0 stays, measured active at 55-61% of solved points on a broad grid (97.5% in the dark, 40% at PPFD >= 100, mostly hot or dry), because asymptotic closure needs a barrier term that moves every result, and solving in log E_s would cost a collar root-find per evaluation; the E_max < r regime is C2.

# phylloptim 0.21.0

## The stem may run past psi_crit, and G is a C2 table (#183, #85)

`psi_crit` no longer bounds the solve for TF24, TF24_floor, ProfitMax, JS22 and CMax, whose costs are defined past it: for these the dry end is the end of the stem curve, where it carries G(inf) to within e^-40. SOX, JW26, CF77 and LeastCost keep `psi_crit` as the end of their own domain (`cost_domain_ends_at_psi_crit<K>()`; the per-curve table is new and is the place to change that). A leaf pinned at the dry end is `BoundaryDry` (`"boundary-dry"`), which replaces `BoundaryCrit`; `BoundaryRootCrit` is gone, because nothing bounds the root either, and the enum after it is renumbered. A diverging curve's leaf shuts when the stem's remaining capacity from the wettest layer falls to `stem_capacity_floor` (1e-11 of G(inf), about 12.4 MPa at the defaults), a resolution floor rather than a modelled bound.

The vulnerability integral G is now a quintic Hermite table in `u = psi/b` with exact value, slope and curvature at every knot, so G is C2 and the slope the solve reads is the derivative of the value it reads. The cubic splines it replaces were C1 and erred by up to 8.7e-11 of G(inf); the table errs by 2e-15 at the default `c` and 2.7e-12 at `c = 1.2`. Because the table depends on `c` alone, a P50 change rebuilds nothing, tables are shared between copies of a `Leaf` and between stem and root when their `c` agree, and `ncontrol_default` is 800 knot intervals (the two tables hold about the doubles the 400-knot splines did). The root's conductivity is its closed form with no clamp (#85), and a fully embolised layer draws nothing. A closed-form G through Boost's incomplete gamma was measured and rejected: plant's TF24 stands ran 3.3-3.9x slower.

Cost: the collar solve is 1-4% slower on rows comparable in both builds (TF24_floor 6.8%), and plant's TF24 stands about 11-12%; a `stem_c` rebuild costs about 17 solves, a P50 change nothing. The goldens are regenerated; what moved, by cause, is in the PR. `operating_points.tsv` gains a 17 MPa soil (672 rows), where every curve shuts.

Migration for plant: `"boundary-crit"` and `"boundary-root-crit"` become `"boundary-dry"`; `supply_psi_crit`, `root_psi_crit_at`, `invert_stem_curve*`, `eval_stem_curve` and `closed_form::PastLastKnot` are removed; `MultiLayerRoots::setup_vulnerability` takes an optional table to share. `Control`'s default follows `ncontrol_default` and needs no change.

# phylloptim 0.20.0

## `collar_at` can hand its coordinates to `outputs_at` (#161)

At an interior point `collar_at` solves the stem potential and intercellular CO2 on its way to the marginal, and `outputs_at` solved them again for profit. Both now take an optional `HeldCoords<S>*`: `collar_at` fills it, `outputs_at` reads it, and refuses one placed at another collar or kind. Values and forward tangents are bit-identical and reverse rows agree to rounding (worst 1.0e-15 on plant's gradient references). plant's long-drought reverse sweep (108 introductions, `ode_tol` 1e-3) is about 4.6% faster once its TF24 strategy passes the struct. Migration: none; both arguments default to null.

## The cost curves' defaults are aligned, so they can be compared at defaults (#156)

Each curve with a level parameter now defaults to the value at which its stem optimum has TF24's conductance at one reference environment (single-potential supply, `psi_soil` 1 MPa, PPFD 900, VPD 2 kPa, 25 C, energy balance off). The new `align_curves()` solves those values, and the defaults are its answer at its own defaults: `JS22_gamma` 1.0 -> 0.4718, `CMax_a` 0.6 -> 0.6386, and `CF77_lambda_` and `LeastCost_beta`, which had no default and refused to solve without one, now 107751 (TF24's own emergent lambda there) and 51.60. Previously the eight curves gave conductances 1.4-fold apart at the same drivers, and that spread said where each curve's parameters sat rather than how it responds. TF24 is the anchor and does not move, so neither do plant's results. `TF24_floor_lambda_o` keeps no default: TF24_floor reaches TF24's conductance only at zero, where it is TF24. SOX, JW26 and ProfitMax have no level parameter and are left where the supply puts them.

The three golden baselines are bit-identical; `test_golden.cpp` now pins JS22 and CMax at their former defaults and clears the two new ones, since it tests the solver rather than the defaults. A value cleared to NA on the leaf is still refused by name. `leaf_gradient_batch()` no longer clears the batch leaf's `CF77_lambda_` and `LeastCost_beta` when run on a route that does not own them, and an NA in an owned column of a supplied `theta` takes the default, as it does in `leaf_gradient()`.

## Documentation (#150, #151, #152)

The developer guide is split into `.claude/CLAUDE.md` (the standing rules), `.claude/hazards.md` and `.claude/measurements.md`, with hazard numbers unchanged (0.19.0's hazard 18 included) and its stale figures and counts corrected against the code (#150). `?phylloptim` is now a package help page, and a pkgdown site (`_pkgdown.yml`) groups the reference by purpose and carries all five vignettes as articles; `docs.yml` builds it beside the Doxygen API and publishes both as one site, still gated on `PUBLISH_DOCS` (#151): pkgdown at the root, the C++ API under `cpp/`, and the resolved registry at `cpp/model-registry.yml`. The cost figures are re-measured at 0.16.2, recorded as a new row of `tools/cost-baseline.tsv`, quoted as ratios in one place each (the README's "Performance from R" for the R side and `set_traits()`, `vignette("the-models")` for the solver) and linked from everywhere else. Two of them changed direction: constructing a `Leaf` from R now costs about 97 solves rather than 45, and a `set_traits()` that rebuilds a vulnerability curve about 16 C++ solves rather than 3.5 (the 400-knot default postdates both old figures), and driving one leaf by hand is now about half the per-row cost of `leaf_solve()` rather than level with it (#152).

## Each golden file has its own cross-platform tolerance (#125)

`test_golden --cross-platform` compared `operating_points.tsv` and `psi_stem_optima.tsv` with one shared pair, set before #173 when the worst disagreements were ~100-1000x larger. Each file now has its own, at ~50x its worst on Linux CI: `operating_points.tsv` profit 1e-9 (measured 2.13e-11) and argmax 5e-6 (8.6e-8); `psi_stem_optima.tsv` profit 1e-10 (2.01e-12) and argmax 1e-7 (1.33e-9). That last class measured 6.63e-6 before 0.19.0, a shut leaf's near-zero `assim` cancelling; since #182 that `A` is an exact zero, and the class sits at the solver floor. Tests only; no number moves.

## `collar_at` refuses a curvature it cannot divide by (#146)

At an interior point `collar_at` divided every collar row by the curvature it was handed, and odelia refuses only an exactly zero or non-finite one, so a curvature of `-1e-300` returned `dcollar/dvcmax = +1.5e296` with nothing flagged. The interior arm now refuses, naming the value, a curvature that is not negative or whose magnitude is below `Leaf::interior_curvature_floor` (1e-3, plant's default `gradient_curvature_floor`; the smallest on `test_transpose`'s grid is 0.14). The double path never reads the curvature and refuses nothing. No number moves.

## The derivative surface's silent failures are loud (#148)

Four ways to get every number finite and some of them wrong now fail by name. `collar_at` and `outputs_at` refuse a curve `K` other than the seated `cost_curve_` (`find_root_collar_psi()`, which solves TF24 whatever is seated, is the one solve this cannot see). `supply_draw_at` refuses a leaf on the single-potential path, since it walks the multi-layer supply whatever the leaf solved on. `PhotoCapacity` and `SupplyDraw` carry a `static_assert` on their layout, so a member added without a matching `for_each_active` entry is a compile error rather than a column of exact zeros. `SupplyDraw::flux` is zero-initialised, where a shutdown draw used to return it indeterminate. And `leaf_control()`'s `vulnerability_curve_ncontrol` default is read from C++ `Leaf::ncontrol_default` rather than restated as a literal. No number moves.

## Tests for what the derivative surface claims (#147)

Six claims the derivative surface made with nothing checking them now have tests. The supplied trait rows are refereed in-tree against a central difference of the solve at every kind `collar_at` places (interior, the wet bound as `boundary-soil` and as `hydraulic-shutdown`, both dry bounds): 387 rows, worst 1.8e-4 relative on a row of size 3.7e-5, at a relative step of 1e-4 read off a step sweep. `test_transpose` now reaches boundary-crit and boundary-root-crit, asserts that every kind was reached, and checks that coordinates handed over give the rows solving them again would. The stem conductivity's closed-form P50 row is bounded against the table's own response, at 9.6e-6 of the largest row at 400 knots and 1.5e-4 at 100. The stem curve's inverse round-trips to within 0.91 of its conditioning `eps * w / G'`. The template `root_network_from_carbon` is bit-identical to the double overload. A forced `method = "ift"` at a pinned row is refused by name in the batch. `leaf_control()`'s knot count against the C++ default is covered by #148. Tests only; no number moves.

# phylloptim 0.19.0

## A residual (cuticular) conductance, and one shut leaf in place of the no-flow family (#182)

A leaf now loses water and exchanges CO2 through its cuticle after its stomata have shut. `g_min` is a new trait, the residual conductance to water vapour at 25 °C in mol m⁻² s⁻¹, default 0, set through `leaf_traits(g_min = )` and refused below zero. It carries CO2 at `rho g_min` with `rho` = 0.05, the upper end of Boyer (2015)'s 1/20–1/40 for the cuticle's CO2:H2O conductance ratio. It rises with leaf temperature along one smooth curve, a softplus blend of two slopes about a wax transition near 40 °C, exactly 1 at 25 °C; its four parameters are placeholders. This is step 1 of #183.

The leaf optimises its stomatal flux `E_s`; the xylem carries the total `E = E_s + r`, which the stem potential fixes through the supply relation as before, so the water budget closes (`E_up = E` on the collar route) and every cost that reads the flux (CF77, TF24_floor, LeastCost and the ProfitMax normalisers) charges the total. The first-order condition keeps its form: with the leaf temperature an explicit function of `E`, `g_min(T_leaf)` needs no linearisation, and its temperature slope joins the energy-balance term. `d E_s / d r` matches the manuscript's `-Theta''/(|B'| + Theta'')` to four figures. Medlyn's `g0` is untouched.

Where the soil cannot supply what the cuticle loses, the realised residual flux is capped at half of what the stem can carry from the collar and, on the collar route, at half of what the roots can deliver from the soil (`MultiLayerRoots::uptake_limit()`, each layer's uptake as the collar dries without bound). Both are soft minima of power 8 rather than `min()`, so no output has a kink in `psi_soil`. Each is inert while `r` is well below it (5e-4 at half the cap) and both are placeholders for the regime #183 leaves to a later step. Without the roots' cap the shut leaf's stem jumped ~10 MPa across 0.005 MPa of soil as a layer approached the root integral's cap. The roots' ceiling is `U (U / (U + r))^2` rather than `U` itself, so the realised flux falls faster than the supply as the last layer stops conducting, and the shut leaf reaches a dead soil's bound, the wettest layer's gravity balance, with its value and its slope in `psi_soil` continuous (stem 1.0000 and profit -0.0866 on both sides). Its second derivative still changes there, which goes with a closed-form root tail (#85).

**Cost.** At `g_min = 0` the collar solve is 1.3% slower than 0.17.0 (3.27 against 3.23 µs on `bench_solve`, interleaved): open leaves are level with 0.17.0, and shut leaves, which are now evaluated where 0.17.0 assigned them, take 0.42 µs against 0.16. The stem routes are level or faster. A residual conductance of 0.5 mmol adds 9%. Two things keep it there: the soft caps take their eighth root by three square roots rather than `pow()`, and below `q^8 = 2^-54` return the demand exactly; and `g_min(T)` is memoised on all five of its inputs.

**One shut leaf.** A leaf with no interval to optimise over sits at its wet bound, where the supply delivers exactly the residual flux, and every output there is the ordinary evaluator's: nothing is written by hand, so the shut leaf and the open leaf beside it are one function. At `g_min = 0` the wet bound is the zero-uptake collar and the stem sits at it, rather than being held at `psi_crit` or the root's critical potential, so it keeps following a drying soil and the kink #175 recorded in the stem and root traits at the default traits is gone. The cuticle's CO2 pathway gives the supply-demand balance a root for every leaf: a leaf in the dark or too hot to fix carbon is an ordinary optimum, a net source whose `ci` sits above `ca`. At `g_min = 0`, with no conductance at all, the shut leaf takes the `g_min -> 0` limit of that: the root of `A(ci) = 0` where one exists, otherwise `ci = +Inf` and the leaf nets `A` as `ci -> Inf` (`-R_d` in the dark). Everything that assigned the no-flow state is deleted: `min(0, A(ca))`, `ci = ca` and `ci_assigned()`, the shade-death kind and its decision at `Tleaf(0)`, and the special values `set_shutdown_state` wrote.

**What plant asked for (#184).** A dry, dark or hot leaf comes back as a value, never an exception: 3168 solves from wet soil to 1000 MPa, dark and in full light, at 10, 30 and 45 °C, with the energy balance off and on, on the collar route and every stem route, and none throws. Every root a reported value rests on is solved to 1e-13 or polished: the collar wet bound to collar precision, the stem wet bound to 1e-15. The surface's outputs name the residual flux beside profit and the per-layer draws, with rows matching its value, so a consumer reads the soil draw, the stem's total and the cuticle's share by name (plant#655). The vulnerability closed form, the curvature floor (#146) and removing `psi_crit` as a bound stay with later steps of #183.

**Smoothness.** `test_leaf` second-differences profit, `A`, `E`, `psi_stem`, `lambda` and the `vcmax_25` and `g_min` derivative rows across soil potential, light, deficit and air temperature, with the energy balance off and on, one and three layers, and the thermal-damage layer, at `g_min` = 0.5 mmol m⁻² s⁻¹. It finds no slope jump above 0.05 of each output's scale away from kind changes, and, bisecting any interval whose change dwarfs its neighbours', no jump above 1e-6 of each output's range. The same holds from `psi_soil` 4 to 12 MPa at a large residual conductance, through both caps. Profit is C1 across the light-compensation point; `A` keeps a kink there, which is the co-limitation's own.

**C++ API:** `set_traits()` takes `g_min` as its sixteenth argument, after `CMax_b`. `LeafOutputs` gains `residual_flux`, and `SupplyDraw` gains `uptake_limit` (visited by `for_each_active` where it is in force). `par_g_min` is 15 and every later `par_*` moves up one, so `gradient::n_theta` is 35 and `n_pars` 36 (after 0.18.0's 34 and 35). `OperatingPointKind::ShadeDeath` is removed, so the later kinds' codes drop by one (`Prescribed` 6, `SolverRefused` 7, `NonFiniteGradient` 8). `HydraulicShutdown` now carries rows, as `BoundarySoil` does. Removed: `ci_assigned()`, `ci_no_flow()`, `no_flow_assim()`, `no_flow_assim_at()`, `no_flow_profit()`, `set_shutdown_state()` and `root_limit_binds_first()`. New: `set_shut_state()`, `find_wet_bound()`, `stem_wet_bound()`, `solve_ci()`, `assim_saturated()`, `residual_flux()`, `g_min_at()` and `supply_ceiling()`, among others in the registry. `optimise_closed()` refuses `g_min != 0`. The uptake walk skips a layer past its root cap rather than dividing by a zero integral.

**R API:** `leaf_traits()` gains `g_min`; `operating_point()` and `leaf_solve()` gain a last column, `residual_flux`; `leaf_infeasible_codes()` gains `stem_wet_bound`. A finite difference in `g_min` at zero is centred one step up, so both of its points are admissible.

**What moves:** nothing at `g_min = 0` on an open leaf. `primitives.tsv` is bit-identical. `operating_points.tsv` moves the 96 shut rows of 576, all at `psi_soil = 6`: the stem goes from `psi_crit` (5.87 MPa) to the zero-uptake collar (6.0 to 6.34), profit by the hydraulic cost there, `ci` on 24 of them to the root of `A(ci) = 0`, and `e_up` and `uptake` on the 64 multi-layer rows to the polished bound's residue (~1e-20, of either sign by platform, so the tolerant comparison now gives the water fluxes an absolute floor of 1e-15). `assim` stays exactly 0: with nothing conducting, `A` at the `ci` root is zero by definition. `psi_stem_optima.tsv` moves 3060 of 5760 rows: 342 collar rows go from `shade-death` to `boundary-soil` (330) or `interior` (12), with the argmax stepped 1e-6 of the bracket inside the bound; 1504 dark or net-source rows move only `ci`, to `+Inf` or a root above `ca`; 1014 hot rows move `A` without moving the argmax, now `A` at their own `ci` rather than at `ca`, among them the reused-leaf pass's ProfitMax degenerate rows, which now write `ci` and `A` instead of leaving the poison (128 fewer unwritten cells); 120 shut collar rows at `psi_soil = 6` move their stem by up to 0.32 MPa as it follows the soil past `psi_crit`, 96 of them with profit lower by the cost of that stem; and 80 stem-route rows at 40 and 50 °C move their optimum, every one with profit higher or equal. `gradient_golden.tsv` moves 6 cells, all on the shut-down case: its `stem_P50` and `stem_c` stem and collar rows go to 0, and its profit rows become the cost's slope at the new stem (1.15 and -0.83). `test-golden.R`'s shut-down row moves in `psi_stem`, `opt_root_psi`, `profit`, `e_up` and `uptake`.

# phylloptim 0.18.0

## Thermal impairment takes Bison & Michaletz's kinetics; every damage rate is one Arrhenius form

The impairment rate was `k_i_` times a saturating switch K/(1+K), so it levelled off at `k_i_` (0.5 per day) above about 40 °C: a leaf at 50 °C took over a day to half-impair, where measured PSII inactivation takes minutes. It is now the first-order heat inactivation of Bison & Michaletz (2026, Nature Ecology & Evolution, doi:10.1038/s41559-026-03196-x), written about the isokinetic temperature their Meyer–Neldel relation makes common to every species:

```
k_imp(T) = k_i_   * exp[(k_i_Ea_/R)(1/T_iso - 1/T)],  T_iso = k_i_Tiso_ + topt_offset_ + dTimpair_max_ * A/(A + K_A_)
k_rec(T) = k_rec_ * exp[(k_rec_Ea_/R)(1/T_25 - 1/T)]
k_mat(T) = k_mat_ * exp[(k_mat_Ea_/R)(1/T_25 - 1/T)]
```

The defaults are their 177-species fit: `k_i_` = 276 per day (3.2e-3 s⁻¹, from the corrected A values, Zenodo 10.5281/zenodo.23752859 v3), `k_i_Tiso_` = 50.5 °C and `k_i_Ea_` = 316 kJ mol⁻¹, so impairment rises tenfold about every 6 K: 0.012 per day at 25 °C, 5.4 at 40 °C, 230 at 50 °C. First-order damage at a constant temperature then halves a leaf in a minute at 54 °C, an hour at 43 °C and a day at 35 °C, as their figures show. Repair no longer shuts off when hot. Impairment outruns the default repair above about 32 °C, so the logistic shut-off had nothing left to do.

**Breaking.** Removed: `impair_H_`, `impair_Tmid_`, `k_i_Tref_`, `m_rep_`, `t_rep_cut_` and `$logistic_()`. Added: `k_i_Tiso_`, `k_rec_Ea_`, `k_mat_Ea_`, all at their defaults inert except `k_i_Tiso_`. `k_i_` and `k_i_baseline_` move from 0.5 to 276. A negative `k_i_Ea_` is refused. In the gradient enumeration `m_rep_` and `t_rep_cut_` give way to `k_rec_Ea_`, so `n_theta` is 34 and `n_pars` 35, and a C++ caller filling `theta` by hand must recount. `damage_par_names()` has fourteen dynamics parameters.

Off by default, and off it reads none of its inputs: every golden file and `gradient_golden.tsv` is bit-identical.

# phylloptim 0.17.0

## A leaf's profit is continuous where it shuts (#175)

A leaf's profit jumped wherever it stopped moving water or gas: down by `R_d` at hydraulic shutdown, up by `R_d` across the light-compensation point, and on entering the narrow-interval `Determined` state, where the stem drew up to 194x what the roots supplied. Hot, bright leaves under the energy balance were refused outright.

A leaf moving no gas now nets `min(0, A(ci = ca))` at its own temperature: zero in light, which is what a transpiring leaf tends to as `E -> 0`, and `A(ca) < 0` in the dark or when too hot to fix carbon. One function, `no_flow_assim_at<T>()`, supplies it to the shutdown and shade-death exits, the zero-transpiration branch and `outputs_at`, so the tape differentiates the value the solve reports. When `psi_stem_to_ci` finds no root in `[gamma*, ca]` an open leaf also sits at `ci = ca`; only a genuine absence of a root reaches that fallback, and a solver failure stops instead. The derivative there is `A(ca)`'s own slope in temperature, which is exact, so `dprofit_dpsi_stem` and `dprofit_at_collar_psi` now report it as feasible and the optimisers refine on it rather than falling back to a bracketing search. Shade death is decided at the shut leaf's temperature, `leaf_temp_from_E(0)`. A narrow collar interval has both ends polished before anything is decided on it, and one that has closed is hydraulic shutdown, so every open leaf takes its stem from continuity. The shut stem holds at the root's critical potential when the drying soil reaches it before the stem's, and at `psi_crit` otherwise.

⚠️ **At the default traits the stem's and the root's critical potentials are equal**, so the held stem has a kink in the stem and root traits exactly there. The rows take the stem's side; a central difference across it records the average of the two. `gradient_golden.tsv`'s shut-down case now sets `root_P50 = 4` to stay off it.

**C++ API:** `OperatingPointKind::Determined` and the kind name `determined` are removed, so the later kinds' integer codes drop by one. `ci_at_compensation_point()` is `ci_assigned()` and `respiration_temp_deriv()` is `temp_deriv()`.

**Not changed:** the Medlyn route, whose closed state still nets `-R_d`, so in light it still steps as it shuts. It is the only route with a residual conductance `g0`, and its output is byte-identical to 0.16.3 at `g0 = 0` and `0.022`.

**What moves:** `operating_points.tsv`, the 96 no-flow rows of 576 (`ci`, `assim`, `profit`), and none of the 480 open rows. `psi_stem_optima.tsv`, 2994 of 5760 rows: 1756 dark rows move `ci` to `ca` with `A` unchanged to 4e-16; 190 shut rows in light go from `A = -R_d` to 0; 904 too-hot rows take `A(ca)`; 106 energy-balance rows on the stem routes move their optimum, by up to 2.6e-6 MPa in the dark, where the refine is now exact (no row loses profit), and by up to 1.85 MPa, with profit up to 1.49 higher, for hot, bright leaves at 40 °C, where the stem optimisers now climb a continuous objective; 30 hot, bright energy-balance rows go from `solver-refused` to `shade-death`; and 8 `TF24_floor` stem-route optima now choose closure, which earns 0 rather than `-R_d`. `gradient_golden.tsv`, 3 cells, all on the shut-down case: `R_d_25`'s `A` and `profit` rows go to zero, and `stem_c`'s `profit` row picks up `-1.7e-10` of rounding. `primitives.tsv` is bit-identical. `test-golden.R`'s shut-down row moves in `ci`, `assim` and `profit`.

# phylloptim 0.16.4

## A reversed gradient is refused, not integrated (#31)

`transpiration(psi_stem, psi_upstream)` now stops when `psi_stem < psi_upstream`. A stem wetter than its upstream end would move water from the leaf to the soil, which the model does not represent, and the integral there was negative: `transpiration()` and `stom_cond_CO2()` returned negative values and `psi_stem_to_ci()` solved on a negative conductance, all as plausible finite numbers. All three now stop with an error naming both potentials. `psi_stem == psi_upstream` is the shut-down state and stays legal, with exactly zero flux.

The operating-point evaluator keeps its rule that `psi_upstream >= psi_stem` is no flow, so `profit_psi_stem_*()` below the boundary is the shut-down state, built on zero flux rather than a negative one. That rule has to stay: a prescribed collar at or clamped to the wet bound (plant's TF24f path) inverts to a stem up to about 1e-9 MPa wetter than the collar. The costs and lambdas that re-read the flux after the evaluator (CF77, TF24_floor, LeastCost, and ProfitMax's thermal term) read it through the new `operating_point_flux()`, which applies the same rule. They used to take a tiny negative flux there; with the guard alone they would have stopped.

**No solve asks `transpiration()` for a reversed pair.** A counter on the branch read zero over all three C++ golden files, the C++ unit suites, the R test suite and every vignette. All three golden files are bit-identical, and `bench_solve` is within 1% on every arm. The refusal's message is built in an out-of-line cold function: inlined into `transpiration()` it cost up to 3.8%.

# phylloptim 0.16.3

## The model registry: equations, parameters and code links from one YAML (#166)

`inst/model-registry.yml` holds the leaf as one record: 114 equations, each with its LaTeX, the C++ symbol that evaluates it, the parameters it reads and the bound fields it writes; and 132 parameters, each with symbol, units, class, meaning, source and a placeholder flag. Fourteen sections run from the drivers to the thermal-damage layer. Documents are generated from it rather than written beside it.

It stores no line number and no default value. `registry_resolve()` reads each symbol's file, line, body span and API anchor from Doxygen's XML, which the `Doxyfile` now emits with the references relation. `registry_values()` reads each value off a live `Leaf`, off `set_drivers()`'s defaults, or off a constexpr's initializer in that XML. `model_registry()` returns the registry as data frames.

CI checks it two ways. `tests/testthat/test-registry.R` asserts that every symbol and field the registry names exists and that every bound field of the `Leaf` is accounted for. `docs.yml` resolves every symbol in the XML through `tools/registry_symbols.R`, checks that what each equation reads is what its function references and that every literal is written in its function, and writes the resolved registry into the site as `model-registry.yml`.

The registry records the model and does not test it: there is no second implementation of the equations. Two tests cover what one would have caught and the suites did not: the closed form's CF77 arm is checked as the identity its USO collapse defines, and `lambda_molar` against its units. No numerics change; the golden baselines are bit-identical. The comment explaining Medlyn's closed state now sits on `set_medlyn_shutdown_if_closed()`, not on `medlyn_beta()`.

# phylloptim 0.16.2

## The Medlyn route reports a CO2-basis conductance (#129)

`medlyn_model_gs` returns the Medlyn (2011) USO expression, which is a **water vapour** conductance -- that is what the 1.6 in it is for, and what published `g1` values are fitted against. Both `solve_medlyn_ci_*` stored it in `stom_cond_CO2_` without converting, so that field came back **1.6x high** on this route. Checked against the model's own outputs, `A = gs*(ca - ci)` was violated by 1.5984-1.6004 across VPD 1-4 and two soil-moisture levels, where the optimality routes satisfy it to 1.00000.

`ci_`, `assim_colimited_` and everything downstream of them are **unchanged**: the solver residual carried a compensating factor, so scaling both sides leaves the root where it was. Only the stored conductance moves. All three golden files are bit-identical and the R suite is unchanged at 1371 passing.

The divisor is 1.6, not `H2O_CO2_stom_diff_ratio_` (1.67): the residual converts with 1.6, so 1.67 would report a conductance inconsistent with the `ci` returned -- a systematic -4.2% violation instead of a 60% one. It is a named constant now, `Leaf::medlyn_H2O_CO2_ratio`, used in all three places that have to agree.

Three smaller defects in the same function, all found while measuring the first:

- **The zero-deficit guard read the wrong variable.** It tested the `atm_vpd` input field while the expression below read the driven `atm_vpd_`, so a leaf driven at zero deficit took the else-branch and returned `Inf` from `sqrt(0)` -- the input is still at its 2.0 default and so never looks like zero.
- **`solve_medlyn_ci_analytical` ignored the soil-moisture factor.** `beta` multiplies `g1` in the coupled route but not in the closed form, so the two entry points disagreed about whether the model has a soil-moisture response at all: chi frozen at 0.645 against 0.681 to 0.446 over beta in [0.2, 1]. At `beta == 1` the value is unchanged, so only a soil-moisture sweep moves.
- **The soil-moisture factor is clamped to [0, 1]**, on both routes, through one function, `medlyn_beta()`. Unclamped, a soil wetter than field capacity scaled `g1` above its fitted value (theta 0.6 against field capacity 0.5 gave beta 1.33) and one drier than the wilting point gave a negative `g1` term. Inside [wilting point, field capacity] every value is bit-identical.
- **Neither route can report a negative conductance now.** With `g0 = 0` near the wilting point the USO expression goes negative, and the numerical route's "closest approach" exit returned it silently (gs = -3.95e-06). It collapses to the class's closed state instead -- no water moving, ci at the compensation point, still respiring.

⚠️ **The two entry points are two forms of ONE model**, which is now asserted rather than assumed: the coupled solve converges on the closed form as `g0 -> 0` (chi agreeing to 1e-6, against a gap of 0.036 at the default `g0`). The remaining gap IS the `g0` term the closed form drops. A consequence worth knowing when reading output: the analytical route's `gs` and `ci` are mutually consistent only in that limit, running 1.084 to 1.316 on the diffusion identity at the default `g0`. That is the approximation, not a units error, and it is pinned so the two cannot be confused.

This route had **no tests and no golden coverage**, which is how a 1.6x error survived in shipped code; `tests/testthat/test-medlyn.R` is the first, and the regression test is the substantive part of this change rather than the arithmetic.

# phylloptim 0.16.1

## The ci root-find is solved to 1e-13

`psi_stem_to_ci`'s fixed tolerance goes from 1e-10 to 1e-13 (#172). Where the root-finder stops inside its tolerance band depends on the state, so at 1e-10 every reported output was noisy along a smooth path through the drivers: up to 3.7e-11 relative along plant height, where plant takes finite differences of the leaf at a step of 1e-7 and so read an error of 1e-4. The IFT trait gradient inherited the same noise through its own central differences of `dprofit`, and on `gradient_golden.tsv`'s one-layer interior case it was wrong by up to 2.3e-3 relative; it now matches a Richardson-extrapolated difference of the solve to 1e-8. A solve and a gradient each cost about 5% more.

**What moves:** the three C++ golden files, at most 3.1e-8 absolute in `profit` and 2.1e-9 Pa in `ci`, with no point changing status, kind or branch; in `primitives.tsv` only the iterative tier moves. `gradient_golden.tsv` moves on its one-layer interior case only, by up to 2.3e-3, toward the reference. `ci_abs_tol` is unchanged and still reaches only the Medlyn route.

# phylloptim 0.16.0

## odelia's reverse-mode names, at 0.10.0

`LinkingTo: odelia (>= 0.10.0)`, the release that renames its reverse-mode API (traitecoevo/odelia#66). The leaf records onto an `odelia::ode::adjoint_scalar` (was `active_scalar`) and puts its supplied rows on the consumer's tape with `odelia::supply_derivatives` (was `record_with_derivatives`). A value and its slope in the collar potential is `odelia::with_slope<T>` (was `value_with_slope<T>`), named directly: `Leaf::pair<T>`, which only renamed it, is gone. A direction is set with `odelia::ode::set_direction` (was `seed_direction`), and "seed" in the derivative sense is gone from the comments: a trait carries a direction, and a sweep starts from an adjoint. A **minor** bump, because the floor moves and `Leaf::pair` is gone. No number moves.

# phylloptim 0.15.4

Three defects the model registry (#166) surfaced, each found by checking an equation's base-R replica against the package.

## The collar route reports the seated curve's lambda_emergent

`lambda_emergent` is the marginal cost of water the operating point implies, `dC/dE`, and it is what puts every cost curve on one axis. The stem routes have always reported each curve's own. The collar route reported `marginal_cost_water_multilayer()` after every solve and every prescribed evaluation: TF24's lambda with the series correction for the soil-to-collar path, whatever curve was seated. CF77 at a price of 1e5, whose lambda is that price by construction, read 1.05e6.

It is now `lambda_collar_for<K>()`, the same chain rule for every curve. The collar is the free variable, so uptake gives `dE/dpsi_r = S`, and the part of a cost that reads `psi_stem` alone takes the series factor `1 + kmax f(psi_r)/S`. The part that reads the potentials only through E (CF77's price, LeastCost's `log W`, the floor's `lambda_o E`, ProfitMax's thermal term) carries over unchanged. JS22 and CMax, which read the collar explicitly, add `dC/dpsi_r` over S. `test_collar_lambda_emergent_is_dA_dE` checks every curve, on one- and three-layer supplies, against `dA/dE` taken along the collar path at an interior optimum, where the first-order condition makes the two equal for every link. 26 interior optima agree to 1.8e-6. On the old reporting 23 of them fail, CF77 by up to 146%.

**What moves:** `lambda_emergent` (and the last element of `operating_point_values()`) on the collar route, for every curve but TF24. TF24's row is the function it always called, so its value is bit-identical, as are all four golden files. Nothing in the solve reads `lambda_emergent`.

## Four bound fields nothing wrote are no longer bound

`count`, `psi_stem` (the member, not `opt_psi_stem_`), `opt_ci_` and `f_r` were bound as R fields and written by nothing in the C++. `count` and `psi_stem` were never initialised, so a read returned whatever the memory held (about 3e-314 measured); `opt_ci_` was always `NA` and `f_r` always empty. They are unbound here. The members stay, because plant's own bindings name them (moving a bound member breaks plant's generated glue); `count` and `psi_stem` are now initialised to `NA`, so a read through plant's bindings is defined.

## The slope of the saturation curve is the exact derivative

`saturation_vapour_pressure_slope()` used FAO-56's quoted numerator 4098, which rounds `17.27 * 237.3 = 4098.171`, so the slope was 4.2e-5 low. The slope reaches the energy-balance first-order condition through `dgc/dT = -gc * Delta(Tleaf) / D_leaf`, and there it made the analytic `dprofit/dpsi_stem` miss a central difference of the profit by 2e-6 to 2e-5 relative at every step from 1e-2 to 1e-6 MPa: a plateau rather than the difference's V. The numerator is now the product of the two Tetens coefficients, and the same comparison falls to the difference's own floor (median 1.5e-8 over 34 points at a step of 3e-4 MPa). `test_leaf.cpp` now makes it with the energy balance on, where only the gate-off case was checked before; on the old constant it fails.

**What moves:** the energy-balance rows of `psi_stem_optima.tsv` and nothing else. 244 of its 5760 rows change, all with the energy balance on, none in status or kind: the argmax by at most 1.2e-5 relative, profit by at most 5.6e-10. `operating_points.tsv`, `primitives.tsv` and `gradient_golden.tsv` are bit-identical. With the gate off the slope is never read.

# phylloptim 0.15.3

## The root curve's derivative answers past its dry end

A soil layer drier than the root curve's last knot (6.89 MPa at the defaults) no longer refuses the derivative path with `stem_curve_domain`. The rows are those of what the table returns there: the last knot's tangent up to G(inf), whose trait rows are G's at the knot, and G(inf) beyond it, with `dG(inf)/db = Gamma(1 + 1/c)` and a digamma row in `c`. The conductivity read past the knot has no row, since the curve there is the 1% floor for every trait. The stem table instead reads the closed form past its last knot (0.14.1), so its rows there are the closed form's; `closed_form_integral()` and `closed_form_curve()` now take a `PastLastKnot` argument naming which tail the table has, and each call site states it. Values are untouched: the golden grid is byte-identical.

## A determined collar answers the derivative path

Where the feasible interval collapses and the collar is held at its midpoint, `collar_at()` no longer refuses with "no operating point to place at one that is determined". The collar's rows are the mean of its two bounds' rows: the wet bound's, and the dry bound's from whichever limit closed it, the continuity root or the root's own critical potential. The stem there sits at its own optimum with the collar held, not at continuity with the draw, so `collar_coords_at()` places it at `psi_crit` with no response to the collar; placed by continuity, profit's row in `root_P50` read 0.053 where a difference of the solve gives -0.83. A stem optimum inside its interval is refused by name. Each bound's row is taken at that bound's own root, not at the midpoint where the draw was taken, up to 5e-4 MPa away: the draw is carried to the root along its slope, and the residual's slope is the supply's at the root. Taken at the midpoint, the collar's row along `psi_soil` read 0.50004 where its bounds give 0.5 exactly; now every collar row matches a difference of the solve to 2e-8.

Both bounds are now polished to collar precision before the midpoint is taken, as a pinned bound already is. Left at the bracket's 1e-4 tolerance, the continuity bound sat up to 5.6e-7 MPa from its root, and a difference of the solve along `root_P50` read 0.014 where the row, and a difference of the polished solve, give -1.6e-6. So a determined collar moves by up to 2.8e-7 MPa, and a stem conductance of 1e-30, refused before as `stem_curve_domain`, is now solved as a determined point. No other point moves: the golden grid has no determined point and is byte-identical.

## `uniroot()` refuses a bisection that spends its budget

`uniroot()` checked `it > max_iterations`, which cannot fire: boost's `bisect` reports the iterations it used and never exceeds its budget. It now checks `>=`, as `uniroot_smooth()` does. Its one caller, the Medlyn ci solve, needs a few tens of halvings to reach `ci_abs_tol` against a default budget of 1000, so no solve that finished before refuses now.

# phylloptim 0.15.2

## LeastCost's cost moves with leaf temperature in its first-order condition (#167 follow-up)

⚠️ **LeastCost results move with the energy balance on; with it off they are bit-identical**, and so is every other curve. LeastCost's capacity term reads `vcmax_` at the leaf temperature, so with the gate on its cost `log W` moves with `psi` through `Tleaf(E(psi))`. Until 0.15.1 the stem routes refused a non-identity link's derivative with the gate on, which hid this. 0.15.1 lifted that refusal and gave the benefit its `h'(A)`, but left out this cost slope, so LeastCost root-found the wrong first-order condition on both routes. Its derivative was 1 to 4% off a finite difference of its own objective, at every step size.

`dprofit_energy_balance_term` now subtracts `beta*umol_to_mol*dVcmax/dT*dT/dpsi / W`, beside ProfitMax's thermal slope. It is cost-side, so it takes no `h'`. The derivative now matches a difference to 4e-6, the same floor as SOX. `lambda_LeastCost()` carries the same slope in `dTheta/dE`, through the new `vcmax_temp_deriv()`, as `lambda_ProfitMax()` carries its thermal term. At an interior optimum it equals `dA/dE` along the stream to 4e-6.

**What moves.** The psi_stem golden file changes on LeastCost's 192 energy-balance rows, in `lambda_emergent` only, with no sign change. Its optima there do not move, because at `beta = 146` every one is pinned at a bound, which is also why no golden file saw the wrong derivative. On the collar route LeastCost pins at most drivers too; at `beta = 2` with three soil layers at 25 °C the old collar sat 0.02 MPa from a dense scan's argmax and now lands on it.

## The stomatal conductance, separated from the boundary layer (#122)

On the energy-balance path `gs_H2O` (and `stom_cond_CO2_`) is the TOTAL leaf-to-air conductance, `E*P/D_leaf`. The transpiration crosses the leaf boundary layer as well as the stoma, so comparing it with a measured stomatal conductance compares different quantities. Two read-only accessors separate them. `gb_H2O` is the boundary layer's conductance to vapour, `1.075 * P/(R*Tair) / ra_`, which is plantecophys's convention for a hypostomatous leaf. `gs_stomatal_H2O` is the stoma's own, `1/(1/gs_H2O - 1/gb_H2O)`.

Off the energy balance there is no boundary layer: `gb_H2O` is `Inf` and `gs_stomatal_H2O` is `gs_H2O` bit for bit. On it, over 216 transpiring driver rows, the total conductance reaches 27% of the boundary layer's at 0.5 m/s wind and 13% at the default 2 m/s. At that worst row the stomatal conductance is 37% above the reported total.

⚠️ **`gs_stomatal_H2O` refuses where `gs_H2O >= gb_H2O`.** Then the transpiration is more than the boundary layer can carry at that deficit, and no finite stomatal conductance exists. Hydraulics set the transpiration, so nothing in the solve prevents it. The refusal sits in the accessor rather than in the solve so that it cannot abort a run; no row of the grid above reaches it.

Nothing in the solve reads either accessor, so every result is unchanged and every golden file is bit-identical. The diffusivity ratio is still the stomatal 1.67 along the whole path, which #122's third part changes separately.

# phylloptim 0.15.1

## SOX, JW26 and ProfitMax have an exact first-order condition with the energy balance on (#124)

⚠️ **Results move on the energy-balance path for these three curves; with the gate off they are bit-identical**, and so is every identity-link curve on either path. With the energy balance on, `dA/dpsi` also runs through `Tleaf(E(psi))`. That term was written for the identity link only. On the stem routes a non-identity link refused it, so those three curves fell back to a bracket-width search. On the collar route nothing refused it: the term was added without the link's `h'(A)`, so the collar solve root-found the wrong first-order condition.

The term now takes the same single `h'(A)` as the rest of `dA/dpsi`. There is no `h''` in a first derivative. ProfitMax's thermal cost `TC(Tleaf)`, the one cost that reads leaf temperature, now moves with `psi` as well. It is cost-side, so it takes no `h'`. The same slope serves `lambda_ProfitMax()`, through the new `thermal_cost_temp_deriv()`. The compensation-point branch carries both. `dprofit_dpsi_stem` no longer throws, and `dprofit_energy_balance_term` is now a template on the curve.

**What moves.** On the collar route, at PPFD 1500 and Tair 25 to 35 °C, the old collar sat 0.1 to 0.16 MPa from a dense scan's argmax for SOX and JW26, and 0.5 MPa or more for ProfitMax. It now lands on the scan argmax. No golden file covers that route for these curves. On the stem routes the optimum moves by at most 5.4e-6 relative, the old search's bracket-width resolution. The psi_stem golden file changes on 92 energy-balance rows of these three curves and nowhere else, with no change of status or kind. Profit moves by at most 1.5e-10 relative, which is the objective's own noise floor from the inner `ci` solve. Argmax second differences in air temperature fall by up to 150 times.

**Cost.** Gate off, the collar and stem solves are unchanged, interleaved ×3.

## The closed form's own solvers refuse the energy balance (#116)

**No numbers move**: all three golden files are bit-identical. #131 already made
every deficit in `closed_form.hpp` read `vpd_leaf_` and made `Leaf::optimise_closed`
refuse `use_energy_balance_`. What was still open was the header's own entry points
-- `closed_form::solve`, `solve_exact_beta2`, `solve_CF77`, `solve_TF24_floor` -- which
a bench, a test or a consumer can call directly, and which with the gate on returned a
number rather than a refusal.

That number is wrong, and **`within_guard` does not see it**. Measured against a
4000-point argmax of the exact TF24 objective with the gate forced through, single
layer, PPFD 900, D_air 1 kPa:

| Tair | A exact | A closed form | error | `within_guard` |
|---|---|---|---|---|
| 25 | 16.37 | 18.76 | **+15%** | pass |
| 30 | 11.86 | 16.37 | **+38%** | pass |
| 35 | 4.19 | 10.28 | **+146%** | pass |
| 40 | −3.76 | 2.94 | **+178%** | pass |

Off the gate the same comparison is 0–4%. The guard tests the solution against three
bounds and reads no driver, so a deficit out by the 3–4x `vpd_leaf_`/`atm_vpd_` ratio
leaves ci/ca plausible. The refusal is therefore a stop in the solvers rather than a
case in the guard, and it is ONE function, `require_prescribed_temperature`, which
`optimise_closed` now calls as well -- so the production route and a direct call
cannot drift on what they refuse. The header records the two repairs considered in
#116 (a Picard pass over `D`; linearising `esat(Tleaf)`) and why neither is taken.

`test_closed_form_is_prescribed_temperature_only` pins all four refusals, that
`vpd_leaf_ == atm_vpd_` bit-for-bit off the gate, and the guard's scope: tripling
the deficit after a solve does not change its verdict.

# phylloptim 0.15.0

## A new cost curve, `LeastCost`: cost per unit carbon (Prentice et al. 2014)

The ninth member, and the one whose objective is not a difference. Least-cost
minimises the cost of acquiring a unit of carbon,

```
C_P = (c_E*E + c_V*Vcmax) / A
```

which is a **ratio**. Minimising a ratio is maximising its reciprocal, and a
logarithm separates the two factors, so it joins the family through the benefit
link rather than through a new kind of objective:

```
max  log A - log(c_E*E + c_V*Vcmax)      h = log,  Theta = log(c_E*E + c_V*Vcmax)
lambda = c_E*A / (c_E*E + c_V*Vcmax)
```

`h = log` is the link `SOX` and `JW26` already use. ⚠️ **It reaches that arm by a
different route, and the sign lands elsewhere.** Those two are products `A*g`,
whose log-cost is `-log g`; this is a ratio, whose log-cost is `+log W`. The link
table reports `h' = 1/A` for all three; `cost_deriv` is where they differ.

⚠️ **The only penalty here that carries something other than water.** Every other
cost in this package is a function of potential alone. This one weighs
photosynthetic capacity against transpiration, which is why its first-order
condition reproduces Prentice et al.'s `chi` and why a light sweep on it measures
acclimation rather than the price of water.

```r
leaf_solve(psi_soil = 0.1, PPFD = 3000, atm_vpd = 2, ca = 40, leaf_temp = 25,
           supply = leaf_supply_singlelayer(),
           model = "LeastCost", LeastCost_beta = 146)   # dimensionless
```

`LeastCost_beta` is `c_V/c_E`, the published unit-cost ratio. **One field is
enough and that is not a simplification:** scaling both unit costs by `k`
multiplies `C_P` by `k`, which moves neither the argmax nor `lambda`. ⚠️ It is
**dimensionless**, so both fluxes are weighed in `mol m^-2 s^-1` — a beta derived
against this package's kg H2O and umol CO2 would be out by `5.6e7`.

⚠️ **Unset is refused, and so is zero** — stricter than `TF24_floor`, which
accepts an explicit zero because zero is `TF24` there. Zero here is not a limiting
case: it deletes the capacity term, leaving `A/E`, which rises without bound as
transpiration falls, so the optimiser returns the wet bound for every driver set.

**Checked against the published analytical solution**, which no other curve here
has. Prentice et al.'s first-order condition is

```
(ci - Gamma*) / (ca - ci) = xi / sqrt(D),    xi = sqrt(beta*(Km + Gamma*)/1.6)
```

and the solved `chi` reproduces it to **1e-11** across `D` = 1, 2 and 4 kPa, on a
leaf driven onto the Rubisco-limited branch with `R_d` zeroed. Doubling beta
scales `(ci - Gamma*)/(ca - ci)` by `sqrt(2)` to eight figures. At this package's
own co-limited defaults the agreement is within 1.3%.

⚠️ **Two reasons the widely-quoted `chi = xi/(xi + sqrt(D))` does not match to
that precision, and neither is a defect.** That form is the condition above with
`Gamma*` dropped from the numerator, so it understates `chi` by
`(Gamma*/ca)/(1 + xi/sqrt(D))` — 4% at 1 kPa, 8% at 4 kPa, growing with deficit.
And the `1.6` in `xi` is the H2O:CO2 stomatal diffusion ratio, where this package
uses **1.67**; testing against 1.6 while the solver diffuses at 1.67 leaves a 0.6%
error. The test suite asserts the size of the first gap rather than tolerating it.

**No closed form.** `method = "closed"` is refused, for the reason `SOX` and
`JW26` are: the log link puts `A` into lambda, and `lambda -> xi -> ci -> A`
closes into a fixed point rather than an inversion. Prentice et al. do publish one,
but they obtain it by solving in `ci` directly, where that `A` cancels.

**The collar route solves the right first-order condition on it.** The cost reads
the upstream potential through `E`, so `cost_reads_upstream<LeastCost>()` is true
and `cost_deriv_upstream` carries `-k_max f(psi_up) / W` (in mol); without that row
the collar solve converges on a point that is not the maximum, with every number
plausible (#141). `test_collar_foc_matches_profit` includes the curve, 24 points.

**What a consumer sees.** `LeastCost_beta` is a new `Leaf` field and a new slot in
the parameter enumeration, `par_LeastCost_beta` = 19, placed after
`TF24_floor_lambda_o` so the three model-owned slots stay contiguous; the fifteen
thermal-damage slots and `par_PPFD` move up one, `gradient::n_theta` is **35** and
`phylloptim::n_pars` **36**. Nothing that reads the pack by name moves; a `theta`
written out by hand must gain the entry. `leaf_solve()`, `leaf_batch()` and
`leaf_gradient()` take `LeastCost_beta`; every other model refuses it by name, and
`leaf_gradient()` differentiates it on `model = "LeastCost"` (the R setter writes it
back per perturbation, as C++'s `gradient::apply()` does).

⚠️ **With the energy balance on, `W` moves with leaf temperature through `vcmax_`
and no derivative here carries that term.** On the stem routes every log-link
derivative is already refused with the gate on and the solver scans instead, so the
optimum is unaffected; the collar route's first-order condition omits it. The same
class of omission the other log-link curves carry there; stated rather than fixed.

# phylloptim 0.14.1

## A soil drier than the stem curve's P99 is a closed leaf, not an error

The stem vulnerability curve's knot grid stops at P99 (6.89 MPa at the default curve), and every read past it threw `stem_curve_domain`. On the single-potential stem routes (`leaf_supply_singlelayer()`), that made a soil drier than P99 unsolvable, although the answer there is a shut leaf, and P95 to P99 already gave it. Past the grid the stem curve now reads its closed form, `G(psi) = (b/c) gamma_lower(1/c, (psi/b)^c)`, with slope `exp(-(psi/b)^c)`, and inverts it there; the knots are seeded from the same expression, so value and slope are continuous at P99. `cumulative_vulnerability_integral_derivatives_at` likewise answers past the grid: its series holds out to `(psi/b)^c = 40`, and beyond that it returns the limit `(b/c) Gamma(1/c)` and the limit's derivatives. A negative potential, a flow below the table, and a flow at or past the curve's total `G(inf)` still raise `stem_curve_domain`. The collar route was not affected. No in-grid number moves: the operating-point, psi_stem and primitives golden files and the gradient golden file are bit-identical (#157).

# phylloptim 0.14.0

## Net longwave depends on air temperature and humidity (#97, #28)

⚠️ **Results move on the energy-balance path; with the gate off they are bit-identical**, which is how plant runs its leaf. The net longwave was a constant, `longwave_net_offset = -40` W m⁻². With an atmospheric emissivity it runs from −144 W m⁻² at Tair 20 °C to 0 at 45 °C (at `atm_vpd` 2 kPa), so the constant was out by up to 104 W m⁻² and wrong in opposite directions at the two ends, crossing over near Tair 37 °C. It is replaced by an isothermal net radiation plus a radiative conductance:

```
ea     = esat(Tair) - atm_vpd                              in PASCALS
ema    = min(0.642 * (ea/Tair_K)^(1/7), 1)                 Brutsaert
Rn_iso = absorbed shortwave - (1 - ema) * sigma * Tair_K^4
         (absorbed shortwave = alpha_PAR_ * PPFD / 4.57 + alpha_NIR_ * NIR_int_, as in 0.13.0)
g_rad  = 4 * eps_leaf * sigma * Tair_K^3
Tleaf  = Tair + (Rn_iso - lambda*E) / (vol_heat_cap_air/ra + g_rad)
```

**The balance is still explicit.** Linearising the leaf's own `Tleaf^4` emission about `Tair` puts it in the denominator, so `Tleaf` stays affine in `E` and `dTleaf/dE = -lambda / (g_sensible + g_radiative)` is still a constant. No inner iteration was added, the first-order conditions and their implicit-function derivatives are structurally unchanged, and `leaf_temp_from_E()` is still the one place the expression is evaluated.

**What moves.** `g_rad` is 13–18% of the total conductance at the default leaf size and wind, so it pulls the leaf toward air temperature, while the longwave term cools a cool leaf much harder than −40 did and a hot one less. Solved with `find_root_collar_psi()` at PPFD 900 with daylight NIR, ψ_soil 0.5, D_air 2, on the 0.13.0 shortwave:

| Tair | net longwave | Tleaf | profit | kind |
|---|---|---|---|---|
| 20 | −40 → −144.4 | 23.06 → 20.24 | 5.786 → 7.229 | interior |
| 25 | −40 → −98.5 | 28.10 → 26.31 | 4.920 → 6.026 | interior |
| 30 | −40 → −69.7 | 33.21 → 32.05 | 3.575 → 4.349 | interior |
| 35 | −40 → −44.5 | 38.45 → 37.82 | 1.591 → 1.958 | interior |
| 40 | −40 → −19.8 | 44.28 → 44.10 | −3.574 → −3.559 | solver-refused |
| 45 | −40 → 0 (clamped) | 49.40 → 49.59 | −3.815 → −3.818 | shade-death |

Below Tair ~40 the leaf is 0.6–2.8 K cooler at the solved point (3.0 K at zero transpiration at Tair 20) and profit is 22–25% higher. Above ~42 the two effects oppose and the longwave wins: at 45 the sky term is 40 W m⁻² warmer than the constant was, `g_rad` pulls back only part of that, and the leaf is 0.2 K warmer. No row changes kind.

⚠️ **`ea` is in Pascals.** The same relation is more often quoted as `1.24*(ea/T)^(1/7)` with `ea` in hPa, and the two agree because `0.642*100^(1/7) = 1.2395`. Using 0.642 with kPa gives ~0.31 where the answer is ~0.82: finite, plausible, and it triples the longwave loss. `test_longwave_and_radiative_conductance` asserts the Pa form against the hPa one.

⚠️ **The emissivity is clamped at 1, and above the clamp `Rn_` no longer responds to air temperature or humidity.** Brutsaert's fit passes unity once `ea` exceeds ~22.25 Pa K⁻¹ × Tair_K: above Tair 38.9 °C in saturated air, 43.9 °C at a 2 kPa deficit and 48.0 °C at 4 kPa. Unclamped, the isothermal longwave would turn positive (+100 W m⁻² by Tair 60 at 2 kPa), a sky heating the leaf. Clamped, the term is exactly 0, which is the correct limit, but a heatwave sweep past the clamp sees only `g_rad`'s temperature dependence. An emissivity valid at those temperatures is not part of this release.

⚠️ **`ra_` is now strictly aerodynamic, and `Rn_` is the isothermal net radiation.** A caller that overwrites `ra_` after `set_drivers()` with a combined convective-plus-radiative resistance (Jones et al.'s `r_HR`, or `1/(Gbh + Grad)`) must also set `l$g_rad_ <- 0`, or the radiative conductance is counted twice: `Tleaf - Tair` and `dTleaf/dE` come out 5–7% small and nothing can detect it. Before this release a combined resistance was the right way to add radiative coupling. `g_rad_ = 0` recovers the single-resistance balance, and a test asserts it. `g_rad_` is a new R-bound field; `longwave_net_offset` is deleted.

**Against Sicangco et al. (2026).** The expression reproduces theirs to 1.4e-14 over Tair 20–60 at two deficits. The range #97 quotes for theirs, −98 to −133 W m⁻², is what it gives at a fixed `ea` of 1.0 kPa (−98.41 at Tair 20, −132.05 at 40), not at a fixed `atm_vpd`: phylloptim's driver is the deficit, and at a fixed deficit `ea = esat(Tair) - atm_vpd` climbs steeply with Tair, so such a sweep is drier than `ea = 1` kPa at the cool end and far more humid at the hot end.

**Cost.** With the gate on a collar solve is unchanged (interleaved ×3). With it off, `set_physiology` gains one `pow` and one `exp`, about 0.026 µs, under 1% of a `set_physiology` plus collar solve, with an identical checksum.

**Golden files.** `operating_points.tsv` and `gradient_golden.tsv` are bit-identical (gate off). `primitives.tsv`: 34 of 544 values move, all `leaf_temp_from_E` and `dleaf_temp_dE` in the arithmetic tier; every other tier is exactly 0. `psi_stem_optima.tsv`: 1708 of 5184 rows move, every one with the energy balance on; all 2592 gate-off rows are byte-identical, and no row changes status or kind.

# phylloptim 0.13.0

## Light arrives intercepted and the leaf absorbs it

The leaf's light input is now named for what it is, and absorption happens inside the leaf. Breaking: the quantum-yield trait is renamed and rebased, and the energy-balance path needs one new input.

* `PPFD` (the `set_drivers()` / `leaf_solve()` driver, `set_physiology()`'s second argument, the `par_PPFD` slot) is the **intercepted** photon flux per unit leaf area: incident on the leaf, before its absorptance. Its value and units are unchanged; the documentation now says which stage of the light chain it is.
* New leaf properties `alpha_PAR_` (0.85) and `alpha_NIR_` (0.20), the PAR and near-infrared absorptances. Settable fields, not fitted traits.
* `a` -> `phi_J`, the quantum yield of electron transport per **absorbed** photon. Electron transport reads `phi_J * alpha_PAR_ * PPFD`. The default `0.30 / 0.85` restates the old incident-basis 0.30, so with the defaults J is unchanged to rounding. Migration:
  * `leaf_traits(a = x)` -> `leaf_traits(phi_J = x / 0.85)`
  * `Leaf(..., a = x, ...)` -> `Leaf(..., phi_J = x / 0.85, ...)`
  * `l$set_traits(..., a, ...)` (8th argument) -> the same position, `phi_J = a / 0.85`
  * `l$a` -> `l$phi_J`
  * `leaf_gradient(pars = "a")` / gradient columns named `a` -> `"phi_J"`; the derivative is per unit `phi_J`, so `d/da = d/dphi_J / 0.85`
  * A gas-exchange fit identifies only `phi_J * alpha_PAR_`: fix `alpha_PAR_` and the fitted `phi_J` is per absorbed photon.
* **Energy-balance path:** net radiation is the absorbed shortwave of two intercepted bands, `Rn = alpha_PAR_ * PPFD / 4.57 + alpha_NIR_ * NIR_int_ - 40`, where it was `2 * PPFD / 4.57 - 40` (all shortwave absorbed). New per-call field `NIR_int_`, the intercepted near-infrared in W m^-2, set before `set_drivers()` like `wind_speed_`. **It has no default**, because the right one depends on the light source (about `PPFD / 4.57` in daylight, 0 under an LED), so `set_physiology()` refuses a non-finite `NIR_int_` with `use_energy_balance_` on. Off the energy-balance path it is never read. Migration: `l$use_energy_balance_ <- TRUE` -> also `l$NIR_int_ <- PPFD / 4.57` for daylight.
* `constants.hpp`: `sw_abs_per_par` (2.0) is removed. It folded the PAR share of shortwave and an absorptance of 1 into one number.

Golden files regenerated. Rows with the energy balance off move by reassociation only (operating points: 342 cells, worst 2.0e-10; psi_stem optima: 836 cells on 133 rows, worst 1.8e-10). Energy-balance rows move by the physics: with daylight NIR at the grid's PPFD of 1500, absorbed shortwave falls from 657 to 345 W m^-2 and the solved leaf temperature by a median 7.9 C (up to 14.2 C); 14 of 2592 rows change kind. The 17 `leaf_temp_from_E` primitives move for the same reason.

# phylloptim 0.12.0

## An optional thermal-damage layer on the leaf, off by default

`Leaf$use_thermal_damage_` turns on a lasting thermal-damage model, ported from the layer plant carried in its own leaf (traitecoevo/plant#566). Three fractions of leaf area, named for what they do — functional F ⇌ impaired I → damaged D, with F = 1 − I − D — plus an acclimation state A, all driven by LEAF temperature:

```
dA/dt = acclim_alpha_ * softplus(T - acclim_T_; acclim_s_) - acclim_beta_ * A
dI/dt = k_imp(T) * F - (k_rec_eff(T) + k_mat_) * I
dD/dt = k_mat_ * I
k_imp(T) = k_i_ * exp[(k_i_Ea_/R)(1/T_ref - 1/T)] * K/(1 + K),  K = exp[(impair_H_/R)(1/Tmid - 1/T)]
```

Rates are per day. The Vcmax and Jmax peaked-Arrhenius curves stay the instantaneous, fully reversible temperature response; impairment is rate-and-state and lasts until repaired, and the two share no parameter. The switch K/(1+K) has its own enthalpy `impair_H_` and midpoint `impair_Tmid_` (defaults 200 kJ mol⁻¹ and 200000/650 K, the Jmax curve's deactivation term), moved by acclimation through its own `dTimpair_max_`.

- **On the solve.** The caller seats `A_`, `I_` and `D_` before each `set_drivers()`. Both capacities are scaled by F (floored at 1e-6), the Vcmax and Jmax optima move by `topt_offset_ + dTopt_max_ * A/(A + K_A_)` through the entropy term, and four maintenance costs, zero by default, are added to `R_d_`, which is otherwise not scaled. `k_imp_eff_`, `k_rec_eff_` and `functional_` are reported at the leaf's temperature; on the energy-balance path that is the operating point's.
- **The dynamics, on the leaf.** `$damage_rates()` evaluates the right-hand side and `leaf_damage_trajectory()` integrates it over a leaf-temperature series (RK4, sub-stepped so the fastest rate × step < 0.1 and ≤ 0.25 K per step), so a leaf experiment can be fitted with this package alone. A caller integrating the states itself — plant's TF24t, which adds dilution of D by leaf turnover and growth — calls `$damage_rates()`.
- **Gradients.** `leaf_gradient()` and `leaf_gradient_batch()` carry the layer: on a leaf with the gate on, its fifteen inputs that reach the solve (`A_`, `I_`, `D_`, `topt_offset_`, `dTopt_max_`, `K_A_`, `k_i_`, `k_i_baseline_`, `k_rec_`, `m_rep_`, `t_rep_cut_` and the four costs) are gradient parameters, appended to the enumeration (`n_theta` 19 → 34) and refused by name with the gate off. `leaf_damage_gradient()` gives the trajectory's exact forward-mode sensitivities to all sixteen dynamics parameters and the three initial states. The derivative surface (`collar_at`, `outputs_at`) and the closed form still refuse the layer.

Off by default, and off it reads none of its inputs: every golden file is bit-identical, `k_imp_eff_` and `k_rec_eff_` are zero and `functional_` one.

⚠️ **`n_theta` is 34 and `n_pars` 35.** A C++ caller filling `theta[gradient::n_theta]` by hand must recount: the layer's slots are read only with the gate on, but a short initialiser shifts nothing else only while it stays short.

Respiration departs from the original layer on purpose: there `R_d` was `0.015 * Vcmax` after the damage discount, so a damaged leaf respired less; here `R_d` stays `R_d_25` on its Q10 (#86) and the layer only adds its costs.

# phylloptim 0.10.3

## odelia is a floor, not a pin, at 0.8.0

`DESCRIPTION` no longer names an odelia commit in `Remotes:`; `LinkingTo: odelia (>= 0.8.0)` states the oldest release this package builds against, and the R check and the C++ workflow both build against that release (tag `v0.8.0`). A consumer such as traitecoevo/plant pins odelia once, in its own `Remotes:`, and this package no longer has to move in step with every odelia it takes: a `Remotes:` pin here was inherited by every consumer, and pak refuses an install whose packages pin one sibling two ways. odelia 0.8.0 changes the ODE step-size default (traitecoevo/odelia#77), which this package does not use; of the headers it does use, only a spline read at an active position changes, and that changes an adjoint only to rounding. Nothing in the package's code changes.

# phylloptim 0.10.2

## `Remotes:` tracks odelia 0.7.3

`LinkingTo: odelia (>= 0.7.3)`; `Remotes:` and the C++ workflow's odelia checkout pin `88cbff7`, odelia 0.7.3 on master: its replay fixes (traitecoevo/odelia#71, #73) and a sweep across steps with no state (traitecoevo/odelia#75), which traitecoevo/plant's FF16 fitness gradient needs. Nothing in this package changes; the pin moves so that a consumer pinning both packages sees one odelia.

# phylloptim 0.10.1

## `Remotes:` tracks odelia 0.7.1

`LinkingTo: odelia (>= 0.7.1)`; `Remotes:` and the C++ workflow's odelia checkout pin `e9a179e`, odelia's fix for a run that starts with no state (traitecoevo/plant's FF16 stand gradient needs it). Nothing in this package changes; the pin moves so that a consumer pinning both packages sees one odelia.

# phylloptim 0.10.0

## `Remotes:` tracks odelia 0.7.0

`LinkingTo: odelia (>= 0.7.0)`; `Remotes:` and the C++ workflow's odelia checkout pin `6d32329`, the merge of traitecoevo/odelia#59 (reverse mode, 0.7.0). The rows below are taken on odelia 0.7.0's `implicit_node.hpp`, `tangent.hpp` and `value_with_slope.hpp`.

## A replayed boundary point hands over its solve's rows

`replay_operating_point(collar, kind)`, which puts a recorded operating point back, lost two things a pinned solve has and a replay does not compute: the polished root of the zero-flux bound, which `bound_at` corrects a `boundary-soil` row back to, and the collar itself where that polished root sits outside the unpolished bracket, which the replay clamped back in. Rows from a replayed `boundary-soil` point therefore differed from the solve's: over 1,048 rows on a grid of dry and hot cases, 258 of the 643 larger than 1e-6 moved by more than 1e-4, worst 25×, with sign flips on `root_c` and `root_P50` at one-layer wet bounds. plant re-solves every recorded leaf through this call during a reverse sweep, so its gradient took those rows. The replay now re-takes the polish from the same bracket and evaluates at the recorded collar; replayed rows are bit-identical to the solve's at interior, `boundary-soil` and `boundary-crit` points, and `test_leaf.cpp` asserts it. Forward values are unchanged.

## The derivative surface is documented in the package

`vignette("derivative-surface")` explains the three calls, the eleven operating-point kinds and what `collar_at` does at each, why the derivative comes from the implicit function theorem, and how it differs from `leaf_gradient()`. `leaf_model.hpp` opens with the same map and the preconditions. The README is retitled and gains a section on it; the developer guide gains a map of the headers.

## ⚠️ Breaking: the leaf supplies its own derivative rows, and three counts moved

`Leaf`'s derivative surface is a set of member templates on the caller's scalar
(`supply_draw_at`, `collar_at`, `outputs_at`), and supplies the derivative of
its operating point from the implicit function theorem at the converged point,
instead of being recorded and differentiated. What that changes for a caller:

* **`vulnerability_curve_ncontrol` defaults to 400, was 100.** The
  pre-integrated vulnerability tables are built on four times the knots, and
  every number read off them moves. `leaf_control()` carries the new default.
  ⚠️ The number is written twice -- `Leaf::ncontrol_default` in
  `leaf_model.hpp`, which plant's `Control` reads, and `leaf_control()`'s formal,
  which R has to pass because only the 15-argument constructor is bound. plant's
  `test-control.R` compares the two with `expect_identical`; nothing in this
  package does, and the golden file cannot see the difference.
* **`n_pars` is 20, and the fitted length is `n_theta` = 19.** The pack the
  kernels read also holds `par_PPFD`, which is seated per observation from the
  drivers and is never fitted, so it has a slot and no name. Read `theta` out to
  `n_pars` and you run past its end.
* **odelia `>= 0.7.0`** is required: the supplied rows are built on its
  `implicit_node.hpp`, `tangent.hpp` and `value_with_slope.hpp`, which 0.6.0
  does not have, and read through its Hermite interpolant (`hermite_spline`,
  odelia 0.5.0).
* **The four golden files are regenerated** (macOS/arm64, odelia 0.7.0), for
  the 400-knot default. Against the 0.9.1 baselines: `operating_points.tsv`
  480 of 576 rows move, 108 beyond the cross-platform tolerance, all near-shut
  rows at `psi_soil = 4, ppfd = 100`; `psi_stem_optima.tsv` 1539 of 5184;
  `primitives.tsv` 129 of 544, with the arithmetic, vulnerability and
  assimilation tiers bit-identical and the spline tier the lowest that moves
  (88 of 110, at 3.5e-8); `gradient_golden.tsv` 15 of 20 rows, worst 1.6e-3
  relative.
* **Two double-path changes a caller can see.** `evaluate_root_collar_psi` and
  `dprofit_droot_collar_psi` dispatch on the seated cost curve, where they had
  been pinned to TF24's, so a CF77 or TF24_floor leaf evaluates its own curve
  there. `ci_at_compensation_point()` is derived from the kernel rather than
  carried as a flag, which decides which branch `dprofit_at_collar_psi` takes at
  a zero-flux probe. The curves read the interpolant directly rather than
  through odelia's refusing front end, so each states its own bound where it
  applies it -- `root_vuln_at` clamps the argument, `root_vuln_integral_at`
  caps the value, `eval_stem_curve` raises.
* **A pinned optimum's bound is solved to collar precision.** The zero-flux and
  continuity bounds were placed to 1e-4 MPa, which is fine for a bracket but not
  where the bound IS the answer. There, the loose root-finder's stopping point
  moved with the parameters differently from the root, so the supplied rows,
  which differentiate the root, disagreed with a finite difference of the solve.
  The worst case was root-trait profit rows 100× off near shutdown. A pinned
  solve now polishes its bound, and the wet-bound row corrects its residual from
  the step-in point back to the root. Against a central-difference referee, bad
  rows fall from 72 to 0 of 2640, and from 183 to 12 of 7784 on a grid across
  deficits. Those 12 are within 5e-8 absolute at one-layer wet bounds, plus one
  interior `curv_elec` row off by 2e-4 relative, which predates this change.
  Pinned and shade-death results move: collar potential by at most 1.3e-4 MPa,
  profit by -1.3e-5 to +5.6e-4, and 112 of 576 golden rows. Interior rows are
  bit-identical. The two bracket root-finds now share their endpoint values, so
  the solve is no slower (3.01 against 3.03 µs).
* **`leaf_gradient(method = "ift")` refuses a pinned optimum by name**, with "a
  pinned optimum" in the message, in R and in the batch. It used to fail because
  a centred step would not fit beside the loose bound, and with the bound solved
  exactly it can.

⚠️ **This retracts two earlier entries.** "`n_pars` is unchanged" under the
`Tleaf` reporting entry, and "`n_pars` is unchanged at 19" under the
`shadow_cost` entry, were true when written and are not now: the count is 20 and
the one a caller wants is `n_theta`.

# phylloptim 0.9.2

## `Remotes:` tracks odelia 0.6.2

`Remotes:` and the C++ workflow's odelia checkout move from `815a54f` (0.6.0) to `64ae092`, the merge of traitecoevo/odelia#69 (0.6.2). plant pins the same odelia as phylloptim and needs 0.6.2, whose solver no longer pays the square of the state length on every resize (traitecoevo/plant#656), which is why this release exists. No code changes. Between 0.6.0 and 0.6.2 odelia added a steady-state solver (0.6.1) and moved when RODAS allocates its buffers (0.6.2); neither changes any arithmetic, so the golden baselines are unchanged. The `LinkingTo: odelia (>= 0.2.0)` floor is unchanged.

# phylloptim 0.9.1

## `Remotes:` tracks odelia 0.6.0

`Remotes:` and the C++ workflow's odelia checkout move from v0.4.0 to `815a54f`, the merge of traitecoevo/odelia#63 (an R-side ODE stepper, which regnans steps the canonical equation with). plant pins the same odelia as phylloptim and must move to 0.6.0 for regnans, which is why this release exists. No code changes.

The golden baselines move with the pin, as they did on the `release/0.8.x` line (0.8.2, #139): odelia 0.5.0 made the spline behind the soil and leaf lookups a local cubic Hermite, which reads the same natural cubic spline as before with values alone, so the goldens differ at rounding only — `operating_points.tsv`, `psi_stem_optima.tsv`, `primitives.tsv` and `gradient_golden.tsv` are the 0.8.x line's, re-verified bit-identical here against the 0.6.0 headers. `test-gradient.R`'s operating point is the 0.8.x line's as well. The `LinkingTo: odelia (>= 0.2.0)` floor is unchanged.

# phylloptim 0.9.0

## Fix: the collar solve maximised the wrong objective for four cost curves

On the collar route the solve chooses the root-collar potential, and `psi_stem` follows from it through the supply. Every cost is billed through `psi_stem`, and the first-order condition carries that channel. Four costs also read the collar itself, because they price the stem path between collar and leaf rather than the state of the leaf end. **CF77** and **TF24_floor** charge for the flow `E = kmax·(G(psi_stem) − G(psi_collar))`. **JS22** charges for the drop `psi_stem − psi_collar`. **CMax** integrates its price from `psi_collar` to `psi_stem`. That direct dependence was missing from the first-order condition, so for these four the collar solve converged, to `|dprofit| ~ 1e-15`, onto a point that was not the maximum of the `profit_` it reported.

The other four price the leaf end alone: **TF24** by the stem conductivity lost at `psi_stem`, **SOX** and **JW26** by a reduction factor of `psi_stem`, **ProfitMax** by `k(psi_stem)` against a fixed soil conductance. They see the collar only through `psi_stem` and are unaffected. So are all eight stem routes, which hold the upstream end at `psi_soil`.

The missing term is `dC/dpsi_upstream`, now a table beside `cost_deriv<K>()` (`cost_reads_upstream<K>()` and `cost_deriv_upstream<K>()`), and appended to the collar FOC so the curves that do not read upstream compile to the same code: `bench_solve` is byte-identical, and all three golden baselines and `gradient_golden.tsv` are bit-identical, because none of them solves an upstream-reading curve on the collar route. `test_collar_foc_matches_profit` is the first test that does.

Measured over 24 driver rows per case (1 and 3 layers), every moved row's `profit_` rose, so the old point was strictly worse by the model's own objective:

| curve | rows moved | worst collar shift, MPa | worst transpiration change |
|---|---|---|---|
| CMax | 24/24 | 0.346 | 34% |
| JS22 | 24/24 | 0.269 | 19% |
| TF24_floor, lambda_o = 1e4 | 24/24 | 0.057 | 5.9% |
| TF24_floor, lambda_o = 1e3 | 24/24 | 0.0054 | 0.5% |
| CF77, lambda = 1e4 | 0/24 (all pinned) | 0 | 0 |

A trait gradient on these curves moves with the solve, since `leaf_gradient()` reaches the FOC through `dprofit_droot_collar_psi`.

## ⚠️ Breaking: root layer thickness is per layer, and it was a 3.7× error

`root_network_from_carbon()`'s C++ signature takes a **vector** of layer
thicknesses where it took one scalar:

```
// was
root_network_from_carbon(carbon, double dz, beta_R_H, beta_R_V [, out])
// now
root_network_from_carbon(carbon, const std::vector<double>& dz, beta_R_H, beta_R_V [, out])
```

`layer_thickness()` → **`layer_thicknesses()`**, returning the per-layer widths
implied by a cumulative depth profile (consecutive differences, implicit 0 at the
surface). Both scalar versions are **deleted** rather than kept as overloads, so a
stale call is a compile error instead of silently binding the old model.

**The R signature does not change.** `root_network_from_carbon(carbon,
soil_depth, ...)` has always taken the profile and derived the widths itself.

**Why this is a fix and not a generalisation.** The vertical resistance is
`beta_R_V * dz^2 / c_r_V`, and one scalar `dz = D/n` is correct only for a profile
of equal layers. What fails otherwise is **discretisation invariance**: slicing a
soil column is a numerical choice, so no property of the plant may depend on it.
With root density uniform over depth the layer-integrated carbon goes as `dz[i]`,
so the per-layer form gives a total `sum(r_R_V) = 3*beta_R_V*D^2/C` for any
slicing — and a single scalar does not. Total vertical root resistance over 1.5 m
carrying 20 kg C m^-2 leaf:

| profile | scalar `dz` (0.8.0) | per-layer `dz[i]` (0.9.0) | ratio |
|---|---|---|---|
| 5 equal layers | 3172.5 | 3172.5 | 1.000 |
| 3 equal layers | 3172.5 | 3172.5 | 1.000 |
| 5 cm surface layer | 6059.5 | 3172.5 | **1.910** |
| 2 cm surface layer | 11688.4 | 3172.5 | **3.684** |

So a graded profile inflated total vertical root resistance up to 3.7×, throttling
uptake, with nothing anywhere reporting it. Unreachable before now only because
nothing built a graded profile; plant #626 makes them reachable, which is what
brought this up.

**Bit-identical for every profile reachable before this change.** All three golden
baselines unchanged (576 operating points, 5184 psi_stem optima, 544 primitives),
and `gradient_golden.tsv` bit-exact. Equal layers give equal widths exactly, and
the golden grid's 1 m boundaries difference exactly.

⚠️ **Two floating-point facts, both asserted in the suite, because they decide how
a caller should reach this function.** Passing widths a caller already holds is
bit-exact against 0.8.0 in 180 of 180 `(depth, n)` pairs. *Differencing* a profile
built as `(i+1)*depth/n` is not — it loses a bit in 122 of those 180, the package
defaults among them. So **a caller holding its own widths must pass them and must
not re-derive them through `layer_thicknesses()`**; the helper is for a caller who
genuinely has only a profile. The one R export therefore is not bit-exact against
0.8.0 for a uniform *non-integer* profile; no golden value moves, because every
baseline here uses integer-metre boundaries.

## `MultiLayerRoots::dz_` is gone

A scalar layer thickness, re-derived from the profile on every `set_soil_state`,
that nothing in this package had read since #33 — and that #626 removed the last
reason to carry, since thickness is now stated by the caller per layer rather than
being something this object holds a second opinion about. Dropped from
`inst/RcppR6_classes.yml`, so `leaf$dz_` no longer exists on the R side; read the
widths off the profile with `layer_thicknesses()`, or off the network via `r_R_V`.

Downstream (plant) must drop the matching `dz_` line from its own
`RcppR6_classes.yml` in the same coupled change.

## Two test-suite bugs the new guards surfaced

`layer_thicknesses()` refuses a profile that is not strictly increasing, and that
immediately caught `test_multi_layer_soil`, which had passed `{0.5, 0.5, 0.5}` as
a cumulative profile — three layers all bottoming at 0.5 m, so two of the three
midpoints sat at the same depth. It went unnoticed because that test asserts
finiteness and sign only. The profile is `{0.5, 1.0, 1.5}` now.

`test-surface.R`'s "reproduces the leaf's own layer thickness" **asserted the
bug**: on a graded profile of 0.4, 0.5 and 0.6 m layers it required all three
`r_R_V` to be equal, computed from a column-average thickness. Replaced by a
per-layer expectation plus the invariance test above.

⚠️ The strictly-increasing guard does **not** catch widths passed where boundaries
were wanted — a width vector that increases with depth passes every check and
yields plausible, different widths. That limit is asserted rather than left
implicit, so the guard is not mistaken for protection it does not give.

## `Remotes:` tracks odelia v0.4.0

A package's `Remotes:` are honoured wherever it is installed from source, so this pin is not private to phylloptim's own CI. Anything depending on both phylloptim and a different odelia ref ends up with two pinned refs for one package in a single graph, which `pkgdepends` refuses with an unexplained `dependency conflict`. That blocked plant's fix for traitecoevo/plant#642, which needs odelia 0.4.0 for step rejection on the pinned-time path.

The `LinkingTo: odelia (>= 0.2.0)` floor is deliberately unchanged: phylloptim uses nothing added since 0.2.0 and should stay installable against older odelia. It is the pin that has to track, not the floor.

Also released as [`v0.8.1`](https://github.com/traitecoevo/phylloptim/tree/release/0.8.x) — the same one-line change off `v0.8.0`, so plant could take odelia 0.4.0 without also taking 0.9.0's breaking root-thickness change.

# phylloptim 0.8.0

## `leaf_solve()` reports the seated curve's lambda

`lambda` is `marginal_cost_water()` — TF24's price at this operating point,
**whatever curve ran** — so on any other curve it answers a question nobody asked.
The seated curve's own `(dC/dpsi)/(dE/dpsi)` is `lambda_emergent`, and it was
reachable only off a `leaf_model()` object. Every `leaf_solve()` caller on a
non-TF24 curve was therefore reconstructing it, and on `TF24_floor` that meant
knowing that `lambda` carries only the hydraulic half and adding the price floor
back by hand:

| curve | `lambda` | `lambda_emergent` |
|---|---|---|
| TF24 | 55757 | 55757 |
| CF77, price 1e5 | 26012 | **100000** |
| SOX | 23354 | 108864 |
| TF24_floor, floor 5e4 | 32996 | 82996 |

CF77's row is the one that pins the field: its emergent lambda is known
independently of the solve, being the constant it was given, and it comes back as
that constant. Appended to `operating_point()` after `shadow_cost`, for the reason
`Tleaf` and `shadow_cost` were appended — those names are positions, so an
insertion beside its namesake would shift every saved output.

## ⚠️ A stale `lambda_emergent` after a shut-down row, found by reporting it

Two terminal exits — hydraulic shutdown and shade death — wrote every other flux
output and left `lambda_emergent_` alone, so a **reused** `Leaf` reported the
previous row's price for a leaf that moves no water. `leaf_solve()` reuses one Leaf
across rows, so a batch containing a dry row returned a finite, plausible price
there — 3.49e+05 — where the same drivers on a fresh leaf give NA.

Both exits now write the NA sentinel, with the reasoning recorded beside the block
in `set_shutdown_state` that already existed for exactly this hazard on the flux
outputs (plant #577, #578). `shadow_cost` needs no equivalent: it is computed from
`transpiration_`, which those exits zero, so it follows to zero on its own.

The bug predates this change and was latent because the field was reachable only
off an object the caller had driven themselves. Reporting it through `leaf_solve()`
is what made it visible, and `test-surface.R`'s row-order test caught it
immediately — the test whose own comment says it crosses into shut-down because
"that is where a last-written-value bug would show up as a duplicated row".

No golden file moves; `n_pars` is unchanged. The C++ assertion that read `Tleaf`
off the end of the reported vector is now two from the end, and the R name/binding
alignment test is what catches that drifting again.

# phylloptim 0.7.0

## A new cost curve, `TF24_floor`: TF24 with a price of water at the wet end

Every conductance-loss cost in this package prices water at **zero** as
transpiration goes to zero. That is not a choice anyone made — it is what "the cost
of water is the conductivity you lose to move it" implies, since no conductivity is
lost when nothing flows — and its consequence is that water is free precisely when
it is abundant. `TF24_floor` is the repair, applied to our own cost.

It takes the split any cost admits, into a part depending on the potential alone
and a part linear in the flux:

```
Theta(E) = Theta~(psi) + lambda_o * E,   lambda = Theta~'(psi)/K(psi) + lambda_o
```

with `K = kmax * f(psi)` the stem conductance, and `Theta~` taken to be **TF24's own
cost at TF24's own traits**. So `lambda_o` is the curve's only parameter, `TF24` is
the same curve at `lambda_o = 0` at identical parameter values, and the comparison
between them is a one-restriction test rather than a fit of two
differently-parameterised models.

```r
leaf_solve(psi_soil = 1.5, PPFD = 1500, model = "TF24_floor",
           supply = leaf_supply_singlelayer(),
           root_network = series_resistance(1e4),
           TF24_floor_lambda_o = 1.5e5)   # umol C (kg H2O)^-1
```

**Two reductions, asserted BIT-FOR-BIT** over 45 driver rows and eight reported
fields, with no fixture and no tolerance: `lambda_o = 0` gives `TF24`, and
`TF24_cost_scale = 0` gives `CF77`. Each term is the parent's own expression, so
zeroing a parameter adds an exact zero to the other curve's exact value — and that
survives fused multiply-add, since the fused addend is then an exact zero.

**`TF24_floor_lambda_o` is an INPUT with no default**, like `CF77_lambda_`, and the
reason is sharper: at zero this curve **is** `TF24`, so a caller who seated it and
never named a price would be running the production model under a new name.
Omitting it is refused, with the units and the routes that can supply one; an
explicit zero is accepted.

**`shadow_cost`, a new reported output, because the two terms are not the same
kind of thing.** `Theta~` is a **realised** cost: carbon actually forgone, in
damaged tissue and in capacity that has to be rebuilt, and it belongs in a carbon
budget. `lambda_o` is not a cost at all. It is a **shadow price** — the value of
water in its best alternative use, which for a leaf is assimilation later — and no
carbon is lost when the plant pays it. It changes the aperture the plant chooses
and nothing else, which is what a multiplier does.

`profit` still reports the objective, so it deducts both. The new column reports
the part that is a price, and the two quantities a caller wants follow:

```
realised cost   = hydraulic_cost - shadow_cost
carbon profit   = profit + shadow_cost
```

⚠️ **A consumer that grows a plant on this leaf wants the second of those, never
`profit`.** Deducting a shadow price from a carbon budget taxes growth by something
the plant never spent. At `lambda_o = 1e5` the objective understates the carbon the
leaf actually kept by 28%.

⚠️ **Zero on every other curve means "this curve does not separate the two", not
"this curve's cost is all realised carbon".** `CF77` is the case worth stating: its
whole cost is `lambda * E`, and reading that as a shadow price is an ordinary thing
to do — but the model supplies one number and nothing to attribute it with, and
behaves identically whether that number is read as carbon lost to water already
taken or as the value of water withheld for later. Reporting it as all-shadow would
present one reading as a fact; as zero-shadow, the other. The field is defined only
where the *curve* defines it, and that is itself an argument for the two-term form:
it is the first curve here in which the question "is `lambda` a cost to deduct, or a
price that only shapes behaviour?" can be posed at all.

The sharpest test of the column is the pair that was already here. `CF77` at
`lambda` and `TF24_floor` at `TF24_cost_scale = 0` are the same model numerically —
every state and flux bit-identical — and `shadow_cost` is the only column that
separates them.

Appended to `operating_point()`, after `Tleaf`, for the reason `Tleaf` was: those
names are positions and an insertion would shift every saved output. It reads the
stored transpiration rather than recomputing it, so `shadow_cost == lambda_o * E`
holds bit-exactly against the reported `E`. No golden file moves, since those are
written from their own struct, and `n_pars` is unchanged at 19 — this is a reported
output, not a parameter.

**A real prediction, and a test of it.** Scaling `kmax` and the vapour deficit
together by 16 leaves stomatal conductance untouched at every potential, so a cost
whose marginal value depends on `psi` alone is invariant and one that prices the
FLUX is not. Measured: `TF24` moves by 0, `CF77` by 64%, and `TF24_floor` in
between and monotonically in `lambda_o` — 4.5e-02, 3.8e-01, 5.7e-01 at `lambda_o`
of 1.5e3, 1.5e4 and 1.5e5.

**The closed form serves it, and its tail is worse than the other two's.** The
wet-end `lambda` IS `lambda_o`, a constant — the `n = 0` case, reached for a reason
belonging to the MODEL rather than to the vulnerability curve — so
`set_model("TF24_floor", "stem", "closed")` works and runs at 1.40 against 2.90
µs/call, 2.11x realised at phi = 0.375. But the leading order discards the whole
TF24 term, and the shared `ci/ca` guard does not see the resulting error: over a
150-row grid the worst SERVED error is 2.9e-01 at `TF24_cost_scale = 0.5` (against
`CF77`'s own 2.4e-01) and **5.2e-01 at the DEFAULT scale**. Read `solve_TF24_floor`'s
table before using that arm. At `lambda_o = 0` it is refused by name and points at
`TF24`, which has its own start.

## `leaf_solve()` and `leaf_batch()` can set the prescribed prices

`$CF77_lambda_` and `$TF24_floor_lambda_o` are fields rather than traits, so
`leaf_traits()` could not carry them — and `leaf_solve()` builds its own `Leaf`
internally. The two priced curves were therefore reachable from `leaf_model()` and
**unreachable from the one-call surface**: `leaf_solve(model = "CF77")` could only
ever raise that curve's own "needs `CF77_lambda_` set" refusal. That was survivable
while CF77 was the only such curve and its documentation said "build the leaf
yourself"; it stopped being survivable for a curve whose only parameter is a price.

`leaf_solve()` takes `CF77_lambda` and `TF24_floor_lambda_o` now, matching
`leaf_batch()`, and applies them where the model is seated — so they survive
`set_drivers()` per row on both the `reuse = TRUE` and `reuse = FALSE` paths.
Passing a price the seated model does not read is **refused rather than ignored**,
naming the curve that does read it: silently ignoring it is how someone spends an
afternoon wondering why their lambda had no effect.

## Smaller things

`.gradient_setter()`'s positional trait call is a `do.call` over the derived trait
vector rather than fifteen subscripts written out. The comment beside it predicted
that adding a trait would break it at run time inside the generated binding,
"naming neither this line nor the count" — which is how #41 broke — and an
intermediate version of this work added one and broke it again, in exactly that
way: 118 gradient-batch rows reporting `error` where they had reported `interior`,
because the R reference threw and the batch did not. The trait was withdrawn, so
the shipped trait vector is unchanged at fifteen; the fix is kept, because the next
trait will not be so easy to withdraw.

⚠️ **`operating_point()$lambda` is TF24's, whatever curve is seated.** It is
`marginal_cost_water()`. The per-curve number is `$lambda_emergent`. That is not
new, but `TF24_floor` is the curve where reading the wrong one is easiest — its
hydraulic half really is TF24's — so the field list now says so, and the R
reduction test pins the one column where the two legitimately disagree.

`gradient_par_names()` goes 18 -> 19: `TF24_floor_lambda_o` is appended after
`CF77_lambda_`, and nothing before it moves. `psi_stem_optima.tsv` grows 4608 ->
5184 rows; every pre-existing row is byte-identical, which is the check the
append-only rule exists for. The 576 operating points, the 544 primitives and the
recorded R gradient baseline are untouched.

## A third axis on `set_model()`: the closed form as a selectable method

`set_model()` took a cost curve and a route. It now takes a **method** as well,
and `"closed"` reaches the same optimum by a different route:

```r
l$set_model("TF24", "stem", "closed")
l$optimise()
l$closed_form_fallback_fraction()   # phi, which sets the realised speedup
```

Given `lambda`, every model here collapses to the Medlyn relation
`ci/ca = xi/(xi + sqrt(D))`; where `lambda` is a wet-end power law in `psi` that
inverts explicitly, so a leading-order potential plus one Newton step on the
supply-minus-demand residual replaces the search over the bracket. `closed_form.hpp`
has existed for a while and nothing called it; it is reachable now, through the same
infrastructure as everything else.

**`method = "exact"` is the default and is bit-identical.** `optimise()` returns
to its previous body before anything new runs: all 4608 golden `psi_stem` optima,
576 operating points and 544 primitives are unchanged.

**Implemented for `TF24` and `CF77`; refused, with reasons, everywhere else.** The
existence criterion is that `h'(A)` — the benefit link's derivative — must not
depend on the solution. `JS22`, `CMax` and `ProfitMax` satisfy it, so a closed
form *exists* for them and only an analytic `dlambda/dpsi` is unwritten. `SOX` and
`JW26` do not: their link is `log`, `h'(A) = 1/A`, so `lambda` carries the
assimilation the solve is for and the dependency closes into a fixed point. Those
two **cannot** have one. Also refused on the collar route, on a layered soil, and
with `use_energy_balance_` on — the last being the resolution of #116, since the
deficit then depends on the transpiration being solved for.

Three things worth knowing before using it:

- **It buys about 2x, not an order of magnitude.** 2.65 -> 1.16 µs for `TF24` and
  2.74 -> 1.37 for `CF77`, both arms in one process. The header's old 10.8x and 47x
  predated the current solver and priced the inversion without the evaluation that
  writes a caller's outputs. The validity guard tests an OUTPUT, so a rejected row
  falls back and the realised figure is `1/[phi + (1-phi)/s]`; `phi` is measurable
  now (`$closed_form_fallback_fraction()`), which it was not before.
- **Accuracy is set by the Medlyn collapse, not by the power law.** The Newton step
  converges — 2, 3 and 5 steps agree — and the residual's fixed point is still not
  the optimum, because `3*Gstar` is the electron-transport-limited `dA/dci`. At
  25 C the error in `A` runs 0.2 % above `ci/ca` 0.9 and 2 % at 0.5-0.6; below
  `ci/ca` about 0.3 the fixed point ceases to exist and the guard is what handles
  it. A fallback row is the exact answer bit for bit.
- ⚠️ **The argmax is discontinuous across the guard boundary, so this must not sit
  under a derivative.** Sweeping `kmax`, the mean second difference of
  `opt_psi_stem_` is 1.2x the exact solve's where the sweep stays on one method and
  **100x** where it crosses. Seating a model resets the method to `"exact"`, so
  `leaf_gradient()` always differentiates the exact solve.

`stem_c` is 2.680147 here, so the no-iteration special case `TF24_beta2 = 1/stem_c`
sits at **0.373**, not the 0.917 the header used to quote. That figure came from a
stem curve built with *E. saligna*'s measured P12 in the P50 slot, and with it goes
the claim that `beta2 = 1/stem_c` was biologically preferred: reaching `beta2 ~ 1`
that way needs `stem_c ~ 1`, a non-threshold vulnerability curve the eucalypt
literature does not support. The case is kept as a valid special case at an
unmotivated parameter value.

## Cowan-Farquhar can have a soil-moisture shutdown

`CF77` is the one cost curve here whose **price** of water does not respond to
soil state. That is the model, not an oversight: `cost_deriv<CF77>` discards
`psi_upstream` and `lambda_for<CF77>` returns the bare constant, so soil water
reaches the optimum only through the feasible bracket and the cost's *value*. The
consequence is that it has no soil-moisture shutdown, which is the one behaviour a
Medlyn-style model is usually reached for.

A new field supplies one, dividing the price by the Medlyn (2011) beta factor:

```r
l$CF77_lambda_ <- 1.5e5
l$CF77_soil_beta_ <- TRUE       # default FALSE
l$theta <- 0.35                 # with l$theta_fc, l$theta_w
l$set_model("CF77", "stem"); l$optimise()
```

`price = CF77_lambda_ / beta`, `beta = (theta - theta_w)/(theta_fc - theta_w)`.
Drying the soil then closes stomata monotonically and shuts the leaf as `theta`
approaches the wilting point; with the option off the same sweep does nothing.

**A field rather than an eighth `CostCurve`.** The enum is append-only with its
ordering propagating into `gradient.hpp`'s parameter indices, `RcppR6_classes.yml`,
R's mirror enumeration and a 4608-row golden file. A field defaulting to the
current behaviour avoids all of it: **all 4608 golden rows, 576 operating points
and 544 primitives are bit-identical**, and `cf77_price()` returns the member
untouched on the off path rather than dividing it by a computed 1.0 — do not
"simplify" that into one expression.

Two decisions worth knowing when reading output:

- **beta is clamped at 1**, so `theta` above field capacity does not price water
  *below* `CF77_lambda_`. That makes the prescribed value a floor rather than
  something the wet end can undercut, and it means that for any
  `theta >= theta_fc` this curve is plain CF77 exactly, not approximately.
- **beta <= 0 is refused rather than optimised.** At or below the wilting point the
  price is infinite or negative, and a negative price pays the leaf to transpire —
  which the optimiser takes to the wet bound and reports as an operating point.
  The check fires only when the option is on.

## `theta`, `theta_fc` and `theta_w` are now settable

⚠️ **Set these, not the trailing-underscore ones.** `theta_`, `theta_fc_` and
`theta_w_` are working copies that `set_physiology` overwrites from the inputs on
every driver set, so a bare `l$theta_ <- 0.25` was silently undone by the next
`set_drivers()` — the same shape as #96, and it would have made the soil-beta
option run every point of a sweep at the default. The three inputs survive
re-driving.

# phylloptim 0.6.0

## Seven optimality models, one optimiser

The package now solves seven stomatal optimality models, and they share a single
templated body rather than one function each. Two things specify a model: a **cost
curve** and a **benefit link**.

| `CostCurve` | objective | link |
|---|---|---|
| `TF24` | `A - C(psi)` | identity |
| `CF77` | `A - lambda*E` | identity |
| `JS22` | `A - gamma*dpsi^2` | identity |
| `CMax` | `A - (a*psibar + b)*dpsi` | identity |
| `ProfitMax` | normalised gain-risk | scaled by `|A|max` |
| `SOX` | `A * g(psi)` | log |
| `JW26` | `A * (1 - psi/psi_crit)` | log |

`JS22`, `CMax`, `SOX` and `JW26` are new. The unification is what makes the
product objectives fit: maximising `A*g` and maximising `log A + log g` share an
argmax, so every model is `max h(A(psi)) - C(psi)` and one derivative serves all
seven --

    d/dpsi [h(A) - C] = h'(A) * dA/dpsi - dC/dpsi

with `h' = 1`, `1/A` or `1/|A|max`. Adding a curve is a row in three dispatch
tables, not a new solver. Energy balance with a non-identity link is refused
rather than approximated.

All 4608 rows of `psi_stem_optima.tsv` are bit-identical across the unification:
the identity arm's derivative is textually unchanged.

## One entry point: `set_model()` then `optimise()`

**BREAKING.** Choosing a model was ten functions: seven per-curve aliases, two
dispatchers taking a curve INDEX, and the TF24 shorthand. It is now

```r
l$set_model("SOX", "stem")     # or ("TF24", "collar"), the default
l$CF77_lambda_ <- 1.5e5        # whatever constants that model needs
l$optimise()                   # takes nothing
l$model_curve(); l$model_route()
```

The model is configuration rather than a call argument, because its **parameters
already were**: `CF77_lambda_` is a field, and `JS22_gamma`, `CMax_a` and `CMax_b`
are traits. Naming the curve at the call site while its constants live on the object
configures half a model in each place. Unknown names are refused with the valid set
listed, so no caller handles an index that silently means a different model when the
enumeration grows.

Removed entirely: `optimise_psi_stem_TF/CF77/JS22/CMax/SOX/JW26/ProfitMax`,
`optimise_psi_stem_by()` and `find_root_collar_psi_by()` — from the R surface first,
and then from C++ too (see *One curve dispatcher* below). A C++ caller that wants
the compile-time form calls `optimise_psi_stem_single<K>()` or
`find_root_collar_psi_for<K>()` directly. `find_root_collar_psi()` is kept: it is
`optimise()` under the default model, and it is what plant's own sources spell.

Defaults are `"TF24"` and `"collar"`, so a caller that never sets a model gets the
production path — and all four golden baselines are bit-identical, which is that
statement checked rather than asserted.

## One curve dispatcher, not five switches

Adding a cost curve meant an arm in six places. Four of them were runtime `switch`es
over `CostCurve` with the same seven arms each — solve-collar, solve-stem,
evaluate-stem, dprofit-stem — and they are now one template, `Leaf::with_curve()`,
which hands a body a compile-time curve tag. Adding a curve touches two places: an
arm there and an arm in `curve_name()`. Neither has a `default:`, deliberately, so
`-Werror=switch` still refuses a forgotten curve.

Deleted with them: the seven `optimise_psi_stem_<curve>()` aliases, each a single
call to `optimise_psi_stem_single<K>()`, none bound to R and none called by the
model. Their per-curve documentation — JS22's always-interior wet end, CMax's
`5.930e-06` step-in tell, SOX's product-versus-log derivation — is one table above
`optimise_psi_stem_single` instead of one paragraph per alias.

⚠️ **plant must move two lines when it bumps its `phylloptim` pin.** Its
`inst/RcppR6_classes.yml` binds `optimise_psi_stem_TF` and `optimise_psi_stem_Sperry`;
the second has not existed here for some time, so plant already cannot build against
`master`, and the bump is already a coordinated edit. Both become `set_model` +
`optimise`.

**BREAKING: no method takes a curve index any more.** `evaluate_psi_stem_by(curve,
psi)` and `dprofit_dpsi_stem_by(curve, psi)` are `evaluate_psi_stem_at(psi)` and
`dprofit_dpsi_stem_checked(psi)`, reading the model `set_model()` seated — the same
rule `optimise()` follows. `vignette("the-models")` had to define its own
`curve_index()` helper to use the old pair, which is the tell that the index was
still API. Both refuse on the collar route rather than answering about a model the
leaf is not seated on, which is a distinction an index could not make.

`leaf_solve()` takes `model =` now, so the one-call R surface reaches all seven
curves rather than only TF24-collar. `leaf_gradient()` has taken it since 0.6.0.

## util::maximise_over_closed_interval is gone

It was `maximise_over_closed_interval_foc` without the root-find. Both routes moved
onto the `_foc` form, leaving it called only by its own test — a second solver in the
file whose guide says not to add one. Its test cases move onto `_foc`, plus one they
could not make: that the returned interior point is *stationary* rather than merely
inside a narrow bracket.

`Leaf::maximise_profit_over_collar` and the two `dprofit` bodies stay as they are,
and now say why at the definition: the first is this algorithm at `n = 0` plus a
five-way `OperatingPointKind` classification that 42 of 240 golden rows read; the
second pair's Identity arm is textually pinned because four golden baselines compare
at the last bit. The compensation branch's `dR_d/dT` — genuinely identical in both —
is now shared.

## Documentation the 0.6.0 changes falsified

`R CMD check` reports **Status: OK**, and `man/` is legible for the first time.

- **`Roxygen: list(markdown = TRUE)` was missing from `DESCRIPTION`**, so every
  `[fn()]` link, `**bold**` and backtick across all of `man/` rendered as literal
  text — `man/leaf_gradient.Rd` contained zero `\code{}` and zero `\link{}`. Nothing
  objected, because generated Rd files are not read and `R CMD check` does not mind a
  `[foo()]` that is only prose. ⚠️ Turning it on means `\%` in a roxygen block
  becomes a literal backslash-percent, which comments out the rest of the Rd line;
  every `\%` is now a bare `%`.
- **`man/Leaf.Rd` documented two R methods that do not exist**, `$optimise_psi_stem_
  ProfitMax()` and `$optimise_psi_stem_TF()`, both removed from the R surface in
  0.6.0. It documents `$set_model()` / `$optimise()` now.
- **Two README examples errored**, both from the `stem_b` → `stem_P50` rename: one
  passed `stem_b` to `leaf_traits()` and one named it in `pars`. `README.md` is plain
  markdown, so nothing in the repo could see them —
  `tests/testthat/test-readme.R` parses and evaluates every ```` ```r ```` block now,
  and catching exactly that regression is what it was checked against.
- **Counts, everywhere.** Fifteen traits, not fourteen; `n_pars` 18, not 16;
  thirteen reported outputs, not twelve; a fifteen-argument constructor, not
  seventeen. `vignette("cpp-interface")` told a C++ consumer to write 16 doubles into
  an 18-element array — the exact under-fill that returns `nan` with no diagnostic.
  `vignette("phylloptim")` computes its counts with inline R now.
- **`vignette("the-models")` described the pre-merge solver**, arguing that the stem
  routes "still scan" and that closing the gap was future work. 0.6.0 closed it. The
  timing table it quoted was pre-merge too: `TF24` 8.69 → **2.66** µs/call, `CF77`
  7.58 → **2.73**, `ProfitMax` 53.70 → **5.24**.
- **`ci_abs_tol` reaches neither family of optimality solver.** The vignette said it
  reached the stem routes; it is read in exactly one place, `solve_medlyn_ci_
  numerical`.
- **Leaf construction is quoted as a ratio now.** Five places said "204 µs" beside
  four different solve counts (33, 55, 70, 73), and its share of a one-parameter
  gradient had gone 40% → **61%** as the solvers got cheaper. Re-measured in one
  process: **~180× a trivial `.Call`, or 45 solves**.
- Dead references removed: `cost_curve_has_derivative()` (deleted in 0.6.0) from two
  `@seealso`; `profitmax_scan_n_` (a member that no longer exists);
  `find_root_collar_psi_by`, whose declaration outlived its definition.
- **`R/gradient.R` had an unreachable error message** telling callers to
  `$optimise_psi_stem_<curve>()`. `.GRADIENT_VERIFIED_LINKS` had come to list every
  route, so the predicate behind it was a constant. Both gone; the model registry is
  derived from `cost_curve_names()`, so it is one list rather than three.

## The collar route's zero-resistance limit converges

A collapsed collar bracket used to DERIVE `psi_stem` from continuity and stop. That
threw away the freedom that survives: feasibility pins the **collar**, not the stem,
so the leaf can still choose how far to let `psi_stem` fall. Measured at `psi_soil`
1.0, PPFD 900: deriving gave `psi_stem` 1.6465 and profit 8.255 where optimising
gives 3.0935 and 13.389 — 38% of the objective.

That branch now optimises `psi_stem` with the upstream potential pinned at the
determined collar, through the same solver everything else uses. The limit is clean
and monotone — the gap to the stem route's answer runs 9.8e-03, 3.5e-05, 3.5e-07,
7.0e-09 as the series resistance goes 1e0 to 1e-6 — where it was 38% off and not even
monotone.

⚠️ **This is not a change of coordinate, and that is what keeps it free.** Optimising
`psi_stem` everywhere would need a nested root-find per objective evaluation, because
the supply maps `psi_collar` to `E` and any other coordinate has to invert it. It is
also unnecessary: `psi_collar` and `psi_stem` are monotone in each other, so
`dJ/dpsi_stem` and `dJ/dpsi_collar` differ by a strictly positive factor and share
their roots. The coordinate was never wrong; only the bracket collapsed. Interleaved
timings are flat (collar 3.06 against 3.05 us/solve) and every golden baseline is
bit-identical, because the branch fires only where the solve previously gave up.

## One solver for every model, and the scan is now an argument

There used to be two optimisers doing the same job with different numerics. The
collar route (`find_root_collar_psi`, what `plant` calls) root-found `dJ/dpsi == 0`
but ran **TF24 only**; the `optimise_psi_stem_*` routes ran all seven curves but
refined on **bracket width**, which has no stationarity guarantee — nothing in Brent
references the derivative. Each had exactly the half the other lacked.

Both now go through one function: evaluate both endpoints, optionally scan for the
basin, then root-find the first-order condition inside the winning cell. The collar
objective is templated on the cost curve, so six of the seven solve on the production
multi-layer topology for the first time (ProfitMax is refused there — `|A|max` has no
multi-layer definition yet, and it says so rather than returning a plausible number).

**Faster and more accurate at once**, because scanning was buying less than it cost:

| | before | after |
|---|---|---|
| `psi_stem:TF` | 9.30 us/call | **2.88 us/call** |
| `psi_stem:CF77` | 8.23 us/call | **2.94 us/call** |
| `psi_stem:ProfitMax` | 58.2 us/call | **5.64 us/call** |
| collar solve | 3.29 us/solve | 3.30 us/solve (interleaved) |
| median `\|dJ/dpsi\|` at psi* | 4.4e-05 | **1.2e-15** |

`operating_points.tsv` is bit-identical, so the production path did not move.
`psi_stem_optima.tsv` was regenerated: the argmax is now stationary rather than
grid-resolved, worst relative movement 1.5e-05 with profit changing at 1e-10.

**Whether to scan is measured, not assumed.** Over a 1728-row sweep (8 air
temperatures x 4 gate combinations x 6 soil potentials x 3 deficits x 3 light levels)
the only objectives carrying two prominent interior basins are `TF24` (21 rows) and
`JS22` (66) — and **every one of those rows has the energy balance on**. The collar
objective has none in 432 rows over the same range. So the scan is skipped with the
gate off, which is what makes the routes ~3x faster. The valleys are 7.0e-01 and
3.5e+00 deep, so this is a cheap scan to skip and an expensive one to skip wrongly:
re-measure before widening the predicate.

### `|A|max` is solved, not scanned

ProfitMax normalises by the largest assimilation on the supply stream, found until
now by a 500-point scan. Measured, that maximum sits at the dry bound on 1318 of 1320
driver rows and interior on 2, beating the bound by at most 1.1e-03 — so it is
endpoints plus a root-find on `dA/dpsi == 0` like everything else. `dA/dpsi` needs no
new algebra: for an identity link it is `dJ/dpsi + dC/dpsi`.

That is the 10x on ProfitMax above, and it also fixes its gradient. A scan argmax is
piecewise constant in the parameters, so differencing through it was step-dependent:
0.014 at `h = 1e-06` against 0.058 at 1e-02. It is now flat across six decades, so
**one finite-difference step (1e-06) serves every route** and the per-route split is
gone. ProfitMax's partial-versus-total gradient distinction remains — the composite
holds `|A|max` fixed — but it is a modelling choice now rather than a staircase.

### Removed

**BREAKING.** All dead or duplicated by the above:

* `cost_curve_has_derivative()` — returned `TRUE` for all seven curves and nothing
  read it; every curve has a derivative through its benefit link
* `leaf_behaviour_fingerprint()`, with `tools/fingerprint.R` and its test — no
  consumer, and it only ever covered two of the four recorded baselines. The
  bit-exact golden comparison is the guard
* `profitmax_scan_n_` — an R-exposed field controlling a scan that no longer exists
* `GSS_tol_abs` is now documented for what it actually reaches: two places, both on
  the collar route. It never reached the `optimise_psi_stem_*` refinement, whose
  tolerance is not settable
* `Imports: tools`

### `leaf_supply_single()` is `leaf_supply_singlelayer()`

**BREAKING**, and purely for symmetry with `leaf_supply_multilayer()`.

⚠️ Worth knowing while renaming callers: the single-layer *supply path* carries a
series resistance, but the `optimise_psi_stem_*` routes pin the upstream potential at
`psi_soil` and **ignore it**. "Single layer" and "no root resistance" are not the same
thing, and a single-layer plant with a real root resistance is not yet representable
in the stem formulation.

## Both vulnerability curves are parameterised on (P50, c)

**BREAKING.** The stem curve's traits are now `stem_P50` and `stem_c`, and the
root curve's `root_P50` and `root_c`. `stem_b`, `psi_crit`, `root_b` and
`root_psi_crit` are **derived**: still readable on a `Leaf`, no longer settable
and no longer arguments to `leaf_traits()`, `set_traits()` or the constructor.

    b        = P50 / (ln 2)^(1/c)
    psi_crit = b * ln(1/0.05)^(1/c)          i.e. P95 of the same curve

`perturb_stem_b()` becomes `perturb_stem_P50()`. See *Model-scoped parameter
names* below for the parameter set as a whole.

Why: `psi_crit` was a free trait describing a curve it was not derived from, and
Sperry's reference conductivity `k_crit = 0.05 * kmax` used the same 0.05 from a
separate hard-coded constant. The three could disagree, and a consistency check
existed to police them. They are now one number, so the check is gone -- along
with the domain check, which becomes vacuous (a derived P95 is inside P99 for
every c > 0).

**The reference values moved**, and were regenerated. The defaults were round in
`stem_b` and `psi_crit` rather than in `P50`, so `psi_crit` was P95 rounded to
seven figures and `f(psi_crit)` was 0.049999968739194864 rather than 0.05.
Adopting a round `P50 = 3.4` shifts `stem_b` by 3.1e-08 and `psi_crit` by 4.7e-08,
which propagates to a worst 4.4e-05 on assimilation over the 576-point grid and
2.2e-05 on profit over the optimiser grid (2304 rows when this was measured; the
grid ships at 4608) -- below the 1e-4 at which this
project calls a difference real. No operating point changed branch. The lowest
primitive tier to move is the vulnerability curve itself, at 3.1e-08; the
arithmetic tier is unchanged.

Two behavioural consequences, both asserted:

* A shut-down leaf's hydraulic cost no longer depends on the curve's scale. It
  sits at `psi_crit`, where remaining conductivity is exactly 0.05 whatever `P50`
  is, so `dprofit/dstem_P50` there is exactly zero rather than negative.
* Stem gradients carry a term through the moving `psi_crit`. It is exactly zero at
  an interior optimum and order 1 at a pinned one, where it is 29% of the total.

## Model-scoped parameter names

**BREAKING.** A parameter that belongs to one model now carries that model's name.

| was | is |
|---|---|
| `beta2` | `TF24_beta2` |
| `cost_scale_TF24` | `TF24_cost_scale` |
| `lambda_` | `CF77_lambda_` |

New parameters follow the same rule: `JS22_gamma`, `CMax_a`, `CMax_b`. Models
without a published name take initials-plus-year.

### What the parameter set looks like against 0.5.2

Everything below landed in one release, so the deltas that matter to a caller are
against **0.5.2**, the last version published:

* the trait vector goes **14 to 15** and `gradient_par_names()` **16 to 18**
* four traits are gone, being derived now: `stem_b`, `psi_crit`, `root_b`,
  `root_psi_crit`
* five are new: `stem_P50`, `root_P50`, `JS22_gamma`, `CMax_a`, `CMax_b`
* two are renamed: `beta2` and `cost_scale_TF24`
* `CF77_lambda_` is differentiable, which is the 18th parameter

Traits are bound **positionally** across the R/C++ boundary, so a caller passing
them by position must be updated with this release. Appending is safe; reordering
differentiates the wrong thing and returns plausible numbers.

## The prescribed-lambda optimiser is removed

**BREAKING.** `optimise_psi_stem_Sperry()`, `profit_psi_stem_Sperry()` and
`hydraulic_cost_Sperry()` are deleted, with no deprecated alias. Nothing in the
plant-family tree called them.

They maximised `A - lambda*(k(psi_soil) - k(psi))` with lambda a prescribed
constant, under Sperry's name. `optimise_psi_stem_ProfitMax()` is Sperry's model,
and the point of it is that the cost scaling is **emergent**: multiplying the
normalised objective by |A|max gives exactly that constant-lambda form with
`lambda = |A|max / (k_soil - k_crit)`, a quantity that moves with the drivers.
Measured here, the implied lambda runs 9.19e4 to 3.15e5 over psi_soil 0.5 to
3 MPa, a 3.4x range. So a fixed lambda is not a variant of the model, and the two
entry points were not two models.

`$lambda_` now has one meaning: the Cowan-Farquhar marginal value of water, in
umol C (kg H2O)^-1, supplied by the caller.

## Gradients for every model, through one entry point

**BREAKING (signature).** `leaf_gradient()` and `leaf_gradient_batch()` take
`model =`, defaulting to `"collar"` -- the production path, unchanged. Naming any
of the seven curves differentiates that model's **stem** optimum instead. All
eight routes are verified; none is refused.

`leaf_gradient_batch()` remains bit-for-bit identical to the R route on every
route, which is what makes a fit free to choose its model.

⚠️ **Each route carries its own finite-difference step.** A difference through a
solve has a noise floor set by the solver being differenced: the collar route
root-finds to ~1e-15 and takes 1e-6, while the stem optimisers scan and refine to
~1e-6 of the bracket and take 1e-3. A 1e-6 step on a stem route returns
quantisation -- measured 0.1855 against a true 0.0551, and the wrong sign on one
curve. Sweep the step and read the floor of the V; do not trust a single step.

New: `cost_curve_names()` and `cost_curve_has_derivative()`.

## `$lambda_` is an input; every curve reports an emergent lambda

**BREAKING for anyone reading `$lambda_` after a ProfitMax solve.**
`optimise_psi_stem_ProfitMax()` used to overwrite `$lambda_`, so solving it and
then Cowan-Farquhar on the same leaf priced water at ProfitMax's number rather
than the caller's — silently, since both are finite and plausible. `$lambda_` is
now an input that no model code touches.

In its place, **every** cost curve now reports the marginal cost of water its
operating point implies, on one shared axis:

    $lambda_emergent = (dC/dpsi) / (dE/dpsi)      umol C (kg H2O)^-1

Same units as the `$lambda_` input, so the two are directly comparable — which is
the point. Cowan-Farquhar prices water at a constant, so its emergent lambda IS
that constant, exactly. TF24's is `lambda_TF24` at the operating point, the collar
solve's carries the series-resistance correction, and ProfitMax's is derived from
its normalised cost and scaled by |A|max to restore carbon units. Verified against
finite differences of each curve's own dC/dE: worst relative error 4.9e-07.

`$lambda_emergent` is read-only, and NA until an optimiser has run.

⚠️ ProfitMax's `|A|max / k_span` is **not** a lambda, whatever it has been called:
its units are carbon per unit CONDUCTANCE lost, because it multiplies
`k(psi_soil) - k(psi)`. It is a normaliser. `$profitmax_A_max` and
`$profitmax_k_span` are now exposed read-only so that ratio is one division away.

Both older golden files stay bit-identical and the fingerprint is unchanged. The
optimiser fixture gains a `lambda_emergent` column; no other column moves.

No numbers move on any surviving path: `operating_points.tsv` and
`primitives.tsv` are bit-identical, the behaviour fingerprint is unchanged, and
`psi_stem_optima.tsv` loses its 576 Sperry rows with no surviving row altered.
The equivalence the deleted entry point used to demonstrate is now asserted
directly as an identity, which holds to 9.5e-16 rather than the 5e-3 two argmaxes
could agree to.

## Constrained optima are named for the bound that binds

**BREAKING for C++ consumers.** `OperatingPointKind::PinnedWet` and `PinnedDry`
are now `BoundarySoil` and `BoundaryCrit`, and `operating_point_kind_name()`
returns `"boundary-soil"` and `"boundary-crit"`. No R surface: the accessor is C++
only.

The old names read backwards. `psi_stem == psi_soil` is the least-tension end of
the feasible interval, so "wet" sounds like the comfortable case -- but it is what
a leaf does under stress, when no interior point earns its water. Warming the
reference grid from 25 to 40 C takes that kind from **24 rows to 80** while the
`psi_crit` end goes from **18 to none**, so the reading inverted with temperature.
Naming each for the bound it sits on cannot invert.

`psi_stem_optima.tsv`'s `kind` column was patched in place rather than
regenerated, so the bit-exact comparison itself proves that the other thirteen columns
did not move: 32 cells, columns 1-9 and 11+ byte-identical.

## The single-layer optimisers reach a maximum at a bound (#94)

⚠️ **Results move on `optimise_psi_stem_TF` and `optimise_psi_stem_Sperry`.** Both
maximised profit with a bare `brent_fmin`, which steps in from the bounds and
follows one basin — so it could return neither endpoint nor the taller of two
humps. Both failure modes are reachable and neither is rare: over a 1728-row
driver sweep (4 gate combinations × 8 air temperatures × 6 soil potentials × 3
deficits × 3 light levels), measured against a 2001-point scan of the same
objective,

| entry point | rows short of the reference | worst profit lost | worst ψ gap |
|---|---|---|---|
| `optimise_psi_stem_TF` | **866 / 1728** | 1.70e-01 | 1.77 MPa |
| `optimise_psi_stem_Sperry` | **837 / 1728** | 4.71e-01 | 4.50 MPa |
| `optimise_psi_stem_ProfitMax` | 29 / 1728 | 2.54e-05 | 2e-04 MPa |

824 and 780 of those rows have their true optimum **at** a bound. All three now
match the reference on **every row**, exactly.

Both now route through `util::maximise_over_closed_interval`, which evaluates the
endpoints, scans `boundary_scan_n_` cells to locate the basin, and refines with
Brent inside the winning cell. `optimise_psi_stem_ProfitMax` keeps reusing
`prepare_profitmax`'s own 500-point scan, so it gains no evaluations.

**One part of this is easy to get wrong and was:** the refinement tolerance has to
scale with the CELL, not with the interval. Brent terminates on bracket width, so a
fixed `GSS_tol_abs` comparable to a cell leaves the answer at essentially the grid
point — and then adding grid points makes the result *worse*, because cells narrow
while the tolerance does not. Measured: with a fixed tolerance the Sperry shortfall
count **rose** from 64 (32 cells) to 104 (64 cells); with the cell-scaled tolerance
64 cells is exact. `optimise_psi_stem_ProfitMax` had the fixed-tolerance form and
now uses the scaled one; that is what closes its remaining 29 rows.

⚠️ **Cost, and it is not small.** `optimise_psi_stem_TF` goes **157 → 395 µs**
(+150%, interleaved ×3, reproducible to ~1%) because its objective is cheap and 65
extra evaluations dominate it. `optimise_psi_stem_Sperry` pays the same ~240 µs on a
larger base (+9%). `optimise_psi_stem_ProfitMax` is unchanged (2443 vs 2445 µs) and
so is **`find_root_collar_psi` — 77.9 µs both ways, checksum identical.** The
production path is untouched: these three entry points require
`supply_is_single_layer()` and plant does not call them.

`boundary_scan_n_` is R-settable for callers who would rather have the speed. 32
cells halves the added cost and is exact for TF24, leaving 6 Sperry rows short by
3.9e-04; 8 cells leaves 17 and 7 rows short by 3.2e-03 and 3.9e-04. Below that the
boundary half still works and the second-hump half stops being reliable.

Both golden files are bit-identical — the grid uses the collar solve, so it never
touched any of this, which is why the defect survived to be measured here.

## Published vulnerability curves convert to (P50, c)

`weibull_p50_c()` takes any two of `P50`, `P88`, `P95`, `P99`, `(px, plc)`, `S50`
and `c`, and returns the pair `leaf_traits()` wants. Also exported:
`psi_at_plc()`, `weibull_s50()` and `weibull_p50_from_b()`.

⚠️ `S50` with a non-P50 quantile has **two** solutions -- `S50(c)` has an analytic
minimum at `c = ln L`, so it is U-shaped. The upper branch is returned, with a
message; the lower branch lands below `c = 1`, where measured angiosperm stems sit
at 1.8-2.6.

Published P50 is negative and is **refused** rather than silently negated: every
potential in this package is a positive magnitude.

## `plot_model_overview()`: the whole model in one figure

New exported function, and the opening figure of `vignette("the-models")`, which
it replaces the bare profit curve with. Four zones: the profit trade-off and the
supply-equals-demand root-find computed from a live `Leaf`, then a schematic of
how the parts are wired and of the soil-to-leaf path they describe.

The wiring panel is what it exists for. One decision variable feeds two branches
that meet again at the objective; the two gates are drawn dashed because both
default off; and the arrows crossing between the columns are the interactions that
are otherwise only stated in prose — the energy balance as the cycle it is, and
the thermal cost reaching past it into the *benefit* column, which is the one gate
that moves both sides of the ledger.

Base graphics on a vector device, so the vignette embeds it as SVG and `pdf()`
gives a manuscript figure whose every line and label stays a separate object.
Shapes are aspect-corrected and box labels shrink to fit, so it holds together at
figure sizes other than the default. No model code is touched and no number moves.

# phylloptim 0.5.2

## `lambda_analytical_` is deleted (#113)

**API change, no numbers move.** The field was declared, bound to R as a settable
field, and read by nothing — no model code, no test, no solve path, no gradient. It
was write-only from R and had no effect on any number. Isaac's reading is that it
belonged to an `optimise_psi_stem_sperry_analytical` that no longer exists.

Deleted rather than populated, and the reason matters because the obvious repair is
wrong. The one λ this class derives is `profitmax_A_max_ / profitmax_k_span_`, and
`profitmax_A_max_` comes from a 500-point scan of the transpiration supply stream —
`prepare_profitmax` explains why the closed form is unavailable once assimilation is
net rather than gross. So the only candidate value was the one value that must not be
called *analytical*, while the genuinely closed-form λs, `lambda_TF24()` and
`marginal_cost_water()`, are accessors by design (hazard 5).

⚠️ **The dual-role problem this field looked like a home for is still open.**
`lambda_` is an input to `optimise_psi_stem_Sperry` and an output of
`optimise_psi_stem_ProfitMax`, so calling the second and then the first silently
optimises a derived λ rather than a prescribed one. That wants a new member with an
honest name; it is #114's third recommendation, not this change.

⚠️ **Hazard 7: plant carries the same generated glue** (`inst/RcppR6_classes.yml`,
`src/RcppR6.cpp`), and `plant::Leaf` *is* `phylloptim::Leaf`, so plant's generated
`obj_->lambda_analytical_` stops compiling against a phylloptim without the field.
plant pins `traitecoevo/phylloptim@037673b7` in `Remotes:` rather than master, so
this merge breaks nothing today — plant's half lands with its next sha bump.
# phylloptim 0.5.1

## The zero-transpiration branch inside the objective (#110, #112)

⚠️ **Results move on the energy-balance path.** `set_leaf_states_rates_from_psi_stem`
has a `psi_upstream >= psi_stem` branch for the no-flow / reversed-gradient case. It
wrote `Tleaf_` but re-derived neither the Farquhar block nor the leaf-to-air deficit
at it — while `ci_` is read off `gamma_`, `vpd_leaf_` is a reported field, and
`assim_colimited_` at the bottom of the function is evaluated against the whole
block. So the branch assembled a state out of two temperatures.

⚠️ **The other temperature was not the air baseline**, which is what separates this
from #105 rather than making it a repeat. The transpiring branch re-derives the block
per candidate *by design*, so this branch inherited whatever ψ the optimiser last
probed. The objective's value at a fixed ψ therefore depended on the **order**
candidates were visited — hazard 3's mechanism, arrived at from the model rather than
from the solver.

Measured at Tair 30 °C, PPFD 900, ψ_soil 1 MPa, D_air 1.5 kPa, where the branch's own
leaf temperature is 39.33 °C:

| reported at the zero-E point | on a cold object | after one transpiring candidate | correct |
|---|---|---|---|
| `assim_colimited_` | −1.9888 | −2.5355 | **−3.0978** |
| `vpd_leaf_` (kPa) | 1.5000 | 2.7368 | **4.3714** |
| `ci_` (Pa) | 5.5700 | 6.9492 | **8.7173** |

Each branch now derives the block and the deficit at the temperature it reports,
before `ci_` is read out of them — the same repair #111 made at the exits *outside*
the solve, and in the same order, so `Tleaf_` moves to the top of the branch.

**Where it changes an answer, not just a report.** `prepare_profitmax` scans from
ψ_soil, which *is* this branch, and records the assimilation there as a profit
candidate — so `optimise_psi_stem_ProfitMax`'s argmax moves. At Tair 40 °C with the
thermal cost on, ψ_soil 0.5 and 1.0 MPa move from full closure to an interior optimum
(1.605 and 1.760 MPa, profit −1.2309 against −1.2825); at Tair 50 and 55 °C the
argmax moves the other way, from barely open (0.5004) to exactly shut. The
prescribed-collar path moves too, where a caller's target clamps onto the wet bound:
profit −2.4329 → −2.5415 at Tair 25 °C.

**`find_root_collar_psi` itself is unmoved** at every PM-path driver tried, because
`dprofit_at_collar_psi` has its own reversed-gradient exit and the collar it returns
is stepped inside the wet bound. That is a statement about the drivers tried, not a
proof: the golden-section fallback and the two-pin failure case both evaluate profit
at that bound.

Cost, interleaved ×4 (hazard 5, since this adds transcendentals to a branch of a
function the collar solve runs ~10³ times per solve): **no signal**. The branch itself
is ≤1% slower per (transpiring + zero-E) pair, within the round-to-round scatter; the
PM collar solve, the ProfitMax optimiser and the gate-off solve are all
indistinguishable. Gate-off is unchanged by construction — the recompute is inside the
gate, and `set_leaf_vpd` returns `atm_vpd_` exactly there. Both golden files are
bit-identical, which here is a statement about the grid (it runs gate-off) rather
than evidence.

# phylloptim 0.5.0

## Documentation: claims that #55 and #93 outlived

No behaviour change. Four places still described the model as it was before two
merged fixes, and one of them was user-facing advice:

- `vignettes/fitting.Rmd` told readers that a bare `l$vcmax_25 <- x` followed by
  re-driving takes a temperature-cache **hit** and silently reports numbers from the
  first `vcmax` it ever saw. #55 widened the cache key to every input of the
  temperature block, so that is no longer true. The advice — use `set_traits()` — is
  unchanged; the reasons are now the splines, the solved operating point and the #25
  checks, which no re-driving repairs.
- `inst/RcppR6_classes.yml` carried the same claim, contradicting the correct
  statement 150 lines above it in the same file.
- `tests/cpp/test_leaf.cpp` claimed a sub-assertion isolated `set_traits`' cache
  invalidation. It cannot any more: with `vcmax_25` in the key the cache misses
  regardless, so the trap is closed twice over and the test cannot say which
  mechanism carried it. The two halves are covered separately, and it now says so.
- `inst/include/phylloptim/closed_form.hpp` divides by `atm_vpd_` throughout, which
  #93 replaced with the leaf-to-air `vpd_leaf_` in the live model. Recorded rather
  than repaired: it is dormant, PM-untested code, and the substitution may not be a
  rename — `vpd_leaf_` depends on Tleaf, which depends on the E these expressions
  solve for.

## `lambda_` is a caller input, so nothing resets it (#96)

`setup_clean_leaf()` cleared `lambda_`, so whether a prescribed Sperry marginal water
cost survived depended on which of two interchangeable-looking calls came next:

```r
l$lambda_ <- 30; set_drivers(l, ...)   # lambda_ is 30
l$lambda_ <- 30; l$set_traits(...)     # lambda_ was NaN -- now 30
```

Neither call touches it now. It is read by `profit_psi_stem_Sperry` and set by no
model code on that path, so wiping it on the caller's behalf was never right — and
clearing it in `set_physiology` too, the other way to make the rule consistent, would
have made a prescribed λ unusable, since the multi-layer path requires
`set_physiology` before every solve. `lambda_analytical_` moves with it, on the same
argument.

Both carry their NA default at the declaration instead, so a fresh `Leaf` still reads
NA and `optimise_psi_stem_Sperry`'s guard against an unset λ still fires. Derived
state and solved outputs are cleared exactly as before.

⚠️ **`lambda_` is now dual-role**, and this change does not fix that: #93's
`optimise_psi_stem_ProfitMax` *writes* it, so calling ProfitMax and then Sperry
optimises ProfitMax's derived λ rather than a prescribed one — finite, plausible, and
past the guard. Documented at both sites and pinned by a test; #114 is the design
question, #113 a possible resolution.

No numbers move: `lambda_` reaches only `profit_psi_stem_Sperry`, which no golden path
calls. Both golden files are bit-identical.

## A shut-down leaf reports one temperature, not two (#105)

⚠️ **Results move on the energy-balance path.** Both zero-transpiration exits of the
collar solve — hydraulic shut-down and shade death — reported a leaf temperature and a
respiration rate that belonged to *different* temperatures. They do not go through
`set_leaf_states_rates_from_psi_stem`, so they inherited the Tair baseline
`set_physiology` derives and then described a state whose transpiration is zero. Both
now re-derive the temperature block at the temperature they report, before profit is
formed.

Measured against `origin/master` at `atm_vpd = 2`, `atm_kpa = 101.3`, Tair 30 °C:

| exit | Tleaf | | `R_d_` | `profit_` | `ci_` | `vpd_leaf_` |
|---|---:|---|---:|---:|---:|---:|
| hydraulic shut-down | 39.325 | before | 1.9888 | −8.9334 | 5.5700 | 2.0000 |
| (`psi_soil = 6.0`, PPFD 900) | | **after** | **3.0978** | **−10.0423** | **8.7173** | **4.8714** |
| shade death | 28.957 | before | 1.9888 | −2.0204 | **NaN** | 2.0000 |
| (`psi_soil = 1.0`, PPFD 1) | | **after** | **1.8679** | **−1.8995** | **5.2889** | **1.7528** |

**The two exits move in opposite directions.** `Rn` is proportional to PPFD with a
fixed longwave offset subtracted, so a hydraulically shut leaf in full sun runs hotter
than the air while a shade-dead one runs cooler.

**`vpd_leaf_` is on the same footing**, and is the half that #93 created. The
leaf-to-air deficit is derived from a leaf temperature and is read by `g1_eff()`, but
these exits left it at the Tair deficit `set_physiology` seated — so a shut-down leaf
reported a 39.3 °C temperature beside the deficit belonging to 30 °C, 2.44× too small.
#93 already makes exactly this argument at its own zero-transpiration branch; the two
exits outside the solve were simply not on its list.

The shade-death exit also never wrote `ci_` at all — it reported the previous solve's
internal CO2, or `NaN` on a cold object. It is now the compensation point, as at the
other exits.

**Nothing moves with the gate off**, which is the default and what both golden files
cover: there `Tleaf_` *is* the `leaf_temp` driver and the block was already derived at
it. `operating_points.tsv` and `primitives.tsv` are bit-identical, and
`leaf_behaviour_fingerprint()` therefore does **not** move — the golden grids run
gate-off, which is the incompleteness its own documentation warns about. The version
bump is what signals this one.

**The feasibility test stays at the Tair baseline, deliberately.** `assim_max_ < 0`
asks whether any *open* operating point is worth taking, and an open leaf transpires
and is therefore cooler than these exits' leaf. Re-taking that test at the
zero-transpiration temperature would judge the coolest option at the hottest
temperature. So the branch is chosen on the transpiring baseline and the state it
reports is self-consistent at the temperature it describes.

# phylloptim 0.4.0

## Infeasibility is a condition class, not a message to match on (#57)

A failed solve can now be caught by class:

```r
tryCatch(leaf_solve(...), phylloptim_infeasible = function(e) e$code)
```

`leaf_infeasible_codes()` documents the seven codes and what each means.
`leaf_solve()` and `leaf_gradient_batch()` classify without the caller doing
anything; `with_phylloptim_conditions()` is for the caller driving a `Leaf` directly,
since `$find_root_collar_psi()` is generated glue that cannot signal on its own.

**Why not message text.** That was the only option before, and it was never stable:
#65 and #79 each rewrote the out-of-domain wording and #92 added a third variant, so
any caller matching on prose was silently broken by an improvement to an error
message.

**Two layers, because it cannot be one.** C++ tags the exit with a stable token —
`[phylloptim:infeasible:collar_bracket] …`, via a new `util::stop_infeasible` and an
`infeasible_error` type — and hand-written R reads the token and re-signals. The join
could not live at the boundary: that is RcppR6-generated and must not be hand-edited,
and `util.hpp` must not reach for Rcpp. A C++ consumer such as plant catches the type
directly.

⚠️ **Seventeen of about fifty `stop` sites are classified, and the conservative
direction was taken everywhere it was arguable.** Every input-validation failure stays
an ordinary error, because classifying one as infeasible would let it be swallowed by
the very `tryCatch` this enables — a fit would then report a plausible likelihood over
whichever rows survived, which is exactly what #39 rejected silent NA rows to prevent.
`?leaf_infeasible_codes` lists the three genuinely arguable sites left ordinary
(`psi_crit` past P99, an explicit `method = "ift"` at a point where the IFT does not
hold, and the trait sign checks) with the reasoning, so the calls can be overruled
deliberately.

⚠️ **The classification cannot drift from its documentation**: a test scans the
*installed* headers for `stop_infeasible("…")` and requires the set to equal
`names(leaf_infeasible_codes())` exactly.

**Still not done, deliberately: `on_error = "na"` with a per-row status column.** The
issue is explicit that it must not be built before the classification is reviewed, or
the opt-in returns `NA` for programming errors too. `leaf_gradient_batch()` already
reports per-row `status`, which is that idea where it was already safe.

⚠️ **`leaf_gradient()` is not routed through the wrapper** — it depends on `missing()`
for two arguments, so delegating to an inner implementation would change argument
matching. Wrap the call, or use `leaf_gradient_batch()`, which a fit should be on
anyway. The gap is documented at `?with_phylloptim_conditions`.

The wrapper is **one `tryCatch` per call, outside the row loop**, so per-row cost is
unchanged.

## A primitive-level baseline, and three dead harnesses removed (#64)

`tests/cpp/golden/primitives.tsv` records **544 values from the functions the solve
calls**, in five call-tree tiers — arithmetic, vulnerability, assimilation, spline,
iterative — read by `tests/cpp/test_primitives.cpp`. Bit-exact on macOS/arm64, with
per-**tier** tolerances elsewhere.

**Why a second baseline.** `operating_points.tsv` detects a difference and cannot
attribute one: the nested solvers amplify a last-bit change up to the loosest
tolerance, so by the time it reaches a reported output it points nowhere. #92 proved
that twice over — the knot grid had dropped its last knot for years and the golden
file recorded the consequences as *correct*, and then the fix's 3496 moved cells could
not be split into "knot displacement" versus "added knot" without a standalone
primitive program, written in a scratch directory and thrown away. Both times.

Demonstrated rather than asserted: a 1e-12 perturbation to
`proportion_of_conductivity` leaves the arithmetic and assimilation tiers at
**exactly 0**, moves vulnerability at 1e-12 and spline at 4.2e-12. **The lowest tier
that moved is the cause.**

The knot grid is now pinned per (b, c) pair — `knot_grid_last_psi` sits beside
`vulnerability_psi_max` for nine pairs, so #92's defect is a row in a file rather
than a test nobody thought to write.

**Deleted:** `tests/validate/compare_with_plant.R`, `compare_primitives.R` and
`primitives.cpp`. The first two compared against plant's own `Leaf`, and plant has no
independent copy of this model any more — it consumes these headers — so that
comparison had become this package against itself. What they established is PLAN
decision 1. `scm_regression.R`, `scm_compare.R` and `tsv_to_hex.c` are kept; the last
is live, in `test-golden.R`'s regeneration recipe.

⚠️ **The enforcement is the deliverable.** Those files rotted through two interface
changes because **no build touched them**, so a revived harness in no build would be
the same defect with a newer date. `test_primitives` is named in `tests/cpp/Makefile`,
`CMakeLists.txt` (as the `leaf_primitives` ctest), `tests/cpp.R` and `.Rbuildignore`,
and CI reaches it through `make` on all three platforms. Adding a program under
`tests/cpp/` now means editing four files, and each of them says so.

⚠️ **The cross-platform tolerances are mechanism-derived ceilings, not measured worst
cases**, and they are labelled as such — the guide's cross-platform table has been
wrong three times from a reading written down as a fact. The real figures come from
the summary line the tolerant mode prints on every run, pass or fail.

## The leaf's own temperature is an output

`Tleaf` joins the reported operating point, as a thirteenth column on
`leaf_solve()` / `operating_point()` and as a `$Tleaf_` field.

**On the energy-balance path it was unreachable.** Leaf temperature is solved there
per operating point, used to re-derive the whole Farquhar block, and was then
discarded — and `leaf_temp` is no substitute, because `set_physiology` reinterprets
that driver as **air** temperature on that path. So the one quantity the
Penman-Monteith path exists to produce was the one thing a caller could not read;
`test_energy_balance_path_runs` had to derive it again from an output to say anything
about it. Measured at `PPFD = 900`, `atm_vpd = 2`, air 30 °C, the leaf runs 8.4–9.3 K
hotter than the air and warms by ~1 K along a drydown as the stomata close.

Off that path `Tleaf` equals the `leaf_temp` driver, deliberately rather than being
NA: a column that is sometimes a temperature and sometimes missing cannot be plotted
against anything.

**Appended, not inserted.** These names are positions, so `profit` stays at index 7
and a saved output does not shift.

⚠️ **Written by every exit from the solve, which is hazard 8 and is the part with
teeth.** The two shut-down exits do not go through
`set_leaf_states_rates_from_psi_stem`, and a leaf that transpires nothing is the
*hottest* one — so leaving them to inherit would have reported the previous
individual's leaf temperature at exactly the operating points where the answer is
most extreme. plant drives every individual in a patch through one persistent `Leaf`.

⚠️ **Exposing it revealed an inconsistency, filed as #105 rather than fixed here.** At
a PM shut-down the reported `A = -R_d` is still computed at *air* temperature while
`Tleaf` is 39.3 °C: respiration 35.8% low, profit −8.93 against −10.04. This release
reports the physically correct E = 0 temperature and leaves the arithmetic alone —
reporting air temperature instead would have made the pair agree by reporting a
temperature the leaf is not at, which is what kept it invisible.

## The `stem_b` shortcut keeps its saving through the batch (#74)

`perturb_stem_b()` exists so that a gradient in `stem_b` needs no vulnerability-spline
rebuild — 24.5× on that perturbation. Through `leaf_gradient_batch()` it was worth
3.4× less than that, because the shortcut leaves `stem_b` displaced from
`stem_b_spline_` and the restore at the end of each observation goes through
`set_traits()`, whose third rebuild clause exists precisely to catch that state. So
every observation paid for one rebuild that no perturbation had asked for.

`gradient::apply()` now undoes the displacement with the shortcut before restoring,
after which `set_traits()` has nothing to rebuild. Measured per observation, three
observations of `d/dstem_b`, interleaved:

| | µs/observation | stem-curve rebuilds/observation |
|---|---:|---:|
| before | 26.5 | 1.00 |
| after | 7.8 | **0.00** |
| `vcmax_25`, the no-spline floor | 7.9 | 0.00 |

**Bit-identical, and now measured rather than argued.** `gradient_golden.tsv` compares
at a worst relative difference of exactly 0, and the C++ golden grid is unmoved. The
reason is that `perturb_stem_b()` writes `stem_b` and nothing else — the splines are
the ones built at `stem_b_spline_`, not a rescaled copy, so putting `stem_b` back
returns the object a rebuild would have produced.

⚠️ **`stem_c` is unchanged and deliberately so.** It has no homogeneity identity, its
3.00 rebuilds per observation are genuine, and only the restore is removable there.

⚠️ **The two figures are different measurements**, which is why this was quoted wrong:
24.5× is per *perturbation* in C++, and what a fit pays is per *observation*. Both are
now printed side by side by `bench_gradient`, and `Leaf::stem_curve_builds_` — a
counter, not state — makes "no rebuild" an integer the test suite asserts rather than
a stopwatch reading.

## Conductance to water, and a settable diffusion ratio (#50, #56)

`gs_H2O` reports stomatal conductance to **water vapour**, alongside the `stom_cond_CO2_`
the model solves for. Every data source records conductance to water, as does the g1
literature, so the conversion was being left to user code — where, as the calibration
study put it, a wrong factor biases `gs` by a constant and the fit launders it into
`kmax` rather than revealing it. An accessor, not stored state (hazard 5).

`H2O_CO2_stom_diff_ratio_` is now a settable field, **defaulting to 1.67 so nothing
moves**: the golden file is bit-identical. It encodes a convention rather than a
measured property of this leaf — Medlyn et al. (2011) and the `g1` literature use 1.6,
from the binary diffusivities — and it was a compile-time constant, so a caller who knew
their comparison wanted 1.6 could not say so.

⚠️ **#50's "~2.2% offset in every reported g1_eff" is wrong, and the correction matters
in both directions.** `g1_eff` does not contain the ratio at all — it is
`chi*sqrt(D)/(1-chi)` — so the effect arrives entirely through `ci` moving in the solve.
Measured at `psi_soil = 2, PPFD = 900`, taking 1.67 to 1.60 (a 4.19% change in the
ratio):

| quantity | 1.67 | 1.60 | rel diff |
|---|---|---|---|
| `g1_eff` | 0.502639 | 0.521100 | **3.67%** |
| `gc` | 0.0192260 | 0.0204187 | 6.20% |
| `gs_H2O` | 0.0321075 | 0.0326700 | 1.75% |
| `A` | 5.60102 | 5.89174 | 5.19% |

So the offset is **1.7x larger than the issue states**, and it is not a unit factor that
could be corrected outside the package — which strengthens the case for exposing it. Note
`gs_H2O` moves *least*, because the ratio partly cancels there.

⚠️ It reaches the solve, not only a reported output (the conductance conversion and
`dprofit`'s `gc_const`), so a solved leaf must be re-solved after changing it.

**#56's second half is documented rather than fixed**, which is what that issue asks for
regardless: `?leaf_supply_singlelayer` now states that `leaf_specific_conductance_max` is
kg-based while `series_resistance()`'s resistance is mol-based, so a caller
parameterising the whole path carries a `molar_mass_h2o` factor between two quantities
presented as two ends of one series. The calibration study recorded dropping that 0.018
as its original error, worth three orders of magnitude. Unifying the bases would move
results *and* change a published argument's units, so it stays a decision.

## The gas constant is the SI value (#51 audit) — MOVES RESULTS, AT 40 C ONLY

`gas_constant` was `8.314`. Since the 2019 SI redefinition R = N_A k_B is **exact**, at
8.314462618153240 J mol^-1 K^-1, so there was no reading on which 8.314 was right — the
only question was whether a 5.56e-05 correction justified a results change.

It does, because the amplification is not 5.56e-05. This constant appears only in
`arrh_curve` and `peak_arrh_curve`, always as `Ea/(R*T)` with `Ea/RT` of order 24 at the
defaults, so a relative change eps in R moves the exponent by ~24 eps and the rate by
~1.3e-03 — an order above the 1e-04 this package calls a real difference.

**And it moves nothing at 25 C, exactly.** `arrh_curve` carries `(leaf_temp - 25)` in its
numerator, so it returns its reference value through `exp(0)` whatever R is;
`peak_arrh_curve`'s `arg2` and `arg3` are the same expression at the same temperature, so
their ratio is exactly 1. Predicted before measuring, and confirmed: **1728 of 1728 moved
golden cells are at 40 C, and the 25 C block is byte-identical.**

Blast radius, all at 40 C: median relative 6.49e-05, largest absolute move 2.04e-03. The
two cells above 1e-02 are `profit` at 5-layer 40 C points where it is ~1e-03, i.e. small
denominators. `gradient_golden.tsv` does not move at all, because `grid_drivers` is 25 C.

⚠️ **The R-side suite was blind to this**: all 1107 tests passed unaltered, and only the
golden grid's second temperature caught it — which is what that second temperature is
for.

The reason was **not** that nothing on the R side runs above 25 C; `test-temperature-
response.R` drives at 35 and 40 C. It is that every assertion up there was
**directional** — `expect_gt`, `expect_lt`, "not `all.equal`" — so a change to a
response curve that preserves the ordering satisfied all of them. There was no *pinned*
value off 25 C anywhere in `tests/testthat/`, because `golden_solve()` hard-coded the
temperature.

**Fixed in the same PR.** `golden_solve()` takes `leaf_temp`, and `test-golden.R` pins
two 40 C points (1-layer and 3-layer, both interior) to hex, alongside a check that they
really are a different operating point from their 25 C counterparts. Measured: reverting
`gas_constant` to 8.314 now fails **18** R-side assertions — exactly the two new rows x
nine fields — against **0** before.

⚠️ plant does **not** re-export `gas_constant` (its `leaf_model.h` says so explicitly),
so nothing over there references it by name; its results still move through the leaf.

## Two constants audited and deliberately left alone

The same audit flagged two more, neither of which has a forced answer, so both are
documented at the constant rather than changed:

- **`gravity_head = 9.8e-3`** is 6.78e-04 below `1000 * 9.80665 / 1e6`. Above the
  real-difference threshold, and it moves results through the per-layer head — but the
  right value depends on a water-density convention nobody has stated (nominal 1000,
  4 C 999.97, 25 C 997.05 span 0.3%, four times the discrepancy). Needs a decision.
- **`latent_heat_vap = 2.45e6`** and its own comment disagree: the comment said "fixed at
  25 deg C", where lambda(25 C) is 2.442e6; 2.45e6 is ~21.5 C. Only reaches the
  default-off Penman-Monteith path, and on that path lambda is arguably the wrong *shape*
  rather than the wrong value, since it depends on the leaf temperature being solved for.
  Left to #28; the comment no longer states a temperature the value does not have.

`C_to_K`, `umol_par_per_joule`, `vol_heat_cap_air` and (since #51) `molar_mass_h2o` check
out. `H2O_CO2_stom_diff_ratio` is #50 and is a convention question, not an error.

## One molar mass of water, so the kg <-> mol conversions are reciprocal (#51) — MOVES RESULTS

`kg_to_mol_h2o` was 55.4939 and `kg_per_mol_h2o` was 0.018015 — two constants naming
the same physical quantity in opposite directions, disagreeing by **0.0277%**, used in
opposite halves of the model. The demand side converted transpiration kg -> mol with
the first; the supply side converted uptake mol -> kg with the second. So a water flux
pushed through both did not come back. The header said the discrepancy was deliberate,
"kept at the historical 0.018015 to preserve results".

**Which one was wrong is not a matter of convention, which is what let this be settled
rather than argued.** 55.4939 is 1/0.018020, i.e. it encodes a molar mass of 18.0200
g/mol. The molar mass of water is 18.015 g/mol — from the standard atomic weights,
2(1.008) + 15.999. So `0.018015` is the physical value and the forward constant was the
odd one. There is now one `molar_mass_h2o = 0.018015` and both old names are derived
from it, reciprocal by construction; both names are kept because plant `using`-declares
both.

Note this is the **larger** of the two possible moves — unifying the other way would
have preserved more digits by adopting a molar mass water does not have.

**Blast radius.** 3710 of 5184 golden cells, median relative move 9.4e-05, and **3377
of the 3710 move by no more than 2x the constant's own 2.77e-04** — i.e. the bulk is
that constant propagating. The largest **absolute** move anywhere in the file is
1.6e-03. The 184 cells whose relative move exceeds 1e-03 are near-zero quantities: the
worst, 2.7e-02, is `profit` = 0.0088 at a 40 C five-layer point where benefit nearly
cancels cost. Recorded gradients: 60 of 100 cells, median 3.2e-04, worst 2.9e-03.
⚠️ That worst is the same order as the gradient file's own cross-platform disagreement
(~2.3e-03), so **off macOS/arm64 this change is not cleanly separable from noise in
that file.**

Directionally: the forward constant rose 0.0277%, so conductance per unit transpiration
rose with it, and the leaf buys slightly more carbon for the same water.

⚠️ **plant's `LinkingTo: phylloptim (>= 0.2.0)` floor is now three minor versions
stale**, and this is a results change it should be able to require. `>= 0.4.0`.

# phylloptim 0.3.0

⚠️ **This section was headed `0.2.1` until now, and the renumbering is the point of
#58's first ask rather than tidying.** `0.2.1` was never tagged or released, and it
accumulated `#41` (dark respiration reallocated, which moved results), `#84`, `#86`,
`#89`, `#90` and `#91` on top of the changes below — so a consumer caching computed
results had no signal that any of it had happened. The rule from here: **a PR that
moves results moves the minor version**, in the same PR that regenerates the golden
file.

⚠️ **Downstream, plant still pins `LinkingTo: phylloptim (>= 0.2.0)`, and that floor
is now two minor versions stale.** #99 has merged, so the bump this note asked for has
happened here and the plant-side half has not: a `>= 0.2.0` floor is satisfied by
every build that predates the vulnerability-domain fix, which is exactly the
"consumer cannot tell" problem one repo over. Raising it to `>= 0.3.0` — together with
the `Remotes: traitecoevo/phylloptim@<sha>` it sits beside — is what makes the pin
mean anything.

## `leaf_behaviour_fingerprint()`, so a consumer can tell when the numbers moved (#58)

The version bump above is the discipline; this is the mechanism, because a promise
about future PRs cannot help a cache built against the history it sits on.
`leaf_behaviour_fingerprint()` returns a 12-character digest of the two recorded
baselines together — `tests/cpp/golden/operating_points.tsv` and
`tests/testthat/gradient_golden.tsv` — so a pipeline can depend on one value instead
of reimplementing "has phylloptim changed?".

Those files already *are* this package's definition of the numbers, and they are
regenerated deliberately, so the fingerprint inherits that discipline rather than
needing new enforcement. `test-fingerprint.R` is what holds it: it recomputes the
digest from the files on disk, so a PR that regenerates a golden file and forgets
`Rscript tools/fingerprint.R` fails.

Two properties chosen against the workaround it replaces — the calibration study
hashed `inst/include`, `R/` and `src/`, which #47's rename would have moved without
moving a number. So **file names are not in the digest** (a rename cannot move it)
and **contents are hashed by line with the trailing `\r` stripped** (a CRLF checkout
cannot either). Both are asserted.

⚠️ It reports what the golden grids *reach*. They build a fresh `Leaf` per point, so
a stale-state bug of hazard 8's kind can be fixed or introduced without moving it —
all three of #15's were golden-bit-identical. `?leaf_behaviour_fingerprint` lists
that and the other two limits.

## The vulnerability spline reaches the domain it claims (#92) — MOVES RESULTS

`cumulative_vulnerability_integral` built its knot grid by accumulating
`psi += step`. Rounding accumulates, so after `resolution` additions `psi` landed a
few ULP either side of `psi_max` and the loop yielded **`resolution` or
`resolution - 1` knots depending on the values of `b` and `c`** — an upper domain
bound that was neither reproducible across parameter values nor equal to
`vulnerability_psi_max`, which the comment there says the inverse "needs the same
bound" as. Both splines built from that grid have extrapolation **disabled**, so the
edge decides whether a lookup at the dry end throws.

**The package defaults were on the losing side of it**: 99 knots, ending at 6.8229
against a `psi_max` of 6.8918 — one full step, 1.0% of the intended domain, silently
missing from both the stem and the root curve.

Indexing the grid fixes the count and makes the last knot exactly
`vulnerability_psi_max(b, c)`. **This moves results**: 3496 of 5184 golden cells,
median 2.0e-15, and the largest **absolute** move anywhere in the file is 1.97e-10.
The 71 cells whose *relative* move exceeds 1e-7 are all at `psi_soil = 3` with 5
layers, where the quantities themselves are 1e-11 to 1e-3. Nothing reaches the
~1e-4 band this package calls a real difference. Recorded gradients move further,
as finite differences of the same solve must: 63 of 100 cells, median 1.1e-08,
worst 5.1e-04 — an order inside the 5e-3 that file is already compared with
cross-platform.

## `psi_crit` is checked against the curve that bounds it (#38)

`psi_crit` looks independent of `stem_b`/`stem_c` and is not: the spline stops at
`P99 = stem_b * log(100)^(1/stem_c)`, `psi_crit` never enters that, and every solve
evaluates the curve *at* `psi_crit`. So `psi_crit > P99` used to fail from inside
the interpolator, in a message naming neither trait — which reads as "the solver
went somewhere strange" when the real answer is a trait combination that was never
valid. Found while calibrating against measured P50/P88 curves, where `stem_b` and
`stem_c` are far from the defaults.

The constructor, `set_traits()` and `perturb_stem_b()` now refuse it by name, and
the message quotes the **P95** that would work — because that is what the defaults
encode: `3.898245 * log(1/0.05)^(1/2.680147) = 5.870283 = psi_crit`, to six decimal
places. `?leaf_traits` states the relationship, which was written down nowhere.

⚠️ **The root curve is deliberately not checked.** #38 assumed `root_psi_crit`
carried the same constraint; since #77 bounded the root curves past their last knot
it does not throw, it **clamps** — a different defect, and #85's question.

This rejects only combinations that already failed, so nothing that worked stops
working. It did catch one inconsistent pair inside this repo's own suite
(`leaf_traits(stem_b = 2.0)` at the default `psi_crit`) and one in `?leaf_traits`'s
example, both now fixed.

## Traits can be read back from the object (#95)

The thirteen `set_traits()` traits were write-only from R: `psi_crit`, `stem_b`,
`stem_c`, `beta2`, `root_*` and the rest could be set and not read. Anything that
had to compute a quantity the model defines in terms of one had to carry it —
Sperry's cost normalises by `k_crit = kmax * proportion_of_conductivity(psi_crit)`,
and a downstream replication probe therefore held a
hard-coded `5.870283`.

All thirteen are now bound **read-only**. `set_traits()` is still the only way to
change one, for the reason hazard 10 gives at length: a bare write leaves up to two
vulnerability splines and the solved operating point describing the old value.
`R_d_25` stays settable, as it already was.

## `R_d_25` is a trait, and respiration rises with temperature (#41)

Dark respiration was `rd_to_vcmax_ratio_ * vcmax_(T)`: a fixed fraction of Vcmax,
unreachable from R, and — because it inherited Vcmax's **peaked** Arrhenius — it
**fell** above the thermal optimum, where real dark respiration rises. No value of
the fraction repairs a function of the wrong shape.

Now `R_d_25` is a `leaf_traits()` trait like any other, and the response is
Tjoelker's declining Q10:

```
Q10(T) = rd_q10_intercept_ - rd_q10_slope_ * (T + 25) / 2      3.09, 0.0430
R_d(T) = R_d_25 * Q10(T)^((T - 25) / 10)
```

The Q10 is evaluated at the mean of T and the reference, which is the form the
land-surface literature implements. Set the slope to zero and the intercept is a
constant Q10, for anyone who wants the conventional form; a constant Q10 of 2 was
tried first and measured too aggressive (R_d 2.8× its 25 °C value over 15 K, −73% on
assimilation at 40 °C, no operating point at all by 45 °C).

**There is no fallback and no sentinel.** `R_d_25` defaults to **1.44** — which is
`0.015 * 96`, the value the old derivation gave at the default `vcmax_25`, and the
same double — so a caller who does not touch it gets exactly what they got before.
An unset or negative value fails rather than deriving something.

⚠️ **RESPIRATION NO LONGER FOLLOWS `vcmax_25`.** That is the point of making it a
trait, and it is the one place "same defaults ⇒ same results" has an exception:
a run that *varies* `vcmax_25` used to get proportionally more or less respiration
and now does not. Measured on plant's TF24 SCM, offspring production:

| | phylloptim master | this branch |
|---|---|---|
| default (`vcmax_25` = 96) | 81.857087216483691 | 81.857087216483691 |
| `vcmax_25` = 60 | 0.40692 | 0.34336 |
| `vcmax_25` = 150 | 73.703 | 92.026 |

The default arm is bit-identical to all 17 digits. If a study wants respiration to
track Vcmax, that coupling belongs in the caller's parameterisation — plant's
`TF24_hyperpar` is where every other derived parameter is already computed.

`R_d_25` is also a `pars` entry for `leaf_gradient()` and `leaf_gradient_batch()`:
`dA/dR_d_25 = -0.24874` against a central difference of the solve at −0.24868. ⚠️ Not
≈ −1, because the optimiser moves the operating point in response. `dY/dvcmax_25` is
correspondingly a **partial** at fixed respiration; a fit that moves both adds the
two columns.

**Blast radius**, and it is zero at 25 °C by construction — the trait *is* the value
there:

| T | R_d old → new | A old → new | |
|---|---|---|---|
| 15 °C | 0.669 → 0.646 | 6.1505 → 6.1572 | +0.11% |
| 25 °C | 1.440 → 1.440 | 5.5992 → 5.5992 | **exact** |
| 35 °C | 1.610 → 2.592 | 3.4924 → 3.1550 | −9.66% |
| 40 °C | 1.012 → 3.171 | 1.8265 → 0.8555 | −53.16% |
| 45 °C | 0.507 → 3.618 | 0.6339 → −3.6176 | shut-down |

Two further consequences:

- **A leaf too hot to gain carbon now shuts down instead of throwing.** A higher R_d
  can drive net assimilation negative across the whole `[gamma*, ca]` bracket, and
  then `psi_stem_to_ci` has no supply-equals-demand root. The model already knew that
  case — `ci_at_compensation_point_` — but only on the energy-balance path; it now
  applies on the default path too, and the temperature at which the model stops
  having an operating point is pinned as a measurement.
- **On plant's Penman-Monteith path the change is large**, because leaf temperature
  there reaches **62 °C, up to +22 K above air**: offspring production −47.2%
  (2.2035 → 1.1638). ⚠️ plant's *default* configuration is unaffected — it sets
  `leaf_temp` as a constant 25 °C driver with PM off — so the exposure is entirely in
  configurations that leave 25 °C.

`R_d_25` is the fourteenth entry of the gradient enumeration, i.e. in `set_traits()`'
own argument order, with `leaf_specific_conductance_max` and `resistance` after it.
Every index in `gradient.hpp` now has a name (`par_vcmax_25`, `par_R_d_25`, ...) so
nothing indexes `theta` with a bare integer.

**API removal:** `rd_to_vcmax_ratio_` is gone, from both C++ and the R bindings. It
existed only as the fraction, and there is no fraction now. plant never used it.

## The golden grid gains a temperature axis: 288 points become 576 (#41)

⚠️ **A grid at one leaf temperature is blind to every temperature response in the
model**, and this one was, at 25 °C. The reference values of Vcmax, Jmax and R_d are
*defined* there, so a change to any response curve is inert **by construction** — the
respiration change above moved results at every temperature except 25 °C and this
file did not move a bit.

The grid is now the same 288 operating points at **25 and 40 °C**, with a `leaf_temp`
column. 40 °C is where the response bites hardest, and the extra pinned rows at the
hot end are a feature.

**Regenerated deliberately, and checkable as an addition**: temperature is the
outermost loop, so the 25 °C block is **byte-identical** to the previous file in all
nine output columns. 288 rows added, none moved.

The classification is counted per temperature now, because it moves:

| T | interior | pinned wet | pinned dry | shutdown |
|---|---|---|---|---|
| 25 °C | 198 | 24 | 18 | 48 |
| 40 °C | 160 | 80 | 0 | 48 |

A single total would let points move between branches and still add up. The direction
is physical: a hot leaf assimilates less and respires more, so it has less to gain
from water and its optimum presses against the **wet** bound instead of the dry one.
The 48 shut-down rows do not move, because there it is hydraulics rather than heat
that forbids transpiration.

`tests/cpp/bench_solve.cpp` deliberately stays at 288 points at 25 °C: a timing
baseline is only useful against its own history (`tools/cost-baseline.tsv`,
`tools/bench_history.sh`).

## The recorded gradients have their own tolerance

`gradient_golden.tsv` used the solved-output baseline's tolerance, and that was
wrong. A value there is a **finite difference**, so cross-platform it carries the
solve's ~1e-9 floor **divided by the step** — one amplification more than the outputs
it is built from. The step is relative, so the smallest-magnitude parameter sets the
tolerance for the whole file: `R_d_25` is 1.44 where `vcmax_25` is 96, a 67× smaller
absolute step, and it disagrees 1.3e-03 on Linux against 1.3e-04 for the other
columns. `gradient_golden_tolerance()` is 5e-03, with the arithmetic.

## The temperature cache now keys on everything it reads, so R_d is genuinely settable (#41)

Binding the temperature-response parameters to R made them *writable* but not
*effective*. `set_physiology()`'s temperature block was cached on
`(leaf_temp_, atm_o2_kpa_)` alone, justified by "same inputs → bit-identical
outputs, so reusing is exact" — a statement that was true while those parameters
were unreachable C++ members and **became false the moment they were bound**. The
block's outputs depend on fourteen further inputs the key never mentioned.

The failure was silent. Setting `rd_to_vcmax_ratio_ <- 0.03` on a solved leaf and
re-supplying the same drivers left `R_d` at **1.44** where a freshly built leaf gave
**2.88**, and assimilation unchanged to every digit. Which is the outcome #41 cared
about: a calibration that cannot move respiration absorbs the mismatch into whatever
it *can* move.

The cache key now covers every scalar `update_temperature_dependent_params()` reads.
Two consequences beyond the reported bug:

- **It also covers `vcmax_25` and `jmax_25`**, which closes the third and least
  visible part of hazard 10 — a bare `l$vcmax_25 <- x` no longer leaves
  `vcmax_`/`jmax_`/`R_d_` describing the old value. ⚠️ `set_traits()` is still the
  correct way to change a trait; the vulnerability splines and the solved operating
  point need clearing too, and a cache key cannot do that.
- ⚠️ **The existing test passed while the feature was broken**, because it pushed
  each change through by calling `update_temperature_dependent_params()` directly —
  a route no caller has. `test_temperature_params_invalidate_cache` uses the real
  one, and fails four ways without this fix.

**No behaviour change at the defaults**: the golden file is bit-identical, and the
solve is 2.99 against 3.00 µs/solve interleaved ×3, inside the ±0.01 within-process
noise.

⚠️ **What this does NOT fix**, recorded because the numbers invite it: `R_d` still
inherits Vcmax's *peaked* Arrhenius and therefore **falls** above the thermal
optimum, where real dark respiration rises. No value of the ratio repairs a function
of the wrong shape — see the note in `update_temperature_dependent_params()`. And
`rd_to_vcmax_ratio` is still not a `leaf_traits()` member, so it cannot yet be
fitted or differentiated; it is set as a field.

## Fick's law uses the LEAF-to-air vapour deficit, not the air's (#7)

Transpiration and stomatal conductance divided by `atm_vpd_` however hot the leaf
got. Diffusion out of a stoma is driven by the deficit *at the leaf*, and on the
energy-balance path the leaf runs above air temperature, so that deficit is the
larger one. Measured at the package defaults it is **3.4× / 3.2× / 4.0× the air's
at Tair 25 / 35 / 45 °C** — the factor by which `gs` and `A` were overstated
([#7](https://github.com/traitecoevo/phylloptim/issues/7), PLAN 13.1).

⚠️ **Off the energy balance nothing moves, and that is a property of the grid
rather than evidence the change is inert.** At air temperature
`vpd_leaf_ == atm_vpd_` exactly, so the prescribed-temperature path is untouched
and `tests/cpp/golden/` is **bit-identical across all five commits** — because
the grid runs with `use_energy_balance_` off. A bit-identical golden run says
nothing here.

⚠️ **The decoupling signature moves from conductance to transpiration.** Over the
acceptance window `E` rises ×1.034 against `D_leaf` ×1.089, so `gs` falls ×0.950.
A test written against the old arithmetic would read the same physics off the
wrong variable.

`constants::vpd_leaf_min = 0.01` kPa floors it. A leaf cooled below the dew point
has a *negative* deficit, and dividing a positive transpiration by it reports a
negative conductance; the floor keeps such points finite and unattractive instead
of sign-inverted. The driver defaults are 1.5–2 kPa, so it binds only where a
one-way diffusion equation has already stopped describing the leaf.

## Sperry's ProfitMax, on the same footing as TF24

`$optimise_psi_stem_ProfitMax()` maximises Sperry (2017)'s profit with **both**
terms normalised, and reports the λ that makes the existing entry point agree.
`$profitmax_curve()` returns the whole cost, gain and profit curve in one call,
and `$optimise_psi_stem_TF()` runs the TF24 cost through the same solver — so the
two formulations can be compared at identical drivers. An optional instantaneous
thermal cost is included, **default off**.

⚠️ **These write `$carbon_gain_`, `$hydraulic_cost_norm_` and `$thermal_cost_`,
which are unitless and are NOT `$hydraulic_cost_`.** Reading the normalised cost
as the TF24 one is a units error the names are chosen to prevent.

The normalised hydraulic cost is invariant to `kmax` by construction, and it is
asserted rather than assumed: worst difference **1.1e-16** across a 3× change.

Motivated by a downstream replication of Sicangco et
al. (2026), which needed ProfitMax measured against TF24 rather than described.

## ⚠️ Four fixes on the single-layer optimisers, and the fourth is a solver bug

These paths are reachable only from `optimise_psi_stem_*`, not from
`find_root_collar_psi()` — so **`plant` is unaffected, provably rather than
probably**: `prepare_collar_solve` has its own `assim_max_ < 0` exit and returns
false before any candidate potential is evaluated.

1. **`optimise_psi_stem_Sperry` searched a NaN objective silently.** `lambda_` is
   an input with no default, cleared by `setup_clean_leaf` and never set by
   `set_physiology`, so a caller who drove the leaf and called this got a
   plausible potential (2.551266 MPa at the defaults) beside `profit_ = NaN`. It
   now refuses.
2. **Hazard 8, live on this path.** `opt_root_psi_`, `E_up_` and
   `soil_consumption_` survived from an earlier collar solve, so the object
   reported `E = 9.216e-5` beside `E_up = 2.626e-5`.
3. **`set_leaf_states_rates_from_psi_stem` zeroed transpiration wherever
   `assim_max_ < 0`.** Transpiration there is the hydraulic supply and does not
   depend on photosynthesis, so the branch made two leaves at the same operating
   point disagree about whether water was moving — purely because one had
   respiration switched on. The shut-down *state* it reached for is unchanged:
   the `ci` root-find has no root in `[gamma*, ca]` and takes its
   compensation-point fallback, which is what the branch set by hand. That
   fallback did not exist when the branch was written.
4. **`optimise_psi_stem_ProfitMax` now scans a grid before refining.**
   `brent_fmin` steps in from the bounds, so it returns neither an endpoint nor
   the global maximum of a multi-modal objective — and at Tair 50 °C with the
   thermal cost on it returned **1.643 MPa where the objective is maximised at
   full closure** (−1.5314 at `psi_soil` against −1.5459 at the interior local
   max), reporting an open stoma where the model says shut. The scan reuses
   `prepare_profitmax`'s own grid, so this costs no extra model evaluations.
   Found by comparing against Sicangco et al.'s Figure 4: their model closed at
   48–56 °C and ours did not, and the whole difference was the search.

⚠️ **`optimise_psi_stem_TF` and `optimise_psi_stem_Sperry` still have fault 4**,
and are documented rather than fixed because neither has a scan to reuse. The
general form is now hazard 11 in the developer guide: *a bracketing optimiser
answers "where is the interior maximum", which is not the same question as "where
is the maximum".*

## The `dgc_dT` term in the energy-balance derivative was re-derived, not sign-flipped

`dprofit_energy_balance_term` carried a named `const double dgc_dT = 0.0` and a
note promising PLAN 13.1 would make it "a one-line change instead of a
re-derivation". It would not have been: the damping factor it multiplied puts the
new term under `A_T` where the derivation puts it under `A_prime`. The two agree
**only at `dgc_dT = 0`**, which is why nothing caught it. The derivation is
written out at the function, and at `dgc_dT = 0` it reduces to exactly the old
expression — so the prescribed-VPD behaviour is unchanged, which is the other
half of why the golden file does not move.

## A comment block may open only one `\verbatim` run, and only CI can see it

Doxygen 1.9.8 — what the runner installs — handles the first `\verbatim` in a
`///` block and drops the second one's OPEN, then reports `unexpected command
endverbatim` at a line in the *filtered* stream that lands in unrelated code: 120
lines past the responsible comment, in the case that found it. Doxygen 1.17
renders the same input in silence, so a local run is not an oracle. `docs.yml`
now counts opens per block and names the offending block.

Established with a probe header of eight isolated constructs in one render: one
run, banner rules, `|` in prose, `|` inside a run, and emoji are all clean; only
two-runs-in-one-block errors. ⚠️ **Two earlier explanations were asserted before
being measured and are both wrong** — a setext heading swallowing the block, and
a bare `\|` opening a Markdown table. A `\|`-escaping change to
`tools/doxygen_filter.awk` had been committed on the second and written into the
filter's header as fact; **that commit is reverted**, and the filter's header now
records the measurement instead.

⚠️ The `/*! \file */` block is exempt from the count, measured rather than
assumed: `closed_form.hpp`'s file block has carried two runs across a long green
master, so counting it would fail the build on code Doxygen renders happily.

C++ suite **486 → 534 checks**, 0 failures.

## A trait gradient at a collar potential the caller supplies

`leaf_gradient()` and `leaf_gradient_batch()` take `psi`, and evaluate there
instead of solving for the profit-maximising collar
([#88](https://github.com/traitecoevo/phylloptim/issues/88)). Until now both
called `find_root_collar_psi()` unconditionally and read `opt_root_psi_` back,
so a model that **tracks** the optimum rather than finding it — plant's TF24f
carries the collar as an ODE state, `dpsi/dt = k * dprofit/dpsi` — could not ask
this package for a trait gradient at the point it was actually operating at. It
got a confident answer about the re-solved optimum instead, with nothing saying
so.

The maths simplifies rather than breaks: `psi` is exogenous, so the indirect term
is whatever the caller says it is, via `dpsi_dtheta` (default zero, the partial
at fixed collar). Nothing is derived from `-M/H`, so nothing needs stationarity,
and `method` is refused — the two routes it chooses between are both about a
solved optimum.

`M`, `H`, `dY_dpsi` and `psi` now come back in the result on both paths. A caller
whose `psi` is *dynamic* cannot supply `dpsi_dtheta` as a constant: for the
gradient-ascent law above it obeys `ds/dt = k(M + H s)`, and those are its
coefficients (traitecoevo/plant#614).

⚠️ **`stationarity` is still computed on the prescribed path, and now means
something better.** It no longer routes anything — it measures how far the collar
you supplied sits from the optimum. It makes exactly one decision, `profit`'s:
at a stationary point the envelope theorem applies and the analytic zero is used;
away from one the *exact* `dprofit/dpsi` is used rather than a difference of it.
One rule, both paths — which is why `psi = <the solved psi*>` with
`dpsi_dtheta = -M/H` reproduces the solving path **bit-for-bit**, asserted with
`identical()` rather than a tolerance.

⚠️ **A clamped `psi` returns no gradient, and this is the case to understand.**
The collar actually used is `psi` clamped into the feasible interval, so it moves
with the *bound* rather than with `dpsi_dtheta` — the active-set problem arriving
through the clamp instead of through the optimiser, where the direct term alone
is plausible and wrong. `status` reports `"clamped"`, the gradient is `NA`, and
`psi` in the result is the collar that was used. Reported rather than thrown
because a tracking model reaches these points routinely: the clamp is how TF24f
pulls an out-of-range state back inside. It also fires for a `psi` within one
step of an end, where `dY/dpsi` cannot be centred.

The solving path is unchanged and bit-identical, including `gradient_golden.tsv`.

⚠️ **An INFEASIBLE prescribed `psi` is `"no-gradient"`, not a sentinel zero.**
`dprofit_droot_collar_psi` returns a hard `0.0` on its shut-down and
reversed-gradient exits, and a bare zero is indistinguishable from a stationary
point. The solving path got away with reading the value alone because `H`
collapses to zero with it and `usable` catches the pair; the prescribed path
never divides by `H`, so it would have adopted the sentinel *as* `dprofit/dpsi`
— silently losing profit's indirect term at exactly the dry points a tracking
model lives in. Most such points are caught as `"clamped"` first, but not the one
where the caller hands back the collar the shut-down state itself seated.

**`Leaf$dprofit_droot_collar_psi_checked()` is new and is what makes that
possible.** The `bool* feasible` out-parameter has been there since #79 and the
C++ vignette has always said a composite ignoring it inherits the bug — but
RcppR6 has no form for a `bool*`, so the generated binding dropped it and every
R-side composite *was* that composite. It returns `{dprofit, feasible}`.

C++ consumers get the same through `gradient::Prescribed` and the new
`psi`/`dpsi_dtheta` arguments to `gradient::batch`. `Status` gains `Prescribed`
and `Clamped` — ⚠️ **appended after `Error`, so no existing integer value
moves**, and `status_name`'s switch is exhaustive with no `default:` so the next
member added is a compiler diagnostic rather than a silent `"error"` label.

## The gradient differentiates `profit`, which is what a demographic caller bills

`leaf_gradient()` and `leaf_gradient_batch()` return a fifth column. The four
that were there — `A`, `gc`, `psi_stem`, `collar` — are what a gas-exchange
calibration observes, and they were chosen for the customer this feature was
built for. They are **disjoint** from what `plant` reads off a solved leaf: its
carbon budget is `leaf.profit_` (not `assim_colimited_`) and its water budget is
`leaf.soil_consumption_`. So no trait gradient this package produced reached a
demographic model at all, at the optimum or anywhere else
([#87](https://github.com/traitecoevo/phylloptim/issues/87)).

⚠️ **`profit` is the one output the envelope theorem reaches, and the only place
this package uses it.** At an interior optimum `dprofit/dpsi = 0`, so the
indirect term `(dprofit/dpsi)(dpsi*/dtheta)` vanishes identically and
`dprofit/dtheta` is the direct partial at fixed ψ — no `dY/dpsi`, no `−M/H`. It
is set from that term rather than computed through the composite, exactly as
`collar` is set from `dpsi*/dtheta` rather than differenced.

**That is a numerical decision, not a tidiness one, and the measurement is the
reason.** The dropped term is *noise*, not an `h²` truncation: `profit` is the
maximum, so it is flat, and a central difference of it divides the solve's ~1e-09
floor by a ~1e-06 step. Over the golden grid's 136 interior rows —

| | median | max |
|---|---|---|
| `\|dprofit/dpsi\|`, exact (forward AD) | 4.8e-15 | 5.4e-10 |
| `\|dprofit/dpsi\|`, central difference | 7.8e-10 | 2.1e-04 |
| relative move in `dprofit/dtheta` if kept | 2.7e-10 | 8.0e-05 |

— eleven orders between the two instruments at the median, and the worst row sits
in the band this repo calls a real difference rather than rounding. The identity
is applied only where `status == "interior"`; at a pinned optimum `psi*` is a
trait-dependent bound, `dprofit/dpsi` is not zero, and `profit` takes the same
finite-difference fallback as the other four.

The four existing columns are **bit-identical** — this is additive, like
appending to `gradient_par_names()`. `tests/testthat/gradient_golden.tsv` gains a
column and no existing cell moved, checked against master rather than against the
branch point.

⚠️ **The shut-down row's profit column is asserted against a closed form, not
only recorded.** It is the one regime where `profit_` is written by a branch that
leaves the other outputs alone (hazard 8), so a hex with nothing saying what it
ought to be would pin a number rather than a fact. There `E = 0`, so `A = -R_d`
exactly and the hydraulic cost does not depend on `R_d_25`: `dprofit/dR_d_25` is
**−1**. And the shut-down collar is pinned at `psi_crit`, so
`dcollar/dpsi_crit` is **1** — which is why `psi_crit` alone carries a non-zero
profit gradient there.

**`gradient_output_names()` is exported**, and R now *reads* the list rather than
keeping a second copy. `gradient_par_names()` has to be duplicated-and-compared
because R builds `theta` before any C++ call; the outputs have no such
constraint, so adding one is a single edit.

⚠️ **`uptake` was considered and is not here.** Every output must be a field R
*copies* out of `operating_point_values()`; `uptake` is one R *computes*, by
summing over the finite soil layers, so adding it means reproducing that
summation and its order on the C++ side too. That is a separate decision from
this one.

## An out-of-domain transport lookup says which spline, and which caller

The stem curve is the only interpolator here built with extrapolation disabled, so
it is the only one a lookup can throw on — and there are **two** of them,
`transpiration_from_psi` and its inverse `psi_from_transpiration`. They are
inverses, they carry different units, and they are read from four places. odelia's
message names the point and the domain but cannot name the spline or the caller,
because it does not know either. Now:

```
Leaf hydraulics: psi_from_transpiration (the INVERSE cumulative xylem
conductivity integral G^-1, argument in E/K_max) evaluated outside its domain:
E/K_max = -3.18471e+07 lies 3.18471e+07 beyond the lower end of [0, 3.46004];
asked by Leaf::transpiration_to_psi_stem, inverting for psi_stem.
```

That lower-end failure is the
[traitecoevo/plant#576](https://github.com/traitecoevo/plant/issues/576)
signature, and reading it is the whole point: a *negative* `E/K_max` says the
collar cannot supply the demanded flux, so the stem potential that would carry it
is wetter than saturation and no widening of the domain can help — the caller
should not have asked. Localising #576 without this meant instrumenting four call
sites by hand to find out which spline was being read and at what value.

Under a `stem_b` rescale the domain is reported in the **caller's** units, not the
spline's, and the rescale is named. The value handed to the spline is `psi / s`, so
quoting the spline's own endpoints would send the reader after a discrepancy that
is not there.

Behaviour is unchanged: 288 golden operating points are bit-identical, and
`find_root_collar_psi()` is unchanged at 3.16 µs. A non-finite argument still falls
through to the spline and returns non-finite rather than throwing — the guard is
`v < min || v > max`, not the negation of an in-range test, and plant documents a
`profit_psi_stem_TF(NA, .) -> NA` contract built on exactly that.

# phylloptim 0.2.0

## A near-embolised root is no longer an ever-stronger pump

The two root vulnerability curves stop at the 1%-conductivity point (6.82 MPa at the
root defaults), and both were left extrapolating. Past the last knot the conductivity
curve crossed zero at 7.3742 MPa and reached -20.35 at 1000 MPa, while the cumulative
integral kept accumulating past a limit it had already reached to 99.83%. The integral
is the denominator of the per-layer mean resistance `r_R_H = r_R_H_min * span /
integral`, so resistance *fell* as a layer dried: the drier a near-embolised layer, the
harder the plant pumped water into it, and whole-plant shutdown keys off the wettest
layer so nothing stopped it.

Both curves are now bounded, and read only through accessors so a bare `.eval` cannot
reintroduce the old behaviour. The integral is capped at its closed form
`G(inf) = (b/c)*Gamma(1/c)`; the conductivity lookup clamps its argument to the last
knot and its spline refuses extrapolation, matching the stem pair. Measured on one
rooted layer with unit horizontal resistance and the collar at 1 MPa, the reverse flux
into a 1000 MPa layer goes from -14.08 to -2.47 mol H2O m^-2 s^-1 — 5.70x too high
before, and now the whole area under the conductivity curve over that interval, which is
the right limit. A tenfold drier layer no longer pumps 4x harder.

**Golden-identical.** Instrumented over the grid, the driest argument the integral ever
sees is 7.0 MPa: past the last knot, short of the 7.3132 MPa where the cap binds. All 288
rows x 9 fields are bit-identical to a baseline generated on the same machine.
`test_root_vulnerability_is_bounded_past_its_grid` is what stands behind the fix
instead. Costs +2.0% on the collar solve (7.33 -> 7.48 us, interleaved x11); a no-op
wrapper around the same `.eval` measures 7.43, so two thirds of that is the extra call
rather than the bound.
## Trait gradients over a batch of observations, composed in C++ -- 22x

`leaf_batch()` and `leaf_gradient_batch()`. The same gradient `leaf_gradient()`
computes, by the same two routes with the same active-set test, but composed in
`inst/include/phylloptim/gradient.hpp` and looped over observations there -- so a
calibration crosses the R boundary **once per likelihood evaluation** instead of 112
times per observation.

```r
b <- leaf_batch(psi_soil = obs$psi_soil, PPFD = obs$PPFD)   # once per fit
g <- leaf_gradient_batch(b, traits, pars = FIT)             # once per draw
```

Measured per observation over 24 gradients at four DIFFERENTIATED parameters --
`length(pars)`, not the number an optimiser moves -- both arms in one
process:

| | us/observation | x a trivial `.Call` |
|---|---:|---:|
| `leaf_gradient()` in a loop | 363.3 | 340 |
| `leaf_gradient(x = )`, reusing one leaf | 235.2 | 220 |
| **`leaf_gradient_batch()`** | **10.59** | **9.9** |

`leaf_batch()` costs 335 us once, for all 24. For a 1,327-observation fit at 30,000
draws that is minutes rather than hours.

⚠️ **Nothing about the gradient got faster.** The C++ model was already 1.5% of a
gradient's cost; this removes the boundary from under it. So the figure needs a batch
to be realised -- calling it with one observation pays the whole per-call overhead for
one row's work.

**Not a likelihood, and not the parameterisation Jacobian.** Both stay in R: the
likelihood is your model, and the chain rule belongs where the win is -- C++ returns
`dY/dtheta` for the four parameters a leaf has, and R applies your `P_fit x P_model`
Jacobian vectorised over observations.

**Per-row status, not an error.** A proposal will reach operating points the solve
cannot handle, and that costs those rows rather than the dataset: every row reports
`status`, a failure is `"error"` with a `message`, and a failed row's gradient is all
`NA` rather than partially filled.

⚠️ **`leaf_batch()` holds C++ pointers, so it does not survive `saveRDS()` or a new
session.** Rebuild it; it says so rather than crashing.

⚠️ **`psi_soil` recycles the way `leaf_solve()`'s does.** A plain numeric vector is N
single-layer observations; a list of numeric vectors is one multi-layer observation
each. This is the easiest mistake to make with `leaf_batch()`, and it fails loudly
rather than silently.

The C++ composite reproduces `leaf_gradient()` **bit-for-bit** -- 0 mismatches over 25
operating points x three methods, including every pinned and shut-down row -- and the
values are also pinned to `tests/testthat/gradient_golden.tsv`, because an equality
test between two implementations cannot see a change applied to both. PLAN item 11g
has the three things that made bit-for-bit possible, including that
`a + b * c` written as one expression compiles to a fused multiply-add and disagrees
with R's two roundings on 28% of random triples.

## `pars` order no longer changes the `stem_b` gradient (#72)

`leaf_gradient()`'s `stem_b` gradient was contaminated by whichever parameter preceded
it in `pars`, by up to **3.4e-5 relative** -- four orders above the ~1e-9 the function
documents as achievable, and on both routes. `perturb_stem_b()` rescales the
vulnerability spline and touches nothing else, while the loops only restored the base
point *after* the whole parameter loop, so a `stem_b` that was not first was
differentiated one step away from base in the preceding parameter.

Fixed by restoring the invariant rather than special-casing the name: **every
parameter's gradient is taken from the base point**. Costs +0.15 us per observation
(0.7%), because `set_traits()` decides its spline rebuilds by comparing the trait pairs
it is given -- so after a parameter that owns no vulnerability curve nothing is
rebuilt.

⚠️ **Gradients move.** 9 of 60 recorded cells in `gradient_golden.tsv`, worst 3.8e-7 --
much smaller than the defect was worth, because every recorded case put `vcmax_25`
immediately before `stem_b` and that is a benign predecessor. `test-gradient.R`'s
arbitrated references are unchanged, having been established with `stem_b` first. If
you have fitted results that differentiated `stem_b` alongside `a`,
`curv_fact_colim`, `root_b` or `stem_c`, they carried the larger error.

The regression guard is a test comparing two `pars` orderings **in the same process**,
not the golden file: the golden file is compared with a 1e-3 tolerance off macOS/arm64,
so it could not have caught a 3.4e-5 recurrence on Linux.

**Also found, and not fixed: [#74](https://github.com/traitecoevo/phylloptim/issues/74).**
The `stem_b` shortcut is undone by a spline rebuild once per observation, so PLAN 11f's
24.5x is **2.4x** through `leaf_gradient_batch()`. Say which figure you mean.

## The gradient's R glue is a third cheaper

Profiling `leaf_gradient()` found that the **C++ model is 1.5% of its cost**: 6 us of
solving against 119 us of `.Call` dispatch over 112 boundary crossings and ~285 us of
R interpreter. `.gradient_setter()` alone was 60% of a gradient. Three changes, no C++:

* **The drivers are resolved once per gradient, not once per perturbation.**
  `set_drivers()` is split into `.resolve_drivers()` (validate and default) and the
  application, and the setter calls the former once and applies the result
  positionally. There is still ONE definition of the defaulting rules -- a second copy
  in the gradient code would have been free to drift from `set_drivers()`.
* **Traits go on positionally**, rather than rebuilding a `leaf_traits` object with
  `structure(as.list(...))` and re-extracting thirteen fields by name per
  perturbation.
* **`.gradient_outputs()` reads all four outputs in one call** instead of four R6
  active bindings: 4.65 -> 0.93 us, eleven times per four-parameter gradient. `$` was
  12.7% of a gradient's self time.

Measured per observation over 24 gradients, interleaved three times against the
commit this lands on:

| `length(pars)` = 4 | before | after |
|---|---|---|
| fresh leaf | 504 us | 366 us (-27%) |
| reused leaf | 400 us | **231 us (-42%)** |

Together with the entry below, a four-parameter calibration goes from 504 to 231 us per
observation, **-54%**.

**Gradients are bit-identical** across both supply paths, both methods and eight
parameters -- asserted rather than assumed, since a change applied to both routes
would otherwise pass every existing test.

⚠️ Two things to know before copying the pattern. The positional trait call **bypasses
`set_traits()`'s invariant checks**, which is sound only because the values came from a
`leaf_traits()` the caller already supplied and the C++ setter asserts #25 itself. And
it is a hard-coded thirteen arguments that would not fail to compile if a trait were
added, so `test-gradient.R` asserts the arity.

## `leaf_gradient()` can reuse a leaf, which is 38% of a call

Closes #52. `leaf_gradient()` built its own `Leaf` every call and offered no way to
pass one in, and construction is **~150 us, half of a one-parameter gradient** -- the
largest single term on the R surface, paid once per observation by a fit that
differentiates per observation. `leaf_gradient(x = l, traits = tr, ...)` reuses one. Measured over 24 gradients,
per observation: **511 -> 409 us at four differentiated parameters (-20%)** and
316 -> 200 us at one (-37%).

⚠️ **The four-parameter figure is the one a calibration gets, and it is the smaller
one.** Construction is a fixed cost per call while the differentiated work is linear
in `pars`, so the saving falls as `pars` grows. Quoting the one-parameter number
overstates what a fit sees.

⚠️ **This does not make the exact gradient the faster arm**, and the issue was right
to say so before the work started. `vignette("fitting")` measures it at ~4.8x the wall
clock of differencing the objective while using ~5x fewer solves; recovering a third
leaves it ~3x. The rest is the R call boundary, which only composition in C++ reaches
(PLAN item 11 stage 2).

Three things about the argument:

* **`x` is a vessel, not the point.** `traits` still says where the gradient is taken
  and is applied to `x` rather than read from it -- a `Leaf` does not expose its
  traits. So `traits` is **required** with `x`: defaulting it would differentiate at
  the package defaults on somebody else's leaf and return a plausible answer for a
  point nobody asked about.
* **`control` and `supply` are refused alongside `x`**, since both were fixed when it
  was built and accepting them invites a silent disagreement. The supply path is read
  off the object, so `resistance` is differentiable exactly when `x` is on the
  single-potential path.
* **`x` comes back solved at the base point**, not at the last perturbation, restored
  on the way out including on an error -- hazard 8. That costs one re-trait and one
  solve, ~18 us against the ~150 saved.

`tests/testthat/test-cost.R` asserts the reuse path constructs **zero** leaves and the
default path exactly one per call, so the argument cannot quietly become decorative.

## `series_resistance()` no longer pays the 58 us R-boundary cost

Bug fix on the entry above. `series_resistance()` built its `RootNetwork` through the
generated constructor, which reaches C++ for a default-constructed struct and pays
`Rcpp::wrap` assembling the five-element named list — **58 us**, against 1.1 us for a
trivial `.Call`. It now overwrites one field of a session-cached prototype instead:
**58 -> 3.1 us**.

This matters because, unlike a fixed root network, a fitted series resistance cannot
be hoisted out of a caller's loop — it changes every proposal. Measured in the
companion calibration study, whose likelihood is 24 rows: the migration to the new
API had made `leaf_predict()` **17% slower** (0.517 -> 0.606 ms) where it was
predicted to get slightly *faster*, because removing one object-resetting
`$set_supply_single()` call saved less than the per-evaluation network build cost.

The prototype is shared, and safe only because R copies on write; `test-surface.R`
pins that a returned network cannot corrupt it and that the field list still comes
from the real constructor rather than a hand-written `structure()` that could drift
from the C++ struct.

## The two supply paths are configured and driven the same way

Follows the entry below, and finishes what it started. `set_physiology()` now takes
the soil-to-collar resistances on **both** supply paths, out of the same
`root_network` argument. Before, the multi-layer path took its resistances per call
and the single-potential path took its resistance at construction — so the same
quantity arrived at a different *time* depending on which path was in force.

* **`leaf_supply_singlelayer()` no longer takes `resistance`.** Migration:
  `leaf_supply_singlelayer(resistance = r)` ->
  `leaf_supply_singlelayer()` plus `root_network = series_resistance(r)` on
  `set_drivers()` / `leaf_solve()` / `leaf_gradient()`. The C++
  `$set_supply_single(resistance, gravity_head)` becomes
  `$set_supply_single(gravity_head)`.
* **`series_resistance(r)` is new**: the single-potential path's counterpart to
  `root_network_from_carbon()`, packaging one series resistance into the
  `RootNetwork` that `set_drivers()` takes. It goes in `r_R_V_sum`, which is exact
  rather than a convention — that field already means "series resistance to the
  collar, with no vulnerability weighting". A network carrying a non-zero
  `r_R_H_min` is **refused** on this path rather than silently reinterpreted.
* **`set_drivers()` accepts `root_network` on the single-potential path** where it
  used to reject it, and defaults it to a nominal 1e3 written out in the body — the
  same treatment the multi-layer default gets. `soil_depth` is still refused there,
  because nothing on that path reads a depth profile.
* **`single_resistance_` is unset until `set_physiology()` runs**, like `psi_soil_`.
  It used to be seated by the constructor.

**What this buys, beyond symmetry.** `resistance` was the only differentiable
parameter whose setter reset the whole object: perturbing it meant calling
`$set_supply_single()`, which calls `setup_clean_leaf()`, wiping the temperature
cache and the solved state. It is a plain driver now, handled exactly like
`leaf_specific_conductance_max`, and `.gradient_setter()` makes one
object-resetting call instead of two.

⚠️ **One asymmetry is left, deliberately: `gravity_head`.** The multi-layer path
derives a per-layer head from the depth profile it is handed
(`gravity_head * z_soil_mid`); the single-potential path has no depth profile to
derive one from, and a bare leaf wants **zero** rather than a geometric default —
which is the caller that path exists for. Making it a driver would mean inventing a
depth for a leaf that has none, or adding a second supply-shaped argument only one
path reads. For the multi-layer rule at one layer of thickness `d`, pass
`gravity_head = 0.00981 * d / 2`.

**Results are unchanged on both paths**, verified bit-identical over 18 operating
points: 3 soil profiles x 2 VPDs on the multi-layer path, and 3 resistances x 2
gravity heads x 2 VPDs on the single-potential one. The 288-point golden file is
bit-identical too.

⚠️ **A caution about how that was verified, because it cost time.** The first
comparison said one row differed by 1.22e-10 and that the *same build* was
nondeterministic run to run. Both were artefacts of a stale `src/*.o`: the two arms
had been installed either side of a `git stash` round-trip, and `R CMD INSTALL`
silently relinks rather than recompiling. Clean-building both arms gave
bit-identical, reproducible results. **Check that an install actually compiled
(`grep -c '^clang++.*-c ' <log>`) before believing a sub-solver-floor difference.**

## `set_physiology()` takes root resistances, not root carbon

**Breaking, in both the C++ and R interfaces**, and coupled with a plant PR.

`Leaf::set_physiology()`'s first argument is now a `RootNetwork` — the per-layer
root hydraulic resistances, per unit leaf area — where it used to be a root carbon
profile. The solve never touched root carbon: `uptake()` and `duptake_dpsi()` read
exactly two vectors, `r_R_H_min` and `r_R_V_sum`, and taking carbon made the leaf
own four things that are not gas exchange — `beta_R_H`, `beta_R_V`, the layer
thickness `dz`, and the 1/3 : 2/3 vertical/horizontal split. It is the move
`leaf_specific_conductance_max` already made: plant computes `kmax` from height and
hands over a scalar, because which conductance-versus-height model is in force is
not this package's business. Which root-architecture model is in force is not
either.

What follows:

* `beta_R_H` and `beta_R_V` are no longer `Leaf` constructor arguments, no longer
  `leaf_traits()` entries, and no longer arguments to `set_traits()`. The C++
  constructor takes 17 arguments rather than 19; `set_traits()` takes 13 rather
  than 15.
* `root_network_from_carbon()` — the architecture model itself — **stays**, is
  public, is tested, and is now exposed to R. It takes the soil-depth profile
  rather than `dz`, so a caller never has to reproduce the layer-thickness rule;
  `phylloptim::layer_thickness()` is the one definition both sides use.
* `set_drivers()` and `leaf_solve()` take `root_network` in place of
  `root_carbon_per_leaf_area`. The default is unchanged numerically — a nominal
  20 kg C m^-2 leaf split evenly — but it is now written out as a call to
  `root_network_from_carbon()` instead of being a number in a signature.
* `RootNetwork()` builds one from R. Only `r_R_H_min` and `r_R_V_sum` matter; the
  other three fields are diagnostics and may be left empty. So a caller with a
  measured or fitted series resistance and no carbon profile can now state what
  they have, which was the point.

**Results are unchanged.** The 288-point golden file is bit-identical, and the R
tie-back to it still passes, because the golden grid now calls
`root_network_from_carbon()` and then `set_physiology()` — the same two steps, with
the boundary moved between them.

**Cost: +0.7% per solve** (2.74 -> 2.76 µs/solve, min-of-2000 interleaved over six
rounds, non-overlapping, identical checksums). `set_physiology()` now copy-assigns
five vectors instead of filling them, where it used to run the architecture model
itself; against that, `set_physiology()` alone got cheaper (0.106 -> 0.083 µs).
`set_root_network()` takes a `const&` and copy-assigns rather than taking by value
and moving: a caller that holds a `RootNetwork` as a member and refills it — which
is what plant does — would have its buffers stolen by a move and reallocate all
five vectors on the next call.

**⚠️ A performance trap the R layer had to be written around.** Materialising a
`RootNetwork` as an R object costs **~58 µs** — the bare `RootNetwork__ctor()`
`.Call`, before any R wrapper, against 1.1 µs for a trivial `.Call` and 4.3 µs for
`$set_physiology()` handed a ready-made network. It is `Rcpp::wrap` building the
five-element named list. Constructing the default one per call made `set_drivers()`
9.4 → 73 µs and `leaf_solve()` 23.7 → 108 µs/row, a 4.5× regression that every test
passed. `set_drivers()` now memoises the default network and the single-potential
path's placeholder, both in size-one memos keyed by `identical()` on `soil_depth`,
which restores 9.6 / 5.8 µs and 24.1 µs/row — master's figures within noise. A
`leaf_solve()` sweep holds `soil_depth` fixed across rows, so it pays the
construction once. **If you build networks in a loop, build them outside it.**

**A new check, because the length agreement is no longer free.** `max_soil_layer`
indexes `psi_soil` and the per-layer gravity head directly, so a network with more
rooted layers than the soil profile has layers is an out-of-bounds read rather than
a wrong number. That used to be guaranteed by validating root carbon against
`soil_depth`; `set_root_network()` now checks it, along with finiteness and
non-negativity of both load-bearing vectors.

**One capability is lost, and it is not dressed up.** `leaf_gradient()` can no
longer differentiate with respect to `beta_R_H` or `beta_R_V` — there is no `pars`
name for them, because the leaf does not have them. Perturbing them is cheap
(`root_network_from_carbon()` is homogeneous of degree 1 in each, so the perturbed
network is a scaling of the base one and needs no rebuild) but the derivative of a
solved output still costs two solves. The previously recorded values are kept in a
comment in `tests/testthat/test-gradient.R` for anyone checking that route.

## Documented when the exact gradient beats differencing, and it is not "at more parameters"

Docs only; no code change. `vignette("fitting")` measured the exact gradient as
*slower* than differencing the objective, from a three-parameter fit, and that
result was over-generalised into "the value of `leaf_gradient()` is exactness, not
speed" — including in this package's own sources.

Both routes are linear in a parameter count, and it is a **different** count for
each:

    T_fd    = 2 * P_fit   * T_objective
    T_exact = N * (t_construct + t_solve)  +  2 * P_model * N * t_eval

An optimiser differences the `P_fit` parameters it is moving; `leaf_gradient()` is
asked for the `P_model` parameters the leaf has, i.e. `length(pars)`. Those are
equal **only when the fit varies traits directly.** Any parameterisation in between
— pooling, a hierarchy, a shared or derived parameter — makes `P_fit > P_model`, and
then the composite's cost does not follow the optimiser's dimension.

⚠️ **Growing the parameter count does not on its own reverse the verdict.** In the
vignette's design the composite's *slope* is the larger of the two, so it loses at
every count up to 13; each extra parameter costs it two more R calls per
observation. What reverses it is `P_fit` exceeding `P_model`.

Both regimes are now measured. The vignette gains a scaling sweep that fits the two
coefficients, and those coefficients — taken from 72 simulated observations —
**predict** a companion calibration study (1,327 observations, 16 species,
`P_fit = 40`, `P_model = 4`) to within a few percent: **638 ms predicted against
679 measured**, break-even 12.3 fitted parameters against 13.1. There the composite
wins **3.4×** and reaches the same optimum as the numerical gradient, and its
57-parameter variant costs the same as its 40-parameter one.

⚠️ Those two figures read 646 ms and 12.4 when first written, i.e. "within about
5%" rather than 6%. They moved because the vignette's 13-parameter candidate list
had to change: `beta_R_H` and `beta_R_V` are not differentiable parameters any more
(see the `set_physiology` entry above), so `psi_crit` and `root_psi_crit` took their
places and the fitted coefficients shifted slightly. The prediction is regenerated
on every build; the *claim* is that two designs differing 18-fold in size agree to a
few percent, and that is unchanged.

`?leaf_gradient` gains a "What it costs" section, and a warning that was missing:
**always pass `pars`**, since it *is* `P_model` and the default of all fourteen is
the most expensive request available.

Also: this file's heading said `leaf`, three renames ago.

## A gradient in `stem_b` no longer rebuilds the vulnerability spline (#4, PLAN 11f)

**24.5× faster in C++** for that parameter — 35.5 → 1.45 µs, which is the same
cost as a trait that owns no spline at all. Results are unchanged: the golden file
is bit-identical and the two routes agree to solver noise.

Perturbing `stem_b` used to rebuild the pre-integrated stem curve, and that rebuild
was essentially the whole cost of its gradient. It is unnecessary, because the
cumulative integral is homogeneous of degree 1 in `stem_b`:

    G(psi; s*b, c) = s * G(psi/s; b, c)

`stem_b` enters only as a scale on both axes, and the knot grid scales with it too,
so this holds for the **spline** and not merely for the integral it approximates —
measured, a rescaled spline reproduces a rebuilt one to 0–3e-16. So the perturbed
curve is the existing spline read at a rescaled argument.

`leaf_gradient()` gains `fast_stem_curve = TRUE`, which is how this is used;
`FALSE` rebuilds and exists so the equivalence is checked rather than assumed.
`Leaf$perturb_stem_b()` is the C++ entry point, for consumers that link the headers.

⚠️ **`stem_c` is unchanged and still rebuilds.** There is no such identity for it —
it changes the curve's shape, not its scale. The obvious alternative, reading the
curve from its closed form, was implemented and rejected: it disagrees with the
rebuild route by a systematic 3.5e-3, because it differentiates the exact integral
while the package evaluates the spline. PLAN 11f has the measurement and the rule
that follows from it.

⚠️ **From R this is 1.23×, not 24.5×**, because every `leaf_gradient()` call pays
204 µs to construct a `Leaf`. Which figure applies depends on which side of the R
boundary you are.

## `leaf_gradient()` covers the two parameters that are not traits (#44)

`pars` now accepts **`leaf_specific_conductance_max`** and — on the
single-potential path — **`resistance`**, alongside the thirteen traits.

They are here because a calibration fits them. Of that study's four free
parameters, two are traits (`cost_scale_TF24`, `beta2`) and two are these
(`K_total` and `f_plant` reach the leaf as a conductance and a resistance), so
restricting `pars` to `leaf_traits()` left half of that fit with no exact
gradient at all.

Nothing in the derivation had to change: the implicit function theorem is applied
to `dprofit/dpsi = 0`, and any parameter the profit function depends on goes
through it identically. What differs per parameter is only which setter applies
the perturbation. Both new parameters agree with a resolved central difference of
the whole solve to **8 and 9 significant figures**, and neither rebuilds a
vulnerability spline, so both sit in the fast class.

`leaf_specific_conductance_max` is also how plant's height reaches the leaf, so
this is a gradient with respect to plant *state* and not only with respect to
traits.

⚠️ **Two consequences worth knowing:**

* **The default `pars` is wider**, so `leaf_gradient()` with no `pars` now
  returns 16 rows on the multi-layer path and 17 on the single-potential one,
  rather than 15.
* **The step rule is per-parameter now.** The traits' rule floors the step at 1,
  which is right for parameters spanning 0.3 to 9.4e3 but catastrophic for
  `leaf_specific_conductance_max`, whose default is 3.14e-05 — flooring there
  would perturb it by 3% and measure a secant rather than a derivative. Those
  two get a plain relative step.

## `leaf_solve()` is 16× faster, and the fast path is now the documented one (#39)

**No results move; every value is bit-identical, including the error messages.**
This is assembly, not arithmetic.

`leaf_solve()` built a one-row `data.frame` of drivers per row, `cbind`ed an
`operating_point()` frame to it, and `rbind`ed the lot at the end. Measured over 32
rows, that cost **344 µs per row** against **2.8 µs** of solving — 96% overhead, on
the function the README and the vignette point people at. It now fills preallocated
storage and assembles once: **21.5 µs per row**.

Two smaller changes make up the rest of it:

* **`operating_point()`: 180 → 4 µs.** It read the twelve outputs through twelve
  separate calls into C++ (~1.1 µs each) and built its one row with
  `data.frame()` (158 µs). `Leaf$operating_point_values()` is new and returns all
  twelve in one call; the row is now built directly, and a test asserts the result
  is `identical()` to what `data.frame()` returned.
* **`Leaf$operating_point_values()`** is public, so C++ and Python consumers get it
  too, though there it is a convenience rather than a speedup — the cost it removes
  is the R boundary's.

⚠️ **The advice that follows from this is the opposite of what it was.**
`leaf_solve()` is now within **6%** of building a `Leaf` and driving it yourself, so
the stateful interface is for access to intermediate state and not for speed. What
is left to optimise is the *number of R calls per row* — about 18 of the 21.5 µs —
which means one vectorised call rather than a loop, and C++ rather than better R if
that is not enough. The README has a "Performance from R" section and `?leaf_solve`
a "Performance" one; both are new, and both replace claims that only ever held on
the C++ side.

## Trait gradients: `leaf_gradient()` and `set_traits()` (#4, PLAN 11e)

**No results move. The golden file is bit-identical** — this adds two entry points
rather than changing one.

**`leaf_gradient()`** returns the derivatives of the solved outputs with respect to
the traits: `dA/dθ`, `dgc/dθ`, `dψ_stem/dθ` and `dψ*/dθ`, for any subset of the
thirteen traits at one operating point. These are not finite differences of the
solve. The outputs are read at the profit-maximising collar potential, so a trait
moves them both directly and by moving that optimum, and differentiating the
optimality condition captures both terms — which matters because for
`cost_scale_TF24`, `beta2`, `stem_b` and `stem_c` the second one is **100%** of the
answer.

⚠️ **That derivation assumes the optimum is interior, and in dry soil it often is
not.** With the optimum pinned to the edge of the feasible collar range the formula
is wrong by up to seven orders of magnitude and does not fail loudly, so the
assumption is **tested at every point** — and `leaf_gradient()` differences the
solve instead where it fails. `$method` and `$status` report which route ran and
why. Across the package's 288-point grid the two populations are five orders of
magnitude apart, so the threshold is measured rather than chosen.

**`set_traits()`** replaces the traits on an existing `Leaf`, which is far cheaper
than building a new one and is the only correct way to do it: a bare trait
assignment would leave the two pre-integrated vulnerability splines, the solved
operating point, and — least visibly — `vcmax_`/`jmax_`/`R_d_` all describing the
old value. The last of those is derived behind a cache keyed on leaf temperature
and O2 alone, so "change the trait, then set the drivers again" silently does not
recompute it. Call `set_drivers()` after `set_traits()`; the object is returned to
its just-constructed state deliberately.

⚠️ **On speed, and the answer depends on which side of the R boundary you are.**
From R this was planned as a 10× speedup and is not one: per parameter the composite
is **6% slower** than differencing the solve, because an R call costs ~1.8 µs against
the 0.26 µs of C++ work it wraps. What *is* worth 4× from R is `set_traits()` —
reusing the object rather than rebuilding it — and that accrues to either method. So
**from R, use `leaf_gradient()` for exactness and the active-set classification, not
for speed.**

In C++ the picture is different, and it is the one a `LinkingTo` consumer gets:
there the composite wins **4.4×** on the eleven traits that do not touch a
vulnerability spline, and only **1.15×** on `stem_b`, `stem_c`, `root_b` and
`root_c`, where rebuilding the spline (21.8 µs) costs more than three solves.
`make -C tests/cpp bench_gradient` reproduces both arms. PLAN 11e has the numbers
and the retraction they corrected.

## The AD replicas are gone, and the model's precision floor moved (#4, PLAN 11b)

Two small results-moving changes, landed in this order because the order matters:

**`psi_stem_to_ci`'s tolerance is 1e-10, was 1e-7** — moves results by
**1.82e-07**, lands **335× closer to a converged solve**, costs **+3.4%**. This
tolerance was invisible while the collar solve's own `GSS_tol_abs` (1e-3) dominated;
once PLAN 11a removed that, it became the model's dominant amplifier and therefore
the floor of what every reported output *means*. The working magnitudes for reading
a diff are now **~1e-16 reassociation, ~1e-9 solver floor, ~1e-4 a real
difference**. Not the same knob as the settable `ci_abs_tol` (default 1e-3), which
reaches only the off-path `optimise_psi_stem_*` solvers.

**`namespace detail`'s hand-maintained AD replicas are deleted** — moves results by
**3.51e-10**. `assim_colimited`, its two components, `hydraulic_cost_TF` and
`proportion_of_conductivity` are now `T = double` instantiations of scalar-generic
member templates, and `dprofit_droot_collar_psi` differentiates the *same code*. So
the forward model and its derivative can no longer be derivatives of different
functions — which they were: the deleted replica associated the electron-limited
term differently, making the AD derivative exact for a function the model did not
evaluate.

**The `double` path is bit-identical**, checked directly across all five entry
points. The whole 3.51e-10 is the derivative changing.

Doing the tolerance first is what makes the unification a 3.51e-10 change rather
than a 4.98e-07 one — a factor of 1400 for free, purely from ordering.

To be clear about what this second change is *not*: it does **not** make the solved
operating point more accurate. Both the old and new builds locate their respective
roots to ~1e-14, and the two derivatives differ so little that their roots differ by
~1e-15 MPa. It is a maintainability and correctness-of-construction fix.

## ⚠️ The collar solve changed, and results moved (#4, PLAN 11a)

**`find_root_collar_psi` no longer maximises profit by golden-section search. It
solves the first-order condition, `dprofit/dpsi == 0`, by a safeguarded
root-find.** Every solved operating point moves by roughly the argmax correction —
worst change over the 288-point golden grid **1.5e-03** relative, 240 rows, with
`profit` itself moving at most 1.4e-03 because it is the maximum and therefore
flat. Shut-down rows are unchanged. **If you have recorded numbers from this
package, they will differ.**

This is a correctness fix, and the evidence is the residual rather than the
values: `|dprofit|` at the returned collar improved on **240 of 240** feasible
grid rows with none worse, from a median of 7.8e-04 to **5.6e-15** on the 198 rows
with an interior optimum. The remaining 42 have a *constrained* optimum pinned to
a bracket bound, where the gradient is genuinely non-zero.

Why it was worth moving results:

- **Trait derivatives were unusable, and worse than noisy.** Golden section
  resolved the argmax only to `GSS_tol_abs`, so a finite difference in a
  photosynthetic trait returned exactly **zero** below a relative step of 1e-4 —
  silently dropping a whole term of `dA/dtheta`. For traits in the hydraulic path
  it was worse: the argmax came back smooth, plausible and **sign-inverted**
  (`root_b` gave −2.6e-03 where the truth is +2.6e-04). A gradient-based
  calibration would have walked those traits the wrong way with nothing to show
  it. Finite differences now agree to ~4 digits at any relative step from 1e-8 to
  1e-2, which unblocks #6.
- **It is faster: 24.5%**, 2.65 against 3.51 µs/solve, interleaved at reps=2000.
- **The argmax got smoother, not rougher** — ~1000× smaller second differences in
  a trait. That is the opposite of what the guide's hazard 3 would lead you to
  expect, and it is measured. ⚠️ In a *trait*; smoothness in plant state, which is
  what hazard 3 is actually about, needs re-measuring on the plant side once
  plant #591 clears.

`GSS_tol_abs` is still a control and still has two jobs — the "interval too narrow
to solve over" threshold and the off-path single-layer optimisers — but it no
longer sets how well the reported operating point is determined.

Also in this change: `dprofit_droot_collar_psi` gained an optional `bool*
feasible` out-parameter, so a caller searching for a zero can tell a real
stationary point from the `0.0` it returns on its shut-down exits. The default
leaves every existing caller, including plant's TF24f, untouched.

## The package is callable from R (#5, stage 1)

`phylloptim::Leaf` now has an R interface, so this is no longer a `LinkingTo`-only
package. `library(phylloptim)` gives you an R6 `Leaf` with the drivers, the solve and
the whole operating point; `vignette("phylloptim")` walks through a solve, a drought
response and a light response. **λ and `g1_eff` are exposed for the first time**
— they existed in C++ and were unreachable from R.

**The C++ headers are unchanged, and still need no R.** The model stays a set of
self-contained headers under `inst/include` that use no R and no Rcpp; the R
layer sits on top of them and is never included by them. So a `LinkingTo: phylloptim`
consumer sees nothing new, and the model remains linkable straight into a C++
program or a Python extension. There is now a **CMake package** (`phylloptim::phylloptim`)
for exactly that, and CI builds and installs it on runners with no R, so the
claim is tested rather than asserted. See PLAN item 6a.

Two things worth knowing if you are working on this package:

- **RcppR6 is not a declared dependency.** The generated glue is committed, so
  only a developer regenerating it needs the generator. CI regenerates and diffs
  to catch a stale commit. PLAN item 6c has the reasoning, including what would
  make it worth swapping for odelia's hand-written style.
- **The golden file's bit-exactness turns out to depend on the optimisation
  level, not only on the platform**: identical at `-O1`/`-O2`/`-O3`, off by
  3.47e-15 at `-O0`, which does not contract `a*b + c` into an FMA. A debug build
  failing `test_golden` by ~1e-15 has found nothing.

## A surface you can type at a console (#5, stage 2)

The bindings above are a faithful translation of the C++, which is what let them
be checked against the golden file — and it means seventeen positional arguments.
On top of them:

- **`leaf_solve()`** — drivers in, operating point out as a data.frame,
  vectorised, so a response curve is one call. This is the entry point for
  someone who would otherwise reach for `plantecophys::Photosyn()`.
- **`leaf_traits()` and `leaf_control()`** — the split the issue asked for. Four
  of the constructor's seventeen arguments are tolerances sitting among the
  physiology, and a trait-calibration loop should not have to know which. A test
  asserts the two functions partition the constructor exactly.
- **`leaf_model()`** — named and defaulted, and the recommended constructor.
  `Leaf()` remains exported as the raw one.
- **`set_drivers()`** and **`operating_point()`** for the stateful path.
- **`vignette("phylloptim")`** — a solve, a drought response, the light/VPD surface,
  and the profit function being maximised.

Two decisions worth knowing:

- **The `psi_soil >= 0` rejection is surfaced, not smoothed over.** A friendly
  wrapper is exactly where someone would be tempted to `abs()` it, and that check
  is the only thing between a script written against the old signed convention
  and a plausible wrong number.
- **The leaf-temperature clamp is NOT in `leaf_control()`**, though the issue
  listed it. It is a guard keeping the Arrhenius block finite on the
  energy-balance path, not a tolerance anyone tunes; making it settable would let
  a caller get NaNs back with no indication why.

## A bare leaf needs no root carbon profile (#5 stage 3, #32)

`leaf_supply_singlelayer()` collapses the whole soil-to-collar path to one series
resistance, so a leaf physiologist with a soil water potential and no root-mass
profile can use the model without going through a plant-shaped one to get at a
leaf. It is also what makes the optimality-model comparison meaningful, since
Medlyn, Prentice least-cost and Cowan-Farquhar are all written against a single
soil potential.

```r
leaf_solve(psi_soil = 1.5, PPFD = 900,
           supply = leaf_supply_singlelayer(resistance = 1e3))
```

The path is chosen when the leaf is built, and **there is no
`leaf$supply_kind <- "single"`.** Flipping a tag would leave the other path's
state configured and silently ignored — and flipping back would make it stale
rather than absent. `leaf_supply_multilayer()` and `leaf_supply_singlelayer()`
reconfigure the object completely instead, so it can never be in a state where
the tag and the supply disagree. `supply_kind`, `single_resistance_` and
`single_gravity_head_` are readable but not settable, for the same reason.

On the multi-layer path `root_carbon_per_leaf_area` still has no good default and
the R layer supplies a stand-in. Taking resistances there too is #33, which is
coupled with plant.

# leaf 0.1.0

**The first version bump since the package was created, and it exists so downstream can
pin.** `0.0.1` has been the version through four merged PRs that changed the API and the
results, so a `leaf (>= 0.0.1)` requirement was satisfied by the original extraction —
which `plant` cannot compile against at all. That is the same trap odelia 0.2.1 was cut
for: a version that names several different header sets is not something a consumer can
depend on.

A minor bump rather than a patch, because the C++ API broke in ways a consumer sees at
compile time, and the numbers moved.

## Breaking API changes

- **`set_physiology` takes 10 arguments, not 14** (#15). `area_leaf`, `rho`, `a_bio` and
  `sapwood_volume_per_leaf_area` were dead stores — assigned, never read — and are gone,
  along with the four members of the same names. `mass_root_prop` became
  `root_carbon_per_leaf_area`, and it is the old value **divided by `area_leaf`**: the
  model is purely intensive, and uptake is exactly homogeneous in that ratio. Passing the
  old absolute carbon compiles and runs, and gives a root system too weak by a factor of
  `1/area_leaf`.
- **`root_collar_psi_` is now `opt_root_psi_`, and it is a positive magnitude** (#25).
  Renamed rather than reused deliberately: keeping the name with a flipped sign is the one
  outcome where an old analysis reads the wrong value in silence, and a rename gives a
  binding error instead.
- **Every ψ in the package is a positive magnitude in MPa** (#25).
  `psi_soil_inverted_`, `psi_soil_inverted_vec_` and `supply_psi_soil_inverted()` are
  deleted. `E_from_Soil_to_Root_Collar`, `find_root_psi`, `find_psi_stem_from_psi_root`,
  `dE_from_soil_dpsi_collar` and `transpiration_to_psi_stem` keep their signatures but
  take magnitudes, and `find_root_psi`'s bracket ends swap (wettest layer first). They
  validate the soil vector and stop on a negative entry, so a pre-#25 caller fails loudly
  rather than returning a wrong number.
- **`dE_from_soil_dpsi_collar` returns a positive conductance** (#25), so callers that
  negated it to recover one must stop.
- **Renames** (#15): `b`/`c` → `stem_b`/`stem_c` (there are two Weibull curves and the
  unmarked pair was the source of a real error), `g1_TF24` → `cost_scale_TF24`.
- **`umol_per_mol_to_Pa` is no longer a namespace-scope constant** (#15). It was
  `0.1013` = 101.3 kPa in disguise; it is now the member `Leaf::umol_per_mol_to_Pa_`,
  derived per call from `atm_kpa_`, so the model is self-consistent away from sea level.

## Behaviour changes

Each landed with its own measured blast radius against the golden file; see the PRs.

- **Four exits no longer leave stale state** (#15, #26). `Leaf` is reused for every
  individual in a patch, so any output a branch declined to write became the previous
  plant's value. Fixed: the three `set_shutdown_state` call sites, the `assim_max_ < 0`
  exit, `soil_consumption_` cleared with `.assign` not `.resize`, and
  `dprofit_droot_collar_psi` (which segfaulted without a prior solve and returned a
  spurious gradient in the reversed-gradient state).
- **The collar bracket is clamped to `root_psi_crit`** (#24, plant #584). The clamp
  compared a magnitude against a signed potential and could never bind. The window is
  empty at this package's defaults and 1.2 MPa wide at plant's.

## Downstream

`plant` requires `leaf (>= 0.1.0)`; anything below it will not compile against plant's
`feature/consume-leaf-package`.
