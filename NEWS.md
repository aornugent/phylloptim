# phylloptim 0.27.4

- For C++ users, each cost curve is one definition, `phylloptim::Curve<K>` in `leaf_cost_curves.hpp`: its name, the parameters it owns and their defaults, its properties, and its cost, the cost's slopes and its marginal cost of water. A curve missing a piece is a compile error naming it. No result changes and every bound method keeps its name ([#226](https://github.com/traitecoevo/phylloptim/issues/226)).

# phylloptim 0.27.3

- `$find_root_collar_psi()` solves the curve the leaf is set to, on the collar route, and refuses what that curve refuses (TF24_floor without its price, for one); it used to solve TF24 whatever curve was set. A leaf on the default model is unchanged.
- `leaf_gradient(model = "collar")` on a reused leaf set to another curve differentiates TF24's collar route, as a fresh leaf does, at a solved or a prescribed `psi`. It used to solve TF24 and test the point against the other curve, reporting it pinned, and to differentiate the other curve at a prescribed `psi`. The leaf is left set to the model the gradient used, as the stem routes already left it, and solved on that route.
- `leaf_gradient_batch()` uses a `step` you name wherever it takes a finite difference (`method = "fd"`, or a pinned row under `"auto"`), as `leaf_gradient()` does; it used to ignore it there ([#219](https://github.com/traitecoevo/phylloptim/issues/219)).
- For C++ users, the derivative surface (`operating_collar()`, `billed_outputs_at()`) refuses a leaf with the energy balance on. It held leaf temperature at the solve's value, so its rows differed from a difference of the solve, by up to the whole row ([#224](https://github.com/traitecoevo/phylloptim/pull/224)).

# phylloptim 0.27.2

- The leaf solve is about 11% faster, and plant's seasonal stand about 15%: the solver no longer recomputes the roots' water uptake three times for each candidate. Results do not change ([#214](https://github.com/traitecoevo/phylloptim/issues/214)).

# phylloptim 0.27.1

- Evaluating a leaf at a stem potential you choose (`evaluate_psi_stem_at()`) now reports that point's marginal cost of water (`lambda_emergent`); it used to report the value left over from the previous solve. Nothing else changes.
- `leaf_gradient_batch()` now refuses prices the way `leaf_solve()` does: a price the chosen model does not use is an error, and so is `model = "TF24_floor"` without `TF24_floor_lambda_o`, instead of an error in every row.
- `leaf_solve()` takes one price for all rows and refuses a vector; for a price per observation, use `leaf_batch()`, which now accepts one ([#207](https://github.com/traitecoevo/phylloptim/issues/207)).

# phylloptim 0.24.0 – 0.27.0

- **0.27.0**: the solver now decides reliably whether a leaf's best opening is at the dry end of its range or just short of it; the decision used to flip under tiny changes of input, so results jumped. Only such leaves move, by about one part in a billion in profit, but a finite difference over a whole plant run notices: plant must recapture its reference gradients (its drought and seasonal sweeps move by about 1%), and its seasonal stand runs 2% slower ([#205](https://github.com/traitecoevo/phylloptim/pull/205)).
- **0.26.2**: a faster solve for the leaf's internal CO₂: the leaf solve is about 21% faster and plant's TF24 stands about 22%; results move only in the last few digits ([#203](https://github.com/traitecoevo/phylloptim/pull/203)).
- **0.26.0 – 0.26.1**: under SOX and JW26 a leaf that cannot fix carbon at any opening (in the dark, say) is reported closed, with no water flow, instead of fully open or refused; a leaf that can photosynthesise is unchanged. For C++ users, the derivatives handed to plant for a closed leaf whose roots outlast its stem are corrected; nothing returned to R changes ([#201](https://github.com/traitecoevo/phylloptim/pull/201), [#202](https://github.com/traitecoevo/phylloptim/pull/202)).
- **0.25.0**: water flow near the dry end of the vulnerability curves is computed accurately where rounding used to lose it, so a leaf no longer shuts early at about 12.4 MPa and roots keep conducting to the end of their curve (15.4 MPa at the default traits, was 14.5). For C++ users, `Leaf::stem_capacity_floor` is removed ([#199](https://github.com/traitecoevo/phylloptim/pull/199)).
- **0.24.0**: with the energy balance on, CO₂ crosses the leaf's boundary layer at that layer's own diffusion rate rather than the stomata's, so photosynthesis rises slightly (a median 0.3–0.8%, at most 1.4%). With the energy balance off nothing changes ([#197](https://github.com/traitecoevo/phylloptim/pull/197)).

Full notes: `git show d972042:NEWS.md`.

# phylloptim 0.21.0 – 0.23.1

- **0.23.0 – 0.23.1**: for C++ users, the derivative functions plant uses are renamed (for example `collar_at` is `operating_collar` and `outputs_at` is `billed_outputs_at`; the full list is in the PRs). Results and the R interface are unchanged ([#195](https://github.com/traitecoevo/phylloptim/pull/195), [#198](https://github.com/traitecoevo/phylloptim/pull/198)).
- **0.22.1 – 0.22.2**: with a single soil potential, a very small series resistance (10 or less) on a dry soil no longer makes the solver report a leaf closed that should be open (its profit was up to 1.14 too low); a resistance below 0.01 is refused ([#194](https://github.com/traitecoevo/phylloptim/pull/194), [#196](https://github.com/traitecoevo/phylloptim/pull/196)).
- **0.22.0**: inputs outside the model's range are refused instead of silently clamped. `set_drivers()` refuses a vapour-pressure deficit at or below zero or above saturation, and drivers that would heat a leaf past what the model can evaluate; the old clamp of leaf temperature to −40–70 °C, with `leaf_temp_min` and `leaf_temp_max`, is gone. A stem potential outside the feasible range and a negative damage state are refused too, and `leaf_gradient()` reports `"refused"` where it said `"clamped"`. Results change only near the old clamps ([#191](https://github.com/traitecoevo/phylloptim/pull/191)).
- **0.21.0**: for TF24, TF24_floor, ProfitMax, JS22 and CMax the search for the best opening now runs to the end of the stem's vulnerability curve instead of stopping at its critical potential `psi_crit`, so a leaf that was held there can open further; SOX, JW26, CF77 and LeastCost still stop at `psi_crit`. The statuses `"boundary-crit"` and `"boundary-root-crit"` become `"boundary-dry"`, and the integer codes of later statuses shift. For C++ users, `supply_psi_crit`, `root_psi_crit_at`, `invert_stem_curve*`, `eval_stem_curve` and `closed_form::PastLastKnot` are removed ([#188](https://github.com/traitecoevo/phylloptim/pull/188)).
- **0.21.0**: the vulnerability curves are read from a smoother, more accurate table that a change of P50 no longer rebuilds; values move only far down the significant digits, and the leaf solve is 1–4% slower (plant's TF24 stands about 11%) ([#188](https://github.com/traitecoevo/phylloptim/pull/188)).

Full notes: `git show v0.26.2:NEWS.md`.

# phylloptim 0.17.0 – 0.20.0

- **0.20.0**: models with a level parameter now default to the value that gives TF24's stomatal conductance at one reference environment (`align_curves()`): `JS22_gamma` 1.0 → 0.4718, `CMax_a` 0.6 → 0.6386, and `CF77_lambda_` and `LeastCost_beta`, which had no default, now 107751 and 51.60. Set the old values to reproduce earlier JS22 and CMax results; TF24's results do not change. The documentation is now a website ([#179](https://github.com/traitecoevo/phylloptim/pull/179)).
- **0.19.0**: a new trait, `g_min`, the cuticular conductance: the water a leaf loses through its surface with the stomata shut (default 0; `leaf_traits(g_min = )`). `leaf_solve()` and `operating_point()` gain a last column, `residual_flux`. For C++ users, `set_traits()` takes `g_min` as its sixteenth argument and the parameter numbering shifts (`n_theta` is 35). At `g_min = 0` an open leaf is unchanged ([#185](https://github.com/traitecoevo/phylloptim/pull/185)).
- **0.19.0**: a leaf too dry, too dark or too hot to open is reported as a closed leaf with ordinary values, never an error. The status `"shade-death"` is removed (later codes shift), as are `set_shutdown_state()`, `ci_assigned()`, the `no_flow_*()` functions and `root_limit_binds_first()`; closed leaves' values change ([#185](https://github.com/traitecoevo/phylloptim/pull/185)).
- **0.18.0**: thermal damage follows Bison & Michaletz (2026), every rate with the same form of temperature response. Removed: `impair_H_`, `impair_Tmid_`, `k_i_Tref_`, `m_rep_`, `t_rep_cut_` and `$logistic_()`; added: `k_i_Tiso_`, `k_rec_Ea_` and `k_mat_Ea_`; `k_i_` and `k_i_baseline_` default to 276 per day (were 0.5). The damage layer is off by default, so nothing else changes ([#186](https://github.com/traitecoevo/phylloptim/pull/186)).
- **0.17.0**: profit no longer jumps as a leaf closes, because respiration no longer switches off abruptly at shutdown or in the dark. The status `"determined"` is removed (later codes shift) and `respiration_temp_deriv()` is `temp_deriv()`; only closed leaves change ([#177](https://github.com/traitecoevo/phylloptim/pull/177)).

Full notes: `git show v0.26.2:NEWS.md`.

# phylloptim 0.16.0 – 0.16.4

- **0.16.4**: `transpiration()`, `stom_cond_CO2()` and `psi_stem_to_ci()` refuse a stem wetter than the soil instead of returning a negative flow. The solver never asks for one, so results are unchanged ([#181](https://github.com/traitecoevo/phylloptim/pull/181)).
- **0.16.3**: `model_registry()` lists every equation and parameter with its units, source and the code that evaluates it ([#168](https://github.com/traitecoevo/phylloptim/pull/168)).
- **0.16.2**: the Medlyn route reported its conductance to water vapour as its conductance to CO₂, 1.6 times too high; it now reports CO₂ (`ci` and `A` are unchanged), its soil-moisture factor stays between 0 and 1, and a leaf at zero deficit no longer returns `Inf` ([#130](https://github.com/traitecoevo/phylloptim/pull/130)).
- **0.16.1**: the leaf's internal CO₂ is solved more tightly, removing small numerical noise that could distort gradients; profit moves by at most 3e-8, and a solve costs about 5% more ([#173](https://github.com/traitecoevo/phylloptim/pull/173)).
- **0.16.0**: requires odelia 0.10.0 or later. For C++ users, odelia's reverse-mode names change (`adjoint_scalar`, `supply_derivatives`, `with_slope<T>`, `set_direction`), and `Leaf::pair<T>` is replaced by `odelia::with_slope<T>`. Results are unchanged ([#164](https://github.com/traitecoevo/phylloptim/pull/164)).

Full notes: `git show v0.26.2:NEWS.md`.

# phylloptim 0.15.0 – 0.15.4

- **0.15.4**: on the collar route, `lambda_emergent` is the chosen model's marginal cost of water; it was always TF24's (wrong by up to 146% for CF77). The unused fields `count`, `psi_stem`, `opt_ci_` and `f_r` are removed from R. With the energy balance on, stem optima move very slightly ([#171](https://github.com/traitecoevo/phylloptim/pull/171)).
- **0.15.3**: gradients are available wherever the solve returns an answer, including very dry soils. For C++ users, `closed_form_integral()` and `closed_form_curve()` gain an argument, removed again in 0.21.0 ([#160](https://github.com/traitecoevo/phylloptim/pull/160)).
- **0.15.2**: new outputs `gb_H2O` (the boundary layer's conductance) and `gs_stomatal_H2O` (the stomata's own); `gs_H2O` is still the total. With the energy balance on, LeastCost's `lambda_emergent` changes slightly ([#169](https://github.com/traitecoevo/phylloptim/pull/169), [#170](https://github.com/traitecoevo/phylloptim/pull/170)).
- **0.15.1**: with the energy balance on, SOX, JW26 and ProfitMax find their true optimum, which could be 0.1 MPa or more away at high light, and the fast closed-form method refuses the energy balance, under which it overstated photosynthesis by 15–178%. With the energy balance off nothing changes ([#121](https://github.com/traitecoevo/phylloptim/pull/121), [#167](https://github.com/traitecoevo/phylloptim/pull/167)).
- **0.15.0**: a ninth model, `LeastCost` (Prentice et al. 2014, which minimises the cost per unit of carbon gained), with a dimensionless `LeastCost_beta`. For C++ users the parameter numbering shifts (`n_theta` is 35), so a hand-written `theta` needs the extra entry ([#137](https://github.com/traitecoevo/phylloptim/pull/137)).

Full notes: `git show v0.15.4:NEWS.md`.

# phylloptim 0.12.0 – 0.14.1

- **0.14.1**: with a single soil potential, a soil drier than the stem can draw from gives a closed leaf instead of an error ([#158](https://github.com/traitecoevo/phylloptim/pull/158)).
- **0.14.0**: with the energy balance on, the longwave radiation balance follows air temperature and humidity instead of a fixed −40 W m⁻², so a leaf below about 40 °C is 0.6–2.8 K cooler and its profit 22–25% higher. `longwave_net_offset` is removed; if you overwrite `ra_` with a combined convective and radiative resistance, set `g_rad_` to 0 as well ([#120](https://github.com/traitecoevo/phylloptim/pull/120)).
- **0.13.0**: the trait `a` is renamed `phi_J`, the quantum yield per absorbed photon: `leaf_traits(a = x)` becomes `leaf_traits(phi_J = x / 0.85)`, and likewise in `Leaf()`, `set_traits()`, `l$a` and `leaf_gradient(pars = "a")`. New absorptances `alpha_PAR_` (0.85) and `alpha_NIR_` (0.20), and `PPFD` is the light the leaf intercepts. At the defaults nothing changes ([#162](https://github.com/traitecoevo/phylloptim/pull/162)).
- **0.13.0**: with the energy balance on, set `l$NIR_int_`, the near-infrared the leaf intercepts (about `PPFD / 4.57` in daylight, 0 under LEDs), or `set_physiology()` refuses. The leaf absorbs less radiation (345 W m⁻² instead of 657 at PPFD 1500) and runs a median 7.9 °C cooler. `sw_abs_per_par` is removed ([#162](https://github.com/traitecoevo/phylloptim/pull/162)).
- **0.12.0**: an optional thermal-damage layer (`use_thermal_damage_`, off by default) tracks how much of the leaf is functional, impaired and damaged as its temperature varies, with `$damage_rates()`, `leaf_damage_trajectory()` and `leaf_damage_gradient()`. For C++ users its inputs join the parameter numbering, so a hand-filled `theta` must recount. Off, nothing changes ([#162](https://github.com/traitecoevo/phylloptim/pull/162)).

Full notes: `git show v0.15.4:NEWS.md`.

# phylloptim 0.10.0 – 0.10.3

- **0.10.1 – 0.10.3**: odelia is a minimum version (`odelia (>= 0.8.0)`) rather than a pinned one, so a package using both pins odelia once; no code changes ([#153](https://github.com/traitecoevo/phylloptim/pull/153), [#154](https://github.com/traitecoevo/phylloptim/pull/154), [#155](https://github.com/traitecoevo/phylloptim/pull/155)).
- **0.10.0**: the leaf hands a model that contains it, such as plant, the exact derivatives of its solved state (`vignette("derivative-surface")`); requires odelia 0.7.0 or later ([#136](https://github.com/traitecoevo/phylloptim/pull/136)).
- **0.10.0**: the vulnerability curves are tabulated more finely (`vulnerability_curve_ncontrol` is 400, was 100), so results move slightly; a leaf held at a bound is placed more precisely, and `leaf_gradient(method = "ift")` refuses one by name ([#136](https://github.com/traitecoevo/phylloptim/pull/136)).
- **0.10.0**: `evaluate_root_collar_psi` and `dprofit_droot_collar_psi` use the chosen model rather than always TF24, and a replayed point hands plant the same derivatives its solve did ([#136](https://github.com/traitecoevo/phylloptim/pull/136)).

Full notes: `git show v0.10.3:NEWS.md`.

# phylloptim 0.6.0 – 0.9.2

- **0.9.0 – 0.9.2**: under CF77, TF24_floor, JS22 and CMax the collar solve maximised the wrong objective; their optima move by up to 0.35 MPa, each to a higher profit. Soil layers can differ in thickness: for C++ users `root_network_from_carbon()` takes a vector of widths and `layer_thickness()` is `layer_thicknesses()`; a graded profile had up to 3.7 times the root resistance it should ([#134](https://github.com/traitecoevo/phylloptim/pull/134), [#135](https://github.com/traitecoevo/phylloptim/pull/135), [#141](https://github.com/traitecoevo/phylloptim/pull/141), [#144](https://github.com/traitecoevo/phylloptim/pull/144)).
- **0.8.0**: `leaf_solve()` reports the chosen model's `lambda_emergent` (its `lambda` column is always TF24's), and a reused leaf no longer carries a price over from a closed leaf before it ([#133](https://github.com/traitecoevo/phylloptim/pull/133)).
- **0.7.0**: a new model, `TF24_floor`, TF24 plus a price of water `TF24_floor_lambda_o`; `leaf_solve()` and `leaf_batch()` take the prices; `set_model()` takes a method, where `"closed"` is a fast approximation for TF24 and CF77; and CF77 can close as the soil dries, set through `theta`, `theta_fc` and `theta_w` ([#128](https://github.com/traitecoevo/phylloptim/pull/128), [#131](https://github.com/traitecoevo/phylloptim/pull/131), [#132](https://github.com/traitecoevo/phylloptim/pull/132)).
- **0.6.0**: seven optimality models behind one solver: choose one with `l$set_model(curve, route)` and solve with `l$optimise()`, which replace the per-model functions; `leaf_solve()` and `leaf_gradient()` take `model =`; every model reports `lambda_emergent`; `plot_model_overview()` draws the model ([#126](https://github.com/traitecoevo/phylloptim/pull/126)).
- **0.6.0** renames: `leaf_supply_single()` is `leaf_supply_singlelayer()`; the vulnerability traits are `stem_P50`, `stem_c`, `root_P50` and `root_c` (convert a published curve with `weibull_p50_c()`), with `stem_b`, `psi_crit`, `root_b` and `root_psi_crit` derived; `beta2`, `cost_scale_TF24` and `lambda_` are `TF24_beta2`, `TF24_cost_scale` and `CF77_lambda_`; `evaluate_psi_stem_by()` and `dprofit_dpsi_stem_by()` are `evaluate_psi_stem_at()` and `dprofit_dpsi_stem_checked()`; `perturb_stem_b()` is `perturb_stem_P50()`; the statuses `PinnedWet`/`PinnedDry` are `"boundary-soil"`/`"boundary-crit"`. Removed: the prescribed-lambda Sperry optimiser, `profitmax_scan_n_` and `leaf_behaviour_fingerprint()` ([#126](https://github.com/traitecoevo/phylloptim/pull/126)).

Full notes: `git show v0.10.3:NEWS.md`.

# phylloptim 0.1.0 – 0.5.2

- **0.5.0 – 0.5.2**: `lambda_` is an input that nothing overwrites; a closed leaf reports one leaf temperature, which changes energy-balance results at shutdown; the unused field `lambda_analytical_` is removed ([#111](https://github.com/traitecoevo/phylloptim/pull/111), [#117](https://github.com/traitecoevo/phylloptim/pull/117), [#118](https://github.com/traitecoevo/phylloptim/pull/118)).
- **0.4.0**: an infeasible input raises a condition class you can catch, rather than a message to match; `Tleaf` and `gs_H2O` are outputs and the H₂O:CO₂ diffusion ratio can be set; the exact gas constant and molar mass of water move results slightly ([#101](https://github.com/traitecoevo/phylloptim/pull/101), [#102](https://github.com/traitecoevo/phylloptim/pull/102), [#103](https://github.com/traitecoevo/phylloptim/pull/103), [#106](https://github.com/traitecoevo/phylloptim/pull/106), [#108](https://github.com/traitecoevo/phylloptim/pull/108)).
- **0.3.0**: `R_d_25` is a trait and respiration rises with temperature (`rd_to_vcmax_ratio_` is removed); diffusion is driven by the leaf-to-air vapour-pressure deficit, which lowers `gs` and `A` with the energy balance on; the vulnerability curves cover their full range, which moves results; Sperry's ProfitMax model; trait gradients; traits can be read back ([#86](https://github.com/traitecoevo/phylloptim/pull/86), [#89](https://github.com/traitecoevo/phylloptim/pull/89), [#90](https://github.com/traitecoevo/phylloptim/pull/90), [#93](https://github.com/traitecoevo/phylloptim/pull/93), [#99](https://github.com/traitecoevo/phylloptim/pull/99)).
- **0.2.0**: the package is phylloptim and callable from R (`Leaf`, `leaf_solve()`, `leaf_supply_singlelayer()`, `leaf_gradient()`, `leaf_gradient_batch()`); `set_physiology()` takes root resistances instead of root carbon; the solver finds the optimum more precisely, which moves results ([#30](https://github.com/traitecoevo/phylloptim/pull/30), [#36](https://github.com/traitecoevo/phylloptim/pull/36), [#47](https://github.com/traitecoevo/phylloptim/pull/47), [#63](https://github.com/traitecoevo/phylloptim/pull/63)).
- **leaf 0.1.0**: the first versioned release, for plant to pin. `set_physiology()` takes 10 arguments, not 14; every water potential is a positive number (a magnitude), so `root_collar_psi_` is `opt_root_psi_` and `dE_from_soil_dpsi_collar` returns a positive conductance; `mass_root_prop` is `root_carbon_per_leaf_area`, per unit leaf area; `b` and `c` are `stem_b` and `stem_c`, and `g1_TF24` is `cost_scale_TF24` ([#29](https://github.com/traitecoevo/phylloptim/pull/29)).

Full notes: `git show v0.10.3:NEWS.md`.
