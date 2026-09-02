# Bringing reverse mode to upstream's leaf: where this is and what is left

Live on `ad/merge-upstream`, cut from `ad/v3-forward` at 33e0048. Your branch is
untouched; nothing is pushed.

## The shape of the work

It is **not a merge**. Upstream's #126 rebuilt this model -- eight optimality
curves behind one solver, a cost curve plus a benefit link, `(P50, c)` on both
vulnerability curves -- and it already does forward-mode AD through `cost_deriv<K>`
and `benefit_link<K>`. This branch's AD internals differentiate the same
quantities for TF24 alone, in odelia's spelling rather than XAD's. They are
superseded.

So: take upstream's model whole, then ADD the surface plant needs. Every file
upstream carries sits at upstream's content, including the finite-difference
gradient surface this branch had deleted -- that surface is the reference an AD
path is checked against, and #89, #90 and #131 built on it.

## Landed

| | |
|---|---|
| `1544b99` | upstream's leaf model whole; FD surface restored; nine unbuilt probes cut. Verified behaviourally identical to upstream |
| `47f7099` | all four vulnerability splines hand their closed-form slopes. Without it the spline tier moves 8.18e-05 against a 1e-10 tolerance and the R suite goes 1484/0 to 1458/36. **This is the ONLY commit that moves an operating point -- see "the golden files" below** |
| `bc91e3d` | the model owns its parameter list; `leaf_pars<S>` carries it as an indexed array, not a struct |
| `20b89b5` | the surface audit -- every coined type tested against what already exists |
| `858e182` | the coordinates (`CollarPoint : CollarCoords`, both `with_slope`), and `stem_integral_at` |
| `ca64fcc` | `VulnerabilityDerivatives` + `rows_in_P50` + `graft.hpp` |
| `3bd1b8b` | the uptake path scalar-generic; double path bit-identical |
| `551f177` | `duptake_dpsi` one body three entries; `supply_draw_at` |
| `bf3c47a` | upstream's kernels take their parameters; `collar_coords_at` with both halves of each pair |
| `6e3fd6a` | the transpose identity, over the surface that exists |
| `ab444bf` | `profit_at` / `outputs_at` / `marginal_at`. profit and the draws bit-identical to the solve's own; the marginal bit-identical to upstream's `dprofit_at_collar_psi` |
| `0b37498` | one spelling of each derivative: the residuals hold the collar, its channel is one supplied slope, and `graft_integral` takes the table's query slope |
| `dc9d4c5` | `marginal_collar_slope`, as a difference of the marginal, with the conditioning measurement that says why |
| `7648dc5` | `collar_at` / `bound_at`; the transpose identity with the collar live AND its control |
| `0acba57` | the four clamp sites counted |
| `bcc81cd` | `replay_operating_point`, `operating_point_kind_count`, `ncontrol_default`, the scalar carbon map |
| `333af55` | `test_supplied_rows` -- the one row this surface supplies |
| `bdff08d` | `outputs_at` / `marginal_at` asserted against the solve's own numbers, at every kind |
| `f2a6dc6` | the curvature falls back to one-sided where a neighbour is outside the interval |

plant carries `c49da429` (onto the rebuilt surface, and four merge breaks that
had never compiled), `6c748e14` (`with_slope` taken from odelia), `c1f4cedd`
(the tests upstream's event-queue constructor left behind, and the 0.9.0 pin),
`6f31d36b` (the pack seated from the leaf before what plant owns overrides it)
and `c6fc6958` (the three FF16 AD probes).

## What the consumer boundary cost, and what it caught

⚠️ **THE MERGE COMMIT HAD NEVER BEEN COMPILED.** Four things in plant were broken
in ways no reader would find and no test could reach: a rename that missed one
site, a lost `template <typename T, typename E>` line, a state read into a
`double` inside a templated environment, and a call carrying an argument this
branch's design had removed. A syntax-only sweep over all 27 translation units is
what surfaced them, and it is the first thing to run after any merge here.

Three more the SUITE caught, none of them about the model:

* upstream's event queue gave `SCM` an `events` argument, and eleven
  three-argument calls across six test files bound `Control()` to `events`.
* the FF16 AD probes named odelia's scalars without including odelia, compiled
  below C++20 where its concepts do not parse, and linked its tape without the
  two XAD storage-class defines -- the hazard that names neither package.
* the parameter pack was filled with zeros and then thirteen of nineteen slots
  set, leaving seven cost curves' parameters at zero rather than the leaf's
  values. plant's gradient ladder named it, by asking that every trait either has
  a column or is refused BY NAME.

Odelia carries `af8b1e4` (compat interpolator, v0.4.0), `90f0dc8`
(RECORDED-DECISIONS.md) and `f9c06de` (`with_slope` moved in). plant carries
`5bde7727`, the merge on your tip.

## What is left

1. **The golden re-bless, on macOS/arm64.** Cannot be done from here -- see "the
   golden files" below.

2. **plant: `k_I` HAS NO ROW, AND THE MODEL DEPENDS ON IT.** ⚠️ A real gradient
   gap, found by plant's own ladder and confirmed against the model.

   The light extinction coefficient's column comes back an EXACT zero, and an
   exact zero is the signature of a missing accumulator rather than of
   insensitivity. Measured: perturbing `k_I` moves the short stand's census at
   **d/dk_I = 16.21**, stable across two step sizes.

   It is read at the active scalar in two places -- `compute_competition_and_slope`
   and `radiation_at` -- so the parameter IS registered and the question is where
   the light field is built. It is plant's competition path, not the leaf.

   Six ladder assertions hang off it: rung3's "every trait the block reads has a
   column", floor's two classification checks, rung4's "no shortlisted trait reads
   exactly zero", declared-zero's, and rung5's newcomer check.

3. **plant: the finite-difference comparisons are over budget, and the AD is the
   side that is right.** Every AD-INTERNAL check passes -- the block Jacobian
   forward against reverse at 2.01e-15, both transposes at 1e-16, the state
   Jacobian along a trajectory at 3.90e-16. Only comparisons against a DIFFERENCE
   of the model fail: the whole Jacobian at 2.25e-04 against 8.08e-07, the
   right-hand-side transpose at 5.78e-04 against 8.82e-06.

   That is this branch's stem-vulnerability finding at stand scale, and the entry
   below says why the AD is right and what closing it costs. It is the same
   decision, arriving where the user will feel it.

4. **`test-stochastic-patch-runner`'s baseline**, which moved on BOTH sides and so
   needs measuring rather than merging. It is the one file whose PASS count varies
   run to run.

5. **TF24f's acclimation rate does not respond to the price.** Everything else
   about TF24_floor does: the solve optimises it, `shadow_cost` reports it, and
   the stomata close. TF24f evaluates at a TRACKED collar rather than optimising
   one, and its tracked state's gradient comes back the same at lambda_o = 0 and
   1e5 -- both about 2e-13. Every phylloptim entry point it calls now names
   TF24_floor, so the question is narrower than the seating: it is about that
   strategy's own differenced gradient.

⚠️ **`test-strategy-ff16 :: Report generation` fails for want of `kableExtra`,
which is an environment gap rather than a defect** -- it should be behind a
`skip_if_not_installed`. Recorded so the next reader does not chase it.

## The gate that has caught everything

**Bit-identity at `T = double`, against upstream's own function.** Not a spot
check -- it is what caught the dropped `kg_per_mol_h2o`, where every per-layer
value stayed exact and only the aggregate moved. A layer-wise check passed. The
tell was the ratio: a clean 1/0.018015 is a missing constant, not drift.

Then, for the active path, the row against whatever upstream computes the same
quantity by -- and where nothing does, the transpose identity, which needs no
reference at all.

⚠️ **COMPARE THE QUANTITY THAT CARRIES THE INFORMATION.** `dsigma/dcollar` is
1 + 9e-6, and the informative part is that 9e-6 -- so a check on the whole number
agrees to 7 digits while the part that matters was 8 per cent wrong.

⚠️ **AN EARLIER REVISION SAID A DIFFERENCE OF THE MODEL COULD NOT REFEREE THESE
SLOPES. IT CAN, AND IT WAS RIGHT WHILE THE ASSEMBLY WAS WRONG.** See
SURFACE-AUDIT.md; the short version is that the solve ran on the TABLE, so the
table's value is the one every partial takes and only the rows come from the
curve. `marginal_at` is bit-identical to upstream's `dprofit_at_collar_psi` now,
and the slopes agree with a difference at the shipped resolution.

## The golden files, correctly attributed

⚠️ **This section was WRONG in an earlier revision, which claimed "32
attributable rows against a CLEAN baseline". The instrument was wrong: it read
the macOS-generated file on Linux, where the whole grid disagrees for reasons
that have nothing to do with this branch.** The measurement below holds odelia
constant and compares our own commits, then checks the baseline separately.

* **Pristine `upstream/master` PASSES `--cross-platform` on Linux, exit 0**, with
  the argmax class at 9.02e-06 against a 5e-03 tolerance -- the figures the
  package's own guide records.
* **This branch gives 147 mismatched values**, argmax class worst 0.0609. That is
  a real, attributable move, not a platform artefact.
* **All of it is `47f7099`.** Regenerating the grid from each commit against one
  odelia: `1544b99` to `47f7099` moves 480 of 576 rows; `47f7099` to HEAD is
  **bit-identical in all 576 x 9**. Every surface commit since is inert at
  double, which is the property the gate exists to keep.
* **The size is uniform, about 1.2e-06 of each column's own scale.** The 0.0609
  is a relative figure at a near-zero row (`assim` = 1.9e-07); read the absolute
  alongside it, as the guide says.

**So the golden files need a deliberate re-bless on macOS/arm64**, with the
measured blast radius in the PR, because a different interpolant is a different
operating point. It cannot be done from here. Until then `make` fails on
`test_golden` and `test_primitives`, and that failure is EXPECTED and attributed
-- which is exactly the state the guide warns rots into being ignored, so it
should not sit here long.

## Standing hazards

* **Every comparison on the active path reads passive.** `std::min` on two active
  values branches on a taped comparison and manufactures a discontinuity.
* **Units are mixed by design.** `E_up` in kg, per-layer draws in mol, converted
  downstream in plant's `compute_rates`. It lives in a comment because it is
  upstream's convention and plant depends on it.
* **The rows are in `(P50, c)`.** `graft.hpp` carries the chain; a row in `b` or
  `psi_crit` reaches nothing, because `set_traits` re-derives them.
* **THE TWO FLUXES ARE DIFFERENT QUANTITIES AND EACH HAS ITS OWN READER.** The
  stem flux `kmax*(G(sigma)-G(collar))` and the soil draw `E_up` agree only to
  the stem splines' round-trip, measured 6.7e-04 apart at a converged interior
  point. `gc`, and so ci, is built from the STEM one; the inverse that PLACED
  sigma is read at the SOIL one. Swapping either is finite and plausible.
* **Every partial takes the TABLE's conductivity with the CURVE's rows**, through
  `graft_curve`, because the solve ran on the table. The closed form's value
  there makes M non-zero at the collar the solve placed.
* **A derivative must have ONE spelling.** Where a quantity is both lifted by a
  residual and carried as an explicit slope, the two are different numbers
  whenever the table and its inverse are not exact mutual inverses -- and both
  are finite, plausible, and wrong in different places. The residuals hold the
  collar for exactly this reason.
* **A nested active is built by assigning into its value.** `nested(inner)`
  strips the inner rows where it compiles at all, giving a slope that is right at
  double and zero in every trait.
