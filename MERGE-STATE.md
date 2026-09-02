# Bringing reverse mode to upstream's leaf: where this is and what is left

Live on `ad/merge-upstream`, cut from `ad/v3-forward` at 33e0048. Your branch is
untouched; nothing is pushed.

## The shape of the work

It is **not a merge**. Upstream's #126 rebuilt this model -- seven optimality
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
| `1544b99` | upstream's leaf model whole; FD surface restored; nine unbuilt probes cut. Verified behaviourally identical to upstream: test_leaf 2168/0, test_golden exit 0 |
| `47f7099` | all four vulnerability splines hand their closed-form slopes. Without it the spline tier moves 8.18e-05 against a 1e-10 tolerance and the R suite goes 1484/0 to 1458/36 |
| `bc91e3d` | the model owns its parameter list; `leaf_pars<S>` carries it as an indexed array, not a struct |
| `20b89b5` | the surface audit -- every coined type tested against what already exists |
| `858e182` | the coordinates (`CollarPoint : CollarCoords`, both `with_slope`), and `stem_integral_at` |
| `ca64fcc` | `VulnerabilityDerivatives` + `rows_in_P50` + `graft.hpp` |
| `3bd1b8b` | the uptake path scalar-generic; double path bit-identical |

Odelia carries `af8b1e4` (compat interpolator, v0.4.0), `90f0dc8`
(RECORDED-DECISIONS.md) and `f9c06de` (`with_slope` moved in). plant carries
`5bde7727`, the merge on your tip.

## What is left, in dependency order

1. **`supply_draw_at`** -- wraps `uptake_at<T>` and `duptake_dpsi`. The per-layer
   collar slopes stay SUPPLIED at double: upstream already has that row, and
   taking the supply again at a live collar records the whole of it twice.
2. **`collar_coords_at`** -- lifts sigma and ci with `implicit_value` on the T1
   and T2 residuals, and fills BOTH halves of each pair. Takes the `SupplyDraw`,
   not a bare flux, so a flux and a slope from different collars cannot be paired.
3. **`test_transpose`** -- port it HERE, not at the end. It needs no reference and
   catches a J that is not self-consistent, which is the first thing that can go
   wrong once (1) and (2) compose.
4. **`profit_at<K,S>`, `outputs_at`, `marginal_assembled<K,S>`** (returning a
   scalar), **`marginal_at`**, **`marginal_collar_slope`**.
5. **`test_supplied_rows`** -- the two rows no difference of the recorded step can
   referee.
6. **`bound_at`** at TWO arms, **`collar_at`** taking the curvature per 33e0048,
   `clamp_sites.hpp`, `clamp_count`, `operating_point_kind_count`.
7. **plant** -- call sites onto the new signatures; the four derived quantities
   stop being passed at all, which is what finally makes the `no_gradient`
   question moot rather than deferred. Also: plant still defines its own
   `with_slope.h` with 11 users, now duplicating odelia's.
8. **The R instruments** -- five files, all addressing the leaf through retired
   trait names.
9. **Verification and re-blessing** -- plant's suites; test_golden's 32
   attributable rows; `test-stochastic-patch-runner`'s baseline, which moved on
   BOTH sides and so needs measuring rather than merging.

## The gate that has caught everything

**Bit-identity at `T = double`, against upstream's own function.** Not a spot
check -- it is what caught the dropped `kg_per_mol_h2o`, where every per-layer
value stayed exact and only the aggregate moved. A layer-wise check passed. The
tell was the ratio: a clean 1/0.018015 is a missing constant, not drift.

Then, for the active path, the row against whatever upstream computes the same
quantity by -- `duptake_dpsi` for the collar, differenced `set_traits` for the
trait rows.

## Standing hazards

* **Every comparison on the active path reads passive.** `std::min` on two active
  values branches on a taped comparison and manufactures a discontinuity.
* **Units are mixed by design.** `E_up` in kg, per-layer draws in mol, converted
  downstream in plant's `compute_rates`. It lives in a comment because it is
  upstream's convention and plant depends on it.
* **The rows are in `(P50, c)`.** `graft.hpp` carries the chain; a row in `b` or
  `psi_crit` reaches nothing, because `set_traits` re-derives them.
* **test_golden carries 32 attributable rows** against a CLEAN baseline -- a
  different interpolant is a different operating point. Worst 0.061, confined to
  `collar/multi3` where fluxes run ~1e-12. A re-bless with a named mechanism, and
  not, as an earlier note here wrongly claimed, fewer failures than the control.
