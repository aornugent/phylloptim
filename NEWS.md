# phylloptim 0.2.1

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

* **`leaf_supply_single()` no longer takes `resistance`.** Migration:
  `leaf_supply_single(resistance = r)` ->
  `leaf_supply_single()` plus `root_network = series_resistance(r)` on
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
**predict** the companion study `leaf-calibration` (1,327 observations, 16 species,
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

They are here because a calibration fits them. Of `leaf-calibration`'s four free
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

`leaf_supply_single()` collapses the whole soil-to-collar path to one series
resistance, so a leaf physiologist with a soil water potential and no root-mass
profile can use the model without going through a plant-shaped one to get at a
leaf. It is also what makes the optimality-model comparison meaningful, since
Medlyn, Prentice least-cost and Cowan-Farquhar are all written against a single
soil potential.

```r
leaf_solve(psi_soil = 1.5, PPFD = 900,
           supply = leaf_supply_single(resistance = 1e3))
```

The path is chosen when the leaf is built, and **there is no
`leaf$supply_kind <- "single"`.** Flipping a tag would leave the other path's
state configured and silently ignored — and flipping back would make it stale
rather than absent. `leaf_supply_multilayer()` and `leaf_supply_single()`
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
