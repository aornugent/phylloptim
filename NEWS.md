# phylloptim 0.9.0

## ⚠️ Breaking: the leaf supplies its own derivative rows, and three counts moved

`Leaf` is templated on its scalar and supplies the derivative of its operating
point from the implicit function theorem at the converged point, instead of
being recorded and differentiated. What that changes for a caller:

* **`vulnerability_curve_ncontrol` defaults to 400, was 100.** The
  pre-integrated vulnerability tables are built on four times the knots, and
  every number read off them moves. `leaf_control()` carries the new default.
  ⚠️ The same number is written in three places -- `Leaf::ncontrol_default` in
  `leaf_model.hpp`, `.leaf_control_defaults` and `leaf_control()`'s formals --
  and nothing checks them against each other.
* **`n_pars` is 20, and the fitted length is `n_theta` = 19.** The pack the
  kernels read also holds `par_PPFD`, which is seated per observation from the
  drivers and is never fitted, so it has a slot and no name. Read `theta` out to
  `n_pars` and you run past its end.
* **odelia `>= 0.5.0`** is required, for the reverse-mode surface and the
  interpolator this leaf's supplied rows are read through. ⚠️ odelia 0.5.0
  retracts its own 0.2.2 promise that an out-of-domain interpolator read is
  refused with a located message: `hermite_interpolator` extends linearly from
  the end knot and there is no switch. Each of this package's three curves now
  states its own bound where it applies it -- `root_vuln_at` clamps the
  argument, `root_vuln_integral_at` caps the value, `eval_stem_curve` raises.

⚠️ **This retracts two earlier entries.** "`n_pars` is unchanged" under the
`Tleaf` reporting entry, and "`n_pars` is unchanged at 19" under the
`shadow_cost` entry, were true when written and are not now: the count is 20 and
the one a caller wants is `n_theta`.

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
