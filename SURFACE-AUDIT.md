# Every type the AD surface coins, tested against what already exists

The same question asked of each: is this a spelling of something upstream, odelia
or this package already has? Four of eight were. The test is worth applying to
anything added later.

## Dissolved -- a spelling of something that exists

| coined | what it already is |
|---|---|
| `SolvedPoint` / `place_solved_point` | upstream's `evaluate_root_collar_psi(target)`, which *"leaves exactly the same outputs as find_root_collar_psi"* and re-derives sigma, ci and every output **from the collar alone**. Four of five fields were redundant. What survives is one double and one enum |
| `MarginalParts{marginal, V, dci_dpsi}` | `V` and `dci_dpsi` are `sigma.slope` and `ci.slope` -- the eighteenth spelling of the value/slope pair `two-paths.md` swept seventeen of. `marginal_assembled` returns a scalar |
| `WhichBound::DryRootPsiCrit` | a bare `return in.root_psi_crit;`, which existed only because root_psi_crit was a free trait. Under (P50, c) it is an ordinary expression AD chains itself. The arm dissolves and `DryBoundArm` folds into the two that remain |
| `with_slope` (in plant) | belongs in odelia: its `for_each_active` exists for `visit_active`, and a model carrying that hook is carrying odelia's problem. Moved |

## Kept, but reshaped

**`CollarPoint`** earns its keep -- its five coordinates *"must be seeded
consistently or the assembly below is a derivative of nothing"*, and a first
attempt that moved only the collar **disagreed by 40 to 100 per cent**. But it
REPEATS `CollarCoords`' two fields instead of deriving from them, which is the
shape odelia already fixed: *"a recording row IS a program row plus its state, so
it derives rather than repeating the two fields."* And its own comments admit a
duplication -- `dEup_dp` IS `transpiration`'s slope. As pairs it is four members,
not five, and consistent seeding stops being a warning and becomes a constraint.

**`SupplyDraw`** earns its keep for the same reason and one more: it holds the
collar it was taken at, so *"the three values cannot come from different ones"*.
Its `flux`/`slope` become a `with_slope`. Its `uptake`/`duptake_dp` do NOT -- the
draws are active and their slopes are SUPPLIED at double, which is
`implicit_value`'s shape, not the pair's. Two different pairings, correctly
distinguished, and worth keeping distinguished.

**`FixedCollarEval{profit, uptake, feasible}`** is `LeafOutputs<double>` plus one
member -- derive, do not repeat. The `feasible` flag itself is legitimate: an
infeasible collar is an EXPECTED outcome on a gradient sweep (1536 of the golden
file's 5184 rows are refused), so a throw would be the wrong shape. Upstream
agrees and has the same convention under a different name -- `_checked` variants
returning value-plus-flag, flattened for the R boundary, exactly as
`profit_at_fixed_collar_values` does. Two spellings of one convention; take
upstream's naming.

## Kept as-is

`LeafOutputs` is the leaf-to-plant contract and is genuinely this package's.
`clamp_count` and `operating_point_kind_count` stay -- the latter counts from the
last enumerator, which is upstream's own `n_cost_curves` pattern arrived at
independently, and convergence is evidence rather than duplication.

## The rule this leaves

Before coining a type on the differentiable surface, ask what it is a pair of, a
record of, or a selection between -- and whether odelia, upstream or this package
already names that. Four of eight did. The two that survived intact both did so
for the same reason: they make a partially-built state unrepresentable, after a
measured bug from building one.

---

# The rule, sharpened by everything since

The test that has actually separated the survivors from the dissolved is not
"does this group things that belong together". It is:

> **Does it make a wrong state unrepresentable, or a wrong pairing checkable?**
> A return type with an invariant earns its keep. A parameter bundle without one
> does not.

Nine for nine now, in both directions.

## Survives -- carries an invariant

| | what it prevents |
|---|---|
| `SupplyDraw` | holds the collar it was taken at; `check_draw` compares it at three call sites, so a flux from one collar and a slope from another cannot be read together |
| `CollarPoint : CollarCoords` | five coordinates that must be seeded together, after a first attempt that moved only the collar and disagreed with a difference by 40 to 100 per cent |
| `SupplyAt` | a view with two genuine owners -- `held_supply()` builds one from this object's members where no owned bundle exists |

## Dissolved -- grouped arguments, or spelled something that existed

| | what it already was |
|---|---|
| `SolvedPoint` | upstream's `evaluate_root_collar_psi`, which re-derives sigma, ci and every output from the collar alone |
| `MarginalParts` | its two extra fields are `sigma.slope` and `ci.slope` -- the eighteenth spelling of the pair `two-paths.md` swept seventeen of |
| `WhichBound::DryRootPsiCrit` | a bare identity, live only while `root_psi_crit` was a free trait. Under `(P50, c)` it is an ordinary expression AD chains itself |
| `with_slope` in plant | odelia's, because `for_each_active` is odelia's obligation |
| `ProfitInputs` / `LeafInputs` | the parameter enumeration upstream already maintains, given a scalar |
| `SupplyValues` | storage only -- every use either declares it or immediately calls `.at()`. plant already owns the vectors |
| `CurveTrait` + two functions | a SELECTOR where a RECORD was wanted. Upstream's integral bundle shows the shape; the integrand now has its sibling |
| `sigma_star` / `ci_star` | threaded members. Every call passes `l.opt_psi_stem_, l.ci_`; the only other is `profit_at` forwarding to `collar_coords_at` |

## What the lens produced rather than deleted

Two pieces that read like they should always have been upstream's:

* **`VulnerabilityDerivatives`** -- the integrand's bundle, in the same header and
  shape as the integral's. Verified against central differences (worst 1.67e-08)
  and against the free cross-check the other bundle gives: **f IS dG/dpsi**, and
  the two agree to 0.00e+00 at every tested potential.
* **`graft.hpp`** -- one graft for both curves. Written twice, one copy keeps a
  partial at fixed `b` and the other carries the chain in `c`; that is a finite,
  plausible, wrong row that nothing reports.

Three refactors in a row -- extracting `rows_in_P50`, then `graft_integral`, then
routing `stem_integral_at` through it -- left the verified numbers unchanged to
the digit: 2.12e-07 and 1.84e-07. Same digits is the check that an abstraction is
the thing that was already there, rather than a rewrite wearing its name.

---

# The lens, applied again on the way through the surface

Same question, same two outcomes. Thirteen for thirteen now.

## Dissolved -- a spelling of something that exists

| coined | what it already is |
|---|---|
| `duptake_dpsi_at` | an OVERLOAD of upstream's `duptake_dpsi`, on the argument that already distinguishes them. Upstream's name plus a `SupplyAt<T>`, exactly as `uptake_at` took one -- so the second name buys nothing |
| `E_from_soil_at` | two statements: size a vector, call `uptake_at`. `supply_draw_at` was its only caller and does both |
| `assim_colimited_slope_kernel` | a TANGENT THROUGH `assim_colimited_kernel`, which is how upstream takes A' at double. A hand-written twin is a second definition of one function, free to disagree while every value stays finite |
| `electron_transport_kernel` | a `leaf_pars<T>` overload of upstream's own `electron_transport()`, whose arguments were already explicit |
| `MarginalParts` | its two extra fields ARE `sigma.slope` and `ci.slope`, and now live there. `marginal_assembled` returns a scalar |

## Kept -- carries an invariant

| | what it prevents |
|---|---|
| `SupplyDraw` + `supply_draw_at` | the only way to make one, so the flux, its collar slope and the per-layer draws cannot come from different collars; `check_draw` refuses a stale one at every read |
| `CollarCoords` with both halves | a consumer seeded from separate places can pair a value from one point with a slope from another. A first attempt that moved only the collar disagreed with a difference by 40 to 100 per cent |
| `transpiration_at` | upstream's `transpiration()` given its parameters, and the name that keeps the STEM flux distinct from the soil draw -- two quantities that agree only to the splines' round-trip, where reading the wrong one was 8 per cent in `dci/dcollar` |

## What the lens produced rather than deleted, again

**Upstream's kernels take their parameters.** The member-reading forms stay and
forward to them, so there is no parallel copy -- and the pack form is what
carries a trait's row. The alternative was a second family of kernels beside
upstream's, which is the shape this rule exists to refuse.

The temperature responses are templated by pure substitution rather than
factored into `reference_value x response`. They ARE exactly linear in the
reference value, so the factoring is algebraically free -- and it reassociates,
which costs the bit-identity that makes the double path checkable.

## The measurement that changed how these are checked -- and the wrong turn in it

`dsigma/dcollar` is 1 + 9e-6. Comparing it against a difference agrees to 7
digits while the informative part -- the 9e-6 -- was 8 per cent out, so the first
lesson stands: **compare the quantity that carries the information, not the one
that is easy to print.**

⚠️ **THE SECOND CONCLUSION DRAWN HERE WAS WRONG AND IS CORRECTED.** It read: a
difference of the model cannot referee these slopes at the shipped resolution,
because sigma comes from an inverse spline whose round-trip dominates them, and
the closed form is right. The evidence was a resolution sweep in which the
difference converged onto the closed form.

That sweep was real and the inference was not. Table and curve coincide as the
resolution rises, so convergence is consistent with EITHER being right; it
cannot separate them. What separates them is that **the solve ran on the table**,
so the operating point is the table's argmax and the closed form is a different
function that is not stationary there. Reading graft.hpp as "take the curve" got
its rule backwards -- the rule is that the VALUE is the table's, precisely
because the solve ran on it, and only the ROWS come from the curve.

The cost of the mistake was visible and specific: `M` read 6.4e-02 at a collar
where the solve had driven upstream's own `dprofit` to 1e-12. The envelope
omission -- dropping `M * dp/dtheta` at an interior point BECAUSE M is zero
there -- was not true of that number.

Fixed by using what upstream already had, in three places:
`stem_curve_integral_inverse_deriv` for sigma's response,
`stem_curve_integral_deriv` through `graft_curve` for every conductivity, and
the SOIL DRAW as the point the inverse is read at. `marginal_at` is then
bit-identical to `dprofit_at_collar_psi`, and a difference of the model referees
the slopes at the shipped resolution: 9.4e-12 and 1.5e-06.

**The lens question caught all three.** Each was a quantity this surface had
spelled for itself -- `1/f(sigma)`, the closed-form conductivity, the stem
flux -- where upstream had the one the model actually uses. Asking "what does
upstream already have for this" is not only a way to delete code; here it was
the difference between a correct gradient and a plausible one.

---

# The lens at the consumer boundary

The same question asked of everything plant reached for. Nothing new was coined
that upstream or odelia could name.

## Dissolved

| coined | what it already is |
|---|---|
| `Leaf::SolvedPoint` (five fields) | a collar and an arm. Four fields were re-derivable from the collar alone, which is exactly what `evaluate_root_collar_psi` does |
| `place_solved_point` / `solved_point` | `replay_operating_point`, below, plus reading `opt_root_psi_` and the kind |
| plant's `with_slope` | odelia's, because `for_each_active` is odelia's obligation. Eleven users, one alias |
| `root_resistances_from_carbon` | an OVERLOAD of upstream's `root_network_from_carbon`, on the output arguments -- the same distinction `duptake_dpsi` and `uptake_at` already make |
| plant's curvature probe | it differenced the marginal to referee an ANALYTIC curvature. The curvature IS that difference now, so the probe compares a difference with itself |
| four of the five R instruments | already answered: the two curves' trait rows by `test_leaf`'s 11292-point check, `set_traits`' stale derived state by `test_set_traits_matches_a_fresh_leaf`, the fixed-collar clamp by upstream's own `_checked` readers |

## Kept -- carries an invariant

| | what it prevents |
|---|---|
| `replay_operating_point` | putting a collar back WITHOUT its arm. `evaluate_root_collar_psi` restores every number bit-for-bit and then tags the point `prescribed`, which is correct from its side and leaves `collar_at` with no condition to place a replayed interior point at. One call, so the pair cannot be split -- and this exact mistake had already cost a throw once, inside `marginal_collar_slope` |
| `dry_bound_is_root_limit_` | which limit closed the dry end. Upstream keeps ONE `BoundaryCrit` kind for two conditions, and re-deriving the comparison would repeat a root-find and differentiate a selector |
| `operating_point_kind_count`, `ncontrol_default` | a consumer sizes an array and seats a spline from them. A literal on that side drifts silently: a kind added here lands outside a tally that stops at a stale count, and a curve built on a different knot number is a different curve with every number plausible |

## What the boundary turned up that no lens would have

Four things in plant that had never compiled, because the phylloptim API changed
in the same commit that broke them: a rename that missed one site, a lost
`template` line, a state read into a `double` inside a templated environment, and
a call carrying an argument this branch's design had removed. A merge that
cannot be built hides all four equally, and the syntax-only sweep over all 27
translation units is what surfaced them.

---

# What only the consumer could find

Three defects survived every check in this package and were caught by plant. Each
is worth knowing as a CLASS, because the same blindness is easy to rebuild.

**A vector sized by the wrong count.** `duptake_dpsi`'s per-layer rows were as
long as the deepest ROOTED layer where the uptake they graft onto is as long as
the soil profile. Every fixture here puts carbon in every layer, so the two
counts coincided and 63 differenced rows agreed perfectly. plant runs
shallow-rooted plants. The test drives each profile twice now -- roots through it,
and roots reaching only the top layer -- and asserts the LENGTH, not only the
entries.

**A pack filled with zeros.** Thirteen of nineteen slots were set and the rest
left at zero rather than the model's values -- a wrong VALUE, not a missing row.
Nothing here reads those slots, because they belong to cost curves this surface
does not cover. plant's ladder asks that every trait either has a column or is
refused BY NAME, and named one.

**A refusal that costs more than it says.** Returning NaN from the curvature at a
collar within a step of a bound reads locally like a careful admission. In plant
it refuses the whole SWEEP: `ranges` comes back zero and the gradient is never
produced, so one awkward point costs a stand's gradient. A one-sided difference
is the answer there.

The common shape: **this package's fixtures are regular where its consumer's
states are not.** Every layer rooted, every collar interior, every parameter
read. The cheapest guard against it is to drive one fixture that is irregular in
each of those ways, which is what the two new sweeps now do.
