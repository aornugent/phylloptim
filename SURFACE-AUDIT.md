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
