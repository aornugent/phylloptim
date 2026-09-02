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
