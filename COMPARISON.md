# How phylloptim compares to other R leaf and canopy models

Three R packages occupy nearby ground. This document says precisely where they
overlap with phylloptim and where they don't, so that anyone deciding what to use —
or reviewing a paper about this one — can see the boundaries.

The short version: **all three assume stomatal conductance and none of them has a
hydraulic architecture.** That is the gap phylloptim fills. In exchange, all three
do things phylloptim does not: ship routines that fit models to measured data,
solve leaf temperature without linearising it, and work at ecosystem scale.

Packages surveyed at: `plantecophys` 1.4-6 (CRAN 2021-03-31), `bigleaf` 0.8.2
(CRAN 2022-08-22; repo 0.8.3 unreleased), `tealeaves` 1.0.7 (CRAN 2026-02-06;
source read at 1.0.6.1). All three are **pure R with no compiled code**.

---

## At a glance

| | **phylloptim** | **plantecophys** | **tealeaves** | **bigleaf** |
|---|---|---|---|---|
| Scale | single leaf | single leaf | single leaf | ecosystem / canopy |
| Direction | forward | forward **and** inverse | forward | almost entirely **inverse** |
| FvCB photosynthesis | yes | yes (C3 + C4) | no | inverted only, from GPP |
| Stomatal conductance | **emergent from hydraulics** | 4 empirical models + Tuzet | fixed input trait | fitted as regression |
| Leaf energy balance | linearised in `T_leaf`, forced convection, solved with the optimum; default off | full, iterated | **full nonlinear root-find** | no leaf; PM at surface |
| Xylem / root vulnerability | **yes** | no | no | no |
| Soil water potential | **yes, multi-layer**, or one potential | prescribed scalar (Tuzet only) | no | no |
| Root resistance | **yes, per layer** | explicitly not implemented | no | no |
| Profit / gain-risk optimisation | **yes** | Cowan-Farquhar with fixed λ | no | no |
| Thermal acclimation or damage | yes — an instantaneous thermal cost (Sicangco et al. 2026) and a damage layer, both default off | no | no | no |
| Fits to measured data | exact trait gradients for a likelihood you write (`vignette("fitting")`); no packaged fit | yes (A-Ci, Ball-Berry) | no | yes (that's the point) |
| Swappable *empirical* gs schemes | one, Medlyn's USO, as a separate route (`medlyn_model_gs`), outside the optimality solver | **yes — 4, plus Tuzet** | n/a | n/a |
| Swappable *hydraulic optimality* schemes | **nine, through one solver at identical drivers** (`vignette("the-models")`) | no (none are hydraulic) | no | no |
| Exact derivatives | **trait gradients through the implicit function theorem; exact rows for plant's TF24 runs (XAD)** | no | no | no |
| Language | C++ header-only | R | R | R |
| Cost per solve | a few µs (README, "Performance from R") | not measured; `mapply` over scalars | one `uniroot` per leaf | vectorised over time series |
| Licence | AGPL-3+ | GPL | MIT | GPL-2+ |

---

## `plantecophys` — the closest neighbour

Duursma's package is the one a leaf physiologist reaches for, and the one
phylloptim would be compared against. It is a steady-state leaf gas-exchange toolkit:
FvCB with all three limitations (Rubisco, RuBP, TPU), optional mesophyll
conductance, and — its distinguishing strength — **A-Ci curve fitting**.

**Where it is ahead of phylloptim:**

- **Inversion.** `fitaci` estimates Vcmax, Jmax and Rd from measured A-Ci curves
  by three methods: `nls` on the full hyperbolic-minimum curve (with a
  brute-force 25×25 grid refinement of starting values, which is what makes it
  robust), a bilinear `lm` fallback that always returns estimates, and De Kauwe's
  one-point method. `fitacis` batches it; `fitBB` fits Ball-Berry-family gs
  models. phylloptim has no packaged fit: it supplies exact gradients of its
  outputs in the traits, and `vignette("fitting")` shows how to write the fit.
- **Interchangeable stomatal models as first-class options.** `gsmodel=` selects
  Medlyn/USO (generalised with Duursma's `gk` exponent), Leuning, Ball-Berry, or
  a user-supplied multiplier. phylloptim has Medlyn's USO as one empirical route
  beside its optimality curves, and reports the Medlyn slope each optimum implies
  (`g1_eff`).
- **A working leaf energy balance.** `PhotosynEB` solves Tleaf with the full
  gs↔Tleaf feedback: an inner `uniroot` on the energy-balance residual for a
  given gs, wrapped in an outer `optimize` that closes the loop because gs itself
  depends on Tleaf through A. Free convection via a Grashof number in
  `abs(Tleaf-Tair)`, radiation conductance in parallel, Brutsaert atmospheric
  emissivity, series stomatal-plus-boundary-layer water path with a
  hypo/amphistomatous switch. phylloptim's energy balance linearises the leaf's
  emission about the air temperature, so the leaf temperature is affine in E and
  is solved with the stomatal optimum; it has forced convection only, the same
  Brutsaert emissivity, and defaults off.
- **C4** (`AciC4`, forward only) and a tidy library of VPD/RH/dewpoint
  conversions.

**Where phylloptim is ahead:**

- **Hydraulic architecture.** `plantecophys` has none. The single concession is
  `PhotosynTuzet`: a *prescribed* soil water potential `psis`, a **constant**
  leaf-specific conductance `kl` (default 2 mmol m⁻² s⁻¹ MPa⁻¹, not a function of
  ψ, so no vulnerability curve and no cavitation), and sigmoidal stomatal closure
  on ψ_leaf. Its own documentation states that "soil-to-root hydraulic
  conductance is not implemented". phylloptim has Weibull vulnerability curves for
  both xylem and roots, multi-layer soil with per-layer root resistance,
  gravitational head, and a root-collar potential that closes the whole
  soil→root→stem→leaf path.
- **Profit maximisation.** `FARAO`/`FARAO2` maximise `A − λE` with λ a fixed
  user-supplied marginal water cost — the Cowan-Farquhar exact solution that
  Medlyn et al. 2011 approximated. There is no hydraulic risk term. phylloptim's cost
  is derived from the vulnerability curve itself, so the optimum moves with the
  plant's actual hydraulic state.
- **Thermal acclimation and damage.** `plantecophys` has none: `EaV`, `EdVC`,
  `delsC` and the rest are fixed inputs. Its `new_T_responses` vignette is about
  updated *default values* from a literature review, not acclimation.
- **Speed and embeddability.** `plantecophys` solves the coupled Ci
  *analytically* by quadratic (`QUADP`/`QUADM`, ported from MAESTRA), which is
  elegant and fast per call — but `Photosyn`'s quadratic solve is dispatched
  through `mapply` because the discriminant test cannot vectorise, and
  `PhotosynEB`, `PhotosynTuzet`, `FARAO` are all `mapply`-over-scalars wrappers.
  `FindTleaf` is documented as not vectorised over `gs`. phylloptim is built to
  be called millions of times from inside a demographic model.

**Note on provenance:** the 2015 PLoS ONE paper predates `PhotosynTuzet`,
`FARAO2`, TPU fitting and the one-point method, so the paper understates the
package. Also: `Photosyn` does *not* apply the atmospheric-pressure correction to
photosynthesis rates (only `fitaci` does) — the docs warn about this in bold.

---

## `tealeaves` — the leaf energy balance done properly

Muir's package does exactly one thing forward, and does it better than anyone
else here: solve the leaf energy budget for equilibrium leaf temperature.
Thirteen exported functions, no photosynthesis, no stomatal model, no hydraulics
— all three deliberately out of scope, deferred to sister packages.

**What it gets right that phylloptim does not:**

- **A full nonlinear root-find**, not a linearisation. `stats::uniroot` on the
  residual of `R_abs − (S_r + H + L) = 0`, bracketed at `T_air ± 30 K`. Nothing
  is linearised, so both the `T_leaf⁴` radiative term and the exponential
  saturation-vapour-pressure term are handled exactly. Non-convergence returns
  `NA` rather than a silently wrong number.
- **Mixed free and forced convection, with surface asymmetry.**
  `Nu = ((a·Re^b)^3.5 + (c·Gr^d)^3.5)^(1/3.5)`, with laminar/turbulent forced
  branches switching at Re = 4000, and free-convection coefficients that differ
  depending on whether a given surface is buoyantly stable or unstable (upper
  surface with T_leaf > T_air, or lower with T_leaf < T_air). Grashof uses a
  *virtual*-temperature difference, assuming the leaf interior is at 100 % RH.
  `Ar()` exposes Gr/Re² so users can check which regime dominates.
- **Amphistomy.** Conductances are computed per surface — stomatal in parallel
  with cuticular, in series with that surface's boundary layer, then the two
  surfaces summed — with the split set by `logit_sr`. This is the package's
  distinctive contribution and nothing else here has it.
- **Units.** `units` is in `Depends`; parameter constructors validate and convert
  on entry. phylloptim documents units and refuses out-of-range inputs by name,
  but does not carry units on its values.

**Where phylloptim differs:** `g_sw` in `tealeaves` is a *fixed input trait*. There is no water supply side at all, so the package cannot tell you when a leaf must close. phylloptim derives gs from hydraulic supply, which is the opposite starting point. It also carries a cuticular conductance in parallel with the stomata (the `g_min` trait, temperature-dependent, drawn through the same hydraulic path), but on one surface: it has no amphistomy.

**The tension worth naming.** `tealeaves` shows what a rigorous leaf energy
balance costs: one `uniroot` per leaf, per condition, scalar internals,
parallelised by `furrr` rather than vectorised. phylloptim evaluates its
objective many times inside every solve, millions of times per model run, and a
nonlinear balance nested there would be a root-find inside each evaluation. The
linearisation keeps the leaf temperature explicit in E; free convection, which
would make it implicit again, is left out for that reason.

`tealeaves` also models **equilibrium, not transient** temperature: no thermal
capacitance, no time-stepping, horizontal leaves only, one homogeneous leaf
temperature.

---

## `bigleaf` — a different problem entirely

Knauer's package overlaps phylloptim almost nowhere, and is included here mainly so
that the comparison is not mistaken for one. It works at **ecosystem scale**,
treating the land surface as one big leaf, and it runs **backwards**: you feed it
half-hourly eddy-covariance fluxes plus met drivers, and it recovers the physical
and physiological properties those fluxes imply. Sixty-one exported functions, of
which sixteen are unit conversions.

What it derives: surface conductance (inverted Penman-Monteith, or flux-gradient),
aerodynamic conductance (log-wind profile plus Monin-Obukhov stability
correction, with Thom / Choudhury / Su boundary-layer options), the
Jarvis-McNaughton decoupling coefficient Ω, potential ET (Priestley-Taylor or PM),
surface temperature and surface CO₂, energy-balance closure, WUE metrics
(including Beer's IWUE and Zhou's uWUE), and Oren-type stomatal sensitivity.

Relevant negatives, all verified against the source:

- **No leaf energy balance.** `surface.conditions()` gets aerodynamic surface
  temperature by inverting the bulk transfer equation, `Tsurf = Tair + H/(ρ·cp·Ga)`
  — H is an *input*, not something solved for. No leaf, no Nusselt/Grashof
  treatment, no nonlinear Tleaf solve.
- **FvCB only backwards.** `photosynthetic.capacity()` inverts the FvCB rate
  expressions to recover bulk canopy Vcmax and Jmax from GPP, Ci and PPFD, with
  the limitation state *assumed* from PPFD windows. There is no forward
  `A(Ci, T, Q)` and no A-gs coupling anywhere. The docs warn the function "should
  be used with care" and that bulk canopy parameters are per ground area and not
  comparable to leaf-level values.
- **No hydraulics, no soil water.** Every `psi_` in the source is a
  Monin-Obukhov stability function. The only water-related data handling is
  dropping hours after rain.
- **No optimisation.** `stomatal.slope()` fits Medlyn's g1 as a regression
  coefficient — a parameter derived from optimality theory, used here without it.
- **No thermal acclimation.** `Arrhenius.temp.response()` *removes* an assumed
  temperature dependence to normalise to 25 °C; it does not model acclimation of
  the parameters. `optimum.temperature()` is an empirical boundary-line Topt of
  GPP.
- **No footprint model.** (Worth stating because it is often assumed: there are
  zero occurrences of `footprint`, `Kljun` or `Kormann` in the source.) Also no
  flux partitioning — GPP and Reco must be supplied — no gap-filling, and no
  sun/shade or multi-layer canopy, since the big-leaf assumption is the premise.

The one genuine point of contact is conceptual: `bigleaf` measures the
*consequences* of stomatal regulation at ecosystem scale, which is the scale
`plant` operates at when it embeds phylloptim. It is a validation target, not an
alternative.

---

## Where that leaves phylloptim

**Every leaf gas-exchange package in R takes stomatal conductance as given — as a fitted empirical model, a regression coefficient, or a fixed trait — and phylloptim derives it from hydraulic supply and carbon cost instead.** The hydraulics behind it (two vulnerability curves, multi-layer soil, root resistance, gravitational head, gain-risk optimisation) have no counterpart in the other three, and it is the only one of the four written to be embedded in a larger model rather than driven from a console.

Two things follow.

**It compares hydraulic optimality theories through one numerical core.** `plantecophys` already compares *empirical* stomatal schemes well (`gsmodel=`, `PhotosynTuzet`, `FARAO`). What none of the three can compare is hydraulically explicit formulations, which is where the live argument sits: gain-risk against least-cost against carbon-maximisation variants. phylloptim runs nine of them through one solver at identical drivers. A companion manuscript (`Falster-stomatal_analytical_analysis`) shows that these models all satisfy `dA/dE = λ` and differ only in λ(state), the marginal cost of water, so a shared core makes the comparison like for like by construction (`vignette("the-models")`).

**It is fast and differentiable.** None of the three has automatic differentiation; all of them difference or grid-search when they need gradients (`fitaci`'s 25×25 start grid exists because `nls` on this objective is fragile). phylloptim returns trait gradients of the operating point through the implicit function theorem at the converged point (`leaf_gradient()`, `leaf_gradient_batch()`; the curvature terms are differences of an analytic first derivative, and a pinned point falls back to a difference of the solve), and exact derivative rows for the TF24 curves, through which plant differentiates whole runs.

The gaps run the other way just as clearly:

1. **No packaged inversion.** `plantecophys::fitaci` is the most-used function in this space. phylloptim supplies the gradients a fit needs, and `vignette("fitting")` writes one, but there is no `fitaci` equivalent.
2. **A simpler energy balance than either leaf-scale neighbour**: linearised in the leaf temperature, forced convection only, one surface, and off by default.
3. **C3 only**, and no mesophyll conductance.
