# How `leaf` compares to other R leaf and canopy models

Three R packages occupy nearby ground. This document says precisely where they
overlap with `leaf` and where they don't, so that anyone deciding what to use —
or reviewing a paper about this one — can see the boundaries.

The short version: **all three assume stomatal conductance and none of them has a
hydraulic architecture.** That is the gap `leaf` fills. In exchange, all three do
things `leaf` does not: fit models to measured data, solve leaf temperature
properly, and work at ecosystem scale.

Packages surveyed at: `plantecophys` 1.4-6 (CRAN 2021-03-31), `bigleaf` 0.8.2
(CRAN 2022-08-22; repo 0.8.3 unreleased), `tealeaves` 1.0.7 (CRAN 2026-02-06;
source read at 1.0.6.1). All three are **pure R with no compiled code**.

---

## At a glance

| | **leaf** | **plantecophys** | **tealeaves** | **bigleaf** |
|---|---|---|---|---|
| Scale | single leaf | single leaf | single leaf | ecosystem / canopy |
| Direction | forward | forward **and** inverse | forward | almost entirely **inverse** |
| FvCB photosynthesis | yes | yes (C3 + C4) | no | inverted only, from GPP |
| Stomatal conductance | **emergent from hydraulics** | 4 empirical models + Tuzet | fixed input trait | fitted as regression |
| Leaf energy balance | minimal, default off | full, iterated | **full nonlinear root-find** | no leaf; PM at surface |
| Xylem / root vulnerability | **yes** | no | no | no |
| Soil water potential | **yes, multi-layer** (separable — PLAN 7b) | prescribed scalar (Tuzet only) | no | no |
| Root resistance | **yes, per layer** | explicitly not implemented | no | no |
| Profit / gain-risk optimisation | **yes** | Cowan-Farquhar with fixed λ | no | no |
| Thermal acclimation or damage | yes — an instantaneous PSII cost here (Sicangco et al. 2026, default off), and a lasting damage ratchet in plant's TF24t | no | no | no |
| Fits to measured data | **no** | yes (A-Ci, Ball-Berry) | no | yes (that's the point) |
| Swappable *empirical* gs schemes | Medlyn present, not dispatched | **yes — 4, plus Tuzet** | n/a | n/a |
| Swappable *hydraulic optimality* schemes | **TF24 and Sperry ProfitMax runnable at identical drivers; Prentice14 planned** | no (none are hydraulic) | no | no |
| Exact derivatives | **yes — forward-mode AD (XAD)** | no | no | no |
| Language | C++ header-only | R | R | R |
| Cost per solve | ~4 µs | not measured; `mapply` over scalars | one `uniroot` per leaf | vectorised over time series |
| Licence | AGPL-3+ | GPL | MIT | GPL-2+ |

---

## `plantecophys` — the closest neighbour

Duursma's package is the one a leaf physiologist reaches for, and the one `leaf`
would be compared against. It is a steady-state leaf gas-exchange toolkit:
FvCB with all three limitations (Rubisco, RuBP, TPU), optional mesophyll
conductance, and — its distinguishing strength — **A-Ci curve fitting**.

**Where it is ahead of `leaf`:**

- **Inversion.** `fitaci` estimates Vcmax, Jmax and Rd from measured A-Ci curves
  by three methods: `nls` on the full hyperbolic-minimum curve (with a
  brute-force 25×25 grid refinement of starting values, which is what makes it
  robust), a bilinear `lm` fallback that always returns estimates, and De Kauwe's
  one-point method. `fitacis` batches it; `fitBB` fits Ball-Berry-family gs
  models. `leaf` cannot invert anything.
- **Interchangeable stomatal models as first-class options.** `gsmodel=` selects
  Medlyn/USO (generalised with Duursma's `gk` exponent), Leuning, Ball-Berry, or
  a user-supplied multiplier. `leaf` has a Medlyn implementation but it is a
  standalone R-callable alternative, not on the solve path.
- **A working leaf energy balance.** `PhotosynEB` solves Tleaf with the full
  gs↔Tleaf feedback: an inner `uniroot` on the energy-balance residual for a
  given gs, wrapped in an outer `optimize` that closes the loop because gs itself
  depends on Tleaf through A. Free convection via a Grashof number in
  `abs(Tleaf-Tair)`, radiation conductance in parallel, Brutsaert atmospheric
  emissivity, series stomatal-plus-boundary-layer water path with a
  hypo/amphistomatous switch. `leaf`'s energy balance is a single explicit pass
  with forced convection only and a fixed longwave offset, and defaults off.
- **C4** (`AciC4`, forward only) and a tidy library of VPD/RH/dewpoint
  conversions.

**Where `leaf` is ahead:**

- **Hydraulic architecture.** `plantecophys` has none. The single concession is
  `PhotosynTuzet`: a *prescribed* soil water potential `psis`, a **constant**
  leaf-specific conductance `kl` (default 2 mmol m⁻² s⁻¹ MPa⁻¹, not a function of
  ψ, so no vulnerability curve and no cavitation), and sigmoidal stomatal closure
  on ψ_leaf. Its own documentation states that "soil-to-root hydraulic
  conductance is not implemented". `leaf` has Weibull vulnerability curves for
  both xylem and roots, multi-layer soil with per-layer root resistance,
  gravitational head, and a root-collar potential that closes the whole
  soil→root→stem→leaf path.
- **Profit maximisation.** `FARAO`/`FARAO2` maximise `A − λE` with λ a fixed
  user-supplied marginal water cost — the Cowan-Farquhar exact solution that
  Medlyn et al. 2011 approximated. There is no hydraulic risk term. `leaf`'s cost
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
  `FindTleaf` is documented as not vectorised over `gs`. `leaf` is built to be
  called millions of times from inside a demographic model.

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

**What it gets right that `leaf` does not:**

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
  on entry. `leaf` documents units in comments and trusts the caller.

**Where `leaf` differs:** `g_sw` in `tealeaves` is a *fixed input trait*. There is
no water supply side at all, so the package cannot tell you when a leaf must
close. `leaf` derives gs from hydraulic supply, which is the opposite starting
point.

**The tension worth naming.** `tealeaves` shows what a rigorous leaf energy
balance costs: one `uniroot` per leaf, per condition, scalar internals,
parallelised by `furrr` rather than vectorised. `leaf` runs ~10³ inner
evaluations inside a golden-section search per solve, millions of times per
model run. Nesting a `tealeaves`-grade energy balance inside that is the trade
discussed in PLAN.md item 13 — and the reason the free-convection term, which is
what makes the balance implicit, is the piece explicitly *not* recommended.

`tealeaves` also models **equilibrium, not transient** temperature: no thermal
capacitance, no time-stepping, horizontal leaves only, one homogeneous leaf
temperature.

---

## `bigleaf` — a different problem entirely

Knauer's package overlaps `leaf` almost nowhere, and is included here mainly so
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
`plant` operates at when it embeds `leaf`. It is a validation target, not an
alternative.

---

## Where that leaves `leaf`

The honest positioning, in one sentence: **every leaf gas-exchange package in R
takes stomatal conductance as given — as a fitted empirical model, a regression
coefficient, or a fixed trait — and `leaf` derives it from hydraulic supply and
carbon cost instead.**

That is a real gap, and the hydraulics behind it (dual vulnerability curves,
multi-layer soil, root resistance, gravitational head, gain-risk optimisation) has
no counterpart in any of the three. So does the C++ implementation: this is the
only one of the four designed to be embedded in a larger model rather than driven
from a console.

Two things follow from that, and they are the reasons this is worth being a
package rather than a file inside plant.

**It can become a platform for comparing hydraulic optimality theories, which
none of these can.** Note the distinction carefully, because it is easy to
overclaim. `plantecophys` is already a good comparison platform for *empirical*
stomatal schemes — `gsmodel=` switches between Ball-Berry, Leuning and Medlyn/USO,
`PhotosynTuzet` adds a water-potential feedback, and `FARAO` adds Cowan-Farquhar
optimal. What it cannot do is compare **hydraulically explicit** formulations,
because none of its schemes has a vulnerability curve. And that is exactly where
the live theoretical argument sits: Sperry-style gain-risk against Prentice
least-cost against the various carbon-maximisation and Wolf-Anderegg-Pacala
variants. This package already contains the TF24 gain-risk formulation and, as
inherited second-class code, both Sperry (2017) and Medlyn (2011).

And there is a sharper version of the claim available, from a companion manuscript
(`Falster-stomatal_analytical_analysis`): those models all maximise a profit, so
they all satisfy `dA/dE = λ` and differ **only** in λ(state), the marginal cost of
water. Given λ, each one collapses to the Medlyn USO functional form — so USO is
the generic solution of the family rather than a competing model, which is why
mutually incompatible schemes fit gas-exchange data equally well. If λ is what
differs, then λ is what should be pluggable, and a package that runs six λ
functions through one shared numerical core makes the comparison apples-to-apples
*by construction*. That is PLAN.md item 7a, and it is a different and more
interesting contribution than a fourth implementation of one scheme.

**It is fast and differentiable, which makes calibration tractable.** None of the
three has automatic differentiation; all three finite-difference or grid-search
when they need gradients (`fitaci`'s 25×25 brute-force start grid is the
tell — it exists because `nls` on this objective is fragile). Forward-mode AD at
4 µs a solve is a different regime. The caveat, stated in full in PLAN.md item 12:
AD currently differentiates with respect to the collar potential, not with respect
to traits, so this is an argument about an architecture until the templated
`Leaf<T>` of item 11 lands and a calibration vignette demonstrates it. Worth noting
that plant already has direct evidence for the underlying claim — the analytic
`dprofit_droot_collar_psi` exists precisely because finite-differencing this
objective was too noisy to drive acclimation tracking.

The gaps run the other way just as clearly, and they are the roadmap:

1. **No inversion.** `plantecophys::fitaci` is the single most-used function in
   this space and `leaf` has no equivalent. Fitting *hydraulic* traits from
   measured A/gs/ψ data would be a new capability rather than a reimplementation
   — see PLAN.md item 12.
2. **A weaker energy balance than either leaf-scale competitor**, and one that is
   off by default. PLAN.md item 13 sets out which parts are worth fixing (the
   leaf-to-air VPD, immediately) and which are not (free convection, on speed
   grounds).
3. **No R interface at all yet.** PLAN.md item 6.
