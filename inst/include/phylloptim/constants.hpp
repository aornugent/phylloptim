// -*-c++-*-
#ifndef PHYLLOPTIM_CONSTANTS_HPP_
#define PHYLLOPTIM_CONSTANTS_HPP_

// Physical and physiological constants for the leaf model.
//
// Split out of leaf_model.h and changed from `static const double` (internal
// linkage, so one copy per translation unit) to `inline constexpr double`,
// which is the right spelling for a header-only library: one entity, usable in
// constant expressions, no per-TU duplication.
//
// Two names changed on the way out of plant, both because they were namespace
// -scope one-letter identifiers in a public header:
//   * `R` (the universal gas constant) -> `gas_constant`. `R` at plant namespace
//     scope was a live collision hazard in a project where `R` also names the
//     language, `R_` prefixes the C API, and several strategies declare local
//     `R`s of their own.
//   * `n` ("number of integration steps") was dead -- nothing read it -- and is
//     dropped. It shadowed local `n`s in every consumer that had one.

#include <cmath>

namespace phylloptim {

// double check best namespace for constants (private vs global)
// kJ mol ^-1
inline constexpr double vcmax_ha = 60000;
// kJ mol ^-1
inline constexpr double vcmax_H_d = 200000;
// kJ mol ^-1
inline constexpr double vcmax_d_S = 650;

// kJ mol ^-1
inline constexpr double jmax_ha = 30000;
// kJ mol ^-1
inline constexpr double jmax_H_d = 200000;
// kJ mol ^-1
inline constexpr double jmax_d_S = 650;

// umol ^ -1 mol ^ 1
inline constexpr double gamma_25 = 42.75;
// kJ mol ^-1
inline constexpr double gamma_ha = 37.83e3;

// umol mol ^-1
inline constexpr double kc_25 = 404.9 ;
// kJ mol ^-1
inline constexpr double kc_ha = 79.43e3;

// umol mol ^-1
inline constexpr double ko_25 = 278400 ;
// kJ mol ^-1
inline constexpr double ko_ha = 36.38e3;

// umol mol^-1 -> Pa is deliberately NOT a constant here: it is 1e-6 * P. See
// Leaf::umol_per_mol_to_Pa_, derived from atm_kpa_. The old `0.1013` was 101.3
// kPa hard-coded; don't reinstate it.

// The molar mass of water -- 18.015 g/mol, from the standard atomic weights -- and
// the ONE place it is written down. Both conversion directions derive from it, so
// they are reciprocal by construction; do not reintroduce a second literal (#51).
//
// Both names below are kept because plant `using`-declares both.
inline constexpr double molar_mass_h2o = 0.018015; // kg mol^-1
// kg mol^-1, for converting a molar water flux back to kg.
inline constexpr double kg_per_mol_h2o = molar_mass_h2o;
// mol H2O kg^-1 -- the reciprocal, now genuinely so.
inline constexpr double kg_to_mol_h2o = 1.0 / molar_mass_h2o;
// mol mol ^-1 / (umol mol ^-1)
inline constexpr double umol_to_mol = 1e-6;
// Pa kPa^-1
inline constexpr double kPa_to_Pa = 1000.0;

// universal gas constant J mol^-1 K^-1. Exact since the 2019 SI redefinition
// (R = N_A k_B), so it is not a value to round.
//
// ⚠️ Read only by arrh_curve and peak_arrh_curve, always as `Ea/(R*T)` with Ea/RT
// of order 24 -- so a relative change here reaches the rates ~24x amplified, and a
// change too small to look like it matters does.
inline constexpr double gas_constant = 8.314462618153240;

//convert deg C to deg K
inline constexpr double C_to_K = 273.15;

//H20:CO2 stomatal diffusion ratio
inline constexpr double H2O_CO2_stom_diff_ratio = 1.67;

// Boundary-layer conductance to water vapour over that to heat, for a
// hypostomatous leaf: the factor plantecophys applies to its two-sided heat
// conductance (`Gbw = 1.075 * Gbh` in LeafEnergyBalance, at StomatalRatio 1).
// Read only by the reported boundary-layer and stomatal conductances (#122);
// nothing in the solve reads it.
inline constexpr double bl_vapour_heat_ratio = 1.075;

// rho_w * g, as a head in MPa per metre of depth.
//
// ⚠️ 9.8e-3 is 6.8e-04 below `1000 * 9.80665 / 1e6`, above this package's 1e-04
// "real difference" threshold, and it reaches results through the per-layer head.
// Deliberately not "corrected": the right value needs a water-density convention
// stated first (nominal 1000, 4 C and 25 C span 0.3%, four times the discrepancy),
// so it is a modelling decision rather than arithmetic.
inline constexpr double gravity_head = 9.8e-3; // MPa / m

// --- Penman-Monteith leaf energy balance (minimal core; #523) -----------------
// See notes/penman-monteith/. These back Leaf::leaf_temp_from_E and the es/Delta
// helpers. Only used on the (default-off) use_energy_balance_ path.
// latent heat of vaporisation of water, J kg^-1, at about 21.5 C.
//
// ⚠️ NOT a 25 C value, which an earlier comment here claimed: lambda(T) ~= 2.501e6
// - 2370*T gives 2.442e6 at 25 C. Left as-is because on the (default-off)
// Penman-Monteith path lambda should become a function of the leaf temperature the
// balance solves for, not a better constant -- see #28.
inline constexpr double latent_heat_vap = 2.45e6;
// volumetric heat capacity of air, J m^-3 K^-1
inline constexpr double vol_heat_cap_air = 1200.0;
// PAR energy conversion for SUNLIGHT: ~4.57 umol photons per J of PAR. The ratio
// is a property of the source spectrum (a red/blue LED is ~5-5.5), so it is only
// correct for a leaf in daylight. It is read on the energy-balance path alone, to
// turn the intercepted photon flux back into the PAR band's energy; photosynthesis
// reads the photon flux directly and never touches it.
inline constexpr double umol_par_per_joule = 4.57;
// --- Longwave radiation (#97, #28) -------------------------------------------
//
// The net longwave is an ISOTHERMAL net radiation (as if the leaf sat at air
// temperature) plus a RADIATIVE CONDUCTANCE for the leaf's own departure from it.
// The absorbed shortwave is each intercepted band times the leaf's own
// absorptance (alpha_PAR_, alpha_NIR_ in leaf_model.hpp):
//
//   Rn_iso = absorbed shortwave - (1 - ema)*sigma*Tair_K^4
//   g_rad  = 4*eps_leaf*sigma*Tair_K^3
//   Tleaf  = Tair + (Rn_iso - lambda*E) / (vol_heat_cap_air/ra + g_rad)
//
// ⚠️ THE LINEARISATION IS WHAT KEEPS THE BALANCE EXPLICIT. A leaf emits
// eps*sigma*Tleaf_K^4, which in Rn would make Rn a function of the temperature
// being solved for: an inner iteration, and no constant dTleaf/dE for the collar
// solve's implicit-function theorem. Linearised about Tair (Monteith & Unsworth;
// the step plantecophys takes) it becomes g_rad in the DENOMINATOR, so Tleaf stays
// affine in E. Do not move any Tleaf dependence into Rn_. The error is of order
// (Tleaf - Tair)/Tair_K, about 1% at 20 K.
inline constexpr double stefan_boltzmann = 5.670374419e-8;  // W m^-2 K^-4
// Longwave emissivity of a leaf. Thermal-infrared leaf emissivities are measured
// at 0.94-0.99; 0.97 is the usual working value.
inline constexpr double leaf_emissivity = 0.97;
// Clear-sky atmospheric emissivity, Brutsaert (1975) as Campbell & Norman write
// it: ema = 0.642*(ea/Tair_K)^(1/7), held below 1 by a soft minimum (see
// set_physiology).
//
// ⚠️ ea MUST BE IN PASCALS. The same relation is more often quoted as
// `1.24*(ea/T)^(1/7)` with ea in hPa; the two agree because 0.642*100^(1/7) =
// 1.2395. Feeding kPa to 0.642 gives ~0.31 where the answer is ~0.82: finite,
// plausible, and it triples the longwave loss. Change the unit and the
// coefficient together, or not at all.
inline constexpr double atmos_emissivity_coef = 0.642;
inline constexpr double atmos_emissivity_exponent = 1.0 / 7.0;
// Width of the soft minimum that holds Brutsaert's emissivity below 1,
// ema - w*softplus((ema - 1)/w). Where the fit passes unity -- hot,
// humid air -- the sky radiates as a blackbody at Tair and the longwave term
// goes to zero smoothly rather than at a corner. Bit-identical to the plain
// fit wherever ema < 0.937, where the correction is below half an ulp.
inline constexpr double atmos_emissivity_softmin_width = 0.002;
// fixed aerodynamic resistance fallback, s m^-1 (doc 6/7.4; used when the wind
// model is unavailable, e.g. a bare Leaf with no wind/d set)
inline constexpr double aerodynamic_resistance_fixed = 50.0;
// leaf boundary-layer coefficient C_ra, s^0.5 m^-1, in ra = C_ra*sqrt(d/U) (doc 4.1)
inline constexpr double aerodynamic_resistance_coef = 200.0;
// Floor on the leaf-to-air vapour pressure deficit (kPa) used by Fick's law.
// A leaf cooled below the dew point has a NEGATIVE deficit, and
// dividing a positive transpiration by it reports a negative stomatal
// conductance. Flooring keeps such points finite and unattractive rather than
// sign-inverted. 0.01 kPa is far below any operating deficit -- the driver
// defaults are 1.5-2 kPa -- so it binds only where the model is already outside
// what a one-way diffusion equation describes.
// ⚠️ THE ONE FLOOR KEPT ON A SOLVED PATH BESIDE THE STOMATAL CORNER, by
// decision. Measured, it binds at candidate points inside the energy-balance
// stem routes' domain and at a solved point only where SOX, JW26 and LeastCost
// pin a dark leaf at psi_crit; never at a collar-route optimum.
inline constexpr double vpd_leaf_min = 0.01;


} // namespace phylloptim

#endif
