// -*-c++-*-
#ifndef PHYLLOPTIM_CONSTANTS_HPP_
#define PHYLLOPTIM_CONSTANTS_HPP_

// Physical and physiological constants for the leaf model, as `inline constexpr`:
// one entity per program, usable in constant expressions. Names are spelled out
// because they sit at namespace scope in a public header (the gas constant is
// `gas_constant`, never `R`). tests/testthat/test-registry.R parses the
// `inline constexpr double NAME = VALUE` lines against inst/model-registry.yml.

#include <cmath>

namespace phylloptim {

// Vcmax and Jmax temperature responses, within the ranges of Medlyn et al. (2002).
// J mol^-1
inline constexpr double vcmax_ha = 60000;
// J mol^-1
inline constexpr double vcmax_H_d = 200000;
// J mol^-1 K^-1
inline constexpr double vcmax_d_S = 650;

// J mol^-1
inline constexpr double jmax_ha = 30000;
// J mol^-1
inline constexpr double jmax_H_d = 200000;
// J mol^-1 K^-1
inline constexpr double jmax_d_S = 650;

// Gamma*, Kc and Ko at 25 C and their activation energies (Bernacchi et al. 2001).
// umol mol^-1
inline constexpr double gamma_25 = 42.75;
// J mol^-1
inline constexpr double gamma_ha = 37.83e3;

// umol mol^-1
inline constexpr double kc_25 = 404.9 ;
// J mol^-1
inline constexpr double kc_ha = 79.43e3;

// umol mol^-1
inline constexpr double ko_25 = 278400 ;
// J mol^-1
inline constexpr double ko_ha = 36.38e3;

// umol mol^-1 -> Pa is not a constant: it is 1e-6 * P, `Leaf::umol_per_mol_to_Pa_`,
// derived from atm_kpa_. Do not hard-code it at 101.3 kPa.

// The molar mass of water, 18.015 g/mol, written once: both conversions derive
// from it, so they are reciprocal by construction. plant `using`-declares both
// names.
inline constexpr double molar_mass_h2o = 0.018015; // kg mol^-1
// kg mol^-1, for converting a molar water flux back to kg.
inline constexpr double kg_per_mol_h2o = molar_mass_h2o;
// mol H2O kg^-1
inline constexpr double kg_to_mol_h2o = 1.0 / molar_mass_h2o;
// mol mol ^-1 / (umol mol ^-1)
inline constexpr double umol_to_mol = 1e-6;
// Pa kPa^-1
inline constexpr double kPa_to_Pa = 1000.0;

// Universal gas constant, J mol^-1 K^-1, exact in the SI (R = N_A k_B). Read as
// Ea/(R*T) with Ea/RT ~ 24, so a relative change reaches the rates ~24x amplified.
inline constexpr double gas_constant = 8.314462618153240;

//convert deg C to deg K
inline constexpr double C_to_K = 273.15;

// H2O:CO2 stomatal diffusion ratio (the g1 literature uses 1.6).
inline constexpr double H2O_CO2_stom_diff_ratio = 1.67;

// Boundary-layer conductance to water vapour over that to heat, for a
// hypostomatous leaf: plantecophys's `Gbw = 1.075 * Gbh` (LeafEnergyBalance, at
// StomatalRatio 1).
inline constexpr double bl_vapour_heat_ratio = 1.075;

// Soft floor, and its width, under the boundary layer's split of the CO2 path,
// d = 1 - (r_s - r_b) g_c,hat / g_bw (Leaf::boundary_layer_split_CO2). Wherever
// the total conductance is below the boundary layer's own, d >= r_b / r_s (0.84 at
// the default ratio), so the split is the plain series there to the last bit; past
// it the floor keeps the conductance finite instead of passing a pole.
inline constexpr double bl_co2_split_floor = 0.1;
inline constexpr double bl_co2_split_floor_width = 0.01;

// rho_w * g, as a head in MPa per metre of depth. 6.8e-04 relative below
// 1000 * 9.80665 / 1e6, which moves results; correcting it needs a water-density
// convention first (1000, 4 C and 25 C span 0.3%), a modelling decision.
inline constexpr double gravity_head = 9.8e-3; // MPa / m

// --- Penman-Monteith leaf energy balance -------------------------------------
// Read by Leaf::leaf_temp_from_E and the es/Delta helpers, only on the
// use_energy_balance_ path.
// Latent heat of vaporisation of water, J kg^-1: the value at about 21.5 C, not
// 25 C (2.442e6). It should become a function of leaf temperature, not a better
// constant.
inline constexpr double latent_heat_vap = 2.45e6;
// volumetric heat capacity of air, J m^-3 K^-1
inline constexpr double vol_heat_cap_air = 1200.0;
// Photons per joule of PAR for sunlight, umol J^-1; a property of the spectrum (a
// red/blue LED is ~5-5.5). Read only by the energy balance, to turn the
// intercepted photon flux back into PAR energy.
inline constexpr double umol_par_per_joule = 4.57;
// --- Longwave radiation ------------------------------------------------------
//
// Net longwave is an isothermal net radiation (the leaf at air temperature) plus a
// radiative conductance for the leaf's departure from it. Absorbed shortwave is
// each intercepted band times the leaf's absorptance (alpha_PAR_, alpha_NIR_):
//
//   Rn_iso = absorbed shortwave - (1 - ema)*sigma*Tair_K^4
//   g_rad  = 4*eps_leaf*sigma*Tair_K^3
//   Tleaf  = Tair + (Rn_iso - lambda*E) / (vol_heat_cap_air/ra + g_rad)
//
// Linearising emission about Tair (Monteith & Unsworth, as plantecophys does) keeps
// Tleaf affine in E, so the balance stays explicit and dTleaf/dE is constant for
// the collar solve. Do not move any Tleaf dependence into Rn_. The error is of
// order (Tleaf - Tair)/Tair_K, about 1% at 20 K.
inline constexpr double stefan_boltzmann = 5.670374419e-8;  // W m^-2 K^-4
// Longwave emissivity of a leaf. Thermal-infrared leaf emissivities are measured
// at 0.94-0.99; 0.97 is the usual working value.
inline constexpr double leaf_emissivity = 0.97;
// Clear-sky atmospheric emissivity, Brutsaert (1975) as Campbell & Norman write
// it: ema = 0.642*(ea/Tair_K)^(1/7), held below 1 by a soft minimum. ea must be
// in Pa: the hPa form is 1.24*(ea/T)^(1/7), and kPa gives ~0.31 where the answer
// is ~0.82, tripling the longwave loss. Change unit and coefficient together.
inline constexpr double atmos_emissivity_coef = 0.642;
inline constexpr double atmos_emissivity_exponent = 1.0 / 7.0;
// Width w of the soft minimum ema - w*softplus((ema - 1)/w) that holds the
// emissivity below 1, so the longwave term reaches zero smoothly in hot, humid
// air. Bit-identical to the plain fit wherever ema < 0.937.
inline constexpr double atmos_emissivity_softmin_width = 0.002;
// Aerodynamic resistance, s m^-1, where no leaf dimension or wind speed is set.
inline constexpr double aerodynamic_resistance_fixed = 50.0;
// leaf boundary-layer coefficient C_ra, s^0.5 m^-1, in ra = C_ra*sqrt(d/U)
inline constexpr double aerodynamic_resistance_coef = 200.0;
// Floor on the leaf-to-air deficit Fick's law divides by, kPa. A leaf cooled below
// the dew point has a negative deficit and would report a negative conductance;
// the floor keeps such points finite and unattractive. It binds at candidate
// points on the energy-balance stem routes and at a solved point only where
// LeastCost pins a dark leaf at psi_crit, never at a collar-route optimum.
inline constexpr double vpd_leaf_min = 0.01;


} // namespace phylloptim

#endif
