# The ATLS thermal-damage layer (traitecoevo/plant#566), behind
# `use_thermal_damage_`.
#
# A lasting-damage model on three fractions of leaf area, functional <=> impaired
# -> damaged (F = 1 - I - D), plus an acclimation state A. The states are seated
# as A_, I_ and D_; the solve scales both capacities by F and reports the
# impairment and repair rates (k_imp_eff_, k_rec_eff_) at the leaf's temperature;
# $damage_rates() and leaf_damage_trajectory() are the dynamics. The structure is
# plant's July layer with its recoverable and permanent damage renamed I and D.
#
# ⚠️ TWO KINDS OF ASSERTION, and both are needed. Off, the layer must be INERT,
# which is asserted with identical() against inputs set to values that would
# matter if they were read. On, it must DO what it says, which is asserted
# against the equations at non-default values -- otherwise an inert gate would
# pass the first half by doing nothing at all.

# The gas constant and Kelvin offset the leaf uses. A reference built on R = 8.314
# would disagree with the impairment rate by more than rounding.
.R_gas <- 8.314462618153240
.T0 <- 273.15

thermal_leaf <- function(supply = leaf_supply_singlelayer()) {
  leaf_model(traits = leaf_traits(), control = leaf_control(), supply = supply)
}

drive_at <- function(l, leaf_temp, psi_soil = 0.5) {
  # Daylight NIR for the PPFD below; read only with the energy balance on.
  l$NIR_int_ <- 1200 / 4.57
  if (identical(l$supply_kind, "single")) {
    set_drivers(l, psi_soil = psi_soil, PPFD = 1200,
                root_network = series_resistance(50),
                leaf_specific_conductance_max = 1e-4, atm_vpd = 1.0, ca = 40,
                leaf_temp = leaf_temp, atm_o2_kpa = 21, atm_kpa = 101.3)
  } else {
    set_drivers(l, psi_soil = psi_soil + c(0, 0.4, 0.8), PPFD = 1200,
                leaf_specific_conductance_max = 1e-4, atm_vpd = 1.0, ca = 40,
                leaf_temp = leaf_temp, atm_o2_kpa = 21, atm_kpa = 101.3)
  }
  invisible(l)
}

# A rate on the layer's centred Arrhenius: one at the reference temperature
# Tref (deg C), rising with activation energy Ea.
arrh_ref <- function(leaf_temp, Ea, Tref) {
  exp((Ea / .R_gas) * (1 / (Tref + .T0) - 1 / (leaf_temp + .T0)))
}

# The default impairment rate, per day: Bison & Michaletz (2026)'s 177-species
# fit, the rate common to every species at the isokinetic temperature and their
# mean activation energy.
k_imp_ref <- function(leaf_temp, k = 276, Ea = 316e3, Tiso = 50.5) {
  k * arrh_ref(leaf_temp, Ea, Tiso)
}

# Every input of the layer, at a value other than its default.
atls_inputs <- list(topt_offset_ = 2.5, k_i_ = 5, k_i_Ea_ = 2.5e5,
                    k_i_Tiso_ = 47, k_i_baseline_ = 8, k_rec_ = 0.4,
                    k_rec_Ea_ = 3e4, dTopt_max_ = 4, dTimpair_max_ = 2.5,
                    K_A_ = 2, A_ = 1.5, I_ = 0.2, D_ = 0.15, k_mat_ = 0.05,
                    k_mat_Ea_ = 2e4, acclim_alpha_ = 0.03, acclim_beta_ = 0.2,
                    acclim_T_ = 28, acclim_s_ = 2, c_acclim_maint_ = 0.05,
                    c_repair_maint_ = 0.3, c_protect_maint_ = 0.4,
                    c_repair_flux_ = 0.8)

set_inputs <- function(l, inputs = atls_inputs) {
  for (nm in names(inputs)) l[[nm]] <- inputs[[nm]]
  invisible(l)
}

# --- ported from plant's July test-leaf-thermal.R ---------------------------

test_that("thermal damage is off by default and writes zero outputs", {
  l <- thermal_leaf()
  expect_false(l$use_thermal_damage_)
  drive_at(l, 30)
  expect_identical(l$k_imp_eff_, 0)
  expect_identical(l$k_rec_eff_, 0)
  expect_identical(l$functional_, 1)

  # On but undamaged and unshifted: the shift is an exact zero in exact
  # arithmetic, and d_S_shifted's round trip H_d/(H_d/(d_S - corr)) + corr
  # returns d_S to rounding only, hence a tolerance rather than identity.
  off <- drive_at(thermal_leaf(), 45)
  on <- thermal_leaf()
  on$use_thermal_damage_ <- TRUE
  drive_at(on, 45)
  expect_equal(on$vcmax_, off$vcmax_, tolerance = 1e-9)
  expect_equal(on$jmax_, off$jmax_, tolerance = 1e-9)
})

test_that("the impairment rate is Bison and Michaletz's centred Arrhenius", {
  for (temp in c(20, 35, 45, 55)) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    drive_at(l, temp)
    expect_equal(l$k_imp_eff_, k_imp_ref(temp), tolerance = 1e-12,
                 info = paste("leaf_temp", temp))
  }
  # At the isokinetic temperature the exponent is an exact zero.
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  drive_at(l, l$k_i_Tiso_)
  expect_identical(l$k_imp_eff_, l$k_i_)
  # The defaults are their fit: exp(b) = 3.2e-3 per second, at 50.5 C, with the
  # mean Ea of 316 kJ/mol.
  expect_equal(l$k_i_ / 86400, 3.2e-3, tolerance = 2e-3)
  expect_identical(c(l$k_i_Ea_, l$k_i_Tiso_), c(316e3, 50.5))
})

# Their Fig. 2d, as first-order damage 1 - exp(-k t) at a constant temperature:
# the temperature that half-damages a leaf in a given exposure. Nicole Bison's
# plot of these curves, from the paper's parameters, puts them at about 54, 47,
# 43 and 35 C for one minute, fifteen minutes, an hour and a day.
test_that("the default impairment rate reproduces their exposure curves", {
  l <- thermal_leaf()
  k_at <- function(temp) l$damage_rates(temp, 0, 0, 0)[2]  # dI/dt at F = 1
  T50 <- function(days) {
    stats::uniroot(function(temp) 1 - exp(-k_at(temp) * days) - 0.5,
                   c(10, 80), tol = 1e-10)$root
  }
  expect_lt(abs(T50(1 / 1440) - 54.0), 0.5)
  expect_lt(abs(T50(15 / 1440) - 46.6), 0.5)
  expect_lt(abs(T50(1 / 24) - 42.9), 0.5)
  expect_lt(abs(T50(1) - 34.8), 0.5)
})

test_that("a negative impairment activation energy is refused", {
  l <- thermal_leaf()
  l$k_i_Ea_ <- -1
  expect_silent(drive_at(l, 30))   # gate off: not read
  l$use_thermal_damage_ <- TRUE
  expect_error(drive_at(l, 30), "needs k_i_Ea_ >= 0")
})

# ⚠️ THE DECISION THIS LAYER ENCODES. The Vcmax and Jmax curves are the
# instantaneous, reversible response and impairment is a separate process, so
# they share no parameter: moving one must not move the other.
test_that("the impairment rate is independent of the Jmax curve", {
  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  drive_at(base, 38)

  jm <- thermal_leaf()
  jm$use_thermal_damage_ <- TRUE
  jm$jmax_H_d_ <- 210000
  jm$jmax_d_S_ <- 670
  drive_at(jm, 38)
  expect_false(isTRUE(all.equal(jm$jmax_, base$jmax_)))  # the curve moved
  expect_identical(jm$k_imp_eff_, base$k_imp_eff_)       # the rate did not

  im <- thermal_leaf()
  im$use_thermal_damage_ <- TRUE
  im$k_i_Tiso_ <- 45
  im$k_i_Ea_ <- 3e5
  drive_at(im, 38)
  expect_false(isTRUE(all.equal(im$k_imp_eff_, base$k_imp_eff_)))  # it moved
  expect_identical(c(im$vcmax_, im$jmax_), c(base$vcmax_, base$jmax_))
})

test_that("each rate is exact at its reference and a centred Arrhenius off it", {
  # k_i_Ea_ = 0: impairment is exactly k_i_ at every temperature.
  for (temp in c(5, 25, 45)) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    l$k_i_Ea_ <- 0
    drive_at(l, temp)
    expect_identical(l$k_imp_eff_, l$k_i_, info = paste("leaf_temp", temp))
  }
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  l$k_i_ <- 0.3
  l$k_i_Ea_ <- 1.5e5
  l$k_i_Tiso_ <- 35
  drive_at(l, 42)
  expect_equal(l$k_imp_eff_, 0.3 * arrh_ref(42, 1.5e5, 35), tolerance = 1e-14)
  # The sharpness Ea sets: tenfold between T1 and the T2 with
  # 1/T2 = 1/T1 - R*ln(10)/Ea, about ln(10)*R*T^2/Ea degrees above it.
  T2 <- 1 / (1 / (35 + .T0) - .R_gas * log(10) / 1.5e5) - .T0
  drive_at(l, T2)
  expect_equal(l$k_imp_eff_ / 0.3, 10, tolerance = 1e-12)
  # Bounded: the exponent is clamped at 30, so an absurd Ea stays finite.
  l$k_i_Ea_ <- 1e8
  drive_at(l, 90)
  expect_equal(l$k_imp_eff_, 0.3 * exp(30), tolerance = 1e-14)
})

test_that("acclimation moves the isokinetic temperature by its own amount", {
  # Acclimated, with the impairment's own maximum shift switched off: the capacity
  # optimum moves, the impairment rate does not.
  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  drive_at(base, 38)
  acc <- thermal_leaf()
  acc$use_thermal_damage_ <- TRUE
  acc$A_ <- 3
  acc$dTimpair_max_ <- 0
  drive_at(acc, 38)
  expect_false(isTRUE(all.equal(acc$vcmax_, base$vcmax_)))
  expect_equal(acc$k_imp_eff_, base$k_imp_eff_, tolerance = 1e-15)

  # And the other way round.
  only <- thermal_leaf()
  only$use_thermal_damage_ <- TRUE
  only$A_ <- 3
  only$dTopt_max_ <- 0
  drive_at(only, 38)
  expect_lt(only$k_imp_eff_, base$k_imp_eff_)
  expect_equal(only$vcmax_, base$vcmax_, tolerance = 1e-12)

  # At the defaults the two maxima are equal, so one shift moves both by the
  # same number of degrees, as in the port.
  expect_identical(base$dTimpair_max_, base$dTopt_max_)
})

test_that("impaired and damaged area discount BOTH capacities by F = 1 - I - D", {
  undamaged <- thermal_leaf()
  undamaged$use_thermal_damage_ <- TRUE
  drive_at(undamaged, 30)
  damaged <- thermal_leaf()
  damaged$use_thermal_damage_ <- TRUE
  damaged$I_ <- 0.3
  damaged$D_ <- 0.1
  drive_at(damaged, 30)
  expect_equal(damaged$vcmax_ / undamaged$vcmax_, 0.6, tolerance = 1e-14)
  expect_equal(damaged$jmax_ / undamaged$jmax_, 0.6, tolerance = 1e-14)

  # Past total damage the functional fraction sits on its floor, so the capacities
  # stay strictly positive and the solve stays well-conditioned.
  crushed <- thermal_leaf()
  crushed$use_thermal_damage_ <- TRUE
  crushed$I_ <- 0.8
  crushed$D_ <- 0.5
  drive_at(crushed, 30)
  expect_equal(crushed$vcmax_ / undamaged$vcmax_, 1e-6, tolerance = 1e-12)
  expect_equal(crushed$jmax_ / undamaged$jmax_, 1e-6, tolerance = 1e-12)
  # A leaf this damaged cannot cover its respiration, so the solve takes the
  # no-flow exit a dark leaf takes -- and completes there, with finite carbon:
  # what its remaining capacity fixes at ambient CO2, less respiration.
  # (lambda_emergent is NaN at that exit by design, as it is for a dark leaf.)
  crushed$find_root_collar_psi()
  expect_true(is.finite(crushed$profit_))
  expect_identical(crushed$assim_colimited_,
                   crushed$assim_colimited(crushed$ca_))
  expect_lt(crushed$assim_colimited_, 0)
})

test_that("a thermostability offset moves the whole response up in temperature", {
  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  drive_at(base, 40)
  tolerant <- thermal_leaf()
  tolerant$use_thermal_damage_ <- TRUE
  tolerant$topt_offset_ <- 8
  drive_at(tolerant, 40)
  expect_gt(tolerant$vcmax_, base$vcmax_)  # the capacity optimum moved up
  expect_lt(tolerant$k_imp_eff_, base$k_imp_eff_)  # and impairment slowed
})

test_that("repair and maturation are centred Arrhenius rates at 25 C", {
  # Default activation energies are zero: the rates are exactly k_rec_ and k_mat_
  # at every temperature, and nothing shuts repair off when hot.
  for (temp in c(20, 45, 60)) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    drive_at(l, temp)
    expect_identical(l$k_rec_eff_, l$k_rec_, info = paste("leaf_temp", temp))
    expect_identical(l$damage_rates(temp, 0, 1, 0)[3], l$k_mat_,
                     info = paste("leaf_temp", temp))
  }
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  l$k_rec_Ea_ <- 3e4
  l$k_mat_Ea_ <- 2e4
  drive_at(l, 25)
  expect_identical(l$k_rec_eff_, l$k_rec_)
  expect_identical(l$damage_rates(25, 0, 1, 0)[3], l$k_mat_)
  drive_at(l, 40)
  expect_equal(l$k_rec_eff_, l$k_rec_ * arrh_ref(40, 3e4, 25), tolerance = 1e-14)
  expect_equal(l$damage_rates(40, 0, 1, 0)[3], l$k_mat_ * arrh_ref(40, 2e4, 25),
               tolerance = 1e-14)
})

test_that("the impairment rate is finite, positive and rising across a wide sweep", {
  temps <- seq(-30, 90, by = 1)
  k <- vapply(temps, function(temp) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    drive_at(l, temp)
    l$k_imp_eff_
  }, numeric(1))
  expect_true(all(is.finite(k)))
  expect_true(all(k > 0))
  # Rising wherever the exponent is inside its clamp, flat where it is not.
  expect_true(all(diff(k) >= 0))
  expect_true(all(diff(k)[temps[-1] > 0] > 0))
})

test_that("outlandish inputs still give finite capacities and outputs", {
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  l$topt_offset_ <- 50
  l$I_ <- 0.4
  l$D_ <- 0.4
  drive_at(l, 99)
  l$k_i_Ea_ <- 5e5
  l$k_i_Tiso_ <- 200
  drive_at(l, 99)
  expect_true(all(is.finite(c(l$vcmax_, l$jmax_, l$k_imp_eff_, l$functional_,
                              l$k_rec_eff_, l$electron_transport_))))
  expect_true(l$k_imp_eff_ >= 0)
  expect_true(all(is.finite(l$damage_rates(99, 50, 0.4, 0.4))))
})

# On the energy-balance path the block is re-seated at each candidate's leaf
# temperature, so the impairment rate follows the operating point: more
# transpiration, a cooler leaf, slower impairment. That is the avoidance axis.
test_that("on the energy-balance path impairment is taken at the operating Tleaf", {
  pm_leaf <- function(gate) {
    l <- thermal_leaf()
    l$use_energy_balance_ <- TRUE
    l$use_thermal_damage_ <- gate
    drive_at(l, 42, psi_soil = 0.2)
    l
  }
  lo <- pm_leaf(TRUE)
  lo$set_leaf_states_rates_from_psi_stem(1.0, 0.2)
  hi <- pm_leaf(TRUE)
  hi$set_leaf_states_rates_from_psi_stem(3.5, 0.2)
  expect_gt(hi$transpiration_, lo$transpiration_)
  expect_lt(hi$Tleaf_, lo$Tleaf_)
  expect_lt(hi$k_imp_eff_, lo$k_imp_eff_)
  # At THAT temperature, not at the air temperature driver.
  expect_equal(hi$k_imp_eff_, k_imp_ref(hi$Tleaf_), tolerance = 1e-12)
  expect_equal(lo$k_imp_eff_, k_imp_ref(lo$Tleaf_), tolerance = 1e-12)

  off <- pm_leaf(FALSE)
  off$set_leaf_states_rates_from_psi_stem(3.5, 0.2)
  expect_identical(off$k_imp_eff_, 0)
})

# --- respiration: phylloptim's R_d, plus the costs --------------------------
#
# ⚠️ THE ONE INTENDED DIVERGENCE FROM THE JULY LAYER. There R_d was 0.015 *
# Vcmax after the damage discount, so a damaged leaf respired less. Here R_d is
# the R_d_25 trait on its Q10, for every model, and the layer only ADDS its
# costs: damage does not lower respiration.
test_that("R_d is not scaled by survival, and the costs add to it", {
  off <- drive_at(thermal_leaf(), 30)
  damaged <- thermal_leaf()
  damaged$use_thermal_damage_ <- TRUE
  damaged$I_ <- 0.3
  damaged$D_ <- 0.1
  drive_at(damaged, 30)
  # Coefficients zero: the four terms add an exact zero.
  expect_identical(damaged$R_d_, off$R_d_)

  cost_of <- function(...) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    set_inputs(l, list(...))
    drive_at(l, 30)
    l
  }
  acc <- cost_of(A_ = 3, c_acclim_maint_ = 0.1)
  expect_equal(acc$R_d_ - off$R_d_, 0.1 * 3, tolerance = 1e-12)
  # The standing repair cost is paid on the CAPACITY, even when cold.
  rep <- cost_of(c_repair_maint_ = 1e-2)
  expect_equal(rep$R_d_ - off$R_d_, 1e-2 * rep$k_rec_, tolerance = 1e-12)
  prot <- cost_of(k_i_ = 0.2, c_protect_maint_ = 0.5)
  expect_equal(prot$R_d_ - off$R_d_, 0.5 * (prot$k_i_baseline_ - 0.2),
               tolerance = 1e-12)
  # Protection above the baseline costs nothing extra: max(0, k_i_baseline - k_i).
  none <- cost_of(k_i_ = 300, c_protect_maint_ = 0.5)
  expect_identical(none$R_d_, off$R_d_)
  flux <- cost_of(I_ = 0.4, c_repair_flux_ = 1.0)
  expect_equal(flux$R_d_ - off$R_d_, flux$k_rec_eff_ * 0.4, tolerance = 1e-12)
  expect_gt(flux$R_d_, off$R_d_)

  # All four at once, with every input at a non-default value.
  all4 <- do.call(cost_of, atls_inputs)
  expected <- with(atls_inputs,
                   c_acclim_maint_ * A_ + c_repair_maint_ * k_rec_ +
                     c_protect_maint_ * (k_i_baseline_ - k_i_) +
                     c_repair_flux_ * all4$k_rec_eff_ * I_)
  expect_equal(all4$R_d_ - off$R_d_, expected, tolerance = 1e-12)
})

# --- off: inert, with every input set -----------------------------------------

test_that("with the gate off, no thermal-damage input is read", {
  # Every input at a value that would move the answer if it were read, against a
  # leaf that never had them set. identical(), on every reported output, across
  # both supply paths, both temperature paths and two temperatures.
  checked <- 0L
  for (supply in list(leaf_supply_singlelayer(), leaf_supply_multilayer())) {
    for (eb in c(FALSE, TRUE)) {
      for (temp in c(25, 40)) {
        clean <- thermal_leaf(supply)
        clean$use_energy_balance_ <- eb
        drive_at(clean, temp)
        clean$find_root_collar_psi()

        set <- set_inputs(thermal_leaf(supply))
        set$use_energy_balance_ <- eb
        drive_at(set, temp)
        set$find_root_collar_psi()

        at <- paste(supply$kind, "eb", eb, "T", temp)
        expect_identical(set$operating_point_values(),
                         clean$operating_point_values(), info = at)
        expect_identical(c(set$vcmax_, set$jmax_, set$R_d_),
                         c(clean$vcmax_, clean$jmax_, clean$R_d_), info = at)
        psi <- clean$opt_root_psi_
        expect_identical(set$dprofit_droot_collar_psi(psi),
                         clean$dprofit_droot_collar_psi(psi), info = at)
        expect_identical(c(set$k_imp_eff_, set$k_rec_eff_, set$functional_),
                         c(0, 0, 1), info = at)

        # ...and the same inputs DO move it with the gate on, or the identity
        # above could be an inert gate rather than an inert layer.
        set$use_thermal_damage_ <- TRUE
        drive_at(set, temp)
        set$find_root_collar_psi()
        expect_false(identical(set$operating_point_values(),
                               clean$operating_point_values()), info = at)
        checked <- checked + 1L
      }
    }
  }
  expect_identical(checked, 8L)
})

# --- on: what it does, at non-default values ----------------------------------

test_that("the entropy shift moves each capacity optimum by the requested amount", {
  # Topt = H_d/(d_S - R*ln(Ea/(H_d - Ea))), so the argmax of vcmax_(T) and of
  # jmax_(T) is known in closed form; locate it numerically through the model.
  argmax_T <- function(l, field) {
    stats::optimize(function(temp) {
      drive_at(l, temp)
      l[[field]]
    }, c(15, 55), maximum = TRUE, tol = 1e-9)$maximum
  }
  topt <- function(ha, H_d, d_S) H_d / (d_S - .R_gas * log(ha / (H_d - ha))) - .T0

  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  shifted <- set_inputs(thermal_leaf())
  shifted$use_thermal_damage_ <- TRUE
  # 2.5 + 4 * 1.5 / (1.5 + 2) = 4.214...
  delta <- with(atls_inputs, topt_offset_ + dTopt_max_ * A_ / (A_ + K_A_))

  for (field in c("vcmax_", "jmax_")) {
    pre <- if (field == "vcmax_") "vcmax" else "jmax"
    T0 <- topt(base[[paste0(pre, "_ha_")]], base[[paste0(pre, "_H_d_")]],
               base[[paste0(pre, "_d_S_")]])
    expect_equal(argmax_T(base, field), T0, tolerance = 1e-5, info = field)
    expect_equal(argmax_T(shifted, field) - argmax_T(base, field), delta,
                 tolerance = 1e-5, info = field)
  }

  # The functional fraction is applied on top of the shifted curve.
  surv <- 1 - atls_inputs$I_ - atls_inputs$D_
  undamaged <- set_inputs(thermal_leaf(),
                          utils::modifyList(atls_inputs,
                                            list(I_ = 0, D_ = 0)))
  undamaged$use_thermal_damage_ <- TRUE
  drive_at(undamaged, 33)
  drive_at(shifted, 33)
  expect_equal(shifted$vcmax_ / undamaged$vcmax_, surv, tolerance = 1e-14)
  expect_equal(shifted$jmax_ / undamaged$jmax_, surv, tolerance = 1e-14)
})

test_that("the isokinetic temperature moves with its own shift", {
  l <- set_inputs(thermal_leaf())
  l$use_thermal_damage_ <- TRUE
  delta <- with(atls_inputs, topt_offset_ + dTimpair_max_ * A_ / (A_ + K_A_))
  # Moved by the constitutive offset and its own acclimated shift, the rate there
  # is k_i_ again.
  drive_at(l, l$k_i_Tiso_ + delta)
  expect_equal(l$k_imp_eff_, l$k_i_, tolerance = 1e-12)
  drive_at(l, 40)
  expect_equal(l$k_imp_eff_,
               l$k_i_ * arrh_ref(40, l$k_i_Ea_, l$k_i_Tiso_ + delta),
               tolerance = 1e-12)
})

test_that("the acclimation forcing is the guarded softplus", {
  l <- thermal_leaf()
  expect_equal(l$softplus_(0, 1), log(2), tolerance = 1e-15)
  expect_equal(l$softplus_(0.7, 2.5), log1p(exp(2.5 * 0.7)) / 2.5,
               tolerance = 1e-15)
  # Clamped at s*x = 30, so it does not overflow.
  expect_true(is.finite(l$softplus_(1e6, 1)))
})

# --- the temperature cache ----------------------------------------------------
#
# ⚠️ THE INVARIANT A LATER EDIT COULD BREAK. The cache key does not hold the
# layer's inputs, because set_physiology bypasses the cache whenever the gate is
# on: A_, I_ and D_ are state a caller moves every step. Remove the bypass and
# the second drive below takes a cache hit and keeps the first damage level.
test_that("with the gate on, every re-drive reads the current state", {
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  drive_at(l, 30)
  l$find_root_collar_psi()
  vcmax_first <- l$vcmax_
  A_first <- l$assim_colimited_

  l$I_ <- 0.25                      # same drivers, new damage state
  drive_at(l, 30)
  l$find_root_collar_psi()
  expect_equal(l$vcmax_ / vcmax_first, 0.75, tolerance = 1e-14)
  expect_lt(l$assim_colimited_, A_first)

  l$A_ <- 2                           # and new acclimation state
  vcmax_damaged <- l$vcmax_
  drive_at(l, 30)
  expect_false(identical(l$vcmax_, vcmax_damaged))

  # Off again at the same drivers: the gate is in the key, so this recomputes
  # rather than serving the gate-on block, and matches a leaf never switched on.
  l$use_thermal_damage_ <- FALSE
  drive_at(l, 30)
  l$find_root_collar_psi()
  never <- drive_at(thermal_leaf(), 30)
  never$find_root_collar_psi()
  expect_identical(l$operating_point_values(), never$operating_point_values())
  expect_identical(c(l$k_imp_eff_, l$k_rec_eff_, l$functional_), c(0, 0, 1))
})

# --- the dynamics ---------------------------------------------------------------
#
# $damage_rates() is the right-hand side everything integrates, and
# leaf_damage_trajectory() integrates it: the route for fitting the layer to a
# leaf experiment with this package alone.

test_that("damage_rates is the three-state model, at every input non-default", {
  l <- set_inputs(thermal_leaf())
  for (temp in c(20, 38, 47)) {
    # A = 1.5 would put dA/dt at an exact cancellation at 38 C.
    A <- 1.2; I <- 0.2; D <- 0.15
    delta_i <- with(atls_inputs, topt_offset_ + dTimpair_max_ * A / (A + K_A_))
    k_imp <- with(atls_inputs, k_i_ * arrh_ref(temp, k_i_Ea_, k_i_Tiso_ + delta_i))
    k_rec <- with(atls_inputs, k_rec_ * arrh_ref(temp, k_rec_Ea_, 25))
    k_mat <- with(atls_inputs, k_mat_ * arrh_ref(temp, k_mat_Ea_, 25))
    # softplus as the leaf computes it: its argument s*x is clamped at 30.
    sp <- with(atls_inputs,
               log1p(exp(min(max(acclim_s_ * (temp - acclim_T_), -30), 30))) / acclim_s_)
    expected <- with(atls_inputs, c(
      acclim_alpha_ * sp - acclim_beta_ * A,
      k_imp * (1 - I - D) - (k_rec + k_mat) * I,
      k_mat * I))
    expect_equal(l$damage_rates(temp, A, I, D), expected, tolerance = 1e-13,
                 info = paste("leaf_temp", temp))
  }
})

test_that("damage_rates reads its arguments, not the seated state or the gate", {
  a <- set_inputs(thermal_leaf())
  b <- set_inputs(thermal_leaf())
  b$A_ <- NaN; b$I_ <- NaN; b$D_ <- NaN
  b$use_thermal_damage_ <- TRUE
  expect_identical(b$damage_rates(40, 1, 0.1, 0.05),
                   a$damage_rates(40, 1, 0.1, 0.05))
  expect_true(all(is.finite(b$damage_rates(40, 1, 0.1, 0.05))))
})

test_that("the trajectory matches the exact solution at a constant temperature", {
  # With A held at its equilibrium, the (I, D) system is linear with constant
  # coefficients, so the matrix exponential is the exact answer.
  l <- set_inputs(thermal_leaf())
  temp <- 41
  sp <- with(atls_inputs, log1p(exp(acclim_s_ * (temp - acclim_T_))) / acclim_s_)
  A_eq <- with(atls_inputs, acclim_alpha_ / acclim_beta_ * sp)
  r0 <- l$damage_rates(temp, A_eq, 0, 0)
  expect_equal(r0[1], 0, tolerance = 1e-14)
  k <- r0[2]                              # k_imp at F = 1
  r1 <- l$damage_rates(temp, A_eq, 1, 0)  # dI = -(k_rec + k_mat), dD = k_mat
  out_rate <- -r1[2]                      # k_rec + k_mat
  kmat <- r1[3]
  # y' = M y + b, y = (I, D)
  M <- matrix(c(-(k + out_rate), kmat, -k, 0), 2, 2)
  b <- c(k, 0)
  times <- c(0, 0.5, 2, 7, 30)
  ev <- eigen(M)
  exact <- t(vapply(times, function(t) {
    E <- Re(ev$vectors %*% diag(exp(ev$values * t)) %*% solve(ev$vectors))
    as.vector((E - diag(2)) %*% solve(M, b))
  }, numeric(2)))
  d <- leaf_damage_trajectory(l, times, temp, A0 = A_eq)
  expect_equal(d$A, rep(A_eq, length(times)), tolerance = 1e-12)
  # RK4 at the integrator's step rule is good to ~2e-8 relative here; a fit to a
  # leaf experiment needs nothing like that, so 1e-7 is a check of the method, not
  # of the last digits.
  expect_equal(d$I, exact[, 1], tolerance = 1e-7)
  expect_equal(d$D, exact[, 2], tolerance = 1e-7)
  expect_equal(d$F, 1 - d$I - d$D)
})

test_that("acclimation relaxes exponentially, on the leaf temperature", {
  l <- set_inputs(thermal_leaf())
  temp <- 36
  sp <- with(atls_inputs, log1p(exp(acclim_s_ * (temp - acclim_T_))) / acclim_s_)
  A_eq <- with(atls_inputs, acclim_alpha_ / acclim_beta_ * sp)
  t <- c(0, 1, 5, 20)
  d <- leaf_damage_trajectory(l, t, temp, A0 = 0)
  expect_equal(d$A, A_eq * (1 - exp(-atls_inputs$acclim_beta_ * t)),
               tolerance = 1e-9)
})

test_that("the output grid does not coarsen the solution", {
  l <- set_inputs(thermal_leaf())
  fine <- seq(0, 3, by = 1 / 96)
  temp_fine <- 30 + 15 * pmax(0, sin(2 * pi * fine))  # daily heat peaks
  coarse_i <- seq(1, length(fine), by = 24)
  a <- leaf_damage_trajectory(l, fine, temp_fine, A0 = 0.5)
  b <- leaf_damage_trajectory(l, fine[coarse_i], temp_fine[coarse_i], A0 = 0.5)
  # The coarse run sees a coarser temperature series, so it agrees only where the
  # series agrees: compare against a fine run of that same piecewise-linear input.
  temp_interp <- approx(fine[coarse_i], temp_fine[coarse_i], xout = fine)$y
  c <- leaf_damage_trajectory(l, fine, temp_interp, A0 = 0.5)
  # Each run is within ~2e-7 of a 50x finer one here, the impairment rate rising
  # tenfold every ~8 K, so the two can differ by that much in either direction.
  expect_equal(b$I, c$I[coarse_i], tolerance = 1e-6)
  expect_equal(b$D, c$D[coarse_i], tolerance = 1e-6)
  expect_true(all(a$I >= 0 & a$D >= 0 & a$F <= 1))
  expect_true(all(diff(a$D) >= 0))                  # damage only accumulates
})

test_that("the trajectory refuses a malformed series", {
  l <- thermal_leaf()
  expect_error(l$damage_trajectory(c(0, 1, 1), c(30, 30, 30), 0, 0, 0),
               "strictly increasing")
  expect_error(l$damage_trajectory(c(0, 1), 30, 0, 0, 0), "same length")
  expect_identical(nrow(leaf_damage_trajectory(l, c(0, 1, 2), 30)), 3L)
})

# --- refusals -----------------------------------------------------------------

# --- gradients ------------------------------------------------------------------
#
# The trait gradient CARRIES the layer: both routes difference the solve's own
# dprofit and outputs, which read the block the layer writes. Its inputs that
# reach the solve are parameters on a leaf with the gate on, and refused by name
# on one with it off.

damage_gradient_pars <- c("A_", "I_", "D_", "topt_offset_", "dTopt_max_", "K_A_",
                          "k_i_", "k_i_baseline_", "k_rec_", "k_rec_Ea_",
                          "c_acclim_maint_", "c_repair_maint_",
                          "c_protect_maint_", "c_repair_flux_")

# The layer's solve-reaching inputs at non-default values. At 35 C rather than 25:
# the peaked Arrhenius is normalised at 25 C, so `topt_offset_` and `dTopt_max_`
# have an exactly zero effect there, and a difference of the solve is a clean
# referee at 35 C (see test_leaf.cpp's note on the 25 C floor).
gradient_leaf <- function(supply = leaf_supply_multilayer()) {
  set_inputs(thermal_leaf(supply), atls_inputs[c(damage_gradient_pars)])
}

test_that("the layer's inputs are gradient parameters only with its gate on", {
  off <- thermal_leaf(leaf_supply_multilayer())
  expect_error(leaf_gradient(psi_soil = c(1, 2, 3), soil_depth = c(.5, 1, 1.5),
                             x = off, traits = leaf_traits(), pars = "I_"),
               "input of the thermal-damage layer, which is off")
  on <- gradient_leaf()
  on$use_thermal_damage_ <- TRUE
  g <- leaf_gradient(psi_soil = c(1, 2, 3), soil_depth = c(.5, 1, 1.5),
                     PPFD = 1500, leaf_temp = 35, x = on, traits = leaf_traits())
  expect_true(all(damage_gradient_pars %in% rownames(g$gradient)))
  # Its DYNAMICS parameters reach no solved output and are not offered.
  expect_false(any(c("k_i_Ea_", "k_mat_", "acclim_beta_") %in%
                     rownames(g$gradient)))
  expect_error(leaf_gradient(psi_soil = 2, x = on, traits = leaf_traits(),
                             pars = "k_mat_"),
               "cannot differentiate")
})

test_that("with the gate on, the composite agrees with differencing the solve", {
  for (supply in list(leaf_supply_multilayer(), leaf_supply_singlelayer())) {
    single <- identical(supply$kind, "single")
    psi <- if (single) 1.0 else c(1, 2, 3)
    depth <- if (single) NULL else c(.5, 1, 1.5)
    l <- gradient_leaf(supply)
    l$use_thermal_damage_ <- TRUE
    pars <- c("vcmax_25", damage_gradient_pars)
    a <- leaf_gradient(psi_soil = psi, soil_depth = depth, PPFD = 1500,
                       leaf_temp = 35, x = l, traits = leaf_traits(), pars = pars,
                       method = "ift")
    b <- leaf_gradient(psi_soil = psi, soil_depth = depth, PPFD = 1500,
                       leaf_temp = 35, x = l, traits = leaf_traits(), pars = pars,
                       method = "fd", step = 1e-5)
    expect_identical(a$status, "interior", info = supply$kind)
    rel <- abs(a$gradient - b$gradient) / pmax(1, abs(a$gradient))
    expect_lt(max(rel), 1e-6, label = paste(supply$kind, "worst", max(rel)))
    # Impaired and damaged area do not photosynthesise: A falls with each.
    expect_lt(a$gradient["I_", "A"], 0)
    expect_lt(a$gradient["D_", "A"], 0)
    # The reused leaf is handed back with its inputs as they were.
    expect_identical(l$I_, atls_inputs$I_)
  }
})

test_that("with the gate on, the batch reproduces leaf_gradient() bit-for-bit", {
  pars <- c("vcmax_25", "I_", "A_", "topt_offset_", "c_repair_flux_")
  psis <- c(0.5, 1.5, 2.5)
  b <- leaf_batch(psi_soil = lapply(psis, function(p) p + c(0, 0.5, 1)),
                  soil_depth = list(c(.5, 1, 1.5)), PPFD = 1500, leaf_temp = 35)
  set_inputs(b$leaf, atls_inputs[damage_gradient_pars])
  b$leaf$use_thermal_damage_ <- TRUE
  for (method in c("ift", "fd")) {
    g <- leaf_gradient_batch(b, pars = pars, method = method)
    for (i in seq_along(psis)) {
      l <- gradient_leaf()
      l$use_thermal_damage_ <- TRUE
      r <- leaf_gradient(psi_soil = psis[i] + c(0, 0.5, 1),
                         soil_depth = c(.5, 1, 1.5), PPFD = 1500, leaf_temp = 35,
                         x = l, traits = leaf_traits(), pars = pars,
                         method = method)
      expect_identical(g$gradient[i, , ], r$gradient,
                       info = paste(method, "row", i))
      expect_identical(g$status[i], r$status, info = paste(method, "row", i))
    }
  }
  # Per-observation state through a theta matrix: each row its own I_.
  th <- matrix(rep(.gradient_theta_matrix(b, b$traits, TRUE), each = 1),
               nrow = 1)
  th <- th[rep(1, 3), , drop = FALSE]
  th[, match("I_", gradient_par_names())] <- c(0.05, 0.2, 0.4)
  g3 <- leaf_gradient_batch(b, pars = "vcmax_25", theta = th)
  expect_true(all(g3$status %in% c("interior", "pinned")))
  expect_false(isTRUE(all.equal(g3$value[1, "A"], g3$value[3, "A"])))
})

test_that("the trajectory gradient is the derivative of the trajectory", {
  l <- set_inputs(thermal_leaf())
  t <- seq(0, 3, by = 1 / 24)
  temp <- 30 + 14 * pmax(0, sin(2 * pi * t))
  g <- leaf_damage_gradient(l, t, temp, A0 = 0.5, I0 = 0.05, D0 = 0.01)
  expect_identical(dimnames(g$gradient)[[3]], l$damage_par_names())
  expect_identical(dim(g$gradient), c(length(t), 3L, 17L))
  base <- leaf_damage_trajectory(l, t, temp, A0 = 0.5, I0 = 0.05, D0 = 0.01)
  expect_identical(g$value, base)
  worst <- 0
  for (nm in l$damage_par_names()) {
    fd <- function(h) {
      shift <- function(sgn) {
        m <- set_inputs(thermal_leaf())
        init <- list(A0 = 0.5, I0 = 0.05, D0 = 0.01)
        if (nm %in% names(init)) {
          init[[nm]] <- init[[nm]] + sgn * h
        } else {
          m[[nm]] <- m[[nm]] + sgn * h
        }
        d <- leaf_damage_trajectory(m, t, temp, init$A0, init$I0, init$D0)
        cbind(d$A, d$I, d$D)
      }
      (shift(1) - shift(-1)) / (2 * h)
    }
    v <- if (nm %in% c("A0", "I0", "D0")) 0 else l[[nm]]
    h <- max(abs(v), 1) * 1e-5
    an <- g$gradient[, , nm]
    rel <- max(abs(an - fd(h)) / pmax(1e-3, abs(an)))
    worst <- max(worst, rel)
    expect_lt(rel, 1e-5, label = paste(nm, "worst relative", signif(rel, 3)))
  }
  # Pure-dynamics inputs genuinely move the states: a nonzero column each.
  expect_true(all(apply(abs(g$gradient), 3, max) > 0))
})

test_that("the trajectory gradient refuses an unknown input", {
  l <- thermal_leaf()
  expect_error(leaf_damage_gradient(l, c(0, 1), 30, pars = "vcmax_25"),
               "does not have")
  expect_error(l$damage_trajectory_gradient(c(0, 1), c(30, 30), 0, 0, 0, 99L),
               "out of range")
})

test_that("a temperature response with no optimum to move is refused", {
  # The activation and deactivation energies are settable, and the peaked
  # Arrhenius has an optimum only for Ea < H_d. With the gate off nothing
  # moves the optimum, so nothing needs one.
  l <- thermal_leaf()
  l$vcmax_ha_ <- 2.5e5
  expect_silent(drive_at(l, 30))
  l$use_thermal_damage_ <- TRUE
  expect_error(drive_at(l, 30), "exists only for 0 < Ea < H_d")
})

test_that("the closed form refuses a leaf with the layer on", {
  l <- thermal_leaf()
  drive_at(l, 25)
  l$set_model("TF24", "stem", "closed")
  l$use_thermal_damage_ <- TRUE
  expect_error(l$optimise(), "closed form is refused with use_thermal_damage_")
  # The exact method on the same leaf solves.
  l$set_model("TF24", "stem")
  l$optimise()
  expect_true(is.finite(l$profit_))
})
