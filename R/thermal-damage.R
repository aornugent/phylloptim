##' Integrate the leaf's thermal-damage states over a temperature series
##'
##' The leaf's lasting damage model, run on its own: three fractions of leaf
##' area, functional ⇌ impaired → damaged, plus an acclimation state, driven by a
##' leaf-temperature time series. This is the route for fitting the layer's
##' parameters to a leaf experiment without a whole-plant model; the comparison
##' with measurements (Fv/Fm, say) belongs to the caller.
##'
##' The dynamics, per day, with \eqn{T} the leaf temperature and \eqn{F = 1 - I - D}:
##'
##' \preformatted{
##'   dA/dt    = acclim_alpha * softplus(T - acclim_T; acclim_s) - acclim_beta * A
##'   dI/dt    = k_imp(T) * F - (k_rec(T) + k_mat(T)) * I
##'   dD/dt    = k_mat(T) * I
##'   k_imp(T) = k_i   * exp[(k_i_Ea/R)(1/T_iso - 1/T)]
##'   k_rec(T) = k_rec * exp[(k_rec_Ea/R)(1/T_25 - 1/T)]
##'   k_mat(T) = k_mat * exp[(k_mat_Ea/R)(1/T_25 - 1/T)]
##' }
##'
##' with \eqn{T_{iso}} = `k_i_Tiso_` shifted by `topt_offset_` and the
##' acclimation, and \eqn{T_{25}} = 25 C. The impairment rate is the first-order
##' PSII heat inactivation of Bison & Michaletz (2026, Nature Ecology & Evolution,
##' \doi{10.1038/s41559-026-03196-x}), and its defaults their 177-species fit:
##' the rate common to every species at the isokinetic temperature (276 per day,
##' at 50.5 C) and the mean activation energy (316 kJ/mol), so it rises tenfold
##' about every 6 K. Nothing switches the rates on or off.
##'
##' Every parameter is a field of `x` (`$k_i_`, `$k_i_Tiso_`, `$acclim_beta_`,
##' ...); set them on the leaf before calling. Neither the gate
##' `use_thermal_damage_` nor the drivers are read. The rates themselves are
##' `x$damage_rates(leaf_temp, A, I, D)`.
##'
##' Temperature is linear between the given times. The integrator is fixed-step
##' RK4 with each interval cut so the fastest rate times the step stays below 0.1
##' and the temperature moves at most 0.25 K per step, so a coarse output grid does
##' not coarsen the solution.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param times days, strictly increasing
##' @param leaf_temp leaf temperature at each time, deg C
##' @param A0,I0,D0 the acclimation state and the impaired and damaged fractions
##'   at `times[1]`
##' @return A data.frame with `time`, `leaf_temp`, `A`, `I`, `D` and the
##'   functional fraction `F = 1 - I - D`.
##' @examples
##' l <- leaf_model()
##' # a 6-hour 46 C heatwave, then two days at 25 C
##' t <- seq(0, 2.25, by = 1 / 24)
##' temp <- ifelse(t < 0.25, 46, 25)
##' d <- leaf_damage_trajectory(l, t, temp)
##' tail(d, 3)
##' @export
leaf_damage_trajectory <- function(x, times, leaf_temp, A0 = 0, I0 = 0, D0 = 0) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  times <- as.numeric(times)
  leaf_temp <- as.numeric(leaf_temp)
  if (length(leaf_temp) == 1L) {
    leaf_temp <- rep(leaf_temp, length(times))
  }
  y <- x$damage_trajectory(times, leaf_temp, A0, I0, D0)
  n <- length(times)
  out <- data.frame(time = times, leaf_temp = leaf_temp,
                    A = y[seq_len(n)], I = y[n + seq_len(n)],
                    D = y[2L * n + seq_len(n)])
  out$F <- 1 - out$I - out$D
  out
}

##' Gradient of the leaf's thermal-damage trajectory
##'
##' The exact sensitivities of [leaf_damage_trajectory()]'s states to the layer's
##' dynamics parameters and initial states: what a fit of the layer to a leaf
##' experiment needs. Computed by forward-mode differentiation of the same RK4 the
##' trajectory runs -- one pass per input -- so it is the derivative of the
##' numbers [leaf_damage_trajectory()] returns, not of an idealised solution.
##'
##' The inputs, named as the leaf's fields: the impairment rate (`k_i_`,
##' `k_i_Ea_`, `k_i_Tiso_`), repair (`k_rec_`, `k_rec_Ea_`), maturation (`k_mat_`,
##' `k_mat_Ea_`), acclimation (`acclim_alpha_`, `acclim_beta_`, `acclim_T_`,
##' `acclim_s_`), the acclimated shift of the isokinetic temperature
##' (`dTimpair_max_`, `K_A_`) and the constitutive offset (`topt_offset_`), then
##' the initial states `A0`, `I0` and `D0`.
##'
##' ⚠️ The sub-step count is chosen from the parameter VALUES, so it is piecewise
##' constant in them: the gradient is exact for the discrete map at that count,
##' and a finite difference straddling a change of count differs from it by the
##' integrator's own error (~1e-8 relative).
##'
##' @inheritParams leaf_damage_trajectory
##' @param pars which inputs to differentiate with respect to; all of them by
##'   default.
##' @return A list: `value`, the [leaf_damage_trajectory()] data.frame, and
##'   `gradient`, an array `[time, state, input]` with states `A`, `I`, `D`.
##' @examples
##' l <- leaf_model()
##' t <- seq(0, 2, by = 1 / 24)
##' g <- leaf_damage_gradient(l, t, ifelse(t < 0.25, 46, 25),
##'                           pars = c("k_i_Tiso_", "k_rec_"))
##' g$gradient[length(t), , ]
##' @export
leaf_damage_gradient <- function(x, times, leaf_temp, A0 = 0, I0 = 0, D0 = 0,
                                 pars = NULL) {
  value <- leaf_damage_trajectory(x, times, leaf_temp, A0, I0, D0)
  names_all <- x$damage_par_names()
  if (is.null(pars)) {
    pars <- names_all
  }
  unknown <- setdiff(pars, names_all)
  if (length(unknown)) {
    stop("`pars` names inputs the damage trajectory does not have: ",
         paste(unknown, collapse = ", "), ". Available: ",
         paste(names_all, collapse = ", "), call. = FALSE)
  }
  n <- length(value$time)
  g <- x$damage_trajectory_gradient(value$time, value$leaf_temp, A0, I0, D0,
                                    match(pars, names_all) - 1L)
  gradient <- array(g, dim = c(n, 3L, length(pars)),
                    dimnames = list(NULL, c("A", "I", "D"), pars))
  list(value = value, gradient = gradient)
}
