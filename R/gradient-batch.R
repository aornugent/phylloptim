# Trait gradients over many observations in one crossing: leaf_batch() and
# leaf_gradient_batch().
#
# The composite is transcribed into inst/include/phylloptim/gradient.hpp and
# looped over observations there; this file gets a parameter matrix in and the
# arrays out. The likelihood and the parameterisation Jacobian stay with the
# caller, in R: a fit with P_fit > P_model applies its own chain rule.
#
# R/gradient.R is the reference, and the C++ composite must agree with it
# bit-for-bit; test-gradient-batch.R holds that, and R/gradient.R's header says
# what it constrains.

##' The observations a batch gradient is taken over
##'
##' Resolves N observations' drivers once and holds them C++-side, together with
##' one `Leaf`. Hand the result to [leaf_gradient_batch()], repeatedly: the
##' drivers do not change across a fit, only the parameters do.
##'
##' @section Why this is a separate call:
##' Because the conversion is expensive and the gradient is not. A [RootNetwork()]
##' costs 60-100 µs to cross the R boundary — some 60 times a trivial `.Call`,
##' because it is an `RcppR6` `list:` class and every crossing rebuilds a
##' five-element named list. Converting 1,327 of them per likelihood evaluation
##' would cost ~80 ms and swamp the ~13 ms of gradient it was carrying. The same
##' argument applies to the `Leaf`: constructing one is the largest single cost in
##' the R layer (the README's "Performance from R" has the measured ratio), so a fit
##' builds it once here rather than once per gradient (#52).
##'
##' ⚠️ **The result holds a C++ pointer, so it does not survive `saveRDS()` or a
##' new session.** Rebuild it rather than restoring it; [leaf_gradient_batch()]
##' says so rather than crashing if you try.
##'
##' @section Recycling:
##' The same rules as [leaf_solve()]. `PPFD`, `atm_vpd`, `ca`, `leaf_temp`,
##' `atm_o2_kpa`, `atm_kpa` and `leaf_specific_conductance_max` are recycled to a
##' common length. `psi_soil` is recycled when it is a plain numeric vector, each
##' element then being a single-layer observation; pass a list of numeric vectors
##' for multi-layer ones, and `soil_depth` follows the same rule. `root_network`
##' is either one [RootNetwork()] used for every observation, or a list of them.
##'
##' @inheritParams leaf_solve
##'
##' @param CF77_lambda Cowan-Farquhar's prescribed marginal value of water, in
##'   umol C (kg H2O)^-1, one per observation. Only `model = "CF77"` reads it,
##'   and [leaf_gradient_batch()] refuses it under any other model, as
##'   [leaf_solve()] does; `NA_real_` (the default) uses its aligned default
##'   (see [align_curves()]).
##' @param TF24_floor_lambda_o the `TF24_floor` curve's price of water as transpiration
##'   goes to zero, same units, one per observation. Only `model = "TF24_floor"`
##'   reads it, and it is refused under any other model; `NA_real_` (the
##'   default) leaves it unset, and [leaf_gradient_batch()] refuses that model
##'   without one, deliberately: at zero it is `TF24`.
##' @param LeastCost_beta least-cost theory's unit-cost ratio (Prentice et al.
##'   2014's beta): the cost of maintaining carboxylation capacity relative to the
##'   cost of transpiration, dimensionless, one per observation. Only
##'   `model = "LeastCost"` reads it, and it is refused under any other model;
##'   `NA_real_` (the default) uses its aligned
##'   default (see [align_curves()]). Zero is refused, since at
##'   zero the capacity term vanishes and what is left is water-use efficiency,
##'   which rises without bound as transpiration falls. The published value is 146.
##'
##' @return A `leaf_batch` object.
##' @seealso [leaf_gradient_batch()], [leaf_gradient()] for one observation.
##' @examples
##' b <- leaf_batch(psi_soil = c(1.0, 1.5, 2.0), PPFD = 900)
##' b
##' leaf_gradient_batch(b, pars = c("vcmax_25", "stem_P50"))$gradient[, , "A"]
##' @export
leaf_batch <- function(psi_soil,
                       PPFD = 900,
                       soil_depth = NULL,
                       root_network = NULL,
                       leaf_specific_conductance_max = 3.14e-5,
                       atm_vpd = 2.0,
                       ca = 40.0,
                       leaf_temp = 25.0,
                       atm_o2_kpa = 21.0,
                       atm_kpa = 101.3,
                       traits = leaf_traits(),
                       control = leaf_control(),
                       supply = leaf_supply_multilayer(),
                       CF77_lambda = NA_real_,
                       TF24_floor_lambda_o = NA_real_,
                       LeastCost_beta = NA_real_) {
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }

  # The same recycling `leaf_solve()` does, through the same helpers, because a
  # second copy of these rules would be free to drift from the first.
  layered <- .as_layer_list(psi_soil, "psi_soil")
  scalars <- list(PPFD = PPFD, atm_vpd = atm_vpd, ca = ca,
                  leaf_temp = leaf_temp, atm_o2_kpa = atm_o2_kpa,
                  atm_kpa = atm_kpa,
                  leaf_specific_conductance_max = leaf_specific_conductance_max)
  n <- max(length(layered), vapply(scalars, length, integer(1)))
  if (n == 0L) {
    stop("nothing to differentiate: every driver has length zero", call. = FALSE)
  }
  layered <- .recycle_to(layered, n, "psi_soil")
  nms <- names(scalars)
  scalars <- lapply(seq_along(scalars), function(i) {
    .recycle_to(scalars[[i]], n, nms[[i]])
  })
  names(scalars) <- nms

  depths <- if (is.null(soil_depth)) NULL else
    .recycle_to(.as_layer_list(soil_depth, "soil_depth"), n, "soil_depth")
  # inherits() first: a RootNetwork *is* a list, so length() on one would be 5
  # and it would be read as five observations' worth of networks.
  networks <- if (is.null(root_network)) {
    NULL
  } else if (inherits(root_network, "RootNetwork")) {
    rep(list(root_network), n)
  } else {
    .recycle_to(root_network, n, "root_network")
  }

  l <- leaf_model(traits, control, supply)

  # `.resolve_drivers()` per row: the only place the driver defaults are applied,
  # so neither this file nor C++ keeps a second copy that could drift.
  resolved <- lapply(seq_len(n), function(i) {
    .resolve_drivers(l, layered[[i]], scalars$PPFD[[i]],
                     if (is.null(depths)) NULL else depths[[i]],
                     if (is.null(networks)) NULL else networks[[i]],
                     scalars$leaf_specific_conductance_max[[i]],
                     scalars$atm_vpd[[i]], scalars$ca[[i]],
                     scalars$leaf_temp[[i]], scalars$atm_o2_kpa[[i]],
                     scalars$atm_kpa[[i]])
  })

  single <- identical(l$supply_kind, "single")
  pick <- function(field) {
    vapply(resolved, function(a) a[[field]], numeric(1))
  }
  # The two parameters that are drivers are columns of `theta`, not part of the
  # prepared bundle, because a fit moves them; their base values come from the
  # resolved drivers, as in `.gradient_theta()`.
  kmax <- pick("leaf_specific_conductance_max")
  resistance <- if (single) {
    vapply(resolved, function(a) a$root_network$r_R_V_sum[[1]], numeric(1))
  } else {
    rep(NA_real_, n)
  }

  prices <- list(CF77_lambda = CF77_lambda,
                 TF24_floor_lambda_o = TF24_floor_lambda_o,
                 LeastCost_beta = LeastCost_beta)
  prices <- lapply(names(prices), function(nm) {
    .recycle_to(as.numeric(prices[[nm]]), n, nm)
  })
  names(prices) <- c("CF77_lambda", "TF24_floor_lambda_o", "LeastCost_beta")

  # One row of `theta` when every observation shares these two and the prices,
  # n rows when they do not. Not an optimisation of the arithmetic -- the C++
  # side handles either -- but of the matrix R builds per likelihood evaluation.
  shared <- n == 1L ||
    (all(kmax == kmax[[1L]]) &&
       (!single || all(resistance == resistance[[1L]])) &&
       all(vapply(prices, function(x) length(unique(x)) == 1L, logical(1))))
  one_or_all <- function(x) if (shared) x[[1L]] else x

  structure(
    list(leaf = l,
         drivers = gradient_batch_prepare(
           root_network = lapply(resolved, `[[`, "root_network"),
           psi_soil = lapply(resolved, `[[`, "psi_soil"),
           soil_depth = lapply(resolved, `[[`, "soil_depth"),
           PPFD = pick("PPFD"), atm_vpd = pick("atm_vpd"), ca = pick("ca"),
           leaf_temp = pick("leaf_temp"), atm_o2_kpa = pick("atm_o2_kpa"),
           atm_kpa = pick("atm_kpa")),
         n = n,
         supply_kind = l$supply_kind,
         traits = traits,
         kmax = one_or_all(kmax),
         resistance = one_or_all(resistance),
         CF77_lambda = one_or_all(prices$CF77_lambda),
         TF24_floor_lambda_o = one_or_all(prices$TF24_floor_lambda_o),
         LeastCost_beta = one_or_all(prices$LeastCost_beta),
         theta_nrow = if (shared) 1L else n),
    class = "leaf_batch")
}

##' @export
print.leaf_batch <- function(x, ...) {
  cat(sprintf("<leaf_batch> %d observation%s, %s supply\n",
              x$n, if (x$n == 1L) "" else "s", x$supply_kind))
  cat(sprintf("  differentiable: %s\n",
              paste(.gradient_available_pars(
                identical(x$supply_kind, "single")), collapse = ", ")))
  invisible(x)
}

##' Trait gradients over a batch of observations
##'
##' [leaf_gradient()], composed in C++ and vectorised over the observations of a
##' [leaf_batch()]. One crossing of the R boundary per call, in place of 112 per
##' observation.
##'
##' @param psi one collar water potential per observation to evaluate at
##'   (positive magnitudes, MPa) instead of solving for the profit-maximising
##'   one. See [leaf_gradient()]'s section on it. Recycling a single value across
##'   the batch is deliberately not offered — rows sharing an imposed collar is a
##'   coincidence rather than a design, and broadcasting one would make a length
##'   mismatch look intentional.
##' @param dpsi_dtheta how each prescribed `psi` responds to each parameter:
##'   an `n` × `length(pars)` matrix, or one vector of `length(pars)` shared by
##'   every observation. Defaults to zero. Must be finite, for the reason
##'   [leaf_gradient()] records. Only meaningful with `psi`.
##'
##' @section What it computes, and where the maths is written down:
##' The same five derivatives at the same solved operating point, by the same two
##' routes, with the same active-set test deciding between them. **Read
##' [leaf_gradient()]'s documentation for all of it** — the implicit function
##' theorem, why the second term is not a correction, why the premise is tested
##' rather than assumed, and what `status` means. Nothing here differs; this is a
##' transcription, and it is required to agree with the R version bit-for-bit.
##'
##' @section Speed, and which side of the boundary the figure belongs to:
##' Through [leaf_gradient()], almost all of a four-parameter gradient's cost per
##' observation is dispatch and the R interpreter; the model work — two solves — is
##' a few percent of it. Composed here in C++, the same gradient costs over an order
##' of magnitude less per observation. The README's "Performance from R" has the
##' measured ratios.
##'
##' ⚠️ **That is a per-observation figure and it needs a batch to be realised.**
##' Calling this with one observation pays the whole R-side call overhead for one
##' row's worth of work and will not be 20× anything. The number of observations
##' is the lever.
##'
##' @section Per-row status, not an error:
##' A proposal during a fit **will** reach operating points the solve cannot
##' handle. That has to cost those rows rather than the whole dataset, so every
##' row reports its own `status` and a failure is `"error"` with a `message`
##' rather than a condition. A failed row's gradient is all `NA`, not partially
##' filled: a parameter loop that failed halfway has valid entries before the
##' failure, and returning them next to `status == "error"` would invite a caller
##' who checks the status per row to read a gradient with a hole in it. `value`,
##' `H` and `stationarity` are kept where they were determined, because they
##' describe the point rather than the derivative.
##'
##' @section Supplying `theta` directly:
##' By default the parameters are `traits` — the batch's, or the ones passed here
##' — with `leaf_specific_conductance_max` and `resistance` taken from the
##' batch's own drivers. Pass `theta` instead when they vary by observation, which
##' is what a hierarchical fit needs: a matrix with one row per observation (or one
##' row shared by all of them) and one column per [gradient_par_names()] entry, in
##' that order. `resistance` is ignored on the multi-layer path and may be `NA`
##' there. An `NA` in the column of the level parameter `model` owns
##' (`CF77_lambda_` or `LeastCost_beta`) takes its aligned default, as it does
##' when `theta` is built for you.
##'
##' @param batch a [leaf_batch()]
##' @param traits a [leaf_traits()] object: the point in trait space to
##'   differentiate at. Defaults to the one `batch` was built with. Ignored when
##'   `theta` is given.
##' @param model which optimality model to differentiate, as in
##'   [leaf_gradient()]: `"collar"` by default, or any [cost_curve_names()] entry.
##' @param pars what to differentiate with respect to, as in [leaf_gradient()].
##'   ⚠️ **Always pass it.** It is `P_model`, so the default — all of them — is
##'   the most expensive thing you can ask for.
##' @param theta the model parameters as a matrix, for the case where they vary
##'   by observation. See the section above.
##' @inheritParams leaf_gradient
##'
##' @return A list with
##'   \describe{
##'     \item{`gradient`}{an `n` × `length(pars)` × 5 array, the last dimension
##'       being `A`, `gc`, `psi_stem`, `collar` and `profit`}
##'     \item{`value`}{an `n` × 5 matrix of the solved outputs}
##'     \item{`method`}{`"ift"`, `"fd"` or `"prescribed"` per observation}
##'     \item{`status`}{`"interior"`, `"pinned"`, `"no-gradient"`, `"error"`,
##'       or with `psi` one of `"prescribed"`, `"refused"`, `"no-gradient"`}
##'     \item{`H`, `stationarity`}{the curvature and the implied Newton step
##'       `method` was decided on, per observation}
##'     \item{`M`}{an `n` × `length(pars)` matrix of mixed partials
##'       `d2profit/dpsi dtheta`}
##'     \item{`dY_dpsi`}{an `n` × 5 matrix of `dY/dpsi` at fixed traits}
##'     \item{`psi`}{the collar each row was evaluated at}
##'     \item{`message`}{`""`, or why that row failed}
##'   }
##'
##' @seealso [leaf_batch()], [leaf_gradient()] for one observation and for the
##'   derivation.
##' @examples
##' b <- leaf_batch(psi_soil = seq(0.5, 4, length.out = 6), PPFD = 900)
##' g <- leaf_gradient_batch(b, pars = c("vcmax_25", "stem_P50"))
##' g$status
##' g$gradient[, "vcmax_25", "A"]
##' @export
leaf_gradient_batch <- function(batch,
                                traits = NULL,
                                pars = NULL,
                                theta = NULL,
                                psi = NULL,
                                dpsi_dtheta = NULL,
                                step = 1e-6,
                                stationarity_tol = 1e-8,
                                method = c("auto", "ift", "fd"),
                                fast_stem_curve = TRUE,
                                model = "collar") {
  if (!inherits(batch, "leaf_batch")) {
    stop("`batch` must come from leaf_batch(); got ", class(batch)[[1]],
         call. = FALSE)
  }
  method <- match.arg(method)
  if (!(is.numeric(step) && length(step) == 1L && step > 0)) {
    stop("`step` must be a single positive number", call. = FALSE)
  }
  # First, before either handle is read: both are external pointers, NULL after
  # serialize(), and the bare Rcpp error gives no hint to rebuild. They are
  # created together, so checking one covers both.
  gradient_batch_check(batch$drivers)

  # The supply path from the batch's record, not `batch$leaf`, which would be a
  # pointer read before the check above. It was taken off the leaf at build time.
  single <- identical(batch$supply_kind, "single")
  # The thermal-damage inputs are parameters only with the batch leaf's gate on.
  damage <- isTRUE(batch$leaf$use_thermal_damage_)
  par_names <- .gradient_par_names()
  if (is.null(pars)) {
    pars <- .gradient_available_pars(single, model, damage)
  }
  .gradient_check_pars(pars, single, model, damage)

  # The same route object leaf_gradient() uses, so the two cannot disagree about
  # what a model name selects, its `fd_step` or the supply path it needs. Built
  # for those; it seats the batch's leaf too, which C++ re-seats for every row
  # (`gradient::route_seat`).
  route <- .gradient_route(batch$leaf, model, list(kind = batch$supply_kind))
  # The curve as an index into cost_curve_names(), -1 for the collar route:
  # C++'s `gradient::Settings::curve`, never seen by an R caller.
  curve <- if (identical(model, .GRADIENT_COLLAR_ROUTE)) {
    -1L
  } else {
    as.integer(match(model, cost_curve_names()) - 1L)
  }
  # The finite-difference arm's step: the caller's if named, as in
  # leaf_gradient(), otherwise the route's, which clears the solver's resolution.
  fd_step <- if (missing(step)) route$fd_step else step
  # ProfitMax's normaliser is cleared at every perturbation, so C++ re-solves it,
  # as `.gradient_route()` does in R. A boolean passed as a double.
  reseat_profitmax <- if (identical(model, "ProfitMax")) 1 else 0

  if (is.null(theta)) {
    if (is.null(traits)) {
      traits <- batch$traits
    } else if (!inherits(traits, "leaf_traits")) {
      stop("`traits` must come from leaf_traits()", call. = FALSE)
    }
    theta <- .gradient_theta_matrix(batch, traits, damage, model)
  } else {
    if (!is.null(traits)) {
      stop("`traits` and `theta` both say where the gradient is taken; pass one ",
           "or the other. `theta` is the general form -- put the traits in its ",
           "columns.", call. = FALSE)
    }
    theta <- .gradient_fill_level_defaults(
      .gradient_check_theta(theta, batch, par_names), model, par_names)
  }

  .gradient_check_theta_prices(theta, model, par_names)

  psi <- .gradient_check_psi_batch(psi, batch$n, method)
  dpsi_dtheta <- .gradient_dpsi_dtheta_batch(dpsi_dtheta, pars, batch$n,
                                             !is.null(psi))

  # One crossing, one condition wrapper: per-row failures come back as `status`,
  # so this catches only what takes the whole call down.
  res <- with_phylloptim_conditions(gradient_batch_run(batch$leaf, batch$drivers, theta,
                            match(pars, par_names) - 1L, step,
                            stationarity_tol, method, fast_stem_curve,
                            psi, dpsi_dtheta, curve, fd_step,
                            reseat_profitmax))

  dimnames(res$gradient) <- list(NULL, pars, .gradient_output_names())
  dimnames(res$value) <- list(NULL, .gradient_output_names())
  dimnames(res$M) <- list(NULL, pars)
  dimnames(res$dY_dpsi) <- list(NULL, .gradient_output_names())
  res
}

# --- internals ---------------------------------------------------------------

# The batch counterparts of `.gradient_check_psi()` and `.gradient_dpsi_dtheta()`.
# A length-1 `psi` is not recycled: rows sharing an imposed collar is a
# coincidence, and broadcasting one would make a length mismatch look intended.
.gradient_check_psi_batch <- function(psi, n, method) {
  if (is.null(psi)) {
    return(NULL)
  }
  if (!identical(method, "auto")) {
    stop("`method` cannot be given with `psi`; see ?leaf_gradient.",
         call. = FALSE)
  }
  if (!(is.numeric(psi) && length(psi) == n)) {
    stop("`psi` must be one collar potential per observation (", n,
         " for this batch); got ", length(psi), call. = FALSE)
  }
  if (!all(is.finite(psi) & psi > 0)) {
    stop("`psi` must be finite and positive: it is a collar water potential ",
         "in MPa, as a positive magnitude.", call. = FALSE)
  }
  as.numeric(psi)
}

.gradient_dpsi_dtheta_batch <- function(dpsi_dtheta, pars, n, prescribed) {
  if (is.null(dpsi_dtheta)) {
    # NULL is zero on the C++ side, so nothing is allocated for the common case
    # of a partial derivative at fixed collar.
    return(NULL)
  }
  if (!prescribed) {
    stop("`dpsi_dtheta` is the trait response of a collar potential YOU ",
         "imposed, so it needs `psi`.", call. = FALSE)
  }
  # Finite, for the reason `.gradient_dpsi_dtheta()` gives; checked here too
  # because the matrix form never reaches it.
  if (!is.numeric(dpsi_dtheta) || !all(is.finite(dpsi_dtheta))) {
    stop("`dpsi_dtheta` must be numeric and finite -- no NA, NaN or Inf. It ",
         "is dpsi/dtheta for a collar YOU imposed; a non-finite value returns ",
         "an all-infinite gradient row with a plausible `profit` in it.",
         call. = FALSE)
  }
  if (is.null(dim(dpsi_dtheta))) {
    if (length(dpsi_dtheta) != length(pars)) {
      stop("`dpsi_dtheta` as a vector is one value per parameter (",
           length(pars), "), shared by every observation; got ",
           length(dpsi_dtheta), ". Pass an ", n, " x ", length(pars),
           " matrix to vary it by observation.", call. = FALSE)
    }
    if (!is.null(names(dpsi_dtheta))) {
      dpsi_dtheta <- .gradient_dpsi_dtheta(dpsi_dtheta, pars, TRUE)
    }
    return(matrix(dpsi_dtheta, nrow = n, ncol = length(pars), byrow = TRUE))
  }
  if (!identical(dim(dpsi_dtheta), c(n, length(pars)))) {
    stop("`dpsi_dtheta` must be ", n, " x ", length(pars),
         " (observations x `pars`); got ",
         paste(dim(dpsi_dtheta), collapse = " x "), call. = FALSE)
  }
  if (!is.null(colnames(dpsi_dtheta))) {
    if (!setequal(colnames(dpsi_dtheta), pars)) {
      stop("`dpsi_dtheta` has column names, so they must be exactly `pars`.",
           call. = FALSE)
    }
    dpsi_dtheta <- dpsi_dtheta[, pars, drop = FALSE]
  }
  matrix(as.numeric(dpsi_dtheta), nrow = n, ncol = length(pars))
}

# `theta` from a leaf_traits() plus the batch's own conductance and resistance:
# the counterpart of `.gradient_theta()`, reading the non-traits back out of the
# resolved drivers so theta holds the value the solve uses.
.gradient_theta_matrix <- function(batch, traits, damage = FALSE,
                                   model = NULL) {
  nr <- batch$theta_nrow
  par_names <- .gradient_par_names()
  m <- matrix(NA_real_, nrow = nr, ncol = length(par_names),
              dimnames = list(NULL, par_names))
  # Every column by name (hazard 13): a positional tail column shifts silently
  # when a parameter is appended.
  trait_names <- setdiff(par_names, .gradient_non_traits)
  tv <- unlist(traits)[trait_names]
  if (anyNA(tv)) {
    stop("`traits` is missing: ",
         paste(trait_names[is.na(tv)], collapse = ", "), call. = FALSE)
  }
  m[, trait_names] <- rep(tv, each = nr)
  m[, "leaf_specific_conductance_max"] <- batch$kmax
  # NA on the multi-layer path, where there is no such parameter and C++ never
  # reads the column.
  m[, "resistance"] <- batch$resistance
  # NA unless the batch was built for the route that owns it; C++ reads every
  # slot, so each holds whatever the solve will use.
  m[, "CF77_lambda_"] <- batch$CF77_lambda
  m[, "TF24_floor_lambda_o"] <- batch$TF24_floor_lambda_o
  m[, "LeastCost_beta"] <- batch$LeastCost_beta
  m <- .gradient_fill_level_defaults(m, model, par_names)
  # The thermal-damage inputs, off the batch's leaf, when its gate is on; NA (and
  # unread) otherwise. One row for every observation: pass `theta` to give each
  # observation its own state.
  if (damage) {
    dv <- .gradient_damage_values(batch$leaf)
    m[, names(dv)] <- rep(dv, each = nr)
  }
  # Handed over WITHOUT names, as `leaf_gradient_batch()` documents: C++ indexes by
  # position, and the check there rejects any colnames that are not this order.
  unname(m)
}

.gradient_check_theta <- function(theta, batch, par_names) {
  if (!is.matrix(theta) || !is.numeric(theta)) {
    stop("`theta` must be a numeric matrix, one column per parameter",
         call. = FALSE)
  }
  if (ncol(theta) != length(par_names)) {
    stop("`theta` must have ", length(par_names), " columns, one per ",
         "gradient_par_names() entry; got ", ncol(theta), call. = FALSE)
  }
  if (!is.null(colnames(theta)) && !identical(colnames(theta), par_names)) {
    stop("`theta`'s columns must be gradient_par_names(), in that order. The ",
         "C++ side indexes them by POSITION, so a reordering would ",
         "differentiate the wrong parameters rather than fail.", call. = FALSE)
  }
  if (!nrow(theta) %in% c(1L, batch$n)) {
    stop("`theta` must have 1 row or one per observation (", batch$n, "); got ",
         nrow(theta), call. = FALSE)
  }
  theta
}

# An unset slot on the route that owns it takes the leaf's own default, whether
# the matrix was built here or supplied, so `theta = NA` means what it means in
# leaf_gradient(). Columns by position: a supplied matrix may have no names.
.gradient_fill_level_defaults <- function(m, model, par_names) {
  if (is.null(model)) {
    return(m)
  }
  for (field in names(.level_defaults)) {
    if (.owns_field(model, field)) {
      j <- match(field, par_names)
      m[is.na(m[, j]), j] <- .level_defaults[[field]]
    }
  }
  m
}

# The owned columns, checked as leaf_solve() checks its arguments: a price the
# model does not read is refused, and the owned field of a model with no default
# must be set on every row. Left to C++, a missing one is a per-row error, which
# a fit reads as observations that failed to solve.
.gradient_check_theta_prices <- function(m, model, par_names) {
  col <- function(field) m[, match(field, par_names)]
  .check_prices(model, col("CF77_lambda_"), col("TF24_floor_lambda_o"),
                col("LeastCost_beta"))
  if (!model %in% names(.solve_price_fields)) {
    return(invisible(m))
  }
  field <- .solve_price_fields[[model]]
  if (anyNA(m[, match(field, par_names)])) {
    stop("the ", model, " cost curve needs ", field, " set on every ",
         "observation: pass `", field, " =` to leaf_batch(), or fill that ",
         "column of `theta`; it has no default.", call. = FALSE)
  }
  invisible(m)
}
