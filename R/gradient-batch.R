# Trait gradients over many observations at once (issue #4, PLAN 11d stage 2).
#
# WHAT THIS FILE IS FOR. `leaf_gradient()` in `gradient.R` computes one
# observation's gradient in R, over primitives that are already C++. Measured
# after #69 and #70, a four-parameter gradient costs ~237 us per observation and
# the model work inside it is 6 us -- 1.5%. The other 98.5% is 112 crossings of
# the R boundary and the interpreter around them. `tests/cpp/bench_gradient.cpp`
# times the same composite in C++ at 1.8 us per trait.
#
# So the fix is not a faster gradient, it is ONE CROSSING. The composite is
# transcribed into `inst/include/phylloptim/gradient.hpp` and looped over
# observations there; this file is the surface that gets it a parameter matrix in
# and three arrays out. For a 1,327-observation MCMC at 30,000 draws that is
# ~2.6 h against ~7 min.
#
# ⚠️ WHAT STAYS IN R, AND WHY THAT IS NOT LAZINESS. The likelihood and the
# parameterisation Jacobian. The likelihood is the CALLER's model -- sigma,
# robustness, a hierarchy -- and putting one in C++ would commit this package to
# it. The Jacobian belongs in R because that is where the win comes from:
# `leaf-calibration` maps 40 fitted parameters onto 4 model ones, so C++ returns
# dY/dtheta for the four and R applies a 40 x 4 chain rule, vectorised over
# observations. That is exactly the `P_fit > P_model` structure
# `vignette("fitting")` identified as where an exact gradient wins at all -- see
# `gradient.R`'s header for the two cost laws and the retraction in each
# direction that produced them.
#
# ⚠️ AND THE THING TO BE MOST CAREFUL ABOUT: the C++ composite is a SECOND
# implementation of `leaf_gradient()`, so R remains the reference and the two are
# required to agree BIT-FOR-BIT. `gradient.hpp`'s header records the three rules
# that makes possible (literal arithmetic order, no fused multiply-add, explicit
# call sequencing) and `test-gradient-batch.R` is what holds them.

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
##' argument applies to the `Leaf`: constructing one from R costs ~204 µs, about
##' 70 solves, so a fit builds it once here rather than once per gradient (#52).
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
##' @return A `leaf_batch` object.
##' @seealso [leaf_gradient_batch()], [leaf_gradient()] for one observation.
##' @examples
##' b <- leaf_batch(psi_soil = c(1.0, 1.5, 2.0), PPFD = 900)
##' b
##' leaf_gradient_batch(b, pars = c("vcmax_25", "stem_b"))$gradient[, , "A"]
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
                       supply = leaf_supply_multilayer()) {
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

  # ⚠️ `.resolve_drivers()` PER ROW, and it is the only place the driver defaults
  # are applied. 1 m layers, the nominal networks, the single path's placeholder
  # depth -- all of them live in `set_drivers()`'s resolver, and reimplementing
  # any of them here or in C++ would create a second copy free to drift silently.
  # It runs once per observation for the whole fit, not once per gradient.
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
  # The two parameters that are DRIVERS: they are columns of `theta` rather than
  # part of the prepared bundle, because a fit moves them. Their base values are
  # read back out of the resolved drivers so that the default `theta` describes
  # the point the batch was built at -- the same thing `.gradient_theta()` does
  # for one observation.
  kmax <- pick("leaf_specific_conductance_max")
  resistance <- if (single) {
    vapply(resolved, function(a) a$root_network$r_R_V_sum[[1]], numeric(1))
  } else {
    rep(NA_real_, n)
  }

  # One row of `theta` when every observation shares these two, n rows when they
  # do not. Not an optimisation of the arithmetic -- the C++ side handles either
  # -- but of the matrix R builds per likelihood evaluation.
  shared <- n == 1L ||
    (all(kmax == kmax[[1L]]) &&
       (!single || all(resistance == resistance[[1L]])))

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
         kmax = if (shared) kmax[[1L]] else kmax,
         resistance = if (shared) resistance[[1L]] else resistance,
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
##' @section What it computes, and where the maths is written down:
##' The same five derivatives at the same solved operating point, by the same two
##' routes, with the same active-set test deciding between them. **Read
##' [leaf_gradient()]'s documentation for all of it** — the implicit function
##' theorem, why the second term is not a correction, why the premise is tested
##' rather than assumed, and what `status` means. Nothing here differs; this is a
##' transcription, and it is required to agree with the R version bit-for-bit.
##'
##' @section Speed, and which side of the boundary the figure belongs to:
##' A four-parameter gradient costs ~237 µs per observation through
##' [leaf_gradient()], of which the model work — two solves — is 6 µs, or 1.5%.
##' Everything else is dispatch and the R interpreter. The C++ composite runs at
##' ~1.8 µs per trait, so the same gradient is order 10 µs here.
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
##' there.
##'
##' @param batch a [leaf_batch()]
##' @param traits a [leaf_traits()] object: the point in trait space to
##'   differentiate at. Defaults to the one `batch` was built with. Ignored when
##'   `theta` is given.
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
##'     \item{`method`}{`"ift"` or `"fd"` per observation}
##'     \item{`status`}{`"interior"`, `"pinned"`, `"no-gradient"` or `"error"`}
##'     \item{`H`, `stationarity`}{the curvature and the implied Newton step
##'       `method` was decided on, per observation}
##'     \item{`message`}{`""`, or why that row failed}
##'   }
##'
##' @seealso [leaf_batch()], [leaf_gradient()] for one observation and for the
##'   derivation.
##' @examples
##' b <- leaf_batch(psi_soil = seq(0.5, 4, length.out = 6), PPFD = 900)
##' g <- leaf_gradient_batch(b, pars = c("vcmax_25", "stem_b"))
##' g$status
##' g$gradient[, "vcmax_25", "A"]
##' @export
leaf_gradient_batch <- function(batch,
                                traits = NULL,
                                pars = NULL,
                                theta = NULL,
                                step = 1e-6,
                                stationarity_tol = 1e-8,
                                method = c("auto", "ift", "fd"),
                                fast_stem_curve = TRUE) {
  if (!inherits(batch, "leaf_batch")) {
    stop("`batch` must come from leaf_batch(); got ", class(batch)[[1]],
         call. = FALSE)
  }
  method <- match.arg(method)
  if (!(is.numeric(step) && length(step) == 1L && step > 0)) {
    stop("`step` must be a single positive number", call. = FALSE)
  }
  # ⚠️ FIRST, before anything touches either handle. Both the drivers and the leaf
  # are external pointers, so a batch that has been through `serialize()` -- via
  # `saveRDS()`, or carried into a new session -- holds NULL ones, and the first
  # read of either reports "external pointer is not valid". That is true and gives
  # no hint that the fix is to rebuild, and it would come from whichever pointer
  # Rcpp unwrapped first rather than from anywhere that could explain itself. The
  # two are created together in `leaf_batch()`, so checking one covers both.
  gradient_batch_check(batch$drivers)

  # The supply path comes from the batch's record rather than from `batch$leaf`,
  # which would be a second pointer read before the check above could speak. It
  # was taken off the leaf when the batch was built, and the C++ side reads it off
  # the leaf again per call -- so the value used to validate `pars` and the value
  # the solve uses still have one source.
  single <- identical(batch$supply_kind, "single")
  par_names <- .gradient_par_names()
  if (is.null(pars)) {
    pars <- .gradient_available_pars(single)
  }
  .gradient_check_pars(pars, single)

  if (is.null(theta)) {
    if (is.null(traits)) {
      traits <- batch$traits
    } else if (!inherits(traits, "leaf_traits")) {
      stop("`traits` must come from leaf_traits()", call. = FALSE)
    }
    theta <- .gradient_theta_matrix(batch, traits)
  } else {
    if (!is.null(traits)) {
      stop("`traits` and `theta` both say where the gradient is taken; pass one ",
           "or the other. `theta` is the general form -- put the traits in its ",
           "columns.", call. = FALSE)
    }
    theta <- .gradient_check_theta(theta, batch, par_names)
  }

  res <- gradient_batch_run(batch$leaf, batch$drivers, theta,
                            match(pars, par_names) - 1L, step,
                            stationarity_tol, method, fast_stem_curve)

  dimnames(res$gradient) <- list(NULL, pars, .gradient_output_names())
  dimnames(res$value) <- list(NULL, .gradient_output_names())
  res
}

# --- internals ---------------------------------------------------------------

# `theta` from a leaf_traits() plus the batch's own conductance and resistance.
# The counterpart of `.gradient_theta()`, which builds the same vector for one
# observation -- and it reads the two non-traits back out of the resolved drivers
# for the same reason that one does: theta has to hold the value the solve will
# actually use.
#
# The two non-traits are `leaf_specific_conductance_max` and, on the
# single-potential path only, `resistance`: both are DRIVERS that a calibration
# fits like traits (#44), so the enumeration carries them after the traits.
.gradient_theta_matrix <- function(batch, traits) {
  nr <- batch$theta_nrow
  par_names <- .gradient_par_names()
  m <- matrix(NA_real_, nrow = nr, ncol = length(par_names))
  # By NAME, not by position. `leaf_traits()` is in `set_traits()`'s order and
  # test-gradient.R asserts it, but this is the one place a silent reordering
  # would produce plausible numbers for the wrong parameters.
  trait_names <- par_names[seq_len(length(par_names) - 2L)]
  tv <- unlist(traits)[trait_names]
  if (anyNA(tv)) {
    stop("`traits` is missing: ",
         paste(trait_names[is.na(tv)], collapse = ", "), call. = FALSE)
  }
  m[, seq_along(tv)] <- rep(tv, each = nr)
  m[, length(par_names) - 1L] <- batch$kmax
  # NA on the multi-layer path, where there is no such parameter and C++ never
  # reads the column.
  m[, length(par_names)] <- batch$resistance
  m
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
