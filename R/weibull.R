# Converting published vulnerability curves to (P50, c).
#
# The model is parameterised on (P50, c) since the curve became
# `f(psi) = 2^-(psi/P50)^c`. Published curves arrive as almost anything else: a
# P88, a P95, a P50 plus a slope. All of it inverts in closed form, so the
# conversion belongs here rather than in each study that needs it -- which is what
# four hand-written copies of `psi_at_plc` downstream had become.
#
# ⚠️ EVERY POTENTIAL HERE IS A POSITIVE MAGNITUDE. Published P50 is negative
# (-3.95 MPa); this package uses magnitudes throughout, so `P50` means `|P50|`.
# A negative input is refused rather than silently negated, because a caller who
# passes a signed value is working in the other convention and everything else they
# pass is likely signed too.

# log2(1/f) for a remaining-conductivity fraction f. This is the only place the
# quantile arithmetic lives.
.weibull_L <- function(f) log(1 / f) / log(2)

# The three named quantiles, as remaining fractions rather than PLC.
.WEIBULL_F <- c(P50 = 0.50, P88 = 0.12, P95 = 0.05, P99 = 0.01)

# S50 = 50*c*ln2/|P50|. The 50 is because the slope is conventionally quoted in
# % PLC per MPa, so it carries the percentage scale.
.WEIBULL_S50_K <- 50 * log(2)

##' The water potential at a given loss of conductivity
##'
##' @description
##' Inverts `f(psi) = 2^-(psi/P50)^c` for the potential at which a given fraction
##' of conductivity has been lost.
##'
##' @param P50 potential at 50% loss, MPa as a positive magnitude
##' @param c shape parameter of the Weibull curve (unitless)
##' @param plc proportion of conductivity LOST, in `(0, 1)`. `plc = 0.5` returns
##'   `P50` by construction; `0.88`, `0.95` and `0.99` return P88, P95 and P99.
##'
##' @return The potential, MPa as a positive magnitude.
##' @seealso [weibull_p50_c()] for the reverse, and [leaf_traits()] for where the
##'   result goes.
##' @examples
##' psi_at_plc(3.4, 2.680147, 0.5)   # returns P50
##' psi_at_plc(3.4, 2.680147, 0.95)  # P95, which is also `psi_crit`
##' @export
psi_at_plc <- function(P50, c, plc) {
  .weibull_check_p50_c(P50, c)
  if (!(is.numeric(plc) && all(plc > 0) && all(plc < 1))) {
    stop("`plc` is a proportion of conductivity LOST and must be in (0, 1); ",
         "0.95 rather than 95.", call. = FALSE)
  }
  P50 * .weibull_L(1 - plc)^(1 / c)
}

##' The slope of the vulnerability curve at P50
##'
##' @description
##' `S50 = 50*c*ln(2)/|P50|`, in % PLC per MPa, which is the other form
##' vulnerability curves are commonly published in.
##'
##' ⚠️ **`S50` is not an independent shape parameter** — it has `P50` built into it
##' by construction, so a curve is not determined by `S50` alone. It is exported
##' rather than kept internal for exactly that reason: a caller reading a published
##' `S50` needs to see what it depends on.
##'
##' @inheritParams psi_at_plc
##' @return The slope at P50, % PLC per MPa.
##' @seealso [weibull_p50_c()], which accepts `S50` as one of its two inputs.
##' @examples
##' weibull_s50(3.4, 2.680147)
##' @export
weibull_s50 <- function(P50, c) {
  .weibull_check_p50_c(P50, c)
  .WEIBULL_S50_K * c / P50
}

##' Convert a Weibull scale parameter to P50
##'
##' @description
##' `b` is the scale of the conventional form `exp(-(psi/b)^c)`; this package is
##' parameterised on `P50` instead. The two are one substitution apart,
##' `b = P50 / (ln 2)^(1/c)`, so
##'
##' \deqn{P_{50} = b\,(\ln 2)^{1/c}}
##'
##' Exported because published fits and older datasets carry `b`, and the
##' alternative is that conversion being pasted at each call site. `Leaf` reports
##' `stem_b` and `root_b` as derived read-only accessors, which is the forward
##' direction of the same identity.
##'
##' ⚠️ **`P50` is always SMALLER than `b`**, by a factor of `(ln 2)^(1/c)` — 0.87
##' at `c = 2.68`. So passing a `b` where a `P50` is wanted, or the reverse, gives
##' a plausible number roughly 15% out rather than an error. Neither function can
##' detect it: both arguments are positive potentials in MPa.
##'
##' @param b scale parameter of the Weibull curve, MPa as a positive magnitude
##' @inheritParams psi_at_plc
##' @return `P50`, MPa as a positive magnitude.
##' @seealso [weibull_p50_c()], [psi_at_plc()]
##' @examples
##' weibull_p50_from_b(3.898245, 2.680147)   # the package default, 3.4
##' @export
weibull_p50_from_b <- function(b, c) {
  .weibull_check_p50_c(b, c)
  b * log(2)^(1 / c)
}

##' Recover (P50, c) from any two published quantities
##'
##' @description
##' A vulnerability curve has two parameters, so **any two** of the quantities
##' below determine it. This solves for `(P50, c)` from whichever pair you have.
##'
##' @section Which pairs are closed-form:
##' All of them except `S50` paired with a quantile other than P50, which is
##' transcendental in `c` and is solved by a root-find over `c` in `(0, 20]`.
##' Everything else follows from
##'
##' \deqn{P_f = P_{50}\,\big(\log_2(1/f)\big)^{1/c}}
##'
##' with `f` the remaining fraction (`f = 1 - PLC`), so `log2(1/f)` is 1, 3.0589,
##' 4.3219 and 6.6439 for P50, P88, P95 and P99. From two quantiles,
##' `c = ln(L2/L1) / ln(P_f2/P_f1)`.
##'
##' @section Guards:
##' `c > 0` is required. A `c` below 1 is **warned about rather than refused**: it
##' is possible, but across measured angiosperm stems it usually flags a
##' measurement artefact — native embolism already present at low tension flattens
##' the curve — rather than a real trait. For reference, `c = 2.0` is a reasonable
##' default for angiosperm stems, measured 1.8–2.6 across three species whose P50
##' spanned 2.9–4.3 MPa, i.e. shape roughly conserved while scale varied freely.
##'
##' @param P50,P88,P95,P99 potentials at 50%, 88%, 95% and 99% loss of
##'   conductivity, MPa as positive magnitudes
##' @param px,plc an arbitrary quantile: the potential `px` at which the
##'   proportion `plc` of conductivity is lost. Supply both or neither; the pair
##'   counts as ONE of the two quantities.
##' @param S50 slope at P50, % PLC per MPa. See [weibull_s50()].
##' @param c shape parameter, if you already have it.
##'
##' @return A list with `P50` and `c`, suitable for splicing into
##'   [leaf_traits()].
##' @seealso [psi_at_plc()], [weibull_s50()]
##' @examples
##' # From a published P50 and P88 pair
##' weibull_p50_c(P50 = 3.4, P88 = 5.0)
##'
##' # From P50 and a slope
##' weibull_p50_c(P50 = 3.4, S50 = 39.4)
##'
##' # From an arbitrary quantile plus a shape
##' weibull_p50_c(px = 5.87, plc = 0.95, c = 2.680147)
##' @export
weibull_p50_c <- function(P50 = NULL, P88 = NULL, P95 = NULL, P99 = NULL,
                          px = NULL, plc = NULL, S50 = NULL, c = NULL) {
  if (is.null(px) != is.null(plc)) {
    stop("`px` and `plc` describe one quantile and must be given together.",
         call. = FALSE)
  }
  # The named quantiles and the arbitrary one reduce to the same thing: a
  # potential plus the remaining fraction it belongs to.
  quant <- list()
  for (nm in names(.WEIBULL_F)) {
    v <- get(nm)
    if (!is.null(v)) quant[[nm]] <- c(psi = v, f = .WEIBULL_F[[nm]])
  }
  if (!is.null(px)) quant[["px"]] <- c(psi = px, f = 1 - plc)

  given <- c(names(quant), if (!is.null(S50)) "S50", if (!is.null(c)) "c")
  if (length(given) != 2L) {
    stop("exactly two quantities are needed to determine a curve; got ",
         if (length(given) == 0L) "none" else paste(given, collapse = " and "),
         ". Available: P50, P88, P95, P99, (px, plc), S50, c.", call. = FALSE)
  }
  for (q in quant) {
    if (!(is.numeric(q[["psi"]]) && length(q[["psi"]]) == 1L &&
          is.finite(q[["psi"]]) && q[["psi"]] > 0)) {
      stop("every potential must be a single finite positive magnitude in MPa; ",
           "published values are negative and this package uses |psi|.",
           call. = FALSE)
    }
    if (q[["f"]] <= 0 || q[["f"]] >= 1) {
      stop("`plc` must be in (0, 1).", call. = FALSE)
    }
  }

  out <- if (length(quant) == 2L) {
    .weibull_from_two_quantiles(quant[[1]], quant[[2]])
  } else if (length(quant) == 1L && !is.null(c)) {
    q <- quant[[1]]
    list(P50 = q[["psi"]] * .weibull_L(q[["f"]])^(-1 / c), c = c)
  } else if (length(quant) == 1L && !is.null(S50)) {
    .weibull_from_quantile_and_s50(quant[[1]], S50)
  } else if (!is.null(S50) && !is.null(c)) {
    # S50 = K*c/P50 rearranges directly.
    list(P50 = .WEIBULL_S50_K * c / S50, c = c)
  } else {
    stop("that pair does not determine a curve.", call. = FALSE)
  }

  .weibull_check_p50_c(out$P50, out$c)
  if (out$c < 1) {
    warning("c = ", format(out$c, digits = 4), " is below 1, which flattens the ",
            "curve. That is possible but across measured angiosperm stems it ",
            "usually indicates native embolism already present at low tension ",
            "rather than a real trait -- check the source curve before using it.",
            call. = FALSE)
  }
  out
}

# --- internals ---------------------------------------------------------------

.weibull_check_p50_c <- function(P50, c) {
  if (!(is.numeric(P50) && length(P50) == 1L && is.finite(P50) && P50 > 0)) {
    stop("`P50` must be a single finite positive magnitude in MPa; got ",
         format(P50), call. = FALSE)
  }
  if (!(is.numeric(c) && length(c) == 1L && is.finite(c) && c > 0)) {
    stop("`c` must be a single finite positive number; got ", format(c),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Two quantiles determine both parameters in closed form. ⚠️ They must be
# DIFFERENT quantiles: two potentials at the same fraction either agree (no
# information) or contradict each other, and the log ratio below is 0/0 either way.
.weibull_from_two_quantiles <- function(a, b) {
  La <- .weibull_L(a[["f"]])
  Lb <- .weibull_L(b[["f"]])
  if (isTRUE(all.equal(La, Lb))) {
    stop("the two quantiles are at the same loss of conductivity, so they carry ",
         "one piece of information rather than two.", call. = FALSE)
  }
  if (isTRUE(all.equal(a[["psi"]], b[["psi"]]))) {
    stop("two different quantiles cannot sit at the same potential: the curve is ",
         "strictly increasing in psi.", call. = FALSE)
  }
  cc <- log(Lb / La) / log(b[["psi"]] / a[["psi"]])
  list(P50 = a[["psi"]] * La^(-1 / cc), c = cc)
}

# S50 with a non-P50 quantile is the one transcendental case: substituting
# `P50 = psi_f * L^(-1/c)` into `S50 = K*c/P50` leaves c on both sides inside an
# exponent.
#
# ⚠️ AND IT DOES NOT HAVE A UNIQUE SOLUTION. Written out,
#
#     S50(c) = K * c * L^(1/c) / psi_f
#
# and `d/dc [c*L^(1/c)] = L^(1/c) * (1 - ln(L)/c)`, which is zero at `c = ln L`.
# So S50 is U-SHAPED in c with a minimum there, and any S50 above that minimum is
# reproduced by TWO curves -- one on each branch. This was found by the round-trip
# test failing with "no root in (0, 20]": both bracket ends were positive because
# the true root sat in a dip between them.
#
# Each branch is searched separately, from the analytic turning point rather than a
# guess. The UPPER branch is returned because that is where measured stems sit --
# `ln L` is 1.46 for P95 and 1.12 for P88, against a measured angiosperm range of
# 1.8-2.6 -- and the lower branch is reported so the caller knows it exists rather
# than being handed one of two answers silently.
#
# When the quantile IS P50 the substitution collapses, `L = 1`, and there is no
# root-find and no ambiguity. That is the common case and it stays exact.
.weibull_from_quantile_and_s50 <- function(q, S50) {
  if (!(is.numeric(S50) && length(S50) == 1L && is.finite(S50) && S50 > 0)) {
    stop("`S50` must be a single finite positive slope in % PLC per MPa; got ",
         format(S50), call. = FALSE)
  }
  L <- .weibull_L(q[["f"]])
  if (isTRUE(all.equal(L, 1))) {
    return(list(P50 = q[["psi"]], c = S50 * q[["psi"]] / .WEIBULL_S50_K))
  }

  resid <- function(cc) {
    .WEIBULL_S50_K * cc * L^(1 / cc) / q[["psi"]] - S50
  }
  turn <- log(L)                      # the analytic minimum of c*L^(1/c)
  if (resid(turn) > 0) {
    stop("no curve reproduces S50 = ", format(S50, digits = 6), " at that ",
         "quantile: the smallest slope any curve can have there is ",
         format(resid(turn) + S50, digits = 6), " % PLC per MPa, at c = ",
         format(turn, digits = 4), ". Check that the slope and the potential ",
         "come from the same curve.", call. = FALSE)
  }

  root_in <- function(a, b) {
    if (resid(a) * resid(b) > 0) return(NULL)
    stats::uniroot(resid, c(a, b), tol = 1e-12)$root
  }
  upper <- root_in(turn, 40)
  lower <- root_in(1e-6, turn)
  if (is.null(upper) && is.null(lower)) {
    stop("no shape parameter in (0, 40] reproduces S50 = ", format(S50),
         " at that quantile.", call. = FALSE)
  }
  cc <- if (!is.null(upper)) upper else lower
  if (!is.null(upper) && !is.null(lower)) {
    other <- q[["psi"]] * L^(-1 / lower)
    message("S50 with a quantile other than P50 has two solutions; returning ",
            "the one on the upper branch (c = ", format(cc, digits = 4),
            "). The other is c = ", format(lower, digits = 4), ", P50 = ",
            format(other, digits = 4),
            ". Supply a second quantile instead if you need it resolved.")
  }
  list(P50 = q[["psi"]] * L^(-1 / cc), c = cc)
}
