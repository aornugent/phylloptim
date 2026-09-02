# Telling infeasibility apart from a caller error (#57).
#
# A parameter proposal during a fit WILL reach operating points the solve cannot
# handle -- that is what an optimiser does -- and those rows have to cost themselves
# rather than the whole likelihood evaluation. Until now every failure arrived as a
# plain R error, so the only way to separate "this operating point is infeasible"
# from "the driver column is all NA" was to match on message text.
#
# ⚠️ Message text is not a stable interface, and that is not hypothetical here: the
# out-of-domain wording was rewritten by #65 and again by #79, and #92 added a third
# variant. Any caller matching on prose was silently broken by each of those, and an
# improvement to an error message should not break a fit.
#
# THE TWO LAYERS, and why it cannot be one. C++ tags the exit with a token in the
# message (`util::stop_infeasible`, see util.hpp for why a derived C++ type is not
# enough); R reads the token and re-signals with a condition class. The join could
# not be at the boundary itself: that is RcppR6-generated and must not be
# hand-edited, and `util.hpp` must not reach for Rcpp (hazard 9, enforced by CI
# building the headers on runners with no R). So the classification lives where the
# knowledge is and the translation lives here.

# The token C++ writes. Kept in one place because the test suite asserts that the
# codes appearing in the installed headers are exactly the ones documented below --
# a code added in C++ and not here would otherwise be an undocumented part of the
# API that a caller discovers from a failure.
.infeasible_pattern <- "^\\[phylloptim:infeasible:([a-z_]+)\\]\\s*"

##' Which failures are reported as infeasible
##'
##' The `code` values a `phylloptim_infeasible` condition can carry, and what each
##' one means. Use these rather than matching on the message text, which changes
##' whenever an error message is improved.
##'
##' @section What is deliberately NOT here:
##' Every input-validation failure. A mismatched `psi_soil`/`soil_depth` length, a
##' non-finite driver, a negative water potential, an unusable `R_d_25`, a
##' misconfigured supply path: those stay ordinary errors, and that is the whole
##' safety property of this list. Classifying one of them as infeasible would let it
##' be swallowed by the very `tryCatch()` this exists to enable, and a fit would
##' report a plausible likelihood over whichever rows survived — which is exactly
##' the failure #39 rejected silent NA rows to prevent.
##'
##' So the test a site has to pass is narrow: it must be reachable from a
##' **well-formed call with in-range data**, because the parameters or the state make
##' the operating point unsolvable. Anything reachable only by handing the model
##' something malformed is not infeasibility, it is a bug in the caller, and it
##' should be loud.
##'
##' @section Three judgement calls left ordinary, and why:
##' Of the 64 `stop` sites in the headers, nineteen are classified. Three
##' of the rest are genuinely arguable, and all three were left ordinary because that
##' is the direction whose failure mode is a loud error rather than a quiet wrong
##' answer. They are listed here so the decision can be overruled on purpose rather
##' than rediscovered:
##'
##' * ~~**`psi_crit` past the stem curve's P99.**~~ No longer reachable, and so no
##'   longer a decision: `psi_crit` is the 5%-conductivity quantile of the stem
##'   curve rather than a trait beside it, so it scales with the curve and lands
##'   inside the domain for every shape parameter. An optimiser moving `stem_P50`
##'   cannot walk out of it.
##'   Constrain the pair instead — `vignette("fitting")` derives both from a P50/P88.
##' * **`method = "ift"` requested at a point where the IFT does not hold**
##'   (`leaf_gradient`). The default `method = "auto"` never reaches it, so getting
##'   here means the caller overrode the choice and the answer is "not at this point".
##' * **Trait sign and magnitude checks** (`psi_crit`, `stem_b`, `root_b`,
##'   `root_psi_crit` must be positive). Reachable by an unbounded proposal, but a
##'   negative `stem_b` is not an infeasible operating point, it is not a stem.
##'
##' Every other unclassified site is input or configuration validation — vector
##' lengths, non-finite drivers, an unusable supply path — where there is no argument
##' for infeasibility at all.
##'
##' @return A named character vector: codes, with a one-line description each.
##' @seealso [with_phylloptim_conditions()]
##' @examples
##' leaf_infeasible_codes()
##' @export
leaf_infeasible_codes <- function() {
  c(collar_bracket =
      paste("the soil-to-collar continuity root-find failed on its bracket, so",
            "no collar potential equilibrates the column"),
    collar_solve =
      paste("profit or stem potential came back non-finite inside the feasible",
            "collar interval, so the maximisation has nothing to maximise"),
    ci_solve =
      paste("the nested ci root-find failed, and not in the way that means",
            "shut-down (which is an answer, not an error)"),
    stem_curve_domain =
      paste("a vulnerability spline was asked for a potential outside the domain",
            "its traits give it"),
    uptake =
      paste("the soil-to-root-collar flux is non-finite, or a layer's",
            "vulnerability weighting is unusable, at this candidate collar"),
    root_find_iterations =
      paste("a one-dimensional solver hit max_iterations without converging"),
    gradient_active_set =
      paste("perturbing a parameter moved the feasible collar interval past the",
            "operating point, so the point sits on an active-set boundary and",
            "cannot be differentiated there"))
}

# The code, or NA if the message carries no token. NA is the safe answer and the
# common one: anything without a token stays an ordinary error.
.infeasible_code <- function(msg) {
  m <- regmatches(msg, regexec(.infeasible_pattern, msg))[[1L]]
  if (length(m) < 2L) NA_character_ else m[[2L]]
}

##' Run leaf code with classified error conditions
##'
##' Evaluates `expr` and, if it fails with an error the model has tagged as
##' infeasible, re-signals it as a `phylloptim_infeasible` condition carrying a
##' machine-readable `code`. Anything else propagates unchanged.
##'
##' [leaf_solve()] and [leaf_gradient_batch()] already do this, so wrap
##' `expr` yourself when you are **driving a `Leaf` directly** —
##' `$find_root_collar_psi()` is generated glue and cannot signal a classified
##' condition on its own — or when you are calling [leaf_gradient()].
##'
##' ⚠️ [leaf_gradient()] is **not** routed through this, and that is a gap rather
##' than a design: it is one long function that depends on `missing()` for two of its
##' arguments, so making it delegate to an inner implementation would change argument
##' matching in ways worth a separate change. Wrap the call yourself until then.
##' [leaf_gradient_batch()] is the path a fit should be on anyway — one boundary
##' crossing rather than one per observation — and it is routed.
##'
##' @section Why this and not a per-row NA mode:
##' Returning `NA` rows instead of throwing is the other half of #57 and is not built
##' yet, deliberately: it must not be built before the classification in
##' [leaf_infeasible_codes()] has been reviewed, or the opt-in would return `NA` for
##' programming errors too and a guard written against "most rows failed" would be
##' applied to bugs as well as to infeasibility. [leaf_gradient_batch()] already
##' reports per-row `status`, which is the same idea where it was safe to have it.
##'
##' @param expr An expression to evaluate.
##' @return The value of `expr`.
##' @seealso [leaf_infeasible_codes()]
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' # An ordinary solve: nothing to catch.
##' with_phylloptim_conditions(l$find_root_collar_psi())
##'
##' # A caller's per-row guard, written against the code rather than the message.
##' tryCatch(
##'   with_phylloptim_conditions(l$find_root_collar_psi()),
##'   phylloptim_infeasible = function(e) message("infeasible: ", e$code))
##' @export
with_phylloptim_conditions <- function(expr) {
  tryCatch(expr, error = function(e) {
    code <- .infeasible_code(conditionMessage(e))
    if (is.na(code)) {
      stop(e)
    }
    stop(.infeasible_condition(e, code))
  })
}

# `phylloptim_error` sits between the specific class and `error` so that a caller can
# catch everything this package classifies without enumerating codes, and so that a
# future second category (there is no second category yet) does not need every
# existing handler rewritten.
.infeasible_condition <- function(e, code) {
  structure(
    class = c("phylloptim_infeasible", "phylloptim_error", "error", "condition"),
    list(message = conditionMessage(e), call = conditionCall(e), code = code))
}
