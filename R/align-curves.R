##' Align the cost curves' level parameters to a common optimum
##'
##' Solves each curve's level parameter so that its stem optimum has the same
##' stomatal conductance as a target, under one reference environment. Aligned
##' curves agree at the reference and differ only in how their optimum moves
##' away from it, so a comparison across curves measures how each RESPONDS to
##' the environment rather than where its parameters happen to put it.
##'
##' The package's own defaults are this function's answer at its defaults:
##' `JS22_gamma` and `CMax_a` in [leaf_traits()], and the `CF77_lambda` and
##' `LeastCost_beta` that [leaf_solve()] and [leaf_batch()] use when those are
##' left `NA`. The defaults' anchor is TF24, which does not move.
##'
##' @section Which curves can be aligned:
##' Four curves have a level parameter: `CF77` (`CF77_lambda`, its price of
##' water), `JS22` (`JS22_gamma`), `CMax` (`CMax_a`, with `CMax_b` held at its
##' trait value) and `LeastCost` (`LeastCost_beta`). Each is solved on a log
##' scale by [stats::uniroot()].
##'
##' The others are not, and this function refuses them:
##'
##' * `SOX`, `JW26` and `ProfitMax` have no level parameter. The shared supply
##'   sets their optimum, and matching them through `kmax` would mean they no
##'   longer share the supply every other curve uses.
##' * `TF24_floor` reads TF24's own cost and adds a price `lambda_o` on top, so
##'   the only `lambda_o` that reaches TF24's conductance is zero, where the
##'   curve is TF24. It keeps no default and still refuses to solve without one.
##'
##' `CF77` aligned to TF24 is an identity rather than a fit: at TF24's emergent
##' `lambda` the two first-order conditions coincide, so the price returned is
##' [operating_point()]'s `lambda_emergent` for TF24 at the reference.
##'
##' @inheritParams leaf_solve
##' @param models which curves to align; any subset of `"CF77"`, `"JS22"`,
##'   `"CMax"` and `"LeastCost"`.
##' @param psi_soil,PPFD,atm_vpd,ca,leaf_temp,atm_o2_kpa,atm_kpa the reference
##'   environment, as in [leaf_solve()]. The defaults are [leaf_solve()]'s own,
##'   with a moderately wet soil.
##' @param supply the supply path. The stem routes need
##'   [leaf_supply_singlelayer()], which is the default.
##' @param target the stomatal conductance to match, mol CO2 m^-2 s^-1. `NULL`
##'   (the default) takes TF24's at the reference, with `traits`.
##' @param tol the tolerance on the log of each level parameter.
##'
##' @return A named numeric vector, one level parameter per curve, named as the
##'   argument or trait that carries it: `CF77_lambda`, `JS22_gamma`, `CMax_a`,
##'   `LeastCost_beta`. Attribute `"target"` holds the conductance matched.
##' @seealso [cost_curve_names()], [leaf_traits()], [leaf_solve()]
##' @examples
##' a <- align_curves()
##' a
##' # The defaults are this answer, so these two agree:
##' leaf_solve(psi_soil = 1, model = "JS22",
##'            supply = leaf_supply_singlelayer())$gc
##' attr(a, "target")
##'
##' # Align to a wetter reference instead, and use the result:
##' w <- align_curves(psi_soil = 0.5)
##' leaf_solve(psi_soil = 0.5, model = "CMax",
##'            supply = leaf_supply_singlelayer(),
##'            traits = leaf_traits(CMax_a = w[["CMax_a"]]))$gc
##' @export
align_curves <- function(models = c("CF77", "JS22", "CMax", "LeastCost"),
                         psi_soil = 1,
                         PPFD = 900,
                         atm_vpd = 2.0,
                         ca = 40.0,
                         leaf_temp = 25.0,
                         atm_o2_kpa = 21.0,
                         atm_kpa = 101.3,
                         leaf_specific_conductance_max = 3.14e-5,
                         traits = leaf_traits(),
                         control = leaf_control(),
                         supply = leaf_supply_singlelayer(),
                         target = NULL,
                         tol = 1e-12) {
  unaligned <- setdiff(models, names(.align_levels))
  if (length(unaligned)) {
    stop("cannot align ", paste0("\"", unaligned, "\"", collapse = ", "),
         ": only ", paste0("\"", names(.align_levels), "\"", collapse = ", "),
         " have a level parameter. SOX, JW26 and ProfitMax have none, and ",
         "TF24_floor reaches TF24's conductance only at lambda_o = 0, where it ",
         "is TF24. See ?align_curves.", call. = FALSE)
  }
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  gc_at <- function(model, extra = list(), tr = traits) {
    args <- c(list(psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd, ca = ca,
                   leaf_temp = leaf_temp, atm_o2_kpa = atm_o2_kpa,
                   atm_kpa = atm_kpa,
                   leaf_specific_conductance_max = leaf_specific_conductance_max,
                   traits = tr, control = control, supply = supply,
                   model = model),
              extra)
    do.call(leaf_solve, args)$gc
  }
  if (is.null(target)) {
    target <- gc_at("TF24")
    if (!(is.finite(target) && target > 0)) {
      stop("TF24 has no positive conductance at this reference (gc = ",
           format(target), "), so there is nothing to align to; choose ",
           "another reference or pass `target`", call. = FALSE)
    }
  }
  if (!(is.numeric(target) && length(target) == 1L && is.finite(target) &&
        target > 0)) {
    stop("`target` must be one positive conductance, mol CO2 m^-2 s^-1",
         call. = FALSE)
  }

  out <- vapply(models, function(model) {
    level <- .align_levels[[model]]
    gc_of <- function(value) {
      if (level$trait) {
        tr <- do.call(leaf_traits,
                      utils::modifyList(unclass(traits),
                                        stats::setNames(list(value), level$name)))
        gc_at(model, tr = tr)
      } else {
        gc_at(model, stats::setNames(list(value), level$name))
      }
    }
    f <- function(log_value) gc_of(exp(log_value)) - target
    bracket <- log(level$bracket)
    ends <- c(f(bracket[[1]]), f(bracket[[2]]))
    if (prod(sign(ends)) > 0) {
      stop(model, " cannot reach a conductance of ", signif(target, 4),
           " at this reference: over ", level$name, " in [",
           paste(signif(level$bracket, 3), collapse = ", "), "] it runs ",
           paste(signif(ends + target, 4), collapse = " to "), ".",
           call. = FALSE)
    }
    exp(stats::uniroot(f, bracket, f.lower = ends[[1]], f.upper = ends[[2]],
                       tol = tol)$root)
  }, numeric(1))
  names(out) <- vapply(models, function(m) .align_levels[[m]]$name,
                       character(1))
  attr(out, "target") <- target
  out
}

# Each alignable curve's level parameter: its name as a leaf_traits() argument
# or a leaf_solve() argument, which of the two it is, and a bracket wide enough
# that the curve's conductance runs from shut to the wet bound across it at
# this package's default leaf. Conductance falls with the first three and
# rises with beta, which weighs capacity against water; uniroot needs only the
# sign change.
.align_levels <- list(
  CF77 = list(name = "CF77_lambda", trait = FALSE, bracket = c(1e3, 1e8)),
  JS22 = list(name = "JS22_gamma", trait = TRUE, bracket = c(1e-4, 1e4)),
  CMax = list(name = "CMax_a", trait = TRUE, bracket = c(1e-4, 1e4)),
  LeastCost = list(name = "LeastCost_beta", trait = FALSE,
                   bracket = c(1e-1, 1e7))
)
