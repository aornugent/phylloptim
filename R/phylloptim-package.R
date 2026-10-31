#' phylloptim: leaf gas exchange with explicit plant hydraulics
#'
#' A leaf gas-exchange model in which stomatal behaviour emerges from
#' hydraulics. Farquhar-von Caemmerer-Berry photosynthesis is coupled to an
#' explicit soil-root-stem-leaf water path, and the operating point is the
#' maximum of an optimality objective over the root-collar water potential.
#' Nine optimality models share one solver; [cost_curve_names()] lists them.
#' The model is a set of R-free C++ headers under `inst/include`, and this
#' package is an R layer on top of them.
#'
#' @section Where to start:
#' * Solve a leaf in one call: [leaf_solve()], with [leaf_traits()] and
#'   [leaf_control()].
#' * Keep a leaf and drive it yourself: [leaf_model()], [set_drivers()],
#'   [set_traits()] and [operating_point()].
#' * Trait gradients: [leaf_gradient()] for one operating point, and
#'   [leaf_batch()] with [leaf_gradient_batch()] for a fit over many
#'   observations.
#' * Thermal damage: [leaf_damage_trajectory()] and [leaf_damage_gradient()].
#' * Vulnerability curves from published quantities: [weibull_p50_c()] and
#'   [psi_at_plc()].
#'
#' @section Vignettes:
#' * `vignette("phylloptim")` -- a tour of the R interface.
#' * `vignette("the-models")` -- the nine optimality models and the maths
#'   behind them.
#' * `vignette("fitting")` -- fitting traits to data.
#' * `vignette("cpp-interface")` -- using the model from C++ with no R.
#' * `vignette("derivative-surface")` -- differentiating the operating point
#'   inside a larger model.
#'
#' @keywords package
"_PACKAGE"
