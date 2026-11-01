## How the two solvers converge, for vignette("the-models"). A search comparing
## values of J stops at a floor, sqrt(4.eps.|J|/k) for curvature k (drawn on
## panel a), because J is flat at its maximum; a root-find on dJ/dpsi, steep
## there, does not. Both get the same bracket, which flatters the search.
##
## Only the outer iterations run in R, since C++ does not return its iterates: J
## and dJ/dpsi are the bound `evaluate_root_collar_psi()` and
## `dprofit_droot_collar_psi()`. The root-find is R's `uniroot` (Brent) traced, and
## the search is textbook golden section.

.sc_col <- list(
  search = "firebrick",
  rootfind = "steelblue4",
  floor = "grey60",
  guide = "grey85"
)

## A sign-changing bracket for dJ/dpsi: the first crossing on a scan, or NULL when
## the optimum is pinned, a real outcome for which the caller picks other drivers.
## Points outside the feasible collar interval, which is narrower than
## [psi_soil, psi_crit], raise a domain condition; the scan drops them.
.sc_bracket <- function(l, n = 60) {
  lo <- l$psi_soil_[1] + 1e-6
  hi <- l$psi_crit - 1e-6
  x <- seq(lo, hi, length.out = n)
  d <- vapply(x, function(p) {
    tryCatch(l$dprofit_droot_collar_psi(p), error = function(e) NA_real_)
  }, numeric(1))
  ok <- which(is.finite(d))
  if (length(ok) < 2) return(NULL)
  x <- x[ok]
  d <- d[ok]
  i <- which(d[-length(d)] > 0 & d[-1] < 0)
  if (!length(i)) return(NULL)
  c(x[i[1]], x[i[1] + 1])
}

## Golden section on J, recording the interval's best point after every
## evaluation. Terminates on interval WIDTH -- which is the whole point: it
## never resolves the location better than the bracket it has left.
.sc_trace_search <- function(l, a, b, n_eval = 60) {
  phi <- (sqrt(5) - 1) / 2
  f <- function(p) {
    l$evaluate_root_collar_psi(p)
    l$profit_
  }
  c1 <- b - phi * (b - a)
  c2 <- a + phi * (b - a)
  f1 <- f(c1)
  f2 <- f(c2)
  out <- data.frame(eval = 2L, x = if (f1 > f2) c1 else c2)
  while (out$eval[nrow(out)] < n_eval) {
    if (f1 > f2) {
      b <- c2
      c2 <- c1
      f2 <- f1
      c1 <- b - phi * (b - a)
      f1 <- f(c1)
    } else {
      a <- c1
      c1 <- c2
      f1 <- f2
      c2 <- a + phi * (b - a)
      f2 <- f(c2)
    }
    out <- rbind(out, data.frame(eval = out$eval[nrow(out)] + 1L,
                                 x = if (f1 > f2) c1 else c2))
  }
  out
}

## Brent's zeroin on dJ/dpsi, traced. The wrapper records every point uniroot
## asks for, so the sequence IS the solver's iterate sequence -- no
## reimplementation, and the evaluation count is exact.
.sc_trace_rootfind <- function(l, a, b, tol = 1e-15) {
  seen <- numeric(0)
  g <- function(p) {
    seen[[length(seen) + 1L]] <<- p
    l$dprofit_droot_collar_psi(p)
  }
  invisible(stats::uniroot(g, c(a, b), tol = tol, maxiter = 200))
  data.frame(eval = seq_along(seen), x = seen)
}

##' Convergence of the two ways to locate an optimum
##'
##' Draws two panels for `vignette("the-models")`: the error in the located
##' operating point, and the first-order residual at it, both against the number
##' of model evaluations spent. One series solves `dJ/dpsi = 0` by Brent's
##' method; the other searches `J` for its maximum by golden section.
##'
##' The iterations run in R against the package's own objective and derivative;
##' only the outer loop is reproduced, because the C++ solvers do not return
##' their iterates. See the file's header comment.
##'
##' @param psi_soil,PPFD,atm_vpd Drivers giving an interior optimum.
##' @param n_eval Evaluations to trace for the search.
##' @return Invisibly, a list with the two traces and the reference optimum.
##' @export
plot_solver_convergence <- function(psi_soil = 1.0, PPFD = 1500,
                                    atm_vpd = 1.5, n_eval = 60) {
  l <- leaf_model()
  set_drivers(l, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd)

  br <- .sc_bracket(l)
  if (is.null(br)) {
    stop("no interior stationary point at these drivers; pick a wetter soil ",
         "or brighter light", call. = FALSE)
  }

  ## The reference. Converged far past either trace, so the error axis measures
  ## the traces rather than the reference.
  star <- stats::uniroot(function(p) l$dprofit_droot_collar_psi(p),
                         br, tol = .Machine$double.eps, maxiter = 200)$root

  search <- .sc_trace_search(l, br[1], br[2], n_eval)
  root <- .sc_trace_rootfind(l, br[1], br[2])

  err <- function(tr) abs(tr$x - star)
  res <- function(tr) abs(vapply(tr$x, function(p)
    l$dprofit_droot_collar_psi(p), numeric(1)))

  ## Where a value-comparing search must stop. k is the curvature at the
  ## optimum, taken from the derivative the root-find already uses.
  h <- 1e-4
  k <- -(l$dprofit_droot_collar_psi(star + h) -
           l$dprofit_droot_collar_psi(star - h)) / (2 * h)
  l$evaluate_root_collar_psi(star)
  floor_x <- sqrt(4 * .Machine$double.eps * abs(l$profit_) / k)

  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1.0),
                      mgp = c(2.5, 0.7, 0), las = 1, bty = "l",
                      cex.axis = 0.85)
  on.exit(graphics::par(op), add = TRUE)

  .sc_panel(search$eval, err(search), root$eval, err(root),
            ylab = expression("|" * psi[n] - psi * "*|" ~ (MPa)),
            main = "(a) error in the located point", floor_x = floor_x)
  .sc_panel(search$eval, res(search), root$eval, res(root),
            ylab = expression("|" * dJ / d * psi * "|" ~ "at" ~ psi[n]),
            main = "(b) first-order residual there")

  invisible(list(search = search, rootfind = root, star = star,
                 bracket = br, curvature = k, search_floor = floor_x))
}

## One panel: both series on a log y, a legend on the first only, and the floor
## drawn where the search stops improving.
.sc_panel <- function(xs, ys, xr, yr, ylab, main, floor_x = NULL) {
  pos <- function(v) ifelse(v <= 0, NA, v)
  ys <- pos(ys)
  yr <- pos(yr)
  ylim <- range(c(ys, yr, floor_x), na.rm = TRUE)
  ylim[1] <- max(ylim[1], 1e-17)
  graphics::plot(NA, xlim = c(0, max(xs)), ylim = ylim, log = "y",
                 xlab = "model evaluations", ylab = ylab, main = main,
                 cex.main = 0.95, font.main = 1)
  graphics::grid(col = .sc_col$guide, lty = 1)
  if (!is.null(floor_x)) {
    graphics::abline(h = floor_x, col = .sc_col$search, lty = 2)
    graphics::text(max(xs), floor_x, adj = c(1, -0.5), cex = 0.75,
                   col = .sc_col$search,
                   labels = expression(sqrt(4 * epsilon * "|J|" / k) ~
                                         "-- the limit of comparing J"))
  }
  graphics::lines(xs, ys, col = .sc_col$search, lwd = 2)
  graphics::points(xs, ys, col = .sc_col$search, pch = 16, cex = 0.4)
  graphics::lines(xr, yr, col = .sc_col$rootfind, lwd = 2)
  graphics::points(xr, yr, col = .sc_col$rootfind, pch = 16, cex = 0.9)
  if (grepl("^\\(a\\)", main)) {
    graphics::legend("topright", bty = "n", cex = 0.85, lwd = 2,
                     col = c(.sc_col$search, .sc_col$rootfind),
                     legend = c("search J for its maximum (golden section)",
                                expression("solve" ~ dJ / d * psi == 0 ~
                                             "(Brent)")))
  }
}
