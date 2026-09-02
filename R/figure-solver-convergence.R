## How the two solvers converge, for vignette("the-models").
##
## The claim being illustrated is stronger than "one converges faster". The
## search does not merely converge slowly -- it STOPS, at a floor no number of
## evaluations gets past, because it compares values of J and J is flat at its
## maximum. Near the optimum J changes by less than a rounding error while psi
## still moves, so the comparison that drives the search becomes noise. The floor
## sits at sqrt(4.eps.|J|/k) for curvature k, which is drawn on panel (a).
##
## A root-find has no such limit: it is looking for a ZERO of a function that is
## steep exactly where J is flat, so the quantity it compares against zero still
## carries information at 1e-15.
##
## Both solvers are given THE SAME sign-changing bracket, which flatters the
## search -- it is handed the basin rather than having to find it, and locating
## the basin is what the 64-point scan in the psi_stem entry points is for.
##
## ⚠️ THE ITERATIONS RUN IN R, THE MODEL DOES NOT. Both traces call the real
## bound methods -- `evaluate_root_collar_psi()` for J and
## `dprofit_droot_collar_psi()` for dJ/dpsi -- so the objective, its derivative
## and every tolerance beneath them are the package's own. Only the outer
## iteration is reproduced here, because C++ does not hand its iterates back.
## The root-find is R's `uniroot` (Brent's zeroin) traced through a counting
## wrapper, so it is a real Brent trace and not a reimplementation; the search is
## textbook golden section, which is what this package used before the collar
## solve moved to the first-order condition.
##
## Base graphics, vector device, everything a separate object if the figure ever
## needs hand-polish -- same rules as figure-overview.R.

.sc_col <- list(
  search = "firebrick",
  rootfind = "steelblue4",
  floor = "grey60",
  guide = "grey85"
)

## A sign-changing bracket for dJ/dpsi: scan the feasible interval and take the
## first crossing. Returns NULL when the optimum is pinned rather than interior,
## which is a real outcome and not a failure -- the caller picks other drivers.
##
## The feasible collar interval is narrower than [psi_soil, psi_crit] -- the
## supply path's own bounds sit inside it -- and asking outside it raises a
## domain condition rather than returning a non-finite number, so the scan
## catches and infeasible points drop out.
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

## One panel. Both series on a log y, a legend only on the first, and the floor
## drawn where the traces stop improving -- because they DO stop, and pretending
## otherwise would be the one dishonest thing this figure could do.
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
