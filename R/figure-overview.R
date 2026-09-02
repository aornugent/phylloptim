## An overview figure for vignette("the-models"): the real curves on the left,
## the causal wiring in the middle, the physical path on the right.
##
## Everything here is base graphics on purpose. The output is a vector device
## (SVG in the vignette, PDF for a manuscript), so every line and label survives
## into Illustrator as a separate object if the figure ever needs hand-polish.
##
## The three zones are drawn by three functions below, each into a 0-100 x 0-100
## user space so the layouts are readable as coordinates rather than as fractions
## of whatever device happens to be open.

## ---------------------------------------------------------------------------
## palette -- shared with the panels, which is the point
##
## The flow diagram is colour-coded to the profit panel: green boxes are the
## green curve, the red box is the red curve. A reader who has understood panel
## (a) can then read the wiring without a legend.

.fig_col <- list(
  benefit = "forestgreen",
  benefit_fill = "#edf5ed",
  cost = "firebrick",
  cost_fill = "#f9ecea",
  profit = "black",
  profit_fill = "#ececec",
  gate = "grey45",
  gate_fill = "#f7f7f7",
  supply = "steelblue4",
  plant = "grey88"
)

## ---------------------------------------------------------------------------
## drawing primitives

## The largest cex up to `cex` at which `txt` fits inside `width` user units.
##
## Boxes are laid out in fixed user coordinates but text is sized in cex, so the
## two part company as soon as the figure is drawn at a different size or the
## caller raises `cex` -- and what that looks like is box labels running out
## through their borders. Measuring is cheap and removes the whole class.
.fig_fit_cex <- function(txt, cex, width, pad = 0.92) {
  lines <- if (is.character(txt)) {
    unlist(strsplit(txt, "\n", fixed = TRUE))
  } else {
    txt
  }
  w <- max(strwidth(lines, cex = cex))
  if (!is.finite(w) || w <= width * pad) cex else cex * width * pad / w
}

## A labelled box. `main` is the equation, `sub` an optional gloss underneath in
## smaller italic, `tag` a panel cross-reference set in the bottom-right corner.
.fig_box <- function(x0, x1, y, h, main, sub = NULL, tag = NULL,
                     border = "black", fill = "white", lty = 1, cex = 1,
                     main_cex = 0.80, sub_cex = 0.66) {
  rect(x0, y - h / 2, x1, y + h / 2, border = border, col = fill, lty = lty,
       lwd = 1.1)
  xm <- (x0 + x1) / 2
  w <- x1 - x0
  mc <- .fig_fit_cex(main, main_cex * cex, w)
  if (is.null(sub)) {
    text(xm, y, main, cex = mc)
  } else {
    text(xm, y + h * 0.20, main, cex = mc)
    text(xm, y - h * 0.24, sub, cex = .fig_fit_cex(sub, sub_cex * cex, w),
         font = 3, col = "grey25")
  }
  if (!is.null(tag)) {
    text(x1 - 1.2, y + h / 2 - 1.4, tag, cex = sub_cex * cex, font = 2,
         col = "grey35", adj = c(1, 1))
  }
  invisible(NULL)
}

## Straight arrow with a consistent head everywhere in the figure.
.fig_arrow <- function(x0, y0, x1, y1, col = "black", lwd = 1.3, lty = 1,
                       len = 0.06) {
  arrows(x0, y0, x1, y1, length = len, angle = 20, col = col, lwd = lwd,
         lty = lty, xpd = NA)
}

## A resistor zigzag along a segment, occupying the middle `frac` of it.
.fig_resistor <- function(x0, y0, x1, y1, n = 5, amp = 1.6, frac = 0.6,
                          col = "black", lwd = 1.2) {
  ## unit vector along the segment, and its normal
  dx <- x1 - x0
  dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len
  uy <- dy / len
  nx <- -uy
  ny <- ux

  pad <- (1 - frac) / 2
  ta <- pad
  tb <- 1 - pad
  ## zigzag parameter: n up-down excursions between ta and tb
  ts <- seq(ta, tb, length.out = 2 * n + 2)
  off <- c(0, rep(c(amp, -amp), n), 0)
  px <- x0 + ts * dx + off * nx
  py <- y0 + ts * dy + off * ny

  segments(x0, y0, x0 + ta * dx, y0 + ta * dy, col = col, lwd = lwd)
  lines(px, py, col = col, lwd = lwd)
  segments(x0 + tb * dx, y0 + tb * dy, x1, y1, col = col, lwd = lwd)
  invisible(NULL)
}

## Open a schematic panel: no axes, a 0-100 x 0-100 space, tight margins.
.fig_panel <- function(mar = c(0.2, 0.2, 0.2, 0.2)) {
  par(mar = mar)
  plot.new()
  plot.window(c(0, 100), c(0, 100), xaxs = "i", yaxs = "i")
  invisible(NULL)
}

## How many y-units are physically as long as one x-unit, in the panel that is
## currently open.
##
## Both schematic panels use a 0-100 x 0-100 space, so a shape specified in user
## units is stretched by whatever aspect the device gives the panel -- and these
## panels are strongly portrait. Multiplying a y-extent by this makes a shape
## come out the way it was drawn on paper: `.fig_aspect() * r` is the y-radius of
## a circle whose x-radius is `r`. Without it the leaf and the resistors distort
## with the figure's dimensions, which is exactly what a reusable figure must not
## do.
.fig_aspect <- function() {
  pin <- par("pin")
  pin[1] / pin[2]
}

## ---------------------------------------------------------------------------
## zone (a) -- the profit panel
##
## Benefit, cost and their difference against the dial. This is the picture every
## other part of the model is in service of, which is why it is top-left.

.fig_profit_panel <- function(l, n = 300, cex = 1) {
  psi <- seq(l$psi_soil_[1], l$psi_crit, length.out = n)
  benefit <- vapply(psi, function(p) {
    l$profit_psi_stem_TF(p, 1.0)
    l$assim_colimited_
  }, numeric(1))
  cost <- vapply(psi, l$hydraulic_cost_TF, numeric(1))
  profit <- benefit - cost

  l$set_model("TF24", "stem"); l$optimise()
  star <- l$opt_psi_stem_
  jstar <- l$profit_

  par(mar = c(3.4, 3.6, 2.0, 0.8), mgp = c(2.1, 0.6, 0), tcl = -0.3,
      cex.axis = 0.78 * cex, cex.lab = 0.85 * cex)
  plot(psi, benefit, type = "n", ylim = range(c(0, benefit, cost, profit)),
       xlab = expression(psi ~ " (MPa)"),
       ylab = expression("C gained or given up  (" * mu * "mol m"^-2 *
                           "s"^-1 * ")"),
       las = 1)

  ## the optimum, marked before the curves so the curves sit on top
  segments(star, par("usr")[3], star, jstar, lty = 3, col = "grey40")
  lines(psi, benefit, lwd = 2, col = .fig_col$benefit)
  lines(psi, cost, lwd = 2, col = .fig_col$cost)
  lines(psi, profit, lwd = 2.6, col = .fig_col$profit)
  points(star, jstar, pch = 21, bg = "white", col = "black", cex = 0.9 * cex,
         lwd = 1.4)

  ## label the curves in place rather than in a legend box
  i <- round(0.86 * n)
  text(psi[i], benefit[i], "benefit  B", col = .fig_col$benefit,
       cex = 0.68 * cex, adj = c(1, -0.6), font = 2)
  text(psi[i], cost[i], "cost  C", col = .fig_col$cost, cex = 0.68 * cex,
       adj = c(1, 1.5), font = 2)
  j <- round(0.72 * n)
  text(psi[j], profit[j], "profit  J = B - C", col = .fig_col$profit,
       cex = 0.68 * cex, adj = c(0.5, -1.0), font = 2)
  text(star, par("usr")[3], expression(psi * "*"), adj = c(-0.3, -0.4),
       cex = 0.78 * cex, col = "grey20")

  mtext("a)  where to operate", side = 3, line = 0.5, adj = 0, font = 2,
        cex = 0.82 * cex)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## zone (b) -- the A/c_i panel
##
## The inner root-find, drawn the way Westoby et al. drew it: biochemical demand
## rising, diffusive supply falling, and the operating point where they cross.

.fig_aci_panel <- function(l, n = 200, cex = 1) {
  l$set_model("TF24", "stem"); l$optimise()
  gc <- l$stom_cond_CO2_
  ca <- l$ca_
  P <- l$atm_kpa_ * 1000
  ci_star <- l$ci_

  supply <- function(ci) gc * (ca - ci) / P * 1e6
  ci <- seq(0.02, ca, length.out = n)
  Ac <- vapply(ci, l$assim_rubisco_limited, numeric(1))
  Aj <- vapply(ci, l$assim_electron_limited, numeric(1))
  A <- vapply(ci, l$assim_colimited, numeric(1))
  sup <- supply(ci)
  Astar <- l$assim_colimited(ci_star)

  par(mar = c(3.4, 3.6, 2.0, 0.8), mgp = c(2.1, 0.6, 0), tcl = -0.3,
      cex.axis = 0.78 * cex, cex.lab = 0.85 * cex)
  plot(ci, A, type = "n", ylim = c(0, max(c(Ac, Aj, sup))),
       xlab = expression(c[i] ~ " (Pa)"),
       ylab = expression("assimilation  (" * mu * "mol m"^-2 ~ "s"^-1 * ")"),
       las = 1)

  segments(ci_star, 0, ci_star, Astar, lty = 3, col = "grey40")
  segments(par("usr")[1], Astar, ci_star, Astar, lty = 3, col = "grey40")
  lines(ci, Ac, lwd = 1.2, col = .fig_col$benefit, lty = 2)
  lines(ci, Aj, lwd = 1.2, col = .fig_col$benefit, lty = 3)
  lines(ci, A, lwd = 2.4, col = .fig_col$benefit)
  lines(ci, sup, lwd = 2, col = .fig_col$supply)
  points(ci_star, Astar, pch = 21, bg = "white", col = "black",
         cex = 0.9 * cex, lwd = 1.4)

  text(ca, supply(ca * 0.80), expression("supply, slope " * -g[c] / P),
       col = .fig_col$supply, cex = 0.66 * cex, adj = c(1, -0.5), font = 2)
  k <- round(0.90 * n)
  text(ci[k], Ac[k], expression(A[c]), col = .fig_col$benefit,
       cex = 0.68 * cex, adj = c(0.5, -0.5))
  text(ci[k], Aj[k], expression(A[j]), col = .fig_col$benefit,
       cex = 0.68 * cex, adj = c(0.5, 1.3))
  m <- round(0.34 * n)
  text(ci[m], A[m], "colimited, net", col = .fig_col$benefit,
       cex = 0.66 * cex, adj = c(0.5, 2.1), font = 2)
  text(ci_star, 0, expression(c[i] * "*"), adj = c(-0.3, -0.4),
       cex = 0.78 * cex, col = "grey20")

  mtext(expression(bold("b)  the inner root-find:  supply = demand")),
        side = 3, line = 0.5, adj = 0, cex = 0.82 * cex)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## zone (c) -- the wiring
##
## One dial, two branches, and the two optional gates that reach across from the
## cost side into the benefit side. The gates are dashed because both default off.

.fig_flow_panel <- function(cex = 1) {
  .fig_panel(c(0.4, 0.6, 0.4, 0.6))

  lx0 <- 2
  lx1 <- 44                              # benefit column
  rx0 <- 56
  rx1 <- 98                              # cost / gate column
  lmid <- (lx0 + lx1) / 2
  rmid <- (rx0 + rx1) / 2
  gmid <- (lx1 + rx0) / 2                # the gutter, where the gates reach across

  text(lx0, 97.5, "c)  how the parts interact", adj = c(0, 0.5), font = 2,
       cex = 0.82 * cex)

  ## the dial, spanning both columns because it feeds both
  .fig_box(lx0, rx1, 88.5, 7,
           expression(psi ~ "-- the single decision variable"),
           "stem potential, or the root collar",
           fill = .fig_col$profit_fill, cex = cex)

  ## benefit branch --------------------------------------------------------
  .fig_arrow(lmid, 84.8, lmid, 80.1)
  .fig_box(lx0, lx1, 75, 9,
           expression(E == k[max] * integral(f(u) * ~ du, psi[up], psi)),
           "hydraulic supply -- saturates, because f falls",
           border = .fig_col$benefit, fill = .fig_col$benefit_fill, cex = cex)

  .fig_arrow(lmid, 70.3, lmid, 64.1)
  .fig_box(lx0, lx1, 59, 9,
           expression(g[c] == P * E * M[w] / (1.6 * D[leaf])),
           "one conductance, two fluxes",
           border = .fig_col$benefit, fill = .fig_col$benefit_fill, cex = cex)

  .fig_arrow(lmid, 54.3, lmid, 47.1)
  .fig_box(lx0, lx1, 41.5, 10.5,
           expression(A(c[i]) == g[c] * (c[a] - c[i]) / P),
           "supply = demand, a root-find at\nevery candidate psi  ->  B",
           border = .fig_col$benefit, fill = .fig_col$benefit_fill,
           tag = "[b]", cex = cex)

  ## cost branch and the two gates ----------------------------------------
  .fig_arrow(rmid, 84.8, rmid, 80.1)
  .fig_box(rx0, rx1, 75, 9,
           "energy balance   (off by default)",
           "T_leaf follows E, then D_leaf follows T_leaf",
           border = .fig_col$gate, fill = .fig_col$gate_fill, lty = 2,
           cex = cex)

  .fig_arrow(rmid, 70.3, rmid, 64.1, col = .fig_col$gate)
  .fig_box(rx0, rx1, 59, 9,
           "thermal cost   (off by default)",
           "a penalty in T_leaf, and x (1 - TC) on J_max",
           border = .fig_col$gate, fill = .fig_col$gate_fill, lty = 2,
           cex = cex)

  .fig_arrow(rmid, 54.3, rmid, 47.1, col = .fig_col$gate)
  .fig_box(rx0, rx1, 41.5, 10.5,
           expression(C(psi) ~ "-- the cost of the tension"),
           paste0("TF24 | Cowan-Farquhar | ProfitMax\n",
                  "the models differ here, and nowhere else"),
           border = .fig_col$cost, fill = .fig_col$cost_fill, cex = cex)

  ## the cross-links: the part the figure exists for -----------------------
  ##
  ## The energy balance is a CYCLE, so it is drawn as the two arrows it is: the
  ## gate reads E, and what it returns re-enters at g_c. The thermal cost gets a
  ## third arrow, reaching past the gutter into the benefit column, because it is
  ## the only gate that moves both sides of the ledger.

  .fig_arrow(lx1 + 0.5, 77, rx0 - 0.5, 77, col = .fig_col$gate, lwd = 1.1)
  .fig_arrow(rx0 - 0.5, 72.5, lx1 + 0.5, 62, col = .fig_col$gate, lwd = 1.1)
  text(gmid, 68, "solved, not evaluated", srt = 90, cex = 0.55 * cex, font = 3,
       col = .fig_col$gate)

  .fig_arrow(rx0 - 0.5, 56, lx1 + 0.5, 45.5, col = .fig_col$gate, lwd = 1.1)
  text(gmid, 50.5, "scales J_max", srt = 90, cex = 0.55 * cex, font = 3,
       col = .fig_col$gate)

  ## converge on the objective -------------------------------------------
  .fig_box(lx0, rx1, 24, 8,
           expression(J(psi) == B(psi) - C(psi)), NULL,
           fill = .fig_col$profit_fill, cex = cex)
  .fig_arrow(lmid, 36.1, lmid, 28.3, col = .fig_col$benefit)
  .fig_arrow(rmid, 36.1, rmid, 28.3, col = .fig_col$cost)

  .fig_box(lx0, rx1, 8, 9,
           expression(psi * "*" == argmax ~ J),
           "solve dJ/dpsi = 0, or report the bound it pinned to",
           fill = .fig_col$profit_fill, cex = cex, tag = "[a]")
  .fig_arrow(lmid, 19.9, lmid, 12.7)

  ## lambda: the quantity that makes the three cost curves comparable
  text(lmid + 5, 16.3,
       "lambda = (dC/dpsi)/(dE/dpsi)  --  the common currency",
       adj = c(0, 0.5), cex = 0.58 * cex, font = 3, col = "grey25")

  invisible(NULL)
}

## ---------------------------------------------------------------------------
## zone (d) -- the physical path
##
## The soil-to-leaf path the maths above is a description of, with numbered notes
## below it separating what is physics from what is a modelling choice.

.fig_plant_panel <- function(cex = 1) {
  .fig_panel(c(0.4, 0.6, 0.4, 0.6))
  k <- .fig_aspect()   # y-units per x-unit at equal physical length

  ## A numbered marker, tying a place on the drawing to a note below it.
  marker <- function(x, y, n) {
    th <- seq(0, 2 * pi, length.out = 60)
    r <- 2.6
    polygon(x + r * cos(th), y + r * k * sin(th), col = "grey35",
            border = "white", lwd = 0.8)
    text(x, y, n, col = "white", cex = 0.56 * cex, font = 2)
  }

  text(2, 99, "d)  the path it describes", adj = c(0, 0.5), font = 2,
       cex = 0.82 * cex)

  ## --- the stem, drawn before the leaf so the leaf sits on top of it -----
  xs <- 46                                # the plant's axis
  polygon(c(xs - 1.5, xs + 1.5, xs + 2.8, xs - 2.8), c(90, 90, 62, 62),
          col = .fig_col$plant, border = "grey55")
  .fig_resistor(xs, 84, xs, 66, n = 8, amp = 2.7, frac = 0.86,
                col = "black", lwd = 1.3)
  text(xs - 8, 78, expression(k(psi) == k[max] * f(psi)), adj = c(1, 0.5),
       cex = 0.64 * cex)
  text(xs - 8, 73.5, "Weibull", adj = c(1, 0.5), cex = 0.56 * cex, font = 3,
       col = "grey35")
  marker(xs + 9, 78, "2")

  ## --- the leaf ---------------------------------------------------------
  th <- seq(0, 2 * pi, length.out = 120)
  polygon(xs + 13 * cos(th), 91 + 13 * 0.40 * k * sin(th), col = "#e8f0e6",
          border = .fig_col$benefit, lwd = 1.2)
  text(xs, 91, expression(psi[stem]), cex = 0.70 * cex)

  ## the fluxes at the leaf surface: water out, CO2 in, radiation in
  .fig_arrow(xs + 7, 93.4, xs + 14, 98, col = .fig_col$supply, lty = 2)
  text(xs + 15, 98, "E", adj = c(0, 0.5), cex = 0.68 * cex,
       col = .fig_col$supply, font = 2)
  .fig_arrow(xs - 1, 98, xs - 5, 93.4, col = .fig_col$benefit)
  text(xs - 0.2, 98, expression(CO[2]), adj = c(0, 0.5), cex = 0.62 * cex,
       col = .fig_col$benefit)
  .fig_arrow(xs - 26, 94.5, xs - 14, 92, col = "grey35")
  text(xs - 27, 94.5, expression(R[n]), adj = c(1, 0.5), cex = 0.64 * cex,
       col = "grey20")
  text(xs + 16, 93, expression(italic(D[leaf])), adj = c(0, 0.5),
       cex = 0.60 * cex, col = "grey35")
  marker(xs + 20, 88, "1")

  ## --- the collar and the root network ----------------------------------
  segments(xs - 6, 62, xs + 6, 62, lwd = 1.5, col = "grey35")
  text(xs - 8, 62, expression(psi[collar]), adj = c(1, 0.5), cex = 0.68 * cex)

  ## the ground, and three soil layers below it. The roots reach out on both
  ## sides so the block reads as one soil column rather than a spray.
  soil_edge <- c(58, 51.5, 45, 38.5)
  soil_fill <- c("#f5f2ec", "#efe9df", "#e7dfd1")
  for (i in 1:3) {
    rect(6, soil_edge[i + 1], 94, soil_edge[i], col = soil_fill[i],
         border = NA)
  }
  segments(6, 58, 94, 58, lwd = 1.3, col = "grey45")
  for (yy in soil_edge[-1]) {
    segments(6, yy, 94, yy, lty = 3, col = "grey70")
  }
  layer_y <- c(54.7, 48.2, 41.7)          # mid-depth of each layer
  root_x <- c(21, 46, 71)                 # where each root reaches
  for (i in seq_along(layer_y)) {
    ## the root out to the layer, then its resistance down into the soil
    lines(c(xs, root_x[i]), c(62, layer_y[i] + 1.6), col = "grey45", lwd = 1.1)
    .fig_resistor(root_x[i], layer_y[i] + 1.6, root_x[i], layer_y[i] - 1.6,
                  n = 3, amp = 1.8, frac = 0.85, col = "black", lwd = 1.1)
    text(root_x[i] + 2.4, layer_y[i], bquote(r[R * .(i)]), adj = c(0, 0.5),
         cex = 0.56 * cex)
    text(root_x[i] - 2.4, layer_y[i], bquote(psi[s * .(i)]), adj = c(1, 0.5),
         cex = 0.58 * cex, col = "grey20")
    text(8, layer_y[i], bquote(z[.(i)]), adj = c(0, 0.5), cex = 0.56 * cex,
         col = "grey45")
  }
  text(94, 35, "depth z, and a gravity head with it", adj = c(1, 0.5),
       cex = 0.56 * cex, font = 3, col = "grey40")
  marker(86, 43, "3")

  ## --- the notes ---------------------------------------------------------
  ##
  ## Numbered rather than set beside the drawing: the panel is portrait, so a
  ## right-hand annotation column would leave the plant two inches wide and six
  ## tall, which is a telegraph pole rather than a tree.
  notes <- list(
    c("The leaf.  Water out and CO2 in through one conductance.",
      "With the energy balance on, T_leaf and D_leaf become",
      "outputs of E rather than drivers of it."),
    c("The stem.  Conductivity falls as tension rises, so E",
      "saturates -- and that saturation is what makes an interior",
      "optimum exist. It is also the tension the cost curve prices:",
      "the three models accept everything else and differ here."),
    c("Supply topology.  Either one psi_soil in series, or this",
      "multi-layer network solved at the collar.  Different",
      "models, not one restricted to the other.")
  )
  y <- 30
  for (i in seq_along(notes)) {
    marker(4.5, y, as.character(i))
    text(9.5, y - (seq_along(notes[[i]]) - 1) * 2.3, notes[[i]],
         adj = c(0, 0.5), cex = 0.58 * cex, col = "grey20")
    y <- y - (length(notes[[i]]) * 2.3 + 1.9)
  }

  invisible(NULL)
}

## ---------------------------------------------------------------------------

##' Overview figure: how the parts of the model fit together
##'
##' Draws the four zones of the model in one figure: the profit panel that
##' everything else serves, the supply-equals-demand root-find inside the
##' benefit, the causal wiring between the parts, and the soil-to-leaf path the
##' equations describe.
##'
##' The two left-hand panels are computed from `model`, so they move if the model
##' does. The two right-hand zones are schematic.
##'
##' The layout follows figure 1 of Westoby et al. (2012), which pairs a data
##' panel with a box-and-arrow flow and a physical cartoon so that the reader can
##' see the same object three ways.
##'
##' @param model A `Leaf`, already driven. Defaults to a leaf at
##'   `psi_soil = 1` MPa and `PPFD = 1500`, which is the reference state used
##'   throughout `vignette("the-models")`.
##' @param cex Overall text scaling, applied on top of the per-element sizes.
##'   Raise it for a slide, lower it for a two-column figure.
##'
##' @return `NULL`, invisibly. Called for the plot.
##'
##' @examples
##' plot_model_overview()
##'
##' # a vector figure for a manuscript
##' \dontrun{
##' pdf("overview.pdf", width = 11, height = 6.6)
##' plot_model_overview()
##' dev.off()
##' }
##'
##' @references
##' Westoby, M., Falster, D.S. and Moles, A.T. (2012) Trade-offs between height
##' growth rate, stem persistence and maximum height among plant species in a
##' post-fire succession. *Journal of Theoretical Biology*.
##'
##' @importFrom graphics arrows layout lines mtext par plot.new plot.window
##' @importFrom graphics points polygon rect segments strwidth text
##' @export
plot_model_overview <- function(model = NULL, cex = 1) {
  if (is.null(model)) {
    model <- leaf_model()
    set_drivers(model, psi_soil = 1.0, PPFD = 1500)
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  layout(matrix(c(1, 3, 4,
                  2, 3, 4), nrow = 2, byrow = TRUE),
         widths = c(1.00, 1.22, 1.16))

  .fig_profit_panel(model, cex = cex)
  .fig_aci_panel(model, cex = cex)
  .fig_flow_panel(cex = cex)
  .fig_plant_panel(cex = cex)

  invisible(NULL)
}
