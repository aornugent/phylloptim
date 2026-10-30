# The model registry's live half: the registry as data frames and the
# parameters' values read off the package. The pure half (reading the YAML,
# Doxygen locations) is R/registry-locations.R.

#' The model registry
#'
#' phylloptim's leaf as one record: every equation (its LaTeX, the C++ symbol
#' that evaluates it, the parameters it reads, the bound fields it writes) and
#' every parameter (symbol, units, class, meaning, source, placeholder flag, the
#' symbol that holds it). The file is
#' `inst/model-registry.yml`; documents are generated from it rather than
#' written beside it, and CI checks it (`tests/testthat/test-registry.R`,
#' `docs.yml`). Its header states the three rules it is kept by: no stored line
#' numbers or defaults, validated in CI, one registry per package.
#'
#' @param path the YAML file; the package's own by default
#' @return an object of class `model_registry`: `package`, `sections` (a data
#'   frame of `id` and `title`), `parameters` and `equations` (data frames, one
#'   row each, with list-columns for `reads`, `writes`, `also` and `aliases`),
#'   `diagnostics` (the bound fields that are neither), `raw` (the
#'   validated list) and `path`
#' @seealso [registry_values()] for the values at render time,
#'   [registry_resolve()] for file, line and API link from Doxygen's XML
#' @examples
#' reg <- model_registry()
#' reg
#' subset(reg$equations, section == "respiration", c(id, title, code))
#' @export
model_registry <- function(path = system.file("model-registry.yml", package = "phylloptim")) {
  raw <- registry_read(path)
  tabs <- registry_tables(raw)
  structure(list(package = raw$package,
                 sections = data.frame(id = names(raw$sections),
                                       title = unname(unlist(raw$sections)),
                                       stringsAsFactors = FALSE),
                 parameters = tabs$parameters,
                 equations = tabs$equations,
                 diagnostics = tabs$diagnostics,
                 raw = raw, path = path),
            class = "model_registry")
}

#' @export
print.model_registry <- function(x, ...) {
  cat("<model_registry> ", x$package, ": ", nrow(x$parameters), " parameters, ",
      nrow(x$equations), " equations\n", sep = "")
  for (i in seq_len(nrow(x$sections))) {
    s <- x$sections$id[i]
    cat(sprintf("  %-16s %3d parameters  %3d equations\n", s,
                sum(x$parameters$section == s), sum(x$equations$section == s)))
  }
  invisible(x)
}

# The bound fields the registry accounts for: parameters' fields, what equations
# write, and the diagnostics.
.registry_covered <- function(registry) {
  unique(c(registry$parameters$field[!is.na(registry$parameters$field)],
           unlist(registry$equations$writes), registry$diagnostics$id))
}

# The constexprs' values from Doxygen's initializers, which may name one another
# (`kg_to_mol_h2o = 1.0 / molar_mass_h2o`): evaluate in passes until nothing new
# resolves. Returns a named numeric vector keyed by qualified symbol.
.registry_initializer_values <- function(locations) {
  if (is.null(locations)) return(numeric())
  vars <- locations[locations$kind == "variable" & !is.na(locations$initializer), , drop = FALSE]
  env <- new.env(parent = baseenv())
  out <- numeric()
  pending <- seq_len(nrow(vars))
  for (pass in 1:4) {
    if (!length(pending)) break
    still <- integer()
    for (i in pending) {
      expr <- gsub("phylloptim::", "", vars$initializer[i], fixed = TRUE)
      v <- tryCatch(eval(parse(text = expr), env), error = function(e) NULL)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
        out[vars$symbol[i]] <- v
        assign(sub("^.*::", "", vars$symbol[i]), v, envir = env)
      } else {
        still <- c(still, i)
      }
    }
    pending <- still
  }
  out
}

#' The parameters' values, read at render time
#'
#' Adds `value` and `value_source` to the registry's parameters: a bound field
#' is read off `leaf`, a driver off [set_drivers()]'s default, a constexpr off
#' its initializer in Doxygen's XML when `locations` is given, and a literal
#' off the registry (the one kind of number it stores). Anything else is `NA`.
#' Nothing is read from a comment or a document, which is rule 1 of the
#' registry's header.
#'
#' @param registry a [model_registry()]
#' @param leaf the `Leaf` to read bound fields off; `NULL` reads none
#' @param locations the output of [registry_locations()], for the constexprs
#' @return the registry's `parameters` data frame with `value` and
#'   `value_source` (`"binding"`, `"binding (vector)"`, `"set_drivers()"`,
#'   `"initializer"`, `"literal"` or `NA`), and the list-column `value_vector`
#'   holding a per-layer field whole (its `value` is then `NA`)
#' @examples
#' v <- registry_values(model_registry())
#' v[v$section == "respiration", c("id", "symbol", "value", "value_source")]
#' @export
registry_values <- function(registry, leaf = leaf_model(), locations = NULL) {
  p <- registry$parameters
  consts <- .registry_initializer_values(locations)
  drivers <- formals(set_drivers)
  bound <- if (is.null(leaf)) character() else names(leaf)
  value <- rep(NA_real_, nrow(p))
  vectors <- vector("list", nrow(p))
  src <- rep(NA_character_, nrow(p))
  for (i in seq_len(nrow(p))) {
    if (!is.na(p$literal[i])) {
      value[i] <- p$literal[i]; src[i] <- "literal"
    } else if (!is.na(p$driver[i]) && p$driver[i] %in% names(drivers)) {
      value[i] <- eval(drivers[[p$driver[i]]]); src[i] <- "set_drivers()"
    } else if (!is.na(f <- p$field[i]) && f %in% bound) {
      v <- leaf[[f]]
      if (is.logical(v)) v <- as.numeric(v)
      # A per-layer field (psi_soil_, r_R_H_min, ...) is a vector, and empty
      # until a solve seats it: its scalar value is NA and the whole is kept.
      if (length(v) == 1L) value[i] <- v else vectors[[i]] <- v
      src[i] <- if (length(v) == 1L) "binding" else "binding (vector)"
    } else if (p$code[i] %in% names(consts)) {
      value[i] <- consts[[p$code[i]]]; src[i] <- "initializer"
    }
  }
  p$value <- value
  p$value_source <- src
  p$value_vector <- vectors
  p
}
