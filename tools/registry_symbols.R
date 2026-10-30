#!/usr/bin/env Rscript
# Resolve the model registry against Doxygen's XML, and fail if it cannot.
#
# Run by .github/workflows/docs.yml after the render, where the package is not
# installed: only R/registry-locations.R is sourced, so that file must stay
# free of anything else in the package (test-registry.R asserts it).
#
#   Rscript tools/registry_symbols.R docs/xml docs/html/model-registry.yml
#
# Four things are checked. Every `code` (and `also`) symbol the registry names
# resolves to a member in the XML, so a rename in the headers that the registry
# did not follow fails here. Every resolved line maps back to the source
# through the filter's line map, which refuses a header whose code lines the
# filter altered. Every parameter an equation `reads` is referenced by the
# equation's function (Doxygen's REFERENCES_RELATION), so an equation that
# quietly stopped reading a parameter, or started reading one not registered,
# fails here. And every `literal:` is a number written in the body of the
# function its entry names, so a constant changed on one side only fails here. The second argument, if given, receives the resolved registry:
# every parameter and equation with its file, line, body span, API anchor and
# initializer, for a document that wants to link without a Doxygen run of its
# own.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: registry_symbols.R <doxygen xml dir> [<output yml>]")
xml_dir <- args[1]
out_path <- if (length(args) >= 2L) args[2] else NA_character_

root <- normalizePath(file.path(dirname(xml_dir), ".."))
source(file.path(root, "R", "registry-locations.R"))

raw <- registry_read(file.path(root, "inst", "model-registry.yml"))
locations <- registry_locations(xml_dir,
                                include_dir = file.path(root, "inst", "include"),
                                filter = file.path(root, "tools", "doxygen_filter.awk"))
res <- registry_resolve(raw, locations = locations, unreferenced = "none")
# What Doxygen records as a function's references differs between versions (CI
# runs 1.9.8), so on a failure print what this one recorded for each function
# involved: the next fix is then a matter of reading the log, not of guessing.
if (nrow(res$unreferenced)) {
  eq <- res$equations
  for (i in seq_len(nrow(res$unreferenced))) {
    e <- res$unreferenced$equation[i]
    p <- res$unreferenced$parameter[i]
    fns <- c(eq$code[eq$id == e], eq$also[[which(eq$id == e)]])
    target <- res$parameters$code[res$parameters$id == p]
    cat(sprintf("\n%s reads %s (%s), which none of these references:\n", e, p, target))
    for (f in fns) {
      refs <- unique(unlist(locations$references[locations$symbol == f]))
      cat("  ", f, ": ", paste(refs, collapse = " "), "\n", sep = "")
    }
  }
  stop("reads not referenced by the equation's function: ",
       paste(res$unreferenced$equation, res$unreferenced$parameter, sep = ":", collapse = ", "),
       call. = FALSE)
}

unwritten <- .registry_literals_unwritten(res, file.path(root, "inst", "include"))
if (length(unwritten)) {
  stop("literals not written in the function their entry names: ",
       paste(unwritten, collapse = ", "), call. = FALSE)
}

cat(sprintf("model registry: %d parameters and %d equations resolved against %d symbols\n",
            nrow(res$parameters), nrow(res$equations), nrow(locations)))

if (!is.na(out_path)) {
  strip <- function(df) {
    lapply(seq_len(nrow(df)), function(i) {
      row <- lapply(df[i, , drop = FALSE], function(col) {
        v <- if (is.list(col)) col[[1]] else col
        if (is.factor(v)) as.character(v) else v
      })
      row <- row[!vapply(row, function(v) length(v) == 1L && is.na(v), TRUE)]
      names(row) <- names(df)[names(df) %in% names(row)]
      row
    })
  }
  out <- list(package = raw$package,
              generated_by = "tools/registry_symbols.R",
              parameters = stats::setNames(strip(res$parameters), res$parameters$id),
              equations = stats::setNames(strip(res$equations), res$equations$id))
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  yaml::write_yaml(out, out_path)
  cat("wrote", out_path, "\n")
}
