# The model registry, inst/model-registry.yml (#166), and the three rules its
# header states. Rule 2 lives here: every symbol and field the registry names
# must exist, and every bound field must appear in it. The Doxygen half of the
# same rule (every `code` resolves in the XML; what each equation reads is what
# its function references; every literal is written in its function) runs in
# .github/workflows/docs.yml through tools/registry_symbols.R, and here too
# when a local render exists.

skip_if_not_installed("yaml")

registry <- model_registry()
leaf <- leaf_model()

# Registry parameters held by a member no R binding exposes: four constants
# (static or constexpr members, valued through Doxygen's initializer) and two
# constructor arguments the Leaf stores but does not bind (leaf_control()
# holds their defaults).
.registry_unbound_members <- c("damage_surv_floor", "collar_root_tol", "k_crit_fraction",
                               "GSS_tol_abs", "vulnerability_curve_ncontrol",
                               "medlyn_H2O_CO2_ratio")

# Inputs with no default by design: an unset one reads NaN and the code refuses
# it where it is used (the three model-owned parameters, the energy balance's
# near-infrared), or it is seated only by a supply network.
.registry_no_default <- c("NIR_int_", "CF77_lambda_", "TF24_floor_lambda_o",
                          "LeastCost_beta", "single_resistance")

test_that("the registry parses and validates", {
  expect_s3_class(registry, "model_registry")
  expect_gt(nrow(registry$parameters), 0)
  expect_gt(nrow(registry$equations), 0)
  expect_output(print(registry), "parameters")
})

test_that("a registry that names what does not exist is refused", {
  broken <- registry$raw
  broken$equations$arrhenius$reads <- c("gas_constant", "no_such_parameter")
  path <- tempfile(fileext = ".yml")
  on.exit(unlink(path), add = TRUE)
  yaml::write_yaml(broken, path)
  expect_error(registry_read(path), "no_such_parameter")
  broken <- registry$raw
  broken$equations$arrhenius$kind <- "approximate"
  yaml::write_yaml(broken, path)
  expect_error(registry_read(path), "unknown kind")
  broken <- registry$raw
  broken$parameters$vcmax_ha_$cite <- "medlyn2002"
  yaml::write_yaml(broken, path)
  expect_error(registry_read(path), "Zotero key")
})

test_that("every symbol the registry names is in the headers, and every field it names is bound", {
  # A textual check on the installed headers; Doxygen's is the authoritative
  # one and runs in docs.yml. This catches a rename with no Doxygen present.
  inc <- system.file("include", "phylloptim", package = "phylloptim")
  src <- unlist(lapply(list.files(inc, pattern = "\\.hpp$", full.names = TRUE),
                       readLines, warn = FALSE))
  symbols <- c(registry$parameters$code, registry$parameters$default_code,
               registry$equations$code, unlist(registry$equations$also),
               unlist(registry$parameters$aliases))
  symbols <- unique(symbols[!is.na(symbols)])
  names_only <- sub("^.*::", "", symbols)
  found <- vapply(names_only, function(n) any(grepl(paste0("\\b", n, "\\b"), src, perl = TRUE)), TRUE)
  expect_true(all(found), info = paste("not in the headers:", paste(symbols[!found], collapse = ", ")))

  # A literal's code names the function it is written in, not a field.
  p <- registry$parameters
  fields <- p$field[is.na(p$literal) & !is.na(p$field)]
  unbound <- sort(setdiff(fields, names(leaf)))
  # Members the registry names that no binding exposes. Their values reach a
  # document through Doxygen's initializer (the constants) or not at all (a
  # constructor argument the leaf keeps no field for); either way the registry
  # records them, and this list says which they are.
  expect_equal(unbound, sort(.registry_unbound_members))

  writes <- unique(unlist(registry$equations$writes))
  expect_true(all(writes %in% names(leaf)),
              info = paste("unbound writes:", paste(setdiff(writes, names(leaf)), collapse = ", ")))
  drivers <- registry$parameters$driver[!is.na(registry$parameters$driver)]
  expect_true(all(drivers %in% names(formals(set_drivers))))
})

# Rule 2: every bound field appears in the registry, as a parameter, as what an
# equation writes, or as a diagnostic (the solver's bookkeeping and the fields
# nothing writes). A field bound later fails here, by name, until it is
# registered.
test_that("every bound field of the Leaf is accounted for", {
  classes <- yaml::read_yaml(system.file("RcppR6_classes.yml", package = "phylloptim"))
  active <- classes$Leaf$active
  bound <- names(active)[vapply(active, function(a) identical(a$access, "field"), TRUE)]
  expect_gt(length(bound), 100)
  missing <- sort(setdiff(bound, .registry_covered(registry)))
  expect_equal(missing, character(), info = paste("unregistered:", paste(missing, collapse = ", ")))
  # And nothing the registry lists as a diagnostic is a parameter or a write.
  expect_length(intersect(registry$diagnostics$id,
                          c(registry$parameters$field, unlist(registry$equations$writes))), 0)
})

test_that("registry_values() reads every parameter live, and the copied defaults equal their constants", {
  v <- registry_values(registry, leaf)
  expect_true(all(v$value_source[!is.na(v$literal)] == "literal"))
  expect_true(all(v$value_source[!is.na(v$driver)] == "set_drivers()"))
  # Every member the leaf binds has a value, except the driver that has no
  # default by design.
  members <- !is.na(v$field) & is.na(v$driver) & is.na(v$literal) &
    !v$field %in% .registry_unbound_members
  expect_true(all(v$value_source[members] %in% c("binding", "binding (vector)")))
  scalar <- members & v$value_source == "binding"
  expect_true(all(is.finite(v$value[scalar & !v$id %in% .registry_no_default])))
  expect_true(all(is.na(v$value[v$id %in% .registry_no_default])))
  # A constexpr has no binding: its value comes from Doxygen's initializer,
  # which this leaf-only call does not have.
  expect_true(all(is.na(v$value_source[grepl("^phylloptim::[a-zA-Z_]+$", v$code) & is.na(v$literal)])))

  # The members whose defaults are copied from constants.hpp equal them. Read
  # the header's `inline constexpr double name = value;` lines, so a constant
  # moved in one place and not the other surfaces.
  hpp <- readLines(system.file("include", "phylloptim", "constants.hpp", package = "phylloptim"), warn = FALSE)
  m <- regmatches(hpp, regexec("^inline constexpr double ([A-Za-z_0-9]+) = ([-0-9.e]+)", hpp))
  consts <- stats::setNames(as.numeric(vapply(Filter(length, m), `[`, "", 3)),
                            paste0("phylloptim::", vapply(Filter(length, m), `[`, "", 2)))
  # Namespace constants only: a default copied from a class member (the knot
  # count) has no constants.hpp line to compare with.
  copied <- v[grepl("^phylloptim::[A-Za-z_0-9]+$", v$default_code), ]
  expect_true(all(copied$default_code %in% names(consts)))
  expect_equal(unname(consts[copied$default_code]), copied$value)
})

test_that("R/registry-locations.R stands alone, as docs.yml sources it", {
  skip_if_not_installed("codetools")
  path <- testthat::test_path("..", "..", "R", "registry-locations.R")
  skip_if_not(file.exists(path), "the source tree is not reachable from here")
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  defined <- ls(env, all.names = TRUE)
  fns <- Filter(is.function, mget(defined, envir = env))
  globals <- unique(unlist(lapply(fns, codetools::findGlobals)))
  leaked <- intersect(globals, setdiff(ls(asNamespace("phylloptim"), all.names = TRUE), defined))
  expect_length(leaked, 0)
})

test_that("the registry resolves against a local Doxygen render", {
  skip_if_not_installed("xml2")
  root <- testthat::test_path("..", "..")
  xml <- file.path(root, "docs", "xml")
  skip_if_not(file.exists(file.path(xml, "index.xml")),
              "no Doxygen XML in docs/xml; run `doxygen` at the package root")
  loc <- registry_locations(xml, include_dir = file.path(root, "inst", "include"),
                            filter = file.path(root, "tools", "doxygen_filter.awk"))
  res <- registry_resolve(registry, locations = loc, unreferenced = "error")
  expect_length(res$unresolved, 0)
  expect_true(all(!is.na(res$parameters$line)))
  expect_true(all(is.na(res$equations$code) | !is.na(res$equations$body_start)))
  # The mapped line holds the symbol it claims to.
  hdr <- function(f) readLines(file.path(root, "inst", "include", f), warn = FALSE)
  for (i in seq_len(nrow(res$parameters))) {
    p <- res$parameters[i, ]
    expect_true(grepl(sub("^.*::", "", p$code), hdr(p$file)[p$line], fixed = TRUE),
                label = paste(p$id, "at", p$file, p$line))
  }
  # Every literal is a number written in the function its entry names.
  expect_equal(.registry_literals_unwritten(res, file.path(root, "inst", "include")),
               character())
  # With the XML, the constants have values too, and the copied defaults match.
  v <- registry_values(registry, leaf, locations = loc)
  expect_true(all(!is.na(v$value[v$class == "constant"])))
  consts <- .registry_initializer_values(loc)
  copied <- v[grepl("^phylloptim::[A-Za-z_0-9]+$", v$default_code), ]
  expect_equal(unname(consts[copied$default_code]), copied$value)
})
