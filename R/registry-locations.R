# The model registry's pure half: read and validate inst/model-registry.yml,
# and resolve the symbols it names against Doxygen's XML.
#
# ⚠️ EVERYTHING IN THIS FILE IS SOURCED STANDALONE BY tools/registry_symbols.R
# in the docs workflow, where the package is not installed (Doxygen runs on a
# runner with no compiler step). So nothing here may call another function of
# this package, read a Leaf, or use a package other than base, yaml and xml2;
# tests/testthat/test-registry.R asserts that with codetools::findGlobals().
# The half that needs a live Leaf is R/registry.R.
#
# Why Doxygen rather than a stored line: a stored line is wrong the day after
# the file above it changes, and nothing can tell. Doxygen's XML carries, for
# every member, the line of its declaration and the span of its body, in the
# stream it parsed. That stream is tools/doxygen_filter.awk's output, which
# inserts lines (a `\file` block, `\verbatim` pairs) and drops held blank
# comment lines, so its line numbers are not the source's. What is identical in
# both streams, byte for byte, is the sequence of lines that are not whole-line
# `//` comments, and docs.yml asserts that on every render. Every line Doxygen
# reports is one of those code lines, so the k-th code line of the filtered
# stream is the k-th code line of the source: registry_line_map() builds that
# map and checks the two sequences are identical before trusting it.

.registry_classes <- c("trait", "hyperparameter", "constant", "driver")
.registry_kinds <- c("closed_form", "identity", "derivative", "integrator")
.registry_cite_pattern <- "^[A-Za-z][A-Za-z-]*-[0-9]{4}[a-z]?$"

#' Read the model registry
#'
#' Reads `inst/model-registry.yml` (or another registry in the same schema) and
#' validates it: every section, class and kind is one of the allowed set, every
#' `reads` names a parameter, and every numeric is finite. The schema is written
#' at the top of the file.
#'
#' This is the raw list. [model_registry()] is the same registry as a pair of
#' data frames, which is what documents want.
#'
#' @param path the YAML file
#' @return the validated list, with `literal` coerced to numeric
#' @seealso [model_registry()], [registry_resolve()]
#' @keywords internal
registry_read <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("registry_read() needs the yaml package", call. = FALSE)
  }
  raw <- yaml::read_yaml(path, eval.expr = FALSE)
  .registry_validate(raw, path)
}

# A number written in the YAML. YAML 1.1 reads `1e-3` as a string, so coerce,
# and refuse anything that does not come out finite.
.registry_number <- function(x, what) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || !is.finite(v)) {
    stop("model registry: ", what, " is not a finite number: ",
         deparse(x), call. = FALSE)
  }
  v
}

.registry_character <- function(x) {
  if (is.null(x)) character() else as.character(unlist(x))
}

.registry_validate <- function(raw, path) {
  fail <- function(...) stop("model registry (", path, "): ", ..., call. = FALSE)
  need <- function(x, keys, what) {
    missing <- setdiff(keys, names(x))
    if (length(missing)) {
      fail(what, " lacks ", paste(missing, collapse = ", "))
    }
  }
  need(raw, c("package", "sections", "parameters", "equations"), "the registry")
  sections <- names(raw$sections)
  pids <- names(raw$parameters)
  eids <- names(raw$equations)
  if (anyDuplicated(pids)) fail("duplicate parameter ids")
  if (anyDuplicated(eids)) fail("duplicate equation ids")
  if (length(intersect(pids, eids))) {
    fail("ids shared by a parameter and an equation: ",
         paste(intersect(pids, eids), collapse = ", "))
  }
  for (id in pids) {
    p <- raw$parameters[[id]]
    what <- paste0("parameter `", id, "`")
    need(p, c("section", "symbol", "units", "class", "meaning", "source",
              "placeholder", "code"), what)
    if (!p$section %in% sections) fail(what, ": unknown section ", p$section)
    if (!p$class %in% .registry_classes) fail(what, ": unknown class ", p$class)
    if (!is.logical(p$placeholder)) fail(what, ": placeholder must be true or false")
    if (!is.character(p$code) || !nzchar(p$code)) fail(what, ": code must name a symbol")
    if (!is.null(p$cite) && !grepl(.registry_cite_pattern, p$cite)) {
      fail(what, ": cite is not a Zotero key (Surname-Year): ", p$cite)
    }
    if (!is.null(p$literal)) {
      raw$parameters[[id]]$literal <- .registry_number(p$literal, paste(what, "literal"))
    }
  }
  for (id in names(raw$diagnostics)) {
    d <- raw$diagnostics[[id]]
    need(d, c("meaning", "code"), paste0("diagnostic `", id, "`"))
  }
  for (id in eids) {
    e <- raw$equations[[id]]
    what <- paste0("equation `", id, "`")
    need(e, c("section", "title", "latex", "reads", "writes", "kind"), what)
    if (!"code" %in% names(e)) fail(what, ": needs code (a symbol, or ~ for none)")
    if (!e$section %in% sections) fail(what, ": unknown section ", e$section)
    if (!e$kind %in% .registry_kinds) fail(what, ": unknown kind ", e$kind)
    unknown <- setdiff(.registry_character(e$reads), pids)
    if (length(unknown)) fail(what, ": reads unknown parameters ", paste(unknown, collapse = ", "))
    if (!is.null(e$cite) && !grepl(.registry_cite_pattern, e$cite)) {
      fail(what, ": cite is not a Zotero key (Surname-Year): ", e$cite)
    }
  }
  raw
}

# The registry as two data frames, one row per parameter and per equation. The
# list-valued fields (reads, writes, also, aliases) are list-columns.
registry_tables <- function(raw) {
  opt <- function(x, key, default = NA_character_) {
    v <- x[[key]]
    if (is.null(v)) default else v
  }
  pids <- names(raw$parameters)
  parameters <- data.frame(
    id = pids,
    section = vapply(raw$parameters, function(p) p$section, ""),
    symbol = vapply(raw$parameters, function(p) p$symbol, ""),
    units = vapply(raw$parameters, function(p) p$units, ""),
    class = vapply(raw$parameters, function(p) p$class, ""),
    meaning = vapply(raw$parameters, function(p) p$meaning, ""),
    source = vapply(raw$parameters, function(p) p$source, ""),
    cite = vapply(raw$parameters, function(p) opt(p, "cite"), ""),
    placeholder = vapply(raw$parameters, function(p) isTRUE(p$placeholder), TRUE),
    code = vapply(raw$parameters, function(p) p$code, ""),
    default_code = vapply(raw$parameters, function(p) opt(p, "default_code"), ""),
    driver = vapply(raw$parameters, function(p) opt(p, "driver"), ""),
    field = vapply(raw$parameters, function(p) {
      if (!is.null(p$field)) return(p$field)
      if (grepl("^phylloptim::Leaf::", p$code)) sub("^phylloptim::Leaf::", "", p$code) else NA_character_
    }, ""),
    literal = vapply(raw$parameters, function(p) opt(p, "literal", NA_real_), 0),
    read_by = vapply(raw$parameters, function(p) opt(p, "read_by"), ""),
    stringsAsFactors = FALSE, row.names = NULL)
  parameters$aliases <- lapply(raw$parameters, function(p) .registry_character(p$aliases))
  eids <- names(raw$equations)
  equations <- data.frame(
    id = eids,
    section = vapply(raw$equations, function(e) e$section, ""),
    title = vapply(raw$equations, function(e) e$title, ""),
    latex = vapply(raw$equations, function(e) trimws(e$latex, "right"), ""),
    code = vapply(raw$equations, function(e) opt(e, "code"), ""),
    kind = vapply(raw$equations, function(e) e$kind, ""),
    cite = vapply(raw$equations, function(e) opt(e, "cite"), ""),
    note = vapply(raw$equations, function(e) opt(e, "note"), ""),
    stringsAsFactors = FALSE, row.names = NULL)
  for (key in c("reads", "writes", "also")) {
    equations[[key]] <- lapply(raw$equations, function(e) .registry_character(e[[key]]))
  }
  dids <- names(raw$diagnostics)
  diagnostics <- data.frame(
    id = as.character(dids),
    meaning = vapply(raw$diagnostics, function(d) d$meaning, character(1)),
    code = vapply(raw$diagnostics, function(d) d$code, character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  list(parameters = parameters, equations = equations, diagnostics = diagnostics)
}

#' Map the lines of a header's filtered stream back to its source
#'
#' Doxygen reads every header through `tools/doxygen_filter.awk`, which inserts
#' and drops comment lines, so the lines in its XML are the filtered stream's.
#' Whole-line `//` comments are the only lines the filter touches: the sequence
#' of all other lines is byte-identical in both streams (docs.yml asserts it),
#' and every line Doxygen reports is one of them. So the k-th such line of the
#' filtered stream is the k-th of the source, and that is the map.
#'
#' @param header path of the header
#' @param filter path of the awk filter
#' @return an integer vector as long as the filtered stream: for each filtered
#'   line, the source line it is, or `NA` for a comment line
#' @keywords internal
registry_line_map <- function(header, filter) {
  src <- readLines(header, warn = FALSE)
  filtered <- system2("awk", c("-f", shQuote(filter), shQuote(header)), stdout = TRUE)
  src_code <- which(!grepl("^[ \t]*//", src))
  # The filtered stream's code lines: outside the `/*! \file` ... `*/` block the
  # filter emits for the first comment, and not a `//` line. The same predicate
  # as docs.yml's byte-identity check.
  in_file <- FALSE
  is_code <- logical(length(filtered))
  for (i in seq_along(filtered)) {
    line <- filtered[i]
    if (!in_file && grepl("^[ \t]*/\\*! \\\\file", line)) {
      in_file <- TRUE
    } else if (in_file) {
      if (grepl("^[ \t]*\\*/[ \t]*$", line)) in_file <- FALSE
    } else if (!grepl("^[ \t]*//", line)) {
      is_code[i] <- TRUE
    }
  }
  f_code <- which(is_code)
  if (length(f_code) != length(src_code) ||
      !identical(filtered[f_code], src[src_code])) {
    stop("registry_line_map(): the filter changed a non-comment line of ",
         header, "; the line map cannot be trusted", call. = FALSE)
  }
  map <- rep(NA_integer_, length(filtered))
  map[f_code] <- src_code
  map
}

#' Where every symbol is, from Doxygen's XML
#'
#' Walks the compounds Doxygen's `index.xml` lists and returns one row per
#' member: its qualified name, kind, file and line (declaration), body span,
#' the anchor of its page in the HTML render, its initializer where it has one
#' (a constexpr's value), and the symbols it references. Lines are the source's,
#' through [registry_line_map()].
#'
#' Needs a render with `GENERATE_XML = YES` and `REFERENCES_RELATION = YES`,
#' which the package's `Doxyfile` sets: `doxygen` at the package root writes
#' `docs/xml`.
#'
#' @param xml_dir Doxygen's XML output directory
#' @param include_dir the directory the XML's file paths are relative to
#'   (`STRIP_FROM_PATH` in the Doxyfile): the package's `inst/include`
#' @param filter the awk filter Doxygen read the headers through
#' @return a data frame with columns `symbol`, `kind`, `id`, `compound`, `file`,
#'   `line`, `body_file`, `body_start`, `body_end`, `url`, `initializer` and the
#'   list-column `references`; the body is in `body_file`, which differs from
#'   `file` for a member defined out of line
#' @export
registry_locations <- function(xml_dir,
                               include_dir = file.path(dirname(dirname(xml_dir)), "inst", "include"),
                               filter = file.path(dirname(dirname(xml_dir)), "tools", "doxygen_filter.awk")) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("registry_locations() needs the xml2 package", call. = FALSE)
  }
  index_path <- file.path(xml_dir, "index.xml")
  if (!file.exists(index_path)) {
    stop("no Doxygen XML at ", xml_dir, " (render with GENERATE_XML = YES)", call. = FALSE)
  }
  index <- xml2::read_xml(index_path)
  compounds <- xml2::xml_find_all(index, "//compound")
  kinds <- xml2::xml_attr(compounds, "kind")
  refids <- xml2::xml_attr(compounds, "refid")
  keep <- kinds %in% c("class", "struct", "namespace")
  text_of <- function(node, xpath) {
    hit <- xml2::xml_find_first(node, xpath)
    if (inherits(hit, "xml_missing")) NA_character_ else xml2::xml_text(hit)
  }
  attr_int <- function(node, name) {
    v <- xml2::xml_attr(node, name)
    if (is.na(v)) NA_integer_ else as.integer(v)
  }
  rows <- list()
  enum_values <- list()
  for (refid in refids[keep]) {
    doc <- xml2::read_xml(file.path(xml_dir, paste0(refid, ".xml")))
    cdef <- xml2::xml_find_first(doc, "//compounddef")
    cname <- text_of(cdef, "./compoundname")
    for (m in xml2::xml_find_all(cdef, ".//memberdef")) {
      name <- text_of(m, "./name")
      qual <- text_of(m, "./qualifiedname")
      if (is.na(qual)) qual <- paste(cname, name, sep = "::")
      loc <- xml2::xml_find_first(m, "./location")
      id <- xml2::xml_attr(m, "id")
      anchor <- sub(paste0("^", refid, "_1"), "", id)
      init <- text_of(m, "./initializer")
      refs <- xml2::xml_find_all(m, "./references")
      # An enum's values are children of its memberdef, and a function that
      # reads one (the damage layer indexes its parameter pack by them) refers
      # to the value's id, so they are named here too.
      for (ev in xml2::xml_find_all(m, "./enumvalue")) {
        enum_values[[length(enum_values) + 1L]] <- c(
          id = xml2::xml_attr(ev, "id"), symbol = paste(cname, text_of(ev, "./name"), sep = "::"))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        symbol = qual, kind = xml2::xml_attr(m, "kind"), id = id, compound = refid,
        file = xml2::xml_attr(loc, "file"),
        line = attr_int(loc, "line"),
        body_file = xml2::xml_attr(loc, "bodyfile"),
        body_start = attr_int(loc, "bodystart"),
        body_end = attr_int(loc, "bodyend"),
        url = paste0(refid, ".html#", anchor),
        initializer = if (is.na(init)) NA_character_ else trimws(sub("^=", "", trimws(init))),
        references = I(list(xml2::xml_attr(refs, "refid"))),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  # References come as member ids; name them. One that is not a member of a
  # compound walked here (a std:: function, say) is dropped.
  by_id <- stats::setNames(out$symbol, out$id)
  for (ev in enum_values) by_id[ev[["id"]]] <- ev[["symbol"]]
  out$references <- I(lapply(out$references, function(ids) unname(by_id[ids[ids %in% names(by_id)]])))
  # Filtered lines to source lines, one map per header.
  maps <- list()
  map_for <- function(file) {
    if (is.null(maps[[file]])) {
      maps[[file]] <<- registry_line_map(file.path(include_dir, file), filter)
    }
    maps[[file]]
  }
  remap <- function(file, line) {
    if (is.na(file) || is.na(line) || line < 1L) return(NA_integer_)
    map <- map_for(file)
    if (line > length(map)) return(NA_integer_)
    map[line]
  }
  out$line <- mapply(remap, out$file, out$line)
  out$body_start <- mapply(remap, out$body_file, out$body_start)
  out$body_end <- mapply(remap, out$body_file, out$body_end)
  rownames(out) <- NULL
  out
}

# The one row for a symbol. Overloads and a declaration/definition pair share a
# qualified name; prefer a function with a body, then the one that references
# the most, which is the one that does the work.
.registry_locate <- function(locations, symbol) {
  hits <- locations[locations$symbol == symbol, , drop = FALSE]
  if (nrow(hits) == 0L) return(NULL)
  if (nrow(hits) > 1L) {
    score <- !is.na(hits$body_start) * 1e6 + lengths(hits$references)
    hits <- hits[order(-score), , drop = FALSE]
  }
  hits[1L, , drop = FALSE]
}

#' Resolve the registry's symbols against Doxygen's XML
#'
#' Joins every parameter and equation to the file, line, body span, page anchor
#' and initializer of the symbol its `code` names (and, for equations, of each
#' symbol in `also`). Fails naming every symbol that does not resolve: a rename
#' in the headers that the registry did not follow.
#'
#' It also checks what each equation reads against what its function
#' references, which Doxygen records: a `reads` entry whose parameter symbol
#' (or one of its `aliases`) is referenced by none of the equation's functions
#' is reported, as a warning by default. Literals are exempt, being written
#' inside the function rather than referenced by it.
#'
#' @param raw a registry list from [registry_read()], or a `model_registry`
#' @param xml_dir Doxygen's XML output directory
#' @param locations the output of [registry_locations()], if already in hand
#' @param unreferenced what to do about a read the function does not reference:
#'   `"warn"`, `"error"` or `"none"`
#' @param ... passed to [registry_locations()]
#' @return a list: `parameters` and `equations` (data frames, with the location
#'   columns added), `locations` (every symbol in the XML), `unresolved`
#'   (character, empty on success) and `unreferenced` (a data frame of equation
#'   id and parameter id)
#' @export
registry_resolve <- function(raw, xml_dir = NULL, locations = NULL,
                             unreferenced = c("warn", "error", "none"), ...) {
  unreferenced <- match.arg(unreferenced)
  if (inherits(raw, "model_registry")) raw <- raw$raw
  if (is.null(locations)) locations <- registry_locations(xml_dir, ...)
  tabs <- registry_tables(raw)
  parameters <- tabs$parameters
  equations <- tabs$equations
  unresolved <- character()
  blank <- function(n) data.frame(file = rep(NA_character_, n), line = NA_integer_,
                                  body_file = NA_character_,
                                  body_start = NA_integer_, body_end = NA_integer_,
                                  url = NA_character_, initializer = NA_character_,
                                  stringsAsFactors = FALSE)
  take <- function(symbol) {
    hit <- .registry_locate(locations, symbol)
    if (is.null(hit)) {
      unresolved <<- c(unresolved, symbol)
      return(blank(1L))
    }
    hit[, c("file", "line", "body_file", "body_start", "body_end", "url", "initializer")]
  }
  parameters <- cbind(parameters, do.call(rbind, lapply(parameters$code, take)))
  eq_loc <- do.call(rbind, lapply(equations$code, function(s) {
    if (is.na(s)) blank(1L) else take(s)
  }))
  equations <- cbind(equations, eq_loc)
  for (s in unique(unlist(equations$also))) take(s)
  for (s in tabs$diagnostics$code) take(s)
  unresolved <- unique(unresolved)
  if (length(unresolved)) {
    stop("model registry names symbols Doxygen does not know: ",
         paste(unresolved, collapse = ", "), call. = FALSE)
  }
  # What each equation reads, against what its functions reference.
  # An overload shares the name; what any of them references counts, since the
  # registry names the symbol and not one signature. An identity has no function
  # (code ~) and is not checked.
  misses <- list()
  for (i in seq_len(nrow(equations))) {
    if (is.na(equations$code[i])) next
    fns <- c(equations$code[i], equations$also[[i]])
    refs <- unique(unlist(locations$references[locations$symbol %in% fns]))
    for (pid in equations$reads[[i]]) {
      p <- parameters[parameters$id == pid, ]
      if (!is.na(p$literal)) next
      if (any(c(p$code, p$aliases[[1]]) %in% refs)) next
      misses[[length(misses) + 1L]] <- data.frame(equation = equations$id[i], parameter = pid,
                                                 stringsAsFactors = FALSE)
    }
  }
  misses <- if (length(misses)) do.call(rbind, misses) else
    data.frame(equation = character(), parameter = character(), stringsAsFactors = FALSE)
  if (nrow(misses) && unreferenced != "none") {
    msg <- paste0("reads not referenced by the equation's function: ",
                  paste(misses$equation, misses$parameter, sep = ":", collapse = ", "))
    if (unreferenced == "error") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  list(parameters = parameters, equations = equations, locations = locations,
       unresolved = unresolved, unreferenced = misses)
}

# Each `literal:` against the function its `code` names: a literal is a number
# written inside that function (Tetens' 17.27), so it has to appear among the
# numbers in the function's body, read from the source over the body span
# Doxygen reported. Returns the ids of the literals whose number does not
# appear; a literal changed in the code and not here, or the reverse, is one.
# `resolved` is registry_resolve()'s output; `include_dir` is the headers' root.
.registry_literals_unwritten <- function(resolved, include_dir) {
  p <- resolved$parameters
  p <- p[!is.na(p$literal), , drop = FALSE]
  number <- "(?<![A-Za-z_0-9.])[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?"
  written <- vapply(seq_len(nrow(p)), function(i) {
    if (is.na(p$body_file[i]) || is.na(p$body_start[i]) || is.na(p$body_end[i])) return(FALSE)
    body <- readLines(file.path(include_dir, p$body_file[i]), warn = FALSE)[p$body_start[i]:p$body_end[i]]
    body <- sub("//.*$", "", body)
    tokens <- unlist(regmatches(body, gregexpr(number, body, perl = TRUE)))
    any(suppressWarnings(as.numeric(sub("\\.$", "", tokens))) == p$literal[i], na.rm = TRUE)
  }, TRUE)
  p$id[!written]
}
