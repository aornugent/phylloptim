# Present leaf's plain `//` comments to Doxygen as documentation, without
# touching a single source file.
#
# WHY THIS EXISTS. Doxygen only treats `///`, `//!`, `/** */` and `/*! */` as
# documentation; an ordinary `//` comment is invisible to it. Every comment in
# this package is an ordinary `//` comment, and they are the substantive part --
# the FMA-contraction note in `hydraulic_cost_Sperry`, the four "do not get this
# wrong" warnings in closed_form.hpp, the two-vulnerability-curves hazard at the
# `stem_b`/`stem_c` declarations. Rendering the API without them would render the
# half a reader can already get from the signatures.
#
# The alternative -- rewriting `//` to `///` across ten headers -- is a large
# diff through exactly the files that feature/api-cleanup (PR #15) rewrites, for
# a cosmetic gain. Doing it at render time costs nothing and conflicts with
# nothing.
#
# WHAT IT DOES. Four rules:
#
#   1. PASS THROUGH the `-*-c++-*-` modeline (an Emacs directive, not prose) and
#      anything already written as `///` or `//!` (a deliberate Doxygen comment,
#      which needs no help and must not be double-marked).
#
#   2. The FIRST comment block in a file becomes a `\file` block, so the "why
#      this header exists" preamble lands on the file's own page instead of being
#      attached to whichever declaration happens to follow it.
#
#   3. Every later line-leading `//` becomes `///`. Trailing comments
#      (`double x; // note`) are left alone: Doxygen would attach them to the
#      wrong entity, and they are terse anyway.
#
#   4. INDENTED RUNS BECOME `\verbatim`. This matters more than it sounds. These
#      comments are full of indented ASCII tables, benchmark numbers and usage
#      examples, and Markdown reflows anything indented by less than four spaces
#      into a paragraph -- turning the usage example at the top of leaf.hpp into
#      one unreadable line. A run of lines indented relative to the surrounding
#      prose is reproduced exactly instead.
#
#      ⚠️ ONE RUN PER COMMENT BLOCK. Doxygen 1.9 -- which is what CI installs --
#      handles the FIRST `\verbatim` in a block and then drops the second one's
#      OPEN, reporting `unexpected command endverbatim` at its close. The line it
#      names is in the filtered stream and lands in unrelated code, so it reads
#      like a parser bug somewhere else entirely: the report that found this was
#      120 lines past the comment responsible. Doxygen 1.17 renders the same input
#      in silence, so a local `doxygen` run says nothing.
#
#      Established by a probe header of eight isolated constructs rendered in CI:
#      one run, banner rules, `|` in prose, `|` inside the run, and an emoji are
#      all clean; every block with two runs errors, and only at the SECOND close.
#      So write one indented display per comment block -- put the equations
#      together rather than one on each side of a paragraph. `docs.yml` asserts
#      it, because nothing else can.
#
#      The `/*! \file */` block rule 2 emits is EXEMPT, and that is measured
#      rather than assumed: closed_form.hpp's file block has carried two runs
#      across a long green master. So the assertion counts `///` blocks only.
#
#      A run only OPENS after a blank line and only on a line that is not a list
#      item. Both conditions are load-bearing. Without the first, the hanging
#      indent under a bullet (`  * POSITIVE magnitudes -- ...` followed by
#      four-space continuations, all over leaf_model.hpp) is mistaken for a code
#      block; that produced a real "'\verbatim' command is not allowed in section
#      title" error before the condition was added. Without the second, the
#      bulleted and numbered lists that DO follow a blank line become
#      preformatted text and stop rendering as lists.
#
# ESCAPING. Text outside a verbatim run is escaped, because none of it was
# written with Doxygen in mind: `\`, `@`, `#`, `%`, `&`, `<` and `>` all mean
# something to Doxygen and here they never do. Without it, `#include <phylloptim.hpp>`
# becomes a broken link to an entity called "include" followed by a swallowed
# HTML tag, and the `\int` in the roots.hpp head-loss note becomes an unknown
# command; both were observed before this was added. Verbatim runs are NOT
# escaped -- Doxygen reproduces them literally, so an escape would show up as a
# stray backslash. If you genuinely want a Doxygen command, write a `///`
# comment and rule 1 will leave it alone.
#
# Deliberately POSIX awk -- no gawk-only `match(s, re, arr)` -- because CI runs
# this on ubuntu, where /usr/bin/awk is mawk.
#
# Usage (see Doxyfile): INPUT_FILTER = "awk -f tools/doxygen_filter.awk"

BEGIN {
  file_block_emitted = 0   # has the \file block already been produced?
  in_file_block = 0        # are we inside it right now?
  in_comment = 0           # inside any run of line-leading `//`?
  in_verbatim = 0          # inside an indented run within that comment?
  pending_blanks = 0       # blank lines held back; see emit()
  emitted_any = 0          # has this block emitted a non-blank line yet?
  prefix = ""              # what each emitted line is prefixed with
}

function indent_of(s,   i, c) {
  i = 1
  while (i <= length(s)) {
    c = substr(s, i, 1)
    if (c != " " && c != "\t") break
    i++
  }
  return substr(s, 1, i - 1)
}

# Neutralise every character Doxygen would read as a command or as markup.
# Backslash first, or it would double the escapes added after it.
function escape(s) {
  gsub(/\\/, "\\\\", s)
  gsub(/@/,  "\\@",  s)
  gsub(/#/,  "\\#",  s)
  gsub(/%/,  "\\%",  s)
  gsub(/&/,  "\\&",  s)
  gsub(/</,  "\\<",  s)
  gsub(/>/,  "\\>",  s)
  return s
}

# A `::` that STARTS a line is a wrapped qualified name -- `MultiLayerRoots` at
# the end of one line and `::set_root_network` at the start of the next -- and
# Doxygen reads it as an explicit link request to a global entity, which cannot
# resolve. Escape only that case: `Leaf::optimise_psi_stem_TF` written inline
# still auto-links, which is worth keeping. Doxygen 1.9 on ubuntu errors on this
# where 1.17 does not, so it only ever showed up in CI.
#
# Written out longhand rather than with a backreference: `\1` in sub() is a GNU
# awk extension and CI runs mawk.
function escape_leading_scope(s,   i) {
  i = 1
  while (i <= length(s) && (substr(s, i, 1) == " " || substr(s, i, 1) == "\t")) i++
  if (substr(s, i, 2) == "::") return substr(s, 1, i - 1) "\\::" substr(s, i + 2)
  return s
}

function put(text) {
  if (length(text)) print prefix " " text
  else print prefix
}

function emit_blanks(   i) {
  for (i = 0; i < pending_blanks; i++) put("")
  pending_blanks = 0
}

function close_verbatim() {
  if (in_verbatim) {
    put("\\endverbatim")
    in_verbatim = 0
  }
}

# One line of comment text, with its `//` already stripped. Blank lines are held
# rather than emitted, so that a blank line at the END of an indented run closes
# the run instead of being swallowed into it.
function emit(text) {
  if (text ~ /^[ \t]*$/) {
    pending_blanks++
    return
  }
  indented = (text ~ /^[ \t][ \t]+[^ \t]/)
  # A list item, at any indentation: `* foo`, `- foo`, `+ foo`, `1. foo`.
  bullet = (text ~ /^[ \t]*([*+-]|[0-9]+[.)])[ \t]/)

  if (in_verbatim && indented) {              # the run continues
    emit_blanks()
    print prefix " " text                     # never escaped, never trimmed
    emitted_any = 1
    return
  }
  if (!in_verbatim && indented && !bullet && (pending_blanks > 0 || !emitted_any)) {
    emit_blanks()                             # blanks BEFORE the run are prose
    put("\\verbatim")
    in_verbatim = 1
    print prefix " " text
    emitted_any = 1
    return
  }
  close_verbatim()
  emit_blanks()
  put(escape_leading_scope(escape(text)))
  emitted_any = 1
}

function end_comment() {
  close_verbatim()
  pending_blanks = 0
  emitted_any = 0
  if (in_file_block) {
    print "*/"
    in_file_block = 0
    file_block_emitted = 1
  }
  in_comment = 0
}

{
  ind = indent_of($0)
  body = substr($0, length(ind) + 1)
  is_comment = (substr(body, 1, 2) == "//")

  # Rule 1: not ours to touch.
  if (is_comment &&
      ($0 ~ /-\*-[ \t]*c\+\+[ \t]*-\*-/ ||
       substr(body, 1, 3) == "///" || substr(body, 1, 3) == "//!")) {
    if (in_comment) end_comment()
    print
    next
  }

  if (!is_comment) {
    if (in_comment) end_comment()
    print
    next
  }

  # Strip the `//` and one following space; the REST of the indentation is
  # meaningful and is what rule 4 keys on.
  text = substr(body, 3)
  if (substr(text, 1, 1) == " ") text = substr(text, 2)

  if (!in_comment) {
    in_comment = 1
    if (!file_block_emitted) {                # Rule 2
      in_file_block = 1
      prefix = ind "/*! \\file"
      print prefix
      prefix = ind
    } else {                                  # Rule 3
      prefix = ind "///"
    }
  }

  emit(text)
}

END {
  if (in_comment) end_comment()
}
