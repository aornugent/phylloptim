#!/usr/bin/env bash
# Run the R-side and C++-side benchmarks against a list of commits and print one
# table, so "have we regressed?" is answerable rather than remembered.
#
# The R script (tools/bench_user_cost.R) is always taken from the CURRENT tree and
# pointed at an OLD install: it is version-aware, and copying an old copy of it into
# each worktree would compare different measurements as if they were the same one.
#
# ⚠️ Compare the RATIO columns, not the microseconds. A bare .Call on this machine
# moved 0.69 -> 1.10 us between two runs an hour apart on the same build, so an
# absolute column is only meaningful within a single invocation of this script.
#
# Usage: tools/bench_history.sh <sha> [<sha> ...]
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
WORK=${TMPDIR:-/tmp}/phylloptim-bench-$$
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

hdr=$(cd "$HERE" && Rscript tools/bench_user_cost.R --tsv --reps 1 2>/dev/null | head -1)
printf 'commit\tsubject\t%s\tcpp_solve_us\tcpp_grad_ift_us\tcpp_grad_fd_us\n' "$hdr"

for sha in "$@"; do
  subj=$(git -C "$HERE" log -1 --format=%s "$sha" | cut -c1-40)
  wt="$WORK/wt-$sha"; lib="$WORK/lib-$sha"; mkdir -p "$lib"
  git -C "$HERE" worktree add --quiet --detach "$wt" "$sha" 2>/dev/null || { echo "$sha: worktree failed" >&2; continue; }

  # ⚠️ Commits before #47 declare `Package: leaf`. They install cleanly, `library(
  # phylloptim)` then falls back to the SITE build, and the row that comes out
  # describes the CURRENT tree while looking like history. Refuse rather than report.
  pkg=$(grep -m1 '^Package:' "$wt/DESCRIPTION" | awk '{print $2}')
  ver=$(grep -m1 '^Version:' "$wt/DESCRIPTION" | awk '{print $2}')
  if [ "$pkg" != "phylloptim" ]; then
    echo "$sha: SKIPPED -- declares Package: $pkg, not phylloptim (pre-#47); needs its own harness" >&2
    git -C "$HERE" worktree remove --force "$wt" 2>/dev/null
    continue
  fi

  R CMD INSTALL --no-docs -l "$lib" "$wt" >/dev/null 2>&1 || { echo "$sha: install failed" >&2; git -C "$HERE" worktree remove --force "$wt"; continue; }

  # R side: current script, old library -- and the script hard-fails if it ends up
  # loading anything other than the version we just installed.
  rrow=$(cd "$HERE" && R_LIBS="$lib" Rscript tools/bench_user_cost.R --tsv --expect-version "$ver" 2>/dev/null | tail -1)
  if [ -z "$rrow" ]; then
    echo "$sha: R bench refused or failed (wrong library loaded?)" >&2
    git -C "$HERE" worktree remove --force "$wt" 2>/dev/null
    continue
  fi

  # C++ side: built from the worktree, so it measures that commit's headers.
  cpp_solve=NA; cpp_ift=NA; cpp_fd=NA
  if make -C "$wt/tests/cpp" bench_solve >/dev/null 2>&1; then
    cpp_solve=$("$wt/tests/cpp/bench_solve" 2>/dev/null | grep -oE '[0-9]+\.[0-9]+ us/solve' | grep -oE '^[0-9.]+')
  fi
  if make -C "$wt/tests/cpp" bench_gradient >/dev/null 2>&1; then
    line=$("$wt/tests/cpp/bench_gradient" 2>/dev/null | grep -E '^vcmax_25' | head -1)
    cpp_fd=$(awk '{print $2}' <<<"$line"); cpp_ift=$(awk '{print $3}' <<<"$line")
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(git -C "$HERE" rev-parse --short "$sha")" "$subj" "$rrow" "${cpp_solve:-NA}" "${cpp_ift:-NA}" "${cpp_fd:-NA}"
  git -C "$HERE" worktree remove --force "$wt" 2>/dev/null
done
