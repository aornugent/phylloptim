Replace this text. Why the change was needed, then what changed in observable
behaviour — two short paragraphs, about 10 lines. It becomes the squash commit
message verbatim, so write it for someone reading `git log` in three years.

Measurements, benchmark tables, suite counts, alternatives rejected, what you
tried first, and stacking notes go in the **first comment** on this PR, not here:
https://github.com/traitecoevo/plant-meta/blob/main/governance/commit-messages.md
That comment quotes `tests/cpp/bench_solve` against master, interleaved, three
runs each, for any change under `inst/include/` (.claude/measurements.md,
"Benchmark method"). A change to any model's equations, parameters, defaults or
operating-point classification updates `inst/model-registry.yml` in this PR.

Closes #
