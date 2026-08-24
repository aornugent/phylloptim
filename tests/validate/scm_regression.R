# The SCM regression, which is the check that mattered there --
# an adaptive stepper and a discrete node-splitting schedule could turn a
# below-tolerance perturbation into a visible difference. It did not: 78/78 nodes.
#
#   R_LIBS=<lib> Rscript tests/validate/scm_regression.R <output.rds>
#
# Runs a TF24 size-structured (SCM) simulation and writes the numeric results to
# an .rds. Run it twice -- once against a plant built with its own leaf, once
# against a plant built on `feature/consume-leaf-package` -- then compare with
# scm_compare.R.
#
# WHY this is worth doing separately from the operating-point grid. The grid showed
# agreement to 1 ULP on a single leaf solve. The SCM integrates that solve millions
# of times through an adaptive ODE solver with a node-splitting schedule, and both
# of those have thresholds: a perturbation far below tolerance can still change
# which side of a refinement decision a step lands on, and that changes the node
# schedule, which is discrete. So the question the grid cannot answer is whether
# 1 ULP stays 1 ULP, or whether it is amplified by the discrete machinery on top.
#
# The scenario is lifted from plant's own tests/testthat/test-strategy-tf24.R so it
# is a configuration known to run and known to be checked. max_patch_lifetime is 5
# there, which keeps this to seconds; raise it for a more demanding check.

suppressMessages(library(plant))

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args)) args[[1]] else "scm_out.rds"

p0 <- scm_base_parameters("TF24")
env <- Environment("TF24")
ctrl <- Control()
p0$max_patch_lifetime <- 5

cat("plant from:", find.package("plant"), "\n")

# --- one species --------------------------------------------------------------
p1 <- add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                     hyperpar = TF24_hyperpar, birth_rate = list(20))
out1 <- run_scm(p1, env, ctrl, collect = TRUE)

# --- two species, so competitive exclusion is exercised too -------------------
p2 <- add_strategies(p0, trait_matrix(c(0.0825, 0.10, 5, 5), c("lma", "hmat")),
                     hyperpar = TF24_hyperpar, birth_rate = list(20, 20))
out2 <- run_scm(p2, env, ctrl, collect = TRUE)

# Record the aggregate outputs, the full ODE time sequence (which is where the
# adaptive stepper's decisions show up), and the collected per-node state -- the
# last being the most sensitive thing available, since it is the raw trajectory
# rather than an integral over it.
result <- list(
  plant_lib = find.package("plant"),
  one = list(
    offspring_production = out1$offspring_production,
    net_reproduction_ratios = out1$net_reproduction_ratios,
    ode_times = out1$ode_times,
    n_ode_times = length(out1$ode_times),
    species = out1$species,
    env = out1$env
  ),
  two = list(
    offspring_production = out2$offspring_production,
    net_reproduction_ratios = out2$net_reproduction_ratios,
    ode_times = out2$ode_times,
    n_ode_times = length(out2$ode_times),
    species = out2$species
  )
)

saveRDS(result, out_path)
cat(sprintf("wrote %s\n", out_path))
cat(sprintf("  1 species: offspring=%.10f, %d ode times\n",
            result$one$offspring_production, result$one$n_ode_times))
cat(sprintf("  2 species: offspring=%.10f / %.6g, %d ode times\n",
            result$two$offspring_production[[1]],
            result$two$offspring_production[[2]],
            result$two$n_ode_times))
