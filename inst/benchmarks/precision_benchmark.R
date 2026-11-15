## Benchmark of fastrel_inv() against nadiv::makeAinv() (run interactively, or
## with Rscript). Needs the suggested package nadiv.
##
## The time taken to compute the inbreeding coefficients, which the inverse
## needs, grows with the number of ancestors of each individual, hence with
## the depth of the pedigree (Meuwissen and Luo 1992). Two kinds of pedigree
## are therefore used:
## - pedigrees simulated by sim_ped(), in a small population over many
##   generations (deep: every individual has many ancestors);
## - pedigrees of a fixed number of discrete generations (shallow), as in
##   most real data sets.
library(fastrel)
library(nadiv)

## number of generations of each individual (founders: 1)
generations <- function(p) {
  m <- match(p$momid, p$id)
  d <- match(p$dadid, p$id)
  g <- integer(nrow(p))
  for (i in seq_len(nrow(p))) g[i] <- 1L + max(0L, g[m[i]], g[d[i]], na.rm = TRUE)
  g
}

## n individuals in n_gen discrete generations of equal size, with 10% of each
## generation used as sires and the rest as dams
shallow_ped <- function(n, n_gen) {
  size <- n %/% n_gen
  p <- data.frame(id = seq_len(size * n_gen), dadid = 0L, momid = 0L)
  for (g in 2:n_gen) {
    kids <- (g - 1) * size + seq_len(size)
    prev <- (g - 2) * size + seq_len(size)
    n_sires <- max(1L, size %/% 10)
    p$dadid[kids] <- sample(prev[seq_len(n_sires)], size, replace = TRUE)
    p$momid[kids] <- sample(prev[-seq_len(n_sires)], size, replace = TRUE)
  }
  p
}

time_it <- function(expr) system.time(expr)[["elapsed"]]

compare <- function(p, label) {
  t_fastrel <- time_it(Ainv <- fastrel_inv(p))
  t_nadiv <- time_it(N <- suppressWarnings(
    makeAinv(data.frame(ID = p$id, Dam = p$momid, Sire = p$dadid))$Ainv))
  ## compared as sparse matrices (makeAinv() names the rows only, and may
  ## reorder them)
  k <- match(as.character(p$id), rownames(N))
  same <- max(abs(N[k, k] - Ainv)) < 1e-8 * max(abs(Ainv))
  data.frame(pedigree = label, n = nrow(p), generations = max(generations(p)),
             fastrel_inv = t_fastrel, nadiv = t_nadiv, same = same)
}

results <- list()
set.seed(1)
for (n in c(1000, 5000, 10000, 20000))
  results[[length(results) + 1]] <- compare(sim_ped(n), "sim_ped (deep)")
for (n in c(10000, 50000, 100000))
  results[[length(results) + 1]] <- compare(shallow_ped(n, 12), "12 generations")
results <- do.call(rbind, results)
print(results, row.names = FALSE)   # times in seconds
