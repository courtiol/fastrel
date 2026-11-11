# fastrel

Fast additive relatedness matrices (A = 2 x kinship, autosomal loci) from pedigrees,
written in C++ with Rcpp, plus an incremental store for simulations in which the
pedigree keeps growing.

## Install

```r
if (!requireNamespace("pak")) {
  install.packages("pak")
}
pak::pkg_install("courtiol/fastrel")
```

## Whole pedigree at once

```r
library(fastrel)
ped <- sim_ped(1000)
A <- relatedness_fast(ped$id, dadid = ped$dadid, momid = ped$momid)
diag(A) - 1   # inbreeding coefficients
```

## Incremental store (simulations)

```r
rel <- rel_new(capacity = 500)                 # about 1.2 x the peak number of individuals
rel_add(rel, data.frame(id = c("A", "B"), mom = NA, dad = NA))     # founders
rel_add(rel, data.frame(id = "C", mom = "A", dad = "B"))           # focal, mother, father
rel_remove(rel, c("A", "B"))                   # dead / emigrated: free their slots
rel_get(rel, "C")
```

Remove every individual that can no longer become a parent (dead, emigrated, ...):
the cost per birth and the memory depend on the number of individuals in the store,
not on the length of the pedigree.

## Checkpointing

```r
rel_save(rel, "store.rds")
rel <- rel_load("store.rds")      # a working store again
```

## Simulating a population with deaths

```r
pop <- sim_pop(n_steps = 100, n_init = 200)    # arrival / exit steps for everyone
for (t in 0:100) {
  rel_add(rel, pop[pop$arrival == t, c("id", "momid", "dadid")])
  rel_remove(rel, pop$id[!is.na(pop$exit) & pop$exit == t])
}
```

## Development

After editing roxygen comments or the C++ attributes:

```r
Rcpp::compileAttributes()
roxygen2::roxygenise()
```

Benchmarks: `system.file("benchmarks/relstore_benchmark.R", package = "fastrel")`.
