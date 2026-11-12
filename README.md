# fastrel

**fast** **rel**atedness (aka `fastrel`) is an R package to compute relatedness 
from a pedigree quickly.

This package offers four main sets of functionalities:

- a set of functions allowing for updating a relatedness matrix, conceived with 
individual based simulations in mind, with no need to store a pedigree, and with
a computation time and memory footprint that do not increase over time as long as
the population is at its demographic equilibrium
(it only requires storing the information about alive individuals)
- a function to compute relatedness from a pedigree more quickly than 
[`ggroups`](https://cran.r-project.org/web/packages/ggroups/index.html) or 
[`Pedixplorer`](https://louislenezet.github.io/Pedixplorer/)
(which now supersedes `kinship2`)
- a function to simulate a pedigree
- a function to simulate a population with its pedigree, deaths, births,
immigrations and emigrations.

The functions have been written with Claude Sonnet 5.5 in C++ and interfaced with R.

The relatedness matrix is here technically defined as the additive relationship matrix, A = 2 × kinship, with diagonal 1 + F.


## Installation

The installation requires the compilation of the C++ code, so you need a C++ 
compiler installed:

- on Linux, those are often installed by default
- on Windows, you will need [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
(same version as your R)
- on MacOS, you will need [Apple Command Line Tools](https://developer.apple.com/download/)
or equivalent.

You may then use [pak](https://pak.r-lib.org/) to install the package simply:

``` r
if (!requireNamespace("pak")) {
  install.packages("pak")
}
pak::pkg_install("courtiol/fastrel")
```

## Basic examples

### Compute the relatedness of a whole pedigree at once


``` r
library(fastrel)
set.seed(1) # for reproducibility
ped <- sim_ped(1000)
tail(ped, n = 3) # last 3 individuals
```

```
##           id   momid   dadid sex
## 998  ind0998 ind0989 ind0959   1
## 999  ind0999 ind0912    <NA>   1
## 1000 ind1000 ind0871 ind0879   1
```

``` r
A <- relatedness_fast(ped$id, dadid = ped$dadid, momid = ped$momid)
A[998:1000, 998:1000] # relatedness matrix between those 3 individuals
```

```
##            ind0998    ind0999    ind1000
## ind0998 1.03676886 0.06467905 0.10907260
## ind0999 0.06467905 1.00000000 0.03918678
## ind1000 0.10907260 0.03918678 1.05214919
```

``` r
max(diag(A) - 1)   # maximum inbreeding coefficients
```

```
## [1] 0.3335501
```

### Compute the relatedness incrementally (for simulations)


``` r
rel <- rel_new() # initialize the store (once per simulation)
rel_add(rel, data.frame(id = c("A", "B"), mom = NA, dad = NA)) # add founders
rel_add(rel, data.frame(id = c("C", "D"), mom = c("A", "A"), dad = c(NA, "B"))) # add births
rel_remove(rel, c("A", "B")) # remove dead / emigrated: free their slots
rel_get(rel) # extract relatedness matrix between remaining individuals
```

```
##      C    D
## C 1.00 0.25
## D 0.25 1.00
```

Removing every individual that can no longer become a parent (dead, emigrated, ...)
implies that the computational cost of adding a birth to the stored relatedness
information and the memory used depend on the number of individuals in the store,
not on the length of the pedigree.

### Saving and loading the relatedness store

Since the store is a C++ object, you need to use specific functions to save and
restore it:


``` r
rel_save(rel, "store.rds")
rel <- rel_load("store.rds") # a working store again
```

### Simulating a population with deaths


``` r
rel <- rel_new()
set.seed(1)
pop <- sim_pop(n_steps = 1000, n_init = 200) # simulate pop (fine tuning the life history is possible)
head(pop) # first 6 simulated individual life histories (no momid and dadid since founder) 
```

```
##            id momid dadid sex arrival exit      fate  origin
## 1 ind00000001  <NA>  <NA>   1       0   20 emigrated founder
## 2 ind00000002  <NA>  <NA>   2       0    7      dead founder
## 3 ind00000003  <NA>  <NA>   1       0    2      dead founder
## 4 ind00000004  <NA>  <NA>   1       0    7      dead founder
## 5 ind00000005  <NA>  <NA>   2       0    3      dead founder
## 6 ind00000006  <NA>  <NA>   1       0   22      dead founder
```

``` r
mean_relatedness <- data.frame(step = 0:1000, mean_rel = NA) # record mean relatedness
for (t in 0:1000) {
  rel_add(rel, pop[pop$arrival == t, c("id", "momid", "dadid")])
  rel_remove(rel, pop$id[!is.na(pop$exit) & pop$exit == t])
  A <- rel_get(rel)
  mean_relatedness[mean_relatedness$step == t, "mean_rel"] <- mean(A[upper.tri(A)])
}
plot(mean_rel ~ step, data = mean_relatedness,
     type = "l", ylab = "Mean relatedness", xlab = "Time step")
```

![plot of chunk simulation](figure/simulation-1.png)

## Benchmarks


``` r
library(microbenchmark)
library(ggroups)

# convert pedigree for ggroups
ID_num <- seq_len(nrow(ped))   # sim_ped() returns individuals in birth order
ped_ggroups <- data.frame(ID = ID_num,
                          SIRE = ID_num[match(ped$dadid, ped$id)],
                          DAM  = ID_num[match(ped$momid, ped$id)])
ped_ggroups$SIRE[is.na(ped_ggroups$SIRE)] <- 0
ped_ggroups$DAM[is.na(ped_ggroups$DAM)]   <- 0

# check that the computed relatedness are the same
all.equal(unname(as.matrix(buildA(ped_ggroups))),
          unname(relatedness_fast(ped$id, dadid = ped$dadid, momid = ped$momid)))
```

```
## [1] TRUE
```

``` r
# benchmark
microbenchmark(ggroups = buildA(ped_ggroups), 
               fastrel = relatedness_fast(ped$id, dadid = ped$dadid, momid = ped$momid))
```

```
## Unit: milliseconds
##     expr        min         lq       mean     median         uq       max neval cld
##  ggroups 199.337294 207.697344 216.820916 213.464191 220.700567 483.82672   100  a 
##  fastrel   2.744712   3.684355   6.263514   4.608985   7.436562  19.44156   100   b
```
As you can see `fastrel::relatedness_fast()` is much faster than `ggroups::buildA()` on our machine, for 1000 individuals.


You can also run the integrated benchmarks as follows:


``` r
source(system.file("benchmarks/relstore_benchmark.R", package = "fastrel"))
source(system.file("benchmarks/sim_pop_benchmark.R", package = "fastrel"))
```


## Development

After editing README.Rmd, do:


``` r
unlink("figure", recursive = TRUE)
knitr::knit("README.Rmd")
```

After editing C++ attributes, do:


``` r
Rcpp::compileAttributes()
```

## Thanks

Nikos Smit, Leonie F. Walter.
