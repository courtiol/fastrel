# Simulate a pedigree

Simulates a pedigree with overlapping generations, some inbreeding, and
some unknown parents. Individuals are numbered in birth order, so that
parents always come before their offspring, and each individual picks
its parents among the `window` individuals born just before it. It is
meant for testing and benchmarking.

## Usage

``` r
sim_ped(
  n = 1000,
  n_founders = 40,
  window = 150,
  p_unknown = 0.05,
  numeric_ids = FALSE
)
```

## Arguments

- n:

  Total number of individuals (larger than `n_founders`).

- n_founders:

  Number of individuals with no parents at the start.

- window:

  Number of preceding individuals among which parents are drawn.

- p_unknown:

  Probability that a parent (mother or father separately) is unknown for
  a non-founder.

- numeric_ids:

  If `TRUE`, individuals are identified by the integers 1 to `n` (in
  birth order) and unknown parents by 0, as many programs of animal
  breeding expect (e.g.
  [`ggroups::buildA()`](https://rdrr.io/pkg/ggroups/man/buildA.html)).
  If `FALSE` (the default), identifiers are character strings and
  unknown parents are `NA`.

## Value

A data frame with columns `id`, `dadid`, `momid` and `sex` (1 = male, 2
= female). Identifiers are character strings such as `"ind0001"`, with
`NA` for unknown parents, or with `numeric_ids = TRUE` integers, with 0
for unknown parents.

Both forms can be given as they are to
[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
[`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
and
[`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md),
which all take 0 as an unknown parent.

## Details

The result depends on the state of the random number generator of the
session: call [`set.seed()`](https://rdrr.io/r/base/Random.html)
beforehand to make it reproducible.

## Examples

``` r
set.seed(1)
ped <- sim_ped(200)
head(ped)
#>        id dadid momid sex
#> 1 ind0001  <NA>  <NA>   1
#> 2 ind0002  <NA>  <NA>   2
#> 3 ind0003  <NA>  <NA>   1
#> 4 ind0004  <NA>  <NA>   1
#> 5 ind0005  <NA>  <NA>   2
#> 6 ind0006  <NA>  <NA>   1
A <- fastrel_ped(ped)
range(diag(A) - 1)
#> [1] 0.00 0.25

## the same pedigree with integer identifiers and 0 for unknown parents
set.seed(1)
ped_num <- sim_ped(200, numeric_ids = TRUE)
head(ped_num)
#>   id dadid momid sex
#> 1  1     0     0   1
#> 2  2     0     0   2
#> 3  3     0     0   1
#> 4  4     0     0   1
#> 5  5     0     0   2
#> 6  6     0     0   1
all.equal(unname(fastrel_ped(ped_num)), unname(A))
#> [1] TRUE

## ready for other packages, e.g. ggroups (columns: individual, sire, dam)
if (requireNamespace("ggroups", quietly = TRUE)) {
  all.equal(unname(ggroups::buildA(ped_num[, c("id", "dadid", "momid")])),
            unname(A))
}
#> [1] TRUE
```
