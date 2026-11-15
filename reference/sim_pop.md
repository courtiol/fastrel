# Simulate a pedigree with deaths, emigration and immigration

Simulates an individual-based population over discrete time steps and
returns the complete pedigree together with the step at which each
individual arrived and left. It is meant for testing and benchmarking
[`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
and
[`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md),
since the number of individuals alive stays roughly constant while the
pedigree keeps growing.

## Usage

``` r
sim_pop(
  n_steps = 50,
  n_init = 100,
  K = 2 * n_init,
  maturity = 2,
  p_birth = 0.6,
  p_death = 0.1,
  p_emigrate = 0.02,
  n_immigrants = 1
)
```

## Arguments

- n_steps:

  Number of time steps after step 0.

- n_init:

  Number of founders at step 0.

- K:

  Carrying capacity: the birth probability decreases linearly to 0 when
  the number of individuals alive, `N`, reaches it.

- maturity:

  Number of steps after birth before an individual can reproduce or
  emigrate. Founders and immigrants are mature on arrival.

- p_birth:

  Maximum probability that a mature female reproduces in a step (one
  offspring at most per step).

- p_death:

  Probability of death per step.

- p_emigrate:

  Probability of emigration per step for mature individuals that did not
  die.

- n_immigrants:

  Mean number of immigrants per step (Poisson).

## Value

A data frame with one row per individual ever present, in order of
arrival, so parents always come before their offspring, and with
columns:

- `id`, `dadid`, `momid`: identifiers (`NA` for unknown parents);

- `sex`: 1 = male, 2 = female;

- `arrival`: step at which the individual enters the pedigree (0 for
  founders);

- `exit`: step at which it leaves, `NA` if it is still present at the
  end;

- `fate`: `"dead"`, `"emigrated"` or `NA`;

- `origin`: `"founder"`, `"born"` or `"immigrant"`.

## Details

Step 0 creates `n_init` adult founders with unknown parents. At each
following step `t`:

1.  Each mature female alive at the start of the step reproduces with
    probability `p_birth * (1 - N / K)` (bounded to `[0, 1]`), where `N`
    is the number of individuals alive. The father is drawn at random
    among the mature males (`NA` if there is none). Newborns arrive at
    step `t`.

2.  A Poisson number (mean `n_immigrants`) of adult immigrants with
    unknown parents arrive at step `t`.

3.  Each individual that was alive at the start of the step dies with
    probability `p_death`, and otherwise, if mature, emigrates with
    probability `p_emigrate`. Its `exit` is `t`.

Newborns and immigrants cannot leave in the step in which they arrive,
and parents are always alive at the start of the step in which their
offspring are born. Therefore, replaying the steps in order with
[`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
for the arrivals followed by
[`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md)
for the exits never refers to an unknown individual or to a removed
parent (see the example).

The result depends on the state of the random number generator of the
session: call [`set.seed()`](https://rdrr.io/r/base/Random.html)
beforehand to make it reproducible.

## See also

[`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md),
[`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md)

## Examples

``` r
set.seed(2)
pop <- sim_pop(n_steps = 40, n_init = 60)
table(pop$origin)
#> 
#>      born   founder immigrant 
#>       271        60        45 
table(pop$fate, useNA = "ifany")
#> 
#>      dead emigrated      <NA> 
#>       271        43        62 

## replay the simulation with the incremental store
rel <- rel_new()
size <- integer(41)
for (t in 0:40) {
  rel_add(rel, pop[pop$arrival == t, ])   # the columns needed are found by name
  rel_remove(rel, pop$id[!is.na(pop$exit) & pop$exit == t])
  size[t + 1] <- rel_size(rel)
}
size
#>  [1] 60 62 58 58 52 55 55 59 60 60 61 65 67 73 78 75 73 78 79 74 72 77 68 70 66
#> [26] 70 70 71 71 69 64 67 58 67 70 65 65 63 61 61 62

## the survivors agree with the relatedness computed from the whole pedigree
A <- fastrel_ped(pop)
alive <- rel_ids(rel)
all.equal(rel_get(rel), A[alive, alive])
#> [1] TRUE
```
