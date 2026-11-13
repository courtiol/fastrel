# fastrel 0.4.0

## Breaking changes

* `relatedness_fast()` is replaced by two functions: `fastrel_vec()`, which is
  the same function under a new name, and `fastrel_ped()`, which takes the
  pedigree as a data frame.
* The argument `carrying` of `sim_pop()` is now called `K`.
* `rel_add()` now finds the columns of `ped` by their names (see below). A data
  frame in which only some of the three columns have a recognised name (e.g.
  `id`, `m`, `d`) used to be read by position and is now refused.

## New features

* `fastrel_ped()` and `rel_add()` find the columns for the focal individual,
  its mother and its father by their names, whatever their order and their
  capitalisation, and ignore the other columns. Usual names are recognised
  (`id`, `focal`, `animal`, ...; `mom`, `mum`, `mother`, `dam`, ...; `dad`,
  `father`, `sire`, ...; with or without `id`, as in `momid` or `mother_id`).
  If no column has a recognised name, the first three columns are used as
  before. The data frames returned by `sim_ped()` and `sim_pop()` can thus be
  passed as they are.
* Tibbles are supported and tested. In particular `fastrel_vec()`,
  `rel_get()` and `rel_remove()` accept a data frame with a single column where
  a vector is expected, which is what `x[, "id"]` returns when `x` is a tibble:
  `fastrel_vec()` used to return a wrong 1 x 1 matrix in that case.
* New benchmark `inst/benchmarks/memory_benchmark.R` measuring the memory used
  against what the documentation says (needs the package lobstr).

## Documentation

* The words "store" and "slot" are now defined, and the help page explains what
  a store is made of: an R environment that points to a C++ object holding the
  matrix.
* Memory figures were corrected after measuring them: `remember_removed = TRUE`
  costs about 30 bytes per individual removed (not 11 to 23), and a store that
  grew in small steps made the R process use about 20% more than its matrix.

# fastrel 0.3.1

* New argument `rel_new(remember_removed = TRUE)`: the store then keeps a trace
  of the individuals removed and `rel_add()` refuses to add them again. Such an
  individual otherwise comes back as an unrelated founder, which is silently
  wrong for its relatives. It is off by default because it costs 11 to 23 bytes
  per individual removed, so the memory used then grows with the length of the
  pedigree; it is recommended for metapopulations. `rel_save()` / `rel_load()`
  keep the setting and the trace.
* The help page of the store and the README explain how to handle a
  metapopulation: a single store for all the patches, in which migrants are
  kept.

# fastrel 0.3.0

## Breaking changes

* `rel_new()` no longer has a `capacity` argument: the store allocates and
  enlarges its matrix by itself. It now grows by half instead of doubling (or
  straight to the size needed by a large `rel_add()`), so the matrix stays
  within about 1.5 times the peak number of individuals held, where the former
  default started at 1000 slots and doubling could overshoot by a factor 2.
  `rel_capacity()` still reports the number of slots.
* `sim_ped()` and `sim_pop()` no longer have a `seed` argument and no longer
  call `set.seed()`: call `set.seed()` yourself beforehand. `set.seed(1)`
  followed by a call gives the same result as the former default.
* `relatedness_fast()` now refuses a `chrtype` other than `"autosome"` (the
  autosomal matrix was returned silently).

## New features

* New `rel_ids()` listing the individuals in a store, to be used instead of
  `rel$ids`.
* The package is checked on GitHub Actions.

## Performance

* `relatedness_fast()` is 2 to 4 times faster: each pair of individuals is now
  computed once instead of twice, and pedigrees whose rows are not in birth
  order no longer pay a cache penalty. Memory use is unchanged.
* `rel_add()` and `rel_remove()` are faster: they only touch the slots that have
  been in use, not all the slots allocated; adding founders or immigrants takes
  constant time; and the store is looked up once per call instead of six times.
  `rel_get()` no longer zero-fills the matrix before writing it.
* `sim_pop()` is about 10 times faster.

## Bug fixes

* `rel_add()` leaves the store unchanged when it fails, whether because of a
  cycle among the new individuals (some of them could already have been
  inserted) or because memory ran out while the store was growing.
* `relatedness_fast()` now refuses `NA` in `id` (unknown parents were matched
  to that individual, which could give a wrong matrix).
* `sim_ped()` refuses `n_founders` or `window` below 1, which could silently
  make an individual its own parent, and `p_unknown` outside `[0, 1]`.
* A store read back with `readRDS()` or `load()` is now reported as such, with
  a pointer to `rel_save()` and `rel_load()`, instead of failing inside Rcpp.
* More input checks: `n_steps`, `n_init` and `maturity` in `sim_pop()`; the
  `format` of a saved store in `rel_load()`; `NA` slots in the C++ object.

## Documentation and packaging

* The help pages of the store and of the package now mention `rel_save()`,
  `rel_load()`, `sim_pop()`, `rel_size()` and `rel_capacity()`.
* Benchmarks: the store is now also timed on a full build, not only on
  extraction.
* `stats` and `utils` are now declared in `Imports`.

# fastrel 0.2.1

* Revised README

# fastrel 0.2.0

* New `rel_save()` and `rel_load()` to checkpoint and restore a store.
* New `sim_pop()` simulating a pedigree with deaths, emigration and immigration,
  to benchmark `rel_remove()` (see `inst/benchmarks/sim_pop_benchmark.R`).

# fastrel 0.1.0

* First version: `relatedness_fast()`, the incremental store (`rel_new()`,
  `rel_add()`, `rel_remove()`, `rel_get()`) and `sim_ped()`.
