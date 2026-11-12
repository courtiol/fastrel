# fastrel 0.3.1

* `rel_add()` now refuses an identifier that was removed from the store
  earlier. Such an individual used to come back as an unrelated founder, which
  was silently wrong for its relatives. The store remembers the identifiers
  removed for this (about 70 bytes each), and `rel_save()` / `rel_load()`
  keep them.
* The help page of the store explains how to handle a metapopulation: a single
  store for all the patches, in which migrants are kept.

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
