# fastrel 0.2.2

* `relatedness_fast()` is 2 to 4 times faster: each pair of individuals is now
  computed once instead of twice, and pedigrees whose rows are not in birth
  order no longer pay a cache penalty. Memory use is unchanged.
* `rel_add()` and `rel_remove()` are faster: they only touch the slots that have
  been in use, not the whole capacity, so an oversized store (e.g. the default
  `rel_new()`) no longer slows them down; adding founders or immigrants takes
  constant time; and the store is looked up once per call instead of six times.
  `rel_get()` no longer zero-fills the matrix before writing it.
* `sim_pop()` is about 10 times faster (same output for a given seed).
* New `rel_ids()` listing the individuals in a store, to be used instead of
  `rel$ids`.
* `sim_ped()` and `sim_pop()` accept `seed = NULL` to leave the random number
  generator of the session as it is.
* `rel_add()` leaves the store unchanged when it fails, whether because of a
  cycle among the new individuals (some of them could already have been
  inserted) or because memory ran out while the store was growing.
* `relatedness_fast()` now refuses `NA` in `id` (unknown parents were matched
  to that individual, which could give a wrong matrix) and any `chrtype` other
  than `"autosome"` (the autosomal matrix was returned silently).
* A store read back with `readRDS()` or `load()` is now reported as such, with
  a pointer to `rel_save()` and `rel_load()`, instead of failing inside Rcpp.
* More input checks: non-numeric or too large `capacity` in `rel_new()`;
  `n_founders` or `window` below 1 in `sim_ped()` (which could silently make an
  individual its own parent) and `p_unknown` outside `[0, 1]`; `n_steps`,
  `n_init` and `maturity` in `sim_pop()`; the `format` of a saved store in
  `rel_load()`; `NA` slots in the C++ object.
* Documentation: the help pages of the store and of the package now mention
  `rel_save()`, `rel_load()`, `sim_pop()`, `rel_size()` and `rel_capacity()`.
* Benchmarks: the store is now also timed on a full build, not only on
  extraction, and a wrong comment on the capacity was corrected.
* `stats` and `utils` are now declared in `Imports`.
* The package is checked on GitHub Actions.

# fastrel 0.2.1

* Revised README

# fastrel 0.2.0

* New `rel_save()` and `rel_load()` to checkpoint and restore a store.
* New `sim_pop()` simulating a pedigree with deaths, emigration and immigration,
  to benchmark `rel_remove()` (see `inst/benchmarks/sim_pop_benchmark.R`).

# fastrel 0.1.0

* First version: `relatedness_fast()`, the incremental store (`rel_new()`,
  `rel_add()`, `rel_remove()`, `rel_get()`) and `sim_ped()`.
