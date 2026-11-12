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
* `rel_add()` leaves the store unchanged when it fails because of a cycle among
  the new individuals, even if some of them had already been inserted.
* More input checks: `NA` in the `id` of `relatedness_fast()`, non-numeric or
  too large `capacity` in `rel_new()`, and `n_founders` or `window` below 1 in
  `sim_ped()` (which could silently make an individual its own parent).
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
