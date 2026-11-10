# fastrel 0.2.0

* New `rel_save()` and `rel_load()` to checkpoint and restore a store.
* New `sim_pop()` simulating a pedigree with deaths, emigration and immigration,
  to benchmark `rel_remove()` (see `inst/benchmarks/sim_pop_benchmark.R`).

# fastrel 0.1.0

* First version: `relatedness_fast()`, the incremental store (`rel_new()`,
  `rel_add()`, `rel_remove()`, `rel_get()`) and `sim_ped()`.
