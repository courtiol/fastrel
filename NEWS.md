# fastrel (development version)

* Breaking: `rel_new(remember_removed = )` now takes `"none"` (the default),
  `"fingerprint"` (formerly `TRUE`) or `"id"`, instead of `FALSE` or `TRUE`.
* Breaking: `rel_load()` only reads stores saved by this version of fastrel.
* When no column of `ped` has a recognised name, `fastrel_ped()` and
  `rel_add()` still take the first three columns as the individual, its mother
  and its father, but now say so in a message, once per R session.
* New `rel_new(remember_removed = "id")`, which keeps the identifiers of the
  individuals removed rather than a fingerprint of each of them. It takes more
  memory (about 70 bytes per individual removed for short identifiers, against
  about 30) but can never mistake a new identifier for a removed one.
  `rel_save()` and `rel_load()` keep the kind of trace used.

# fastrel 0.7.0

* New argument `sim_ped(numeric_ids = TRUE)` to simulate a pedigree with
integer identifiers and 0 for unknown parents, which can be given as it is to
`fastrel_ped()`, `rel_add()`, `kinship2::kinship()` or `ggroups::buildA()`.
* 0 can now be used to define an unknown parent everywhere (synonymous to `NA`).

# fastrel 0.6.2

* The results are now validated in the tests against `kinship2` and `ggroups`,
  two independent implementations, for all four measures and on pedigrees of
  many shapes (inbreeding, selfing, unknown parents, unsorted rows, removals
  from a store, ...).
* The section "Range of values" moves from the help page of the package to
  `?relatedness_metrics`.

# fastrel 0.6.1

* Documentation only. In `?relatedness_metrics`, the table is split in two (it
  was too wide for the PDF manual) and the link to `cov2cor()` is fixed. The
  README is refreshed and explains how to build the manual.

# fastrel 0.6.0

* `fastrel_ped()`, `fastrel_vec()` and `rel_get()` gain an argument `metric`
  to return, besides the additive relationship (`"additive"`, the default),
  the kinship coefficient (`"kinship"`), Wright's coefficient of relationship
  (`"wright"`) or Hamilton's regression coefficient of relatedness
  (`"hamilton"`). The conversion uses no extra memory and adds little time.
* New help page `?relatedness_metrics` defining these measures, with
  references.
* Breaking: `fastrel_vec()` has the arguments `sex` and `chrtype` of
  `kinship2::kinship()` instead of `...`; any other argument is now an error.

# fastrel 0.5.1

* The documentation states the range of relatedness values (0 to 2, above 1
  only with inbreeding).

# fastrel 0.5.0

* Tibbles are supported. This fixes `fastrel_vec()`, which returned a wrong
  1 x 1 matrix when each argument was a tibble with one column instead of a
  vector.
* Inputs that cannot be a vector of identifiers now raise a clear error.

# fastrel 0.4.0

* Breaking: `relatedness_fast()` is replaced by `fastrel_vec()` (same
  function) and `fastrel_ped()`, which takes the pedigree as a data frame; the
  argument `carrying` of `sim_pop()` is now `K`.
* `fastrel_ped()` and `rel_add()` find the columns for the individual, its
  mother and its father by their names (`id`, `mom`, `dad`, `sire`, `dam`,
  ...), whatever their order; the other columns are ignored.
* New benchmark of memory use, `inst/benchmarks/memory_benchmark.R`.
* The documentation defines "store" and "slot", and its memory figures were
  corrected after measuring them.

# fastrel 0.3.1

* New argument `rel_new(remember_removed = TRUE)` to refuse adding again an
  individual that was removed from a store, which would otherwise come back
  as unrelated to its relatives. Recommended for metapopulations, whose
  handling is now documented.

# fastrel 0.3.0

* Breaking: `rel_new()` no longer has a `capacity` argument, as the store now
  sizes itself; `sim_ped()` and `sim_pop()` no longer have a `seed` argument
  (call `set.seed()` beforehand).
* Much faster `relatedness_fast()`, `rel_add()`, `rel_remove()` and
  `sim_pop()`, without using more memory.
* New `rel_ids()` listing the individuals in a store.
* Bug fixes: `rel_add()` leaves the store unchanged when it fails; `NA`
  identifiers, invalid arguments of the simulators and stores read back with
  `readRDS()` now raise clear errors instead of giving wrong results.
* The package is checked on GitHub Actions.

# fastrel 0.2.1

* Revised README.

# fastrel 0.2.0

* New `rel_save()` and `rel_load()` to save and restore a store.
* New `sim_pop()` simulating a population with births, deaths, emigration and
  immigration.

# fastrel 0.1.0

* First version: `relatedness_fast()`, the incremental store (`rel_new()`,
  `rel_add()`, `rel_remove()`, `rel_get()`) and `sim_ped()`.
