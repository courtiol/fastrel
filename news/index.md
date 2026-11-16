# Changelog

## fastrel 1.0.0

- Breaking:
  [`fastrel_inv()`](https://courtiol.github.io/fastrel/reference/fastrel_inv.md)
  now returns a `dgCMatrix`, which stores both triangles of the matrix,
  instead of a `dsCMatrix`, which stores only one:
  [`MCMCglmm::MCMCglmm()`](https://rdrr.io/pkg/MCMCglmm/man/MCMCglmm.html)
  (version 2.36) silently gives wrong results with the latter. The new
  argument `matrix_class = "dsCMatrix"` still returns the leaner
  symmetric class, for software that reads it correctly, such as spaMM.
- [`fastrel_inv()`](https://courtiol.github.io/fastrel/reference/fastrel_inv.md)
  is tested in animal models fitted with spaMM and MCMCglmm.
- The documentation calls the additive relationship the genetic
  covariance (Emik and Terrill 1949), and says between what it is a
  covariance.
- New hex logo.

## fastrel 0.9.0

- Breaking: when no column of `ped` has a recognised name,
  [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  and
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  now take the first three columns as the individual, its father and its
  mother (instead of mother then father). This aligns with the order
  used by
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
  and standard animal breeding pedigrees.
  [`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md)
  and
  [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md)
  return their columns in the same order (`id`, `dadid`, `momid`).
- New website, built with pkgdown:
  <https://courtiol.github.io/fastrel/>.
- New
  [`rel_copy()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  which returns an independent copy of a store (a store is otherwise
  shared, not copied, when it is assigned).
- New
  [`fastrel_inv()`](https://courtiol.github.io/fastrel/reference/fastrel_inv.md),
  which returns the inverse of the relatedness matrix (the precision
  matrix used by animal models) as a sparse matrix, built directly from
  the pedigree. Matrix is now imported.
- The results are now also validated in the tests against
  [`nadiv::makeA()`](https://rdrr.io/pkg/nadiv/man/makeA.html), and the
  README compares the speed of fastrel with that of nadiv.
- Documentation slightly reworked.

## fastrel 0.8.0

- Breaking: `rel_new(remember_removed = )` now takes `"none"` (the
  default), `"fingerprint"` (formerly `TRUE`) or `"id"`, instead of
  `FALSE` or `TRUE`.
- Breaking:
  [`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  only reads stores saved by this version of fastrel.
- When no column of `ped` has a recognised name,
  [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  and
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  still take the first three columns as the individual, its mother and
  its father, but now say so in a message, once per R session.
- New `rel_new(remember_removed = "id")`, which keeps the identifiers of
  the individuals removed rather than a fingerprint of each of them. It
  takes more memory but can never mistake a new identifier for a removed
  one.
  [`rel_save()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  and
  [`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  keep the kind of trace used.
- Preparation for CRAN: references with DOIs in the description of the
  package, stricter checks of the arguments of
  [`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md)
  and
  [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md)
  (single whole numbers or probabilities, no `NA`), canonical CRAN links
  in the README, and unused packages removed from `Suggests`.

## fastrel 0.7.0

- New argument `sim_ped(numeric_ids = TRUE)` to simulate a pedigree with
  integer identifiers and 0 for unknown parents, which can be given as
  it is to
  [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
  or [`ggroups::buildA()`](https://rdrr.io/pkg/ggroups/man/buildA.html).
- 0 can now be used to define an unknown parent everywhere (synonymous
  with `NA`).

## fastrel 0.6.2

- The results are now validated in the tests against `kinship2` and
  `ggroups`, two independent implementations, for all four measures and
  on pedigrees of many shapes (inbreeding, selfing, unknown parents,
  unsorted rows, removals from a store, …).
- The section “Range of values” moves from the help page of the package
  to
  [`?relatedness_metrics`](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md).

## fastrel 0.6.1

- Documentation only. In
  [`?relatedness_metrics`](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md),
  the table is split in two (it was too wide for the PDF manual) and the
  link to [`cov2cor()`](https://rdrr.io/r/stats/cor.html) is fixed. The
  README is refreshed and explains how to build the manual.

## fastrel 0.6.0

- [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  and
  [`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  gain an argument `metric` to return, besides the additive relationship
  (`"additive"`, the default), the kinship coefficient (`"kinship"`),
  Wright’s coefficient of relationship (`"wright"`) or Hamilton’s
  regression coefficient of relatedness (`"hamilton"`). The conversion
  uses no extra memory and adds little time.
- New help page
  [`?relatedness_metrics`](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
  defining these measures, with references.
- Breaking:
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  has the arguments `sex` and `chrtype` of
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
  instead of `...`; any other argument is now an error.

## fastrel 0.5.1

- The documentation states the range of relatedness values (0 to 2,
  above 1 only with inbreeding).

## fastrel 0.5.0

- Tibbles are supported. This fixes
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
  which returned a wrong 1 x 1 matrix when each argument was a tibble
  with one column instead of a vector.
- Inputs that cannot be a vector of identifiers now raise a clear error.

## fastrel 0.4.0

- Breaking: `relatedness_fast()` is replaced by
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  (same function) and
  [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
  which takes the pedigree as a data frame; the argument `carrying` of
  [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md)
  is now `K`.
- [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  and
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  find the columns for the individual, its mother and its father by
  their names (`id`, `mom`, `dad`, `sire`, `dam`, …), whatever their
  order; the other columns are ignored.
- New benchmark of memory use, `inst/benchmarks/memory_benchmark.R`.
- The documentation defines “store” and “slot”, and its memory figures
  were corrected after measuring them.

## fastrel 0.3.1

- New argument `rel_new(remember_removed = TRUE)` to refuse adding again
  an individual that was removed from a store, which would otherwise
  come back as unrelated to its relatives. Recommended for
  metapopulations, whose handling is now documented.

## fastrel 0.3.0

- Breaking:
  [`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  no longer has a `capacity` argument, as the store now sizes itself;
  [`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md)
  and
  [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md)
  no longer have a `seed` argument (call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand).
- Much faster `relatedness_fast()`,
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  and
  [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md),
  without using more memory.
- New
  [`rel_ids()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  listing the individuals in a store.
- Bug fixes:
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  leaves the store unchanged when it fails; `NA` identifiers, invalid
  arguments of the simulators and stores read back with
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) now raise clear
  errors instead of giving wrong results.
- The package is checked on GitHub Actions.

## fastrel 0.2.1

- Revised README.

## fastrel 0.2.0

- New
  [`rel_save()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  and
  [`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  to save and restore a store.
- New
  [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md)
  simulating a population with births, deaths, emigration and
  immigration.

## fastrel 0.1.0

- First version: `relatedness_fast()`, the incremental store
  ([`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md))
  and
  [`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md).
