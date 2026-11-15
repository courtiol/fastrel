# Package index

## Relatedness of a whole pedigree

- [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  : Relatedness matrix from a pedigree

## Inverse of the relatedness matrix

The precision matrix used by animal models, sparse and built directly
from the pedigree.

- [`fastrel_inv()`](https://courtiol.github.io/fastrel/reference/fastrel_inv.md)
  : Inverse of the relatedness matrix (precision matrix)

## Relatedness updated step by step (the store)

A store keeps the relatedness matrix among the individuals currently in
it, and is updated as individuals are added and removed.

- [`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_copy()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_ids()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_size()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`rel_capacity()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  [`print(`*`<relatedness_store>`*`)`](https://courtiol.github.io/fastrel/reference/relstore.md)
  : Incremental store of the relatedness matrix
- [`rel_save()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  [`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  : Save and restore a store

## Measures of relatedness

- [`relatedness_metrics`](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
  : Measures of relatedness

## Simulating data

- [`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md)
  : Simulate a pedigree
- [`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md)
  : Simulate a pedigree with deaths, emigration and immigration

## Package

- [`fastrel`](https://courtiol.github.io/fastrel/reference/fastrel-package.md)
  [`fastrel-package`](https://courtiol.github.io/fastrel/reference/fastrel-package.md)
  : fastrel: fast relatedness, relationship and kinship matrices from
  pedigrees
