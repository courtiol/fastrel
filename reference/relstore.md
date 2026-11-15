# Incremental store of the relatedness matrix

A *store* is an object that keeps the additive relationship matrix \\A =
2\phi\\ among the individuals currently in it, and that is updated as
individuals are added and removed. Because relatedness between existing
individuals never changes when new individuals are born, a newborn's row
is computed once from its parents' rows, at a cost proportional to the
number of individuals the store has held at once instead of recomputing
the whole matrix. Individuals who can no longer become parents (dead,
emigrated, ...) should be removed so that their slots are recycled. The
cost per birth and the memory then depend on the number of individuals
in the store, not on the length of the pedigree.

## Usage

``` r
rel_new(remember_removed = c("none", "fingerprint", "id"))

rel_copy(rel)

rel_add(rel, ped)

rel_remove(rel, ids)

rel_get(
  rel,
  ids = rel_ids(rel),
  metric = c("additive", "kinship", "wright", "hamilton", "Wright", "Hamilton")
)

rel_ids(rel)

rel_size(rel)

rel_capacity(rel)

# S3 method for class 'relatedness_store'
print(x, ...)
```

## Arguments

- remember_removed:

  Should the store keep a trace of the individuals removed, so that
  `rel_add()` refuses to add them again? `"none"` (the default) for no
  trace, `"fingerprint"` for a compact trace that is wrong with a
  negligible probability, or `"id"` for the identifiers themselves,
  which takes more memory but is never wrong (see the section Removed
  individuals).

- rel:

  A store created by `rel_new()`.

- ped:

  A data frame (or a tibble) with one row per individual to add and
  three columns giving the identifier of the individual, of its father
  and of its mother (`NA` or 0 for unknown parents). The columns are
  found by their names (e.g. `id`, `dad`, `mom` or `Animal`, `Sire`,
  `Dam`), whatever their order, and other columns are ignored: see the
  section Columns of `ped` in
  [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md).
  Only individuals not yet in the store can be added. Known parents must
  be in the store or in `ped` itself: any other identifier raises an
  error rather than being treated as unknown, since it usually reveals a
  mistake (a typo, a parent never added or removed too early). 0 is
  taken as an unknown parent unless an individual is called 0, in which
  case it is refused as ambiguous. The rows of `ped` may be in any
  order.

- ids:

  Character vector (or coercible to it) of individual identifiers. A
  data frame with a single column is taken as that column. For
  `rel_get()`, the default is every individual in the store.

- metric:

  The measure of relatedness returned by `rel_get()`: `"additive"` (the
  default), `"kinship"`, `"wright"` (or `"Wright"`) or `"hamilton"` (or
  `"Hamilton"`). Abbreviations are accepted. See
  [relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md).
  The store always holds the additive relationship, and the conversion
  is done on the matrix extracted only.

- x:

  A store created by `rel_new()`.

- ...:

  Ignored.

## Value

`rel_new()` returns the store. `rel_copy()` returns a new store,
identical to `rel` (same individuals, slots, capacity and trace of the
removed individuals) but independent of it. `rel_add()` and
`rel_remove()` return the store invisibly. `rel_get()` returns a matrix
with `ids` as dimnames: with the default `metric = "additive"`, it is
symmetric with diagonal `1 + F` (see
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
for the other measures). `rel_ids()` returns a character vector, in
order of insertion. `rel_size()` and `rel_capacity()` return integers.

## Details

The functions work as follows:

- `rel_new()` creates an empty store.

- `rel_add()` adds individuals, given by a data frame with columns for
  the focal individual, its father and its mother.

- `rel_remove()` removes individuals.

- `rel_get()` extracts the relatedness matrix among chosen individuals,
  as any of the measures described in
  [relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md).

- `rel_copy()` returns an independent copy of a store.

- `rel_ids()` lists the individuals in the store.

- `rel_size()` and `rel_capacity()` report the number of individuals in
  the store and the number of slots allocated.

## Vocabulary

*Store* and *slot* are everyday words of computer science, not terms of
R or of C++.

- A **store** is a container that keeps data in memory between calls,
  here the relatedness matrix. You create it once and then pass it to
  the other functions, which read it or update it.

- A **slot** is one position in the matrix of the store, i.e. one row
  and the matching column, reserved for one individual. An individual
  occupies a slot from the moment it is added until it is removed; the
  slot is then free and is given to an individual added later. The
  matrix is allocated with more slots than needed so that it does not
  have to be enlarged at each birth: its number of slots is the
  *capacity* of the store.

## What a store is made of

The object returned by `rel_new()` is an R environment (the kind of
object created by
[`new.env()`](https://rdrr.io/r/base/environment.html)) of class
`"relatedness_store"`. It holds the identifiers of the individuals, the
slot of each of them, and a pointer to a C++ object which holds the
matrix itself. Two consequences follow, which differ from the usual
behaviour of R objects:

- Like any R environment, a store is modified in place and is not copied
  when it is assigned or passed to a function: `rel_add()` and
  `rel_remove()` change the store they are given, and after
  `rel2 <- rel` both names refer to the same store. Use `rel_copy()` to
  get a store that can be modified without changing the original.

- The matrix lives in memory managed by C++, outside the memory managed
  by R. It is therefore not written by
  [`save()`](https://rdrr.io/r/base/save.html) or
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html), and a store read
  back from such a file is unusable: use
  [`rel_save()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  and
  [`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  instead. For the same reason
  [`utils::object.size()`](https://rdrr.io/r/utils/object.size.html)
  does not see it: see the section Memory management.

The content of the environment is internal: use the functions above
rather than reading or modifying it directly.

## What to keep in the store

Keep an individual for as long as it can still become a parent, or for
as long as you want to query its relatedness. A removed individual
cannot be used as a parent later (`rel_add()` raises an error), and a
newborn gets no relatedness value with individuals removed before its
birth. Relatedness values among the remaining individuals are not
affected by removals.

Removal is final. The relatedness between a removed individual and those
born afterwards is not tracked, so if it is added again it comes back as
unrelated to everybody, which is wrong for its relatives. By default
nothing prevents this: it is up to you never to add again an individual
that was removed. With `rel_new(remember_removed = "fingerprint")` or
`rel_new(remember_removed = "id")`, `rel_add()` raises an error instead
(see the section Removed individuals).

## Removed individuals

With `remember_removed` set to `"fingerprint"` or `"id"` (rather than
`"none"`, the default), the store keeps a trace of every individual
removed and `rel_add()` refuses to add it again. This protects against a
silent error, at a cost in memory that grows slowly with the length of
the pedigree, instead of depending only on the number of individuals in
the store. This is why it is off by default. Two kinds of trace are
offered:

- `"fingerprint"` keeps a 64-bit fingerprint of each identifier, about
  30 bytes per individual removed (30 MB per million). A new identifier
  can be mistaken for a removed one, but with a negligible probability:
  the number of individuals removed divided by `2^64`, i.e. less than 1
  in 10 million millions per addition after a million removals.

- `"id"` keeps the identifiers themselves, so it can never be wrong, but
  it takes more memory, which increases with the length of the
  identifiers: about 70 bytes per individual removed for identifiers of
  up to 15 characters, and about 115 bytes for identifiers of 34
  characters.

Both add little to the time taken by `rel_add()` and `rel_remove()`,
`"id"` a little more than `"fingerprint"` since it handles whole
identifiers. The cost only becomes noticeable when individuals are added
and removed in large numbers without much else being computed.

## Metapopulations

To follow several patches connected by migration, use a single store for
the whole metapopulation and do not remove the individuals that move to
another patch: remove them only when they die or leave the
metapopulation for good. Keep track of who is in which patch in your
simulation, and extract the matrix of a patch with
`rel_get(rel, ids_of_the_patch)`. Migrants then keep their exact
relatedness with everybody, including the relatives they left behind and
the offspring they produce elsewhere. The time per birth and the memory
depend on the number of individuals alive over all the patches (see the
section Memory management).

Creating the store with `rel_new(remember_removed = "fingerprint")` or
`rel_new(remember_removed = "id")` is recommended here: with individuals
moving between patches, removing one that left a patch and adding it
again when it comes back is an easy mistake to make, and it would go
unnoticed.

## Memory management

The store allocates its matrix by itself: there is nothing to set. The
matrix has one row and one column per slot, hence uses
`8 * rel_capacity(rel)^2` bytes. When all the slots are taken, the
number of slots grows by half (or more if a single call to `rel_add()`
requires it), which costs a one-off copy of the matrix. Beyond the
initial 16 slots, the number of slots therefore stays within about 1.5
times the largest number of individuals the store has held at once. It
never shrinks.

R does not see this memory:
[`utils::object.size()`](https://rdrr.io/r/utils/object.size.html) or
[`lobstr::obj_size()`](https://lobstr.r-lib.org/reference/obj_size.html)
only report the small R side of a store. Also, each time the matrix is
enlarged, the old and the new matrix coexist for a moment, and the
memory of the old one is not always returned to the operating system at
once: in our benchmark, a store that grew in small steps made the R
process use about 20% more memory than its matrix. All the figures given
here can be checked on your machine with the script
`system.file("benchmarks/memory_benchmark.R", package = "fastrel")`.

## See also

[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
to compute the whole matrix from a pedigree;
[`rel_save()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
and
[`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
to save a store to disk.

## Examples

``` r
rel <- rel_new()

## founders: unknown parents
rel_add(rel, data.frame(id = c("A", "B", "C", "D"), dad = NA, mom = NA))

## births: the columns are found by their names, in any order
rel_add(rel, data.frame(sire = c("B", "D"), dam = c("A", "C"), animal = c("E", "F")))
rel_add(rel, data.frame(id = "G", dad = "F", mom = "E"))

## deaths (or emigration): free the slots
rel_remove(rel, c("A", "B", "C", "D"))
rel_size(rel)
#> [1] 3
rel_ids(rel)
#> [1] "E" "F" "G"

## a mother-son mating
rel_add(rel, data.frame(id = "H", dad = "G", mom = "E"))

A <- rel_get(rel)
A
#>      E    F    G    H
#> E 1.00 0.00 0.50 0.75
#> F 0.00 1.00 0.50 0.25
#> G 0.50 0.50 1.00 0.75
#> H 0.75 0.25 0.75 1.25
diag(A) - 1   # inbreeding coefficients
#>    E    F    G    H 
#> 0.00 0.00 0.00 0.25 
rel_get(rel, metric = "kinship")   # other measures: see ?relatedness_metrics
#>       E     F     G     H
#> E 0.500 0.000 0.250 0.375
#> F 0.000 0.500 0.250 0.125
#> G 0.250 0.250 0.500 0.375
#> H 0.375 0.125 0.375 0.625

## an independent copy, e.g. to explore two scenarios from the same state
rel2 <- rel_copy(rel)
rel_remove(rel2, "H")
rel_ids(rel)    # unchanged
#> [1] "E" "F" "G" "H"
rel_ids(rel2)
#> [1] "E" "F" "G"

## two patches: one store, and the patch of each individual kept aside
rel <- rel_new(remember_removed = "id")
rel_add(rel, data.frame(id = c("a1", "a2", "b1", "b2"), dad = NA, mom = NA))
patch <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
rel_add(rel, data.frame(id = "a3", dad = "a2", mom = "a1"))   # born in A
patch["a3"] <- "A"
patch["a3"] <- "B"                                            # a3 moves to B
rel_add(rel, data.frame(id = "b3", dad = "a3", mom = "b1"))   # and breeds there
patch["b3"] <- "B"
rel_get(rel, names(patch)[patch == "B"])
#>     b1 b2  a3  b3
#> b1 1.0  0 0.0 0.5
#> b2 0.0  1 0.0 0.0
#> a3 0.0  0 1.0 0.5
#> b3 0.5  0 0.5 1.0
rel_get(rel)["b3", "a1"]   # b3 is related to its grandparent in patch A
#> [1] 0.25

## with remember_removed = "id", a removed individual cannot come back
rel_remove(rel, "a3")
try(rel_add(rel, data.frame(id = "a3", dad = NA, mom = NA)))
#> Error : some focal IDs were removed from the store earlier and cannot be added again: a3
```
