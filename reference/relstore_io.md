# Save and restore a store

A store lives in C++ memory and cannot be saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) directly.
`rel_save()` extracts its content (the identifiers and the relatedness
matrix among them) and `rel_load()` rebuilds a working store from it, so
that a simulation can be checkpointed and resumed.

## Usage

``` r
rel_save(rel, file = NULL)

rel_load(x)
```

## Arguments

- rel:

  A store created by
  [`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md).

- file:

  Path of an `.rds` file. If `NULL` (default), `rel_save()` returns the
  saved object instead of writing it.

- x:

  For `rel_load()`, either the path of an `.rds` file written by
  `rel_save()` or the object returned by `rel_save(rel)`.

## Value

`rel_save()` returns the saved object (a list of class
`"relatedness_store_saved"`), or, when `file` is given, the file path
invisibly. `rel_load()` returns a new store.

## Details

The saved object contains the full dense relatedness matrix among the
individuals in the store, so its size is `8 * n^2` bytes for `n`
individuals, plus the trace of the individuals removed so far if the
store keeps one: 8 bytes per individual removed with
`remember_removed = "fingerprint"`, or their identifiers with
`remember_removed = "id"`. After loading, individuals have new internal
slots, which is invisible when you work through the identifiers. The
restored store behaves exactly like the original one: you can keep
adding and removing individuals.

The saved object also records the version of fastrel that saved it
(`fastrel_version`) and the format of the object (`format`, which
changes only when the structure of saved objects changes). `rel_load()`
reads objects in the format of the installed version of fastrel. If an
object is in another format, it says which version of fastrel saved it.

## See also

[`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md)

## Examples

``` r
rel <- rel_new()
rel_add(rel, data.frame(id = c("A", "B"), dad = NA, mom = NA))
rel_add(rel, data.frame(id = "C", dad = "B", mom = "A"))

f <- tempfile(fileext = ".rds")
rel_save(rel, f)
rel2 <- rel_load(f)
all.equal(rel_get(rel), rel_get(rel2))
#> [1] TRUE

## the restored store can be used as usual
rel_add(rel2, data.frame(id = "D", dad = "A", mom = "C"))
```
