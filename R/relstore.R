#' Incremental store of the relatedness matrix
#'
#' A *store* is an object that keeps the additive relationship matrix
#' \eqn{A = 2\phi}{A = 2 * phi} among the individuals currently in it, and
#' that is updated as individuals are added and removed. Because
#' relatedness between existing individuals never changes when new individuals
#' are born, a newborn's row is computed once from its parents' rows, at a cost
#' proportional to the number of individuals the store has held at once instead
#' of recomputing the whole matrix. Individuals who can no longer become parents
#' (dead, emigrated, ...) should be removed so that their slots are recycled.
#' The cost per birth and the memory then depend on the number of individuals
#' in the store, not on the length of the pedigree.
#'
#' The functions work as follows:
#' * `rel_new()` creates an empty store.
#' * `rel_add()` adds individuals, given by a data frame with columns for the
#'   focal individual, its father and its mother.
#' * `rel_remove()` removes individuals.
#' * `rel_get()` extracts the relatedness matrix among chosen individuals, as
#'   any of the measures described in [relatedness_metrics].
#' * `rel_copy()` returns an independent copy of a store.
#' * `rel_ids()` lists the individuals in the store.
#' * `rel_size()` and `rel_capacity()` report the number of individuals in the
#'   store and the number of slots allocated.
#'
#' @section Vocabulary:
#' *Store* and *slot* are everyday words of computer science, not terms of R or
#' of C++.
#' * A **store** is a container that keeps data in memory between calls, here
#'   the relatedness matrix. You create it once and then pass it to the other
#'   functions, which read it or update it.
#' * A **slot** is one position in the matrix of the store, i.e. one row and
#'   the matching column, reserved for one individual. An individual occupies a
#'   slot from the moment it is added until it is removed; the slot is then
#'   free and is given to an individual added later. The matrix is allocated
#'   with more slots than needed so that it does not have to be enlarged at
#'   each birth: its number of slots is the *capacity* of the store.
#'
#' @section What a store is made of:
#' The object returned by `rel_new()` is an R environment (the kind of object
#' created by [new.env()]) of class `"relatedness_store"`. It holds the
#' identifiers of the individuals, the slot of each of them, and a pointer to a
#' C++ object which holds the matrix itself. Two consequences follow, which
#' differ from the usual behaviour of R objects:
#' * Like any R environment, a store is modified in place and is not copied
#'   when it is assigned or passed to a function: `rel_add()` and `rel_remove()`
#'   change the store they are given, and after `rel2 <- rel` both names refer
#'   to the same store. Use `rel_copy()` to get a store that can be modified
#'   without changing the original.
#' * The matrix lives in memory managed by C++, outside the memory managed by
#'   R. It is therefore not written by [save()] or [saveRDS()], and a store
#'   read back from such a file is unusable: use [rel_save()] and [rel_load()]
#'   instead. For the same reason [utils::object.size()] does not see it: see
#'   the section Memory management.
#'
#' The content of the environment is internal: use the functions above rather
#' than reading or modifying it directly.
#'
#' @section What to keep in the store:
#' Keep an individual for as long as it can still become a parent, or for as long
#' as you want to query its relatedness. A removed individual cannot be used as a
#' parent later (`rel_add()` raises an error), and a newborn gets no relatedness
#' value with individuals removed before its birth. Relatedness values among the
#' remaining individuals are not affected by removals.
#'
#' Removal is final. The relatedness between a removed individual and those
#' born afterwards is not tracked, so if it is added again it comes back as
#' unrelated to everybody, which is wrong for its relatives. By default
#' nothing prevents this: it is up to you never to add again an individual
#' that was removed. With `rel_new(remember_removed = "fingerprint")` or
#' `rel_new(remember_removed = "id")`, `rel_add()` raises an error instead (see
#' the section Removed individuals).
#'
#' @section Removed individuals:
#' With `remember_removed` set to `"fingerprint"` or `"id"` (rather than
#' `"none"`, the default), the store keeps a trace of every individual removed
#' and `rel_add()` refuses to add it again. This protects against a silent error, at a cost in memory that grows
#' slowly with the length of the pedigree, instead of depending only on the
#' number of individuals in the store. This is why it is off by default. Two
#' kinds of trace are offered:
#' * `"fingerprint"` keeps a 64-bit fingerprint of each
#'   identifier, about 30 bytes per individual removed (30 MB per million). A
#'   new identifier can be mistaken for a removed one, but with a negligible
#'   probability: the number of individuals removed divided by `2^64`, i.e.
#'   less than 1 in 10 million millions per addition after a million removals.
#' * `"id"` keeps the identifiers themselves, so it can never be wrong, but it
#'   takes more memory, which increases with the length of the identifiers:
#'   about 70 bytes per individual removed for identifiers of up to 15
#'   characters, and about 115 bytes for identifiers of 34 characters.
#'
#' Both add little to the time taken by `rel_add()` and `rel_remove()`,
#' `"id"` a little more than `"fingerprint"` since it handles whole
#' identifiers. The cost only becomes noticeable when individuals are added
#' and removed in large numbers without much else being computed.
#'
#' @section Metapopulations:
#' To follow several patches connected by migration, use a single store for
#' the whole metapopulation and do not remove the individuals that move to
#' another patch: remove them only when they die or leave the metapopulation
#' for good. Keep track of who is in which patch in your simulation, and
#' extract the matrix of a patch with `rel_get(rel, ids_of_the_patch)`.
#' Migrants then keep their exact relatedness with everybody, including the
#' relatives they left behind and the offspring they produce elsewhere. The
#' time per birth and the memory depend on the number of individuals alive
#' over all the patches (see the section Memory management).
#'
#' Creating the store with `rel_new(remember_removed = "fingerprint")` or
#' `rel_new(remember_removed = "id")` is recommended here: with individuals
#' moving between patches, removing one that left a patch and adding it again
#' when it comes back is an easy mistake to make, and it would go unnoticed.
#'
#' @section Memory management:
#' The store allocates its matrix by itself: there is nothing to set. The
#' matrix has one row and one column per slot, hence uses `8 * rel_capacity(rel)^2`
#' bytes. When all the slots are taken, the number of slots grows by half (or
#' more if a single call to `rel_add()` requires it), which costs a one-off
#' copy of the matrix. Beyond the initial 16 slots, the number of slots therefore
#' stays within about 1.5 times the largest number of individuals the store has
#' held at once. It never shrinks.
#'
#' R does not see this memory: [utils::object.size()] or `lobstr::obj_size()`
#' only report the small R side of a store. Also, each time the matrix is
#' enlarged, the old and the new matrix coexist for a moment, and the memory of
#' the old one is not always returned to the operating system at once: in our
#' benchmark, a store that grew in small steps made the R process use about
#' 20% more memory than its matrix. All the figures given here can be checked
#' on your machine with the script
#' `system.file("benchmarks/memory_benchmark.R", package = "fastrel")`.
#'
#' @param remember_removed Should the store keep a trace of the individuals
#'   removed, so that `rel_add()` refuses to add them again? `"none"` (the
#'   default) for no trace, `"fingerprint"` for a compact trace that is wrong
#'   with a negligible probability, or `"id"` for the identifiers themselves,
#'   which takes more memory but is never wrong (see the section Removed
#'   individuals).
#' @param rel A store created by `rel_new()`.
#' @param ped A data frame (or a tibble) with one row per individual to add and three columns
#'   giving the identifier of the individual, of its father and of its mother
#'   (`NA` or 0 for unknown parents). The columns are found by their names (e.g.
#'   `id`, `dad`, `mom` or `Animal`, `Sire`, `Dam`), whatever their order, and
#'   other columns are ignored: see the section Columns of `ped` in
#'   [fastrel_ped()]. Only individuals not yet in the store can be added.
#'   Known parents must be in the store or in `ped` itself: any other
#'   identifier raises an error rather than being treated as unknown, since it
#'   usually reveals a mistake (a typo, a parent never added or removed too
#'   early). 0 is taken as an unknown parent unless an individual is called 0,
#'   in which case it is refused as ambiguous. The rows of `ped` may be in any
#'   order.
#' @param ids Character vector (or coercible to it) of individual identifiers.
#'   A data frame with a single column is taken as that column.
#'   For `rel_get()`, the default is every individual in the store.
#' @param metric The measure of relatedness returned by `rel_get()`:
#'   `"additive"` (the default), `"kinship"`, `"wright"` (or `"Wright"`) or
#'   `"hamilton"` (or `"Hamilton"`). Abbreviations are accepted. See
#'   [relatedness_metrics]. The store always holds the additive relationship,
#'   and the conversion is done on the matrix extracted only.
#'
#' @return `rel_new()` returns the store. `rel_copy()` returns a new store,
#'   identical to `rel` (same individuals, slots, capacity and trace of the
#'   removed individuals) but independent of it. `rel_add()` and `rel_remove()` return
#'   the store invisibly. `rel_get()` returns a matrix with `ids` as dimnames:
#'   with the default `metric = "additive"`, it is symmetric with diagonal
#'   `1 + F` (see [relatedness_metrics] for the other measures). `rel_ids()`
#'   returns a character vector, in order of insertion. `rel_size()` and
#'   `rel_capacity()` return integers.
#'
#' @seealso [fastrel_ped()] to compute the whole matrix from a pedigree;
#'   [rel_save()] and [rel_load()] to save a store to disk.
#'
#' @examples
#' rel <- rel_new()
#'
#' ## founders: unknown parents
#' rel_add(rel, data.frame(id = c("A", "B", "C", "D"), dad = NA, mom = NA))
#'
#' ## births: the columns are found by their names, in any order
#' rel_add(rel, data.frame(sire = c("B", "D"), dam = c("A", "C"), animal = c("E", "F")))
#' rel_add(rel, data.frame(id = "G", dad = "F", mom = "E"))
#'
#' ## deaths (or emigration): free the slots
#' rel_remove(rel, c("A", "B", "C", "D"))
#' rel_size(rel)
#' rel_ids(rel)
#'
#' ## a mother-son mating
#' rel_add(rel, data.frame(id = "H", dad = "G", mom = "E"))
#'
#' A <- rel_get(rel)
#' A
#' diag(A) - 1   # inbreeding coefficients
#' rel_get(rel, metric = "kinship")   # other measures: see ?relatedness_metrics
#'
#' ## an independent copy, e.g. to explore two scenarios from the same state
#' rel2 <- rel_copy(rel)
#' rel_remove(rel2, "H")
#' rel_ids(rel)    # unchanged
#' rel_ids(rel2)
#'
#' ## two patches: one store, and the patch of each individual kept aside
#' rel <- rel_new(remember_removed = "id")
#' rel_add(rel, data.frame(id = c("a1", "a2", "b1", "b2"), dad = NA, mom = NA))
#' patch <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
#' rel_add(rel, data.frame(id = "a3", dad = "a2", mom = "a1"))   # born in A
#' patch["a3"] <- "A"
#' patch["a3"] <- "B"                                            # a3 moves to B
#' rel_add(rel, data.frame(id = "b3", dad = "a3", mom = "b1"))   # and breeds there
#' patch["b3"] <- "B"
#' rel_get(rel, names(patch)[patch == "B"])
#' rel_get(rel)["b3", "a1"]   # b3 is related to its grandparent in patch A
#'
#' ## with remember_removed = "id", a removed individual cannot come back
#' rel_remove(rel, "a3")
#' try(rel_add(rel, data.frame(id = "a3", dad = NA, mom = NA)))
#' @name relstore
NULL

#' @rdname relstore
#' @export
rel_new <- function(remember_removed = c("none", "fingerprint", "id")) {
  remember_removed <- match.arg(remember_removed)
  ## code for the C++ side: 0 = no trace, 1 = fingerprints, 2 = identifiers
  mode <- match(remember_removed, c("none", "fingerprint", "id")) - 1L
  self <- new.env()
  ## 16: initial number of slots; grows as needed
  self$ks    <- new(RelStore, 16L, mode)
  self$ptr   <- self$ks$.pointer   # kept aside: quick to test in check_store()
  self$ids   <- character(0)   # IDs currently in the store
  self$slots <- integer(0)     # their slots, same order
  class(self) <- "relatedness_store"
  self
}

#' @rdname relstore
#' @export
rel_copy <- function(rel) {
  check_store(rel)
  out <- rel_new(c("none", "fingerprint", "id")[rel$ks$remember() + 1L])
  out$ks$copy_from(rel$ptr)
  out$ids   <- rel$ids
  out$slots <- rel$slots
  out
}

#' @rdname relstore
#' @export
rel_add <- function(rel, ped) {
  check_store(rel)
  cols <- ped_columns(ped)                        # found by their names
  id  <- as.character(ped[[cols[1L]]])
  dad <- as.character(ped[[cols[2L]]])
  mom <- as.character(ped[[cols[3L]]])
  n <- length(id)
  if (n == 0L) return(invisible(rel))             # nothing to add

  if (anyNA(id))                stop("focal IDs cannot be NA")
  if (anyDuplicated(id))        stop("duplicated focal IDs in ped")
  mom <- zero_as_unknown(mom, c(id, rel$ids))    # 0 = unknown parent
  dad <- zero_as_unknown(dad, c(id, rel$ids))

  ## a single lookup in the store for the focal individuals and their parents
  k <- match(c(id, mom, dad), rel$ids)
  if (!all(is.na(k[seq_len(n)]))) stop("some focal IDs are already in the store")
  sm <- rel$slots[k[n + seq_len(n)]]              # slots of the parents; NA if
  sd <- rel$slots[k[2L * n + seq_len(n)]]         # unknown or not in the store
  m_new <- is.na(sm) & !is.na(mom)                # known parents that are not in
  d_new <- is.na(sd) & !is.na(dad)                # the store: must be in ped

  rel$ks$prepare(id)       # refuse removed individuals; make room once for them all
  if (!any(m_new) && !any(d_new)) {
    ## Usual case: all the parents are already in the store.
    ord <- seq_len(n)
    new <- rel$ks$add(sm, sd)
  } else {
    pm <- match(mom, id)                          # row of ped holding the parent
    pd <- match(dad, id)
    missing_par <- unique(c(mom[m_new & is.na(pm)], dad[d_new & is.na(pd)]))
    if (length(missing_par))
      stop("parents not found in the store or in ped: ",
           paste(head(missing_par, 5), collapse = ", "))

    ## Parents must be inserted before their children: insert, round after
    ## round, the individuals whose parents are all in the store.
    new  <- rep(NA_integer_, n)                   # slots given to the rows of ped
    done <- logical(n)
    ord  <- integer(0)                            # rows of ped, in insertion order
    while (length(ord) < n) {
      todo  <- which(!done)
      ready <- todo[(!m_new[todo] | done[pm[todo]]) & (!d_new[todo] | done[pd[todo]])]
      if (!length(ready)) {
        rel$ks$remove(new[ord], character(0))     # leave the store unchanged
        stop("impossible pedigree: cycle among the new individuals")
      }
      s1 <- sm[ready]; i <- m_new[ready]; s1[i] <- new[pm[ready[i]]]
      s2 <- sd[ready]; i <- d_new[ready]; s2[i] <- new[pd[ready[i]]]
      new[ready]  <- rel$ks$add(s1, s2)
      done[ready] <- TRUE
      ord <- c(ord, ready)
    }
  }
  rel$ids   <- c(rel$ids, id[ord])
  rel$slots <- c(rel$slots, new[ord])
  invisible(rel)
}

#' @rdname relstore
#' @export
rel_remove <- function(rel, ids) {
  check_store(rel)
  ids <- unique(as.character(id_vector(ids, "ids")))
  if (!length(ids)) return(invisible(rel))          # nothing to remove
  i <- match(ids, rel$ids)
  if (anyNA(i)) stop("unknown IDs: ", paste(head(ids[is.na(i)], 5), collapse = ", "))
  rel$ks$remove(rel$slots[i], ids)                # ids: they cannot be added again
  rel$ids   <- rel$ids[-i]
  rel$slots <- rel$slots[-i]
  invisible(rel)
}

#' @rdname relstore
#' @export
rel_get <- function(rel, ids = rel_ids(rel),
                    metric = c("additive", "kinship", "wright", "hamilton", "Wright", "Hamilton")) {
  metric <- match.arg(metric)
  code <- metric_code(metric)
  check_store(rel)
  ids <- as.character(id_vector(ids, "ids"))
  i <- match(ids, rel$ids)
  if (anyNA(i)) stop("unknown IDs: ", paste(head(ids[is.na(i)], 5), collapse = ", "))
  A <- rel$ks$get(rel$slots[i], code)
  dimnames(A) <- list(ids, ids)
  A
}

#' @rdname relstore
#' @export
rel_ids <- function(rel) {
  check_store(rel)
  rel$ids
}

#' @rdname relstore
#' @export
rel_size <- function(rel) {
  check_store(rel)
  rel$ks$size()
}

#' @rdname relstore
#' @export
rel_capacity <- function(rel) {
  check_store(rel)
  rel$ks$capacity()
}

#' @param x A store created by `rel_new()`.
#' @param ... Ignored.
#' @rdname relstore
#' @export
print.relatedness_store <- function(x, ...) {
  check_store(x)
  cat("<relatedness_store>", length(x$ids), "individuals,",
      x$ks$capacity(), "slots\n")
  invisible(x)
}

check_store <- function(rel) {
  if (!inherits(rel, "relatedness_store"))
    stop("`rel` must be a store created by rel_new()")
  ## a store read back with readRDS() or load() has lost its C++ content
  if (identical(rel$ptr, nil_pointer))
    stop("this store is no longer valid: a store cannot be saved with save() ",
         "or saveRDS(), use rel_save() and rel_load() instead")
  invisible(TRUE)
}

## what the pointer to the C++ object becomes once saved and read back
nil_pointer <- methods::new("externalptr")

#' Save and restore a store
#'
#' A store lives in C++ memory and cannot be saved with [saveRDS()] directly.
#' `rel_save()` extracts its content (the identifiers and the relatedness matrix
#' among them) and `rel_load()` rebuilds a working store from
#' it, so that a simulation can be checkpointed and resumed.
#'
#' The saved object contains the full dense relatedness matrix among the
#' individuals in the store, so its size is `8 * n^2` bytes for `n` individuals,
#' plus the trace of the individuals removed so far if the store keeps one: 8
#' bytes per individual removed with `remember_removed = "fingerprint"`, or
#' their identifiers with `remember_removed = "id"`.
#' After loading, individuals have new internal slots, which is invisible when
#' you work through the identifiers. The restored store behaves exactly like the
#' original one: you can keep adding and removing individuals.
#'
#' The saved object also records the version of fastrel that saved it
#' (`fastrel_version`) and the format of the object (`format`, which changes
#' only when the structure of saved objects changes). `rel_load()` reads
#' objects in the format of the installed version of fastrel. If an object is
#' in another format, it says which version of fastrel saved it.
#'
#' @param rel A store created by [rel_new()].
#' @param file Path of an `.rds` file. If `NULL` (default), `rel_save()` returns
#'   the saved object instead of writing it.
#' @param x For `rel_load()`, either the path of an `.rds` file written by
#'   `rel_save()` or the object returned by `rel_save(rel)`.
#'
#' @return `rel_save()` returns the saved object (a list of class
#'   `"relatedness_store_saved"`), or, when `file` is given, the file path
#'   invisibly. `rel_load()` returns a new store.
#'
#' @examples
#' rel <- rel_new()
#' rel_add(rel, data.frame(id = c("A", "B"), dad = NA, mom = NA))
#' rel_add(rel, data.frame(id = "C", dad = "B", mom = "A"))
#'
#' f <- tempfile(fileext = ".rds")
#' rel_save(rel, f)
#' rel2 <- rel_load(f)
#' all.equal(rel_get(rel), rel_get(rel2))
#'
#' ## the restored store can be used as usual
#' rel_add(rel2, data.frame(id = "D", dad = "A", mom = "C"))
#' @seealso [rel_new()]
#' @name relstore_io
NULL

#' @rdname relstore_io
#' @export
rel_save <- function(rel, file = NULL) {
  check_store(rel)
  out <- structure(
    list(ids = rel$ids, A = rel_get(rel),
         remember_removed = c("none", "fingerprint", "id")[rel$ks$remember() + 1L],
         removed = rel$ks$gone(),
         format = saved_format,
         fastrel_version = as.character(utils::packageVersion("fastrel"))),
    class = "relatedness_store_saved")
  if (is.null(file)) return(out)
  saveRDS(out, file)
  invisible(file)
}

#' @rdname relstore_io
#' @export
rel_load <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- readRDS(x)
  if (!inherits(x, "relatedness_store_saved") || !is.list(x))
    stop("`x` must be a file or an object produced by rel_save()")
  ## the format is checked first: an object in another format may have other
  ## fields
  if (!identical(x$format, saved_format))
    stop("this object was saved ",
         if (is.character(x$fastrel_version)) paste0("by fastrel ", x$fastrel_version, " "),
         "in a format that this version of fastrel (",
         as.character(utils::packageVersion("fastrel")), ") cannot read",
         if (is.numeric(x$format) && isTRUE(x$format > saved_format))
           "; update fastrel to read it")
  fields <- c("ids", "A", "remember_removed", "removed", "format", "fastrel_version")
  if (!all(fields %in% names(x)))
    stop("`x` must be a file or an object produced by rel_save()")
  k <- length(x$ids)
  if (!is.matrix(x$A) || nrow(x$A) != k || ncol(x$A) != k)
    stop("the saved matrix does not match the saved identifiers")
  if (anyNA(x$ids) || anyDuplicated(x$ids))
    stop("the saved identifiers must be unique and not NA")
  remember <- x$remember_removed
  if ((identical(remember, "fingerprint") && !is.raw(x$removed)) ||
      (identical(remember, "id") && !is.character(x$removed)))
    stop("the saved list of removed individuals is not valid")

  rel <- rel_new(remember_removed = remember)
  if (remember != "none") rel$ks$load_gone(x$removed)
  if (k > 0L) {
    rel$slots <- rel$ks$restore(unname(x$A))
    rel$ids   <- as.character(x$ids)
  }
  rel
}

## Format of the objects written by rel_save(). It changes only when their
## structure changes; rel_load() reads this format (and, in future versions,
## may convert the earlier ones).
saved_format <- 1L
