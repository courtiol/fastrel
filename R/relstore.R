#' Incremental store of the relatedness matrix
#'
#' A persistent store holds the additive relationship matrix
#' \eqn{A = 2\phi}{A = 2 * phi} among the individuals currently in it. Because
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
#' * `rel_add()` adds individuals, given by a data frame with the columns
#'   focal, mother and father.
#' * `rel_remove()` removes individuals.
#' * `rel_get()` extracts the relatedness matrix among chosen individuals.
#' * `rel_ids()` lists the individuals in the store.
#' * `rel_size()` and `rel_capacity()` report the number of individuals in the
#'   store and the number of slots allocated.
#'
#' A store is an environment: `rel_add()` and `rel_remove()` modify it in place,
#' and copying the object does not copy the data. A store cannot be saved to disk
#' with [save()] or [saveRDS()] and restored in another session: use
#' [rel_save()] and [rel_load()] instead. The content of the environment is
#' internal: use the functions above rather than reading or modifying it.
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
#' that was removed. With `rel_new(remember_removed = TRUE)`, `rel_add()`
#' raises an error instead (see the section Removed individuals).
#'
#' @section Removed individuals:
#' With `remember_removed = TRUE`, the store keeps a trace of every individual
#' removed and `rel_add()` refuses to add it again. This protects against a
#' silent error, at the cost of 11 to 23 bytes per individual removed: the
#' memory used then grows slowly with the length of the pedigree (11 to 23 MB
#' per million individuals removed), instead of depending only on the number
#' of individuals in the store. This is why it is off by default. The time
#' taken by `rel_add()` and `rel_remove()` is not affected.
#'
#' The trace is a 64-bit fingerprint of each identifier, not the identifier
#' itself. A new identifier can thus be mistaken for a removed one, but with a
#' negligible probability: the number of individuals removed divided by
#' `2^64`, i.e. less than 1 in 10 million millions per addition after a million
#' removals.
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
#' over all the patches (see the section Memory).
#'
#' Creating the store with `rel_new(remember_removed = TRUE)` is recommended
#' here: with individuals moving between patches, removing one that left a
#' patch and adding it again when it comes back is an easy mistake to make,
#' and it would go unnoticed.
#'
#' @section Memory:
#' The store allocates its matrix by itself: there is nothing to set. The
#' matrix has one row and one column per slot, hence uses `8 * rel_capacity(rel)^2`
#' bytes. When all the slots are taken, the number of slots grows by half (or
#' more if a single call to `rel_add()` requires it), which costs a one-off
#' copy of the matrix. Beyond the initial 16 slots, the number of slots therefore
#' stays within about 1.5 times the largest number of individuals the store has
#' held at once. It never shrinks.
#'
#' @param remember_removed Should the store keep a trace of the individuals
#'   removed, so that `rel_add()` refuses to add them again? `FALSE` by default
#'   (see the section Removed individuals).
#' @param rel A store created by `rel_new()`.
#' @param ped A data frame whose first three columns are the identifier of the
#'   focal individual, of its mother and of its father (`NA` for unknown
#'   parents). Only individuals not yet in the store can be added. Parents must
#'   be in the store or in `ped` itself; the rows of `ped` may be in any order.
#' @param ids Character vector (or coercible to it) of individual identifiers.
#'   For `rel_get()`, the default is every individual in the store.
#'
#' @return `rel_new()` returns the store. `rel_add()` and `rel_remove()` return
#'   the store invisibly. `rel_get()` returns a symmetric matrix with `ids` as
#'   dimnames (diagonal `1 + F`). `rel_ids()` returns a character vector, in
#'   order of insertion. `rel_size()` and `rel_capacity()` return integers.
#'
#' @seealso [relatedness_fast()] to compute the whole matrix from a pedigree;
#'   [rel_save()] and [rel_load()] to save a store to disk.
#'
#' @examples
#' rel <- rel_new()
#'
#' ## founders: unknown parents
#' rel_add(rel, data.frame(id = c("A", "B", "C", "D"), mom = NA, dad = NA))
#'
#' ## births: columns are focal, mother, father
#' rel_add(rel, data.frame(id = c("E", "F"), mom = c("A", "C"), dad = c("B", "D")))
#' rel_add(rel, data.frame(id = "G", mom = "E", dad = "F"))
#'
#' ## deaths (or emigration): free the slots
#' rel_remove(rel, c("A", "B", "C", "D"))
#' rel_size(rel)
#' rel_ids(rel)
#'
#' ## a mother-son mating
#' rel_add(rel, data.frame(id = "H", mom = "E", dad = "G"))
#'
#' A <- rel_get(rel)
#' A
#' diag(A) - 1   # inbreeding coefficients
#'
#' ## two patches: one store, and the patch of each individual kept aside
#' rel <- rel_new(remember_removed = TRUE)
#' rel_add(rel, data.frame(id = c("a1", "a2", "b1", "b2"), mom = NA, dad = NA))
#' patch <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
#' rel_add(rel, data.frame(id = "a3", mom = "a1", dad = "a2"))   # born in A
#' patch["a3"] <- "A"
#' patch["a3"] <- "B"                                            # a3 moves to B
#' rel_add(rel, data.frame(id = "b3", mom = "b1", dad = "a3"))   # and breeds there
#' patch["b3"] <- "B"
#' rel_get(rel, names(patch)[patch == "B"])
#' rel_get(rel)["b3", "a1"]   # b3 is related to its grandparent in patch A
#'
#' ## with remember_removed = TRUE, a removed individual cannot come back
#' rel_remove(rel, "a3")
#' try(rel_add(rel, data.frame(id = "a3", mom = NA, dad = NA)))
#' @name relstore
NULL

#' @rdname relstore
#' @export
rel_new <- function(remember_removed = FALSE) {
  if (!is.logical(remember_removed) || length(remember_removed) != 1L ||
      is.na(remember_removed))
    stop("`remember_removed` must be TRUE or FALSE")
  self <- new.env()
  ## 16: initial number of slots; grows as needed
  self$ks    <- new(RelStore, 16L, remember_removed)
  self$ptr   <- self$ks$.pointer   # kept aside: quick to test in check_store()
  self$ids   <- character(0)   # IDs currently in the store
  self$slots <- integer(0)     # their slots, same order
  class(self) <- "relatedness_store"
  self
}

#' @rdname relstore
#' @export
rel_add <- function(rel, ped) {
  check_store(rel)
  if (!is.data.frame(ped) || ncol(ped) < 3L)
    stop("`ped` must be a data frame with columns focal, mother and father")
  id  <- as.character(ped[[1]])
  mom <- as.character(ped[[2]])
  dad <- as.character(ped[[3]])
  n <- length(id)
  if (n == 0L) return(invisible(rel))             # nothing to add

  if (anyNA(id))                stop("focal IDs cannot be NA")
  if (anyDuplicated(id))        stop("duplicated focal IDs in ped")

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
  ids <- unique(as.character(ids))
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
rel_get <- function(rel, ids = rel$ids) {
  check_store(rel)
  ids <- as.character(ids)
  i <- match(ids, rel$ids)
  if (anyNA(i)) stop("unknown IDs: ", paste(head(ids[is.na(i)], 5), collapse = ", "))
  A <- rel$ks$get(rel$slots[i])
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
#' plus, for a store created with `remember_removed = TRUE`, 8 bytes per
#' individual removed so far.
#' After loading, individuals have new internal slots, which is invisible when
#' you work through the identifiers. The restored store behaves exactly like the
#' original one: you can keep adding and removing individuals.
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
#' rel_add(rel, data.frame(id = c("A", "B"), mom = NA, dad = NA))
#' rel_add(rel, data.frame(id = "C", mom = "A", dad = "B"))
#'
#' f <- tempfile(fileext = ".rds")
#' rel_save(rel, f)
#' rel2 <- rel_load(f)
#' all.equal(rel_get(rel), rel_get(rel2))
#'
#' ## the restored store can be used as usual
#' rel_add(rel2, data.frame(id = "D", mom = "C", dad = "A"))
#' @seealso [rel_new()]
#' @name relstore_io
NULL

#' @rdname relstore_io
#' @export
rel_save <- function(rel, file = NULL) {
  check_store(rel)
  out <- structure(
    list(ids = rel$ids, A = rel_get(rel), remember_removed = rel$ks$remember(),
         removed = rel$ks$gone(), format = 1L),
    class = "relatedness_store_saved")
  if (is.null(file)) return(out)
  saveRDS(out, file)
  invisible(file)
}

#' @rdname relstore_io
#' @export
rel_load <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- readRDS(x)
  if (!inherits(x, "relatedness_store_saved") ||
      !all(c("ids", "A") %in% names(x)))
    stop("`x` must be a file or an object produced by rel_save()")
  if (!is.null(x$format) && !identical(as.integer(x$format), 1L))
    stop("this file was written in a format (", x$format, ") that this ",
         "version of fastrel cannot read")
  k <- length(x$ids)
  if (!is.matrix(x$A) || nrow(x$A) != k || ncol(x$A) != k)
    stop("the saved matrix does not match the saved identifiers")
  if (anyNA(x$ids) || anyDuplicated(x$ids))
    stop("the saved identifiers must be unique and not NA")

  remember <- isTRUE(x$remember_removed)
  if (remember && !is.raw(x$removed))
    stop("the saved list of removed individuals is not valid")

  rel <- rel_new(remember_removed = remember)
  if (remember) rel$ks$load_gone(x$removed)
  if (k > 0L) {
    rel$slots <- rel$ks$restore(unname(x$A))
    rel$ids   <- as.character(x$ids)
  }
  rel
}
