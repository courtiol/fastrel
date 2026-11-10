#' Incremental store of the relatedness matrix
#'
#' A persistent store holds the additive relationship matrix
#' \eqn{A = 2\phi}{A = 2 * phi} among the individuals currently in it. Because
#' relatedness between existing individuals never changes when new individuals
#' are born, a newborn's row is computed once from its parents' rows, at a cost
#' proportional to the capacity of the store instead of recomputing the whole
#' matrix. Individuals who can no longer become parents (dead, emigrated, ...)
#' should be removed so that their slots are recycled. The cost per birth and
#' the memory (`8 * capacity^2` bytes) then depend on the number of individuals
#' in the store, not on the length of the pedigree.
#'
#' The functions work as follows:
#' * `rel_new()` creates an empty store.
#' * `rel_add()` adds individuals, given by a data frame with the columns
#'   focal, mother and father.
#' * `rel_remove()` removes individuals.
#' * `rel_get()` extracts the relatedness matrix among chosen individuals.
#' * `rel_size()` and `rel_capacity()` report the number of individuals in the
#'   store and the number of slots allocated.
#'
#' A store is an environment: `rel_add()` and `rel_remove()` modify it in place,
#' and copying the object does not copy the data. A store cannot be saved to disk
#' with [save()] or [saveRDS()] and restored in another session.
#'
#' @section What to keep in the store:
#' Keep an individual for as long as it can still become a parent, or for as long
#' as you want to query its relatedness. A removed individual cannot be used as a
#' parent later (`rel_add()` raises an error), and a newborn gets no relatedness
#' value with individuals removed before its birth. Relatedness values among the
#' remaining individuals are not affected by removals.
#'
#' @section Capacity:
#' The matrix has `capacity` rows and columns. It doubles automatically when it
#' is full, which costs a one-off copy and may overshoot. If you know roughly the
#' peak number of individuals in the store, start with about 1.2 times that
#' number. Capacity never shrinks.
#'
#' @param capacity Initial number of slots (a positive integer).
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
#'   dimnames (diagonal `1 + F`). `rel_size()` and `rel_capacity()` return
#'   integers.
#'
#' @seealso [relatedness_fast()] to compute the whole matrix from a pedigree.
#'
#' @examples
#' rel <- rel_new(16)
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
#'
#' ## a mother-son mating
#' rel_add(rel, data.frame(id = "H", mom = "E", dad = "G"))
#'
#' A <- rel_get(rel)
#' A
#' diag(A) - 1   # inbreeding coefficients
#' @name relstore
NULL

#' @rdname relstore
#' @export
rel_new <- function(capacity = 1000L) {
  if (length(capacity) != 1L || is.na(capacity) || capacity < 1)
    stop("`capacity` must be a single positive number")
  self <- new.env()
  self$ks    <- new(RelStore, as.integer(capacity))
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

  if (anyNA(id))                stop("focal IDs cannot be NA")
  if (anyDuplicated(id))        stop("duplicated focal IDs in ped")
  if (any(id %in% rel$ids))     stop("some focal IDs are already in the store")

  missing_par <- setdiff(c(mom, dad), c(NA, rel$ids, id))
  if (length(missing_par))
    stop("parents not found in the store or in ped: ",
         paste(head(missing_par, 5), collapse = ", "))

  ## Parents must be inserted before their children. Normally all parents
  ## are already in the store, so this loop runs once.
  done <- logical(length(id))
  while (!all(done)) {
    todo  <- which(!done)
    ready <- todo[(is.na(mom[todo]) | mom[todo] %in% rel$ids) &
                  (is.na(dad[todo]) | dad[todo] %in% rel$ids)]
    if (!length(ready)) stop("impossible pedigree: cycle among the new individuals")

    sm  <- rel$slots[match(mom[ready], rel$ids)]   # NA stays NA
    sd  <- rel$slots[match(dad[ready], rel$ids)]
    new <- rel$ks$add(sm, sd)

    rel$ids   <- c(rel$ids, id[ready])
    rel$slots <- c(rel$slots, new)
    done[ready] <- TRUE
  }
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
  rel$ks$remove(rel$slots[i])
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
  cat("<relatedness_store>", length(x$ids), "individuals,",
      x$ks$capacity(), "slots\n")
  invisible(x)
}

check_store <- function(rel) {
  if (!inherits(rel, "relatedness_store"))
    stop("`rel` must be a store created by rel_new()")
  invisible(TRUE)
}
