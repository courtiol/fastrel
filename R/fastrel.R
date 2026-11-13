#' Relatedness matrix from a pedigree
#'
#' Compute the additive relationship matrix \eqn{A = 2\phi}{A = 2 * phi} of
#' autosomal loci for all individuals of a pedigree, using the tabular method
#' implemented in C++. The two functions do the same thing and differ only in
#' how the pedigree is given:
#' * `fastrel_ped()` takes a data frame;
#' * `fastrel_vec()` takes three vectors, as [kinship2::kinship()] does. It is a
#'   fast replacement for `2 * kinship2::kinship()`.
#'
#' Rows of the pedigree can be in any order: parents do not need to appear
#' before their offspring. Parents that are `NA` or that are not found among
#' the individuals are treated as unknown. The algorithm uses a topological
#' sort of the pedigree, so an impossible pedigree (someone being their own
#' ancestor) raises an error.
#'
#' Note that the argument order of `fastrel_vec()` (`dadid` before `momid`)
#' follows [kinship2::kinship()].
#'
#' @section Columns of `ped`:
#' `fastrel_ped()` and [rel_add()] find the three columns they need by their
#' names, so the columns can be in any order and `ped` can contain other
#' columns, which are ignored. Names are compared without regard to
#' capitalisation, to punctuation (`mother_id`, `Mother.ID` and `motherid` are
#' the same) and to a leading or trailing `id` (`momid`, `id_mom` and `mom` are
#' the same). The names recognised are:
#' * for the focal individual: `id`, `iid`, `focal`, `ind`, `indiv`,
#'   `individual`, `animal`, `ego`, `offspring`, `child`;
#' * for its mother: `mom`, `mum`, `mother`, `dam`, `mat`, `maternal`;
#' * for its father: `dad`, `father`, `sire`, `pat`, `paternal`.
#'
#' If none of the columns has a recognised name, the first three columns are
#' taken to be the focal individual, its mother and its father, in this order.
#' If only some of the three are recognised, or if several columns could play
#' the same role, an error is raised rather than guessing.
#'
#' @param ped A data frame (or a tibble) with one row per individual and three columns giving
#'   the identifier of the individual, of its mother and of its father (`NA` if
#'   unknown). See the section Columns of `ped`.
#' @param id Vector of unique individual identifiers. For this argument and
#'   the next two, a data frame with a single column (what `x[, "id"]` returns
#'   when `x` is a tibble) is taken as that column.
#' @param dadid Vector of the same length as `id` with the identifier of each
#'   individual's father (`NA` if unknown).
#' @param momid Vector of the same length as `id` with the identifier of each
#'   individual's mother (`NA` if unknown).
#' @param ... Ignored. Present so that calls written for [kinship2::kinship()]
#'   (e.g. with `sex` or `chrtype = "autosome"`) do not fail. Only autosomes are
#'   supported: any other `chrtype` raises an error.
#'
#' @return A symmetric numeric matrix with one row and one column per
#'   individual, in the order given, and the identifiers as row and column
#'   names. The diagonal is `1 + F`, where `F` is the inbreeding coefficient.
#'
#' @seealso [rel_new()] for the incremental version.
#'
#' @examples
#' ped <- data.frame(
#'   id  = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
#'   mom = c(NA, NA, NA, NA, "A", "C", "E", "E", "E"),
#'   dad = c(NA, NA, NA, NA, "B", "D", "F", "G", "G"))
#' A <- fastrel_ped(ped)
#' A
#'
#' ## inbreeding coefficients (H and I are the offspring of a mother-son mating)
#' diag(A) - 1
#'
#' ## the same from vectors
#' all.equal(A, fastrel_vec(ped$id, dadid = ped$dad, momid = ped$mom))
#'
#' ## columns are found by their names, whatever their order
#' ped2 <- data.frame(Sire = ped$dad, Dam = ped$mom, Animal = ped$id, weight = 1:9)
#' all.equal(A, fastrel_ped(ped2))
#'
#' ## the same as twice the kinship matrix of the kinship2 package
#' if (requireNamespace("kinship2", quietly = TRUE)) {
#'   all.equal(A, 2 * kinship2::kinship(ped$id, dadid = ped$dad, momid = ped$mom))
#' }
#' @export
fastrel_ped <- function(ped) {
  cols <- ped_columns(ped)
  fastrel_vec(ped[[cols[1L]]], dadid = ped[[cols[3L]]], momid = ped[[cols[2L]]])
}

#' @rdname fastrel_ped
#' @export
fastrel_vec <- function(id, dadid, momid, ...) {
  chrtype <- list(...)[["chrtype"]]
  if (!is.null(chrtype) && !identical(tolower(chrtype), "autosome"))
    stop("only autosomes are supported (`chrtype = \"autosome\"`)")
  id    <- id_vector(id, "id")
  dadid <- id_vector(dadid, "dadid")
  momid <- id_vector(momid, "momid")
  if (anyNA(id)) stop("id values cannot be NA")
  if (anyDuplicated(id)) stop("All id values must be unique")
  n <- length(id)
  if (length(dadid) != n) stop("Invalid father id")
  if (length(momid) != n) stop("Invalid mother id")
  if (n == 0L) return(matrix(numeric(0), 0, 0))

  A <- relatedness_cpp(mom = match(momid, id), dad = match(dadid, id))
  dimnames(A) <- list(id, id)
  A
}
