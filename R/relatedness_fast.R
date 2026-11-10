#' Relatedness matrix from a pedigree
#'
#' Computes the additive relationship matrix \eqn{A = 2\phi}{A = 2 * phi} of
#' autosomal loci for all individuals of a pedigree, using the tabular method
#' implemented in C++. It is a fast replacement for `2 * kinship2::kinship()`.
#'
#' Rows of the pedigree can be in any order: parents do not need to appear
#' before their offspring. Parents that are `NA` or that are not found in `id`
#' are treated as unknown. The algorithm uses a topological sort of the pedigree,
#' so an impossible pedigree (someone being their own ancestor) raises an error.
#'
#' Note that the argument order (`dadid` before `momid`) follows
#' [kinship2::kinship()], whereas [rel_add()] takes focal, mother, father.
#'
#' @param id Vector of unique individual identifiers.
#' @param dadid Vector of the same length as `id` with the identifier of each
#'   individual's father (`NA` if unknown).
#' @param momid Vector of the same length as `id` with the identifier of each
#'   individual's mother (`NA` if unknown).
#' @param ... Ignored. Present so that calls written for [kinship2::kinship()]
#'   (e.g. with `sex` or `chrtype = "autosome"`) do not fail. Only autosomes are
#'   supported.
#'
#' @return A symmetric numeric matrix with `length(id)` rows and columns and
#'   `id` as row and column names. The diagonal is `1 + F`, where `F` is the
#'   inbreeding coefficient.
#'
#' @seealso [rel_new()] for the incremental version.
#'
#' @examples
#' ped <- data.frame(
#'   id  = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
#'   mom = c(NA, NA, NA, NA, "A", "C", "E", "E", "E"),
#'   dad = c(NA, NA, NA, NA, "B", "D", "F", "G", "G"))
#' A <- relatedness_fast(ped$id, dadid = ped$dad, momid = ped$mom)
#' A
#'
#' ## inbreeding coefficients (H and I are the offspring of a mother-son mating)
#' diag(A) - 1
#'
#' ## the same as twice the kinship matrix of the kinship2 package
#' if (requireNamespace("kinship2", quietly = TRUE)) {
#'   all.equal(A, 2 * kinship2::kinship(ped$id, dadid = ped$dad, momid = ped$mom))
#' }
#' @export
relatedness_fast <- function(id, dadid, momid, ...) {
  if (anyDuplicated(id)) stop("All id values must be unique")
  n <- length(id)
  if (length(dadid) != n) stop("Invalid father id")
  if (length(momid) != n) stop("Invalid mother id")
  if (n == 0L) return(matrix(numeric(0), 0, 0))

  A <- relatedness_cpp(mom = match(momid, id), dad = match(dadid, id))
  dimnames(A) <- list(id, id)
  A
}
