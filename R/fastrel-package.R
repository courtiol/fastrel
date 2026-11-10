#' fastrel: fast additive relatedness matrices from pedigrees
#'
#' Tools to compute the additive relationship ("relatedness") matrix
#' \eqn{A = 2\phi}{A = 2 * phi}, where \eqn{\phi}{phi} is the kinship
#' coefficient, for autosomal loci.
#'
#' Two interfaces are provided:
#' * [relatedness_fast()] computes the whole matrix from a pedigree in one go
#'   (a fast replacement for `2 * kinship2::kinship()`).
#' * [rel_new()], [rel_add()], [rel_remove()] and [rel_get()] maintain the matrix
#'   incrementally, which is much cheaper when a pedigree keeps growing, as in
#'   an individual-based simulation.
#'
#' [sim_ped()] simulates a pedigree for testing and benchmarking.
#'
#' @section Scale:
#' All functions return relatedness, not kinship: \eqn{A_{ij} = 2\phi_{ij}}{A[i,j] = 2 * phi[i,j]}.
#' The diagonal is \eqn{1 + F_i}{1 + F[i]}, where \eqn{F_i}{F[i]} is the
#' inbreeding coefficient, so `diag(A) - 1` gives the inbreeding coefficients.
#' Founders (individuals with unknown parents) are assumed to be unrelated and
#' non-inbred, and an unknown parent is treated as unrelated to everybody.
#'
#' @useDynLib fastrel, .registration = TRUE
#' @importFrom Rcpp evalCpp loadModule
#' @importFrom methods new
#' @importFrom utils head globalVariables
"_PACKAGE"

## The class generator RelStore is created at load time by loadModule()
utils::globalVariables("RelStore")

Rcpp::loadModule("relstore", TRUE)
