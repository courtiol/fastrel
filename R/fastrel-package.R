#' fastrel: fast additive relatedness matrices from pedigrees
#'
#' Tools to compute the additive relationship ("relatedness") matrix
#' \eqn{A = 2\phi}{A = 2 * phi}, where \eqn{\phi}{phi} is the kinship
#' coefficient, for autosomal loci.
#'
#' Two interfaces are provided:
#' * [fastrel_ped()] and [fastrel_vec()] compute the whole matrix from a
#'   pedigree in one go, given as a data frame or as vectors (the latter is a
#'   fast replacement for `2 * kinship2::kinship()`).
#' * [rel_new()], [rel_add()], [rel_remove()] and [rel_get()] maintain the matrix
#'   incrementally, which is computationally much cheaper when a pedigree keeps growing, as in
#'   an individual-based simulation. [rel_ids()], [rel_size()] and
#'   [rel_capacity()] describe the content of such a store, and [rel_save()] and
#'   [rel_load()] save it to disk and restore it.
#'
#' [sim_ped()] simulates a pedigree, and [sim_pop()] a population with births,
#' deaths, emigration and immigration, for testing and benchmarking.
#'
#' @section Scale:
#' All functions return relatedness, not kinship: \eqn{A_{ij} = 2\phi_{ij}}{A[i,j] = 2 * phi[i,j]}.
#' The diagonal is \eqn{1 + F_i}{1 + F[i]}, where \eqn{F_i}{F[i]} is the
#' inbreeding coefficient, so `diag(A) - 1` gives the inbreeding coefficients.
#' Founders (individuals with unknown parents) are assumed to be unrelated and
#' non-inbred, and an unknown parent is treated as unrelated to everybody.
#'
#' @section Range of values:
#' Relatedness values lie between 0 and 2, since a kinship coefficient lies
#' between 0 and 1.
#' * The diagonal goes from 1 (not inbred) to 2 (fully inbred).
#' * The relatedness between two individuals is at most the geometric mean of
#'   their diagonal values, \eqn{A_{ij} \le \sqrt{A_{ii} A_{jj}}}{A[i,j] <= sqrt(A[i,i] * A[j,j])}.
#'   Between individuals that are not inbred it is therefore at most 1, the
#'   value for identical twins or clones; usual values are 0.5 for parent and
#'   offspring or for full siblings, and 0.25 for half siblings.
#' * Values above 1 thus only occur with inbreeding. For example, after `n`
#'   generations of selfing the inbreeding coefficient is `1 - (1/2)^n`, so the
#'   diagonal is 1.5, 1.75, 1.875, ... The maximum of 2 is approached but
#'   never reached, since founders are assumed not to be inbred.
#'
#' For a measure that cannot exceed 1 whatever the inbreeding, divide each
#' value by the geometric mean of the two diagonal values, which is what
#' [stats::cov2cor()] does: `cov2cor(A)`.
#'
#' @useDynLib fastrel, .registration = TRUE
#' @importFrom Rcpp evalCpp loadModule
#' @importFrom methods new
#' @importFrom utils head
"_PACKAGE"

## The class generator RelStore is created at load time by loadModule()
utils::globalVariables("RelStore")

Rcpp::loadModule("relstore", TRUE)
