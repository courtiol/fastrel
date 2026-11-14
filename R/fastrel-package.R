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
#' @section Measures of relatedness:
#' By default, all functions return the additive relationship
#' \eqn{A_{ij} = 2\phi_{ij}}{A[i,j] = 2 * phi[i,j]}, where \eqn{\phi_{ij}}{phi[i,j]}
#' is the kinship coefficient. The kinship coefficient is a probability: that
#' of an allele drawn at random from \eqn{i} and an allele drawn at random from
#' \eqn{j}, at the same locus, being identical by descent. The additive
#' relationship is instead an expected count of alleles: \eqn{A_{ij}}{A[i,j]}
#' is the expected number of the two alleles of \eqn{j} that are identical by
#' descent with an allele drawn at random from \eqn{i}. This is why it lies
#' between 0 and 2, and can exceed 1 (see the section Range of values). For
#' example, an allele drawn from a parent is one of the two alleles of its
#' offspring with probability 1/2, so their additive relationship is 0.5
#' without inbreeding. Its diagonal is \eqn{1 + F_i}{1 + F[i]}, where
#' \eqn{F_i}{F[i]} is the inbreeding coefficient, so `diag(A) - 1` gives the
#' inbreeding coefficients. Founders (individuals with unknown parents) are
#' assumed to be unrelated and non-inbred, and an unknown parent is treated as
#' unrelated to everybody.
#'
#' [fastrel_ped()], [fastrel_vec()] and [rel_get()] can instead return the
#' kinship coefficient, Wright's coefficient of relationship or Hamilton's
#' regression coefficient of relatedness, with their argument `metric`. These
#' measures, what they mean and when to use them are described in
#' [relatedness_metrics].
#'
#' @section Range of values:
#' The additive relationship lies between 0 and 2, since it counts how many of
#' the two alleles of an individual are identical by descent with an allele of
#' another (equivalently, since a kinship coefficient is a probability and lies
#' between 0 and 1).
#' * The diagonal goes from 1 (not inbred) to 2 (fully inbred).
#' * The relationship between two individuals is at most the smaller of their
#'   two diagonal values, \eqn{A_{ij} \le \min(A_{ii}, A_{jj})}{A[i,j] <= min(A[i,i], A[j,j])}:
#'   nobody is more related to another individual than to itself.
#'   Between individuals that are not inbred it is therefore at most 1, the
#'   value for identical twins or clones; usual values are 0.5 for parent and
#'   offspring or for full siblings, and 0.25 for half siblings.
#' * Values above 1 thus only occur with inbreeding. For example, after `n`
#'   generations of selfing the inbreeding coefficient is `1 - (1/2)^n`, so the
#'   diagonal is 1.5, 1.75, 1.875, ... The maximum of 2 is approached but
#'   never reached, since founders are assumed not to be inbred.
#'
#' For a measure that cannot exceed 1 whatever the inbreeding, use
#' `metric = "kinship"` (a probability) or `metric = "wright"` (a correlation,
#' the same as `cov2cor(A)`): see [relatedness_metrics].
#'
#' @useDynLib fastrel, .registration = TRUE
#' @importFrom Rcpp evalCpp loadModule
#' @importFrom methods new
#' @importFrom utils head
"_PACKAGE"

## The class generator RelStore is created at load time by loadModule()
utils::globalVariables("RelStore")

Rcpp::loadModule("relstore", TRUE)
