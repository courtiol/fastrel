#' fastrel: fast relatedness, relationship and kinship matrices from pedigrees
#'
#' Tools to compute, for autosomal loci of diploid individuals, the additive
#' relationship ("relatedness") matrix \eqn{A = 2\phi}{A = 2 * phi}, where
#' \eqn{\phi}{phi} is the kinship coefficient (Malécot 1948), as well as the
#' kinship coefficient itself, Wright's coefficient of relationship (Wright
#' 1922) and Hamilton's regression coefficient of relatedness (Hamilton 1972).
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
#' between 0 and 2, and can exceed 1 (see [relatedness_metrics]). For
#' example, an allele drawn from a parent is one of the two alleles of its
#' offspring with probability 1/2, so their additive relationship is 0.5
#' without inbreeding. Its diagonal is \eqn{1 + F_i}{1 + F[i]}, where
#' \eqn{F_i}{F[i]} is the inbreeding coefficient, so `diag(A) - 1` gives the
#' inbreeding coefficients. Founders (individuals with unknown parents) are
#' assumed to be unrelated and non-inbred, and an unknown parent is treated as
#' unrelated to everybody.
#'
#' [fastrel_ped()], [fastrel_vec()] and [rel_get()] can instead return the
#' kinship coefficient (Malécot 1948), Wright's coefficient of relationship
#' (Wright 1922) or Hamilton's regression coefficient of relatedness
#' (Hamilton 1972), with their argument `metric`. These measures, what they
#' mean and when to use them are described in [relatedness_metrics].
#'
#' The matrices are computed with the tabular method (Emik and Terrill 1949),
#' implemented in C++. See [relatedness_metrics] for the full list of
#' references.
#'
#' @references
#' Emik, L. O., & Terrill, C. E. (1949). Systematic procedures for calculating
#' inbreeding coefficients. *Journal of Heredity*, 40(2), 51–55.
#' \doi{10.1093/oxfordjournals.jhered.a105986}
#'
#' Hamilton, W. D. (1972). Altruism and related phenomena, mainly in social
#' insects. *Annual Review of Ecology and Systematics*, 3, 193–232.
#' \doi{10.1146/annurev.es.03.110172.001205}
#'
#' Malécot, G. (1948). *Les mathématiques de l'hérédité*. Masson et Cie,
#' Paris.
#'
#' Wright, S. (1922). Coefficients of inbreeding and relationship. *The
#' American Naturalist*, 56(645), 330–338. \doi{10.1086/279872}
#'
#' @useDynLib fastrel, .registration = TRUE
#' @importFrom Rcpp evalCpp loadModule
#' @importFrom methods new
#' @importFrom utils head
"_PACKAGE"

## The class generator RelStore is created at load time by loadModule()
utils::globalVariables("RelStore")

Rcpp::loadModule("relstore", TRUE)
