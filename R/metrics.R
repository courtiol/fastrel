#' Measures of relatedness
#'
#' @description
#' [fastrel_ped()], [fastrel_vec()] and [rel_get()] can return four measures
#' of the genetic resemblance between individuals, chosen with their argument
#' `metric`. All four are derived from the same quantity, the kinship
#' coefficient, so choosing one or the other costs a single pass over the
#' matrix at most.
#'
#' | `metric` | Measure | Formula | Diagonal | Off-diagonal | Symmetric |
#' |---|---|---|---|---|---|
#' | `"additive"` (default) | additive (or numerator) relationship | \eqn{A_{ij} = 2\phi_{ij}}{A[i,j] = 2 phi[i,j]} | \eqn{1 + F_i}{1 + F[i]}, in \[1, 2\] | \[0, 2\] | yes |
#' | `"kinship"` | kinship (or coancestry) coefficient | \eqn{\phi_{ij}}{phi[i,j]} | \eqn{(1 + F_i)/2}{(1 + F[i]) / 2}, in \[0.5, 1\] | \[0, 1\] | yes |
#' | `"wright"` | Wright's coefficient of relationship | \eqn{A_{ij} / \sqrt{A_{ii} A_{jj}}}{A[i,j] / sqrt(A[i,i] A[j,j])} | 1 | \[0, 1\] | yes |
#' | `"hamilton"` | regression coefficient of relatedness, as in Hamilton's rule | \eqn{A_{ij} / A_{ii}}{A[i,j] / A[i,i]} | 1 | \[0, 1\] | no |
#'
#' Here \eqn{\phi_{ij}}{phi[i,j]} is the kinship coefficient between
#' individuals \eqn{i} and \eqn{j}, and \eqn{F_i}{F[i]} is the inbreeding
#' coefficient of \eqn{i}. Without inbreeding (all \eqn{F = 0}) the measures
#' `"additive"`, `"wright"` and `"hamilton"` are identical, and equal twice the
#' kinship coefficient: 0.5 for parent and offspring or full siblings, 0.25 for
#' half siblings, 0.125 for first cousins. They differ only with inbreeding.
#'
#' @section Kinship coefficient (`"kinship"`):
#' The kinship coefficient \eqn{\phi_{ij}}{phi[i,j]} (Malécot 1948, who called
#' it the coefficient of coancestry) is the probability that an allele drawn at
#' random from \eqn{i} and an allele drawn at random from \eqn{j}, at the same
#' autosomal locus, are identical by descent, i.e. copies of a single allele of
#' a common ancestor. Being a probability, it lies between 0 and 1. This is
#' what [kinship2::kinship()] returns (Sinnwell et al. 2014).
#'
#' Two related quantities follow. The inbreeding coefficient \eqn{F_i}{F[i]}
#' (Wright 1922) is the probability that the two alleles of \eqn{i} are
#' identical by descent; it equals the kinship coefficient between the parents
#' of \eqn{i}. The kinship of an individual with itself is
#' \eqn{\phi_{ii} = (1 + F_i)/2}{phi[i,i] = (1 + F[i]) / 2}: when two alleles
#' are drawn from \eqn{i} with replacement, the same allele is drawn twice with
#' probability 1/2, and otherwise the two alleles are identical by descent with
#' probability \eqn{F_i}{F[i]}.
#'
#' @section Additive relationship (`"additive"`):
#' The additive relationship \eqn{A_{ij} = 2\phi_{ij}}{A[i,j] = 2 phi[i,j]} is
#' what animal and plant breeders call the numerator relationship matrix
#' (Henderson 1976), computed here with the tabular method (Emik and Terrill
#' 1949). It is not a probability but an expected count: \eqn{A_{ij}}{A[i,j]}
#' is the expected number of the two alleles of \eqn{j} that are identical by
#' descent with an allele drawn at random from \eqn{i}. It therefore lies
#' between 0 and 2, and exceeds 1 only with inbreeding. For example, two
#' offspring obtained by selfing an individual that was itself obtained by
#' selfing have an additive relationship of 1.5 (and a kinship of 0.75).
#'
#' Its main use is quantitative genetics: under an additive model, the
#' covariance between the breeding values of \eqn{i} and \eqn{j} is
#' \eqn{A_{ij} \sigma^2_A}{A[i,j] * sigma2_A}, where \eqn{\sigma^2_A}{sigma2_A}
#' is the additive genetic variance of the base population (Lynch and Walsh
#' 1998, chapter 7). The diagonal, \eqn{1 + F_i}{1 + F[i]}, shows that the
#' breeding values of inbred individuals vary more, since inbred individuals
#' are more homozygous. This matrix is the one used in animal models and to
#' estimate heritability.
#'
#' @section Wright's coefficient of relationship (`"wright"`):
#' Wright (1922) defined the coefficient of relationship as the correlation
#' between the genetic values of two individuals, which is obtained by dividing
#' the covariance \eqn{A_{ij}}{A[i,j]} by the product of the two standard
#' deviations:
#' \deqn{r_{ij} = \frac{A_{ij}}{\sqrt{A_{ii} A_{jj}}} = \frac{2\phi_{ij}}{\sqrt{(1 + F_i)(1 + F_j)}}}{r[i,j] = A[i,j] / sqrt(A[i,i] * A[j,j]) = 2 phi[i,j] / sqrt((1 + F[i]) * (1 + F[j]))}
#' This is what [stats::cov2cor()] computes from the additive relationship
#' matrix. Being a correlation, it lies between 0 and 1 (the Cauchy-Schwarz
#' inequality guarantees that
#' \eqn{A_{ij} \le \sqrt{A_{ii} A_{jj}}}{A[i,j] <= sqrt(A[i,i] * A[j,j])}), and
#' every individual has a relationship of 1 with itself. It measures genetic
#' resemblance relative to the variability of each individual: a value of 1
#' means that two individuals are as alike as an individual is to itself, as
#' with members of a line selfed for many generations. It does not remove the
#' effect of inbreeding on resemblance, it puts it on a common scale. In the
#' selfing example above, the two offspring have \eqn{A = 1.5} but
#' \eqn{r = 1.5 / 1.75 \approx 0.86}{r = 1.5 / 1.75 = 0.86}.
#'
#' @section Hamilton's regression coefficient of relatedness (`"hamilton"`):
#' Hamilton (1964) used Wright's coefficient of relationship in the first
#' formulation of inclusive fitness. He later showed that with inbreeding the
#' relevant measure for kin selection, the \eqn{r} of Hamilton's rule
#' \eqn{rb > c}, is a regression coefficient rather than a correlation
#' (Hamilton 1972; Michod and Hamilton 1980): the regression of the genetic
#' value of the recipient of a behaviour on that of the actor,
#' \deqn{b_{ij} = \frac{\mathrm{Cov}(g_i, g_j)}{\mathrm{Var}(g_i)} = \frac{A_{ij}}{A_{ii}} = \frac{2\phi_{ij}}{1 + F_i}}{b[i,j] = Cov(g[i], g[j]) / Var(g[i]) = A[i,j] / A[i,i] = 2 phi[i,j] / (1 + F[i])}
#' where \eqn{i} is the actor and \eqn{j} the recipient. In the matrix
#' returned, **rows are actors and columns are recipients**: `R[i, j]` is the
#' relatedness of recipient `j` from the point of view of actor `i`. The
#' matrix is not symmetric when individuals differ in their inbreeding
#' coefficients. Its values lie between 0 and 1 (see the section Bounds), and
#' all individuals have a relatedness of 1 with themselves.
#'
#' @section Bounds:
#' In a pedigree, nobody is more related to another individual than to itself:
#' \deqn{A_{ij} \le \min(A_{ii}, A_{jj})}{A[i,j] <= min(A[i,i], A[j,j])}
#' This holds because, when \eqn{i} is added after its parents \eqn{m} and
#' \eqn{d}, \eqn{A_{ij} = (A_{mj} + A_{dj})/2}{A[i,j] = (A[m,j] + A[d,j]) / 2}
#' and \eqn{A_{ii} = 1 + A_{md}/2}{A[i,i] = 1 + A[m,d] / 2}. Each allele of
#' \eqn{j} can be identical by descent with an allele drawn from \eqn{m} and
#' with an allele drawn from \eqn{d} only if these two are themselves
#' identical by descent, so \eqn{A_{mj} + A_{dj} \le 2 + A_{md}}{A[m,j] + A[d,j] <= 2 + A[m,d]},
#' which gives \eqn{A_{ij} \le A_{ii}}{A[i,j] <= A[i,i]}; and
#' \eqn{A_{ij} \le A_{jj}}{A[i,j] <= A[j,j]} follows, by induction, from the
#' same property for the parents. The bounds of the table follow: the
#' off-diagonal values of the additive relationship are at most 2 and those of
#' the kinship at most 1, Wright's coefficient (bounded by the geometric mean
#' of the two diagonal values) and Hamilton's coefficient (bounded by the
#' diagonal value of the actor) are at most 1.
#'
#' @section Assumptions:
#' Whatever the measure, founders (individuals with unknown parents) are
#' assumed to be unrelated to each other and not inbred, an unknown parent is
#' treated as unrelated to everybody, and only autosomal loci are considered.
#' All values are therefore relative to the founders of the pedigree.
#'
#' @references
#' Emik, L. O., & Terrill, C. E. (1949). Systematic procedures for calculating
#' inbreeding coefficients. *Journal of Heredity*, 40(2), 51–55.
#' \doi{10.1093/oxfordjournals.jhered.a105986}
#'
#' Hamilton, W. D. (1964). The genetical evolution of social behaviour. I.
#' *Journal of Theoretical Biology*, 7(1), 1–16.
#' \doi{10.1016/0022-5193(64)90038-4}
#'
#' Hamilton, W. D. (1972). Altruism and related phenomena, mainly in social
#' insects. *Annual Review of Ecology and Systematics*, 3, 193–232.
#'
#' Henderson, C. R. (1976). A simple method for computing the inverse of a
#' numerator relationship matrix used in prediction of breeding values.
#' *Biometrics*, 32(1), 69–83. \doi{10.2307/2529339}
#'
#' Lynch, M., & Walsh, B. (1998). *Genetics and Analysis of Quantitative
#' Traits*. Sinauer Associates, Sunderland, MA. ISBN 978-0-87893-481-2.
#'
#' Malécot, G. (1948). *Les mathématiques de l'hérédité*. Masson et Cie,
#' Paris.
#'
#' Michod, R. E., & Hamilton, W. D. (1980). Coefficients of relatedness in
#' sociobiology. *Nature*, 288, 694–697. \doi{10.1038/288694a0}
#'
#' Sinnwell, J. P., Therneau, T. M., & Schaid, D. J. (2014). The kinship2 R
#' package for pedigree data. *Human Heredity*, 78(2), 91–93.
#' \doi{10.1159/000363105}
#'
#' Wright, S. (1922). Coefficients of inbreeding and relationship. *The
#' American Naturalist*, 56(645), 330–338. \doi{10.1086/279872}
#'
#' @seealso [fastrel_ped()], [fastrel_vec()] and [rel_get()], which take the
#'   argument `metric`.
#'
#' @examples
#' ## G is a founder; P is obtained by selfing G; O1 and O2 by selfing P
#' ped <- data.frame(id  = c("G", "P", "O1", "O2"),
#'                   mom = c(NA, "G", "P", "P"),
#'                   dad = c(NA, "G", "P", "P"))
#' fastrel_ped(ped)                       # additive relationship, the default
#' fastrel_ped(ped, metric = "kinship")   # probabilities
#' fastrel_ped(ped, metric = "wright")    # correlations
#' fastrel_ped(ped, metric = "hamilton")  # regressions: row = actor
#'
#' ## the conversions from the additive relationship matrix
#' A <- fastrel_ped(ped)
#' all.equal(fastrel_ped(ped, metric = "kinship"), A / 2)
#' all.equal(fastrel_ped(ped, metric = "wright"), cov2cor(A))
#' all.equal(fastrel_ped(ped, metric = "hamilton"), A / diag(A))
#'
#' ## without inbreeding, all but the kinship are the same
#' ped2 <- data.frame(id = c("a", "b", "c", "d"), mom = c(NA, NA, "a", "a"),
#'                    dad = c(NA, NA, "b", "b"))
#' all.equal(fastrel_ped(ped2), fastrel_ped(ped2, metric = "wright"))
#' all.equal(fastrel_ped(ped2), fastrel_ped(ped2, metric = "hamilton"))
#' @name relatedness_metrics
NULL

## names of the measures, in the order of their codes in src/metrics.h
metric_names <- c("additive", "kinship", "wright", "hamilton")

## integer code of a measure, for the C++ side
metric_code <- function(metric) {
  if (!is.character(metric) || length(metric) != 1L || is.na(metric))
    stop("`metric` must be one of: ", paste0('"', metric_names, '"', collapse = ", "))
  code <- pmatch(tolower(metric), metric_names)
  if (is.na(code))
    stop("unknown `metric` \"", metric, "\": must be one of ",
         paste0('"', metric_names, '"', collapse = ", "),
         ". See ?relatedness_metrics.")
  code - 1L
}
