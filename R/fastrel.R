#' Relatedness matrix from a pedigree
#'
#' Compute a measure of relatedness of autosomal loci between all individuals
#' of a pedigree, using the tabular method (Emik and Terrill 1949; Henderson
#' 1976) implemented in C++. By default, the
#' measure is the additive relationship \eqn{A = 2\phi}{A = 2 * phi}, where
#' \eqn{\phi}{phi} is the kinship coefficient; the argument `metric` gives the
#' kinship coefficient itself, Wright's coefficient of relationship or
#' Hamilton's regression coefficient of relatedness instead (see
#' [relatedness_metrics]). The two functions do the same thing and differ only
#' in how the pedigree is given:
#' * `fastrel_ped()` takes a data frame;
#' * `fastrel_vec()` takes three vectors, as [kinship2::kinship()] does. It is a
#'   fast replacement for `2 * kinship2::kinship()`, or for
#'   `kinship2::kinship()` when `metric = "kinship"`.
#'
#' Rows of the pedigree can be in any order: parents do not need to appear
#' before their offspring. Parents that are `NA` or 0 (as in many programs of
#' animal breeding), or that are not found among the individuals, are treated
#' as unknown; 0 is refused if an individual is itself called 0, as it would
#' then be ambiguous. The algorithm uses a topological
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
#' taken to be the focal individual, its father and its mother, in this order,
#' as for the arguments of `fastrel_vec()` and [kinship2::kinship()] and as in
#' pedigrees used in animal breeding (animal, sire, dam); a message says so the
#' first time this happens in an R session.
#' If only some of the three are recognised, or if several columns could play
#' the same role, an error is raised rather than guessing.
#'
#' @param ped A data frame (or a tibble) with one row per individual and three
#'   columns giving the identifier of the individual, of its father and of its
#'   mother (`NA` or 0 if unknown). See the section Columns of `ped`.
#' @param id Vector of unique individual identifiers. For this argument and
#'   the next two, a data frame with a single column (what `x[, "id"]` returns
#'   when `x` is a tibble) also works.
#' @param dadid Vector of the same length as `id` with the identifier of each
#'   individual's father (`NA` or 0 if unknown).
#' @param momid Vector of the same length as `id` with the identifier of each
#'   individual's mother (`NA` or 0 if unknown).
#' @param sex,chrtype Arguments of [kinship2::kinship()], present so that
#'   calls written for that function work unchanged. `sex` is not used, as the
#'   sex of individuals does not matter for autosomes; if given, it must have
#'   one value per individual. `chrtype` must be `"autosome"`, the only type of
#'   chromosome supported: anything else raises an error.
#' @param metric The measure of relatedness to return: `"additive"` (the
#'   default), `"kinship"`, `"wright"` (or `"Wright"`) or `"hamilton"` (or
#'   `"Hamilton"`). Abbreviations are accepted (e.g. `"kin"`). See
#'   [relatedness_metrics] for their definitions. In `fastrel_vec()` it comes
#'   after the arguments used for compatibility with [kinship2::kinship()], so
#'   it is best given by name.
#'
#' @return A numeric matrix with one row and one column per individual, in the
#'   order given, and the identifiers as row and column names. With the default
#'   `metric = "additive"`, the matrix is symmetric, its diagonal is `1 + F`,
#'   where `F` is the inbreeding coefficient, and its values lie between 0 and
#'   2, exceeding 1 only with inbreeding. For the other measures, see
#'   [relatedness_metrics]; note that with `metric = "hamilton"` the matrix is
#'   not symmetric (rows are actors, columns are recipients).
#'
#' @references
#' Emik, L. O., & Terrill, C. E. (1949). Systematic procedures for calculating
#' inbreeding coefficients. *Journal of Heredity*, 40(2), 51–55.
#' \doi{10.1093/oxfordjournals.jhered.a105986}
#'
#' Henderson, C. R. (1976). A simple method for computing the inverse of a
#' numerator relationship matrix used in prediction of breeding values.
#' *Biometrics*, 32(1), 69–83. \doi{10.2307/2529339}
#'
#' @seealso [relatedness_metrics] for the measures of relatedness;
#'   [rel_new()] for the incremental version.
#'
#' @examples
#' ped <- data.frame(
#'   id  = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
#'   dad = c(NA, NA, NA, NA, "B", "D", "F", "G", "G"),
#'   mom = c(NA, NA, NA, NA, "A", "C", "E", "E", "E"))
#' A <- fastrel_ped(ped)
#' A
#'
#' ## inbreeding coefficients (H and I are the offspring of a mother-son mating)
#' diag(A) - 1
#'
#' ## other measures of relatedness (see ?relatedness_metrics)
#' fastrel_ped(ped, metric = "kinship")
#' fastrel_ped(ped, metric = "wright")
#' fastrel_ped(ped, metric = "hamilton")
#'
#' ## the same from vectors
#' all.equal(A, fastrel_vec(ped$id, dadid = ped$dad, momid = ped$mom))
#'
#' ## columns are found by their names, whatever their order
#' ped2 <- data.frame(Sire = ped$dad, Dam = ped$mom, Animal = ped$id)
#' all.equal(A, fastrel_ped(ped2))
#'
#' ## the same as twice the kinship matrix of the kinship2 package
#' if (requireNamespace("kinship2", quietly = TRUE)) {
#'   all.equal(A, 2 * kinship2::kinship(ped$id, dadid = ped$dad, momid = ped$mom))
#'   all.equal(fastrel_ped(ped, metric = "kinship"),
#'             kinship2::kinship(ped$id, dadid = ped$dad, momid = ped$mom))
#' }
#' @export
fastrel_ped <- function(ped,
                        metric = c("additive", "kinship", "wright", "hamilton", "Wright", "Hamilton")) {
  metric <- match.arg(metric)
  cols <- ped_columns(ped)
  fastrel_vec(ped[[cols[1L]]], dadid = ped[[cols[3L]]], momid = ped[[cols[2L]]],
              metric = metric)
}

#' @rdname fastrel_ped
#' @export
fastrel_vec <- function(id, dadid, momid, sex, chrtype = "autosome",
                        metric = c("additive", "kinship", "wright", "hamilton", "Wright", "Hamilton")) {
  metric <- match.arg(metric)
  code <- metric_code(metric)
  if (!is.character(chrtype) || length(chrtype) != 1L ||
      !identical(tolower(chrtype), "autosome"))
    stop("only autosomes are supported (`chrtype = \"autosome\"`)")
  id    <- id_vector(id, "id")
  dadid <- id_vector(dadid, "dadid")
  momid <- id_vector(momid, "momid")
  if (anyNA(id)) stop("id values cannot be NA")
  if (anyDuplicated(id)) stop("All id values must be unique")
  n <- length(id)
  if (length(dadid) != n) stop("Invalid father id")
  if (length(momid) != n) stop("Invalid mother id")
  momid <- zero_as_unknown(momid, id)             # 0 = unknown parent
  dadid <- zero_as_unknown(dadid, id)
  if (!missing(sex) && length(sex) != n)
    stop("`sex` must have one value per individual (it is not used, since ",
         "only autosomes are supported). To choose the measure returned, use ",
         "`metric = `.")
  if (n == 0L) return(matrix(numeric(0), 0, 0))

  A <- relatedness_cpp(mom = match(momid, id), dad = match(dadid, id), metric = code)
  dimnames(A) <- list(id, id)
  A
}
