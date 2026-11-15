#' Simulate a pedigree
#'
#' Simulates a pedigree with overlapping generations, some inbreeding, and some
#' unknown parents. Individuals are numbered in birth order, so that parents
#' always come before their offspring, and each individual picks its parents among
#' the `window` individuals born just before it. It is meant for testing and
#' benchmarking.
#'
#' The result depends on the state of the random number generator of the
#' session: call [set.seed()] beforehand to make it reproducible.
#'
#' @param n Total number of individuals (larger than `n_founders`).
#' @param n_founders Number of individuals with no parents at the start.
#' @param window Number of preceding individuals among which parents are drawn.
#' @param p_unknown Probability that a parent (mother or father separately) is
#'   unknown for a non-founder.
#' @param numeric_ids If `TRUE`, individuals are identified by the integers 1
#'   to `n` (in birth order) and unknown parents by 0, as many programs of
#'   animal breeding expect (e.g. [ggroups::buildA()]). If `FALSE` (the
#'   default), identifiers are character strings and unknown parents are `NA`.
#'
#' @return A data frame with columns `id`, `momid`, `dadid` and `sex` (1 =
#'   male, 2 = female). Identifiers are character strings such as `"ind0001"`,
#'   with `NA` for unknown parents, or with `numeric_ids = TRUE` integers, with
#'   0 for unknown parents.
#'
#'   Both forms can be given as they are to [fastrel_ped()], [fastrel_vec()]
#'   and [rel_add()], which all take 0 as an unknown parent.
#'
#' @examples
#' set.seed(1)
#' ped <- sim_ped(200)
#' head(ped)
#' A <- fastrel_ped(ped)
#' range(diag(A) - 1)
#'
#' ## the same pedigree with integer identifiers and 0 for unknown parents
#' set.seed(1)
#' ped_num <- sim_ped(200, numeric_ids = TRUE)
#' head(ped_num)
#' all.equal(unname(fastrel_ped(ped_num)), unname(A))
#'
#' ## ready for other packages, e.g. ggroups (columns: individual, sire, dam)
#' if (requireNamespace("ggroups", quietly = TRUE)) {
#'   all.equal(unname(ggroups::buildA(ped_num[, c("id", "dadid", "momid")])),
#'             unname(A))
#' }
#' @export
sim_ped <- function(n = 1000, n_founders = 40, window = 150,
                    p_unknown = 0.05, numeric_ids = FALSE) {
  if (!is.logical(numeric_ids) || length(numeric_ids) != 1L || is.na(numeric_ids))
    stop("`numeric_ids` must be TRUE or FALSE")
  if (n <= n_founders) stop("`n` must be larger than `n_founders`")
  if (n_founders < 1) stop("`n_founders` must be at least 1")
  if (window < 1) stop("`window` must be at least 1")
  if (p_unknown < 0 || p_unknown > 1) stop("`p_unknown` must be a probability")
  sex <- sample(1:2, n, replace = TRUE)       # 1 = male, 2 = female
  mom <- dad <- rep(NA_integer_, n)

  for (i in (n_founders + 1):n) {
    pool <- max(1, i - window):(i - 1)        # recent individuals -> overlapping
    f <- pool[sex[pool] == 2]                 # generations and inbreeding
    m <- pool[sex[pool] == 1]
    if (length(f) && stats::runif(1) > p_unknown) mom[i] <- f[sample.int(length(f), 1)]
    if (length(m) && stats::runif(1) > p_unknown) dad[i] <- m[sample.int(length(m), 1)]
  }

  if (numeric_ids) {
    mom[is.na(mom)] <- 0L
    dad[is.na(dad)] <- 0L
    return(data.frame(id = seq_len(n), momid = mom, dadid = dad, sex = sex))
  }

  width <- max(4L, nchar(as.character(n)))
  ped <- data.frame(id = sprintf(paste0("ind%0", width, "d"), 1:n), sex = sex)
  ped$momid <- ped$id[mom]
  ped$dadid <- ped$id[dad]
  ped[, c("id", "momid", "dadid", "sex")]
}
