#' Simulate a pedigree
#'
#' Simulates a pedigree with overlapping generations, some inbreeding, and some
#' unknown parents. Individuals are numbered in birth order, so that parents
#' always come before their offspring, and each individual picks its parents among
#' the `window` individuals born just before it. It is meant for testing and
#' benchmarking.
#'
#' Unless `seed = NULL`, the function calls [set.seed()], which changes the
#' state of the random number generator of the session.
#'
#' @param n Total number of individuals (larger than `n_founders`).
#' @param n_founders Number of individuals with no parents at the start.
#' @param window Number of preceding individuals among which parents are drawn.
#' @param p_unknown Probability that a parent (mother or father separately) is
#'   unknown for a non-founder.
#' @param seed Seed of the random number generator, or `NULL` to use the
#'   generator of the session as it is.
#'
#' @return A data frame with columns `id` (e.g. `"ind0001"`), `momid`, `dadid`
#'   and `sex` (1 = male, 2 = female).
#'
#' @examples
#' ped <- sim_ped(200)
#' head(ped)
#' A <- relatedness_fast(ped$id, ped$dadid, ped$momid)
#' range(diag(A) - 1)
#' @export
sim_ped <- function(n = 1000, n_founders = 40, window = 150,
                    p_unknown = 0.05, seed = 1) {
  if (n <= n_founders) stop("`n` must be larger than `n_founders`")
  if (n_founders < 1) stop("`n_founders` must be at least 1")
  if (window < 1) stop("`window` must be at least 1")
  if (p_unknown < 0 || p_unknown > 1) stop("`p_unknown` must be a probability")
  if (!is.null(seed)) set.seed(seed)
  sex <- sample(1:2, n, replace = TRUE)       # 1 = male, 2 = female
  mom <- dad <- rep(NA_integer_, n)

  for (i in (n_founders + 1):n) {
    pool <- max(1, i - window):(i - 1)        # recent individuals -> overlapping
    f <- pool[sex[pool] == 2]                 # generations and inbreeding
    m <- pool[sex[pool] == 1]
    if (length(f) && stats::runif(1) > p_unknown) mom[i] <- f[sample.int(length(f), 1)]
    if (length(m) && stats::runif(1) > p_unknown) dad[i] <- m[sample.int(length(m), 1)]
  }

  width <- max(4L, nchar(as.character(n)))
  ped <- data.frame(id = sprintf(paste0("ind%0", width, "d"), 1:n), sex = sex)
  ped$momid <- ped$id[mom]
  ped$dadid <- ped$id[dad]
  ped[, c("id", "momid", "dadid", "sex")]
}
