#' Simulate a pedigree with deaths, emigration and immigration
#'
#' Simulates an individual-based population over discrete time steps and returns
#' the complete pedigree together with the step at which each individual
#' arrived and left. It is meant for testing and benchmarking [rel_add()] and
#' [rel_remove()], since the number of individuals alive stays roughly constant
#' while the pedigree keeps growing.
#'
#' Step 0 creates `n_init` adult founders with unknown parents. At each
#' following step `t`:
#' 1. Each mature female alive at the start of the step reproduces with
#'    probability `p_birth * (1 - N / carrying)` (bounded to `[0, 1]`), where `N`
#'    is the number of individuals alive. The father is drawn at random among
#'    the mature males (`NA` if there is none). Newborns arrive at step `t`.
#' 2. A Poisson number (mean `n_immigrants`) of adult immigrants with unknown
#'    parents arrive at step `t`.
#' 3. Each individual that was alive at the start of the step dies with
#'    probability `p_death`, and otherwise, if mature, emigrates with probability
#'    `p_emigrate`. Its `exit` is `t`.
#'
#' Newborns and immigrants cannot leave in the step in which they arrive, and
#' parents are always alive at the start of the step in which their offspring are
#' born. Therefore, replaying the steps in order with `rel_add()` for the
#' arrivals followed by `rel_remove()` for the exits never refers to an unknown
#' individual or to a removed parent (see the example).
#'
#' Unless `seed = NULL`, the function calls [set.seed()], which changes the
#' state of the random number generator of the session.
#'
#' @param n_steps Number of time steps after step 0.
#' @param n_init Number of founders at step 0.
#' @param carrying Carrying capacity: the birth probability decreases linearly to
#'   0 when `N` reaches it.
#' @param maturity Number of steps after birth before an individual can
#'   reproduce or emigrate. Founders and immigrants are mature on arrival.
#' @param p_birth Maximum probability that a mature female reproduces in a step
#'   (one offspring at most per step).
#' @param p_death Probability of death per step.
#' @param p_emigrate Probability of emigration per step for mature individuals
#'   that did not die.
#' @param n_immigrants Mean number of immigrants per step (Poisson).
#' @param seed Seed of the random number generator, or `NULL` to use the
#'   generator of the session as it is.
#'
#' @return A data frame with one row per individual ever present, in order of
#'   arrival, so parents always come before their offspring, and with columns:
#'   * `id`, `momid`, `dadid`: identifiers (`NA` for unknown parents);
#'   * `sex`: 1 = male, 2 = female;
#'   * `arrival`: step at which the individual enters the pedigree (0 for
#'     founders);
#'   * `exit`: step at which it leaves, `NA` if it is still present at the end;
#'   * `fate`: `"dead"`, `"emigrated"` or `NA`;
#'   * `origin`: `"founder"`, `"born"` or `"immigrant"`.
#'
#' @examples
#' pop <- sim_pop(n_steps = 40, n_init = 60, seed = 2)
#' table(pop$origin)
#' table(pop$fate, useNA = "ifany")
#'
#' ## replay the simulation with the incremental store
#' rel <- rel_new(300)
#' size <- integer(41)
#' for (t in 0:40) {
#'   rel_add(rel, pop[pop$arrival == t, c("id", "momid", "dadid")])
#'   rel_remove(rel, pop$id[!is.na(pop$exit) & pop$exit == t])
#'   size[t + 1] <- rel_size(rel)
#' }
#' size
#'
#' ## the survivors agree with the relatedness computed from the whole pedigree
#' A <- relatedness_fast(pop$id, pop$dadid, pop$momid)
#' alive <- rel_ids(rel)
#' all.equal(rel_get(rel), A[alive, alive])
#' @seealso [sim_ped()], [rel_remove()]
#' @export
sim_pop <- function(n_steps = 50, n_init = 100, carrying = 2 * n_init,
                    maturity = 2, p_birth = 0.6, p_death = 0.1,
                    p_emigrate = 0.02, n_immigrants = 1, seed = 1) {
  whole <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && x == round(x)
  if (!whole(n_steps) || n_steps < 0) stop("`n_steps` must be a non-negative whole number")
  if (!whole(n_init) || n_init < 1) stop("`n_init` must be a whole number, at least 1")
  if (!whole(maturity) || maturity < 0) stop("`maturity` must be a non-negative whole number")
  if (carrying <= 0) stop("`carrying` must be positive")
  if (any(c(p_birth, p_death, p_emigrate) < 0 | c(p_birth, p_death, p_emigrate) > 1))
    stop("`p_birth`, `p_death` and `p_emigrate` must be probabilities")
  if (n_immigrants < 0) stop("`n_immigrants` must be non-negative")
  if (!is.null(seed)) set.seed(seed)

  ## Individuals are numbered 1, 2, ... in order of arrival; the numbers are
  ## turned into identifiers at the end. Records are kept as plain vectors, one
  ## chunk per step, and assembled into a data frame only once.
  n_total <- n_init                               # individuals ever present

  ## state of the individuals currently present
  a_id  <- seq_len(n_init)
  a_sex <- sample(1:2, n_init, replace = TRUE)
  a_mat <- rep(0L, n_init)                        # step from which it is mature

  ## records of the arrivals (founders first) and of the exits
  r_mom <- r_dad <- r_sex <- r_arr <- r_ori <- vector("list", n_steps + 1L)
  r_mom[[1L]] <- r_dad[[1L]] <- rep(NA_integer_, n_init)
  r_sex[[1L]] <- a_sex
  r_arr[[1L]] <- rep(0L, n_init)
  r_ori[[1L]] <- rep("founder", n_init)
  x_id <- x_step <- x_fate <- vector("list", n_steps)

  for (t in seq_len(n_steps)) {
    n_alive <- length(a_id)
    if (n_alive == 0L) {
      warning("the population went extinct at step ", t - 1L)
      break
    }
    mature <- a_mat <= t

    ## 1. births
    mothers <- a_id[mature & a_sex == 2L]
    fathers <- a_id[mature & a_sex == 1L]
    p_b <- min(1, max(0, p_birth * (1 - n_alive / carrying)))
    mothers <- mothers[stats::runif(length(mothers)) < p_b]
    nb <- length(mothers)
    dads <- sex_b <- integer(0)
    if (nb > 0L) {
      dads <- if (length(fathers)) fathers[sample.int(length(fathers), nb, replace = TRUE)]
              else rep(NA_integer_, nb)
      sex_b <- sample(1:2, nb, replace = TRUE)
    }

    ## 2. immigrants
    ni <- stats::rpois(1, n_immigrants)
    sex_i <- if (ni > 0L) sample(1:2, ni, replace = TRUE) else integer(0)

    ## 3. exits among those present at the start of the step
    dead <- stats::runif(n_alive) < p_death
    emig <- !dead & mature & stats::runif(n_alive) < p_emigrate
    gone <- dead | emig
    if (any(gone)) {
      x_id[[t]]   <- a_id[gone]
      x_step[[t]] <- rep(t, sum(gone))
      x_fate[[t]] <- ifelse(dead[gone], "dead", "emigrated")
    }

    ## record the arrivals and update the individuals present
    na <- nb + ni
    if (na > 0L) {
      r_mom[[t + 1L]] <- c(mothers, rep(NA_integer_, ni))
      r_dad[[t + 1L]] <- c(dads, rep(NA_integer_, ni))
      r_sex[[t + 1L]] <- c(sex_b, sex_i)
      r_arr[[t + 1L]] <- rep(t, na)
      r_ori[[t + 1L]] <- rep(c("born", "immigrant"), c(nb, ni))
    }
    a_id  <- c(a_id[!gone],  n_total + seq_len(na))
    a_sex <- c(a_sex[!gone], sex_b, sex_i)
    a_mat <- c(a_mat[!gone],
               c(rep(t + as.integer(maturity), nb), rep(t, ni)))
    n_total <- n_total + na
  }

  ids <- sprintf("ind%08d", seq_len(n_total))
  exit <- rep(NA_integer_, n_total)
  fate <- rep(NA_character_, n_total)
  gone <- unlist(x_id)
  exit[gone] <- unlist(x_step)
  fate[gone] <- unlist(x_fate)
  data.frame(id = ids, momid = ids[unlist(r_mom)], dadid = ids[unlist(r_dad)],
             sex = unlist(r_sex), arrival = unlist(r_arr), exit = exit,
             fate = fate, origin = unlist(r_ori), stringsAsFactors = FALSE)
}
