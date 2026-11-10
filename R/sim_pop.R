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
#' The function calls [set.seed()], which changes the state of the random number
#' generator of the session.
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
#' @param seed Seed of the random number generator.
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
#' alive <- rel$ids
#' all.equal(rel_get(rel), A[alive, alive])
#' @seealso [sim_ped()], [rel_remove()]
#' @export
sim_pop <- function(n_steps = 50, n_init = 100, carrying = 2 * n_init,
                    maturity = 2, p_birth = 0.6, p_death = 0.1,
                    p_emigrate = 0.02, n_immigrants = 1, seed = 1) {
  if (n_init < 1) stop("`n_init` must be at least 1")
  if (carrying <= 0) stop("`carrying` must be positive")
  if (any(c(p_birth, p_death, p_emigrate) < 0 | c(p_birth, p_death, p_emigrate) > 1))
    stop("`p_birth`, `p_death` and `p_emigrate` must be probabilities")
  if (n_immigrants < 0) stop("`n_immigrants` must be non-negative")
  set.seed(seed)

  next_id <- 0L
  new_ids <- function(k) {
    ids <- sprintf("ind%08d", next_id + seq_len(k))
    next_id <<- next_id + k
    ids
  }

  ## state of the individuals currently present
  a_id <- new_ids(n_init)
  a_sex <- sample(1:2, n_init, replace = TRUE)
  a_mat <- rep(0L, n_init)                      # step from which it is mature

  ## records
  chunks <- list(data.frame(id = a_id, momid = NA_character_, dadid = NA_character_,
                            sex = a_sex, arrival = 0L, origin = "founder",
                            stringsAsFactors = FALSE))
  exits <- vector("list", n_steps)

  for (t in seq_len(n_steps)) {
    n_alive <- length(a_id)
    if (n_alive == 0L) {
      warning("the population went extinct at step ", t - 1L)
      exits <- exits[seq_len(t - 1L)]
      break
    }
    mature <- a_mat <= t

    ## 1. births
    mothers <- a_id[mature & a_sex == 2L]
    fathers <- a_id[mature & a_sex == 1L]
    p_b <- min(1, max(0, p_birth * (1 - n_alive / carrying)))
    mothers <- mothers[stats::runif(length(mothers)) < p_b]
    nb <- length(mothers)
    born <- NULL
    if (nb > 0L) {
      dads <- if (length(fathers)) fathers[sample.int(length(fathers), nb, replace = TRUE)]
              else rep(NA_character_, nb)
      born <- data.frame(id = new_ids(nb), momid = mothers, dadid = dads,
                         sex = sample(1:2, nb, replace = TRUE), arrival = t,
                         origin = "born", stringsAsFactors = FALSE)
    }

    ## 2. immigrants
    ni <- stats::rpois(1, n_immigrants)
    imm <- NULL
    if (ni > 0L) {
      imm <- data.frame(id = new_ids(ni), momid = NA_character_, dadid = NA_character_,
                        sex = sample(1:2, ni, replace = TRUE), arrival = t,
                        origin = "immigrant", stringsAsFactors = FALSE)
    }

    ## 3. exits among those present at the start of the step
    dead <- stats::runif(n_alive) < p_death
    emig <- !dead & mature & stats::runif(n_alive) < p_emigrate
    gone <- dead | emig
    if (any(gone))
      exits[[t]] <- data.frame(id = a_id[gone], exit = t,
                               fate = ifelse(dead[gone], "dead", "emigrated"),
                               stringsAsFactors = FALSE)

    ## update the individuals present
    arrivals <- rbind(born, imm)
    a_id  <- c(a_id[!gone],  arrivals$id)
    a_sex <- c(a_sex[!gone], arrivals$sex)
    a_mat <- c(a_mat[!gone],
               c(rep(t + as.integer(maturity), nb), rep(t, ni)))
    if (!is.null(arrivals)) chunks[[length(chunks) + 1L]] <- arrivals
  }

  ped <- do.call(rbind, chunks)
  ex <- do.call(rbind, exits)
  k <- if (is.null(ex)) rep(NA_integer_, nrow(ped)) else match(ped$id, ex$id)
  ped$exit <- if (is.null(ex)) NA_integer_ else as.integer(ex$exit[k])
  ped$fate <- if (is.null(ex)) NA_character_ else ex$fate[k]
  rownames(ped) <- NULL
  ped[, c("id", "momid", "dadid", "sex", "arrival", "exit", "fate", "origin")]
}
