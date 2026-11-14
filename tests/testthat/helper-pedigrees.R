## Pedigrees of various shapes used to validate fastrel against independent
## implementations (test-kinship2.R, test-ggroups.R). Each is a data frame with
## columns id, dad, mom, in birth order (parents before offspring), with NA for
## unknown parents.
validation_pedigrees <- function() {
  peds <- list()

  peds$founders <- data.frame(id = c("a", "b", "c"), dad = NA, mom = NA)
  peds$single <- data.frame(id = "a", dad = NA, mom = NA)

  ## full siblings, half siblings, a single known parent
  peds$nuclear <- data.frame(
    id  = c("F", "M1", "M2", "s1", "s2", "h1", "u1", "u2"),
    dad = c(NA,  NA,   NA,   "F",  "F",  "F",  "F",  NA),
    mom = c(NA,  NA,   NA,   "M1", "M1", "M2", NA,   "M2"))

  ## grandparents G1 x G2 have sons S1, S2 and daughter D; S1 x W1 and S2 x W2
  ## have cousins; two brothers marrying two sisters give double first
  ## cousins; g is the offspring of first cousins
  peds$cousins <- data.frame(
    id  = c("G1", "G2", "S1", "S2", "D", "W1", "W2", "c1", "c2", "c3",
            "B1", "B2", "Z1", "Z2", "d1", "d2", "g"),
    dad = c(NA, NA, "G1", "G1", "G1", NA, NA, "S1", "S2", NA,
            NA, NA, NA, NA, "B1", "B2", "c1"),
    mom = c(NA, NA, "G2", "G2", "G2", NA, NA, "W1", "W2", "D",
            NA, NA, NA, NA, "Z1", "Z2", "c2"))

  ## mother-son and father-daughter matings
  peds$parent_offspring <- data.frame(
    id  = c("A", "B", "C", "D", "E", "F"),
    dad = c(NA, NA, "A", "C", "A", "A"),
    mom = c(NA, NA, "B", "B", "B", "E"))

  ## ten generations of full-sibling mating (F reaches about 0.89)
  g <- 10
  peds$sib_mating <- data.frame(
    id  = c("m0", "f0", as.vector(rbind(paste0("m", 1:g), paste0("f", 1:g)))),
    dad = c(NA, NA, rep(paste0("m", 0:(g - 1)), each = 2)),
    mom = c(NA, NA, rep(paste0("f", 0:(g - 1)), each = 2)))

  ## fifteen generations of selfing (an individual is both parents)
  id <- paste0("s", 0:15)
  peds$selfing <- data.frame(id = id, dad = c(NA, id[-16]), mom = c(NA, id[-16]))

  ## unconnected families
  peds$families <- data.frame(
    id  = c("a1", "a2", "a3", "b1", "b2", "b3", "b4"),
    dad = c(NA, NA, "a1", NA, NA, "b1", "b1"),
    mom = c(NA, NA, "a2", NA, NA, "b2", "b2"))

  ## numeric identifiers
  peds$numeric <- data.frame(id = 1:6, dad = c(NA, NA, 1, 1, 3, 3),
                             mom = c(NA, NA, 2, 2, 4, 4))

  ## simulated pedigrees
  settings <- list(
    default_like  = list(n = 500,  n_founders = 40, window = 150, p_unknown = 0.05),
    heavy_inbreed = list(n = 500,  n_founders = 3,  window = 6,   p_unknown = 0),
    many_unknown  = list(n = 500,  n_founders = 20, window = 30,  p_unknown = 0.5),
    tiny_pop      = list(n = 800,  n_founders = 2,  window = 3,   p_unknown = 0.1),
    larger        = list(n = 1500, n_founders = 50, window = 300, p_unknown = 0.02))
  for (s in names(settings)) {
    set.seed(match(s, names(settings)))
    p <- do.call(sim_ped, settings[[s]])
    peds[[paste0("sim_ped_", s)]] <- data.frame(id = p$id, dad = p$dadid, mom = p$momid)
  }
  set.seed(11)
  pop <- sim_pop(n_steps = 60, n_init = 40, K = 80)
  peds$sim_pop <- data.frame(id = pop$id, dad = pop$dadid, mom = pop$momid)

  peds
}

## fastrel's four measures derived from an additive relationship matrix A, as
## documented in ?relatedness_metrics
metrics_from_additive <- function(A) {
  list(additive = A, kinship = A / 2, wright = cov2cor(A), hamilton = A / diag(A))
}

## compare all four measures of fastrel_vec(), fastrel_ped() and of a store
## (fed with the rows shuffled) with those derived from a reference additive
## relationship matrix `A_ref`, computed on the pedigree `p`
expect_matches_reference <- function(p, A_ref, what, tolerance = 1e-12) {
  ref <- metrics_from_additive(A_ref)
  ped <- data.frame(id = p$id, mom = p$mom, dad = p$dad)
  shuffled <- ped[sample(nrow(ped)), ]
  rel <- rel_new()
  rel_add(rel, shuffled)
  for (m in names(ref)) {
    lab <- function(f) paste0(what, ": ", f, "(metric = \"", m, "\")")
    expect_equal(unname(fastrel_vec(p$id, p$dad, p$mom, metric = m)), unname(ref[[m]]),
                 tolerance = tolerance, label = lab("fastrel_vec"))
    expect_equal(unname(fastrel_ped(shuffled, metric = m)[as.character(p$id), as.character(p$id), drop = FALSE]),
                 unname(ref[[m]]), tolerance = tolerance, label = lab("fastrel_ped, shuffled rows,"))
    expect_equal(unname(rel_get(rel, p$id, metric = m)), unname(ref[[m]]),
                 tolerance = tolerance, label = lab("rel_get"))
  }
}
