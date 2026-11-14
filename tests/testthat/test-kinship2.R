## Validation against kinship2::kinship(), an independent implementation of
## the kinship coefficient (Sinnwell, Therneau & Schaid 2014), on pedigrees of
## various shapes (see helper-pedigrees.R). All four measures of fastrel are
## checked, each being a documented function of the kinship coefficient.

skip_if_not_installed("kinship2")

kinship2_additive <- function(p) 2 * kinship2::kinship(p$id, dadid = p$dad, momid = p$mom)

test_that("all measures match kinship2 on pedigrees of various shapes", {
  peds <- validation_pedigrees()
  set.seed(1)
  for (nm in names(peds))
    expect_matches_reference(peds[[nm]], kinship2_additive(peds[[nm]]), what = nm)
})

test_that("kinship2 and fastrel agree on special inputs", {
  ## a parent absent from id is an unknown parent for both
  p <- data.frame(id = c("a1", "a2", "x"), dad = c(NA, NA, "ghost"), mom = c(NA, NA, "a1"))
  expect_equal(fastrel_ped(p, metric = "kinship"),
               kinship2::kinship(p$id, dadid = p$dad, momid = p$mom), tolerance = 1e-12)
  ## factor identifiers, and rows not in birth order (kinship2 sorts them too)
  f <- data.frame(id = factor(c("r", "p", "q")), dad = factor(c("p", NA, NA)),
                  mom = factor(c("q", NA, NA)))
  expect_equal(fastrel_ped(f, metric = "kinship"),
               kinship2::kinship(f$id, dadid = f$dad, momid = f$mom), tolerance = 1e-12)
  ## kinship2's own example pedigree
  ex <- get(utils::data("sample.ped", package = "kinship2", envir = environment()))
  expect_equal(fastrel_vec(ex$id, ex$father, ex$mother, metric = "kinship"),
               kinship2::kinship(ex$id, dadid = ex$father, momid = ex$mother),
               tolerance = 1e-12)
})

test_that("a store after removals keeps the values kinship2 gives on the whole pedigree", {
  set.seed(11)
  pop <- sim_pop(n_steps = 60, n_init = 40, K = 80)
  A_ref <- kinship2_additive(data.frame(id = pop$id, dad = pop$dadid, mom = pop$momid))
  ref <- metrics_from_additive(A_ref)
  rel <- rel_new(remember_removed = TRUE)
  for (t in 0:60) {
    rel_add(rel, pop[pop$arrival == t, ])
    rel_remove(rel, pop$id[!is.na(pop$exit) & pop$exit == t])
  }
  alive <- rel_ids(rel)
  expect_gt(length(alive), 10)
  expect_lt(length(alive), nrow(pop))
  for (m in names(ref))
    expect_equal(rel_get(rel, alive, metric = m), ref[[m]][alive, alive], tolerance = 1e-12)
})
