## Validation against nadiv::makeA(), an independent implementation of the
## additive (numerator) relationship matrix used in evolutionary ecology, on
## the same pedigrees as for kinship2 and ggroups (see helper-pedigrees.R).

skip_if_not_installed("nadiv")

## makeA() takes the columns ID, Dam, Sire, and returns a sparse matrix whose
## rows may be reordered; it warns about selfing, which it handles correctly
nadiv_additive <- function(p) {
  A <- suppressWarnings(nadiv::makeA(data.frame(ID = p$id, Dam = p$mom, Sire = p$dad)))
  A <- as.matrix(A)
  id <- as.character(p$id)
  A[id, id, drop = FALSE]
}

## makeA() fails when no parent is known at all (a limitation of nadiv), so
## those pedigrees are checked against kinship2 and ggroups only
nadiv_pedigrees <- function() {
  peds <- validation_pedigrees()
  peds[!names(peds) %in% c("founders", "single")]
}

test_that("all measures match nadiv on pedigrees of various shapes", {
  peds <- nadiv_pedigrees()
  set.seed(3)
  for (nm in names(peds))
    expect_matches_reference(peds[[nm]], nadiv_additive(peds[[nm]]), what = nm)
})

test_that("nadiv and kinship2 agree with each other on these pedigrees", {
  skip_if_not_installed("kinship2")
  peds <- nadiv_pedigrees()
  for (nm in names(peds)) {
    p <- peds[[nm]]
    expect_equal(unname(nadiv_additive(p)),
                 unname(2 * kinship2::kinship(p$id, dadid = p$dad, momid = p$mom)),
                 tolerance = 1e-12, label = nm)
  }
})
