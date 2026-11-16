## Validation against nadiv::makeA(), an independent implementation of the
## additive (numerator) relationship matrix used in evolutionary ecology, on
## the same pedigrees as for kinship2 and ggroups (see helper-pedigrees.R).
## Only the additive relationship is compared: the other measures are derived
## from it in fastrel, and are checked against kinship2 (test-kinship2.R).

skip_if_not_installed("nadiv")

## makeA() takes the columns ID, Dam, Sire, and returns a sparse matrix whose
## rows may be reordered; it warns about selfing, which it handles correctly
nadiv_additive <- function(p) {
  A <- suppressWarnings(nadiv::makeA(data.frame(ID = p$id, Dam = p$mom, Sire = p$dad)))
  A <- as.matrix(A)
  id <- as.character(p$id)
  A <- A[id, id, drop = FALSE]
  dimnames(A) <- list(id, id)
  A
}

test_that("the additive relationship matches nadiv on pedigrees of various shapes", {
  ## makeA() fails when no parent is known at all (a limitation of nadiv), so
  ## those pedigrees are checked against kinship2 and ggroups only
  peds <- validation_pedigrees()
  peds <- peds[!names(peds) %in% c("founders", "single")]
  for (nm in names(peds))
    expect_equal(fastrel_ped(peds[[nm]]), nadiv_additive(peds[[nm]]),
                 tolerance = 1e-12, label = nm)
})
