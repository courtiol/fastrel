## Validation against ggroups::buildA(), an independent implementation of the
## additive (numerator) relationship matrix used in animal breeding, on the
## same pedigrees as for kinship2 (see helper-pedigrees.R). Only the additive
## relationship is compared: the other measures are derived from it in
## fastrel, and are checked against kinship2 (test-kinship2.R). ggroups needs
## the parents to come before their offspring, and codes unknown parents as 0.

skip_if_not_installed("ggroups")

## ggroups wants integer identifiers, with 0 for an unknown parent: the
## individuals are numbered in birth order
ggroups_additive <- function(p) {
  id  <- as.character(p$id)
  num <- function(x) { k <- match(as.character(x), id); ifelse(is.na(k), 0L, k) }
  A <- as.matrix(ggroups::buildA(data.frame(ID = seq_along(id), SIRE = num(p$dad),
                                            DAM = num(p$mom))))
  dimnames(A) <- list(id, id)
  A
}

test_that("the additive relationship matches ggroups on pedigrees of various shapes", {
  ## buildA() fails on a pedigree of a single individual (a limitation of
  ## ggroups), so that one is checked against kinship2 only
  peds <- validation_pedigrees()
  peds <- peds[names(peds) != "single"]
  for (nm in names(peds))
    expect_equal(fastrel_ped(peds[[nm]]), ggroups_additive(peds[[nm]]),
                 tolerance = 1e-12, label = nm)
})
