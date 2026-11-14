## Validation against ggroups::buildA(), an independent implementation of the
## additive (numerator) relationship matrix used in animal breeding, on the
## same pedigrees as for kinship2 (see helper-pedigrees.R). ggroups needs the
## parents to come before their offspring, and codes unknown parents as 0.

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

## buildA() fails on a pedigree of a single individual (a limitation of
## ggroups), so that one is checked against kinship2 only
ggroups_pedigrees <- function() {
  peds <- validation_pedigrees()
  peds[names(peds) != "single"]
}

test_that("all measures match ggroups on pedigrees of various shapes", {
  peds <- ggroups_pedigrees()
  set.seed(2)
  for (nm in names(peds))
    expect_matches_reference(peds[[nm]], ggroups_additive(peds[[nm]]), what = nm)
})

test_that("ggroups and kinship2 agree with each other on these pedigrees", {
  ## a check of the references themselves: if they disagreed, a pass against
  ## one of them would mean little
  skip_if_not_installed("kinship2")
  peds <- ggroups_pedigrees()
  for (nm in names(peds)) {
    p <- peds[[nm]]
    expect_equal(unname(ggroups_additive(p)),
                 unname(2 * kinship2::kinship(p$id, dadid = p$dad, momid = p$mom)),
                 tolerance = 1e-12, label = nm)
  }
})
