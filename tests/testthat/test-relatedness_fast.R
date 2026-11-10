test_that("known values for the mother-son pedigree", {
  p <- test_ped()
  A <- relatedness_fast(p$id, dadid = p$dad, momid = p$mom)
  expect_equal(unname(diag(A) - 1), c(rep(0, 7), 0.25, 0.25))
  expect_equal(A["H", "E"], 0.75)
  expect_equal(A["H", "I"], 0.75)
  expect_equal(A["E", "G"], 0.5)
  expect_equal(A["H", "F"], 0.25)
  expect_true(isSymmetric(A))
})

test_that("matches the reference implementation, whatever the row order", {
  ped <- sim_ped(300, n_founders = 20, window = 60)
  A0 <- ref_relatedness(ped$id, ped$dadid, ped$momid)
  expect_equal(relatedness_fast(ped$id, ped$dadid, ped$momid), A0)

  set.seed(10)
  shuffled <- ped[sample(nrow(ped)), ]
  A1 <- relatedness_fast(shuffled$id, shuffled$dadid, shuffled$momid)
  expect_equal(A1[ped$id, ped$id], A0)
})

test_that("matches kinship2 when available", {
  skip_if_not_installed("kinship2")
  ped <- sim_ped(200, n_founders = 20, window = 50)
  expect_equal(relatedness_fast(ped$id, ped$dadid, ped$momid),
               2 * kinship2::kinship(ped$id, dadid = ped$dadid, momid = ped$momid))
})

test_that("unknown parents and unrelated individuals", {
  A <- relatedness_fast(c("a", "b", "c"), dadid = c(NA, NA, "zzz"), momid = c(NA, NA, "a"))
  expect_equal(unname(A["c", "a"]), 0.5)   # one known parent
  expect_equal(unname(A["c", "b"]), 0)
  expect_equal(unname(diag(A)), c(1, 1, 1))
})

test_that("input checks and edge cases", {
  expect_error(relatedness_fast(c("a", "a"), c(NA, NA), c(NA, NA)), "unique")
  expect_error(relatedness_fast(c("a", "b"), NA, c(NA, NA)), "father")
  expect_error(relatedness_fast(c("a", "b"), c(NA, NA), NA), "mother")
  expect_equal(dim(relatedness_fast(character(0), character(0), character(0))), c(0L, 0L))
  expect_equal(unname(relatedness_fast("a", NA, NA)), matrix(1, 1, 1))
  ## a cycle: a is the child of b, b is the child of a
  expect_error(relatedness_fast(c("a", "b"), dadid = c("b", "a"), momid = c(NA, NA)),
               "own ancestor")
})
