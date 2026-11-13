test_that("known values for the mother-son pedigree", {
  p <- test_ped()
  A <- fastrel_vec(p$id, dadid = p$dad, momid = p$mom)
  expect_equal(unname(diag(A) - 1), c(rep(0, 7), 0.25, 0.25))
  expect_equal(A["H", "E"], 0.75)
  expect_equal(A["H", "I"], 0.75)
  expect_equal(A["E", "G"], 0.5)
  expect_equal(A["H", "F"], 0.25)
  expect_true(isSymmetric(A))
})

test_that("matches the reference implementation, whatever the row order", {
  set.seed(1)
  ped <- sim_ped(300, n_founders = 20, window = 60)
  A0 <- ref_relatedness(ped$id, ped$dadid, ped$momid)
  expect_equal(fastrel_vec(ped$id, ped$dadid, ped$momid), A0)

  set.seed(10)
  shuffled <- ped[sample(nrow(ped)), ]
  A1 <- fastrel_vec(shuffled$id, shuffled$dadid, shuffled$momid)
  expect_equal(A1[ped$id, ped$id], A0)
})

test_that("matches kinship2 when available", {
  skip_if_not_installed("kinship2")
  set.seed(1)
  ped <- sim_ped(200, n_founders = 20, window = 50)
  expect_equal(fastrel_vec(ped$id, ped$dadid, ped$momid),
               2 * kinship2::kinship(ped$id, dadid = ped$dadid, momid = ped$momid))
})

test_that("unknown parents and unrelated individuals", {
  A <- fastrel_vec(c("a", "b", "c"), dadid = c(NA, NA, "zzz"), momid = c(NA, NA, "a"))
  expect_equal(unname(A["c", "a"]), 0.5)   # one known parent
  expect_equal(unname(A["c", "b"]), 0)
  expect_equal(unname(diag(A)), c(1, 1, 1))
})

test_that("input checks and edge cases", {
  expect_error(fastrel_vec(c("a", "a"), c(NA, NA), c(NA, NA)), "unique")
  expect_error(fastrel_vec(c("a", "b"), NA, c(NA, NA)), "father")
  expect_error(fastrel_vec(c("a", "b"), c(NA, NA), NA), "mother")
  expect_equal(dim(fastrel_vec(character(0), character(0), character(0))), c(0L, 0L))
  expect_equal(unname(fastrel_vec("a", NA, NA)), matrix(1, 1, 1))
  ## a cycle: a is the child of b, b is the child of a
  expect_error(fastrel_vec(c("a", "b"), dadid = c("b", "a"), momid = c(NA, NA)),
               "own ancestor")
})

test_that("rows in reverse birth order, and NA in id", {
  set.seed(1)
  ped <- sim_ped(150, n_founders = 5, window = 20, p_unknown = 0.3)
  A0 <- ref_relatedness(ped$id, ped$dadid, ped$momid)
  rev_ped <- ped[nrow(ped):1, ]
  A1 <- fastrel_vec(rev_ped$id, rev_ped$dadid, rev_ped$momid)
  expect_identical(rownames(A1), rev_ped$id)
  expect_equal(A1[ped$id, ped$id], A0)
  expect_true(isSymmetric(A1))
  expect_error(fastrel_vec(c("a", NA), c(NA, NA), c(NA, NA)), "NA")
  expect_error(fastrel_vec(c("a", NA, "c"), dadid = c(NA, "zz", NA),
                                momid = c(NA, "zz", NA)), "NA")
})

test_that("selfing, a self-parent, factor or numeric ids, and chrtype", {
  ## 2 is the offspring of 1 selfing; 3 is a backcross of 2 to 1
  A <- fastrel_vec(c(1, 2, 3), dadid = c(NA, 1, 1), momid = c(NA, 1, 2))
  expect_equal(unname(A), matrix(c(1, 1, 1, 1, 1.5, 1.25, 1, 1.25, 1.5), 3))
  expect_identical(rownames(A), c("1", "2", "3"))
  p <- test_ped()
  A0 <- fastrel_vec(p$id, dadid = p$dad, momid = p$mom)
  expect_equal(fastrel_vec(factor(p$id), dadid = factor(p$dad), momid = factor(p$mom)), A0)
  expect_error(fastrel_vec("a", dadid = "a", momid = NA), "own ancestor")
  expect_equal(fastrel_vec(p$id, dadid = p$dad, momid = p$mom, sex = 1,
                                chrtype = "autosome"), A0)
  expect_error(fastrel_vec(p$id, dadid = p$dad, momid = p$mom, chrtype = "X"),
               "autosome")
})

test_that("fastrel_ped gives the same as fastrel_vec, whatever the columns", {
  p <- test_ped()
  A0 <- fastrel_vec(p$id, dadid = p$dad, momid = p$mom)
  expect_equal(fastrel_ped(p), A0)
  ## any order, any capitalisation, usual synonyms, extra columns
  expect_equal(fastrel_ped(data.frame(weight = 1:9, Sire = p$dad, DAM = p$mom, Animal = p$id)), A0)
  expect_equal(fastrel_ped(data.frame(Father = p$dad, focal = p$id, Mother = p$mom)), A0)
  expect_equal(fastrel_ped(data.frame(mum = p$mom, dad = p$dad, ID = p$id)), A0)
  expect_equal(fastrel_ped(data.frame(dad_id = p$dad, Mom.ID = p$mom, ind_id = p$id)), A0)
  expect_equal(fastrel_ped(data.frame(IID = p$id, PAT = p$dad, MAT = p$mom)), A0)
  set.seed(1)
  ped <- sim_ped(100)                            # id, momid, dadid, sex
  expect_equal(fastrel_ped(ped), fastrel_vec(ped$id, ped$dadid, ped$momid))
  ## no recognised name at all: the first three columns, in this order
  expect_equal(fastrel_ped(data.frame(x = p$id, y = p$mom, z = p$dad)), A0)
  expect_equal(dim(fastrel_ped(p[0, ])), c(0L, 0L))
})

test_that("fastrel_ped refuses ambiguous or incomplete columns", {
  p <- test_ped()
  expect_error(fastrel_ped(data.frame(id = p$id, mom = p$mom, x = p$dad)), "father")
  expect_error(fastrel_ped(data.frame(id = p$id, m = p$mom, d = p$dad)),
               "mother or as the father")
  expect_error(fastrel_ped(data.frame(x = p$id, mom = p$mom, dad = p$dad)), "focal")
  expect_error(fastrel_ped(data.frame(id = p$id, mom = p$mom, mother = p$mom, dad = p$dad)),
               "several columns.*mother.*mom, mother")
  expect_error(fastrel_ped(data.frame(id = p$id, animal = p$id, mom = p$mom, dad = p$dad)),
               "several columns.*focal")
  expect_error(fastrel_ped(p[, 1:2]), "data frame")
  expect_error(fastrel_ped(as.matrix(p)), "data frame")
  expect_error(relatedness_fast(p$id, p$dad, p$mom), "could not find function")
})

test_that("range of values stated in the documentation", {
  ## repeated selfing: the diagonal is 2 - (1/2)^n and approaches 2
  n <- 12
  id <- seq_len(n + 1)
  A <- fastrel_vec(id, dadid = c(NA, id[-(n + 1)]), momid = c(NA, id[-(n + 1)]))
  expect_equal(unname(diag(A)), 2 - 0.5^(0:n))
  expect_lt(max(A), 2)
  ## two offspring of the same selfed parent are related by its diagonal value
  B <- fastrel_vec(1:4, dadid = c(NA, 1, 2, 2), momid = c(NA, 1, 2, 2))
  expect_equal(B["3", "4"], B["2", "2"])
  expect_gt(B["3", "4"], 1)
  ## in general: 0 <= A[i,j] <= sqrt(A[i,i] * A[j,j]) <= 2, so cov2cor(A) <= 1
  set.seed(4)
  ped <- sim_ped(400, n_founders = 4, window = 6, p_unknown = 0)   # much inbreeding
  A <- fastrel_ped(ped)
  expect_gt(max(diag(A)), 1.5)
  expect_gt(max(A[upper.tri(A)]), 1)
  expect_true(all(A >= 0 & A <= 2))
  expect_true(all(A <= sqrt(outer(diag(A), diag(A))) + 1e-12))
  expect_lte(max(cov2cor(A)), 1 + 1e-12)
  ## without inbreeding, relatedness between two individuals is at most 1
  p <- test_ped()[1:7, ]
  A <- fastrel_ped(p)
  expect_equal(unname(diag(A)), rep(1, 7))
  expect_lte(max(A[upper.tri(A)]), 0.5)
})
