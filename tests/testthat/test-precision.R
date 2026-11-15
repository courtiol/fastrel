## fastrel_inv(): the inverse of the relatedness matrix, built directly

test_that("fastrel_inv() is the inverse of fastrel_ped() on pedigrees of various shapes", {
  peds <- validation_pedigrees()
  set.seed(4)
  for (nm in names(peds)) {
    p <- peds[[nm]]
    A <- fastrel_ped(p)
    Ainv <- fastrel_inv(p)
    expect_s4_class(Ainv, "dsCMatrix")
    expect_identical(dimnames(Ainv), dimnames(A), label = nm)
    ## relative to the largest entry, which is large in inbred pedigrees
    S <- solve(A)
    expect_lt(max(abs(as.matrix(Ainv) - S)) / max(abs(S)), 1e-10, label = nm)
    expect_lt(max(abs(as.matrix(Ainv) %*% A - diag(nrow(A)))), 1e-8, label = nm)
    ## whatever the order of the rows
    shuffled <- p[sample(nrow(p)), ]
    expect_equal(as.matrix(fastrel_inv(shuffled))[rownames(A), rownames(A), drop = FALSE],
                 as.matrix(Ainv), tolerance = 1e-10, label = paste(nm, "shuffled"))
  }
})

test_that("fastrel_inv() has at most 6 entries per individual in its upper triangle", {
  set.seed(5)
  p <- sim_ped(3000)
  Ainv <- fastrel_inv(p)
  expect_lte(length(Ainv@x), 6 * nrow(p))
  expect_equal(as.matrix(Ainv), solve(fastrel_ped(p)), tolerance = 1e-8)
})

test_that("fastrel_inv() handles selfing, unknown parents given as 0, and edge cases", {
  ## selfing: the parent appears twice
  p <- data.frame(id = c("G", "P", "O"), dad = c(NA, "G", "P"), mom = c(NA, "G", "P"))
  expect_equal(as.matrix(fastrel_inv(p)), solve(fastrel_ped(p)), tolerance = 1e-12)
  ## 0 for unknown parents, integer identifiers
  q <- data.frame(id = 1:4, dad = c(0, 0, 1, 1), mom = c(0, 0, 2, 2))
  expect_equal(as.matrix(fastrel_inv(q)), solve(fastrel_ped(q)), tolerance = 1e-12)
  ## founders only: the identity
  f <- fastrel_inv(data.frame(id = c("a", "b"), dad = NA, mom = NA))
  expect_equal(as.matrix(f), diag(2), ignore_attr = TRUE)
  ## empty pedigree
  expect_equal(dim(fastrel_inv(p[0, ])), c(0L, 0L))
  ## errors as for fastrel_ped()
  expect_error(fastrel_inv(data.frame(id = c("a", "a"), dad = NA, mom = NA)), "unique")
  expect_error(fastrel_inv(data.frame(id = c("a", "b"), dad = c("b", "a"), mom = NA)),
               "own ancestor")
})

test_that("fastrel_inv() matches nadiv::makeAinv()", {
  skip_if_not_installed("nadiv")
  set.seed(6)
  p <- sim_ped(1000)
  N <- suppressWarnings(nadiv::makeAinv(data.frame(ID = p$id, Dam = p$momid, Sire = p$dadid))$Ainv)
  N <- as.matrix(N)
  colnames(N) <- rownames(N)          # makeAinv() names the rows only
  expect_equal(N[p$id, p$id], as.matrix(fastrel_inv(p)), tolerance = 1e-10)
})
