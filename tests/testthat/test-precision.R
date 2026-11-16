## fastrel_inv(): the inverse of the relatedness matrix, built directly

test_that("fastrel_inv() is the inverse of fastrel_ped() on pedigrees of various shapes", {
  peds <- validation_pedigrees()
  set.seed(4)
  for (nm in names(peds)) {
    p <- peds[[nm]]
    A <- fastrel_ped(p)
    Ainv <- fastrel_inv(p)
    ## both triangles stored (MCMCglmm misreads a matrix storing only one)
    expect_s4_class(Ainv, "dgCMatrix")
    expect_true(Matrix::isSymmetric(Ainv))
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

test_that("fastrel_inv(matrix_class = \"dsCMatrix\") gives the same matrix, one triangle stored", {
  set.seed(7)
  p <- sim_ped(500)
  general <- fastrel_inv(p)
  symmetric <- fastrel_inv(p, matrix_class = "dsCMatrix")
  expect_s4_class(symmetric, "dsCMatrix")
  expect_identical(dimnames(symmetric), dimnames(general))
  expect_equal(as.matrix(symmetric), as.matrix(general), tolerance = 1e-12)
  entries <- function(m) as.numeric(object.size(m@x) + object.size(m@i) + object.size(m@p))
  expect_lt(entries(symmetric), 0.7 * entries(general))
  expect_s4_class(fastrel_inv(p[0, ], matrix_class = "dsCMatrix"), "dsCMatrix")
  expect_error(fastrel_inv(p, matrix_class = "dense"))
})

test_that("fastrel_inv() has at most 6 entries per individual in its upper triangle", {
  set.seed(5)
  p <- sim_ped(3000)
  Ainv <- fastrel_inv(p)
  expect_lte(Matrix::nnzero(Matrix::triu(Ainv)), 6 * nrow(p))
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
