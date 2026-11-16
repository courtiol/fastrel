## fastrel_inv() used as it is in MCMCglmm, to fit an animal model: from the
## same random seed, the chains must be the same as with the inverse computed
## by MCMCglmm itself, MCMCglmm::inverseA().

skip_if_not_installed("MCMCglmm")

test_that("MCMCglmm gives the same chains with fastrel_inv() as with inverseA()", {
  set.seed(1)
  p <- sim_ped(300)
  A <- fastrel_ped(p)
  Ainv <- fastrel_inv(p)
  ## breeding values with additive variance 2, residual variance 1
  a <- drop(crossprod(chol(A), rnorm(nrow(p))) * sqrt(2))
  d <- data.frame(animal = p$id, y = 10 + a + rnorm(nrow(p)))

  ## the two inverses are the same matrix
  ref <- MCMCglmm::inverseA(data.frame(animal = p$id, dam = p$momid,
                                       sire = p$dadid))$Ainv
  expect_s4_class(ref, "dgCMatrix")
  k <- match(p$id, rownames(ref))     # inverseA() names the rows only
  expect_equal(unname(as.matrix(ref[k, k])), unname(as.matrix(Ainv)), tolerance = 1e-10)

  prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                R = list(V = 1, nu = 0.002))
  fit <- function(Ginv, data = d) {
    set.seed(42)
    MCMCglmm::MCMCglmm(y ~ 1, random = ~animal, ginverse = list(animal = Ginv),
                       data = data, prior = prior, nitt = 3000, burnin = 1000,
                       thin = 10, verbose = FALSE)
  }
  m_fastrel <- fit(Ainv)
  m_ref <- fit(ref)
  expect_equal(unclass(m_fastrel$VCV), unclass(m_ref$VCV), tolerance = 1e-8)
  expect_equal(unclass(m_fastrel$Sol), unclass(m_ref$Sol), tolerance = 1e-8)

  ## and the estimates are sensible (true values: 2 and 1)
  means <- colMeans(m_fastrel$VCV)
  expect_gt(means[["animal"]], 0.5)
  expect_gt(means[["units"]], 0.3)
})
