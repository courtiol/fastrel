## fastrel_inv() used as it is in spaMM, to fit an animal model: the fit with
## the precision matrix must be the same as with the relationship matrix.

skip_if_not_installed("spaMM")

test_that("spaMM fits the same animal model with fastrel_inv() as with fastrel_ped()", {
  set.seed(1)
  p <- sim_ped(300)
  A <- fastrel_ped(p)
  Ainv <- fastrel_inv(p)
  ## breeding values with additive variance 2, residual variance 1
  a <- drop(crossprod(chol(A), rnorm(nrow(p)))) * sqrt(2)
  d <- data.frame(id = p$id, y = 10 + a + rnorm(nrow(p)))

  ## the precision matrix is matched to the data by its names: no warning
  expect_no_warning(
    fit_prec <- spaMM::fitme(y ~ 1 + corrMatrix(1 | id),
                             covStruct = list(precision = Ainv), data = d))
  fit_dense <- suppressMessages(
    spaMM::fitme(y ~ 1 + corrMatrix(1 | id), covStruct = list(corrMatrix = A), data = d))

  expect_equal(logLik(fit_prec), logLik(fit_dense), tolerance = 1e-8)
  expect_equal(unlist(spaMM::get_ranPars(fit_prec, "lambda")),
               unlist(spaMM::get_ranPars(fit_dense, "lambda")), tolerance = 1e-4)
  expect_equal(fit_prec$phi[1], fit_dense$phi[1], tolerance = 1e-4)
  expect_equal(spaMM::fixef(fit_prec), spaMM::fixef(fit_dense), tolerance = 1e-6)
  u_prec <- spaMM::ranef(fit_prec)[[1]]
  u_dense <- spaMM::ranef(fit_dense)[[1]]
  expect_setequal(names(u_prec), p$id)
  expect_equal(u_prec[p$id], u_dense[p$id], tolerance = 1e-6)

  ## the rows of the data can be in any order
  fit_shuffled <- spaMM::fitme(y ~ 1 + corrMatrix(1 | id),
                               covStruct = list(precision = Ainv),
                               data = d[sample(nrow(d)), ])
  expect_equal(logLik(fit_shuffled), logLik(fit_prec), tolerance = 1e-8)

  ## spaMM also reads the symmetric class, which stores one triangle only
  fit_sym <- spaMM::fitme(y ~ 1 + corrMatrix(1 | id),
                          covStruct = list(precision = fastrel_inv(p, matrix_class = "dsCMatrix")),
                          data = d)
  expect_equal(logLik(fit_sym), logLik(fit_prec), tolerance = 1e-8)
})
