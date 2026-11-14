## The four measures of relatedness (see ?relatedness_metrics)

## G is a founder; P is obtained by selfing G; O1 and O2 by selfing P
selfing_ped <- function() {
  data.frame(id  = c("G", "P", "O1", "O2"),
             mom = c(NA, "G", "P", "P"),
             dad = c(NA, "G", "P", "P"))
}

test_that("each measure is the documented function of the additive relationship", {
  set.seed(4)
  ped <- sim_ped(500, n_founders = 4, window = 6, p_unknown = 0.1)   # inbred
  A <- fastrel_ped(ped)
  expect_equal(fastrel_ped(ped, metric = "additive"), A)
  expect_equal(fastrel_ped(ped, metric = "kinship"), A / 2)
  expect_equal(fastrel_ped(ped, metric = "wright"), cov2cor(A))
  expect_equal(fastrel_ped(ped, metric = "hamilton"), A / diag(A))   # row i / A[i,i]
  ## the same from vectors (metric must be named), and in any row order
  expect_equal(fastrel_vec(ped$id, ped$dadid, ped$momid, metric = "wright"), cov2cor(A))
  sh <- ped[sample(nrow(ped)), ]
  expect_equal(fastrel_ped(sh, metric = "hamilton")[ped$id, ped$id], A / diag(A))
  expect_equal(fastrel_ped(sh, metric = "wright")[ped$id, ped$id], cov2cor(A))
})

test_that("known values for repeated selfing", {
  ped <- selfing_ped()
  expect_equal(fastrel_ped(ped)["O1", "O2"], 1.5)
  expect_equal(fastrel_ped(ped, metric = "kinship")["O1", "O2"], 0.75)
  expect_equal(fastrel_ped(ped, metric = "wright")["O1", "O2"], 1.5 / 1.75)
  H <- fastrel_ped(ped, metric = "hamilton")
  expect_equal(H["G", "P"], 1)              # actor G (outbred), recipient P
  expect_equal(H["P", "G"], 1 / 1.5)        # actor P (inbred), recipient G
  expect_false(isSymmetric(H))
  expect_equal(unname(diag(H)), rep(1, 4))
  expect_equal(unname(diag(fastrel_ped(ped, metric = "wright"))), rep(1, 4))
  expect_equal(unname(diag(fastrel_ped(ped, metric = "kinship"))), c(1, 1.5, 1.75, 1.75) / 2)
})

test_that("without inbreeding, additive, wright and hamilton coincide", {
  p <- test_ped()[1:7, ]                    # no inbreeding among the first 7
  A <- fastrel_ped(p)
  expect_equal(fastrel_ped(p, metric = "wright"), A)
  expect_equal(fastrel_ped(p, metric = "hamilton"), A)
  expect_equal(A["E", "A"], 0.5)            # parent and offspring
  expect_equal(fastrel_ped(p, metric = "kinship")["E", "A"], 0.25)
})

test_that("bounds: nobody is more related to another than to itself", {
  for (s in 1:5) {
    set.seed(s)
    ped <- sim_ped(300, n_founders = 2, window = sample(3:8, 1), p_unknown = 0.2)
    A <- fastrel_ped(ped)
    d <- diag(A)
    expect_true(all(A <= pmin(outer(d, rep(1, length(d))), outer(rep(1, length(d)), d)) + 1e-12))
    expect_true(all(A >= 0 & A <= 2))
    K <- fastrel_ped(ped, metric = "kinship")
    expect_true(all(K >= 0 & K <= 1))
    for (m in c("wright", "hamilton")) {
      R <- fastrel_ped(ped, metric = m)
      expect_true(all(R >= 0 & R <= 1 + 1e-12))
      expect_equal(unname(diag(R)), rep(1, nrow(ped)))
    }
  }
})

test_that("rel_get returns every measure, for any subset and order", {
  set.seed(2)
  ped <- sim_ped(400, n_founders = 5, window = 10)
  rel <- rel_new()
  rel_add(rel, ped)
  sel <- sample(ped$id, 150)
  for (m in c("additive", "kinship", "wright", "hamilton"))
    expect_equal(rel_get(rel, sel, metric = m), fastrel_ped(ped, metric = m)[sel, sel])
  expect_equal(dim(rel_get(rel, character(0), metric = "wright")), c(0L, 0L))
  ## the store keeps the additive relationship whatever was extracted
  expect_equal(rel_get(rel, ped$id), fastrel_ped(ped))
  expect_equal(rel_get(rel_load(rel_save(rel)), ped$id, metric = "wright"),
               fastrel_ped(ped, metric = "wright"))
})

test_that("metric: partial matching, case and invalid values", {
  ped <- selfing_ped()
  A <- fastrel_ped(ped)
  expect_equal(fastrel_ped(ped, metric = "kin"), A / 2)
  expect_equal(fastrel_ped(ped, metric = "Wright"), cov2cor(A))
  expect_error(fastrel_ped(ped, metric = "foo"), "unknown `metric`")
  expect_error(fastrel_ped(ped, metric = c("kinship", "wright")), "must be one of")
  expect_error(fastrel_ped(ped, metric = NA), "must be one of")
  rel <- rel_new()
  rel_add(rel, ped)
  expect_error(rel_get(rel, metric = "foo"), "unknown `metric`")
  expect_equal(dim(fastrel_ped(ped[0, ], metric = "wright")), c(0L, 0L))
  ## a fourth positional argument, as for kinship2::kinship(id, dadid, momid, sex),
  ## is ignored rather than taken as the metric
  expect_equal(fastrel_vec(ped$id, ped$dad, ped$mom, c(1, 2, 1, 2)), A)
  expect_equal(fastrel_vec(ped$id, ped$dad, ped$mom, c("m", "f", "m", "f"), "autosome"), A)
  ## ... but a measure given without its name is refused, not ignored
  expect_error(fastrel_vec(ped$id, ped$dad, ped$mom, "kinship"), "must be named")
  expect_error(fastrel_vec(ped$id, ped$dad, ped$mom, c(1, 2, 1, 2), "Wright"), "must be named")
  expect_error(fastrel_vec(ped$id, ped$dad, ped$mom, "ham"), "must be named")
  expect_equal(fastrel_vec(ped$id, ped$dad, ped$mom, sex = "kinship"), A)   # named: not a measure
  ## internal calls pass it by name
  expect_equal(fastrel_ped(ped, "kinship"), A / 2)
})

test_that("kinship matches kinship2 when available", {
  skip_if_not_installed("kinship2")
  set.seed(1)
  ped <- sim_ped(200, n_founders = 20, window = 50)
  expect_equal(fastrel_ped(ped, metric = "kinship"),
               kinship2::kinship(ped$id, dadid = ped$dadid, momid = ped$momid))
})
