## A tibble differs from a data frame in ways that matter here: x[, "id"]
## returns a tibble with one column instead of a vector, and so does x[1, "id"].

test_that("fastrel_vec and fastrel_ped work with tibbles", {
  skip_if_not_installed("tibble")
  p <- test_ped()
  A0 <- fastrel_vec(p$id, dadid = p$dad, momid = p$mom)
  tb <- tibble::as_tibble(p)

  ## columns extracted from a tibble as vectors
  expect_equal(fastrel_vec(tb$id, dadid = tb$dad, momid = tb$mom), A0)
  expect_equal(fastrel_vec(tb[["id"]], tb[["dad"]], tb[["mom"]]), A0)
  ## columns extracted with [, ]: tibbles with one column
  expect_equal(fastrel_vec(tb[, "id"], dadid = tb[, "dad"], momid = tb[, "mom"]), A0)
  expect_equal(fastrel_vec(tb["id"], tb["dad"], tb["mom"]), A0)
  ## a mix of both
  expect_equal(fastrel_vec(tb$id, dadid = tb[, "dad"], momid = p$mom), A0)
  ## the same goes for data frames with one column
  expect_equal(fastrel_vec(p["id"], p["dad"], p["mom"]), A0)

  ## the whole tibble
  expect_equal(fastrel_ped(tb), A0)
  expect_equal(fastrel_ped(tb[, c("dad", "id", "mom")]), A0)
  expect_equal(fastrel_ped(tibble::tibble(Sire = p$dad, Dam = p$mom, Animal = p$id)), A0)
  set.seed(1)
  ped <- tibble::as_tibble(sim_ped(200))
  expect_equal(fastrel_ped(ped), fastrel_vec(ped$id, ped$dadid, ped$momid))
  expect_equal(dim(fastrel_ped(tb[0, ])), c(0L, 0L))

  ## what cannot be a vector of identifiers is refused, not misread
  expect_error(fastrel_vec(tb, tb$dad, tb$mom), "`id` must be a vector")
  expect_error(fastrel_vec(tb$id, tb[, c("dad", "mom")], tb$mom), "`dadid` must be a vector")
  expect_error(fastrel_vec(as.list(p$id), p$dad, p$mom), "`id` must be a vector")
})

test_that("the store works with tibbles", {
  skip_if_not_installed("tibble")
  p <- test_ped()
  A0 <- fastrel_ped(p)
  tb <- tibble::as_tibble(p)
  rel <- rel_new()
  rel_add(rel, tb[1:4, ])
  rel_add(rel, tb[5:9, c("dad", "mom", "id")])
  expect_equal(rel_get(rel, tb$id), A0)
  expect_equal(rel_get(rel, tb[, "id"]), A0)       # a tibble with one column
  keep <- tb[tb$id %in% c("E", "H", "I"), "id"]
  expect_equal(rel_get(rel, keep), A0[keep$id, keep$id])
  rel_remove(rel, tb[1:2, "id"])
  expect_identical(rel_ids(rel), p$id[3:9])
  expect_error(rel_get(rel, tb), "`ids` must be a vector")
  expect_error(rel_remove(rel, tb[, c("id", "mom")]), "`ids` must be a vector")
  expect_identical(rel_ids(rel), p$id[3:9])
})
