test_that("save/load round trip, in memory and via file", {
  set.seed(1)
  ped <- sim_ped(200, n_founders = 20, window = 50)
  rel <- rel_new()
  rel_add(rel, ped[, c("id", "dadid", "momid")])
  rel_remove(rel, ped$id[1:60])                  # slots are recycled / holes exist

  obj <- rel_save(rel)
  expect_s3_class(obj, "relatedness_store_saved")
  rel2 <- rel_load(obj)
  expect_equal(rel_get(rel2), rel_get(rel))
  expect_equal(rel_size(rel2), rel_size(rel))
  expect_gte(rel_capacity(rel2), rel_size(rel2))

  f <- tempfile(fileext = ".rds")
  expect_equal(rel_save(rel, f), f)
  rel3 <- rel_load(f)
  expect_equal(rel_get(rel3), rel_get(rel))
})

test_that("a restored store keeps working and stays exact", {
  set.seed(1)
  ped <- sim_ped(300, n_founders = 20, window = 60)
  full <- fastrel_vec(ped$id, ped$dadid, ped$momid)
  rel <- rel_new()
  rel_add(rel, ped[1:150, c("id", "dadid", "momid")])
  rel2 <- rel_load(rel_save(rel))
  rel_add(rel2, ped[151:300, c("id", "dadid", "momid")])
  expect_equal(rel_get(rel2, ped$id), full)
  rel_remove(rel2, ped$id[1:100])
  keep <- ped$id[101:300]      # rel_get() default order is insertion order
  expect_equal(rel_get(rel2, keep), full[keep, keep])
  ## the original is independent of the copy
  expect_equal(rel_size(rel), 150L)
})

test_that("empty stores and invalid input", {
  rel <- rel_new()
  rel2 <- rel_load(rel_save(rel))
  expect_equal(rel_size(rel2), 0L)
  rel_add(rel2, data.frame(id = "a", mom = NA, dad = NA))
  expect_equal(unname(rel_get(rel2)), matrix(1, 1, 1))

  expect_error(rel_load(list()), "rel_save")
  expect_error(rel_load(1:3), "rel_save")
  bad <- rel_save(rel2); bad$ids <- c("a", "b")
  expect_error(rel_load(bad), "match")
  bad <- rel_save(rel2); bad$A[1, 1] <- NA
  expect_error(rel_load(bad), "non-finite")
  expect_error(rel_save(list()), "store")
  bad <- rel_save(rel2); bad$remember_removed <- NULL
  expect_error(rel_load(bad), "rel_save")
})

test_that("the saved object records its format and the version of fastrel", {
  obj <- rel_save(rel_new())
  expect_identical(obj$format, 1L)
  expect_identical(obj$fastrel_version, as.character(utils::packageVersion("fastrel")))
  ## another version of fastrel in the same format: read as usual
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a", "b", "c"), dad = c(NA, NA, "a"), mom = c(NA, NA, "b")))
  other <- rel_save(rel)
  other$fastrel_version <- "9.9.9"
  expect_equal(rel_get(rel_load(other)), rel_get(rel))
  ## an object in another format is refused, naming the version that saved it
  future <- obj
  future$format <- 2L
  future$fastrel_version <- "9.9.9"
  expect_error(rel_load(future), "saved by fastrel 9.9.9 in a format.*update fastrel")
  ## ... even if its structure differs entirely
  future$ids <- NULL
  expect_error(rel_load(future), "saved by fastrel 9.9.9 in a format")
  bad <- obj
  bad$format <- NULL
  expect_error(rel_load(bad), "cannot read")
  bad <- obj
  bad$fastrel_version <- NULL
  expect_error(rel_load(bad), "rel_save")
})
