test_that("store reproduces relatedness_fast, any row order, with growth", {
  ped <- sim_ped(400, n_founders = 20, window = 60)
  rel <- rel_new(8)                          # forces several doublings
  set.seed(3)
  rel_add(rel, ped[sample(nrow(ped)), c("id", "momid", "dadid")])
  expect_equal(rel_size(rel), 400L)
  expect_gte(rel_capacity(rel), 400L)
  expect_equal(rel_get(rel, ped$id),
               relatedness_fast(ped$id, ped$dadid, ped$momid))
})

test_that("incremental additions and removals stay exact (slot reuse)", {
  set.seed(42)
  ped <- sim_ped(600, n_founders = 20, window = 60)
  full <- relatedness_fast(ped$id, ped$dadid, ped$momid)
  rel <- rel_new(40)
  rel_add(rel, ped[1:30, c("id", "momid", "dadid")])
  batch <- 10
  for (start in seq(31, 591, by = batch)) {
    rows <- start:(start + batch - 1)
    rel_add(rel, ped[rows, c("id", "momid", "dadid")])
    ## remove only individuals that never parent again
    later <- ped[-seq_len(max(rows)), ]
    future_parents <- unique(stats::na.omit(c(later$momid, later$dadid)))
    candidates <- setdiff(rel$ids, future_parents)
    die <- head(sample(candidates), 8)
    rel_remove(rel, die)
  }
  alive <- rel$ids
  expect_lt(length(alive), 600)
  expect_equal(rel_get(rel), full[alive, alive])
})

test_that("known values from the documentation example", {
  rel <- rel_new(16)
  rel_add(rel, data.frame(id = c("A", "B", "C", "D"), mom = NA, dad = NA))
  rel_add(rel, data.frame(id = c("E", "F"), mom = c("A", "C"), dad = c("B", "D")))
  rel_add(rel, data.frame(id = "G", mom = "E", dad = "F"))
  rel_remove(rel, c("A", "B", "C", "D"))
  rel_add(rel, data.frame(id = "H", mom = "E", dad = "G"))
  A <- rel_get(rel)
  expect_equal(unname(diag(A) - 1), c(0, 0, 0, 0.25))
  expect_equal(rel_size(rel), 4L)
})

test_that("error handling", {
  rel <- rel_new(4)
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  expect_error(rel_add(rel, data.frame(id = "a", mom = NA, dad = NA)), "already")
  expect_error(rel_add(rel, data.frame(id = c("x", "x"), mom = NA, dad = NA)), "duplicated")
  expect_error(rel_add(rel, data.frame(id = "c", mom = "zz", dad = "a")), "not found")
  expect_error(rel_add(rel, data.frame(id = NA_character_, mom = NA, dad = NA)), "NA")
  expect_error(rel_add(rel, data.frame(id = c("p", "q"), mom = c("q", "p"), dad = NA)), "cycle")
  expect_error(rel_remove(rel, "nobody"), "unknown")
  expect_error(rel_get(rel, "nobody"), "unknown")
  expect_error(rel_new(0), "positive")
  expect_error(rel_add(list(), data.frame(id = "z", mom = NA, dad = NA)), "store")
  ## a failed call leaves the store unchanged
  expect_equal(rel_size(rel), 2L)
  expect_equal(rel$ids, c("a", "b"))
})

test_that("removing nothing, or a parent, behaves", {
  rel <- rel_new(4)
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  rel_remove(rel, character(0))
  expect_equal(rel$ids, c("a", "b"))
  rel_remove(rel, "a")
  expect_error(rel_add(rel, data.frame(id = "c", mom = "a", dad = "b")), "not found")
  expect_equal(rel_size(rel), 1L)
})

test_that("print method works", {
  rel <- rel_new(4)
  expect_output(print(rel), "relatedness_store")
})

test_that("a cycle found after some insertions leaves the store unchanged", {
  rel <- rel_new(2)
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  before <- rel_get(rel)
  ## f and c can be inserted, p and q form a cycle
  expect_error(rel_add(rel, data.frame(id  = c("f", "c", "p", "q"),
                                       mom = c(NA, "f", "q", "p"),
                                       dad = c(NA, "a", NA, NA))), "cycle")
  expect_equal(rel_size(rel), 2L)
  expect_equal(rel$ids, c("a", "b"))
  expect_equal(rel_get(rel), before)
  ## the store still works, and c is inserted after its mother f
  rel_add(rel, data.frame(id = c("c", "f"), mom = c("f", "a"), dad = c("b", NA)))
  expect_equal(rel$ids, c("a", "b", "f", "c"))
  expect_equal(rel_get(rel, c("a", "b", "f", "c")),
               relatedness_fast(c("a", "b", "f", "c"), dadid = c(NA, NA, NA, "b"),
                                momid = c(NA, NA, "a", "f")))
})

test_that("an oversized or emptied store gives the same values", {
  ped <- sim_ped(120, n_founders = 10, window = 30)
  full <- relatedness_fast(ped$id, ped$dadid, ped$momid)
  rel <- rel_new(5000)
  rel_add(rel, ped[, c("id", "momid", "dadid")])
  expect_equal(rel_get(rel, ped$id), full)
  expect_equal(rel_capacity(rel), 5000L)
  ## empty it completely and start again: no value must be left behind
  rel_remove(rel, ped$id)
  expect_equal(rel_size(rel), 0L)
  expect_equal(dim(rel_get(rel)), c(0L, 0L))
  rel_add(rel, ped[, c("id", "momid", "dadid")])
  expect_equal(rel_get(rel, ped$id), full)
})

test_that("adding nothing and invalid capacities", {
  rel <- rel_new(4)
  expect_identical(rel_add(rel, data.frame(id = character(0), mom = character(0),
                                           dad = character(0))), rel)
  expect_equal(rel_size(rel), 0L)
  expect_error(rel_new("a"), "positive")
  expect_error(rel_new(1e10), "positive")
  expect_error(rel_new(c(2, 3)), "positive")
})
