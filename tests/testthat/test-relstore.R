test_that("store reproduces relatedness_fast, any row order, with growth", {
  set.seed(1)
  ped <- sim_ped(400, n_founders = 20, window = 60)
  rel <- rel_new()
  set.seed(3)
  rel_add(rel, ped[sample(nrow(ped)), c("id", "momid", "dadid")])
  expect_equal(rel_size(rel), 400L)
  expect_gte(rel_capacity(rel), 400L)
  expect_equal(rel_get(rel, ped$id),
               relatedness_fast(ped$id, ped$dadid, ped$momid))
})

test_that("incremental additions and removals stay exact (slot reuse)", {
  set.seed(1)
  ped <- sim_ped(600, n_founders = 20, window = 60)
  full <- relatedness_fast(ped$id, ped$dadid, ped$momid)
  rel <- rel_new()
  rel_add(rel, ped[1:30, c("id", "momid", "dadid")])
  batch <- 10
  for (start in seq(31, 591, by = batch)) {
    rows <- start:(start + batch - 1)
    rel_add(rel, ped[rows, c("id", "momid", "dadid")])
    ## remove only individuals that never parent again
    later <- ped[-seq_len(max(rows)), ]
    future_parents <- unique(stats::na.omit(c(later$momid, later$dadid)))
    candidates <- setdiff(rel_ids(rel), future_parents)
    die <- head(sample(candidates), 8)       # random: hence the set.seed() above
    rel_remove(rel, die)
  }
  alive <- rel_ids(rel)
  expect_lt(length(alive), 600)
  expect_equal(rel_get(rel), full[alive, alive])
})

test_that("known values from the documentation example", {
  rel <- rel_new()
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
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  expect_error(rel_add(rel, data.frame(id = "a", mom = NA, dad = NA)), "already")
  expect_error(rel_add(rel, data.frame(id = c("x", "x"), mom = NA, dad = NA)), "duplicated")
  expect_error(rel_add(rel, data.frame(id = "c", mom = "zz", dad = "a")), "not found")
  expect_error(rel_add(rel, data.frame(id = NA_character_, mom = NA, dad = NA)), "NA")
  expect_error(rel_add(rel, data.frame(id = c("p", "q"), mom = c("q", "p"), dad = NA)), "cycle")
  expect_error(rel_remove(rel, "nobody"), "unknown")
  expect_error(rel_get(rel, "nobody"), "unknown")
  expect_error(rel_add(list(), data.frame(id = "z", mom = NA, dad = NA)), "store")
  ## a failed call leaves the store unchanged
  expect_equal(rel_size(rel), 2L)
  expect_equal(rel_ids(rel), c("a", "b"))
})

test_that("removing nothing, or a parent, behaves", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  rel_remove(rel, character(0))
  expect_equal(rel_ids(rel), c("a", "b"))
  rel_remove(rel, "a")
  expect_error(rel_add(rel, data.frame(id = "c", mom = "a", dad = "b")), "not found")
  expect_equal(rel_size(rel), 1L)
})

test_that("print method works", {
  rel <- rel_new()
  expect_output(print(rel), "relatedness_store")
})

test_that("a cycle found after some insertions leaves the store unchanged", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  before <- rel_get(rel)
  ## f and c can be inserted, p and q form a cycle
  expect_error(rel_add(rel, data.frame(id  = c("f", "c", "p", "q"),
                                       mom = c(NA, "f", "q", "p"),
                                       dad = c(NA, "a", NA, NA))), "cycle")
  expect_equal(rel_size(rel), 2L)
  expect_equal(rel_ids(rel), c("a", "b"))
  expect_equal(rel_get(rel), before)
  ## the store still works, and c is inserted after its mother f
  rel_add(rel, data.frame(id = c("c", "f"), mom = c("f", "a"), dad = c("b", NA)))
  expect_equal(rel_ids(rel), c("a", "b", "f", "c"))
  expect_equal(rel_get(rel, c("a", "b", "f", "c")),
               relatedness_fast(c("a", "b", "f", "c"), dadid = c(NA, NA, NA, "b"),
                                momid = c(NA, NA, "a", "f")))
})

test_that("an oversized or emptied store gives the same values", {
  set.seed(1)
  ped <- sim_ped(120, n_founders = 10, window = 30)
  full <- relatedness_fast(ped$id, ped$dadid, ped$momid)
  rel <- rel_new()
  rel_add(rel, ped[, c("id", "momid", "dadid")])
  expect_equal(rel_get(rel, ped$id), full)
  cap <- rel_capacity(rel)
  ## empty it completely and start again: no value must be left behind
  rel_remove(rel, ped$id)
  expect_equal(rel_size(rel), 0L)
  expect_equal(dim(rel_get(rel)), c(0L, 0L))
  ## (under new identifiers: a removed individual cannot be added again)
  ped2 <- ped[, c("id", "momid", "dadid")]
  ped2[] <- lapply(ped2, function(x) ifelse(is.na(x), NA, paste0("new_", x)))
  rel_add(rel, ped2)
  expect_equal(unname(rel_get(rel, ped2$id)), unname(full))
  expect_equal(rel_capacity(rel), cap)         # the slots were reused
})

test_that("a removed individual cannot be added again", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  rel_add(rel, data.frame(id = "c", mom = "a", dad = "b"))
  rel_remove(rel, "a")
  ## it would come back unrelated to its offspring c
  expect_error(rel_add(rel, data.frame(id = "a", mom = NA, dad = NA)), "removed")
  expect_error(rel_add(rel, data.frame(id = c("d", "a"), mom = NA, dad = NA)), "removed")
  expect_identical(rel_ids(rel), c("b", "c"))  # nothing was inserted
  rel_add(rel, data.frame(id = "d", mom = "c", dad = "b"))
  expect_identical(rel_ids(rel), c("b", "c", "d"))
  ## a failed rel_add() does not count as a removal: x can be added later
  expect_error(rel_add(rel, data.frame(id = c("x", "p", "q"), mom = c(NA, "q", "p"),
                                       dad = NA)), "cycle")
  rel_add(rel, data.frame(id = "x", mom = NA, dad = NA))
  expect_true("x" %in% rel_ids(rel))
  ## the memory survives a save and load
  rel2 <- rel_load(rel_save(rel))
  expect_error(rel_add(rel2, data.frame(id = "a", mom = NA, dad = NA)), "removed")
  rel_add(rel2, data.frame(id = "e", mom = "d", dad = NA))
  expect_equal(rel_size(rel2), 5L)
  ## a file saved without the list of removed individuals still loads
  old <- rel_save(rel); old$removed <- NULL
  expect_equal(rel_get(rel_load(old)), rel_get(rel))
})

test_that("a metapopulation in one store: migrants keep their relatives", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a1", "a2", "b1", "b2"), mom = NA, dad = NA))
  patch <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
  rel_add(rel, data.frame(id = "a3", mom = "a1", dad = "a2"))
  patch["a3"] <- "B"                           # a3 was born in A and moves to B
  rel_add(rel, data.frame(id = "b3", mom = "b1", dad = "a3"))
  patch["b3"] <- "B"
  rel_add(rel, data.frame(id = "a4", mom = "a1", dad = "a2"))   # a sib born later in A
  patch["a4"] <- "A"
  B <- rel_get(rel, names(patch)[patch == "B"])
  expect_identical(rownames(B), c("b1", "b2", "a3", "b3"))
  expect_equal(B["b3", "a3"], 0.5)
  A <- rel_get(rel)
  expect_equal(A["a3", "a4"], 0.5)             # full sibs in different patches
  expect_equal(A["b3", "a1"], 0.25)            # grandparent left behind
  expect_equal(A["b3", "a4"], 0.25)            # uncle or aunt born after the move
})

test_that("the capacity is handled automatically and stays close to the peak", {
  expect_error(rel_new(100), "unused argument")
  rel <- rel_new()
  set.seed(1)
  ped <- sim_ped(900, n_founders = 20, window = 60)
  full <- relatedness_fast(ped$id, ped$dadid, ped$momid)
  ## births one at a time; from the 301st on, the oldest individual is removed
  for (i in seq_len(nrow(ped))) {
    rel_add(rel, ped[i, c("id", "momid", "dadid")])
    if (i > 300) rel_remove(rel, ped$id[i - 300])
    expect_gte(rel_capacity(rel), rel_size(rel))
  }
  expect_equal(rel_size(rel), 300L)
  expect_lte(rel_capacity(rel), 1.5 * 301 + 1)
  keep <- ped$id[601:900]
  expect_equal(rel_get(rel, keep), full[keep, keep])
  ## a large batch makes room in one step, without overshooting
  rel2 <- rel_new()
  rel_add(rel2, ped[, c("id", "momid", "dadid")])
  expect_equal(rel_capacity(rel2), 900L)
  expect_equal(rel_get(rel2, ped$id), full)
})

test_that("adding nothing", {
  rel <- rel_new()
  expect_identical(rel_add(rel, data.frame(id = character(0), mom = character(0),
                                           dad = character(0))), rel)
  expect_equal(rel_size(rel), 0L)
})

test_that("rel_ids, and rel_get with a subset, another order or duplicates", {
  set.seed(1)
  ped <- sim_ped(80, n_founders = 10, window = 20)
  full <- relatedness_fast(ped$id, ped$dadid, ped$momid)
  rel <- rel_new()
  for (i in seq_len(nrow(ped))) rel_add(rel, ped[i, c("id", "momid", "dadid")])
  expect_identical(rel_ids(rel), ped$id)       # order of insertion
  expect_error(rel_ids(list()), "store")
  sel <- ped$id[c(70, 3, 41, 3, 80)]
  expect_equal(rel_get(rel, sel), full[sel, sel])
  expect_equal(rel_get(rel, rev(ped$id)), full[80:1, 80:1])
  expect_equal(dim(rel_get(rel, character(0))), c(0L, 0L))
  rel_remove(rel, ped$id[c(80, 80, 79)])       # duplicates are removed once
  expect_identical(rel_ids(rel), ped$id[1:78])
})

test_that("selfing, a self-parent, and factor or numeric ids", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = 1, mom = NA, dad = NA))
  rel_add(rel, data.frame(id = factor("s"), mom = factor("1"), dad = factor("1")))
  rel_add(rel, data.frame(id = 3, mom = "s", dad = 1))
  A <- rel_get(rel, c(1, "s", 3))
  expect_equal(unname(A), unname(relatedness_fast(c(1, 2, 3), dadid = c(NA, 1, 1),
                                                   momid = c(NA, 1, 2))))
  expect_equal(A["s", "s"], 1.5)               # selfing: F = 1/2
  expect_equal(A["s", "1"], 1)
  expect_error(rel_add(rel, data.frame(id = "z", mom = "z", dad = NA)), "cycle")
  expect_identical(rel_ids(rel), c("1", "s", "3"))
})

test_that("a store read back with readRDS is rejected with a clear message", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = "a", mom = NA, dad = NA))
  f <- tempfile(fileext = ".rds")
  saveRDS(rel, f)
  dead <- readRDS(f)
  expect_error(rel_size(dead), "rel_save")
  expect_error(rel_add(dead, data.frame(id = "b", mom = NA, dad = NA)), "rel_save")
  expect_error(print(dead), "rel_save")
  expect_equal(rel_size(rel), 1L)              # the original still works
})

test_that("the C++ object rejects invalid slots and stays unchanged", {
  rel <- rel_new()
  rel_add(rel, data.frame(id = c("a", "b"), mom = NA, dad = NA))
  expect_error(rel$ks$get(NA_integer_), "slot")
  expect_error(rel$ks$get(3L), "slot")
  ## the second newborn has an invalid parent: nothing must be inserted
  expect_error(rel$ks$add(c(1L, 1L), c(2L, 7L)), "slot")
  expect_equal(rel$ks$size(), 2L)
  expect_equal(rel_ids(rel), c("a", "b"))
})
