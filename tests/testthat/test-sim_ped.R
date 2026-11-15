test_that("sim_ped gives a valid pedigree in birth order", {
  set.seed(1)
  ped <- sim_ped(300)
  expect_named(ped, c("id", "momid", "dadid", "sex"))
  expect_equal(nrow(ped), 300)
  ## parents always come before their offspring
  expect_true(all(match(ped$momid, ped$id) < seq_len(300), na.rm = TRUE))
  expect_true(all(match(ped$dadid, ped$id) < seq_len(300), na.rm = TRUE))
  ## ids sort the same way alphabetically and numerically
  set.seed(1)
  big <- sim_ped(12000, n_founders = 40, window = 50)
  expect_identical(big$id, sort(big$id))
  expect_error(sim_ped(10, n_founders = 20), "larger")
  expect_error(sim_ped(10, n_founders = 0), "n_founders")
  expect_error(sim_ped(10, n_founders = 2, window = 0), "window")
  expect_error(sim_ped(10, n_founders = 2, p_unknown = 1.5), "probability")
})

test_that("the simulators follow the random number generator of the session", {
  set.seed(1); a <- sim_ped(100)
  set.seed(1); expect_identical(sim_ped(100), a)
  expect_false(identical(sim_ped(100), a))     # the generator has moved on
  set.seed(5); p1 <- sim_pop(10, n_init = 20)
  set.seed(5); expect_identical(sim_pop(10, n_init = 20), p1)
  expect_error(sim_ped(100, seed = 1), "unused argument")
})

test_that("numeric_ids gives the same pedigree with integers and 0 for unknown parents", {
  set.seed(3); a <- sim_ped(300, p_unknown = 0.2)
  set.seed(3); b <- sim_ped(300, p_unknown = 0.2, numeric_ids = TRUE)
  expect_named(b, c("id", "momid", "dadid", "sex"))
  expect_identical(b$id, 1:300)
  expect_type(b$momid, "integer")
  expect_type(b$dadid, "integer")
  expect_identical(b$sex, a$sex)
  num <- function(x) { k <- match(x, a$id); k[is.na(k)] <- 0L; k }
  expect_identical(b$momid, num(a$momid))
  expect_identical(b$dadid, num(a$dadid))
  expect_true(any(b$momid == 0L & b$id > 40))            # unknown parents occur
  ## parents come before their offspring
  expect_true(all(b$momid < b$id & b$dadid < b$id))
  ## fastrel treats 0 as an unknown parent: same matrix as the character version
  expect_equal(unname(fastrel_ped(b)), unname(fastrel_ped(a)))
  expect_equal(unname(fastrel_vec(b$id, b$dadid, b$momid, metric = "kinship")),
               unname(fastrel_ped(a, metric = "kinship")))
  ## rel_add() takes 0 as an unknown parent too
  rel <- rel_new()
  rel_add(rel, b)
  expect_equal(unname(rel_get(rel, b$id)), unname(fastrel_ped(a)))
  expect_error(sim_ped(10, n_founders = 2, numeric_ids = NA), "TRUE or FALSE")
})

test_that("numeric_ids pedigrees can be given as they are to kinship2 and ggroups", {
  set.seed(4)
  p <- sim_ped(250, n_founders = 5, window = 15, numeric_ids = TRUE)
  A <- fastrel_ped(p)
  skip_if_not_installed("kinship2")
  expect_equal(unname(2 * kinship2::kinship(p$id, p$dadid, p$momid)), unname(A),
               tolerance = 1e-12)
  skip_if_not_installed("ggroups")
  expect_equal(unname(ggroups::buildA(p[, c("id", "dadid", "momid")])), unname(A),
               tolerance = 1e-12)
})
