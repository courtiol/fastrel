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
