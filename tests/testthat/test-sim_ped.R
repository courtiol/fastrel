test_that("sim_ped gives a valid pedigree in birth order", {
  ped <- sim_ped(300)
  expect_named(ped, c("id", "momid", "dadid", "sex"))
  expect_equal(nrow(ped), 300)
  ## parents always come before their offspring
  expect_true(all(match(ped$momid, ped$id) < seq_len(300), na.rm = TRUE))
  expect_true(all(match(ped$dadid, ped$id) < seq_len(300), na.rm = TRUE))
  ## ids sort the same way alphabetically and numerically
  big <- sim_ped(12000, n_founders = 40, window = 50)
  expect_identical(big$id, sort(big$id))
  expect_error(sim_ped(10, n_founders = 20), "larger")
  expect_error(sim_ped(10, n_founders = 0), "n_founders")
  expect_error(sim_ped(10, n_founders = 2, window = 0), "window")
})
