test_that("sim_pop returns a consistent pedigree", {
  set.seed(5)
  pop <- sim_pop(30, n_init = 40)
  expect_named(pop, c("id", "momid", "dadid", "sex", "arrival", "exit", "fate", "origin"))
  expect_false(anyDuplicated(pop$id) > 0)
  expect_equal(sum(pop$origin == "founder"), 40)
  ## parents come first and arrived earlier than their offspring
  expect_true(all(match(pop$momid, pop$id) < seq_len(nrow(pop)), na.rm = TRUE))
  born <- pop[pop$origin == "born", ]
  expect_true(all(pop$arrival[match(born$momid, pop$id)] < born$arrival))
  expect_true(all(pop$sex[match(born$momid, pop$id)] == 2))
  expect_true(all(pop$sex[match(born$dadid, pop$id)] == 1, na.rm = TRUE))
  ## exits are after arrivals, and fate is set exactly when exit is
  left <- !is.na(pop$exit)
  expect_true(all(pop$exit[left] > pop$arrival[left]))
  expect_identical(left, !is.na(pop$fate))
  expect_true(all(pop$fate[left] %in% c("dead", "emigrated")))
  ## parents were still present when the offspring were born
  m <- match(born$momid, pop$id)
  expect_true(all(is.na(pop$exit[m]) | pop$exit[m] >= born$arrival))
})

test_that("population is roughly stable and replay with the store is exact", {
  set.seed(7)
  pop <- sim_pop(40, n_init = 50)
  full <- fastrel_vec(pop$id, pop$dadid, pop$momid)
  rel <- rel_new()
  for (t in 0:40) {
    rel_add(rel, pop[pop$arrival == t, c("id", "momid", "dadid")])
    rel_remove(rel, pop$id[!is.na(pop$exit) & pop$exit == t])
    n_alive <- sum(pop$arrival <= t & (is.na(pop$exit) | pop$exit > t))
    expect_equal(rel_size(rel), n_alive)
  }
  expect_gt(rel_size(rel), 10)
  expect_lt(rel_size(rel), nrow(pop))            # removals actually happened
  alive <- rel_ids(rel)
  expect_equal(rel_get(rel), full[alive, alive])
})

test_that("edge cases", {
  expect_error(sim_pop(5, n_init = 0), "n_init")
  expect_error(sim_pop(5, p_death = 2), "probabilities")
  expect_error(sim_pop(-1), "n_steps")
  expect_error(sim_pop(5, K = 0), "K")
  expect_error(sim_pop(5, carrying = 100), "unused argument")
  set.seed(2); small <- sim_pop(60, n_init = 40, K = 50)
  set.seed(2); large <- sim_pop(60, n_init = 40, K = 400)
  expect_lt(nrow(small), nrow(large))            # K limits the births
  expect_error(sim_pop(5, n_init = 2.5), "n_init")
  expect_error(sim_pop(5, maturity = -1), "maturity")
  set.seed(1)
  expect_warning(pop <- sim_pop(10, n_init = 3, p_death = 1, n_immigrants = 0), "extinct")
  expect_true(all(!is.na(pop$exit)))
  set.seed(1)
  none <- sim_pop(0, n_init = 5)
  expect_equal(nrow(none), 5)
  expect_true(all(is.na(none$exit)))
})
