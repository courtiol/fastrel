## Benchmark of the incremental store in a simulation with deaths, emigration and
## immigration (uses rel_remove). Compares, at every step, with recomputing the
## relatedness matrix from scratch on the whole pedigree to date.
library(fastrel)

n_steps <- 150
set.seed(1)
pop <- sim_pop(n_steps = n_steps, n_init = 300, K = 600)
cat("individuals ever present:", nrow(pop), "\n")

## the store only needs to hold the individuals that are present
rel <- rel_new()                       # allocates and grows by itself (see rel_capacity())
rec <- data.frame(step = 0:n_steps, present = 0L, ever = 0L,
                  store = 0, scratch = 0, capacity = 0L)

for (t in 0:n_steps) {
  arrivals <- pop[pop$arrival == t, c("id", "momid", "dadid")]
  leaving  <- pop$id[!is.na(pop$exit) & pop$exit == t]

  ## incremental store: add the arrivals, then remove the individuals that left
  ## (in this order, since an individual can leave at the step it reproduces)
  rec$store[t + 1] <- system.time({
    rel_add(rel, arrivals)
    rel_remove(rel, leaving)
  })[["elapsed"]]

  ## from scratch: the whole pedigree to date (ancestors are needed, so the dead
  ## cannot be dropped). Done every 10 steps to keep the benchmark short.
  seen <- pop[pop$arrival <= t, ]
  if (t %% 10 == 0)
    rec$scratch[t + 1] <- system.time(
      fastrel_vec(seen$id, seen$dadid, seen$momid))[["elapsed"]]

  rec$present[t + 1]  <- rel_size(rel)
  rec$ever[t + 1]     <- nrow(seen)
  rec$capacity[t + 1] <- rel_capacity(rel)
}

## correctness: the survivors agree with the from-scratch computation
A <- fastrel_vec(pop$id, pop$dadid, pop$momid)
alive <- rel_ids(rel)
stopifnot(isTRUE(all.equal(rel_get(rel), A[alive, alive])))

print(tail(rec[rec$step %% 10 == 0, ], 5))
cat("total store time (add + remove, all steps):", sum(rec$store), "s\n")

op <- par(mfrow = c(1, 2))
plot(present ~ step, data = rec, type = "l", ylim = c(0, max(rec$ever)),
     ylab = "individuals")
lines(ever ~ step, data = rec, lty = 2)
legend("topleft", lty = 1:2, bty = "n", legend = c("in the store", "in the pedigree"))
sc <- rec[rec$scratch > 0, ]
plot(store ~ step, data = rec, type = "l", ylab = "seconds per step",
     ylim = c(0, max(sc$scratch, rec$store)))
lines(scratch ~ step, data = sc, lty = 2, type = "b")
legend("topleft", lty = 1:2, bty = "n", legend = c("store (add + remove)", "from scratch"))
par(op)
