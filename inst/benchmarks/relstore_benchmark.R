## Benchmarks of fastrel against kinship2 (run interactively).
## Needs the suggested packages kinship2 and microbenchmark.
library(fastrel)
library(kinship2)
library(microbenchmark)

## kinship2::kinship returns kinship, so double it to get relatedness
source_orig <- function(id, dadid, momid, ...) {
  2 * kinship2::kinship(id, dadid = dadid, momid = momid, ...)
}

## comparison ----
ped <- sim_ped(1000)                                # columns: id, momid, dadid, sex
rel <- rel_new(64)                                  # small start, to exercise growth
rel_add(rel, ped[, c("id", "momid", "dadid")])      # any row order works
A1 <- rel_get(rel, ped$id)
A2 <- relatedness_fast(ped$id, ped$dadid, ped$momid)
A3 <- source_orig(ped$id, ped$dadid, ped$momid)

all.equal(A1, A2)
all.equal(A1, A3)

## speed ----
## orig and fast compute the whole matrix from the pedigree. store_build does
## the same with a store (creation, insertion and extraction), whereas
## store_get only extracts the matrix from a store that is already filled.
microbenchmark(orig  = source_orig(ped$id, ped$dadid, ped$momid),
               fast  = relatedness_fast(ped$id, ped$dadid, ped$momid),
               store_build = {
                 r <- rel_new(1000)
                 rel_add(r, ped[, c("id", "momid", "dadid")])
                 rel_get(r, ped$id)
               },
               store_get = rel_get(rel, ped$id),
               times = 200)

## scaling: a pedigree that keeps growing ----
## sim_ped() returns individuals in birth order, so the first n rows always
## form a closed pedigree (parents come before their offspring)
ped <- sim_ped(10000)

batch <- 20
steps <- seq(1000, 4000, by = batch)

rel <- rel_new(10)
rel_add(rel, ped[1:1000, c("id", "momid", "dadid")])

ped_ori <- ped[1:1000, ]
rec_time <- data.frame(ped_length = steps, relstore = 0, orig = 0, cap = 0)

for (n in steps[-1]) {
  new <- ped[(n - batch + 1):n, ]
  k   <- rec_time$ped_length == n

  ## incremental store: add the newborns and extract the whole matrix
  rec_time$relstore[k] <- system.time({
    rel_add(rel, new[, c("id", "momid", "dadid")])
    rel_get(rel)
  })[["elapsed"]]

  ## from scratch on the whole pedigree
  ped_ori <- rbind(ped_ori, new)
  rec_time$orig[k] <- system.time(
    source_orig(id = ped_ori$id, dadid = ped_ori$dadid, momid = ped_ori$momid)
  )[["elapsed"]]

  rec_time$cap[k] <- rel_capacity(rel)
}

head(rec_time)

## correctness at the end of the loop
all.equal(rel_get(rel, ped_ori$id),
          source_orig(ped_ori$id, ped_ori$dadid, ped_ori$momid))

rt <- rec_time[-1, ]   # the first row is the initial batch (not timed)
plot(relstore ~ ped_length, data = rt, type = "l", ylim = c(0, 0.4))
lines(orig ~ ped_length, data = rt, lty = 2)
lines(0.3 * (cap / max(cap)) ~ ped_length, data = rt, col = 2)
legend("topleft", lty = c(1, 2, 1), col = c(1, 1, 2), bty = "n",
       legend = c("store (add + get)", "kinship2 from scratch", "capacity (rescaled)"))
