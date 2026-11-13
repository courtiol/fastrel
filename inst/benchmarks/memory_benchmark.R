## Benchmark of the memory used by fastrel, to check the claims made in the
## documentation (run interactively, or with Rscript).
## Needs the suggested package lobstr.
##
## Two kinds of measurements are used:
## - lobstr::obj_size() gives the exact size of an R object (a matrix returned
##   by fastrel, a saved store, ...);
## - the matrix of a store lives in memory managed by C++, which R (hence
##   lobstr) does not see. For this, we read the memory used by the whole R
##   process as reported by the operating system (the "resident set size"),
##   before and after. This is approximate (the system counts in pages, and R
##   itself allocates and frees memory in between), so sizes are chosen large
##   enough for the signal to dominate, and the comparisons are made with a
##   tolerance. It is not available on Windows, where those lines show NA.
library(fastrel)
library(lobstr)

MB <- 1e6        # as lobstr does

## memory used by the R process, in bytes
process_memory <- function() {
  if (file.exists("/proc/self/status")) {                      # Linux
    x <- grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)
    return(1024 * as.numeric(gsub("[^0-9]", "", x)))
  }
  if (.Platform$OS.type == "unix")                             # macOS, BSD
    return(1024 * as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)))
  NA_real_
}

## memory taken by the process while evaluating `expr` (and kept afterwards)
memory_taken <- function(expr) {
  invisible(gc())
  before <- process_memory()
  force(expr)
  invisible(gc())
  process_memory() - before
}

results <- list()
record <- function(claim, expected, measured, tolerance = 0.10, upper_only = FALSE) {
  ok <- if (is.na(measured)) NA
        else if (upper_only) measured <= expected * (1 + tolerance)
        else abs(measured - expected) <= tolerance * expected
  results[[length(results) + 1L]] <<- data.frame(
    claim = claim, expected = expected, measured = measured, ok = ok)
}

## a store holding n unrelated individuals
store_of <- function(n, remember_removed = FALSE) {
  rel <- rel_new(remember_removed)
  rel_add(rel, data.frame(id = seq_len(n), mom = NA, dad = NA))
  rel
}

## 1. The matrix returned for a whole pedigree uses 8 * n^2 bytes -------------
n <- 3000
set.seed(1)
ped <- sim_ped(n)
A <- fastrel_ped(ped)
## (plus the row and column names: about 64 bytes per individual)
record("fastrel_ped(): result is about 8 * n^2 bytes (MB)",
       8 * n^2 / MB, as.numeric(obj_size(A)) / MB, tolerance = 0.01)

## 2. Computing it needs no second matrix, whatever the order of the rows -----
## (peak memory used by R during the call, compared with the size of the result)
peak_during <- function(expr) {
  invisible(gc(reset = TRUE))
  base <- sum(gc()[, 2])                       # in use before (gc() reports MiB)
  force(expr)
  (sum(gc()[, 6]) - base) * 1024^2 / MB        # max used minus the baseline
}
record("fastrel_ped(), rows in birth order: peak / result",
       1, peak_during(fastrel_ped(ped)) / (8 * n^2 / MB), tolerance = 0.05)
shuffled <- ped[sample(n), ]
record("fastrel_ped(), rows shuffled: peak / result",
       1, peak_during(fastrel_ped(shuffled)) / (8 * n^2 / MB), tolerance = 0.05)
rm(A)

## 3. A store uses 8 * capacity^2 bytes --------------------------------------
n <- 4000
taken <- memory_taken(rel <- store_of(n))
record("store: memory is 8 * rel_capacity(rel)^2 bytes (MB)",
       8 * rel_capacity(rel)^2 / MB, taken / MB)
## ... which R does not see: the R side of the store is small
cat("R side of this store, as seen by lobstr::obj_size():",
    format(obj_size(rel$ids, rel$slots), units = "auto"), "\n")

## 4. rel_get() returns 8 * k^2 bytes, and a saved store 8 * n^2 bytes --------
k <- 1000
record("rel_get() of k individuals: about 8 * k^2 bytes (MB)",
       8 * k^2 / MB, as.numeric(obj_size(rel_get(rel, seq_len(k)))) / MB,
       tolerance = 0.01)
record("rel_save(): saved object is about 8 * n^2 bytes (MB)",
       8 * n^2 / MB, as.numeric(obj_size(rel_save(rel))) / MB, tolerance = 0.02)
rm(rel)

## 5. The capacity stays within about 1.5 times the peak number held ----------
## individuals are added in small batches, the worst case for the growth
worst <- 0
taken <- memory_taken({
  rel <- rel_new()
  for (i in 1:300) {
    rel_add(rel, data.frame(id = paste(i, 1:10), mom = NA, dad = NA))
    if (rel_capacity(rel) > 16)
      worst <- max(worst, rel_capacity(rel) / rel_size(rel))
  }
})
record("store: capacity / largest number held (at most about 1.5)",
       1.5, worst, tolerance = 0.02, upper_only = TRUE)
## Each time the matrix is enlarged, the old one is released, but the system
## does not always get that memory back at once: the process can hold more
## than the matrix in use.
record("store grown step by step: memory of the process / (8 * capacity^2)",
       1, taken / (8 * rel_capacity(rel)^2), tolerance = 1, upper_only = TRUE)
## a single large batch makes room exactly
rel <- store_of(2500)
record("store: capacity after one batch of 2500", 2500, rel_capacity(rel), tolerance = 0)
rm(rel)

## 6. In a simulation, memory depends on the number of individuals present, ---
##    not on the length of the pedigree
n_steps <- 600
set.seed(1)
pop <- sim_pop(n_steps = n_steps, n_init = 1000, K = 2000)
arrivals <- split(pop[, c("id", "momid", "dadid")], factor(pop$arrival, 0:n_steps))
exits    <- split(pop$id, factor(pop$exit, 0:n_steps))
rel <- rel_new()
track <- data.frame(step = 0:n_steps, pedigree = NA, present = NA, capacity = NA, memory = NA)
invisible(gc())
start <- process_memory()
for (t in 0:n_steps) {
  rel_add(rel, arrivals[[t + 1]])
  track$present[t + 1] <- rel_size(rel)          # before the removals: the peak
  rel_remove(rel, exits[[t + 1]])
  track$pedigree[t + 1] <- sum(pop$arrival <= t)
  track$capacity[t + 1] <- rel_capacity(rel)
  if (t %% 50 == 0) { invisible(gc()); track$memory[t + 1] <- (process_memory() - start) / MB }
}
print(track[track$step %% 100 == 0, ], row.names = FALSE)
second_half <- track[track$step >= n_steps / 2, ]
record("simulation: pedigree at the end / at mid-run (it keeps growing)",
       NA, max(second_half$pedigree) / min(second_half$pedigree))
record("simulation: capacity at the end / at mid-run (it does not)",
       1, max(second_half$capacity) / min(second_half$capacity), tolerance = 0)
record("simulation: capacity / largest number present (at most about 1.5)",
       1.5, max(track$capacity) / max(track$present), tolerance = 0.02, upper_only = TRUE)
mem <- second_half$memory[!is.na(second_half$memory)]
record("simulation: memory at the end / at mid-run",
       1, mem[length(mem)] / mem[1], tolerance = 0.10)
rm(rel)

## 7. remember_removed = TRUE costs 11 to 23 bytes per individual removed -----
##    for the table of traces itself, and up to twice that for the process
##    (as in 5: the memory released when the table is enlarged is not always
##    returned to the system at once)
## Two million individuals are added and removed by batches, with and without
## the option. Each run is done in a fresh R process, so that the two are
## strictly comparable (within one process, memory that R has already obtained
## from the system for the first run would hide part of the second).
n_removed <- 2e6
run <- function(remember_removed) {
  code <- sprintf('
    library(fastrel)
    rel <- rel_new(%s)
    ids <- split(sprintf("ind%%08d", seq_len(%d)), rep(1:2000, each = %d))
    for (b in ids) {
      rel_add(rel, data.frame(id = b, mom = NA, dad = NA))
      rel_remove(rel, b)
    }
    invisible(gc())
    x <- if (file.exists("/proc/self/status")) {
      gsub("[^0-9]", "", grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE))
    } else if (.Platform$OS.type == "unix") {
      system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)
    } else NA
    cat(1024 * as.numeric(x))', remember_removed, n_removed, n_removed / 2000)
  script <- tempfile(fileext = ".R")
  writeLines(code, script)
  on.exit(unlink(script))
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script),
                 stdout = TRUE, env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)))
  as.numeric(out[length(out)])
}
per_removed <- (run(TRUE) - run(FALSE)) / n_removed
record("remember_removed = TRUE: bytes per individual removed, for the process (at most 2 * 23)",
       46, per_removed, tolerance = 0, upper_only = TRUE)
rel <- rel_new(remember_removed = TRUE)
rel_add(rel, data.frame(id = 1:1000, mom = NA, dad = NA))
rel_remove(rel, 1:1000)
record("rel_save() with remember_removed: bytes saved per individual removed",
       8, as.numeric(obj_size(rel_save(rel)$removed) - obj_size(raw(0))) / 1000,
       tolerance = 0.01)

## summary ---------------------------------------------------------------------
results <- do.call(rbind, results)
results$expected <- signif(results$expected, 4)
results$measured <- signif(results$measured, 4)
print(results, row.names = FALSE, right = FALSE)
if (any(!results$ok, na.rm = TRUE)) warning("some claims are not met, see the table")
