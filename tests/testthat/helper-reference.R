## Slow but obviously correct reference: tabular method in plain R.
## Processes individuals in an order where parents come first.
ref_relatedness <- function(id, dadid, momid) {
  n <- length(id)
  m <- match(momid, id); d <- match(dadid, id)
  A <- matrix(0, n, n)
  done <- logical(n)
  while (!all(done)) {
    for (j in which(!done)) {
      pm <- m[j]; pd <- d[j]
      if ((is.na(pm) || done[pm]) && (is.na(pd) || done[pd])) {
        row <- numeric(n)
        if (!is.na(pm)) row <- row + A[pm, ]
        if (!is.na(pd)) row <- row + A[pd, ]
        row <- row / 2
        row[j] <- 1 + (if (!is.na(pm) && !is.na(pd)) A[pm, pd] else 0) / 2
        A[j, ] <- row; A[, j] <- row
        done[j] <- TRUE
      }
    }
  }
  dimnames(A) <- list(id, id)
  A
}

## The pedigree used in the documentation: H and I are offspring of a mother-son mating
test_ped <- function() {
  data.frame(
    id  = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    mom = c(NA, NA, NA, NA, "A", "C", "E", "E", "E"),
    dad = c(NA, NA, NA, NA, "B", "D", "F", "G", "G"),
    stringsAsFactors = FALSE)
}
