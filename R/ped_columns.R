## Names recognised for the three columns of a pedigree, once in lower case,
## without punctuation and without a leading or trailing "id".
ped_names <- list(
  focal  = c("id", "iid", "focal", "ind", "indiv", "individual", "animal",
             "ego", "offspring", "child"),
  mother = c("mom", "mum", "mother", "dam", "mat", "maternal"),
  father = c("dad", "father", "sire", "pat", "paternal"))

## Role (1 = focal, 2 = mother, 3 = father, NA = none) of each column name.
ped_roles <- function(nm) {
  nm   <- gsub("[^a-z0-9]", "", tolower(nm))
  tail <- sub("^(.+)id$", "\\1", nm)             # "momid" -> "mom"; "id" stays
  head <- sub("^id(.+)$", "\\1", nm)             # "idmom" -> "mom"
  role <- rep(NA_integer_, length(nm))
  for (r in 1:3) {
    known <- ped_names[[r]]
    role[is.na(role) & (nm %in% known | tail %in% known | head %in% known)] <- r
  }
  role
}

## Positions of the columns focal, mother and father of a pedigree `ped`
## (see the section "Columns of ped" of ?fastrel_ped).
## The answer for the last set of column names seen is kept, since a
## simulation calls rel_add() over and over with the same columns.
ped_columns <- function(ped) {
  if (!is.data.frame(ped) || ncol(ped) < 3L)
    stop("`ped` must be a data frame with columns focal, mother and father")
  nm <- names(ped)
  if (identical(nm, ped_cache$names)) return(ped_cache$cols)

  role <- ped_roles(nm)
  if (all(is.na(role))) {
    cols <- 1:3                                   # no known name: by position
  } else {
    what <- c("the focal individual", "the mother", "the father")
    n_found <- tabulate(role, 3L)
    if (any(n_found > 1L)) {
      r <- which(n_found > 1L)[1L]
      stop("several columns of `ped` could be ", what[r], ": ",
           paste(nm[which(role == r)], collapse = ", "),
           ". Keep only one of them or rename the others.")
    }
    if (any(n_found == 0L))
      stop("no column of `ped` was recognised as ",
           paste(what[n_found == 0L], collapse = " or as "),
           " (columns: ", paste(nm, collapse = ", "),
           "). See the names recognised in ?fastrel_ped.")
    cols <- match(1:3, role)
  }
  ped_cache$names <- nm
  ped_cache$cols  <- cols
  cols
}

ped_cache <- new.env()
