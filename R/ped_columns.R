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
## (in this order, which is the order of the roles, not of the columns)
## (see the section "Columns of ped" of ?fastrel_ped).
## The answer for the last set of column names seen is kept, since a
## simulation calls rel_add() over and over with the same columns.
ped_columns <- function(ped) {
  if (!is.data.frame(ped) || ncol(ped) < 3L)
    stop("`ped` must be a data frame with columns for the individual, its ",
         "father and its mother")
  nm <- names(ped)
  if (identical(nm, ped_cache$names)) return(ped_cache$cols)

  role <- ped_roles(nm)
  if (all(is.na(role))) {
    ## no known name: by position, in the order individual, father, mother
    ## (that of kinship2::kinship() and of pedigrees in animal breeding)
    cols <- c(1L, 3L, 2L)
    if (!isTRUE(ped_cache$told_position)) {       # say so once per session
      ped_cache$told_position <- TRUE
      message("No column of `ped` has a recognised name, so its first three ",
              "columns (", paste(nm[1:3], collapse = ", "), ") are taken as the ",
              "individual, its father and its mother, in this order. Name them ",
              "(e.g. id, dad, mom) to avoid any ambiguity; see ?fastrel_ped. ",
              "This message is shown once per session.")
    }
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

## Memory of the session: the last column names seen (see above), and whether
## the message about columns taken by position was shown.
ped_cache <- new.env()

## Identifiers given as a vector. A data frame with a single column is
## accepted too and replaced by that column: this is what `x[, "id"]` returns
## when `x` is a tibble, where a base data frame returns the vector. Anything
## else that is not a plain vector is refused, as it would otherwise be
## treated as a single identifier.
id_vector <- function(x, arg) {
  if (is.data.frame(x)) {
    if (ncol(x) != 1L)
      stop("`", arg, "` must be a vector, not a data frame with ", ncol(x), " columns")
    x <- x[[1L]]
  }
  if (!is.atomic(x) && !is.null(x))
    stop("`", arg, "` must be a vector")
  x
}

## The parent code 0 means "unknown", as in many programs of animal breeding,
## unless an individual is itself called 0, which would make it ambiguous.
## `id`: identifiers of the individuals that 0 could refer to; `parents`:
## identifiers of parents (any type). Returns `parents`, of the same type, with
## the 0 (or "0") replaced by NA.
zero_as_unknown <- function(parents, id) {
  zero <- !is.na(parents) & as.character(parents) == "0"
  if (!any(zero)) return(parents)
  if (any(as.character(id) == "0", na.rm = TRUE))
    stop("0 is used as a parent, but it is also the identifier of an individual: ",
         "it cannot also mean an unknown parent. Use NA for unknown parents, or ",
         "give that individual another identifier.")
  parents[zero] <- NA
  parents
}
