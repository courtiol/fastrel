# Relatedness matrix from a pedigree

Compute a measure of relatedness of autosomal loci between all
individuals of a pedigree, using the tabular method (Emik and Terrill
1949; Henderson 1976) implemented in C++. By default, the measure is the
additive relationship \\A = 2\phi\\, where \\\phi\\ is the kinship
coefficient; the argument `metric` gives the kinship coefficient itself,
Wright's coefficient of relationship or Hamilton's regression
coefficient of relatedness instead (see
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)).
The two functions do the same thing and differ only in how the pedigree
is given:

- `fastrel_ped()` takes a data frame;

- `fastrel_vec()` takes three vectors, as
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
  does. It is a fast replacement for `2 * kinship2::kinship()`, or for
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
  when `metric = "kinship"`.

## Usage

``` r
fastrel_ped(
  ped,
  metric = c("additive", "kinship", "wright", "hamilton", "Wright", "Hamilton")
)

fastrel_vec(
  id,
  dadid,
  momid,
  sex,
  chrtype = "autosome",
  metric = c("additive", "kinship", "wright", "hamilton", "Wright", "Hamilton")
)
```

## Arguments

- ped:

  A data frame (or a tibble) with one row per individual and three
  columns giving the identifier of the individual, of its father and of
  its mother (`NA` or 0 if unknown). See the section Columns of `ped`.

- metric:

  The measure of relatedness to return: `"additive"` (the default),
  `"kinship"`, `"wright"` (or `"Wright"`) or `"hamilton"` (or
  `"Hamilton"`). Abbreviations are accepted (e.g. `"kin"`). See
  [relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
  for their definitions. In `fastrel_vec()` it comes after the arguments
  used for compatibility with
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html),
  so it is best given by name.

- id:

  Vector of unique individual identifiers. For this argument and the
  next two, a data frame with a single column (what `x[, "id"]` returns
  when `x` is a tibble) also works.

- dadid:

  Vector of the same length as `id` with the identifier of each
  individual's father (`NA` or 0 if unknown).

- momid:

  Vector of the same length as `id` with the identifier of each
  individual's mother (`NA` or 0 if unknown).

- sex, chrtype:

  Arguments of
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html),
  present so that calls written for that function work unchanged. `sex`
  is not used, as the sex of individuals does not matter for autosomes;
  if given, it must have one value per individual. `chrtype` must be
  `"autosome"`, the only type of chromosome supported: anything else
  raises an error.

## Value

A numeric matrix with one row and one column per individual, in the
order given, and the identifiers as row and column names. With the
default `metric = "additive"`, the matrix is symmetric, its diagonal is
`1 + F`, where `F` is the inbreeding coefficient, and its values lie
between 0 and 2, exceeding 1 only with inbreeding. For the other
measures, see
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md);
note that with `metric = "hamilton"` the matrix is not symmetric (rows
are actors, columns are recipients).

## Details

Rows of the pedigree can be in any order: parents do not need to appear
before their offspring. Parents that are `NA` or 0 (as in many programs
of animal breeding), or that are not found among the individuals, are
treated as unknown; 0 is refused if an individual is itself called 0, as
it would then be ambiguous. The algorithm uses a topological sort of the
pedigree, so an impossible pedigree (someone being their own ancestor)
raises an error.

Note that the argument order of `fastrel_vec()` (`dadid` before `momid`)
follows
[`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html).

## Columns of `ped`

`fastrel_ped()` and
[`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md)
find the three columns they need by their names, so the columns can be
in any order and `ped` can contain other columns, which are ignored.
Names are compared without regard to capitalisation, to punctuation
(`mother_id`, `Mother.ID` and `motherid` are the same) and to a leading
or trailing `id` (`momid`, `id_mom` and `mom` are the same). The names
recognised are:

- for the focal individual: `id`, `iid`, `focal`, `ind`, `indiv`,
  `individual`, `animal`, `ego`, `offspring`, `child`;

- for its father: `dad`, `father`, `sire`, `pat`, `paternal`;

- for its mother: `mom`, `mum`, `mother`, `dam`, `mat`, `maternal`.

If none of the columns has a recognised name, the first three columns
are taken to be the focal individual, its father and its mother, in this
order, as for the arguments of `fastrel_vec()` and
[`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
and as in pedigrees used in animal breeding (animal, sire, dam); a
message says so the first time this happens in an R session. If only
some of the three are recognised, or if several columns could play the
same role, an error is raised rather than guessing.

## References

Emik, L. O., & Terrill, C. E. (1949). Systematic procedures for
calculating inbreeding coefficients. *Journal of Heredity*, 40(2),
51–55.
[doi:10.1093/oxfordjournals.jhered.a105986](https://doi.org/10.1093/oxfordjournals.jhered.a105986)

Henderson, C. R. (1976). A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values.
*Biometrics*, 32(1), 69–83.
[doi:10.2307/2529339](https://doi.org/10.2307/2529339)

## See also

[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
for the measures of relatedness;
[`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md)
for the incremental version.

## Examples

``` r
ped <- data.frame(
  id  = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  dad = c(NA, NA, NA, NA, "B", "D", "F", "G", "G"),
  mom = c(NA, NA, NA, NA, "A", "C", "E", "E", "E"))
A <- fastrel_ped(ped)
A
#>       A     B     C     D    E    F    G     H     I
#> A 1.000 0.000 0.000 0.000 0.50 0.00 0.25 0.375 0.375
#> B 0.000 1.000 0.000 0.000 0.50 0.00 0.25 0.375 0.375
#> C 0.000 0.000 1.000 0.000 0.00 0.50 0.25 0.125 0.125
#> D 0.000 0.000 0.000 1.000 0.00 0.50 0.25 0.125 0.125
#> E 0.500 0.500 0.000 0.000 1.00 0.00 0.50 0.750 0.750
#> F 0.000 0.000 0.500 0.500 0.00 1.00 0.50 0.250 0.250
#> G 0.250 0.250 0.250 0.250 0.50 0.50 1.00 0.750 0.750
#> H 0.375 0.375 0.125 0.125 0.75 0.25 0.75 1.250 0.750
#> I 0.375 0.375 0.125 0.125 0.75 0.25 0.75 0.750 1.250

## inbreeding coefficients (H and I are the offspring of a mother-son mating)
diag(A) - 1
#>    A    B    C    D    E    F    G    H    I 
#> 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.25 0.25 

## other measures of relatedness (see ?relatedness_metrics)
fastrel_ped(ped, metric = "kinship")
#>        A      B      C      D     E     F     G      H      I
#> A 0.5000 0.0000 0.0000 0.0000 0.250 0.000 0.125 0.1875 0.1875
#> B 0.0000 0.5000 0.0000 0.0000 0.250 0.000 0.125 0.1875 0.1875
#> C 0.0000 0.0000 0.5000 0.0000 0.000 0.250 0.125 0.0625 0.0625
#> D 0.0000 0.0000 0.0000 0.5000 0.000 0.250 0.125 0.0625 0.0625
#> E 0.2500 0.2500 0.0000 0.0000 0.500 0.000 0.250 0.3750 0.3750
#> F 0.0000 0.0000 0.2500 0.2500 0.000 0.500 0.250 0.1250 0.1250
#> G 0.1250 0.1250 0.1250 0.1250 0.250 0.250 0.500 0.3750 0.3750
#> H 0.1875 0.1875 0.0625 0.0625 0.375 0.125 0.375 0.6250 0.3750
#> I 0.1875 0.1875 0.0625 0.0625 0.375 0.125 0.375 0.3750 0.6250
fastrel_ped(ped, metric = "wright")
#>           A         B         C         D         E         F         G
#> A 1.0000000 0.0000000 0.0000000 0.0000000 0.5000000 0.0000000 0.2500000
#> B 0.0000000 1.0000000 0.0000000 0.0000000 0.5000000 0.0000000 0.2500000
#> C 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000 0.5000000 0.2500000
#> D 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.5000000 0.2500000
#> E 0.5000000 0.5000000 0.0000000 0.0000000 1.0000000 0.0000000 0.5000000
#> F 0.0000000 0.0000000 0.5000000 0.5000000 0.0000000 1.0000000 0.5000000
#> G 0.2500000 0.2500000 0.2500000 0.2500000 0.5000000 0.5000000 1.0000000
#> H 0.3354102 0.3354102 0.1118034 0.1118034 0.6708204 0.2236068 0.6708204
#> I 0.3354102 0.3354102 0.1118034 0.1118034 0.6708204 0.2236068 0.6708204
#>           H         I
#> A 0.3354102 0.3354102
#> B 0.3354102 0.3354102
#> C 0.1118034 0.1118034
#> D 0.1118034 0.1118034
#> E 0.6708204 0.6708204
#> F 0.2236068 0.2236068
#> G 0.6708204 0.6708204
#> H 1.0000000 0.6000000
#> I 0.6000000 1.0000000
fastrel_ped(ped, metric = "hamilton")
#>      A    B    C    D   E   F    G     H     I
#> A 1.00 0.00 0.00 0.00 0.5 0.0 0.25 0.375 0.375
#> B 0.00 1.00 0.00 0.00 0.5 0.0 0.25 0.375 0.375
#> C 0.00 0.00 1.00 0.00 0.0 0.5 0.25 0.125 0.125
#> D 0.00 0.00 0.00 1.00 0.0 0.5 0.25 0.125 0.125
#> E 0.50 0.50 0.00 0.00 1.0 0.0 0.50 0.750 0.750
#> F 0.00 0.00 0.50 0.50 0.0 1.0 0.50 0.250 0.250
#> G 0.25 0.25 0.25 0.25 0.5 0.5 1.00 0.750 0.750
#> H 0.30 0.30 0.10 0.10 0.6 0.2 0.60 1.000 0.600
#> I 0.30 0.30 0.10 0.10 0.6 0.2 0.60 0.600 1.000

## the same from vectors
all.equal(A, fastrel_vec(ped$id, dadid = ped$dad, momid = ped$mom))
#> [1] TRUE

## columns are found by their names, whatever their order
ped2 <- data.frame(Sire = ped$dad, Dam = ped$mom, Animal = ped$id)
all.equal(A, fastrel_ped(ped2))
#> [1] TRUE

## the same as twice the kinship matrix of the kinship2 package
if (requireNamespace("kinship2", quietly = TRUE)) {
  all.equal(A, 2 * kinship2::kinship(ped$id, dadid = ped$dad, momid = ped$mom))
  all.equal(fastrel_ped(ped, metric = "kinship"),
            kinship2::kinship(ped$id, dadid = ped$dad, momid = ped$mom))
}
#> [1] TRUE
```
