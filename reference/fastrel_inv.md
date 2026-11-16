# Inverse of the relatedness matrix (precision matrix)

Computes the inverse \\A^{-1}\\ of the additive relationship matrix \\A
= 2\phi\\ of a pedigree, i.e. the precision matrix of the additive
genetic values, as a sparse matrix. It is built directly from the
pedigree with the rules of Henderson (1976) and Quaas (1976), without
ever computing \\A\\, so it can be obtained for pedigrees far too large
for
[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md).

## Usage

``` r
fastrel_inv(ped, matrix_class = c("dgCMatrix", "dsCMatrix"))
```

## Arguments

- ped:

  A data frame (or a tibble) with one row per individual and three
  columns giving the identifier of the individual, of its father and of
  its mother (`NA` or 0 if unknown). See the section Columns of `ped`.

- matrix_class:

  The class of sparse matrix returned, from the package Matrix:

  - `"dgCMatrix"` (the default) stores both triangles of the matrix, as
    [`MCMCglmm::inverseA()`](https://rdrr.io/pkg/MCMCglmm/man/inverseA.html)
    does. This is the safe choice: it is read correctly by any software
    that accepts sparse matrices.

  - `"dsCMatrix"` stores only the upper triangle, and the matrix is
    known to be symmetric. Its entries take about 40% less memory (at
    most 6 per individual instead of 11; the diagonal is stored in both
    cases, as are the row and column names), which matters only for very
    large pedigrees. Use it only with software that reads symmetric
    sparse matrices correctly, such as
    [`spaMM::fitme()`](https://rdrr.io/pkg/spaMM/man/fitme.html):
    [`MCMCglmm::MCMCglmm()`](https://rdrr.io/pkg/MCMCglmm/man/MCMCglmm.html)
    (version 2.36) silently gives wrong results with it.

## Value

A sparse matrix of class `matrix_class`, with the identifiers of the
individuals as row and column names, in the order of the rows of `ped`.
Use [`as.matrix()`](https://rdrr.io/r/base/matrix.html) to get a base
(dense) matrix.

## Details

Whereas \\A\\ is dense in most pedigrees (nearly all individuals of a
population end up sharing an ancestor), its inverse is always sparse:
each individual adds at most six distinct entries to its upper triangle,
in the cells of itself and of its parents. This is what makes the
inverse useful in mixed models of genetic values (animal models), which
can be fitted with \\A^{-1}\\ rather than with \\A\\: e.g. as
`covStruct=list(precision=Ainv)` in
[`spaMM::fitme()`](https://rdrr.io/pkg/spaMM/man/fitme.html) or
`ginverse=list(var=Ainv)` in
[`MCMCglmm::MCMCglmm()`](https://rdrr.io/pkg/MCMCglmm/man/MCMCglmm.html).

The inverse is \\A^{-1} = \sum_i q_i q_i^\top / D_i\\, where \\q_i\\ has
1 for individual \\i\\ and \\-1/2\\ for each of its known parents, and
\\D_i\\ is the variance of the Mendelian sampling of \\i\\, relative to
the additive genetic variance: \\1/2 - (F_m + F_d)/4\\ if both parents
are known, \\3/4 - F_p/4\\ if only parent \\p\\ is known, and 1 if no
parent is known, \\F\\ being the inbreeding coefficients. These are
computed with the algorithm of Meuwissen and Luo (1992), which needs
memory proportional to the number of individuals only.

The inverse of the kinship matrix is `2 * fastrel_inv(ped)`. The other
measures of
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
are not offered here.

## References

Henderson, C. R. (1976). A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values.
*Biometrics*, 32(1), 69–83.
[doi:10.2307/2529339](https://doi.org/10.2307/2529339)

Meuwissen, T. H. E., & Luo, Z. (1992). Computing inbreeding coefficients
in large populations. *Genetics Selection Evolution*, 24(4), 305–313.
[doi:10.1186/1297-9686-24-4-305](https://doi.org/10.1186/1297-9686-24-4-305)

Quaas, R. L. (1976). Computing the diagonal elements and inverse of a
large numerator relationship matrix. *Biometrics*, 32(4), 949–953.
[doi:10.2307/2529279](https://doi.org/10.2307/2529279)

## See also

[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
for the relatedness matrix itself.

## Examples

``` r
ped <- data.frame(id  = c("A", "B", "C", "D", "E"),
                  dad = c(NA,  NA,  "A", "A", "C"),
                  mom = c(NA,  NA,  "B", "B", "D"))   # E: offspring of siblings
Ainv <- fastrel_inv(ped)
Ainv
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>    A  B    C    D  E
#> A  2  1 -1.0 -1.0  .
#> B  1  2 -1.0 -1.0  .
#> C -1 -1  2.5  0.5 -1
#> D -1 -1  0.5  2.5 -1
#> E  .  . -1.0 -1.0  2

## it is the inverse of the relatedness matrix
all.equal(as.matrix(Ainv), solve(fastrel_ped(ped)))
#> [1] TRUE

## a large pedigree: the inverse is sparse
set.seed(1)
big <- sim_ped(2000)
Ainv <- fastrel_inv(big)
Matrix::nnzero(Ainv) / prod(dim(Ainv))   # share of nonzero entries
#> [1] 0.003229

## less memory, for software that reads symmetric sparse matrices
Ainv_sym <- fastrel_inv(big, matrix_class = "dsCMatrix")
c(object.size(Ainv), object.size(Ainv_sym))
#> [1] 420592 355320
```
