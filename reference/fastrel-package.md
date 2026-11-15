# fastrel: fast relatedness, relationship and kinship matrices from pedigrees

Tools to compute, for autosomal loci of diploid individuals, the
additive relationship ("relatedness") matrix \\A = 2\phi\\ (Emik and
Terrill 1949; Henderson 1976), where \\\phi\\ is the kinship coefficient
(Malécot 1948), as well as the kinship coefficient itself, Wright's
coefficient of relationship (Wright 1922) and Hamilton's regression
coefficient of relatedness (Hamilton 1972; Michod and Hamilton 1980).

## Details

Two interfaces compute the relatedness matrix itself:

- [`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  and
  [`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
  compute the whole matrix from a pedigree in one go, given as a data
  frame or as vectors (the latter is a fast replacement for
  `2 * kinship2::kinship()`).

- [`rel_new()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_add()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_remove()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  and
  [`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  maintain the matrix incrementally, which is computationally much
  cheaper when a pedigree keeps growing, as in an individual-based
  simulation.
  [`rel_ids()`](https://courtiol.github.io/fastrel/reference/relstore.md),
  [`rel_size()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  and
  [`rel_capacity()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  describe the content of such a store, and
  [`rel_save()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  and
  [`rel_load()`](https://courtiol.github.io/fastrel/reference/relstore_io.md)
  save it to disk and restore it.
  [`rel_copy()`](https://courtiol.github.io/fastrel/reference/relstore.md)
  makes an independent copy of it.

[`fastrel_inv()`](https://courtiol.github.io/fastrel/reference/fastrel_inv.md)
computes the inverse of the relatedness matrix (the precision matrix
used by animal models) as a sparse matrix, directly from the pedigree,
so that it can be obtained even for very large pedigrees.

[`sim_ped()`](https://courtiol.github.io/fastrel/reference/sim_ped.md)
simulates a pedigree, and
[`sim_pop()`](https://courtiol.github.io/fastrel/reference/sim_pop.md) a
population with births, deaths, emigration and immigration, for testing
and benchmarking.

The README of the package, at <https://github.com/courtiol/fastrel> and
on the website of the package, <https://courtiol.github.io/fastrel/>,
gives an overview with worked examples: relatedness of a whole pedigree,
incremental updates in a simulation with births and deaths,
metapopulations, and benchmarks against other packages.

## Measures of relatedness

By default, all functions return the additive relationship \\A\_{ij} =
2\phi\_{ij}\\, where \\\phi\_{ij}\\ is the kinship coefficient. The
kinship coefficient is a probability: that of an allele drawn at random
from \\i\\ and an allele drawn at random from \\j\\, at the same locus,
being identical by descent. The additive relationship is instead an
expected count of alleles: \\A\_{ij}\\ is the expected number of the two
alleles of \\j\\ that are identical by descent with an allele drawn at
random from \\i\\. This is why it lies between 0 and 2, and can exceed 1
(see
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)).
For example, an allele drawn from a parent is one of the two alleles of
its offspring with probability 1/2, so their additive relationship is
0.5 without inbreeding. Its diagonal is \\1 + F_i\\, where \\F_i\\ is
the inbreeding coefficient, so `diag(A) - 1` gives the inbreeding
coefficients. Founders (individuals with unknown parents) are assumed to
be unrelated and non-inbred, and an unknown parent is treated as
unrelated to everybody.

[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
[`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
and
[`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md)
can instead return the kinship coefficient (Malécot 1948), Wright's
coefficient of relationship (Wright 1922) or Hamilton's regression
coefficient of relatedness (Hamilton 1972), with their argument
`metric`. These measures, what they mean and when to use them are
described in
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md).

The matrices are computed with the tabular method (Emik and Terrill
1949), implemented in C++. See
[relatedness_metrics](https://courtiol.github.io/fastrel/reference/relatedness_metrics.md)
for the full list of references.

## References

Emik, L. O., & Terrill, C. E. (1949). Systematic procedures for
calculating inbreeding coefficients. *Journal of Heredity*, 40(2),
51–55.
[doi:10.1093/oxfordjournals.jhered.a105986](https://doi.org/10.1093/oxfordjournals.jhered.a105986)

Hamilton, W. D. (1972). Altruism and related phenomena, mainly in social
insects. *Annual Review of Ecology and Systematics*, 3, 193–232.
[doi:10.1146/annurev.es.03.110172.001205](https://doi.org/10.1146/annurev.es.03.110172.001205)

Henderson, C. R. (1976). A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values.
*Biometrics*, 32(1), 69–83.
[doi:10.2307/2529339](https://doi.org/10.2307/2529339)

Malécot, G. (1948). *Les mathématiques de l'hérédité*. Masson et Cie,
Paris.

Michod, R. E., & Hamilton, W. D. (1980). Coefficients of relatedness in
sociobiology. *Nature*, 288, 694–697.
[doi:10.1038/288694a0](https://doi.org/10.1038/288694a0)

Wright, S. (1922). Coefficients of inbreeding and relationship. *The
American Naturalist*, 56(645), 330–338.
[doi:10.1086/279872](https://doi.org/10.1086/279872)

## See also

Useful links:

- <https://courtiol.github.io/fastrel/>

- <https://github.com/courtiol/fastrel>

- Report bugs at <https://github.com/courtiol/fastrel/issues>

## Author

**Maintainer**: Alexandre Courtiol <alexandre.courtiol@gmail.com>
([ORCID](https://orcid.org/0000-0003-0637-2959))

Authors:

- Alexandre Courtiol <alexandre.courtiol@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0637-2959))
