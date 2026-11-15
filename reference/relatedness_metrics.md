# Measures of relatedness

[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
[`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
and
[`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md)
can return four measures of the genetic resemblance between individuals,
chosen with their argument `metric`. All four are derived from the same
quantity, the kinship coefficient, so choosing one or the other costs a
single pass over the matrix at most.

|  |  |  |
|----|----|----|
| `metric` | Measure | Formula |
| `"additive"` (default) | additive relationship | \\A\_{ij} = 2\phi\_{ij}\\ |
| `"kinship"` | kinship coefficient | \\\phi\_{ij}\\ |
| `"wright"` | Wright's coefficient of relationship | \\A\_{ij} / \sqrt{A\_{ii} A\_{jj}}\\ |
| `"hamilton"` | Hamilton's coefficient of relatedness | \\A\_{ij} / A\_{ii}\\ |

Their values are as follows:

|  |  |  |  |  |
|----|----|----|----|----|
| `metric` | Diagonal | Range of the diagonal | Range off the diagonal | Symmetric |
| `"additive"` | \\1 + F_i\\ | 1 to 2 | 0 to 2 | yes |
| `"kinship"` | \\(1 + F_i)/2\\ | 0.5 to 1 | 0 to 1 | yes |
| `"wright"` | 1 | 1 | 0 to 1 | yes |
| `"hamilton"` | 1 | 1 | 0 to 1 | no |

Here \\\phi\_{ij}\\ is the kinship coefficient between individuals \\i\\
and \\j\\, and \\F_i\\ is the inbreeding coefficient of \\i\\. Without
inbreeding (all \\F = 0\\) the measures `"additive"`, `"wright"` and
`"hamilton"` are identical, and equal twice the kinship coefficient: 0.5
for parent and offspring or full siblings, 0.25 for half siblings, 0.125
for first cousins. They differ only with inbreeding.

## Kinship coefficient (`"kinship"`)

The kinship coefficient \\\phi\_{ij}\\ (Malécot 1948, who called it the
*coefficient de parenté*; it is also known as the coefficient of
coancestry) is the probability that an allele drawn at random from \\i\\
and an allele drawn at random from \\j\\, at the same autosomal locus,
are identical by descent, i.e. copies of a single allele of a common
ancestor. Being a probability, it lies between 0 and 1. This is what
[`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
returns (Sinnwell et al. 2014).

Two related quantities follow. The inbreeding coefficient \\F_i\\,
introduced by Wright (1922) as the correlation between the egg and the
sperm that unite to form \\i\\, is the probability that the two alleles
of \\i\\ are identical by descent: Malécot (1948) defined it so, as his
*coefficient de consanguinité*, and noted that it equals the kinship
coefficient between the parents of \\i\\. The kinship of an individual
with itself is \\\phi\_{ii} = (1 + F_i)/2\\: when two alleles are drawn
from \\i\\ with replacement, the same allele is drawn twice with
probability 1/2, and otherwise the two alleles are identical by descent
with probability \\F_i\\.

## Additive relationship (`"additive"`)

The additive relationship \\A\_{ij} = 2\phi\_{ij}\\ is what animal and
plant breeders call the numerator relationship, as it is the numerator
of Wright's coefficient of relationship (see below). Emik and Terrill
(1949) introduced this name and described the tabular method used here
to compute it. Individuals are added one by one, parents before
offspring. When individual \\i\\, with mother \\m\\ and father \\d\\, is
added:

- its additive relationship with each individual \\j\\ already present
  is the average of those of its two parents with \\j\\, \\A\_{ij} =
  (A\_{mj} + A\_{dj}) / 2\\ (an unknown parent counting as 0);

- its additive relationship with itself is \\A\_{ii} = 1 + A\_{md} /
  2\\, where \\A\_{md}\\ is the additive relationship between its two
  parents, since its inbreeding coefficient is half of that; if a parent
  is unknown, \\A\_{ii} = 1\\.

Henderson (1976) gives the same method for the whole numerator
relationship matrix, as used to predict breeding values.

The additive relationship is not a probability but an expected count:
\\A\_{ij}\\ is the expected number of the two alleles of \\j\\ that are
identical by descent with an allele drawn at random from \\i\\, which is
why it can exceed 1 (see the section Range of values). For example, two
offspring obtained by selfing an individual that was itself obtained by
selfing have an additive relationship of 1.5 (and a kinship of 0.75).

Its main use is quantitative genetics: under an additive model, the
covariance between the breeding values of \\i\\ and \\j\\ is \\A\_{ij}
\sigma^2_A\\, where \\\sigma^2_A\\ is the additive genetic variance of
the base population: the variance matrix of breeding values is \\A
\sigma^2_A\\ (Henderson 1976; Lynch and Walsh 1998, chapter 7). In a
population without inbreeding, the covariance is thus half the additive
variance between parent and offspring and a quarter of it between half
siblings, as Falconer (1960, chapter 9) derives directly. The diagonal,
\\1 + F_i\\, shows that the breeding values of inbred individuals vary
more, since inbred individuals are more homozygous. This matrix is the
one used in animal models and to estimate heritability.

## Wright's coefficient of relationship (`"wright"`)

Wright (1922) defined the coefficient of relationship as the correlation
between two individuals, i.e. between their genotypes (Michod and
Hamilton 1980), which is obtained by dividing the covariance \\A\_{ij}\\
by the product of the two standard deviations, as in Wright's own
formula: \$\$r\_{ij} = \frac{A\_{ij}}{\sqrt{A\_{ii} A\_{jj}}} =
\frac{2\phi\_{ij}}{\sqrt{(1 + F_i)(1 + F_j)}}\$\$ Conversely, \\A\_{ij}
= r\_{ij} \sqrt{A\_{ii} A\_{jj}}\\ (Henderson 1976). This is what
[`stats::cov2cor()`](https://rdrr.io/r/stats/cor.html) computes from the
additive relationship matrix. Being a correlation, it lies between 0 and
1 (the Cauchy-Schwarz inequality guarantees that \\A\_{ij} \le
\sqrt{A\_{ii} A\_{jj}}\\), and every individual has a relationship of 1
with itself.

Its interpretation follows: it measures how similar two individuals are
compared with how similar an individual is to itself. A value of 1 means
that two individuals are **as similar as an individual is to itself**,
as with members of a line selfed for many generations; a value of 0.5
means that they are half as similar. It does not remove the effect of
inbreeding on resemblance, it measures resemblance against that
reference, which is why it stays between 0 and 1 whatever the
inbreeding. In the selfing example above, the two offspring have \\A =
1.5\\ but \\r = 1.5 / 1.75 \approx 0.86\\: they are 86 percent as
similar to each other as each of them is to itself.

## Hamilton's regression coefficient of relatedness (`"hamilton"`)

Hamilton (1964a, b) used Wright's coefficient of relationship, in a
population without inbreeding, in the first formulation of inclusive
fitness. To allow for inbreeding, Hamilton (1972) then defined the
regression coefficient of relatedness of a relative \\B\\ to an actor
\\A\\ as \\b\_{AB} = 2 r\_{AB} / (1 + F_A)\\, where \\r\_{AB}\\, the
probability that a gene of \\A\\ is present in a random gamete of \\B\\,
is their kinship coefficient. In his words, it "measures the regression
of B's additive genotype on that of A". With inbreeding, the \\r\\ of
Hamilton's rule \\rb \> c\\ is this regression coefficient rather than a
correlation. Michod and Hamilton (1980) showed that the various
coefficients of relatedness proposed for Hamilton's rule are equivalent,
and that they reduce, in particular when gene effects are additive, to
this regression of the genotype of the recipient of a behaviour on that
of the actor, \$\$b\_{ij} = \frac{\mathrm{Cov}(g_i,
g_j)}{\mathrm{Var}(g_i)} = \frac{A\_{ij}}{A\_{ii}} =
\frac{2\phi\_{ij}}{1 + F_i}\$\$ where \\i\\ is the actor and \\j\\ the
recipient. In the matrix returned, **rows are actors and columns are
recipients**: `R[i, j]` is the relatedness of recipient `j` from the
point of view of actor `i`. The matrix is not symmetric when individuals
differ in their inbreeding coefficients, which is also when the
regression and the correlation (Wright's coefficient) differ (Michod and
Hamilton 1980). Its values lie between 0 and 1 (see the section Range of
values), and all individuals have a relatedness of 1 with themselves.
Its interpretation is similar to that of Wright's coefficient, but from
the point of view of the actor only: `R[i, j] = 1` means that, as far as
the genes of the actor `i` are concerned, the recipient `j` is as
similar to `i` as `i` is to itself.

## Range of values

The ranges given in the second table above all follow from one property:
in a pedigree of diploid individuals, nobody is more related to another
individual than to itself, \$\$A\_{ij} \le \min(A\_{ii}, A\_{jj})\$\$
Consequently:

- The off-diagonal values of the additive relationship are at most 2 and
  those of the kinship at most 1, while Wright's coefficient (the
  additive relationship divided by the geometric mean of the two
  diagonal values) and Hamilton's coefficient (divided by the diagonal
  value of the actor) are at most 1.

- Without inbreeding, every diagonal value of the additive relationship
  is 1, so the additive relationship between two individuals is then at
  most 1, the value for identical twins or clones.

- Values of the additive relationship above 1 therefore only occur with
  inbreeding. For example, after `n` generations of selfing the
  inbreeding coefficient is `1 - (1/2)^n`, so the diagonal is 1.5, 1.75,
  1.875, ... The maximum of 2 (a kinship of 1) is approached but never
  reached, since founders are assumed not to be inbred.

Here is why the property \\A\_{ij} \le \min(A\_{ii}, A\_{jj})\\ holds.
Take two individuals, and call \\i\\ the one that is not an ancestor of
the other, \\j\\; let \\m\\ and \\d\\ be the parents of \\i\\, so that
\\A\_{ij} = (A\_{mj} + A\_{dj})/2\\ and \\A\_{ii} = 1 + A\_{md}/2\\.

- Draw one allele at random from each of \\j\\, \\m\\ and \\d\\. If the
  allele of \\j\\ is identical by descent both with that of \\m\\ and
  with that of \\d\\, then these two are identical by descent with each
  other. In terms of probabilities, \\\phi\_{mj} + \phi\_{dj} \le 1 +
  \phi\_{md}\\, i.e. \\A\_{mj} + A\_{dj} \le 2 + A\_{md}\\, which is
  \\A\_{ij} \le A\_{ii}\\.

- \\A\_{ij}\\ is the mean of \\A\_{mj}\\ and \\A\_{dj}\\. If neither
  exceeds \\A\_{jj}\\, their mean does not either. As \\m\\ and \\d\\
  come before \\i\\ in the pedigree, the result follows step by step
  from the founders, which are related to nobody but themselves.

An unknown parent counts as a founder that appears nowhere else in the
pedigree, so the same reasoning applies.

## Assumptions

Whatever the measure, founders (individuals with unknown parents) are
assumed to be unrelated to each other and not inbred, an unknown parent
is treated as unrelated to everybody, and only autosomal loci are
considered. All values are therefore relative to the founders of the
pedigree.

## References

Emik, L. O., & Terrill, C. E. (1949). Systematic procedures for
calculating inbreeding coefficients. *Journal of Heredity*, 40(2),
51–55.
[doi:10.1093/oxfordjournals.jhered.a105986](https://doi.org/10.1093/oxfordjournals.jhered.a105986)

Falconer, D. S. (1960). *Introduction to Quantitative Genetics*. Oliver
and Boyd, Edinburgh.

Hamilton, W. D. (1964a). The genetical evolution of social behaviour. I.
*Journal of Theoretical Biology*, 7(1), 1–16.
[doi:10.1016/0022-5193(64)90038-4](https://doi.org/10.1016/0022-5193%2864%2990038-4)

Hamilton, W. D. (1964b). The genetical evolution of social behaviour.
II. *Journal of Theoretical Biology*, 7(1), 17–52.
[doi:10.1016/0022-5193(64)90039-6](https://doi.org/10.1016/0022-5193%2864%2990039-6)

Hamilton, W. D. (1972). Altruism and related phenomena, mainly in social
insects. *Annual Review of Ecology and Systematics*, 3, 193–232.
[doi:10.1146/annurev.es.03.110172.001205](https://doi.org/10.1146/annurev.es.03.110172.001205)

Henderson, C. R. (1976). A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values.
*Biometrics*, 32(1), 69–83.
[doi:10.2307/2529339](https://doi.org/10.2307/2529339)

Lynch, M., & Walsh, B. (1998). *Genetics and Analysis of Quantitative
Traits*. Sinauer Associates, Sunderland, MA. ISBN 978-0-87893-481-2.

Malécot, G. (1948). *Les mathématiques de l'hérédité*. Masson et Cie,
Paris.

Michod, R. E., & Hamilton, W. D. (1980). Coefficients of relatedness in
sociobiology. *Nature*, 288, 694–697.
[doi:10.1038/288694a0](https://doi.org/10.1038/288694a0)

Sinnwell, J. P., Therneau, T. M., & Schaid, D. J. (2014). The kinship2 R
package for pedigree data. *Human Heredity*, 78(2), 91–93.
[doi:10.1159/000363105](https://doi.org/10.1159/000363105)

Wright, S. (1922). Coefficients of inbreeding and relationship. *The
American Naturalist*, 56(645), 330–338.
[doi:10.1086/279872](https://doi.org/10.1086/279872)

## See also

[`fastrel_ped()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md),
[`fastrel_vec()`](https://courtiol.github.io/fastrel/reference/fastrel_ped.md)
and
[`rel_get()`](https://courtiol.github.io/fastrel/reference/relstore.md),
which take the argument `metric`.

## Examples

``` r
## G is a founder; P is obtained by selfing G; O1 and O2 by selfing P
ped <- data.frame(id  = c("G", "P", "O1", "O2"),
                  dad = c(NA, "G", "P", "P"),
                  mom = c(NA, "G", "P", "P"))
fastrel_ped(ped)                       # additive relationship, the default
#>    G   P   O1   O2
#> G  1 1.0 1.00 1.00
#> P  1 1.5 1.50 1.50
#> O1 1 1.5 1.75 1.50
#> O2 1 1.5 1.50 1.75
fastrel_ped(ped, metric = "kinship")   # probabilities
#>      G    P    O1    O2
#> G  0.5 0.50 0.500 0.500
#> P  0.5 0.75 0.750 0.750
#> O1 0.5 0.75 0.875 0.750
#> O2 0.5 0.75 0.750 0.875
fastrel_ped(ped, metric = "wright")    # correlations
#>            G         P        O1        O2
#> G  1.0000000 0.8164966 0.7559289 0.7559289
#> P  0.8164966 1.0000000 0.9258201 0.9258201
#> O1 0.7559289 0.9258201 1.0000000 0.8571429
#> O2 0.7559289 0.9258201 0.8571429 1.0000000
fastrel_ped(ped, metric = "hamilton")  # regressions: row = actor
#>            G         P        O1        O2
#> G  1.0000000 1.0000000 1.0000000 1.0000000
#> P  0.6666667 1.0000000 1.0000000 1.0000000
#> O1 0.5714286 0.8571429 1.0000000 0.8571429
#> O2 0.5714286 0.8571429 0.8571429 1.0000000

## the conversions from the additive relationship matrix
A <- fastrel_ped(ped)
all.equal(fastrel_ped(ped, metric = "kinship"), A / 2)
#> [1] TRUE
all.equal(fastrel_ped(ped, metric = "wright"), cov2cor(A))
#> [1] TRUE
all.equal(fastrel_ped(ped, metric = "hamilton"), A / diag(A))
#> [1] TRUE

## without inbreeding, all but the kinship are the same
ped2 <- data.frame(id = c("a", "b", "c", "d"), dad = c(NA, NA, "b", "b"),
                   mom = c(NA, NA, "a", "a"))
all.equal(fastrel_ped(ped2), fastrel_ped(ped2, metric = "wright"))
#> [1] TRUE
all.equal(fastrel_ped(ped2), fastrel_ped(ped2, metric = "hamilton"))
#> [1] TRUE
```
