// Created by Claude Sonnet 5.5 from Alex's prompts

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Additive relationship ("relatedness") matrix A = 2 * kinship, computed with
// the tabular (Henderson) method.
//
// Definitions, with phi = kinship coefficient and F = inbreeding coefficient:
//   A[i,j] = 2 * phi(i,j)
//   A[i,i] = 1 + F_i               (so founders have 1, inbred individuals > 1)
//
// Recursion for a non-founder j with mother m and father d, for any i != j:
//   A[j,i] = ( A[m,i] + A[d,i] ) / 2      -> j's alleles come half from each parent
//   A[j,j] = 1 + A[m,d] / 2               -> 1 + F_j, where F_j = A[m,d] / 2
// A missing parent counts as unrelated to everyone (contributes 0).
// Founders are unrelated to each other: A[i,i] = 1, off-diagonal = 0.
//
// Arguments:
//   mom, dad: 1-based row indices of each individual's parents, NA = unknown
// [[Rcpp::export]]
NumericMatrix relatedness_cpp(IntegerVector mom, IntegerVector dad) {
  const int n = mom.size();

  // m[i], d[i]: 0-based index of mother / father of i, or -1 if unknown
  // indeg[i]  : number of known parents of i (0, 1 or 2)
  // kids[p]   : list of children of p (used to walk the pedigree downwards)
  std::vector<int> m(n), d(n), indeg(n, 0);
  std::vector<std::vector<int>> kids(n);

  for (int i = 0; i < n; i++) {
    m[i] = (mom[i] == NA_INTEGER) ? -1 : mom[i] - 1;   // R is 1-based, C++ 0-based
    d[i] = (dad[i] == NA_INTEGER) ? -1 : dad[i] - 1;
    if (m[i] >= 0) { indeg[i]++; kids[m[i]].push_back(i); }
    if (d[i] >= 0) { indeg[i]++; kids[d[i]].push_back(i); }
  }

  // Step 1: order the individuals so that parents always come before their
  // children (topological sort, Kahn's algorithm). This replaces e.g. kinship2::kindepth():
  // we only need a valid processing order, not the generation numbers.
  std::vector<int> order;
  order.reserve(n);
  // start with the founders (no known parents)
  for (int i = 0; i < n; i++) if (indeg[i] == 0) order.push_back(i);
  // "order" grows while we iterate: when all parents of a child have been
  // placed (indeg drops to 0), the child is appended
  for (size_t h = 0; h < order.size(); h++) {
    for (int c : kids[order[h]])
      if (--indeg[c] == 0) order.push_back(c);
  }
  // if some individuals were never placed, they are part of a cycle
  if ((int)order.size() < n)
    stop("Impossible pedigree: someone is their own ancestor");

  // Step 2: fill the matrix.
  NumericMatrix A(n, n);               // zero-initialised: everyone unrelated
  double* p = A.begin();               // column-major: element (i,j) is p[i + j*n]
  // Diagonal starts at 1 (no inbreeding). Founders keep this value; it is
  // overwritten below for non-founders, whose diagonal depends on their parents.
  for (int i = 0; i < n; i++) p[i + (size_t)i * n] = 1.0;

  std::vector<double> v(n);            // scratch buffer for the new row/column
  for (int j : order) {                // parents are always processed before j
    // Founder: row and column are already 0, diagonal already 1. Nothing to do.
    if (m[j] < 0 && d[j] < 0) continue;

    // Pointers to the columns of the mother and the father. Since A is
    // symmetric, a column is also a row: A[m, i] = pm[i].
    // nullptr means the parent is unknown (contributes 0).
    const double* pm = (m[j] >= 0) ? p + (size_t)m[j] * n : nullptr;
    const double* pd = (d[j] >= 0) ? p + (size_t)d[j] * n : nullptr;

    // Relatedness between the two parents, A[m,d]. Needed for j's diagonal.
    // If either parent is unknown, they count as unrelated: 0.
    const double amd = (m[j] >= 0 && d[j] >= 0) ? p[m[j] + (size_t)d[j] * n] : 0.0;

    // Relatedness of j with every individual i: the average of the mother's
    // and the father's relatedness with i.
    for (int i = 0; i < n; i++)
      v[i] = ((pm ? pm[i] : 0.0) + (pd ? pd[i] : 0.0)) / 2.0;

    // Diagonal: 1 + F_j, where F_j = A[m,d] / 2 is the inbreeding coefficient
    // (probability that j's two alleles are identical by descent). This
    // overrides the value that the averaging gave for i == j.
    v[j] = 1.0 + amd / 2.0;

    // Write the result as column j and, by symmetry, as row j.
    // The column write is contiguous; the row write is strided (the slow part).
    double* colj = p + (size_t)j * n;
    for (int i = 0; i < n; i++) {
      colj[i] = v[i];                  // A[i,j]
      p[j + (size_t)i * n] = v[i];     // A[j,i]
    }
  }
  return A;
}
