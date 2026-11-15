#include <Rcpp.h>
#include "metrics.h"
#include "pedigree_order.h"
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
//   metric: measure to return, see metrics.h (0 = A itself)
//
// Implementation notes (speed):
//   * Individuals are processed in an order where parents come first, so when
//     j is processed only its relatedness with the individuals processed
//     before it needs computing: the others fill that cell in when their own
//     turn comes. This halves the work compared with filling whole rows.
//   * The matrix is filled in "processing order" labels (individual number k
//     is the k-th to be processed). Then the individuals processed before k
//     are simply 0..k-1, a contiguous block, and consecutive individuals write
//     to neighbouring memory. If the pedigree was not already sorted, the
//     matrix is relabelled in place at the end, which needs no second matrix.
// [[Rcpp::export]]
NumericMatrix relatedness_cpp(IntegerVector mom, IntegerVector dad, int metric = 0) {
  const int n = (int)mom.size();
  if (dad.size() != n) stop("mom and dad must have the same length");

  // m[i], d[i]: 0-based index of mother / father of i, or -1 if unknown
  std::vector<int> m(n), d(n);
  for (int i = 0; i < n; i++) {
    m[i] = (mom[i] == NA_INTEGER) ? -1 : mom[i] - 1;   // R is 1-based, C++ 0-based
    d[i] = (dad[i] == NA_INTEGER) ? -1 : dad[i] - 1;
    if (m[i] >= n || d[i] >= n || m[i] < -1 || d[i] < -1)
      stop("parent index out of range");
  }

  // Step 1: if needed, order the individuals so that parents always come
  // before their children (see pedigree_order.h). Afterwards m and d are
  // expressed in processing order, and pos[i] is the position at which
  // individual i is processed.
  std::vector<int> pos;
  const bool sorted = parents_first(m, d, pos);

  // Step 2: fill the matrix (in processing order: parents of j are < j).
  NumericMatrix A(n, n);               // zero-initialised: everyone unrelated
  double* p = A.begin();               // column-major: element (i,j) is p[i + j*n]

  for (int j = 0; j < n; j++) {
    double* colj = p + (size_t)j * n;

    // Founder: row and column are already 0, the diagonal is 1 (no inbreeding).
    if (m[j] < 0 && d[j] < 0) { colj[j] = 1.0; continue; }

    // Relatedness of j with every individual i processed before it: the
    // average of the mother's and the father's relatedness with i. Since A is
    // symmetric, a parent's column is also its row: A[m, i] = pm[i]. An
    // unknown parent contributes 0.
    if (m[j] >= 0 && d[j] >= 0) {
      const double* pm = p + (size_t)m[j] * n;
      const double* pd = p + (size_t)d[j] * n;
      for (int i = 0; i < j; i++) colj[i] = (pm[i] + pd[i]) * 0.5;
      // Diagonal: 1 + F_j, where F_j = A[m,d] / 2 is the inbreeding
      // coefficient (probability that j's two alleles are identical by descent).
      colj[j] = 1.0 + pm[d[j]] * 0.5;
    } else {
      const double* pp = p + (size_t)(m[j] >= 0 ? m[j] : d[j]) * n;
      for (int i = 0; i < j; i++) colj[i] = pp[i] * 0.5;
      colj[j] = 1.0;                   // a single known parent: not inbred
    }

    // Copy the result into row j, by symmetry. The column write above is
    // contiguous; this one is strided (the slow part).
    double* rowj = p + j;
    for (int i = 0; i < j; i++) rowj[(size_t)i * n] = colj[i];
  }

  // Step 3: back to the labels of the input, A[i, j] = B[pos[i], pos[j]],
  // done in place with a buffer of one column.
  if (!sorted) {
    std::vector<double> buf(n);
    // rows, one column at a time
    for (int j = 0; j < n; j++) {
      double* colj = p + (size_t)j * n;
      for (int i = 0; i < n; i++) buf[i] = colj[pos[i]];
      std::copy(buf.begin(), buf.end(), colj);
    }
    // columns: new column j is old column pos[j]; follow each cycle of the
    // permutation so that every column is moved exactly once
    std::vector<char> moved(n, 0);
    for (int s = 0; s < n; s++) {
      if (moved[s]) continue;
      moved[s] = 1;
      if (pos[s] == s) continue;
      std::copy(p + (size_t)s * n, p + (size_t)(s + 1) * n, buf.begin());
      int cur = s;
      while (pos[cur] != s) {
        const int src = pos[cur];
        std::copy(p + (size_t)src * n, p + (size_t)(src + 1) * n, p + (size_t)cur * n);
        moved[src] = 1;
        cur = src;
      }
      std::copy(buf.begin(), buf.end(), p + (size_t)cur * n);
    }
  }
  // Step 4: the measure requested, in place
  apply_metric(p, n, metric);
  return A;
}
