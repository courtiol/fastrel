#include <Rcpp.h>
#include "pedigree_order.h"
#include <vector>
#include <queue>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
using namespace Rcpp;

// Inverse of the additive relationship matrix A = 2 * kinship (the precision
// matrix of the additive genetic values), built directly from the pedigree
// with the rules of Henderson (1976) and Quaas (1976), without ever forming A.
//
// A = L D L', where L is lower triangular (in an order where parents come
// first) and D is diagonal: D[i] is the variance of the Mendelian sampling
// term of i, relative to the additive genetic variance:
//   both parents known:  D[i] = 1/2 - (F_m + F_d) / 4
//   one parent p known:  D[i] = 3/4 - F_p / 4
//   no parent known:     D[i] = 1
// Hence A^-1 = sum over i of q_i q_i' / D[i], with q_i = e_i - e_m / 2 - e_d / 2
// (an unknown parent dropping out): each individual adds at most 6 distinct
// entries to the upper triangle, so A^-1 is sparse whatever the pedigree.
//
// The inbreeding coefficients F are computed with the algorithm of Meuwissen
// and Luo (1992): A[i,i] = sum over the ancestors j of i (i included) of
// L[i,j]^2 D[j], where L[i,j] is obtained by going up the pedigree from i,
// each individual passing half of its coefficient to each of its parents.
// Ancestors are visited from the youngest to the oldest (decreasing position
// in the processing order), so that the coefficient of each is complete when
// it is visited. The memory used is proportional to the number of individuals.
//
// Arguments: mom, dad: 1-based row indices of each individual's parents,
// NA = unknown. Returns the nonzero entries of the upper triangle of A^-1 as
// triplets (i, j, x), 1-based in the order of the input, possibly with
// repeated (i, j) pairs whose values are to be summed, and the inbreeding
// coefficients f (in the order of the input).
// [[Rcpp::export]]
List precision_cpp(IntegerVector mom, IntegerVector dad) {
  const int n = (int)mom.size();
  if (dad.size() != n) stop("mom and dad must have the same length");

  std::vector<int> m(n), d(n);
  for (int i = 0; i < n; i++) {
    m[i] = (mom[i] == NA_INTEGER) ? -1 : mom[i] - 1;
    d[i] = (dad[i] == NA_INTEGER) ? -1 : dad[i] - 1;
    if (m[i] >= n || d[i] >= n || m[i] < -1 || d[i] < -1)
      stop("parent index out of range");
  }
  std::vector<int> pos;
  const bool sorted = parents_first(m, d, pos);
  // orig[k]: index in the input of the k-th individual processed
  std::vector<int> orig(n);
  if (sorted) for (int k = 0; k < n; k++) orig[k] = k;
  else        for (int i = 0; i < n; i++) orig[pos[i]] = i;

  // Step 1: inbreeding coefficients F and Mendelian sampling variances D,
  // in processing order (parents first).
  std::vector<double> F(n, 0.0), D(n, 1.0);
  std::vector<double> L(n, 0.0);         // coefficients L[i, j], reset after use
  std::vector<char> queued(n, 0);
  std::priority_queue<int> todo;         // ancestors to visit, youngest first
  // F only depends on the pair of parents: full siblings share it
  std::unordered_map<std::int64_t, double> F_of_parents;

  for (int i = 0; i < n; i++) {
    const int a = m[i], b = d[i];
    if (a >= 0 && b >= 0) D[i] = 0.5 - 0.25 * (F[a] + F[b]);
    else if (a >= 0 || b >= 0) D[i] = 0.75 - 0.25 * F[a >= 0 ? a : b];
    if (a < 0 || b < 0) continue;        // F = 0 without two known parents

    const std::int64_t key = (std::int64_t)std::min(a, b) * n + std::max(a, b);
    auto found = F_of_parents.find(key);
    if (found != F_of_parents.end()) { F[i] = found->second; continue; }

    // A[i,i] = sum of L[i,j]^2 D[j] over the ancestors j of i
    double Aii = 0.0;
    L[i] = 1.0; todo.push(i); queued[i] = 1;
    while (!todo.empty()) {
      const int j = todo.top(); todo.pop();
      queued[j] = 0;
      const double l = L[j];
      L[j] = 0.0;
      Aii += l * l * D[j];
      const int pj[2] = {m[j], d[j]};
      for (int p : pj) {
        if (p < 0) continue;
        L[p] += 0.5 * l;                 // twice for a selfed individual
        if (!queued[p]) { queued[p] = 1; todo.push(p); }
      }
    }
    F[i] = Aii - 1.0;
    F_of_parents.emplace(key, F[i]);
  }

  // Step 2: A^-1 = sum of q_i q_i' / D[i], upper triangle in input labels.
  std::vector<int> ri, ci;
  std::vector<double> x;
  ri.reserve((size_t)6 * n); ci.reserve((size_t)6 * n); x.reserve((size_t)6 * n);
  int who[3]; double coef[3];
  for (int i = 0; i < n; i++) {
    const double w = 1.0 / D[i];
    int k = 0;
    who[k] = i; coef[k] = 1.0; k++;
    if (m[i] >= 0) { who[k] = m[i]; coef[k] = -0.5; k++; }
    if (d[i] >= 0) {
      if (k == 2 && who[1] == d[i]) coef[1] = -1.0;   // selfing: same parent twice
      else { who[k] = d[i]; coef[k] = -0.5; k++; }
    }
    for (int u = 0; u < k; u++)
      for (int v = u; v < k; v++) {
        int r = orig[who[u]] + 1, c = orig[who[v]] + 1;
        if (r > c) std::swap(r, c);
        ri.push_back(r); ci.push_back(c);
        // an off-diagonal cell is reached once per (u, v) pair, and stands
        // for both (u, v) and (v, u) of the outer product: same value
        x.push_back(w * coef[u] * coef[v]);
      }
  }

  NumericVector f(n);
  for (int k = 0; k < n; k++) f[orig[k]] = F[k];
  return List::create(_["i"] = wrap(ri), _["j"] = wrap(ci), _["x"] = wrap(x),
                      _["f"] = f);
}
