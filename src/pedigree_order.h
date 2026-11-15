// Ordering of a pedigree so that parents always come before their children,
// shared by relatedness_cpp() and precision_cpp().

#ifndef FASTREL_PEDIGREE_ORDER_H
#define FASTREL_PEDIGREE_ORDER_H

#include <Rcpp.h>
#include <vector>

// m[i], d[i]: 0-based index of the mother / father of individual i, or -1 if
// unknown. If some parent comes after one of its children, the individuals
// are put in an order where parents come first (topological sort, Kahn's
// algorithm), and m and d are rewritten in that order:
//   pos[i]: position at which individual i is processed.
// Returns false in that case, and true (leaving m, d and pos unchanged) if
// the pedigree was already in such an order. Stops with an error on a cycle.
// This replaces e.g. kinship2::kindepth(): we only need a valid processing
// order, not the generation numbers.
inline bool parents_first(std::vector<int>& m, std::vector<int>& d,
                          std::vector<int>& pos) {
  const int n = (int)m.size();
  bool sorted = true;
  for (int i = 0; i < n && sorted; i++)
    if (m[i] >= i || d[i] >= i) sorted = false;
  if (sorted) return true;

  // indeg[i]: number of known parents of i (0, 1 or 2)
  // children of p are kid[start[p]] ... kid[start[p + 1] - 1]
  std::vector<int> indeg(n, 0), start(n + 1, 0);
  for (int i = 0; i < n; i++) {
    if (m[i] >= 0) { indeg[i]++; start[m[i] + 1]++; }
    if (d[i] >= 0) { indeg[i]++; start[d[i] + 1]++; }
  }
  for (int i = 0; i < n; i++) start[i + 1] += start[i];
  std::vector<int> kid(start[n]), fill(start.begin(), start.end() - 1);
  for (int i = 0; i < n; i++) {
    if (m[i] >= 0) kid[fill[m[i]]++] = i;
    if (d[i] >= 0) kid[fill[d[i]]++] = i;
  }

  std::vector<int> order;
  order.reserve(n);
  // start with the founders (no known parents)
  for (int i = 0; i < n; i++) if (indeg[i] == 0) order.push_back(i);
  // "order" grows while we iterate: when all parents of a child have been
  // placed (indeg drops to 0), the child is appended
  for (size_t h = 0; h < order.size(); h++) {
    const int q = order[h];
    for (int c = start[q]; c < start[q + 1]; c++)
      if (--indeg[kid[c]] == 0) order.push_back(kid[c]);
  }
  // if some individuals were never placed, they are part of a cycle
  if ((int)order.size() < n)
    Rcpp::stop("Impossible pedigree: someone is their own ancestor");

  pos.resize(n);
  for (int k = 0; k < n; k++) pos[order[k]] = k;
  std::vector<int> m2(n), d2(n);
  for (int k = 0; k < n; k++) {
    const int i = order[k];
    m2[k] = (m[i] >= 0) ? pos[m[i]] : -1;
    d2[k] = (d[i] >= 0) ? pos[d[i]] : -1;
  }
  m.swap(m2);
  d.swap(d2);
  return false;
}

#endif
