// Conversion of an additive relationship matrix A = 2 * kinship into the
// other measures of relatedness the package can return (see ?relatedness_metrics).
// The conversion is done in place, in a single pass over the matrix, so that
// no second matrix is ever allocated.

#ifndef FASTREL_METRICS_H
#define FASTREL_METRICS_H

#include <vector>
#include <cmath>

// codes shared with R (see metric_code() in R/metrics.R)
enum Metric { ADDITIVE = 0, KINSHIP = 1, WRIGHT = 2, HAMILTON = 3 };

// p: k x k matrix in column-major order (element (i,j) is p[i + j*k]),
// holding A on entry and the requested measure on exit.
inline void apply_metric(double* p, int k, int metric) {
  if (metric == ADDITIVE || k == 0) return;
  const size_t K = (size_t)k;

  if (metric == KINSHIP) {                 // phi = A / 2
    for (size_t x = 0; x < K * K; x++) p[x] *= 0.5;
    return;
  }

  // the diagonal, A[i,i] = 1 + F_i >= 1, is read before it is overwritten
  std::vector<double> s(k);
  if (metric == WRIGHT) {                  // r = A[i,j] / sqrt(A[i,i] * A[j,j])
    for (int i = 0; i < k; i++) s[i] = 1.0 / std::sqrt(p[i + i * K]);
    for (int j = 0; j < k; j++) {
      double* cj = p + j * K;
      const double sj = s[j];
      for (int i = 0; i < k; i++) cj[i] *= s[i] * sj;
      cj[j] = 1.0;                         // exactly 1, without rounding error
    }
  } else if (metric == HAMILTON) {         // b = A[i,j] / A[i,i], row i = actor
    for (int i = 0; i < k; i++) s[i] = 1.0 / p[i + i * K];
    for (int j = 0; j < k; j++) {
      double* cj = p + j * K;
      for (int i = 0; i < k; i++) cj[i] *= s[i];
      cj[j] = 1.0;
    }
  }
}

#endif
