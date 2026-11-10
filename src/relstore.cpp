// Created by Claude Sonnet 5.5 from Alex's prompts

#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Persistent store of the additive relationship ("relatedness") matrix
// A = 2 * kinship among the individuals currently in the store.
//
//   A[i,j] = 2 * phi(i,j)       (phi = kinship coefficient)
//   A[i,i] = 1 + F_i            (F_i = inbreeding coefficient)
//
// Key property: A between existing individuals never changes when new ones
// are born. A newborn's row is computed once from its parents' rows:
//   A[c,i] = ( A[m,i] + A[d,i] ) / 2     for every other individual i
//   A[c,c] = 1 + A[m,d] / 2              (F_c = A[m,d] / 2)
// so each birth costs O(capacity) instead of recomputing the whole matrix.
//
// Layout: a dense capacity x capacity matrix in column-major order, element
// (i,j) is K[i + j*cap]. Each individual occupies a "slot" (row/column index).
// Removed individuals free their slot for reuse.
// Class invariant: rows/columns of free slots are all zero.
class RelStore {
public:
  explicit RelStore(int capacity) { init(capacity); }

  // Add individuals. mom/dad are the SLOTS of the parents (1-based, NA =
  // unknown parent), not the slots of the newborns. Returns the slots that
  // were assigned to the newborns, in order.
  // Newborns are processed one after the other in the order given, so a
  // newborn cannot use another newborn of the same call as a parent.
  IntegerVector add(IntegerVector mom, IntegerVector dad) {
    const int nb = mom.size();
    if (dad.size() != nb) stop("mom and dad must have the same length");
    IntegerVector out(nb);
    std::vector<double> v;                 // scratch buffer for the new row/column

    for (int b = 0; b < nb; b++) {
      // 0-based parent slots (-1 if unknown); slot() also validates them
      const int m = slot(mom[b]), d = slot(dad[b]);

      // take a free slot for the newborn, enlarging the matrix if none is left
      if (free_.empty()) grow();
      const int j = free_.back(); free_.pop_back();
      used[j] = true;
      v.assign(cap, 0.0);                  // unrelated to everyone by default

      // Parent columns (A is symmetric, so a column is also a row).
      // nullptr = unknown parent, which contributes 0.
      const double* pm = (m >= 0) ? col(m) : nullptr;
      const double* pd = (d >= 0) ? col(d) : nullptr;

      // Relatedness between the two parents, A[m,d]; 0 if either is unknown.
      const double amd = (m >= 0 && d >= 0) ? K[m + (size_t)d * cap] : 0.0;

      // Relatedness of the newborn with every slot i: average of the mother's
      // and the father's relatedness with i. Free slots hold zeros, so they
      // stay zero. A founder (no parents) keeps v = 0 everywhere.
      if (pm || pd)
        for (int i = 0; i < cap; i++)
          v[i] = ((pm ? pm[i] : 0.0) + (pd ? pd[i] : 0.0)) / 2.0;

      // Diagonal: 1 + F_j, with F_j = A[m,d] / 2. For a founder or an
      // individual with a single known parent, amd = 0, so the diagonal is 1.
      // This overrides the value the averaging produced for i == j.
      v[j] = 1.0 + amd / 2.0;

      // Write v as column j and, by symmetry, as row j.
      // The column write is contiguous; the row write is strided (the slow part).
      double* cj = col(j);
      for (int i = 0; i < cap; i++) {
        cj[i] = v[i];                      // A[i,j]
        K[j + (size_t)i * cap] = v[i];     // A[j,i]
      }
      out[b] = j + 1;                      // back to 1-based for R
    }
    return out;
  }

  // Remove individuals (they can no longer be used as parents). Their row and
  // column are zeroed to keep the invariant, and the slot is recycled.
  // Relatedness values of the remaining individuals are unaffected.
  void remove(IntegerVector slots) {
    for (int s : slots) {
      const int j = slot(s);
      if (j < 0) continue;                 // NA: nothing to remove
      double* cj = col(j);
      for (int i = 0; i < cap; i++) {
        cj[i] = 0.0;
        K[j + (size_t)i * cap] = 0.0;
      }
      used[j] = false;
      free_.push_back(j);
    }
  }

  // Relatedness submatrix for the given slots (already on the A = 2*phi scale).
  NumericMatrix get(IntegerVector slots) {
    const int k = slots.size();
    // convert and validate the slots once, not inside the double loop
    std::vector<int> s(k);
    for (int a = 0; a < k; a++) s[a] = slot(slots[a]);

    NumericMatrix out(k, k);
    for (int a = 0; a < k; a++)
      for (int b = 0; b < k; b++)
        out(b, a) = K[s[b] + (size_t)s[a] * cap];
    return out;
  }

  // 1-based slots currently in use
  IntegerVector alive() {
    std::vector<int> s;
    for (int i = 0; i < cap; i++) if (used[i]) s.push_back(i + 1);
    return wrap(s);
  }

  int size() { return (int)(cap - free_.size()); }   // individuals in the store
  int capacity() { return cap; }                     // slots available (used + free)

private:
  int cap;                    // current number of slots
  std::vector<double> K;      // cap x cap relatedness matrix, column-major
  std::vector<char> used;     // used[i] = slot i is occupied
  std::vector<int> free_;     // stack of free slots

  void init(int c) {
    cap = c;
    K.assign((size_t)c * c, 0.0);
    used.assign(c, false);
    free_.clear();
    for (int i = c - 1; i >= 0; i--) free_.push_back(i);   // slot 0 is handed out first
  }

  double* col(int j) { return K.data() + (size_t)j * cap; }

  // 1-based R slot -> 0-based index; NA -> -1 (unknown parent);
  // an out-of-range or removed slot raises an error.
  int slot(int s) {
    if (s == NA_INTEGER) return -1;
    if (s < 1 || s > cap || !used[s - 1]) stop("invalid or removed slot: %d", s);
    return s - 1;
  }

  // Double the capacity: copy the old matrix into the top-left corner of a
  // larger, zero-filled one and register the new slots as free.
  void grow() {
    const int nc = cap * 2;
    std::vector<double> NK((size_t)nc * nc, 0.0);
    for (int j = 0; j < cap; j++)
      std::copy(K.begin() + (size_t)j * cap, K.begin() + (size_t)(j + 1) * cap,
                NK.begin() + (size_t)j * nc);
    K.swap(NK);
    used.resize(nc, false);
    for (int i = nc - 1; i >= cap; i--) free_.push_back(i);
    cap = nc;
  }
};

RCPP_MODULE(relstore) {
  class_<RelStore>("RelStore")
  .constructor<int>()
  .method("add", &RelStore::add)
  .method("remove", &RelStore::remove)
  .method("get", &RelStore::get)
  .method("alive", &RelStore::alive)
  .method("size", &RelStore::size)
  .method("capacity", &RelStore::capacity);
}
