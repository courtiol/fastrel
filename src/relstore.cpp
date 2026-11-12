// Created by Claude Sonnet 5.5 from Alex's prompts

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
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
// so each birth costs O(individuals held) instead of recomputing the whole
// matrix.
//
// Layout: a dense capacity x capacity matrix in column-major order, element
// (i,j) is K[i + j*cap]. Each individual occupies a "slot" (row/column index).
// Removed individuals free their slot for reuse. The most recently freed slot
// is reused first and fresh slots are handed out in increasing order, so the
// occupied slots stay packed at the low indices: the work per birth depends on
// the peak number of individuals held, not on the capacity.
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
    // Everything that can fail is done before the first slot is taken, so
    // that an error leaves the store unchanged (and in step with the R side):
    // validate all the parents, then make room for all the newborns.
    for (int b = 0; b < nb; b++) { slot(mom[b]); slot(dad[b]); }
    IntegerVector out(nb);
    if ((int)free_.size() < nb) grow(nb);

    for (int b = 0; b < nb; b++) {
      const int m = slot(mom[b]), d = slot(dad[b]);   // 0-based, -1 if unknown

      // take a free slot for the newborn
      const int j = free_.back(); free_.pop_back();
      used[j] = true;
      if (j >= hi) hi = j + 1;
      double* cj = col(j);                 // all zero: the slot was free

      // Relatedness of the newborn with every slot i: average of the mother's
      // and the father's relatedness with i (A is symmetric, so a parent's
      // column is also its row). An unknown parent contributes 0. Free slots
      // hold zeros, so they stay zero. A founder (no parents) is unrelated to
      // everyone: its column and row are already zero, nothing to compute.
      if (m >= 0 || d >= 0) {
        if (m >= 0 && d >= 0) {
          const double* pm = col(m);
          const double* pd = col(d);
          for (int i = 0; i < hi; i++) cj[i] = (pm[i] + pd[i]) * 0.5;
        } else {
          const double* pp = col(m >= 0 ? m : d);
          for (int i = 0; i < hi; i++) cj[i] = pp[i] * 0.5;
        }
        // Copy column j into row j, by symmetry. The column write above is
        // contiguous; this one is strided (the slow part).
        double* rj = K.data() + j;
        for (int i = 0; i < hi; i++) rj[(size_t)i * cap] = cj[i];
      }

      // Diagonal: 1 + F_j, with F_j = A[m,d] / 2. For a founder or an
      // individual with a single known parent the diagonal is 1.
      // This overrides the value the averaging produced for i == j.
      cj[j] = (m >= 0 && d >= 0) ? 1.0 + K[m + (size_t)d * cap] * 0.5 : 1.0;

      out[b] = j + 1;                      // back to 1-based for R
    }
    return out;
  }

  // Remove individuals (they can no longer be used as parents). Their row and
  // column are zeroed to keep the invariant, and the slot is recycled.
  // Relatedness values of the remaining individuals are unaffected.
  void remove(IntegerVector slots) {
    for (int s : slots) {
      if (s == NA_INTEGER) continue;       // NA: nothing to remove
      const int j = slot(s);
      double* cj = col(j);
      std::fill(cj, cj + hi, 0.0);
      double* rj = K.data() + j;
      for (int i = 0; i < hi; i++) rj[(size_t)i * cap] = 0.0;
      used[j] = false;
      free_.push_back(j);
    }
    while (hi > 0 && !used[hi - 1]) hi--;
  }

  // Relatedness submatrix for the given slots (already on the A = 2*phi scale).
  NumericMatrix get(IntegerVector slots) {
    const int k = slots.size();
    // convert and validate the slots once, not inside the double loop
    std::vector<int> s(k);
    for (int a = 0; a < k; a++) {
      if (slots[a] == NA_INTEGER) stop("invalid or removed slot: NA");
      s[a] = slot(slots[a]);
    }

    NumericMatrix out(no_init(k, k));      // every cell is written below
    double* o = out.begin();
    for (int a = 0; a < k; a++, o += k) {
      const double* ca = col(s[a]);
      for (int b = 0; b < k; b++) o[b] = ca[s[b]];
    }
    return out;
  }

  // Fill an EMPTY store from a complete relatedness matrix (used to restore a
  // saved store). The k individuals of A (rows/columns, in this order) get
  // fresh slots; the capacity grows if k exceeds it. Returns their 1-based
  // slots, in the order of the rows of A.
  IntegerVector restore(NumericMatrix A) {
    const int k = A.nrow();
    if (A.ncol() != k) stop("A must be a square matrix");
    if (size() != 0) stop("restore() requires an empty store");
    for (int a = 0; a < k; a++)
      for (int b = 0; b <= a; b++) {
        if (!R_finite(A(a, b))) stop("A contains missing or non-finite values");
        if (std::abs(A(a, b) - A(b, a)) > 1e-10) stop("A must be symmetric");
      }
    if (cap < k) grow(k);

    std::vector<int> js(k);
    IntegerVector out(k);
    for (int a = 0; a < k; a++) {
      const int j = free_.back(); free_.pop_back();
      used[j] = true;
      if (j >= hi) hi = j + 1;
      js[a] = j;
      out[a] = j + 1;
    }
    for (int a = 0; a < k; a++)
      for (int b = 0; b < k; b++)
        K[js[b] + (size_t)js[a] * cap] = A(b, a);
    return out;
  }

  // Make room for `need` more individuals in a single step.
  void reserve(int need) {
    if (need > 0 && (int)free_.size() < need) grow(need);
  }

  int size() { return (int)(cap - free_.size()); }   // individuals in the store
  int capacity() { return cap; }                     // slots available (used + free)

private:
  int cap;                    // current number of slots
  std::vector<double> K;      // cap x cap relatedness matrix, column-major
  std::vector<char> used;     // used[i] = slot i is occupied
  std::vector<int> free_;     // stack of free slots
  int hi;                     // all the slots >= hi are free, so rows/columns
                              // only need to be read or written below hi

  void init(int c) {
    if (c < 1) stop("capacity must be at least 1");
    cap = c;
    hi = 0;
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

  // Enlarge the matrix so that at least `need` slots are free: copy the old
  // matrix into the top-left corner of a larger, zero-filled one and register
  // the new slots as free. The capacity grows by half (not by doubling, which
  // could leave a matrix four times larger than needed), or straight to the
  // size required if that is more, so a large batch costs a single copy.
  // The new slots go to the bottom of the stack, so that the free slots with
  // a low index are still used first. All the allocations are done before the
  // store is modified: if memory runs out, the store is left as it was.
  void grow(int need) {
    const double want = std::max(cap * 1.5 + 1.0, (double)size() + need);
    if (want > INT_MAX) stop("the store cannot grow any further");
    const int nc = (int)want;
    std::vector<double> NK((size_t)nc * nc, 0.0);
    std::vector<int> nf;
    nf.reserve(free_.size() + (nc - cap));
    used.reserve(nc);
    for (int i = nc - 1; i >= cap; i--) nf.push_back(i);
    nf.insert(nf.end(), free_.begin(), free_.end());
    for (int j = 0; j < hi; j++)
      std::copy(K.begin() + (size_t)j * cap, K.begin() + (size_t)j * cap + hi,
                NK.begin() + (size_t)j * nc);
    K.swap(NK);
    free_.swap(nf);
    used.resize(nc, false);
    cap = nc;
  }
};

RCPP_MODULE(relstore) {
  class_<RelStore>("RelStore")
  .constructor<int>()
  .method("add", &RelStore::add)
  .method("remove", &RelStore::remove)
  .method("get", &RelStore::get)
  .method("restore", &RelStore::restore)
  .method("reserve", &RelStore::reserve)
  .method("size", &RelStore::size)
  .method("capacity", &RelStore::capacity);
}
