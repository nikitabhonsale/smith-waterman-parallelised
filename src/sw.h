// sw.h - shared types and declarations for the wavefront Smith-Waterman project.
//
// Conventions used everywhere in the project:
//   * Sequence A (length m) runs down the ROWS of the DP matrix, sequence B (length n) across the COLUMNS.
//   * H(i, j) is the best local-alignment score ending at A[i-1] and B[j-1] (1-based matrix indices).
//   * Row 0 and column 0 are all zeros (a local alignment may start anywhere).
//   * Linear gap penalty: every gap position costs `gap` (a negative number).
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
#endif

namespace sw {

inline int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// Wall-clock time in seconds (omp_get_wtime when built with OpenMP). Only differences of two
// readings taken on the same thread are meaningful.
inline double wtime() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
#endif
}

// ---------------------------------------------------------------------------------------------
// Scoring
// ---------------------------------------------------------------------------------------------
struct Scoring {
    int match = 2;      // reward for identical bases (positive)
    int mismatch = -1;  // penalty for different bases (negative)
    int gap = -2;       // penalty per gap position (negative)
};

// 'N' (unknown base) never counts as a match, even against another 'N'.
inline int substitution(char x, char y, const Scoring& sc) {
    return (x == y && x != 'N') ? sc.match : sc.mismatch;
}

// ---------------------------------------------------------------------------------------------
// The recurrence for one cell. Every implementation calls this same function, so the sequential
// and parallel versions cannot drift apart.
//     H(i,j) = max( 0, H(i-1,j-1) + s(a_i,b_j), H(i-1,j) + gap, H(i,j-1) + gap )
// ---------------------------------------------------------------------------------------------
inline int32_t cell_score(int32_t diag, int32_t up, int32_t left, int sub, int gap) {
    int32_t v = diag + sub;
    const int32_t u = up + gap;
    const int32_t l = left + gap;
    if (u > v) v = u;
    if (l > v) v = l;
    return v > 0 ? v : 0;
}

// ---------------------------------------------------------------------------------------------
// Best cell found so far. Ties on score are broken by the smallest row, then the smallest column.
// Because this rule does not depend on the order cells are visited, every implementation reports
// exactly the same cell, no matter how the work is split across threads.
// ---------------------------------------------------------------------------------------------
struct Best {
    int32_t score = 0;
    int i = 0;
    int j = 0;
};

inline void update_best(Best& b, int32_t score, int i, int j) {
    if (score > b.score || (score == b.score && (i < b.i || (i == b.i && j < b.j)))) {
        b.score = score;
        b.i = i;
        b.j = j;
    }
}

inline bool same_cell(const Best& x, const Best& y) {
    return x.score == y.score && x.i == y.i && x.j == y.j;
}

// ---------------------------------------------------------------------------------------------
// Dense (m+1) x (n+1) matrix of 32-bit scores, row-major. The buffer is reused between runs when it
// is big enough, and its pages are touched once on allocation so that first-use page faults are
// not charged to whichever implementation happens to run first.
// ---------------------------------------------------------------------------------------------
class DPMatrix {
public:
    static size_t bytes_needed(int m, int n) {
        return (size_t)(m + 1) * (size_t)(n + 1) * sizeof(int32_t);
    }

    void resize(int m, int n) {
        const size_t cells = (size_t)(m + 1) * (size_t)(n + 1);
        if (cells > capacity_) {
            data_.reset();
            capacity_ = 0;
            data_.reset(new int32_t[cells]);  // may throw std::bad_alloc
            capacity_ = cells;
            std::memset(data_.get(), 0, cells * sizeof(int32_t));
        }
        m_ = m;
        n_ = n;
        stride_ = (size_t)n + 1;
    }

    // Row 0 and column 0 must be zero before a fill; every fill function calls this first.
    void init_borders() {
        std::memset(row(0), 0, stride_ * sizeof(int32_t));
        for (int i = 1; i <= m_; ++i) row(i)[0] = 0;
    }

    int rows() const { return m_; }
    int cols() const { return n_; }
    int32_t* row(int i) { return data_.get() + (size_t)i * stride_; }
    const int32_t* row(int i) const { return data_.get() + (size_t)i * stride_; }
    int32_t at(int i, int j) const { return row(i)[j]; }

private:
    std::unique_ptr<int32_t[]> data_;
    size_t capacity_ = 0;  // in cells
    int m_ = 0;
    int n_ = 0;
    size_t stride_ = 1;
};

// ---------------------------------------------------------------------------------------------
// Phase 1: fill the matrix and return the best cell. Three interchangeable implementations.
// ---------------------------------------------------------------------------------------------

// Reference implementation: plain row-by-row loops, no OpenMP.
Best fill_sequential(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H);

// Wavefront: one anti-diagonal at a time, cells of a diagonal split across threads, barrier between
// diagonals.
Best fill_wavefront(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H,
                    int threads);

// Tiled wavefront: the matrix is cut into tile x tile blocks; tiles on the same anti-diagonal of the
// tile grid run in parallel, each tile is computed sequentially, barrier between tile diagonals.
Best fill_tiled(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H,
                int threads, int tile);

// ---------------------------------------------------------------------------------------------
// Phase 2: traceback and presentation
// ---------------------------------------------------------------------------------------------
struct Alignment {
    std::string a_aligned;  // sequence A with '-' for gaps
    std::string b_aligned;  // sequence B with '-' for gaps
    std::string markers;    // '|' match, '.' mismatch, ' ' gap
    int32_t score = 0;
    int a_start = 0, a_end = 0;  // 1-based inclusive positions in A (0 when the alignment is empty)
    int b_start = 0, b_end = 0;
    int matches = 0, mismatches = 0, gaps = 0;

    int length() const { return (int)a_aligned.size(); }
    double identity() const { return length() ? 100.0 * matches / length() : 0.0; }
};

// Walks back from the best cell until a zero is reached. Fixed tie-break order when several
// predecessors give the same score: diagonal, then up, then left.
Alignment traceback(const std::string& a, const std::string& b, const Scoring& sc, const DPMatrix& H,
                    const Best& best);

// Recomputes an alignment's score from its aligned strings (independent check used by the tests).
int rescore(const Alignment& aln, const Scoring& sc);

bool same_alignment(const Alignment& x, const Alignment& y);

void print_alignment_summary(std::ostream& os, const Alignment& aln);

// Prints the alignment in blocks of `width` columns. max_blocks <= 0 means "print everything".
void print_alignment(std::ostream& os, const Alignment& aln, const std::string& name_a,
                     const std::string& name_b, int width, int max_blocks);

}  // namespace sw
