// sw_tiled.cpp - tiled (blocked) wavefront with OpenMP.
//
// The basic wavefront synchronises once per cell-diagonal (m + n barriers) and walks memory along
// anti-diagonals, which is unfriendly to the cache. Here the matrix is cut into tile x tile blocks:
//   * tiles on the same anti-diagonal of the TILE grid are independent -> run in parallel,
//   * each tile is computed sequentially, row by row (cache friendly),
//   * there is one barrier per tile diagonal, i.e. about (m + n) / tile barriers.
//
// A tile (ti, tj) needs the bottom row of tile (ti-1, tj), the right column of tile (ti, tj-1) and
// the corner of tile (ti-1, tj-1). All three live on earlier tile diagonals, so they are finished.
#include <algorithm>

#include "sw.h"

namespace sw {

namespace {

// Computes the cells with rows i0..i1 and columns j0..j1 (1-based, inclusive), row by row.
inline void compute_tile(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H,
                         int i0, int i1, int j0, int j1, Best& local) {
    for (int i = i0; i <= i1; ++i) {
        const int32_t* prev = H.row(i - 1);
        int32_t* cur = H.row(i);
        const char ai = a[i - 1];
        for (int j = j0; j <= j1; ++j) {
            const int32_t v = cell_score(prev[j - 1], prev[j], cur[j - 1],
                                         substitution(ai, b[j - 1], sc), sc.gap);
            cur[j] = v;
            update_best(local, v, i, j);
        }
    }
}

}  // namespace

Best fill_tiled(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H,
                int threads, int tile) {
    const int m = (int)a.size();
    const int n = (int)b.size();

    H.init_borders();
    Best best;
    if (m == 0 || n == 0) return best;

    const int tiles_m = (m + tile - 1) / tile;  // tile rows
    const int tiles_n = (n + tile - 1) / tile;  // tile columns

#pragma omp parallel num_threads(threads)
    {
        Best local;

        for (int td = 0; td < tiles_m + tiles_n - 1; ++td) {  // td = ti + tj
            const int ti_lo = std::max(0, td - (tiles_n - 1));
            const int ti_hi = std::min(tiles_m - 1, td);

#pragma omp for schedule(static)
            for (int ti = ti_lo; ti <= ti_hi; ++ti) {
                const int tj = td - ti;
                const int i0 = ti * tile + 1;
                const int i1 = std::min(m, (ti + 1) * tile);
                const int j0 = tj * tile + 1;
                const int j1 = std::min(n, (tj + 1) * tile);
                compute_tile(a, b, sc, H, i0, i1, j0, j1, local);
            }
            // <- implicit barrier: tile diagonal td is complete before td + 1 starts
        }

#pragma omp critical
        update_best(best, local.score, local.i, local.j);
    }
    return best;
}

}  // namespace sw
