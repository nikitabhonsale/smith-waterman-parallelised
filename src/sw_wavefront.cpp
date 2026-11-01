// sw_wavefront.cpp - anti-diagonal (wavefront) parallelisation with OpenMP.
//
// Cells with the same i + j lie on one anti-diagonal. Each cell needs only cells from the two
// previous anti-diagonals, so all cells of one diagonal are independent and can be computed at the
// same time. The diagonals themselves must be processed in order, which is what the implicit
// barrier at the end of each `omp for` guarantees.
//
//        j ->  0  1  2  3  4
//   i = 0      .  .  .  .  .
//   i = 1      .  1  2  3  4      numbers = anti-diagonal index d - 1; cells with equal numbers
//   i = 2      .  2  3  4  5      are computed in parallel
//   i = 3      .  3  4  5  6
#include <algorithm>

#include "sw.h"

namespace sw {

Best fill_wavefront(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H,
                    int threads) {
    const int m = (int)a.size();
    const int n = (int)b.size();

    H.init_borders();
    Best best;
    if (m == 0 || n == 0) return best;

#pragma omp parallel num_threads(threads)
    {
        Best local;  // each thread tracks its own best cell; merged once at the end

        // Every thread walks the same sequence of diagonals d = i + j (2 .. m+n).
        for (int d = 2; d <= m + n; ++d) {
            const int i_lo = std::max(1, d - n);  // first row that has a cell on this diagonal
            const int i_hi = std::min(m, d - 1);  // last row

#pragma omp for schedule(static)
            for (int i = i_lo; i <= i_hi; ++i) {
                const int j = d - i;
                const int32_t v = cell_score(H.row(i - 1)[j - 1], H.row(i - 1)[j], H.row(i)[j - 1],
                                             substitution(a[i - 1], b[j - 1], sc), sc.gap);
                H.row(i)[j] = v;
                update_best(local, v, i, j);
            }
            // <- implicit barrier: diagonal d is complete before any thread starts d + 1
        }

#pragma omp critical
        update_best(best, local.score, local.i, local.j);
    }
    return best;
}

}  // namespace sw
