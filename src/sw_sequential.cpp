// sw_sequential.cpp - reference implementation (no OpenMP).
//
// Every cell depends on its upper, left and upper-left neighbours, so filling the matrix row by row
// is the simplest correct order. This is the baseline all parallel versions are checked against.
#include "sw.h"

namespace sw {

Best fill_sequential(const std::string& a, const std::string& b, const Scoring& sc, DPMatrix& H) {
    const int m = (int)a.size();
    const int n = (int)b.size();

    H.init_borders();
    Best best;

    for (int i = 1; i <= m; ++i) {
        const int32_t* prev = H.row(i - 1);
        int32_t* cur = H.row(i);
        const char ai = a[i - 1];
        for (int j = 1; j <= n; ++j) {
            const int32_t v = cell_score(prev[j - 1], prev[j], cur[j - 1],
                                         substitution(ai, b[j - 1], sc), sc.gap);
            cur[j] = v;
            update_best(best, v, i, j);
        }
    }
    return best;
}

}  // namespace sw
