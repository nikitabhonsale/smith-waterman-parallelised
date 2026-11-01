// alignment.cpp - traceback, alignment checking and pretty-printing.
#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include "sw.h"

namespace sw {

// ---------------------------------------------------------------------------------------------
// Traceback. Starts at the best cell and follows, at each step, the neighbour the score came from,
// until a cell with score 0 is reached (that is where the local alignment begins).
// When more than one neighbour explains the score, the first match in this order wins:
//     1. diagonal (A[i] aligned with B[j])
//     2. up       (A[i] aligned with a gap)
//     3. left     (B[j] aligned with a gap)
// The matrix is identical for every implementation, so the alignment is too.
// ---------------------------------------------------------------------------------------------
Alignment traceback(const std::string& a, const std::string& b, const Scoring& sc, const DPMatrix& H,
                    const Best& best) {
    Alignment aln;
    aln.score = best.score;
    if (best.score <= 0) return aln;  // nothing aligns: empty alignment

    std::string ra, rb;  // built back to front, reversed at the end
    int i = best.i;
    int j = best.j;

    while (i > 0 && j > 0 && H.at(i, j) > 0) {
        const int32_t h = H.at(i, j);
        if (h == H.at(i - 1, j - 1) + substitution(a[i - 1], b[j - 1], sc)) {
            ra.push_back(a[i - 1]);
            rb.push_back(b[j - 1]);
            --i;
            --j;
        } else if (h == H.at(i - 1, j) + sc.gap) {
            ra.push_back(a[i - 1]);
            rb.push_back('-');
            --i;
        } else if (h == H.at(i, j - 1) + sc.gap) {
            ra.push_back('-');
            rb.push_back(b[j - 1]);
            --j;
        } else {
            throw std::logic_error("traceback: matrix is inconsistent at a non-zero cell");
        }
    }
    std::reverse(ra.begin(), ra.end());
    std::reverse(rb.begin(), rb.end());

    aln.a_aligned = ra;
    aln.b_aligned = rb;
    aln.a_start = i + 1;  // (i, j) is the zero cell just before the first aligned pair
    aln.b_start = j + 1;
    aln.a_end = best.i;
    aln.b_end = best.j;

    aln.markers.resize(ra.size());
    for (size_t k = 0; k < ra.size(); ++k) {
        if (ra[k] == '-' || rb[k] == '-') {
            aln.markers[k] = ' ';
            ++aln.gaps;
        } else if (ra[k] == rb[k] && ra[k] != 'N') {
            aln.markers[k] = '|';
            ++aln.matches;
        } else {
            aln.markers[k] = '.';
            ++aln.mismatches;
        }
    }
    return aln;
}

int rescore(const Alignment& aln, const Scoring& sc) {
    int total = 0;
    for (size_t k = 0; k < aln.a_aligned.size(); ++k) {
        const char x = aln.a_aligned[k];
        const char y = aln.b_aligned[k];
        total += (x == '-' || y == '-') ? sc.gap : substitution(x, y, sc);
    }
    return total;
}

bool same_alignment(const Alignment& x, const Alignment& y) {
    return x.score == y.score && x.a_aligned == y.a_aligned && x.b_aligned == y.b_aligned &&
           x.a_start == y.a_start && x.a_end == y.a_end && x.b_start == y.b_start &&
           x.b_end == y.b_end;
}

void print_alignment_summary(std::ostream& os, const Alignment& aln) {
    os << "  Score        : " << aln.score << "\n";
    if (aln.length() == 0) {
        os << "  Alignment    : none (no positive-scoring local alignment)\n";
        return;
    }
    os << "  Seq A range  : " << aln.a_start << " - " << aln.a_end << "\n";
    os << "  Seq B range  : " << aln.b_start << " - " << aln.b_end << "\n";
    os << "  Length       : " << aln.length() << " columns\n";
    os << "  Matches      : " << aln.matches << " (" << std::fixed << std::setprecision(1)
       << aln.identity() << "% identity)\n";
    os << "  Mismatches   : " << aln.mismatches << "\n";
    os << "  Gap columns  : " << aln.gaps << "\n";
}

void print_alignment(std::ostream& os, const Alignment& aln, const std::string& name_a,
                     const std::string& name_b, int width, int max_blocks) {
    if (aln.length() == 0) return;

    const int len = aln.length();
    const int blocks = (len + width - 1) / width;
    const int shown = (max_blocks > 0) ? std::min(blocks, max_blocks) : blocks;
    const int label_w = 10;

    const std::string la = name_a.substr(0, label_w);
    const std::string lb = name_b.substr(0, label_w);

    auto count_bases = [](const std::string& s) {
        return (int)std::count_if(s.begin(), s.end(), [](char c) { return c != '-'; });
    };

    int pos_a = aln.a_start;
    int pos_b = aln.b_start;
    for (int blk = 0; blk < shown; ++blk) {
        const int off = blk * width;
        const int cols = std::min(width, len - off);
        const std::string sa = aln.a_aligned.substr(off, cols);
        const std::string sb = aln.b_aligned.substr(off, cols);
        const std::string mk = aln.markers.substr(off, cols);
        const int ca = count_bases(sa);
        const int cb = count_bases(sb);

        os << std::left << std::setw(label_w) << la << ' ' << std::right << std::setw(7) << pos_a << ' '
           << sa << ' ' << (pos_a + ca - 1) << "\n";
        os << std::string(label_w + 1 + 7 + 1, ' ') << mk << "\n";
        os << std::left << std::setw(label_w) << lb << ' ' << std::right << std::setw(7) << pos_b << ' '
           << sb << ' ' << (pos_b + cb - 1) << "\n\n";

        pos_a += ca;
        pos_b += cb;
    }
    os << std::left;
    if (shown < blocks) {
        os << "  ... " << (blocks - shown) << " more block(s) of " << width
           << " columns not shown (the output file has the full alignment)\n\n";
    }
}

}  // namespace sw
