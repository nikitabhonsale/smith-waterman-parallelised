// main.cpp - sw_align: align two FASTA sequences with Smith-Waterman.
//
// Runs the selected implementation(s), checks that all of them agree with the sequential reference,
// prints a report to the console and writes the full report to a file.
#include <algorithm>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cli.h"
#include "fasta.h"
#include "sw.h"

using namespace sw;

namespace {

double ms_since(double t0) { return (wtime() - t0) * 1000.0; }

struct RunResult {
    std::string impl;
    int threads = 1;
    double fill_ms = 0;
    double traceback_ms = 0;
    Best best;
    Alignment aln;
    bool agrees = true;  // matches the first (reference) run
};

struct Report {
    std::string file_a, file_b, name_a, name_b;
    int m = 0, n = 0;
    Scoring sc;
    int tile = 0;
    std::vector<RunResult> runs;
};

void usage(std::ostream& os) {
    os << "Usage: sw_align <seqA.fasta> <seqB.fasta> [options]\n"
          "\n"
          "Smith-Waterman local alignment of two DNA sequences (first record of each FASTA file).\n"
          "\n"
          "Options:\n"
          "  --impl NAME        sequential | wavefront | tiled | all        (default: all)\n"
          "  --threads N        OpenMP threads for the parallel versions    (default: all cores)\n"
          "  --tile N           tile size for the tiled version             (default: 64)\n"
          "  --match N          score for a match, > 0                      (default: 2)\n"
          "  --mismatch N       score for a mismatch, <= 0                  (default: -1)\n"
          "  --gap N            score per gap position, <= 0                (default: -2)\n"
          "  --out FILE         file for the full report                    (default: alignment_result.txt)\n"
          "  --max-mem-mb N     refuse inputs whose matrix needs more MB    (default: 2048)\n"
          "  --help             show this message\n";
}

void write_report(std::ostream& os, const Report& r, int max_blocks) {
    const RunResult& ref = r.runs.front();

    os << "=== Smith-Waterman local alignment ===\n";
    os << "Sequence A : " << r.name_a << " (" << r.m << " bases) from " << r.file_a << "\n";
    os << "Sequence B : " << r.name_b << " (" << r.n << " bases) from " << r.file_b << "\n";
    os << "Scoring    : match " << r.sc.match << ", mismatch " << r.sc.mismatch << ", gap " << r.sc.gap
       << " (linear)\n";
    os << "Matrix     : " << (r.m + 1) << " x " << (r.n + 1) << " cells ("
       << std::fixed << std::setprecision(1) << DPMatrix::bytes_needed(r.m, r.n) / (1024.0 * 1024.0) << " MB)\n\n";

    os << "--- Best local alignment (traced back from the '" << ref.impl << "' matrix) ---\n";
    print_alignment_summary(os, ref.aln);
    os << "\n";
    print_alignment(os, ref.aln, r.name_a, r.name_b, 60, max_blocks);

    os << "--- Timing (matrix fill only; traceback timed separately) ---\n";
    os << std::left << std::setw(12) << "Impl" << std::right << std::setw(8) << "Threads" << std::setw(14)
       << "Fill (ms)" << std::setw(10) << "Speedup" << std::setw(9) << "GCUPS" << std::setw(15)
       << "Traceback (ms)" << "  " << std::left << "Agrees\n";

    const bool have_seq = ref.impl == "sequential";
    const double cells = (double)r.m * (double)r.n;
    for (const RunResult& run : r.runs) {
        os << std::left << std::setw(12) << run.impl << std::right << std::setw(8) << run.threads
           << std::setw(14) << std::fixed << std::setprecision(2) << run.fill_ms;
        if (have_seq && run.fill_ms > 0)
            os << std::setw(9) << std::setprecision(2) << ref.fill_ms / run.fill_ms << "x";
        else
            os << std::setw(10) << "-";
        os << std::setw(9) << std::setprecision(3) << (run.fill_ms > 0 ? cells / (run.fill_ms * 1e6) : 0.0)
           << std::setw(15) << std::setprecision(3) << run.traceback_ms << "  " << std::left
           << (&run == &ref ? "(reference)" : (run.agrees ? "yes" : "NO - MISMATCH")) << "\n";
    }
    os << std::left;
    os << "\nGCUPS = giga cell updates per second (cells filled / fill time).\n";
}

}  // namespace

int main(int argc, char** argv) try {
    Args args(argc, argv, {"--impl", "--threads", "--tile", "--match", "--mismatch", "--gap", "--out",
                           "--max-mem-mb"});
    if (args.help()) {
        usage(std::cout);
        return 0;
    }
    if (args.positional().size() != 2) {
        usage(std::cerr);
        return 1;
    }

    // ---- options -------------------------------------------------------------------------
    Report rep;
    rep.sc.match = args.get_int("--match", 2);
    rep.sc.mismatch = args.get_int("--mismatch", -1);
    rep.sc.gap = args.get_int("--gap", -2);
    if (rep.sc.match <= 0) throw std::runtime_error("--match must be positive");
    if (rep.sc.mismatch > 0) throw std::runtime_error("--mismatch must be zero or negative");
    if (rep.sc.gap > 0) throw std::runtime_error("--gap must be zero or negative");

    const int threads = args.get_int("--threads", max_threads());
    rep.tile = args.get_int("--tile", 64);
    if (threads < 1) throw std::runtime_error("--threads must be at least 1");
    if (rep.tile < 1) throw std::runtime_error("--tile must be at least 1");

    std::string impl = args.get("--impl", "all");
    if (impl == "seq") impl = "sequential";
    std::vector<std::string> impls;
    if (impl == "all")
        impls = {"sequential", "wavefront", "tiled"};
    else if (impl == "sequential" || impl == "wavefront" || impl == "tiled")
        impls = {impl};
    else
        throw std::runtime_error("--impl must be sequential, wavefront, tiled or all");

    const std::string out_path = args.get("--out", "alignment_result.txt");
    const long long max_mb = args.get_int("--max-mem-mb", 2048);

    // ---- input ---------------------------------------------------------------------------
    rep.file_a = args.positional()[0];
    rep.file_b = args.positional()[1];
    const FastaRecord fa = read_fasta(rep.file_a);
    const FastaRecord fb = read_fasta(rep.file_b);
    rep.name_a = fa.name;
    rep.name_b = fb.name;
    rep.m = (int)fa.sequence.size();
    rep.n = (int)fb.sequence.size();

    if ((long long)rep.sc.match * std::min(rep.m, rep.n) > INT_MAX)
        throw std::runtime_error("sequences are too long for 32-bit scores with this match score");

    const size_t need = DPMatrix::bytes_needed(rep.m, rep.n);
    if (need > (size_t)max_mb * 1024 * 1024) {
        throw std::runtime_error(
            "the " + std::to_string(rep.m + 1) + " x " + std::to_string(rep.n + 1) + " matrix needs " +
            std::to_string(need / (1024 * 1024)) + " MB, above the " + std::to_string(max_mb) +
            " MB limit. Use shorter sequences or raise --max-mem-mb if you have the RAM.");
    }

    DPMatrix H;
    try {
        H.resize(rep.m, rep.n);
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("out of memory allocating the score matrix");
    }

    // ---- run -----------------------------------------------------------------------------
    for (const std::string& name : impls) {
        RunResult run;
        run.impl = name;
        run.threads = (name == "sequential") ? 1 : threads;

        double t0 = wtime();
        if (name == "sequential")
            run.best = fill_sequential(fa.sequence, fb.sequence, rep.sc, H);
        else if (name == "wavefront")
            run.best = fill_wavefront(fa.sequence, fb.sequence, rep.sc, H, threads);
        else
            run.best = fill_tiled(fa.sequence, fb.sequence, rep.sc, H, threads, rep.tile);
        run.fill_ms = ms_since(t0);

        // The matrix still holds this implementation's result, so trace back now before the next
        // implementation overwrites it.
        t0 = wtime();
        run.aln = traceback(fa.sequence, fb.sequence, rep.sc, H, run.best);
        run.traceback_ms = ms_since(t0);

        if (!rep.runs.empty()) {
            const RunResult& ref = rep.runs.front();
            run.agrees = same_cell(run.best, ref.best) && same_alignment(run.aln, ref.aln);
        }
        rep.runs.push_back(std::move(run));
    }

    // ---- report --------------------------------------------------------------------------
    write_report(std::cout, rep, /*max_blocks=*/10);

    std::ofstream out(out_path);
    if (!out) throw std::runtime_error("cannot write '" + out_path + "'");
    write_report(out, rep, /*max_blocks=*/0);
    std::cout << "Full report written to " << out_path << "\n";

    for (const RunResult& run : rep.runs) {
        if (!run.agrees) {
            std::cerr << "ERROR: '" << run.impl << "' disagrees with the reference result\n";
            return 2;
        }
    }
    return 0;
} catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
}
