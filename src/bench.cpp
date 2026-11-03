// bench.cpp - sw_bench: benchmark harness.
//
// For each problem size N the first N bases of both FASTA sequences are aligned (N x N matrix).
// Each configuration is run once untimed (warm-up) and then `--reps` times; the median is reported.
// Only the matrix fill is timed (the traceback is sequential and takes negligible time).
// Every run's best cell is compared with the sequential result, and results go to a CSV file.
#include <algorithm>
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

struct Timing {
    double median = 0, min = 0, max = 0;
    Best best;            // result of the warm-up run
    bool stable = true;   // every timed run produced the same result as the warm-up run
};

template <class F>
Timing measure(int reps, F&& run_once) {
    Timing t;
    t.best = run_once();  // warm-up, not timed

    std::vector<double> ms;
    for (int r = 0; r < reps; ++r) {
        const double t0 = wtime();
        const Best b = run_once();
        ms.push_back(ms_since(t0));
        if (!same_cell(b, t.best)) t.stable = false;
    }
    std::sort(ms.begin(), ms.end());
    t.min = ms.front();
    t.max = ms.back();
    const size_t k = ms.size();
    t.median = (k % 2) ? ms[k / 2] : 0.5 * (ms[k / 2 - 1] + ms[k / 2]);
    return t;
}

std::vector<int> default_threads() {
    const int maxt = max_threads();
    std::vector<int> v;
    for (int t = 1; t < maxt; t *= 2) v.push_back(t);
    v.push_back(maxt);
    return v;
}

void usage(std::ostream& os) {
    os << "Usage: sw_bench <seqA.fasta> <seqB.fasta> [options]\n"
          "\n"
          "Times the sequential, wavefront and tiled implementations on the first N bases of both\n"
          "sequences, for several N and several thread counts, and writes the results to a CSV file.\n"
          "\n"
          "Options:\n"
          "  --sizes LIST       comma-separated problem sizes N      (default: 1000,2000,4000,8000)\n"
          "  --threads LIST     comma-separated thread counts        (default: 1,2,4,... up to all cores)\n"
          "  --reps N           timed repetitions per configuration  (default: 5)\n"
          "  --tile N           tile size for the tiled version      (default: 64)\n"
          "  --csv FILE         output CSV                           (default: bench_results.csv)\n"
          "  --max-mem-mb N     skip sizes whose matrix exceeds this (default: 2048)\n"
          "  --help             show this message\n";
}

}  // namespace

int main(int argc, char** argv) try {
    Args args(argc, argv, {"--sizes", "--threads", "--reps", "--tile", "--csv", "--max-mem-mb"});
    if (args.help()) {
        usage(std::cout);
        return 0;
    }
    if (args.positional().size() != 2) {
        usage(std::cerr);
        return 1;
    }

    const int reps = args.get_int("--reps", 5);
    const int tile = args.get_int("--tile", 64);
    const long long max_mb = args.get_int("--max-mem-mb", 2048);
    const std::string csv_path = args.get("--csv", "bench_results.csv");
    std::vector<int> sizes = args.get_int_list("--sizes", {1000, 2000, 4000, 8000});
    std::vector<int> threads = args.get_int_list("--threads", default_threads());
    if (reps < 1 || tile < 1) throw std::runtime_error("--reps and --tile must be at least 1");
    for (int t : threads)
        if (t < 1) throw std::runtime_error("--threads values must be at least 1");

    const FastaRecord fa = read_fasta(args.positional()[0]);
    const FastaRecord fb = read_fasta(args.positional()[1]);
    const int available = (int)std::min(fa.sequence.size(), fb.sequence.size());

    // Drop sizes that do not fit the input or the memory budget.
    std::vector<int> usable;
    for (int n : sizes) {
        if (n < 1) continue;
        if (n > available) {
            std::cerr << "skipping size " << n << ": the shorter sequence has only " << available
                      << " bases\n";
        } else if (DPMatrix::bytes_needed(n, n) > (size_t)max_mb * 1024 * 1024) {
            std::cerr << "skipping size " << n << ": matrix needs "
                      << DPMatrix::bytes_needed(n, n) / (1024 * 1024) << " MB, above the " << max_mb
                      << " MB limit\n";
        } else {
            usable.push_back(n);
        }
    }
    if (usable.empty()) throw std::runtime_error("no usable problem sizes");
    std::sort(usable.begin(), usable.end());

    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot write '" + csv_path + "'");
    csv << "size,impl,threads,tile,reps,median_ms,min_ms,max_ms,speedup_vs_seq,efficiency,gcups,verified\n";

    std::cout << "Sequences : " << fa.name << ", " << fb.name << " (" << available << " usable bases)\n";
    std::cout << "Max threads on this machine: " << max_threads() << ", reps per configuration: " << reps
              << ", tile: " << tile << "\n\n";

    const Scoring sc;  // default scoring: +2 / -1 / -2
    DPMatrix H;
    bool all_ok = true;

    auto emit = [&](int n, const std::string& impl, int t, int tile_col, const Timing& tm, double seq_median,
                    bool verified) {
        const double speedup = seq_median / tm.median;
        const double eff = speedup / t;
        const double gcups = (double)n * (double)n / (tm.median * 1e6);
        csv << n << ',' << impl << ',' << t << ',' << tile_col << ',' << reps << ',' << std::fixed
            << std::setprecision(4) << tm.median << ',' << tm.min << ',' << tm.max << ',' << speedup << ','
            << eff << ',' << gcups << ',' << (verified ? 1 : 0) << std::defaultfloat << std::endl;
        std::cout << "N=" << std::setw(6) << n << "  " << std::left << std::setw(11) << impl << std::right
                  << " threads=" << std::setw(2) << t << "  median " << std::setw(10) << std::fixed
                  << std::setprecision(2) << tm.median << " ms  speedup " << std::setw(6)
                  << std::setprecision(2) << speedup << "x  " << std::setw(6) << std::setprecision(3)
                  << gcups << " GCUPS  " << (verified ? "ok" : "MISMATCH") << std::defaultfloat << "\n"
                  << std::flush;
    };

    for (int n : usable) {
        const std::string a = fa.sequence.substr(0, n);
        const std::string b = fb.sequence.substr(0, n);
        H.resize(n, n);

        // Sequential baseline (also provides the reference result).
        const Timing seq = measure(reps, [&] { return fill_sequential(a, b, sc, H); });
        const bool seq_ok = seq.stable;
        all_ok = all_ok && seq_ok;
        emit(n, "sequential", 1, 0, seq, seq.median, seq_ok);

        for (int t : threads) {
            const Timing wf = measure(reps, [&] { return fill_wavefront(a, b, sc, H, t); });
            const bool wf_ok = wf.stable && same_cell(wf.best, seq.best);
            all_ok = all_ok && wf_ok;
            emit(n, "wavefront", t, 0, wf, seq.median, wf_ok);
        }
        for (int t : threads) {
            const Timing tl = measure(reps, [&] { return fill_tiled(a, b, sc, H, t, tile); });
            const bool tl_ok = tl.stable && same_cell(tl.best, seq.best);
            all_ok = all_ok && tl_ok;
            emit(n, "tiled", t, tile, tl, seq.median, tl_ok);
        }
        std::cout << "\n";
    }

    std::cout << "Results written to " << csv_path << "\n";
    if (!all_ok) {
        std::cerr << "ERROR: at least one run disagreed with the sequential result\n";
        return 2;
    }
    return 0;
} catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
}
