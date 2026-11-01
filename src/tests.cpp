// tests.cpp - sw_tests: correctness tests.
//
// Every test runs all implementations (sequential, wavefront and tiled, with several thread counts
// and tile sizes, including awkward ones such as tile = 1) and checks that
//   * each alignment is valid: its gap-free letters are the stated substrings of the inputs and its
//     score, recomputed independently from the aligned strings, equals the reported score,
//   * all implementations give exactly the same best cell and the same alignment as the sequential
//     reference.
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "sw.h"

using namespace sw;

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::cerr << "    FAIL: " << what << "\n";
    }
}

// ---- running one implementation ----------------------------------------------------------------

struct Config {
    std::string impl;
    int threads;
    int tile;
    std::string label() const {
        return impl + "(threads=" + std::to_string(threads) + (impl == "tiled" ? ", tile=" + std::to_string(tile) : "") + ")";
    }
};

struct Outcome {
    Best best;
    Alignment aln;
};

Outcome run(const Config& c, const std::string& a, const std::string& b, const Scoring& sc) {
    DPMatrix H;
    H.resize((int)a.size(), (int)b.size());
    Outcome o;
    if (c.impl == "sequential")
        o.best = fill_sequential(a, b, sc, H);
    else if (c.impl == "wavefront")
        o.best = fill_wavefront(a, b, sc, H, c.threads);
    else
        o.best = fill_tiled(a, b, sc, H, c.threads, c.tile);
    o.aln = traceback(a, b, sc, H, o.best);
    return o;
}

std::vector<Config> all_configs() {
    std::vector<Config> v = {{"sequential", 1, 0}};
    for (int t : {1, 2, 4}) v.push_back({"wavefront", t, 0});
    for (int t : {1, 2, 4})
        for (int tile : {1, 3, 8, 64}) v.push_back({"tiled", t, tile});
    return v;
}

// Independent validity check of one outcome.
void validate(const std::string& a, const std::string& b, const Scoring& sc, const Outcome& o,
              const std::string& label) {
    const Alignment& aln = o.aln;
    check(aln.score == o.best.score, label + ": alignment score equals best-cell score");
    if (aln.length() == 0) {
        check(o.best.score == 0, label + ": empty alignment only when the score is 0");
        return;
    }
    check(rescore(aln, sc) == aln.score, label + ": rescoring the aligned strings gives the same score");

    std::string ga, gb;  // gap-free letters
    bool double_gap = false;
    for (size_t k = 0; k < aln.a_aligned.size(); ++k) {
        if (aln.a_aligned[k] == '-' && aln.b_aligned[k] == '-') double_gap = true;
        if (aln.a_aligned[k] != '-') ga.push_back(aln.a_aligned[k]);
        if (aln.b_aligned[k] != '-') gb.push_back(aln.b_aligned[k]);
    }
    check(!double_gap, label + ": no column is a gap in both sequences");
    check(aln.a_start >= 1 && aln.a_end <= (int)a.size() && ga == a.substr(aln.a_start - 1, aln.a_end - aln.a_start + 1),
          label + ": A letters match the stated range");
    check(aln.b_start >= 1 && aln.b_end <= (int)b.size() && gb == b.substr(aln.b_start - 1, aln.b_end - aln.b_start + 1),
          label + ": B letters match the stated range");
}

// Runs every configuration, validates each, and compares them all with the sequential result.
void run_all_and_compare(const std::string& a, const std::string& b, const Scoring& sc, const std::string& name) {
    const std::vector<Config> configs = all_configs();
    const Outcome ref = run(configs.front(), a, b, sc);
    validate(a, b, sc, ref, name + " / " + configs.front().label());
    for (size_t k = 1; k < configs.size(); ++k) {
        const Outcome o = run(configs[k], a, b, sc);
        const std::string label = name + " / " + configs[k].label();
        validate(a, b, sc, o, label);
        check(same_cell(o.best, ref.best), label + ": same best cell as sequential");
        check(same_alignment(o.aln, ref.aln), label + ": same alignment as sequential");
    }
}

// ---- random data -------------------------------------------------------------------------------

std::mt19937 g_rng(20261101);

std::string random_dna(int len) {
    static const char bases[] = "ACGT";
    std::string s;
    for (int k = 0; k < len; ++k) s.push_back(bases[g_rng() % 4]);
    return s;
}

// Copy of `s` with substitutions, insertions and deletions, so the pair has a long, gappy alignment.
std::string mutate(const std::string& s, double rate) {
    static const char bases[] = "ACGT";
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::string out;
    for (char c : s) {
        const double r = u(g_rng);
        if (r < rate / 4) continue;                              // deletion
        if (r < rate / 2) out.push_back(bases[g_rng() % 4]);     // insertion before c
        out.push_back(r < rate ? bases[g_rng() % 4] : c);        // possible substitution
    }
    return out;
}

// ---- tests -------------------------------------------------------------------------------------

void test_empty_sequences() {
    const Scoring sc;
    const std::vector<std::pair<std::string, std::string>> cases = {{"", ""}, {"", "ACGT"}, {"ACGT", ""}};
    for (const auto& [a, b] : cases) {
        for (const Config& c : all_configs()) {
            const std::string label = "A='" + a + "' B='" + b + "' / " + c.label();
            const Outcome o = run(c, a, b, sc);
            check(o.best.score == 0, label + ": score is 0");
            check(o.aln.length() == 0, label + ": alignment is empty");
            check(o.aln.a_aligned.empty() && o.aln.b_aligned.empty(), label + ": no aligned letters");
        }
    }
}

void test_identical_sequences() {
    const Scoring sc;
    for (int len : {1, 2, 17, 100, 1000}) {
        const std::string s = random_dna(len);
        for (const Config& c : all_configs()) {
            const std::string label = "identical length " + std::to_string(len) + " / " + c.label();
            const Outcome o = run(c, s, s, sc);
            check(o.best.score == sc.match * len, label + ": score is match * length");
            check(o.aln.a_aligned == s && o.aln.b_aligned == s, label + ": alignment is the whole sequence");
            check(o.aln.matches == len && o.aln.mismatches == 0 && o.aln.gaps == 0,
                  label + ": 100% identity, no gaps");
            check(o.aln.a_start == 1 && o.aln.a_end == len && o.aln.b_start == 1 && o.aln.b_end == len,
                  label + ": covers positions 1.." + std::to_string(len));
        }
    }
}

void test_no_similarity() {
    const Scoring sc;
    for (const Config& c : all_configs()) {
        Outcome o = run(c, "AAAAAA", "CCCCCC", sc);
        check(o.best.score == 0 && o.aln.length() == 0, "no shared bases / " + c.label() + ": score 0, empty alignment");
        o = run(c, "NNNN", "NNNN", sc);
        check(o.best.score == 0 && o.aln.length() == 0, "N never matches / " + c.label() + ": score 0, empty alignment");
    }
}

// Worked example from the Wikipedia article on the Smith-Waterman algorithm:
// match +3, mismatch -3, gap -2 per position; best local alignment scores 13 as GTT-AC / GTTGAC.
void test_known_example() {
    Scoring sc;
    sc.match = 3;
    sc.mismatch = -3;
    sc.gap = -2;
    const std::string a = "TGTTACGG";
    const std::string b = "GGTTGACTA";
    for (const Config& c : all_configs()) {
        const Outcome o = run(c, a, b, sc);
        const std::string label = "textbook example / " + c.label();
        check(o.best.score == 13, label + ": score is 13");
        check(o.aln.a_aligned == "GTT-AC" && o.aln.b_aligned == "GTTGAC", label + ": alignment is GTT-AC / GTTGAC");
    }
}

void test_random_pairs() {
    const std::vector<Scoring> schemes = {{2, -1, -2}, {1, -1, -1}, {5, -4, -8}, {3, -3, -2}};

    // Hand-picked shapes: very thin and very wide matrices, and a square one.
    const std::vector<std::pair<int, int>> shapes = {{1, 200}, {200, 1}, {300, 5}, {5, 300}, {120, 120}};
    int idx = 0;
    for (const auto& [m, n] : shapes) {
        const std::string a = random_dna(m);
        const std::string b = random_dna(n);
        run_all_and_compare(a, b, schemes[idx++ % schemes.size()], "shape " + std::to_string(m) + "x" + std::to_string(n));
    }

    // Random lengths, alternating between unrelated pairs and mutated copies.
    for (int iter = 0; iter < 120; ++iter) {
        const int m = 1 + (int)(g_rng() % 250);
        const std::string a = random_dna(m);
        const std::string b = (iter % 2 == 0) ? mutate(a, 0.05 + 0.25 * (g_rng() % 100) / 100.0)
                                              : random_dna(1 + (int)(g_rng() % 250));
        if (b.empty()) continue;
        run_all_and_compare(a, b, schemes[iter % schemes.size()], "random pair #" + std::to_string(iter));
    }
}

void run_test(const char* name, void (*fn)()) {
    const int before = g_failures;
    fn();
    std::cout << (g_failures == before ? "[ ok ] " : "[FAIL] ") << name << "\n";
}

}  // namespace

int main() {
    std::cout << "Running correctness tests (" << all_configs().size()
              << " implementation configurations per case)\n\n";

    run_test("empty sequences", test_empty_sequences);
    run_test("identical sequences", test_identical_sequences);
    run_test("sequences with nothing in common", test_no_similarity);
    run_test("textbook example (score 13)", test_known_example);
    run_test("random pairs vs sequential reference", test_random_pairs);

    std::cout << "\n" << (g_checks - g_failures) << " / " << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
