# Wavefront Smith-Waterman DNA Alignment (OpenMP)

Local alignment of two DNA sequences (Smith-Waterman), with three interchangeable matrix-fill implementations:

| Implementation | File | Idea |
|---|---|---|
| `sequential` | `src/sw_sequential.cpp` | Plain row-by-row loops. The reference everything is checked against. |
| `wavefront` | `src/sw_wavefront.cpp` | Cells on one anti-diagonal are independent: split them across threads, barrier between diagonals. |
| `tiled` | `src/sw_tiled.cpp` | Same idea on tile x tile blocks: tiles on one tile-diagonal run in parallel, each tile is computed sequentially. Far fewer barriers and better cache use. |

Scoring is linear-gap: match +2, mismatch -1, gap -2 per position (all configurable). Input is FASTA only.

## Build

**macOS** (Apple's clang has no built-in OpenMP):

```
brew install libomp
make
```

Alternative with GCC: `brew install gcc`, then `make CXX=g++-14` (use whatever version you installed).

**Linux:** `make`

Programs are built into `bin/`: `sw_align`, `sw_bench`, `sw_tests`.

## Run

```
make test     # correctness tests
make demo     # align the small example pair
make bench    # benchmark the large example pair -> bench_results.csv
```

Align your own files:

```
./bin/sw_align a.fasta b.fasta                       # runs all three, checks they agree
./bin/sw_align a.fasta b.fasta --impl tiled --threads 4 --tile 128
./bin/sw_align a.fasta b.fasta --match 1 --mismatch -1 --gap -1 --out result.txt
```

The console shows the first 10 blocks of the alignment; the output file (`alignment_result.txt` by default) has all of it.

Benchmark (median of `--reps` runs, one CSV row per configuration):

```
./bin/sw_bench examples/large_a.fasta examples/large_b.fasta --sizes 2000,4000,8000 --threads 1,2,4,8 --reps 5
python3 tools/plot_results.py bench_results.csv --outdir plots     # needs matplotlib
```

The benchmark aligns the first N bases of each file, so one pair of long FASTA files covers every size.
`tools/generate_fasta.py` creates related test pairs of any length.

## Memory

The score matrix is (m+1) x (n+1) 32-bit integers. The programs refuse inputs needing more than 2048 MB
(`--max-mem-mb`), which is roughly 23,000 x 23,000 bases. That is a safe limit on an 8 GB machine.

## How correctness is guaranteed

* All implementations call the same `cell_score()` function (`src/sw.h`).
* The best cell is chosen by a fixed rule (highest score, then smallest row, then smallest column) and traceback prefers
  diagonal, then up, then left, so all implementations produce byte-identical alignments.
* `sw_tests` checks empty sequences, identical sequences, sequences with nothing in common, a textbook example, and
  random pairs (including thin/wide shapes, tile size 1, and several scoring schemes) across 16 implementation/thread/tile
  configurations. Each alignment is also re-scored independently from its aligned strings.

## Layout

```
src/sw.h              shared types: scoring, DP matrix, cell recurrence, declarations
src/sw_sequential.cpp / sw_wavefront.cpp / sw_tiled.cpp
src/alignment.cpp     traceback, alignment checking and printing
src/fasta.cpp         FASTA reader
src/main.cpp          sw_align (command-line aligner)
src/bench.cpp         sw_bench (benchmark harness -> CSV)
src/tests.cpp         sw_tests
tools/                generate_fasta.py, plot_results.py
examples/             small pair (demo) and large pair (benchmark)
```
