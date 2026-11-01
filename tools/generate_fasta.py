#!/usr/bin/env python3
"""Generate a pair of related DNA sequences as FASTA files (for demos and benchmarks).

The first sequence is random. The second is a mutated copy (substitutions plus a few insertions
and deletions), so the two have a long, gappy local alignment.

Example:
    python3 tools/generate_fasta.py --length 400 --out-a examples/small_a.fasta \
        --out-b examples/small_b.fasta --name small --seed 1
"""
import argparse
import random

BASES = "ACGT"


def write_fasta(path, header, seq, width=70):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for k in range(0, len(seq), width):
            f.write(seq[k:k + width] + "\n")


def mutate(seq, sub_rate, indel_rate, rng):
    out = []
    for base in seq:
        r = rng.random()
        if r < indel_rate / 2:
            continue                                  # deletion of this base
        if r < indel_rate:
            out.append(rng.choice(BASES))             # insertion before this base
        if rng.random() < sub_rate:
            base = rng.choice([b for b in BASES if b != base])  # substitution
        out.append(base)
    return "".join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--length", type=int, default=10000, help="length of the first sequence")
    p.add_argument("--sub-rate", type=float, default=0.05, help="per-base substitution probability")
    p.add_argument("--indel-rate", type=float, default=0.01, help="per-base insertion/deletion probability")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--name", default="seq", help="prefix for the FASTA header names")
    p.add_argument("--out-a", default="seq_a.fasta")
    p.add_argument("--out-b", default="seq_b.fasta")
    args = p.parse_args()

    rng = random.Random(args.seed)
    a = "".join(rng.choice(BASES) for _ in range(args.length))
    b = mutate(a, args.sub_rate, args.indel_rate, rng)

    write_fasta(args.out_a, f"{args.name}_a random sequence, seed {args.seed}", a)
    write_fasta(args.out_b, f"{args.name}_b mutated copy of {args.name}_a", b)
    print(f"wrote {args.out_a} ({len(a)} bases) and {args.out_b} ({len(b)} bases)")


if __name__ == "__main__":
    main()
