#!/usr/bin/env python3
"""
Script 2: Clump Finder – (L, k, t)-Clumps

Identifies k-mers that appear at least t times within
a window of length L using an efficient sliding window approach.

Input:
    python script2_clump_finder.py genomic.fa

Output:
    results/script2_clumps.tsv
"""

import sys
import os
from collections import defaultdict
from Bio import SeqIO

# ---------------- PARAMETERS ----------------
K = 8
T = 3
L = 1000
OUT_DIR = "results"
OUT_FILE = "script2_clumps.tsv"
# --------------------------------------------


def get_kmer(seq, i, k):
    """Extract k-mer starting at position i."""
    return seq[i:i + k]


def main(fasta_file):
    os.makedirs(OUT_DIR, exist_ok=True)

    genome = str(SeqIO.read(fasta_file, "fasta").seq).upper()
    genome_len = len(genome)

    freq = defaultdict(int)
    clumps = set()

    # Initialize first window
    for i in range(L - K + 1):
        kmer = get_kmer(genome, i, K)
        if "N" not in kmer:
            freq[kmer] += 1
            if freq[kmer] >= T:
                clumps.add(kmer)

    # Slide window
    for i in range(1, genome_len - L + 1):
        # Remove outgoing k-mer
        out_kmer = get_kmer(genome, i - 1, K)
        if "N" not in out_kmer:
            freq[out_kmer] -= 1

        # Add incoming k-mer
        in_kmer = get_kmer(genome, i + L - K, K)
        if "N" not in in_kmer:
            freq[in_kmer] += 1
            if freq[in_kmer] >= T:
                clumps.add(in_kmer)

    # Write results
    output_path = os.path.join(OUT_DIR, OUT_FILE)
    with open(output_path, "w") as f:
        f.write("K-mer\n")
        for kmer in sorted(clumps):
            f.write(f"{kmer}\n")

    print(f"[INFO] (L, k, t)-clumps saved → {output_path}")
    print(f"[INFO] Total clump-forming k-mers: {len(clumps)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script2_clump_finder.py genomic.fa")
        sys.exit(1)

    main(sys.argv[1])
