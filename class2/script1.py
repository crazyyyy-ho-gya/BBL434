#!/usr/bin/env python3
"""
Script 1: ORI Signal Checker – K-mer Enrichment Analysis

Identifies enriched 8-mers across a genome using sliding windows
and generates an enrichment plot.

Input:
    python script1_kmer_enrichment.py genomic.fa

Output:
    results/plots/script1_plot1.png
"""

import sys
import os
from collections import Counter, defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt

# ---------------- PARAMETERS ----------------
K = 8
WINDOW_SIZE = 5000
STEP_SIZE = 500
TOP_KMERS = 5
OUT_DIR = "results/plots"
PLOT_NAME = "script1_plot1.png"
# --------------------------------------------


def count_kmers(seq, k):
    """Return k-mer counts for a sequence."""
    return Counter(
        seq[i:i + k]
        for i in range(len(seq) - k + 1)
        if "N" not in seq[i:i + k]
    )


def sliding_windows(genome):
    """Generate window center positions and k-mer counts."""
    positions, counts = [], []
    for start in range(0, len(genome) - WINDOW_SIZE + 1, STEP_SIZE):
        window = genome[start:start + WINDOW_SIZE]
        positions.append(start + WINDOW_SIZE // 2)
        counts.append(count_kmers(window, K))
    return positions, counts


def main(fasta_file):
    os.makedirs(OUT_DIR, exist_ok=True)

    genome = str(SeqIO.read(fasta_file, "fasta").seq).upper()

    positions, window_counts = sliding_windows(genome)

    global_counts = Counter()
    for wc in window_counts:
        global_counts.update(wc)

    top_kmers = [k for k, _ in global_counts.most_common(TOP_KMERS)]

    enrichment = defaultdict(list)
    for wc in window_counts:
        total = sum(wc.values())
        for kmer in top_kmers:
            enrichment[kmer].append(wc[kmer] / total if total else 0)

    # Plot
    plt.figure(figsize=(12, 6))
    for kmer in top_kmers:
        plt.plot(positions, enrichment[kmer], label=kmer)

    plt.xlabel("Genome position (bp)")
    plt.ylabel("K-mer frequency")
    plt.title("ORI Signal Checker: 8-mer Enrichment")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, PLOT_NAME), dpi=300)
    plt.close()

    print(f"[INFO] Plot saved → {OUT_DIR}/{PLOT_NAME}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script1_kmer_enrichment.py genomic.fa")
        sys.exit(1)
    main(sys.argv[1])

