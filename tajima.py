#!/usr/bin/env python

import os
import subprocess
from itertools import combinations

"""
This is a protein level implementation
of Tajima's D for the universal peptide
project. It aligns peptide sequences
and removes indels first.
"""


def _run_msa(file):
    """Run clustal from executable."""
    subprocess.call(["MUSCLE", "-in", file, "-out", "MUSCLE.out"])

def _calculate_pairwise(sequences):
    """Calculate pi, number of pairwise differences."""
    for seq in sequences:
        if len(seq) != len(sequences[0]):
            raise("All sequences must have the same length.")

    numseqs = len(sequences)

    num = float(numseqs * (numseqs - 1)) / float(2)

    combos = combinations(sequences, 2)
    counts = []
    for pair in combos:
        seqA = pair[0]
        seqB = pair[1]
        count = sum(1 for a, b in zip(seqA, seqB) if a != b)
        counts.append(count)

    return(float(sum(counts)) / float(num))


def _calculate_segregating_sites(sequences):
    """Calculate S, number of segregation sites)."""
    # Assume if we're in here seqs have already been checked
    combos = combinations(sequences, 2)
    indexes = []
    for pair in combos:
        seqA = pair[0]
        seqB = pair[1]
        print seqA, seqB
        for idx, (i, j) in enumerate(zip(seqA, seqB)):
            if i != j:
                indexes.append(idx)

    indexes = list(set(indexes))

    S = len(indexes)
    n = len(sequences)

    denom = 0
    for i in range(1, n):
        denom += (float(1) / float(i))

    return float(S / denom)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def _read_sequences(file):
    sequences = []
    with open(file) as f:
        for name, seq in read_fasta(f):
            sequences.append(seq)
    return sequences

def tajimas(file):
    """Return number of pairwise mismatches against number of segreating sites."""

    _run_msa(file)

    sequences = _read_sequences("MUSCLE.out")
    
    os.remove("MUSCLE.out")

    pi = _calculate_pairwise(sequences)
    S = _calculate_segregating_sites(sequences)

    return pi - S
