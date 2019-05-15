#!/usr/bin/env python3

import argparse

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


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--infile", "-I", type=str, required=True,
    help="FASTA file to flatten")
    parser.add_argument("--outfile", "-O", type=str, required=True,
    help="FASTA file to write to")

    return parser.parse_args()

def main():
    args = parse_args()

    outfile = open(args.outfile, "w")

    with open(args.infile) as f:
        for name, seq in read_fasta(f):
            outfile.write(name.strip() + "\n" + seq.strip() + "\n")

if __name__ == "__main__":
    main()