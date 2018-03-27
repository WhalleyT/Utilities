#!/usr/bin/env python

import argparse

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("--motifs", "-M", dest="mot", type=str, help="A file containing motifs one would like to search for."
		                "For example, if one wanted to look for a 9-mer with a W at position 3 and V at position 5 "
		                "they would have XXWXVXXXX",
		                required=True)
	parser.add_argument("--fasta", "-F", dest="fasta", type=str, help="FASTA file containing sequences to search against.",
		                required=True)
	return parser.parse_args()


def read_fasta(fp):
	#https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def print_motif_lengths(motifs):
	lens = list(set([len(x) for x in motifs]))
	out = "".join(str(lens))
	out = out.replace("[", "")
	out = out.replace("]", "")
	print "There are motifs found of length %s in your motif file" %out


def _get_motif_amino_acids(motif):
	idxs = []
	lett = []

	for i, m in enumerate(motif):
		if m != "X":
			idxs.append(i)
			lett.append(m)
	return idxs, lett


def count_motifs(motifs, fasta):
	for item in fasta:
		names= fasta[0]
		seqs = fasta[1]

		for motif in motifs:
			motiflen = len(motif)
			indexes, letters = _get_motif_amino_acids(motif)

			for seq, name in zip(seqs, names):
				if len(seq) >= motiflen:
					upper = (len(seq) - motiflen) + 1
					for i in range(0, upper):
						peptide = seq[i:i+motiflen]
						matches = []

						for idx, lett in zip(indexes, letters):
							matches.append(lett == peptide[idx])

						if len(matches) == sum(matches):
							print "Motif '%s' found matching pattern %s in %s" %(peptide, motif, name)




def main():

	args = parse_args()

	#Read in fasta file
	names = []
	sequences = []

	with open(args.fasta) as f:
		for name, seq in read_fasta(f):
			names.append(name)
			sequences.append(seq)

	fasta_obj = (names, sequences)

	#read in motifs
	with open(args.mot) as f:
		motifs = [line.strip() for line in f]

	print_motif_lengths(motifs)
	count_motifs(motifs, fasta_obj)

if __name__ == '__main__':
	main()


