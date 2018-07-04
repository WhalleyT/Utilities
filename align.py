import argparse
import sys
import progressbar
import multiprocessing

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import Bio.SubsMat.MatrixInfo as mat_info

from joblib import Parallel, delayed

#Parse arguments from command line.
def parse_args():
	parser = argparse.ArgumentParser("Parser for proteome aligment")

	parser.add_argument("-R", "--reference", help="Reference file you want to align your peptides \
		                against. Must be a fasta file containing amino acids", required=True, type=str,
		                dest="reffile")
	parser.add_argument("-Q", "--query", help="peptide file you want to find alignemnts for \
		                Must be a fasta file containing amino acids", required=True, type=str,
		                dest="pepfile")
	parser.add_argument("-M", "--matrix", help="Protein alignment matrix", default="pam30",
		                dest="matrix", required=False, type=str)
	parser.add_argument("-L", "--length", help="Peptide length of interest", default=9,
	                    dest="peplen", required=False, type=int)
	parser.add_argument("-O", "--outfile", help="Output file", default="alignments.csv",
		                dest="outfile", required=False, type=str)
	parser.add_argument("-N", "--nthreads", help="Number of parallel threads. Pick an integer for a \
		                specific amount, leave blank or 'max' for the maximum,", default="max", type=str, 
		                dest="nthreads")

	return parser.parse_args()


#generator for reading fasta
def _read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


#use generator to pass to dictionary
def fasta_dict(infile):
	d = {}

	with open(infile) as f:
		for name, seq in _read_fasta(f):
			d[name] = seq

	return d


def check_matrix(mat):

	available = mat_info.available_matrices
	if mat not in available:
		print "Could not find correct alignment matrix." 
		print "Please check from the following list to make sure you have selected correctly:"

		for  i in available:
			print "\t-" + i

		sys.exit()
	else:
		return getattr(mat_info, mat)


def alignment(reference, viral, peplen, matrix):
	out = []
	widgets = [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()]
	pbar = progressbar.ProgressBar(widgets=widgets)

	for viral_name in pbar(viral):
		viral_protein = viral[viral_name]
		v_endrange = len(viral_protein) - (peplen - 1)

		for reference_name in reference:
			reference_protein = reference[reference_name]
			r_endrange = len(reference_protein) - (peplen - 1)


			if viral_name != reference_name:

				for i in range(0, v_endrange):
					viral_pep = viral_protein[i:i+peplen]
					viral_scores = []
					matches = 0

					for j in range(0, r_endrange):
						ref_pep = reference_protein[j:j+peplen]
						#put endwsith R or K in pos13 (end)

						score = 0

						if ref_pep == viral_pep:
							matches += 1

						for vAA, rAA in zip(viral_pep, ref_pep):
							AA = (vAA, rAA)
							rev_AA = (rAA, vAA)

							if AA in matrix:
								score += matrix[AA]
							elif rev_AA in matrix:
								score += matrix[rev_AA]
							else:
								score += 0

						viral_scores.append(float(score))

					if len(viral_scores) != 0 and sum(viral_scores) != 0:
						mean = sum(viral_scores) / float(len(viral_scores))
					else:
						mean = 0.0
					out.append([viral_name, viral_pep, mean, matches])



	return pd.DataFrame(out, columns = ["viral_name", "viral_peptide", "mean_score", "matches"])



def parallel_alignment(reference, viral, peplen, matrix, viral_name):
	out = []

	viral_protein = viral[viral_name]
	v_endrange = len(viral_protein) - (peplen - 1)

	for reference_name in reference:
		reference_protein = reference[reference_name]
		r_endrange = len(reference_protein) - (peplen - 1)


		if viral_name != reference_name:

			for i in range(0, v_endrange):
				viral_pep = viral_protein[i:i+peplen]
				viral_scores = []
				matches = 0

				for j in range(0, r_endrange):
					ref_pep = reference_protein[j:j+peplen]
					#put endwsith R or K in pos13 (end)

					score = 0

					if ref_pep == viral_pep:
						matches += 1

					for vAA, rAA in zip(viral_pep, ref_pep):
						AA = (vAA, rAA)
						rev_AA = (rAA, vAA)

						if AA in matrix:
							score += matrix[AA]
						elif rev_AA in matrix:
							score += matrix[rev_AA]
						else:
							score += 0

					viral_scores.append(float(score))

				if len(viral_scores) != 0 and sum(viral_scores) != 0:
					mean = sum(viral_scores) / float(len(viral_scores))
				else:
					mean = 0.0

				out.append([viral_name, viral_pep, mean, matches])

	return pd.DataFrame(out, columns = ["viral_name", "viral_peptide", "mean_score", "matches"])


def mpl_clear():
	plt.cla()
	plt.clf()
	plt.close()


def main():
	args = parse_args()

	if args.nthreads == "max":
		num_cores = multiprocessing.cpu_count()
	else:
		num_cores = int(args.nthreads)

	reference = fasta_dict(args.reffile)
	virus = fasta_dict(args.pepfile)
	alignment_matrix = check_matrix(args.matrix)

	if num_cores == 1:
		print "Starting alignment in serial"
		df = alignment(reference, virus, args.peplen, alignment_matrix)
		print "Saving dataframe output to " + args.outfile
		df.to_csv(args.outfile)
	else:
		print "Starting alignment in parallel"
		list_of_df = Parallel(n_jobs=num_cores, verbose=5)(delayed(parallel_alignment)
			(reference, virus, args.peplen, alignment_matrix, viral_name)
			for viral_name in virus)
		df = pd.concat(list_of_df)
		print "Saving dataframe output to " + args.outfile
		df.to_csv(args.outfile)



if __name__ == '__main__':
	main()