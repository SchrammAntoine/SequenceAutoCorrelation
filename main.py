

import argparse
import numpy as np

from matplotlib import pyplot as plt

from SequenceAutoCorrelation import RunSequenceAutocorrelation
from SequenceSequenceComparison import ComparisonScoringMatrices



def ParseCommandLineArguments() :
	argument_parser = argparse.ArgumentParser()

	argument_parser.add_argument("-f","--fasta",
		type=str,
		required=True,
		help="input fasta file"
	) # -f		fasta
	argument_parser.add_argument("-m","--matrix",
		type=str,
		default="identity"
	) # -m		matrix
	argument_parser.add_argument("-l","--length",
		type=int,
		default = 64,
		help="motif length"
	) # -l		length
	argument_parser.add_argument("-r1","--scan_range",
		type=int,
		nargs="+",
		default=[0,-1,1]
	) # -r1	range iteration 1
	argument_parser.add_argument("-r2","--correl_range",
		type=int,
		nargs="+",
		default=[1,24,1]
	) # -r2	range iteration 2
	argument_parser.add_argument("-n","--norm",
		type=float,
		nargs=2,
		default=[0.5,8]
	) # -n		normalization

	args = argument_parser.parse_args()
	args.fasta = ExtractSequenceFromFile( args.fasta )
	args.matrix = GetComparisonMatrix( args.matrix )
	if   args.scan_range[1] < 0 : args.scan_range[1] = len(args.fasta) +1 - args.scan_range[1]
	if args.correl_range[1] < 0 : args.correl_range[1] = args.scan_range[1] +1 - args.correl_range[1]
	return args

def PrintParameters(args) :
	print("-- Sequence Analysis by AutoCorrelation --")
	print("-- Author : Antoine Schramm --")
	print()
	print("Input Sequence : ")
	for i in range(1+len(args.fasta)//100) : print(args.fasta[i*100:i*100+100])
	print()
	print(f"Auto-correlation will be carried out on motifs from position {args.scan_range[0]} to {args.scan_range[1]} with step of {args.scan_range[2]}")
	print(f"Motifs will be compared with downstream motifs with offset ranging from {args.correl_range[0]} to {args.correl_range[1]} with step {args.correl_range[2]}")
	print(f"Motifs used for Motif-Motif comparison are {args.length} letters long")
	print(f"Motif-Motif comparison will be carried out by sequence superposition using the following matrix")
	print()
	print(" " + "".join( list(args.matrix.keys()) ) )
	for l1 in args.matrix : 
		string = l1
		for l2 in args.matrix : 
			value = args.matrix[l1][l2]
			if value : string += "█"
			else : string += "░"
		print(string)
	print()
	print(f"Scores will be normalized using a Hill Sigmoid Function parametrized with inflexion at {args.norm[0]} and steepness power of {args.norm[1]}")
	print()


def ExtractSequenceFromFile(file_name ) :
	f = open(file_name,"r")
	sequence = ""
	for line in f :
		if line.startswith(">") : sequence = ""
		else :
			sequence += line.strip().upper()
	return sequence

def GetComparisonMatrix(matrix_name ) :
	if matrix_name == "identity" : return ComparisonScoringMatrices.identity_all_letters
	if matrix_name == "nucleotides" : return ComparisonScoringMatrices.identity_nucleotides
	if matrix_name == "aminoacids" : return ComparisonScoringMatrices.identity_aminoacids

def PlotResults( dataframe ) :
	dataframe = np.transpose(dataframe)
	plt.matshow( dataframe , aspect="auto" )
	plt.show()

def main() :
	args = ParseCommandLineArguments()
	PrintParameters(args)
	results = RunSequenceAutocorrelation(
		sequence = args.fasta,
		motif_length = args.length,
		scanning_range = args.scan_range,
		correlation_range = args.correl_range,
		scoring_matrix = args.matrix,
		normalization_parameters = args.norm
	)
	PlotResults(results)


if __name__ == "__main__" :
	main()