
from functools import lru_cache



class ComparisonScoringMatrices :

	identity_nucleotides = {
		"A" : { "A":1, "C":0, "G":0, "T":0, "U":0, "N":1 },
		"C" : { "A":0, "C":1, "G":0, "T":0, "U":0, "N":1 },
		"G" : { "A":0, "C":0, "G":1, "T":0, "U":0, "N":1 },
		"T" : { "A":0, "C":0, "G":0, "T":1, "U":1, "N":1 },
		"U" : { "A":0, "C":0, "G":0, "T":1, "U":1, "N":1 },
		"N" : { "A":1, "C":1, "G":1, "T":1, "U":1, "N":1 }
	}

	identity_aminoacids = {
		"A" : {"A":1,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"C" : {"A":0,"C":1,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"D" : {"A":0,"C":0,"D":1,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"E" : {"A":0,"C":0,"D":0,"E":1,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"F" : {"A":0,"C":0,"D":0,"E":0,"F":1,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"G" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":1,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"H" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":1,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"I" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":1,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"K" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":1,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"L" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":1,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"M" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":1,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"N" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":1,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"P" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":1,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"Q" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":1,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"R" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":1,"S":0,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"S" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":1,"T":0,"V":0,"W":0,"Y":0, "X":1},
		"T" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":1,"V":0,"W":0,"Y":0, "X":1},
		"V" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":1,"W":0,"Y":0, "X":1},
		"W" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":1,"Y":0, "X":1},
		"Y" : {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":1, "X":1},
		"X" : {"A":1,"C":1,"D":1,"E":1,"F":1,"G":1,"H":1,"I":1,"K":1,"L":1,"M":1,"N":1,"P":1,"Q":1,"R":1,"S":1,"T":1,"V":1,"W":1,"Y":1, "X":1}
	}

	identity_all_letters = dict()
	for l1 in "AZERTYUIOPQSDFGHJKLMWXCVBN" :
		identity_all_letters[l1] = dict()
		for l2 in "AZERTYUIOPQSDFGHJKLMWXCVBN" :
			identity_all_letters[l1][l2] = l1==l2


class SequenceComparison :

	def GetSequenceComparisonFunction( scoring_matrix, sequence_length ) :

		def SequenceComparisonFunction( sequence_1, sequence_2 ) :
			score_array = [ scoring_matrix[l1][l2] for l1,l2 in zip( sequence_1, sequence_2) ]
			score = sum( score_array )
			return score / sequence_length

		return SequenceComparisonFunction


class NormalizationFunctions :

	def GetHillSigmoid_Caching(inflexion, steepness, cache_maxsize=128) :

		@lru_cache(cache_maxsize)
		def HillSigmoid(value) :
			return value**steepness / (value**steepness + inflexion**steepness)

		return HillSigmoid

	def GetHillSigmoid_Uncaching(inflexion, steepness) :

		def HillSigmoid(value) :
			return value**steepness / (value**steepness + inflexion**steepness)

		return HillSigmoid




if __name__ == "__main__" :

	import numpy as np
	from matplotlib import pyplot as plt
	import time
	import random


	### TESTING LRU_CACHE ###

	print("  len        caching        uncaching            diff")
	for sequence_length in range(20,500,20) :

		ScoreNormalizer = NormalizationFunctions.GetHillSigmoid_Uncaching( 0.5, 8 )
		x = np.random.choice(sequence_length,100000)
		timer_zero = time.time()
		y = [ ScoreNormalizer(i) for i in x ]
		timer_1 = time.time()
		uncaching_timer = timer_1 - timer_zero

		ScoreNormalizer = NormalizationFunctions.GetHillSigmoid_Caching( 0.5, 8 , sequence_length)
		x = np.random.choice(sequence_length,100000)
		timer_zero = time.time()
		y = [ ScoreNormalizer(i) for i in x ]
		timer_1 = time.time()
		caching_timer = timer_1 - timer_zero

		delta = uncaching_timer - caching_timer

		print(f"{sequence_length:>5} {caching_timer:>15.8f} {uncaching_timer:>15.8f} {delta:>15.8f}")

