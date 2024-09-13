
import numpy as np

from SequenceMotifGenerators import GenerateMotifsPairwise
from SequenceSequenceComparison import NormalizationFunctions, SequenceComparison

def InstanciateResultTable(scan_range, correl_range ) :
	ncols = (scan_range[1] - scan_range[0]) // scan_range[2]
	nrows = (correl_range[1] - correl_range[0]) // correl_range[2]
	return np.zeros( (ncols, nrows) )


def RunSequenceAutocorrelation(sequence, motif_length, scanning_range, correlation_range, scoring_matrix, normalization_parameters ) :

	motifs_generator = GenerateMotifsPairwise(
		sequence = sequence, 
		motif_length = motif_length, 
		main_motif_range_parameters = scanning_range, 
		paired_motif_range_parameters = correlation_range
	)

	MotifMotifComparisonFunc = SequenceComparison.GetSequenceComparisonFunction( 
		scoring_matrix = scoring_matrix, 
		sequence_length = motif_length
	)

	ScoringNormalizerFunc = NormalizationFunctions.GetHillSigmoid_Caching( 
		inflexion = normalization_parameters[0], 
		steepness = normalization_parameters[1],
		cache_maxsize = motif_length
	)

	result = InstanciateResultTable(scanning_range, correlation_range)
	for motif_1, motif_2 in motifs_generator :
		comparison_score_raw  = MotifMotifComparisonFunc( motif_1.sequence, motif_2.sequence)
		comparison_score_norm = ScoringNormalizerFunc(comparison_score_raw)
		index_1 = (motif_1.position - scanning_range[0])//scanning_range[2]
		index_2 = motif_2.position - motif_1.position - correlation_range[0]
		result[index_1,index_2] = comparison_score_norm
	return result
