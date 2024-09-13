

from dataclasses import dataclass


@dataclass
class SequenceMotif :
	sequence:str
	position:int



def GenerateMotifs( sequence, motif_length, iteration_start=0, iteration_end=-1, iteration_step=1 ) :
	if iteration_end == -1 : iteration_end = len(sequence) - motif_length
	for index in range( iteration_start, iteration_end, iteration_step) :
		motif_seq = sequence[index : index+motif_length]
		if len(motif_seq) != motif_length : continue
		yield SequenceMotif( position=index, sequence=motif_seq )

def GenerateMotifsPairwise( sequence, motif_length, main_motif_range_parameters = [1,-1,1], paired_motif_range_parameters = [1,24,1] ) :
	motif_1_generator = GenerateMotifs(
			sequence = sequence,
			motif_length = motif_length,
			iteration_start = main_motif_range_parameters[0],
			iteration_end = main_motif_range_parameters[1],
			iteration_step = main_motif_range_parameters[2]
		)

	for motif_1 in motif_1_generator :

		motif_1_position = motif_1.position
		motif_2_iteration_start = motif_1.position + paired_motif_range_parameters[0]
		motif_2_iteration_end = motif_1.position + paired_motif_range_parameters[1]
		motif_2_iteration_step = paired_motif_range_parameters[2]

		motif_2_generator = GenerateMotifs(
			sequence = sequence,
			motif_length = motif_length,
			iteration_start = motif_2_iteration_start,
			iteration_end = motif_2_iteration_end,
			iteration_step = motif_2_iteration_step
		)

		for motif_2 in motif_2_generator :
			yield motif_1, motif_2



if __name__ == "__main__" :

	sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

	motifs_generator =  GenerateMotifsPairwise(
		sequence = sequence, 
		motif_length=20,
	)

	for motif1, motif2 in motifs_generator :
		print(motif1, motif2)
