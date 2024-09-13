# SequenceAutoCorrelation

small python3 script for sequence complexity analysis.
SequenceAutoCorrelation can help in the identification of low complexity regions and/or repeated regions within biological gene/protein sequences.

requires python3 libraries :
  - numpy
  - matplotlib
  - argparse
  - functools
  - dataclasses

usage: main.py [-h] -f FASTA [-m MATRIX] [-l LENGTH] [-r1 SCAN_RANGE [SCAN_RANGE ...]] [-r2 CORREL_RANGE [CORREL_RANGE ...]] [-n NORM NORM]

help :
  -f FASTA, --fasta FASTA
                        input fasta file
  -m {identity,nucleotides,aminoacids}, --matrix {identity,nucleotides,aminoacids}
                        name of the sequence comparison matrix to use for motif-motif comparison
  -l LENGTH, --length LENGTH
                        motif length used for motif-motif comparison
  -r1 SCAN_RANGE SCAN_RANGE SCAN_RANGE, --scan_range SCAN_RANGE SCAN_RANGE SCAN_RANGE
                        sequence scanning range parameters. 3 values : start, end, step
  -r2 CORREL_RANGE [CORREL_RANGE ...], --correl_range CORREL_RANGE [CORREL_RANGE ...]
                        motifs scanning range parameters for autocorrelation. 3 values : start, end, step
  -n NORM NORM, --norm NORM NORM
                        Hill Sigmoid function parameters used for normalization
