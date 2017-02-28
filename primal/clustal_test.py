from Bio.Align.Applications import ClustalOmegaCommandline
import sys

in_file = '../4_HAV-genomes-multi.fa'
out_file = 'aligned.fasta'
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, percentid=True, distmat_out="dists.txt", distmat_full=True)
print clustalomega_cline
clustalomega_cline()
