from Bio.Align.Applications import ClustalOmegaCommandline
import sys

in_file = sys.argv[1]
out_file = in_file + '.align'
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, percentid=True, distmat_out=in_file + '.dist', distmat_full=True)
print clustalomega_cline
clustalomega_cline()
