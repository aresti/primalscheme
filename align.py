import sys
import re
from Bio import pairwise2, Seq, SeqIO

def align():
	ref = list(SeqIO.parse(open('rabies.fasta', 'r'), 'fasta'))[1].seq[0:101]
	query = Seq.Seq('GCAGAAAGTGTCAGTTGTAAAGCA')
	for a in pairwise2.align.localms(query, ref, 2, -1, -1, -1):
		end = a[4]
		start = a[3]
		aln_len = end - start
		print start, end, aln_len
		print a[0][start-5:end+5]
		print a[1][start-5:end+5]

if __name__ == '__main__':
	align()
	
