import sys
from Bio import pairwise2, Seq, SeqIO

def align():
	query = Seq.Seq('AGCAGAAAGTGTCAGTTGTAAAGCA')
	reference = list(SeqIO.parse(open('../rabies.fasta', 'r'), 'fasta'))[0]
	print pairwise2.align.localms(query, reference.seq[10:60], 2, -1, -1, -1)[0]

if __name__ == '__main__':
	align()
	
