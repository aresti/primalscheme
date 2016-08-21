from itertools import product
from Bio import pairwise2
from Bio import Seq

def pairwise(pair):
	primer1 = Seq.Seq(pair[0].left_0_seq)
	primer2 = Seq.Seq(pair[1].right_0_seq)
	aln = pairwise2.align.localms(primer1, primer2, 2, -1, -2, -1)[0]
	if aln[2] >= 20:
		print pair[0].left_0_name, pair[1].right_0_name, aln
	
def interactions(args, scheme):
	for each in product(scheme, repeat=2):
		if each[0] == each[1]:
			continue
		else:
			pairwise(each)
