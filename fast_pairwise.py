from Bio import pairwise2
from Bio import Seq
import sys

def mismatches(record, pair):

	matches = [
		set(['A', 'T']), 
		set(['C', 'G']), 
		set(['G', 'T'])
		]

	mismatches = [
		set(['A', 'A']), 
		set(['A', 'C']),
		set(['C', 'C']), 
		set(['C', 'T']),
		set(['G', 'A']),
		set(['G', 'G']),
		set(['T', 'T'])
		]

	left = Seq.Seq(pair.left_0_seq)
	aln = pairwise2.align.localms(left, record, 2, -1, -1, -1)[0]
	score = aln[2]
	start = aln[3]
	end = aln[4]
	aligned_length = end - start
	aligned_query = aln[0][start:end]
	aligned_reference = aln[1].seq[start:end]
	primer_3prime = left[-1]
	reference_complement = str(Seq.Seq(str(aligned_reference)).complement())
	template_3prime = reference_complement[-1]
	print
	print 'record           ', record.id
	print 'length           ', len(left), start, end
	print 'score            ', score
	print 'primer           ', left
	print 'aligned_query    ', aligned_query
	print 'aligned_reference', aligned_reference
	print 'primer_3prime    ', primer_3prime
	print 'template_3prime  ', template_3prime
	for mismatch in mismatches:
		if set([primer_3prime, template_3prime]) == mismatch:
			pass
			#print '3\' mismatch left ', '5\'-%s-3\'' %left 
			#print '                 ', '3\'-%s-5\'' %reference_complement
	
#	right = Seq.Seq(pair.right_0_seq).reverse_complement()
#	aln = pairwise2.align.localms(right, record, 2, -1, -2, -1)[0]
#	start = aln[3]
#	end = start+len(right)
#	aligned_reference = aln[1].seq[start:end]
#	primer_3prime = left[-1]
#	template_3prime = str(Seq.Seq(aligned_reference[-1]).complement())
#	for mismatch in mismatches:
#		if set([primer_3prime, template_3prime]) == mismatch:
#			print '3\' mismatch right primer'
#			print '5\'-%s-3\'' %left
#			print '3\'-%s-5\'' %(Seq.Seq(str(aligned_reference)).complement())
#	
	return None
	
# G/T OK
# G/A, G/G NOT OK
