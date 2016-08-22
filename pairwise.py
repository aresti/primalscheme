from Bio import pairwise2
from Bio import Seq
import sys

def fast_pairwise(record, pairs, number):

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

	left_name = getattr(pairs,'left_%i_name' %number)
	left_seq = Seq.Seq(getattr(pairs,'left_%i_seq' %number))
	left_start = getattr(pairs,'left_%i_start' %number)
	left_end = getattr(pairs,'left_%i_end' %number)
	search_start = (left_start-50 if left_start-50 >= 0 else 0)
	search_end = (left_end+50 if left_end+50 <= len(record) else len(record))
	aln = pairwise2.align.localms(left_seq, record.seq[search_start:search_end], 2, -1, -1, -1)[0]
	left_aln_score = aln[2]
	aln_start = aln[3]
	aln_end = aln[4]
	aln_length = aln_end - aln_start
	if aln_length < len(left_seq): #in some cases this fails
		print 'unexpected alignment length'
		aln_end += 1
	aln_query = aln[0][aln_start:aln_end]
	aln_reference = aln[1][aln_start:aln_end]
	primer_3prime = aln_query[-1]
	reference_complement = str(Seq.Seq(str(aln_reference)).complement())
	template_3prime = reference_complement[-1]
	setattr(pairs,'left_%i_3prime_mm' %number, False)
	setattr(pairs,'left_%i_aln_score' %number, left_aln_score)

	print 'primer           ', left_name
	print 'sequence         ', left_seq
	print 'reference        ', record.id
	print 'length           ', len(left_seq)
	print 'start            ', search_start + aln_start
	print 'end              ', search_start + aln_end
	print 'aligned_length   ', aln_length
	print 'score            ', left_aln_score
	print 'aligned_query    ', aln_query
	print 'aligned_reference', reference_complement
	print '3\' pairing       ', '%s/%s' %(primer_3prime, template_3prime)

	for mismatch in mismatches:
		if set([primer_3prime, template_3prime]) == mismatch:
			print '3\' mismatch left ', '5\' %s 3\'' %aln_query 
			print '                 ', '3\' %s 5\'' %reference_complement
			setattr(pairs,'left_%i_3prime_mm' %number, True)
			
	print

	right_name = getattr(pairs,'right_%i_name' %number)
	right_seq = Seq.Seq(getattr(pairs,'right_%i_seq' %number))
	right_start = getattr(pairs,'right_%i_start' %number)
	right_end = getattr(pairs,'right_%i_end' %number)
	search_start = (right_start-50 if right_start-50 >= 0 else 0)
	search_end = (right_end+50 if right_end+50 <= len(record) else len(record))
	aln = pairwise2.align.localms(right_seq, record.seq[search_start:search_end].reverse_complement(), 2, -1, -1, -1)[0]
	right_aln_score = aln[2] #probably wrong
	aln_start = aln[3]
	aln_end = aln[4]
	aln_length = aln_end - aln_start
	if aln_length < len(right_seq):
		print 'unexpected alignment length'
		aln_start -= 1
	aln_query = aln[0][aln_start:aln_end]
	aln_reference = aln[1][aln_start:aln_end]
	primer_3prime = aln_query[-1]
	reference_complement = str(Seq.Seq(str(aln_reference)).complement())
	template_3prime = reference_complement[-1]
	setattr(pairs,'right_%i_3prime_mm' %number, False)
	setattr(pairs,'right_%i_aln_score' %number, right_aln_score)

	print 'primer           ', right_name
	print 'sequence         ', right_seq
	print 'reference        ', record.id
	print 'length           ', len(right_seq)
	print 'start            ', search_start + aln_start
	print 'end              ', search_start + aln_end
	print 'aligned_length   ', aln_length
	print 'score            ', right_aln_score
	print 'aligned_query    ', aln_query
	print 'aligned_reference', reference_complement
	print '3\' pairing       ', '%s/%s' %(primer_3prime, template_3prime)

	for mismatch in mismatches:
		if set([primer_3prime, template_3prime]) == mismatch:
			print '3\' mismatch right ', '5\' %s 3\'' %aln_query 
			print '                 ', '3\' %s 5\'' %reference_complement
			setattr(pairs,'right_%i_3prime_mm' %number, True)
	print

	return None

# G/T OK
# G/A, G/G NOT OK

def slow_pairwise(record, pair):

	print
	print 'record           ', record.id

	left = Seq.Seq(pair.left_0_seq)
	aln = pairwise2.align.localms(left, record, 2, -1, -1, -1)[0]
	left_start = aln[3]
	left_end = aln[4]
	left_score = aln[2]

	right = Seq.Seq(pair.right_0_seq)
	aln = pairwise2.align.localms(right, record, 2, -1, -1, -1)[0]
	right_start = aln[3]
	right_end = aln[4]
	right_score = aln[2]

	print 'primer           ', left
	print 'start            ', left_start
	print 'end              ', left_end
	print 'score            ', left_score	
	print 'primer           ', right
	print 'start            ', right_start
	print 'end              ', right_end
	print 'score            ', right_score

	return None
