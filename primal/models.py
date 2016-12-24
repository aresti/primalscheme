import settings
import re
import sys

from Bio import pairwise2
from Bio import Seq
from pprint import pprint


class Explain():
	def __init__(self, string):
		self.considered = 0
		self.GC_content_failed = 0
		self.low_tm = 0
		self.high_tm = 0
		self.not_in_any_ok_left_region = 0
		self.not_in_any_ok_right_region = 0
		self.ok = 0
		self.explain(string)

	def explain(self, string):
		fields = ['considered', 'GC content failed', 'low tm', 'high tm', \
			  'not in any ok left region', 'not in any ok right region', 'ok']
		for field in fields:
			pattern = field.replace(' ', '\s') + '\s(\d+)'
			matches = re.search(pattern, string)
			if matches:
				setattr(self, field.replace('\s','_'), int(matches.group(1)))


class Primitive_primer():
	def __init__(self, direction, name, seq):
		self.direction = direction
		self.name = name
		self.seq = seq
		self.start = 0
		self.end = 10440 #For zika
		self.alignment = primal_models.Alignment(self, ref)


class PrimerPair():
	def __init__(self, left, right):
		self.left = left
		self.right = right
		self.total = left.sub_total + right.sub_total


class Primer():
	def __init__(self, scheme, region, output, n, direction, references):
		self.direction = direction
		self.name = '%i_%i_%s_%i' %(scheme, region, direction, n)
		self.length = int(output['PRIMER_%s_%i' %(direction, n)][1])
		if self.direction == 'LEFT':
			self.start = int(output['PRIMER_%s_%i' %(direction, n)][0])
			self.end = int(self.start + self.length)
		elif self.direction == 'RIGHT':
			self.start = int(output['PRIMER_%s_%i' %(direction, n)][0] + 1)
			self.end = int(self.start - self.length)
		self.seq = str(output['PRIMER_%s_%i_SEQUENCE' %(direction, n)])
		self.gc = float(output['PRIMER_%s_%i_GC_PERCENT' %(direction, n)])
		self.tm = float(output['PRIMER_%s_%i_TM' %(direction, n)])
		self.sub_total = 0
		self.alignments = []
		#pprint(vars(self), width=1)

		for ref in references:
			alignment = Alignment(self, ref)
			#pprint(vars(alignment), width=1)
			self.alignments.append(alignment)
			self.sub_total += alignment.score


class Region():
	def __init__(self, scheme, region_num, output, references, position):
		self.scheme = scheme
		self.region_num = region_num
		self.pool = '2' if self.region_num % 2 == 0 else '1'
		self.pairs = []
		#should pass number of candidates to function
		for cand_num in range(5):
			lenkey = 'PRIMER_LEFT_%s' %cand_num
			if lenkey not in output:
				break
			output[lenkey] = (output[lenkey][0] + position, output[lenkey][1])
			#print output[lenkey]
			left = Primer(scheme, region_num, output, cand_num, 'LEFT', references)
			right = Primer(scheme, region_num, output, cand_num, 'RIGHT', references)
			self.pairs.append(PrimerPair(left, right))
		self.pairs.sort(key=lambda x: x.total, reverse=True)

class Alignment():
	def __init__(self, primer, reference):
		self.start = 0
		self.end = 0
		self.length = 0
		self.score = 0
		self.aln_query = ''
		self.aln_ref = ''
		self.aln_ref_comp = ''
		self.template_3prime = ''
		self.primer_3prime = ''
		self.mm_3prime = False
		self.fast_pairwise(primer, reference)

	def fast_pairwise(self, primer, ref):
		self.primer = primer.seq
		self.primer_length = len(self.primer)
		if primer.direction == 'LEFT':
			search_start = primer.start-50 if primer.start >= 0 else 0
			search_end = primer.end+50 if primer.end+50 <= len(ref) else len(ref)
			alns = pairwise2.align.localms(primer.seq, ref.seq[search_start:search_end], 2, -1, -1, -1, penalize_end_gaps=True)
		elif primer.direction == 'RIGHT':
			search_start = primer.start-50
			search_end = primer.end+50 if primer.end+50 <= len(ref) else len(ref)
			alns = pairwise2.align.localms(primer.seq, ref.seq[search_start:search_end].reverse_complement(), 2, -1, -1, -1, penalize_end_gaps=True)
		if alns:
			aln = alns[0]
			p = re.compile('(-*)([ACGTN][ACGTN\-]*[ACGTN])(-*)')
			m = list(re.finditer(p, str(aln[0])))[0]
			self.score = aln[2]
			self.start = (search_start + m.span(2)[0] if primer.direction == 'LEFT' else search_end - m.span(2)[0])
			self.end = (search_start + m.span(2)[1] if primer.direction == 'LEFT' else search_end - m.span(2)[1])
			self.length = (self.end - self.start if primer.direction == 'LEFT' else self.start - self.end)
			#print self.start, self.end, self.length, self.score
			self.aln_query = aln[0][m.span(2)[0]:m.span(2)[1]]
			self.aln_ref = aln[1][m.span(2)[0]:m.span(2)[1]]
			self.primer_3prime = self.aln_query[-1]
			self.aln_ref_comp = Seq.Seq(str(self.aln_ref)).complement()
			self.template_3prime = self.aln_ref_comp[-1]
			self.mm_3prime = False
			#pprint(vars(self), width=1)

			#check this is left
			for mismatch in settings.MISMATCHES:
				if set([self.primer_3prime, self.template_3prime]) == mismatch:
					print '3\' mismatch left ', '5\' %s 3\'' %self.aln_query
					print '                 ', '3\' %s 5\'' %self.aln_ref_comp
					self.mm_3prime = True
					self.score = 0
		else:
			self.score = 0
