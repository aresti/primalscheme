import settings
import re
import sys

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Seq
from pprint import pprint


class Explain(object):
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


class Primer(object):
	def __init__(self, direction, name, seq):
		# TODO: Validate direction is LEFT or RIGHT
		self.direction = direction
		self.name = name
		self.seq = seq

	@property
	def length(self):
		return len(self.seq)


class CandidatePrimer(Primer):
	def __init__(self, direction, name, seq, start, gc, tm, references):
		super(CandidatePrimer, self).__init__(direction, name, seq)
		self.start = start
		self.gc = gc
		self.tm = tm

		self.sub_total = 0
		self.alignments = []

		for ref in references:
			alignment = Alignment(self, ref)
			self.alignments.append(alignment)
			self.sub_total += alignment.score

	@property
	def end(self):
		if self.direction == 'LEFT':
			return self.start + self.length
		else:
			return self.start - self.length


class CandidatePrimerPair(object):
	def __init__(self, left, right):
		self.left = left
		self.right = right
		self.total = left.sub_total + right.sub_total

	@property
	def product_length(self):
		return self.right.start - self.left.start + 1


class Region(object):
	def __init__(self, prefix, region_num, start_limits, primer3_output, references):
		self.region_num = region_num
		self.pool = '2' if self.region_num % 2 == 0 else '1'
		self.candidate_pairs = []

		for cand_num in range(5):
			lenkey = 'PRIMER_LEFT_%s' % (cand_num)
			left_name = '%s_%i_%s' % (prefix, region_num, 'LEFT')
			right_name = '%s_%i_%s' % (prefix, region_num, 'RIGHT')
			if lenkey not in primer3_output:
				print 'Only %s candidate primers' % (cand_num)
				break

			left_seq = str(primer3_output['PRIMER_LEFT_%i_SEQUENCE' % (cand_num)])
			right_seq = str(primer3_output['PRIMER_RIGHT_%i_SEQUENCE' % (cand_num)])

			left_start = int(primer3_output['PRIMER_LEFT_%i' % (cand_num)][0] + start_limits[0])
			right_start = int(primer3_output['PRIMER_RIGHT_%i' % (cand_num)][0] + start_limits[0] + 1)

			left_gc = float(primer3_output['PRIMER_LEFT_%i_GC_PERCENT' % (cand_num)])
			right_gc = float(primer3_output['PRIMER_RIGHT_%i_GC_PERCENT' % (cand_num)])

			left_tm = float(primer3_output['PRIMER_LEFT_%i_TM' % (cand_num)])
			right_tm = float(primer3_output['PRIMER_RIGHT_%i_TM' % (cand_num)])

			left = CandidatePrimer('LEFT', left_name, left_seq, left_start, left_gc, left_tm, references)
			right = CandidatePrimer('RIGHT', right_name, right_seq, right_start, right_gc, right_tm, references)

			self.candidate_pairs.append(CandidatePrimerPair(left, right))
		self.candidate_pairs.sort(key=lambda x: x.total, reverse=True)


class Alignment():
	def __init__(self, primer, ref):
		#get alignments
		if primer.direction == 'LEFT':
			search_start = primer.start - 50 if primer.start > 50 else 0
			search_end = primer.end + 50 if primer.end + 50 <= len(ref) else len(ref)
			alns = pairwise2.align.globalms(str(primer.seq), str(ref.seq[search_start:search_end]), 2, -1, -1, -1, penalize_end_gaps=False, one_alignment_only=True)
		elif primer.direction == 'RIGHT':
			search_start = primer.end - 50
			search_end = primer.start + 50 if primer.start + 50 <= len(ref) else len(ref)
			alns = pairwise2.align.globalms(primer.seq, ref.seq[search_start:search_end].reverse_complement(), 2, -1, -1, -1, penalize_end_gaps=False, one_alignment_only=True)
		if alns:
			#print(format_alignment(*alns[0]))
			aln = alns[0]
			p = re.compile('(-*)([ACGTN][ACGTN\-]*[ACGTN])(-*)')
			m = list(re.finditer(p, str(aln[0])))[0]
			self.score = aln[2]

			if primer.direction == 'LEFT':
				self.start = search_start + m.span(2)[0]
				self.end = search_start + m.span(2)[1]
				self.length = self.end - self.start
			else:
				self.start = search_end - m.span(2)[0]
				self.end = search_end - m.span(2)[1]
				self.length = self.start - self.end

			self.aln_query = aln[0][m.span(2)[0]:m.span(2)[1]]
			self.aln_ref = aln[1][m.span(2)[0]:m.span(2)[1]]
			self.aln_ref_comp = Seq.Seq(str(self.aln_ref)).complement()
			self.cigar = ''
			self.mm_3prime = False

			#make cigar
			for a, b in zip(self.aln_query, self.aln_ref):
				if a == '-' or b == '-':
					self.cigar += ' '
					continue
				if a != b:
					self.cigar += '*'
					continue
				else:
					self.cigar += '|'

			#check 3' mismatches
			if set([self.aln_query[-1], self.aln_ref_comp[-1]]) in settings.MISMATCHES:
				self.mm_3prime = True
				self.score = 0

				#pprint(vars(self), width=1)
				#print
				#print '3\' mismatch'
				#short_primer = primer.name[:50] if len(primer.name) > 50 else primer.name
				#short_ref = ref.id[:50] if len(ref.id) > 50 else ref.id
				#print "{: <50}".format(short_primer), '5\'-%s-3\'' %self.aln_query
				#print "{: <50}".format(''), '   %s' %self.cigar
				#print "{: <50}".format(short_ref), '3\'-%s-5\'' %self.aln_ref_comp
		else:
			self.score = 0
