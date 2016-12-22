from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio import pairwise2
import primer3
import re

class Primer():

	def __init__(self, n, s, e, l, p, gc, tm, reference):
		self.name = n
		self.start = s
		self.end = e
		self.length = l
		self.seq = p
		self.gc1 = gc
		self.Tm1 = tm
		self.reference = reference
		self.calcTm2()
		self.calcHP()
		self.homopolymer()
		self.gc2()
		self.ambiguous()
		self.align()

	def calcTm2(self):
		self.Tm2 = primer3.calcTm(self.seq, dv_conc=2.0, mv_conc=50, dntp_conc=0.8)

	def calcHP(self):
		self.HP = primer3.calcHairpin(self.seq, dv_conc=2.0, mv_conc=50, dntp_conc=0.8).tm

	def homopolymer(self):
		self.homopolymer = max([len(m.group()) for m in re.finditer(r'([ACGT])\1{1,}', self.seq)])
	def gc2(self):
		self.gc2 = GC(self.seq)

	def ambiguous(self):
		self.ambiguous_bases = self.length-(self.seq.count('A') + self.seq.count('C') + self.seq.count('G') + self.seq.count('T'))

	def align(self):
		#query = (self.seq if self.strand == +1 else Seq(self.seq).reverse_complement())
		query = Seq(self.seq)
		revcomp = query.reverse_complement()
		aln_query = pairwise2.align.localms(query, self.reference.seq, 2, -1, -2, -1)[0]
		aln_revcomp = pairwise2.align.localms(revcomp, self.reference.seq, 2, -1, -2, -1)[0]
		if aln_query[2] > aln_revcomp[2]:
			self.aln_strand = +1
			self.score = aln_query[2]
			self.aln_start, self.aln_end = aln_query[3], aln_query[3]+self.length
		else:
			self.aln_strand = -1
			self.score = aln_revcomp[2]
			self.aln_start, self.aln_end = aln_revcomp[4]-1, aln_revcomp[4]-self.length-1
