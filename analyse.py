import sys
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO
import primer3
import re


fwd_primers = []
rev_primers = []
for line in open(sys.argv[1], 'r'):
	if line.startswith('#'):
		continue
	cols = line.strip().split(',')
	if cols[4] == '.':
		continue
	fwd_primers.append(cols[4])
	rev_primers.append(cols[5])

print 'primer, length, Tm_primer3, GC, homopolymer, ambigious'
for primer in fwd_primers + rev_primers:
	print primer, len(primer), '%0.2f' %(primer3.calcTm(primer, dv_conc=2.0, mv_conc=50, dntp_conc=0.8)), '%0.2f' %(GC(primer)), max([len(m.group()) for m in re.finditer(r'([ACGT])\1{1,}', primer)]), len(primer)-(primer.count('A') + primer.count('C') + primer.count('G') + primer.count('T')), '%0.2f' %(primer3.calcHairpin(primer, dv_conc=2.0, mv_conc=50, dntp_conc=0.8).tm)
