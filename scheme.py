import sys
from Bio import SeqIO
import primer3
from pprint import pprint
import settings
import interactions
import pairwise
from primal_models import Region, Explain
from operator import itemgetter

def run(args):
	#Use first record for primer picking
	references = list(SeqIO.parse(open(args.g, 'r'), 'fasta'))
	seq = str(references[0].seq)
	regions, remain = divmod(len(seq)-(2*args.overlap), args.length-(2*args.overlap))
	#substract the overlaps from the ends then divide the rest by the non-overlapping region 
	step = args.length-(2*args.overlap) + int(round(remain/float(regions)))
	#divide the remainder across the regions
	print 'Reference size %i, step size %i, regions %i, remain %i' %(len(seq), step, regions, remain)
	seq_params = {
		     'SEQUENCE_ID': 'blah',
		     'SEQUENCE_TEMPLATE': seq,
		     'SEQUENCE_INCLUDED_REGION': [0, len(seq)],
		     }
	return references, seq, step, seq_params, regions 

def multiplex(args, parser=None):
	references, seq, step, seq_params, regions = run(args)
	outer_params = settings.outer_pair(args)
	result = []
	key = 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'
	print 'Running primer3'

	for i, region in enumerate(range(0, len(seq), step)[:-1]):
		print 'Region %i, position %i' %(i+1, region)
		if i == 0:
			max_left_end = region
		else:
			max_left_end = result[-1].pairs[0].right.end-(args.overlap+20)
		left_start = max_left_end
		left_end = args.overlap
		right_start = region+args.length-args.overlap
		right_end = args.overlap
		if i+1 == regions:
			right_end = 25
			right_start = len(seq)-right_end
		for j in range(5):
			seq_params[key] = [left_start, left_end, right_start, right_end]
			print seq_params[key]
			output = primer3.bindings.designPrimers(seq_params, outer_params)
			#pprint(output, width=1)
			left_ok = Explain(output['PRIMER_LEFT_EXPLAIN']).ok
			right_ok = Explain(output['PRIMER_RIGHT_EXPLAIN']).ok
			if left_ok < 3:
				print 'Found %i left and %i right primers on iteration %i' %(left_ok, right_ok, j+1)
				left_start -= 10
				left_end += 10
			elif right_ok < 3:
				print 'Found %i left and %i right primers on iteration %i' %(left_ok, right_ok, j+1)
				if i+1 == regions:
					#Work back from right for last region
					right_end += 10
					right_start = len(seq)-right_end
				else:
					right_start -= 10
					right_end += 10 
			else:
				print 'Found %i left and %i right primers on iteration %i' %(left_ok, right_ok, j+1)
				break
		

		region = Region(args.length, i+1, output, references)
		#print vars(region.pairs[0].left)
		result.append(region)
		#pprint(output, width=1)
		#pprint(vars(region), width=1)
		#pprint(vars(region.pairs[0]), width=1)

	return result

	
def nested_multiplex(parser, args):
        seq, seq_params = run(args)
	outer = []
	inner = []
