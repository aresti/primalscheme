import sys
from Bio import SeqIO
import primer3
from pprint import pprint
import settings
import interactions
import pairwise
from primer_pair import region_object
#from openpyxl import load_workbook
#from openpyxl.worksheet.datavalidation import DataValidation
from operator import itemgetter

def run(args):
	#Use first record for primer picking
	records = list(SeqIO.parse(open(args.g, 'r'), 'fasta'))
	seq = str(records[0].seq)
	step = args.length - 2*(args.overlap)
	seq_params = {
		     'SEQUENCE_ID': 'blah',
		     'SEQUENCE_TEMPLATE': seq,
		     'SEQUENCE_INCLUDED_REGION': [0, len(seq)],
		     }
	return records, seq, step, seq_params 

def multiplex(args, parser=None):
	records, seq, step, seq_params = run(args)
	outer_params = settings.outer_pair(args)
	outer_pairs = []
	print 'Running primer3'
	for i, region in enumerate(range(args.length/2, len(seq), step)):
		sys.stdout.write('.')
		sys.stdout.flush()
		seq_params['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [region-(args.length/2), args.overlap, region+(args.length/2)-args.overlap, args.overlap]
		if region > len(seq)-step:
			seq_params['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [region-200, args.overlap, len(seq)-50, 50]
			output = primer3.bindings.designPrimers(seq_params, outer_params)
		else:
			output = primer3.bindings.designPrimers(seq_params, outer_params)
			try:
				output['PRIMER_LEFT_0']
				output['PRIMER_RIGHT_0']
			except KeyError:
				seq_params['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [region-(args.length/2), args.overlap+25, region+(args.length/2)-(args.overlap+25), args.overlap+25]
				output = primer3.bindings.designPrimers(seq_params, outer_params)

		#pprint(output, width=1)
		region = region_object(args.length, i+1, output)
		#pprint(vars(region), width=1)
		outer_pairs.append(region)

	#check for mismatches
	print
	print 'Checking for 3\' mismatches'
	totals = []
	for region_pairs in outer_pairs:
		sub_totals = []
		for pair in region_pairs.pairs:
			scores = 0
			for record in records:
				pairwise.fast_pairwise(record, pair)
				scores += (pair.left_aln_score + pair.right_aln_score)
				if (pair.left_3prime_mm == True or pair.right_3prime_mm == True):
					print 'mismatch'
					scores = 0 #should always be 0?
			pair.sub_total = scores
			#sub_totals.append((pair, scores))
		#region_pairs.totals = sub_totals
		#totals.append(sub_totals)
	#pprint(vars(outer_pairs[0].pairs[0]), width=1)

	#check for interactions
	#print interactions.interactions(args, outer)

	result = []
	for i, region_pairs in enumerate(outer_pairs):
		pool = str(['1' if i%2==0 else '2'][0])
		region_pairs.pool = pool
		result.append(sorted(region_pairs.pairs, key=lambda x: x.sub_total, reverse=True)[0])
	return result
	

def nested_multiplex(parser, args):
        seq, seq_params = run(args)
	outer = []
	inner = []
