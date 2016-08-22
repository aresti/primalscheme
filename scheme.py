import sys
from Bio import SeqIO
import primer3
from pprint import pprint
import settings
import interactions
import pairwise
from primerpair import primerpair
#from openpyxl import load_workbook
#from openpyxl.worksheet.datavalidation import DataValidation
#from operator import itemgetter

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

def multiplex(parser, args):
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
		pair = primerpair(args.length, i+1, output)
		#pprint(vars(pair), width=1)
		outer_pairs.append(pair)

	#check for mismatches
	print
	print 'Checking for 3\' mismatches'
	pair_scores = []
	totals = []
	for region_pairs in outer_pairs:
		sub_totals = []
		for pair in range(4):
			scores = 0
			for record in records:
				pairwise.fast_pairwise(record, region_pairs, pair)
				left_score = getattr(region_pairs,'left_%i_aln_score' %pair) 
				right_score = getattr(region_pairs,'right_%i_aln_score' %pair)
				left_3prime_mm = getattr(region_pairs,'left_%i_3prime_mm' %pair)
				right_3prime_mm = getattr(region_pairs,'right_%i_3prime_mm' %pair)
				scores += (left_score + right_score)
				if (left_3prime_mm == True or right_3prime_mm == True):
					print 'mismatch'
					scores = 0 #may not always be 0
			sub_totals.append((pair, scores))
		totals.append(sub_totals)

	picked = []
	for each in totals:
		ordered = sorted(each, key=itemgetter(1), reverse=True)
		print ordered
		picked.append(ordered[0])

	#check for interactions
	#print interactions.interactions(args, outer)

	#return results as list
	result = []
	for i, each in enumerate(outer_pairs):
		pool = str(['1' if i%2==0 else '2'][0])
		result.append(getattr(each, 'left_%s_name' %picked[i][0]), '\t', getattr(each, 'left_%s_seq' %picked[i][0]), '\t', pool)
		result.append(getattr(each, 'right_%s_name' %picked[i][0]), '\t', getattr(each, 'right_%s_seq' %picked[i][0]), '\t', pool)
	return result

def nested_multiplex(parser, args):
        seq, seq_params = run(args)
	outer = []
	inner = []
