import sys
from Bio import SeqIO
import primer3
from pprint import pprint
import settings
import interactions
import pairwise
from primerpair import primerpair
from openpyxl import load_workbook
from openpyxl.worksheet.datavalidation import DataValidation

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
	best_pair = []
	for region_pairs in outer_pairs:
		for pair in range(3):
			total = 0
			for record in records:
				left_name, left_score, left_3prime_mm, right_name, right_score, right_3prime_mm = pairwise.fast_pairwise(record, region_pairs, pair)
				pair_scores.append((left_name, left_score, left_3prime_mm, right_name, right_score, right_3prime_mm))
				total += (left_score + right_score)
				if (left_3prime_mm or right_3prime_mm):
					total = 0
			best_pair.append((region_pairs, pair, total))
	for each in best_pair:
		print each
	print
	for each in pair_scores:
		print each
				

	#check for interactions
	#print interactions.interactions(args, outer)

	with open(args.o + '_sigma.csv', 'w') as fileout:
		for i, each in enumerate(outer_pairs):
			pool = str(['1' if i%2==0 else '2'][0])
			print each.left_0_name, '\t', each.left_0_seq, '\t', pool
			print each.right_0_name, '\t', each.right_0_seq, '\t', pool

def nested_multiplex(parser, args):
        seq, seq_params = run(args)
	outer = []
	inner = []
