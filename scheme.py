import sys
from Bio import SeqIO
import primer3
from pprint import pprint
import settings
import interactions
import pairwise
from primal_models import Region
from operator import itemgetter

def run(args):
	#Use first record for primer picking
	references = list(SeqIO.parse(open(args.g, 'r'), 'fasta'))
	seq = str(references[0].seq)
	step = args.length - 2*(args.overlap)
	seq_params = {
		     'SEQUENCE_ID': 'blah',
		     'SEQUENCE_TEMPLATE': seq,
		     'SEQUENCE_INCLUDED_REGION': [0, len(seq)],
		     }
	return references, seq, step, seq_params 

def multiplex(args, parser=None):
	references, seq, step, seq_params = run(args)
	outer_params = settings.outer_pair(args)
	result = []
	print 'Running primer3'
	for i, region in enumerate(range(args.length/2, len(seq), step)):
		#sys.stdout.write('.')
		#sys.stdout.flush()
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

		region = Region(args.length, i+1, output, references)
		result.append(region)
		#pprint(output, width=1)
		#pprint(vars(region), width=1)
		#pprint(vars(region.pairs[0]), width=1)

	return result

	
def nested_multiplex(parser, args):
        seq, seq_params = run(args)
	outer = []
	inner = []
