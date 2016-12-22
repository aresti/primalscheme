import sys
import primer3

from Bio import SeqIO
from pprint import pprint

import primal.settings
from primal.models import Region, Explain


def find_primers(length, overlap, window_size, references, seq, p3_global_args, results=[], verbose=False):
	p3_seq_args = {
		'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0, window_size, length-window_size, length],
		'SEQUENCE_TEMPLATE': seq,
		'SEQUENCE_INCLUDED_REGION': [0, len(seq)],
	}
	p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[length-(2*window_size), length]]
	if verbose:
		print pprint(p3_seq_args, width=1)
		print pprint(p3_global_args, width=1)
	#make some decision based on remaining seq length
	while True:
		output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
		#if verbose:
		#	pprint(output, width=1)
		left_ok = Explain(output['PRIMER_LEFT_EXPLAIN']).ok
		right_ok = Explain(output['PRIMER_RIGHT_EXPLAIN']).ok
		#print 'OK', left_ok, right_ok
		if left_ok > 0 and right_ok > 0:
			break
		#adjustments to window
		if left_ok == 0 or right_ok == 0:
			print 'Failed to find any primers'
			sys.exit()
	position = len(references[0])-len(seq)
	region = Region(length, len(results)+1, output, references, position)
	results.append(region)
	seq = seq[region.pairs[0].right.end-overlap:]
	if verbose:
		print 'Results len', len(results)
		print 'Position', position
		print 'Remaining', len(seq)
	if len(seq) < length:
		return results
	find_primers(length, overlap, window_size, references, seq, p3_global_args, results, verbose=verbose)


def multiplex(args, parser=None):
	references = list(SeqIO.parse(open(args.g, 'r'), 'fasta'))
	return find_primers(args.length, args.overlap, 50, references, str(references[0].seq), settings.outer_params, verbose=args.verbose)
