import sys
import primer3

from Bio import SeqIO
from pprint import pprint

from primal.models import Region, Explain
import primal.settings


def find_primers(amplicon_length, overlap, window_size, references, seq, region_num, start, verbose=False):
	"""
	Given a list of biopython SeqRecords (references), and a string representation
	of the pimary reference (seq), return a list of Region objects containing candidate
	primer pairs sorted by an alignment score summed over all references.
	"""

	# Primer3 setup
	p3_global_args = settings.outer_params
	region_key = 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'
	p3_seq_args = {
		region_key: [
			start,
			start + window_size,
			start + amplicon_length - window_size,
			start + amplicon_length],
		'SEQUENCE_TEMPLATE': seq,
		'SEQUENCE_INCLUDED_REGION': [0, len(seq) - 1],
	}
	p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[amplicon_length-(2*window_size), amplicon_length]]

	if verbose:
		print "Primer3 Settings:"
		pprint(p3_seq_args, width=1)
		pprint(p3_global_args, width=1)

	# make some decision based on remaining seq length

	while True:
		p3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)

		left_ok = Explain(p3_output['PRIMER_LEFT_EXPLAIN']).ok
		right_ok = Explain(p3_output['PRIMER_RIGHT_EXPLAIN']).ok
		if left_ok > 0 and right_ok > 0:
			break

		# Step right and try again
		print "Stepping right +10"
		p3_seq_args[region_key] = [min(x+10, len(seq)) for x in p3_seq_args[region_key]]

		if p3_seq_args[region_key][0] >= overlap:
			print "Total failure!"
			sys.exit()

	position = 0
	return Region(amplicon_length, region_num, p3_output, references, position)


def multiplex(args, parser=None):
	references = list(SeqIO.parse(open(args.g, 'r'), 'fasta'))
	results =[]
	start = 0
	region_num = 0

	while True:
		region_num += 1
		print "Region {}".format(region_num)
		print "Start: {}".format(start)
		region = find_primers(args.amplicon_length, args.overlap, 50, references,
					   str(references[0].seq), region_num, start, verbose=args.verbose)
		print "Right end = {}".format(region.pairs[0].right.end)
		print "Right start = {}".format(region.pairs[0].right.start)
		start = region.pairs[0].right.end - args.overlap
		results.append(region)

	return results
