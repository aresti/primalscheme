import sys
import primer3

from Bio import SeqIO
from pprint import pprint

from primal.models import Region, Explain
import primal.settings


def find_primers(amplicon_length, overlap, window_size, references, seq, region_num, start, same_pool_limit, verbose=False, very_verbose=False):
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
			min(start, len(seq)-amplicon_length),
			window_size,
			min(start+amplicon_length-window_size, len(seq)-window_size),
			window_size],
		'SEQUENCE_TEMPLATE': seq,
		'SEQUENCE_INCLUDED_REGION': [0, len(seq)-1],
	}
	# Product size does not include primer sequences
	p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[amplicon_length-(2*window_size)+(2*22), amplicon_length+(2*22)]]

	if very_verbose:
		print "\nPrimer3 Settings:"
		pprint(p3_seq_args, width=1)
		pprint(p3_global_args, width=1)

	while True:
		p3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
		if very_verbose:
			pprint(p3_output, width=1)
	
		left_ok = Explain(p3_output['PRIMER_LEFT_EXPLAIN']).ok
		right_ok = Explain(p3_output['PRIMER_RIGHT_EXPLAIN']).ok
		if verbose:
			print 'ok', left_ok, right_ok
		if left_ok > 0 and right_ok > 0:
			break

		# Step right and try again
		print "Stepping left -10"
		p3_seq_args[region_key][0] -= 10
		p3_seq_args[region_key][2] -= 10

		# Check if overlap is too large
		if p3_seq_args[region_key][0] < same_pool_limit:
			print "Failure! Please try different settings"
			sys.exit()

	return Region(amplicon_length, region_num, p3_output, references)


def multiplex(args, parser=None):
	references = list(SeqIO.parse(open(args.g, 'r'), 'fasta'))
	results =[]
	start = 0
	region_num = 0
	window_size = 50
	overlap = 0
	same_pool_limit = 0

	while True:
		region_num += 1
		if len(results) > 2:
			same_pool_limit = results[-2].pairs[0].right.start
		region = find_primers(args.amplicon_length, args.overlap, window_size, 
			references, str(references[0].seq), region_num, start, same_pool_limit, 
			verbose=args.verbose, very_verbose=args.very_verbose)
		print "Region %i, %i:%i" %( region_num, region.pairs[0].left.start, 
			region.pairs[0].right.start)
		if results:
			overlap = results[-1].pairs[0].right.start - region.pairs[0].left.end
		print "Product length %i, overlap %i" %(region.pairs[0].right.start - 
			region.pairs[0].left.start, overlap)
		if args.very_verbose:
			pprint(vars(region.pairs[0].left))
			pprint(vars(region.pairs[0].right))
		start = region.pairs[0].right.end - args.overlap
		results.append(region)
		if region.pairs[0].right.start >= len(references[0].seq)-window_size:
			print 'Finished!'
			break
		

	return results
