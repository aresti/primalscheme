import os
import sys
import primer3
import primal.settings

from Bio import SeqIO
from pprint import pprint
from .models import Region, Explain
from plot import plot_schemeadelica


class PoolOverlapException(Exception):
	"""No suitable primer found between --overlap and right edge of last primer
	primer in the same pool"""
	pass


def find_primers(prefix, amplicon_length, overlap, window_size, references, seq, region_num, start, same_pool_limit, v=False, vvv=False):
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
	p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[amplicon_length-(2*window_size)+(2*p3_global_args['PRIMER_OPT_SIZE']), amplicon_length+(2*p3_global_args['PRIMER_OPT_SIZE'])]]

	if vvv:
		print "\nPrimer3 Settings:"
		pprint(p3_seq_args, width=1)
		pprint(p3_global_args, width=1)

	while True:
		primer3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
		if vvv:
			pprint(primer3_output, width=1)

		left_ok = Explain(primer3_output['PRIMER_LEFT_EXPLAIN']).ok
		right_ok = Explain(primer3_output['PRIMER_RIGHT_EXPLAIN']).ok
		if v:
			print 'ok', left_ok, right_ok
		if left_ok > 0 and right_ok > 0:
			break

		# Step right and try again
		print "Stepping left -10"
		p3_seq_args[region_key][0] -= 10
		p3_seq_args[region_key][2] -= 10

		# Check if overlap is too large
		if p3_seq_args[region_key][0] < same_pool_limit:
			raise PoolOverlapException("No suitable primers found for region {} with current parameters. Try adjusting --overlap and/or --amplicon-length.".format(region_num))

	return Region(prefix, region_num, primer3_output, references)


def write_bed(prefix, results, reference_id, path='./'):
	filepath = os.path.join(path, '{}.bed'.format(prefix))
	with open(filepath, 'w') as bedhandle:
		for r in results:
			print >>bedhandle, '\t'.join(map(str, [reference_id, r.candidate_pairs[0].left.start, r.candidate_pairs[0].left.end, r.candidate_pairs[0].left.name, r.pool]))
			print >>bedhandle, '\t'.join(map(str, [reference_id, r.candidate_pairs[0].right.end, r.candidate_pairs[0].right.start, r.candidate_pairs[0].right.name, r.pool]))


def multiplex(args, parser=None):
	# Check for sensible parameters
	if args.amplicon_length < 100 or args.amplicon_length > 2000:
		raise ValueError("--amplicon-length must be within the range 100 to 2000")

	references = list(SeqIO.parse(open(args.f, 'r'), 'fasta'))
	results =[]
	start = 0
	region_num = 0
	window_size = 50
	same_pool_limit = 0

	while True:
		region_num += 1

		# Same pool limit based on the previous region in this pool (-2)
		if region_num > 2:
			same_pool_limit = results[-2].candidate_pairs[0].right.start

		# Find primers for this region
		region = find_primers(args.p, args.amplicon_length, args.overlap, window_size,
				references, str(references[0].seq), region_num, start, same_pool_limit,
				v=args.v, vvv=args.vvv)
		results.append(region)

		# Reporting
		print "\nRegion %i, %i:%i" %( region_num, region.candidate_pairs[0].left.start,
			region.candidate_pairs[0].right.start)
		product_length = region.candidate_pairs[0].right.start - region.candidate_pairs[0].left.start
		if region_num > 1:
			# Remember, results now include this one, so -2 is the other pool
			overlap = results[-2].candidate_pairs[0].right.start - region.candidate_pairs[0].left.end
			print "Product length %i, overlap %i" % (product_length, overlap)
		else:
			print "Product length %i" % (product_length)
		if args.v:
			pprint(vars(region.candidate_pairs[0].left))
			pprint(vars(region.candidate_pairs[0].right))
		start = region.candidate_pairs[0].right.end - args.overlap
		if region.candidate_pairs[0].right.start >= len(references[0].seq)-window_size:
			print 'Finished!'
			break

	# write bed and image files
	write_bed(args.p, results, references[0].id, args.output_path)
	plot_schemeadelica(args.p, results, references[0], args.output_path)

	return results
