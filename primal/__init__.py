import os
import sys
import primer3
import primal.settings

from Bio import SeqIO
from pprint import pprint
from .models import Region, Explain
from plot import plot_schemeadelica


class PoolOverlapException(Exception):
	"""No suitable primer found between --min-overlap and right edge of last primer
	primer in the same pool"""
	pass


def find_primers(prefix, amplicon_length, min_overlap, search_space, references, region_num, start_limits, v=False, vvv=False):
	"""
	Given a list of biopython SeqRecords (references), and a string representation
	of the pimary reference (seq), return a list of Region objects containing candidate
	primer pairs sorted by an alignment score summed over all references.
	"""

	# Slice primary reference to speed up Primer3 on long sequences
	chunk_end = len(references[0]) if region_num == 1 else start_limits[1] + amplicon_length
	seq = str(references[0].seq)[start_limits[0]:chunk_end]


	# Primer3 setup
	p3_global_args = settings.outer_params
	region_key = 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'
	start = start_limits[1] if region_num == 1 else start_limits[1] - start_limits[0] - search_space
	p3_seq_args = {
		region_key: [min(start, len(seq) - amplicon_length), search_space, -1, -1],
		'SEQUENCE_TEMPLATE': seq,
		'SEQUENCE_INCLUDED_REGION': [0, len(seq) - 1],
	}


	allowed_variation = amplicon_length * 0.1
	p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[amplicon_length - allowed_variation, amplicon_length + allowed_variation]]


	while True:
		if vvv:
			print "\nPrimer3 Settings:"
			pprint(p3_seq_args, width=1)
			pprint(p3_global_args, width=1)

		primer3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
		#if vvv:
		#	pprint(primer3_output, width=1)

		num_pairs_returned = primer3_output['PRIMER_PAIR_NUM_RETURNED']
		if num_pairs_returned:
			break


		# Step right only if no primers in region 1
		if region_num == 1:
			p3_seq_args[region_key][1] += settings.STEP_DISTANCE
			print "Stepping right, position %i, limit %s" %(p3_seq_args[region_key][1], 'none')
		# Step left (increase overlap)
		else:
			p3_seq_args[region_key][0] -= settings.STEP_DISTANCE
			p3_seq_args[region_key][1] += settings.STEP_DISTANCE
			print "Stepping left, position %i, limit %i" %(p3_seq_args[region_key][0] + start_limits[0], start_limits[0])


		# Check if we've run into the left limit
		if p3_seq_args[region_key][0] < 0:
			raise PoolOverlapException("No suitable primers found for region {} with current parameters. Try adjusting --min-overlap and/or --amplicon-length.".format(region_num))

	return Region(prefix, region_num, start_limits, primer3_output, references)


def write_bed(prefix, results, reference_id, path='./'):
	filepath = os.path.join(path, '{}.bed'.format(prefix))
	with open(filepath, 'w') as bedhandle:
		for r in results:
			print >>bedhandle, '\t'.join(map(str, [reference_id, r.candidate_pairs[0].left.start, r.candidate_pairs[0].left.end, r.candidate_pairs[0].left.name, r.pool]))
			print >>bedhandle, '\t'.join(map(str, [reference_id, r.candidate_pairs[0].right.end, r.candidate_pairs[0].right.start, r.candidate_pairs[0].right.name, r.pool]))


def write_tsv(prefix, results, path='./'):
	filepath = os.path.join(path, '{}.tsv'.format(prefix))
	with open(filepath, 'w') as tsvhandle:
		print >>tsvhandle, '\t'.join(['name', 'seq', 'length', '%gc', 'tm (use 65)'])
		for r in results:
			left = r.candidate_pairs[0].left
			right = r.candidate_pairs[0].right
			print >>tsvhandle, '\t'.join(map(str, [left.name, left.seq, left.length, left.gc, left.tm]))
			print >>tsvhandle, '\t'.join(map(str, [right.name, right.seq, right.length, right.gc, right.tm]))


def multiplex(args, parser=None):
	# Check for sensible parameters
	if args.amplicon_length < 100 or args.amplicon_length > 2000:
		raise ValueError("--amplicon-length must be within the range 100 to 2000")

	references = args.references
	results =[]
	region_num = 0
	window_size = 50

	# Check there are enough regions to allow limits and end logic to work
	if len(references[0]) < (3 * args.amplicon_length - 2 * args.min_overlap):
		raise ValueError("length of reference must be at least ")

	while True:
		region_num += 1

		# Left limit prevents crashing into the previous primer in this pool
		left_start_limit = results[-2].candidate_pairs[0].right.start + 1 if region_num > 2 else 0

		# Right limit maintains a minimum overlap of 0 (no gap)
		right_start_limit = results[-1].candidate_pairs[0].right.end - args.min_overlap - 1 if region_num > 1 else 0

		# Find primers for this region
		region = find_primers(args.p, args.amplicon_length, args.min_overlap, args.search_space, references, region_num, (left_start_limit, right_start_limit), v=args.v, vvv=args.vvv)
		results.append(region)

		# Reporting
		print "\nRegion %i, %i:%i" %( region_num, results[-1].candidate_pairs[0].left.start,
			results[-1].candidate_pairs[0].right.start)
		if region_num > 1:
			# Remember, results now include this one, so -2 is the other pool
			trimmed_overlap = results[-2].candidate_pairs[0].right.end - results[-1].candidate_pairs[0].left.end - 1
 			print "Product length %i, trimmed overlap %i" % (results[-1].candidate_pairs[0].product_length, trimmed_overlap)
		else:
			print "Product length %i" % (results[-1].candidate_pairs[0].product_length)
		if args.vvv:
			pprint(vars(results[-1].candidate_pairs[0].left))
			pprint(vars(results[-1].candidate_pairs[0].right))

		# Handle the end so maximum uncovered genome is one overlaps length
		if region_num > 2:
			if len(references[0]) - results[-2].candidate_pairs[0].right.start < args.amplicon_length:
				print 'Finished!'
				break

	# write bed and image files
	write_bed(args.p, results, references[0].id, args.output_path)
	write_tsv(args.p, results, args.output_path)
	plot_schemeadelica(args.p, results, references[0], args.output_path)

	return results
