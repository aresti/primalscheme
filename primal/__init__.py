import os
import sys
import pickle
import primer3
import primal.settings

from Bio import SeqIO
from .models import Region, Explain
from plot import plot_schemeadelica


class NoSuitableException(Exception):
    """No suitable primer found"""
    pass

class MaxGapException(Exception):
    """Maximum gap exceeded, increase amplicon length"""
    pass

def step(logger, p3_seq_args, step_type, search_space, start_limits, region_key, region_num, keep_right, step_size=settings.STEP_DISTANCE):
    print 'stepping %s' %(step_type)
    if step_type == 'left':
        if p3_seq_args[region_key][0] < step_size:
            p3_seq_args[region_key][0] = 0
            p3_seq_args[region_key][1] = search_space
        else:
            p3_seq_args[region_key][0] -= step_size
            p3_seq_args[region_key][1] += step_size
    if step_type == 'right':
        p3_seq_args[region_key][0] = 0
        p3_seq_args[region_key][1] += step_size

    logger.info("Region %i: step type %s, range %i:%i, limit %s, keep right=%s" %(region_num, step_type, p3_seq_args[region_key][0] + start_limits[0], p3_seq_args[region_key][0] + start_limits[0] + p3_seq_args[region_key][1], str(start_limits[0]) if step_type == 'left' else 'none', keep_right))
    return p3_seq_args

def find_primers(prefix, logger, amplicon_length, min_overlap, max_gap, search_space, max_candidates, references, region_num, start_limits, debug, is_last_region):
    """
    Given a list of biopython SeqRecords (references), and a string representation
    of the pimary reference (seq), return a list of Region objects containing candidate
    primer pairs sorted by an alignment score summed over all references.
    """
    logger.info("Region %i: forward primer limits %i:%i" %(region_num, start_limits[0], start_limits[1]))

    # Slice primary reference to speed up Primer3 on long sequences
    chunk_end = min(len(references[0]), int(1.1 * \
        (amplicon_length + max_gap) if region_num == 1 else start_limits[1] + 1.1 * \
        (amplicon_length + max_gap)))
    seq = str(references[0].seq)[start_limits[0]:chunk_end]
    logger.info("Region %i: reference slice %i:%i, length %i" %(region_num, start_limits[0], chunk_end, len(seq)))

    # Primer3 setup
    p3_global_args = settings.outer_params
    region_key = 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'

    #Reset to 0 to prevent invalid region key
    region_end = search_space
    if region_num == 1:
        region_start = 0
    elif is_last_region:
        region_start = len(seq) - amplicon_length
    else:
        if start_limits[1] - start_limits[0] - search_space < 0:
            region_start = 0
            region_end = search_space + (start_limits[1] - start_limits[0] - search_space)
        else:
            region_start = start_limits[1] - start_limits[0] - search_space

    p3_seq_args = {
        region_key: [region_start, region_end, -1, -1],
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_INCLUDED_REGION': [0, len(seq) - 1]
    }

    p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
        [int(amplicon_length * 0.9), int(amplicon_length * 1.1)]]
    p3_global_args['PRIMER_NUM_RETURN'] = max_candidates

    if debug:
        for a, b in [(b, p3_seq_args[b]) for b in p3_seq_args.keys()]:
            logger.debug('%s:%s' %(a, b))

    #print start_limits
    keep_right = False
    while True:
        print p3_seq_args[region_key]
        primer3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
        #if debug:
        #    logger.debug(primer3_output)

        num_returned = primer3_output['PRIMER_PAIR_NUM_RETURNED']
        logger.info("Region %i: current settings returned %i pairs" %(region_num, num_returned))
        if num_returned:
            break

        #if region_num == 1 or keep_right == True:
        if p3_seq_args[region_key][0] == 0 or keep_right:
            p3_seq_args = step(logger, p3_seq_args, 'right', search_space, start_limits, region_key, region_num, keep_right)
            keep_right = True

        else:
            p3_seq_args = step(logger, p3_seq_args, 'left', search_space, start_limits, region_key, region_num, keep_right)
            if p3_seq_args[region_key][0] < 0:
                keep_right = True

        print start_limits[0] + p3_seq_args[region_key][0] + p3_seq_args[region_key][1], len(references[0])
        if start_limits[0] + p3_seq_args[region_key][0] + p3_seq_args[region_key][1] > len(references[0]):
            return

    return Region(prefix, region_num, max_candidates, start_limits, primer3_output, references)


def write_bed(prefix, results, reference_id, path='./'):
    filepath = os.path.join(path, '{}.bed'.format(prefix))
    with open(filepath, 'w') as bedhandle:
        for r in results:
            print >>bedhandle, '\t'.join(map(
                str, [reference_id, r.candidate_pairs[0].left.start, r.candidate_pairs[0].left.end, r.candidate_pairs[0].left.name, r.pool]))
            print >>bedhandle, '\t'.join(map(str, [reference_id, r.candidate_pairs[0].right.end,
                                                   r.candidate_pairs[0].right.start, r.candidate_pairs[0].right.name, r.pool]))


def write_tsv(prefix, results, path='./'):
    filepath = os.path.join(path, '{}.tsv'.format(prefix))
    with open(filepath, 'w') as tsvhandle:
        print >>tsvhandle, '\t'.join(
            ['name', 'seq', 'length', '%gc', 'tm (use 65)'])
        for r in results:
            left = r.candidate_pairs[0].left
            right = r.candidate_pairs[0].right
            print >>tsvhandle, '\t'.join(
                map(str, [left.name, left.seq, left.length, left.gc, left.tm]))
            print >>tsvhandle, '\t'.join(
                map(str, [right.name, right.seq, right.length, right.gc, right.tm]))

def write_pickle(prefix, results, path='./'):
    filepath = os.path.join(path, '{}.pickle'.format(prefix))
    with open(filepath, 'wb') as pickleobj:
        pickle.dump(results, pickleobj)

def write_refs(prefix, refs, path='./'):
    filepath = os.path.join(path, '{}.fasta'.format(prefix))
    with open(filepath, 'w') as refhandle:
        SeqIO.write(refs, filepath, 'fasta')


def smart(args, logger, parser=None):
    return


def multiplex(args, logger, parser=None):

    references = args.references
    results = []
    region_num = 0
    window_size = 50

    # Check for sensible amplicon length
    if args.amplicon_length < 100 or args.amplicon_length > 2000:
        logger.error(
            "--amplicon-length set to %s, must be between 100 and 2000" % (args.amplicon_length))
        raise ValueError(
            "--amplicon-length set to %s, must be between 100 and 2000" % (args.amplicon_length))

    # Check there are enough regions to allow limits and end logic to work
    if len(references[0]) < 2 * args.amplicon_length:
        # This needs changing to length > reference for sinlge region genomes
                # See other pool overlap exception
        logger.error("Reference length %s, must be twice amplicon length" % (
            len(references[0])))
        raise ValueError(
            "Reference length %s, must be twice amplicon length" % (len(references[0])))

    while True:
        region_num += 1
        print "Region %i" %(region_num)
        if region_num > 2:
        # Left limit prevents crashing into the previous primer in this pool
            if results[-1].candidate_pairs[0].left.start > results[-2].candidate_pairs[0].right.start:
                # If gap we need to change left_start_limit
                left_start_limit = results[-1].candidate_pairs[0].left.end + 1
            else:
                left_start_limit = results[-2].candidate_pairs[0].right.start + 1
        else :
            left_start_limit = 0

        # Right limit maintains a minimum overlap of 0 (no gap)
        right_start_limit = results[-1].candidate_pairs[0].right.end - \
            args.min_overlap - 1 if region_num > 1 else 0
        #print (left_start_limit, right_start_limit)

        if region_num > 2:
            if right_start_limit <= left_start_limit:
                raise ValueError("Amplicon length too short for specified overlap")

        # Determine if is last region
        #if region_num > 1:
        #    print len(references[0]) - results[-1].candidate_pairs[0].right.start, args.amplicon_length
        is_last_region = (region_num > 1 and len(
            references[0]) - results[-1].candidate_pairs[0].right.start < args.amplicon_length)
        print is_last_region
        try:
            region = find_primers(args.prefix, logger, args.amplicon_length, args.min_overlap, args.max_gap, args.search_space,
                                  args.max_candidates, references, region_num, (left_start_limit, right_start_limit), args.debug, is_last_region)
            results.append(region)
            print results

            #Report scores and alignments
            if args.debug:
                #Don't report aligment to primary reference
                for i in range(1, len(args.references)):
                    logger.debug(
                        results[-1].candidate_pairs[0].left.alignments[i].formatted_alignment)
                    logger.debug('Subtotal: %i' % (
                        results[-1].candidate_pairs[0].left.alignments[i].score))
                    logger.debug(
                        results[-1].candidate_pairs[0].right.alignments[i].formatted_alignment)
                    logger.debug('Subtotal: %i' % (
                        results[-1].candidate_pairs[0].right.alignments[i].score))
                logger.debug('Totals for sorted candidates: ' + ','.join(['%.2f' %each.total for each in results[-1].candidate_pairs]))
                logger.debug('Right start for sorted candidates: ' + ','.join(['%i' %each.right.start for each in results[-1].candidate_pairs]))
                logger.debug('Right end for sorted candidates: ' + ','.join(['%i' %each.right.end for each in results[-1].candidate_pairs]))
                logger.debug('Length for sorted candidates: ' + ','.join(['%i' %each.right.length for each in results[-1].candidate_pairs]))

            #Reporting region
            if region_num > 1:
                # Results now include this one, so -2 is the other pool
                trimmed_overlap = results[-2].candidate_pairs[0].right.end - \
                    results[-1].candidate_pairs[0].left.end - 1
                logger.info("Region %i: highest scoring product %i:%i, length %i, trimmed overlap %i" % (region_num, results[-1].candidate_pairs[0].left.start, results[-1].candidate_pairs[0].right.start, results[-1].candidate_pairs[0].product_length, trimmed_overlap))
            else:
                logger.info("Region %i: highest scoring product %i:%i, length %i" %
                            (region_num, results[-1].candidate_pairs[0].left.start, results[-1].candidate_pairs[0].right.start, results[-1].candidate_pairs[0].product_length))
            if args.debug:
                for a in vars(results[-1].candidate_pairs[0].left):
                    logger.info('%s:%s' %(a, str(vars(results[-1].candidate_pairs[0].left)[a])))
                for b in vars(results[-1].candidate_pairs[0].right):
                    logger.info('%s:%s' %(b, str(vars(results[-1].candidate_pairs[0].right)[b])))
        except NoSuitableException:
            logger.warning("No suitable primer found for region")
            pass

        #Handle the end so maximum uncovered genome is one overlaps length
        if region_num > 1:
            #print len(references[0]) - results[-1].candidate_pairs[0].right.start, args.amplicon_length
            if len(references[0]) - results[-2].candidate_pairs[0].right.start < args.amplicon_length:
                break

    #Clean up
    print results
    if results[-1] == None:
            del results[-1]

    # Write bed and image files
    write_pickle(args.prefix, results, args.output_path)
    logger.info('pickled results')
    write_bed(args.prefix, results, references[0].id, args.output_path)
    logger.info('BED written')
    write_tsv(args.prefix, results, args.output_path)
    logger.info('TSV written')
    write_refs(args.prefix, references, args.output_path)
    logger.info('FASTA written')
    plot_schemeadelica(args.prefix, results, references[0], args.output_path)
    logger.info('PDF written')
    logger.info('Complete')

    return results
