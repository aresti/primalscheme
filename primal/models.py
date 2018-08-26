import logging
import os
import sys
import pickle
import primer3
import re
import settings

from Bio import pairwise2, Seq, SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors

from exceptions import NoSuitableException

sys.path.append('Porechop/porechop')
from cpp_function_wrappers import adapter_alignment


logger = logging.getLogger('Primal Log')


class Primer(object):
    """A simple primer."""

    def __init__(self, direction, name, seq):
        # TODO: Validate direction is LEFT or RIGHT
        self.direction = direction
        self.name = name
        self.seq = seq

    @property
    def length(self):
        return len(self.seq)


class CandidatePrimer(Primer):
    """A candidate primer for a region."""

    def __init__(self, direction, name, seq, start, gc, tm, references):
        super(CandidatePrimer, self).__init__(direction, name, seq)
        self.start = start
        self.gc = gc
        self.tm = tm

        self.sub_total = 0
        self.alignments = []

        for ref in references:
            alignment = CAlignment(self, ref)
            self.alignments.append(alignment)
            self.sub_total += alignment.score

    @property
    def end(self):
        if self.direction == 'LEFT':
            return self.start + self.length
        else:
            return self.start - self.length


class CandidatePrimerPair(object):
    """A pair of candidate primers for a region."""

    def __init__(self, left, right):
        self.left = left
        self.right = right
        self.total = left.sub_total + right.sub_total

    @property
    def product_length(self):
        return self.right.start - self.left.start + 1


class Region(object):
    """A region that forms part of a scheme."""

    def __init__(self, prefix, region_num, max_candidates, chunk_start, primer3_output, references):
        self.region_num = region_num
        self.pool = '2' if self.region_num % 2 == 0 else '1'
        self.candidate_pairs = []

        for cand_num in range(max_candidates):
            lenkey = 'PRIMER_LEFT_%s' % (cand_num)
            left_name = '%s_%i_%s' % (prefix, region_num, 'LEFT')
            right_name = '%s_%i_%s' % (prefix, region_num, 'RIGHT')
            if lenkey not in primer3_output:
                break

            left_seq = str(primer3_output['PRIMER_LEFT_%i_SEQUENCE' % (cand_num)])
            right_seq = str(primer3_output['PRIMER_RIGHT_%i_SEQUENCE' % (cand_num)])

            left_start = int(primer3_output['PRIMER_LEFT_%i' % (cand_num)][0] + chunk_start)
            right_start = int(primer3_output['PRIMER_RIGHT_%i' % (cand_num)][0] + chunk_start + 1)

            left_gc = float(primer3_output['PRIMER_LEFT_%i_GC_PERCENT' % (cand_num)])
            right_gc = float(primer3_output['PRIMER_RIGHT_%i_GC_PERCENT' % (cand_num)])

            left_tm = float(primer3_output['PRIMER_LEFT_%i_TM' % (cand_num)])
            right_tm = float(primer3_output['PRIMER_RIGHT_%i_TM' % (cand_num)])

            left = CandidatePrimer('LEFT', left_name, left_seq, left_start, left_gc, left_tm, references)
            right = CandidatePrimer('RIGHT', right_name, right_seq, right_start, right_gc, right_tm, references)

            self.candidate_pairs.append(CandidatePrimerPair(left, right))
        # Select the highest scoring pair with the rightmost position
        self.candidate_pairs.sort(key=lambda x: (x.total, x.right.end), reverse=True)

    @property
    def top_pair(self):
        return self.candidate_pairs[0]

class CAlignment(object):
    """An alignment of a primer against a reference."""
    def __init__(self, primer, ref):
        if primer.direction == 'LEFT':
            alignment_result = adapter_alignment(str(ref.seq), str(primer.seq), [2, -1, -2, -1])
        elif primer.direction == 'RIGHT':
            alignment_result = adapter_alignment(str(ref.seq.reverse_complement()), str(primer.seq), [2, -1, -2, -1])
        result_parts = alignment_result.split(',')
        ref_start = int(result_parts[0])

        # If the read start is -1, that indicates that the alignment failed completely.
        if ref_start == -1:
            self.score = 0.0
            self.formatted_alignment = 'None'
        else:
            ref_end = int(result_parts[1]) + 1
            #aligned_region_percent_identity = float(result_parts[5])
            full_primer_percent_identity = float(result_parts[6])

            if primer.direction == 'LEFT':
                self.start = ref_start
                self.end = ref_end
                self.length = self.end - self.start
            else:
                self.start = len(ref) - ref_start
                self.end = len(ref) - (int(result_parts[1]) + 1)
                self.length = self.start - self.end

            # Percentage identity for glocal alignment
            self.score = full_primer_percent_identity

            # Get alignment strings
            self.aln_query = result_parts[8][ref_start:ref_end]
            self.aln_ref = result_parts[7][ref_start:ref_end]
            self.aln_ref_comp = Seq.Seq(str(self.aln_ref)).complement()
            self.ref_id = ref.id
            self.mm_3prime = False

            # Make cigar
            self.cigar = ''
            for a, b in zip(self.aln_query, self.aln_ref):
                if a == '-' or b == '-':
                    self.cigar += ' '
                    continue
                if a != b:
                    self.cigar += '*'
                    continue
                else:
                    self.cigar += '|'

            # Format alignment
            short_primer = primer.name[:30] if len(primer.name) > 30 else primer.name
            short_ref = ref.id[:30] if len(ref.id) > 30 else ref.id
            self.formatted_alignment = "\n{: <30}5\'-{}-3\'\n{: <33}{}\n{: <30}3\'-{}-5\'".format(short_primer, self.aln_query, '', self.cigar, short_ref, self.aln_ref_comp)

            # Check 3' mismatches
            if set([self.aln_query[-1], self.aln_ref_comp[-1]]) in settings.MISMATCHES:
                self.mm_3prime = True
                self.score = 0


class Alignment(object):
    """An alignment of a primer against a reference."""

    def __init__(self, primer, ref):
        # Do alignments
        if primer.direction == 'LEFT':
            search_start = primer.start - 100 if primer.start > 100 else 0
            search_end = primer.end + 100 if primer.end + 100 <= len(ref) else len(ref)
            alns = pairwise2.align.globalms(str(primer.seq), str(ref.seq[search_start:search_end]), 2, -1, -2, -1, penalize_end_gaps=False, one_alignment_only=True)
        elif primer.direction == 'RIGHT':
            search_start = primer.end - 100 if primer.start > 100 else 0
            search_end = primer.start + 100 if primer.start + 100 <= len(ref) else len(ref)
            alns = pairwise2.align.globalms(str(primer.seq), str(ref.seq[search_start:search_end].reverse_complement()), 2, -1, -2, -1, penalize_end_gaps=False, one_alignment_only=True)
        if alns:
            aln = alns[0]

            p = re.compile('(-*)([ACGTN][ACGTN\-]*[ACGTN])(-*)')
            m = list(re.finditer(p, str(aln[0])))[0]

            if primer.direction == 'LEFT':
                self.start = search_start + m.span(2)[0]
                self.end = search_start + m.span(2)[1]
                self.length = self.end - self.start
            else:
                self.start = search_end - m.span(2)[0]
                self.end = search_end - m.span(2)[1]
                self.length = self.start - self.end

            # Normalise alignment score by length
            self.score = aln[2] / self.length

            # Get alignment strings
            self.aln_query = aln[0][m.span(2)[0]:m.span(2)[1]]
            self.aln_ref = aln[1][m.span(2)[0]:m.span(2)[1]]
            self.aln_ref_comp = Seq.Seq(str(self.aln_ref)).complement()
            self.ref_id = ref.id
            self.mm_3prime = False

            # Make cigar
            self.cigar = ''
            for a, b in zip(self.aln_query, self.aln_ref):
                if a == '-' or b == '-':
                    self.cigar += ' '
                    continue
                if a != b:
                    self.cigar += '*'
                    continue
                else:
                    self.cigar += '|'

            # Format alignment
            short_primer = primer.name[:30] if len(primer.name) > 30 else primer.name
            short_ref = ref.id[:30] if len(ref.id) > 30 else ref.id
            self.formatted_alignment = "\n{: <30}5\'-{}-3\'\n{: <33}{}\n{: <30}3\'-{}-5\'".format(short_primer, self.aln_query, '', self.cigar, short_ref, self.aln_ref_comp)

            # Check 3' mismatches
            if set([self.aln_query[-1], self.aln_ref_comp[-1]]) in settings.MISMATCHES:
                self.mm_3prime = True
                self.score = 0

        else:
            self.score = 0
            self.formatted_alignment = 'None found'


class MultiplexScheme(object):
    """A complete multiplex primer scheme."""

    def __init__(self, references, amplicon_length, min_overlap=20, max_gap=100, window_size=50, search_space=30,
                 max_candidates=10, step_size=1, prefix='PRIMAL_SCHEME'):
        self.references = references
        self.amplicon_length = amplicon_length
        self.min_overlap = min_overlap
        self.max_gap = max_gap
        self.window_size = window_size
        self.search_space = search_space
        self.max_candidates = max_candidates
        self.step_size = step_size
        self.prefix = prefix
        self.regions = []

        self.run()

    @property
    def primary_reference(self):
        return self.references[0]

    def run(self):
        regions = []
        region_num = 0

        while True:
            region_num += 1
            prev_pair = regions[-1].candidate_pairs[0] if region_num > 1 else None
            prev_pair_same_pool = regions[-2].candidate_pairs[0] if region_num > 2 else None

            # Left start limit
            if prev_pair_same_pool:
                # Left limit prevents crashing into the previous primer in this pool
                if prev_pair.left.start > prev_pair_same_pool.right.start:
                    # Where there is a gap, we need to change left_primer_left_limit
                    left_primer_left_limit = prev_pair.left.end + 1
                else:
                    left_primer_left_limit = prev_pair_same_pool.right.start + 1
            else:
                left_primer_left_limit = 0

            # Right start limit; maintains a minimum overlap of 0 (no gap)
            left_primer_right_limit = prev_pair.right.end - self.min_overlap - 1 if prev_pair else 0

            #if prev_pair_same_pool and left_primer_right_limit <= left_primer_left_limit:
            #    raise ValueError("Amplicon length too short for specified overlap")

            is_last_region = (region_num > 1 and len(self.primary_reference) - prev_pair.right.start < self.amplicon_length)

            #Log limits
            logger.debug('Region {}: forward primer limits {}:{}'.format(region_num, left_primer_left_limit, left_primer_right_limit))

            # Slice primary reference to constrain primer positions to where we want them
            max_variation = 0.1
            if region_num == 1:
                chunk_start = 0
                chunk_end = (1 + max_variation/2) * self.amplicon_length
            elif is_last_region:
                # Last time work backwards
                chunk_end = len(self.primary_reference)
                chunk_start = len(self.primary_reference) - ((1 + max_variation/2) * self.amplicon_length)
            else:
                # right limit - min overlap - diff max min product length - max primer length
                chunk_start = left_primer_right_limit - (0.1 * self.amplicon_length) - settings.global_args['PRIMER_MAX_SIZE']
                chunk_end = chunk_start + ((1 + max_variation/2) * self.amplicon_length)

            # Find primers
            try:
                region = self._find_primers_only(region_num, left_primer_left_limit, left_primer_right_limit, chunk_start, chunk_end, self.min_overlap)
                regions.append(region)
            except NoSuitableException:
                pass

            # Handle the end; maximum uncovered genome is one overlap's length
            if prev_pair and len(self.primary_reference) - prev_pair.right.start < self.amplicon_length:
                break

            #Report scores and alignments
            #Don't report aligment to primary reference
            for i in range(1, len(self.references)):
                logger.debug(regions[-1].candidate_pairs[0].left.alignments[i].formatted_alignment)
                logger.debug('Subtotal: %i' % (regions[-1].candidate_pairs[0].left.alignments[i].score))
                logger.debug(regions[-1].candidate_pairs[0].right.alignments[i].formatted_alignment)
                logger.debug('Subtotal: %i' % (regions[-1].candidate_pairs[0].right.alignments[i].score))
                logger.debug('Totals for sorted candidates: ' + ','.join(['%.2f' %each.total for each in regions[-1].candidate_pairs]))
                logger.debug('Right start for sorted candidates: ' + ','.join(['%i' %each.right.start for each in regions[-1].candidate_pairs]))
                logger.debug('Right end for sorted candidates: ' + ','.join(['%i' %each.right.end for each in regions[-1].candidate_pairs]))
                logger.debug('Length for sorted candidates: ' + ','.join(['%i' %each.right.length for each in regions[-1].candidate_pairs]))

            if region_num > 1:
            # Remember, results now include this one, so -2 is the other pool
                trimmed_overlap = regions[-2].candidate_pairs[0].right.end - regions[-1].candidate_pairs[0].left.end - 1
                logger.info("Region %i: highest scoring product %i:%i, length %i, trimmed overlap %i" % (region_num, regions[-1].candidate_pairs[0].left.start, regions[-1].candidate_pairs[0].right.start, regions[-1].candidate_pairs[0].product_length, trimmed_overlap))
            else:
                logger.info("Region %i: highest scoring product %i:%i, length %i" % (region_num, regions[-1].candidate_pairs[0].left.start, regions[-1].candidate_pairs[0].right.start, regions[-1].candidate_pairs[0].product_length))

        self.regions = regions
        print regions
        if regions[-1] == None:
            del regions[-1]

    def write_bed(self, path='./'):
        logger.info('Writing BED')
        filepath = os.path.join(path, '{}.bed'.format(self.prefix))
        with open(filepath, 'w') as bedhandle:
            for r in self.regions:
                print >>bedhandle, '\t'.join(map(
                    str, [self.primary_reference.id, r.top_pair.left.start, r.top_pair.left.end, r.top_pair.left.name, r.pool]))
                print >>bedhandle, '\t'.join(map(str, [self.primary_reference.id, r.top_pair.right.end,
                                                       r.top_pair.right.start, r.top_pair.right.name, r.pool]))

    def write_tsv(self, path='./'):
        logger.info('Writing TSV')
        filepath = os.path.join(path, '{}.tsv'.format(self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print >>tsvhandle, '\t'.join(
                ['name', 'seq', 'length', '%gc', 'tm (use 65)'])
            for r in self.regions:
                left = r.top_pair.left
                right = r.top_pair.right
                print >>tsvhandle, '\t'.join(
                    map(str, [left.name, left.seq, left.length, left.gc, left.tm]))
                print >>tsvhandle, '\t'.join(
                    map(str, [right.name, right.seq, right.length, right.gc, right.tm]))

    def write_pickle(self, path='./'):
        logger.info('Writing pickles')
        filepath = os.path.join(path, '{}.pickle'.format(self.prefix))
        with open(filepath, 'wb') as pickleobj:
            pickle.dump(self.regions, pickleobj)

    def write_refs(self, path='./'):
        logger.info('Writing references')
        filepath = os.path.join(path, '{}.fasta'.format(self.prefix))
        with open(filepath, 'w') as refhandle:
            SeqIO.write(self.references, filepath, 'fasta')

    def write_schemadelica_plot(self, path='./'):
        logger.info('Writing plot')
        gd_diagram = GenomeDiagram.Diagram("Primer Scheme", track_size=1)
        scale_track = GenomeDiagram.Track(
            name='scale', scale=True, scale_fontsize=10, scale_largetick_interval=1000, height=0.1)
        gd_diagram.add_track(scale_track, 2)

        primer_feature_set_1 = GenomeDiagram.FeatureSet()
        primer_feature_set_2 = GenomeDiagram.FeatureSet()

        for r in self.regions:
            cols1 = [self.primary_reference.id, r.top_pair.left.start,
                     r.top_pair.left.end, r.top_pair.left.name, r.pool]
            cols2 = [self.primary_reference.id, r.top_pair.right.end,
                     r.top_pair.right.start, r.top_pair.right.name, r.pool]
            region = str(r.region_num)
            fwd_feature = SeqFeature(FeatureLocation(
                int(cols1[1]), int(cols1[2]), strand=0))
            rev_feature = SeqFeature(FeatureLocation(
                int(cols2[1]), int(cols2[2]), strand=0))
            region_feature = SeqFeature(FeatureLocation(
                int(cols1[1]), int(cols2[2]), strand=0))
            if int(region) % 2 == 0:
                primer_feature_set_1.add_feature(region_feature, color=colors.palevioletred,
                                                 name=region, label=True, label_size=10, label_position="middle", label_angle=0)
                primer_feature_set_1.add_feature(
                    fwd_feature, color=colors.red, name=region, label=False)
                primer_feature_set_1.add_feature(
                    rev_feature, color=colors.red, name=region, label=False)
            else:
                primer_feature_set_2.add_feature(region_feature, color=colors.palevioletred,
                                                 name=region, label=True, label_size=10, label_position="middle", label_angle=0)
                primer_feature_set_2.add_feature(
                    fwd_feature, color=colors.red, name=region, label=False)
                primer_feature_set_2.add_feature(
                    rev_feature, color=colors.red, name=region, label=False)

        primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
        primer_track.add_set(primer_feature_set_1)
        gd_diagram.add_track(primer_track, 4)

        primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
        primer_track.add_set(primer_feature_set_2)
        gd_diagram.add_track(primer_track, 6)

        rows = max(2, int(round(len(self.primary_reference) / 10000.0)))
        gd_diagram.draw(format='linear', pagesize=(300 * rows, 200 * rows),
                        fragments=rows, start=0, end=len(self.primary_reference))

        pdf_filepath = os.path.join(path, '{}.pdf'.format(self.prefix))
        svg_filepath = os.path.join(path, '{}.svg'.format(self.prefix))
        gd_diagram.write(pdf_filepath, 'PDF', dpi=300)
        gd_diagram.write(svg_filepath, 'SVG', dpi=300)

    def _find_primers_only(self, region_num, left_limit, right_limit, chunk_start, chunk_end, overlap):
        """
        Find primers for a given region.

        Given a list of biopython SeqRecords (references), and a string representation
        of the pimary reference (seq), return a list of Region objects containing candidate
        primer pairs sorted by an alignment score summed over all references.
        """

        # Primer3 setup
        p3_global_args = settings.global_args
        p3_seq_args = settings.seq_args
        max_variation = 0.1
        p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
            [int(self.amplicon_length * (1-max_variation/2)), int(self.amplicon_length * (1+max_variation/2))]]
        p3_global_args['PRIMER_NUM_RETURN'] = self.max_candidates


        # Run Primer3 until primers are found
        hit_left_limit = False
        while True:
            chunk_end = int(chunk_end)
            chunk_start = int(chunk_start)
            seq = str(self.primary_reference.seq[chunk_start:chunk_end])
            p3_seq_args['SEQUENCE_TEMPLATE'] = seq
            p3_seq_args['SEQUENCE_INCLUDED_REGION'] = [0, len(seq) - 1]
            logger.debug("Region %i: reference chunk %i:%i, length %i" %(region_num, chunk_start, chunk_end, len(seq)))
            primer3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
            num_returned = primer3_output['PRIMER_PAIR_NUM_RETURNED']
            logger.debug("Region %i: current settings returned %i pairs" %(region_num, num_returned))
            if num_returned:
                break

            #Move right if first region or if hit left limit
            if region_num == 1 or hit_left_limit:
                logger.debug("Region %i: stepping right, position %i" %(region_num, chunk_start))
                chunk_start += self.step_size
                chunk_end += self.step_size
            #Move left for all other regions
            else:
                logger.debug("Region %i: stepping left, position %i, limit %s" %(region_num, chunk_start, left_limit))
                chunk_start -= self.step_size
                chunk_end -= self.step_size
                if chunk_start <= left_limit:
                    #This is correct for a step size of 1
                    logger.debug("Region %i: hit left limit" %(region_num))
                    chunk_start = right_limit - (0.1 * self.amplicon_length) - settings.global_args['PRIMER_MAX_SIZE']
                    chunk_end = chunk_start + ((1 + max_variation/2) * self.amplicon_length)
                    chunk_end = int(chunk_end)
                    hit_left_limit = True


        # right_limit not used by region so could be removed
        return Region(self.prefix, region_num, self.max_candidates, chunk_start, primer3_output,
                      self.references)


    def _find_primers(self, region_num, left_limit, right_limit, is_last_region):
        """
        Find primers for a given region.

        Given a list of biopython SeqRecords (references), and a string representation
        of the pimary reference (seq), return a list of Region objects containing candidate
        primer pairs sorted by an alignment score summed over all references.
        """
        logger.info('Processing region {}'.format(region_num))
        logger.debug('Region {}: forward primer limits {}:{}'.format(region_num, left_limit, right_limit))

        # Slice primary reference to speed up Primer3 on long sequences
        if region_num == 1:
            chunk_end = min(len(self.primary_reference), 1.1 * (self.amplicon_length + self.max_gap))
        else:
            chunk_end = min(len(self.primary_reference), right_limit + 1.1 * (self.amplicon_length + self.max_gap))
        chunk_end = int(chunk_end)
        seq = str(self.primary_reference.seq)[left_limit:chunk_end]

        # Primer3 setup
        p3_global_args = settings.outer_params
        region_key = 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'

        # Reset to 0 to prevent invalid region key
        region_end = self.search_space
        if region_num == 1:
            region_start = 0
        elif is_last_region:
            region_start = len(seq) - self.amplicon_length
        else:
            if right_limit - left_limit - self.search_space < 0:
                region_start = 0
                region_end = self.search_space + right_limit - left_limit - self.search_space
            else:
                region_start = right_limit - left_limit - self.search_space

        p3_seq_args = {
            region_key: [region_start, region_end, -1, -1],
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_INCLUDED_REGION': [0, len(seq) - 1]
        }
        p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
            [int(self.amplicon_length * 0.9), int(self.amplicon_length * 1.1)]]
        p3_global_args['PRIMER_NUM_RETURN'] = self.max_candidates
        keep_right = False

        while True:
            primer3_output = primer3.bindings.designPrimers(p3_seq_args, p3_global_args)
            num_returned = primer3_output['PRIMER_PAIR_NUM_RETURNED']
            if num_returned:
                break

            if p3_seq_args[region_key][0] == 0 or keep_right:
                step_type = 'right'
                p3_seq_args[region_key][0] = 0
                p3_seq_args[region_key][1] += self.step_size
                keep_right = True
            else:
                step_type = 'left'
                p3_seq_args[region_key][0] -= self.step_size
                p3_seq_args[region_key][1] += self.step_size
                if p3_seq_args[region_key][0] < 0:
                    keep_right = True

            logger.debug('Region {}: step type {}, range {}:{}, limit {}, keep right={}'.format(region_num, step_type, p3_seq_args[region_key][0] + left_limit, p3_seq_args[region_key][0] + left_limit + p3_seq_args[region_key][1], str(left_limit) if step_type == 'left' else 'none', keep_right))

            if left_limit + p3_seq_args[region_key][0] + p3_seq_args[region_key][1] > len(self.primary_reference):
                raise NoSuitableException

        return Region(self.prefix, region_num, self.max_candidates, (left_limit, right_limit), primer3_output,
                      self.references)
