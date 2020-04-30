"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the MutilplexScheme object. Instantiation of this object
will result in a complete multiplex primer scheme.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""

import logging
import primer3

from primalscheme import settings
from primalscheme.exceptions import NoSuitablePrimers
from primalscheme.models import CandidatePrimer, CandidatePrimerPair, Region

logger = logging.getLogger('Primal Log')


class MultiplexScheme(object):
    """A complete multiplex primer scheme."""

    def __init__(self, references, amplicon_length, min_overlap, max_gap,
                 max_alts, max_candidates, step_size, max_variation,
                 prefix='primalscheme'):
        self.references = references
        self.amplicon_length = amplicon_length
        self.min_overlap = min_overlap
        self.max_gap = max_gap
        self.max_alts = max_alts
        self.max_candidates = max_candidates
        self.step_size = step_size
        self.max_variation = max_variation
        self.prefix = prefix
        self.regions = []

        self.run()

    @property
    def primary_reference(self):
        return self.references[0]

    def run(self):
        regions = []
        region_num = 0
        is_last_region = False

        while True:
            region_num += 1

            # Get the previous region in each pool
            prev_pair = (regions[-1].candidate_pairs[0]
                         if len(regions) >= 1 else None)
            prev_pair_same_pool = (regions[-2].candidate_pairs[0]
                                   if len(regions) > 2 else None)

            # If there are two regions or more
            if prev_pair_same_pool:
                # Gap opened between -1 and -2 regions
                if prev_pair.left.start > prev_pair_same_pool.right.start:
                    # If gap, left primer cannot overlap with -1 region
                    left_primer_left_limit = prev_pair.left.end + 1
                else:
                    # Left primer cannot overlap -2 region
                    left_primer_left_limit = (prev_pair_same_pool.right.start
                                              + 1)
            # If there is more than one region
            elif prev_pair:
                # Left primer cannot overlap with -1 region or you don't move
                left_primer_left_limit = prev_pair.left.end + 1
            else:
                # Region one only limit is 0
                left_primer_left_limit = 0

            # Right start limit maintains the minimum_overlap
            left_primer_right_limit = (prev_pair.right.end
                                       - self.min_overlap - 1
                                       if prev_pair else self.max_gap)

            # Last region if less than one amplicon length remaining
            if prev_pair:
                if ((len(self.primary_reference) - prev_pair.right.end)
                        < self.amplicon_length):
                    is_last_region = True
                    logger.debug(
                        'Region {}: is last region'.format(region_num))

            # Log limits
            logger.debug('Region {}: forward primer limits {}:{}'.format(
                region_num, left_primer_left_limit, left_primer_right_limit))

            # Find primers or handle no suitable error
            try:
                region = self._find_primers(region_num, left_primer_left_limit,
                                            left_primer_right_limit,
                                            is_last_region)
                regions.append(region)
            except NoSuitablePrimers:
                logger.debug('Region {}: no suitable primer error'.format(
                    region_num))
                break

            # Handle the end
            if is_last_region:
                logger.debug('Region {}: ending normally'.format(region_num))
                break

            # Report scores and alignments
            for i in range(0, len(self.references)):
                # Don't display alignment to reference
                logger.debug(regions[-1].candidate_pairs[0].left.alignments[i]
                             .formatted_alignment)
            logger.debug('Identities for sorted left candidates: ' + ','.join(
                ['%.2f' % each.left.percent_identity
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Left start for sorted candidates: ' + ','.join(
                ['%i' % each.left.start
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Left end for sorted candidates: '
                         + ','.join(['%i' % each.left.end
                                     for each in regions[-1].candidate_pairs]))
            logger.debug('Left length for sorted candidates: ' + ','.join(
                ['%i' % each.left.length
                 for each in regions[-1].candidate_pairs]))

            for i in range(0, len(self.references)):
                logger.debug(regions[-1].candidate_pairs[0].right.alignments[i]
                             .formatted_alignment)
            logger.debug('Identities for sorted right candidates: ' + ','.join(
                ['%.2f' % each.right.percent_identity
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Right start for sorted candidates: ' + ','.join(
                ['%i' % each.right.start
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Right end for sorted candidates: ' + ','.join(
                ['%i' % each.right.end
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Right length for sorted candidates: ' + ','.join(
                ['%i' % each.right.length
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Totals for sorted pairs: ' + ','.join(
                ['%.2f' % each.mean_percent_identity
                 for each in regions[-1].candidate_pairs]))

            if len(regions) > 1:
                # Remember, results now include this one, so -2
                # is the other pool
                trimmed_overlap = (regions[-2].candidate_pairs[0].right.end
                                   - regions[-1].candidate_pairs[0].left.end
                                   - 1)
                logger.info(
                    "Region %i: highest scoring product %i:%i, length %i, "
                    "trimmed overlap %i" % (
                        region_num, regions[-1].candidate_pairs[0].left.start,
                        regions[-1].candidate_pairs[0].right.start,
                        regions[-1].candidate_pairs[0].product_length,
                        trimmed_overlap))
            else:
                logger.info(
                    "Region %i: highest scoring product %i:%i, length %i" % (
                        region_num, regions[-1].candidate_pairs[0].left.start,
                        regions[-1].candidate_pairs[0].right.start,
                        regions[-1].candidate_pairs[0].product_length))

        # Return regions
        self.regions = regions

    def _find_primers(self, region_num, left_primer_left_limit,
                      left_primer_right_limit, is_last_region):
        """
        Find primers for a given region.

        Return a list of Region objects containing candidate
        primer pairs sorted by mean percent identity against all references.
        """

        # Calculate where to slice the reference
        if region_num == 1:
            chunk_start = 0
            chunk_end = (int((1 + self.max_variation / 2)
                             * self.amplicon_length))
        elif is_last_region:
            # Last time work backwards
            chunk_start = (int(len(self.primary_reference)
                               - ((1 + self.max_variation / 2)
                               * self.amplicon_length)))
            chunk_end = len(self.primary_reference)
        else:
            # right limit
            # - min overlap
            # - diff max min product length
            # - max primer length
            chunk_start = (int(left_primer_right_limit
                               - (self.max_variation * self.amplicon_length)
                               - settings.global_args['PRIMER_MAX_SIZE']))
            chunk_end = (int(chunk_start
                             + ((1 + self.max_variation/2)
                                 * self.amplicon_length)))
        initial_chunk_start = chunk_start
        initial_chunk_end = chunk_end

        # Primer3 setup
        p3_global_args = settings.global_args
        p3_seq_args = settings.seq_args
        p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
            [int(self.amplicon_length * (1 - self.max_variation / 2)),
             int(self.amplicon_length * (1 + self.max_variation / 2))]]
        p3_global_args['PRIMER_NUM_RETURN'] = self.max_candidates

        # Run primer3 until primers are found
        hit_left_limit = False
        while True:
            # Slice primary reference
            seq = str(self.primary_reference.seq[chunk_start:chunk_end])
            p3_seq_args['SEQUENCE_TEMPLATE'] = seq
            p3_seq_args['SEQUENCE_INCLUDED_REGION'] = [0, len(seq) - 1]
            logger.debug("Region %i: reference chunk %i:%i, length %i" % (
                region_num, chunk_start, chunk_end, len(seq)))
            primer3_output = primer3.bindings.designPrimers(p3_seq_args,
                                                            p3_global_args)

            candidate_pairs = []

            for cand_num in range(self.max_candidates):
                lenkey = 'PRIMER_LEFT_%i' % (cand_num)
                left_name = '%s_%i_%s' % (self.prefix, region_num, 'LEFT')
                right_name = '%s_%i_%s' % (self.prefix, region_num, 'RIGHT')

                if lenkey not in primer3_output:
                    break

                left_seq = str(
                    primer3_output['PRIMER_LEFT_%i_SEQUENCE' % (cand_num)])
                right_seq = str(
                    primer3_output['PRIMER_RIGHT_%i_SEQUENCE' % (cand_num)])

                left_start = int(
                    primer3_output['PRIMER_LEFT_%i' % (cand_num)][0]
                    + chunk_start)
                right_start = int(
                    primer3_output['PRIMER_RIGHT_%i' % (cand_num)][0]
                    + chunk_start + 1)

                left = CandidatePrimer('LEFT', left_name, left_seq, left_start)
                right = CandidatePrimer('RIGHT', right_name, right_seq,
                                        right_start)

                candidate_pairs.append(CandidatePrimerPair(left, right))

            set_left = set(pair.left.seq for pair in candidate_pairs)
            set_right = set(pair.right.seq for pair in candidate_pairs)

            logger.info(
                "Region %i: current position returned %i left and %i "
                "right unique" % (region_num, len(set_left), len(set_right)))

            if len(set_left) > 2 and len(set_right) > 2:
                return Region(region_num, chunk_start, candidate_pairs,
                              self.references, self.prefix, self.max_alts)

            # Move right if first region or to open gap
            if region_num == 1 or hit_left_limit:
                logger.debug("Region %i: stepping right, position %i"
                             % (region_num, chunk_start))
                chunk_start += self.step_size
                chunk_end += self.step_size
                # Hit end of regerence
                if chunk_end > len(self.primary_reference):
                    logger.debug("Region %i: hit right limit %i"
                                 % (region_num, len(self.primary_reference)))
                    raise NoSuitablePrimers("No suitable primers in region")
            else:
                # Move left for all other regions
                logger.debug(
                    "Region %i: stepping left, position %i, limit %s"
                    % (region_num, chunk_start, left_primer_left_limit))
                chunk_start -= self.step_size
                chunk_end -= self.step_size
                if chunk_start <= left_primer_left_limit:
                    # Switch direction to open gap
                    logger.debug("Region %i: hit left limit" % (region_num))
                    chunk_start = initial_chunk_start
                    chunk_end = initial_chunk_end
                    hit_left_limit = True
