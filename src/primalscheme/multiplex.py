"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains MultiplexScheme and supporting classes.

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

from primalscheme.wrapper import design_primers, InsufficientPrimersError
from primalscheme.components import (
    CandidatePrimer,
    CandidatePrimerPair,
    reversor,
)

logger = logging.getLogger("primalscheme")


class MultiplexScheme:
    """A complete multiplex primer scheme."""

    def __init__(
        self,
        references,
        p3_global,
        target_overlap,
        step_distance,
        min_unique,
        prefix,
        progress_func=None,
    ):

        self.references = references
        self.p3_global = p3_global
        self.target_overlap = target_overlap
        self.step_distance = step_distance
        self.min_unique = min_unique
        self.prefix = prefix
        self.progress_func = progress_func

        # derived
        self.primary_ref = references[0]
        self.ref_len = len(self.primary_ref)
        self.amplicon_size_min = p3_global["PRIMER_PRODUCT_SIZE_RANGE"][0][0]
        self.amplicon_size_max = p3_global["PRIMER_PRODUCT_SIZE_RANGE"][0][1]
        self.amplicon_max_variation = self.amplicon_size_max - self.amplicon_size_min
        self.primer_max_size = p3_global["PRIMER_MAX_SIZE"]

        self.regions = []

    def design_scheme(self):
        """Design a multiplex primer scheme"""

        regions = []
        region_num = 0
        is_last = False

        while not is_last:
            region_num += 1
            prev = regions[-1].top_pair if regions else None
            prev_in_pool = regions[-2].top_pair if region_num > 2 else None

            # determine left limit
            if region_num == 1:
                left_limit = 0
            elif region_num == 2 or prev.left.start > prev_in_pool.right.start:
                # first region in second pool, or we have a gap. Force progress.
                left_limit = prev.left.end + 1
            else:
                # default case (just don't crash into prev in pool)
                left_limit = prev_in_pool.right.start + 1

            # determine slice start
            if region_num == 1:
                slice_start = 0
            else:
                insert_start = prev.right.end - self.target_overlap - 1
                slice_start = (
                    insert_start - self.primer_max_size - self.amplicon_max_variation
                )

                # if target overlap is impossible, take left_limit
                slice_start = max(slice_start, left_limit)

                # handle last region
                remaining_distance = self.ref_len - prev.right.start
                is_last = remaining_distance <= self.amplicon_size_max
                if is_last:
                    right_aligned = self.ref_len - self.amplicon_size_max
                    # if forced to choose, take a gap at the end
                    slice_start = min(slice_start, right_aligned)

            # create region, find primers
            region = Region(region_num, self, left_limit, slice_start)
            try:
                region.find_primers()
            except NoSuitablePrimersError as e:
                if region_num == 1:
                    raise e  # we've got nothing
                break
            regions.append(region)
            if self.progress_func:
                self.progress_func(region.top_pair.right.start, self.ref_len)

        self.regions = regions
        if self.progress_func:
            self.progress_func(self.ref_len, self.ref_len)


class Window:
    """
    A sliding window representing slice coordinates
    against the primary reference.
    """

    def __init__(self, scheme, left_limit, slice_start, right_limit=None):
        self.scheme = scheme
        self.left_limit = left_limit
        self.slice_start = slice_start
        self.right_limit = right_limit or len(scheme.primary_ref.seq)

        self._initial_slice_start = self.slice_start

        logger.debug(
            f"Window: left_limit {left_limit}, slice_start {slice_start}, "
            f"right_limit {self.right_limit}"
        )

        # check bounds
        if slice_start < left_limit or self.slice_end > self.right_limit:
            raise SliceOutOfBoundsError("The window slice is out of bounds.")

    def step_left(self):
        distance = self.scheme.step_distance
        if (self.slice_start - distance) < self.left_limit:
            logger.debug("Left limit reached")
            raise SliceOutOfBoundsError("Left window limit reached.")
        self.slice_start -= distance
        logger.debug(f"Stepping left to {self.slice_start}")

    def step_right(self):
        distance = self.scheme.step_distance
        if (self.slice_end + distance) > self.right_limit:
            logger.debug("Right limit reached")
            raise SliceOutOfBoundsError("Right window limit reached.")
        self.slice_start += distance
        logger.debug(f"Stepping right to {self.slice_start}")

    def reset_slice(self):
        self.slice_start = self._initial_slice_start

    @property
    def slice_end(self):
        return self.slice_start + self.scheme.amplicon_size_max

    @property
    def ref_slice(self):
        primary_ref = self.scheme.primary_ref.seq  # Seq object
        return str(primary_ref[self.slice_start : self.slice_end])


class Region(Window):
    """
    Region forming part of a tiling amplicon scheme.
    """

    def __init__(self, region_num, *args, **kwargs):
        self.region_num = region_num
        self.pool = 1 if self.region_num % 2 == 1 else 2
        self.candidate_pairs = []
        self.top_pair = None
        super().__init__(*args, **kwargs)  # init Window

        logger.debug(f"Region {region_num}, pool {self.pool}")

    def find_primers(self):
        """
        Try <-> step window logic:
        """
        exhausted_left = False
        while True:
            try:
                self._find_primers_for_slice()
                break
            except InsufficientPrimersError:
                if exhausted_left:
                    try:
                        # step right, retry
                        self.step_right()
                    except SliceOutOfBoundsError:
                        raise NoSuitablePrimersError("Right limit reached.")
                else:
                    try:
                        # step left, retry
                        self.step_left()
                    except SliceOutOfBoundsError:
                        self.reset_slice()
                        exhausted_left = True
            except NoSuitablePrimersError as e:
                # out of options
                raise e

    def _find_primers_for_slice(self):
        """Try to find sufficient primers for the current slice"""

        logger.debug(f"Finding primers for slice [{self.slice_start}:{self.slice_end}]")

        pairs = design_primers(
            self.ref_slice,
            self.scheme.p3_global,
            self.scheme.min_unique,
            offset=self.slice_start,
        )

        name = f"{self.scheme.prefix}_{self.region_num}"
        for i in range(len(pairs[0])):
            left = CandidatePrimer(
                pairs[0][i].seq,
                pairs[0][i].start,
                CandidatePrimer.Direction.left,
                f"{name}_LEFT",
                penalty=pairs[0][i].penalty,
            )
            right = CandidatePrimer(
                pairs[1][i].seq,
                pairs[1][i].start,
                CandidatePrimer.Direction.right,
                f"{name}_RIGHT",
                penalty=pairs[1][i].penalty,
            )
            left.align(self.scheme.references)
            right.align(self.scheme.references)
            self.candidate_pairs.append(CandidatePrimerPair(left, right))

        self._pick_pair()

    def _sort_candidate_pairs(self):
        """Sort the list of candidate pairs in place"""

        self.candidate_pairs.sort(
            key=lambda x: (reversor(x.mean_identity), x.combined_penalty)
        )

    def _pick_pair(self):
        self._sort_candidate_pairs()
        self.top_pair = self.candidate_pairs[0]

        logger.debug(f"Picked: {self.top_pair}")
        for alignment in self.top_pair.left.alignments:
            logger.debug(alignment[1])
        for alignment in self.top_pair.right.alignments:
            logger.debug(alignment[1])
        for i, pair in enumerate(self.candidate_pairs):
            logger.debug(
                f"Candidate pair {i}: "
                f"{pair.mean_identity} identity, "
                f"penalty {pair.combined_penalty}"
            )


class NoSuitablePrimersError(Exception):
    """No suitable primers found."""

    pass


class SliceOutOfBoundsError(Exception):
    """The requested start position would put the slice out of bounds."""

    pass
