"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the Window and Region classes.

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

from operator import attrgetter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from primalscheme import config
from primalscheme.align import align_secondary_reference, FailedAlignmentError
from primalscheme.primer import design_primers, primer_thermo_filter, Direction

logger = logging.getLogger("primalscheme")


class Window:
    """
    A sliding window representing slice coordinates
    against the primary reference.
    """

    def __init__(self, scheme, left_limit, slice_start, right_limit=None):
        """Init Window."""
        self.scheme = scheme
        self.left_limit = left_limit
        self.slice_start = slice_start
        self.right_limit = right_limit or len(scheme.primary_ref.seq)

        self._initial_slice_start = self.slice_start
        self.left_flank_msa = None
        self.right_flank_msa = None

        logger.debug(
            f"Window: left_limit {left_limit}, slice_start {slice_start}, "
            f"right_limit {self.right_limit}"
        )

        # check bounds
        if slice_start < left_limit or self.slice_end > self.right_limit:
            raise SliceOutOfBoundsError("The window slice is out of bounds.")

    def step_left(self):
        """Step the slice start position left, raise if limit reached"""
        distance = config.STEP_DISTANCE
        if (self.slice_start - distance) < self.left_limit:
            logger.debug("Left limit reached")
            raise SliceOutOfBoundsError("Left window limit reached.")
        self.slice_start -= distance
        self.left_flank_msa, self.right_flank_msa = None, None
        logger.debug(f"Stepping left to {self.slice_start}")

    def step_right(self):
        """Step the slice start position right, raise if limit reached"""
        distance = config.STEP_DISTANCE
        if (self.slice_end + distance) > self.right_limit:
            logger.debug("Right limit reached")
            raise SliceOutOfBoundsError("Right window limit reached.")
        self.slice_start += distance
        self.left_flank_msa, self.right_flank_msa = None, None
        logger.debug(f"Stepping right to {self.slice_start}")

    def reset_slice(self):
        """Reset the slice start position to its initial value"""
        self.slice_start = self._initial_slice_start

    @property
    def slice_end(self):
        """The slice end position"""
        return self.slice_start + self.scheme.amplicon_size_max - 1

    @property
    def ref_slice(self):
        """The reference sequence for the slice"""
        return self.scheme.primary_ref[self.slice_start : self.slice_end + 1]

    @property
    def flank_size(self):
        """The size of the slice flanks, where possible primers could be located"""
        return (
            int((self.scheme.amplicon_size_max - self.scheme.amplicon_size_min) / 2)
            + config.PRIMER_SIZE_RANGE.max
        )

    @property
    def left_flank(self):
        """The reference sequence for the left flank"""
        return self.ref_slice[: self.flank_size]

    @property
    def right_flank(self):
        """The reference sequence for the right flank"""
        return self.ref_slice[-self.flank_size :]

    def align_flanks(self):
        """Perform multiple seq alignment of all refs for flanks."""
        self.left_flank_msa = MultipleSeqAlignment([self.left_flank])
        right_rev = SeqRecord(
            self.right_flank.reverse_complement().seq, id=self.scheme.primary_ref.id
        )
        self.right_flank_msa = MultipleSeqAlignment([right_rev])
        for ref in self.scheme.secondary_refs:
            self.left_flank_msa.append(align_secondary_reference(self.left_flank, ref))
            ra = align_secondary_reference(self.right_flank, ref)
            self.right_flank_msa.append(
                SeqRecord(ra.reverse_complement().seq, id=ra.id)
            )


class Region(Window):
    """
    A region forming part of a tiling amplicon scheme.
    """

    def __init__(self, region_num, *args, **kwargs):
        """Init Region."""
        self.region_num = region_num
        self.pool = 1 if self.region_num % 2 == 1 else 2
        super().__init__(*args, **kwargs)

        self.left = None
        self.right = None
        self.left_candidates = []
        self.right_candidates = []
        self.exhausted_left_stepping = False
        self.failed_alignment_count = 0

        logger.debug(f"Region {region_num}, pool {self.pool}")

    @property
    def product_size(self):
        """The product size for the picked left/right primers"""
        if self.left and self.right:
            return self.right.start - self.left.start + 1
        return None

    def find_primers(self):
        """Find primers for this region"""
        while True:
            try:
                self._find_primers_for_slice()
                return
            except FailedAlignmentError as exc:
                self.failed_alignment_count += 1
                if self.failed_alignment_count > self.scheme._max_failed_aln:
                    raise exc
                self._try_stepping()
            except NoSuitablePrimersError:
                self._try_stepping()

    def _try_stepping(self):
        """Step region left, or right (from initial start) once left limit reached"""
        if self.exhausted_left_stepping:
            try:
                self.step_right()
            except SliceOutOfBoundsError:
                raise NoSuitablePrimersError("Right limit reached.")
        else:
            try:
                self.step_left()
            except SliceOutOfBoundsError:
                self.reset_slice()
                self.exhausted_left_stepping = True

    def _find_primers_for_slice(self):
        """Try to find suitable primers for the current slice"""
        logger.debug(
            f"Finding primers for slice {self.slice_start} to {self.slice_end}"
        )

        # Align flanks at this position
        self.align_flanks()

        # Design primers for the left and right flanks
        candidates = design_primers(
            self.left_flank_msa,
            Direction.LEFT,
            self.pool,
            offset=self.slice_start,
            primary_only=self.scheme.primary_only,
        )
        candidates.extend(
            design_primers(
                self.right_flank_msa,
                Direction.RIGHT,
                self.pool,
                offset=self.slice_end - self.flank_size + 1,
                primary_only=self.scheme.primary_only,
            )
        )

        # Update considered (trigger progress update)
        self.scheme.add_considered(len(candidates))

        # Hard thermo filter
        filtered_candidates = [
            primer for primer in candidates if primer_thermo_filter(primer)
        ]

        # Check mismatch threshold
        passing_candidates = []
        for primer in filtered_candidates:
            if max(primer.mismatch_counts) <= config.PRIMER_MAX_MISMATCHES:
                passing_candidates.append(primer)
        candidates = passing_candidates

        # Pull out left and right from passing candidates
        self.left_candidates = [
            primer for primer in candidates if primer.direction == Direction.LEFT
        ]
        self.right_candidates = [
            primer for primer in candidates if primer.direction == Direction.RIGHT
        ]

        if not (self.left_candidates and self.right_candidates):
            raise NoSuitablePrimersError(
                "Unable to find at least one suitable pair of candidate primers."
            )

        # Pick best-scoring left and right candidates
        self._pick_pair()

    def _sort_candidates(self, candidates):
        """
        Sort candidates by penalty, start (higher), seq.
        seq is necessary to maintain deterministic output.
        """
        candidates.sort(key=attrgetter("seq"))
        candidates.sort(key=attrgetter("start"), reverse=True)
        candidates.sort(key=attrgetter("combined_penalty"))

    def _pick_pair(self):
        """Pick the best scoring left and right primer for the region"""
        self._sort_candidates(self.left_candidates)
        self._sort_candidates(self.right_candidates)
        self.left = self._pick_candidate(self.left_candidates)
        self.right = self._pick_candidate(self.right_candidates)

        if logger.level >= logging.DEBUG:
            self._log_debug(Direction.LEFT)
            self._log_debug(Direction.RIGHT)

    def _pick_candidate(self, candidates):
        """Pick the best scoring candidate that passes a same-pool heterodimer check"""
        for candidate in candidates:
            if not self._check_for_heterodimers(candidate):
                return candidate
        raise NoSuitablePrimersError(
            "All candidates form stable heterodimers with existing primers "
            "in this pool."
        )

    def _check_for_heterodimers(self, candidate):
        """
        Return True if candidate primer forms stable heterodimer with
        an existing primer in the same pool.
        """
        for existing in self.scheme.primers_in_pool(self.pool):
            thermo_end = primer3.bindings.calcEndStability(
                candidate.seq,
                existing.seq,
                mv_conc=config.MV_CONC,
                dv_conc=config.DV_CONC,
                dna_conc=config.DNA_CONC,
                dntp_conc=config.DNTP_CONC,
                temp_c=config.TEMP_C,
            )
            if thermo_end.dg / 1000 < config.HETERODIMER_DG_THRESHOLD:
                logger.debug(
                    f"Primer interaction between {candidate.seq} and {existing.seq} "
                    f"predicted with a ∆G of {thermo_end.dg / 1000:.2f} kcal/mol"
                )
                candidate.interacts_with = existing
                return True
        return False

    def _log_debug(self, direction):
        """Log detailed debug info for the region"""
        if direction == Direction.LEFT:
            logger.debug(f"Left region flank MSA: {self.left_flank_msa}")
            candidates = self.left_candidates
            picked = self.left
        else:
            logger.debug(f"Right region flank MSA: {self.right_flank_msa}")
            candidates = self.right_candidates
            picked = self.right

        logger.debug(f"Picked primer {picked}")
        logger.debug(f"Picked MSA slice: {picked.annotated_msa}")
        for i, primer in enumerate(candidates[:10]):
            max_mis = max(primer.mismatch_counts)
            total_mis = sum(primer.mismatch_counts)
            logger.debug(
                f"{direction.name} candidate {i}: "
                f"{'*interaction, ' if primer.interacts_with else ''}"
                f"base pen {primer.base_penalty:.3f}, "
                f"combined pen {primer.combined_penalty:.3f}, "
                f"{max_mis} max mismatch{'' if max_mis == 1 else 'es'}, "
                f"{total_mis} total mismatch{'' if total_mis == 1 else 'es'}."
            )


class NoSuitablePrimersError(Exception):
    """Unable to find suitable primers at the current cursor position."""

    pass


class SliceOutOfBoundsError(Exception):
    """The requested start position would put the slice out of bounds."""

    pass
