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

from operator import attrgetter

from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from primalscheme import config
from primalscheme.align import align_secondary_reference, FailedAlignmentError
from primalscheme.primer import (
    calc_heterdodimer,
    calc_tm,
    design_primers,
    primer_thermo_filter,
    Dimer,
    Direction,
)

logger = logging.getLogger("primalscheme")


class Region:
    """
    A region forming part of a tiling amplicon scheme.
    """

    def __init__(self, scheme, region_num, left_limit, slice_start, right_limit=None):
        """Init Region."""
        self.region_num = region_num
        self.scheme = scheme
        self.left_limit = left_limit
        self.slice_start = slice_start
        self.right_limit = right_limit or len(scheme.primary_ref.seq)

        self.initial_slice_start = self.slice_start
        self.pool = 1 if self.region_num % 2 == 1 else 2

        self.left = None
        self.right = None
        self.left_flank_msa = None
        self.right_flank_msa = None
        self.left_candidates = []
        self.right_candidates = []
        self.exhausted_left_stepping = False
        self.failed_aln_ref_ids = []

        logger.debug(
            f"Region {region_num}, pool {self.pool}, left_limit {left_limit}, "
            f"slice_start {slice_start}, right_limit {self.right_limit}"
        )

        # check bounds
        if slice_start < left_limit or self.slice_end > self.right_limit:
            raise SliceOutOfBoundsError("The window slice is out of bounds.")

    @property
    def slice_end(self):
        """Slice end position."""
        return self.slice_start + self.scheme.amplicon_size_max - 1

    @property
    def ref_slice(self):
        """Primary reference sequence for the slice."""
        return self.scheme.primary_ref[self.slice_start : self.slice_end + 1]

    @property
    def flank_size(self):
        """Size of the slice flanks, where possible primers could be located."""
        return (
            int((self.scheme.amplicon_size_max - self.scheme.amplicon_size_min) / 2)
            + config.PRIMER_SIZE_RANGE.max
        )

    @property
    def left_flank(self):
        """Primary reference sequence for the left flank."""
        return self.ref_slice[: self.flank_size]

    @property
    def right_flank(self):
        """Primary reference sequence for the right flank."""
        return self.ref_slice[-self.flank_size :]

    @property
    def product_size(self):
        """Product size for the picked left/right primers."""
        if self.left and self.right:
            return self.right.start - self.left.start + 1
        return None

    def step_left(self, reason=None):
        """Step the slice start position left; raise if limit reached."""
        distance = config.STEP_DISTANCE
        if (self.slice_start - distance) < self.left_limit:
            logger.debug("Left limit reached")
            raise SliceOutOfBoundsError("Left window limit reached.")
        self.slice_start -= distance
        self.left_flank_msa, self.right_flank_msa = None, None
        logger.debug(
            f"Stepping left to {self.slice_start}"
            f"{' due to ' + reason if reason else ''}"
        )

    def step_right(self, reason=None):
        """Step the slice start position right; raise if limit reached."""
        distance = config.STEP_DISTANCE
        if (self.slice_end + distance) > self.right_limit:
            logger.debug("Right limit reached")
            raise SliceOutOfBoundsError("Right window limit reached.")
        self.slice_start += distance
        self.left_flank_msa, self.right_flank_msa = None, None
        logger.debug(
            f"Stepping right to {self.slice_start}"
            f"{' due to ' + reason if reason else ''}"
        )

    def reset_slice(self):
        """Reset the slice start position to its initial value."""
        self.slice_start = self.initial_slice_start

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

    def find_primers(self):
        """Find primers for this region."""
        while True:
            try:
                self._find_primers_for_slice()
                return
            except FailedAlignmentError as exc:
                if self.exhausted_left_stepping:
                    self._handle_failed_alignments(exc)
                self._try_stepping(reason="failed alignment")
            except NoSuitablePrimersError:
                self._try_stepping(reason="no suitable primers")

    def _align_flanks(self):
        """Perform multiple seq alignment of all refs for flanks."""
        self.left_flank_msa = MultipleSeqAlignment([self.left_flank])
        right_rev = SeqRecord(
            self.right_flank.reverse_complement().seq, id=self.scheme.primary_ref.id
        )
        self.right_flank_msa = MultipleSeqAlignment([right_rev])
        failed_ref_ids = []
        for ref in self.scheme.secondary_refs:
            try:
                self.left_flank_msa.append(
                    align_secondary_reference(self.left_flank, ref)
                )
                ra = align_secondary_reference(self.right_flank, ref)
                self.right_flank_msa.append(
                    SeqRecord(ra.reverse_complement().seq, id=ra.id)
                )
            except FailedAlignmentError:
                failed_ref_ids.append(ref.id)
        if failed_ref_ids:
            raise FailedAlignmentError(
                "Failed alignment between primary flank and secondary references.",
                ref_ids=failed_ref_ids,
            )

    def _try_stepping(self, reason=None):
        """Step region left, or right (from initial) on reaching left limit."""
        if self.exhausted_left_stepping:
            try:
                self.step_right(reason=reason)
            except SliceOutOfBoundsError:
                raise NoSuitablePrimersError("Right limit reached.")
        else:
            try:
                self.step_left(reason=reason)
            except SliceOutOfBoundsError:
                self.reset_slice()
                self.exhausted_left_stepping = True

    def _find_primers_for_slice(self):
        """Try to find suitable primers for the current slice."""
        logger.debug(
            f"Finding primers for slice {self.slice_start} to {self.slice_end}"
        )

        # Align flanks at this position
        self._align_flanks()

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

        # Check mismatch threshold & max penalty
        passing_candidates = []
        for primer in filtered_candidates:
            if (
                max(primer.mismatch_counts) <= config.PRIMER_MAX_MISMATCHES
                and primer.base_penalty <= config.PRIMER_MAX_BASE_PENALTY
            ):
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

    def _sort_candidates(self):
        """
        Sort candidates by penalty, start (higher), seq.
        seq is necessary to maintain deterministic output.
        """
        for candidates in [self.left_candidates, self.right_candidates]:
            candidates.sort(key=attrgetter("seq"))
            candidates.sort(key=attrgetter("start"), reverse=True)
            candidates.sort(key=attrgetter("combined_penalty"))

    def _pick_pair(self):
        """Pick the best scoring left and right primer for this region."""
        self._sort_candidates()

        no_suitable_exc = None
        try:
            self.left = self._pick_candidate(self.left_candidates)
            self.right = self._pick_candidate(self.right_candidates)
        except NoSuitablePrimersError as e:
            no_suitable_exc = e

        if logger.level >= logging.DEBUG:
            for direction in Direction:
                self._log_debug(direction)

        if no_suitable_exc:
            raise no_suitable_exc

    def _pick_candidate(self, candidates):
        """
        Pick the best scoring candidate that passes interaction check.
        """
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
        an existing primer in the same pool (or itself).
        """

        same_pool_primers = self.scheme.primers_in_pool(self.pool)

        for existing in same_pool_primers + [candidate]:  # also check for self-dimer
            thermo_het = calc_heterdodimer(candidate.seq, existing.seq)

            if thermo_het.structure_found:
                # Only exact 3' matches
                if (
                    len(thermo_het.ascii_structure_lines[0].rstrip("-").split()) == 2
                ) or (
                    len(thermo_het.ascii_structure_lines[3].rstrip("-").split()) == 2
                ):
                    # Slice overlap sequence
                    ol = thermo_het.ascii_structure_lines[1].split()[1]
                    ol_tm = calc_tm(ol)
                    if ol_tm > config.PRIMER_MAX_HETERODIMER_TH:
                        candidate.flagged_dimer = Dimer(existing, ol_tm)
                        return True
        return False

    def _handle_failed_alignments(self, exc):
        self.failed_aln_ref_ids.extend(exc.ref_ids)
        unique_ids = list(set(self.failed_aln_ref_ids))
        exclude_ids = list(
            filter(
                lambda x: self.failed_aln_ref_ids.count(x)
                >= self.scheme._max_failed_aln,
                unique_ids,
            )
        )
        if exclude_ids:
            self.scheme.exclude_references(exclude_ids)
            self.reset_slice()
        else:
            self._try_stepping()

    def _log_debug(self, direction):
        """Log detailed debug info for this region."""
        if direction == Direction.LEFT:
            candidates = self.left_candidates
            picked = self.left
        else:
            candidates = self.right_candidates
            picked = self.right

        if picked:
            logger.debug(f"Picked MSA slice: {picked.annotated_msa}")
        for i, primer in enumerate(candidates):
            max_mis = max(primer.mismatch_counts)
            total_mis = sum(primer.mismatch_counts)
            dimer_str = ""
            if primer.flagged_dimer:
                tm_str = f"{primer.flagged_dimer[1]:.2f}"
                if primer.flagged_dimer.existing == primer:
                    dimer_str = f"*homodimer tm {tm_str}, "
                else:
                    dimer_str = f"*heterodimer tm {tm_str}, "
            logger.debug(
                f"{direction.name} candidate {i}: {dimer_str}"
                f"{'*picked, ' if primer == picked else ''}"
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
