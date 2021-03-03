"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains MultiplexScheme.

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

from collections import namedtuple

from primalscheme import config
from primalscheme.align import FailedAlignmentError
from primalscheme.region import Region, NoSuitablePrimersError

logger = logging.getLogger("primalscheme")


Insert = namedtuple("Insert", "start end pool")


class MultiplexPanel:
    """
    A panel of one or more schemes, where the primers are checked for interactions
    between schemes
    """

    def __init__(
        self,
        reference_collections,
        prefix=config.PREFIX,
        amplicon_size_min=config.AMPLICON_SIZE_MIN,
        amplicon_size_max=config.AMPLICON_SIZE_MAX,
        target_overlap=config.TARGET_OVERLAP,
        primary_only=False,
        high_gc=False,
        progress_tracker=None,
    ):
        """Init MultiplexScheme."""

        self.reference_collections = reference_collections
        self.prefix = prefix
        self.amplicon_size_min = amplicon_size_min
        self.amplicon_size_max = amplicon_size_max
        self.target_overlap = target_overlap
        self.primary_only = primary_only
        self.high_gc = high_gc
        self.progress_tracker = progress_tracker
        self.combined_ref_len = sum(len(refs[0]) for refs in reference_collections)

        # config switching
        config_key = "HIGH_GC" if high_gc else "DEFAULT"
        config.PRIMER_SIZE_RANGE = config.PRIMER_SIZE_RANGES[config_key]
        config.PRIMER_GC_RANGE = config.PRIMER_GC_RANGES[config_key]

        # progress & logging setup
        if self.progress_tracker:
            self.combined_progress = 0
            self.progress_tracker.end = self.combined_ref_len

        # schemes
        self.schemes = []

        logger.debug(str(self))

    def __str__(self):
        lines = [
            "MuliplexPanel",
            f"Prefix: {self.prefix}",
            f"Amplicon size range: {self.amplicon_size_min} - {self.amplicon_size_max}",
            f"Target overlap: {self.target_overlap}",
            f"Primary only (pinned): {self.primary_only}",
        ]
        return "\n".join(lines)

    @property
    def primers(self):
        """Return a list of all primers in all schemes"""
        primers = []
        for scheme in self.schemes:
            primers.extend(scheme.primers)
        return primers

    def primers_in_pool(self, pool):
        """Return a list of all primers a pool across all schemes."""
        primers = []
        for scheme in self.schemes:
            primers.extend(scheme.primers_in_pool(pool))
        return primers

    @property
    def regions(self):
        regions = []
        for scheme in self.schemes:
            regions.extend(scheme.regions)
        return regions

    @property
    def gap_count(self):
        return sum(scheme.gap_count for scheme in self.schemes)

    @property
    def percent_coverage(self):
        return 100

    def design_panel(self):
        """Design a panel of multiplex primer schemes."""
        for references in self.reference_collections:
            scheme = MultiplexScheme(references, self)
            try:
                scheme.design_scheme()
            except NoSuitablePrimersError as e:
                if len(self.reference_collections) == 1:
                    raise e
            self.schemes.append(scheme)

            if self.progress_tracker:
                self.combined_progress += scheme.ref_len

        if self.progress_tracker:
            self.progress_tracker.goto(self.progress_tracker.end)
            self.progress_tracker.finish()


class MultiplexScheme:
    """A multiplex primer scheme, as part of a panel."""

    def __init__(self, references, panel):
        """Init MultiplexScheme."""
        self.regions = []
        self.panel = panel
        self.references = references
        self.primary_ref = references[0]
        self.secondary_refs = references[1:]
        self.excluded_refs = []
        self.ref_len = len(self.primary_ref)
        self._max_failed_aln = int(
            config.MAX_ALN_GAP_PERCENT * self.ref_len / config.STEP_DISTANCE
        )
        self.considered = 0
        self.no_suitable_primers = False

        logger.debug(str(self))

    def __str__(self):
        lines = [
            "MuliplexScheme",
            "References:",
        ]
        lines.extend(f" - {ref.id}" for ref in self.references)
        return "\n".join(lines)

    @property
    def region_count(self):
        """The number of regions in the scheme."""
        return len(self.regions)

    @property
    def _region_num(self):
        """The next/pending region number."""
        return len(self.regions) + 1

    @property
    def _prev(self):
        """The previous region."""
        return self.regions[-1] if self.regions else None

    @property
    def _prev_in_pool(self):
        """The previous region in the same pool as the next."""
        return self.regions[-2] if self.region_count > 1 else None

    @property
    def _remaining_distance(self):
        """The distance between the last region and the primary reference end."""
        if self.regions:
            return self.ref_len - self._prev.right.start
        return self.ref_len

    @property
    def _is_last_region(self):
        """Is there space for only one more region?"""
        return self._remaining_distance <= self.panel.amplicon_size_max

    @property
    def _left_limit(self):
        """The left limit for the next region."""
        if self._region_num == 1:
            return 0
        elif self._region_num == 2 or (
            self._prev.left.start > self._prev_in_pool.right.start
        ):
            # first region in second pool, or we have a gap. Force progress.
            return self._prev.left.end + 1
        else:
            # default case (just don't crash into prev in pool)
            return self._prev_in_pool.right.start + 1

    @property
    def _slice_start(self):
        """The ideal slice start position for the next region."""

        if self._region_num == 1:
            return 0

        insert_start = self._prev.right.end - self.panel.target_overlap - 1
        desired_slice_start = (
            insert_start
            - config.PRIMER_SIZE_RANGE.max
            - (self.panel.amplicon_size_max - self.panel.amplicon_size_min)
        )

        # if target overlap is impossible, take left_limit
        slice_start = max(desired_slice_start, self._left_limit)

        if self._is_last_region:
            # if forced to choose, take a gap at the end
            right_aligned = self.ref_len - self.panel.amplicon_size_max
            return min(slice_start, right_aligned)

        return slice_start

    @property
    def primers(self):
        """Return a list of all primers in the scheme."""
        primers = []
        for region in [region for region in self.regions]:
            primers.append(region.left)
            primers.append(region.right)
        return primers

    def primers_in_pool(self, pool):
        """Return a list of all primers in a pool."""
        return [primer for primer in self.primers if primer.pool == pool]

    def design_scheme(self):
        """Design a multiplex primer scheme."""
        is_last = False
        while not is_last:
            is_last = self._is_last_region
            region = Region(self, self._region_num, self._left_limit, self._slice_start)
            try:
                region.find_primers()
            except NoSuitablePrimersError as exc:
                if self._region_num == 1:
                    self.no_suitable_primers = True
                    raise exc
                break
            except FailedAlignmentError as exc:
                self._exclude_reference(exc.reference)
                continue

            self.regions.append(region)

            if self.panel.progress_tracker:
                self.panel.progress_tracker.region_num = self._region_num
                self.panel.progress_tracker.goto(
                    region.right.start + self.panel.combined_progress
                )

            if self._is_last_region and self._slice_start < self._left_limit:
                # additional region no longer possible due to small ref, relative to
                # amplicon size, or significant right stepping of prev region.
                break

    def _exclude_reference(self, reference):
        """Exclude a secondary reference."""
        for i, ref in enumerate(self.secondary_refs):
            if ref.id == reference.id:
                self.excluded_refs.append(self.secondary_refs.pop(i))
        self.panel.progress_tracker.interrupt()
        logger.info(f"Excluding reference {reference.id}: unable to align.")

    def add_considered(self, num):
        """Increment number of considered primers; update progress tracker."""
        self.considered += num
        if self.panel.progress_tracker:
            self.panel.progress_tracker.considered += self.considered

    @property
    def inserts(self):
        """A list of insert (start, end) tuples."""
        return [Insert(r.left.end + 1, r.right.end - 1, r.pool) for r in self.regions]

    @property
    def gap_count(self):
        """The total number of gaps in the scheme."""
        gaps = 0
        last_covered = self.inserts[0].end
        for insert in self.inserts[1:]:
            if insert.start > last_covered:
                gaps += 1
            last_covered = insert.end
        return gaps

    @property
    def percent_coverage(self):
        """Coverage %, with respect to the primary reference."""
        covered_coords = set(
            [x for insert in self.inserts for x in range(insert.start, insert.end + 1)]
        )
        return round(len(covered_coords) / self.ref_len * 100, 2)
