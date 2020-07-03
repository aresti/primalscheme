"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the MutilplexReporter object.
This object extends a MultiplexScheme object to provide reporting methods.

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

import csv
import json
import logging
import pickle

from abc import ABC, abstractmethod
from collections import namedtuple

from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors

from primalscheme import config, __version__ as version
from primalscheme.multiplex import MultiplexScheme
from primalscheme.primer import Direction

logger = logging.getLogger("primalscheme")


Insert = namedtuple("Insert", "start end pool")


class MultiplexReporter(MultiplexScheme):
    """Reporting methods to extend MultiplexScheme."""

    def __init__(self, outpath, *args, **kwargs):
        """Init MultiplexReporter."""
        self.outpath = outpath
        super().__init__(*args, **kwargs)

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

    def calc_gc(self, sequence):
        """Return gc content for a sequence as fraction."""
        g = sequence.count("G") + sequence.count("g")
        c = sequence.count("C") + sequence.count("c")
        return (g + c) / len(sequence)

    def write_default_outputs(self):
        """Write all default output files."""
        self.write_primer_bed()
        self.write_insert_bed()
        self.write_primer_tsv()
        self.write_refs()
        self.write_schemadelica_plot()
        self.write_run_report_json()

    def write_primer_bed(self):
        """Write primer BED file."""
        filepath = self.outpath / f"{self.prefix}.primer.bed"
        logger.info(f"Writing {filepath}")

        rows = []
        ref_id = self.primary_ref.id

        for i, p in enumerate(self.primers):
            region_num = int(i / 2) + 1
            start = p.start if p.direction == Direction.LEFT else p.end
            end = p.end if p.direction == Direction.LEFT else p.start
            name = f"{self.prefix}_{region_num}_{p.direction.name}"
            rows.append(
                [
                    ref_id,
                    start,
                    end + 1,  # BED end is 1-based
                    name,
                    p.pool,
                    p.direction.value,
                ]
            )

        with open(filepath, "w") as fh:
            cw = csv.writer(fh, delimiter="\t")
            cw.writerows(rows)

    def write_insert_bed(self):
        """Write insert BED file."""
        filepath = self.outpath / f"{self.prefix}.insert.bed"
        logger.info(f"Writing {filepath}")

        rows = []
        ref_id = self.primary_ref.id

        for insert_num, insert in enumerate(self.inserts):
            rows.append(
                [
                    ref_id,
                    insert.start,
                    insert.end + 1,  # BED end is 1-based
                    f"{self.prefix}_INSERT_{insert_num + 1}",
                    insert.pool,
                    "+",
                ]
            )

        with open(filepath, "w") as fh:
            cw = csv.writer(fh, delimiter="\t")
            cw.writerows(rows)

    def write_primer_tsv(self):
        """Write primer TSV file."""
        filepath = self.outpath / f"{self.prefix}.primer.tsv"
        logger.info(f"Writing {filepath}")

        rows = [["name", "pool", "seq", "size", "%gc", "tm (use 65)"]]

        for i, p in enumerate(self.primers):
            region_num = int(i / 2) + 1
            rows.append(
                [
                    f"{self.prefix}_{region_num}_{p.direction.name}",
                    p.pool,
                    p.seq,
                    p.size,
                    f"{p.gc:.2f}",
                    f"{p.tm:.2f}",
                ]
            )

        with open(filepath, "w") as fh:
            cw = csv.writer(fh, delimiter="\t")
            cw.writerows(rows)

    def write_pickle(self):
        """Write pickle file."""
        filepath = self.outpath / f"{self.prefix}.pickle"
        logger.info(f"Writing {filepath}")
        with open(filepath, "wb") as pickleobj:
            pickle.dump(self, pickleobj)

    def write_refs(self):
        """Write reference FASTA."""
        filepath = self.outpath / f"{self.prefix}.reference.fasta"
        logger.info(f"Writing {filepath}")
        with open(filepath, "w"):
            SeqIO.write(self.references, filepath, "fasta")

    def apply_to_window(self, sequence, window_size, function, step=None):
        """
        Modified from
        https://github.com/biopython/biopython/blob/master/Tests/test_GenomeDiagram.py
        Apply function to windows of the given sequence.
        Returns a list of (position, value) tuples for fragments of the passed
        sequence of length window_size (stepped by step), calculated by the
        passed function.  Returned positions are the midpoint of each window.
        Arguments:
        - sequence - Bio.Seq.Seq object.
        - window_size - an integer describing the length of sequence
          to consider.
        - step - an integer describing the step to take between windows
          (default = window_size//2).
        - function - Method or function that accepts a Bio.Seq.Seq object
          as its sole argument and returns a single value.
        apply_to_window(sequence, window_size, function) ->
            [(int, float),(int, float),...]
        """
        seqlen = len(sequence)  # Total length of sequence to be used
        if step is None:  # Use half window-width or 1 if larger
            step = max(window_size // 2, 1)
        else:  # Use specified step, or 1 if greater
            step = max(step, 1)

        results = []  # Holds (position, value) results

        # Perform the passed function on as many windows as possible, short of
        # overrunning the sequence
        pos = 0
        while pos < seqlen - window_size + 1:
            # Obtain sequence fragment
            start = pos
            middle = (pos + window_size + pos) // 2
            end = pos + window_size
            fragment = sequence[start:end]
            # Apply function to the sequence fragment
            value = function(fragment)
            results.append((middle, value))  # Add results to list
            # Advance to next fragment
            pos += step

        # Use the last available window on the sequence, even if it means
        # re-covering old ground
        if pos != seqlen - window_size:
            # Obtain sequence fragment
            pos = seqlen - window_size
            start = pos
            middle = (pos + window_size + pos) // 2
            end = pos + window_size
            fragment = sequence[start:end]
            # Apply function to sequence fragment
            value = function(fragment)
            results.append((middle, value))  # Add results to list

        return results  # Return the list of (position, value) results

    def write_schemadelica_plot(self):
        """Write schemadelica plot as SVG and PDF."""
        gd_diagram = GenomeDiagram.Diagram("Primer Scheme", track_size=0.15)
        primer_feature_set = GenomeDiagram.FeatureSet()

        # make the gc track
        window = 50
        gc_set = GenomeDiagram.GraphSet("GC content")
        graphdata1 = self.apply_to_window(self.primary_ref.seq, window, self.calc_gc)
        gc_set.new_graph(
            graphdata1,
            "GC content",
            style="line",
            color=colors.violet,
            altcolor=colors.purple,
        )
        gc_track = GenomeDiagram.Track(
            "GC content", height=1.5, greytrack=0, scale_largetick_interval=1e3
        )
        gc_track.add_set(gc_set)

        # make the primer track
        for r in self.regions:
            region = str(r.region_num)
            strand = 1 if r.region_num % 2 else -1

            fwd_feature = SeqFeature(
                FeatureLocation(r.left.start, r.left.end, strand=strand)
            )
            rev_feature = SeqFeature(
                FeatureLocation(r.right.end, r.right.start, strand=strand)
            )
            region_feature = SeqFeature(
                FeatureLocation(r.left.start, r.right.start, strand=strand)
            )

            primer_color = colors.red
            region_color = colors.palevioletred

            primer_feature_set.add_feature(
                region_feature,
                color=region_color,
                name=region,
                label=True,
                label_position="middle",
                label_angle=0 if strand == 1 else -180,
            )
            primer_feature_set.add_feature(fwd_feature, color=primer_color, name=region)
            primer_feature_set.add_feature(rev_feature, color=primer_color, name=region)

        primer_track = GenomeDiagram.Track(name="Annotated Features", height=1)
        primer_track.add_set(primer_feature_set)

        gd_diagram.add_track(primer_track, 2)
        gd_diagram.add_track(gc_track, 1)

        rows = max(2, int(round(len(self.primary_ref) / 10000.0)))
        gd_diagram.draw(
            format="linear",
            pagesize=(300 * rows, 200 * rows),
            fragments=rows,
            start=0,
            end=len(self.primary_ref),
        )

        pdf_filepath = self.outpath / f"{self.prefix}.plot.pdf"
        svg_filepath = self.outpath / f"{self.prefix}.plot.svg"
        logger.info(f"Writing {pdf_filepath}")
        logger.info(f"Writing {svg_filepath}")
        gd_diagram.write(str(pdf_filepath), "PDF", dpi=300)
        gd_diagram.write(str(svg_filepath), "SVG", dpi=300)

    def write_run_report_json(self):
        """Write run report json."""
        filepath = self.outpath / f"{self.prefix}.report.json"
        logger.info(f"Writing {filepath}")
        data = {
            "references": [ref.id for ref in self.references],
            "primary_ref": self.references[0].id,
            "secondary_refs": [ref.id for ref in self.secondary_refs],
            "secondary_ref_count": len(self.secondary_refs),
            "excluded_refs": [ref.id for ref in self.excluded_refs],
            "excluded_ref_count": len(self.excluded_refs),
            "regions": len(self.regions),
            "percent_coverage": self.percent_coverage,
            "gaps": self.gap_count,
            "config": {
                "amplicon_size_min": self.amplicon_size_min,
                "amplicon_size_max": self.amplicon_size_max,
                "target_overlap": self.target_overlap,
                "high_gc_mode": self.high_gc,
                "step_distance": config.STEP_DISTANCE,
                "primer_size_range": {
                    "min": config.PRIMER_SIZE_RANGE.min,
                    "max": config.PRIMER_SIZE_RANGE.max,
                    "opt": config.PRIMER_SIZE_RANGE.opt,
                },
                "primer_gc_range": {
                    "min": config.PRIMER_GC_RANGE.min,
                    "max": config.PRIMER_GC_RANGE.max,
                    "opt": config.PRIMER_GC_RANGE.opt,
                },
                "primalscheme_version": version,
            },
        }
        filepath.write_text(json.dumps(data))


class ProgressTracker(ABC):
    """Abstract base class for ProgressTracker."""

    @abstractmethod
    def goto(self, val):
        """Update progress to val."""
        ...

    @property
    def end(self):
        """Progress end value."""
        ...

    @end.setter
    @abstractmethod
    def end(self, val):
        """Set progress end value."""
        ...

    @property
    def considered(self):
        """Count of considered primers."""
        ...

    @considered.setter
    @abstractmethod
    def considered(self, considered):
        """Set count of considered primers."""
        ...

    @property
    def region_num(self):
        """The current region num."""
        ...

    @region_num.setter
    @abstractmethod
    def region_num(self, region_num):
        """Set current region num."""
        ...

    @abstractmethod
    def interrupt(self):
        """Prepare to be interrupted by a log message."""
        ...

    @abstractmethod
    def finish(self):
        """Finish tracking progress."""
        ...
