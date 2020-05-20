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

import logging
import os
import pickle

from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors

from primalscheme.multiplex import MultiplexScheme

logger = logging.getLogger("primalscheme")


class MultiplexReporter(MultiplexScheme):
    """Reporting methods to extend MultiplexScheme"""

    def write_all(self, path):
        self.write_bed(path)
        self.write_pickle(path)
        self.write_tsv(path)
        self.write_refs(path)
        self.write_schemadelica_plot(path)

    @property
    def inserts(self):
        return [(r.top_pair.left.end + 1, r.top_pair.right.end) for r in self.regions]

    @property
    def gap_count(self):
        gaps = 0
        last_covered = self.inserts[0][1]
        for insert in self.inserts:
            if insert[0] > last_covered:
                gaps += 1
            last_covered = insert[1]
        return gaps

    @property
    def percent_coverage(self):
        covered_coords = set([x for insert in self.inserts for x in range(*insert)])
        return round(len(covered_coords) / self.ref_len * 100, 2)

    def write_bed(self, path):
        """
        Write BED format to file.
        """
        filepath = path / f"{self.prefix}.scheme.bed"
        logger.info(f"Writing {filepath}")

        with open(filepath, "w") as bedhandle:
            for r in self.regions:
                print(
                    *map(
                        str,
                        [
                            self.primary_ref.id,
                            r.top_pair.left.start,
                            r.top_pair.left.end + 1,  # BED end is 1-based
                            r.top_pair.left.name,
                            r.pool,
                            "+",
                        ],
                    ),
                    sep="\t",
                    file=bedhandle,
                )
                print(
                    *map(
                        str,
                        [
                            self.primary_ref.id,
                            r.top_pair.right.end,
                            r.top_pair.right.start + 1,
                            r.top_pair.right.name,
                            r.pool,
                            "-",
                        ],
                    ),
                    sep="\t",
                    file=bedhandle,
                )

    def write_tsv(self, path):
        filepath = path / f"{self.prefix}.tsv"
        logger.info(f"Writing {filepath}")
        with open(filepath, "w") as tsvhandle:
            print(
                *["name", "pool", "seq", "size", "%gc", "tm (use 65)"],
                sep="\t",
                file=tsvhandle,
            )
            for r in self.regions:
                left = r.top_pair.left
                right = r.top_pair.right
                print(
                    *map(
                        str,
                        [
                            left.name,
                            r.pool,
                            left.seq,
                            left.size,
                            "%.2f" % left.gc,
                            "%.2f" % left.tm,
                        ],
                    ),
                    sep="\t",
                    file=tsvhandle,
                )
                print(
                    *map(
                        str,
                        [
                            right.name,
                            r.pool,
                            right.seq,
                            right.size,
                            "%.2f" % right.gc,
                            "%.2f" % right.tm,
                        ],
                    ),
                    sep="\t",
                    file=tsvhandle,
                )

    def write_pickle(self, path):
        filepath = path / f"{self.prefix}.pickle"
        logger.info(f"Writing {filepath}")
        with open(filepath, "wb") as pickleobj:
            pickle.dump(self, pickleobj)

    def write_refs(self, path):
        filepath = path / f"{self.prefix}.reference.fasta"
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

    def calc_gc(self, sequence):
        """
        Return gc content as fraction
        """
        g = sequence.count("G") + sequence.count("g")
        c = sequence.count("C") + sequence.count("c")
        return (g + c) / len(sequence)

    def write_schemadelica_plot(self, path):
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
                FeatureLocation(
                    r.top_pair.left.start, r.top_pair.left.end, strand=strand
                )
            )
            rev_feature = SeqFeature(
                FeatureLocation(
                    r.top_pair.right.end, r.top_pair.right.start, strand=strand
                )
            )
            region_feature = SeqFeature(
                FeatureLocation(
                    r.top_pair.left.start, r.top_pair.right.start, strand=strand
                )
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

        pdf_filepath = os.path.join(path, "{}.pdf".format(self.prefix))
        svg_filepath = os.path.join(path, "{}.svg".format(self.prefix))
        logger.info(f"Writing {pdf_filepath}")
        logger.info(f"Writing {svg_filepath}")
        gd_diagram.write(pdf_filepath, "PDF", dpi=300)
        gd_diagram.write(svg_filepath, "SVG", dpi=300)
