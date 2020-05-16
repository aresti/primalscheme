"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains supporting scheme components.

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
import parasail

from Bio.Seq import Seq
from primer3 import calcTm

logger = logging.getLogger("primalscheme")


class Primer(object):
    """A primer."""

    def __init__(self, seq, start, direction):
        self.direction = direction
        self.seq = seq
        self.start = start

    def __str__(self):
        return f"{self.direction}:{self.seq}:{self.start}"

    @property
    def length(self):
        return len(self.seq)


class CandidatePrimer(Primer):
    """A candidate primer."""

    def __init__(self, seq, start, direction, name="", penalty=None):
        super().__init__(seq, start, direction)
        self.name = name
        self.penalty = penalty
        self.identity = 0
        self.tm = calcTm(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)
        self.gc = 100.0 * (seq.count("G") + seq.count("C")) / len(seq)
        self.alignments = []

    def align(self, references):
        for ref in references[1:]:
            # align against non-primary references
            alignment = get_alignment(self, ref)
            if alignment:
                self.alignments.append(alignment)

        # Calculate average percent identity
        if self.alignments:
            scores = [a[0] for a in self.alignments]
            self.identity = sum(scores) / len(scores)
            return

    @property
    def end(self):
        if self.direction == "LEFT":
            return self.start + self.length
        else:
            return self.start - self.length


class CandidatePrimerPair(object):
    """A pair of candidate primers."""

    def __init__(self, left, right):
        self.left = left
        self.right = right

    @property
    def mean_identity(self):
        return (self.left.identity + self.right.identity) / 2

    @property
    def product_length(self):
        return self.right.start - self.left.start


def get_alignment(primer, reference):
    """An seqan alignment of a primer against a reference."""

    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 2
    EXTEND = 1
    IDENTITY_THRESHOLD = 0.7

    if primer.direction == "LEFT":
        query = primer.seq
    elif primer.direction == "RIGHT":
        query = str(Seq(primer.seq).reverse_complement())

    # Semi-Global, do not penalize gaps at beginning and end of s2/database
    trace = parasail.sg_dx_trace_striped_sat(
        query, str(reference.seq), OPEN, EXTEND, MATRIX
    )
    traceback = trace.get_traceback()

    query_end = trace.end_query
    ref_end = trace.end_ref

    # Get alignment strings
    aln_query = traceback.query[ref_end - query_end : ref_end + 1]
    cigar = traceback.comp[ref_end - query_end : ref_end + 1]
    aln_ref = traceback.ref[ref_end - query_end : ref_end + 1]

    # Identity for glocal alignment
    identity = cigar.count("|") / len(cigar)

    # Alignment failed
    if identity < IDENTITY_THRESHOLD:
        return None

    # Format alignment
    refid = reference.id[:30]
    name = primer.name[:30]
    formatted_query = f"\n{name: <30} {1: >6} {aln_query} {query_end}"
    if primer.direction == "LEFT":
        formatted_ref = f"\n{refid: <30} {ref_end - query_end: >6} {aln_ref} {ref_end}"
    elif primer.direction == "RIGHT":
        formatted_ref = f"\n{refid: <30} {ref_end: >6} {aln_ref} {ref_end - query_end}"

    formatted_alignment = (
        formatted_query + f"\n{'': <30} {'': >6} {cigar}" + formatted_ref
    )

    del trace

    return (identity, formatted_alignment)
