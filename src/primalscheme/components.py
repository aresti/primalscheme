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

from enum import Enum

from primer3 import calcTm

logger = logging.getLogger("primalscheme")


class Primer(object):
    """A primer."""

    class Direction(Enum):
        left = 1
        right = -1

    def __init__(self, seq, start, direction):
        self.direction = direction
        self.seq = seq
        self.start = start

    def __str__(self):
        return f"{self.direction}:{self.seq}:{self.start}"

    @property
    def size(self):
        return len(self.seq)

    @property
    def end(self):
        if self.direction == Primer.Direction.left:
            return self.start + self.size - 1
        elif self.direction == Primer.Direction.right:
            return self.start - self.size + 1


class CandidatePrimer(Primer):
    """A candidate primer."""

    def __init__(self, seq, start, direction, name="", penalty=0):
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


class CandidatePrimerPair(object):
    """A pair of candidate primers."""

    def __init__(self, left, right):
        self.left = left
        self.right = right

    def __str__(self):
        return f"{self.left} -:- {self.right}"

    @property
    def mean_identity(self):
        return (self.left.identity + self.right.identity) / 2

    @property
    def product_size(self):
        return self.right.start - self.left.start + 1

    @property
    def combined_penalty(self):
        return self.left.penalty + self.right.penalty


def get_alignment(primer, reference):
    """An seqan alignment of a primer against a reference."""

    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 2
    EXTEND = 1
    IDENTITY_THRESHOLD = 0.7

    if primer.direction == Primer.Direction.left:
        ref = reference.seq
    elif primer.direction == Primer.Direction.right:
        ref = reference.reverse_complement().seq

    # Semi-Global, do not penalize gaps at beginning and end of s2/database
    trace = parasail.sg_dx_trace_striped_sat(
        str(primer.seq), str(ref), OPEN, EXTEND, MATRIX
    )
    traceback = trace.get_traceback()

    query_end = trace.end_query + 1
    ref_end = trace.end_ref + 1

    # Get alignment strings
    aln_query = traceback.query[ref_end - query_end : ref_end]
    cigar = traceback.comp[ref_end - query_end : ref_end]
    aln_ref = traceback.ref[ref_end - query_end : ref_end]

    # Identity for glocal alignment
    identity = cigar.count("|") / len(cigar)

    # Alignment failed
    if identity < IDENTITY_THRESHOLD:
        return None

    # Format alignment
    refid = reference.id[:30]
    name = primer.name[:30]
    formatted_query = f"{name: <30} {1: >6} {aln_query} {query_end}"
    if primer.direction == Primer.Direction.left:
        formatted_ref = (
            f"{refid: <30} {ref_end - query_end + 1: >6} {aln_ref} {ref_end}"
        )
    elif primer.direction == Primer.Direction.right:
        rev_start = len(reference) - ref_end
        formatted_ref = (
            f"{refid: <30} {rev_start + query_end: >6} {aln_ref} {rev_start + 1}"
        )
    formatted_cigar = f"{'': <30} {'': >6} {cigar}"
    formatted_alignment = "\n".join(
        ["", formatted_query, formatted_cigar, formatted_ref]
    )

    del trace

    return (identity, formatted_alignment)


class reversor:
    """Decorator to reverse sort comparisons"""

    def __init__(self, obj):
        self.obj = obj

    def __eq__(self, other):
        return other.obj == self.obj

    def __lt__(self, other):
        return other.obj < self.obj

    def __gt__(self, other):
        return other.obj > self.obj
