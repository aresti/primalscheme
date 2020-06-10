"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains alignment functions.

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

import parasail

from collections import namedtuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Alignment = namedtuple("Alignment", "mismatches formatted_alignment")


class FailedAlignmentError(Exception):
    """No aligment between the primer and the refernce"""

    pass


def align_primer(primer, reference):
    """An seqan alignment of a primer against a reference."""

    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 2
    EXTEND = 1

    if primer.direction.name == "LEFT":
        ref = reference.seq
    elif primer.direction.name == "RIGHT":
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

    mismatches = cigar.count(".") + cigar.count(" ")

    # Format alignment
    refid = reference.id[:30]
    name = primer.direction.name
    formatted_query = f"{name: <30} {1: >6} {aln_query} {query_end}"
    if primer.direction.name == "LEFT":
        formatted_ref = (
            f"{refid: <30} {ref_end - query_end + 1: >6} {aln_ref} {ref_end}"
        )
    elif primer.direction.name == "RIGHT":
        rev_start = len(reference) - ref_end
        formatted_ref = (
            f"{refid: <30} {rev_start + query_end: >6} {aln_ref} {rev_start + 1}"
        )
    formatted_cigar = f"{'': <30} {'': >6} {cigar}"
    formatted_alignment = "\n".join(
        ["", formatted_query, formatted_cigar, formatted_ref]
    )

    del trace

    return Alignment(mismatches, formatted_alignment)


def align_secondary_reference(primary_flank, secondary_ref):
    """
    An seqan alignment of a primary reference slice against a complete
    secondary reference.
    """

    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 2
    EXTEND = 1

    # Semi-Global, do not penalize gaps at beginning and end of s2/database
    trace = parasail.sg_dx_trace_striped_sat(
        str(primary_flank.seq), str(secondary_ref.seq), OPEN, EXTEND, MATRIX
    )
    traceback = trace.get_traceback()
    query_end = trace.end_query
    ref_end = trace.end_ref
    aligned_query = traceback.query[ref_end - query_end : ref_end + 1]
    aligned_ref = traceback.ref[ref_end - query_end : ref_end + 1]

    del trace

    # Alignment failed (indels)
    if "-" in aligned_query + aligned_ref or len(primary_flank) != len(aligned_ref):
        raise FailedAlignmentError

    return SeqRecord(Seq(aligned_ref), id=secondary_ref.id)
