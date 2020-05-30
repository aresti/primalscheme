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
from collections import namedtuple
from itertools import groupby
from enum import Enum
from primer3 import calcTm, calcHairpin

logger = logging.getLogger("primalscheme")

Kmer = namedtuple("Kmer", "seq start")
Alignment = namedtuple("Alignment", "mismatches formatted_alignment")
SimplePrimer = namedtuple("SimplePrimer", "seq start direction gc tm hairpin max_homo")


class FailedAlignmentError(Exception):
    """No aligment between the primer and the refernce"""

    pass


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

    MISMATCH_PENALTY = 1

    def __init__(self, simple_primer, base_penalty, name=""):
        super().__init__(
            simple_primer.seq, simple_primer.start, simple_primer.direction
        )
        self.name = name
        self.tm = simple_primer.tm
        self.gc = simple_primer.gc
        self.alignments = []
        self.alignment_cigars = []
        self.base_penalty = base_penalty
        self.identity = 0

    def align(self, references):
        for ref in references[1:]:
            # align against non-primary references
            alignment = get_alignment(self, ref)
            self.alignments.append(alignment)

        # Calculate average
        if self.alignments:
            scores = [
                (len(self.seq) - a.mismatches) / len(self.seq) for a in self.alignments
            ]
            self.identity = sum(scores) / len(scores)
            return

    @property
    def mismatch_counts(self):
        return [cigar.count(".") for cigar in self.alignment_cigars]

    @property
    def mismatch_penalty_matrix(self):
        PENALTIES = [3, 2, 1]
        matrix = []
        for cigar in self.alignment_cigars:
            row = []
            for i, base in enumerate(cigar):
                if base == "|":
                    row.append(0)
                    continue
                dist_3p = self.size - i
                if dist_3p > len(PENALTIES) - 1:
                    row.append(PENALTIES[-1])
                else:
                    row.append(PENALTIES[dist_3p])
            matrix.append(row)
        return matrix

    @property
    def mismatch_penalties(self):
        return [sum(row) for row in self.mismatch_penalty_matrix]

    @property
    def total_penalty(self):
        return self.base_penalty + sum(self.mismatch_penalties)


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
        return self.left.total_penalty + self.right.total_penalty


def get_alignment(primer, reference):
    """An seqan alignment of a primer against a reference."""

    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 2
    EXTEND = 1
    MISMATCH_THRESHOLD = 10

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

    mismatches = cigar.count(".")

    # Alignment failed
    if mismatches > MISMATCH_THRESHOLD:
        raise FailedAlignmentError

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

    return Alignment(mismatches, formatted_alignment)


def align_secondary_reference(primary_ref_slice, ref, limit=None):
    """
    An seqan alignment of a primary reference slice against a complete
    secondary reference.
    """

    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 2
    EXTEND = 1

    # Semi-Global, do not penalize gaps at beginning and end of s2/database
    if limit and limit > 0:
        primary_ref_slice = primary_ref_slice[: limit + 1]
    elif limit and limit < 0:
        primary_ref_slice = primary_ref_slice[len(primary_ref_slice) + limit :]
    trace = parasail.sg_dx_trace_striped_sat(
        str(primary_ref_slice), str(ref.seq), OPEN, EXTEND, MATRIX
    )
    traceback = trace.get_traceback()
    query_end = trace.end_query
    ref_end = trace.end_ref
    cigar = traceback.comp[ref_end - query_end : ref_end]

    # Alignment failed (indels)
    if " " in cigar:
        raise FailedAlignmentError

    del trace

    return cigar


class InsufficientPrimersError(Exception):
    """Unable to find sufficient unique primers at the current cursor position."""

    pass


def design_primers(seq, offset, p3_global, reverse=False):
    min_size = p3_global["PRIMER_MIN_SIZE"]
    max_size = p3_global["PRIMER_MAX_SIZE"]
    variation = max_size - min_size

    # Digest seq into k-mers
    all_kmers = set()
    for kmer_size in range(min_size, max_size + 1):
        all_kmers.update(digest_seq(str(seq[:-variation]), kmer_size))

    # Generate SimplePrimers
    simple_primers = [
        SimplePrimer(
            kmer.seq,
            offset + len(seq) - 1 - kmer.start if reverse else offset + kmer.start,
            Primer.Direction.right if reverse else Primer.Direction.left,
            calc_gc(kmer.seq),
            calc_tm(kmer.seq),
            calc_hairpin(kmer.seq),
            calc_max_homo(kmer.seq),
        )
        for kmer in all_kmers
    ]

    # Hard filter
    filtered_candidates = [
        CandidatePrimer(primer, calc_base_penalty(primer.seq, p3_global))
        for primer in simple_primers
        if hard_filter(primer)
    ]
    return filtered_candidates


def hard_filter(primer):
    return (
        (30 <= primer.gc <= 55)
        and (60 <= primer.tm <= 63)
        and (primer.hairpin <= 50.0)
        and (primer.max_homo <= 5)
    )


def calc_base_penalty(seq, p3_global):
    """
    Calculate primer penalty score.
    As per proutine described in http://primer3.ut.ee/primer3web_help.htm
    """

    penalty = 0
    tm = calc_tm(seq)
    gc = calc_gc(seq)
    length = calc_length(seq)

    # High Tm
    if tm > p3_global["PRIMER_OPT_TM"]:
        penalty += p3_global["PRIMER_WT_TM_GT"] * (tm - p3_global["PRIMER_OPT_TM"])

    # Low Tm
    if tm < p3_global["PRIMER_OPT_TM"]:
        penalty += p3_global["PRIMER_WT_TM_LT"] * (p3_global["PRIMER_OPT_TM"] - tm)

    # High GC
    if gc > p3_global["PRIMER_OPT_GC_PERCENT"]:
        penalty += p3_global["PRIMER_WT_GC_PERCENT_GT"] * (
            gc - p3_global["PRIMER_OPT_GC_PERCENT"]
        )

    # Low GC
    if gc < p3_global["PRIMER_OPT_GC_PERCENT"]:
        penalty += p3_global["PRIMER_WT_GC_PERCENT_LT"] * (
            p3_global["PRIMER_OPT_GC_PERCENT"] - gc
        )

    # High size
    if length > p3_global["PRIMER_OPT_SIZE"]:
        penalty += p3_global["PRIMER_WT_SIZE_GT"] * (
            length - p3_global["PRIMER_OPT_SIZE"]
        )
    # Low size
    if length < p3_global["PRIMER_OPT_SIZE"]:
        penalty += p3_global["PRIMER_WT_SIZE_LT"] * (
            p3_global["PRIMER_OPT_SIZE"] - length
        )

    return penalty


# Calc Tm
def calc_tm(seq):
    return calcTm(seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)


# Calc GC content
def calc_gc(seq):
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


# Calc length of primer
def calc_length(seq):
    return len(seq)


# Calc max homopolymer length using itertools
def calc_max_homo(seq):
    return sorted([(len(list(g))) for k, g in groupby(seq)], reverse=True)[0]


# Calc hairpin stability
def calc_hairpin(seq):
    return calcHairpin(seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm


# Digest seq into k-mers
def digest_seq(seq, kmer_size):
    return [Kmer(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]
