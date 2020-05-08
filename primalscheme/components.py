"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains classes that constitute the various components of a
multiplex primer scheme.

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

from Bio.Seq import Seq
from primer3 import calcTm, calcHairpin, calcHomodimer
from porechop.cpp_function_wrappers import adapter_alignment

logger = logging.getLogger("primalscheme")


class Primer(object):
    """A simple primer."""

    def __init__(self, seq, start, direction):
        self.direction = direction
        self.seq = seq
        self.start = start

    @property
    def length(self):
        return len(self.seq)


class CandidatePrimer(Primer):
    """A candidate primer."""

    def __init__(self, seq, start, direction, name="", penalty=None):
        super().__init__(seq, start, direction)  # TODO tidy up
        self.name = name
        self.penalty = penalty

        self.percent_identity = 0
        self.tm = calcTm(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)
        self.homodimer = calcHomodimer(
            self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6
        ).tm
        self.hairpin = calcHairpin(self.seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm
        self.gc = 100.0 * (seq.count("G") + seq.count("C")) / len(seq)
        self.alignments = []

    def align(self, references):
        for ref in references:
            alignment = CAlignment(self, ref)
            self.alignments.append(alignment)
        # Calculate average percent identity
        idents = [i.percent_identity for i in self.alignments if i.percent_identity]
        if idents:
            self.percent_identity = sum(idents) / len(idents)
        return self

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
    def mean_percent_identity(self):
        return (self.left.percent_identity + self.right.percent_identity) / 2

    @property
    def product_length(self):
        return self.right.start - self.left.start + 1


class CAlignment(object):
    """An seqan alignment of a primer against a reference."""

    MISMATCHES = [
        set(["A", "A"]),
        set(["A", "C"]),
        set(["C", "C"]),
        set(["G", "A"]),
        set(["G", "G"]),
    ]

    def __init__(self, primer, ref):
        self.start = None
        self.end = None
        self.length = None
        self.percent_identity = None
        self.aln_query = None
        self.aln_ref = None
        self.aln_ref_comp = None
        self.ref_id = None
        self.mm_3prime = None
        self.cigar = None
        self.formatted_alignment = None

        if primer.direction == "LEFT":
            alignment_result = adapter_alignment(
                str(ref.seq), str(primer.seq), [2, -1, -2, -1]
            )
        elif primer.direction == "RIGHT":
            alignment_result = adapter_alignment(
                str(ref.seq.reverse_complement()), str(primer.seq), [2, -1, -2, -1]
            )
        result_parts = alignment_result.split(",")
        ref_start = int(result_parts[0])
        full_primer_percent_identity = float(result_parts[6])

        # If the read start is -1, that indicates that the alignment
        # failed completely.
        if ref_start == -1 or full_primer_percent_identity < 70:
            return
        else:
            ref_end = int(result_parts[1]) + 1

            if primer.direction == "LEFT":
                self.start = ref_start
                self.end = ref_end
                self.length = self.end - self.start
            else:
                self.start = len(ref) - ref_start
                self.end = len(ref) - (int(result_parts[1]) + 1)
                self.length = self.start - self.end

            # Percentage identity for glocal alignment
            self.percent_identity = full_primer_percent_identity

            # Get alignment strings
            self.aln_query = result_parts[8][ref_start:ref_end]
            self.aln_ref = result_parts[7][ref_start:ref_end]
            self.aln_ref_comp = Seq(str(self.aln_ref)).complement()
            self.ref_id = ref.id
            self.mm_3prime = False

            # Make cigar
            self.cigar = ""
            for a, b in zip(self.aln_query, self.aln_ref):
                if a == "-" or b == "-":
                    self.cigar += " "
                    continue
                if a != b:
                    self.cigar += "*"
                    continue
                else:
                    self.cigar += "|"

            # Format alignment
            short_primer = primer.name[:30] if len(primer.name) > 30 else primer.name
            short_ref = ref.id[:30] if len(ref.id) > 30 else ref.id
            self.formatted_alignment = "\n{: <30}5'-{}-3'\n{: <33}{}\n{: <30}3'-{}-5'".format(
                short_primer,
                self.aln_query,
                "",
                self.cigar,
                short_ref,
                self.aln_ref_comp,
            )

            # Check 3' mismatches
            if (
                set([self.aln_query[-1], self.aln_ref_comp[-1]])
                in CAlignment.MISMATCHES
            ):
                self.mm_3prime = True
                self.percent_identity = 0
