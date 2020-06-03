"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains classes and functions related to primer design.

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

from collections import namedtuple
from enum import Enum
from itertools import groupby

from primer3 import calcTm as p3_calcTm, calcHairpin as p3_calcHairpin
from primalscheme import config
from primalscheme.align import align_primer

Kmer = namedtuple("Kmer", "seq start")


class Direction(Enum):
    """Primer direction."""

    left = "left"
    right = "right"


class Primer:
    """A primer."""

    def __init__(self, seq, start, direction):
        self.direction = direction
        self.seq = seq
        self.start = start
        self.alignments = []
        self.reference_msa = None

    def __str__(self):
        return f"{self.direction.value}:{self.seq}:{self.start}"

    @property
    def size(self):
        return len(self.seq)

    @property
    def end(self):
        """The (inclusive) end position of the primer"""
        if self.direction == Direction.left:
            return self.start + self.size - 1
        elif self.direction == Direction.right:
            return self.start - self.size + 1

    @property
    def seq(self):
        return self.__seq

    @seq.setter
    def seq(self, seq):
        """Set seq and calculate derived characteristics"""
        self.__seq = seq
        self.gc = calc_gc(self.seq)
        self.tm = calc_tm(self.seq)
        self.hairpin = calc_hairpin(self.seq)
        self.max_homo = calc_max_homo(self.seq)
        self.__base_penalty = None

    @property
    def base_penalty(self):
        """
        Calculate intrinsic primer penalty.
        As per routine described in http://primer3.ut.ee/primer3web_help.htm
        Lazy-evaluation.
        """
        if self.__base_penalty is not None:
            return self.__base_penalty

        penalty = 0

        # High Tm
        if self.tm > config.PRIMER_OPT_TM:
            penalty += config.PRIMER_WT_TM_GT * (self.tm - config.PRIMER_OPT_TM)

        # Low Tm
        if self.tm < config.PRIMER_OPT_TM:
            penalty += config.PRIMER_WT_TM_LT * (config.PRIMER_OPT_TM - self.tm)

        # High GC
        if self.gc > config.PRIMER_OPT_GC_PERCENT:
            penalty += config.PRIMER_WT_GC_PERCENT_GT * (
                self.gc - config.PRIMER_OPT_GC_PERCENT
            )

        # Low GC
        if self.gc < config.PRIMER_OPT_GC_PERCENT:
            penalty += config.PRIMER_WT_GC_PERCENT_LT * (
                config.PRIMER_OPT_GC_PERCENT - self.gc
            )

        # High size
        if self.size > config.PRIMER_OPT_SIZE:
            penalty += config.PRIMER_WT_SIZE_GT * (self.size - config.PRIMER_OPT_SIZE)

        # Low size
        if self.size < config.PRIMER_OPT_SIZE:
            penalty += config.PRIMER_WT_SIZE_LT * (config.PRIMER_OPT_SIZE - self.size)

        self.__base_penalty = penalty
        return penalty

    def _align(self, references):  # TODO, derive from msa
        """Align primer against secondary references for debug purposes"""
        for ref in references[1:]:
            alignment = align_primer(self, ref)
            self.alignments.append(alignment)

    def _mismatch_count_for_ref(self, ref):
        count = 0
        for i, base in enumerate(self.seq):
            if base != str(ref.seq)[i]:
                count += 1
        return count

    def _mismatch_penalties_for_ref(self, ref):
        positions = []
        penalties = config.PRIMER_MISMATCH_PENALTIES
        ref_seq = str(ref.seq)
        for i, base in enumerate(self.seq):
            if base == ref_seq[i]:
                positions.append(0)
                continue
            dist_3p = self.size - i
            if dist_3p > len(penalties) - 1:
                positions.append(penalties[-1])
            else:
                positions.append(penalties[dist_3p])
        return positions

    @property
    def mismatch_counts(self):
        return [self._mismatch_count_for_ref(ref) for ref in self.reference_msa]

    @property
    def mismatch_penalty_matrix(self):
        return [self._mismatch_penalties_for_ref(ref) for ref in self.reference_msa]

    @property
    def mismatch_penalties(self):
        return [sum(row) for row in self.mismatch_penalty_matrix]

    @property
    def combined_penalty(self):
        return self.base_penalty + sum(self.mismatch_penalties)


def calc_gc(seq):
    """Calculate percent GC for a sequence"""
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def calc_tm(seq):
    """Calculate Tm for a sequence"""
    return p3_calcTm(
        seq, mv_conc=config.MV_CONC, dv_conc=config.DV_CONC, dntp_conc=config.DNTP_CONC
    )


def calc_max_homo(seq):
    """Calculate max homopolymer length for a sequence"""
    return sorted([(len(list(g))) for k, g in groupby(seq)], reverse=True)[0]


def calc_hairpin(seq):
    """Calculate hairpin Tm for a sequence"""
    return p3_calcHairpin(
        seq, mv_conc=config.MV_CONC, dv_conc=config.DV_CONC, dntp_conc=config.DNTP_CONC
    ).tm


def design_primers(msa, offset=0, reverse=False):
    min_size = config.PRIMER_SIZE_MIN
    max_size = config.PRIMER_SIZE_MAX
    variation = max_size - min_size

    sequences = [record.seq for record in msa]

    # Digest all sequences into k-mers
    all_kmers = set()
    for seq in sequences:
        for kmer_size in range(min_size, max_size + 1):
            all_kmers.update(digest_seq(str(seq[:-variation]), kmer_size))

    # Generate primers
    primers = [
        Primer(
            kmer.seq,
            offset + len(msa[0].seq) - 1 - kmer.start
            if reverse
            else offset + kmer.start,
            Direction.right if reverse else Direction.left,
        )
        for kmer in all_kmers
    ]

    # Hard filter
    filtered_candidates = [primer for primer in primers if hard_filter(primer)]
    return filtered_candidates


def hard_filter(primer):
    return (
        (config.PRIMER_MIN_GC <= primer.gc <= config.PRIMER_MAX_GC)
        and (config.PRIMER_MIN_TM <= primer.tm <= config.PRIMER_MAX_TM)
        and (primer.hairpin <= config.PRIMER_MAX_HAIRPIN_TH)
        and (primer.max_homo <= config.PRIMER_MAX_HOMO)
    )


def digest_seq(seq, kmer_size):
    return [Kmer(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]
