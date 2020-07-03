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

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3 import calcTm as p3_calcTm, calcHairpin as p3_calcHairpin
from primalscheme import config

Kmer = namedtuple("Kmer", "seq start")


class Direction(Enum):
    """Primer direction."""

    LEFT = "+"
    RIGHT = "-"


class Primer:
    """A primer."""

    def __init__(self, seq, start, direction, pool, reference_msa):
        """Init Primer."""
        self.seq = seq
        self.start = start
        self.direction = direction
        self.pool = pool
        self.reference_msa = reference_msa
        self.interacts_with = None

    def __str__(self):
        """Primer string representation."""
        return f"{self.direction.name}:{self.seq}:{self.start}"

    @property
    def size(self):
        """Primer size (length)."""
        return len(self.seq)

    @property
    def end(self):
        """Primer (inclusive) end position."""
        if self.direction == Direction.LEFT:
            return self.start + self.size - 1
        elif self.direction == Direction.RIGHT:
            return self.start - self.size + 1

    @property
    def seq(self):
        """Primer sequence."""
        return self.__seq

    @seq.setter
    def seq(self, seq):
        """Set seq and calculate derived characteristics."""
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
        if self.gc > config.PRIMER_GC_RANGE.opt:
            penalty += config.PRIMER_WT_GC_PERCENT_GT * (
                self.gc - config.PRIMER_GC_RANGE.opt
            )

        # Low GC
        if self.gc < config.PRIMER_GC_RANGE.opt:
            penalty += config.PRIMER_WT_GC_PERCENT_LT * (
                config.PRIMER_GC_RANGE.opt - self.gc
            )

        # High size
        if self.size > config.PRIMER_SIZE_RANGE.opt:
            penalty += config.PRIMER_WT_SIZE_GT * (
                self.size - config.PRIMER_SIZE_RANGE.opt
            )

        # Low size
        if self.size < config.PRIMER_SIZE_RANGE.opt:
            penalty += config.PRIMER_WT_SIZE_LT * (
                config.PRIMER_SIZE_RANGE.opt - self.size
            )

        self.__base_penalty = penalty
        return penalty

    def _mismatch_count_for_ref(self, ref):
        """Mismatch count for the primer against a given reference."""
        count = 0
        for i, base in enumerate(self.seq):
            if base != str(ref.seq)[i]:
                count += 1
        return count

    def _mismatch_penalties_for_ref(self, ref):
        """Mismatch penalties for each position against a given reference."""
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
        """List of mismatch counts between the primer and all references."""
        return [self._mismatch_count_for_ref(ref) for ref in self.reference_msa]

    @property
    def mismatch_penalty_matrix(self):
        """Matrix of mistmatch penalties for each position, against all references."""
        return [self._mismatch_penalties_for_ref(ref) for ref in self.reference_msa]

    @property
    def mismatch_penalties(self):
        """List of total mismatch penalties for the primer against all references."""
        return [sum(row) for row in self.mismatch_penalty_matrix]

    @property
    def combined_penalty(self):
        """Combined penalty (base + mistmatch)."""
        return self.base_penalty + sum(self.mismatch_penalties)

    @property
    def reference_msa(self):
        """Reference MSA, sliced to match primer position."""
        return self.__reference_msa

    @reference_msa.setter
    def reference_msa(self, msa):
        """Set reference MSA."""
        self.__reference_msa = msa
        self.__annotated_msa = None

    @property
    def annotated_msa(self):
        """Reference MSA with annotated cigar for logging (lazy eval)."""
        if self.__annotated_msa:
            return self.__annotated_msa

        # Calculate cigar
        cigar = ""
        for i in range(self.__reference_msa.get_alignment_length()):
            num = len(set(self.__reference_msa[::-1, i]))
            if num == 1:
                cigar += "*"
            else:
                cigar += " "
        cigar = SeqRecord(Seq(cigar), id="Cigar")
        self.__annotated_msa = self.__reference_msa[:]
        self.__annotated_msa.append(cigar)
        return self.__annotated_msa


def calc_gc(seq):
    """Calculate percent GC for a sequence."""
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def calc_tm(seq):
    """Calculate Tm for a sequence."""
    return p3_calcTm(
        seq,
        mv_conc=config.MV_CONC,
        dv_conc=config.DV_CONC,
        dntp_conc=config.DNTP_CONC,
        dna_conc=config.DNA_CONC,
    )


def calc_max_homo(seq):
    """Calculate max homopolymer length for a sequence."""
    return sorted([(len(list(g))) for k, g in groupby(seq)], reverse=True)[0]


def calc_hairpin(seq):
    """Calculate hairpin Tm for a sequence."""
    return p3_calcHairpin(
        seq,
        mv_conc=config.MV_CONC,
        dv_conc=config.DV_CONC,
        dntp_conc=config.DNTP_CONC,
        dna_conc=config.DNA_CONC,
    ).tm


def design_primers(msa, direction, pool, offset=0, primary_only=False):
    """Design primers against a reference MSA."""
    min_size = config.PRIMER_SIZE_RANGE.min
    max_size = config.PRIMER_SIZE_RANGE.max
    variation = max_size - min_size

    if primary_only:
        sequences = [msa[0].seq]
    else:
        sequences = [record.seq for record in msa]

    # Digest all sequences into k-mers
    all_kmers = set()
    for seq in sequences:
        for kmer_size in range(min_size, max_size + 1):
            all_kmers.update(digest_seq(str(seq[:-variation]), kmer_size))

    # Generate primers
    primers = []
    for kmer in all_kmers:
        sliced_msa = msa[:, kmer.start : kmer.start + len(kmer.seq)]
        if direction == Direction.LEFT:
            primer_start = offset + kmer.start
        else:
            primer_start = offset + len(msa[0].seq) - 1 - kmer.start
        primers.append(Primer(kmer.seq, primer_start, direction, pool, sliced_msa))

    return primers


def primer_thermo_filter(primer):
    """Hard filter for candidate primers."""
    return (
        (config.PRIMER_GC_RANGE.min <= primer.gc <= config.PRIMER_GC_RANGE.max)
        and (config.PRIMER_MIN_TM <= primer.tm <= config.PRIMER_MAX_TM)
        and (primer.hairpin <= config.PRIMER_MAX_HAIRPIN_TH)
        and (primer.max_homo <= config.PRIMER_MAX_HOMO)
    )


def digest_seq(seq, kmer_size):
    """Digest a sequence into kmers."""
    return [Kmer(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]
