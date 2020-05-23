import logging
import primer3
from Bio import Seq

from collections import namedtuple
from itertools import groupby


logger = logging.getLogger("primalscheme")


class InsufficientPrimersError(Exception):
    """Unable to find sufficient unique primers at the current cursor position."""

    pass


def generate(seq, offset, p3_global, direction):

    SimplePrimer = namedtuple(
        "SimplePrimer", "seq start gc tm hairpin max_homo penalty"
    )

    # Digest seq into k-mers
    all_kmers = set()
    for kmer_size in range(
        p3_global["PRIMER_MIN_SIZE"], p3_global["PRIMER_MAX_SIZE"] + 1
    ):
        all_kmers.update(digest_seq(seq, kmer_size))

    # Filter for valid start position only
    valid_positions = [k for k in all_kmers if 0 <= k[1] < 40]

    # Generate SimplePrimers
    simple_primers = [
        SimplePrimer(
            k[0],
            offset + k[1] if direction == "left" else offset + len(seq) - k[1],
            calc_gc(k[0]),
            calc_tm(k[0]),
            calc_hairpin(k[0]),
            calc_max_homo(k[0]),
            calc_penalty(k[0], p3_global),
        )
        for k in valid_positions
    ]

    # Filter function
    def hard_filter(p):
        return (
            (30 <= p.gc <= 55)
            and (60 <= p.tm <= 63)
            and (p.hairpin <= 50.0)
            and (p.max_homo <= 5)
        )

    # Perform the hard filtering
    filtered_primers = [p for p in simple_primers if hard_filter(p)]
    return filtered_primers


def design_primers(seq, p3_global, min_unique, offset=0):
    """Find primer pairs for a sequence slice."""

    left_primers = generate(seq, offset, p3_global, "left")
    right_primers = generate(
        str(Seq.Seq(seq).reverse_complement()), offset, p3_global, "right"
    )

    filt_pairs = [(f, r) for f in left_primers for r in right_primers]
    pairs = ([], [])
    for d in range(2):
        for i in sorted(filt_pairs, key=lambda x: (x[0].penalty + x[1].penalty))[:10]:
            pairs[d].append(i[d])

    if len(filt_pairs) < min_unique:
        logger.debug(f"Does not satisfy {min_unique} min unique")
        raise InsufficientPrimersError(
            f"Failed to find {min_unique} unique left or right primers."
        )

    return pairs


def calc_penalty(seq, p3_global):
    """
    Calculate primer penalty score.
    As per proutine described in http://primer3.ut.ee/primer3web_help.htm
    """

    # High Tm
    penalty = 0
    tm = calc_tm(seq)
    gc = calc_gc(seq)
    length = calc_length(seq)

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
    return primer3.calcTm(seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6)


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
    return primer3.calcHairpin(seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.6).tm


# Digest seq into k-mers
def digest_seq(seq, kmer_size):
    return [(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]
