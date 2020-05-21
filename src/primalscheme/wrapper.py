import logging
import primer3

from collections import namedtuple


logger = logging.getLogger("primalscheme")


class InsufficientPrimersError(Exception):
    """Unable to find sufficient unique primers at the current cursor position."""

    pass


def design_primers(seq, p3_global, min_unique, offset=0):
    """Find primer pairs for a sequence slice."""

    SimplePrimer = namedtuple("SimplePrimer", "seq start penalty")

    p3_seq = {
        "SEQUENCE_TEMPLATE": seq,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [-1, -1, -1, -1],
        "SEQUENCE_INCLUDED_REGION": [-1, -1],
    }

    # Call primer3
    p3_output = primer3.bindings.designPrimers(p3_seq, p3_global)

    # Parse result
    pairs = ([], [])
    text_dir = ("LEFT", "RIGHT")

    for d in range(2):
        num_returned = p3_output[f"PRIMER_{text_dir[d]}_NUM_RETURNED"]
        for i in range(num_returned):
            seq = str(p3_output[f"PRIMER_{text_dir[d]}_{i}_SEQUENCE"])
            penalty = calc_penalty(seq, p3_global)
            start = offset + int(p3_output[f"PRIMER_{text_dir[d]}_{i}"][0])
            pairs[d].append(SimplePrimer(seq, start, penalty))

    # If we don't have min unique left and right, then raise
    unique_count = min(len(set(pairs[0])), len(set(pairs[1])))
    logger.debug(f"Primer3 returned {unique_count} unique pairs")

    if unique_count < min_unique:
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
