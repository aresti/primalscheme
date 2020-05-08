import primer3

from collections import namedtuple


class InsufficientPrimersError(Exception):
    """Unable to find sufficient unique primers at the current cursor position."""
    pass


def design_primers(seq, p3_global, min_unique, offset=0):
    """Find primer pairs for a sequence slice."""

    SimplePrimer = namedtuple('SimplePrimer', 'seq start penalty')

    p3_seq = {
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [-1, -1, -1, -1],
        'SEQUENCE_INCLUDED_REGION': [-1, -1],
    }

    # Call primer3
    p3_output = primer3.bindings.designPrimers(p3_seq, p3_global)

    # Parse result
    pairs = ([], [])
    text_dir = ('LEFT', 'RIGHT')

    for d in range(2):
        num_returned = p3_output[f'PRIMER_{text_dir[d]}_NUM_RETURNED']
        for i in range(num_returned):
            seq = str(p3_output[f'PRIMER_{text_dir[d]}_{i}_SEQUENCE'])
            penalty = float(p3_output[f'PRIMER_{text_dir[d]}_{i}_PENALTY'])
            start = offset + int(p3_output[f'PRIMER_{text_dir[d]}_{i}'][0])
            if d == 1:
                """
                This is a deliberate re-introduction of bug,
                to compare refactored output to old working version
                """
                start += 1
            pairs[d].append(SimplePrimer(seq, start, penalty))
    
    # If we don't have min unique left and right, then raise
    if any(len(set(side)) < min_unique for side in pairs):
        raise InsufficientPrimersError(
            f'Failed to find {min_unique} unique left or right primers.')

    return pairs
