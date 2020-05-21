import primer3

from primalscheme.wrapper import calc_penalty


def test_calc_penalty_matches_primer3_penalty(random_reference_slice, default_config):
    """Check internal penalty calculation matches primer3 output"""

    p3_global = default_config["primer3"]
    p3_seq = {
        "SEQUENCE_TEMPLATE": random_reference_slice,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [-1, -1, -1, -1],
        "SEQUENCE_INCLUDED_REGION": [-1, -1],
    }

    p3_output = primer3.bindings.designPrimers(p3_seq, p3_global)
    text_dir = ("LEFT", "RIGHT")

    for d in range(2):
        num_returned = p3_output[f"PRIMER_{text_dir[d]}_NUM_RETURNED"]
        for i in range(num_returned):
            seq = str(p3_output[f"PRIMER_{text_dir[d]}_{i}_SEQUENCE"])
            p3_penalty = float(p3_output[f"PRIMER_{text_dir[d]}_{i}_PENALTY"])
            penalty = calc_penalty(seq, p3_global)

            assert penalty == p3_penalty
