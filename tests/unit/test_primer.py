from primalscheme.primer import calc_max_homo


def test_calc_max_homo():
    assert calc_max_homo("G") == 1
    assert calc_max_homo("GCAT") == 1
    assert calc_max_homo("AAAAA") == 5
    assert calc_max_homo("AAGGGTTAAAAAA") == 6
    assert calc_max_homo("GTCAGGGGAAAAC") == 4
