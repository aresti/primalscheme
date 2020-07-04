import pytest

from Bio import SeqIO

from primalscheme.cli import process_fasta
from primalscheme.multiplex import MultiplexScheme
from primalscheme.primer import Direction


@pytest.fixture(scope="session")
def chikv_scheme(chikv_input):
    return get_scheme(chikv_input)


@pytest.fixture(scope="session")
def ebola_scheme(ebola_input):
    return get_scheme(ebola_input)


def get_scheme(fasta, **kwargs):
    references = process_fasta(fasta)
    scheme = MultiplexScheme(references, **kwargs)
    scheme.design_scheme()
    return scheme


def no_collisions(regions):
    pool_ids = (1, 2)
    for pool_id in pool_ids:
        pool = [r for r in regions if r.pool == pool_id]
        last_coord = -1
        for r in pool:
            if r.left.start <= last_coord:
                return False
            last_coord = r.right.start
    return True


def test_chikv_scheme_has_no_gaps(chikv_scheme):
    regions = chikv_scheme.regions
    inserts = [(r.left.end + 1, r.right.end) for r in regions]
    covered_coords = set([x for insert in inserts for x in range(*insert)])
    all_coords = set(range(regions[0].left.end + 1, regions[-1].right.end))

    assert all_coords.issubset(covered_coords)


def test_no_collisions_in_any_test_scheme(all_stored_inputs):
    scheme = get_scheme(all_stored_inputs)
    assert no_collisions(scheme.regions)


@pytest.mark.parametrize("amplicon_size", range(200, 801, 200))
def test_scheme_varying_amplicon_sizes(amplicon_size, chikv_input):
    variation = int(0.1 * amplicon_size / 2)
    amplicon_size_min = amplicon_size - variation
    amplicon_size_max = amplicon_size + variation
    scheme = get_scheme(
        chikv_input,
        amplicon_size_min=amplicon_size_min,
        amplicon_size_max=amplicon_size_max,
    )

    regions = scheme.regions
    for region in regions[:4]:
        product_size = region.product_size
        assert product_size <= amplicon_size_max
        assert product_size >= amplicon_size_min


def test_large_target_overlap_does_not_result_in_collision(chikv_input):
    scheme = get_scheme(chikv_input, target_overlap=300)

    assert no_collisions(scheme.regions)


def test_candidates_are_sorted_first_by_penalty(chikv_scheme):
    regions = chikv_scheme.regions
    for region in regions:
        candidates = region.left_candidates + region.right_candidates
        prev = None
        direction = "LEFT"
        for candidate in candidates:
            if candidate.direction.name != direction:
                direction = "RIGHT"
                prev = None
            if prev:
                if candidate.combined_penalty < prev:
                    pytest.fail("Candidates are not sorted first by combined_penalty.")
            prev = candidate.combined_penalty


def test_candidates_are_sorted_second_by_start(chikv_scheme):
    regions = chikv_scheme.regions
    for region in regions:
        candidates = region.left_candidates + region.right_candidates
        prev = None
        direction = "LEFT"
        for candidate in candidates:
            attrs = (candidate.combined_penalty, candidate.start)
            if candidate.direction.name != direction:
                direction = "RIGHT"
                prev = None
            if prev:
                if attrs[0] == prev[0] and attrs[1] > prev[1]:
                    pytest.fail("Candidates are not sorted second by start position.")
            prev = attrs


def test_candidates_are_sorted_third_by_seq(ebola_scheme):
    regions = ebola_scheme.regions
    for region in regions:
        candidates = region.left_candidates + region.right_candidates
        prev = None
        direction = "LEFT"
        for candidate in candidates:
            attrs = (candidate.combined_penalty, candidate.start, candidate.seq)
            if candidate.direction.name != direction:
                direction = "RIGHT"
                prev = None
            if prev:
                if attrs[0] == prev[0] and attrs[1] == prev[1] and attrs[2] < prev[2]:
                    pytest.fail("Candidates are not sorted third by seq.")
            prev = attrs


def test_left_primer_seq_matches_some_ref_slice(chikv_scheme):
    left = chikv_scheme.regions[0].left
    ref_slices = left.reference_msa

    assert any(ref_slice.seq == left.seq for ref_slice in ref_slices)


def test_right_primer_seq_matches_some_ref_slice(chikv_scheme):
    right = chikv_scheme.regions[0].right
    ref_slices = right.reference_msa

    assert any(ref_slice.seq == right.seq for ref_slice in ref_slices)


def test_multiple_references_used_for_primer_design(chikv_scheme):
    for dir in Direction:
        assert any(
            map(
                lambda x: x.seq != x.reference_msa[0].seq,
                [p for p in chikv_scheme.primers if p.direction == dir],
            )
        )


def test_scheme_runs_with_single_reference(ncov_single_ref_input):
    references = process_fasta(ncov_single_ref_input)
    scheme = MultiplexScheme(references)
    scheme.design_scheme()

    assert len(scheme.regions) > 30


def test_first_reference_is_primary(ebola_input, ebola_scheme):
    records = list(SeqIO.parse(ebola_input, "fasta"))
    assert records[0].id == ebola_scheme.primary_ref.id


def test_first_only_option_has_no_mistmatches_against_primary(chikv_input):
    scheme = get_scheme(chikv_input, primary_only=True)
    for dir in Direction:
        assert all(
            map(
                lambda x: x.seq == x.reference_msa[0].seq,
                [p for p in scheme.primers if p.direction == dir],
            )
        )
