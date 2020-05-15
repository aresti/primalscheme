
import pytest

from primalscheme.cli import get_arguments, process_fasta
from primalscheme.multiplex import MultiplexScheme


@pytest.fixture
def default_chikv_scheme(chikv_input):
    args = get_arguments(test=['multiplex', str(chikv_input)])
    return scheme_for_args(args)


def scheme_for_args(args):
    references = process_fasta(args.fasta)
    scheme = MultiplexScheme(
        references, args.amplicon_size, args.amplicon_max_variation,
        args.target_overlap, args.step_distance, args.min_unique, args.prefix,
        args.primer3)
    scheme.design_scheme()
    return scheme


def no_collisions(regions):
    pool_ids = (1, 2)
    for pool_id in pool_ids:
        pool = [r for r in regions if r.pool == pool_id]
        last_coord = -1
        for r in pool:
            if r.top_pair.left.start <= last_coord:
                return False
            last_coord = r.top_pair.right.start
    return True


def test_chikv_scheme_has_no_gaps(default_chikv_scheme):
    regions = default_chikv_scheme.regions
    inserts = [(r.top_pair.left.end + 1, r.top_pair.right.end)
               for r in regions]
    covered_coords = set()
    for insert in inserts:
        covered_coords.update(range(insert[0], insert[1]))
    all_coords = set(range(regions[0].top_pair.left.end + 1,
                           regions[-1].top_pair.right.end))

    assert all_coords - covered_coords == set()


def test_scheme_pools_do_not_have_collisions(all_stored_inputs):
    args = get_arguments(test=['multiplex', str(all_stored_inputs)])
    scheme = scheme_for_args(args)

    assert no_collisions(scheme.regions)


@pytest.mark.parametrize('amplicon_size', range(200, 900, 200))
def test_scheme_varying_amplicon_sizes(amplicon_size, chikv_input):
    max_variation = int(0.1 * amplicon_size / 2)
    args = get_arguments(test=[
        'multiplex', str(chikv_input),
        '--amplicon-size', str(amplicon_size),
        '--amplicon-max-variation', str(max_variation)])
    scheme = scheme_for_args(args)

    regions = scheme.regions
    for region in regions[:4]:
        product_size = region.top_pair.product_length
        assert product_size <= amplicon_size + max_variation
        assert product_size >= amplicon_size - max_variation


def test_large_target_overlap_does_not_result_in_collision(chikv_input):
    args = get_arguments(test=[
        'multiplex', str(chikv_input),
        '--target-overlap', str(200)])
    scheme = scheme_for_args(args)

    assert no_collisions(scheme.regions)
