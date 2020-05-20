import pytest
import warnings

from Bio.SeqRecord import SeqRecord

with warnings.catch_warnings():
    warnings.simplefilter("ignore", PendingDeprecationWarning)
    from primalscheme.cli import process_fasta


def test_process_fasta_empty_input(input_fasta_empty):
    with pytest.raises(ValueError, match="does not contain any valid references"):
        process_fasta(input_fasta_empty)


def test_process_fasta_invalid_alphabet(input_fasta_invalid_alphabet):
    with pytest.raises(ValueError, match="invalid nucleotide codes"):
        process_fasta(input_fasta_invalid_alphabet)


def test_process_fasta_too_many_records(input_fasta_101_random_valid):
    with pytest.raises(ValueError, match="A maximum of 100"):
        process_fasta(input_fasta_101_random_valid)


def test_process_fasta_size_difference_over_500(input_fasta_size_difference_over_500,):
    with pytest.raises(ValueError, match="too different in size"):
        process_fasta(input_fasta_size_difference_over_500)


def test_process_fasta_returns_list_of_seq_records(input_fasta_5_random_valid):
    references = process_fasta(input_fasta_5_random_valid)
    assert len(references) == 5
    assert all(isinstance(r, SeqRecord) for r in references)


def test_process_fasta_chikv_demo(chikv_input):
    references = process_fasta(chikv_input)
    assert len(references) == 2


def test_process_fasta_gaps_removed(input_fasta_valid_with_gaps):
    references = process_fasta(input_fasta_valid_with_gaps)
    assert "-" not in references[0]


def test_process_fasta_too_short_input(input_fasta_short_500):
    with pytest.raises(ValueError, match="too short"):
        process_fasta(input_fasta_short_500, min_ref_size=900)
