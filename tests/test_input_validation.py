import os
import pytest

from Bio.SeqRecord import SeqRecord
from primalscheme.primalscheme import process_fasta


def input_file_path(file_name):
   return os.path.join(os.path.dirname(__file__), 'inputs/{}'.format(file_name))

def test_empty_input_file():
    with pytest.raises(ValueError):
        process_fasta(input_file_path('empty_file.fa'))

def test_invalid_alphabet_input_file():
    with pytest.raises(ValueError):
        process_fasta(input_file_path('10_acgtrn_seq.fa'))

def test_valid_fasta_returns_list_of_seq_records():
    references = process_fasta(input_file_path('10_acgt_seq.fa'))
    assert len(references) == 10
    assert all(isinstance(r, SeqRecord) for r in references)

