import pytest
import random
import uuid

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path


STORED_TEST_INPUTS = {
    'chikv': 'CHIKV_demo.50ca2db6b3ff.fa',
}

def seq_record_factory(seq_len=500, alphabet='acgt', id=''):
    """Generate a random SeqRecord for testing purposes"""
    id = id or f'random_seq_{uuid.uuid4()}'
    seq = ''.join([random.choice(alphabet) for i in range(seq_len)])
    return SeqRecord(Seq(seq), id=id)

def multi_seq_generator(n, **kwargs):
    """Generate multiple random SeqRecords for testing purposes"""
    return [seq_record_factory(**kwargs) for i in range(n)]

@pytest.fixture(scope='session')
def stored_inputs_path():
    """Return path to stored test inputs dir"""
    return Path(__file__).parent / 'inputs'

@pytest.fixture(scope='session')
def temp_inputs_path(tmp_path_factory):
    """Return a temp_path for generating input files"""
    return tmp_path_factory.mktemp('inputs')

@pytest.fixture(scope='session')
def input_fasta_empty(temp_inputs_path):
    """Generate an empty fasta input file"""
    fh = temp_inputs_path / "empty_input.fa"
    fh.write_text('')
    return fh

@pytest.fixture(scope='session')
def input_fasta_invalid_alphabet(temp_inputs_path):
    """Generate a fasta file with one valid record and one with an invalid alphabet"""
    fh = temp_inputs_path / "invalid_alphabet.fa"
    records = [seq_record_factory(), seq_record_factory(alphabet='acgtRN')]
    SeqIO.write(records, fh, 'fasta')
    return fh

@pytest.fixture(scope='session')
def input_fasta_5_random_valid(temp_inputs_path):
    """Generate a random multi-FASTA with 5 valid records"""
    fh = temp_inputs_path / "valid_random_5_fasta.fa"
    SeqIO.write(multi_seq_generator(5), fh, 'fasta')
    return fh

@pytest.fixture(scope='session')
def input_fasta_101_random_valid(temp_inputs_path):
    """Generate a random multi-FASTA with 101 records"""
    fh = temp_inputs_path / "valid_random_101_fasta.fa"
    SeqIO.write(multi_seq_generator(101, seq_len=10), fh, 'fasta')
    return fh

@pytest.fixture(scope='session')
def input_fasta_valid_with_gaps(temp_inputs_path):
    """Generate a FASTA that includes gaps - """
    fh = temp_inputs_path / "valid_random_1_fasta_with_gaps.fa"
    SeqIO.write(seq_record_factory(seq_len=500, alphabet='acgt-'), fh, 'fasta')
    return fh

@pytest.fixture(scope='session')
def input_fasta_length_difference_over_500(temp_inputs_path):
    """Generate a multi-FASTA where the second record is over 500nt longer than the first"""
    fh = temp_inputs_path / "valid_random_1_fasta_with_gaps.fa"
    records = [seq_record_factory(seq_len=10), seq_record_factory(seq_len=511)]
    SeqIO.write(records, fh, 'fasta')
    return fh

@pytest.fixture(scope='session')
def input_fasta_chikv(stored_inputs_path):
    """Return CHIKV test fasta"""
    return stored_inputs_path / STORED_TEST_INPUTS['chikv']