import pytest
import random
import uuid

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from primalscheme.cli import get_config


def seq_record_factory(seq_len=5000, alphabet="acgt", id=""):
    """Generate a random SeqRecord for testing purposes"""
    id = id or f"random_seq_{uuid.uuid4()}"
    seq = "".join([random.choice(alphabet) for i in range(seq_len)])
    return SeqRecord(Seq(seq), id=id)


def multi_seq_generator(n, **kwargs):
    """Generate multiple random SeqRecords for testing purposes"""
    return [seq_record_factory(**kwargs) for i in range(n)]


@pytest.fixture(scope="session")
def default_config():
    return get_config()


@pytest.fixture(scope="session")
def stored_inputs_path():
    """Return path to stored test inputs dir"""
    return Path(__file__).parent / "inputs"


@pytest.fixture(scope="session")
def temp_inputs_path(tmp_path_factory):
    """Return a temp_path for generating input files"""
    return tmp_path_factory.mktemp("inputs")


@pytest.fixture(scope="session")
def input_fasta_empty(temp_inputs_path):
    """Generate an empty fasta input file"""
    fh = temp_inputs_path / "empty_input.fa"
    fh.write_text("")
    return fh


@pytest.fixture(scope="session")
def input_fasta_invalid_alphabet(temp_inputs_path):
    """Generate a fasta file with one valid record and one with an invalid alphabet"""
    fh = temp_inputs_path / "invalid_alphabet.fa"
    records = [seq_record_factory(), seq_record_factory(alphabet="acgtRN")]
    SeqIO.write(records, fh, "fasta")
    return fh


@pytest.fixture(scope="session")
def input_fasta_5_random_valid(temp_inputs_path):
    """Generate a random multi-FASTA with 5 valid records"""
    fh = temp_inputs_path / "valid_random_5_fasta.fa"
    SeqIO.write(multi_seq_generator(5), fh, "fasta")
    return fh


@pytest.fixture(scope="session")
def input_fasta_101_random_valid(temp_inputs_path):
    """Generate a random multi-FASTA with 101 records"""
    fh = temp_inputs_path / "valid_random_101_fasta.fa"
    SeqIO.write(multi_seq_generator(101), fh, "fasta")
    return fh


@pytest.fixture(scope="session")
def input_fasta_valid_with_gaps(temp_inputs_path):
    """Generate a FASTA that includes gaps - """
    fh = temp_inputs_path / "valid_fasta_with_gaps.fa"
    SeqIO.write(seq_record_factory(alphabet="acgt-"), fh, "fasta")
    return fh


@pytest.fixture(scope="session")
def input_fasta_size_difference_over_500(temp_inputs_path):
    """Generate a FASTA where the  record is over 500nt longer than the first"""
    fh = temp_inputs_path / "invalid_different.fa"
    records = [seq_record_factory(seq_len=5000), seq_record_factory(seq_len=5501)]
    SeqIO.write(records, fh, "fasta")
    return fh


@pytest.fixture(scope="session")
def input_fasta_short_500(temp_inputs_path):
    """Generate a FASTA where the record is short"""
    fh = temp_inputs_path / "invalid_short.fa"
    SeqIO.write(seq_record_factory(seq_len=500), fh, "fasta")
    return fh


STORED_INPUTS = [
    "CHIKV_demo.fa",
    "Ebov-10-Pan.fasta",
]


@pytest.fixture(params=STORED_INPUTS)
def all_stored_inputs(request, stored_inputs_path):
    """Return all stored input paths"""
    return stored_inputs_path / request.param


@pytest.fixture(scope="session")
def chikv_input(stored_inputs_path):
    """Return all stored input paths"""
    return stored_inputs_path / STORED_INPUTS[0]
