"""
Tests for command line interface (CLI)
"""
import os
import pytest
import re
import warnings
import primalscheme.cli

from argparse import Namespace
from Bio.SeqRecord import SeqRecord

with warnings.catch_warnings():
    warnings.simplefilter("ignore", PendingDeprecationWarning)
    from primalscheme.cli import process_fasta


def test_runas_module():
    """
    Can this package be run as a Python module?
    """
    exit_status = os.system("python -m primalscheme --help")
    assert exit_status == 0


def test_entrypoint():
    """
    Is entrypoint script installed? (setup.py)
    """
    exit_status = os.system("primalscheme --help")
    assert exit_status == 0


def test_mutliplex_command():
    """
    Is multiplex command available?
    """
    exit_status = os.system("primalscheme multiplex --help")
    assert exit_status == 0


def test_cli_fails_without_command():
    """
    Does CLI stop execution w/o a command argument?
    """
    with pytest.raises(SystemExit):
        primalscheme.cli.main()
        pytest.fail("CLI doesn't abort asking for a command argument")


def test_cli_fails_without_fasta():
    """
    Does CLI stop execution w/o a positional fasta arg?
    """
    with pytest.raises(SystemExit):
        args = ["multiplex"]
        primalscheme.cli.parse_arguments(args)


@pytest.mark.parametrize(
    "option, value",
    [
        ("--prefix", "myScheme"),
        ("--amplicon-size-min", "380"),
        ("--amplicon-size-max", "420"),
        ("--target-overlap", "0"),
        ("--output-path", "./scheme"),
        ("--step-distance", "11"),
        ("--debug", None),
        ("--no-sort", None),
        ("--force", None),
    ],
)
def test_availability_of_cli_options(option, value):
    """
    Are options available?
    """
    args = ["multiplex", "some.fa", option]
    if value:
        args.append(value)
    parsed = primalscheme.cli.parse_arguments(args)
    assert isinstance(parsed, Namespace)


def test_cli_fails_with_invalid_option(chikv_input):
    """
    Does CLI stop execution with an invalid option?
    """
    with pytest.raises(SystemExit):
        args = ["multiplex", f"{str(chikv_input)}", "--force", "--badOption"]
        primalscheme.cli.parse_arguments(args)


@pytest.mark.parametrize(
    "option", ["-V", "--version"],
)
def test_cli_version_output(option, capsys):
    """
    Does CLI output a sensible version number for -V and --version?
    """
    try:
        primalscheme.cli.parse_arguments([option])
    except SystemExit:
        pass
    out, err = capsys.readouterr()
    assert re.match(r"^primalscheme \d{1,2}.\d{1,2}.\d{1,2}[a-zA-Z]{0,3}\d*$", out)


def test_force_not_required_when_output_path_does_not_exist(tmp_path):
    path = tmp_path / "output"
    primalscheme.cli.get_output_path(output_path=path)


def test_force_required_when_output_path_does_exist(tmp_path):
    path = tmp_path / "output"
    path.mkdir()
    with pytest.raises(IOError):
        primalscheme.cli.get_output_path(output_path=path)


def test_force_allows_existing_output_path(tmp_path):
    path = tmp_path / "output"
    path.mkdir()
    primalscheme.cli.get_output_path(output_path=path, force=True)


def test_existing_output_path_not_dir(tmp_path):
    path = tmp_path / "output"
    path.write_text("text")
    with pytest.raises(IOError):
        primalscheme.cli.get_output_path(output_path=path, force=True)


def test_too_short_fasta_size_vs_amplicon_size(input_fasta_short_500, tmp_path, caplog):
    args = ["multiplex", str(input_fasta_short_500)]
    parsed = primalscheme.cli.parse_arguments(args)
    outpath = tmp_path / "output"
    output_path = primalscheme.cli.get_output_path(outpath)

    with pytest.raises(SystemExit):
        primalscheme.cli.multiplex(parsed, output_path)
    messages = [log.message for log in caplog.records]
    assert any("too short" in msg for msg in messages)


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


def test_process_fasta_returns_longest_reference_first(input_fasta_shortest_first):
    references = process_fasta(input_fasta_shortest_first, sort=True)
    first_ref_len = len(references[0])
    for ref in references[1:]:
        assert len(ref) <= first_ref_len


def test_process_fasta_no_sort(input_fasta_shortest_first):
    references = process_fasta(input_fasta_shortest_first, sort=False)
    first_ref_len = len(references[0])
    for ref in references[1:]:
        assert len(ref) >= first_ref_len
