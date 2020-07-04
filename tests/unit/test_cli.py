"""
Tests for command line interface (CLI)
"""

import os
import pytest
import re
import warnings

from Bio.SeqRecord import SeqRecord
from primalscheme.cli import cli

with warnings.catch_warnings():
    warnings.simplefilter("ignore", PendingDeprecationWarning)
    from primalscheme.cli import process_fasta


def test_runas_module():
    """Can this package be run as a Python module?"""
    exit_status = os.system("python -m primalscheme --help")

    assert exit_status == 0


def test_entrypoint():
    """Is entrypoint script installed?"""
    exit_status = os.system("primalscheme --help")

    assert exit_status == 0


def test_mutliplex_command(cli_runner):
    """Is multiplex command available, with help?"""
    result = cli_runner.invoke(cli, ["multiplex", "-h"])

    assert result.exit_code == 0


def test_cli_without_command_shows_help(cli_runner):
    """Does CLI display help when run w/o a command argument?"""
    result = cli_runner.invoke(cli)

    assert "--help" in result.output


def test_multiplex_command_without_fasta_returns_exit_code_2(cli_runner):
    """Does CLI stop execution w/o a positional fasta arg, returning exit code 2?"""
    result = cli_runner.invoke(cli, ["multiplex"])

    assert result.exit_code == 2


@pytest.mark.parametrize(
    "option, value",
    [
        ("--amplicon-size", "380"),
        ("-a", "420"),
        ("--outpath", "./scheme"),
        ("-o", "./scheme"),
        ("--name", "myScheme"),
        ("-n", "myScheme"),
        ("--target-overlap", "0"),
        ("-t", "0"),
        ("--debug", None),
        ("-d", None),
        ("--force", None),
        ("-f", None),
        ("--pinned", None),
        ("-p", None),
        ("--high-gc", None),
        ("-g", None),
        ("--help", None),
        ("-h", None),
    ],
)
def test_availability_of_cli_options(option, value, cli_runner):
    """Are options available?"""
    args = ["multiplex", "-h", option]
    if value:
        args.append(value)
    result = cli_runner.invoke(cli, args)

    assert "Error: no such option" not in result.output


def test_cli_fails_with_invalid_option(chikv_input, cli_runner):
    """Does CLI stop execution with an invalid option?"""
    args = ["multiplex", f"{str(chikv_input)}", "--force", "--badOption"]
    result = cli_runner.invoke(cli, args)

    assert result.exit_code == 2
    assert "Error: no such option" in result.output


@pytest.mark.parametrize(
    "option", ["-V", "--version"],
)
def test_cli_version_output(option, cli_runner):
    """Does CLI output a sensible version number for -V and --version?"""
    result = cli_runner.invoke(cli, option)

    assert re.match(
        r"^cli, version \d{1,2}.\d{1,2}.\d{1,2}[a-zA-Z]{0,3}\d*$", result.output,
    )


def test_force_not_required_when_outpath_does_not_exist(
    tmp_path, cli_runner, chikv_input
):
    """Does CLI proceed without --force when outpath does not exist?"""
    path = str(tmp_path / "output")
    result = cli_runner.invoke(cli, ["multiplex", str(chikv_input), "-o", path, "-p"])

    assert result.exit_code == 0


def test_force_required_when_output_path_does_exist(tmp_path, cli_runner, chikv_input):
    """Does CLI require --force when outpath does exist?"""
    path = tmp_path / "output"
    path.mkdir()
    result = cli_runner.invoke(
        cli, ["multiplex", str(chikv_input), "-o", str(path), "-p"]
    )

    assert result.exit_code != 0
    assert "--force" in result.output


def test_force_allows_existing_output_path(tmp_path, cli_runner, chikv_input):
    """Does CLI proceed with an existing outpath when --force is set?"""
    path = tmp_path / "output"
    path.mkdir()
    result = cli_runner.invoke(
        cli, ["multiplex", str(chikv_input), "-o", str(path), "-p", "-f"]
    )

    assert result.exit_code == 0


def test_existing_output_path_not_dir(tmp_path, cli_runner, chikv_input):
    """Does CLI require outpath to be a directory"""
    path = tmp_path / "output"
    path.write_text("text")
    result = cli_runner.invoke(
        cli, ["multiplex", str(chikv_input), "-o", str(path), "-p"]
    )

    assert result.exit_code != 0
    assert "file" in result.output


def test_warning_on_high_gc_input(tmp_path, cli_runner, high_gc_input):
    """Does CLI warn user when using a high-GC input without the --high-gc option?"""
    path = tmp_path / "output"
    result = cli_runner.invoke(cli, ["multiplex", str(high_gc_input), "-o", str(path)])

    assert "warning:" in result.output.lower()
    assert "high-gc" in result.output.lower()


def test_warning_silenced_when_using_high_gc(tmp_path, cli_runner, high_gc_input):
    """Does CLI silence the warning when using the --high-gc option?"""
    path = tmp_path / "output"
    result = cli_runner.invoke(
        cli, ["multiplex", str(high_gc_input), "-o", str(path), "-g"]
    )

    assert "warning:" not in result.output.lower()
    assert "high-gc" not in result.output.lower()


def test_no_high_gc_warning_on_average_input(
    tmp_path, cli_runner, ncov_single_ref_input
):
    """Does CLI stay quiet when using a regular input wihout the --high-gc option?"""
    path = tmp_path / "output"
    result = cli_runner.invoke(
        cli, ["multiplex", str(ncov_single_ref_input), "-o", str(path)]
    )

    assert "warning:" not in result.output.lower()
    assert "high-gc" not in result.output.lower()


def test_too_short_fasta_size_vs_amplicon_size(
    input_fasta_short_500, tmp_path, cli_runner
):
    """Does CLI exit where amplicon size is greater than reference size?"""
    path = str(tmp_path / "output")
    args = [
        "multiplex",
        str(input_fasta_short_500),
        "-a",
        "500",
        "-a",
        "550",
        "-o",
        path,
    ]
    result = cli_runner.invoke(cli, args)

    assert result.exit_code != 0
    assert "too short" in result.output


def test_process_fasta_empty_input(input_fasta_empty):
    """Does process_fasta raise for empty input?"""
    with pytest.raises(ValueError, match="does not contain any valid references"):
        process_fasta(input_fasta_empty)


def test_process_fasta_invalid_alphabet(input_fasta_invalid_alphabet):
    """Does process_fasta raise for invalid alphabet?"""
    with pytest.raises(ValueError, match="invalid nucleotide codes"):
        process_fasta(input_fasta_invalid_alphabet)


def test_process_fasta_too_many_records(input_fasta_101_random_valid):
    """Does process_fasta raise for too many references?"""
    with pytest.raises(ValueError, match="A maximum of 100"):
        process_fasta(input_fasta_101_random_valid)


def test_process_fasta_size_difference_over_500(input_fasta_size_difference_over_500,):
    """Does process_fasta raise for size difference over 500?"""
    with pytest.raises(ValueError, match="too different in size"):
        process_fasta(input_fasta_size_difference_over_500)


def test_process_fasta_returns_list_of_seq_records(input_fasta_5_random_valid):
    """Does process_fasta return a list of SeqRecord objects?"""
    references = process_fasta(input_fasta_5_random_valid)
    assert len(references) == 5
    assert all(isinstance(r, SeqRecord) for r in references)


def test_process_fasta_chikv_demo(chikv_input):
    """Does process_fasta return 2 SeqRecord objects for demo CHIKV input?"""
    references = process_fasta(chikv_input)
    assert len(references) == 2


def test_process_fasta_gaps_removed(input_fasta_valid_with_gaps):
    """Does process_fasta remove gaps?"""
    references = process_fasta(input_fasta_valid_with_gaps)
    assert "-" not in references[0]


def test_process_fasta_too_short_input(input_fasta_short_500):
    """Does process_fasta raise for too short input?"""
    with pytest.raises(ValueError, match="too short"):
        process_fasta(input_fasta_short_500, min_ref_size=900)
