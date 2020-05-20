"""
Tests for command line interface (CLI)
"""
import os
import pytest
import re

import primalscheme.cli

from argparse import Namespace


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


def test_cli_fails_without_fasta(default_config):
    """
    Does CLI stop execution w/o a positional fasta arg?
    """
    with pytest.raises(SystemExit):
        args = ["multiplex"]
        primalscheme.cli.parse_arguments(args, default_config)


@pytest.mark.parametrize(
    "option, value",
    [
        ("--prefix", "myScheme"),
        ("--amplicon-size-min", "380"),
        ("--amplicon-size-max", "420"),
        ("--target-overlap", "0"),
        ("--min-unique", "3"),
        ("--max-candidates", "10"),
        ("--output-path", "./scheme"),
        ("--step-distance", "11"),
        ("--debug", None),
        ("--force", None),
    ],
)
def test_availability_of_cli_options(option, value, default_config):
    """
    Are options available?
    """
    args = ["multiplex", "some.fa", option]
    if value:
        args.append(value)
    parsed = primalscheme.cli.parse_arguments(args, default_config)
    assert isinstance(parsed, Namespace)


def test_cli_fails_with_invalid_option(chikv_input, default_config):
    """
    Does CLI stop execution with an invalid option?
    """
    with pytest.raises(SystemExit):
        args = ["multiplex", f"{str(chikv_input)}", "--force", "--badOption"]
        primalscheme.cli.parse_arguments(args, default_config)


@pytest.mark.parametrize(
    "option", ["-V", "--version"],
)
def test_cli_version_output(option, capsys, default_config):
    """
    Does CLI output a sensible version number for -V and --version?
    """
    try:
        primalscheme.cli.parse_arguments([option], default_config)
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


def test_too_short_fasta_size_vs_amplicon_size(
    input_fasta_short_500, default_config, tmp_path, caplog
):
    args = ["multiplex", str(input_fasta_short_500)]
    parsed = primalscheme.cli.parse_arguments(args, default_config)
    outpath = tmp_path / "output"
    output_path = primalscheme.cli.get_output_path(outpath)

    with pytest.raises(SystemExit):
        primalscheme.cli.multiplex(parsed, output_path)
    messages = [log.message for log in caplog.records]
    assert any("too short" in msg for msg in messages)
