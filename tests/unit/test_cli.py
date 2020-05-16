"""
Tests for command line interface (CLI)
"""
import os
import pytest

import primalscheme.cli

from argparse import Namespace


def test_runas_module():
    """
    Can this package be run as a Python module?
    """
    exit_status = os.system('python3 -m primalscheme --help')
    assert exit_status == 0


def test_entrypoint():
    """
    Is entrypoint script installed? (setup.py)
    """
    exit_status = os.system('primalscheme --help')
    assert exit_status == 0


def test_mutliplex_command():
    """
    Is multiplex command available?
    """
    exit_status = os.system('primalscheme multiplex --help')
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
    exit_status = os.system('primalscheme multiplex')
    assert exit_status != 0


@pytest.mark.parametrize('option, value', [
    ('--prefix', 'myScheme'),
    ('--amplicon-size-min', '380'),
    ('--amplicon-size-max', '420'),
    ('--target-overlap', '0'),
    ('--min-unique', '3'),
    ('--max-candidates', '10'),
    ('--output-path', './scheme'),
    ('--step-distance', '11'),
    ('--debug', None),
    ('--force', None),
])
def test_availability_of_options(option, value, default_config):
    """
    Are options available?
    """
    args = ['multiplex', 'some.fa', option]
    if value:
        args.append(value)
    parsed = primalscheme.cli.parse_arguments(args, default_config)
    assert isinstance(parsed, Namespace)


def test_cli_fails_with_invalid_option(chikv_input):
    """
    Does CLI stop execution with an invalid option?
    """
    exit_status = os.system(
        f'primalscheme multiplex {str(chikv_input)} --force --badOption')
    assert exit_status != 0
