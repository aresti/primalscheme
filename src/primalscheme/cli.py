"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the CLI for PrimalScheme.
It is executed when the user runs 'primalscheme' after installation,
or 'python3 -m primalscheme'.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""

import json
import sys
import os
import argparse
import logging

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from primalscheme.reporting import MultiplexReporter

logger = logging.getLogger('primalscheme')


def main():
    """
    Entry point for primalscheme.
    Proccess args, set up logging and call the command function.
    """

    args = get_arguments()
    output_path = args.output_path

    try:
        check_output_dir(output_path, force=args.force)
        setup_logging(output_path, debug=args.debug, prefix=args.prefix)

        for arg in vars(args):
            logger.debug('{}: {}'.format(arg, str(vars(args)[arg])))

        # Run
        args.func(args)
    except Exception as e:
        logger.error('Error: {}'.format(e))
        raise e
        sys.exit()


def multiplex(args):
    """
    Multipex scheme command.
    """

    references = process_fasta(args.fasta)

    logger.info(f'Designing primers using reference: {references[0].id}')
    for ref in references[1:]:
        logger.info(f'Checking alignments against reference: {ref.id}')

    scheme = MultiplexReporter(
        references, args.amplicon_size, args.amplicon_max_variation,
        args.target_overlap, args.step_distance, args.min_unique, args.prefix,
        args.primer3, progress_func=stdout_progress)
    scheme.design_scheme()

    scheme.write_all(args.output_path)


def process_fasta(file_path):
    """
    Parse and validate the fasta file.
    """

    references = []
    records = SeqIO.parse(file_path, 'fasta')  # may raise

    # Remove gaps
    for record in records:
        ref = SeqRecord(
            Seq(str(record.seq).replace('-', '').upper()),
            id=record.id, description=record.id
        )
        references.append(ref)

    # Check for no references
    if not references:
        raise ValueError(
            'The input FASTA file does not contain any valid references.')

    # Check for too many references
    if len(references) > 100:
        raise ValueError(
            'A maximum of 100 reference genomes is currently supported.')

    # Check for max difference in length between references
    primary_ref = references[0]
    primary_ref_len = len(primary_ref)

    if any(abs(len(r) - primary_ref_len) > 500 for r in references):
        raise ValueError(
            'One or more of your references is too different in length to '
            'the primary (first) reference. The maximum difference is '
            '500 nt.'
        )

    # Check for a valid alphabet
    VALID_ALPHABET = "AGCTacgt"
    for r in references:
        if any(l not in VALID_ALPHABET for l in set(r.seq)):
            raise ValueError(
                'One or more of your fasta sequences contain invalid '
                "nucleotide codes. The supported alphabet is '{}'. "
                'Ambiguity codes and gaps are not currently supported.'
                .format(VALID_ALPHABET)
            )

    return references


def setup_logging(output_path, debug=False, prefix='primalscheme'):
    """
    Setup logging output and verbosity.
    """

    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # File handler
    log_filepath = os.path.join(output_path, '{}.log'.format(prefix))
    fh = logging.FileHandler(log_filepath)
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    # Stream handler STDOUT
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh_formatter = logging.Formatter('%(message)s')
    sh.setFormatter(sh_formatter)
    logger.addHandler(sh)

    logger.info(f'Writing log to {log_filepath}')


def check_output_dir(output_path, force=False):
    """
    Check for an existing output dir, require --force to overwrite.
    Otherwise, create a new one.
    """

    if os.path.isdir(output_path) and not force:
        raise IOError('Directory exists add --force to overwrite')
    if not os.path.isdir(output_path):
        os.mkdir(output_path)


def get_arguments(test=[]):
    """
    Parse command line arguments, update config
    """
    # Setup parsers
    parser = argparse.ArgumentParser(
        prog='primalscheme',
        description=(
            'A primer3 wrapper for designing multiplex primer schemes.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title='[subcommands]', dest='command', required=True)

    parser_scheme = subparsers.add_parser(
        'multiplex', help='Multiplex PCR scheme')

    # Set func and defaults from config file
    config_file = Path(__file__).parent / 'config.json'
    config = json.loads(config_file.read_text())
    parser_scheme.set_defaults(func=multiplex, **config)

    parser_scheme.add_argument(
        'fasta', help='FASTA file')
    parser_scheme.add_argument(
        '--prefix', 
        help='Prefix used for primer names and output files '
             '(default: %(default)s)')
    parser_scheme.add_argument(
        '--amplicon-size', type=int,
        help='Desired amplicon size (default: %(default)i)')
    parser_scheme.add_argument(
        '--amplicon-max-variation', type=int,
        help='Maximum amplicon variation (default: %(default)i)')
    parser_scheme.add_argument(
        '--debug', action='store_true', help=f'Verbose logging')
    parser_scheme.add_argument(
        '--force', action='store_true',
        help='Force output to an existing directory and overwrite output '
             'files.')
    parser_scheme.add_argument(
        '--target-overlap', type=int,
        help='Target overlap size (default: %(default)i)')
    parser_scheme.add_argument(
        '--min-unique', type=int,
        help='Minimum unique candiate pairs (default: %(default)i)')
    parser_scheme.add_argument(
        '--max-candidates', type=int,
        help='Maximum candidate pairs (default: %(default)i)')
    parser_scheme.add_argument(
        '--output-path',
        help='Output directory (default: %(default)s)')
    parser_scheme.add_argument(
        '--step-distance', type=int,
        help='Distance to step between find attempts (default: %(default)i)')

    # Generate args
    args = parser.parse_args(test if test else None)

    # update primer3 with derivatives
    args.__dict__['primer3'].update({
        'PRIMER_NUM_RETURN': args.max_candidates,
        'PRIMER_PRODUCT_SIZE_RANGE': [[
            args.amplicon_size - args.amplicon_max_variation,
            args.amplicon_size + args.amplicon_max_variation
        ]],
    })

    return args


def stdout_progress(count, total):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write(f'\r[{bar}] {percents}%')
    if count == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


if __name__ == '__main__':
    main()
