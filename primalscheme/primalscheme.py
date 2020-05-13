"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes
Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the main script for PrimalScheme.
It is executed when the user runs 'primalscheme' after installation,
or 'primal.py' (directly from the source directory).

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

from primalscheme.multiplex_reporting import MultiplexReporter

logger = logging.getLogger('primalscheme')


def main():
    """
    Entry point for primalscheme. Proccesses args, sets up logging and calls
    the appropriate scheme function.
    """

    config = get_config()
    args = get_arguments(config['scheme'])
    output_path = args.output_path

    try:
        check_output_dir(output_path, force=args.force)
        setup_logging(output_path, debug=args.debug, prefix=args.prefix)
        logger.info('PrimalScheme started...')

        for arg in vars(args):
            logger.debug('{}: {}'.format(arg, str(vars(args)[arg])))

        # Run
        args.func(config, args)  # TODO get rid of subparser
    except Exception as e:
        logger.error('Error: {}'.format(e))
        raise e
        sys.exit()


def multiplex(config, args):
    """
    Multipex scheme runner.
    """

    references = process_fasta(args.fasta)

    for i, record in enumerate(references):
        logger.info(
            f'{"PRIMARY " if i == 0 else ""}Reference {i + 1}: {record.id}')

    # update p3_global from args
    p3_global = config['primer3']
    p3_global['PRIMER_NUM_RETURN'] = args.max_candidates
    p3_global['PRIMER_PRODUCT_SIZE_RANGE'] = [[
        args.amplicon_size - args.amplicon_max_variation,
        args.amplicon_size + args.amplicon_max_variation
    ]]

    scheme = MultiplexReporter(
        references, args.amplicon_size, args.amplicon_max_variation,
        args.target_overlap, args.step_distance, args.min_unique, args.prefix,
        p3_global, progress_func=stdout_progress)
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
        raise ValueError('The input FASTA file does not contain any valid references.')

    # Check for too many references
    if len(references) > 100:
        raise ValueError('A maximum of 100 reference genomes is currently supported.')

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

    fh = logging.FileHandler(
        os.path.join(output_path, '{}.log'.format(prefix)))
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh_formatter = logging.Formatter('%(message)s')
    sh.setFormatter(sh_formatter)
    logger.addHandler(sh)


def check_output_dir(output_path, force=False):
    """
    Check for an existing output dir, require --force to overwrite.
    Otherwise, create a new one.
    """

    if os.path.isdir(output_path) and not force:
        raise IOError('Directory exists add --force to overwrite')
    if not os.path.isdir(output_path):
        os.mkdir(output_path)


def get_config():
    config_file = Path(__file__).parent / 'config.json'
    return json.loads(config_file.read_text())

    
def get_arguments(defaults):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        prog='primalscheme',
        description='A primer3 wrapper for designing multiplex primer schemes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='[subcommands]', dest='command', required=True)

    # Standard multiplex scheme
    parser_scheme = subparsers.add_parser(
        'multiplex', help='Multiplex PCR scheme')
    parser_scheme.add_argument(
        'fasta', help='FASTA file')
    parser_scheme.add_argument(
        '--prefix', default=defaults['prefix'],
        help='Prefix for primer names (default: %(default)s)')
    parser_scheme.add_argument(
        '--amplicon-size', type=int, default=defaults['amplicon_size'],
        help='Desired amplicon size (default: %(default)i)')
    parser_scheme.add_argument(
        '--amplicon-max-variation', type=int,
        default=defaults['amplicon_max_variation'],
        help='Maximum amplicon variation (default: %(default)i')
    parser_scheme.add_argument(
        '--debug', action='store_true', help=f'Verbose logging')
    parser_scheme.add_argument(
        '--force', action='store_true', help=f'Force overwrite')
    parser_scheme.add_argument(
        '--target-overlap', type=int, default=defaults['target_overlap'],
        help='Target overlap size (default: %(default)i)')
    parser_scheme.add_argument(
        '--min-unique', type=int, default=defaults['min_unique'],
        help='Minimum unique candiate pairs (default: %(default)i)')
    parser_scheme.add_argument(
        '--max-candidates', type=int, default=defaults['max_candidates'],
        help='Maximum candidate pairs (default: %(default)i)')
    parser_scheme.add_argument(
        '--output-path', default='./output',
        help='Output directory (default: %(default)s)')
    parser_scheme.add_argument(
        '--step-distance', type=int, default=defaults['step_distance'],
        help='Distance to step between find attempts (default: %(default)i)')
    parser_scheme.set_defaults(func=multiplex)

    # Generate args
    args = parser.parse_args()
    return args


def stdout_progress(count, total):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write(f'[{bar}] {percents}%' +
                     ('\n' if count == total else '\r'))
    sys.stdout.flush()


if __name__ == '__main__':
    main()
