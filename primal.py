#!/usr/bin/env python2.7
# Primal scheme by Josh Quick and Andy Smith 2016
# www.github.com/aresti/primalrefactor.git

import sys
import os
import argparse
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from primal.models import MultiplexScheme


def multiplex(args):
    scheme = MultiplexScheme(args.references, args.amplicon_length, min_overlap=args.min_overlap, max_gap=args.max_gap,
                             search_space=args.search_space, max_candidates=args.max_candidates, prefix=args.prefix)
    scheme.write_bed(args.output_path)
    scheme.write_pickle(args.output_path)
    scheme.write_tsv(args.output_path)
    scheme.write_refs(args.output_path)
    scheme.write_schemadelica_plot(args.output_path)


def main():
    logger = logging.getLogger('Primal Log')
    parser = argparse.ArgumentParser(prog='primal', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    # Standard scheme
    parser_scheme = subparsers.add_parser('scheme', help='Tiling amplicons designer')
    parser_scheme.add_argument('fasta', help='FASTA file')
    parser_scheme.add_argument('prefix', help='Prefix')
    parser_scheme.add_argument('--amplicon-length', help='Amplicon length (default: %(default)i)', type=int, default=400)
    parser_scheme.add_argument('--min-overlap', help='Minimum overlap length (default: %(default)i)', type=int, default=0)
    parser_scheme.add_argument('--max-gap', help='Maximum gap to introduce before failing (default: %(default)i)', type=int, default=100)
    parser_scheme.add_argument('--max-candidates', help='Maximum candidate primers (default: %(default)i)', type=int, default=10)
    parser_scheme.add_argument('--search-space', help='Initial primer search space (default: %(default)i)', type=int, default=40)
    parser_scheme.add_argument('--output-path', help='Output directory to save files (default: %(default)s)', default='./')
    parser_scheme.add_argument('--force', help='Force overwrite', action="store_true")
    parser_scheme.add_argument('--debug', help='Verbose logging', action="store_true")
    parser_scheme.set_defaults(func=multiplex)

    # Generate args
    args = parser.parse_args()
    args.references = []
    for record in SeqIO.parse(open(args.fasta, 'r'), 'fasta'):
        args.references.append(SeqRecord(Seq(str(record.seq).replace('-', '').upper()), id=record.id, description=record.id))

    # Check directory exists
    if os.path.isdir(args.output_path) and not args.force:
        logger.error('Directory exists add --force to overwrite')
        raise IOError('Directory exists add --force to overwrite')
        sys.exit()
    if not os.path.isdir(args.output_path):
        os.mkdir(args.output_path)

    # Logging
    logger.setLevel(logging.DEBUG if args.debug else logging.INFO)

    fh = logging.FileHandler(os.path.join(args.output_path, '{}.log'.format(args.prefix)))
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG)
    sh_formatter = logging.Formatter('%(message)s')
    sh.setFormatter(sh_formatter)
    logger.addHandler(sh)

    logger.info('Primal scheme started...)')
    for arg in vars(args):
        logger.debug('{}: {}'.format(arg, str(vars(args)[arg])))

    for r in args.references:
        logger.info('Reference: {}'.format(r.id))

    # Run
    args.func(args)


if __name__ == '__main__':
    main()
