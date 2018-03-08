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

from pprint import pprint
from primal import multiplex
from primal import smart


def main():
    parser = argparse.ArgumentParser(prog='primal', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    #Standard scheme
    parser_scheme = subparsers.add_parser('scheme', help='Tiling amplicons designer')
    parser_scheme.add_argument('fasta', help='FASTA file')
    parser_scheme.add_argument('prefix', help='Prefix')
    parser_scheme.add_argument('--amplicon-length', help='Amplicon length', type=int, default=400)
    parser_scheme.add_argument('--min-overlap', help='Minimum overlap length', type=int, default=20)
    parser_scheme.add_argument('--max-gap', help='Maximum gap to introduce before failing', type=int, default=100)
    parser_scheme.add_argument('--max-candidates', help='Maximum candidate primers', type=int, default=10)
    parser_scheme.add_argument('--search-space', help='Initial primer search space', type=int, default=40)
    parser_scheme.add_argument('--output-path', help='Output directory to save files', default='./')
    parser_scheme.add_argument('--force', help='Force overwrite', action="store_true")
    parser_scheme.add_argument('--debug', help='Verbose logging', action="store_true")
    parser_scheme.set_defaults(func=multiplex)

    #Smart scheme
    parser_smart = subparsers.add_parser('smart', help='Multiplex Smart designer')
    parser_smart.add_argument('fasta', help='FASTA file')
    parser_smart.add_argument('prefix', help='Prefix')
    parser_smart.set_defaults(func=smart)

    #Generate args
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

    #Log
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(os.path.join(args.output_path, '{}.log'.format(args.prefix)))
    logger.addHandler(handler)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.setLevel(logging.DEBUG)
    logger.debug('Logging started...')
    for arg in vars(args):
        logger.info('%s:%s' %(arg, str(vars(args)[arg])))

    #Run
    args.func(args, logger)

if __name__ == '__main__':
    main()
