#!/usr/bin/env python2.7

import sys
import argparse
from Bio import SeqIO

from pprint import pprint
from primal import multiplex


def main():
	parser = argparse.ArgumentParser(prog='primal', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

	#scheme
	parser_scheme = subparsers.add_parser('scheme', help='Tiling amplicons designer')
	parser_scheme.add_argument('-f', help='FASTA file', required=True)
	parser_scheme.add_argument('-p', help='Prefix', required=True)
	parser_scheme.add_argument('--amplicon-length', help='Amplicon length', type=int, default=400)
	parser_scheme.add_argument('--min-overlap', help='Minimum overlap length', type=int, default=20)
	parser_scheme.add_argument('--search-space', help='Initial primer search space', type=int, default=40)
	parser_scheme.add_argument('--output-path', help='Output path (dir) for bed and image files', default='./')
	parser_scheme.add_argument('--v', help='Verbose mode', action="store_true")
	parser_scheme.add_argument('--vvv', help='Very verbose mode', action="store_true")
	parser_scheme.set_defaults(func=multiplex)

	#run
	args = parser.parse_args()
	args.references = list(SeqIO.parse(open(args.f, 'r'), 'fasta'))
	args.func(args)

if __name__ == '__main__':
	main()
