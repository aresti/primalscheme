#!/usr/bin/env python2.7

import sys
import argparse
import pprint

def run_command(parser, args):
	if args.command == 'scheme':
		import scheme as command
	return command.multiplex(args, parser=parser)

def main():
	parser = argparse.ArgumentParser(prog='primal', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

	#scheme
	parser_scheme = subparsers.add_parser('scheme', help='Tiling amplicons designer')
	parser_scheme.add_argument('-g', help='FASTA file', metavar='STRING', required=True)
	parser_scheme.add_argument('-o', help='Prefix', metavar='STRING', required=True)
	parser_scheme.add_argument('--length', help='Amplicon length', type=int, default=400)
	parser_scheme.add_argument('--overlap', help='Overlap length', type=int, default=75)
	parser_scheme.set_defaults(func=run_command)

	#run
	args = parser.parse_args()
	result = args.func(parser, args)


if __name__ == '__main__':
	main()
