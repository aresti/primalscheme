#!/usr/bin/env python2.7

import sys
import argparse
from pprint import pprint

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
	parser_scheme.add_argument('--verbose', help='Verbose mode', action="store_true")
	parser_scheme.set_defaults(func=run_command)

	#run
	args = parser.parse_args()
	result = args.func(parser, args)
	
if __name__ == '__main__':
	main()


"""
for each in result:
	left = each.pairs[0].left
	right = each.pairs[0].right
print '\t'.join(['gi|992324757|gb|KU707826.1|', str(left.start), str(left.end), left.name, '0', '+'])
print '\t'.join(['gi|992324757|gb|KU707826.1|', str(right.end), str(right.start), right.name, '0', '-'])
print '\t'.join(['CHIK_', left.name, left.seq, '100nm', 'STD'])
print '\t'.join(['CHIK_', right.name, right.seq, '100nm', 'STD'])
print '\t'.join(['ZIKA_'+ left.name, left.seq, str(each.pool), str(left.length), str(left.tm), str(left.gc)])
print '\t'.join(['ZIKA_'+ right.name, right.seq, str(each.pool), str(right.length), str(right.tm), str(right.gc)])
"""
