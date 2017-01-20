
outer_params = {
	'PRIMER_OPT_SIZE': 22,
	'PRIMER_MIN_SIZE': 22,
	'PRIMER_MAX_SIZE': 30,
	'PRIMER_OPT_TM': 61.5,
	'PRIMER_MIN_TM': 60.0,
	'PRIMER_MAX_TM': 63.0,
	'PRIMER_MIN_GC': 30.0,
	'PRIMER_MAX_GC': 55.0,
	'PRIMER_MAX_POLY_X': 5,
	'PRIMER_SALT_MONOVALENT': 50.0,
	'PRIMER_DNA_CONC': 50.0,
	'PRIMER_MAX_NS_ACCEPTED': 0,
	'PRIMER_MAX_SELF_ANY': 8,
	'PRIMER_MAX_SELF_END': 47,
	'PRIMER_PAIR_MAX_COMPL_ANY': 8,
	'PRIMER_PAIR_MAX_COMPL_END': 47,
}

MATCHES = [
	set(['A', 'T']),
	set(['C', 'G']),
	set(['G', 'T']),
	set(['C', 'T']),
	set(['T', 'T'])
	 ],

MISMATCHES = [
	set(['A', 'A']),
	set(['A', 'C']),
	set(['C', 'C']),
	set(['G', 'A']),
	set(['G', 'G']),
	]

STEP_DISTANCE = 10

# Length 22, 30, 24.43
# Tm 59.96, 62.72, 61.33
# GC 33.33, 54.55, 45.71
# Longest homopolymer 5
# Ambiguous 0
