import subprocess
import os

def make_2bit(file, out):
	try:	
		p = subprocess.Popen(['faToTwoBit', file, out])
	except Exception:
		'Print add faToTwoBit to path'

if __name__ == '__main__':
	make_2bit('x', 'y')
