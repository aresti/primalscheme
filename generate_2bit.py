import subprocess
import os

def make_2bit(file, out):
	p = subprocess.Popen(['faToTwoBit', file, out)
