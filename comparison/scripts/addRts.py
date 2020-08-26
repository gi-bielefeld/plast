#!/usr/bin/env python3

import sys
from glob import glob

compRt = 0.0
maxMem = 0

#Get all input file names
rtfilenames = glob(sys.argv[-1].split(".txt")[0] + "_part*")

for n in rtfilenames:
	rtfile = open(n, 'r')
	
	for l in rtfile:
		if l.find("User") >= 0:
			compRt += float(l.split(' ')[3])
		elif l.find("Maximum") >= 0:
			maxMem = max(maxMem, int(l.split(' ')[5]))

	rtfile.close()

outfile = open(sys.argv[-1], 'w')
outfile.write("User time (seconds): " + str(compRt) + '\n')
outfile.write("Maximum resident set size (kbytes): " + str(maxMem))
outfile.close()