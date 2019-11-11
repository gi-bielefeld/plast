#!/usr/bin/env python3

import sys
import random as r

#Parameter 1: Length of random sequence
#Parameter 2: Number of sequences

#For each sequence to be generated
for i in range(int(sys.argv[2])):
	#Initialize it
	seq = ""

	#Draw a base for each sequence position
	for i in range(int(sys.argv[1])):
		seq += r.choice("ACGT")

	#Output sequence
	print(seq)
