#!/usr/bin/env python

import sys
import random as r

seq = ""

for i in range(int(sys.argv[1])):
	seq += r.choice("ACGT")

print seq
