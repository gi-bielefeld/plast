#!/usr/bin/env python3

import argparse as args
from numpy import mean, std

#Setting up the argument parser
parser = args.ArgumentParser(description="Show statistics for PLAST results")
parser.add_argument('R', metavar='ResultFile', type=args.FileType('r'), help="PLAST result file")
# parser.add_argument('-H', metavar='PlotHist', type=str, choices=["nres", "lens", "evals"], help="Plot histogram of specified property")
arguments = parser.parse_args()

#Parse data
nbRes = 0
scores = []
lengths = []
evals = []

for l in arguments.R:
	if l.startswith("Score:"):
		nbRes += 1
		elems = l.split('\t')
		scores.append(int(elems[0].split(' ')[1]))
		lengths.append(int(elems[1].split(' ')[1]))
		evals.append(float(elems[2].split(' ')[1]))

#Output descriptive stats
print("Number of found results:", nbRes)
print(f"Results scores: max: {max(scores)} min: {min(scores)} avg.: {mean(scores)} std.: {std(scores)}")
print(f"Results lengths: max: {max(lengths)} min: {min(lengths)} avg.: {mean(lengths)} std.: {std(lengths)}")
print(f"Results e-values: max: {max(evals)} min: {min(evals)} avg.: {mean(evals)} std.: {std(evals)}")
