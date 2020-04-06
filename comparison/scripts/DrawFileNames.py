#!/usr/bin/python3

import sys
from random import choice

#This script takes as input a number n and a bunch of files and randomly outputs the names of n of these files
#Parameters:
#	1. Number of files to select
#	2... File names

#Get number of input parameters
numInp = len(sys.argv)

#Check if we have enough input parameters
if numInp < 3:
	print("Too few input arguments", file=sys.stderr)

#Get n
n = int(sys.argv[1])

#Check if we have enough files to draw from
if numInp - 2 < n:
	print("Number of files to be selected is not smaller than the number of file names to select from! Is that really what you want?", file=sys.stderr)
	exit(-1)

#Get list of all file names
fNames = sys.argv[2:]

#Draw n times
for i in range(n):
	#Draw a name
	ch = choice(fNames)
	#Delete file name from list
	fNames.remove(ch)
	#Output file name
	print(ch, end=' ')
	print(ch, end=' ', file=sys.stderr)