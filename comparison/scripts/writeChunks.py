#!/usr/bin/env python3

import sys

CHUNK_SIZE = 100

colorNamesFile = open(sys.argv[1], 'r')
colorNamesLine = colorNamesFile.readline()
colorNamesFile.close()

colorCounter = 0

for n in colorNamesLine.split(' ')[:-1]:
	if colorCounter % CHUNK_SIZE == 0:
		if colorCounter > 0:
			outfile.close()

		outfileName = '/'.join(sys.argv[1].split('/')[:-1]) + "/merged_" + sys.argv[1].split(".txt")[0].split("chosen_")[1] + "_part" + str(int(colorCounter / CHUNK_SIZE) + 1) + ".fasta"
		outfile = open(outfileName, 'w')

	fastaFile = open(n.replace("fasta", "fastatmp").replace('\n',''), 'r')

	for l in fastaFile:
		outfile.write(l)

	fastaFile.close()
	colorCounter += 1