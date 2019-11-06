#!/usr/bin/env python

import sys

#This function adds a new score count to a dictionary
def pushScore(score, sDict):
	#Check if score value is already inside the dictionary
	if score in sDict:
		#Increment counter
		sDict[score] += 1
	else:
		#Add new entry
		sDict[score] = 1

#This function outputs found scores and their abundances
def reportHit(sDict):
	#Get entries in ungapped dictionary
	entries = sDict.keys()
	#Sort scores decreasingly
	entries.sort()

	#Go through all dictionary entries
	for e in entries:
		#Output entry
		print e, sDict[e]

if __name__ == '__main__':
	#Dictionary with max ungapped score-count key-value pairs
	ugCnts = {}
	#Dictionary with max gapped score-count key-value pairs
	gCnts = {}

	#Read each file
	for f in sys.argv[1:]:
		#Open file
		sFile = open(f, 'r')
		#Initialize maximum scores
		maxScoreUg = 0
		maxScoreG = 0

		#Read each line
		for l in sFile:
			#Check if line introduces results for new query which is not the first
			if l.startswith("Query") and maxScoreUg != 0:
				#Add maximum scores to dictionaries
				pushScore(maxScoreUg, ugCnts)
				pushScore(maxScoreG, gCnts)
				#Reset maximum scores
				maxScoreUg = 0
				maxScoreG = 0
			#Check if current line decodes an ungapped score
			elif "ungapped" in l:
				#Update maximum score
				maxScoreUg = max(maxScoreUg, int(l.split(' ')[2]))
			#Check if current line decodes a gapped score
			elif "(gapped" in l:
				#Update maximum score
				maxScoreG = max(maxScoreG, int(l.split(' ')[2]))
			#Every other line can be disregarded
			else:
				continue

		#Add maximum scores to dictionaries a last time
		pushScore(maxScoreUg, ugCnts)
		pushScore(maxScoreG, gCnts)
		#Close file
		sFile.close()

	#Output scores
	print "Ungapped results:"
	reportHit(ugCnts)
	print "Gapped results:"
	reportHit(gCnts)