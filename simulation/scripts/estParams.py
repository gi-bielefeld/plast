#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from matplotlib.patches import Rectangle

#This function changes counts to normalized cumulative counts. It is expected that values are decreasingly sorted)
def chNormCum(counts):
	#Get total number of counts
	nbCnts = sum(counts[1])

	for c in range(1,len(counts[1])):
		#Add the number of larger scores to the next smaller one
		counts[1][c] += counts[1][c - 1]
		#Normalize previous entry
		counts[1][c - 1] /= nbCnts

	#Normalize lowest score count as well
	counts[1][len(counts[1]) - 1] /= nbCnts

#The line function we want to fit
def lineFunc(x, lmd, c):
	return -lmd * x + c

#The function describing the expected random distribution
def func(x, a, c):
	return c * np.exp(-a * x)

if __name__ == '__main__':
	#A flag to distinguish gapped and ungapped results read in
	isGapped = False
	#A list for ungapped results 1. score value 2. count
	ugRes = [[], []]
	#A list for gapped results 1. score value 2. count
	gRes = [[], []]
	#A dictionary harboring a results line shape mapping
	lineShapes = {"ungapped": "--", "gapped": "-."}
	#A dictionary harboring a results dot shape mapping
	dotstyles = {"ungapped": 'x', "gapped": '+'}

	#Open input file
	file = open(sys.argv[1], 'r')
	#The first line in the file is not of interest for us
	file.readline()

	#Read each remaining line
	for l in file:
		#Check if we have reached the gapped results
		if l.startswith("Gapped"):
			#Mark that gapped results will follow
			isGapped = True
			continue

		#Check if we still handle ungapped results
		if not isGapped:	
			#Parse line
			elems = [float(e) for e in l.split(' ')]
			#Add score value pair
			ugRes[0].append(elems[0])
			ugRes[1].append(elems[1])
		else:
			#Parse line
			elems = [int(e) for e in l.split(' ')]
			#Add score value pair
			gRes[0].append(elems[0])
			gRes[1].append(elems[1])

	#Change counts to normalized cumulative counts
	chNormCum(ugRes)
	chNormCum(gRes)

	fig, ax = plt.subplots()

	#Repeat for both results
	for r in [["ungapped", ugRes], ["gapped", gRes]]:
		#Flags indicating if considerable interval for line fitting have already been found
		smlValFound = False
		lrgValFound = False

		#Go through counts
		for i in range(len(r[1][1])):
			#If the first score outside of our fitting interval is found we also know the last one inside
			if r[1][1][i] > 0.01 and not smlValFound:
				#Save last score inside the interval
				smlVal = r[1][0][i - 1]
				#Mark smallest score as found
				smlValFound = True
			#Check if the first score inside our fitting interval is reached
			elif r[1][1][i] >= 0.0001 and not lrgValFound:
				#Save first score inside the interval
				lrgVal = r[1][0][i]
				#Mark largest score as found
				lrgValFound = True

		#Get smallest index (corresponding to largest score as values are decreasingly sorted)
		smlst = r[1][0].index(lrgVal)
		#Get largest index
		lrgst = r[1][0].index(smlVal)
		#Fit line
		popt, pcov = curve_fit(lineFunc, r[1][0][smlst:lrgst], np.log(r[1][1][smlst:lrgst]))
		#Output parameters
		print("Lambda (%s):" %r[0], popt[0], "C (%s):" %r[0], popt[1])
		#Calculate fitting errors
		perr = np.sqrt(np.diag(pcov))
		#Output errors
		print("Errors:", perr[0], perr[1])
		#Calculate lines y-values
		ys = [func(i, popt[0], np.exp(popt[1])) for i in r[1][0]]
		#Plot line
		plt.plot(r[1][0], ys, color = 'tab:gray', linestyle = lineShapes[r[0]])
		#Plot score value pairs
		plt.semilogy(r[1][0], r[1][1], 'k' + dotstyles[r[0]], label = r[0])

	#Define highlight for fitted region
	rect = Rectangle((10,0.0001),80,0.01,linewidth=1,edgecolor='none',facecolor='lightgray')
	#Enable a grid
	plt.grid(True)
	#Add the highlight
	ax.add_patch(rect)
	#Add a legend
	plt.legend()
	#Label x-axis
	plt.xlabel("Score")
	#Save figure
	plt.savefig("results/dist.pdf", format="pdf")