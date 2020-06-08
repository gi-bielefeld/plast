#!/usr/bin/env python

import sys

#All infos are saved in a dictionary (key=(#colors, subsetID) value=[userTime, memory])
data = {}

#Read in files
for f in sys.argv[2:]:
	rtFile = open(f, 'r')

	#Get infos from file name
	if f.find("_s") >= 0:
		nbColors = int(f.split('c')[1].split("_s")[0])
	else:
		nbColors = int(f.split('c')[1].split('.')[0])

	subsetInfo = f.split('.')[0].split('sub')

	if len(subsetInfo) < 2:
		subsetID = 0
	else:
		subsetID = int(subsetInfo[1])

	#Parse relevant information
	for l in rtFile:
		if l.find("User") >= 0:
			data[(nbColors, subsetID)] = [float(l.split(' ')[3])]
		elif l.find("Maximum r") >= 0:
			data[(nbColors, subsetID)].append(int(l.split(' ')[5]))

	rtFile.close()

#Output data
keys = sorted(data.keys())

#Print head line
print("Prog m Quorum Colors UserTime Memory subset")

for k in keys:
	print(sys.argv[1] + " 0 0 " + str(k[0]) + " " + ' '.join([str(v) for v in data[k]]) + " " + str(k[1]))