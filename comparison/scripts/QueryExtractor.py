import random as r
import argparse as args
import sys
from Bio import SeqIO

if __name__ == '__main__':

	#Setting up the argument parser
	parser = args.ArgumentParser(description="This function iterates over the contigs of an assembly and cuts out a query of a specified length")
	parser.add_argument('Q', metavar='QueryLength', type=int, help="Length of query to be cut out")
	parser.add_argument('A', metavar='AssemblyFiles', type=str, nargs='+', help="Paths to fasta files containing the assembled sequences")

	arguments = parser.parse_args()

	#List to store all contig sequences whose length is at least the required query length
	posContSeqs = []
	#List to store the number of possible start positions for each conitg and all contig's in the list in front of it
	startPos = []
	#The total number of possible start positions found so far
	posStartNum = 0

	#Iterate over input files
	for f in arguments.A:
		#Parse the fasta file
		records = SeqIO.parse(f, 'fasta')

		#Go through the file
		for record in records:
			#Make sure that the contig sequence is long enough to extract a query of required length
			if len(record) > arguments.Q:
				#Add the sequence to the list of possible sequences (without its line break at the end)
				posContSeqs.append(record.seq)
				#Update the total number of possible start positions
				posStartNum += len(record) - arguments.Q
				#Create an entry for this contig in startPos
				startPos.append(posStartNum)

	#Check if there is no contig that is long enough to extract a query sequence
	if posStartNum == 0:
		#Put an error message
		print("Assembly does not provide a contig large enough to extract a query of required length from it", file=sys.stderr)
		#Terminate the program
		exit(1)

	#Decide where to extract the query sequence
	extrPos = r.randint(0, posStartNum - 1)

	#Testing
	#print "Extraction position is", extrPos

	#Go through startPos
	for i in range(len(startPos)):
		#Check whether the selected start position lies within the range of the current contig
		if startPos[i] > extrPos:
			#Check whether the start position lies in the list's first contig
			if i == 0:
				#In this case the relative start position is also the absolute start position
				relStartPos = extrPos
			else:
				#Calculate the relative start position inside the unitig
				relStartPos = extrPos - startPos[i - 1]

			#Output the query sequence
			print(posContSeqs[i][relStartPos:relStartPos + arguments.Q])
			#Terminate the program
			exit(0)
				



