#!/usr/bin/env bash

#This script estimates alignment statistic parameters for a specified pangenome graph

#Constants

#Path to the pangenome graph
GRAPH_PATH=$1
#k-mer length used for graph
K=31
#Minimizer length used for graph
MINI_LEN=23
#Minimal seed length used for index
SEED_LEN=11
#Number of queries used for simulation
NB_QUERIES=1000000
#Query lengths
Q_LEN=1000
#Number of queries per file
BATCH_SIZE=10000
#PLAST binary path
PLAST_BIN="../src/PLAST"

#Functions

usage()
{
	echo "usage: estParams.sh graphFile.gfa [-k k-mer-length] [-g minimizer-length] [-w seed-length] [-h]"
}

#Main

#Check if path to graph was given
if [ "$1" == "" ]; then
	#Print help
	usage
	exit 1
fi

#Parse input arguments
while [ "$2" != "" ]; do
	case $2 in
		-k )	shift
				K=$2
				;;
		-g )	shift
				MINI_LEN=$2
				;;
		-w )	shift
				SEED_LEN=$2
				;;
		-h )	usage
				exit
				;;
		* )		usage
				exit 1
	esac
	shift
done

#Simulate queries#
echo "Simulate query sequences"

#Calculate number of batches
NB_BATCHES=$(($NB_QUERIES / $BATCH_SIZE))
#Create query directory
mkdir -p queries

#For each batch...
for ((i=1; i<=$NB_BATCHES; i++))
do
	#Generate random query file
	scripts/GenRandDNA_Seq.py $Q_LEN $BATCH_SIZE > queries/Random_"$i".q
done

#Query the graph#
echo "Query the graph"

#Create result directory
mkdir -p results

#Query each batch
for ((i=1; i<=$NB_BATCHES; i++))
do
	$PLAST_BIN Search -i $(echo $GRAPH_PATH | sed 's/.gfa//g')  -q queries/Random_"$i".q -k $K -g $MINI_LEN -w $SEED_LEN -e 100 -u > results/maxScores_"$i".txt 2> results/maxScores_"$i".log
done

#Aggregate scores#
echo "Aggregate scores"

scripts/aggScores.py results/maxScores_*.txt > results/aggMaxScores.txt

#Estimate parameters#
echo "Estimate parameters"

scripts/estParams.py results/aggMaxScores.txt > results/parameters.txt

#Clean up#
echo "Clean up"

rm -r queries/* results/maxScores_* results/aggMaxScores.txt
