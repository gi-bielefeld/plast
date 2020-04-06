#!/usr/bin/env bash

#Parse input parameters
#	1. blat executable (has to be an absolute path)
#	1. database list
#	2. query sequences file
#	3. output file name
#	4. runtime file

BLAT_BIN=$1
DB_LIST=../"$2"
QUERY=../"$3"
RES_FILE=../"$4"
RT_FILE=../"$5"

MAX_DB_NB=750

#Check if run has to be splitted
NB_DBS=$(wc -l $DB_LIST | cut -d' ' -f1)

if [[ $NB_DBS -lt $MAX_DB_NB ]]; then
	#Simply run BLAT once
	/usr/bin/time -v $BLAT_BIN $DB_LIST $QUERY $RES_FILE 2> $RT_FILE
else
	#Split database list into managable chuncks and perform runs
	cat $DB_LIST > "$DB_LIST".tmp
	PART=1

	while
		head -n $MAX_DB_NB "$DB_LIST".tmp > "$DB_LIST".part"$PART"
		tail -n +$(($MAX_DB_NB+1)) "$DB_LIST".tmp > "$DB_LIST".tmp2
		mv "$DB_LIST".tmp2 "$DB_LIST".tmp
		/usr/bin/time -v $BLAT_BIN "$DB_LIST".part"$PART" $QUERY "$RES_FILE".part"$PART" 2> "$RT_FILE".part"$PART"
		PART=$(($PART+1))
		NB_DBS=$(wc -l "$DB_LIST".tmp | cut -d' ' -f1)
		[[ $NB_DBS != 0 ]]
	do true; done

	#Merge results
	cat "$RES_FILE".part* > $RES_FILE
	RT=0
	MAX_MEM=0

	for i in "$RT_FILE".part*; do
		#Time...
		TM_LINE=$(grep User $i)
		RT=$(echo $RT + $(echo $TM_LINE | cut -d' ' -f4) | bc)
		#...and memory
		MEM_LINE=$(grep Max $i)
		CUR_MEM=$(echo $MEM_LINE | cut -d' ' -f6)

		if [[ $MAX_MEM -lt $CUR_MEM ]]; then
			MAX_MEM=$CUR_MEM
		fi
	done
	echo $(echo $TM_LINE | cut -d' ' -f-3)" "$(echo $RT) > $RT_FILE
	echo $(echo $MEM_LINE | cut -d' ' -f-5)" "$(echo $MAX_MEM) >> $RT_FILE

	#Clean up
	rm "$DB_LIST".part* "$RES_FILE".part* "$DB_LIST".tmp
fi