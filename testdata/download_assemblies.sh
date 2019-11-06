#!/usr/bin/env bash

BARCODES="$1"
TOKEN="$2"
mkdir -p fa

tok=$(cat $TOKEN)

while read p;
do
curl -L --user "${tok}:" "https://enterobase.warwick.ac.uk/upload/download?assembly_barcode=${p}&database=senterica" > fa/${p}.fasta;
done < ${BARCODES}
