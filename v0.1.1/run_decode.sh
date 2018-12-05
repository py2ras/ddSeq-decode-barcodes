#!/bin/bash

# Script to run all the decoding functions
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of last modification - 11/26/2018

POSITIONAL=()

while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-r1|--reads1)
			reads1_file="$2"
			shift
			shift
			;;
		-r2|--reads2)
			reads2_file="$2"
			shift
			shift
			;;
		*)
			POSITIONAL+=("$1")
			shift
			;;
	esac
done

set -- "${POSITIONAL[@]}"

./decoder.py $reads1_file
./filterReads.py out.tsv
./correctReads.py fout.tsv
./formatReads.py cout.tsv $reads1_file $reads2_file


