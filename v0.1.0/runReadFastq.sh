#!/bin/bash

# Shell script for running the python script
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of last modification: 05/11/2018

mkdir formattedSeqs
for f in ./*R1.fastq.split*;
do
	file1=$(basename ${f})
	file1Part1=${file1:0:13}
	file1Part2=R2${file1:15:${#file1}}
	file2=$file1Part1$file1Part2
	echo $file1
	echo $file2
	time python readFastq.py ${file1} ${file2}
done
