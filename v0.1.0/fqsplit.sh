#!/bin/bash

# Shell script for splitting fastq files
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of last modification: 05/10/2018

for f in ./*.fastq;
do
	filename=$(basename ${f})
	split -l$((`wc -l < ${filename}`/14)) ${filename} ${filename}.split -da 4
done
