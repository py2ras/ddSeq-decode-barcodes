#!/usr/bin/python

# Program to filter the reads based
# the ACG and PolyT sequences and positions
# @author - Sarthak Sharma
# Date of last modification - 09/06/2018

import pandas as pd
import sys
import distance
import time

def removeBadBlocks(df):
    # to remove ACG and GACTTT (polyT) blocks allowing 1ED
    rows = df.values.tolist()
    columns = ["ReadNumber","BC1","BC2","BC3","ACG","UMI","GACTTT","IndicesString"]
    filteredRows = []
    for row in rows:
        acg = row[4]
        #polyT = str(row[6])[0:6].ljust(6,'T')
        polyT = str(row[6])
        if isBlockCorrect(acg, "ACG") and isBlockCorrect(polyT, "GACTTT"):
            filteredRows.append(row)
        else:
            continue
    filteredDataFrame = pd.DataFrame(filteredRows, columns = columns)
    return filteredDataFrame

def isBlockCorrect(block,correctSequence):
    try:
        if len(block) == len(correctSequence):
            if distance.hamming(block,correctSequence) <= 1:
                return True
        else:
            return False
    except TypeError:
        print block
        return False

def main():
    inFile = sys.argv[1]
    start_time = time.time()
    df = pd.read_csv(inFile, sep='\t', header=0)
    filteredDf = removeBadBlocks(df)
    filteredDf.to_csv("fout.tsv", sep="\t", index=False, encoding="utf-8")
    print "Time to filter reads:", time.time() - start_time

if __name__ == '__main__':
    main()
