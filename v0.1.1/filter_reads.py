#!/usr/bin/python

# Program to remove reads with incorrect positions
# of ACG or polyT blocks
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of last modification - 09/04/2018

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
    # Sometimes I get a "nan" at the ACG and polyT blocks 
    # when the linkers are too far from the starting position 
    # of the read. In such a case, finding hamming distance
    # throws a TypeError.
    # Nonetheless, such reads should be filtered and hence,
    # I am returning such a block as an incorrect block
    try:
        if len(block) == len(correctSequence):
            if distance.hamming(block,correctSequence) <= 1:
                return True
        else:
            return False
    except TypeError:
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
