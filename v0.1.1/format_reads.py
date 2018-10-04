#!/usr/bin/python

# Program to format the read 2 files
# in correspondence with the read 1 files
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of Last Modification - 09/17/2018

import pandas as pd
import sys
import time

def formatRead1(reads1_file, corrected_reads1_df):
    # [TODO]
    # Divide this into further smaller functions
    # 1. i == ((accepted_read_number+1)*4 can be used in another function
    # 2. read_count < len(accepted_read_numbers)-1 can be used in another function
    concat_bc_df = corrected_reads1_df.iloc[:,1:4].apply(lambda x: ''.join(x), axis=1)
    accepted_read_numbers = corrected_reads1_df["ReadNumber"].tolist()
    output = []
    with open(reads1_file) as r1:
        i = 0
        read_count = 0
        accepted_read_number = accepted_read_numbers[read_count]
        for line in r1:
            if isValidRead(i,accepted_read_number):
                print "read Count: %s; i: %s" % (accepted_read_number, i)
                output.append(line)
                i += 1
            else:
                i += 1
            if i == ((accepted_read_number+1)*4) and read_count < len(accepted_read_numbers)-1:
                read_count += 1
                accepted_read_number = accepted_read_numbers[read_count]
    return output

def isValidRead(i,accepted_read_number):
    # the four lines corresponding to the 
    # accepted read number in the file
    # should be accepted
    # print "i: %s; read: %s; range: %s" % (i, accepted_read_number, range((accepted_read_number*4), (accepted_read_number+1)*4))
    if i in range((accepted_read_number*4), (accepted_read_number+1)*4):
        return True
    else:
        return False

def filterRead2():
    pass

def writeFormattedRead1(formatted_reads):
    with open("formatted_read1.txt","w") as outFile:
        for line in formatted_reads:
            outFile.write(line)

def main():
    start_time = time.time()
    corrected_reads1_file = sys.argv[1]
    reads1_file = sys.argv[2]
    reads2_file = sys.argv[3]
    corrected_reads1_df = pd.read_csv(corrected_reads1_file, sep='\t', header=0)
    formatted_reads = formatRead1(reads1_file, corrected_reads1_df)
    writeFormattedRead1(formatted_reads)
    print "Total time:", time.time() - start_time
    #filterRead2(reads2_file, accepted_read_numbers)

    pass

if __name__ == '__main__':
    main()
