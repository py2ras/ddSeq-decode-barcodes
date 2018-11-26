#!/usr/bin/python

# Program to format the read 2 files
# in correspondence with the read 1 files
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of Last Modification - 09/17/2018

import pandas as pd
import sys
import time
from itertools import izip
import os

def formatReads(reads1_file, reads2_file, corrected_reads1_df):
    '''
    [TODO]
    Divide this into further smaller functions
    1. i == ((accepted_read_number+1)*4 can be used in another function
    2. read_count < len(accepted_read_numbers)-1 can be used in another function
    '''

    print "Formatting Reads ..."

    start_time = time.time()

    # if one would want to create the new files
    # in the same directory as the original
    #r1_dirname = os.path.dirname(reads1_file)
    #r2_dirname = os.path.dirname(reads2_file)

    r1_filename = os.path.basename(reads1_file)
    r2_filename = os.path.basename(reads2_file)

    # table with the concatenated reads
    concat_bc_df = corrected_reads1_df.iloc[:,1:7].apply(lambda x: ''.join(x), axis=1)
    accepted_read_numbers = corrected_reads1_df["ReadNumber"].tolist()

    r1_in = open(reads1_file,"r")
    r2_in = open(reads2_file,"r")
    r1_out = open("formatted_" + r1_filename,"w")
    r2_out = open("formatted_" + r2_filename,"w")

    #r1_out = open(r1_dirname + "/formatted_" + reads1_file,"w")
    #r2_out = open(r2_dirname + "/formatted_" + reads2_file,"w")

    line_num = 0
    read_count = 0
    accepted_read_number = accepted_read_numbers[read_count]
    for line_r1,line_r2 in izip(r1_in,r2_in):
        # check if the read is a valid read 
        if isValidRead(line_num,accepted_read_number):
            #print "read Count: %s; i: %s" % (accepted_read_number, i)

            # write the corresponding read to the new read2 file
            r2_out.write(line_r2)

            # check if the line is a sequence
            if isSequence(line_num):
                # write the formatted sequence instead of the original sequence
                r1_out.write(concat_bc_df[read_count]+'\n')
            else:
                # otherwise simply copy the original file contents
                r1_out.write(line_r1)
            line_num += 1
        else:
            line_num += 1
        # check if the line number corresponds to an accepted read
        if line_num == ((accepted_read_number+1)*4) and read_count < len(accepted_read_numbers)-1:
            read_count += 1
            accepted_read_number = accepted_read_numbers[read_count]

    r1_in.close()
    r2_in.close()
    r1_out.close()
    r2_out.close()

    print "Reads formatted!"

    print "Time:", time.time() - start_time


def isSequence(line_num):
    if line_num%4 == 1:
        return True
    else:
        return False


def formatRead1_test(reads2_df, corrected_reads1_df):
    '''
    This is a test function which uses pandas tables for
    formatting reads 1 file.
    But my tests show that this method is really slow.
    This could be because I am trying to do a lot of 
    formatting after loading the files as tables.
    This method could be faster for filtering reads 2
    since that requires only merging.
    '''
    reads1_df["LineNumber"] = reads1_df.index
    reads1_df["ReadNumber"] = reads1_df["LineNumber"].apply(lambda x: int(x/4))
    comb_df = pd.merge(reads1_df, corrected_reads1_df, on="ReadNumber")
    comb_df["ConcatRead"] = comb_df[comb_df.columns[3:9]].apply(lambda x: ''.join(x), axis=1)
    comb_df = comb_df.drop(comb_df.columns[[1,3,4,5,6,7,8,9]], axis=1)
    comb_df.to_csv("test_pd_read1.fastq",sep="\t",index=False)
    print "Done"


def isValidRead(i,accepted_read_number):
    # the four lines corresponding to the 
    # accepted read number in the file
    # should be accepted
    # print "i: %s; read: %s; range: %s" % (i, accepted_read_number, range((accepted_read_number*4), (accepted_read_number+1)*4))
    if i in range((accepted_read_number*4), (accepted_read_number+1)*4):
        return True
    else:
        return False


def filterRead2(reads2_file, corrected_reads1_df):
    print "Formatting Reads 2 ... "

    start_time = time.time()
    
    reads2_df = pd.read_csv(reads2_file,sep="\n",names=["Lines"])
    reads2_df["LineNumber"] = reads2_df.index
    reads2_df["ReadNumber"] = reads2_df["LineNumber"].apply(lambda x: int(x/4))
    comb_df = pd.merge(reads2_df, corrected_reads1_df, on="ReadNumber")
    comb_df.to_csv("test_pd_read2.fastq",sep="\t",index=False, header=None)

    print "Reads 2 formatted!"

    print "Time:", time.time() - start_time


def writeFormattedRead1(formatted_reads):
    reads2_df = pd.read_csv(reads2_file,sep="\n",names=["Lines"])
    with open("formatted_read1.txt","w") as outFile:
        for line in formatted_reads:
            outFile.write(line)


def main():
    start_time = time.time()
    if len(sys.argv) < 4:
        print "Usage: python " + sys.argv[0] + " corrected_reads_1 reads_1 reads_2"
        quit()
    corrected_reads1_file = sys.argv[1]
    reads1_file = sys.argv[2]
    reads2_file = sys.argv[3]
    corrected_reads1_df = pd.read_csv(corrected_reads1_file, sep='\t', header=0)
    formatReads(reads1_file, reads2_file, corrected_reads1_df)

    print "Total time:", time.time() - start_time

if __name__ == '__main__':
    main()
