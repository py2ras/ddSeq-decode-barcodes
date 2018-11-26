#!/usr/bin/python

# Program to correct the filtered reads
# [TODO]
# 1. Add more detailed description
# 2. bktree is still in the test stage - the modified version would be added soon
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of Last Modification - 09/06/2018

import sys
import pandas as pd
import time
import ddBKtree

columns = ["ReadNumber","BC1","BC2","BC3","ACG","UMI","GACTTT","IndicesString"]
possible_bar_codes = ['AAAGAA', 'AACAGC', 'AACGTG', 'AAGCCA', 'AAGTAT', 'AATTGG', 'ACAAGG', 'ACCCAA', 'ACCTTC', 'ACGGAC', 'ACTGCA', 'AGACCC', 'AGATGT', 'AGCACG', 'AGGTTA', 'AGTAAA', 'AGTCTG', 'ATACTT', 'ATAGCG', 'ATATAC', 'ATCCGG', 'ATGAAG', 'ATTAGT', 'CAACCG', 'CAAGTC', 'CACCAC', 'CACTGT', 'CAGACT', 'CAGGAG', 'CATAGA', 'CCACGC', 'CCGATG', 'CCGTAA', 'CCTCTA', 'CGAAAG', 'CGAGCA', 'CGCATA', 'CGGCGT', 'CGGTCC', 'CGTTAT', 'CTAGGT', 'CTATTA', 'CTCAAT', 'CTGTGG', 'CTTACG', 'CTTGAA', 'GAAATA', 'GAAGGG', 'GACTCG', 'GAGCTT', 'GAGGCC', 'GAGTGA', 'GATCAA', 'GCCAGA', 'GCCGTT', 'GCGAAT', 'GCGCGG', 'GCTCCC', 'GCTGAG', 'GCTTGT', 'GGACGA', 'GGATTG', 'GGCCAT', 'GGGATC', 'GGTAGG', 'GGTGCT', 'GTACAG', 'GTCCTA', 'GTCGGC', 'GTGGTG', 'GTTAAC', 'GTTTCA', 'TAAGCT', 'TAATAG', 'TACCGA', 'TAGAGG', 'TATTTC', 'TCAGTG', 'TCATCA', 'TCCAAG', 'TCGCCT', 'TCGGGA', 'TCTAGC', 'TGAATT', 'TGAGAC', 'TGCGGT', 'TGCTAA', 'TGGCAG', 'TGTGTA', 'TGTTCG', 'TTAAGA', 'TTCGCA', 'TTCTTG', 'TTGCTC', 'TTGGAT', 'TTTGGG']

def simpleCorrection(filteredDf):
    # exhaustive search algorithm - very slow
    rows = filteredDf.values.tolist()
    corrected_rows = []
    corrected_reads = 0
    for row in rows:
        bc1,bc2,bc3 = row[1],row[2],row[3]
        old_row = [x for x in row]
        for possible_bar_code in possible_bar_codes:
            bc1_ed = distance.hamming(bc1,possible_bar_code)
            bc2_ed = distance.hamming(bc2,possible_bar_code)
            bc3_ed = distance.hamming(bc3,possible_bar_code)
            if bc1_ed == 1:
                row[1] = possible_bar_code
            if bc2_ed == 1:
                row[2] = possible_bar_code
            if bc3_ed == 1:
                row[3] = possible_bar_code
        if old_row != row:
            corrected_reads = corrected_reads + 1
    print "No. of reads corrected:",corrected_reads
    return rows

def hamming_distance(s1,s2):
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def BKTreeCorrection(filteredDf):
    # implements BK Tree for searching similar words - very fast
    #tree = ddBKtree.Tree(possible_bar_codes[0])
    #for possible_bar_code in possible_bar_codes[1:]:
    #    tree.add_child(possible_bar_code)
    tree = ddBKtree.Tree()
    tree.add_word_list(possible_bar_codes)
    rows = filteredDf.values.tolist()
    corrected_reads = 0
    for row in rows:
        read_number = row[0]
        bc1,bc2,bc3 = row[1],row[2],row[3]
        old_row = [x for x in row]
        row[1] = getNewBC(bc1,tree)
        row[2] = getNewBC(bc2,tree)
        row[3] = getNewBC(bc3,tree)
        if old_row != row:
            corrected_reads = corrected_reads + 1
    print "No. of reads corrected:",corrected_reads
    return rows

def getNewBC(bc, tree):
    bc_new = tree.search(bc,1)
    if len(bc_new) == 0:
        return bc
    else:
        return bc_new[0][1]

def main():
    start_time = time.time()
    if len(sys.argv) < 2:
        print "Usage: python " + sys.argv[0] + " filtered_reads_table"
        quit()
    inFile = sys.argv[1]
    filteredDf = pd.read_csv(inFile, sep='\t', header=0)
    corrected_rows = BKTreeCorrection(filteredDf)
    corrected_df = pd.DataFrame(corrected_rows, columns=columns)
    corrected_df.to_csv("cout.tsv", sep='\t',index=False, encoding="utf-8")
    print "Time taken:", time.time() - start_time

if __name__ == '__main__':
    main()

