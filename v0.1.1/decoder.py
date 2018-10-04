#!/usr/bin/python

# Program to decode ddSeq barcodes
# for custom secondary analysis
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of Last Modification - 08/31/2018

import time
import distance
import sys
import pandas as pd

# abbreviations: bc - barcode, umi - Unique Molecular Identifier, acg - sequence ACG, polyT - GACTTTTT... block

linker1 = "TAGCCATCGCATTGC"
linker2 = "TACCTCTGAGCTGAA"
bc_length = 6
linker_length = 15
acg_length = 3
umi_length = 8
polyT_length = 6

rejected_reads = []
accept_index = []
k = 15
columns = ["ReadNumber","BC1","BC2","BC3","ACG","UMI","GACTTT","IndicesString"]

def getAllBlocks(read, read_number):
    if findLinkerIndices(read, read_number):
        [l1_ind, l2_ind] = findLinkerIndices(read, read_number)
        [bc1_ind,bc2_ind,bc3_ind,acg_ind,umi_ind,polyT_ind] = findAllOtherIndices(l1_ind, l2_ind)
        # all indices are the starting points of the respective blocks
        bc1 = read[bc1_ind : bc1_ind+bc_length]
        bc2 = read[bc2_ind : bc2_ind+bc_length]
        bc3 = read[bc3_ind : bc3_ind+bc_length]
        acg = read[acg_ind : acg_ind+acg_length]
        umi = read[umi_ind : umi_ind+umi_length]
        polyT = read[polyT_ind : polyT_ind+polyT_length]
        indices_string = "-".join(str(ind) for ind in [bc1_ind,bc2_ind,bc3_ind,acg_ind,umi_ind,polyT_ind])
        all_blocks = [read_number,bc1,bc2,bc3,acg,umi,polyT,indices_string]
        return all_blocks
    else:
        return None

def findLinkerIndices(read, read_number):
    l2_ind = 0
    l1_ind = 0
    for i in range(0,len(read)-k):
        kmer = read[i:i+k]
        if isLinker(kmer,linker1):
            l1_ind = i
        elif isLinker(kmer,linker2):
            l2_ind = i
            break
    if ((l2_ind - l1_ind) != 21) or (l1_ind < 6):
        rejectRead(read_number)
        return None
    return [l1_ind, l2_ind]

def findAllOtherIndices(l1_ind, l2_ind):
    bc1_ind = l1_ind - bc_length
    bc2_ind = l2_ind - bc_length
    bc3_ind = l2_ind + linker_length
    acg_ind = bc3_ind + bc_length
    umi_ind = acg_ind + acg_length
    polyT_ind = umi_ind + umi_length
    return [bc1_ind,bc2_ind,bc3_ind,acg_ind,umi_ind,polyT_ind]

def isLinker(kmer,linker):
    if distance.hamming(kmer,linker) <= 1:
        return True
    else:
        return False

def rejectRead(read_number):
    rejected_reads.append(read_number)
    return

def main():
    #read = "NCAGGATTGTAGCCATCGCATTGCGAAGGGTACCTCTGAGCTGAAATTAGTACGCATATAAAGACTTG"
    #read_number = 1
    start_time = time.time()
    rows = []
    inFile = sys.argv[1]
    with open(inFile,'r') as fIn:
        for i,line in enumerate(fIn):
            if i%4 == 1:
                read = line.strip()
                read_number = i/4
                all_blocks = getAllBlocks(read, read_number)
                if all_blocks != None:
                    rows.append(all_blocks)
    df = pd.DataFrame(rows, columns=columns)
    print "Time taken to decode:", time.time() - start_time
    start_time = time.time()
    df.to_csv("out.tsv", sep="\t", index=False, encoding="utf-8")
    print "Time to write:", time.time() - start_time

if __name__ == '__main__':
    main()
