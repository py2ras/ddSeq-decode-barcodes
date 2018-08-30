#!/usr/bin/python

# Program to read FASTQ files from single cell sequencing
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of Last Modification - 05/07/2018

import sys
import distance
#from multiprocessing.dummy import Pool as ThreadPool
#import threading

class FastqReader:
    def __init__(self):
        self.stringPairs = []
        self.sequences = []
        self.readWiseKmers = []
        self.linker1 = "TAGCCATCGCATTGC"
        self.linker2 = "TACCTCTGAGCTGAA"
        self.linkerIndices = []
        self.filterIndices = []
        self.blocks = []
        self.rejectedReads = []
        self.qualityScores = []
        self.possibleBarCodes = ["AAAGAA","AGTCTG","CCGTAA","GACTCG","GGTAGG","TCGCCT","AACAGC","ATACTT","CCTCTA","GAGCTT","GGTGCT","TCGGGA","AACGTG","ATAGCG","CGAAAG","GAGGCC","GTACAG","TCTAGC","AAGCCA","ATATAC","CGAGCA","GAGTGA","GTCCTA","TGAATT","AAGTAT","ATCCGG","CGCATA","GATCAA","GTCGGC","TGAGAC","AATTGG","ATGAAG","CGGCGT","GCCAGA","GTGGTG","TGCGGT","ACAAGG","ATTAGT","CGGTCC","GCCGTT","GTTAAC","TGCTAA","ACCCAA","CAACCG","CGTTAT","GCGAAT","GTTTCA","TGGCAG","ACCTTC","CAAGTC","CTAGGT","GCGCGG","TAAGCT","TGTGTA","ACGGAC","CACCAC","CTATTA","GCTCCC","TAATAG","TGTTCG","ACTGCA","CACTGT","CTCAAT","GCTGAG","TACCGA","TTAAGA","AGACCC","CAGACT","CTGTGG","GCTTGT","TAGAGG","TTCGCA","AGATGT","CAGGAG","CTTACG","GGACGA","TATTTC","TTCTTG","AGCACG","CATAGA","CTTGAA","GGATTG","TCAGTG","TTGCTC","AGGTTA","CCACGC","GAAATA","GGCCAT","TCATCA","TTGGAT","AGTAAA","CCGATG","GAAGGG","GGGATC","TCCAAG","TTTGGG"]
        #self.numCorrectedBarCodes = 0
        self.sequenceNames = []
        self.sequencePlaceHolders = []
        self.blockSequenceNames = []
        self.blockPlaceHolders = []
        self.blockQualityScores = []
        self.read2Lines = []

    def getSequences(self, inFile1,inFile2):
        lines = []
        print "Reading Sequences ..."
        with open(inFile1) as fIn:
            lines = fIn.readlines()
            #self.sequences = fIn.readlines()
        self.sequenceNames = [read.strip() for read in lines[0::4]]
        self.sequences = [read.strip() for read in lines[1::4]]
        self.sequencePlaceHolders = [read.strip() for read in lines[2::4]]
        self.qualityScores = [read.strip() for read in lines[3::4]]
        with open(inFile2) as fIn2:
            lines = fIn2.readlines()
        self.R2sequenceNames = [read.strip() for read in lines[0::4]]
        self.R2sequences = [read.strip() for read in lines[1::4]]
        self.R2sequencePlaceHolders = [read.strip() for read in lines[2::4]]
        self.R2qualityScores = [read.strip() for read in lines[3::4]]
        del lines

    def printSequences(self):
        for sequence in self.sequences:
            print sequence

    def printReadKmers(self):
        for readKmers in self.readWiseKmers:
            print readKmers

    def findLinkerIndices(self):
        print "Finding Linker Indices ..."
        for i,readKmers in enumerate(self.readWiseKmers):
            linker1Index, linker2Index = self._getLinkerIndices(readKmers)
            self.linkerIndices.append((linker1Index, linker2Index))

    def _getLinkerIndices(self, readKmers):
        minLinker1Ed = 2
        minLinker2Ed = 2
        linker1Index = None
        linker2Index = None
        for i in range(len(readKmers)):
            kmer = readKmers[i]
            linker1Ed = distance.hamming(kmer,self.linker1)
            if linker1Ed < minLinker1Ed:
                minLinker1Ed = linker1Ed
                linker1Index = i
                break
        if linker1Index != None:
            linker2Start = linker1Index + 15
        else:
            linker2Start = 0

        for i in range(linker2Start,len(readKmers)):
            kmer = readKmers[i]
            linker2Ed = distance.hamming(kmer,self.linker2)
            if linker2Ed < minLinker2Ed:
                minLinker2Ed = linker2Ed
                linker2Index = i
                break
        return linker1Index, linker2Index

    def checkSequenceLength(self):
        lengths = {}
        for sequence in self.sequences:
            if len(sequence) in lengths:
                lengths[len(sequence)] += 1
            else:
                lengths[len(sequence)] = 0
        print lengths

    def findOtherBlocks(self):
        print "Finding other blocks ..."
        for i in range(len(self.linkerIndices)):
            indexPair = self.linkerIndices[i]
            if not self.isInvalidIndexPair(indexPair):
                linker1Index = indexPair[0]
                linker2Index = indexPair[1]
                sequenceName = self.sequenceNames[i]
                sequence = self.sequences[i]
                placeHolder = self.sequencePlaceHolders[i]
                qualityScore = self.qualityScores[i]
                R2sequenceName = self.R2sequenceNames[i]
                R2sequence = self.R2sequences[i]
                R2placeHolder = self.R2sequencePlaceHolders[i]
                R2qualityScore = self.R2qualityScores[i]

                barCode1 = sequence[linker1Index-6:linker1Index]
                barCode2 = sequence[linker1Index+15:linker2Index]
                barCode3 = sequence[linker2Index+15:linker2Index+15+6]
                acgBlock = sequence[linker2Index+15+6:linker2Index+15+6+3]
                umiBlock = sequence[linker2Index+15+6+3:linker2Index+15+6+3+8]
                polyTBlock = sequence[linker2Index+15+6+3+8:linker2Index+15+6+3+8+4]

                barCode1QS = qualityScore[linker1Index-6:linker1Index]
                barCode2QS = qualityScore[linker1Index+15:linker2Index]
                barCode3QS = qualityScore[linker2Index+15:linker2Index+15+6]
                umiBlockQS = qualityScore[linker2Index+15+6+3:linker2Index+15+6+3+8]
                if self.barCodesAreCorrect(barCode1,barCode2,barCode3):
                    for possibleBarCode in self.possibleBarCodes:
                        bc1Ed = distance.hamming(barCode1,possibleBarCode)
                        bc2Ed = distance.hamming(barCode2,possibleBarCode)
                        bc3Ed = distance.hamming(barCode3,possibleBarCode)
                        if bc1Ed == 1:
                            barCode1 = possibleBarCode
                            #self.numCorrectedBarCodes = self.numCorrectedBarCodes + 1
                        if bc2Ed == 1:
                            barCode2 = possibleBarCode
                            #self.numCorrectedBarCodes = self.numCorrectedBarCodes + 1
                        if bc3Ed == 1:
                            barCode3 = possibleBarCode
                            #self.numCorrectedBarCodes = self.numCorrectedBarCodes + 1
                if self.blocksAreCorrect(acgBlock,polyTBlock):
                    self.blockSequenceNames.append(sequenceName)
                    self.blocks.append([barCode1,barCode2,barCode3,acgBlock,umiBlock,polyTBlock])
                    self.blockPlaceHolders.append(placeHolder)
                    self.blockQualityScores.append([barCode1QS,barCode2QS,barCode3QS,umiBlockQS])
                    self.read2Lines.extend([R2sequenceName,R2sequence,R2placeHolder,R2qualityScore])
                else:
                    self.filterIndices.append(i)
        #print "Number of corrected Barcodes: " + str(self.numCorrectedBarCodes)

    def barCodesAreCorrect(self,bc1,bc2,bc3):
        if (len(bc1) != 6) or (len(bc2) != 6) or (len(bc3) != 6):
            return False
        else:
            return True

    def blocksAreCorrect(self,block1,block2):
        if (len(block2) < 4):
            return False
        block1Ed = distance.hamming(block1,"ACG")
        block2Ed = distance.hamming(block2,"GACT")
        if (block1Ed > 1) or (block2Ed > 1):
            return False
        else:
            return True

    def correctBarcodes(self):
        print "Correcting barcodes ..."
        for i in range(len(self.blocks)):
            bc1 = self.blocks[i][0]
            bc2 = self.blocks[i][1]
            bc3 = self.blocks[i][2]
            for possibleBarCode in self.possibleBarCodes:
                bc1Ed = distance.hamming(bc1,possibleBarCode)
                bc2Ed = distance.hamming(bc2,possibleBarCode)
                bc3Ed = distance.hamming(bc3,possibleBarCode)
                if bc1Ed == 1:
                    self.blocks[i][0] = possibleBarCode
                if bc2Ed == 1:
                    self.blocks[i][1] = possibleBarCode
                if bc3Ed == 1:
                    self.blocks[i][2] = possibleBarCode

    def isInvalidIndexPair(self, indexPair):
        if (indexPair[1] > 50) or (indexPair[0] < 6) or (indexPair[0] == None) or (indexPair[1] == None) or (indexPair[1] - indexPair[0] != 21):
            return True
        else:
            return False

    def filterReads(self):
        print "Filering Reads ..."
        preFilterLength = len(self.sequences)
        for i, indexPair in enumerate(self.linkerIndices):
            if self.isInvalidIndexPair(indexPair):
                self.filterIndices.append(i)
        #for i in self.filterIndices:
        #    self.rejectedReads.append(self.sequences[i])
            #self.sequences[i] = None
        postFilterLength = len(self.sequences) - len(self.filterIndices)
        print "Started with " + str(preFilterLength) + " reads. Now remaining " + str(postFilterLength) + " reads"


    def getKmers(self):
        print "Getting k-mers ..."
        # pool = ThreadPool(4)
        # self.readWiseKmers = pool.map(self._buildKmers, self.sequences)
        # pool.close()
        # pool.join()
        for sequence in self.sequences:
            readKmers = self._buildKmers(sequence,15)
            self.readWiseKmers.append(readKmers)

    def _buildKmers(self,sequence,k=15):
        kmers = []
        for i in range(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmers.append(kmer)
        return kmers

    def printLinkerIndices(self):
        for linkerIndex in self.linkerIndices:
            print linkerIndex

    def printBlocks(self):
        for block in self.blocks:
            print block

    def writeFormattedSequences(self,outFile1,outFile2):
        print "Writing to Output files ..."
        with open("./formattedSeqs/"+outFile1+".formatted",'w') as fOut:
            for i in range(len(self.blocks)):
                sequenceName = self.blockSequenceNames[i]
                block = self.blocks[i]
                placeHolder = self.blockPlaceHolders[i]
                blockQualityScore = self.blockQualityScores[i]
                formattedSeq = block[0] + block[1] + block[2] + block[4]
                formattedQS = ''.join(blockQualityScore)
                fOut.write(sequenceName)
                fOut.write("\n")
                fOut.write(formattedSeq)
                fOut.write("\n")
                fOut.write(placeHolder)
                fOut.write("\n")
                fOut.write(formattedQS)
                fOut.write("\n")
        with open("./formattedSeqs/"+outFile2+".formatted","w") as fOut2:
            for line in self.read2Lines:
                fOut2.write(line+"\n")
        #with open(outFile+".filterIndices.log",'w') as log:
        #    for i in self.filterIndices:
        #        log.write(str(i) + "\n")
    def checkLengths(self):
        print len(self.sequences)*4
        print len(self.filterIndices)*4
        print len(self.read2Lines)

def main():
    if len(sys.argv)<3:
        print "Usage: python", sys.argv[0], "<input fastq>"
        quit()
    inputFastq1 = sys.argv[1]
    inputFastq2 = sys.argv[2]
    outFile1 = inputFastq1
    outFile2 = inputFastq2
    fr = FastqReader()
    fr.getSequences(inputFastq1,inputFastq2)
    #fr.checkSequenceLength()
    fr.getKmers()
    fr.findLinkerIndices()
    fr.findOtherBlocks()
    #fr.correctBarcodes()
    fr.filterReads()
    fr.writeFormattedSequences(outFile1,outFile2)
    fr.checkLengths()

if __name__ == '__main__':
    main()
