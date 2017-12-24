# import numpy as np
# import sys
import os
import argparse
import re
# import commands

#todo search for genome sizes and comlete
genomeSizeDic = {'hg19' : 3200000000, 'hg18' : 3200000000, 'hg38' : 3200000000, 'mm8' : '1', 'mm9' : '1', 'rn6' : '1'}


# this file is used in order to downsample the trio's alignment files
# it gets a sam file input to downsample from, the requested coverage, and the genome version (for effective size)
def main():
    print "This is the begining!"

    parser = argparse.ArgumentParser()
    parser.add_argument('bamFile', type=str, help='a path to the bam file', default='no_file_given')
    parser.add_argument('requestedCoverage', type=float, help='a float that is the requested coverage', default=1)
    parser.add_argument('genome',type=str,help='the genome version')

    args = parser.parse_args()
    if (args.genome in genomeSizeDic):
        genomeSize = genomeSizeDic[args.genome]
    else:
        exit("Genome version " + args.genome + " is unknown to the program!")

    if (os.path.isfile(args.samFile) is False):
        exit(args.samFile + " does not exists")

    print "Counting lines of file: " + args.samFile + ".\nMay take some time..."
    numOfSamLines = sum(1 for line in open(args.samFile))

    SamFile = open(args.samFile, 'r')

    totalLen = 0
    if (numOfSamLines > 100000):
        for i in range(100000):
            line = SamFile.readline()
            seq = line.split("\t")[9]
            totalLen += len(seq)
        avgReadLen = float(totalLen)/100000
    else:
        for i in range(numOfSamLines):
            line = SamFile.readline()
            seq = line.split("\t")[9]
            totalLen += len(seq)
        avgReadLen = float(totalLen)/numOfSamLines
    SamFile.close()

    estimatedCoverage = (avgReadLen * numOfSamLines)/(float(genomeSize))

    # factor = estimatedCoverage/args.requestedCoverage

    factor = (args.requestedCoverage/estimatedCoverage) + 1

    # numOfLinesNeeded = int(numOfSamLines/factor)
    # roundFactor = int(factor)

    # SamFile = open(args.samFile, 'r')
    m = re.search('/?(.+).sam', args.samFile)
    outputFileName =  m.group(1) + "_" + str(args.requestedCoverage) + ".sam"

    # os.system("samtools view -s " + str(factor) + " > " + outputFileName)
    os.system("""sbatch -c 10 --mem-per-cpu=6000 --time=3-0 --wrap="samtools view -s """ + str(factor) + """ """ + args.samFile + """ > """ + outputFileName + """ " &""")

    # newSamFile = open(args.samFile + "_" + str(args.requestedCoverage), 'w')

    # fivePer = int(float(numOfLinesNeeded)/20)

    # progress = 5

    # for i in range(numOfLinesNeeded):
    #     if ((i+1)%fivePer == 0):
    #         print (str(progress) + "% is done")
    #         progress += 5
    #     for j in range(roundFactor):
    #         line = SamFile.readline()
    #         if (j == 1):
    #             newSamFile.write(line)


    # newSamFile.close()
    # SamFile.close()

    return

main()
