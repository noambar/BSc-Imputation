import random
import argparse
import os
import re
import datetime

GENOME_SIZE = 3200000000
READ_LENGTH = 148

# this file is used in order to downsample the trio's alignment files
# samFile - the alignment file to downsample from
# numOfFileRows - the number of rows/reads in the input files
# resquestedCoverage - the requested down sampled coverage
# readLength - the average read length of the input file
def main():

    # print str(datetime.datetime.now())
    parser = argparse.ArgumentParser()
    parser.add_argument('samFile', type=str, help='a path to the sam file', default='no_file_given')
    parser.add_argument('numOfFileRows', type=int, help='the number of the files rows', default=-1)
    parser.add_argument('requestedCoverage', type=float, help='the requested genome coverage', default=1)
    parser.add_argument('readLength', type=int, help='the average read length', default=-1)

    args = parser.parse_args()

    if (os.path.isfile(args.samFile) is False):
        exit(args.samFile + " does not exists")

    SamFile = open(args.samFile, 'r')

    numberOfLinesNeeded = int((GENOME_SIZE * args.requestedCoverage)/float(args.readLength))
    print "number of lines needed: " + str(numberOfLinesNeeded)

    if (numberOfLinesNeeded > args.numOfFileRows):
        exit("requested coverage is higher than the given sam file!")

    print "Generating random list of integers"
    linesNumbers = random.sample(range(args.numOfFileRows), numberOfLinesNeeded)

    print "Sorting list of integers, may take a while..."
    linesNumbers.sort()
    print str(len(linesNumbers))

    counter = 0
    posInFile = 0

    m = re.search('/?(.+).sam', args.samFile)
    outputFileName =  m.group(1) + "_down_sampled_" + str(args.requestedCoverage) + "x_" + str(datetime.datetime.now()).split(" ")[1] + ".sam"
    outputFile = open(outputFileName, 'w')

    print "writing to " + outputFileName + ", this will take a while..."

    for line in SamFile:
        if (counter == linesNumbers[posInFile]):
            outputFile.write(line)
            posInFile += 1
        if (posInFile == numberOfLinesNeeded):
            break
        counter += 1

    outputFile.close()
    SamFile.close()

    SamFile.close()

    # os.system("""sbatch -c 10 --mem-per-cpu=6000 --time=3-0 --wrap="samtools view -s """ + str(factor) + """ """ + args.samFile + """ > """ + outputFileName + """ " &""")

    print "succsess"

    return

main()
