import re
import os

# this function gets the directory of the proband output files, and parse them in order to extract
# only the coordinates that had some evidence
# param probandDir - a path to the directory of proband files (as outputed by createVarMatrix())
# returns a list of lists of tuples. where each tuple is (origPos, major, minor)
def extractProbandVariations(probandDir):
    maxChr = 0
    # get number of chromosome files
    for chrFile in os.listdir(probandDir):
        if (chrFile.startswith("chr") and chrFile.endswith(".out")):
            m = re.search('chr([0-9]+)\.out', chrFile)
            chromosome =  m.group(1)
            if (int(chromosome) > maxChr):
                maxChr = int(chromosome)
    for chrFile in os.listdir(probandDir):
        # print chrFile ## todo delete
        if (chrFile.startswith("chr") and chrFile.endswith(".out")):
            m = re.search('chr([0-9]+)\.out', chrFile)
            chromosome =  int(m.group(1))
            inputFile = open(probandDir + "/" + chrFile, 'r')
            outputFile = open(probandDir + "/" + "chr" + str(chromosome) + "_probandMat.txt", 'w')
            counter = 0

            # adding only rows that had a non 0 value
            for line in inputFile:
                major = line.split(" ")[0]
                minor = line.split(" ")[1]
                position = (line.split(" ")[2]).split("\n")[0]
                if (int(minor) > 0 or int(major) > 0):
                    outputFile.write(str(counter) + " " + position + " " + major + " " + minor + "\n")
                counter += 1
            outputFile.close()
            inputFile.close()

    return