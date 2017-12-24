import os
import re

# this function calculates the allele frequency for each variation in the haplotip files
# param hapDir - a path to a dir holding the haplotip files.
# the function creates new files that holds the allele frequencies of each variation
# in the position of the original variation at the haplotip file
def calculateAlleleFrequency(hapDir):
    for hapFile in os.listdir(hapDir):
        if (hapFile.startswith("tagc128.shapeit.") and hapFile.endswith(".hap")):
            m = re.search('tagc128\.shapeit\.([0-9]+)\.hap', hapFile)
            m = re.search('(tagc128\.shapeit\.[0-9]+)\.hap', hapFile)
            outPath = m.group(1)
            inputFile = open(hapDir + "/" + hapFile, 'r')
            outputFile = open(hapDir + "/" + outPath + ".freq", 'w')
            for line in inputFile:
                alleles = line.split(" ")
                minor = 0
                for i in range(len(alleles) - 1):
                    minor += int(alleles[i])
                minor += int((alleles[-1]).split("\n")[0])
                numOfAl = float(len(alleles))
                outputFile.write(str(minor/numOfAl) + "\n")
            outputFile.close()
            inputFile.close()
    return


