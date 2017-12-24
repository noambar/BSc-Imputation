import os
import re

# this function reads the reference and writes new files in the format:
# "tagc128.shapeit.[chromosome].as_pair"
# where each reference is represented by the sum of its haplotypes.
# @param hapDir - a path to the reference's directory.
def sumRefPairs(hapDir):
    for hapFile in os.listdir(hapDir):
        if (hapFile.startswith("tagc128.shapeit.") and hapFile.endswith(".hap")):
            m = re.search('(tagc128\.shapeit\.[0-9]+)\.hap', hapFile)
            outPath = m.group(1)
            inputFile = open(hapDir + "/" + hapFile, 'r')
            outputFile = open(hapDir + "/" + outPath + ".as_pair", 'w')
            for line in inputFile:
                # only non empty lines
                if (line.strip() != ""):
                    alleles = (line.strip()).split(" ")
                    for i in range(0, (len(alleles)) - 2, 2):
                        outputFile.write(str(int(alleles[i]) + int(alleles[i + 1])) + " ")
                    outputFile.write(str(int(alleles[-2]) + int(alleles[-1])) + "\n")
            outputFile.close()
            inputFile.close()
    return