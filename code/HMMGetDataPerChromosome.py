import numpy as np
from GetBestTwoRefs import *

START_IDX = 0
END_IDX = 1
NUMBER_OF_VARIANTS_IDX = 2
REF_HOMO = 0
HETERO = 1
ALT_HOMO = 2
PROBAND_LINE_IDX = 0
PROBAND_POS_IDX = 1
REF_COUNT_IDX = 2
ALT_COUNT_IDX = 3
HOMO_THRESHOLD = 10


# this helper function divides the variants of a chromosome to the windows.
# and for each window returns the best two haplotypes.
# the variants are returned as a numpy array with values 0,1,2 for:
# 0 - variant is homozygus to the reference allele.
# 1 - variant is heterozygus.
# 2 - variant is homozygus to the alternative allele.
def getWindowsData(chromosome, probandFile, naivScoresFile, refs):
    refFullList = []
    reFile = open("impute2/tagc128.shapeit." + str(chromosome) + ".hap", 'r')
    # gather all haplotypes
    for line in reFile:
        refFullList.append((line.strip()).split(" "))
    reFile.close()

    windowStats = []
    scoreFile = open(naivScoresFile, 'r')
    # get all stats about each window
    for line in scoreFile:
        windowStats.append(line.split(" ")[:3] + [0, 0])
    scoreFile.close()

    window = 0
    # the data structures that will hold the variants and ref's haplotypes
    refList = np.zeros((len(refFullList), 4))
    variants = np.zeros((len(refFullList), 3))


    counter = 0
    pFile = open(probandFile, 'r')
    for line in pFile:
        fields = (line.strip()).split(" ")
        probandPos = int(fields[2])
        refCount = float(fields[0])
        altCount = float(fields[1])
        # variant belongs to the bext window
        if (probandPos > int(windowStats[window][END_IDX])):
            if (window != len(windowStats) - 1):            
                windowStats[window][4] = counter - 1
                windowStats[window + 1][3] = counter
                window += 1
        # getting the counts of reference alleles, alternative alleles and the SNP position
        variants[counter][0] = refCount
        variants[counter][1] = altCount
        variants[counter][2] = probandPos

        # get the reference's haplotypes
        for i in range(2):
            refList[counter][i] = refFullList[counter][int(refs[window][0])*2-2+i]
            refList[counter][i+2] = refFullList[counter][int(refs[window][1])*2-2+i]
        counter += 1
    pFile.close()
    windowStats[window][4] = len(refFullList) - 1
    
    return variants, refList, windowStats
