import numpy as np
from GetHapMapPositions import *
import sys
import os

SIZE_OF_REFERENCE = 128
PROBAND_COUNTER = 0
PROBAND_POS = 1
PROBAND_REF = 2
PROBAND_ALT = 3
REF_ALLELE = 0
ALT_ALLELE = 2
HET_ALLELE = 1
HOMO_THRESHOLD = 10
ALL_ALLELES_KNOWN = 0
ONE_ALLELE_MASKED = 1
CONSTANT_MISSMATCH_ALLELES = -1


countMat = np.zeros(SIZE_OF_REFERENCE) ###

# using the naiv algorithm
def main():
    naivImputeAlgo(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6], sys.argv[7])
    return

# this function implements the naiv scoring imputing algorithm. going over windows of size
# min(centiMorgan, fixed_window_size) and calculates scores for each corresponding window
# in the reference panel.
# @param probandMatPath - a path to a file holding only relevant data about the proband variants.
# @param hapDir - a path to the reference haplotype directory.
# @param centiMorgan - size in centiMorgan to divide the chromosome by.
# @param hapMapDir - a path to the hapMap directory.
# @param chromosome - the chromosome number.
# @param outputDir - a path to the output directory.
def naivImputeAlgo(probandMatPath, hapDir, centiMorgan, hapMapDir, chromosome, outputDir, mode):
    # open relevant file for reading and writing, ## open files for writing will override existing files with the same name ##
    asPairFile = open(hapDir + "/tagc128.shapeit." + str(chromosome) + ".as_pair", 'r')
    freqFile = open(hapDir + "/tagc128.shapeit." + str(chromosome) + ".freq", 'r')
    outputFile = open(outputDir + "/tagc128.shapeit." + str(chromosome) + ".naiv", 'w')
    countFile = open(outputDir + "/tagc128.shapeit." + str(chromosome) + ".count", 'w') ###
    methodFile = open(outputDir + "/tagc128.shapeit." + str(chromosome) + ".method", 'w') ### ###
    methodCounter = np.zeros(5) ### ###

    probandMatFile = open(probandMatPath, 'r')

    # create the proband matrix by reading it from file
    probandMat = []
    for line in probandMatFile:
        fields = (line.strip()).split(" ")
        probandMat.append((int(fields[PROBAND_COUNTER]), int(fields[PROBAND_POS]), int(fields[PROBAND_REF]), int(fields[PROBAND_ALT])))
    probandMatFile.close()
    coordinates = probandMat[0][PROBAND_POS]
    freq_pair_counter = 0
    positionInProbandMat = 0
    pairLine = asPairFile.readline()
    freqLine = freqFile.readline()
    variantCounterPerWindow = 0
    # get next break of chromosome
    for nextBreak in getHapMapPosition(hapMapDir, chromosome, centiMorgan, mode): ### mode added
        # writing starting coordinates of window
        outputFile.write(str(coordinates) + " ")

        # creating a vector for reference scores for each window
        referenceScoreMat = np.zeros(SIZE_OF_REFERENCE)
        coordinates = probandMat[positionInProbandMat][PROBAND_POS]

        # taking all variations before the next break
        while (coordinates <= nextBreak):
            # finding next existing variant in the proband
            # getting to the relevant rows in the attached files
            while (freq_pair_counter < probandMat[positionInProbandMat][PROBAND_COUNTER]):
                pairLine = asPairFile.readline()
                freqLine = freqFile.readline()
                freq_pair_counter += 1

            ref = probandMat[positionInProbandMat][PROBAND_REF]
            alter = probandMat[positionInProbandMat][PROBAND_ALT]

            # using a helper function in order to calculate score for each reference sample in the window
            referenceScoreMat, methodCounter = calculateAlleleScore(ref, alter, pairLine, freqLine, referenceScoreMat, methodCounter) ### ###
            positionInProbandMat += 1
            variantCounterPerWindow += 1

            # in case last existing variant in proband matrix appears before the next break, moving on to next window
            if (positionInProbandMat == len(probandMat)):
                break

            coordinates = probandMat[positionInProbandMat][PROBAND_POS]

        # writing end of window coordinates
        outputFile.write(str(nextBreak) + " ")
        outputFile.write(str(variantCounterPerWindow) + " ")
        variantCounterPerWindow = 0
        # print window scores to file
        for i in range(SIZE_OF_REFERENCE - 1):
            outputFile.write(str(referenceScoreMat[i]) + " ")
            countFile.write(str(countMat[i]) + " ") ###
            countMat[i] = 0 ###
        outputFile.write(str(referenceScoreMat[-1]) + "\n")
        countFile.write(str(countMat[-1]) + "\n") ###
        countMat[-1] = 0 ###

    # in case last existing variant in last window in proband matrix appears before the last break, finish
    if (positionInProbandMat == len(probandMat)):
        outputFile.close()
        freqFile.close()
        asPairFile.close()
        # should write method file because otherwise it wont be wrriten
        ########
        methodFile.write("all_alleles_known_reference	all_alleles_known_alt	all_alleles_known_hetro	one_allele_masked_ref	one_allele_masked_alt" + "\n") ### ###
        methodFile.write(str(methodCounter[0]) + "\t" + str(methodCounter[1]) + "\t" + str(methodCounter[2]) + "\t" + str(methodCounter[3]) + "\t" + str(methodCounter[4]) + "\n") ### ### 
        methodFile.close() ### ###

        ### add a sorting code ###
        sortResults(outputDir, chromosome)
        ###########################
        print ("algorithm run complete successfully")
        ########
        return

    outputFile.write(str(coordinates) + " ")
    referenceScoreMat = np.zeros(SIZE_OF_REFERENCE)

    # going over last partial window in chromosome
    for prob in range(positionInProbandMat, len(probandMat)):
        # getting to the relevant rows in the attached files
        while (freq_pair_counter < probandMat[prob][PROBAND_COUNTER]):
            pairLine = asPairFile.readline()
            freqLine = freqFile.readline()
            freq_pair_counter += 1
        ref = probandMat[prob][PROBAND_REF]
        alter = probandMat[prob][PROBAND_ALT]
        # using a helper function in order to calculate score for each reference sample in last window
        referenceScoreMat, methodCounter = calculateAlleleScore(ref, alter, pairLine, freqLine, referenceScoreMat, methodCounter)  ### ###
        variantCounterPerWindow += 1
    outputFile.write(str(probandMat[-1][PROBAND_POS]) + " ")
    outputFile.write(str(variantCounterPerWindow) + " ")
    for i in range(SIZE_OF_REFERENCE - 1):
        outputFile.write(str(referenceScoreMat[i]) + " ")
        countFile.write(str(countMat[i]) + " ") ###
        countMat[i] = 0 ###
    outputFile.write(str(referenceScoreMat[-1]) + "\n")
    countFile.write(str(countMat[-1]) + "\n") ###
    countMat[-1] = 0 ###
    outputFile.close()
    countFile.close() ###
    methodFile.write("all_alleles_known_reference	all_alleles_known_alt	all_alleles_known_hetro	one_allele_masked_ref	one_allele_masked_alt" + "\n") ### ###
    methodFile.write(str(methodCounter[0]) + "\t" + str(methodCounter[1]) + "\t" + str(methodCounter[2]) + "\t" + str(methodCounter[3]) + "\t" + str(methodCounter[4]) + "\n") ### ### 
    methodFile.close() ### ###
    freqFile.close()
    asPairFile.close()

    ### add a sorting code ###
    sortResults(outputDir, chromosome)
    ###########################
    print ("algorithm run complete successfully")

    return



# this function determines if a variant is as the reference allele or as the alternative allele,
# and by looking at the reference haplotype, determines what scoring method to use.
# @param refCount - the number of times the reference allele appeared in the proband.
# @param altCount - the number of times the alternative allele appeared in the proband.
# @param refPairLine - the appropriate line in the pair file.
# @param refFreqLine - the appropriate line in the freq file.
# @param referenceScoreMat - the reference's scoring matrix to be updated
# returns the reference scoring matrix after update.
def calculateAlleleScore(refCount, altCount, refPairLine, refFreqLine, referenceScoreMat, methodCounter): ### ###
    # it is clear that the proband is homozygous for the alternative allele
    if (refCount == 0 and altCount > HOMO_THRESHOLD or (refCount > 0 and altCount/float(refCount) > HOMO_THRESHOLD)):
        probandAllele = ALT_ALLELE
        method = ALL_ALLELES_KNOWN
        methodCounter[1] += 1 ### ###
    # it is clear that the proband is homozygous for the reference allele
    elif (altCount == 0 and refCount > HOMO_THRESHOLD or (altCount > 0 and refCount/float(altCount) > HOMO_THRESHOLD)):
        probandAllele = REF_ALLELE 
        method = ALL_ALLELES_KNOWN
        methodCounter[0] += 1 ### ###
    # probably one allele is missing and there are only alternative reads
    elif (altCount == 0 and refCount <= HOMO_THRESHOLD):
        probandAllele = REF_ALLELE ############ major bug might happend - original: probandAllele = ALT_ALLELE ############
        method = ONE_ALLELE_MASKED
        methodCounter[3] += 1 ### ###
    # probably one allele is missing and there are only reference reads
    elif (refCount == 0 and altCount <= HOMO_THRESHOLD):
        probandAllele = ALT_ALLELE ############ major bug might happend - original: probandAllele = REF_ALLELE ############
        method = ONE_ALLELE_MASKED
        methodCounter[4] += 1 ### ###
    # proband is likely to be heterozygous
    else:
        probandAllele = HET_ALLELE
        method = ALL_ALLELES_KNOWN
        methodCounter[2] += 1 ### ###
    referenceFields = (refPairLine.strip()).split(" ")
    alterAlleleFreq = float(refFreqLine.strip())
    for i in range(len(referenceFields)):
        referenceScoreMat[i] += getAlleleProbability(probandAllele, method, alterAlleleFreq, int(referenceFields[i]), i) ###
    return referenceScoreMat, methodCounter ### ###



# this function returns a score for each variant for each sample in the reference.
# @param probandAllele - the proband allele (either REF_ALLELE|ALT_ALLELE|HET_ALLELE)
# @param method - the method by which to calculate the score (either ALL_ALLELES_KNOWN|ONE_ALLELE_MASKED)
# @param refAlterAlleleFreq - the reference's alternative allele frequency.
# @param referenceAllele - the reference allele (either REF_ALLELE|ALT_ALLELE|HET_ALLELE)
# returns the score calculated by the parameters.
def getAlleleProbability(probandAllele, method, refAlterAlleleFreq, referenceAllele, i): ###
    alleleScore = 0
    Q = refAlterAlleleFreq
    P = 1 - Q
    if (method == ALL_ALLELES_KNOWN):
        if (probandAllele == REF_ALLELE and referenceAllele == REF_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./P)
        elif (probandAllele == REF_ALLELE and referenceAllele == HET_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./(2 * P))
        elif (probandAllele == ALT_ALLELE and referenceAllele == ALT_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./Q)
        elif (probandAllele == ALT_ALLELE and referenceAllele == HET_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./(2 * Q))
        elif (probandAllele == HET_ALLELE and referenceAllele == REF_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./(2 * P))
        elif (probandAllele == HET_ALLELE and referenceAllele == ALT_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./(2 * Q))
        elif (probandAllele == HET_ALLELE and referenceAllele == HET_ALLELE):
            countMat[i] += 1 ###
            return np.log(1./(4 * Q * P))
        elif (probandAllele == REF_ALLELE and referenceAllele == ALT_ALLELE or probandAllele == ALT_ALLELE and referenceAllele == REF_ALLELE):
            return CONSTANT_MISSMATCH_ALLELES

    elif (method == ONE_ALLELE_MASKED):
        if (probandAllele == REF_ALLELE and referenceAllele == REF_ALLELE):
            countMat[i] += 1 ###
            return np.log((P**3 + ((P**2)*Q)/2.)/(P**3))
        elif (probandAllele == REF_ALLELE and referenceAllele == ALT_ALLELE):
            return np.log(1./2)
        elif (probandAllele == REF_ALLELE and referenceAllele == HET_ALLELE):
            countMat[i] += 1 ###
            return np.log(((P**2)*Q + P*Q/2)/(2*(P**2)*Q))
        elif (probandAllele == ALT_ALLELE and referenceAllele == REF_ALLELE):
            return np.log(1./2)
        elif (probandAllele == ALT_ALLELE and referenceAllele == ALT_ALLELE):
            countMat[i] += 1 ###
            return np.log((Q**3 + ((Q**2)*P)/2.)/(Q**3))
        elif (probandAllele == ALT_ALLELE and referenceAllele == HET_ALLELE):
            countMat[i] += 1 ###
            return np.log(((Q**2)*P + P*Q/2)/(2*(Q**2)*P))

    return alleleScore

# sorting function
def sortResults(outputDir, chromosome):
    fileToSort = open(outputDir + "/tagc128.shapeit." + str(chromosome) + ".naiv", 'r')
    sortFile = open(outputDir + "/tagc128.shapeit." + str(chromosome) + ".naiv_sorted", 'w')
    for line in fileToSort:
        values = (line.strip()).split(" ")[3:]
        npval = np.zeros((len(values)))
        for val in range(len(values)):
            npval[val] = float(values[val])   
        s = npval.argsort()
        k = []
        for i in range(len(npval)):
            k.append(str(npval[s[i]]) + "_" + str(s[i]+1))
        k.reverse()

        for i in range(len(k) - 1):
            sortFile.write(k[i] + " ")
        sortFile.write(k[-1] + "\n")
    sortFile.close()
    fileToSort.close()
    return

main()
