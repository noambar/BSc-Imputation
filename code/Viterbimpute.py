from HMMGetDataPerChromosome import *
from ViterbiPerWindow import *
import sys

STATES = [[0,2],[0,3],[1,2],[1,3]]

# using viterbi algorithm infer the final haplotypes

def main():#sampleDir, outputDir, method):

    sampleDir = sys.argv[1]
    
    for chromosome in range(1,23):#3):
        print ("chromosome: " + str(chromosome))
        viterbiList = []
        allRefLists = []
        for mode in range(3):
            print ("mode: " + str(mode))
            naivScoresFile = sampleDir + "/mode" + str(mode) + "/tagc128.shapeit." + str(chromosome) + ".naiv"
            naivScoresFileSorted = naivScoresFile + "_sorted"
            probandFile = sampleDir + "/chr" + str(chromosome) + ".out"

            # get pairs of best references for each chromosome in each mode
            refs = getPair(naivScoresFileSorted, "refferenceAsProband/", chromosome, mode)
            print ("done finding ref pairs")
            variants, refList, windowStats = getWindowsData(chromosome, probandFile, naivScoresFile, refs)
            print ("done getting data")


            viterbiStates = []
            for window in range(len(windowStats)):
                print (window)
                # adding the viterbi output foreach of the windows
                startIdx = windowStats[window][3]
                endIdx = windowStats[window][4]
                vs, vl = viterbi(variants[startIdx:endIdx+1], refList[startIdx:endIdx+1])
                viterbiStates += vs
            viterbiList.append(viterbiStates)
            allRefLists.append(refList)
            inferHapFromViterbiDirect(viterbiStates, refList, 'tests/' + sampleDir.split("/")[0] + '_chr' + str(chromosome) + '_viterbi_mode' + str(mode) + '.out')
        inferHapFromViterbiMajorVote(viterbiList, allRefLists, 'tests/' + sampleDir.split("/")[0] + '_chr' + str(chromosome) + '_viterbi_major.out')
    return






# infering the proband haplotypes for a chromosome using viterbi states as is. means that
# foreach mode their will be a different output file.
def inferHapFromViterbiDirect(viterbiStates, refList, outputPath):
    outputFile = open(outputPath, 'w')
    for var in range(len(viterbiStates)):
        firstHap = STATES[int(viterbiStates[var])][0]
        secondHap = STATES[int(viterbiStates[var])][1]
        outputFile.write(str(refList[var][firstHap]) + " " + str(refList[var][secondHap]) + "\n")
    outputFile.close()
    return


# infering the proband haplotypes for a chromosome using majority vote viterbi states.
# each variant is overlapped by n modes (n should be odd), and the haplotypes are determined
# by the majority votes.
def inferHapFromViterbiMajorVote(viterbiStates, refList, outputPath):
    outputFile = open(outputPath, 'w')
    for var in range(len(viterbiStates[0])):
        hapCount = np.zeros((4))
        for mode in range(3):
            firstHap = STATES[int(viterbiStates[mode][var])][0]
            secondHap = STATES[int(viterbiStates[mode][var])][1]
            firstHap = refList[mode][var][firstHap]
            secondHap = refList[mode][var][secondHap]
            hapCount[firstHap] += 1
            hapCount[secondHap + 2] += 1
        if (hapCount[0] > hapCount[1]):
            firstHap = 0
        else:
            firstHap = 1
        if (hapCount[2] > hapCount[3]):
            secondHap = 0
        else:
            secondHap = 1
        outputFile.write(str(firstHap) + " " + str(secondHap) + "\n")
    outputFile.close()
    return

# infering the proband haplotypes for a chromosome using weighted majority vote viterbi states.
# each variant is overlapped by n modes (n should be odd), and the haplotypes are determined
# by the weighted majority votes.
def inferHapFromViterbiMajorVoteWeighted(viterbiStates, refList, outputPath): ## not implemented yet ##
    outputFile = open(outputPath, 'w')
    for var in range(len(viterbiStates)):
        hapCount = np.zeros((4))
        for mode in range(3):
            firstHap = STATES[int(viterbiStates[mode][var])][0]
            secondHap = STATES[int(viterbiStates[mode][var])][1]
            firstHap = refList[mode][var][firstHap]
            secondHap = refList[mode][var][secondHap]
            hapCount[firstHap] += 1
            hapCount[secondHap + 2] += 1
        if (hapCount[0] > hapCount[1]):
            firstHap = 0
        else:
            firstHap = 1
        if (hapCount[2] > hapCount[3]):
            secondHap = 0
        else:
            secondHap = 1
        outputFile.write(str(firstHap) + " " + str(secondHap) + "\n")
    outputFile.close()
    return


main()
