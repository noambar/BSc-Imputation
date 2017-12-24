import numpy as np
import scipy.special

from HMMGetDataPerChromosome import *

# a viterbi algorithm implementaion for each window of each chromosome

NUMBER_OF_STATES = 4
REF_HOMO = 0
HETERO = 1
ALT_HOMO = 2
VAR_POS_IDX = 2
MEAN_DISTANCE_BETWEEN_VARIANTS = 260.

PEP = 0.01 # PROBABILITY_FOR_PHASING_ERROR = 300/30000 # temporary value

SEP = 0.01 # PROBABILITY_FOR_SEQUENCE_ERROR = 0.01

# define a transition matrix - # also consult with Shai
TRANSITIONS = np.array([[0, 1, 1, 2],
                        [1, 0, 2, 1],
                        [1, 2, 0, 1],
                        [2, 1, 1, 0]])

# there are four states (numbered from 0..3):
# 0 - R1H1,R2H1
# 1 - R1H1,R2H2
# 2 - R1H2,R2H1
# 3 - R1H2,R2H2
# Rx - x(1 or 2) represent the number of reference among the two references.
# Hy - y(1 or 2) stands for the number of haplotype out of the reference two haplotypes.


def viterbi(variants, refs):
    numberOfVariants = len(variants)

    viterbiTable = np.zeros((NUMBER_OF_STATES, numberOfVariants))
    traceBack = np.zeros((NUMBER_OF_STATES, numberOfVariants))
    viterbiTable += -np.inf

    # fill the first column as with just using the emissions as score
    for state in range(NUMBER_OF_STATES):
        viterbiTable[state][0] = np.log(calculateEmission(variants[0], state, refs[0]))

    for pos in range(1,numberOfVariants):
        transitions = calculateTransition(variants[pos][VAR_POS_IDX] - variants[pos-1][VAR_POS_IDX])
        for state in range(NUMBER_OF_STATES):
            emission = calculateEmission(variants[pos], state, refs[pos])
            options = []
            for lastState in range(NUMBER_OF_STATES):
                options.append(viterbiTable[lastState][pos-1] + np.log(transitions[TRANSITIONS[lastState][state]]))
            viterbiTable[state][pos] = max(options) + np.log(emission)
            traceBack[state][pos] = np.argmax(options)

    options = []
    for state in range(NUMBER_OF_STATES):
        options.append(viterbiTable[state][-1])
    bestEnding = np.argmax(options)
    viterbiLikelihood = options[bestEnding]
    
    viterbiStates = [bestEnding]
    for reversePos in range(numberOfVariants-1, 0, -1):
        bestEnding = traceBack[bestEnding][reversePos]
        viterbiStates.insert(0, bestEnding)

    #for i in viterbiTable:
    #    print (i)
    #for i in traceBack:
    #    print (i)
    #for i in viterbiStates:
    #    print (i)

    return viterbiStates, viterbiLikelihood








def forwardBackward(variants, refs):
    numberOfVariants = len(variants)

    # initialize tables to hold the forward and backward values
    forwardTable = np.zeros((NUMBER_OF_STATES, numberOfVariants))
    backwardTable = np.zeros((NUMBER_OF_STATES, numberOfVariants))
    forwardTable += -np.inf
    backwardTable += -np.inf

    # fill the first column as with just using the emissions as score
    for state in range(NUMBER_OF_STATES):
        forwardTable[state][0] = np.log(calculateEmission(variants[0], state, refs[0]))

    for pos in range(1,numberOfVariants):
        # transitions are calculated once for all options, using only the distance
        transitions = calculateTransition(variants[pos][VAR_POS_IDX] - variants[pos-1][VAR_POS_IDX])
        for state in range(NUMBER_OF_STATES):
            nextVal = -np.inf
            emission = calculateEmission(variants[pos], state, refs[pos])
            # foreach position go over all states in the previous position add up the probabilities
            for lastState in range(NUMBER_OF_STATES):
                if (nextVal == -np.inf):
                    nextVal = forwardTable[lastState][pos-1] + np.log(transitions[TRANSITIONS[lastState][state]])
                else:
                    nextVal = np.logaddexp(nextVal, forwardTable[lastState][pos-1] + np.log(transitions[TRANSITIONS[lastState][state]]))
            forwardTable[state][pos] = nextVal + np.log(emission)

    # at the last position, assign 0 to all states
    for state in range(NUMBER_OF_STATES):
        backwardTable[state][-1] = 0
    # go over backward table from end to beginning and fill it
    for pos in range(numberOfVariants-2,-1,-1):
        transitions = calculateTransition(variants[pos][VAR_POS_IDX] - variants[pos-1][VAR_POS_IDX])
        for state in range(NUMBER_OF_STATES):
            nextVal = -np.inf
            for lastState in range(NUMBER_OF_STATES):
                # calculate emission for the next position
                emission = calculateEmission(variants[pos+1], lastState, refs[pos+1])
                # foreach position go over all states in the next position and add up probabilities
                if (nextVal == -np.inf):
                    nextVal = backwardTable[lastState][pos+1] + np.log(transitions[TRANSITIONS[state][lastState]]) + np.log(emission)
                else:
                    nextVal = np.logaddexp(nextVal, backwardTable[lastState][pos+1] + np.log(transitions[TRANSITIONS[state][lastState]]) + np.log(emission))
            backwardTable[state][pos] = nextVal


    #for i in forwardTable:
    #    print (i)
    #for i in backwardTable:
    #    print (i)

    return forwardTable, backwardTable






# emissions are calculated using binomial distribution
def calculateEmission(proband, state, refs):
    if (proband[0] == 0 and proband[1] == 0): ## added in order to deal with variants that has no calls
        return 1 ##
    binom = scipy.special.binom(proband[0]+proband[1],proband[0])
    firstHap = refs[0:2][int(state/2)]
    secondHap = refs[2:4][state%2]
    error = 0
    if (firstHap + secondHap == 0):
        error = ((1-SEP)**proband[0]) * (SEP**proband[1])
    elif (firstHap + secondHap == 1):
        error = 2**(-(proband[0]+proband[1]))
    elif (firstHap + secondHap == 2):
        error = ((1-SEP)**proband[1]) * (SEP**proband[0])
    return binom*error

# transitions are calculated using the distance between the two variants
def calculateTransition(distance):
    # calculating the transitions according to the distance between the two variants
    p = (1 - (1-PEP)**(distance / MEAN_DISTANCE_BETWEEN_VARIANTS))
    # transition value greater than 0.5 should be taken as random (0.5)
    if (p >= 0.5):
        return [0.25, 0.25, 0.25]
    stay = (1-p)**2
    oneSwitch = p*(1-p)
    doubleSwitch = p**2
    return [stay, oneSwitch, doubleSwitch]



    





def main():
    #r = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,1,0,1],[0,1,0,1],[0,1,0,1],[0,1,0,1],[0,1,0,1],[0,1,0,1],[0,1,0,1]]
    #v = [[1,0,300],[2,0,700],[1,0,1500],[1,0,2000],[2,0,2100],[2,0,2500],[2,0,2900],[2,0,3400],[2,0,5000],[2,0,8000],[0,0,0]]
    #viterbi(v,r)
    # print (calculateEmission((1,1,0), 0, np.array([0,0,1,0])))
    #return
    sampleDir = "HG002.hs37d5.60x.1_down_sampled_0.2x_15:53:16.35632/"
    chromosome = 1
    mode = 0
    naivScoresFile = sampleDir + "/mode" + str(mode) + "/tagc128.shapeit." + str(chromosome) + ".naiv"
    naivScoresFileSorted = naivScoresFile + "_sorted"
    probandFile = sampleDir + "/chr" + str(chromosome) + ".out"
    # get pairs of best references for each chromosome in each mode
    refs = getPair(naivScoresFileSorted, "refferenceAsProband/", chromosome, mode)
    print ("done finding ref pairs")
    variants, refList, windowStats = getWindowsData(chromosome, probandFile, naivScoresFile, refs)
    print ("done getting window data")
    startIdx = windowStats[0][3]
    endIdx = windowStats[0][4]
    print (refs[0])
    
    v,l = viterbi(variants[startIdx:endIdx+1], refList[startIdx:endIdx+1])
    f,b = forwardBackward(variants[startIdx:endIdx+1], refList[startIdx:endIdx+1])
    return
    print (v[0:30])
    print (refList[startIdx:startIdx+30])
    print (v[2000:2030])
    print (refList[2000:2030])
    #startIdx = windowStats[1][3]
    #endIdx = windowStats[1][4]
    #print (refs[1])
    #viterbi(variants[startIdx:endIdx+1], refList[startIdx:endIdx+1])
    #startIdx = windowStats[2][3]
    #endIdx = windowStats[2][4]
    #print (refs[2])
    #viterbi(variants[startIdx:endIdx+1], refList[startIdx:endIdx+1])
    return
    for window in range(len(windowStats)):
        #print (window)
        #print (refs[window])
        startIdx = windowStats[window][3]
        endIdx = windowStats[window][4]
        #print (str(startIdx) + " - " + str(endIdx))
        #print (refList[startIdx])
        viterbi(variants[startIdx:endIdx+1], refList[startIdx:endIdx+1])
    return


main()
