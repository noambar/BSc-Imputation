import numpy as np
import scipy.special
from HMMGetDataPerChromosome import *

# this is an implementation of the forward-backward algorithm


NUMBER_OF_STATES = 4
REF_HOMO = 0
HETERO = 1
ALT_HOMO = 2
VAR_POS_IDX = 2
MEAN_DISTANCE_BETWEEN_VARIANTS = 300.

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


def forwardBackward(variants, refs):
    # count number of variants in window (in the future know that number in advance)
    numberOfVariants = 0
    for line in variants:
        if (line[VAR_POS_IDX] == 0):
            break
        numberOfVariants += 1

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
    # go over backward table from end to begining and fill it
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


    for i in forwardTable:
        print (i)
    for i in backwardTable:
        print (i)

    return backwardTable, backwardTable







# emissions are calculated using binomial distribution
def calculateEmission(proband, state, refs):
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



    



### something doesn't work, viterbi table ends with -inf !!! also look at the transitions and check if they make sense



def main():
    r = [[0,1,0,0],[0,1,1,0],[0,1,0,1]]
    v = [[1,0,300],[0,1,2000],[0,2,2100],[0,0,0]]
    #viterbi(v,r)
    forwardBackward(v,r)
    # print (calculateEmission((1,1,0), 0, np.array([0,0,1,0])))
    return
    sampleDir = "HG002.hs37d5.60x.1_down_sampled_0.5x_19:25:42.07399/"
    chromosome = 1
    mode = 0
    naivScoresFile = sampleDir + "/mode" + str(mode) + "/tagc128.shapeit." + str(chromosome) + ".naiv"
    naivScoresFileSorted = naivScoresFile + "_sorted"
    probandFile = sampleDir + "/chr" + str(chromosome) + "_probandMat.txt"
    # get pairs of best references for each chromosome in each mode
    refs = getPair(naivScoresFileSorted, "refferenceAsProband/", chromosome, mode)
    print ("done finding ref pairs")
    variants, refList = getWindowsData(chromosome, probandFile, naivScoresFile, refs)
    print ("done getting window data")
    #viterbi(variants[0], refList[0])

    forwardBackward(variants[0], refList[0])
    #for window in range(len(variants)):
    #    print (window)
    #    viterbiTable, traceBack, viterbiStates = viterbi(variants[window], refList[window])
    return


main()
