import re
import os
import argparse
import numpy as np

NUMMBER_OF_CHROMOSOMES = 22
NUM_OF_WINDOWS_IDX = 0
FATHER_COUNTER_IDX = 1
MOTHER_COUNTER_IDX = 2
BOTH_COUNTER_IDX = 3

# this file generates the trio's stats table
# sonDir - the trio's son output directory created by the naive algorithm
# fatherDir - the trio's father output directory created by the naive algorithm
# motherDir - the trio's mother output directory created by the naive algorithm
# mode - a number 0-2, the mode of the windows
# outputDir - a selected output directory for the output
# window - (optional) a positive number, thw window size in terms of the number of top scores at parents
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sonDir', type=str, help='a path to the son directory', default='no_file_given')
    parser.add_argument('fatherDir', type=str, help='a path to the father directory', default='no_file_given')
    parser.add_argument('motherDir', type=str, help='a path to the mother directory', default='no_file_given')
    parser.add_argument('mode', type=str, help='the mode of the naive algorithm', default='0')
    parser.add_argument('outputDir', type=str, help='a diretory to hold the output')
    parser.add_argument('--window', type=int, help='the size of window to search', default=5)
    args = parser.parse_args()
    mode = args.mode
    outputFileName = args.outputDir + "/tripleStats_" + args.sonDir.split("/")[0] + "_" + args.fatherDir.split("/")[0] + "_" + args.motherDir.split("/")[0] + "_mode" + mode + ".txt"
    if (os.path.exists(outputFileName)):
        print outputFileName + " allready exists!!!"
        exit(1)

    outputFile = open(outputFileName, 'w')
    outputFile.write("chromosome    number_of_windows   son-father  son-mother  son-both	%of_son_windows_covered_by_parents" + "\n")

    # go over each chromosome file of son, save best index in a matrix
    for chromosome in range(1, NUMMBER_OF_CHROMOSOMES + 1):
        sonFile = open(args.sonDir + "/mode" + mode + "/tagc128.shapeit." + str(chromosome) + ".naiv_sorted", 'r')
        sonScoreMatrix = []
        for line in sonFile:
            sonScoreMatrix.append(int((line.split(" ")[0]).split("_")[1]))
        sonFile.close()

        fatherFile = open(args.fatherDir + "/mode" + mode + "/tagc128.shapeit." + str(chromosome) + ".naiv_sorted", 'r')
        fatherScoreMatrix = []
        for line in fatherFile:
            values = line.split(" ")[0:args.window]
            ltemp = []
            for i in range(args.window):
                ltemp.append(int(values[i].split("_")[1]))
            fatherScoreMatrix.append(ltemp)

        motherFile = open(args.motherDir + "/mode" + mode + "/tagc128.shapeit." + str(chromosome) + ".naiv_sorted", 'r')
        motherScoreMatrix = []
        for line in motherFile:
            values = line.split(" ")[0:args.window]
            ltemp = []
            for i in range(args.window):
                ltemp.append(int(values[i].split("_")[1]))
            motherScoreMatrix.append(ltemp)

        sonFile.close()
        fatherFile.close()
        motherFile.close()

        statsMatrix = np.zeros(4)
        for row in range(len(sonScoreMatrix)):
            temp = 0
            if (sonScoreMatrix[row] in fatherScoreMatrix[row]):
                statsMatrix[FATHER_COUNTER_IDX] += 1
                temp = 1
            if (sonScoreMatrix[row] in motherScoreMatrix[row]):
                statsMatrix[MOTHER_COUNTER_IDX] += 1
                if (temp == 1):
                    statsMatrix[BOTH_COUNTER_IDX] += 1
            # outputFile.write(str(sonScoreMatrix[row]) + "\t" + str(fatherScoreMatrix[row]) + "\t" + str(motherScoreMatrix[row]) + "\n")
        statsMatrix[NUM_OF_WINDOWS_IDX] = len(sonScoreMatrix)
        numOfWin = int(statsMatrix[NUM_OF_WINDOWS_IDX])
        fatherWin = int(statsMatrix[FATHER_COUNTER_IDX])
        motherWin = int(statsMatrix[MOTHER_COUNTER_IDX])
        bothWin = int(statsMatrix[BOTH_COUNTER_IDX])
        percentWin = (fatherWin + motherWin - bothWin) / float(numOfWin)
        outputFile.write(str(chromosome) + "\t" + str(numOfWin) + "\t" + str(fatherWin) + "\t" + str(motherWin) + "\t" + str(bothWin)  + "\t" + str(percentWin) + "\n")

    outputFile.close()

    return



main()
