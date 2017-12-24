__author__ = 'noam.bar'

HAP_CM_COL = 3
HAP_POS_COL = 1
MIN_WINDOW_SIZE = 3000000
MIN_CENTROMER_SIZE = 3000000

MODE_NORM = 0
MODE_ONE_THIRD = 1
MODE_TWO_THIRDS = 2
#### should be changed to - max between constant and centiMorgan.
#### and take in account the region of N's.

# this function iterates over a hapMap file and divides the chromosome into windows
# of the minimum between the centimorgan parameter or a constant window size.
# param hapMapDir - a directory that holds the hap map files in the format:
# genetic_map_GRCh37_chr$.txt ($ - the chromosome)
# param chromosome - the chromosome.
# param centiMorgan - a value of centi Morgan where the chromosome will be divided by.
# param mode - the offset in which to return the breaks: 0 - normal, 1 - one third in, 2 - two thirds in
def getHapMapPosition(hapMapDir, chromosome, centiMorgan, mode):
    gapsDir = 'HG19_Gaps_21_1_2016/'
    hapFilePath = "genetic_map_GRCh37_chr" + str(chromosome) + ".txt"
    hapFile = open(hapMapDir + "/" + hapFilePath, 'r')
    gapFilePath = gapsDir + "/" + "chr" + str(chromosome) + "_gaps_over_300000.gap"
    gapFile = open(gapFilePath, 'r')
    gaps = []
    nextGapExists = False
    for line in gapFile:
        fields = (line.strip()).split("\t")
        if (fields[0] != "0"):
            gaps.append((int(fields[0]), int(fields[1])))
            nextGapExists = True
	#print gaps
    positionInGaps = 0
    nextCMBreak = centiMorgan
    nextMinWinBreak = -1


    ###

    breaks = []

    ###
    for line in hapFile:
        if (line.startswith("chr")):
            fields = (line.strip()).split("\t")
            cm = float(fields[HAP_CM_COL])
            pos = int(fields[HAP_POS_COL])
            if (nextMinWinBreak == -1):
                nextMinWinBreak = pos + MIN_WINDOW_SIZE
            if (cm >= nextCMBreak and pos >= nextMinWinBreak):
                if (nextGapExists):
                    if (pos <= gaps[positionInGaps][0] and gaps[positionInGaps][0] - pos <= MIN_WINDOW_SIZE):
                        breaks.append(gaps[positionInGaps][0])
                        nextMinWinBreak = gaps[positionInGaps][1] + MIN_WINDOW_SIZE
                        nextCMBreak = cm + centiMorgan #nextCMBreak += centiMorgan 23.3.2016
                        positionInGaps += 1
                        if (len(gaps) == positionInGaps):
                            nextGapExists = False
                    else:
                        #print (cm) #todo delete
                        #print (pos) #todo delete
                        breaks.append(pos)
                        nextCMBreak = cm + centiMorgan #nextCMBreak += centiMorgan 23.3.2016
                        nextMinWinBreak = pos + MIN_WINDOW_SIZE #nextMinWinBreak += MIN_WINDOW_SIZE 23.3.2016
                else:
                    #print (cm) #todo delete
                    #print (pos) #todo delete
                    breaks.append(pos)
                    nextCMBreak = cm + centiMorgan #nextCMBreak += centiMorgan 23.3.2016
                    nextMinWinBreak = pos + MIN_WINDOW_SIZE #nextMinWinBreak += MIN_WINDOW_SIZE 23.3.2016
            else:
                if (nextGapExists):
                    if (pos <= gaps[positionInGaps][0] and gaps[positionInGaps][0] != 0 and gaps[positionInGaps][0]-pos < MIN_WINDOW_SIZE):
                        #print gaps[positionInGaps][0]
                        #print ("here")
                        breaks.append(gaps[positionInGaps][0])
                        nextMinWinBreak = gaps[positionInGaps][1] + MIN_WINDOW_SIZE
                        nextCMBreak = cm + centiMorgan #nextCMBreak += centiMorgan 23.3.2016
                        positionInGaps += 1
                        if (len(gaps) == positionInGaps):
                            nextGapExists = False
    ######                  

    if (int(mode) == MODE_NORM):
        for b in breaks:
            yield b

    allGaps = []
    for g in gaps:
        allGaps.append(g[0])
        allGaps.append(g[1])

    if (int(mode) == MODE_ONE_THIRD):
        for i in range(len(breaks) - 1):
            if (breaks[i] in allGaps):
                yield breaks[i]
            else:
                b = breaks[i+1]-breaks[i]
                yield int(breaks[i] + b/3)
        yield breaks[-1]

    if (int(mode) == MODE_TWO_THIRDS):
        for i in range(len(breaks) - 1):
            if (breaks[i] in allGaps):
                yield breaks[i]
            else:
                b = breaks[i+1]-breaks[i]
                yield int(breaks[i] + 2*b/3)
        yield breaks[-1]


    ######
    gapFile.close()            
    hapFile.close()
    return



