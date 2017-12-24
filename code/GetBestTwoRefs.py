SECOND_WITH_IBD = -1
SIMILAR_RANGE = 5

# this funcion gets the best scoring pair of reference individuals for a chromosome in a giving mode
# sorted file - the sorted scores file
# referenceDir - the reference directory holding the scores of individulas between themselfs
# chromosome - a number between 1-22
# mode - currently a numebr between 0-2

def getPair(sortedFile, referenceDir, chromosome, mode):
    secondRef = SECOND_WITH_IBD
    sortedFile = open(sortedFile, 'r')
    windowCounter = 1
    pairs = []
    for line in sortedFile:
        fields = line.strip().split(" ")
        firstRef = fields[0].split("_")[1]
        compareRefFile = open(referenceDir + "/" + "ref" + firstRef + "/mode" + str(mode) + "/tagc128.shapeit." + str(chromosome) + ".naiv_sorted", 'r')
        lines = compareRefFile.readlines()
        compareRefFile.close()
        best = []
        for i in range(1,SIMILAR_RANGE+1):
            best.append(lines[windowCounter-1].split(" ")[i].split("_")[1])
        #print (best)
        secondNotFound = True
        posInSorted = 1
        while (secondNotFound):
            nextRef = fields[posInSorted].split("_")[1]
            if nextRef in best:
                posInSorted += 1
                #print (posInSorted)
                continue
            secondRef = nextRef
            break
        if (posInSorted >= SIMILAR_RANGE):
            print ("all ten best scores are similar")
            secondRef = fields[1].split("_")[1]
        pairs.append((firstRef, secondRef))
        #print (pairs[-1])
    return pairs
