import numpy as np

# This function reads a sam file and according to variations in given files, creates new haplotip
# files in a given directory.
# param samInput - a path to a sam file (file must be of format sam, and sorted by coordinates)
# param ouDir - an existing directory, output files will be written to that directory with names:
# chr$.out ($ is chromosome number)
# param hapDir - a directory that holds legend files of variation in the folowing format:
# var_chr$_% % T G ($ is the chromosome, % is the coordinates, T is the major allele and G is minor)
def createVarMatrix(samInput, outDir, hapDir):

    samFile = open(samInput, 'r')

    chromosome = ""
    legend = []
    posInLegend = 0
    newHap = np.zeros((1,1))
    lineCounter = 0

    for line in samFile:
        if (line.startswith("@")):
            continue
        fields = line.split("\t")
        # we have reached a new chromosome
        if (chromosome != fields[2]):
            # finished with chromosome - print values to output file
            if (chromosome != ""):
                # print "createVarMatrix chr" + chromosome ##todo delete
                outFile = open(outDir + "/chr" + chromosome + ".out", 'w')
                for i in range(lineCounter):
                    outFile.write(str(int(newHap[i][0])) + " " + str(int(newHap[i][1])) + " " + str(legend[i][0]) + "\n")
                outFile.close()

            chromosome = fields[2]
            legend = []
            lineCounter = 0

            if (chromosome.isdigit() is False): ## todo think about fixing this
                samFile.close()
                return

            # create a list of all variations in a chromosome
            for leg in open(hapDir + "/tagc128.shapeit." + chromosome + ".legend", 'r' ):
                if (leg.startswith("var")):
                    legend.append((int(leg.split(" ")[1]), leg.split(" ")[2], (leg.split(" ")[3]).split("\n")[0]))
                    lineCounter += 1
            newHap = np.zeros((lineCounter,2))
            posInLegend = 0
        if (posInLegend == lineCounter):
            continue

        # get start and end coordinates
        start = int(fields[3])
        readLen = len(fields[9])

        # read ends before variant
        if ((start + readLen) < legend[posInLegend][0]):
            continue

        # read overlaps variation
        elif (start <= legend[posInLegend][0] and (start + readLen - 1) >= legend[posInLegend][0]):
            currentLeg = 0
            # check more variations downstream
            while ((start + readLen - 1) >= legend[posInLegend + currentLeg][0]):
                seq = fields[9]
                SNPpos = (legend[posInLegend + currentLeg][0] - start)
                if (seq[SNPpos] == legend[posInLegend + currentLeg][1]):
                    newHap[posInLegend][0] += 1
                elif (seq[SNPpos] == legend[posInLegend + currentLeg][2]):
                    newHap[posInLegend][1] += 1
                currentLeg += 1
                if (posInLegend + currentLeg == lineCounter):
                    break

        # read starts after the variation
        elif (start > legend[posInLegend][0]):
            # promote the variation
            while (start > legend[posInLegend][0]):
                posInLegend += 1
                if (posInLegend == lineCounter):
                    break
            if (posInLegend == lineCounter):
                continue
            # read ends before the variation
            if ((start + readLen - 1) < legend[posInLegend][0]):
                continue
            # read overlaps variation
            elif (start <= legend[posInLegend][0] and (start + readLen - 1) >= legend[posInLegend][0]):
                currentLeg = 0
                while ((start + readLen - 1) >= legend[posInLegend + currentLeg][0]):
                    seq = fields[9]
                    SNPpos = (legend[posInLegend + currentLeg][0] - start)
                    if (seq[SNPpos] == legend[posInLegend + currentLeg][1]):
                        newHap[posInLegend][0] += 1
                    elif (seq[SNPpos] == legend[posInLegend + currentLeg][2]):
                        newHap[posInLegend][1] += 1
                    currentLeg += 1
                    if (posInLegend + currentLeg == lineCounter):
                        break

    samFile.close()

    # write last chromosome to output file
    # print "createVarMatrix chr" + chromosome ##todo delete
    outFile = open(outDir + "/chr" + chromosome + ".out", 'w')
    for i in range(lineCounter):
        outFile.write(str(int(newHap[i][0])) + " " + str(int(newHap[i][1])) + " " + str(legend[i][0]) + "\n")
    outFile.close()
    return