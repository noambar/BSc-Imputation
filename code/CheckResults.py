import sys
import numpy as np

HAPLO = ("0/0","0/1","1/1","1/2")

def main():
    inputFile = open(sys.argv[1], 'r')
    results = np.zeros((4,4))
    for line in inputFile:
        fields = line.strip().split(" ")
        vcf = int(fields[1])
        proband = int(fields[2])
        results[vcf][proband] += 1
    inputFile.close()
    outputFile = open(sys.argv[2], 'w')
    outputFile.write("vcf/proband\t0/0\t0/1\t1/1\t1/2\n")
    counter = 0
    for line in results:
        outputFile.write(HAPLO[counter] + "\t")
        for value in line:
            outputFile.write(str(int(value)) + "\t")
        outputFile.write(str(float(results[counter][counter])/sum(results[counter])) + "\n")
        counter += 1
    totCor = results[0][0] + results[1][1] + results[2][2]
    percent = totCor/results.sum()
    outputFile.write("total correct sights " + str(percent) + "\n")
    outputFile.close()
    return


main()
