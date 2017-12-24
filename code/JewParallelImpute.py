from ExtractProbandVariations import *
from CreateVarMatrix import *
import argparse

# this main function is used in order to run the navie algorithm. lines that are with commetns (#) should be uncomment, in order for
# the algorithm to run fully.
# sum - the alignment file
# hapDir - the impute2 dir
# outPutDir - a selected uotput directory
# hapMapDir - the directory holding the HapMap data
# centimorgan - (optional) a positive number for the window size
def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('sam', type=str, help='a path to the sam file', default='no_file_given')
    parser.add_argument('hapDir', type=str, help='a path to the Dir holding the refference haplotypes', default='no_file_given')
    parser.add_argument('outPutDir', type=str, help='a path to wanted output Dir', default='no_file_given')
    parser.add_argument("hapMapDir", type=str,  help="a path to the directory of HapMap files")
    parser.add_argument("--centiMorgan", type=str,  help="a parameter to be used when dividing the chromosome to windows", default=3)
    args = parser.parse_args()

    #if (os.path.exists(args.outPutDir)):       # this line should be uncomment
    #    print args.outPutDir + " allready exists!!!"       # this line should be uncomment
    #    exit(1)       # this line should be uncomment

    #os.makedirs(args.outPutDir)       # this line should be uncomment

    for i in range(3):
        os.makedirs(args.outPutDir + "/mode" + str(i))

    print "creating var files"
    # create the SNP calling for the proband and store the files in an output directory
    #createVarMatrix(args.sam, args.outPutDir, args.hapDir)       # this line should be uncomment

    print "creating proband matix"
    # create a proband matrix holding only SNPs that had calls in the proband
    #extractProbandVariations(args.outPutDir)       # this line should be uncomment

    # for each chromosome run the naiv imputing algorithm in sbatch
    for probandFile in os.listdir(args.outPutDir):
        if (probandFile.startswith("chr") and probandFile.endswith("_probandMat.txt")):
            m = re.search('chr([0-9]+)_probandMat\.txt', probandFile)
            chromosome =  m.group(1)
            for i in range(3):
                os.system("sbatch -c 2 --mem-per-cpu=8000 --time=1-0 -o naiv-algo_" + args.outPutDir.split("/")[0] + "_chr" + chromosome + "-%j.out --wrap=\"python NaivImputeAlgorithm.py " + args.outPutDir + "/" + probandFile + " " + args.hapDir + " " + str(args.centiMorgan) + " " +  args.hapMapDir + " " + chromosome + " " + args.outPutDir + "/mode" + str(i) + " " + str(i) + "\" ")
    return


main()
