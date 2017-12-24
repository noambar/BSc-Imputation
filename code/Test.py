import numpy as np
from GetHapMapPositions import *
from CreateVarMatrix import *
import os
from GetBestTwoRefs import *
import matplotlib.pyplot as plt

# this file was only used in order to test python functions
def main():

    # for it in iteratorfunc():
    #     print it

    #if (5 in [1,2,3,4]):
    #    print "yes"

    # m = min(7,9,8,6,4,3,2,5,6,7,8)
    #
    # print m
    
    a = getHapMapPosition('HapMap_GRCh37/', 1, 3, 2)
    for b in a:
        print (b)
    #   print (nextBreak)
    # createVarMatrix('HG002.hs37d5.down_sampled_1x.sam', 'test', 'impute2/')
    # if (os.path.exists("tempdir")):
    #     print "dir allready exists"
    #     exit(1)
    # os.makedirs("tempdir")
    #
    # output = open("tempdir/file.txt", 'w')
    #
    # output.write("aaaa\n")
    #
    # os.system("sbatch --wrap=\"echo noambar\"")
    #p = getPair("HG002.hs37d5.60x.1_down_sampled_1.0x_19:25:30.88795/mode0/tagc128.shapeit.8.naiv_sorted", "refferenceAsProband/", 8, 0)
    #for a in p:
    #    print (a)



    return


def iteratorfunc():
    l = [1, 2, 3, 4]
    for y in l:
        yield y


main()
