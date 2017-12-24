import matplotlib.pyplot as plt

# this program when running it creates a figure that show the distribution of
# window sizes across chromosomes, modes and overall, in directory windowSizes
def main():
    total = []
    for m in range(3):
        mode = []
        for c in range(1,23):
            chrom = []
            f = open("HG004.hs37d5.60x.1_down_sampled_0.1x/mode" + str(m) + "/tagc128.shapeit." + str(c) + ".naiv", 'r')
            for line in f:
                fields = line.split(" ")
                chrom.append(int(fields[1]) - int(fields[0]))
            f.close()
            plt.figure()
            plt.hist(chrom, 100, facecolor='green')
            plt.xlabel("Length of window [bps]")
            plt.ylabel("Number of windows")
            plt.grid(True)
            #print (c)
            #print (m)
            #print (chrom)
            plt.savefig("windowSizes/WLD_chr" + str(c) + "_mode" + str(m) + ".png")
            mode += chrom
            chrom = []
        plt.figure()
        plt.hist(mode, 1000, facecolor='blue')
        plt.xlabel("Length of window [bps]")
        plt.ylabel("Number of windows")
        plt.grid(True)
        plt.savefig("windowSizes/WLD_mode" + str(m) + ".png")
        total += mode
        mode = []
    plt.figure()
    plt.hist(total, 1000, facecolor='red')
    plt.xlabel("Length of window [bps]")
    plt.ylabel("Number of windows")
    plt.grid(True)
    plt.savefig("windowSizes/WLD_all.png")
    return




main()
