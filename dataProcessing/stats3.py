import params
from utils import utils
import numpy as np
from utils import utils
import params
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Process, Value, Queue
import time

P_THRESHOLD = 0.05


def merge3():
    fin = open("%s/FSUBTEST/3/FileMap.txt" % params.FADER_OUT)
    fout = open("%s/FSUBTEST/3/3.txt" % params.FADER_OUT, "w")
    dDrug2Se = dict()

    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        hashFile = parts[0]
        f = open("%s/FSUBTEST/3/%s" % (params.FADER_OUT, hashFile))
        while True:
            l = f.readline()
            if l == "":
                break
            parts = l.strip().split("_")
            drug = parts[0]
            se = parts[1].split("\t")[0]
            ses = utils.get_insert_key_dict(dDrug2Se, drug, [])
            ses.append(se)

        f.close()
    for k, v in dDrug2Se.items():
        fout.write("%s\t%s\n" % (k, ",".join(v)))
    fout.close()


def stats3():
    fin = open("%s/FSUBTEST/3/FileMap.txt" % params.FADER_OUT)
    dComCount = dict()

    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        hashFile = parts[0]
        f = open("%s/FSUBTEST/3/%s" % (params.FADER_OUT, hashFile))
        while True:
            l = f.readline()
            if l == "":
                break
            parts = l.strip().split("_")
            drug = parts[0]
            se = parts[1].split("\t")[0]
            utils.add_dict_counter(dComCount, drug)


        f.close()

    counts = dComCount.values()
    maxCount = max(counts)
    dCountFreq = dict()
    for c in counts:
        utils.add_dict_counter(dCountFreq, c)

    arCountFreq = np.zeros(maxCount)
    for c,v in dCountFreq.items():
        arCountFreq[c-1] = v
    x = np.arange(1, maxCount + 1)

    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,2,1)
    ax.scatter(x, arCountFreq)
    plt.xlabel("3-Drug Combinations with Frequencies")
    plt.ylabel("Number of Combinations")
    sum = np.sum(arCountFreq)
    prop = arCountFreq / sum

    ax2 = ax.twinx()

    ax2.scatter(x, prop)
    plt.ylabel("Percentage")
    plt.xlim(1,5)
    # fig.add_subplot(1,2,2)
    # plt.scatter(x, prop)
    # plt.xlabel("Freq of 3 Drug Combinations")
    # plt.ylabel("Proportion")
    plt.title("%s Combinations (%s)" % (np.sum(arCountFreq[:5]), round(np.sum(prop[:5]), 2)))
    plt.tight_layout()

    plt.savefig("%s/D3_Freq.png" % params.FIG_DIR)



if __name__ == "__main__":
    # merge3()
    stats3()