import params
from utils import utils

OUT_DIR = params.CAD_OUT
PREF = "SCAD"

def writeSortedDictC(dl, p):
     f = open(p, "w")
     for kv in dl:
         k, v = kv
         f.write("%s\t%s\n" % (k, v))
     f.close()

def exportPair():
    fin = open("%s/CADER.txt" % OUT_DIR)
    # fout = open("%s/JADERIndPair.txt" % params.JADER_OUT, "w")
    validDrugs = dict()
    validPairs = dict()
    validIndicates = dict()
    validSes = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("$")
        drugComb = parts[1]
        indications = parts[2]
        ses = parts[3]
        drugs = drugComb.split(",")
        # print(drugs)
        for drug in drugs:
            utils.add_dict_counter(validDrugs, drug)
        for ind in indications.split(","):
            utils.add_dict_counter(validIndicates, ind)
        for se in ses.split(","):
            utils.add_dict_counter(validSes, se)
        if len(drugs) >= 2 and len(drugs) <= 20:
            drugs = sorted(drugs)
            for i in range(len(drugs)):
                for j in range(i + 1, len(drugs)):
                    d1, d2 = drugs[i], drugs[j]
                    pair = "%s,%s" % (d1, d2)
                    utils.add_dict_counter(validPairs, pair)


    cDrug = utils.sort_dict(validDrugs)
    cInd = utils.sort_dict(validIndicates)
    cSe = utils.sort_dict(validSes)
    cPair = utils.sort_dict(validPairs)
    print(len(cPair))
    writeSortedDictC(cDrug, "%s/%sADrugs.txt" % (OUT_DIR,PREF))
    writeSortedDictC(cInd, "%s/%sAInd.txt" % (OUT_DIR,PREF))
    writeSortedDictC(cSe, "%s/%sASe.txt" % (OUT_DIR,PREF))
    writeSortedDictC(cPair, "%s/%sPairs.txt" % (OUT_DIR,PREF))

if __name__ == "__main__":
    exportPair()