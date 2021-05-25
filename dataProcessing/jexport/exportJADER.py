from utils import utils

import params
import glob
import os




def loadValidDrugMap():
    path = "%s/finalMap/DrugFMap.txt" % params.OUTPUT_DIR
    d = dict()
    if os.path.isfile(path):
        lines = open(path).readlines()
        for line in lines:
            line = line.strip()
            parts = line.split("||")
            d[parts[0]] = parts[1]
    return d

def exportJADER():
    d = loadValidDrugMap()
    fin = open("%s/FinalJADER.txt" % params.JADER_OUT)
    fout = open("%s/JADER.txt" % params.JADER_OUT, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        drugs = parts[1].split(",")
        isValid = True
        for drug in drugs:
            if drug not in d:
                isValid = False
                break
        if isValid:
            fout.write("%s$%s\n" % (parts[1], parts[-1]))
    fout.close()


def exportJADER2():
    d = loadValidDrugMap()
    fin = open("%s/FinalJADER.txt" % params.JADER_OUT)
    fout = open("%s/JADERInd.txt" % params.JADER_OUT, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        drugs = parts[1].split(",")
        # indicates = parts[2].split(",")
        isValid = True
        for drug in drugs:
            if drug not in d:
                isValid = False
                break
        if isValid:
            fout.write("%s$%s$%s$%s\n" % (parts[1], parts[2], parts[3], parts[-1]))
    fout.close()

def exportJADERPair():
    fin = open("%s/JADERInd.txt" % params.JADER_OUT)
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
        drugComb = parts[0]
        indications = parts[2]
        ses = parts[3]
        drugs = drugComb.split(",")
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
    writeSortedDictC(cDrug, "%s/SJADERADrugs.txt" % params.JADER_OUT)
    writeSortedDictC(cInd, "%s/SJADERAInd.txt" % params.JADER_OUT)
    writeSortedDictC(cSe, "%s/SJADERASe.txt" % params.JADER_OUT)
    writeSortedDictC(cPair, "%s/SJADERPairs.txt" % params.JADER_OUT)

def writeSortedDictC(dl, p):
     f = open(p, "w")
     for kv in dl:
         k, v = kv
         f.write("%s\t%s\n" % (k, v))
     f.close()


def checkX():
    fin = open("%s/SJADERPairs.txt" % params.JADER_OUT)
    inpPairs = set()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        dPair = parts[0]
        inpPairs.add(dPair)


    fin = open("%s/GJADERDDI.txt" % params.JADER_OUT)
    outPairs = set()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("|")
        pair = "%s,%s" % (parts[0], parts[1])
        outPairs.add(pair)
    cc1 = 0
    for pair1 in inpPairs:
        if pair1 in outPairs:
            cc1 += 1
    cc2 = 0
    for pair2 in outPairs:
        if pair2 in inpPairs:
            cc2 += 1
    print("NInp, NOut, NInpOout, NOutInp", cc1, cc2, len(inpPairs), len(outPairs))
if __name__ == "__main__":
    # exportJADER()
    # exportJADER2()
    # exportJADERPair()
    checkX()