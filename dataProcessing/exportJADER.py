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

if __name__ == "__main__":
    exportJADER()