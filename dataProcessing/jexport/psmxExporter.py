from utils import utils
import params
import numpy as np


def loadDictName2Id(path):
    fin = open(path)
    d = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        name = parts[0]
        utils.get_update_dict_index(d, name)

    return d


def loadData():
    dDrug2Id = loadDictName2Id("%s/SJADERADrugs.txt" % params.JADER_OUT)
    dInd2Id = loadDictName2Id("%s/SJADERAInd.txt" % params.JADER_OUT)
    dSe2Id = loadDictName2Id("%s/SJADERASe.txt" % params.JADER_OUT)
    dDrugPair2Id = loadDictName2Id("%s/SJADERPairs.txt" % params.JADER_OUT)

    fin = open("%s/JADERInd.txt" % params.JADER_OUT)
    dList = []



def exportADrugPair(drugPair):
    pass
