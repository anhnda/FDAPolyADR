from utils import utils

import params
import glob
import os


def stripDrugNameO(dName):
    if len(dName) > 0:
        if dName[-1] == "." or dName[-1] == ",":
            dName = dName[:-1]
    return dName


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


def getDrugSet(path, dDrugSet, dDrugCombSet, dMap=dict()):
    fin = open(path, encoding="utf8", errors='ignore')
    fin.readline()

    currentId = -1
    currentDrugSet = set()
    print("Loading: ...", path)
    skipCase = False
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().lower()
        parts = line.split("$")
        caseId = parts[1]
        drugName = parts[4]
        drugName = stripDrugNameO(drugName)
        if len(drugName) == 0:
            skipCase = True
            currentId = caseId
            currentDrugSet = set()
            continue
        if len(dMap) == 0:
            utils.add_dict_counter(dDrugSet, drugName)
        else:
            drugName = utils.get_dict(dMap, drugName, -1)
            if drugName == -1:
                skipCase = True

        if caseId != currentId:
            if currentId != -1 and not skipCase:
                utils.add_dict_counter(dDrugCombSet, tuple(currentDrugSet), 1)
                for dName in currentDrugSet:
                    utils.add_dict_counter(dDrugSet, dName)
            currentId = caseId
            currentDrugSet = set()
            if drugName != -1:
                skipCase = False

        if not skipCase:
            if type(drugName) == int:
                print(currentId, caseId)
                print(line)
                exit(-1)
            currentDrugSet.add(drugName)
    fin.close()


def getSideEffectSet(path, seCounter, dValidSes=set()):
    fin = open(path, encoding="utf8", errors='ignore')
    fin.readline()

    currentId = -1
    currentSESet = set()
    print("Loading: ...", path)
    dCase2Se = dict()
    skipCase = False
    assert 'medication error' not in dValidSes
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().lower()
        parts = line.split("$")
        caseId = parts[1]
        seName = parts[2]

        assert len(seName) > 0

        if caseId != currentId:
            if currentId != -1:
                dCase2Se[currentId] = currentSESet
                for se in currentSESet:
                    utils.add_dict_counter(seCounter, se)
            currentId = caseId
            currentSESet = set()

        if seName in dValidSes:
            currentSESet.add(seName)
    fin.close()

    return dCase2Se


def getSideEffectSet2(path, seCounter, dValidSes=set()):
    fin = open(path, encoding="utf8", errors='ignore')
    fin.readline()

    currentId = -1
    currentSESet = set()
    print("Loading: ...", path)
    dCase2Se = dict()
    skipCase = False
    assert 'medication error' not in dValidSes
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().lower()
        parts = line.split("$")
        caseId = parts[1]
        seName = parts[2]

        assert len(seName) > 0

        if len(dValidSes) > 0:
            if seName not in dValidSes:
                skipCase = True

        if caseId != currentId:
            if currentId != -1 and not skipCase:
                dCase2Se[currentId] = currentSESet
                for se in currentSESet:
                    if se == 'medication error':
                        print("???")
                        exit(-1)

                    utils.add_dict_counter(seCounter, se)
            currentId = caseId
            currentSESet = set()
            if seName in dValidSes:
                skipCase = False
        if not skipCase:
            currentSESet.add(seName)
    fin.close()

    return dCase2Se
def getDrugSEMappingFile(path, fout, dMap, dCaseSE):
    fin = open(path, encoding="utf8", errors='ignore')
    fin.readline()

    currentId = -1
    currentDrugSet = set()
    print("Loading: ...", path)
    skipCase = False
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().lower()
        parts = line.split("$")
        caseId = parts[1]
        drugName = parts[4]
        drugName = stripDrugNameO(drugName)
        if len(drugName) == 0:
            skipCase = True
            currentId = caseId
            currentDrugSet = set()
            continue

        drugName = utils.get_dict(dMap, drugName, -1)

        if drugName == -1:
            skipCase = True

        if caseId != currentId:
            if currentId != -1 and not skipCase:
                seSet = utils.get_dict(dCaseSE, currentId, set())
                if len(seSet) > 0:
                    fout.write("%s$%s\n" % (",".join(list(currentDrugSet)), ",".join(list(seSet))))

            currentId = caseId
            currentDrugSet = set()
            if drugName != -1:
                skipCase = False

        if not skipCase:
            if type(drugName) == int:
                print(currentId, caseId)
                print(line)
                exit(-1)
            currentDrugSet.add(drugName)
    fin.close()


def getDrugFile(dir):
    suff = dir[-4:]
    path = "%s/ascii/DRUG%s.txt" % (dir, suff.upper())
    if not os.path.isfile(path):
        path = "%s/ASCII/DRUG%s.txt" % (dir, suff.upper())
    if not os.path.isfile(path):
        path = "%s/ascii/drug%s.txt" % (dir, suff)
    return path


def getSEFile(dir):
    suff = dir[-4:]
    path = "%s/ascii/REAC%s.txt" % (dir, suff.upper())
    if not os.path.isfile(path):
        path = "%s/ASCII/REAC%s.txt" % (dir, suff.upper())
    if not os.path.isfile(path):
        path = "%s/ascii/reac%s.txt" % (dir, suff)
    return path


def getAllDrugSet():
    dirs = glob.glob("%s/*" % params.FADER_DIR)
    drugNameSet = dict()
    drugCombSet = dict()
    dMap = loadValidDrugMap()
    nSize = len(dMap)
    print("DMAP SIZE: ", nSize)
    for dir in dirs:
        path = getDrugFile(dir)
        assert os.path.isfile(path)
        getDrugSet(path, drugNameSet, drugCombSet, dMap)

    print("Saving...")
    utils.save_obj(drugNameSet, "%s/FDrugNameCount_%s" % (params.FADER_OUT, nSize))
    utils.save_obj(drugCombSet, "%s/FDrugCombCount_%s" % (params.FADER_OUT, nSize))


def getAllDrugSEMap():
    dirs = glob.glob("%s/*" % params.FADER_DIR)
    validSes = loadValidSEs()
    nSE = len(validSes)
    fout = open("%s/FDrug2SeList_%s.txt" % (params.FADER_OUT, nSE), "w")

    dMap = loadValidDrugMap()
    assert len(dMap) > 0
    nSize = len(dMap)
    print("DMAP SIZE: ", nSize)
    seCount = dict()
    for dir in dirs:
        pathDrug = getDrugFile(dir)
        assert os.path.isfile(pathDrug)
        pathSE = getSEFile(dir)
        caseSEMap = getSideEffectSet(pathSE, seCount, validSes)
        getDrugSEMappingFile(pathDrug, fout, dMap, caseSEMap)

    print("Saving...")
    utils.save_obj(seCount, "%s/FSECount_%s_%s" % (params.FADER_OUT, nSize, nSE))

    fout.close()


def stats1(nSize=0):
    print("Loading...")
    drugComb = utils.load_obj("%s/FDrugNameCount_%s" % (params.FADER_OUT, nSize))
    print("Sorting..")
    kvs = utils.sort_dict(drugComb)

    fout = open("%s/FDrugNamesSort_%s" % (params.FADER_OUT, nSize), "w")
    print("Saving...")
    for kv in kvs:
        k, v = kv
        if len(k) <= 1:
            continue

        fout.write("%s$%s\n" % (k, v))

    fout.close()


def stats2(nSize=0):
    print("Loading...")
    drugComb = utils.load_obj("%s/FDrugCombCount_%s" % (params.FADER_OUT, nSize))

    print("Sorting..")
    kvs = utils.sort_dict(drugComb)

    fout = open("%s/FDrugCombSort_%s" % (params.FADER_OUT, nSize), "w")

    print("Saving...")
    cc = 0
    for kv in kvs:
        k, v = kv
        # print(k, v)
        cc += v
        fout.write("%s$%s\n" % (",".join(k), v))
    fout.close()
    print("Total: %s cases" % cc)
    from plotLib import plotCul2

    plotCul2(kvs[::-1], 200, 1, "SelectedCombDrugCutOff", xLabel="ThreshHold: Freq >=", yLabel="Number of Combs")


def exportSeCount(nSize=9210):
    d = utils.load_obj("%s/FSECount_%s_0" % (params.FADER_OUT, nSize))
    kvs = utils.sort_dict(d)
    fout = open("%s/FSECountSorted_%s_0" % (params.FADER_OUT, nSize), "w")
    for kv in kvs:
        k, v = kv
        fout.write("%s\t%s\n" % (k, v))
    fout.close()


def exportValidSEs(nSize=9210):
    def loadException(path="%s/InValidSEs.txt" % params.FADER_OUT):
        lines = open(path).readlines()
        invalidSes = set()
        invalidTokens = list()
        for line in lines:
            line = line.strip()
            if line[0] == '#':
                invalidTokens.append(line[1:])
            else:
                invalidSes.add(line)
        return invalidSes, invalidTokens

    invalidSes, invalidTokens = loadException()

    fout = open("%s/ValidSes.txt" % params.FADER_OUT, "w")
    d = utils.load_obj("%s/FSECount_%s_0" % (params.FADER_OUT, nSize))
    kvs = utils.sort_dict(d)

    for kv in kvs:
        k, v = kv
        if k in invalidSes:
            continue
        isInvalid = False

        for token in invalidTokens:
            if k.__contains__(token):
                isInvalid = True
                break
        if isInvalid:
            continue
        fout.write("%s\t%s\n" % (k, v))
    fout.close()


def loadValidSEs():
    fin = open("%s/ValidSes.txt" % params.FADER_OUT)
    try:
        lines = fin.readlines()
    except:
        print("All SEs allowed")
        return set()
    validSes = set([line.strip().split("\t")[0] for line in lines])
    fin.close()
    print("NSE: ", len(validSes))
    return validSes

def exportDrugCom2Side():
    fin = open("%s/JADER.txt" % params.JADER_OUT)
    fout = open("%s/JADER2AllSeList.txt" % params.JADER_OUT, "w")
    dDrugComb2Se = dict()
    dDrugCombCount = dict()
    dDrugCom2Lenght = dict()
    drugCont = dict()
    seCount = dict()
    cc = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        cc += 1
        line = line.strip()
        parts = line.split("$")
        drugCom = parts[0]
        dDrugCom2Lenght[drugCom] = len(drugCom.split(","))

        ses = parts[1].split(",")
        utils.add_dict_counter(dDrugCombCount, drugCom, 1)
        for drug in drugCom.split(","):
            utils.add_dict_counter(drugCont, drug, 1)
        sesComb = utils.get_insert_key_dict(dDrugComb2Se, drugCom, dict())
        for se in ses:
            utils.add_dict_counter(sesComb, se, 1)
            utils.add_dict_counter(seCount, se)

    kvs = utils.sort_dict(dDrugCombCount)
    for kv in kvs:
        k, v = kv
        seCountKv = utils.sort_dict(dDrugComb2Se[k])
        sString = []
        for seCountx in seCountKv:
            se,count = seCountx
            sString.append("%s:%s"% (se, count))

        fout.write("%s:%s$%s$%s\n" % (k, v, len(sString), ",".join(sString)))
    fout.close()
    utils.save_obj(seCount, "%s/JADERSeCountFX" % params.JADER_OUT)
    utils.save_obj(dDrugCom2Lenght, "%s/DrugCombLength" % params.JADER_OUT)
    print(len(drugCont), len(seCount))
def plotDrugCombCount():
    fin = open("%s/JADER2AllSeList.txt" % params.JADER_OUT)
    dLength2NReports = dict()
    kv = []
    vs = []
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().split("$")
        parts = line[0].split(":")
        c = int(parts[1])
        drugCombLenght = len(parts[0].split(","))
        utils.add_dict_counter(dLength2NReports, drugCombLenght, c)
        vs.append(c)
        kv.append([parts[0], c])


    import matplotlib.pyplot as plt
    import numpy as np
    maxX = max(dLength2NReports.keys())
    x = [ i for i in range(1, maxX+1)]
    y = np.zeros(maxX)
    for k, v in dLength2NReports.items():
        y[k-1] = v

    plt.scatter(x, y)
    plt.xlabel("DrugComb length")
    plt.ylabel("Num Reports")
    plt.tight_layout()
    plt.savefig("%s/JADERReportsOnDrugLength.png" % params.FIG_DIR)

    from dataProcessing.plotLib import plotCul2, plotCul, plotHistD
    print(len(kv), kv[-1])
    print(kv[0])
    print(max(vs), min(vs))
    plotCul(kv[::-1], 10, 1, "JADERDrugCombFreq", xLabel="Threshold of DrugComb Frequency", yLabel="Num DrugComb")
    plotCul2(kv[::-1], 10, 1, "JADERDrugCombReports",  xLabel="Threshold of DrugComb Frequency", yLabel="Num Reports" )

    # plotHistD(vs, 100, "HistDrugCombFrequency")

def plotSeCount():
    seCount = utils.load_obj( "%s/JADERSeCountFX" % params.JADER_OUT)
    kvs = utils.sort_dict(seCount)


    from dataProcessing.plotLib import plotCul2, plotCul, plotHistD

    plotCul(kvs[::-1], 50, 1, "JADERSEFreq", xLabel="Thresholds of SE Frequency", yLabel="Num. SEs")
    # plotCul2(kvs[::-1], 50, 1, "SeReports")

def plotDrugCombLength():
    dLength = utils.load_obj("%s/DrugCombLength" % params.JADER_OUT)

    kvs = utils.sort_dict(dLength)
    dCount = dict()
    for kv in kvs:
        _, v = kv
        utils.add_dict_counter(dCount, v)

    maxLength = max(dCount.keys())
    x = [i for i in range(1, maxLength+1)]
    import numpy as np

    y = np.zeros(maxLength)
    for k, v in dCount.items():
        y[k-1] = v

    import matplotlib.pyplot as plt
    plt.scatter(x,y)
    plt.xlabel("DrugComb length")
    plt.ylabel("Num DrugComb")
    plt.tight_layout()
    plt.savefig("%s/%s.png" % (params.FIG_DIR, "JADERDrugLength"))


def plotDrugLength2NSEs():
    import numpy as np

    dDrugLength2NSes = dict()
    fin = open("%s/JADER2AllSeList.txt" % params.JADER_OUT)


    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("$")
        drugCombs = parts[0].split(":")[0]
        nDrug = len(drugCombs.split(","))
        nSe = int(parts[1])
        seLengths = utils.get_insert_key_dict(dDrugLength2NSes, nDrug, [])
        seLengths.append(nSe)


    x = dDrugLength2NSes.keys()
    xmax = max(x)
    x  = [i for i in range(1, xmax + 1)]
    y = np.zeros(xmax)
    for k, v in dDrugLength2NSes.items():
        # avg = sum(v) / len(v)
        y[k-1] = np.median(v)
    import matplotlib.pyplot as plt
    plt.scatter(x,y)
    plt.xlabel("DrugComb length")
    plt.ylabel("Median SEs")
    plt.tight_layout()
    plt.savefig("%s/%s.png" % (params.FIG_DIR, "JADERAvgSEDrugCombLength"))
if __name__ == "__main__":


    # exportDrugCom2Side()
    # plotDrugCombCount()
    # plotSeCount()
    # plotDrugCombLength()
    plotDrugLength2NSEs()
    pass
