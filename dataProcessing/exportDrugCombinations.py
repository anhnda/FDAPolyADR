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


if __name__ == "__main__":
    # nSize = len(loadValidDrugMap())

    # getAllDrugSet()
    # stats1(nSize)
    # stats2(nSize)
    exportValidSEs()
    getAllDrugSEMap()
    # exportSeCount()
    pass
