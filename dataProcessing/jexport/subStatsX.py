import params
from utils import utils
import numpy as np
from utils import utils
import params
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Process, Value, Queue
from utils.banard import barnard_exact
import time


P_THRESHOLD = 0.05

_exact = fisher_exact

def exportSub():
    fin = open("%s/JADER.txt" % params.JADER_OUT)
    foutDict = dict()
    dlen2SeCount = dict()
    nA = 0
    print("Reading...")

    while True:
        line = fin.readline()
        if line == "":
            break
        nA += 1
        print("\r%s" % nA, end="")
        parts = line.strip().split("$")
        drugCmb = parts[0]
        ses = parts[1]
        drugs = drugCmb.split(",")
        nD = len(drugs)
        sortNames = ",".join(sorted(drugs))

        fO = utils.get_dict(foutDict, nD, -1)
        if fO == -1:
            fO = open("%s/SUB/%s" % (params.JADER_OUT, nD), "w")
            foutDict[nD] = fO
        fO.write("%s$%s\n" % (sortNames, ses))
        len2SeCount = utils.get_insert_key_dict(dlen2SeCount, nD, dict())
        sess = ses.split(",")
        for se in sess:
            utils.add_dict_counter(len2SeCount, se)

    for k, v in foutDict.items():
        v.close()

    d2 = dict()
    for k, v in dlen2SeCount.items():
        kvs = utils.sort_dict(v)
        ks = []
        for kv in kvs:
            kk, _ = kv
            ks.append(kk)
        d2[k] = ks
    utils.save_obj(d2, "%s/SUB/drugSize2CommonSEs" % params.JADER_OUT)




def exportSubG2():
    fin = open("%s/JADER.txt" % params.JADER_OUT)
    foutDict = dict()
    dlen2SeCount = dict()
    nA = 0
    print("Reading...")

    while True:
        line = fin.readline()
        if line == "":
            break
        nA += 1
        print("\r%s" % nA, end="")
        parts = line.strip().split("$")
        drugCmb = parts[0]
        ses = parts[1]
        drugs = drugCmb.split(",")
        nD = len(drugs)
        drugs = sorted(drugs)
        sortNames = ",".join(drugs)

        fO = utils.get_dict(foutDict, nD, -1)
        if fO == -1:
            fO = open("%s/SUB/G%s" % (params.JADER_OUT, nD), "w")
            foutDict[nD] = fO
        fO.write("%s$%s\n" % (sortNames, ses))
        if len(drugs) > 2 and len(drugs) <= 20:
            for i in range(len(drugs)):
                for j in range(i+1, len(drugs)):
                    d1 = drugs[i]
                    d2 = drugs[j]
                    pair = "%s,%s" % (d1, d2)
                    try:
                        f2 = foutDict[2]
                    except:
                        f2 = open("%s/SUB/G%s" % (params.JADER_OUT, 2), "w")
                        foutDict[2] = f2
                    f2.write("%s$%s\n" % (pair, ses))
        len2SeCount = utils.get_insert_key_dict(dlen2SeCount, nD, dict())
        sess = ses.split(",")
        for se in sess:
            utils.add_dict_counter(len2SeCount, se)

    for k, v in foutDict.items():
        v.close()

    d2 = dict()
    for k, v in dlen2SeCount.items():
        kvs = utils.sort_dict(v)
        ks = []
        for kv in kvs:
            kk, _ = kv
            ks.append(kk)
        d2[k] = ks
    utils.save_obj(d2, "%s/SUB/drugSize2CommonSEs" % params.JADER_OUT)


def producer(queue, arrs):
    for comAr in arrs:
        com, ar = comAr
        p = _exact(ar, 'greater')

        queue.put([com, ar[0, 0], p])


def consumer(queue, counter, counter2, fout=None, caches=None):
    while True:
        data = queue.get()
        if data is None:
            print("Receive terminate signal")
            with counter.get_lock():
                counter.value = 1
            if caches is not None:
                for line in caches:
                    fout.write("%s" % line)

            fout.flush()
            break
        com, cc, p = data
        with counter2.get_lock():
            counter2.value += 1

        # print(drugJader,">>", drugBankName)
        if fout is not None:
            ord, pv = p
            if pv <= P_THRESHOLD:
                if caches is None:
                    fout.write("%s\t%s\t%s\t%s\n" % (com, cc, ord, pv))
                else:
                    caches.append("%s\t%s\t%s\t%s\n" % (com, cc, ord, pv))


def exportBySE(seNames, pathIn, dirOut, pathInfo):
    fin = open(pathIn)
    dCombCount = dict()
    dCombSe = dict()
    dSe = dict()
    nA = 0
    print("Reading...")
    if not type(seNames) == set:
        seNames = set(seNames)
    print(seNames)
    while True:
        line = fin.readline()
        if line == "":
            break
        nA += 1
        parts = line.strip().split("$")
        drugCmb = parts[0]
        ses = parts[1]

        ses = set(ses.split(","))

        for se in seNames:
            dCombCountx = utils.get_insert_key_dict(dCombCount, se, dict())
            utils.add_dict_counter(dCombCountx, drugCmb)
            if se in ses:
                dComSEx = utils.get_insert_key_dict(dCombSe, se, dict())
                utils.add_dict_counter(dSe, se)
                utils.add_dict_counter(dComSEx, drugCmb)

    fin.close()
    print("Cal Contingency table...")
    dContigenTable = dict()

    for se in seNames:
        dCombCountx = dCombCount[se]
        dComSEx = utils.get_dict(dCombSe, se, dict())
        nSe = utils.get_dict(dSe, se, 0)
        if nSe == 0:
            continue
        for drugComb, nComb in dCombCountx.items():
            ar = np.zeros((2, 2))
            nCombSe = utils.get_dict(dComSEx, drugComb, 0)
            if nCombSe == 0:
                # print("SKIP")
                continue
            ar[0, 0] = nCombSe
            ar[1, 0] = nComb - nCombSe
            ar[0, 1] = nSe - nCombSe
            ar[1, 1] = nA - (nComb + nSe - nCombSe)
            nName = "%s_%s" % (drugComb, se)
            dContigenTable[nName] = ar

    producers = []
    consumers = []
    queue = Queue(params.K_FOLD)
    counter = Value('i', 0)
    counter2 = Value('i', 0)

    inputList = list(dContigenTable.items())
    nInputList = len(inputList)
    nDPerWorker = int(nInputList / params.N_DATA_WORKER)
    # assert 'g-csf' in allDrugNames
    for i in range(params.N_DATA_WORKER):
        startInd = i * nDPerWorker
        endInd = (i + 1) * nDPerWorker
        endInd = min(endInd, nInputList)
        if i == params.N_DATA_WORKER - 1:
            endInd = nInputList
        data = inputList[startInd:endInd]
        producers.append(Process(target=producer, args=(queue, data)))

    sname = "__".join(list(seNames))
    seNameString = "%s" % hash(sname)

    fFileNameMap = open(pathInfo, "a")
    fFileNameMap.write("%s\t%s\n" % (seNameString, sname))
    fFileNameMap.close()
    fout = open("%s/%s" % (dirOut, seNameString), "w")
    p = Process(target=consumer, args=(queue, counter, counter2, fout, []))
    p.daemon = True
    consumers.append(p)

    print("Start Producers...")
    for p in producers:
        p.start()
    print("Start Consumers...")
    for p in consumers:
        p.start()

    for p in producers:
        p.join()
    print("Finish Producers")

    queue.put(None)

    while True:
        if counter.value == 0:
            time.sleep(0.01)
            continue
        else:
            break
    fout.flush()
    fout.close()


def merge1():
    fin = open("%s/FSUBTEST/1/FileMap.txt" % params.JADER_OUT)
    fout = open("%s/FSUBTEST/1/1.txt" % params.JADER_OUT, "w")
    dDrug2Se = dict()

    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        hashFile = parts[0]
        f = open("%s/FSUBTEST/1/%s" % (params.JADER_OUT, hashFile))
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

def merge2():
    fin = open("%s/FSUBTEST/2/FileMap.txt" % params.JADER_OUT)
    fout = open("%s/FSUBTEST/2/2.txt" % params.JADER_OUT, "w")
    dDrug2Se = dict()

    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        hashFile = parts[0]
        f = open("%s/FSUBTEST/2/%s" % (params.JADER_OUT, hashFile))
        while True:
            l = f.readline()
            if l == "":
                break
            parts = l.strip().split("_")

            drug = parts[0]
            se = parts[1].split("\t")[0]
            ses = utils.get_insert_key_dict(dDrug2Se, drug, [])
            ses.append(se)
            # print(drug, ses)

        f.close()
    for k, v in dDrug2Se.items():
        fout.write("%s\t%s\n" % (k, ",".join(v)))
    fout.close()


def mergeG2():
    fin = open("%s/FSUBTEST/2/GFileMap.txt" % params.JADER_OUT)
    fout = open("%s/FSUBTEST/2/G2.txt" % params.JADER_OUT, "w")
    dDrug2Se = dict()

    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        hashFile = parts[0]
        f = open("%s/FSUBTEST/2/%s" % (params.JADER_OUT, hashFile))
        while True:
            l = f.readline()
            if l == "":
                break
            parts = l.strip().split("_")

            drug = parts[0]
            se = parts[1].split("\t")[0]
            ses = utils.get_insert_key_dict(dDrug2Se, drug, [])
            ses.append(se)
            # print(drug, ses)

        f.close()
    for k, v in dDrug2Se.items():
        fout.write("%s\t%s\n" % (k, ",".join(v)))
    fout.close()


def filter2():
    dir2 = "%s/FSUBTEST/2" % params.JADER_OUT
    utils.ensure_dir(dir2)

    dDrug1Se = dict()
    fin = open("%s/FSUBTEST/1/1.txt" % params.JADER_OUT)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        drug = parts[0]
        ses = set(parts[1].split(","))
        dDrug1Se[drug] = ses
    fin.close()
    fin = open("%s/SUB/2" % params.JADER_OUT)
    fout = open("%s/SUB/F2" % params.JADER_OUT, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("$")
        dDrug = parts[0].split(",")
        ses = parts[1].split(",")
        invalidSes = set()
        for drug in dDrug:
            sD = utils.get_dict(dDrug1Se, drug, set())
            for s in sD:
                invalidSes.add(s)
        validSes = []
        for se in ses:
            if se not in invalidSes:
                validSes.append(se)
        fout.write("%s$%s\n" % (parts[0], ",".join(validSes)))
    fout.close()



def filterg2():
    dir2 = "%s/FSUBTEST/2" % params.JADER_OUT
    utils.ensure_dir(dir2)

    dDrug1Se = dict()
    fin = open("%s/FSUBTEST/1/1.txt" % params.JADER_OUT)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        drug = parts[0]
        ses = set(parts[1].split(","))
        dDrug1Se[drug] = ses
    fin.close()
    fin = open("%s/SUB/G2" % params.JADER_OUT)
    fout = open("%s/SUB/GF2" % params.JADER_OUT, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("$")
        dDrug = parts[0].split(",")
        ses = parts[1].split(",")
        invalidSes = set()
        for drug in dDrug:
            sD = utils.get_dict(dDrug1Se, drug, set())
            for s in sD:
                invalidSes.add(s)
        validSes = []
        for se in ses:
            if se not in invalidSes:
                validSes.append(se)
        fout.write("%s$%s\n" % (parts[0], ",".join(validSes)))
    fout.close()


def filter3():
    dir3 = "%s/FSUBTEST/3" % params.JADER_OUT
    utils.ensure_dir(dir3)

    dDrug1Se = dict()
    dDrug2Se = dict()
    fin1 = open("%s/FSUBTEST/1/1.txt" % params.JADER_OUT)
    fin2 = open("%s/FSUBTEST/2/2.txt" % params.JADER_OUT)

    while True:
        line = fin1.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        drug = parts[0]
        ses = set(parts[1].split(","))
        dDrug1Se[drug] = ses
    fin1.close()

    while True:
        line = fin2.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        drug = parts[0]
        ses = set(parts[1].split(","))
        dDrug2Se[drug] = ses
    fin1.close()

    fin = open("%s/SUB/3" % params.JADER_OUT)
    fout = open("%s/SUB/F3" % params.JADER_OUT, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("$")
        dDrug = parts[0].split(",")
        ses = parts[1].split(",")
        invalidSes = set()
        for drug in dDrug:
            sD = utils.get_dict(dDrug1Se, drug, set())
            for s in sD:
                invalidSes.add(s)
        drugS = sorted(dDrug)
        drugPairs = []
        for i in range(len(drugS)):
            for j in range(i+1, len(drugS)):
                pair = "%s,%s" % (drugS[i], drugS[j])
                drugPairs.append(pair)
        for pair in drugPairs:
            sD = utils.get_dict(dDrug2Se, pair, set())
            for s in sD:
                invalidSes.add(s)
        validSes = []
        for se in ses:
            # if se not in invalidSes:
            validSes.append(se)
        fout.write("%s$%s\n" % (parts[0], ",".join(validSes)))
    fout.close()

def exportAllSes1():
    seList = utils.load_obj("%s/SUB/drugSize2CommonSEs" % params.JADER_OUT)[1]
    # seList = ['product dose omission']
    nSize = 50
    import os
    p = "%s/FSUBTEST/1/*" % params.JADER_OUT
    p = p.replace(" ", '\ ')

    cmd = "rm %s" % p
    try:
        os.system(cmd)
    except:
        pass
    pathInfo1 = "%s/FSUBTEST/1/FileMap.txt" % params.JADER_OUT
    pathIn1 = "%s/SUB/1" % (params.JADER_OUT)
    dirOut1 = "%s/FSUBTEST/1" % params.JADER_OUT
    fFileNameMap = open(pathInfo1, "w")
    fFileNameMap.close()

    nSeg = max(int(len(seList) / nSize), 1)

    for i in range(nSeg):
        start = i * nSize
        end = min((i + 1) * nSize, len(seList))
        exportBySE(seList[start:end], pathIn1, dirOut1, pathInfo1)


def exportAllSes2():
    seList = utils.load_obj("%s/SUB/drugSize2CommonSEs" % params.JADER_OUT)[2]
    # seList = ['product dose omission']
    nSize = 50
    import os
    p = "%s/FSUBTEST/2/*" % params.JADER_OUT
    p = p.replace(" ", "\ ")

    cmd = "rm %s" % p
    try:
        os.system(cmd)
    except:
        pass
    pathInfo1 = "%s/FSUBTEST/2/FileMap.txt" % params.JADER_OUT
    pathIn1 = "%s/SUB/F2" % (params.JADER_OUT)
    dirOut1 = "%s/FSUBTEST/2" % params.JADER_OUT
    fFileNameMap = open(pathInfo1, "w")
    fFileNameMap.close()

    nSeg = max(int(len(seList) / nSize), 1)

    for i in range(nSeg):
        start = i * nSize
        end = min((i + 1) * nSize, len(seList))
        exportBySE(seList[start:end], pathIn1, dirOut1, pathInfo1)


def exportAllSesG2():
    seList = utils.load_obj("%s/SUB/drugSize2CommonSEs" % params.JADER_OUT)[2]
    # seList = ['product dose omission']
    nSize = 50
    import os
    p = "%s/FSUBTEST/2/*" % params.JADER_OUT
    p = p.replace(" ", "\ ")

    cmd = "rm %s" % p
    try:
        os.system(cmd)
    except:
        pass
    pathInfo1 = "%s/FSUBTEST/2/GFileMap.txt" % params.JADER_OUT
    pathIn1 = "%s/SUB/GF2" % (params.JADER_OUT)
    dirOut1 = "%s/FSUBTEST/2" % params.JADER_OUT
    fFileNameMap = open(pathInfo1, "w")
    fFileNameMap.close()

    nSeg = max(int(len(seList) / nSize), 1)

    for i in range(nSeg):
        start = i * nSize
        end = min((i + 1) * nSize, len(seList))
        exportBySE(seList[start:end], pathIn1, dirOut1, pathInfo1)
def exportAllSes3():
    seList = utils.load_obj("%s/SUB/drugSize2CommonSEs" % params.JADER_OUT)[2]
    # seList = ['product dose omission']
    nSize = 50
    import os
    p = "%s/FSUBTEST/3/*" % params.JADER_OUT
    p = p.replace(" ", "\ ")

    cmd = "rm %s" % p
    try:
        os.system(cmd)
    except:
        pass
    pathInfo1 = "%s/FSUBTEST/3/FileMap.txt" % params.JADER_OUT
    pathIn1 = "%s/SUB/F3" % (params.JADER_OUT)
    dirOut1 = "%s/FSUBTEST/3" % params.JADER_OUT
    fFileNameMap = open(pathInfo1, "w")
    fFileNameMap.close()

    nSeg = max(int(len(seList) / nSize), 1)

    for i in range(nSeg):
        start = i * nSize
        end = min((i + 1) * nSize, len(seList))
        exportBySE(seList[start:end], pathIn1, dirOut1, pathInfo1)
def ensureDIR():
    utils.ensure_dir("%s/FSUBTEST" % params.JADER_OUT)
    utils.ensure_dir("%s/FSUBTEST/1" % params.JADER_OUT)
    utils.ensure_dir("%s/FSUBTEST/2" % params.JADER_OUT)

    utils.ensure_dir("%s/SUB" % params.JADER_OUT)

def exportPolyJADER():
    fin = open("%s/Data/DrugBank/DrugBankNames.txt" % params.C_DIR)
    dName2Inchi = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("||")
        drugName = parts[0]
        inchi = parts[3]
        dName2Inchi[drugName] = inchi
    fin.close()

    fin = open("%s/FSUBTEST/2/2.txt" % params.JADER_OUT)
    fout = open("%s/FSUBTEST/2/JADERDDI.txt" % params.JADER_OUT, "w")


    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        d1, d2 = parts[0].split(",")
        i1, i2 = utils.get_dict(dName2Inchi, d1, -1), utils.get_dict(dName2Inchi, d2, -1)
        if i1 == -1 or i2 == -1:
            continue
        if len(i1) < 2 or len(i2) < 2:
            continue
        fout.write("%s|%s|%s|%s|%s\n" % (d1, d2, i1, i2, parts[1]))
    fout.close()
    fin.close()


def exportPolyGJADER():
    fin = open("%s/Data/DrugBank/DrugBankNames.txt" % params.C_DIR)
    dName2Inchi = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("||")
        drugName = parts[0]
        inchi = parts[3]
        dName2Inchi[drugName] = inchi
    fin.close()

    fin = open("%s/FSUBTEST/2/G2.txt" % params.JADER_OUT)
    fout = open("%s/FSUBTEST/2/GJADERDDI.txt" % params.JADER_OUT, "w")


    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        d1, d2 = parts[0].split(",")
        i1, i2 = utils.get_dict(dName2Inchi, d1, -1), utils.get_dict(dName2Inchi, d2, -1)
        if i1 == -1 or i2 == -1:
            continue
        if len(i1) < 2 or len(i2) < 2:
            continue
        fout.write("%s|%s|%s|%s|%s\n" % (d1, d2, i1, i2, parts[1]))
    fout.close()
    fin.close()
if __name__ == "__main__":
    ensureDIR()
    # exportSub()
    # exportSubG2()
    # exportAllSes1()
    # merge1()
    # filter2()
    # filterg2()
    exportAllSesG2()
    # merge2()
    mergeG2()

    # filter3()
    # exportAllSes3()
    # exportPolyJADER()
    exportPolyGJADER()

    pass