from utils import utils
import params
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Process, Value, Queue
import time

P_THRESHOLD = 0.05


def producer(queue, arrs):
    for comAr in arrs:
        com, ar = comAr
        p = fisher_exact(ar, 'greater')
        queue.put([com, p])


def consumer(queue, counter, counter2, fout=None):
    while True:
        data = queue.get()
        if data is None:
            print("Receive terminate signal")
            with counter.get_lock():
                counter.value = 1
            fout.flush()
            break
        com, p = data
        with counter2.get_lock():
            counter2.value += 1

        # print(drugJader,">>", drugBankName)
        if fout is not None:
            ord, pv = p
            if pv <= P_THRESHOLD:
                fout.write("%s\t%s\t%s\n" % (com, ord, pv))


def exportBySE(seNames):
    fin = open("%s/FDrug2SeList_20050.txt" % params.FADER_OUT)
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
    print("Cal Contigence table...")
    dContigenTable = dict()

    for se in seNames:
        dCombCountx = dCombCount[se]
        dComSEx = utils.get_dict(dCombSe, se, dict())
        nSe = dSe[se]
        for drugComb, nComb in dCombCountx.items():
            ar = np.zeros((2, 2))
            nCombSe = utils.get_dict(dComSEx, drugComb, 0)

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

    fFileNameMap = open("%s/FTest/FileMap.txt" % params.FADER_OUT, "a")
    fFileNameMap.write("%s\t%s\n" % (seNameString, sname))
    fFileNameMap.close()

    fout = open("%s/FTest/%s" % (params.FADER_OUT, seNameString), "w")
    p = Process(target=consumer, args=(queue, counter, counter2, fout))
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


def loadValidSEs(nT=0):
    lines = open("%s/ValidSes.txt" % params.FADER_OUT).readlines()
    validSes = list()
    for line in lines:
        parts = line.strip().split("\t")
        se, v = parts[0], int(parts[1])
        if v > nT:
            validSes.append(se)
        else:
            break

    return validSes


def exportAllSes():
    validSEs = loadValidSEs(nT=100)
    seList = list(validSEs)

    # seList = ['product dose omission']
    nSize = 10
    import os
    os.system("rm '%s/FTest/*'" % params.FADER_OUT)
    fFileNameMap = open("%s/FTest/FileMap.txt" % params.FADER_OUT, "w")
    fFileNameMap.close()

    nSeg = max(int(len(seList) / nSize), 1)

    for i in range(nSeg):
        start = i * nSize
        end = min((i+1) * nSize, len(seList))

        exportBySE(seList[start:end])


if __name__ == "__main__":
    # seNames = "dizziness"
    # exportBySE(seNames)
    exportAllSes()
