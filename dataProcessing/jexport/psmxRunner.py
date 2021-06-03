import math
import random

from utils import utils
import params
import numpy as np
import time
from multiprocessing import Process, Value, Queue
from dataProcessing.jexport.psmxExporter import loadDictName2Id
import statsmodels.api as sm

from statsmodels.genmod.generalized_linear_model import GLM

OUT_DIR = params.JADER_OUT
PREF = "SJADER"


def calPerson(matX, Y):
    Y = Y[:, np.newaxis]
    xy = matX * Y
    Exy = np.mean(xy, axis=0)
    Ex = np.mean(matX, axis=0)
    Ey = np.mean(Y, axis=0)
    dxy = Exy - Ex * Ey
    x2 = matX * matX
    Ex2 = np.mean(x2, axis=0)
    dx = np.sqrt(Ex2 - Ex * Ex)
    y2 = Y * Y
    Ey2 = np.mean(y2, axis=0)
    dy = np.sqrt(Ey2 - Ey * Ey)
    p = dxy / (dx * dy + 1e-10)
    return p


def tP():
    matX = np.random.uniform(0, 1, (1000, 2000))
    matY = np.random.uniform(0, 1, 1000)

    p = calPerson(matX, matY)
    print(p)


def matching(score1, score2):
    # n1 = score1.shape[-1]
    # n2 = score2.shape[-1]
    min1 = np.min(score1)
    max1 = np.max(score1)
    nBin = 20
    binSize = (max1 - min1) / nBin + 1e-10

    bin1s = [[] for _ in range(nBin)]
    bin2s = [[] for _ in range(nBin)]

    sbin1s = np.floor((score1 - min1) / binSize)
    sbin2s = np.floor((score2 - min1) / binSize)

    for i, v in enumerate(sbin1s):
        v = int(v)
        if v == nBin:
            v = nBin - 1
        if v >= nBin or v < 0:
            print(v)
        bin1s[v].append(i)
    for i, v in enumerate(sbin2s):
        v = int(v)
        if v == nBin:
            v = nBin - 1
        if nBin > v >= 0:
            bin2s[v].append(i)
    r1 = []
    r2 = []
    for i in range(nBin):
        bin1 = bin1s[i]
        bin2 = bin2s[i]

        nSample = 10 * len(bin1)
        if nSample > 0 and len(bin2) > 0:
            rbin2 = np.random.choice(bin2, nSample)
            for j in bin1:
                r1.append(j)
            for j in rbin2:
                r2.append(j)

    return r1, r2


def getSubList(list, ids):
    subList = []
    for id in ids:
        subList.append(list[id])
    return subList


def producer(queue, datas):
    oRs, drugPairList, dDrug2Id, dInd2Id, dList = datas
    nD = len(dDrug2Id)
    nInd = len(dInd2Id)
    fSize = nD + nInd
    for oR in oRs:
        pId, rExposeIds, rNonExposeIds = oR
        rExpose = getSubList(dList, rExposeIds)
        rNoneExpose = getSubList(dList, rNonExposeIds)
        nExpose = len(rExpose)
        nNonExpose = len(rNoneExpose)
        drugPair = drugPairList[pId]
        d1, d2 = drugPair.split(",")
        d1, d2 = utils.get_dict(dDrug2Id, d1, -1), utils.get_dict(dDrug2Id, d2, -1)
        matX1 = np.zeros((nExpose, fSize))
        Y1 = np.ones(nExpose)
        matX2 = np.zeros((nNonExpose, fSize))
        Y2 = np.zeros(nNonExpose)
        for i, v in enumerate(rExpose):
            drugIds, indIds, _ = v
            matX1[i, drugIds] = 1
            matX1[i, d1] = 0
            matX1[i, d2] = 0
            matX1[i, indIds] = 1
        for i, v in enumerate(rNoneExpose):
            drugIds, indIds, _ = v
            matX2[i, drugIds] = 1
            matX2[i, indIds] = 1

        matX = np.vstack((matX1, matX2))
        Y = np.concatenate((Y1, Y2))
        p = calPerson(matX, Y)
        args = np.argsort(p)[::-1][:200]
        matX = matX[:, args]

        glm = GLM(Y, matX, family=sm.families.Binomial())
        res = glm.fit()
        scores = res.predict(matX)
        scores1 = scores[:nExpose]
        scores2 = scores[nExpose:]
        r1, r2 = matching(scores1, scores2)
        r1 = getSubList(rExposeIds, r1)
        r2 = getSubList(rNonExposeIds, r2)
        del glm
        del res
        queue.put([pId, r1, r2])
        # print("\rPut...%s" % queue.size(), end="")


def consumer(queue, counter, counter2, fout=None, caches=None, maxCache=20):
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
        i, expose, nonExpose = data
        expose = ["%s" % i for i in expose]
        nonExpose = ["%s" % i for i in nonExpose]
        with counter2.get_lock():
            counter2.value += 1
            # if counter2.value % 1000 == 0:
            print("\r Processed %s" % counter2.value, end="")

        if fout is not None:

            if caches is None:
                fout.write("%s\t%s\t%s\n" % (i, ",".join(expose), ",".join(nonExpose)))
            else:
                caches.append("%s\t%s\t%s\n" % (i, ",".join(expose), ",".join(nonExpose)))
                if len(caches) >= maxCache:
                    for line in caches:
                        fout.write("%s" % line)
                    fout.flush()
                    print("Flush...", end="")
                    caches.clear()


def loadRawExpose():
    fin = open("%s/%s" % (OUT_DIR, "rawExpose"))
    ors = []
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        pid = int(parts[0])
        exposeIds = []
        nonExposeIds = []
        try:
            for ep in parts[1].split(","):
                exposeIds.append(int(ep))
            for ne in parts[2].split(","):
                nonExposeIds.append(int(ne))
        except:
            print(line)
        ors.append([pid, exposeIds, nonExposeIds])
    fin.close()
    return ors


def pExport():
    producers = []
    consumers = []
    queue = Queue(params.K_FOLD)
    counter = Value('i', 0)
    counter2 = Value('i', 0)

    dList = utils.load_obj("%s/DataDump.o" % OUT_DIR)
    dDrugPair2Id, drugPairList = loadDictName2Id("%s/%sPairs.txt" % (OUT_DIR, PREF), nMax=-1, min=1)
    dDrug2Id, _ = loadDictName2Id("%s/%sADrugs.txt" % (OUT_DIR, PREF))
    dInd2Id, _ = loadDictName2Id("%s/%sAInd.txt" % (OUT_DIR, PREF))

    inputList = loadRawExpose()
    print("inputLen: ", len(inputList))
    nInputList = len(drugPairList)
    nDPerWorker = int(nInputList / params.N_DATA_WORKER)
    # assert 'g-csf' in allDrugNames
    for i in range(params.N_DATA_WORKER):
        startInd = i * nDPerWorker
        endInd = (i + 1) * nDPerWorker
        endInd = min(endInd, nInputList)
        if i == params.N_DATA_WORKER - 1:
            endInd = nInputList
        data = inputList[startInd:endInd], drugPairList, dDrug2Id, dInd2Id, dList
        producers.append(Process(target=producer, args=(queue, data)))

    fout = open("%s/%s" % (OUT_DIR, "psmExpose"), "w")
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


if __name__ == "__main__":
    pExport()
