import math
import random

from utils import utils
import params
import numpy as np
import time
from multiprocessing import Process, Value, Queue

NEG_SAMPLE_RATE = 0.01
MAX_NEG = 5000

OUT_DIR = params.JADER_OUT
PREF = "SJADER"

def loadDictName2Id(path, nMax=1000, min=1):
    fin = open(path)
    d = dict()
    ic = 0
    keyList = []
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        name = parts[0]
        cout = float(parts[1])
        if cout <= min:
            break
        utils.get_update_dict_index(d, name)
        ic += 1
        keyList.append(name)
        if ic == nMax:
            break
    fin.close()
    return d, keyList


def exportOData():
    dDrug2Id, _ = loadDictName2Id("%s/%sADrugs.txt" % (OUT_DIR, PREF))
    dInd2Id, _ = loadDictName2Id("%s/%sAInd.txt" % (OUT_DIR, PREF))
    dSe2Id, _ = loadDictName2Id("%s/%sASe.txt" % (OUT_DIR, PREF))
    fin = open("%s/JADERInd.txt" % OUT_DIR)

    dList = []
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("$")
        drugs = parts[0].split(",")
        inds = parts[2].split(",")
        ses = parts[-1].split(",")
        drugIds = []
        indIds = []
        seIds = []
        if len(drugs) > 20:
            continue
        for drug in drugs:
            drugId = utils.get_dict(dDrug2Id, drug, -1)
            # print(drug, drugId)
            if drugId != -1:
                drugIds.append(drugId)
        for ind in inds:
            indId = utils.get_dict(dInd2Id, ind, -1)
            if indId != -1:
                indIds.append(indId)
        for se in ses:
            seId = utils.get_dict(dSe2Id, se, -1)
            if seId != -1:
                seIds.append(seId)
        # print(drugIds, indIds, seIds)
        dList.append([drugIds, indIds, seIds])

    utils.save_obj(dList, "%s/DataDump.o" % OUT_DIR)


def producer(queue, datas):
    js, drugPairList, dDrug2Id, dList = datas
    nReport = len(dList)
    arReport = [i for i in range(nReport)]
    for j in js:
        drugPair = drugPairList[j]
        d1, d2 = drugPair.split(",")
        d1, d2 = utils.get_dict(dDrug2Id, d1, -1), utils.get_dict(dDrug2Id, d2, -1)
        exposeD12 = []
        nonExposeD12 = []
        random.shuffle(arReport)
        for ir in arReport:
            record = dList[ir]
            drugIds, _, _ = record
            # print(drugIds)
            if d1 in drugIds and d2 in drugIds:
                exposeD12.append("%s" % ir)
            else:
                nonExposeD12.append("%s" % ir)
        # print(d1, d2, len(exposeD12), len(nonExposeD12))
        nonExposeD12 = np.random.choice(nonExposeD12, max(10 * len(exposeD12), 5000), replace=False)
        queue.put([j, exposeD12, nonExposeD12])


def consumer(queue, counter, counter2, fout=None, caches=None, maxCache=5000):
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
        with counter2.get_lock():
            counter2.value += 1
            if counter2.value % 1000 == 0:
                print("\r%s"% counter2.value, end="")

        # print(drugJader,">>", drugBankName)
        if fout is not None:

            if caches is None:
                fout.write("%s\t%s\t%s\n" % (i, ",".join(expose), ",".join(nonExpose)))
            else:
                caches.append("%s\t%s\t%s\n" % (i, ",".join(expose), ",".join(nonExpose)))
                if len(caches) >= maxCache:
                    for line in caches:
                        fout.write("%s" % line)
                    fout.flush()
                    caches.clear()


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

    nInputList = len(drugPairList)
    inputList = [i for i in range(nInputList)]
    print("Len input: ", len(inputList))
    nDPerWorker = int(nInputList / params.N_DATA_WORKER)
    # assert 'g-csf' in allDrugNames
    for i in range(params.N_DATA_WORKER):
        startInd = i * nDPerWorker
        endInd = (i + 1) * nDPerWorker
        endInd = min(endInd, nInputList)
        if i == params.N_DATA_WORKER - 1:
            endInd = nInputList
        data = inputList[startInd:endInd], drugPairList, dDrug2Id, dList
        producers.append(Process(target=producer, args=(queue, data)))

    fout = open("%s/%s" % (OUT_DIR, "rawExpose"), "w")
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
    # exportOData()
    pExport()
