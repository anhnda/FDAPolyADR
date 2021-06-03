from utils import utils
import params
import numpy as np

from scipy.stats import ttest_ind

from multiprocessing import Process, Value, Queue

OUT_DIR = params.JADER_OUT
PREF = "SJADER"
P_THRESHOLD = 0.05
import time

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

def getSubList(list, ids):
    subList = []
    for id in ids:
        subList.append(list[id])
    return subList
def producer(queue, datas):
    oRs, drugPairList, dDrug2Id, dId2Se, dList = datas
    for oR in oRs:
        pId, rExposeIds, rNonExposeIds = oR
        dPair = drugPairList[pId]
        rExpose = getSubList(dList,rExposeIds)
        # rNoneExpose = dList[rNonExposeIds]

        seSet = set()
        for r in rExpose:
            _, _, ses = r
            for se in ses:
                seSet.add(se)
        n1 = max(int(len(rExposeIds) / 10), 1)
        n2 = max(int(len(rNonExposeIds) / 10), 1)
        # print(n1, n2, len(rExposeIds), len(rNonExposeIds))
        # ar1 = np.random.choice(rExposeIds, (1000, n1), replace=False)
        # ar2 = np.random.choice(rNonExposeIds, (1000, n2), replace=False)
        nSe = len(seSet)
        dOldSeId2NewId = dict()
        for se in seSet:
            dOldSeId2NewId[se] = len(dOldSeId2NewId)
        dId2NewSeIdOld = utils.reverse_dict(dOldSeId2NewId)
        def calRatio(dList, ar, nSe, nCount):
            appears = np.zeros((1000, nSe))
            for i in range(1000):
                rIds = np.random.choice(ar, nCount, replace=False)
                rs = getSubList(dList, rIds)

                for r in rs:
                    _, _, seIds = r
                    see = []
                    for seId in seIds:
                        newSeId = utils.get_dict(dOldSeId2NewId, seId, -1)
                        if newSeId != -1:
                            see.append(newSeId)
                        appears[i, see] += 1

            notAppear = nCount - appears + 1e-10
            ratio = appears / notAppear
            return ratio

        ratioExpose = calRatio(dList, rExposeIds, nSe, n1)
        ratioNonExpose = calRatio(dList, rNonExposeIds, nSe, n2)
        sigSes = []
        for i in range(nSe):
            _, p = ttest_ind(ratioExpose[ :, i], ratioNonExpose[ :, i],  alternative="greater")
            if p <= P_THRESHOLD:
                sigSes.append([dId2NewSeIdOld[i], p])
        for v in sigSes:
            se, p = v
            seName = dId2Se[se]
            queue.put([dPair, seName, p])


def consumer(queue, counter, counter2, fout=None, caches=None, maxCache=100):
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
        drugPair, seName, p = data
        with counter2.get_lock():
            counter2.value += 1
            if counter2.value % 10 == 0:
                print("\r%s"% counter2.value, end="")

        # print(drugJader,">>", drugBankName)
        if fout is not None:

            if caches is None:
                fout.write("%s\t%s\t%s\n" % (drugPair, seName, p))
            else:
                caches.append("%s\t%s\t%s\n" % (drugPair, seName, p))
                if len(caches) >= maxCache:
                    for line in caches:
                        fout.write("%s" % line)
                    fout.flush()
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


def loadPSMExpose():
    fin = open("%s/%s" % (OUT_DIR, "psmExpose"))
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

def runTTest():
    producers = []
    consumers = []
    queue = Queue(params.K_FOLD)
    counter = Value('i', 0)
    counter2 = Value('i', 0)

    dList = utils.load_obj("%s/DataDump.o" % OUT_DIR)
    dDrugPair2Id, drugPairList = loadDictName2Id("%s/%sPairs.txt" % (OUT_DIR, PREF), nMax=-1, min=1)
    dDrug2Id, _ = loadDictName2Id("%s/%sADrugs.txt" % (OUT_DIR, PREF))
    dInd2Id, _ = loadDictName2Id("%s/%sAInd.txt" % (OUT_DIR, PREF))
    dSe2Id, _ = loadDictName2Id("%s/%sASe.txt" % (OUT_DIR, PREF))
    dId2Se = utils.reverse_dict(dSe2Id)

    inputList = loadPSMExpose() # loadRawExpose()
    nInputList = len(inputList)

    nDPerWorker = int(nInputList / params.N_DATA_WORKER)
    # assert 'g-csf' in allDrugNames
    for i in range(params.N_DATA_WORKER):
        startInd = i * nDPerWorker
        endInd = (i + 1) * nDPerWorker
        endInd = min(endInd, nInputList)
        if i == params.N_DATA_WORKER - 1:
            endInd = nInputList
        data = inputList[startInd:endInd], drugPairList, dDrug2Id, dId2Se, dList
        producers.append(Process(target=producer, args=(queue, data)))

    fout = open("%s/%s" % (OUT_DIR, "ttStatsRe"), "w")
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


def exportPolySE():
    fin = open("%s/%s" % (OUT_DIR, "ttStatsRe"))
    dDrugPair2Se = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        drugPairs = parts[0]
        se = parts[1]
        seList = utils.get_insert_key_dict(dDrugPair2Se, drugPairs, [])
        seList.append(se)
    fin.close()

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

    fout = open("%s/%s" % (OUT_DIR, "JPolySE"), "w")
    for dp, ses in dDrugPair2Se.items():
        d1, d2 = dp.split(",")
        i1, i2 = utils.get_dict(dName2Inchi, d1, -1), utils.get_dict(dName2Inchi, d2, -1)
        if i1 == -1 or i2 == -1:
            continue
        if len(i1) < 2 or len(i2) < 2:
            continue
        fout.write("%s|%s|%s|%s|%s\n" % (d1, d2, i1, i2, ",".join(ses)))
    fout.close()

if __name__ == "__main__":
    runTTest()
    exportPolySE()