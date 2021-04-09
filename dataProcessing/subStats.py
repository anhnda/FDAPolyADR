import params
from utils import utils
import numpy as np
from utils import utils
import params
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Process, Value, Queue
import time
P_THRESHOLD = 0.05
def exportSub():
    fin = open("%s/FDrug2SeList_19814.txt" % params.FADER_OUT)
    foutDict = dict()
    dlen2SeCount = dict()
    nA = 0
    print("Reading...")


    while True:
        line = fin.readline()
        if line == "":
            break
        nA += 1
        parts = line.strip().split("$")
        drugCmb = parts[0]
        ses = parts[1]
        drugs = drugCmb.split(",")
        nD = len(drugs)
        sortNames = ",".join(sorted(drugs))
        fO = utils.get_insert_key_dict(foutDict, nD, open("%s/SUB/%s" % (params.FADER_OUT, nD), "w"))
        fO.write("%s$%s\n" % (sortNames, ses))

        len2SeCount = utils.get_insert_key_dict(dlen2SeCount, nD, dict())
        sess = ses.split(",")
        for se in sess:
            utils.add_dict_counter(len2SeCount, se)

    for k, v in foutDict.items():
        v.close()


    d2 = dict()
    for k, v in dlen2SeCount:
        kvs = utils.sort_dict(v)
        ks = []
        for kv in kvs:
            kk, _ = kv
            ks.append(kk)
        d2[k] = ks
    utils.save_obj(d2, "%s/SUB/drugSize2CommonSEs" % params.FADER_OUT)




def producer(queue, arrs):
    for comAr in arrs:
        com, ar = comAr
        p = fisher_exact(ar, 'greater')
        queue.put([com, ar[0,0], p])


def consumer(queue, counter, counter2, fout=None, caches = None):
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
                    caches.append("%s\t%s\t%s\t%s\n" % (com,cc, ord, pv))


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
        nSe = utils.get_dict(dSe,se, 0)
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




def exportAllSes1():

    seList = utils.load_obj( "%s/SUB/drugSize2CommonSEs" % params.FADER_OUT)[1]
    # seList = ['product dose omission']
    nSize = 50
    import os
    p = "%s/FSUBTEST/1/*" % params.FADER_OUT
    p = p.replace(" ", '\ ')

    cmd = "rm %s" % p

    os.system(cmd)
    pathInfo1 = "%s/FSUBTEST/1/FileMap.txt" % params.FADER_OUT
    pathIn1 = "%s/SUB/%s" % (params.FADER_OUT, 1)
    dirOut1 = "%s/FSUBTEST/1" % params.FADER_OUT
    fFileNameMap = open(pathInfo1, "w")
    fFileNameMap.close()

    nSeg = max(int(len(seList) / nSize), 1)

    for i in range(nSeg):
        start = i * nSize
        end = min((i + 1) * nSize, len(seList))
        exportBySE(seList[start:end], pathIn1, dirOut1, pathInfo1)



if __name__ == "__main__":
    exportSub()
