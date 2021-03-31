import params
from multiprocessing import Process, Value, Queue
import time
from utils import utils
import editdistance as ed
import pylcs
import numpy as np
import os

T_DRUG_JADER = -1

MAX_DRUG_NAME_SPACE = 6
MIN_FAERS_DRUG_COUNT = 100


def loadDrugBankNames(sName=None, onlySmallMolecule=True):
    path = "%s/DrugBank/DrugBankNames.txt" % params.DATA_DIR
    fin = open(path)
    dHardDrug = dict()
    dSyn = dict()
    allNames = set()
    while True:
        line = fin.readline()
        if line == "":
            break
        name = line.strip()
        parts = name.split("||")
        syns = parts[-1].split("|")
        dtype = parts[1]
        if dtype != 'small molecule' and onlySmallMolecule:
            continue
        print(dtype)
        hardDrug = parts[0]
        if len(hardDrug) < 2:
            continue
        rsyns = set()
        for syn in syns:
            if len(syn) < 5:
                continue
            syn = syn.strip()
            rsyns.add(syn)
            allNames.add(syn)
            dHardDrug[syn] = hardDrug

        dSyn[hardDrug] = rsyns
        dHardDrug[hardDrug] = hardDrug

        allNames.add(hardDrug)

        # if sName is not None:
        #     if name == sName:
        #         print("Found %s"%sName)
        #         exit(-1)
    fin.close()
    # print(allNames)
    print(len(dHardDrug), len(dSyn), len(allNames))
    return dHardDrug, dSyn, allNames


def loadFDrugNames(exclusivePath=None, sName=None):
    exclusiveTokens = list()
    if exclusivePath is not None:
        fin = open(exclusivePath)
        while True:
            line = fin.readline()
            if line == "":
                break
            exclusiveTokens.append(line.strip())
        fin.close()

    print(exclusiveTokens)

    path = "%s/FDrugNamesSort_0" % params.FADER_OUT
    names = set()
    fin = open(path)
    cc = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("$")
        if len(parts) < 2:
            print(line)
            print(parts[1])
            print(parts[0])
            print(len(parts[0].strip()))
            exit(-1)
        name = parts[0].lower().strip()
        count = int(parts[1])
        if count < MIN_FAERS_DRUG_COUNT:
            continue

        isCont = False
        for token in exclusiveTokens:
            if token in name:
                isCont = True
                break
        if isCont:
            continue

        if sName is not None:
            if name.__contains__(sName):
                print("?? %s " % sName)
                exit(-1)

        names.add(name)
        cc += 1
        if cc == T_DRUG_JADER:
            break
    fin.close()
    return names


def isValidMatch(sCore, sAll):
    if sAll.__contains__(","):
        return False

    idx = sAll.find(sCore)
    if idx == -1:
        return False
    if idx > 0:
        if sAll[idx - 1].isalpha():
            return False
    return True


def producer(queue, datum):
    segdrugJader, drugBank = datum
    for drugJader in segdrugJader:
        for drugBankName in drugBank:
            if isValidMatch(drugBankName, drugJader):  # drugJader.__contains__(drugBankName):
                queue.put([drugJader, drugBankName])


def consumer(queue, counter, counter2, dDrugBank, fout=None):
    while True:
        data = queue.get()
        if data is None:
            print("Receive terminate signal")
            with counter.get_lock():
                counter.value = 1
            fout.flush()
            break
        drugJader, drugBankName = data
        with counter2.get_lock():
            counter2.value += 1

        # print(drugJader,">>", drugBankName)
        if fout is not None:
            fout.write("%s||%s||%s\n" % (drugJader, drugBankName, dDrugBank[drugBankName]))


def match():
    dHardDrug, dSyns, allDrugNames = loadDrugBankNames()
    exclusivePath = "%s/Exceptions/ExclusiveTokens" % params.OUTPUT_DIR
    faersNames = loadFDrugNames(exclusivePath)  # exclusivePath="%s/Exceptions/ExclusiveTokens" % params.OUTPUT_DIR
    print("FAERS Drugs: ", len(faersNames))
    nMatch1 = 0
    mathcingSet = set()
    for jName in faersNames:
        if jName in allDrugNames:
            nMatch1 += 1
            mathcingSet.add(jName)

    f1 = open("%s/rawMatching/MatchingDrug1.txt" % params.OUTPUT_DIR, "w")
    for drug in mathcingSet:
        f1.write("%s||%s\n" % (drug, dHardDrug[drug]))
    f1.close()
    remainJader = set()

    for name in faersNames:
        if name not in mathcingSet:
            remainJader.add(name)

    producers = []
    consumers = []
    queue = Queue(params.K_FOLD)
    counter = Value('i', 0)
    counter2 = Value('i', 0)

    jaderList = sorted(list(remainJader))
    nJader = len(jaderList)
    nDPerWorker = int(nJader / params.N_DATA_WORKER)
    # assert 'g-csf' in allDrugNames
    for i in range(params.N_DATA_WORKER):
        startInd = i * nDPerWorker
        endInd = (i + 1) * nDPerWorker
        endInd = min(endInd, nJader)
        if i == params.N_DATA_WORKER - 1:
            endInd = nJader
        data = jaderList[startInd:endInd], list(allDrugNames)
        producers.append(Process(target=producer, args=(queue, data)))

    fout = open("%s/rawMatching/MatchingDrug2.txt" % params.OUTPUT_DIR, "w")
    p = Process(target=consumer, args=(queue, counter, counter2, dHardDrug, fout))
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
    print(nMatch1, nMatch1 * 1.0 / len(faersNames))

    print("New match: ", counter2.value)
    print("Exporting Salt Candidates...")
    exportCanSaltFreq()


def loadMatchingFiles(path):
    d = dict()
    fin = open(path)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.lower().strip().split("||")
        d[parts[0]] = parts[1]
    fin.close()
    return d


def loadFreq(path):
    fin = open(path)
    kvs = []
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        kv = [parts[1].strip().lower(), int(parts[0])]
        kvs.append(kv)
    fin.close()
    return kvs


def exportNoMatching():
    d1 = loadMatchingFiles("%s/rawMatching/MatchingDrug1.txt" % params.OUTPUT_DIR)
    d2 = loadMatchingFiles("%s/rawMatching/MatchingDrug2.txt" % params.OUTPUT_DIR)
    # d3 = loadMatchingFiles("%s/MatchingDrugH.txt" % params.OUTPUT_DIR)
    print("Matching targets: ", len(d1.values()), len(set(d2.values())))
    for k, v in d2.items():
        d1[k] = v
    # for k, v in d3.items():
    #     d1[k] = v
    selectedDrugs = d1.keys()
    fin = open("%s/Tmp/DrugFreq2.txt" % params.OUTPUT_DIR)
    fout2 = open("%s/rawMatching/NoMatchingDrugFreq.txt" % params.OUTPUT_DIR, "w")

    dMatchCout = dict()
    noMatchingList = set()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        name = parts[1].strip().lower()
        cout = int(parts[0])
        # if line.__contains__("theophyllline"):
        #     print(name, name in selectedDrugs)
        if name in selectedDrugs:
            targetName = d1[name]
            utils.add_dict_counter(dMatchCout, targetName, cout)
        else:
            fout2.write("%s" % line)
            noMatchingList.add(name)

    fin.close()
    fout2.close()


def producer2(queue, datum):
    segdrugJader, drugBank = datum
    for drugJader in segdrugJader:
        scores = []
        for drugBankName in drugBank:
            # edScore = ed.eval(drugJader, drugBankName)
            edScore = (-pylcs.lcs(drugJader, drugBankName) + ed.eval(drugJader, drugBankName)) / np.log(
                len(drugBankName))
            scores.append(edScore)
        scores = np.asarray(scores)
        args = np.argsort(scores)
        assert scores[args[0]] <= scores[args[1]]
        candidateIds = args[:10]
        candidates = [drugBank[i] for i in candidateIds]

        queue.put([drugJader, candidates])


def consumer2(queue, counter, counter2, fout=None):
    while True:
        data = queue.get()
        if data is None:
            print("Receive terminate signal")
            with counter.get_lock():
                counter.value = 1
            fout.flush()
            break
        drugJader, drugBankName = data
        with counter2.get_lock():
            counter2.value += 1

        # print(drugJader,">>", drugBankName)
        if fout is not None:
            fout.write("%s||%s\n" % (drugJader, "|".join(drugBankName)))


def exportCandidateNonMatching():
    dHardDrug, dSyn, allDrugNames = loadDrugBankNames()
    allDrugNames = list(allDrugNames)

    noMatchDrugs = list()
    fin = open("%s/rawMatching/NoMatchingDrugFreq.txt" % params.OUTPUT_DIR)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        name = parts[1].strip().lower()
        # cout = int(parts[0])
        noMatchDrugs.append(name)

    producers = []
    consumers = []
    queue = Queue(params.K_FOLD)
    counter = Value('i', 0)
    counter2 = Value('i', 0)

    jaderList = noMatchDrugs
    nJader = len(jaderList)
    nDPerWorker = int(nJader / params.N_DATA_WORKER)
    for i in range(params.N_DATA_WORKER):
        startInd = i * nDPerWorker
        endInd = (i + 1) * nDPerWorker
        endInd = min(endInd, nJader)
        if i == params.N_DATA_WORKER - 1:
            endInd = nJader
        data = jaderList[startInd:endInd], allDrugNames
        producers.append(Process(target=producer2, args=(queue, data)))

    fout = open("%s/typosMatching/CandidateMatchingDrug.txt" % params.OUTPUT_DIR, "w")
    p = Process(target=consumer2, args=(queue, counter, counter2, fout))
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


def filterCandidates1():
    fin = open("%s/typosMatching/CandidateMatchingDrug.txt" % params.OUTPUT_DIR)
    fout = open("%s/typosMatching/FilterCandidateMatchingDrug.txt" % params.OUTPUT_DIR, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        if line.__contains__("/") or line.__contains__("+") or line.__contains__(",") or \
                line[0].isnumeric() or line.startswith("unk"):
            continue
        fout.write("%s" % line)
    fin.close()
    fout.close()


def filterCandidates2():
    filterCandidates1()
    fin = open("%s/typosMatching/FilterCandidateMatchingDrug.txt" % params.OUTPUT_DIR)
    fout2 = open("%s/typosMatching/FilterCandidateMatchingDrug2.txt" % params.OUTPUT_DIR, "w")
    fout3 = open("%s/typosMatching/ReFilterCandidateMatchingDrug2.txt" % params.OUTPUT_DIR, "w")
    dHardDrug, _, _ = loadDrugBankNames()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        jaderName = parts[0]
        drugBankNames = parts[1].split("|")
        firstCandidate = drugBankNames[0]
        if jaderName[0] == firstCandidate[0] and (ed.eval(jaderName, firstCandidate) < 3 or
                                                  jaderName.startswith(firstCandidate) or
                                                  firstCandidate.startswith(jaderName)):
            fout2.write("%s||%s||%s\n" % (jaderName, firstCandidate, dHardDrug[firstCandidate]))
        else:
            fout3.write("%s" % line)

    fin.close()
    fout2.close()
    fout3.close()


def filterCandidates3():
    fin1 = open("%s/typosMatching/MatchingDrugTypos.txt" % params.OUTPUT_DIR)
    fin2 = open("%s/typosMatching/FilterCandidateMatchingDrug2.txt" % params.OUTPUT_DIR)

    fout2 = open("%s/typosMatching/ReFilterCandidateMatchingDrug3.txt" % params.OUTPUT_DIR, "w")

    dHardDrug, _, _ = loadDrugBankNames()

    line1 = set(fin1.readlines())
    line2 = fin2.readlines()
    for line in line2:
        if line in line1:
            continue
        fout2.write("%s" % line)

    fin1.close()
    fin2.close()
    fout2.close()


def exportCanSaltFreq():
    fin = open("%s/rawMatching/MatchingDrug2.txt" % params.OUTPUT_DIR)
    wordFreqs = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        # words = parts[0].split(" ")
        # for word in words:
        #     if "(" not in word and ")" not in word:
        #         utils.add_dict_counter(wordFreqs, word)
        utils.add_dict_counter(wordFreqs, parts[1])
    kvs = utils.sort_dict(wordFreqs)

    fout = open("%s/rawMatching/CandSaltFreq.txt" % params.OUTPUT_DIR, "w")
    for kv in kvs:
        k, v = kv
        if v <= 2:
            continue
        fout.write("%s\n" % (k))
    fout.close()


def exportAllDict1():
    dDict1 = dict()
    # Perfect matching:
    fin = open("%s/rawMatching/MatchingDrug1.txt" % params.OUTPUT_DIR, "r")
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        t = utils.get_insert_key_dict(dDict1, parts[0], set())
        t.add(parts[1])

    fin.close()

    # Salt:
    fin = open("%s/rawMatching/Salt.txt" % params.OUTPUT_DIR)
    lines = fin.readlines()
    salts = set()
    for salt in lines:
        salt = salt.strip()
        if salt.__contains__("#"):
            continue
        salts.add(salt)

    # Partial matching:
    dHardDrug, _, _ = loadDrugBankNames()
    fin = open("%s/rawMatching/MatchingDrug2.txt" % params.OUTPUT_DIR)
    lines = fin.readlines()
    fin.close()
    for line in lines:
        line = line.strip()
        parts = line.split("||")
        if parts[1] in salts:
            continue

        drugBankName = dHardDrug[parts[1]]
        jaderName = parts[0]

        t = utils.get_insert_key_dict(dDict1, jaderName, set())
        t.add(drugBankName)

    # Typos

    fin = open("%s/typosMatching/MatchingDrugTypos.txt" % params.OUTPUT_DIR)
    lines = fin.readlines()
    fin.close()
    for line in lines:

        line = line.strip()
        if line.__contains__("#"):
            continue
        parts = line.split("||")
        t = utils.get_insert_key_dict(dDict1, parts[0], set())
        t.add(parts[-1])

    fout = open("%s/finalMap/DrugMap1.txt" % params.OUTPUT_DIR, "w")
    fout2 = open("%s/finalMap/DrugMap2.txt" % params.OUTPUT_DIR, "w")
    for k, v in dDict1.items():
        v = list(v)
        if len(v) == 1:
            fout.write("%s||%s\n" % (k, v[0]))
        else:
            fout2.write("%s||%s\n" % (k, "|".join(v)))
    fout.close()
    fout2.close()


def finalStats():
    fin = open("%s/finalMap/FinalMap.txt" % params.OUTPUT_DIR)
    lines = fin.readlines()
    lines = [line.strip() for line in lines]
    dMap = dict()
    for line in lines:
        parts = line.split("||")
        dMap[parts[0]] = parts[1]
    fin.close()

    fin = open("%s/finalMap/FinalMapH.txt" % params.OUTPUT_DIR)
    lines = fin.readlines()
    lines = [line.strip() for line in lines]
    dMapH = dict()
    for line in lines:
        parts = line.split("||")
        dMapH[parts[0]] = parts[1]
    fin.close()

    dFreq = dict()
    fin = open("%s/Tmp/DrugFreq2.txt" % params.OUTPUT_DIR)
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        drugJader = parts[1]
        c = int(parts[0])
        dDrugBank = utils.get_dict(dMap, drugJader, -1)
        d2 = utils.get_dict(dMapH, drugJader, -1)
        if dDrugBank != -1:
            utils.add_dict_counter(dFreq, dDrugBank, c)
        elif d2 != -1:
            utils.add_dict_counter(dFreq, drugJader, c)

    kvs = utils.sort_dict(dFreq)
    fout = open("%s/FinalDrugFreq.txt" % params.OUTPUT_DIR, "w")
    for kv in kvs:
        k, v = kv
        fout.write("%.6s\t%s\n" % (v, k))
    from plotLib import plotHistD, plotCul
    plotCul(kvs[::-1], 50, 2, "SelectedDrugCutOff", xLabel="ThreshHold: Freq >=", yLabel="Number of Drugs")

    fout.close()
    from plotLib import plotHistD, plotCul
    plotCul(kvs[::-1], 20, 1, "SelectedDrugCutOff", xLabel="ThreshHold: Freq >=", yLabel="Number of Drugs")


def statsSelectMatching():
    from plotLib import plotHistD, plotCul
    kvs = loadFreq("%s/SelectedDrugFreq.txt" % params.OUTPUT_DIR)

    plotCul(kvs[::-1], 20, 1, "SelectedDrugCutOff", xLabel="ThreshHold: Freq >=", yLabel="Number of Drugs")


def delFile(path):
    try:
        os.remove(path)
    except OSError:
        pass


def clean():
    delFile("%s/typosMatching/FilterCandidateMatchingDrug.txt" % params.OUTPUT_DIR)


def filterDrugMatching2():
    salts = set()
    saltStrings = open("%s/rawMatching/Salt.txt" % params.OUTPUT_DIR).readlines()
    for line in saltStrings:
        salts.add(line.strip())
    fin = open("%s/rawMatching/MatchingDrug2.txt" % params.OUTPUT_DIR)
    d = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("||")
        match2 = utils.get_insert_key_dict(d, parts[0], set())
        mD = parts[-1]
        if mD in salts and not parts[0].startswith('hydro'):
            continue
        if parts[0].__contains__('sodium chloride'):
            if mD == 'chloride ion' or mD == 'chlorine':
                continue
        match2.add(parts[-1])

    fin.close()
    fout1 = open("%s/rawMatching/MatchingDrug2_1.txt" % params.OUTPUT_DIR, "w")
    fout2 = open("%s/rawMatching/MatchingDrug2_2.txt" % params.OUTPUT_DIR, "w")

    for k, v in d.items():
        v = list(v)
        if len(v) == 1:
            fout1.write("%s||%s\n" % (k, v[0]))
        else:

            fout2.write("%s||%s\n" % (k, ",".join(v)))

    fout1.close()
    fout2.close()


def exportFinalMap():
    cmd = "cat '%s/rawMatching/MatchingDrug1.txt' '%s/rawMatching/MatchingDrug2_1.txt' >>" \
          " '%s/finalMap/DrugFMap.txt'" % (params.OUTPUT_DIR, params.OUTPUT_DIR, params.OUTPUT_DIR)
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    match()
    filterDrugMatching2()
    exportFinalMap()

    pass
