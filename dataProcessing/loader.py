import params
from utils import utils
import random
import numpy as np
import torch

def exportPolySes():
    dDrug = dict()
    dSe = dict()
    dDrugComb2Ses = dict()
    fin = open("%s/PolySes.txt" % params.FADER_OUT)
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("\t")
        drugCom = parts[0]
        ses = parts[1]
        drugs = drugCom.split(",")
        if len(drugs) > params.MAX_N_DRUG:
            continue

        for drug in drugs:
            utils.get_update_dict_index(dDrug, drug)

        ses = ses.split(",")
        for se in ses:
            utils.get_update_dict_index(dSe, se)

        dDrugComb2Ses[drugCom] = ses

    nDrug = len(dDrug)
    nSe = len(dSe)
    nComb = len(dDrugComb2Ses)

    print("Drugs, Ses, Comb: ", nDrug, nSe, nComb)

    fout = open("%s/PolySe_%s" % (params.FADER_OUT, params.MAX_N_DRUG), "w")
    kvs = []
    for drugCom, ses in dDrugComb2Ses.items():
        fout.write("%s\t%s\n" % (drugCom, ",".join(ses)))
        kvs.append([drugCom.split(","), ses])

    fout.close()

    random.shuffle(kvs)

    SEG_SIZE = int(nComb / params.K_FOLD)

    for iFold in range(params.K_FOLD):
        print("Generating fold...", iFold)
        tests = []
        validates = []
        trains = []
        startTest = iFold * SEG_SIZE
        endTest = (iFold + 1) * SEG_SIZE
        if endTest > nComb:
            endTest = nComb

        startValid = endTest
        if iFold == params.K_FOLD - 1:
            startValid = 0

        endValid = startValid + SEG_SIZE

        if endValid > nComb:
            endValid = nComb

        for j, kv in enumerate(kvs):
            if startTest <= j < endTest:
                seg = tests
            elif startValid <= j < endValid:
                seg = validates
            else:
                seg = trains
            seg.append(kv)

        utils.save_obj((dDrug, dSe, trains, tests, validates), "%s/_%s" % (params.FADER_KFOLD, iFold))


def loadFold(iFold):
    path = "%s/_%s" % (params.FADER_KFOLD, iFold)
    dDes = utils.load_obj("%s/DrugBank/DrugMorganDes" % params.DATA_DIR)
    dDrug, dSe, trains, tests, validates = utils.load_obj(path)
    return dDrug, dSe, trains, tests, validates, dDes


class PolySEData:
    def __init__(self, iFold):
        self.iFold = iFold
        self.__loadRawFold(iFold)
        pass

    def __loadRawFold(self, iFold):
        self.dDrug, self.dSe, self.trains, self.tests, self.validates, self.dDes = loadFold(iFold)
        self.nD = len(self.dDrug)
        self.nSe = len(self.dSe)
        self.currentTrainIdx = 0
        self.currentTestIdx = 0
        self.currentValidIdx = 0

    def __getNextRaw(self, inps, batchSize, currentIdx, shuffle=True, onePass=False):
        l = len(inps)
        data = []
        if batchSize == -1:
            batchSize = l
        end = currentIdx + batchSize
        remain = 0
        if end > l:
            remain = end - l
            end = l
        for i in range(currentIdx, end):
            data.append(inps[i])
        if remain > 0:
            if onePass:
                return data, -1

            for i in range(0, remain):
                data.append(inps[i])

            if shuffle:
                random.shuffle(inps)
            currentIdx = remain

        return data, currentIdx

    def __convertSegData2Binary(self, data):
        matInp = []
        matOut = []

        for kv in data:
            # print(kv)
            k, v = kv
            inp = np.zeros(self.nD)
            for drug in k:

                inp[self.dDrug[drug]] = 1
            out = np.zeros(self.nSe)
            for se in v:
                out[self.dSe[se]] = 1

            matInp.append(inp)
            matOut.append(out)

        matInp = np.vstack(matInp)
        matOut = np.vstack(matOut)

        return matInp, matOut

    def resetOnePassIndx(self):
        self.currentTestIdx = 0
        self.currentValidIdx = 0

    def getNextMinibatchTrain(self, batchSize, totorch=False, device=None):
        trainMinibatch, self.currentTrainIdx = self.__getNextRaw(self.trains, batchSize, self.currentTrainIdx)
        matInp, matOut = self.__convertSegData2Binary(trainMinibatch)
        if totorch:
            matInp = torch.from_numpy(matInp).float()
            matOut = torch.from_numpy(matOut).float()
            if device is not None:
                matInp = matInp.to(device)
                matOut = matOut.to(device)
        return matInp, matOut

    def getNextMinibatchTest(self, batchSize, totorch=False, device=None):

        testMinibatch, self.currentTestIdx = self.__getNextRaw(self.tests, batchSize, self.currentTestIdx, shuffle=False, onePass=True)
        matInp, matOut = self.__convertSegData2Binary(testMinibatch)
        if totorch:
            matInp = torch.from_numpy(matInp).float()
            matOut = torch.from_numpy(matOut).float()
            if device is not None:
                matInp = matInp.to(device)
                # matOut = matOut.to(device)
        return matInp, matOut, self.currentTestIdx

    def getNextMinibatchValid(self, batchSize, totorch=False, device=None):

        validMinibatch, self.currentValidIdx = self.__getNextRaw(self.validates, batchSize, self.currentValidIdx,
                                                                 shuffle=False, onePass=True)
        matInp, matOut = self.__convertSegData2Binary(validMinibatch)
        if totorch:
            matInp = torch.from_numpy(matInp).float()
            matOut = torch.from_numpy(matOut).float()
            if device is not None:
                matInp = matInp.to(device)
                # matOut = matOut.to(device)
        return matInp, matOut, self.currentTestIdx


if __name__ == "__main__":
    np.random.seed(params.TORCH_SEED)
    random.seed(params.TORCH_SEED)
    exportPolySes()
