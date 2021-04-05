from models.ffnn import FFNN
import params
import torch
import numpy as np
from dataProcessing.dataFactory import PolySEData
from sklearn.metrics import roc_auc_score, average_precision_score
import inspect

from utils import utils


def getMSE(a1, a2):
    v = a1 - a2
    v = np.multiply(v, v)
    return np.sqrt(np.sum(v) / (v.shape[0] * v.shape[0]))


def weightMSELoss(target, pred, device=None, avg=True):
    mask = torch.ones(target.size())
    if device is not None:
        mask = mask.to(device)
    mask[target == 0] = params.WEIGHT_ZERO
    loss = torch.sum(mask * (target - pred) ** 2)
    if avg:
        loss /= target.size()[0]
    return loss


class FFNNModel:
    def __init__(self):
        if params.FORCE_CPU:
            self.device = torch.device('cpu')
        else:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.name = "FFNN"
        self.isFitAndPredict = True

    def setLogger(self, logger):
        self.logger = logger
        self.logger.infoAll(inspect.getsource(FFNN))

    def train(self, iFold):
        print("Fold: ", iFold)
        polySeData = PolySEData(iFold, nmaxDrug=params.MAX_N_DRUG, nChanel=params.N_CHANEL)

        model = FFNN(params.EMBEDDING_SIZE, polySeData.nSe, polySeData.nD, nLayer=params.N_LAYER,
                     device=self.device)
        self.model = model.to(self.device)

        lossFunc = torch.nn.MSELoss()
        topks = [i  for i in range(1, 20)]

        if params.OPTIMIZER == "Adam":
            optimizer = torch.optim.Adam(self.model.parameters(), lr=0.01)
        else:
            optimizer = torch.optim.Adagrad(self.model.parameters(), lr=0.01)

        for i in range(params.N_ITER):
            optimizer.zero_grad()
            trainInp, trainOut, _ = polySeData.getNextMinibatchTrain(params.BATCH_SIZE, totorch=True, device=self.device)
            # self.db(polySeData)
            trainPred = self.model(trainInp)
            lss = weightMSELoss(trainOut, trainPred, self.device)  # lossFunc(trainOut, trainPred)
            lss.backward()
            optimizer.step()
            # print(lss)

            if i % 20 == 0:
                with torch.no_grad():
                    self.logger.infoAll(("ITER: ", i))
                    # testPred = []
                    # validPred = []
                    polySeData.resetOnePassIndx()
                    # print("DB BEFORE")
                    # db(polySeData)
                    testInp, testOut, _ = polySeData.getNextMinibatchTest(-1, totorch=True, device=self.device)
                    validInp, validOut, _ = polySeData.getNextMinibatchValid(-1, totorch=True, device=self.device)

                    polySeData.resetOnePassIndx()
                    testInp2, testOut2, _ = polySeData.getNextMinibatchTest(-1)

                    # testInpNumpy = testInp.cpu().detach().numpy()
                    # testOutNumpy = testOut.cpu().numpy()

                    # print(testInpNumpy.shape, testOutNumpy.shape)
                    # print(torch.nonzero(testInp[0]).squeeze(), torch.sum(testInp))
                    # print(np.nonzero(testInpNumpy[0]))
                    # print(np.nonzero(testInp2[0]), np.sum(testInp2))

                    # print("DB AFTER")

                    # db(polySeData)
                    # print("DB2")
                    # db2(testInpNumpy, testOutNumpy, polySeData)
                    # print(testInp.shape)
                    testPred = self.model(testInp)
                    # validPred = self.model(validInp)

                    # aucTest, auprTest = evalAUCAUPR1(testOut, testPred.cpu().detach())
                    # aucValid, auprValid = evalAUCAUPR1(validOut, validPred.cpu().detach())
                    # self.logger.infoAll(("Test: ", aucTest, auprTest))
                    # self.logger.infoAll(("Valid: ", aucValid, auprValid))

                    # precTests, recallTests = evalX(testOut, testPred.cpu().detach(), topks)
                    # precVals, recallValis = evalX(validOut, validPred.cpu().detach(), topks)
                    #
                    # self.logger.infoAll(("Test: ", precTests, recallTests))
                    # self.logger.infoAll(("Valid: ", precVals, recallValis))

                    # evalX(testOut[:10], testPred.cpu().detach()[:10], topks, testInp.cpu().detach()[:10], polySeData)
                    # evalX(testOut, testPred.cpu().detach(), topks)
                    # evalX(testOut, testPred.cpu().detach(), topks, testInp.cpu().detach(), polySeData)
                    # evalX(testOut, testPred.cpu().detach(), topks,None, None)

                    prec, recall = evalX(testOut.cpu().detach(), testPred.cpu().detach(), topks, None, None)
                    self.logger.infoAll(("Prec: ", prec))
                    self.logger.infoAll(("Recal: ", recall))

                    # exit(-1)


def evalX(target, pred, topks, input=None, polySE=None):
    if input is not None:
        # print(input.shape, target.shape, pred.shape)

        # db2(input.numpy(), target.numpy(), polySE)
        _, predTopK = torch.topk(pred, 5)
        predTopK = predTopK.numpy()
        ids = [i for i in range(2)]
        # print(drugs, ses, predTopK)
        for id in ids:
            drugx = torch.nonzero(input[id]).squeeze().numpy()
            sesx = torch.nonzero(target[id]).squeeze().numpy()
            psesx = predTopK[id]
            drugx = np.atleast_1d(drugx)
            sesx = np.atleast_1d(sesx)
            # print(drugx, sesx, drugx.shape, sesx.shape)
            drugNames = [polySE.id2Drug[d] for d in drugx]
            seNames = [polySE.id2Se[s] for s in sesx]
            pses = [polySE.id2Se[s] for s in psesx]
            print(id)
            print("Drugs: ", ",".join(drugNames))
            print("Target Ses: ", ",".join(seNames))
            print("Predicted ses: ", ",".join(pses))
        # return

    precisionKs, recallKs = [], []
    for topk in topks:
        _, predTopKIndices = torch.topk(pred, topk)

        prec, recall = getPrecisionRecall(target, predTopKIndices)
        precisionKs.append(prec)
        recallKs.append(recall)
    # print(precisionKs)
    # print(recallKs)
    return precisionKs, recallKs


def db2(testInpDrug, testOutSe, polySE):
    id0 = 1
    print("DB222222")
    t3 = testInpDrug[id0]
    to3 = testOutSe[id0]
    nzd = np.nonzero(t3)[0]
    nzs = np.nonzero(to3)[0]
    print(nzd, nzs)
    dId2Drug = polySE.id2Drug
    dId2Se = polySE.id2Se

    drugNames = [dId2Drug[i] for i in nzd]
    seNames = [dId2Se[i] for i in nzs]
    print(",".join(drugNames))
    print(",".join(seNames))
    print("__________________________")


def db(polySE):
    print("___________DB_______________")
    from utils import utils
    polySE.resetOnePassIndx()
    matInp, matOut, _ = polySE.getNextMinibatchTest(-1)

    for ii in range(2):
        print(ii)
        t3 = matInp[ii]
        to3 = matOut[ii]
        nzd = np.nonzero(t3)[0]
        nzs = np.nonzero(to3)[0]
        dId2Drug = utils.reverse_dict(polySE.dDrug)
        dId2Se = utils.reverse_dict(polySE.dSe)

        drugNames = [dId2Drug[i] for i in nzd]
        seNames = [dId2Se[i] for i in nzs]
        print(",".join(drugNames))
        print(",".join(seNames))
    print("__________________________")


def getPrecisionRecall(target, predictIndices):
    pred = torch.zeros(target.size())
    nP, topk = predictIndices.shape[0], predictIndices.shape[-1]
    # print(predictIndices.shape)
    # print(predictIndices)
    # print(nP, topk, target.shape, predictIndices.shape)
    setValue(pred, predictIndices, 1)
    # pred[predictIndices] = 1
    pred[target == 0] = 0

    s0 = torch.sum(pred, dim=-1)
    ln = torch.sum(target, dim=-1)
    # print(s0, s0.shape, ln, ln.shape)
    # print(torch.sum(s0), torch.sum(ln))
    prec = s0 / topk
    recall = s0 / ln
    prec = torch.sum(prec) / nP
    recall = torch.sum(recall) / nP

    return prec, recall


def evalX2(target, pred, topks):
    precisionKs, recallKs = [], []
    for topk in topks:
        _, predTopKIndices = torch.topk(pred, topk)
        prec, recall = getPrecisionRecall(target, predTopKIndices)
        precisionKs.append(prec)
        recallKs.append(recall)
    return precisionKs, recallKs


def evalAUCAUPR1(target, pred):
    target = target.reshape(-1)
    pred = pred.reshape(-1)

    # trueOut = 1 - trueOut
    # predicted = 1 - predicted
    aupr = average_precision_score(target, pred)
    auc = roc_auc_score(target, pred)
    return auc, aupr



def convertToIndices(indices):
    nD, nC = indices.shape
    ar = torch.arange(0, nD)
    ar = ar.repeat(nC, 1).t()
    ar = torch.vstack([ar.reshape(-1), indices.reshape(-1)])
    return ar

def setValue(tensor, indices, value):
    indices = convertToIndices(indices)
    tensor[tuple(indices)] = value

if __name__ == "__main__":
    pass
