from models.ffnn import FFNN
import params
import torch
import numpy as np
from dataProcessing.loader import PolySEData
from sklearn.metrics import roc_auc_score, average_precision_score
import inspect


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
        polySeData = PolySEData(iFold)

        model = FFNN(params.EMBEDDING_SIZE, polySeData.nSe, polySeData.nD, nLayer=params.N_LAYER,
                     device=self.device)
        self.model = model.to(self.device)

        lossFunc = torch.nn.MSELoss()
        topks = [i * 5 for i in range(1, 10)]

        if params.OPTIMIZER == "Adam":
            optimizer = torch.optim.Adam(self.model.parameters(), lr=0.01)
        else:
            optimizer = torch.optim.Adagrad(self.model.parameters(), lr=0.01)

        for i in range(params.N_ITER):
            optimizer.zero_grad()
            trainInp, trainOut = polySeData.getNextMinibatchTrain(params.BATCH_SIZE, totorch=True, device=self.device)

            trainPred = self.model(trainInp)
            lss = weightMSELoss(trainOut, trainPred, self.device)  # lossFunc(trainOut, trainPred)
            lss.backward()
            optimizer.step()
            # print(lss)

            if i % 20 == 0:
                with torch.no_grad():
                    print("ITER: ", i)
                    # testPred = []
                    # validPred = []
                    polySeData.resetOnePassIndx()
                    testInp, testOut, _ = polySeData.getNextMinibatchTest(-1, True, device=self.device)
                    validInp, validOut, _ = polySeData.getNextMinibatchValid(-1, True, device=self.device)

                    testPred = self.model(testInp)
                    validPred = self.model(validInp)

                    # aucTest, auprTest = evalAUCAUPR1(testOut, testPred.cpu().detach())
                    # aucValid, auprValid = evalAUCAUPR1(validOut, validPred.cpu().detach())
                    # self.logger.infoAll(("Test: ", aucTest, auprTest))
                    # self.logger.infoAll(("Valid: ", aucValid, auprValid))

                    precTests, recallTests = evalX(testOut, testPred.cpu().detach())
                    precVals, recallValis = evalX(validOut, validPred.cpu().detach())

                    self.logger.infoAll(("Test: ", precTests, recallTests))
                    self.logger.infoAll(("Valid: ", precVals, recallValis))
def getPrecisionRecall(target, predictIndices):
    pred = torch.zeros(target.size())
    topk = predictIndices.shape[-1]
    pred[predictIndices] = 1
    pred[pred == target] = 1
    s0 = torch.sum(pred, dim=-1)
    ln = topk.sum(target, dim=-1)
    prec = s0 / topk
    recall = s0 / ln
    return prec, recall

def evalX(target, pred, topks):
    precisionKs, recallKs = [], []
    for topk in topks:
        _, predTopKIndices = torch.topk(topk, pred)
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


if __name__ == "__main__":
    pass
