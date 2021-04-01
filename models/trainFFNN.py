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


def weightMSELoss(target, pred):
    mask = torch.ones(target.size())
    mask[target == 0] = params.WEIGHT_ZERO
    return torch.sum(mask * (target - pred) ** 2)


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

        if params.OPTIMIZER == "Adam":
            optimizer = torch.optim.Adam(self.model.parameters(), lr=0.01)
        else:
            optimizer = torch.optim.Adagrad(self.model.parameters(), lr=0.01)

        for i in range(params.N_ITER):
            optimizer.zero_grad()
            trainInp, trainOut = polySeData.getNextMinibatchTrain(params.BATCH_SIZE, totorch=True)
            trainPred = self.model(trainInp)
            lss = weightMSELoss(trainOut, trainPred)  # lossFunc(trainOut, trainPred)
            lss.backward()
            optimizer.step()

            if i % 20 == 0:
                with torch.no_grad():
                    # testPred = []
                    # validPred = []

                    testInp, testOut, _ = polySeData.getNextMinibatchTest(-1, True)
                    validInp, validOut, _ = polySeData.getNextMinibatchValid(-1, True)

                    testPred = self.model(testInp)
                    validPred = self.model(validInp)

                    aucTest, auprTest = evalAUCAUPR1(testOut, testPred)
                    aucValid, auprValid = evalAUCAUPR1(validOut, validPred)

                    self.logger.infoAll(("Test: ", aucTest, auprTest))
                    self.logger.infoAll(("Valid: ", aucValid, auprValid))


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
