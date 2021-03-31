from utils.logger.logger2 import MyLogger
from models.trainFFNN import FFNNModel
import params
import random
import torch
import numpy as np
from utils import utils

def resetRandomSeed():
    random.seed(params.TORCH_SEED)
    torch.manual_seed(params.TORCH_SEED)
    np.random.seed(params.TORCH_SEED)


class ModelRunner:
    def __init__(self):
        resetRandomSeed()

        self.data = None

        utils.ensure_dir("%s/logs" % params.C_DIR)

        PREX = "FDA"
        logPath = "%s/logs/%s_%s_%s" % (params.C_DIR, PREX,"FFNN", utils.getCurrentTimeString())
        self.model = FFNNModel()
        self.logger = MyLogger(logPath)
        self.model.setLogger(self.logger)

    def run(self):
        folds = np.arange(0, params.K_FOLD)
        for iFold in folds:
            self.model.train(iFold)



if __name__ == "__main__":
    modelRunner = ModelRunner()
    modelRunner.run()