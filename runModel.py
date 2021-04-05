from utils.logger.logger2 import MyLogger
from models.trainFFNN import FFNNModel
from models.trainMili import MILIModel
from models.trainMil import MILModel
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
    def __init__(self, model = "FFNN"):
        resetRandomSeed()

        self.data = None

        utils.ensure_dir("%s/logs" % params.C_DIR)

        PREX = "FDA"
        if model == "FFNN":
            self.model = FFNNModel()
        elif model == "MILI":
            self.model = MILIModel()
        elif model == "MIL":
            self.model = MILModel()
        logPath = "%s/logs/%s_%s_%s" % (params.C_DIR, PREX, self.model.name, utils.getCurrentTimeString())

        self.logger = MyLogger(logPath)
        self.model.setLogger(self.logger)
        self.logger.infoAll((self.model.name))
        # self.logger.infoAll(self.model.model.named_parameters())
        self.logger.infoAll(("LAYERS, EMBEDDING_SIZE, WEIGHT, BATCH_SIZE", params.N_LAYER, params.EMBEDDING_SIZE, params.WEIGHT_ZERO, params.BATCH_SIZE))
        self.logger.infoAll(("NCHANELS, DK, MAX DRUG: ",params.N_CHANEL, params.DK, params.MAX_N_DRUG))
    def run(self):
        folds = np.arange(1)
        for iFold in folds:
            self.model.train(iFold)


if __name__ == "__main__":
    model = params.METHOD
    modelRunner = ModelRunner(model)
    modelRunner.run()
