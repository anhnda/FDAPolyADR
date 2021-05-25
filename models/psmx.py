from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from statsmodels.distributions.empirical_distribution import ECDF
from scipy import stats
from collections import Counter
from itertools import chain
import statsmodels.api as sm
import patsy
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys


class PSMX:
    def __init__(self, test, control, norm=True):
        self.test = test
        self.control = control

    def normalized(self):
        pass
