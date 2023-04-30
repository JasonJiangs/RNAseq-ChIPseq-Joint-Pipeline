import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso

class Regression():
    def __init__(self):
        self.regression_models['linear'] = LinearRegression()
        self.regression_models['ridge'] = Ridge()
        self.regression_models['lasso'] = Lasso()