import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler

class DR_LDA():
    def __init__(self, ctrl):
        self.lda = LinearDiscriminantAnalysis()
        self.scaler = StandardScaler()
        self.lda_result_dir = None