import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class DR_PCA():
    def __init__(self, ctrl):
        self.pca = PCA(n_components=2)
        self.scaler = StandardScaler()
        self.pca_result_dir = None
