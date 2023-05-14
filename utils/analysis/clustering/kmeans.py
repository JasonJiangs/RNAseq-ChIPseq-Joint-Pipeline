import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

class clustering():
    def __init__(self):
        self.clustering = KMeans(n_clusters=2)
        self.scaler = StandardScaler()
        self.clustering_result_dir = None