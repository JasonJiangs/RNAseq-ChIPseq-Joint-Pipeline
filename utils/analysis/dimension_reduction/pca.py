import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')


class DR_PCA():
    def __init__(self, data):
        self.pca = PCA(n_components=2)
        self.data = data
        self.result = None
        self.pca_result_dir = 'result/fpkm_pca'

    def scale(self):
        # log2 transformation
        self.data = np.log2(self.data + 1)

    def fit(self):
        pca_result = self.pca.fit_transform(self.data.T)
        self.result = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
        # save result with row names
        self.result.index = self.data.columns
        self.result.to_csv(self.pca_result_dir + '/pca_result.csv')

    def plot(self):
        # plot
        plt.figure(figsize=(10, 10))
        sns.scatterplot(
            x="PC1", y="PC2",
            data=self.result,
            legend="full",
            alpha=0.3
        )
        plt.title('PCA result')
        plt.savefig(self.pca_result_dir + '/pca_result.png')
        plt.close()

    def find_load(self):
        # the loadings of PC2
        loadings = self.pca.components_[1]
        # the indices of the top 10 genes with the highest loadings
        top10 = np.argsort(loadings)[-10:]
        # the indices of the top 10 genes with the lowest loadings
        bottom10 = np.argsort(loadings)[:10]
        # the names of the top 10 genes with the highest loadings
        top10_genes = self.data.index[top10]
        # the names of the top 10 genes with the lowest loadings
        bottom10_genes = self.data.index[bottom10]
        # the loadings of the top 10 genes with the highest loadings
        top10_loadings = loadings[top10]
        # the loadings of the top 10 genes with the lowest loadings
        bottom10_loadings = loadings[bottom10]
        # the names of the top 10 genes with the highest loadings and their loadings
        top10_genes_loadings = pd.DataFrame(
            data={'gene': top10_genes, 'loading': top10_loadings})
        # the names of the top 10 genes with the lowest loadings and their loadings
        bottom10_genes_loadings = pd.DataFrame(
            data={'gene': bottom10_genes, 'loading': bottom10_loadings})
        # save result
        top10_genes_loadings.to_csv(
            self.pca_result_dir + '/top10_genes_loadings.csv')
        bottom10_genes_loadings.to_csv(
            self.pca_result_dir + '/bottom10_genes_loadings.csv')



if __name__ == '__main__':
    # read data
    data = pd.read_csv('../data/rnaseq/gene_fpkm_matrix.csv', index_col=0)
    # drop the last two columns
    # data = data.iloc[:, :-2]
    # check correlation between first two columns
    print(data.iloc[:, 0].corr(data.iloc[:, 1]))
    print(data.iloc[:, 2].corr(data.iloc[:, 3]))
    model = DR_PCA(data)
    model.scale()
    model.fit()
    model.plot()
    model.find_load()
