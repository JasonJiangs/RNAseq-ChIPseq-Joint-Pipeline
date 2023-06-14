import matplotlib.pyplot as plt
import numpy as np

class regression_model():
    def __init__(self):
        self.model = None
        self.result = None
        self.result_dir = 'result/regression_model'

    def fit(self, X, y):
        self.model.fit(X, y)
        self.result = self.model.predict(X)
        # save result
        np.savetxt(self.result_dir + '/result.csv', self.result, delimiter=',')

    def plot(self):
        plt.figure(figsize=(10, 10))
        plt.plot(self.result, label='result')
        plt.legend()
        plt.savefig(self.result_dir + '/result.png')
        plt.close()

    def score(self, X, y):
        return self.model.score(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def get_params(self):
        return self.model.get_params()

    def set_params(self, **params):
        self.model.set_params(**params)

    def get_coef(self):
        return self.model.coef_

    def get_intercept(self):
        return self.model.intercept_

    def get_feature_importance(self):
        return self.model.feature_importances_


if __name__ == '__main__':
    print('This is a module for regression model.')
    # MACS2'output




