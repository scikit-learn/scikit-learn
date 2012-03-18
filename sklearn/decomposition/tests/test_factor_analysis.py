import numpy as np
from nose.tools import assert_equal

from sklearn import datasets
from sklearn.decomposition import FactorAnalysis

iris = datasets.load_iris()

def test_fa():
    """Factor analysis on iris."""
    X, y = iris.data, iris.target
    n_samples, n_features = X.shape

    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    X_ = PCA(n_components=2).fit_transform(X)
    color = np.array(['r', 'g', 'b'])
    plt.scatter(X_[:,0], X[:, 1], c=color[y])

    fa = FactorAnalysis(n_components=2, tol=10)
    fa.fit(X)
    X_transformed = fa.transform(X)
    assert_equal(X_transformed.shape, (n_samples, 2))

    plt.figure()
    plt.scatter(X_transformed[:,0], X_transformed[:,1], c=color[y])
    plt.show()

if __name__ == "__main__":
    test_fa()
