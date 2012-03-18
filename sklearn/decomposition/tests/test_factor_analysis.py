import numpy as np
from nose.tools import assert_equal

from sklearn import datasets
from sklearn.decomposition import FactorAnalysis

iris = datasets.load_iris()

def test_fa():
    """Factor analysis on iris."""
    X, y = iris.data, iris.target
    n_samples, n_features = X.shape

    fa = FactorAnalysis(n_components=2)
    fa.fit(X)
    X_transformed = fa.transform(X)
    assert_equal(X_transformed.shape, (n_samples, 2))

def _test_fa_against_mdp():
    digits = datasets.fetch_mldata("MNIST original")
    X, y = digits.data, digits.target
    X = X.astype(np.float) + np.random.normal(size=X.shape)
    from time import time
    n_samples, n_features = X.shape

    import matplotlib.pyplot as plt
    from mdp.nodes import FANode
    fa_mdp = FANode(output_dim=2)
    start = time()
    fa_mdp.train(X)
    X_mdp = fa_mdp.execute(X)
    print("mdp: %f" % (time() - start))

    color = np.array(['r', 'g', 'b']*10)
    plt.scatter(X_mdp[:,0], X_mdp[:, 1], c=color[y])

    fa = FactorAnalysis(n_components=2, random_state=10)
    start = time()
    fa.fit(X)
    X_transformed = fa.transform(X)
    print("sklearn: %f" % (time() - start))
    assert_equal(X_transformed.shape, (n_samples, 2))

    plt.figure()
    plt.scatter(X_transformed[:,0], X_transformed[:,1], c=color[y])
    plt.show()

if __name__ == "__main__":
    test_fa()
    #test_fa_against_mdp()
