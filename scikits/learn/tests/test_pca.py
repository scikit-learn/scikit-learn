import numpy as np

from .. import datasets
from ..pca import PCA

iris = datasets.load_iris()

X = iris.data

def test_pca():
    """
    PCA
    """

    pca = PCA(k=2)
    X_r = pca.fit(X).transform(X)
    np.testing.assert_equal(X_r.shape[1], 2)

    pca = PCA()
    pca.fit(X)
    np.testing.assert_almost_equal(pca.explained_variance_.sum(),  1.0, 3)
