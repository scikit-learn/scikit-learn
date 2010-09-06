from nose.tools import assert_true
from scikits.learn import datasets
from scikits.learn.pca import PCA

iris = datasets.load_iris()

X = iris.data

def test_pca():
    """
    PCA
    """

    pca = PCA(k=2)
    X_r = pca.fit(X).transform(X)
    assert_true(X_r.shape[1] == 2)

    pca = PCA()
    pca.fit(X)
    assert_true(pca.explained_variance_.sum() == 1.0)
