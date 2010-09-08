import numpy as np
from nose.tools import assert_true


from .. import datasets
from ..pca import PCA,  _assess_dimension_, _infer_dimension_

iris = datasets.load_iris()

X = iris.data

def test_pca():
    """
    PCA
    """
    pca = PCA(n_comp=2)
    X_r = pca.fit(X).transform(X)
    np.testing.assert_equal(X_r.shape[1], 2)

    pca = PCA()
    pca.fit(X)
    np.testing.assert_almost_equal(pca.explained_variance_.sum(),  1.0, 3)

def test_pca_check_project():
    """test that the projection of data is correct
    """
    n, p = 100, 3
    X = np.random.randn(n, p)*.1
    X[:10] += np.array([3, 4, 5])
    pca = PCA(n_comp=2)
    X_r = pca.fit(X)
    Xt = 0.1* np.random.randn(1, p) + np.array([3, 4, 5])
    Yt = pca.transform(Xt)
    Yt /= np.sqrt((Yt**2).sum())
    np.testing.assert_almost_equal(Yt[0][0], 1., 1)

def test_pca_dim():
    """
    """
    n, p = 100, 5
    X = np.random.randn(n, p)*.1
    X[:10] += np.array([3, 4, 5, 1, 2])
    pca = PCA(n_comp='mle')
    pca.fit(X)
    assert_true(pca.n_comp==1)

def test_infer_dim_1():
    """
    """
    n, p = 100, 5
    X = np.random.randn(n, p)*.1
    X[:10] += np.array([3, 4, 5, 1, 2])
    pca = PCA(n_comp=p)
    pca.fit(X)
    spect = pca.explained_variance_
    ll = []
    for k in range(p):
         ll.append(_assess_dimension_(spect, k, n, p))
    ll = np.array(ll)
    assert_true(ll.argmax()==1)

def test_infer_dim_2():
    """
    """
    n, p = 100, 5
    X = np.random.randn(n, p)*.1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    pca = PCA(n_comp=p)
    pca.fit(X)
    spect = pca.explained_variance_
    print _infer_dimension_(spect, n, p)
    assert_true(_infer_dimension_(spect, n, p)==2)

def test_infer_dim_3():
    """
    """
    n, p = 100, 5
    X = np.random.randn(n, p)*.1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    X[30:40] += 2*np.array([-1, 1, -1, 1, -1])
    pca = PCA(n_comp=p)
    pca.fit(X)
    spect = pca.explained_variance_
    assert_true(_infer_dimension_(spect, n, p)==3)

