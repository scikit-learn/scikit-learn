import numpy as np
from scipy.sparse import csr_matrix

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_less, assert_greater

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA, CCIPCA

iris = datasets.load_iris()

def test_incremental_pca():
    """IPCA on dense arrays"""
    X = iris.data
    
    # test shape
    ipca = IncrementalPCA(n_components=2)
    X_r1 = ipca.fit(X).transform(X)
    np.testing.assert_equal(X_r1.shape[1], 2)
    
    # test equivalence of separate fits
    ipca = IncrementalPCA(n_components=2)
    X_r2 = ipca.fit(X).transform(X)
    assert_array_almost_equal(X_r1, X_r2)
    
    # check that the subspace represents the variance well
    ipca = IncrementalPCA(n_components=4)
    ipca.fit(X)
    assert_almost_equal(ipca.explained_variance_ratio_.sum(), 1.0, 3)
    
    # test equivalence with pca
    pca = PCA(n_components=4)
    pca.fit(X)
    
    # ipca and pca will be close but not identical, check to 3 decimal
    assert_array_almost_equal(pca.explained_variance_ratio_, ipca.explained_variance_ratio_,3)
    

def test_incremental_pca_check_projection():
    """Test that the projection of data is correct"""
    rng = np.random.RandomState(0)
    n, p = 100, 3
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * rng.randn(1, p) + np.array([3, 4, 5])

    Yt = IncrementalPCA(n_components=2).fit(X).transform(Xt)

    Yt /= np.sqrt((Yt * Yt).sum())
    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_incremental_pca_inverse():
    """Test that the projection of data can be inverted"""
    rng = np.random.RandomState(0)
    n, p = 50, 3
    X = rng.randn(n, p)  # spherical data
    X[:, 1] *= .00001  # make middle component relatively small
    X += [5, 4, 3]  # make a large mean

    # same check that we can find the original data from the transformed
    # signal (since the data is almost of rank n_components)
    ipca = IncrementalPCA(n_components=2).fit(X)
    Y = ipca.transform(X)
    Y_inverse = ipca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=3)


def test_ccipca():
    """CCIPCA on dense arrays"""
    X = iris.data
    
    # test shape
    ccipca = CCIPCA(n_components=2, amnesic=1.0)
    X_r1 = ccipca.fit(X).transform(X)
    np.testing.assert_equal(X_r1.shape[1], 2)
    
    # test equivalence of separate fits
    ccipca = CCIPCA(n_components=2, amnesic=1.0)
    X_r2 = ccipca.fit(X).transform(X)
    assert_array_almost_equal(X_r1, X_r2)
    
    # check that the subspace represents the variance well
    ccipca = CCIPCA(n_components=4, amnesic=1.0)
    ccipca.fit(X)
    assert_almost_equal(ccipca.explained_variance_ratio_.sum(), 1.0, 3)
    
    # test equivalence with pca
    pca = PCA(n_components=4)
    pca.fit(X)
    
    # ccipca and pca will be close but not identical, check to 1 decimal
    assert_array_almost_equal(pca.explained_variance_ratio_, ccipca.explained_variance_ratio_,1)
    

def test_ccipca_check_projection():
    """Test that the projection of data is correct"""
    rng = np.random.RandomState(0)
    n, p = 100, 3
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * rng.randn(1, p) + np.array([3, 4, 5])

    Yt = CCIPCA(n_components=2).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_ccipca_inverse():
    """Test that the projection of data can be inverted"""
    rng = np.random.RandomState(0)
    n, p = 50, 3
    X = rng.randn(n, p)  # spherical data
    X[:, 1] *= .00001  # make middle component relatively small
    X += [5, 4, 3]  # make a large mean

    # same check that we can find the original data from the transformed
    # signal (since the data is almost of rank n_components)
    ccipca = CCIPCA(n_components=2).fit(X)
    Y = ccipca.transform(X)
    Y_inverse = ccipca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=0) # again not as precise as PCA without iteration

if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
