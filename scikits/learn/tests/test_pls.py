import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from scikits.learn.datasets import load_linnerud
from scikits.learn import pls

d=load_linnerud()
X = d['data_exercise']
Y = d['data_physiological']

def test_pls():
    n_components = 2
    
    # 1) Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
    # -----------------------------------------------------------    
    # Compare 2 algo.: nipals vs svd
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pls_bynipals = pls.PLS(deflation_mode="canonical")
    pls_bynipals.fit(X,Y, n_components=n_components)
    pls_bysvd = pls.PLS(deflation_mode="canonical",algorithm="svd")
    pls_bysvd.fit(X,Y, n_components=n_components)

    assert_array_almost_equal(pls_bynipals.x_loadings_,pls_bysvd.x_loadings_, 5)
    assert_array_almost_equal(pls_bynipals.y_loadings_,pls_bysvd.y_loadings_, 5)
    assert_array_almost_equal(pls_bynipals.x_scores_,pls_bysvd.x_scores_, 5)
    assert_array_almost_equal(pls_bynipals.y_scores_,pls_bysvd.y_scores_, 5)

    # Non regression test
    # ~~~~~~~~~~~~~~~~~~~
    x_loadings = np.array(
        [[-0.58989155, -0.78900503],
        [-0.77134037,  0.61351764],
        [ 0.23887653,  0.03266757]])
    assert_array_almost_equal(pls_bynipals.x_loadings_, x_loadings)
    
    y_loadings = np.array(
        [[ 0.61330741, -0.25616063],
        [ 0.7469717,  -0.11930623],
        [ 0.25668522,  0.95924333]])
    assert_array_almost_equal(pls_bynipals.y_loadings_, y_loadings)

    
    # check orthogonality of latent scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    assert_array_almost_equal(np.corrcoef(pls_bynipals.x_scores_,rowvar=0),
                              np.eye(n_components))
    assert_array_almost_equal(np.corrcoef(pls_bynipals.y_scores_,rowvar=0),
                              np.eye(n_components))
    
    ## Regression PLS (PLS 2 blocks regression mode A known as PLS2)
    ## -------------------------------------------------------------
    pls2 = pls.PLS(deflation_mode="regression")
    pls2.fit(X, Y, n_components=n_components)

    x_loadings = np.array(
        [[-0.58989155,  0.46874883],
        [-0.77134037, -0.56798901],
        [ 0.23887653, -0.67650796]])
    assert_array_almost_equal(pls2.x_loadings_, x_loadings)
    
    y_loadings = np.array(
        [[ 0.61330741,  0.74851609],
        [ 0.7469717,   0.64704492],
        [ 0.25668522,  0.14510869]])
    assert_array_almost_equal(pls2.y_loadings_, y_loadings)

