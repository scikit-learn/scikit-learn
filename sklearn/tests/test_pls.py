import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.datasets import load_linnerud
from sklearn import pls

d = load_linnerud()
X = d.data
Y = d.target


def test_pls():
    n_components = 2
    # 1) Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
    # ===========================================================
    # Compare 2 algo.: nipals vs. svd
    # ------------------------------
    pls_bynipals = pls.PLSCanonical(n_components=n_components)
    pls_bynipals.fit(X, Y)
    pls_bysvd = pls.PLSCanonical(algorithm="svd", n_components=n_components)
    pls_bysvd.fit(X, Y)
    # check that the loading vectors are highly correlated
    assert_array_almost_equal(
        [np.abs(np.corrcoef(pls_bynipals.x_loadings_[:, k],
                               pls_bysvd.x_loadings_[:, k])[1, 0])
                           for k in xrange(n_components)],
        np.ones(n_components),
        err_msg="nipals and svd implementation lead to different x loadings")

    assert_array_almost_equal(
        [np.abs(np.corrcoef(pls_bynipals.y_loadings_[:, k],
                               pls_bysvd.y_loadings_[:, k])[1, 0])
                           for k in xrange(n_components)],
        np.ones(n_components),
        err_msg="nipals and svd implementation lead to different y loadings")

    # Check PLS properties (with n_components=X.shape[1])
    # ---------------------------------------------------
    plsca = pls.PLSCanonical(n_components=X.shape[1])
    plsca.fit(X, Y)
    T = plsca.x_scores_
    P = plsca.x_loadings_
    Wx = plsca.x_weights_
    U = plsca.y_scores_
    Q = plsca.y_loadings_
    Wy = plsca.y_weights_

    def check_ortho(M, err_msg):
        K = np.dot(M.T, M)
        assert_array_almost_equal(K, np.diag(np.diag(K)), err_msg=err_msg)

    # Orthogonality of weights
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(Wx, "x weights are not orthogonal")
    check_ortho(Wy, "y weights are not orthogonal")

    # Orthogonality of latent scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(T, "x scores are not orthogonal")
    check_ortho(U, "y scores are not orthogonal")

    # Check X = TP' and Y = UQ' (with (p == q) components)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # center scale X, Y
    Xc, Yc, x_mean, y_mean, x_std, y_std =\
         pls._center_scale_xy(X.copy(), Y.copy(), scale=True)
    assert_array_almost_equal(Xc, np.dot(T, P.T),
        err_msg="X != TP'")
    assert_array_almost_equal(Yc, np.dot(U, Q.T),
        err_msg="Y != UQ'")

    # Check that rotations on training data lead to scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Xr = plsca.transform(X)
    assert_array_almost_equal(Xr, plsca.x_scores_,
        err_msg="rotation on X failed")
    Xr, Yr = plsca.transform(X, Y)
    assert_array_almost_equal(Xr, plsca.x_scores_,
        err_msg="rotation on X failed")
    assert_array_almost_equal(Yr, plsca.y_scores_,
        err_msg="rotation on Y failed")

    # "Non regression test" on canonical PLS
    # --------------------------------------
    pls_bynipals = pls.PLSCanonical(n_components=n_components)
    pls_bynipals.fit(X, Y)

    pls_ca = pls_bynipals
    x_loadings = np.array(
       [[-0.61470416,  0.37877695],
        [-0.65625755,  0.01196893],
        [-0.51733059, -0.93984954]])
    assert_array_almost_equal(pls_ca.x_loadings_, x_loadings)

    y_loadings = np.array(
        [[0.66591533,  0.77358148],
         [0.67602364, -0.62871191],
         [-0.35892128, -0.11981924]])
    assert_array_almost_equal(pls_ca.y_loadings_, y_loadings)

    x_weights = np.array(
        [[-0.61330704,  0.25616119],
         [-0.74697144,  0.11930791],
         [-0.25668686, -0.95924297]])
    assert_array_almost_equal(pls_ca.x_weights_, x_weights)

    y_weights = np.array(
        [[0.58989127,  0.7890047],
         [0.77134053, -0.61351791],
         [-0.2388767, -0.03267062]])
    assert_array_almost_equal(pls_ca.y_weights_, y_weights)

    # 2) Regression PLS (PLS2): "Non regression test"
    # ===============================================
    pls2 = pls.PLSRegression(n_components=n_components)
    pls2.fit(X, Y)

    x_loadings = np.array(
        [[-0.61470416, -0.24574278],
         [-0.65625755, -0.14396183],
         [-0.51733059,  1.00609417]])
    assert_array_almost_equal(pls2.x_loadings_, x_loadings)

    y_loadings = np.array(
        [[0.32456184,  0.29892183],
         [0.42439636,  0.61970543],
         [-0.13143144, -0.26348971]])
    assert_array_almost_equal(pls2.y_loadings_, y_loadings)

    x_weights = np.array(
        [[-0.61330704, -0.00443647],
         [-0.74697144, -0.32172099],
         [-0.25668686,  0.94682413]])
    assert_array_almost_equal(pls2.x_weights_, x_weights)

    y_weights = np.array(
        [[0.58989127,  0.40572461],
         [0.77134053,  0.84112205],
         [-0.2388767, -0.35763282]])
    assert_array_almost_equal(pls2.y_weights_, y_weights)

    ypred_2 = np.array(
        [[180.33278555,   35.57034871,   56.06817703],
         [192.06235219,   37.95306771,   54.12925192]])

    assert_array_almost_equal(pls2.predict(X[:2]), ypred_2)
