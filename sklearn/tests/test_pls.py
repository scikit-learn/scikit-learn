import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.datasets import load_linnerud
from sklearn import pls

d = load_linnerud()
X = d.data
Y = d.target


def test_pls():
    # 1) Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
    # ===========================================================
    # Compare 2 algo.: nipals vs. svd
    # ------------------------------
    pls_bynipals = pls.PLSCanonical(n_components=X.shape[1])
    pls_bynipals.fit(X, Y)
    pls_bysvd = pls.PLSCanonical(algorithm="svd", n_components=X.shape[1])
    pls_bysvd.fit(X, Y)
    # check equalities of loading (up to the sign of the second column)
    assert_array_almost_equal(
        pls_bynipals.x_loadings_,
        np.multiply(pls_bysvd.x_loadings_, np.array([1, -1, 1])), decimal=5,
        err_msg="nipals and svd implementation lead to different x loadings")

    assert_array_almost_equal(
        pls_bynipals.y_loadings_,
        np.multiply(pls_bysvd.y_loadings_, np.array([1, -1, 1])), decimal=5,
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
    pls_ca = pls.PLSCanonical(n_components=X.shape[1])
    pls_ca.fit(X, Y)

    x_weights = np.array(
        [[-0.61330704,  0.25616119, -0.74715187],
         [-0.74697144,  0.11930791,  0.65406368],
         [-0.25668686, -0.95924297, -0.11817271]])
    assert_array_almost_equal(pls_ca.x_weights_, x_weights)

    x_rotations = np.array(
        [[-0.61330704,  0.41591889, -0.62297525],
         [-0.74697144,  0.31388326,  0.77368233],
         [-0.25668686, -0.89237972, -0.24121788]])
    assert_array_almost_equal(pls_ca.x_rotations_, x_rotations)

    y_weights = np.array(
        [[+0.58989127,  0.7890047,   0.1717553],
         [+0.77134053, -0.61351791,  0.16920272],
         [-0.23887670, -0.03267062,  0.97050016]])
    assert_array_almost_equal(pls_ca.y_weights_, y_weights)

    y_rotations = np.array(
        [[+0.58989127,  0.7168115,  0.30665872],
         [+0.77134053, -0.70791757,  0.19786539],
         [-0.23887670, -0.00343595,  0.94162826]])
    assert_array_almost_equal(pls_ca.y_rotations_, y_rotations)

    # 2) Regression PLS (PLS2): "Non regression test"
    # ===============================================
    pls_2 = pls.PLSRegression(n_components=X.shape[1])
    pls_2.fit(X, Y)

    x_weights = np.array(
        [[-0.61330704, -0.00443647,  0.78983213],
         [-0.74697144, -0.32172099, -0.58183269],
         [-0.25668686,  0.94682413, -0.19399983]])
    assert_array_almost_equal(pls_2.x_weights_, x_weights)

    x_loadings = np.array(
        [[-0.61470416, -0.24574278,  0.78983213],
         [-0.65625755, -0.14396183, -0.58183269],
         [-0.51733059,  1.00609417, -0.19399983]])
    assert_array_almost_equal(pls_2.x_loadings_, x_loadings)

    y_weights = np.array(
        [[+0.32456184,  0.29892183,  0.20316322],
         [+0.42439636,  0.61970543,  0.19320542],
         [-0.13143144, -0.26348971, -0.17092916]])
    assert_array_almost_equal(pls_2.y_weights_, y_weights)

    y_loadings = np.array(
        [[+0.32456184,  0.29892183,  0.20316322],
         [+0.42439636,  0.61970543,  0.19320542],
         [-0.13143144, -0.26348971, -0.17092916]])
    assert_array_almost_equal(pls_2.y_loadings_, y_loadings)
