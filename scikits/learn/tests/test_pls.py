import numpy as np
from numpy.testing import assert_array_almost_equal
from scikits.learn.datasets import load_linnerud
from scikits.learn import pls

d = load_linnerud()
X = d['data_exercise']
Y = d['data_physiological']


def test_pls():
    n_components = 2
    # 1) Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
    # ===========================================================
    # Compare 2 algo.: nipals vs. svd
    # ------------------------------
    pls_bynipals = pls.PLSCanonical()
    pls_bynipals.fit(X, Y, n_components=n_components)
    pls_bysvd = pls.PLSCanonical(algorithm="svd")
    pls_bysvd.fit(X, Y, n_components=n_components)
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
    plsca = pls.PLSCanonical()
    plsca.fit(X, Y, n_components=X.shape[1])
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
    pls_bynipals = pls.PLSCanonical()
    pls_bynipals.fit(X, Y, n_components=n_components)

    pls_ca = pls_bynipals
    x_loadings = np.array(
          [[-0.66591531,  0.77356014],
           [-0.67602366, -0.62873035],
           [ 0.35892139, -0.11993352]])
    assert_array_almost_equal(pls_ca.x_loadings_, x_loadings)

    y_loadings = np.array(
          [[ 0.6147046 ,  0.37877736],
           [ 0.65625787,  0.01196607],
           [ 0.51733043, -0.93984895]])
    assert_array_almost_equal(pls_ca.y_loadings_, y_loadings)

    x_weights = np.array(
          [[-0.58989082,  0.78900159],
           [-0.77134081, -0.61352087],
           [ 0.23887693, -0.03269003]])
    assert_array_almost_equal(pls_ca.x_weights_, x_weights)

    y_weights = np.array(
          [[ 0.61330742,  0.25616374],
           [ 0.74697171,  0.11930342],
           [ 0.25668516, -0.95924284]])
    assert_array_almost_equal(pls_ca.y_weights_, y_weights)


    # 2) Regression PLS (PLS2): "Non regression test"
    # ===============================================
    pls2 = pls.PLSRegression()
    pls2.fit(X, Y, n_components=n_components)

    x_loadings = np.array(
          [[-0.66591531, -0.01976472],
           [-0.67602366, -0.35471175],
           [ 0.35892139, -1.19418263]])
    assert_array_almost_equal(pls2.x_loadings_, x_loadings)

    y_loadings = np.array(
          [[ 0.34163079,  0.33638156],
           [ 0.41608584,  0.29078321],
           [ 0.1429814 ,  0.0651935 ]])
    assert_array_almost_equal(pls2.y_loadings_, y_loadings)

    x_weights = np.array(
          [[-0.58989082,  0.46882134],
           [-0.77134081, -0.56802033],
           [ 0.23887693, -0.67643141]])
    assert_array_almost_equal(pls2.x_weights_, x_weights)

    y_weights = np.array(
          [[ 0.61330742,  0.74851767],
           [ 0.74697171,  0.64705202],
           [ 0.25668516,  0.14506885]])
    assert_array_almost_equal(pls2.y_weights_, y_weights)

    ypred = np.array(
           [[ 9.34053037e+00,   1.39572471e+02,   6.75634346e+01],
           [  8.35625179e+00,   1.28543567e+02,   6.52136587e+01],
           [  6.88442847e+00,   1.12444632e+02,   6.19046849e+01],
           [  9.48402047e+00,   1.51099673e+02,   7.30739923e+01],
           [  1.05852278e+01,   1.53487226e+02,   7.05181552e+01],
           [  8.90885030e+00,   1.38282450e+02,   6.83808597e+01],
           [  6.77862431e+00,   1.07612118e+02,   5.97520250e+01],
           [  1.04183596e+01,   1.61076073e+02,   7.50477312e+01],
           [  1.13829272e+01,   1.78285171e+02,   8.06854024e+01],
           [  1.18461685e+01,   1.78850010e+02,   7.93812540e+01],
           [  1.13542309e+01,   1.67786960e+02,   7.53148497e+01],
           [  1.20340762e+01,   1.77605016e+02,   7.80842300e+01],
           [  1.02470983e+01,   1.62671902e+02,   7.64701009e+01],
           [  1.70939023e-01,   1.74451012e+01,   3.55750249e+01],
           [  9.69490309e+00,   1.41835825e+02,   6.75198672e+01],
           [  7.17068297e+00,   1.16570160e+02,   6.30663554e+01],
           [  8.37882344e+00,   1.31621907e+02,   6.67395903e+01],
           [  1.30076335e+01,   1.91010591e+02,   8.17092462e+01],
           [  1.20067084e+01,   1.79630851e+02,   7.92341246e+01],
           [  1.09495154e+01,   1.75568297e+02,   8.07654128e+01]])

    assert_array_almost_equal(pls2.predict(X), ypred)
    

