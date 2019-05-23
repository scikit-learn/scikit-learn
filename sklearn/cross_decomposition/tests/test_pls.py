import numpy as np
from numpy.testing import assert_approx_equal

from sklearn.utils.testing import (assert_equal, assert_array_almost_equal,
                                   assert_array_equal, assert_raise_message,
                                   assert_warns)
from sklearn.datasets import load_linnerud
from sklearn.cross_decomposition import pls_, CCA
from sklearn.preprocessing import StandardScaler
from sklearn.utils import check_random_state
from sklearn.exceptions import ConvergenceWarning


def test_pls():
    d = load_linnerud()
    X = d.data
    Y = d.target
    # 1) Canonical (symmetric) PLS (PLS 2 blocks canonical mode A)
    # ===========================================================
    # Compare 2 algo.: nipals vs. svd
    # ------------------------------
    pls_bynipals = pls_.PLSCanonical(n_components=X.shape[1])
    pls_bynipals.fit(X, Y)
    pls_bysvd = pls_.PLSCanonical(algorithm="svd", n_components=X.shape[1])
    pls_bysvd.fit(X, Y)
    # check equalities of loading (up to the sign of the second column)
    assert_array_almost_equal(
        pls_bynipals.x_loadings_,
        pls_bysvd.x_loadings_, decimal=5,
        err_msg="nipals and svd implementations lead to different x loadings")

    assert_array_almost_equal(
        pls_bynipals.y_loadings_,
        pls_bysvd.y_loadings_, decimal=5,
        err_msg="nipals and svd implementations lead to different y loadings")

    # Check PLS properties (with n_components=X.shape[1])
    # ---------------------------------------------------
    plsca = pls_.PLSCanonical(n_components=X.shape[1])
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
        pls_._center_scale_xy(X.copy(), Y.copy(), scale=True)
    assert_array_almost_equal(Xc, np.dot(T, P.T), err_msg="X != TP'")
    assert_array_almost_equal(Yc, np.dot(U, Q.T), err_msg="Y != UQ'")

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
    # The results were checked against the R-package plspm
    pls_ca = pls_.PLSCanonical(n_components=X.shape[1])
    pls_ca.fit(X, Y)

    x_weights = np.array(
        [[-0.61330704,  0.25616119, -0.74715187],
         [-0.74697144,  0.11930791,  0.65406368],
         [-0.25668686, -0.95924297, -0.11817271]])
    # x_weights_sign_flip holds columns of 1 or -1, depending on sign flip
    # between R and python
    x_weights_sign_flip = pls_ca.x_weights_ / x_weights

    x_rotations = np.array(
        [[-0.61330704,  0.41591889, -0.62297525],
         [-0.74697144,  0.31388326,  0.77368233],
         [-0.25668686, -0.89237972, -0.24121788]])
    x_rotations_sign_flip = pls_ca.x_rotations_ / x_rotations

    y_weights = np.array(
        [[+0.58989127,  0.7890047,   0.1717553],
         [+0.77134053, -0.61351791,  0.16920272],
         [-0.23887670, -0.03267062,  0.97050016]])
    y_weights_sign_flip = pls_ca.y_weights_ / y_weights

    y_rotations = np.array(
        [[+0.58989127,  0.7168115,  0.30665872],
         [+0.77134053, -0.70791757,  0.19786539],
         [-0.23887670, -0.00343595,  0.94162826]])
    y_rotations_sign_flip = pls_ca.y_rotations_ / y_rotations

    # x_weights = X.dot(x_rotation)
    # Hence R/python sign flip should be the same in x_weight and x_rotation
    assert_array_almost_equal(x_rotations_sign_flip, x_weights_sign_flip)
    # This test that R / python give the same result up to column
    # sign indeterminacy
    assert_array_almost_equal(np.abs(x_rotations_sign_flip), 1, 4)
    assert_array_almost_equal(np.abs(x_weights_sign_flip), 1, 4)


    assert_array_almost_equal(y_rotations_sign_flip, y_weights_sign_flip)
    assert_array_almost_equal(np.abs(y_rotations_sign_flip), 1, 4)
    assert_array_almost_equal(np.abs(y_weights_sign_flip), 1, 4)

    # 2) Regression PLS (PLS2): "Non regression test"
    # ===============================================
    # The results were checked against the R-packages plspm, misOmics and pls
    pls_2 = pls_.PLSRegression(n_components=X.shape[1])
    pls_2.fit(X, Y)

    x_weights = np.array(
        [[-0.61330704, -0.00443647,  0.78983213],
         [-0.74697144, -0.32172099, -0.58183269],
         [-0.25668686,  0.94682413, -0.19399983]])
    x_weights_sign_flip = pls_2.x_weights_ / x_weights

    x_loadings = np.array(
        [[-0.61470416, -0.24574278,  0.78983213],
         [-0.65625755, -0.14396183, -0.58183269],
         [-0.51733059,  1.00609417, -0.19399983]])
    x_loadings_sign_flip = pls_2.x_loadings_ / x_loadings

    y_weights = np.array(
        [[+0.32456184,  0.29892183,  0.20316322],
         [+0.42439636,  0.61970543,  0.19320542],
         [-0.13143144, -0.26348971, -0.17092916]])
    y_weights_sign_flip = pls_2.y_weights_ / y_weights

    y_loadings = np.array(
        [[+0.32456184,  0.29892183,  0.20316322],
         [+0.42439636,  0.61970543,  0.19320542],
         [-0.13143144, -0.26348971, -0.17092916]])
    y_loadings_sign_flip = pls_2.y_loadings_ / y_loadings

    # x_loadings[:, i] = Xi.dot(x_weights[:, i]) \forall i
    assert_array_almost_equal(x_loadings_sign_flip, x_weights_sign_flip, 4)
    assert_array_almost_equal(np.abs(x_loadings_sign_flip), 1, 4)
    assert_array_almost_equal(np.abs(x_weights_sign_flip), 1, 4)

    assert_array_almost_equal(y_loadings_sign_flip, y_weights_sign_flip, 4)
    assert_array_almost_equal(np.abs(y_loadings_sign_flip), 1, 4)
    assert_array_almost_equal(np.abs(y_weights_sign_flip), 1, 4)

    # 3) Another non-regression test of Canonical PLS on random dataset
    # =================================================================
    # The results were checked against the R-package plspm
    n = 500
    p_noise = 10
    q_noise = 5
    # 2 latents vars:
    rng = check_random_state(11)
    l1 = rng.normal(size=n)
    l2 = rng.normal(size=n)
    latents = np.array([l1, l1, l2, l2]).T
    X = latents + rng.normal(size=4 * n).reshape((n, 4))
    Y = latents + rng.normal(size=4 * n).reshape((n, 4))
    X = np.concatenate(
        (X, rng.normal(size=p_noise * n).reshape(n, p_noise)), axis=1)
    Y = np.concatenate(
        (Y, rng.normal(size=q_noise * n).reshape(n, q_noise)), axis=1)

    pls_ca = pls_.PLSCanonical(n_components=3)
    pls_ca.fit(X, Y)

    x_weights = np.array(
        [[0.65803719,  0.19197924,  0.21769083],
         [0.7009113,  0.13303969, -0.15376699],
         [0.13528197, -0.68636408,  0.13856546],
         [0.16854574, -0.66788088, -0.12485304],
         [-0.03232333, -0.04189855,  0.40690153],
         [0.1148816, -0.09643158,  0.1613305],
         [0.04792138, -0.02384992,  0.17175319],
         [-0.06781, -0.01666137, -0.18556747],
         [-0.00266945, -0.00160224,  0.11893098],
         [-0.00849528, -0.07706095,  0.1570547],
         [-0.00949471, -0.02964127,  0.34657036],
         [-0.03572177,  0.0945091,  0.3414855],
         [0.05584937, -0.02028961, -0.57682568],
         [0.05744254, -0.01482333, -0.17431274]])
    x_weights_sign_flip = pls_ca.x_weights_ / x_weights


    x_loadings = np.array(
        [[0.65649254,  0.1847647,  0.15270699],
         [0.67554234,  0.15237508, -0.09182247],
         [0.19219925, -0.67750975,  0.08673128],
         [0.2133631, -0.67034809, -0.08835483],
         [-0.03178912, -0.06668336,  0.43395268],
         [0.15684588, -0.13350241,  0.20578984],
         [0.03337736, -0.03807306,  0.09871553],
         [-0.06199844,  0.01559854, -0.1881785],
         [0.00406146, -0.00587025,  0.16413253],
         [-0.00374239, -0.05848466,  0.19140336],
         [0.00139214, -0.01033161,  0.32239136],
         [-0.05292828,  0.0953533,  0.31916881],
         [0.04031924, -0.01961045, -0.65174036],
         [0.06172484, -0.06597366, -0.1244497]])
    x_loadings_sign_flip = pls_ca.x_loadings_ / x_loadings

    y_weights = np.array(
        [[0.66101097,  0.18672553,  0.22826092],
         [0.69347861,  0.18463471, -0.23995597],
         [0.14462724, -0.66504085,  0.17082434],
         [0.22247955, -0.6932605, -0.09832993],
         [0.07035859,  0.00714283,  0.67810124],
         [0.07765351, -0.0105204, -0.44108074],
         [-0.00917056,  0.04322147,  0.10062478],
         [-0.01909512,  0.06182718,  0.28830475],
         [0.01756709,  0.04797666,  0.32225745]])
    y_weights_sign_flip = pls_ca.y_weights_ / y_weights

    y_loadings = np.array(
        [[0.68568625,  0.1674376,  0.0969508],
         [0.68782064,  0.20375837, -0.1164448],
         [0.11712173, -0.68046903,  0.12001505],
         [0.17860457, -0.6798319, -0.05089681],
         [0.06265739, -0.0277703,  0.74729584],
         [0.0914178,  0.00403751, -0.5135078],
         [-0.02196918, -0.01377169,  0.09564505],
         [-0.03288952,  0.09039729,  0.31858973],
         [0.04287624,  0.05254676,  0.27836841]])
    y_loadings_sign_flip = pls_ca.y_loadings_ / y_loadings

    assert_array_almost_equal(x_loadings_sign_flip, x_weights_sign_flip, 4)
    assert_array_almost_equal(np.abs(x_weights_sign_flip), 1, 4)
    assert_array_almost_equal(np.abs(x_loadings_sign_flip), 1, 4)

    assert_array_almost_equal(y_loadings_sign_flip, y_weights_sign_flip, 4)
    assert_array_almost_equal(np.abs(y_weights_sign_flip), 1, 4)
    assert_array_almost_equal(np.abs(y_loadings_sign_flip), 1, 4)

    # Orthogonality of weights
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(pls_ca.x_weights_, "x weights are not orthogonal")
    check_ortho(pls_ca.y_weights_, "y weights are not orthogonal")

    # Orthogonality of latent scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(pls_ca.x_scores_, "x scores are not orthogonal")
    check_ortho(pls_ca.y_scores_, "y scores are not orthogonal")


def test_convergence_fail():
    d = load_linnerud()
    X = d.data
    Y = d.target
    pls_bynipals = pls_.PLSCanonical(n_components=X.shape[1],
                                     max_iter=2, tol=1e-10)
    assert_warns(ConvergenceWarning, pls_bynipals.fit, X, Y)


def test_PLSSVD():
    # Let's check the PLSSVD doesn't return all possible component but just
    # the specified number
    d = load_linnerud()
    X = d.data
    Y = d.target
    n_components = 2
    for clf in [pls_.PLSSVD, pls_.PLSRegression, pls_.PLSCanonical]:
        pls = clf(n_components=n_components)
        pls.fit(X, Y)
        assert_equal(n_components, pls.y_scores_.shape[1])


def test_univariate_pls_regression():
    # Ensure 1d Y is correctly interpreted
    d = load_linnerud()
    X = d.data
    Y = d.target

    clf = pls_.PLSRegression()
    # Compare 1d to column vector
    model1 = clf.fit(X, Y[:, 0]).coef_
    model2 = clf.fit(X, Y[:, :1]).coef_
    assert_array_almost_equal(model1, model2)


def test_predict_transform_copy():
    # check that the "copy" keyword works
    d = load_linnerud()
    X = d.data
    Y = d.target
    clf = pls_.PLSCanonical()
    X_copy = X.copy()
    Y_copy = Y.copy()
    clf.fit(X, Y)
    # check that results are identical with copy
    assert_array_almost_equal(clf.predict(X), clf.predict(X.copy(), copy=False))
    assert_array_almost_equal(clf.transform(X), clf.transform(X.copy(), copy=False))

    # check also if passing Y
    assert_array_almost_equal(clf.transform(X, Y),
                              clf.transform(X.copy(), Y.copy(), copy=False))
    # check that copy doesn't destroy
    # we do want to check exact equality here
    assert_array_equal(X_copy, X)
    assert_array_equal(Y_copy, Y)
    # also check that mean wasn't zero before (to make sure we didn't touch it)
    assert np.all(X.mean(axis=0) != 0)


def test_scale_and_stability():
    # We test scale=True parameter
    # This allows to check numerical stability over platforms as well

    d = load_linnerud()
    X1 = d.data
    Y1 = d.target
    # causes X[:, -1].std() to be zero
    X1[:, -1] = 1.0

    # From bug #2821
    # Test with X2, T2 s.t. clf.x_score[:, 1] == 0, clf.y_score[:, 1] == 0
    # This test robustness of algorithm when dealing with value close to 0
    X2 = np.array([[0., 0., 1.],
                   [1., 0., 0.],
                   [2., 2., 2.],
                   [3., 5., 4.]])
    Y2 = np.array([[0.1, -0.2],
                   [0.9, 1.1],
                   [6.2, 5.9],
                   [11.9, 12.3]])

    for (X, Y) in [(X1, Y1), (X2, Y2)]:
        X_std = X.std(axis=0, ddof=1)
        X_std[X_std == 0] = 1
        Y_std = Y.std(axis=0, ddof=1)
        Y_std[Y_std == 0] = 1

        X_s = (X - X.mean(axis=0)) / X_std
        Y_s = (Y - Y.mean(axis=0)) / Y_std

        for clf in [CCA(), pls_.PLSCanonical(), pls_.PLSRegression(),
                    pls_.PLSSVD()]:
            clf.set_params(scale=True)
            X_score, Y_score = clf.fit_transform(X, Y)
            clf.set_params(scale=False)
            X_s_score, Y_s_score = clf.fit_transform(X_s, Y_s)
            assert_array_almost_equal(X_s_score, X_score, decimal=4)
            assert_array_almost_equal(Y_s_score, Y_score, decimal=4)
            # Scaling should be idempotent
            clf.set_params(scale=True)
            X_score, Y_score = clf.fit_transform(X_s, Y_s)
            assert_array_almost_equal(X_s_score, X_score, decimal=4)
            assert_array_almost_equal(Y_s_score, Y_score, decimal=4)


def test_pls_errors():
    d = load_linnerud()
    X = d.data
    Y = d.target
    for clf in [pls_.PLSCanonical(), pls_.PLSRegression(),
                pls_.PLSSVD()]:
        clf.n_components = 4
        assert_raise_message(ValueError, "Invalid number of components",
                             clf.fit, X, Y)


def test_pls_scaling():
    # sanity check for scale=True
    n_samples = 1000
    n_targets = 5
    n_features = 10

    rng = check_random_state(0)

    Q = rng.randn(n_targets, n_features)
    Y = rng.randn(n_samples, n_targets)
    X = np.dot(Y, Q) + 2 * rng.randn(n_samples, n_features) + 1
    X *= 1000
    X_scaled = StandardScaler().fit_transform(X)

    pls = pls_.PLSRegression(n_components=5, scale=True)

    pls.fit(X, Y)
    score = pls.score(X, Y)

    pls.fit(X_scaled, Y)
    score_scaled = pls.score(X_scaled, Y)

    assert_approx_equal(score, score_scaled)
