import numpy as np
from numpy import dot
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.datasets import load_linnerud

import sklearn.NIPALS as pls
from math import log
from sklearn.pls import PLSRegression
from sklearn.pls import PLSCanonical
from sklearn.pls import CCA
from sklearn.pls import PLSSVD
from sklearn.pls import _center_scale_xy

def test_NIPALS():

    d = load_linnerud()
    X = d.data
    Y = d.target
    tol = 5e-12
    Xorig = X.copy()
    Yorig = Y.copy()

    # 1) Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
    # ===========================================================
    # Compare 2 algo.: nipals vs. svd
    # ------------------------------
    pls_bynipals = PLSCanonical(n_components = X.shape[1])
    pls_bynipals.fit(X, Y)
    pls_bysvd = PLSCanonical(algorithm = "svd", n_components = X.shape[1])
    pls_bysvd.fit(X, Y)
    pls_byNIPALS = pls.PLSC(num_comp=X.shape[1], tolerance=tol, max_iter=1000)
    pls_byNIPALS.fit(X, Y)

    # check equalities of loading (up to the sign of the second column)
    pls_bynipals.x_loadings_, pls_bysvd.x_loadings_ = \
            pls.direct(pls_bynipals.x_loadings_, pls_bysvd.x_loadings_, compare = True)
    assert_array_almost_equal(
        pls_bynipals.x_loadings_, pls_bysvd.x_loadings_, decimal=5,
        err_msg="NIPALS and svd implementations lead to different X loadings")
    pls_bynipals.x_loadings_, pls_byNIPALS.P = \
            pls.direct(pls_bynipals.x_loadings_, pls_byNIPALS.P, compare = True)
    assert_array_almost_equal(
        pls_bynipals.x_loadings_, pls_byNIPALS.P, decimal=5,
        err_msg="The two NIPALS implementations lead to different X loadings")
    print "Testing equality of X loadings ... OK!"

    pls_bynipals.y_loadings_, pls_bysvd.y_loadings_ = \
            pls.direct(pls_bynipals.y_loadings_, pls_bysvd.y_loadings_, compare = True)
    assert_array_almost_equal(
        pls_bynipals.y_loadings_, pls_bysvd.y_loadings_, decimal = 5,
        err_msg="nipals and svd implementation lead to different Y loadings")
    pls_bynipals.y_loadings_, pls_byNIPALS.Q = \
            pls.direct(pls_bynipals.y_loadings_, pls_byNIPALS.Q, compare = True)
    assert_array_almost_equal(
        pls_bynipals.y_loadings_, pls_byNIPALS.Q, decimal = 5,
        err_msg="The two NIPALS implementations lead to different Y loadings")
    print "Testing equality of Y loadings ... OK!"

    # Check PLS properties (with n_components=X.shape[1])
    # ---------------------------------------------------
    plsca = PLSCanonical(n_components=X.shape[1])
    plsca.fit(X, Y)
    pls_byNIPALS = pls.PLSC(num_comp=X.shape[1], tolerance=tol, max_iter=1000)
    pls_byNIPALS.fit(X, Y)

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
    check_ortho(Wx, "X weights are not orthogonal")
    check_ortho(Wy, "Y weights are not orthogonal")
    check_ortho(pls_byNIPALS.W, "X weights are not orthogonal")
    check_ortho(pls_byNIPALS.C, "Y weights are not orthogonal")
    print "Testing orthogonality of weights ... OK!"

    # Orthogonality of latent scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(T, "X scores are not orthogonal")
    check_ortho(U, "Y scores are not orthogonal")
    check_ortho(pls_byNIPALS.T, "X scores are not orthogonal")
    check_ortho(pls_byNIPALS.U, "Y scores are not orthogonal")
    print "Testing orthogonality of scores ... OK!"


    # Check X = TP' and Y = UQ' (with (p == q) components)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # center scale X, Y
    Xc, Yc, x_mean, y_mean, x_std, y_std = \
        _center_scale_xy(X.copy(), Y.copy(), scale=True)
    assert_array_almost_equal(Xc, dot(T, P.T), err_msg="X != TP'")
    assert_array_almost_equal(Yc, dot(U, Q.T), err_msg="Y != UQ'")
    print "Testing equality of matriecs and their models ... OK!"

    Xc, mX = pls.center(X, return_means = True)
    Xc, sX = pls.scale(Xc, return_stds = True)
    Yc, mY = pls.center(Y, return_means = True)
    Yc, sY = pls.scale(Yc, return_stds = True)

    assert_array_almost_equal(Xc, dot(pls_byNIPALS.T, pls_byNIPALS.P.T), err_msg="X != TP'")
    assert_array_almost_equal(Yc, dot(pls_byNIPALS.U, pls_byNIPALS.Q.T), err_msg="Y != UQ'")
    print "Testing equality of matriecs and their models ... OK!"


    # Check that rotations on training data lead to scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plsca = PLSCanonical(n_components=X.shape[1])
    plsca.fit(X, Y)
    pls_byNIPALS = pls.PLSC(num_comp=X.shape[1], tolerance=tol, max_iter=1000)
    pls_byNIPALS.fit(X, Y)

    Xr = plsca.transform(X)
    assert_array_almost_equal(Xr, plsca.x_scores_,
                              err_msg="Rotation of X failed")
    Xr = pls_byNIPALS.transform(X)
    assert_array_almost_equal(Xr, pls_byNIPALS.T,
                              err_msg="Rotation of X failed")
    print "Testing equality of computed X scores and the transform ... OK!"
    Xr, Yr = plsca.transform(X, Y)
    assert_array_almost_equal(Xr, plsca.x_scores_,
                              err_msg="Rotation of X failed")
    assert_array_almost_equal(Yr, plsca.y_scores_,
                              err_msg="Rotation of Y failed")
    print "Testing equality of computed X and Y scores and the transform ... OK!"


    # "Non regression test" on canonical PLS
    # --------------------------------------
    # The results were checked against the R-package plspm
    pls_ca = PLSCanonical(n_components=X.shape[1])
    pls_ca.fit(X, Y)
    pls_byNIPALS = pls.PLSC(num_comp=X.shape[1], tolerance=tol, max_iter=1000)
    pls_byNIPALS.fit(X, Y)

    x_weights = np.array(
        [[-0.61330704,  0.25616119, -0.74715187],
         [-0.74697144,  0.11930791,  0.65406368],
         [-0.25668686, -0.95924297, -0.11817271]])
    x_weights, pls_ca.x_weights_ = \
            pls.direct(x_weights, pls_ca.x_weights_, compare = True)
    assert_array_almost_equal(x_weights, pls_ca.x_weights_, decimal = 5)
    x_weights, pls_byNIPALS.W = \
            pls.direct(x_weights, pls_byNIPALS.W, compare = True)
    assert_array_almost_equal(x_weights, pls_byNIPALS.W, decimal = 5)
    print "Testing equality of X weights ... OK!"

    x_rotations = np.array(
        [[-0.61330704,  0.41591889, -0.62297525],
         [-0.74697144,  0.31388326,  0.77368233],
         [-0.25668686, -0.89237972, -0.24121788]])
    x_rotations, pls_ca.x_rotations_ = \
            pls.direct(x_rotations, pls_ca.x_rotations_, compare = True)
    assert_array_almost_equal(x_rotations, pls_ca.x_rotations_, decimal = 5)
    x_rotations, pls_byNIPALS.Ws = \
            pls.direct(x_rotations, pls_byNIPALS.Ws, compare = True)
    assert_array_almost_equal(x_rotations, pls_byNIPALS.Ws, decimal = 5)
    print "Testing equality of X loadings weights ... OK!"

    y_weights = np.array(
        [[+0.58989127,  0.7890047,   0.1717553],
         [+0.77134053, -0.61351791,  0.16920272],
         [-0.23887670, -0.03267062,  0.97050016]])
    y_weights, pls_ca.y_weights_ = \
            pls.direct(y_weights, pls_ca.y_weights_, compare = True)
    assert_array_almost_equal(y_weights, pls_ca.y_weights_, decimal = 5)
    y_weights, pls_byNIPALS.C = \
            pls.direct(y_weights, pls_byNIPALS.C, compare = True)
    assert_array_almost_equal(y_weights, pls_byNIPALS.C, decimal = 5)
    print "Testing equality of Y weights ... OK!"

    y_rotations = np.array(
        [[+0.58989127,  0.7168115,  0.30665872],
         [+0.77134053, -0.70791757,  0.19786539],
         [-0.23887670, -0.00343595,  0.94162826]])
    pls_ca.y_rotations_, y_rotations = \
            pls.direct(pls_ca.y_rotations_, y_rotations, compare = True)
    assert_array_almost_equal(pls_ca.y_rotations_, y_rotations)
    y_rotations, pls_byNIPALS.Cs = \
            pls.direct(y_rotations, pls_byNIPALS.Cs, compare = True)
    assert_array_almost_equal(y_rotations, pls_byNIPALS.Cs, decimal = 5)
    print "Testing equality of Y loadings weights ... OK!"

    assert_array_almost_equal(X, Xorig, decimal = 5, err_msg = "X and Xorig are not equal!!")
    assert_array_almost_equal(Y, Yorig, decimal = 5, err_msg = "Y and Yorig are not equal!!")


    # 2) Regression PLS (PLS2): "Non regression test"
    # ===============================================
    # The results were checked against the R-packages plspm, misOmics and pls
    pls_2 = PLSRegression(n_components=X.shape[1])
    pls_2.fit(X, Y)

    pls_NIPALS = pls.PLSR(num_comp = X.shape[1],
                          center = True, scale = True,
                          tolerance=tol, max_iter=1000)
    pls_NIPALS.fit(X, Y)

    x_weights = np.array(
        [[-0.61330704, -0.00443647,  0.78983213],
         [-0.74697144, -0.32172099, -0.58183269],
         [-0.25668686,  0.94682413, -0.19399983]])
    x_weights, pls_NIPALS.W = pls.direct(x_weights, pls_NIPALS.W, compare = True)
    assert_array_almost_equal(pls_NIPALS.W, x_weights, decimal=5,
            err_msg="sklearn.NIPALS.PLSR and sklearn.pls.PLSRegression " \
                    "implementations lead to different X weights")
    print "Comparing X weights of sklearn.NIPALS.PLSR and " \
            "sklearn.pls.PLSRegression ... OK!"

    x_loadings = np.array(
        [[-0.61470416, -0.24574278,  0.78983213],
         [-0.65625755, -0.14396183, -0.58183269],
         [-0.51733059,  1.00609417, -0.19399983]])
    x_loadings, pls_NIPALS.P = pls.direct(x_loadings, pls_NIPALS.P, compare = True)
    assert_array_almost_equal(pls_NIPALS.P, x_loadings, decimal=5,
            err_msg="sklearn.NIPALS.PLSR and sklearn.pls.PLSRegression " \
                    "implementations lead to different X loadings")
    print "Comparing X loadings of sklearn.NIPALS.PLSR and " \
            "sklearn.pls.PLSRegression ... OK!"

    y_weights = np.array(
        [[+0.32456184,  0.29892183,  0.20316322],
         [+0.42439636,  0.61970543,  0.19320542],
         [-0.13143144, -0.26348971, -0.17092916]])
    y_weights, pls_NIPALS.C = pls.direct(y_weights, pls_NIPALS.C, compare = True)
    assert_array_almost_equal(pls_NIPALS.C, y_weights, decimal=5,
            err_msg="sklearn.NIPALS.PLSR and sklearn.pls.PLSRegression " \
                    "implementations lead to different Y weights")
    print "Comparing Y weights of sklearn.NIPALS.PLSR and " \
            "sklearn.pls.PLSRegression ... OK!"

    X_, m = pls.center(X, return_means = True)
    X_, s = pls.scale(X_, return_stds = True, centered = True)
    t1 = dot(X_, x_weights[:,[0]])
    t1, pls_NIPALS.T[:,[0]] = pls.direct(t1, pls_NIPALS.T[:,[0]], compare = True)
    assert_array_almost_equal(t1, pls_NIPALS.T[:,[0]], decimal=5,
            err_msg="sklearn.NIPALS.PLSR and sklearn.pls.PLSRegression " \
                    "implementations lead to different X scores")
    print "Comparing scores of sklearn.NIPALS.PLSR and " \
            "sklearn.pls.PLSRegression ... OK!"

    y_loadings = np.array(
        [[+0.32456184,  0.29892183,  0.20316322],
         [+0.42439636,  0.61970543,  0.19320542],
         [-0.13143144, -0.26348971, -0.17092916]])
    y_loadings, pls_NIPALS.C = pls.direct(y_loadings, pls_NIPALS.C, compare = True)
    assert_array_almost_equal(pls_NIPALS.C, y_loadings)


    # 3) Another non-regression test of Canonical PLS on random dataset
    # =================================================================
    # The results were checked against the R-package plspm
    #
    # Warning! This example is not stable, and the reference weights have
    # not converged properly!
    #
    n = 500
    p_noise = 10
    q_noise = 5
    # 2 latents vars:
    np.random.seed(11)
    l1 = np.random.normal(size=n)
    l2 = np.random.normal(size=n)
    latents = np.array([l1, l1, l2, l2]).T
    X = latents + np.random.normal(size=4 * n).reshape((n, 4))
    Y = latents + np.random.normal(size=4 * n).reshape((n, 4))
    X = np.concatenate(
        (X, np.random.normal(size=p_noise * n).reshape(n, p_noise)), axis=1)
    Y = np.concatenate(
        (Y, np.random.normal(size=q_noise * n).reshape(n, q_noise)), axis=1)
    np.random.seed(None)
    x_weights = np.array(
        [[ 0.65803719,  0.19197924,  0.21769083],
         [ 0.7009113,   0.13303969, -0.15376699],
         [ 0.13528197, -0.68636408,  0.13856546],
         [ 0.16854574, -0.66788088, -0.12485304],
         [-0.03232333, -0.04189855,  0.40690153],
         [ 0.1148816,  -0.09643158,  0.1613305 ],
         [ 0.04792138, -0.02384992,  0.17175319],
         [-0.06781,    -0.01666137, -0.18556747],
         [-0.00266945, -0.00160224,  0.11893098],
         [-0.00849528, -0.07706095,  0.1570547 ],
         [-0.00949471, -0.02964127,  0.34657036],
         [-0.03572177,  0.0945091,   0.3414855 ],
         [ 0.05584937, -0.02028961, -0.57682568],
         [ 0.05744254, -0.01482333, -0.17431274]])
    tols = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]#, 5e-7, 5e-8, 5e-9, 5e-10, 5e-11]
    for i in tols:
        pls_ca = PLSCanonical(n_components=3, max_iter=1000, tol=i)
        pls_ca.fit(X, Y)

        x_weights, pls_ca.x_weights_ = pls.direct(x_weights, pls_ca.x_weights_, compare = True)
        print "tolerance: "+str(i).rjust(6)+", error:", np.max(pls_ca.x_weights_ - x_weights)

    assert_array_almost_equal(pls_ca.x_weights_, x_weights, decimal = 4,
            err_msg="sklearn.pls.PLSCanonical does not give the same " \
                    "X weights as the reference model")
    print "Comparing X weights of sklearn.pls.PLSCanonical ... OK!"

    for i in tols:
        pls_NIPALS = pls.PLSC(num_comp=3, tolerance=i, max_iter=1000)
        pls_NIPALS.fit(X, Y)

        x_weights, pls_NIPALS.W = pls.direct(x_weights, pls_NIPALS.W, compare = True)
        print "tolerance: "+str(i).rjust(6)+", error:", np.max(x_weights - pls_NIPALS.W)

    assert_array_almost_equal(pls_NIPALS.W, x_weights, decimal = 4,
            err_msg="sklearn.NIPALS.PLSC does not give the same " \
                    "X weights as the reference model")
    print "Comparing X weights of sklearn.NIPALS.PLSC ... OK! "

    x_loadings = np.array(
        [[ 0.65649254,  0.1847647,   0.15270699],
         [ 0.67554234,  0.15237508, -0.09182247],
         [ 0.19219925, -0.67750975,  0.08673128],
         [ 0.2133631,  -0.67034809, -0.08835483],
         [-0.03178912, -0.06668336,  0.43395268],
         [ 0.15684588, -0.13350241,  0.20578984],
         [ 0.03337736, -0.03807306,  0.09871553],
         [-0.06199844,  0.01559854, -0.1881785 ],
         [ 0.00406146, -0.00587025,  0.16413253],
         [-0.00374239, -0.05848466,  0.19140336],
         [ 0.00139214, -0.01033161,  0.32239136],
         [-0.05292828,  0.0953533,   0.31916881],
         [ 0.04031924, -0.01961045, -0.65174036],
         [ 0.06172484, -0.06597366, -0.1244497]])
    pls_ca.x_loadings_, x_loadings = pls.direct(pls_ca.x_loadings_, x_loadings, compare = True)
    assert_array_almost_equal(pls_ca.x_loadings_, x_loadings, decimal = 4,
            err_msg="sklearn.pls.PLSCanonical does not give the same " \
                    "X loadings as the reference model")
    print "Comparing X loadings of sklearn.pls.PLSCanonical ... OK!"

    pls_NIPALS.P, x_loadings = pls.direct(pls_NIPALS.P, x_loadings, compare = True)
    assert_array_almost_equal(pls_NIPALS.P, x_loadings, decimal = 4,
            err_msg="sklearn.NIPALS.PLSC does not give the same " \
                    "loadings as the reference model")
    print "Comparing X loadings of sklearn.NIPALS.PLSC ... OK! "

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
    pls_ca.y_weights_, y_weights = pls.direct(pls_ca.y_weights_, y_weights, compare = True)
    assert_array_almost_equal(pls_ca.y_weights_, y_weights, decimal = 4,
            err_msg="sklearn.pls.PLSCanonical does not give the same " \
                    "Y weights as the reference model")
    print "Comparing Y weights of sklearn.pls.PLSCanonical ... OK!"

    pls_NIPALS.C, y_weights = pls.direct(pls_NIPALS.C, y_weights, compare = True)
    assert_array_almost_equal(pls_NIPALS.C, y_weights, decimal = 4,
            err_msg="sklearn.NIPALS.PLSC does not give the same " \
                    "loadings as the reference model")
    print "Comparing Y weights of sklearn.NIPALS.PLSC ... OK! "

    y_loadings = np.array(
        [[0.68568625,   0.1674376,   0.0969508 ],
         [0.68782064,   0.20375837, -0.1164448 ],
         [0.11712173,  -0.68046903,  0.12001505],
         [0.17860457,  -0.6798319,  -0.05089681],
         [0.06265739,  -0.0277703,   0.74729584],
         [0.0914178,    0.00403751, -0.5135078 ],
         [-0.02196918, -0.01377169,  0.09564505],
         [-0.03288952,  0.09039729,  0.31858973],
         [0.04287624,   0.05254676,  0.27836841]])
    pls_ca.y_loadings_, y_loadings = pls.direct(pls_ca.y_loadings_, y_loadings, compare = True)
    assert_array_almost_equal(pls_ca.y_loadings_, y_loadings, decimal = 4,
            err_msg="sklearn.pls.PLSCanonical does not give the same " \
                    "Y loadings as the reference model")
    print "Comparing Y loadings of sklearn.pls.PLSCanonical ... OK!"

    pls_NIPALS.Q, y_loadings = pls.direct(pls_NIPALS.Q, y_loadings, compare = True)
    assert_array_almost_equal(pls_NIPALS.Q, y_loadings, decimal = 4,
            err_msg="sklearn.NIPALS.PLSC does not give the same " \
                    "Y loadings as the reference model")
    print "Comparing Y loadings of sklearn.NIPALS.PLSC ... OK!"

    # Orthogonality of weights
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(pls_ca.x_weights_, "X weights are not orthogonal in sklearn.pls.PLSCanonical")
    check_ortho(pls_ca.y_weights_, "Y weights are not orthogonal in sklearn.pls.PLSCanonical")
    print "Confirming orthogonality of weights in sklearn.pls.PLSCanonical ... OK!"
    check_ortho(pls_NIPALS.W, "X weights are not orthogonal in sklearn.NIPALS.PLSC")
    check_ortho(pls_NIPALS.C, "Y weights are not orthogonal in sklearn.NIPALS.PLSC")
    print "Confirming orthogonality of weights in sklearn.NIPALS.PLSC ... OK!"

    # Orthogonality of latent scores
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    check_ortho(pls_ca.x_scores_, "X scores are not orthogonal in sklearn.pls.PLSCanonical")
    check_ortho(pls_ca.y_scores_, "Y scores are not orthogonal in sklearn.pls.PLSCanonical")
    print "Confirming orthogonality of scores in sklearn.pls.PLSCanonical ... OK!"
    check_ortho(pls_NIPALS.T, "X scores are not orthogonal in sklearn.NIPALS.PLSC")
    check_ortho(pls_NIPALS.U, "Y scores are not orthogonal in sklearn.NIPALS.PLSC")
    print "Confirming orthogonality of scores in sklearn.NIPALS.PLSC ... OK!"


    # Compare sparse sklearn.NIPALS.PLSR and sklearn.pls.PLSRegression

    d = load_linnerud()
    X = np.asarray(d.data)
    Y = d.target
    num_comp = 3
    tol = 5e-12

    for st in [0.1, 0.01, 0.001, 0.0001, 0]:

        plsr = pls.PLSR(num_comp = num_comp, center = True, scale = True,
                        tolerance=tol, max_iter=1000, soft_threshold = st)
        plsr.fit(X, Y)
        Yhat = plsr.predict(X)

        pls2 = PLSRegression(n_components=num_comp, scale=True,
                     max_iter=1000, tol=tol, copy=True)
        pls2.fit(X, Y)
        Yhat_ = pls2.predict(X)

        SSY     = np.sum(Y**2)
        SSYdiff = np.sum((Y-Yhat)**2)
#        print np.sum(1 - (SSYdiff/SSY))
        SSY     = np.sum(Y**2)
        SSYdiff = np.sum((Y-Yhat_)**2)
#        print np.sum(1 - (SSYdiff/SSY))

        if st < tol:
            num_decimals = 5
        else:
            num_decimals = int(log(1./st, 10) + 0.5)
        assert_array_almost_equal(Yhat, Yhat_, decimal=num_decimals-2,
                err_msg="NIPALS SVD and numpy.linalg.svd implementations " \
                "lead to different loadings")
        print "Comparing loadings of PLSR and sklearn.pls ... OK!" \
                " (err=%.4f, threshold=%0.4f)" % (np.sum((Yhat-Yhat_)**2), st)


    # Compare PLSR in sklearn.NIPALS and sklearn.pls

    d = load_linnerud()
    X = np.asarray(d.data)
    Y = d.target
    num_comp = 3
    tol = 5e-12

    plsr = pls.PLSR(num_comp = num_comp, center = True, scale = True,
                    tolerance=tol, max_iter=1000, soft_threshold = 0)
    plsr.fit(X, Y)
    Yhat = plsr.predict(X)

    Xorig = X.copy()
    Yorig = Y.copy()
    X, mX = pls.center(X, return_means = True)
    X, sX = pls.scale(X, return_stds = True, centered = True)
    Y, mY = pls.center(Y, return_means = True)
    Y, sY = pls.scale(Y, return_stds = True, centered = True)
    W = np.zeros((X.shape[1],num_comp))
    P = np.zeros((X.shape[1],num_comp))
    C = np.zeros((Y.shape[1],num_comp))
    Xuv = X.copy()
#    Yuv = Y.copy()

    XY = dot(X.T, Y)
    W_, S_, C_ = np.linalg.svd(XY)
    w1 = W_[:,[0]]
    W[:,0] = w1.ravel()
    t1 = dot(X, w1)
    p1 = dot(X.T, t1) / dot(t1.T, t1)
    P[:,0] = p1.ravel()
    c1 = dot(Y.T, t1) / dot(t1.T, t1)
    C[:,0] = c1.ravel()
    X -= dot(t1,p1.T)

    if num_comp >= 2:
        XY = dot(X.T, Y)
        W_, S_, C_ = np.linalg.svd(XY)
        w2 = W_[:,[0]]
        W[:,1] = w2.ravel()
        t2 = dot(X, w2)
        p2 = dot(X.T, t2) / dot(t2.T, t2)
        P[:,1] = p2.ravel()
        c2 = dot(Y.T, t2) / dot(t2.T, t2)
        C[:,1] = c2.ravel()
        X -= dot(t2,p2.T)

    if num_comp >= 3:
        XY = dot(X.T, Y)
        W_, S_, C_ = np.linalg.svd(XY)
        w3 = W_[:,[0]]
        W[:,2] = w3.ravel()
        t3 = dot(X, w3)
        p3 = dot(X.T, t3) / dot(t3.T, t3)
        P[:,2] = p3.ravel()
        c3 = dot(Y.T, t3) / dot(t3.T, t3)
        C[:,2] = c3.ravel()
        X -= dot(t3,p3.T)

    Ws  = dot(W, np.linalg.inv(dot(P.T, W)))
    B   = dot(Ws, C.T)
    Yhat_ = (dot(Xuv, B) * sY) + mY

    pls2 = PLSRegression(n_components=num_comp, scale=True,
                 max_iter=1000, tol=tol, copy=True)
    pls2.fit(Xorig, Yorig)
    Yhat__ = pls2.predict(Xorig)

    SSY     = np.sum(Yorig**2)
    SSYdiff = np.sum((Yorig-Yhat)**2)
    print "sklearn.NIPALS: Comparing original and predicted Y ... OK!"\
            " (R2Yhat = %.4f)" % (1 - (SSYdiff / SSY))
    SSYdiff = np.sum((Yorig-Yhat_)**2)
    print "Manual        : Comparing original and predicted Y ... OK!"\
            " (R2Yhat = %.4f)" % (1 - (SSYdiff / SSY))
    SSYdiff = np.sum((Yorig-Yhat__)**2)
    print "sklearn.pls   : Comparing original and predicted Y ... OK!"\
            " (R2Yhat = %.4f)" % (1 - (SSYdiff / SSY))
    assert_array_almost_equal(Yhat__, Yhat, decimal=5,
            err_msg="NIPALS and pls implementations lead to different" \
            " predictions")


    # Compare SVD with and without sparsity constraint to numpy.linalg.svd

    Xtr = np.random.rand(6,6)
    #Xte = np.random.rand(2,6)
    num_comp = 3
    tol = 5e-12

    Xtr, m = pls.center(Xtr, return_means = True)
    Xtr, s = pls.scale(Xtr, return_stds = True)
    #Xte = (Xte - m) / s

    for st in [0.1, 0.01, 0.001, 0.0001, 0]:
        # NIPALS.SVD
        svd = pls.SVD(num_comp = num_comp,
                      tolerance=tol, max_iter=1000, soft_threshold = st)
        svd.fit(Xtr)
#        svd.V, svd.U = pls.direct(svd.V, svd.U)

        # numpy.lialg.svd
        U, S, V = np.linalg.svd(Xtr)
        V = V.T
        S = np.diag(S)
        U = U[:,0:num_comp]
        S = S[:,0:num_comp]
        V = V[:,0:num_comp]
#        V, U = pls.direct(V, U)
        #SVDte = dot(Xte, V)

        if st < tol:
            num_decimals = 5
        else:
            num_decimals = int(log(1./st, 10) + 0.5)
        svd.V, V = pls.direct(svd.V, V, compare = True)
        assert_array_almost_equal(svd.V, V, decimal=num_decimals-2,
                err_msg="sklearn.NIPALS.SVD and numpy.linalg.svd implementations " \
                "lead to different loadings")
        print "Comparing loadings of sklearn.NIPALS.SVD and numpy.linalg.svd... OK!"\
                " (diff=%.4f, threshold=%0.4f)" % (np.max(np.abs(V - svd.V)), st)


    # Compare PCA with sparsity constraint to numpy.linalg.svd

    Xtr = np.random.rand(5,5)
    num_comp = 5
    tol = 5e-9

    Xtr, m = pls.center(Xtr, return_means = True)
    Xtr, s = pls.scale(Xtr, return_stds = True)

    for st in [0.1, 0.01, 0.001, 0.0001, 0]:
        pca = pls.PCA(center = False, scale = False, num_comp = num_comp,
                  tolerance=tol, max_iter=500, soft_threshold = st)
        pca.fit(Xtr)
        Tte = pca.transform(Xtr)
        U, S, V = np.linalg.svd(Xtr)
        V = V.T
        US = dot(U,np.diag(S))
        US = US[:,0:num_comp]
        V  = V[:,0:num_comp]
    
        err = np.max(np.abs(Xtr - dot(pca.T,pca.P.T)))

        if st < tol:
            num_decimals = 5
        else:
            num_decimals = int(log(1./st, 10) + 0.5)
        assert_array_almost_equal(Xtr, dot(pca.T,pca.P.T), decimal = num_decimals-1,
                err_msg="Model does not equal the matrices")
        print "PCA: Testing equality of model and matrix ... OK! (err=%.5f, threshold=%.4f)" % (err, st)


    return


    # Compare PCA without the sparsity constraint to numpy.linalg.svd

    Xtr = np.random.rand(50,50)
    Xte = np.random.rand(20,50)
    num_comp = 3

    Xtr, m = pls.center(Xtr, return_means = True)
    Xtr, s = pls.scale(Xtr, return_stds = True)
    Xte = (Xte - m) / s

    pca = pls.PCA(center = False, scale = False, num_comp = num_comp,
                  tolerance=5e-12, max_iter=1000)
    pca.fit(Xtr)
    pca.P, pca.T = pls.direct(pca.P, pca.T)
    Tte = pca.transform(Xte)

    U, S, V = np.linalg.svd(Xtr)
    V = V.T
    US = dot(U,np.diag(S))
    US = US[:,0:num_comp]
    V  = V[:,0:num_comp]
    V, US = pls.direct(V, US)
    SVDte = dot(Xte, V)

    assert_array_almost_equal(pca.P, V, decimal=2, err_msg="NIPALS PCA and "
            "numpy.linalg.svd implementations lead to different loadings")
    print "Comparing loadings of NIPALS PCA and numpy.linalg.svd... OK! "\
          "(max diff = %.4f)" % np.max(np.abs(V - pca.P))

    assert_array_almost_equal(pca.T, US, decimal=2, err_msg="NIPALS PCA and "
            "numpy.linalg.svd implementations lead to different scores")
    print "Comparing scores of NIPALS PCA and numpy.linalg.svd...   OK! "\
          "(max diff = %.4f)" % np.max(np.abs(US - pca.T))

    assert_array_almost_equal(Tte, SVDte, decimal=2, err_msg="NIPALS PCA and "
            "numpy.linalg.svd implementations lead to different scores")
    print "Comparing test set of NIPALS PCA and numpy.linalg.svd... OK! "\
          "(max diff = %.4f)" % np.max(np.abs(Tte - SVDte))


    # Test CCA
    
    num_comp = 3
    tol = 5e-9
    
    d = load_linnerud()
    X = np.asarray(d.data)
    Y = d.target

    #print X
    #print Y



def test_scale():

    d = load_linnerud()
    X = d.data
    Y = d.target

    # causes X[:, -1].std() to be zero
    X[:, -1] = 1.0

    methods = [PLSCanonical(), PLSRegression(), CCA(), PLSSVD(),
               pls.PCA(), pls.SVD(), pls.PLSR(), pls.PLSC()]
    names   = ["PLSCanonical", "PLSRegression", "CCA", "PLSSVD",
               "pls.PCA", "pls.SVD", "pls.PLSR", "pls.PLSC"]
    for i in xrange(len(methods)):
        clf = methods[i]
        print "Testing scale of "+names[i]
#        clf.set_params(scale=True)
        clf.scale = True
        clf.fit(X, Y)


if __name__ == "__main__":
    test_NIPALS()
    test_scale()
