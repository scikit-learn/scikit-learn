import numpy as np
from numpy import dot
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.datasets import load_linnerud

import sklearn.NIPALS as pls
from math import log
from numpy import dot

def test_NIPALS():

    pass
#    d = load_linnerud()
#    X = d.data
#    Y = d.target
#    # 1) Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
#    # ===========================================================
#    # Compare 2 algo.: nipals vs. svd
#    # ------------------------------
#    pls_bynipals = pls.PLSCanonical(n_components=X.shape[1])
#    pls_bynipals.fit(X, Y)
#    pls_bysvd = pls.PLSCanonical(algorithm="svd", n_components=X.shape[1])
#    pls_bysvd.fit(X, Y)
#    # check equalities of loading (up to the sign of the second column)
#    assert_array_almost_equal(
#        pls_bynipals.x_loadings_,
#        np.multiply(pls_bysvd.x_loadings_, np.array([1, -1, 1])), decimal=5,
#        err_msg="nipals and svd implementation lead to different x loadings")
#
#    assert_array_almost_equal(
#        pls_bynipals.y_loadings_,
#        np.multiply(pls_bysvd.y_loadings_, np.array([1, -1, 1])), decimal=5,
#        err_msg="nipals and svd implementation lead to different y loadings")
#
#    # Check PLS properties (with n_components=X.shape[1])
#    # ---------------------------------------------------
#    plsca = pls.PLSCanonical(n_components=X.shape[1])
#    plsca.fit(X, Y)
#    T = plsca.x_scores_
#    P = plsca.x_loadings_
#    Wx = plsca.x_weights_
#    U = plsca.y_scores_
#    Q = plsca.y_loadings_
#    Wy = plsca.y_weights_
#
#    def check_ortho(M, err_msg):
#        K = np.dot(M.T, M)
#        assert_array_almost_equal(K, np.diag(np.diag(K)), err_msg=err_msg)
#
#    # Orthogonality of weights
#    # ~~~~~~~~~~~~~~~~~~~~~~~~
#    check_ortho(Wx, "x weights are not orthogonal")
#    check_ortho(Wy, "y weights are not orthogonal")
#
#    # Orthogonality of latent scores
#    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    check_ortho(T, "x scores are not orthogonal")
#    check_ortho(U, "y scores are not orthogonal")
#
#    # Check X = TP' and Y = UQ' (with (p == q) components)
#    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    # center scale X, Y
#    Xc, Yc, x_mean, y_mean, x_std, y_std =\
#        pls._center_scale_xy(X.copy(), Y.copy(), scale=True)
#    assert_array_almost_equal(Xc, np.dot(T, P.T), err_msg="X != TP'")
#    assert_array_almost_equal(Yc, np.dot(U, Q.T), err_msg="Y != UQ'")
#
#    # Check that rotations on training data lead to scores
#    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    Xr = plsca.transform(X)
#    assert_array_almost_equal(Xr, plsca.x_scores_,
#                              err_msg="rotation on X failed")
#    Xr, Yr = plsca.transform(X, Y)
#    assert_array_almost_equal(Xr, plsca.x_scores_,
#                              err_msg="rotation on X failed")
#    assert_array_almost_equal(Yr, plsca.y_scores_,
#                              err_msg="rotation on Y failed")
#
#    # "Non regression test" on canonical PLS
#    # --------------------------------------
#    # The results were checked against the R-package plspm
#    pls_ca = pls.PLSCanonical(n_components=X.shape[1])
#    pls_ca.fit(X, Y)
#
#    x_weights = np.array(
#        [[-0.61330704,  0.25616119, -0.74715187],
#         [-0.74697144,  0.11930791,  0.65406368],
#         [-0.25668686, -0.95924297, -0.11817271]])
#    assert_array_almost_equal(pls_ca.x_weights_, x_weights)
#
#    x_rotations = np.array(
#        [[-0.61330704,  0.41591889, -0.62297525],
#         [-0.74697144,  0.31388326,  0.77368233],
#         [-0.25668686, -0.89237972, -0.24121788]])
#    assert_array_almost_equal(pls_ca.x_rotations_, x_rotations)
#
#    y_weights = np.array(
#        [[+0.58989127,  0.7890047,   0.1717553],
#         [+0.77134053, -0.61351791,  0.16920272],
#         [-0.23887670, -0.03267062,  0.97050016]])
#    assert_array_almost_equal(pls_ca.y_weights_, y_weights)
#
#    y_rotations = np.array(
#        [[+0.58989127,  0.7168115,  0.30665872],
#         [+0.77134053, -0.70791757,  0.19786539],
#         [-0.23887670, -0.00343595,  0.94162826]])
#    assert_array_almost_equal(pls_ca.y_rotations_, y_rotations)
#
#    # 2) Regression PLS (PLS2): "Non regression test"
#    # ===============================================
#    # The results were checked against the R-packages plspm, misOmics and pls
#    pls_2 = pls.PLSRegression(n_components=X.shape[1])
#    pls_2.fit(X, Y)
#
#    x_weights = np.array(
#        [[-0.61330704, -0.00443647,  0.78983213],
#         [-0.74697144, -0.32172099, -0.58183269],
#         [-0.25668686,  0.94682413, -0.19399983]])
#    assert_array_almost_equal(pls_2.x_weights_, x_weights)
#
#    x_loadings = np.array(
#        [[-0.61470416, -0.24574278,  0.78983213],
#         [-0.65625755, -0.14396183, -0.58183269],
#         [-0.51733059,  1.00609417, -0.19399983]])
#    assert_array_almost_equal(pls_2.x_loadings_, x_loadings)
#
#    y_weights = np.array(
#        [[+0.32456184,  0.29892183,  0.20316322],
#         [+0.42439636,  0.61970543,  0.19320542],
#         [-0.13143144, -0.26348971, -0.17092916]])
#    assert_array_almost_equal(pls_2.y_weights_, y_weights)
#
#    y_loadings = np.array(
#        [[+0.32456184,  0.29892183,  0.20316322],
#         [+0.42439636,  0.61970543,  0.19320542],
#         [-0.13143144, -0.26348971, -0.17092916]])
#    assert_array_almost_equal(pls_2.y_loadings_, y_loadings)
#
#    # 3) Another non-regression test of Canonical PLS on random dataset
#    # =================================================================
#    # The results were checked against the R-package plspm
#    n = 500
#    p_noise = 10
#    q_noise = 5
#    # 2 latents vars:
#    np.random.seed(11)
#    l1 = np.random.normal(size=n)
#    l2 = np.random.normal(size=n)
#    latents = np.array([l1, l1, l2, l2]).T
#    X = latents + np.random.normal(size=4 * n).reshape((n, 4))
#    Y = latents + np.random.normal(size=4 * n).reshape((n, 4))
#    X = np.concatenate(
#        (X, np.random.normal(size=p_noise * n).reshape(n, p_noise)), axis=1)
#    Y = np.concatenate(
#        (Y, np.random.normal(size=q_noise * n).reshape(n, q_noise)), axis=1)
#    np.random.seed(None)
#    pls_ca = pls.PLSCanonical(n_components=3)
#    pls_ca.fit(X, Y)
#
#    x_weights = np.array(
#        [[0.65803719,  0.19197924,  0.21769083],
#         [0.7009113,  0.13303969, -0.15376699],
#         [0.13528197, -0.68636408,  0.13856546],
#         [0.16854574, -0.66788088, -0.12485304],
#         [-0.03232333, -0.04189855,  0.40690153],
#         [0.1148816, -0.09643158,  0.1613305],
#         [0.04792138, -0.02384992,  0.17175319],
#         [-0.06781, -0.01666137, -0.18556747],
#         [-0.00266945, -0.00160224,  0.11893098],
#         [-0.00849528, -0.07706095,  0.1570547],
#         [-0.00949471, -0.02964127,  0.34657036],
#         [-0.03572177,  0.0945091,  0.3414855],
#         [0.05584937, -0.02028961, -0.57682568],
#         [0.05744254, -0.01482333, -0.17431274]])
#    assert_array_almost_equal(pls_ca.x_weights_, x_weights)
#
#    x_loadings = np.array(
#        [[0.65649254,  0.1847647,  0.15270699],
#         [0.67554234,  0.15237508, -0.09182247],
#         [0.19219925, -0.67750975,  0.08673128],
#         [0.2133631, -0.67034809, -0.08835483],
#         [-0.03178912, -0.06668336,  0.43395268],
#         [0.15684588, -0.13350241,  0.20578984],
#         [0.03337736, -0.03807306,  0.09871553],
#         [-0.06199844,  0.01559854, -0.1881785],
#         [0.00406146, -0.00587025,  0.16413253],
#         [-0.00374239, -0.05848466,  0.19140336],
#         [0.00139214, -0.01033161,  0.32239136],
#         [-0.05292828,  0.0953533,  0.31916881],
#         [0.04031924, -0.01961045, -0.65174036],
#         [0.06172484, -0.06597366, -0.1244497]])
#    assert_array_almost_equal(pls_ca.x_loadings_, x_loadings)
#
#    y_weights = np.array(
#        [[0.66101097,  0.18672553,  0.22826092],
#         [0.69347861,  0.18463471, -0.23995597],
#         [0.14462724, -0.66504085,  0.17082434],
#         [0.22247955, -0.6932605, -0.09832993],
#         [0.07035859,  0.00714283,  0.67810124],
#         [0.07765351, -0.0105204, -0.44108074],
#         [-0.00917056,  0.04322147,  0.10062478],
#         [-0.01909512,  0.06182718,  0.28830475],
#         [0.01756709,  0.04797666,  0.32225745]])
#    assert_array_almost_equal(pls_ca.y_weights_, y_weights)
#
#    y_loadings = np.array(
#        [[0.68568625,  0.1674376,  0.0969508],
#         [0.68782064,  0.20375837, -0.1164448],
#         [0.11712173, -0.68046903,  0.12001505],
#         [0.17860457, -0.6798319, -0.05089681],
#         [0.06265739, -0.0277703,  0.74729584],
#         [0.0914178,  0.00403751, -0.5135078],
#         [-0.02196918, -0.01377169,  0.09564505],
#         [-0.03288952,  0.09039729,  0.31858973],
#         [0.04287624,  0.05254676,  0.27836841]])
#    assert_array_almost_equal(pls_ca.y_loadings_, y_loadings)
#
#    # Orthogonality of weights
#    # ~~~~~~~~~~~~~~~~~~~~~~~~
#    check_ortho(pls_ca.x_weights_, "x weights are not orthogonal")
#    check_ortho(pls_ca.y_weights_, "y weights are not orthogonal")
#
#    # Orthogonality of latent scores
#    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    check_ortho(pls_ca.x_scores_, "x scores are not orthogonal")
#    check_ortho(pls_ca.y_scores_, "y scores are not orthogonal")
#
#
def test_scale():
    pass
#    d = load_linnerud()
#    X = d.data
#    Y = d.target
#
#    # causes X[:, -1].std() to be zero
#    X[:, -1] = 1.0
#
#    for clf in [pls.PLSCanonical(), pls.PLSRegression(), pls.CCA(),
#                pls.PLSSVD()]:
#        clf.set_params(scale=True)
#        clf.fit(X, Y)

def main():

    d = load_linnerud()
    X = np.asarray(d.data)
    Y = d.target
    num_comp = 3
    tol = 5e-12

    X, mX = pls.center(X, return_means = True)
    X, sX = pls.scale(X, return_stds = True)
    Y, mY = pls.center(Y, return_means = True)
    Y, sY = pls.scale(Y, return_stds = True)
    Xorig = X.copy()

    plsr = pls.PLSR(num_comp = num_comp, center = False, scale = False,
                    tolerance=tol, max_iter=1000, soft_threshold = 0)
    plsr.fit(X, Y)
    Yhat = np.dot(Xorig, plsr.B)

#    plsr.W, plsr.T, plsr.P = pls.direct(plsr.W, plsr.T, plsr.P)
#    plsr.C, plsr.U, plsr.Q = pls.direct(plsr.C, plsr.U, plsr.Q)
#    plsr.Ws = dot(plsr.W, np.linalg.inv(dot(plsr.P.T,plsr.W)))

    XY = dot(X.T, Y)
    W, S, C = np.linalg.svd(XY)
    w1 = W[:,[0]]
    t1 = dot(X, w1)
    p1 = dot(X.T, t1) / dot(t1.T, t1)
    c1 = dot(Y.T, t1) / dot(t1.T, t1)
    X -= dot(t1,p1.T)
#    w1, t1, p1 = pls.direct(w1, t1, p1)

    if num_comp >= 2:
        XY = dot(X.T, Y)
        W, S, C = np.linalg.svd(XY)
        w2 = W[:,[0]]
        t2 = dot(X, w2)
        p2 = dot(X.T, t2) / dot(t2.T, t2)
        c2 = dot(Y.T, t2) / dot(t2.T, t2)
        X -= dot(t2,p2.T)
#        w2, t2, p2 = pls.direct(w2, t2, p2)

    if num_comp >= 3:
        XY = dot(X.T, Y)
        W, S, C = np.linalg.svd(XY)
        w3 = W[:,[0]]
        t3 = dot(X, w3)
        p3 = dot(X.T, t3) / dot(t3.T, t3)
        c3 = dot(Y.T, t3) / dot(t3.T, t3)
        X -= dot(t3,p3.T)
#        w3, t3, p3 = pls.direct(w3, t3, p3)

    W = np.zeros((X.shape[1],num_comp))
    W[:,0] = w1.ravel()
    if num_comp >= 2:
        W[:,1] = w2.ravel()
    if num_comp >= 3:
        W[:,2] = w3.ravel()
    T = np.zeros((X.shape[0],num_comp))
    T[:,0] = t1.ravel()
    if num_comp >= 2:
        T[:,1] = t2.ravel()
    if num_comp >= 3:
        T[:,2] = t3.ravel()
    P = np.zeros((X.shape[1],num_comp))
    P[:,0] = p1.ravel()
    if num_comp >= 2:
        P[:,1] = p2.ravel()
    if num_comp >= 3:
        P[:,2] = p3.ravel()
    C = np.zeros((Y.shape[1],num_comp))
    C[:,0] = c1.ravel()
    if num_comp >= 2:
        C[:,1] = c2.ravel()
    if num_comp >= 3:
        C[:,2] = c3.ravel()

    Ws  = dot(W, np.linalg.inv(dot(P.T, W)))
    ws1 = Ws[:,[0]]
    ws2 = Ws[:,[1]]
    ws3 = Ws[:,[2]]
    B   = dot(Ws, C.T)
    Yhat_ = dot(Xorig, B)

    print Y
    print Yhat
    print Yhat_

    print "diff w1:", np.sum((w1 - plsr.W[:,[0]])**2)
    if num_comp >= 2:
        if dot(w2.T, plsr.W[:,[1]]) < 0:
            w2 *= -1
        print "diff w2:", np.sum((w2 - plsr.W[:,[1]])**2)
    if num_comp >= 3:
        if dot(w3.T, plsr.W[:,[2]]) < 0:
            w3 *= -1
        print "diff w3:", np.sum((w3 - plsr.W[:,[2]])**2)

    print "diff t1:", np.sum((t1 - plsr.T[:,[0]])**2)
    if num_comp >= 2:
        if dot(t2.T, plsr.T[:,[1]]) < 0:
            t2 *= -1
        print "diff t2:", np.sum((t2 - plsr.T[:,[1]])**2)
    if num_comp >= 3:
        if dot(t3.T, plsr.T[:,[2]]) < 0:
            t3 *= -1
        print "diff t3:", np.sum((t3 - plsr.T[:,[2]])**2)

    print "diff p1:", np.sum((p1 - plsr.P[:,[0]])**2)
    if num_comp >= 2:
        if dot(p2.T, plsr.P[:,[1]]) < 0:
            p2 *= -1
        print "diff p2:", np.sum((p2 - plsr.P[:,[1]])**2)
    if num_comp >= 3:
        if dot(p3.T, plsr.P[:,[2]]) < 0:
            p3 *= -1
        print "diff p3:", np.sum((p3 - plsr.P[:,[2]])**2)

    print "diff ws1:", np.sum((ws1 - plsr.Ws[:,[0]])**2)
    if num_comp >= 2:
        if dot(ws2.T, plsr.Ws[:,[1]]) < 0:
            ws2 *= -1
        print "diff ws2:", np.sum((ws2 - plsr.Ws[:,[1]])**2)
    if num_comp >= 3:
        if dot(ws3.T, plsr.Ws[:,[2]]) < 0:
            ws3 *= -1
        print "diff ws3:", np.sum((ws3 - plsr.Ws[:,[2]])**2)

    print "diff c1:", np.sum((c1 - plsr.C[:,[0]])**2), np.linalg.norm(c1), np.linalg.norm(plsr.C[:,[0]])
    if num_comp >= 2:
        if dot(c2.T, plsr.C[:,[1]]) < 0:
            c2 *= -1
        print "diff c2:", np.sum((c2 - plsr.C[:,[1]])**2)
    if num_comp >= 3:
        if dot(c3.T, plsr.C[:,[2]]) < 0:
            c3 *= -1
        print "diff c3:", np.sum((c3 - plsr.C[:,[2]])**2)

    print "diff b1:", np.sum((B[:,[0]] - plsr.B[:,[0]])**2)
    if num_comp >= 2:
        print "diff b1:", np.sum((B[:,[1]] - plsr.B[:,[1]])**2)
    if num_comp >= 3:
        print "diff b1:", np.sum((B[:,[2]] - plsr.B[:,[2]])**2)

    print B
    print plsr.B

    print "ss left:", np.sum(X**2)

    SSY     = np.sum(Y**2)
    SSYdiff = np.sum((Y-Yhat_)**2)
    print "Comparing original and predicted Y... OK!"\
            " (R2Yhat = %.4f)" % (1 - (SSYdiff / SSY))
    SSYdiff = np.sum((Y-Yhat)**2)
    print "Comparing original and predicted Y... OK!"\
            " (R2Yhat = %.4f)" % (1 - (SSYdiff / SSY))

    return

    # Compare SVD with and without sparsity constraint to numpy.linalg.svd

    Xtr = np.random.rand(6,6)
    Xte = np.random.rand(2,6)
    num_comp = 3

    Xtr, m = pls.center(Xtr, return_means = True)
    Xtr, s = pls.scale(Xtr, return_stds = True)
    Xte = (Xte - m) / s

    tol = 5e-12
    for st in [0.1, 0.01, 0.001, 0.0001, 0]:
        # NIPALS.SVD
        svd = pls.SVD(num_comp = num_comp,
                      tolerance=tol, max_iter=1000, soft_threshold = st)
        svd.fit(Xtr)
        svd.V, svd.U = pls.direct(svd.V, svd.U)

        # numpy.lialg.svd
        U, S, V = np.linalg.svd(Xtr)
        V = V.T
        S = np.diag(S)
        U = U[:,0:num_comp]
        S = S[:,0:num_comp]
        V = V[:,0:num_comp]
        V, U = pls.direct(V, U)
        SVDte = dot(Xte, V)

        if st < tol:
            num_decimals = 5
        else:
            num_decimals = int(log(1./st, 10) + 0.5)
        assert_array_almost_equal(svd.V, V, decimal=num_decimals-1,
                err_msg="NIPALS SVD and numpy.linalg.svd implementations " \
                "lead to different loadings")
        print "Comparing loadings of NIPALS SVD and numpy.linalg.svd... OK!"\
                " (max diff = %.4f, threshold = %0.4f)" % (np.max(np.abs(V - svd.V)), st)

    return


    # Compare PCA with sparsity constraint to numpy.linalg.svd

    Xtr = np.random.rand(5,5)
    Xte = np.random.rand(2,5)
    num_comp = 3

    Xtr, m = pls.center(Xtr, return_means = True)
    Xtr, s = pls.scale(Xtr, return_stds = True)
    Xte = (Xte - m) / s

    pca = pls.PCA(center = False, scale = False, num_comp = num_comp,
              tolerance=5e-12, max_iter=1000, soft_threshold = 0.1)
    pca.fit(Xtr)
    pca.P, pca.T = pls.direct(pca.P, pca.T)
    Tte = pca.transform(Xte)

    print pca.W

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

if __name__ == "__main__":
    main()
