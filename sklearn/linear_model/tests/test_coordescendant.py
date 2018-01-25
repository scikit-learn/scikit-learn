# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

import itertools
import numpy as np
import scipy
from sklearn.utils import check_random_state
from sklearn.utils.testing import (assert_array_almost_equal, assert_less,
                                   assert_equal, assert_less_equal,
                                   assert_array_equal, assert_false,
                                   assert_true, assert_array_less,
                                   assert_almost_equal, ignore_warnings)
from sklearn.linear_model import ElasticNet, MultiTaskElasticNet
from sklearn.linear_model.coordescendant import (coordescendant, L11_PENALTY,
                                                 L2INF_CONSTRAINT, L21_PENALTY,
                                                 L1INF_CONSTRAINT, NOP)
from .coordescendant_slow import coordescendant_slow
from .dual_gap_slow import compute_dual_gap_slow
from .prox_slow import prox_l1_slow, prox_l2_slow, proj_l1_slow, proj_l2_slow
from sklearn.linear_model.python_wrappers import (
    _py_prox_l1, _py_prox_l2, _py_proj_l1, _py_proj_l2, _py_compute_dual_gap)


def as_complex(X, Y):
    imag_unit = scipy.sqrt(-1)
    return X + imag_unit * Y


def _make_test_data(n_samples, n_features, n_targets, single_precision=False,
                    real=True, random_state=0, fake_complex=False,
                    sparsity=0., noise_std=0., double_precision=True):
    rng = check_random_state(random_state)
    if single_precision:
        real_dtype = np.float32
        complex_dtype = np.complex64
    else:
        real_dtype = np.float64
        complex_dtype = np.complex128
    if real and not fake_complex:
        complex_dtype = real_dtype
    noise_std = real_dtype(noise_std)
    X = rng.randn(n_samples, n_features)
    X = real_dtype(X)
    if X.ndim == 0:
        X = np.array([[X]])
    W = rng.randn(n_features, n_targets)
    W[rng.rand(*W.shape) < sparsity] = 0.
    W = real_dtype(W)
    if W.ndim == 0:
        W = np.array([[W]])
    if not real:
        X = as_complex(X, real_dtype(rng.randn(n_samples, n_features)))
        W = as_complex(W, real_dtype(rng.randn(n_features, n_targets)))
    elif fake_complex:
        X = complex_dtype(X)
        W = complex_dtype(W)
    Y = X.dot(W)
    if noise_std:
        if real and not fake_complex:
            Y += real_dtype(rng.randn(n_samples, n_targets))
        else:
            Y += noise_std * as_complex(
                real_dtype(rng.randn(n_samples, n_targets)),
                real_dtype(rng.randn(n_samples, n_targets)))
    Gram = np.dot(X.conjugate().T, X)
    Cov = np.dot(X.conjugate().T, Y)
    assert_array_almost_equal(
        Gram, Gram.conjugate().T)  # Gram must be conj symm
    X = np.asarray(X, order="F")
    Y = np.asarray(Y, order="F")
    Gram = np.asarray(Gram, order="F")
    Cov = np.asarray(Cov, order="F")
    W = np.asarray(W, order="C")
    assert_array_less(0., np.diag(Gram))
    assert_array_equal(X.shape, (n_samples, n_features))
    assert_array_equal(Y.shape, (n_samples, n_targets))
    assert_array_equal(Cov.shape, (n_features, n_targets))

    return X, Y, Gram, Cov, W


@ignore_warnings
def test_lstsq(random_state=0):
    for precompute, real in itertools.product([True, False],
                                              [True, False]):
        imag_unit = scipy.sqrt(-1)
        rng = check_random_state(random_state)
        X = np.diag([1., 10. + 2. * imag_unit])
        X[0, 1] = -4.
        X = np.vstack((X, np.ones(2)))
        W = rng.randn(2, 3) + imag_unit * rng.randn(2, 3)
        X = X.astype(np.complex128)
        W = W.astype(np.complex128)
        if real:
            W = W.real
            X = X.real
        Y = np.dot(X, W)
        Y_flat = Y.ravel()
        Y_norm2 = np.dot(Y_flat.conjugate(), Y_flat)
        if precompute:
            X_or_Gram = np.dot(X.conjugate().T, X)
            Y_or_Cov = np.dot(X.conjugate().T, Y)
        else:
            X_or_Gram = X
            Y_or_Cov = Y
        X_or_Gram, Y_or_Cov = map(np.asfortranarray, (X_or_Gram, Y_or_Cov))
        for backend in [coordescendant_slow, coordescendant]:
            W_ = np.ascontiguousarray(np.zeros_like(W, dtype=W.dtype))
            W_ = backend(
                W_, 0., 0., X_or_Gram, Y_or_Cov, precomputed=precompute,
                max_iter=10, penalty_model=NOP, Y_norm2=Y_norm2)[0]
            assert_array_almost_equal(W, W_, decimal=12)
            assert_array_almost_equal(Y, np.dot(X, W_), decimal=12)


@ignore_warnings
def test_cmp_with_sklearn(max_iter=10, X=None, Y=None, random_state=0):
    rng = check_random_state(random_state)
    for (n_samples, n_features, n_targets, penalty_model, alpha, l1_ratio,
         precompute, positive) in itertools.product(
             range(1, 3), range(1, 3), range(1, 3), [L11_PENALTY, L21_PENALTY],
             [.1, 1., 10.], [.5, 1.], [True, False], [True, False]):
        X, Y, Gram, Cov, _ = _make_test_data(
            n_samples, n_features, n_targets, random_state=rng)
        if positive and penalty_model != L11_PENALTY:
            continue
        kwargs = {}
        if penalty_model == L11_PENALTY:
            kwargs = dict(precompute=precompute, positive=positive)
            sk_cls = ElasticNet
        elif penalty_model == L21_PENALTY:
            if precompute:
                continue
            sk_cls = MultiTaskElasticNet

        alpha_ = alpha * n_samples
        l1_reg = alpha_ * l1_ratio
        l2_reg = alpha_ - l1_reg
        if n_targets > 1 and penalty_model == L11_PENALTY:
            tol = 0.
        else:
            tol = 1e-2
        assert_array_almost_equal(Gram, Gram.conjugate().T, decimal=14)
        sk_model = sk_cls(
            alpha=alpha, l1_ratio=l1_ratio, fit_intercept=False, tol=tol,
            random_state=rng, max_iter=max_iter, **kwargs)
        W1 = sk_model.fit(X.real, Y.real).coef_.T
        Y_flat = Y.ravel()
        Y_norm2 = np.dot(Y_flat.conjugate(), Y_flat)
        if precompute:
            X_or_Gram = Gram
            Y_or_Cov = Cov
        else:
            X_or_Gram = X
            Y_or_Cov = Y
        for backend in [coordescendant_slow, coordescendant]:
            W2 = np.zeros_like(Cov, order="C")
            W2, gap, _, n_iter = backend(
                W2, l1_reg, l2_reg, X_or_Gram, Y_or_Cov, Y_norm2=Y_norm2,
                precomputed=precompute, max_iter=max_iter, tol=tol,
                penalty_model=penalty_model, positive=positive)
            if W1.ndim == 1:
                W2 = W2.ravel()
            if n_targets == 1:
                assert_equal(n_iter, sk_model.n_iter_)
            assert_array_almost_equal(W1, W2.real, decimal=13)


@ignore_warnings
def test_user_prox(random_state=0, max_iter=10):
    rng = check_random_state(random_state)
    n_samples, n_features, n_targets = 4, 2, 3
    for ((penalty_model, proxes), real, fake_complex,
         alpha, l1_ratio, emulate_sklearn_dl) in itertools.product(
            zip([L11_PENALTY, L21_PENALTY, L2INF_CONSTRAINT, L1INF_CONSTRAINT],
                [[prox_l1_slow, _py_prox_l1], [prox_l2_slow, _py_prox_l2],
                 [proj_l2_slow, _py_proj_l2], [proj_l1_slow, _py_proj_l1]]),
             [True, False], [True, False], [0., 1., 10.], [0., .5, 1.],
             [True, False]):
        if penalty_model == L1INF_CONSTRAINT and (not real
                                                  or fake_complex):
            continue
        if emulate_sklearn_dl and penalty_model not in [L1INF_CONSTRAINT,
                                                        L2INF_CONSTRAINT]:
            continue
        for prox in proxes:
            reg = alpha * l1_ratio
            l2_reg = alpha - reg
            X, Y, Gram, Cov, _ = _make_test_data(
                n_samples, n_features, n_targets, random_state=rng, real=real,
                fake_complex=fake_complex)
            stuff = dict(W=[], tol=[], gap=[], n_iter=[])
            for precompute in [True, False][1:]:
                if precompute:
                    X_or_Gram = Gram
                    Y_or_Cov = Cov
                    Y_flat = Y.ravel()
                    Y_norm2 = np.dot(Y_flat.conjugate(), Y_flat)
                else:
                    X_or_Gram = X
                    Y_or_Cov = Y
                    Y_norm2 = np.nan
                X_or_Gram = np.asfortranarray(X_or_Gram)
                Y_or_Cov = np.asfortranarray(Y_or_Cov)
                for backend in [coordescendant_slow, coordescendant]:
                    W = np.zeros_like(Cov, order="C")
                    W, gap, tol, n_iter = backend(
                        W, reg, l2_reg, X_or_Gram, Y_or_Cov, max_iter=max_iter,
                        penalty_model=penalty_model, precomputed=precompute,
                        emulate_sklearn_dl=emulate_sklearn_dl, tol=0,
                        Y_norm2=Y_norm2)
                    assert_true(isinstance(gap, float), msg=backend)
                    stuff["W"].append(W)
                    stuff["gap"].append(gap)
                    stuff["tol"].append(tol)
                    stuff["n_iter"].append(n_iter)
                for name, vals in stuff.items():
                    for a, b in itertools.combinations(vals, 2):
                        if name == "W":
                            assert_false(np.any(np.isnan(a)))
                            assert_false(np.any(np.isnan(b)))
                            assert_array_almost_equal(a, b, decimal=13)
                        elif name in ["tol", "n_iter"]:
                            assert_equal(a, b)
                        else:
                            assert_almost_equal(a, b, decimal=13)


@ignore_warnings
def test_coordescendant_cython_equals_python(random_state=0, max_iter=10):
    rng = check_random_state(random_state)
    n_samples, n_features, n_targets = 4, 2, 3
    for (penalty_model, real, fake_complex, alpha, l1_ratio, positive,
         emulate_sklearn_dl) in itertools.product(
            [NOP, L11_PENALTY, L21_PENALTY, L2INF_CONSTRAINT,
             L1INF_CONSTRAINT], [True, False], [True, False],
             [0., 1., 10.], [0., .5, 1.], [True, False][:1], [True, False]):
        if penalty_model == L1INF_CONSTRAINT and (not real or
                                                  fake_complex):
            continue
        if emulate_sklearn_dl and penalty_model not \
           in [L1INF_CONSTRAINT, L2INF_CONSTRAINT]:
            continue
        if positive:
            if penalty_model != L11_PENALTY:
                continue
            if (fake_complex or not real):
                continue
            if emulate_sklearn_dl:
                continue
        reg = alpha * l1_ratio
        l2_reg = alpha - reg
        X, Y, Gram, Cov, _ = _make_test_data(
            n_samples, n_features, n_targets, random_state=rng,
            real=real, fake_complex=fake_complex)
        for precompute in [True, False]:
            if precompute:
                X_or_Gram = Gram
                Y_or_Cov = Cov
                Y_flat = Y.ravel()
                Y_norm2 = np.dot(Y_flat.conjugate(), Y_flat)
            else:
                X_or_Gram = X
                Y_or_Cov = Y
                Y_norm2 = np.nan
            W1 = np.zeros_like(Cov, order="C")
            W1, gap1, tol1, n_iter1 = coordescendant(
                W1, reg, l2_reg, X_or_Gram, Y_or_Cov,
                penalty_model=penalty_model, max_iter=max_iter,
                emulate_sklearn_dl=emulate_sklearn_dl,
                precomputed=precompute, Y_norm2=Y_norm2, positive=positive)
            assert_false(np.any(np.isnan(W1)))
            W2 = np.zeros_like(Cov, order="C")
            W2, gap2, tol2, n_iter2 = coordescendant_slow(
                W2, reg, l2_reg, X_or_Gram, Y_or_Cov,
                penalty_model=penalty_model, max_iter=max_iter,
                emulate_sklearn_dl=emulate_sklearn_dl,
                precomputed=precompute, Y_norm2=Y_norm2, positive=positive)
            assert_false(np.any(np.isnan(W2)))
            assert_equal(n_iter1, n_iter2)
            assert_equal(tol1, tol2)
            assert_almost_equal(gap1, gap2, decimal=1)
            assert_array_almost_equal(W1, W2, decimal=13)


@ignore_warnings
def test_precompute(random_state=0, max_iter=1):
    for (tol, sp, real, positive, n_samples, n_features, n_targets, alpha,
         l1_ratio, penalty_model) in itertools.product(
             [0, 1e-2], [True, False][1:], [True, False], [True, False][1:],
             range(1, 3), range(1, 3), range(1, 3), [0., 1., 10.],
             [1e-4, .5, 1.], [L11_PENALTY, L21_PENALTY]):
        if positive:
            if penalty_model != L11_PENALTY or not real:
                continue
        X, Y, Gram, Cov, _ = _make_test_data(
            n_samples, n_features, n_targets, real=real,
            random_state=random_state)
        reg = alpha * l1_ratio
        l2_reg = alpha - reg
        if penalty_model == L1INF_CONSTRAINT and not real:
            continue
        stuff = dict(gap=[], n_iter=[], W=[])
        Y_flat = Y.ravel()
        Y_norm2 = np.dot(Y_flat.conjugate(), Y_flat).real
        for precompute in [True, False]:
            if precompute:
                X_or_Gram = Gram
                Y_or_Cov = Cov
            else:
                X_or_Gram = X
                Y_or_Cov = Y
            for backend in [coordescendant_slow, coordescendant]:
                W = np.zeros_like(Cov, order="C")
                W, gap, tol, n_iter = backend(
                    W, reg, l2_reg, X_or_Gram, Y_or_Cov, max_iter=max_iter,
                    precomputed=precompute, penalty_model=penalty_model,
                    tol=tol, Y_norm2=Y_norm2, positive=positive)
                stuff["W"].append(W)
                stuff["gap"].append(gap)
                stuff["n_iter"].append(n_iter)

        for name, vals in stuff.items():
            for a, b in itertools.combinations(vals, 2):
                if name == "W":
                    assert_array_almost_equal(a, b, decimal=13)
                if name == "n_iter":
                    assert_equal(a, b)
                elif name == "gap":
                    assert_almost_equal(a, b, decimal=12)


@ignore_warnings
def test_compute_dual_gap(random_state=0):
    for (sp, real, positive, n_samples, n_features, n_targets, alpha, l1_ratio,
         penalty_model) in itertools.product(
             [True, False][1:], [True, False], [True, False], range(1, 3),
             range(1, 3), range(1, 3), [0., 1., 10.], [0., .5, 1.],
             [L11_PENALTY, L21_PENALTY]):
        if positive:
            if not real or penalty_model != L11_PENALTY:
                continue
        reg = l1_ratio * alpha
        l2_reg = alpha - reg
        if sp:
            reg = np.float32(reg)
            l2_reg = np.float32(reg)
        if sp:
            decimal = 4
        else:
            decimal = 12
        X, Y, Gram, Cov, W = _make_test_data(
            n_samples, n_features, n_targets, real=real, single_precision=sp,
            random_state=random_state)
        Grad = np.zeros_like(W, order="F")
        gaps = []
        for precompute in [True, False]:
            R = Y - np.dot(X, W)
            Y_norm2 = np.nan
            if precompute:
                X_or_Gram = Gram
                Y_or_Cov = Cov
                Y_flat = Y.ravel()
                Y_norm2 = np.dot(Y_flat.conjugate(), Y_flat).real
                R = np.dot(X.conjugate().T, R)
            else:
                X_or_Gram = X
                Y_or_Cov = Y
            R = np.asfortranarray(R)
            for dgap in [compute_dual_gap_slow, _py_compute_dual_gap][1:]:
                print(dgap)
                gap = dgap(
                    W, reg, l2_reg, X_or_Gram.conjugate(), Y_or_Cov, R, Grad,
                    Y_norm2=Y_norm2, precomputed=precompute,
                    penalty_model=penalty_model, positive=positive)
                assert_true(np.isreal(gap), msg="gap=%s" % gap)
                gaps.append(gap)
                assert_less_equal(0., gap, msg="%g (%s)" % (gap, dgap))
        for a, b in itertools.combinations(gaps, 2):
            assert_almost_equal(a, b, decimal=decimal)
