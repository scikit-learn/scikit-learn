# Author: Elvis Dohmatob <gmdopp@gmail.com>

import numpy as np
from sklearn.utils import check_random_state
from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_array_almost_equal)
from .prox_slow import prox_l1_slow, prox_l2_slow, proj_l1_slow, proj_l2_slow
from sklearn.linear_model.python_wrappers import (_py_proj_l2, _py_prox_l2,
                                                  _py_prox_l1, _py_proj_l1)
from sklearn.linear_model.coordescendant import (
    L1_PENALTY, L11_PENALTY, L2INF_CONSTRAINT, L2_CONSTRAINT, L1_CONSTRAINT,
    L1INF_CONSTRAINT)
from .test_coordescendant import as_complex


def test_penalty_model_pseudos():
    assert_equal(L1_PENALTY, L11_PENALTY)
    assert_equal(L2INF_CONSTRAINT, L2_CONSTRAINT)
    assert_equal(L1INF_CONSTRAINT, L1_CONSTRAINT)


def test_prox_l2_basic():
    w = np.array([2., 0.])
    _py_prox_l2(w, 1., 1.)
    assert_array_equal(w, [1., 0.])
    w = np.array([3., 4.])
    _py_prox_l2(w, 5., 1.)
    assert_array_equal(w, 0.)

    w = np.array([3., 4.])
    old_w = w.copy()
    _py_prox_l2(w, 0., 1.)
    assert_array_equal(w, old_w)


def test_prox_l1_basic():
    w = np.array([2., 1., 0., .5, -4.])
    _py_prox_l1(w, 1., 1.)
    assert_array_equal(w, [1., 0., 0., 0., -3.])
    w = np.array([3., 4.])
    old_w = w.copy()
    _py_prox_l1(w, 0., 1.)
    assert_array_equal(w, old_w)


def test_prox_cython_equals_python(random_state=0):
    rng = check_random_state(random_state)
    n = 5
    for _ in range(6):
        for real in [True, False]:
            w = rng.randn(n)
            if not real:
                w = as_complex(w, rng.randn(n))

            for cython_prox, python_prox in zip(
                    [_py_prox_l1, _py_prox_l2, _py_proj_l1, _py_proj_l2],
                    [prox_l1_slow, prox_l2_slow, proj_l1_slow, proj_l2_slow]):
                if python_prox == proj_l1_slow and not real:
                    continue
                for reg in [0., .1, 1., 10.]:
                    for ajj in [0., 1e-4, 1., 1e4]:
                        w_ = w.copy()
                        w__ = w.copy()
                        cython_prox(w_, reg, ajj)
                        python_prox(w__, reg, ajj)
                        assert_array_almost_equal(w__, w_, decimal=10)


def test_prox_scaling_trick(random_state=0):
    rng = check_random_state(random_state)
    n = 5
    for _ in range(6):
        for real in [True, False]:
            w = rng.randn(n)
            if not real:
                w = as_complex(w, rng.randn(n))

            for prox in [_py_prox_l1, _py_prox_l2, prox_l1_slow, prox_l2_slow]:
                for reg in [0., .1, 1., 10.]:
                    for ajj in [1e-4, 1., 1e4]:
                        old_w = w.copy()
                        w_ = w / ajj
                        prox(w, reg, ajj)
                        prox(w_, reg / ajj, 1.)
                        assert_array_almost_equal(w_, w, decimal=10)
                        w = old_w.copy()


def test_proj_zero_radius(random_state=0):
    rng = check_random_state(random_state)
    n = 5
    for _ in range(6):
        for real in [True, False]:
            w = rng.randn(n)
            if not real:
                w = as_complex(w, rng.randn(n))
            for proj in [_py_proj_l1, proj_l1_slow, _py_proj_l2, proj_l2_slow][2:3]:
                if not real and proj in [_py_proj_l1, proj_l1_slow]:
                    continue
                for ajj in [0., .1, 1., 1e4]:
                    w_ = w.copy()
                    proj(w_, 0., ajj)
                    assert_array_equal(w_, 0.,
                                       err_msg="w=%s, ajj=%g" % (w, ajj))


def test_proj_zero():
    n = 5
    for real in [True, False]:
        w = np.zeros(n)
        if not real:
            w = as_complex(w, w)
    for proj in [_py_proj_l1, proj_l1_slow, _py_proj_l2, proj_l2_slow]:
        if not real and proj in [_py_proj_l1, proj_l1_slow]:
            continue
        for reg in [.1, 1., 10.]:
            for ajj in [.1, 1., 10.]:
                w_ = w.copy()
                proj(w_, reg, ajj)
                assert_array_equal(w_, 0.)


def test_prox_zero_penalty(random_state=0):
    rng = check_random_state(random_state)
    n = 5
    for _ in range(6):
        for real in [True, False]:
            w = rng.randn(n)
            if not real:
                w = as_complex(w, rng.randn(n))
            for prox in [_py_prox_l1, _py_prox_l2, prox_l1_slow, prox_l2_slow]:
                for ajj in [0., .1, 10.]:
                    w_ = w.copy()
                    prox(w_, 0., ajj)
                    if ajj == 0.:
                        assert_array_equal(w_, 0.)
                    else:
                        assert_array_almost_equal(w_, w / ajj, decimal=10)
