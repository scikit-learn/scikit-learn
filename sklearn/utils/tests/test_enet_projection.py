# Author: Arthur Mensch
# License: BSD 3 clause

import numpy as np
from numpy import sqrt
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.enet_projection import enet_norm, enet_projection


def _enet_norm_for_projection(v, gamma):
    return np.sum(v * (1 + gamma / 2 * v))


def enet_norm_slow(v, l1_ratio=0.1):
    if l1_ratio == 0:
        return sqrt(np.sum(v ** 2))
    b_abs = np.abs(v)
    return np.sum(b_abs * (l1_ratio + b_abs * (1 - l1_ratio)))


def enet_projection_slow(v, radius=1, l1_ratio=0.1):
    """Projection on the elastic-net ball
    **References:**

    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)
    """
    random_state = check_random_state(None)
    if l1_ratio == 0:
        return v / sqrt(np.sum(v ** 2))
    gamma = 2 / l1_ratio - 2
    radius /= l1_ratio
    m = v.shape[0]
    b_abs = np.abs(v)
    norm = _enet_norm_for_projection(b_abs, gamma)
    if norm <= radius:
        return v
    else:
        s = 0
        rho = 0
        U = np.arange(m)
        mask = np.ones(m, dtype=np.bool)
        mask_non_zero = mask.nonzero()[0]
        while mask_non_zero.shape[0] != 0:
            k = random_state.randint(mask_non_zero.shape[0])
            idx = mask_non_zero[k]
            k = U[idx]
            sel = b_abs < b_abs[k]
            G = U[~sel * mask]
            d_rho = G.shape[0]
            d_s = _enet_norm_for_projection(b_abs[G], gamma)
            if s + d_s - (rho + d_rho) * (1 + gamma / 2 * b_abs[k]) * b_abs[k]\
                    < radius * (1 + gamma * b_abs[k]) ** 2:
                s += d_s
                rho += d_rho
                mask *= sel
            else:
                mask *= ~sel
                mask[idx] = False
            mask_non_zero = mask.nonzero()[0]
        if gamma != 0:
            a = gamma ** 2 * radius + gamma * rho * 0.5
            b_ = 2 * radius * gamma + rho
            c = radius - s
            l = (-b_ + np.sqrt(b_ ** 2 - 4 * a * c)) / (2*a)
        else:
            l = (s - radius) / rho
        b_sign = np.sign(v)
        b_sign[b_sign == 0] = 1
        return b_sign * np.maximum(np.zeros_like(b_abs), b_abs - l)\
               / (1 + l * gamma)


def test_slow_enet_norm():
    norms = np.zeros(10)
    norms2 = np.zeros(10)
    random_state = check_random_state(0)

    for i in range(10):
        a = random_state.randn(10000)
        norms[i] = enet_norm_slow(a, l1_ratio=0.1)
        norms2[i] = (1 - 0.1) * (a ** 2).sum() + 0.1 * np.abs(a).sum()
    assert_array_almost_equal(norms, norms2)

def test_slow_enet_projection_norm():
    norms = np.zeros(10)
    random_state = check_random_state(0)

    for i in range(10):
        a = random_state.randn(10000)
        b = np.asarray(enet_projection_slow(a, radius=1, l1_ratio=0.1))
        norms[i] = enet_norm_slow(b, l1_ratio=0.1)
    assert_array_almost_equal(norms, np.ones(10))


def test_fast_enet_projection_norm():
    random_state = check_random_state(0)
    norms = np.zeros(10)
    for i in range(10):
        a = random_state.randn(20000)
        a /= np.sqrt(np.sum(a ** 2))
        c = np.zeros(20000)
        c[:] = enet_projection(a, 1, 0.15)
        norms[i] = enet_norm(c, l1_ratio=0.15)
    assert_array_almost_equal(norms, np.ones(10))


def test_fast_enet_projection():
    c = np.empty((10, 100))
    b = np.empty((10, 100))
    random_state = check_random_state(0)
    for i in range(10):
        a = random_state.randn(100)
        b[i, :] = enet_projection_slow(a, radius=1, l1_ratio=0.1)
        c[i] = enet_projection(a, radius=1, l1_ratio=0.1)
    assert_array_almost_equal(c, b, 4)


def test_fast_enet_l2_ball():
    random_state = check_random_state(0)
    norms = np.zeros(10)
    for i in range(10):
        a = random_state.randn(100)
        c = np.zeros(100)
        c[:] = enet_projection(a, 1, 0.0)
        norms[i] = np.sqrt(np.sum(c ** 2))
    assert_array_almost_equal(norms, np.ones(10))


def test_fast_enet_l1_ball():
    random_state = check_random_state(0)
    norms = np.zeros(10)
    for i in range(10):
        a = random_state.randn(100)
        b = np.zeros(100)
        b[:] = enet_projection(a, 1, 1.0)
        norms[i] = np.sum(np.abs(b))
    assert_array_almost_equal(norms, np.ones(10))
