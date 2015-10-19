"""Projection on the elastic-net ball
**References:**

J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)
"""
import numpy as np
from .enet_proj_fast import enet_projection_inplace
from .enet_proj_fast import enet_norm as c_enet_norm
from . import check_array


def enet_projection(v, radius=1., l1_ratio=0.1, check_input=True):
    if check_input:
        v = check_array(v, dtype=np.float64, order='C', copy=False)
    b = np.zeros(v.shape, dtype=np.float64)
    if v.ndim == 1:
        enet_projection_inplace(v, b, radius, l1_ratio)
    else:
        for i in range(v.shape[0]):
            enet_projection_inplace(v[i], b[i], radius, l1_ratio)
    return b


def enet_norm(v, l1_ratio=0.1):
    v = check_array(v, dtype=np.float64, order='C', copy=False)
    if v.ndim == 1:
        return c_enet_norm(v, l1_ratio)
    else:
        m = v.shape[0]
        norms = np.zeros(m, dtype=np.float64)
        for i in range(m):
            norms[i] = c_enet_norm(v[i], l1_ratio)
    return norms

def enet_scale(v, l1_ratio=0.1, radius=1, inplace=False):
    if not inplace:
        v = v.copy()
    l1_v = np.sum(np.abs(v), axis=1) * l1_ratio
    l1_v[l1_v == 0] = 1
    if l1_ratio != 1:
        l2_v = np.sum(v ** 2, axis=1) * (1 - l1_ratio)
        S = (- l1_v + np.sqrt(l1_v ** 2 + 4 * radius * l2_v))
        l2_v[l2_v == 0] = 1
        S /= (2 * l2_v)
        v *= S[:, np.newaxis]
    else:
        v /= l1_v[:, np.newaxis] / radius
    return v


def enet_threshold(v, l1_ratio=0.1, radius=1, inplace=False):
    if not inplace:
        v = v.copy()
    Sv = np.sqrt(np.sum(v ** 2, axis=1)) / radius
    Sv[Sv == 0] = 1
    v[:] = enet_projection(v / Sv[:, np.newaxis], l1_ratio=l1_ratio,
                        radius=radius)
    Sb = np.sqrt(np.sum(v ** 2, axis=1)) / radius
    Sb[Sb == 0] = 1
    v *= (Sv / Sb)[:, np.newaxis]
    return v


def l2_sphere_projection(v, radius=1, inplace=False):
    if not inplace:
        v = v.copy()
    Sv = np.sqrt(np.sum(v ** 2, axis=1)) / radius
    Sv[Sv == 0] = 1
    v /= Sv[:, np.newaxis]
    return v