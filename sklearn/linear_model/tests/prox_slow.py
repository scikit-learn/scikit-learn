# Synopsis: Pure-python (slow!) implementation of proximal operators.
#           This is only for testing purposes and is not meant to be
#           used otherwise. See cython implementations.
#
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD
#
# What's a proximal operator ?
#
#          := argmin .5 * ||z - w / ajj||_2^2 + (reg / ajj) * pen(z)
#               z
#
# where ajj > 0 is a constant. ajj = 0 corresponds to argmin_z pen(z)
# while ajj = 1 corresponds to classical definition of the prox.

from math import sqrt
import numpy as np


def nrm2(w, squared=False):
    """L2 norm of real or complex vector."""
    ww = np.dot(w.conjugate(), w)
    if w[0].dtype in [np.complex64, np.complex128]:
        ww = ww.real
    return ww if squared else sqrt(ww)


def prox_l1_slow(w, reg, ajj):
    """Soft-thresholding operator"""
    if ajj == 0.:
        w[:] = 0.
    else:
        w_nz = w.nonzero()
        shrink = np.maximum(1. - reg / np.abs(w[w_nz]), 0.)
        if ajj != 1.:
            shrink /= ajj
        w[w_nz] *= shrink
    return w


def proj_l1_slow(w, reg, ajj, l1_ratio=1.):
    """Projection onto L1 ball"""
    try:
        from modl.utils.math.enet import enet_projection
    except ImportError:
        raise NotImplementedError("proj_l1")
    if ajj == 0. or reg == 0.:
        w[:] = 0.
    else:
        if ajj != 1.:
            w /= ajj
        out = np.zeros_like(w)
        enet_projection(w, out, reg, l1_ratio)
        w[:] = out


def proj_l2_slow(w, reg, ajj):
    """Projection onto L1 ball"""
    if ajj == 0.:
        w[:] = 0.
    else:
        norm = nrm2(w)
        if norm > ajj * reg:
            scaling = reg / norm
        else:
            scaling = 1. / ajj
        if scaling != 1.:
            w *= scaling
    return w


def prox_l2_slow(w, reg, ajj):
    """Proximal operator of L2 norm"""
    if ajj == 0.:
        w[:] = 0.
    else:
        norm = nrm2(w)
        if norm > 0.:
            shrink = np.maximum(1 - reg / norm, 0) / ajj
            if shrink != 1.:
                w *= shrink
    return w
