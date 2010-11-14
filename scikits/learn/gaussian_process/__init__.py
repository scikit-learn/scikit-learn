#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    This module contains a contribution to the scikit-learn project that
    implements Gaussian Process based prediction (also known as Kriging).

    The present implementation is based on a transliteration of the DACE
    Matlab toolbox <http://www2.imm.dtu.dk/~hbn/dace/>.
"""

from .gaussian_process import GaussianProcess
from .correlation import corrlin, corrcubic, correxp1, \
                         correxp2, correxpg, corriid
from .regression import regpoly0, regpoly1, regpoly2
