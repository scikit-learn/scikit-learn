#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    This module contains a contribution to the scikit-learn project that
    implements Gaussian Process based prediction (also known as Kriging).
    
    The present implementation is based on a transliteration of the DACE
    Matlab toolbox <http://www2.imm.dtu.dk/~hbn/dace/>.
"""

from .gaussian_process_model import GaussianProcessModel
from .correlation_models import *
from .regression_models import *