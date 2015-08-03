# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# Licence: BSD 3 clause

"""
The :mod:`sklearn.gaussian_process` module implements scalar Gaussian Process
based predictions.
"""

from .gpr import GaussianProcessRegressor
from .gpc import (GaussianProcessClassifier,
    BinaryGaussianProcessClassifierLaplace)
from . import kernels

from .gaussian_process import GaussianProcess
from . import correlation_models
from . import regression_models

__all__ = ['GaussianProcess', 'correlation_models', 'regression_models']
