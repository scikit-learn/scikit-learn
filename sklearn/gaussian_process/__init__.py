#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# License: BSD style

"""
A module that implements scalar Gaussian Process based prediction (also
known as Kriging).

Contains
--------
GaussianProcess: The main class of the module that implements the Gaussian
                 Process prediction theory.

regression_models: A submodule that contains the built-in regression models.

correlation_models: A submodule that contains the built-in correlation models.
"""

from .gaussian_process import GaussianProcess
from . import correlation_models
from . import regression_models
