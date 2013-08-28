"""Bagging meta-estimator"""

#Author: Maheshakya Wijewardena <pmaheshakya4@gmail.com>
#License: BSD 3 clause

from __future__ import division

import itertools
import numbers
import numpy as np
from warnings import warn
from abc import ABCMeta, abstractmethod
from inspect import getargspec

"""
if __name__ == "__main__" and __package__ is None:
    __package__ = "sklearn.ensemble"
"""

from .base import BaseEnsemble

class BaseBagging():

