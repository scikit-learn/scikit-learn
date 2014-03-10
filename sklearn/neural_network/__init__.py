"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

from .rbm import BernoulliRBM
from .sgvb import SGVB

__all__ = ['BernoulliRBM', 'SGVB']
