"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

from .rbm import BernoulliRBM
from .random_basis_function import RandomBasisFunction

__all__ = ['BernoulliRBM', 'RandomBasisFunction']
