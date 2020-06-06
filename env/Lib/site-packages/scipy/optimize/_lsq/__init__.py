"""This module contains least-squares algorithms."""
from __future__ import division, print_function, absolute_import

from .least_squares import least_squares
from .lsq_linear import lsq_linear

__all__ = ['least_squares', 'lsq_linear']
