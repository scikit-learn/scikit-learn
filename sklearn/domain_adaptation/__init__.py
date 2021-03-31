"""
The :mod:`sklearn.domain_adaptation` module implements domain adaptation
algorithms. These algorithms minimize the domain gaps between domains
to improve the classification performance. This module includes Transfer 
Component Analysis (TCA) and Balanced Distribution Adaptation (BDA).
"""

from .tca import TCA
from .bda import BDA

__all__ = ['TCA', 'BDA']