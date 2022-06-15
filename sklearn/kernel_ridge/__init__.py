"""
The :mod:`sklearn.kernel_model` module implements a variety of kernel models.
"""

# http://scikit-learn.sourceforge.net/modules/kernel_model.html for
# complete documentation.

from .krr import KernelRidge
from .krc import KernelRidgeClassifier

__all__ = [
    "KernelRidge",
    "KernelRidgeClassifier",
]
