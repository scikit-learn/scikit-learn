"""
The :mod:`sklearn.kernel_ridge` module includes models based on
Kernel Ridge.
"""

# License: BSD 3 clause

from .krc import KernelRidgeClassification
from .krr import KernelRidgeRegression
from .kernel_ridge import KernelRidge

KRC = KernelRidgeClassification
KRR = KernelRidgeRegression

__all__ = ["KernelRidge"
           "KernelRidgeRegression",
           "KRR",
           "KernelRidgeClassification",
           "KRC"]
