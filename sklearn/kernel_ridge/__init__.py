"""
The :mod:`sklearn.kernel_model` module implements a variety of kernel models.
"""

# http://scikit-learn.sourceforge.net/modules/kernel_model.html for
# complete documentation.

from .kernel_ridge import KernelRidge

__all__ = [
    "KernelRidge",
    "KernelRidgeClassification",
    ]
