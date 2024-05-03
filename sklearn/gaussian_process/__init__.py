# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Vincent Dubourg <vincent.dubourg@gmail.com>
#         Conrad Stevens <conrad.stevens@sydney.uni.com>
#         (mostly translation, see implementation details)
# License: BSD 3 clause

"""
The :mod:`sklearn.gaussian_process` module implements Gaussian Process
based regression and classification.
"""

from . import kernels
from ._gpc import GaussianProcessClassifier
from ._gpr import GaussianProcessRegressor
from ._tpr import TProcessRegressor

__all__ = ["GaussianProcessRegressor",
           "GaussianProcessClassifier",
           "TProcessRegressor",
           "kernels"]
