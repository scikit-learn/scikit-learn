# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# License: BSD 3 clause

"""
The :mod:`sklearn.gaussian_process` module implements Gaussian Process
based regression and classification.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    ["kernels"],
    submod_attrs={
        "_gpc": ["GaussianProcessClassifier"],
        "_gpr": ["GaussianProcessRegressor"],
    },
)
