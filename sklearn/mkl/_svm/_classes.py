# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ...svm import _classes
from ._base import BaseLibSVMforMKL


class SVC(BaseLibSVMforMKL, _classes.SVC):
    pass


class SVR(BaseLibSVMforMKL, _classes.SVR):
    pass


class OneClassSVM(BaseLibSVMforMKL, _classes.OneClassSVM):
    pass


class NuSVC(BaseLibSVMforMKL, _classes.NuSVC):
    pass


class NuSVR(BaseLibSVMforMKL, _classes.NuSVR):
    pass
