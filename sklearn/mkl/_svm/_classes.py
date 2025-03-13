# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ...svm import _classes
from ._base import BaseLibSVMforMKL


class SVC(BaseLibSVMforMKL, _classes.SVC):
    def __init__(self, shrinking=False, **kwargs):
        super().__init__(shrinking=shrinking, **kwargs)


class SVR(BaseLibSVMforMKL, _classes.SVR):
    def __init__(self, shrinking=False, **kwargs):
        super().__init__(shrinking=shrinking, **kwargs)


class OneClassSVM(BaseLibSVMforMKL, _classes.OneClassSVM):
    def __init__(self, shrinking=False, **kwargs):
        super().__init__(shrinking=shrinking, **kwargs)


class NuSVC(BaseLibSVMforMKL, _classes.NuSVC):
    def __init__(self, shrinking=False, **kwargs):
        super().__init__(shrinking=shrinking, **kwargs)


class NuSVR(BaseLibSVMforMKL, _classes.NuSVR):
    def __init__(self, shrinking=False, **kwargs):
        super().__init__(shrinking=shrinking, **kwargs)
