# -*- coding: utf-8 -*-

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#          Carlos Perales <sir.perales@gmail.com>
# License: BSD 3 clause

from .krc import KernelRidgeClassification
from .krr import KernelRidgeRegression

KernelRidge = KernelRidgeRegression
KRR = KernelRidgeRegression
KRC = KernelRidgeClassification

__all__ = ["KernelRidge",
           "KernelRidgeRegression",
           "KRR",
           "KernelRidgeClassification",
           "KRC"]
