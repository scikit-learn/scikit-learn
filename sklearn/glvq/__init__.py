# -*- coding: utf-8 -*-

# Author: Joris Jensen <jjensen@techfak.uni-bielefeld.de>
#
# License: BSD 3 clause

"""
The :mod:`sklearn.glvq` module implements Generalized Learning Vector
Quantization based classification.
"""

from .glvq import GlvqModel
from .grlvq import GrlvqModel
from .gmlvq import GmlvqModel
from .lgmlvq import LgmlvqModel

__all__ = ['GlvqModel', 'GrlvqModel', 'GmlvqModel', 'LgmlvqModel']
