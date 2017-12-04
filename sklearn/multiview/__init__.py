"""
The :mod:`sklearn.multivew` module implements multiview data techniques.
"""

from .cpcmv import MVCPC
from .mvmds import MVMDS
from .mvsc import MVSC
from .mvtsne import MvtSNE

__all__ = ['MVCPC', 'MVMDS', 'MVSC', 'MvtSNE']
