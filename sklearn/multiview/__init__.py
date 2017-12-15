"""
The :mod:`sklearn.multivew` module implements multiview data techniques.
"""

from .cpcmv import MVCPC
from .mvsc import MVSC
from .mvmds import MVMDS
from .utils import x2p, whiten
from .mvtsne import MvtSNE

__all__ = ['MVCPC', 'MVMDS', 'MVSC', 'MvtSNE']
