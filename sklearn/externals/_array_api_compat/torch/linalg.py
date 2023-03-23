from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import torch
    array = torch.Tensor

from torch.linalg import *

# torch.linalg doesn't define __all__
# from torch.linalg import __all__ as linalg_all
from torch import linalg as torch_linalg
linalg_all = [i for i in dir(torch_linalg) if not i.startswith('_')]

# These are implemented in torch but aren't in the linalg namespace
from torch import outer, trace
from ._aliases import _fix_promotion, matrix_transpose, tensordot

# Note: torch.linalg.cross does not default to axis=-1 (it defaults to the
# first axis with size 3), see https://github.com/pytorch/pytorch/issues/58743
def cross(x1: array, x2: array, /, *, axis: int = -1) -> array:
    x1, x2 = _fix_promotion(x1, x2, only_scalar=False)
    return torch_linalg.cross(x1, x2, dim=axis)

__all__ = linalg_all + ['outer', 'trace', 'matrix_transpose', 'tensordot']

del linalg_all
