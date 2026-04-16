import numpy as np
from numpy.typing import NDArray
from typing import TypeVar

# this mimicks the `ctypedef fused lapack_t`
_LapackT = TypeVar("_LapackT", np.float32, np.float64, np.complex64, np.complex128)

def lu_dispatcher(a: NDArray[_LapackT], u: NDArray[_LapackT], piv: NDArray[np.integer],
                  permute_l: bool) -> None: ...
