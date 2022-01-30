from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Callable,
)

import numpy as np

from pandas._typing import Scalar
from pandas.compat._optional import import_optional_dependency

from pandas.core.util.numba_ import (
    NUMBA_FUNC_CACHE,
    get_jit_arguments,
)


def generate_shared_aggregator(
    func: Callable[..., Scalar],
    engine_kwargs: dict[str, bool] | None,
    cache_key_str: str,
):
    """
    Generate a Numba function that loops over the columns 2D object and applies
    a 1D numba kernel over each column.

    Parameters
    ----------
    func : function
        aggregation function to be applied to each column
    engine_kwargs : dict
        dictionary of arguments to be passed into numba.jit
    cache_key_str: str
        string to access the compiled function of the form
        <caller_type>_<aggregation_type> e.g. rolling_mean, groupby_mean

    Returns
    -------
    Numba function
    """
    nopython, nogil, parallel = get_jit_arguments(engine_kwargs, None)

    cache_key = (func, cache_key_str)
    if cache_key in NUMBA_FUNC_CACHE:
        return NUMBA_FUNC_CACHE[cache_key]

    if TYPE_CHECKING:
        import numba
    else:
        numba = import_optional_dependency("numba")

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def column_looper(
        values: np.ndarray,
        start: np.ndarray,
        end: np.ndarray,
        min_periods: int,
        *args,
    ):
        result = np.empty((len(start), values.shape[1]), dtype=np.float64)
        for i in numba.prange(values.shape[1]):
            result[:, i] = func(values[:, i], start, end, min_periods, *args)
        return result

    return column_looper
