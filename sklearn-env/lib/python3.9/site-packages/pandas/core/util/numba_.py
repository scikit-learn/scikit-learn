"""Common utilities for Numba operations"""
from __future__ import annotations

import types
from typing import (
    TYPE_CHECKING,
    Callable,
)

import numpy as np

from pandas.compat._optional import import_optional_dependency
from pandas.errors import NumbaUtilError

GLOBAL_USE_NUMBA: bool = False
NUMBA_FUNC_CACHE: dict[tuple[Callable, str], Callable] = {}


def maybe_use_numba(engine: str | None) -> bool:
    """Signal whether to use numba routines."""
    return engine == "numba" or (engine is None and GLOBAL_USE_NUMBA)


def set_use_numba(enable: bool = False) -> None:
    global GLOBAL_USE_NUMBA
    if enable:
        import_optional_dependency("numba")
    GLOBAL_USE_NUMBA = enable


def get_jit_arguments(
    engine_kwargs: dict[str, bool] | None = None, kwargs: dict | None = None
) -> tuple[bool, bool, bool]:
    """
    Return arguments to pass to numba.JIT, falling back on pandas default JIT settings.

    Parameters
    ----------
    engine_kwargs : dict, default None
        user passed keyword arguments for numba.JIT
    kwargs : dict, default None
        user passed keyword arguments to pass into the JITed function

    Returns
    -------
    (bool, bool, bool)
        nopython, nogil, parallel

    Raises
    ------
    NumbaUtilError
    """
    if engine_kwargs is None:
        engine_kwargs = {}

    nopython = engine_kwargs.get("nopython", True)
    if kwargs and nopython:
        raise NumbaUtilError(
            "numba does not support kwargs with nopython=True: "
            "https://github.com/numba/numba/issues/2916"
        )
    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)
    return nopython, nogil, parallel


def jit_user_function(
    func: Callable, nopython: bool, nogil: bool, parallel: bool
) -> Callable:
    """
    JIT the user's function given the configurable arguments.

    Parameters
    ----------
    func : function
        user defined function
    nopython : bool
        nopython parameter for numba.JIT
    nogil : bool
        nogil parameter for numba.JIT
    parallel : bool
        parallel parameter for numba.JIT

    Returns
    -------
    function
        Numba JITed function
    """
    if TYPE_CHECKING:
        import numba
    else:
        numba = import_optional_dependency("numba")

    if numba.extending.is_jitted(func):
        # Don't jit a user passed jitted function
        numba_func = func
    else:

        @numba.generated_jit(nopython=nopython, nogil=nogil, parallel=parallel)
        def numba_func(data, *_args):
            if getattr(np, func.__name__, False) is func or isinstance(
                func, types.BuiltinFunctionType
            ):
                jf = func
            else:
                jf = numba.jit(func, nopython=nopython, nogil=nogil)

            def impl(data, *_args):
                return jf(data, *_args)

            return impl

    return numba_func
