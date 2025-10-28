from __future__ import annotations

import re
import sys
from collections.abc import Hashable
from functools import cache
from importlib import import_module
from importlib.util import find_spec
from types import ModuleType
from typing import TYPE_CHECKING, Any, ClassVar, cast

_ALTAIR_AVAILABLE = True
_DELTALAKE_AVAILABLE = True
_FSSPEC_AVAILABLE = True
_GEVENT_AVAILABLE = True
_GREAT_TABLES_AVAILABLE = True
_HYPOTHESIS_AVAILABLE = True
_NUMPY_AVAILABLE = True
_PANDAS_AVAILABLE = True
_POLARS_CLOUD_AVAILABLE = True
_PYARROW_AVAILABLE = True
_PYDANTIC_AVAILABLE = True
_PYICEBERG_AVAILABLE = True
_TORCH_AVAILABLE = True
_PYTZ_AVAILABLE = True


class _LazyModule(ModuleType):
    """
    Module that can act both as a lazy-loader and as a proxy.

    Notes
    -----
    We do NOT register this module with `sys.modules` so as not to cause
    confusion in the global environment. This way we have a valid proxy
    module for our own use, but it lives *exclusively* within polars.
    """

    __lazy__ = True

    _mod_pfx: ClassVar[dict[str, str]] = {
        "numpy": "np.",
        "pandas": "pd.",
        "pyarrow": "pa.",
        "polars_cloud": "pc.",
    }

    def __init__(
        self,
        module_name: str,
        *,
        module_available: bool,
    ) -> None:
        """
        Initialise lazy-loading proxy module.

        Parameters
        ----------
        module_name : str
            the name of the module to lazy-load (if available).

        module_available : bool
            indicate if the referenced module is actually available (we will proxy it
            in both cases, but raise a helpful error when invoked if it doesn't exist).
        """
        self._module_available = module_available
        self._module_name = module_name
        self._globals = globals()
        super().__init__(module_name)

    def _import(self) -> ModuleType:
        # import the referenced module, replacing the proxy in this module's globals
        module = import_module(self.__name__)
        self._globals[self._module_name] = module
        self.__dict__.update(module.__dict__)
        return module

    def __getattr__(self, name: str) -> Any:
        # have "hasattr('__wrapped__')" return False without triggering import
        # (it's for decorators, not modules, but keeps "make doctest" happy)
        if name == "__wrapped__":
            msg = f"{self._module_name!r} object has no attribute {name!r}"
            raise AttributeError(msg)

        # accessing the proxy module's attributes triggers import of the real thing
        if self._module_available:
            # import the module and return the requested attribute
            module = self._import()
            return getattr(module, name)

        # user has not installed the proxied/lazy module
        elif name == "__name__":
            return self._module_name
        elif re.match(r"^__\w+__$", name) and name != "__version__":
            # allow some minimal introspection on private module
            # attrs to avoid unnecessary error-handling elsewhere
            return None
        else:
            # all other attribute access raises a helpful exception
            pfx = self._mod_pfx.get(self._module_name, "")
            msg = f"{pfx}{name} requires {self._module_name!r} module to be installed"
            raise ModuleNotFoundError(msg) from None


def _lazy_import(module_name: str) -> tuple[ModuleType, bool]:
    """
    Lazy import the given module; avoids up-front import costs.

    Parameters
    ----------
    module_name : str
        name of the module to import, eg: "pyarrow".

    Notes
    -----
    If the requested module is not available (eg: has not been installed), a proxy
    module is created in its place, which raises an exception on any attribute
    access. This allows for import and use as normal, without requiring explicit
    guard conditions - if the module is never used, no exception occurs; if it
    is, then a helpful exception is raised.

    Returns
    -------
    tuple of (Module, bool)
        A lazy-loading module and a boolean indicating if the requested/underlying
        module exists (if not, the returned module is a proxy).
    """
    # check if module is LOADED
    if module_name in sys.modules:
        return sys.modules[module_name], True

    # check if module is AVAILABLE
    try:
        module_spec = find_spec(module_name)
        module_available = not (module_spec is None or module_spec.loader is None)
    except ModuleNotFoundError:
        module_available = False

    # create lazy/proxy module that imports the real one on first use
    # (or raises an explanatory ModuleNotFoundError if not available)
    return (
        _LazyModule(
            module_name=module_name,
            module_available=module_available,
        ),
        module_available,
    )


if TYPE_CHECKING:
    import dataclasses
    import html
    import json
    import pickle
    import subprocess

    import altair
    import boto3
    import deltalake
    import fsspec
    import gevent
    import great_tables
    import hypothesis
    import numpy
    import pandas
    import polars_cloud
    import pyarrow
    import pydantic
    import pyiceberg
    import pyiceberg.schema
    import pytz
    import torch

else:
    # infrequently-used builtins
    dataclasses, _ = _lazy_import("dataclasses")
    html, _ = _lazy_import("html")
    json, _ = _lazy_import("json")
    pickle, _ = _lazy_import("pickle")
    subprocess, _ = _lazy_import("subprocess")

    # heavy/optional third party libs
    altair, _ALTAIR_AVAILABLE = _lazy_import("altair")
    boto3, _BOTO3_AVAILABLE = _lazy_import("boto3")
    deltalake, _DELTALAKE_AVAILABLE = _lazy_import("deltalake")
    fsspec, _FSSPEC_AVAILABLE = _lazy_import("fsspec")
    gevent, _GEVENT_AVAILABLE = _lazy_import("gevent")
    great_tables, _GREAT_TABLES_AVAILABLE = _lazy_import("great_tables")
    hypothesis, _HYPOTHESIS_AVAILABLE = _lazy_import("hypothesis")
    numpy, _NUMPY_AVAILABLE = _lazy_import("numpy")
    pandas, _PANDAS_AVAILABLE = _lazy_import("pandas")
    polars_cloud, _POLARS_CLOUD_AVAILABLE = _lazy_import("polars_cloud")
    pyarrow, _PYARROW_AVAILABLE = _lazy_import("pyarrow")
    pydantic, _PYDANTIC_AVAILABLE = _lazy_import("pydantic")
    pyiceberg, _PYICEBERG_AVAILABLE = _lazy_import("pyiceberg")
    torch, _TORCH_AVAILABLE = _lazy_import("torch")
    pytz, _PYTZ_AVAILABLE = _lazy_import("pytz")


@cache
def _might_be(cls: type, type_: str) -> bool:
    # infer whether the given class "might" be associated with the given
    # module (in which case it's reasonable to do a real isinstance check;
    # we defer that so as not to unnecessarily trigger module import)
    try:
        return any(f"{type_}." in str(o) for o in cls.mro())
    except TypeError:
        return False


def _check_for_numpy(obj: Any, *, check_type: bool = True) -> bool:
    return _NUMPY_AVAILABLE and _might_be(
        cast(Hashable, type(obj) if check_type else obj), "numpy"
    )


def _check_for_pandas(obj: Any, *, check_type: bool = True) -> bool:
    return _PANDAS_AVAILABLE and _might_be(
        cast(Hashable, type(obj) if check_type else obj), "pandas"
    )


def _check_for_pyarrow(obj: Any, *, check_type: bool = True) -> bool:
    return _PYARROW_AVAILABLE and _might_be(
        cast(Hashable, type(obj) if check_type else obj), "pyarrow"
    )


def _check_for_pydantic(obj: Any, *, check_type: bool = True) -> bool:
    return _PYDANTIC_AVAILABLE and _might_be(
        cast(Hashable, type(obj) if check_type else obj), "pydantic"
    )


def _check_for_torch(obj: Any, *, check_type: bool = True) -> bool:
    return _TORCH_AVAILABLE and _might_be(
        cast(Hashable, type(obj) if check_type else obj), "torch"
    )


def _check_for_pytz(obj: Any, *, check_type: bool = True) -> bool:
    return _PYTZ_AVAILABLE and _might_be(
        cast(Hashable, type(obj) if check_type else obj), "pytz"
    )


def import_optional(
    module_name: str,
    err_prefix: str = "required package",
    err_suffix: str = "not found",
    min_version: str | tuple[int, ...] | None = None,
    min_err_prefix: str = "requires",
    install_message: str | None = None,
) -> Any:
    """
    Import an optional dependency, returning the module.

    Parameters
    ----------
    module_name : str
        Name of the dependency to import.
    err_prefix : str, optional
        Error prefix to use in the raised exception (appears before the module name).
    err_suffix: str, optional
        Error suffix to use in the raised exception (follows the module name).
    min_version : {str, tuple[int]}, optional
        If a minimum module version is required, specify it here.
    min_err_prefix : str, optional
        Override the standard "requires" prefix for the minimum version error message.
    install_message : str, optional
        Override the standard "Please install it using..." exception message fragment.

    Examples
    --------
    >>> from polars._dependencies import import_optional
    >>> import_optional(
    ...     "definitely_a_real_module",
    ...     err_prefix="super-important package",
    ... )  # doctest: +SKIP
    ImportError: super-important package 'definitely_a_real_module' not installed.
    Please install it using the command `pip install definitely_a_real_module`.
    """
    from polars._utils.various import parse_version
    from polars.exceptions import ModuleUpgradeRequiredError

    module_root = module_name.split(".", 1)[0]
    try:
        module = import_module(module_name)
    except ImportError:
        prefix = f"{err_prefix.strip(' ')} " if err_prefix else ""
        suffix = f" {err_suffix.strip(' ')}" if err_suffix else ""
        err_message = f"{prefix}'{module_name}'{suffix}.\n" + (
            install_message
            or f"Please install using the command `pip install {module_root}`."
        )
        raise ModuleNotFoundError(err_message) from None

    if min_version:
        min_version = parse_version(min_version)
        mod_version = parse_version(module.__version__)
        if mod_version < min_version:
            msg = (
                f"{min_err_prefix} {module_root} "
                f"{'.'.join(str(v) for v in min_version)} or higher"
                f" (found {'.'.join(str(v) for v in mod_version)})"
            )
            raise ModuleUpgradeRequiredError(msg)

    return module


__all__ = [
    # lazy-load rarely-used/heavy builtins (for fast startup)
    "dataclasses",
    "html",
    "json",
    "pickle",
    "subprocess",
    # lazy-load third party libs
    "altair",
    "boto3",
    "deltalake",
    "fsspec",
    "gevent",
    "great_tables",
    "numpy",
    "pandas",
    "polars_cloud",
    "pydantic",
    "pyiceberg",
    "pyarrow",
    "torch",
    "pytz",
    # lazy utilities
    "_check_for_numpy",
    "_check_for_pandas",
    "_check_for_pyarrow",
    "_check_for_pydantic",
    "_check_for_torch",
    "_check_for_pytz",
    # exported flags/guards
    "_ALTAIR_AVAILABLE",
    "_DELTALAKE_AVAILABLE",
    "_FSSPEC_AVAILABLE",
    "_GEVENT_AVAILABLE",
    "_GREAT_TABLES_AVAILABLE",
    "_HYPOTHESIS_AVAILABLE",
    "_NUMPY_AVAILABLE",
    "_PANDAS_AVAILABLE",
    "_POLARS_CLOUD_AVAILABLE",
    "_PYARROW_AVAILABLE",
    "_PYDANTIC_AVAILABLE",
    "_PYICEBERG_AVAILABLE",
    "_TORCH_AVAILABLE",
]
