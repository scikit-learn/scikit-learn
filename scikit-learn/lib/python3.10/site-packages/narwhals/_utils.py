from __future__ import annotations

import os
import re
import sys
from collections.abc import Collection, Container, Iterable, Iterator, Mapping, Sequence
from datetime import timezone
from enum import Enum, auto
from functools import cache, lru_cache, partial, wraps
from importlib.util import find_spec
from inspect import getattr_static, getdoc
from itertools import chain
from operator import attrgetter
from secrets import token_hex
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Final,
    Generic,
    Literal,
    Protocol,
    TypeVar,
    Union,
    cast,
    overload,
)

from narwhals._enum import NoAutoEnum
from narwhals._exceptions import issue_deprecation_warning
from narwhals._typing_compat import assert_never, deprecated
from narwhals.dependencies import (
    get_cudf,
    get_dask_dataframe,
    get_duckdb,
    get_ibis,
    get_modin,
    get_pandas,
    get_polars,
    get_pyarrow,
    get_pyspark_connect,
    get_pyspark_sql,
    get_sqlframe,
    is_narwhals_series,
    is_narwhals_series_int,
    is_numpy_array_1d,
    is_numpy_array_1d_int,
    is_pandas_like_dataframe,
    is_pandas_like_series,
)
from narwhals.exceptions import ColumnNotFoundError, DuplicateError, InvalidOperationError

if TYPE_CHECKING:
    from collections.abc import Set  # noqa: PYI025
    from types import ModuleType

    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import (
        Concatenate,
        LiteralString,
        ParamSpec,
        Self,
        TypeAlias,
        TypeIs,
    )

    from narwhals._compliant import (
        CompliantExpr,
        CompliantExprT,
        CompliantFrameT,
        CompliantSeriesOrNativeExprT_co,
        CompliantSeriesT,
        NativeSeriesT_co,
    )
    from narwhals._compliant.any_namespace import NamespaceAccessor
    from narwhals._compliant.typing import (
        Accessor,
        EvalNames,
        NativeDataFrameT,
        NativeLazyFrameT,
    )
    from narwhals._namespace import Namespace
    from narwhals._native import (
        NativeArrow,
        NativeCuDF,
        NativeDask,
        NativeDuckDB,
        NativeIbis,
        NativeModin,
        NativePandas,
        NativePandasLike,
        NativePolars,
        NativePySpark,
        NativePySparkConnect,
        NativeSQLFrame,
    )
    from narwhals._translate import ArrowStreamExportable, IntoArrowTable, ToNarwhalsT_co
    from narwhals._typing import (
        Backend,
        IntoBackend,
        _ArrowImpl,
        _CuDFImpl,
        _DaskImpl,
        _DuckDBImpl,
        _EagerAllowedImpl,
        _IbisImpl,
        _LazyAllowedImpl,
        _LazyFrameCollectImpl,
        _ModinImpl,
        _PandasImpl,
        _PandasLikeImpl,
        _PolarsImpl,
        _PySparkConnectImpl,
        _PySparkImpl,
        _SQLFrameImpl,
    )
    from narwhals.dataframe import DataFrame, LazyFrame
    from narwhals.dtypes import DType
    from narwhals.series import Series
    from narwhals.typing import (
        CompliantDataFrame,
        CompliantLazyFrame,
        CompliantSeries,
        DTypes,
        FileSource,
        IntoSeriesT,
        MultiIndexSelector,
        SingleIndexSelector,
        SizedMultiIndexSelector,
        SizeUnit,
        SupportsNativeNamespace,
        TimeUnit,
        _1DArray,
        _SliceIndex,
        _SliceName,
        _SliceNone,
    )

    UnknownBackendName: TypeAlias = str

    FrameOrSeriesT = TypeVar(
        "FrameOrSeriesT", bound=Union[LazyFrame[Any], DataFrame[Any], Series[Any]]
    )

    _T1 = TypeVar("_T1")
    _T2 = TypeVar("_T2")
    _T3 = TypeVar("_T3")
    _Fn = TypeVar("_Fn", bound="Callable[..., Any]")
    P = ParamSpec("P")
    R = TypeVar("R")
    R1 = TypeVar("R1")
    R2 = TypeVar("R2")

    class _SupportsVersion(Protocol):
        __version__: str

    class _SupportsGet(Protocol):  # noqa: PYI046
        def __get__(self, instance: Any, owner: Any | None = None, /) -> Any: ...

    class _StoresColumns(Protocol):
        @property
        def columns(self) -> Sequence[str]: ...


_T = TypeVar("_T")
NativeT_co = TypeVar("NativeT_co", covariant=True)
CompliantT_co = TypeVar("CompliantT_co", covariant=True)
_IntoContext: TypeAlias = "_FullContext | NamespaceAccessor[_FullContext]"
_IntoContextT = TypeVar("_IntoContextT", bound=_IntoContext)
_Method: TypeAlias = "Callable[Concatenate[_IntoContextT, P], R]"
_Constructor: TypeAlias = "Callable[Concatenate[_T, P], R2]"


class _StoresNative(Protocol[NativeT_co]):
    """Provides access to a native object.

    Native objects have types like:

    >>> from pandas import Series
    >>> from pyarrow import Table
    """

    @property
    def native(self) -> NativeT_co:
        """Return the native object."""
        ...


class _StoresCompliant(Protocol[CompliantT_co]):  # noqa: PYI046
    """Provides access to a compliant object.

    Compliant objects have types like:

    >>> from narwhals._pandas_like.series import PandasLikeSeries
    >>> from narwhals._arrow.dataframe import ArrowDataFrame
    """

    @property
    def compliant(self) -> CompliantT_co:
        """Return the compliant object."""
        ...


class _StoresBackendVersion(Protocol):
    @property
    def _backend_version(self) -> tuple[int, ...]:
        """Version tuple for a native package."""
        ...


class _StoresVersion(Protocol):
    _version: Version
    """Narwhals API version (V1 or MAIN)."""


class _StoresImplementation(Protocol):
    _implementation: Implementation
    """Implementation of native object (pandas, Polars, PyArrow, ...)."""


class _LimitedContext(_StoresImplementation, _StoresVersion, Protocol):
    """Provides 2 attributes.

    - `_implementation`
    - `_version`
    """


class _FullContext(_StoresImplementation, _StoresBackendVersion, Protocol):
    """Provides 2 attributes.

    - `_implementation`
    - `_backend_version`
    """


class ValidateBackendVersion(_StoresImplementation, Protocol):
    """Ensure the target `Implementation` is on a supported version."""

    def _validate_backend_version(self) -> None:
        """Raise if installed version below `nw._utils.MIN_VERSIONS`.

        **Only use this when moving between backends.**
        Otherwise, the validation will have taken place already.
        """
        _ = self._implementation._backend_version()


class Version(Enum):
    V1 = auto()
    V2 = auto()
    MAIN = auto()

    @property
    def namespace(self) -> type[Namespace[Any]]:
        if self is Version.V1:
            from narwhals.stable.v1._namespace import Namespace as NamespaceV1

            return NamespaceV1
        if self is Version.V2:
            from narwhals.stable.v2._namespace import Namespace as NamespaceV2

            return NamespaceV2
        from narwhals._namespace import Namespace

        return Namespace

    @property
    def dtypes(self) -> DTypes:
        if self is Version.V1:
            from narwhals.stable.v1 import dtypes as dtypes_v1

            return dtypes_v1
        if self is Version.V2:
            from narwhals.stable.v2 import dtypes as dtypes_v2

            return dtypes_v2
        from narwhals import dtypes

        return dtypes

    @property
    def dataframe(self) -> type[DataFrame[Any]]:
        if self is Version.V1:
            from narwhals.stable.v1 import DataFrame as DataFrameV1

            return DataFrameV1
        if self is Version.V2:
            from narwhals.stable.v2 import DataFrame as DataFrameV2

            return DataFrameV2
        from narwhals.dataframe import DataFrame

        return DataFrame

    @property
    def lazyframe(self) -> type[LazyFrame[Any]]:
        if self is Version.V1:
            from narwhals.stable.v1 import LazyFrame as LazyFrameV1

            return LazyFrameV1
        if self is Version.V2:
            from narwhals.stable.v2 import LazyFrame as LazyFrameV2

            return LazyFrameV2
        from narwhals.dataframe import LazyFrame

        return LazyFrame

    @property
    def series(self) -> type[Series[Any]]:
        if self is Version.V1:
            from narwhals.stable.v1 import Series as SeriesV1

            return SeriesV1
        if self is Version.V2:
            from narwhals.stable.v2 import Series as SeriesV2

            return SeriesV2
        from narwhals.series import Series

        return Series


class Implementation(NoAutoEnum):
    """Implementation of native object (pandas, Polars, PyArrow, ...)."""

    PANDAS = "pandas"
    """pandas implementation."""
    MODIN = "modin"
    """Modin implementation."""
    CUDF = "cudf"
    """cuDF implementation."""
    PYARROW = "pyarrow"
    """PyArrow implementation."""
    PYSPARK = "pyspark"
    """PySpark implementation."""
    POLARS = "polars"
    """Polars implementation."""
    DASK = "dask"
    """Dask implementation."""
    DUCKDB = "duckdb"
    """DuckDB implementation."""
    IBIS = "ibis"
    """Ibis implementation."""
    SQLFRAME = "sqlframe"
    """SQLFrame implementation."""
    PYSPARK_CONNECT = "pyspark[connect]"
    """PySpark Connect implementation."""
    UNKNOWN = "unknown"
    """Unknown implementation."""

    def __str__(self) -> str:
        return str(self.value)

    @classmethod
    def from_native_namespace(
        cls: type[Self], native_namespace: ModuleType
    ) -> Implementation:  # pragma: no cover
        """Instantiate Implementation object from a native namespace module.

        Arguments:
            native_namespace: Native namespace.
        """
        mapping = {
            get_pandas(): Implementation.PANDAS,
            get_modin(): Implementation.MODIN,
            get_cudf(): Implementation.CUDF,
            get_pyarrow(): Implementation.PYARROW,
            get_pyspark_sql(): Implementation.PYSPARK,
            get_polars(): Implementation.POLARS,
            get_dask_dataframe(): Implementation.DASK,
            get_duckdb(): Implementation.DUCKDB,
            get_ibis(): Implementation.IBIS,
            get_sqlframe(): Implementation.SQLFRAME,
            get_pyspark_connect(): Implementation.PYSPARK_CONNECT,
        }
        return mapping.get(native_namespace, Implementation.UNKNOWN)

    @classmethod
    def from_string(cls: type[Self], backend_name: str) -> Implementation:
        """Instantiate Implementation object from a native namespace module.

        Arguments:
            backend_name: Name of backend, expressed as string.
        """
        try:
            return cls(backend_name)
        except ValueError:
            return Implementation.UNKNOWN

    @classmethod
    def from_backend(
        cls: type[Self], backend: IntoBackend[Backend] | UnknownBackendName
    ) -> Implementation:
        """Instantiate from native namespace module, string, or Implementation.

        Arguments:
            backend: Backend to instantiate Implementation from.
        """
        return (
            cls.from_string(backend)
            if isinstance(backend, str)
            else backend
            if isinstance(backend, Implementation)
            else cls.from_native_namespace(backend)
        )

    def to_native_namespace(self) -> ModuleType:
        """Return the native namespace module corresponding to Implementation."""
        if self is Implementation.UNKNOWN:
            msg = "Cannot return native namespace from UNKNOWN Implementation"
            raise AssertionError(msg)

        self._backend_version()
        module_name = _IMPLEMENTATION_TO_MODULE_NAME.get(self, self.value)
        return _import_native_namespace(module_name)

    def is_pandas(self) -> bool:
        """Return whether implementation is pandas.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_pandas()
            True
        """
        return self is Implementation.PANDAS

    def is_pandas_like(self) -> bool:
        """Return whether implementation is pandas, Modin, or cuDF.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_pandas_like()
            True
        """
        return self in {Implementation.PANDAS, Implementation.MODIN, Implementation.CUDF}

    def is_spark_like(self) -> bool:
        """Return whether implementation is pyspark or sqlframe.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_spark_like()
            False
        """
        return self in {
            Implementation.PYSPARK,
            Implementation.SQLFRAME,
            Implementation.PYSPARK_CONNECT,
        }

    def is_polars(self) -> bool:
        """Return whether implementation is Polars.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_polars()
            True
        """
        return self is Implementation.POLARS

    def is_cudf(self) -> bool:
        """Return whether implementation is cuDF.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_cudf()
            False
        """
        return self is Implementation.CUDF  # pragma: no cover

    def is_modin(self) -> bool:
        """Return whether implementation is Modin.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_modin()
            False
        """
        return self is Implementation.MODIN  # pragma: no cover

    def is_pyspark(self) -> bool:
        """Return whether implementation is PySpark.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_pyspark()
            False
        """
        return self is Implementation.PYSPARK  # pragma: no cover

    def is_pyspark_connect(self) -> bool:
        """Return whether implementation is PySpark.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_pyspark_connect()
            False
        """
        return self is Implementation.PYSPARK_CONNECT  # pragma: no cover

    def is_pyarrow(self) -> bool:
        """Return whether implementation is PyArrow.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_pyarrow()
            False
        """
        return self is Implementation.PYARROW  # pragma: no cover

    def is_dask(self) -> bool:
        """Return whether implementation is Dask.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_dask()
            False
        """
        return self is Implementation.DASK  # pragma: no cover

    def is_duckdb(self) -> bool:
        """Return whether implementation is DuckDB.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_duckdb()
            False
        """
        return self is Implementation.DUCKDB  # pragma: no cover

    def is_ibis(self) -> bool:
        """Return whether implementation is Ibis.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_ibis()
            False
        """
        return self is Implementation.IBIS  # pragma: no cover

    def is_sqlframe(self) -> bool:
        """Return whether implementation is SQLFrame.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.implementation.is_sqlframe()
            False
        """
        return self is Implementation.SQLFRAME  # pragma: no cover

    def _backend_version(self) -> tuple[int, ...]:
        """Returns backend version."""
        return backend_version(self)


MIN_VERSIONS: Mapping[Implementation, tuple[int, ...]] = {
    Implementation.PANDAS: (1, 1, 3),
    Implementation.MODIN: (0, 8, 2),
    Implementation.CUDF: (24, 10),
    Implementation.PYARROW: (13,),
    Implementation.PYSPARK: (3, 5),
    Implementation.PYSPARK_CONNECT: (3, 5),
    Implementation.POLARS: (0, 20, 4),
    Implementation.DASK: (2024, 8),
    Implementation.DUCKDB: (1, 1),
    Implementation.IBIS: (6,),
    Implementation.SQLFRAME: (3, 22, 0),
}

_IMPLEMENTATION_TO_MODULE_NAME: Mapping[Implementation, str] = {
    Implementation.DASK: "dask.dataframe",
    Implementation.MODIN: "modin.pandas",
    Implementation.PYSPARK: "pyspark.sql",
    Implementation.PYSPARK_CONNECT: "pyspark.sql.connect",
}
"""Stores non default mapping from Implementation to module name"""


@lru_cache(maxsize=16)
def _import_native_namespace(module_name: str) -> ModuleType:
    from importlib import import_module

    return import_module(module_name)


# NOTE: We can safely use an unbounded cache, the size is constrained by `len(Implementation._member_names_)`
# Faster than `lru_cache`
# https://docs.python.org/3/library/functools.html#functools.cache
@cache
def backend_version(implementation: Implementation, /) -> tuple[int, ...]:
    if not isinstance(implementation, Implementation):
        assert_never(implementation)
    if implementation is Implementation.UNKNOWN:
        return (0, 0, 0)
    into_version: ModuleType | str
    impl = implementation
    module_name = _IMPLEMENTATION_TO_MODULE_NAME.get(impl, impl.value)
    native_namespace = _import_native_namespace(module_name)
    if impl.is_sqlframe():
        import sqlframe._version

        into_version = sqlframe._version
    elif impl.is_pyspark() or impl.is_pyspark_connect():  # pragma: no cover
        import pyspark  # ignore-banned-import

        into_version = pyspark
    elif impl.is_dask():
        import dask  # ignore-banned-import

        into_version = dask
    else:
        into_version = native_namespace
    version = parse_version(into_version)
    if version < (min_version := MIN_VERSIONS[impl]):
        msg = f"Minimum version of {impl} supported by Narwhals is {min_version}, found: {version}"
        raise ValueError(msg)
    return version


def flatten(args: Any) -> list[Any]:
    return list(args[0] if (len(args) == 1 and _is_iterable(args[0])) else args)


def tupleify(arg: Any) -> Any:
    if not isinstance(arg, (list, tuple)):  # pragma: no cover
        return (arg,)
    return arg


def _is_iterable(arg: Any | Iterable[Any]) -> bool:
    from narwhals.series import Series

    if (
        (pd := get_pandas()) is not None and isinstance(arg, (pd.Series, pd.DataFrame))
    ) or (
        (pl := get_polars()) is not None
        and isinstance(arg, (pl.Series, pl.Expr, pl.DataFrame, pl.LazyFrame))
    ):
        # Non-exhaustive check for common potential mistakes.
        msg = (
            f"Expected Narwhals class or scalar, got: {qualified_type_name(arg)!r}.\n\n"
            "Hint: Perhaps you\n"
            "- forgot a `nw.from_native` somewhere?\n"
            "- used `pl.col` instead of `nw.col`?"
        )
        raise TypeError(msg)

    return isinstance(arg, Iterable) and not isinstance(arg, (str, bytes, Series))


def is_iterator(val: Iterable[_T] | Any) -> TypeIs[Iterator[_T]]:
    return isinstance(val, Iterator)


def parse_version(version: str | ModuleType | _SupportsVersion) -> tuple[int, ...]:
    """Simple version parser; split into a tuple of ints for comparison.

    Arguments:
        version: Version string, or object with one, to parse.
    """
    # lifted from Polars
    # [marco]: Take care of DuckDB pre-releases which end with e.g. `-dev4108`
    # and pandas pre-releases which end with e.g. .dev0+618.gb552dc95c9
    version_str = version if isinstance(version, str) else version.__version__
    version_str = re.sub(r"(\D?dev.*$)", "", version_str)
    return tuple(int(re.sub(r"\D", "", v)) for v in version_str.split("."))


@overload
def isinstance_or_issubclass(
    obj_or_cls: type, cls_or_tuple: type[_T]
) -> TypeIs[type[_T]]: ...


@overload
def isinstance_or_issubclass(
    obj_or_cls: object | type, cls_or_tuple: type[_T]
) -> TypeIs[_T | type[_T]]: ...


@overload
def isinstance_or_issubclass(
    obj_or_cls: type, cls_or_tuple: tuple[type[_T1], type[_T2]]
) -> TypeIs[type[_T1 | _T2]]: ...


@overload
def isinstance_or_issubclass(
    obj_or_cls: object | type, cls_or_tuple: tuple[type[_T1], type[_T2]]
) -> TypeIs[_T1 | _T2 | type[_T1 | _T2]]: ...


@overload
def isinstance_or_issubclass(
    obj_or_cls: type, cls_or_tuple: tuple[type[_T1], type[_T2], type[_T3]]
) -> TypeIs[type[_T1 | _T2 | _T3]]: ...


@overload
def isinstance_or_issubclass(
    obj_or_cls: object | type, cls_or_tuple: tuple[type[_T1], type[_T2], type[_T3]]
) -> TypeIs[_T1 | _T2 | _T3 | type[_T1 | _T2 | _T3]]: ...


@overload
def isinstance_or_issubclass(
    obj_or_cls: Any, cls_or_tuple: tuple[type, ...]
) -> TypeIs[Any]: ...


def isinstance_or_issubclass(obj_or_cls: Any, cls_or_tuple: Any) -> bool:
    from narwhals.dtypes import DType

    if isinstance(obj_or_cls, DType):
        return isinstance(obj_or_cls, cls_or_tuple)
    return isinstance(obj_or_cls, cls_or_tuple) or (
        isinstance(obj_or_cls, type) and issubclass(obj_or_cls, cls_or_tuple)
    )


def validate_laziness(items: Iterable[Any]) -> None:
    from narwhals.dataframe import DataFrame, LazyFrame

    if all(isinstance(item, DataFrame) for item in items) or (
        all(isinstance(item, LazyFrame) for item in items)
    ):
        return
    msg = f"The items to concatenate should either all be eager, or all lazy, got: {[type(item) for item in items]}"
    raise TypeError(msg)


def maybe_align_index(
    lhs: FrameOrSeriesT, rhs: Series[Any] | DataFrame[Any] | LazyFrame[Any]
) -> FrameOrSeriesT:
    """Align `lhs` to the Index of `rhs`, if they're both pandas-like.

    Arguments:
        lhs: Dataframe or Series.
        rhs: Dataframe or Series to align with.

    Notes:
        This is only really intended for backwards-compatibility purposes,
        for example if your library already aligns indices for users.
        If you're designing a new library, we highly encourage you to not
        rely on the Index.
        For non-pandas-like inputs, this only checks that `lhs` and `rhs`
        are the same length.

    Examples:
        >>> import pandas as pd
        >>> import polars as pl
        >>> import narwhals as nw
        >>> df_pd = pd.DataFrame({"a": [1, 2]}, index=[3, 4])
        >>> s_pd = pd.Series([6, 7], index=[4, 3])
        >>> df = nw.from_native(df_pd)
        >>> s = nw.from_native(s_pd, series_only=True)
        >>> nw.to_native(nw.maybe_align_index(df, s))
           a
        4  2
        3  1
    """
    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals._pandas_like.series import PandasLikeSeries

    def _validate_index(index: Any) -> None:
        if not index.is_unique:
            msg = "given index doesn't have a unique index"
            raise ValueError(msg)

    lhs_any = cast("Any", lhs)
    rhs_any = cast("Any", rhs)
    if isinstance(
        getattr(lhs_any, "_compliant_frame", None), PandasLikeDataFrame
    ) and isinstance(getattr(rhs_any, "_compliant_frame", None), PandasLikeDataFrame):
        _validate_index(lhs_any._compliant_frame.native.index)
        _validate_index(rhs_any._compliant_frame.native.index)
        return lhs_any._with_compliant(
            lhs_any._compliant_frame._with_native(
                lhs_any._compliant_frame.native.loc[rhs_any._compliant_frame.native.index]
            )
        )
    if isinstance(
        getattr(lhs_any, "_compliant_frame", None), PandasLikeDataFrame
    ) and isinstance(getattr(rhs_any, "_compliant_series", None), PandasLikeSeries):
        _validate_index(lhs_any._compliant_frame.native.index)
        _validate_index(rhs_any._compliant_series.native.index)
        return lhs_any._with_compliant(
            lhs_any._compliant_frame._with_native(
                lhs_any._compliant_frame.native.loc[
                    rhs_any._compliant_series.native.index
                ]
            )
        )
    if isinstance(
        getattr(lhs_any, "_compliant_series", None), PandasLikeSeries
    ) and isinstance(getattr(rhs_any, "_compliant_frame", None), PandasLikeDataFrame):
        _validate_index(lhs_any._compliant_series.native.index)
        _validate_index(rhs_any._compliant_frame.native.index)
        return lhs_any._with_compliant(
            lhs_any._compliant_series._with_native(
                lhs_any._compliant_series.native.loc[
                    rhs_any._compliant_frame.native.index
                ]
            )
        )
    if isinstance(
        getattr(lhs_any, "_compliant_series", None), PandasLikeSeries
    ) and isinstance(getattr(rhs_any, "_compliant_series", None), PandasLikeSeries):
        _validate_index(lhs_any._compliant_series.native.index)
        _validate_index(rhs_any._compliant_series.native.index)
        return lhs_any._with_compliant(
            lhs_any._compliant_series._with_native(
                lhs_any._compliant_series.native.loc[
                    rhs_any._compliant_series.native.index
                ]
            )
        )
    if len(lhs_any) != len(rhs_any):
        msg = f"Expected `lhs` and `rhs` to have the same length, got {len(lhs_any)} and {len(rhs_any)}"
        raise ValueError(msg)
    return lhs


def maybe_get_index(obj: DataFrame[Any] | LazyFrame[Any] | Series[Any]) -> Any | None:
    """Get the index of a DataFrame or a Series, if it's pandas-like.

    Arguments:
        obj: Dataframe or Series.

    Notes:
        This is only really intended for backwards-compatibility purposes,
        for example if your library already aligns indices for users.
        If you're designing a new library, we highly encourage you to not
        rely on the Index.
        For non-pandas-like inputs, this returns `None`.

    Examples:
        >>> import pandas as pd
        >>> import polars as pl
        >>> import narwhals as nw
        >>> df_pd = pd.DataFrame({"a": [1, 2], "b": [4, 5]})
        >>> df = nw.from_native(df_pd)
        >>> nw.maybe_get_index(df)
        RangeIndex(start=0, stop=2, step=1)
        >>> series_pd = pd.Series([1, 2])
        >>> series = nw.from_native(series_pd, series_only=True)
        >>> nw.maybe_get_index(series)
        RangeIndex(start=0, stop=2, step=1)
    """
    obj_any = cast("Any", obj)
    native_obj = obj_any.to_native()
    if is_pandas_like_dataframe(native_obj) or is_pandas_like_series(native_obj):
        return native_obj.index
    return None


def maybe_set_index(
    obj: FrameOrSeriesT,
    column_names: str | list[str] | None = None,
    *,
    index: Series[IntoSeriesT] | list[Series[IntoSeriesT]] | None = None,
) -> FrameOrSeriesT:
    """Set the index of a DataFrame or a Series, if it's pandas-like.

    Arguments:
        obj: object for which maybe set the index (can be either a Narwhals `DataFrame`
            or `Series`).
        column_names: name or list of names of the columns to set as index.
            For dataframes, only one of `column_names` and `index` can be specified but
            not both. If `column_names` is passed and `df` is a Series, then a
            `ValueError` is raised.
        index: series or list of series to set as index.

    Raises:
        ValueError: If one of the following conditions happens

            - none of `column_names` and `index` are provided
            - both `column_names` and `index` are provided
            - `column_names` is provided and `df` is a Series

    Notes:
        This is only really intended for backwards-compatibility purposes, for example if
        your library already aligns indices for users.
        If you're designing a new library, we highly encourage you to not
        rely on the Index.

        For non-pandas-like inputs, this is a no-op.

    Examples:
        >>> import pandas as pd
        >>> import polars as pl
        >>> import narwhals as nw
        >>> df_pd = pd.DataFrame({"a": [1, 2], "b": [4, 5]})
        >>> df = nw.from_native(df_pd)
        >>> nw.to_native(nw.maybe_set_index(df, "b"))  # doctest: +NORMALIZE_WHITESPACE
           a
        b
        4  1
        5  2
    """
    from narwhals.translate import to_native

    df_any = cast("Any", obj)
    native_obj = df_any.to_native()

    if column_names is not None and index is not None:
        msg = "Only one of `column_names` or `index` should be provided"
        raise ValueError(msg)

    if not column_names and index is None:
        msg = "Either `column_names` or `index` should be provided"
        raise ValueError(msg)

    if index is not None:
        keys = (
            [to_native(idx, pass_through=True) for idx in index]
            if _is_iterable(index)
            else to_native(index, pass_through=True)
        )
    else:
        keys = column_names

    if is_pandas_like_dataframe(native_obj):
        return df_any._with_compliant(
            df_any._compliant_frame._with_native(native_obj.set_index(keys))
        )
    if is_pandas_like_series(native_obj):
        from narwhals._pandas_like.utils import set_index

        if column_names:
            msg = "Cannot set index using column names on a Series"
            raise ValueError(msg)

        native_obj = set_index(
            native_obj,
            keys,
            implementation=obj._compliant_series._implementation,  # type: ignore[union-attr]
        )
        return df_any._with_compliant(df_any._compliant_series._with_native(native_obj))
    return df_any


def maybe_reset_index(obj: FrameOrSeriesT) -> FrameOrSeriesT:
    """Reset the index to the default integer index of a DataFrame or a Series, if it's pandas-like.

    Arguments:
        obj: Dataframe or Series.

    Notes:
        This is only really intended for backwards-compatibility purposes,
        for example if your library already resets the index for users.
        If you're designing a new library, we highly encourage you to not
        rely on the Index.
        For non-pandas-like inputs, this is a no-op.

    Examples:
        >>> import pandas as pd
        >>> import polars as pl
        >>> import narwhals as nw
        >>> df_pd = pd.DataFrame({"a": [1, 2], "b": [4, 5]}, index=([6, 7]))
        >>> df = nw.from_native(df_pd)
        >>> nw.to_native(nw.maybe_reset_index(df))
           a  b
        0  1  4
        1  2  5
        >>> series_pd = pd.Series([1, 2])
        >>> series = nw.from_native(series_pd, series_only=True)
        >>> nw.maybe_get_index(series)
        RangeIndex(start=0, stop=2, step=1)
    """
    obj_any = cast("Any", obj)
    native_obj = obj_any.to_native()
    if is_pandas_like_dataframe(native_obj):
        native_namespace = obj_any.__native_namespace__()
        if _has_default_index(native_obj, native_namespace):
            return obj_any
        return obj_any._with_compliant(
            obj_any._compliant_frame._with_native(native_obj.reset_index(drop=True))
        )
    if is_pandas_like_series(native_obj):
        native_namespace = obj_any.__native_namespace__()
        if _has_default_index(native_obj, native_namespace):
            return obj_any
        return obj_any._with_compliant(
            obj_any._compliant_series._with_native(native_obj.reset_index(drop=True))
        )
    return obj_any


if TYPE_CHECKING:
    zip_strict = partial(zip, strict=True)
else:
    import sys

    if sys.version_info >= (3, 10):
        zip_strict = partial(zip, strict=True)
    else:  # pragma: no cover
        # https://stackoverflow.com/questions/32954486/zip-iterators-asserting-for-equal-length-in-python/69485272#69485272

        def zip_strict(*iterables: Iterable[Any]) -> Iterable[tuple[Any, ...]]:
            # For trivial cases, use pure zip.
            if len(iterables) < 2:
                return zip(*iterables)
            # Tail for the first iterable
            first_stopped = False

            def first_tail() -> Any:
                nonlocal first_stopped
                first_stopped = True
                return
                yield

            # Tail for the zip
            def zip_tail() -> Any:
                if not first_stopped:  # pragma: no cover
                    msg = "zip_strict: first iterable is longer"
                    raise ValueError(msg)
                for _ in chain.from_iterable(rest):  # pragma: no cover
                    msg = "zip_strict: first iterable is shorter"
                    raise ValueError(msg)
                    yield

            # Put the pieces together
            iterables_it = iter(iterables)
            first = chain(next(iterables_it), first_tail())
            rest = list(map(iter, iterables_it))
            return chain(zip(first, *rest), zip_tail())


def _is_range_index(obj: Any, native_namespace: Any) -> TypeIs[pd.RangeIndex]:
    return isinstance(obj, native_namespace.RangeIndex)


def _has_default_index(
    native_frame_or_series: pd.Series[Any] | pd.DataFrame, native_namespace: Any
) -> bool:
    index = native_frame_or_series.index
    return (
        _is_range_index(index, native_namespace)
        and index.start == 0
        and index.stop == len(index)
        and index.step == 1
    )


def maybe_convert_dtypes(
    obj: FrameOrSeriesT, *args: bool, **kwargs: bool | str
) -> FrameOrSeriesT:
    """Convert columns or series to the best possible dtypes using dtypes supporting ``pd.NA``, if df is pandas-like.

    Arguments:
        obj: DataFrame or Series.
        *args: Additional arguments which gets passed through.
        **kwargs: Additional arguments which gets passed through.

    Notes:
        For non-pandas-like inputs, this is a no-op.
        Also, `args` and `kwargs` just get passed down to the underlying library as-is.

    Examples:
        >>> import pandas as pd
        >>> import polars as pl
        >>> import narwhals as nw
        >>> import numpy as np
        >>> df_pd = pd.DataFrame(
        ...     {
        ...         "a": pd.Series([1, 2, 3], dtype=np.dtype("int32")),
        ...         "b": pd.Series([True, False, np.nan], dtype=np.dtype("O")),
        ...     }
        ... )
        >>> df = nw.from_native(df_pd)
        >>> nw.to_native(
        ...     nw.maybe_convert_dtypes(df)
        ... ).dtypes  # doctest: +NORMALIZE_WHITESPACE
        a             Int32
        b           boolean
        dtype: object
    """
    if not obj.implementation.is_pandas_like():
        return obj
    result = obj._with_compliant(
        obj._compliant._with_native(obj.to_native().convert_dtypes(*args, **kwargs))
    )
    return cast("FrameOrSeriesT", result)


def scale_bytes(sz: int, unit: SizeUnit) -> int | float:
    """Scale size in bytes to other size units (eg: "kb", "mb", "gb", "tb").

    Arguments:
        sz: original size in bytes
        unit: size unit to convert into
    """
    if unit in {"b", "bytes"}:
        return sz
    if unit in {"kb", "kilobytes"}:
        return sz / 1024
    if unit in {"mb", "megabytes"}:
        return sz / 1024**2
    if unit in {"gb", "gigabytes"}:
        return sz / 1024**3
    if unit in {"tb", "terabytes"}:
        return sz / 1024**4
    msg = f"`unit` must be one of {{'b', 'kb', 'mb', 'gb', 'tb'}}, got {unit!r}"
    raise ValueError(msg)


def is_ordered_categorical(series: Series[Any]) -> bool:
    """Return whether indices of categories are semantically meaningful.

    This is a convenience function to accessing what would otherwise be
    the `is_ordered` property from the DataFrame Interchange Protocol,
    see https://data-apis.org/dataframe-protocol/latest/API.html.

    - For Polars:
      - Enums are always ordered.
      - Categoricals are ordered if `dtype.ordering == "physical"`.
    - For pandas-like APIs:
      - Categoricals are ordered if `dtype.cat.ordered == True`.
    - For PyArrow table:
      - Categoricals are ordered if `dtype.type.ordered == True`.

    Arguments:
        series: Input Series.

    Examples:
        >>> import narwhals as nw
        >>> import pandas as pd
        >>> import polars as pl
        >>> data = ["x", "y"]
        >>> s_pd = pd.Series(data, dtype=pd.CategoricalDtype(ordered=True))
        >>> s_pl = pl.Series(data, dtype=pl.Categorical(ordering="lexical"))

        Let's define a library-agnostic function:

        >>> @nw.narwhalify
        ... def func(s):
        ...     return nw.is_ordered_categorical(s)

        Then, we can pass any supported library to `func`:

        >>> func(s_pd)
        True
        >>> func(s_pl)
        False
    """
    from narwhals._interchange.series import InterchangeSeries

    dtypes = series._compliant_series._version.dtypes
    compliant = series._compliant_series
    # If it doesn't match any branches, let's just play it safe and return False.
    result: bool = False
    if isinstance(compliant, InterchangeSeries) and isinstance(
        series.dtype, dtypes.Categorical
    ):
        result = compliant.native.describe_categorical["is_ordered"]
    elif series.dtype == dtypes.Enum:
        result = True
    elif series.dtype != dtypes.Categorical:
        result = False
    else:
        native = series.to_native()
        impl = series.implementation
        if impl.is_polars() and impl._backend_version() < (1, 32):
            # NOTE: Deprecated https://github.com/pola-rs/polars/pull/23779
            # Since version 1.32.0, ordering parameter is ignored and
            # it always behaves as if 'lexical' was passed.
            result = cast("pl.Categorical", native.dtype).ordering == "physical"
        elif impl.is_pandas_like():
            result = bool(native.cat.ordered)
        elif impl.is_pyarrow():
            from narwhals._arrow.utils import is_dictionary

            result = is_dictionary(native.type) and native.type.ordered
    return result


def generate_unique_token(
    n_bytes: int, columns: Container[str], prefix: str = "nw"
) -> str:  # pragma: no cover
    msg = (
        "Use `generate_temporary_column_name` instead. `generate_unique_token` is "
        "deprecated and it will be removed in future versions"
    )
    issue_deprecation_warning(msg, _version="1.13.0")
    return generate_temporary_column_name(n_bytes=n_bytes, columns=columns, prefix=prefix)


def generate_temporary_column_name(
    n_bytes: int, columns: Container[str], prefix: str = "nw"
) -> str:
    """Generates a unique column name that is not present in the given list of columns.

    It relies on [python secrets token_hex](https://docs.python.org/3/library/secrets.html#secrets.token_hex)
    function to return a string nbytes random bytes.

    Arguments:
        n_bytes: The number of bytes to generate for the token.
        columns: The list of columns to check for uniqueness.
        prefix: prefix with which the temporary column name should start with.

    Returns:
        A unique token that is not present in the given list of columns.

    Raises:
        AssertionError: If a unique token cannot be generated after 100 attempts.

    Examples:
        >>> import narwhals as nw
        >>> columns = ["abc", "xyz"]
        >>> nw.generate_temporary_column_name(n_bytes=8, columns=columns) not in columns
        True
        >>> temp_name = nw.generate_temporary_column_name(
        ...     n_bytes=8, columns=columns, prefix="foo"
        ... )
        >>> temp_name not in columns and temp_name.startswith("foo")
        True
    """
    counter = 0
    while True:
        token = f"{prefix}{token_hex(n_bytes - 1)}"
        if token not in columns:
            return token

        counter += 1
        if counter > 100:
            msg = (
                "Internal Error: Narwhals was not able to generate a column name with "
                f"{n_bytes=} and not in {columns}"
            )
            raise AssertionError(msg)


def parse_columns_to_drop(
    frame: _StoresColumns, subset: Iterable[str], /, *, strict: bool
) -> list[str]:
    if not strict:
        return list(set(frame.columns).intersection(subset))
    to_drop = list(subset)
    if error := check_columns_exist(to_drop, available=frame.columns):
        raise error
    return to_drop


def is_sequence_but_not_str(sequence: Sequence[_T] | Any) -> TypeIs[Sequence[_T]]:
    return isinstance(sequence, Sequence) and not isinstance(sequence, str)


def is_slice_none(obj: Any) -> TypeIs[_SliceNone]:
    return isinstance(obj, slice) and obj == slice(None)


def is_sized_multi_index_selector(
    obj: Any,
) -> TypeIs[SizedMultiIndexSelector[Series[Any] | CompliantSeries[Any]]]:
    return (
        (
            is_sequence_but_not_str(obj)
            and ((len(obj) > 0 and isinstance(obj[0], int)) or (len(obj) == 0))
        )
        or is_numpy_array_1d_int(obj)
        or is_narwhals_series_int(obj)
        or is_compliant_series_int(obj)
    )


def is_sequence_like(
    obj: Sequence[_T] | Any,
) -> TypeIs[Sequence[_T] | Series[Any] | _1DArray]:
    return (
        is_sequence_but_not_str(obj)
        or is_numpy_array_1d(obj)
        or is_narwhals_series(obj)
        or is_compliant_series(obj)
    )


def is_slice_index(obj: Any) -> TypeIs[_SliceIndex]:
    return isinstance(obj, slice) and (
        isinstance(obj.start, int)
        or isinstance(obj.stop, int)
        or (isinstance(obj.step, int) and obj.start is None and obj.stop is None)
    )


def is_range(obj: Any) -> TypeIs[range]:
    return isinstance(obj, range)


def is_single_index_selector(obj: Any) -> TypeIs[SingleIndexSelector]:
    return bool(isinstance(obj, int) and not isinstance(obj, bool))


def is_index_selector(
    obj: Any,
) -> TypeIs[SingleIndexSelector | MultiIndexSelector[Series[Any] | CompliantSeries[Any]]]:
    return (
        is_single_index_selector(obj)
        or is_sized_multi_index_selector(obj)
        or is_slice_index(obj)
    )


def is_list_of(obj: Any, tp: type[_T]) -> TypeIs[list[_T]]:
    # Check if an object is a list of `tp`, only sniffing the first element.
    return bool(isinstance(obj, list) and obj and isinstance(obj[0], tp))


def predicates_contains_list_of_bool(
    predicates: Collection[Any],
) -> TypeIs[Collection[list[bool]]]:
    return any(is_list_of(pred, bool) for pred in predicates)


def is_sequence_of(obj: Any, tp: type[_T]) -> TypeIs[Sequence[_T]]:
    # Check if an object is a sequence of `tp`, only sniffing the first element.
    return bool(
        is_sequence_but_not_str(obj)
        and (first := next(iter(obj), None))
        and isinstance(first, tp)
    )


def validate_strict_and_pass_though(
    strict: bool | None,  # noqa: FBT001
    pass_through: bool | None,  # noqa: FBT001
    *,
    pass_through_default: bool,
) -> bool:
    if strict is None and pass_through is None:
        pass_through = pass_through_default
    elif strict is not None and pass_through is None:
        pass_through = not strict
    elif strict is None and pass_through is not None:
        pass
    else:
        msg = "Cannot pass both `strict` and `pass_through`"
        raise ValueError(msg)
    return pass_through


def deprecate_native_namespace(
    *, warn_version: str = "", required: bool = False
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Decorator to transition from `native_namespace` to `backend` argument.

    Arguments:
        warn_version: Emit a deprecation warning from this version.
        required: Raise when both `native_namespace`, `backend` are `None`.

    Returns:
        Wrapped function, with `native_namespace` **removed**.
    """

    def decorate(fn: Callable[P, R], /) -> Callable[P, R]:
        @wraps(fn)
        def wrapper(*args: P.args, **kwds: P.kwargs) -> R:
            backend = kwds.pop("backend", None)
            native_namespace = kwds.pop("native_namespace", None)
            if native_namespace is not None and backend is None:
                if warn_version:
                    msg = (
                        "`native_namespace` is deprecated, please use `backend` instead.\n\n"
                        "Note: `native_namespace` will remain available in `narwhals.stable.v1`.\n"
                        "See https://narwhals-dev.github.io/narwhals/backcompat/ for more information.\n"
                    )
                    issue_deprecation_warning(msg, _version=warn_version)
                backend = native_namespace
            elif native_namespace is not None and backend is not None:
                msg = "Can't pass both `native_namespace` and `backend`"
                raise ValueError(msg)
            elif native_namespace is None and backend is None and required:
                msg = f"`backend` must be specified in `{fn.__name__}`."
                raise ValueError(msg)
            kwds["backend"] = backend
            return fn(*args, **kwds)

        return wrapper

    return decorate


def _validate_rolling_arguments(
    window_size: int, min_samples: int | None
) -> tuple[int, int]:
    ensure_type(window_size, int, param_name="window_size")
    ensure_type(min_samples, int, type(None), param_name="min_samples")

    if window_size < 1:
        msg = "window_size must be greater or equal than 1"
        raise ValueError(msg)

    if min_samples is not None:
        if min_samples < 1:
            msg = "min_samples must be greater or equal than 1"
            raise ValueError(msg)

        if min_samples > window_size:
            msg = "`min_samples` must be less or equal than `window_size`"
            raise InvalidOperationError(msg)
    else:
        min_samples = window_size

    return window_size, min_samples


def generate_repr(header: str, native_repr: str) -> str:
    try:
        terminal_width = os.get_terminal_size().columns
    except OSError:
        terminal_width = int(os.getenv("COLUMNS", 80))  # noqa: PLW1508
    native_lines = native_repr.expandtabs().splitlines()
    max_native_width = max(len(line) for line in native_lines)

    if max_native_width + 2 <= terminal_width:
        length = max(max_native_width, len(header))
        output = f"┌{'─' * length}┐\n"
        header_extra = length - len(header)
        output += f"|{' ' * (header_extra // 2)}{header}{' ' * (header_extra // 2 + header_extra % 2)}|\n"
        output += f"|{'-' * (length)}|\n"
        start_extra = (length - max_native_width) // 2
        end_extra = (length - max_native_width) // 2 + (length - max_native_width) % 2
        for line in native_lines:
            output += f"|{' ' * (start_extra)}{line}{' ' * (end_extra + max_native_width - len(line))}|\n"
        output += f"└{'─' * length}┘"
        return output

    diff = 39 - len(header)
    return (
        f"┌{'─' * (39)}┐\n"
        f"|{' ' * (diff // 2)}{header}{' ' * (diff // 2 + diff % 2)}|\n"
        "| Use `.to_native` to see native output |\n└"
        f"{'─' * 39}┘"
    )


def check_columns_exist(
    subset: Collection[str], /, *, available: Collection[str]
) -> ColumnNotFoundError | None:
    if missing := set(subset).difference(available):
        return ColumnNotFoundError.from_missing_and_available_column_names(
            missing, available
        )
    return None


def check_column_names_are_unique(columns: Collection[str]) -> None:
    if len(columns) != len(set(columns)):
        from collections import Counter

        counter = Counter(columns)
        duplicates = {k: v for k, v in counter.items() if v > 1}
        msg = "".join(f"\n- '{k}' {v} times" for k, v in duplicates.items())
        msg = f"Expected unique column names, got:{msg}"
        raise DuplicateError(msg)


def _parse_time_unit_and_time_zone(
    time_unit: TimeUnit | Iterable[TimeUnit] | None,
    time_zone: str | timezone | Iterable[str | timezone | None] | None,
) -> tuple[Set[TimeUnit], Set[str | None]]:
    time_units: Set[TimeUnit] = (
        {"ms", "us", "ns", "s"}
        if time_unit is None
        else {time_unit}
        if isinstance(time_unit, str)
        else set(time_unit)
    )
    time_zones: Set[str | None] = (
        {None}
        if time_zone is None
        else {str(time_zone)}
        if isinstance(time_zone, (str, timezone))
        else {str(tz) if tz is not None else None for tz in time_zone}
    )
    return time_units, time_zones


def dtype_matches_time_unit_and_time_zone(
    dtype: DType, dtypes: DTypes, time_units: Set[TimeUnit], time_zones: Set[str | None]
) -> bool:
    return (
        isinstance(dtype, dtypes.Datetime)
        and (dtype.time_unit in time_units)
        and (
            dtype.time_zone in time_zones
            or ("*" in time_zones and dtype.time_zone is not None)
        )
    )


def get_column_names(frame: _StoresColumns, /) -> Sequence[str]:
    return frame.columns


def exclude_column_names(frame: _StoresColumns, names: Container[str]) -> Sequence[str]:
    return [col_name for col_name in frame.columns if col_name not in names]


def passthrough_column_names(names: Sequence[str], /) -> EvalNames[Any]:
    def fn(_frame: Any, /) -> Sequence[str]:
        return names

    return fn


_SENTINEL: Final = object()


def _hasattr_static(obj: Any, attr: str) -> bool:
    return getattr_static(obj, attr, _SENTINEL) is not _SENTINEL


def is_compliant_dataframe(
    obj: CompliantDataFrame[
        CompliantSeriesT, CompliantExprT, NativeDataFrameT, ToNarwhalsT_co
    ]
    | Any,
) -> TypeIs[
    CompliantDataFrame[CompliantSeriesT, CompliantExprT, NativeDataFrameT, ToNarwhalsT_co]
]:
    return _hasattr_static(obj, "__narwhals_dataframe__")


def is_compliant_lazyframe(
    obj: CompliantLazyFrame[CompliantExprT, NativeLazyFrameT, ToNarwhalsT_co] | Any,
) -> TypeIs[CompliantLazyFrame[CompliantExprT, NativeLazyFrameT, ToNarwhalsT_co]]:
    return _hasattr_static(obj, "__narwhals_lazyframe__")


def is_compliant_series(
    obj: CompliantSeries[NativeSeriesT_co] | Any,
) -> TypeIs[CompliantSeries[NativeSeriesT_co]]:
    return _hasattr_static(obj, "__narwhals_series__")


def is_compliant_series_int(
    obj: CompliantSeries[NativeSeriesT_co] | Any,
) -> TypeIs[CompliantSeries[NativeSeriesT_co]]:
    return is_compliant_series(obj) and obj.dtype.is_integer()


def is_compliant_expr(
    obj: CompliantExpr[CompliantFrameT, CompliantSeriesOrNativeExprT_co] | Any,
) -> TypeIs[CompliantExpr[CompliantFrameT, CompliantSeriesOrNativeExprT_co]]:
    return hasattr(obj, "__narwhals_expr__")


def _is_namespace_accessor(obj: _IntoContext) -> TypeIs[NamespaceAccessor[_FullContext]]:
    # NOTE: Only `compliant` has false positives **internally**
    # - https://github.com/narwhals-dev/narwhals/blob/cc69bac35eb8c81a1106969c49bfba9fd569b856/narwhals/_compliant/group_by.py#L44-L49
    # - https://github.com/narwhals-dev/narwhals/blob/cc69bac35eb8c81a1106969c49bfba9fd569b856/narwhals/_namespace.py#L166-L168
    # NOTE: Only `_accessor` has false positives **upstream**
    # - https://github.com/pandas-dev/pandas/blob/e209a35403f8835bbcff97636b83d2fc39b51e68/pandas/core/accessor.py#L200-L233
    # - https://github.com/pola-rs/polars/blob/a60c5019f7b694c97009ef9208d25aaa4cc1d8a6/py-polars/polars/api.py#L29-L42
    return _hasattr_static(obj, "compliant") and _hasattr_static(obj, "_accessor")


def is_eager_allowed(impl: Implementation, /) -> TypeIs[_EagerAllowedImpl]:
    """Return True if `impl` allows eager operations."""
    return impl in {
        Implementation.CUDF,
        Implementation.MODIN,
        Implementation.PANDAS,
        Implementation.POLARS,
        Implementation.PYARROW,
    }


def can_lazyframe_collect(impl: Implementation, /) -> TypeIs[_LazyFrameCollectImpl]:
    """Return True if `LazyFrame.collect(impl)` is allowed."""
    return impl in {Implementation.PANDAS, Implementation.POLARS, Implementation.PYARROW}


def is_lazy_allowed(impl: Implementation, /) -> TypeIs[_LazyAllowedImpl]:
    """Return True if `DataFrame.lazy(impl)` is allowed."""
    return impl in {
        Implementation.DASK,
        Implementation.DUCKDB,
        Implementation.IBIS,
        Implementation.POLARS,
        Implementation.PYSPARK,
        Implementation.PYSPARK_CONNECT,
        Implementation.SQLFRAME,
    }


def has_native_namespace(obj: Any) -> TypeIs[SupportsNativeNamespace]:
    return _hasattr_static(obj, "__native_namespace__")


def supports_arrow_c_stream(obj: Any) -> TypeIs[ArrowStreamExportable]:
    return _hasattr_static(obj, "__arrow_c_stream__")


def _remap_full_join_keys(
    left_on: Collection[str], right_on: Collection[str], suffix: str
) -> dict[str, str]:
    """Remap join keys to avoid collisions.

    If left keys collide with the right keys, append the suffix.
    If there's no collision, let the right keys be.

    Arguments:
        left_on: Left keys.
        right_on: Right keys.
        suffix: Suffix to append to right keys.

    Returns:
        A map of old to new right keys.
    """
    right_keys_suffixed = (
        f"{key}{suffix}" if key in left_on else key for key in right_on
    )
    return dict(zip(right_on, right_keys_suffixed))


def _into_arrow_table(data: IntoArrowTable, context: _LimitedContext, /) -> pa.Table:
    """Guards `ArrowDataFrame.from_arrow` w/ safer imports.

    Arguments:
        data: Object which implements `__arrow_c_stream__`.
        context: Initialized compliant object.
    """
    if find_spec("pyarrow"):
        ns = context._version.namespace.from_backend("pyarrow").compliant
        return ns._dataframe.from_arrow(data, context=ns).native
    msg = f"'pyarrow>=14.0.0' is required for `from_arrow` for object of type {qualified_type_name(data)!r}."  # pragma: no cover
    raise ModuleNotFoundError(msg)  # pragma: no cover


# TODO @dangotbanned: Extend with runtime behavior for `v1.*`
# See `narwhals.exceptions.NarwhalsUnstableWarning`
def unstable(fn: _Fn, /) -> _Fn:
    """Visual-only marker for unstable functionality.

    Arguments:
        fn: Function to decorate.

    Returns:
        Decorated function (unchanged).

    Examples:
        >>> from narwhals._utils import unstable
        >>> @unstable
        ... def a_work_in_progress_feature(*args):
        ...     return args
        >>>
        >>> a_work_in_progress_feature.__name__
        'a_work_in_progress_feature'
        >>> a_work_in_progress_feature(1, 2, 3)
        (1, 2, 3)
    """
    return fn


def _is_naive_format(format: str) -> bool:
    """Determines if a datetime format string is 'naive', i.e., does not include timezone information.

    A format is considered naive if it does not contain any of the following

    - '%s': Unix timestamp
    - '%z': UTC offset
    - 'Z' : UTC timezone designator

    Arguments:
        format: The datetime format string to check.

    Returns:
        bool: True if the format is naive (does not include timezone info), False otherwise.
    """
    return not any(x in format for x in ("%s", "%z", "Z"))


class not_implemented:  # noqa: N801
    """Mark some functionality as unsupported.

    Arguments:
        alias: optional name used instead of the data model hook [`__set_name__`].

    Returns:
        An exception-raising [descriptor].

    Notes:
        - Attribute/method name *doesn't* need to be declared twice
        - Allows different behavior when looked up on the class vs instance
        - Allows us to use `isinstance(...)` instead of monkeypatching an attribute to the function

    Examples:
        >>> from narwhals._utils import not_implemented
        >>> class Thing:
        ...     def totally_ready(self) -> str:
        ...         return "I'm ready!"
        ...
        ...     not_ready_yet = not_implemented()
        >>>
        >>> thing = Thing()
        >>> thing.totally_ready()
        "I'm ready!"
        >>> thing.not_ready_yet()
        Traceback (most recent call last):
            ...
        NotImplementedError: 'not_ready_yet' is not implemented for: 'Thing'.
        ...
        >>> isinstance(Thing.not_ready_yet, not_implemented)
        True

    [`__set_name__`]: https://docs.python.org/3/reference/datamodel.html#object.__set_name__
    [descriptor]: https://docs.python.org/3/howto/descriptor.html
    """

    def __init__(self, alias: str | None = None, /) -> None:
        # NOTE: Don't like this
        # Trying to workaround `mypy` requiring `@property` everywhere
        self._alias: str | None = alias

    def __repr__(self) -> str:
        return f"<{type(self).__name__}>: {self._name_owner}.{self._name}"

    def __set_name__(self, owner: type[_T], name: str) -> None:
        # https://docs.python.org/3/howto/descriptor.html#customized-names
        self._name_owner: str = owner.__name__
        self._name: str = self._alias or name

    def __get__(
        self, instance: _T | Literal["raise"] | None, owner: type[_T] | None = None, /
    ) -> Any:
        if instance is None:
            # NOTE: Branch for `cls._name`
            # We can check that to see if an instance of `type(self)` for
            # https://narwhals-dev.github.io/narwhals/api-completeness/expr/
            return self
        # NOTE: Prefer not exposing the actual class we're defining in
        # `_implementation` may not be available everywhere
        who = getattr(instance, "_implementation", self._name_owner)
        _raise_not_implemented_error(self._name, who)
        return None  # pragma: no cover

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        # NOTE: Purely to duck-type as assignable to **any** instance method
        # Wouldn't be reachable through *regular* attribute access
        return self.__get__("raise")

    @classmethod
    def deprecated(cls, message: LiteralString, /) -> Self:
        """Alt constructor, wraps with `@deprecated`.

        Arguments:
            message: **Static-only** deprecation message, emitted in an IDE.

        [descriptor]: https://docs.python.org/3/howto/descriptor.html
        """
        obj = cls()
        return deprecated(message)(obj)


def _raise_not_implemented_error(what: str, who: str, /) -> NotImplementedError:
    msg = (
        f"{what!r} is not implemented for: {who!r}.\n\n"
        "If you would like to see this functionality in `narwhals`, "
        "please open an issue at: https://github.com/narwhals-dev/narwhals/issues"
    )
    raise NotImplementedError(msg)


class requires:  # noqa: N801
    """Method decorator for raising under certain constraints.

    Attributes:
        _min_version: Minimum backend version.
        _hint: Optional suggested alternative.

    Examples:
        >>> from narwhals._utils import requires, Implementation
        >>> class SomeBackend:
        ...     _implementation = Implementation.PYARROW
        ...     _backend_version = 20, 0, 0
        ...
        ...     @requires.backend_version((9000, 0, 0))
        ...     def really_complex_feature(self) -> str:
        ...         return "hello"
        >>> backend = SomeBackend()
        >>> backend.really_complex_feature()
        Traceback (most recent call last):
            ...
        NotImplementedError: `really_complex_feature` is only available in 'pyarrow>=9000.0.0', found version '20.0.0'.
    """

    _min_version: tuple[int, ...]
    _hint: str
    _wrapped_name: str
    """(Unqualified) decorated method name.

    When used in a namespace accessor, it will be prefixed by the property name.
    """

    @classmethod
    def backend_version(cls, minimum: tuple[int, ...], /, hint: str = "") -> Self:
        """Method decorator for raising below a minimum `_backend_version`.

        Arguments:
            minimum: Minimum backend version.
            hint: Optional suggested alternative.
        """
        obj = cls.__new__(cls)
        obj._min_version = minimum
        obj._hint = hint
        return obj

    @staticmethod
    def _unparse_version(backend_version: tuple[int, ...], /) -> str:
        return ".".join(f"{d}" for d in backend_version)

    def _qualify_accessor_name(self, prefix: Accessor, /) -> None:
        # NOTE: Should only need to do this once per class (the first time the method is called)
        if "." not in self._wrapped_name:
            self._wrapped_name = f"{prefix}.{self._wrapped_name}"

    def _unwrap_context(self, instance: _IntoContext, /) -> tuple[tuple[int, ...], str]:
        if _is_namespace_accessor(instance):
            self._qualify_accessor_name(instance._accessor)
            compliant = instance.compliant
        else:
            compliant = instance
        return compliant._backend_version, str(compliant._implementation)

    def _ensure_version(self, instance: _IntoContext, /) -> None:
        version, backend = self._unwrap_context(instance)
        if version >= self._min_version:
            return
        minimum = self._unparse_version(self._min_version)
        found = self._unparse_version(version)
        msg = f"`{self._wrapped_name}` is only available in '{backend}>={minimum}', found version {found!r}."
        if self._hint:
            msg = f"{msg}\n{self._hint}"
        raise NotImplementedError(msg)

    def __call__(
        self, fn: _Method[_IntoContextT, P, R], /
    ) -> _Method[_IntoContextT, P, R]:
        self._wrapped_name = fn.__name__

        @wraps(fn)
        def wrapper(instance: _IntoContextT, *args: P.args, **kwds: P.kwargs) -> R:
            self._ensure_version(instance)
            return fn(instance, *args, **kwds)

        # NOTE: Only getting a complaint from `mypy`
        return wrapper  # type: ignore[return-value]


def convert_str_slice_to_int_slice(
    str_slice: _SliceName, columns: Sequence[str]
) -> tuple[int | None, int | None, Any]:
    start = columns.index(str_slice.start) if str_slice.start is not None else None
    stop = columns.index(str_slice.stop) + 1 if str_slice.stop is not None else None
    step = str_slice.step
    return (start, stop, step)


def inherit_doc(
    tp_parent: Callable[P, R1], /
) -> Callable[[_Constructor[_T, P, R2]], _Constructor[_T, P, R2]]:
    """Steal the class-level docstring from parent and attach to child `__init__`.

    Returns:
        Decorated constructor.

    Notes:
        - Passes static typing (mostly)
        - Passes at runtime
    """

    def decorate(init_child: _Constructor[_T, P, R2], /) -> _Constructor[_T, P, R2]:
        if init_child.__name__ == "__init__" and issubclass(type(tp_parent), type):
            init_child.__doc__ = getdoc(tp_parent)
            return init_child
        msg = (  # pragma: no cover
            f"`@{inherit_doc.__name__}` is only allowed to decorate an `__init__` with a class-level doc.\n"
            f"Method: {init_child.__qualname__!r}\n"
            f"Parent: {tp_parent!r}"
        )
        raise TypeError(msg)  # pragma: no cover

    return decorate


def qualified_type_name(obj: object | type[Any], /) -> str:
    tp = obj if isinstance(obj, type) else type(obj)
    module = tp.__module__ if tp.__module__ != "builtins" else ""
    return f"{module}.{tp.__name__}".lstrip(".")


def ensure_type(obj: Any, /, *valid_types: type[Any], param_name: str = "") -> None:
    """Validate that an object is an instance of one or more specified types.

    Parameters:
        obj: The object to validate.
        *valid_types: One or more valid types that `obj` is expected to match.
        param_name: The name of the parameter being validated.
            Used to improve error message clarity.

    Raises:
        TypeError: If `obj` is not an instance of any of the provided `valid_types`.

    Examples:
        >>> from narwhals._utils import ensure_type
        >>> ensure_type(42, int, float)
        >>> ensure_type("hello", str)

        >>> ensure_type("hello", int, param_name="test")
        Traceback (most recent call last):
            ...
        TypeError: Expected 'int', got: 'str'
            test='hello'
                 ^^^^^^^
        >>> import polars as pl
        >>> import pandas as pd
        >>> df = pl.DataFrame([[1], [2], [3], [4], [5]], schema=[*"abcde"])
        >>> ensure_type(df, pd.DataFrame, param_name="df")
        Traceback (most recent call last):
            ...
        TypeError: Expected 'pandas.core.frame.DataFrame', got: 'polars.dataframe.frame.DataFrame'
            df=polars.dataframe.frame.DataFrame(...)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    """
    if not isinstance(obj, valid_types):  # pragma: no cover
        tp_names = " | ".join(qualified_type_name(tp) for tp in valid_types)
        msg = f"Expected {tp_names!r}, got: {qualified_type_name(obj)!r}"
        if param_name:
            left_pad = " " * 4
            val = repr(obj)
            if len(val) > 40:  # truncate long reprs
                val = f"{qualified_type_name(obj)}(...)"
            assign = f"{left_pad}{param_name}="
            underline = (" " * len(assign)) + ("^" * len(val))
            msg = f"{msg}\n{assign}{val}\n{underline}"
        raise TypeError(msg)


class _DeferredIterable(Generic[_T]):
    """Store a callable producing an iterable to defer collection until we need it."""

    def __init__(self, into_iter: Callable[[], Iterable[_T]], /) -> None:
        self._into_iter: Callable[[], Iterable[_T]] = into_iter

    def __iter__(self) -> Iterator[_T]:
        yield from self._into_iter()

    def to_tuple(self) -> tuple[_T, ...]:
        # Collect and return as a `tuple`.
        it = self._into_iter()
        return it if isinstance(it, tuple) else tuple(it)


@lru_cache(maxsize=64)
def deep_attrgetter(attr: str, *nested: str) -> attrgetter[Any]:
    name = ".".join((attr, *nested)) if nested else attr
    return attrgetter(name)


def deep_getattr(obj: Any, name_1: str, *nested: str) -> Any:
    """Perform a nested attribute lookup on `obj`."""
    return deep_attrgetter(name_1, *nested)(obj)


class Compliant(
    _StoresNative[NativeT_co], _StoresImplementation, Protocol[NativeT_co]
): ...


class Narwhals(Protocol[NativeT_co]):
    """Minimal *Narwhals-level* protocol.

    Provides access to a compliant object:

        obj: Narwhals[NativeT_co]]
        compliant: Compliant[NativeT_co] = obj._compliant

    Which itself exposes:

        implementation: Implementation = compliant.implementation
        native: NativeT_co = compliant.native

    This interface is used for revealing which `Implementation` member is associated with **either**:
    - One or more [nominal] native type(s)
    - One or more [structural] type(s)
      - where the true native type(s) are [assignable to] *at least* one of them

    These relationships are defined in the `@overload`s of `_Implementation.__get__(...)`.

    [nominal]: https://typing.python.org/en/latest/spec/glossary.html#term-nominal
    [structural]: https://typing.python.org/en/latest/spec/glossary.html#term-structural
    [assignable to]: https://typing.python.org/en/latest/spec/glossary.html#term-assignable
    """

    @property
    def _compliant(self) -> Compliant[NativeT_co]: ...


class _Implementation:
    """Descriptor for matching an opaque `Implementation` on a generic class.

    Based on [pyright comment](https://github.com/microsoft/pyright/issues/3071#issuecomment-1043978070)
    """

    def __set_name__(self, owner: type[Any], name: str) -> None:
        self.__name__: str = name

    @overload
    def __get__(self, instance: Narwhals[NativePolars], owner: Any) -> _PolarsImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativePandas], owner: Any) -> _PandasImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativeModin], owner: Any) -> _ModinImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativeCuDF], owner: Any) -> _CuDFImpl: ...
    @overload
    def __get__(
        self, instance: Narwhals[NativePandasLike], owner: Any
    ) -> _PandasLikeImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativeArrow], owner: Any) -> _ArrowImpl: ...
    @overload
    def __get__(
        self, instance: Narwhals[NativePolars | NativeArrow | NativePandas], owner: Any
    ) -> _PolarsImpl | _PandasImpl | _ArrowImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativeDuckDB], owner: Any) -> _DuckDBImpl: ...
    @overload
    def __get__(
        self, instance: Narwhals[NativeSQLFrame], owner: Any
    ) -> _SQLFrameImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativeDask], owner: Any) -> _DaskImpl: ...
    @overload
    def __get__(self, instance: Narwhals[NativeIbis], owner: Any) -> _IbisImpl: ...
    @overload
    def __get__(
        self, instance: Narwhals[NativePySpark | NativePySparkConnect], owner: Any
    ) -> _PySparkImpl | _PySparkConnectImpl: ...
    # NOTE: https://docs.python.org/3/howto/descriptor.html#invocation-from-a-class
    @overload
    def __get__(self, instance: None, owner: type[Narwhals[Any]]) -> Self: ...
    @overload
    def __get__(
        self, instance: DataFrame[Any] | Series[Any], owner: Any
    ) -> _EagerAllowedImpl: ...
    @overload
    def __get__(self, instance: LazyFrame[Any], owner: Any) -> _LazyAllowedImpl: ...
    def __get__(self, instance: Narwhals[Any] | None, owner: Any) -> Any:
        return self if instance is None else instance._compliant._implementation


def to_pyarrow_table(tbl: pa.Table | pa.RecordBatchReader) -> pa.Table:
    import pyarrow as pa  # ignore-banned-import

    if isinstance(tbl, pa.RecordBatchReader):  # pragma: no cover
        return pa.Table.from_batches(tbl)
    return tbl


def normalize_path(source: FileSource, /) -> str:
    if isinstance(source, str):
        return source
    from pathlib import Path

    return str(Path(source))


def extend_bool(
    value: bool | Iterable[bool],  # noqa: FBT001
    n_match: int,
) -> Sequence[bool]:
    """Ensure the given bool or sequence of bools is the correct length.

    Stolen from https://github.com/pola-rs/polars/blob/b8bfb07a4a37a8d449d6d1841e345817431142df/py-polars/polars/_utils/various.py#L580-L594
    """
    return (value,) * n_match if isinstance(value, bool) else tuple(value)
