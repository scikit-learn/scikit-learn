from __future__ import annotations

import datetime as dt
from decimal import Decimal
from functools import wraps
from typing import TYPE_CHECKING, Any, Callable, Literal, TypeVar, overload

from narwhals._constants import EPOCH, MS_PER_SECOND
from narwhals._native import (
    is_native_arrow,
    is_native_pandas_like,
    is_native_polars,
    is_native_spark_like,
)
from narwhals._utils import Implementation, Version, has_native_namespace
from narwhals.dependencies import (
    get_dask_expr,
    get_numpy,
    get_pandas,
    is_cupy_scalar,
    is_dask_dataframe,
    is_duckdb_relation,
    is_ibis_table,
    is_numpy_scalar,
    is_pandas_like_dataframe,
    is_polars_lazyframe,
    is_polars_series,
    is_pyarrow_scalar,
    is_pyarrow_table,
)

if TYPE_CHECKING:
    from narwhals.dataframe import DataFrame, LazyFrame
    from narwhals.series import Series
    from narwhals.typing import (
        DataFrameT,
        Frame,
        IntoDataFrameT,
        IntoFrame,
        IntoLazyFrameT,
        IntoSeries,
        IntoSeriesT,
        LazyFrameT,
        SeriesT,
    )

T = TypeVar("T")

NON_TEMPORAL_SCALAR_TYPES = (bool, bytes, str, int, float, complex, Decimal)
TEMPORAL_SCALAR_TYPES = (dt.date, dt.timedelta, dt.time)


@overload
def to_native(
    narwhals_object: DataFrame[IntoDataFrameT], *, pass_through: Literal[False] = ...
) -> IntoDataFrameT: ...
@overload
def to_native(
    narwhals_object: LazyFrame[IntoLazyFrameT], *, pass_through: Literal[False] = ...
) -> IntoLazyFrameT: ...
@overload
def to_native(
    narwhals_object: Series[IntoSeriesT], *, pass_through: Literal[False] = ...
) -> IntoSeriesT: ...
@overload
def to_native(narwhals_object: Any, *, pass_through: bool) -> Any: ...


def to_native(
    narwhals_object: DataFrame[IntoDataFrameT]
    | LazyFrame[IntoLazyFrameT]
    | Series[IntoSeriesT],
    *,
    pass_through: bool = False,
) -> IntoDataFrameT | IntoLazyFrameT | IntoSeriesT | Any:
    """Convert Narwhals object to native one.

    Arguments:
        narwhals_object: Narwhals object.
        pass_through: Determine what happens if `narwhals_object` isn't a Narwhals class

            - `False` (default): raise an error
            - `True`: pass object through as-is

    Returns:
        Object of class that user started with.
    """
    from narwhals.dataframe import BaseFrame
    from narwhals.series import Series

    if isinstance(narwhals_object, BaseFrame):
        return narwhals_object._compliant_frame._native_frame
    if isinstance(narwhals_object, Series):
        return narwhals_object._compliant_series.native

    if not pass_through:
        msg = f"Expected Narwhals object, got {type(narwhals_object)}."
        raise TypeError(msg)
    return narwhals_object


@overload
def from_native(native_object: SeriesT, **kwds: Any) -> SeriesT: ...


@overload
def from_native(native_object: DataFrameT, **kwds: Any) -> DataFrameT: ...


@overload
def from_native(native_object: LazyFrameT, **kwds: Any) -> LazyFrameT: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    pass_through: Literal[True],
    eager_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoLazyFrameT | IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | LazyFrame[IntoLazyFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoLazyFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> LazyFrame[IntoLazyFrameT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoFrame | IntoSeries,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[Any] | LazyFrame[Any] | Series[Any]: ...


@overload
def from_native(
    native_object: IntoSeriesT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


# All params passed in as variables
@overload
def from_native(
    native_object: Any,
    *,
    pass_through: bool,
    eager_only: bool,
    series_only: bool,
    allow_series: bool | None,
) -> Any: ...


def from_native(  # noqa: D417
    native_object: IntoLazyFrameT
    | IntoDataFrameT
    | IntoSeriesT
    | IntoFrame
    | IntoSeries
    | T,
    *,
    pass_through: bool = False,
    eager_only: bool = False,
    series_only: bool = False,
    allow_series: bool | None = None,
    **kwds: Any,
) -> LazyFrame[IntoLazyFrameT] | DataFrame[IntoDataFrameT] | Series[IntoSeriesT] | T:
    """Convert `native_object` to Narwhals Dataframe, Lazyframe, or Series.

    Arguments:
        native_object: Raw object from user.
            Depending on the other arguments, input object can be

            - a Dataframe / Lazyframe / Series supported by Narwhals (pandas, Polars, PyArrow, ...)
            - an object which implements `__narwhals_dataframe__`, `__narwhals_lazyframe__`,
              or `__narwhals_series__`
        pass_through: Determine what happens if the object can't be converted to Narwhals

            - `False` (default): raise an error
            - `True`: pass object through as-is
        eager_only: Whether to only allow eager objects

            - `False` (default): don't require `native_object` to be eager
            - `True`: only convert to Narwhals if `native_object` is eager
        series_only: Whether to only allow Series

            - `False` (default): don't require `native_object` to be a Series
            - `True`: only convert to Narwhals if `native_object` is a Series
        allow_series: Whether to allow Series (default is only Dataframe / Lazyframe)

            - `False` or `None` (default): don't convert to Narwhals if `native_object` is a Series
            - `True`: allow `native_object` to be a Series

    Returns:
        DataFrame, LazyFrame, Series, or original object, depending
            on which combination of parameters was passed.
    """
    if kwds:
        msg = f"from_native() got an unexpected keyword argument {next(iter(kwds))!r}"
        raise TypeError(msg)

    return _from_native_impl(  # type: ignore[no-any-return]
        native_object,
        pass_through=pass_through,
        eager_only=eager_only,
        eager_or_interchange_only=False,
        series_only=series_only,
        allow_series=allow_series,
        version=Version.MAIN,
    )


def _from_native_impl(  # noqa: C901, PLR0911, PLR0912, PLR0915
    native_object: Any,
    *,
    pass_through: bool = False,
    eager_only: bool = False,
    # Interchange-level was removed after v1
    eager_or_interchange_only: bool = False,
    series_only: bool = False,
    allow_series: bool | None = None,
    version: Version,
) -> Any:
    from narwhals._interchange.dataframe import supports_dataframe_interchange
    from narwhals._utils import (
        is_compliant_dataframe,
        is_compliant_lazyframe,
        is_compliant_series,
    )
    from narwhals.dataframe import DataFrame, LazyFrame
    from narwhals.series import Series

    # Early returns
    if isinstance(native_object, (DataFrame, LazyFrame)) and not series_only:
        return native_object
    if isinstance(native_object, Series) and (series_only or allow_series):
        return native_object

    if series_only:
        if allow_series is False:
            msg = "Invalid parameter combination: `series_only=True` and `allow_series=False`"
            raise ValueError(msg)
        allow_series = True
    if eager_only and eager_or_interchange_only:
        msg = "Invalid parameter combination: `eager_only=True` and `eager_or_interchange_only=True`"
        raise ValueError(msg)

    # Extensions
    if is_compliant_dataframe(native_object):
        if series_only:
            if not pass_through:
                msg = "Cannot only use `series_only` with dataframe"
                raise TypeError(msg)
            return native_object
        return version.dataframe(
            native_object.__narwhals_dataframe__()._with_version(version), level="full"
        )
    if is_compliant_lazyframe(native_object):
        if series_only:
            if not pass_through:
                msg = "Cannot only use `series_only` with lazyframe"
                raise TypeError(msg)
            return native_object
        if eager_only or eager_or_interchange_only:
            if not pass_through:
                msg = "Cannot only use `eager_only` or `eager_or_interchange_only` with lazyframe"
                raise TypeError(msg)
            return native_object
        return version.lazyframe(
            native_object.__narwhals_lazyframe__()._with_version(version), level="full"
        )
    if is_compliant_series(native_object):
        if not allow_series:
            if not pass_through:
                msg = "Please set `allow_series=True` or `series_only=True`"
                raise TypeError(msg)
            return native_object
        return version.series(
            native_object.__narwhals_series__()._with_version(version), level="full"
        )

    # Polars
    if is_native_polars(native_object):
        if series_only and not is_polars_series(native_object):
            if not pass_through:
                msg = f"Cannot only use `series_only` with {type(native_object).__qualname__}"
                raise TypeError(msg)
            return native_object
        if (eager_only or eager_or_interchange_only) and is_polars_lazyframe(
            native_object
        ):
            if not pass_through:
                msg = "Cannot only use `eager_only` or `eager_or_interchange_only` with polars.LazyFrame"
                raise TypeError(msg)
            return native_object
        if (not allow_series) and is_polars_series(native_object):
            if not pass_through:
                msg = "Please set `allow_series=True` or `series_only=True`"
                raise TypeError(msg)
            return native_object
        return (
            version.namespace.from_native_object(native_object)
            .compliant.from_native(native_object)
            .to_narwhals()
        )

    # PandasLike
    if is_native_pandas_like(native_object):
        if is_pandas_like_dataframe(native_object):
            if series_only:
                if not pass_through:
                    msg = f"Cannot only use `series_only` with {type(native_object).__qualname__}"
                    raise TypeError(msg)
                return native_object
        elif not allow_series:
            if not pass_through:
                msg = "Please set `allow_series=True` or `series_only=True`"
                raise TypeError(msg)
            return native_object
        return (
            version.namespace.from_native_object(native_object)
            .compliant.from_native(native_object)
            .to_narwhals()
        )

    # PyArrow
    if is_native_arrow(native_object):
        if is_pyarrow_table(native_object):
            if series_only:
                if not pass_through:
                    msg = f"Cannot only use `series_only` with {type(native_object).__qualname__}"
                    raise TypeError(msg)
                return native_object
        elif not allow_series:
            if not pass_through:
                msg = "Please set `allow_series=True` or `series_only=True`"
                raise TypeError(msg)
            return native_object
        return (
            version.namespace.from_native_object(native_object)
            .compliant.from_native(native_object)
            .to_narwhals()
        )

    # Dask
    if is_dask_dataframe(native_object):
        if series_only:
            if not pass_through:
                msg = "Cannot only use `series_only` with dask DataFrame"
                raise TypeError(msg)
            return native_object
        if eager_only or eager_or_interchange_only:
            if not pass_through:
                msg = "Cannot only use `eager_only` or `eager_or_interchange_only` with dask DataFrame"
                raise TypeError(msg)
            return native_object
        if (
            Implementation.DASK._backend_version() <= (2024, 12, 1)
            and get_dask_expr() is None
        ):  # pragma: no cover
            msg = "Please install dask-expr"
            raise ImportError(msg)
        return (
            version.namespace.from_backend(Implementation.DASK)
            .compliant.from_native(native_object)
            .to_narwhals()
        )

    # DuckDB
    if is_duckdb_relation(native_object):
        if eager_only or series_only:  # pragma: no cover
            if not pass_through:
                msg = "Cannot only use `series_only=True` or `eager_only=False` with DuckDBPyRelation"
                raise TypeError(msg)
            return native_object
        return (
            version.namespace.from_native_object(native_object)
            .compliant.from_native(native_object)
            .to_narwhals()
        )

    # Ibis
    if is_ibis_table(native_object):
        if eager_only or series_only:  # pragma: no cover
            if not pass_through:
                msg = "Cannot only use `series_only=True` or `eager_only=False` with ibis.Table"
                raise TypeError(msg)
            return native_object
        return (
            version.namespace.from_native_object(native_object)
            .compliant.from_native(native_object)
            .to_narwhals()
        )

    # PySpark
    if is_native_spark_like(native_object):  # pragma: no cover
        ns_spark = version.namespace.from_native_object(native_object)
        if series_only or eager_only or eager_or_interchange_only:
            if not pass_through:
                msg = (
                    "Cannot only use `series_only`, `eager_only` or `eager_or_interchange_only` "
                    f"with {ns_spark.implementation} DataFrame"
                )
                raise TypeError(msg)
            return native_object
        return ns_spark.compliant.from_native(native_object).to_narwhals()

    # Interchange protocol
    if supports_dataframe_interchange(native_object):
        from narwhals._interchange.dataframe import InterchangeFrame

        if eager_only or series_only:
            if not pass_through:
                msg = (
                    "Cannot only use `series_only=True` or `eager_only=False` "
                    "with object which only implements __dataframe__"
                )
                raise TypeError(msg)
            return native_object
        if version is not Version.V1:
            if pass_through:
                return native_object
            msg = (
                "The Dataframe Interchange Protocol is no longer supported in the main `narwhals` namespace.\n\n"
                "You may want to:\n"
                " - Use `narwhals.stable.v1`, where it is still supported.\n"
                "    - See https://narwhals-dev.github.io/narwhals/backcompat\n"
                " - Use `pass_through=True` to pass the object through without raising."
            )
            raise TypeError(msg)
        return Version.V1.dataframe(InterchangeFrame(native_object), level="interchange")

    if not pass_through:
        msg = f"Expected pandas-like dataframe, Polars dataframe, or Polars lazyframe, got: {type(native_object)}"
        raise TypeError(msg)
    return native_object


def get_native_namespace(*obj: Frame | Series[Any] | IntoFrame | IntoSeries) -> Any:
    """Get native namespace from object.

    Arguments:
        obj: Dataframe, Lazyframe, or Series. Multiple objects can be
            passed positionally, in which case they must all have the
            same native namespace (else an error is raised).

    Returns:
        Native module.

    Examples:
        >>> import polars as pl
        >>> import pandas as pd
        >>> import narwhals as nw
        >>> df = nw.from_native(pd.DataFrame({"a": [1, 2, 3]}))
        >>> nw.get_native_namespace(df)
        <module 'pandas'...>
        >>> df = nw.from_native(pl.DataFrame({"a": [1, 2, 3]}))
        >>> nw.get_native_namespace(df)
        <module 'polars'...>
    """
    if not obj:
        msg = "At least one object must be passed to `get_native_namespace`."
        raise ValueError(msg)
    result = {_get_native_namespace_single_obj(x) for x in obj}
    if len(result) != 1:
        msg = f"Found objects with different native namespaces: {result}."
        raise ValueError(msg)
    return result.pop()


def _get_native_namespace_single_obj(
    obj: Frame | Series[Any] | IntoFrame | IntoSeries,
) -> Any:
    if has_native_namespace(obj):
        return obj.__native_namespace__()
    return Version.MAIN.namespace.from_native_object(
        obj
    ).implementation.to_native_namespace()


def narwhalify(
    func: Callable[..., Any] | None = None,
    *,
    pass_through: bool = True,
    eager_only: bool = False,
    series_only: bool = False,
    allow_series: bool | None = True,
) -> Callable[..., Any]:
    """Decorate function so it becomes dataframe-agnostic.

    This will try to convert any dataframe/series-like object into the Narwhals
    respective DataFrame/Series, while leaving the other parameters as they are.
    Similarly, if the output of the function is a Narwhals DataFrame or Series, it will be
    converted back to the original dataframe/series type, while if the output is another
    type it will be left as is.
    By setting `pass_through=False`, then every input and every output will be required to be a
    dataframe/series-like object.

    Arguments:
        func: Function to wrap in a `from_native`-`to_native` block.
        pass_through: Determine what happens if the object can't be converted to Narwhals

            - `False`: raise an error
            - `True` (default): pass object through as-is
        eager_only: Whether to only allow eager objects

            - `False` (default): don't require `native_object` to be eager
            - `True`: only convert to Narwhals if `native_object` is eager
        series_only: Whether to only allow Series

            - `False` (default): don't require `native_object` to be a Series
            - `True`: only convert to Narwhals if `native_object` is a Series
        allow_series: Whether to allow Series (default is only Dataframe / Lazyframe)

            - `False` or `None`: don't convert to Narwhals if `native_object` is a Series
            - `True` (default): allow `native_object` to be a Series

    Returns:
        Decorated function.

    Examples:
        Instead of writing

        >>> import narwhals as nw
        >>> def agnostic_group_by_sum(df):
        ...     df = nw.from_native(df, pass_through=True)
        ...     df = df.group_by("a").agg(nw.col("b").sum())
        ...     return nw.to_native(df)

        you can just write

        >>> @nw.narwhalify
        ... def agnostic_group_by_sum(df):
        ...     return df.group_by("a").agg(nw.col("b").sum())
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            args = [
                from_native(
                    arg,
                    pass_through=pass_through,
                    eager_only=eager_only,
                    series_only=series_only,
                    allow_series=allow_series,
                )
                for arg in args
            ]  # type: ignore[assignment]

            kwargs = {
                name: from_native(
                    value,
                    pass_through=pass_through,
                    eager_only=eager_only,
                    series_only=series_only,
                    allow_series=allow_series,
                )
                for name, value in kwargs.items()
            }

            backends = {
                b()
                for v in (*args, *kwargs.values())
                if (b := getattr(v, "__native_namespace__", None))
            }

            if len(backends) > 1:
                msg = "Found multiple backends. Make sure that all dataframe/series inputs come from the same backend."
                raise ValueError(msg)

            result = func(*args, **kwargs)

            return to_native(result, pass_through=pass_through)

        return wrapper

    if func is None:
        return decorator
    # If func is not None, it means the decorator is used without arguments
    return decorator(func)


def to_py_scalar(scalar_like: Any) -> Any:
    """If a scalar is not Python native, converts it to Python native.

    Arguments:
        scalar_like: Scalar-like value.

    Raises:
        ValueError: If the object is not convertible to a scalar.

    Examples:
        >>> import narwhals as nw
        >>> import pandas as pd
        >>> df = nw.from_native(pd.DataFrame({"a": [1, 2, 3]}))
        >>> nw.to_py_scalar(df["a"].item(0))
        1
        >>> import pyarrow as pa
        >>> df = nw.from_native(pa.table({"a": [1, 2, 3]}))
        >>> nw.to_py_scalar(df["a"].item(0))
        1
        >>> nw.to_py_scalar(1)
        1
    """
    scalar: Any
    pd = get_pandas()
    if scalar_like is None or isinstance(scalar_like, NON_TEMPORAL_SCALAR_TYPES):
        scalar = scalar_like
    elif (
        (np := get_numpy())
        and isinstance(scalar_like, np.datetime64)
        and scalar_like.dtype == "datetime64[ns]"
    ):
        ms = scalar_like.item() // MS_PER_SECOND
        scalar = EPOCH + dt.timedelta(microseconds=ms)
    elif is_numpy_scalar(scalar_like) or is_cupy_scalar(scalar_like):
        scalar = scalar_like.item()
    elif pd and isinstance(scalar_like, pd.Timestamp):
        scalar = scalar_like.to_pydatetime()
    elif pd and isinstance(scalar_like, pd.Timedelta):
        scalar = scalar_like.to_pytimedelta()
    # pd.Timestamp and pd.Timedelta subclass datetime and timedelta,
    # so we need to check this separately
    elif isinstance(scalar_like, TEMPORAL_SCALAR_TYPES):
        scalar = scalar_like
    elif _is_pandas_na(scalar_like):
        scalar = None
    elif is_pyarrow_scalar(scalar_like):
        scalar = scalar_like.as_py()
    else:
        msg = (
            f"Expected object convertible to a scalar, found {type(scalar_like)}.\n"
            f"{scalar_like!r}"
        )
        raise ValueError(msg)
    return scalar


def _is_pandas_na(obj: Any) -> bool:
    return bool((pd := get_pandas()) and pd.api.types.is_scalar(obj) and pd.isna(obj))


__all__ = ["get_native_namespace", "narwhalify", "to_native", "to_py_scalar"]
