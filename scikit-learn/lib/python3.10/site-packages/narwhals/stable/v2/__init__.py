from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING, Any, Callable, Final, Literal, cast, overload

import narwhals as nw
from narwhals import exceptions, functions as nw_f
from narwhals._typing_compat import TypeVar, assert_never
from narwhals._utils import (
    Implementation,
    Version,
    generate_temporary_column_name,
    inherit_doc,
    is_ordered_categorical,
    maybe_align_index,
    maybe_convert_dtypes,
    maybe_get_index,
    maybe_reset_index,
    maybe_set_index,
    not_implemented,
)
from narwhals.dataframe import DataFrame as NwDataFrame, LazyFrame as NwLazyFrame
from narwhals.dtypes import (
    Array,
    Binary,
    Boolean,
    Categorical,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Field,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    List,
    Object,
    String,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    Unknown,
)
from narwhals.expr import Expr as NwExpr
from narwhals.functions import _new_series_impl, concat, show_versions
from narwhals.schema import Schema as NwSchema
from narwhals.series import Series as NwSeries
from narwhals.stable.v2 import dependencies, dtypes, selectors
from narwhals.translate import _from_native_impl, get_native_namespace, to_py_scalar
from narwhals.typing import IntoDataFrameT, IntoLazyFrameT

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

    from typing_extensions import ParamSpec, Self

    from narwhals._translate import IntoArrowTable
    from narwhals._typing import (
        Arrow,
        Backend,
        EagerAllowed,
        IntoBackend,
        LazyAllowed,
        Pandas,
        Polars,
    )
    from narwhals.dataframe import MultiColSelector, MultiIndexSelector
    from narwhals.dtypes import DType
    from narwhals.typing import (
        IntoDType,
        IntoExpr,
        IntoFrame,
        IntoSchema,
        IntoSeries,
        NonNestedLiteral,
        SingleColSelector,
        SingleIndexSelector,
        _1DArray,
        _2DArray,
    )

    DataFrameT = TypeVar("DataFrameT", bound="DataFrame[Any]")
    LazyFrameT = TypeVar("LazyFrameT", bound="LazyFrame[Any]")
    SeriesT = TypeVar("SeriesT", bound="Series[Any]")
    T = TypeVar("T", default=Any)
    P = ParamSpec("P")
    R = TypeVar("R")

IntoSeriesT = TypeVar("IntoSeriesT", bound="IntoSeries", default=Any)


class DataFrame(NwDataFrame[IntoDataFrameT]):
    _version = Version.V2

    @inherit_doc(NwDataFrame)
    def __init__(self, df: Any, *, level: Literal["full", "lazy", "interchange"]) -> None:
        assert df._version is Version.V2  # noqa: S101
        super().__init__(df, level=level)

    # We need to override any method which don't return Self so that type
    # annotations are correct.

    @classmethod
    def from_arrow(
        cls, native_frame: IntoArrowTable, *, backend: IntoBackend[EagerAllowed]
    ) -> DataFrame[Any]:
        result = super().from_arrow(native_frame, backend=backend)
        return cast("DataFrame[Any]", result)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        schema: Mapping[str, DType] | Schema | None = None,
        *,
        backend: IntoBackend[EagerAllowed] | None = None,
    ) -> DataFrame[Any]:
        result = super().from_dict(data, schema, backend=backend)
        return cast("DataFrame[Any]", result)

    @classmethod
    def from_dicts(
        cls,
        data: Sequence[Mapping[str, Any]],
        schema: IntoSchema | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> DataFrame[Any]:
        result = super().from_dicts(data, schema, backend=backend)
        return cast("DataFrame[Any]", result)

    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        schema: Mapping[str, DType] | Schema | Sequence[str] | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> DataFrame[Any]:
        result = super().from_numpy(data, schema, backend=backend)
        return cast("DataFrame[Any]", result)

    @property
    def _series(self) -> type[Series[Any]]:
        return cast("type[Series[Any]]", Series)

    @property
    def _lazyframe(self) -> type[LazyFrame[Any]]:
        return cast("type[LazyFrame[Any]]", LazyFrame)

    @overload
    def __getitem__(self, item: tuple[SingleIndexSelector, SingleColSelector]) -> Any: ...

    @overload
    def __getitem__(  # type: ignore[overload-overlap]
        self, item: str | tuple[MultiIndexSelector, SingleColSelector]
    ) -> Series[Any]: ...

    @overload
    def __getitem__(
        self,
        item: (
            SingleIndexSelector
            | MultiIndexSelector
            | MultiColSelector
            | tuple[SingleIndexSelector, MultiColSelector]
            | tuple[MultiIndexSelector, MultiColSelector]
        ),
    ) -> Self: ...
    def __getitem__(
        self,
        item: (
            SingleIndexSelector
            | SingleColSelector
            | MultiColSelector
            | MultiIndexSelector
            | tuple[SingleIndexSelector, SingleColSelector]
            | tuple[SingleIndexSelector, MultiColSelector]
            | tuple[MultiIndexSelector, SingleColSelector]
            | tuple[MultiIndexSelector, MultiColSelector]
        ),
    ) -> Series[Any] | Self | Any:
        return super().__getitem__(item)

    def get_column(self, name: str) -> Series:
        # Type checkers complain that `nw.Series` is not assignable to `nw.v2.stable.Series`.
        # However the return type actually is `nw.v2.stable.Series`, check `tests/v2_test.py`.
        return super().get_column(name)  # type: ignore[return-value]

    def lazy(
        self,
        backend: IntoBackend[LazyAllowed] | None = None,
        *,
        session: Any | None = None,
    ) -> LazyFrame[Any]:
        return _stableify(super().lazy(backend=backend, session=session))

    @overload  # type: ignore[override]
    def to_dict(self, *, as_series: Literal[True] = ...) -> dict[str, Series[Any]]: ...
    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...
    @overload
    def to_dict(
        self, *, as_series: bool = True
    ) -> dict[str, Series[Any]] | dict[str, list[Any]]: ...
    def to_dict(  # pyright: ignore[reportIncompatibleMethodOverride]
        self, *, as_series: bool = True
    ) -> dict[str, Series[Any]] | dict[str, list[Any]]:
        # Type checkers complain that `nw.Series` is not assignable to `nw.v2.stable.Series`.
        # However the return type actually is `nw.v2.stable.Series`, check `tests/v2_test.py::test_to_dict_as_series`.
        return super().to_dict(as_series=as_series)  # type: ignore[return-value]

    def is_duplicated(self) -> Series[Any]:
        return _stableify(super().is_duplicated())

    def is_unique(self) -> Series[Any]:
        return _stableify(super().is_unique())


class LazyFrame(NwLazyFrame[IntoLazyFrameT]):
    @inherit_doc(NwLazyFrame)
    def __init__(self, df: Any, *, level: Literal["full", "lazy", "interchange"]) -> None:
        assert df._version is Version.V2  # noqa: S101
        super().__init__(df, level=level)

    @property
    def _dataframe(self) -> type[DataFrame[Any]]:
        return DataFrame

    def collect(
        self, backend: IntoBackend[Polars | Pandas | Arrow] | None = None, **kwargs: Any
    ) -> DataFrame[Any]:
        return _stableify(super().collect(backend=backend, **kwargs))


class Series(NwSeries[IntoSeriesT]):
    _version = Version.V2

    @inherit_doc(NwSeries)
    def __init__(
        self, series: Any, *, level: Literal["full", "lazy", "interchange"]
    ) -> None:
        assert series._version is Version.V2  # noqa: S101
        super().__init__(series, level=level)

    # We need to override any method which don't return Self so that type
    # annotations are correct.

    @property
    def _dataframe(self) -> type[DataFrame[Any]]:
        return DataFrame

    @classmethod
    def from_numpy(
        cls,
        name: str,
        values: _1DArray,
        dtype: IntoDType | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> Series[Any]:
        result = super().from_numpy(name, values, dtype, backend=backend)
        return cast("Series[Any]", result)

    @classmethod
    def from_iterable(
        cls,
        name: str,
        values: Iterable[Any],
        dtype: IntoDType | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> Series[Any]:
        result = super().from_iterable(name, values, dtype, backend=backend)
        return cast("Series[Any]", result)

    def to_frame(self) -> DataFrame[Any]:
        return _stableify(super().to_frame())

    def value_counts(
        self,
        *,
        sort: bool = False,
        parallel: bool = False,
        name: str | None = None,
        normalize: bool = False,
    ) -> DataFrame[Any]:
        return _stableify(
            super().value_counts(
                sort=sort, parallel=parallel, name=name, normalize=normalize
            )
        )

    # Too unstable to consider including here.
    hist: Any = not_implemented()


class Expr(NwExpr): ...


class Schema(NwSchema):
    _version = Version.V2

    @inherit_doc(NwSchema)
    def __init__(
        self, schema: Mapping[str, DType] | Iterable[tuple[str, DType]] | None = None
    ) -> None:
        super().__init__(schema)


@overload
def _stableify(obj: NwDataFrame[IntoDataFrameT]) -> DataFrame[IntoDataFrameT]: ...
@overload
def _stableify(obj: NwLazyFrame[IntoLazyFrameT]) -> LazyFrame[IntoLazyFrameT]: ...
@overload
def _stableify(obj: NwSeries[IntoSeriesT]) -> Series[IntoSeriesT]: ...
@overload
def _stableify(obj: NwExpr) -> Expr: ...


def _stableify(
    obj: NwDataFrame[IntoDataFrameT]
    | NwLazyFrame[IntoLazyFrameT]
    | NwSeries[IntoSeriesT]
    | NwExpr,
) -> DataFrame[IntoDataFrameT] | LazyFrame[IntoLazyFrameT] | Series[IntoSeriesT] | Expr:
    if isinstance(obj, NwDataFrame):
        return DataFrame(obj._compliant_frame._with_version(Version.V2), level=obj._level)
    if isinstance(obj, NwLazyFrame):
        return LazyFrame(obj._compliant_frame._with_version(Version.V2), level=obj._level)
    if isinstance(obj, NwSeries):
        return Series(obj._compliant_series._with_version(Version.V2), level=obj._level)
    if isinstance(obj, NwExpr):
        return Expr(obj._to_compliant_expr, obj._metadata)
    assert_never(obj)


@overload
def from_native(native_object: SeriesT, **kwds: Any) -> SeriesT: ...


@overload
def from_native(native_object: DataFrameT, **kwds: Any) -> DataFrameT: ...


@overload
def from_native(native_object: LazyFrameT, **kwds: Any) -> LazyFrameT: ...


@overload
def from_native(
    native_object: DataFrameT | LazyFrameT, **kwds: Any
) -> DataFrameT | LazyFrameT: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoSeries,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT]: ...


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
    native_object: IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


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
    native_object: IntoDataFrameT
    | IntoLazyFrameT
    | IntoFrame
    | IntoSeriesT
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
    # Early returns
    if isinstance(native_object, (DataFrame, LazyFrame)) and not series_only:
        return native_object
    if isinstance(native_object, Series) and (series_only or allow_series):
        return native_object

    if kwds:
        msg = f"from_native() got an unexpected keyword argument {next(iter(kwds))!r}"
        raise TypeError(msg)

    return _from_native_impl(  # type: ignore[no-any-return]
        native_object,
        pass_through=pass_through,
        eager_only=eager_only,
        series_only=series_only,
        allow_series=allow_series,
        version=Version.V2,
    )


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
    return nw.to_native(narwhals_object, pass_through=pass_through)


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

            if backends.__len__() > 1:
                msg = "Found multiple backends. Make sure that all dataframe/series inputs come from the same backend."
                raise ValueError(msg)

            result = func(*args, **kwargs)

            return to_native(result, pass_through=pass_through)

        return wrapper

    if func is None:
        return decorator
    # If func is not None, it means the decorator is used without arguments
    return decorator(func)


def all() -> Expr:
    """Instantiate an expression representing all columns."""
    return _stableify(nw.all())


def col(*names: str | Iterable[str]) -> Expr:
    """Creates an expression that references one or more columns by their name(s).

    Arguments:
        names: Name(s) of the columns to use.
    """
    return _stableify(nw.col(*names))


def exclude(*names: str | Iterable[str]) -> Expr:
    """Creates an expression that excludes columns by their name(s).

    Arguments:
        names: Name(s) of the columns to exclude.
    """
    return _stableify(nw.exclude(*names))


def nth(*indices: int | Sequence[int]) -> Expr:
    """Creates an expression that references one or more columns by their index(es).

    Notes:
        `nth` is not supported for Polars version<1.0.0. Please use
        [`narwhals.col`][] instead.

    Arguments:
        indices: One or more indices representing the columns to retrieve.
    """
    return _stableify(nw.nth(*indices))


def len() -> Expr:
    """Return the number of rows."""
    return _stableify(nw.len())


def lit(value: NonNestedLiteral, dtype: IntoDType | None = None) -> Expr:
    """Return an expression representing a literal value.

    Arguments:
        value: The value to use as literal.
        dtype: The data type of the literal value. If not provided, the data type will
            be inferred by the native library.
    """
    return _stableify(nw.lit(value, dtype))


def min(*columns: str) -> Expr:
    """Return the minimum value.

    Note:
       Syntactic sugar for ``nw.col(columns).min()``.

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function.
    """
    return _stableify(nw.min(*columns))


def max(*columns: str) -> Expr:
    """Return the maximum value.

    Note:
       Syntactic sugar for ``nw.col(columns).max()``.

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function.
    """
    return _stableify(nw.max(*columns))


def mean(*columns: str) -> Expr:
    """Get the mean value.

    Note:
        Syntactic sugar for ``nw.col(columns).mean()``

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function
    """
    return _stableify(nw.mean(*columns))


def median(*columns: str) -> Expr:
    """Get the median value.

    Notes:
        - Syntactic sugar for ``nw.col(columns).median()``
        - Results might slightly differ across backends due to differences in the
            underlying algorithms used to compute the median.

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function
    """
    return _stableify(nw.median(*columns))


def sum(*columns: str) -> Expr:
    """Sum all values.

    Note:
        Syntactic sugar for ``nw.col(columns).sum()``

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function
    """
    return _stableify(nw.sum(*columns))


def sum_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Sum all values horizontally across columns.

    Warning:
        Unlike Polars, we support horizontal sum over numeric columns only.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
    """
    return _stableify(nw.sum_horizontal(*exprs))


def all_horizontal(*exprs: IntoExpr | Iterable[IntoExpr], ignore_nulls: bool) -> Expr:
    r"""Compute the bitwise AND horizontally across columns.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
        ignore_nulls: Whether to ignore nulls:

            - If `True`, null values are ignored. If there are no elements, the result
              is `True`.
            - If `False`, Kleene logic is followed. Note that this is not allowed for
              pandas with classical NumPy dtypes when null values are present.
    """
    return _stableify(nw.all_horizontal(*exprs, ignore_nulls=ignore_nulls))


def any_horizontal(*exprs: IntoExpr | Iterable[IntoExpr], ignore_nulls: bool) -> Expr:
    r"""Compute the bitwise OR horizontally across columns.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
        ignore_nulls: Whether to ignore nulls:

            - If `True`, null values are ignored. If there are no elements, the result
              is `False`.
            - If `False`, Kleene logic is followed. Note that this is not allowed for
              pandas with classical NumPy dtypes when null values are present.
    """
    return _stableify(nw.any_horizontal(*exprs, ignore_nulls=ignore_nulls))


def mean_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Compute the mean of all values horizontally across columns.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
    """
    return _stableify(nw.mean_horizontal(*exprs))


def min_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Get the minimum value horizontally across columns.

    Notes:
        We support `min_horizontal` over numeric columns only.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
    """
    return _stableify(nw.min_horizontal(*exprs))


def max_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Get the maximum value horizontally across columns.

    Notes:
        We support `max_horizontal` over numeric columns only.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
    """
    return _stableify(nw.max_horizontal(*exprs))


def concat_str(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    separator: str = "",
    ignore_nulls: bool = False,
) -> Expr:
    r"""Horizontally concatenate columns into a single string column.

    Arguments:
        exprs: Columns to concatenate into a single string column. Accepts expression
            input. Strings are parsed as column names, other non-expression inputs are
            parsed as literals. Non-`String` columns are cast to `String`.
        *more_exprs: Additional columns to concatenate into a single string column,
            specified as positional arguments.
        separator: String that will be used to separate the values of each column.
        ignore_nulls: Ignore null values (default is `False`).
            If set to `False`, null values will be propagated and if the row contains any
            null values, the output is null.
    """
    return _stableify(
        nw.concat_str(exprs, *more_exprs, separator=separator, ignore_nulls=ignore_nulls)
    )


def format(f_string: str, *args: IntoExpr) -> Expr:
    """Format expressions as a string.

    Arguments:
        f_string: A string that with placeholders.
        args: Expression(s) that fill the placeholders.
    """
    return _stableify(nw.format(f_string, *args))


def coalesce(exprs: IntoExpr | Iterable[IntoExpr], *more_exprs: IntoExpr) -> Expr:
    """Folds the columns from left to right, keeping the first non-null value.

    Arguments:
        exprs: Columns to coalesce, must be a str, nw.Expr, or nw.Series
            where strings are parsed as column names and both nw.Expr/nw.Series
            are passed through as-is. Scalar values must be wrapped in `nw.lit`.

        *more_exprs: Additional columns to coalesce, specified as positional arguments.

    Raises:
        TypeError: If any of the inputs are not a str, nw.Expr, or nw.Series.
    """
    return _stableify(nw.coalesce(exprs, *more_exprs))


class When(nw_f.When):
    @classmethod
    def from_when(cls, when: nw_f.When) -> When:
        return cls(when._predicate)

    def then(self, value: IntoExpr | NonNestedLiteral | _1DArray) -> Then:
        return Then.from_then(super().then(value))


class Then(nw_f.Then, Expr):
    @classmethod
    def from_then(cls, then: nw_f.Then) -> Then:
        return cls(then._to_compliant_expr, then._metadata)

    def otherwise(self, value: IntoExpr | NonNestedLiteral | _1DArray) -> Expr:
        return _stableify(super().otherwise(value))


def when(*predicates: IntoExpr | Iterable[IntoExpr]) -> When:
    """Start a `when-then-otherwise` expression.

    Expression similar to an `if-else` statement in Python. Always initiated by a
    `pl.when(<condition>).then(<value if condition>)`, and optionally followed by a
    `.otherwise(<value if condition is false>)` can be appended at the end. If not
    appended, and the condition is not `True`, `None` will be returned.

    Info:
        Chaining multiple `.when(<condition>).then(<value>)` statements is currently
        not supported.
        See [Narwhals#668](https://github.com/narwhals-dev/narwhals/issues/668).

    Arguments:
        predicates: Condition(s) that must be met in order to apply the subsequent
            statement. Accepts one or more boolean expressions, which are implicitly
            combined with `&`. String input is parsed as a column name.

    Returns:
        A "when" object, which `.then` can be called on.
    """
    return When.from_when(nw_f.when(*predicates))


def new_series(
    name: str,
    values: Any,
    dtype: IntoDType | None = None,
    *,
    backend: IntoBackend[EagerAllowed],
) -> Series[Any]:
    """Instantiate Narwhals Series from iterable (e.g. list or array).

    Arguments:
        name: Name of resulting Series.
        values: Values of make Series from.
        dtype: (Narwhals) dtype. If not provided, the native library
            may auto-infer it from `values`.
        backend: specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
    """
    return _stableify(_new_series_impl(name, values, dtype, backend=backend))


def from_arrow(
    native_frame: IntoArrowTable, *, backend: IntoBackend[EagerAllowed]
) -> DataFrame[Any]:
    """Construct a DataFrame from an object which supports the PyCapsule Interface.

    Arguments:
        native_frame: Object which implements `__arrow_c_stream__`.
        backend: specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
    """
    return _stableify(nw_f.from_arrow(native_frame, backend=backend))


def from_dict(
    data: Mapping[str, Any],
    schema: Mapping[str, DType] | Schema | None = None,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
) -> DataFrame[Any]:
    """Instantiate DataFrame from dictionary.

    Indexes (if present, for pandas-like backends) are aligned following
    the [left-hand-rule](../concepts/pandas_index.md/).

    Notes:
        For pandas-like dataframes, conversion to schema is applied after dataframe
        creation.

    Arguments:
        data: Dictionary to create DataFrame from.
        schema: The DataFrame schema as Schema or dict of {name: type}. If not
            specified, the schema will be inferred by the native library.
        backend: specifies which eager backend instantiate to. Only
            necessary if inputs are not Narwhals Series.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
    """
    return _stableify(nw_f.from_dict(data, schema, backend=backend))


from_dicts: Final = DataFrame.from_dicts


def from_numpy(
    data: _2DArray,
    schema: Mapping[str, DType] | Schema | Sequence[str] | None = None,
    *,
    backend: IntoBackend[EagerAllowed],
) -> DataFrame[Any]:
    """Construct a DataFrame from a NumPy ndarray.

    Notes:
        Only row orientation is currently supported.

        For pandas-like dataframes, conversion to schema is applied after dataframe
        creation.

    Arguments:
        data: Two-dimensional data represented as a NumPy ndarray.
        schema: The DataFrame schema as Schema, dict of {name: type}, or a sequence of str.
        backend: specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
    """
    return _stableify(nw_f.from_numpy(data, schema, backend=backend))


def read_csv(
    source: str, *, backend: IntoBackend[EagerAllowed], **kwargs: Any
) -> DataFrame[Any]:
    """Read a CSV file into a DataFrame.

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        kwargs: Extra keyword arguments which are passed to the native CSV reader.
            For example, you could use
            `nw.read_csv('file.csv', backend='pandas', engine='pyarrow')`.
    """
    return _stableify(nw_f.read_csv(source, backend=backend, **kwargs))


def scan_csv(
    source: str, *, backend: IntoBackend[Backend], **kwargs: Any
) -> LazyFrame[Any]:
    """Lazily read from a CSV file.

    For the libraries that do not support lazy dataframes, the function reads
    a csv file eagerly and then converts the resulting dataframe to a lazyframe.

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        kwargs: Extra keyword arguments which are passed to the native CSV reader.
            For example, you could use
            `nw.scan_csv('file.csv', backend=pd, engine='pyarrow')`.
    """
    return _stableify(nw_f.scan_csv(source, backend=backend, **kwargs))


def read_parquet(
    source: str, *, backend: IntoBackend[EagerAllowed], **kwargs: Any
) -> DataFrame[Any]:
    """Read into a DataFrame from a parquet file.

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        kwargs: Extra keyword arguments which are passed to the native parquet reader.
            For example, you could use
            `nw.read_parquet('file.parquet', backend=pd, engine='pyarrow')`.
    """
    return _stableify(nw_f.read_parquet(source, backend=backend, **kwargs))


def scan_parquet(
    source: str, *, backend: IntoBackend[Backend], **kwargs: Any
) -> LazyFrame[Any]:
    """Lazily read from a parquet file.

    For the libraries that do not support lazy dataframes, the function reads
    a parquet file eagerly and then converts the resulting dataframe to a lazyframe.

    Note:
        Spark like backends require a session object to be passed in `kwargs`.

        For instance:

        ```py
        import narwhals as nw
        from sqlframe.duckdb import DuckDBSession

        nw.scan_parquet(source, backend="sqlframe", session=DuckDBSession())
        ```

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN`, `CUDF`, `PYSPARK` or `SQLFRAME`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"`, `"cudf"`,
                `"pyspark"` or `"sqlframe"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin`, `cudf`,
                `pyspark.sql` or `sqlframe`.
        kwargs: Extra keyword arguments which are passed to the native parquet reader.
            For example, you could use
            `nw.scan_parquet('file.parquet', backend=pd, engine='pyarrow')`.
    """
    return _stableify(nw_f.scan_parquet(source, backend=backend, **kwargs))


__all__ = [
    "Array",
    "Binary",
    "Boolean",
    "Categorical",
    "DataFrame",
    "Date",
    "Datetime",
    "Decimal",
    "Duration",
    "Enum",
    "Expr",
    "Field",
    "Float32",
    "Float64",
    "Implementation",
    "Int8",
    "Int16",
    "Int32",
    "Int64",
    "Int128",
    "LazyFrame",
    "List",
    "Object",
    "Schema",
    "Series",
    "String",
    "Struct",
    "Time",
    "UInt8",
    "UInt16",
    "UInt32",
    "UInt64",
    "UInt128",
    "Unknown",
    "all",
    "all_horizontal",
    "any_horizontal",
    "coalesce",
    "col",
    "concat",
    "concat_str",
    "dependencies",
    "dtypes",
    "dtypes",
    "exceptions",
    "exclude",
    "format",
    "from_arrow",
    "from_dict",
    "from_dicts",
    "from_native",
    "from_numpy",
    "generate_temporary_column_name",
    "get_native_namespace",
    "is_ordered_categorical",
    "len",
    "lit",
    "max",
    "max_horizontal",
    "maybe_align_index",
    "maybe_convert_dtypes",
    "maybe_get_index",
    "maybe_reset_index",
    "maybe_set_index",
    "mean",
    "mean_horizontal",
    "median",
    "min",
    "min_horizontal",
    "narwhalify",
    "new_series",
    "nth",
    "read_csv",
    "read_parquet",
    "scan_csv",
    "scan_parquet",
    "selectors",
    "selectors",
    "show_versions",
    "sum",
    "sum_horizontal",
    "to_native",
    "to_py_scalar",
    "when",
]
