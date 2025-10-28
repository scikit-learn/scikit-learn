from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Generic, Literal, overload

from narwhals._utils import (
    Implementation,
    Version,
    _Implementation,
    _validate_rolling_arguments,
    ensure_type,
    generate_repr,
    is_compliant_series,
    is_eager_allowed,
    is_index_selector,
    qualified_type_name,
    supports_arrow_c_stream,
)
from narwhals.dependencies import is_numpy_array, is_numpy_array_1d, is_numpy_scalar
from narwhals.dtypes import _validate_dtype, _validate_into_dtype
from narwhals.exceptions import ComputeError, InvalidOperationError
from narwhals.series_cat import SeriesCatNamespace
from narwhals.series_dt import SeriesDateTimeNamespace
from narwhals.series_list import SeriesListNamespace
from narwhals.series_str import SeriesStringNamespace
from narwhals.series_struct import SeriesStructNamespace
from narwhals.translate import to_native
from narwhals.typing import IntoSeriesT

if TYPE_CHECKING:
    from types import ModuleType
    from typing import NoReturn

    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import Self

    from narwhals._compliant import CompliantSeries
    from narwhals._typing import EagerAllowed, IntoBackend
    from narwhals.dataframe import DataFrame, MultiIndexSelector
    from narwhals.dtypes import DType
    from narwhals.typing import (
        ClosedInterval,
        FillNullStrategy,
        IntoDType,
        ModeKeepStrategy,
        NonNestedLiteral,
        NumericLiteral,
        PythonLiteral,
        RankMethod,
        RollingInterpolationMethod,
        SingleIndexSelector,
        TemporalLiteral,
        _1DArray,
    )


class Series(Generic[IntoSeriesT]):
    """Narwhals Series, backed by a native series.

    Warning:
        This class is not meant to be instantiated directly - instead:

        - If the native object is a series from one of the supported backend (e.g.
            pandas.Series, polars.Series, pyarrow.ChunkedArray), you can use
            [`narwhals.from_native`][]:
            ```py
            narwhals.from_native(native_series, allow_series=True)
            narwhals.from_native(native_series, series_only=True)
            ```

        - If the object is a generic sequence (e.g. a list or a tuple of values), you can
            create a series via [`narwhals.new_series`][], e.g.:
            ```py
            narwhals.new_series(name="price", values=[10.5, 9.4, 1.2], backend="pandas")
            ```
    """

    _version: ClassVar[Version] = Version.MAIN

    @property
    def _compliant(self) -> CompliantSeries[IntoSeriesT]:
        return self._compliant_series

    @property
    def _dataframe(self) -> type[DataFrame[Any]]:
        from narwhals.dataframe import DataFrame

        return DataFrame

    def __init__(
        self, series: Any, *, level: Literal["full", "lazy", "interchange"]
    ) -> None:
        self._level: Literal["full", "lazy", "interchange"] = level
        if is_compliant_series(series):
            self._compliant_series: CompliantSeries[IntoSeriesT] = (
                series.__narwhals_series__()
            )
        else:  # pragma: no cover
            msg = f"Expected Polars Series or an object which implements `__narwhals_series__`, got: {type(series)}."
            raise AssertionError(msg)

    @classmethod
    def from_numpy(
        cls,
        name: str,
        values: _1DArray,
        dtype: IntoDType | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> Series[Any]:
        """Construct a Series from a NumPy ndarray.

        Arguments:
            name: Name of resulting Series.
            values: One-dimensional data represented as a NumPy ndarray.
            dtype: (Narwhals) dtype. If not provided, the native library
                may auto-infer it from `values`.
            backend: specifies which eager backend instantiate to.

                `backend` can be specified in various ways

                - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                    `POLARS`, `MODIN` or `CUDF`.
                - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
                - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

        Examples:
            >>> import numpy as np
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> arr = np.arange(5, 10)
            >>> nw.Series.from_numpy("arr", arr, dtype=nw.Int8, backend="polars")
            ┌──────────────────┐
            | Narwhals Series  |
            |------------------|
            |shape: (5,)       |
            |Series: 'arr' [i8]|
            |[                 |
            |        5         |
            |        6         |
            |        7         |
            |        8         |
            |        9         |
            |]                 |
            └──────────────────┘
        """
        if not is_numpy_array_1d(values):
            msg = "`from_numpy` only accepts 1D numpy arrays"
            raise ValueError(msg)
        if dtype:
            _validate_into_dtype(dtype)
        implementation = Implementation.from_backend(backend)
        if is_eager_allowed(implementation):
            ns = cls._version.namespace.from_backend(implementation).compliant
            compliant = ns.from_numpy(values).alias(name)
            if dtype:
                return cls(compliant.cast(dtype), level="full")
            return cls(compliant, level="full")
        msg = (
            f"{implementation} support in Narwhals is lazy-only, but `Series.from_numpy` is an eager-only function.\n\n"
            "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
            f"    nw.Series.from_numpy(arr, backend='pyarrow').to_frame().lazy('{implementation}')"
        )
        raise ValueError(msg)

    @classmethod
    def from_iterable(
        cls,
        name: str,
        values: Iterable[Any],
        dtype: IntoDType | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> Series[Any]:
        """Construct a Series from an iterable.

        Arguments:
            name: Name of resulting Series.
            values: One-dimensional data represented as an iterable.
            dtype: (Narwhals) dtype. If not provided, the native library
                may auto-infer it from `values`.
            backend: specifies which eager backend instantiate to.

                `backend` can be specified in various ways

                - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                    `POLARS`, `MODIN` or `CUDF`.
                - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
                - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> values = [4, 1, 3, 2]
            >>> nw.Series.from_iterable("a", values, dtype=nw.UInt32, backend="pandas")
            ┌──────────────────────┐
            |   Narwhals Series    |
            |----------------------|
            |0    4                |
            |1    1                |
            |2    3                |
            |3    2                |
            |Name: a, dtype: uint32|
            └──────────────────────┘
        """
        if is_numpy_array(values):
            return cls.from_numpy(name, values, dtype, backend=backend)
        if dtype:
            _validate_into_dtype(dtype)
        if not isinstance(values, Iterable):
            msg = f"Expected values to be an iterable, got: {qualified_type_name(values)!r}."
            raise TypeError(msg)
        implementation = Implementation.from_backend(backend)
        if is_eager_allowed(implementation):
            ns = cls._version.namespace.from_backend(implementation).compliant
            compliant = ns._series.from_iterable(
                values, context=ns, name=name, dtype=dtype
            )
            return cls(compliant, level="full")
        msg = (
            f"{implementation} support in Narwhals is lazy-only, but `Series.from_iterable` is an eager-only function.\n\n"
            "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
            f"    nw.Series.from_iterable('a', [1,2,3], backend='pyarrow').to_frame().lazy('{implementation}')"
        )
        raise ValueError(msg)

    implementation: _Implementation = _Implementation()
    """Return [`narwhals.Implementation`][] of native Series.

    This can be useful when you need to use special-casing for features outside of
    Narwhals' scope - for example, when dealing with pandas' Period Dtype.

    Examples:
        >>> import narwhals as nw
        >>> import pandas as pd
        >>> s_native = pd.Series([1, 2, 3])
        >>> s = nw.from_native(s_native, series_only=True)
        >>> s.implementation
        <Implementation.PANDAS: 'pandas'>
        >>> s.implementation.is_pandas()
        True
        >>> s.implementation.is_pandas_like()
        True
        >>> s.implementation.is_polars()
        False
    """

    def __bool__(self) -> NoReturn:
        msg = (
            "the truth value of a Series is ambiguous"
            "\n\n"
            "Here are some things you might want to try:\n"
            "- instead of `if s`, use `if not s.is_empty()`\n"
            "- instead of `s1 and s2`, use `s1 & s2`\n"
            "- instead of `s1 or s2`, use `s1 | s2`\n"
            "- instead of `s in [y, z]`, use `s.is_in([y, z])`\n"
        )
        raise TypeError(msg)

    def __array__(self, dtype: Any = None, copy: bool | None = None) -> _1DArray:  # noqa: FBT001
        return self._compliant_series.__array__(dtype=dtype, copy=copy)

    @overload
    def __getitem__(self, idx: SingleIndexSelector) -> Any: ...

    @overload
    def __getitem__(self, idx: MultiIndexSelector) -> Self: ...

    def __getitem__(self, idx: SingleIndexSelector | MultiIndexSelector) -> Any | Self:
        """Retrieve elements from the object using integer indexing or slicing.

        Arguments:
            idx: The index, slice, or sequence of indices to retrieve.

                - If `idx` is an integer, a single element is returned.
                - If `idx` is a slice, a sequence of integers, or another Series
                    (with integer values) a subset of the Series is returned.

        Returns:
            A single element if `idx` is an integer, else a subset of the Series.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> nw.from_native(s_native, series_only=True)[0]
            1

            >>> nw.from_native(s_native, series_only=True)[
            ...     :2
            ... ].to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                2
              ]
            ]
        """
        if isinstance(idx, int) or (
            is_numpy_scalar(idx) and idx.dtype.kind in {"i", "u"}
        ):
            idx = int(idx) if not isinstance(idx, int) else idx
            return self._compliant_series.item(idx)

        if isinstance(idx, self.to_native().__class__):
            idx = self._with_compliant(self._compliant_series._with_native(idx))

        if not is_index_selector(idx):
            msg = (
                f"Unexpected type for `Series.__getitem__`: {type(idx)}.\n\n"
                "Hints:\n"
                "- use `s.item` to select a single item.\n"
                "- Use `s[indices]` to select rows positionally.\n"
                "- Use `s.filter(mask)` to filter rows based on a boolean mask."
            )
            raise TypeError(msg)
        if isinstance(idx, Series):
            return self._with_compliant(self._compliant_series[idx._compliant_series])
        assert not isinstance(idx, int)  # noqa: S101  # help mypy
        return self._with_compliant(self._compliant_series[idx])

    def __native_namespace__(self) -> ModuleType:
        return self._compliant_series.__native_namespace__()

    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object:
        """Export a Series via the Arrow PyCapsule Interface.

        Narwhals doesn't implement anything itself here:

        - if the underlying series implements the interface, it'll return that
        - else, it'll call `to_arrow` and then defer to PyArrow's implementation

        See [PyCapsule Interface](https://arrow.apache.org/docs/dev/format/CDataInterface/PyCapsuleInterface.html)
        for more.
        """
        native_series = self._compliant_series.native
        if supports_arrow_c_stream(native_series):
            return native_series.__arrow_c_stream__(requested_schema=requested_schema)
        try:
            pa_version = Implementation.PYARROW._backend_version()
        except ModuleNotFoundError as exc:
            msg = f"'pyarrow>=16.0.0' is required for `Series.__arrow_c_stream__` for object of type {type(native_series)}"
            raise ModuleNotFoundError(msg) from exc
        if pa_version < (16, 0):  # pragma: no cover
            msg = f"'pyarrow>=16.0.0' is required for `Series.__arrow_c_stream__` for object of type {type(native_series)}"
            raise ModuleNotFoundError(msg)
        from narwhals._arrow.utils import chunked_array

        ca = chunked_array(self.to_arrow())
        return ca.__arrow_c_stream__(requested_schema=requested_schema)

    def to_native(self) -> IntoSeriesT:
        """Convert Narwhals series to native series.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i64]
            [
              1
              2
            ]
        """
        return self._compliant_series.native

    def scatter(self, indices: int | Sequence[int], values: Any) -> Self:
        """Set value(s) at given position(s).

        Arguments:
            indices: Position(s) to set items at.
            values: Values to set.

        Note:
            This method always returns a new Series, without modifying the original one.
            Using this function in a for-loop is an anti-pattern, we recommend building
            up your positions and values beforehand and doing an update in one go.

            For example, instead of

            ```python
            for i in [1, 3, 2]:
                value = some_function(i)
                s = s.scatter(i, value)
            ```

            prefer

            ```python
            positions = [1, 3, 2]
            values = [some_function(x) for x in positions]
            s = s.scatter(positions, values)
            ```

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> df_native = pa.table({"a": [1, 2, 3], "b": [4, 5, 6]})
            >>> df_nw = nw.from_native(df_native)
            >>> df_nw.with_columns(df_nw["a"].scatter([0, 1], [999, 888])).to_native()
            pyarrow.Table
            a: int64
            b: int64
            ----
            a: [[999,888,3]]
            b: [[4,5,6]]
        """
        return self._with_compliant(
            self._compliant_series.scatter(indices, self._extract_native(values))
        )

    @property
    def shape(self) -> tuple[int]:
        """Get the shape of the Series.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).shape
            (3,)
        """
        return (self._compliant_series.len(),)

    def _extract_native(self, arg: Any) -> Any:
        from narwhals.series import Series

        if isinstance(arg, Series):
            return arg._compliant_series
        return arg

    def _with_compliant(self, series: Any) -> Self:
        return self.__class__(series, level=self._level)

    def pipe(self, function: Callable[[Any], Self], *args: Any, **kwargs: Any) -> Self:
        """Pipe function call.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([1, 2, 3])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.pipe(lambda x: x + 2).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [i64]
            [
                3
                4
                5
            ]
        """
        return function(self, *args, **kwargs)

    def __repr__(self) -> str:  # pragma: no cover
        return generate_repr("Narwhals Series", self.to_native().__repr__())

    def __len__(self) -> int:
        return len(self._compliant_series)

    def len(self) -> int:
        r"""Return the number of elements in the Series.

        Null values count towards the total.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, None]])
            >>> nw.from_native(s_native, series_only=True).len()
            3
        """
        return len(self._compliant_series)

    @property
    def dtype(self) -> DType:
        """Get the data type of the Series.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).dtype
            Int64
        """
        return self._compliant_series.dtype

    @property
    def name(self) -> str:
        """Get the name of the Series.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series("foo", [1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).name
            'foo'
        """
        return self._compliant_series.name

    def ewm_mean(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Self:
        r"""Compute exponentially-weighted moving average.

        Arguments:
            com: Specify decay in terms of center of mass, $\gamma$, with <br> $\alpha = \frac{1}{1+\gamma}\forall\gamma\geq0$
            span: Specify decay in terms of span, $\theta$, with <br> $\alpha = \frac{2}{\theta + 1} \forall \theta \geq 1$
            half_life: Specify decay in terms of half-life, $\tau$, with <br> $\alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \tau } \right\} \forall \tau > 0$
            alpha: Specify smoothing factor alpha directly, $0 < \alpha \leq 1$.
            adjust: Divide by decaying adjustment factor in beginning periods to account for imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights $w_i = (1 - \alpha)^i$
                - When `adjust=False` the EW function is calculated recursively by
                  $$
                  y_0=x_0
                  $$
                  $$
                  y_t = (1 - \alpha)y_{t - 1} + \alpha x_t
                  $$
            min_samples: Minimum number of observations in window required to have a value (otherwise result is null).
            ignore_nulls: Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of $x_0$ and $x_2$ used in
                  calculating the final weighted average of $[x_0, None, x_2]$ are
                  $(1-\alpha)^2$ and $1$ if `adjust=True`, and
                  $(1-\alpha)^2$ and $\alpha$ if `adjust=False`.
                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  $x_0$ and $x_2$ used in calculating the final weighted
                  average of $[x_0, None, x_2]$ are
                  $1-\alpha$ and $1$ if `adjust=True`,
                  and $1-\alpha$ and $\alpha$ if `adjust=False`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series(name="a", data=[1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).ewm_mean(
            ...     com=1, ignore_nulls=False
            ... ).to_native()
            0    1.000000
            1    1.666667
            2    2.428571
            Name: a, dtype: float64
        """
        return self._with_compliant(
            self._compliant_series.ewm_mean(
                com=com,
                span=span,
                half_life=half_life,
                alpha=alpha,
                adjust=adjust,
                min_samples=min_samples,
                ignore_nulls=ignore_nulls,
            )
        )

    def cast(self, dtype: IntoDType) -> Self:
        """Cast between data types.

        Arguments:
            dtype: Data type that the object will be cast into.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[True, False, True]])
            >>> nw.from_native(s_native, series_only=True).cast(nw.Int64).to_native()
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                0,
                1
              ]
            ]
        """
        _validate_dtype(dtype)
        return self._with_compliant(self._compliant_series.cast(dtype))

    def to_frame(self) -> DataFrame[Any]:
        """Convert to dataframe.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series("a", [1, 2])
            >>> nw.from_native(s_native, series_only=True).to_frame().to_native()
            shape: (2, 1)
            ┌─────┐
            │ a   │
            │ --- │
            │ i64 │
            ╞═════╡
            │ 1   │
            │ 2   │
            └─────┘
        """
        return self._dataframe(self._compliant_series.to_frame(), level=self._level)

    def to_list(self) -> list[Any]:
        """Convert to list.

        Notes:
            This function converts to Python scalars. It's typically
            more efficient to keep your data in the format native to
            your original dataframe, so we recommend only calling this
            when you absolutely need to.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> nw.from_native(s_native, series_only=True).to_list()
            [1, 2, 3]
        """
        return self._compliant_series.to_list()

    def mean(self) -> float:
        """Reduce this Series to the mean value.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1.2, 4.2])
            >>> nw.from_native(s_native, series_only=True).mean()
            np.float64(2.7)
        """
        return self._compliant_series.mean()

    def median(self) -> float:
        """Reduce this Series to the median value.

        Notes:
            Results might slightly differ across backends due to differences in the underlying algorithms used to compute the median.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[5, 3, 8]])
            >>> nw.from_native(s_native, series_only=True).median()
            5.0
        """
        return self._compliant_series.median()

    def skew(self) -> float | None:
        """Calculate the sample skewness of the Series.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 1, 2, 10, 100])
            >>> nw.from_native(s_native, series_only=True).skew()
            1.4724267269058975

        Notes:
            The skewness is a measure of the asymmetry of the probability distribution.
            A perfectly symmetric distribution has a skewness of 0.
        """
        return self._compliant_series.skew()

    def kurtosis(self) -> float | None:
        """Compute the kurtosis (Fisher's definition) without bias correction.

        Kurtosis is the fourth central moment divided by the square of the variance.
        The Fisher's definition is used where 3.0 is subtracted from the result to give 0.0 for a normal distribution.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 1, 2, 10, 100])
            >>> nw.from_native(s_native, series_only=True).kurtosis()
            0.2106571340718002
        """
        return self._compliant_series.kurtosis()

    def count(self) -> int:
        """Returns the number of non-null elements in the Series.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, None]])
            >>> nw.from_native(s_native, series_only=True).count()
            2
        """
        return self._compliant_series.count()

    def any(self) -> bool:
        """Return whether any of the values in the Series are True.

        If there are no non-null elements, the result is `False`.

        Notes:
            Only works on Series of data type Boolean.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([False, True, False])
            >>> nw.from_native(s_native, series_only=True).any()
            np.True_
        """
        return self._compliant_series.any()

    def all(self) -> bool:
        """Return whether all values in the Series are True.

        If there are no non-null elements, the result is `True`.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[False, True, False]])
            >>> nw.from_native(s_native, series_only=True).all()
            False
        """
        return self._compliant_series.all()

    def min(self) -> Any:
        """Get the minimal value in this Series.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).min()
            1
        """
        return self._compliant_series.min()

    def max(self) -> Any:
        """Get the maximum value in this Series.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).max()
            np.int64(3)
        """
        return self._compliant_series.max()

    def arg_min(self) -> int:
        """Returns the index of the minimum value.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> nw.from_native(s_native, series_only=True).arg_min()
            0
        """
        return self._compliant_series.arg_min()

    def arg_max(self) -> int:
        """Returns the index of the maximum value.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).arg_max()
            2
        """
        return self._compliant_series.arg_max()

    def sum(self) -> float:
        """Reduce this Series to the sum value.

        If there are no non-null elements, the result is zero.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> nw.from_native(s_native, series_only=True).sum()
            6
        """
        return self._compliant_series.sum()

    def std(self, *, ddof: int = 1) -> float:
        """Get the standard deviation of this Series.

        Arguments:
            ddof: "Delta Degrees of Freedom": the divisor used in the calculation is N - ddof,
                    where N represents the number of elements.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).std()
            1.0
        """
        return self._compliant_series.std(ddof=ddof)

    def var(self, *, ddof: int = 1) -> float:
        """Get the variance of this Series.

        Arguments:
            ddof: "Delta Degrees of Freedom": the divisor used in the calculation is N - ddof,
                    where N represents the number of elements.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> nw.from_native(s_native, series_only=True).var()
            1.0
        """
        return self._compliant_series.var(ddof=ddof)

    def clip(
        self,
        lower_bound: Self | NumericLiteral | TemporalLiteral | None = None,
        upper_bound: Self | NumericLiteral | TemporalLiteral | None = None,
    ) -> Self:
        r"""Clip values in the Series.

        Arguments:
            lower_bound: Lower bound value.
            upper_bound: Upper bound value.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([-1, 1, -3, 3, -5, 5])
            >>> nw.from_native(s_native, series_only=True).clip(-1, 3).to_native()
            0   -1
            1    1
            2   -1
            3    3
            4   -1
            5    3
            dtype: int64
        """
        return self._with_compliant(
            self._compliant_series.clip(
                lower_bound=self._extract_native(lower_bound),
                upper_bound=self._extract_native(upper_bound),
            )
        )

    def first(self) -> PythonLiteral:
        """Get the first element of the Series.

        Returns:
            A scalar value or `None` if the Series is empty.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3])
            >>> s_nw = nw.from_native(s_native, series_only=True)
            >>> s_nw.first()
            1
            >>> s_nw.filter(s_nw > 5).first() is None
            True
        """
        return self._compliant_series.first()

    def last(self) -> PythonLiteral:
        """Get the last element of the Series.

        Returns:
            A scalar value or `None` if the Series is empty.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> s_nw = nw.from_native(s_native, series_only=True)
            >>> s_nw.last()
            3
            >>> s_nw.filter(s_nw > 5).last() is None
            True
        """
        return self._compliant_series.last()

    def is_in(self, other: Any) -> Self:
        """Check if the elements of this Series are in the other sequence.

        Arguments:
            other: Sequence of primitive type.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.is_in([3, 2, 8]).to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                false,
                true,
                true
              ]
            ]
        """
        return self._with_compliant(
            self._compliant_series.is_in(to_native(other, pass_through=True))
        )

    def arg_true(self) -> Self:
        """Find elements where boolean Series is True.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, None, None, 2])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).is_null().arg_true().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [u32]
            [
               1
               2
            ]
        """
        return self._with_compliant(self._compliant_series.arg_true())

    def drop_nulls(self) -> Self:
        """Drop null values.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([2, 4, None, 3, 5])
            >>> nw.from_native(s_native, series_only=True).drop_nulls().to_native()
            0    2.0
            1    4.0
            3    3.0
            4    5.0
            dtype: float64
        """
        return self._with_compliant(self._compliant_series.drop_nulls())

    def abs(self) -> Self:
        """Calculate the absolute value of each element.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[2, -4, 3]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).abs().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                2,
                4,
                3
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.abs())

    def cum_sum(self, *, reverse: bool = False) -> Self:
        """Calculate the cumulative sum.

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([2, 4, 3])
            >>> nw.from_native(s_native, series_only=True).cum_sum().to_native()
            0    2
            1    6
            2    9
            dtype: int64
        """
        return self._with_compliant(self._compliant_series.cum_sum(reverse=reverse))

    def unique(self, *, maintain_order: bool = False) -> Self:
        """Returns unique values of the series.

        Arguments:
            maintain_order: Keep the same order as the original series. This may be more
                expensive to compute. Settings this to `True` blocks the possibility
                to run on the streaming engine for Polars.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([2, 4, 4, 6])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.unique(
            ...     maintain_order=True
            ... ).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [i64]
            [
               2
               4
               6
            ]
        """
        return self._with_compliant(
            self._compliant_series.unique(maintain_order=maintain_order)
        )

    def diff(self) -> Self:
        """Calculate the difference with the previous element, for each element.

        Notes:
            pandas may change the dtype here, for example when introducing missing
            values in an integer column. To ensure, that the dtype doesn't change,
            you may want to use `fill_null` and `cast`. For example, to calculate
            the diff and fill missing values with `0` in a Int64 column, you could
            do:

                s.diff().fill_null(0).cast(nw.Int64)

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[2, 4, 3]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).diff().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                null,
                2,
                -1
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.diff())

    def shift(self, n: int) -> Self:
        """Shift values by `n` positions.

        Arguments:
            n: Number of indices to shift forward. If a negative value is passed,
                values are shifted in the opposite direction instead.

        Notes:
            pandas may change the dtype here, for example when introducing missing
            values in an integer column. To ensure, that the dtype doesn't change,
            you may want to use `fill_null` and `cast`. For example, to shift
            and fill missing values with `0` in a Int64 column, you could
            do:

                s.shift(1).fill_null(0).cast(nw.Int64)

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([2, 4, 3])
            >>> nw.from_native(s_native, series_only=True).shift(1).to_native()
            0    NaN
            1    2.0
            2    4.0
            dtype: float64
        """
        ensure_type(n, int, param_name="n")

        return self._with_compliant(self._compliant_series.shift(n))

    def sample(
        self,
        n: int | None = None,
        *,
        fraction: float | None = None,
        with_replacement: bool = False,
        seed: int | None = None,
    ) -> Self:
        """Sample randomly from this Series.

        Arguments:
            n: Number of items to return. Cannot be used with fraction.
            fraction: Fraction of items to return. Cannot be used with n.
            with_replacement: Allow values to be sampled more than once.
            seed: Seed for the random number generator. If set to None (default), a random
                seed is generated for each sample operation.

        Notes:
            The `sample` method returns a Series with a specified number of
            randomly selected items chosen from this Series.
            The results are not consistent across libraries.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3, 4])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.sample(
            ...     fraction=1.0, with_replacement=True
            ... ).to_native()  # doctest: +SKIP
            shape: (4,)
            Series: '' [i64]
            [
               1
               4
               3
               4
            ]
        """
        return self._with_compliant(
            self._compliant_series.sample(
                n=n, fraction=fraction, with_replacement=with_replacement, seed=seed
            )
        )

    def alias(self, name: str) -> Self:
        """Rename the Series.

        Notes:
            This method is very cheap, but does not guarantee that data
            will be copied. For example:

            ```python
            s1: nw.Series
            s2 = s1.alias("foo")
            arr = s2.to_numpy()
            arr[0] = 999
            ```

            may (depending on the backend, and on the version) result in
            `s1`'s data being modified. We recommend:

                - if you need to alias an object and don't need the original
                  one around any more, just use `alias` without worrying about it.
                - if you were expecting `alias` to copy data, then explicitly call
                  `.clone` before calling `alias`.

        Arguments:
            name: The new name.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3], name="foo")
            >>> nw.from_native(s_native, series_only=True).alias("bar").to_native()
            0    1
            1    2
            2    3
            Name: bar, dtype: int64
        """
        return self._with_compliant(self._compliant_series.alias(name=name))

    def rename(self, name: str) -> Self:
        """Rename the Series.

        Alias for `Series.alias()`.

        Notes:
            This method is very cheap, but does not guarantee that data
            will be copied. For example:

            ```python
            s1: nw.Series
            s2 = s1.rename("foo")
            arr = s2.to_numpy()
            arr[0] = 999
            ```

            may (depending on the backend, and on the version) result in
            `s1`'s data being modified. We recommend:

                - if you need to rename an object and don't need the original
                  one around any more, just use `rename` without worrying about it.
                - if you were expecting `rename` to copy data, then explicitly call
                  `.clone` before calling `rename`.

        Arguments:
            name: The new name.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series("foo", [1, 2, 3])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.rename("bar").to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: 'bar' [i64]
            [
               1
               2
               3
            ]
        """
        return self.alias(name=name)

    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any] | None = None,
        *,
        return_dtype: IntoDType | None = None,
    ) -> Self:
        """Replace all values by different values.

        This function must replace all non-null input values (else it raises an error).

        Arguments:
            old: Sequence of values to replace. It also accepts a mapping of values to
                their replacement as syntactic sugar for
                `replace_strict(old=list(mapping.keys()), new=list(mapping.values()))`.
            new: Sequence of values to replace by. Length must match the length of `old`.
            return_dtype: The data type of the resulting expression. If set to `None`
                (default), the data type is determined automatically based on the other
                inputs.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([3, 0, 1, 2], name="a")
            >>> nw.from_native(s_native, series_only=True).replace_strict(
            ...     [0, 1, 2, 3], ["zero", "one", "two", "three"], return_dtype=nw.String
            ... ).to_native()
            0    three
            1     zero
            2      one
            3      two
            Name: a, dtype: object
        """
        if new is None:
            if not isinstance(old, Mapping):
                msg = "`new` argument is required if `old` argument is not a Mapping type"
                raise TypeError(msg)

            new = list(old.values())
            old = list(old.keys())

        return self._with_compliant(
            self._compliant_series.replace_strict(old, new, return_dtype=return_dtype)
        )

    def sort(self, *, descending: bool = False, nulls_last: bool = False) -> Self:
        """Sort this Series. Place null values first.

        Arguments:
            descending: Sort in descending order.
            nulls_last: Place null values last instead of first.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([5, None, 1, 2])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.sort(descending=True).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (4,)
            Series: '' [i64]
            [
               null
               5
               2
               1
            ]
        """
        return self._with_compliant(
            self._compliant_series.sort(descending=descending, nulls_last=nulls_last)
        )

    def is_null(self) -> Self:
        """Returns a boolean Series indicating which values are null.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, None]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).is_null().to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                false,
                false,
                true
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.is_null())

    def is_nan(self) -> Self:
        """Returns a boolean Series indicating which values are NaN.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([0.0, None, 2.0], dtype="Float64")
            >>> nw.from_native(s_native, series_only=True).is_nan().to_native()
            0    False
            1     <NA>
            2    False
            dtype: boolean
        """
        return self._with_compliant(self._compliant_series.is_nan())

    def fill_null(
        self,
        value: Self | NonNestedLiteral = None,
        strategy: FillNullStrategy | None = None,
        limit: int | None = None,
    ) -> Self:
        """Fill null values using the specified value.

        Arguments:
            value: Value used to fill null values.
            strategy: Strategy used to fill null values.
            limit: Number of consecutive null values to fill when using the 'forward' or 'backward' strategy.

        Notes:
            - pandas handles null values differently from other libraries.
              See [null_handling](../concepts/null_handling.md/)
              for reference.
            - For pandas Series of `object` dtype, `fill_null` will not automatically change the
              Series' dtype as pandas used to do. Explicitly call `cast` if you want the dtype to change.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, None])
            >>>
            >>> nw.from_native(s_native, series_only=True).fill_null(5).to_native()
            0    1.0
            1    2.0
            2    5.0
            dtype: float64

            Or using a strategy:

            >>> nw.from_native(s_native, series_only=True).fill_null(
            ...     strategy="forward", limit=1
            ... ).to_native()
            0    1.0
            1    2.0
            2    2.0
            dtype: float64
        """
        if value is not None and strategy is not None:
            msg = "cannot specify both `value` and `strategy`"
            raise ValueError(msg)
        if value is None and strategy is None:
            msg = "must specify either a fill `value` or `strategy`"
            raise ValueError(msg)
        if strategy is not None and strategy not in {"forward", "backward"}:
            msg = f"strategy not supported: {strategy}"
            raise ValueError(msg)
        return self._with_compliant(
            self._compliant_series.fill_null(
                value=self._extract_native(value), strategy=strategy, limit=limit
            )
        )

    def fill_nan(self, value: float | None) -> Self:
        """Fill floating point NaN values with given value.

        Arguments:
            value: Value used to fill NaN values.

        Notes:
            This function only fills `'NaN'` values, not null ones, except for pandas
            which doesn't distinguish between them.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([1.0, 2.0, float("nan"), None])
            >>> result = nw.from_native(s_native, series_only=True).fill_nan(0)
            >>> result.to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (4,)
            Series: '' [f64]
            [
               1.0
               2.0
               0.0
               null
            ]
        """
        return self._with_compliant(
            self._compliant_series.fill_nan(value=self._extract_native(value))
        )

    def is_between(
        self,
        lower_bound: Any | Self,
        upper_bound: Any | Self,
        closed: ClosedInterval = "both",
    ) -> Self:
        """Get a boolean mask of the values that are between the given lower/upper bounds.

        Arguments:
            lower_bound: Lower bound value.
            upper_bound: Upper bound value.
            closed: Define which sides of the interval are closed (inclusive).

        Notes:
            If the value of the `lower_bound` is greater than that of the `upper_bound`,
            then the values will be False, as no value can satisfy the condition.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3, 4, 5]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.is_between(2, 4, "right").to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                false,
                false,
                true,
                true,
                false
              ]
            ]
        """
        return self._with_compliant(
            self._compliant_series.is_between(
                self._extract_native(lower_bound),
                self._extract_native(upper_bound),
                closed=closed,
            )
        )

    def n_unique(self) -> int:
        """Count the number of unique values.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 2, 3])
            >>> nw.from_native(s_native, series_only=True).n_unique()
            3
        """
        return self._compliant_series.n_unique()

    def to_numpy(self) -> _1DArray:
        """Convert to numpy.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3], name="a")
            >>> nw.from_native(s_native, series_only=True).to_numpy()
            array([1, 2, 3]...)
        """
        return self._compliant_series.to_numpy(None, copy=None)

    def to_pandas(self) -> pd.Series[Any]:
        """Convert to pandas Series.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series("a", [1, 2, 3])
            >>> nw.from_native(s_native, series_only=True).to_pandas()
            0    1
            1    2
            2    3
            Name: a, dtype: int64
        """
        return self._compliant_series.to_pandas()

    def to_polars(self) -> pl.Series:
        """Convert to polars Series.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).to_polars()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [i64]
            [
                1
                2
                3
            ]
        """
        return self._compliant_series.to_polars()

    def __add__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__add__(self._extract_native(other))
        )

    def __radd__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__radd__(self._extract_native(other))
        )

    def __sub__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__sub__(self._extract_native(other))
        )

    def __rsub__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__rsub__(self._extract_native(other))
        )

    def __mul__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__mul__(self._extract_native(other))
        )

    def __rmul__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__rmul__(self._extract_native(other))
        )

    def __truediv__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__truediv__(self._extract_native(other))
        )

    def __rtruediv__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__rtruediv__(self._extract_native(other))
        )

    def __floordiv__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__floordiv__(self._extract_native(other))
        )

    def __rfloordiv__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__rfloordiv__(self._extract_native(other))
        )

    def __pow__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__pow__(self._extract_native(other))
        )

    def __rpow__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__rpow__(self._extract_native(other))
        )

    def __mod__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__mod__(self._extract_native(other))
        )

    def __rmod__(self, other: object) -> Self:
        return self._with_compliant(
            self._compliant_series.__rmod__(self._extract_native(other))
        )

    def __eq__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_compliant(
            self._compliant_series.__eq__(self._extract_native(other))
        )

    def __ne__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_compliant(
            self._compliant_series.__ne__(self._extract_native(other))
        )

    def __gt__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__gt__(self._extract_native(other))
        )

    def __ge__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__ge__(self._extract_native(other))
        )

    def __lt__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__lt__(self._extract_native(other))
        )

    def __le__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__le__(self._extract_native(other))
        )

    def __and__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__and__(self._extract_native(other))
        )

    def __rand__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__rand__(self._extract_native(other))
        )

    def __or__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__or__(self._extract_native(other))
        )

    def __ror__(self, other: Any) -> Self:
        return self._with_compliant(
            self._compliant_series.__ror__(self._extract_native(other))
        )

    # unary
    def __invert__(self) -> Self:
        return self._with_compliant(self._compliant_series.__invert__())

    def filter(self, predicate: Any) -> Self:
        """Filter elements in the Series based on a condition.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([4, 10, 15, 34, 50])
            >>> s_nw = nw.from_native(s_native, series_only=True)
            >>> s_nw.filter(s_nw > 10).to_native()
            2    15
            3    34
            4    50
            dtype: int64
        """
        return self._with_compliant(
            self._compliant_series.filter(self._extract_native(predicate))
        )

    # --- descriptive ---
    def is_duplicated(self) -> Self:
        r"""Get a mask of all duplicated rows in the Series.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3, 1]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).is_duplicated().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                true,
                false,
                false,
                true
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.is_duplicated())

    def is_empty(self) -> bool:
        r"""Check if the series is empty.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3])
            >>> s_nw = nw.from_native(s_native, series_only=True)

            >>> s_nw.is_empty()
            False
            >>> s_nw.filter(s_nw > 10).is_empty()
            True
        """
        return self._compliant_series.is_empty()

    def is_unique(self) -> Self:
        r"""Get a mask of all unique rows in the Series.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3, 1])
            >>> nw.from_native(s_native, series_only=True).is_unique().to_native()
            0    False
            1     True
            2     True
            3    False
            dtype: bool
        """
        return self._with_compliant(self._compliant_series.is_unique())

    def null_count(self) -> int:
        r"""Count the number of null values.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, None, None]])
            >>> nw.from_native(s_native, series_only=True).null_count()
            2
        """
        return self._compliant_series.null_count()

    def is_first_distinct(self) -> Self:
        r"""Return a boolean mask indicating the first occurrence of each distinct value.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 1, 2, 3, 2])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).is_first_distinct().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (5,)
            Series: '' [bool]
            [
                true
                false
                true
                true
                false
            ]
        """
        return self._with_compliant(self._compliant_series.is_first_distinct())

    def is_last_distinct(self) -> Self:
        r"""Return a boolean mask indicating the last occurrence of each distinct value.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 1, 2, 3, 2])
            >>> nw.from_native(s_native, series_only=True).is_last_distinct().to_native()
            0    False
            1     True
            2    False
            3     True
            4     True
            dtype: bool
        """
        return self._with_compliant(self._compliant_series.is_last_distinct())

    def is_sorted(self, *, descending: bool = False) -> bool:
        r"""Check if the Series is sorted.

        Arguments:
            descending: Check if the Series is sorted in descending order.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[3, 2, 1]])
            >>> s_nw = nw.from_native(s_native, series_only=True)

            >>> s_nw.is_sorted(descending=False)
            False

            >>> s_nw.is_sorted(descending=True)
            True
        """
        return self._compliant_series.is_sorted(descending=descending)

    def value_counts(
        self,
        *,
        sort: bool = False,
        parallel: bool = False,
        name: str | None = None,
        normalize: bool = False,
    ) -> DataFrame[Any]:
        r"""Count the occurrences of unique values.

        Arguments:
            sort: Sort the output by count in descending order. If set to False (default),
                the order of the output is random.
            parallel: Execute the computation in parallel. Used for Polars only.
            name: Give the resulting count column a specific name; if `normalize` is True
                defaults to "proportion", otherwise defaults to "count".
            normalize: If true gives relative frequencies of the unique values

                - The original values as first column
                - Either count or proportion as second column, depending on normalize parameter.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 1, 2, 3, 2], name="s")
            >>> nw.from_native(s_native, series_only=True).value_counts(
            ...     sort=True
            ... ).to_native()
               s  count
            0  1      2
            1  2      2
            2  3      1
        """
        return self._dataframe(
            self._compliant_series.value_counts(
                sort=sort, parallel=parallel, name=name, normalize=normalize
            ),
            level=self._level,
        )

    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> float:
        """Get quantile value of the series.

        Note:
            pandas and Polars may have implementation differences for a given interpolation method.

        Arguments:
            quantile: Quantile between 0.0 and 1.0.
            interpolation: Interpolation method.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series(list(range(50)))
            >>> s_nw = nw.from_native(s_native, series_only=True)
            >>> [
            ...     s_nw.quantile(quantile=q, interpolation="nearest")
            ...     for q in (0.1, 0.25, 0.5, 0.75, 0.9)
            ... ]
            [5.0, 12.0, 25.0, 37.0, 44.0]
        """
        return self._compliant_series.quantile(
            quantile=quantile, interpolation=interpolation
        )

    def zip_with(self, mask: Self, other: Self) -> Self:
        """Take values from self or other based on the given mask.

        Where mask evaluates true, take values from self. Where mask evaluates false,
        take values from other.

        Arguments:
            mask: Boolean Series
            other: Series of same type.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> data_native = pa.chunked_array([[1, 2, 3, 4, 5]])
            >>> other_native = pa.chunked_array([[5, 4, 3, 2, 1]])
            >>> mask_native = pa.chunked_array([[True, False, True, False, True]])
            >>>
            >>> data_nw = nw.from_native(data_native, series_only=True)
            >>> other_nw = nw.from_native(other_native, series_only=True)
            >>> mask_nw = nw.from_native(mask_native, series_only=True)
            >>>
            >>> data_nw.zip_with(mask_nw, other_nw).to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                4,
                3,
                2,
                5
              ]
            ]
        """
        return self._with_compliant(
            self._compliant_series.zip_with(
                self._extract_native(mask), self._extract_native(other)
            )
        )

    def item(self, index: int | None = None) -> Any:
        r"""Return the Series as a scalar, or return the element at the given index.

        If no index is provided, this is equivalent to `s[0]`, with a check
        that the shape is (1,). With an index, this is equivalent to `s[index]`.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> nw.from_native(pl.Series("a", [1]), series_only=True).item()
            1

            >>> nw.from_native(pl.Series("a", [9, 8, 7]), series_only=True).item(-1)
            7
        """
        return self._compliant_series.item(index=index)

    def head(self, n: int = 10) -> Self:
        r"""Get the first `n` rows.

        Arguments:
            n: Number of rows to return.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series(list(range(10)))
            >>> nw.from_native(s_native, series_only=True).head(3).to_native()
            0    0
            1    1
            2    2
            dtype: int64
        """
        return self._with_compliant(self._compliant_series.head(n))

    def tail(self, n: int = 10) -> Self:
        r"""Get the last `n` rows.

        Arguments:
            n: Number of rows to return.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([list(range(10))])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.tail(3).to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                7,
                8,
                9
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.tail(n))

    def round(self, decimals: int = 0) -> Self:
        r"""Round underlying floating point data by `decimals` digits.

        Arguments:
            decimals: Number of decimals to round by.

        Notes:
            For values exactly halfway between rounded decimal values pandas behaves differently than Polars and Arrow.

            pandas rounds to the nearest even value (e.g. -0.5 and 0.5 round to 0.0, 1.5 and 2.5 round to 2.0, 3.5 and
            4.5 to 4.0, etc..).

            Polars and Arrow round away from 0 (e.g. -0.5 to -1.0, 0.5 to 1.0, 1.5 to 2.0, 2.5 to 3.0, etc..).

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1.12345, 2.56789, 3.901234])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.round(1).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [f64]
            [
               1.1
               2.6
               3.9
            ]
        """
        return self._with_compliant(self._compliant_series.round(decimals))

    def floor(self) -> Self:
        r"""Compute the numerical floor.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([[1.1, 4.3, -1.3]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.floor().to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                4,
                -2
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.floor())

    def ceil(self) -> Self:
        r"""Compute the numerical ceiling.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([[1.1, 4.3, -1.3]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.ceil().to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                2,
                5,
                -1
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.ceil())

    def to_dummies(
        self, *, separator: str = "_", drop_first: bool = False
    ) -> DataFrame[Any]:
        r"""Get dummy/indicator variables.

        Arguments:
            separator: Separator/delimiter used when generating column names.
            drop_first: Remove the first category from the variable being encoded.

        Notes:
            pandas and Polars handle null values differently. Polars distinguishes
            between NaN and Null, whereas pandas doesn't.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1, 2, 3], name="a")
            >>> s_nw = nw.from_native(s_native, series_only=True)

            >>> s_nw.to_dummies(drop_first=False).to_native()
               a_1  a_2  a_3
            0    1    0    0
            1    0    1    0
            2    0    0    1

            >>> s_nw.to_dummies(drop_first=True).to_native()
               a_2  a_3
            0    0    0
            1    1    0
            2    0    1
        """
        return self._dataframe(
            self._compliant_series.to_dummies(separator=separator, drop_first=drop_first),
            level=self._level,
        )

    def gather_every(self, n: int, offset: int = 0) -> Self:
        r"""Take every nth value in the Series and return as new Series.

        Arguments:
            n: Gather every *n*-th row.
            offset: Starting index.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 2, 3, 4]])
            >>> nw.from_native(s_native, series_only=True).gather_every(
            ...     n=2, offset=1
            ... ).to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                2,
                4
              ]
            ]
        """
        return self._with_compliant(
            self._compliant_series.gather_every(n=n, offset=offset)
        )

    def to_arrow(self) -> pa.Array[Any]:
        r"""Convert to arrow.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 2, 3, 4])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).to_arrow()  # doctest:+NORMALIZE_WHITESPACE
            <pyarrow.lib.Int64Array object at ...>
            [
                1,
                2,
                3,
                4
            ]
        """
        return self._compliant_series.to_arrow()

    @overload
    def mode(self, *, keep: Literal["all"] = "all") -> Self: ...

    @overload
    def mode(self, *, keep: Literal["any"]) -> NonNestedLiteral: ...

    def mode(self, *, keep: ModeKeepStrategy = "all") -> Self | NonNestedLiteral:
        r"""Compute the most occurring value(s).

        Can return multiple values.

        Note:
            For `keep="any"` a scalar is returned, while for `keep="all"` a Series in
            returned even in the case of unimodal values.

        Arguments:
            keep: Whether to keep all modes or any mode found. Remark that `keep='any'`
                is not deterministic for multimodal values.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([1, 1, 2, 2, 3])
            >>> nw.from_native(s_native, series_only=True).mode().sort().to_native()
            0    1
            1    2
            dtype: int64
        """
        _supported_keep_values = ("all", "any")
        if keep not in _supported_keep_values:  # pragma: no cover
            msg = f"`keep` must be one of {_supported_keep_values}, found '{keep}'"
            raise ValueError(msg)

        result = self._with_compliant(self._compliant_series.mode(keep=keep))
        return result.item(0) if keep == "any" else result

    def is_finite(self) -> Self:
        """Returns a boolean Series indicating which values are finite.

        Warning:
            Different backend handle null values differently. `is_finite` will return
            False for NaN and Null's in the Dask and pandas non-nullable backend, while
            for Polars, PyArrow and pandas nullable backends null values are kept as such.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[float("nan"), float("inf"), 2.0, None]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).is_finite().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                false,
                false,
                true,
                null
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.is_finite())

    def cum_count(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative count of the non-null values in the series.

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series(["x", "k", None, "d"])
            >>> nw.from_native(s_native, series_only=True).cum_count(
            ...     reverse=True
            ... ).to_native()  # doctest:+NORMALIZE_WHITESPACE
            shape: (4,)
            Series: '' [u32]
            [
                3
                2
                1
                1
            ]
        """
        return self._with_compliant(self._compliant_series.cum_count(reverse=reverse))

    def cum_min(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative min of the non-null values in the series.

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([3, 1, None, 2])
            >>> nw.from_native(s_native, series_only=True).cum_min().to_native()
            0    3.0
            1    1.0
            2    NaN
            3    1.0
            dtype: float64
        """
        return self._with_compliant(self._compliant_series.cum_min(reverse=reverse))

    def cum_max(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative max of the non-null values in the series.

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1, 3, None, 2]])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).cum_max().to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                3,
                null,
                3
              ]
            ]
        """
        return self._with_compliant(self._compliant_series.cum_max(reverse=reverse))

    def cum_prod(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative product of the non-null values in the series.

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1, 3, None, 2])
            >>> nw.from_native(
            ...     s_native, series_only=True
            ... ).cum_prod().to_native()  # doctest:+NORMALIZE_WHITESPACE
            shape: (4,)
            Series: '' [i64]
            [
               1
               3
               null
               6
            ]
        """
        return self._with_compliant(self._compliant_series.cum_prod(reverse=reverse))

    def rolling_sum(
        self, window_size: int, *, min_samples: int | None = None, center: bool = False
    ) -> Self:
        """Apply a rolling sum (moving sum) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their sum.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`
            center: Set the labels at the center of the window.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1.0, 2.0, 3.0, 4.0])
            >>> nw.from_native(s_native, series_only=True).rolling_sum(
            ...     window_size=2
            ... ).to_native()
            0    NaN
            1    3.0
            2    5.0
            3    7.0
            dtype: float64
        """
        window_size, min_samples_int = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        if len(self) == 0:  # pragma: no cover
            return self

        return self._with_compliant(
            self._compliant_series.rolling_sum(
                window_size=window_size, min_samples=min_samples_int, center=center
            )
        )

    def rolling_mean(
        self, window_size: int, *, min_samples: int | None = None, center: bool = False
    ) -> Self:
        """Apply a rolling mean (moving mean) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their mean.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`
            center: Set the labels at the center of the window.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[1.0, 2.0, 3.0, 4.0]])
            >>> nw.from_native(s_native, series_only=True).rolling_mean(
            ...     window_size=2
            ... ).to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                null,
                1.5,
                2.5,
                3.5
              ]
            ]
        """
        window_size, min_samples = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        if len(self) == 0:  # pragma: no cover
            return self

        return self._with_compliant(
            self._compliant_series.rolling_mean(
                window_size=window_size, min_samples=min_samples, center=center
            )
        )

    def rolling_var(
        self,
        window_size: int,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Self:
        """Apply a rolling variance (moving variance) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their variance.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`.
            center: Set the labels at the center of the window.
            ddof: Delta Degrees of Freedom; the divisor for a length N window is N - ddof.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> s_native = pl.Series([1.0, 3.0, 1.0, 4.0])
            >>> nw.from_native(s_native, series_only=True).rolling_var(
            ...     window_size=2, min_samples=1
            ... ).to_native()  # doctest:+NORMALIZE_WHITESPACE
            shape: (4,)
            Series: '' [f64]
            [
               null
               2.0
               2.0
               4.5
            ]
        """
        window_size, min_samples = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        if len(self) == 0:  # pragma: no cover
            return self

        return self._with_compliant(
            self._compliant_series.rolling_var(
                window_size=window_size, min_samples=min_samples, center=center, ddof=ddof
            )
        )

    def rolling_std(
        self,
        window_size: int,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Self:
        """Apply a rolling standard deviation (moving standard deviation) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their standard deviation.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`.
            center: Set the labels at the center of the window.
            ddof: Delta Degrees of Freedom; the divisor for a length N window is N - ddof.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>>
            >>> s_native = pd.Series([1.0, 3.0, 1.0, 4.0])
            >>> nw.from_native(s_native, series_only=True).rolling_std(
            ...     window_size=2, min_samples=1
            ... ).to_native()
            0         NaN
            1    1.414214
            2    1.414214
            3    2.121320
            dtype: float64
        """
        window_size, min_samples = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        if len(self) == 0:  # pragma: no cover
            return self

        return self._with_compliant(
            self._compliant_series.rolling_std(
                window_size=window_size, min_samples=min_samples, center=center, ddof=ddof
            )
        )

    def __iter__(self) -> Iterator[Any]:
        yield from self._compliant_series.__iter__()

    def __contains__(self, other: Any) -> bool:
        return self._compliant_series.__contains__(other)

    def rank(self, method: RankMethod = "average", *, descending: bool = False) -> Self:
        """Assign ranks to data, dealing with ties appropriately.

        Notes:
            The resulting dtype may differ between backends.

        Arguments:
            method: The method used to assign ranks to tied elements.
                The following methods are available (default is 'average')

                - *"average"*: The average of the ranks that would have been assigned to
                    all the tied values is assigned to each value.
                - *"min"*: The minimum of the ranks that would have been assigned to all
                    the tied values is assigned to each value. (This is also referred to
                    as "competition" ranking.)
                - *"max"*: The maximum of the ranks that would have been assigned to all
                    the tied values is assigned to each value.
                - *"dense"*: Like "min", but the rank of the next highest element is
                    assigned the rank immediately after those assigned to the tied elements.
                - *"ordinal"*: All values are given a distinct rank, corresponding to the
                    order that the values occur in the Series.

            descending: Rank in descending order.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> s_native = pa.chunked_array([[3, 6, 1, 1, 6]])
            >>> nw.from_native(s_native, series_only=True).rank(
            ...     method="dense"
            ... ).to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                2,
                3,
                1,
                1,
                3
              ]
            ]
        """
        supported_rank_methods = {"average", "min", "max", "dense", "ordinal"}
        if method not in supported_rank_methods:
            msg = (
                "Ranking method must be one of {'average', 'min', 'max', 'dense', 'ordinal'}. "
                f"Found '{method}'"
            )
            raise ValueError(msg)

        return self._with_compliant(
            self._compliant_series.rank(method=method, descending=descending)
        )

    def hist(
        self,
        bins: list[float] | None = None,
        *,
        bin_count: int | None = None,
        include_breakpoint: bool = True,
    ) -> DataFrame[Any]:
        """Bin values into buckets and count their occurrences.

        Warning:
            This functionality is considered **unstable**. It may be changed at any point
            without it being considered a breaking change.

        Arguments:
            bins: A monotonically increasing sequence of values.
            bin_count: If no bins provided, this will be used to determine the distance of the bins.
            include_breakpoint: Include a column that shows the intervals as categories.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([1, 3, 8, 8, 2, 1, 3], name="a")
            >>> nw.from_native(s_native, series_only=True).hist(bin_count=4)
            ┌────────────────────┐
            | Narwhals DataFrame |
            |--------------------|
            |   breakpoint  count|
            |0        2.75      3|
            |1        4.50      2|
            |2        6.25      0|
            |3        8.00      2|
            └────────────────────┘
        """
        if bins is not None:
            if any(bins[i - 1] >= bins[i] for i in range(1, len(bins))):
                msg = "bins must increase monotonically"
                raise ComputeError(msg)
            if bin_count is not None:
                msg = f"can only provide one of `bin_count` or `bins`, got: {bin_count=}, {bins=}"
                raise ComputeError(msg)
            result = self._compliant_series.hist_from_bins(
                bins=bins, include_breakpoint=include_breakpoint
            )
        else:
            # polars (v1.20) sets bin=10 if neither are provided.
            default = 10
            bin_count = default if bin_count is None else bin_count
            result = self._compliant_series.hist_from_bin_count(
                bin_count=bin_count, include_breakpoint=include_breakpoint
            )

        return self._dataframe(result, level=self._level)

    def log(self, base: float = math.e) -> Self:
        r"""Compute the logarithm to a given base.

        Arguments:
            base: Given base, defaults to `e`

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([1, 2, 4], name="a")
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.log(base=2)
            ┌───────────────────────┐
            |    Narwhals Series    |
            |-----------------------|
            |0    0.0               |
            |1    1.0               |
            |2    2.0               |
            |Name: a, dtype: float64|
            └───────────────────────┘
        """
        return self._with_compliant(self._compliant_series.log(base=base))

    def exp(self) -> Self:
        r"""Compute the exponent.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([-1, 0, 1], name="a")
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.exp()
            ┌───────────────────────┐
            |    Narwhals Series    |
            |-----------------------|
            |0    0.367879          |
            |1    1.000000          |
            |2    2.718282          |
            |Name: a, dtype: float64|
            └───────────────────────┘
        """
        return self._with_compliant(self._compliant_series.exp())

    def sqrt(self) -> Self:
        r"""Compute the square root.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([1, 4, 9], name="a")
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.sqrt()
            ┌───────────────────────┐
            |    Narwhals Series    |
            |-----------------------|
            |0    1.0               |
            |1    2.0               |
            |2    3.0               |
            |Name: a, dtype: float64|
            └───────────────────────┘
        """
        return self._with_compliant(self._compliant_series.sqrt())

    def is_close(
        self,
        other: Self | NumericLiteral,
        *,
        abs_tol: float = 0.0,
        rel_tol: float = 1e-09,
        nans_equal: bool = False,
    ) -> Self:
        r"""Get a boolean mask of the values being close to the other values.

        Two values `a` and `b` are considered close if the following condition holds:

        $$
        |a-b| \le max \{ \text{rel\_tol} \cdot max \{ |a|, |b| \}, \text{abs\_tol} \}
        $$

        Arguments:
            other: Values to compare with.
            abs_tol: Absolute tolerance. This is the maximum allowed absolute difference
                between two values. Must be non-negative.
            rel_tol: Relative tolerance. This is the maximum allowed difference between
                two values, relative to the larger absolute value. Must be in the range
                [0, 1).
            nans_equal: Whether NaN values should be considered equal.

        Notes:
            The implementation of this method is symmetric and mirrors the behavior of
            `math.isclose`. Specifically note that this behavior is different to
            `numpy.isclose`.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> data = [1.0, float("inf"), 1.41, None, float("nan")]
            >>> s_native = pa.chunked_array([data])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.is_close(1.4, abs_tol=0.1).to_native()  # doctest:+ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                false,
                false,
                true,
                null,
                false
              ]
            ]
        """
        if not self.dtype.is_numeric():
            msg = (
                f"is_close operation not supported for dtype `{self.dtype}`\n\n"
                "Hint: `is_close` is only supported for numeric types"
            )
            raise InvalidOperationError(msg)

        if abs_tol < 0:
            msg = f"`abs_tol` must be non-negative but got {abs_tol}"
            raise ComputeError(msg)

        if not (0 <= rel_tol < 1):
            msg = f"`rel_tol` must be in the range [0, 1) but got {rel_tol}"
            raise ComputeError(msg)

        return self._with_compliant(
            self._compliant_series.is_close(
                self._extract_native(other),
                abs_tol=abs_tol,
                rel_tol=rel_tol,
                nans_equal=nans_equal,
            )
        )

    @property
    def str(self) -> SeriesStringNamespace[Self]:
        return SeriesStringNamespace(self)

    @property
    def dt(self) -> SeriesDateTimeNamespace[Self]:
        return SeriesDateTimeNamespace(self)

    @property
    def cat(self) -> SeriesCatNamespace[Self]:
        return SeriesCatNamespace(self)

    @property
    def list(self) -> SeriesListNamespace[Self]:
        return SeriesListNamespace(self)

    @property
    def struct(self) -> SeriesStructNamespace[Self]:
        return SeriesStructNamespace(self)
