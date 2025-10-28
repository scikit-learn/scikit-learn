from __future__ import annotations

from abc import abstractmethod
from itertools import chain
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    Generic,
    Literal,
    NoReturn,
    TypeVar,
    get_args,
    overload,
)

from narwhals._exceptions import issue_warning
from narwhals._expression_parsing import (
    ExprKind,
    check_expressions_preserve_length,
    is_into_expr_eager,
    is_scalar_like,
)
from narwhals._typing import Arrow, Pandas, _LazyAllowedImpl, _LazyFrameCollectImpl
from narwhals._utils import (
    Implementation,
    Version,
    _Implementation,
    can_lazyframe_collect,
    check_columns_exist,
    flatten,
    generate_repr,
    is_compliant_dataframe,
    is_compliant_lazyframe,
    is_eager_allowed,
    is_index_selector,
    is_iterator,
    is_lazy_allowed,
    is_list_of,
    is_sequence_like,
    is_slice_none,
    predicates_contains_list_of_bool,
    qualified_type_name,
    supports_arrow_c_stream,
    zip_strict,
)
from narwhals.dependencies import is_numpy_array_2d, is_pyarrow_table
from narwhals.exceptions import (
    ColumnNotFoundError,
    InvalidIntoExprError,
    InvalidOperationError,
    PerformanceWarning,
)
from narwhals.functions import _from_dict_no_backend, _is_into_schema
from narwhals.schema import Schema
from narwhals.series import Series
from narwhals.translate import to_native

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import Concatenate, ParamSpec, Self, TypeAlias

    from narwhals._compliant import CompliantDataFrame, CompliantLazyFrame
    from narwhals._compliant.typing import CompliantExprAny, EagerNamespaceAny
    from narwhals._translate import IntoArrowTable
    from narwhals._typing import EagerAllowed, IntoBackend, LazyAllowed, Polars
    from narwhals.group_by import GroupBy, LazyGroupBy
    from narwhals.typing import (
        AsofJoinStrategy,
        IntoDataFrame,
        IntoExpr,
        IntoFrame,
        IntoLazyFrame,
        IntoSchema,
        JoinStrategy,
        MultiColSelector as _MultiColSelector,
        MultiIndexSelector as _MultiIndexSelector,
        PivotAgg,
        SingleColSelector,
        SingleIndexSelector,
        SizeUnit,
        UniqueKeepStrategy,
        _2DArray,
    )

    PS = ParamSpec("PS")
    Incomplete: TypeAlias = Any

_FrameT = TypeVar("_FrameT", bound="IntoFrame")
LazyFrameT = TypeVar("LazyFrameT", bound="IntoLazyFrame")
DataFrameT = TypeVar("DataFrameT", bound="IntoDataFrame")
R = TypeVar("R")

MultiColSelector: TypeAlias = "_MultiColSelector[Series[Any]]"
MultiIndexSelector: TypeAlias = "_MultiIndexSelector[Series[Any]]"


class BaseFrame(Generic[_FrameT]):
    _compliant_frame: Any
    _level: Literal["full", "lazy", "interchange"]

    implementation: _Implementation = _Implementation()
    """Return [`narwhals.Implementation`][] of native frame.

    This can be useful when you need to use special-casing for features outside of
    Narwhals' scope - for example, when dealing with pandas' Period Dtype.

    Examples:
        >>> import narwhals as nw
        >>> import pandas as pd
        >>> df_native = pd.DataFrame({"a": [1, 2, 3]})
        >>> df = nw.from_native(df_native)
        >>> df.implementation
        <Implementation.PANDAS: 'pandas'>
        >>> df.implementation.is_pandas()
        True
        >>> df.implementation.is_pandas_like()
        True
        >>> df.implementation.is_polars()
        False
    """

    @property
    @abstractmethod
    def _compliant(self) -> Any: ...

    def __native_namespace__(self) -> ModuleType:
        return self._compliant_frame.__native_namespace__()  # type: ignore[no-any-return]

    def __narwhals_namespace__(self) -> Any:
        return self._compliant_frame.__narwhals_namespace__()

    def _with_compliant(self, df: Any) -> Self:
        # construct, preserving properties
        return self.__class__(df, level=self._level)  # type: ignore[call-arg]

    def _flatten_and_extract(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> tuple[list[CompliantExprAny], list[ExprKind]]:
        # Process `args` and `kwargs`, extracting underlying objects as we go.
        # NOTE: Strings are interpreted as column names.
        out_exprs = []
        out_kinds = []
        for expr in flatten(exprs):
            compliant_expr = self._extract_compliant(expr)
            out_exprs.append(compliant_expr)
            out_kinds.append(ExprKind.from_into_expr(expr, str_as_lit=False))
        for alias, expr in named_exprs.items():
            compliant_expr = self._extract_compliant(expr).alias(alias)
            out_exprs.append(compliant_expr)
            out_kinds.append(ExprKind.from_into_expr(expr, str_as_lit=False))
        return out_exprs, out_kinds

    @abstractmethod
    def _extract_compliant(self, arg: Any) -> Any:
        raise NotImplementedError

    def _extract_compliant_frame(self, other: Self | Any, /) -> Any:
        if isinstance(other, type(self)):
            return other._compliant_frame
        msg = f"Expected `other` to be a {qualified_type_name(self)!r}, got: {qualified_type_name(other)!r}"
        raise TypeError(msg)

    def _check_columns_exist(self, subset: Sequence[str]) -> ColumnNotFoundError | None:
        return check_columns_exist(subset, available=self.columns)

    @property
    def schema(self) -> Schema:
        return Schema(self._compliant_frame.schema.items())

    def collect_schema(self) -> Schema:
        native_schema = dict(self._compliant_frame.collect_schema())

        return Schema(native_schema)

    def pipe(
        self,
        function: Callable[Concatenate[Self, PS], R],
        *args: PS.args,
        **kwargs: PS.kwargs,
    ) -> R:
        return function(self, *args, **kwargs)

    def drop_nulls(self, subset: str | list[str] | None) -> Self:
        subset = [subset] if isinstance(subset, str) else subset
        return self._with_compliant(self._compliant_frame.drop_nulls(subset=subset))

    @property
    def columns(self) -> list[str]:
        return self._compliant_frame.columns  # type: ignore[no-any-return]

    def with_columns(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> Self:
        compliant_exprs, kinds = self._flatten_and_extract(*exprs, **named_exprs)
        compliant_exprs = [
            compliant_expr.broadcast(kind) if is_scalar_like(kind) else compliant_expr
            for compliant_expr, kind in zip_strict(compliant_exprs, kinds)
        ]
        return self._with_compliant(self._compliant_frame.with_columns(*compliant_exprs))

    def select(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> Self:
        flat_exprs = tuple(flatten(exprs))
        if flat_exprs and all(isinstance(x, str) for x in flat_exprs) and not named_exprs:
            # fast path!
            try:
                return self._with_compliant(
                    self._compliant_frame.simple_select(*flat_exprs)
                )
            except Exception as e:
                # Column not found is the only thing that can realistically be raised here.
                if error := self._check_columns_exist(flat_exprs):
                    raise error from e
                raise
        compliant_exprs, kinds = self._flatten_and_extract(*flat_exprs, **named_exprs)
        if compliant_exprs and all(is_scalar_like(kind) for kind in kinds):
            return self._with_compliant(self._compliant_frame.aggregate(*compliant_exprs))
        compliant_exprs = [
            compliant_expr.broadcast(kind) if is_scalar_like(kind) else compliant_expr
            for compliant_expr, kind in zip_strict(compliant_exprs, kinds)
        ]
        return self._with_compliant(self._compliant_frame.select(*compliant_exprs))

    def rename(self, mapping: dict[str, str]) -> Self:
        return self._with_compliant(self._compliant_frame.rename(mapping))

    def head(self, n: int) -> Self:
        return self._with_compliant(self._compliant_frame.head(n))

    def tail(self, n: int) -> Self:
        return self._with_compliant(self._compliant_frame.tail(n))

    def drop(self, *columns: Iterable[str], strict: bool) -> Self:
        return self._with_compliant(self._compliant_frame.drop(columns, strict=strict))

    def filter(
        self, *predicates: IntoExpr | Iterable[IntoExpr], **constraints: Any
    ) -> Self:
        from narwhals.functions import col

        flat_predicates = flatten(predicates)
        check_expressions_preserve_length(*flat_predicates, function_name="filter")
        plx = self.__narwhals_namespace__()
        compliant_predicates, _kinds = self._flatten_and_extract(*flat_predicates)
        compliant_constraints = (
            (col(name) == v)._to_compliant_expr(plx) for name, v in constraints.items()
        )
        predicate = plx.all_horizontal(
            *chain(compliant_predicates, compliant_constraints), ignore_nulls=False
        )
        return self._with_compliant(self._compliant_frame.filter(predicate))

    def sort(
        self,
        by: str | Iterable[str],
        *more_by: str,
        descending: bool | Sequence[bool] = False,
        nulls_last: bool = False,
    ) -> Self:
        by = flatten([*flatten([by]), *more_by])
        return self._with_compliant(
            self._compliant_frame.sort(*by, descending=descending, nulls_last=nulls_last)
        )

    def top_k(
        self, k: int, *, by: str | Iterable[str], reverse: bool | Sequence[bool] = False
    ) -> Self:
        flatten_by = flatten([by])
        return self._with_compliant(
            self._compliant_frame.top_k(k, by=flatten_by, reverse=reverse)
        )

    def join(
        self,
        other: Incomplete,
        on: str | list[str] | None,
        how: JoinStrategy,
        *,
        left_on: str | list[str] | None,
        right_on: str | list[str] | None,
        suffix: str,
    ) -> Self:
        _supported_joins = ("inner", "left", "full", "cross", "anti", "semi")
        on = [on] if isinstance(on, str) else on
        left_on = [left_on] if isinstance(left_on, str) else left_on
        right_on = [right_on] if isinstance(right_on, str) else right_on
        compliant = self._compliant_frame
        other = self._extract_compliant_frame(other)

        if how not in _supported_joins:
            msg = f"Only the following join strategies are supported: {_supported_joins}; found '{how}'."
            raise NotImplementedError(msg)
        if how == "cross":
            if left_on is not None or right_on is not None or on is not None:
                msg = "Can not pass `left_on`, `right_on` or `on` keys for cross join"
                raise ValueError(msg)
            result = compliant.join(
                other, how=how, left_on=None, right_on=None, suffix=suffix
            )
        elif on is None:
            if left_on is None or right_on is None:
                msg = f"Either (`left_on` and `right_on`) or `on` keys should be specified for {how}."
                raise ValueError(msg)
            if len(left_on) != len(right_on):
                msg = "`left_on` and `right_on` must have the same length."
                raise ValueError(msg)
            result = compliant.join(
                other, how=how, left_on=left_on, right_on=right_on, suffix=suffix
            )
        else:
            if left_on is not None or right_on is not None:
                msg = f"If `on` is specified, `left_on` and `right_on` should be None for {how}."
                raise ValueError(msg)
            result = compliant.join(
                other, how=how, left_on=on, right_on=on, suffix=suffix
            )
        return self._with_compliant(result)

    def gather_every(self, n: int, offset: int = 0) -> Self:
        return self._with_compliant(
            self._compliant_frame.gather_every(n=n, offset=offset)
        )

    def join_asof(
        self,
        other: Incomplete,
        *,
        left_on: str | None,
        right_on: str | None,
        on: str | None,
        by_left: str | list[str] | None,
        by_right: str | list[str] | None,
        by: str | list[str] | None,
        strategy: AsofJoinStrategy,
        suffix: str,
    ) -> Self:
        _supported_strategies = ("backward", "forward", "nearest")

        if strategy not in _supported_strategies:
            msg = f"Only the following strategies are supported: {_supported_strategies}; found '{strategy}'."
            raise NotImplementedError(msg)

        if (on is None) and (left_on is None or right_on is None):
            msg = "Either (`left_on` and `right_on`) or `on` keys should be specified."
            raise ValueError(msg)
        if (on is not None) and (left_on is not None or right_on is not None):
            msg = "If `on` is specified, `left_on` and `right_on` should be None."
            raise ValueError(msg)
        if (by is None) and (
            (by_left is None and by_right is not None)
            or (by_left is not None and by_right is None)
        ):
            msg = (
                "Can not specify only `by_left` or `by_right`, you need to specify both."
            )
            raise ValueError(msg)
        if (by is not None) and (by_left is not None or by_right is not None):
            msg = "If `by` is specified, `by_left` and `by_right` should be None."
            raise ValueError(msg)
        if on is not None:
            left_on = right_on = on
        if by is not None:
            by_left = by_right = by

        by_left = [by_left] if isinstance(by_left, str) else by_left
        by_right = [by_right] if isinstance(by_right, str) else by_right

        if (isinstance(by_left, list) and isinstance(by_right, list)) and (
            len(by_left) != len(by_right)
        ):
            msg = "`by_left` and `by_right` must have the same length."
            raise ValueError(msg)

        return self._with_compliant(
            self._compliant_frame.join_asof(
                self._extract_compliant_frame(other),
                left_on=left_on,
                right_on=right_on,
                by_left=by_left,
                by_right=by_right,
                strategy=strategy,
                suffix=suffix,
            )
        )

    def unpivot(
        self,
        on: str | list[str] | None,
        *,
        index: str | list[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        on = [on] if isinstance(on, str) else on
        index = [index] if isinstance(index, str) else index

        return self._with_compliant(
            self._compliant_frame.unpivot(
                on=on, index=index, variable_name=variable_name, value_name=value_name
            )
        )

    def __neq__(self, other: object) -> NoReturn:
        msg = (
            "DataFrame.__neq__ and LazyFrame.__neq__ are not implemented, please "
            "use expressions instead.\n\n"
            "Hint: instead of\n"
            "    df != 0\n"
            "you may want to use\n"
            "    df.select(nw.all() != 0)"
        )
        raise NotImplementedError(msg)

    def __eq__(self, other: object) -> NoReturn:
        msg = (
            "DataFrame.__eq__ and LazyFrame.__eq__ are not implemented, please "
            "use expressions instead.\n\n"
            "Hint: instead of\n"
            "    df == 0\n"
            "you may want to use\n"
            "    df.select(nw.all() == 0)"
        )
        raise NotImplementedError(msg)

    def explode(self, columns: str | Sequence[str], *more_columns: str) -> Self:
        to_explode = (
            [columns, *more_columns]
            if isinstance(columns, str)
            else [*columns, *more_columns]
        )

        return self._with_compliant(self._compliant_frame.explode(columns=to_explode))


class DataFrame(BaseFrame[DataFrameT]):
    """Narwhals DataFrame, backed by a native eager dataframe.

    Warning:
        This class is not meant to be instantiated directly - instead:

        - If the native object is a eager dataframe from one of the supported
            backend (e.g. pandas.DataFrame, polars.DataFrame, pyarrow.Table),
            you can use [`narwhals.from_native`][]:
            ```py
            narwhals.from_native(native_dataframe)
            narwhals.from_native(native_dataframe, eager_only=True)
            ```

        - If the object is a dictionary of column names and generic sequences mapping
            (e.g. `dict[str, list]`), you can create a DataFrame via
            [`narwhals.from_dict`][]:
            ```py
            narwhals.from_dict(
                data={"a": [1, 2, 3]},
                backend=narwhals.get_native_namespace(another_object),
            )
            ```
    """

    _version: ClassVar[Version] = Version.MAIN

    @property
    def _compliant(self) -> CompliantDataFrame[Any, Any, DataFrameT, Self]:
        return self._compliant_frame

    def _extract_compliant(self, arg: Any) -> Any:
        if is_into_expr_eager(arg):
            plx: EagerNamespaceAny = self.__narwhals_namespace__()
            return plx.parse_into_expr(arg, str_as_lit=False)
        raise InvalidIntoExprError.from_invalid_type(type(arg))

    @property
    def _series(self) -> type[Series[Any]]:
        return Series

    @property
    def _lazyframe(self) -> type[LazyFrame[Any]]:
        return LazyFrame

    def __init__(self, df: Any, *, level: Literal["full", "lazy", "interchange"]) -> None:
        self._level: Literal["full", "lazy", "interchange"] = level
        self._compliant_frame: CompliantDataFrame[Any, Any, DataFrameT, Self]
        if is_compliant_dataframe(df):
            self._compliant_frame = df.__narwhals_dataframe__()
        else:  # pragma: no cover
            msg = f"Expected an object which implements `__narwhals_dataframe__`, got: {type(df)}"
            raise AssertionError(msg)

    @classmethod
    def from_arrow(
        cls, native_frame: IntoArrowTable, *, backend: IntoBackend[EagerAllowed]
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

        Examples:
            >>> import pandas as pd
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> df_native = pd.DataFrame({"a": [1, 2], "b": [4.2, 5.1]})
            >>> nw.DataFrame.from_arrow(df_native, backend="polars")
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  shape: (2, 2)   |
            |  ┌─────┬─────┐   |
            |  │ a   ┆ b   │   |
            |  │ --- ┆ --- │   |
            |  │ i64 ┆ f64 │   |
            |  ╞═════╪═════╡   |
            |  │ 1   ┆ 4.2 │   |
            |  │ 2   ┆ 5.1 │   |
            |  └─────┴─────┘   |
            └──────────────────┘
        """
        if not (supports_arrow_c_stream(native_frame) or is_pyarrow_table(native_frame)):
            msg = f"Given object of type {type(native_frame)} does not support PyCapsule interface"
            raise TypeError(msg)
        implementation = Implementation.from_backend(backend)
        if is_eager_allowed(implementation):
            ns = cls._version.namespace.from_backend(implementation).compliant
            compliant = ns._dataframe.from_arrow(native_frame, context=ns)
            return cls(compliant, level="full")
        msg = (
            f"{implementation} support in Narwhals is lazy-only, but `DataFrame.from_arrow` is an eager-only function.\n\n"
            "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
            f"    nw.DataFrame.from_arrow(df, backend='pyarrow').lazy('{implementation}')"
        )
        raise ValueError(msg)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        schema: IntoSchema | None = None,
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

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> data = {"c": [5, 2], "d": [1, 4]}
            >>> nw.DataFrame.from_dict(data, backend="pandas")
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        c  d      |
            |     0  5  1      |
            |     1  2  4      |
            └──────────────────┘
        """
        if backend is None:
            data, backend = _from_dict_no_backend(data)
        implementation = Implementation.from_backend(backend)
        if is_eager_allowed(implementation):
            ns = cls._version.namespace.from_backend(implementation).compliant
            compliant = ns._dataframe.from_dict(data, schema=schema, context=ns)
            return cls(compliant, level="full")
        # NOTE: (#2786) needs resolving for extensions
        msg = (
            f"{implementation} support in Narwhals is lazy-only, but `DataFrame.from_dict` is an eager-only function.\n\n"
            "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
            f"    nw.DataFrame.from_dict({{'a': [1, 2]}}, backend='pyarrow').lazy('{implementation}')"
        )
        raise ValueError(msg)

    @classmethod
    def from_dicts(
        cls,
        data: Sequence[Mapping[str, Any]],
        schema: IntoSchema | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> DataFrame[Any]:
        """Instantiate DataFrame from a sequence of dictionaries representing rows.

        Notes:
            For pandas-like dataframes, conversion to schema is applied after dataframe
            creation.

        Arguments:
            data: Sequence with dictionaries mapping column name to value.
            schema: The DataFrame schema as Schema or dict of {name: type}. If not
                specified, the schema will be inferred by the native library.
            backend: Specifies which eager backend instantiate to.

                `backend` can be specified in various ways

                - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                    `POLARS`, `MODIN` or `CUDF`.
                - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
                - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

        Tip:
            If you expect non-uniform keys in `data`, consider passing `schema` for
            more consistent results, as **inference varies between backends**:

            - pandas uses all rows
            - polars uses the first 100 rows
            - pyarrow uses only the first row

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> data = [
            ...     {"item": "apple", "weight": 80, "price": 0.60},
            ...     {"item": "egg", "weight": 55, "price": 0.40},
            ... ]
            >>> nw.DataFrame.from_dicts(data, backend="polars")
            ┌──────────────────────────┐
            |    Narwhals DataFrame    |
            |--------------------------|
            |shape: (2, 3)             |
            |┌───────┬────────┬───────┐|
            |│ item  ┆ weight ┆ price │|
            |│ ---   ┆ ---    ┆ ---   │|
            |│ str   ┆ i64    ┆ f64   │|
            |╞═══════╪════════╪═══════╡|
            |│ apple ┆ 80     ┆ 0.6   │|
            |│ egg   ┆ 55     ┆ 0.4   │|
            |└───────┴────────┴───────┘|
            └──────────────────────────┘
        """
        implementation = Implementation.from_backend(backend)
        if is_eager_allowed(implementation):
            ns = cls._version.namespace.from_backend(implementation).compliant
            compliant = ns._dataframe.from_dicts(data, schema=schema, context=ns)
            return cls(compliant, level="full")
        # NOTE: (#2786) needs resolving for extensions
        msg = (
            f"{implementation} support in Narwhals is lazy-only, but `DataFrame.from_dicts` is an eager-only function.\n\n"
            "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
            f"    nw.DataFrame.from_dicts([{{'a': 1}}, {{'a': 2}}], backend='pyarrow').lazy('{implementation}')"
        )
        raise ValueError(msg)

    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        schema: IntoSchema | Sequence[str] | None = None,
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

        Examples:
            >>> import numpy as np
            >>> import polars as pl
            >>> import narwhals as nw
            >>>
            >>> arr = np.array([[5, 2, 1], [1, 4, 3]])
            >>> schema = {"c": nw.Int16(), "d": nw.Float32(), "e": nw.Int8()}
            >>> nw.DataFrame.from_numpy(arr, schema=schema, backend="polars")
            ┌───────────────────┐
            |Narwhals DataFrame |
            |-------------------|
            |shape: (2, 3)      |
            |┌─────┬─────┬─────┐|
            |│ c   ┆ d   ┆ e   │|
            |│ --- ┆ --- ┆ --- │|
            |│ i16 ┆ f32 ┆ i8  │|
            |╞═════╪═════╪═════╡|
            |│ 5   ┆ 2.0 ┆ 1   │|
            |│ 1   ┆ 4.0 ┆ 3   │|
            |└─────┴─────┴─────┘|
            └───────────────────┘
        """
        if not is_numpy_array_2d(data):
            msg = "`from_numpy` only accepts 2D numpy arrays"
            raise ValueError(msg)
        if not _is_into_schema(schema):
            msg = (
                "`schema` is expected to be one of the following types: "
                "Mapping[str, DType] | Schema | Sequence[str]. "
                f"Got {type(schema)}."
            )
            raise TypeError(msg)
        implementation = Implementation.from_backend(backend)
        if is_eager_allowed(implementation):
            ns = cls._version.namespace.from_backend(implementation).compliant
            return cls(ns.from_numpy(data, schema), level="full")
        msg = (
            f"{implementation} support in Narwhals is lazy-only, but `DataFrame.from_numpy` is an eager-only function.\n\n"
            "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
            f"    nw.DataFrame.from_numpy(arr, backend='pyarrow').lazy('{implementation}')"
        )
        raise ValueError(msg)

    def __len__(self) -> int:
        return self._compliant_frame.__len__()

    def __array__(self, dtype: Any = None, copy: bool | None = None) -> _2DArray:  # noqa: FBT001
        return self._compliant_frame.__array__(dtype, copy=copy)

    def __repr__(self) -> str:  # pragma: no cover
        return generate_repr("Narwhals DataFrame", self.to_native().__repr__())

    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object:
        """Export a DataFrame via the Arrow PyCapsule Interface.

        - if the underlying dataframe implements the interface, it'll return that
        - else, it'll call `to_arrow` and then defer to PyArrow's implementation

        See [PyCapsule Interface](https://arrow.apache.org/docs/dev/format/CDataInterface/PyCapsuleInterface.html)
        for more.
        """
        native_frame = self._compliant_frame._native_frame
        if supports_arrow_c_stream(native_frame):
            return native_frame.__arrow_c_stream__(requested_schema=requested_schema)
        try:
            pa_version = Implementation.PYARROW._backend_version()
        except ModuleNotFoundError as exc:
            msg = f"'pyarrow>=14.0.0' is required for `DataFrame.__arrow_c_stream__` for object of type {type(native_frame)}"
            raise ModuleNotFoundError(msg) from exc
        if pa_version < (14, 0):  # pragma: no cover
            msg = f"'pyarrow>=14.0.0' is required for `DataFrame.__arrow_c_stream__` for object of type {type(native_frame)}"
            raise ModuleNotFoundError(msg) from None
        pa_table = self.to_arrow()
        return pa_table.__arrow_c_stream__(requested_schema=requested_schema)  # type: ignore[no-untyped-call]

    def lazy(
        self,
        backend: IntoBackend[LazyAllowed] | None = None,
        *,
        session: Any | None = None,
    ) -> LazyFrame[Any]:
        """Restrict available API methods to lazy-only ones.

        If `backend` is specified, then a conversion between different backends
        might be triggered.

        If a library does not support lazy execution and `backend` is not specified,
        then this is will only restrict the API to lazy-only operations. This is useful
        if you want to ensure that you write dataframe-agnostic code which all has
        the possibility of running entirely lazily.

        Note:
            If `backend` is spark-like, then a valid `session` is required.

            For instance:

            ```py
            import narwhals as nw
            from sqlframe.duckdb import DuckDBSession

            df.lazy(backend=nw.Implementation.SQLFRAME, session=DuckDBSession())
            ```

        Arguments:
            backend: Which lazy backend collect to. This will be the underlying
                backend for the resulting Narwhals LazyFrame. If not specified, and the
                given library does not support lazy execution, then this will restrict
                the API to lazy-only operations.

                `backend` can be specified in various ways

                - As `Implementation.<BACKEND>` with `BACKEND` being `DASK`, `DUCKDB`,
                    `IBIS` or `POLARS`.
                - As a string: `"dask"`, `"duckdb"`, `"ibis"` or `"polars"`
                - Directly as a module `dask.dataframe`, `duckdb`, `ibis` or `polars`.
            session: Session to be used if backend is spark-like.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 2], "b": [4, 6]})
            >>> df = nw.from_native(df_native)

            If we call `df.lazy`, we get a `narwhals.LazyFrame` backed by a Polars
            LazyFrame.

            >>> df.lazy()  # doctest: +SKIP
            ┌─────────────────────────────┐
            |     Narwhals LazyFrame      |
            |-----------------------------|
            |<LazyFrame at 0x7F52B9937230>|
            └─────────────────────────────┘

            We can also pass DuckDB as the backend, and then we'll get a
            `narwhals.LazyFrame` backed by a `duckdb.DuckDBPyRelation`.

            >>> df.lazy(backend=nw.Implementation.DUCKDB)
            ┌──────────────────┐
            |Narwhals LazyFrame|
            |------------------|
            |┌───────┬───────┐ |
            |│   a   │   b   │ |
            |│ int64 │ int64 │ |
            |├───────┼───────┤ |
            |│     1 │     4 │ |
            |│     2 │     6 │ |
            |└───────┴───────┘ |
            └──────────────────┘
        """
        lazy = self._compliant_frame.lazy
        if backend is None:
            return self._lazyframe(lazy(None, session=session), level="lazy")
        lazy_backend = Implementation.from_backend(backend)
        if is_lazy_allowed(lazy_backend):
            return self._lazyframe(lazy(lazy_backend, session=session), level="lazy")
        msg = f"Not-supported backend.\n\nExpected one of {get_args(_LazyAllowedImpl)} or `None`, got {lazy_backend}"
        raise ValueError(msg)

    def to_native(self) -> DataFrameT:
        """Convert Narwhals DataFrame to native one.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"foo": [1, 2, 3], "bar": [6.0, 7.0, 8.0], "ham": ["a", "b", "c"]}
            ... )

            Calling `to_native` on a Narwhals DataFrame returns the native object:

            >>> nw.from_native(df_native).to_native()
               foo  bar ham
            0    1  6.0   a
            1    2  7.0   b
            2    3  8.0   c
        """
        return self._compliant_frame._native_frame

    def to_pandas(self) -> pd.DataFrame:
        """Convert this DataFrame to a pandas DataFrame.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"foo": [1, 2, 3], "bar": [6.0, 7.0, 8.0], "ham": ["a", "b", "c"]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.to_pandas()
               foo  bar ham
            0    1  6.0   a
            1    2  7.0   b
            2    3  8.0   c
        """
        return self._compliant_frame.to_pandas()

    def to_polars(self) -> pl.DataFrame:
        """Convert this DataFrame to a polars DataFrame.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> df = nw.from_native(df_native)
            >>> df.to_polars()
            shape: (2, 2)
            ┌─────┬─────┐
            │ foo ┆ bar │
            │ --- ┆ --- │
            │ i64 ┆ f64 │
            ╞═════╪═════╡
            │ 1   ┆ 6.0 │
            │ 2   ┆ 7.0 │
            └─────┴─────┘
        """
        return self._compliant_frame.to_polars()

    @overload
    def write_csv(self, file: None = None) -> str: ...

    @overload
    def write_csv(self, file: str | Path | BytesIO) -> None: ...

    def write_csv(self, file: str | Path | BytesIO | None = None) -> str | None:
        r"""Write dataframe to comma-separated values (CSV) file.

        Arguments:
            file: String, path object or file-like object to which the dataframe will be
                written. If None, the resulting csv format is returned as a string.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"foo": [1, 2, 3], "bar": [6.0, 7.0, 8.0], "ham": ["a", "b", "c"]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.write_csv()  # doctest: +SKIP
            'foo,bar,ham\n1,6.0,a\n2,7.0,b\n3,8.0,c\n'

            If we had passed a file name to `write_csv`, it would have been
            written to that file.
        """
        return self._compliant_frame.write_csv(file)

    def write_parquet(self, file: str | Path | BytesIO) -> None:
        """Write dataframe to parquet file.

        Arguments:
            file: String, path object or file-like object to which the dataframe will be
                written.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> df = nw.from_native(df_native)
            >>> df.write_parquet("out.parquet")  # doctest:+SKIP
        """
        self._compliant_frame.write_parquet(file)

    def to_numpy(self) -> _2DArray:
        """Convert this DataFrame to a NumPy ndarray.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [1, 2], "bar": [6.5, 7.0]})
            >>> df = nw.from_native(df_native)
            >>> df.to_numpy()
            array([[1. , 6.5],
                   [2. , 7. ]])
        """
        return self._compliant_frame.to_numpy(None, copy=None)

    @property
    def shape(self) -> tuple[int, int]:
        """Get the shape of the DataFrame.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [1, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.shape
            (2, 1)
        """
        return self._compliant_frame.shape

    def get_column(self, name: str) -> Series[Any]:
        """Get a single column by name.

        Arguments:
            name: The column name as a string.

        Notes:
            Although `name` is typed as `str`, pandas does allow non-string column
            names, and they will work when passed to this function if the
            `narwhals.DataFrame` is backed by a pandas dataframe with non-string
            columns. This function can only be used to extract a column by name, so
            there is no risk of ambiguity.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.get_column("a").to_native()
            0    1
            1    2
            Name: a, dtype: int64
        """
        return self._series(self._compliant_frame.get_column(name), level=self._level)

    def estimated_size(self, unit: SizeUnit = "b") -> int | float:
        """Return an estimation of the total (heap) allocated size of the `DataFrame`.

        Estimated size is given in the specified unit (bytes by default).

        Arguments:
            unit: 'b', 'kb', 'mb', 'gb', 'tb', 'bytes', 'kilobytes', 'megabytes',
                'gigabytes', or 'terabytes'.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> df = nw.from_native(df_native)
            >>> df.estimated_size()
            32
        """
        return self._compliant_frame.estimated_size(unit=unit)

    # `str` overlaps with `Sequence[str]`
    # We can ignore this but we must keep this overload ordering
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
    def __getitem__(  # noqa: C901, PLR0912
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
        """Extract column or slice of DataFrame.

        Arguments:
            item: How to slice dataframe. What happens depends on what is passed. It's easiest
                to explain by example. Suppose we have a Dataframe `df`

                - `df['a']` extracts column `'a'` and returns a `Series`.
                - `df[0:2]` extracts the first two rows and returns a `DataFrame`.
                - `df[0:2, 'a']` extracts the first two rows from column `'a'` and returns
                    a `Series`.
                - `df[0:2, 0]` extracts the first two rows from the first column and returns
                    a `Series`.
                - `df[[0, 1], [0, 1, 2]]` extracts the first two rows and the first three columns
                    and returns a `DataFrame`
                - `df[:, [0, 1, 2]]` extracts all rows from the first three columns and returns a
                  `DataFrame`.
                - `df[:, ['a', 'c']]` extracts all rows and columns `'a'` and `'c'` and returns a
                  `DataFrame`.
                - `df[['a', 'c']]` extracts all rows and columns `'a'` and `'c'` and returns a
                  `DataFrame`.
                - `df[0: 2, ['a', 'c']]` extracts the first two rows and columns `'a'` and `'c'` and
                    returns a `DataFrame`
                - `df[:, 0: 2]` extracts all rows from the first two columns and returns a `DataFrame`
                - `df[:, 'a': 'c']` extracts all rows and all columns positioned between `'a'` and `'c'`
                    _inclusive_ and returns a `DataFrame`. For example, if the columns are
                    `'a', 'd', 'c', 'b'`, then that would extract columns `'a'`, `'d'`, and `'c'`.

        Notes:
            - Integers are always interpreted as positions
            - Strings are always interpreted as column names.

            In contrast with Polars, pandas allows non-string column names.
            If you don't know whether the column name you're trying to extract
            is definitely a string (e.g. `df[df.columns[0]]`) then you should
            use `DataFrame.get_column` instead.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2]})
            >>> df = nw.from_native(df_native)
            >>> df["a"].to_native()
            0    1
            1    2
            Name: a, dtype: int64
        """
        from narwhals.series import Series

        msg = (
            f"Unexpected type for `DataFrame.__getitem__`, got: {type(item)}.\n\n"
            "Hints:\n"
            "- use `df.item` to select a single item.\n"
            "- Use `df[indices, :]` to select rows positionally.\n"
            "- Use `df.filter(mask)` to filter rows based on a boolean mask."
        )

        if isinstance(item, tuple):
            if len(item) > 2:
                tuple_msg = (
                    "Tuples cannot be passed to DataFrame.__getitem__ directly.\n\n"
                    "Hint: instead of `df[indices]`, did you mean `df[indices, :]`?"
                )
                raise TypeError(tuple_msg)
            rows = None if not item or is_slice_none(item[0]) else item[0]
            columns = None if len(item) < 2 or is_slice_none(item[1]) else item[1]
            if rows is None and columns is None:
                return self
        elif is_index_selector(item):
            rows = item
            columns = None
        elif is_sequence_like(item) or isinstance(item, (slice, str)):
            rows = None
            columns = item
        else:
            raise TypeError(msg)

        if isinstance(rows, str):
            raise TypeError(msg)

        compliant = self._compliant_frame

        if isinstance(columns, (int, str)):
            if isinstance(rows, int):
                return self.item(rows, columns)
            col_name = columns if isinstance(columns, str) else self.columns[columns]
            series = self.get_column(col_name)
            return series[rows] if rows is not None else series
        if isinstance(rows, Series):
            rows = rows._compliant_series
        if isinstance(columns, Series):
            columns = columns._compliant_series
        if rows is None:
            return self._with_compliant(compliant[:, columns])
        if columns is None:
            return self._with_compliant(compliant[rows, :])
        return self._with_compliant(compliant[rows, columns])

    def __contains__(self, key: str) -> bool:
        return key in self.columns

    @overload
    def to_dict(self, *, as_series: Literal[True] = ...) -> dict[str, Series[Any]]: ...
    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...
    @overload
    def to_dict(
        self, *, as_series: bool = True
    ) -> dict[str, Series[Any]] | dict[str, list[Any]]: ...
    def to_dict(
        self, *, as_series: bool = True
    ) -> dict[str, Series[Any]] | dict[str, list[Any]]:
        """Convert DataFrame to a dictionary mapping column name to values.

        Arguments:
            as_series: If set to true ``True``, then the values are Narwhals Series,
                    otherwise the values are Any.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"A": [1, 2], "fruits": ["banana", "apple"]})
            >>> df = nw.from_native(df_native)
            >>> df.to_dict(as_series=False)
            {'A': [1, 2], 'fruits': ['banana', 'apple']}
        """
        if as_series:
            return {
                key: self._series(value, level=self._level)
                for key, value in self._compliant_frame.to_dict(
                    as_series=as_series
                ).items()
            }
        return self._compliant_frame.to_dict(as_series=as_series)

    def row(self, index: int) -> tuple[Any, ...]:
        """Get values at given row.

        Warning:
            You should NEVER use this method to iterate over a DataFrame;
            if you require row-iteration you should strongly prefer use of iter_rows()
            instead.

        Arguments:
            index: Row number.

        Notes:
            cuDF doesn't support this method.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [1, 2], "b": [4, 5]})
            >>> nw.from_native(df_native).row(1)
            (<pyarrow.Int64Scalar: 2>, <pyarrow.Int64Scalar: 5>)
        """
        return self._compliant_frame.row(index)

    # inherited
    def pipe(
        self,
        function: Callable[Concatenate[Self, PS], R],
        *args: PS.args,
        **kwargs: PS.kwargs,
    ) -> R:
        """Pipe function call.

        Arguments:
            function: Function to apply.
            args: Positional arguments to pass to function.
            kwargs: Keyword arguments to pass to function.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2], "ba": [4, 5]})
            >>> nw.from_native(df_native).pipe(
            ...     lambda _df: _df.select(
            ...         [x for x in _df.columns if len(x) == 1]
            ...     ).to_native()
            ... )
               a
            0  1
            1  2
        """
        return super().pipe(function, *args, **kwargs)

    def drop_nulls(self, subset: str | list[str] | None = None) -> Self:
        """Drop rows that contain null values.

        Arguments:
            subset: Column name(s) for which null values are considered. If set to None
                (default), use all columns.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md)
            for reference.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [1.0, None], "ba": [1.0, 2.0]})
            >>> nw.from_native(df_native).drop_nulls().to_native()
            pyarrow.Table
            a: double
            ba: double
            ----
            a: [[1]]
            ba: [[1]]
        """
        return super().drop_nulls(subset=subset)

    def with_row_index(
        self, name: str = "index", *, order_by: str | Sequence[str] | None = None
    ) -> Self:
        """Insert column which enumerates rows.

        Arguments:
            name: The name of the column as a string. The default is "index".
            order_by: Column(s) to order by when computing the row index.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [1, 2], "b": [4, 5]})
            >>> nw.from_native(df_native).with_row_index().to_native()
            pyarrow.Table
            index: int64
            a: int64
            b: int64
            ----
            index: [[0,1]]
            a: [[1,2]]
            b: [[4,5]]
        """
        order_by_ = [order_by] if isinstance(order_by, str) else order_by
        return self._with_compliant(
            self._compliant_frame.with_row_index(name, order_by=order_by_)
        )

    @property
    def schema(self) -> Schema:
        r"""Get an ordered mapping of column names to their data type.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> nw.from_native(df_native).schema
            Schema({'foo': Int64, 'bar': Float64})
        """
        return super().schema

    def collect_schema(self) -> Schema:
        r"""Get an ordered mapping of column names to their data type.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> nw.from_native(df_native).collect_schema()
            Schema({'foo': Int64, 'bar': Float64})
        """
        return super().collect_schema()

    @property
    def columns(self) -> list[str]:
        """Get column names.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> nw.from_native(df_native).columns
            ['foo', 'bar']
        """
        return super().columns

    @overload
    def rows(self, *, named: Literal[False] = False) -> list[tuple[Any, ...]]: ...

    @overload
    def rows(self, *, named: Literal[True]) -> list[dict[str, Any]]: ...

    @overload
    def rows(self, *, named: bool) -> list[tuple[Any, ...]] | list[dict[str, Any]]: ...

    def rows(
        self, *, named: bool = False
    ) -> list[tuple[Any, ...]] | list[dict[str, Any]]:
        """Returns all data in the DataFrame as a list of rows of python-native values.

        Arguments:
            named: By default, each row is returned as a tuple of values given
                in the same order as the frame columns. Setting named=True will
                return rows of dictionaries instead.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> nw.from_native(df_native).rows()
            [(1, 6.0), (2, 7.0)]
        """
        return self._compliant_frame.rows(named=named)  # type: ignore[return-value]

    def iter_columns(self) -> Iterator[Series[Any]]:
        """Returns an iterator over the columns of this DataFrame.

        Yields:
            A Narwhals Series, backed by a native series.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> iter_columns = nw.from_native(df_native).iter_columns()
            >>> next(iter_columns)
            ┌───────────────────────┐
            |    Narwhals Series    |
            |-----------------------|
            |0    1                 |
            |1    2                 |
            |Name: foo, dtype: int64|
            └───────────────────────┘
            >>> next(iter_columns)
            ┌─────────────────────────┐
            |     Narwhals Series     |
            |-------------------------|
            |0    6.0                 |
            |1    7.0                 |
            |Name: bar, dtype: float64|
            └─────────────────────────┘
        """
        for series in self._compliant_frame.iter_columns():
            yield self._series(series, level=self._level)

    @overload
    def iter_rows(
        self, *, named: Literal[False] = ..., buffer_size: int = ...
    ) -> Iterator[tuple[Any, ...]]: ...

    @overload
    def iter_rows(
        self, *, named: Literal[True], buffer_size: int = ...
    ) -> Iterator[dict[str, Any]]: ...

    @overload
    def iter_rows(
        self, *, named: bool, buffer_size: int = ...
    ) -> Iterator[tuple[Any, ...]] | Iterator[dict[str, Any]]: ...

    def iter_rows(
        self, *, named: bool = False, buffer_size: int = 512
    ) -> Iterator[tuple[Any, ...]] | Iterator[dict[str, Any]]:
        """Returns an iterator over the DataFrame of rows of python-native values.

        Arguments:
            named: By default, each row is returned as a tuple of values given
                in the same order as the frame columns. Setting named=True will
                return rows of dictionaries instead.
            buffer_size: Determines the number of rows that are buffered
                internally while iterating over the data.
                See https://docs.pola.rs/api/python/stable/reference/dataframe/api/polars.DataFrame.iter_rows.html

        Notes:
            cuDF doesn't support this method.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> iter_rows = nw.from_native(df_native).iter_rows()
            >>> next(iter_rows)
            (1, 6.0)
            >>> next(iter_rows)
            (2, 7.0)
        """
        return self._compliant_frame.iter_rows(named=named, buffer_size=buffer_size)  # type: ignore[return-value]

    def with_columns(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> Self:
        r"""Add columns to this DataFrame.

        Added columns will replace existing columns with the same name.

        Arguments:
            *exprs: Column(s) to add, specified as positional arguments.
                     Accepts expression input. Strings are parsed as column names, other
                     non-expression inputs are parsed as literals.

            **named_exprs: Additional columns to add, specified as keyword arguments.
                            The columns will be renamed to the keyword used.

        Note:
            Creating a new DataFrame using this method does not create a new copy of
            existing data.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2], "b": [0.5, 4.0]})
            >>> (
            ...     nw.from_native(df_native)
            ...     .with_columns((nw.col("a") * 2).alias("a*2"))
            ...     .to_native()
            ... )
               a    b  a*2
            0  1  0.5    2
            1  2  4.0    4
        """
        return super().with_columns(*exprs, **named_exprs)

    def select(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> Self:
        r"""Select columns from this DataFrame.

        Arguments:
            *exprs: Column(s) to select, specified as positional arguments.
                     Accepts expression input. Strings are parsed as column names,
                     other non-expression inputs are parsed as literals.

            **named_exprs: Additional columns to select, specified as keyword arguments.
                            The columns will be renamed to the keyword used.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [1, 2], "b": [3, 4]})
            >>> nw.from_native(df_native).select("a", a_plus_1=nw.col("a") + 1)
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |pyarrow.Table     |
            |a: int64          |
            |a_plus_1: int64   |
            |----              |
            |a: [[1,2]]        |
            |a_plus_1: [[2,3]] |
            └──────────────────┘
        """
        return super().select(*exprs, **named_exprs)

    def rename(self, mapping: dict[str, str]) -> Self:
        """Rename column names.

        Arguments:
            mapping: Key value pairs that map from old name to new name.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, 2], "bar": [6, 7]})
            >>> nw.from_native(df_native).rename({"foo": "apple"}).to_native()
            pyarrow.Table
            apple: int64
            bar: int64
            ----
            apple: [[1,2]]
            bar: [[6,7]]
        """
        return super().rename(mapping)

    def head(self, n: int = 5) -> Self:
        """Get the first `n` rows.

        Arguments:
            n: Number of rows to return. If a negative value is passed, return all rows
                except the last `abs(n)`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2], "b": [0.5, 4.0]})
            >>> nw.from_native(df_native).head(1).to_native()
               a    b
            0  1  0.5
        """
        return super().head(n)

    def tail(self, n: int = 5) -> Self:
        """Get the last `n` rows.

        Arguments:
            n: Number of rows to return. If a negative value is passed, return all rows
                except the first `abs(n)`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2], "b": [0.5, 4.0]})
            >>> nw.from_native(df_native).tail(1)
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |       a    b     |
            |    1  2  4.0     |
            └──────────────────┘
        """
        return super().tail(n)

    def drop(self, *columns: str | Iterable[str], strict: bool = True) -> Self:
        """Remove columns from the dataframe.

        Arguments:
            *columns: Names of the columns that should be removed from the dataframe.
            strict: Validate that all column names exist in the schema and throw an
                exception if a column name does not exist in the schema.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"foo": [1, 2], "bar": [6.0, 7.0], "ham": ["a", "b"]}
            ... )
            >>> nw.from_native(df_native).drop("ham").to_native()
               foo  bar
            0    1  6.0
            1    2  7.0
        """
        return super().drop(*flatten(columns), strict=strict)

    def unique(
        self,
        subset: str | list[str] | None = None,
        *,
        keep: UniqueKeepStrategy = "any",
        maintain_order: bool = False,
        order_by: str | Sequence[str] | None = None,
    ) -> Self:
        """Drop duplicate rows from this dataframe.

        Arguments:
            subset: Column name(s) to consider when identifying duplicate rows.
            keep: {'first', 'last', 'any', 'none'}
                Which of the duplicate rows to keep.

                * 'any': Does not give any guarantee of which row is kept.
                        This allows more optimizations.
                * 'none': Don't keep duplicate rows.
                * 'first': Keep first unique row.
                * 'last': Keep last unique row.
            maintain_order: Keep the same order as the original DataFrame. This may be more
                expensive to compute.
            order_by: Column(s) to order by when computing the row index.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"foo": [1, 2], "bar": ["a", "a"], "ham": ["b", "b"]}
            ... )
            >>> nw.from_native(df_native).unique(["bar", "ham"]).to_native()
               foo bar ham
            0    1   a   b
        """
        if keep not in {"any", "none", "first", "last"}:
            msg = f"Expected {'any', 'none', 'first', 'last'}, got: {keep}"
            raise ValueError(msg)
        if isinstance(subset, str):
            subset = [subset]
        return self._with_compliant(
            self._compliant_frame.unique(
                subset, keep=keep, maintain_order=maintain_order, order_by=order_by
            )
        )

    def filter(
        self, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool], **constraints: Any
    ) -> Self:
        r"""Filter the rows in the DataFrame based on one or more predicate expressions.

        The original order of the remaining rows is preserved.

        Arguments:
            *predicates: Expression(s) that evaluates to a boolean Series. Can
                also be a boolean list(s).
            **constraints: Column filters; use `name = value` to filter columns by the supplied value.
                Each constraint will behave the same as `nw.col(name).eq(value)`, and will be implicitly
                joined with the other filter conditions using &.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"foo": [1, 2, 3], "bar": [6, 7, 8], "ham": ["a", "b", "c"]}
            ... )

            Filter on one condition

            >>> nw.from_native(df_native).filter(nw.col("foo") > 1).to_native()
               foo  bar ham
            1    2    7   b
            2    3    8   c

            Filter on multiple conditions with implicit `&`

            >>> nw.from_native(df_native).filter(
            ...     nw.col("foo") < 3, nw.col("ham") == "a"
            ... ).to_native()
               foo  bar ham
            0    1    6   a

            Filter on multiple conditions with `|`

            >>> nw.from_native(df_native).filter(
            ...     (nw.col("foo") == 1) | (nw.col("ham") == "c")
            ... ).to_native()
               foo  bar ham
            0    1    6   a
            2    3    8   c

            Filter using `**kwargs` syntax

            >>> nw.from_native(df_native).filter(foo=2, ham="b").to_native()
               foo  bar ham
            1    2    7   b
        """
        impl = self.implementation
        parsed_predicates = (
            self._series.from_iterable("", p, backend=impl) if is_list_of(p, bool) else p
            for p in predicates
        )
        return super().filter(*parsed_predicates, **constraints)

    @overload
    def group_by(
        self, *keys: IntoExpr | Iterable[IntoExpr], drop_null_keys: Literal[False] = ...
    ) -> GroupBy[Self]: ...

    @overload
    def group_by(
        self, *keys: str | Iterable[str], drop_null_keys: Literal[True]
    ) -> GroupBy[Self]: ...

    def group_by(
        self, *keys: IntoExpr | Iterable[IntoExpr], drop_null_keys: bool = False
    ) -> GroupBy[Self]:
        r"""Start a group by operation.

        Arguments:
            *keys: Column(s) to group by. Accepts expression input. Strings are parsed as
                column names.
            drop_null_keys: if True, then groups where any key is null won't be included
                in the result.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {
            ...         "a": ["a", "b", "a", "b", "c"],
            ...         "b": [1, 2, 1, 3, 3],
            ...         "c": [5, 4, 3, 2, 1],
            ...     }
            ... )

            Group by one column and compute the sum of another column

            >>> nw.from_native(df_native, eager_only=True).group_by("a").agg(
            ...     nw.col("b").sum()
            ... ).sort("a").to_native()
               a  b
            0  a  2
            1  b  5
            2  c  3

            Group by multiple columns and compute the max of another column

            >>> (
            ...     nw.from_native(df_native, eager_only=True)
            ...     .group_by(["a", "b"])
            ...     .agg(nw.max("c"))
            ...     .sort("a", "b")
            ...     .to_native()
            ... )
               a  b  c
            0  a  1  5
            1  b  2  4
            2  b  3  2
            3  c  3  1

            Expressions are also accepted.

            >>> nw.from_native(df_native, eager_only=True).group_by(
            ...     "a", nw.col("b") // 2
            ... ).agg(nw.col("c").mean()).to_native()
               a  b    c
            0  a  0  4.0
            1  b  1  3.0
            2  c  1  1.0
        """
        from narwhals.group_by import GroupBy

        flat_keys = flatten(keys)

        if all(isinstance(key, str) for key in flat_keys):
            return GroupBy(self, flat_keys, drop_null_keys=drop_null_keys)

        from narwhals import col
        from narwhals.expr import Expr
        from narwhals.series import Series

        key_is_expr_or_series = tuple(isinstance(k, (Expr, Series)) for k in flat_keys)

        if drop_null_keys and any(key_is_expr_or_series):
            msg = "drop_null_keys cannot be True when keys contains Expr or Series"
            raise NotImplementedError(msg)

        _keys = [
            k if is_expr else col(k)
            for k, is_expr in zip_strict(flat_keys, key_is_expr_or_series)
        ]
        expr_flat_keys, _kinds = self._flatten_and_extract(*_keys)
        check_expressions_preserve_length(*_keys, function_name="DataFrame.group_by")
        return GroupBy(self, expr_flat_keys, drop_null_keys=drop_null_keys)

    def sort(
        self,
        by: str | Iterable[str],
        *more_by: str,
        descending: bool | Sequence[bool] = False,
        nulls_last: bool = False,
    ) -> Self:
        r"""Sort the dataframe by the given columns.

        Arguments:
            by: Column(s) names to sort by.
            *more_by: Additional columns to sort by, specified as positional arguments.
            descending: Sort in descending order. When sorting by multiple columns, can be
                specified per column by passing a sequence of booleans.
            nulls_last: Place null values last.

        Note:
            Unlike Polars, it is not possible to specify a sequence of booleans for
            `nulls_last` in order to control per-column behaviour. Instead a single
            boolean is applied for all `by` columns.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"foo": [2, 1], "bar": [6.0, 7.0], "ham": ["a", "b"]}
            ... )
            >>> nw.from_native(df_native).sort("foo")
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |    foo  bar ham  |
            | 1    1  7.0   b  |
            | 0    2  6.0   a  |
            └──────────────────┘
        """
        return super().sort(by, *more_by, descending=descending, nulls_last=nulls_last)

    def top_k(
        self, k: int, *, by: str | Iterable[str], reverse: bool | Sequence[bool] = False
    ) -> Self:
        r"""Return the `k` largest rows.

        Non-null elements are always preferred over null elements,
        regardless of the value of reverse. The output is not guaranteed
        to be in any particular order, sort the outputs afterwards if you wish the output to be sorted.

        Arguments:
            k: Number of rows to return.
            by: Column(s) used to determine the top rows. Accepts expression input. Strings are parsed as column names.
            reverse: Consider the k smallest elements of the by column(s) (instead of the k largest).
                This can be specified per column by passing a sequence of booleans.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": ["a", "b", "a", "b", None, "c"], "b": [2, 1, 1, 3, 2, 1]}
            ... )
            >>> nw.from_native(df_native).top_k(4, by=["b", "a"])
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |          a  b    |
            |    3     b  3    |
            |    0     a  2    |
            |    4  None  2    |
            |    5     c  1    |
            └──────────────────┘
        """
        return super().top_k(k, by=by, reverse=reverse)

    def join(
        self,
        other: Self,
        on: str | list[str] | None = None,
        how: JoinStrategy = "inner",
        *,
        left_on: str | list[str] | None = None,
        right_on: str | list[str] | None = None,
        suffix: str = "_right",
    ) -> Self:
        r"""Join in SQL-like fashion.

        Arguments:
            other: DataFrame to join with.
            on: Name(s) of the join columns in both DataFrames. If set, `left_on` and
                `right_on` should be None.
            how: Join strategy.

                  * *inner*: Returns rows that have matching values in both tables.
                  * *left*: Returns all rows from the left table, and the matched rows from the right table.
                  * *full*: Returns all rows in both dataframes, with the suffix appended to the right join keys.
                  * *cross*: Returns the Cartesian product of rows from both tables.
                  * *semi*: Filter rows that have a match in the right table.
                  * *anti*: Filter rows that do not have a match in the right table.
            left_on: Join column of the left DataFrame.
            right_on: Join column of the right DataFrame.
            suffix: Suffix to append to columns with a duplicate name.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_1_native = pd.DataFrame({"id": ["a", "b"], "price": [6.0, 7.0]})
            >>> df_2_native = pd.DataFrame({"id": ["a", "b", "c"], "qty": [1, 2, 3]})
            >>> nw.from_native(df_1_native).join(nw.from_native(df_2_native), on="id")
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |   id  price  qty |
            | 0  a    6.0    1 |
            | 1  b    7.0    2 |
            └──────────────────┘
        """
        return super().join(
            other, how=how, left_on=left_on, right_on=right_on, on=on, suffix=suffix
        )

    def join_asof(
        self,
        other: Self,
        *,
        left_on: str | None = None,
        right_on: str | None = None,
        on: str | None = None,
        by_left: str | list[str] | None = None,
        by_right: str | list[str] | None = None,
        by: str | list[str] | None = None,
        strategy: AsofJoinStrategy = "backward",
        suffix: str = "_right",
    ) -> Self:
        """Perform an asof join.

        This is similar to a left-join except that we match on nearest key rather than equal keys.

        For Polars, both DataFrames must be sorted by the `on` key (within each `by` group
        if specified).

        Arguments:
            other: DataFrame to join with.
            left_on: Name(s) of the left join column(s).
            right_on: Name(s) of the right join column(s).
            on: Join column of both DataFrames. If set, left_on and right_on should be None.
            by_left: join on these columns before doing asof join.
            by_right: join on these columns before doing asof join.
            by: join on these columns before doing asof join.
            strategy: Join strategy. The default is "backward".
            suffix: Suffix to append to columns with a duplicate name.

                  * *backward*: selects the last row in the right DataFrame whose "on" key is less than or equal to the left's key.
                  * *forward*: selects the first row in the right DataFrame whose "on" key is greater than or equal to the left's key.
                  * *nearest*: search selects the last row in the right DataFrame whose value is nearest to the left's key.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> data_gdp = {
            ...     "datetime": [
            ...         datetime(2016, 1, 1),
            ...         datetime(2017, 1, 1),
            ...         datetime(2018, 1, 1),
            ...         datetime(2019, 1, 1),
            ...         datetime(2020, 1, 1),
            ...     ],
            ...     "gdp": [4164, 4411, 4566, 4696, 4827],
            ... }
            >>> data_population = {
            ...     "datetime": [
            ...         datetime(2016, 3, 1),
            ...         datetime(2018, 8, 1),
            ...         datetime(2019, 1, 1),
            ...     ],
            ...     "population": [82.19, 82.66, 83.12],
            ... }
            >>> gdp_native = pd.DataFrame(data_gdp)
            >>> population_native = pd.DataFrame(data_population)
            >>> gdp = nw.from_native(gdp_native)
            >>> population = nw.from_native(population_native)
            >>> population.join_asof(gdp, on="datetime", strategy="backward")
            ┌──────────────────────────────┐
            |      Narwhals DataFrame      |
            |------------------------------|
            |    datetime  population   gdp|
            |0 2016-03-01       82.19  4164|
            |1 2018-08-01       82.66  4566|
            |2 2019-01-01       83.12  4696|
            └──────────────────────────────┘
        """
        return super().join_asof(
            other,
            left_on=left_on,
            right_on=right_on,
            on=on,
            by_left=by_left,
            by_right=by_right,
            by=by,
            strategy=strategy,
            suffix=suffix,
        )

    # --- descriptive ---
    def is_duplicated(self) -> Series[Any]:
        r"""Get a mask of all duplicated rows in this DataFrame.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [2, 2, 2], "bar": [6.0, 6.0, 7.0]})
            >>> nw.from_native(df_native).is_duplicated()
            ┌───────────────┐
            |Narwhals Series|
            |---------------|
            |  0     True   |
            |  1     True   |
            |  2    False   |
            |  dtype: bool  |
            └───────────────┘
        """
        return ~self.is_unique()

    def is_empty(self) -> bool:
        r"""Check if the dataframe is empty.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [2, 2, 2], "bar": [6.0, 6.0, 7.0]})
            >>> nw.from_native(df_native).is_empty()
            False
        """
        return len(self) == 0

    def is_unique(self) -> Series[Any]:
        r"""Get a mask of all unique rows in this DataFrame.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [2, 2, 2], "bar": [6.0, 6.0, 7.0]})
            >>> nw.from_native(df_native).is_unique()
            ┌───────────────┐
            |Narwhals Series|
            |---------------|
            |  0    False   |
            |  1    False   |
            |  2     True   |
            |  dtype: bool  |
            └───────────────┘
        """
        return self._series(self._compliant_frame.is_unique(), level=self._level)

    def null_count(self) -> Self:
        r"""Create a new DataFrame that shows the null counts per column.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, None], "bar": [2, 3]})
            >>> nw.from_native(df_native).null_count()
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  pyarrow.Table   |
            |  foo: int64      |
            |  bar: int64      |
            |  ----            |
            |  foo: [[1]]      |
            |  bar: [[0]]      |
            └──────────────────┘
        """
        plx = self._compliant_frame.__narwhals_namespace__()
        result = self._compliant_frame.select(plx.all().null_count())
        return self._with_compliant(result)

    def item(self, row: int | None = None, column: int | str | None = None) -> Any:
        r"""Return the DataFrame as a scalar, or return the element at the given row/column.

        Arguments:
            row: The *n*-th row.
            column: The column selected via an integer or a string (column name).

        Notes:
            If row/col not provided, this is equivalent to df[0,0], with a check that the shape is (1,1).
            With row/col, this is equivalent to df[row,col].

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, None], "bar": [2, 3]})
            >>> nw.from_native(df_native).item(0, 1)
            2
        """
        return self._compliant_frame.item(row=row, column=column)

    def clone(self) -> Self:
        r"""Create a copy of this DataFrame."""
        return self._with_compliant(self._compliant_frame.clone())

    def gather_every(self, n: int, offset: int = 0) -> Self:
        r"""Take every nth row in the DataFrame and return as a new DataFrame.

        Arguments:
            n: Gather every *n*-th row.
            offset: Starting index.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"foo": [1, None, 2, 3]})
            >>> nw.from_native(df_native).gather_every(2)
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  pyarrow.Table   |
            |  foo: int64      |
            |  ----            |
            |  foo: [[1,2]]    |
            └──────────────────┘
        """
        return super().gather_every(n=n, offset=offset)

    def pivot(
        self,
        on: str | list[str],
        *,
        index: str | list[str] | None = None,
        values: str | list[str] | None = None,
        aggregate_function: PivotAgg | None = None,
        maintain_order: bool | None = None,
        sort_columns: bool = False,
        separator: str = "_",
    ) -> Self:
        r"""Create a spreadsheet-style pivot table as a DataFrame.

        Arguments:
            on: Name of the column(s) whose values will be used as the header of the
                output DataFrame.
            index: One or multiple keys to group by. If None, all remaining columns not
                specified on `on` and `values` will be used. At least one of `index` and
                `values` must be specified.
            values: One or multiple keys to group by. If None, all remaining columns not
                specified on `on` and `index` will be used. At least one of `index` and
                `values` must be specified.
            aggregate_function: Choose from

                - None: no aggregation takes place, will raise error if multiple values
                    are in group.
                - A predefined aggregate function string, one of
                    {'min', 'max', 'first', 'last', 'sum', 'mean', 'median', 'len'}
            maintain_order: Has no effect and is kept around only for backwards-compatibility.
            sort_columns: Sort the transposed columns by name. Default is by order of
                discovery.
            separator: Used as separator/delimiter in generated column names in case of
                multiple `values` columns.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> data = {
            ...     "ix": [1, 1, 2, 2, 1, 2],
            ...     "col": ["a", "a", "a", "a", "b", "b"],
            ...     "foo": [0, 1, 2, 2, 7, 1],
            ...     "bar": [0, 2, 0, 0, 9, 4],
            ... }
            >>> df_native = pd.DataFrame(data)
            >>> nw.from_native(df_native).pivot(
            ...     "col", index="ix", aggregate_function="sum"
            ... )
            ┌─────────────────────────────────┐
            |       Narwhals DataFrame        |
            |---------------------------------|
            |   ix  foo_a  foo_b  bar_a  bar_b|
            |0   1      1      7      2      9|
            |1   2      4      1      0      4|
            └─────────────────────────────────┘
        """
        if values is None and index is None:
            msg = "At least one of `values` and `index` must be passed"
            raise ValueError(msg)
        if maintain_order is not None:
            msg = (
                "`maintain_order` has no effect and is only kept around for backwards-compatibility. "
                "You can safely remove this argument."
            )
            issue_warning(msg, UserWarning)
        on = [on] if isinstance(on, str) else on
        values = [values] if isinstance(values, str) else values
        index = [index] if isinstance(index, str) else index

        return self._with_compliant(
            self._compliant_frame.pivot(
                on=on,
                index=index,
                values=values,
                aggregate_function=aggregate_function,
                sort_columns=sort_columns,
                separator=separator,
            )
        )

    def to_arrow(self) -> pa.Table:
        r"""Convert to arrow table.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [1, None], "bar": [2, 3]})
            >>> nw.from_native(df_native).to_arrow()
            pyarrow.Table
            foo: double
            bar: int64
            ----
            foo: [[1,null]]
            bar: [[2,3]]
        """
        return self._compliant_frame.to_arrow()

    def sample(
        self,
        n: int | None = None,
        *,
        fraction: float | None = None,
        with_replacement: bool = False,
        seed: int | None = None,
    ) -> Self:
        r"""Sample from this DataFrame.

        Arguments:
            n: Number of items to return. Cannot be used with fraction.
            fraction: Fraction of items to return. Cannot be used with n.
            with_replacement: Allow values to be sampled more than once.
            seed: Seed for the random number generator. If set to None (default), a random
                seed is generated for each sample operation.

        Notes:
            The results may not be consistent across libraries.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [1, 2, 3], "bar": [19, 32, 4]})
            >>> nw.from_native(df_native).sample(n=2)  # doctest:+SKIP
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |      foo  bar    |
            |   2    3    4    |
            |   1    2   32    |
            └──────────────────┘
        """
        return self._with_compliant(
            self._compliant_frame.sample(
                n=n, fraction=fraction, with_replacement=with_replacement, seed=seed
            )
        )

    def unpivot(
        self,
        on: str | list[str] | None = None,
        *,
        index: str | list[str] | None = None,
        variable_name: str = "variable",
        value_name: str = "value",
    ) -> Self:
        r"""Unpivot a DataFrame from wide to long format.

        Optionally leaves identifiers set.

        This function is useful to massage a DataFrame into a format where one or more
        columns are identifier variables (index) while all other columns, considered
        measured variables (on), are "unpivoted" to the row axis leaving just
        two non-identifier columns, 'variable' and 'value'.

        Arguments:
            on: Column(s) to use as values variables; if `on` is empty all columns that
                are not in `index` will be used.
            index: Column(s) to use as identifier variables.
            variable_name: Name to give to the `variable` column. Defaults to "variable".
            value_name: Name to give to the `value` column. Defaults to "value".

        Notes:
            If you're coming from pandas, this is similar to `pandas.DataFrame.melt`,
            but with `index` replacing `id_vars` and `on` replacing `value_vars`.
            In other frameworks, you might know this operation as `pivot_longer`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> data = {"a": ["x", "y", "z"], "b": [1, 3, 5], "c": [2, 4, 6]}
            >>> df_native = pd.DataFrame(data)
            >>> nw.from_native(df_native).unpivot(["b", "c"], index="a")
            ┌────────────────────┐
            | Narwhals DataFrame |
            |--------------------|
            |   a variable  value|
            |0  x        b      1|
            |1  y        b      3|
            |2  z        b      5|
            |3  x        c      2|
            |4  y        c      4|
            |5  z        c      6|
            └────────────────────┘
        """
        return super().unpivot(
            on=on, index=index, variable_name=variable_name, value_name=value_name
        )

    def explode(self, columns: str | Sequence[str], *more_columns: str) -> Self:
        """Explode the dataframe to long format by exploding the given columns.

        Notes:
            It is possible to explode multiple columns only if these columns must have
            matching element counts.

        Arguments:
            columns: Column names. The underlying columns being exploded must be of the `List` data type.
            *more_columns: Additional names of columns to explode, specified as positional arguments.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> data = {"a": ["x", "y"], "b": [[1, 2], [3]]}
            >>> df_native = pl.DataFrame(data)
            >>> nw.from_native(df_native).explode("b").to_native()
            shape: (3, 2)
            ┌─────┬─────┐
            │ a   ┆ b   │
            │ --- ┆ --- │
            │ str ┆ i64 │
            ╞═════╪═════╡
            │ x   ┆ 1   │
            │ x   ┆ 2   │
            │ y   ┆ 3   │
            └─────┴─────┘
        """
        return super().explode(columns, *more_columns)


class LazyFrame(BaseFrame[LazyFrameT]):
    """Narwhals LazyFrame, backed by a native lazyframe.

    Warning:
        This class is not meant to be instantiated directly - instead use
        [`narwhals.from_native`][] with a native
        object that is a lazy dataframe from one of the supported
        backend (e.g. polars.LazyFrame, dask_expr._collection.DataFrame):
        ```py
        narwhals.from_native(native_lazyframe)
        ```
    """

    @property
    def _compliant(self) -> CompliantLazyFrame[Any, LazyFrameT, Self]:
        return self._compliant_frame

    def _extract_compliant(self, arg: Any) -> Any:
        from narwhals.expr import Expr
        from narwhals.series import Series

        if isinstance(arg, Series):  # pragma: no cover
            msg = "Binary operations between Series and LazyFrame are not supported."
            raise TypeError(msg)
        if isinstance(arg, (Expr, str)):
            if isinstance(arg, Expr):
                if arg._metadata.n_orderable_ops:
                    msg = (
                        "Order-dependent expressions are not supported for use in LazyFrame.\n\n"
                        "Hint: To make the expression valid, use `.over` with `order_by` specified.\n\n"
                        "For example, if you wrote `nw.col('price').cum_sum()` and you have a column\n"
                        "`'date'` which orders your data, then replace:\n\n"
                        "   nw.col('price').cum_sum()\n\n"
                        " with:\n\n"
                        "   nw.col('price').cum_sum().over(order_by='date')\n"
                        "                            ^^^^^^^^^^^^^^^^^^^^^^\n\n"
                        "See https://narwhals-dev.github.io/narwhals/concepts/order_dependence/."
                    )
                    raise InvalidOperationError(msg)
                if arg._metadata.is_filtration:
                    msg = (
                        "Length-changing expressions are not supported for use in LazyFrame, unless\n"
                        "followed by an aggregation.\n\n"
                        "Hints:\n"
                        "- Instead of `lf.select(nw.col('a').head())`, use `lf.select('a').head()\n"
                        "- Instead of `lf.select(nw.col('a').drop_nulls()).select(nw.sum('a'))`,\n"
                        "  use `lf.select(nw.col('a').drop_nulls().sum())\n"
                    )
                    raise InvalidOperationError(msg)
            return self.__narwhals_namespace__().parse_into_expr(arg, str_as_lit=False)
        raise InvalidIntoExprError.from_invalid_type(type(arg))

    @property
    def _dataframe(self) -> type[DataFrame[Any]]:
        return DataFrame

    def __init__(self, df: Any, *, level: Literal["full", "lazy", "interchange"]) -> None:
        self._level = level
        self._compliant_frame: CompliantLazyFrame[Any, LazyFrameT, Self]
        if is_compliant_lazyframe(df):
            self._compliant_frame = df.__narwhals_lazyframe__()
        else:  # pragma: no cover
            msg = f"Expected Polars LazyFrame or an object that implements `__narwhals_lazyframe__`, got: {type(df)}"
            raise AssertionError(msg)

    def __repr__(self) -> str:  # pragma: no cover
        return generate_repr("Narwhals LazyFrame", self.to_native().__repr__())

    def __getitem__(self, item: str | slice) -> NoReturn:
        msg = "Slicing is not supported on LazyFrame"
        raise TypeError(msg)

    def collect(
        self, backend: IntoBackend[Polars | Pandas | Arrow] | None = None, **kwargs: Any
    ) -> DataFrame[Any]:
        r"""Materialize this LazyFrame into a DataFrame.

        As each underlying lazyframe has different arguments to set when materializing
        the lazyframe into a dataframe, we allow to pass them as kwargs (see examples
        below for how to generalize the specification).

        Arguments:
            backend: specifies which eager backend collect to. This will be the underlying
                backend for the resulting Narwhals DataFrame. If None, then the following
                default conversions will be applied

                - `polars.LazyFrame` -> `polars.DataFrame`
                - `dask.DataFrame` -> `pandas.DataFrame`
                - `duckdb.PyRelation` -> `pyarrow.Table`
                - `pyspark.DataFrame` -> `pyarrow.Table`

                `backend` can be specified in various ways

                - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`
                    or `POLARS`.
                - As a string: `"pandas"`, `"pyarrow"` or `"polars"`
                - Directly as a module `pandas`, `pyarrow` or `polars`.
            kwargs: backend specific kwargs to pass along. To know more please check the
                backend specific documentation

                - [polars.LazyFrame.collect](https://docs.pola.rs/api/python/dev/reference/lazyframe/api/polars.LazyFrame.collect.html)
                - [dask.dataframe.DataFrame.compute](https://docs.dask.org/en/stable/generated/dask.dataframe.DataFrame.compute.html)

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 2), (3, 4) df(a, b)")
            >>> lf = nw.from_native(lf_native)
            >>> lf
            ┌──────────────────┐
            |Narwhals LazyFrame|
            |------------------|
            |┌───────┬───────┐ |
            |│   a   │   b   │ |
            |│ int32 │ int32 │ |
            |├───────┼───────┤ |
            |│     1 │     2 │ |
            |│     3 │     4 │ |
            |└───────┴───────┘ |
            └──────────────────┘
            >>> lf.collect()
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  pyarrow.Table   |
            |  a: int32        |
            |  b: int32        |
            |  ----            |
            |  a: [[1,3]]      |
            |  b: [[2,4]]      |
            └──────────────────┘
        """
        collect = self._compliant_frame.collect
        if backend is None:
            return self._dataframe(collect(None, **kwargs), level="full")
        eager_backend = Implementation.from_backend(backend)
        if can_lazyframe_collect(eager_backend):
            return self._dataframe(collect(eager_backend, **kwargs), level="full")
        msg = f"Unsupported `backend` value.\nExpected one of {get_args(_LazyFrameCollectImpl)} or None, got: {eager_backend}."
        raise ValueError(msg)

    def to_native(self) -> LazyFrameT:
        """Convert Narwhals LazyFrame to native one.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 2), (3, 4) df(a, b)")
            >>> nw.from_native(lf_native).to_native()
            ┌───────┬───────┐
            │   a   │   b   │
            │ int32 │ int32 │
            ├───────┼───────┤
            │     1 │     2 │
            │     3 │     4 │
            └───────┴───────┘
            <BLANKLINE>
        """
        return to_native(narwhals_object=self, pass_through=False)

    # inherited
    def pipe(
        self,
        function: Callable[Concatenate[Self, PS], R],
        *args: PS.args,
        **kwargs: PS.kwargs,
    ) -> R:
        """Pipe function call.

        Arguments:
            function: Function to apply.
            args: Positional arguments to pass to function.
            kwargs: Keyword arguments to pass to function.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 2), (3, 4) df(a, b)")
            >>> nw.from_native(lf_native).pipe(lambda x: x.select("a")).to_native()
            ┌───────┐
            │   a   │
            │ int32 │
            ├───────┤
            │     1 │
            │     3 │
            └───────┘
            <BLANKLINE>
        """
        return super().pipe(function, *args, **kwargs)

    def drop_nulls(self, subset: str | list[str] | None = None) -> Self:
        """Drop rows that contain null values.

        Arguments:
            subset: Column name(s) for which null values are considered. If set to None
                (default), use all columns.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, NULL), (3, 4) df(a, b)")
            >>> nw.from_native(lf_native).drop_nulls()
            ┌──────────────────┐
            |Narwhals LazyFrame|
            |------------------|
            |┌───────┬───────┐ |
            |│   a   │   b   │ |
            |│ int32 │ int32 │ |
            |├───────┼───────┤ |
            |│     3 │     4 │ |
            |└───────┴───────┘ |
            └──────────────────┘
        """
        return super().drop_nulls(subset=subset)

    def with_row_index(
        self, name: str = "index", *, order_by: str | Sequence[str]
    ) -> Self:
        """Insert column which enumerates rows.

        Arguments:
            name: The name of the column as a string. The default is "index".
            order_by: Column(s) to order by when computing the row index.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 5), (2, 4) df(a, b)")
            >>> nw.from_native(lf_native).with_row_index(order_by="a").sort("a").collect()
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  pyarrow.Table   |
            |  index: int64    |
            |  a: int32        |
            |  b: int32        |
            |  ----            |
            |  index: [[0,1]]  |
            |  a: [[1,2]]      |
            |  b: [[5,4]]      |
            └──────────────────┘
            >>> nw.from_native(lf_native).with_row_index(order_by="b").sort("a").collect()
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  pyarrow.Table   |
            |  index: int64    |
            |  a: int32        |
            |  b: int32        |
            |  ----            |
            |  index: [[1,0]]  |
            |  a: [[1,2]]      |
            |  b: [[5,4]]      |
            └──────────────────┘
        """
        order_by_ = [order_by] if isinstance(order_by, str) else order_by
        return self._with_compliant(
            self._compliant_frame.with_row_index(name, order_by=order_by_)
        )

    @property
    def schema(self) -> Schema:
        r"""Get an ordered mapping of column names to their data type.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 4.5), (3, 2.) df(a, b)")
            >>> nw.from_native(lf_native).schema  # doctest:+SKIP
            Schema({'a': Int32, 'b': Decimal})
        """
        if self._compliant_frame._version is not Version.V1:
            msg = (
                "Resolving the schema of a LazyFrame is a potentially expensive operation. "
                "Use `LazyFrame.collect_schema()` to get the schema without this warning."
            )
            issue_warning(msg, PerformanceWarning)
        return super().schema

    def collect_schema(self) -> Schema:
        r"""Get an ordered mapping of column names to their data type.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 4.5), (3, 2.) df(a, b)")
            >>> nw.from_native(lf_native).collect_schema()
            Schema({'a': Int32, 'b': Decimal})
        """
        return super().collect_schema()

    @property
    def columns(self) -> list[str]:
        r"""Get column names.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 4.5), (3, 2.) df(a, b)")
            >>> nw.from_native(lf_native).columns
            ['a', 'b']
        """
        return super().columns

    def with_columns(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> Self:
        r"""Add columns to this LazyFrame.

        Added columns will replace existing columns with the same name.

        Arguments:
            *exprs: Column(s) to add, specified as positional arguments.
                     Accepts expression input. Strings are parsed as column names, other
                     non-expression inputs are parsed as literals.

            **named_exprs: Additional columns to add, specified as keyword arguments.
                            The columns will be renamed to the keyword used.

        Note:
            Creating a new LazyFrame using this method does not create a new copy of
            existing data.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 4.5), (3, 2.) df(a, b)")
            >>> nw.from_native(lf_native).with_columns(c=nw.col("a") + 1)
            ┌────────────────────────────────┐
            |       Narwhals LazyFrame       |
            |--------------------------------|
            |┌───────┬──────────────┬───────┐|
            |│   a   │      b       │   c   │|
            |│ int32 │ decimal(2,1) │ int32 │|
            |├───────┼──────────────┼───────┤|
            |│     1 │          4.5 │     2 │|
            |│     3 │          2.0 │     4 │|
            |└───────┴──────────────┴───────┘|
            └────────────────────────────────┘
        """
        if not exprs and not named_exprs:
            msg = "At least one expression must be passed to LazyFrame.with_columns"
            raise ValueError(msg)
        return super().with_columns(*exprs, **named_exprs)

    def select(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> Self:
        r"""Select columns from this LazyFrame.

        Arguments:
            *exprs: Column(s) to select, specified as positional arguments.
                Accepts expression input. Strings are parsed as column names.
            **named_exprs: Additional columns to select, specified as keyword arguments.
                The columns will be renamed to the keyword used.

        Notes:
            If you'd like to select a column whose name isn't a string (for example,
            if you're working with pandas) then you should explicitly use `nw.col` instead
            of just passing the column name. For example, to select a column named
            `0` use `df.select(nw.col(0))`, not `df.select(0)`.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 4.5), (3, 2.) df(a, b)")
            >>> nw.from_native(lf_native).select("a", a_plus_1=nw.col("a") + 1)
            ┌────────────────────┐
            | Narwhals LazyFrame |
            |--------------------|
            |┌───────┬──────────┐|
            |│   a   │ a_plus_1 │|
            |│ int32 │  int32   │|
            |├───────┼──────────┤|
            |│     1 │        2 │|
            |│     3 │        4 │|
            |└───────┴──────────┘|
            └────────────────────┘
        """
        if not exprs and not named_exprs:
            msg = "At least one expression must be passed to LazyFrame.select"
            raise ValueError(msg)
        return super().select(*exprs, **named_exprs)

    def rename(self, mapping: dict[str, str]) -> Self:
        r"""Rename column names.

        Arguments:
            mapping: Key value pairs that map from old name to new name, or a
                      function that takes the old name as input and returns the
                      new name.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 4.5), (3, 2.) df(a, b)")
            >>> nw.from_native(lf_native).rename({"a": "c"})
            ┌────────────────────────┐
            |   Narwhals LazyFrame   |
            |------------------------|
            |┌───────┬──────────────┐|
            |│   c   │      b       │|
            |│ int32 │ decimal(2,1) │|
            |├───────┼──────────────┤|
            |│     1 │          4.5 │|
            |│     3 │          2.0 │|
            |└───────┴──────────────┘|
            └────────────────────────┘
        """
        return super().rename(mapping)

    def head(self, n: int = 5) -> Self:
        r"""Get `n` rows.

        Arguments:
            n: Number of rows to return.

        Examples:
            >>> import dask.dataframe as dd
            >>> import narwhals as nw
            >>> lf_native = dd.from_dict({"a": [1, 2, 3], "b": [4, 5, 6]}, npartitions=1)
            >>> nw.from_native(lf_native).head(2).collect()
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a  b      |
            |     0  1  4      |
            |     1  2  5      |
            └──────────────────┘
        """
        return super().head(n)

    def drop(self, *columns: str | Iterable[str], strict: bool = True) -> Self:
        r"""Remove columns from the LazyFrame.

        Arguments:
            *columns: Names of the columns that should be removed from the dataframe.
            strict: Validate that all column names exist in the schema and throw an
                exception if a column name does not exist in the schema.

        Warning:
            `strict` argument is ignored for `polars<1.0.0`.

            Please consider upgrading to a newer version or pass to eager mode.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 2), (3, 4) df(a, b)")
            >>> nw.from_native(lf_native).drop("a").to_native()
            ┌───────┐
            │   b   │
            │ int32 │
            ├───────┤
            │     2 │
            │     4 │
            └───────┘
            <BLANKLINE>
        """
        return super().drop(*flatten(columns), strict=strict)

    def unique(
        self,
        subset: str | list[str] | None = None,
        *,
        keep: UniqueKeepStrategy = "any",
        order_by: str | Sequence[str] | None = None,
    ) -> Self:
        """Drop duplicate rows from this LazyFrame.

        Arguments:
            subset: Column name(s) to consider when identifying duplicate rows.
                     If set to `None`, use all columns.
            keep: {'any', 'none', 'first', 'last}
                Which of the duplicate rows to keep.

                * 'any': Does not give any guarantee of which row is kept.
                * 'none': Don't keep duplicate rows.
                * 'first': Keep the first row. Requires `order_by` to be specified.
                * 'last': Keep the last row. Requires `order_by` to be specified.
            order_by: Column(s) to order by when computing the row index.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> lf_native = duckdb.sql("SELECT * FROM VALUES (1, 3), (1, 4) df(a, b)")
            >>> nw.from_native(lf_native).unique("a").sort("a", descending=True)
            ┌──────────────────┐
            |Narwhals LazyFrame|
            |------------------|
            |┌───────┬───────┐ |
            |│   a   │   b   │ |
            |│ int32 │ int32 │ |
            |├───────┼───────┤ |
            |│     1 │     3 │ |
            |└───────┴───────┘ |
            └──────────────────┘
        """
        if keep not in {"any", "none", "first", "last"}:
            msg = f"Expected {'any', 'none', 'first', 'last'}, got: {keep}"
            raise ValueError(msg)
        if keep in {"first", "last"} and not order_by:
            msg = (
                "narwhals.LazyFrame makes no assumptions about row order, so only "
                "'first' and 'last' are only supported if `order_by` is passed."
            )
            raise InvalidOperationError(msg)
        if isinstance(subset, str):
            subset = [subset]
        return self._with_compliant(
            self._compliant_frame.unique(subset=subset, keep=keep, order_by=order_by)
        )

    def filter(
        self, *predicates: IntoExpr | Iterable[IntoExpr], **constraints: Any
    ) -> Self:
        r"""Filter the rows in the LazyFrame based on a predicate expression.

        The original order of the remaining rows is preserved.

        Arguments:
            *predicates: Expression(s) that evaluates to a boolean Series.
            **constraints: Column filters; use `name = value` to filter columns by the supplied value.
                Each constraint will behave the same as `nw.col(name).eq(value)`, and will be implicitly
                joined with the other filter conditions using &.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql('''
            ...     SELECT * FROM VALUES
            ...         (1, 6, 'a'),
            ...         (2, 7, 'b'),
            ...         (3, 8, 'c')
            ...     df(foo, bar, ham)
            ... ''')

            Filter on one condition

            >>> nw.from_native(df_native).filter(nw.col("foo") > 1).to_native()
            ┌───────┬───────┬─────────┐
            │  foo  │  bar  │   ham   │
            │ int32 │ int32 │ varchar │
            ├───────┼───────┼─────────┤
            │     2 │     7 │ b       │
            │     3 │     8 │ c       │
            └───────┴───────┴─────────┘
            <BLANKLINE>

            Filter on multiple conditions with implicit `&`

            >>> nw.from_native(df_native).filter(
            ...     nw.col("foo") < 3, nw.col("ham") == "a"
            ... ).to_native()
            ┌───────┬───────┬─────────┐
            │  foo  │  bar  │   ham   │
            │ int32 │ int32 │ varchar │
            ├───────┼───────┼─────────┤
            │     1 │     6 │ a       │
            └───────┴───────┴─────────┘
            <BLANKLINE>

            Filter on multiple conditions with `|`

            >>> nw.from_native(df_native).filter(
            ...     (nw.col("foo") == 1) | (nw.col("ham") == "c")
            ... ).to_native()
            ┌───────┬───────┬─────────┐
            │  foo  │  bar  │   ham   │
            │ int32 │ int32 │ varchar │
            ├───────┼───────┼─────────┤
            │     1 │     6 │ a       │
            │     3 │     8 │ c       │
            └───────┴───────┴─────────┘
            <BLANKLINE>

            Filter using `**kwargs` syntax

            >>> nw.from_native(df_native).filter(foo=2, ham="b").to_native()
            ┌───────┬───────┬─────────┐
            │  foo  │  bar  │   ham   │
            │ int32 │ int32 │ varchar │
            ├───────┼───────┼─────────┤
            │     2 │     7 │ b       │
            └───────┴───────┴─────────┘
            <BLANKLINE>
        """
        predicates_ = tuple(tuple(p) if is_iterator(p) else p for p in predicates)
        if predicates_contains_list_of_bool(predicates_):
            msg = "`LazyFrame.filter` is not supported with Python boolean masks - use expressions instead."
            raise TypeError(msg)

        return super().filter(*predicates_, **constraints)

    def sink_parquet(self, file: str | Path | BytesIO) -> None:
        """Write LazyFrame to Parquet file.

        This may allow larger-than-RAM datasets to be written to disk.

        Arguments:
            file: String, path object or file-like object to which the dataframe will be
                written.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.LazyFrame({"foo": [1, 2], "bar": [6.0, 7.0]})
            >>> df = nw.from_native(df_native)
            >>> df.sink_parquet("out.parquet")  # doctest:+SKIP
        """
        self._compliant_frame.sink_parquet(file)

    @overload
    def group_by(
        self, *keys: IntoExpr | Iterable[IntoExpr], drop_null_keys: Literal[False] = ...
    ) -> LazyGroupBy[Self]: ...

    @overload
    def group_by(
        self, *keys: str | Iterable[str], drop_null_keys: Literal[True]
    ) -> LazyGroupBy[Self]: ...

    def group_by(
        self, *keys: IntoExpr | Iterable[IntoExpr], drop_null_keys: bool = False
    ) -> LazyGroupBy[Self]:
        r"""Start a group by operation.

        Arguments:
            *keys: Column(s) to group by. Accepts expression input. Strings are parsed as
                column names.
            drop_null_keys: if True, then groups where any key is null won't be
                included in the result.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES (1, 'a'), (2, 'b'), (3, 'a') df(a, b)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.group_by("b").agg(nw.col("a").sum()).sort("b").to_native()
            ┌─────────┬────────┐
            │    b    │   a    │
            │ varchar │ int128 │
            ├─────────┼────────┤
            │ a       │      4 │
            │ b       │      2 │
            └─────────┴────────┘
            <BLANKLINE>

            Expressions are also accepted.

            >>> df.group_by(nw.col("b").str.len_chars()).agg(
            ...     nw.col("a").sum()
            ... ).to_native()
            ┌───────┬────────┐
            │   b   │   a    │
            │ int64 │ int128 │
            ├───────┼────────┤
            │     1 │      6 │
            └───────┴────────┘
            <BLANKLINE>
        """
        from narwhals.group_by import LazyGroupBy

        flat_keys = flatten(keys)

        if all(isinstance(key, str) for key in flat_keys):
            return LazyGroupBy(self, flat_keys, drop_null_keys=drop_null_keys)

        from narwhals import col
        from narwhals.expr import Expr

        key_is_expr = tuple(isinstance(k, Expr) for k in flat_keys)

        if drop_null_keys and any(key_is_expr):
            msg = "drop_null_keys cannot be True when keys contains Expr"
            raise NotImplementedError(msg)

        _keys = [
            k if is_expr else col(k) for k, is_expr in zip_strict(flat_keys, key_is_expr)
        ]
        expr_flat_keys, _kinds = self._flatten_and_extract(*_keys)
        check_expressions_preserve_length(*_keys, function_name="LazyFrame.group_by")

        return LazyGroupBy(self, expr_flat_keys, drop_null_keys=drop_null_keys)

    def sort(
        self,
        by: str | Iterable[str],
        *more_by: str,
        descending: bool | Sequence[bool] = False,
        nulls_last: bool = False,
    ) -> Self:
        r"""Sort the LazyFrame by the given columns.

        Arguments:
            by: Column(s) names to sort by.
            *more_by: Additional columns to sort by, specified as positional arguments.
            descending: Sort in descending order. When sorting by multiple columns, can be
                specified per column by passing a sequence of booleans.
            nulls_last: Place null values last; can specify a single boolean applying to
                all columns or a sequence of booleans for per-column control.

        Warning:
            Unlike Polars, it is not possible to specify a sequence of booleans for
            `nulls_last` in order to control per-column behaviour. Instead a single
            boolean is applied for all `by` columns.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES (1, 6.0, 'a'), (2, 5.0, 'c'), (NULL, 4.0, 'b') df(a, b, c)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.sort("a")
            ┌──────────────────────────────────┐
            |        Narwhals LazyFrame        |
            |----------------------------------|
            |┌───────┬──────────────┬─────────┐|
            |│   a   │      b       │    c    │|
            |│ int32 │ decimal(2,1) │ varchar │|
            |├───────┼──────────────┼─────────┤|
            |│  NULL │          4.0 │ b       │|
            |│     1 │          6.0 │ a       │|
            |│     2 │          5.0 │ c       │|
            |└───────┴──────────────┴─────────┘|
            └──────────────────────────────────┘
        """
        return super().sort(by, *more_by, descending=descending, nulls_last=nulls_last)

    def top_k(
        self, k: int, *, by: str | Iterable[str], reverse: bool | Sequence[bool] = False
    ) -> Self:
        r"""Return the `k` largest rows.

        Non-null elements are always preferred over null elements,
        regardless of the value of reverse. The output is not guaranteed
        to be in any particular order, sort the outputs afterwards if you wish the output to be sorted.

        Arguments:
            k: Number of rows to return.
            by: Column(s) used to determine the top rows. Accepts expression input. Strings are parsed as column names.
            reverse: Consider the k smallest elements of the by column(s) (instead of the k largest).
                This can be specified per column by passing a sequence of booleans.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES ('a', 2), ('b', 1), ('a', 1), ('b', 3), (NULL, 2), ('c', 1) df(a, b)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.top_k(4, by=["b", "a"])
            ┌───────────────────┐
            |Narwhals LazyFrame |
            |-------------------|
            |┌─────────┬───────┐|
            |│    a    │   b   │|
            |│ varchar │ int32 │|
            |├─────────┼───────┤|
            |│ b       │     3 │|
            |│ a       │     2 │|
            |│ NULL    │     2 │|
            |│ c       │     1 │|
            |└─────────┴───────┘|
            └───────────────────┘
        """
        return super().top_k(k, by=by, reverse=reverse)

    def join(
        self,
        other: Self,
        on: str | list[str] | None = None,
        how: JoinStrategy = "inner",
        *,
        left_on: str | list[str] | None = None,
        right_on: str | list[str] | None = None,
        suffix: str = "_right",
    ) -> Self:
        r"""Add a join operation to the Logical Plan.

        Arguments:
            other: Lazy DataFrame to join with.
            on: Name(s) of the join columns in both DataFrames. If set, `left_on` and
                `right_on` should be None.
            how: Join strategy.

                  * *inner*: Returns rows that have matching values in both tables.
                  * *left*: Returns all rows from the left table, and the matched rows from the right table.
                  * *full*: Returns all rows in both dataframes, with the suffix appended to the right join keys.
                  * *cross*: Returns the Cartesian product of rows from both tables.
                  * *semi*: Filter rows that have a match in the right table.
                  * *anti*: Filter rows that do not have a match in the right table.
            left_on: Join column of the left DataFrame.
            right_on: Join column of the right DataFrame.
            suffix: Suffix to append to columns with a duplicate name.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native1 = duckdb.sql(
            ...     "SELECT * FROM VALUES (1, 'a'), (2, 'b') df(a, b)"
            ... )
            >>> df_native2 = duckdb.sql(
            ...     "SELECT * FROM VALUES (1, 'x'), (3, 'y') df(a, c)"
            ... )
            >>> df1 = nw.from_native(df_native1)
            >>> df2 = nw.from_native(df_native2)
            >>> df1.join(df2, on="a")
            ┌─────────────────────────────┐
            |     Narwhals LazyFrame      |
            |-----------------------------|
            |┌───────┬─────────┬─────────┐|
            |│   a   │    b    │    c    │|
            |│ int32 │ varchar │ varchar │|
            |├───────┼─────────┼─────────┤|
            |│     1 │ a       │ x       │|
            |└───────┴─────────┴─────────┘|
            └─────────────────────────────┘
        """
        return super().join(
            other, how=how, left_on=left_on, right_on=right_on, on=on, suffix=suffix
        )

    def join_asof(
        self,
        other: Self,
        *,
        left_on: str | None = None,
        right_on: str | None = None,
        on: str | None = None,
        by_left: str | list[str] | None = None,
        by_right: str | list[str] | None = None,
        by: str | list[str] | None = None,
        strategy: AsofJoinStrategy = "backward",
        suffix: str = "_right",
    ) -> Self:
        """Perform an asof join.

        This is similar to a left-join except that we match on nearest key rather than equal keys.

        For Polars, both DataFrames must be sorted by the `on` key (within each `by` group
        if specified).

        Arguments:
            other: DataFrame to join with.
            left_on: Name(s) of the left join column(s).
            right_on: Name(s) of the right join column(s).
            on: Join column of both DataFrames. If set, left_on and right_on should be None.
            by_left: join on these columns before doing asof join
            by_right: join on these columns before doing asof join
            by: join on these columns before doing asof join
            strategy: Join strategy. The default is "backward".

                  * *backward*: selects the last row in the right DataFrame whose "on" key is less than or equal to the left's key.
                  * *forward*: selects the first row in the right DataFrame whose "on" key is greater than or equal to the left's key.
                  * *nearest*: search selects the last row in the right DataFrame whose value is nearest to the left's key.

            suffix: Suffix to append to columns with a duplicate name.

        Examples:
            >>> from datetime import datetime
            >>> import polars as pl
            >>> import narwhals as nw
            >>> data_gdp = {
            ...     "datetime": [
            ...         datetime(2016, 1, 1),
            ...         datetime(2017, 1, 1),
            ...         datetime(2018, 1, 1),
            ...         datetime(2019, 1, 1),
            ...         datetime(2020, 1, 1),
            ...     ],
            ...     "gdp": [4164, 4411, 4566, 4696, 4827],
            ... }
            >>> data_population = {
            ...     "datetime": [
            ...         datetime(2016, 3, 1),
            ...         datetime(2018, 8, 1),
            ...         datetime(2019, 1, 1),
            ...     ],
            ...     "population": [82.19, 82.66, 83.12],
            ... }
            >>> gdp_native = pl.DataFrame(data_gdp)
            >>> population_native = pl.DataFrame(data_population)
            >>> gdp = nw.from_native(gdp_native)
            >>> population = nw.from_native(population_native)
            >>> population.join_asof(gdp, on="datetime", strategy="backward").to_native()
            shape: (3, 3)
            ┌─────────────────────┬────────────┬──────┐
            │ datetime            ┆ population ┆ gdp  │
            │ ---                 ┆ ---        ┆ ---  │
            │ datetime[μs]        ┆ f64        ┆ i64  │
            ╞═════════════════════╪════════════╪══════╡
            │ 2016-03-01 00:00:00 ┆ 82.19      ┆ 4164 │
            │ 2018-08-01 00:00:00 ┆ 82.66      ┆ 4566 │
            │ 2019-01-01 00:00:00 ┆ 83.12      ┆ 4696 │
            └─────────────────────┴────────────┴──────┘
        """
        return super().join_asof(
            other,
            left_on=left_on,
            right_on=right_on,
            on=on,
            by_left=by_left,
            by_right=by_right,
            by=by,
            strategy=strategy,
            suffix=suffix,
        )

    def lazy(self) -> Self:
        """Restrict available API methods to lazy-only ones.

        This is a no-op, and exists only for compatibility with `DataFrame.lazy`.
        """
        return self

    def unpivot(
        self,
        on: str | list[str] | None = None,
        *,
        index: str | list[str] | None = None,
        variable_name: str = "variable",
        value_name: str = "value",
    ) -> Self:
        r"""Unpivot a DataFrame from wide to long format.

        Optionally leaves identifiers set.

        This function is useful to massage a DataFrame into a format where one or more
        columns are identifier variables (index) while all other columns, considered
        measured variables (on), are "unpivoted" to the row axis leaving just
        two non-identifier columns, 'variable' and 'value'.

        Arguments:
            on: Column(s) to use as values variables; if `on` is empty all columns that
                are not in `index` will be used.
            index: Column(s) to use as identifier variables.
            variable_name: Name to give to the `variable` column. Defaults to "variable".
            value_name: Name to give to the `value` column. Defaults to "value".

        Notes:
            If you're coming from pandas, this is similar to `pandas.DataFrame.melt`,
            but with `index` replacing `id_vars` and `on` replacing `value_vars`.
            In other frameworks, you might know this operation as `pivot_longer`.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES ('x', 1, 2), ('y', 3, 4), ('z', 5, 6) df(a, b, c)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.unpivot(on=["b", "c"], index="a").sort("a", "variable").to_native()
            ┌─────────┬──────────┬───────┐
            │    a    │ variable │ value │
            │ varchar │ varchar  │ int32 │
            ├─────────┼──────────┼───────┤
            │ x       │ b        │     1 │
            │ x       │ c        │     2 │
            │ y       │ b        │     3 │
            │ y       │ c        │     4 │
            │ z       │ b        │     5 │
            │ z       │ c        │     6 │
            └─────────┴──────────┴───────┘
            <BLANKLINE>
        """
        return super().unpivot(
            on=on, index=index, variable_name=variable_name, value_name=value_name
        )

    def explode(self, columns: str | Sequence[str], *more_columns: str) -> Self:
        """Explode the dataframe to long format by exploding the given columns.

        Notes:
            It is possible to explode multiple columns only if these columns have
            matching element counts.

        Arguments:
            columns: Column names. The underlying columns being exploded must be of the `List` data type.
            *more_columns: Additional names of columns to explode, specified as positional arguments.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES ('x', [1, 2]), ('y', [3, 4]), ('z', [5, 6]) df(a, b)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.explode("b").to_native()
            ┌─────────┬───────┐
            │    a    │   b   │
            │ varchar │ int32 │
            ├─────────┼───────┤
            │ x       │     1 │
            │ x       │     2 │
            │ y       │     3 │
            │ y       │     4 │
            │ z       │     5 │
            │ z       │     6 │
            └─────────┴───────┘
            <BLANKLINE>
        """
        return super().explode(columns, *more_columns)
