from __future__ import annotations

import contextlib
import math
import operator
import sys
import warnings
from collections.abc import Collection, Mapping, Sequence
from datetime import timedelta
from functools import reduce
from io import BytesIO, StringIO
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    NoReturn,
    TypeVar,
)

import polars._reexport as pl
from polars import functions as F
from polars._dependencies import _check_for_numpy
from polars._dependencies import numpy as np
from polars._utils.convert import negate_duration_string, parse_as_duration_string
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    deprecated,
    issue_deprecation_warning,
)
from polars._utils.parse import (
    parse_into_expression,
    parse_into_list_of_expressions,
    parse_predicates_constraints_into_expression,
)
from polars._utils.unstable import issue_unstable_warning, unstable
from polars._utils.various import (
    BUILDING_SPHINX_DOCS,
    extend_bool,
    find_stacklevel,
    no_default,
    normalize_filepath,
    sphinx_accessor,
    warn_null_comparison,
)
from polars._utils.wrap import wrap_expr, wrap_s
from polars.datatypes import (
    Int64,
    parse_into_datatype_expr,
)
from polars.exceptions import (
    CustomUFuncWarning,
    OutOfBoundsError,
    PolarsInefficientMapWarning,
)
from polars.expr.array import ExprArrayNameSpace
from polars.expr.binary import ExprBinaryNameSpace
from polars.expr.categorical import ExprCatNameSpace
from polars.expr.datetime import ExprDateTimeNameSpace
from polars.expr.list import ExprListNameSpace
from polars.expr.meta import ExprMetaNameSpace
from polars.expr.name import ExprNameNameSpace
from polars.expr.string import ExprStringNameSpace
from polars.expr.struct import ExprStructNameSpace
from polars.meta import thread_pool_size

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import arg_where as py_arg_where

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyExpr

if TYPE_CHECKING:
    with contextlib.suppress(ImportError):  # Module not available when building docs
        from polars._plr import PySeries

    with contextlib.suppress(ImportError):  # Module not available when building docs
        import polars._plr as plr

    from collections.abc import Iterable
    from io import IOBase

    from polars import DataFrame, LazyFrame, Series
    from polars._typing import (
        ClosedInterval,
        FillNullStrategy,
        InterpolationMethod,
        IntoExpr,
        IntoExprColumn,
        MapElementsStrategy,
        NullBehavior,
        NumericLiteral,
        PolarsDataType,
        QuantileMethod,
        RankMethod,
        RoundMode,
        SchemaDict,
        SearchSortedSide,
        SerializationFormat,
        TemporalLiteral,
        WindowMappingStrategy,
    )
    from polars._utils.various import NoDefault

    if sys.version_info >= (3, 11):
        from typing import Concatenate, ParamSpec
    else:
        from typing_extensions import Concatenate, ParamSpec

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004

    T = TypeVar("T")
    P = ParamSpec("P")

elif BUILDING_SPHINX_DOCS:
    # note: we assign this way to work around an autocomplete issue in ipython/jedi
    # (ref: https://github.com/davidhalter/jedi/issues/2057)
    current_module = sys.modules[__name__]
    current_module.property = sphinx_accessor


class Expr:
    """Expressions that can be used in various contexts."""

    # NOTE: This `= None` is needed to generate the docs with sphinx_accessor.
    _pyexpr: PyExpr = None  # type: ignore[assignment]
    _accessors: ClassVar[set[str]] = {
        "arr",
        "bin",
        "cat",
        "dt",
        "list",
        "meta",
        "name",
        "str",
        "struct",
    }

    @classmethod
    def _from_pyexpr(cls, pyexpr: PyExpr) -> Expr:
        expr = cls.__new__(cls)
        expr._pyexpr = pyexpr
        return expr

    def _repr_html_(self) -> str:
        return self._pyexpr.to_str()

    def __repr__(self) -> str:
        if len(expr_str := self._pyexpr.to_str()) > 30:
            expr_str = f"{expr_str[:30]}…"
        return f"<{self.__class__.__name__} [{expr_str!r}] at 0x{id(self):X}>"

    def __str__(self) -> str:
        return self._pyexpr.to_str()

    def __bool__(self) -> NoReturn:
        msg = (
            "the truth value of an Expr is ambiguous"
            "\n\n"
            "You probably got here by using a Python standard library function instead "
            "of the native expressions API.\n"
            "Here are some things you might want to try:\n"
            "- instead of `pl.col('a') and pl.col('b')`, use `pl.col('a') & pl.col('b')`\n"
            "- instead of `pl.col('a') in [y, z]`, use `pl.col('a').is_in([y, z])`\n"
            "- instead of `max(pl.col('a'), pl.col('b'))`, use `pl.max_horizontal(pl.col('a'), pl.col('b'))`\n"
        )
        raise TypeError(msg)

    def __abs__(self) -> Expr:
        return self.abs()

    # operators
    def __add__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr + other_pyexpr)

    def __radd__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(other_pyexpr + self._pyexpr)

    def __and__(self, other: IntoExprColumn | int | bool) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr.and_(other_pyexpr))

    def __rand__(self, other: IntoExprColumn | int | bool) -> Expr:
        other_expr = parse_into_expression(other)
        return wrap_expr(other_expr.and_(self._pyexpr))

    def __eq__(self, other: IntoExpr) -> Expr:  # type: ignore[override]
        warn_null_comparison(other)
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.eq(other_pyexpr))

    def __floordiv__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr // other_pyexpr)

    def __rfloordiv__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(other_pyexpr // self._pyexpr)

    def __ge__(self, other: IntoExpr) -> Expr:
        warn_null_comparison(other)
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.gt_eq(other_pyexpr))

    def __gt__(self, other: IntoExpr) -> Expr:
        warn_null_comparison(other)
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.gt(other_pyexpr))

    def __invert__(self) -> Expr:
        return self.not_()

    def __le__(self, other: IntoExpr) -> Expr:
        warn_null_comparison(other)
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.lt_eq(other_pyexpr))

    def __lt__(self, other: IntoExpr) -> Expr:
        warn_null_comparison(other)
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.lt(other_pyexpr))

    def __mod__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr % other_pyexpr)

    def __rmod__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(other_pyexpr % self._pyexpr)

    def __mul__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr * other_pyexpr)

    def __rmul__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(other_pyexpr * self._pyexpr)

    def __ne__(self, other: IntoExpr) -> Expr:  # type: ignore[override]
        warn_null_comparison(other)
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.neq(other_pyexpr))

    def __neg__(self) -> Expr:
        return wrap_expr(-self._pyexpr)

    def __or__(self, other: IntoExprColumn | int | bool) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr.or_(other_pyexpr))

    def __ror__(self, other: IntoExprColumn | int | bool) -> Expr:
        other_expr = parse_into_expression(other)
        return wrap_expr(other_expr.or_(self._pyexpr))

    def __pos__(self) -> Expr:
        return self

    def __pow__(self, exponent: IntoExprColumn | int | float) -> Expr:
        exponent_pyexpr = parse_into_expression(exponent)
        return wrap_expr(self._pyexpr.pow(exponent_pyexpr))

    def __rpow__(self, base: IntoExprColumn | int | float) -> Expr:
        base_pyexpr = parse_into_expression(base)
        return wrap_expr(base_pyexpr) ** self

    def __sub__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr - other_pyexpr)

    def __rsub__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(other_pyexpr - self._pyexpr)

    def __truediv__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr / other_pyexpr)

    def __rtruediv__(self, other: IntoExpr) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(other_pyexpr / self._pyexpr)

    def __xor__(self, other: IntoExprColumn | int | bool) -> Expr:
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr.xor_(other_pyexpr))

    def __rxor__(self, other: IntoExprColumn | int | bool) -> Expr:
        other_expr = parse_into_expression(other)
        return wrap_expr(other_expr.xor_(self._pyexpr))

    def __getstate__(self) -> bytes:
        return self._pyexpr.__getstate__()

    def __setstate__(self, state: bytes) -> None:
        self._pyexpr = F.lit(0)._pyexpr  # Initialize with a dummy
        self._pyexpr.__setstate__(state)

    def __array_ufunc__(
        self, ufunc: Callable[..., Any], method: str, *inputs: Any, **kwargs: Any
    ) -> Expr:
        """Numpy universal functions."""
        if method != "__call__":
            msg = f"Only call is implemented not {method}"
            raise NotImplementedError(msg)
        # Numpy/Scipy ufuncs have signature None but numba signatures always exists.
        is_custom_ufunc = getattr(ufunc, "signature") is not None  # noqa: B009
        if is_custom_ufunc is True:
            msg = (
                "Native numpy ufuncs are dispatched using `map_batches(ufunc, is_elementwise=True)` which "
                "is safe for native Numpy and Scipy ufuncs but custom ufuncs in a group_by "
                "context won't be properly grouped. Custom ufuncs are dispatched with is_elementwise=False. "
                f"If {ufunc.__name__} needs elementwise then please use map_batches directly."
            )
            warnings.warn(
                msg,
                CustomUFuncWarning,
                stacklevel=find_stacklevel(),
            )
        if len(inputs) == 1 and len(kwargs) == 0:
            # if there is only 1 input then it must be an Expr for this func to
            # have been called. If there are no kwargs then call map_batches
            # directly on the ufunc
            if not isinstance(inputs[0], Expr):
                msg = "Input must be expression."
                raise OutOfBoundsError(msg)
            return inputs[0].map_batches(ufunc, is_elementwise=not is_custom_ufunc)
        num_expr = sum(isinstance(inp, Expr) for inp in inputs)
        exprs = [
            (inp, True, i) if isinstance(inp, Expr) else (inp, False, i)
            for i, inp in enumerate(inputs)
        ]

        if num_expr == 1:
            root_expr = next(expr[0] for expr in exprs if expr[1])
        else:
            # We rename all but the first expression in case someone did e.g.
            # np.divide(pl.col("a"), pl.col("a")); we'll be creating a struct
            # below, and structs can't have duplicate names.
            first_renameable_expr = True
            actual_exprs = []
            for inp, is_actual_expr, index in exprs:
                if is_actual_expr:
                    if first_renameable_expr:
                        first_renameable_expr = False
                    else:
                        inp = inp.alias(f"argument_{index}")
                    actual_exprs.append(inp)
            root_expr = F.struct(actual_exprs)

        def function(s: Series) -> Series:  # pragma: no cover
            args: list[Any] = []
            for i, expr in enumerate(exprs):
                if expr[1] and num_expr > 1:
                    args.append(s.struct[i])
                elif expr[1]:
                    args.append(s)
                else:
                    args.append(expr[0])
            return ufunc(*args, **kwargs)

        return root_expr.map_batches(function, is_elementwise=not is_custom_ufunc)

    @classmethod
    def deserialize(
        cls,
        source: str | Path | IOBase | bytes,
        *,
        format: SerializationFormat = "binary",
    ) -> Expr:
        """
        Read a serialized expression from a file.

        Parameters
        ----------
        source
            Path to a file or a file-like object (by file-like object, we refer to
            objects that have a `read()` method, such as a file handler (e.g.
            via builtin `open` function) or `BytesIO`).
        format
            The format with which the Expr was serialized. Options:

            - `"binary"`: Deserialize from binary format (bytes). This is the default.
            - `"json"`: Deserialize from JSON format (string).

        Warnings
        --------
        This function uses :mod:`pickle` if the logical plan contains Python UDFs,
        and as such inherits the security implications. Deserializing can execute
        arbitrary code, so it should only be attempted on trusted data.

        See Also
        --------
        Expr.meta.serialize

        Notes
        -----
        Serialization is not stable across Polars versions: a LazyFrame serialized
        in one Polars version may not be deserializable in another Polars version.

        Examples
        --------
        >>> import io
        >>> expr = pl.col("foo").sum().over("bar")
        >>> bytes = expr.meta.serialize()
        >>> pl.Expr.deserialize(io.BytesIO(bytes))
        <Expr ['col("foo").sum().over([col("ba…'] at ...>
        """
        if isinstance(source, StringIO):
            source = BytesIO(source.getvalue().encode())
        elif isinstance(source, (str, Path)):
            source = normalize_filepath(source)
        elif isinstance(source, bytes):
            source = BytesIO(source)

        if format == "binary":
            deserializer = PyExpr.deserialize_binary
        elif format == "json":
            deserializer = PyExpr.deserialize_json
        else:
            msg = f"`format` must be one of {{'binary', 'json'}}, got {format!r}"
            raise ValueError(msg)

        return cls._from_pyexpr(deserializer(source))

    def to_physical(self) -> Expr:
        """
        Cast to physical representation of the logical dtype.

        - :func:`polars.datatypes.Date` -> :func:`polars.datatypes.Int32`
        - :func:`polars.datatypes.Datetime` -> :func:`polars.datatypes.Int64`
        - :func:`polars.datatypes.Time` -> :func:`polars.datatypes.Int64`
        - :func:`polars.datatypes.Duration` -> :func:`polars.datatypes.Int64`
        - :func:`polars.datatypes.Categorical` -> :func:`polars.datatypes.UInt32`
        - `List(inner)` -> `List(physical of inner)`
        - `Array(inner)` -> `Struct(physical of inner)`
        - `Struct(fields)` -> `Array(physical of fields)`

        Other data types will be left unchanged.

        Warnings
        --------
        The physical representations are an implementation detail
        and not guaranteed to be stable.

        Examples
        --------
        Replicating the pandas
        `pd.factorize
        <https://pandas.pydata.org/docs/reference/api/pandas.factorize.html>`_
        function.

        >>> pl.DataFrame({"vals": ["a", "x", None, "a"]}).with_columns(
        ...     pl.col("vals").cast(pl.Categorical),
        ...     pl.col("vals")
        ...     .cast(pl.Categorical)
        ...     .to_physical()
        ...     .alias("vals_physical"),
        ... )
        shape: (4, 2)
        ┌──────┬───────────────┐
        │ vals ┆ vals_physical │
        │ ---  ┆ ---           │
        │ cat  ┆ u32           │
        ╞══════╪═══════════════╡
        │ a    ┆ 0             │
        │ x    ┆ 1             │
        │ null ┆ null          │
        │ a    ┆ 0             │
        └──────┴───────────────┘
        """
        return wrap_expr(self._pyexpr.to_physical())

    def any(self, *, ignore_nulls: bool = True) -> Expr:
        """
        Return whether any of the values in the column are `True`.

        Only works on columns of data type :class:`Boolean`.

        Parameters
        ----------
        ignore_nulls
            * If set to `True` (default), null values are ignored. If there
              are no non-null values, the output is `False`.
            * If set to `False`, `Kleene logic`_ is used to deal with nulls:
              if the column contains any null values and no `True` values,
              the output is null.

            .. _Kleene logic: https://en.wikipedia.org/wiki/Three-valued_logic

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [True, False],
        ...         "b": [False, False],
        ...         "c": [None, False],
        ...     }
        ... )
        >>> df.select(pl.col("*").any())
        shape: (1, 3)
        ┌──────┬───────┬───────┐
        │ a    ┆ b     ┆ c     │
        │ ---  ┆ ---   ┆ ---   │
        │ bool ┆ bool  ┆ bool  │
        ╞══════╪═══════╪═══════╡
        │ true ┆ false ┆ false │
        └──────┴───────┴───────┘

        Enable Kleene logic by setting `ignore_nulls=False`.

        >>> df.select(pl.col("*").any(ignore_nulls=False))
        shape: (1, 3)
        ┌──────┬───────┬──────┐
        │ a    ┆ b     ┆ c    │
        │ ---  ┆ ---   ┆ ---  │
        │ bool ┆ bool  ┆ bool │
        ╞══════╪═══════╪══════╡
        │ true ┆ false ┆ null │
        └──────┴───────┴──────┘
        """
        return wrap_expr(self._pyexpr.any(ignore_nulls))

    def all(self, *, ignore_nulls: bool = True) -> Expr:
        """
        Return whether all values in the column are `True`.

        Only works on columns of data type :class:`Boolean`.

        .. note::
            This method is not to be confused with the function :func:`polars.all`,
            which can be used to select all columns.

        Parameters
        ----------
        ignore_nulls
            * If set to `True` (default), null values are ignored. If there
              are no non-null values, the output is `True`.
            * If set to `False`, `Kleene logic`_ is used to deal with nulls:
              if the column contains any null values and no `False` values,
              the output is null.

            .. _Kleene logic: https://en.wikipedia.org/wiki/Three-valued_logic

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [True, True],
        ...         "b": [False, True],
        ...         "c": [None, True],
        ...     }
        ... )
        >>> df.select(pl.col("*").all())
        shape: (1, 3)
        ┌──────┬───────┬──────┐
        │ a    ┆ b     ┆ c    │
        │ ---  ┆ ---   ┆ ---  │
        │ bool ┆ bool  ┆ bool │
        ╞══════╪═══════╪══════╡
        │ true ┆ false ┆ true │
        └──────┴───────┴──────┘

        Enable Kleene logic by setting `ignore_nulls=False`.

        >>> df.select(pl.col("*").all(ignore_nulls=False))
        shape: (1, 3)
        ┌──────┬───────┬──────┐
        │ a    ┆ b     ┆ c    │
        │ ---  ┆ ---   ┆ ---  │
        │ bool ┆ bool  ┆ bool │
        ╞══════╪═══════╪══════╡
        │ true ┆ false ┆ null │
        └──────┴───────┴──────┘
        """
        return wrap_expr(self._pyexpr.all(ignore_nulls))

    def arg_true(self) -> Expr:
        """
        Return indices where expression evaluates `True`.

        .. warning::
            Modifies number of rows returned, so will fail in combination with other
            expressions. Use as only expression in `select` / `with_columns`.

        See Also
        --------
        Series.arg_true : Return indices where Series is True
        polars.arg_where

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2, 1]})
        >>> df.select((pl.col("a") == 1).arg_true())
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 0   │
        │ 1   │
        │ 3   │
        └─────┘
        """
        return wrap_expr(py_arg_where(self._pyexpr))

    def sqrt(self) -> Expr:
        """
        Compute the square root of the elements.

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1.0, 2.0, 4.0]})
        >>> df.select(pl.col("values").sqrt())
        shape: (3, 1)
        ┌──────────┐
        │ values   │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.0      │
        │ 1.414214 │
        │ 2.0      │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.sqrt())

    def cbrt(self) -> Expr:
        """
        Compute the cube root of the elements.

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1.0, 2.0, 4.0]})
        >>> df.select(pl.col("values").cbrt())
        shape: (3, 1)
        ┌──────────┐
        │ values   │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.0      │
        │ 1.259921 │
        │ 1.587401 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.cbrt())

    def log10(self) -> Expr:
        """
        Compute the base 10 logarithm of the input array, element-wise.

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1.0, 2.0, 4.0]})
        >>> df.select(pl.col("values").log10())
        shape: (3, 1)
        ┌─────────┐
        │ values  │
        │ ---     │
        │ f64     │
        ╞═════════╡
        │ 0.0     │
        │ 0.30103 │
        │ 0.60206 │
        └─────────┘
        """
        return self.log(10.0)

    def exp(self) -> Expr:
        """
        Compute the exponential, element-wise.

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1.0, 2.0, 4.0]})
        >>> df.select(pl.col("values").exp())
        shape: (3, 1)
        ┌──────────┐
        │ values   │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 2.718282 │
        │ 7.389056 │
        │ 54.59815 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.exp())

    def alias(self, name: str) -> Expr:
        """
        Rename the expression.

        Parameters
        ----------
        name
            The new name.

        See Also
        --------
        name.map
        name.prefix
        name.suffix

        Examples
        --------
        Rename an expression to avoid overwriting an existing column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": ["x", "y", "z"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("a") + 10,
        ...     pl.col("b").str.to_uppercase().alias("c"),
        ... )
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ str ┆ str │
        ╞═════╪═════╪═════╡
        │ 11  ┆ x   ┆ X   │
        │ 12  ┆ y   ┆ Y   │
        │ 13  ┆ z   ┆ Z   │
        └─────┴─────┴─────┘

        Overwrite the default name of literal columns to prevent errors due to duplicate
        column names.

        >>> df.with_columns(
        ...     pl.lit(True).alias("c"),
        ...     pl.lit(4.0).alias("d"),
        ... )
        shape: (3, 4)
        ┌─────┬─────┬──────┬─────┐
        │ a   ┆ b   ┆ c    ┆ d   │
        │ --- ┆ --- ┆ ---  ┆ --- │
        │ i64 ┆ str ┆ bool ┆ f64 │
        ╞═════╪═════╪══════╪═════╡
        │ 1   ┆ x   ┆ true ┆ 4.0 │
        │ 2   ┆ y   ┆ true ┆ 4.0 │
        │ 3   ┆ z   ┆ true ┆ 4.0 │
        └─────┴─────┴──────┴─────┘
        """
        return wrap_expr(self._pyexpr.alias(name))

    def exclude(
        self,
        columns: str | PolarsDataType | Collection[str] | Collection[PolarsDataType],
        *more_columns: str | PolarsDataType,
    ) -> Expr:
        """
        Exclude columns from a multi-column expression.

        Only works after a wildcard or regex column selection, and you cannot provide
        both string column names *and* dtypes (you may prefer to use selectors instead).

        Parameters
        ----------
        columns
            The name or datatype of the column(s) to exclude. Accepts regular expression
            input. Regular expressions should start with `^` and end with `$`.
        *more_columns
            Additional names or datatypes of columns to exclude, specified as positional
            arguments.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "aa": [1, 2, 3],
        ...         "ba": ["a", "b", None],
        ...         "cc": [None, 2.5, 1.5],
        ...     }
        ... )
        >>> df
        shape: (3, 3)
        ┌─────┬──────┬──────┐
        │ aa  ┆ ba   ┆ cc   │
        │ --- ┆ ---  ┆ ---  │
        │ i64 ┆ str  ┆ f64  │
        ╞═════╪══════╪══════╡
        │ 1   ┆ a    ┆ null │
        │ 2   ┆ b    ┆ 2.5  │
        │ 3   ┆ null ┆ 1.5  │
        └─────┴──────┴──────┘

        Exclude by column name(s):

        >>> df.select(pl.all().exclude("ba"))
        shape: (3, 2)
        ┌─────┬──────┐
        │ aa  ┆ cc   │
        │ --- ┆ ---  │
        │ i64 ┆ f64  │
        ╞═════╪══════╡
        │ 1   ┆ null │
        │ 2   ┆ 2.5  │
        │ 3   ┆ 1.5  │
        └─────┴──────┘

        Exclude by regex, e.g. removing all columns whose names end with the letter "a":

        >>> df.select(pl.all().exclude("^.*a$"))
        shape: (3, 1)
        ┌──────┐
        │ cc   │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ null │
        │ 2.5  │
        │ 1.5  │
        └──────┘

        Exclude by dtype(s), e.g. removing all columns of type Int64 or Float64:

        >>> df.select(pl.all().exclude([pl.Int64, pl.Float64]))
        shape: (3, 1)
        ┌──────┐
        │ ba   │
        │ ---  │
        │ str  │
        ╞══════╡
        │ a    │
        │ b    │
        │ null │
        └──────┘
        """
        return self.meta.as_selector().exclude(columns, *more_columns).as_expr()

    def pipe(
        self,
        function: Callable[Concatenate[Expr, P], T],
        *args: P.args,
        **kwargs: P.kwargs,
    ) -> T:
        r'''
        Offers a structured way to apply a sequence of user-defined functions (UDFs).

        Parameters
        ----------
        function
            Callable; will receive the expression as the first parameter,
            followed by any given args/kwargs.
        *args
            Arguments to pass to the UDF.
        **kwargs
            Keyword arguments to pass to the UDF.

        Examples
        --------
        >>> def extract_number(expr: pl.Expr) -> pl.Expr:
        ...     """Extract the digits from a string."""
        ...     return expr.str.extract(r"\d+", 0).cast(pl.Int64)
        >>>
        >>> def scale_negative_even(expr: pl.Expr, *, n: int = 1) -> pl.Expr:
        ...     """Set even numbers negative, and scale by a user-supplied value."""
        ...     expr = pl.when(expr % 2 == 0).then(-expr).otherwise(expr)
        ...     return expr * n
        >>>
        >>> df = pl.DataFrame({"val": ["a: 1", "b: 2", "c: 3", "d: 4"]})
        >>> df.with_columns(
        ...     udfs=(
        ...         pl.col("val").pipe(extract_number).pipe(scale_negative_even, n=5)
        ...     ),
        ... )
        shape: (4, 2)
        ┌──────┬──────┐
        │ val  ┆ udfs │
        │ ---  ┆ ---  │
        │ str  ┆ i64  │
        ╞══════╪══════╡
        │ a: 1 ┆ 5    │
        │ b: 2 ┆ -10  │
        │ c: 3 ┆ 15   │
        │ d: 4 ┆ -20  │
        └──────┴──────┘

        '''
        return function(self, *args, **kwargs)

    def not_(self) -> Expr:
        """
        Negate a boolean expression.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [True, False, False],
        ...         "b": ["a", "b", None],
        ...     }
        ... )
        >>> df
        shape: (3, 2)
        ┌───────┬──────┐
        │ a     ┆ b    │
        │ ---   ┆ ---  │
        │ bool  ┆ str  │
        ╞═══════╪══════╡
        │ true  ┆ a    │
        │ false ┆ b    │
        │ false ┆ null │
        └───────┴──────┘
        >>> df.select(pl.col("a").not_())
        shape: (3, 1)
        ┌───────┐
        │ a     │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ false │
        │ true  │
        │ true  │
        └───────┘
        """
        return wrap_expr(self._pyexpr.not_())

    def is_null(self) -> Expr:
        """
        Returns a boolean Series indicating which values are null.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None, 1, 5],
        ...         "b": [1.0, 2.0, float("nan"), 1.0, 5.0],
        ...     }
        ... )
        >>> df.with_columns(pl.all().is_null().name.suffix("_isnull"))  # nan != null
        shape: (5, 4)
        ┌──────┬─────┬──────────┬──────────┐
        │ a    ┆ b   ┆ a_isnull ┆ b_isnull │
        │ ---  ┆ --- ┆ ---      ┆ ---      │
        │ i64  ┆ f64 ┆ bool     ┆ bool     │
        ╞══════╪═════╪══════════╪══════════╡
        │ 1    ┆ 1.0 ┆ false    ┆ false    │
        │ 2    ┆ 2.0 ┆ false    ┆ false    │
        │ null ┆ NaN ┆ true     ┆ false    │
        │ 1    ┆ 1.0 ┆ false    ┆ false    │
        │ 5    ┆ 5.0 ┆ false    ┆ false    │
        └──────┴─────┴──────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.is_null())

    def is_not_null(self) -> Expr:
        """
        Returns a boolean Series indicating which values are not null.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None, 1, 5],
        ...         "b": [1.0, 2.0, float("nan"), 1.0, 5.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.all().is_not_null().name.suffix("_not_null")  # nan != null
        ... )
        shape: (5, 4)
        ┌──────┬─────┬────────────┬────────────┐
        │ a    ┆ b   ┆ a_not_null ┆ b_not_null │
        │ ---  ┆ --- ┆ ---        ┆ ---        │
        │ i64  ┆ f64 ┆ bool       ┆ bool       │
        ╞══════╪═════╪════════════╪════════════╡
        │ 1    ┆ 1.0 ┆ true       ┆ true       │
        │ 2    ┆ 2.0 ┆ true       ┆ true       │
        │ null ┆ NaN ┆ false      ┆ true       │
        │ 1    ┆ 1.0 ┆ true       ┆ true       │
        │ 5    ┆ 5.0 ┆ true       ┆ true       │
        └──────┴─────┴────────────┴────────────┘
        """
        return wrap_expr(self._pyexpr.is_not_null())

    def is_finite(self) -> Expr:
        """
        Returns a boolean Series indicating which values are finite.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "A": [1.0, 2],
        ...         "B": [3.0, float("inf")],
        ...     }
        ... )
        >>> df.select(pl.all().is_finite())
        shape: (2, 2)
        ┌──────┬───────┐
        │ A    ┆ B     │
        │ ---  ┆ ---   │
        │ bool ┆ bool  │
        ╞══════╪═══════╡
        │ true ┆ true  │
        │ true ┆ false │
        └──────┴───────┘
        """
        return wrap_expr(self._pyexpr.is_finite())

    def is_infinite(self) -> Expr:
        """
        Returns a boolean Series indicating which values are infinite.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "A": [1.0, 2],
        ...         "B": [3.0, float("inf")],
        ...     }
        ... )
        >>> df.select(pl.all().is_infinite())
        shape: (2, 2)
        ┌───────┬───────┐
        │ A     ┆ B     │
        │ ---   ┆ ---   │
        │ bool  ┆ bool  │
        ╞═══════╪═══════╡
        │ false ┆ false │
        │ false ┆ true  │
        └───────┴───────┘
        """
        return wrap_expr(self._pyexpr.is_infinite())

    def is_nan(self) -> Expr:
        """
        Returns a boolean Series indicating which values are NaN.

        Notes
        -----
        Floating point `NaN` (Not A Number) should not be confused
        with missing data represented as `Null/None`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None, 1, 5],
        ...         "b": [1.0, 2.0, float("nan"), 1.0, 5.0],
        ...     }
        ... )
        >>> df.with_columns(pl.col(pl.Float64).is_nan().name.suffix("_isnan"))
        shape: (5, 3)
        ┌──────┬─────┬─────────┐
        │ a    ┆ b   ┆ b_isnan │
        │ ---  ┆ --- ┆ ---     │
        │ i64  ┆ f64 ┆ bool    │
        ╞══════╪═════╪═════════╡
        │ 1    ┆ 1.0 ┆ false   │
        │ 2    ┆ 2.0 ┆ false   │
        │ null ┆ NaN ┆ true    │
        │ 1    ┆ 1.0 ┆ false   │
        │ 5    ┆ 5.0 ┆ false   │
        └──────┴─────┴─────────┘
        """
        return wrap_expr(self._pyexpr.is_nan())

    def is_not_nan(self) -> Expr:
        """
        Returns a boolean Series indicating which values are not NaN.

        Notes
        -----
        Floating point `NaN` (Not A Number) should not be confused
        with missing data represented as `Null/None`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None, 1, 5],
        ...         "b": [1.0, 2.0, float("nan"), 1.0, 5.0],
        ...     }
        ... )
        >>> df.with_columns(pl.col(pl.Float64).is_not_nan().name.suffix("_is_not_nan"))
        shape: (5, 3)
        ┌──────┬─────┬──────────────┐
        │ a    ┆ b   ┆ b_is_not_nan │
        │ ---  ┆ --- ┆ ---          │
        │ i64  ┆ f64 ┆ bool         │
        ╞══════╪═════╪══════════════╡
        │ 1    ┆ 1.0 ┆ true         │
        │ 2    ┆ 2.0 ┆ true         │
        │ null ┆ NaN ┆ false        │
        │ 1    ┆ 1.0 ┆ true         │
        │ 5    ┆ 5.0 ┆ true         │
        └──────┴─────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.is_not_nan())

    def agg_groups(self) -> Expr:
        """
        Get the group indexes of the group by operation.

        .. deprecated:: 1.35
            use `df.with_row_index().group_by(...).agg(pl.col('index'))` instead.
            This method will be removed in Polars 2.0.

        Should be used in aggregation context only.

        Examples
        --------
        >>> import warnings
        >>> warnings.filterwarnings("ignore", category=DeprecationWarning)
        >>> df = pl.DataFrame(
        ...     {
        ...         "group": [
        ...             "one",
        ...             "one",
        ...             "one",
        ...             "two",
        ...             "two",
        ...             "two",
        ...         ],
        ...         "value": [94, 95, 96, 97, 97, 99],
        ...     }
        ... )
        >>> df.group_by("group", maintain_order=True).agg(pl.col("value").agg_groups())
        shape: (2, 2)
        ┌───────┬───────────┐
        │ group ┆ value     │
        │ ---   ┆ ---       │
        │ str   ┆ list[u32] │
        ╞═══════╪═══════════╡
        │ one   ┆ [0, 1, 2] │
        │ two   ┆ [3, 4, 5] │
        └───────┴───────────┘

        New recommended approach:
        >>> (
        ...     df.with_row_index()
        ...     .group_by("group", maintain_order=True)
        ...     .agg(pl.col("index"))
        ... )
        shape: (2, 2)
        ┌───────┬───────────┐
        │ group ┆ index     │
        │ ---   ┆ ---       │
        │ str   ┆ list[u32] │
        ╞═══════╪═══════════╡
        │ one   ┆ [0, 1, 2] │
        │ two   ┆ [3, 4, 5] │
        └───────┴───────────┘
        """
        warnings.warn(
            "agg_groups() is deprecated and will be removed in Polars 2.0. "
            "Use df.with_row_index().group_by(...).agg(pl.col('index')) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return wrap_expr(self._pyexpr.agg_groups())

    def count(self) -> Expr:
        """
        Return the number of non-null elements in the column.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32`.

        See Also
        --------
        len

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": [None, 4, 4]})
        >>> df.select(pl.all().count())
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ u32 ┆ u32 │
        ╞═════╪═════╡
        │ 3   ┆ 2   │
        └─────┴─────┘
        """
        return wrap_expr(self._pyexpr.count())

    def len(self) -> Expr:
        """
        Return the number of elements in the column.

        Null values count towards the total.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32`.

        See Also
        --------
        count

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": [None, 4, 4]})
        >>> df.select(pl.all().len())
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ u32 ┆ u32 │
        ╞═════╪═════╡
        │ 3   ┆ 3   │
        └─────┴─────┘
        """
        return wrap_expr(self._pyexpr.len())

    def slice(self, offset: int | Expr, length: int | Expr | None = None) -> Expr:
        """
        Get a slice of this expression.

        Parameters
        ----------
        offset
            Start index. Negative indexing is supported.
        length
            Length of the slice. If set to `None`, all rows starting at the offset
            will be selected.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [8, 9, 10, 11],
        ...         "b": [None, 4, 4, 4],
        ...     }
        ... )
        >>> df.select(pl.all().slice(1, 2))
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 9   ┆ 4   │
        │ 10  ┆ 4   │
        └─────┴─────┘
        """
        if not isinstance(offset, Expr):
            offset = F.lit(offset)
        if not isinstance(length, Expr):
            length = F.lit(length)
        return wrap_expr(self._pyexpr.slice(offset._pyexpr, length._pyexpr))

    def append(self, other: IntoExpr, *, upcast: bool = True) -> Expr:
        """
        Append expressions.

        This is done by adding the chunks of `other` to this `Series`.

        Parameters
        ----------
        other
            Expression to append.
        upcast
            Cast both `Series` to the same supertype.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [8, 9, 10],
        ...         "b": [None, 4, 4],
        ...     }
        ... )
        >>> df.select(pl.all().head(1).append(pl.all().tail(1)))
        shape: (2, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ i64 ┆ i64  │
        ╞═════╪══════╡
        │ 8   ┆ null │
        │ 10  ┆ 4    │
        └─────┴──────┘
        """
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr.append(other_pyexpr, upcast))

    def rechunk(self) -> Expr:
        """
        Create a single chunk of memory for this Series.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2]})

        Create a Series with 3 nulls, append column `a`, then rechunk.

        >>> df.select(pl.repeat(None, 3).append(pl.col("a")).rechunk())
        shape: (6, 1)
        ┌────────┐
        │ repeat │
        │ ---    │
        │ i64    │
        ╞════════╡
        │ null   │
        │ null   │
        │ null   │
        │ 1      │
        │ 1      │
        │ 2      │
        └────────┘
        """
        return wrap_expr(self._pyexpr.rechunk())

    def drop_nulls(self) -> Expr:
        """
        Drop all null values.

        The original order of the remaining elements is preserved.

        See Also
        --------
        drop_nans

        Notes
        -----
        A null value is not the same as a NaN value.
        To drop NaN values, use :func:`drop_nans`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0, None, 3.0, float("nan")]})
        >>> df.select(pl.col("a").drop_nulls())
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        │ 3.0 │
        │ NaN │
        └─────┘
        """
        return wrap_expr(self._pyexpr.drop_nulls())

    def drop_nans(self) -> Expr:
        """
        Drop all floating point NaN values.

        The original order of the remaining elements is preserved.

        See Also
        --------
        drop_nulls

        Notes
        -----
        A NaN value is not the same as a null value.
        To drop null values, use :func:`drop_nulls`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0, None, 3.0, float("nan")]})
        >>> df.select(pl.col("a").drop_nans())
        shape: (3, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ 1.0  │
        │ null │
        │ 3.0  │
        └──────┘
        """
        return wrap_expr(self._pyexpr.drop_nans())

    def cum_sum(self, *, reverse: bool = False) -> Expr:
        """
        Get an array with the cumulative sum computed at every element.

        Parameters
        ----------
        reverse
            Reverse the operation.

        Notes
        -----
        Dtypes in {Int8, UInt8, Int16, UInt16} are cast to
        Int64 before summing to prevent overflow issues.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 4]})
        >>> df.with_columns(
        ...     pl.col("a").cum_sum().alias("cum_sum"),
        ...     pl.col("a").cum_sum(reverse=True).alias("cum_sum_reverse"),
        ... )
        shape: (4, 3)
        ┌─────┬─────────┬─────────────────┐
        │ a   ┆ cum_sum ┆ cum_sum_reverse │
        │ --- ┆ ---     ┆ ---             │
        │ i64 ┆ i64     ┆ i64             │
        ╞═════╪═════════╪═════════════════╡
        │ 1   ┆ 1       ┆ 10              │
        │ 2   ┆ 3       ┆ 9               │
        │ 3   ┆ 6       ┆ 7               │
        │ 4   ┆ 10      ┆ 4               │
        └─────┴─────────┴─────────────────┘

        Null values are excluded, but can also be filled by calling
        `fill_null(strategy="forward")`.

        >>> df = pl.DataFrame({"values": [None, 10, None, 8, 9, None, 16, None]})
        >>> df.with_columns(
        ...     pl.col("values").cum_sum().alias("value_cum_sum"),
        ...     pl.col("values")
        ...     .cum_sum()
        ...     .fill_null(strategy="forward")
        ...     .alias("value_cum_sum_all_filled"),
        ... )
        shape: (8, 3)
        ┌────────┬───────────────┬──────────────────────────┐
        │ values ┆ value_cum_sum ┆ value_cum_sum_all_filled │
        │ ---    ┆ ---           ┆ ---                      │
        │ i64    ┆ i64           ┆ i64                      │
        ╞════════╪═══════════════╪══════════════════════════╡
        │ null   ┆ null          ┆ null                     │
        │ 10     ┆ 10            ┆ 10                       │
        │ null   ┆ null          ┆ 10                       │
        │ 8      ┆ 18            ┆ 18                       │
        │ 9      ┆ 27            ┆ 27                       │
        │ null   ┆ null          ┆ 27                       │
        │ 16     ┆ 43            ┆ 43                       │
        │ null   ┆ null          ┆ 43                       │
        └────────┴───────────────┴──────────────────────────┘
        """
        return wrap_expr(self._pyexpr.cum_sum(reverse))

    def cum_prod(self, *, reverse: bool = False) -> Expr:
        """
        Get an array with the cumulative product computed at every element.

        Parameters
        ----------
        reverse
            Reverse the operation.

        Notes
        -----
        Dtypes in {Int8, UInt8, Int16, UInt16} are cast to
        Int64 before summing to prevent overflow issues.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 4]})
        >>> df.with_columns(
        ...     pl.col("a").cum_prod().alias("cum_prod"),
        ...     pl.col("a").cum_prod(reverse=True).alias("cum_prod_reverse"),
        ... )
        shape: (4, 3)
        ┌─────┬──────────┬──────────────────┐
        │ a   ┆ cum_prod ┆ cum_prod_reverse │
        │ --- ┆ ---      ┆ ---              │
        │ i64 ┆ i64      ┆ i64              │
        ╞═════╪══════════╪══════════════════╡
        │ 1   ┆ 1        ┆ 24               │
        │ 2   ┆ 2        ┆ 24               │
        │ 3   ┆ 6        ┆ 12               │
        │ 4   ┆ 24       ┆ 4                │
        └─────┴──────────┴──────────────────┘
        """
        return wrap_expr(self._pyexpr.cum_prod(reverse))

    def cum_min(self, *, reverse: bool = False) -> Expr:
        """
        Get an array with the cumulative min computed at every element.

        Parameters
        ----------
        reverse
            Reverse the operation.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [3, 1, 2]})
        >>> df.with_columns(
        ...     pl.col("a").cum_min().alias("cum_min"),
        ...     pl.col("a").cum_min(reverse=True).alias("cum_min_reverse"),
        ... )
        shape: (3, 3)
        ┌─────┬─────────┬─────────────────┐
        │ a   ┆ cum_min ┆ cum_min_reverse │
        │ --- ┆ ---     ┆ ---             │
        │ i64 ┆ i64     ┆ i64             │
        ╞═════╪═════════╪═════════════════╡
        │ 3   ┆ 3       ┆ 1               │
        │ 1   ┆ 1       ┆ 1               │
        │ 2   ┆ 1       ┆ 2               │
        └─────┴─────────┴─────────────────┘
        """
        return wrap_expr(self._pyexpr.cum_min(reverse))

    def cum_max(self, *, reverse: bool = False) -> Expr:
        """
        Get an array with the cumulative max computed at every element.

        Parameters
        ----------
        reverse
            Reverse the operation.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 3, 2]})
        >>> df.with_columns(
        ...     pl.col("a").cum_max().alias("cum_max"),
        ...     pl.col("a").cum_max(reverse=True).alias("cum_max_reverse"),
        ... )
        shape: (3, 3)
        ┌─────┬─────────┬─────────────────┐
        │ a   ┆ cum_max ┆ cum_max_reverse │
        │ --- ┆ ---     ┆ ---             │
        │ i64 ┆ i64     ┆ i64             │
        ╞═════╪═════════╪═════════════════╡
        │ 1   ┆ 1       ┆ 3               │
        │ 3   ┆ 3       ┆ 3               │
        │ 2   ┆ 3       ┆ 2               │
        └─────┴─────────┴─────────────────┘


        Null values are excluded, but can also be filled by calling
        `fill_null(strategy="forward")`.

        >>> df = pl.DataFrame({"values": [None, 10, None, 8, 9, None, 16, None]})
        >>> df.with_columns(
        ...     pl.col("values").cum_max().alias("cum_max"),
        ...     pl.col("values")
        ...     .cum_max()
        ...     .fill_null(strategy="forward")
        ...     .alias("cum_max_all_filled"),
        ... )
        shape: (8, 3)
        ┌────────┬─────────┬────────────────────┐
        │ values ┆ cum_max ┆ cum_max_all_filled │
        │ ---    ┆ ---     ┆ ---                │
        │ i64    ┆ i64     ┆ i64                │
        ╞════════╪═════════╪════════════════════╡
        │ null   ┆ null    ┆ null               │
        │ 10     ┆ 10      ┆ 10                 │
        │ null   ┆ null    ┆ 10                 │
        │ 8      ┆ 10      ┆ 10                 │
        │ 9      ┆ 10      ┆ 10                 │
        │ null   ┆ null    ┆ 10                 │
        │ 16     ┆ 16      ┆ 16                 │
        │ null   ┆ null    ┆ 16                 │
        └────────┴─────────┴────────────────────┘
        """
        return wrap_expr(self._pyexpr.cum_max(reverse))

    def cum_count(self, *, reverse: bool = False) -> Expr:
        """
        Return the cumulative count of the non-null values in the column.

        Parameters
        ----------
        reverse
            Reverse the operation.

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["x", "k", None, "d"]})
        >>> df.with_columns(
        ...     pl.col("a").cum_count().alias("cum_count"),
        ...     pl.col("a").cum_count(reverse=True).alias("cum_count_reverse"),
        ... )
        shape: (4, 3)
        ┌──────┬───────────┬───────────────────┐
        │ a    ┆ cum_count ┆ cum_count_reverse │
        │ ---  ┆ ---       ┆ ---               │
        │ str  ┆ u32       ┆ u32               │
        ╞══════╪═══════════╪═══════════════════╡
        │ x    ┆ 1         ┆ 3                 │
        │ k    ┆ 2         ┆ 2                 │
        │ null ┆ 2         ┆ 1                 │
        │ d    ┆ 3         ┆ 1                 │
        └──────┴───────────┴───────────────────┘
        """
        return wrap_expr(self._pyexpr.cum_count(reverse))

    def floor(self) -> Expr:
        """
        Rounds down to the nearest integer value.

        Only works on floating point Series.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.3, 0.5, 1.0, 1.1]})
        >>> df.select(pl.col("a").floor())
        shape: (4, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.0 │
        │ 0.0 │
        │ 1.0 │
        │ 1.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.floor())

    def ceil(self) -> Expr:
        """
        Rounds up to the nearest integer value.

        Only works on floating point Series.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.3, 0.5, 1.0, 1.1]})
        >>> df.select(pl.col("a").ceil())
        shape: (4, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        │ 1.0 │
        │ 1.0 │
        │ 2.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.ceil())

    def round(self, decimals: int = 0, mode: RoundMode = "half_to_even") -> Expr:
        """
        Round underlying floating point data by `decimals` digits.

        The default rounding mode is "half to even" (also known as "bankers' rounding").

        Parameters
        ----------
        decimals
            Number of decimals to round by.
        mode : {'half_to_even', 'half_away_from_zero'}
            RoundMode.

            * *half_to_even*
                round to the nearest even number
            * *half_away_from_zero*
                round to the nearest number away from zero

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.33, 0.52, 1.02, 1.17]})
        >>> df.select(pl.col("a").round(1))
        shape: (4, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.3 │
        │ 0.5 │
        │ 1.0 │
        │ 1.2 │
        └─────┘

        >>> df = pl.DataFrame(
        ...     {
        ...         "f64": [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5],
        ...         "d": ["-3.5", "-2.5", "-1.5", "-0.5", "0.5", "1.5", "2.5", "3.5"],
        ...     },
        ...     schema_overrides={"d": pl.Decimal(scale=1)},
        ... )
        >>> df.with_columns(
        ...     pl.all().round(mode="half_away_from_zero").name.suffix("_away"),
        ...     pl.all().round(mode="half_to_even").name.suffix("_to_even"),
        ... )
        shape: (8, 6)
        ┌──────┬───────────────┬──────────┬───────────────┬─────────────┬───────────────┐
        │ f64  ┆ d             ┆ f64_away ┆ d_away        ┆ f64_to_even ┆ d_to_even     │
        │ ---  ┆ ---           ┆ ---      ┆ ---           ┆ ---         ┆ ---           │
        │ f64  ┆ decimal[38,1] ┆ f64      ┆ decimal[38,1] ┆ f64         ┆ decimal[38,1] │
        ╞══════╪═══════════════╪══════════╪═══════════════╪═════════════╪═══════════════╡
        │ -3.5 ┆ -3.5          ┆ -4.0     ┆ -4.0          ┆ -4.0        ┆ -4.0          │
        │ -2.5 ┆ -2.5          ┆ -3.0     ┆ -3.0          ┆ -2.0        ┆ -2.0          │
        │ -1.5 ┆ -1.5          ┆ -2.0     ┆ -2.0          ┆ -2.0        ┆ -2.0          │
        │ -0.5 ┆ -0.5          ┆ -1.0     ┆ -1.0          ┆ -0.0        ┆ 0.0           │
        │ 0.5  ┆ 0.5           ┆ 1.0      ┆ 1.0           ┆ 0.0         ┆ 0.0           │
        │ 1.5  ┆ 1.5           ┆ 2.0      ┆ 2.0           ┆ 2.0         ┆ 2.0           │
        │ 2.5  ┆ 2.5           ┆ 3.0      ┆ 3.0           ┆ 2.0         ┆ 2.0           │
        │ 3.5  ┆ 3.5           ┆ 4.0      ┆ 4.0           ┆ 4.0         ┆ 4.0           │
        └──────┴───────────────┴──────────┴───────────────┴─────────────┴───────────────┘
        """  # noqa: W505
        return wrap_expr(self._pyexpr.round(decimals, mode))

    def round_sig_figs(self, digits: int) -> Expr:
        """
        Round to a number of significant figures.

        Parameters
        ----------
        digits
            Number of significant figures to round to.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.01234, 3.333, 1234.0]})
        >>> df.with_columns(pl.col("a").round_sig_figs(2).alias("round_sig_figs"))
        shape: (3, 2)
        ┌─────────┬────────────────┐
        │ a       ┆ round_sig_figs │
        │ ---     ┆ ---            │
        │ f64     ┆ f64            │
        ╞═════════╪════════════════╡
        │ 0.01234 ┆ 0.012          │
        │ 3.333   ┆ 3.3            │
        │ 1234.0  ┆ 1200.0         │
        └─────────┴────────────────┘
        """
        return wrap_expr(self._pyexpr.round_sig_figs(digits))

    def dot(self, other: Expr | str) -> Expr:
        """
        Compute the dot/inner product between two Expressions.

        Parameters
        ----------
        other
            Expression to compute dot product with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> df.select(pl.col("a").dot(pl.col("b")))
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 44  │
        └─────┘
        """
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr.dot(other_pyexpr))

    def mode(self) -> Expr:
        """
        Compute the most occurring value(s).

        Can return multiple Values.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 1, 2, 3],
        ...         "b": [1, 1, 2, 2],
        ...     }
        ... )
        >>> df.select(pl.all().mode().first())  # doctest: +IGNORE_RESULT
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 1   │
        └─────┴─────┘
        """
        return wrap_expr(self._pyexpr.mode())

    def cast(
        self,
        dtype: PolarsDataType | pl.DataTypeExpr | type[Any],
        *,
        strict: bool = True,
        wrap_numerical: bool = False,
    ) -> Expr:
        r"""
        Cast between data types.

        Parameters
        ----------
        dtype
            DataType to cast to.
        strict
            Raise if cast is invalid on rows after predicates are pushed down.
            If `False`, invalid casts will produce null values.
        wrap_numerical
            If True numeric casts wrap overflowing values instead of
            marking the cast as invalid.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": ["4", "5", "6"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("a").cast(pl.Float64),
        ...     pl.col("b").cast(pl.Int32),
        ... )
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ f64 ┆ i32 │
        ╞═════╪═════╡
        │ 1.0 ┆ 4   │
        │ 2.0 ┆ 5   │
        │ 3.0 ┆ 6   │
        └─────┴─────┘
        """
        dtype = parse_into_datatype_expr(dtype)
        return wrap_expr(
            self._pyexpr.cast(dtype._pydatatype_expr, strict, wrap_numerical)
        )

    def sort(self, *, descending: bool = False, nulls_last: bool = False) -> Expr:
        """
        Sort this column.

        When used in a projection/selection context, the whole column is sorted.
        When used in a group by context, the groups are sorted.

        Parameters
        ----------
        descending
            Sort in descending order.
        nulls_last
            Place null values last.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, None, 3, 2],
        ...     }
        ... )
        >>> df.select(pl.col("a").sort())
        shape: (4, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ i64  │
        ╞══════╡
        │ null │
        │ 1    │
        │ 2    │
        │ 3    │
        └──────┘
        >>> df.select(pl.col("a").sort(descending=True))
        shape: (4, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ i64  │
        ╞══════╡
        │ null │
        │ 3    │
        │ 2    │
        │ 1    │
        └──────┘
        >>> df.select(pl.col("a").sort(nulls_last=True))
        shape: (4, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ i64  │
        ╞══════╡
        │ 1    │
        │ 2    │
        │ 3    │
        │ null │
        └──────┘

        When sorting in a group by context, the groups are sorted.

        >>> df = pl.DataFrame(
        ...     {
        ...         "group": ["one", "one", "one", "two", "two", "two"],
        ...         "value": [1, 98, 2, 3, 99, 4],
        ...     }
        ... )
        >>> df.group_by("group").agg(pl.col("value").sort())  # doctest: +IGNORE_RESULT
        shape: (2, 2)
        ┌───────┬────────────┐
        │ group ┆ value      │
        │ ---   ┆ ---        │
        │ str   ┆ list[i64]  │
        ╞═══════╪════════════╡
        │ two   ┆ [3, 4, 99] │
        │ one   ┆ [1, 2, 98] │
        └───────┴────────────┘
        """
        return wrap_expr(self._pyexpr.sort_with(descending, nulls_last))

    def top_k(self, k: int | IntoExprColumn = 5) -> Expr:
        r"""
        Return the `k` largest elements.

        Non-null elements are always preferred over null elements. The output
        is not guaranteed to be in any particular order, call :func:`sort`
        after this function if you wish the output to be sorted.

        This has time complexity:

        .. math:: O(n)

        Parameters
        ----------
        k
            Number of elements to return.

        See Also
        --------
        top_k_by
        bottom_k
        bottom_k_by

        Examples
        --------
        Get the 5 largest values in series.

        >>> df = pl.DataFrame({"value": [1, 98, 2, 3, 99, 4]})
        >>> df.select(
        ...     pl.col("value").top_k().alias("top_k"),
        ...     pl.col("value").bottom_k().alias("bottom_k"),
        ... )
        shape: (5, 2)
        ┌───────┬──────────┐
        │ top_k ┆ bottom_k │
        │ ---   ┆ ---      │
        │ i64   ┆ i64      │
        ╞═══════╪══════════╡
        │ 4     ┆ 1        │
        │ 98    ┆ 98       │
        │ 2     ┆ 2        │
        │ 3     ┆ 3        │
        │ 99    ┆ 4        │
        └───────┴──────────┘
        """
        k_pyexpr = parse_into_expression(k)
        return wrap_expr(self._pyexpr.top_k(k_pyexpr))

    @deprecate_renamed_parameter("descending", "reverse", version="1.0.0")
    def top_k_by(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        k: int | IntoExprColumn = 5,
        *,
        reverse: bool | Sequence[bool] = False,
    ) -> Expr:
        r"""
        Return the elements corresponding to the `k` largest elements of the `by` column(s).

        Non-null elements are always preferred over null elements, regardless of
        the value of `reverse`. The output is not guaranteed to be in any
        particular order, call :func:`sort` after this function if you wish the
        output to be sorted.

        This has time complexity:

        .. math:: O(n \log{n})

        .. versionchanged:: 1.0.0
            The `descending` parameter was renamed to `reverse`.

        Parameters
        ----------
        by
            Column(s) used to determine the largest elements.
            Accepts expression input. Strings are parsed as column names.
        k
            Number of elements to return.
        reverse
            Consider the `k` smallest elements of the `by` column(s) (instead of the `k`
            largest). This can be specified per column by passing a sequence of
            booleans.

        See Also
        --------
        top_k
        bottom_k
        bottom_k_by

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4, 5, 6],
        ...         "b": [6, 5, 4, 3, 2, 1],
        ...         "c": ["Apple", "Orange", "Apple", "Apple", "Banana", "Banana"],
        ...     }
        ... )
        >>> df
        shape: (6, 3)
        ┌─────┬─────┬────────┐
        │ a   ┆ b   ┆ c      │
        │ --- ┆ --- ┆ ---    │
        │ i64 ┆ i64 ┆ str    │
        ╞═════╪═════╪════════╡
        │ 1   ┆ 6   ┆ Apple  │
        │ 2   ┆ 5   ┆ Orange │
        │ 3   ┆ 4   ┆ Apple  │
        │ 4   ┆ 3   ┆ Apple  │
        │ 5   ┆ 2   ┆ Banana │
        │ 6   ┆ 1   ┆ Banana │
        └─────┴─────┴────────┘

        Get the top 2 rows by column `a` or `b`.

        >>> df.select(
        ...     pl.all().top_k_by("a", 2).name.suffix("_top_by_a"),
        ...     pl.all().top_k_by("b", 2).name.suffix("_top_by_b"),
        ... )
        shape: (2, 6)
        ┌────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐
        │ a_top_by_a ┆ b_top_by_a ┆ c_top_by_a ┆ a_top_by_b ┆ b_top_by_b ┆ c_top_by_b │
        │ ---        ┆ ---        ┆ ---        ┆ ---        ┆ ---        ┆ ---        │
        │ i64        ┆ i64        ┆ str        ┆ i64        ┆ i64        ┆ str        │
        ╞════════════╪════════════╪════════════╪════════════╪════════════╪════════════╡
        │ 6          ┆ 1          ┆ Banana     ┆ 1          ┆ 6          ┆ Apple      │
        │ 5          ┆ 2          ┆ Banana     ┆ 2          ┆ 5          ┆ Orange     │
        └────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘

        Get the top 2 rows by multiple columns with given order.

        >>> df.select(
        ...     pl.all()
        ...     .top_k_by(["c", "a"], 2, reverse=[False, True])
        ...     .name.suffix("_by_ca"),
        ...     pl.all()
        ...     .top_k_by(["c", "b"], 2, reverse=[False, True])
        ...     .name.suffix("_by_cb"),
        ... )
        shape: (2, 6)
        ┌─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
        │ a_by_ca ┆ b_by_ca ┆ c_by_ca ┆ a_by_cb ┆ b_by_cb ┆ c_by_cb │
        │ ---     ┆ ---     ┆ ---     ┆ ---     ┆ ---     ┆ ---     │
        │ i64     ┆ i64     ┆ str     ┆ i64     ┆ i64     ┆ str     │
        ╞═════════╪═════════╪═════════╪═════════╪═════════╪═════════╡
        │ 2       ┆ 5       ┆ Orange  ┆ 2       ┆ 5       ┆ Orange  │
        │ 5       ┆ 2       ┆ Banana  ┆ 6       ┆ 1       ┆ Banana  │
        └─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘

        Get the top 2 rows by column `a` in each group.

        >>> (
        ...     df.group_by("c", maintain_order=True)
        ...     .agg(pl.all().top_k_by("a", 2))
        ...     .explode(pl.all().exclude("c"))
        ... )
        shape: (5, 3)
        ┌────────┬─────┬─────┐
        │ c      ┆ a   ┆ b   │
        │ ---    ┆ --- ┆ --- │
        │ str    ┆ i64 ┆ i64 │
        ╞════════╪═════╪═════╡
        │ Apple  ┆ 4   ┆ 3   │
        │ Apple  ┆ 3   ┆ 4   │
        │ Orange ┆ 2   ┆ 5   │
        │ Banana ┆ 6   ┆ 1   │
        │ Banana ┆ 5   ┆ 2   │
        └────────┴─────┴─────┘
        """  # noqa: W505
        k_pyexpr = parse_into_expression(k)
        by_pyexprs = parse_into_list_of_expressions(by)

        reverse = extend_bool(reverse, len(by_pyexprs), "reverse", "by")

        return wrap_expr(self._pyexpr.top_k_by(by_pyexprs, k=k_pyexpr, reverse=reverse))

    def bottom_k(self, k: int | IntoExprColumn = 5) -> Expr:
        r"""
        Return the `k` smallest elements.

        Non-null elements are always preferred over null elements. The output is
        not guaranteed to be in any particular order, call :func:`sort` after
        this function if you wish the output to be sorted.

        This has time complexity:

        .. math:: O(n)

        Parameters
        ----------
        k
            Number of elements to return.

        See Also
        --------
        top_k
        top_k_by
        bottom_k_by

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "value": [1, 98, 2, 3, 99, 4],
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("value").top_k().alias("top_k"),
        ...     pl.col("value").bottom_k().alias("bottom_k"),
        ... )
        shape: (5, 2)
        ┌───────┬──────────┐
        │ top_k ┆ bottom_k │
        │ ---   ┆ ---      │
        │ i64   ┆ i64      │
        ╞═══════╪══════════╡
        │ 4     ┆ 1        │
        │ 98    ┆ 98       │
        │ 2     ┆ 2        │
        │ 3     ┆ 3        │
        │ 99    ┆ 4        │
        └───────┴──────────┘
        """
        k_pyexpr = parse_into_expression(k)
        return wrap_expr(self._pyexpr.bottom_k(k_pyexpr))

    @deprecate_renamed_parameter("descending", "reverse", version="1.0.0")
    def bottom_k_by(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        k: int | IntoExprColumn = 5,
        *,
        reverse: bool | Sequence[bool] = False,
    ) -> Expr:
        r"""
        Return the elements corresponding to the `k` smallest elements of the `by` column(s).

        Non-null elements are always preferred over null elements, regardless of
        the value of `reverse`. The output is not guaranteed to be in any
        particular order, call :func:`sort` after this function if you wish the
        output to be sorted.

        This has time complexity:

        .. math:: O(n \log{n})

        .. versionchanged:: 1.0.0
            The `descending` parameter was renamed `reverse`.

        Parameters
        ----------
        by
            Column(s) used to determine the smallest elements.
            Accepts expression input. Strings are parsed as column names.
        k
            Number of elements to return.
        reverse
            Consider the `k` largest elements of the `by` column(s) (instead of the `k`
            smallest). This can be specified per column by passing a sequence of
            booleans.

        See Also
        --------
        top_k
        top_k_by
        bottom_k

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4, 5, 6],
        ...         "b": [6, 5, 4, 3, 2, 1],
        ...         "c": ["Apple", "Orange", "Apple", "Apple", "Banana", "Banana"],
        ...     }
        ... )
        >>> df
        shape: (6, 3)
        ┌─────┬─────┬────────┐
        │ a   ┆ b   ┆ c      │
        │ --- ┆ --- ┆ ---    │
        │ i64 ┆ i64 ┆ str    │
        ╞═════╪═════╪════════╡
        │ 1   ┆ 6   ┆ Apple  │
        │ 2   ┆ 5   ┆ Orange │
        │ 3   ┆ 4   ┆ Apple  │
        │ 4   ┆ 3   ┆ Apple  │
        │ 5   ┆ 2   ┆ Banana │
        │ 6   ┆ 1   ┆ Banana │
        └─────┴─────┴────────┘

        Get the bottom 2 rows by column `a` or `b`.

        >>> df.select(
        ...     pl.all().bottom_k_by("a", 2).name.suffix("_btm_by_a"),
        ...     pl.all().bottom_k_by("b", 2).name.suffix("_btm_by_b"),
        ... )
        shape: (2, 6)
        ┌────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐
        │ a_btm_by_a ┆ b_btm_by_a ┆ c_btm_by_a ┆ a_btm_by_b ┆ b_btm_by_b ┆ c_btm_by_b │
        │ ---        ┆ ---        ┆ ---        ┆ ---        ┆ ---        ┆ ---        │
        │ i64        ┆ i64        ┆ str        ┆ i64        ┆ i64        ┆ str        │
        ╞════════════╪════════════╪════════════╪════════════╪════════════╪════════════╡
        │ 1          ┆ 6          ┆ Apple      ┆ 6          ┆ 1          ┆ Banana     │
        │ 2          ┆ 5          ┆ Orange     ┆ 5          ┆ 2          ┆ Banana     │
        └────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘

        Get the bottom 2 rows by multiple columns with given order.

        >>> df.select(
        ...     pl.all()
        ...     .bottom_k_by(["c", "a"], 2, reverse=[False, True])
        ...     .name.suffix("_by_ca"),
        ...     pl.all()
        ...     .bottom_k_by(["c", "b"], 2, reverse=[False, True])
        ...     .name.suffix("_by_cb"),
        ... )
        shape: (2, 6)
        ┌─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
        │ a_by_ca ┆ b_by_ca ┆ c_by_ca ┆ a_by_cb ┆ b_by_cb ┆ c_by_cb │
        │ ---     ┆ ---     ┆ ---     ┆ ---     ┆ ---     ┆ ---     │
        │ i64     ┆ i64     ┆ str     ┆ i64     ┆ i64     ┆ str     │
        ╞═════════╪═════════╪═════════╪═════════╪═════════╪═════════╡
        │ 4       ┆ 3       ┆ Apple   ┆ 1       ┆ 6       ┆ Apple   │
        │ 3       ┆ 4       ┆ Apple   ┆ 3       ┆ 4       ┆ Apple   │
        └─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘

        Get the bottom 2 rows by column `a` in each group.

        >>> (
        ...     df.group_by("c", maintain_order=True)
        ...     .agg(pl.all().bottom_k_by("a", 2))
        ...     .explode(pl.all().exclude("c"))
        ... )
        shape: (5, 3)
        ┌────────┬─────┬─────┐
        │ c      ┆ a   ┆ b   │
        │ ---    ┆ --- ┆ --- │
        │ str    ┆ i64 ┆ i64 │
        ╞════════╪═════╪═════╡
        │ Apple  ┆ 1   ┆ 6   │
        │ Apple  ┆ 3   ┆ 4   │
        │ Orange ┆ 2   ┆ 5   │
        │ Banana ┆ 5   ┆ 2   │
        │ Banana ┆ 6   ┆ 1   │
        └────────┴─────┴─────┘
        """  # noqa: W505
        k_pyexpr = parse_into_expression(k)
        by_pyexpr = parse_into_list_of_expressions(by)
        reverse = extend_bool(reverse, len(by_pyexpr), "reverse", "by")
        return wrap_expr(
            self._pyexpr.bottom_k_by(by_pyexpr, k=k_pyexpr, reverse=reverse)
        )

    def arg_sort(self, *, descending: bool = False, nulls_last: bool = False) -> Expr:
        """
        Get the index values that would sort this column.

        Parameters
        ----------
        descending
            Sort in descending (descending) order.
        nulls_last
            Place null values last instead of first.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32`.

        See Also
        --------
        Expr.gather: Take values by index.
        Expr.rank : Get the rank of each row.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [20, 10, 30],
        ...         "b": [1, 2, 3],
        ...     }
        ... )
        >>> df.select(pl.col("a").arg_sort())
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 1   │
        │ 0   │
        │ 2   │
        └─────┘

        Use gather to apply the arg sort to other columns.

        >>> df.select(pl.col("b").gather(pl.col("a").arg_sort()))
        shape: (3, 1)
        ┌─────┐
        │ b   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        │ 1   │
        │ 3   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arg_sort(descending, nulls_last))

    def arg_max(self) -> Expr:
        """
        Get the index of the maximal value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [20, 10, 30],
        ...     }
        ... )
        >>> df.select(pl.col("a").arg_max())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 2   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arg_max())

    def arg_min(self) -> Expr:
        """
        Get the index of the minimal value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [20, 10, 30],
        ...     }
        ... )
        >>> df.select(pl.col("a").arg_min())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 1   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arg_min())

    def index_of(self, element: IntoExpr) -> Expr:
        """
        Get the index of the first occurrence of a value, or ``None`` if it's not found.

        Parameters
        ----------
        element
            Value to find.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, None, 17]})
        >>> df.select(
        ...     [
        ...         pl.col("a").index_of(17).alias("seventeen"),
        ...         pl.col("a").index_of(None).alias("null"),
        ...         pl.col("a").index_of(55).alias("fiftyfive"),
        ...     ]
        ... )
        shape: (1, 3)
        ┌───────────┬──────┬───────────┐
        │ seventeen ┆ null ┆ fiftyfive │
        │ ---       ┆ ---  ┆ ---       │
        │ u32       ┆ u32  ┆ u32       │
        ╞═══════════╪══════╪═══════════╡
        │ 2         ┆ 1    ┆ null      │
        └───────────┴──────┴───────────┘
        """
        element_pyexpr = parse_into_expression(element, str_as_lit=True)
        return wrap_expr(self._pyexpr.index_of(element_pyexpr))

    def search_sorted(
        self,
        element: IntoExpr | np.ndarray[Any, Any],
        side: SearchSortedSide = "any",
        *,
        descending: bool = False,
    ) -> Expr:
        """
        Find indices where elements should be inserted to maintain order.

        .. math:: a[i-1] < v <= a[i]

        Parameters
        ----------
        element
            Expression or scalar value.
        side : {'any', 'left', 'right'}
            If 'any', the index of the first suitable location found is given.
            If 'left', the index of the leftmost suitable location found is given.
            If 'right', return the rightmost suitable location found is given.
        descending
            Boolean indicating whether the values are descending or not (they
            are required to be sorted either way).

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "values": [1, 2, 3, 5],
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("values").search_sorted(0).alias("zero"),
        ...         pl.col("values").search_sorted(3).alias("three"),
        ...         pl.col("values").search_sorted(6).alias("six"),
        ...     ]
        ... )
        shape: (1, 3)
        ┌──────┬───────┬─────┐
        │ zero ┆ three ┆ six │
        │ ---  ┆ ---   ┆ --- │
        │ u32  ┆ u32   ┆ u32 │
        ╞══════╪═══════╪═════╡
        │ 0    ┆ 2     ┆ 4   │
        └──────┴───────┴─────┘
        """
        element_pyexpr = parse_into_expression(
            element, str_as_lit=True, list_as_series=True
        )
        return wrap_expr(self._pyexpr.search_sorted(element_pyexpr, side, descending))

    def sort_by(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        *more_by: IntoExpr,
        descending: bool | Sequence[bool] = False,
        nulls_last: bool | Sequence[bool] = False,
        multithreaded: bool = True,
        maintain_order: bool = False,
    ) -> Expr:
        """
        Sort this column by the ordering of other columns.

        When used in a projection/selection context, the whole column is sorted.
        When used in a group by context, the groups are sorted.

        Parameters
        ----------
        by
            Column(s) to sort by. Accepts expression input. Strings are parsed as column
            names.
        *more_by
            Additional columns to sort by, specified as positional arguments.
        descending
            Sort in descending order. When sorting by multiple columns, can be specified
            per column by passing a sequence of booleans.
        nulls_last
            Place null values last; can specify a single boolean applying to all columns
            or a sequence of booleans for per-column control.
        multithreaded
            Sort using multiple threads.
        maintain_order
            Whether the order should be maintained if elements are equal.

        Examples
        --------
        Pass a single column name to sort by that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "group": ["a", "a", "b", "b"],
        ...         "value1": [1, 3, 4, 2],
        ...         "value2": [8, 7, 6, 5],
        ...     }
        ... )
        >>> df.select(pl.col("group").sort_by("value1"))
        shape: (4, 1)
        ┌───────┐
        │ group │
        │ ---   │
        │ str   │
        ╞═══════╡
        │ a     │
        │ b     │
        │ a     │
        │ b     │
        └───────┘

        Sorting by expressions is also supported.

        >>> df.select(pl.col("group").sort_by(pl.col("value1") + pl.col("value2")))
        shape: (4, 1)
        ┌───────┐
        │ group │
        │ ---   │
        │ str   │
        ╞═══════╡
        │ b     │
        │ a     │
        │ a     │
        │ b     │
        └───────┘

        Sort by multiple columns by passing a list of columns.

        >>> df.select(pl.col("group").sort_by(["value1", "value2"], descending=True))
        shape: (4, 1)
        ┌───────┐
        │ group │
        │ ---   │
        │ str   │
        ╞═══════╡
        │ b     │
        │ a     │
        │ b     │
        │ a     │
        └───────┘

        Or use positional arguments to sort by multiple columns in the same way.

        >>> df.select(pl.col("group").sort_by("value1", "value2"))
        shape: (4, 1)
        ┌───────┐
        │ group │
        │ ---   │
        │ str   │
        ╞═══════╡
        │ a     │
        │ b     │
        │ a     │
        │ b     │
        └───────┘

        When sorting in a group by context, the groups are sorted.

        >>> df.group_by("group").agg(
        ...     pl.col("value1").sort_by("value2")
        ... )  # doctest: +IGNORE_RESULT
        shape: (2, 2)
        ┌───────┬───────────┐
        │ group ┆ value1    │
        │ ---   ┆ ---       │
        │ str   ┆ list[i64] │
        ╞═══════╪═══════════╡
        │ a     ┆ [3, 1]    │
        │ b     ┆ [2, 4]    │
        └───────┴───────────┘

        Take a single row from each group where a column attains its minimal value
        within that group.

        >>> df.group_by("group").agg(
        ...     pl.all().sort_by("value2").first()
        ... )  # doctest: +IGNORE_RESULT
        shape: (2, 3)
        ┌───────┬────────┬────────┐
        │ group ┆ value1 ┆ value2 |
        │ ---   ┆ ---    ┆ ---    │
        │ str   ┆ i64    ┆ i64    |
        ╞═══════╪════════╪════════╡
        │ a     ┆ 3      ┆ 7      |
        │ b     ┆ 2      ┆ 5      |
        └───────┴────────┴────────┘
        """
        by_pyexprs = parse_into_list_of_expressions(by, *more_by)
        descending = extend_bool(descending, len(by_pyexprs), "descending", "by")
        nulls_last = extend_bool(nulls_last, len(by_pyexprs), "nulls_last", "by")
        return wrap_expr(
            self._pyexpr.sort_by(
                by_pyexprs, descending, nulls_last, multithreaded, maintain_order
            )
        )

    def gather(
        self, indices: int | Sequence[int] | IntoExpr | Series | np.ndarray[Any, Any]
    ) -> Expr:
        """
        Take values by index.

        Parameters
        ----------
        indices
            An expression that leads to a UInt32 dtyped Series.

        Returns
        -------
        Expr
            Expression of the same data type.

        See Also
        --------
        Expr.get : Take a single value

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "group": [
        ...             "one",
        ...             "one",
        ...             "one",
        ...             "two",
        ...             "two",
        ...             "two",
        ...         ],
        ...         "value": [1, 98, 2, 3, 99, 4],
        ...     }
        ... )
        >>> df.group_by("group", maintain_order=True).agg(
        ...     pl.col("value").gather([2, 1])
        ... )
        shape: (2, 2)
        ┌───────┬───────────┐
        │ group ┆ value     │
        │ ---   ┆ ---       │
        │ str   ┆ list[i64] │
        ╞═══════╪═══════════╡
        │ one   ┆ [2, 98]   │
        │ two   ┆ [4, 99]   │
        └───────┴───────────┘
        """
        if (isinstance(indices, Sequence) and not isinstance(indices, str)) or (
            _check_for_numpy(indices) and isinstance(indices, np.ndarray)
        ):
            indices_lit_pyexpr = F.lit(pl.Series("", indices, dtype=Int64))._pyexpr
        else:
            indices_lit_pyexpr = parse_into_expression(indices)
        return wrap_expr(self._pyexpr.gather(indices_lit_pyexpr))

    def get(self, index: int | Expr) -> Expr:
        """
        Return a single value by index.

        Parameters
        ----------
        index
            An expression that leads to a UInt32 index.

        Returns
        -------
        Expr
            Expression of the same data type.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "group": [
        ...             "one",
        ...             "one",
        ...             "one",
        ...             "two",
        ...             "two",
        ...             "two",
        ...         ],
        ...         "value": [1, 98, 2, 3, 99, 4],
        ...     }
        ... )
        >>> df.group_by("group", maintain_order=True).agg(pl.col("value").get(1))
        shape: (2, 2)
        ┌───────┬───────┐
        │ group ┆ value │
        │ ---   ┆ ---   │
        │ str   ┆ i64   │
        ╞═══════╪═══════╡
        │ one   ┆ 98    │
        │ two   ┆ 99    │
        └───────┴───────┘
        """
        index_lit_pyexpr = parse_into_expression(index)
        return wrap_expr(self._pyexpr.get(index_lit_pyexpr))

    def shift(
        self, n: int | IntoExprColumn = 1, *, fill_value: IntoExpr | None = None
    ) -> Expr:
        """
        Shift values by the given number of indices.

        Parameters
        ----------
        n
            Number of indices to shift forward. If a negative value is passed, values
            are shifted in the opposite direction instead.
        fill_value
            Fill the resulting null values with this scalar value.

        Notes
        -----
        This method is similar to the `LAG` operation in SQL when the value for `n`
        is positive. With a negative value for `n`, it is similar to `LEAD`.

        See Also
        --------
        fill_null

        Examples
        --------
        By default, values are shifted forward by one index.

        >>> df = pl.DataFrame({"a": [1, 2, 3, 4]})
        >>> df.with_columns(shift=pl.col("a").shift())
        shape: (4, 2)
        ┌─────┬───────┐
        │ a   ┆ shift │
        │ --- ┆ ---   │
        │ i64 ┆ i64   │
        ╞═════╪═══════╡
        │ 1   ┆ null  │
        │ 2   ┆ 1     │
        │ 3   ┆ 2     │
        │ 4   ┆ 3     │
        └─────┴───────┘

        Pass a negative value to shift in the opposite direction instead.

        >>> df.with_columns(shift=pl.col("a").shift(-2))
        shape: (4, 2)
        ┌─────┬───────┐
        │ a   ┆ shift │
        │ --- ┆ ---   │
        │ i64 ┆ i64   │
        ╞═════╪═══════╡
        │ 1   ┆ 3     │
        │ 2   ┆ 4     │
        │ 3   ┆ null  │
        │ 4   ┆ null  │
        └─────┴───────┘

        Specify `fill_value` to fill the resulting null values.

        >>> df.with_columns(shift=pl.col("a").shift(-2, fill_value=100))
        shape: (4, 2)
        ┌─────┬───────┐
        │ a   ┆ shift │
        │ --- ┆ ---   │
        │ i64 ┆ i64   │
        ╞═════╪═══════╡
        │ 1   ┆ 3     │
        │ 2   ┆ 4     │
        │ 3   ┆ 100   │
        │ 4   ┆ 100   │
        └─────┴───────┘
        """
        if fill_value is not None:
            fill_value_pyexpr = parse_into_expression(fill_value, str_as_lit=True)
        else:
            fill_value_pyexpr = None
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.shift(n_pyexpr, fill_value_pyexpr))

    def fill_null(
        self,
        value: Any | Expr | None = None,
        strategy: FillNullStrategy | None = None,
        limit: int | None = None,
    ) -> Expr:
        """
        Fill null values using the specified value or strategy.

        To interpolate over null values see interpolate.
        See the examples below to fill nulls with an expression.

        Parameters
        ----------
        value
            Value used to fill null values.
        strategy : {None, 'forward', 'backward', 'min', 'max', 'mean', 'zero', 'one'}
            Strategy used to fill null values.
        limit
            Number of consecutive null values to fill when using the 'forward' or
            'backward' strategy.

        See Also
        --------
        backward_fill
        fill_nan
        forward_fill

        Notes
        -----
        A null value is not the same as a NaN value.
        To fill NaN values, use :func:`fill_nan`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None],
        ...         "b": [4, None, 6],
        ...     }
        ... )
        >>> df.with_columns(pl.col("b").fill_null(strategy="zero"))
        shape: (3, 2)
        ┌──────┬─────┐
        │ a    ┆ b   │
        │ ---  ┆ --- │
        │ i64  ┆ i64 │
        ╞══════╪═════╡
        │ 1    ┆ 4   │
        │ 2    ┆ 0   │
        │ null ┆ 6   │
        └──────┴─────┘
        >>> df.with_columns(pl.col("b").fill_null(99))
        shape: (3, 2)
        ┌──────┬─────┐
        │ a    ┆ b   │
        │ ---  ┆ --- │
        │ i64  ┆ i64 │
        ╞══════╪═════╡
        │ 1    ┆ 4   │
        │ 2    ┆ 99  │
        │ null ┆ 6   │
        └──────┴─────┘
        >>> df.with_columns(pl.col("b").fill_null(strategy="forward"))
        shape: (3, 2)
        ┌──────┬─────┐
        │ a    ┆ b   │
        │ ---  ┆ --- │
        │ i64  ┆ i64 │
        ╞══════╪═════╡
        │ 1    ┆ 4   │
        │ 2    ┆ 4   │
        │ null ┆ 6   │
        └──────┴─────┘
        >>> df.with_columns(pl.col("b").fill_null(pl.col("b").median()))
        shape: (3, 2)
        ┌──────┬─────┐
        │ a    ┆ b   │
        │ ---  ┆ --- │
        │ i64  ┆ f64 │
        ╞══════╪═════╡
        │ 1    ┆ 4.0 │
        │ 2    ┆ 5.0 │
        │ null ┆ 6.0 │
        └──────┴─────┘
        >>> df.with_columns(pl.all().fill_null(pl.all().median()))
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ f64 ┆ f64 │
        ╞═════╪═════╡
        │ 1.0 ┆ 4.0 │
        │ 2.0 ┆ 5.0 │
        │ 1.5 ┆ 6.0 │
        └─────┴─────┘
        """
        if value is not None and strategy is not None:
            msg = "cannot specify both `value` and `strategy`"
            raise ValueError(msg)
        elif value is None and strategy is None:
            msg = "must specify either a fill `value` or `strategy`"
            raise ValueError(msg)
        elif strategy not in ("forward", "backward") and limit is not None:
            msg = "can only specify `limit` when strategy is set to 'backward' or 'forward'"
            raise ValueError(msg)

        if value is not None:
            value_pyexpr = parse_into_expression(value, str_as_lit=True)
            return wrap_expr(self._pyexpr.fill_null(value_pyexpr))
        else:
            assert strategy is not None
            return wrap_expr(self._pyexpr.fill_null_with_strategy(strategy, limit))

    def fill_nan(self, value: int | float | Expr | None) -> Expr:
        """
        Fill floating point NaN value with a fill value.

        Parameters
        ----------
        value
            Value used to fill NaN values.

        See Also
        --------
        fill_null

        Notes
        -----
        A NaN value is not the same as a null value.
        To fill null values, use :func:`fill_null`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1.0, None, float("nan")],
        ...         "b": [4.0, float("nan"), 6],
        ...     }
        ... )
        >>> df.with_columns(pl.col("b").fill_nan(0))
        shape: (3, 2)
        ┌──────┬─────┐
        │ a    ┆ b   │
        │ ---  ┆ --- │
        │ f64  ┆ f64 │
        ╞══════╪═════╡
        │ 1.0  ┆ 4.0 │
        │ null ┆ 0.0 │
        │ NaN  ┆ 6.0 │
        └──────┴─────┘
        """
        fill_value_pyexpr = parse_into_expression(value, str_as_lit=True)
        return wrap_expr(self._pyexpr.fill_nan(fill_value_pyexpr))

    def forward_fill(self, limit: int | None = None) -> Expr:
        """
        Fill missing values with the last non-null value.

        This is an alias of `.fill_null(strategy="forward")`.

        Parameters
        ----------
        limit
            The number of consecutive null values to forward fill.

        See Also
        --------
        backward_fill
        fill_null
        shift
        """
        return self.fill_null(strategy="forward", limit=limit)

    def backward_fill(self, limit: int | None = None) -> Expr:
        """
        Fill missing values with the next non-null value.

        This is an alias of `.fill_null(strategy="backward")`.

        Parameters
        ----------
        limit
            The number of consecutive null values to backward fill.

        See Also
        --------
        fill_null
        forward_fill
        shift
        """
        return self.fill_null(strategy="backward", limit=limit)

    def reverse(self) -> Expr:
        """
        Reverse the selection.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "A": [1, 2, 3, 4, 5],
        ...         "fruits": ["banana", "banana", "apple", "apple", "banana"],
        ...         "B": [5, 4, 3, 2, 1],
        ...         "cars": ["beetle", "audi", "beetle", "beetle", "beetle"],
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.all(),
        ...         pl.all().reverse().name.suffix("_reverse"),
        ...     ]
        ... )
        shape: (5, 8)
        ┌─────┬────────┬─────┬────────┬───────────┬────────────────┬───────────┬──────────────┐
        │ A   ┆ fruits ┆ B   ┆ cars   ┆ A_reverse ┆ fruits_reverse ┆ B_reverse ┆ cars_reverse │
        │ --- ┆ ---    ┆ --- ┆ ---    ┆ ---       ┆ ---            ┆ ---       ┆ ---          │
        │ i64 ┆ str    ┆ i64 ┆ str    ┆ i64       ┆ str            ┆ i64       ┆ str          │
        ╞═════╪════════╪═════╪════════╪═══════════╪════════════════╪═══════════╪══════════════╡
        │ 1   ┆ banana ┆ 5   ┆ beetle ┆ 5         ┆ banana         ┆ 1         ┆ beetle       │
        │ 2   ┆ banana ┆ 4   ┆ audi   ┆ 4         ┆ apple          ┆ 2         ┆ beetle       │
        │ 3   ┆ apple  ┆ 3   ┆ beetle ┆ 3         ┆ apple          ┆ 3         ┆ beetle       │
        │ 4   ┆ apple  ┆ 2   ┆ beetle ┆ 2         ┆ banana         ┆ 4         ┆ audi         │
        │ 5   ┆ banana ┆ 1   ┆ beetle ┆ 1         ┆ banana         ┆ 5         ┆ beetle       │
        └─────┴────────┴─────┴────────┴───────────┴────────────────┴───────────┴──────────────┘
        """  # noqa: W505
        return wrap_expr(self._pyexpr.reverse())

    def std(self, ddof: int = 1) -> Expr:
        """
        Get standard deviation.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 0, 1]})
        >>> df.select(pl.col("a").std())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.std(ddof))

    def var(self, ddof: int = 1) -> Expr:
        """
        Get variance.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 0, 1]})
        >>> df.select(pl.col("a").var())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.var(ddof))

    def max(self) -> Expr:
        """
        Get maximum value.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1.0, float("nan"), 1.0]})
        >>> df.select(pl.col("a").max())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.max())

    def min(self) -> Expr:
        """
        Get minimum value.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1.0, float("nan"), 1.0]})
        >>> df.select(pl.col("a").min())
        shape: (1, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ -1.0 │
        └──────┘
        """
        return wrap_expr(self._pyexpr.min())

    def nan_max(self) -> Expr:
        """
        Get maximum value, but propagate/poison encountered NaN values.

        This differs from numpy's `nanmax` as numpy defaults to propagating NaN values,
        whereas polars defaults to ignoring them.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.0, float("nan")]})
        >>> df.select(pl.col("a").nan_max())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ NaN │
        └─────┘
        """
        return wrap_expr(self._pyexpr.nan_max())

    def nan_min(self) -> Expr:
        """
        Get minimum value, but propagate/poison encountered NaN values.

        This differs from numpy's `nanmax` as numpy defaults to propagating NaN values,
        whereas polars defaults to ignoring them.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.0, float("nan")]})
        >>> df.select(pl.col("a").nan_min())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ NaN │
        └─────┘
        """
        return wrap_expr(self._pyexpr.nan_min())

    def sum(self) -> Expr:
        """
        Get sum value.

        Notes
        -----
        * Dtypes in {Int8, UInt8, Int16, UInt16} are cast to
          Int64 before summing to prevent overflow issues.
        * If there are no non-null values, then the output is `0`.
          If you would prefer empty sums to return `None`, you can
          use `pl.when(expr.count()>0).then(expr.sum())` instead
          of `expr.sum()`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 0, 1]})
        >>> df.select(pl.col("a").sum())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │  0  │
        └─────┘
        """
        return wrap_expr(self._pyexpr.sum())

    def mean(self) -> Expr:
        """
        Get mean value.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 0, 1]})
        >>> df.select(pl.col("a").mean())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.mean())

    def median(self) -> Expr:
        """
        Get median value using linear interpolation.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 0, 1]})
        >>> df.select(pl.col("a").median())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.median())

    def product(self) -> Expr:
        """
        Compute the product of an expression.

        Notes
        -----
        If there are no non-null values, then the output is `1`.
        If you would prefer empty products to return `None`, you can
        use `pl.when(expr.count()>0).then(expr.product())` instead
        of `expr.product()`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").product())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 6   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.product())

    def n_unique(self) -> Expr:
        """
        Count unique values.

        Notes
        -----
        `null` is considered to be a unique value for the purposes of this operation.

        Examples
        --------
        >>> df = pl.DataFrame({"x": [1, 1, 2, 2, 3], "y": [1, 1, 1, None, None]})
        >>> df.select(
        ...     x_unique=pl.col("x").n_unique(),
        ...     y_unique=pl.col("y").n_unique(),
        ... )
        shape: (1, 2)
        ┌──────────┬──────────┐
        │ x_unique ┆ y_unique │
        │ ---      ┆ ---      │
        │ u32      ┆ u32      │
        ╞══════════╪══════════╡
        │ 3        ┆ 2        │
        └──────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.n_unique())

    def approx_n_unique(self) -> Expr:
        """
        Approximate count of unique values.

        This is done using the HyperLogLog++ algorithm for cardinality estimation.

        Examples
        --------
        >>> df = pl.DataFrame({"n": [1, 1, 2]})
        >>> df.select(pl.col("n").approx_n_unique())
        shape: (1, 1)
        ┌─────┐
        │ n   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 2   │
        └─────┘
        >>> df = pl.DataFrame({"n": range(1000)})
        >>> df.select(
        ...     exact=pl.col("n").n_unique(),
        ...     approx=pl.col("n").approx_n_unique(),
        ... )  # doctest: +SKIP
        shape: (1, 2)
        ┌───────┬────────┐
        │ exact ┆ approx │
        │ ---   ┆ ---    │
        │ u32   ┆ u32    │
        ╞═══════╪════════╡
        │ 1000  ┆ 1005   │
        └───────┴────────┘
        """
        return wrap_expr(self._pyexpr.approx_n_unique())

    def null_count(self) -> Expr:
        """
        Count null values.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [None, 1, None],
        ...         "b": [10, None, 300],
        ...         "c": [350, 650, 850],
        ...     }
        ... )
        >>> df.select(pl.all().null_count())
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ u32 ┆ u32 ┆ u32 │
        ╞═════╪═════╪═════╡
        │ 2   ┆ 1   ┆ 0   │
        └─────┴─────┴─────┘
        """
        return wrap_expr(self._pyexpr.null_count())

    def has_nulls(self) -> Expr:
        """
        Check whether the expression contains one or more null values.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [None, 1, None],
        ...         "b": [10, None, 300],
        ...         "c": [350, 650, 850],
        ...     }
        ... )
        >>> df.select(pl.all().has_nulls())
        shape: (1, 3)
        ┌──────┬──────┬───────┐
        │ a    ┆ b    ┆ c     │
        │ ---  ┆ ---  ┆ ---   │
        │ bool ┆ bool ┆ bool  │
        ╞══════╪══════╪═══════╡
        │ true ┆ true ┆ false │
        └──────┴──────┴───────┘
        """
        return self.null_count() > 0

    def arg_unique(self) -> Expr:
        """
        Get index of first unique value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [8, 9, 10],
        ...         "b": [None, 4, 4],
        ...     }
        ... )
        >>> df.select(pl.col("a").arg_unique())
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 0   │
        │ 1   │
        │ 2   │
        └─────┘
        >>> df.select(pl.col("b").arg_unique())
        shape: (2, 1)
        ┌─────┐
        │ b   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 0   │
        │ 1   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arg_unique())

    def unique(self, *, maintain_order: bool = False) -> Expr:
        """
        Get unique values of this expression.

        Parameters
        ----------
        maintain_order
            Maintain order of data. This requires more work.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2]})
        >>> df.select(pl.col("a").unique())  # doctest: +IGNORE_RESULT
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        │ 1   │
        └─────┘
        >>> df.select(pl.col("a").unique(maintain_order=True))
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
        if maintain_order:
            return wrap_expr(self._pyexpr.unique_stable())
        return wrap_expr(self._pyexpr.unique())

    def first(self) -> Expr:
        """
        Get the first value.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2]})
        >>> df.select(pl.col("a").first())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.first())

    def last(self) -> Expr:
        """
        Get the last value.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 3, 2]})
        >>> df.select(pl.col("a").last())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.last())

    @unstable()
    def item(self) -> Expr:
        """
        Get the single value.

        This raises an error if there is not exactly one value.

        See Also
        --------
        :meth:`Expr.get` : Get a single value by index.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1]})
        >>> df.select(pl.col("a").item())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        └─────┘
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").item())
        Traceback (most recent call last):
        ...
        polars.exceptions.ComputeError: aggregation 'item' expected a single value, got 3 values
        """  # noqa: W505
        return wrap_expr(self._pyexpr.item())

    def over(
        self,
        partition_by: IntoExpr | Iterable[IntoExpr] | None = None,
        *more_exprs: IntoExpr,
        order_by: IntoExpr | Iterable[IntoExpr] | None = None,
        descending: bool = False,
        nulls_last: bool = False,
        mapping_strategy: WindowMappingStrategy = "group_to_rows",
    ) -> Expr:
        """
        Compute expressions over the given groups.

        This expression is similar to performing a group by aggregation and joining the
        result back into the original DataFrame.

        The outcome is similar to how `window functions
        <https://www.postgresql.org/docs/current/tutorial-window.html>`_
        work in PostgreSQL.

        Parameters
        ----------
        partition_by
            Column(s) to group by. Accepts expression input. Strings are parsed as
            column names.
        *more_exprs
            Additional columns to group by, specified as positional arguments.
        order_by
            Order the window functions/aggregations with the partitioned groups by the
            result of the expression passed to `order_by`.
        descending
            In case 'order_by' is given, indicate whether to order in
            ascending or descending order.
        nulls_last
            In case 'order_by' is given, indicate whether to order
            the nulls in last position.
        mapping_strategy: {'group_to_rows', 'join', 'explode'}
            - group_to_rows
                If the aggregation results in multiple values per group, map them back
                to their row position in the DataFrame. This can only be done if each
                group yields the same elements before aggregation as after. If the
                aggregation results in one scalar value per group, this value will be
                mapped to every row.
            - join
                If the aggregation may result in multiple values per group, join the
                values as 'List<group_dtype>' to each row position. Warning: this can be
                memory intensive. If the aggregation always results in one scalar value
                per group, join this value as '<group_dtype>' to each row position.
            - explode
                If the aggregation may result in multiple values per group, map each
                value to a new row, similar to the results of `group_by` + `agg` +
                `explode`. If the aggregation always results in one scalar value per
                group, map this value to one row position. Sorting of the given groups
                is required if the groups are not part of the window operation for the
                operation, otherwise the result would not make sense. This operation
                changes the number of rows.

        Examples
        --------
        Pass the name of a column to compute the expression over that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["a", "a", "b", "b", "b"],
        ...         "b": [1, 2, 3, 5, 3],
        ...         "c": [5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> df.with_columns(c_max=pl.col("c").max().over("a"))
        shape: (5, 4)
        ┌─────┬─────┬─────┬───────┐
        │ a   ┆ b   ┆ c   ┆ c_max │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ str ┆ i64 ┆ i64 ┆ i64   │
        ╞═════╪═════╪═════╪═══════╡
        │ a   ┆ 1   ┆ 5   ┆ 5     │
        │ a   ┆ 2   ┆ 4   ┆ 5     │
        │ b   ┆ 3   ┆ 3   ┆ 3     │
        │ b   ┆ 5   ┆ 2   ┆ 3     │
        │ b   ┆ 3   ┆ 1   ┆ 3     │
        └─────┴─────┴─────┴───────┘

        Expression input is also supported.

        >>> df.with_columns(c_max=pl.col("c").max().over(pl.col("b") // 2))
        shape: (5, 4)
        ┌─────┬─────┬─────┬───────┐
        │ a   ┆ b   ┆ c   ┆ c_max │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ str ┆ i64 ┆ i64 ┆ i64   │
        ╞═════╪═════╪═════╪═══════╡
        │ a   ┆ 1   ┆ 5   ┆ 5     │
        │ a   ┆ 2   ┆ 4   ┆ 4     │
        │ b   ┆ 3   ┆ 3   ┆ 4     │
        │ b   ┆ 5   ┆ 2   ┆ 2     │
        │ b   ┆ 3   ┆ 1   ┆ 4     │
        └─────┴─────┴─────┴───────┘

        Group by multiple columns by passing multiple column names or expressions.

        >>> df.with_columns(c_min=pl.col("c").min().over("a", pl.col("b") % 2))
        shape: (5, 4)
        ┌─────┬─────┬─────┬───────┐
        │ a   ┆ b   ┆ c   ┆ c_min │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ str ┆ i64 ┆ i64 ┆ i64   │
        ╞═════╪═════╪═════╪═══════╡
        │ a   ┆ 1   ┆ 5   ┆ 5     │
        │ a   ┆ 2   ┆ 4   ┆ 4     │
        │ b   ┆ 3   ┆ 3   ┆ 1     │
        │ b   ┆ 5   ┆ 2   ┆ 1     │
        │ b   ┆ 3   ┆ 1   ┆ 1     │
        └─────┴─────┴─────┴───────┘

        Mapping strategy `join` joins the values by group.

        >>> df.with_columns(
        ...     c_pairs=pl.col("c").head(2).over("a", mapping_strategy="join")
        ... )
        shape: (5, 4)
        ┌─────┬─────┬─────┬───────────┐
        │ a   ┆ b   ┆ c   ┆ c_pairs   │
        │ --- ┆ --- ┆ --- ┆ ---       │
        │ str ┆ i64 ┆ i64 ┆ list[i64] │
        ╞═════╪═════╪═════╪═══════════╡
        │ a   ┆ 1   ┆ 5   ┆ [5, 4]    │
        │ a   ┆ 2   ┆ 4   ┆ [5, 4]    │
        │ b   ┆ 3   ┆ 3   ┆ [3, 2]    │
        │ b   ┆ 5   ┆ 2   ┆ [3, 2]    │
        │ b   ┆ 3   ┆ 1   ┆ [3, 2]    │
        └─────┴─────┴─────┴───────────┘

        Mapping strategy `explode` maps the values to new rows, changing the shape.

        >>> df.select(
        ...     c_first_2=pl.col("c").head(2).over("a", mapping_strategy="explode")
        ... )
        shape: (4, 1)
        ┌───────────┐
        │ c_first_2 │
        │ ---       │
        │ i64       │
        ╞═══════════╡
        │ 5         │
        │ 4         │
        │ 3         │
        │ 2         │
        └───────────┘

        You can use non-elementwise expressions with `over` too. By default they are
        evaluated using row-order, but you can specify a different one using `order_by`.

        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "store_id": ["a", "a", "b", "b"],
        ...         "date": [
        ...             date(2024, 9, 18),
        ...             date(2024, 9, 17),
        ...             date(2024, 9, 18),
        ...             date(2024, 9, 16),
        ...         ],
        ...         "sales": [7, 9, 8, 10],
        ...     }
        ... )
        >>> df.with_columns(
        ...     cumulative_sales=pl.col("sales")
        ...     .cum_sum()
        ...     .over("store_id", order_by="date")
        ... )
        shape: (4, 4)
        ┌──────────┬────────────┬───────┬──────────────────┐
        │ store_id ┆ date       ┆ sales ┆ cumulative_sales │
        │ ---      ┆ ---        ┆ ---   ┆ ---              │
        │ str      ┆ date       ┆ i64   ┆ i64              │
        ╞══════════╪════════════╪═══════╪══════════════════╡
        │ a        ┆ 2024-09-18 ┆ 7     ┆ 16               │
        │ a        ┆ 2024-09-17 ┆ 9     ┆ 9                │
        │ b        ┆ 2024-09-18 ┆ 8     ┆ 18               │
        │ b        ┆ 2024-09-16 ┆ 10    ┆ 10               │
        └──────────┴────────────┴───────┴──────────────────┘

        If you don't require that the group order be preserved, then the more performant
        option is to use `mapping_strategy='explode'` - be careful however to only ever
        use this in a `select` statement, not a `with_columns` one.

        >>> window = {
        ...     "partition_by": "store_id",
        ...     "order_by": "date",
        ...     "mapping_strategy": "explode",
        ... }
        >>> df.select(
        ...     pl.all().over(**window),
        ...     cumulative_sales=pl.col("sales").cum_sum().over(**window),
        ... )
        shape: (4, 4)
        ┌──────────┬────────────┬───────┬──────────────────┐
        │ store_id ┆ date       ┆ sales ┆ cumulative_sales │
        │ ---      ┆ ---        ┆ ---   ┆ ---              │
        │ str      ┆ date       ┆ i64   ┆ i64              │
        ╞══════════╪════════════╪═══════╪══════════════════╡
        │ a        ┆ 2024-09-17 ┆ 9     ┆ 9                │
        │ a        ┆ 2024-09-18 ┆ 7     ┆ 16               │
        │ b        ┆ 2024-09-16 ┆ 10    ┆ 10               │
        │ b        ┆ 2024-09-18 ┆ 8     ┆ 18               │
        └──────────┴────────────┴───────┴──────────────────┘
        """
        if partition_by is not None:
            partition_by_pyexprs = parse_into_list_of_expressions(
                partition_by, *more_exprs
            )
        else:
            partition_by_pyexprs = None

        if order_by is not None:
            order_by_pyexprs = parse_into_list_of_expressions(order_by)
        else:
            order_by_pyexprs = None

        return wrap_expr(
            self._pyexpr.over(
                partition_by_pyexprs,
                order_by=order_by_pyexprs,
                order_by_descending=descending,
                order_by_nulls_last=False,  # does not work yet
                mapping_strategy=mapping_strategy,
            )
        )

    def rolling(
        self,
        index_column: str,
        *,
        period: str | timedelta,
        offset: str | timedelta | None = None,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Create rolling groups based on a temporal or integer column.

        If you have a time series `<t_0, t_1, ..., t_n>`, then by default the
        windows created will be

            * (t_0 - period, t_0]
            * (t_1 - period, t_1]
            * ...
            * (t_n - period, t_n]

        whereas if you pass a non-default `offset`, then the windows will be

            * (t_0 + offset, t_0 + offset + period]
            * (t_1 + offset, t_1 + offset + period]
            * ...
            * (t_n + offset, t_n + offset + period]

        The `period` and `offset` arguments are created either from a timedelta, or
        by using the following string language:

        - 1ns   (1 nanosecond)
        - 1us   (1 microsecond)
        - 1ms   (1 millisecond)
        - 1s    (1 second)
        - 1m    (1 minute)
        - 1h    (1 hour)
        - 1d    (1 calendar day)
        - 1w    (1 calendar week)
        - 1mo   (1 calendar month)
        - 1q    (1 calendar quarter)
        - 1y    (1 calendar year)
        - 1i    (1 index count)

        Or combine them:
        "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

        By "calendar day", we mean the corresponding time on the next day (which may
        not be 24 hours, due to daylight savings). Similarly for "calendar week",
        "calendar month", "calendar quarter", and "calendar year".

        Parameters
        ----------
        index_column
            Column used to group based on the time window.
            Often of type Date/Datetime.
            This column must be sorted in ascending order.
            In case of a rolling group by on indices, dtype needs to be one of
            {UInt32, UInt64, Int32, Int64}. Note that the first three get temporarily
            cast to Int64, so if performance matters use an Int64 column.
        period
            Length of the window - must be non-negative.
        offset
            Offset of the window. Default is `-period`.
        closed : {'right', 'left', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive).

        Examples
        --------
        >>> dates = [
        ...     "2020-01-01 13:45:48",
        ...     "2020-01-01 16:42:13",
        ...     "2020-01-01 16:45:09",
        ...     "2020-01-02 18:12:48",
        ...     "2020-01-03 19:45:32",
        ...     "2020-01-08 23:16:43",
        ... ]
        >>> df = pl.DataFrame({"dt": dates, "a": [3, 7, 5, 9, 2, 1]}).with_columns(
        ...     pl.col("dt").str.strptime(pl.Datetime).set_sorted()
        ... )
        >>> df.with_columns(
        ...     sum_a=pl.sum("a").rolling(index_column="dt", period="2d"),
        ...     min_a=pl.min("a").rolling(index_column="dt", period="2d"),
        ...     max_a=pl.max("a").rolling(index_column="dt", period="2d"),
        ... )
        shape: (6, 5)
        ┌─────────────────────┬─────┬───────┬───────┬───────┐
        │ dt                  ┆ a   ┆ sum_a ┆ min_a ┆ max_a │
        │ ---                 ┆ --- ┆ ---   ┆ ---   ┆ ---   │
        │ datetime[μs]        ┆ i64 ┆ i64   ┆ i64   ┆ i64   │
        ╞═════════════════════╪═════╪═══════╪═══════╪═══════╡
        │ 2020-01-01 13:45:48 ┆ 3   ┆ 3     ┆ 3     ┆ 3     │
        │ 2020-01-01 16:42:13 ┆ 7   ┆ 10    ┆ 3     ┆ 7     │
        │ 2020-01-01 16:45:09 ┆ 5   ┆ 15    ┆ 3     ┆ 7     │
        │ 2020-01-02 18:12:48 ┆ 9   ┆ 24    ┆ 3     ┆ 9     │
        │ 2020-01-03 19:45:32 ┆ 2   ┆ 11    ┆ 2     ┆ 9     │
        │ 2020-01-08 23:16:43 ┆ 1   ┆ 1     ┆ 1     ┆ 1     │
        └─────────────────────┴─────┴───────┴───────┴───────┘
        """
        if offset is None:
            offset = negate_duration_string(parse_as_duration_string(period))

        period = parse_as_duration_string(period)
        offset = parse_as_duration_string(offset)

        return wrap_expr(self._pyexpr.rolling(index_column, period, offset, closed))

    def is_unique(self) -> Expr:
        """
        Get mask of unique values.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2]})
        >>> df.select(pl.col("a").is_unique())
        shape: (3, 1)
        ┌───────┐
        │ a     │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ false │
        │ false │
        │ true  │
        └───────┘
        """
        return wrap_expr(self._pyexpr.is_unique())

    def is_first_distinct(self) -> Expr:
        """
        Return a boolean mask indicating the first occurrence of each distinct value.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2, 3, 2]})
        >>> df.with_columns(pl.col("a").is_first_distinct().alias("first"))
        shape: (5, 2)
        ┌─────┬───────┐
        │ a   ┆ first │
        │ --- ┆ ---   │
        │ i64 ┆ bool  │
        ╞═════╪═══════╡
        │ 1   ┆ true  │
        │ 1   ┆ false │
        │ 2   ┆ true  │
        │ 3   ┆ true  │
        │ 2   ┆ false │
        └─────┴───────┘
        """
        return wrap_expr(self._pyexpr.is_first_distinct())

    def is_last_distinct(self) -> Expr:
        """
        Return a boolean mask indicating the last occurrence of each distinct value.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2, 3, 2]})
        >>> df.with_columns(pl.col("a").is_last_distinct().alias("last"))
        shape: (5, 2)
        ┌─────┬───────┐
        │ a   ┆ last  │
        │ --- ┆ ---   │
        │ i64 ┆ bool  │
        ╞═════╪═══════╡
        │ 1   ┆ false │
        │ 1   ┆ true  │
        │ 2   ┆ false │
        │ 3   ┆ true  │
        │ 2   ┆ true  │
        └─────┴───────┘
        """
        return wrap_expr(self._pyexpr.is_last_distinct())

    def is_duplicated(self) -> Expr:
        """
        Return a boolean mask indicating duplicated values.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2]})
        >>> df.select(pl.col("a").is_duplicated())
        shape: (3, 1)
        ┌───────┐
        │ a     │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ true  │
        │ true  │
        │ false │
        └───────┘
        """
        return wrap_expr(self._pyexpr.is_duplicated())

    def peak_max(self) -> Expr:
        """
        Get a boolean mask of the local maximum peaks.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 4, 5]})
        >>> df.select(pl.col("a").peak_max())
        shape: (5, 1)
        ┌───────┐
        │ a     │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ false │
        │ false │
        │ false │
        │ false │
        │ true  │
        └───────┘
        """
        return wrap_expr(self._pyexpr.peak_max())

    def peak_min(self) -> Expr:
        """
        Get a boolean mask of the local minimum peaks.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [4, 1, 3, 2, 5]})
        >>> df.select(pl.col("a").peak_min())
        shape: (5, 1)
        ┌───────┐
        │ a     │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ false │
        │ true  │
        │ false │
        │ true  │
        │ false │
        └───────┘
        """
        return wrap_expr(self._pyexpr.peak_min())

    def quantile(
        self,
        quantile: float | Expr,
        interpolation: QuantileMethod = "nearest",
    ) -> Expr:
        """
        Get quantile value.

        Parameters
        ----------
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0, 1, 2, 3, 4, 5]})
        >>> df.select(pl.col("a").quantile(0.3))
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 2.0 │
        └─────┘
        >>> df.select(pl.col("a").quantile(0.3, interpolation="higher"))
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 2.0 │
        └─────┘
        >>> df.select(pl.col("a").quantile(0.3, interpolation="lower"))
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        └─────┘
        >>> df.select(pl.col("a").quantile(0.3, interpolation="midpoint"))
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.5 │
        └─────┘
        >>> df.select(pl.col("a").quantile(0.3, interpolation="linear"))
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.5 │
        └─────┘
        """  # noqa: W505
        quantile_pyexpr = parse_into_expression(quantile)
        return wrap_expr(self._pyexpr.quantile(quantile_pyexpr, interpolation))

    @unstable()
    def cut(
        self,
        breaks: Sequence[float],
        *,
        labels: Sequence[str] | None = None,
        left_closed: bool = False,
        include_breaks: bool = False,
    ) -> Expr:
        """
        Bin continuous values into discrete categories.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        breaks
            List of unique cut points.
        labels
            Names of the categories. The number of labels must be equal to the number
            of cut points plus one.
        left_closed
            Set the intervals to be left-closed instead of right-closed.
        include_breaks
            Include a column with the right endpoint of the bin each observation falls
            in. This will change the data type of the output from a
            :class:`Categorical` to a :class:`Struct`.

        Returns
        -------
        Expr
            Expression of data type :class:`Categorical` if `include_breaks` is set to
            `False` (default), otherwise an expression of data type :class:`Struct`.

        See Also
        --------
        qcut

        Examples
        --------
        Divide a column into three categories.

        >>> df = pl.DataFrame({"foo": [-2, -1, 0, 1, 2]})
        >>> df.with_columns(
        ...     pl.col("foo").cut([-1, 1], labels=["a", "b", "c"]).alias("cut")
        ... )
        shape: (5, 2)
        ┌─────┬─────┐
        │ foo ┆ cut │
        │ --- ┆ --- │
        │ i64 ┆ cat │
        ╞═════╪═════╡
        │ -2  ┆ a   │
        │ -1  ┆ a   │
        │ 0   ┆ b   │
        │ 1   ┆ b   │
        │ 2   ┆ c   │
        └─────┴─────┘

        Add both the category and the breakpoint.

        >>> df.with_columns(
        ...     pl.col("foo").cut([-1, 1], include_breaks=True).alias("cut")
        ... ).unnest("cut")
        shape: (5, 3)
        ┌─────┬────────────┬────────────┐
        │ foo ┆ breakpoint ┆ category   │
        │ --- ┆ ---        ┆ ---        │
        │ i64 ┆ f64        ┆ cat        │
        ╞═════╪════════════╪════════════╡
        │ -2  ┆ -1.0       ┆ (-inf, -1] │
        │ -1  ┆ -1.0       ┆ (-inf, -1] │
        │ 0   ┆ 1.0        ┆ (-1, 1]    │
        │ 1   ┆ 1.0        ┆ (-1, 1]    │
        │ 2   ┆ inf        ┆ (1, inf]   │
        └─────┴────────────┴────────────┘
        """
        return wrap_expr(self._pyexpr.cut(breaks, labels, left_closed, include_breaks))

    @unstable()
    def qcut(
        self,
        quantiles: Sequence[float] | int,
        *,
        labels: Sequence[str] | None = None,
        left_closed: bool = False,
        allow_duplicates: bool = False,
        include_breaks: bool = False,
    ) -> Expr:
        """
        Bin continuous values into discrete categories based on their quantiles.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        quantiles
            Either a list of quantile probabilities between 0 and 1 or a positive
            integer determining the number of bins with uniform probability.
        labels
            Names of the categories. The number of labels must be equal to the number
            of categories.
        left_closed
            Set the intervals to be left-closed instead of right-closed.
        allow_duplicates
            If set to `True`, duplicates in the resulting quantiles are dropped,
            rather than raising a `DuplicateError`. This can happen even with unique
            probabilities, depending on the data.
        include_breaks
            Include a column with the right endpoint of the bin each observation falls
            in. This will change the data type of the output from a
            :class:`Categorical` to a :class:`Struct`.

        Returns
        -------
        Expr
            Expression of data type :class:`Categorical` if `include_breaks` is set to
            `False` (default), otherwise an expression of data type :class:`Struct`.

        See Also
        --------
        cut

        Examples
        --------
        Divide a column into three categories according to pre-defined quantile
        probabilities.

        >>> df = pl.DataFrame({"foo": [-2, -1, 0, 1, 2]})
        >>> df.with_columns(
        ...     pl.col("foo").qcut([0.25, 0.75], labels=["a", "b", "c"]).alias("qcut")
        ... )
        shape: (5, 2)
        ┌─────┬──────┐
        │ foo ┆ qcut │
        │ --- ┆ ---  │
        │ i64 ┆ cat  │
        ╞═════╪══════╡
        │ -2  ┆ a    │
        │ -1  ┆ a    │
        │ 0   ┆ b    │
        │ 1   ┆ b    │
        │ 2   ┆ c    │
        └─────┴──────┘

        Divide a column into two categories using uniform quantile probabilities.

        >>> df.with_columns(
        ...     pl.col("foo")
        ...     .qcut(2, labels=["low", "high"], left_closed=True)
        ...     .alias("qcut")
        ... )
        shape: (5, 2)
        ┌─────┬──────┐
        │ foo ┆ qcut │
        │ --- ┆ ---  │
        │ i64 ┆ cat  │
        ╞═════╪══════╡
        │ -2  ┆ low  │
        │ -1  ┆ low  │
        │ 0   ┆ high │
        │ 1   ┆ high │
        │ 2   ┆ high │
        └─────┴──────┘

        Add both the category and the breakpoint.

        >>> df.with_columns(
        ...     pl.col("foo").qcut([0.25, 0.75], include_breaks=True).alias("qcut")
        ... ).unnest("qcut")
        shape: (5, 3)
        ┌─────┬────────────┬────────────┐
        │ foo ┆ breakpoint ┆ category   │
        │ --- ┆ ---        ┆ ---        │
        │ i64 ┆ f64        ┆ cat        │
        ╞═════╪════════════╪════════════╡
        │ -2  ┆ -1.0       ┆ (-inf, -1] │
        │ -1  ┆ -1.0       ┆ (-inf, -1] │
        │ 0   ┆ 1.0        ┆ (-1, 1]    │
        │ 1   ┆ 1.0        ┆ (-1, 1]    │
        │ 2   ┆ inf        ┆ (1, inf]   │
        └─────┴────────────┴────────────┘
        """
        if isinstance(quantiles, int):
            pyexpr = self._pyexpr.qcut_uniform(
                quantiles, labels, left_closed, allow_duplicates, include_breaks
            )
        else:
            pyexpr = self._pyexpr.qcut(
                quantiles, labels, left_closed, allow_duplicates, include_breaks
            )

        return wrap_expr(pyexpr)

    def rle(self) -> Expr:
        """
        Compress the column data using run-length encoding.

        Run-length encoding (RLE) encodes data by storing each *run* of identical values
        as a single value and its length.

        Returns
        -------
        Expr
            Expression of data type `Struct` with fields `len` of data type `UInt32`
            and `value` of the original data type.

        See Also
        --------
        rle_id

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 2, 1, None, 1, 3, 3]})
        >>> df.select(pl.col("a").rle()).unnest("a")
        shape: (6, 2)
        ┌─────┬───────┐
        │ len ┆ value │
        │ --- ┆ ---   │
        │ u32 ┆ i64   │
        ╞═════╪═══════╡
        │ 2   ┆ 1     │
        │ 1   ┆ 2     │
        │ 1   ┆ 1     │
        │ 1   ┆ null  │
        │ 1   ┆ 1     │
        │ 2   ┆ 3     │
        └─────┴───────┘
        """
        return wrap_expr(self._pyexpr.rle())

    def rle_id(self) -> Expr:
        """
        Get a distinct integer ID for each run of identical values.

        The ID starts at 0 and increases by one each time the value of the column
        changes.

        Returns
        -------
        Expr
            Expression of data type `UInt32`.

        See Also
        --------
        rle

        Notes
        -----
        This functionality is especially useful for defining a new group for every time
        a column's value changes, rather than for every distinct value of that column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 1, 1, 1],
        ...         "b": ["x", "x", None, "y", "y"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     rle_id_a=pl.col("a").rle_id(),
        ...     rle_id_ab=pl.struct("a", "b").rle_id(),
        ... )
        shape: (5, 4)
        ┌─────┬──────┬──────────┬───────────┐
        │ a   ┆ b    ┆ rle_id_a ┆ rle_id_ab │
        │ --- ┆ ---  ┆ ---      ┆ ---       │
        │ i64 ┆ str  ┆ u32      ┆ u32       │
        ╞═════╪══════╪══════════╪═══════════╡
        │ 1   ┆ x    ┆ 0        ┆ 0         │
        │ 2   ┆ x    ┆ 1        ┆ 1         │
        │ 1   ┆ null ┆ 2        ┆ 2         │
        │ 1   ┆ y    ┆ 2        ┆ 3         │
        │ 1   ┆ y    ┆ 2        ┆ 3         │
        └─────┴──────┴──────────┴───────────┘
        """
        return wrap_expr(self._pyexpr.rle_id())

    def filter(
        self,
        *predicates: IntoExprColumn | Iterable[IntoExprColumn],
        **constraints: Any,
    ) -> Expr:
        """
        Filter the expression based on one or more predicate expressions.

        The original order of the remaining elements is preserved.

        Elements where the filter does not evaluate to True are discarded, including
        nulls.

        Mostly useful in an aggregation context. If you want to filter on a DataFrame
        level, use `LazyFrame.filter`.

        Parameters
        ----------
        predicates
            Expression(s) that evaluates to a boolean Series.
        constraints
            Column filters; use `name = value` to filter columns by the supplied value.
            Each constraint will behave the same as `pl.col(name).eq(value)`, and
            be implicitly joined with the other filter conditions using `&`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "group_col": ["g1", "g1", "g2"],
        ...         "b": [1, 2, 3],
        ...     }
        ... )
        >>> df.group_by("group_col").agg(
        ...     lt=pl.col("b").filter(pl.col("b") < 2).sum(),
        ...     gte=pl.col("b").filter(pl.col("b") >= 2).sum(),
        ... ).sort("group_col")
        shape: (2, 3)
        ┌───────────┬─────┬─────┐
        │ group_col ┆ lt  ┆ gte │
        │ ---       ┆ --- ┆ --- │
        │ str       ┆ i64 ┆ i64 │
        ╞═══════════╪═════╪═════╡
        │ g1        ┆ 1   ┆ 2   │
        │ g2        ┆ 0   ┆ 3   │
        └───────────┴─────┴─────┘

        Filter expressions can also take constraints as keyword arguments.

        >>> df = pl.DataFrame(
        ...     {
        ...         "key": ["a", "a", "a", "a", "b", "b", "b", "b", "b"],
        ...         "n": [1, 2, 2, 3, 1, 3, 3, 2, 3],
        ...     },
        ... )
        >>> df.group_by("key").agg(
        ...     n_1=pl.col("n").filter(n=1).sum(),
        ...     n_2=pl.col("n").filter(n=2).sum(),
        ...     n_3=pl.col("n").filter(n=3).sum(),
        ... ).sort(by="key")
        shape: (2, 4)
        ┌─────┬─────┬─────┬─────┐
        │ key ┆ n_1 ┆ n_2 ┆ n_3 │
        │ --- ┆ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 4   ┆ 3   │
        │ b   ┆ 1   ┆ 2   ┆ 9   │
        └─────┴─────┴─────┴─────┘
        """
        predicate = parse_predicates_constraints_into_expression(
            *predicates, **constraints
        )
        return wrap_expr(self._pyexpr.filter(predicate))

    @deprecated("`where` is deprecated; use `filter` instead.")
    def where(self, predicate: Expr) -> Expr:
        """
        Filter a single column.

        .. deprecated:: 0.20.4
            Use the :func:`filter` method instead.

        Alias for :func:`filter`.

        Parameters
        ----------
        predicate
            Boolean expression.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "group_col": ["g1", "g1", "g2"],
        ...         "b": [1, 2, 3],
        ...     }
        ... )
        >>> df.group_by("group_col").agg(  # doctest: +SKIP
        ...     [
        ...         pl.col("b").where(pl.col("b") < 2).sum().alias("lt"),
        ...         pl.col("b").where(pl.col("b") >= 2).sum().alias("gte"),
        ...     ]
        ... ).sort("group_col")
        shape: (2, 3)
        ┌───────────┬─────┬─────┐
        │ group_col ┆ lt  ┆ gte │
        │ ---       ┆ --- ┆ --- │
        │ str       ┆ i64 ┆ i64 │
        ╞═══════════╪═════╪═════╡
        │ g1        ┆ 1   ┆ 2   │
        │ g2        ┆ 0   ┆ 3   │
        └───────────┴─────┴─────┘
        """
        return self.filter(predicate)

    def map_batches(
        self,
        function: Callable[[Series], Series | Any],
        return_dtype: PolarsDataType | pl.DataTypeExpr | None = None,
        *,
        agg_list: bool = False,
        is_elementwise: bool = False,
        returns_scalar: bool = False,
    ) -> Expr:
        """
        Apply a custom python function to a whole Series or sequence of Series.

        The output of this custom function is presumed to be either a Series,
        or a NumPy array (in which case it will be automatically converted into
        a Series), or a scalar that will be converted into a Series. If the
        result is a scalar and you want it to stay as a scalar, pass in
        ``returns_scalar=True``. If you want to apply a
        custom function elementwise over single values, see :func:`map_elements`.
        A reasonable use case for `map` functions is transforming the values
        represented by an expression using a third-party library.

        Parameters
        ----------
        function
            Lambda/function to apply.
        return_dtype
            Datatype of the output Series.

            It is recommended to set this whenever possible. If this is `None`, it tries
            to infer the datatype by calling the function with dummy data and looking at
            the output.
        agg_list
            First implode when in a group-by aggregation.

            .. deprecated:: 1.32.0

                Use `expr.implode().map_batches(..)` instead.
        is_elementwise
            Set to true if the operations is elementwise for better performance
            and optimization.

            An elementwise operations has unit or equal length for all inputs
            and can be ran sequentially on slices without results being affected.
        returns_scalar
            If the function returns a scalar, by default it will be wrapped in
            a list in the output, since the assumption is that the function
            always returns something Series-like. If you want to keep the
            result as a scalar, set this argument to True.

        Notes
        -----
        A UDF passed to `map_batches` must be pure, meaning that it cannot modify
        or depend on state other than its arguments. Polars may call the function
        with arbitrary input data.

        See Also
        --------
        map_elements
        replace

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "sine": [0.0, 1.0, 0.0, -1.0],
        ...         "cosine": [1.0, 0.0, -1.0, 0.0],
        ...     }
        ... )
        >>> df.select(
        ...     pl.all().map_batches(
        ...         lambda x: x.to_numpy().argmax(),
        ...         returns_scalar=True,
        ...     )
        ... )
        shape: (1, 2)
        ┌──────┬────────┐
        │ sine ┆ cosine │
        │ ---  ┆ ---    │
        │ i64  ┆ i64    │
        ╞══════╪════════╡
        │ 1    ┆ 0      │
        └──────┴────────┘

        Here's an example of a function that returns a scalar, where we want it
        to stay as a scalar:

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [0, 1, 0, 1],
        ...         "b": [1, 2, 3, 4],
        ...     }
        ... )
        >>> df.group_by("a").agg(
        ...     pl.col("b").map_batches(
        ...         lambda x: x.max(), returns_scalar=True, return_dtype=pl.self_dtype()
        ...     )
        ... )  # doctest: +IGNORE_RESULT
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 4   │
        │ 0   ┆ 3   │
        └─────┴─────┘

        Call a function that takes multiple arguments by creating a `struct` and
        referencing its fields inside the function call.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [5, 1, 0, 3],
        ...         "b": [4, 2, 3, 4],
        ...     }
        ... )
        >>> df.with_columns(
        ...     a_times_b=pl.struct("a", "b").map_batches(
        ...         lambda x: np.multiply(x.struct.field("a"), x.struct.field("b")),
        ...         return_dtype=pl.Int64,
        ...     )
        ... )
        shape: (4, 3)
        ┌─────┬─────┬───────────┐
        │ a   ┆ b   ┆ a_times_b │
        │ --- ┆ --- ┆ ---       │
        │ i64 ┆ i64 ┆ i64       │
        ╞═════╪═════╪═══════════╡
        │ 5   ┆ 4   ┆ 20        │
        │ 1   ┆ 2   ┆ 2         │
        │ 0   ┆ 3   ┆ 0         │
        │ 3   ┆ 4   ┆ 12        │
        └─────┴─────┴───────────┘
        """
        if agg_list:
            msg = f"""using 'agg_list=True' is deprecated and will be removed in 2.0

Consider using {self}.implode() instead"""
            raise DeprecationWarning(msg)
            self = self.implode()

        def _wrap(sl: Sequence[pl.Series], *args: Any, **kwargs: Any) -> pl.Series:
            return function(sl[0], *args, **kwargs)

        return F.map_batches(
            [self],
            _wrap,
            return_dtype,
            is_elementwise=is_elementwise,
            returns_scalar=returns_scalar,
        )

    def map_elements(
        self,
        function: Callable[[Any], Any],
        return_dtype: PolarsDataType | pl.DataTypeExpr | None = None,
        *,
        skip_nulls: bool = True,
        pass_name: bool = False,
        strategy: MapElementsStrategy = "thread_local",
        returns_scalar: bool = False,
    ) -> Expr:
        """
        Map a custom/user-defined function (UDF) to each element of a column.

        .. warning::
            This method is much slower than the native expressions API.
            Only use it if you cannot implement your logic otherwise.

            Suppose that the function is: `x ↦ sqrt(x)`:

            - For mapping elements of a series, consider:
              `pl.col("col_name").sqrt()`.
            - For mapping inner elements of lists, consider:
              `pl.col("col_name").list.eval(pl.element().sqrt())`.
            - For mapping elements of struct fields, consider:
              `pl.col("col_name").struct.field("field_name").sqrt()`.

            If you want to replace the original column or field,
            consider :meth:`.with_columns <polars.DataFrame.with_columns>`
            and :meth:`.with_fields <polars.Expr.struct.with_fields>`.

        Parameters
        ----------
        function
            Lambda/function to map.
        return_dtype
            Datatype of the output Series.

            It is recommended to set this whenever possible. If this is `None`, it tries
            to infer the datatype by calling the function with dummy data and looking at
            the output.
        skip_nulls
            Don't map the function over values that contain nulls (this is faster).
        pass_name
            Pass the Series name to the custom function (this is more expensive).
        returns_scalar
            .. deprecated:: 1.32.0
                Is ignored and will be removed in 2.0.
        strategy : {'thread_local', 'threading'}
            The threading strategy to use.

            - 'thread_local': run the python function on a single thread.
            - 'threading': run the python function on separate threads. Use with
              care as this can slow performance. This might only speed up
              your code if the amount of work per element is significant
              and the python function releases the GIL (e.g. via calling
              a c function)

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Notes
        -----
        * Using `map_elements` is strongly discouraged as you will be effectively
          running python "for" loops, which will be very slow. Wherever possible you
          should prefer the native expression API to achieve the best performance.

        * If your function is expensive and you don't want it to be called more than
          once for a given input, consider applying an `@lru_cache` decorator to it.
          If your data is suitable you may achieve *significant* speedups.

        * Window function application using `over` is considered a GroupBy context
          here, so `map_elements` can be used to map functions over window groups.

        * A UDF passed to `map_elements` must be pure, meaning that it cannot modify or
          depend on state other than its arguments. Polars may call the function
          with arbitrary input data.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 1],
        ...         "b": ["a", "b", "c", "c"],
        ...     }
        ... )

        The function is applied to each element of column `'a'`:

        >>> df.with_columns(  # doctest: +SKIP
        ...     pl.col("a")
        ...     .map_elements(lambda x: x * 2, return_dtype=pl.self_dtype())
        ...     .alias("a_times_2"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬───────────┐
        │ a   ┆ b   ┆ a_times_2 │
        │ --- ┆ --- ┆ ---       │
        │ i64 ┆ str ┆ i64       │
        ╞═════╪═════╪═══════════╡
        │ 1   ┆ a   ┆ 2         │
        │ 2   ┆ b   ┆ 4         │
        │ 3   ┆ c   ┆ 6         │
        │ 1   ┆ c   ┆ 2         │
        └─────┴─────┴───────────┘

        Tip: it is better to implement this with an expression:

        >>> df.with_columns(
        ...     (pl.col("a") * 2).alias("a_times_2"),
        ... )  # doctest: +IGNORE_RESULT

        >>> (
        ...     df.lazy()
        ...     .group_by("b")
        ...     .agg(
        ...         pl.col("a")
        ...         .implode()
        ...         .map_elements(lambda x: x.sum(), return_dtype=pl.Int64)
        ...     )
        ...     .collect()
        ... )  # doctest: +IGNORE_RESULT
        shape: (3, 2)
        ┌─────┬─────┐
        │ b   ┆ a   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ a   ┆ 1   │
        │ b   ┆ 2   │
        │ c   ┆ 4   │
        └─────┴─────┘

        Tip: again, it is better to implement this with an expression:

        >>> (
        ...     df.lazy()
        ...     .group_by("b", maintain_order=True)
        ...     .agg(pl.col("a").sum())
        ...     .collect()
        ... )  # doctest: +IGNORE_RESULT

        Window function application using `over` will behave as a GroupBy
        context, with your function receiving individual window groups:

        >>> df = pl.DataFrame(
        ...     {
        ...         "key": ["x", "x", "y", "x", "y", "z"],
        ...         "val": [1, 1, 1, 1, 1, 1],
        ...     }
        ... )
        >>> df.with_columns(
        ...     scaled=pl.col("val")
        ...     .implode()
        ...     .map_elements(lambda s: s * len(s), return_dtype=pl.List(pl.Int64))
        ...     .explode()
        ...     .over("key"),
        ... ).sort("key")
        shape: (6, 3)
        ┌─────┬─────┬────────┐
        │ key ┆ val ┆ scaled │
        │ --- ┆ --- ┆ ---    │
        │ str ┆ i64 ┆ i64    │
        ╞═════╪═════╪════════╡
        │ x   ┆ 1   ┆ 3      │
        │ x   ┆ 1   ┆ 3      │
        │ x   ┆ 1   ┆ 3      │
        │ y   ┆ 1   ┆ 2      │
        │ y   ┆ 1   ┆ 2      │
        │ z   ┆ 1   ┆ 1      │
        └─────┴─────┴────────┘

        Note that this function would *also* be better-implemented natively:

        >>> df.with_columns(
        ...     scaled=(pl.col("val") * pl.col("val").count()).over("key"),
        ... ).sort("key")  # doctest: +IGNORE_RESULT

        """
        if strategy == "threading":
            issue_unstable_warning(
                "the 'threading' strategy for `map_elements` is considered unstable."
            )

        # input x: Series of type list containing the group values
        from polars._utils.udfs import warn_on_inefficient_map

        root_names = self.meta.root_names()
        if len(root_names) > 0:
            warn_on_inefficient_map(function, columns=root_names, map_target="expr")

        if pass_name:

            def wrap_f(x: Series, **kwargs: Any) -> Series:  # pragma: no cover
                return_dtype = kwargs["return_dtype"]

                def inner(s: Series | Any) -> Series:  # pragma: no cover
                    if isinstance(s, pl.Series):
                        s = s.alias(x.name)
                    return function(s)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", PolarsInefficientMapWarning)
                    return x.map_elements(
                        inner, return_dtype=return_dtype, skip_nulls=skip_nulls
                    )

        else:

            def wrap_f(x: Series, **kwargs: Any) -> Series:  # pragma: no cover
                return_dtype = kwargs["return_dtype"]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", PolarsInefficientMapWarning)

                    return x.map_elements(
                        function, return_dtype=return_dtype, skip_nulls=skip_nulls
                    )

        if strategy == "thread_local":
            return self.map_batches(
                wrap_f,
                agg_list=False,
                return_dtype=return_dtype,
                returns_scalar=False,
                is_elementwise=True,
            )
        elif strategy == "threading":

            def wrap_threading(x: Series) -> Series:
                def get_lazy_promise(df: DataFrame) -> LazyFrame:
                    return df.lazy().select(
                        F.col("x").map_batches(
                            wrap_f,
                            agg_list=False,
                            return_dtype=return_dtype,
                            returns_scalar=False,
                        )
                    )

                df = x.to_frame("x")

                if x.len() == 0:
                    return get_lazy_promise(df).collect().to_series()

                n_threads = thread_pool_size()
                chunk_size = x.len() // n_threads
                remainder = x.len() % n_threads
                if chunk_size == 0:
                    chunk_sizes = [1 for _ in range(remainder)]
                else:
                    chunk_sizes = [
                        chunk_size + 1 if i < remainder else chunk_size
                        for i in range(n_threads)
                    ]

                # create partitions with LazyFrames
                # these are promises on a computation
                partitions = []
                b = 0
                for step in chunk_sizes:
                    a = b
                    b = b + step
                    partition_df = df[a:b, :]
                    partitions.append(get_lazy_promise(partition_df))

                out = [df.to_series() for df in F.collect_all(partitions)]
                return F.concat(out, rechunk=False)

            return self.map_batches(
                wrap_threading,
                agg_list=False,
                return_dtype=return_dtype,
                returns_scalar=False,
                is_elementwise=True,
            )
        else:
            msg = f"strategy {strategy!r} is not supported"
            raise ValueError(msg)

    def flatten(self) -> Expr:
        """
        Flatten a list or string column.

        Alias for :func:`Expr.list.explode`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "group": ["a", "b", "b"],
        ...         "values": [[1, 2], [2, 3], [4]],
        ...     }
        ... )
        >>> df.group_by("group").agg(pl.col("values").flatten())  # doctest: +SKIP
        shape: (2, 2)
        ┌───────┬───────────┐
        │ group ┆ values    │
        │ ---   ┆ ---       │
        │ str   ┆ list[i64] │
        ╞═══════╪═══════════╡
        │ a     ┆ [1, 2]    │
        │ b     ┆ [2, 3, 4] │
        └───────┴───────────┘
        """
        return wrap_expr(self._pyexpr.explode())

    def explode(self) -> Expr:
        """
        Explode a list expression.

        This means that every item is expanded to a new row.

        Returns
        -------
        Expr
            Expression with the data type of the list elements.

        See Also
        --------
        Expr.list.explode : Explode a list column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "group": ["a", "b"],
        ...         "values": [
        ...             [1, 2],
        ...             [3, 4],
        ...         ],
        ...     }
        ... )
        >>> df.select(pl.col("values").explode())
        shape: (4, 1)
        ┌────────┐
        │ values │
        │ ---    │
        │ i64    │
        ╞════════╡
        │ 1      │
        │ 2      │
        │ 3      │
        │ 4      │
        └────────┘
        """
        return wrap_expr(self._pyexpr.explode())

    def implode(self) -> Expr:
        """
        Aggregate values into a list.

        The returned list itself is a scalar value of `list` dtype.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": [4, 5, 6],
        ...     }
        ... )
        >>> df.select(pl.all().implode())
        shape: (1, 2)
        ┌───────────┬───────────┐
        │ a         ┆ b         │
        │ ---       ┆ ---       │
        │ list[i64] ┆ list[i64] │
        ╞═══════════╪═══════════╡
        │ [1, 2, 3] ┆ [4, 5, 6] │
        └───────────┴───────────┘
        """
        return wrap_expr(self._pyexpr.implode())

    def gather_every(self, n: int, offset: int = 0) -> Expr:
        """
        Take every nth value in the Series and return as a new Series.

        Parameters
        ----------
        n
            Gather every *n*-th row.
        offset
            Starting index.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5, 6, 7, 8, 9]})
        >>> df.select(pl.col("foo").gather_every(3))
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 4   │
        │ 7   │
        └─────┘

        >>> df.select(pl.col("foo").gather_every(3, offset=1))
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        │ 5   │
        │ 8   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.gather_every(n, offset))

    def head(self, n: int | Expr = 10) -> Expr:
        """
        Get the first `n` rows.

        Parameters
        ----------
        n
            Number of rows to return.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5, 6, 7]})
        >>> df.select(pl.col("foo").head(3))
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘
        """
        return self.slice(0, n)

    def tail(self, n: int | Expr = 10) -> Expr:
        """
        Get the last `n` rows.

        Parameters
        ----------
        n
            Number of rows to return.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5, 6, 7]})
        >>> df.select(pl.col("foo").tail(3))
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 5   │
        │ 6   │
        │ 7   │
        └─────┘
        """
        # This cast enables tail with expressions that return unsigned integers,
        # for which negate otherwise raises InvalidOperationError.
        offset = -(
            wrap_expr(parse_into_expression(n)).cast(
                Int64, strict=False, wrap_numerical=True
            )
        )
        return self.slice(offset, n)

    def limit(self, n: int | Expr = 10) -> Expr:
        """
        Get the first `n` rows (alias for :func:`Expr.head`).

        Parameters
        ----------
        n
            Number of rows to return.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5, 6, 7]})
        >>> df.select(pl.col("foo").limit(3))
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘
        """
        return self.head(n)

    def and_(self, *others: Any) -> Expr:
        """
        Method equivalent of bitwise "and" operator `expr & other & ...`.

        Parameters
        ----------
        *others
            One or more integer or boolean expressions to evaluate/combine.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [5, 6, 7, 4, 8],
        ...         "y": [1.5, 2.5, 1.0, 4.0, -5.75],
        ...         "z": [-9, 2, -1, 4, 8],
        ...     }
        ... )
        >>> df.select(
        ...     (pl.col("x") >= pl.col("z"))
        ...     .and_(
        ...         pl.col("y") >= pl.col("z"),
        ...         pl.col("y") == pl.col("y"),
        ...         pl.col("z") <= pl.col("x"),
        ...         pl.col("y") != pl.col("x"),
        ...     )
        ...     .alias("all")
        ... )
        shape: (5, 1)
        ┌───────┐
        │ all   │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ true  │
        │ true  │
        │ true  │
        │ false │
        │ false │
        └───────┘
        """
        return reduce(operator.and_, (self, *others))

    def or_(self, *others: Any) -> Expr:
        """
        Method equivalent of bitwise "or" operator `expr | other | ...`.

        Parameters
        ----------
        *others
            One or more integer or boolean expressions to evaluate/combine.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [5, 6, 7, 4, 8],
        ...         "y": [1.5, 2.5, 1.0, 4.0, -5.75],
        ...         "z": [-9, 2, -1, 4, 8],
        ...     }
        ... )
        >>> df.select(
        ...     (pl.col("x") == pl.col("y"))
        ...     .or_(
        ...         pl.col("x") == pl.col("y"),
        ...         pl.col("y") == pl.col("z"),
        ...         pl.col("y").cast(int) == pl.col("z"),
        ...     )
        ...     .alias("any")
        ... )
        shape: (5, 1)
        ┌───────┐
        │ any   │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ false │
        │ true  │
        │ false │
        │ true  │
        │ false │
        └───────┘
        """
        return reduce(operator.or_, (self,) + others)

    def eq(self, other: Any) -> Expr:
        """
        Method equivalent of equality operator `expr == other`.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [1.0, 2.0, float("nan"), 4.0],
        ...         "y": [2.0, 2.0, float("nan"), 4.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").eq(pl.col("y")).alias("x == y"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬────────┐
        │ x   ┆ y   ┆ x == y │
        │ --- ┆ --- ┆ ---    │
        │ f64 ┆ f64 ┆ bool   │
        ╞═════╪═════╪════════╡
        │ 1.0 ┆ 2.0 ┆ false  │
        │ 2.0 ┆ 2.0 ┆ true   │
        │ NaN ┆ NaN ┆ true   │
        │ 4.0 ┆ 4.0 ┆ true   │
        └─────┴─────┴────────┘
        """
        return self.__eq__(other)

    def eq_missing(self, other: Any) -> Expr:
        """
        Method equivalent of equality operator `expr == other` where `None == None`.

        This differs from default `eq` where null values are propagated.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [1.0, 2.0, float("nan"), 4.0, None, None],
        ...         "y": [2.0, 2.0, float("nan"), 4.0, 5.0, None],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").eq(pl.col("y")).alias("x eq y"),
        ...     pl.col("x").eq_missing(pl.col("y")).alias("x eq_missing y"),
        ... )
        shape: (6, 4)
        ┌──────┬──────┬────────┬────────────────┐
        │ x    ┆ y    ┆ x eq y ┆ x eq_missing y │
        │ ---  ┆ ---  ┆ ---    ┆ ---            │
        │ f64  ┆ f64  ┆ bool   ┆ bool           │
        ╞══════╪══════╪════════╪════════════════╡
        │ 1.0  ┆ 2.0  ┆ false  ┆ false          │
        │ 2.0  ┆ 2.0  ┆ true   ┆ true           │
        │ NaN  ┆ NaN  ┆ true   ┆ true           │
        │ 4.0  ┆ 4.0  ┆ true   ┆ true           │
        │ null ┆ 5.0  ┆ null   ┆ false          │
        │ null ┆ null ┆ null   ┆ true           │
        └──────┴──────┴────────┴────────────────┘
        """
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.eq_missing(other_pyexpr))

    def ge(self, other: Any) -> Expr:
        """
        Method equivalent of "greater than or equal" operator `expr >= other`.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [5.0, 4.0, float("nan"), 2.0],
        ...         "y": [5.0, 3.0, float("nan"), 1.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").ge(pl.col("y")).alias("x >= y"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬────────┐
        │ x   ┆ y   ┆ x >= y │
        │ --- ┆ --- ┆ ---    │
        │ f64 ┆ f64 ┆ bool   │
        ╞═════╪═════╪════════╡
        │ 5.0 ┆ 5.0 ┆ true   │
        │ 4.0 ┆ 3.0 ┆ true   │
        │ NaN ┆ NaN ┆ true   │
        │ 2.0 ┆ 1.0 ┆ true   │
        └─────┴─────┴────────┘
        """
        return self.__ge__(other)

    def gt(self, other: Any) -> Expr:
        """
        Method equivalent of "greater than" operator `expr > other`.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [5.0, 4.0, float("nan"), 2.0],
        ...         "y": [5.0, 3.0, float("nan"), 1.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").gt(pl.col("y")).alias("x > y"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬───────┐
        │ x   ┆ y   ┆ x > y │
        │ --- ┆ --- ┆ ---   │
        │ f64 ┆ f64 ┆ bool  │
        ╞═════╪═════╪═══════╡
        │ 5.0 ┆ 5.0 ┆ false │
        │ 4.0 ┆ 3.0 ┆ true  │
        │ NaN ┆ NaN ┆ false │
        │ 2.0 ┆ 1.0 ┆ true  │
        └─────┴─────┴───────┘
        """
        return self.__gt__(other)

    def le(self, other: Any) -> Expr:
        """
        Method equivalent of "less than or equal" operator `expr <= other`.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [5.0, 4.0, float("nan"), 0.5],
        ...         "y": [5.0, 3.5, float("nan"), 2.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").le(pl.col("y")).alias("x <= y"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬────────┐
        │ x   ┆ y   ┆ x <= y │
        │ --- ┆ --- ┆ ---    │
        │ f64 ┆ f64 ┆ bool   │
        ╞═════╪═════╪════════╡
        │ 5.0 ┆ 5.0 ┆ true   │
        │ 4.0 ┆ 3.5 ┆ false  │
        │ NaN ┆ NaN ┆ true   │
        │ 0.5 ┆ 2.0 ┆ true   │
        └─────┴─────┴────────┘
        """
        return self.__le__(other)

    def lt(self, other: Any) -> Expr:
        """
        Method equivalent of "less than" operator `expr < other`.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [1.0, 2.0, float("nan"), 3.0],
        ...         "y": [2.0, 2.0, float("nan"), 4.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").lt(pl.col("y")).alias("x < y"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬───────┐
        │ x   ┆ y   ┆ x < y │
        │ --- ┆ --- ┆ ---   │
        │ f64 ┆ f64 ┆ bool  │
        ╞═════╪═════╪═══════╡
        │ 1.0 ┆ 2.0 ┆ true  │
        │ 2.0 ┆ 2.0 ┆ false │
        │ NaN ┆ NaN ┆ false │
        │ 3.0 ┆ 4.0 ┆ true  │
        └─────┴─────┴───────┘
        """
        return self.__lt__(other)

    def ne(self, other: Any) -> Expr:
        """
        Method equivalent of inequality operator `expr != other`.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [1.0, 2.0, float("nan"), 4.0],
        ...         "y": [2.0, 2.0, float("nan"), 4.0],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").ne(pl.col("y")).alias("x != y"),
        ... )
        shape: (4, 3)
        ┌─────┬─────┬────────┐
        │ x   ┆ y   ┆ x != y │
        │ --- ┆ --- ┆ ---    │
        │ f64 ┆ f64 ┆ bool   │
        ╞═════╪═════╪════════╡
        │ 1.0 ┆ 2.0 ┆ true   │
        │ 2.0 ┆ 2.0 ┆ false  │
        │ NaN ┆ NaN ┆ false  │
        │ 4.0 ┆ 4.0 ┆ false  │
        └─────┴─────┴────────┘
        """
        return self.__ne__(other)

    def ne_missing(self, other: Any) -> Expr:
        """
        Method equivalent of equality operator `expr != other` where `None == None`.

        This differs from default `ne` where null values are propagated.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "x": [1.0, 2.0, float("nan"), 4.0, None, None],
        ...         "y": [2.0, 2.0, float("nan"), 4.0, 5.0, None],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("x").ne(pl.col("y")).alias("x ne y"),
        ...     pl.col("x").ne_missing(pl.col("y")).alias("x ne_missing y"),
        ... )
        shape: (6, 4)
        ┌──────┬──────┬────────┬────────────────┐
        │ x    ┆ y    ┆ x ne y ┆ x ne_missing y │
        │ ---  ┆ ---  ┆ ---    ┆ ---            │
        │ f64  ┆ f64  ┆ bool   ┆ bool           │
        ╞══════╪══════╪════════╪════════════════╡
        │ 1.0  ┆ 2.0  ┆ true   ┆ true           │
        │ 2.0  ┆ 2.0  ┆ false  ┆ false          │
        │ NaN  ┆ NaN  ┆ false  ┆ false          │
        │ 4.0  ┆ 4.0  ┆ false  ┆ false          │
        │ null ┆ 5.0  ┆ null   ┆ true           │
        │ null ┆ null ┆ null   ┆ false          │
        └──────┴──────┴────────┴────────────────┘
        """
        other_pyexpr = parse_into_expression(other, str_as_lit=True)
        return wrap_expr(self._pyexpr.neq_missing(other_pyexpr))

    def add(self, other: Any) -> Expr:
        """
        Method equivalent of addition operator `expr + other`.

        Parameters
        ----------
        other
            numeric or string value; accepts expression input.

        Examples
        --------
        >>> df = pl.DataFrame({"x": [1, 2, 3, 4, 5]})
        >>> df.with_columns(
        ...     pl.col("x").add(2).alias("x+int"),
        ...     pl.col("x").add(pl.col("x").cum_prod()).alias("x+expr"),
        ... )
        shape: (5, 3)
        ┌─────┬───────┬────────┐
        │ x   ┆ x+int ┆ x+expr │
        │ --- ┆ ---   ┆ ---    │
        │ i64 ┆ i64   ┆ i64    │
        ╞═════╪═══════╪════════╡
        │ 1   ┆ 3     ┆ 2      │
        │ 2   ┆ 4     ┆ 4      │
        │ 3   ┆ 5     ┆ 9      │
        │ 4   ┆ 6     ┆ 28     │
        │ 5   ┆ 7     ┆ 125    │
        └─────┴───────┴────────┘

        >>> df = pl.DataFrame(
        ...     {"x": ["a", "d", "g"], "y": ["b", "e", "h"], "z": ["c", "f", "i"]}
        ... )
        >>> df.with_columns(pl.col("x").add(pl.col("y")).add(pl.col("z")).alias("xyz"))
        shape: (3, 4)
        ┌─────┬─────┬─────┬─────┐
        │ x   ┆ y   ┆ z   ┆ xyz │
        │ --- ┆ --- ┆ --- ┆ --- │
        │ str ┆ str ┆ str ┆ str │
        ╞═════╪═════╪═════╪═════╡
        │ a   ┆ b   ┆ c   ┆ abc │
        │ d   ┆ e   ┆ f   ┆ def │
        │ g   ┆ h   ┆ i   ┆ ghi │
        └─────┴─────┴─────┴─────┘
        """
        return self.__add__(other)

    def floordiv(self, other: Any) -> Expr:
        """
        Method equivalent of integer division operator `expr // other`.

        Parameters
        ----------
        other
            Numeric literal or expression value.

        See Also
        --------
        truediv

        Examples
        --------
        >>> df = pl.DataFrame({"x": [1, 2, 3, 4, 5]})
        >>> df.with_columns(
        ...     pl.col("x").truediv(2).alias("x/2"),
        ...     pl.col("x").floordiv(2).alias("x//2"),
        ... )
        shape: (5, 3)
        ┌─────┬─────┬──────┐
        │ x   ┆ x/2 ┆ x//2 │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ f64 ┆ i64  │
        ╞═════╪═════╪══════╡
        │ 1   ┆ 0.5 ┆ 0    │
        │ 2   ┆ 1.0 ┆ 1    │
        │ 3   ┆ 1.5 ┆ 1    │
        │ 4   ┆ 2.0 ┆ 2    │
        │ 5   ┆ 2.5 ┆ 2    │
        └─────┴─────┴──────┘

        Note that Polars' `floordiv` is subtly different from Python's floor division.
        For example, consider 6.0 floor-divided by 0.1.
        Python gives:

        >>> 6.0 // 0.1
        59.0

        because `0.1` is not represented internally as that exact value,
        but a slightly larger value.
        So the result of the division is slightly less than 60,
        meaning the flooring operation returns 59.0.

        Polars instead first does the floating-point division,
        resulting in a floating-point value of 60.0,
        and then performs the flooring operation using :any:`floor`:

        >>> df = pl.DataFrame({"x": [6.0, 6.03]})
        >>> df.with_columns(
        ...     pl.col("x").truediv(0.1).alias("x/0.1"),
        ... ).with_columns(
        ...     pl.col("x/0.1").floor().alias("x/0.1 floor"),
        ... )
        shape: (2, 3)
        ┌──────┬───────┬─────────────┐
        │ x    ┆ x/0.1 ┆ x/0.1 floor │
        │ ---  ┆ ---   ┆ ---         │
        │ f64  ┆ f64   ┆ f64         │
        ╞══════╪═══════╪═════════════╡
        │ 6.0  ┆ 60.0  ┆ 60.0        │
        │ 6.03 ┆ 60.3  ┆ 60.0        │
        └──────┴───────┴─────────────┘

        yielding the more intuitive result 60.0.
        The row with x = 6.03 is included to demonstrate
        the effect of the flooring operation.

        `floordiv` combines those two steps
        to give the same result with one expression:

        >>> df.with_columns(
        ...     pl.col("x").floordiv(0.1).alias("x//0.1"),
        ... )
        shape: (2, 2)
        ┌──────┬────────┐
        │ x    ┆ x//0.1 │
        │ ---  ┆ ---    │
        │ f64  ┆ f64    │
        ╞══════╪════════╡
        │ 6.0  ┆ 60.0   │
        │ 6.03 ┆ 60.0   │
        └──────┴────────┘
        """
        return self.__floordiv__(other)

    def mod(self, other: Any) -> Expr:
        """
        Method equivalent of modulus operator `expr % other`.

        Parameters
        ----------
        other
            Numeric literal or expression value.

        Examples
        --------
        >>> df = pl.DataFrame({"x": [0, 1, 2, 3, 4]})
        >>> df.with_columns(pl.col("x").mod(2).alias("x%2"))
        shape: (5, 2)
        ┌─────┬─────┐
        │ x   ┆ x%2 │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 0   ┆ 0   │
        │ 1   ┆ 1   │
        │ 2   ┆ 0   │
        │ 3   ┆ 1   │
        │ 4   ┆ 0   │
        └─────┴─────┘
        """
        return self.__mod__(other)

    def mul(self, other: Any) -> Expr:
        """
        Method equivalent of multiplication operator `expr * other`.

        Parameters
        ----------
        other
            Numeric literal or expression value.

        Examples
        --------
        >>> df = pl.DataFrame({"x": [1, 2, 4, 8, 16]})
        >>> df.with_columns(
        ...     pl.col("x").mul(2).alias("x*2"),
        ...     pl.col("x").mul(pl.col("x").log(2)).alias("x * xlog2"),
        ... )
        shape: (5, 3)
        ┌─────┬─────┬───────────┐
        │ x   ┆ x*2 ┆ x * xlog2 │
        │ --- ┆ --- ┆ ---       │
        │ i64 ┆ i64 ┆ f64       │
        ╞═════╪═════╪═══════════╡
        │ 1   ┆ 2   ┆ 0.0       │
        │ 2   ┆ 4   ┆ 2.0       │
        │ 4   ┆ 8   ┆ 8.0       │
        │ 8   ┆ 16  ┆ 24.0      │
        │ 16  ┆ 32  ┆ 64.0      │
        └─────┴─────┴───────────┘
        """
        return self.__mul__(other)

    def sub(self, other: Any) -> Expr:
        """
        Method equivalent of subtraction operator `expr - other`.

        Parameters
        ----------
        other
            Numeric literal or expression value.

        Examples
        --------
        >>> df = pl.DataFrame({"x": [0, 1, 2, 3, 4]})
        >>> df.with_columns(
        ...     pl.col("x").sub(2).alias("x-2"),
        ...     pl.col("x").sub(pl.col("x").cum_sum()).alias("x-expr"),
        ... )
        shape: (5, 3)
        ┌─────┬─────┬────────┐
        │ x   ┆ x-2 ┆ x-expr │
        │ --- ┆ --- ┆ ---    │
        │ i64 ┆ i64 ┆ i64    │
        ╞═════╪═════╪════════╡
        │ 0   ┆ -2  ┆ 0      │
        │ 1   ┆ -1  ┆ 0      │
        │ 2   ┆ 0   ┆ -1     │
        │ 3   ┆ 1   ┆ -3     │
        │ 4   ┆ 2   ┆ -6     │
        └─────┴─────┴────────┘
        """
        return self.__sub__(other)

    def neg(self) -> Expr:
        """
        Method equivalent of unary minus operator `-expr`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 0, 2, None]})
        >>> df.with_columns(pl.col("a").neg())
        shape: (4, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ i64  │
        ╞══════╡
        │ 1    │
        │ 0    │
        │ -2   │
        │ null │
        └──────┘
        """
        return self.__neg__()

    def truediv(self, other: Any) -> Expr:
        """
        Method equivalent of float division operator `expr / other`.

        Parameters
        ----------
        other
            Numeric literal or expression value.

        Notes
        -----
        Zero-division behaviour follows IEEE-754:

        0/0: Invalid operation - mathematically undefined, returns NaN.
        n/0: On finite operands gives an exact infinite result, eg: ±infinity.

        See Also
        --------
        floordiv

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"x": [-2, -1, 0, 1, 2], "y": [0.5, 0.0, 0.0, -4.0, -0.5]}
        ... )
        >>> df.with_columns(
        ...     pl.col("x").truediv(2).alias("x/2"),
        ...     pl.col("x").truediv(pl.col("y")).alias("x/y"),
        ... )
        shape: (5, 4)
        ┌─────┬──────┬──────┬───────┐
        │ x   ┆ y    ┆ x/2  ┆ x/y   │
        │ --- ┆ ---  ┆ ---  ┆ ---   │
        │ i64 ┆ f64  ┆ f64  ┆ f64   │
        ╞═════╪══════╪══════╪═══════╡
        │ -2  ┆ 0.5  ┆ -1.0 ┆ -4.0  │
        │ -1  ┆ 0.0  ┆ -0.5 ┆ -inf  │
        │ 0   ┆ 0.0  ┆ 0.0  ┆ NaN   │
        │ 1   ┆ -4.0 ┆ 0.5  ┆ -0.25 │
        │ 2   ┆ -0.5 ┆ 1.0  ┆ -4.0  │
        └─────┴──────┴──────┴───────┘
        """
        return self.__truediv__(other)

    def pow(self, exponent: IntoExprColumn | int | float) -> Expr:
        """
        Method equivalent of exponentiation operator `expr ** exponent`.

        If the exponent is float, the result follows the dtype of exponent.
        Otherwise, it follows dtype of base.

        Parameters
        ----------
        exponent
            Numeric literal or expression exponent value.

        Examples
        --------
        >>> df = pl.DataFrame({"x": [1, 2, 4, 8]})
        >>> df.with_columns(
        ...     pl.col("x").pow(3).alias("cube"),
        ...     pl.col("x").pow(pl.col("x").log(2)).alias("x ** xlog2"),
        ... )
        shape: (4, 3)
        ┌─────┬──────┬────────────┐
        │ x   ┆ cube ┆ x ** xlog2 │
        │ --- ┆ ---  ┆ ---        │
        │ i64 ┆ i64  ┆ f64        │
        ╞═════╪══════╪════════════╡
        │ 1   ┆ 1    ┆ 1.0        │
        │ 2   ┆ 8    ┆ 2.0        │
        │ 4   ┆ 64   ┆ 16.0       │
        │ 8   ┆ 512  ┆ 512.0      │
        └─────┴──────┴────────────┘

        Raising an integer to a positive integer results in an integer - in order
        to raise to a negative integer, you can cast either the base or the exponent
        to float first:

        >>> df.with_columns(
        ...     x_squared=pl.col("x").pow(2),
        ...     x_inverse=pl.col("x").pow(-1.0),
        ... )
        shape: (4, 3)
        ┌─────┬───────────┬───────────┐
        │ x   ┆ x_squared ┆ x_inverse │
        │ --- ┆ ---       ┆ ---       │
        │ i64 ┆ i64       ┆ f64       │
        ╞═════╪═══════════╪═══════════╡
        │ 1   ┆ 1         ┆ 1.0       │
        │ 2   ┆ 4         ┆ 0.5       │
        │ 4   ┆ 16        ┆ 0.25      │
        │ 8   ┆ 64        ┆ 0.125     │
        └─────┴───────────┴───────────┘
        """
        return self.__pow__(exponent)

    def xor(self, other: Any) -> Expr:
        """
        Method equivalent of bitwise exclusive-or operator `expr ^ other`.

        Parameters
        ----------
        other
            Integer or boolean value; accepts expression input.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"x": [True, False, True, False], "y": [True, True, False, False]}
        ... )
        >>> df.with_columns(pl.col("x").xor(pl.col("y")).alias("x ^ y"))
        shape: (4, 3)
        ┌───────┬───────┬───────┐
        │ x     ┆ y     ┆ x ^ y │
        │ ---   ┆ ---   ┆ ---   │
        │ bool  ┆ bool  ┆ bool  │
        ╞═══════╪═══════╪═══════╡
        │ true  ┆ true  ┆ false │
        │ false ┆ true  ┆ true  │
        │ true  ┆ false ┆ true  │
        │ false ┆ false ┆ false │
        └───────┴───────┴───────┘

        >>> def binary_string(n: int) -> str:
        ...     return bin(n)[2:].zfill(8)
        >>>
        >>> df = pl.DataFrame(
        ...     data={"x": [10, 8, 250, 66], "y": [1, 2, 3, 4]},
        ...     schema={"x": pl.UInt8, "y": pl.UInt8},
        ... )
        >>> df.with_columns(
        ...     pl.col("x")
        ...     .map_elements(binary_string, return_dtype=pl.String)
        ...     .alias("bin_x"),
        ...     pl.col("y")
        ...     .map_elements(binary_string, return_dtype=pl.String)
        ...     .alias("bin_y"),
        ...     pl.col("x").xor(pl.col("y")).alias("xor_xy"),
        ...     pl.col("x")
        ...     .xor(pl.col("y"))
        ...     .map_elements(binary_string, return_dtype=pl.String)
        ...     .alias("bin_xor_xy"),
        ... )
        shape: (4, 6)
        ┌─────┬─────┬──────────┬──────────┬────────┬────────────┐
        │ x   ┆ y   ┆ bin_x    ┆ bin_y    ┆ xor_xy ┆ bin_xor_xy │
        │ --- ┆ --- ┆ ---      ┆ ---      ┆ ---    ┆ ---        │
        │ u8  ┆ u8  ┆ str      ┆ str      ┆ u8     ┆ str        │
        ╞═════╪═════╪══════════╪══════════╪════════╪════════════╡
        │ 10  ┆ 1   ┆ 00001010 ┆ 00000001 ┆ 11     ┆ 00001011   │
        │ 8   ┆ 2   ┆ 00001000 ┆ 00000010 ┆ 10     ┆ 00001010   │
        │ 250 ┆ 3   ┆ 11111010 ┆ 00000011 ┆ 249    ┆ 11111001   │
        │ 66  ┆ 4   ┆ 01000010 ┆ 00000100 ┆ 70     ┆ 01000110   │
        └─────┴─────┴──────────┴──────────┴────────┴────────────┘
        """
        return self.__xor__(other)

    def is_in(
        self,
        other: Expr | Collection[Any] | Series,
        *,
        nulls_equal: bool = False,
    ) -> Expr:
        """
        Check if elements of this expression are present in the other Series.

        Parameters
        ----------
        other
            Series or sequence of primitive type.
        nulls_equal : bool, default False
            If True, treat null as a distinct value. Null values will not propagate.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"sets": [[1, 2, 3], [1, 2], [9, 10]], "optional_members": [1, 2, 3]}
        ... )
        >>> df.with_columns(contains=pl.col("optional_members").is_in("sets"))
        shape: (3, 3)
        ┌───────────┬──────────────────┬──────────┐
        │ sets      ┆ optional_members ┆ contains │
        │ ---       ┆ ---              ┆ ---      │
        │ list[i64] ┆ i64              ┆ bool     │
        ╞═══════════╪══════════════════╪══════════╡
        │ [1, 2, 3] ┆ 1                ┆ true     │
        │ [1, 2]    ┆ 2                ┆ true     │
        │ [9, 10]   ┆ 3                ┆ false    │
        └───────────┴──────────────────┴──────────┘
        """
        if isinstance(other, Collection) and not isinstance(other, (str, pl.Series)):
            other = list(other)  # eg: set, frozenset, etc

        other_pyexpr = parse_into_expression(other)
        return wrap_expr(self._pyexpr.is_in(other_pyexpr, nulls_equal))

    def repeat_by(self, by: pl.Series | Expr | str | int) -> Expr:
        """
        Repeat the elements in this Series as specified in the given expression.

        The repeated elements are expanded into a `List`.

        Parameters
        ----------
        by
            Numeric column that determines how often the values will be repeated.
            The column will be coerced to UInt32. Give this dtype to make the coercion a
            no-op.

        Returns
        -------
        Expr
            Expression of data type :class:`List`, where the inner data type is equal
            to the original data type.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["x", "y", "z"],
        ...         "n": [1, 2, 3],
        ...     }
        ... )
        >>> df.select(pl.col("a").repeat_by("n"))
        shape: (3, 1)
        ┌─────────────────┐
        │ a               │
        │ ---             │
        │ list[str]       │
        ╞═════════════════╡
        │ ["x"]           │
        │ ["y", "y"]      │
        │ ["z", "z", "z"] │
        └─────────────────┘
        """
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(self._pyexpr.repeat_by(by_pyexpr))

    def is_between(
        self,
        lower_bound: IntoExpr,
        upper_bound: IntoExpr,
        closed: ClosedInterval = "both",
    ) -> Expr:
        """
        Check if this expression is between the given lower and upper bounds.

        Parameters
        ----------
        lower_bound
            Lower bound value. Accepts expression input. Strings are parsed as column
            names, other non-expression inputs are parsed as literals.
        upper_bound
            Upper bound value. Accepts expression input. Strings are parsed as column
            names, other non-expression inputs are parsed as literals.
        closed : {'both', 'left', 'right', 'none'}
            Define which sides of the interval are closed (inclusive).

        Notes
        -----
        If the value of the `lower_bound` is greater than that of the `upper_bound`
        then the result will be False, as no value can satisfy the condition.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame({"num": [1, 2, 3, 4, 5]})
        >>> df.with_columns(pl.col("num").is_between(2, 4).alias("is_between"))
        shape: (5, 2)
        ┌─────┬────────────┐
        │ num ┆ is_between │
        │ --- ┆ ---        │
        │ i64 ┆ bool       │
        ╞═════╪════════════╡
        │ 1   ┆ false      │
        │ 2   ┆ true       │
        │ 3   ┆ true       │
        │ 4   ┆ true       │
        │ 5   ┆ false      │
        └─────┴────────────┘

        Use the `closed` argument to include or exclude the values at the bounds:

        >>> df.with_columns(
        ...     pl.col("num").is_between(2, 4, closed="left").alias("is_between")
        ... )
        shape: (5, 2)
        ┌─────┬────────────┐
        │ num ┆ is_between │
        │ --- ┆ ---        │
        │ i64 ┆ bool       │
        ╞═════╪════════════╡
        │ 1   ┆ false      │
        │ 2   ┆ true       │
        │ 3   ┆ true       │
        │ 4   ┆ false      │
        │ 5   ┆ false      │
        └─────┴────────────┘

        You can also use strings as well as numeric/temporal values (note: ensure that
        string literals are wrapped with `lit` so as not to conflate them with
        column names):

        >>> df = pl.DataFrame({"a": ["a", "b", "c", "d", "e"]})
        >>> df.with_columns(
        ...     pl.col("a")
        ...     .is_between(pl.lit("a"), pl.lit("c"), closed="both")
        ...     .alias("is_between")
        ... )
        shape: (5, 2)
        ┌─────┬────────────┐
        │ a   ┆ is_between │
        │ --- ┆ ---        │
        │ str ┆ bool       │
        ╞═════╪════════════╡
        │ a   ┆ true       │
        │ b   ┆ true       │
        │ c   ┆ true       │
        │ d   ┆ false      │
        │ e   ┆ false      │
        └─────┴────────────┘

        Use column expressions as lower/upper bounds, comparing to a literal value:

        >>> df = pl.DataFrame({"a": [1, 2, 3, 4, 5], "b": [5, 4, 3, 2, 1]})
        >>> df.with_columns(
        ...     pl.lit(3).is_between(pl.col("a"), pl.col("b")).alias("between_ab")
        ... )
        shape: (5, 3)
        ┌─────┬─────┬────────────┐
        │ a   ┆ b   ┆ between_ab │
        │ --- ┆ --- ┆ ---        │
        │ i64 ┆ i64 ┆ bool       │
        ╞═════╪═════╪════════════╡
        │ 1   ┆ 5   ┆ true       │
        │ 2   ┆ 4   ┆ true       │
        │ 3   ┆ 3   ┆ true       │
        │ 4   ┆ 2   ┆ false      │
        │ 5   ┆ 1   ┆ false      │
        └─────┴─────┴────────────┘
        """
        lower_bound_pyexpr = parse_into_expression(lower_bound)
        upper_bound_pyexpr = parse_into_expression(upper_bound)

        return wrap_expr(
            self._pyexpr.is_between(lower_bound_pyexpr, upper_bound_pyexpr, closed)
        )

    def is_close(
        self,
        other: IntoExpr,
        *,
        abs_tol: float = 0.0,
        rel_tol: float = 1e-09,
        nans_equal: bool = False,
    ) -> Expr:
        r"""
        Check if this expression is close, i.e. almost equal, to the other expression.

        Two values `a` and `b` are considered close if the following condition holds:

        .. math::
            |a-b| \le max \{ \text{rel_tol} \cdot max \{ |a|, |b| \}, \text{abs_tol} \}

        Parameters
        ----------
        other
            A literal or expression value to compare with.
        abs_tol
            Absolute tolerance. This is the maximum allowed absolute difference between
            two values. Must be non-negative.
        rel_tol
            Relative tolerance. This is the maximum allowed difference between two
            values, relative to the larger absolute value. Must be non-negative.
        nans_equal
            Whether NaN values should be considered equal.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Notes
        -----
            The implementation of this method is symmetric and mirrors the behavior of
            :meth:`math.isclose`. Specifically note that this behavior is different to
            :meth:`numpy.isclose`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.5, 2.0, 2.5], "b": [1.55, 2.2, 3.0]})
        >>> df.with_columns(pl.col("a").is_close("b", abs_tol=0.1).alias("is_close"))
        shape: (3, 3)
        ┌─────┬──────┬──────────┐
        │ a   ┆ b    ┆ is_close │
        │ --- ┆ ---  ┆ ---      │
        │ f64 ┆ f64  ┆ bool     │
        ╞═════╪══════╪══════════╡
        │ 1.5 ┆ 1.55 ┆ true     │
        │ 2.0 ┆ 2.2  ┆ false    │
        │ 2.5 ┆ 3.0  ┆ false    │
        └─────┴──────┴──────────┘
        """
        other_pyexpr = parse_into_expression(other)
        return wrap_expr(
            self._pyexpr.is_close(other_pyexpr, abs_tol, rel_tol, nans_equal)
        )

    def hash(
        self,
        seed: int = 0,
        seed_1: int | None = None,
        seed_2: int | None = None,
        seed_3: int | None = None,
    ) -> Expr:
        """
        Hash the elements in the selection.

        The hash value is of type `UInt64`.

        Parameters
        ----------
        seed
            Random seed parameter. Defaults to 0.
        seed_1
            Random seed parameter. Defaults to `seed` if not set.
        seed_2
            Random seed parameter. Defaults to `seed` if not set.
        seed_3
            Random seed parameter. Defaults to `seed` if not set.

        Notes
        -----
        This implementation of `hash` does not guarantee stable results
        across different Polars versions. Its stability is only guaranteed within a
        single version.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None],
        ...         "b": ["x", None, "z"],
        ...     }
        ... )
        >>> df.with_columns(pl.all().hash(10, 20, 30, 40))  # doctest: +IGNORE_RESULT
        shape: (3, 2)
        ┌──────────────────────┬──────────────────────┐
        │ a                    ┆ b                    │
        │ ---                  ┆ ---                  │
        │ u64                  ┆ u64                  │
        ╞══════════════════════╪══════════════════════╡
        │ 9774092659964970114  ┆ 13614470193936745724 │
        │ 1101441246220388612  ┆ 11638928888656214026 │
        │ 11638928888656214026 ┆ 13382926553367784577 │
        └──────────────────────┴──────────────────────┘
        """
        k0 = seed
        k1 = seed_1 if seed_1 is not None else seed
        k2 = seed_2 if seed_2 is not None else seed
        k3 = seed_3 if seed_3 is not None else seed
        return wrap_expr(self._pyexpr.hash(k0, k1, k2, k3))

    def reinterpret(self, *, signed: bool = True) -> Expr:
        """
        Reinterpret the underlying bits as a signed/unsigned integer.

        This operation is only allowed for 64bit integers. For lower bits integers,
        you can safely use that cast operation.

        Parameters
        ----------
        signed
            If True, reinterpret as `pl.Int64`. Otherwise, reinterpret as `pl.UInt64`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 1, 2], dtype=pl.UInt64)
        >>> df = pl.DataFrame([s])
        >>> df.select(
        ...     [
        ...         pl.col("a").reinterpret(signed=True).alias("reinterpreted"),
        ...         pl.col("a").alias("original"),
        ...     ]
        ... )
        shape: (3, 2)
        ┌───────────────┬──────────┐
        │ reinterpreted ┆ original │
        │ ---           ┆ ---      │
        │ i64           ┆ u64      │
        ╞═══════════════╪══════════╡
        │ 1             ┆ 1        │
        │ 1             ┆ 1        │
        │ 2             ┆ 2        │
        └───────────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.reinterpret(signed))

    def inspect(self, fmt: str = "{}") -> Expr:
        """
        Print the value that this expression evaluates to and pass on the value.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 1, 2]})
        >>> df.select(pl.col("foo").cum_sum().inspect("value is: {}").alias("bar"))
        value is: shape: (3,)
        Series: 'foo' [i64]
        [
            1
            2
            4
        ]
        shape: (3, 1)
        ┌─────┐
        │ bar │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 4   │
        └─────┘
        """

        def inspect(s: Series) -> Series:  # pragma: no cover
            print(fmt.format(s))
            return s

        return self.map_batches(inspect, return_dtype=F.dtype_of(self))

    def interpolate(self, method: InterpolationMethod = "linear") -> Expr:
        """
        Interpolate intermediate values.

        Nulls at the beginning and end of the series remain null.

        Parameters
        ----------
        method : {'linear', 'nearest'}
            Interpolation method.

        Examples
        --------
        Fill null values using linear interpolation.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, None, 3],
        ...         "b": [1.0, float("nan"), 3.0],
        ...     }
        ... )
        >>> df.select(pl.all().interpolate())
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ f64 ┆ f64 │
        ╞═════╪═════╡
        │ 1.0 ┆ 1.0 │
        │ 2.0 ┆ NaN │
        │ 3.0 ┆ 3.0 │
        └─────┴─────┘

        Fill null values using nearest interpolation.

        >>> df.select(pl.all().interpolate("nearest"))
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ f64 │
        ╞═════╪═════╡
        │ 1   ┆ 1.0 │
        │ 3   ┆ NaN │
        │ 3   ┆ 3.0 │
        └─────┴─────┘

        Regrid data to a new grid.

        >>> df_original_grid = pl.DataFrame(
        ...     {
        ...         "grid_points": [1, 3, 10],
        ...         "values": [2.0, 6.0, 20.0],
        ...     }
        ... )  # Interpolate from this to the new grid
        >>> df_new_grid = pl.DataFrame({"grid_points": range(1, 11)})
        >>> df_new_grid.join(
        ...     df_original_grid, on="grid_points", how="left", coalesce=True
        ... ).with_columns(pl.col("values").interpolate())
        shape: (10, 2)
        ┌─────────────┬────────┐
        │ grid_points ┆ values │
        │ ---         ┆ ---    │
        │ i64         ┆ f64    │
        ╞═════════════╪════════╡
        │ 1           ┆ 2.0    │
        │ 2           ┆ 4.0    │
        │ 3           ┆ 6.0    │
        │ 4           ┆ 8.0    │
        │ 5           ┆ 10.0   │
        │ 6           ┆ 12.0   │
        │ 7           ┆ 14.0   │
        │ 8           ┆ 16.0   │
        │ 9           ┆ 18.0   │
        │ 10          ┆ 20.0   │
        └─────────────┴────────┘
        """
        return wrap_expr(self._pyexpr.interpolate(method))

    def interpolate_by(self, by: IntoExpr) -> Expr:
        """
        Fill null values using interpolation based on another column.

        Nulls at the beginning and end of the series remain null.

        Parameters
        ----------
        by
            Column to interpolate values based on.

        Examples
        --------
        Fill null values using linear interpolation.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, None, None, 3],
        ...         "b": [1, 2, 7, 8],
        ...     }
        ... )
        >>> df.with_columns(a_interpolated=pl.col("a").interpolate_by("b"))
        shape: (4, 3)
        ┌──────┬─────┬────────────────┐
        │ a    ┆ b   ┆ a_interpolated │
        │ ---  ┆ --- ┆ ---            │
        │ i64  ┆ i64 ┆ f64            │
        ╞══════╪═════╪════════════════╡
        │ 1    ┆ 1   ┆ 1.0            │
        │ null ┆ 2   ┆ 1.285714       │
        │ null ┆ 7   ┆ 2.714286       │
        │ 3    ┆ 8   ┆ 3.0            │
        └──────┴─────┴────────────────┘
        """
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(self._pyexpr.interpolate_by(by_pyexpr))

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_min_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Apply a rolling min based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling min with the temporal windows closed on the right (default)

        >>> df_temporal.with_columns(
        ...     rolling_row_min=pl.col("index").rolling_min_by("date", window_size="2h")
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_min │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ u32             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0               │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0               │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1               │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 2               │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 3               │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 19              │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 20              │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 21              │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 22              │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 23              │
        └───────┴─────────────────────┴─────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_min_by(by_pyexpr, window_size, min_samples, closed)
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_max_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Apply a rolling max based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling max with the temporal windows closed on the right (default)

        >>> df_temporal.with_columns(
        ...     rolling_row_max=pl.col("index").rolling_max_by("date", window_size="2h")
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_max │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ u32             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0               │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 1               │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 2               │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 3               │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 4               │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 20              │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 21              │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 22              │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 23              │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 24              │
        └───────┴─────────────────────┴─────────────────┘

        Compute the rolling max with the closure of windows on both sides

        >>> df_temporal.with_columns(
        ...     rolling_row_max=pl.col("index").rolling_max_by(
        ...         "date", window_size="2h", closed="both"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_max │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ u32             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0               │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 1               │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 2               │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 3               │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 4               │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 20              │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 21              │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 22              │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 23              │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 24              │
        └───────┴─────────────────────┴─────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_max_by(by_pyexpr, window_size, min_samples, closed)
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_mean_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Apply a rolling mean based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling mean with the temporal windows closed on the right (default)

        >>> df_temporal.with_columns(
        ...     rolling_row_mean=pl.col("index").rolling_mean_by(
        ...         "date", window_size="2h"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬──────────────────┐
        │ index ┆ date                ┆ rolling_row_mean │
        │ ---   ┆ ---                 ┆ ---              │
        │ u32   ┆ datetime[μs]        ┆ f64              │
        ╞═══════╪═════════════════════╪══════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0.0              │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.5              │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1.5              │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 2.5              │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 3.5              │
        │ …     ┆ …                   ┆ …                │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 19.5             │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 20.5             │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 21.5             │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 22.5             │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 23.5             │
        └───────┴─────────────────────┴──────────────────┘

        Compute the rolling mean with the closure of windows on both sides

        >>> df_temporal.with_columns(
        ...     rolling_row_mean=pl.col("index").rolling_mean_by(
        ...         "date", window_size="2h", closed="both"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬──────────────────┐
        │ index ┆ date                ┆ rolling_row_mean │
        │ ---   ┆ ---                 ┆ ---              │
        │ u32   ┆ datetime[μs]        ┆ f64              │
        ╞═══════╪═════════════════════╪══════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0.0              │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.5              │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1.0              │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 2.0              │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 3.0              │
        │ …     ┆ …                   ┆ …                │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 19.0             │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 20.0             │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 21.0             │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 22.0             │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 23.0             │
        └───────┴─────────────────────┴──────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_mean_by(
                by_pyexpr,
                window_size,
                min_samples,
                closed,
            )
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_sum_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Apply a rolling sum based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling sum with the temporal windows closed on the right (default)

        >>> df_temporal.with_columns(
        ...     rolling_row_sum=pl.col("index").rolling_sum_by("date", window_size="2h")
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_sum │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ u32             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0               │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 1               │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 3               │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 5               │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 7               │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 39              │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 41              │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 43              │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 45              │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 47              │
        └───────┴─────────────────────┴─────────────────┘

        Compute the rolling sum with the closure of windows on both sides

        >>> df_temporal.with_columns(
        ...     rolling_row_sum=pl.col("index").rolling_sum_by(
        ...         "date", window_size="2h", closed="both"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_sum │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ u32             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0               │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 1               │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 3               │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 6               │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 9               │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 57              │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 60              │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 63              │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 66              │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 69              │
        └───────┴─────────────────────┴─────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_sum_by(by_pyexpr, window_size, min_samples, closed)
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_std_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
        ddof: int = 1,
    ) -> Expr:
        """
        Compute a rolling standard deviation based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling std with the temporal windows closed on the right (default)

        >>> df_temporal.with_columns(
        ...     rolling_row_std=pl.col("index").rolling_std_by("date", window_size="2h")
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_std │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ f64             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ null            │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.707107        │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 0.707107        │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 0.707107        │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 0.707107        │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 0.707107        │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 0.707107        │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 0.707107        │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 0.707107        │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 0.707107        │
        └───────┴─────────────────────┴─────────────────┘

        Compute the rolling std with the closure of windows on both sides

        >>> df_temporal.with_columns(
        ...     rolling_row_std=pl.col("index").rolling_std_by(
        ...         "date", window_size="2h", closed="both"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_std │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ f64             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ null            │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.707107        │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1.0             │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 1.0             │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 1.0             │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 1.0             │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 1.0             │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 1.0             │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 1.0             │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 1.0             │
        └───────┴─────────────────────┴─────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_std_by(
                by_pyexpr,
                window_size,
                min_samples,
                closed,
                ddof,
            )
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_var_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
        ddof: int = 1,
    ) -> Expr:
        """
        Compute a rolling variance based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling var with the temporal windows closed on the right (default)

        >>> df_temporal.with_columns(
        ...     rolling_row_var=pl.col("index").rolling_var_by("date", window_size="2h")
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_var │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ f64             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ null            │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.5             │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 0.5             │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 0.5             │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 0.5             │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 0.5             │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 0.5             │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 0.5             │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 0.5             │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 0.5             │
        └───────┴─────────────────────┴─────────────────┘

        Compute the rolling var with the closure of windows on both sides

        >>> df_temporal.with_columns(
        ...     rolling_row_var=pl.col("index").rolling_var_by(
        ...         "date", window_size="2h", closed="both"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬─────────────────┐
        │ index ┆ date                ┆ rolling_row_var │
        │ ---   ┆ ---                 ┆ ---             │
        │ u32   ┆ datetime[μs]        ┆ f64             │
        ╞═══════╪═════════════════════╪═════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ null            │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.5             │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1.0             │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 1.0             │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 1.0             │
        │ …     ┆ …                   ┆ …               │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 1.0             │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 1.0             │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 1.0             │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 1.0             │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 1.0             │
        └───────┴─────────────────────┴─────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_var_by(
                by_pyexpr,
                window_size,
                min_samples,
                closed,
                ddof,
            )
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_median_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Compute a rolling median based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling median with the temporal windows closed on the right:

        >>> df_temporal.with_columns(
        ...     rolling_row_median=pl.col("index").rolling_median_by(
        ...         "date", window_size="2h"
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬────────────────────┐
        │ index ┆ date                ┆ rolling_row_median │
        │ ---   ┆ ---                 ┆ ---                │
        │ u32   ┆ datetime[μs]        ┆ f64                │
        ╞═══════╪═════════════════════╪════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0.0                │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.5                │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1.5                │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 2.5                │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 3.5                │
        │ …     ┆ …                   ┆ …                  │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 19.5               │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 20.5               │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 21.5               │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 22.5               │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 23.5               │
        └───────┴─────────────────────┴────────────────────┘
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_median_by(by_pyexpr, window_size, min_samples, closed)
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_quantile_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        quantile: float,
        interpolation: QuantileMethod = "nearest",
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Compute a rolling quantile based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.
        window_size
            The length of the window. Can be a dynamic
            temporal size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a DataFrame with a datetime column and a row number column

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> df_temporal = pl.DataFrame(
        ...     {"date": pl.datetime_range(start, stop, "1h", eager=True)}
        ... ).with_row_index()
        >>> df_temporal
        shape: (25, 2)
        ┌───────┬─────────────────────┐
        │ index ┆ date                │
        │ ---   ┆ ---                 │
        │ u32   ┆ datetime[μs]        │
        ╞═══════╪═════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 │
        │ 1     ┆ 2001-01-01 01:00:00 │
        │ 2     ┆ 2001-01-01 02:00:00 │
        │ 3     ┆ 2001-01-01 03:00:00 │
        │ 4     ┆ 2001-01-01 04:00:00 │
        │ …     ┆ …                   │
        │ 20    ┆ 2001-01-01 20:00:00 │
        │ 21    ┆ 2001-01-01 21:00:00 │
        │ 22    ┆ 2001-01-01 22:00:00 │
        │ 23    ┆ 2001-01-01 23:00:00 │
        │ 24    ┆ 2001-01-02 00:00:00 │
        └───────┴─────────────────────┘

        Compute the rolling quantile with the temporal windows closed on the right:

        >>> df_temporal.with_columns(
        ...     rolling_row_quantile=pl.col("index").rolling_quantile_by(
        ...         "date", window_size="2h", quantile=0.3
        ...     )
        ... )
        shape: (25, 3)
        ┌───────┬─────────────────────┬──────────────────────┐
        │ index ┆ date                ┆ rolling_row_quantile │
        │ ---   ┆ ---                 ┆ ---                  │
        │ u32   ┆ datetime[μs]        ┆ f64                  │
        ╞═══════╪═════════════════════╪══════════════════════╡
        │ 0     ┆ 2001-01-01 00:00:00 ┆ 0.0                  │
        │ 1     ┆ 2001-01-01 01:00:00 ┆ 0.0                  │
        │ 2     ┆ 2001-01-01 02:00:00 ┆ 1.0                  │
        │ 3     ┆ 2001-01-01 03:00:00 ┆ 2.0                  │
        │ 4     ┆ 2001-01-01 04:00:00 ┆ 3.0                  │
        │ …     ┆ …                   ┆ …                    │
        │ 20    ┆ 2001-01-01 20:00:00 ┆ 19.0                 │
        │ 21    ┆ 2001-01-01 21:00:00 ┆ 20.0                 │
        │ 22    ┆ 2001-01-01 22:00:00 ┆ 21.0                 │
        │ 23    ┆ 2001-01-01 23:00:00 ┆ 22.0                 │
        │ 24    ┆ 2001-01-02 00:00:00 ┆ 23.0                 │
        └───────┴─────────────────────┴──────────────────────┘
        """  # noqa: W505
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_quantile_by(
                by_pyexpr,
                quantile,
                interpolation,
                window_size,
                min_samples,
                closed,
            )
        )

    @unstable()
    def rolling_rank_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        method: RankMethod = "average",
        *,
        seed: int | None = None,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Expr:
        """
        Compute a rolling rank based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic
            temporal size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        method : {'average', 'min', 'max', 'dense', 'random'}
            The method used to assign ranks to tied elements.
            The following methods are available (default is 'average'):

            - 'average' : The average of the ranks that would have been assigned to
              all the tied values is assigned to each value.
            - 'min' : The minimum of the ranks that would have been assigned to all
              the tied values is assigned to each value. (This is also referred to
              as "competition" ranking.)
            - 'max' : The maximum of the ranks that would have been assigned to all
              the tied values is assigned to each value.
            - 'dense' : Like 'min', but the rank of the next highest element is
              assigned the rank immediately after those assigned to the tied
              elements.
            - 'random' : Choose a random rank for each value in a tie.
        seed
            Random seed used when `method='random'`. If set to None (default), a
            random seed is generated for each rolling rank operation.
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Returns
        -------
        Expr
            An Expr of data :class:`.Float64` if `method` is `"average"` or,
            the index size (see :func:`.get_index_type()`) otherwise.
        """
        window_size = _prepare_rolling_by_window_args(window_size)
        by_pyexpr = parse_into_expression(by)
        return wrap_expr(
            self._pyexpr.rolling_rank_by(
                by_pyexpr,
                window_size,
                method,
                seed,
                min_samples,
                closed,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_min(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Apply a rolling min (moving min) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their min.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_min=pl.col("A").rolling_min(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_min │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 1.0         │
        │ 3.0 ┆ 2.0         │
        │ 4.0 ┆ 3.0         │
        │ 5.0 ┆ 4.0         │
        │ 6.0 ┆ 5.0         │
        └─────┴─────────────┘

        Specify weights to multiply the values in the window with:

        >>> df.with_columns(
        ...     rolling_min=pl.col("A").rolling_min(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_min │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 0.25        │
        │ 3.0 ┆ 0.5         │
        │ 4.0 ┆ 0.75        │
        │ 5.0 ┆ 1.0         │
        │ 6.0 ┆ 1.25        │
        └─────┴─────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_min=pl.col("A").rolling_min(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_min │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 1.0         │
        │ 3.0 ┆ 2.0         │
        │ 4.0 ┆ 3.0         │
        │ 5.0 ┆ 4.0         │
        │ 6.0 ┆ null        │
        └─────┴─────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_min(
                window_size,
                weights,
                min_samples,
                center=center,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_max(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Apply a rolling max (moving max) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their max.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_max=pl.col("A").rolling_max(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_max │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 2.0         │
        │ 3.0 ┆ 3.0         │
        │ 4.0 ┆ 4.0         │
        │ 5.0 ┆ 5.0         │
        │ 6.0 ┆ 6.0         │
        └─────┴─────────────┘

        Specify weights to multiply the values in the window with:

        >>> df.with_columns(
        ...     rolling_max=pl.col("A").rolling_max(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_max │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 1.5         │
        │ 3.0 ┆ 2.25        │
        │ 4.0 ┆ 3.0         │
        │ 5.0 ┆ 3.75        │
        │ 6.0 ┆ 4.5         │
        └─────┴─────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_max=pl.col("A").rolling_max(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_max │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 3.0         │
        │ 3.0 ┆ 4.0         │
        │ 4.0 ┆ 5.0         │
        │ 5.0 ┆ 6.0         │
        │ 6.0 ┆ null        │
        └─────┴─────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_max(
                window_size,
                weights,
                min_samples,
                center,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_mean(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Apply a rolling mean (moving mean) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their mean. Weights
        are normalized to sum to 1.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window, after being normalized to sum to
            1.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_mean=pl.col("A").rolling_mean(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────┐
        │ A   ┆ rolling_mean │
        │ --- ┆ ---          │
        │ f64 ┆ f64          │
        ╞═════╪══════════════╡
        │ 1.0 ┆ null         │
        │ 2.0 ┆ 1.5          │
        │ 3.0 ┆ 2.5          │
        │ 4.0 ┆ 3.5          │
        │ 5.0 ┆ 4.5          │
        │ 6.0 ┆ 5.5          │
        └─────┴──────────────┘

        Specify weights to multiply the values in the window with:

        >>> df.with_columns(
        ...     rolling_mean=pl.col("A").rolling_mean(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────┐
        │ A   ┆ rolling_mean │
        │ --- ┆ ---          │
        │ f64 ┆ f64          │
        ╞═════╪══════════════╡
        │ 1.0 ┆ null         │
        │ 2.0 ┆ 1.75         │
        │ 3.0 ┆ 2.75         │
        │ 4.0 ┆ 3.75         │
        │ 5.0 ┆ 4.75         │
        │ 6.0 ┆ 5.75         │
        └─────┴──────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_mean=pl.col("A").rolling_mean(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────┐
        │ A   ┆ rolling_mean │
        │ --- ┆ ---          │
        │ f64 ┆ f64          │
        ╞═════╪══════════════╡
        │ 1.0 ┆ null         │
        │ 2.0 ┆ 2.0          │
        │ 3.0 ┆ 3.0          │
        │ 4.0 ┆ 4.0          │
        │ 5.0 ┆ 5.0          │
        │ 6.0 ┆ null         │
        └─────┴──────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_mean(
                window_size,
                weights,
                min_samples,
                center,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_sum(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Apply a rolling sum (moving sum) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their sum.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_sum=pl.col("A").rolling_sum(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_sum │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 3.0         │
        │ 3.0 ┆ 5.0         │
        │ 4.0 ┆ 7.0         │
        │ 5.0 ┆ 9.0         │
        │ 6.0 ┆ 11.0        │
        └─────┴─────────────┘

        Specify weights to multiply the values in the window with:

        >>> df.with_columns(
        ...     rolling_sum=pl.col("A").rolling_sum(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_sum │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 1.75        │
        │ 3.0 ┆ 2.75        │
        │ 4.0 ┆ 3.75        │
        │ 5.0 ┆ 4.75        │
        │ 6.0 ┆ 5.75        │
        └─────┴─────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_sum=pl.col("A").rolling_sum(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_sum │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 6.0         │
        │ 3.0 ┆ 9.0         │
        │ 4.0 ┆ 12.0        │
        │ 5.0 ┆ 15.0        │
        │ 6.0 ┆ null        │
        └─────┴─────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_sum(
                window_size,
                weights,
                min_samples,
                center,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_std(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Expr:
        """
        Compute a rolling standard deviation.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their std. Weights
        are normalized to sum to 1.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window after being normalized to sum to
            1.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_std=pl.col("A").rolling_std(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_std │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 0.707107    │
        │ 3.0 ┆ 0.707107    │
        │ 4.0 ┆ 0.707107    │
        │ 5.0 ┆ 0.707107    │
        │ 6.0 ┆ 0.707107    │
        └─────┴─────────────┘

        Specify weights to multiply the values in the window with:

        >>> df.with_columns(
        ...     rolling_std=pl.col("A").rolling_std(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_std │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 0.433013    │
        │ 3.0 ┆ 0.433013    │
        │ 4.0 ┆ 0.433013    │
        │ 5.0 ┆ 0.433013    │
        │ 6.0 ┆ 0.433013    │
        └─────┴─────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_std=pl.col("A").rolling_std(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_std │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 1.0         │
        │ 3.0 ┆ 1.0         │
        │ 4.0 ┆ 1.0         │
        │ 5.0 ┆ 1.0         │
        │ 6.0 ┆ null        │
        └─────┴─────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_std(
                window_size,
                weights,
                min_samples,
                center=center,
                ddof=ddof,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_var(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Expr:
        """
        Compute a rolling variance.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their var. Weights
        are normalized to sum to 1.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window after being normalized to sum to
            1.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_var=pl.col("A").rolling_var(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_var │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 0.5         │
        │ 3.0 ┆ 0.5         │
        │ 4.0 ┆ 0.5         │
        │ 5.0 ┆ 0.5         │
        │ 6.0 ┆ 0.5         │
        └─────┴─────────────┘

        Specify weights to multiply the values in the window with:

        >>> df.with_columns(
        ...     rolling_var=pl.col("A").rolling_var(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_var │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 0.1875      │
        │ 3.0 ┆ 0.1875      │
        │ 4.0 ┆ 0.1875      │
        │ 5.0 ┆ 0.1875      │
        │ 6.0 ┆ 0.1875      │
        └─────┴─────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_var=pl.col("A").rolling_var(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬─────────────┐
        │ A   ┆ rolling_var │
        │ --- ┆ ---         │
        │ f64 ┆ f64         │
        ╞═════╪═════════════╡
        │ 1.0 ┆ null        │
        │ 2.0 ┆ 1.0         │
        │ 3.0 ┆ 1.0         │
        │ 4.0 ┆ 1.0         │
        │ 5.0 ┆ 1.0         │
        │ 6.0 ┆ null        │
        └─────┴─────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_var(
                window_size,
                weights,
                min_samples,
                center=center,
                ddof=ddof,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_median(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Compute a rolling median.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their median.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_median=pl.col("A").rolling_median(window_size=2),
        ... )
        shape: (6, 2)
        ┌─────┬────────────────┐
        │ A   ┆ rolling_median │
        │ --- ┆ ---            │
        │ f64 ┆ f64            │
        ╞═════╪════════════════╡
        │ 1.0 ┆ null           │
        │ 2.0 ┆ 1.5            │
        │ 3.0 ┆ 2.5            │
        │ 4.0 ┆ 3.5            │
        │ 5.0 ┆ 4.5            │
        │ 6.0 ┆ 5.5            │
        └─────┴────────────────┘

        Specify weights for the values in each window:

        >>> df.with_columns(
        ...     rolling_median=pl.col("A").rolling_median(
        ...         window_size=2, weights=[0.25, 0.75]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬────────────────┐
        │ A   ┆ rolling_median │
        │ --- ┆ ---            │
        │ f64 ┆ f64            │
        ╞═════╪════════════════╡
        │ 1.0 ┆ null           │
        │ 2.0 ┆ 1.5            │
        │ 3.0 ┆ 2.5            │
        │ 4.0 ┆ 3.5            │
        │ 5.0 ┆ 4.5            │
        │ 6.0 ┆ 5.5            │
        └─────┴────────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_median=pl.col("A").rolling_median(window_size=3, center=True),
        ... )
        shape: (6, 2)
        ┌─────┬────────────────┐
        │ A   ┆ rolling_median │
        │ --- ┆ ---            │
        │ f64 ┆ f64            │
        ╞═════╪════════════════╡
        │ 1.0 ┆ null           │
        │ 2.0 ┆ 2.0            │
        │ 3.0 ┆ 3.0            │
        │ 4.0 ┆ 4.0            │
        │ 5.0 ┆ 5.0            │
        │ 6.0 ┆ null           │
        └─────┴────────────────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_median(
                window_size,
                weights,
                min_samples,
                center=center,
            )
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_quantile(
        self,
        quantile: float,
        interpolation: QuantileMethod = "nearest",
        window_size: int = 2,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Compute a rolling quantile.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weights` vector. The resulting values will be aggregated to their quantile.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        >>> df = pl.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
        >>> df.with_columns(
        ...     rolling_quantile=pl.col("A").rolling_quantile(
        ...         quantile=0.25, window_size=4
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────────┐
        │ A   ┆ rolling_quantile │
        │ --- ┆ ---              │
        │ f64 ┆ f64              │
        ╞═════╪══════════════════╡
        │ 1.0 ┆ null             │
        │ 2.0 ┆ null             │
        │ 3.0 ┆ null             │
        │ 4.0 ┆ 2.0              │
        │ 5.0 ┆ 3.0              │
        │ 6.0 ┆ 4.0              │
        └─────┴──────────────────┘

        Specify weights for the values in each window:

        >>> df.with_columns(
        ...     rolling_quantile=pl.col("A").rolling_quantile(
        ...         quantile=0.25, window_size=4, weights=[0.2, 0.4, 0.4, 0.2]
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────────┐
        │ A   ┆ rolling_quantile │
        │ --- ┆ ---              │
        │ f64 ┆ f64              │
        ╞═════╪══════════════════╡
        │ 1.0 ┆ null             │
        │ 2.0 ┆ null             │
        │ 3.0 ┆ null             │
        │ 4.0 ┆ 2.0              │
        │ 5.0 ┆ 3.0              │
        │ 6.0 ┆ 4.0              │
        └─────┴──────────────────┘

        Specify weights and interpolation method

        >>> df.with_columns(
        ...     rolling_quantile=pl.col("A").rolling_quantile(
        ...         quantile=0.25,
        ...         window_size=4,
        ...         weights=[0.2, 0.4, 0.4, 0.2],
        ...         interpolation="linear",
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────────┐
        │ A   ┆ rolling_quantile │
        │ --- ┆ ---              │
        │ f64 ┆ f64              │
        ╞═════╪══════════════════╡
        │ 1.0 ┆ null             │
        │ 2.0 ┆ null             │
        │ 3.0 ┆ null             │
        │ 4.0 ┆ 1.625            │
        │ 5.0 ┆ 2.625            │
        │ 6.0 ┆ 3.625            │
        └─────┴──────────────────┘

        Center the values in the window

        >>> df.with_columns(
        ...     rolling_quantile=pl.col("A").rolling_quantile(
        ...         quantile=0.2, window_size=5, center=True
        ...     ),
        ... )
        shape: (6, 2)
        ┌─────┬──────────────────┐
        │ A   ┆ rolling_quantile │
        │ --- ┆ ---              │
        │ f64 ┆ f64              │
        ╞═════╪══════════════════╡
        │ 1.0 ┆ null             │
        │ 2.0 ┆ null             │
        │ 3.0 ┆ 2.0              │
        │ 4.0 ┆ 3.0              │
        │ 5.0 ┆ null             │
        │ 6.0 ┆ null             │
        └─────┴──────────────────┘
        """  # noqa: W505
        return wrap_expr(
            self._pyexpr.rolling_quantile(
                quantile,
                interpolation,
                window_size,
                weights,
                min_samples,
                center=center,
            )
        )

    @unstable()
    def rolling_rank(
        self,
        window_size: int,
        method: RankMethod = "average",
        *,
        seed: int | None = None,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Compute a rolling rank.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        A window of length `window_size` will traverse the array. The values
        that fill this window will be ranked according to the `method`
        parameter. The resulting values will be the rank of the value that is
        at the end of the sliding window.

        Parameters
        ----------
        window_size
            Integer size of the rolling window.
        method : {'average', 'min', 'max', 'dense', 'random'}
            The method used to assign ranks to tied elements.
            The following methods are available (default is 'average'):

            - 'average' : The average of the ranks that would have been assigned to
              all the tied values is assigned to each value.
            - 'min' : The minimum of the ranks that would have been assigned to all
              the tied values is assigned to each value. (This is also referred to
              as "competition" ranking.)
            - 'max' : The maximum of the ranks that would have been assigned to all
              the tied values is assigned to each value.
            - 'dense' : Like 'min', but the rank of the next highest element is
              assigned the rank immediately after those assigned to the tied
              elements.
            - 'random' : Choose a random rank for each value in a tie.
        seed
            Random seed used when `method='random'`. If set to None (default), a
            random seed is generated for each rolling rank operation.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Returns
        -------
        Expr
            An Expr of data :class:`.Float64` if `method` is `"average"` or,
            the index size (see :func:`.get_index_type()`) otherwise.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 4, 4, 1, 9]})
        >>> df.select(pl.col("a").rolling_rank(3, method="average"))
            shape: (5, 1)
            ┌──────┐
            │ a    │
            │ ---  │
            │ f64  │
            ╞══════╡
            │ null │
            │ null │
            │ 2.5  │
            │ 1.0  │
            │ 3.0  │
            └──────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_rank(
                window_size,
                method,
                seed,
                min_samples,
                center,
            )
        )

    @unstable()
    def rolling_skew(
        self,
        window_size: int,
        *,
        bias: bool = True,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Compute a rolling skew.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        Parameters
        ----------
        window_size
            Integer size of the rolling window.
        bias
            If False, the calculations are corrected for statistical bias.
                     bias: bool = True,
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        See Also
        --------
        Expr.skew

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 4, 2, 9]})
        >>> df.select(pl.col("a").rolling_skew(3))
        shape: (4, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ null     │
        │ null     │
        │ 0.381802 │
        │ 0.47033  │
        └──────────┘

        Note how the values match the following:

        >>> pl.Series([1, 4, 2]).skew(), pl.Series([4, 2, 9]).skew()
        (0.38180177416060584, 0.47033046033698594)
        """
        return wrap_expr(
            self._pyexpr.rolling_skew(
                window_size, bias=bias, min_periods=min_samples, center=center
            )
        )

    @unstable()
    def rolling_kurtosis(
        self,
        window_size: int,
        *,
        fisher: bool = True,
        bias: bool = True,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Compute a rolling kurtosis.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        Parameters
        ----------
        window_size
            Integer size of the rolling window.
        fisher : bool, optional
            If True, Fisher's definition is used (normal ==> 0.0). If False,
            Pearson's definition is used (normal ==> 3.0).
        bias : bool, optional
            If False, the calculations are corrected for statistical bias.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        See Also
        --------
        Expr.kurtosis

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 4, 2, 9]})
        >>> df.select(pl.col("a").rolling_kurtosis(3))
        shape: (4, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ null │
        │ null │
        │ -1.5 │
        │ -1.5 │
        └──────┘
        """
        return wrap_expr(
            self._pyexpr.rolling_kurtosis(
                window_size,
                fisher=fisher,
                bias=bias,
                min_periods=min_samples,
                center=center,
            )
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_map(
        self,
        function: Callable[[Series], Any],
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Expr:
        """
        Compute a custom rolling window function.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        function
            Custom aggregation function.
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Warnings
        --------
        Computing custom functions is extremely slow. Use specialized rolling
        functions such as :func:`Expr.rolling_sum` if at all possible.

        Examples
        --------
        >>> from numpy import nansum
        >>> df = pl.DataFrame({"a": [11.0, 2.0, 9.0, float("nan"), 8.0]})
        >>> df.select(pl.col("a").rolling_map(nansum, window_size=3))
        shape: (5, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ null │
        │ null │
        │ 22.0 │
        │ 11.0 │
        │ 17.0 │
        └──────┘
        """
        if min_samples is None:
            min_samples = window_size

        def _wrap(pys: PySeries) -> PySeries:
            s = wrap_s(pys)
            rv = function(s)
            if isinstance(rv, pl.Series):
                return rv._s
            return pl.Series([rv])._s

        return wrap_expr(
            self._pyexpr.rolling_map(_wrap, window_size, weights, min_samples, center)
        )

    def abs(self) -> Expr:
        """
        Compute absolute values.

        Same as `abs(expr)`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "A": [-1.0, 0.0, 1.0, 2.0],
        ...     }
        ... )
        >>> df.select(pl.col("A").abs())
        shape: (4, 1)
        ┌─────┐
        │ A   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        │ 0.0 │
        │ 1.0 │
        │ 2.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.abs())

    def rank(
        self,
        method: RankMethod = "average",
        *,
        descending: bool = False,
        seed: int | None = None,
    ) -> Expr:
        """
        Assign ranks to data, dealing with ties appropriately.

        Parameters
        ----------
        method : {'average', 'min', 'max', 'dense', 'ordinal', 'random'}
            The method used to assign ranks to tied elements.
            The following methods are available (default is 'average'):

            - 'average' : The average of the ranks that would have been assigned to
              all the tied values is assigned to each value.
            - 'min' : The minimum of the ranks that would have been assigned to all
              the tied values is assigned to each value. (This is also referred to
              as "competition" ranking.)
            - 'max' : The maximum of the ranks that would have been assigned to all
              the tied values is assigned to each value.
            - 'dense' : Like 'min', but the rank of the next highest element is
              assigned the rank immediately after those assigned to the tied
              elements.
            - 'ordinal' : All values are given a distinct rank, corresponding to
              the order that the values occur in the Series.
            - 'random' : Like 'ordinal', but the rank for ties is not dependent
              on the order that the values occur in the Series.
        descending
            Rank in descending order.
        seed
            If `method="random"`, use this as seed.

        Examples
        --------
        The 'average' method:

        >>> df = pl.DataFrame({"a": [3, 6, 1, 1, 6]})
        >>> df.select(pl.col("a").rank())
        shape: (5, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 3.0 │
        │ 4.5 │
        │ 1.5 │
        │ 1.5 │
        │ 4.5 │
        └─────┘

        The 'ordinal' method:

        >>> df = pl.DataFrame({"a": [3, 6, 1, 1, 6]})
        >>> df.select(pl.col("a").rank("ordinal"))
        shape: (5, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 3   │
        │ 4   │
        │ 1   │
        │ 2   │
        │ 5   │
        └─────┘

        Use 'rank' with 'over' to rank within groups:

        >>> df = pl.DataFrame({"a": [1, 1, 2, 2, 2], "b": [6, 7, 5, 14, 11]})
        >>> df.with_columns(pl.col("b").rank().over("a").alias("rank"))
        shape: (5, 3)
        ┌─────┬─────┬──────┐
        │ a   ┆ b   ┆ rank │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ i64 ┆ f64  │
        ╞═════╪═════╪══════╡
        │ 1   ┆ 6   ┆ 1.0  │
        │ 1   ┆ 7   ┆ 2.0  │
        │ 2   ┆ 5   ┆ 1.0  │
        │ 2   ┆ 14  ┆ 3.0  │
        │ 2   ┆ 11  ┆ 2.0  │
        └─────┴─────┴──────┘

        Divide by the length or number of non-null values
        to compute the percentile rank.

        >>> df = pl.DataFrame({"a": [6, 7, None, 14, 11]})
        >>> df.with_columns(
        ...     pct=pl.col("a").rank() / pl.len(),
        ...     pct_valid=pl.col("a").rank() / pl.count("a"),
        ... )
        shape: (5, 3)
        ┌──────┬──────┬───────────┐
        │ a    ┆ pct  ┆ pct_valid │
        │ ---  ┆ ---  ┆ ---       │
        │ i64  ┆ f64  ┆ f64       │
        ╞══════╪══════╪═══════════╡
        │ 6    ┆ 0.2  ┆ 0.25      │
        │ 7    ┆ 0.4  ┆ 0.5       │
        │ null ┆ null ┆ null      │
        │ 14   ┆ 0.8  ┆ 1.0       │
        │ 11   ┆ 0.6  ┆ 0.75      │
        └──────┴──────┴───────────┘

        """
        return wrap_expr(self._pyexpr.rank(method, descending, seed))

    def diff(
        self, n: int | IntoExpr = 1, null_behavior: NullBehavior = "ignore"
    ) -> Expr:
        """
        Calculate the first discrete difference between shifted items.

        Parameters
        ----------
        n
            Number of slots to shift.
        null_behavior : {'ignore', 'drop'}
            How to handle null values.

        Examples
        --------
        >>> df = pl.DataFrame({"int": [20, 10, 30, 25, 35]})
        >>> df.with_columns(change=pl.col("int").diff())
        shape: (5, 2)
        ┌─────┬────────┐
        │ int ┆ change │
        │ --- ┆ ---    │
        │ i64 ┆ i64    │
        ╞═════╪════════╡
        │ 20  ┆ null   │
        │ 10  ┆ -10    │
        │ 30  ┆ 20     │
        │ 25  ┆ -5     │
        │ 35  ┆ 10     │
        └─────┴────────┘

        >>> df.with_columns(change=pl.col("int").diff(n=2))
        shape: (5, 2)
        ┌─────┬────────┐
        │ int ┆ change │
        │ --- ┆ ---    │
        │ i64 ┆ i64    │
        ╞═════╪════════╡
        │ 20  ┆ null   │
        │ 10  ┆ null   │
        │ 30  ┆ 10     │
        │ 25  ┆ 15     │
        │ 35  ┆ 5      │
        └─────┴────────┘

        >>> df.select(pl.col("int").diff(n=2, null_behavior="drop").alias("diff"))
        shape: (3, 1)
        ┌──────┐
        │ diff │
        │ ---  │
        │ i64  │
        ╞══════╡
        │ 10   │
        │ 15   │
        │ 5    │
        └──────┘
        """
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.diff(n_pyexpr, null_behavior))

    def pct_change(self, n: int | IntoExprColumn = 1) -> Expr:
        """
        Computes percentage change between values.

        Percentage change (as fraction) between current element and most-recent
        non-null element at least `n` period(s) before the current element.

        Computes the change from the previous row by default.

        Parameters
        ----------
        n
            periods to shift for forming percent change.

        Notes
        -----
        Null values are preserved. If you're coming from pandas, this matches
        their ``fill_method=None`` behaviour.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [10, 11, 12, None, 12],
        ...     }
        ... )
        >>> df.with_columns(pl.col("a").pct_change().alias("pct_change"))
        shape: (5, 2)
        ┌──────┬────────────┐
        │ a    ┆ pct_change │
        │ ---  ┆ ---        │
        │ i64  ┆ f64        │
        ╞══════╪════════════╡
        │ 10   ┆ null       │
        │ 11   ┆ 0.1        │
        │ 12   ┆ 0.090909   │
        │ null ┆ null       │
        │ 12   ┆ null       │
        └──────┴────────────┘
        """
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.pct_change(n_pyexpr))

    def skew(self, *, bias: bool = True) -> Expr:
        r"""
        Compute the sample skewness of a data set.

        For normally distributed data, the skewness should be about zero. For
        unimodal continuous distributions, a skewness value greater than zero means
        that there is more weight in the right tail of the distribution. The
        function `skewtest` can be used to determine if the skewness value
        is close enough to zero, statistically speaking.


        See scipy.stats for more information.

        Parameters
        ----------
        bias : bool, optional
            If False, the calculations are corrected for statistical bias.

        Notes
        -----
        The sample skewness is computed as the Fisher-Pearson coefficient
        of skewness, i.e.

        .. math:: g_1=\frac{m_3}{m_2^{3/2}}

        where

        .. math:: m_i=\frac{1}{N}\sum_{n=1}^N(x[n]-\bar{x})^i

        is the biased sample :math:`i\texttt{th}` central moment, and
        :math:`\bar{x}` is
        the sample mean. If `bias` is False, the calculations are
        corrected for bias and the value computed is the adjusted
        Fisher-Pearson standardized moment coefficient, i.e.

        .. math::
            G_1 = \frac{k_3}{k_2^{3/2}} = \frac{\sqrt{N(N-1)}}{N-2}\frac{m_3}{m_2^{3/2}}

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 2, 1]})
        >>> df.select(pl.col("a").skew())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.343622 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.skew(bias))

    def kurtosis(self, *, fisher: bool = True, bias: bool = True) -> Expr:
        """
        Compute the kurtosis (Fisher or Pearson) of a dataset.

        Kurtosis is the fourth central moment divided by the square of the
        variance. If Fisher's definition is used, then 3.0 is subtracted from
        the result to give 0.0 for a normal distribution.
        If bias is False then the kurtosis is calculated using k statistics to
        eliminate bias coming from biased moment estimators.

        See scipy.stats for more information

        Parameters
        ----------
        fisher : bool, optional
            If True, Fisher's definition is used (normal ==> 0.0). If False,
            Pearson's definition is used (normal ==> 3.0).
        bias : bool, optional
            If False, the calculations are corrected for statistical bias.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 2, 1]})
        >>> df.select(pl.col("a").kurtosis())
        shape: (1, 1)
        ┌───────────┐
        │ a         │
        │ ---       │
        │ f64       │
        ╞═══════════╡
        │ -1.153061 │
        └───────────┘
        """
        return wrap_expr(self._pyexpr.kurtosis(fisher, bias))

    def clip(
        self,
        lower_bound: NumericLiteral | TemporalLiteral | IntoExprColumn | None = None,
        upper_bound: NumericLiteral | TemporalLiteral | IntoExprColumn | None = None,
    ) -> Expr:
        """
        Set values outside the given boundaries to the boundary value.

        Parameters
        ----------
        lower_bound
            Lower bound. Accepts expression input. Non-expression inputs are
            parsed as literals. Strings are parsed as column names.
        upper_bound
            Upper bound. Accepts expression input. Non-expression inputs are
            parsed as literals. Strings are parsed as column names.

        See Also
        --------
        when

        Notes
        -----
        This method only works for numeric and temporal columns. To clip other data
        types, consider writing a `when-then-otherwise` expression. See :func:`when`.

        Examples
        --------
        Specifying both a lower and upper bound:

        >>> df = pl.DataFrame({"a": [-50, 5, 50, None]})
        >>> df.with_columns(clip=pl.col("a").clip(1, 10))
        shape: (4, 2)
        ┌──────┬──────┐
        │ a    ┆ clip │
        │ ---  ┆ ---  │
        │ i64  ┆ i64  │
        ╞══════╪══════╡
        │ -50  ┆ 1    │
        │ 5    ┆ 5    │
        │ 50   ┆ 10   │
        │ null ┆ null │
        └──────┴──────┘

        Specifying only a single bound:

        >>> df.with_columns(clip=pl.col("a").clip(upper_bound=10))
        shape: (4, 2)
        ┌──────┬──────┐
        │ a    ┆ clip │
        │ ---  ┆ ---  │
        │ i64  ┆ i64  │
        ╞══════╪══════╡
        │ -50  ┆ -50  │
        │ 5    ┆ 5    │
        │ 50   ┆ 10   │
        │ null ┆ null │
        └──────┴──────┘

        Using columns as bounds:

        >>> df = pl.DataFrame(
        ...     {"a": [-50, 5, 50, None], "low": [10, 1, 0, 0], "up": [20, 4, 3, 2]}
        ... )
        >>> df.with_columns(clip=pl.col("a").clip("low", "up"))
        shape: (4, 4)
        ┌──────┬─────┬─────┬──────┐
        │ a    ┆ low ┆ up  ┆ clip │
        │ ---  ┆ --- ┆ --- ┆ ---  │
        │ i64  ┆ i64 ┆ i64 ┆ i64  │
        ╞══════╪═════╪═════╪══════╡
        │ -50  ┆ 10  ┆ 20  ┆ 10   │
        │ 5    ┆ 1   ┆ 4   ┆ 4    │
        │ 50   ┆ 0   ┆ 3   ┆ 3    │
        │ null ┆ 0   ┆ 2   ┆ null │
        └──────┴─────┴─────┴──────┘
        """
        if lower_bound is not None:
            lower_bound_pyexpr = parse_into_expression(lower_bound)
        else:
            lower_bound_pyexpr = None
        if upper_bound is not None:
            upper_bound_pyexpr = parse_into_expression(upper_bound)
        else:
            upper_bound_pyexpr = None
        return wrap_expr(self._pyexpr.clip(lower_bound_pyexpr, upper_bound_pyexpr))

    def lower_bound(self) -> Expr:
        """
        Calculate the lower bound.

        Returns a unit Series with the lowest value possible for the dtype of this
        expression.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 2, 1]})
        >>> df.select(pl.col("a").lower_bound())
        shape: (1, 1)
        ┌──────────────────────┐
        │ a                    │
        │ ---                  │
        │ i64                  │
        ╞══════════════════════╡
        │ -9223372036854775808 │
        └──────────────────────┘
        """
        return wrap_expr(self._pyexpr.lower_bound())

    def upper_bound(self) -> Expr:
        """
        Calculate the upper bound.

        Returns a unit Series with the highest value possible for the dtype of this
        expression.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 2, 1]})
        >>> df.select(pl.col("a").upper_bound())
        shape: (1, 1)
        ┌─────────────────────┐
        │ a                   │
        │ ---                 │
        │ i64                 │
        ╞═════════════════════╡
        │ 9223372036854775807 │
        └─────────────────────┘
        """
        return wrap_expr(self._pyexpr.upper_bound())

    def sign(self) -> Expr:
        """
        Compute the element-wise sign function on numeric types.

        The returned value is computed as follows:

        * -1 if x < 0.
        *  1 if x > 0.
        *  x otherwise (typically 0, but could be NaN if the input is).

        Null values are preserved as-is, and the dtype of the input is preserved.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-9.0, -0.0, 0.0, 4.0, float("nan"), None]})
        >>> df.select(pl.col.a.sign())
        shape: (6, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ -1.0 │
        │ -0.0 │
        │ 0.0  │
        │ 1.0  │
        │ NaN  │
        │ null │
        └──────┘
        """
        return wrap_expr(self._pyexpr.sign())

    def sin(self) -> Expr:
        """
        Compute the element-wise value for the sine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.0]})
        >>> df.select(pl.col("a").sin())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.sin())

    def cos(self) -> Expr:
        """
        Compute the element-wise value for the cosine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.0]})
        >>> df.select(pl.col("a").cos())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.cos())

    def tan(self) -> Expr:
        """
        Compute the element-wise value for the tangent.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").tan().round(2))
        shape: (1, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ 1.56 │
        └──────┘
        """
        return wrap_expr(self._pyexpr.tan())

    def cot(self) -> Expr:
        """
        Compute the element-wise value for the cotangent.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").cot().round(2))
        shape: (1, 1)
        ┌──────┐
        │ a    │
        │ ---  │
        │ f64  │
        ╞══════╡
        │ 0.64 │
        └──────┘
        """
        return wrap_expr(self._pyexpr.cot())

    def arcsin(self) -> Expr:
        """
        Compute the element-wise value for the inverse sine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").arcsin())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.570796 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.arcsin())

    def arccos(self) -> Expr:
        """
        Compute the element-wise value for the inverse cosine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [0.0]})
        >>> df.select(pl.col("a").arccos())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.570796 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.arccos())

    def arctan(self) -> Expr:
        """
        Compute the element-wise value for the inverse tangent.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").arctan())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.785398 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.arctan())

    def sinh(self) -> Expr:
        """
        Compute the element-wise value for the hyperbolic sine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").sinh())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.175201 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.sinh())

    def cosh(self) -> Expr:
        """
        Compute the element-wise value for the hyperbolic cosine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").cosh())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.543081 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.cosh())

    def tanh(self) -> Expr:
        """
        Compute the element-wise value for the hyperbolic tangent.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").tanh())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.761594 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.tanh())

    def arcsinh(self) -> Expr:
        """
        Compute the element-wise value for the inverse hyperbolic sine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").arcsinh())
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.881374 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.arcsinh())

    def arccosh(self) -> Expr:
        """
        Compute the element-wise value for the inverse hyperbolic cosine.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").arccosh())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arccosh())

    def arctanh(self) -> Expr:
        """
        Compute the element-wise value for the inverse hyperbolic tangent.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1.0]})
        >>> df.select(pl.col("a").arctanh())
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ inf │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arctanh())

    def degrees(self) -> Expr:
        """
        Convert from radians to degrees.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> import math
        >>> df = pl.DataFrame({"a": [x * math.pi for x in range(-4, 5)]})
        >>> df.select(pl.col("a").degrees())
        shape: (9, 1)
        ┌────────┐
        │ a      │
        │ ---    │
        │ f64    │
        ╞════════╡
        │ -720.0 │
        │ -540.0 │
        │ -360.0 │
        │ -180.0 │
        │ 0.0    │
        │ 180.0  │
        │ 360.0  │
        │ 540.0  │
        │ 720.0  │
        └────────┘
        """
        return wrap_expr(self._pyexpr.degrees())

    def radians(self) -> Expr:
        """
        Convert from degrees to radians.

        Returns
        -------
        Expr
            Expression of data type :class:`Float64`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-720, -540, -360, -180, 0, 180, 360, 540, 720]})
        >>> df.select(pl.col("a").radians())
        shape: (9, 1)
        ┌────────────┐
        │ a          │
        │ ---        │
        │ f64        │
        ╞════════════╡
        │ -12.566371 │
        │ -9.424778  │
        │ -6.283185  │
        │ -3.141593  │
        │ 0.0        │
        │ 3.141593   │
        │ 6.283185   │
        │ 9.424778   │
        │ 12.566371  │
        └────────────┘
        """
        return wrap_expr(self._pyexpr.radians())

    def reshape(self, dimensions: tuple[int, ...]) -> Expr:
        """
        Reshape this Expr to a flat column or an Array column.

        Parameters
        ----------
        dimensions
            Tuple of the dimension sizes. If -1 is used as the value for the
            first dimension, that dimension is inferred.
            Because the size of the Column may not be known in advance, it is
            only possible to use -1 for the first dimension.

        Returns
        -------
        Expr
            If a single dimension is given, results in an expression of the original
            data type.
            If a multiple dimensions are given, results in an expression of data type
            :class:`Array` with shape `dimensions`.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5, 6, 7, 8, 9]})
        >>> square = df.select(pl.col("foo").reshape((3, 3)))
        >>> square
        shape: (3, 1)
        ┌───────────────┐
        │ foo           │
        │ ---           │
        │ array[i64, 3] │
        ╞═══════════════╡
        │ [1, 2, 3]     │
        │ [4, 5, 6]     │
        │ [7, 8, 9]     │
        └───────────────┘
        >>> square.select(pl.col("foo").reshape((9,)))
        shape: (9, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        │ 4   │
        │ 5   │
        │ 6   │
        │ 7   │
        │ 8   │
        │ 9   │
        └─────┘

        See Also
        --------
        Expr.list.explode : Explode a list column.
        """
        return wrap_expr(self._pyexpr.reshape(dimensions))

    def shuffle(self, seed: int | None = None) -> Expr:
        """
        Shuffle the contents of this expression.

        Note this is shuffled independently of any other column or Expression. If you
        want each row to stay the same use df.sample(shuffle=True)

        Parameters
        ----------
        seed
            Seed for the random number generator. If set to None (default), a
            random seed is generated each time the shuffle is called.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").shuffle(seed=1))
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        │ 3   │
        │ 1   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.shuffle(seed))

    def sample(
        self,
        n: int | IntoExprColumn | None = None,
        *,
        fraction: float | IntoExprColumn | None = None,
        with_replacement: bool = False,
        shuffle: bool = False,
        seed: int | None = None,
    ) -> Expr:
        """
        Sample from this expression.

        Parameters
        ----------
        n
            Number of items to return. Cannot be used with `fraction`. Defaults to 1 if
            `fraction` is None.
        fraction
            Fraction of items to return. Cannot be used with `n`.
        with_replacement
            Allow values to be sampled more than once.
        shuffle
            Shuffle the order of sampled data points.
        seed
            Seed for the random number generator. If set to None (default), a
            random seed is generated for each sample operation.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").sample(fraction=1.0, with_replacement=True, seed=1))
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 3   │
        │ 3   │
        │ 1   │
        └─────┘
        """
        if n is not None and fraction is not None:
            msg = "cannot specify both `n` and `fraction`"
            raise ValueError(msg)

        if fraction is not None:
            fraction_pyexpr = parse_into_expression(fraction)
            return wrap_expr(
                self._pyexpr.sample_frac(
                    fraction_pyexpr, with_replacement, shuffle, seed
                )
            )

        if n is None:
            n = 1
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(
            self._pyexpr.sample_n(n_pyexpr, with_replacement, shuffle, seed)
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
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
    ) -> Expr:
        r"""
        Compute exponentially-weighted moving average.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        com
            Specify decay in terms of center of mass, :math:`\gamma`, with

                .. math::
                    \alpha = \frac{1}{1 + \gamma} \; \forall \; \gamma \geq 0
        span
            Specify decay in terms of span, :math:`\theta`, with

                .. math::
                    \alpha = \frac{2}{\theta + 1} \; \forall \; \theta \geq 1
        half_life
            Specify decay in terms of half-life, :math:`\tau`, with

                .. math::
                    \alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \tau } \right\} \;
                    \forall \; \tau > 0
        alpha
            Specify smoothing factor alpha directly, :math:`0 < \alpha \leq 1`.
        adjust
            Divide by decaying adjustment factor in beginning periods to account for
            imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights :math:`w_i = (1 - \alpha)^i`
                - When `adjust=False` the EW function is calculated
                  recursively by

                  .. math::
                    y_0 &= x_0 \\
                    y_t &= (1 - \alpha)y_{t - 1} + \alpha x_t
        min_samples
            Minimum number of observations in window required to have a value
            (otherwise result is null).
        ignore_nulls
            Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of :math:`x_0` and :math:`x_2` used in
                  calculating the final weighted average of
                  [:math:`x_0`, None, :math:`x_2`] are
                  :math:`(1-\alpha)^2` and :math:`1` if `adjust=True`, and
                  :math:`(1-\alpha)^2` and :math:`\alpha` if `adjust=False`.

                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  :math:`x_0` and :math:`x_2` used in calculating the final weighted
                  average of [:math:`x_0`, None, :math:`x_2`] are
                  :math:`1-\alpha` and :math:`1` if `adjust=True`,
                  and :math:`1-\alpha` and :math:`\alpha` if `adjust=False`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").ewm_mean(com=1, ignore_nulls=False))
        shape: (3, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.0      │
        │ 1.666667 │
        │ 2.428571 │
        └──────────┘
        """
        alpha = _prepare_alpha(com, span, half_life, alpha)
        return wrap_expr(
            self._pyexpr.ewm_mean(alpha, adjust, min_samples, ignore_nulls)
        )

    def ewm_mean_by(
        self,
        by: str | IntoExpr,
        *,
        half_life: str | timedelta,
    ) -> Expr:
        r"""
        Compute time-based exponentially weighted moving average.

        Given observations :math:`x_0, x_1, \ldots, x_{n-1}` at times
        :math:`t_0, t_1, \ldots, t_{n-1}`, the EWMA is calculated as

            .. math::

                y_0 &= x_0

                \alpha_i &= 1 - \exp \left\{ \frac{ -\ln(2)(t_i-t_{i-1}) }
                    { \tau } \right\}

                y_i &= \alpha_i x_i + (1 - \alpha_i) y_{i-1}; \quad i > 0

        where :math:`\tau` is the `half_life`.

        Parameters
        ----------
        by
            Times to calculate average by. Should be ``DateTime``, ``Date``, ``UInt64``,
            ``UInt32``, ``Int64``, or ``Int32`` data type.
        half_life
            Unit over which observation decays to half its value.

            Can be created either from a timedelta, or
            by using the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 day)
            - 1w    (1 week)
            - 1i    (1 index count)

            Or combine them:
            "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

            Note that `half_life` is treated as a constant duration - calendar
            durations such as months (or even days in the time-zone-aware case)
            are not supported, please express your duration in an approximately
            equivalent number of hours (e.g. '370h' instead of '1mo').

        Returns
        -------
        Expr
            Float32 if input is Float32, otherwise Float64.

        Examples
        --------
        >>> from datetime import date, timedelta
        >>> df = pl.DataFrame(
        ...     {
        ...         "values": [0, 1, 2, None, 4],
        ...         "times": [
        ...             date(2020, 1, 1),
        ...             date(2020, 1, 3),
        ...             date(2020, 1, 10),
        ...             date(2020, 1, 15),
        ...             date(2020, 1, 17),
        ...         ],
        ...     }
        ... ).sort("times")
        >>> df.with_columns(
        ...     result=pl.col("values").ewm_mean_by("times", half_life="4d"),
        ... )
        shape: (5, 3)
        ┌────────┬────────────┬──────────┐
        │ values ┆ times      ┆ result   │
        │ ---    ┆ ---        ┆ ---      │
        │ i64    ┆ date       ┆ f64      │
        ╞════════╪════════════╪══════════╡
        │ 0      ┆ 2020-01-01 ┆ 0.0      │
        │ 1      ┆ 2020-01-03 ┆ 0.292893 │
        │ 2      ┆ 2020-01-10 ┆ 1.492474 │
        │ null   ┆ 2020-01-15 ┆ null     │
        │ 4      ┆ 2020-01-17 ┆ 3.254508 │
        └────────┴────────────┴──────────┘
        """
        by_pyexpr = parse_into_expression(by)
        half_life = parse_as_duration_string(half_life)
        return wrap_expr(self._pyexpr.ewm_mean_by(by_pyexpr, half_life))

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def ewm_std(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        bias: bool = False,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Expr:
        r"""
        Compute exponentially-weighted moving standard deviation.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        com
            Specify decay in terms of center of mass, :math:`\gamma`, with

                .. math::
                    \alpha = \frac{1}{1 + \gamma} \; \forall \; \gamma \geq 0
        span
            Specify decay in terms of span, :math:`\theta`, with

                .. math::
                    \alpha = \frac{2}{\theta + 1} \; \forall \; \theta \geq 1
        half_life
            Specify decay in terms of half-life, :math:`\lambda`, with

                .. math::
                    \alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \lambda } \right\} \;
                    \forall \; \lambda > 0
        alpha
            Specify smoothing factor alpha directly, :math:`0 < \alpha \leq 1`.
        adjust
            Divide by decaying adjustment factor in beginning periods to account for
            imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights :math:`w_i = (1 - \alpha)^i`
                - When `adjust=False` the EW function is calculated
                  recursively by

                  .. math::
                    y_0 &= x_0 \\
                    y_t &= (1 - \alpha)y_{t - 1} + \alpha x_t
        bias
            When `bias=False`, apply a correction to make the estimate statistically
            unbiased.
        min_samples
            Minimum number of observations in window required to have a value
            (otherwise result is null).
        ignore_nulls
            Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of :math:`x_0` and :math:`x_2` used in
                  calculating the final weighted average of
                  [:math:`x_0`, None, :math:`x_2`] are
                  :math:`(1-\alpha)^2` and :math:`1` if `adjust=True`, and
                  :math:`(1-\alpha)^2` and :math:`\alpha` if `adjust=False`.

                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  :math:`x_0` and :math:`x_2` used in calculating the final weighted
                  average of [:math:`x_0`, None, :math:`x_2`] are
                  :math:`1-\alpha` and :math:`1` if `adjust=True`,
                  and :math:`1-\alpha` and :math:`\alpha` if `adjust=False`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").ewm_std(com=1, ignore_nulls=False))
        shape: (3, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.0      │
        │ 0.707107 │
        │ 0.963624 │
        └──────────┘
        """
        alpha = _prepare_alpha(com, span, half_life, alpha)
        return wrap_expr(
            self._pyexpr.ewm_std(alpha, adjust, bias, min_samples, ignore_nulls)
        )

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def ewm_var(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        bias: bool = False,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Expr:
        r"""
        Compute exponentially-weighted moving variance.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        com
            Specify decay in terms of center of mass, :math:`\gamma`, with

                .. math::
                    \alpha = \frac{1}{1 + \gamma} \; \forall \; \gamma \geq 0
        span
            Specify decay in terms of span, :math:`\theta`, with

                .. math::
                    \alpha = \frac{2}{\theta + 1} \; \forall \; \theta \geq 1
        half_life
            Specify decay in terms of half-life, :math:`\lambda`, with

                .. math::
                    \alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \lambda } \right\} \;
                    \forall \; \lambda > 0
        alpha
            Specify smoothing factor alpha directly, :math:`0 < \alpha \leq 1`.
        adjust
            Divide by decaying adjustment factor in beginning periods to account for
            imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights :math:`w_i = (1 - \alpha)^i`
                - When `adjust=False` the EW function is calculated
                  recursively by

                  .. math::
                    y_0 &= x_0 \\
                    y_t &= (1 - \alpha)y_{t - 1} + \alpha x_t
        bias
            When `bias=False`, apply a correction to make the estimate statistically
            unbiased.
        min_samples
            Minimum number of observations in window required to have a value
            (otherwise result is null).
        ignore_nulls
            Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of :math:`x_0` and :math:`x_2` used in
                  calculating the final weighted average of
                  [:math:`x_0`, None, :math:`x_2`] are
                  :math:`(1-\alpha)^2` and :math:`1` if `adjust=True`, and
                  :math:`(1-\alpha)^2` and :math:`\alpha` if `adjust=False`.

                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  :math:`x_0` and :math:`x_2` used in calculating the final weighted
                  average of [:math:`x_0`, None, :math:`x_2`] are
                  :math:`1-\alpha` and :math:`1` if `adjust=True`,
                  and :math:`1-\alpha` and :math:`\alpha` if `adjust=False`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").ewm_var(com=1, ignore_nulls=False))
        shape: (3, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.0      │
        │ 0.5      │
        │ 0.928571 │
        └──────────┘
        """
        alpha = _prepare_alpha(com, span, half_life, alpha)
        return wrap_expr(
            self._pyexpr.ewm_var(alpha, adjust, bias, min_samples, ignore_nulls)
        )

    def extend_constant(self, value: IntoExpr, n: int | IntoExprColumn) -> Expr:
        """
        Extremely fast method for extending the Series with 'n' copies of a value.

        Parameters
        ----------
        value
            A constant literal value or a unit expression with which to extend the
            expression result Series; can pass None to extend with nulls.
        n
            The number of additional values that will be added.

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1, 2, 3]})
        >>> df.select((pl.col("values") - 1).extend_constant(99, n=2))
        shape: (5, 1)
        ┌────────┐
        │ values │
        │ ---    │
        │ i64    │
        ╞════════╡
        │ 0      │
        │ 1      │
        │ 2      │
        │ 99     │
        │ 99     │
        └────────┘
        """
        value_pyexpr = parse_into_expression(value, str_as_lit=True)
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.extend_constant(value_pyexpr, n_pyexpr))

    def value_counts(
        self,
        *,
        sort: bool = False,
        parallel: bool = False,
        name: str | None = None,
        normalize: bool = False,
    ) -> Expr:
        """
        Count the occurrence of unique values.

        Parameters
        ----------
        sort
            Sort the output by count, in descending order.
            If set to `False` (default), the order is non-deterministic.
        parallel
            Execute the computation in parallel.

            .. note::
                This option should likely *not* be enabled in a `group_by` context,
                as the computation will already be parallelized per group.
        name
            Give the resulting count column a specific name; if `normalize` is
            True this defaults to "proportion", otherwise defaults to "count".
        normalize
            If True, the count is returned as the relative frequency of unique
            values normalized to 1.0.

        Returns
        -------
        Expr
            Expression of type :class:`Struct`, mapping unique values to their
            count (or proportion).

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"color": ["red", "blue", "red", "green", "blue", "blue"]}
        ... )
        >>> df_count = df.select(pl.col("color").value_counts())
        >>> df_count  # doctest: +IGNORE_RESULT
        shape: (3, 1)
        ┌─────────────┐
        │ color       │
        │ ---         │
        │ struct[2]   │
        ╞═════════════╡
        │ {"green",1} │
        │ {"blue",3}  │
        │ {"red",2}   │
        └─────────────┘

        >>> df_count.unnest("color")  # doctest: +IGNORE_RESULT
        shape: (3, 2)
        ┌───────┬───────┐
        │ color ┆ count │
        │ ---   ┆ ---   │
        │ str   ┆ u32   │
        ╞═══════╪═══════╡
        │ green ┆ 1     │
        │ blue  ┆ 3     │
        │ red   ┆ 2     │
        └───────┴───────┘

        Sort the output by (descending) count, customize the field name,
        and normalize the count to its relative proportion (of 1.0).

        >>> df_count = df.select(
        ...     pl.col("color").value_counts(
        ...         name="fraction",
        ...         normalize=True,
        ...         sort=True,
        ...     )
        ... )
        >>> df_count
        shape: (3, 1)
        ┌────────────────────┐
        │ color              │
        │ ---                │
        │ struct[2]          │
        ╞════════════════════╡
        │ {"blue",0.5}       │
        │ {"red",0.333333}   │
        │ {"green",0.166667} │
        └────────────────────┘

        >>> df_count.unnest("color")
        shape: (3, 2)
        ┌───────┬──────────┐
        │ color ┆ fraction │
        │ ---   ┆ ---      │
        │ str   ┆ f64      │
        ╞═══════╪══════════╡
        │ blue  ┆ 0.5      │
        │ red   ┆ 0.333333 │
        │ green ┆ 0.166667 │
        └───────┴──────────┘

        Note that `group_by` can be used to generate counts.

        >>> df.group_by("color").len()  # doctest: +IGNORE_RESULT
        shape: (3, 2)
        ┌───────┬─────┐
        │ color ┆ len │
        │ ---   ┆ --- │
        │ str   ┆ u32 │
        ╞═══════╪═════╡
        │ red   ┆ 2   │
        │ green ┆ 1   │
        │ blue  ┆ 3   │
        └───────┴─────┘

        To add counts as a new column `pl.len()` can be used as a window function.

        >>> df.with_columns(pl.len().over("color"))
        shape: (6, 2)
        ┌───────┬─────┐
        │ color ┆ len │
        │ ---   ┆ --- │
        │ str   ┆ u32 │
        ╞═══════╪═════╡
        │ red   ┆ 2   │
        │ blue  ┆ 3   │
        │ red   ┆ 2   │
        │ green ┆ 1   │
        │ blue  ┆ 3   │
        │ blue  ┆ 3   │
        └───────┴─────┘

        >>> df.with_columns((pl.len().over("color") / pl.len()).alias("fraction"))
        shape: (6, 2)
        ┌───────┬──────────┐
        │ color ┆ fraction │
        │ ---   ┆ ---      │
        │ str   ┆ f64      │
        ╞═══════╪══════════╡
        │ red   ┆ 0.333333 │
        │ blue  ┆ 0.5      │
        │ red   ┆ 0.333333 │
        │ green ┆ 0.166667 │
        │ blue  ┆ 0.5      │
        │ blue  ┆ 0.5      │
        └───────┴──────────┘
        """
        name = name or ("proportion" if normalize else "count")
        return wrap_expr(self._pyexpr.value_counts(sort, parallel, name, normalize))

    def unique_counts(self) -> Expr:
        """
        Return a count of the unique values in the order of appearance.

        This method differs from `value_counts` in that it does not return the
        values, only the counts and might be faster

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "id": ["a", "b", "b", "c", "c", "c"],
        ...     }
        ... )
        >>> df.select(pl.col("id").unique_counts())
        shape: (3, 1)
        ┌─────┐
        │ id  │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘

        Note that `group_by` can be used to generate counts.

        >>> df.group_by("id", maintain_order=True).len().select("len")
        shape: (3, 1)
        ┌─────┐
        │ len │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘

        To add counts as a new column `pl.len()` can be used as a window function.

        >>> df.with_columns(pl.len().over("id"))
        shape: (6, 2)
        ┌─────┬─────┐
        │ id  ┆ len │
        │ --- ┆ --- │
        │ str ┆ u32 │
        ╞═════╪═════╡
        │ a   ┆ 1   │
        │ b   ┆ 2   │
        │ b   ┆ 2   │
        │ c   ┆ 3   │
        │ c   ┆ 3   │
        │ c   ┆ 3   │
        └─────┴─────┘
        """
        return wrap_expr(self._pyexpr.unique_counts())

    def log(self, base: float | IntoExpr = math.e) -> Expr:
        """
        Compute the logarithm to a given base.

        Parameters
        ----------
        base
            Given base, defaults to `e`

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").log(base=2))
        shape: (3, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.0      │
        │ 1.0      │
        │ 1.584963 │
        └──────────┘
        """
        base_pyexpr = parse_into_expression(base)
        return wrap_expr(self._pyexpr.log(base_pyexpr))

    def log1p(self) -> Expr:
        """
        Compute the natural logarithm of each element plus one.

        This computes `log(1 + x)` but is more numerically stable for `x` close to zero.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").log1p())
        shape: (3, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.693147 │
        │ 1.098612 │
        │ 1.386294 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.log1p())

    def entropy(self, base: float = math.e, *, normalize: bool = True) -> Expr:
        """
        Computes the entropy.

        Uses the formula `-sum(pk * log(pk))` where `pk` are discrete probabilities.

        Parameters
        ----------
        base
            Given base, defaults to `e`
        normalize
            Normalize pk if it doesn't sum to 1.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select(pl.col("a").entropy(base=2))
        shape: (1, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 1.459148 │
        └──────────┘
        >>> df.select(pl.col("a").entropy(base=2, normalize=False))
        shape: (1, 1)
        ┌───────────┐
        │ a         │
        │ ---       │
        │ f64       │
        ╞═══════════╡
        │ -6.754888 │
        └───────────┘
        """
        return wrap_expr(self._pyexpr.entropy(base, normalize))

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def cumulative_eval(self, expr: Expr, *, min_samples: int = 1) -> Expr:
        """
        Run an expression over a sliding window that increases `1` slot every iteration.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        expr
            Expression to evaluate
        min_samples
            Number of valid values there should be in the window before the expression
            is evaluated. valid values = `length - null_count`

        Warnings
        --------
        This can be really slow as it can have `O(n^2)` complexity. Don't use this
        for operations that visit all elements.

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1, 2, 3, 4, 5]})
        >>> df.select(
        ...     [
        ...         pl.col("values").cumulative_eval(
        ...             pl.element().first() - pl.element().last() ** 2
        ...         )
        ...     ]
        ... )
        shape: (5, 1)
        ┌────────┐
        │ values │
        │ ---    │
        │ i64    │
        ╞════════╡
        │ 0      │
        │ -3     │
        │ -8     │
        │ -15    │
        │ -24    │
        └────────┘
        """
        return wrap_expr(self._pyexpr.cumulative_eval(expr._pyexpr, min_samples))

    def set_sorted(self, *, descending: bool = False) -> Expr:
        """
        Flags the expression as 'sorted'.

        Enables downstream code to user fast paths for sorted arrays.

        Parameters
        ----------
        descending
            Whether the `Series` order is descending.

        Warnings
        --------
        This can lead to incorrect results if the data is NOT sorted!!
        Use with care!

        Examples
        --------
        >>> df = pl.DataFrame({"values": [1, 2, 3]})
        >>> df.select(pl.col("values").set_sorted().max())
        shape: (1, 1)
        ┌────────┐
        │ values │
        │ ---    │
        │ i64    │
        ╞════════╡
        │ 3      │
        └────────┘
        """
        return wrap_expr(self._pyexpr.set_sorted_flag(descending))

    @deprecated(
        "`Expr.shrink_dtype` is deprecated and is a no-op; use `Series.shrink_dtype` instead."
    )
    def shrink_dtype(self) -> Expr:
        """
        Shrink numeric columns to the minimal required datatype.

        Shrink to the dtype needed to fit the extrema of this [`Series`].
        This can be used to reduce memory pressure.

        .. versionchanged:: 1.33.0
            Deprecated and turned into a no-op. The operation does not match the
            Polars data-model during lazy execution since the output datatype
            cannot be known without inspecting the data.

            Use `Series.shrink_dtype` instead.

        Examples
        --------
        >>> pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": [1, 2, 2 << 32],
        ...         "c": [-1, 2, 1 << 30],
        ...         "d": [-112, 2, 112],
        ...         "e": [-112, 2, 129],
        ...         "f": ["a", "b", "c"],
        ...         "g": [0.1, 1.32, 0.12],
        ...         "h": [True, None, False],
        ...     }
        ... ).select(pl.all().shrink_dtype())  # doctest: +SKIP
        shape: (3, 8)
        ┌─────┬────────────┬────────────┬──────┬──────┬─────┬──────┬───────┐
        │ a   ┆ b          ┆ c          ┆ d    ┆ e    ┆ f   ┆ g    ┆ h     │
        │ --- ┆ ---        ┆ ---        ┆ ---  ┆ ---  ┆ --- ┆ ---  ┆ ---   │
        │ i8  ┆ i64        ┆ i32        ┆ i8   ┆ i16  ┆ str ┆ f32  ┆ bool  │
        ╞═════╪════════════╪════════════╪══════╪══════╪═════╪══════╪═══════╡
        │ 1   ┆ 1          ┆ -1         ┆ -112 ┆ -112 ┆ a   ┆ 0.1  ┆ true  │
        │ 2   ┆ 2          ┆ 2          ┆ 2    ┆ 2    ┆ b   ┆ 1.32 ┆ null  │
        │ 3   ┆ 8589934592 ┆ 1073741824 ┆ 112  ┆ 129  ┆ c   ┆ 0.12 ┆ false │
        └─────┴────────────┴────────────┴──────┴──────┴─────┴──────┴───────┘
        """
        return self

    @unstable()
    def hist(
        self,
        bins: IntoExpr | None = None,
        *,
        bin_count: int | None = None,
        include_category: bool = False,
        include_breakpoint: bool = False,
    ) -> Expr:
        """
        Bin values into buckets and count their occurrences.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        bins
            Bin edges. If None given, we determine the edges based on the data.
        bin_count
            If `bins` is not provided, `bin_count` uniform bins are created that fully
            encompass the data.
        include_breakpoint
            Include a column that indicates the upper breakpoint.
        include_category
            Include a column that shows the intervals as categories.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 3, 8, 8, 2, 1, 3]})
        >>> df.select(pl.col("a").hist(bins=[1, 2, 3]))
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 3   │
        │ 2   │
        └─────┘
        >>> df.select(
        ...     pl.col("a").hist(
        ...         bins=[1, 2, 3], include_breakpoint=True, include_category=True
        ...     )
        ... )
        shape: (2, 1)
        ┌──────────────────────┐
        │ a                    │
        │ ---                  │
        │ struct[3]            │
        ╞══════════════════════╡
        │ {2.0,"[1.0, 2.0]",3} │
        │ {3.0,"(2.0, 3.0]",2} │
        └──────────────────────┘
        """
        if bins is not None:
            if isinstance(bins, list):
                bins = pl.Series(bins)
            bins_pyexpr = parse_into_expression(bins)
        else:
            bins_pyexpr = None
        return wrap_expr(
            self._pyexpr.hist(
                bins_pyexpr, bin_count, include_category, include_breakpoint
            )
        )

    def replace(
        self,
        old: IntoExpr | Sequence[Any] | Mapping[Any, Any],
        new: IntoExpr | Sequence[Any] | NoDefault = no_default,
        *,
        default: IntoExpr | NoDefault = no_default,
        return_dtype: PolarsDataType | None = None,
    ) -> Expr:
        """
        Replace the given values by different values of the same data type.

        Parameters
        ----------
        old
            Value or sequence of values to replace.
            Accepts expression input. Sequences are parsed as Series,
            other non-expression inputs are parsed as literals.
            Also accepts a mapping of values to their replacement as syntactic sugar for
            `replace(old=Series(mapping.keys()), new=Series(mapping.values()))`.
        new
            Value or sequence of values to replace by.
            Accepts expression input. Sequences are parsed as Series,
            other non-expression inputs are parsed as literals.
            Length must match the length of `old` or have length 1.

        default
            Set values that were not replaced to this value.
            Defaults to keeping the original value.
            Accepts expression input. Non-expression inputs are parsed as literals.

            .. deprecated:: 1.0.0
                Use :meth:`replace_strict` instead to set a default while replacing
                values.

        return_dtype
            The data type of the resulting expression. If set to `None` (default),
            the data type of the original column is preserved.

            .. deprecated:: 1.0.0
                Use :meth:`replace_strict` instead to set a return data type while
                replacing values, or explicitly call :meth:`cast` on the output.

        See Also
        --------
        replace_strict
        str.replace

        Notes
        -----
        The global string cache must be enabled when replacing categorical values.

        Examples
        --------
        Replace a single value by another value. Values that were not replaced remain
        unchanged.

        >>> df = pl.DataFrame({"a": [1, 2, 2, 3]})
        >>> df.with_columns(replaced=pl.col("a").replace(2, 100))
        shape: (4, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ i64 ┆ i64      │
        ╞═════╪══════════╡
        │ 1   ┆ 1        │
        │ 2   ┆ 100      │
        │ 2   ┆ 100      │
        │ 3   ┆ 3        │
        └─────┴──────────┘

        Replace multiple values by passing sequences to the `old` and `new` parameters.

        >>> df.with_columns(replaced=pl.col("a").replace([2, 3], [100, 200]))
        shape: (4, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ i64 ┆ i64      │
        ╞═════╪══════════╡
        │ 1   ┆ 1        │
        │ 2   ┆ 100      │
        │ 2   ┆ 100      │
        │ 3   ┆ 200      │
        └─────┴──────────┘

        Passing a mapping with replacements is also supported as syntactic sugar.

        >>> mapping = {2: 100, 3: 200}
        >>> df.with_columns(replaced=pl.col("a").replace(mapping))
        shape: (4, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ i64 ┆ i64      │
        ╞═════╪══════════╡
        │ 1   ┆ 1        │
        │ 2   ┆ 100      │
        │ 2   ┆ 100      │
        │ 3   ┆ 200      │
        └─────┴──────────┘

        The original data type is preserved when replacing by values of a different
        data type. Use :meth:`replace_strict` to replace and change the return data
        type.

        >>> df = pl.DataFrame({"a": ["x", "y", "z"]})
        >>> mapping = {"x": 1, "y": 2, "z": 3}
        >>> df.with_columns(replaced=pl.col("a").replace(mapping))
        shape: (3, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ str ┆ str      │
        ╞═════╪══════════╡
        │ x   ┆ 1        │
        │ y   ┆ 2        │
        │ z   ┆ 3        │
        └─────┴──────────┘

        Expression input is supported.

        >>> df = pl.DataFrame({"a": [1, 2, 2, 3], "b": [1.5, 2.5, 5.0, 1.0]})
        >>> df.with_columns(
        ...     replaced=pl.col("a").replace(
        ...         old=pl.col("a").max(),
        ...         new=pl.col("b").sum(),
        ...     )
        ... )
        shape: (4, 3)
        ┌─────┬─────┬──────────┐
        │ a   ┆ b   ┆ replaced │
        │ --- ┆ --- ┆ ---      │
        │ i64 ┆ f64 ┆ i64      │
        ╞═════╪═════╪══════════╡
        │ 1   ┆ 1.5 ┆ 1        │
        │ 2   ┆ 2.5 ┆ 2        │
        │ 2   ┆ 5.0 ┆ 2        │
        │ 3   ┆ 1.0 ┆ 10       │
        └─────┴─────┴──────────┘
        """
        if return_dtype is not None:
            issue_deprecation_warning(
                "the `return_dtype` parameter for `replace` is deprecated."
                " Use `replace_strict` instead to set a return data type while replacing values.",
                version="1.0.0",
            )
        if default is not no_default:
            issue_deprecation_warning(
                "the `default` parameter for `replace` is deprecated."
                " Use `replace_strict` instead to set a default while replacing values.",
                version="1.0.0",
            )
            return self.replace_strict(
                old, new, default=default, return_dtype=return_dtype
            )

        if new is no_default:
            if not isinstance(old, Mapping):
                msg = (
                    "`new` argument is required if `old` argument is not a Mapping type"
                )
                raise TypeError(msg)
            new = list(old.values())
            old = list(old.keys())
        else:
            if isinstance(old, Sequence) and not isinstance(old, (str, pl.Series)):
                old = pl.Series(old)
            if isinstance(new, Sequence) and not isinstance(new, (str, pl.Series)):
                new = pl.Series(new)

        old_pyexpr = parse_into_expression(old, str_as_lit=True)  # type: ignore[arg-type]
        new_pyexpr = parse_into_expression(new, str_as_lit=True)

        result = wrap_expr(self._pyexpr.replace(old_pyexpr, new_pyexpr))

        if return_dtype is not None:
            result = result.cast(return_dtype)

        return result

    def replace_strict(
        self,
        old: IntoExpr | Sequence[Any] | Mapping[Any, Any],
        new: IntoExpr | Sequence[Any] | NoDefault = no_default,
        *,
        default: IntoExpr | NoDefault = no_default,
        return_dtype: PolarsDataType | pl.DataTypeExpr | None = None,
    ) -> Expr:
        """
        Replace all values by different values.

        Parameters
        ----------
        old
            Value or sequence of values to replace.
            Accepts expression input. Sequences are parsed as Series,
            other non-expression inputs are parsed as literals.
            Also accepts a mapping of values to their replacement as syntactic sugar for
            `replace_strict(old=Series(mapping.keys()), new=Series(mapping.values()))`.
        new
            Value or sequence of values to replace by.
            Accepts expression input. Sequences are parsed as Series,
            other non-expression inputs are parsed as literals.
            Length must match the length of `old` or have length 1.
        default
            Set values that were not replaced to this value. If no default is specified,
            (default), an error is raised if any values were not replaced.
            Accepts expression input. Non-expression inputs are parsed as literals.
        return_dtype
            The data type of the resulting expression. If set to `None` (default),
            the data type is determined automatically based on the other inputs.

        Raises
        ------
        InvalidOperationError
            If any non-null values in the original column were not replaced, and no
            `default` was specified.

        See Also
        --------
        replace
        str.replace

        Notes
        -----
        The global string cache must be enabled when replacing categorical values.

        Examples
        --------
        Replace values by passing sequences to the `old` and `new` parameters.

        >>> df = pl.DataFrame({"a": [1, 2, 2, 3]})
        >>> df.with_columns(
        ...     replaced=pl.col("a").replace_strict([1, 2, 3], [100, 200, 300])
        ... )
        shape: (4, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ i64 ┆ i64      │
        ╞═════╪══════════╡
        │ 1   ┆ 100      │
        │ 2   ┆ 200      │
        │ 2   ┆ 200      │
        │ 3   ┆ 300      │
        └─────┴──────────┘

        Passing a mapping with replacements is also supported as syntactic sugar.

        >>> mapping = {1: 100, 2: 200, 3: 300}
        >>> df.with_columns(replaced=pl.col("a").replace_strict(mapping))
        shape: (4, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ i64 ┆ i64      │
        ╞═════╪══════════╡
        │ 1   ┆ 100      │
        │ 2   ┆ 200      │
        │ 2   ┆ 200      │
        │ 3   ┆ 300      │
        └─────┴──────────┘

        By default, an error is raised if any non-null values were not replaced.
        Specify a default to set all values that were not matched.

        >>> mapping = {2: 200, 3: 300}
        >>> df.with_columns(
        ...     replaced=pl.col("a").replace_strict(mapping)
        ... )  # doctest: +SKIP
        Traceback (most recent call last):
        ...
        polars.exceptions.InvalidOperationError: incomplete mapping specified for `replace_strict`
        >>> df.with_columns(replaced=pl.col("a").replace_strict(mapping, default=-1))
        shape: (4, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ i64 ┆ i64      │
        ╞═════╪══════════╡
        │ 1   ┆ -1       │
        │ 2   ┆ 200      │
        │ 2   ┆ 200      │
        │ 3   ┆ 300      │
        └─────┴──────────┘

        Replacing by values of a different data type sets the return type based on
        a combination of the `new` data type and the `default` data type.

        >>> df = pl.DataFrame({"a": ["x", "y", "z"]})
        >>> mapping = {"x": 1, "y": 2, "z": 3}
        >>> df.with_columns(replaced=pl.col("a").replace_strict(mapping))
        shape: (3, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ str ┆ i64      │
        ╞═════╪══════════╡
        │ x   ┆ 1        │
        │ y   ┆ 2        │
        │ z   ┆ 3        │
        └─────┴──────────┘
        >>> df.with_columns(replaced=pl.col("a").replace_strict(mapping, default="x"))
        shape: (3, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ str ┆ str      │
        ╞═════╪══════════╡
        │ x   ┆ 1        │
        │ y   ┆ 2        │
        │ z   ┆ 3        │
        └─────┴──────────┘

        Set the `return_dtype` parameter to control the resulting data type directly.

        >>> df.with_columns(
        ...     replaced=pl.col("a").replace_strict(mapping, return_dtype=pl.UInt8)
        ... )
        shape: (3, 2)
        ┌─────┬──────────┐
        │ a   ┆ replaced │
        │ --- ┆ ---      │
        │ str ┆ u8       │
        ╞═════╪══════════╡
        │ x   ┆ 1        │
        │ y   ┆ 2        │
        │ z   ┆ 3        │
        └─────┴──────────┘

        Expression input is supported for all parameters.

        >>> df = pl.DataFrame({"a": [1, 2, 2, 3], "b": [1.5, 2.5, 5.0, 1.0]})
        >>> df.with_columns(
        ...     replaced=pl.col("a").replace_strict(
        ...         old=pl.col("a").max(),
        ...         new=pl.col("b").sum(),
        ...         default=pl.col("b"),
        ...     )
        ... )
        shape: (4, 3)
        ┌─────┬─────┬──────────┐
        │ a   ┆ b   ┆ replaced │
        │ --- ┆ --- ┆ ---      │
        │ i64 ┆ f64 ┆ f64      │
        ╞═════╪═════╪══════════╡
        │ 1   ┆ 1.5 ┆ 1.5      │
        │ 2   ┆ 2.5 ┆ 2.5      │
        │ 2   ┆ 5.0 ┆ 5.0      │
        │ 3   ┆ 1.0 ┆ 10.0     │
        └─────┴─────┴──────────┘
        """  # noqa: W505
        if new is no_default:
            if not isinstance(old, Mapping):
                msg = (
                    "`new` argument is required if `old` argument is not a Mapping type"
                )
                raise TypeError(msg)
            new = list(old.values())
            old = list(old.keys())

        old_pyexpr = parse_into_expression(old, str_as_lit=True)  # type: ignore[arg-type]
        new_pyexpr = parse_into_expression(new, str_as_lit=True)  # type: ignore[arg-type]

        dtype_pyexpr: plr.PyDataTypeExpr | None = None
        if return_dtype is not None:
            dtype_pyexpr = parse_into_datatype_expr(return_dtype)._pydatatype_expr
        else:
            dtype_pyexpr = None

        default_pyexpr = (
            None
            if default is no_default
            else parse_into_expression(default, str_as_lit=True)
        )

        return wrap_expr(
            self._pyexpr.replace_strict(
                old_pyexpr, new_pyexpr, default_pyexpr, dtype_pyexpr
            )
        )

    def bitwise_count_ones(self) -> Expr:
        """Evaluate the number of set bits."""
        return wrap_expr(self._pyexpr.bitwise_count_ones())

    def bitwise_count_zeros(self) -> Expr:
        """Evaluate the number of unset bits."""
        return wrap_expr(self._pyexpr.bitwise_count_zeros())

    def bitwise_leading_ones(self) -> Expr:
        """Evaluate the number most-significant set bits before seeing an unset bit."""
        return wrap_expr(self._pyexpr.bitwise_leading_ones())

    def bitwise_leading_zeros(self) -> Expr:
        """Evaluate the number most-significant unset bits before seeing a set bit."""
        return wrap_expr(self._pyexpr.bitwise_leading_zeros())

    def bitwise_trailing_ones(self) -> Expr:
        """Evaluate the number least-significant set bits before seeing an unset bit."""
        return wrap_expr(self._pyexpr.bitwise_trailing_ones())

    def bitwise_trailing_zeros(self) -> Expr:
        """Evaluate the number least-significant unset bits before seeing a set bit."""
        return wrap_expr(self._pyexpr.bitwise_trailing_zeros())

    def bitwise_and(self) -> Expr:
        """Perform an aggregation of bitwise ANDs.

        Examples
        --------
        >>> df = pl.DataFrame({"n": [-1, 0, 1]})
        >>> df.select(pl.col("n").bitwise_and())
        shape: (1, 1)
        ┌─────┐
        │ n   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 0   │
        └─────┘
        >>> df = pl.DataFrame(
        ...     {"grouper": ["a", "a", "a", "b", "b"], "n": [-1, 0, 1, -1, 1]}
        ... )
        >>> df.group_by("grouper", maintain_order=True).agg(pl.col("n").bitwise_and())
        shape: (2, 2)
        ┌─────────┬─────┐
        │ grouper ┆ n   │
        │ ---     ┆ --- │
        │ str     ┆ i64 │
        ╞═════════╪═════╡
        │ a       ┆ 0   │
        │ b       ┆ 1   │
        └─────────┴─────┘
        """
        return wrap_expr(self._pyexpr.bitwise_and())

    def bitwise_or(self) -> Expr:
        """Perform an aggregation of bitwise ORs.

        Examples
        --------
        >>> df = pl.DataFrame({"n": [-1, 0, 1]})
        >>> df.select(pl.col("n").bitwise_or())
        shape: (1, 1)
        ┌─────┐
        │ n   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ -1  │
        └─────┘
        >>> df = pl.DataFrame(
        ...     {"grouper": ["a", "a", "a", "b", "b"], "n": [-1, 0, 1, -1, 1]}
        ... )
        >>> df.group_by("grouper", maintain_order=True).agg(pl.col("n").bitwise_or())
        shape: (2, 2)
        ┌─────────┬─────┐
        │ grouper ┆ n   │
        │ ---     ┆ --- │
        │ str     ┆ i64 │
        ╞═════════╪═════╡
        │ a       ┆ -1  │
        │ b       ┆ -1  │
        └─────────┴─────┘
        """
        return wrap_expr(self._pyexpr.bitwise_or())

    def bitwise_xor(self) -> Expr:
        """Perform an aggregation of bitwise XORs.

        Examples
        --------
        >>> df = pl.DataFrame({"n": [-1, 0, 1]})
        >>> df.select(pl.col("n").bitwise_xor())
        shape: (1, 1)
        ┌─────┐
        │ n   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ -2  │
        └─────┘
        >>> df = pl.DataFrame(
        ...     {"grouper": ["a", "a", "a", "b", "b"], "n": [-1, 0, 1, -1, 1]}
        ... )
        >>> df.group_by("grouper", maintain_order=True).agg(pl.col("n").bitwise_xor())
        shape: (2, 2)
        ┌─────────┬─────┐
        │ grouper ┆ n   │
        │ ---     ┆ --- │
        │ str     ┆ i64 │
        ╞═════════╪═════╡
        │ a       ┆ -2  │
        │ b       ┆ -2  │
        └─────────┴─────┘
        """
        return wrap_expr(self._pyexpr.bitwise_xor())

    @deprecated(
        "`register_plugin` is deprecated; "
        "use `polars.plugins.register_plugin_function` instead."
    )
    def register_plugin(
        self,
        *,
        lib: str,
        symbol: str,
        args: list[IntoExpr] | None = None,
        kwargs: dict[Any, Any] | None = None,
        is_elementwise: bool = False,
        input_wildcard_expansion: bool = False,
        returns_scalar: bool = False,
        cast_to_supertypes: bool = False,
        pass_name_to_apply: bool = False,
        changes_length: bool = False,
    ) -> Expr:
        """
        Register a plugin function.

        .. deprecated:: 0.20.16
            Use :func:`polars.plugins.register_plugin_function` instead.

        See the `user guide <https://docs.pola.rs/user-guide/plugins/>`_
        for more information about plugins.

        Warnings
        --------
        This method is deprecated. Use the new `polars.plugins.register_plugin_function`
        function instead.

        This is highly unsafe as this will call the C function loaded by
        `lib::symbol`.

        The parameters you set dictate how Polars will handle the function.
        Make sure they are correct!

        Parameters
        ----------
        lib
            Library to load.
        symbol
            Function to load.
        args
            Arguments (other than self) passed to this function.
            These arguments have to be of type Expression.
        kwargs
            Non-expression arguments. They must be JSON serializable.
        is_elementwise
            If the function only operates on scalars
            this will trigger fast paths.
        input_wildcard_expansion
            Expand expressions as input of this function.
        returns_scalar
            Automatically explode on unit length if it ran as final aggregation.
            this is the case for aggregations like `sum`, `min`, `covariance` etc.
        cast_to_supertypes
            Cast the input datatypes to their supertype.
        pass_name_to_apply
            if set, then the `Series` passed to the function in the group_by operation
            will ensure the name is set. This is an extra heap allocation per group.
        changes_length
            For example a `unique` or a `slice`
        """
        from polars.plugins import register_plugin_function

        if args is None:
            args = [self]
        else:
            args = [self, *list(args)]

        return register_plugin_function(
            plugin_path=lib,
            function_name=symbol,
            args=args,
            kwargs=kwargs,
            is_elementwise=is_elementwise,
            changes_length=changes_length,
            returns_scalar=returns_scalar,
            cast_to_supertype=cast_to_supertypes,
            input_wildcard_expansion=input_wildcard_expansion,
            pass_name_to_apply=pass_name_to_apply,
        )

    def _row_encode(
        self,
        *,
        unordered: bool = False,
        descending: bool | None = None,
        nulls_last: bool | None = None,
    ) -> Expr:
        return F._row_encode(
            [self],
            unordered=unordered,
            descending=None if descending is None else [descending],
            nulls_last=None if nulls_last is None else [nulls_last],
        )

    def _row_decode(
        self,
        names: Sequence[str],
        dtypes: Sequence[pl.DataTypeExpr | PolarsDataType],
        *,
        unordered: bool = False,
        descending: Sequence[bool] | None = None,
        nulls_last: Sequence[bool] | None = None,
    ) -> Expr:
        dtypes_pyexprs = [
            parse_into_datatype_expr(dtype)._pydatatype_expr for dtype in dtypes
        ]

        if unordered:
            assert descending is None
            assert nulls_last is None

            result = self._pyexpr.row_decode_unordered(names, dtypes_pyexprs)
        else:
            result = self._pyexpr.row_decode_ordered(
                names, dtypes_pyexprs, descending, nulls_last
            )

        return wrap_expr(result)

    @classmethod
    def from_json(cls, value: str) -> Expr:
        """
        Read an expression from a JSON encoded string to construct an Expression.

        .. deprecated:: 0.20.11
            This method has been renamed to :meth:`deserialize`.
            Note that the new method operates on file-like inputs rather than strings.
            Enclose your input in `io.StringIO` to keep the same behavior.

        Parameters
        ----------
        value
            JSON encoded string value
        """
        issue_deprecation_warning(
            "`Expr.from_json` is deprecated. It has been renamed to `Expr.deserialize`."
            " Note that the new method operates on file-like inputs rather than strings."
            " Enclose your input in `io.StringIO` to keep the same behavior.",
            version="0.20.11",
        )
        return cls.deserialize(StringIO(value), format="json")

    @property
    def bin(self) -> ExprBinaryNameSpace:
        """
        Create an object namespace of all binary related methods.

        See the individual method pages for full details
        """
        return ExprBinaryNameSpace(self)

    @property
    def cat(self) -> ExprCatNameSpace:
        """
        Create an object namespace of all categorical related methods.

        See the individual method pages for full details

        Examples
        --------
        >>> df = pl.DataFrame({"values": ["a", "b"]}).select(
        ...     pl.col("values").cast(pl.Categorical)
        ... )
        >>> df.select(pl.col("values").cat.get_categories())
        shape: (2, 1)
        ┌────────┐
        │ values │
        │ ---    │
        │ str    │
        ╞════════╡
        │ a      │
        │ b      │
        └────────┘
        """
        return ExprCatNameSpace(self)

    @property
    def dt(self) -> ExprDateTimeNameSpace:
        """Create an object namespace of all datetime related methods."""
        return ExprDateTimeNameSpace(self)

    # Keep the `list` and `str` properties below at the end of the definition of Expr,
    # as to not confuse mypy with the type annotation `str` and `list`

    @property
    def list(self) -> ExprListNameSpace:
        """
        Create an object namespace of all list related methods.

        See the individual method pages for full details.
        """
        return ExprListNameSpace(self)

    @property
    def arr(self) -> ExprArrayNameSpace:
        """
        Create an object namespace of all array related methods.

        See the individual method pages for full details.
        """
        return ExprArrayNameSpace(self)

    @property
    def meta(self) -> ExprMetaNameSpace:
        """
        Create an object namespace of all meta related expression methods.

        This can be used to modify and traverse existing expressions.
        """
        return ExprMetaNameSpace(self)

    @property
    def name(self) -> ExprNameNameSpace:
        """
        Create an object namespace of all expressions that modify expression names.

        See the individual method pages for full details.
        """
        return ExprNameNameSpace(self)

    @property
    def str(self) -> ExprStringNameSpace:
        """
        Create an object namespace of all string related methods.

        See the individual method pages for full details.

        Examples
        --------
        >>> df = pl.DataFrame({"letters": ["a", "b"]})
        >>> df.select(pl.col("letters").str.to_uppercase())
        shape: (2, 1)
        ┌─────────┐
        │ letters │
        │ ---     │
        │ str     │
        ╞═════════╡
        │ A       │
        │ B       │
        └─────────┘
        """
        return ExprStringNameSpace(self)

    @property
    def struct(self) -> ExprStructNameSpace:
        """
        Create an object namespace of all struct related methods.

        See the individual method pages for full details.

        Examples
        --------
        >>> df = (
        ...     pl.DataFrame(
        ...         {
        ...             "int": [1, 2],
        ...             "str": ["a", "b"],
        ...             "bool": [True, None],
        ...             "list": [[1, 2], [3]],
        ...         }
        ...     )
        ...     .to_struct("my_struct")
        ...     .to_frame()
        ... )
        >>> df.select(pl.col("my_struct").struct.field("str"))
        shape: (2, 1)
        ┌─────┐
        │ str │
        │ --- │
        │ str │
        ╞═════╡
        │ a   │
        │ b   │
        └─────┘
        """
        return ExprStructNameSpace(self)

    def _skip_batch_predicate(self, schema: SchemaDict) -> Expr | None:
        result = self._pyexpr.skip_batch_predicate(schema)
        if result is None:
            return None
        return wrap_expr(result)


def _prepare_alpha(
    com: float | int | None = None,
    span: float | int | None = None,
    half_life: float | int | None = None,
    alpha: float | int | None = None,
) -> float:
    """Normalise EWM decay specification in terms of smoothing factor 'alpha'."""
    if sum((param is not None) for param in (com, span, half_life, alpha)) > 1:
        msg = (
            "parameters `com`, `span`, `half_life`, and `alpha` are mutually exclusive"
        )
        raise ValueError(msg)
    if com is not None:
        if com < 0.0:
            msg = f"require `com` >= 0 (found {com!r})"
            raise ValueError(msg)
        alpha = 1.0 / (1.0 + com)

    elif span is not None:
        if span < 1.0:
            msg = f"require `span` >= 1 (found {span!r})"
            raise ValueError(msg)
        alpha = 2.0 / (span + 1.0)

    elif half_life is not None:
        if half_life <= 0.0:
            msg = f"require `half_life` > 0 (found {half_life!r})"
            raise ValueError(msg)
        alpha = 1.0 - math.exp(-math.log(2.0) / half_life)

    elif alpha is None:
        msg = "one of `com`, `span`, `half_life`, or `alpha` must be set"
        raise ValueError(msg)

    elif not (0 < alpha <= 1):
        msg = f"require 0 < `alpha` <= 1 (found {alpha!r})"
        raise ValueError(msg)

    return alpha


def _prepare_rolling_by_window_args(window_size: timedelta | str) -> str:
    if isinstance(window_size, timedelta):
        window_size = parse_as_duration_string(window_size)
    return window_size
