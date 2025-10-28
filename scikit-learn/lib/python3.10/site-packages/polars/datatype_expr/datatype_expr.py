from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import polars._reexport as pl
from polars._utils.various import BUILDING_SPHINX_DOCS, sphinx_accessor
from polars.datatype_expr.array import DataTypeExprArrNameSpace
from polars.datatype_expr.list import DataTypeExprListNameSpace
from polars.datatype_expr.struct import DataTypeExprStructNameSpace

if TYPE_CHECKING:
    import contextlib
    from typing import ClassVar

    from polars import DataType
    from polars._typing import PolarsDataType, SchemaDict

    with contextlib.suppress(ImportError):  # Module not available when building docs
        from polars._plr import PyDataTypeExpr
elif BUILDING_SPHINX_DOCS:
    import sys

    # note: we assign this way to work around an autocomplete issue in ipython/jedi
    # (ref: https://github.com/davidhalter/jedi/issues/2057)
    current_module = sys.modules[__name__]
    current_module.property = sphinx_accessor


class DataTypeExpr:
    """
    A lazily instantiated :class:`DataType` that can be used in an :class:`Expr`.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    This expression is made to represent a :class:`DataType` that can be used to
    reference a datatype in a lazy context.

    Examples
    --------
    >>> lf = pl.LazyFrame({"a": [1, 2, 3]})
    >>> lf.with_columns(
    ...     pl.col.a.map_batches(lambda x: x * 2, return_dtype=pl.dtype_of("a"))
    ... ).collect()
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 2   │
    │ 4   │
    │ 6   │
    └─────┘
    """

    # NOTE: This `= None` is needed to generate the docs with sphinx_accessor.
    _pydatatype_expr: PyDataTypeExpr = None  # type: ignore[assignment]
    _accessors: ClassVar[set[str]] = {
        "arr",
        "enum",
        "list",
        "struct",
    }

    def __eq__(self, value: PolarsDataType | DataTypeExpr) -> pl.Expr:  # type: ignore[override]
        cmp_with: DataTypeExpr
        if isinstance(value, pl.DataType):
            cmp_with = value.to_dtype_expr()
        elif isinstance(value, pl.DataTypeClass):
            cmp_with = value.to_dtype_expr()
        elif isinstance(value, DataTypeExpr):
            cmp_with = value
        else:
            msg = f"cannot compare {self!r} to {value!r}"
            raise TypeError(msg) from None

        return pl.Expr._from_pyexpr(
            self._pydatatype_expr.equals(cmp_with._pydatatype_expr)
        )

    def __ne__(self, value: PolarsDataType | DataTypeExpr) -> pl.Expr:  # type: ignore[override]
        return (self == value).not_()

    @classmethod
    def _from_pydatatype_expr(cls, pydatatype_expr: PyDataTypeExpr) -> DataTypeExpr:
        slf = cls()
        slf._pydatatype_expr = pydatatype_expr
        return slf

    def inner_dtype(self) -> DataTypeExpr:
        """Get the inner DataType of a List or Array."""
        return DataTypeExpr._from_pydatatype_expr(self._pydatatype_expr.inner_dtype())

    def display(self) -> pl.Expr:
        """
        Get a formatted version of the output DataType.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": ["X", "Y", "Z"],
        ...         "c": [1.3, 3.7, 4.2],
        ...     }
        ... )
        >>> df.select(
        ...     a=pl.dtype_of("a").display(),
        ...     b=pl.dtype_of("b").display(),
        ...     c=pl.dtype_of("c").display(),
        ... ).transpose(include_header=True, column_names=["dtype"])
        shape: (3, 2)
        ┌────────┬───────┐
        │ column ┆ dtype │
        │ ---    ┆ ---   │
        │ str    ┆ str   │
        ╞════════╪═══════╡
        │ a      ┆ i64   │
        │ b      ┆ str   │
        │ c      ┆ f64   │
        └────────┴───────┘
        """
        return pl.Expr._from_pyexpr(self._pydatatype_expr.display())

    def matches(self, selector: pl.Selector) -> pl.Expr:
        """
        Get whether the output DataType is matches a certain selector.

        Examples
        --------
        >>> import polars.selectors as cs
        >>> pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...     }
        ... ).select(
        ...     a_is_string=pl.dtype_of("a").matches(cs.string()),
        ...     a_is_integer=pl.dtype_of("a").matches(cs.integer()),
        ... )
        shape: (1, 2)
        ┌─────────────┬──────────────┐
        │ a_is_string ┆ a_is_integer │
        │ ---         ┆ ---          │
        │ bool        ┆ bool         │
        ╞═════════════╪══════════════╡
        │ false       ┆ true         │
        └─────────────┴──────────────┘
        """
        return pl.Expr._from_pyexpr(self._pydatatype_expr.matches(selector._pyselector))

    def wrap_in_list(self) -> DataTypeExpr:
        """
        Get the DataType wrapped in a list.

        Examples
        --------
        >>> pl.Int32.to_dtype_expr().wrap_in_list().collect_dtype({})
        List(Int32)

        """
        return DataTypeExpr._from_pydatatype_expr(self._pydatatype_expr.wrap_in_list())

    def wrap_in_array(self, *, width: int) -> DataTypeExpr:
        """
        Get the DataType wrapped in an array.

        Examples
        --------
        >>> pl.Int32.to_dtype_expr().wrap_in_array(width=5).collect_dtype({})
        Array(Int32, shape=(5,))
        """
        return DataTypeExpr._from_pydatatype_expr(
            self._pydatatype_expr.wrap_in_array(width)
        )

    def to_unsigned_integer(self) -> pl.DataTypeExpr:
        """
        Get the unsigned integer version of the same bitsize.

        Examples
        --------
        >>> int32 = pl.Int32.to_dtype_expr()
        >>> int32.to_unsigned_integer().collect_dtype({})
        UInt32
        """
        return pl.DataTypeExpr._from_pydatatype_expr(
            self._pydatatype_expr.to_unsigned_integer()
        )

    def to_signed_integer(self) -> pl.DataTypeExpr:
        """
        Get the signed integer version of the same bitsize.

        Examples
        --------
        >>> uint32 = pl.UInt32.to_dtype_expr()
        >>> uint32.to_signed_integer().collect_dtype({})
        Int32
        """
        return pl.DataTypeExpr._from_pydatatype_expr(
            self._pydatatype_expr.to_signed_integer()
        )

    def default_value(
        self,
        n: int = 1,
        *,
        numeric_to_one: bool = False,
        num_list_values: int = 0,
    ) -> pl.Expr:
        """
        Get a default value of a specific type.

        - Integers and floats are their zero value as default, unless otherwise
          specified
        - Temporals are a physical zero as default
        - `pl.Decimal` is zero as default
        - `pl.String` and `pl.Binary` are an empty string
        - `pl.List` is an empty list, unless otherwise specified
        - `pl.Array` is the inner default value repeated over the shape
        - `pl.Struct` is the inner default value for all fields
        - `pl.Enum` is the first category if it exists
        - `pl.Null`, `pl.Object` and `pl.Categorical` are `null`.

        Parameters
        ----------
        n
            Number of types you want the value
        numeric_to_one
            Use `1` instead of `0` as the default value for numeric types
        num_list_values
            The amount of values a list contains

        Examples
        --------
        >>> uint32 = pl.UInt32.to_dtype_expr()
        >>> pl.select(default=uint32.default_value())
        shape: (1, 1)
        ┌─────────┐
        │ default │
        │ ---     │
        │ u32     │
        ╞═════════╡
        │ 0       │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(
            self._pydatatype_expr.default_value(
                n=n, numeric_to_one=numeric_to_one, num_list_values=num_list_values
            )
        )

    @property
    def list(self) -> DataTypeExprListNameSpace:
        """Create an object namespace of all list related methods."""
        return DataTypeExprListNameSpace(self)

    @property
    def arr(self) -> DataTypeExprArrNameSpace:
        """Create an object namespace of all array related methods."""
        return DataTypeExprArrNameSpace(self)

    @property
    def struct(self) -> DataTypeExprStructNameSpace:
        """Create an object namespace of all struct related methods."""
        return DataTypeExprStructNameSpace(self)

    def collect_dtype(
        self, context: SchemaDict | pl.Schema | pl.DataFrame | pl.LazyFrame
    ) -> DataType:
        """
        Materialize the :class:`DataTypeExpr` in a specific context.

        This is a useful function when debugging datatype expressions.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...     }
        ... )
        >>> pl.dtype_of("a").collect_dtype(lf)
        Int64
        >>> pl.dtype_of("a").collect_dtype({"a": pl.String})
        String
        """
        schema: pl.Schema
        if isinstance(context, pl.Schema):
            schema = context
        elif isinstance(context, Mapping):
            schema = pl.Schema(context)
        elif isinstance(context, pl.DataFrame):
            schema = context.schema
        elif isinstance(context, pl.LazyFrame):
            schema = context.collect_schema()
        else:
            msg = f"DataTypeExpr.collect_dtype did not expect {context!r}"
            raise TypeError(msg)

        return self._pydatatype_expr.collect_dtype(schema)
