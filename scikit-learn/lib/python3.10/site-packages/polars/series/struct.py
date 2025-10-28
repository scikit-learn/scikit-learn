from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from polars._utils.various import (
    BUILDING_SPHINX_DOCS,
    qualified_type_name,
    sphinx_accessor,
)
from polars._utils.wrap import wrap_df
from polars.schema import Schema
from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars import DataFrame, Series
    from polars._plr import PySeries
elif BUILDING_SPHINX_DOCS:
    # note: we assign this way to work around an autocomplete issue in ipython/jedi
    # (ref: https://github.com/davidhalter/jedi/issues/2057)
    current_module = sys.modules[__name__]
    current_module.property = sphinx_accessor


@expr_dispatch
class StructNameSpace:
    """Series.struct namespace."""

    _accessor = "struct"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def __getitem__(self, item: int | str) -> Series:
        if isinstance(item, int):
            return self.field(self.fields[item])
        elif isinstance(item, str):
            return self.field(item)
        else:
            msg = f"expected type 'int | str', got {qualified_type_name(item)!r}"
            raise TypeError(msg)

    def _ipython_key_completions_(self) -> list[str]:
        return self.fields

    @property
    def fields(self) -> list[str]:
        """
        Get the names of the fields.

        Examples
        --------
        >>> s = pl.Series([{"a": 1, "b": 2}, {"a": 3, "b": 4}])
        >>> s.struct.fields
        ['a', 'b']
        """
        if getattr(self, "_s", None) is None:
            return []
        return self._s.struct_fields()

    def field(self, name: str) -> Series:
        """
        Retrieve one of the fields of this `Struct` as a new Series.

        Parameters
        ----------
        name
            Name of the field.

        Examples
        --------
        >>> s = pl.Series([{"a": 1, "b": 2}, {"a": 3, "b": 4}])
        >>> s.struct.field("a")
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            3
        ]
        """

    def rename_fields(self, names: Sequence[str]) -> Series:
        """
        Rename the fields of the struct.

        Parameters
        ----------
        names
            New names in the order of the struct's fields.

        Examples
        --------
        >>> s = pl.Series([{"a": 1, "b": 2}, {"a": 3, "b": 4}])
        >>> s.struct.fields
        ['a', 'b']
        >>> s = s.struct.rename_fields(["c", "d"])
        >>> s.struct.fields
        ['c', 'd']
        """

    @property
    def schema(self) -> Schema:
        """
        Get the struct definition as a name/dtype schema dict.

        Examples
        --------
        >>> s = pl.Series([{"a": 1, "b": 2}, {"a": 3, "b": 4}])
        >>> s.struct.schema
        Schema({'a': Int64, 'b': Int64})
        """
        if getattr(self, "_s", None) is None:
            return Schema({})

        schema = self._s.dtype().to_schema()
        return Schema(schema, check_dtypes=False)

    def unnest(self) -> DataFrame:
        """
        Convert this struct Series to a DataFrame with a separate column for each field.

        Examples
        --------
        >>> s = pl.Series([{"a": 1, "b": 2}, {"a": 3, "b": 4}])
        >>> s.struct.unnest()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 2   │
        │ 3   ┆ 4   │
        └─────┴─────┘
        """
        return wrap_df(self._s.struct_unnest())

    def json_encode(self) -> Series:
        """
        Convert this struct to a string column with json values.

        Examples
        --------
        >>> s = pl.Series("a", [{"a": [1, 2], "b": [45]}, {"a": [9, 1, 3], "b": None}])
        >>> s.struct.json_encode()
        shape: (2,)
        Series: 'a' [str]
        [
            "{"a":[1,2],"b":[45]}"
            "{"a":[9,1,3],"b":null}"
        ]
        """
