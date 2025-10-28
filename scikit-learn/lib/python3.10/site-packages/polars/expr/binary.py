from __future__ import annotations

from typing import TYPE_CHECKING

from polars._utils.parse import parse_into_expression
from polars._utils.various import scale_bytes
from polars._utils.wrap import wrap_expr
from polars.datatypes import parse_into_datatype_expr

if TYPE_CHECKING:
    from polars import DataTypeExpr, Expr
    from polars._typing import (
        Endianness,
        IntoExpr,
        PolarsDataType,
        SizeUnit,
        TransferEncoding,
    )


class ExprBinaryNameSpace:
    """Namespace for bin related expressions."""

    _accessor = "bin"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def contains(self, literal: IntoExpr) -> Expr:
        r"""
        Check if binaries in Series contain a binary substring.

        Parameters
        ----------
        literal
            The binary substring to look for

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        See Also
        --------
        starts_with : Check if the binary substring exists at the start
        ends_with : Check if the binary substring exists at the end

        Examples
        --------
        >>> colors = pl.DataFrame(
        ...     {
        ...         "name": ["black", "yellow", "blue"],
        ...         "code": [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"],
        ...         "lit": [b"\x00", b"\xff\x00", b"\xff\xff"],
        ...     }
        ... )
        >>> colors.select(
        ...     "name",
        ...     pl.col("code").bin.contains(b"\xff").alias("contains_with_lit"),
        ...     pl.col("code").bin.contains(pl.col("lit")).alias("contains_with_expr"),
        ... )
        shape: (3, 3)
        ┌────────┬───────────────────┬────────────────────┐
        │ name   ┆ contains_with_lit ┆ contains_with_expr │
        │ ---    ┆ ---               ┆ ---                │
        │ str    ┆ bool              ┆ bool               │
        ╞════════╪═══════════════════╪════════════════════╡
        │ black  ┆ false             ┆ true               │
        │ yellow ┆ true              ┆ true               │
        │ blue   ┆ true              ┆ false              │
        └────────┴───────────────────┴────────────────────┘
        """
        literal_pyexpr = parse_into_expression(literal, str_as_lit=True)
        return wrap_expr(self._pyexpr.bin_contains(literal_pyexpr))

    def ends_with(self, suffix: IntoExpr) -> Expr:
        r"""
        Check if string values end with a binary substring.

        Parameters
        ----------
        suffix
            Suffix substring.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        See Also
        --------
        starts_with : Check if the binary substring exists at the start
        contains : Check if the binary substring exists anywhere

        Examples
        --------
        >>> colors = pl.DataFrame(
        ...     {
        ...         "name": ["black", "yellow", "blue"],
        ...         "code": [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"],
        ...         "suffix": [b"\x00", b"\xff\x00", b"\x00\x00"],
        ...     }
        ... )
        >>> colors.select(
        ...     "name",
        ...     pl.col("code").bin.ends_with(b"\xff").alias("ends_with_lit"),
        ...     pl.col("code").bin.ends_with(pl.col("suffix")).alias("ends_with_expr"),
        ... )
        shape: (3, 3)
        ┌────────┬───────────────┬────────────────┐
        │ name   ┆ ends_with_lit ┆ ends_with_expr │
        │ ---    ┆ ---           ┆ ---            │
        │ str    ┆ bool          ┆ bool           │
        ╞════════╪═══════════════╪════════════════╡
        │ black  ┆ false         ┆ true           │
        │ yellow ┆ false         ┆ true           │
        │ blue   ┆ true          ┆ false          │
        └────────┴───────────────┴────────────────┘
        """
        suffix_pyexpr = parse_into_expression(suffix, str_as_lit=True)
        return wrap_expr(self._pyexpr.bin_ends_with(suffix_pyexpr))

    def starts_with(self, prefix: IntoExpr) -> Expr:
        r"""
        Check if values start with a binary substring.

        Parameters
        ----------
        prefix
            Prefix substring.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        See Also
        --------
        ends_with : Check if the binary substring exists at the end
        contains : Check if the binary substring exists anywhere

        Examples
        --------
        >>> colors = pl.DataFrame(
        ...     {
        ...         "name": ["black", "yellow", "blue"],
        ...         "code": [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"],
        ...         "prefix": [b"\x00", b"\xff\x00", b"\x00\x00"],
        ...     }
        ... )
        >>> colors.select(
        ...     "name",
        ...     pl.col("code").bin.starts_with(b"\xff").alias("starts_with_lit"),
        ...     pl.col("code")
        ...     .bin.starts_with(pl.col("prefix"))
        ...     .alias("starts_with_expr"),
        ... )
        shape: (3, 3)
        ┌────────┬─────────────────┬──────────────────┐
        │ name   ┆ starts_with_lit ┆ starts_with_expr │
        │ ---    ┆ ---             ┆ ---              │
        │ str    ┆ bool            ┆ bool             │
        ╞════════╪═════════════════╪══════════════════╡
        │ black  ┆ false           ┆ true             │
        │ yellow ┆ true            ┆ false            │
        │ blue   ┆ false           ┆ true             │
        └────────┴─────────────────┴──────────────────┘
        """
        prefix_pyexpr = parse_into_expression(prefix, str_as_lit=True)
        return wrap_expr(self._pyexpr.bin_starts_with(prefix_pyexpr))

    def decode(self, encoding: TransferEncoding, *, strict: bool = True) -> Expr:
        r"""
        Decode values using the provided encoding.

        Parameters
        ----------
        encoding : {'hex', 'base64'}
            The encoding to use.
        strict
            Raise an error if the underlying value cannot be decoded,
            otherwise mask out with a null value.

        Returns
        -------
        Expr
            Expression of data type :class:`Binary`.

        Examples
        --------
        >>> colors = pl.DataFrame(
        ...     {
        ...         "name": ["black", "yellow", "blue"],
        ...         "encoded": [b"000000", b"ffff00", b"0000ff"],
        ...     }
        ... )
        >>> colors.with_columns(
        ...     pl.col("encoded").bin.decode("hex").alias("code"),
        ... )
        shape: (3, 3)
        ┌────────┬───────────┬─────────────────┐
        │ name   ┆ encoded   ┆ code            │
        │ ---    ┆ ---       ┆ ---             │
        │ str    ┆ binary    ┆ binary          │
        ╞════════╪═══════════╪═════════════════╡
        │ black  ┆ b"000000" ┆ b"\x00\x00\x00" │
        │ yellow ┆ b"ffff00" ┆ b"\xff\xff\x00" │
        │ blue   ┆ b"0000ff" ┆ b"\x00\x00\xff" │
        └────────┴───────────┴─────────────────┘
        """
        if encoding == "hex":
            return wrap_expr(self._pyexpr.bin_hex_decode(strict))
        elif encoding == "base64":
            return wrap_expr(self._pyexpr.bin_base64_decode(strict))
        else:
            msg = f"`encoding` must be one of {{'hex', 'base64'}}, got {encoding!r}"
            raise ValueError(msg)

    def encode(self, encoding: TransferEncoding) -> Expr:
        r"""
        Encode a value using the provided encoding.

        Parameters
        ----------
        encoding : {'hex', 'base64'}
            The encoding to use.

        Returns
        -------
        Expr
            Expression of data type :class:`Binary`.

        Examples
        --------
        >>> colors = pl.DataFrame(
        ...     {
        ...         "color": ["black", "yellow", "blue"],
        ...         "code": [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"],
        ...     }
        ... )
        >>> colors.with_columns(
        ...     pl.col("code").bin.encode("hex").alias("encoded"),
        ... )
        shape: (3, 3)
        ┌────────┬─────────────────┬─────────┐
        │ color  ┆ code            ┆ encoded │
        │ ---    ┆ ---             ┆ ---     │
        │ str    ┆ binary          ┆ str     │
        ╞════════╪═════════════════╪═════════╡
        │ black  ┆ b"\x00\x00\x00" ┆ 000000  │
        │ yellow ┆ b"\xff\xff\x00" ┆ ffff00  │
        │ blue   ┆ b"\x00\x00\xff" ┆ 0000ff  │
        └────────┴─────────────────┴─────────┘
        """
        if encoding == "hex":
            return wrap_expr(self._pyexpr.bin_hex_encode())
        elif encoding == "base64":
            return wrap_expr(self._pyexpr.bin_base64_encode())
        else:
            msg = f"`encoding` must be one of {{'hex', 'base64'}}, got {encoding!r}"
            raise ValueError(msg)

    def size(self, unit: SizeUnit = "b") -> Expr:
        r"""
        Get the size of binary values in the given unit.

        Parameters
        ----------
        unit : {'b', 'kb', 'mb', 'gb', 'tb'}
            Scale the returned size to the given unit.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32` or `Float64`.

        Examples
        --------
        >>> from os import urandom
        >>> df = pl.DataFrame({"data": [urandom(n) for n in (512, 256, 1024)]})
        >>> df.with_columns(  # doctest: +IGNORE_RESULT
        ...     n_bytes=pl.col("data").bin.size(),
        ...     n_kilobytes=pl.col("data").bin.size("kb"),
        ... )
        shape: (4, 3)
        ┌─────────────────────────────────┬─────────┬─────────────┐
        │ data                            ┆ n_bytes ┆ n_kilobytes │
        │ ---                             ┆ ---     ┆ ---         │
        │ binary                          ┆ u32     ┆ f64         │
        ╞═════════════════════════════════╪═════════╪═════════════╡
        │ b"y?~B\x83\xf4V\x07\xd3\xfb\xb… ┆ 512     ┆ 0.5         │
        │ b"\xee$4@f\xc14\x07\x8e\x88\x1… ┆ 256     ┆ 0.25        │
        │ b"\x80\xbd\xb9nEq;2\x99$\xf9\x… ┆ 1024    ┆ 1.0         │
        └─────────────────────────────────┴─────────┴─────────────┘
        """
        sz = wrap_expr(self._pyexpr.bin_size_bytes())
        sz = scale_bytes(sz, unit)
        return sz

    def reinterpret(
        self, *, dtype: PolarsDataType | DataTypeExpr, endianness: Endianness = "little"
    ) -> Expr:
        r"""
        Interpret bytes as another type.

        Supported types are numerical or temporal dtypes, or an ``Array`` of
        these dtypes.

        Parameters
        ----------
        dtype : PolarsDataType
            Which type to interpret binary column into.
        endianness : {"big", "little"}, optional
            Which endianness to use when interpreting bytes, by default "little".

        Returns
        -------
        Expr
            Expression of data type `dtype`.
            Note that rows of the binary array where the length does not match
            the size in bytes of the output array (number of items * byte size
            of item) will become NULL.

        Examples
        --------
        >>> df = pl.DataFrame({"data": [b"\x05\x00\x00\x00", b"\x10\x00\x01\x00"]})
        >>> df.with_columns(  # doctest: +IGNORE_RESULT
        ...     bin2int=pl.col("data").bin.reinterpret(
        ...         dtype=pl.Int32, endianness="little"
        ...     ),
        ... )
        shape: (2, 2)
        ┌─────────────────────┬─────────┐
        │ data                ┆ bin2int │
        │ ---                 ┆ ---     │
        │ binary              ┆ i32     │
        ╞═════════════════════╪═════════╡
        │ b"\x05\x00\x00\x00" ┆ 5       │
        │ b"\x10\x00\x01\x00" ┆ 65552   │
        └─────────────────────┴─────────┘
        """
        dtype = parse_into_datatype_expr(dtype)

        return wrap_expr(
            self._pyexpr.bin_reinterpret(dtype._pydatatype_expr, endianness)
        )
