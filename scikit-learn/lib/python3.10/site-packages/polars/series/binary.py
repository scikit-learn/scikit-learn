from __future__ import annotations

from typing import TYPE_CHECKING

from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    from polars import Series
    from polars._plr import PySeries
    from polars._typing import (
        Endianness,
        IntoExpr,
        PolarsDataType,
        SizeUnit,
        TransferEncoding,
    )


@expr_dispatch
class BinaryNameSpace:
    """Series.bin namespace."""

    _accessor = "bin"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def contains(self, literal: IntoExpr) -> Series:
        r"""
        Check if binaries in Series contain a binary substring.

        Parameters
        ----------
        literal
            The binary substring to look for

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("colors", [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"])
        >>> s.bin.contains(b"\xff")
        shape: (3,)
        Series: 'colors' [bool]
        [
            false
            true
            true
        ]
        """

    def ends_with(self, suffix: IntoExpr) -> Series:
        r"""
        Check if string values end with a binary substring.

        Parameters
        ----------
        suffix
            Suffix substring.

        Examples
        --------
        >>> s = pl.Series("colors", [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"])
        >>> s.bin.ends_with(b"\x00")
        shape: (3,)
        Series: 'colors' [bool]
        [
            true
            true
            false
        ]
        """

    def starts_with(self, prefix: IntoExpr) -> Series:
        r"""
        Check if values start with a binary substring.

        Parameters
        ----------
        prefix
            Prefix substring.

        Examples
        --------
        >>> s = pl.Series("colors", [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"])
        >>> s.bin.starts_with(b"\x00")
        shape: (3,)
        Series: 'colors' [bool]
        [
            true
            false
            true
        ]
        """

    def decode(self, encoding: TransferEncoding, *, strict: bool = True) -> Series:
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
        Series
            Series of data type :class:`String`.

        Examples
        --------
        Decode values using hexadecimal encoding.

        >>> s = pl.Series("colors", [b"000000", b"ffff00", b"0000ff"])
        >>> s.bin.decode("hex")
        shape: (3,)
        Series: 'colors' [binary]
        [
            b"\x00\x00\x00"
            b"\xff\xff\x00"
            b"\x00\x00\xff"
        ]

        Decode values using Base64 encoding.

        >>> s = pl.Series("colors", [b"AAAA", b"//8A", b"AAD/"])
        >>> s.bin.decode("base64")
        shape: (3,)
        Series: 'colors' [binary]
        [
            b"\x00\x00\x00"
            b"\xff\xff\x00"
            b"\x00\x00\xff"
        ]

        Set `strict=False` to set invalid values to null instead of raising an error.

        >>> s = pl.Series("colors", [b"000000", b"ffff00", b"invalid_value"])
        >>> s.bin.decode("hex", strict=False)
        shape: (3,)
        Series: 'colors' [binary]
        [
            b"\x00\x00\x00"
            b"\xff\xff\x00"
            null
        ]
        """

    def encode(self, encoding: TransferEncoding) -> Series:
        r"""
        Encode values using the provided encoding.

        Parameters
        ----------
        encoding : {'hex', 'base64'}
            The encoding to use.

        Returns
        -------
        Series
            Series of data type :class:`String`.

        Examples
        --------
        Encode values using hexadecimal encoding.

        >>> s = pl.Series("colors", [b"\x00\x00\x00", b"\xff\xff\x00", b"\x00\x00\xff"])
        >>> s.bin.encode("hex")
        shape: (3,)
        Series: 'colors' [str]
        [
            "000000"
            "ffff00"
            "0000ff"
        ]

        Encode values using Base64 encoding.

        >>> s.bin.encode("base64")
        shape: (3,)
        Series: 'colors' [str]
        [
            "AAAA"
            "//8A"
            "AAD/"
        ]
        """

    def size(self, unit: SizeUnit = "b") -> Series:
        r"""
        Get the size of the binary values in a Series in the given unit.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`.

        Examples
        --------
        >>> from os import urandom
        >>> s = pl.Series("data", [urandom(n) for n in (512, 256, 2560, 1024)])
        >>> s.bin.size("kb")
        shape: (4,)
        Series: 'data' [f64]
        [
            0.5
            0.25
            2.5
            1.0
        ]
        """

    def reinterpret(
        self, *, dtype: PolarsDataType, endianness: Endianness = "little"
    ) -> Series:
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
        Series
            Series of data type `dtype`.
            Note that rows of the binary array where the length does not match
            the size in bytes of the output array (number of items * byte size
            of item) will become NULL.

        Examples
        --------
        >>> s = pl.Series("data", [b"\x05\x00\x00\x00", b"\x10\x00\x01\x00"])
        >>> s.bin.reinterpret(dtype=pl.Int32, endianness="little")
        shape: (2,)
        Series: 'data' [i32]
        [
            5
            65552
        ]

        """
