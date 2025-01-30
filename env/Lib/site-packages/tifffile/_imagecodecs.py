# tifffile/_imagecodecs.py

# Copyright (c) 2008-2024, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Fallback imagecodecs codecs.

This module provides alternative, pure Python and NumPy implementations of
some functions of the `imagecodecs`_ package. The functions may raise
`NotImplementedError`.

.. _imagecodecs: https://github.com/cgohlke/imagecodecs

"""

from __future__ import annotations

__all__ = [
    'bitorder_decode',
    'delta_decode',
    'delta_encode',
    'float24_decode',
    'lzma_decode',
    'lzma_encode',
    'packbits_decode',
    'packints_decode',
    'packints_encode',
    'zlib_decode',
    'zlib_encode',
]

from typing import TYPE_CHECKING, overload

import numpy

if TYPE_CHECKING:
    from typing import Any, Literal

    from numpy.typing import DTypeLike, NDArray

try:
    import lzma

    def lzma_encode(
        data: bytes | NDArray[Any],
        /,
        level: int | None = None,
        *,
        out: Any = None,
    ) -> bytes:
        """Compress LZMA."""
        if isinstance(data, numpy.ndarray):
            data = data.tobytes()
        return lzma.compress(data)

    def lzma_decode(data: bytes, /, *, out: Any = None) -> bytes:
        """Decompress LZMA."""
        return lzma.decompress(data)

except ImportError:
    # Python was built without lzma
    def lzma_encode(
        data: bytes | NDArray[Any],
        /,
        level: int | None = None,
        *,
        out: Any = None,
    ) -> bytes:
        """Raise ImportError."""
        import lzma  # noqa

        return b''

    def lzma_decode(data: bytes, /, *, out: Any = None) -> bytes:
        """Raise ImportError."""
        import lzma  # noqa

        return b''


try:
    import zlib

    def zlib_encode(
        data: bytes | NDArray[Any],
        /,
        level: int | None = None,
        *,
        out: Any = None,
    ) -> bytes:
        """Compress Zlib DEFLATE."""
        if isinstance(data, numpy.ndarray):
            data = data.tobytes()
        return zlib.compress(data, 6 if level is None else level)

    def zlib_decode(data: bytes, /, *, out: Any = None) -> bytes:
        """Decompress Zlib DEFLATE."""
        return zlib.decompress(data)

except ImportError:
    # Python was built without zlib

    def zlib_encode(
        data: bytes | NDArray[Any],
        /,
        level: int | None = None,
        *,
        out: Any = None,
    ) -> bytes:
        """Raise ImportError."""
        import zlib  # noqa

        return b''

    def zlib_decode(data: bytes, /, *, out: Any = None) -> bytes:
        """Raise ImportError."""
        import zlib  # noqa

        return b''


def packbits_decode(encoded: bytes, /, *, out: Any = None) -> bytes:
    r"""Decompress PackBits encoded byte string.

    >>> packbits_decode(b'\x80\x80')  # NOP
    b''
    >>> packbits_decode(b'\x02123')
    b'123'
    >>> packbits_decode(
    ...     b'\xfe\xaa\x02\x80\x00\x2a\xfd\xaa\x03\x80\x00\x2a\x22\xf7\xaa'
    ... )[:-5]
    b'\xaa\xaa\xaa\x80\x00*\xaa\xaa\xaa\xaa\x80\x00*"\xaa\xaa\xaa\xaa\xaa'

    """
    out = []
    out_extend = out.extend
    i = 0
    try:
        while True:
            n = ord(encoded[i : i + 1]) + 1
            i += 1
            if n > 129:
                # replicate
                out_extend(encoded[i : i + 1] * (258 - n))
                i += 1
            elif n < 129:
                # literal
                out_extend(encoded[i : i + n])
                i += n
    except TypeError:
        pass
    return bytes(out)


@overload
def delta_encode(
    data: bytes | bytearray,
    /,
    axis: int = -1,
    dist: int = 1,
    *,
    out: Any = None,
) -> bytes: ...


@overload
def delta_encode(
    data: NDArray[Any], /, axis: int = -1, dist: int = 1, *, out: Any = None
) -> NDArray[Any]: ...


def delta_encode(
    data: bytes | bytearray | NDArray[Any],
    /,
    axis: int = -1,
    dist: int = 1,
    *,
    out: Any = None,
) -> bytes | NDArray[Any]:
    """Encode Delta."""
    if dist != 1:
        raise NotImplementedError(
            f"delta_encode with {dist=} requires the 'imagecodecs' package"
        )
    if isinstance(data, (bytes, bytearray)):
        data = numpy.frombuffer(data, dtype=numpy.uint8)
        diff = numpy.diff(data, axis=0)
        return numpy.insert(diff, 0, data[0]).tobytes()

    dtype = data.dtype
    if dtype.kind == 'f':
        data = data.view(f'{dtype.byteorder}u{dtype.itemsize}')
    diff = numpy.diff(data, axis=axis)
    key: list[int | slice] = [slice(None)] * data.ndim
    key[axis] = 0
    diff = numpy.insert(diff, 0, data[tuple(key)], axis=axis)
    if not data.dtype.isnative:
        diff = diff.byteswap(True)
        diff = diff.view(diff.dtype.newbyteorder())
    if dtype.kind == 'f':
        return diff.view(dtype)
    return diff


@overload
def delta_decode(
    data: bytes | bytearray, /, axis: int, dist: int, *, out: Any
) -> bytes: ...


@overload
def delta_decode(
    data: NDArray[Any], /, axis: int, dist: int, *, out: Any
) -> NDArray[Any]: ...


def delta_decode(
    data: bytes | bytearray | NDArray[Any],
    /,
    axis: int = -1,
    dist: int = 1,
    *,
    out: Any = None,
) -> bytes | NDArray[Any]:
    """Decode Delta."""
    if dist != 1:
        raise NotImplementedError(
            f"delta_decode with {dist=} requires the 'imagecodecs' package"
        )
    if out is not None and not out.flags.writeable:
        out = None
    if isinstance(data, (bytes, bytearray)):
        data = numpy.frombuffer(data, dtype=numpy.uint8)
        return numpy.cumsum(  # type: ignore[no-any-return]
            data, axis=0, dtype=numpy.uint8, out=out
        ).tobytes()
    if data.dtype.kind == 'f':
        if not data.dtype.isnative:
            raise NotImplementedError(
                f'delta_decode with {data.dtype!r} '
                "requires the 'imagecodecs' package"
            )
        view = data.view(f'{data.dtype.byteorder}u{data.dtype.itemsize}')
        view = numpy.cumsum(view, axis=axis, dtype=view.dtype)
        return view.view(data.dtype)
    return numpy.cumsum(  # type: ignore[no-any-return]
        data, axis=axis, dtype=data.dtype, out=out
    )


@overload
def bitorder_decode(
    data: bytes | bytearray, /, *, out: Any = None, _bitorder: list[Any] = []
) -> bytes: ...


@overload
def bitorder_decode(
    data: NDArray[Any], /, *, out: Any = None, _bitorder: list[Any] = []
) -> NDArray[Any]: ...


def bitorder_decode(
    data: bytes | bytearray | NDArray[Any],
    /,
    *,
    out: Any = None,
    _bitorder: list[Any] = [],
) -> bytes | NDArray[Any]:
    r"""Reverse bits in each byte of bytes or numpy array.

    Decode data where pixels with lower column values are stored in the
    lower-order bits of the bytes (TIFF FillOrder is LSB2MSB).

    Parameters:
        data:
            Data to bit-reversed. If bytes type, a new bit-reversed
            bytes is returned. NumPy arrays are bit-reversed in-place.

    Examples:
        >>> bitorder_decode(b'\x01\x64')
        b'\x80&'
        >>> data = numpy.array([1, 666], dtype='uint16')
        >>> bitorder_decode(data)
        >>> data
        array([  128, 16473], dtype=uint16)

    """
    if not _bitorder:
        _bitorder.append(
            b'\x00\x80@\xc0 \xa0`\xe0\x10\x90P\xd00\xb0p\xf0\x08\x88H'
            b'\xc8(\xa8h\xe8\x18\x98X\xd88\xb8x\xf8\x04\x84D\xc4$\xa4d'
            b'\xe4\x14\x94T\xd44\xb4t\xf4\x0c\x8cL\xcc,\xacl\xec\x1c\x9c'
            b'\\\xdc<\xbc|\xfc\x02\x82B\xc2"\xa2b\xe2\x12\x92R\xd22'
            b'\xb2r\xf2\n\x8aJ\xca*\xaaj\xea\x1a\x9aZ\xda:\xbaz\xfa'
            b'\x06\x86F\xc6&\xa6f\xe6\x16\x96V\xd66\xb6v\xf6\x0e\x8eN'
            b'\xce.\xaen\xee\x1e\x9e^\xde>\xbe~\xfe\x01\x81A\xc1!\xa1a'
            b'\xe1\x11\x91Q\xd11\xb1q\xf1\t\x89I\xc9)\xa9i\xe9\x19'
            b'\x99Y\xd99\xb9y\xf9\x05\x85E\xc5%\xa5e\xe5\x15\x95U\xd55'
            b'\xb5u\xf5\r\x8dM\xcd-\xadm\xed\x1d\x9d]\xdd=\xbd}\xfd'
            b'\x03\x83C\xc3#\xa3c\xe3\x13\x93S\xd33\xb3s\xf3\x0b\x8bK'
            b'\xcb+\xabk\xeb\x1b\x9b[\xdb;\xbb{\xfb\x07\x87G\xc7\'\xa7g'
            b'\xe7\x17\x97W\xd77\xb7w\xf7\x0f\x8fO\xcf/\xafo\xef\x1f\x9f_'
            b'\xdf?\xbf\x7f\xff'
        )
        _bitorder.append(numpy.frombuffer(_bitorder[0], dtype=numpy.uint8))
    if isinstance(data, (bytes, bytearray)):
        return data.translate(_bitorder[0])
    try:
        view = data.view('uint8')
        numpy.take(_bitorder[1], view, out=view)
        return data
    except ValueError as exc:
        raise NotImplementedError(
            "bitorder_decode of slices requires the 'imagecodecs' package"
        ) from exc
    return None  # type: ignore[unreachable]


def packints_decode(
    data: bytes,
    /,
    dtype: DTypeLike,
    bitspersample: int,
    runlen: int = 0,
    *,
    out: Any = None,
) -> NDArray[Any]:
    """Decompress bytes to array of integers.

    This implementation only handles itemsizes 1, 8, 16, 32, and 64 bits.
    Install the Imagecodecs package for decoding other integer sizes.

    Parameters:
        data:
            Data to decompress.
        dtype:
            Numpy boolean or integer type.
        bitspersample:
            Number of bits per integer.
        runlen:
            Number of consecutive integers after which to start at next byte.

    Examples:
        >>> packints_decode(b'a', 'B', 1)
        array([0, 1, 1, 0, 0, 0, 0, 1], dtype=uint8)

    """
    if bitspersample == 1:  # bitarray
        data_array = numpy.frombuffer(data, '|B')
        data_array = numpy.unpackbits(data_array)
        if runlen % 8:
            data_array = data_array.reshape(-1, runlen + (8 - runlen % 8))
            data_array = data_array[:, :runlen].reshape(-1)
        return data_array.astype(dtype)
    if bitspersample in (8, 16, 32, 64):
        return numpy.frombuffer(data, dtype)
    raise NotImplementedError(
        f'packints_decode of {bitspersample}-bit integers '
        "requires the 'imagecodecs' package"
    )


def packints_encode(
    data: NDArray[Any],
    /,
    bitspersample: int,
    axis: int = -1,
    *,
    out: Any = None,
) -> bytes:
    """Tightly pack integers."""
    raise NotImplementedError(
        "packints_encode requires the 'imagecodecs' package"
    )


def float24_decode(
    data: bytes, /, byteorder: Literal['>', '<']
) -> NDArray[Any]:
    """Return float32 array from float24."""
    raise NotImplementedError(
        "float24_decode requires the 'imagecodecs' package"
    )
