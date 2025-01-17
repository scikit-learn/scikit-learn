from __future__ import annotations

from typing import TYPE_CHECKING

from polars.interchange.protocol import (
    Buffer,
    CopyNotAllowedError,
    DlpackDeviceType,
    DtypeKind,
)
from polars.interchange.utils import polars_dtype_to_dtype

if TYPE_CHECKING:
    from typing import NoReturn

    from polars import Series


class PolarsBuffer(Buffer):
    """
    A buffer object backed by a Polars Series consisting of a single chunk.

    Parameters
    ----------
    data
        The Polars Series backing the buffer object.
    allow_copy
        Allow data to be copied during operations on this column. If set to `False`,
        a RuntimeError will be raised if data would be copied.
    """

    def __init__(self, data: Series, *, allow_copy: bool = True) -> None:
        if data.n_chunks() > 1:
            if not allow_copy:
                msg = "non-contiguous buffer must be made contiguous"
                raise CopyNotAllowedError(msg)
            data = data.rechunk()

        self._data = data

    @property
    def bufsize(self) -> int:
        """Buffer size in bytes."""
        dtype = polars_dtype_to_dtype(self._data.dtype)

        if dtype[0] == DtypeKind.BOOL:
            _, offset, length = self._data._get_buffer_info()
            n_bits = offset + length
            n_bytes, rest = divmod(n_bits, 8)
            # Round up to the nearest byte
            if rest == 0:
                return n_bytes
            else:
                return n_bytes + 1

        return self._data.len() * (dtype[1] // 8)

    @property
    def ptr(self) -> int:
        """Pointer to start of the buffer as an integer."""
        pointer, _, _ = self._data._get_buffer_info()
        return pointer

    def __dlpack__(self) -> NoReturn:
        """Represent this structure as DLPack interface."""
        msg = "__dlpack__"
        raise NotImplementedError(msg)

    def __dlpack_device__(self) -> tuple[DlpackDeviceType, None]:
        """Device type and device ID for where the data in the buffer resides."""
        return (DlpackDeviceType.CPU, None)

    def __repr__(self) -> str:
        bufsize = self.bufsize
        ptr = self.ptr
        device = self.__dlpack_device__()[0].name
        return f"PolarsBuffer(bufsize={bufsize}, ptr={ptr}, device={device!r})"
