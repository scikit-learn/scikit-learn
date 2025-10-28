from __future__ import annotations

from typing import TYPE_CHECKING

from polars.datatypes import Boolean, Categorical, Enum, String
from polars.interchange.buffer import PolarsBuffer
from polars.interchange.protocol import (
    Column,
    ColumnNullType,
    CopyNotAllowedError,
    DtypeKind,
    Endianness,
)
from polars.interchange.utils import polars_dtype_to_dtype

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Any

    from polars import Series
    from polars.interchange.protocol import CategoricalDescription, ColumnBuffers, Dtype


class PolarsColumn(Column):
    """
    A column object backed by a Polars Series.

    Parameters
    ----------
    column
        The Polars Series backing the column object.
    allow_copy
        Allow data to be copied during operations on this column. If set to `False`,
        a RuntimeError will be raised if data would be copied.
    """

    def __init__(self, column: Series, *, allow_copy: bool = True) -> None:
        self._col = column
        self._allow_copy = allow_copy

    def size(self) -> int:
        """Size of the column in elements."""
        return self._col.len()

    @property
    def offset(self) -> int:
        """Offset of the first element with respect to the start of the underlying buffer."""  # noqa: W505
        if self._col.dtype == Boolean:
            return self._col._get_buffer_info()[1]
        else:
            return 0

    @property
    def dtype(self) -> Dtype:
        """Data type of the column."""
        pl_dtype = self._col.dtype
        return polars_dtype_to_dtype(pl_dtype)

    @property
    def describe_categorical(self) -> CategoricalDescription:
        """
        Description of the categorical data type of the column.

        Raises
        ------
        TypeError
            If the data type of the column is not categorical.
        """
        dtype = self._col.dtype
        if dtype == Categorical:
            categories = self._col.cat.get_categories()
            is_ordered = False
        elif dtype == Enum:
            categories = dtype.categories  # type: ignore[attr-defined]
            is_ordered = True
        else:
            msg = "`describe_categorical` only works on categorical columns"
            raise TypeError(msg)

        return {
            "is_ordered": is_ordered,
            "is_dictionary": True,
            "categories": PolarsColumn(categories, allow_copy=self._allow_copy),
        }

    @property
    def describe_null(self) -> tuple[ColumnNullType, int | None]:
        """Description of the null representation the column uses."""
        if self.null_count == 0:
            return ColumnNullType.NON_NULLABLE, None
        else:
            return ColumnNullType.USE_BITMASK, 0

    @property
    def null_count(self) -> int:
        """The number of null elements."""
        return self._col.null_count()

    @property
    def metadata(self) -> dict[str, Any]:
        """The metadata for the column."""
        return {}

    def num_chunks(self) -> int:
        """Return the number of chunks the column consists of."""
        return self._col.n_chunks()

    def get_chunks(self, n_chunks: int | None = None) -> Iterator[PolarsColumn]:
        """
        Return an iterator yielding the column chunks.

        Parameters
        ----------
        n_chunks
            The number of chunks to return. Must be a multiple of the number of chunks
            in the column.

        Notes
        -----
        When `n_chunks` is higher than the number of chunks in the column, a slice
        must be performed that is not on the chunk boundary. This will trigger some
        compute if the column contains null values or if the column is of data type
        boolean.
        """
        total_n_chunks = self.num_chunks()
        chunks = self._col.get_chunks()

        if (n_chunks is None) or (n_chunks == total_n_chunks):
            for chunk in chunks:
                yield PolarsColumn(chunk, allow_copy=self._allow_copy)

        elif (n_chunks <= 0) or (n_chunks % total_n_chunks != 0):
            msg = (
                "`n_chunks` must be a multiple of the number of chunks of this column"
                f" ({total_n_chunks})"
            )
            raise ValueError(msg)

        else:
            subchunks_per_chunk = n_chunks // total_n_chunks
            for chunk in chunks:
                size = len(chunk)
                step = size // subchunks_per_chunk
                if size % subchunks_per_chunk != 0:
                    step += 1
                for start in range(0, step * subchunks_per_chunk, step):
                    yield PolarsColumn(
                        chunk[start : start + step], allow_copy=self._allow_copy
                    )

    def get_buffers(self) -> ColumnBuffers:
        """Return a dictionary containing the underlying buffers."""
        dtype = self._col.dtype

        if dtype == String and not self._allow_copy:
            msg = "string buffers must be converted"
            raise CopyNotAllowedError(msg)

        buffers = self._col._get_buffers()

        return {
            "data": self._wrap_data_buffer(buffers["values"]),
            "validity": self._wrap_validity_buffer(buffers["validity"]),
            "offsets": self._wrap_offsets_buffer(buffers["offsets"]),
        }

    def _wrap_data_buffer(self, buffer: Series) -> tuple[PolarsBuffer, Dtype]:
        interchange_buffer = PolarsBuffer(buffer, allow_copy=self._allow_copy)
        dtype = polars_dtype_to_dtype(buffer.dtype)
        return interchange_buffer, dtype

    def _wrap_validity_buffer(
        self, buffer: Series | None
    ) -> tuple[PolarsBuffer, Dtype] | None:
        if buffer is None:
            return None

        interchange_buffer = PolarsBuffer(buffer, allow_copy=self._allow_copy)
        dtype = (DtypeKind.BOOL, 1, "b", Endianness.NATIVE)
        return interchange_buffer, dtype

    def _wrap_offsets_buffer(
        self, buffer: Series | None
    ) -> tuple[PolarsBuffer, Dtype] | None:
        if buffer is None:
            return None

        interchange_buffer = PolarsBuffer(buffer, allow_copy=self._allow_copy)
        dtype = (DtypeKind.INT, 64, "l", Endianness.NATIVE)
        return interchange_buffer, dtype
