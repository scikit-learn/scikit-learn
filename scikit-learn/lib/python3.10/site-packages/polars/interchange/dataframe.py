from __future__ import annotations

from collections.abc import Sequence
from itertools import accumulate
from typing import TYPE_CHECKING

from polars.interchange.column import PolarsColumn
from polars.interchange.protocol import CopyNotAllowedError
from polars.interchange.protocol import DataFrame as InterchangeDataFrame

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Any

    from polars import DataFrame


class PolarsDataFrame(InterchangeDataFrame):
    """
    A dataframe object backed by a Polars DataFrame.

    Parameters
    ----------
    df
        The Polars DataFrame backing the dataframe object.
    allow_copy
        Allow data to be copied during operations on this column. If set to `False`,
        a RuntimeError is raised if data would be copied.
    """

    version = 0

    def __init__(self, df: DataFrame, *, allow_copy: bool = True) -> None:
        self._df = df
        self._allow_copy = allow_copy

    def __dataframe__(
        self,
        nan_as_null: bool = False,  # noqa: FBT001
        allow_copy: bool = True,  # noqa: FBT001
    ) -> PolarsDataFrame:
        """
        Construct a new dataframe object, potentially changing the parameters.

        Parameters
        ----------
        nan_as_null
            Overwrite null values in the data with `NaN`.

            .. warning::
                This functionality has not been implemented and the parameter will be
                removed in a future version.
                Setting this to `True` will raise a `NotImplementedError`.
        allow_copy
            Allow memory to be copied to perform the conversion. If set to `False`,
            causes conversions that are not zero-copy to fail.
        """
        if nan_as_null:
            msg = (
                "functionality for `nan_as_null` has not been implemented and the"
                " parameter will be removed in a future version"
                "\n\nUse the default `nan_as_null=False`."
            )
            raise NotImplementedError(msg)
        return PolarsDataFrame(self._df, allow_copy=allow_copy)

    @property
    def metadata(self) -> dict[str, Any]:
        """The metadata for the dataframe."""
        return {}

    def num_columns(self) -> int:
        """Return the number of columns in the dataframe."""
        return self._df.width

    def num_rows(self) -> int:
        """Return the number of rows in the dataframe."""
        return self._df.height

    def num_chunks(self) -> int:
        """
        Return the number of chunks the dataframe consists of.

        It is possible for a Polars DataFrame to consist of columns with a varying
        number of chunks. This method returns the number of chunks of the first
        column.

        See Also
        --------
        polars.dataframe.frame.DataFrame.n_chunks
        """
        return self._df.n_chunks("first")

    def column_names(self) -> list[str]:
        """Return the column names."""
        return self._df.columns

    def get_column(self, i: int) -> PolarsColumn:
        """
        Return the column at the indicated position.

        Parameters
        ----------
        i
            Index of the column.
        """
        s = self._df.to_series(i)
        return PolarsColumn(s, allow_copy=self._allow_copy)

    def get_column_by_name(self, name: str) -> PolarsColumn:
        """
        Return the column with the given name.

        Parameters
        ----------
        name
            Name of the column.
        """
        s = self._df.get_column(name)
        return PolarsColumn(s, allow_copy=self._allow_copy)

    def get_columns(self) -> Iterator[PolarsColumn]:
        """Return an iterator yielding the columns."""
        for column in self._df.get_columns():
            yield PolarsColumn(column, allow_copy=self._allow_copy)

    def select_columns(self, indices: Sequence[int]) -> PolarsDataFrame:
        """
        Create a new dataframe by selecting a subset of columns by index.

        Parameters
        ----------
        indices
            Column indices
        """
        if not isinstance(indices, Sequence):
            msg = "`indices` is not a sequence"
            raise TypeError(msg)
        if not isinstance(indices, list):
            indices = list(indices)

        return PolarsDataFrame(
            self._df[:, indices],
            allow_copy=self._allow_copy,
        )

    def select_columns_by_name(self, names: Sequence[str]) -> PolarsDataFrame:
        """
        Create a new dataframe by selecting a subset of columns by name.

        Parameters
        ----------
        names
            Column names.
        """
        if not isinstance(names, Sequence):
            msg = "`names` is not a sequence"
            raise TypeError(msg)

        return PolarsDataFrame(
            self._df.select(names),
            allow_copy=self._allow_copy,
        )

    def get_chunks(self, n_chunks: int | None = None) -> Iterator[PolarsDataFrame]:
        """
        Return an iterator yielding the chunks of the dataframe.

        Parameters
        ----------
        n_chunks
            The number of chunks to return. Must be a multiple of the number of chunks
            in the dataframe. If set to `None` (default), returns all chunks.

        Notes
        -----
        When the columns in the dataframe are chunked unevenly, or when `n_chunks` is
        higher than the number of chunks in the dataframe, a slice must be performed
        that is not on the chunk boundary. This will trigger some compute for columns
        that contain null values and boolean columns.
        """
        total_n_chunks = self.num_chunks()
        chunks = self._get_chunks_from_col_chunks()

        if (n_chunks is None) or (n_chunks == total_n_chunks):
            for chunk in chunks:
                yield PolarsDataFrame(chunk, allow_copy=self._allow_copy)

        elif (n_chunks <= 0) or (n_chunks % total_n_chunks != 0):
            msg = (
                "`n_chunks` must be a multiple of the number of chunks of this"
                f" dataframe ({total_n_chunks})"
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
                    yield PolarsDataFrame(
                        chunk[start : start + step, :],
                        allow_copy=self._allow_copy,
                    )

    def _get_chunks_from_col_chunks(self) -> Iterator[DataFrame]:
        """
        Return chunks of this dataframe according to the chunks of the first column.

        If columns are not all chunked identically, they will be rechunked like the
        first column. If copy is not allowed, this raises a RuntimeError.
        """
        col_chunks = self.get_column(0).get_chunks()
        chunk_sizes = [chunk.size() for chunk in col_chunks]
        starts = [0] + list(accumulate(chunk_sizes))

        for i in range(len(starts) - 1):
            start, end = starts[i : i + 2]
            chunk = self._df[start:end, :]

            if not all(x == 1 for x in chunk.n_chunks("all")):
                if not self._allow_copy:
                    msg = "unevenly chunked columns must be rechunked"
                    raise CopyNotAllowedError(msg)
                chunk = chunk.rechunk()

            yield chunk
