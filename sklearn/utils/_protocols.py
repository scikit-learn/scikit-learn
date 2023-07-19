"""Defines Python Protocol for the DataFrame InterChange Protocol.

See more in https://data-apis.org/dataframe-protocol/latest/API.html
"""
from typing import Protocol, runtime_checkable


@runtime_checkable
class DataFrameInterchangeProtocol(Protocol):
    def __dataframe__(self, nan_as_null, allow_copy):
        """Construct a new exchange object, potentially changing the parameters."""

    @property
    def metadata(self):
        """Metadata for data frame."""

    def num_columns(self):
        """Return the number of columns in the DataFrame."""

    def num_rows(self):
        """Return the number of rows in the DataFrame, if available."""

    def num_chunks(self):
        """Return the number of chunks the DataFrame consists of."""

    def column_names(self):
        """Return an iterator yielding the column names."""

    def get_column(self, i):
        """Return the column at the indicated position."""

    def get_column_by_name(self, name):
        """Return the column whose name is the indicated name."""

    def get_columns(self):
        """Return an iterator yielding the columns."""

    def select_columns(self, indices):
        """Create a new DataFrame by selecting a subset of columns by index."""

    def select_columns_by_name(self, names):
        """Create a new DataFrame by selecting a subset of columns by name."""

    def get_chunks(self, n_chunks):
        """Return an iterator yielding the chunks."""
