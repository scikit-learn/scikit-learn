from typing import Protocol, runtime_checkable


@runtime_checkable
class DataFrameInterchangeProtocol(Protocol):
    def __dataframe__(self, nan_as_null, allow_copy):
        ...

    @property
    def metadata(self):
        ...

    def num_columns(self):
        ...

    def num_rows(self):
        ...

    def num_chunks(self):
        ...

    def column_names(self):
        ...

    def get_column(self, i):
        ...

    def get_column_by_name(self, name):
        ...

    def get_columns(self):
        ...

    def select_columns(self, indices):
        ...

    def select_columns_by_name(self, names):
        ...

    def get_chunks(self, n_chunks):
        ...
