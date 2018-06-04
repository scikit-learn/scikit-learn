""" parquet compat """

from warnings import catch_warnings
from distutils.version import LooseVersion
from pandas import DataFrame, RangeIndex, Int64Index, get_option
from pandas.compat import string_types
import pandas.core.common as com
from pandas.io.common import get_filepath_or_buffer, is_s3_url


def get_engine(engine):
    """ return our implementation """

    if engine == 'auto':
        engine = get_option('io.parquet.engine')

    if engine == 'auto':
        # try engines in this order
        try:
            return PyArrowImpl()
        except ImportError:
            pass

        try:
            return FastParquetImpl()
        except ImportError:
            pass

        raise ImportError("Unable to find a usable engine; "
                          "tried using: 'pyarrow', 'fastparquet'.\n"
                          "pyarrow or fastparquet is required for parquet "
                          "support")

    if engine not in ['pyarrow', 'fastparquet']:
        raise ValueError("engine must be one of 'pyarrow', 'fastparquet'")

    if engine == 'pyarrow':
        return PyArrowImpl()
    elif engine == 'fastparquet':
        return FastParquetImpl()


class BaseImpl(object):

    api = None  # module

    @staticmethod
    def validate_dataframe(df):

        if not isinstance(df, DataFrame):
            raise ValueError("to_parquet only supports IO with DataFrames")

        # must have value column names (strings only)
        if df.columns.inferred_type not in {'string', 'unicode'}:
            raise ValueError("parquet must have string column names")

        # index level names must be strings
        valid_names = all(
            isinstance(name, string_types)
            for name in df.index.names
            if name is not None
        )
        if not valid_names:
            raise ValueError("Index level names must be strings")

    def write(self, df, path, compression, **kwargs):
        raise com.AbstractMethodError(self)

    def read(self, path, columns=None, **kwargs):
        raise com.AbstractMethodError(self)


class PyArrowImpl(BaseImpl):

    def __init__(self):
        # since pandas is a dependency of pyarrow
        # we need to import on first use
        try:
            import pyarrow
            import pyarrow.parquet
        except ImportError:
            raise ImportError(
                "pyarrow is required for parquet support\n\n"
                "you can install via conda\n"
                "conda install pyarrow -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U pyarrow\n"
            )
        if LooseVersion(pyarrow.__version__) < '0.4.1':
            raise ImportError(
                "pyarrow >= 0.4.1 is required for parquet support\n\n"
                "you can install via conda\n"
                "conda install pyarrow -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U pyarrow\n"
            )

        self._pyarrow_lt_060 = (
            LooseVersion(pyarrow.__version__) < LooseVersion('0.6.0'))
        self._pyarrow_lt_070 = (
            LooseVersion(pyarrow.__version__) < LooseVersion('0.7.0'))

        self.api = pyarrow

    def write(self, df, path, compression='snappy',
              coerce_timestamps='ms', **kwargs):
        self.validate_dataframe(df)
        if self._pyarrow_lt_070:
            self._validate_write_lt_070(df)
        path, _, _, _ = get_filepath_or_buffer(path, mode='wb')

        if self._pyarrow_lt_060:
            table = self.api.Table.from_pandas(df, timestamps_to_ms=True)
            self.api.parquet.write_table(
                table, path, compression=compression, **kwargs)

        else:
            table = self.api.Table.from_pandas(df)
            self.api.parquet.write_table(
                table, path, compression=compression,
                coerce_timestamps=coerce_timestamps, **kwargs)

    def read(self, path, columns=None, **kwargs):
        path, _, _, should_close = get_filepath_or_buffer(path)
        if self._pyarrow_lt_070:
            result = self.api.parquet.read_pandas(path, columns=columns,
                                                  **kwargs).to_pandas()
        else:
            kwargs['use_pandas_metadata'] = True
            result = self.api.parquet.read_table(path, columns=columns,
                                                 **kwargs).to_pandas()
        if should_close:
            try:
                path.close()
            except:  # noqa: flake8
                pass

        return result

    def _validate_write_lt_070(self, df):
        # Compatibility shim for pyarrow < 0.7.0
        # TODO: Remove in pandas 0.23.0
        from pandas.core.indexes.multi import MultiIndex
        if isinstance(df.index, MultiIndex):
            msg = (
                "Multi-index DataFrames are only supported "
                "with pyarrow >= 0.7.0"
            )
            raise ValueError(msg)
        # Validate index
        if not isinstance(df.index, Int64Index):
            msg = (
                "pyarrow < 0.7.0 does not support serializing {} for the "
                "index; you can .reset_index() to make the index into "
                "column(s), or install the latest version of pyarrow or "
                "fastparquet."
            )
            raise ValueError(msg.format(type(df.index)))
        if not df.index.equals(RangeIndex(len(df))):
            raise ValueError(
                "pyarrow < 0.7.0 does not support serializing a non-default "
                "index; you can .reset_index() to make the index into "
                "column(s), or install the latest version of pyarrow or "
                "fastparquet."
            )
        if df.index.name is not None:
            raise ValueError(
                "pyarrow < 0.7.0 does not serialize indexes with a name; you "
                "can set the index.name to None or install the latest version "
                "of pyarrow or fastparquet."
            )


class FastParquetImpl(BaseImpl):

    def __init__(self):
        # since pandas is a dependency of fastparquet
        # we need to import on first use
        try:
            import fastparquet
        except ImportError:
            raise ImportError(
                "fastparquet is required for parquet support\n\n"
                "you can install via conda\n"
                "conda install fastparquet -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U fastparquet"
            )
        if LooseVersion(fastparquet.__version__) < '0.1.0':
            raise ImportError(
                "fastparquet >= 0.1.0 is required for parquet "
                "support\n\n"
                "you can install via conda\n"
                "conda install fastparquet -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U fastparquet"
            )
        self.api = fastparquet

    def write(self, df, path, compression='snappy', **kwargs):
        self.validate_dataframe(df)
        # thriftpy/protocol/compact.py:339:
        # DeprecationWarning: tostring() is deprecated.
        # Use tobytes() instead.

        if is_s3_url(path):
            # path is s3:// so we need to open the s3file in 'wb' mode.
            # TODO: Support 'ab'

            path, _, _, _ = get_filepath_or_buffer(path, mode='wb')
            # And pass the opened s3file to the fastparquet internal impl.
            kwargs['open_with'] = lambda path, _: path
        else:
            path, _, _, _ = get_filepath_or_buffer(path)

        with catch_warnings(record=True):
            self.api.write(path, df,
                           compression=compression, **kwargs)

    def read(self, path, columns=None, **kwargs):
        if is_s3_url(path):
            # When path is s3:// an S3File is returned.
            # We need to retain the original path(str) while also
            # pass the S3File().open function to fsatparquet impl.
            s3, _, _, should_close = get_filepath_or_buffer(path)
            try:
                parquet_file = self.api.ParquetFile(path, open_with=s3.s3.open)
            finally:
                s3.close()
        else:
            path, _, _, _ = get_filepath_or_buffer(path)
            parquet_file = self.api.ParquetFile(path)

        return parquet_file.to_pandas(columns=columns, **kwargs)


def to_parquet(df, path, engine='auto', compression='snappy', **kwargs):
    """
    Write a DataFrame to the parquet format.

    Parameters
    ----------
    df : DataFrame
    path : string
        File path
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.
    compression : {'snappy', 'gzip', 'brotli', None}, default 'snappy'
        Name of the compression to use. Use ``None`` for no compression.
    kwargs
        Additional keyword arguments passed to the engine
    """
    impl = get_engine(engine)
    return impl.write(df, path, compression=compression, **kwargs)


def read_parquet(path, engine='auto', columns=None, **kwargs):
    """
    Load a parquet object from the file path, returning a DataFrame.

    .. versionadded 0.21.0

    Parameters
    ----------
    path : string
        File path
    columns: list, default=None
        If not None, only these columns will be read from the file.

        .. versionadded 0.21.1
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.
    kwargs are passed to the engine

    Returns
    -------
    DataFrame

    """

    impl = get_engine(engine)
    return impl.read(path, columns=columns, **kwargs)
