""" pickle compat """
import warnings

import numpy as np
from numpy.lib.format import read_array, write_array
from pandas.compat import BytesIO, cPickle as pkl, pickle_compat as pc, PY3
from pandas.core.dtypes.common import is_datetime64_dtype, _NS_DTYPE
from pandas.io.common import _get_handle, _infer_compression, _stringify_path


def to_pickle(obj, path, compression='infer', protocol=pkl.HIGHEST_PROTOCOL):
    """
    Pickle (serialize) object to file.

    Parameters
    ----------
    obj : any object
        Any python object.
    path : str
        File path where the pickled object will be stored.
    compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default 'infer'
        A string representing the compression to use in the output file. By
        default, infers from the file extension in specified path.

        .. versionadded:: 0.20.0
    protocol : int
        Int which indicates which protocol should be used by the pickler,
        default HIGHEST_PROTOCOL (see [1], paragraph 12.1.2). The possible
        values for this parameter depend on the version of Python. For Python
        2.x, possible values are 0, 1, 2. For Python>=3.0, 3 is a valid value.
        For Python >= 3.4, 4 is a valid value. A negative value for the
        protocol parameter is equivalent to setting its value to
        HIGHEST_PROTOCOL.

        .. [1] https://docs.python.org/3/library/pickle.html
        .. versionadded:: 0.21.0

    See Also
    --------
    read_pickle : Load pickled pandas object (or any object) from file.
    DataFrame.to_hdf : Write DataFrame to an HDF5 file.
    DataFrame.to_sql : Write DataFrame to a SQL database.
    DataFrame.to_parquet : Write a DataFrame to the binary parquet format.

    Examples
    --------
    >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
    >>> original_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9
    >>> pd.to_pickle(original_df, "./dummy.pkl")

    >>> unpickled_df = pd.read_pickle("./dummy.pkl")
    >>> unpickled_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9

    >>> import os
    >>> os.remove("./dummy.pkl")
    """
    path = _stringify_path(path)
    inferred_compression = _infer_compression(path, compression)
    f, fh = _get_handle(path, 'wb',
                        compression=inferred_compression,
                        is_text=False)
    if protocol < 0:
        protocol = pkl.HIGHEST_PROTOCOL
    try:
        f.write(pkl.dumps(obj, protocol=protocol))
    finally:
        for _f in fh:
            _f.close()


def read_pickle(path, compression='infer'):
    """
    Load pickled pandas object (or any object) from file.

    .. warning::

       Loading pickled data received from untrusted sources can be
       unsafe. See `here <https://docs.python.org/3/library/pickle.html>`__.

    Parameters
    ----------
    path : str
        File path where the pickled object will be loaded.
    compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default 'infer'
        For on-the-fly decompression of on-disk data. If 'infer', then use
        gzip, bz2, xz or zip if path ends in '.gz', '.bz2', '.xz',
        or '.zip' respectively, and no decompression otherwise.
        Set to None for no decompression.

        .. versionadded:: 0.20.0

    Returns
    -------
    unpickled : type of object stored in file

    See Also
    --------
    DataFrame.to_pickle : Pickle (serialize) DataFrame object to file.
    Series.to_pickle : Pickle (serialize) Series object to file.
    read_hdf : Read HDF5 file into a DataFrame.
    read_sql : Read SQL query or database table into a DataFrame.
    read_parquet : Load a parquet object, returning a DataFrame.

    Examples
    --------
    >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
    >>> original_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9
    >>> pd.to_pickle(original_df, "./dummy.pkl")

    >>> unpickled_df = pd.read_pickle("./dummy.pkl")
    >>> unpickled_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9

    >>> import os
    >>> os.remove("./dummy.pkl")
    """
    path = _stringify_path(path)
    inferred_compression = _infer_compression(path, compression)

    def read_wrapper(func):
        # wrapper file handle open/close operation
        f, fh = _get_handle(path, 'rb',
                            compression=inferred_compression,
                            is_text=False)
        try:
            return func(f)
        finally:
            for _f in fh:
                _f.close()

    def try_read(path, encoding=None):
        # try with cPickle
        # try with current pickle, if we have a Type Error then
        # try with the compat pickle to handle subclass changes
        # pass encoding only if its not None as py2 doesn't handle
        # the param

        # cpickle
        # GH 6899
        try:
            with warnings.catch_warnings(record=True):
                # We want to silencce any warnings about, e.g. moved modules.
                return read_wrapper(lambda f: pkl.load(f))
        except Exception:
            # reg/patched pickle
            try:
                return read_wrapper(
                    lambda f: pc.load(f, encoding=encoding, compat=False))
            # compat pickle
            except:
                return read_wrapper(
                    lambda f: pc.load(f, encoding=encoding, compat=True))
    try:
        return try_read(path)
    except:
        if PY3:
            return try_read(path, encoding='latin1')
        raise


# compat with sparse pickle / unpickle


def _pickle_array(arr):
    arr = arr.view(np.ndarray)

    buf = BytesIO()
    write_array(buf, arr)

    return buf.getvalue()


def _unpickle_array(bytes):
    arr = read_array(BytesIO(bytes))

    # All datetimes should be stored as M8[ns].  When unpickling with
    # numpy1.6, it will read these as M8[us].  So this ensures all
    # datetime64 types are read as MS[ns]
    if is_datetime64_dtype(arr):
        arr = arr.view(_NS_DTYPE)

    return arr
