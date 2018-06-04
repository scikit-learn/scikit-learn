from __future__ import absolute_import, division, print_function

from ..core import DataFrame
from ...base import tokenize
from ...bytes.core import get_fs_token_paths
from ...utils import import_required
from .utils import _get_pyarrow_dtypes, _meta_from_dtypes

__all__ = ('read_orc',)


def _read_orc_stripe(fs, path, stripe, columns=None):
    """Pull out specific data from specific part of ORC file"""
    orc = import_required('pyarrow.orc', 'Please install pyarrow >= 0.9.0')
    with fs.open(path, 'rb') as f:
        o = orc.ORCFile(f)
        return o.read_stripe(stripe, columns).to_pandas()


def read_orc(path, columns=None, storage_options=None):
    """Read dataframe from ORC file(s)

    Parameters
    ----------
    path: str or list(str)
        Location of file(s), which can be a full URL with protocol specifier,
        and may include glob character if a single string.
    columns: None or list(str)
        Columns to load. If None, loads all.
    storage_options: None or dict
        Further parameters to pass to the bytes backend.

    Returns
    -------
    Dask.DataFrame (even if there is only one column)

    Examples
    --------
    >>> df = dd.read_orc('https://github.com/apache/orc/raw/'
    ...                  'master/examples/demo-11-zlib.orc')  # doctest: +SKIP
    """
    orc = import_required('pyarrow.orc', 'Please install pyarrow >= 0.9.0')
    storage_options = storage_options or {}
    fs, fs_token, paths = get_fs_token_paths(path, mode='rb',
                                             storage_options=storage_options)
    schema = None
    nstripes_per_file = []
    for path in paths:
        with fs.open(path, 'rb') as f:
            o = orc.ORCFile(f)
            if schema is None:
                schema = o.schema
            elif schema != o.schema:
                raise ValueError("Incompatible schemas while parsing ORC files")
            nstripes_per_file.append(o.nstripes)
    schema = _get_pyarrow_dtypes(schema, categories=None)
    if columns is not None:
        ex = set(columns) - set(schema)
        if ex:
            raise ValueError("Requested columns (%s) not in schema (%s)" % (
                ex, set(schema)
            ))
    else:
        columns = list(schema)
    meta = _meta_from_dtypes(columns, schema, [], [])

    name = 'read-orc-' + tokenize(fs_token, path, columns)
    dsk = {}
    N = 0
    for path, n in zip(paths, nstripes_per_file):
        for stripe in range(n):
            dsk[(name, N)] = (_read_orc_stripe, fs, path, stripe, columns)
            N += 1

    return DataFrame(dsk, name, meta, [None] * (len(dsk) + 1))
