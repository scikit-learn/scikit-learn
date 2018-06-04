from __future__ import absolute_import, division, print_function

from fnmatch import fnmatch
from glob import glob
import os
import uuid
from warnings import warn

import pandas as pd
from toolz import merge

from .io import _link
from ..core import DataFrame, new_dd_object
from ... import multiprocessing
from ...base import tokenize, compute_as_if_collection
from ...bytes.utils import build_name_function
from ...compatibility import PY3
from ...context import _globals
from ...delayed import Delayed, delayed
from ...local import get_sync
from ...utils import effective_get, get_scheduler_lock


def _pd_to_hdf(pd_to_hdf, lock, args, kwargs=None):
    """ A wrapper function around pd_to_hdf that enables locking"""

    if lock:
        lock.acquire()
    try:
        pd_to_hdf(*args, **kwargs)
    finally:
        if lock:
            lock.release()

    return None


def to_hdf(df, path, key, mode='a', append=False, get=None,
           name_function=None, compute=True, lock=None, dask_kwargs={},
           **kwargs):
    """ Store Dask Dataframe to Hierarchical Data Format (HDF) files

    This is a parallel version of the Pandas function of the same name.  Please
    see the Pandas docstring for more detailed information about shared keyword
    arguments.

    This function differs from the Pandas version by saving the many partitions
    of a Dask DataFrame in parallel, either to many files, or to many datasets
    within the same file.  You may specify this parallelism with an asterix
    ``*`` within the filename or datapath, and an optional ``name_function``.
    The asterix will be replaced with an increasing sequence of integers
    starting from ``0`` or with the result of calling ``name_function`` on each
    of those integers.

    This function only supports the Pandas ``'table'`` format, not the more
    specialized ``'fixed'`` format.

    Parameters
    ----------
    path: string
        Path to a target filename.  May contain a ``*`` to denote many filenames
    key: string
        Datapath within the files.  May contain a ``*`` to denote many locations
    name_function: function
        A function to convert the ``*`` in the above options to a string.
        Should take in a number from 0 to the number of partitions and return a
        string. (see examples below)
    compute: bool
        Whether or not to execute immediately.  If False then this returns a
        ``dask.Delayed`` value.
    lock: Lock, optional
        Lock to use to prevent concurrency issues.  By default a
        ``threading.Lock``, ``multiprocessing.Lock`` or ``SerializableLock``
        will be used depending on your scheduler if a lock is required. See
        dask.utils.get_scheduler_lock for more information about lock
        selection.
    **other:
        See pandas.to_hdf for more information

    Examples
    --------
    Save Data to a single file

    >>> df.to_hdf('output.hdf', '/data')            # doctest: +SKIP

    Save data to multiple datapaths within the same file:

    >>> df.to_hdf('output.hdf', '/data-*')          # doctest: +SKIP

    Save data to multiple files:

    >>> df.to_hdf('output-*.hdf', '/data')          # doctest: +SKIP

    Save data to multiple files, using the multiprocessing scheduler:

    >>> df.to_hdf('output-*.hdf', '/data', get=dask.multiprocessing.get) # doctest: +SKIP

    Specify custom naming scheme.  This writes files as
    '2000-01-01.hdf', '2000-01-02.hdf', '2000-01-03.hdf', etc..

    >>> from datetime import date, timedelta
    >>> base = date(year=2000, month=1, day=1)
    >>> def name_function(i):
    ...     ''' Convert integer 0 to n to a string '''
    ...     return base + timedelta(days=i)

    >>> df.to_hdf('*.hdf', '/data', name_function=name_function) # doctest: +SKIP

    Returns
    -------
    None: if compute == True
    delayed value: if compute == False

    See Also
    --------
    read_hdf:
    to_parquet:
    """
    name = 'to-hdf-' + uuid.uuid1().hex

    pd_to_hdf = getattr(df._partition_type, 'to_hdf')

    single_file = True
    single_node = True

    # if path is string, format using i_name
    if isinstance(path, str):
        if path.count('*') + key.count('*') > 1:
            raise ValueError("A maximum of one asterisk is accepted in file "
                             "path and dataset key")

        fmt_obj = lambda path, i_name: path.replace('*', i_name)

        if '*' in path:
            single_file = False
    else:
        if key.count('*') > 1:
            raise ValueError("A maximum of one asterisk is accepted in "
                             "dataset key")

        fmt_obj = lambda path, _: path

    if '*' in key:
        single_node = False

    if 'format' in kwargs and kwargs['format'] not in ['t', 'table']:
        raise ValueError("Dask only support 'table' format in hdf files.")

    if mode not in ('a', 'w', 'r+'):
        raise ValueError("Mode must be one of 'a', 'w' or 'r+'")

    if name_function is None:
        name_function = build_name_function(df.npartitions - 1)

    # we guarantee partition order is preserved when its saved and read
    # so we enforce name_function to maintain the order of its input.
    if not (single_file and single_node):
        formatted_names = [name_function(i) for i in range(df.npartitions)]
        if formatted_names != sorted(formatted_names):
            warn("To preserve order between partitions name_function "
                 "must preserve the order of its input")

    # If user did not specify scheduler and write is sequential default to the
    # sequential scheduler. otherwise let the _get method choose the scheduler
    if get is None and 'get' not in _globals and single_node and single_file:
        get = get_sync

    # handle lock default based on whether we're writing to a single entity
    _actual_get = effective_get(get, df)
    if lock is None:
        if not single_node:
            lock = True
        elif not single_file and _actual_get is not multiprocessing.get:
            # if we're writing to multiple files with the multiprocessing
            # scheduler we don't need to lock
            lock = True
        else:
            lock = False
    if lock:
        lock = get_scheduler_lock(get, df)

    kwargs.update({'format': 'table', 'mode': mode, 'append': append})

    dsk = dict()

    i_name = name_function(0)
    dsk[(name, 0)] = (_pd_to_hdf, pd_to_hdf, lock,
                      [(df._name, 0), fmt_obj(path, i_name),
                       key.replace('*', i_name)], kwargs)

    kwargs2 = kwargs.copy()
    if single_file:
        kwargs2['mode'] = 'a'
    if single_node:
        kwargs2['append'] = True

    filenames = []
    for i in range(0,df.npartitions):
        i_name = name_function(i)
        filenames.append(fmt_obj(path, i_name))

    for i in range(1, df.npartitions):
        i_name = name_function(i)
        task = (_pd_to_hdf, pd_to_hdf, lock,
                [(df._name, i), fmt_obj(path, i_name),
                    key.replace('*', i_name)], kwargs2)
        if single_file:
            link_dep = i - 1 if single_node else 0
            task = (_link, (name, link_dep), task)
        dsk[(name, i)] = task

    dsk = merge(df.dask, dsk)
    if single_file and single_node:
        keys = [(name, df.npartitions - 1)]
    else:
        keys = [(name, i) for i in range(df.npartitions)]

    if compute:
        compute_as_if_collection(DataFrame, dsk, keys, get=get, **dask_kwargs)
        return filenames
    else:
        return delayed([Delayed(k, dsk) for k in keys])


dont_use_fixed_error_message = """
This HDFStore is not partitionable and can only be use monolithically with
pandas.  In the future when creating HDFStores use the ``format='table'``
option to ensure that your dataset can be parallelized"""

read_hdf_error_msg = """
The start and stop keywords are not supported when reading from more than
one file/dataset.

The combination is ambiguous because it could be interpreted as the starting
and stopping index per file, or starting and stopping index of the global
dataset."""


def _read_single_hdf(path, key, start=0, stop=None, columns=None,
                     chunksize=int(1e6), sorted_index=False, lock=None,
                     mode='a'):
    """
    Read a single hdf file into a dask.dataframe. Used for each file in
    read_hdf.
    """
    def get_keys_stops_divisions(path, key, stop, sorted_index, chunksize):
        """
        Get the "keys" or group identifiers which match the given key, which
        can contain wildcards. This uses the hdf file identified by the
        given path. Also get the index of the last row of data for each matched
        key.
        """
        with pd.HDFStore(path, mode=mode) as hdf:
            keys = [k for k in hdf.keys() if fnmatch(k, key)]
            stops = []
            divisions = []
            for k in keys:
                storer = hdf.get_storer(k)
                if storer.format_type != 'table':
                    raise TypeError(dont_use_fixed_error_message)
                if stop is None:
                    stops.append(storer.nrows)
                elif stop > storer.nrows:
                    raise ValueError("Stop keyword exceeds dataset number "
                                     "of rows ({})".format(storer.nrows))
                else:
                    stops.append(stop)
                if sorted_index:
                    division = [storer.read_column('index', start=start, stop=start + 1)[0]
                                for start in range(0, storer.nrows, chunksize)]
                    division_end = storer.read_column('index',
                                                      start=storer.nrows - 1,
                                                      stop=storer.nrows)[0]

                    division.append(division_end)
                    divisions.append(division)
                else:
                    divisions.append(None)

        return keys, stops, divisions

    def one_path_one_key(path, key, start, stop, columns, chunksize, division, lock):
        """
        Get the data frame corresponding to one path and one key (which should
        not contain any wildcards).
        """
        empty = pd.read_hdf(path, key, mode=mode, stop=0)
        if columns is not None:
            empty = empty[columns]

        token = tokenize((path, os.path.getmtime(path), key, start,
                          stop, empty, chunksize, division))
        name = 'read-hdf-' + token
        if empty.ndim == 1:
            base = {'name': empty.name, 'mode': mode}
        else:
            base = {'columns': empty.columns, 'mode': mode}

        if start >= stop:
            raise ValueError("Start row number ({}) is above or equal to stop "
                             "row number ({})".format(start, stop))

        def update(s):
            new = base.copy()
            new.update({'start': s, 'stop': s + chunksize})
            return new

        dsk = dict(((name, i), (_pd_read_hdf, path, key, lock,
                                update(s)))
                   for i, s in enumerate(range(start, stop, chunksize)))

        if division:
            divisions = division
        else:
            divisions = [None] * (len(dsk) + 1)

        return new_dd_object(dsk, name, empty, divisions)

    keys, stops, divisions = get_keys_stops_divisions(path, key, stop, sorted_index, chunksize)
    if (start != 0 or stop is not None) and len(keys) > 1:
        raise NotImplementedError(read_hdf_error_msg)
    from ..multi import concat
    return concat([one_path_one_key(path, k, start, s, columns, chunksize, d, lock)
                   for k, s, d in zip(keys, stops, divisions)])


def _pd_read_hdf(path, key, lock, kwargs):
    """ Read from hdf5 file with a lock """
    if lock:
        lock.acquire()
    try:
        result = pd.read_hdf(path, key, **kwargs)
    finally:
        if lock:
            lock.release()
    return result


def read_hdf(pattern, key, start=0, stop=None, columns=None,
             chunksize=1000000, sorted_index=False, lock=True, mode='a'):
    """
    Read HDF files into a Dask DataFrame

    Read hdf files into a dask dataframe. This function is like
    ``pandas.read_hdf``, except it can read from a single large file, or from
    multiple files, or from multiple keys from the same file.

    Parameters
    ----------
    pattern : string, list
        File pattern (string), buffer to read from, or list of file
        paths. Can contain wildcards.
    key : group identifier in the store. Can contain wildcards
    start : optional, integer (defaults to 0), row number to start at
    stop : optional, integer (defaults to None, the last row), row number to
        stop at
    columns : list of columns, optional
        A list of columns that if not None, will limit the return
        columns (default is None)
    chunksize : positive integer, optional
        Maximal number of rows per partition (default is 1000000).
    sorted_index : boolean, optional
        Option to specify whether or not the input hdf files have a sorted
        index (default is False).
    lock : boolean, optional
        Option to use a lock to prevent concurrency issues (default is True).
    mode : {'a', 'r', 'r+'}, default 'a'. Mode to use when opening file(s).
        'r'
            Read-only; no data can be modified.
        'a'
            Append; an existing file is opened for reading and writing,
            and if the file does not exist it is created.
        'r+'
            It is similar to 'a', but the file must already exist.

    Returns
    -------
    dask.DataFrame

    Examples
    --------
    Load single file

    >>> dd.read_hdf('myfile.1.hdf5', '/x')  # doctest: +SKIP

    Load multiple files

    >>> dd.read_hdf('myfile.*.hdf5', '/x')  # doctest: +SKIP

    >>> dd.read_hdf(['myfile.1.hdf5', 'myfile.2.hdf5'], '/x')  # doctest: +SKIP

    Load multiple datasets

    >>> dd.read_hdf('myfile.1.hdf5', '/*')  # doctest: +SKIP
    """
    if lock is True:
        lock = get_scheduler_lock()

    key = key if key.startswith('/') else '/' + key
    if isinstance(pattern, str):
        paths = sorted(glob(pattern))
    else:
        paths = pattern
    if (start != 0 or stop is not None) and len(paths) > 1:
        raise NotImplementedError(read_hdf_error_msg)
    if chunksize <= 0:
        raise ValueError("Chunksize must be a positive integer")
    if (start != 0 or stop is not None) and sorted_index:
        raise ValueError("When assuming pre-partitioned data, data must be "
                         "read in its entirety using the same chunksizes")
    from ..multi import concat
    return concat([_read_single_hdf(path, key, start=start, stop=stop,
                                    columns=columns, chunksize=chunksize,
                                    sorted_index=sorted_index,
                                    lock=lock, mode=mode)
                   for path in paths])


if PY3:
    from ..core import _Frame
    _Frame.to_hdf.__doc__ = to_hdf.__doc__
