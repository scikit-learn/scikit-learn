from __future__ import absolute_import, division, print_function

import re
import copy
import json
import warnings
import os
from distutils.version import LooseVersion

import numpy as np
import pandas as pd

from ..core import DataFrame, Series
from ..utils import (clear_known_categories, strip_unknown_categories,
                     UNKNOWN_CATEGORIES)
from ...bytes.compression import compress
from ...base import tokenize
from ...compatibility import PY3, string_types
from ...delayed import delayed
from ...bytes.core import get_fs_token_paths
from ...bytes.utils import infer_storage_options
from ...utils import import_required, natural_sort_key
from .utils import _get_pyarrow_dtypes, _meta_from_dtypes

__all__ = ('read_parquet', 'to_parquet')


def _parse_pandas_metadata(pandas_metadata):
    """Get the set of names from the pandas metadata section

    Parameters
    ----------
    pandas_metadata : dict
        Should conform to the pandas parquet metadata spec

    Returns
    -------
    index_names : list
        List of strings indicating the actual index names
    column_names : list
        List of strings indicating the actual column names
    storage_name_mapping : dict
        Pairs of storage names (e.g. the field names for
        PyArrow) and actual names. The storage and field names will
        differ for index names for certain writers (pyarrow > 0.8).
    column_indexes_names : list
        The names for ``df.columns.name`` or ``df.columns.names`` for
        a MultiIndex in the columns

    Notes
    -----
    This should support metadata written by at least

    * fastparquet>=0.1.3
    * pyarrow>=0.7.0
    """
    index_storage_names = pandas_metadata['index_columns']
    index_name_xpr = re.compile('__index_level_\d+__')

    # older metadatas will not have a 'field_name' field so we fall back
    # to the 'name' field
    pairs = [(x.get('field_name', x['name']), x['name'])
             for x in pandas_metadata['columns']]

    # Need to reconcile storage and real names. These will differ for
    # pyarrow, which uses __index_leveL_d__ for the storage name of indexes.
    # The real name may be None (e.g. `df.index.name` is None).
    pairs2 = []
    for storage_name, real_name in pairs:
        if real_name and index_name_xpr.match(real_name):
            real_name = None
        pairs2.append((storage_name, real_name))
    index_names = [name for (storage_name, name) in pairs2
                   if name != storage_name]

    # column_indexes represents df.columns.name
    # It was added to the spec after pandas 0.21.0+, and implemented
    # in PyArrow 0.8. It's not currently impelmented in fastparquet.
    column_index_names = pandas_metadata.get("column_indexes", [{'name': None}])
    column_index_names = [x['name'] for x in column_index_names]

    # Now we need to disambiguate between columns and index names. PyArrow
    # 0.8.0+ allows for duplicates between df.index.names and df.columns
    if not index_names:
        # For PyArrow < 0.8, Any fastparquet. This relies on the facts that
        # 1. Those versions used the real index name as the index storage name
        # 2. Those versions did not allow for duplicate index / column names
        # So we know that if a name is in index_storage_names, it must be an
        # index name
        index_names = list(index_storage_names)  # make a copy
        index_storage_names2 = set(index_storage_names)
        column_names = [name for (storage_name, name)
                        in pairs if name not in index_storage_names2]
    else:
        # For newer PyArrows the storage names differ from the index names
        # iff it's an index level. Though this is a fragile assumption for
        # other systems...
        column_names = [name for (storage_name, name) in pairs2
                        if name == storage_name]

    storage_name_mapping = dict(pairs2)   # TODO: handle duplicates gracefully

    return index_names, column_names, storage_name_mapping, column_index_names


def _normalize_index_columns(user_columns, data_columns, user_index, data_index):
    """Normalize user and file-provided column and index names

    Parameters
    ----------
    user_columns : None, str or list of str
    data_columns : list of str
    user_index : None, str, or list of str
    data_index : list of str

    Returns
    -------
    column_names : list of str
    index_names : list of str
    out_type : {pd.Series, pd.DataFrame}
    """
    specified_columns = user_columns is not None
    specified_index = user_index is not None

    out_type = DataFrame

    if user_columns is None:
        user_columns = list(data_columns)
    elif isinstance(user_columns, string_types):
        user_columns = [user_columns]
        out_type = Series
    else:
        user_columns = list(user_columns)

    if user_index is None:
        user_index = data_index
    elif user_index is False:
        # When index is False, use no index and all fields should be treated as
        # columns (unless `columns` provided).
        user_index = []
        data_columns = data_index + data_columns
    elif isinstance(user_index, string_types):
        user_index = [user_index]
    else:
        user_index = list(user_index)

    if specified_index and not specified_columns:
        # Only `index` provided. Use specified index, and all column fields
        # that weren't specified as indices
        index_names = user_index
        column_names = [x for x in data_columns if x not in index_names]
    elif specified_columns and not specified_index:
        # Only `columns` provided. Use specified columns, and all index fields
        # that weren't specified as columns
        column_names = user_columns
        index_names = [x for x in data_index if x not in column_names]
    elif specified_index and specified_columns:
        # Both `index` and `columns` provided. Use as specified, but error if
        # they intersect.
        column_names = user_columns
        index_names = user_index
        if set(column_names).intersection(index_names):
            raise ValueError("Specified index and column names must not "
                             "intersect")
    else:
        # Use default columns and index from the metadata
        column_names = data_columns
        index_names = data_index

    return column_names, index_names, out_type


# ----------------------------------------------------------------------
# Fastparquet interface


def _read_fastparquet(fs, fs_token, paths, columns=None, filters=None,
                      categories=None, index=None, infer_divisions=None):
    import fastparquet
    from fastparquet.util import check_column_names

    if len(paths) > 1:
        pf = fastparquet.ParquetFile(paths, open_with=fs.open, sep=fs.sep)
    else:
        try:
            pf = fastparquet.ParquetFile(paths[0] + fs.sep + '_metadata',
                                         open_with=fs.open,
                                         sep=fs.sep)
        except Exception:
            pf = fastparquet.ParquetFile(paths[0], open_with=fs.open, sep=fs.sep)

    # Validate infer_divisions
    if os.path.split(pf.fn)[-1] != '_metadata' and infer_divisions is True:
        raise NotImplementedError("infer_divisions=True is not supported by the fastparquet engine for datasets "
                                  "that do not contain a global '_metadata' file")

    check_column_names(pf.columns, categories)
    if isinstance(columns, tuple):
        # ensure they tokenize the same
        columns = list(columns)

    if pf.fmd.key_value_metadata:
        pandas_md = [x.value for x in pf.fmd.key_value_metadata if x.key == 'pandas']
    else:
        pandas_md = []

    if len(pandas_md) == 0:
        # Fall back to the storage information
        index_names = pf._get_index()
        if not isinstance(index_names, list):
            index_names = [index_names]
        column_names = pf.columns + list(pf.cats)
        storage_name_mapping = {k: k for k in column_names}
    elif len(pandas_md) == 1:
        index_names, column_names, storage_name_mapping, column_index_names = (
            _parse_pandas_metadata(json.loads(pandas_md[0]))
        )
        column_names.extend(pf.cats)
    else:
        raise ValueError("File has multiple entries for 'pandas' metadata")

    # Normalize user inputs

    if filters is None:
        filters = []

    column_names, index_names, out_type = _normalize_index_columns(
        columns, column_names, index, index_names)

    if categories is None:
        categories = pf.categories
    elif isinstance(categories, string_types):
        categories = [categories]
    else:
        categories = list(categories)

    # TODO: write partition_on to pandas metadata...
    all_columns = list(column_names)
    all_columns.extend(x for x in index_names if x not in column_names)

    rgs = [rg for rg in pf.row_groups if
           not (fastparquet.api.filter_out_stats(rg, filters, pf.schema)) and
           not (fastparquet.api.filter_out_cats(rg, filters))]

    dtypes = pf._dtypes(categories)
    dtypes = {storage_name_mapping.get(k, k): v for k, v in dtypes.items()}

    meta = _meta_from_dtypes(all_columns, dtypes, index_names, [None])
    # fastparquet doesn't handle multiindex
    if len(index_names) > 1:
        raise ValueError("Cannot read DataFrame with MultiIndex.")
    elif len(index_names) == 0:
        index_names = None

    for cat in categories:
        if cat in meta:
            meta[cat] = pd.Series(pd.Categorical([],
                                                 categories=[UNKNOWN_CATEGORIES]),
                                  index=meta.index)

    for catcol in pf.cats:
        if catcol in meta.columns:
            meta[catcol] = meta[catcol].cat.set_categories(pf.cats[catcol])
        elif meta.index.name == catcol:
            meta.index = meta.index.set_categories(pf.cats[catcol])

    if out_type == Series:
        assert len(meta.columns) == 1
        meta = meta[meta.columns[0]]

    name = 'read-parquet-' + tokenize(fs_token, paths, all_columns, filters,
                                      categories)

    dsk = {(name, i): (_read_parquet_row_group, fs, pf.row_group_filename(rg),
                       index_names, all_columns, rg, out_type == Series,
                       categories, pf.schema, pf.cats, pf.dtypes,
                       pf.file_scheme, storage_name_mapping)
           for i, rg in enumerate(rgs)}
    if not dsk:
        # empty dataframe
        dsk = {(name, 0): meta}
        divisions = (None, None)
        return out_type(dsk, name, meta, divisions)

    if index_names and infer_divisions is not False:
        index_name = meta.index.name
        minmax = fastparquet.api.sorted_partitioned_columns(pf)
        if index_name in minmax:
            divisions = (list(minmax[index_name]['min']) +
                         [minmax[index_name]['max'][-1]])
            divisions = [divisions[i] for i, rg in enumerate(pf.row_groups)
                         if rg in rgs] + [divisions[-1]]
        else:
            if infer_divisions is True:
                raise ValueError(
                    ("Unable to infer divisions for index of '{index_name}' because it is not known to be "
                     "sorted across partitions").format(index_name=index_name))

            divisions = (None,) * (len(rgs) + 1)
    else:
        if infer_divisions is True:
            raise ValueError(
                'Unable to infer divisions for because no index column was discovered')

        divisions = (None,) * (len(rgs) + 1)

    if isinstance(divisions[0], np.datetime64):
        divisions = [pd.Timestamp(d) for d in divisions]

    return out_type(dsk, name, meta, divisions)


def _read_parquet_row_group(fs, fn, index, columns, rg, series, categories,
                            schema, cs, dt, scheme, storage_name_mapping, *args):
    from fastparquet.api import _pre_allocate
    from fastparquet.core import read_row_group_file
    from collections import OrderedDict

    name_storage_mapping = {v: k for k, v in storage_name_mapping.items()}
    if not isinstance(columns, (tuple, list)):
        columns = [columns,]
        series = True
    if index:
        index, = index
        if index not in columns:
            columns = columns + [index]

    columns = [name_storage_mapping.get(col, col) for col in columns]
    index = name_storage_mapping.get(index, index)
    cs = OrderedDict([(k, v) for k, v in cs.items() if k in columns])

    df, views = _pre_allocate(rg.num_rows, columns, categories, index, cs, dt)
    read_row_group_file(fn, rg, columns, categories, schema, cs,
                        open=fs.open, assign=views, scheme=scheme)

    if df.index.nlevels == 1:
        if index:
            df.index.name = storage_name_mapping.get(index, index)
    else:
        if index:
            df.index.names = [storage_name_mapping.get(name, name)
                              for name in index]
    df.columns = [storage_name_mapping.get(col, col)
                  for col in columns
                  if col != index]

    if series:
        return df[df.columns[0]]
    else:
        return df


def _write_partition_fastparquet(df, fs, path, filename, fmd, compression,
                                 partition_on):
    from fastparquet.writer import partition_on_columns, make_part_file
    import fastparquet
    # Fastparquet mutates this in a non-threadsafe manner. For now we just copy
    # it before forwarding to fastparquet.
    fmd = copy.copy(fmd)
    if not len(df):
        # Write nothing for empty partitions
        rgs = None
    elif partition_on:
        if LooseVersion(fastparquet.__version__) >= '0.1.4':
            rgs = partition_on_columns(df, partition_on, path, filename, fmd,
                                       compression, fs.open, fs.mkdirs)
        else:
            rgs = partition_on_columns(df, partition_on, path, filename, fmd,
                                       fs.sep, compression, fs.open, fs.mkdirs)
    else:
        # Fastparquet current doesn't properly set `num_rows` in the output
        # metadata. Set it here to fix that.
        fmd.num_rows = len(df)
        with fs.open(fs.sep.join([path, filename]), 'wb') as fil:
            rgs = make_part_file(fil, df, fmd.schema, compression=compression,
                                 fmd=fmd)
    return rgs


def _write_fastparquet(df, fs, fs_token, path, write_index=None, append=False,
                       ignore_divisions=False, partition_on=None,
                       compression=None, **kwargs):
    import fastparquet

    fs.mkdirs(path)
    sep = fs.sep

    object_encoding = kwargs.pop('object_encoding', 'utf8')
    if object_encoding == 'infer' or (isinstance(object_encoding, dict) and 'infer' in object_encoding.values()):
        raise ValueError('"infer" not allowed as object encoding, '
                         'because this required data in memory.')

    divisions = df.divisions
    if write_index is True or write_index is None and df.known_divisions:
        df = df.reset_index()
        index_cols = [df.columns[0]]
    else:
        ignore_divisions = True
        index_cols = []

    if append:
        try:
            pf = fastparquet.api.ParquetFile(path, open_with=fs.open, sep=sep)
        except (IOError, ValueError):
            # append for create
            append = False
    if append:
        if pf.file_scheme not in ['hive', 'empty', 'flat']:
            raise ValueError('Requested file scheme is hive, '
                             'but existing file scheme is not.')
        elif ((set(pf.columns) != set(df.columns) - set(partition_on)) or (set(partition_on) != set(pf.cats))):
            raise ValueError('Appended columns not the same.\n'
                             'New: {} | Previous: {}'
                             .format(pf.columns, list(df.columns)))
        elif set(pf.dtypes[c] for c in pf.columns) != set(df[pf.columns].dtypes):
            raise ValueError('Appended dtypes differ.\n{}'
                             .format(set(pf.dtypes.items()) ^
                                     set(df.dtypes.iteritems())))
        else:
            df = df[pf.columns + partition_on]

        fmd = pf.fmd
        i_offset = fastparquet.writer.find_max_part(fmd.row_groups)

        if not ignore_divisions:
            minmax = fastparquet.api.sorted_partitioned_columns(pf)
            old_end = minmax[index_cols[0]]['max'][-1]
            if divisions[0] < old_end:
                raise ValueError(
                    'Appended divisions overlapping with the previous ones.\n'
                    'New: {} | Previous: {}'.format(old_end, divisions[0]))
    else:
        fmd = fastparquet.writer.make_metadata(df._meta,
                                               object_encoding=object_encoding,
                                               index_cols=index_cols,
                                               ignore_columns=partition_on,
                                               **kwargs)
        i_offset = 0

    filenames = ['part.%i.parquet' % (i + i_offset)
                 for i in range(df.npartitions)]

    write = delayed(_write_partition_fastparquet, pure=False)
    writes = [write(part, fs, path, filename, fmd, compression, partition_on)
              for filename, part in zip(filenames, df.to_delayed())]

    return delayed(_write_metadata)(writes, filenames, fmd, path, fs, sep)


def _write_metadata(writes, filenames, fmd, path, fs, sep):
    """ Write Parquet metadata after writing all row groups

    See Also
    --------
    to_parquet
    """
    import fastparquet
    fmd = copy.copy(fmd)
    for fn, rg in zip(filenames, writes):
        if rg is not None:
            if isinstance(rg, list):
                for r in rg:
                    fmd.row_groups.append(r)
            else:
                for chunk in rg.columns:
                    chunk.file_path = fn
                fmd.row_groups.append(rg)

    fn = sep.join([path, '_metadata'])
    fastparquet.writer.write_common_metadata(fn, fmd, open_with=fs.open,
                                             no_row_groups=False)

    fn = sep.join([path, '_common_metadata'])
    fastparquet.writer.write_common_metadata(fn, fmd, open_with=fs.open)


# ----------------------------------------------------------------------
# PyArrow interface

def _read_pyarrow(fs, fs_token, paths, columns=None, filters=None,
                  categories=None, index=None, infer_divisions=None):
    from ...bytes.core import get_pyarrow_filesystem
    import pyarrow.parquet as pq
    import pyarrow as pa

    # In pyarrow, the physical storage field names may differ from
    # the actual dataframe names. This is true for Index names when
    # PyArrow >= 0.8.
    # We would like to resolve these to the correct dataframe names
    # as soon as possible.

    if filters is not None:
        raise NotImplementedError("Predicate pushdown not implemented")

    if isinstance(categories, string_types):
        categories = [categories]
    elif categories is None:
        categories = []
    else:
        categories = list(categories)

    if isinstance(columns, tuple):
        columns = list(columns)

    dataset = pq.ParquetDataset(paths, filesystem=get_pyarrow_filesystem(fs))
    schema = dataset.schema.to_arrow_schema()
    has_pandas_metadata = schema.metadata is not None and b'pandas' in schema.metadata

    if has_pandas_metadata:
        pandas_metadata = json.loads(schema.metadata[b'pandas'].decode('utf8'))
        index_names, column_names, storage_name_mapping, column_index_names = (
            _parse_pandas_metadata(pandas_metadata)
        )
    else:
        index_names = []
        column_names = schema.names
        storage_name_mapping = {k: k for k in column_names}
        column_index_names = [None]

    if pa.__version__ < LooseVersion('0.8.0'):
        # the pyarrow 0.7.0 *reader* expects the storage names for index names
        # that are None.
        if any(x is None for x in index_names):
            name_storage_mapping = {v: k for
                                    k, v in storage_name_mapping.items()}
            index_names = [name_storage_mapping.get(name, name)
                           for name in index_names]

    column_names, index_names, out_type = _normalize_index_columns(
        columns, column_names, index, index_names)

    all_columns = index_names + column_names

    # Find non-empty pieces
    non_empty_pieces = []
    # Determine valid pieces
    _open = lambda fn: pq.ParquetFile(fs.open(fn, mode='rb'))
    for piece in dataset.pieces:
        pf = piece.get_metadata(_open)
        # non_empty_pieces.append(piece)
        if pf.num_row_groups > 0:
            non_empty_pieces.append(piece)

    # Sort pieces naturally
    # If a single input path resulted in multiple dataset pieces, then sort
    # the pieces naturally. If multiple paths were supplied then we leave
    # the order of the resulting pieces unmodified
    if len(paths) == 1 and len(dataset.pieces) > 1:
        non_empty_pieces = sorted(
            non_empty_pieces, key=lambda piece: natural_sort_key(piece.path))

    # Determine divisions
    if len(index_names) == 1:

        # Look up storage name of the single index column
        divisions_names = [storage_name for storage_name, name
                           in storage_name_mapping.items()
                           if index_names[0] == name]

        if divisions_names:
            divisions_name = divisions_names[0]
        else:
            divisions_name = None
    else:
        divisions_name = None

    divisions = _get_pyarrow_divisions(non_empty_pieces, divisions_name,
                                       schema, infer_divisions)

    # Build task
    dtypes = _get_pyarrow_dtypes(schema, categories)
    dtypes = {storage_name_mapping.get(k, k): v for k, v in dtypes.items()}

    meta = _meta_from_dtypes(all_columns, dtypes, index_names,
                             column_index_names)
    meta = clear_known_categories(meta, cols=categories)

    if out_type == Series:
        assert len(meta.columns) == 1
        meta = meta[meta.columns[0]]

    task_name = 'read-parquet-' + tokenize(fs_token, paths, all_columns)

    if non_empty_pieces:
        task_plan = {
            (task_name, i): (_read_pyarrow_parquet_piece,
                             fs,
                             piece,
                             column_names,
                             index_names,
                             out_type == Series,
                             dataset.partitions,
                             categories)
            for i, piece in enumerate(non_empty_pieces)
        }
    else:
        meta = strip_unknown_categories(meta)
        task_plan = {(task_name, 0): meta}

    return out_type(task_plan, task_name, meta, divisions)


def _to_ns(val, unit):
    """
    Convert an input time in the specified units to nanoseconds

    Parameters
    ----------
    val: int
        Input time value
    unit : str
        Time units of `val`.
        One of 's', 'ms', 'us', 'ns'

    Returns
    -------
    int
        Time val in nanoseconds
    """
    factors = {'s': int(1e9), 'ms': int(1e6), 'us': int(1e3), 'ns': 1}
    try:
        factor = factors.get(unit)
    except KeyError:
        raise ValueError("Unsupported time unit '{unit}'".format(unit=unit))

    return val * factor


def _get_pyarrow_divisions(pa_pieces, divisions_name, pa_schema, infer_divisions):
    """
    Compute DataFrame divisions from a list of pyarrow dataset pieces

    Parameters
    ----------
    pa_pieces : list[pyarrow.parquet.ParquetDatasetPiece]
        List of dataset pieces. Each piece corresponds to a single partition in the eventual dask DataFrame
    divisions_name : str|None
        The name of the column to compute divisions for
    pa_schema : pyarrow.lib.Schema
        The pyarrow schema for the dataset
    infer_divisions : bool or None
        If True divisions must be inferred (otherwise an exception is raised). If False or None divisions are not
        inferred
    Returns
    -------
    list
    """
    # Local imports
    import pyarrow as pa
    import pyarrow.parquet as pq

    if infer_divisions is True and pa.__version__ < LooseVersion('0.9.0'):
        raise NotImplementedError('infer_divisions=True requires pyarrow >=0.9.0')

    # Check whether divisions_name is in the schema
    # Note: get_field_index returns -1 if not found, but it does not accept None
    if infer_divisions is True:
        divisions_name_in_schema = divisions_name is not None and pa_schema.get_field_index(divisions_name) >= 0

        if divisions_name_in_schema is False and infer_divisions is True:
            raise ValueError(
                'Unable to infer divisions for because no index column was discovered')
    else:
        divisions_name_in_schema = None

    if pa_pieces and divisions_name_in_schema:
        # We have pieces and a valid division column.
        # Compute min/max for column in each row group
        min_maxs = []
        last_max = None

        # Initialize index of divisions column within the row groups.
        # To be computed during while processing the first piece below
        for piece in pa_pieces:
            pf = piece.get_metadata(pq.ParquetFile)
            rg = pf.row_group(0)

            # Compute division column index
            rg_paths = [rg.column(i).path_in_schema for i in range(rg.num_columns)]
            try:
                divisions_col_index = rg_paths.index(divisions_name)
            except ValueError:
                # Divisions not valid
                min_maxs = None
                break

            col_meta = rg.column(divisions_col_index)
            stats = col_meta.statistics
            if stats.has_min_max and (last_max is None or last_max < stats.min):
                min_maxs.append((stats.min, stats.max))
                last_max = stats.max
            else:
                # Divisions not valid
                min_maxs = None
                break

        if min_maxs:
            # We have min/max pairs
            divisions = [mn for mn, mx in min_maxs] + [min_maxs[-1][1]]

            # Handle conversion to pandas timestamp divisions
            index_field = pa_schema.field_by_name(divisions_name)
            if isinstance(index_field.type, pa.TimestampType):
                time_unit = index_field.type.unit
                divisions_ns = [_to_ns(d, time_unit) for d in
                                divisions]
                divisions = [pd.Timestamp(ns) for ns in divisions_ns]

            # Handle encoding of bytes string
            if index_field.type == pa.string():
                # Parquet strings are always encoded as utf-8
                encoding = 'utf-8'
                divisions = [d.decode(encoding).strip() for d in divisions]

        else:
            if infer_divisions is True:
                raise ValueError(
                    ("Unable to infer divisions for index of '{index_name}' because it is not known to be "
                     "sorted across partitions").format(index_name=divisions_name_in_schema))

            divisions = (None,) * (len(pa_pieces) + 1)
    elif pa_pieces:
        divisions = (None,) * (len(pa_pieces) + 1)
    else:
        divisions = (None, None)
    return divisions


def _read_pyarrow_parquet_piece(fs, piece, columns, index_cols, is_series,
                                partitions, categories):
    import pyarrow as pa

    with fs.open(piece.path, mode='rb') as f:
        table = piece.read(columns=index_cols + columns,
                           partitions=partitions,
                           use_pandas_metadata=True,
                           file=f)

    if pa.__version__ < LooseVersion('0.9.0'):
        df = table.to_pandas()
        for cat in categories:
            df[cat] = df[cat].astype('category')
    else:
        df = table.to_pandas(categories=categories)
    has_index = not isinstance(df.index, pd.RangeIndex)

    if not has_index and index_cols:
        # Index should be set, but it isn't
        df = df.set_index(index_cols)
    elif has_index and df.index.names != index_cols:
        # Index is set, but isn't correct
        # This can happen when reading in not every column in a multi-index
        df = df.reset_index(drop=False)
        if index_cols:
            df = df.set_index(index_cols)
        drop = list(set(df.columns).difference(columns))
        if drop:
            df = df.drop(drop, axis=1)
        # Ensure proper ordering
        df = df.reindex(columns=columns, copy=False)

    if is_series:
        return df[df.columns[0]]
    else:
        return df


_pyarrow_write_table_kwargs = {'row_group_size', 'version', 'use_dictionary',
                               'compression', 'use_deprecated_int96_timestamps',
                               'coerce_timestamps', 'flavor', 'chunk_size'}

_pyarrow_write_metadata_kwargs = {'version', 'use_deprecated_int96_timestamps',
                                  'coerce_timestamps'}


def _write_pyarrow(df, fs, fs_token, path, write_index=None, append=False,
                   ignore_divisions=False, partition_on=None, **kwargs):
    if append:
        raise NotImplementedError("`append` not implemented for "
                                  "`engine='pyarrow'`")

    if ignore_divisions:
        raise NotImplementedError("`ignore_divisions` not implemented for "
                                  "`engine='pyarrow'`")

    # We can check only write_table kwargs, as it is a superset of kwargs for write functions
    if set(kwargs).difference(_pyarrow_write_table_kwargs):
        msg = ("Unexpected keyword arguments: " +
               "%r" % list(set(kwargs).difference(_pyarrow_write_table_kwargs)))
        raise TypeError(msg)

    if write_index is None and df.known_divisions:
        write_index = True

    fs.mkdirs(path)

    template = fs.sep.join([path, 'part.%i.parquet'])

    write = delayed(_write_partition_pyarrow, pure=False)
    first_kwargs = kwargs.copy()
    first_kwargs['metadata_path'] = fs.sep.join([path, '_common_metadata'])
    writes = [write(part, path, fs, template % i, write_index, partition_on,
                    **(kwargs if i else first_kwargs))
              for i, part in enumerate(df.to_delayed())]
    return delayed(writes)


def _write_partition_pyarrow(df, path, fs, filename, write_index,
                             partition_on, metadata_path=None, **kwargs):
    import pyarrow as pa
    from pyarrow import parquet
    t = pa.Table.from_pandas(df, preserve_index=write_index)

    if partition_on:
        parquet.write_to_dataset(t, path, partition_cols=partition_on,
                                 preserve_index=write_index,
                                 filesystem=fs, **kwargs)
    else:
        with fs.open(filename, 'wb') as fil:
            parquet.write_table(t, fil, **kwargs)

    if metadata_path is not None:
        with fs.open(metadata_path, 'wb') as fil:
            # Get only arguments specified in the function
            kwargs_meta = {k: v for k, v in kwargs.items()
                           if k in _pyarrow_write_metadata_kwargs}
            parquet.write_metadata(t.schema, fil, **kwargs_meta)


# ----------------------------------------------------------------------
# User API


_ENGINES = {}


def get_engine(engine):
    """Get the parquet engine backend implementation.

    Parameters
    ----------
    engine : {'auto', 'fastparquet', 'pyarrow'}, default 'auto'
        Parquet reader library to use. Default is first installed in this list.

    Returns
    -------
    A dict containing a ``'read'`` and ``'write'`` function.
    """
    if engine in _ENGINES:
        return _ENGINES[engine]

    if engine == 'auto':
        for eng in ['fastparquet', 'pyarrow']:
            try:
                return get_engine(eng)
            except RuntimeError:
                pass
        else:
            raise RuntimeError("Please install either fastparquet or pyarrow")

    elif engine == 'fastparquet':
        import_required('fastparquet', "`fastparquet` not installed")

        _ENGINES['fastparquet'] = eng = {'read': _read_fastparquet,
                                         'write': _write_fastparquet}
        return eng

    elif engine == 'pyarrow':
        import_required('pyarrow', "`pyarrow` not installed")

        _ENGINES['pyarrow'] = eng = {'read': _read_pyarrow,
                                     'write': _write_pyarrow}
        return eng

    elif engine == 'arrow':
        warnings.warn("parquet with `engine='arrow'` is deprecated, "
                      "use `engine='pyarrow'` instead")
        return get_engine('pyarrow')

    else:
        raise ValueError('Unsupported engine type: {0}'.format(engine))


def read_parquet(path, columns=None, filters=None, categories=None, index=None,
                 storage_options=None, engine='auto', infer_divisions=None):
    """
    Read ParquetFile into a Dask DataFrame

    This reads a directory of Parquet data into a Dask.dataframe, one file per
    partition.  It selects the index among the sorted columns if any exist.

    Parameters
    ----------
    path : string or list
        Source directory for data, or path(s) to individual parquet files.
        Prefix with a protocol like ``s3://`` to read from alternative
        filesystems. To read from multiple files you can pass a globstring or a
        list of paths, with the caveat that they must all have the same
        protocol.
    columns : string, list or None (default)
        Field name(s) to read in as columns in the output. By default all
        non-index fields will be read (as determined by the pandas parquet
        metadata, if present). Provide a single field name instead of a list to
        read in the data as a Series.
    filters : list
        List of filters to apply, like ``[('x', '>', 0), ...]``. This implements
        row-group (partition) -level filtering only, i.e., to prevent the
        loading of some chunks of the data, and only if relevant statistics
        have been included in the metadata.
    index : string, list, False or None (default)
        Field name(s) to use as the output frame index. By default will be
        inferred from the pandas parquet file metadata (if present). Use False
        to read all fields as columns.
    categories : list, dict or None
        For any fields listed here, if the parquet encoding is Dictionary,
        the column will be created with dtype category. Use only if it is
        guaranteed that the column is encoded as dictionary in all row-groups.
        If a list, assumes up to 2**16-1 labels; if a dict, specify the number
        of labels expected; if None, will load categories automatically for
        data written by dask/fastparquet, not otherwise.
    storage_options : dict
        Key/value pairs to be passed on to the file-system backend, if any.
    engine : {'auto', 'fastparquet', 'pyarrow'}, default 'auto'
        Parquet reader library to use. If only one library is installed, it
        will use that one; if both, it will use 'fastparquet'
    infer_divisions : bool or None (default).
        By default, divisions are inferred if the read `engine` supports
        doing so efficiently and the `index` of the underlying dataset is
        sorted across the individual parquet files. Set to ``True`` to
        force divisions to be inferred in all cases. Note that this may
        require reading metadata from each file in the dataset, which may
        be expensive. Set to ``False`` to never infer divisions.

    Examples
    --------
    >>> df = dd.read_parquet('s3://bucket/my-parquet-data')  # doctest: +SKIP

    See Also
    --------
    to_parquet
    """
    read = get_engine(engine)['read']

    fs, fs_token, paths = get_fs_token_paths(path, mode='rb',
                                             storage_options=storage_options)

    if isinstance(path, string_types) and len(paths) > 1:
        # Sort paths naturally if multiple paths resulted from a single
        # specification (by '*' globbing)
        paths = sorted(paths, key=natural_sort_key)

    return read(fs, fs_token, paths, columns=columns, filters=filters,
                categories=categories, index=index, infer_divisions=infer_divisions)


def to_parquet(df, path, engine='auto', compression='default', write_index=None,
               append=False, ignore_divisions=False, partition_on=None,
               storage_options=None, compute=True, **kwargs):
    """Store Dask.dataframe to Parquet files

    Notes
    -----
    Each partition will be written to a separate file.

    Parameters
    ----------
    df : dask.dataframe.DataFrame
    path : string
        Destination directory for data.  Prepend with protocol like ``s3://``
        or ``hdfs://`` for remote data.
    engine : {'auto', 'fastparquet', 'pyarrow'}, default 'auto'
        Parquet library to use. If only one library is installed, it will use
        that one; if both, it will use 'fastparquet'.
    compression : string or dict, optional
        Either a string like ``"snappy"`` or a dictionary mapping column names
        to compressors like ``{"name": "gzip", "values": "snappy"}``. The
        default is ``"default"``, which uses the default compression for
        whichever engine is selected.
    write_index : boolean, optional
        Whether or not to write the index. Defaults to True *if* divisions are
        known.
    append : bool, optional
        If False (default), construct data-set from scratch. If True, add new
        row-group(s) to an existing data-set. In the latter case, the data-set
        must exist, and the schema must match the input data.
    ignore_divisions : bool, optional
        If False (default) raises error when previous divisions overlap with
        the new appended divisions. Ignored if append=False.
    partition_on : list, optional
        Construct directory-based partitioning by splitting on these fields'
        values. Each dask partition will result in one or more datafiles,
        there will be no global groupby.
    storage_options : dict, optional
        Key/value pairs to be passed on to the file-system backend, if any.
    compute : bool, optional
        If True (default) then the result is computed immediately. If False
        then a ``dask.delayed`` object is returned for future computation.
    **kwargs
        Extra options to be passed on to the specific backend.

    Examples
    --------
    >>> df = dd.read_csv(...)  # doctest: +SKIP
    >>> to_parquet('/path/to/output/', df, compression='snappy')  # doctest: +SKIP

    See Also
    --------
    read_parquet: Read parquet data to dask.dataframe
    """
    # TODO: remove once deprecation cycle is finished
    if isinstance(path, DataFrame):
        warnings.warn("DeprecationWarning: The order of `df` and `path` in "
                      "`dd.to_parquet` has switched, please update your code")
        df, path = path, df

    partition_on = partition_on or []

    if set(partition_on) - set(df.columns):
        raise ValueError('Partitioning on non-existent column')

    if compression != 'default':
        kwargs['compression'] = compression
    elif 'snappy' in compress:
        kwargs['compression'] = 'snappy'

    write = get_engine(engine)['write']

    fs, fs_token, _ = get_fs_token_paths(path, mode='wb',
                                         storage_options=storage_options)
    # Trim any protocol information from the path before forwarding
    path = infer_storage_options(path)['path']

    out = write(df, fs, fs_token, path, write_index=write_index, append=append,
                ignore_divisions=ignore_divisions, partition_on=partition_on,
                **kwargs)

    if compute:
        out.compute()
        return None
    return out


if PY3:
    DataFrame.to_parquet.__doc__ = to_parquet.__doc__
