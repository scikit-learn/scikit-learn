import numpy as np
import pandas as pd

from ... import delayed
from ...compatibility import string_types
from .io import from_delayed, from_pandas


def read_sql_table(table, uri, index_col, divisions=None, npartitions=None,
                   limits=None, columns=None, bytes_per_chunk=256 * 2**20,
                   head_rows=5, schema=None, meta=None, **kwargs):
    """
    Create dataframe from an SQL table.

    If neither divisions or npartitions is given, the memory footprint of the
    first few rows will be determined, and partitions of size ~256MB will
    be used.

    Parameters
    ----------
    table : string or sqlalchemy expression
        Select columns from here.
    uri : string
        Full sqlalchemy URI for the database connection
    index_col : string
        Column which becomes the index, and defines the partitioning. Should
        be a indexed column in the SQL server, and any orderable type. If the
        type is number or time, then partition boundaries can be inferred from
        npartitions or bytes_per_chunk; otherwide must supply explicit
        ``divisions=``.
        ``index_col`` could be a function to return a value, e.g.,
        ``sql.func.abs(sql.column('value')).label('abs(value)')``.
        Labeling columns created by functions or arithmetic operations is
        required.
    divisions: sequence
        Values of the index column to split the table by. If given, this will
        override npartitions and bytes_per_chunk. The divisions are the value
        boundaries of the index column used to define the partitions. For
        example, ``divisions=list('acegikmoqsuwz')`` could be used to partition
        a string column lexographically into 12 partitions, with the implicit
        assumption that each partition contains similar numbers of records.
    npartitions : int
        Number of partitions, if divisions is not given. Will split the values
        of the index column linearly between limits, if given, or the column
        max/min. The index column must be numeric or time for this to work
    limits: 2-tuple or None
        Manually give upper and lower range of values for use with npartitions;
        if None, first fetches max/min from the DB. Upper limit, if
        given, is inclusive.
    columns : list of strings or None
        Which columns to select; if None, gets all; can include sqlalchemy
        functions, e.g.,
        ``sql.func.abs(sql.column('value')).label('abs(value)')``.
        Labeling columns created by functions or arithmetic operations is
        recommended.
    bytes_per_chunk: int
        If both divisions and npartitions is None, this is the target size of
        each partition, in bytes
    head_rows: int
        How many rows to load for inferring the data-types, unless passing meta
    meta: empty DataFrame or None
        If provided, do not attempt to infer dtypes, but use these, coercing
        all chunks on load
    schema: str or None
        If using a table name, pass this to sqlalchemy to select which DB
        schema to use within the URI connection
    kwargs : dict
        Additional parameters to pass to `pd.read_sql()`

    Returns
    -------
    dask.dataframe

    Examples
    --------
    >>> df = dd.read_sql('accounts', 'sqlite:///path/to/bank.db',
    ...                  npartitions=10, index_col='id')  # doctest: +SKIP
    """
    import sqlalchemy as sa
    from sqlalchemy import sql
    from sqlalchemy.sql import elements
    if index_col is None:
        raise ValueError("Must specify index column to partition on")
    engine = sa.create_engine(uri)
    m = sa.MetaData()
    if isinstance(table, string_types):
        table = sa.Table(table, m, autoload=True, autoload_with=engine,
                         schema=schema)

    index = (table.columns[index_col] if isinstance(index_col, string_types)
             else index_col)
    if not isinstance(index_col, string_types + (elements.Label,)):
        raise ValueError('Use label when passing an SQLAlchemy instance'
                         ' as the index (%s)' % index)
    if divisions and npartitions:
        raise TypeError('Must supply either divisions or npartitions, not both')

    columns = ([(table.columns[c] if isinstance(c, string_types) else c)
                for c in columns]
               if columns else list(table.columns))
    if index_col not in columns:
        columns.append(table.columns[index_col]
                       if isinstance(index_col, string_types)
                       else index_col)

    if isinstance(index_col, string_types):
        kwargs['index_col'] = index_col
    else:
        # function names get pandas auto-named
        kwargs['index_col'] = index_col.name

    if meta is None:
        # derrive metadata from first few rows
        q = sql.select(columns).limit(head_rows).select_from(table)
        head = pd.read_sql(q, engine, **kwargs)

        if head.empty:
            # no results at all
            name = table.name
            head = pd.read_sql_table(name, uri, index_col=index_col)
            return from_pandas(head, npartitions=1)

        bytes_per_row = (head.memory_usage(deep=True, index=True)).sum() / 5
        meta = head[:0]
    else:
        if divisions is None and npartitions is None:
            raise ValueError('Must provide divisions or npartitions when'
                             'using explicit meta.')

    if divisions is None:
        if limits is None:
            # calculate max and min for given index
            q = sql.select([sql.func.max(index), sql.func.min(index)]
                           ).select_from(table)
            minmax = pd.read_sql(q, engine)
            maxi, mini = minmax.iloc[0]
            dtype = minmax.dtypes['max_1']
        else:
            mini, maxi = limits
            dtype = pd.Series(limits).dtype
        if npartitions is None:
            q = sql.select([sql.func.count(index)]).select_from(table)
            count = pd.read_sql(q, engine)['count_1'][0]
            npartitions = round(count * bytes_per_row / bytes_per_chunk) or 1
        if dtype.kind == "M":
            divisions = pd.date_range(
                start=mini, end=maxi, freq='%iS' % (
                    (maxi - mini) / npartitions).total_seconds()).tolist()
            divisions[0] = mini
            divisions[-1] = maxi
        else:
            divisions = np.linspace(mini, maxi, npartitions + 1).tolist()

    parts = []
    lowers, uppers = divisions[:-1], divisions[1:]
    for i, (lower, upper) in enumerate(zip(lowers, uppers)):
        cond = index <= upper if i == len(lowers) - 1 else index < upper
        q = sql.select(columns).where(sql.and_(index >= lower, cond)
                                      ).select_from(table)
        parts.append(delayed(_read_sql_chunk)(q, uri, meta, **kwargs))

    return from_delayed(parts, meta, divisions=divisions)


def _read_sql_chunk(q, uri, meta, **kwargs):
    df = pd.read_sql(q, uri, **kwargs)
    if df.empty:
        return meta
    else:
        return df.astype(meta.dtypes.to_dict(), copy=False)
