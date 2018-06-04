"""
Table Schema builders

http://specs.frictionlessdata.io/json-table-schema/
"""
import warnings

import pandas._libs.json as json
from pandas import DataFrame
from pandas.api.types import CategoricalDtype
import pandas.core.common as com
from pandas.core.dtypes.common import (
    is_integer_dtype, is_timedelta64_dtype, is_numeric_dtype,
    is_bool_dtype, is_datetime64_dtype, is_datetime64tz_dtype,
    is_categorical_dtype, is_period_dtype, is_string_dtype
)

loads = json.loads


def as_json_table_type(x):
    """
    Convert a NumPy / pandas type to its corresponding json_table.

    Parameters
    ----------
    x : array or dtype

    Returns
    -------
    t : str
        the Table Schema data types

    Notes
    -----
    This table shows the relationship between NumPy / pandas dtypes,
    and Table Schema dtypes.

    ==============  =================
    Pandas type     Table Schema type
    ==============  =================
    int64           integer
    float64         number
    bool            boolean
    datetime64[ns]  datetime
    timedelta64[ns] duration
    object          str
    categorical     any
    =============== =================
    """
    if is_integer_dtype(x):
        return 'integer'
    elif is_bool_dtype(x):
        return 'boolean'
    elif is_numeric_dtype(x):
        return 'number'
    elif (is_datetime64_dtype(x) or is_datetime64tz_dtype(x) or
          is_period_dtype(x)):
        return 'datetime'
    elif is_timedelta64_dtype(x):
        return 'duration'
    elif is_categorical_dtype(x):
        return 'any'
    elif is_string_dtype(x):
        return 'string'
    else:
        return 'any'


def set_default_names(data):
    """Sets index names to 'index' for regular, or 'level_x' for Multi"""
    if com._all_not_none(*data.index.names):
        nms = data.index.names
        if len(nms) == 1 and data.index.name == 'index':
            warnings.warn("Index name of 'index' is not round-trippable")
        elif len(nms) > 1 and any(x.startswith('level_') for x in nms):
            warnings.warn("Index names beginning with 'level_' are not "
                          "round-trippable")
        return data

    data = data.copy()
    if data.index.nlevels > 1:
        names = [name if name is not None else 'level_{}'.format(i)
                 for i, name in enumerate(data.index.names)]
        data.index.names = names
    else:
        data.index.name = data.index.name or 'index'
    return data


def convert_pandas_type_to_json_field(arr, dtype=None):
    dtype = dtype or arr.dtype
    if arr.name is None:
        name = 'values'
    else:
        name = arr.name
    field = {'name': name,
             'type': as_json_table_type(dtype)}

    if is_categorical_dtype(arr):
        if hasattr(arr, 'categories'):
            cats = arr.categories
            ordered = arr.ordered
        else:
            cats = arr.cat.categories
            ordered = arr.cat.ordered
        field['constraints'] = {"enum": list(cats)}
        field['ordered'] = ordered
    elif is_period_dtype(arr):
        field['freq'] = arr.freqstr
    elif is_datetime64tz_dtype(arr):
        if hasattr(arr, 'dt'):
            field['tz'] = arr.dt.tz.zone
        else:
            field['tz'] = arr.tz.zone
    return field


def convert_json_field_to_pandas_type(field):
    """
    Converts a JSON field descriptor into its corresponding NumPy / pandas type

    Parameters
    ----------
    field
        A JSON field descriptor

    Returns
    -------
    dtype

    Raises
    -----
    ValueError
        If the type of the provided field is unknown or currently unsupported

    Examples
    --------
    >>> convert_json_field_to_pandas_type({'name': 'an_int',
                                           'type': 'integer'})
    'int64'
    >>> convert_json_field_to_pandas_type({'name': 'a_categorical',
                                           'type': 'any',
                                           'contraints': {'enum': [
                                                          'a', 'b', 'c']},
                                           'ordered': True})
    'CategoricalDtype(categories=['a', 'b', 'c'], ordered=True)'
    >>> convert_json_field_to_pandas_type({'name': 'a_datetime',
                                           'type': 'datetime'})
    'datetime64[ns]'
    >>> convert_json_field_to_pandas_type({'name': 'a_datetime_with_tz',
                                           'type': 'datetime',
                                           'tz': 'US/Central'})
    'datetime64[ns, US/Central]'
    """
    typ = field['type']
    if typ == 'string':
        return 'object'
    elif typ == 'integer':
        return 'int64'
    elif typ == 'number':
        return 'float64'
    elif typ == 'boolean':
        return 'bool'
    elif typ == 'duration':
        return 'timedelta64'
    elif typ == 'datetime':
        if field.get('tz'):
            return 'datetime64[ns, {tz}]'.format(tz=field['tz'])
        else:
            return 'datetime64[ns]'
    elif typ == 'any':
        if 'constraints' in field and 'ordered' in field:
            return CategoricalDtype(categories=field['constraints']['enum'],
                                    ordered=field['ordered'])
        else:
            return 'object'

    raise ValueError("Unsupported or invalid field type: {}".format(typ))


def build_table_schema(data, index=True, primary_key=None, version=True):
    """
    Create a Table schema from ``data``.

    Parameters
    ----------
    data : Series, DataFrame
    index : bool, default True
        Whether to include ``data.index`` in the schema.
    primary_key : bool or None, default True
        column names to designate as the primary key.
        The default `None` will set `'primaryKey'` to the index
        level or levels if the index is unique.
    version : bool, default True
        Whether to include a field `pandas_version` with the version
        of pandas that generated the schema.

    Returns
    -------
    schema : dict

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {'A': [1, 2, 3],
    ...      'B': ['a', 'b', 'c'],
    ...      'C': pd.date_range('2016-01-01', freq='d', periods=3),
    ...     }, index=pd.Index(range(3), name='idx'))
    >>> build_table_schema(df)
    {'fields': [{'name': 'idx', 'type': 'integer'},
    {'name': 'A', 'type': 'integer'},
    {'name': 'B', 'type': 'string'},
    {'name': 'C', 'type': 'datetime'}],
    'pandas_version': '0.20.0',
    'primaryKey': ['idx']}

    Notes
    -----
    See `_as_json_table_type` for conversion types.
    Timedeltas as converted to ISO8601 duration format with
    9 decimal places after the secnods field for nanosecond precision.

    Categoricals are converted to the `any` dtype, and use the `enum` field
    constraint to list the allowed values. The `ordered` attribute is included
    in an `ordered` field.
    """
    if index is True:
        data = set_default_names(data)

    schema = {}
    fields = []

    if index:
        if data.index.nlevels > 1:
            for level in data.index.levels:
                fields.append(convert_pandas_type_to_json_field(level))
        else:
            fields.append(convert_pandas_type_to_json_field(data.index))

    if data.ndim > 1:
        for column, s in data.iteritems():
            fields.append(convert_pandas_type_to_json_field(s))
    else:
        fields.append(convert_pandas_type_to_json_field(data))

    schema['fields'] = fields
    if index and data.index.is_unique and primary_key is None:
        if data.index.nlevels == 1:
            schema['primaryKey'] = [data.index.name]
        else:
            schema['primaryKey'] = data.index.names
    elif primary_key is not None:
        schema['primaryKey'] = primary_key

    if version:
        schema['pandas_version'] = '0.20.0'
    return schema


def parse_table_schema(json, precise_float):
    """
    Builds a DataFrame from a given schema

    Parameters
    ----------
    json :
        A JSON table schema
    precise_float : boolean
        Flag controlling precision when decoding string to double values, as
        dictated by ``read_json``

    Returns
    -------
    df : DataFrame

    Raises
    ------
    NotImplementedError
        If the JSON table schema contains either timezone or timedelta data

    Notes
    -----
        Because :func:`DataFrame.to_json` uses the string 'index' to denote a
        name-less :class:`Index`, this function sets the name of the returned
        :class:`DataFrame` to ``None`` when said string is encountered with a
        normal :class:`Index`. For a :class:`MultiIndex`, the same limitation
        applies to any strings beginning with 'level_'. Therefore, an
        :class:`Index` name of 'index'  and :class:`MultiIndex` names starting
        with 'level_' are not supported.

    See also
    --------
    build_table_schema : inverse function
    pandas.read_json
    """
    table = loads(json, precise_float=precise_float)
    col_order = [field['name'] for field in table['schema']['fields']]
    df = DataFrame(table['data'])[col_order]

    dtypes = {field['name']: convert_json_field_to_pandas_type(field)
              for field in table['schema']['fields']}

    # Cannot directly use as_type with timezone data on object; raise for now
    if any(str(x).startswith('datetime64[ns, ') for x in dtypes.values()):
        raise NotImplementedError('table="orient" can not yet read timezone '
                                  'data')

    # No ISO constructor for Timedelta as of yet, so need to raise
    if 'timedelta64' in dtypes.values():
        raise NotImplementedError('table="orient" can not yet read '
                                  'ISO-formatted Timedelta data')

    df = df.astype(dtypes)

    df = df.set_index(table['schema']['primaryKey'])
    if len(df.index.names) == 1:
        if df.index.name == 'index':
            df.index.name = None
    else:
        df.index.names = [None if x.startswith('level_') else x for x in
                          df.index.names]

    return df
