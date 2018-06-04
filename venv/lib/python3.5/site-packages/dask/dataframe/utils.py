from __future__ import absolute_import, division, print_function

import re
import textwrap
from distutils.version import LooseVersion

from collections import Iterator
import sys
import traceback
from contextlib import contextmanager

import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas.api.types import is_categorical_dtype, is_scalar
try:
    from pandas.api.types import is_datetime64tz_dtype
except ImportError:
    # pandas < 0.19.2
    from pandas.core.common import is_datetime64tz_dtype

from ..core import get_deps
from ..local import get_sync
from ..utils import asciitable, is_arraylike


PANDAS_VERSION = LooseVersion(pd.__version__)


def shard_df_on_index(df, divisions):
    """ Shard a DataFrame by ranges on its index

    Examples
    --------

    >>> df = pd.DataFrame({'a': [0, 10, 20, 30, 40], 'b': [5, 4 ,3, 2, 1]})
    >>> df
        a  b
    0   0  5
    1  10  4
    2  20  3
    3  30  2
    4  40  1

    >>> shards = list(shard_df_on_index(df, [2, 4]))
    >>> shards[0]
        a  b
    0   0  5
    1  10  4

    >>> shards[1]
        a  b
    2  20  3
    3  30  2

    >>> shards[2]
        a  b
    4  40  1

    >>> list(shard_df_on_index(df, []))[0]  # empty case
        a  b
    0   0  5
    1  10  4
    2  20  3
    3  30  2
    4  40  1
    """

    if isinstance(divisions, Iterator):
        divisions = list(divisions)
    if not len(divisions):
        yield df
    else:
        divisions = np.array(divisions)
        df = df.sort_index()
        index = df.index
        if is_categorical_dtype(index):
            index = index.as_ordered()
        indices = index.searchsorted(divisions)
        yield df.iloc[:indices[0]]
        for i in range(len(indices) - 1):
            yield df.iloc[indices[i]: indices[i + 1]]
        yield df.iloc[indices[-1]:]


_META_TYPES = "meta : pd.DataFrame, pd.Series, dict, iterable, tuple, optional"
_META_DESCRIPTION = """\
An empty ``pd.DataFrame`` or ``pd.Series`` that matches the dtypes and
column names of the output. This metadata is necessary for many algorithms
in dask dataframe to work.  For ease of use, some alternative inputs are
also available. Instead of a ``DataFrame``, a ``dict`` of ``{name: dtype}``
or iterable of ``(name, dtype)`` can be provided. Instead of a series, a
tuple of ``(name, dtype)`` can be used. If not provided, dask will try to
infer the metadata. This may lead to unexpected results, so providing
``meta`` is recommended. For more information, see
``dask.dataframe.utils.make_meta``.
"""


def insert_meta_param_description(*args, **kwargs):
    """Replace `$META` in docstring with param description.

    If pad keyword is provided, will pad description by that number of
    spaces (default is 8)."""
    if not args:
        return lambda f: insert_meta_param_description(f, **kwargs)
    f = args[0]
    indent = " " * kwargs.get('pad', 8)
    body = textwrap.wrap(_META_DESCRIPTION, initial_indent=indent,
                         subsequent_indent=indent, width=78)
    descr = '{0}\n{1}'.format(_META_TYPES, '\n'.join(body))
    if f.__doc__:
        if '$META' in f.__doc__:
            f.__doc__ = f.__doc__.replace('$META', descr)
        else:
            # Put it at the end of the parameters section
            parameter_header = 'Parameters\n%s----------' % indent[4:]
            first, last = re.split('Parameters\\n[ ]*----------', f.__doc__)
            parameters, rest = last.split('\n\n', 1)
            f.__doc__ = '{0}{1}{2}\n{3}{4}\n\n{5}'.format(first, parameter_header,
                                                          parameters, indent[4:],
                                                          descr, rest)
    return f


@contextmanager
def raise_on_meta_error(funcname=None, udf=False):
    """Reraise errors in this block to show metadata inference failure.

    Parameters
    ----------
    funcname : str, optional
        If provided, will be added to the error message to indicate the
        name of the method that failed.
    """
    try:
        yield
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = ''.join(traceback.format_tb(exc_traceback))
        msg = "Metadata inference failed{0}.\n\n"
        if udf:
            msg += ("You have supplied a custom function and Dask is unable to \n"
                    "determine the type of output that that function returns. \n\n"
                    "To resolve this please provide a meta= keyword.\n"
                    "The docstring of the Dask function you ran should have more information.\n\n")
        msg += ("Original error is below:\n"
                "------------------------\n"
                "{1}\n\n"
                "Traceback:\n"
                "---------\n"
                "{2}")
        msg = msg.format(" in `{0}`".format(funcname) if funcname else "", repr(e), tb)
        raise ValueError(msg)


UNKNOWN_CATEGORIES = '__UNKNOWN_CATEGORIES__'


def has_known_categories(x):
    """Returns whether the categories in `x` are known.

    Parameters
    ----------
    x : Series or CategoricalIndex
    """
    x = getattr(x, '_meta', x)
    if isinstance(x, pd.Series):
        return UNKNOWN_CATEGORIES not in x.cat.categories
    elif isinstance(x, pd.CategoricalIndex):
        return UNKNOWN_CATEGORIES not in x.categories
    raise TypeError("Expected Series or CategoricalIndex")


def strip_unknown_categories(x):
    """Replace any unknown categoricals with empty categoricals.

    Useful for preventing ``UNKNOWN_CATEGORIES`` from leaking into results.
    """
    if isinstance(x, (pd.Series, pd.DataFrame)):
        x = x.copy()
        if isinstance(x, pd.DataFrame):
            cat_mask = x.dtypes == 'category'
            if cat_mask.any():
                cats = cat_mask[cat_mask].index
                for c in cats:
                    if not has_known_categories(x[c]):
                        x[c].cat.set_categories([], inplace=True)
        elif isinstance(x, pd.Series):
            if is_categorical_dtype(x.dtype) and not has_known_categories(x):
                x.cat.set_categories([], inplace=True)
        if (isinstance(x.index, pd.CategoricalIndex) and not
                has_known_categories(x.index)):
            x.index = x.index.set_categories([])
    elif isinstance(x, pd.CategoricalIndex) and not has_known_categories(x):
        x = x.set_categories([])
    return x


def clear_known_categories(x, cols=None, index=True):
    """Set categories to be unknown.

    Parameters
    ----------
    x : DataFrame, Series, Index
    cols : iterable, optional
        If x is a DataFrame, set only categoricals in these columns to unknown.
        By default, all categorical columns are set to unknown categoricals
    index : bool, optional
        If True and x is a Series or DataFrame, set the clear known categories
        in the index as well.
    """
    if isinstance(x, (pd.Series, pd.DataFrame)):
        x = x.copy()
        if isinstance(x, pd.DataFrame):
            mask = x.dtypes == 'category'
            if cols is None:
                cols = mask[mask].index
            elif not mask.loc[cols].all():
                raise ValueError("Not all columns are categoricals")
            for c in cols:
                x[c].cat.set_categories([UNKNOWN_CATEGORIES], inplace=True)
        elif isinstance(x, pd.Series):
            if is_categorical_dtype(x.dtype):
                x.cat.set_categories([UNKNOWN_CATEGORIES], inplace=True)
        if index and isinstance(x.index, pd.CategoricalIndex):
            x.index = x.index.set_categories([UNKNOWN_CATEGORIES])
    elif isinstance(x, pd.CategoricalIndex):
        x = x.set_categories([UNKNOWN_CATEGORIES])
    return x


def _empty_series(name, dtype, index=None):
    if isinstance(dtype, str) and dtype == 'category':
        return pd.Series(pd.Categorical([UNKNOWN_CATEGORIES]),
                         name=name, index=index).iloc[:0]
    return pd.Series([], dtype=dtype, name=name, index=index)


def make_meta(x, index=None):
    """Create an empty pandas object containing the desired metadata.

    Parameters
    ----------
    x : dict, tuple, list, pd.Series, pd.DataFrame, pd.Index, dtype, scalar
        To create a DataFrame, provide a `dict` mapping of `{name: dtype}`, or
        an iterable of `(name, dtype)` tuples. To create a `Series`, provide a
        tuple of `(name, dtype)`. If a pandas object, names, dtypes, and index
        should match the desired output. If a dtype or scalar, a scalar of the
        same dtype is returned.
    index :  pd.Index, optional
        Any pandas index to use in the metadata. If none provided, a
        `RangeIndex` will be used.

    Examples
    --------
    >>> make_meta([('a', 'i8'), ('b', 'O')])
    Empty DataFrame
    Columns: [a, b]
    Index: []
    >>> make_meta(('a', 'f8'))
    Series([], Name: a, dtype: float64)
    >>> make_meta('i8')
    1
    """
    if hasattr(x, '_meta'):
        return x._meta
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.iloc[0:0]
    elif isinstance(x, pd.Index):
        return x[0:0]
    elif is_arraylike(x):
        return x[:0]
    index = index if index is None else index[0:0]

    if isinstance(x, dict):
        return pd.DataFrame({c: _empty_series(c, d, index=index)
                             for (c, d) in x.items()}, index=index)
    if isinstance(x, tuple) and len(x) == 2:
        return _empty_series(x[0], x[1], index=index)
    elif isinstance(x, (list, tuple)):
        if not all(isinstance(i, tuple) and len(i) == 2 for i in x):
            raise ValueError("Expected iterable of tuples of (name, dtype), "
                             "got {0}".format(x))
        return pd.DataFrame({c: _empty_series(c, d, index=index) for (c, d) in x},
                            columns=[c for c, d in x], index=index)
    elif not hasattr(x, 'dtype') and x is not None:
        # could be a string, a dtype object, or a python type. Skip `None`,
        # because it is implictly converted to `dtype('f8')`, which we don't
        # want here.
        try:
            dtype = np.dtype(x)
            return _scalar_from_dtype(dtype)
        except Exception:
            # Continue on to next check
            pass

    if is_scalar(x):
        return _nonempty_scalar(x)

    raise TypeError("Don't know how to create metadata from {0}".format(x))


if PANDAS_VERSION >= "0.20.0":
    _numeric_index_types = (pd.Int64Index, pd.Float64Index, pd.UInt64Index)
else:
    _numeric_index_types = (pd.Int64Index, pd.Float64Index)


def _nonempty_index(idx):
    typ = type(idx)
    if typ is pd.RangeIndex:
        return pd.RangeIndex(2, name=idx.name)
    elif typ in _numeric_index_types:
        return typ([1, 2], name=idx.name)
    elif typ is pd.Index:
        return pd.Index(['a', 'b'], name=idx.name)
    elif typ is pd.DatetimeIndex:
        start = '1970-01-01'
        # Need a non-monotonic decreasing index to avoid issues with
        # partial string indexing see https://github.com/dask/dask/issues/2389
        # and https://github.com/pandas-dev/pandas/issues/16515
        # This doesn't mean `_meta_nonempty` should ever rely on
        # `self.monotonic_increasing` or `self.monotonic_decreasing`
        data = [start, '1970-01-02'] if idx.freq is None else None
        return pd.DatetimeIndex(data, start=start, periods=2, freq=idx.freq,
                                tz=idx.tz, name=idx.name)
    elif typ is pd.PeriodIndex:
        return pd.PeriodIndex(start='1970-01-01', periods=2, freq=idx.freq,
                              name=idx.name)
    elif typ is pd.TimedeltaIndex:
        start = np.timedelta64(1, 'D')
        data = [start, start + 1] if idx.freq is None else None
        return pd.TimedeltaIndex(data, start=start, periods=2, freq=idx.freq,
                                 name=idx.name)
    elif typ is pd.CategoricalIndex:
        if len(idx.categories) == 0:
            data = _nonempty_index(idx.categories)
            cats = None
        else:
            data = _nonempty_index(_nonempty_index(idx.categories))
            cats = idx.categories
        return pd.CategoricalIndex(data, categories=cats,
                                   ordered=idx.ordered, name=idx.name)
    elif typ is pd.MultiIndex:
        levels = [_nonempty_index(l) for l in idx.levels]
        labels = [[0, 0] for i in idx.levels]
        return pd.MultiIndex(levels=levels, labels=labels, names=idx.names)
    raise TypeError("Don't know how to handle index of "
                    "type {0}".format(type(idx).__name__))


_simple_fake_mapping = {
    'b': np.bool_(True),
    'V': np.void(b' '),
    'M': np.datetime64('1970-01-01'),
    'm': np.timedelta64(1),
    'S': np.str_('foo'),
    'a': np.str_('foo'),
    'U': np.unicode_('foo'),
    'O': 'foo'
}


def _scalar_from_dtype(dtype):
    if dtype.kind in ('i', 'f', 'u'):
        return dtype.type(1)
    elif dtype.kind == 'c':
        return dtype.type(complex(1, 0))
    elif dtype.kind in _simple_fake_mapping:
        o = _simple_fake_mapping[dtype.kind]
        return o.astype(dtype) if dtype.kind in ('m', 'M') else o
    else:
        raise TypeError("Can't handle dtype: {0}".format(dtype))


def _nonempty_scalar(x):
    if isinstance(x, (pd.Timestamp, pd.Timedelta, pd.Period)):
        return x
    elif np.isscalar(x):
        dtype = x.dtype if hasattr(x, 'dtype') else np.dtype(type(x))
        return _scalar_from_dtype(dtype)
    else:
        raise TypeError("Can't handle meta of type "
                        "'{0}'".format(type(x).__name__))


def _nonempty_series(s, idx):

    dtype = s.dtype
    if is_datetime64tz_dtype(dtype):
        entry = pd.Timestamp('1970-01-01', tz=dtype.tz)
        data = [entry, entry]
    elif is_categorical_dtype(dtype):
        if len(s.cat.categories):
            data = [s.cat.categories[0]] * 2
            cats = s.cat.categories
        else:
            data = _nonempty_index(s.cat.categories)
            cats = None
        data = pd.Categorical(data, categories=cats,
                              ordered=s.cat.ordered)
    else:
        entry = _scalar_from_dtype(dtype)
        data = np.array([entry, entry], dtype=dtype)

    return pd.Series(data, name=s.name, index=idx)


def meta_nonempty(x):
    """Create a nonempty pandas object from the given metadata.

    Returns a pandas DataFrame, Series, or Index that contains two rows
    of fake data.
    """
    if isinstance(x, pd.Index):
        return _nonempty_index(x)
    elif isinstance(x, pd.Series):
        idx = _nonempty_index(x.index)
        return _nonempty_series(x, idx)
    elif isinstance(x, pd.DataFrame):
        idx = _nonempty_index(x.index)
        data = {i: _nonempty_series(x.iloc[:, i], idx)
                for i, c in enumerate(x.columns)}
        res = pd.DataFrame(data, index=idx,
                           columns=np.arange(len(x.columns)))
        res.columns = x.columns
        return res
    elif is_scalar(x):
        return _nonempty_scalar(x)
    else:
        raise TypeError("Expected Index, Series, DataFrame, or scalar, "
                        "got {0}".format(type(x).__name__))


def check_meta(x, meta, funcname=None, numeric_equal=True):
    """Check that the dask metadata matches the result.

    If metadata matches, ``x`` is passed through unchanged. A nice error is
    raised if metadata doesn't match.

    Parameters
    ----------
    x : DataFrame, Series, or Index
    meta : DataFrame, Series, or Index
        The expected metadata that ``x`` should match
    funcname : str, optional
        The name of the function in which the metadata was specified. If
        provided, the function name will be included in the error message to be
        more helpful to users.
    numeric_equal : bool, optionl
        If True, integer and floating dtypes compare equal. This is useful due
        to panda's implicit conversion of integer to floating upon encountering
        missingness, which is hard to infer statically.
    """
    eq_types = {'i', 'f'} if numeric_equal else {}

    def equal_dtypes(a, b):
        if is_categorical_dtype(a) != is_categorical_dtype(b):
            return False
        if (a is '-' or b is '-'):
            return False
        if is_categorical_dtype(a) and is_categorical_dtype(b):
            # Pandas 0.21 CategoricalDtype compat
            if (PANDAS_VERSION >= '0.21.0' and
                    (UNKNOWN_CATEGORIES in a.categories or
                     UNKNOWN_CATEGORIES in b.categories)):
                return True
            return a == b
        return (a.kind in eq_types and b.kind in eq_types) or (a == b)

    if not isinstance(meta, (pd.Series, pd.Index, pd.DataFrame)):
        raise TypeError("Expected partition to be DataFrame, Series, or "
                        "Index, got `%s`" % type(meta).__name__)

    if type(x) != type(meta):
        errmsg = ("Expected partition of type `%s` but got "
                  "`%s`" % (type(meta).__name__, type(x).__name__))
    elif isinstance(meta, pd.DataFrame):
        dtypes = pd.concat([x.dtypes, meta.dtypes], axis=1)
        bad = [(col, a, b) for col, a, b in dtypes.fillna('-').itertuples()
               if not equal_dtypes(a, b)]
        if not bad:
            return x
        errmsg = ("Partition type: `%s`\n%s" %
                  (type(meta).__name__,
                   asciitable(['Column', 'Found', 'Expected'], bad)))
    else:
        if equal_dtypes(x.dtype, meta.dtype):
            return x
        errmsg = ("Partition type: `%s`\n%s" %
                  (type(meta).__name__,
                   asciitable(['', 'dtype'], [('Found', x.dtype),
                                              ('Expected', meta.dtype)])))

    raise ValueError("Metadata mismatch found%s.\n\n"
                     "%s" % ((" in `%s`" % funcname if funcname else ""),
                             errmsg))


def index_summary(idx, name=None):
    """Summarized representation of an Index.
    """
    n = len(idx)
    if name is None:
        name = idx.__class__.__name__
    if n:
        head = idx[0]
        tail = idx[-1]
        summary = ', {} to {}'.format(head, tail)
    else:
        summary = ''

    return "{}: {} entries{}".format(name, n, summary)


###############################################################
# Testing
###############################################################


def _check_dask(dsk, check_names=True, check_dtypes=True, result=None):
    import dask.dataframe as dd
    if hasattr(dsk, 'dask'):
        if result is None:
            result = dsk.compute(get=get_sync)
        if isinstance(dsk, dd.Index):
            assert isinstance(result, pd.Index), type(result)
            assert isinstance(dsk._meta, pd.Index), type(dsk._meta)
            if check_names:
                assert dsk.name == result.name
                assert dsk._meta.name == result.name
                if isinstance(result, pd.MultiIndex):
                    assert result.names == dsk._meta.names
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
        elif isinstance(dsk, dd.Series):
            assert isinstance(result, pd.Series), type(result)
            assert isinstance(dsk._meta, pd.Series), type(dsk._meta)
            if check_names:
                assert dsk.name == result.name, (dsk.name, result.name)
                assert dsk._meta.name == result.name
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
            _check_dask(dsk.index, check_names=check_names,
                        check_dtypes=check_dtypes, result=result.index)
        elif isinstance(dsk, dd.DataFrame):
            assert isinstance(result, pd.DataFrame), type(result)
            assert isinstance(dsk.columns, pd.Index), type(dsk.columns)
            assert isinstance(dsk._meta, pd.DataFrame), type(dsk._meta)
            if check_names:
                tm.assert_index_equal(dsk.columns, result.columns)
                tm.assert_index_equal(dsk._meta.columns, result.columns)
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
            _check_dask(dsk.index, check_names=check_names,
                        check_dtypes=check_dtypes, result=result.index)
        elif isinstance(dsk, dd.core.Scalar):
            assert (np.isscalar(result) or
                    isinstance(result, (pd.Timestamp, pd.Timedelta)))
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
        else:
            msg = 'Unsupported dask instance {0} found'.format(type(dsk))
            raise AssertionError(msg)
        return result
    return dsk


def _maybe_sort(a):
    # sort by value, then index
    try:
        if isinstance(a, pd.DataFrame):
            a = a.sort_values(by=a.columns.tolist())
        else:
            a = a.sort_values()
    except (TypeError, IndexError, ValueError):
        pass
    return a.sort_index()


def assert_eq(a, b, check_names=True, check_dtypes=True,
              check_divisions=True, check_index=True, **kwargs):
    if check_divisions:
        assert_divisions(a)
        assert_divisions(b)
        if hasattr(a, 'divisions') and hasattr(b, 'divisions'):
            at = type(np.asarray(a.divisions).tolist()[0])  # numpy to python
            bt = type(np.asarray(b.divisions).tolist()[0])  # scalar conversion
            assert at == bt, (at, bt)
    assert_sane_keynames(a)
    assert_sane_keynames(b)
    a = _check_dask(a, check_names=check_names, check_dtypes=check_dtypes)
    b = _check_dask(b, check_names=check_names, check_dtypes=check_dtypes)
    if not check_index:
        a = a.reset_index(drop=True)
        b = b.reset_index(drop=True)
    if isinstance(a, pd.DataFrame):
        a = _maybe_sort(a)
        b = _maybe_sort(b)
        tm.assert_frame_equal(a, b, **kwargs)
    elif isinstance(a, pd.Series):
        a = _maybe_sort(a)
        b = _maybe_sort(b)
        tm.assert_series_equal(a, b, check_names=check_names, **kwargs)
    elif isinstance(a, pd.Index):
        tm.assert_index_equal(a, b, **kwargs)
    else:
        if a == b:
            return True
        else:
            if np.isnan(a):
                assert np.isnan(b)
            else:
                assert np.allclose(a, b)
    return True


def assert_dask_graph(dask, label):
    if hasattr(dask, 'dask'):
        dask = dask.dask
    assert isinstance(dask, dict)
    for k in dask:
        if isinstance(k, tuple):
            k = k[0]
        if k.startswith(label):
            return True
    raise AssertionError("given dask graph doesn't contain label: {label}"
                         .format(label=label))


def assert_divisions(ddf):
    if not hasattr(ddf, 'divisions'):
        return
    if not hasattr(ddf, 'index'):
        return
    if not ddf.known_divisions:
        return

    def index(x):
        return (x if isinstance(x, pd.Index)
                else x.index.get_level_values(0))

    results = get_sync(ddf.dask, ddf.__dask_keys__())
    for i, df in enumerate(results[:-1]):
        if len(df):
            assert index(df).min() >= ddf.divisions[i]
            assert index(df).max() < ddf.divisions[i + 1]

    if len(results[-1]):
        assert index(results[-1]).min() >= ddf.divisions[-2]
        assert index(results[-1]).max() <= ddf.divisions[-1]


def assert_sane_keynames(ddf):
    if not hasattr(ddf, 'dask'):
        return
    for k in ddf.dask.keys():
        while isinstance(k, tuple):
            k = k[0]
        assert isinstance(k, (str, bytes))
        assert len(k) < 100
        assert ' ' not in k
        if sys.version_info[0] >= 3:
            assert k.split('-')[0].isidentifier()


def assert_dask_dtypes(ddf, res, numeric_equal=True):
    """Check that the dask metadata matches the result.

    If `numeric_equal`, integer and floating dtypes compare equal. This is
    useful due to the implicit conversion of integer to floating upon
    encountering missingness, which is hard to infer statically."""

    eq_types = {'O', 'S', 'U', 'a'}     # treat object and strings alike
    if numeric_equal:
        eq_types.update(('i', 'f'))

    if isinstance(res, pd.DataFrame):
        for col, a, b in pd.concat([ddf._meta.dtypes, res.dtypes],
                                   axis=1).itertuples():
            assert (a.kind in eq_types and b.kind in eq_types) or (a == b)
    elif isinstance(res, (pd.Series, pd.Index)):
        a = ddf._meta.dtype
        b = res.dtype
        assert (a.kind in eq_types and b.kind in eq_types) or (a == b)
    else:
        if hasattr(ddf._meta, 'dtype'):
            a = ddf._meta.dtype
            if not hasattr(res, 'dtype'):
                assert np.isscalar(res)
                b = np.dtype(type(res))
            else:
                b = res.dtype
            assert (a.kind in eq_types and b.kind in eq_types) or (a == b)
        else:
            assert type(ddf._meta) == type(res)


def assert_max_deps(x, n, eq=True):
    dependencies, dependents = get_deps(x.dask)
    if eq:
        assert max(map(len, dependencies.values())) == n
    else:
        assert max(map(len, dependencies.values())) <= n
