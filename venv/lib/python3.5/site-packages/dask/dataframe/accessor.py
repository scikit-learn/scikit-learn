from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd
from toolz import partial

from ..utils import derived_from
from .utils import is_categorical_dtype, PANDAS_VERSION


def maybe_wrap_pandas(obj, x):
    if isinstance(x, np.ndarray):
        if isinstance(obj, pd.Series):
            return pd.Series(x, index=obj.index, dtype=x.dtype)
        return pd.Index(x)
    return x


class Accessor(object):
    """
    Base class for pandas Accessor objects cat, dt, and str.

    Notes
    -----
    Subclasses should define the following attributes:

    * _accessor
    * _accessor_name
    """
    _not_implemented = set()

    def __init__(self, series):
        from .core import Series
        if not isinstance(series, Series):
            raise ValueError('Accessor cannot be initialized')
        self._validate(series)
        self._series = series

    def _validate(self, series):
        pass

    @staticmethod
    def _delegate_property(obj, accessor, attr):
        out = getattr(getattr(obj, accessor, obj), attr)
        return maybe_wrap_pandas(obj, out)

    @staticmethod
    def _delegate_method(obj, accessor, attr, args, kwargs):
        out = getattr(getattr(obj, accessor, obj), attr)(*args, **kwargs)
        return maybe_wrap_pandas(obj, out)

    def _property_map(self, attr):
        meta = self._delegate_property(self._series._meta,
                                       self._accessor_name, attr)
        token = '%s-%s' % (self._accessor_name, attr)
        return self._series.map_partitions(self._delegate_property,
                                           self._accessor_name, attr,
                                           token=token, meta=meta)

    def _function_map(self, attr, *args, **kwargs):
        meta = self._delegate_method(self._series._meta_nonempty,
                                     self._accessor_name, attr, args, kwargs)
        token = '%s-%s' % (self._accessor_name, attr)
        return self._series.map_partitions(self._delegate_method,
                                           self._accessor_name, attr, args,
                                           kwargs, meta=meta, token=token)

    @property
    def _delegates(self):
        return set(dir(self._accessor)).difference(self._not_implemented)

    def __dir__(self):
        o = self._delegates
        o.update(self.__dict__)
        o.update(dir(type(self)))
        return list(o)

    def __getattr__(self, key):
        if key in self._delegates:
            if isinstance(getattr(self._accessor, key), property):
                return self._property_map(key)
            else:
                return partial(self._function_map, key)
        else:
            raise AttributeError(key)


class DatetimeAccessor(Accessor):
    """ Accessor object for datetimelike properties of the Series values.

    Examples
    --------

    >>> s.dt.microsecond  # doctest: +SKIP
    """
    _accessor = pd.Series.dt
    _accessor_name = 'dt'


class StringAccessor(Accessor):
    """ Accessor object for string properties of the Series values.

    Examples
    --------

    >>> s.str.lower()  # doctest: +SKIP
    """
    _accessor = pd.Series.str
    _accessor_name = 'str'
    _not_implemented = {'get_dummies'}

    def _validate(self, series):
        if not (series.dtype == 'object' or (
                is_categorical_dtype(series) and
                series.cat.categories.dtype == 'object')):
            raise AttributeError("Can only use .str accessor with object dtype")

    @derived_from(pd.core.strings.StringMethods)
    def split(self, pat=None, n=-1):
        return self._function_map('split', pat=pat, n=n)

    @derived_from(pd.core.strings.StringMethods)
    def cat(self, others=None, sep=None, na_rep=None):
        from .core import Series, Index
        if others is None:
            raise NotImplementedError("x.str.cat() with `others == None`")

        valid_types = (Series, Index, pd.Series, pd.Index)
        if isinstance(others, valid_types):
            others = [others]
        elif not all(isinstance(a, valid_types) for a in others):
            raise TypeError("others must be Series/Index")

        return self._series.map_partitions(str_cat, *others, sep=sep,
                                           na_rep=na_rep, meta=self._series._meta)

    @derived_from(pd.core.strings.StringMethods)
    def extractall(self, pat, flags=0):
        # TODO: metadata inference here won't be necessary for pandas >= 0.23.0
        meta = self._series._meta.str.extractall(pat, flags=flags)
        if PANDAS_VERSION < '0.23.0':
            index_name = self._series.index.name
            meta.index = pd.MultiIndex(levels=[[], []],
                                       labels=[[], []],
                                       names=[index_name, 'match'])
        return self._series.map_partitions(str_extractall, pat, flags,
                                           meta=meta, token='str-extractall')

    def __getitem__(self, index):
        return self._series.map_partitions(str_get, index,
                                           meta=self._series._meta)


def str_extractall(series, pat, flags):
    return series.str.extractall(pat, flags=flags)


def str_get(series, index):
    """ Implements series.str[index] """
    return series.str[index]


def str_cat(self, *others, **kwargs):
    return self.str.cat(others=others, **kwargs)
