# -*- coding: utf-8 -*-
# Authors: Adrin Jalali <adrin.jalali@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin
from .validation import check_array, column_or_1d


class NamedArray(NDArrayOperatorsMixin):
    _feature_names = None
    _data = None

    def __init__(self, data, feature_names=None):
        self.data = data
        self.feature_names = feature_names

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        value = check_array(value)
        self._data = value

        if self.feature_names is None:
            return

        if len(self.feature_names) != self._col_count(value):
            self._feature_names = None

    @property
    def feature_names(self):
        return self._feature_names

    @feature_names.setter
    def feature_names(self, value):
        if value is None:
            self._feature_names = None
            return

        value = column_or_1d(value)
        col_count = self._col_count(self.data)
        if len(value) != col_count:
            raise ValueError("{} column names provided, but data has {} "
                             "columns".format(len(value), col_count))

        self._feature_names = value

    def _col_count(self, value):
        if value.ndim == 1:
            return 1
        else:
            return value.shape[1]

    def __getattr__(self, name):
        return getattr(self._data, name)

    def __getitem__(self, slice):
        return self.data[slice]

    def __array__(self, *args, **kwargs):
        return self.data.__array__(*args, **kwargs)

    def __repr__(self):
        prefix = self.__class__.__name__ + '('
        base_repr = np.array2string(self.data,
                                    prefix=prefix)
        return (prefix + base_repr
                + f',\n           feature_names={self.feature_names!r})')
