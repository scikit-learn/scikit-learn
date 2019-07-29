# -*- coding: utf-8 -*-
# Authors: Adrin Jalali <adrin.jalali@gmail.com>
#
# License: BSD 3 clause

import numpy as np

# NDArrayOperatorsMixin was added in numpy 1.13
# TODO: cleanup once we support numpy 1.13+
try:
    from numpy.lib.mixins import NDArrayOperatorsMixin
except ImportError:
    raise NotImplementedError("In order to use NamedAraay, please upgrade your "
                              "numpy to 1.13+!")

from .validation import check_array, column_or_1d


class NamedArray(NDArrayOperatorsMixin):
    """A wrapper to a numpy ndarray holding some metadata about the data.

    Instances of this object behave like a numpy array, and loose all metadata
    information in numerical operations.

    Parameters
    ----------
    data: array-like
        A one or two dimensional array like data.

    feature_names: list or array of strings, or None, default=None
        Feature names associated with the columns of the data. The number of
        columns should always be the same as the number of feature names.
        Setting the `data` of an instance, would result in `feature_names` to
        be `None` if the number of columns do not match the number of stored
        feature names.
    """

    def __init__(self, data, feature_names=None):
        data = check_array(data, ensure_2d=False)
        self._data = data
        self.feature_names = feature_names

    @property
    def feature_names(self):
        return self._feature_names

    @feature_names.setter
    def feature_names(self, value):
        if value is None:
            self._feature_names = None
            return

        if np.isscalar(value):
            value = [value]
        value = column_or_1d(value)
        col_count = self._col_count(self._data)
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

    def __dir__(self):
        return list(set(dir(NamedArray)).union(set(dir(self._data))))

    def __getitem__(self, slice):
        return self._data[slice]

    def __repr__(self):
        prefix = self.__class__.__name__ + '('
        base_repr = np.array2string(self._data,
                                    prefix=prefix)
        return (prefix + base_repr
                + ',\n           feature_names={})'.format(
                    str(self.feature_names)))

    def todataframe(self):
        """Returns a `pandas.DataFrame` with set column names."""
        import pandas as pd
        return pd.DataFrame(self._data, columns=self.feature_names)
