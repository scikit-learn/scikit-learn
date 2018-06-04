import warnings
import numpy as np
from pandas import compat
from pandas._libs import reduction
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.common import (
    is_extension_type,
    is_sequence)
from pandas.util._decorators import cache_readonly

from pandas.io.formats.printing import pprint_thing


def frame_apply(obj, func, axis=0, broadcast=None,
                raw=False, reduce=None, result_type=None,
                ignore_failures=False,
                args=None, kwds=None):
    """ construct and return a row or column based frame apply object """

    axis = obj._get_axis_number(axis)
    if axis == 0:
        klass = FrameRowApply
    elif axis == 1:
        klass = FrameColumnApply

    return klass(obj, func, broadcast=broadcast,
                 raw=raw, reduce=reduce, result_type=result_type,
                 ignore_failures=ignore_failures,
                 args=args, kwds=kwds)


class FrameApply(object):

    def __init__(self, obj, func, broadcast, raw, reduce, result_type,
                 ignore_failures, args, kwds):
        self.obj = obj
        self.raw = raw
        self.ignore_failures = ignore_failures
        self.args = args or ()
        self.kwds = kwds or {}

        if result_type not in [None, 'reduce', 'broadcast', 'expand']:
            raise ValueError("invalid value for result_type, must be one "
                             "of {None, 'reduce', 'broadcast', 'expand'}")

        if broadcast is not None:
            warnings.warn("The broadcast argument is deprecated and will "
                          "be removed in a future version. You can specify "
                          "result_type='broadcast' to broadcast the result "
                          "to the original dimensions",
                          FutureWarning, stacklevel=4)
            if broadcast:
                result_type = 'broadcast'

        if reduce is not None:
            warnings.warn("The reduce argument is deprecated and will "
                          "be removed in a future version. You can specify "
                          "result_type='reduce' to try to reduce the result "
                          "to the original dimensions",
                          FutureWarning, stacklevel=4)
            if reduce:

                if result_type is not None:
                    raise ValueError(
                        "cannot pass both reduce=True and result_type")

                result_type = 'reduce'

        self.result_type = result_type

        # curry if needed
        if kwds or args and not isinstance(func, np.ufunc):
            def f(x):
                return func(x, *args, **kwds)
        else:
            f = func

        self.f = f

        # results
        self.result = None
        self.res_index = None
        self.res_columns = None

    @property
    def columns(self):
        return self.obj.columns

    @property
    def index(self):
        return self.obj.index

    @cache_readonly
    def values(self):
        return self.obj.values

    @cache_readonly
    def dtypes(self):
        return self.obj.dtypes

    @property
    def agg_axis(self):
        return self.obj._get_agg_axis(self.axis)

    def get_result(self):
        """ compute the results """

        # all empty
        if len(self.columns) == 0 and len(self.index) == 0:
            return self.apply_empty_result()

        # string dispatch
        if isinstance(self.f, compat.string_types):
            # Support for `frame.transform('method')`
            # Some methods (shift, etc.) require the axis argument, others
            # don't, so inspect and insert if nescessary.
            func = getattr(self.obj, self.f)
            sig = compat.signature(func)
            if 'axis' in sig.args:
                self.kwds['axis'] = self.axis
            return func(*self.args, **self.kwds)

        # ufunc
        elif isinstance(self.f, np.ufunc):
            with np.errstate(all='ignore'):
                results = self.f(self.values)
            return self.obj._constructor(data=results, index=self.index,
                                         columns=self.columns, copy=False)

        # broadcasting
        if self.result_type == 'broadcast':
            return self.apply_broadcast()

        # one axis empty
        elif not all(self.obj.shape):
            return self.apply_empty_result()

        # raw
        elif self.raw and not self.obj._is_mixed_type:
            return self.apply_raw()

        return self.apply_standard()

    def apply_empty_result(self):
        """
        we have an empty result; at least 1 axis is 0

        we will try to apply the function to an empty
        series in order to see if this is a reduction function
        """

        # we are not asked to reduce or infer reduction
        # so just return a copy of the existing object
        if self.result_type not in ['reduce', None]:
            return self.obj.copy()

        # we may need to infer
        reduce = self.result_type == 'reduce'

        from pandas import Series
        if not reduce:

            EMPTY_SERIES = Series([])
            try:
                r = self.f(EMPTY_SERIES, *self.args, **self.kwds)
                reduce = not isinstance(r, Series)
            except Exception:
                pass

        if reduce:
            return self.obj._constructor_sliced(np.nan, index=self.agg_axis)
        else:
            return self.obj.copy()

    def apply_raw(self):
        """ apply to the values as a numpy array """

        try:
            result = reduction.reduce(self.values, self.f, axis=self.axis)
        except Exception:
            result = np.apply_along_axis(self.f, self.axis, self.values)

        # TODO: mixed type case
        if result.ndim == 2:
            return self.obj._constructor(result,
                                         index=self.index,
                                         columns=self.columns)
        else:
            return self.obj._constructor_sliced(result,
                                                index=self.agg_axis)

    def apply_broadcast(self, target):
        result_values = np.empty_like(target.values)

        # axis which we want to compare compliance
        result_compare = target.shape[0]

        for i, col in enumerate(target.columns):
            res = self.f(target[col])
            ares = np.asarray(res).ndim

            # must be a scalar or 1d
            if ares > 1:
                raise ValueError("too many dims to broadcast")
            elif ares == 1:

                # must match return dim
                if result_compare != len(res):
                    raise ValueError("cannot broadcast result")

            result_values[:, i] = res

        # we *always* preserve the original index / columns
        result = self.obj._constructor(result_values,
                                       index=target.index,
                                       columns=target.columns)
        return result

    def apply_standard(self):

        # try to reduce first (by default)
        # this only matters if the reduction in values is of different dtype
        # e.g. if we want to apply to a SparseFrame, then can't directly reduce

        # we cannot reduce using non-numpy dtypes,
        # as demonstrated in gh-12244
        if (self.result_type in ['reduce', None] and
                not self.dtypes.apply(is_extension_type).any()):

            # Create a dummy Series from an empty array
            from pandas import Series
            values = self.values
            index = self.obj._get_axis(self.axis)
            labels = self.agg_axis
            empty_arr = np.empty(len(index), dtype=values.dtype)
            dummy = Series(empty_arr, index=index, dtype=values.dtype)

            try:
                result = reduction.reduce(values, self.f,
                                          axis=self.axis,
                                          dummy=dummy,
                                          labels=labels)
                return self.obj._constructor_sliced(result, index=labels)
            except Exception:
                pass

        # compute the result using the series generator
        self.apply_series_generator()

        # wrap results
        return self.wrap_results()

    def apply_series_generator(self):
        series_gen = self.series_generator
        res_index = self.result_index

        i = None
        keys = []
        results = {}
        if self.ignore_failures:
            successes = []
            for i, v in enumerate(series_gen):
                try:
                    results[i] = self.f(v)
                    keys.append(v.name)
                    successes.append(i)
                except Exception:
                    pass

            # so will work with MultiIndex
            if len(successes) < len(res_index):
                res_index = res_index.take(successes)

        else:
            try:
                for i, v in enumerate(series_gen):
                    results[i] = self.f(v)
                    keys.append(v.name)
            except Exception as e:
                if hasattr(e, 'args'):

                    # make sure i is defined
                    if i is not None:
                        k = res_index[i]
                        e.args = e.args + ('occurred at index %s' %
                                           pprint_thing(k), )
                raise

        self.results = results
        self.res_index = res_index
        self.res_columns = self.result_columns

    def wrap_results(self):
        results = self.results

        # see if we can infer the results
        if len(results) > 0 and is_sequence(results[0]):

            return self.wrap_results_for_axis()

        # dict of scalars
        result = self.obj._constructor_sliced(results)
        result.index = self.res_index

        return result


class FrameRowApply(FrameApply):
    axis = 0

    def get_result(self):

        # dispatch to agg
        if isinstance(self.f, (list, dict)):
            return self.obj.aggregate(self.f, axis=self.axis,
                                      *self.args, **self.kwds)

        return super(FrameRowApply, self).get_result()

    def apply_broadcast(self):
        return super(FrameRowApply, self).apply_broadcast(self.obj)

    @property
    def series_generator(self):
        return (self.obj._ixs(i, axis=1)
                for i in range(len(self.columns)))

    @property
    def result_index(self):
        return self.columns

    @property
    def result_columns(self):
        return self.index

    def wrap_results_for_axis(self):
        """ return the results for the rows """

        results = self.results
        result = self.obj._constructor(data=results)

        if not isinstance(results[0], ABCSeries):
            try:
                result.index = self.res_columns
            except ValueError:
                pass

        try:
            result.columns = self.res_index
        except ValueError:
            pass

        return result


class FrameColumnApply(FrameApply):
    axis = 1

    def apply_broadcast(self):
        result = super(FrameColumnApply, self).apply_broadcast(self.obj.T)
        return result.T

    @property
    def series_generator(self):
        constructor = self.obj._constructor_sliced
        return (constructor(arr, index=self.columns, name=name)
                for i, (arr, name) in enumerate(zip(self.values,
                                                    self.index)))

    @property
    def result_index(self):
        return self.index

    @property
    def result_columns(self):
        return self.columns

    def wrap_results_for_axis(self):
        """ return the results for the columns """
        results = self.results

        # we have requested to expand
        if self.result_type == 'expand':
            result = self.infer_to_same_shape()

        # we have a non-series and don't want inference
        elif not isinstance(results[0], ABCSeries):
            from pandas import Series
            result = Series(results)
            result.index = self.res_index

        # we may want to infer results
        else:
            result = self.infer_to_same_shape()

        return result

    def infer_to_same_shape(self):
        """ infer the results to the same shape as the input object """
        results = self.results

        result = self.obj._constructor(data=results)
        result = result.T

        # set the index
        result.index = self.res_index

        # infer dtypes
        result = result.infer_objects()

        return result
