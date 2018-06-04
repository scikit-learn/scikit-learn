from sys import getsizeof
import operator
from datetime import timedelta

import numpy as np
from pandas._libs import index as libindex

from pandas.core.dtypes.common import (
    is_integer,
    is_scalar,
    is_int64_dtype)
from pandas.core.dtypes.generic import ABCSeries, ABCTimedeltaIndex

from pandas import compat
from pandas.compat import lrange, range, get_range_parameters
from pandas.compat.numpy import function as nv

import pandas.core.common as com
from pandas.core import ops
from pandas.core.indexes.base import Index, _index_shared_docs
from pandas.util._decorators import Appender, cache_readonly
import pandas.core.dtypes.concat as _concat
import pandas.core.indexes.base as ibase

from pandas.core.indexes.numeric import Int64Index


class RangeIndex(Int64Index):

    """
    Immutable Index implementing a monotonic integer range.

    RangeIndex is a memory-saving special case of Int64Index limited to
    representing monotonic ranges. Using RangeIndex may in some instances
    improve computing speed.

    This is the default index type used
    by DataFrame and Series when no explicit index is provided by the user.

    Parameters
    ----------
    start : int (default: 0), or other RangeIndex instance.
        If int and "stop" is not given, interpreted as "stop" instead.
    stop : int (default: 0)
    step : int (default: 1)
    name : object, optional
        Name to be stored in the index
    copy : bool, default False
        Unused, accepted for homogeneity with other index types.

    See Also
    --------
    Index : The base pandas Index type
    Int64Index : Index of int64 data

    Attributes
    ----------
    None

    Methods
    -------
    from_range
    """

    _typ = 'rangeindex'
    _engine_type = libindex.Int64Engine

    def __new__(cls, start=None, stop=None, step=None,
                dtype=None, copy=False, name=None, fastpath=False):

        if fastpath:
            return cls._simple_new(start, stop, step, name=name)

        cls._validate_dtype(dtype)

        # RangeIndex
        if isinstance(start, RangeIndex):
            if name is None:
                name = start.name
            return cls._simple_new(name=name,
                                   **dict(start._get_data_as_items()))

        # validate the arguments
        def _ensure_int(value, field):
            msg = ("RangeIndex(...) must be called with integers,"
                   " {value} was passed for {field}")
            if not is_scalar(value):
                raise TypeError(msg.format(value=type(value).__name__,
                                           field=field))
            try:
                new_value = int(value)
                assert(new_value == value)
            except (TypeError, ValueError, AssertionError):
                raise TypeError(msg.format(value=type(value).__name__,
                                           field=field))

            return new_value

        if com._all_none(start, stop, step):
            msg = "RangeIndex(...) must be called with integers"
            raise TypeError(msg)
        elif start is None:
            start = 0
        else:
            start = _ensure_int(start, 'start')
        if stop is None:
            stop = start
            start = 0
        else:
            stop = _ensure_int(stop, 'stop')
        if step is None:
            step = 1
        elif step == 0:
            raise ValueError("Step must not be zero")
        else:
            step = _ensure_int(step, 'step')

        return cls._simple_new(start, stop, step, name)

    @classmethod
    def from_range(cls, data, name=None, dtype=None, **kwargs):
        """ create RangeIndex from a range (py3), or xrange (py2) object """
        if not isinstance(data, range):
            raise TypeError(
                '{0}(...) must be called with object coercible to a '
                'range, {1} was passed'.format(cls.__name__, repr(data)))

        start, stop, step = get_range_parameters(data)
        return RangeIndex(start, stop, step, dtype=dtype, name=name, **kwargs)

    @classmethod
    def _simple_new(cls, start, stop=None, step=None, name=None,
                    dtype=None, **kwargs):
        result = object.__new__(cls)

        # handle passed None, non-integers
        if start is None and stop is None:
            # empty
            start, stop, step = 0, 0, 1

        if start is None or not is_integer(start):
            try:

                return RangeIndex(start, stop, step, name=name, **kwargs)
            except TypeError:
                return Index(start, stop, step, name=name, **kwargs)

        result._start = start
        result._stop = stop or 0
        result._step = step or 1
        result.name = name
        for k, v in compat.iteritems(kwargs):
            setattr(result, k, v)

        result._reset_identity()
        return result

    @staticmethod
    def _validate_dtype(dtype):
        """ require dtype to be None or int64 """
        if not (dtype is None or is_int64_dtype(dtype)):
            raise TypeError('Invalid to pass a non-int64 dtype to RangeIndex')

    @cache_readonly
    def _constructor(self):
        """ return the class to use for construction """
        return Int64Index

    @cache_readonly
    def _data(self):
        return np.arange(self._start, self._stop, self._step, dtype=np.int64)

    @cache_readonly
    def _int64index(self):
        return Int64Index(self._data, name=self.name, fastpath=True)

    def _get_data_as_items(self):
        """ return a list of tuples of start, stop, step """
        return [('start', self._start),
                ('stop', self._stop),
                ('step', self._step)]

    def __reduce__(self):
        d = self._get_attributes_dict()
        d.update(dict(self._get_data_as_items()))
        return ibase._new_Index, (self.__class__, d), None

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr, formatted_value)
        """
        attrs = self._get_data_as_items()
        if self.name is not None:
            attrs.append(('name', ibase.default_pprint(self.name)))
        return attrs

    def _format_data(self, name=None):
        # we are formatting thru the attributes
        return None

    @cache_readonly
    def nbytes(self):
        """
        Return the number of bytes in the underlying data
        On implementations where this is undetermined (PyPy)
        assume 24 bytes for each value
        """
        return sum(getsizeof(getattr(self, v), 24) for v in
                   ['_start', '_stop', '_step'])

    def memory_usage(self, deep=False):
        """
        Memory usage of my values

        Parameters
        ----------
        deep : bool
            Introspect the data deeply, interrogate
            `object` dtypes for system-level memory consumption

        Returns
        -------
        bytes used

        Notes
        -----
        Memory usage does not include memory consumed by elements that
        are not components of the array if deep=False

        See Also
        --------
        numpy.ndarray.nbytes
        """
        return self.nbytes

    @property
    def dtype(self):
        return np.dtype(np.int64)

    @property
    def is_unique(self):
        """ return if the index has unique values """
        return True

    @cache_readonly
    def is_monotonic_increasing(self):
        return self._step > 0 or len(self) <= 1

    @cache_readonly
    def is_monotonic_decreasing(self):
        return self._step < 0 or len(self) <= 1

    @property
    def has_duplicates(self):
        return False

    def tolist(self):
        return lrange(self._start, self._stop, self._step)

    @Appender(_index_shared_docs['_shallow_copy'])
    def _shallow_copy(self, values=None, **kwargs):
        if values is None:
            return RangeIndex(name=self.name, fastpath=True,
                              **dict(self._get_data_as_items()))
        else:
            kwargs.setdefault('name', self.name)
            return self._int64index._shallow_copy(values, **kwargs)

    @Appender(ibase._index_shared_docs['copy'])
    def copy(self, name=None, deep=False, dtype=None, **kwargs):
        self._validate_dtype(dtype)
        if name is None:
            name = self.name
        return RangeIndex(name=name, fastpath=True,
                          **dict(self._get_data_as_items()))

    def _minmax(self, meth):
        no_steps = len(self) - 1
        if no_steps == -1:
            return np.nan
        elif ((meth == 'min' and self._step > 0) or
              (meth == 'max' and self._step < 0)):
            return self._start

        return self._start + self._step * no_steps

    def min(self):
        """The minimum value of the RangeIndex"""
        return self._minmax('min')

    def max(self):
        """The maximum value of the RangeIndex"""
        return self._minmax('max')

    def argsort(self, *args, **kwargs):
        """
        Returns the indices that would sort the index and its
        underlying data.

        Returns
        -------
        argsorted : numpy array

        See also
        --------
        numpy.ndarray.argsort
        """
        nv.validate_argsort(args, kwargs)

        if self._step > 0:
            return np.arange(len(self))
        else:
            return np.arange(len(self) - 1, -1, -1)

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if isinstance(other, RangeIndex):
            ls = len(self)
            lo = len(other)
            return (ls == lo == 0 or
                    ls == lo == 1 and
                    self._start == other._start or
                    ls == lo and
                    self._start == other._start and
                    self._step == other._step)

        return super(RangeIndex, self).equals(other)

    def intersection(self, other):
        """
        Form the intersection of two Index objects. Sortedness of the result is
        not guaranteed

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        intersection : Index
        """
        if not isinstance(other, RangeIndex):
            return super(RangeIndex, self).intersection(other)

        if not len(self) or not len(other):
            return RangeIndex._simple_new(None)

        first = self[::-1] if self._step < 0 else self
        second = other[::-1] if other._step < 0 else other

        # check whether intervals intersect
        # deals with in- and decreasing ranges
        int_low = max(first._start, second._start)
        int_high = min(first._stop, second._stop)
        if int_high <= int_low:
            return RangeIndex._simple_new(None)

        # Method hint: linear Diophantine equation
        # solve intersection problem
        # performance hint: for identical step sizes, could use
        # cheaper alternative
        gcd, s, t = first._extended_gcd(first._step, second._step)

        # check whether element sets intersect
        if (first._start - second._start) % gcd:
            return RangeIndex._simple_new(None)

        # calculate parameters for the RangeIndex describing the
        # intersection disregarding the lower bounds
        tmp_start = first._start + (second._start - first._start) * \
            first._step // gcd * s
        new_step = first._step * second._step // gcd
        new_index = RangeIndex(tmp_start, int_high, new_step, fastpath=True)

        # adjust index to limiting interval
        new_index._start = new_index._min_fitting_element(int_low)

        if (self._step < 0 and other._step < 0) is not (new_index._step < 0):
            new_index = new_index[::-1]
        return new_index

    def _min_fitting_element(self, lower_limit):
        """Returns the smallest element greater than or equal to the limit"""
        no_steps = -(-(lower_limit - self._start) // abs(self._step))
        return self._start + abs(self._step) * no_steps

    def _max_fitting_element(self, upper_limit):
        """Returns the largest element smaller than or equal to the limit"""
        no_steps = (upper_limit - self._start) // abs(self._step)
        return self._start + abs(self._step) * no_steps

    def _extended_gcd(self, a, b):
        """
        Extended Euclidean algorithms to solve Bezout's identity:
           a*x + b*y = gcd(x, y)
        Finds one particular solution for x, y: s, t
        Returns: gcd, s, t
        """
        s, old_s = 0, 1
        t, old_t = 1, 0
        r, old_r = b, a
        while r:
            quotient = old_r // r
            old_r, r = r, old_r - quotient * r
            old_s, s = s, old_s - quotient * s
            old_t, t = t, old_t - quotient * t
        return old_r, old_s, old_t

    def union(self, other):
        """
        Form the union of two Index objects and sorts if possible

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        union : Index
        """
        self._assert_can_do_setop(other)
        if len(other) == 0 or self.equals(other):
            return self
        if len(self) == 0:
            return other
        if isinstance(other, RangeIndex):
            start_s, step_s = self._start, self._step
            end_s = self._start + self._step * (len(self) - 1)
            start_o, step_o = other._start, other._step
            end_o = other._start + other._step * (len(other) - 1)
            if self._step < 0:
                start_s, step_s, end_s = end_s, -step_s, start_s
            if other._step < 0:
                start_o, step_o, end_o = end_o, -step_o, start_o
            if len(self) == 1 and len(other) == 1:
                step_s = step_o = abs(self._start - other._start)
            elif len(self) == 1:
                step_s = step_o
            elif len(other) == 1:
                step_o = step_s
            start_r = min(start_s, start_o)
            end_r = max(end_s, end_o)
            if step_o == step_s:
                if ((start_s - start_o) % step_s == 0 and
                        (start_s - end_o) <= step_s and
                        (start_o - end_s) <= step_s):
                    return RangeIndex(start_r, end_r + step_s, step_s)
                if ((step_s % 2 == 0) and
                        (abs(start_s - start_o) <= step_s / 2) and
                        (abs(end_s - end_o) <= step_s / 2)):
                    return RangeIndex(start_r, end_r + step_s / 2, step_s / 2)
            elif step_o % step_s == 0:
                if ((start_o - start_s) % step_s == 0 and
                        (start_o + step_s >= start_s) and
                        (end_o - step_s <= end_s)):
                    return RangeIndex(start_r, end_r + step_s, step_s)
            elif step_s % step_o == 0:
                if ((start_s - start_o) % step_o == 0 and
                        (start_s + step_o >= start_o) and
                        (end_s - step_o <= end_o)):
                    return RangeIndex(start_r, end_r + step_o, step_o)

        return self._int64index.union(other)

    @Appender(_index_shared_docs['join'])
    def join(self, other, how='left', level=None, return_indexers=False,
             sort=False):
        if how == 'outer' and self is not other:
            # note: could return RangeIndex in more circumstances
            return self._int64index.join(other, how, level, return_indexers,
                                         sort)

        return super(RangeIndex, self).join(other, how, level, return_indexers,
                                            sort)

    def _concat_same_dtype(self, indexes, name):
        return _concat._concat_rangeindex_same_dtype(indexes).rename(name)

    def __len__(self):
        """
        return the length of the RangeIndex
        """
        return max(0, -(-(self._stop - self._start) // self._step))

    @property
    def size(self):
        return len(self)

    def __getitem__(self, key):
        """
        Conserve RangeIndex type for scalar and slice keys.
        """
        super_getitem = super(RangeIndex, self).__getitem__

        if is_scalar(key):
            n = int(key)
            if n != key:
                return super_getitem(key)
            if n < 0:
                n = len(self) + key
            if n < 0 or n > len(self) - 1:
                raise IndexError("index {key} is out of bounds for axis 0 "
                                 "with size {size}".format(key=key,
                                                           size=len(self)))
            return self._start + n * self._step

        if isinstance(key, slice):

            # This is basically PySlice_GetIndicesEx, but delegation to our
            # super routines if we don't have integers

            l = len(self)

            # complete missing slice information
            step = 1 if key.step is None else key.step
            if key.start is None:
                start = l - 1 if step < 0 else 0
            else:
                start = key.start

                if start < 0:
                    start += l
                if start < 0:
                    start = -1 if step < 0 else 0
                if start >= l:
                    start = l - 1 if step < 0 else l

            if key.stop is None:
                stop = -1 if step < 0 else l
            else:
                stop = key.stop

                if stop < 0:
                    stop += l
                if stop < 0:
                    stop = -1
                if stop > l:
                    stop = l

            # delegate non-integer slices
            if (start != int(start) or
                    stop != int(stop) or
                    step != int(step)):
                return super_getitem(key)

            # convert indexes to values
            start = self._start + self._step * start
            stop = self._start + self._step * stop
            step = self._step * step

            return RangeIndex(start, stop, step, name=self.name, fastpath=True)

        # fall back to Int64Index
        return super_getitem(key)

    def __floordiv__(self, other):
        if is_integer(other) and other != 0:
            if (len(self) == 0 or
                    self._start % other == 0 and
                    self._step % other == 0):
                start = self._start // other
                step = self._step // other
                stop = start + len(self) * step
                return RangeIndex(start, stop, step, name=self.name,
                                  fastpath=True)
            if len(self) == 1:
                start = self._start // other
                return RangeIndex(start, start + 1, 1, name=self.name,
                                  fastpath=True)
        return self._int64index // other

    @classmethod
    def _add_numeric_methods_binary(cls):
        """ add in numeric methods, specialized to RangeIndex """

        def _make_evaluate_binop(op, step=False):
            """
            Parameters
            ----------
            op : callable that accepts 2 parms
                perform the binary op
            step : callable, optional, default to False
                op to apply to the step parm if not None
                if False, use the existing step
            """

            def _evaluate_numeric_binop(self, other):
                if isinstance(other, ABCSeries):
                    return NotImplemented
                elif isinstance(other, ABCTimedeltaIndex):
                    # Defer to TimedeltaIndex implementation
                    return NotImplemented
                elif isinstance(other, (timedelta, np.timedelta64)):
                    # GH#19333 is_integer evaluated True on timedelta64,
                    # so we need to catch these explicitly
                    return op(self._int64index, other)

                other = self._validate_for_numeric_binop(other, op)
                attrs = self._get_attributes_dict()
                attrs = self._maybe_update_attributes(attrs)

                left, right = self, other

                try:
                    # apply if we have an override
                    if step:
                        with np.errstate(all='ignore'):
                            rstep = step(left._step, right)

                        # we don't have a representable op
                        # so return a base index
                        if not is_integer(rstep) or not rstep:
                            raise ValueError

                    else:
                        rstep = left._step

                    with np.errstate(all='ignore'):
                        rstart = op(left._start, right)
                        rstop = op(left._stop, right)

                    result = RangeIndex(rstart,
                                        rstop,
                                        rstep,
                                        **attrs)

                    # for compat with numpy / Int64Index
                    # even if we can represent as a RangeIndex, return
                    # as a Float64Index if we have float-like descriptors
                    if not all(is_integer(x) for x in
                               [rstart, rstop, rstep]):
                        result = result.astype('float64')

                    return result

                except (ValueError, TypeError, ZeroDivisionError):
                    # Defer to Int64Index implementation
                    return op(self._int64index, other)
                    # TODO: Do attrs get handled reliably?

            return _evaluate_numeric_binop

        cls.__add__ = _make_evaluate_binop(operator.add)
        cls.__radd__ = _make_evaluate_binop(ops.radd)
        cls.__sub__ = _make_evaluate_binop(operator.sub)
        cls.__rsub__ = _make_evaluate_binop(ops.rsub)
        cls.__mul__ = _make_evaluate_binop(operator.mul, step=operator.mul)
        cls.__rmul__ = _make_evaluate_binop(ops.rmul, step=ops.rmul)
        cls.__truediv__ = _make_evaluate_binop(operator.truediv,
                                               step=operator.truediv)
        cls.__rtruediv__ = _make_evaluate_binop(ops.rtruediv,
                                                step=ops.rtruediv)
        if not compat.PY3:
            cls.__div__ = _make_evaluate_binop(operator.div, step=operator.div)
            cls.__rdiv__ = _make_evaluate_binop(ops.rdiv, step=ops.rdiv)


RangeIndex._add_numeric_methods()
RangeIndex._add_logical_methods()
