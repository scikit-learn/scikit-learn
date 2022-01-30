"""
Base and utility classes for tseries type pandas objects.
"""
from __future__ import annotations

from datetime import datetime
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Sequence,
    TypeVar,
    cast,
    final,
)
import warnings

import numpy as np

from pandas._libs import (
    NaT,
    Timedelta,
    lib,
)
from pandas._libs.tslibs import (
    BaseOffset,
    Resolution,
    Tick,
    parsing,
    to_offset,
)
from pandas.compat.numpy import function as nv
from pandas.util._decorators import (
    Appender,
    cache_readonly,
    doc,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_dtype_equal,
    is_integer,
    is_list_like,
)
from pandas.core.dtypes.concat import concat_compat

from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    PeriodArray,
    TimedeltaArray,
)
from pandas.core.arrays.datetimelike import DatetimeLikeArrayMixin
import pandas.core.common as com
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import (
    Index,
    _index_shared_docs,
)
from pandas.core.indexes.extension import (
    NDArrayBackedExtensionIndex,
    inherit_names,
)
from pandas.core.indexes.range import RangeIndex
from pandas.core.tools.timedeltas import to_timedelta

if TYPE_CHECKING:
    from pandas import CategoricalIndex

_index_doc_kwargs = dict(ibase._index_doc_kwargs)

_T = TypeVar("_T", bound="DatetimeIndexOpsMixin")
_TDT = TypeVar("_TDT", bound="DatetimeTimedeltaMixin")


@inherit_names(
    ["inferred_freq", "_resolution_obj", "resolution"],
    DatetimeLikeArrayMixin,
    cache=True,
)
@inherit_names(["mean", "asi8", "freq", "freqstr"], DatetimeLikeArrayMixin)
class DatetimeIndexOpsMixin(NDArrayBackedExtensionIndex):
    """
    Common ops mixin to support a unified interface datetimelike Index.
    """

    _is_numeric_dtype = False
    _can_hold_strings = False
    _data: DatetimeArray | TimedeltaArray | PeriodArray
    freq: BaseOffset | None
    freqstr: str | None
    _resolution_obj: Resolution

    # error: "Callable[[Any], Any]" has no attribute "fget"
    hasnans = cast(
        bool,
        cache_readonly(
            DatetimeLikeArrayMixin._hasna.fget  # type: ignore[attr-defined]
        ),
    )

    @property
    def _is_all_dates(self) -> bool:
        return True

    # ------------------------------------------------------------------------

    def equals(self, other: Any) -> bool:
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if not isinstance(other, Index):
            return False
        elif other.dtype.kind in ["f", "i", "u", "c"]:
            return False
        elif not isinstance(other, type(self)):
            should_try = False
            inferable = self._data._infer_matches
            if other.dtype == object:
                should_try = other.inferred_type in inferable
            elif is_categorical_dtype(other.dtype):
                other = cast("CategoricalIndex", other)
                should_try = other.categories.inferred_type in inferable

            if should_try:
                try:
                    other = type(self)(other)
                except (ValueError, TypeError, OverflowError):
                    # e.g.
                    #  ValueError -> cannot parse str entry, or OutOfBoundsDatetime
                    #  TypeError  -> trying to convert IntervalIndex to DatetimeIndex
                    #  OverflowError -> Index([very_large_timedeltas])
                    return False

        if not is_dtype_equal(self.dtype, other.dtype):
            # have different timezone
            return False

        return np.array_equal(self.asi8, other.asi8)

    @Appender(Index.__contains__.__doc__)
    def __contains__(self, key: Any) -> bool:
        hash(key)
        try:
            self.get_loc(key)
        except (KeyError, TypeError, ValueError):
            return False
        return True

    _can_hold_na = True

    def _convert_tolerance(self, tolerance, target):
        tolerance = np.asarray(to_timedelta(tolerance).to_numpy())
        return super()._convert_tolerance(tolerance, target)

    # --------------------------------------------------------------------
    # Rendering Methods

    def format(
        self,
        name: bool = False,
        formatter: Callable | None = None,
        na_rep: str = "NaT",
        date_format: str | None = None,
    ) -> list[str]:
        """
        Render a string representation of the Index.
        """
        header = []
        if name:
            header.append(
                ibase.pprint_thing(self.name, escape_chars=("\t", "\r", "\n"))
                if self.name is not None
                else ""
            )

        if formatter is not None:
            return header + list(self.map(formatter))

        return self._format_with_header(header, na_rep=na_rep, date_format=date_format)

    def _format_with_header(
        self, header: list[str], na_rep: str = "NaT", date_format: str | None = None
    ) -> list[str]:
        # matches base class except for whitespace padding and date_format
        return header + list(
            self._format_native_types(na_rep=na_rep, date_format=date_format)
        )

    @property
    def _formatter_func(self):
        return self._data._formatter()

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value).
        """
        attrs = super()._format_attrs()
        for attrib in self._attributes:
            # iterating over _attributes prevents us from doing this for PeriodIndex
            if attrib == "freq":
                freq = self.freqstr
                if freq is not None:
                    freq = repr(freq)  # e.g. D -> 'D'
                attrs.append(("freq", freq))
        return attrs

    @Appender(Index._summary.__doc__)
    def _summary(self, name=None) -> str:
        result = super()._summary(name=name)
        if self.freq:
            result += f"\nFreq: {self.freqstr}"

        return result

    # --------------------------------------------------------------------
    # Indexing Methods

    def _can_partial_date_slice(self, reso: Resolution) -> bool:
        raise NotImplementedError

    def _parsed_string_to_bounds(self, reso: Resolution, parsed):
        raise NotImplementedError

    def _parse_with_reso(self, label: str):
        # overridden by TimedeltaIndex
        parsed, reso_str = parsing.parse_time_string(label, self.freq)
        reso = Resolution.from_attrname(reso_str)
        return parsed, reso

    def _get_string_slice(self, key: str):
        parsed, reso = self._parse_with_reso(key)
        try:
            return self._partial_date_slice(reso, parsed)
        except KeyError as err:
            raise KeyError(key) from err

    @final
    def _partial_date_slice(
        self,
        reso: Resolution,
        parsed: datetime,
    ):
        """
        Parameters
        ----------
        reso : Resolution
        parsed : datetime

        Returns
        -------
        slice or ndarray[intp]
        """
        if not self._can_partial_date_slice(reso):
            raise ValueError

        t1, t2 = self._parsed_string_to_bounds(reso, parsed)
        vals = self._data._ndarray
        unbox = self._data._unbox

        if self.is_monotonic_increasing:

            if len(self) and (
                (t1 < self[0] and t2 < self[0]) or (t1 > self[-1] and t2 > self[-1])
            ):
                # we are out of range
                raise KeyError

            # TODO: does this depend on being monotonic _increasing_?

            # a monotonic (sorted) series can be sliced
            left = vals.searchsorted(unbox(t1), side="left")
            right = vals.searchsorted(unbox(t2), side="right")
            return slice(left, right)

        else:
            lhs_mask = vals >= unbox(t1)
            rhs_mask = vals <= unbox(t2)

            # try to find the dates
            return (lhs_mask & rhs_mask).nonzero()[0]

    def _maybe_cast_slice_bound(self, label, side: str, kind=lib.no_default):
        """
        If label is a string, cast it to scalar type according to resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : {'loc', 'getitem'} or None

        Returns
        -------
        label : object

        Notes
        -----
        Value of `side` parameter should be validated in caller.
        """
        assert kind in ["loc", "getitem", None, lib.no_default]
        self._deprecated_arg(kind, "kind", "_maybe_cast_slice_bound")

        if isinstance(label, str):
            try:
                parsed, reso = self._parse_with_reso(label)
            except ValueError as err:
                # DTI -> parsing.DateParseError
                # TDI -> 'unit abbreviation w/o a number'
                # PI -> string cannot be parsed as datetime-like
                raise self._invalid_indexer("slice", label) from err

            lower, upper = self._parsed_string_to_bounds(reso, parsed)
            return lower if side == "left" else upper
        elif not isinstance(label, self._data._recognized_scalars):
            raise self._invalid_indexer("slice", label)

        return label

    # --------------------------------------------------------------------
    # Arithmetic Methods

    def shift(self: _T, periods: int = 1, freq=None) -> _T:
        """
        Shift index by desired number of time frequency increments.

        This method is for shifting the values of datetime-like indexes
        by a specified time increment a given number of times.

        Parameters
        ----------
        periods : int, default 1
            Number of periods (or increments) to shift by,
            can be positive or negative.
        freq : pandas.DateOffset, pandas.Timedelta or string, optional
            Frequency increment to shift by.
            If None, the index is shifted by its own `freq` attribute.
            Offset aliases are valid strings, e.g., 'D', 'W', 'M' etc.

        Returns
        -------
        pandas.DatetimeIndex
            Shifted index.

        See Also
        --------
        Index.shift : Shift values of Index.
        PeriodIndex.shift : Shift values of PeriodIndex.
        """
        arr = self._data.view()
        arr._freq = self.freq
        result = arr._time_shift(periods, freq=freq)
        return type(self)._simple_new(result, name=self.name)

    # --------------------------------------------------------------------

    @doc(Index._maybe_cast_listlike_indexer)
    def _maybe_cast_listlike_indexer(self, keyarr):
        try:
            res = self._data._validate_listlike(keyarr, allow_object=True)
        except (ValueError, TypeError):
            if not isinstance(keyarr, ExtensionArray):
                # e.g. we don't want to cast DTA to ndarray[object]
                res = com.asarray_tuplesafe(keyarr)
                # TODO: com.asarray_tuplesafe shouldn't cast e.g. DatetimeArray
            else:
                res = keyarr
        return Index(res, dtype=res.dtype)


class DatetimeTimedeltaMixin(DatetimeIndexOpsMixin):
    """
    Mixin class for methods shared by DatetimeIndex and TimedeltaIndex,
    but not PeriodIndex
    """

    _data: DatetimeArray | TimedeltaArray
    _comparables = ["name", "freq"]
    _attributes = ["name", "freq"]

    # Compat for frequency inference, see GH#23789
    _is_monotonic_increasing = Index.is_monotonic_increasing
    _is_monotonic_decreasing = Index.is_monotonic_decreasing
    _is_unique = Index.is_unique

    _join_precedence = 10

    def _with_freq(self, freq):
        arr = self._data._with_freq(freq)
        return type(self)._simple_new(arr, name=self._name)

    def is_type_compatible(self, kind: str) -> bool:
        warnings.warn(
            f"{type(self).__name__}.is_type_compatible is deprecated and will be "
            "removed in a future version.",
            FutureWarning,
            stacklevel=find_stack_level(),
        )
        return kind in self._data._infer_matches

    @property
    def values(self) -> np.ndarray:
        # NB: For Datetime64TZ this is lossy
        return self._data._ndarray

    # --------------------------------------------------------------------
    # Set Operation Methods

    @cache_readonly
    def _as_range_index(self) -> RangeIndex:
        # Convert our i8 representations to RangeIndex
        # Caller is responsible for checking isinstance(self.freq, Tick)
        freq = cast(Tick, self.freq)
        tick = freq.delta.value
        rng = range(self[0].value, self[-1].value + tick, tick)
        return RangeIndex(rng)

    def _can_range_setop(self, other):
        return isinstance(self.freq, Tick) and isinstance(other.freq, Tick)

    def _wrap_range_setop(self, other, res_i8):
        new_freq = None
        if not len(res_i8):
            # RangeIndex defaults to step=1, which we don't want.
            new_freq = self.freq
        elif isinstance(res_i8, RangeIndex):
            new_freq = to_offset(Timedelta(res_i8.step))
            res_i8 = res_i8

        # TODO: we cannot just do
        #  type(self._data)(res_i8.values, dtype=self.dtype, freq=new_freq)
        # because test_setops_preserve_freq fails with _validate_frequency raising.
        # This raising is incorrect, as 'on_freq' is incorrect. This will
        # be fixed by GH#41493
        res_values = res_i8.values.view(self._data._ndarray.dtype)
        result = type(self._data)._simple_new(
            res_values, dtype=self.dtype, freq=new_freq
        )
        return self._wrap_setop_result(other, result)

    def _range_intersect(self, other, sort):
        # Dispatch to RangeIndex intersection logic.
        left = self._as_range_index
        right = other._as_range_index
        res_i8 = left.intersection(right, sort=sort)
        return self._wrap_range_setop(other, res_i8)

    def _range_union(self, other, sort):
        # Dispatch to RangeIndex union logic.
        left = self._as_range_index
        right = other._as_range_index
        res_i8 = left.union(right, sort=sort)
        return self._wrap_range_setop(other, res_i8)

    def _intersection(self, other: Index, sort=False) -> Index:
        """
        intersection specialized to the case with matching dtypes and both non-empty.
        """
        other = cast("DatetimeTimedeltaMixin", other)

        if self._can_range_setop(other):
            return self._range_intersect(other, sort=sort)

        if not self._can_fast_intersect(other):
            result = Index._intersection(self, other, sort=sort)
            # We need to invalidate the freq because Index._intersection
            #  uses _shallow_copy on a view of self._data, which will preserve
            #  self.freq if we're not careful.
            # At this point we should have result.dtype == self.dtype
            #  and type(result) is type(self._data)
            result = self._wrap_setop_result(other, result)
            return result._with_freq(None)._with_freq("infer")

        else:
            return self._fast_intersect(other, sort)

    def _fast_intersect(self, other, sort):
        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        # after sorting, the intersection always starts with the right index
        # and ends with the index of which the last elements is smallest
        end = min(left[-1], right[-1])
        start = right[0]

        if end < start:
            result = self[:0]
        else:
            lslice = slice(*left.slice_locs(start, end))
            result = left._values[lslice]

        return result

    def _can_fast_intersect(self: _T, other: _T) -> bool:
        # Note: we only get here with len(self) > 0 and len(other) > 0
        if self.freq is None:
            return False

        elif other.freq != self.freq:
            return False

        elif not self.is_monotonic_increasing:
            # Because freq is not None, we must then be monotonic decreasing
            return False

        # this along with matching freqs ensure that we "line up",
        #  so intersection will preserve freq
        # Note we are assuming away Ticks, as those go through _range_intersect
        # GH#42104
        return self.freq.n == 1

    def _can_fast_union(self: _T, other: _T) -> bool:
        # Assumes that type(self) == type(other), as per the annotation
        # The ability to fast_union also implies that `freq` should be
        #  retained on union.
        freq = self.freq

        if freq is None or freq != other.freq:
            return False

        if not self.is_monotonic_increasing:
            # Because freq is not None, we must then be monotonic decreasing
            # TODO: do union on the reversed indexes?
            return False

        if len(self) == 0 or len(other) == 0:
            # only reached via union_many
            return True

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        right_start = right[0]
        left_end = left[-1]

        # Only need to "adjoin", not overlap
        return (right_start == left_end + freq) or right_start in left

    def _fast_union(self: _TDT, other: _TDT, sort=None) -> _TDT:
        # Caller is responsible for ensuring self and other are non-empty

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        elif sort is False:
            # TDIs are not in the "correct" order and we don't want
            #  to sort but want to remove overlaps
            left, right = self, other
            left_start = left[0]
            loc = right.searchsorted(left_start, side="left")
            right_chunk = right._values[:loc]
            dates = concat_compat((left._values, right_chunk))
            result = type(self)._simple_new(dates, name=self.name)
            return result
        else:
            left, right = other, self

        left_end = left[-1]
        right_end = right[-1]

        # concatenate
        if left_end < right_end:
            loc = right.searchsorted(left_end, side="right")
            right_chunk = right._values[loc:]
            dates = concat_compat([left._values, right_chunk])
            # The can_fast_union check ensures that the result.freq
            #  should match self.freq
            dates = type(self._data)(dates, freq=self.freq)
            result = type(self)._simple_new(dates)
            return result
        else:
            return left

    def _union(self, other, sort):
        # We are called by `union`, which is responsible for this validation
        assert isinstance(other, type(self))
        assert self.dtype == other.dtype

        if self._can_range_setop(other):
            return self._range_union(other, sort=sort)

        if self._can_fast_union(other):
            result = self._fast_union(other, sort=sort)
            # in the case with sort=None, the _can_fast_union check ensures
            #  that result.freq == self.freq
            return result
        else:
            return super()._union(other, sort)._with_freq("infer")

    # --------------------------------------------------------------------
    # Join Methods

    def _get_join_freq(self, other):
        """
        Get the freq to attach to the result of a join operation.
        """
        freq = None
        if self._can_fast_union(other):
            freq = self.freq
        return freq

    def _wrap_joined_index(self, joined, other):
        assert other.dtype == self.dtype, (other.dtype, self.dtype)
        result = super()._wrap_joined_index(joined, other)
        result._data._freq = self._get_join_freq(other)
        return result

    def _get_engine_target(self) -> np.ndarray:
        # engine methods and libjoin methods need dt64/td64 values cast to i8
        return self._data._ndarray.view("i8")

    def _from_join_target(self, result: np.ndarray):
        # view e.g. i8 back to M8[ns]
        result = result.view(self._data._ndarray.dtype)
        return self._data._from_backing_data(result)

    # --------------------------------------------------------------------
    # List-like Methods

    def _get_delete_freq(self, loc: int | slice | Sequence[int]):
        """
        Find the `freq` for self.delete(loc).
        """
        freq = None
        if self.freq is not None:
            if is_integer(loc):
                if loc in (0, -len(self), -1, len(self) - 1):
                    freq = self.freq
            else:
                if is_list_like(loc):
                    # error: Incompatible types in assignment (expression has
                    # type "Union[slice, ndarray]", variable has type
                    # "Union[int, slice, Sequence[int]]")
                    loc = lib.maybe_indices_to_slice(  # type: ignore[assignment]
                        np.asarray(loc, dtype=np.intp), len(self)
                    )
                if isinstance(loc, slice) and loc.step in (1, None):
                    if loc.start in (0, None) or loc.stop in (len(self), None):
                        freq = self.freq
        return freq

    def _get_insert_freq(self, loc: int, item):
        """
        Find the `freq` for self.insert(loc, item).
        """
        value = self._data._validate_scalar(item)
        item = self._data._box_func(value)

        freq = None
        if self.freq is not None:
            # freq can be preserved on edge cases
            if self.size:
                if item is NaT:
                    pass
                elif (loc == 0 or loc == -len(self)) and item + self.freq == self[0]:
                    freq = self.freq
                elif (loc == len(self)) and item - self.freq == self[-1]:
                    freq = self.freq
            else:
                # Adding a single item to an empty index may preserve freq
                if isinstance(self.freq, Tick):
                    # all TimedeltaIndex cases go through here; is_on_offset
                    #  would raise TypeError
                    freq = self.freq
                elif self.freq.is_on_offset(item):
                    freq = self.freq
        return freq

    @doc(NDArrayBackedExtensionIndex.delete)
    def delete(self, loc):
        result = super().delete(loc)
        result._data._freq = self._get_delete_freq(loc)
        return result

    @doc(NDArrayBackedExtensionIndex.insert)
    def insert(self, loc: int, item):
        result = super().insert(loc, item)
        if isinstance(result, type(self)):
            # i.e. parent class method did not cast
            result._data._freq = self._get_insert_freq(loc, item)
        return result

    # --------------------------------------------------------------------
    # NDArray-Like Methods

    @Appender(_index_shared_docs["take"] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True, fill_value=None, **kwargs):
        nv.validate_take((), kwargs)
        indices = np.asarray(indices, dtype=np.intp)

        result = NDArrayBackedExtensionIndex.take(
            self, indices, axis, allow_fill, fill_value, **kwargs
        )

        maybe_slice = lib.maybe_indices_to_slice(indices, len(self))
        if isinstance(maybe_slice, slice):
            freq = self._data._get_getitem_freq(maybe_slice)
            result._data._freq = freq
        return result
