"""
Provide classes to perform the groupby aggregate operations.

These are not exposed to the user and provide implementations of the grouping
operations, primarily in cython. These classes (BaseGrouper and BinGrouper)
are contained *in* the SeriesGroupBy and DataFrameGroupBy objects.
"""
from __future__ import annotations

import collections
import functools
from typing import (
    Callable,
    Generic,
    Hashable,
    Iterator,
    Sequence,
    final,
    overload,
)

import numpy as np

from pandas._libs import (
    NaT,
    lib,
)
import pandas._libs.groupby as libgroupby
import pandas._libs.reduction as libreduction
from pandas._typing import (
    ArrayLike,
    DtypeObj,
    NDFrameT,
    Shape,
    npt,
)
from pandas.errors import AbstractMethodError
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.cast import (
    maybe_cast_pointwise_result,
    maybe_downcast_to_dtype,
)
from pandas.core.dtypes.common import (
    ensure_float64,
    ensure_int64,
    ensure_platform_int,
    is_1d_only_ea_obj,
    is_bool_dtype,
    is_categorical_dtype,
    is_complex_dtype,
    is_datetime64_any_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_sparse,
    is_timedelta64_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.missing import (
    isna,
    maybe_fill,
)

from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    PeriodArray,
    TimedeltaArray,
)
from pandas.core.arrays.boolean import BooleanDtype
from pandas.core.arrays.floating import (
    Float64Dtype,
    FloatingDtype,
)
from pandas.core.arrays.integer import (
    Int64Dtype,
    _IntegerDtype,
)
from pandas.core.arrays.masked import (
    BaseMaskedArray,
    BaseMaskedDtype,
)
from pandas.core.arrays.string_ import StringDtype
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.groupby import grouper
from pandas.core.indexes.api import (
    CategoricalIndex,
    Index,
    MultiIndex,
    ensure_index,
)
from pandas.core.series import Series
from pandas.core.sorting import (
    compress_group_index,
    decons_obs_group_ids,
    get_flattened_list,
    get_group_index,
    get_group_index_sorter,
    get_indexer_dict,
)


class WrappedCythonOp:
    """
    Dispatch logic for functions defined in _libs.groupby
    """

    # Functions for which we do _not_ attempt to cast the cython result
    #  back to the original dtype.
    cast_blocklist = frozenset(["rank", "count", "size", "idxmin", "idxmax"])

    def __init__(self, kind: str, how: str):
        self.kind = kind
        self.how = how

    _CYTHON_FUNCTIONS = {
        "aggregate": {
            "add": "group_add",
            "prod": "group_prod",
            "min": "group_min",
            "max": "group_max",
            "mean": "group_mean",
            "median": "group_median",
            "var": "group_var",
            "first": "group_nth",
            "last": "group_last",
            "ohlc": "group_ohlc",
        },
        "transform": {
            "cumprod": "group_cumprod",
            "cumsum": "group_cumsum",
            "cummin": "group_cummin",
            "cummax": "group_cummax",
            "rank": "group_rank",
        },
    }

    _MASKED_CYTHON_FUNCTIONS = {"cummin", "cummax", "min", "max"}

    _cython_arity = {"ohlc": 4}  # OHLC

    # Note: we make this a classmethod and pass kind+how so that caching
    #  works at the class level and not the instance level
    @classmethod
    @functools.lru_cache(maxsize=None)
    def _get_cython_function(
        cls, kind: str, how: str, dtype: np.dtype, is_numeric: bool
    ):

        dtype_str = dtype.name
        ftype = cls._CYTHON_FUNCTIONS[kind][how]

        # see if there is a fused-type version of function
        # only valid for numeric
        f = getattr(libgroupby, ftype)
        if is_numeric:
            return f
        elif dtype == object:
            if "object" not in f.__signatures__:
                # raise NotImplementedError here rather than TypeError later
                raise NotImplementedError(
                    f"function is not implemented for this dtype: "
                    f"[how->{how},dtype->{dtype_str}]"
                )
            return f

    def get_cython_func_and_vals(self, values: np.ndarray, is_numeric: bool):
        """
        Find the appropriate cython function, casting if necessary.

        Parameters
        ----------
        values : np.ndarray
        is_numeric : bool

        Returns
        -------
        func : callable
        values : np.ndarray
        """
        how = self.how
        kind = self.kind

        if how in ["median", "cumprod"]:
            # these two only have float64 implementations
            if is_numeric:
                values = ensure_float64(values)
            else:
                raise NotImplementedError(
                    f"function is not implemented for this dtype: "
                    f"[how->{how},dtype->{values.dtype.name}]"
                )
            func = getattr(libgroupby, f"group_{how}_float64")
            return func, values

        func = self._get_cython_function(kind, how, values.dtype, is_numeric)

        if values.dtype.kind in ["i", "u"]:
            if how in ["add", "var", "prod", "mean", "ohlc"]:
                # result may still include NaN, so we have to cast
                values = ensure_float64(values)

        return func, values

    def _disallow_invalid_ops(self, dtype: DtypeObj, is_numeric: bool = False):
        """
        Check if we can do this operation with our cython functions.

        Raises
        ------
        NotImplementedError
            This is either not a valid function for this dtype, or
            valid but not implemented in cython.
        """
        how = self.how

        if is_numeric:
            # never an invalid op for those dtypes, so return early as fastpath
            return

        if is_categorical_dtype(dtype):
            # NotImplementedError for methods that can fall back to a
            #  non-cython implementation.
            if how in ["add", "prod", "cumsum", "cumprod"]:
                raise TypeError(f"{dtype} type does not support {how} operations")
            raise NotImplementedError(f"{dtype} dtype not supported")

        elif is_sparse(dtype):
            # categoricals are only 1d, so we
            #  are not setup for dim transforming
            raise NotImplementedError(f"{dtype} dtype not supported")
        elif is_datetime64_any_dtype(dtype):
            # we raise NotImplemented if this is an invalid operation
            #  entirely, e.g. adding datetimes
            if how in ["add", "prod", "cumsum", "cumprod"]:
                raise TypeError(f"datetime64 type does not support {how} operations")
        elif is_timedelta64_dtype(dtype):
            if how in ["prod", "cumprod"]:
                raise TypeError(f"timedelta64 type does not support {how} operations")

    def _get_output_shape(self, ngroups: int, values: np.ndarray) -> Shape:
        how = self.how
        kind = self.kind

        arity = self._cython_arity.get(how, 1)

        out_shape: Shape
        if how == "ohlc":
            out_shape = (ngroups, 4)
        elif arity > 1:
            raise NotImplementedError(
                "arity of more than 1 is not supported for the 'how' argument"
            )
        elif kind == "transform":
            out_shape = values.shape
        else:
            out_shape = (ngroups,) + values.shape[1:]
        return out_shape

    def get_out_dtype(self, dtype: np.dtype) -> np.dtype:
        how = self.how

        if how == "rank":
            out_dtype = "float64"
        else:
            if is_numeric_dtype(dtype):
                out_dtype = f"{dtype.kind}{dtype.itemsize}"
            else:
                out_dtype = "object"
        return np.dtype(out_dtype)

    @overload
    def _get_result_dtype(self, dtype: np.dtype) -> np.dtype:
        ...  # pragma: no cover

    @overload
    def _get_result_dtype(self, dtype: ExtensionDtype) -> ExtensionDtype:
        ...  # pragma: no cover

    def _get_result_dtype(self, dtype: DtypeObj) -> DtypeObj:
        """
        Get the desired dtype of a result based on the
        input dtype and how it was computed.

        Parameters
        ----------
        dtype : np.dtype or ExtensionDtype
            Input dtype.

        Returns
        -------
        np.dtype or ExtensionDtype
            The desired dtype of the result.
        """
        how = self.how

        if how in ["add", "cumsum", "sum", "prod"]:
            if dtype == np.dtype(bool):
                return np.dtype(np.int64)
            elif isinstance(dtype, (BooleanDtype, _IntegerDtype)):
                return Int64Dtype()
        elif how in ["mean", "median", "var"]:
            if isinstance(dtype, (BooleanDtype, _IntegerDtype)):
                return Float64Dtype()
            elif is_float_dtype(dtype) or is_complex_dtype(dtype):
                return dtype
            elif is_numeric_dtype(dtype):
                return np.dtype(np.float64)
        return dtype

    def uses_mask(self) -> bool:
        return self.how in self._MASKED_CYTHON_FUNCTIONS

    @final
    def _ea_wrap_cython_operation(
        self,
        values: ExtensionArray,
        min_count: int,
        ngroups: int,
        comp_ids: np.ndarray,
        **kwargs,
    ) -> ArrayLike:
        """
        If we have an ExtensionArray, unwrap, call _cython_operation, and
        re-wrap if appropriate.
        """
        # TODO: general case implementation overridable by EAs.
        if isinstance(values, BaseMaskedArray) and self.uses_mask():
            return self._masked_ea_wrap_cython_operation(
                values,
                min_count=min_count,
                ngroups=ngroups,
                comp_ids=comp_ids,
                **kwargs,
            )

        if isinstance(values, (DatetimeArray, PeriodArray, TimedeltaArray)):
            # All of the functions implemented here are ordinal, so we can
            #  operate on the tz-naive equivalents
            npvalues = values._ndarray.view("M8[ns]")
        elif isinstance(values.dtype, (BooleanDtype, _IntegerDtype)):
            # IntegerArray or BooleanArray
            npvalues = values.to_numpy("float64", na_value=np.nan)
        elif isinstance(values.dtype, FloatingDtype):
            # FloatingArray
            npvalues = values.to_numpy(values.dtype.numpy_dtype, na_value=np.nan)
        elif isinstance(values.dtype, StringDtype):
            # StringArray
            npvalues = values.to_numpy(object, na_value=np.nan)
        else:
            raise NotImplementedError(
                f"function is not implemented for this dtype: {values.dtype}"
            )

        res_values = self._cython_op_ndim_compat(
            npvalues,
            min_count=min_count,
            ngroups=ngroups,
            comp_ids=comp_ids,
            mask=None,
            **kwargs,
        )

        if self.how in ["rank"]:
            # i.e. how in WrappedCythonOp.cast_blocklist, since
            #  other cast_blocklist methods dont go through cython_operation
            return res_values

        return self._reconstruct_ea_result(values, res_values)

    def _reconstruct_ea_result(self, values, res_values):
        """
        Construct an ExtensionArray result from an ndarray result.
        """
        # TODO: allow EAs to override this logic

        if isinstance(
            values.dtype, (BooleanDtype, _IntegerDtype, FloatingDtype, StringDtype)
        ):
            dtype = self._get_result_dtype(values.dtype)
            cls = dtype.construct_array_type()
            return cls._from_sequence(res_values, dtype=dtype)

        elif needs_i8_conversion(values.dtype):
            i8values = res_values.view("i8")
            return type(values)(i8values, dtype=values.dtype)

        raise NotImplementedError

    @final
    def _masked_ea_wrap_cython_operation(
        self,
        values: BaseMaskedArray,
        min_count: int,
        ngroups: int,
        comp_ids: np.ndarray,
        **kwargs,
    ) -> BaseMaskedArray:
        """
        Equivalent of `_ea_wrap_cython_operation`, but optimized for masked EA's
        and cython algorithms which accept a mask.
        """
        orig_values = values

        # Copy to ensure input and result masks don't end up shared
        mask = values._mask.copy()
        result_mask = np.zeros(ngroups, dtype=bool)
        arr = values._data

        res_values = self._cython_op_ndim_compat(
            arr,
            min_count=min_count,
            ngroups=ngroups,
            comp_ids=comp_ids,
            mask=mask,
            result_mask=result_mask,
            **kwargs,
        )

        dtype = self._get_result_dtype(orig_values.dtype)
        assert isinstance(dtype, BaseMaskedDtype)
        cls = dtype.construct_array_type()

        if self.kind != "aggregate":
            return cls(res_values.astype(dtype.type, copy=False), mask)
        else:
            return cls(res_values.astype(dtype.type, copy=False), result_mask)

    @final
    def _cython_op_ndim_compat(
        self,
        values: np.ndarray,
        *,
        min_count: int,
        ngroups: int,
        comp_ids: np.ndarray,
        mask: np.ndarray | None = None,
        result_mask: np.ndarray | None = None,
        **kwargs,
    ) -> np.ndarray:
        if values.ndim == 1:
            # expand to 2d, dispatch, then squeeze if appropriate
            values2d = values[None, :]
            if mask is not None:
                mask = mask[None, :]
            if result_mask is not None:
                result_mask = result_mask[None, :]
            res = self._call_cython_op(
                values2d,
                min_count=min_count,
                ngroups=ngroups,
                comp_ids=comp_ids,
                mask=mask,
                result_mask=result_mask,
                **kwargs,
            )
            if res.shape[0] == 1:
                return res[0]

            # otherwise we have OHLC
            return res.T

        return self._call_cython_op(
            values,
            min_count=min_count,
            ngroups=ngroups,
            comp_ids=comp_ids,
            mask=mask,
            result_mask=result_mask,
            **kwargs,
        )

    @final
    def _call_cython_op(
        self,
        values: np.ndarray,  # np.ndarray[ndim=2]
        *,
        min_count: int,
        ngroups: int,
        comp_ids: np.ndarray,
        mask: np.ndarray | None,
        result_mask: np.ndarray | None,
        **kwargs,
    ) -> np.ndarray:  # np.ndarray[ndim=2]
        orig_values = values

        dtype = values.dtype
        is_numeric = is_numeric_dtype(dtype)

        is_datetimelike = needs_i8_conversion(dtype)

        if is_datetimelike:
            values = values.view("int64")
            is_numeric = True
        elif is_bool_dtype(dtype):
            values = values.astype("int64")
        elif is_integer_dtype(dtype):
            # GH#43329 If the dtype is explicitly of type uint64 the type is not
            # changed to prevent overflow.
            if dtype != np.uint64:
                values = values.astype(np.int64, copy=False)
        elif is_numeric:
            if not is_complex_dtype(dtype):
                values = ensure_float64(values)

        values = values.T
        if mask is not None:
            mask = mask.T
            if result_mask is not None:
                result_mask = result_mask.T

        out_shape = self._get_output_shape(ngroups, values)
        func, values = self.get_cython_func_and_vals(values, is_numeric)
        out_dtype = self.get_out_dtype(values.dtype)

        result = maybe_fill(np.empty(out_shape, dtype=out_dtype))
        if self.kind == "aggregate":
            counts = np.zeros(ngroups, dtype=np.int64)
            if self.how in ["min", "max", "mean"]:
                func(
                    result,
                    counts,
                    values,
                    comp_ids,
                    min_count,
                    mask=mask,
                    result_mask=result_mask,
                    is_datetimelike=is_datetimelike,
                )
            elif self.how in ["add"]:
                # We support datetimelike
                func(
                    result,
                    counts,
                    values,
                    comp_ids,
                    min_count,
                    datetimelike=is_datetimelike,
                )
            else:
                func(result, counts, values, comp_ids, min_count)
        else:
            # TODO: min_count
            if self.uses_mask():
                func(
                    result,
                    values,
                    comp_ids,
                    ngroups,
                    is_datetimelike,
                    mask=mask,
                    **kwargs,
                )
            else:
                func(result, values, comp_ids, ngroups, is_datetimelike, **kwargs)

        if self.kind == "aggregate":
            # i.e. counts is defined.  Locations where count<min_count
            # need to have the result set to np.nan, which may require casting,
            # see GH#40767
            if is_integer_dtype(result.dtype) and not is_datetimelike:
                cutoff = max(1, min_count)
                empty_groups = counts < cutoff
                if empty_groups.any():
                    # Note: this conversion could be lossy, see GH#40767
                    result = result.astype("float64")
                    result[empty_groups] = np.nan

        result = result.T

        if self.how not in self.cast_blocklist:
            # e.g. if we are int64 and need to restore to datetime64/timedelta64
            # "rank" is the only member of cast_blocklist we get here
            res_dtype = self._get_result_dtype(orig_values.dtype)
            op_result = maybe_downcast_to_dtype(result, res_dtype)
        else:
            op_result = result

        # error: Incompatible return value type (got "Union[ExtensionArray, ndarray]",
        # expected "ndarray")
        return op_result  # type: ignore[return-value]

    @final
    def cython_operation(
        self,
        *,
        values: ArrayLike,
        axis: int,
        min_count: int = -1,
        comp_ids: np.ndarray,
        ngroups: int,
        **kwargs,
    ) -> ArrayLike:
        """
        Call our cython function, with appropriate pre- and post- processing.
        """
        if values.ndim > 2:
            raise NotImplementedError("number of dimensions is currently limited to 2")
        elif values.ndim == 2:
            assert axis == 1, axis
        elif not is_1d_only_ea_obj(values):
            # Note: it is *not* the case that axis is always 0 for 1-dim values,
            #  as we can have 1D ExtensionArrays that we need to treat as 2D
            assert axis == 0

        dtype = values.dtype
        is_numeric = is_numeric_dtype(dtype)

        # can we do this operation with our cython functions
        # if not raise NotImplementedError
        self._disallow_invalid_ops(dtype, is_numeric)

        if not isinstance(values, np.ndarray):
            # i.e. ExtensionArray
            return self._ea_wrap_cython_operation(
                values,
                min_count=min_count,
                ngroups=ngroups,
                comp_ids=comp_ids,
                **kwargs,
            )

        return self._cython_op_ndim_compat(
            values,
            min_count=min_count,
            ngroups=ngroups,
            comp_ids=comp_ids,
            mask=None,
            **kwargs,
        )


class BaseGrouper:
    """
    This is an internal Grouper class, which actually holds
    the generated groups

    Parameters
    ----------
    axis : Index
    groupings : Sequence[Grouping]
        all the grouping instances to handle in this grouper
        for example for grouper list to groupby, need to pass the list
    sort : bool, default True
        whether this grouper will give sorted result or not
    group_keys : bool, default True
    mutated : bool, default False
    indexer : np.ndarray[np.intp], optional
        the indexer created by Grouper
        some groupers (TimeGrouper) will sort its axis and its
        group_info is also sorted, so need the indexer to reorder

    """

    axis: Index

    def __init__(
        self,
        axis: Index,
        groupings: Sequence[grouper.Grouping],
        sort: bool = True,
        group_keys: bool = True,
        mutated: bool = False,
        indexer: npt.NDArray[np.intp] | None = None,
        dropna: bool = True,
    ):
        assert isinstance(axis, Index), axis

        self.axis = axis
        self._groupings: list[grouper.Grouping] = list(groupings)
        self._sort = sort
        self.group_keys = group_keys
        self.mutated = mutated
        self.indexer = indexer
        self.dropna = dropna

    @property
    def groupings(self) -> list[grouper.Grouping]:
        return self._groupings

    @property
    def shape(self) -> Shape:
        return tuple(ping.ngroups for ping in self.groupings)

    def __iter__(self):
        return iter(self.indices)

    @property
    def nkeys(self) -> int:
        return len(self.groupings)

    def get_iterator(
        self, data: NDFrameT, axis: int = 0
    ) -> Iterator[tuple[Hashable, NDFrameT]]:
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        splitter = self._get_splitter(data, axis=axis)
        keys = self.group_keys_seq
        for key, group in zip(keys, splitter):
            yield key, group.__finalize__(data, method="groupby")

    @final
    def _get_splitter(self, data: NDFrame, axis: int = 0) -> DataSplitter:
        """
        Returns
        -------
        Generator yielding subsetted objects

        __finalize__ has not been called for the subsetted objects returned.
        """
        ids, _, ngroups = self.group_info
        return get_splitter(data, ids, ngroups, axis=axis)

    def _get_grouper(self):
        """
        We are a grouper as part of another's groupings.

        We have a specific method of grouping, so cannot
        convert to a Index for our grouper.
        """
        return self.groupings[0].grouping_vector

    @final
    @cache_readonly
    def group_keys_seq(self):
        if len(self.groupings) == 1:
            return self.levels[0]
        else:
            ids, _, ngroups = self.group_info

            # provide "flattened" iterator for multi-group setting
            return get_flattened_list(ids, ngroups, self.levels, self.codes)

    @final
    def apply(
        self, f: Callable, data: DataFrame | Series, axis: int = 0
    ) -> tuple[list, bool]:
        mutated = self.mutated
        splitter = self._get_splitter(data, axis=axis)
        group_keys = self.group_keys_seq
        result_values = []

        # This calls DataSplitter.__iter__
        zipped = zip(group_keys, splitter)

        for key, group in zipped:
            group = group.__finalize__(data, method="groupby")
            object.__setattr__(group, "name", key)

            # group might be modified
            group_axes = group.axes
            res = f(group)
            if not mutated and not _is_indexed_like(res, group_axes, axis):
                mutated = True
            result_values.append(res)

        # getattr pattern for __name__ is needed for functools.partial objects
        if len(group_keys) == 0 and getattr(f, "__name__", None) not in [
            "idxmin",
            "idxmax",
            "nanargmin",
            "nanargmax",
        ]:
            # If group_keys is empty, then no function calls have been made,
            #  so we will not have raised even if this is an invalid dtype.
            #  So do one dummy call here to raise appropriate TypeError.
            f(data.iloc[:0])

        return result_values, mutated

    @cache_readonly
    def indices(self) -> dict[Hashable, npt.NDArray[np.intp]]:
        """dict {group name -> group indices}"""
        if len(self.groupings) == 1 and isinstance(self.result_index, CategoricalIndex):
            # This shows unused categories in indices GH#38642
            return self.groupings[0].indices
        codes_list = [ping.codes for ping in self.groupings]
        keys = [ping.group_index for ping in self.groupings]
        return get_indexer_dict(codes_list, keys)

    @final
    @property
    def codes(self) -> list[np.ndarray]:
        return [ping.codes for ping in self.groupings]

    @property
    def levels(self) -> list[Index]:
        return [ping.group_index for ping in self.groupings]

    @property
    def names(self) -> list[Hashable]:
        return [ping.name for ping in self.groupings]

    @final
    def size(self) -> Series:
        """
        Compute group sizes.
        """
        ids, _, ngroups = self.group_info
        out: np.ndarray | list
        if ngroups:
            out = np.bincount(ids[ids != -1], minlength=ngroups)
        else:
            out = []
        return Series(out, index=self.result_index, dtype="int64")

    @cache_readonly
    def groups(self) -> dict[Hashable, np.ndarray]:
        """dict {group name -> group labels}"""
        if len(self.groupings) == 1:
            return self.groupings[0].groups
        else:
            to_groupby = zip(*(ping.grouping_vector for ping in self.groupings))
            index = Index(to_groupby)
            return self.axis.groupby(index)

    @final
    @cache_readonly
    def is_monotonic(self) -> bool:
        # return if my group orderings are monotonic
        return Index(self.group_info[0]).is_monotonic

    @cache_readonly
    def group_info(self) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.intp], int]:
        comp_ids, obs_group_ids = self._get_compressed_codes()

        ngroups = len(obs_group_ids)
        comp_ids = ensure_platform_int(comp_ids)

        return comp_ids, obs_group_ids, ngroups

    @final
    @cache_readonly
    def codes_info(self) -> npt.NDArray[np.intp]:
        # return the codes of items in original grouped axis
        ids, _, _ = self.group_info
        if self.indexer is not None:
            sorter = np.lexsort((ids, self.indexer))
            ids = ids[sorter]
            ids = ensure_platform_int(ids)
            # TODO: if numpy annotates np.lexsort, this ensure_platform_int
            #  may become unnecessary
        return ids

    @final
    def _get_compressed_codes(self) -> tuple[np.ndarray, npt.NDArray[np.intp]]:
        # The first returned ndarray may have any signed integer dtype
        if len(self.groupings) > 1:
            group_index = get_group_index(self.codes, self.shape, sort=True, xnull=True)
            return compress_group_index(group_index, sort=self._sort)

        ping = self.groupings[0]
        return ping.codes, np.arange(len(ping.group_index), dtype=np.intp)

    @final
    @cache_readonly
    def ngroups(self) -> int:
        return len(self.result_index)

    @property
    def reconstructed_codes(self) -> list[np.ndarray]:
        codes = self.codes
        ids, obs_ids, _ = self.group_info
        return decons_obs_group_ids(ids, obs_ids, self.shape, codes, xnull=True)

    @final
    @cache_readonly
    def result_arraylike(self) -> ArrayLike:
        """
        Analogous to result_index, but returning an ndarray/ExtensionArray
        allowing us to retain ExtensionDtypes not supported by Index.
        """
        # TODO(ExtensionIndex): once Index supports arbitrary EAs, this can
        #  be removed in favor of result_index
        if len(self.groupings) == 1:
            return self.groupings[0].group_arraylike

        # result_index is MultiIndex
        return self.result_index._values

    @cache_readonly
    def result_index(self) -> Index:
        if len(self.groupings) == 1:
            return self.groupings[0].result_index.rename(self.names[0])

        codes = self.reconstructed_codes
        levels = [ping.result_index for ping in self.groupings]
        return MultiIndex(
            levels=levels, codes=codes, verify_integrity=False, names=self.names
        )

    @final
    def get_group_levels(self) -> list[ArrayLike]:
        # Note: only called from _insert_inaxis_grouper_inplace, which
        #  is only called for BaseGrouper, never for BinGrouper
        if len(self.groupings) == 1:
            return [self.groupings[0].group_arraylike]

        name_list = []
        for ping, codes in zip(self.groupings, self.reconstructed_codes):
            codes = ensure_platform_int(codes)
            levels = ping.group_arraylike.take(codes)

            name_list.append(levels)

        return name_list

    # ------------------------------------------------------------
    # Aggregation functions

    @final
    def _cython_operation(
        self,
        kind: str,
        values,
        how: str,
        axis: int,
        min_count: int = -1,
        **kwargs,
    ) -> ArrayLike:
        """
        Returns the values of a cython operation.
        """
        assert kind in ["transform", "aggregate"]

        cy_op = WrappedCythonOp(kind=kind, how=how)

        ids, _, _ = self.group_info
        ngroups = self.ngroups
        return cy_op.cython_operation(
            values=values,
            axis=axis,
            min_count=min_count,
            comp_ids=ids,
            ngroups=ngroups,
            **kwargs,
        )

    @final
    def agg_series(
        self, obj: Series, func: Callable, preserve_dtype: bool = False
    ) -> ArrayLike:
        """
        Parameters
        ----------
        obj : Series
        func : function taking a Series and returning a scalar-like
        preserve_dtype : bool
            Whether the aggregation is known to be dtype-preserving.

        Returns
        -------
        np.ndarray or ExtensionArray
        """
        # test_groupby_empty_with_category gets here with self.ngroups == 0
        #  and len(obj) > 0

        if len(obj) == 0:
            # SeriesGrouper would raise if we were to call _aggregate_series_fast
            result = self._aggregate_series_pure_python(obj, func)

        elif not isinstance(obj._values, np.ndarray):
            result = self._aggregate_series_pure_python(obj, func)

            # we can preserve a little bit more aggressively with EA dtype
            #  because maybe_cast_pointwise_result will do a try/except
            #  with _from_sequence.  NB we are assuming here that _from_sequence
            #  is sufficiently strict that it casts appropriately.
            preserve_dtype = True

        else:
            result = self._aggregate_series_pure_python(obj, func)

        npvalues = lib.maybe_convert_objects(result, try_float=False)
        if preserve_dtype:
            out = maybe_cast_pointwise_result(npvalues, obj.dtype, numeric_only=True)
        else:
            out = npvalues
        return out

    @final
    def _aggregate_series_pure_python(
        self, obj: Series, func: Callable
    ) -> npt.NDArray[np.object_]:
        ids, _, ngroups = self.group_info

        counts = np.zeros(ngroups, dtype=int)
        result = np.empty(ngroups, dtype="O")
        initialized = False

        # equiv: splitter = self._get_splitter(obj, axis=0)
        splitter = get_splitter(obj, ids, ngroups, axis=0)

        for i, group in enumerate(splitter):
            group = group.__finalize__(obj, method="groupby")
            res = func(group)
            res = libreduction.extract_result(res)

            if not initialized:
                # We only do this validation on the first iteration
                libreduction.check_result_array(res, group.dtype)
                initialized = True

            counts[i] = group.shape[0]
            result[i] = res

        return result


class BinGrouper(BaseGrouper):
    """
    This is an internal Grouper class

    Parameters
    ----------
    bins : the split index of binlabels to group the item of axis
    binlabels : the label list
    mutated : bool, default False
    indexer : np.ndarray[np.intp]

    Examples
    --------
    bins: [2, 4, 6, 8, 10]
    binlabels: DatetimeIndex(['2005-01-01', '2005-01-03',
        '2005-01-05', '2005-01-07', '2005-01-09'],
        dtype='datetime64[ns]', freq='2D')

    the group_info, which contains the label of each item in grouped
    axis, the index of label in label list, group number, is

    (array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4]), array([0, 1, 2, 3, 4]), 5)

    means that, the grouped axis has 10 items, can be grouped into 5
    labels, the first and second items belong to the first label, the
    third and forth items belong to the second label, and so on

    """

    bins: npt.NDArray[np.int64]
    binlabels: Index
    mutated: bool

    def __init__(
        self,
        bins,
        binlabels,
        mutated: bool = False,
        indexer=None,
    ):
        self.bins = ensure_int64(bins)
        self.binlabels = ensure_index(binlabels)
        self.mutated = mutated
        self.indexer = indexer

        # These lengths must match, otherwise we could call agg_series
        #  with empty self.bins, which would raise in libreduction.
        assert len(self.binlabels) == len(self.bins)

    @cache_readonly
    def groups(self):
        """dict {group name -> group labels}"""
        # this is mainly for compat
        # GH 3881
        result = {
            key: value
            for key, value in zip(self.binlabels, self.bins)
            if key is not NaT
        }
        return result

    @property
    def nkeys(self) -> int:
        # still matches len(self.groupings), but we can hard-code
        return 1

    def _get_grouper(self):
        """
        We are a grouper as part of another's groupings.

        We have a specific method of grouping, so cannot
        convert to a Index for our grouper.
        """
        return self

    def get_iterator(self, data: NDFrame, axis: int = 0):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        if axis == 0:
            slicer = lambda start, edge: data.iloc[start:edge]
        else:
            slicer = lambda start, edge: data.iloc[:, start:edge]

        length = len(data.axes[axis])

        start = 0
        for edge, label in zip(self.bins, self.binlabels):
            if label is not NaT:
                yield label, slicer(start, edge)
            start = edge

        if start < length:
            yield self.binlabels[-1], slicer(start, None)

    @cache_readonly
    def indices(self):
        indices = collections.defaultdict(list)

        i = 0
        for label, bin in zip(self.binlabels, self.bins):
            if i < bin:
                if label is not NaT:
                    indices[label] = list(range(i, bin))
                i = bin
        return indices

    @cache_readonly
    def group_info(self) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.intp], int]:
        ngroups = self.ngroups
        obs_group_ids = np.arange(ngroups, dtype=np.intp)
        rep = np.diff(np.r_[0, self.bins])

        rep = ensure_platform_int(rep)
        if ngroups == len(self.bins):
            comp_ids = np.repeat(np.arange(ngroups), rep)
        else:
            comp_ids = np.repeat(np.r_[-1, np.arange(ngroups)], rep)

        return (
            ensure_platform_int(comp_ids),
            obs_group_ids,
            ngroups,
        )

    @cache_readonly
    def reconstructed_codes(self) -> list[np.ndarray]:
        # get unique result indices, and prepend 0 as groupby starts from the first
        return [np.r_[0, np.flatnonzero(self.bins[1:] != self.bins[:-1]) + 1]]

    @cache_readonly
    def result_index(self):
        if len(self.binlabels) != 0 and isna(self.binlabels[0]):
            return self.binlabels[1:]

        return self.binlabels

    @property
    def levels(self) -> list[Index]:
        return [self.binlabels]

    @property
    def names(self) -> list[Hashable]:
        return [self.binlabels.name]

    @property
    def groupings(self) -> list[grouper.Grouping]:
        lev = self.binlabels
        ping = grouper.Grouping(lev, lev, in_axis=False, level=None)
        return [ping]

    def _aggregate_series_fast(self, obj: Series, func: Callable) -> np.ndarray:
        # -> np.ndarray[object]
        raise NotImplementedError(
            "This should not be reached; use _aggregate_series_pure_python"
        )


def _is_indexed_like(obj, axes, axis: int) -> bool:
    if isinstance(obj, Series):
        if len(axes) > 1:
            return False
        return obj.axes[axis].equals(axes[axis])
    elif isinstance(obj, DataFrame):
        return obj.axes[axis].equals(axes[axis])

    return False


# ----------------------------------------------------------------------
# Splitting / application


class DataSplitter(Generic[NDFrameT]):
    def __init__(
        self,
        data: NDFrameT,
        labels: npt.NDArray[np.intp],
        ngroups: int,
        axis: int = 0,
    ):
        self.data = data
        self.labels = ensure_platform_int(labels)  # _should_ already be np.intp
        self.ngroups = ngroups

        self.axis = axis
        assert isinstance(axis, int), axis

    @cache_readonly
    def slabels(self) -> npt.NDArray[np.intp]:
        # Sorted labels
        return self.labels.take(self._sort_idx)

    @cache_readonly
    def _sort_idx(self) -> npt.NDArray[np.intp]:
        # Counting sort indexer
        return get_group_index_sorter(self.labels, self.ngroups)

    def __iter__(self):
        sdata = self.sorted_data

        if self.ngroups == 0:
            # we are inside a generator, rather than raise StopIteration
            # we merely return signal the end
            return

        starts, ends = lib.generate_slices(self.slabels, self.ngroups)

        for start, end in zip(starts, ends):
            yield self._chop(sdata, slice(start, end))

    @cache_readonly
    def sorted_data(self) -> NDFrameT:
        return self.data.take(self._sort_idx, axis=self.axis)

    def _chop(self, sdata, slice_obj: slice) -> NDFrame:
        raise AbstractMethodError(self)


class SeriesSplitter(DataSplitter):
    def _chop(self, sdata: Series, slice_obj: slice) -> Series:
        # fastpath equivalent to `sdata.iloc[slice_obj]`
        mgr = sdata._mgr.get_slice(slice_obj)
        # __finalize__ not called here, must be applied by caller if applicable
        return sdata._constructor(mgr, name=sdata.name, fastpath=True)


class FrameSplitter(DataSplitter):
    def _chop(self, sdata: DataFrame, slice_obj: slice) -> DataFrame:
        # Fastpath equivalent to:
        # if self.axis == 0:
        #     return sdata.iloc[slice_obj]
        # else:
        #     return sdata.iloc[:, slice_obj]
        mgr = sdata._mgr.get_slice(slice_obj, axis=1 - self.axis)
        # __finalize__ not called here, must be applied by caller if applicable
        return sdata._constructor(mgr)


def get_splitter(
    data: NDFrame, labels: np.ndarray, ngroups: int, axis: int = 0
) -> DataSplitter:
    if isinstance(data, Series):
        klass: type[DataSplitter] = SeriesSplitter
    else:
        # i.e. DataFrame
        klass = FrameSplitter

    return klass(data, labels, ngroups, axis)
