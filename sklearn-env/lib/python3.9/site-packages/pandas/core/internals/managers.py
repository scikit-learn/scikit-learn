from __future__ import annotations

import itertools
from typing import (
    Any,
    Callable,
    Hashable,
    Sequence,
    TypeVar,
    cast,
)
import warnings

import numpy as np

from pandas._libs import (
    internals as libinternals,
    lib,
)
from pandas._libs.internals import BlockPlacement
from pandas._typing import (
    ArrayLike,
    DtypeObj,
    Shape,
    npt,
    type_t,
)
from pandas.errors import PerformanceWarning
from pandas.util._decorators import cache_readonly
from pandas.util._exceptions import find_stack_level
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import infer_dtype_from_scalar
from pandas.core.dtypes.common import (
    ensure_platform_int,
    is_1d_only_ea_dtype,
    is_dtype_equal,
    is_list_like,
    needs_i8_conversion,
)
from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCSeries,
)
from pandas.core.dtypes.missing import (
    array_equals,
    isna,
)

import pandas.core.algorithms as algos
from pandas.core.arrays._mixins import NDArrayBackedExtensionArray
from pandas.core.arrays.sparse import SparseDtype
from pandas.core.construction import (
    ensure_wrapped_if_datetimelike,
    extract_array,
)
from pandas.core.indexers import maybe_convert_indices
from pandas.core.indexes.api import (
    Float64Index,
    Index,
    ensure_index,
)
from pandas.core.internals.base import (
    DataManager,
    SingleDataManager,
    interleaved_dtype,
)
from pandas.core.internals.blocks import (
    Block,
    DatetimeTZBlock,
    NumpyBlock,
    ensure_block_shape,
    extend_blocks,
    get_block_type,
    new_block,
    new_block_2d,
)
from pandas.core.internals.ops import (
    blockwise_all,
    operate_blockwise,
)

T = TypeVar("T", bound="BaseBlockManager")


class BaseBlockManager(DataManager):
    """
    Core internal data structure to implement DataFrame, Series, etc.

    Manage a bunch of labeled 2D mixed-type ndarrays. Essentially it's a
    lightweight blocked set of labeled data to be manipulated by the DataFrame
    public API class

    Attributes
    ----------
    shape
    ndim
    axes
    values
    items

    Methods
    -------
    set_axis(axis, new_labels)
    copy(deep=True)

    get_dtypes

    apply(func, axes, block_filter_fn)

    get_bool_data
    get_numeric_data

    get_slice(slice_like, axis)
    get(label)
    iget(loc)

    take(indexer, axis)
    reindex_axis(new_labels, axis)
    reindex_indexer(new_labels, indexer, axis)

    delete(label)
    insert(loc, label, value)
    set(label, value)

    Parameters
    ----------
    blocks: Sequence of Block
    axes: Sequence of Index
    verify_integrity: bool, default True

    Notes
    -----
    This is *not* a public API class
    """

    __slots__ = ()

    _blknos: npt.NDArray[np.intp]
    _blklocs: npt.NDArray[np.intp]
    blocks: tuple[Block, ...]
    axes: list[Index]

    ndim: int
    _known_consolidated: bool
    _is_consolidated: bool

    def __init__(self, blocks, axes, verify_integrity: bool = True):
        raise NotImplementedError

    @classmethod
    def from_blocks(cls: type_t[T], blocks: list[Block], axes: list[Index]) -> T:
        raise NotImplementedError

    @property
    def blknos(self) -> npt.NDArray[np.intp]:
        """
        Suppose we want to find the array corresponding to our i'th column.

        blknos[i] identifies the block from self.blocks that contains this column.

        blklocs[i] identifies the column of interest within
        self.blocks[self.blknos[i]]
        """
        if self._blknos is None:
            # Note: these can be altered by other BlockManager methods.
            self._rebuild_blknos_and_blklocs()

        return self._blknos

    @property
    def blklocs(self) -> npt.NDArray[np.intp]:
        """
        See blknos.__doc__
        """
        if self._blklocs is None:
            # Note: these can be altered by other BlockManager methods.
            self._rebuild_blknos_and_blklocs()

        return self._blklocs

    def make_empty(self: T, axes=None) -> T:
        """return an empty BlockManager with the items axis of len 0"""
        if axes is None:
            axes = [Index([])] + self.axes[1:]

        # preserve dtype if possible
        if self.ndim == 1:
            assert isinstance(self, SingleBlockManager)  # for mypy
            blk = self.blocks[0]
            arr = blk.values[:0]
            bp = BlockPlacement(slice(0, 0))
            nb = blk.make_block_same_class(arr, placement=bp)
            blocks = [nb]
        else:
            blocks = []
        return type(self).from_blocks(blocks, axes)

    def __nonzero__(self) -> bool:
        return True

    # Python3 compat
    __bool__ = __nonzero__

    def _normalize_axis(self, axis: int) -> int:
        # switch axis to follow BlockManager logic
        if self.ndim == 2:
            axis = 1 if axis == 0 else 0
        return axis

    def set_axis(self, axis: int, new_labels: Index) -> None:
        # Caller is responsible for ensuring we have an Index object.
        self._validate_set_axis(axis, new_labels)
        self.axes[axis] = new_labels

    @property
    def is_single_block(self) -> bool:
        # Assumes we are 2D; overridden by SingleBlockManager
        return len(self.blocks) == 1

    @property
    def items(self) -> Index:
        return self.axes[0]

    def get_dtypes(self):
        dtypes = np.array([blk.dtype for blk in self.blocks])
        return dtypes.take(self.blknos)

    @property
    def arrays(self) -> list[ArrayLike]:
        """
        Quick access to the backing arrays of the Blocks.

        Only for compatibility with ArrayManager for testing convenience.
        Not to be used in actual code, and return value is not the same as the
        ArrayManager method (list of 1D arrays vs iterator of 2D ndarrays / 1D EAs).
        """
        return [blk.values for blk in self.blocks]

    def __repr__(self) -> str:
        output = type(self).__name__
        for i, ax in enumerate(self.axes):
            if i == 0:
                output += f"\nItems: {ax}"
            else:
                output += f"\nAxis {i}: {ax}"

        for block in self.blocks:
            output += f"\n{block}"
        return output

    def apply(
        self: T,
        f,
        align_keys: list[str] | None = None,
        ignore_failures: bool = False,
        **kwargs,
    ) -> T:
        """
        Iterate over the blocks, collect and create a new BlockManager.

        Parameters
        ----------
        f : str or callable
            Name of the Block method to apply.
        align_keys: List[str] or None, default None
        ignore_failures: bool, default False
        **kwargs
            Keywords to pass to `f`

        Returns
        -------
        BlockManager
        """
        assert "filter" not in kwargs

        align_keys = align_keys or []
        result_blocks: list[Block] = []
        # fillna: Series/DataFrame is responsible for making sure value is aligned

        aligned_args = {k: kwargs[k] for k in align_keys}

        for b in self.blocks:

            if aligned_args:

                for k, obj in aligned_args.items():
                    if isinstance(obj, (ABCSeries, ABCDataFrame)):
                        # The caller is responsible for ensuring that
                        #  obj.axes[-1].equals(self.items)
                        if obj.ndim == 1:
                            kwargs[k] = obj.iloc[b.mgr_locs.indexer]._values
                        else:
                            kwargs[k] = obj.iloc[:, b.mgr_locs.indexer]._values
                    else:
                        # otherwise we have an ndarray
                        kwargs[k] = obj[b.mgr_locs.indexer]

            try:
                if callable(f):
                    applied = b.apply(f, **kwargs)
                else:
                    applied = getattr(b, f)(**kwargs)
            except (TypeError, NotImplementedError):
                if not ignore_failures:
                    raise
                continue
            result_blocks = extend_blocks(applied, result_blocks)

        if ignore_failures:
            return self._combine(result_blocks)

        out = type(self).from_blocks(result_blocks, self.axes)
        return out

    def where(self: T, other, cond, align: bool) -> T:
        if align:
            align_keys = ["other", "cond"]
        else:
            align_keys = ["cond"]
            other = extract_array(other, extract_numpy=True)

        return self.apply(
            "where",
            align_keys=align_keys,
            other=other,
            cond=cond,
        )

    def setitem(self: T, indexer, value) -> T:
        """
        Set values with indexer.

        For SingleBlockManager, this backs s[indexer] = value
        """
        return self.apply("setitem", indexer=indexer, value=value)

    def putmask(self, mask, new, align: bool = True):

        if align:
            align_keys = ["new", "mask"]
        else:
            align_keys = ["mask"]
            new = extract_array(new, extract_numpy=True)

        return self.apply(
            "putmask",
            align_keys=align_keys,
            mask=mask,
            new=new,
        )

    def diff(self: T, n: int, axis: int) -> T:
        axis = self._normalize_axis(axis)
        return self.apply("diff", n=n, axis=axis)

    def interpolate(self: T, **kwargs) -> T:
        return self.apply("interpolate", **kwargs)

    def shift(self: T, periods: int, axis: int, fill_value) -> T:
        axis = self._normalize_axis(axis)
        if fill_value is lib.no_default:
            fill_value = None

        if (
            axis == 0
            and self.ndim == 2
            and (
                self.nblocks > 1
                or (
                    # If we only have one block and we know that we can't
                    #  keep the same dtype (i.e. the _can_hold_element check)
                    #  then we can go through the reindex_indexer path
                    #  (and avoid casting logic in the Block method).
                    #  The exception to this (until 2.0) is datetimelike
                    #  dtypes with integers, which cast.
                    not self.blocks[0]._can_hold_element(fill_value)
                    # TODO(2.0): remove special case for integer-with-datetimelike
                    #  once deprecation is enforced
                    and not (
                        lib.is_integer(fill_value)
                        and needs_i8_conversion(self.blocks[0].dtype)
                    )
                )
            )
        ):
            # GH#35488 we need to watch out for multi-block cases
            # We only get here with fill_value not-lib.no_default
            ncols = self.shape[0]
            nper = abs(periods)
            nper = min(nper, ncols)
            if periods > 0:
                indexer = np.array(
                    [-1] * nper + list(range(ncols - periods)), dtype=np.intp
                )
            else:
                indexer = np.array(
                    list(range(nper, ncols)) + [-1] * nper, dtype=np.intp
                )
            result = self.reindex_indexer(
                self.items,
                indexer,
                axis=0,
                fill_value=fill_value,
                allow_dups=True,
                consolidate=False,
            )
            return result

        return self.apply("shift", periods=periods, axis=axis, fill_value=fill_value)

    def fillna(self: T, value, limit, inplace: bool, downcast) -> T:
        return self.apply(
            "fillna", value=value, limit=limit, inplace=inplace, downcast=downcast
        )

    def astype(self: T, dtype, copy: bool = False, errors: str = "raise") -> T:
        return self.apply("astype", dtype=dtype, copy=copy, errors=errors)

    def convert(
        self: T,
        copy: bool = True,
        datetime: bool = True,
        numeric: bool = True,
        timedelta: bool = True,
    ) -> T:
        return self.apply(
            "convert",
            copy=copy,
            datetime=datetime,
            numeric=numeric,
            timedelta=timedelta,
        )

    def replace(self: T, to_replace, value, inplace: bool) -> T:
        inplace = validate_bool_kwarg(inplace, "inplace")
        # NDFrame.replace ensures the not-is_list_likes here
        assert not is_list_like(to_replace)
        assert not is_list_like(value)
        return self.apply(
            "replace", to_replace=to_replace, value=value, inplace=inplace
        )

    def replace_regex(self, **kwargs):
        return self.apply("_replace_regex", **kwargs)

    def replace_list(
        self: T,
        src_list: list[Any],
        dest_list: list[Any],
        inplace: bool = False,
        regex: bool = False,
    ) -> T:
        """do a list replace"""
        inplace = validate_bool_kwarg(inplace, "inplace")

        bm = self.apply(
            "replace_list",
            src_list=src_list,
            dest_list=dest_list,
            inplace=inplace,
            regex=regex,
        )
        bm._consolidate_inplace()
        return bm

    def to_native_types(self: T, **kwargs) -> T:
        """
        Convert values to native types (strings / python objects) that are used
        in formatting (repr / csv).
        """
        return self.apply("to_native_types", **kwargs)

    @property
    def is_numeric_mixed_type(self) -> bool:
        return all(block.is_numeric for block in self.blocks)

    @property
    def any_extension_types(self) -> bool:
        """Whether any of the blocks in this manager are extension blocks"""
        return any(block.is_extension for block in self.blocks)

    @property
    def is_view(self) -> bool:
        """return a boolean if we are a single block and are a view"""
        if len(self.blocks) == 1:
            return self.blocks[0].is_view

        # It is technically possible to figure out which blocks are views
        # e.g. [ b.values.base is not None for b in self.blocks ]
        # but then we have the case of possibly some blocks being a view
        # and some blocks not. setting in theory is possible on the non-view
        # blocks w/o causing a SettingWithCopy raise/warn. But this is a bit
        # complicated

        return False

    def _get_data_subset(self: T, predicate: Callable) -> T:
        blocks = [blk for blk in self.blocks if predicate(blk.values)]
        return self._combine(blocks, copy=False)

    def get_bool_data(self: T, copy: bool = False) -> T:
        """
        Select blocks that are bool-dtype and columns from object-dtype blocks
        that are all-bool.

        Parameters
        ----------
        copy : bool, default False
            Whether to copy the blocks
        """

        new_blocks = []

        for blk in self.blocks:
            if blk.dtype == bool:
                new_blocks.append(blk)

            elif blk.is_object:
                nbs = blk._split()
                for nb in nbs:
                    if nb.is_bool:
                        new_blocks.append(nb)

        return self._combine(new_blocks, copy)

    def get_numeric_data(self: T, copy: bool = False) -> T:
        """
        Parameters
        ----------
        copy : bool, default False
            Whether to copy the blocks
        """
        numeric_blocks = [blk for blk in self.blocks if blk.is_numeric]
        if len(numeric_blocks) == len(self.blocks):
            # Avoid somewhat expensive _combine
            if copy:
                return self.copy(deep=True)
            return self
        return self._combine(numeric_blocks, copy)

    def _combine(
        self: T, blocks: list[Block], copy: bool = True, index: Index | None = None
    ) -> T:
        """return a new manager with the blocks"""
        if len(blocks) == 0:
            if self.ndim == 2:
                # retain our own Index dtype
                if index is not None:
                    axes = [self.items[:0], index]
                else:
                    axes = [self.items[:0]] + self.axes[1:]
                return self.make_empty(axes)
            return self.make_empty()

        # FIXME: optimization potential
        indexer = np.sort(np.concatenate([b.mgr_locs.as_array for b in blocks]))
        inv_indexer = lib.get_reverse_indexer(indexer, self.shape[0])

        new_blocks: list[Block] = []
        for b in blocks:
            b = b.copy(deep=copy)
            b.mgr_locs = BlockPlacement(inv_indexer[b.mgr_locs.indexer])
            new_blocks.append(b)

        axes = list(self.axes)
        if index is not None:
            axes[-1] = index
        axes[0] = self.items.take(indexer)

        return type(self).from_blocks(new_blocks, axes)

    @property
    def nblocks(self) -> int:
        return len(self.blocks)

    def copy(self: T, deep=True) -> T:
        """
        Make deep or shallow copy of BlockManager

        Parameters
        ----------
        deep : bool or string, default True
            If False, return shallow copy (do not copy data)
            If 'all', copy data and a deep copy of the index

        Returns
        -------
        BlockManager
        """
        # this preserves the notion of view copying of axes
        if deep:
            # hit in e.g. tests.io.json.test_pandas

            def copy_func(ax):
                return ax.copy(deep=True) if deep == "all" else ax.view()

            new_axes = [copy_func(ax) for ax in self.axes]
        else:
            new_axes = list(self.axes)

        res = self.apply("copy", deep=deep)

        res.axes = new_axes

        if self.ndim > 1:
            # Avoid needing to re-compute these
            blknos = self._blknos
            if blknos is not None:
                res._blknos = blknos.copy()
                res._blklocs = self._blklocs.copy()

        if deep:
            res._consolidate_inplace()
        return res

    def consolidate(self: T) -> T:
        """
        Join together blocks having same dtype

        Returns
        -------
        y : BlockManager
        """
        if self.is_consolidated():
            return self

        bm = type(self)(self.blocks, self.axes, verify_integrity=False)
        bm._is_consolidated = False
        bm._consolidate_inplace()
        return bm

    def reindex_indexer(
        self: T,
        new_axis: Index,
        indexer,
        axis: int,
        fill_value=None,
        allow_dups: bool = False,
        copy: bool = True,
        consolidate: bool = True,
        only_slice: bool = False,
        *,
        use_na_proxy: bool = False,
    ) -> T:
        """
        Parameters
        ----------
        new_axis : Index
        indexer : ndarray of int64 or None
        axis : int
        fill_value : object, default None
        allow_dups : bool, default False
        copy : bool, default True
        consolidate: bool, default True
            Whether to consolidate inplace before reindexing.
        only_slice : bool, default False
            Whether to take views, not copies, along columns.
        use_na_proxy : bool, default False
            Whether to use a np.void ndarray for newly introduced columns.

        pandas-indexer with -1's only.
        """
        if indexer is None:
            if new_axis is self.axes[axis] and not copy:
                return self

            result = self.copy(deep=copy)
            result.axes = list(self.axes)
            result.axes[axis] = new_axis
            return result

        if consolidate:
            self._consolidate_inplace()

        # some axes don't allow reindexing with dups
        if not allow_dups:
            self.axes[axis]._validate_can_reindex(indexer)

        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

        if axis == 0:
            new_blocks = self._slice_take_blocks_ax0(
                indexer,
                fill_value=fill_value,
                only_slice=only_slice,
                use_na_proxy=use_na_proxy,
            )
        else:
            new_blocks = [
                blk.take_nd(
                    indexer,
                    axis=1,
                    fill_value=(
                        fill_value if fill_value is not None else blk.fill_value
                    ),
                )
                for blk in self.blocks
            ]

        new_axes = list(self.axes)
        new_axes[axis] = new_axis

        new_mgr = type(self).from_blocks(new_blocks, new_axes)
        if axis == 1:
            # We can avoid the need to rebuild these
            new_mgr._blknos = self.blknos.copy()
            new_mgr._blklocs = self.blklocs.copy()
        return new_mgr

    def _slice_take_blocks_ax0(
        self,
        slice_or_indexer: slice | np.ndarray,
        fill_value=lib.no_default,
        only_slice: bool = False,
        *,
        use_na_proxy: bool = False,
    ) -> list[Block]:
        """
        Slice/take blocks along axis=0.

        Overloaded for SingleBlock

        Parameters
        ----------
        slice_or_indexer : slice or np.ndarray[int64]
        fill_value : scalar, default lib.no_default
        only_slice : bool, default False
            If True, we always return views on existing arrays, never copies.
            This is used when called from ops.blockwise.operate_blockwise.
        use_na_proxy : bool, default False
            Whether to use a np.void ndarray for newly introduced columns.

        Returns
        -------
        new_blocks : list of Block
        """
        allow_fill = fill_value is not lib.no_default

        sl_type, slobj, sllen = _preprocess_slice_or_indexer(
            slice_or_indexer, self.shape[0], allow_fill=allow_fill
        )

        if self.is_single_block:
            blk = self.blocks[0]

            if sl_type == "slice":
                # GH#32959 EABlock would fail since we can't make 0-width
                # TODO(EA2D): special casing unnecessary with 2D EAs
                if sllen == 0:
                    return []
                bp = BlockPlacement(slice(0, sllen))
                return [blk.getitem_block_columns(slobj, new_mgr_locs=bp)]
            elif not allow_fill or self.ndim == 1:
                if allow_fill and fill_value is None:
                    fill_value = blk.fill_value

                if not allow_fill and only_slice:
                    # GH#33597 slice instead of take, so we get
                    #  views instead of copies
                    blocks = [
                        blk.getitem_block_columns(
                            slice(ml, ml + 1), new_mgr_locs=BlockPlacement(i)
                        )
                        for i, ml in enumerate(slobj)
                    ]
                    # We have
                    #  all(np.shares_memory(nb.values, blk.values) for nb in blocks)
                    return blocks
                else:
                    bp = BlockPlacement(slice(0, sllen))
                    return [
                        blk.take_nd(
                            slobj,
                            axis=0,
                            new_mgr_locs=bp,
                            fill_value=fill_value,
                        )
                    ]

        if sl_type == "slice":
            blknos = self.blknos[slobj]
            blklocs = self.blklocs[slobj]
        else:
            blknos = algos.take_nd(
                self.blknos, slobj, fill_value=-1, allow_fill=allow_fill
            )
            blklocs = algos.take_nd(
                self.blklocs, slobj, fill_value=-1, allow_fill=allow_fill
            )

        # When filling blknos, make sure blknos is updated before appending to
        # blocks list, that way new blkno is exactly len(blocks).
        blocks = []
        group = not only_slice
        for blkno, mgr_locs in libinternals.get_blkno_placements(blknos, group=group):
            if blkno == -1:
                # If we've got here, fill_value was not lib.no_default

                blocks.append(
                    self._make_na_block(
                        placement=mgr_locs,
                        fill_value=fill_value,
                        use_na_proxy=use_na_proxy,
                    )
                )
            else:
                blk = self.blocks[blkno]

                # Otherwise, slicing along items axis is necessary.
                if not blk._can_consolidate and not blk._validate_ndim:
                    # i.e. we dont go through here for DatetimeTZBlock
                    # A non-consolidatable block, it's easy, because there's
                    # only one item and each mgr loc is a copy of that single
                    # item.
                    for mgr_loc in mgr_locs:
                        newblk = blk.copy(deep=False)
                        newblk.mgr_locs = BlockPlacement(slice(mgr_loc, mgr_loc + 1))
                        blocks.append(newblk)

                else:
                    # GH#32779 to avoid the performance penalty of copying,
                    #  we may try to only slice
                    taker = blklocs[mgr_locs.indexer]
                    max_len = max(len(mgr_locs), taker.max() + 1)
                    if only_slice:
                        taker = lib.maybe_indices_to_slice(taker, max_len)

                    if isinstance(taker, slice):
                        nb = blk.getitem_block_columns(taker, new_mgr_locs=mgr_locs)
                        blocks.append(nb)
                    elif only_slice:
                        # GH#33597 slice instead of take, so we get
                        #  views instead of copies
                        for i, ml in zip(taker, mgr_locs):
                            slc = slice(i, i + 1)
                            bp = BlockPlacement(ml)
                            nb = blk.getitem_block_columns(slc, new_mgr_locs=bp)
                            # We have np.shares_memory(nb.values, blk.values)
                            blocks.append(nb)
                    else:
                        nb = blk.take_nd(taker, axis=0, new_mgr_locs=mgr_locs)
                        blocks.append(nb)

        return blocks

    def _make_na_block(
        self, placement: BlockPlacement, fill_value=None, use_na_proxy: bool = False
    ) -> Block:
        # Note: we only get here with self.ndim == 2

        if use_na_proxy:
            assert fill_value is None
            shape = (len(placement), self.shape[1])
            vals = np.empty(shape, dtype=np.void)
            nb = NumpyBlock(vals, placement, ndim=2)
            return nb

        if fill_value is None:
            fill_value = np.nan
        block_shape = list(self.shape)
        block_shape[0] = len(placement)

        dtype, fill_value = infer_dtype_from_scalar(fill_value)
        # error: Argument "dtype" to "empty" has incompatible type "Union[dtype,
        # ExtensionDtype]"; expected "Union[dtype, None, type, _SupportsDtype, str,
        # Tuple[Any, int], Tuple[Any, Union[int, Sequence[int]]], List[Any], _DtypeDict,
        # Tuple[Any, Any]]"
        block_values = np.empty(block_shape, dtype=dtype)  # type: ignore[arg-type]
        block_values.fill(fill_value)
        return new_block_2d(block_values, placement=placement)

    def take(self: T, indexer, axis: int = 1, verify: bool = True) -> T:
        """
        Take items along any axis.

        indexer : np.ndarray or slice
        axis : int, default 1
        verify : bool, default True
            Check that all entries are between 0 and len(self) - 1, inclusive.
            Pass verify=False if this check has been done by the caller.

        Returns
        -------
        BlockManager
        """
        # We have 6 tests that get here with a slice
        indexer = (
            np.arange(indexer.start, indexer.stop, indexer.step, dtype=np.intp)
            if isinstance(indexer, slice)
            else np.asanyarray(indexer, dtype=np.intp)
        )

        n = self.shape[axis]
        indexer = maybe_convert_indices(indexer, n, verify=verify)

        new_labels = self.axes[axis].take(indexer)
        return self.reindex_indexer(
            new_axis=new_labels,
            indexer=indexer,
            axis=axis,
            allow_dups=True,
            consolidate=False,
        )


class BlockManager(libinternals.BlockManager, BaseBlockManager):
    """
    BaseBlockManager that holds 2D blocks.
    """

    ndim = 2

    # ----------------------------------------------------------------
    # Constructors

    def __init__(
        self,
        blocks: Sequence[Block],
        axes: Sequence[Index],
        verify_integrity: bool = True,
    ):

        if verify_integrity:
            # Assertion disabled for performance
            # assert all(isinstance(x, Index) for x in axes)

            for block in blocks:
                if self.ndim != block.ndim:
                    raise AssertionError(
                        f"Number of Block dimensions ({block.ndim}) must equal "
                        f"number of axes ({self.ndim})"
                    )
                if isinstance(block, DatetimeTZBlock) and block.values.ndim == 1:
                    # TODO(2.0): remove once fastparquet no longer needs this
                    warnings.warn(
                        "In a future version, the BlockManager constructor "
                        "will assume that a DatetimeTZBlock with block.ndim==2 "
                        "has block.values.ndim == 2.",
                        DeprecationWarning,
                        stacklevel=find_stack_level(),
                    )

                    # error: Incompatible types in assignment (expression has type
                    # "Union[ExtensionArray, ndarray]", variable has type
                    # "DatetimeArray")
                    block.values = ensure_block_shape(  # type: ignore[assignment]
                        block.values, self.ndim
                    )
                    try:
                        block._cache.clear()
                    except AttributeError:
                        # _cache not initialized
                        pass

            self._verify_integrity()

    def _verify_integrity(self) -> None:
        mgr_shape = self.shape
        tot_items = sum(len(x.mgr_locs) for x in self.blocks)
        for block in self.blocks:
            if block.shape[1:] != mgr_shape[1:]:
                raise construction_error(tot_items, block.shape[1:], self.axes)
        if len(self.items) != tot_items:
            raise AssertionError(
                "Number of manager items must equal union of "
                f"block items\n# manager items: {len(self.items)}, # "
                f"tot_items: {tot_items}"
            )

    @classmethod
    def from_blocks(cls, blocks: list[Block], axes: list[Index]) -> BlockManager:
        """
        Constructor for BlockManager and SingleBlockManager with same signature.
        """
        return cls(blocks, axes, verify_integrity=False)

    # ----------------------------------------------------------------
    # Indexing

    def fast_xs(self, loc: int) -> ArrayLike:
        """
        Return the array corresponding to `frame.iloc[loc]`.

        Parameters
        ----------
        loc : int

        Returns
        -------
        np.ndarray or ExtensionArray
        """
        if len(self.blocks) == 1:
            return self.blocks[0].iget((slice(None), loc))

        dtype = interleaved_dtype([blk.dtype for blk in self.blocks])

        n = len(self)
        if isinstance(dtype, ExtensionDtype):
            cls = dtype.construct_array_type()
            result = cls._empty((n,), dtype=dtype)
        else:
            result = np.empty(n, dtype=dtype)
            result = ensure_wrapped_if_datetimelike(result)

        for blk in self.blocks:
            # Such assignment may incorrectly coerce NaT to None
            # result[blk.mgr_locs] = blk._slice((slice(None), loc))
            for i, rl in enumerate(blk.mgr_locs):
                result[rl] = blk.iget((i, loc))

        return result

    def iget(self, i: int) -> SingleBlockManager:
        """
        Return the data as a SingleBlockManager.
        """
        block = self.blocks[self.blknos[i]]
        values = block.iget(self.blklocs[i])

        # shortcut for select a single-dim from a 2-dim BM
        bp = BlockPlacement(slice(0, len(values)))
        nb = type(block)(values, placement=bp, ndim=1)
        return SingleBlockManager(nb, self.axes[1])

    def iget_values(self, i: int) -> ArrayLike:
        """
        Return the data for column i as the values (ndarray or ExtensionArray).
        """
        block = self.blocks[self.blknos[i]]
        values = block.iget(self.blklocs[i])
        return values

    @property
    def column_arrays(self) -> list[np.ndarray]:
        """
        Used in the JSON C code to access column arrays.
        This optimizes compared to using `iget_values` by converting each
        """
        # This is an optimized equivalent to
        #  result = [self.iget_values(i) for i in range(len(self.items))]
        result: list[np.ndarray | None] = [None] * len(self.items)

        for blk in self.blocks:
            mgr_locs = blk._mgr_locs
            values = blk.values_for_json()
            if values.ndim == 1:
                # TODO(EA2D): special casing not needed with 2D EAs
                result[mgr_locs[0]] = values

            else:
                for i, loc in enumerate(mgr_locs):
                    result[loc] = values[i]

        # error: Incompatible return value type (got "List[None]",
        # expected "List[ndarray[Any, Any]]")
        return result  # type: ignore[return-value]

    def iset(
        self, loc: int | slice | np.ndarray, value: ArrayLike, inplace: bool = False
    ):
        """
        Set new item in-place. Does not consolidate. Adds new Block if not
        contained in the current set of items
        """

        # FIXME: refactor, clearly separate broadcasting & zip-like assignment
        #        can prob also fix the various if tests for sparse/categorical
        if self._blklocs is None and self.ndim > 1:
            self._rebuild_blknos_and_blklocs()

        # Note: we exclude DTA/TDA here
        value_is_extension_type = is_1d_only_ea_dtype(value.dtype)
        if not value_is_extension_type:
            if value.ndim == 2:
                value = value.T
            else:
                value = ensure_block_shape(value, ndim=2)

            if value.shape[1:] != self.shape[1:]:
                raise AssertionError(
                    "Shape of new values must be compatible with manager shape"
                )

        if lib.is_integer(loc):
            # We have 6 tests where loc is _not_ an int.
            # In this case, get_blkno_placements will yield only one tuple,
            #  containing (self._blknos[loc], BlockPlacement(slice(0, 1, 1)))

            # Check if we can use _iset_single fastpath
            loc = cast(int, loc)
            blkno = self.blknos[loc]
            blk = self.blocks[blkno]
            if len(blk._mgr_locs) == 1:  # TODO: fastest way to check this?
                return self._iset_single(
                    loc,
                    value,
                    inplace=inplace,
                    blkno=blkno,
                    blk=blk,
                )

            # error: Incompatible types in assignment (expression has type
            # "List[Union[int, slice, ndarray]]", variable has type "Union[int,
            # slice, ndarray]")
            loc = [loc]  # type: ignore[assignment]

        # categorical/sparse/datetimetz
        if value_is_extension_type:

            def value_getitem(placement):
                return value

        else:

            def value_getitem(placement):
                return value[placement.indexer]

        # Accessing public blknos ensures the public versions are initialized
        blknos = self.blknos[loc]
        blklocs = self.blklocs[loc].copy()

        unfit_mgr_locs = []
        unfit_val_locs = []
        removed_blknos = []
        for blkno_l, val_locs in libinternals.get_blkno_placements(blknos, group=True):
            blk = self.blocks[blkno_l]
            blk_locs = blklocs[val_locs.indexer]
            if inplace and blk.should_store(value):
                blk.set_inplace(blk_locs, value_getitem(val_locs))
            else:
                unfit_mgr_locs.append(blk.mgr_locs.as_array[blk_locs])
                unfit_val_locs.append(val_locs)

                # If all block items are unfit, schedule the block for removal.
                if len(val_locs) == len(blk.mgr_locs):
                    removed_blknos.append(blkno_l)
                else:
                    blk.delete(blk_locs)
                    self._blklocs[blk.mgr_locs.indexer] = np.arange(len(blk))

        if len(removed_blknos):
            # Remove blocks & update blknos accordingly
            is_deleted = np.zeros(self.nblocks, dtype=np.bool_)
            is_deleted[removed_blknos] = True

            new_blknos = np.empty(self.nblocks, dtype=np.intp)
            new_blknos.fill(-1)
            new_blknos[~is_deleted] = np.arange(self.nblocks - len(removed_blknos))
            self._blknos = new_blknos[self._blknos]
            self.blocks = tuple(
                blk for i, blk in enumerate(self.blocks) if i not in set(removed_blknos)
            )

        if unfit_val_locs:
            unfit_idxr = np.concatenate(unfit_mgr_locs)
            unfit_count = len(unfit_idxr)

            new_blocks: list[Block] = []
            if value_is_extension_type:
                # This code (ab-)uses the fact that EA blocks contain only
                # one item.
                # TODO(EA2D): special casing unnecessary with 2D EAs
                new_blocks.extend(
                    new_block_2d(
                        values=value,
                        placement=BlockPlacement(slice(mgr_loc, mgr_loc + 1)),
                    )
                    for mgr_loc in unfit_idxr
                )

                self._blknos[unfit_idxr] = np.arange(unfit_count) + len(self.blocks)
                self._blklocs[unfit_idxr] = 0

            else:
                # unfit_val_locs contains BlockPlacement objects
                unfit_val_items = unfit_val_locs[0].append(unfit_val_locs[1:])

                new_blocks.append(
                    new_block_2d(
                        values=value_getitem(unfit_val_items),
                        placement=BlockPlacement(unfit_idxr),
                    )
                )

                self._blknos[unfit_idxr] = len(self.blocks)
                self._blklocs[unfit_idxr] = np.arange(unfit_count)

            self.blocks += tuple(new_blocks)

            # Newly created block's dtype may already be present.
            self._known_consolidated = False

    def _iset_single(
        self, loc: int, value: ArrayLike, inplace: bool, blkno: int, blk: Block
    ) -> None:
        """
        Fastpath for iset when we are only setting a single position and
        the Block currently in that position is itself single-column.

        In this case we can swap out the entire Block and blklocs and blknos
        are unaffected.
        """
        # Caller is responsible for verifying value.shape

        if inplace and blk.should_store(value):
            iloc = self.blklocs[loc]
            blk.set_inplace(slice(iloc, iloc + 1), value)
            return

        nb = new_block_2d(value, placement=blk._mgr_locs)
        old_blocks = self.blocks
        new_blocks = old_blocks[:blkno] + (nb,) + old_blocks[blkno + 1 :]
        self.blocks = new_blocks
        return

    def insert(self, loc: int, item: Hashable, value: ArrayLike) -> None:
        """
        Insert item at selected position.

        Parameters
        ----------
        loc : int
        item : hashable
        value : np.ndarray or ExtensionArray
        """
        # insert to the axis; this could possibly raise a TypeError
        new_axis = self.items.insert(loc, item)

        if value.ndim == 2:
            value = value.T
            if len(value) > 1:
                raise ValueError(
                    f"Expected a 1D array, got an array with shape {value.T.shape}"
                )
        else:
            value = ensure_block_shape(value, ndim=self.ndim)

        bp = BlockPlacement(slice(loc, loc + 1))
        block = new_block_2d(values=value, placement=bp)

        if not len(self.blocks):
            # Fastpath
            self._blklocs = np.array([0], dtype=np.intp)
            self._blknos = np.array([0], dtype=np.intp)
        else:
            self._insert_update_mgr_locs(loc)
            self._insert_update_blklocs_and_blknos(loc)

        self.axes[0] = new_axis
        self.blocks += (block,)

        self._known_consolidated = False

        if sum(not block.is_extension for block in self.blocks) > 100:
            warnings.warn(
                "DataFrame is highly fragmented.  This is usually the result "
                "of calling `frame.insert` many times, which has poor performance.  "
                "Consider joining all columns at once using pd.concat(axis=1) "
                "instead. To get a de-fragmented frame, use `newframe = frame.copy()`",
                PerformanceWarning,
                stacklevel=find_stack_level(),
            )

    def _insert_update_mgr_locs(self, loc) -> None:
        """
        When inserting a new Block at location 'loc', we increment
        all of the mgr_locs of blocks above that by one.
        """
        for blkno, count in _fast_count_smallints(self.blknos[loc:]):
            # .620 this way, .326 of which is in increment_above
            blk = self.blocks[blkno]
            blk._mgr_locs = blk._mgr_locs.increment_above(loc)

    def _insert_update_blklocs_and_blknos(self, loc) -> None:
        """
        When inserting a new Block at location 'loc', we update our
        _blklocs and _blknos.
        """

        # Accessing public blklocs ensures the public versions are initialized
        if loc == self.blklocs.shape[0]:
            # np.append is a lot faster, let's use it if we can.
            self._blklocs = np.append(self._blklocs, 0)
            self._blknos = np.append(self._blknos, len(self.blocks))
        elif loc == 0:
            # np.append is a lot faster, let's use it if we can.
            self._blklocs = np.append(self._blklocs[::-1], 0)[::-1]
            self._blknos = np.append(self._blknos[::-1], len(self.blocks))[::-1]
        else:
            new_blklocs, new_blknos = libinternals.update_blklocs_and_blknos(
                self.blklocs, self.blknos, loc, len(self.blocks)
            )
            self._blklocs = new_blklocs
            self._blknos = new_blknos

    def idelete(self, indexer) -> BlockManager:
        """
        Delete selected locations, returning a new BlockManager.
        """
        is_deleted = np.zeros(self.shape[0], dtype=np.bool_)
        is_deleted[indexer] = True
        taker = (~is_deleted).nonzero()[0]

        nbs = self._slice_take_blocks_ax0(taker, only_slice=True)
        new_columns = self.items[~is_deleted]
        axes = [new_columns, self.axes[1]]
        return type(self)(tuple(nbs), axes, verify_integrity=False)

    # ----------------------------------------------------------------
    # Block-wise Operation

    def grouped_reduce(self: T, func: Callable, ignore_failures: bool = False) -> T:
        """
        Apply grouped reduction function blockwise, returning a new BlockManager.

        Parameters
        ----------
        func : grouped reduction function
        ignore_failures : bool, default False
            Whether to drop blocks where func raises TypeError.

        Returns
        -------
        BlockManager
        """
        result_blocks: list[Block] = []
        dropped_any = False

        for blk in self.blocks:
            if blk.is_object:
                # split on object-dtype blocks bc some columns may raise
                #  while others do not.
                for sb in blk._split():
                    try:
                        applied = sb.apply(func)
                    except (TypeError, NotImplementedError):
                        if not ignore_failures:
                            raise
                        dropped_any = True
                        continue
                    result_blocks = extend_blocks(applied, result_blocks)
            else:
                try:
                    applied = blk.apply(func)
                except (TypeError, NotImplementedError):
                    if not ignore_failures:
                        raise
                    dropped_any = True
                    continue
                result_blocks = extend_blocks(applied, result_blocks)

        if len(result_blocks) == 0:
            index = Index([None])  # placeholder
        else:
            index = Index(range(result_blocks[0].values.shape[-1]))

        if dropped_any:
            # faster to skip _combine if we haven't dropped any blocks
            return self._combine(result_blocks, copy=False, index=index)

        return type(self).from_blocks(result_blocks, [self.axes[0], index])

    def reduce(
        self: T, func: Callable, ignore_failures: bool = False
    ) -> tuple[T, np.ndarray]:
        """
        Apply reduction function blockwise, returning a single-row BlockManager.

        Parameters
        ----------
        func : reduction function
        ignore_failures : bool, default False
            Whether to drop blocks where func raises TypeError.

        Returns
        -------
        BlockManager
        np.ndarray
            Indexer of mgr_locs that are retained.
        """
        # If 2D, we assume that we're operating column-wise
        assert self.ndim == 2

        res_blocks: list[Block] = []
        for blk in self.blocks:
            nbs = blk.reduce(func, ignore_failures)
            res_blocks.extend(nbs)

        index = Index([None])  # placeholder
        if ignore_failures:
            if res_blocks:
                indexer = np.concatenate([blk.mgr_locs.as_array for blk in res_blocks])
                new_mgr = self._combine(res_blocks, copy=False, index=index)
            else:
                indexer = []
                new_mgr = type(self).from_blocks([], [self.items[:0], index])
        else:
            indexer = np.arange(self.shape[0])
            new_mgr = type(self).from_blocks(res_blocks, [self.items, index])
        return new_mgr, indexer

    def operate_blockwise(self, other: BlockManager, array_op) -> BlockManager:
        """
        Apply array_op blockwise with another (aligned) BlockManager.
        """
        return operate_blockwise(self, other, array_op)

    def _equal_values(self: BlockManager, other: BlockManager) -> bool:
        """
        Used in .equals defined in base class. Only check the column values
        assuming shape and indexes have already been checked.
        """
        return blockwise_all(self, other, array_equals)

    def quantile(
        self: T,
        *,
        qs: Float64Index,
        axis: int = 0,
        interpolation="linear",
    ) -> T:
        """
        Iterate over blocks applying quantile reduction.
        This routine is intended for reduction type operations and
        will do inference on the generated blocks.

        Parameters
        ----------
        axis: reduction axis, default 0
        consolidate: bool, default True. Join together blocks having same
            dtype
        interpolation : type of interpolation, default 'linear'
        qs : list of the quantiles to be computed

        Returns
        -------
        BlockManager
        """
        # Series dispatches to DataFrame for quantile, which allows us to
        #  simplify some of the code here and in the blocks
        assert self.ndim >= 2
        assert is_list_like(qs)  # caller is responsible for this
        assert axis == 1  # only ever called this way

        new_axes = list(self.axes)
        new_axes[1] = Float64Index(qs)

        blocks = [
            blk.quantile(axis=axis, qs=qs, interpolation=interpolation)
            for blk in self.blocks
        ]

        return type(self)(blocks, new_axes)

    # ----------------------------------------------------------------

    def unstack(self, unstacker, fill_value) -> BlockManager:
        """
        Return a BlockManager with all blocks unstacked.

        Parameters
        ----------
        unstacker : reshape._Unstacker
        fill_value : Any
            fill_value for newly introduced missing values.

        Returns
        -------
        unstacked : BlockManager
        """
        new_columns = unstacker.get_new_columns(self.items)
        new_index = unstacker.new_index

        allow_fill = not unstacker.mask_all
        if allow_fill:
            # calculating the full mask once and passing it to Block._unstack is
            #  faster than letting calculating it in each repeated call
            new_mask2D = (~unstacker.mask).reshape(*unstacker.full_shape)
            needs_masking = new_mask2D.any(axis=0)
        else:
            needs_masking = np.zeros(unstacker.full_shape[1], dtype=bool)

        new_blocks: list[Block] = []
        columns_mask: list[np.ndarray] = []

        if len(self.items) == 0:
            factor = 1
        else:
            fac = len(new_columns) / len(self.items)
            assert fac == int(fac)
            factor = int(fac)

        for blk in self.blocks:
            mgr_locs = blk.mgr_locs
            new_placement = mgr_locs.tile_for_unstack(factor)

            blocks, mask = blk._unstack(
                unstacker,
                fill_value,
                new_placement=new_placement,
                needs_masking=needs_masking,
            )

            new_blocks.extend(blocks)
            columns_mask.extend(mask)

            # Block._unstack should ensure this holds,
            assert mask.sum() == sum(len(nb._mgr_locs) for nb in blocks)
            # In turn this ensures that in the BlockManager call below
            #  we have len(new_columns) == sum(x.shape[0] for x in new_blocks)
            #  which suffices to allow us to pass verify_inegrity=False

        new_columns = new_columns[columns_mask]

        bm = BlockManager(new_blocks, [new_columns, new_index], verify_integrity=False)
        return bm

    def to_dict(self, copy: bool = True):
        """
        Return a dict of str(dtype) -> BlockManager

        Parameters
        ----------
        copy : bool, default True

        Returns
        -------
        values : a dict of dtype -> BlockManager
        """

        bd: dict[str, list[Block]] = {}
        for b in self.blocks:
            bd.setdefault(str(b.dtype), []).append(b)

        # TODO(EA2D): the combine will be unnecessary with 2D EAs
        return {dtype: self._combine(blocks, copy=copy) for dtype, blocks in bd.items()}

    def as_array(
        self,
        dtype: np.dtype | None = None,
        copy: bool = False,
        na_value=lib.no_default,
    ) -> np.ndarray:
        """
        Convert the blockmanager data into an numpy array.

        Parameters
        ----------
        dtype : np.dtype or None, default None
            Data type of the return array.
        copy : bool, default False
            If True then guarantee that a copy is returned. A value of
            False does not guarantee that the underlying data is not
            copied.
        na_value : object, default lib.no_default
            Value to be used as the missing value sentinel.

        Returns
        -------
        arr : ndarray
        """
        if len(self.blocks) == 0:
            arr = np.empty(self.shape, dtype=float)
            return arr.transpose()

        # We want to copy when na_value is provided to avoid
        # mutating the original object
        copy = copy or na_value is not lib.no_default

        if self.is_single_block:
            blk = self.blocks[0]
            if blk.is_extension:
                # Avoid implicit conversion of extension blocks to object

                # error: Item "ndarray" of "Union[ndarray, ExtensionArray]" has no
                # attribute "to_numpy"
                arr = blk.values.to_numpy(  # type: ignore[union-attr]
                    dtype=dtype,
                    na_value=na_value,
                ).reshape(blk.shape)
            else:
                arr = np.asarray(blk.get_values())
                if dtype:
                    arr = arr.astype(dtype, copy=False)
        else:
            arr = self._interleave(dtype=dtype, na_value=na_value)
            # The underlying data was copied within _interleave
            copy = False

        if copy:
            arr = arr.copy()

        if na_value is not lib.no_default:
            arr[isna(arr)] = na_value

        return arr.transpose()

    def _interleave(
        self,
        dtype: np.dtype | None = None,
        na_value=lib.no_default,
    ) -> np.ndarray:
        """
        Return ndarray from blocks with specified item order
        Items must be contained in the blocks
        """
        if not dtype:
            # Incompatible types in assignment (expression has type
            # "Optional[Union[dtype[Any], ExtensionDtype]]", variable has
            # type "Optional[dtype[Any]]")
            dtype = interleaved_dtype(  # type: ignore[assignment]
                [blk.dtype for blk in self.blocks]
            )

        # TODO: https://github.com/pandas-dev/pandas/issues/22791
        # Give EAs some input on what happens here. Sparse needs this.
        if isinstance(dtype, SparseDtype):
            dtype = dtype.subtype
            dtype = cast(np.dtype, dtype)
        elif isinstance(dtype, ExtensionDtype):
            dtype = np.dtype("object")
        elif is_dtype_equal(dtype, str):
            dtype = np.dtype("object")

        result = np.empty(self.shape, dtype=dtype)

        itemmask = np.zeros(self.shape[0])

        if dtype == np.dtype("object") and na_value is lib.no_default:
            # much more performant than using to_numpy below
            for blk in self.blocks:
                rl = blk.mgr_locs
                arr = blk.get_values(dtype)
                result[rl.indexer] = arr
                itemmask[rl.indexer] = 1
            return result

        for blk in self.blocks:
            rl = blk.mgr_locs
            if blk.is_extension:
                # Avoid implicit conversion of extension blocks to object

                # error: Item "ndarray" of "Union[ndarray, ExtensionArray]" has no
                # attribute "to_numpy"
                arr = blk.values.to_numpy(  # type: ignore[union-attr]
                    dtype=dtype,
                    na_value=na_value,
                )
            else:
                arr = blk.get_values(dtype)
            result[rl.indexer] = arr
            itemmask[rl.indexer] = 1

        if not itemmask.all():
            raise AssertionError("Some items were not contained in blocks")

        return result

    # ----------------------------------------------------------------
    # Consolidation

    def is_consolidated(self) -> bool:
        """
        Return True if more than one block with the same dtype
        """
        if not self._known_consolidated:
            self._consolidate_check()
        return self._is_consolidated

    def _consolidate_check(self) -> None:
        if len(self.blocks) == 1:
            # fastpath
            self._is_consolidated = True
            self._known_consolidated = True
            return
        dtypes = [blk.dtype for blk in self.blocks if blk._can_consolidate]
        self._is_consolidated = len(dtypes) == len(set(dtypes))
        self._known_consolidated = True

    def _consolidate_inplace(self) -> None:
        if not self.is_consolidated():
            self.blocks = tuple(_consolidate(self.blocks))
            self._is_consolidated = True
            self._known_consolidated = True
            self._rebuild_blknos_and_blklocs()


class SingleBlockManager(BaseBlockManager, SingleDataManager):
    """manage a single block with"""

    ndim = 1
    _is_consolidated = True
    _known_consolidated = True
    __slots__ = ()
    is_single_block = True

    def __init__(
        self,
        block: Block,
        axis: Index,
        verify_integrity: bool = False,
        fastpath=lib.no_default,
    ):
        # Assertions disabled for performance
        # assert isinstance(block, Block), type(block)
        # assert isinstance(axis, Index), type(axis)

        if fastpath is not lib.no_default:
            warnings.warn(
                "The `fastpath` keyword is deprecated and will be removed "
                "in a future version.",
                FutureWarning,
                stacklevel=find_stack_level(),
            )

        self.axes = [axis]
        self.blocks = (block,)

    @classmethod
    def from_blocks(cls, blocks: list[Block], axes: list[Index]) -> SingleBlockManager:
        """
        Constructor for BlockManager and SingleBlockManager with same signature.
        """
        assert len(blocks) == 1
        assert len(axes) == 1
        return cls(blocks[0], axes[0], verify_integrity=False)

    @classmethod
    def from_array(cls, array: ArrayLike, index: Index) -> SingleBlockManager:
        """
        Constructor for if we have an array that is not yet a Block.
        """
        block = new_block(array, placement=slice(0, len(index)), ndim=1)
        return cls(block, index)

    def to_2d_mgr(self, columns: Index) -> BlockManager:
        """
        Manager analogue of Series.to_frame
        """
        blk = self.blocks[0]
        arr = ensure_block_shape(blk.values, ndim=2)
        bp = BlockPlacement(0)
        new_blk = type(blk)(arr, placement=bp, ndim=2)
        axes = [columns, self.axes[0]]
        return BlockManager([new_blk], axes=axes, verify_integrity=False)

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_items = [self.items[b.mgr_locs.indexer] for b in self.blocks]
        axes_array = list(self.axes)

        extra_state = {
            "0.14.1": {
                "axes": axes_array,
                "blocks": [
                    {"values": b.values, "mgr_locs": b.mgr_locs.indexer}
                    for b in self.blocks
                ],
            }
        }

        # First three elements of the state are to maintain forward
        # compatibility with 0.13.1.
        return axes_array, block_values, block_items, extra_state

    def __setstate__(self, state):
        def unpickle_block(values, mgr_locs, ndim: int) -> Block:
            # TODO(EA2D): ndim would be unnecessary with 2D EAs
            # older pickles may store e.g. DatetimeIndex instead of DatetimeArray
            values = extract_array(values, extract_numpy=True)
            return new_block(values, placement=mgr_locs, ndim=ndim)

        if isinstance(state, tuple) and len(state) >= 4 and "0.14.1" in state[3]:
            state = state[3]["0.14.1"]
            self.axes = [ensure_index(ax) for ax in state["axes"]]
            ndim = len(self.axes)
            self.blocks = tuple(
                unpickle_block(b["values"], b["mgr_locs"], ndim=ndim)
                for b in state["blocks"]
            )
        else:
            raise NotImplementedError("pre-0.14.1 pickles are no longer supported")

        self._post_setstate()

    def _post_setstate(self):
        pass

    @cache_readonly
    def _block(self) -> Block:
        return self.blocks[0]

    @property
    def _blknos(self):
        """compat with BlockManager"""
        return None

    @property
    def _blklocs(self):
        """compat with BlockManager"""
        return None

    def getitem_mgr(self, indexer) -> SingleBlockManager:
        # similar to get_slice, but not restricted to slice indexer
        blk = self._block
        array = blk._slice(indexer)
        if array.ndim > 1:
            # This will be caught by Series._get_values
            raise ValueError("dimension-expanding indexing not allowed")

        bp = BlockPlacement(slice(0, len(array)))
        block = type(blk)(array, placement=bp, ndim=1)

        new_idx = self.index[indexer]
        return type(self)(block, new_idx)

    def get_slice(self, slobj: slice, axis: int = 0) -> SingleBlockManager:
        # Assertion disabled for performance
        # assert isinstance(slobj, slice), type(slobj)
        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

        blk = self._block
        array = blk._slice(slobj)
        bp = BlockPlacement(slice(0, len(array)))
        block = type(blk)(array, placement=bp, ndim=1)
        new_index = self.index._getitem_slice(slobj)
        return type(self)(block, new_index)

    @property
    def index(self) -> Index:
        return self.axes[0]

    @property
    def dtype(self) -> DtypeObj:
        return self._block.dtype

    def get_dtypes(self) -> np.ndarray:
        return np.array([self._block.dtype])

    def external_values(self):
        """The array that Series.values returns"""
        return self._block.external_values()

    def internal_values(self):
        """The array that Series._values returns"""
        return self._block.values

    def array_values(self):
        """The array that Series.array returns"""
        return self._block.array_values

    def get_numeric_data(self, copy: bool = False):
        if self._block.is_numeric:
            if copy:
                return self.copy()
            return self
        return self.make_empty()

    @property
    def _can_hold_na(self) -> bool:
        return self._block._can_hold_na

    def idelete(self, indexer) -> SingleBlockManager:
        """
        Delete single location from SingleBlockManager.

        Ensures that self.blocks doesn't become empty.
        """
        self._block.delete(indexer)
        self.axes[0] = self.axes[0].delete(indexer)
        return self

    def fast_xs(self, loc):
        """
        fast path for getting a cross-section
        return a view of the data
        """
        raise NotImplementedError("Use series._values[loc] instead")

    def set_values(self, values: ArrayLike):
        """
        Set the values of the single block in place.

        Use at your own risk! This does not check if the passed values are
        valid for the current Block/SingleBlockManager (length, dtype, etc).
        """
        self.blocks[0].values = values
        self.blocks[0]._mgr_locs = BlockPlacement(slice(len(values)))

    def _equal_values(self: T, other: T) -> bool:
        """
        Used in .equals defined in base class. Only check the column values
        assuming shape and indexes have already been checked.
        """
        # For SingleBlockManager (i.e.Series)
        if other.ndim != 1:
            return False
        left = self.blocks[0].values
        right = other.blocks[0].values
        return array_equals(left, right)


# --------------------------------------------------------------------
# Constructor Helpers


def create_block_manager_from_blocks(
    blocks: list[Block],
    axes: list[Index],
    consolidate: bool = True,
    verify_integrity: bool = True,
) -> BlockManager:
    # If verify_integrity=False, then caller is responsible for checking
    #  all(x.shape[-1] == len(axes[1]) for x in blocks)
    #  sum(x.shape[0] for x in blocks) == len(axes[0])
    #  set(x for for blk in blocks for x in blk.mgr_locs) == set(range(len(axes[0])))
    #  all(blk.ndim == 2 for blk in blocks)
    # This allows us to safely pass verify_integrity=False

    try:
        mgr = BlockManager(blocks, axes, verify_integrity=verify_integrity)

    except ValueError as err:
        arrays = [blk.values for blk in blocks]
        tot_items = sum(arr.shape[0] for arr in arrays)
        raise construction_error(tot_items, arrays[0].shape[1:], axes, err)

    if consolidate:
        mgr._consolidate_inplace()
    return mgr


def create_block_manager_from_column_arrays(
    arrays: list[ArrayLike],
    axes: list[Index],
    consolidate: bool = True,
) -> BlockManager:
    # Assertions disabled for performance (caller is responsible for verifying)
    # assert isinstance(axes, list)
    # assert all(isinstance(x, Index) for x in axes)
    # assert all(isinstance(x, (np.ndarray, ExtensionArray)) for x in arrays)
    # assert all(type(x) is not PandasArray for x in arrays)
    # assert all(x.ndim == 1 for x in arrays)
    # assert all(len(x) == len(axes[1]) for x in arrays)
    # assert len(arrays) == len(axes[0])
    # These last three are sufficient to allow us to safely pass
    #  verify_integrity=False below.

    try:
        blocks = _form_blocks(arrays, consolidate)
        mgr = BlockManager(blocks, axes, verify_integrity=False)
    except ValueError as e:
        raise construction_error(len(arrays), arrays[0].shape, axes, e)
    if consolidate:
        mgr._consolidate_inplace()
    return mgr


def construction_error(
    tot_items: int,
    block_shape: Shape,
    axes: list[Index],
    e: ValueError | None = None,
):
    """raise a helpful message about our construction"""
    passed = tuple(map(int, [tot_items] + list(block_shape)))
    # Correcting the user facing error message during dataframe construction
    if len(passed) <= 2:
        passed = passed[::-1]

    implied = tuple(len(ax) for ax in axes)
    # Correcting the user facing error message during dataframe construction
    if len(implied) <= 2:
        implied = implied[::-1]

    # We return the exception object instead of raising it so that we
    #  can raise it in the caller; mypy plays better with that
    if passed == implied and e is not None:
        return e
    if block_shape[0] == 0:
        return ValueError("Empty data passed with indices specified.")
    return ValueError(f"Shape of passed values is {passed}, indices imply {implied}")


# -----------------------------------------------------------------------


def _grouping_func(tup: tuple[int, ArrayLike]) -> tuple[int, bool, DtypeObj]:
    # compat for numpy<1.21, in which comparing a np.dtype with an ExtensionDtype
    # raises instead of returning False. Once earlier numpy versions are dropped,
    # this can be simplified to `return tup[1].dtype`
    dtype = tup[1].dtype

    if is_1d_only_ea_dtype(dtype):
        # We know these won't be consolidated, so don't need to group these.
        # This avoids expensive comparisons of CategoricalDtype objects
        sep = id(dtype)
    else:
        sep = 0

    return sep, isinstance(dtype, np.dtype), dtype


def _form_blocks(arrays: list[ArrayLike], consolidate: bool) -> list[Block]:
    tuples = list(enumerate(arrays))

    if not consolidate:
        nbs = _tuples_to_blocks_no_consolidate(tuples)
        return nbs

    # group by dtype
    grouper = itertools.groupby(tuples, _grouping_func)

    nbs = []
    for (_, _, dtype), tup_block in grouper:
        block_type = get_block_type(dtype)

        if isinstance(dtype, np.dtype):
            is_dtlike = dtype.kind in ["m", "M"]

            if issubclass(dtype.type, (str, bytes)):
                dtype = np.dtype(object)

            values, placement = _stack_arrays(list(tup_block), dtype)
            if is_dtlike:
                values = ensure_wrapped_if_datetimelike(values)
            blk = block_type(values, placement=BlockPlacement(placement), ndim=2)
            nbs.append(blk)

        elif is_1d_only_ea_dtype(dtype):
            dtype_blocks = [
                block_type(x[1], placement=BlockPlacement(x[0]), ndim=2)
                for x in tup_block
            ]
            nbs.extend(dtype_blocks)

        else:
            dtype_blocks = [
                block_type(
                    ensure_block_shape(x[1], 2), placement=BlockPlacement(x[0]), ndim=2
                )
                for x in tup_block
            ]
            nbs.extend(dtype_blocks)
    return nbs


def _tuples_to_blocks_no_consolidate(tuples) -> list[Block]:
    # tuples produced within _form_blocks are of the form (placement, array)
    return [
        new_block_2d(ensure_block_shape(x[1], ndim=2), placement=BlockPlacement(x[0]))
        for x in tuples
    ]


def _stack_arrays(tuples, dtype: np.dtype):

    placement, arrays = zip(*tuples)

    first = arrays[0]
    shape = (len(arrays),) + first.shape

    stacked = np.empty(shape, dtype=dtype)
    for i, arr in enumerate(arrays):
        stacked[i] = arr

    return stacked, placement


def _consolidate(blocks: tuple[Block, ...]) -> list[Block]:
    """
    Merge blocks having same dtype, exclude non-consolidating blocks
    """
    # sort by _can_consolidate, dtype
    gkey = lambda x: x._consolidate_key
    grouper = itertools.groupby(sorted(blocks, key=gkey), gkey)

    new_blocks: list[Block] = []
    for (_can_consolidate, dtype), group_blocks in grouper:
        merged_blocks = _merge_blocks(
            list(group_blocks), dtype=dtype, can_consolidate=_can_consolidate
        )
        new_blocks = extend_blocks(merged_blocks, new_blocks)
    return new_blocks


def _merge_blocks(
    blocks: list[Block], dtype: DtypeObj, can_consolidate: bool
) -> list[Block]:

    if len(blocks) == 1:
        return blocks

    if can_consolidate:

        # TODO: optimization potential in case all mgrs contain slices and
        # combination of those slices is a slice, too.
        new_mgr_locs = np.concatenate([b.mgr_locs.as_array for b in blocks])

        new_values: ArrayLike

        if isinstance(blocks[0].dtype, np.dtype):
            # error: List comprehension has incompatible type List[Union[ndarray,
            # ExtensionArray]]; expected List[Union[complex, generic,
            # Sequence[Union[int, float, complex, str, bytes, generic]],
            # Sequence[Sequence[Any]], SupportsArray]]
            new_values = np.vstack([b.values for b in blocks])  # type: ignore[misc]
        else:
            bvals = [blk.values for blk in blocks]
            bvals2 = cast(Sequence[NDArrayBackedExtensionArray], bvals)
            new_values = bvals2[0]._concat_same_type(bvals2, axis=0)

        argsort = np.argsort(new_mgr_locs)
        new_values = new_values[argsort]
        new_mgr_locs = new_mgr_locs[argsort]

        bp = BlockPlacement(new_mgr_locs)
        return [new_block_2d(new_values, placement=bp)]

    # can't consolidate --> no merge
    return blocks


def _fast_count_smallints(arr: npt.NDArray[np.intp]):
    """Faster version of set(arr) for sequences of small numbers."""
    counts = np.bincount(arr)
    nz = counts.nonzero()[0]
    # Note: list(zip(...) outperforms list(np.c_[nz, counts[nz]]) here,
    #  in one benchmark by a factor of 11
    return zip(nz, counts[nz])


def _preprocess_slice_or_indexer(
    slice_or_indexer: slice | np.ndarray, length: int, allow_fill: bool
):
    if isinstance(slice_or_indexer, slice):
        return (
            "slice",
            slice_or_indexer,
            libinternals.slice_len(slice_or_indexer, length),
        )
    else:
        if (
            not isinstance(slice_or_indexer, np.ndarray)
            or slice_or_indexer.dtype.kind != "i"
        ):
            dtype = getattr(slice_or_indexer, "dtype", None)
            raise TypeError(type(slice_or_indexer), dtype)

        indexer = ensure_platform_int(slice_or_indexer)
        if not allow_fill:
            indexer = maybe_convert_indices(indexer, length)
        return "fancy", indexer, len(indexer)
