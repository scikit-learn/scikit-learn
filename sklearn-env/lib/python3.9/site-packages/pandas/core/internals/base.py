"""
Base class for the internal managers. Both BlockManager and ArrayManager
inherit from this class.
"""
from __future__ import annotations

from typing import (
    TypeVar,
    final,
)

from pandas._typing import (
    ArrayLike,
    DtypeObj,
    Shape,
)
from pandas.errors import AbstractMethodError

from pandas.core.dtypes.cast import find_common_type

from pandas.core.base import PandasObject
from pandas.core.indexes.api import (
    Index,
    default_index,
)

T = TypeVar("T", bound="DataManager")


class DataManager(PandasObject):

    # TODO share more methods/attributes

    axes: list[Index]

    @property
    def items(self) -> Index:
        raise AbstractMethodError(self)

    @final
    def __len__(self) -> int:
        return len(self.items)

    @property
    def ndim(self) -> int:
        return len(self.axes)

    @property
    def shape(self) -> Shape:
        return tuple(len(ax) for ax in self.axes)

    @final
    def _validate_set_axis(self, axis: int, new_labels: Index) -> None:
        # Caller is responsible for ensuring we have an Index object.
        old_len = len(self.axes[axis])
        new_len = len(new_labels)

        if axis == 1 and len(self.items) == 0:
            # If we are setting the index on a DataFrame with no columns,
            #  it is OK to change the length.
            pass

        elif new_len != old_len:
            raise ValueError(
                f"Length mismatch: Expected axis has {old_len} elements, new "
                f"values have {new_len} elements"
            )

    def reindex_indexer(
        self: T,
        new_axis,
        indexer,
        axis: int,
        fill_value=None,
        allow_dups: bool = False,
        copy: bool = True,
        consolidate: bool = True,
        only_slice: bool = False,
    ) -> T:
        raise AbstractMethodError(self)

    @final
    def reindex_axis(
        self: T,
        new_index: Index,
        axis: int,
        fill_value=None,
        consolidate: bool = True,
        only_slice: bool = False,
    ) -> T:
        """
        Conform data manager to new index.
        """
        new_index, indexer = self.axes[axis].reindex(new_index)

        return self.reindex_indexer(
            new_index,
            indexer,
            axis=axis,
            fill_value=fill_value,
            copy=False,
            consolidate=consolidate,
            only_slice=only_slice,
        )

    def _equal_values(self: T, other: T) -> bool:
        """
        To be implemented by the subclasses. Only check the column values
        assuming shape and indexes have already been checked.
        """
        raise AbstractMethodError(self)

    @final
    def equals(self, other: object) -> bool:
        """
        Implementation for DataFrame.equals
        """
        if not isinstance(other, DataManager):
            return False

        self_axes, other_axes = self.axes, other.axes
        if len(self_axes) != len(other_axes):
            return False
        if not all(ax1.equals(ax2) for ax1, ax2 in zip(self_axes, other_axes)):
            return False

        return self._equal_values(other)

    def apply(
        self: T,
        f,
        align_keys: list[str] | None = None,
        ignore_failures: bool = False,
        **kwargs,
    ) -> T:
        raise AbstractMethodError(self)

    @final
    def isna(self: T, func) -> T:
        return self.apply("apply", func=func)

    # --------------------------------------------------------------------
    # Consolidation: No-ops for all but BlockManager

    def is_consolidated(self) -> bool:
        return True

    def consolidate(self: T) -> T:
        return self

    def _consolidate_inplace(self) -> None:
        return


class SingleDataManager(DataManager):
    ndim = 1

    @final
    @property
    def array(self) -> ArrayLike:
        """
        Quick access to the backing array of the Block or SingleArrayManager.
        """
        # error: "SingleDataManager" has no attribute "arrays"; maybe "array"
        return self.arrays[0]  # type: ignore[attr-defined]

    def setitem_inplace(self, indexer, value) -> None:
        """
        Set values with indexer.

        For Single[Block/Array]Manager, this backs s[indexer] = value

        This is an inplace version of `setitem()`, mutating the manager/values
        in place, not returning a new Manager (and Block), and thus never changing
        the dtype.
        """
        self.array[indexer] = value

    def grouped_reduce(self, func, ignore_failures: bool = False):
        """
        ignore_failures : bool, default False
            Not used; for compatibility with ArrayManager/BlockManager.
        """

        arr = self.array
        res = func(arr)
        index = default_index(len(res))

        mgr = type(self).from_array(res, index)
        return mgr

    @classmethod
    def from_array(cls, arr: ArrayLike, index: Index):
        raise AbstractMethodError(cls)


def interleaved_dtype(dtypes: list[DtypeObj]) -> DtypeObj | None:
    """
    Find the common dtype for `blocks`.

    Parameters
    ----------
    blocks : List[DtypeObj]

    Returns
    -------
    dtype : np.dtype, ExtensionDtype, or None
        None is returned when `blocks` is empty.
    """
    if not len(dtypes):
        return None

    return find_common_type(dtypes)
