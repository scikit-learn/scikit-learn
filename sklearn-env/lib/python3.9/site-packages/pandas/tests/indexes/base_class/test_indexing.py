import numpy as np
import pytest

import pandas as pd
from pandas import (
    Index,
    NaT,
)
import pandas._testing as tm


class TestGetSliceBounds:
    @pytest.mark.parametrize("kind", ["getitem", "loc", None])
    @pytest.mark.parametrize("side, expected", [("left", 4), ("right", 5)])
    def test_get_slice_bounds_within(self, kind, side, expected):
        index = Index(list("abcdef"))
        with tm.assert_produces_warning(FutureWarning, match="'kind' argument"):
            result = index.get_slice_bound("e", kind=kind, side=side)
        assert result == expected

    @pytest.mark.parametrize("kind", ["getitem", "loc", None])
    @pytest.mark.parametrize("side", ["left", "right"])
    @pytest.mark.parametrize(
        "data, bound, expected", [(list("abcdef"), "x", 6), (list("bcdefg"), "a", 0)]
    )
    def test_get_slice_bounds_outside(self, kind, side, expected, data, bound):
        index = Index(data)
        with tm.assert_produces_warning(FutureWarning, match="'kind' argument"):
            result = index.get_slice_bound(bound, kind=kind, side=side)
        assert result == expected

    def test_get_slice_bounds_invalid_side(self):
        with pytest.raises(ValueError, match="Invalid value for side kwarg"):
            Index([]).get_slice_bound("a", side="middle")


class TestGetIndexerNonUnique:
    def test_get_indexer_non_unique_dtype_mismatch(self):
        # GH#25459
        indexes, missing = Index(["A", "B"]).get_indexer_non_unique(Index([0]))
        tm.assert_numpy_array_equal(np.array([-1], dtype=np.intp), indexes)
        tm.assert_numpy_array_equal(np.array([0], dtype=np.intp), missing)


class TestGetLoc:
    @pytest.mark.slow  # to_flat_index takes a while
    def test_get_loc_tuple_monotonic_above_size_cutoff(self):
        # Go through the libindex path for which using
        # _bin_search vs ndarray.searchsorted makes a difference

        lev = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        dti = pd.date_range("2016-01-01", periods=100)

        mi = pd.MultiIndex.from_product([lev, range(10 ** 3), dti])
        oidx = mi.to_flat_index()

        loc = len(oidx) // 2
        tup = oidx[loc]

        res = oidx.get_loc(tup)
        assert res == loc

    def test_get_loc_nan_object_dtype_nonmonotonic_nonunique(self):
        # case that goes through _maybe_get_bool_indexer
        idx = Index(["foo", np.nan, None, "foo", 1.0, None], dtype=object)

        # we dont raise KeyError on nan
        res = idx.get_loc(np.nan)
        assert res == 1

        # we only match on None, not on np.nan
        res = idx.get_loc(None)
        expected = np.array([False, False, True, False, False, True])
        tm.assert_numpy_array_equal(res, expected)

        # we don't match at all on mismatched NA
        with pytest.raises(KeyError, match="NaT"):
            idx.get_loc(NaT)
