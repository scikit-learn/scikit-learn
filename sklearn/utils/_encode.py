# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import Counter
from collections.abc import Iterable
from numbers import Real
from typing import NamedTuple

import numpy as np

from sklearn.utils._array_api import device, get_namespace, size
from sklearn.utils._missing import is_scalar_nan


def _unique(values, *, return_inverse=False, return_counts=False):
    """Helper function to find unique values with support for python objects.

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.

    Parameters
    ----------
    values : ndarray
        Values to find uniques from.

    return_inverse : bool, default=False
        If True, also return the indices of the unique values.

    return_counts : bool, default=False
        If True, also return the number of times each unique item appears in
        values.

    Returns
    -------
    unique : ndarray
        The sorted unique values.

    unique_inverse : ndarray
        The indices to reconstruct the original array from the unique array.
        Only provided if `return_inverse` is True.

    unique_counts : ndarray
        The number of times each of the unique values comes up in the original
        array. Only provided if `return_counts` is True.
    """
    if values.dtype == object:
        return _unique_python(
            values, return_inverse=return_inverse, return_counts=return_counts
        )
    # numerical
    return _unique_np(
        values, return_inverse=return_inverse, return_counts=return_counts
    )


def _unique_np(values, return_inverse=False, return_counts=False):
    """Helper function to find unique values for numpy arrays that correctly
    accounts for nans. See `_unique` documentation for details."""
    xp, _ = get_namespace(values)

    inverse, counts = None, None

    if return_inverse and return_counts:
        uniques, _, inverse, counts = xp.unique_all(values)
    elif return_inverse:
        uniques, inverse = xp.unique_inverse(values)
    elif return_counts:
        uniques, counts = xp.unique_counts(values)
    else:
        uniques = xp.unique_values(values)

    # np.unique will have duplicate missing values at the end of `uniques`
    # here we clip the nans and remove it from uniques
    if size(uniques) and is_scalar_nan(uniques[-1]):
        nan_idx = xp.searchsorted(uniques, xp.nan)
        uniques = uniques[: nan_idx + 1]
        if return_inverse:
            inverse[inverse > nan_idx] = nan_idx

        if return_counts:
            counts[nan_idx] = xp.sum(counts[nan_idx:])
            counts = counts[: nan_idx + 1]

    ret = (uniques,)

    if return_inverse:
        ret += (inverse,)

    if return_counts:
        ret += (counts,)

    return ret[0] if len(ret) == 1 else ret


class MissingValues(NamedTuple):
    """Data class for missing data information"""

    nan: bool
    none: bool
    # All the different nan instances seen, useful because they end up merged
    # into just np.nan:
    all_nans: Iterable[Real] = (np.nan,)

    def to_list(self):
        """Convert tuple to a list where None is always first."""
        output = []
        if self.none:
            output.append(None)
        if self.nan:
            output.append(np.nan)
        return output


def _extract_missing(values):
    """Extract missing values from `values`.

    Parameters
    ----------
    values: set
        Set of values to extract missing from.

    Returns
    -------
    output: set
        Set with missing values extracted.

    missing_values: MissingValues
        Object with missing value information.
    """
    missing_values_set = {
        value for value in values if value is None or is_scalar_nan(value)
    }

    if not missing_values_set:
        return values, MissingValues(nan=False, none=False)

    all_nans = missing_values_set - {None}
    if None in missing_values_set:
        if len(missing_values_set) == 1:
            output_missing_values = MissingValues(
                nan=False, none=True, all_nans=all_nans
            )
        else:
            # If there is more than one missing value, then it has to be
            # float('nan') or np.nan
            output_missing_values = MissingValues(
                nan=True, none=True, all_nans=all_nans
            )
    else:
        output_missing_values = MissingValues(nan=True, none=False, all_nans=all_nans)

    # create set without the missing values
    output = values - missing_values_set
    return output, output_missing_values


class _nandict(dict):
    """Dictionary with support for nans."""

    def __init__(self, mapping):
        super().__init__(mapping)
        for key, value in mapping.items():
            if is_scalar_nan(key):
                self.nan_value = value
                break

    def __missing__(self, key):
        if hasattr(self, "nan_value") and is_scalar_nan(key):
            return self.nan_value
        return -1


def _map_to_integer(values, uniques):
    """Map values based on their position in uniques.

    Values not present in `uniques` are encoded as -1.
    """
    xp, _ = get_namespace(values, uniques)
    table = _nandict({val: i for i, val in enumerate(uniques)})
    return xp.asarray([table[v] for v in values], device=device(values))


def _unique_python(values, *, return_inverse, return_counts):
    # Only used in `_unique`, see docstring there for details
    try:
        uniques_set = set(values)
        uniques_set, missing_values = _extract_missing(uniques_set)

        uniques = sorted(uniques_set)
        uniques.extend(missing_values.to_list())
        uniques = np.array(uniques, dtype=values.dtype)
    except TypeError:
        types = sorted(t.__qualname__ for t in set(type(v) for v in values))
        raise TypeError(
            "Encoders require their input argument must be uniformly "
            f"strings or numbers. Got {types}"
        )
    ret = (uniques,)

    if return_inverse:
        ret += (_map_to_integer(values, uniques),)

    if return_counts:
        ret += (_get_counts(values, uniques, missing_values.all_nans),)

    return ret[0] if len(ret) == 1 else ret


def _unique_categorical(
    values, *, return_inverse=False, return_counts=False, filter_present=True
):
    """Find uniques for categorical pandas Series without materializing values."""
    uniques = values.cat.categories.to_numpy()
    codes = values.cat.codes.to_numpy(copy=True)
    isna = codes == -1
    is_numeric = np.issubdtype(uniques.dtype, np.number)
    if is_numeric:
        is_sorted = (uniques[:-1] < uniques[1:]).all()
        # Pandas sorts numerical categories but user-specified categorical dtypes
        # can still be unordered. We reject unordered numerical categories because
        # downstream numerical encoding might rely on sorted categories.
        if not is_sorted:
            raise TypeError(
                "Unsorted categories are not supported for numerical categorical data"
            )

    codes[isna] = uniques.size
    if isna.any():
        if is_numeric:
            uniques = np.r_[uniques, np.nan]
        else:
            uniques = np.r_[uniques.astype(object, copy=False), np.nan]

    if filter_present or return_counts:
        counts = np.bincount(codes, minlength=uniques.size)
        is_present = counts > 0
    if filter_present and not is_present.all():
        codes_mapping = np.empty_like(is_present, dtype=np.intp)
        codes_mapping[is_present] = np.arange(is_present.sum())
        codes = codes_mapping[codes]
        uniques = uniques[is_present]
        counts = counts[is_present]

    ret = (uniques,)

    if return_inverse:
        ret += (codes,)

    if return_counts:
        ret += (counts,)

    return ret[0] if len(ret) == 1 else ret


def _encode_labels(values, *, uniques):
    """Encode labels into [0, n_uniques - 1].

    Relies on `_encode`, see docstring there for more details.

    Unknown values raise a ValueError, i.e. when `_encode` returns
    a non-empty diff.

    Parameters
    ----------
    values : ndarray
        Values to encode.
    uniques : ndarray
        The unique values in `values`. If the dtype is not object, then
        `uniques` needs to be sorted.

    Returns
    -------
    encoded : ndarray
        Encoded values.

    Raises
    ------
    ValueError
        If `values` contains labels that are not in `uniques`.
    """
    encoded, diff = _encode(values, uniques=uniques, return_diff=True)
    if size(diff):
        raise ValueError(f"y contains previously unseen labels: {diff}")
    return encoded


def _encode(values, *, uniques, return_diff=False):
    """Encode values into [0, n_uniques - 1].

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.
    The numpy method has the limitation that the `uniques` need to
    be sorted. Importantly, this is not checked but assumed to already be
    the case. The calling method needs to ensure this for all non-object
    values.

    Values that are not present in `uniques` are encoded as -1.

    Parameters
    ----------
    values : ndarray
        Values to encode.
    uniques : ndarray
        The unique values in `values`. If the dtype is not object, then
        `uniques` needs to be sorted.
    return_diff : bool, default=False
        If True, also return the unique values in `values` that are not
        present in `uniques`.

    Returns
    -------
    encoded : ndarray
        Encoded values.
    diff : ndarray
        The unique values present in `values` and not in `uniques`. Only
        returned if ``return_diff=True``.
    """
    xp, _ = get_namespace(values, uniques)
    if not xp.isdtype(values.dtype, "numeric"):
        encoded = _map_to_integer(values, uniques)
    else:
        encoded = xp.searchsorted(uniques, values)
        if size(uniques):
            max_idx = xp.asarray(
                uniques.shape[0] - 1, dtype=encoded.dtype, device=device(encoded)
            )
            encoded_safe = xp.minimum(encoded, max_idx)
            matches = (encoded < uniques.shape[0]) & (uniques[encoded_safe] == values)

            if xp.any(xp.isnan(uniques)):
                matches |= (
                    (encoded < uniques.shape[0])
                    & xp.isnan(uniques[encoded_safe])
                    & xp.isnan(values)
                )
        else:
            matches = xp.zeros_like(encoded, dtype=xp.bool)
        encoded[~matches] = -1

    if return_diff:
        diff = _unique(values[encoded == -1])
        return encoded, diff

    return encoded


def _get_counts(values, uniques, nan_values=(np.nan,)):
    """Get the count of each of the `uniques` in `values`.

    The counts will use the order passed in by `uniques`.  `np.nan` is assumed
    to be the last item in `uniques`, if it was one of the values.  For
    non-object dtypes, `uniques` is assumed to be sorted.
    """
    if values.dtype.kind in "OU":
        counter = Counter(values)
        output = np.zeros(len(uniques), dtype=np.int64)
        for i, item in enumerate(uniques):
            output[i] = counter[item]
        if len(uniques) > 0 and is_scalar_nan(uniques[-1]):
            # Should be the sum of all nans:
            output[-1] = sum(counter[nan] for nan in nan_values)
        return output

    unique_values, counts = _unique_np(values, return_counts=True)

    # Recorder unique_values based on input: `uniques`
    uniques_in_values = np.isin(uniques, unique_values, assume_unique=True)
    if np.isnan(unique_values[-1]) and np.isnan(uniques[-1]):
        uniques_in_values[-1] = True

    unique_valid_indices = np.searchsorted(unique_values, uniques[uniques_in_values])
    output = np.zeros_like(uniques, dtype=np.int64)
    output[uniques_in_values] = counts[unique_valid_indices]
    return output
