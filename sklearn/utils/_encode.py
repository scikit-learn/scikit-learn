from typing import NamedTuple

import numpy as np
from . import is_scalar_nan
from collections import Counter


def _unique(values, *, return_inverse=False, return_counts=False, cats=None):
    """Helper function to find unique values with support for python objects.

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.

    Parameters
    ----------
    values : ndarray
        Values to check for unknowns.

    return_inverse : bool, default=False
        If True, also return the indices of the unique values.

    return_counts: bool, default=False
        If True, also return count of the unique values.

    cats: list of array-like, default=False
        List of specified categories. Will add zero count categories.
    Returns
    -------
    unique : ndarray
        The sorted unique values or cats if `cats` specified.

    unique_inverse : ndarray
        The indices to reconstruct the original array from the unique array.
        Only provided if `return_inverse` is True.

    counts : ndarray
        Count of each value. Indicies match unique.
        Only provided if `return_counts` is True.
    """

    if cats is not None:
        cat_zero_count = []
        for cat in cats:
            if cat not in values:
                # add values to be counted to be set to 0 later
                values = np.append(values, cat)
                cat_zero_count.append(cat)

    if values.dtype == object and cats is not None:
        uniques, counts = _unique_python(values,
                                         return_inverse=return_inverse,
                                         return_counts=return_counts)
        uniques, counts = _zero_counts(uniques, counts, cat_zero_count)
        counts = _reindex_specified_cats(uniques, cats, counts)
        return cats, counts
    elif values.dtype == object:
        return _unique_python(values,
                              return_inverse=return_inverse,
                              return_counts=return_counts)
    # numerical
    out = np.unique(values,
                    return_inverse=return_inverse,
                    return_counts=return_counts)

    if return_inverse:
        uniques, inverse = out
    elif return_counts:
        uniques, counts = out
    else:
        uniques = out

    # np.unique will have duplicate missing values at the end of `uniques`
    # here we clip the nans and remove it from uniques
    if uniques.size and is_scalar_nan(uniques[-1]):
        nan_idx = np.searchsorted(uniques, np.nan)
        uniques = uniques[:nan_idx + 1]
        if return_inverse:
            inverse[inverse > nan_idx] = nan_idx
        if return_counts:
            counts[counts > nan_idx] = nan_idx

    if return_inverse:
        return uniques, inverse
    elif return_counts and cats is not None:
        uniques, counts = out
        uniques, counts = _zero_counts(uniques, counts, cat_zero_count)
        counts = _reindex_specified_cats(uniques, cats, counts)
        return cats, counts
    elif return_counts:
        return uniques, counts
    return uniques


def _reindex_specified_cats(uniques, cats, counts):
    """Reindex counts by specified categories"""
    nan_idx_sorted = [idx for idx, cat in enumerate(uniques)
                      if is_scalar_nan(cat) or cat is None]
    nan_idx_cats = [idx for idx, cat in enumerate(cats)
                    if is_scalar_nan(cat) or cat is None]

    if len(nan_idx_cats) != 0:
        nan_cats = cats[nan_idx_cats]
        nan_counts = counts[nan_idx_sorted]

        cats = np.delete(cats, nan_idx_cats)
        counts = np.delete(counts, nan_idx_sorted)

    idx = np.argsort(cats)
    idx_2 = np.argsort(idx)
    counts = counts[idx_2]

    if len(nan_idx_cats) != 0:
        cats = np.append(cats, nan_cats)
        counts = np.insert(counts, nan_idx_cats, nan_counts)

    return counts


def _zero_counts(uniques, counts, cat_zero_count):
    """Set added specified categories count to 0"""
    for cat in cat_zero_count:
        idx = np.where(uniques == cat)
        if len(idx) != 0:
            counts[idx[0]] = 0
    return uniques, counts


class MissingValues(NamedTuple):
    """Data class for missing data information"""
    nan: bool
    none: bool

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
    missing_values_set = {value for value in values
                          if value is None or is_scalar_nan(value)}

    if not missing_values_set:
        return values, MissingValues(nan=False, none=False)

    if None in missing_values_set:
        if len(missing_values_set) == 1:
            output_missing_values = MissingValues(nan=False, none=True)
        else:
            # If there is more than one missing value, then it has to be
            # float('nan') or np.nan
            output_missing_values = MissingValues(nan=True, none=True)
    else:
        output_missing_values = MissingValues(nan=True, none=False)

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
        if hasattr(self, 'nan_value') and is_scalar_nan(key):
            return self.nan_value
        raise KeyError(key)


def _map_to_integer(values, uniques):
    """Map values based on its position in uniques."""
    table = _nandict({val: i for i, val in enumerate(uniques)})
    return np.array([table[v] for v in values])


def _unique_python(values, *, return_inverse, return_counts):
    # Only used in `_uniques`, see docstring there for details
    try:
        uniques_set = set(values)
        uniques_set, missing_values = _extract_missing(uniques_set)

        uniques = sorted(uniques_set)
        uniques.extend(missing_values.to_list())
        uniques = np.array(uniques, dtype=values.dtype)
    except TypeError:
        types = sorted(t.__qualname__
                       for t in set(type(v) for v in values))
        raise TypeError("Encoders require their input to be uniformly "
                        f"strings or numbers. Got {types}")
    if return_inverse and not return_counts:
        return uniques, _map_to_integer(values, uniques)
    if return_counts:
        ctr_all = Counter(values)
        ctr = {k: ctr_all[k] for k in uniques_set}
        ctr = dict(sorted(ctr.items()))
        if len(missing_values.to_list()) != 0:
            ctr_missing = {k: ctr_all[k] for k in missing_values.to_list()}
            ctr.update(ctr_missing)
        uniques = np.array(list(ctr.keys()), dtype=values.dtype)
        counts = np.array(list(ctr.values()), dtype=int)
        return uniques, counts
    return uniques


def _encode(values, *, uniques, check_unknown=True):
    """Helper function to encode values into [0, n_uniques - 1].

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.
    The numpy method has the limitation that the `uniques` need to
    be sorted. Importantly, this is not checked but assumed to already be
    the case. The calling method needs to ensure this for all non-object
    values.

    Parameters
    ----------
    values : ndarray
        Values to encode.
    uniques : ndarray
        The unique values in `values`. If the dtype is not object, then
        `uniques` needs to be sorted.
    check_unknown : bool, default=True
        If True, check for values in `values` that are not in `unique`
        and raise an error. This is ignored for object dtype, and treated as
        True in this case. This parameter is useful for
        _BaseEncoder._transform() to avoid calling _check_unknown()
        twice.

    Returns
    -------
    encoded : ndarray
        Encoded values
    """
    if values.dtype == object:
        try:
            return _map_to_integer(values, uniques)
        except KeyError as e:
            raise ValueError(f"y contains previously unseen labels: {str(e)}")
    else:
        if check_unknown:
            diff = _check_unknown(values, uniques)
            if diff:
                raise ValueError(f"y contains previously unseen labels: "
                                 f"{str(diff)}")
        return np.searchsorted(uniques, values)


def _check_unknown(values, known_values, return_mask=False):
    """
    Helper function to check for unknowns in values to be encoded.

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.

    Parameters
    ----------
    values : array
        Values to check for unknowns.
    known_values : array
        Known values. Must be unique.
    return_mask : bool, default=False
        If True, return a mask of the same shape as `values` indicating
        the valid values.

    Returns
    -------
    diff : list
        The unique values present in `values` and not in `know_values`.
    valid_mask : boolean array
        Additionally returned if ``return_mask=True``.

    """
    valid_mask = None

    if values.dtype.kind in 'UO':
        values_set = set(values)
        values_set, missing_in_values = _extract_missing(values_set)

        uniques_set = set(known_values)
        uniques_set, missing_in_uniques = _extract_missing(uniques_set)
        diff = values_set - uniques_set

        nan_in_diff = missing_in_values.nan and not missing_in_uniques.nan
        none_in_diff = missing_in_values.none and not missing_in_uniques.none

        def is_valid(value):
            return (value in uniques_set or
                    missing_in_uniques.none and value is None or
                    missing_in_uniques.nan and is_scalar_nan(value))

        if return_mask:
            if diff or nan_in_diff or none_in_diff:
                valid_mask = np.array([is_valid(value) for value in values])
            else:
                valid_mask = np.ones(len(values), dtype=bool)

        diff = list(diff)
        if none_in_diff:
            diff.append(None)
        if nan_in_diff:
            diff.append(np.nan)
    else:
        unique_values = np.unique(values)
        diff = np.setdiff1d(unique_values, known_values,
                            assume_unique=True)
        if return_mask:
            if diff.size:
                valid_mask = np.in1d(values, known_values)
            else:
                valid_mask = np.ones(len(values), dtype=bool)

        # check for nans in the known_values
        if np.isnan(known_values).any():
            diff_is_nan = np.isnan(diff)
            if diff_is_nan.any():
                # removes nan from valid_mask
                if diff.size and return_mask:
                    is_nan = np.isnan(values)
                    valid_mask[is_nan] = 1

                # remove nan from diff
                diff = diff[~diff_is_nan]
        diff = list(diff)

    if return_mask:
        return diff, valid_mask
    return diff
