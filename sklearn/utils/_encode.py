import numpy as np
from . import is_scalar_nan


def _unique(values, *, return_inverse=False):
    """Helper function to find unique values with support for python objects.

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.

    Parameters
    ----------
    values : ndarray
        Values to check for unknowns.

    return_inverse : bool, default=False
        If True, also return the indices of the unique values.

    Returns
    -------
    unique : ndarray
        The sorted unique values.

    unique_inverse : ndarray
        The indices to reconstruct the original array from the unique array.
        Only provided if `return_inverse` is True.
    """
    if values.dtype == object:
        return _unique_python(values, return_inverse=return_inverse)
    # numerical
    out = np.unique(values, return_inverse=return_inverse)

    if return_inverse:
        uniques, inverse = out
    else:
        uniques, inverse = out, None

    # np.unique will have duplicate missing values at the end of `uniques`
    # here we clip the nans and remove it from uniques
    if uniques.size and is_scalar_nan(uniques[-1]):
        nan_idx = np.searchsorted(uniques, np.nan)
        uniques = uniques[:nan_idx + 1]
        if return_inverse:
            inverse[inverse > nan_idx] = nan_idx

    if return_inverse:
        return uniques, inverse
    return uniques


def _split_missing(items, none_is_missing=True):
    """Split items into a set and a list of missing values. If none_is_missing
    is True, then None and nan are considered missing. If none_is_missing
    is False, then only nans are considered missing.
    """
    def is_missing(value):
        if none_is_missing:
            return value is None or is_scalar_nan(value)
        return is_scalar_nan(value)

    # Return items without missing items
    missing_values = [value for value in items if is_missing(value)]

    if not missing_values:
        return items, []

    # Enforces an order where None always comes first
    if None in missing_values:
        if len(missing_values) == 1:
            output_missing_values = [None]
        else:
            output_missing_values = [None, np.nan]
    else:
        output_missing_values = [np.nan]

    # create set without the missing values
    output = set(value for value in items if not is_missing(value))
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


def _unique_python(values, *, return_inverse):
    # Only used in `_uniques`, see docstring there for details
    try:
        uniques_set = set(values)
        uniques_set, missing_values = _split_missing(uniques_set)

        uniques = sorted(uniques_set)
        uniques.extend(missing_values)
        uniques = np.array(uniques, dtype=values.dtype)
    except TypeError:
        types = sorted(t.__qualname__
                       for t in set(type(v) for v in values))
        raise TypeError("Encoders require their input to be uniformly "
                        f"strings or numbers. Got {types}")

    if return_inverse:
        return uniques, _map_to_integer(values, uniques)

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
        values_set, nan_in_values = _split_missing(values_set,
                                                   none_is_missing=False)

        uniques_set = set(known_values)
        uniques_set, nan_in_uniques = _split_missing(uniques_set,
                                                     none_is_missing=False)
        diff = values_set - uniques_set
        is_missing_in_uniques = bool(nan_in_uniques)

        nan_in_diff = nan_in_values and not nan_in_uniques

        def is_valid(value):
            value_in_set = value in uniques_set
            if nan_in_uniques and is_scalar_nan(value):
                return is_missing_in_uniques
            return value_in_set

        if return_mask:
            if diff or nan_in_diff:
                valid_mask = np.array([is_valid(value) for value in values])
            else:
                valid_mask = np.ones(len(values), dtype=bool)

        # ensure that None is at the end
        none_in_diff = None in diff
        if none_in_diff:
            diff.remove(None)
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
