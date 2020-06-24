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


def _unique_python(values, *, return_inverse):
    # Only used in `_uniques`, see docstring there for details
    try:
        uniques_set = set(values)
        missing_values_in_set = [value for value in (None, np.nan)
                                 if value in uniques_set]

        uniques_set -= set(missing_values_in_set)
        uniques = sorted(uniques_set)
        uniques.extend(missing_values_in_set)
        uniques = np.array(uniques, dtype=values.dtype)
    except TypeError:
        types = sorted(t.__qualname__
                       for t in set(type(v) for v in values))
        raise TypeError("Encoders require their input to be uniformly "
                        f"strings or numbers. Got {types}")

    if return_inverse:
        table = {val: i for i, val in enumerate(uniques)}
        inverse = np.array([table[v] for v in values])
        return uniques, inverse

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
    check_unknown : bool, default True
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
        table = {val: i for i, val in enumerate(uniques)}
        try:
            return np.array([table[v] for v in values])
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
    return_mask : bool, default False
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
        uniques_set = set(known_values)
        diff = set(values) - uniques_set
        if return_mask:
            if diff:
                valid_mask = np.array([val in uniques_set for val in values])
            else:
                valid_mask = np.ones(len(values), dtype=bool)

        missing_values_in_diff = [value for value in (None, np.nan)
                                  if value in diff]

        diff -= set(missing_values_in_diff)
        diff = list(diff)
        diff.extend(missing_values_in_diff)

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
    # valid_mask is
    return diff
