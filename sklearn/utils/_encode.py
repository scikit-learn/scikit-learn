from collections import Counter
import numpy as np


def _unique(values, *, return_inverse=False, return_counts=False):
    """Helper function to find unique values with support for python objects.

    Uses pure python method for object dtype, and numpy method for
    all other dtypes.

    Parameters
    ----------
    values : ndarray
        Values to check for unknowns.

    return_inverse : bool, default=False
        If True, also return the indices of the unique values.

    return_count : bool, default=False
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
        return _unique_python(values, return_inverse=return_inverse,
                              return_counts=return_counts)
    # numerical
    return np.unique(values, return_inverse=return_inverse,
                     return_counts=return_counts)


def _unique_python(values, *, return_inverse, return_counts):
    # Only used in `_uniques`, see docstring there for details
    try:
        uniques = sorted(set(values))
        uniques = np.array(uniques, dtype=values.dtype)
    except TypeError:
        types = sorted(t.__qualname__
                       for t in set(type(v) for v in values))
        raise TypeError("Encoders require their input to be uniformly "
                        f"strings or numbers. Got {types}")
    ret = (uniques, )

    if return_inverse:
        table = {val: i for i, val in enumerate(uniques)}
        inverse = np.array([table[v] for v in values])
        ret += (inverse, )

    if return_counts:
        counts = _get_counts(values, uniques)
        ret += (counts, )

    if len(ret) == 1:
        ret = ret[0]

    return ret


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
    if values.dtype == object:
        uniques_set = set(known_values)
        diff = list(set(values) - uniques_set)
        if return_mask:
            if diff:
                valid_mask = np.array([val in uniques_set for val in values])
            else:
                valid_mask = np.ones(len(values), dtype=bool)
            return diff, valid_mask
        else:
            return diff
    else:
        unique_values = np.unique(values)
        diff = list(np.setdiff1d(unique_values, known_values,
                                 assume_unique=True))
        if return_mask:
            if diff:
                valid_mask = np.in1d(values, known_values)
            else:
                valid_mask = np.ones(len(values), dtype=bool)
            return diff, valid_mask
        else:
            return diff


def _get_counts(values, uniques):
    """Get the count of each of the `uniques` in `values`. The counts will use
    the order passed in by `uniques`.

    For non-object dtypes, `uniques` is assumed to be sorted.
    """
    if values.dtype == object:
        counter = Counter(values)
        counts = np.array([counter[item] for item in uniques],
                          dtype=int)
        return counts

    unique_values, counts = np.unique(values, return_counts=True)
    uniques_in_values = np.isin(uniques, unique_values, assume_unique=True)
    unique_valid_indices = np.searchsorted(unique_values,
                                           uniques[uniques_in_values])

    output = np.zeros_like(uniques)
    output[uniques_in_values] = counts[unique_valid_indices]
    return output
