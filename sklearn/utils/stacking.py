import numpy as np
from scipy import sparse

def flexible_vstack(it, final_len=None):
    """Helper that concatenates the elements of an iterable along axis=0.

    Supports iterables of arrays, lists, sparse matrices or tuples thereof.

    Parameters
    ----------
    it :

    final_len :

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import sparse
    >>>
    >>> def make_example(typ):
    ...     yield typ([1, 2])
    ...     yield typ([3])
    ...     yield typ([4, 5, 6])
    ...
    >>> flexible_concatenate(make_example(list))
    [1, 2, 3, 4, 5, 6]
    >>> flexible_concatenate(make_example(np.array))
    array([1, 2, 3, 4, 5, 6])
    >>> flexible_concatenate(zip(make_example(list), make_example(np.array)))
    ([1, 2, 3, 4, 5, 6], array([1, 2, 3, 4, 5, 6]))
    >>> flexible_concatenate(make_example(np.array))
    array([1, 2, 3, 4, 5, 6])
    >>> flexible_concatenate(make_example(np.array), final_len=6)
    array([1, 2, 3, 4, 5, 6])
    >>> flexible_concatenate(make_example(
    ...     lambda x: np.array(x).reshape(-1, 1)))
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([[1], [2], [3], [4], [5], [6]])
    """

    def make_accumulator(prototype):
        if isinstance(prototype, tuple):
            return tuple(make_accumulator(y_proto) for y_proto in prototype)
        if isinstance(prototype, np.ndarray) and final_len is not None:
            return np.empty((final_len,) + prototype.shape[1:],
                            dtype=prototype.dtype)
        else:
            return []

    def accumulate(x, accumulator, prototype):
        if isinstance(prototype, tuple):
            for y, y_acc, y_prototype in zip(x, accumulator, prototype):
                n_rows = accumulate(y, y_acc, y_prototype)
                # XXX: could assert all n_rows are identical
            return n_rows
        elif isinstance(prototype, np.ndarray) and final_len is not None:
            accumulator[offset:offset + len(x)] = x
            return len(x)
        elif isinstance(prototype, list):
            accumulator.extend(x)
            return len(x)
        else:
            accumulator.append(x)
            if hasattr(x, 'shape'):
                return x.shape[0]
            return len(x)

    def finalize(accumulator, prototype):
        if isinstance(prototype, tuple):
            return tuple(finalize(y_acc, y_prototype)
                         for y_acc, y_prototype in zip(accumulator, prototype))
        elif isinstance(prototype, list):
            return accumulator
        elif isinstance(prototype, np.ndarray) and final_len is not None:
            return accumulator
        elif isinstance(prototype, np.ndarray):
            return np.concatenate(accumulator, axis=0)
        elif sparse.isspmatrix(prototype):
            return sparse.vstack(accumulator).asformat(prototype.format)
        else:
            raise NotImplementedError('No finalizing for accumulation of %s'
                                      % type(prototype))

    it = iter(it)
    try:
        first = next(it)
    except StopIteration:
        raise ValueError('Require at least one output from the iterator')

    accumulator = make_accumulator(first)
    offset = 0
    offset = accumulate(first, accumulator, first)
    for x in it:
        offset += accumulate(x, accumulator, first)

    if final_len is not None:
        assert offset == final_len, 'Expected %d, got %d' % (final_len, offset)

    return finalize(accumulator, first)
