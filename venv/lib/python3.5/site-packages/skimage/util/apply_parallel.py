from math import ceil
from multiprocessing import cpu_count

__all__ = ['apply_parallel']


try:
    import dask.array as da
    dask_available = True
except ImportError:
    dask_available = False


def _get_chunks(shape, ncpu):
    """Split the array into equal sized chunks based on the number of
    available processors. The last chunk in each dimension absorbs the
    remainder array elements if the number of CPUs does not divide evenly into
    the number of array elements.

    Examples
    --------
    >>> _get_chunks((4, 4), 4)
    ((2, 2), (2, 2))
    >>> _get_chunks((4, 4), 2)
    ((2, 2), (4,))
    >>> _get_chunks((5, 5), 2)
    ((2, 3), (5,))
    >>> _get_chunks((2, 4), 2)
    ((1, 1), (4,))
    """
    chunks = []
    nchunks_per_dim = int(ceil(ncpu ** (1./len(shape))))

    used_chunks = 1
    for i in shape:
        if used_chunks < ncpu:
            regular_chunk = i // nchunks_per_dim
            remainder_chunk = regular_chunk + (i % nchunks_per_dim)

            if regular_chunk == 0:
                chunk_lens = (remainder_chunk,)
            else:
                chunk_lens = ((regular_chunk,) * (nchunks_per_dim - 1) +
                              (remainder_chunk,))
        else:
            chunk_lens = (i,)

        chunks.append(chunk_lens)
        used_chunks *= nchunks_per_dim
    return tuple(chunks)


def apply_parallel(function, array, chunks=None, depth=0, mode=None,
                 extra_arguments=(), extra_keywords={}):
    """Map a function in parallel across an array.

    Split an array into possibly overlapping chunks of a given depth and
    boundary type, call the given function in parallel on the chunks, combine
    the chunks and return the resulting array.

    Parameters
    ----------
    function : function
        Function to be mapped which takes an array as an argument.
    array : numpy array
        Array which the function will be applied to.
    chunks : int, tuple, or tuple of tuples, optional
        A single integer is interpreted as the length of one side of a square
        chunk that should be tiled across the array.  One tuple of length
        ``array.ndim`` represents the shape of a chunk, and it is tiled across
        the array.  A list of tuples of length ``ndim``, where each sub-tuple
        is a sequence of chunk sizes along the corresponding dimension. If
        None, the array is broken up into chunks based on the number of
        available cpus. More information about chunks is in the documentation
        `here <https://dask.pydata.org/en/latest/array-design.html>`_.
    depth : int, optional
        Integer equal to the depth of the added boundary cells. Defaults to
        zero.
    mode : {'reflect', 'symmetric', 'periodic', 'wrap', 'nearest', 'edge'}, optional
        type of external boundary padding.
    extra_arguments : tuple, optional
        Tuple of arguments to be passed to the function.
    extra_keywords : dictionary, optional
        Dictionary of keyword arguments to be passed to the function.

    Notes
    -----
    Numpy edge modes 'symmetric', 'wrap', and 'edge' are converted to the
    equivalent `dask` boundary modes 'reflect', 'periodic' and 'nearest',
    respectively.
    """
    if not dask_available:
        raise RuntimeError("Could not import 'dask'.  Please install "
                           "using 'pip install dask'")

    if chunks is None:
        shape = array.shape
        try:
            ncpu = cpu_count()
        except NotImplementedError:
            ncpu = 4
        chunks = _get_chunks(shape, ncpu)

    if mode == 'wrap':
        mode = 'periodic'
    elif mode == 'symmetric':
        mode = 'reflect'
    elif mode == 'edge':
        mode = 'nearest'

    def wrapped_func(arr):
        return function(arr, *extra_arguments, **extra_keywords)

    darr = da.from_array(array, chunks=chunks)
    return darr.map_overlap(wrapped_func, depth, boundary=mode).compute()
