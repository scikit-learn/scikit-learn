import warnings
import functools
import sys
import numpy as np
import types
import numbers

import six

from ..util import img_as_float
from ._warnings import all_warnings, warn

__all__ = ['deprecated', 'get_bound_method_class', 'all_warnings',
           'safe_as_int', 'assert_nD', 'warn']


class skimage_deprecation(Warning):
    """Create our own deprecation class, since Python >= 2.7
    silences deprecations by default.

    """
    pass


class deprecated(object):
    """Decorator to mark deprecated functions with warning.

    Adapted from <http://wiki.python.org/moin/PythonDecoratorLibrary>.

    Parameters
    ----------
    alt_func : str
        If given, tell user what function to use instead.
    behavior : {'warn', 'raise'}
        Behavior during call to deprecated function: 'warn' = warn user that
        function is deprecated; 'raise' = raise error.
    removed_version : str
        The package version in which the deprecated function will be removed.
    """

    def __init__(self, alt_func=None, behavior='warn', removed_version=None):
        self.alt_func = alt_func
        self.behavior = behavior
        self.removed_version = removed_version

    def __call__(self, func):

        alt_msg = ''
        if self.alt_func is not None:
            alt_msg = ' Use ``%s`` instead.' % self.alt_func
        rmv_msg = ''
        if self.removed_version is not None:
            rmv_msg = (' and will be removed in version %s' %
                       self.removed_version)

        msg = ('Function ``%s`` is deprecated' % func.__name__ +
               rmv_msg + '.' + alt_msg)

        @functools.wraps(func)
        def wrapped(*args, **kwargs):
            if self.behavior == 'warn':
                func_code = six.get_function_code(func)
                warnings.simplefilter('always', skimage_deprecation)
                warnings.warn_explicit(msg,
                                       category=skimage_deprecation,
                                       filename=func_code.co_filename,
                                       lineno=func_code.co_firstlineno + 1)
            elif self.behavior == 'raise':
                raise skimage_deprecation(msg)
            return func(*args, **kwargs)

        # modify doc string to display deprecation warning
        doc = '**Deprecated function**.' + alt_msg
        if wrapped.__doc__ is None:
            wrapped.__doc__ = doc
        else:
            wrapped.__doc__ = doc + '\n\n    ' + wrapped.__doc__

        return wrapped


def get_bound_method_class(m):
    """Return the class for a bound method.

    """
    return m.im_class if sys.version < '3' else m.__self__.__class__


def safe_as_int(val, atol=1e-3):
    """
    Attempt to safely cast values to integer format.

    Parameters
    ----------
    val : scalar or iterable of scalars
        Number or container of numbers which are intended to be interpreted as
        integers, e.g., for indexing purposes, but which may not carry integer
        type.
    atol : float
        Absolute tolerance away from nearest integer to consider values in
        ``val`` functionally integers.

    Returns
    -------
    val_int : NumPy scalar or ndarray of dtype `np.int64`
        Returns the input value(s) coerced to dtype `np.int64` assuming all
        were within ``atol`` of the nearest integer.

    Notes
    -----
    This operation calculates ``val`` modulo 1, which returns the mantissa of
    all values. Then all mantissas greater than 0.5 are subtracted from one.
    Finally, the absolute tolerance from zero is calculated. If it is less
    than ``atol`` for all value(s) in ``val``, they are rounded and returned
    in an integer array. Or, if ``val`` was a scalar, a NumPy scalar type is
    returned.

    If any value(s) are outside the specified tolerance, an informative error
    is raised.

    Examples
    --------
    >>> safe_as_int(7.0)
    7

    >>> safe_as_int([9, 4, 2.9999999999])
    array([9, 4, 3])

    >>> safe_as_int(53.1)
    Traceback (most recent call last):
        ...
    ValueError: Integer argument required but received 53.1, check inputs.

    >>> safe_as_int(53.01, atol=0.01)
    53

    """
    mod = np.asarray(val) % 1                # Extract mantissa

    # Check for and subtract any mod values > 0.5 from 1
    if mod.ndim == 0:                        # Scalar input, cannot be indexed
        if mod > 0.5:
            mod = 1 - mod
    else:                                    # Iterable input, now ndarray
        mod[mod > 0.5] = 1 - mod[mod > 0.5]  # Test on each side of nearest int

    try:
        np.testing.assert_allclose(mod, 0, atol=atol)
    except AssertionError:
        raise ValueError("Integer argument required but received "
                         "{0}, check inputs.".format(val))

    return np.round(val).astype(np.int64)


def assert_nD(array, ndim, arg_name='image'):
    """
    Verify an array meets the desired ndims and array isn't empty.

    Parameters
    ----------
    array : array-like
        Input array to be validated
    ndim : int or iterable of ints
        Allowable ndim or ndims for the array.
    arg_name : str, optional
        The name of the array in the original function.

    """
    array = np.asanyarray(array)
    msg_incorrect_dim = "The parameter `%s` must be a %s-dimensional array"
    msg_empty_array = "The parameter `%s` cannot be an empty array"
    if isinstance(ndim, int):
        ndim = [ndim]
    if array.size == 0:
        raise ValueError(msg_empty_array % (arg_name))
    if not array.ndim in ndim:
        raise ValueError(msg_incorrect_dim % (arg_name, '-or-'.join([str(n) for n in ndim])))


def copy_func(f, name=None):
    """Create a copy of a function.

    Parameters
    ----------
    f : function
        Function to copy.
    name : str, optional
        Name of new function.

    """
    return types.FunctionType(six.get_function_code(f),
                              six.get_function_globals(f), name or f.__name__,
                              six.get_function_defaults(f), six.get_function_closure(f))


def check_random_state(seed):
    """Turn seed into a `np.random.RandomState` instance.

    Parameters
    ----------
    seed : None, int or np.random.RandomState
           If `seed` is None, return the RandomState singleton used by `np.random`.
           If `seed` is an int, return a new RandomState instance seeded with `seed`.
           If `seed` is already a RandomState instance, return it.

    Raises
    ------
    ValueError
        If `seed` is of the wrong type.

    """
    # Function originally from scikit-learn's module sklearn.utils.validation
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


def convert_to_float(image, preserve_range):
    """Convert input image to double image with the appropriate range.

    Parameters
    ----------
    image : ndarray
        Input image.
    preserve_range : bool
        Determines if the range of the image should be kept or transformed
        using img_as_float.

    Returns
    -------
    image : ndarray
        Transformed version of the input.
    """
    if preserve_range:
        image = image.astype(np.double)
    else:
        image = img_as_float(image)
    return image

