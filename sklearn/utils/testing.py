"""Testing utilities."""

# Copyright (c) 2011, 2012
# Authors: Pietro Berkes,
#          Andreas Muller
#          Mathieu Blondel
#          Olivier Grisel
#          Arnaud Joly
#          Denis Engemann
#          Giorgio Patrini
#          Thierry Guillemot
# License: BSD 3 clause
import os
import os.path as op
import inspect
import pkgutil
import warnings
import sys
import functools
import tempfile
from subprocess import check_output, STDOUT, CalledProcessError
from subprocess import TimeoutExpired


import scipy as sp
import scipy.io
from functools import wraps
from operator import itemgetter
from inspect import signature
from urllib.request import urlopen
from urllib.error import HTTPError

import tempfile
import shutil
import atexit
import unittest

# WindowsError only exist on Windows
try:
    WindowsError
except NameError:
    WindowsError = None

from numpy.testing import assert_allclose
from numpy.testing import assert_almost_equal
from numpy.testing import assert_approx_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_less
import numpy as np

import sklearn
from sklearn.base import (BaseEstimator, ClassifierMixin, ClusterMixin,
                          RegressorMixin, TransformerMixin)
from sklearn.utils import deprecated, IS_PYPY, _IS_32BIT
from sklearn.utils._joblib import joblib
from sklearn.utils._unittest_backport import TestCase

additional_names_in_all = []
try:
    from nose.tools import raises as _nose_raises
    deprecation_message = (
        'sklearn.utils.testing.raises has been deprecated in version 0.20 '
        'and will be removed in 0.22. Please use '
        'sklearn.utils.testing.assert_raises instead.')
    raises = deprecated(deprecation_message)(_nose_raises)
    additional_names_in_all.append('raises')
except ImportError:
    pass

try:
    from nose.tools import with_setup as _with_setup
    deprecation_message = (
        'sklearn.utils.testing.with_setup has been deprecated in version 0.20 '
        'and will be removed in 0.22.'
        'If your code relies on with_setup, please use'
        ' nose.tools.with_setup instead.')
    with_setup = deprecated(deprecation_message)(_with_setup)
    additional_names_in_all.append('with_setup')
except ImportError:
    pass

__all__ = ["assert_equal", "assert_not_equal", "assert_raises",
           "assert_raises_regexp", "assert_true",
           "assert_false", "assert_almost_equal", "assert_array_equal",
           "assert_array_almost_equal", "assert_array_less",
           "assert_less", "assert_less_equal",
           "assert_greater", "assert_greater_equal",
           "assert_approx_equal", "assert_allclose",
           "assert_run_python_script", "SkipTest"]
__all__.extend(additional_names_in_all)

_dummy = TestCase('__init__')
assert_equal = _dummy.assertEqual
assert_not_equal = _dummy.assertNotEqual
assert_raises = _dummy.assertRaises
SkipTest = unittest.case.SkipTest
assert_dict_equal = _dummy.assertDictEqual
assert_in = _dummy.assertIn
assert_not_in = _dummy.assertNotIn
assert_less = _dummy.assertLess
assert_greater = _dummy.assertGreater
assert_less_equal = _dummy.assertLessEqual
assert_greater_equal = _dummy.assertGreaterEqual

assert_raises_regex = _dummy.assertRaisesRegex
# assert_raises_regexp is deprecated in Python 3.4 in favor of
# assert_raises_regex but lets keep the backward compat in scikit-learn with
# the old name for now
assert_raises_regexp = assert_raises_regex

deprecation_message = "'assert_true' is deprecated in version 0.21 " \
                      "and will be removed in version 0.23. " \
                      "Please use 'assert' instead."
assert_true = deprecated(deprecation_message)(_dummy.assertTrue)

deprecation_message = "'assert_false' is deprecated in version 0.21 " \
                      "and will be removed in version 0.23. " \
                      "Please use 'assert' instead."
assert_false = deprecated(deprecation_message)(_dummy.assertFalse)


def assert_warns(warning_class, func, *args, **kw):
    """Test that a certain warning occurs.

    Parameters
    ----------
    warning_class : the warning class
        The class to test for, e.g. UserWarning.

    func : callable
        Callable object to trigger warnings.

    *args : the positional arguments to `func`.

    **kw : the keyword arguments to `func`

    Returns
    -------

    result : the return value of `func`

    """
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        result = func(*args, **kw)
        if hasattr(np, 'VisibleDeprecationWarning'):
            # Filter out numpy-specific warnings in numpy >= 1.9
            w = [e for e in w
                 if e.category is not np.VisibleDeprecationWarning]

        # Verify some things
        if not len(w) > 0:
            raise AssertionError("No warning raised when calling %s"
                                 % func.__name__)

        found = any(warning.category is warning_class for warning in w)
        if not found:
            raise AssertionError("%s did not give warning: %s( is %s)"
                                 % (func.__name__, warning_class, w))
    return result


def assert_warns_message(warning_class, message, func, *args, **kw):
    # very important to avoid uncontrolled state propagation
    """Test that a certain warning occurs and with a certain message.

    Parameters
    ----------
    warning_class : the warning class
        The class to test for, e.g. UserWarning.

    message : str | callable
        The message or a substring of the message to test for. If callable,
        it takes a string as the argument and will trigger an AssertionError
        if the callable returns `False`.

    func : callable
        Callable object to trigger warnings.

    *args : the positional arguments to `func`.

    **kw : the keyword arguments to `func`.

    Returns
    -------
    result : the return value of `func`

    """
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        if hasattr(np, 'VisibleDeprecationWarning'):
            # Let's not catch the numpy internal DeprecationWarnings
            warnings.simplefilter('ignore', np.VisibleDeprecationWarning)
        # Trigger a warning.
        result = func(*args, **kw)
        # Verify some things
        if not len(w) > 0:
            raise AssertionError("No warning raised when calling %s"
                                 % func.__name__)

        found = [issubclass(warning.category, warning_class) for warning in w]
        if not any(found):
            raise AssertionError("No warning raised for %s with class "
                                 "%s"
                                 % (func.__name__, warning_class))

        message_found = False
        # Checks the message of all warnings belong to warning_class
        for index in [i for i, x in enumerate(found) if x]:
            # substring will match, the entire message with typo won't
            msg = w[index].message  # For Python 3 compatibility
            msg = str(msg.args[0] if hasattr(msg, 'args') else msg)
            if callable(message):  # add support for certain tests
                check_in_message = message
            else:
                check_in_message = lambda msg: message in msg

            if check_in_message(msg):
                message_found = True
                break

        if not message_found:
            raise AssertionError("Did not receive the message you expected "
                                 "('%s') for <%s>, got: '%s'"
                                 % (message, func.__name__, msg))

    return result


def assert_warns_div0(func, *args, **kw):
    """Assume that numpy's warning for divide by zero is raised

    Handles the case of platforms that do not support warning on divide by zero

    Parameters
    ----------
    func
    *args
    **kw
    """

    with np.errstate(divide='warn', invalid='warn'):
        try:
            assert_warns(RuntimeWarning, np.divide, 1, np.zeros(1))
        except AssertionError:
            # This platform does not report numpy divide by zeros
            return func(*args, **kw)
        return assert_warns_message(RuntimeWarning,
                                    'invalid value encountered',
                                    func, *args, **kw)


# To remove when we support numpy 1.7
def assert_no_warnings(func, *args, **kw):
    """
    Parameters
    ----------
    func
    *args
    **kw
    """
    # very important to avoid uncontrolled state propagation
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')

        result = func(*args, **kw)
        if hasattr(np, 'VisibleDeprecationWarning'):
            # Filter out numpy-specific warnings in numpy >= 1.9
            w = [e for e in w
                 if e.category is not np.VisibleDeprecationWarning]

        if len(w) > 0:
            raise AssertionError("Got warnings when calling %s: [%s]"
                                 % (func.__name__,
                                    ', '.join(str(warning) for warning in w)))
    return result


def ignore_warnings(obj=None, category=Warning):
    """Context manager and decorator to ignore warnings.

    Note: Using this (in both variants) will clear all warnings
    from all python modules loaded. In case you need to test
    cross-module-warning-logging, this is not your tool of choice.

    Parameters
    ----------
    obj : callable or None
        callable where you want to ignore the warnings.
    category : warning class, defaults to Warning.
        The category to filter. If Warning, all categories will be muted.

    Examples
    --------
    >>> with ignore_warnings():
    ...     warnings.warn('buhuhuhu')

    >>> def nasty_warn():
    ...    warnings.warn('buhuhuhu')
    ...    print(42)

    >>> ignore_warnings(nasty_warn)()
    42
    """
    if isinstance(obj, type) and issubclass(obj, Warning):
        # Avoid common pitfall of passing category as the first positional
        # argument which result in the test not being run
        warning_name = obj.__name__
        raise ValueError(
            "'obj' should be a callable where you want to ignore warnings. "
            "You passed a warning class instead: 'obj={warning_name}'. "
            "If you want to pass a warning class to ignore_warnings, "
            "you should use 'category={warning_name}'".format(
                warning_name=warning_name))
    elif callable(obj):
        return _IgnoreWarnings(category=category)(obj)
    else:
        return _IgnoreWarnings(category=category)


class _IgnoreWarnings:
    """Improved and simplified Python warnings context manager and decorator.

    This class allows the user to ignore the warnings raised by a function.
    Copied from Python 2.7.5 and modified as required.

    Parameters
    ----------
    category : tuple of warning class, default to Warning
        The category to filter. By default, all the categories will be muted.

    """

    def __init__(self, category):
        self._record = True
        self._module = sys.modules['warnings']
        self._entered = False
        self.log = []
        self.category = category

    def __call__(self, fn):
        """Decorator to catch and hide warnings without visual nesting."""
        @wraps(fn)
        def wrapper(*args, **kwargs):
            clean_warning_registry()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", self.category)
                return fn(*args, **kwargs)

        return wrapper

    def __repr__(self):
        args = []
        if self._record:
            args.append("record=True")
        if self._module is not sys.modules['warnings']:
            args.append("module=%r" % self._module)
        name = type(self).__name__
        return "%s(%s)" % (name, ", ".join(args))

    def __enter__(self):
        if self._entered:
            raise RuntimeError("Cannot enter %r twice" % self)
        self._entered = True
        self._filters = self._module.filters
        self._module.filters = self._filters[:]
        self._showwarning = self._module.showwarning
        clean_warning_registry()
        warnings.simplefilter("ignore", self.category)

    def __exit__(self, *exc_info):
        if not self._entered:
            raise RuntimeError("Cannot exit %r without entering first" % self)
        self._module.filters = self._filters
        self._module.showwarning = self._showwarning
        self.log[:] = []
        clean_warning_registry()


def assert_raise_message(exceptions, message, function, *args, **kwargs):
    """Helper function to test the message raised in an exception.

    Given an exception, a callable to raise the exception, and
    a message string, tests that the correct exception is raised and
    that the message is a substring of the error thrown. Used to test
    that the specific message thrown during an exception is correct.

    Parameters
    ----------
    exceptions : exception or tuple of exception
        An Exception object.

    message : str
        The error message or a substring of the error message.

    function : callable
        Callable object to raise error.

    *args : the positional arguments to `function`.

    **kwargs : the keyword arguments to `function`.
    """
    try:
        function(*args, **kwargs)
    except exceptions as e:
        error_message = str(e)
        if message not in error_message:
            raise AssertionError("Error message does not include the expected"
                                 " string: %r. Observed error message: %r" %
                                 (message, error_message))
    else:
        # concatenate exception names
        if isinstance(exceptions, tuple):
            names = " or ".join(e.__name__ for e in exceptions)
        else:
            names = exceptions.__name__

        raise AssertionError("%s not raised by %s" %
                             (names, function.__name__))


def assert_allclose_dense_sparse(x, y, rtol=1e-07, atol=1e-9, err_msg=''):
    """Assert allclose for sparse and dense data.

    Both x and y need to be either sparse or dense, they
    can't be mixed.

    Parameters
    ----------
    x : array-like or sparse matrix
        First array to compare.

    y : array-like or sparse matrix
        Second array to compare.

    rtol : float, optional
        relative tolerance; see numpy.allclose

    atol : float, optional
        absolute tolerance; see numpy.allclose. Note that the default here is
        more tolerant than the default for numpy.testing.assert_allclose, where
        atol=0.

    err_msg : string, default=''
        Error message to raise.
    """
    if sp.sparse.issparse(x) and sp.sparse.issparse(y):
        x = x.tocsr()
        y = y.tocsr()
        x.sum_duplicates()
        y.sum_duplicates()
        assert_array_equal(x.indices, y.indices, err_msg=err_msg)
        assert_array_equal(x.indptr, y.indptr, err_msg=err_msg)
        assert_allclose(x.data, y.data, rtol=rtol, atol=atol, err_msg=err_msg)
    elif not sp.sparse.issparse(x) and not sp.sparse.issparse(y):
        # both dense
        assert_allclose(x, y, rtol=rtol, atol=atol, err_msg=err_msg)
    else:
        raise ValueError("Can only compare two sparse matrices,"
                         " not a sparse matrix and an array.")


@deprecated('deprecated in version 0.20 to be removed in version 0.22')
def fake_mldata(columns_dict, dataname, matfile, ordering=None):
    """Create a fake mldata data set.

    .. deprecated:: 0.20
        Will be removed in version 0.22

    Parameters
    ----------
    columns_dict : dict, keys=str, values=ndarray
        Contains data as columns_dict[column_name] = array of data.

    dataname : string
        Name of data set.

    matfile : string or file object
        The file name string or the file-like object of the output file.

    ordering : list, default None
        List of column_names, determines the ordering in the data set.

    Notes
    -----
    This function transposes all arrays, while fetch_mldata only transposes
    'data', keep that into account in the tests.
    """
    datasets = dict(columns_dict)

    # transpose all variables
    for name in datasets:
        datasets[name] = datasets[name].T

    if ordering is None:
        ordering = sorted(list(datasets.keys()))
    # NOTE: setting up this array is tricky, because of the way Matlab
    # re-packages 1D arrays
    datasets['mldata_descr_ordering'] = sp.empty((1, len(ordering)),
                                                 dtype='object')
    for i, name in enumerate(ordering):
        datasets['mldata_descr_ordering'][0, i] = name

    scipy.io.savemat(matfile, datasets, oned_as='column')


@deprecated('deprecated in version 0.20 to be removed in version 0.22')
class mock_mldata_urlopen:
    """Object that mocks the urlopen function to fake requests to mldata.

    When requesting a dataset with a name that is in mock_datasets, this object
    creates a fake dataset in a StringIO object and returns it. Otherwise, it
    raises an HTTPError.

    .. deprecated:: 0.20
        Will be removed in version 0.22

    Parameters
    ----------
    mock_datasets : dict
        A dictionary of {dataset_name: data_dict}, or
        {dataset_name: (data_dict, ordering). `data_dict` itself is a
        dictionary of {column_name: data_array}, and `ordering` is a list of
        column_names to determine the ordering in the data set (see
        :func:`fake_mldata` for details).
    """
    def __init__(self, mock_datasets):
        self.mock_datasets = mock_datasets

    def __call__(self, urlname):
        """
        Parameters
        ----------
        urlname : string
            The url
        """
        dataset_name = urlname.split('/')[-1]
        if dataset_name in self.mock_datasets:
            resource_name = '_' + dataset_name
            from io import BytesIO
            matfile = BytesIO()

            dataset = self.mock_datasets[dataset_name]
            ordering = None
            if isinstance(dataset, tuple):
                dataset, ordering = dataset
            fake_mldata(dataset, resource_name, matfile, ordering)

            matfile.seek(0)
            return matfile
        else:
            raise HTTPError(urlname, 404, dataset_name + " is not available",
                            [], None)


def install_mldata_mock(mock_datasets):
    """
    Parameters
    ----------
    mock_datasets : dict
        A dictionary of {dataset_name: data_dict}, or
        {dataset_name: (data_dict, ordering). `data_dict` itself is a
        dictionary of {column_name: data_array}, and `ordering` is a list of
        column_names to determine the ordering in the data set (see
        :func:`fake_mldata` for details).
    """
    # Lazy import to avoid mutually recursive imports
    from sklearn import datasets
    datasets.mldata.urlopen = mock_mldata_urlopen(mock_datasets)


def uninstall_mldata_mock():
    # Lazy import to avoid mutually recursive imports
    from sklearn import datasets
    datasets.mldata.urlopen = urlopen


def all_estimators(include_meta_estimators=None,
                   include_other=None, type_filter=None,
                   include_dont_test=None):
    """Get a list of all estimators from sklearn.

    This function crawls the module and gets all classes that inherit
    from BaseEstimator. Classes that are defined in test-modules are not
    included.
    By default meta_estimators such as GridSearchCV are also not included.

    Parameters
    ----------
    include_meta_estimators : boolean, default=False
        Deprecated, ignored.
        .. deprecated:: 0.21
           ``include_meta_estimators`` has been deprecated and has no effect in
           0.21 and will be removed in 0.23.

    include_other : boolean, default=False
        Deprecated, ignored.
        .. deprecated:: 0.21
           ``include_other`` has been deprecated and has not effect in 0.21 and
           will be removed in 0.23.

    type_filter : string, list of string,  or None, default=None
        Which kind of estimators should be returned. If None, no filter is
        applied and all estimators are returned.  Possible values are
        'classifier', 'regressor', 'cluster' and 'transformer' to get
        estimators only of these specific types, or a list of these to
        get the estimators that fit at least one of the types.

    include_dont_test : boolean, default=False
        Deprecated, ignored.
        .. deprecated:: 0.21
           ``include_dont_test`` has been deprecated and has no effect in 0.21
           and will be removed in 0.23.

    Returns
    -------
    estimators : list of tuples
        List of (name, class), where ``name`` is the class name as string
        and ``class`` is the actuall type of the class.
    """
    def is_abstract(c):
        if not(hasattr(c, '__abstractmethods__')):
            return False
        if not len(c.__abstractmethods__):
            return False
        return True

    if include_other is not None:
        warnings.warn("include_other was deprecated in version 0.21,"
                      " has no effect and will be removed in 0.23",
                      DeprecationWarning)

    if include_dont_test is not None:
        warnings.warn("include_dont_test was deprecated in version 0.21,"
                      " has no effect and will be removed in 0.23",
                      DeprecationWarning)

    if include_meta_estimators is not None:
        warnings.warn("include_meta_estimators was deprecated in version 0.21,"
                      " has no effect and will be removed in 0.23",
                      DeprecationWarning)

    all_classes = []
    # get parent folder
    path = sklearn.__path__
    for importer, modname, ispkg in pkgutil.walk_packages(
            path=path, prefix='sklearn.', onerror=lambda x: None):
        if ".tests." in modname or "externals" in modname:
            continue
        if IS_PYPY and ('_svmlight_format' in modname or
                        'feature_extraction._hashing' in modname):
            continue
        module = __import__(modname, fromlist="dummy")
        classes = inspect.getmembers(module, inspect.isclass)
        all_classes.extend(classes)

    all_classes = set(all_classes)

    estimators = [c for c in all_classes
                  if (issubclass(c[1], BaseEstimator) and
                      c[0] != 'BaseEstimator')]
    # get rid of abstract base classes
    estimators = [c for c in estimators if not is_abstract(c[1])]

    if type_filter is not None:
        if not isinstance(type_filter, list):
            type_filter = [type_filter]
        else:
            type_filter = list(type_filter)  # copy
        filtered_estimators = []
        filters = {'classifier': ClassifierMixin,
                   'regressor': RegressorMixin,
                   'transformer': TransformerMixin,
                   'cluster': ClusterMixin}
        for name, mixin in filters.items():
            if name in type_filter:
                type_filter.remove(name)
                filtered_estimators.extend([est for est in estimators
                                            if issubclass(est[1], mixin)])
        estimators = filtered_estimators
        if type_filter:
            raise ValueError("Parameter type_filter must be 'classifier', "
                             "'regressor', 'transformer', 'cluster' or "
                             "None, got"
                             " %s." % repr(type_filter))

    # drop duplicates, sort for reproducibility
    # itemgetter is used to ensure the sort does not extend to the 2nd item of
    # the tuple
    return sorted(set(estimators), key=itemgetter(0))


def set_random_state(estimator, random_state=0):
    """Set random state of an estimator if it has the `random_state` param.

    Parameters
    ----------
    estimator : object
        The estimator
    random_state : int, RandomState instance or None, optional, default=0
        Pseudo random number generator state.  If int, random_state is the seed
        used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.
    """
    if "random_state" in estimator.get_params():
        estimator.set_params(random_state=random_state)


try:
    import pytest

    skip_if_32bit = pytest.mark.skipif(_IS_32BIT,
                                       reason='skipped on 32bit platforms')
    skip_travis = pytest.mark.skipif(os.environ.get('TRAVIS') == 'true',
                                     reason='skip on travis')
    fails_if_pypy = pytest.mark.xfail(IS_PYPY, raises=NotImplementedError,
                                      reason='not compatible with PyPy')
    skip_if_no_parallel = pytest.mark.skipif(not joblib.parallel.mp,
                                             reason="joblib is in serial mode")

    #  Decorator for tests involving both BLAS calls and multiprocessing.
    #
    #  Under POSIX (e.g. Linux or OSX), using multiprocessing in conjunction
    #  with some implementation of BLAS (or other libraries that manage an
    #  internal posix thread pool) can cause a crash or a freeze of the Python
    #  process.
    #
    #  In practice all known packaged distributions (from Linux distros or
    #  Anaconda) of BLAS under Linux seems to be safe. So we this problem seems
    #  to only impact OSX users.
    #
    #  This wrapper makes it possible to skip tests that can possibly cause
    #  this crash under OS X with.
    #
    #  Under Python 3.4+ it is possible to use the `forkserver` start method
    #  for multiprocessing to avoid this issue. However it can cause pickling
    #  errors on interactively defined functions. It therefore not enabled by
    #  default.

    if_safe_multiprocessing_with_blas = pytest.mark.skipif(
            sys.platform == 'darwin',
            reason="Possible multi-process bug with some BLAS")
except ImportError:
    pass


def clean_warning_registry():
    """Clean Python warning registry for easier testing of warning messages.

    We may not need to do this any more when getting rid of Python 2, not
    entirely sure. See https://bugs.python.org/issue4180 and
    https://bugs.python.org/issue21724 for more details.

    """
    reg = "__warningregistry__"
    for mod_name, mod in list(sys.modules.items()):
        if hasattr(mod, reg):
            getattr(mod, reg).clear()


def check_skip_network():
    if int(os.environ.get('SKLEARN_SKIP_NETWORK_TESTS', 0)):
        raise SkipTest("Text tutorial requires large dataset download")


def _delete_folder(folder_path, warn=False):
    """Utility function to cleanup a temporary folder if still existing.

    Copy from joblib.pool (for independence).
    """
    try:
        if os.path.exists(folder_path):
            # This can fail under windows,
            #  but will succeed when called by atexit
            shutil.rmtree(folder_path)
    except WindowsError:
        if warn:
            warnings.warn("Could not delete temporary folder %s" % folder_path)


class TempMemmap:
    """
    Parameters
    ----------
    data
    mmap_mode
    """
    def __init__(self, data, mmap_mode='r'):
        self.mmap_mode = mmap_mode
        self.data = data

    def __enter__(self):
        data_read_only, self.temp_folder = create_memmap_backed_data(
            self.data, mmap_mode=self.mmap_mode, return_folder=True)
        return data_read_only

    def __exit__(self, exc_type, exc_val, exc_tb):
        _delete_folder(self.temp_folder)


def create_memmap_backed_data(data, mmap_mode='r', return_folder=False):
    """
    Parameters
    ----------
    data
    mmap_mode
    return_folder
    """
    temp_folder = tempfile.mkdtemp(prefix='sklearn_testing_')
    atexit.register(functools.partial(_delete_folder, temp_folder, warn=True))
    filename = op.join(temp_folder, 'data.pkl')
    joblib.dump(data, filename)
    memmap_backed_data = joblib.load(filename, mmap_mode=mmap_mode)
    result = (memmap_backed_data if not return_folder
              else (memmap_backed_data, temp_folder))
    return result


# Utils to test docstrings


def _get_args(function, varargs=False):
    """Helper to get function arguments"""

    try:
        params = signature(function).parameters
    except ValueError:
        # Error on builtin C function
        return []
    args = [key for key, param in params.items()
            if param.kind not in (param.VAR_POSITIONAL, param.VAR_KEYWORD)]
    if varargs:
        varargs = [param.name for param in params.values()
                   if param.kind == param.VAR_POSITIONAL]
        if len(varargs) == 0:
            varargs = None
        return args, varargs
    else:
        return args


def _get_func_name(func, class_name=None):
    """Get function full name

    Parameters
    ----------
    func : callable
        The function object.
    class_name : string, optional (default: None)
       If ``func`` is a class method and the class name is known specify
       class_name for the error message.

    Returns
    -------
    name : str
        The function name.
    """
    parts = []
    module = inspect.getmodule(func)
    if module:
        parts.append(module.__name__)
    if class_name is not None:
        parts.append(class_name)
    elif hasattr(func, 'im_class'):
        parts.append(func.im_class.__name__)

    parts.append(func.__name__)
    return '.'.join(parts)


def check_docstring_parameters(func, doc=None, ignore=None, class_name=None):
    """Helper to check docstring

    Parameters
    ----------
    func : callable
        The function object to test.
    doc : str, optional (default: None)
        Docstring if it is passed manually to the test.
    ignore : None | list
        Parameters to ignore.
    class_name : string, optional (default: None)
       If ``func`` is a class method and the class name is known specify
       class_name for the error message.

    Returns
    -------
    incorrect : list
        A list of string describing the incorrect results.
    """
    from numpydoc import docscrape
    incorrect = []
    ignore = [] if ignore is None else ignore

    func_name = _get_func_name(func, class_name=class_name)
    if (not func_name.startswith('sklearn.') or
            func_name.startswith('sklearn.externals')):
        return incorrect
    # Don't check docstring for property-functions
    if inspect.isdatadescriptor(func):
        return incorrect
    # Don't check docstring for setup / teardown pytest functions
    if func_name.split('.')[-1] in ('setup_module', 'teardown_module'):
        return incorrect
    # Dont check estimator_checks module
    if func_name.split('.')[2] == 'estimator_checks':
        return incorrect
    args = list(filter(lambda x: x not in ignore, _get_args(func)))
    # drop self
    if len(args) > 0 and args[0] == 'self':
        args.remove('self')

    if doc is None:
        with warnings.catch_warnings(record=True) as w:
            try:
                doc = docscrape.FunctionDoc(func)
            except Exception as exp:
                incorrect += [func_name + ' parsing error: ' + str(exp)]
                return incorrect
        if len(w):
            raise RuntimeError('Error for %s:\n%s' % (func_name, w[0]))

    param_names = []
    for name, type_definition, param_doc in doc['Parameters']:
        if not type_definition.strip():
            if ':' in name and name[:name.index(':')][-1:].strip():
                incorrect += [func_name +
                              ' There was no space between the param name and '
                              'colon (%r)' % name]
            elif name.rstrip().endswith(':'):
                incorrect += [func_name +
                              ' Parameter %r has an empty type spec. '
                              'Remove the colon' % (name.lstrip())]

        if '*' not in name:
            param_names.append(name.split(':')[0].strip('` '))

    param_names = list(filter(lambda x: x not in ignore, param_names))

    if len(param_names) != len(args):
        bad = str(sorted(list(set(param_names) ^ set(args))))
        incorrect += [func_name + ' arg mismatch: ' + bad]
    else:
        for n1, n2 in zip(param_names, args):
            if n1 != n2:
                incorrect += [func_name + ' ' + n1 + ' != ' + n2]
    return incorrect


def assert_run_python_script(source_code, timeout=60):
    """Utility to check assertions in an independent Python subprocess.

    The script provided in the source code should return 0 and not print
    anything on stderr or stdout.

    This is a port from cloudpickle https://github.com/cloudpipe/cloudpickle

    Parameters
    ----------
    source_code : str
        The Python source code to execute.
    timeout : int
        Time in seconds before timeout.
    """
    fd, source_file = tempfile.mkstemp(suffix='_src_test_sklearn.py')
    os.close(fd)
    try:
        with open(source_file, 'wb') as f:
            f.write(source_code.encode('utf-8'))
        cmd = [sys.executable, source_file]
        cwd = op.normpath(op.join(op.dirname(sklearn.__file__), '..'))
        env = os.environ.copy()
        try:
            env["PYTHONPATH"] = os.pathsep.join([cwd, env["PYTHONPATH"]])
        except KeyError:
            env["PYTHONPATH"] = cwd
        kwargs = {
            'cwd': cwd,
            'stderr': STDOUT,
            'env': env
        }
        # If coverage is running, pass the config file to the subprocess
        coverage_rc = os.environ.get("COVERAGE_PROCESS_START")
        if coverage_rc:
            kwargs['env']['COVERAGE_PROCESS_START'] = coverage_rc

        kwargs['timeout'] = timeout
        try:
            try:
                out = check_output(cmd, **kwargs)
            except CalledProcessError as e:
                raise RuntimeError(u"script errored with output:\n%s"
                                   % e.output.decode('utf-8'))
            if out != b"":
                raise AssertionError(out.decode('utf-8'))
        except TimeoutExpired as e:
            raise RuntimeError(u"script timeout, output so far:\n%s"
                               % e.output.decode('utf-8'))
    finally:
        os.unlink(source_file)
