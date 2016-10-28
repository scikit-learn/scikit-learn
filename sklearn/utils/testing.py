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
import inspect
import pkgutil
import warnings
import sys
import re
import platform
import struct

import scipy as sp
import scipy.io
from functools import wraps
from operator import itemgetter
try:
    # Python 2
    from urllib2 import urlopen
    from urllib2 import HTTPError
except ImportError:
    # Python 3+
    from urllib.request import urlopen
    from urllib.error import HTTPError

import tempfile
import shutil
import os.path as op
import atexit
import unittest

# WindowsError only exist on Windows
try:
    WindowsError
except NameError:
    WindowsError = None

import sklearn
from sklearn.base import BaseEstimator
from sklearn.externals import joblib

from nose.tools import raises
from nose import with_setup

from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_less
from numpy.testing import assert_approx_equal
import numpy as np

from sklearn.base import (ClassifierMixin, RegressorMixin, TransformerMixin,
                          ClusterMixin)
from sklearn.cluster import DBSCAN

__all__ = ["assert_equal", "assert_not_equal", "assert_raises",
           "assert_raises_regexp", "raises", "with_setup", "assert_true",
           "assert_false", "assert_almost_equal", "assert_array_equal",
           "assert_array_almost_equal", "assert_array_less",
           "assert_less", "assert_less_equal",
           "assert_greater", "assert_greater_equal",
           "assert_approx_equal", "SkipTest"]


_dummy = unittest.TestCase('__init__')
assert_equal = _dummy.assertEqual
assert_not_equal = _dummy.assertNotEqual
assert_true = _dummy.assertTrue
assert_false = _dummy.assertFalse
assert_raises = _dummy.assertRaises

try:
    SkipTest = unittest.case.SkipTest
except AttributeError:
    # Python <= 2.6, we stil need nose here
    from nose import SkipTest


try:
    assert_dict_equal = _dummy.assertDictEqual
    assert_in = _dummy.assertIn
    assert_not_in = _dummy.assertNotIn
except AttributeError:
    # Python <= 2.6

    assert_dict_equal = assert_equal

    def assert_in(x, container):
        assert_true(x in container, msg="%r in %r" % (x, container))

    def assert_not_in(x, container):
        assert_false(x in container, msg="%r in %r" % (x, container))

try:
    assert_raises_regex = _dummy.assertRaisesRegex
except AttributeError:
    # for Python 2.6
    def assert_raises_regex(expected_exception, expected_regexp,
                            callable_obj=None, *args, **kwargs):
        """Helper function to check for message patterns in exceptions."""
        not_raised = False
        try:
            callable_obj(*args, **kwargs)
            not_raised = True
        except expected_exception as e:
            error_message = str(e)
            if not re.compile(expected_regexp).search(error_message):
                raise AssertionError("Error message should match pattern "
                                     "%r. %r does not." %
                                     (expected_regexp, error_message))
        if not_raised:
            raise AssertionError("%s not raised by %s" %
                                 (expected_exception.__name__,
                                  callable_obj.__name__))

# assert_raises_regexp is deprecated in Python 3.4 in favor of
# assert_raises_regex but lets keep the backward compat in scikit-learn with
# the old name for now
assert_raises_regexp = assert_raises_regex


def _assert_less(a, b, msg=None):
    message = "%r is not lower than %r" % (a, b)
    if msg is not None:
        message += ": " + msg
    assert a < b, message


def _assert_greater(a, b, msg=None):
    message = "%r is not greater than %r" % (a, b)
    if msg is not None:
        message += ": " + msg
    assert a > b, message


def assert_less_equal(a, b, msg=None):
    message = "%r is not lower than or equal to %r" % (a, b)
    if msg is not None:
        message += ": " + msg
    assert a <= b, message


def assert_greater_equal(a, b, msg=None):
    message = "%r is not greater than or equal to %r" % (a, b)
    if msg is not None:
        message += ": " + msg
    assert a >= b, message


def assert_warns(warning_class, func, *args, **kw):
    """Test that a certain warning occurs.

    Parameters
    ----------
    warning_class : the warning class
        The class to test for, e.g. UserWarning.

    func : callable
        Calable object to trigger warnings.

    *args : the positional arguments to `func`.

    **kw : the keyword arguments to `func`

    Returns
    -------

    result : the return value of `func`

    """
    # very important to avoid uncontrolled state propagation
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
        The entire message or a substring to  test for. If callable,
        it takes a string as argument and will trigger an assertion error
        if it returns `False`.

    func : callable
        Calable object to trigger warnings.

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


# To remove when we support numpy 1.7
def assert_no_warnings(func, *args, **kw):
    # XXX: once we may depend on python >= 2.6, this can be replaced by the

    # warnings module context manager.
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

    Note. Using this (in both variants) will clear all warnings
    from all python modules loaded. In case you need to test
    cross-module-warning-logging this is not your tool of choice.

    Parameters
    ----------
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
    if callable(obj):
        return _IgnoreWarnings(category=category)(obj)
    else:
        return _IgnoreWarnings(category=category)


class _IgnoreWarnings(object):
    """Improved and simplified Python warnings context manager and decorator.

    This class allows to ignore the warnings raise by a function.
    Copied from Python 2.7.5 and modified as required.

    Parameters
    ----------
    category : tuple of warning class, defaut to Warning
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
            # very important to avoid uncontrolled state propagation
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
        clean_warning_registry()  # be safe and not propagate state + chaos
        warnings.simplefilter("ignore", self.category)
        if self._entered:
            raise RuntimeError("Cannot enter %r twice" % self)
        self._entered = True
        self._filters = self._module.filters
        self._module.filters = self._filters[:]
        self._showwarning = self._module.showwarning

    def __exit__(self, *exc_info):
        if not self._entered:
            raise RuntimeError("Cannot exit %r without entering first" % self)
        self._module.filters = self._filters
        self._module.showwarning = self._showwarning
        self.log[:] = []
        clean_warning_registry()  # be safe and not propagate state + chaos


try:
    assert_less = _dummy.assertLess
    assert_greater = _dummy.assertGreater
except AttributeError:
    assert_less = _assert_less
    assert_greater = _assert_greater


def _assert_allclose(actual, desired, rtol=1e-7, atol=0,
                     err_msg='', verbose=True):
    actual, desired = np.asanyarray(actual), np.asanyarray(desired)
    if np.allclose(actual, desired, rtol=rtol, atol=atol):
        return
    msg = ('Array not equal to tolerance rtol=%g, atol=%g: '
           'actual %s, desired %s') % (rtol, atol, actual, desired)
    raise AssertionError(msg)


if hasattr(np.testing, 'assert_allclose'):
    assert_allclose = np.testing.assert_allclose
else:
    assert_allclose = _assert_allclose


def assert_raise_message(exceptions, message, function, *args, **kwargs):
    """Helper function to test error messages in exceptions.

    Parameters
    ----------
    exceptions : exception or tuple of exception
        Name of the estimator

    function : callable
        Calable object to raise error

    *args : the positional arguments to `function`.

    **kw : the keyword arguments to `function`
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


def fake_mldata(columns_dict, dataname, matfile, ordering=None):
    """Create a fake mldata data set.

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


class mock_mldata_urlopen(object):

    def __init__(self, mock_datasets):
        """Object that mocks the urlopen function to fake requests to mldata.

        `mock_datasets` is a dictionary of {dataset_name: data_dict}, or
        {dataset_name: (data_dict, ordering).
        `data_dict` itself is a dictionary of {column_name: data_array},
        and `ordering` is a list of column_names to determine the ordering
        in the data set (see `fake_mldata` for details).

        When requesting a dataset with a name that is in mock_datasets,
        this object creates a fake dataset in a StringIO object and
        returns it. Otherwise, it raises an HTTPError.
        """
        self.mock_datasets = mock_datasets

    def __call__(self, urlname):
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
    # Lazy import to avoid mutually recursive imports
    from sklearn import datasets
    datasets.mldata.urlopen = mock_mldata_urlopen(mock_datasets)


def uninstall_mldata_mock():
    # Lazy import to avoid mutually recursive imports
    from sklearn import datasets
    datasets.mldata.urlopen = urlopen


# Meta estimators need another estimator to be instantiated.
META_ESTIMATORS = ["OneVsOneClassifier", "MultiOutputEstimator",
                   "MultiOutputRegressor", "MultiOutputClassifier",
                   "OutputCodeClassifier", "OneVsRestClassifier",
                   "RFE", "RFECV", "BaseEnsemble"]
# estimators that there is no way to default-construct sensibly
OTHER = ["Pipeline", "FeatureUnion", "GridSearchCV", "RandomizedSearchCV",
         "SelectFromModel"]

# some trange ones
DONT_TEST = ['SparseCoder', 'EllipticEnvelope', 'DictVectorizer',
             'LabelBinarizer', 'LabelEncoder',
             'MultiLabelBinarizer', 'TfidfTransformer',
             'TfidfVectorizer', 'IsotonicRegression',
             'OneHotEncoder', 'RandomTreesEmbedding',
             'FeatureHasher', 'DummyClassifier', 'DummyRegressor',
             'TruncatedSVD', 'PolynomialFeatures',
             'GaussianRandomProjectionHash', 'HashingVectorizer',
             'CheckingClassifier', 'PatchExtractor', 'CountVectorizer',
             # GradientBoosting base estimators, maybe should
             # exclude them in another way
             'ZeroEstimator', 'ScaledLogOddsEstimator',
             'QuantileEstimator', 'MeanEstimator',
             'LogOddsEstimator', 'PriorProbabilityEstimator',
             '_SigmoidCalibration', 'VotingClassifier']


def all_estimators(include_meta_estimators=False,
                   include_other=False, type_filter=None,
                   include_dont_test=False):
    """Get a list of all estimators from sklearn.

    This function crawls the module and gets all classes that inherit
    from BaseEstimator. Classes that are defined in test-modules are not
    included.
    By default meta_estimators such as GridSearchCV are also not included.

    Parameters
    ----------
    include_meta_estimators : boolean, default=False
        Whether to include meta-estimators that can be constructed using
        an estimator as their first argument. These are currently
        BaseEnsemble, OneVsOneClassifier, OutputCodeClassifier,
        OneVsRestClassifier, RFE, RFECV.

    include_other : boolean, default=False
        Wether to include meta-estimators that are somehow special and can
        not be default-constructed sensibly. These are currently
        Pipeline, FeatureUnion and GridSearchCV

    include_dont_test : boolean, default=False
        Whether to include "special" label estimator or test processors.

    type_filter : string, list of string,  or None, default=None
        Which kind of estimators should be returned. If None, no filter is
        applied and all estimators are returned.  Possible values are
        'classifier', 'regressor', 'cluster' and 'transformer' to get
        estimators only of these specific types, or a list of these to
        get the estimators that fit at least one of the types.

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

    all_classes = []
    # get parent folder
    path = sklearn.__path__
    for importer, modname, ispkg in pkgutil.walk_packages(
            path=path, prefix='sklearn.', onerror=lambda x: None):
        if (".tests." in modname):
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

    if not include_dont_test:
        estimators = [c for c in estimators if not c[0] in DONT_TEST]

    if not include_other:
        estimators = [c for c in estimators if not c[0] in OTHER]
    # possibly get rid of meta estimators
    if not include_meta_estimators:
        estimators = [c for c in estimators if not c[0] in META_ESTIMATORS]
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

    Classes for whom random_state is deprecated are ignored. Currently DBSCAN
    is one such class.
    """
    if isinstance(estimator, DBSCAN):
        return

    if "random_state" in estimator.get_params():
        estimator.set_params(random_state=random_state)


def if_matplotlib(func):
    """Test decorator that skips test if matplotlib not installed."""
    @wraps(func)
    def run_test(*args, **kwargs):
        try:
            import matplotlib
            matplotlib.use('Agg', warn=False)
            # this fails if no $DISPLAY specified
            import matplotlib.pyplot as plt
            plt.figure()
        except ImportError:
            raise SkipTest('Matplotlib not available.')
        else:
            return func(*args, **kwargs)
    return run_test


def skip_if_32bit(func):
    """Test decorator that skips tests on 32bit platforms."""
    @wraps(func)
    def run_test(*args, **kwargs):
        bits = 8 * struct.calcsize("P")
        if bits == 32:
            raise SkipTest('Test skipped on 32bit platforms.')
        else:
            return func(*args, **kwargs)
    return run_test


def if_not_mac_os(versions=('10.7', '10.8', '10.9'),
                  message='Multi-process bug in Mac OS X >= 10.7 '
                          '(see issue #636)'):
    """Test decorator that skips test if OS is Mac OS X and its
    major version is one of ``versions``.
    """
    warnings.warn("if_not_mac_os is deprecated in 0.17 and will be removed"
                  " in 0.19: use the safer and more generic"
                  " if_safe_multiprocessing_with_blas instead",
                  DeprecationWarning)
    mac_version, _, _ = platform.mac_ver()
    skip = '.'.join(mac_version.split('.')[:2]) in versions

    def decorator(func):
        if skip:
            @wraps(func)
            def func(*args, **kwargs):
                raise SkipTest(message)
        return func
    return decorator


def if_safe_multiprocessing_with_blas(func):
    """Decorator for tests involving both BLAS calls and multiprocessing.

    Under POSIX (e.g. Linux or OSX), using multiprocessing in conjunction with
    some implementation of BLAS (or other libraries that manage an internal
    posix thread pool) can cause a crash or a freeze of the Python process.

    In practice all known packaged distributions (from Linux distros or
    Anaconda) of BLAS under Linux seems to be safe. So we this problem seems to
    only impact OSX users.

    This wrapper makes it possible to skip tests that can possibly cause
    this crash under OS X with.

    Under Python 3.4+ it is possible to use the `forkserver` start method
    for multiprocessing to avoid this issue. However it can cause pickling
    errors on interactively defined functions. It therefore not enabled by
    default.
    """
    @wraps(func)
    def run_test(*args, **kwargs):
        if sys.platform == 'darwin':
            raise SkipTest(
                "Possible multi-process bug with some BLAS")
        return func(*args, **kwargs)
    return run_test


def clean_warning_registry():
    """Safe way to reset warnings."""
    warnings.resetwarnings()
    reg = "__warningregistry__"
    for mod_name, mod in list(sys.modules.items()):
        if 'six.moves' in mod_name:
            continue
        if hasattr(mod, reg):
            getattr(mod, reg).clear()


def check_skip_network():
    if int(os.environ.get('SKLEARN_SKIP_NETWORK_TESTS', 0)):
        raise SkipTest("Text tutorial requires large dataset download")


def check_skip_travis():
    """Skip test if being run on Travis."""
    if os.environ.get('TRAVIS') == "true":
        raise SkipTest("This test needs to be skipped on Travis")


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


class TempMemmap(object):
    def __init__(self, data, mmap_mode='r'):
        self.temp_folder = tempfile.mkdtemp(prefix='sklearn_testing_')
        self.mmap_mode = mmap_mode
        self.data = data

    def __enter__(self):
        fpath = op.join(self.temp_folder, 'data.pkl')
        joblib.dump(self.data, fpath)
        data_read_only = joblib.load(fpath, mmap_mode=self.mmap_mode)
        atexit.register(lambda: _delete_folder(self.temp_folder, warn=True))
        return data_read_only

    def __exit__(self, exc_type, exc_val, exc_tb):
        _delete_folder(self.temp_folder)


with_network = with_setup(check_skip_network)
with_travis = with_setup(check_skip_travis)


class _named_check(object):
    """Wraps a check to show a useful description

    Parameters
    ----------
    check : function
        Must have ``__name__`` and ``__call__``
    arg_text : str
        A summary of arguments to the check
    """
    # Setting the description on the function itself can give incorrect results
    # in failing tests
    def __init__(self, check, arg_text):
        self.check = check
        self.description = ("{0[1]}.{0[3]}:{1.__name__}({2})".format(
            inspect.stack()[1], check, arg_text))

    def __call__(self, *args, **kwargs):
        return self.check(*args, **kwargs)
