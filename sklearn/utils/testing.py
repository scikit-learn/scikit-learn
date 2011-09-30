"""Testing utilities."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD
import sys
import warnings
from .fixes import savemat
import urllib2
from StringIO import StringIO
import scipy as sp
import numpy.testing


def assert_in(obj, in_=None, out_=None):
    """Checks that all names in `in_` as in `obj`, but no name
    in `out_` is."""
    if in_ is not None:
        for name in in_:
            assert name in obj
    if out_ is not None:
        for name in out_:
            assert name not in obj


def fake_mldata_cache(columns_dict, dataname, matfile, ordering=None):
    """Create a fake mldata data set in the cache_path.

    Parameters
    ----------
    columns_dict: contains data as
                  columns_dict[column_name] = array of data
    dataname: name of data set
    matfile: file-like object or file name
    ordering: list of column_names, determines the ordering in the data set

    Note: this function transposes all arrays, while fetch_mldata only
    transposes 'data', keep that into account in the tests.
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

    savemat(matfile, datasets, oned_as='column')


class mock_urllib2(object):

    def __init__(self, mock_datasets):
        """Object that mocks the urllib2 module to fake requests to mldata.

        `mock_datasets` is a dictionary of {dataset_name: data_dict}, or
        {dataset_name: (data_dict, ordering).
        `data_dict` itself is a dictionary of {column_name: data_array},
        and `ordering` is a list of column_names to determine the ordering
        in the data set (see `fake_mldata_cache` for details).

        When requesting a dataset with a name that is in mock_datasets,
        this object creates a fake dataset in a StringIO object and
        returns it. Otherwise, it raises an URLError.
        """
        self.mock_datasets = mock_datasets

    class HTTPError(urllib2.URLError):
        code = 404

    def urlopen(self, urlname):
        dataset_name = urlname.split('/')[-1]
        if dataset_name in self.mock_datasets:
            resource_name = '_' + dataset_name
            matfile = StringIO()

            dataset = self.mock_datasets[dataset_name]
            ordering = None
            if isinstance(dataset, tuple):
                dataset, ordering = dataset
            fake_mldata_cache(dataset, resource_name, matfile, ordering)

            matfile.seek(0)
            return matfile
        else:
            raise mock_urllib2.HTTPError('%s not found.' % urlname)

    def quote(self, string, safe='/'):
        return urllib2.quote(string, safe)


#----------------------------------------------------------------------
# These are copied from python 2.6 warnings module
# It is also available in numpy > 1.3

class WarningMessage(object):

    """
    Holds the result of a single showwarning() call.

    Notes
    -----
    `WarningMessage` is copied from the Python 2.6 warnings module,
    so it can be used in NumPy with older Python versions.

    """

    _WARNING_DETAILS = ("message", "category", "filename", "lineno", "file",
                        "line")

    def __init__(self, message, category, filename, lineno, file=None,
                    line=None):
        local_values = locals()
        for attr in self._WARNING_DETAILS:
            setattr(self, attr, local_values[attr])
        if category:
            self._category_name = category.__name__
        else:
            self._category_name = None

    def __str__(self):
        return ("{message : %r, category : %r, filename : %r, lineno : %s, "
                    "line : %r}" % (self.message, self._category_name,
                                    self.filename, self.lineno, self.line))


class WarningManager(object):
    """
    A context manager that copies and restores the warnings filter upon
    exiting the context.

    The 'record' argument specifies whether warnings should be captured by a
    custom implementation of ``warnings.showwarning()`` and be appended to a
    list returned by the context manager. Otherwise None is returned by the
    context manager. The objects appended to the list are arguments whose
    attributes mirror the arguments to ``showwarning()``.

    The 'module' argument is to specify an alternative module to the module
    named 'warnings' and imported under that name. This argument is only useful
    when testing the warnings module itself.

    Notes
    -----
    `WarningManager` is a copy of the ``catch_warnings`` context manager
    from the Python 2.6 warnings module, with slight modifications.
    It is copied so it can be used in NumPy with older Python versions.

    """
    def __init__(self, record=False, module=None):
        self._record = record
        if module is None:
            self._module = sys.modules['warnings']
        else:
            self._module = module
        self._entered = False

    def __enter__(self):
        if self._entered:
            raise RuntimeError("Cannot enter %r twice" % self)
        self._entered = True
        self._filters = self._module.filters
        self._module.filters = self._filters[:]
        self._showwarning = self._module.showwarning
        if self._record:
            log = []

            def showwarning(*args, **kwargs):
                log.append(WarningMessage(*args, **kwargs))
            self._module.showwarning = showwarning
            return log
        else:
            return None

    def __exit__(self):
        if not self._entered:
            raise RuntimeError("Cannot exit %r without entering first" % self)
        self._module.filters = self._filters
        self._module.showwarning = self._showwarning


def assert_warns(warning_class, func, *args, **kw):
    """
    Fail unless the given callable throws the specified warning.

    A warning of class warning_class should be thrown by the callable when
    invoked with arguments args and keyword arguments kwargs.
    If a different type of warning is thrown, it will not be caught, and the
    test case will be deemed to have suffered an error.

    Parameters
    ----------
    warning_class : class
        The class defining the warning that `func` is expected to throw.
    func : callable
        The callable to test.
    \\*args : Arguments
        Arguments passed to `func`.
    \\*\\*kwargs : Kwargs
        Keyword arguments passed to `func`.

    Returns
    -------
    None

    """
    ctx = WarningManager(record=True)
    l = ctx.__enter__()
    warnings.simplefilter('always')
    try:
        func(*args, **kw)
        if not len(l) > 0:
            raise AssertionError("No warning raised when calling %s"
                    % func.__name__)
        if not l[0].category is warning_class:
            raise AssertionError("First warning for %s is not a " \
                    "%s( is %s)" % (func.__name__, warning_class, l[0]))
    finally:
        ctx.__exit__()
