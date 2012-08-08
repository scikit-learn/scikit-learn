"""Testing utilities."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD
import inspect
import pkgutil

import urllib2
from StringIO import StringIO
import scipy as sp
import sklearn
from sklearn.base import BaseEstimator
from .fixes import savemat


try:
    from nose.tools import assert_in, assert_not_in
except ImportError:
    # Nose < 1.0.0
    from nose.tools import assert_true, assert_false

    def assert_in(x, container):
        assert_true(x in container, msg="%r in %r" % (x, container))

    def assert_not_in(x, container):
        assert_false(x in container, msg="%r in %r" % (x, container))


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


try:
    from nose.tools import assert_less
except ImportError:
    assert_less = _assert_less

try:
    from nose.tools import assert_greater
except ImportError:
    assert_greater = _assert_greater


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


def all_estimators():
    def is_abstract(c):
        if not(hasattr(c, '__abstractmethods__')):
            return False
        if not len(c.__abstractmethods__):
            return False
        return True

    all_classes = []
    # get parent folder
    path = sklearn.__path__
    for importer, modname, ispkg in pkgutil.walk_packages(path=path,
                            prefix='sklearn.', onerror=lambda x: None):
        module = __import__(modname, fromlist="dummy")
        classes = inspect.getmembers(module, inspect.isclass)
        # get rid of abstract base classes
        all_classes.extend(classes)

    all_classes = set(all_classes)

    estimators = [c for c in all_classes if issubclass(c[1], BaseEstimator)]
    estimators = [c for c in estimators if not is_abstract(c[1])]
    # We sort in order to have reproducible test failures
    return sorted(estimators)
