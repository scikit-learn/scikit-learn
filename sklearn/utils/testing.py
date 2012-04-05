"""Testing utilities."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD

from .fixes import savemat
import urllib2
from StringIO import StringIO
import scipy as sp

import nose.tools

_NOSE_TOOLS = dir(nose.tools)

if 'assert_in' in _NOSE_TOOLS:
    assert_in = nose.tools.assert_in
else:
    def assert_in(member, container, msg=None):
        """
        Checks member is present in container
        definition of assert_in if missing in nose.tools

        Parameters
        ----------
        member: object
        container: container or iterator
        msg: message if assertion fails
        """
        if not msg:
            msg = "%s not in %s" % (member, container)
        nose.tools.assert_true((member in container), msg=msg)


if 'assert_not_in' in _NOSE_TOOLS:
    assert_in = nose.tools.assert_in
else:
    def assert_not_in(member, container, msg=None):
        """
        Checks member is NOT in container
        definition of assert_not_in if missing in nose.tools

        Parameters
        ----------
        member: object
        container: container or iterator
        msg: message if assertion fails
        """
        if not msg:
            msg = "%s in %s" % (member, container)
        nose.tools.assert_false((member in container), msg=msg)


def assert_contained(obj, in_=None, out_=None):
    """Checks
    * that all objects in `in_` are in `obj`,
    * that all objects in `out` are NOT in `object`

    Parameters
    ----------
    obj: container or iterable (implements __contains__ or __iter__).
    in_: iterable or None (default)
        if not None, the items in this iterable must be in obj
        for the assert to succeed
    out_: iterable or None (default)
        if not None the items in this iterable must NOT be in obj
        for the assert to succeed
    """
    if in_ is not None:
        for item in in_:
            assert_in(item, obj)
    if out_ is not None:
        for item in out_:
            assert_not_in(item, obj)


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


# objects that are imported from this module by default
__all__ = assert_contained, assert_in, assert_not_in, fake_mldata_cache, \
    mock_urllib2
