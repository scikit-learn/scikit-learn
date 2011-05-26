"""Testing utilities."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD

import os
import scipy as sp
from scipy.io import savemat


def assert_in(obj, in_=None, out_=None):
    """Checks that all names in `in_` as in `obj`, but no name
    in `out_` is."""
    if in_ is not None:
        for name in in_:
            assert name in obj
    if out_ is not None:
        for name in out_:
            assert name not in obj


def fake_mldata_cache(columns_dict, dataname, cachedir, ordering=None):
    """Create a fake mldata data set in the cache_path.

    Parameters
    ----------
    columns_dict: contains data as
                  columns_dict[column_name] = array of data
    dataname: name of data set
    cachedir: cache path
    ordering: list of column_names, determines the ordering in the data set

    Note: this function transposes all arrays, while fetch_mldata only
    transposes 'data', keep that into account in the tests.
    """
    dset = dict(columns_dict)

    # transpose all variables
    for name in dset:
        dset[name] = dset[name].T

    if ordering is None:
        ordering = dset.keys()
    # NOTE: setting up this array is tricky, because of the way Matlab
    # re-packages 1D arrays
    dset['mldata_descr_ordering'] = sp.empty((1, len(ordering)),
                                             dtype='object')
    for i, name in enumerate(ordering):
        dset['mldata_descr_ordering'][0, i] = name

    savemat(os.path.join(cachedir, 'mldata', dataname), dset,
            oned_as='column')


class mock_urllib2(object):

    def __init__(self, mock_datasets, cachedir):
        """Object that mocks the urllib2 module to fake requests to mldata.

        `mock_datasets` is a dictionary of {dataset_name: data_dict}, and
        `data_dict` itself is a dictionary of {column_name: data_array}

        When requesting a dataset with a name that is in mock_datasets,
        this object creates a fake dataset in cachedir,
        then returns a file handle to it. Otherwise, it raises
        an URLError.
        """
        self.mock_datasets = mock_datasets
        self.cachedir = cachedir

    class URLError(Exception):
        pass

    def urlopen(self, urlname):
        dataset_name = urlname.split('/')[-1]
        if dataset_name in self.mock_datasets:
            resource_name = '_' + dataset_name
            fake_mldata_cache(self.mock_datasets[dataset_name],
                              resource_name, self.cachedir)
            return open(os.path.join(self.cachedir, 'mldata',
                                     resource_name + '.mat'), 'rb')
        else:
            raise mock_urllib2.URLError('%s not found.', urlname)
