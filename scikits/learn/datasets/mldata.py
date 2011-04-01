"""Automatically download MLdata datasets."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD

from scipy import io

from os.path import join, exists
from os import makedirs
import urllib2

from .base import get_data_home, Bunch

# TODO: test: it download the first time, it loads it the second
# XXX: is the assignment of data and labels the mldata.org standard?

MLDATA_BASE_URL = "http://mldata.org/repository/data/download/matlab/%s"

def fetch_mldata(dataname, data_home=None):
    """
    Fecth an mldata.org data set.
    
    If the file does not exist yet, it is downloaded from mldata.org .

    Parameters
    ----------
    dataname:
        Name of the dataset on mldata.org,
        e.g.: "leukemia", "Whistler Daily Snowfall", etc.
    
    data_home: optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels,
        and 'DESCR', the full description of the dataset.
    
    Example
    -------
    
    Load the 'leukemia' dataset from mldata.org, and print the number of
    data points:
    
    >>> data = fetch_mldata('leukemia')
    >>> print data.data.shape[0]
    7129
    """

    # normalize dataset name
    dataname = dataname.lower().replace(' ', '-').replace('.', '')

    # check if this data set has been already downloaded
    data_home = get_data_home(data_home=data_home)
    data_home = join(data_home, 'mldata')
    if not exists(data_home):
        makedirs(data_home)

    matlab_name = dataname + '.mat'
    filename = join(data_home, matlab_name)
    print matlab_name, filename

    # if the file does not exist, download it
    if not exists(filename):
        print 'download'
        urlname = MLDATA_BASE_URL % (dataname)
        try:
            mldata_url = urllib2.urlopen(urlname)
        except urllib2.URLError as e:
            msg = "Dataset '%s' not found on mldata.org." % dataname
            raise IOError(msg)
        # store Matlab file
        with open(filename, 'w+b') as matlab_file:
            matlab_file.write(mldata_url.read())
        mldata_url.close()

    # load dataset matlab file
    with open(filename, 'rb') as matlab_file:
        matlab_dict = io.loadmat(matlab_file, struct_as_record=True)

    # extract data from matlab_dict
    ordering = matlab_dict['mldata_descr_ordering']
    # no labels
    if ordering.shape[1] == 1:
        data_name = ordering[0, 0][0]
        return Bunch(data=matlab_dict[data_name],
                     DESCR='mldata.org dataset: %s' % dataname)
    # with labels
    label_name = ordering[0, 0][0]
    data_name = ordering[0, 1][0]
    return Bunch(data=matlab_dict[data_name],
                 target=matlab_dict[label_name],
                 DESCR='mldata.org dataset: %s' % dataname)
