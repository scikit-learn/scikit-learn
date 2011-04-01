"""Automatically download MLdata datasets."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD

import scipy as sp
import scipy.io
from os.path import join, exists
import urllib2
from .base import get_data_home, Bunch

# TODO: test: it download the first time, it loads it the second
# TODO: function to search dataset name and return best match
# TODO: error checking!

MLDATA_BASE_URL = "http://mldata.org/repository/data/download/matlab/%s"

def fetch_mldata(dataname, data_home=None):
    # check if this data set has been already downloaded
    data_home = get_data_home(data_home=data_home)
    matlab_name = dataname + '.mat'
    filename = join(data_home, matlab_name)
    print matlab_name, filename

    # if the file does not exist, download it
    if not exists(filename):
        print 'download'
        urlname = MLDATA_BASE_URL % (dataname)
        mldata_url = urllib2.urlopen(urlname)
        # store Matlab file
        with open(filename, 'w+b') as matlab_file:
            matlab_file.write(mldata_url.read())
        mldata_url.close()

    # load dataset
    with open(filename, 'rb') as matlab_file:
        matlab_dict = sp.io.loadmat(matlab_file)

    # create Bunch object
    return Bunch(data=matlab_dict['data'],
                 target=matlab_dict.get('label', None),
                 mldata_descr_ordering=matlab_dict['mldata_descr_ordering'],
                 DESCR='mldata.org dataset: %s' % dataname)
