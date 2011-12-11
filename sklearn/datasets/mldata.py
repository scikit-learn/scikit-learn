"""Automatically download MLdata datasets."""

# Copyright (c) 2011 Pietro Berkes
# License: Simplified BSD

import os
from os.path import join, exists
import re
from scipy import io
from shutil import copyfileobj
import urllib2

from .base import get_data_home, Bunch

MLDATA_BASE_URL = "http://mldata.org/repository/data/download/matlab/%s"


def mldata_filename(dataname):
    """Convert a raw name for a data set in a mldata.org filename."""
    dataname = dataname.lower().replace(' ', '-')
    return re.sub(r'[().]', '', dataname)


def fetch_mldata(dataname, target_name='label', data_name='data',
                 transpose_data=True, data_home=None):
    """Fetch an mldata.org data set

    If the file does not exist yet, it is downloaded from mldata.org .

    mldata.org does not have an enforced convention for storing data or
    naming the columns in a data set. The default behavior of this function
    works well with the most common cases:
      1) data values are stored in the column 'data', and target values in the
         column 'label'
      2) alternatively, the first column stores target values, and the second
         data values
      3) the data array is stored as `n_features x n_samples` , and thus needs
         to be transposed to match the `sklearn` standard

    Keyword arguments allow to adapt these defaults to specific data sets
    (see parameters `target_name`, `data_name`, `transpose_data`, and
    the examples below).

    mldata.org data sets may have multiple columns, which are stored in the
    Bunch object with their original name.

    Parameters
    ----------

    dataname:
        Name of the data set on mldata.org,
        e.g.: "leukemia", "Whistler Daily Snowfall", etc.
        The raw name is automatically converted to a mldata.org URL .

    target_name: optional, default: 'label'
        Name or index of the column containing the target values.

    data_name: optional, default: 'data'
        Name or index of the column containing the data.

    transpose_data: optional, default: True
        If True, transpose the downloaded data array.

    data_home: optional, default: None
        Specify another download and cache folder for the data sets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    Returns
    -------

    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels,
        'DESCR', the full description of the dataset, and
        'COL_NAMES', the original names of the dataset columns.

    Examples
    --------
    Load the 'iris' dataset from mldata.org:
    >>> from sklearn.datasets.mldata import fetch_mldata
    >>> iris = fetch_mldata('iris')
    >>> iris.target[0]
    1
    >>> print iris.data[0]
    [-0.555556  0.25     -0.864407 -0.916667]

    Load the 'leukemia' dataset from mldata.org, which respects the
    sklearn axes convention:
    >>> leuk = fetch_mldata('leukemia', transpose_data=False)
    >>> print leuk.data.shape[0]
    7129

    Load an alternative 'iris' dataset, which has different names for the
    columns:
    >>> iris2 = fetch_mldata('datasets-UCI iris', target_name=1,
    ...                      data_name=0)
    >>> iris3 = fetch_mldata('datasets-UCI iris',
    ...                      target_name='class', data_name='double0')
    """

    # normalize dataset name
    dataname = mldata_filename(dataname)

    # check if this data set has been already downloaded
    data_home = get_data_home(data_home=data_home)
    data_home = join(data_home, 'mldata')
    if not exists(data_home):
        os.makedirs(data_home)

    matlab_name = dataname + '.mat'
    filename = join(data_home, matlab_name)

    # if the file does not exist, download it
    if not exists(filename):
        urlname = MLDATA_BASE_URL % urllib2.quote(dataname)
        try:
            mldata_url = urllib2.urlopen(urlname)
        except urllib2.HTTPError as e:
            if e.code == 404:
                e.msg = "Dataset '%s' not found on mldata.org." % dataname
            raise
        # store Matlab file
        try:
            with open(filename, 'w+b') as matlab_file:
                copyfileobj(mldata_url, matlab_file)
        except:
            os.remove(filename)
            raise
        mldata_url.close()

    # load dataset matlab file
    with open(filename, 'rb') as matlab_file:
        matlab_dict = io.loadmat(matlab_file, struct_as_record=True)

    # -- extract data from matlab_dict

    # flatten column names
    col_names = [str(descr[0])
                 for descr in matlab_dict['mldata_descr_ordering'][0]]

    # if target or data names are indices, transform then into names
    if isinstance(target_name, int):
        target_name = col_names[target_name]
    if isinstance(data_name, int):
        data_name = col_names[data_name]

    # rules for making sense of the mldata.org data format
    # (earlier ones have priority):
    # 1) there is only one array => it is "data"
    # 2) there are multiple arrays
    #    a) copy all columns in the bunch, using their column name
    #    b) if there is a column called `target_name`, set "target" to it,
    #        otherwise set "target" to first column
    #    c) if there is a column called `data_name`, set "data" to it,
    #        otherwise set "data" to second column

    dataset = {'DESCR': 'mldata.org dataset: %s' % dataname,
               'COL_NAMES': col_names}

    # 1) there is only one array => it is considered data
    if len(col_names) == 1:
        data_name = col_names[0]
        dataset['data'] = matlab_dict[data_name]
    # 2) there are multiple arrays
    else:
        for name in col_names:
            dataset[name] = matlab_dict[name]

        if target_name in col_names:
            del dataset[target_name]
            dataset['target'] = matlab_dict[target_name]
        else:
            del dataset[col_names[0]]
            dataset['target'] = matlab_dict[col_names[0]]

        if data_name in col_names:
            del dataset[data_name]
            dataset['data'] = matlab_dict[data_name]
        else:
            del dataset[col_names[1]]
            dataset['data'] = matlab_dict[col_names[1]]

    # set axes to sklearn conventions
    if transpose_data:
        dataset['data'] = dataset['data'].T
    if 'target' in dataset:
        dataset['target'] = dataset['target'].squeeze()

    return Bunch(**dataset)
