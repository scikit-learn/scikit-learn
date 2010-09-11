"""
Base IO code for all datasets
"""

# Copyright (c) 2007 David Cournapeau <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
# License: Simplified BSD

import csv
import os

import numpy as np


class Bunch(dict):
    """ Container object for datasets: dictionnary-like object that
        exposes its keys as attributes.
    """

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self


################################################################################

def load_iris():
    """load the iris dataset and returns it.
    
    Returns
    -------
    data : Bunch
        Dictionnary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels, 
        'target_names', the meaning of the labels, and 'DESCR', the
        full description of the dataset.

    Example
    -------
    Let's say you are interested in the samples 10, 25, and 50, and want to
    know their class name.

    >>> from scikits.learn.datasets import load_iris
    >>> data = load_iris()
    >>> data.target[[10, 25, 50]]
    array([0, 0, 1])
    >>> data.target_names
    array(['setosa', 'versicolor', 'virginica'],
          dtype='|S10')

    """
    
    data_file = csv.reader(open(os.path.dirname(__file__) 
                        + '/data/iris.csv'))
    fdescr = open(os.path.dirname(__file__) 
                        + '/descr/iris.rst')
    temp = data_file.next()
    n_samples = int(temp[0])
    n_features = int(temp[1])
    target_names = np.array(temp[2:])
    data = np.empty((n_samples, n_features))
    target = np.empty((n_samples,), dtype=np.int)
    for i, ir in enumerate(data_file):
        data[i] = np.asanyarray(ir[:-1], dtype=np.float)
        target[i] = np.asanyarray(ir[-1], dtype=np.int)
    return Bunch(data=data, target=target, target_names=target_names, 
                 DESCR=fdescr.read())


def load_digits():
    """load the digits dataset and returns it.
    
    Returns
    -------
    data : Bunch
        Dictionnary-like object, the interesting attributes are:
        'data', the data to learn, `images`, the images corresponding
        to each sample, 'target', the classification labels for each
        sample, 'target_names', the meaning of the labels, and 'DESCR', 
        the full description of the dataset.

    Example
    -------
    To load the data and visualize the images::

        import pylab as pl
        digits = datasets.load_digits()
        pl.gray()
        # Visualize the first image:
        pl.matshow(digits.raw_data[0])

    """
    
    data = np.loadtxt(os.path.join(os.path.dirname(__file__) 
                        + '/data/digits.csv.gz'), delimiter=',')
    fdescr = open(os.path.join(os.path.dirname(__file__) 
                        + '/descr/digits.rst'))
    target = data[:, -1]
    flat_data = data[:, :-1]
    images = flat_data.view()
    images.shape = (-1, 8, 8)
    return Bunch(data=flat_data, target=target.astype(np.int), 
                 target_names=np.arange(10), 
                 images=images,
                 DESCR=fdescr.read())



def load_diabetes():
    data = np.loadtxt(os.path.join(os.path.dirname(__file__) +
                                   '/data/diabetes.csv'))
    target = data[:, -1]
    data   = data[:, :-1]
    return Bunch (data=data, target=target)
