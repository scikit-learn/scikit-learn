"""
Base object for all datasets
"""

# Copyright (c) 2007 David Cournapeau <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
#


import numpy as np

class Bunch(dict):
    """
    Container for dataset.

    Members
    -------
    - data : a record array with the actual data
    - label : label[i] = label index of data[i]
    - class : class[i] is the string corresponding to label index i.
    - COPYRIGHT, TITLE, SOURCE, DESCRSHORT, DESCRLONG,
      NOTE. Information about the dataset.
    """

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self

def load(dataset):
    """load the data and returns them.
    
    Returns
    -------
    data : Bunch
        See docstring of bunch for a complete description of its members.

    Available datasets
     - iris

    Example
    -------
    Let's say you are interested in the samples 10, 25, and 50, and want to
    know their class name.

    >>> data = load()
    >>> print data.label #doctest: +ELLIPSIS
    [ 0.  0. ...]
    """
    import csv
    import os
    
    firis = csv.reader(open(os.path.dirname(__file__) 
                        + '/data/%s.csv' % dataset))
    fdescr = open(os.path.dirname(__file__) 
                        + '/descr/%s.rst' % dataset)
    temp = firis.next()
    nsamples = int(temp[0])
    nfeat = int(temp[1])
    targetnames = temp[2:]
    data = np.empty((nsamples, nfeat))
    target = np.empty((nsamples,))
    for i, ir in enumerate(firis):
        data[i] = np.asanyarray(ir[:-1], dtype=np.float)
        target[i] = np.asanyarray(ir[-1], dtype=np.float)
    return Bunch(data=data, target=target, targetnames=targetnames, 
                 DESCR=fdescr.read())

