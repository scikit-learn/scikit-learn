"""
=============================
Species distribution dataset
=============================

This dataset represents the geographic distribution of species.
The dataset is provided by Phillips et. al. (2006).

The two species are:

 - `"Bradypus variegatus"
   <http://www.iucnredlist.org/apps/redlist/details/3038/0>`_ ,
   the Brown-throated Sloth.

 - `"Microryzomys minutus"
   <http://www.iucnredlist.org/apps/redlist/details/13408/0>`_ ,
   also known as the Forest Small Rice Rat, a rodent that lives in Peru,
   Colombia, Ecuador, Peru, and Venezuela.

References:

 * `"Maximum entropy modeling of species geographic distributions"
   <http://www.cs.princeton.edu/~schapire/papers/ecolmod.pdf>`_
   S. J. Phillips, R. P. Anderson, R. E. Schapire - Ecological Modelling,
   190:231-259, 2006.

Notes:

 * See examples/applications/plot_species_distribution_modeling.py
"""

# Authors: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
#
# License: Simplified BSD

from cStringIO import StringIO

import os
from os import makedirs
from os.path import join, normpath, split, exists

import urllib2

import pylab as pl
import numpy as np

from sklearn.datasets.base import get_data_home, Bunch

DIRECTORY_URL = "http://www.cs.princeton.edu/~schapire/maxent/datasets/" 

SAMPLES_URL = join(DIRECTORY_URL, "samples.zip")
COVERAGES_URL = join(DIRECTORY_URL, "coverages.zip")

DATA_ARCHIVE_NAME = "species_coverage.npz"

def load_coverage(F, header_length=6,
                  dtype=np.int16):
    """
    load a coverage file.
    This will return a numpy array of the given dtype
    """
    try:
        header = [F.readline() for i in range(header_length)]
    except:
        F = open(F)
        header = [F.readline() for i in range(header_length)]

    make_tuple = lambda t: (t.split()[0], float(t.split()[1]))
    header = dict([make_tuple(line) for line in header])

    M = np.loadtxt(F, dtype=dtype)
    nodata = header['NODATA_value']
    if nodata != -9999:
        M[nodata] = -9999
    
    return M


def load_csv(F):
    """Load csv file.
    
    Paramters
    ---------
    F : string or file object
        file object or name of file

    Returns
    -------
    rec : np.ndarray
        record array representing the data
    """
    try:
        names = F.readline().strip().split(',')
    except:
        F = open(F)
        names = F.readline().strip().split(',')

    rec = np.loadtxt(F, skiprows=1, delimiter=',',
                     dtype = 'a22,f4,f4')
    rec.dtype.names = names

    return rec

def save_compressed(filename,
                    x_left_lower_corner = -94.8, Nx=1212,
                    y_left_lower_corner = -56.05, Ny = 1592,
                    grid_size = 0.05, dtype=np.int16):
    print 80 * '_'
    print "loading coverages"
    coverages = [load_coverage(os.path.join('coverages', f), dtype=dtype)
                 for f in os.listdir('coverages')]
    test = load_csv('samples/alltrain.csv')
    train = load_csv('samples/alltest.csv')

    coverages = np.asarray(coverages,
                           dtype=dtype)

    np.savez(filename,
             coverages=coverages,
             test=test,
             train=train,
             x_left_lower_corner = x_left_lower_corner,
             y_left_lower_corner = y_left_lower_corner,
             Nx=Nx, Ny=Ny, grid_size=grid_size)


def load_compressed(filename):
    if not(os.path.exists(filename)):
        save_compressed(filename)
    X = np.load(filename)
    return Bunch(**dict([(f,X[f]) for f in X.files]))


def fetch_species_distributions(data_home=None,
                                download_if_missing=True):
    """Loader for species distribution dataset from Phillips et. al. (2006)

    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing: optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    Notes
    ------

    This dataset represents the geographic distribution of species.
    The dataset is provided by Phillips et. al. (2006).

    The two species are:

    - `"Bradypus variegatus"
      <http://www.iucnredlist.org/apps/redlist/details/3038/0>`_ ,
      the Brown-throated Sloth.

    - `"Microryzomys minutus"
      <http://www.iucnredlist.org/apps/redlist/details/13408/0>`_ ,
      also known as the Forest Small Rice Rat, a rodent that lives in Peru,
      Colombia, Ecuador, Peru, and Venezuela.

    References
    ----------

    * `"Maximum entropy modeling of species geographic distributions"
      <http://www.cs.princeton.edu/~schapire/papers/ecolmod.pdf>`_
      S. J. Phillips, R. P. Anderson, R. E. Schapire - Ecological Modelling,
      190:231-259, 2006.
    """
    # parameters for the data files.  These should not be changed unless
    # the data model changes.  They will be saved in the npz file with
    # the downloaded data.
    x_left_lower_corner = -94.8
    Nx=1212
    y_left_lower_corner = -56.05
    Ny = 1592
    grid_size = 0.05
    dtype = np.int16
    
    data_home = get_data_home(data_home)
    if not exists(data_home):
        makedirs(data_home)

    if not exists(join(data_home, DATA_ARCHIVE_NAME)):
        print 'Downloading species data from %s to %s' % (SAMPLES_URL,
                                                          data_home)
        
        X = np.load(StringIO(urllib2.urlopen(SAMPLES_URL).read()))

        for f in X.files:
            fhandle = StringIO(X[f])
            
            if 'train' in f:
                train = load_csv(fhandle)
            if 'test' in f:
                test = load_csv(fhandle)

        print 'Downloading coverage data from %s to %s' % (COVERAGES_URL,
                                                           data_home)

        X = np.load(StringIO(urllib2.urlopen(COVERAGES_URL).read()))

        coverages = []

        for f in sorted(X.files):
            fhandle = StringIO(X[f])
            print '-', f
            coverages.append(load_coverage(fhandle))

        coverages = np.asarray(coverages,
                               dtype=dtype)

        np.savez(join(data_home, DATA_ARCHIVE_NAME),
                 coverages=coverages,
                 test=test,
                 train=train,
                 x_left_lower_corner = x_left_lower_corner,
                 y_left_lower_corner = y_left_lower_corner,
                 Nx=Nx, Ny=Ny, grid_size=grid_size)

        bunch = Bunch(coverages=coverages,
                      test=test,
                      train=train,
                      x_left_lower_corner = x_left_lower_corner,
                      y_left_lower_corner = y_left_lower_corner,
                      Nx=Nx, Ny=Ny, grid_size=grid_size)
    else:
        X = np.load(join(data_home, DATA_ARCHIVE_NAME))
        return Bunch(**dict([(f,X[f]) for f in X.files]))
        

if __name__ == '__main__':
    fetch_species_distributions()
