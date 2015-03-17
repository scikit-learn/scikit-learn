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
   for an example of using this dataset
"""

# Authors: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
#
# License: BSD 3 clause

from io import BytesIO
from os import makedirs
from os.path import join
from os.path import exists

try:
    # Python 2
    from urllib2 import urlopen
    PY2 = True
except ImportError:
    # Python 3
    from urllib.request import urlopen
    PY2 = False

import numpy as np

from sklearn.datasets.base import get_data_home, Bunch
from sklearn.externals import joblib

DIRECTORY_URL = "http://www.cs.princeton.edu/~schapire/maxent/datasets/"

SAMPLES_URL = join(DIRECTORY_URL, "samples.zip")
COVERAGES_URL = join(DIRECTORY_URL, "coverages.zip")

DATA_ARCHIVE_NAME = "species_coverage.pkz"


def _load_coverage(F, header_length=6, dtype=np.int16):
    """Load a coverage file from an open file object.

    This will return a numpy array of the given dtype
    """
    header = [F.readline() for i in range(header_length)]
    make_tuple = lambda t: (t.split()[0], float(t.split()[1]))
    header = dict([make_tuple(line) for line in header])

    M = np.loadtxt(F, dtype=dtype)
    nodata = header[b'NODATA_value']
    if nodata != -9999:
        print(nodata)
        M[nodata] = -9999
    return M


def _load_csv(F):
    """Load csv file.

    Parameters
    ----------
    F : file object
        CSV file open in byte mode.

    Returns
    -------
    rec : np.ndarray
        record array representing the data
    """
    if PY2:
        # Numpy recarray wants Python 2 str but not unicode
        names = F.readline().strip().split(',')
    else:
        # Numpy recarray wants Python 3 str but not bytes...
        names = F.readline().decode('ascii').strip().split(',')
    rec = np.loadtxt(F, skiprows=0, delimiter=',', dtype='a22,f4,f4')
    rec.dtype.names = names
    return rec


def construct_grids(batch):
    """Construct the map grid from the batch object

    Parameters
    ----------
    batch : Batch object
        The object returned by :func:`fetch_species_distributions`

    Returns
    -------
    (xgrid, ygrid) : 1-D arrays
        The grid corresponding to the values in batch.coverages
    """
    # x,y coordinates for corner cells
    xmin = batch.x_left_lower_corner + batch.grid_size
    xmax = xmin + (batch.Nx * batch.grid_size)
    ymin = batch.y_left_lower_corner + batch.grid_size
    ymax = ymin + (batch.Ny * batch.grid_size)

    # x coordinates of the grid cells
    xgrid = np.arange(xmin, xmax, batch.grid_size)
    # y coordinates of the grid cells
    ygrid = np.arange(ymin, ymax, batch.grid_size)

    return (xgrid, ygrid)


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

    Returns
    --------
    The data is returned as a Bunch object with the following attributes:

    coverages : array, shape = [14, 1592, 1212]
        These represent the 14 features measured at each point of the map grid.
        The latitude/longitude values for the grid are discussed below.
        Missing data is represented by the value -9999.

    train : record array, shape = (1623,)
        The training points for the data.  Each point has three fields:

        - train['species'] is the species name
        - train['dd long'] is the longitude, in degrees
        - train['dd lat'] is the latitude, in degrees

    test : record array, shape = (619,)
        The test points for the data.  Same format as the training data.

    Nx, Ny : integers
        The number of longitudes (x) and latitudes (y) in the grid

    x_left_lower_corner, y_left_lower_corner : floats
        The (x,y) position of the lower-left corner, in degrees

    grid_size : float
        The spacing between points of the grid, in degrees

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

    Notes
    -----

    * See examples/applications/plot_species_distribution_modeling.py
      for an example of using this dataset with scikit-learn

    """
    data_home = get_data_home(data_home)
    if not exists(data_home):
        makedirs(data_home)

    # Define parameters for the data files.  These should not be changed
    # unless the data model changes.  They will be saved in the npz file
    # with the downloaded data.
    extra_params = dict(x_left_lower_corner=-94.8,
                        Nx=1212,
                        y_left_lower_corner=-56.05,
                        Ny=1592,
                        grid_size=0.05)
    dtype = np.int16

    if not exists(join(data_home, DATA_ARCHIVE_NAME)):
        print('Downloading species data from %s to %s' % (SAMPLES_URL,
                                                          data_home))
        X = np.load(BytesIO(urlopen(SAMPLES_URL).read()))

        for f in X.files:
            fhandle = BytesIO(X[f])
            if 'train' in f:
                train = _load_csv(fhandle)
            if 'test' in f:
                test = _load_csv(fhandle)

        print('Downloading coverage data from %s to %s' % (COVERAGES_URL,
                                                           data_home))

        X = np.load(BytesIO(urlopen(COVERAGES_URL).read()))

        coverages = []
        for f in X.files:
            fhandle = BytesIO(X[f])
            print(' - converting', f)
            coverages.append(_load_coverage(fhandle))
        coverages = np.asarray(coverages, dtype=dtype)

        bunch = Bunch(coverages=coverages,
                      test=test,
                      train=train,
                      **extra_params)
        joblib.dump(bunch, join(data_home, DATA_ARCHIVE_NAME), compress=9)
    else:
        bunch = joblib.load(join(data_home, DATA_ARCHIVE_NAME))

    return bunch
