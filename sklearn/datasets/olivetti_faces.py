"""Modified Olivetti faces dataset.

The original database was available from (now defunct)

    http://www.uk.research.att.com/facedatabase.html

The version retrieved here comes in MATLAB format from the personal
web page of Sam Roweis:

    http://www.cs.nyu.edu/~roweis/

There are ten different images of each of 40 distinct subjects. For some
subjects, the images were taken at different times, varying the lighting,
facial expressions (open / closed eyes, smiling / not smiling) and facial
details (glasses / no glasses). All the images were taken against a dark
homogeneous background with the subjects in an upright, frontal position (with
tolerance for some side movement).

The original dataset consisted of 92 x 112, while the Roweis version
consists of 64x64 images.
"""
# Copyright (c) 2011 David Warde-Farley <wardefar at iro dot umontreal dot ca>
# License: BSD 3 clause

from io import BytesIO
from os.path import join, exists
from os import makedirs
try:
    # Python 2
    import urllib2
    urlopen = urllib2.urlopen
except ImportError:
    # Python 3
    import urllib.request
    urlopen = urllib.request.urlopen

import numpy as np
from scipy.io.matlab import loadmat

from .base import get_data_home, Bunch
from ..utils import check_random_state
from ..externals import joblib


DATA_URL = "http://cs.nyu.edu/~roweis/data/olivettifaces.mat"
TARGET_FILENAME = "olivetti.pkz"

# Grab the module-level docstring to use as a description of the
# dataset
MODULE_DOCS = __doc__


def fetch_olivetti_faces(data_home=None, shuffle=False, random_state=0,
                         download_if_missing=True):
    """Loader for the Olivetti faces data-set from AT&T.

    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    shuffle : boolean, optional
        If True the order of the dataset is shuffled to avoid having
        images of the same person grouped.

    download_if_missing: optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    random_state : optional, integer or RandomState object
        The seed or the random number generator used to shuffle the
        data.

    Returns
    -------
    An object with the following attributes:

    data : numpy array of shape (400, 4096)
        Each row corresponds to a ravelled face image of original size 64 x 64 pixels.

    images : numpy array of shape (400, 64, 64)
        Each row is a face image corresponding to one of the 40 subjects of the dataset.

    target : numpy array of shape (400, )
        Labels associated to each face image. Those labels are ranging from 0-39 and correspond to the Subject IDs.

    DESCR : string
        Description of the modified Olivetti Faces Dataset.

    Notes
    ------

    This dataset consists of 10 pictures each of 40 individuals. The original
    database was available from (now defunct)

        http://www.uk.research.att.com/facedatabase.html

    The version retrieved here comes in MATLAB format from the personal
    web page of Sam Roweis:

        http://www.cs.nyu.edu/~roweis/

    """
    data_home = get_data_home(data_home=data_home)
    if not exists(data_home):
        makedirs(data_home)
    if not exists(join(data_home, TARGET_FILENAME)):
        print('downloading Olivetti faces from %s to %s'
              % (DATA_URL, data_home))
        fhandle = urlopen(DATA_URL)
        buf = BytesIO(fhandle.read())
        mfile = loadmat(buf)
        faces = mfile['faces'].T.copy()
        joblib.dump(faces, join(data_home, TARGET_FILENAME), compress=6)
        del mfile
    else:
        faces = joblib.load(join(data_home, TARGET_FILENAME))
    # We want floating point data, but float32 is enough (there is only
    # one byte of precision in the original uint8s anyway)
    faces = np.float32(faces)
    faces = faces - faces.min()
    faces /= faces.max()
    faces = faces.reshape((400, 64, 64)).transpose(0, 2, 1)
    # 10 images per class, 400 images total, each class is contiguous.
    target = np.array([i // 10 for i in range(400)])
    if shuffle:
        random_state = check_random_state(random_state)
        order = random_state.permutation(len(faces))
        faces = faces[order]
        target = target[order]
    return Bunch(data=faces.reshape(len(faces), -1),
                 images=faces,
                 target=target,
                 DESCR=MODULE_DOCS)
