"""Loader for the modified Olivetti faces dataset.

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
# License: Simplified BSD

from os.path import join, exists
from os import makedirs
from cStringIO import StringIO
import urllib2

import numpy as np
from scipy.io.matlab import loadmat

from .base import get_data_home
from .base import Bunch

logger = logging.getLogger(__name__)


DATA_URL = "http://cs.nyu.edu/~roweis/data/olivettifaces.mat"
TARGET_FILENAME = "olivetti.npy"

def fetch_olivetti_faces(data_home=None, download_if_missing=True):
    """Loader for the Olivetti faces data-set from AT&T.

    This dataset consists of 10 pictures each of 40 individuals. The original
    database was available from (now defunct)

        http://www.uk.research.att.com/facedatabase.html

    The version retrieved here comes in MATLAB format from the personal
    web page of Sam Roweis:

        http://www.cs.nyu.edu/~roweis/

    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing: optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.
    """
    data_home = get_data_home(data_home=data_home)
    if not exists(data_home):
        makedirs(data_home)
    if not exists(join(data_home, TARGET_FILENAME)):
        print 'Not found; downloading Olivetti faces from %s' % DATA_URL
        fhandle = urllib2.urlopen(DATA_URL)
        buf = StringIO(fhandle.read())
        mfile = loadmat(buf)
        np.save(join(data_home, TARGET_FILENAME), mfile['faces'].T)
        faces = mfile['faces'].T.copy()
        del mfile
    else:
        faces = np.load(join(data_home, TARGET_FILENAME))
    # We want floating point data, but float32 is enough (there is only
    # one byte of precision in the original uint8s anyway)
    faces = np.float32(faces)
    faces = faces - faces.min()
    faces /= faces.max()
    faces = faces.reshape((400, 64, 64)).transpose(0, 2, 1)
    # 10 images per class, 400 images total, each class is contiguous.
    target = np.array([i // 10 for i in range(400)])
    return Bunch(data=faces.reshape(len(faces), -1),
                 images=faces,
                 target=target,
                 DESCR="64x64 rescaled Olivetti faces dataset")
