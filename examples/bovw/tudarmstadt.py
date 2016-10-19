"""TU Darmstadt dataset.
The original database was available from
    http://host.robots.ox.ac.uk/pascal/VOC/download/tud.tar.gz
"""

import os
from os.path import join, exists
try:
    # Python 2
    from urllib2 import HTTPError
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.error import HTTPError
    from urllib.request import urlopen

import numpy as np
import tarfile

from sklearn.datasets import get_data_home

DATA_URL = "http://host.robots.ox.ac.uk/pascal/VOC/download/tud.tar.gz"
TARGET_FILENAME = "tud.pkz"

# Grab the module-level docstring to use as a description of the
# dataset
MODULE_DOCS = __doc__


def fetch_tu_darmstadt(data_home=None):
    """Loader for the TU Darmstadt dataset.

    Read more in the :ref:`User Guide <datasets>`.


    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    Returns
    -------
    images_list : list
        Python list with the path of each image to consider during the
        classification.

    labels : array-like, shape (n_images, )
        An array with the different label corresponding to the categories.
        0: motorbikes - 1: cars - 2: cows.

    Notes
    ------
    The dataset is composed of 124 motorbikes images, 100 cars, and 112 cows.

    Examples
    --------
    Load the 'tu-darmstadt' dataset:

    >>> from sklearn.datasets.tudarmstadt import fetch_tu_darmstadt
    >>> import tempfile
    >>> test_data_home = tempfile.mkdtemp()
    >>> im_list, labels = fetch_tu_darmstadt(data_home=test_data_home)
    """

    # check if the data has been already downloaded
    data_home = get_data_home(data_home=data_home)
    data_home = join(data_home, 'tu_darmstadt')
    if not exists(data_home):
        os.makedirs(data_home)

    # dataset tar file
    filename = join(data_home, 'tud.tar.gz')

    # if the file does not exist, download it
    if not exists(filename):
        try:
            db_url = urlopen(DATA_URL)
            with open(filename, 'wb') as f:
                f.write(db_url.read())
            db_url.close()
        except HTTPError as e:
            if e.code == 404:
                e.msg = 'TU Darmstadt dataset not found.'
            raise
    # Try to extract the complete archieve
    try:
        tarfile.open(filename, "r:gz").extractall(path=data_home)
    except:
        os.remove(filename)
        raise

    # the file 'motorbikes023' is a gray scale image and need to be removed
    file_removal = [
        join(data_home,
             'TUDarmstadt/PNGImages/motorbike-testset/motorbikes023.png'),
        join(data_home,
             'TUDarmstadt/Annotations/motorbike-testset/motorbikes023.txt'),
    ]
    for f in file_removal:
        os.remove(f)

    # list the different images
    data_path = join(data_home, 'TUDarmstadt/PNGImages')
    images_list = [os.path.join(root, name)
                   for root, dirs, files in os.walk(data_path)
                   for name in files
                   if name.endswith((".png"))]

    # create the label array
    labels = []
    for imf in images_list:
        if 'motorbike' in imf:
            labels.append(0)
        elif 'cars' in imf:
            labels.append(1)
        elif 'cows' in imf:
            labels.append(2)

    # Return these information
    return images_list, np.array(labels)
