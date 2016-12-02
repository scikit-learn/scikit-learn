"""Loader for the Labeled Faces in the Wild (LFW) dataset

This dataset is a collection of JPEG pictures of famous people collected
over the internet, all details are available on the official website:

    http://vis-www.cs.umass.edu/lfw/

Each picture is centered on a single face. The typical task is called
Face Verification: given a pair of two pictures, a binary classifier
must predict whether the two images are from the same person.

An alternative task, Face Recognition or Face Identification is:
given the picture of the face of an unknown person, identify the name
of the person by referring to a gallery of previously seen pictures of
identified persons.

Both Face Verification and Face Recognition are tasks that are typically
performed on the output of a model trained to perform Face Detection. The
most popular model for Face Detection is called Viola-Johns and is
implemented in the OpenCV library. The LFW faces were extracted by this face
detector from various online websites.
"""
# Copyright (c) 2011 Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

from os import listdir, makedirs, remove
from os.path import join, exists, isdir

from sklearn.utils import deprecated

import logging
import numpy as np

try:
    import urllib.request as urllib  # for backwards compatibility
except ImportError:
    import urllib

from .base import get_data_home, Bunch
from ..externals.joblib import Memory

from ..externals.six import b

logger = logging.getLogger(__name__)


BASE_URL = "http://vis-www.cs.umass.edu/lfw/"
ARCHIVE_NAME = "lfw.tgz"
FUNNELED_ARCHIVE_NAME = "lfw-funneled.tgz"
TARGET_FILENAMES = [
    'pairsDevTrain.txt',
    'pairsDevTest.txt',
    'pairs.txt',
]


def scale_face(face):
    """Scale back to 0-1 range in case of normalization for plotting"""
    scaled = face - face.min()
    scaled /= scaled.max()
    return scaled


#
# Common private utilities for data fetching from the original LFW website
# local disk caching, and image decoding.
#


def check_fetch_lfw(data_home=None, funneled=True, download_if_missing=True):
    """Helper function to download any missing LFW data"""
    data_home = get_data_home(data_home=data_home)
    lfw_home = join(data_home, "lfw_home")

    if funneled:
        archive_path = join(lfw_home, FUNNELED_ARCHIVE_NAME)
        data_folder_path = join(lfw_home, "lfw_funneled")
        archive_url = BASE_URL + FUNNELED_ARCHIVE_NAME
    else:
        archive_path = join(lfw_home, ARCHIVE_NAME)
        data_folder_path = join(lfw_home, "lfw")
        archive_url = BASE_URL + ARCHIVE_NAME

    if not exists(lfw_home):
        makedirs(lfw_home)

    for target_filename in TARGET_FILENAMES:
        target_filepath = join(lfw_home, target_filename)
        if not exists(target_filepath):
            if download_if_missing:
                url = BASE_URL + target_filename
                logger.warning("Downloading LFW metadata: %s", url)
                urllib.urlretrieve(url, target_filepath)
            else:
                raise IOError("%s is missing" % target_filepath)

    if not exists(data_folder_path):

        if not exists(archive_path):
            if download_if_missing:
                logger.warning("Downloading LFW data (~200MB): %s",
                               archive_url)
                urllib.urlretrieve(archive_url, archive_path)
            else:
                raise IOError("%s is missing" % target_filepath)

        import tarfile
        logger.info("Decompressing the data archive to %s", data_folder_path)
        tarfile.open(archive_path, "r:gz").extractall(path=lfw_home)
        remove(archive_path)

    return lfw_home, data_folder_path


def _load_imgs(file_paths, slice_, color, resize):
    """Internally used to load images"""

    # Try to import imread and imresize from PIL. We do this here to prevent
    # the whole sklearn.datasets module from depending on PIL.
    try:
        try:
            from scipy.misc import imread
        except ImportError:
            from scipy.misc.pilutil import imread
        from scipy.misc import imresize
    except ImportError:
        raise ImportError("The Python Imaging Library (PIL)"
                          " is required to load data from jpeg files")

    # compute the portion of the images to load to respect the slice_ parameter
    # given by the caller
    default_slice = (slice(0, 250), slice(0, 250))
    if slice_ is None:
        slice_ = default_slice
    else:
        slice_ = tuple(s or ds for s, ds in zip(slice_, default_slice))

    h_slice, w_slice = slice_
    h = (h_slice.stop - h_slice.start) // (h_slice.step or 1)
    w = (w_slice.stop - w_slice.start) // (w_slice.step or 1)

    if resize is not None:
        resize = float(resize)
        h = int(resize * h)
        w = int(resize * w)

    # allocate some contiguous memory to host the decoded image slices
    n_faces = len(file_paths)
    if not color:
        faces = np.zeros((n_faces, h, w), dtype=np.float32)
    else:
        faces = np.zeros((n_faces, h, w, 3), dtype=np.float32)

    # iterate over the collected file path to load the jpeg files as numpy
    # arrays
    for i, file_path in enumerate(file_paths):
        if i % 1000 == 0:
            logger.info("Loading face #%05d / %05d", i + 1, n_faces)

        # Checks if jpeg reading worked. Refer to issue #3594 for more
        # details.
        img = imread(file_path)
        if img.ndim is 0:
            raise RuntimeError("Failed to read the image file %s, "
                               "Please make sure that libjpeg is installed"
                               % file_path)

        face = np.asarray(img[slice_], dtype=np.float32)
        face /= 255.0  # scale uint8 coded colors to the [0.0, 1.0] floats
        if resize is not None:
            face = imresize(face, resize)
        if not color:
            # average the color channels to compute a gray levels
            # representation
            face = face.mean(axis=2)

        faces[i, ...] = face

    return faces


#
# Task #1:  Face Identification on picture with names
#

def _fetch_lfw_people(data_folder_path, slice_=None, color=False, resize=None,
                      min_faces_per_person=0):
    """Perform the actual data loading for the lfw people dataset

    This operation is meant to be cached by a joblib wrapper.
    """
    # scan the data folder content to retain people with more that
    # `min_faces_per_person` face pictures
    person_names, file_paths = [], []
    for person_name in sorted(listdir(data_folder_path)):
        folder_path = join(data_folder_path, person_name)
        if not isdir(folder_path):
            continue
        paths = [join(folder_path, f) for f in listdir(folder_path)]
        n_pictures = len(paths)
        if n_pictures >= min_faces_per_person:
            person_name = person_name.replace('_', ' ')
            person_names.extend([person_name] * n_pictures)
            file_paths.extend(paths)

    n_faces = len(file_paths)
    if n_faces == 0:
        raise ValueError("min_faces_per_person=%d is too restrictive" %
                         min_faces_per_person)

    target_names = np.unique(person_names)
    target = np.searchsorted(target_names, person_names)

    faces = _load_imgs(file_paths, slice_, color, resize)

    # shuffle the faces with a deterministic RNG scheme to avoid having
    # all faces of the same person in a row, as it would break some
    # cross validation and learning algorithms such as SGD and online
    # k-means that make an IID assumption

    indices = np.arange(n_faces)
    np.random.RandomState(42).shuffle(indices)
    faces, target = faces[indices], target[indices]
    return faces, target, target_names


def fetch_lfw_people(data_home=None, funneled=True, resize=0.5,
                     min_faces_per_person=0, color=False,
                     slice_=(slice(70, 195), slice(78, 172)),
                     download_if_missing=True):
    """Loader for the Labeled Faces in the Wild (LFW) people dataset

    This dataset is a collection of JPEG pictures of famous people
    collected on the internet, all details are available on the
    official website:

        http://vis-www.cs.umass.edu/lfw/

    Each picture is centered on a single face. Each pixel of each channel
    (color in RGB) is encoded by a float in range 0.0 - 1.0.

    The task is called Face Recognition (or Identification): given the
    picture of a face, find the name of the person given a training set
    (gallery).

    The original images are 250 x 250 pixels, but the default slice and resize
    arguments reduce them to 62 x 74.

    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    funneled : boolean, optional, default: True
        Download and use the funneled variant of the dataset.

    resize : float, optional, default 0.5
        Ratio used to resize the each face picture.

    min_faces_per_person : int, optional, default None
        The extracted dataset will only retain pictures of people that have at
        least `min_faces_per_person` different pictures.

    color : boolean, optional, default False
        Keep the 3 RGB channels instead of averaging them to a single
        gray level channel. If color is True the shape of the data has
        one more dimension than than the shape with color = False.

    slice_ : optional
        Provide a custom 2D slice (height, width) to extract the
        'interesting' part of the jpeg files and avoid use statistical
        correlation from the background

    download_if_missing : optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    Returns
    -------
    dataset : dict-like object with the following attributes:

    dataset.data : numpy array of shape (13233, 2914)
        Each row corresponds to a ravelled face image of original size 62 x 47
        pixels. Changing the ``slice_`` or resize parameters will change the
        shape of the output.

    dataset.images : numpy array of shape (13233, 62, 47)
        Each row is a face image corresponding to one of the 5749 people in
        the dataset. Changing the ``slice_`` or resize parameters will change
        the shape of the output.

    dataset.target : numpy array of shape (13233,)
        Labels associated to each face image. Those labels range from 0-5748
        and correspond to the person IDs.

    dataset.DESCR : string
        Description of the Labeled Faces in the Wild (LFW) dataset.
    """
    lfw_home, data_folder_path = check_fetch_lfw(
        data_home=data_home, funneled=funneled,
        download_if_missing=download_if_missing)
    logger.info('Loading LFW people faces from %s', lfw_home)

    # wrap the loader in a memoizing function that will return memmaped data
    # arrays for optimal memory usage
    m = Memory(cachedir=lfw_home, compress=6, verbose=0)
    load_func = m.cache(_fetch_lfw_people)

    # load and memoize the pairs as np arrays
    faces, target, target_names = load_func(
        data_folder_path, resize=resize,
        min_faces_per_person=min_faces_per_person, color=color, slice_=slice_)

    # pack the results as a Bunch instance
    return Bunch(data=faces.reshape(len(faces), -1), images=faces,
                 target=target, target_names=target_names,
                 DESCR="LFW faces dataset")


#
# Task #2:  Face Verification on pairs of face pictures
#


def _fetch_lfw_pairs(index_file_path, data_folder_path, slice_=None,
                     color=False, resize=None):
    """Perform the actual data loading for the LFW pairs dataset

    This operation is meant to be cached by a joblib wrapper.
    """
    # parse the index file to find the number of pairs to be able to allocate
    # the right amount of memory before starting to decode the jpeg files
    with open(index_file_path, 'rb') as index_file:
        split_lines = [ln.strip().split(b('\t')) for ln in index_file]
    pair_specs = [sl for sl in split_lines if len(sl) > 2]
    n_pairs = len(pair_specs)

    # iterating over the metadata lines for each pair to find the filename to
    # decode and load in memory
    target = np.zeros(n_pairs, dtype=np.int)
    file_paths = list()
    for i, components in enumerate(pair_specs):
        if len(components) == 3:
            target[i] = 1
            pair = (
                (components[0], int(components[1]) - 1),
                (components[0], int(components[2]) - 1),
            )
        elif len(components) == 4:
            target[i] = 0
            pair = (
                (components[0], int(components[1]) - 1),
                (components[2], int(components[3]) - 1),
            )
        else:
            raise ValueError("invalid line %d: %r" % (i + 1, components))
        for j, (name, idx) in enumerate(pair):
            try:
                person_folder = join(data_folder_path, name)
            except TypeError:
                person_folder = join(data_folder_path, str(name, 'UTF-8'))
            filenames = list(sorted(listdir(person_folder)))
            file_path = join(person_folder, filenames[idx])
            file_paths.append(file_path)

    pairs = _load_imgs(file_paths, slice_, color, resize)
    shape = list(pairs.shape)
    n_faces = shape.pop(0)
    shape.insert(0, 2)
    shape.insert(0, n_faces // 2)
    pairs.shape = shape

    return pairs, target, np.array(['Different persons', 'Same person'])


@deprecated("Function 'load_lfw_people' has been deprecated in 0.17 and will "
            "be removed in 0.19."
            "Use fetch_lfw_people(download_if_missing=False) instead.")
def load_lfw_people(download_if_missing=False, **kwargs):
    """Alias for fetch_lfw_people(download_if_missing=False)

    Check fetch_lfw_people.__doc__ for the documentation and parameter list.
    """
    return fetch_lfw_people(download_if_missing=download_if_missing, **kwargs)


def fetch_lfw_pairs(subset='train', data_home=None, funneled=True, resize=0.5,
                    color=False, slice_=(slice(70, 195), slice(78, 172)),
                    download_if_missing=True):
    """Loader for the Labeled Faces in the Wild (LFW) pairs dataset

    This dataset is a collection of JPEG pictures of famous people
    collected on the internet, all details are available on the
    official website:

        http://vis-www.cs.umass.edu/lfw/

    Each picture is centered on a single face. Each pixel of each channel
    (color in RGB) is encoded by a float in range 0.0 - 1.0.

    The task is called Face Verification: given a pair of two pictures,
    a binary classifier must predict whether the two images are from
    the same person.

    In the official `README.txt`_ this task is described as the
    "Restricted" task.  As I am not sure as to implement the
    "Unrestricted" variant correctly, I left it as unsupported for now.

      .. _`README.txt`: http://vis-www.cs.umass.edu/lfw/README.txt

    The original images are 250 x 250 pixels, but the default slice and resize
    arguments reduce them to 62 x 74.

    Read more in the :ref:`User Guide <labeled_faces_in_the_wild>`.

    Parameters
    ----------
    subset : optional, default: 'train'
        Select the dataset to load: 'train' for the development training
        set, 'test' for the development test set, and '10_folds' for the
        official evaluation set that is meant to be used with a 10-folds
        cross validation.

    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By
        default all scikit learn data is stored in '~/scikit_learn_data'
        subfolders.

    funneled : boolean, optional, default: True
        Download and use the funneled variant of the dataset.

    resize : float, optional, default 0.5
        Ratio used to resize the each face picture.

    color : boolean, optional, default False
        Keep the 3 RGB channels instead of averaging them to a single
        gray level channel. If color is True the shape of the data has
        one more dimension than than the shape with color = False.

    slice_ : optional
        Provide a custom 2D slice (height, width) to extract the
        'interesting' part of the jpeg files and avoid use statistical
        correlation from the background

    download_if_missing : optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    Returns
    -------
    The data is returned as a Bunch object with the following attributes:

    data : numpy array of shape (2200, 5828). Shape depends on ``subset``.
        Each row corresponds to 2 ravel'd face images of original size 62 x 47
        pixels. Changing the ``slice_``, ``resize`` or ``subset`` parameters
        will change the shape of the output.

    pairs : numpy array of shape (2200, 2, 62, 47). Shape depends on
            ``subset``.
        Each row has 2 face images corresponding to same or different person
        from the dataset containing 5749 people. Changing the ``slice_``,
        ``resize`` or ``subset`` parameters will change the shape of the
        output.

    target : numpy array of shape (2200,). Shape depends on ``subset``.
        Labels associated to each pair of images. The two label values being
        different persons or the same person.

    DESCR : string
        Description of the Labeled Faces in the Wild (LFW) dataset.

    """
    lfw_home, data_folder_path = check_fetch_lfw(
        data_home=data_home, funneled=funneled,
        download_if_missing=download_if_missing)
    logger.info('Loading %s LFW pairs from %s', subset, lfw_home)

    # wrap the loader in a memoizing function that will return memmaped data
    # arrays for optimal memory usage
    m = Memory(cachedir=lfw_home, compress=6, verbose=0)
    load_func = m.cache(_fetch_lfw_pairs)

    # select the right metadata file according to the requested subset
    label_filenames = {
        'train': 'pairsDevTrain.txt',
        'test': 'pairsDevTest.txt',
        '10_folds': 'pairs.txt',
    }
    if subset not in label_filenames:
        raise ValueError("subset='%s' is invalid: should be one of %r" % (
            subset, list(sorted(label_filenames.keys()))))
    index_file_path = join(lfw_home, label_filenames[subset])

    # load and memoize the pairs as np arrays
    pairs, target, target_names = load_func(
        index_file_path, data_folder_path, resize=resize, color=color,
        slice_=slice_)

    # pack the results as a Bunch instance
    return Bunch(data=pairs.reshape(len(pairs), -1), pairs=pairs,
                 target=target, target_names=target_names,
                 DESCR="'%s' segment of the LFW pairs dataset" % subset)


@deprecated("Function 'load_lfw_pairs' has been deprecated in 0.17 and will "
            "be removed in 0.19."
            "Use fetch_lfw_pairs(download_if_missing=False) instead.")
def load_lfw_pairs(download_if_missing=False, **kwargs):
    """Alias for fetch_lfw_pairs(download_if_missing=False)

    Check fetch_lfw_pairs.__doc__ for the documentation and parameter list.
    """
    return fetch_lfw_pairs(download_if_missing=download_if_missing, **kwargs)
