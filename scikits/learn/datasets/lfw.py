"""Loader for the Labeled Faces in the Wild (LFW) dataset

This dataset is a collection of JPEG pictures of famous people collected
over the internet, all details are available on the official website:

    http://vis-www.cs.umass.edu/lfw/

Each picture is centered on a single face. The typical task is called
Face Verification: given a pair of two pictures, a binary classifier
must predict whether the two images are from the same person.

An alternative task, Face Recognition or Face Identification is:
given the picture of the face of an unknown person, identify the name
of the person by refering to a gallery of previously seen pictures of
identified persons.

Both Face Verification and Face Recognition are tasks that are typically
performed on the output of a model trained to perform Face Detection. The
most popular model for Face Detection is called Viola-Johns and is
implemented in the OpenCV library. The LFW faces were extracted by this face
detector from various online websites.
"""
# Copyright (c) 2011 Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

from os.path import join
from os.path import exists
from os.path import isdir
from os import listdir
from os import makedirs

import urllib
import logging

from scipy.misc import imread
from scipy.misc import imresize
import numpy as np

from scikits.learn.externals.joblib import Memory
from .base import get_data_home
from .base import Bunch


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
# Common data fetching from the original LFW website and local disk caching
#

def check_fetch_lfw(data_home=None, funneled=True):
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

    if not exists(archive_path):
        logging.info("Downloading LFW data: %s", archive_url)
        downloader = urllib.urlopen(archive_url)
        open(archive_path, 'wb').write(downloader.read())

    for target_filename in TARGET_FILENAMES:
        target_filepath = join(lfw_home, target_filename)
        if not exists(target_filepath):
            url = BASE_URL + target_filename
            logging.info("Downloading LFW metadata: %s", url)
            downloader = urllib.urlopen(BASE_URL + target_filename)
            open(target_filepath, 'wb').write(downloader.read())

    if not exists(data_folder_path):
        import tarfile
        logging.info("Decompressing the data archive to %s", data_folder_path)
        tarfile.open(archive_path, "r:gz").extractall(path=lfw_home)

    return lfw_home, data_folder_path


#
# Task #1:  Face Identification on picture with names
#

def _load_lfw_people(data_folder_path, slice_=None, color=False, resize=None,
                     min_faces_per_person=10):
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
    class_names = np.unique(person_names)
    target = np.searchsorted(class_names, person_names)

    # compute the portion of the images to load to respect the slice_ parameter
    # given by the caller
    default_slice = (slice(0, 250), slice(0, 250))
    if slice_ is None:
        slice_ = default_slice
    else:
        slice_ = tuple(s or ds for s, ds in zip(slice_, default_slice))

    h_slice, w_slice = slice_
    h = (h_slice.stop - h_slice.start) / (h_slice.step or 1)
    w = (w_slice.stop - w_slice.start) / (w_slice.step or 1)

    if resize is not None:
        h = int(resize * h)
        w = int(resize * w)

    # allocate some contiguous memory to host the decoded image slices
    if not color:
        faces = np.zeros((n_faces, h, w), dtype=np.float32)
    else:
        faces = np.zeros((n_faces, h, w, 3), dtype=np.float32)

    # iterate over the collected file path to load the jpeg files as numpy
    # arrays
    for i, file_path in enumerate(file_paths):
        if i % 1000 == 0:
            logging.info("Loading face #%05d / %05d", i + 1, n_faces)

        face = np.asarray(imread(file_path)[slice_], dtype=np.float32)
        face /= 255.0  # scale uint8 coded colors to the [0.0, 1.0] floats
        if resize is not None:
            face = imresize(face, resize)
        if not color:
            # average the color channels to compute a gray levels
            # representaion
            faces[i, :, :] = face.mean(axis=2)
        else:
            faces[i, :, :, :] = face

    # shuffle the faces with a deterministic RNG scheme to avoid having
    # all faces of the same person in a row, as it would break some
    # cross validation and learning algorithms such as SGD and online
    # k-means that make an IID assumption

    indices = np.arange(n_faces)
    np.random.RandomState(42).shuffle(indices)
    faces, target = faces[indices], target[indices]
    return faces, target, class_names


def load_lfw_people(data_home=None, funneled=True, resize=0.5,
                    min_faces_per_person=10, color=False,
                    slice_=(slice(70, 195), slice(78, 172))):
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

    Parameters
    ----------
    data_home: optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    funneled: boolean, optional, default: True
        Download and use the funneled variant of the dataset.

    resize: float, optional, default 0.5
        Ratio used to resize the each face picture.

    min_faces_per_person=10: int, optional, default 10
        The extracted dataset will only retain pictures of people that have at
        least `min_faces_per_person` different pictures.

    color: boolean, optional, default False
        Keep the 3 RGB channels instead of averaging them to a single
        gray level channel. If color is True the shape of the data has
        one more dimension than than the shape with color = False.

    slice_: optional
        Provide a custom 2D slice (height, width) to extract the
        'interesting' part of the jpeg files and avoid use statistical
        correlation from the background
    """

    parameters = locals().copy()
    del parameters['data_home']
    del parameters['funneled']

    lfw_home, data_folder_path = check_fetch_lfw(data_home=data_home,
                                                 funneled=funneled)
    logging.info('Loading LFW people faces from %s', lfw_home)

    # wrap the loader in a memoizing function that will return memmaped data
    # arrays for optimal memory usage
    m = Memory(cachedir=join(lfw_home, 'joblib'), mmap_mode='c', verbose=0)
    load_func = m.cache(_load_lfw_people)

    # load and memoize the pairs as np arrays
    faces, target, class_names = load_func(data_folder_path, **parameters)

    # pack the results as a Bunch instance
    return Bunch(data=faces, target=target, class_names=class_names,
                 DESCR="LFW faces dataset")


#
# Task #2:  Face Verification on pairs of face pictures
#


def _load_lfw_pairs(index_file_path, data_folder_path, slice_=None,
                    color=False, resize=None):
    """Perform the actual data loading for the LFW pairs dataset

    This operation is meant to be cached by a joblib wrapper.
    """
    # parse the index file to find the number of pairs to be able to allocate
    # the right amount of memory before starting to decode the jpeg files
    with open(index_file_path) as f:
        splitted_lines = [l.strip().split('\t') for l in f.readlines()]
    pair_specs = [l for l in splitted_lines if len(l) > 2]
    n_pairs = len(pair_specs)

    # compute the portion of the images to load to respect the slice_ parameter
    # given by the caller
    default_slice = (slice(0, 250), slice(0, 250))
    if slice_ is None:
        slice_ = default_slice
    else:
        slice_ = tuple(s or ds for s, ds in zip(slice_, default_slice))

    h_slice, w_slice = slice_
    h = (h_slice.stop - h_slice.start) / (h_slice.step or 1)
    w = (w_slice.stop - w_slice.start) / (w_slice.step or 1)

    if resize is not None:
        h = int(resize * h)
        w = int(resize * w)

    # allocate some contiguous memory to host the decoded image slices
    target = np.zeros(n_pairs, dtype=np.int)
    if not color:
        pairs = np.zeros((n_pairs, 2, h, w), dtype=np.float32)
    else:
        pairs = np.zeros((n_pairs, 2, h, w, 3), dtype=np.float32)

    # interating over the metadata lines for each pair to find the filename to
    # decode and load in memory
    for i, components in enumerate(pair_specs):
        if i % 1000 == 0:
            logging.info("Loading pair #%05d / %05d", i + 1, n_pairs)

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
            person_folder = join(data_folder_path, name)
            filenames = list(sorted(listdir(person_folder)))
            filepath = join(person_folder, filenames[idx])
            face = np.asarray(imread(filepath)[slice_], dtype=np.float32)
            face /= 255.0  # scale uint8 coded colors to the [0.0, 1.0] floats
            if resize is not None:
                face = imresize(face, resize)
            if not color:
                # average the color channels to compute a gray levels
                # representaion
                pairs[i, j, :, :] = face.mean(axis=2)
            else:
                pairs[i, j, :, :, :] = face

    return pairs, target, np.array(['Different persons', 'Same person'])


def load_lfw_pairs(subset='train', data_home=None, funneled=True, resize=0.5,
                   color=False, slice_=(slice(70, 195), slice(78, 172))):
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

    Parameters
    ----------
    subset: optional, default: 'train'
        Select the dataset to load: 'train' for the development training
        set, 'test' for the development test set, and '10_folds' for the
        official evaluation set that is meant to be used with a 10-folds
        cross validation.

    data_home: optional, default: None
        Specify another download and cache folder for the datasets. By
        default all scikit learn data is stored in '~/scikit_learn_data'
        subfolders.

    funneled: boolean, optional, default: True
        Download and use the funneled variant of the dataset.

    resize: float, optional, default 0.5
        Ratio used to resize the each face picture.

    color: boolean, optional, default False
        Keep the 3 RGB channels instead of averaging them to a single
        gray level channel. If color is True the shape of the data has
        one more dimension than than the shape with color = False.

    slice_: optional
        Provide a custom 2D slice (height, width) to extract the
        'interesting' part of the jpeg files and avoid use statistical
        correlation from the background
    """

    parameters = locals().copy()
    del parameters['subset']
    del parameters['data_home']
    del parameters['funneled']

    lfw_home, data_folder_path = check_fetch_lfw(data_home=data_home,
                                                 funneled=funneled)
    logging.info('Loading %s LFW pairs from %s', subset, lfw_home)

    # wrap the loader in a memoizing function that will return memmaped data
    # arrays for optimal memory usage
    m = Memory(cachedir=join(lfw_home, 'joblib'), mmap_mode='c', verbose=0)
    load_func = m.cache(_load_lfw_pairs)

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
    pairs, target, class_names = load_func(index_file_path, data_folder_path,
                                           **parameters)

    # pack the results as a Bunch instance
    return Bunch(data=pairs, target=target, class_names=class_names,
                 DESCR="'%s' segment of the LFW pairs dataset" % subset)
