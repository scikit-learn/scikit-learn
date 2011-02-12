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
implemented in the OpenCV library.
"""
# Copyright (c) 2011 Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

from os.path import join
from os.path import exists
from os.path import expanduser
from os import listdir
from os import makedirs

import logging

from .base import get_data_home


BASE_URL = "http://vis-www.cs.umass.edu/lfw/"
ARCHIVE_NAME = "lfw.tgz"
TARGET_FILENAMES = [
    'pairsDevTrain.txt',
    'pairsDevTest.txt',
    'pairs.txt',
    'people.txt',
]


def check_fetch_lfw(data_home=None):
    """Helper function to download any missing LFW data"""
    data_home = get_data_home(data_home=data_home)
    lfw_home = join(data_home, "lfw_home")
    archive_path = join(lfw_home, ARCHIVE_NAME)
    data_folder_path = join(lfw_home, "lfw")

    if not exists(lfw_home):
        makedirs(lfw_home)

    logging.info("LFW parent data folder: %s", lfw_home)

    if not exists(archive_path):
        import urllib
        url = BASE_URL + ARCHIVE_NAME
        logging.info("Downloading LFW data (234MB): %s", url)
        downloader = urllib.urlopen(url)
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

    return lfw_home


def load_lfw_pairs(subset='train', data_home=None):
    """Loader for the Labeled Faces in the Wild (LFW) dataset

    This dataset is a collection of JPEG pictures of famous people
    collected on the internet, all details are available on the
    official website:

        http://vis-www.cs.umass.edu/lfw/

    Each picture is centered on a single face. The task is called Face
    Verification: given a pair of two pictures, a binary classifier must
    predict whether the two images are from the same person.

    Parameters
    ----------
    subset: optional, default: 'train'
        Select the dataset to load: 'train' for the development
        training set, 'test' for the development test set and '0_train',
        '1_train'...,  '9_train' and '0_test', '1_test'..., '9_test'
        for the 10-fold cross-validation used for comparative evaluations.

    """

    lfw_home = check_fetch_lfw(data_home=data_home)


