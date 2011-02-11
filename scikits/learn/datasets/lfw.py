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

from .base import get_data_dir

def load_lfw_pairs(subset='train', data_dir=None):
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

    data_dir = get_data_dir(data_dir=data_dir)
    base_url = "http://vis-www.cs.umass.edu/lfw/"
    archive_name = "lfw.tgz"
    archive_path = join(data_dir, archive_name)
    folder_path = join(data_dir, "lfw")
    target_filenames = ['pairsDevTrain.txt', 'pairsDevTest.txt',
                        'pairs.txt', 'people.txt']

    if exists(folder_path):
        os.makedirs(folder_path)
        if not exists(archive_path):
            import urllib
            url = base_url + archive_name
            print "Downloading LFW data, please wait (234MB)..."
            print url
            downloader = urllib.urlopen(url)
            open(archive_path, 'wb').write(downloader.read())
            print

        for target_filename in target_filename:
            if not exists(target_filepath):
                url = base_url + target_filepath
                print "Downloading " + url
                target_filepath = join(folder_path, target_filepath)
                downloader = urllib.urlopen(base_url + target_filename)
                open(target_filepath, 'wb').write(downloader.read())

        import tarfile
        print "Decompressiong the archive: " + archive_path
        tarfile.open(archive_path, "r:gz").extractall()
        print


