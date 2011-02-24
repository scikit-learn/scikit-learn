"""This test for the LFW require medium-size data dowloading and processing

If the data has not been already downloaded by runnning the examples, the tests
won't run (skipped).

If the test are run, the first execution will be long (typically a bit more than
a couple of minutes) but as the dataset loader is leveraging joblib, successive
runs will be fast (less than 200ms).
"""

import os
from nose import SkipTest

from scikits.learn.datasets import load_lfw_pairs
from scikits.learn.datasets import load_lfw_people
from scikits.learn.datasets import get_data_home

from numpy.testing import assert_array_equal
from numpy.testing import assert_equal


def test_load_lfw_people():
    if not os.path.exists(os.path.join(get_data_home(), 'lfw_home')):
        # skip this test is the data has not already been previously
        # downloaded to avoid having tests rely on the availability of a
        # fast internet connection

        # to download the data, run the face recognition / verification examples
        # or call load_lfw_people function from an interactive shell for
        # instance
        raise SkipTest

    lfw_people = load_lfw_people(min_faces_per_person=100)

    # only 5 person have more than 100 pictures each in the dataset
    top_classes = ['Colin Powell', 'Donald Rumsfeld', 'George W Bush',
                   'Gerhard Schroeder', 'Tony Blair']
    assert_array_equal(lfw_people.class_names, top_classes)

    # default slice is a rectangular shape around the face, removing
    # most of the background
    assert_equal(lfw_people.data.shape, (1140, 62, 47))

    # person ids have been shuffled to avoid having the photo ordered by
    # alphabetical ordering as in the default tarball layout
    assert_equal(lfw_people.target.shape, (1140,))
    assert_array_equal(lfw_people.target[:5], [2, 3, 1, 4, 1])

    # it is possible to slice the data in different ways and to resize the
    # outpout without changing the width / heigh ratio
    lfw_people = load_lfw_people(min_faces_per_person=100,
                                 slice_=(slice(50, 200), slice(50, 200)),
                                 resize=0.1)
    assert_equal(lfw_people.data.shape, (1140, 15, 15))

    # it is also possible to load the color version of the data, in that
    # case the color channels are stored in the last dimension of the data
    lfw_people = load_lfw_people(min_faces_per_person=100, color=True)
    assert_equal(lfw_people.data.shape, (1140, 62, 47, 3))


def test_load_lfw_pairs():
    if not os.path.exists(os.path.join(get_data_home(), 'lfw_home')):
        raise SkipTest

    lfw_pairs_train = load_lfw_pairs(subset='train')

    # this dataset is used for training supervised face verification models,
    # this is a binary classification task
    top_classes = ['Different persons', 'Same person']
    assert_array_equal(lfw_pairs_train.class_names, top_classes)

    # default slice is a rectangular shape around the face, removing
    # most of the background, for each of the 2 face pictures
    assert_equal(lfw_pairs_train.data.shape, (2200, 2, 62, 47))

    # the ordering is respecting the metadata text file of the official LFW
    # tasks
    assert_equal(lfw_pairs_train.target.shape, (2200,))
    assert_array_equal(lfw_pairs_train.target[:5], [1, 1, 1, 1, 1])
    assert_array_equal(lfw_pairs_train.target[-5:], [0, 0, 0, 0, 0])

    # as for the people loader it is also possible to load the color channels
    # in the last dimension
    lfw_pairs_train = load_lfw_pairs(subset='train', color=True)
    assert_equal(lfw_pairs_train.data.shape, (2200, 2, 62, 47, 3))

    # the data also has a test development set and a 10-fold CV dataset for
    # final evaluation
    lfw_pairs_test = load_lfw_pairs(subset='test')
    assert_equal(lfw_pairs_test.data.shape, (1000, 2, 62, 47))

    lfw_pairs_10_folds = load_lfw_pairs(subset='10_folds')
    assert_equal(lfw_pairs_10_folds.data.shape, (6000, 2, 62, 47))

