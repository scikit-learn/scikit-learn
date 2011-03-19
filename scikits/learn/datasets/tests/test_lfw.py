"""This test for the LFW require medium-size data dowloading and processing

If the data has not been already downloaded by runnning the examples,
the tests won't run (skipped).

If the test are run, the first execution will be long (typically a bit
more than a couple of minutes) but as the dataset loader is leveraging
joblib, successive runs will be fast (less than 200ms).
"""

import random
import os
import shutil
import tempfile
import numpy as np
from scipy.misc import imsave

from scikits.learn.datasets import fetch_lfw_pairs
from scikits.learn.datasets import fetch_lfw_people
from scikits.learn.datasets import get_data_home

from numpy.testing import assert_array_equal
from numpy.testing import assert_equal
from nose import SkipTest
from nose.tools import raises


SCIKIT_LEARN_DATA = tempfile.mkdtemp(prefix="scikit_learn_lfw_test_")
LFW_HOME = os.path.join(SCIKIT_LEARN_DATA, 'lfw_home')
FAKE_NAMES = [
    'Abdelatif_Smith',
    'Abhati_Kepler',
    'Camara_Alvaro',
    'Chen_Dupont',
    'John_Lee',
    'Lin_Bauman',
    'Onur_Lopez',
]


def setup_module():
    """Test fixture run once and common to all tests of this module"""
    if not os.path.exists(LFW_HOME):
        os.makedirs(LFW_HOME)

    rng = random.Random(42)
    np_rng = np.random.RandomState(42)

    # generate some random jpeg files for each person
    counts = {}
    for name in FAKE_NAMES:
        folder_name = os.path.join(LFW_HOME, 'lfw_funneled', name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        n_faces = np_rng.randint(1, 5)
        counts[name] = n_faces
        for i in range(n_faces):
            file_path = os.path.join(folder_name, name + '_%04d.jpg' % i)
            uniface = np_rng.randint(0, 255, size=(250, 250, 3))
            imsave(file_path, uniface)

    # add some random file pollution to test robustness
    f = open(os.path.join(LFW_HOME, 'lfw_funneled', '.test.swp'), 'wb')
    f.write('Text file to be ignored by the dataset loader.')
    f.close()

    # generate some pairing metadata files using the same format as LFW
    f = open(os.path.join(LFW_HOME, 'pairsDevTrain.txt'), 'wb')
    f.write("10\n")
    more_than_two = [name for name, count in counts.iteritems()
                     if count >= 2]
    for i in range(5):
        name = rng.choice(more_than_two)
        first, second = rng.sample(range(counts[name]), 2)
        f.write('%s\t%d\t%d\n' % (name, first, second))

    for i in range(5):
        first_name, second_name = rng.sample(FAKE_NAMES, 2)
        first_index = rng.choice(range(counts[first_name]))
        second_index = rng.choice(range(counts[second_name]))
        f.write('%s\t%d\t%s\t%d\n' % (first_name, first_index,
                                      second_name, second_index))
    f.close()

    f = open(os.path.join(LFW_HOME, 'pairsDevTest.txt'), 'wb')
    f.write("Fake place holder that won't be tested")
    f.close()

    f = open(os.path.join(LFW_HOME, 'pairs.txt'), 'wb')
    f.write("Fake place holder that won't be tested")
    f.close()


def teardown_module():
    """Test fixture (clean up) run once after all tests of this module"""
    if os.path.isdir(SCIKIT_LEARN_DATA):
        shutil.rmtree(SCIKIT_LEARN_DATA)


def test_load_fake_lfw_people():
    lfw_people = fetch_lfw_people(data_home=SCIKIT_LEARN_DATA,
                                 min_faces_per_person=3)

    # The data is croped around the center as a rectangular bounding box
    # arounthe the face. Colors are converted to gray levels:
    assert_equal(lfw_people.data.shape, (10, 62, 47))

    # the target is array of person integer ids
    assert_array_equal(lfw_people.target, [2, 0, 1, 0, 2, 0, 2, 1, 1, 2])

    # names of the persons can be found using the target_names array
    expected_classes = ['Abdelatif Smith', 'Abhati Kepler', 'Onur Lopez']
    assert_array_equal(lfw_people.target_names, expected_classes)

    # It is possible to ask for the original data without any croping or color
    # conversion
    lfw_people = fetch_lfw_people(data_home=SCIKIT_LEARN_DATA,
                                 min_faces_per_person=3,
                                 resize=None, slice_=None, color=True)
    assert_equal(lfw_people.data.shape, (10, 250, 250, 3))

    # the ids and class names are the same as previously
    assert_array_equal(lfw_people.target, [2, 0, 1, 0, 2, 0, 2, 1, 1, 2])
    assert_array_equal(lfw_people.target_names, expected_classes)


@raises(ValueError)
def test_load_fake_lfw_people_too_restrictive():
    fetch_lfw_people(data_home=SCIKIT_LEARN_DATA, min_faces_per_person=100)


def test_load_fake_lfw_pairs():
    lfw_pairs_train = fetch_lfw_pairs(data_home=SCIKIT_LEARN_DATA)

    # The data is croped around the center as a rectangular bounding box
    # arounthe the face. Colors are converted to gray levels:
    assert_equal(lfw_pairs_train.data.shape, (10, 2, 62, 47))

    # the target is whether the person is the same or not
    assert_array_equal(lfw_pairs_train.target, [1, 1, 1, 1, 1, 0, 0, 0, 0, 0])

    # names of the persons can be found using the target_names array
    expected_classes = ['Different persons', 'Same person']
    assert_array_equal(lfw_pairs_train.target_names, expected_classes)

    # It is possible to ask for the original data without any croping or color
    # conversion
    lfw_pairs_train = fetch_lfw_pairs(data_home=SCIKIT_LEARN_DATA,
                                     resize=None, slice_=None, color=True)
    assert_equal(lfw_pairs_train.data.shape, (10, 2, 250, 250, 3))

    # the ids and class names are the same as previously
    assert_array_equal(lfw_pairs_train.target, [1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
    assert_array_equal(lfw_pairs_train.target_names, expected_classes)


def test_fetch_lfw_people():
    if not os.path.exists(os.path.join(get_data_home(), 'lfw_home')):
        # skip this test is the data has not already been previously
        # downloaded to avoid having tests rely on the availability of a
        # fast internet connection

        # to download the data, run the face recognition / verification
        # examples or call fetch_lfw_people function from an interactive shell
        # for instance
        raise SkipTest

    lfw_people = fetch_lfw_people(min_faces_per_person=100)

    # only 5 person have more than 100 pictures each in the dataset
    top_classes = ['Colin Powell', 'Donald Rumsfeld', 'George W Bush',
                   'Gerhard Schroeder', 'Tony Blair']
    assert_array_equal(lfw_people.target_names, top_classes)

    # default slice is a rectangular shape around the face, removing
    # most of the background
    assert_equal(lfw_people.data.shape, (1140, 62, 47))

    # person ids have been shuffled to avoid having the photo ordered by
    # alphabetical ordering as in the default tarball layout
    assert_equal(lfw_people.target.shape, (1140,))
    assert_array_equal(lfw_people.target[:5], [2, 3, 1, 4, 1])

    # it is possible to slice the data in different ways and to resize the
    # outpout without changing the width / heigh ratio
    lfw_people = fetch_lfw_people(min_faces_per_person=100,
                                 slice_=(slice(50, 200), slice(50, 200)),
                                 resize=0.1)
    assert_equal(lfw_people.data.shape, (1140, 15, 15))

    # it is also possible to load the color version of the data, in that
    # case the color channels are stored in the last dimension of the data
    lfw_people = fetch_lfw_people(min_faces_per_person=100, color=True)
    assert_equal(lfw_people.data.shape, (1140, 62, 47, 3))


def test_fetch_lfw_pairs():
    if not os.path.exists(os.path.join(get_data_home(), 'lfw_home')):
        raise SkipTest

    lfw_pairs_train = fetch_lfw_pairs(subset='train')

    # this dataset is used for training supervised face verification models,
    # this is a binary classification task
    top_classes = ['Different persons', 'Same person']
    assert_array_equal(lfw_pairs_train.target_names, top_classes)

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
    lfw_pairs_train = fetch_lfw_pairs(subset='train', color=True)
    assert_equal(lfw_pairs_train.data.shape, (2200, 2, 62, 47, 3))

    # the data also has a test development set and a 10-fold CV dataset for
    # final evaluation
    lfw_pairs_test = fetch_lfw_pairs(subset='test')
    assert_equal(lfw_pairs_test.data.shape, (1000, 2, 62, 47))

    lfw_pairs_10_folds = fetch_lfw_pairs(subset='10_folds')
    assert_equal(lfw_pairs_10_folds.data.shape, (6000, 2, 62, 47))
