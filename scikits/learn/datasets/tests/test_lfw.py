
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

        # to dowload the data, run the face recognition / verification examples
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
