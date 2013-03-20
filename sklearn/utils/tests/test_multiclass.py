import numpy as np

from sklearn.externals.six.moves import xrange

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises

from sklearn.utils.multiclass import unique_labels
from sklearn.utils.multiclass import is_label_indicator_matrix
from sklearn.utils.multiclass import is_multilabel


def test_unique_labels():
    # Empty iterable
    assert_raises(ValueError, unique_labels)

    # Multiclass problem
    assert_array_equal(unique_labels(xrange(10)), np.arange(10))
    assert_array_equal(unique_labels(np.arange(10)), np.arange(10))
    assert_array_equal(unique_labels([4, 0, 2]), np.array([0, 2, 4]))

    # Multilabels
    assert_array_equal(unique_labels([(0, 1, 2), (0,), tuple(), (2, 1)]),
                       np.arange(3))
    assert_array_equal(unique_labels([[0, 1, 2], [0], list(), [2, 1]]),
                       np.arange(3))
    assert_array_equal(unique_labels(np.array([[0, 0, 1],
                                               [1, 0, 1],
                                               [0, 0, 0]])),
                       np.arange(3))

    # Several arrays passed
    assert_array_equal(unique_labels([4, 0, 2], xrange(5)),
                       np.arange(5))
    assert_array_equal(unique_labels((0, 1, 2), (0,), (2, 1)),
                       np.arange(3))


def test_is_multilabel():
    assert_true(is_multilabel([[1], [2], [0, 1]]))
    assert_true(is_multilabel([[1], [2]]))
    assert_true(is_multilabel([[1], [2], []]))
    assert_true(is_multilabel([[1], [0, 2], []]))
    assert_true(is_multilabel(np.random.randint(2, size=(10, 10))))

    assert_false(is_multilabel(range(10)))
    assert_false(is_multilabel(np.arange(10)))
    assert_false(is_multilabel(np.reshape(np.arange(10), (1, -1))))
    assert_false(is_multilabel(np.reshape(np.arange(10), (-1, 1))))
    assert_false(is_multilabel(np.random.randint(2, size=(10, ))))
    assert_false(is_multilabel(np.random.randint(2, size=(10, 1))))


def test_is_label_indicator_matrix():
    assert_true(is_label_indicator_matrix(np.random.randint(2, size=(10, 10))))

    assert_false(is_label_indicator_matrix([[1], [2], [0, 1]]))
    assert_false(is_label_indicator_matrix([[1], [2]]))
    assert_false(is_label_indicator_matrix([[1], [2], []]))
    assert_false(is_label_indicator_matrix([[1], [0, 2], []]))
    assert_false(is_label_indicator_matrix(range(10)))
    assert_false(is_label_indicator_matrix(np.arange(10)))
    assert_false(is_label_indicator_matrix(np.reshape(np.arange(9), (3, 3))))
    assert_false(is_label_indicator_matrix(np.reshape(np.arange(10), (-1, 1))))
    assert_false(is_label_indicator_matrix(np.reshape(np.arange(10), (1, -1))))
    assert_false(is_label_indicator_matrix(np.random.randint(2, size=(10, ))))
    assert_false(is_label_indicator_matrix(np.random.randint(2, size=(10, 1))))
