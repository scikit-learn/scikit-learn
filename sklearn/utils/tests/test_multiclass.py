import numpy as np

from sklearn.externals.six.moves import xrange

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises

from sklearn.utils.multiclass import unique_labels
from sklearn.utils.multiclass import is_label_indicator_matrix
from sklearn.utils.multiclass import is_multilabel
from sklearn.utils.multiclass import is_sequence_of_sequences
from sklearn.utils.multiclass import type_of_target

INDICATOR_EXAMPLES = [
    np.random.randint(2, size=(10, 10)),
    np.array([[0, 1], [1, 0]]),
    np.array([[0, 0], [0, 0]]),
    np.array([[-1, 1], [1, -1]]),
    np.array([[-3, 3], [3, -3]]),
]

SEQUENCES_EXAMPLES = [
    [[0, 1]],
    [[0], [1]],
    [[1, 2, 3]],
    [[1], [2], [0, 1]],
    [[1], [2]],
    [[]],
    [()],
    np.array([[], [1, 2]], dtype='object'),
]

MULTICLASS_EXAMPLES = [
    [1, 0, 2, 2, 1, 4, 2, 4, 4, 4],
    np.array([1, 0, 2]),
    np.array([[1], [0], [2]]),
    np.array([[1, 0, 2]]),
    [],
    [0],
    [0, 1, 2],
    ['a', 'b', 'c'],
]

BINARY_EXAMPLES = [
    [0, 1],
    np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1]),
    np.array([[0], [1]]),
    np.array([[0, 1]]),
    [1, -1],
    [3, 5],
    ['a', 'b'],
    ['abc', 'def'],
]

CONTINUOUS_EXAMPLES = [
    [1e-5],
    [0, .5],
]

BAD_EXAMPLES = [
    # Could be confused for indicator or other:
    np.array([[0.1, 0.2], [0.2, 0.1]]),
    np.array([[1, 2], [3, 1]]),
    # not currently supported sequence of sequences
    np.array([np.array([]), np.array([1, 2, 3])], dtype=object),
]


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
    for example in SEQUENCES_EXAMPLES + INDICATOR_EXAMPLES:
        assert_true(is_multilabel(example),
                    msg='is_multilabel(%r) should be True' % example)

    for example in (MULTICLASS_EXAMPLES + BINARY_EXAMPLES +
                    CONTINUOUS_EXAMPLES + BAD_EXAMPLES):
        assert_false(is_multilabel(example),
                     msg='is_multilabel(%r) should be False' % example)


def test_is_label_indicator_matrix():
    for example in INDICATOR_EXAMPLES:
        assert_true(is_label_indicator_matrix(example),
                    msg='is_label_indicator_matrix(%r) should be True'
                    % example)

    for example in (SEQUENCES_EXAMPLES + MULTICLASS_EXAMPLES +
                    BINARY_EXAMPLES + CONTINUOUS_EXAMPLES + BAD_EXAMPLES):
        assert_false(is_label_indicator_matrix(example),
                     msg='is_label_indicator_matrix(%r) should be False'
                     % example)


def test_is_sequence_of_sequences():
    for example in SEQUENCES_EXAMPLES:
        assert_true(is_sequence_of_sequences(example),
                    msg='is_sequence_of_sequences(%r) should be True'
                    % example)

    for example in (INDICATOR_EXAMPLES + MULTICLASS_EXAMPLES +
                    BINARY_EXAMPLES + CONTINUOUS_EXAMPLES + BAD_EXAMPLES):
        assert_false(is_sequence_of_sequences(example),
                     msg='is_sequence_of_sequences(%r) should be False'
                     % example)


def test_type_of_target():
    for example in SEQUENCES_EXAMPLES:
        assert_equal(type_of_target(example), 'multilabel-sequences')
    for example in INDICATOR_EXAMPLES:
        assert_equal(type_of_target(example), 'multilabel-indicator')
    for example in MULTICLASS_EXAMPLES:
        assert_equal(type_of_target(example), 'multiclass')
    for example in BINARY_EXAMPLES:
        assert_equal(type_of_target(example), 'binary')
    for example in CONTINUOUS_EXAMPLES:
        assert_equal(type_of_target(example), 'continuous')
#    for example in BAD_EXAMPLES:
#        assert_raises(ValueError, type_of_target, example)
