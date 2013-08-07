import numpy as np
from itertools import product
from sklearn.externals.six.moves import xrange
from sklearn.externals.six import iteritems

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

EXAMPLES = {
    'multilabel-indicator': [
        np.random.RandomState(42).randint(2, size=(10, 10)),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, 1], [1, 0]], dtype=np.bool),
        np.array([[0, 1], [1, 0]], dtype=np.int8),
        np.array([[0, 1], [1, 0]], dtype=np.uint8),
        np.array([[0, 1], [1, 0]], dtype=np.float),
        np.array([[0, 1], [1, 0]], dtype=np.float32),
        np.array([[0, 0], [0, 0]]),
        np.array([[-1, 1], [1, -1]]),
        np.array([[-3, 3], [3, -3]]),
        np.array([[0, 1]]),
    ],
    'multilabel-sequences': [
        [[0, 1]],
        [[0], [1]],
        [[1, 2, 3]],
        [[1, 2, 1]],  # duplicate values, why not?
        [[1], [2], [0, 1]],
        [[1], [2]],
        [[]],
        [()],
        np.array([[], [1, 2]], dtype='object'),
    ],
    'multiclass': [
        [1, 0, 2, 2, 1, 4, 2, 4, 4, 4],
        np.array([1, 0, 2]),
        np.array([1, 0, 2], dtype=np.int8),
        np.array([1, 0, 2], dtype=np.uint8),
        np.array([1, 0, 2], dtype=np.float),
        np.array([1, 0, 2], dtype=np.float32),
        np.array([[1], [0], [2]]),
        [0, 1, 2],
        ['a', 'b', 'c'],
        np.array([u'a', u'b', u'c']),
    ],
    'multiclass-multioutput': [
        np.array([[1, 0, 2, 2], [1, 4, 2, 4]]),
        np.array([[1, 0, 2, 2], [1, 4, 2, 4]], dtype=np.int8),
        np.array([[1, 0, 2, 2], [1, 4, 2, 4]], dtype=np.uint8),
        np.array([[1, 0, 2, 2], [1, 4, 2, 4]], dtype=np.float),
        np.array([[1, 0, 2, 2], [1, 4, 2, 4]], dtype=np.float32),
        np.array([['a', 'b'], ['c', 'd']]),
        np.array([[u'a', u'b'], [u'c', u'd']]),
        np.array([[1, 0, 2]]),
    ],
    'binary': [
        [0, 1],
        [1, 1],
        [],
        [0],
        np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1]),
        np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1], dtype=np.bool),
        np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1], dtype=np.int8),
        np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1], dtype=np.uint8),
        np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1], dtype=np.float),
        np.array([0, 1, 1, 1, 0, 0, 0, 1, 1, 1], dtype=np.float32),
        np.array([[0], [1]]),
        [1, -1],
        [3, 5],
        ['a'],
        ['a', 'b'],
        ['abc', 'def'],
        [u'a', u'b'],
    ],
    'continuous': [
        [1e-5],
        [0, .5],
        np.array([[0], [.5]]),
        np.array([[0], [.5]], dtype=np.float32),
    ],
    'continuous-multioutput': [
        np.array([[0, .5], [.5, 0]]),
        np.array([[0, .5], [.5, 0]], dtype=np.float32),
        np.array([[0, .5]]),
    ],
    'unknown': [
        # empty second dimension
        np.array([[], []]),
        # 3d
        np.array([[[0, 1], [2, 3]], [[4, 5], [6, 7]]]),
        # not currently supported sequence of sequences
        np.array([np.array([]), np.array([1, 2, 3])], dtype=object),
        [np.array([]), np.array([1, 2, 3])],
        [set([1, 2, 3]), set([1, 2])],
        [frozenset([1, 2, 3]), frozenset([1, 2])],
        # and also confusable as sequences of sequences
        [{0: 'a', 1: 'b'}, {0: 'a'}],
    ]
}

NON_ARRAY_LIKE_EXAMPLES = [
    set([1, 2, 3]),
    {0: 'a', 1: 'b'},
    {0: [5], 1: [5]},
    'abc',
    frozenset([1, 2, 3]),
    None,
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

    assert_array_equal(unique_labels(np.array([[0, 0, 1],
                                               [0, 0, 0]])),
                       np.arange(3))

    # Several arrays passed
    assert_array_equal(unique_labels([4, 0, 2], xrange(5)),
                       np.arange(5))
    assert_array_equal(unique_labels((0, 1, 2), (0,), (2, 1)),
                       np.arange(3))

    # Border line case with binary indicator matrix
    assert_raises(ValueError, unique_labels, [4, 0, 2], np.ones((5, 5)))
    assert_raises(ValueError, unique_labels, np.ones((5, 4)), np.ones((5, 5)))
    assert_array_equal(unique_labels(np.ones((4, 5)), np.ones((5, 5))),
                       np.arange(5))

    # Some tests with strings input
    assert_array_equal(unique_labels(["a", "b", "c"], ["d"]),
                       ["a", "b", "c", "d"])
    assert_array_equal(unique_labels([["a", "b"], ["c"]], [["d"]]),
                       ["a", "b", "c", "d"])

    # Smoke test for all supported format
    for format in ["binary", "multiclass", "multilabel-sequences",
                   "multilabel-indicator"]:
        for y in EXAMPLES[format]:
            unique_labels(y)

    # We don't support those format at the moment
    for example in NON_ARRAY_LIKE_EXAMPLES:
        assert_raises(ValueError, unique_labels, example)

    for y_type in ["unknown", "continuous", 'continuous-multioutput',
                   'multiclass-multioutput']:
        for example in EXAMPLES[y_type]:
            assert_raises(ValueError, unique_labels, example)

    #Mix of multilabel-indicator and multilabel-sequences
    mix_multilabel_format = product(EXAMPLES["multilabel-indicator"],
                                    EXAMPLES["multilabel-sequences"])
    for y_multilabel, y_multiclass in mix_multilabel_format:
        assert_raises(ValueError, unique_labels, y_multiclass, y_multilabel)
        assert_raises(ValueError, unique_labels, y_multilabel, y_multiclass)

    #Mix with binary or multiclass and multilabel
    mix_clf_format = product(EXAMPLES["multilabel-indicator"] +
                             EXAMPLES["multilabel-sequences"],
                             EXAMPLES["multiclass"] +
                             EXAMPLES["binary"])

    for y_multilabel, y_multiclass in mix_clf_format:
        assert_raises(ValueError, unique_labels, y_multiclass, y_multilabel)
        assert_raises(ValueError, unique_labels, y_multilabel, y_multiclass)

    # Mix string and number input type
    assert_raises(ValueError, unique_labels, [[1, 2], [3]],
                  [["a", "d"]])
    assert_raises(ValueError, unique_labels, ["1", 2])
    assert_raises(ValueError, unique_labels, [["1", 2], [3]])
    assert_raises(ValueError, unique_labels, [["1", "2"], [3]])

    assert_array_equal(unique_labels([(2,), (0, 2,)], [(), ()]), [0, 2])
    assert_array_equal(unique_labels([("2",), ("0", "2",)], [(), ()]),
                       ["0", "2"])


def test_is_multilabel():
    for group, group_examples in iteritems(EXAMPLES):
        if group.startswith('multilabel'):
            assert_, exp = assert_true, 'True'
        else:
            assert_, exp = assert_false, 'False'
        for example in group_examples:
            assert_(is_multilabel(example),
                    msg='is_multilabel(%r) should be %s' % (example, exp))


def test_is_label_indicator_matrix():
    for group, group_examples in iteritems(EXAMPLES):
        if group == 'multilabel-indicator':
            assert_, exp = assert_true, 'True'
        else:
            assert_, exp = assert_false, 'False'
        for example in group_examples:
            assert_(is_label_indicator_matrix(example),
                    msg='is_label_indicator_matrix(%r) should be %s'
                    % (example, exp))


def test_is_sequence_of_sequences():
    for group, group_examples in iteritems(EXAMPLES):
        if group == 'multilabel-sequences':
            assert_, exp = assert_true, 'True'
        else:
            assert_, exp = assert_false, 'False'
        for example in group_examples:
            assert_(is_sequence_of_sequences(example),
                    msg='is_sequence_of_sequences(%r) should be %s'
                    % (example, exp))


def test_type_of_target():
    for group, group_examples in iteritems(EXAMPLES):
        for example in group_examples:
            assert_equal(type_of_target(example), group,
                         msg='type_of_target(%r) should be %r, got %r'
                         % (example, group, type_of_target(example)))

    for example in NON_ARRAY_LIKE_EXAMPLES:
        assert_raises(ValueError, type_of_target, example)
