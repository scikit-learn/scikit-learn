import warnings

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises

from sklearn.feature_selection.by_score import mask_by_score
from sklearn.feature_selection.by_score import SelectByScore

SCORES = np.array([2, 1, np.nan, 5, 1])
SCORES_NO_NAN = np.array([2, 1, 0, 5, 1])


def test_no_args():
    assert_equal(mask_by_score(SCORES).dtype, np.bool)
    assert_array_equal(mask_by_score(SCORES), [1, 1, 1, 1, 1])


def test_max_min_number():
    assert_equal(mask_by_score(SCORES, minimum=2).dtype, np.bool)
    assert_array_equal(mask_by_score(SCORES, minimum=2), [1, 0, 0, 1, 0])
    assert_array_equal(mask_by_score(SCORES, minimum=2.1), [0, 0, 0, 1, 0])
    assert_array_equal(mask_by_score(SCORES, maximum=2), [1, 1, 0, 0, 1])
    assert_array_equal(mask_by_score(SCORES, maximum=1.9), [0, 1, 0, 0, 1])
    assert_array_equal(mask_by_score(SCORES, minimum=2, maximum=4),
                       [1, 0, 0, 0, 0])


def test_max_min_string():
    tests = [
        ('mean', 1.8),
        ('.5 * mean', 0.9),
        ('median', 1.0),
        ('min', 0.0),
        ('max', 5.0),
        ('sum * .5', 4.5),
        ('length * .5', 2.5),
    ]
    for arg, val in tests:
        print arg, val
        assert_array_equal(mask_by_score(SCORES_NO_NAN, minimum=arg),
                           mask_by_score(SCORES_NO_NAN, minimum=val))
        assert_array_equal(mask_by_score(SCORES_NO_NAN, maximum=arg),
                           mask_by_score(SCORES_NO_NAN, maximum=val))


def test_max_min_callable():
    val = 2

    def fn(scores):
        assert_array_equal(scores, SCORES)
        return 2

    assert_array_equal(mask_by_score(SCORES, minimum=fn),
                       mask_by_score(SCORES, minimum=val))
    assert_array_equal(mask_by_score(SCORES, maximum=fn),
                       mask_by_score(SCORES, maximum=val))


def test_limit_int():
    assert_equal(mask_by_score(SCORES, limit=2).dtype, np.bool)
    assert_array_equal(mask_by_score(SCORES, limit=2), [1, 0, 0, 1, 0])
    assert_array_equal(mask_by_score(SCORES, limit=-2), [0, 1, 0, 0, 1])

    # Bounds
    assert_array_equal(mask_by_score(SCORES, limit=0), [0, 0, 0, 0, 0])
    assert_array_equal(mask_by_score(SCORES, limit=5), [1, 1, 1, 1, 1])
    assert_array_equal(mask_by_score(SCORES, limit=10), [1, 1, 1, 1, 1])

    # Tie-breaking
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        top3 = mask_by_score(SCORES, limit=3)
        assert_equal(top3.sum(), 3)
        assert_equal((SCORES[top3] == 1).sum(), 1)
        assert_equal(len(w), 1)
        del w[:]
        bottom1 = mask_by_score(SCORES, limit=-1)
        assert_equal(bottom1.sum(), 1)
        assert_equal((SCORES[bottom1] == 1).sum(), 1)
        assert_equal(len(w), 1)

    # Combination with min/maximum
    assert_array_equal(mask_by_score(SCORES, limit=1, maximum=3),
                       [1, 0, 0, 0, 0])
    assert_array_equal(mask_by_score(SCORES, limit=-1, minimum=3),
                       [0, 0, 0, 1, 0])

    # Checking nan support
    assert_array_equal(mask_by_score(SCORES, limit=4), [1, 1, 0, 1, 1])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        assert_equal(mask_by_score([np.nan, np.nan], limit=1).sum(), 1)
        assert_equal(len(w), 1)
        del w[:]


def test_limit_float():
    assert_raises(ValueError, mask_by_score, SCORES, limit=1.1)
    assert_raises(ValueError, mask_by_score, SCORES, limit=-1.1)
    assert_raises(ValueError, mask_by_score, SCORES, limit=2.0)
    assert_raises(ValueError, mask_by_score, SCORES, limit=-2.0)
    assert_array_equal(mask_by_score(SCORES, limit=0.4), [1, 0, 0, 1, 0])
    assert_array_equal(mask_by_score(SCORES, limit=-0.4), [0, 1, 0, 0, 1])
    assert_array_equal(mask_by_score(SCORES, limit=0.5), [1, 0, 0, 1, 0])
    assert_array_equal(mask_by_score(SCORES, limit=-0.5), [0, 1, 0, 0, 1])

    # Bounds
    assert_array_equal(mask_by_score(SCORES, limit=0.0), [0, 0, 0, 0, 0])
    assert_array_equal(mask_by_score(SCORES, limit=1.0), [1, 1, 1, 1, 1])

    # Tie-breaking
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        top3 = mask_by_score(SCORES, limit=0.6)
        assert_equal(top3.sum(), 3)
        assert_equal((SCORES[top3] == 1).sum(), 1)
        assert_equal(len(w), 1)
        del w[:]
        bottom1 = mask_by_score(SCORES, limit=-0.2)
        assert_equal(bottom1.sum(), 1)
        assert_equal((SCORES[bottom1] == 1).sum(), 1)
        assert_equal(len(w), 1)

    # Combination with min/maximum
    assert_array_equal(mask_by_score(SCORES, limit=0.2, maximum=3),
                       [1, 0, 0, 0, 0])
    assert_array_equal(mask_by_score(SCORES, limit=-0.2, minimum=3),
                       [0, 0, 0, 1, 0])

    # Checking nan support
    assert_array_equal(mask_by_score(SCORES, limit=0.8), [1, 1, 0, 1, 1])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        assert_equal(mask_by_score([np.nan, np.nan], limit=0.5).sum(), 1)
        assert_equal(len(w), 1)
        del w[:]


def test_transformer():
    transformer = SelectByScore(lambda X, y: [5, 3, 2, 4, 1],
                                minimum=2, maximum=4, limit=2)
    X = [[10, 20, 30, 40, 50]]
    assert_array_equal(transformer.fit_transform(X, [1]), [[20, 40]])
