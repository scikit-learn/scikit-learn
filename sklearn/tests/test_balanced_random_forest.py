from sklearn.ensemble.forest import\
    _get_class_balance_data, _generate_balanced_sample_indices
import numpy as np
from numpy.testing import assert_array_equal


def test_get_class_balance_data():
    y = np.array([0, 1, 0, 1, 1, 2])
    classes, class_counts, class_indices = _get_class_balance_data(y)
    assert_array_equal(classes, [0, 1, 2])
    assert_array_equal(class_counts, [2, 3, 1])
    assert_array_equal(class_indices[0], [0, 2])
    assert_array_equal(class_indices[1], [1, 3, 4])
    assert_array_equal(class_indices[2], [5])


def test_generate_balanced_sample_indices():
    y = np.array([0, 1, 0, 1, 1, 2])
    random_state = 0
    balance_data = _get_class_balance_data(y)
    sample_indices = _generate_balanced_sample_indices(random_state,
                                                       balance_data)
    assert_array_equal(sample_indices, [0, 3, 5])
