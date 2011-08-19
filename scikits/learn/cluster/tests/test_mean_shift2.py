import numpy as np
from nose.tools import assert_equal, assert_true

from scikits.learn.cluster import mean_shift, get_bin_seeds

def test_bin_seeds():
    # Data is just 6 points in the plane
    X = np.array([[1., 1.], [1.5, 1.5], [1.8, 1.2],
                  [2., 1.], [2.1, 1.1], [0., 0.]])
    
    # With a bin coarseness of 1.0 and min_bin_freq of 1, 3 bins should be found
    ground_truth = set([(1.,1.), (2.,1.), (0.,0.)])
    test_bins = get_bin_seeds(X, 1, 1)
    test_result = set([tuple(p) for p in test_bins])
    assert_true(len(ground_truth.symmetric_difference(test_result)) == 0)
    
    # With a bin coarseness of 1.0 and min_bin_freq of 2, 2 bins should be found
    ground_truth = set([(1.,1.), (2.,1.)])
    test_bins = get_bin_seeds(X, 1, 2)
    test_result = set([tuple(p) for p in test_bins])
    assert_true(len(ground_truth.symmetric_difference(test_result)) == 0)

    # With a bin size of 0.01 and min_bin_freq of 1, 6 bins should be found
    test_bins = get_bin_seeds(X, 0.01, 1)
    test_result = set([tuple(p) for p in test_bins])
    assert_true(len(test_result) == 6)

    
