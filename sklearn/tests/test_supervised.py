import numpy as np

from sklearn.utils._testing import (assert_array_equal,
                                    assert_array_almost_equal)

from sklearn.metrics import pair_confusion_matrix, rand_score


def test_pair_confusion_matrix():
    n = 10
    N = n * n

    # edge case: every element is its own cluster
    clustering1 = list(range(N))
    clustering2 = clustering1
    expected = np.array([[N * (N - 1), 0], [0, 0]])
    assert_array_equal(expected, pair_confusion_matrix(clustering1,
                                                       clustering2))

    # edge case: only one cluster
    clustering1 = np.full((N,), fill_value=0)
    clustering2 = clustering1
    expected = np.array([[0, 0], [0, N * (N - 1)]])
    assert_array_equal(expected, pair_confusion_matrix(clustering1,
                                                       clustering2))

    # regular case: different non-trivial clusterings
    clustering1 = np.array([i+1 for i in range(n) for j in range(n)])
    clustering2 = np.array([i+1 for i in range(n) for j in range(n+1)][:N])
    # basic quadratic implementation
    expected = np.full(shape=(2, 2), fill_value=0, dtype=np.int64)
    for i in range(len(clustering1)):
        for j in range(len(clustering1)):
            if i != j:
                same_cluster_1 = int(clustering1[i] == clustering1[j])
                same_cluster_2 = int(clustering2[i] == clustering2[j])
                expected[same_cluster_1, same_cluster_2] += 1
    assert_array_equal(expected, pair_confusion_matrix(clustering1,
                                                       clustering2))


def test_rand_score():
    n = 10
    N = n * n

    # edge case: every element is its own cluster
    clustering1 = list(range(N))
    clustering2 = clustering1
    expected = np.array([1.])
    assert_array_almost_equal(expected, rand_score(clustering1,
                                                   clustering2))

    # edge case: only one cluster
    clustering1 = np.full((N,), fill_value=0)
    clustering2 = clustering1
    expected = np.array([1.])
    assert_array_almost_equal(expected, rand_score(clustering1,
                                                   clustering2))

    # regular case: different non-trivial clusterings
    clustering1 = np.array([i+1 for i in range(n) for j in range(n)])
    clustering2 = np.array([i+1 for i in range(n) for j in range(n+1)][:N])
    # basic quadratic implementation
    expected = np.full(shape=(2, 2), fill_value=0, dtype=np.int64)
    for i in range(len(clustering1)):
        for j in range(len(clustering1)):
            if i != j:
                same_cluster_1 = int(clustering1[i] == clustering1[j])
                same_cluster_2 = int(clustering2[i] == clustering2[j])
                expected[same_cluster_1, same_cluster_2] += 1
        expected_numerator = expected[0, 0] + expected[1, 1]
    expected_denominator = (expected[0, 0] + expected[1, 1] + expected[0, 1] +
                            expected[1, 0])
    expected_score = expected_numerator / expected_denominator
    assert_array_almost_equal(expected_score, rand_score(clustering1,
                                                         clustering2))
