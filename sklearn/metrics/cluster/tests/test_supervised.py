import numpy as np

from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics.cluster import entropy
from sklearn.metrics.cluster import expected_mutual_information
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_completeness_v_measure
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import mutual_info_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import v_measure_score
from sklearn.metrics.cluster.supervised import _generalized_average

from sklearn.utils import assert_all_finite
from sklearn.utils.testing import (
        assert_equal, assert_almost_equal, assert_raise_message,
        assert_warns_message, ignore_warnings
)
from numpy.testing import assert_array_almost_equal


score_funcs = [
    adjusted_rand_score,
    homogeneity_score,
    completeness_score,
    v_measure_score,
    adjusted_mutual_info_score,
    normalized_mutual_info_score,
]


def test_future_warning():
    score_funcs_with_changing_means = [
        normalized_mutual_info_score,
        adjusted_mutual_info_score,
    ]
    warning_msg = "The behavior of "
    args = [0, 0, 0], [0, 0, 0]
    for score_func in score_funcs_with_changing_means:
        assert_warns_message(FutureWarning, warning_msg, score_func, *args)


@ignore_warnings(category=FutureWarning)
def test_error_messages_on_wrong_input():
    for score_func in score_funcs:
        expected = ('labels_true and labels_pred must have same size,'
                    ' got 2 and 3')
        assert_raise_message(ValueError, expected, score_func,
                             [0, 1], [1, 1, 1])

        expected = "labels_true must be 1D: shape is (2"
        assert_raise_message(ValueError, expected, score_func,
                             [[0, 1], [1, 0]], [1, 1, 1])

        expected = "labels_pred must be 1D: shape is (2"
        assert_raise_message(ValueError, expected, score_func,
                             [0, 1, 0], [[1, 1], [0, 0]])


def test_generalized_average():
    a, b = 1, 2
    methods = ["min", "geometric", "arithmetic", "max"]
    means = [_generalized_average(a, b, method) for method in methods]
    assert means[0] <= means[1] <= means[2] <= means[3]
    c, d = 12, 12
    means = [_generalized_average(c, d, method) for method in methods]
    assert means[0] == means[1] == means[2] == means[3]


@ignore_warnings(category=FutureWarning)
def test_perfect_matches():
    for score_func in score_funcs:
        assert_equal(score_func([], []), 1.0)
        assert_equal(score_func([0], [1]), 1.0)
        assert_equal(score_func([0, 0, 0], [0, 0, 0]), 1.0)
        assert_equal(score_func([0, 1, 0], [42, 7, 42]), 1.0)
        assert_equal(score_func([0., 1., 0.], [42., 7., 42.]), 1.0)
        assert_equal(score_func([0., 1., 2.], [42., 7., 2.]), 1.0)
        assert_equal(score_func([0, 1, 2], [42, 7, 2]), 1.0)
    score_funcs_with_changing_means = [
        normalized_mutual_info_score,
        adjusted_mutual_info_score,
    ]
    means = {"min", "geometric", "arithmetic", "max"}
    for score_func in score_funcs_with_changing_means:
        for mean in means:
            assert score_func([], [], mean) == 1.0
            assert score_func([0], [1], mean) == 1.0
            assert score_func([0, 0, 0], [0, 0, 0], mean) == 1.0
            assert score_func([0, 1, 0], [42, 7, 42], mean) == 1.0
            assert score_func([0., 1., 0.], [42., 7., 42.], mean) == 1.0
            assert score_func([0., 1., 2.], [42., 7., 2.], mean) == 1.0
            assert score_func([0, 1, 2], [42, 7, 2], mean) == 1.0


def test_homogeneous_but_not_complete_labeling():
    # homogeneous but not complete clustering
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, 1, 2, 2])
    assert_almost_equal(h, 1.00, 2)
    assert_almost_equal(c, 0.69, 2)
    assert_almost_equal(v, 0.81, 2)


def test_complete_but_not_homogeneous_labeling():
    # complete but not homogeneous clustering
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 1, 1, 2, 2],
        [0, 0, 1, 1, 1, 1])
    assert_almost_equal(h, 0.58, 2)
    assert_almost_equal(c, 1.00, 2)
    assert_almost_equal(v, 0.73, 2)


def test_not_complete_and_not_homogeneous_labeling():
    # neither complete nor homogeneous but not so bad either
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 2, 2])
    assert_almost_equal(h, 0.67, 2)
    assert_almost_equal(c, 0.42, 2)
    assert_almost_equal(v, 0.52, 2)


def test_beta_parameter():
    # test for when beta passed to
    # homogeneity_completeness_v_measure
    # and v_measure_score
    beta_test = 0.2
    h_test = 0.67
    c_test = 0.42
    v_test = ((1 + beta_test) * h_test * c_test
              / (beta_test * h_test + c_test))

    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 2, 2],
        beta=beta_test)
    assert_almost_equal(h, h_test, 2)
    assert_almost_equal(c, c_test, 2)
    assert_almost_equal(v, v_test, 2)

    v = v_measure_score(
        [0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 2, 2],
        beta=beta_test)
    assert_almost_equal(v, v_test, 2)


def test_non_consecutive_labels():
    # regression tests for labels with gaps
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 2, 2, 2],
        [0, 1, 0, 1, 2, 2])
    assert_almost_equal(h, 0.67, 2)
    assert_almost_equal(c, 0.42, 2)
    assert_almost_equal(v, 0.52, 2)

    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 4, 0, 4, 2, 2])
    assert_almost_equal(h, 0.67, 2)
    assert_almost_equal(c, 0.42, 2)
    assert_almost_equal(v, 0.52, 2)

    ari_1 = adjusted_rand_score([0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 2, 2])
    ari_2 = adjusted_rand_score([0, 0, 0, 1, 1, 1], [0, 4, 0, 4, 2, 2])
    assert_almost_equal(ari_1, 0.24, 2)
    assert_almost_equal(ari_2, 0.24, 2)


@ignore_warnings(category=FutureWarning)
def uniform_labelings_scores(score_func, n_samples, k_range, n_runs=10,
                             seed=42):
    # Compute score for random uniform cluster labelings
    random_labels = np.random.RandomState(seed).randint
    scores = np.zeros((len(k_range), n_runs))
    for i, k in enumerate(k_range):
        for j in range(n_runs):
            labels_a = random_labels(low=0, high=k, size=n_samples)
            labels_b = random_labels(low=0, high=k, size=n_samples)
            scores[i, j] = score_func(labels_a, labels_b)
    return scores


@ignore_warnings(category=FutureWarning)
def test_adjustment_for_chance():
    # Check that adjusted scores are almost zero on random labels
    n_clusters_range = [2, 10, 50, 90]
    n_samples = 100
    n_runs = 10

    scores = uniform_labelings_scores(
        adjusted_rand_score, n_samples, n_clusters_range, n_runs)

    max_abs_scores = np.abs(scores).max(axis=1)
    assert_array_almost_equal(max_abs_scores, [0.02, 0.03, 0.03, 0.02], 2)


@ignore_warnings(category=FutureWarning)
def test_adjusted_mutual_info_score():
    # Compute the Adjusted Mutual Information and test against known values
    labels_a = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
    labels_b = np.array([1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 3, 3, 3, 2, 2])
    # Mutual information
    mi = mutual_info_score(labels_a, labels_b)
    assert_almost_equal(mi, 0.41022, 5)
    # with provided sparse contingency
    C = contingency_matrix(labels_a, labels_b, sparse=True)
    mi = mutual_info_score(labels_a, labels_b, contingency=C)
    assert_almost_equal(mi, 0.41022, 5)
    # with provided dense contingency
    C = contingency_matrix(labels_a, labels_b)
    mi = mutual_info_score(labels_a, labels_b, contingency=C)
    assert_almost_equal(mi, 0.41022, 5)
    # Expected mutual information
    n_samples = C.sum()
    emi = expected_mutual_information(C, n_samples)
    assert_almost_equal(emi, 0.15042, 5)
    # Adjusted mutual information
    ami = adjusted_mutual_info_score(labels_a, labels_b)
    assert_almost_equal(ami, 0.27502, 5)
    ami = adjusted_mutual_info_score([1, 1, 2, 2], [2, 2, 3, 3])
    assert_equal(ami, 1.0)
    # Test with a very large array
    a110 = np.array([list(labels_a) * 110]).flatten()
    b110 = np.array([list(labels_b) * 110]).flatten()
    ami = adjusted_mutual_info_score(a110, b110)
    # This is not accurate to more than 2 places
    assert_almost_equal(ami, 0.37, 2)


def test_expected_mutual_info_overflow():
    # Test for regression where contingency cell exceeds 2**16
    # leading to overflow in np.outer, resulting in EMI > 1
    assert expected_mutual_information(np.array([[70000]]), 70000) <= 1


def test_int_overflow_mutual_info_fowlkes_mallows_score():
    # Test overflow in mutual_info_classif and fowlkes_mallows_score
    x = np.array([1] * (52632 + 2529) + [2] * (14660 + 793) + [3] * (3271 +
                 204) + [4] * (814 + 39) + [5] * (316 + 20))
    y = np.array([0] * 52632 + [1] * 2529 + [0] * 14660 + [1] * 793 +
                 [0] * 3271 + [1] * 204 + [0] * 814 + [1] * 39 + [0] * 316 +
                 [1] * 20)

    assert_all_finite(mutual_info_score(x, y))
    assert_all_finite(fowlkes_mallows_score(x, y))


def test_entropy():
    ent = entropy([0, 0, 42.])
    assert_almost_equal(ent, 0.6365141, 5)
    assert_almost_equal(entropy([]), 1)


def test_contingency_matrix():
    labels_a = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
    labels_b = np.array([1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 3, 3, 3, 2, 2])
    C = contingency_matrix(labels_a, labels_b)
    C2 = np.histogram2d(labels_a, labels_b,
                        bins=(np.arange(1, 5),
                              np.arange(1, 5)))[0]
    assert_array_almost_equal(C, C2)
    C = contingency_matrix(labels_a, labels_b, eps=.1)
    assert_array_almost_equal(C, C2 + .1)


def test_contingency_matrix_sparse():
    labels_a = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
    labels_b = np.array([1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 3, 3, 3, 2, 2])
    C = contingency_matrix(labels_a, labels_b)
    C_sparse = contingency_matrix(labels_a, labels_b, sparse=True).toarray()
    assert_array_almost_equal(C, C_sparse)
    C_sparse = assert_raise_message(ValueError,
                                    "Cannot set 'eps' when sparse=True",
                                    contingency_matrix, labels_a, labels_b,
                                    eps=1e-10, sparse=True)


@ignore_warnings(category=FutureWarning)
def test_exactly_zero_info_score():
    # Check numerical stability when information is exactly zero
    for i in np.logspace(1, 4, 4).astype(np.int):
        labels_a, labels_b = (np.ones(i, dtype=np.int),
                              np.arange(i, dtype=np.int))
        assert_equal(normalized_mutual_info_score(labels_a, labels_b), 0.0)
        assert_equal(v_measure_score(labels_a, labels_b), 0.0)
        assert_equal(adjusted_mutual_info_score(labels_a, labels_b), 0.0)
        assert_equal(normalized_mutual_info_score(labels_a, labels_b), 0.0)
        for method in ["min", "geometric", "arithmetic", "max"]:
            assert adjusted_mutual_info_score(labels_a, labels_b,
                                              method) == 0.0
            assert normalized_mutual_info_score(labels_a, labels_b,
                                                method) == 0.0


def test_v_measure_and_mutual_information(seed=36):
    # Check relation between v_measure, entropy and mutual information
    for i in np.logspace(1, 4, 4).astype(np.int):
        random_state = np.random.RandomState(seed)
        labels_a, labels_b = (random_state.randint(0, 10, i),
                              random_state.randint(0, 10, i))
        assert_almost_equal(v_measure_score(labels_a, labels_b),
                            2.0 * mutual_info_score(labels_a, labels_b) /
                            (entropy(labels_a) + entropy(labels_b)), 0)
        avg = 'arithmetic'
        assert_almost_equal(v_measure_score(labels_a, labels_b),
                            normalized_mutual_info_score(labels_a, labels_b,
                                                         average_method=avg)
                            )


def test_fowlkes_mallows_score():
    # General case
    score = fowlkes_mallows_score([0, 0, 0, 1, 1, 1],
                                  [0, 0, 1, 1, 2, 2])
    assert_almost_equal(score, 4. / np.sqrt(12. * 6.))

    # Perfect match but where the label names changed
    perfect_score = fowlkes_mallows_score([0, 0, 0, 1, 1, 1],
                                          [1, 1, 1, 0, 0, 0])
    assert_almost_equal(perfect_score, 1.)

    # Worst case
    worst_score = fowlkes_mallows_score([0, 0, 0, 0, 0, 0],
                                        [0, 1, 2, 3, 4, 5])
    assert_almost_equal(worst_score, 0.)


def test_fowlkes_mallows_score_properties():
    # handcrafted example
    labels_a = np.array([0, 0, 0, 1, 1, 2])
    labels_b = np.array([1, 1, 2, 2, 0, 0])
    expected = 1. / np.sqrt((1. + 3.) * (1. + 2.))
    # FMI = TP / sqrt((TP + FP) * (TP + FN))

    score_original = fowlkes_mallows_score(labels_a, labels_b)
    assert_almost_equal(score_original, expected)

    # symmetric property
    score_symmetric = fowlkes_mallows_score(labels_b, labels_a)
    assert_almost_equal(score_symmetric, expected)

    # permutation property
    score_permuted = fowlkes_mallows_score((labels_a + 1) % 3, labels_b)
    assert_almost_equal(score_permuted, expected)

    # symmetric and permutation(both together)
    score_both = fowlkes_mallows_score(labels_b, (labels_a + 2) % 3)
    assert_almost_equal(score_both, expected)
