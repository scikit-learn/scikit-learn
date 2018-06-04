import sys

import numpy as np
from scipy import sparse as sp

import pytest

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import if_safe_multiprocessing_with_blas
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.validation import _num_samples
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning

from sklearn.utils.extmath import row_norms
from sklearn.metrics.cluster import v_measure_score
from sklearn.cluster import KMeans, k_means
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster.k_means_ import _labels_inertia
from sklearn.cluster.k_means_ import _mini_batch_step
from sklearn.datasets.samples_generator import make_blobs
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.metrics.cluster import homogeneity_score


# non centered, sparse centers to check the
centers = np.array([
    [1.0, 2.0, 0.0, 1.0, 0.0],
    [1.0, 1.0, 3.0, 0.0, 0.0],
    [1.0, 0.0, 2.0, 1.0, 1.0],
])
n_samples = 100
n_clusters, n_features = centers.shape
X, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                            cluster_std=1., random_state=42)
X_csr = sp.csr_matrix(X)

def _sort_centers(centers):
    return np.sort(centers, axis=0)

data_args = {
            "init":"k-means++",
            "n_clusters":n_clusters,
            "random_state":42
    }

def weight_invariance_kmeans(model1, model2):
    """
    :param model1:
    :param model2:
    :return:
    """
    rng = np.random.RandomState(42)

    # Test for weighting of ones vs None

    est_ones = clone(estimator).fit(X, sample_weight=np.ones(len(X)))
    est_none = clone(estimator).fit(X, sample_weight=None)

    assert_almost_equal(v_measure_score(est_ones.labels_, est_none.labels_), 1.0)

    # Test that random weightings produces a different model to None

    est_rand = clone(estimator).fit(X, sample_weight=rng.randint(1, 5, len(X)))

    assert_true(np.any(np.not_equal(est_rand.labels_, est_none.labels_)))

    # Test that scalar multiples of sample_weights are the same
    int_array = np.ones(len(X)) * 5
    const_k = 3
    est_s = clone(estimator).fit(X, sample_weight=int_array)
    est_sxk = clone(estimator).fit(X, sample_weight=int_array * const_k)

    assert_almost_equal(v_measure_score(est_s.labels_, est_sxk.labels_), 1.0)

    return True


def check_sample_weight_invariance(data_args, to_fit, is_equal=False):
    """
    Parameters
    ----------

    data_args : dict
        Keyword arguments to pass to fit, and which would need to be repeated
        to test equivalence to integer sample weights.
    fit : callable
        e.g. estimator
        Passed data args, returns a model that can be compared with is_equal
    is_equal : callable
        Passed two models returned from fit, returns a bool to indicate equality
        between models
    """

    estimator = to_fit(**data_args)

    is_equal(clone(estimator),clone(estimator))



def test_weighted_vs_repeated():
    # a sample weight of N should yield the same result as an N-fold
    # repetition of the sample
    r = np.random.RandomState(42)
    sample_weight = r.randint(1, 5, size=n_samples)
    X_repeat = np.repeat(X, sample_weight, axis=0)
    estimators = [KMeans(init="k-means++", n_clusters=n_clusters,
                         random_state=42),
                  KMeans(init="random", n_clusters=n_clusters,
                         random_state=42),
                  KMeans(init=centers.copy(), n_clusters=n_clusters,
                         random_state=42),
                  MiniBatchKMeans(n_clusters=n_clusters, batch_size=25, tol=1e-9,
                                  random_state=42)]
    for estimator in estimators:
        # print(estimator)
        est_weighted = clone(estimator).fit(X, sample_weight=sample_weight)
        # print(est_weighted)
        est_repeated = clone(estimator).fit(X_repeat)
        # print(est_repeated)
        repeated_labels = np.repeat(est_weighted.labels_, sample_weight)
        # print(repeated_labels)
        assert_almost_equal(v_measure_score(est_repeated.labels_,
                                            repeated_labels), 1.0)
        if not isinstance(estimator, MiniBatchKMeans):
            assert_almost_equal(_sort_centers(est_weighted.cluster_centers_),
                                _sort_centers(est_repeated.cluster_centers_))
