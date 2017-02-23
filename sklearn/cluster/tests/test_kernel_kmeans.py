""" Testing for Kernel K means"""
import sys

import numpy as np
from contextlib import contextmanager

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_not_equal

from sklearn.cluster import KernelKMeans
from sklearn.datasets.samples_generator import make_blobs
from sklearn.externals.six.moves import cStringIO as StringIO

centers = np.array([
    [0.0, 5.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 4.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 5.0, 1.0],
])
n_samples = 100
n_clusters, n_features = centers.shape
X, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                            cluster_std=1., random_state=42)


def test_kernel_kmeans_dtype():
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 2))
    X = (X * 10).astype(np.uint8)
    clf = KernelKMeans().fit(X)
    pred_x = clf.predict(X)
    assert_array_equal(clf.labels_, pred_x)


def test_predict():
    clf = KernelKMeans(n_clusters=n_clusters, random_state=42)
    clf.fit(X)

    # sanity check: re-predict labeling for training set samples
    pred = clf.predict(X)
    assert_array_equal(pred, clf.labels_)

    # re-predict labels for training set using fit_predict
    pred = clf.fit_predict(X)
    assert_array_equal(pred, clf.labels_)

    centers_test = np.array([
        [0.0, 5.0, 0.0, 0.0, 3.0],
        [0.0, 1.0, 4.0, 1.0, 0.0],
        [4.0, 0.0, 0.0, 5.0, 1.0],
    ])
    n_samples_test = 100
    X_test, _ = make_blobs(n_samples=n_samples_test, centers=centers_test,
                           cluster_std=1., random_state=51)
    pred = clf.predict(X_test)
    assert_raises(AssertionError, assert_array_equal, pred, clf.labels_)


def test_kernels():
    def check_kernel(kernel):
        clf = KernelKMeans(n_clusters=n_clusters, kernel=kernel, gamma=0.5,
                           degree=3, coef0=2.0, random_state=42).fit(X)

        # check that clusters formed are equal to n_clusters
        labels = clf.labels_
        assert_equal(np.unique(labels).shape[0], n_clusters)

        # check error on dataset being too small
        assert_raises(ValueError, clf.fit, [[0., 1.]])

    for kernel in ['linear', 'rbf', 'polynomial', 'sigmoid']:
        yield check_kernel, kernel


def test_kernel_kmeans_n_init():
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 2))

    # two tests on bad n_init argument
    assert_raises_regex(ValueError, "n_init", KernelKMeans(n_init=0).fit, X)
    assert_raises_regex(ValueError, "n_init", KernelKMeans(n_init=-1).fit, X)


def test_kernel_kmeans_max_iter():
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 2))

    # two tests on bad max_iter argument
    assert_raises_regex(ValueError, "max_iter", KernelKMeans(max_iter=0).fit,
                        X)
    assert_raises_regex(ValueError, "max_iter", KernelKMeans(max_iter=-1).fit,
                        X)


def test_kernel_kmeans_verbose():
    clf = KernelKMeans(n_clusters=n_clusters, random_state=42, verbose=1)

    @contextmanager
    def capture_output():
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            yield sys.stdout
        finally:
            sys.stdout = old_stdout

    with capture_output() as out:
        clf.fit(X)
        output = out.getvalue().strip()
        assert_not_equal(output, '')
