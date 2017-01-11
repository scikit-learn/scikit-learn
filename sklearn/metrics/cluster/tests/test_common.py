from functools import partial

import numpy as np

from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import mutual_info_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import v_measure_score
from sklearn.metrics.cluster import silhouette_score
from sklearn.metrics.cluster import calinski_harabaz_score

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_contains
from sklearn.utils.testing import assert_almost_equal


# Dictionaries of metrics
# ------------------------
# The goal of having those dictionaries is to have an easy way to call a
# particular metric and associate a name to each function:
#   - SUPERVISED_METRICS: all supervised cluster metrics - (when given a
# ground truth value)
#   - UNSUPERVISED_METRICS: all unsupervised cluster metrics
#
# Those dictionaries will be used to test systematically some invariance
# properties, e.g. invariance toward several input layout.
#

# Metrics used to test similarity between bicluster
SUPERVISED_METRICS = {
    "adjusted_mutual_info_score": adjusted_mutual_info_score,
    "adjusted_rand_score": adjusted_rand_score,
    "completeness_score": completeness_score,
    "homogeneity_score": homogeneity_score,
    "mutual_info_score": mutual_info_score,
    "normalized_mutual_info_score": normalized_mutual_info_score,
    "v_measure_score": v_measure_score,
    "fowlkes_mallows_score": fowlkes_mallows_score
}

UNSUPERVISED_METRICS = {
    "silhouette_score": silhouette_score,
    "silhouette_manhattan": partial(silhouette_score, metric='manhattan'),
    "calinski_harabaz_score": calinski_harabaz_score
}

# Lists of metrics with common properties
# ---------------------------------------
# Lists of metrics with common properties are used to test systematically some
# functionalities and invariance, e.g. SYMMETRIC_METRICS lists all metrics
# that are symmetric with respect to their input argument y_true and y_pred.
#
# --------------------------------------------------------------------
# Symmetric with respect to their input arguments y_true and y_pred.
# Symmetric metrics only apply to supervised clusters.
SYMMETRIC_METRICS = [
    "adjusted_rand_score", "v_measure_score",
    "mutual_info_score", "adjusted_mutual_info_score",
    "normalized_mutual_info_score", "fowlkes_mallows_score"
]

NON_SYMMETRIC_METRICS = ["homogeneity_score", "completeness_score"]


# Metrics whose upper bound is 1
NORMALIZED_METRICS = [
    "adjusted_rand_score", "homogeneity_score", "completeness_score",
    "v_measure_score", "adjusted_mutual_info_score", "fowlkes_mallows_score",
    "normalized_mutual_info_score"
]


def test_symmetry():
    rng = np.random.RandomState(0)
    y1 = rng.randint(3, size=30)
    y2 = rng.randint(3, size=30)
    for name in SYMMETRIC_METRICS:
        metric = SUPERVISED_METRICS[name]
        assert_almost_equal(metric(y1, y2), metric(y2, y1),
                            err_msg="%s is not symmetric" % name)

    for name in NON_SYMMETRIC_METRICS:
        metric = SUPERVISED_METRICS[name]
        assert_not_equal(metric(y1, y2), metric(y2, y1),
                         msg="%s is symmetric" % name)


def test_normalized_output():
    upper_bound_1 = [0, 0, 0, 1, 1, 1]
    upper_bound_2 = [0, 0, 0, 1, 1, 1]
    for name in NORMALIZED_METRICS:
        metric = SUPERVISED_METRICS[name]
        assert_greater(metric([0, 0, 0, 1, 1], [0, 0, 0, 1, 2]), 0.0)
        assert_greater(metric([0, 0, 1, 1, 2], [0, 0, 1, 1, 1]), 0.0)
        assert_less(metric([0, 0, 0, 1, 2], [0, 1, 1, 1, 1]), 1.0)
        assert_less(metric([0, 0, 0, 1, 2], [0, 1, 1, 1, 1]), 1.0)
        assert_equal(metric(upper_bound_1, upper_bound_2), 1.0,
                     msg="%s has upper_bound greater than 1" % name)

    lower_bound_1 = [0, 0, 0, 0, 0, 0]
    lower_bound_2 = [0, 1, 2, 3, 4, 5]
    for name in NORMALIZED_METRICS:
        metric = SUPERVISED_METRICS[name]
        assert_contains(0, [metric(lower_bound_1, lower_bound_2), 
                            metric(lower_bound_2, lower_bound_2)])


# All clustering metrics do not change score due to permutations of labels
# that is when 0 and 1 exchchanged.
def test_permute_labels():
    y_label = np.array([0, 0, 0, 1, 1, 0, 1])
    y_pred = np.array([1, 0, 1, 0, 1, 1, 0])
    for name in SUPERVISED_METRICS:
        metric = SUPERVISED_METRICS[name]
        score_1 = metric(y_pred, y_label)
        assert_equal(score_1, metric(1 - y_pred, y_label),
                     msg="%s failed labels permutation" % name)
        assert_equal(score_1, metric(1 - y_pred, 1 - y_label),
                     msg="%s failed labels permutation" % name)
        assert_almost_equal(score_1, metric(y_pred, 1 - y_label),
                            err_msg="%s failed labels permutation" % name)

    for name in UNSUPERVISED_METRICS:
        metric = UNSUPERVISED_METRICS[name]
        X = np.random.randint(10, size=(7, 10))
        score_1 = metric(X, y_pred)
        assert_almost_equal(score_1, metric(X, 1 - y_pred),
                            err_msg="%s failed labels permutation" % name)


# For ALL clustering metrics Input parameters can be both
# in the form of arrays lists , positive , negetive or string
def test_format_invariance():
    labels = [0, 0, 0, 0, 1, 1, 1, 1]
    y_pred = [0, 1, 2, 3, 4, 5, 6, 7]

    def generate_formats(y):
        y = np.array(y)
        y = np.array(y)
        yield y, 'array of ints'
        yield y.tolist(), 'list of ints'
        yield [str(x) for x in y.tolist()], 'list of strs'
        yield y - 1, 'including negative ints'
        yield y + 1, 'strictly positive ints'

    for name in SUPERVISED_METRICS:
        metric = SUPERVISED_METRICS[name]
        score_1 = metric(labels, y_pred)
        mygenerator_1 = generate_formats(labels)
        mygenerator_2 = generate_formats(y_pred)
        for i, j in zip(mygenerator_1, mygenerator_2):
            assert_equal(score_1, metric(i[0], j[0]),
                         msg="%s failed %s format invariance" % (name, j[1]))

    for name in UNSUPERVISED_METRICS:
        metric = UNSUPERVISED_METRICS[name]
        X = np.random.randint(10, size=(8, 10))
        score_1 = metric(X, labels)
        assert_equal(score_1, metric(X.astype(float), labels))
        mygenerator_1 = generate_formats(labels)
        for j in mygenerator_1:
            assert_equal(score_1, metric(X, j[0]),
                         msg="%s failed %s format invariance" % (name, j[1]))
