from __future__ import division, print_function

from functools import partial
from itertools import product

import numpy as np
import scipy.sparse as sp

from sklearn.datasets import make_multilabel_classification
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils.multiclass import type_of_target
from sklearn.utils.validation import check_random_state
from sklearn.utils import shuffle

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import _named_check

from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import brier_score_loss
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import coverage_error
from sklearn.metrics import explained_variance_score
from sklearn.metrics import f1_score
from sklearn.metrics import fbeta_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import hinge_loss
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics import label_ranking_average_precision_score
from sklearn.metrics import label_ranking_loss
from sklearn.metrics import log_loss
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import median_absolute_error
from sklearn.metrics import precision_score
from sklearn.metrics import r2_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import zero_one_loss

# TODO Curve are currently not covered by invariance test
# from sklearn.metrics import precision_recall_curve
# from sklearn.metrics import roc_curve


from sklearn.metrics.base import _average_binary_score


# Note toward developers about metric testing
# -------------------------------------------
# It is often possible to write one general test for several metrics:
#
#   - invariance properties, e.g. invariance to sample order
#   - common behavior for an argument, e.g. the "normalize" with value True
#     will return the mean of the metrics and with value False will return
#     the sum of the metrics.
#
# In order to improve the overall metric testing, it is a good idea to write
# first a specific test for the given metric and then add a general test for
# all metrics that have the same behavior.
#
# Two types of datastructures are used in order to implement this system:
# dictionaries of metrics and lists of metrics wit common properties.
#
# Dictionaries of metrics
# ------------------------
# The goal of having those dictionaries is to have an easy way to call a
# particular metric and associate a name to each function:
#
#   - REGRESSION_METRICS: all regression metrics.
#   - CLASSIFICATION_METRICS: all classification metrics
#     which compare a ground truth and the estimated targets as returned by a
#     classifier.
#   - THRESHOLDED_METRICS: all classification metrics which
#     compare a ground truth and a score, e.g. estimated probabilities or
#     decision function (format might vary)
#
# Those dictionaries will be used to test systematically some invariance
# properties, e.g. invariance toward several input layout.
#

REGRESSION_METRICS = {
    "mean_absolute_error": mean_absolute_error,
    "mean_squared_error": mean_squared_error,
    "median_absolute_error": median_absolute_error,
    "explained_variance_score": explained_variance_score,
    "r2_score": partial(r2_score, multioutput='variance_weighted'),
}

CLASSIFICATION_METRICS = {
    "accuracy_score": accuracy_score,
    "unnormalized_accuracy_score": partial(accuracy_score, normalize=False),
    "confusion_matrix": confusion_matrix,
    "hamming_loss": hamming_loss,

    "jaccard_similarity_score": jaccard_similarity_score,
    "unnormalized_jaccard_similarity_score":
    partial(jaccard_similarity_score, normalize=False),

    "zero_one_loss": zero_one_loss,
    "unnormalized_zero_one_loss": partial(zero_one_loss, normalize=False),

    # These are needed to test averaging
    "precision_score": precision_score,
    "recall_score": recall_score,
    "f1_score": f1_score,
    "f2_score": partial(fbeta_score, beta=2),
    "f0.5_score": partial(fbeta_score, beta=0.5),
    "matthews_corrcoef_score": matthews_corrcoef,

    "weighted_f0.5_score": partial(fbeta_score, average="weighted", beta=0.5),
    "weighted_f1_score": partial(f1_score, average="weighted"),
    "weighted_f2_score": partial(fbeta_score, average="weighted", beta=2),
    "weighted_precision_score": partial(precision_score, average="weighted"),
    "weighted_recall_score": partial(recall_score, average="weighted"),

    "micro_f0.5_score": partial(fbeta_score, average="micro", beta=0.5),
    "micro_f1_score": partial(f1_score, average="micro"),
    "micro_f2_score": partial(fbeta_score, average="micro", beta=2),
    "micro_precision_score": partial(precision_score, average="micro"),
    "micro_recall_score": partial(recall_score, average="micro"),

    "macro_f0.5_score": partial(fbeta_score, average="macro", beta=0.5),
    "macro_f1_score": partial(f1_score, average="macro"),
    "macro_f2_score": partial(fbeta_score, average="macro", beta=2),
    "macro_precision_score": partial(precision_score, average="macro"),
    "macro_recall_score": partial(recall_score, average="macro"),

    "samples_f0.5_score": partial(fbeta_score, average="samples", beta=0.5),
    "samples_f1_score": partial(f1_score, average="samples"),
    "samples_f2_score": partial(fbeta_score, average="samples", beta=2),
    "samples_precision_score": partial(precision_score, average="samples"),
    "samples_recall_score": partial(recall_score, average="samples"),

    "cohen_kappa_score": cohen_kappa_score,
}

THRESHOLDED_METRICS = {
    "coverage_error": coverage_error,
    "label_ranking_loss": label_ranking_loss,
    "log_loss": log_loss,
    "unnormalized_log_loss": partial(log_loss, normalize=False),

    "hinge_loss": hinge_loss,

    "brier_score_loss": brier_score_loss,

    "roc_auc_score": roc_auc_score,
    "weighted_roc_auc": partial(roc_auc_score, average="weighted"),
    "samples_roc_auc": partial(roc_auc_score, average="samples"),
    "micro_roc_auc": partial(roc_auc_score, average="micro"),
    "macro_roc_auc": partial(roc_auc_score, average="macro"),

    "average_precision_score": average_precision_score,
    "weighted_average_precision_score":
    partial(average_precision_score, average="weighted"),
    "samples_average_precision_score":
    partial(average_precision_score, average="samples"),
    "micro_average_precision_score":
    partial(average_precision_score, average="micro"),
    "macro_average_precision_score":
    partial(average_precision_score, average="macro"),
    "label_ranking_average_precision_score":
    label_ranking_average_precision_score,
}

ALL_METRICS = dict()
ALL_METRICS.update(THRESHOLDED_METRICS)
ALL_METRICS.update(CLASSIFICATION_METRICS)
ALL_METRICS.update(REGRESSION_METRICS)

# Lists of metrics with common properties
# ---------------------------------------
# Lists of metrics with common properties are used to test systematically some
# functionalities and invariance, e.g. SYMMETRIC_METRICS lists all metrics that
# are symmetric with respect to their input argument y_true and y_pred.
#
# When you add a new metric or functionality, check if a general test
# is already written.

# Those metrics don't support binary inputs
METRIC_UNDEFINED_BINARY = [
    "samples_f0.5_score",
    "samples_f1_score",
    "samples_f2_score",
    "samples_precision_score",
    "samples_recall_score",
    "coverage_error",

    "roc_auc_score",
    "micro_roc_auc",
    "weighted_roc_auc",
    "macro_roc_auc",
    "samples_roc_auc",

    "average_precision_score",
    "weighted_average_precision_score",
    "micro_average_precision_score",
    "macro_average_precision_score",
    "samples_average_precision_score",

    "label_ranking_loss",
    "label_ranking_average_precision_score",
]

# Those metrics don't support multiclass inputs
METRIC_UNDEFINED_MULTICLASS = [
    "brier_score_loss",

    # with default average='binary', multiclass is prohibited
    "precision_score",
    "recall_score",
    "f1_score",
    "f2_score",
    "f0.5_score",
]

# Metric undefined with "binary" or "multiclass" input
METRIC_UNDEFINED_BINARY_MULTICLASS = set(METRIC_UNDEFINED_BINARY).union(
    set(METRIC_UNDEFINED_MULTICLASS))

# Metrics with an "average" argument
METRICS_WITH_AVERAGING = [
    "precision_score", "recall_score", "f1_score", "f2_score", "f0.5_score"
]

# Threshold-based metrics with an "average" argument
THRESHOLDED_METRICS_WITH_AVERAGING = [
    "roc_auc_score", "average_precision_score",
]

# Metrics with a "pos_label" argument
METRICS_WITH_POS_LABEL = [
    "roc_curve",

    "brier_score_loss",

    "precision_score", "recall_score", "f1_score", "f2_score", "f0.5_score",

    # pos_label support deprecated; to be removed in 0.18:
    "weighted_f0.5_score", "weighted_f1_score", "weighted_f2_score",
    "weighted_precision_score", "weighted_recall_score",

    "micro_f0.5_score", "micro_f1_score", "micro_f2_score",
    "micro_precision_score", "micro_recall_score",

    "macro_f0.5_score", "macro_f1_score", "macro_f2_score",
    "macro_precision_score", "macro_recall_score",
]

# Metrics with a "labels" argument
# TODO: Handle multi_class metrics that has a labels argument as well as a
# decision function argument. e.g hinge_loss
METRICS_WITH_LABELS = [
    "confusion_matrix",

    "hamming_loss",

    "precision_score", "recall_score", "f1_score", "f2_score", "f0.5_score",

    "weighted_f0.5_score", "weighted_f1_score", "weighted_f2_score",
    "weighted_precision_score", "weighted_recall_score",

    "micro_f0.5_score", "micro_f1_score", "micro_f2_score",
    "micro_precision_score", "micro_recall_score",

    "macro_f0.5_score", "macro_f1_score", "macro_f2_score",
    "macro_precision_score", "macro_recall_score",

    "cohen_kappa_score",
]

# Metrics with a "normalize" option
METRICS_WITH_NORMALIZE_OPTION = [
    "accuracy_score",
    "jaccard_similarity_score",
    "zero_one_loss",
]

# Threshold-based metrics with "multilabel-indicator" format support
THRESHOLDED_MULTILABEL_METRICS = [
    "log_loss",
    "unnormalized_log_loss",

    "roc_auc_score", "weighted_roc_auc", "samples_roc_auc",
    "micro_roc_auc", "macro_roc_auc",

    "average_precision_score", "weighted_average_precision_score",
    "samples_average_precision_score", "micro_average_precision_score",
    "macro_average_precision_score",

    "coverage_error", "label_ranking_loss",
]

# Classification metrics with  "multilabel-indicator" format
MULTILABELS_METRICS = [
    "accuracy_score", "unnormalized_accuracy_score",
    "hamming_loss",
    "jaccard_similarity_score", "unnormalized_jaccard_similarity_score",
    "zero_one_loss", "unnormalized_zero_one_loss",

    "weighted_f0.5_score", "weighted_f1_score", "weighted_f2_score",
    "weighted_precision_score", "weighted_recall_score",

    "macro_f0.5_score", "macro_f1_score", "macro_f2_score",
    "macro_precision_score", "macro_recall_score",

    "micro_f0.5_score", "micro_f1_score", "micro_f2_score",
    "micro_precision_score", "micro_recall_score",

    "samples_f0.5_score", "samples_f1_score", "samples_f2_score",
    "samples_precision_score", "samples_recall_score",
]

# Regression metrics with "multioutput-continuous" format support
MULTIOUTPUT_METRICS = [
    "mean_absolute_error", "mean_squared_error", "r2_score",
    "explained_variance_score"
]

# Symmetric with respect to their input arguments y_true and y_pred
# metric(y_true, y_pred) == metric(y_pred, y_true).
SYMMETRIC_METRICS = [
    "accuracy_score", "unnormalized_accuracy_score",
    "hamming_loss",
    "jaccard_similarity_score", "unnormalized_jaccard_similarity_score",
    "zero_one_loss", "unnormalized_zero_one_loss",

    "f1_score", "micro_f1_score", "macro_f1_score",
    "weighted_recall_score",
    # P = R = F = accuracy in multiclass case
    "micro_f0.5_score", "micro_f1_score", "micro_f2_score",
    "micro_precision_score", "micro_recall_score",

    "matthews_corrcoef_score", "mean_absolute_error", "mean_squared_error",
    "median_absolute_error",

    "cohen_kappa_score",
]

# Asymmetric with respect to their input arguments y_true and y_pred
# metric(y_true, y_pred) != metric(y_pred, y_true).
NOT_SYMMETRIC_METRICS = [
    "explained_variance_score",
    "r2_score",
    "confusion_matrix",

    "precision_score", "recall_score", "f2_score", "f0.5_score",

    "weighted_f0.5_score", "weighted_f1_score", "weighted_f2_score",
    "weighted_precision_score",

    "macro_f0.5_score", "macro_f2_score", "macro_precision_score",
    "macro_recall_score", "log_loss", "hinge_loss"
]


# No Sample weight support
METRICS_WITHOUT_SAMPLE_WEIGHT = [
    "confusion_matrix", # Left this one here because the tests in this file do
                        # not work for confusion_matrix, as its output is a
                        # matrix instead of a number. Testing of
                        # confusion_matrix with sample_weight is in
                        # test_classification.py
    "median_absolute_error",
]


@ignore_warnings
def test_symmetry():
    # Test the symmetry of score and loss functions
    random_state = check_random_state(0)
    y_true = random_state.randint(0, 2, size=(20, ))
    y_pred = random_state.randint(0, 2, size=(20, ))

    # We shouldn't forget any metrics
    assert_equal(set(SYMMETRIC_METRICS).union(
        NOT_SYMMETRIC_METRICS, THRESHOLDED_METRICS,
        METRIC_UNDEFINED_BINARY_MULTICLASS),
        set(ALL_METRICS))

    assert_equal(
        set(SYMMETRIC_METRICS).intersection(set(NOT_SYMMETRIC_METRICS)),
        set([]))

    # Symmetric metric
    for name in SYMMETRIC_METRICS:
        metric = ALL_METRICS[name]
        assert_almost_equal(metric(y_true, y_pred),
                            metric(y_pred, y_true),
                            err_msg="%s is not symmetric" % name)

    # Not symmetric metrics
    for name in NOT_SYMMETRIC_METRICS:
        metric = ALL_METRICS[name]
        assert_true(np.any(metric(y_true, y_pred) != metric(y_pred, y_true)),
                    msg="%s seems to be symmetric" % name)


@ignore_warnings
def test_sample_order_invariance():
    random_state = check_random_state(0)
    y_true = random_state.randint(0, 2, size=(20, ))
    y_pred = random_state.randint(0, 2, size=(20, ))
    y_true_shuffle, y_pred_shuffle = shuffle(y_true, y_pred, random_state=0)

    for name, metric in ALL_METRICS.items():
        if name in METRIC_UNDEFINED_BINARY_MULTICLASS:
            continue

        assert_almost_equal(metric(y_true, y_pred),
                            metric(y_true_shuffle, y_pred_shuffle),
                            err_msg="%s is not sample order invariant"
                                    % name)


@ignore_warnings
def test_sample_order_invariance_multilabel_and_multioutput():
    random_state = check_random_state(0)

    # Generate some data
    y_true = random_state.randint(0, 2, size=(20, 25))
    y_pred = random_state.randint(0, 2, size=(20, 25))
    y_score = random_state.normal(size=y_true.shape)

    y_true_shuffle, y_pred_shuffle, y_score_shuffle = shuffle(y_true,
                                                              y_pred,
                                                              y_score,
                                                              random_state=0)

    for name in MULTILABELS_METRICS:
        metric = ALL_METRICS[name]
        assert_almost_equal(metric(y_true, y_pred),
                            metric(y_true_shuffle, y_pred_shuffle),
                            err_msg="%s is not sample order invariant"
                                    % name)

    for name in THRESHOLDED_MULTILABEL_METRICS:
        metric = ALL_METRICS[name]
        assert_almost_equal(metric(y_true, y_score),
                            metric(y_true_shuffle, y_score_shuffle),
                            err_msg="%s is not sample order invariant"
                                    % name)

    for name in MULTIOUTPUT_METRICS:
        metric = ALL_METRICS[name]
        assert_almost_equal(metric(y_true, y_score),
                            metric(y_true_shuffle, y_score_shuffle),
                            err_msg="%s is not sample order invariant"
                                    % name)
        assert_almost_equal(metric(y_true, y_pred),
                            metric(y_true_shuffle, y_pred_shuffle),
                            err_msg="%s is not sample order invariant"
                                    % name)


@ignore_warnings
def test_format_invariance_with_1d_vectors():
    random_state = check_random_state(0)
    y1 = random_state.randint(0, 2, size=(20, ))
    y2 = random_state.randint(0, 2, size=(20, ))

    y1_list = list(y1)
    y2_list = list(y2)

    y1_1d, y2_1d = np.array(y1), np.array(y2)
    assert_equal(y1_1d.ndim, 1)
    assert_equal(y2_1d.ndim, 1)
    y1_column = np.reshape(y1_1d, (-1, 1))
    y2_column = np.reshape(y2_1d, (-1, 1))
    y1_row = np.reshape(y1_1d, (1, -1))
    y2_row = np.reshape(y2_1d, (1, -1))

    for name, metric in ALL_METRICS.items():
        if name in METRIC_UNDEFINED_BINARY_MULTICLASS:
            continue

        measure = metric(y1, y2)

        assert_almost_equal(metric(y1_list, y2_list), measure,
                            err_msg="%s is not representation invariant "
                                    "with list" % name)

        assert_almost_equal(metric(y1_1d, y2_1d), measure,
                            err_msg="%s is not representation invariant "
                                    "with np-array-1d" % name)

        assert_almost_equal(metric(y1_column, y2_column), measure,
                            err_msg="%s is not representation invariant "
                                    "with np-array-column" % name)

        # Mix format support
        assert_almost_equal(metric(y1_1d, y2_list), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and list" % name)

        assert_almost_equal(metric(y1_list, y2_1d), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and list" % name)

        assert_almost_equal(metric(y1_1d, y2_column), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and np-array-column"
                                    % name)

        assert_almost_equal(metric(y1_column, y2_1d), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and np-array-column"
                                    % name)

        assert_almost_equal(metric(y1_list, y2_column), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix list and np-array-column"
                                    % name)

        assert_almost_equal(metric(y1_column, y2_list), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix list and np-array-column"
                                    % name)

        # These mix representations aren't allowed
        assert_raises(ValueError, metric, y1_1d, y2_row)
        assert_raises(ValueError, metric, y1_row, y2_1d)
        assert_raises(ValueError, metric, y1_list, y2_row)
        assert_raises(ValueError, metric, y1_row, y2_list)
        assert_raises(ValueError, metric, y1_column, y2_row)
        assert_raises(ValueError, metric, y1_row, y2_column)

        # NB: We do not test for y1_row, y2_row as these may be
        # interpreted as multilabel or multioutput data.
        if (name not in (MULTIOUTPUT_METRICS + THRESHOLDED_MULTILABEL_METRICS +
                         MULTILABELS_METRICS)):
            assert_raises(ValueError, metric, y1_row, y2_row)


@ignore_warnings
def test_invariance_string_vs_numbers_labels():
    # Ensure that classification metrics with string labels
    random_state = check_random_state(0)
    y1 = random_state.randint(0, 2, size=(20, ))
    y2 = random_state.randint(0, 2, size=(20, ))

    y1_str = np.array(["eggs", "spam"])[y1]
    y2_str = np.array(["eggs", "spam"])[y2]

    pos_label_str = "spam"
    labels_str = ["eggs", "spam"]

    for name, metric in CLASSIFICATION_METRICS.items():
        if name in METRIC_UNDEFINED_BINARY_MULTICLASS:
            continue

        measure_with_number = metric(y1, y2)

        # Ugly, but handle case with a pos_label and label
        metric_str = metric
        if name in METRICS_WITH_POS_LABEL:
            metric_str = partial(metric_str, pos_label=pos_label_str)

        measure_with_str = metric_str(y1_str, y2_str)

        assert_array_equal(measure_with_number, measure_with_str,
                           err_msg="{0} failed string vs number invariance "
                                   "test".format(name))

        measure_with_strobj = metric_str(y1_str.astype('O'),
                                         y2_str.astype('O'))
        assert_array_equal(measure_with_number, measure_with_strobj,
                           err_msg="{0} failed string object vs number "
                                   "invariance test".format(name))

        if name in METRICS_WITH_LABELS:
            metric_str = partial(metric_str, labels=labels_str)
            measure_with_str = metric_str(y1_str, y2_str)
            assert_array_equal(measure_with_number, measure_with_str,
                               err_msg="{0} failed string vs number  "
                                       "invariance test".format(name))

            measure_with_strobj = metric_str(y1_str.astype('O'),
                                             y2_str.astype('O'))
            assert_array_equal(measure_with_number, measure_with_strobj,
                               err_msg="{0} failed string vs number  "
                                       "invariance test".format(name))

    for name, metric in THRESHOLDED_METRICS.items():
        if name in ("log_loss", "hinge_loss", "unnormalized_log_loss",
                    "brier_score_loss"):
            # Ugly, but handle case with a pos_label and label
            metric_str = metric
            if name in METRICS_WITH_POS_LABEL:
                metric_str = partial(metric_str, pos_label=pos_label_str)

            measure_with_number = metric(y1, y2)
            measure_with_str = metric_str(y1_str, y2)
            assert_array_equal(measure_with_number, measure_with_str,
                               err_msg="{0} failed string vs number "
                                       "invariance test".format(name))

            measure_with_strobj = metric(y1_str.astype('O'), y2)
            assert_array_equal(measure_with_number, measure_with_strobj,
                               err_msg="{0} failed string object vs number "
                                       "invariance test".format(name))
        else:
            # TODO those metrics doesn't support string label yet
            assert_raises(ValueError, metric, y1_str, y2)
            assert_raises(ValueError, metric, y1_str.astype('O'), y2)


def test_inf_nan_input():
    invalids =[([0, 1], [np.inf, np.inf]),
               ([0, 1], [np.nan, np.nan]),
               ([0, 1], [np.nan, np.inf])]

    METRICS = dict()
    METRICS.update(THRESHOLDED_METRICS)
    METRICS.update(REGRESSION_METRICS)

    for metric in METRICS.values():
        for y_true, y_score in invalids:
            assert_raise_message(ValueError,
                                 "contains NaN, infinity",
                                 metric, y_true, y_score)

    # Classification metrics all raise a mixed input exception
    for metric in CLASSIFICATION_METRICS.values():
        for y_true, y_score in invalids:
            assert_raise_message(ValueError,
                                 "Classification metrics can't handle a mix "
                                 "of binary and continuous targets",
                                 metric, y_true, y_score)


@ignore_warnings
def check_single_sample(name):
    # Non-regression test: scores should work with a single sample.
    # This is important for leave-one-out cross validation.
    # Score functions tested are those that formerly called np.squeeze,
    # which turns an array of size 1 into a 0-d array (!).
    metric = ALL_METRICS[name]

    # assert that no exception is thrown
    for i, j in product([0, 1], repeat=2):
        metric([i], [j])


@ignore_warnings
def check_single_sample_multioutput(name):
    metric = ALL_METRICS[name]
    for i, j, k, l in product([0, 1], repeat=4):
        metric(np.array([[i, j]]), np.array([[k, l]]))


def test_single_sample():
    for name in ALL_METRICS:
        if (name in METRIC_UNDEFINED_BINARY_MULTICLASS or
                name in THRESHOLDED_METRICS):
            # Those metrics are not always defined with one sample
            # or in multiclass classification
            continue

        yield check_single_sample, name

    for name in MULTIOUTPUT_METRICS + MULTILABELS_METRICS:
        yield check_single_sample_multioutput, name


def test_multioutput_number_of_output_differ():
    y_true = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
    y_pred = np.array([[0, 0], [1, 0], [0, 0]])

    for name in MULTIOUTPUT_METRICS:
        metric = ALL_METRICS[name]
        assert_raises(ValueError, metric, y_true, y_pred)


def test_multioutput_regression_invariance_to_dimension_shuffling():
    # test invariance to dimension shuffling
    random_state = check_random_state(0)
    y_true = random_state.uniform(0, 2, size=(20, 5))
    y_pred = random_state.uniform(0, 2, size=(20, 5))

    for name in MULTIOUTPUT_METRICS:
        metric = ALL_METRICS[name]
        error = metric(y_true, y_pred)

        for _ in range(3):
            perm = random_state.permutation(y_true.shape[1])
            assert_almost_equal(metric(y_true[:, perm], y_pred[:, perm]),
                                error,
                                err_msg="%s is not dimension shuffling "
                                        "invariant" % name)


@ignore_warnings
def test_multilabel_representation_invariance():
    # Generate some data
    n_classes = 4
    n_samples = 50

    _, y1 = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                           random_state=0, n_samples=n_samples,
                                           allow_unlabeled=True)
    _, y2 = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                           random_state=1, n_samples=n_samples,
                                           allow_unlabeled=True)

    # To make sure at least one empty label is present
    y1 = np.vstack([y1, [[0] * n_classes]])
    y2 = np.vstack([y2, [[0] * n_classes]])

    y1_sparse_indicator = sp.coo_matrix(y1)
    y2_sparse_indicator = sp.coo_matrix(y2)

    for name in MULTILABELS_METRICS:
        metric = ALL_METRICS[name]

        # XXX cruel hack to work with partial functions
        if isinstance(metric, partial):
            metric.__module__ = 'tmp'
            metric.__name__ = name

        measure = metric(y1, y2)

        # Check representation invariance
        assert_almost_equal(metric(y1_sparse_indicator,
                                   y2_sparse_indicator),
                            measure,
                            err_msg="%s failed representation invariance  "
                                    "between dense and sparse indicator "
                                    "formats." % name)


def test_raise_value_error_multilabel_sequences():
    # make sure the multilabel-sequence format raises ValueError
    multilabel_sequences = [
        [[0, 1]],
        [[1], [2], [0, 1]],
        [(), (2), (0, 1)],
        [[]],
        [()],
        np.array([[], [1, 2]], dtype='object')]

    for name in MULTILABELS_METRICS:
        metric = ALL_METRICS[name]
        for seq in multilabel_sequences:
            assert_raises(ValueError, metric, seq, seq)


def test_normalize_option_binary_classification(n_samples=20):
    # Test in the binary case
    random_state = check_random_state(0)
    y_true = random_state.randint(0, 2, size=(n_samples, ))
    y_pred = random_state.randint(0, 2, size=(n_samples, ))

    for name in METRICS_WITH_NORMALIZE_OPTION:
        metrics = ALL_METRICS[name]
        measure = metrics(y_true, y_pred, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true, y_pred, normalize=False)
                            / n_samples, measure)


def test_normalize_option_multiclass_classification():
    # Test in the multiclass case
    random_state = check_random_state(0)
    y_true = random_state.randint(0, 4, size=(20, ))
    y_pred = random_state.randint(0, 4, size=(20, ))
    n_samples = y_true.shape[0]

    for name in METRICS_WITH_NORMALIZE_OPTION:
        metrics = ALL_METRICS[name]
        measure = metrics(y_true, y_pred, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true, y_pred, normalize=False)
                            / n_samples, measure)


def test_normalize_option_multilabel_classification():
    # Test in the multilabel case
    n_classes = 4
    n_samples = 100

    # for both random_state 0 and 1, y_true and y_pred has at least one
    # unlabelled entry
    _, y_true = make_multilabel_classification(n_features=1,
                                               n_classes=n_classes,
                                               random_state=0,
                                               allow_unlabeled=True,
                                               n_samples=n_samples)
    _, y_pred = make_multilabel_classification(n_features=1,
                                               n_classes=n_classes,
                                               random_state=1,
                                               allow_unlabeled=True,
                                               n_samples=n_samples)

    # To make sure at least one empty label is present
    y_true += [0]*n_classes
    y_pred += [0]*n_classes

    for name in METRICS_WITH_NORMALIZE_OPTION:
        metrics = ALL_METRICS[name]
        measure = metrics(y_true, y_pred, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true, y_pred, normalize=False)
                            / n_samples, measure,
                            err_msg="Failed with %s" % name)


@ignore_warnings
def _check_averaging(metric, y_true, y_pred, y_true_binarize, y_pred_binarize,
                     is_multilabel):
    n_samples, n_classes = y_true_binarize.shape

    # No averaging
    label_measure = metric(y_true, y_pred, average=None)
    assert_array_almost_equal(label_measure,
                              [metric(y_true_binarize[:, i],
                                      y_pred_binarize[:, i])
                               for i in range(n_classes)])

    # Micro measure
    micro_measure = metric(y_true, y_pred, average="micro")
    assert_almost_equal(micro_measure, metric(y_true_binarize.ravel(),
                                              y_pred_binarize.ravel()))

    # Macro measure
    macro_measure = metric(y_true, y_pred, average="macro")
    assert_almost_equal(macro_measure, np.mean(label_measure))

    # Weighted measure
    weights = np.sum(y_true_binarize, axis=0, dtype=int)

    if np.sum(weights) != 0:
        weighted_measure = metric(y_true, y_pred, average="weighted")
        assert_almost_equal(weighted_measure, np.average(label_measure,
                                                         weights=weights))
    else:
        weighted_measure = metric(y_true, y_pred, average="weighted")
        assert_almost_equal(weighted_measure, 0)

    # Sample measure
    if is_multilabel:
        sample_measure = metric(y_true, y_pred, average="samples")
        assert_almost_equal(sample_measure,
                            np.mean([metric(y_true_binarize[i],
                                            y_pred_binarize[i])
                                     for i in range(n_samples)]))

    assert_raises(ValueError, metric, y_true, y_pred, average="unknown")
    assert_raises(ValueError, metric, y_true, y_pred, average="garbage")


def check_averaging(name, y_true, y_true_binarize, y_pred, y_pred_binarize,
                    y_score):
    is_multilabel = type_of_target(y_true).startswith("multilabel")

    metric = ALL_METRICS[name]

    if name in METRICS_WITH_AVERAGING:
        _check_averaging(metric, y_true, y_pred, y_true_binarize,
                         y_pred_binarize, is_multilabel)
    elif name in THRESHOLDED_METRICS_WITH_AVERAGING:
        _check_averaging(metric, y_true, y_score, y_true_binarize,
                         y_score, is_multilabel)
    else:
        raise ValueError("Metric is not recorded as having an average option")


def test_averaging_multiclass(n_samples=50, n_classes=3):
    random_state = check_random_state(0)
    y_true = random_state.randint(0, n_classes, size=(n_samples, ))
    y_pred = random_state.randint(0, n_classes, size=(n_samples, ))
    y_score = random_state.uniform(size=(n_samples, n_classes))

    lb = LabelBinarizer().fit(y_true)
    y_true_binarize = lb.transform(y_true)
    y_pred_binarize = lb.transform(y_pred)

    for name in METRICS_WITH_AVERAGING:
        yield (_named_check(check_averaging, name), name, y_true,
               y_true_binarize, y_pred, y_pred_binarize, y_score)


def test_averaging_multilabel(n_classes=5, n_samples=40):
    _, y = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                          random_state=5, n_samples=n_samples,
                                          allow_unlabeled=False)
    y_true = y[:20]
    y_pred = y[20:]
    y_score = check_random_state(0).normal(size=(20, n_classes))
    y_true_binarize = y_true
    y_pred_binarize = y_pred

    for name in METRICS_WITH_AVERAGING + THRESHOLDED_METRICS_WITH_AVERAGING:
        yield (_named_check(check_averaging, name), name, y_true,
               y_true_binarize, y_pred, y_pred_binarize, y_score)


def test_averaging_multilabel_all_zeroes():
    y_true = np.zeros((20, 3))
    y_pred = np.zeros((20, 3))
    y_score = np.zeros((20, 3))
    y_true_binarize = y_true
    y_pred_binarize = y_pred

    for name in METRICS_WITH_AVERAGING:
        yield (_named_check(check_averaging, name), name, y_true,
               y_true_binarize, y_pred, y_pred_binarize, y_score)

    # Test _average_binary_score for weight.sum() == 0
    binary_metric = (lambda y_true, y_score, average="macro":
                     _average_binary_score(
                         precision_score, y_true, y_score, average))
    _check_averaging(binary_metric, y_true, y_pred, y_true_binarize,
                     y_pred_binarize, is_multilabel=True)


def test_averaging_multilabel_all_ones():
    y_true = np.ones((20, 3))
    y_pred = np.ones((20, 3))
    y_score = np.ones((20, 3))
    y_true_binarize = y_true
    y_pred_binarize = y_pred

    for name in METRICS_WITH_AVERAGING:
        yield (_named_check(check_averaging, name), name, y_true,
               y_true_binarize, y_pred, y_pred_binarize, y_score)


@ignore_warnings
def check_sample_weight_invariance(name, metric, y1, y2):
    rng = np.random.RandomState(0)
    sample_weight = rng.randint(1, 10, size=len(y1))

    # check that unit weights gives the same score as no weight
    unweighted_score = metric(y1, y2, sample_weight=None)
    assert_almost_equal(
        unweighted_score,
        metric(y1, y2, sample_weight=np.ones(shape=len(y1))),
        err_msg="For %s sample_weight=None is not equivalent to "
                "sample_weight=ones" % name)

    # check that the weighted and unweighted scores are unequal
    weighted_score = metric(y1, y2, sample_weight=sample_weight)
    assert_not_equal(
        unweighted_score, weighted_score,
        msg="Unweighted and weighted scores are unexpectedly "
            "equal (%f) for %s" % (weighted_score, name))

    # check that sample_weight can be a list
    weighted_score_list = metric(y1, y2,
                                 sample_weight=sample_weight.tolist())
    assert_almost_equal(
        weighted_score, weighted_score_list,
        err_msg=("Weighted scores for array and list "
                 "sample_weight input are not equal (%f != %f) for %s") % (
                     weighted_score, weighted_score_list, name))

    # check that integer weights is the same as repeated samples
    repeat_weighted_score = metric(
        np.repeat(y1, sample_weight, axis=0),
        np.repeat(y2, sample_weight, axis=0), sample_weight=None)
    assert_almost_equal(
        weighted_score, repeat_weighted_score,
        err_msg="Weighting %s is not equal to repeating samples" % name)

    # check that ignoring a fraction of the samples is equivalent to setting
    # the corresponding weights to zero
    sample_weight_subset = sample_weight[1::2]
    sample_weight_zeroed = np.copy(sample_weight)
    sample_weight_zeroed[::2] = 0
    y1_subset = y1[1::2]
    y2_subset = y2[1::2]
    weighted_score_subset = metric(y1_subset, y2_subset,
                                   sample_weight=sample_weight_subset)
    weighted_score_zeroed = metric(y1, y2,
                                   sample_weight=sample_weight_zeroed)
    assert_almost_equal(
        weighted_score_subset, weighted_score_zeroed,
        err_msg=("Zeroing weights does not give the same result as "
                 "removing the corresponding samples (%f != %f) for %s" %
                 (weighted_score_zeroed, weighted_score_subset, name)))

    if not name.startswith('unnormalized'):
        # check that the score is invariant under scaling of the weights by a
        # common factor
        for scaling in [2, 0.3]:
            assert_almost_equal(
                weighted_score,
                metric(y1, y2, sample_weight=sample_weight * scaling),
                err_msg="%s sample_weight is not invariant "
                        "under scaling" % name)

    # Check that if sample_weight.shape[0] != y_true.shape[0], it raised an
    # error
    assert_raises(Exception, metric, y1, y2,
                  sample_weight=np.hstack([sample_weight, sample_weight]))


def test_sample_weight_invariance(n_samples=50):
    random_state = check_random_state(0)
    # regression
    y_true = random_state.random_sample(size=(n_samples,))
    y_pred = random_state.random_sample(size=(n_samples,))
    for name in ALL_METRICS:
        if name not in REGRESSION_METRICS:
            continue
        if name in METRICS_WITHOUT_SAMPLE_WEIGHT:
            continue
        metric = ALL_METRICS[name]
        yield _named_check(check_sample_weight_invariance, name), name,\
            metric, y_true, y_pred

    # binary
    random_state = check_random_state(0)
    y_true = random_state.randint(0, 2, size=(n_samples, ))
    y_pred = random_state.randint(0, 2, size=(n_samples, ))
    y_score = random_state.random_sample(size=(n_samples,))
    for name in ALL_METRICS:
        if name in REGRESSION_METRICS:
            continue
        if (name in METRICS_WITHOUT_SAMPLE_WEIGHT or
                name in METRIC_UNDEFINED_BINARY):
            continue
        metric = ALL_METRICS[name]
        if name in THRESHOLDED_METRICS:
            yield _named_check(check_sample_weight_invariance, name), name,\
                  metric, y_true, y_score
        else:
            yield _named_check(check_sample_weight_invariance, name), name,\
                  metric, y_true, y_pred

    # multiclass
    random_state = check_random_state(0)
    y_true = random_state.randint(0, 5, size=(n_samples, ))
    y_pred = random_state.randint(0, 5, size=(n_samples, ))
    y_score = random_state.random_sample(size=(n_samples, 5))
    for name in ALL_METRICS:
        if name in REGRESSION_METRICS:
            continue
        if (name in METRICS_WITHOUT_SAMPLE_WEIGHT or
                name in METRIC_UNDEFINED_BINARY_MULTICLASS):
            continue
        metric = ALL_METRICS[name]
        if name in THRESHOLDED_METRICS:
            yield _named_check(check_sample_weight_invariance, name), name,\
                  metric, y_true, y_score
        else:
            yield _named_check(check_sample_weight_invariance, name), name,\
                  metric, y_true, y_pred

    # multilabel indicator
    _, ya = make_multilabel_classification(n_features=1, n_classes=20,
                                           random_state=0, n_samples=100,
                                           allow_unlabeled=False)
    _, yb = make_multilabel_classification(n_features=1, n_classes=20,
                                           random_state=1, n_samples=100,
                                           allow_unlabeled=False)
    y_true = np.vstack([ya, yb])
    y_pred = np.vstack([ya, ya])
    y_score = random_state.randint(1, 4, size=y_true.shape)

    for name in (MULTILABELS_METRICS + THRESHOLDED_MULTILABEL_METRICS +
                 MULTIOUTPUT_METRICS):
        if name in METRICS_WITHOUT_SAMPLE_WEIGHT:
            continue

        metric = ALL_METRICS[name]
        if name in THRESHOLDED_METRICS:
            yield (_named_check(check_sample_weight_invariance, name), name,
                   metric, y_true, y_score)
        else:
            yield (_named_check(check_sample_weight_invariance, name), name,
                   metric, y_true, y_pred)


@ignore_warnings
def test_no_averaging_labels():
    # test labels argument when not using averaging
    # in multi-class and multi-label cases
    y_true_multilabel = np.array([[1, 1, 0, 0], [1, 1, 0, 0]])
    y_pred_multilabel = np.array([[0, 0, 1, 1], [0, 1, 1, 0]])
    y_true_multiclass = np.array([0, 1, 2])
    y_pred_multiclass = np.array([0, 2, 3])
    labels = np.array([3, 0, 1, 2])
    _, inverse_labels = np.unique(labels, return_inverse=True)

    for name in METRICS_WITH_AVERAGING:
        for y_true, y_pred in [[y_true_multiclass, y_pred_multiclass],
                               [y_true_multilabel, y_pred_multilabel]]:
            if name not in MULTILABELS_METRICS and y_pred.ndim > 1:
                continue

            metric = ALL_METRICS[name]

            score_labels = metric(y_true, y_pred, labels=labels, average=None)
            score = metric(y_true, y_pred, average=None)
            assert_array_equal(score_labels, score[inverse_labels])
