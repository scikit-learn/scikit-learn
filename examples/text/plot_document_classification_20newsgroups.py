"""
======================================================
Classification of text documents using sparse features
======================================================

This is an example showing how scikit-learn can be used to classify documents by
topics using a bag-of-words approach. This example uses a Tf-idf-weighted
document-term sparse matrix to encode the features and demonstrates various
classifiers that can efficiently handle sparse matrices.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Mathieu Blondel <mathieu@mblondel.org>
#         Lars Buitinck
# License: BSD 3 clause


# %%
# Load data
# ---------
# We define a function to load data from :ref:`20newsgroups_dataset`, which
# comprises around 18000 newsgroups posts on 20 topics split in two subsets: one
# for training (or development) and the other one for testing (or for
# performance evaluation). Note that, by default, the text samples contain
# some message metadata such as 'headers', 'footers' (signatures) and 'quotes'
# to other posts. The `fetch_20newsgroups` function therefore accepts a
# parameter named `remove` to make it possible to attempt stripping those
# extra information that can make the classification problem "too easy". This
# is achieved using simple heuristics that are neither perfect nor standard, hence
# disabled by default.

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from time import time

categories = [
    "alt.atheism",
    "talk.religion.misc",
    "comp.graphics",
    "sci.space",
]


def size_mb(docs):
    return sum(len(s.encode("utf-8")) for s in docs) / 1e6


def load_dataset(verbose=False, remove=()):
    """Load and vectorize the 20 newsgroups dataset."""

    data_train = fetch_20newsgroups(
        subset="train",
        categories=categories,
        shuffle=True,
        random_state=42,
        remove=remove,
    )

    data_test = fetch_20newsgroups(
        subset="test",
        categories=categories,
        shuffle=True,
        random_state=42,
        remove=remove,
    )

    # order of labels in `target_names` can be different from `categories`
    target_names = data_train.target_names

    # split target in a training set and a test set
    y_train, y_test = data_train.target, data_test.target

    # Extracting features from the training data using a sparse vectorizer
    t0 = time()
    vectorizer = TfidfVectorizer(
        sublinear_tf=True, max_df=0.5, min_df=5, stop_words="english"
    )
    X_train = vectorizer.fit_transform(data_train.data)
    duration_train = time() - t0

    # Extracting features from the test data using the same vectorizer
    t0 = time()
    X_test = vectorizer.transform(data_test.data)
    duration_test = time() - t0

    feature_names = vectorizer.get_feature_names_out()

    if verbose:

        # compute size of loaded data
        data_train_size_mb = size_mb(data_train.data)
        data_test_size_mb = size_mb(data_test.data)

        print(
            "%d documents - %0.3fMB (training set)"
            % (len(data_train.data), data_train_size_mb)
        )
        print(
            "%d documents - %0.3fMB (test set)"
            % (len(data_test.data), data_test_size_mb)
        )
        print("%d categories" % len(target_names))
        print(
            "vectorize training done in %fs at %0.3fMB/s"
            % (duration_train, data_train_size_mb / duration_train)
        )
        print("n_samples: %d, n_features: %d" % X_train.shape)
        print(
            "vectorize testing done in %fs at %0.3fMB/s"
            % (duration_test, data_test_size_mb / duration_test)
        )
        print("n_samples: %d, n_features: %d" % X_test.shape)

    return X_train, X_test, y_train, y_test, feature_names, target_names


# %%
# Compare feature effects
# -----------------------
# We train a first classification model without attempting to strip the metadata of
# the dataset.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import ConfusionMatrixDisplay

X_train, X_test, y_train, y_test, feature_names, target_names = load_dataset(
    verbose=True
)

# %%
# Our first model is an instance of the
# :class:`~sklearn.linear_model.RidgeClassifier` class. This is a
# linear classification model that uses the mean squared error on
# {-1, 1} encoded targets, one for each possible class. Contrary to
# :class:`~sklearn.linear_model.LogisticRegression`, `RidgeClassifier`
# does not provide probabilistic predictions (no `predict_proba` method),
# but it is often faster to train.
clf = RidgeClassifier(tol=1e-2, solver="sparse_cg")
clf.fit(X_train, y_train)
pred = clf.predict(X_test)

fig, ax = plt.subplots(figsize=(10, 5))
ConfusionMatrixDisplay.from_predictions(y_test, pred, ax=ax)
ax.xaxis.set_ticklabels(target_names)
ax.yaxis.set_ticklabels(target_names)
_ = ax.set_title("Confusion Matrix with headers for the\n" + str(clf).split("(")[0])


def plot_feature_effects():
    # learned coefficients weighted by frequency of appearance
    feature_effects = clf.coef_ * np.asarray(X_train.mean(axis=0)).ravel()

    for i, label in enumerate(target_names):
        top5 = np.argsort(feature_effects[i])[-5:][::-1]
        if i == 0:
            top = pd.DataFrame(feature_names[top5], columns=[label])
            top_indices = top5
        else:
            top[label] = feature_names[top5]
            top_indices = np.concatenate((top_indices, top5), axis=None)
    top_indices = np.unique(top_indices)
    predictive_words = feature_names[top_indices]

    # plot feature effects
    bar_size = 0.25
    padding = 0.75
    y_locs = np.arange(len(top_indices)) * (4 * bar_size + padding)

    fig, ax = plt.subplots(figsize=(10, 8))
    for i, label in enumerate(target_names):
        ax.barh(
            y_locs + (i - 2) * bar_size,
            feature_effects[i, top_indices],
            height=bar_size,
            label=label,
        )
    ax.set(
        yticks=y_locs,
        yticklabels=predictive_words,
        ylim=[
            0 - 4 * bar_size,
            len(top_indices) * (4 * bar_size + padding) - 4 * bar_size,
        ],
    )
    ax.legend(loc="lower right")

    print("top 5 keywords per class:")
    print(top)

    return ax


_ = plot_feature_effects().set_title("Feature impact with headers")

# %%
# The word "caltech" is one of the top predictive features for atheism due to
# pollution in the dataset introduced by the metadata
data_train = fetch_20newsgroups(
    subset="train", categories=categories, shuffle=True, random_state=42
)

for doc in data_train.data:
    if "caltech" in doc:
        print(doc)
        break

# %%
# We can remove the metadata, but other sources of pollution may still exist
(
    X_train,
    X_test,
    y_train,
    y_test,
    feature_names,
    target_names,
) = load_dataset(remove=("headers", "footers", "quotes"))

clf = RidgeClassifier(tol=1e-2, solver="sparse_cg")
clf.fit(X_train, y_train)
pred = clf.predict(X_test)

fig, ax = plt.subplots(figsize=(10, 5))
ConfusionMatrixDisplay.from_predictions(y_test, pred, ax=ax)
ax.xaxis.set_ticklabels(target_names)
ax.yaxis.set_ticklabels(target_names)
_ = ax.set_title("Confusion Matrix without headers for the\n" + str(clf).split("(")[0])

# %%
_ = plot_feature_effects().set_title("Feature impact without headers")

# %%
# By looking at the confusion matrix, it is more evident that the scores of the
# model trained with metadata were overoptimistic. The classification problem
# without access to the metadata is harder (more errors) but more representative
# of the intended text classification problem.

# %%
# Benchmark classifiers
# ---------------------
#
# First we define small benchmarking utilities
from sklearn.utils.extmath import density


def trim(s):
    """Trim string to fit on terminal (assuming 80-column display)"""
    return s if len(s) <= 80 else s[:77] + "..."


def benchmark(clf, custom_name=False):
    print("_" * 80)
    print("Training: ")
    print(clf)
    t0 = time()
    clf.fit(X_train, y_train)
    train_time = time() - t0
    print("train time: %0.3fs" % train_time)

    t0 = time()
    pred = clf.predict(X_test)
    test_time = time() - t0
    print("test time:  %0.3fs" % test_time)

    score = metrics.accuracy_score(y_test, pred)
    print("accuracy:   %0.3f" % score)

    if hasattr(clf, "coef_"):
        print("dimensionality: %d" % clf.coef_.shape[1])
        print("density: %f" % density(clf.coef_))
        print("top 10 keywords per class:")
        for i, label in enumerate(target_names):
            top10 = np.argsort(clf.coef_[i])[-10:]
            print(trim("%s: %s" % (label, " ".join(feature_names[top10]))))
        print()

    print("classification report:")
    print(metrics.classification_report(y_test, pred, target_names=target_names))

    print("confusion matrix:")
    print(metrics.confusion_matrix(y_test, pred))

    print()
    if custom_name:
        clf_descr = str(custom_name)
    else:
        clf_descr = str(clf).split("(")[0]
    return clf_descr, score, train_time, test_time


# %%
# We now train and test the datasets with 8 different classification
# models and get performance results for each model. The goal of this study
# is to highlight the computation/accuracy tradeoffs of the choice of the
# type of classifier for such a multi-class text classification problem.
#
# Note that the most important hyper-parameters values where tuned with
# a simple grid search procedure not shown in this notebook for the sake
# of simplicity. 
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import ComplementNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import RandomForestClassifier


results = []
for clf, name in (
    (LogisticRegression(C=5, max_iter=1000), "Logistic Regression"),
    (RidgeClassifier(alpha=1.0, solver="sparse_cg"), "Ridge Classifier"),
    (KNeighborsClassifier(n_neighbors=100), "kNN"),
    (RandomForestClassifier(), "Random Forest"),
):
    print("=" * 80)
    print(name)
    results.append(benchmark(clf, name))

# Train Liblinear model
print("=" * 80)
print("L2 penalty Linear SVC")
results.append(benchmark(LinearSVC(C=0.1, dual=False, max_iter=1000)))

# Train SGD model
print("=" * 80)
print("L2 penalty Linear SGD")
results.append(
    benchmark(
        SGDClassifier(loss="log", alpha=1e-4, n_iter_no_change=3, early_stopping=True)
    )
)

# Train NearestCentroid without threshold
print("=" * 80)
print("NearestCentroid (aka Rocchio classifier)")
results.append(benchmark(NearestCentroid()))

# Train sparse Naive Bayes classifiers
print("=" * 80)
print("Naive Bayes")
results.append(benchmark(ComplementNB(alpha=0.1)))

# %%
# Plot accuracy, training and test time of each classifier
# --------------------------------------------------------
# The scatter plots show the trade-off between the test accuracy and the
# training and testing time of each classifier.

indices = np.arange(len(results))

results = [[x[i] for x in results] for i in range(4)]

clf_names, score, training_time, test_time = results
training_time = np.array(training_time)
test_time = np.array(test_time)

fig, ax1 = plt.subplots(figsize=(10, 8))
ax1.scatter(score, training_time, s=60)
ax1.set(
    title="Score-training time trade-off",
    yscale="log",
    xlabel="test accuracy",
    ylabel="training time (s)",
)
fig, ax2 = plt.subplots(figsize=(10, 8))
ax2.scatter(score, test_time, s=60)
ax2.set(
    title="Score-test time trade-off",
    yscale="log",
    xlabel="test accuracy",
    ylabel="test time (s)",
)

for i, txt in enumerate(clf_names):
    ax1.annotate(txt, (score[i], training_time[i]))
    ax2.annotate(txt, (score[i], test_time[i]))

# %%
# The Naive Bayes model has the best trade-off between score and
# training/testing time, while Random Forest is both slow to train
# expensive to predict and has a comparatively bad accuracy. This is
# expected: for high-dimensional prediction problems, linear models are
# often better suited as all problem became linearly separable as the
# dimensionality of the feature space increases to 10,000 dimensions or
# more.
# KNeighborsClassifier has relatively low accuracy and has the highest testing
# time. The long prediction time is expected: for each prediction the model
# has to compute the pairwise distances with each document in the training
# set which is very expensive.
# Furthermore, the curse of dimensionality harms the ability of this model to
# yield competitive accuracy in the high dimensional feature space of text
# classification problems.
