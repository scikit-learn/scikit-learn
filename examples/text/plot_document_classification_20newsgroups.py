"""
======================================================
Classification of text documents using sparse features
======================================================

This is an example showing how scikit-learn can be used to classify documents
by topics using a bag-of-words approach. This example uses a scipy.sparse
matrix to store the features and demonstrates various classifiers that can
efficiently handle sparse matrices.

The dataset used in this example is the 20 newsgroups dataset. It will be
automatically downloaded, then cached.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Mathieu Blondel <mathieu@mblondel.org>
#         Lars Buitinck
# License: BSD 3 clause


# %%
# Configuration options for the analysis
# --------------------------------------

# If True, we use `HashingVectorizer`, otherwise we use a `TfidfVectorizer`
USE_HASHING = False

# Number of features used by `HashingVectorizer`
N_FEATURES = 2**16

# Optional feature selection: either False, or an integer: the number of
# features to select
SELECT_CHI2 = False

# If True, remove headers, footers and quotes when fetching dataset
REMOVE_ALL = False

# %%
# Load data from the training set
# -------------------------------
# Let's load data from the newsgroups dataset which comprises around 18000
# newsgroups posts on 20 topics split in two subsets: one for training (or
# development) and the other one for testing (or for performance evaluation).
from sklearn.datasets import fetch_20newsgroups

categories = [
    "alt.atheism",
    "talk.religion.misc",
    "comp.graphics",
    "sci.space",
]

if REMOVE_ALL:
    remove = ("headers", "footers", "quotes")
else:
    remove = ()

data_train = fetch_20newsgroups(
    subset="train", categories=categories, shuffle=True, random_state=42, remove=remove
)

data_test = fetch_20newsgroups(
    subset="test", categories=categories, shuffle=True, random_state=42, remove=remove
)
print("data loaded")

# order of labels in `target_names` can be different from `categories`
target_names = data_train.target_names

# split target in a training set and a test set
y_train, y_test = data_train.target, data_test.target


def size_mb(docs):
    return sum(len(s.encode("utf-8")) for s in docs) / 1e6


data_train_size_mb = size_mb(data_train.data)
data_test_size_mb = size_mb(data_test.data)

print(
    "%d documents - %0.3fMB (training set)" % (len(data_train.data), data_train_size_mb)
)
print("%d documents - %0.3fMB (test set)" % (len(data_test.data), data_test_size_mb))
print("%d categories" % len(target_names))

# %%
# Vectorize the training and test data
# -------------------------------------
#
# Extracting features from the training data using a sparse vectorizer
from time import time

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import HashingVectorizer

t0 = time()

if USE_HASHING:
    vectorizer = HashingVectorizer(
        stop_words="english", alternate_sign=False, n_features=N_FEATURES
    )
    X_train = vectorizer.transform(data_train.data)
else:
    vectorizer = TfidfVectorizer(sublinear_tf=True, max_df=0.5, stop_words="english")
    X_train = vectorizer.fit_transform(data_train.data)
duration = time() - t0
print("done in %fs at %0.3fMB/s" % (duration, data_train_size_mb / duration))
print("n_samples: %d, n_features: %d" % X_train.shape)

# %%
# Extracting features from the test data using the same vectorizer
t0 = time()
X_test = vectorizer.transform(data_test.data)
duration = time() - t0
print("done in %fs at %0.3fMB/s" % (duration, data_test_size_mb / duration))
print("n_samples: %d, n_features: %d" % X_test.shape)

# %%
# Mapping from integer feature name to original token string
if USE_HASHING:
    feature_names = None
else:
    feature_names = vectorizer.get_feature_names_out()

# %%
# Keeping only the best features if SELECT_CHI2 == True
from sklearn.feature_selection import SelectKBest, chi2

if SELECT_CHI2:
    print("Extracting %d best features by a chi-squared test" % SELECT_CHI2)
    t0 = time()
    ch2 = SelectKBest(chi2, k=SELECT_CHI2)
    X_train = ch2.fit_transform(X_train, y_train)
    X_test = ch2.transform(X_test)
    if feature_names is not None:
        # keep selected feature names
        feature_names = feature_names[ch2.get_support()]
    print("done in %fs" % (time() - t0))
    print()

# %%
# Compare feature effects
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.linear_model import LogisticRegression

clf = LogisticRegression(C=0.2, max_iter=200)
clf.fit(X_train, y_train)
pred = clf.predict(X_test)
score = metrics.accuracy_score(y_test, pred)

if feature_names is not None:
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
    padding = 0.5
    y_locs = np.arange(len(top_indices)) * (bar_size * 3 + padding)

    fig, ax = plt.subplots(figsize=(10, 8))
    for i, label in enumerate(target_names):
        ax.barh(
            y_locs + i * bar_size,
            feature_effects[i, top_indices],
            height=bar_size,
            label=label,
        )
    ax.set(
        title="Feature impact",
        yticks=y_locs,
        yticklabels=predictive_words,
        ylim=[0 - padding, len(y_locs) + 4.5],
    )
    _ = ax.legend(loc="lower right")

    print("top 5 keywords per class:")
    print(top)

# %%
# The word "caltech" as one of the top predictive features is the result of the
# pollution in the dataset introduced by the headers. This effect can be fixed
# by setting REMOVE_ALL == True, but other sources of pollution may still exist
if feature_names is not None and not REMOVE_ALL:
    for doc in data_train.data:
        if "caltech" in doc:
            print(doc)
            break

# %%
# Setting REMOVE_ALL == True increases the overlap between classes
import seaborn as sns

cf_matrix = metrics.confusion_matrix(y_test, pred)
ax = sns.heatmap(cf_matrix / np.sum(cf_matrix), annot=True, fmt=".2%", cmap="Blues")

ax.set_xlabel("\nPredicted Category")
ax.set_ylabel("Actual Category ")
ax.xaxis.set_ticklabels(target_names)
ax.yaxis.set_ticklabels(target_names)
plt.tight_layout(rect=(0, 0, 1, 0.92))
_ = ax.set_title("Confusion Matrix for the\n" + str(clf).split("(")[0])

# %%
# Benchmark classifiers
# ---------------------
#
# First we define small benchmarking utilities
from sklearn.utils.extmath import density


def trim(s):
    """Trim string to fit on terminal (assuming 80-column display)"""
    return s if len(s) <= 80 else s[:77] + "..."


def benchmark(clf):
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

        if feature_names is not None:
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
    clf_descr = str(clf).split("(")[0]
    return clf_descr, score, train_time, test_time


# %%
# We now train and test the datasets with 9 different classification
# models and get performance results for each model.
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import RidgeClassifier
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import ComplementNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import RandomForestClassifier


results = []
for clf, name in (
    (LogisticRegression(), "Logistic Regression"),
    (RidgeClassifier(tol=1e-2, solver="sag"), "Ridge Classifier"),
    (KNeighborsClassifier(n_neighbors=10), "kNN"),
    (RandomForestClassifier(), "Random forest"),
):
    print("=" * 80)
    print(name)
    results.append(benchmark(clf))

# Train Liblinear model
results.append(benchmark(LinearSVC(C=10, dual=False, tol=1e-3)))

# Train SGD model
results.append(benchmark(SGDClassifier(alpha=1e-5, early_stopping=True)))

# Train NearestCentroid without threshold
print("=" * 80)
print("NearestCentroid (aka Rocchio classifier)")
results.append(benchmark(NearestCentroid()))

# Train sparse Naive Bayes classifiers
print("=" * 80)
print("Naive Bayes")
results.append(benchmark(ComplementNB(alpha=0.1)))

print("=" * 80)
print("LinearSVC with L1-based feature selection")
# The smaller C, the stronger the regularization.
# The more regularization, the more sparsity.
results.append(
    benchmark(
        Pipeline(
            [
                (
                    "feature_selection",
                    SelectFromModel(LinearSVC(penalty="l1", dual=False, tol=1e-3)),
                ),
                ("classification", LinearSVC(penalty="l2")),
            ]
        )
    )
)


# %%
# Plot accuracy, training and test time of each classifier
# --------------------------------------------------------
# The bar plot indicates the accuracy, training time (normalized) and test time
# (normalized) of each classifier.
import matplotlib.pyplot as plt

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
# The Naive Bayesian models appear to have the best trade-off between score and
# training/testing time. LogisticRegression and KNeighborsClassifier have
# relatively low accuracy and high training and testing time, respectively.
