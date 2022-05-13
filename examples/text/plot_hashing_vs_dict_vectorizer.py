"""
===========================================
FeatureHasher and DictVectorizer Comparison
===========================================

Compares :func:`~sklearn.feature_extraction.FeatureHasher` and
:func:`~sklearn.feature_extraction.DictVectorizer` by using both to vectorize
text documents.

The example demonstrates syntax and speed only; it doesn't actually do anything
useful with the extracted vectors. See the example scripts
:ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py` and
:ref:`sphx_glr_auto_examples_text_plot_document_clustering.py` for actual
learning on text documents.

A discrepancy between the number of terms reported for
:func:`~sklearn.feature_extraction.DictVectorizer` and for
:func:`~sklearn.feature_extraction.FeatureHasher` is to be expected due to hash
collisions.

"""

# Author: Lars Buitinck
# License: BSD 3 clause

# %%
# Load Data
# ---------
# We load data from the :ref:`20newsgroups_dataset`, which comprises around
# 18000 newsgroups posts on 20 topics split in two subsets: one for training and
# one for testing. For the sake of simplicity and for reducing the computational
# cost, here we select a subset of 7 topics and use the training set only.

from sklearn.datasets import fetch_20newsgroups

categories = [
    "alt.atheism",
    "comp.graphics",
    "comp.sys.ibm.pc.hardware",
    "misc.forsale",
    "rec.autos",
    "sci.space",
    "talk.religion.misc",
]

print("Loading 20 newsgroups training data")
raw_data, _ = fetch_20newsgroups(subset="train", categories=categories, return_X_y=True)
data_size_mb = sum(len(s.encode("utf-8")) for s in raw_data) / 1e6
print("%d documents - %0.3fMB" % (len(raw_data), data_size_mb))
print()

# %%
# Define preprocessing functions
# ------------------------------
# This functions are later used for passing either a list of lists or a
# dictionary to the vectorizers.

import re
from collections import defaultdict


def tokens(doc):
    """Extract tokens from doc.

    This uses a simple regex to break strings into tokens. For a more
    principled approach, see CountVectorizer or TfidfVectorizer.
    """
    return (tok.lower() for tok in re.findall(r"\w+", doc))


def token_freqs(doc):
    """Extract a dict mapping tokens from doc to their frequencies."""
    freq = defaultdict(int)
    for tok in tokens(doc):
        freq[tok] += 1
    return freq


# %%
# Benchmarking
# ------------
# First we benchmark the `DictVectorizer`, then we compare it
# to `CountVectorizer`.

import matplotlib.pyplot as plt
from time import time
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer

dict_count_vectorizers = defaultdict(list)

print("DictVectorizer")
t0 = time()
vectorizer = DictVectorizer()
vectorizer.fit_transform(token_freqs(d) for d in raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique tokens" % len(vectorizer.get_feature_names_out()))
print()

print("CountVectorizer")
t0 = time()
vectorizer = CountVectorizer()
vectorizer.fit_transform(raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique terms" % len(vectorizer.get_feature_names_out()))
print()

print("TfidfVectorizer")
t0 = time()
vectorizer = TfidfVectorizer()
vectorizer.fit_transform(raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique terms" % len(vectorizer.get_feature_names_out()))
print()

# %%
# We plot the speed of both methods for vectorizing.

import numpy as np

fig, ax = plt.subplots()

y_pos = np.arange(len(dict_count_vectorizers["vectorizer"]))
ax.barh(y_pos, dict_count_vectorizers["speed"], align="center")
ax.set_yticks(y_pos, labels=dict_count_vectorizers["vectorizer"])
ax.invert_yaxis()
_ = ax.set_xlabel("speed (MB/s)")

# %%
# We see that using the `CountVectorizer` implementation is approximately twice
# as fast as using the `DictVectorizer` along with the simple function we
# defined for mapping the tokens.
#
# Now we make a similar experiment with `FeatureHasher` and the `HashingVectorizer`.
# In this case we need an additional function to count active features.


def n_nonzero_columns(X):
    """Number of columns with at least one non-zero value in a CSR matrix.

    This is useful to count the number of features columns that are effectively
    active when using the FeatureHasher.
    """
    return len(np.unique(X.nonzero()[1]))


# %%
# The default number of features for the `FeatureHasher` is 2**20. Here we set
# `n_features = 2**18` to illustrate hash collisions. Notice that a large number
# of features causes larger coefficient dimensions and then requires more
# training time.

from sklearn.feature_extraction import FeatureHasher
from sklearn.feature_extraction.text import HashingVectorizer

hasher_hashing_vectorizer = defaultdict(list)

print("FeatureHasher on frequency dicts")
t0 = time()
hasher = FeatureHasher(n_features=2**18)
X = hasher.transform(token_freqs(d) for d in raw_data)
duration = time() - t0
hasher_hashing_vectorizer["vectorizer"].append(hasher.__class__.__name__)
hasher_hashing_vectorizer["speed"].append(data_size_mb / duration)
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique tokens" % n_nonzero_columns(X))
print()

print("FeatureHasher on raw tokens")
t0 = time()
hasher = FeatureHasher(n_features=2**18, input_type="string")
X = hasher.transform(tokens(d) for d in raw_data)
duration = time() - t0
hasher_hashing_vectorizer["vectorizer"].append(hasher.__class__.__name__)
hasher_hashing_vectorizer["speed"].append(data_size_mb / duration)
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique tokens" % n_nonzero_columns(X))
print()

print("HashingVectorizer")
t0 = time()
vectorizer = HashingVectorizer()
vectorizer.fit_transform(raw_data)
duration = time() - t0
hasher_hashing_vectorizer["vectorizer"].append(vectorizer.__class__.__name__)
hasher_hashing_vectorizer["speed"].append(data_size_mb / duration)
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print()

fig, ax = plt.subplots()

y_pos = np.arange(len(hasher_hashing_vectorizer["vectorizer"]))
ax.barh(y_pos, hasher_hashing_vectorizer["speed"], align="center")
ax.set_yticks(y_pos, labels=hasher_hashing_vectorizer["vectorizer"])
ax.invert_yaxis()
_ = ax.set_xlabel("speed (MB/s)")

# %%
# The number of unique tokens when using the `FeatureHasher` is lower than those
# obtained using the `DictVectorizer` and `CountVectorizer`. This is due to hash
# collisions.
