"""
===========================================
FeatureHasher and DictVectorizer Comparison
===========================================

Compares FeatureHasher and DictVectorizer by using both to vectorize text
documents.

The example demonstrates syntax and speed only; it doesn't actually do anything
useful with the extracted vectors. See the example scripts
:ref:`sphx_glr_auto_examples_text_plot_document_classification.py` and
:ref:`sphx_glr_auto_examples_text_plot_document_clustering.py` for actual
learning on text documents.

A discrepancy between the number of terms reported for DictVectorizer and for
FeatureHasher is to be expected due to hash collisions.

"""

# Author: Lars Buitinck
# License: BSD 3 clause

# %%
import sys

print(__doc__)
print("Usage: %s [n_features_for_hashing]" % sys.argv[0])

# %%
# Load Data
# ---------

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
# Uncomment the following line to use a larger set (11k+ documents)
# categories = None

print("Loading 20 newsgroups training data")
raw_data, _ = fetch_20newsgroups(subset="train", categories=categories, return_X_y=True)
data_size_mb = sum(len(s.encode("utf-8")) for s in raw_data) / 1e6
print("%d documents - %0.3fMB" % (len(raw_data), data_size_mb))
print()

# %%
# Define preprocessing functions
# ------------------------------

import re
import numpy as np
from collections import defaultdict


def n_nonzero_columns(X):
    """Returns the number of non-zero columns in a CSR matrix X."""
    return len(np.unique(X.nonzero()[1]))


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

from time import time
from sklearn.feature_extraction import DictVectorizer, FeatureHasher

print("DictVectorizer")
t0 = time()
vectorizer = DictVectorizer()
vectorizer.fit_transform(token_freqs(d) for d in raw_data)
duration = time() - t0
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique terms" % len(vectorizer.get_feature_names_out()))
print()

# The default number of features is 2**20.
try:
    n_features = int(sys.argv[1])
except IndexError:
    n_features = 2**18
except ValueError:
    n_features = 2**18

print("FeatureHasher on frequency dicts")
t0 = time()
hasher = FeatureHasher(n_features=n_features)
X = hasher.transform(token_freqs(d) for d in raw_data)
duration = time() - t0
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique terms" % n_nonzero_columns(X))
print()

print("FeatureHasher on raw tokens")
t0 = time()
hasher = FeatureHasher(n_features=n_features, input_type="string")
X = hasher.transform(tokens(d) for d in raw_data)
duration = time() - t0
print("done in %fs at %0.3fMB/s" % (duration, data_size_mb / duration))
print("Found %d unique terms" % n_nonzero_columns(X))
