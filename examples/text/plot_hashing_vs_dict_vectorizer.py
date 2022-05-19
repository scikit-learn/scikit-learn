"""
===========================================
FeatureHasher and DictVectorizer Comparison
===========================================

In this example we first compare
:func:`~sklearn.feature_extraction.FeatureHasher` and
:func:`~sklearn.feature_extraction.DictVectorizer` by using both to vectorize
text documents. Then we introduce and further analyze the text specific
vectorizers :func:`~sklearn.feature_extraction.text.HashingVectorizer`,
:func:`~sklearn.feature_extraction.text.CountVectorizer` and
:func:`~sklearn.feature_extraction.text.TfidfVectorizer`.

The objective of the example is to demonstrate syntax and speed only; it doesn't
actually do anything useful with the extracted vectors. See the example scripts
:ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
and :ref:`sphx_glr_auto_examples_text_plot_document_clustering.py` for actual
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
#
# We load data from the :ref:`20newsgroups_dataset`, which comprises around
# 18000 newsgroups posts on 20 topics split in two subsets: one for training and
# one for testing. For the sake of simplicity and reducing the computational
# cost, we select a subset of 7 topics and use the training set only.

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
print(f"{len(raw_data)} documents - {data_size_mb:.3f}MB")

# %%
# Define preprocessing functions
# ------------------------------
#
# A token may be a word, part of a word or anything comprised between spaces or
# symbols in a string. Here we define a function that extracts the tokens using
# a simple regular expression (regex) that matches Unicode word characters. This
# includes most characters that can be part of a word in any language, as well
# as numbers and the underscore:

import re


def tokens(doc):
    """Extract tokens from doc.

    This uses a simple regex that matches word characters to break strings
    into tokens. For a more principled approach, see CountVectorizer or
    TfidfVectorizer.
    """
    return (tok.lower() for tok in re.findall(r"\w+", doc))


list(tokens("This is a simple example, isn't it?"))

# %%
# We define an additional function that transforms a document into a frequency
# dictionary to be used by the vectorizers.

from collections import defaultdict


def token_freqs(doc):
    """Extract a dict mapping tokens from doc to their frequencies.

    Count the number of occurrences of all tokens of each document.
    """

    freq = defaultdict(int)
    for tok in tokens(doc):
        freq[tok] += 1
    return freq


token_freqs("That is one example, but this is another one")

# %%
# DictVectorizer
# --------------
#
# First we benchmark the `DictVectorizer`, then we compare it to `FeatureHasher`
# as both of them receive dictionaries as input.

from time import time
from sklearn.feature_extraction import DictVectorizer

dict_count_vectorizers = defaultdict(list)

t0 = time()
vectorizer = DictVectorizer()
vectorizer.fit_transform(token_freqs(d) for d in raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print(f"done in {duration:.3f}s at {data_size_mb / duration:.3f}MB/s")
print(f"Found {len(vectorizer.get_feature_names_out())} unique terms")

# %%
# FeatureHasher
# -------------
#
# Dictionaries take up a large amount of storage space and grow in size as the
# training set grows. Instead of growing the vectors along with a dictionary,
# feature hashing builds a vector of pre-defined length by applying a hash
# function `h` to the features (e.g., tokens), then using the hash values
# directly as feature indices and updating the resulting vector at those
# indices. Hashing functions tend to map distinct values to the same hash code
# (hash collisions). As a result, it is impossible to determine what object
# generated any particular hash code.
#
# Because of the above it is impossible to recover the original tokens from the
# feature matrix and the best approach to estimate the number of unique terms in
# the original dictionary is to count the number of active columns in the
# encoded feature matrix. For such purpose we define the following function:

import numpy as np


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
#
# **FeatureHasher on frequency dictionaries**

from sklearn.feature_extraction import FeatureHasher

t0 = time()
hasher = FeatureHasher(n_features=2**18)
X = hasher.transform(token_freqs(d) for d in raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(
    hasher.__class__.__name__ + "\non freq dicts"
)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print(f"done in {duration:.3f}s at {data_size_mb / duration:.3f}MB/s")
print(f"Found {n_nonzero_columns(X)} unique tokens")

# %%
# **FeatureHasher on raw tokens**

t0 = time()
hasher = FeatureHasher(n_features=2**18, input_type="string")
X = hasher.transform(tokens(d) for d in raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(
    hasher.__class__.__name__ + "\non raw tokens"
)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print(f"done in {duration:.3f}s at {data_size_mb / duration:.3f}MB/s")
print(f"Found {n_nonzero_columns(X)} unique tokens")

# %%
# The number of unique tokens when using the `FeatureHasher` is lower than those
# obtained using the `DictVectorizer`. This is due to hash collisions.
#
# We now plot the speed of the above methods for vectorizing.

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 6))

y_pos = np.arange(len(dict_count_vectorizers["vectorizer"]))
ax.barh(y_pos, dict_count_vectorizers["speed"], align="center")
ax.set_yticks(y_pos)
ax.set_yticklabels(dict_count_vectorizers["vectorizer"])
ax.invert_yaxis()
_ = ax.set_xlabel("speed (MB/s)")

# %%
# In both cases `FeatureHasher` is approximately twice as fast as
# `DictVectorizer`. This is handy when dealing with large amounts of data, with
# the downside of losing the invertibility of the transformation, which in turn
# makes the interpretation of a model a more complex task.
#
# Comparison with special purpose text vectorizers
# ------------------------------------------------
#
# `CountVectorizer` accepts raw data as it internally implements tokenization
# and occurrence counting. It is similar to the `DictVectorizer` when used along
# with the customized function `token_freqs` as done in the previous section.
# The difference being that `CountVectorizer` is more flexible. In particular it
# accepts various regex patterns through the `token_pattern` parameter.

from sklearn.feature_extraction.text import CountVectorizer

t0 = time()
vectorizer = CountVectorizer()
vectorizer.fit_transform(raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print(f"done in {duration:.3f}s at {data_size_mb / duration:.3f}MB/s")
print(f"Found {len(vectorizer.get_feature_names_out())} unique terms")

# %%
# We see that using the `CountVectorizer` implementation is approximately twice
# as fast as using the `DictVectorizer` along with the simple function we
# defined for mapping the tokens. The reason is that `CountVectorizer` produces
# a sparse representation of the counts using `scipy.sparse.csr_matrix`.
#
# Now we make a similar experiment with the `HashingVectorizer`, which is
# equivalent to combining the “hashing trick” implemented by the `FeatureHasher`
# class and the text preprocessing and tokenization of the `CountVectorizer`.

from sklearn.feature_extraction.text import HashingVectorizer

t0 = time()
vectorizer = HashingVectorizer(n_features=2**18)
vectorizer.fit_transform(raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print(f"done in {duration:.3f}s at {data_size_mb / duration:.3f}MB/s")

# %%
# TfidfVectorizer
# ---------------
#
# In a large text corpus, some words appear with higher frequency (e.g. “the”,
# “a”, “is” in English) which do not carry meaningful information about the
# actual contents of a document. If we were to feed the word count data directly
# to a classifier those very common terms would shadow the frequencies of rarer
# yet more interesting terms. In order to re-weight the count features into
# floating point values suitable for usage by a classifier it is very common to
# use the tf–idf transform as implemented by the `TfidfTransformer`. Tf stands
# for "term-frequency" while "tf–idf" means term-frequency times inverse
# document-frequency.
#
# We now benchmark the `TfidfVectorizer`, which is equivalent to combining the
# tokenization and occurrence counting of the `CountVectorizer` along with the
# normalizing and weighting from a `TfidfTransformer`.

from sklearn.feature_extraction.text import TfidfVectorizer

t0 = time()
vectorizer = TfidfVectorizer()
vectorizer.fit_transform(raw_data)
duration = time() - t0
dict_count_vectorizers["vectorizer"].append(vectorizer.__class__.__name__)
dict_count_vectorizers["speed"].append(data_size_mb / duration)
print(f"done in {duration:.3f}s at {data_size_mb / duration:.3f}MB/s")
print(f"Found {len(vectorizer.get_feature_names_out())} unique terms")

fig, ax = plt.subplots(figsize=(12, 6))

y_pos = np.arange(len(dict_count_vectorizers["vectorizer"]))
ax.barh(y_pos, dict_count_vectorizers["speed"], align="center")
ax.set_yticks(y_pos)
ax.set_yticklabels(dict_count_vectorizers["vectorizer"])
ax.invert_yaxis()
_ = ax.set_xlabel("speed (MB/s)")

# %%
# Notice from the plot that `TfidfVectorizer` is slightly slower than
# `CountVectorizer` because of the extra operation induced by the
# `TfidfTransformer`. Also notice that, by setting the number of features
# `n_features = 2**18`, the `HashingVectorizer` performs much faster than the
# `CountVectorizer`. Keep in mind that this may cause hash collisions and no
# inverse transform can be used to recover the original dictionary.
