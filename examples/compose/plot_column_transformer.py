"""
==================================================
Column Transformer with Heterogeneous Data Sources
==================================================

Datasets can often contain components that require different feature
extraction and processing pipelines. This scenario might occur when:

1. your dataset consists of heterogeneous data types (e.g. raster images and
   text captions),
2. your dataset is stored in a :class:`pandas.DataFrame` and different columns
   require different processing pipelines.

This example demonstrates how to use
:class:`~sklearn.compose.ColumnTransformer` on a dataset containing
different types of features. The choice of features is not particularly
helpful, but serves to illustrate the technique.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_20newsgroups
from sklearn.decomposition import PCA
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import classification_report
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer
from sklearn.svm import LinearSVC

##############################################################################
# 20 newsgroups dataset
# ---------------------
#
# We will use the :ref:`20 newsgroups dataset <20newsgroups_dataset>`, which
# comprises posts from newsgroups on 20 topics. This dataset is split
# into train and test subsets based on messages posted before and after
# a specific date. We will only use posts from 2 categories to speed up running
# time.

categories = ["sci.med", "sci.space"]
X_train, y_train = fetch_20newsgroups(
    random_state=1,
    subset="train",
    categories=categories,
    remove=("footers", "quotes"),
    return_X_y=True,
)
X_test, y_test = fetch_20newsgroups(
    random_state=1,
    subset="test",
    categories=categories,
    remove=("footers", "quotes"),
    return_X_y=True,
)

##############################################################################
# Each feature comprises meta information about that post, such as the subject,
# and the body of the news post.

print(X_train[0])

##############################################################################
# Creating transformers
# ---------------------
#
# First, we would like a transformer that extracts the subject and
# body of each post. Since this is a stateless transformation (does not
# require state information from training data), we can define a function that
# performs the data transformation then use
# :class:`~sklearn.preprocessing.FunctionTransformer` to create a scikit-learn
# transformer.


def subject_body_extractor(posts):
    # construct object dtype array with two columns
    # first column = 'subject' and second column = 'body'
    features = np.empty(shape=(len(posts), 2), dtype=object)
    for i, text in enumerate(posts):
        # temporary variable `_` stores '\n\n'
        headers, _, body = text.partition("\n\n")
        # store body text in second column
        features[i, 1] = body

        prefix = "Subject:"
        sub = ""
        # save text after 'Subject:' in first column
        for line in headers.split("\n"):
            if line.startswith(prefix):
                sub = line[len(prefix) :]
                break
        features[i, 0] = sub

    return features


subject_body_transformer = FunctionTransformer(subject_body_extractor)

##############################################################################
# We will also create a transformer that extracts the
# length of the text and the number of sentences.


def text_stats(posts):
    return [{"length": len(text), "num_sentences": text.count(".")} for text in posts]


text_stats_transformer = FunctionTransformer(text_stats)

##############################################################################
# Classification pipeline
# -----------------------
#
# The pipeline below extracts the subject and body from each post using
# ``SubjectBodyExtractor``, producing a (n_samples, 2) array. This array is
# then used to compute standard bag-of-words features for the subject and body
# as well as text length and number of sentences on the body, using
# ``ColumnTransformer``. We combine them, with weights, then train a
# classifier on the combined set of features.

pipeline = Pipeline(
    [
        # Extract subject & body
        ("subjectbody", subject_body_transformer),
        # Use ColumnTransformer to combine the subject and body features
        (
            "union",
            ColumnTransformer(
                [
                    # bag-of-words for subject (col 0)
                    ("subject", TfidfVectorizer(min_df=50), 0),
                    # bag-of-words with decomposition for body (col 1)
                    (
                        "body_bow",
                        Pipeline(
                            [
                                ("tfidf", TfidfVectorizer()),
                                ("best", PCA(n_components=50, svd_solver="arpack")),
                            ]
                        ),
                        1,
                    ),
                    # Pipeline for pulling text stats from post's body
                    (
                        "body_stats",
                        Pipeline(
                            [
                                (
                                    "stats",
                                    text_stats_transformer,
                                ),  # returns a list of dicts
                                (
                                    "vect",
                                    DictVectorizer(),
                                ),  # list of dicts -> feature matrix
                            ]
                        ),
                        1,
                    ),
                ],
                # weight above ColumnTransformer features
                transformer_weights={
                    "subject": 0.8,
                    "body_bow": 0.5,
                    "body_stats": 1.0,
                },
            ),
        ),
        # Use a SVC classifier on the combined features
        ("svc", LinearSVC(dual=False)),
    ],
    verbose=True,
)

##############################################################################
# Finally, we fit our pipeline on the training data and use it to predict
# topics for ``X_test``. Performance metrics of our pipeline are then printed.

pipeline.fit(X_train, y_train)
y_pred = pipeline.predict(X_test)
print("Classification report:\n\n{}".format(classification_report(y_test, y_pred)))
